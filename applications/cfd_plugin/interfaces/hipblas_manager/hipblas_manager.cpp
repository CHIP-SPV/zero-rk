#include "hip/hip_runtime.h"
#include "hip/hip_complex.h"
#include "../../gpu_err_check.h"
#include "hipblas_manager.h"

#include <chrono>
#include <cstdio>


static const char* hipblas_status_str(hipblasStatus_t s) {
  switch (s) {
    case HIPBLAS_STATUS_SUCCESS:          return "HIPBLAS_STATUS_SUCCESS";
    case HIPBLAS_STATUS_NOT_INITIALIZED:  return "HIPBLAS_STATUS_NOT_INITIALIZED";
    case HIPBLAS_STATUS_ALLOC_FAILED:     return "HIPBLAS_STATUS_ALLOC_FAILED";
    case HIPBLAS_STATUS_INVALID_VALUE:    return "HIPBLAS_STATUS_INVALID_VALUE";
    case HIPBLAS_STATUS_MAPPING_ERROR:    return "HIPBLAS_STATUS_MAPPING_ERROR";
    case HIPBLAS_STATUS_EXECUTION_FAILED: return "HIPBLAS_STATUS_EXECUTION_FAILED";
    case HIPBLAS_STATUS_INTERNAL_ERROR:   return "HIPBLAS_STATUS_INTERNAL_ERROR";
    case HIPBLAS_STATUS_NOT_SUPPORTED:    return "HIPBLAS_STATUS_NOT_SUPPORTED";
    case HIPBLAS_STATUS_ARCH_MISMATCH:    return "HIPBLAS_STATUS_ARCH_MISMATCH";
    case HIPBLAS_STATUS_HANDLE_IS_NULLPTR:return "HIPBLAS_STATUS_HANDLE_IS_NULLPTR";
    case HIPBLAS_STATUS_INVALID_ENUM:     return "HIPBLAS_STATUS_INVALID_ENUM";
    default:                              return "HIPBLAS_STATUS_UNKNOWN";
  }
}

#define HIP_CHECK(cmd) do {                                             \
  hipblasStatus_t s = (cmd); \
if (s != HIPBLAS_STATUS_SUCCESS) {				       \
    fprintf(stderr, "HIP error at %s:%d\n",                         \
            __FILE__, __LINE__);                 \
    exit(1);                                                            \
  }                                                                     \
} while (0)



template<typename T>
hipblas_manager<T>::hipblas_manager() :
  n_(-1),
  num_batches_(-1),
  factored_(false),
  info_dev_(),
  tmp_dev_(),
  matrix_inverse_dev_(),
  matrix_inverse_pointers_dev_(),
  matrix_pointers_dev_(),
  tmp_pointers_dev_()
{
  hipblasCreate(&hipblas_handle_);
}

template<typename T>
hipblas_manager<T>::~hipblas_manager()
{
  if(factored_) {
    FreeDeviceMemory();
    hipblasDestroy(hipblas_handle_);
  }
}

template<typename T>
void hipblas_manager<T>::setup_memory()
{
  if(factored_) {
    FreeDeviceMemory();
  }
  AllocateDeviceMemory();
}


template<typename T>
void hipblas_manager<T>::AllocateDeviceMemory()
{
  hipDeviceSynchronize();
  gpu_err_check(hipGetLastError());

  matrix_inverse_dev_.resize(n_*n_*num_batches_);
  matrix_inverse_pointers_dev_.resize(num_batches_);
  matrix_pointers_dev_.resize(num_batches_);
  info_dev_.resize(num_batches_);
  tmp_dev_.resize(num_batches_*n_);
  tmp_pointers_dev_.resize(num_batches_);

  data_ptrs_.resize(num_batches_);
  tmp_ptrs_.resize(num_batches_);

  for(int j = 0; j < num_batches_; ++j) {
    data_ptrs_[j] = thrust::raw_pointer_cast(matrix_inverse_dev_.data()) + j*n_*n_;
  }
  hipMemcpy(thrust::raw_pointer_cast(matrix_inverse_pointers_dev_.data()), data_ptrs_.data(), sizeof(T*)*num_batches_, hipMemcpyHostToDevice);
  gpu_err_check(hipGetLastError());

  for(int j = 0; j < num_batches_; ++j) {
    tmp_ptrs_[j] = thrust::raw_pointer_cast(tmp_dev_.data()) + j*n_;
  }
  hipMemcpy(thrust::raw_pointer_cast(tmp_pointers_dev_.data()), tmp_ptrs_.data(), sizeof(T*)*num_batches_, hipMemcpyHostToDevice);
  gpu_err_check(hipGetLastError());

  // Identity pivots for the batched getrs (which requires a non-NULL ipiv even
  // though the factorization is no-pivot). ipiv[i] = i+1 (1-based) means "no row
  // swap". The batched getrs reads a contiguous n*num_batches pivot array, so
  // build one identity block per batch. Built once here for a given shape.
  ipiv_dev_.resize(n_*num_batches_);
  std::vector<int> ipiv_host(n_*num_batches_);
  for(int j = 0; j < num_batches_; ++j) {
    for(int i = 0; i < n_; ++i) {
      ipiv_host[j*n_ + i] = i + 1;
    }
  }
  hipMemcpy(thrust::raw_pointer_cast(ipiv_dev_.data()), ipiv_host.data(), sizeof(int)*n_*num_batches_, hipMemcpyHostToDevice);
  gpu_err_check(hipGetLastError());
}

template<typename T>
void hipblas_manager<T>::FreeDeviceMemory()
{
}

template<>
void hipblas_manager<double>::getrf_batched() {
  int lda = n_;
  int* ipiv = NULL; //Turns off pivoting
  printf("[hipblas_manager] Dgetrf_batched: n=%d lda=%d num_batches=%d\n",
         n_, lda, num_batches_);
  hipDeviceSynchronize();
  auto t0 = std::chrono::high_resolution_clock::now();
  HIP_CHECK(hipblasDgetrfBatched(hipblas_handle_, n_,
                       thrust::raw_pointer_cast(matrix_pointers_dev_.data()), lda,
                       ipiv, thrust::raw_pointer_cast(info_dev_.data()), num_batches_));
  hipDeviceSynchronize();
  auto t1 = std::chrono::high_resolution_clock::now();
  printf("[hipblas_manager] Dgetrf_batched: %.6f ms\n",
         std::chrono::duration<double, std::milli>(t1 - t0).count());
}

template<>
void hipblas_manager<double>::getri_batched() {
  int lda = n_;
  int* ipiv = NULL; //Turns off pivoting
  int ldc = n_;
  double* const* const_matrix_pointers_dev = (double* const*) thrust::raw_pointer_cast(matrix_pointers_dev_.data());
  printf("[hipblas_manager] Dgetri_batched: n=%d lda=%d ldc=%d num_batches=%d\n",
         n_, lda, ldc, num_batches_);
  hipDeviceSynchronize();
  auto t0 = std::chrono::high_resolution_clock::now();
  HIP_CHECK(hipblasDgetriBatched(hipblas_handle_, n_, const_matrix_pointers_dev,
                       lda, ipiv, thrust::raw_pointer_cast(matrix_inverse_pointers_dev_.data()),
				 ldc, thrust::raw_pointer_cast(info_dev_.data()), num_batches_));
  hipDeviceSynchronize();
  auto t1 = std::chrono::high_resolution_clock::now();
  printf("[hipblas_manager] Dgetri_batched: %.6f ms\n",
         std::chrono::duration<double, std::milli>(t1 - t0).count());
}

template<>
void hipblas_manager<hipDoubleComplex>::getrf_batched() {
  int lda = n_;
  int* ipiv = NULL; //Turns off pivoting
  hipblasDoubleComplex* const* const_matrix_pointers_dev = (hipblasDoubleComplex* const*) thrust::raw_pointer_cast(matrix_pointers_dev_.data());
  printf("[hipblas_manager] Zgetrf_batched: n=%d lda=%d num_batches=%d\n",
         n_, lda, num_batches_);
  hipDeviceSynchronize();
  auto t0 = std::chrono::high_resolution_clock::now();
  hipblasZgetrfBatched(hipblas_handle_, n_,
                       const_matrix_pointers_dev, lda,
                       ipiv, thrust::raw_pointer_cast(info_dev_.data()), num_batches_);
  hipDeviceSynchronize();
  auto t1 = std::chrono::high_resolution_clock::now();
  printf("[hipblas_manager] Zgetrf_batched: %.6f ms\n",
         std::chrono::duration<double, std::milli>(t1 - t0).count());
}

template<>
void hipblas_manager<hipDoubleComplex>::getri_batched() {
  int lda = n_;
  int* ipiv = NULL; //Turns off pivoting
  int ldc = n_;
  hipblasDoubleComplex* const* const_matrix_pointers_dev = (hipblasDoubleComplex* const*) thrust::raw_pointer_cast(matrix_pointers_dev_.data());
  hipblasDoubleComplex* const* const_matrix_inverse_pointers_dev = (hipblasDoubleComplex* const*) thrust::raw_pointer_cast(matrix_inverse_pointers_dev_.data());
  printf("[hipblas_manager] Zgetri_batched: n=%d lda=%d ldc=%d num_batches=%d\n",
         n_, lda, ldc, num_batches_);
  hipDeviceSynchronize();
  auto t0 = std::chrono::high_resolution_clock::now();
  hipblasZgetriBatched(hipblas_handle_, n_, const_matrix_pointers_dev,
                       lda, ipiv, const_matrix_inverse_pointers_dev,
                       ldc, thrust::raw_pointer_cast(info_dev_.data()), num_batches_);
  hipDeviceSynchronize();
  auto t1 = std::chrono::high_resolution_clock::now();
  printf("[hipblas_manager] Zgetri_batched: %.6f ms\n",
         std::chrono::duration<double, std::milli>(t1 - t0).count());
}


template<typename T>
int hipblas_manager<T>::factor_invert(int num_batches, int n, T* values) {
  if(n != n_ || num_batches != num_batches_) {
    n_ = n;
    num_batches_ = num_batches;
    setup_memory();
  }
  if(values == NULL) {
    return 1;
  }

  bool need_tx = false;
  for(int j = 0; j < num_batches_; ++j) {
    if(data_ptrs_[j] != values + j*n_*n_) {
      data_ptrs_[j] = values + j*n_*n_;
      need_tx = true;
    }
  }
  if(need_tx) {
    hipMemcpy(thrust::raw_pointer_cast(matrix_pointers_dev_.data()), data_ptrs_.data(), sizeof(T*)*num_batches_, hipMemcpyHostToDevice);
  }

  this->getrf_batched();
  this->getri_batched();

  int ierr = 0;
#ifdef ZERORK_FULL_DEBUG
  info_.resize(num_batches_);
  gpu_err_check(hipMemcpy(info_.data(), thrust::raw_pointer_cast(info_dev_.data()), num_batches_*sizeof(int), hipMemcpyDeviceToHost));
  //Check for errors
  // factor_error > 0, singular matrix, zero diagonal at row,col = factor_error
  // factor_error = 0, success
  // factor_error < 0, illegal input
  for(int i=0; i < num_batches_; ++i) {
    if(info_[i]!=0) {
      ierr = info_[i];
      break;
    }
  }
#endif

  factored_ = true;
  return ierr;
}

template<typename T>
int hipblas_manager<T>::factor_lu(int num_batches, int n, T* values) {
  if(n != n_ || num_batches != num_batches_) {
    n_ = n;
    num_batches_ = num_batches;
    setup_memory();
  }
  if(values == NULL) {
    return 1;
  }

  bool need_tx = false;
  for(int j = 0; j < num_batches_; ++j) {
    if(data_ptrs_[j] != values + j*n_*n_) {
      data_ptrs_[j] = values + j*n_*n_;
      need_tx = true;
    }
  }
  if(need_tx) {
    hipMemcpy(thrust::raw_pointer_cast(matrix_pointers_dev_.data()), data_ptrs_.data(), sizeof(T*)*num_batches_, hipMemcpyHostToDevice);
  }

  this->getrf_batched();

  int ierr = 0;
#ifdef ZERORK_FULL_DEBUG
  info_.resize(num_batches_);
  gpu_err_check(hipMemcpy(info_.data(), thrust::raw_pointer_cast(info_dev_.data()), num_batches_*sizeof(int), hipMemcpyDeviceToHost));
  //Check for errors
  // factor_error > 0, singular matrix, zero diagonal at row,col = factor_error
  // factor_error = 0, success
  // factor_error < 0, illegal input
  for(int i=0; i < num_batches_; ++i) {
    if(info_[i]!=0) {
      ierr = info_[i];
      break;
    }
  }
#endif

  factored_ = true;
  return ierr;
}


//The following modified from cuda sdk-5.0
#define TRANSPOSE_TILE_DIM    32
#define TRANSPOSE_BLOCK_ROWS  8

template<typename T>
static __global__ void HIPBLAS_MANAGER_TransposeNoBankConflicts(T *odata, const T *idata, const int width, const int height)
{
    __shared__ T tile[TRANSPOSE_TILE_DIM][TRANSPOSE_TILE_DIM+1];
    int xIndex,yIndex,index_in,index_out;

    xIndex = blockIdx.x * TRANSPOSE_TILE_DIM + threadIdx.x;
    yIndex = blockIdx.y * TRANSPOSE_TILE_DIM + threadIdx.y;
    index_in = xIndex + (yIndex)*width;

    for (int i=0; i<TRANSPOSE_TILE_DIM; i+=TRANSPOSE_BLOCK_ROWS)
    {
        if(xIndex < width && yIndex+i < height){
        tile[threadIdx.y+i][threadIdx.x] = idata[index_in+i*width];}
    }

    __syncthreads();

    xIndex = blockIdx.y * TRANSPOSE_TILE_DIM + threadIdx.x;
    yIndex = blockIdx.x * TRANSPOSE_TILE_DIM + threadIdx.y;
    index_out = xIndex + (yIndex)*height;

    for (int i=0; i<TRANSPOSE_TILE_DIM; i+=TRANSPOSE_BLOCK_ROWS)
    {
        if(yIndex+i < width && xIndex < height){
        odata[index_out+i*height] = tile[threadIdx.x][threadIdx.y+i];}
    }
}

template<typename T>
void hipblas_manager<T>::gpu_transpose(T* odata, const T* idata, const int width, const int height)
{
    // Put df/dy in "normal" order
    dim3 nBlocks2D,nThreads2D;
    nThreads2D.x = TRANSPOSE_TILE_DIM;
    nThreads2D.y = TRANSPOSE_BLOCK_ROWS;
    nBlocks2D.x = (width+TRANSPOSE_TILE_DIM-1)/TRANSPOSE_TILE_DIM;
    nBlocks2D.y = (height+TRANSPOSE_TILE_DIM-1)/TRANSPOSE_TILE_DIM;
    HIPBLAS_MANAGER_TransposeNoBankConflicts<T><<<nBlocks2D,nThreads2D>>>(odata,idata,width,height);
#ifdef ZERORK_FULL_DEBUG
    gpu_err_check( hipPeekAtLastError() );
    gpu_err_check( hipDeviceSynchronize() );
#endif
}

namespace {
template<typename T>
void __global__ HIPBLAS_MANAGER_gpu_bdmv_kernel
(
    const int mtx_block_size,
    const int num_mtx_blocks,
    const T* A_dev,
    const T* X_dev ,
    T * Y_dev
)
{
  int tidx = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  for( ; tidx < num_mtx_blocks*mtx_block_size; tidx += stride)
  {
    int local_row   = tidx % mtx_block_size;
    int local_block = tidx / mtx_block_size;
    T Y_dev_accum = 0.0;
    for(int i = 0; i < mtx_block_size; ++i) //columns
    {
      int data_idx = mtx_block_size*mtx_block_size*local_block + mtx_block_size*i + local_row;
      Y_dev_accum += A_dev[data_idx]*X_dev[i+local_block*mtx_block_size];
    }
    Y_dev[local_row+local_block*mtx_block_size] = Y_dev_accum;
  }
}

template<>
void __global__ HIPBLAS_MANAGER_gpu_bdmv_kernel
(
    const int mtx_block_size,
    const int num_mtx_blocks,
    const hipDoubleComplex* A_dev,
    const hipDoubleComplex* X_dev ,
    hipDoubleComplex * Y_dev
)
{
  int tidx = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  for( ; tidx < num_mtx_blocks*mtx_block_size; tidx += stride)
  {
    int local_row   = tidx % mtx_block_size;
    int local_block = tidx / mtx_block_size;
    hipDoubleComplex Y_dev_accum = make_hipDoubleComplex(0.0,0.0);
    for(int i = 0; i < mtx_block_size; ++i) //columns
    {
      int data_idx = mtx_block_size*mtx_block_size*local_block + mtx_block_size*i + local_row;
      //Y_dev_accum += A_dev[data_idx]*X_dev[i+local_block*mtx_block_size];
      Y_dev_accum = hipCadd(Y_dev_accum, hipCmul(A_dev[data_idx],X_dev[i+local_block*mtx_block_size]));
    }
    Y_dev[local_row+local_block*mtx_block_size] = Y_dev_accum;
  }
}
} //anonymous namespace

template<typename T>
int hipblas_manager<T>::gpu_bdmv(int n, int nbatch, T* A_dev, T* B_dev, T* Y_dev)
{
  int threads = std::min(n*nbatch,1024);
  int blocks=(nbatch*n+threads-1)/threads;
  HIPBLAS_MANAGER_gpu_bdmv_kernel<T><<<blocks,threads>>>(n, nbatch, A_dev, B_dev, Y_dev);
#ifdef ZERORK_FULL_DEBUG
  gpu_err_check(hipPeekAtLastError());
  gpu_err_check(hipDeviceSynchronize());
#endif
  return 0;  
}

template<typename T>
int hipblas_manager<T>::solve_invert(int num_batches, int n, const T* rhs, T* soln) {
  if(n != n_ || num_batches != num_batches_) {
    return 1;
  }

  // Transpose rhs into soln
  this->gpu_transpose(soln,rhs,num_batches_,n_);

  // Block-diagonal matrix vector multiplication
  this->gpu_bdmv(n_, num_batches_, thrust::raw_pointer_cast(matrix_inverse_dev_.data()), soln, thrust::raw_pointer_cast(tmp_dev_.data()));

  // Put tmp back into block order
  this->gpu_transpose(soln,thrust::raw_pointer_cast(tmp_dev_.data()),n_,num_batches_);

  return(0);
}

template<>
void hipblas_manager<double>::getrs_batched() {
  int lda = n_;
  int ldb = n_;
  int info = 0;
  // ipiv is the identity pivot array (getrs rejects NULL) laid out as one
  // n-element block per batch; the factorization is no-pivot.
  int* ipiv = thrust::raw_pointer_cast(ipiv_dev_.data());
  double* const* const_matrix_pointers_dev = (double* const*) thrust::raw_pointer_cast(matrix_pointers_dev_.data());
  printf("[hipblas_manager] Dgetrs_batched: n=%d nrhs=%d lda=%d ldb=%d num_batches=%d\n",
         n_, 1, lda, ldb, num_batches_);
  hipDeviceSynchronize();
  auto t0 = std::chrono::high_resolution_clock::now();
  hipblasDgetrsBatched(hipblas_handle_, HIPBLAS_OP_N, n_, 1,
                       const_matrix_pointers_dev, lda,
                       ipiv, thrust::raw_pointer_cast(tmp_pointers_dev_.data()), ldb, &info, num_batches_);
  hipDeviceSynchronize();
  auto t1 = std::chrono::high_resolution_clock::now();
  printf("[hipblas_manager] Dgetrs_batched: %.6f ms\n",
         std::chrono::duration<double, std::milli>(t1 - t0).count());
}

template<>
void hipblas_manager<hipDoubleComplex>::getrs_batched() {
  int* ipiv = NULL; //Turns off pivoting
  int lda = n_;
  int ldb = n_;
  int info = 0;
  hipblasDoubleComplex* const* const_matrix_pointers_dev = (hipblasDoubleComplex* const*) thrust::raw_pointer_cast(matrix_pointers_dev_.data());
  hipblasDoubleComplex* const* const_tmp_pointers_dev = (hipblasDoubleComplex* const*) thrust::raw_pointer_cast(tmp_pointers_dev_.data());
  printf("[hipblas_manager] Zgetrs_batched: n=%d nrhs=%d lda=%d ldb=%d num_batches=%d\n",
         n_, 1, lda, ldb, num_batches_);
  hipDeviceSynchronize();
  auto t0 = std::chrono::high_resolution_clock::now();
  hipblasZgetrsBatched(hipblas_handle_, HIPBLAS_OP_N, n_, 1,
                       const_matrix_pointers_dev, lda,
                       ipiv, const_tmp_pointers_dev, ldb, &info, num_batches_);
  hipDeviceSynchronize();
  auto t1 = std::chrono::high_resolution_clock::now();
  printf("[hipblas_manager] Zgetrs_batched: %.6f ms\n",
         std::chrono::duration<double, std::milli>(t1 - t0).count());
}

template<typename T>
int hipblas_manager<T>::solve_lu(int num_batches, int n, const T* rhs, T* soln) {
  if(n != n_ || num_batches != num_batches_) {
    return 1;
  }

  // Transpose rhs into tmp_dev_
  this->gpu_transpose(thrust::raw_pointer_cast(tmp_dev_.data()),rhs,num_batches_,n_);

  // HIPBLAS forward and back substitution
  this->getrs_batched();

  // Put tmp back into block order
  this->gpu_transpose(soln,thrust::raw_pointer_cast(tmp_dev_.data()),n_,num_batches_);

  return(0);
}

template class hipblas_manager<double>;
template class hipblas_manager<hipDoubleComplex>;

