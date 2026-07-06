#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <cassert>

#include "mpi.h"
#include "zerork_cfd_plugin.h"

struct Block2D {
  int xi, xj; // Start and end X indices (inclusive)
  int yi, yj; // Start and end Y indices (inclusive)
  int size; // number of elements
};

// Helper function to find the best 2D factors for N closest to the grid aspect ratio
std::pair<int, int>
findGridFactors(int N, int Nx, int Ny) {
  int fx = 1, fy = N;
  double targetRatio = static_cast<double>(Nx) / Ny;
  double bestDiff = std::abs((static_cast<double>(fx) / fy) - targetRatio);

  for (int i = 1; i <= std::sqrt(N); ++i) {
    if (N % i == 0) {
      int f1 = i;
      int f2 = N / i;

      // Test f1 as X, f2 as Y
      double diff1 = std::abs((static_cast<double>(f1) / f2) - targetRatio);
      if (diff1 < bestDiff) {
	bestDiff = diff1;
	fx = f1;
	fy = f2;
      }

      // Test f2 as X, f1 as Y
      double diff2 = std::abs((static_cast<double>(f2) / f1) - targetRatio);
      if (diff2 < bestDiff) {
	bestDiff = diff2;
	fx = f2;
	fy = f1;
      }
    }
  }
  return {fx, fy};
}

std::vector<Block2D>
sliceGrid2D(int Nx, int Ny, int N) {
  std::vector<Block2D> blocks;
  if (N <= 0 || Nx <= 0 || Ny <= 0) return blocks;

  // 1. Determine how many splits to make along X and Y
  auto [Fx, Fy] = findGridFactors(N, Nx, Ny);
    
  int baseWidth = Nx / Fx;
  int remX = Nx % Fx;

  int baseHeight = Ny / Fy;
  int remY = Ny % Fy;

  // 2. Compute indices and attach leftovers to the last blocks
  int currentX = 0;
  for (int i = 0; i < Fx; ++i) {
    int width = baseWidth + (i == Fx - 1 ? remX : 0);
    int nxa = currentX;
    int nxb = currentX + width - 1;

    int currentY = 0;
    for (int j = 0; j < Fy; ++j) {
      int height = baseHeight + (j == Fy - 1 ? remY : 0);
      int nya = currentY;
      int nyb = currentY + height - 1;

      int const size = (nxb - nxa + 1)*(nyb - nya + 1);
      blocks.push_back({nxa, nxb, nya, nyb, size});
      currentY += height;
    }
    currentX += width;
  }
  return blocks;
}

template <typename T>
void
printVector(std::vector<T> const &values)
{
  int const n = values.size();
  std::cout<<"[";
  for (int i = 0; i < n; ++i)
    std::cout<<values[i]<<" ";
  std::cout<<"]"<<std::endl;
}


std::vector<double>
getBoxTemperature(Block2D block, unsigned int Nx, unsigned int Ny, unsigned int Nz,
		  double const Tmin, double const Tmax)
{
  size_t const N = block.size*Nz;
  std::vector<double> T(N, Tmin);
  if (Tmin == Tmax) return T;
  
  double alpha = 0;
  size_t index = 0;
  for (int i = block.xi; i <= block.xj; ++i) {
    for (int j = block.yi; j <= block.yj; ++j) {
      alpha = (double)j/(Ny-1.0);
      for (int k = 0; k < Nz; ++k) {
	T[index++] = Tmin + (Tmax - Tmin)*alpha;
      }
    }
  }
  return T;
}

std::vector<double>
getBoxPressure(Block2D block, unsigned int Nx, unsigned int Ny, unsigned int Nz,
	       double const Pmin, double const Pmax)
{
  size_t const N = block.size*Nz;
  std::vector<double> P(N, Pmin);
  if (Pmin == Pmax) return P;

  double alpha = 0;
  size_t index = 0;
  for (int i = block.xi; i <= block.xj; ++i) {
    for (int j = block.yi; j <= block.yj; ++j) {
      alpha = (double)j/(Ny-1.0);
      for (int k = 0; k < Nz; ++k) {
	P[index++] = Pmin + (Pmax - Pmin)*alpha;
      }
    }
  }
  return P;
}
      
std::vector<double>
getBoxMassFractions(Block2D block, unsigned int Nx, unsigned int Ny, unsigned int Nz,
		    std::vector<double> const &Yu, std::vector<double> const &Yb)
{
  size_t const N = block.size*Nz;  
  unsigned int ns = Yu.size();
  std::vector<double> Y(N*ns, 0.0);

  std::vector<double> Ymix(ns, 0.0);
  double phi = 0;
  size_t index = 0;
  for (int i = block.xi; i <= block.xj; ++i) {
    phi = (double)i/(Nx-1.0);    
    for (int j = block.yi; j <= block.yj; ++j) {
      for (int k = 0; k < Nz; ++k) {
	for (int s = 0; s < ns; ++s) {
	  Ymix[s] = (1.0-phi)*Yu[s] + phi*Yb[s];
	  Y[index*ns + s] = Ymix[s];
	}
	index++;
      }
    }
  }
  
  return Y;
}

void initialize_zerork(zerork_handle &zrm_handle)
{
  std::string const zerork_file = "zerork.yml";
  std::string const chem_file = "chem.inp";
  std::string const therm_file = "therm.dat";

  zrm_handle = zerork_reactor_init();

  zerork_reactor_read_options_file(zerork_file.c_str(), zrm_handle);
  zerork_reactor_set_mechanism_files(chem_file.c_str(), therm_file.c_str(), zrm_handle);
  zerork_reactor_load_mechanism(zrm_handle);

  zerork_reactor_set_int_option("constant_volume", 1, zrm_handle);
}

int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  int rank = 0;
  int nranks = 1;
  
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nranks);

  zerork_handle zrm_handle;
  initialize_zerork(zrm_handle);

  // unsigned int const ns = 6; //2-step propane: C3H8 O2 H2O CO CO2 N2
  // std::vector<double> const Yf{1.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  // std::vector<double> const Yo{0.00, 0.23, 0.00, 0.00, 0.00, 0.77};

  unsigned int const ns = 9; //H2 mech
  std::vector<double> const Yf{1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  std::vector<double> const Yo{0.00, 0.23, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.77};
  
  // unsigned int const ns = 20;
  // std::vector<double> const Yf{0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  //   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  // std::vector<double> const Yo{0.23, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  //   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.77};

  size_t const Po = 0;
  size_t const Nx = 100*(Po+1);
  size_t const Ny = 100*(Po+1);
  size_t const Nz = 5*(Po+1);
  size_t const N_total = Nx*Ny*Nz;

  auto const blocks = sliceGrid2D(Nx, Ny, nranks);
  if (rank == 0) {
    std::cout<<"Box Reactor Size = [" << Nx << "x" << Ny << "x" << Nz <<"]" << std::endl;
    std::cout<<"Number of Ranks  = " << nranks <<std::endl;    
    std::cout<<"Box Partitioning: [";
    size_t sum = 0;
    for (int i = 0; i < nranks; ++i) {
      //std::cout<<"Block2D[" << i << "] = " << blocks[i].xi << ":" << blocks[i].xj << " " << blocks[i].yi << ":" << blocks[i].yj << std::endl;
      size_t const N = blocks[i].size*Nz;      
      std::cout<<i<<"="<<N<<" ";
      sum += N;
    }
    std::cout<<"]"<<std::endl;
    
    std::cout<<"Sum of Partitions = " << sum << std::endl;
    std::cout<<"Total Reactors Nx*Ny*Nz = " << N_total <<std::endl;
  }
  
  auto T = getBoxTemperature(blocks[rank], Nx, Ny, Nz, 600.0, 1800.0);
  auto P = getBoxPressure(blocks[rank], Nx, Ny, Nz, 101325.0, 101325.0);
  auto Y = getBoxMassFractions(blocks[rank], Nx, Ny, Nz, Yf, Yo);
  
  double t = 0;
  double const dt = 1e-8; //2e-7; //1e-8;
  unsigned int const n_steps = 2000; //200000;

  std::default_random_engine gen;
  gen.seed(std::chrono::system_clock::now().time_since_epoch().count());

  if (rank == 0) std::cout<<"Begin Box Reactor computation for n_steps = " << n_steps <<std::endl;
  if (rank == 0) std::cout<<"Number of species = " << ns <<std::endl;
  
  double T_min_g, T_max_g;

  std::ofstream results;
  if (rank == 0) {
    results.open("box_reactor_results.csv");
    results << "Step,Time,T_min,T_max" << std::endl;
  }
  
  for (int step = 0; step < n_steps; ++step) {

    zerork_reactor_solve(step, t, dt, T.size(), T.data(), P.data(), Y.data(), zrm_handle);

    t += dt;
    
    auto const [min_it, max_it] = std::minmax_element(T.begin(), T.end());
    double Tmin = *min_it;
    double Tmax = *max_it;
    MPI_Allreduce(&Tmin, &T_min_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&Tmax, &T_max_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    if (rank == 0) std::cout<<"step = "<<step<<", t = " <<t<<"s : N_total = "<<N_total<<", Tmin = " << T_min_g << ", Tmax = " << T_max_g << std::endl;
    if (rank == 0) results << step+1 << "," << t << "," << T_min_g << "," << T_max_g << std::endl;
  }

  if (rank == 0) results.close();
  
  MPI_Finalize();

  if (rank == 0) std::cout<<"End of computation!"<<std::endl;
  
  return 1;
}
