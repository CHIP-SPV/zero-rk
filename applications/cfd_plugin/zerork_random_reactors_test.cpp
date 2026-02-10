#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

#include "mpi.h"
#include "zerork_cfd_plugin.h"

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
getRandomDoubles(unsigned int N, double const min, double const max)
{
  std::vector<double> x(N, 0.0);
  std::default_random_engine gen;
  gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<> rand(min, max);

  for (int i = 0; i < N; ++i)
    x[i] = rand(gen);
  return x;
}

std::vector<double>
getRandomY(unsigned int N,
	   std::vector<double> const &Yu,
	   std::vector<double> const &Yb)
{
  unsigned int ns = Yu.size();
  std::vector<double> Y(N*ns, 0.0);

  std::default_random_engine gen;
  gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<> rand(0, 1.0);

  std::vector<double> Ymix(ns, 0.0);
  for (int i = 0; i < N; ++i) {
    double const phi = rand(gen);
    for (int k = 0l; k < ns; ++k) {
      Ymix[k] = (1.0-phi)*Yu[k] + phi*Yb[k];
      Y[i*ns + k] = Ymix[k];
    }    
  }

  return Y;
}

void initialize_zerork(zerork_handle &zrm_handle)
{
  std::string const zerork_file = "zerork.yml";
  std::string const chem_file = "chem.inp";
  std::string const therm_file = "therm.dat";

  //std::cout<<"Reading chemistry using files: "<<chem_file<<" and "<<therm_file<<std::endl;
  
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

  unsigned int const ns = 20;
  std::vector<double> const Yf{0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  std::vector<double> const Yo{0.23, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.77};

  size_t const N_total = 2e9;  
  unsigned int const N = N_total/nranks;

  if (rank == 0) std::cout<<"Total Reactors, N_total = " << N_total <<std::endl;
  if (rank == 0) std::cout<<"No. of Ranks, nranks = " << nranks <<std::endl;  
  if (rank == 0) std::cout<<"Reactors per Rank, N = " << N <<std::endl;
  
  auto T = getRandomDoubles(N, 300.0, 2100.0);
  auto P = getRandomDoubles(N, 101325.0, 101325.0);
  auto Y = getRandomY(N, Yf, Yo);

  double t = 0;
  double const dt = 1e-8;
  unsigned int const n_steps = 1000;

  std::default_random_engine gen;
  gen.seed(std::chrono::system_clock::now().time_since_epoch().count());

  if (rank == 0) std::cout<<"Begin computation for n_steps = " << n_steps <<std::endl;

  double T_min_g, T_max_g;
  
  for (int step = 0; step < n_steps; ++step) {

    zerork_reactor_solve(step, t, dt, N, T.data(), P.data(), Y.data(), zrm_handle);

    t += dt;
    
    auto const [min_it, max_it] = std::minmax_element(T.begin(), T.end());
    double Tmin = *min_it;
    double Tmax = *max_it;
    MPI_Allreduce(&Tmin, &T_min_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&Tmax, &T_max_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    if (rank == 0) std::cout<<"step = "<<step<<", t = " <<t<<"s : N = "<<N<<", Tmin = " << T_min_g << ", Tmax = " << T_max_g << std::endl;
  }

  MPI_Finalize();

  if (rank == 0) std::cout<<"End of computation!"<<std::endl;
  
  return 1;
}
