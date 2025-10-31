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


std::vector<unsigned int>
getRandomInts(unsigned int N, unsigned int const min, unsigned int const max)
{
  std::vector<unsigned int> I(N, 0.0);
  std::default_random_engine gen;
  gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<> rand(min, max-1);

  for (int k = 0; k < N; ++k)
    I[k] = rand(gen);
  return I;
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
getRandomSamples(std::vector<unsigned int> const &indices,
		 std::vector<double> const &values)
{
  unsigned int const n = indices.size();
  std::vector<double> x(n);
  
  for (int i = 0; i < n; ++i)
    x[i] = values[indices[i]];
  return x;
}

std::vector<double>
getRandomSamplesY(std::vector<unsigned int> const &indices,
		  unsigned int const ns,
		  std::vector<double> const &Y)
{
  unsigned int const n = indices.size();
  std::vector<double> Ys(n*ns);
  
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < ns; ++k)
      Ys[i*ns + k] = Y[indices[i]*ns + k];
  }
  return Ys;
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

void
update_solution_vectors(std::vector<unsigned int> const &indices, unsigned int const ns,
			std::vector<double> &T, std::vector<double> &P, std::vector<double> &Y,
			std::vector<double> const &Ts, std::vector<double> const &Ps, std::vector<double> const &Ys)
{
  unsigned int const n = indices.size();

  for (int i = 0; i < n; ++i) {
    T[indices[i]] = Ts[i];
    P[indices[i]] = Ps[i];
    for (int k = 0; k < ns; ++k) {
      Y[indices[i]*ns + k ] = Ys[i*ns + k];
    }
  }
}

void initialize_zerork(zerork_handle &zrm_handle)
{
  std::string const zerork_file = "zerork.yml";
  std::string const chem_file = "chem.inp";
  std::string const therm_file = "therm.dat";

  std::cout<<"Reading chemistry using files: "<<chem_file<<" and "<<therm_file<<std::endl;
  
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

  zerork_handle zrm_handle;
  initialize_zerork(zrm_handle);

  // unsigned int const ns = 10; //H2 mech species order: h,h2,o,o2,oh,h2o,n2,ho2,h2o2,ar
  // std::vector<double> const Yu{0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // std::vector<double> const Yb{0.0, 0.0, 0.0, 0.23, 0.0, 0.0, 0.77, 0.0, 0.0, 0.0};

  unsigned int const ns = 6; //CH4 BFER mechanism: O2, H2O, CH4, CO, CO2, N2
  std::vector<double> const Yu{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
  std::vector<double> const Yb{0.23, 0.0, 0.0, 0.0, 0.0, 0.77};
  
  unsigned int const N = 3e6;
  auto T = getRandomDoubles(N, 300.0, 2000.0);
  auto P = getRandomDoubles(N, 101325.0, 101325.0);
  auto Y = getRandomY(N, Yu, Yb);

  double t = 0;
  double const dt = 5.207e-9;
  unsigned int const n_steps = 100;

  std::default_random_engine gen;
  gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<> samples(2e6,N);

  std::cout<<"Begin computation for n_steps = " << n_steps <<std::endl;
  
  for (int step = 0; step < n_steps; ++step) {

    unsigned int M = samples(gen);

    // select M random indices
    auto const indices = getRandomInts(M, 0, N);

    // extract T,P,Y
    auto Ts = getRandomSamples(indices, T);
    auto Ps = getRandomSamples(indices, P);
    auto Ys = getRandomSamplesY(indices, ns, Y);
   
    zerork_reactor_solve(step, t, dt, M, Ts.data(), Ps.data(), Ys.data(), zrm_handle);
    t += dt;
    
    // update solution back into the main arrays
    update_solution_vectors(indices, ns, T, P, Y, Ts, Ps, Ys);
    
    auto const [min_it, max_it] = std::minmax_element(T.begin(), T.end());

    std::cout<<"step = "<<step<<", t = " <<t<<"s : M = "<<M<<", Tmin = " << *min_it << ", Tmax = " << *max_it << std::endl;
  }

  MPI_Finalize();

  std::cout<<"End of computation!"<<std::endl;
  
  return 1;
}
