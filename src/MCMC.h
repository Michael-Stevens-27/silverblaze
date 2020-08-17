
#pragma once

#include <Rcpp.h>
#include "Data.h"
#include "Parameters.h"
#include "Particle.h"

//------------------------------------------------
// class defining MCMC
class MCMC {
  
public:
  
  // pointers to data and parameters
  Parameters * p;
  Data * d;
  
  // thermodynamic parameters
  std::vector<int> rung_order;
  int cold_rung;
  
  // vector of particles
  std::vector<Particle> particle_vec;
  
  // ordering of labels
  std::vector<int> label_order;
  
  // qmatrices
  std::vector<std::vector<double>> log_qmatrix_running;
  std::vector<std::vector<double>> qmatrix_final;
  
  // objects for storing results
  std::vector<std::vector<double>> loglike_burnin;
  std::vector<std::vector<double>> source_lon_burnin;
  std::vector<std::vector<double>> source_lat_burnin;
  std::vector<std::vector<bool>> source_realised_burnin;
  std::vector<std::vector<double>> sigma_burnin;
  std::vector<std::vector<double>> ep_burnin;
  std::vector<double> alpha_burnin;
  
  std::vector<std::vector<double>> loglike_sampling;
  std::vector<std::vector<double>> source_lon_sampling;
  std::vector<std::vector<double>> source_lat_sampling;
  std::vector<std::vector<bool>> source_realised_sampling;
  std::vector<std::vector<double>> sigma_sampling;
  std::vector<std::vector<double>> ep_sampling;
  std::vector<double> alpha_sampling;
  
  // objects for storing acceptance rates
  std::vector<std::vector<int>> source_accept_burnin;
  std::vector<std::vector<int>> source_accept_sampling;
  std::vector<std::vector<int>> sigma_accept_burnin;
  std::vector<std::vector<int>> sigma_accept_sampling;
  std::vector<std::vector<int>> ep_accept_burnin;
  std::vector<std::vector<int>> ep_accept_sampling;
  std::vector<int> alpha_accept_burnin;
  std::vector<int> alpha_accept_sampling;
  std::vector<int> coupling_accept_burnin;
  std::vector<int> coupling_accept_sampling;
  
  // store convergence
  std::vector<bool> rung_converged;
  int convergence_iteration;  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC() {};
  MCMC(Data &data, Parameters &params, Lookup &lookup, Spatial_prior &spatprior);
  
  // other functions
  void burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void metropolis_coupling(bool burnin_phase);
  void sample_realised_sources(std::vector<std::vector<double>> &qmatrix, std::vector<int> const &label_order, std::vector<bool> &source_realised);
  
};
