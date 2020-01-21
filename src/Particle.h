
#pragma once

#include <Rcpp.h>

#include "Lookup.h"
#include "Spatial_prior.h"
#include "Data.h"
#include "Parameters.h"

//------------------------------------------------
// class defining particle
class Particle {

public:
  // PUBLIC OBJECTS
  
  // pointers
  Data * d;
  Parameters * p;
  Lookup * l;
  Spatial_prior * sp;
  
  // value of the thermodynamic power
  double beta;
  
  // source locations
  std::vector<double> source_lon;
  std::vector<double> source_lat;
  
  // standard deviation of sources (km)
  std::vector<double> sigma;
  
  // scaling factor on hazard surface
  double expected_popsize;
  double log_expected_popsize;
  
  // qmatrices
  std::vector<std::vector<double>> log_qmatrix;
  std::vector<std::vector<double>> qmatrix;
  
  // proposal standard deviations
  std::vector<double> source_propSD;
  std::vector<double> sigma_propSD;
  
  // Robbins-Monro stepsize
  double source_rm_stepsize;
  double sigma_rm_stepsize;
  
  // misc constants
  double log_sentinel_area;
  int counts_total;
  double log_K;
  
  // likelihood
  std::vector<std::vector<double>> dist_source_data;
  std::vector<double> dist_source_data_prop;
  std::vector<std::vector<double>> log_hazard_height;
  std::vector<double> log_hazard_height_prop;
  std::vector<std::vector<double>> log_hazard_height_prop2;
  double logprior;
  double loglike;
  
  // initialise ordering of labels
  std::vector<int> label_order;
  std::vector<int> label_order_new;
  
  // objects for solving label switching problem
  std::vector<std::vector<double>> cost_mat;
  std::vector<int> best_perm;
  std::vector<int> best_perm_order;
  std::vector<int> edges_left;
  std::vector<int> edges_right;
  std::vector<int> blocked_left;
  std::vector<int> blocked_right;
  
  // store acceptance rates
  std::vector<int> source_accept_burnin;
  std::vector<int> source_accept_sampling;
  std::vector<int> sigma_accept_burnin;
  std::vector<int> sigma_accept_sampling;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  Particle(Data &data, Parameters &params, Lookup &lookup, Spatial_prior &spatprior, double beta);
  
  // other functions
  void reset(double beta);
  double calculate_logprior_source(double source_lon_prop, double source_lat_prop);
  double calculate_loglike_source(double source_lon_prop, double source_lat_prop, int k);
  void update_sources(bool robbins_monro_on, int iteration);
  void update_sigma(bool robbins_monro_on, int iteration);
  void update_sigma_single(bool robbins_monro_on, int iteration);
  void update_sigma_independent(bool robbins_monro_on, int iteration);
  void update_expected_popsize();
  void update_qmatrix();
  void solve_label_switching(const std::vector<std::vector<double>> &log_qmatrix_running);
  
};
