
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
  std::vector<double> expected_popsize;
  std::vector<double> pop_size_domain;
  std::vector<double> cum_sum_density; 
  std::vector<double> cum_sum_normalised;
  
  // parameter controlling the nbinom variance = mean + alpha*mean^2 
  double alpha; 
  
  // qmatrices
  std::vector<std::vector<double>> log_qmatrix;
  std::vector<std::vector<double>> qmatrix;
  
  // proposal standard deviations
  std::vector<double> source_propSD;
  std::vector<double> sigma_propSD;
  std::vector<double> ep_propSD;
  double alpha_propSD;  
    
  // Robbins-Monro stepsize
  double source_rm_stepsize;
  double sigma_rm_stepsize;
  double ep_rm_stepsize;
  double alpha_rm_stepsize;
  
  // misc constants
  double log_sentinel_area;
  int counts_total;
  int tested_total;
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
  std::vector<int> ep_accept_burnin;
  std::vector<int> ep_accept_sampling;
  int alpha_accept_burnin;
  int alpha_accept_sampling;
    
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  Particle(Data &data, Parameters &params, Lookup &lookup, Spatial_prior &spatprior, double beta);
  
  // other functions
  void reset(double beta);
  
  // main update switches
  void update_sources(bool robbins_monro_on, int iteration);
  void update_sigma(bool robbins_monro_on, int iteration);
  void update_expected_popsize(bool robbins_monro_on, int iteration);
  void update_alpha(bool robbins_monro_on, int iteration);
  
  // calculate hazrad based on dispersal kernel
  double calculate_hazard(double dist, double single_scale);
  
  // switch for sources
  double calculate_loglike_source(double source_lon_prop, double source_lat_prop, int k);
  
  // loglikelihood functions for sources
  double calculate_loglike_source_pois(double source_lon_prop, double source_lat_prop, int k);
  double calculate_loglike_source_ind_exp_pop(double source_lon_prop, double source_lat_prop, int k);
  double calculate_loglike_source_binom(double source_lon_prop, double source_lat_prop, int k);
  double calculate_loglike_source_points(double source_lon_prop, double source_lat_prop, int k);
  double calculate_loglike_source_negative_binomial_indpendent_lambda(double source_lon_prop, double source_lat_prop, int k);
  
  double calculate_logprior_source(double source_lon_prop, double source_lat_prop);
    
  // loglikelihood and update functions for sigmas  
  void update_sigma_pois(bool robbins_monro_on, int iteration);
  void update_sigma_pois_ind_exp_pop(bool robbins_monro_on, int iteration);
  void update_sigma_binom(bool robbins_monro_on, int iteration);
  void update_sigma_points(bool robbins_monro_on, int iteration); 
  void update_sigma_negative_binomial_ind_exp_pop(bool robbins_monro_on, int iteration);

  // loglikelihood and update functions for expected population size
  void update_expected_popsize_pois_single();
  void update_expected_popsize_pois_independent(bool robbins_monro_on, int iteration);
  // void update_expected_popsize_pois_independent(bool robbins_monro_on);
  void update_expected_popsize_binom(bool robbins_monro_on, int iteration);
  void update_expected_popsize_negative_binomial_independent(bool robbins_monro_on, int iteration);
  
  // loglikelihood and update for alpha parameter controlling the nbinom variance
  void update_alpha_negative_binomial(bool robbins_monro_on, int iteration);

  // misc  
  void update_qmatrix();
  void solve_label_switching(const std::vector<std::vector<double>> &log_qmatrix_running);
  
};
