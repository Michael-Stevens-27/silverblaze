
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing all parameters
class Parameters {
  
public:
  
  // MCMC parameters
  int burnin;
  int samples;
  std::vector<double> beta_vec;
  int rungs;
  bool auto_converge;
  int converge_test;
  bool pb_markdown;
  bool silent;
  bool coupling_on;
  
  // model parameters
  int K;
  double sentinel_radius;
  double min_lon;
  double max_lon;
  double res_lon;
  int n_lon;
  double min_lat;
  double max_lat;
  double res_lat;
  int n_lat;
  double search_area;
  double cell_area;
  std::vector<double> source_init_lon;
  std::vector<double> source_init_lat;
  int sigma_model;
  double sigma_prior_meanlog;
  double sigma_prior_sdlog;
  int ep_model;
  int dispersal_model;
  double ep_prior_mean;
  double ep_prior_sd;
  double ep_prior_shape;
  double ep_prior_rate;
  double ep_prior_meanlog;
  double ep_prior_sdlog;
  int count_type;
  double alpha_prior_meanlog;
  double alpha_prior_sdlog;
  
  bool bugged;
  int chosen_rung;
  double dirichlet_scale;
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
};
