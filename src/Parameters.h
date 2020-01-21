
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
  std::vector<double> source_init;
  int sigma_model;
  double sigma_prior_meanlog;
  double sigma_prior_sdlog;
  double expected_popsize_prior_mean;
  double expected_popsize_prior_sd;
  double expected_popsize_prior_shape;
  double expected_popsize_prior_rate;
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
};
