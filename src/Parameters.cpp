
#include "Parameters.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// constructor for Parameters class
Parameters::Parameters(const Rcpp::List &args) {
  
  // MCMC parameters
  burnin = rcpp_to_int(args["burnin"]);
  samples = rcpp_to_int(args["samples"]);
  beta_vec = rcpp_to_vector_double(args["beta_vec"]);
  rungs = int(beta_vec.size());
  auto_converge = rcpp_to_bool(args["auto_converge"]);
  converge_test = rcpp_to_int(args["converge_test"]);
  pb_markdown = rcpp_to_bool(args["pb_markdown"]);
  silent = rcpp_to_bool(args["silent"]);
  coupling_on = rcpp_to_bool(args["coupling_on"]);
  chosen_rung = rcpp_to_int(args["rung_store"]);
  
  // model parameters
  K = rcpp_to_int(args["K"]);
  sentinel_radius = rcpp_to_double(args["sentinel_radius"]);
  min_lon = rcpp_to_double(args["min_lon"]);
  max_lon = rcpp_to_double(args["max_lon"]);
  res_lon = rcpp_to_double(args["res_lon"]);
  n_lon = rcpp_to_int(args["n_lon"]);
  min_lat = rcpp_to_double(args["min_lat"]);
  max_lat = rcpp_to_double(args["max_lat"]);
  res_lat = rcpp_to_double(args["res_lat"]);
  n_lat = rcpp_to_int(args["n_lat"]);
  search_area = rcpp_to_double(args["study_area"]);
  cell_area = rcpp_to_double(args["cell_area"]);
  source_init_lon = rcpp_to_vector_double(args["source_init_lon"]);
  source_init_lat = rcpp_to_vector_double(args["source_init_lat"]);
  sigma_model = rcpp_to_int(args["sigma_model_numeric"]);
  ep_model = rcpp_to_int(args["expected_popsize_model_numeric"]);
  dispersal_model = rcpp_to_int(args["dispersal_model_numeric"]);
  count_type = rcpp_to_int(args["count_type_numeric"]);  // negative binomial = 1, non-negative binomial (currently Poisson) = 2
  
    // get sigma prior mean and sd in log space from raw inputs
  double sigma_prior_mean = rcpp_to_double(args["sigma_prior_mean"]);
  double sigma_prior_sd = rcpp_to_double(args["sigma_prior_sd"]);
  double sigma_prior_varlog = log(pow(sigma_prior_sd, 2) / pow(sigma_prior_mean, 2) + 1);
  sigma_prior_sdlog = sqrt(sigma_prior_varlog);
  sigma_prior_meanlog = log(sigma_prior_mean) - sigma_prior_varlog / 2.0;
  
  // get expected_popsize shape and rate parameters from raw inputs
  ep_prior_mean = rcpp_to_double(args["expected_popsize_prior_mean"]);
  ep_prior_sd = rcpp_to_double(args["expected_popsize_prior_sd"]);
  ep_prior_shape = pow(ep_prior_mean, 2) / pow(ep_prior_sd, 2);
  ep_prior_rate = ep_prior_mean / pow(ep_prior_sd, 2);
  double ep_prior_varlog = log(pow(ep_prior_sd, 2) / pow(ep_prior_mean, 2) + 1);
  ep_prior_sdlog = sqrt(ep_prior_varlog);
  ep_prior_meanlog = log(ep_prior_mean) - ep_prior_varlog / 2.0;
  
  // get alpha prior mean and sd in log space from raw inputs
  double alpha_prior_mean = rcpp_to_double(args["alpha_prior_mean"]);
  double alpha_prior_sd = rcpp_to_double(args["alpha_prior_sd"]);
  double alpha_prior_varlog = log(pow(alpha_prior_sd, 2) / pow(alpha_prior_mean, 2) + 1);
  alpha_prior_sdlog = sqrt(alpha_prior_varlog);
  alpha_prior_meanlog = log(alpha_prior_mean) - alpha_prior_varlog / 2.0;
  
  // get parameter for dirichlet prior on weights
  weight_prior = rcpp_to_double(args["weight_prior"]);
}
