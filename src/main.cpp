
#include <chrono>

#include "main.h"
#include "Data.h"
#include "Parameters.h"
#include "Spatial_prior.h"
#include "Lookup.h"
#include "MCMC.h"
#include "misc.h"
//#include "probability.h"
//#include "hungarian.h"

using namespace std;

//------------------------------------------------
// run main MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_data = args["args_data"];
  Rcpp::List args_model = args["args_model"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  
  // read in data into separate class
  Data data(args_data);
  
  // read in parameters into separate class
  Parameters params(args_model);
  
  // define spatial prior
  Spatial_prior spatprior(args_model, params);
  
  // define lookup table from data and parameters. The lookup table 
  // precalculates the distance in km of every data point from every possible 
  // source location in a lon/lat grid. The grid is defined from the lon/lat min
  // and max values ()in the parameter set), and the lon/lat precision (in the
  // MCMC arguments).
  Lookup lookup(data, params);
  lookup.recalc();
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create MCMC object
  MCMC mcmc(data, params, lookup, spatprior);
  
  // run MCMC
  mcmc.burnin_mcmc(args_functions, args_progress);
  mcmc.sampling_mcmc(args_functions, args_progress);
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!params.silent) {
    print("   completed in", time_span.count(), "seconds\n");
  }
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( mcmc.loglike_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.source_lon_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.source_lat_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.source_realised_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.sigma_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.ep_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.alpha_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.loglike_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.source_lon_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.source_lat_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.source_realised_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.sigma_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.ep_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.alpha_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.qmatrix_final ));
  ret.push_back(Rcpp::wrap( mcmc.source_accept_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.source_accept_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.sigma_accept_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.sigma_accept_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.ep_accept_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.ep_accept_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.alpha_accept_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.alpha_accept_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.coupling_accept_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.coupling_accept_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.rung_converged ));
  ret.push_back(Rcpp::wrap( mcmc.convergence_iteration ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("loglike_burnin");
  ret_names.push_back("source_lon_burnin");
  ret_names.push_back("source_lat_burnin");
  ret_names.push_back("source_realised_burnin");
  ret_names.push_back("sigma_burnin");
  ret_names.push_back("ep_burnin");
  ret_names.push_back("alpha_burnin");
  ret_names.push_back("loglike_sampling");
  ret_names.push_back("source_lon_sampling");
  ret_names.push_back("source_lat_sampling");
  ret_names.push_back("source_realised_sampling");
  ret_names.push_back("sigma_sampling");
  ret_names.push_back("ep_sampling");
  ret_names.push_back("alpha_sampling");
  ret_names.push_back("qmatrix");
  ret_names.push_back("source_accept_burnin");
  ret_names.push_back("source_accept_sampling");
  ret_names.push_back("sigma_accept_burnin");
  ret_names.push_back("sigma_accept_sampling");
  ret_names.push_back("ep_accept_burnin");
  ret_names.push_back("ep_accept_sampling");
  ret_names.push_back("alpha_accept_burnin");
  ret_names.push_back("alpha_accept_sampling");
  ret_names.push_back("coupling_accept_burnin");
  ret_names.push_back("coupling_accept_sampling");
  ret_names.push_back("rung_converged");
  ret_names.push_back("convergence_iteration");
  
  ret.names() = ret_names;
  return ret;
}
