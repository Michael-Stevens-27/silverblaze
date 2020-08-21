
#include "MCMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor for MCMC class
MCMC::MCMC(Data &data, Parameters &params, Lookup &lookup, Spatial_prior &spatprior) {
  
  // store pointer to data and parameters
  p = &params;
  d = &data;
  
  // initialise rung order
  rung_order = seq_int(0, p->rungs - 1);
  cold_rung = rung_order[p->rungs - 1];
  
  // vector of particles
  particle_vec = vector<Particle>(p->rungs);
  for (int r = 0; r < p->rungs; r++) {
    particle_vec[r] = Particle(data, params, lookup, spatprior, p->beta_vec[r]);
  }
  
  // initialise ordering of labels
  label_order = seq_int(0, p->K - 1);
  
  // qmatrices
  log_qmatrix_running = vector<vector<double>>(d->n, vector<double>(p->K));
  if (p->K == 1) {
    qmatrix_final = vector<vector<double>>(d->n, vector<double>(p->K, p->samples));
  } else {
    qmatrix_final = vector<vector<double>>(d->n, vector<double>(p->K));
  }
  
  // objects for storing results
  loglike_burnin = vector<vector<double>>(p->rungs, vector<double>(p->burnin));
  source_lon_burnin = vector<vector<double>>(p->burnin, vector<double>(p->K));
  source_lat_burnin = vector<vector<double>>(p->burnin, vector<double>(p->K));
  source_realised_burnin = vector<vector<bool>>(p->burnin, vector<bool>(p->K, false));
  sigma_burnin = vector<vector<double>>(p->burnin, vector<double>(p->K));
  ep_burnin =  vector<vector<double>>(p->burnin, vector<double>(p->K));
  alpha_burnin = vector<double>(p->burnin);
  
  loglike_sampling = vector<vector<double>>(p->rungs, vector<double>(p->samples));
  source_lon_sampling = vector<vector<double>>(p->samples, vector<double>(p->K));
  source_lat_sampling = vector<vector<double>>(p->samples, vector<double>(p->K));
  source_realised_sampling = vector<vector<bool>>(p->samples, vector<bool>(p->K, false));
  sigma_sampling = vector<vector<double>>(p->samples, vector<double>(p->K));
  ep_sampling = vector<vector<double>>(p->samples, vector<double>(p->K));
  alpha_sampling = vector<double>(p->samples);
  
  // objects for storing acceptance rates
  source_accept_burnin = vector<vector<int>>(p->rungs, vector<int>(p->K));
  source_accept_sampling = vector<vector<int>>(p->rungs, vector<int>(p->K));
  sigma_accept_burnin = vector<vector<int>>(p->rungs, vector<int>(p->K));
  sigma_accept_sampling = vector<vector<int>>(p->rungs, vector<int>(p->K));
  ep_accept_burnin = vector<vector<int>>(p->rungs, vector<int>(p->K));
  ep_accept_sampling = vector<vector<int>>(p->rungs, vector<int>(p->K));
  alpha_accept_burnin = vector<int>(p->rungs);
  alpha_accept_sampling = vector<int>(p->rungs);
  
  coupling_accept_burnin = vector<int>(p->rungs - 1);
  coupling_accept_sampling = vector<int>(p->rungs - 1);
  
  // store convergence
  rung_converged = vector<bool>(p->rungs, false);
  convergence_iteration = p->burnin;
  
}

//------------------------------------------------
// run burn-in phase of MCMC
void MCMC::burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // print header
  if (!p->silent) {
    print("Running MCMC for K =", p->K);
    print("Burn-in phase");
  }
  
  // read in R functions
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // reset particles
  for (int r = 0; r < p->rungs; r++) {
    particle_vec[r].reset(p->beta_vec[r]);
  }
  rung_order = seq_int(0, p->rungs-1);
  
  // loop through burnin iterations
  bool all_convergence_reached = false;
  for (int rep = 0; rep < p->burnin; rep++) {
    
    // update particles
    for (int r = 0; r < p->rungs; r++) {
      int rung = rung_order[r];
      
      // update sources
      particle_vec[rung].update_sources(true, rep + 1);
      
      // update sigma
      particle_vec[rung].update_sigma(true, rep + 1);
      
      // update expected population size
      particle_vec[rung].update_expected_popsize(true, rep + 1);
      
      // update alpha
      particle_vec[rung].update_alpha(true, rep + 1);
      
    } // end loop over rungs
    
    // apply Metropolis-coupling
    if (p->coupling_on) {
      metropolis_coupling(true);
    }
    
    // focus on cold rung
    cold_rung = rung_order[p->rungs - 1];
    
    // methods that only apply when K>1
    if (p->K > 1) {
      
      // update qmatrix of cold rung
      particle_vec[cold_rung].update_qmatrix();
      
      // fix labels
      particle_vec[cold_rung].solve_label_switching(log_qmatrix_running);
      label_order = particle_vec[cold_rung].label_order;
      
      // add particle log_qmatrix to log_qmatrix_running
      for (int i = 0; i < d->n; ++i) {
        for (int k = 0; k < p->K; ++k) {
          log_qmatrix_running[i][k] = log_sum(log_qmatrix_running[i][k], particle_vec[cold_rung].log_qmatrix[i][label_order[k]]);
        }
      }
      
      // draw realised sources
      sample_realised_sources(particle_vec[cold_rung].qmatrix, label_order, source_realised_burnin[rep]);
      
    }
    
    // store loglikelihood
    for (int r = 0; r < p->rungs; ++r) {
      int rung = rung_order[r];
      loglike_burnin[r][rep] = particle_vec[rung].loglike;
    }
    
    // store source locations
    for (int k = 0; k < p->K; ++k) {
      source_lon_burnin[rep][k] = particle_vec[cold_rung].source_lon[label_order[k]];
      source_lat_burnin[rep][k] = particle_vec[cold_rung].source_lat[label_order[k]];
      
      // store sigma
      sigma_burnin[rep][k] = particle_vec[cold_rung].sigma[label_order[k]];
      
      // store expected population size or source weights
      if (d->data_type == 3){
        ep_burnin[rep][k] = particle_vec[cold_rung].source_weights[label_order[k]];
      } else {
        ep_burnin[rep][k] = particle_vec[cold_rung].expected_popsize[label_order[k]];
      }
      
      // store alpha
      alpha_burnin[rep] = particle_vec[cold_rung].alpha;
    }
    
    // update progress bars
    if (!p->silent) {
      if ((rep+1) == p->burnin) {
        update_progress(args_progress, "pb_burnin", rep + 1, p->burnin);
      } else {
        int remainder = rep % int(ceil(double(p->burnin)/100));
        if (remainder == 0 && !p->pb_markdown) {
          update_progress(args_progress, "pb_burnin", rep + 1, p->burnin);
        }
      }
    }
    
    // check for convergence
    if ((p->auto_converge && ((rep + 1) % p->converge_test) == 0) || (rep + 1) == p->burnin) {
      
      // check for convergence of each chain
      for (int r = 0; r < p->rungs; ++r) {
        rung_converged[r] = rcpp_to_bool(test_convergence(loglike_burnin[r], rep + 1));
      }
      
      // break if convergence reached
      all_convergence_reached = true;
      for (int r = 0; r < p->rungs; ++r) {
        if (!rung_converged[r]) {
          all_convergence_reached = false;
          break;
        }
      }
      
      // end if all reached convergence
      if (all_convergence_reached) {
        convergence_iteration = rep + 1;
        if (!p->silent) {
          update_progress(args_progress, "pb_burnin", p->burnin, p->burnin);
          print("   converged within", convergence_iteration, "iterations");
        }
        break;
      }
      
    }  // end if auto_converge
    
  }  // end burn-in iterations
  
  // store acceptance rates
  for (int r = 0; r < p->rungs; ++r) {
    int rung = rung_order[r];
    source_accept_burnin[r] = particle_vec[rung].source_accept_burnin;
    sigma_accept_burnin[r] = particle_vec[rung].sigma_accept_burnin;
    ep_accept_burnin[r] = particle_vec[rung].ep_accept_burnin;
    alpha_accept_burnin[r] = particle_vec[rung].alpha_accept_burnin;
  }
  
  // warning if still not converged
  if (!all_convergence_reached && !p->silent) {
    print("   Warning: convergence still not reached within", p->burnin, "iterations");
  }
  
}

//------------------------------------------------
// run sampling phase of MCMC
void MCMC::sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // print header
  if (!p->silent) {
    print("Sampling phase");
  }
  
  // read in R functions
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // loop through sampling iterations
  for (int rep = 0; rep < p->samples; ++rep) {
    
    // update particles
    for (int r = 0; r < p->rungs; ++r) {
      int rung = rung_order[r];
      
      // update sources
      particle_vec[rung].update_sources(false, 0);
      
      // update sigma
      particle_vec[rung].update_sigma(false, 0);
      
      // update expected popsize
      particle_vec[rung].update_expected_popsize(false, 0);
      
      // update alpha
      particle_vec[rung].update_alpha(false, 0);
            
    } // end loop over rungs
    
    // apply Metropolis-coupling
    if (p->coupling_on) {
      metropolis_coupling(false);
    }
    
    // focus on cold rung
    cold_rung = rung_order[p->rungs - 1];
    
    // methods that only apply when K>1
    if (p->K > 1) {
      
      // update qmatrix of cold rung
      particle_vec[cold_rung].update_qmatrix();
      
      // fix labels
      particle_vec[cold_rung].solve_label_switching(log_qmatrix_running);
      label_order = particle_vec[cold_rung].label_order;
      
      // add particle log_qmatrix to log_qmatrix_running
      for (int i = 0; i < d->n; ++i) {
        for (int k = 0; k < p->K; ++k) {
          log_qmatrix_running[i][k] = log_sum(log_qmatrix_running[i][k], particle_vec[cold_rung].log_qmatrix[i][label_order[k]]);
        }
      }
      
      // add particle qmatrix to qmatrix_final
      for (int i = 0; i < d->n; ++i) {
        for (int k = 0; k < p->K; ++k) {
          qmatrix_final[i][k] += particle_vec[cold_rung].qmatrix[i][label_order[k]];
        }
      }
      
      // draw realised sources
      sample_realised_sources(particle_vec[cold_rung].qmatrix, label_order, source_realised_sampling[rep]);
      
    }
    
    // store loglikelihood
    for (int r = 0; r < p->rungs; ++r) {
      int rung = rung_order[r];
      loglike_sampling[r][rep] = particle_vec[rung].loglike;
    }
    
    // store source locations
    for (int k = 0; k < p->K; ++k) {
      source_lon_sampling[rep][k] = particle_vec[cold_rung].source_lon[label_order[k]];
      source_lat_sampling[rep][k] = particle_vec[cold_rung].source_lat[label_order[k]];
      
      // store sigma
      sigma_sampling[rep][k] = particle_vec[cold_rung].sigma[label_order[k]];
      
      // store expected population size or source weights
      if (d->data_type == 3){
        ep_sampling[rep][k] = particle_vec[cold_rung].source_weights[label_order[k]];
      } else {
        ep_sampling[rep][k] = particle_vec[cold_rung].expected_popsize[label_order[k]];
      }
            
      // store alpha
      alpha_sampling[rep] = particle_vec[cold_rung].alpha;
    }
    
    // update progress bars
    if (!p->silent) {
      if ((rep+1) == p->samples) {
        update_progress(args_progress, "pb_samples", rep+1, p->samples);
      } else {
        int remainder = rep % int(ceil(double(p->samples)/100));
        if (remainder == 0 && !p->pb_markdown) {
          update_progress(args_progress, "pb_samples", rep+1, p->samples);
        }
      }
    }
    
  } // end sampling iterations
  
  // store acceptance rates
  for (int r = 0; r < p->rungs; ++r) {
    int rung = rung_order[r];
    
    source_accept_sampling[r] = particle_vec[rung].source_accept_sampling;
    sigma_accept_sampling[r] = particle_vec[rung].sigma_accept_sampling;
    ep_accept_sampling[r] = particle_vec[rung].ep_accept_sampling;
    alpha_accept_sampling[r] = particle_vec[rung].alpha_accept_sampling;
  }

}

//------------------------------------------------
// Metropolis-coupling to propose swaps between temperature rungs
void MCMC::metropolis_coupling(bool burnin_phase) {
  
  // loop over rungs, starting with the hottest chain and moving to the cold
  // chain. Each time propose a swap with the next rung up
  for (int i = 0; i < (p->rungs - 1); i++) {
    
    // define rungs of interest
    int rung1 = rung_order[i];
    int rung2 = rung_order[i + 1];
    
    // get log-likelihoods and beta values of two chains in the comparison
    double loglike1 = particle_vec[rung1].loglike;
    double loglike2 = particle_vec[rung2].loglike;
    
    double beta1 = particle_vec[rung1].beta;
    double beta2 = particle_vec[rung2].beta;
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (loglike2*beta1 + loglike1*beta2) - (loglike1*beta1 + loglike2*beta2);
    
    // accept or reject move
    double rand1 = runif1();
    if (log(rand1) < acceptance) {
      
      // swap beta values
      particle_vec[rung1].beta = beta2;
      particle_vec[rung2].beta = beta1;
      
      // swap source proposal SD
      vector<double> store_source_propSD1 = particle_vec[rung1].source_propSD;
      particle_vec[rung1].source_propSD = particle_vec[rung2].source_propSD;
      particle_vec[rung2].source_propSD = store_source_propSD1;
      
      // swap sigma proposal SD
      vector<double> store_sigma_propSD1 = particle_vec[rung1].sigma_propSD;
      particle_vec[rung1].sigma_propSD = particle_vec[rung2].sigma_propSD;
      particle_vec[rung2].sigma_propSD = store_sigma_propSD1;
      
      // swap expected_pop proposal SD
      vector<double> store_ep_propSD1 = particle_vec[rung1].ep_propSD;
      particle_vec[rung1].ep_propSD = particle_vec[rung2].ep_propSD;
      particle_vec[rung2].ep_propSD = store_ep_propSD1;
      
      // swap alpha proposal SD
      double store_alpha_propSD1 = particle_vec[rung1].alpha_propSD;
      particle_vec[rung1].alpha_propSD = particle_vec[rung2].alpha_propSD;
      particle_vec[rung2].alpha_propSD = store_alpha_propSD1;
      
      // swap rung order
      rung_order[i] = rung2;
      rung_order[i + 1] = rung1;
      
      if (burnin_phase){
        // update burnin coupling acceptance rates
        coupling_accept_burnin[i]++;
      } else {
        // update sampling coupling acceptance rates
        coupling_accept_sampling[i]++;
      }
    }
  }

}

//------------------------------------------------
// draw from the posterior allocation of observations to sources
void MCMC::sample_realised_sources(vector<vector<double>> &qmatrix, vector<int> const &label_order, vector<bool> &source_realised) {
  
  // currently only applies to prevalence model
  if (d->data_type != 2) {
    return;
  }
  
  // sample allocation for each row of qmatrix
  vector<int> rand_multinom(p->K);
  for (int i = 0; i < d->n; ++i) {
    
    // skip if no observation at this site
    if (d->positive[i] == 0) {
      continue;
    }
    
    // sample a source for each positive observation
    rmultinom1(d->positive[i], qmatrix[i], 1.0, rand_multinom);
    
    // source is realised if any observation is allocated to ii
    for (int j = 0; j < p->K; ++j) {
      if (rand_multinom[label_order[j]] > 0) {
        source_realised[j] = true;
      }
    }
    
  }  // end i loop
  
}
