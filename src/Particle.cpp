
#include "Particle.h"
#include "misc.h"
#include "probability.h"
#include "hungarian.h"

using namespace std;

//------------------------------------------------
// constructor for Particle class
Particle::Particle(Data &data, Parameters &params, Lookup &lookup, Spatial_prior &spatprior, double beta) {
  
  // load pointers
  d = &data;
  p = &params;
  l = &lookup;
  sp = &spatprior;
  
  // value of the thermodynamic power
  this -> beta = beta;
  
  // source locations
  source_lon = vector<double>(p->K);
  source_lat = vector<double>(p->K);
  
  // standard deviation of sources (km)
  sigma = vector<double>(p->K, 1);
  
  // expected population size for each soure
  expected_popsize = vector<double>(p->K, 100);
  ep_total = 0;
    
  // weights for each soure
  source_weights = vector<double>(p->K, 1/(p->K));
  source_weight_prop = vector<double>(p->K, 1/(p->K));
  
  // alpha param for nbinom variance
  alpha = 1;
  
  // qmatrices
  log_qmatrix = vector<vector<double>>(d->n, vector<double>(p->K));
  qmatrix = vector<vector<double>>(d->n, vector<double>(p->K));
  
  // proposal standard deviations
  source_propSD = vector<double>(p->K, 0.01);
  sigma_propSD = vector<double>(p->K, 1.0);
  ep_propSD = vector<double>(p->K, 100.0);
  alpha_propSD = 1;
  
  // Robbins-Monro stepsize constants
  source_rm_stepsize = 5.0;
  sigma_rm_stepsize = 4.0;
  ep_rm_stepsize = 4.0;
  alpha_rm_stepsize = 4.0;
  
  // misc constants
  // area around a sentinel site
  log_sentinel_area = LOG_PI + 2*log(p->sentinel_radius);
  
  // sum of counts over all sentinel sites
  counts_total = sum(d->counts);
  tested_total = sum(d->tested);
    
  // log of K
  log_K = log(p->K);
  
  // likelihood
  dist_source_data = vector<vector<double>>(d->n, vector<double>(p->K));
  dist_source_data_prop = vector<double>(d->n);
  log_hazard_height = vector<vector<double>>(d->n, vector<double>(p->K));
  log_hazard_height_prop = vector<double>(d->n);
  log_hazard_height_prop2 = vector<vector<double>>(d->n, vector<double>(p->K));
  logprior = 0;
  loglike = 0;
  
  // initialise ordering of labels
  label_order = seq_int(0, p->K - 1);
  label_order_new = vector<int>(p->K);
  
  // objects for solving label switching problem
  cost_mat = vector<vector<double>>(p->K, vector<double>(p->K));
  best_perm = vector<int>(p->K);
  best_perm_order = vector<int>(p->K);
  edges_left = vector<int>(p->K);
  edges_right = vector<int>(p->K);
  blocked_left = vector<int>(p->K);
  blocked_right = vector<int>(p->K);
  
  // store acceptance rates
  source_accept_burnin = vector<int>(p->K);
  source_accept_sampling = vector<int>(p->K);
  sigma_accept_burnin = vector<int>(p->K);
  sigma_accept_sampling = vector<int>(p->K);
  ep_accept_burnin = vector<int>(p->K);
  ep_accept_sampling = vector<int>(p->K);
  alpha_accept_burnin = 0;
  alpha_accept_sampling = 0;
}

//------------------------------------------------
// reset particle
void Particle::reset(double beta) {
  
  // reset beta value
  this->beta = beta;
  
  // initialise source locations
  for (int k = 0; k < p->K; ++k) {
    // source_lon[k] = p->source_init_lon[k];
    // source_lat[k] = p->source_init_lat[k];
    
    // randomly draw from spatial extent of prior
    source_lon[k] = runif1(p->min_lon, p->max_lon);
    source_lat[k] = runif1(p->min_lat, p->max_lat);
  }
  
  // draw sigma from prior
  if (p->sigma_prior_sdlog == 0) {
    sigma = vector<double>(p->K, exp(p->sigma_prior_meanlog));
  } else {
    if (p->sigma_model == 1) {
      sigma = vector<double>(p->K, exp(rnorm1(p->sigma_prior_meanlog, p->sigma_prior_sdlog)) );
    } else if (p->sigma_model == 2) {
      for (int k = 0; k < p->K; ++k) {
        sigma[k] = exp(rnorm1(p->sigma_prior_meanlog, p->sigma_prior_sdlog));
      }
    }
  }
  
  // draw expected popsize from prior
  if (p->ep_prior_sd <= 0) {
    expected_popsize = vector<double>(p->K, p->ep_prior_mean);
  } else {
    if(d->data_type == 1){
      if(p->ep_model == 1){
        expected_popsize = vector<double>(p->K, rgamma1(p->ep_prior_shape, p->ep_prior_rate));
      } else if(p->ep_model == 2){
        for (int k = 0; k < p->K; ++k) {
          expected_popsize[k] = rgamma1(p->ep_prior_shape, p->ep_prior_rate);
        }
      }
    } else if(d->data_type == 2){
      for (int k = 0; k < p->K; ++k) {
        expected_popsize[k] = rgamma1(p->ep_prior_shape / double(p->K), p->ep_prior_rate);
      } 
    } else if(d->data_type == 3){
      for (int k = 0; k < p->K; ++k) {
        expected_popsize[k] = rgamma1(p->ep_prior_shape / double(p->K), p->ep_prior_rate);
        ep_total += expected_popsize[k];
      }
      for (int k = 0; k < p->K; ++k) { // build source weights based on expected pop size
        source_weights[k] = expected_popsize[k]/ep_total;
      }
    }
  }
  
  // draw alpha (negative binomial parameter) from prior
  alpha = exp(rnorm1(p->alpha_prior_meanlog, p->alpha_prior_sdlog));
  
  // initialise proposal standard deviations
  source_propSD = vector<double>(p->K, 0.01);
  sigma_propSD = vector<double>(p->K, 1.0);
  ep_propSD = vector<double>(p->K, 1.0);
  alpha_propSD = 1; 
  
  // calculate initial likelihood. Calling calculate_loglike_source() on each
  // source updates the dist_source_data_prop and log_hazard_height_prop
  // objects, which can then be stored back into the final matrices. This is
  // equivalent to running a Metropolis-Hastings step in which the move is
  // guaranteed to be accepted
  for (int k = 0; k < p->K; ++k) {
    logprior = calculate_logprior_source(source_lon[k], source_lat[k]);
    loglike = calculate_loglike_source(source_lon[k], source_lat[k], k);
    
    for (int i = 0; i < d->n; ++i) {
      dist_source_data[i][k] = dist_source_data_prop[i];
      log_hazard_height[i][k] = log_hazard_height_prop[i];
    }
  }
  log_hazard_height_prop2 = log_hazard_height; 
  
  // reset acceptance rates
  source_accept_burnin = vector<int>(p->K);
  source_accept_sampling = vector<int>(p->K);
  sigma_accept_burnin = vector<int>(p->K);
  sigma_accept_sampling = vector<int>(p->K);
  ep_accept_burnin = vector<int>(p->K);
  ep_accept_sampling = vector<int>(p->K);
  alpha_accept_burnin = 0;
  alpha_accept_sampling = 0;
}

//------------------------------------------------
// calculate log-likelihood given new proposed source
double Particle::calculate_loglike_source(double source_lon_prop, double source_lat_prop, int k) {
  
  // update source based on binomial or poisson model
  double ret = 0.0;
  if (d->data_type == 1) { // count data
    if(p->count_type == 1) { // negative binomial count data
      ret = calculate_loglike_source_negative_binomial_indpendent_lambda(source_lon_prop, source_lat_prop, k);
    } else {
      if(p->ep_model == 1){ // single expected pop size
      ret = calculate_loglike_source_pois(source_lon_prop, source_lat_prop, k);
    } else if(p->ep_model == 2){ // independent expected pop size
      ret = calculate_loglike_source_ind_exp_pop(source_lon_prop, source_lat_prop, k);
    } 
  }
  } else if (d->data_type == 2) { // prevalence data
    ret = calculate_loglike_source_binom(source_lon_prop, source_lat_prop, k);
  } else if (d->data_type == 3){ // point pattern data
    ret = calculate_loglike_source_points(source_lon_prop, source_lat_prop, k);
  }
  
  return ret;
}

//------------------------------------------------
// update sigma
void Particle::update_sigma(bool robbins_monro_on, int iteration) {
  
  // return if prior is exact
  if (p->sigma_prior_sdlog == 0) {
    return;
  }
  
  if (d->data_type == 1) { // count data
    if(p->count_type == 1) { // negative binomial count data
      update_sigma_negative_binomial_ind_exp_pop(robbins_monro_on, iteration);
    } else{
      if(p->ep_model == 1){ // single expected pop size
        update_sigma_pois(robbins_monro_on, iteration);
      } else if(p->ep_model == 2){ // independent expected pop size
        update_sigma_pois_ind_exp_pop(robbins_monro_on, iteration);
      } 
    }    
  } else if (d->data_type == 2) { // prevalence data
    update_sigma_binom(robbins_monro_on, iteration);
  } else if (d->data_type == 3){ // point pattern data
    update_sigma_points(robbins_monro_on, iteration);
  }
}

//------------------------------------------------
// update expected popsize
void Particle::update_expected_popsize(bool robbins_monro_on, int iteration) {
  
  // return if prior is exact
  // return if using point pattern data
  if (p->ep_prior_sd == 0) {
    return;
  } 
  
  // update expected_popsize based on binomial or poisson model
  if (d->data_type == 1) { // count data
    if (p->count_type == 1) { // negative binomial count data
      update_expected_popsize_negative_binomial_independent(robbins_monro_on, iteration);
    } else {
      if (p->ep_model == 1) { // single expected pop size
        update_expected_popsize_pois_single();
      } else if (p->ep_model == 2) { // independent expected pop size
        update_expected_popsize_pois_independent(robbins_monro_on, iteration);
      }
    }
  } else if (d->data_type == 2) { // prevalence data
    update_expected_popsize_binom(robbins_monro_on, iteration);
  } else if (d->data_type == 3) { // point-pattern data
    update_weights_point_pattern(robbins_monro_on, iteration);
  }  
}

//------------------------------------------------
// update alpha (negative binomial parameter)
void Particle::update_alpha(bool robbins_monro_on, int iteration) {
  
  // update if using count data and negative binomial model
  if (d->data_type == 1 && p->count_type == 1) {
    update_alpha_negative_binomial(robbins_monro_on, iteration);
  }
}

//------------------------------------------------
// calculate log-prior given new proposed source
double Particle::calculate_logprior_source(double source_lon_prop, double source_lat_prop) {
  
  // get logprior probability
  double logprior_prob = sp->get_value(source_lon_prop, source_lat_prop);
  
  // catch values with zero prior probability
  if (logprior_prob == 0) {
    logprior_prob = -OVERFLO;
  }
  
  return logprior_prob;
}

//------------------------------------------------
// calculate log-likelihood under Poisson model given new proposed source
double Particle::calculate_loglike_source_pois(double source_lon_prop, double source_lat_prop, int k) {
  
  // initialise running values
  double loglike_prop = 0;
  double theta_sum = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // get distance from proposed source to data point i
    double dist = l->get_data_dist(source_lon_prop, source_lat_prop, i);
    dist_source_data_prop[i] = dist;
    
    // calculate bivariate height of data point i from proposed source.
    log_hazard_height_prop[i] = calculate_hazard(dist, sigma[k]); 
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop[i];
    for (int j = 0; j < p->K; ++j) {
      if (j == k) {
        continue;
      }
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // define theta_i as the sentinel area * the mean hazard. Calculate
    // log(theta_i) and add theta_i to running sum
    double log_theta_i = log_sentinel_area + log_hazard_sum - log_K;
    theta_sum += exp(log_theta_i);
    
    // add necessary terms to loglikelihood
    loglike_prop += d->counts[i]*log_theta_i - lgamma(d->counts[i] + 1);
  }
  
  // complete loglikelihood
  if (p->ep_prior_sd == 0) {
    loglike_prop += counts_total*log(p->ep_prior_mean) - p->ep_prior_mean*theta_sum;
  } else {
    double gamma_shape = p->ep_prior_shape;
    double gamma_rate = p->ep_prior_rate;
    
    loglike_prop += gamma_shape*log(gamma_rate) - (gamma_shape + counts_total)*log(gamma_rate + theta_sum) + lgamma(gamma_shape + counts_total) - lgamma(gamma_shape);
  }
  
  return loglike_prop;
}

//------------------------------------------------
// calculate log-likelihood given new proposed source
double Particle::calculate_loglike_source_ind_exp_pop(double source_lon_prop, double source_lat_prop, int k) {
  
  // initialise new likelihood
  double loglike_prop = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // get distance from proposed source to data point i
    double dist = l->get_data_dist(source_lon_prop, source_lat_prop, i);
    dist_source_data_prop[i] = dist;
    
    // calculate bivariate height of data point i from proposed source.
    log_hazard_height_prop[i] =  log(expected_popsize[k]) + calculate_hazard(dist, sigma[k]);
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop[i];
    for (int j = 0; j < p->K; ++j) {
      if (j == k) {
        continue;
      }
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // define the rate lambda of the Poisson process at this sentinel site,
    // while remaining in log space
    double log_lambda = log_sentinel_area + log_hazard_sum;
        
    // calculate the Poisson log-probability of the counts at this sentinel site
    loglike_prop += d->counts[i]*log_lambda - exp(log_lambda) - lgamma(d->counts[i] + 1);
  }
  
  return loglike_prop;
}

//------------------------------------------------
// calculate log-likelihood under a binomial model given new proposed source
double Particle::calculate_loglike_source_binom(double source_lon_prop, double source_lat_prop, int k) {
  
  // initialise running values
  double loglike_prop = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // get distance from proposed source to data point i
    double dist = l->get_data_dist(source_lon_prop, source_lat_prop, i);
    dist_source_data_prop[i] = dist;
    
    // calculate bivariate height of data point i from proposed source.
    log_hazard_height_prop[i] = log(expected_popsize[k]) + calculate_hazard(dist, sigma[k]); 
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop[i];
    for (int j = 0; j < p->K; ++j) {
      if (j == k) {
        continue;
      }
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // define theta_i as the sentinel area * the mean hazard. Calculate
    // log(theta_i) and add theta_i to running sum
    double log_theta_i = log_hazard_sum;
    
    // add necessary terms to loglikelihood
    loglike_prop += lgamma(d->tested[i] + 1) - lgamma(d->positive[i] + 1) - lgamma(d->tested[i] - d->positive[i]  + 1)
                    + d->positive[i]*(log_theta_i) - d->tested[i]*log(1 + exp(log_theta_i));
    
  }
  
  return loglike_prop;
}

//------------------------------------------------
// calculate log-likelihood under point pattern model given new proposed source
double Particle::calculate_loglike_source_points(double source_lon_prop, double source_lat_prop, int k){

  // initialise running values
  double loglike_prop = 0;

  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // get distance from proposed source to data point i
    double dist = l->get_data_dist(source_lon_prop, source_lat_prop, i);
    dist_source_data_prop[i] = dist;
    
    // calculate bivariate height of data point i from proposed source.
    log_hazard_height_prop[i] = log(source_weights[k]) + calculate_hazard(dist, sigma[k]); 
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop[i];
    for (int j = 0; j < p->K; ++j) {
      if (j == k) {
        continue;
      }
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // add necessary terms to loglikelihood
    loglike_prop += log_hazard_sum;
  }

  return loglike_prop;
}

//------------------------------------------------
// update source locations
void Particle::update_sources(bool robbins_monro_on, int iteration) {

  // loop through all sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new source location
    double source_lon_prop = rnorm1(source_lon[k], source_propSD[k]);
    double source_lat_prop = rnorm1(source_lat[k], source_propSD[k]);
    
    // check proposed source within defined range
    if (source_lon_prop <= p->min_lon || source_lon_prop >= p->max_lon ||
        source_lat_prop <= p->min_lat || source_lat_prop >= p->max_lat) {
      
      // auto-reject proposed move
      if (robbins_monro_on) {
        source_propSD[k] = exp(log(source_propSD[k]) - source_rm_stepsize*0.234/sqrt(iteration));
      }
      continue;
    }
    
    // calculate new logprior and loglikelihood
    double logprior_prop = calculate_logprior_source(source_lon_prop, source_lat_prop);
    double loglike_prop;
    
    // update source
    loglike_prop = calculate_loglike_source(source_lon_prop, source_lat_prop, k);
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update source
      source_lon[k] = source_lon_prop;
      source_lat[k] = source_lat_prop;
      
      // update stored distances and hazard values
      for (int i = 0; i < d->n; ++i) {
        dist_source_data[i][k] = dist_source_data_prop[i];
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood and prior
      logprior = logprior_prop;
      loglike = loglike_prop;

      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        source_propSD[k] = exp(log(source_propSD[k]) + source_rm_stepsize*(1 - 0.234)/sqrt(iteration));
        source_accept_burnin[k]++;
      } else {
        source_accept_sampling[k]++;
      }
    
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        source_propSD[k] = exp(log(source_propSD[k]) - source_rm_stepsize*0.234/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
  
  }  // end loop through sources
  
}

//------------------------------------------------
// update independent sigma under a Poisson model for each source
void Particle::update_sigma_pois(bool robbins_monro_on, int iteration) {
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double sigma_prop = rnorm1(sigma[k], sigma_propSD[k]);
    sigma_prop = (sigma_prop < 0) ? -sigma_prop : sigma_prop;
    sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
    
    // initialise running values
    double loglike_prop = 0;
    double theta_sum = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = calculate_hazard(dist, sigma_prop);       
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log_hazard_height_prop[i];
      for (int j = 0; j < p->K; ++j) {
        if (j == k) {
          continue;
        }
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      
      
      // define theta_i as the sentinel area * the mean hazard. Calculate
      // log(theta_i) and add theta_i to running sum
      double log_theta_i = log_sentinel_area + log_hazard_sum - log_K;
      theta_sum += exp(log_theta_i);
      
      // add necessary terms to loglikelihood
      loglike_prop += d->counts[i]*log_theta_i - lgamma(d->counts[i] + 1);
      
    }
    
    
    // complete loglikelihood
    if (p->ep_prior_sd == 0) {
      loglike_prop += counts_total*log(p->ep_prior_mean) - p->ep_prior_mean*theta_sum;
    } else {
      double gamma_shape = p->ep_prior_shape;
      double gamma_rate = p->ep_prior_rate;
      
      loglike_prop += gamma_shape*log(gamma_rate) - (gamma_shape + counts_total)*log(gamma_rate + theta_sum) 
                      + lgamma(gamma_shape + counts_total) - lgamma(gamma_shape);
    }
    //-----------------------------------------------------------------------------------------------------------------------
    
    // calculate priors
    double logprior = dlnorm1(sigma[k], p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    double logprior_prop = dlnorm1(sigma_prop, p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update sigma for this source
      sigma[k] = sigma_prop;
      
      // if sigma model is single, stop looping over K and update all sigma
      if(p->sigma_model == 1){
        // set same sigma for all sources
        for (int k = 1; k < p->K; ++k) {
          sigma[k] = sigma_prop;
        }
      }
      
      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) + sigma_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        sigma_accept_burnin[k]++;
      } else {
        sigma_accept_sampling[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) - sigma_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
    if(p->sigma_model == 1){
      // break loop over sources
      break;
    }
  
  }  // end loop over sources
  
}

//------------------------------------------------
// update independent sigma under a Poisson model for each source
void Particle::update_sigma_pois_ind_exp_pop(bool robbins_monro_on, int iteration) {
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double sigma_prop = rnorm1(sigma[k], sigma_propSD[k]);
    sigma_prop = (sigma_prop < 0) ? -sigma_prop : sigma_prop;
    sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
    
    // initialise new likelihood
    double loglike_prop = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = log(expected_popsize[k]) + calculate_hazard(dist, sigma_prop);
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log_hazard_height_prop[i];
      for (int j = 0; j < p->K; ++j) {
        if (j == k) {
          continue;
        }
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      
      // define the rate lambda of the Poisson process at this sentinel site,
      // while remaining in log space
      double log_lambda = log_sentinel_area + log_hazard_sum;
      
      // calculate the Poisson log-probability of the counts at this sentinel site
      loglike_prop += d->counts[i]*log_lambda - exp(log_lambda) - lgamma(d->counts[i] + 1);
    }
    //-----------------------------------------------------------------------------------------------------------------------
    
    // calculate priors
    double logprior = dlnorm1(sigma[k], p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    double logprior_prop = dlnorm1(sigma_prop, p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update sigma for this source
      sigma[k] = sigma_prop;
      
      if(p->sigma_model == 1){
        // update sigma for all sources
        for (int k = 0; k < p->K; ++k) {
          sigma[k] = sigma_prop;
        }
      }
  
      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) + sigma_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        sigma_accept_burnin[k]++;
      } else {
        sigma_accept_sampling[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) - sigma_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
    if(p->sigma_model == 1){
      break;
    }
  }  // end loop over sources
  
}

//------------------------------------------------
// update independent sigma under a binomial model for each source
void Particle::update_sigma_binom(bool robbins_monro_on, int iteration) {
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double sigma_prop = rnorm1(sigma[k], sigma_propSD[k]);
    sigma_prop = (sigma_prop < 0) ? -sigma_prop : sigma_prop;
    sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
    
    // initialise running values
    double loglike_prop = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = log(expected_popsize[k]) + calculate_hazard(dist, sigma_prop);
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log_hazard_height_prop[i];
      for (int j = 0; j < p->K; ++j) {
        if (j == k) {
          continue;
        }
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      
      // define theta_i as the sentinel area * the mean hazard. Calculate
      // log(theta_i) and add theta_i to running sum
      double log_theta_i = log_hazard_sum;

      // add necessary terms to loglikelihood
      loglike_prop += lgamma(d->tested[i] + 1) - lgamma(d->positive[i] + 1) - lgamma(d->tested[i] - d->positive[i]  + 1)
        + d->positive[i]*(log_theta_i) - d->tested[i]*log(1 + exp(log_theta_i));
      
    }
    
    // calculate priors
    double logprior = dlnorm1(sigma[k], p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    double logprior_prop = dlnorm1(sigma_prop, p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update sigma for this source
      sigma[k] = sigma_prop;
      
      if(p->sigma_model==1){
        // update sigma for all sources
        for (int k = 0; k < p->K; ++k) {
          sigma[k] = sigma_prop;
        }
      }
    
      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) + sigma_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        sigma_accept_burnin[k]++;
      } else {
        sigma_accept_sampling[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) - sigma_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    if(p->sigma_model==1){
        break;
    }
  }  // end loop over sources
  
}

//------------------------------------------------
// update sigma under a point pattern model 
void Particle::update_sigma_points(bool robbins_monro_on, int iteration){
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double sigma_prop = rnorm1(sigma[k], sigma_propSD[k]);
    sigma_prop = (sigma_prop < 0) ? -sigma_prop : sigma_prop;
    sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
    
    // initialise running values
    double loglike_prop = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = log(source_weights[k]) + calculate_hazard(dist, sigma_prop);       
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log_hazard_height_prop[i];
      for (int j = 0; j < p->K; ++j) {
        if (j == k) {
          continue;
        }
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      
      // add necessary terms to loglikelihood
      loglike_prop += log_hazard_sum;
      
    }
    
    //--------------------------------------------------------------------------
    
    // calculate priors
    double logprior = dlnorm1(sigma[k], p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    double logprior_prop = dlnorm1(sigma_prop, p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update sigma for this source
      sigma[k] = sigma_prop;
      
      // if sigma model is single, stop looping over K and update all sigma
      if(p->sigma_model == 1){
        // set same sigma for all sources
        for (int k = 1; k < p->K; ++k) {
          sigma[k] = sigma_prop;
        }
      }
      
      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) + sigma_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        sigma_accept_burnin[k]++;
      } else {
        sigma_accept_sampling[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) - sigma_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
    if(p->sigma_model == 1){
      // break loop over sources
      break;
    }

  }  // end loop over sources
}

//------------------------------------------------
// update expected popsize under a Poisson model
void Particle::update_expected_popsize_pois_single() {
  
  // sum of Poisson rate over sentinel sites
  double lambda_total = 0;
  for (int i = 0; i < d->n; ++i) {
    for (int k = 0; k < p->K; ++k) {
      lambda_total += exp(log_sentinel_area + log_hazard_height[i][k] - log_K);
    }
  }
  
  // draw new expected population size
  double posterior_shape = p->ep_prior_shape + beta*counts_total;
  double posterior_rate = p->ep_prior_rate + beta*lambda_total;
  expected_popsize[0] = rgamma1(posterior_shape, posterior_rate);
  
  // update ep for all sources
  for (int k = 0; k < p->K; ++k) {
  expected_popsize[k] = expected_popsize[0];
  }
}

//------------------------------------------------
// update independent ep under a Poisson model for each source
void Particle::update_expected_popsize_pois_independent(bool robbins_monro_on, int iteration) {

  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
      double ep_prop = rnorm1(expected_popsize[k], ep_propSD[k]);
      ep_prop = (ep_prop < 0) ? -ep_prop : ep_prop;
      ep_prop = (ep_prop < UNDERFLO) ? UNDERFLO : ep_prop;
      
      // initialise new likelihood
      double loglike_prop = 0;
      
      // loop through sentinel sites
      for (int i = 0; i < d->n; ++i) {
        
        // recalculate hazard given new ep
        double dist = dist_source_data[i][k];
        log_hazard_height_prop[i] = log(ep_prop) + calculate_hazard(dist, sigma[k]);
        
        // sum hazard over sources while remaining in log space
        double log_hazard_sum = log_hazard_height_prop[i];
        for (int j = 0; j < p->K; ++j) {
          if (j == k) {
            continue;
          }
          if (log_hazard_sum < log_hazard_height[i][j]) {
            log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
          } else {
            log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
          }
        }
        
        // define the rate lambda of the Poisson process at this sentinel site,
        // while remaining in log space
        double log_lambda = log_sentinel_area + log_hazard_sum; 
        
        // calculate the Poisson log-probability of the counts at this sentinel site
        loglike_prop += d->counts[i]*log_lambda - exp(log_lambda) - lgamma(d->counts[i] + 1);
      }
    //-----------------------------------------------------------------------------------------------------------------------
    
    // calculate priors
    double logprior = log(dgamma1(expected_popsize[k], p->ep_prior_shape, p->ep_prior_rate));
    double logprior_prop = log(dgamma1(ep_prop, p->ep_prior_shape, p->ep_prior_rate));
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      // update ep for this source
      expected_popsize[k] = ep_prop;

      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }

      // update likelihood
      loglike = loglike_prop;

      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD[k] = exp(log(ep_propSD[k]) + ep_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        ep_accept_burnin[k]++;
      } else {
        ep_accept_sampling[k]++;
      }

    } else {
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD[k] = exp(log(ep_propSD[k]) - ep_rm_stepsize*0.44/sqrt(iteration));
      }

    }  // end Metropolis-Hastings step

  }  // end loop over sources

}

//------------------------------------------------
// update expected popsize under a binomial model
void Particle::update_expected_popsize_binom(bool robbins_monro_on, int iteration) {
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double ep_prop = rnorm1(expected_popsize[k], ep_propSD[k]);
    if (ep_prop < 0) {
      ep_prop *= -1;
    }
    if (ep_prop < UNDERFLO) {
      ep_prop = UNDERFLO;
    }
    
    // initialise running values
    double loglike_prop = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      // recalculate hazard given new ep
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = log(ep_prop) + calculate_hazard(dist, sigma[k]);
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log_hazard_height_prop[i];
      for (int j = 0; j < p->K; ++j) {
        
        if (j == k) {
          continue;
        } 
        
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      
      // define theta_i as the sentinel area * the mean hazard. Calculate
      // log(theta_i) and add theta_i to running sum
      double log_theta_i = log_hazard_sum;
      
      // add necessary terms to loglikelihood
      loglike_prop += lgamma(d->tested[i] + 1) - lgamma(d->positive[i] + 1) - lgamma(d->tested[i] - d->positive[i]  + 1)
        + d->positive[i]*(log_theta_i) - d->tested[i]*log(1 + exp(log_theta_i));
    }
    
    // calculate priors
    double logprior = log(dgamma1(expected_popsize[k], p->ep_prior_shape / double(p->K), p->ep_prior_rate));
    double logprior_prop = log(dgamma1(ep_prop, p->ep_prior_shape / double(p->K), p->ep_prior_rate));
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update sigma for this source
      expected_popsize[k] = ep_prop;
      
      if (p->ep_model==1){
        // update expected popsize for all sources
        for (int k = 0; k < p->K; ++k) {
          expected_popsize[k] = ep_prop;
        }
      }
      
      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD[k] = exp(log(ep_propSD[k]) + ep_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        ep_accept_burnin[k]++;
      } else {
        ep_accept_sampling[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD[k] = exp(log(ep_propSD[k]) - ep_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
    // if single expected popsize model stop here
    if (p->ep_model == 1){
      break;
    }
  
  }  // end loop over sources
  
}

//------------------------------------------------
// update weigths for a finite mixture using point pattern data
void Particle::update_weights_point_pattern(bool robbins_monro_on, int iteration) {
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double ep_prop = rnorm1(expected_popsize[k], ep_propSD[k]);
    if (ep_prop < 0) {
      ep_prop *= -1;
    }
    if (ep_prop < UNDERFLO) {
      ep_prop = UNDERFLO;
    }
    
    double ep_prop_total = ep_total - expected_popsize[k] + ep_prop; 
    
    // update all source weights
    for (int l = 0; l < p->K; ++l) {
      if (l == k) {
        source_weight_prop[l] = ep_prop/ep_prop_total; 
      } else{
        source_weight_prop[l] = expected_popsize[l]/ep_prop_total; 
      } 
    }
    
    // initialise running values
    double loglike_prop = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      for (int l = 0; l < p->K; ++l) {
        // recalculate hazard of each source given new weights
        double dist = dist_source_data[i][l];
        log_hazard_height_prop2[i][l] = log(source_weight_prop[l]) + calculate_hazard(dist, sigma[l]);       
      } 
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log(0);
      
      for (int j = 0; j < p->K; ++j) {
        if (log_hazard_sum < log_hazard_height_prop2[i][j]) {
          log_hazard_sum = log_hazard_height_prop2[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height_prop2[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height_prop2[i][j] - log_hazard_sum));
        }
      }
      
      // add additional term to log likelihood 
      loglike_prop += log_hazard_sum; 

    }

    // calculate priors (uniform prior on weights)
    double logprior = 0;
    double logprior_prop = 0;
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {

      // update expected pop size for this source
      expected_popsize[k] = ep_prop;
      ep_total = ep_prop_total;
        
      // update the weights for each source
      for (int j = 0; j < p->K; ++j) {
        source_weights[j] = source_weight_prop[j];
      }
      
      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        for (int j = 0; j < p->K; ++j) {
          log_hazard_height[i][j] = log_hazard_height_prop2[i][j];
        }
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD[k] = exp(log(ep_propSD[k]) + ep_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        ep_accept_burnin[k]++;
      } else {
        ep_accept_sampling[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD[k] = exp(log(ep_propSD[k]) - ep_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
  }  // end loop over sources
  
}

//------------------------------------------------
// update qmatrix
void Particle::update_qmatrix() {
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // dependent on data type - skip if no observation at this site
    if(d->data_type == 1){
      if (d->counts[i] == 0) {
        continue;
      }     
    } else if(d->data_type == 2){
      if (d->positive[i] == 0) {
        continue;
      }
    }
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = -OVERFLO;
    for (int j = 0; j < p->K; ++j) {
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // qmatrix equal to normalised hazard
    for (int j = 0; j < p->K; ++j) {
      log_qmatrix[i][j] = log_hazard_height[i][j] - log_hazard_sum;
      qmatrix[i][j] = exp(log_hazard_height[i][j] - log_hazard_sum);
    }
    
  }  // end loop through sentinel sites
  
}

//------------------------------------------------
// solve label switching problem
void Particle::solve_label_switching(const vector<vector<double>> &log_qmatrix_running) {
  
  // recalculate cost matrix
  for (int k1 = 0; k1 < p->K; ++k1) {
    fill(cost_mat[k1].begin(), cost_mat[k1].end(), 0);
    for (int k2 = 0; k2 < p->K; ++k2) {
      for (int i = 0; i < d->n; ++i) {
        // update cost matrix based on data type and assuming data is positive
        if (d->data_type == 1 && d->counts[i] > 0) { 
          // count data
          cost_mat[k1][k2] += qmatrix[i][label_order[k1]]*(log_qmatrix[i][label_order[k1]] 
                              - log_qmatrix_running[i][k2]);
        } else if(d->data_type == 2 && d->positive[i] > 0){ 
          // prevalence data
          cost_mat[k1][k2] += qmatrix[i][label_order[k1]]*(log_qmatrix[i][label_order[k1]] 
                              - log_qmatrix_running[i][k2]);
        } else if(d->data_type == 3){ 
          // point pattern  
          cost_mat[k1][k2] += qmatrix[i][label_order[k1]]*(log_qmatrix[i][label_order[k1]] 
                              - log_qmatrix_running[i][k2]);
        }
      
      }
    }
  }
  
  // find best permutation of current labels using Hungarian algorithm
  best_perm = hungarian(cost_mat, edges_left, edges_right, blocked_left, blocked_right);
  
  // define best_perm_order
  for (int k = 0; k < p->K; ++k) {
    best_perm_order[best_perm[k]] = k;
  }
  
  // replace old label order with new
  for (int k = 0; k < p->K; ++k) {
    label_order_new[k] = label_order[best_perm_order[k]];
  }
  label_order = label_order_new;
  
}

//------------------------------------------------
// calculate log-likelihood under Poisson model given new proposed source
double Particle::calculate_loglike_source_negative_binomial_indpendent_lambda(double source_lon_prop, double source_lat_prop, int k) {
  
  // initialise running values
  double loglike_prop = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // get distance from proposed source to data point i
    double dist = l->get_data_dist(source_lon_prop, source_lat_prop, i);
    dist_source_data_prop[i] = dist;
    
    // calculate bivariate height of data point i from proposed source.
    log_hazard_height_prop[i] = log(expected_popsize[k]) + calculate_hazard(dist, sigma[k]); 
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop[i];
    for (int j = 0; j < p->K; ++j) {
      if (j == k) {
        continue;
      }
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // define theta_i as the sentinel area * the mean hazard. Calculate log(theta_i) 
    double log_theta_i = log_sentinel_area + log_hazard_sum;
    loglike_prop += lgamma(1/alpha + d->counts[i]) - lgamma(d->counts[i] + 1) -
                    lgamma(1/alpha) + (1/alpha)*log(1/alpha) + 
                    d->counts[i]*log_theta_i - 
                    (1/alpha + d->counts[i])*log(1/alpha + exp(log_theta_i));
  }
  
  return loglike_prop;
}

//------------------------------------------------
// update independent sigma under a Poisson model for each source
void Particle::update_sigma_negative_binomial_ind_exp_pop(bool robbins_monro_on, int iteration) {
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double sigma_prop = rnorm1(sigma[k], sigma_propSD[k]);
    sigma_prop = (sigma_prop < 0) ? -sigma_prop : sigma_prop;
    sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
    
    // initialise new likelihood
    double loglike_prop = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = log(expected_popsize[k]) + calculate_hazard(dist, sigma_prop);
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log_hazard_height_prop[i];
      for (int j = 0; j < p->K; ++j) {
        if (j == k) {
          continue;
        }
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      
      // define theta_i as the sentinel area * the mean hazard. Calculate log(theta_i) 
      double log_theta_i = log_sentinel_area + log_hazard_sum;
      loglike_prop += lgamma(1/alpha + d->counts[i]) - lgamma(d->counts[i] + 1) -
                      lgamma(1/alpha) + (1/alpha)*log(1/alpha) + 
                      d->counts[i]*log_theta_i - 
                      (1/alpha + d->counts[i])*log(1/alpha + exp(log_theta_i));
    }
    
    //-----------------------------------------------------------------------------------------------------------------------
    
    // calculate priors
    double logprior = dlnorm1(sigma[k], p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    double logprior_prop = dlnorm1(sigma_prop, p->sigma_prior_meanlog, p->sigma_prior_sdlog);
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update sigma for this source
      sigma[k] = sigma_prop;
      
      if(p->sigma_model == 1){
        // update sigma for all sources
        for (int k = 0; k < p->K; ++k) {
          sigma[k] = sigma_prop;
        }
      }
      
      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) + sigma_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        sigma_accept_burnin[k]++;
      } else {
        sigma_accept_sampling[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) - sigma_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
    if(p->sigma_model == 1){
      break;
    }
  }  // end loop over sources
  
}

//------------------------------------------------
// update independent ep under a negative binomial model for each source
void Particle::update_expected_popsize_negative_binomial_independent(bool robbins_monro_on, int iteration) {
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double ep_prop = rnorm1(expected_popsize[k], ep_propSD[k]);
    ep_prop = (ep_prop < 0) ? -ep_prop : ep_prop;
    ep_prop = (ep_prop < UNDERFLO) ? UNDERFLO : ep_prop;
    
    // initialise new likelihood
    double loglike_prop = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      // recalculate hazard given new ep
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = log(ep_prop) + calculate_hazard(dist, sigma[k]);
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log_hazard_height_prop[i];
      for (int j = 0; j < p->K; ++j) {
        if (j == k) {
          continue;
        }
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      
      // define theta_i as the sentinel area * the mean hazard. Calculate log(theta_i) 
      double log_theta_i = log_sentinel_area + log_hazard_sum;
      loglike_prop += lgamma(1/alpha + d->counts[i]) - lgamma(d->counts[i] + 1) -
                      lgamma(1/alpha) + (1/alpha)*log(1/alpha) + 
                      d->counts[i]*log_theta_i - 
                      (1/alpha + d->counts[i])*log(1/alpha + exp(log_theta_i));
    }
    
    //-----------------------------------------------------------------------------------------------------------------------
    
    // calculate priors
    double logprior = log(dgamma1(expected_popsize[k], p->ep_prior_shape, p->ep_prior_rate));
    double logprior_prop = log(dgamma1(ep_prop, p->ep_prior_shape, p->ep_prior_rate));
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      // update ep for this source
      expected_popsize[k] = ep_prop;
      
      if(p->ep_model == 1){
        // update sigma for all sources
        for (int k = 0; k < p->K; ++k) {
          expected_popsize[k] = ep_prop;
        }
      }
      
      // update stored hazard values
      for (int i = 0; i < d->n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD[k] = exp(log(ep_propSD[k]) + ep_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        ep_accept_burnin[k]++;
      } else {
        ep_accept_sampling[k]++;
      }
      
    } else {
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD[k] = exp(log(ep_propSD[k]) - ep_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
    if(p->ep_model == 1){
      break;
    }
  }  // end loop over sources
  
}

//------------------------------------------------
// update independent ep under a negative binomial model for each source
void Particle::update_alpha_negative_binomial(bool robbins_monro_on, int iteration) {
  
  // propose new value
  double alpha_prop = rnorm1(alpha, alpha_propSD);
  alpha_prop = (alpha_prop < 0) ? -alpha_prop : alpha_prop;
  alpha_prop = (alpha_prop < UNDERFLO) ? UNDERFLO : alpha_prop;
  
  // initialise new likelihood
  double loglike_prop = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log(0);
    
    for (int j = 0; j < p->K; ++j) {
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // define theta_i as the sentinel area * the mean hazard. Calculate log(theta_i) 
    double log_theta_i = log_sentinel_area + log_hazard_sum;
    loglike_prop += lgamma(1/alpha_prop + d->counts[i]) - lgamma(d->counts[i] + 1) -
                    lgamma(1/alpha_prop) + (1/alpha_prop)*log(1/alpha_prop) + 
                    d->counts[i]*log_theta_i - 
                    (1/alpha_prop + d->counts[i])*log(1/alpha_prop + exp(log_theta_i));
  }
  
  //-----------------------------------------------------------------------------------------------------------------------
  // calculate priors
  double logprior = dlnorm1(alpha, p->alpha_prior_meanlog, p->alpha_prior_sdlog);
  double logprior_prop = dlnorm1(alpha_prop, p->alpha_prior_meanlog, p->alpha_prior_sdlog);
  
  // Metropolis-Hastings ratio
  double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
  
  // Metropolis-Hastings step
  if (log(runif_0_1()) < MH_ratio) {
    
    // update alpha for this source
    alpha = alpha_prop;
    
    // update likelihood
    loglike = loglike_prop;
    
    // Robbins-Monro positive update (on the log scale)
    if (robbins_monro_on) {
      alpha_propSD = exp(log(alpha_propSD) + alpha_rm_stepsize*(1 - 0.44)/sqrt(iteration));
      alpha_accept_burnin++;
    } else {
      alpha_accept_sampling++;
    }
    
  } else {
    // Robbins-Monro negative update (on the log scale)
    if (robbins_monro_on) {
      alpha_propSD = exp(log(alpha_propSD) - alpha_rm_stepsize*0.44/sqrt(iteration));
    }
    
  }  // end Metropolis-Hastings step
  
}

//------------------------------------------------
// calculate the log hazard height dependent on the kernel
double Particle::calculate_hazard(double dist, double single_scale) {
  
  double hazard_height = 0.0;  
  
  // update source based on dispersal kernel model
  if (p->dispersal_model == 1) {
    
    // calculate bivariate NORMAL height 
    // equivalent to the univariate normal density of the euclidian distance 
    // between source and data multiplied by the univariate normal density of 
    // 0 from 0 (the latter is needed to make it a bivariate density). Finally, 
    // in log space densities are summed not multiplied.
    hazard_height = dnorm1(dist, 0, single_scale, true) + dnorm1(0, 0, single_scale, true);  
    
  } else if (p->dispersal_model == 2) {
    
    // calculate bivariate CAUCHY height 
    hazard_height = - log(2*M_PI) + 0.5*log(single_scale) - 1.5*log(pow(dist, 2) + single_scale);
    
  } else if (p->dispersal_model == 3) {
    
    // calculate bivariate LAPLACE height
    // An accurate version of the laplace density is obtained using the bessel 
    // function below - it comes with underflow issues
    // hazard_height = -log(M_PI) - 2*log(single_scale) + log(bessel(sqrt(2)*dist*pow(single_scale, -1), 0));
    hazard_height = - 0.5*log(M_PI) - 0.75*log(2*single_scale) - 0.5*log(dist) - dist*sqrt(2*pow(single_scale, -1));
    
  } else if (p->dispersal_model == 4) {
    
    // calculate bivariate Student's t height
    // hazard_height = log(cell_area);
  }
  
  return hazard_height;
}

// // //------------------------------------------------
// // // update expected popsize under a Poisson model for each source TODO COMMENTS
// void Particle::update_expected_popsize_pois_independent(bool robbins_monro_on) {
// 
//   int domain_max = 1e3;
//   int domain_size = 100;
//   pop_size_domain = vector<double>(domain_max, 0);
//   vector<double> current_eps = expected_popsize;
// 
//   // generate domain 
//   for (int domain = 0; domain < domain_size; ++domain) {
//     pop_size_domain[domain] = domain*(domain_max/domain_size);
//   }
// 
//   // // loop through sources
//   for (int k = 0; k < p->K; ++k) {
// 
//    cum_sum_density = vector<double>(domain_size, 0);
//    cum_sum_normalised = vector<double>(domain_size, 0);
//    double density_sum = 0;
// 
//   // calculate the density of a single expected pop size
//    for (int i = 1; i < domain_size; ++i) {
//        // initiate values that make up loglikelihood
//        double exp_component = 0;
//        double product_component = 0;
// 
//        // loop over sentinel sites to sum over their individual loglikelihoods
//        for (int j = 0; j < d->n; ++j) {
// 
//          // TODO, deal with underflow here
//          exp_component += exp(log_hazard_height[j][k] - log(current_eps[k]));
// 
//          // prep to sum over hazards with new expected pop size 
//          double sentinel_component = log(pop_size_domain[i]) + log_hazard_height[j][k] - log(current_eps[k]);;
// 
//          for (int b = 0; b < p->K; ++b) {
//            if(b == k){
//              continue;
//            } 
//            if (sentinel_component < log_hazard_height[j][b]) {
//              sentinel_component = log_hazard_height[j][b] + log(1 + exp(sentinel_component - log_hazard_height[j][b]));
//            } else {
//              sentinel_component = sentinel_component + log(1 + exp(log_hazard_height[j][b] - sentinel_component));
//            }
//          }
// 
//          product_component += sentinel_component*d->counts[j];
// 
//       } // finish looping over sentinel sites and combine loglikelihood components
//       exp_component = pop_size_domain[i]*exp(log_sentinel_area)*exp_component;
// 
//       double single_lambda_loglikelihood = beta*(product_component - exp_component);
//       double single_lambda_logprior = log(dgamma1(pop_size_domain[i],  p->ep_prior_shape,  p->ep_prior_rate)); 
//       double single_lambda_logposterior = single_lambda_loglikelihood + single_lambda_logprior;
//       cum_sum_density[i] = cum_sum_density[i - 1] + exp(single_lambda_logposterior);
//       density_sum += exp(single_lambda_logposterior);
//       if(robbins_monro_on){
//           print(product_component);
//           print(i);
//         }
//       }
// 
//       // normalise
//       for(int i = 0; i < domain_size; ++i) {
//            cum_sum_normalised[i] = cum_sum_density[i]/density_sum;
//         }
// 
//         // draw new expected population size
//         double random_val = runif_0_1();
// 
//         // the index of the value nearest to r 
//         int lo = 0;
// 
//         for(int i = 0; i < domain_size; ++i) {
//             if(cum_sum_normalised[i] < random_val){
//                 lo++; 
//               } else {
//                   break;
//                 }
//               }
// 
//               // how far away from the nearest value is r? 
//               double gap = R::imax2(random_val - cum_sum_normalised[lo], 0); 
// 
//               // linear interpolation: y = c + mx from the lower value
//               double lif = gap/(cum_sum_normalised[lo + 1] - cum_sum_normalised[lo]);
//               expected_popsize[k] = pop_size_domain[lo] + lif * (pop_size_domain[lo + 1] - pop_size_domain[lo]);
// 
//               // update log hazard height given new expected pop size
//               for (int j = 0; j < d->n; ++j) {
//                     log_hazard_height[j][k] = log(expected_popsize[k]) + log_hazard_height[j][k] - log(current_eps[k]);
//                 }
//               }
//             }
