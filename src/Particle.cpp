
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
  
  // scaling factor on hazard surface, equivalent to the expected total
  // population size (both observed and unobserved) in unit time
  if (p->ep_prior_sd <= 0) {
    expected_popsize = p->ep_prior_mean;
  } else {
    if(p->model_type == 1){
      expected_popsize = rgamma1(p->ep_prior_shape, p->ep_prior_rate);
    } else if(p->model_type == 2){
      expected_popsize = abs(rnorm1(p->ep_prior_mean, p->ep_prior_sd));
    }
  }
  log_expected_popsize = log(expected_popsize);
  
  // qmatrices
  log_qmatrix = vector<vector<double>>(d->n, vector<double>(p->K));
  qmatrix = vector<vector<double>>(d->n, vector<double>(p->K));
  
  // proposal standard deviations
  source_propSD = vector<double>(p->K, 0.01);
  sigma_propSD = vector<double>(p->K, 1.0);
  ep_propSD = 1.0;
  
  // Robbins-Monro stepsize constants
  source_rm_stepsize = 5.0;
  sigma_rm_stepsize = 4.0;
  ep_rm_stepsize = 4.0;

  // misc constants
  // area around a sentinel site
  log_sentinel_area = LOG_PI + 2*log(p->sentinel_radius);
  // sum of counts over all sentinel sites
  counts_total = sum(d->sentinel_counts);
  tested_total = sum(d->total_counts);
    
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
  int ep_accept_burnin;
  int ep_accept_sampling;
  
}

//------------------------------------------------
// reset particle
void Particle::reset(double beta) {
  
  // reset beta value
  this->beta = beta;
  
  // initialise source locations
  for (int k = 0; k < p->K; ++k) {
    source_lon[k] = p->source_init[0];
    source_lat[k] = p->source_init[1];
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
    expected_popsize = p->ep_prior_mean;
  } else {
    if(p->model_type == 1){
      expected_popsize = rgamma1(p->ep_prior_shape, p->ep_prior_rate);
    } else if(p->model_type == 2){
      expected_popsize = abs(rnorm1(p->ep_prior_mean, p->ep_prior_sd));
    }
  }
  log_expected_popsize = log(expected_popsize);
  
  // initialise proposal standard deviations
  source_propSD = vector<double>(p->K, 0.01);
  sigma_propSD = vector<double>(p->K, 1.0);
  double ep_propSD = 1.0;
  
  // calculate initial likelihood. Calling calculate_loglike_source() on each
  // source updates the dist_source_data_prop and log_hazard_height_prop
  // objects, which can then be stored back into the final matrices. This is
  // equivalent to running a Metropolis-Hastings step in which the move is
  // guaranteed to be accepted
  for (int k = 0; k < p->K; ++k) {
    logprior = calculate_logprior_source(source_lon[k], source_lat[k]);
    // update source based on binomial or poisson model
    if (p->model_type == 1) {
      loglike = calculate_loglike_source(source_lon[k], source_lat[k], k);
    } else if (p->model_type == 2) {
      loglike = BINOM_calculate_loglike_source(source_lon[k], source_lat[k], k);
    }
    
    for (int i = 0; i < d->n; ++i) {
      dist_source_data[i][k] = dist_source_data_prop[i];
      log_hazard_height[i][k] = log_hazard_height_prop[i];
    }
  }
  
  // reset acceptance rates
  source_accept_burnin = vector<int>(p->K);
  source_accept_sampling = vector<int>(p->K);
  sigma_accept_burnin = vector<int>(p->K);
  sigma_accept_sampling = vector<int>(p->K);
  int ep_accept_burnin;
  int ep_accept_sampling;  
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
// calculate log-likelihood given new proposed source
double Particle::calculate_loglike_source(double source_lon_prop, double source_lat_prop, int k) {
  
  // initialise running values
  double loglike_prop = 0;
  double theta_sum = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // get distance from proposed source to data point i
    double dist = l->get_data_dist(source_lon_prop, source_lat_prop, i);
    dist_source_data_prop[i] = dist;
    
    // calculate bivariate normal height of data point i from proposed source.
    // This is equivalent to the univariate normal density of the distance
    // between data and source in lon, multiplied by the univariate normal
    // density of the distance between data and source in lat. Due to the
    // properties of normal distributions this is equivalent to the univariate
    // normal density of the euclidian distance between source and data,
    // multiplied by the univariate normal density of 0 from 0 (the latter is
    // needed to make it a bivariate density). Finally, in log space densities
    // are summed not multiplied.
    log_hazard_height_prop[i] = dnorm1(dist, 0, sigma[k], true) + dnorm1(0, 0, sigma[k], true);
    
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
    loglike_prop += d->sentinel_counts[i]*log_theta_i - lgamma(d->sentinel_counts[i] + 1);
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
// calculate log-likelihood given new proposed source under BIOMIAL model
double Particle::BINOM_calculate_loglike_source(double source_lon_prop, double source_lat_prop, int k) {
  
  // initialise running values
  double loglike_prop = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // get distance from proposed source to data point i
    double dist = l->get_data_dist(source_lon_prop, source_lat_prop, i);
    dist_source_data_prop[i] = dist;
    
    // calculate bivariate normal height of data point i from proposed source.
    // This is equivalent to the univariate normal density of the distance
    // between data and source in lon, multiplied by the univariate normal
    // density of the distance between data and source in lat. Due to the
    // properties of normal distributions this is equivalent to the univariate
    // normal density of the euclidian distance between source and data,
    // multiplied by the univariate normal density of 0 from 0 (the latter is
    // needed to make it a bivariate density). Finally, in log space densities
    // are summed not multiplied.
    log_hazard_height_prop[i] = dnorm1(dist, 0, sigma[k], true) + dnorm1(0, 0, sigma[k], true);
    
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
    
    // add necessary terms to loglikelihood
    loglike_prop += lgamma(d->total_counts[i] + 1) - lgamma(d->sentinel_counts[i] + 1) - lgamma(d->total_counts[i] - d->sentinel_counts[i]  + 1)
                    + d->sentinel_counts[i]*log(expected_popsize*exp(log_theta_i)) - d->total_counts[i]*log(1 + expected_popsize*exp(log_theta_i));
    
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
    
    // update source based on binomial or poisson model
    if (p->model_type == 1) {
      loglike_prop = calculate_loglike_source(source_lon_prop, source_lat_prop, k);
    } else if (p->model_type == 2) {
      loglike_prop = BINOM_calculate_loglike_source(source_lon_prop, source_lat_prop, k);
    }
    
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
// update sigma
void Particle::update_sigma(bool robbins_monro_on, int iteration) {
  
  // return if prior is exact
  if (p->sigma_prior_sdlog == 0) {
    return;
  }
  // update source based on binomial or poisson model
  if (p->model_type == 1) { // Poisson
    // update single sigma or separately for each source
    if (p->sigma_model == 1) {
      update_sigma_single(robbins_monro_on, iteration);
    } else if (p->sigma_model == 2) {
      update_sigma_independent(robbins_monro_on, iteration);
    }
  } else if (p->model_type == 2) { // binomial
    // update single sigma or separately for each source
    if (p->sigma_model == 1) {
      BINOM_update_sigma_single(robbins_monro_on, iteration);
    } else if (p->sigma_model == 2) {
      BINOM_update_sigma_independent(robbins_monro_on, iteration);
    }
  }
}

//------------------------------------------------
// update single sigma for all sources
void Particle::update_sigma_single(bool robbins_monro_on, int iteration) {
  
  // propose new value
  double sigma_prop = rnorm1(sigma[0], sigma_propSD[0]);
  sigma_prop = (sigma_prop < 0) ? -sigma_prop : sigma_prop;
  sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;

  // initialise running values
  double loglike_prop = 0;
  double theta_sum = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // loop through sources
    for (int k = 0; k < p->K; ++k) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop2[i][k] = dnorm1(dist, 0, sigma_prop, true) + dnorm1(0, 0, sigma_prop, true);
    }
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop2[i][0];
    for (int k = 1; k < p->K; ++k) {
      if (log_hazard_sum < log_hazard_height_prop2[i][k]) {
        log_hazard_sum = log_hazard_height_prop2[i][k] + log(1 + exp(log_hazard_sum - log_hazard_height_prop2[i][k]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height_prop2[i][k] - log_hazard_sum));
      }
    }
    
    // define theta_i as the sentinel area * the mean hazard. Calculate
    // log(theta_i) and add theta_i to running sum
    double log_theta_i = log_sentinel_area + log_hazard_sum - log_K;
    theta_sum += exp(log_theta_i);
    
    // add necessary terms to loglikelihood
    loglike_prop += d->sentinel_counts[i]*log_theta_i - lgamma(d->sentinel_counts[i] + 1);
    
  }
  
  // complete loglikelihood
  if (p->ep_prior_sd == 0) {
    loglike_prop += counts_total*log(p->ep_prior_mean) - p->ep_prior_mean*theta_sum;
  } else {
    double gamma_shape = p->ep_prior_shape;
    double gamma_rate = p->ep_prior_rate;
    
    loglike_prop += gamma_shape*log(gamma_rate) - (gamma_shape + counts_total)*log(gamma_rate + theta_sum) + lgamma(gamma_shape + counts_total) - lgamma(gamma_shape);
  }
  
  // calculate priors
  double logprior = dlnorm1(sigma[0], p->sigma_prior_meanlog, p->sigma_prior_sdlog);
  double logprior_prop = dlnorm1(sigma_prop, p->sigma_prior_meanlog, p->sigma_prior_sdlog);
  
  // Metropolis-Hastings ratio
  double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
  
  // Metropolis-Hastings step
  if (log(runif_0_1()) < MH_ratio) {
    
    // update sigma for all sources
    for (int k = 0; k < p->K; ++k) {
      sigma[k] = sigma_prop;
    }
    
    // update stored hazard values
    log_hazard_height = log_hazard_height_prop2;
    
    // update likelihood
    loglike = loglike_prop;
    
    // Robbins-Monro positive update (on the log scale)
    if (robbins_monro_on) {
      sigma_propSD[0] = exp(log(sigma_propSD[0]) + sigma_rm_stepsize*(1 - 0.44)/sqrt(iteration));
      sigma_accept_burnin[0]++;
    } else {
      sigma_accept_sampling[0]++;
    }
    
  } else {
    
    // Robbins-Monro negative update (on the log scale)
    if (robbins_monro_on) {
      sigma_propSD[0] = exp(log(sigma_propSD[0]) - sigma_rm_stepsize*0.44/sqrt(iteration));
    }
    
  }  // end Metropolis-Hastings step
  
}

//------------------------------------------------
// update single sigma for all sources
void Particle::BINOM_update_sigma_single(bool robbins_monro_on, int iteration) {
  
  // propose new value
  double sigma_prop = rnorm1(sigma[0], sigma_propSD[0]);
  sigma_prop = (sigma_prop < 0) ? -sigma_prop : sigma_prop;
  sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
  
  // initialise running values
  double loglike_prop = 0;
  double theta_sum = 0;
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // loop through sources
    for (int k = 0; k < p->K; ++k) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop2[i][k] = dnorm1(dist, 0, sigma_prop, true) + dnorm1(0, 0, sigma_prop, true);
    }
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop2[i][0];
    for (int k = 1; k < p->K; ++k) {
      if (log_hazard_sum < log_hazard_height_prop2[i][k]) {
        log_hazard_sum = log_hazard_height_prop2[i][k] + log(1 + exp(log_hazard_sum - log_hazard_height_prop2[i][k]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height_prop2[i][k] - log_hazard_sum));
      }
    }
    
    // define theta_i as the sentinel area * the mean hazard. Calculate
    // log(theta_i) and add theta_i to running sum
    double log_theta_i = log_sentinel_area + log_hazard_sum - log_K;

    // add necessary terms to loglikelihood
    loglike_prop += lgamma(d->total_counts[i] + 1) - lgamma(d->sentinel_counts[i] + 1) - lgamma(d->total_counts[i] - d->sentinel_counts[i]  + 1)
                    + d->sentinel_counts[i]*log(expected_popsize*exp(log_theta_i)) - d->total_counts[i]*log(1 + expected_popsize*exp(log_theta_i));

  }
  
  // calculate priors
  double logprior = dlnorm1(sigma[0], p->sigma_prior_meanlog, p->sigma_prior_sdlog);
  double logprior_prop = dlnorm1(sigma_prop, p->sigma_prior_meanlog, p->sigma_prior_sdlog);
  
  // Metropolis-Hastings ratio
  double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
  
  // Metropolis-Hastings step
  if (log(runif_0_1()) < MH_ratio) {
    
    // update sigma for all sources
    for (int k = 0; k < p->K; ++k) {
      sigma[k] = sigma_prop;
    }
    
    // update stored hazard values
    log_hazard_height = log_hazard_height_prop2;
    
    // update likelihood
    loglike = loglike_prop;
    
    // Robbins-Monro positive update (on the log scale)
    if (robbins_monro_on) {
      sigma_propSD[0] = exp(log(sigma_propSD[0]) + sigma_rm_stepsize*(1 - 0.44)/sqrt(iteration));
      sigma_accept_burnin[0]++;
    } else {
      sigma_accept_sampling[0]++;
    }
    
  } else {
    
    // Robbins-Monro negative update (on the log scale)
    if (robbins_monro_on) {
      sigma_propSD[0] = exp(log(sigma_propSD[0]) - sigma_rm_stepsize*0.44/sqrt(iteration));
    }
    
  }  // end Metropolis-Hastings step
  
}


//------------------------------------------------
// update independent sigma for each source
void Particle::update_sigma_independent(bool robbins_monro_on, int iteration) {
  
  // loop through sources
  for (int k = 0; k < p->K; ++k) {
    
    // propose new value
    double sigma_prop = rnorm1(sigma[k], sigma_propSD[k]);
    sigma_prop = (sigma_prop < 0) ? -sigma_prop : sigma_prop;
    sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
    
    //-----------------------------------------------------------------------------------------------------------------------
    // initialise running values
    double loglike_prop = 0;
    double theta_sum = 0;
    
    // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = dnorm1(dist, 0, sigma_prop, true) + dnorm1(0, 0, sigma_prop, true);
      
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
      loglike_prop += d->sentinel_counts[i]*log_theta_i - lgamma(d->sentinel_counts[i] + 1);
      
    }
    
    // complete loglikelihood
    if (p->ep_prior_sd == 0) {
      loglike_prop += counts_total*log(p->ep_prior_mean) - p->ep_prior_mean*theta_sum;
    } else {
      double gamma_shape = p->ep_prior_shape;
      double gamma_rate = p->ep_prior_rate;
      
      loglike_prop += gamma_shape*log(gamma_rate) - (gamma_shape + counts_total)*log(gamma_rate + theta_sum) + lgamma(gamma_shape + counts_total) - lgamma(gamma_shape);
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
    
  }  // end loop over sources
  
}

//------------------------------------------------
// update independent sigma for each source
void Particle::BINOM_update_sigma_independent(bool robbins_monro_on, int iteration) {
  
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
      log_hazard_height_prop[i] = dnorm1(dist, 0, sigma_prop, true) + dnorm1(0, 0, sigma_prop, true);
      
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

      // add necessary terms to loglikelihood
      loglike_prop += lgamma(d->total_counts[i] + 1) - lgamma(d->sentinel_counts[i] + 1) - lgamma(d->total_counts[i] - d->sentinel_counts[i]  + 1)
                      + d->sentinel_counts[i]*log(expected_popsize*exp(log_theta_i)) - d->total_counts[i]*log(1 + expected_popsize*exp(log_theta_i));
      
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
    
  }  // end loop over sources
  
}

//------------------------------------------------
// update sigma
void Particle::update_expected_popsize(bool robbins_monro_on, int iteration) {
  
  // return if prior is exact
  if (p->ep_prior_sd == 0) {
    return;
  }
  
  // update expected_pop_size based on binomial or poisson model
  if (p->model_type == 1) { // Poisson 
    update_expected_popsize_gibbs();
  } else if (p->model_type == 2) { // Binomial 
    BINOM_update_expected_popsize(robbins_monro_on, iteration);
  }
}

//------------------------------------------------
// update separate sigma for each source
void Particle::update_expected_popsize_gibbs() {
  
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
  expected_popsize = rgamma1(posterior_shape, posterior_rate);
  log_expected_popsize = log(expected_popsize);
}

//------------------------------------------------
// update expected pop size using a MH step 
void Particle::BINOM_update_expected_popsize(bool robbins_monro_on, int iteration) {
  
    // propose new value
    double ep_prop = rnorm1(expected_popsize, ep_propSD);
    ep_prop = (ep_prop < 0) ? -ep_prop : ep_prop;
    ep_prop = (ep_prop < UNDERFLO) ? UNDERFLO : ep_prop;
    
    // initialise running values
    double loglike_prop = 0;
    
    // // loop through sentinel sites
    for (int i = 0; i < d->n; ++i) {
    
      // sum hazard over sources while remaining in log space
      double log_hazard_sum;
      for (int j = 0; j < p->K; ++j) {
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      // define theta_i as the sentinel area * the mean hazard. Calculate
      // log(theta_i) and add theta_i to running sum
      double log_theta_i = log_sentinel_area + log_hazard_sum - log_K;

      // add necessary terms to loglikelihood
      loglike_prop += lgamma(d->total_counts[i] + 1) - lgamma(d->sentinel_counts[i] + 1) - lgamma(d->total_counts[i] - d->sentinel_counts[i] + 1)
                      + d->sentinel_counts[i]*log(expected_popsize*exp(log_theta_i)) - d->total_counts[i]*log(1 + expected_popsize*exp(log_theta_i));
      }
    
    // calculate priors
    double logprior = dlnorm1(expected_popsize, p->ep_prior_meanlog, p->ep_prior_sdlog);
    double logprior_prop = dlnorm1(ep_prop, p->ep_prior_meanlog, p->ep_prior_sdlog);
    
    // Metropolis-Hastings ratio
    double MH_ratio = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update sigma for this source
      expected_popsize = ep_prop;
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD = exp(log(ep_propSD) + ep_rm_stepsize*(1 - 0.44)/sqrt(iteration));
        ep_accept_burnin++;
      } else {
        ep_accept_sampling++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        ep_propSD = exp(log(ep_propSD) - ep_rm_stepsize*0.44/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
}

//------------------------------------------------
// update qmatrix
void Particle::update_qmatrix() {
  
  // loop through sentinel sites
  for (int i = 0; i < d->n; ++i) {
    
    // skip if no observations at this site
    if (d->sentinel_counts[i] == 0) {
      continue;
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
        if (d->sentinel_counts[i] > 0) {
          cost_mat[k1][k2] += qmatrix[i][label_order[k1]]*(log_qmatrix[i][label_order[k1]] - log_qmatrix_running[i][k2]);
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
