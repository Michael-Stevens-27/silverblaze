
#pragma once

#include <Rcpp.h>
#include "Parameters.h"

//------------------------------------------------
// class containing spatial prior
class Spatial_prior {

public:
  
  // pointer to parameters
  Parameters * p;
  
  // spatial prior
  static std::vector<bool> spatial_prior_mask;

  // constructors
  Spatial_prior() {};
  Spatial_prior(const Rcpp::List &args, Parameters &params);

  // other methods
  double get_value(double lon, double lat);

};
