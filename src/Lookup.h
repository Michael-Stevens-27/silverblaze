
#pragma once

#include <Rcpp.h>

#include "Data.h"
#include "Parameters.h"

//------------------------------------------------
// class containing lookup tables
class Lookup {
  
public:
  
  // pointers
  Data * d;
  Parameters * p;
  
  // lookup table
  std::vector<double> lookup_dist;
  
  // constructors
  Lookup() {};
  Lookup(Data &data, Parameters &params);
  
  // other methods
  void recalc();
  double get_data_dist(double source_lon, double source_lat, int data_index);
  
};
