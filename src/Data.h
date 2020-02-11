
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing data
class Data {
  
public:
  
  std::vector<double> sentinel_lon;
  std::vector<double> sentinel_lat;
  std::vector<int> sentinel_counts;
  std::vector<int> total_counts;
  int n;
  
  // constructors
  Data() {};
  Data(const Rcpp::List &args);
  
};
