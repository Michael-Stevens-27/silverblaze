
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing data
class Data {
  
public:
  
  std::vector<double> sentinel_lon;
  std::vector<double> sentinel_lat;
  std::vector<int> counts;
  std::vector<int> tested;
  std::vector<int> positive;
  int n;
  int data_type;
  
  // constructors
  Data() {};
  Data(const Rcpp::List &args);
  
};
