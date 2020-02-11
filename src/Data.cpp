
#include "Data.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// constructor for Data class
Data::Data(const Rcpp::List &args) {
  
  sentinel_lon = rcpp_to_vector_double(args["longitude"]);
  sentinel_lat = rcpp_to_vector_double(args["latitude"]);
  counts = rcpp_to_vector_int(args["counts"]);
  tested = rcpp_to_vector_int(args["tested"]);
  positive = rcpp_to_vector_int(args["positive"]);
  n = int(counts.size());
  data_type = rcpp_to_int(args["data_type"]);
}
