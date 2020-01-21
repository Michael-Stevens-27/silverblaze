
#include "Data.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// constructor for Data class
Data::Data(const Rcpp::List &args) {
  
  sentinel_lon = rcpp_to_vector_double(args["longitude"]);
  sentinel_lat = rcpp_to_vector_double(args["latitude"]);
  sentinel_counts = rcpp_to_vector_int(args["counts"]);
  n = int(sentinel_counts.size());
}
