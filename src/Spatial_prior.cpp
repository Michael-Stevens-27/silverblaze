
#include "Spatial_prior.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// constructor for Parameters class
Spatial_prior::Spatial_prior(const Rcpp::List &args, Parameters &params) {
  
  // pointer to parameters
  p = &params;
  
  // get number of cells in each dimension
  vector<double> spatial_prior_values = rcpp_to_vector_double(args["spatial_prior_values"]);
  
  // populate mask 
  spatial_prior_mask = vector<double>(p->n_lat*p->n_lon, -OVERFLO);
  for (int i = 0; i < (p->n_lat*p->n_lon); ++i) {
    spatial_prior_mask[i] = spatial_prior_values[i];
  }
  
}

//------------------------------------------------
// get value
double Spatial_prior::get_value(double lon, double lat) {
  
  // convert lon/lat to index using precision
  int lon_index = floor((lon - p->min_lon)/p->res_lon);
  int lat_index;
  
  // pull the bugged index with the floor function 
  // or the proposed solution using ceil instead
  if(p->bugged){
    lat_index = floor((lat - p->min_lat)/p->res_lat);
    
  } else {
    lat_index = ceil((lat - p->min_lat)/p->res_lat);
  }

  int transform = (p->n_lat - lat_index)*p->n_lon + lon_index;

  // lookup value
  return spatial_prior_mask[transform];
  
}
