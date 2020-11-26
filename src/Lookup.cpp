
#include "Lookup.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// constructor
Lookup::Lookup(Data &data, Parameters &params) {
  
  // define pointers
  d = &data;
  p = &params;
  
}

//------------------------------------------------
// recalculate lookup table
void Lookup::recalc() {
  
  // print size of lookup table
  //if (!silent) {
  //  print("lookup table", sizeof(double)*n_lon*n_lat*n/1e6, "megabytes in size");
  //}
  
  // populate lookup table
  lookup_dist = vector<double>(p->n_lon * p->n_lat * d->n);
  int j = 0;
  for (int i_lon = 0; i_lon < p->n_lon; ++i_lon) {
    for (int i_lat = 0; i_lat < p->n_lat; ++i_lat) {
      for (int i = 0; i < d->n; ++i) {
        
        // store great circle distance between this lon/lat point and the data
        double lon = p->min_lon + (i_lon + 0.5)*p->res_lon;
        double lat = p->min_lat + (i_lat + 0.5)*p->res_lat;
        lookup_dist[j++] = gc_dist(lon, lat, d->sentinel_lon[i], d->sentinel_lat[i]);
      }
    }
  }
  
}

//------------------------------------------------
// get distance in km of a data point (indexed by data_index) from a given 
// source lon/lat
double Lookup::get_data_dist(std::vector<double> source_prop, int data_index) {
  
  // convert lon/lat to index
  int lon_index = floor((source_prop[0] - p->min_lon)/p->res_lon);
  int lat_index = floor((source_prop[1] - p->min_lat)/p->res_lat);

  double val = lookup_dist[(lon_index*p->n_lat + lat_index)*d->n + data_index]; 
  
  // lookup value
  return val;
}
