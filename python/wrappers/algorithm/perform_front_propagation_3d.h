#ifndef PERFORM_FRONT_PROPAGATION_3D_H
#define PERFORM_FRONT_PROPAGATION_3D_H

#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "perform_front_propagation_3d_config.h"
#include "fib.h"


struct point {
  point( int ii, int jj, int kk )
  { i = ii; j = jj; k = kk; }
  int i,j,k;
  double* pval;
};


struct Propagator
{
  typedef std::vector<point*> point_list;
  typedef bool (*T_callback_intert_node)(int i, int j, int k, int ii, int jj, int kk);

  unsigned int n, p, q;
  double* D = NULL;
  double* S = NULL;
  const double* W = NULL;
  double* Q = NULL;
  double* L = NULL;
  double* H = NULL;
  const double* start_points = NULL;
  double* end_points = NULL;
  const double* values = NULL;
  int nb_iter_max = 100000;
  int nb_start_points = 0;
  int nb_end_points = 0;
  fibheap_el** heap_pool = NULL;

  void perform_front_propagation_3d( T_callback_intert_node callback_insert_node = NULL );
  bool end_points_reached(const int i, const int j, const int k );
};

#endif
