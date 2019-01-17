/*

  perform_front_propagation_3d - perform a fast marching front propagation.
  
    'D' is a 2D array containing the value of the distance function to seed.
    'S' is a 2D array containing the state of each point : 
        -1 : dead, distance have been computed.
         0 : open, distance is being computed but not set.
         1 : far, distance not already computed.
    'W' is the weight matrix (inverse of the speed).
    'start_points' is a 3 x num_start_points matrix where k is the number of starting points.
    'H' is an heuristic (distance that remains to goal). This is a 2D matrix.
    
  Copyright (c) 2004 Gabriel Peyre
  
  Ported to python by Andrew

*/

#include <iostream>
#include <exception>

#include "labelfusion.h"
#include "perform_front_propagation_3d.h"


#define kDead -1
#define kOpen 0
#define kFar 1
#define ACCESS_ARRAY(a,i,j,k) a[(i)+n*(j)+n*p*(k)]
#define D_(i,j,k) ACCESS_ARRAY(D,i,j,k)
#define S_(i,j,k) ACCESS_ARRAY(S,i,j,k)
#define W_(i,j,k) ACCESS_ARRAY(W,i,j,k)
#define H_(i,j,k) ACCESS_ARRAY(H,i,j,k)
#define L_(i,j,k) ACCESS_ARRAY(L,i,j,k)
#define Q_(i,j,k) ACCESS_ARRAY(Q,i,j,k)
#define heap_pool_(i,j,k) ACCESS_ARRAY(heap_pool,i,j,k)
#define start_points_(i,s) start_points[(i)+3*(s)]
#define end_points_(i,s) end_points[(i)+3*(s)]


/*
  Performs a fast marching front propagation and computes the distance function from given
  starting points. This is a python port of the perform_front_propagation_3d tool written by Gabriel Peyre.
*/
pyarrayd performFrontPropagation3D(const pyarrayd &weights, const pyarrayd &start, int max_iters, const pyarrayd &values)
{
  Propagator prop;

  // weights
  prop.W = weights.data(0);
  prop.n = weights.shape(0);
  prop.p = weights.shape(1);
  prop.q = weights.shape(2);

  // start points
  prop.start_points = start.data(0);
  prop.nb_start_points = start.shape(1);
  if (start.shape(0) != 3) throw std::runtime_error("start points must be of shape 3 x N");

  // max iters
  prop.nb_iter_max = max_iters;

  // values
  prop.values = values.data(0);
  if (values.shape(0) != prop.nb_start_points) throw std::runtime_error("len of values does not equal number of start points");

  // init arrays
  prop.D = new double[weights.size()];
  prop.S = new double[weights.size()];
  prop.Q = new double[weights.size()];

  // propagate
  prop.perform_front_propagation_3d();

  // free S and Q since we won't be returning them
  GW_DELETEARRAY(prop.S)
  GW_DELETEARRAY(prop.Q)

  return makeArray({prop.n, prop.p, prop.q}, MemoryOrder::Fortran, prop.D);
}



int compare_points(void *x, void *y)
{
  point& a = *( (point*) x );
  point& b = *( (point*) y );
  return cmp( *a.pval, *b.pval );
}


bool Propagator::end_points_reached(const int i, const int j, const int k )
{
  for( int s=0; s<nb_end_points; ++s ) {
    if( i==((int)end_points_(0,s)) && j==((int)end_points_(1,s)) && k==((int)end_points_(2,s)) ) {
      return true;
    }
  }
  return false;
}


void Propagator::perform_front_propagation_3d( T_callback_intert_node callback_insert_node ) 
{
  // create the Fibonacci heap
  struct fibheap* open_heap = fh_makeheap();
  fh_setcmp(open_heap, compare_points);

  double h = 1.0/n;

  // initialize points
  for( unsigned int i=0; i<n; ++i )
  for( unsigned int j=0; j<p; ++j )
  for( unsigned int k=0; k<q; ++k ) {
    D_(i,j,k) = GW_INFINITE;
    S_(i,j,k) = kFar;
    Q_(i,j,k) = -1;
  }

  // record all the points
  heap_pool = new fibheap_el*[n*p*q];
  memset( heap_pool, 0, n*p*q*sizeof(fibheap_el*) );

  // initalize open list
  point_list existing_points;
  for( int s=0; s<nb_start_points; ++s ) {
    int i = (int) start_points_(0,s);
    int j = (int) start_points_(1,s);
    int k = (int) start_points_(2,s);

    if( D_( i,j,k )==0 ) throw std::runtime_error("start_points should not contain duplicates");

    point* pt = new point( i,j,k );
    pt->pval = &D_( i,j,k );
    existing_points.push_back( pt );      // for deleting at the end
    heap_pool_(i,j,k) = fh_insert( open_heap, pt );     // add to heap
    if( values==NULL ) {
      D_( i,j,k ) = 0;
    } else {
      D_( i,j,k ) = values[s];      
    }
    S_( i,j,k ) = kOpen;
    Q_( i,j,k ) = s;
  }

  // perform the front propagation
  int num_iter = 0;
  bool stop_iteration = GW_False;
  while( !fh_isempty(open_heap) && num_iter<nb_iter_max && !stop_iteration ) {
    num_iter++;

    // remove from open list and set up state to dead
    point& cur_point = * ((point*) fh_extractmin( open_heap )); // current point
    int i = cur_point.i;
    int j = cur_point.j;
    int k = cur_point.k;
    heap_pool_(i,j,k) = NULL;
    S_(i,j,k) = kDead;
    stop_iteration = end_points_reached(i,j,k);

    // recurse on each neighbor
    int nei_i[6] = {i+1,i,i-1,i,i,i};
    int nei_j[6] = {j,j+1,j,j-1,j,j};
    int nei_k[6] = {k,k,k,k,k-1,k+1};
    for( int s=0; s<6; ++s ) {
      int ii = nei_i[s];
      int jj = nei_j[s];
      int kk = nei_k[s];
      
      bool bInsert = true;
      if( callback_insert_node!=NULL ) {
        bInsert = callback_insert_node(i,j,k,ii,jj,kk);
      }
        
      if( ii>=0 && jj>=0 && ii<(int)n && jj<(int)p && kk>=0 && kk<(int)q && bInsert ) {
        double P = h/W_(ii,jj,kk);
        // compute its neighboring values
        double a1 = GW_INFINITE;
        if( ii<(int)n-1 )
          a1 = D_(ii+1,jj,kk);
        if( ii>0 )
          a1 = GW_MIN( a1, D_(ii-1,jj,kk) );
        double a2 = GW_INFINITE;
        if( jj<(int)p-1 )
          a2 = D_(ii,jj+1,kk);
        if( jj>0 )
          a2 = GW_MIN( a2, D_(ii,jj-1,kk) );
        double a3 = GW_INFINITE;
        if( kk<(int)q-1 )
          a3 = D_(ii,jj,kk+1);
        if( kk>0 )
          a3 = GW_MIN( a3, D_(ii,jj,kk-1) );
        // order so that a1<a2<a3
        double tmp = 0;
        #define SWAP(a,b) tmp = a; a = b; b = tmp
        #define SWAPIF(a,b) if(a>b) { SWAP(a,b); }
        SWAPIF(a2,a3)
        SWAPIF(a1,a2)
        SWAPIF(a2,a3)
        // update its distance
        // now the equation is   (a-a1)^2+(a-a2)^2+(a-a3)^2 - P^2 = 0, with a >= a3 >= a2 >= a1.
        // =>    3*a^2 - 2*(a2+a1+a3)*a - P^2 + a1^2 + a3^2 + a2^2
        // => delta = (a2+a1+a3)^2 - 3*(a1^2 + a3^2 + a2^2 - P^2)
        double delta = (a2+a1+a3)*(a2+a1+a3) - 3*(a1*a1 + a2*a2 + a3*a3 - P*P);
        double A1 = 0;
        if( delta>=0 ) {
          A1 = ( a2+a1+a3 + sqrt(delta) )/3.0;
        }
        if( A1<=a3 ) {
          // at least a3 is too large, so we have
          // a >= a2 >= a1  and  a<a3 so the equation is 
          //    (a-a1)^2+(a-a2)^2 - P^2 = 0
          //=> 2*a^2 - 2*(a1+a2)*a + a1^2+a2^2-P^2
          // delta = (a2+a1)^2 - 2*(a1^2 + a2^2 - P^2)
          delta = (a2+a1)*(a2+a1) - 2*(a1*a1 + a2*a2 - P*P);
          A1 = 0;
          if( delta>=0 )
            A1 = 0.5 * ( a2+a1 +sqrt(delta) );
          if( A1<=a2 )
            A1 = a1 + P;
        }
        // update the value
        if( ((int) S_(ii,jj,kk)) == kDead ) {
          // check if action has change. Should not appen for FM
          // if( A1<D_(ii,jj,kk) )
          if( A1<D_(ii,jj,kk) )  { // should not happen for FM
            D_(ii,jj,kk) = A1;
            Q_(ii,jj,kk) = Q_(i,j,k);
          }
        }
        else if( ((int) S_(ii,jj,kk)) == kOpen ) {
          // check if action has change.
          if( A1<D_(ii,jj,kk) ) {
            D_(ii,jj,kk) = A1;
            Q_(ii,jj,kk) = Q_(i,j,k);
            // Modify the value in the heap
            fibheap_el* cur_el = heap_pool_(ii,jj,kk);
            if( cur_el!=NULL ) {
              fh_replacedata( open_heap, cur_el, cur_el->fhe_data );  // use same data for update
            } else {
              throw std::runtime_error("error in heap pool allocation");
            }
          }
        }
        else if( ((int) S_(ii,jj,kk)) == kFar ) {
          if( D_(ii,jj,kk)!=GW_INFINITE ) {
            std::cerr << "warning: distance must be initialized to Inf." << std::endl;
          }
          if( L==NULL || A1<=L_(ii,jj,kk) ) {
            S_(ii,jj,kk) = kOpen;
            // distance must have change.
            D_(ii,jj,kk) = A1;
            Q_(ii,jj,kk) = Q_(i,j,k);
            // add to open list
            point* pt = new point(ii,jj,kk);
            pt->pval = &D_( ii,jj,kk );
            existing_points.push_back( pt );
            heap_pool_(ii,jj,kk) = fh_insert( open_heap, pt );      // add to heap  
          }
        }
        else {
          std::cerr << "warning: unknown state." << std::endl;
        }
      }  // end swich
    }  // end for
  }  // end while

  // free heap
  fh_deleteheap(open_heap);

  // free point pool
  for( point_list::iterator it = existing_points.begin(); it!=existing_points.end(); ++it ) {
    GW_DELETE( *it );
  }

  // free fibheap pool
  GW_DELETEARRAY(heap_pool);

  return;
}
