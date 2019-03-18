
#include "surf_powell.h"



#include "numerics.h"
#include "matrix.h"
;

static float
compute_powell_cost(float *p);


//------------------------------

// global variables
//

PointsContainerType* g_container;
Energy_pointer g_fpointer;
float         *g_transf; // data buffer for the transformation

float
powell_minimize( PointsContainerType& container,
                 float *iotransform,
                 Energy_pointer fun_pointer)
{
  float *p;
  float **xi;
  float fret;
  float fstart;
  int   r,c;
  int   iter;

  // assign global vars
  g_container = &container;
  g_fpointer = fun_pointer;

  int npars = 12;

  g_transf = new float[npars];

  p = vector(1, npars);
  xi = matrix(1,npars, 1, npars);

  for (int i=0; i<npars; i++) p[i+1] = iotransform[i];

  for ( r=1; r<=npars; r++)
    for ( c=1; c<=npars; c++)
      xi[r][c] = (r==c)?1:0;

  fstart = compute_powell_cost(p);

  OpenPowell(p, xi,npars, 1.0e-4, &iter, &fret, compute_powell_cost);
  if ( fret < fstart )
  {
    for (int i=0; i<npars; i++) iotransform[i]=p[i+1];
  }
  // iterate greedy
  do
  {
    // reinit directions
    for ( r=1; r<=npars; r++)
      for ( c=1; c<=npars; c++)
        xi[r][c] = (r==c)?1:0;

    fstart = fret;
    OpenPowell(p, xi, npars, 1.0e-4, &iter, &fret, compute_powell_cost);
    if ( fret < fstart )
    {
      for (int i=0; i<npars; i++) iotransform[i]=p[i+1];
    }
  }
  while ( fret < fstart && fret/255>.10  );

  free_matrix(xi, 1, npars, 1, npars);
  free_vector(p,1,npars);
  delete[] g_transf;

  return fstart;

}


float
compute_powell_cost(float *p)
{
  int npars = 12;
  // copy p in g_transf
  for (int i=0; i<npars; i++) g_transf[i]=p[i+1];

  return g_fpointer(*g_container, g_transf);
}

