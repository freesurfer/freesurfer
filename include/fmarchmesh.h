/* Header function for computing geodesic distance on a triangular 
 * surfmesh
 */

#ifndef FMARCHING_H
#define FMARCHING_H

#include <math.h>
#include "heap.h"
#include "mrisurf.h"

#define ALIVE 1
#define NBAND 2
#define FAWAY 3 
#define INFINITY 1000000000 

typedef struct
{
  float x;
  float y;
  float z;
} MyVector; 

/* function declarations */
float *FastMarchMesh(MRI_SURFACE *mesh, int *contour, int numinitvert, float thred);


#endif

