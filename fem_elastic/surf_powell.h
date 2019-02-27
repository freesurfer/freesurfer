
#ifndef H_LIN_POWELL_H
#define H_LIN_POWELL_H

#include "surf_types.h"



// both /usr/pubsw/packages/petsc/current/include/petscvec.h and mrisurf.h
// define NORM_MAX, so undef it to reveal who uses it, if anybody.
#undef NORM_MAX
#include "mrisurf.h"
#undef NORM_MAX


//
// the energy defined will be the MAX of the norm of the error
//
// the purpose is to come close to an optimal linear registration
// and thus to deal with small displacements
typedef float (*Energy_pointer)(const PointsContainerType&,
                                float* transform);


// the surface structure will include 2 sets of points
//
// it will try to map the vertex x,y,z
// to the origx,y,z + vertex x,y,z values
//
// the transform structure will be a 4x3 matrix
// which will be FULLY LINEAR
float
powell_minimize( PointsContainerType&,
                 float *iotransform,
                 Energy_pointer fun_pointer);

#endif
