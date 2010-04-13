/*

Gheorghe Postelnicu, 2006

Computes the euclidian distance between two surfaces with
1-2-1 vertex correspondance

used by Powell

*/

#ifndef _h_surf_energy_h_
#define _h_surf_energy_h_

#include "surf_powell.h"

// the contents of the container will not be cleared
void add_to_container(const MRI_SURFACE* mris_x,
                      const MRI_SURFACE* mris_fx,
                      PointsContainerType& container);

float energy(const PointsContainerType& container,
             float *transform);

void apply_lin_transform(MRI_SURFACE* mris, float* transform);

#endif
