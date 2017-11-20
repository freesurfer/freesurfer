#ifndef GEOS_H
#define GEOS_H

#include "mri.h"
#include "voxlist.h"

typedef struct
{
  double  lambda ;
  int     wsize ;
  double  max_dist ;
} GEOS_PARMS, GP ;


VOXEL_LIST *GEOSproposeSegmentation(MRI *image, 
				    MRI *mri_mask,
				    MRI_REGION *box,
				    VOXEL_LIST *vl_interior, 
				    VOXEL_LIST *vl_exterior,
				    GEOS_PARMS *parms);


#endif
