#include "geos.h"


VOXEL_LIST *
GEOSproposeSegmentation(MRI *image, 
			MRI *mri_mask,
			MRI_REGION *box,
			VOXEL_LIST **interior, int num_interior, 
			VOXEL_LIST **exterior, int num_exterior,
			GEOS_PARMS *parms)
{
  VOXEL_LIST  *vl_dst ;
  
  vl_dst = VLSTalloc(10) ;
  return(vl_dst) ;
}

