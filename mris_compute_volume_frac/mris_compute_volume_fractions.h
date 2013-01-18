#include "utils.h"
#include "mrisurf.h"


typedef struct {
  double frac;  /* volume fraction */
  double err;  /* error upper bound */
} volFraction;

typedef struct {
  double cent[3]; /* center in mm */
  double corn[8][3]; /* corners in mm */
  double vsize[3]; /* size in mm */
  double vox[3]; /* origin in mm */
} octTreeVoxel;
MRI* MRIcomputeVolumeFractionsFromSurface(MRI_SURFACE *mris, double accuracy, MRI *mri_src, MRI *mri_fractions);
volFraction MRIcomputeVoxelFractions(octTreeVoxel V, VERTEX *v, double acc, int max_depth, MRI_SURFACE *mris);
octTreeVoxel octTreeVoxelCreate (double *vox, double *vsize);
octTreeVoxel octTreeVoxelDivide (int type, octTreeVoxel v);

