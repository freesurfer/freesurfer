#ifndef mri_transform_h
#define mri_transform_h

#include "volume_io/basic.h" /* defines Real */


                                   /* need to set the size of the volume
                                      and its resolution */
void trans_SetBounds ( float ifStartX, float ifEndX,
                    float ifStartY, float ifEndY,
                    float ifStartZ, float ifEndZ );

void trans_SetResolution ( float ifSizeX, float ifSizeY, float ifSizeZ );

                                  /* converts ras coords to
                                     voxel coords and back */
void trans_RASToVoxel
          ( Real irX, Real irY, Real irZ,    /* incoming ras coords */
            Real *onX, Real *onY, Real *onZ );  /* outgoing voxel coords */
void trans_VoxelToRAS
          ( Real   inX, Real   inY, Real   inZ,   /* incoming voxel coords */
            Real *orX, Real *orY, Real *orZ ); /* outgoing RAS coords */
void trans_RASToVoxelIndex
          ( Real irX, Real irY, Real irZ,    /* incoming ras coords */
            int *onX, int *onY, int *onZ );  /* outgoing voxel coords */
void trans_VoxelIndexToRAS ( int inVoxX, int  inVoxY, int  inVoxZ,  
                        Real *orRASX, Real *orRASY, Real *orRASZ ) ;


#endif
