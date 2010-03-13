/**
 * @file  mri_transform.h
 * @brief transform utils
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/03/13 01:32:40 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef mri_transform_h
#define mri_transform_h

/* need to set the size of the volume
   and its resolution */
void trans_SetBounds ( float ifStartX, float ifEndX,
                       float ifStartY, float ifEndY,
                       float ifStartZ, float ifEndZ );

void trans_SetResolution ( float ifSizeX, float ifSizeY, float ifSizeZ );

/* converts ras coords to
   voxel coords and back */
void trans_RASToVoxel
( double irX, double irY, double irZ,    /* incoming ras coords */
  double *onX, double *onY, double *onZ );  /* outgoing voxel coords */
void trans_VoxelToRAS
( double   inX, double   inY, double   inZ,   /* incoming voxel coords */
  double *orX, double *orY, double *orZ ); /* outgoing RAS coords */
void trans_RASToVoxelIndex
( double irX, double irY, double irZ,    /* incoming ras coords */
  int *onX, int *onY, int *onZ );  /* outgoing voxel coords */
void trans_VoxelIndexToRAS ( int inVoxX, int  inVoxY, int  inVoxZ,
                             double *orRASX, double *orRASY, double *orRASZ ) ;

#endif
