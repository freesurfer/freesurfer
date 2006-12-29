/**
 * @file  mri_transform.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


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
