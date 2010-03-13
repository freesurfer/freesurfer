/**
 * @file  mri_transform.c
 * @brief transform utils
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/03/13 01:32:44 $
 *    $Revision: 1.3 $
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

#include "mri_transform.h"

static float mfStartX, mfEndX;
static float mfStartY, mfEndY;
static float mfStartZ, mfEndZ;
static float mfSizeX, mfSizeY, mfSizeZ;

void trans_SetBounds ( float ifStartX, float ifEndX,
                       float ifStartY, float ifEndY,
                       float ifStartZ, float ifEndZ )
{
  mfStartX = ifStartX;
  mfEndX   = ifEndX;
  mfStartY = ifStartY;
  mfEndY   = ifEndY;
  mfStartZ = ifStartZ;
  mfEndZ   = ifEndZ;
}

void trans_SetResolution ( float ifSizeX, float ifSizeY, float ifSizeZ )
{
  mfSizeX = ifSizeX;
  mfSizeY = ifSizeY;
  mfSizeZ = ifSizeZ;
}

void trans_RASToVoxelIndex ( double irRASX, double irRASY, double irRASZ,
                             int *onVoxX, int *onVoxY, int *onVoxZ )
{
  double orVoxX, orVoxY, orVoxZ ;

  trans_RASToVoxel(irRASX, irRASY, irRASZ, &orVoxX, &orVoxY, &orVoxZ) ;
  *onVoxX = (int)orVoxX ;
  *onVoxY = (int)orVoxY ;
  *onVoxZ = (int)orVoxZ ;
}

void trans_RASToVoxel ( double irRASX, double irRASY, double irRASZ,
                        double *orVoxX, double *orVoxY, double *orVoxZ )
{
  /* we need to stay in the same type so
     as not to typecast thee conversion to
     int until right at the very end. */
  float fRASX, fRASY, fRASZ;
  fRASX = (float) irRASX;
  fRASY = (float) irRASY;
  fRASZ = (float) irRASZ;

  *orVoxX = ( (  mfEndX - fRASX    ) / mfSizeX );
  *orVoxY = ( ( -fRASZ  - mfStartY ) / mfSizeY );
  *orVoxZ = ( (  fRASY  - mfStartZ ) / mfSizeZ );
}

void trans_VoxelToRAS ( double irVoxX, double  irVoxY, double  irVoxZ,
                        double *orRASX, double *orRASY, double *orRASZ )
{
  *orRASX = (double) (     mfEndX   - ( irVoxX * mfSizeX )   );
  *orRASY = (double) (     mfStartZ + ( irVoxZ * mfSizeZ )   );
  *orRASZ = (double) ( - ( mfStartY + ( irVoxY * mfSizeY ) ) );
}

void trans_VoxelIndexToRAS ( int inVoxX, int  inVoxY, int  inVoxZ,
                             double *orRASX, double *orRASY, double *orRASZ )
{
  double irVoxX, irVoxY, irVoxZ ;

  irVoxX = (double)inVoxX ;
  irVoxY = (double)inVoxY ;
  irVoxZ = (double)inVoxZ ;
  trans_VoxelToRAS(inVoxX, inVoxY, inVoxZ, orRASX, orRASY, orRASZ) ;
}
