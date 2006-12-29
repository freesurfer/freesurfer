/**
 * @file  mri_transform.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:35 $
 *    $Revision: 1.2 $
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

void trans_RASToVoxelIndex ( Real irRASX, Real irRASY, Real irRASZ,
                             int *onVoxX, int *onVoxY, int *onVoxZ )
{
  Real orVoxX, orVoxY, orVoxZ ;

  trans_RASToVoxel(irRASX, irRASY, irRASZ, &orVoxX, &orVoxY, &orVoxZ) ;
  *onVoxX = (int)orVoxX ;
  *onVoxY = (int)orVoxY ;
  *onVoxZ = (int)orVoxZ ;
}

void trans_RASToVoxel ( Real irRASX, Real irRASY, Real irRASZ,
                        Real *orVoxX, Real *orVoxY, Real *orVoxZ )
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

void trans_VoxelToRAS ( Real irVoxX, Real  irVoxY, Real  irVoxZ,
                        Real *orRASX, Real *orRASY, Real *orRASZ )
{


  *orRASX = (Real) (     mfEndX   - ( irVoxX * mfSizeX )   );
  *orRASY = (Real) (     mfStartZ + ( irVoxZ * mfSizeZ )   );
  *orRASZ = (Real) ( - ( mfStartY + ( irVoxY * mfSizeY ) ) );

}

void trans_VoxelIndexToRAS ( int inVoxX, int  inVoxY, int  inVoxZ,
                             Real *orRASX, Real *orRASY, Real *orRASZ )
{

  Real irVoxX, irVoxY, irVoxZ ;

  irVoxX = (Real)inVoxX ;
  irVoxY = (Real)inVoxY ;
  irVoxZ = (Real)inVoxZ ;
  trans_VoxelToRAS(inVoxX, inVoxY, inVoxZ, orRASX, orRASY, orRASZ) ;
}
