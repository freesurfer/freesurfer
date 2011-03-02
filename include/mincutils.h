/**
 * @file  mincutils.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef MINCUTILS_H
#define MINCUTILS_H

#include <stdio.h>
#include "mri.h"

/*-----------------------------------------------------*/
typedef struct
{
  int   VoxAxisId;     /*   0,   1,     2 */
  char *VoxAxisName;   /* col, row, slice */
  int   MINCAxisId;    /*   0,       1      2    */
  char *MINCAxisName;  /* xspace, yspace, zspace */
  int   Len;           /* Length of axis */
  float Res;           /* Resolution */
  float DirCos[3];     /* Direction Cosine*/
}
MINCAXIS;

/*-----------------------------------------------------*/
typedef struct
{
  int VoxAxisStorageOrder[3];
  float VolCenterVox[3];
  float VolCenterWorld[3];
  MINCAXIS Axis[3];
}
MINCAXES;

/*-----------------------------------------------------*/
int DumpMINCAxes(FILE *fp, MINCAXES *MA);
MINCAXES *ConfigMINCAxes(MRI *mri);
int NameMINCAxes(MINCAXES *MA);
int MINCAxesStorageOrder(MINCAXES *MA);















#endif /* #ifndef MINCUTILS_H */
