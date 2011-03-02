/**
 * @file  mincutils.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
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


#include <stdio.h>
#include "mri.h"
#include "mincutils.h"

/*-----------------------------------------------------*/
int DumpMINCAxes(FILE *fp, MINCAXES *MA)
{
  int i;

  fprintf(fp,"Storage Order: %d %d %d\n",
          MA->VoxAxisStorageOrder[0],
          MA->VoxAxisStorageOrder[1],
          MA->VoxAxisStorageOrder[2]);
  fprintf(fp,"Volume Voxel Center: %g %g %g\n",
          MA->VolCenterVox[0],
          MA->VolCenterVox[1],
          MA->VolCenterVox[2]);
  fprintf(fp,"Volume World Center: %g %g %g\n",
          MA->VolCenterWorld[0],
          MA->VolCenterWorld[1],
          MA->VolCenterWorld[2]);

  for (i=0;i<3;i++)
  {
    fprintf(fp,"VoxAxisName:  %s ------ \n",MA->Axis[i].VoxAxisName);
    fprintf(fp,"VoxAxisId:    %d\n",MA->Axis[i].VoxAxisId);
    fprintf(fp,"MINCAxisId:   %d\n",MA->Axis[i].MINCAxisId);
    fprintf(fp,"MINCAxisName: %s\n",MA->Axis[i].MINCAxisName);
    fprintf(fp,"Length:       %d\n",MA->Axis[i].Len);
    fprintf(fp,"Resolution:   %g\n",MA->Axis[i].Res);
    fprintf(fp,"DC:           %g %g %g \n",MA->Axis[i].DirCos[0],
            MA->Axis[i].DirCos[1],MA->Axis[i].DirCos[2]);
  }

  return(0);
}

/*-----------------------------------------------------*/
MINCAXES *ConfigMINCAxes(MRI *mri)
{
  MINCAXES *MA;
  int va;

  MA = (MINCAXES *) calloc(1,sizeof(MINCAXES));

  /* column voxel axis */
  va = 0;
  MA->Axis[va].VoxAxisId = va;
  MA->Axis[va].VoxAxisName = "col";
  MA->Axis[va].VoxAxisId = va;
  MA->Axis[va].Len = mri->width;
  MA->Axis[va].Res = mri->xsize;
  MA->Axis[va].DirCos[0] = mri->x_r;
  MA->Axis[va].DirCos[1] = mri->x_a;
  MA->Axis[va].DirCos[2] = mri->x_s;

  /* row voxel axis */
  va = 1;
  MA->Axis[va].VoxAxisId = va;
  MA->Axis[va].VoxAxisName = "row";
  MA->Axis[va].VoxAxisId = va;
  MA->Axis[va].Len = mri->height;
  MA->Axis[va].Res = mri->ysize;
  MA->Axis[va].DirCos[0] = mri->y_r;
  MA->Axis[va].DirCos[1] = mri->y_a;
  MA->Axis[va].DirCos[2] = mri->y_s;

  /* slice voxel axis */
  va = 2;
  MA->Axis[va].VoxAxisId = va;
  MA->Axis[va].VoxAxisName = "slice";
  MA->Axis[va].VoxAxisId = va;
  MA->Axis[va].Len = mri->depth;
  MA->Axis[va].Res = mri->zsize;
  MA->Axis[va].DirCos[0] = mri->z_r;
  MA->Axis[va].DirCos[1] = mri->z_a;
  MA->Axis[va].DirCos[2] = mri->z_s;

  MA->VolCenterWorld[0] = mri->c_r;
  MA->VolCenterWorld[1] = mri->c_a;
  MA->VolCenterWorld[2] = mri->c_s;

  NameMINCAxes(MA);
  MINCAxesStorageOrder(MA);

  return(MA);
}

/*-----------------------------------------------------*/
int NameMINCAxes(MINCAXES *MA)
{
  int xspacehit, yspacehit, zspacehit;
  int xspaceid, yspaceid, zspaceid;
  char *space[3];
  int   colspaceid=0, rowspaceid=0, slcspaceid=0;
  float col_dc_x, col_dc_y, col_dc_z;
  float row_dc_x, row_dc_y, row_dc_z;
  float slc_dc_x, slc_dc_y, slc_dc_z;
  int err;

  col_dc_x = fabs(MA->Axis[0].DirCos[0]);
  col_dc_y = fabs(MA->Axis[0].DirCos[1]);
  col_dc_z = fabs(MA->Axis[0].DirCos[2]);

  row_dc_x = fabs(MA->Axis[1].DirCos[0]);
  row_dc_y = fabs(MA->Axis[1].DirCos[1]);
  row_dc_z = fabs(MA->Axis[1].DirCos[2]);

  slc_dc_x = fabs(MA->Axis[2].DirCos[0]);
  slc_dc_y = fabs(MA->Axis[2].DirCos[1]);
  slc_dc_z = fabs(MA->Axis[2].DirCos[2]);

  xspaceid = 0;
  yspaceid = 1;
  zspaceid = 2;

  xspacehit = 0;
  yspacehit = 0;
  zspacehit = 0;

  /* Name the Column Axis */
  if (col_dc_x >= row_dc_x && col_dc_x >= slc_dc_x)
  {
    colspaceid = xspaceid;
    MA->VolCenterVox[0] = (MA->Axis[0].Len-1)/2;
    xspacehit = 1;
  }
  else if (col_dc_y >= row_dc_y && col_dc_y >= slc_dc_y)
  {
    colspaceid = yspaceid;
    MA->VolCenterVox[1] = (MA->Axis[0].Len-1)/2;
    yspacehit = 1;
  }
  else if (col_dc_z >= row_dc_z && col_dc_z >= slc_dc_z)
  {
    colspaceid = zspaceid;
    MA->VolCenterVox[2] = (MA->Axis[0].Len-1)/2;
    zspacehit = 1;
  }

  /* Name the Row Axis */
  if (!xspacehit && row_dc_x >= slc_dc_x)
  {
    rowspaceid = xspaceid;
    MA->VolCenterVox[0] = (MA->Axis[1].Len-1)/2;
    xspacehit = 1;
  }
  else if (!yspacehit && row_dc_y >= slc_dc_y)
  {
    rowspaceid = yspaceid;
    MA->VolCenterVox[1] = (MA->Axis[1].Len-1)/2;
    yspacehit = 1;
  }
  else if (!zspacehit && row_dc_z >= slc_dc_z)
  {
    rowspaceid = zspaceid;
    MA->VolCenterVox[2] = (MA->Axis[1].Len-1)/2;
    zspacehit = 1;
  }

  /* Name the Slice Axis */
  if (!xspacehit)
  {
    slcspaceid = xspaceid;
    MA->VolCenterVox[0] = (MA->Axis[2].Len-1)/2;
    xspacehit = 1;
  }
  else if (!yspacehit)
  {
    slcspaceid = yspaceid;
    MA->VolCenterVox[1] = (MA->Axis[2].Len-1)/2;
    yspacehit = 1;
  }
  if (!zspacehit)
  {
    slcspaceid = zspaceid;
    MA->VolCenterVox[2] = (MA->Axis[2].Len-1)/2;
    zspacehit = 1;
  }

  /* Check for errors */
  err = 0;
  if (!xspacehit)
  {
    printf("ERROR: could not assign xspace\n");
    err = 1;
  }
  if (!yspacehit)
  {
    printf("ERROR: could not assign yspace\n");
    err = 1;
  }
  if (!zspacehit)
  {
    printf("ERROR: could not assign zspace\n");
    err = 1;
  }
  if (err) return(1);

  /* Make final assignments */
  space[xspaceid] = MIxspace;
  space[yspaceid] = MIyspace;
  space[zspaceid] = MIzspace;

  MA->Axis[0].MINCAxisId   = colspaceid;
  MA->Axis[1].MINCAxisId   = rowspaceid;
  MA->Axis[2].MINCAxisId   = slcspaceid;

  MA->Axis[0].MINCAxisName = space[colspaceid];
  MA->Axis[1].MINCAxisName = space[rowspaceid];
  MA->Axis[2].MINCAxisName = space[slcspaceid];

  return(0);
}
/*-----------------------------------------------------*/
int MINCAxesStorageOrder(MINCAXES *MA)
{
  int ncols, nrows, nslcs;
  int col, row, slc;

  col = 0;
  row = 1;
  slc = 2;

  ncols = MA->Axis[0].Len;
  nrows = MA->Axis[1].Len;
  nslcs = MA->Axis[2].Len;

  if (ncols >= nrows && ncols >= nslcs)
  {
    MA->VoxAxisStorageOrder[0] = col;
    if (nrows >= nslcs)
    {
      MA->VoxAxisStorageOrder[1] = row;
      MA->VoxAxisStorageOrder[2] = slc;
    }
    else
    {
      MA->VoxAxisStorageOrder[1] = slc;
      MA->VoxAxisStorageOrder[2] = row;
    }
  }
  else if (nrows >= nslcs)
  {
    MA->VoxAxisStorageOrder[0] = row;
    if (ncols >= nslcs)
    {
      MA->VoxAxisStorageOrder[1] = col;
      MA->VoxAxisStorageOrder[2] = slc;
    }
    else
    {
      MA->VoxAxisStorageOrder[1] = slc;
      MA->VoxAxisStorageOrder[2] = col;
    }
  }
  else
  {
    MA->VoxAxisStorageOrder[0] = slc;
    if (ncols >= nrows)
    {
      MA->VoxAxisStorageOrder[1] = col;
      MA->VoxAxisStorageOrder[2] = row;
    }
    else
    {
      MA->VoxAxisStorageOrder[1] = row;
      MA->VoxAxisStorageOrder[2] = col;
    }
  }

  return(0);
}


