#ifndef MINCUTILS_H
#define MINCUTILS_H

#include <stdio.h>
#include "mri.h"

/*-----------------------------------------------------*/
typedef struct {
  int   VoxAxisId;     /*   0,   1,     2 */
  char *VoxAxisName;   /* col, row, slice */
  int   MINCAxisId;    /*   0,       1      2    */
  char *MINCAxisName;  /* xspace, yspace, zspace */
  int   Len;           /* Length of axis */
  float Res;           /* Resolution */
  float DirCos[3];     /* Direction Cosine*/
} MINCAXIS;

/*-----------------------------------------------------*/
typedef struct {
  int VoxAxisStorageOrder[3];
  float VolCenterVox[3];
  float VolCenterWorld[3];
  MINCAXIS Axis[3];
} MINCAXES;

/*-----------------------------------------------------*/
int DumpMINCAxes(FILE *fp, MINCAXES *MA);
MINCAXES *ConfigMINCAxes(MRI *mri);
int NameMINCAxes(MINCAXES *MA);
int MINCAxesStorageOrder(MINCAXES *MA);















#endif /* #ifndef MINCUTILS_H */
