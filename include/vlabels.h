/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef VOXEL_LABELS_H
#define VOXEL_LABELS_H


typedef struct
{
  unsigned short nlabels ;
  unsigned char *labels ;
  unsigned short *counts ;
}
VOXEL_LABELS, VL ;

typedef struct
{
  int          width ;
  int          height ;
  int          depth ;
  float        resolution ;
  VOXEL_LABELS ***vl ;
}
VOXEL_LABELS_IMAGE, VLI ;

VOXEL_LABELS_IMAGE  *VLalloc(int width, int height,int depth,float resolution);
int                 VLfree(VLI **pvli) ;
int                 VLwrite(VLI *vli, char *fname) ;
VLI                 *VLread(char *fname) ;
VL                  *VLreadVoxel(char *fname, int x, int y, int z,  VL *vl) ;
int                 VLnormalize(VLI *vli) ;


#define VL_MAGIC 0xaefcdae

#endif
