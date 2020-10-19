/*
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


#ifndef REGION_H
#define REGION_H

#include <stdio.h>
#include "mri.h"


MRI_REGION *REGIONsubtract(MRI_REGION *reg1, MRI_REGION *reg2,
                           MRI_REGION *rdst) ;
MRI_REGION *REGIONadd(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst);
MRI_REGION *REGIONclear(MRI_REGION *r) ;
MRI_REGION *REGIONintersect(MRI_REGION *reg1, MRI_REGION *reg2,
                            MRI_REGION *rdst) ;
MRI_REGION *REGIONunion(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst);
MRI_REGION *REGIONalloc(void) ;
MRI_REGION *REGIONcopy(MRI_REGION *rsrc, MRI_REGION *rdst) ;
int        REGIONinside(MRI_REGION *reg, int x, int y, int z) ;
MRI_REGION *REGIONexpand(MRI_REGION *rsrc, MRI_REGION *rdst, int n) ;
float      REGIONminCornerDistance(MRI_REGION *r1, MRI_REGION *r2) ;
int REGIONprint(FILE *fp, MRI_REGION *r);
MRI_REGION *REGIONgetBoundingBox(MRI *mask, int npad);
MRI_REGION *REGIONgetBoundingBoxM(const MRI *mask, const int npad[6]);
MRI_REGION *REGIONgetBoundingBoxEqOdd(MRI *mask, int npad);

#define REGION_INSIDE      1
#define REGION_ON_BORDER   -1
#define REGION_OUTSIDE     0

#endif
