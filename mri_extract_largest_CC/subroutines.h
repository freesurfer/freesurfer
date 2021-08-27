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


#ifndef SUBROUTINES_H
#define SUBROUTINES_H

#include "mri.h"

typedef struct {
  int x,y,z;
}
POINTI;


void RemoveHoles(MRI *orivol);
void  GrassFire(MRI *orivol, MRI *Label, int label, POINTI *Pt,
                int *curSize, int *minX, int *maxX, int *minY, int *maxY,
                int *minZ, int *maxZ);
void  GrassFire6(MRI *orivol, MRI *Label, int label, POINTI *Pt,
                 int *curSize, int *minX, int *maxX, int *minY, int *maxY,
                 int *minZ, int *maxZ);

void  GrassFire18(MRI *orivol, MRI *Label, int label, POINTI *Pt,
                  int *curSize, int *minX, int *maxX, int *minY, int *maxY,
                  int *minZ, int *maxZ);
void GetLargestCC18(MRI *orivol);
void GetLargestCC6(MRI *orivol);
MRI * Dilation6(MRI *ori, MRI *out, int);
MRI * Erosion6(MRI *ori, MRI *out, int);
MRI * Dilation26(MRI *ori, MRI *out, int);
MRI * Erosion26(MRI *ori, MRI *out, int);
MRI * BinaryOpen6(MRI *ori, MRI *out, int R);
MRI * BinaryOpen26(MRI *ori, MRI *out, int R);
MRI * BinaryClose6(MRI *ori, MRI *out, int R);
MRI * BinaryClose26(MRI *ori, MRI *out, int R);


#endif
