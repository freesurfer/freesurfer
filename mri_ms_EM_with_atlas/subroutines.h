/**
 * @file  subroutines.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:07 $
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
