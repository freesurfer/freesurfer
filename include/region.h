/**
 * @file  region.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.5 $
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


#define REGION_INSIDE      1
#define REGION_ON_BORDER   -1
#define REGION_OUTSIDE     0

#endif
