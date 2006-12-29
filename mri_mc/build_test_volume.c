/**
 * @file  build_test_volume.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:07 $
 *    $Revision: 1.3 $
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


#include <stdio.h>
#include "mri.h"

char *Progname;

int main(int argc, char *argv[]) {
  MRI *mri;
  Progname=argv[0];
  printf("Generating test_volume.mgz...\n");

  mri=MRIalloc(20,20,20,MRI_UCHAR);

  MRIvox(mri,10,10,10)=255;
  MRIvox(mri,11,10,10)=255;
  MRIvox(mri,11,11,10)=255;
  MRIvox(mri,11,11,11)=255;
  MRIvox(mri,10,11,11)=255;
  MRIvox(mri,10,10,11)=255;

  MRIwrite(mri,"test_volume.mgz");

  return 0;
}
