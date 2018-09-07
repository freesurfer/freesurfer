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
 *    $Date: 2011/03/02 00:04:23 $
 *    $Revision: 1.4 $
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

const char *Progname;

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
