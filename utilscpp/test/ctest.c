/**
 * @file  ctest.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:56 $
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


//
// ctest.c
//
// purpose: usability of cpp library exported C function
//
#include <stdio.h>
#include <stdlib.h>
#include "mri.h"
#include "fastmarching.h"

char *Progname = "ctest";

int main(int argc, char *argv[]) {
  MRI *src;
  MRI *dst;
  int label=2;
  float max_distance=10.;
  int mode=1;

  if (argc < 1) {
    fprintf(stderr, "Usage: ctest <volume>\n");
    return -1;
  }
  src = MRIread(argv[1]);
  dst = MRIextractDistanceMap(src, NULL, label, max_distance, mode,NULL);
  return 0;
}

