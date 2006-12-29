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
 *    $Date: 2006/12/29 02:09:19 $
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
  dst = MRIextractDistanceMap(src, NULL, label, max_distance, mode);
  return 0;
}

