/**
 * @file  globals.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:19 $
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


#include "topology/globals.h"

#include <iostream>
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
  //#include <stdio.h>
  //#include <stdlib.h>
  //#include <string.h>
  //#include <math.h>
  //#include <ctype.h>
  //#include "macros.h"
#include "error.h"
  //#include "tags.h"
  //#include "diag.h"
  //#include "proto.h"
  //#include "timer.h"
  //#include "mrisurf.h"
  //#include "mri.h"
#include "macros.h"
  //#include "icosahedron.h"
  //#include "mrishash.h"
  //#include "version.h"
#ifdef __cplusplus
}
#endif


void check(bool exp) {
  if (exp==false) cout << "e";
}

void ErrorExit(string s) {
  cout << endl << "ERROR: " << s << endl;
  exit(-1);
}

int Random(int nmax) {
  //  return rand()*nmax/RAND_MAX;
  return nint(randomNumber(0.0, (double)nmax-1));
}

