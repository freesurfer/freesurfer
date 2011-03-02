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
 *    $Date: 2011/03/02 00:04:56 $
 *    $Revision: 1.5 $
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
#ifdef Windows_NT
#include "proto.h"
#endif
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

