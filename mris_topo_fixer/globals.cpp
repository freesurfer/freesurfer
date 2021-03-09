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


#include "globals.h"

#include <iostream>
using namespace std;

#include "error.h"
#include "macros.h"

#ifdef Windows_NT
#include "proto.h"
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

