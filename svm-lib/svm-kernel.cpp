/**
 * @file  svm-kernel.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
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


#include "svm-kernel.h"

using namespace std;


bool Kernel::parse(const char *paramString) {
  if ( _p != NULL )
    delete _p;

  int type, t;
  sscanf(paramString,"%d", &type);


  switch (type) {
  case LINEAR_KERNEL:
    _p = new LinearKernelParam();
    break;
  case POLY_KERNEL:
    int d;
    sscanf(paramString,"%d %d", &t, &d);
    _p = new PolyKernelParam(d);
    break;
  case RBF_KERNEL:
    double gamma;
    sscanf(paramString,"%d %lf", &t, &gamma);
    _p = new RbfKernelParam(gamma);
    break;
  default:
    cerr << "Kernel error: Unknown kernel type " << type << ".\n";
    _p = NULL;
    return false;
  }

  return true;
}



bool Kernel::read(FILE* f, bool binary) {
  char header[300];
  fgets(header,300,f);

  return parse(header);
}

