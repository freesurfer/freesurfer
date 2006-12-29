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
 *    $Date: 2006/12/29 02:09:17 $
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

