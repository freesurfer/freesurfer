/**
 * @file  vltest.cpp
 * @brief test routines
 *
 */
/*
 * Original Author: Y. Tosa
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/12/10 19:26:12 $
 *    $Revision: 1.4 $
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

#include <iostream>
#include <fstream>
#include <mcheck.h>

extern "C"
{
#include "label.h"
#include "vlabels.h"
#include "transform.h"
const char *Progname = "vltest";
}

using namespace std;

int main(int, char **)
{
  VLI *pVli = 0;
  VL ***voxel_labels;
  VOXEL_LABELS *vl;

  mtrace();
  for (size_t i=0; i < 10; ++i)
  {
    pVli = VLread(const_cast<char *>("vltest.dat"));
    voxel_labels = pVli->vl;
    vl = &voxel_labels[0][0][0];
    VLfree(&pVli);
  }
}
