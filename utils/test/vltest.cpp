/**
 * @brief test routines
 *
 */
/*
 * Original Author: Y. Tosa
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
