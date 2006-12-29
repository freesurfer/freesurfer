/**
 * @file  mri_copy_params.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:07 $
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
// mri_copy_params.cpp
//
// purpose: copy volume parameters from template and write out the volume
//
// written by: tosa@nmr.mgh.harvard.edu
// date      : Feb. 23rd, 2005
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2006/12/29 02:09:07 $
// Revision       : $Revision: 1.2 $

#include <iostream>
#include <iomanip>

extern "C" {
#include "error.h"
#include "mri.h"
#include "version.h"
  char *Progname = "mri_copy_params";
}

using namespace std;

void print_usage() {
  cout << "Usage: mri_copy_params <in_vol> <template_vol> <out_vol>" << endl;
  cout << "     : where all volume parameters of in_vol are replaced with those of template_vol." << endl;
}

int main(int argc, char *argv[]) {
  bool bVolumeDifferent = false;
  bool bSizeDifferent = false;
  int nargs;
  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_copy_params.cpp,v 1.2 2006/12/29 02:09:07 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  if (argc < 3) {
    print_usage();
    return -1;
  }

  MRI *in = MRIread(argv[1]);
  if (!in) {
    cerr << "could not open " << argv[1] << endl;
    return -1;
  }
  MRI *temp = MRIreadHeader(argv[2], MRI_UNDEFINED);
  if (!temp) {
    cerr << "could not open " << argv[2] << endl;
    return -1;
  }
  MRI *dst = MRIcopy(in, NULL);

  // check few things
  if ((temp->width != in->width)
      || (temp->height != in->height)
      || (temp->depth != in->depth)) {
    cerr << "WARNING: volume sizes are different" << endl;
    cerr << "    in_vol : " << in->width << ", " << in->height << ", " << in->depth << endl;
    cerr << "  temp_vol : " << temp->width << ", " << temp->height << ", " << temp->depth << endl;
    bVolumeDifferent = true;
  }
  if ((temp->xsize != in->xsize)
      || (temp->ysize != in->ysize)
      || (temp->zsize != in->zsize)) {
    cerr << "WARNING: voxel sizes are different" << endl;
    cerr << "    in_vol : " << in->xsize << ", " << in->ysize << ", " << in->zsize << endl;
    cerr << "  temp_vol : " << temp->xsize << ", " << temp->ysize << ", " << temp->zsize << endl;
    bSizeDifferent= true;
  }
  // copy everything in the header from template
  MRIcopyHeader(temp, dst);
  // just few things restored
  if (bVolumeDifferent) {
    dst->width = in->width;
    dst->height = in->height;
    dst->depth = in->depth;
  }
  if (bSizeDifferent);
  {
    dst->xsize = in->xsize;
    dst->ysize = in->ysize;
    dst->zsize = in->zsize;
  }
  //
  MRIwrite(dst, argv[3]);

  MRIfree(&in);
  MRIfree(&temp);
  MRIfree(&dst);

  return (NO_ERROR);
}
