/**
 * @file  mri_uchar.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:46 $
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
// mri_uchar.cpp
//
// purpose: convert data to uchar value without scaling
//

#include <iostream>
#include <iomanip>

extern "C"
{
#include "mri.h"

  char *Progname = "mri_uchar";
}

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 1)
  {
    cerr << "Usage: mri_uchar <involume> <outvolume>" << endl;
    cerr << "       any involume val is set to uchar volume " << endl;
    cerr << "       i.e. -0.5 <= val < 0.5 becomes 0. " << endl;
    cerr << "       make sure that involume does not have values more than 255 " << endl;
    return -1;
  }
  MRI *src = MRIread(argv[1]);
  if (!src)
  {
    cerr << "could not open " << argv[1] << endl;
    return -1;
  }
  MRI *dst = MRIalloc(src->width, src->height, src->depth, MRI_UCHAR);
  if (!dst)
  {
    cerr << "could not allocate memory for the target" << endl;
    return -1;
  }
  // copy geometry information
  MRIcopyHeader(src, dst);

  int count = 0;
  for (int f=0; f < src->nframes; f++)
    for (int k=0; k < src->depth; k++)
      for (int j=0; j < src->height; j++)
        for (int i=0; i < src->width; i++)
        {
          float val = MRIgetVoxVal(src, i, j, k, f);
          // -0.5 up to 0.5(not including) becomes 0
          float fapp = floorf(val+0.5);
          if (fapp > 0)
            count++;
          MRIsetVoxVal(dst, i, j, k, f, fapp);
        }
  cout << "non-zero value count = " << count << endl;
  MRIwrite(dst, argv[2]);
}
