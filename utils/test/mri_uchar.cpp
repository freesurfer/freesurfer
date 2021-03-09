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

  const char *Progname = "mri_uchar";
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
