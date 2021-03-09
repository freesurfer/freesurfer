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
// difftool.cpp
//

#include <iostream>
#include <iomanip>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>

extern "C"
{
#include "mri.h"

  char *Progname="difftool";
}

using namespace std;

int main(int argc, char *argv[])
{
  string f1;
  string f2;
  string outfile;

  if (argc < 3)
  {
    // tie in and out
    cin.tie(&cout);
    cout << "1st filename :";

    cin >> f1;
    cout << "2nd filename :";

    cin >> f2;

    cout << "output filename :";

    cin >> outfile;
  }
  else
  {
    f1 = argv[1];
    f2 = argv[2];
    if (argv[3])
      outfile = argv[3];
    else
      outfile = "";
  }

  MRI *mri1 = MRIread(const_cast<char *> (f1.c_str()));
  if (mri1==0)
  {
    cout << "could not read " << f1.c_str() << endl;
    return -1;
  }

  MRI *mri2 = MRIread(const_cast<char *> (f2.c_str()));
  if (mri2==0)
  {
    cout << "could not read " << f2.c_str() << endl;
    return -1;
  }
  // verify width, height, depth and type are the same
  if ( (mri1->width != mri2->width)
       || (mri1->height != mri2->height)
       || (mri1->depth != mri2->depth)
       || (mri1->type != mri2->type))
  {
    cout << "size or type did not match" << endl;
    return -1;
  }
  // now diff
  MRI *res = MRIalloc(mri1->width, mri1->height, mri1->depth, mri1->type);
  MRIcopyHeader(mri1, res);

  unsigned char u1, u2, udiff;

  for (int k=0; k < mri1->depth; ++k)
    for (int j=0; j < mri1->height; ++j)
      for (int i=0; i < mri1->width; ++i)
      {
        switch (mri1->type)
        {
        case MRI_UCHAR:
          u1 = MRIvox(mri1, i,j,k);
          u2 = MRIvox(mri2, i,j,k);
          udiff = (u1-u2)/2 + 127;
          MRIvox(res, i,j,k) = udiff;
          if (u1 != u2)
            cout << "(" << i << ", " << j << ", " << k << ") = (" << (int) u1 << ", " << (int) u2 << ")" << endl;
          break;
        case MRI_SHORT:
        case MRI_FLOAT:
        default:
          cout << "Sorry.  currently supports only uchar, short, or float" << endl;
          goto cleanup;
        }
      }

  if (outfile.length())
    MRIwrite(res, const_cast<char *>(outfile.c_str()));

cleanup:
  cout << "Cleanup" << endl;
  MRIfree(&mri1);
  MRIfree(&mri2);
  MRIfree(&res);
  cout << "Done" << endl;
}
