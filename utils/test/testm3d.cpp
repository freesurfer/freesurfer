/**
 * @file  testm3d.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:46 $
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


//
// testm3d.cpp
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
#include "gcamorph.h"
#include "transform.h"
#include "error.h"
#include "utils.h"
  char *Progname="testm3d";
}

using namespace std;

int main(int argc, char *argv[])
{
  cout << "m3d file: ";
  string infile;
  cin >> infile;

  cout << "mri file: ";
  string mrifile;
  cin >> mrifile;

  printf("Try 1 *************************************************\n");
  printf("before loading the transform : heap usage %d\n", getMemoryUsed());
  TRANSFORM *transform=0;
  cout << "reading transform file " << infile.c_str() << endl;
  transform = TransformRead(const_cast<char *> (infile.c_str()));
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
              Progname, const_cast<char *> (infile.c_str()));
  printf("after loading the transform : heap usage %d\n", getMemoryUsed());
  TransformFree(&transform);
  printf("after freeing transform : heap usage %d\n", getMemoryUsed());

  printf("Try 2 *************************************************\n");
  printf("before loading the transform : heap usage %d\n", getMemoryUsed());
  transform=0;
  cout << "reading transform file " << infile.c_str() << endl;
  transform = TransformRead(const_cast<char *> (infile.c_str()));
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
              Progname, const_cast<char *> (infile.c_str()));

  printf("before loading mri : heap usage %d\n", getMemoryUsed());
  cout << "reading mri file " << mrifile.c_str() << endl;
  MRI *mri = MRIread(const_cast<char *> (mrifile.c_str()));

  printf("after  loading mri : heap usage %d\n", getMemoryUsed());

  // modify transform to store inverse also
  cout << "TransformInvert processing ..." << endl;
  TransformInvert(transform, mri);

  printf("after inverting the transform  heap usage %d\n", getMemoryUsed());

  cout << "Free memory..." << endl;
  MRIfree(&mri);
  printf("after freeing mri : heap usage %d\n", getMemoryUsed());
  TransformFree(&transform);
  printf("after freeing transform : heap usage %d\n", getMemoryUsed());

  cout << "Done" << endl;
}
