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

#include <iomanip>
#include <iostream>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>

#include "error.h"
#include "gcamorph.h"
#include "mri.h"
#include "transform.h"
#include "utils.h"
const char *Progname = "testm3d";

using namespace std;

int main(int argc, char *argv[]) {
  std::cout << "m3d file: ";
  std::string infile;
  cin >> infile;

  std::cout << "mri file: ";
  std::string mrifile;
  cin >> mrifile;

  printf("Try 1 *************************************************\n");
  printf("before loading the transform : heap usage %d\n", getMemoryUsed());
  TRANSFORM *transform = 0;
  std::cout << "reading transform file " << infile.c_str() << std::endl;
  transform = TransformRead(const_cast<char *>(infile.c_str()));
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
              Progname, const_cast<char *>(infile.c_str()));
  printf("after loading the transform : heap usage %d\n", getMemoryUsed());
  TransformFree(&transform);
  printf("after freeing transform : heap usage %d\n", getMemoryUsed());

  printf("Try 2 *************************************************\n");
  printf("before loading the transform : heap usage %d\n", getMemoryUsed());
  transform = 0;
  std::cout << "reading transform file " << infile.c_str() << std::endl;
  transform = TransformRead(const_cast<char *>(infile.c_str()));
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
              Progname, const_cast<char *>(infile.c_str()));

  printf("before loading mri : heap usage %d\n", getMemoryUsed());
  std::cout << "reading mri file " << mrifile.c_str() << std::endl;
  MRI *mri = MRIread(const_cast<char *>(mrifile.c_str()));

  printf("after  loading mri : heap usage %d\n", getMemoryUsed());

  // modify transform to store inverse also
  std::cout << "TransformInvert processing ..." << std::endl;
  TransformInvert(transform, mri);

  printf("after inverting the transform  heap usage %d\n", getMemoryUsed());

  std::cout << "Free memory..." << std::endl;
  MRIfree(&mri);
  printf("after freeing mri : heap usage %d\n", getMemoryUsed());
  TransformFree(&transform);
  printf("after freeing transform : heap usage %d\n", getMemoryUsed());

  std::cout << "Done" << std::endl;
}
