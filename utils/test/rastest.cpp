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

#include "mri.h"
const char *Progname = "rastest";

using namespace std;

double t = 0.000001;

bool equal(double &x, double &y, double &tol) {
  if (fabs(x - y) < tol)
    return true;
  else
    return false;
}
int main(int argc, char *argv[]) {
  // verify various vox, ras, surfaceRAS, relations
  if (argc < 2) {
    std::cout << "Usage: rastest <vol>" << std::endl;
    exit(1);
  }
  MRI *  mri = MRIread(argv[1]);
  double xr, yr, zr, xr1, yr1, zr1;
  double xsr, ysr, zsr, xsr1, ysr1, zsr1;
  double xv, yv, zv, xv1, yv1, zv1, xv2, yv2, zv2;
  char   more = 'y';
  while (more != 'n') {
    std::cout << "vox x: ";
    cin >> xv;
    std::cout << "vox y: ";
    cin >> yv;
    std::cout << "vox z: ";
    cin >> zv;
    std::cout << "Going from Voxel to SurfaceRAS " << std::endl;
    // one way to get to surface RAS
    MRIvoxelToWorld(mri, xv, yv, zv, &xr, &yr, &zr);
    MRIRASToSurfaceRAS(mri, xr, yr, zr, &xsr, &ysr, &zsr);
    // another way to get to surface RAS
    MRIvoxelToSurfaceRAS(mri, xv, yv, zv, &xsr1, &ysr1, &zsr1);
    // compare
    if (!equal(xsr, xsr1, t) || !equal(ysr, ysr1, t) || !equal(zsr, zsr1, t))
      std::cout << "Error" << std::endl;
    else
      std::cout << "Fine" << std::endl;
    std::cout << '\t' << xsr << " <-> " << xsr1 << std::endl;
    std::cout << '\t' << ysr << " <-> " << ysr1 << std::endl;
    std::cout << '\t' << zsr << " <-> " << zsr1 << std::endl;
    //////////////////////////////////////////////////////////////
    std::cout << "Coming back from SurfaceRAS to Voxel " << std::endl;
    // one way
    MRIsurfaceRASToRAS(mri, xsr, ysr, zsr, &xr1, &yr1, &zr1);
    MRIworldToVoxel(mri, xr1, yr1, zr1, &xv1, &yv1, &zv1);
    // another way
    MRIsurfaceRASToVoxel(mri, xsr, ysr, zsr, &xv2, &yv2, &zv2);
    // compare
    if (!equal(xv1, xv2, t) || !equal(yv1, yv2, t) || !equal(zv1, zv2, t))
      std::cout << "Error *****************************" << std::endl;
    else
      std::cout << "Fine" << std::endl;
    std::cout << '\t' << xv1 << " <-> " << xv2 << std::endl;
    std::cout << '\t' << yv1 << " <-> " << yv2 << std::endl;
    std::cout << '\t' << zv1 << " <-> " << zv2 << std::endl;

    std::cout << "RAS location comparision going and coming back" << std::endl;
    // RAS comparison
    if (!equal(xr, xr1, t) || !equal(yr, yr1, t) || !equal(zr, zr1, t))
      std::cout << "Error ****************************" << std::endl;
    else
      std::cout << "Fine" << std::endl;
    std::cout << '\t' << xr << " <-> " << xr1 << std::endl;
    std::cout << '\t' << yr << " <-> " << yr1 << std::endl;
    std::cout << '\t' << zr << " <-> " << zr1 << std::endl;

    std::cout << "More? (y/n) ";
    cin >> more;
  }
}
