/**
 * @brief 
 *
 */
/*
 * Original Author: Doug Greve
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef MRIS_SPHSHAPEPVF
#define MRIS_SPHSHAPEPVF

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <float.h>
#include "mrisurf.h"
#include "mrisutils.h"
#include "utils.h"
#include "mri2.h"
#include "icosahedron.h"
#include "error.h"
#include "dmatrix.h"
#include "matrix.h"
#include "diag.h"
#include <vector>
#include <array>

#undef private

class SphShapePVFSim {
  /* This is an abstract class used to manage simulations of partial
     volume (PV) effects in shapes defined by a spherical
     function. This function takes theta,phi spherical coordinates and
     returns the radius. The most simple example is that of a sphere
     where the function always returns a constant (the radius), but
     more complicated functions can exist. All coordinates are assumed
     to be in tkreg space. The partial volume fraction (PVF) in a
     given voxel is computed by dividing up the voxel into subvoxels
     (fsubsample) and determining whether the subvoxel is inside or
     outside the shape.
   */
public:
  MRI  *vol=NULL;
  MRIS *surf=NULL;
  int icoOrder=5;
  double fsubsample=0.1;
  double VolPVF=0; // volume of shape as computed from the PVF image
  int debug = 0;
  std::string ShapeName;

  // RadiusFunc(theta,phi): radius of the point on the shape at the coordinate theta,phi
  virtual double RadiusFunc(double theta, double phi)=0;
  // NormalFunc(theta,phi): normal vector at the point on the shape at the coordinate theta,phi
  virtual std::array<double,3> NormalFunc(double theta, double phi)=0;
  // ShapeVolume() Ideal volume inside the shape (good for error checking)
  virtual double ShapeVolume(void)=0; 

  // Fills the voxels in vol with the partial volume fraction (PVF)
  int ComputePVF(void);

  // RadiusFunc(xyz): radius of the point on the shape that the ray through xyz passes
  double RadiusFunc(const std::array<double,3>xyz);
  // NormalFunc(xyz): normal vector at the point on the shape at the
  // that the ray through xyz passes
  std::array<double,3> NormalFunc(const std::array<double,3>xyz);
  // XYZFunc(): return the xyz of the point on the shape at the coordinate theta,phi
  // This is used to get the vertex coords of the shape (eg, to make a surf)
  std::array<double,3> XYZFunc(const double thetaRad, const double phiRad);
  std::array<double,3> XYZ2RTP(std::array<double,3> xyz);
  std::array<double,3> RTP2XYZ(std::array<double,3> rtp);
  double NormDot(MRIS *surf2);
  int SetSurfXYZ(void);
  int MakeMRI(int dim, double res);
  int LoadIcoSurf(void);
};

//===================================================================
class BasicSpherePVF : public SphShapePVFSim {
public:
  double radius = 25; // could add center
  BasicSpherePVF(){
    ShapeName = "BasicSphere";
  }
  double RadiusFunc(double theta, double phi){
    return(radius);
  }
  std::array<double,3> NormalFunc(double theta, double phi){
    std::array<double,3> rtp = {radius,theta,phi};
    std::array<double,3> nxyz = RTP2XYZ(rtp);
    nxyz[0] /= rtp[0];
    nxyz[1] /= rtp[0];
    nxyz[2] /= rtp[0];
    return(nxyz);
  };
  double ShapeVolume(void){
    return(4.0*M_PI*pow(radius,3.0)/3);
  }
};

#endif
