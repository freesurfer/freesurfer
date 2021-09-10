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

#include "mris_sphshapepvf.h"

int SphShapePVFSim::MakeMRI(int dim, double res){
  // just a convenient way to make an MRI, 
  // might just want to load one from file
  vol = MRIalloc(dim,dim,dim,MRI_FLOAT);
  if(vol==NULL) return(1);
  vol->xsize = res;
  vol->ysize = res;
  vol->zsize = res;
  MATRIX *M = MRIxfmCRS2XYZtkreg(vol);
  // Set the vox2ras to the tkreg vox2ras
  MRIsetVox2RASFromMatrix(vol, M);
  MatrixFree(&M);
  if(surf) getVolGeom(vol, &(surf->vg));
  return(0);
}

int SphShapePVFSim::LoadIcoSurf(void){
  // Load the icosahedron. Radius=100 is meaningless
  // Note: you can load your own surface instead
  surf = ReadIcoByOrder(icoOrder,100);
  if(surf==NULL) return(1);
  if(vol) getVolGeom(vol, &(surf->vg));
  MRIScomputeMetricProperties(surf);
  return(0);
}

double SphShapePVFSim::RadiusFunc(const std::array<double,3>xyz){
  // Compute theta and phi from xyz, then run the virtual RadiusFunc()
  std::array<double,3> rtp = XYZ2RTP(xyz);
  double radius = RadiusFunc(rtp[1],rtp[2]);
  return(radius);
}
std::array<double,3> SphShapePVFSim::NormalFunc(const std::array<double,3>xyz){
  // Compute theta and phi from xyz, then run the virtual NormalFunc()
  std::array<double,3> rtp = XYZ2RTP(xyz);
  std::array<double,3> nxyz = NormalFunc(rtp[1],rtp[2]);
  return(nxyz);
}

std::array<double,3> SphShapePVFSim::XYZFunc(const double thetaRad, const double phiRad){
  // Return the xyz of the point on the shape at the coordinate
  // theta,phi Radius already determined form abstract function. This
  // can be used to set the vertex coords of the surface
  double radius = RadiusFunc(thetaRad,phiRad);
  std::array<double,3> rtp = {radius,thetaRad,phiRad};
  std::array<double,3> xyz = RTP2XYZ(rtp);
  return(xyz);
};

std::array<double,3> SphShapePVFSim::XYZ2RTP(std::array<double,3> xyz){
  // This just computes spherical coords from xyz
  std::array<double,3> rtp;
  rtp[0] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
  rtp[1] = acos(xyz[2]/rtp[0]);
  rtp[2] = atan2(xyz[1],xyz[0]);
  return(rtp);
}
std::array<double,3> SphShapePVFSim::RTP2XYZ(std::array<double,3> rtp){
  // This just computes xyz from spherical coords
  std::array<double,3> xyz;
  xyz[0] = rtp[0]*cos(rtp[2])*sin(rtp[1]);
  xyz[1] = rtp[0]*sin(rtp[2])*sin(rtp[1]);
  xyz[2] = rtp[0]*cos(rtp[1]);
  return(xyz);
}

int SphShapePVFSim::SetSurfXYZ(void){
  if(surf==NULL) {
    printf("ERROR: SetSurfXYZ(): must load a surface first\n");
    return(1);
  }
  int n;
  VERTEX *v;
  for(n=0; n < surf->nvertices; n++){
    v = &(surf->vertices[n]);
    std::array<double,3> xyz = {v->x,v->y,v->z};
    std::array<double,3> rtp = SphShapePVFSim::XYZ2RTP(xyz);
    //std::array<double,3> rtp = {1,2,3};
    xyz = XYZFunc(rtp[1],rtp[2]);
    v->x = xyz[0];
    v->y = xyz[1];
    v->z = xyz[2];
  }
  if(vol) getVolGeom(vol, &(surf->vg));
  MRIScomputeMetricProperties(surf);
  return(0);
}

/*!
  double SphShapePVFSim::NormDot(MRIS *testsurf)
  \brief Compute the dot product between the normal vector at vertex
  and the theoretical normal (should be 1.0).
 */

double SphShapePVFSim::NormDot(MRIS *testsurf)
{
  double doterr=0;
  int n;
  VERTEX *v;
  for(n=0; n < testsurf->nvertices; n++){
    v = &(surf->vertices[n]);
    std::array<double,3> xyz = {v->x,v->y,v->z};
    std::array<double,3> nxyz = NormalFunc(xyz);
    double dot = v->nx*nxyz[0] + v->ny*nxyz[1] + v->nz*nxyz[2];
    doterr += (1-dot)*(1-dot);
    v->val = dot;
    if(debug && n < 10){
      printf("%d   %6.4lf %6.4lf %6.4lf   %6.4lf %6.4lf %6.4lf   %g\n",
	     n,v->nx,v->ny,v->nz,nxyz[0],nxyz[1],nxyz[2],dot);
    }
  }
  doterr = sqrt(doterr);
  printf("dotsum = %g\n",doterr);
  return(doterr);
}

/*!
  int SphShapePVFSim::ComputePVF(void)
  Compute partial volume fraction at each voxel given the shape
*/
int SphShapePVFSim::ComputePVF(void)
{
  double voxR,VolTrue;
  double res=fsubsample;
  MATRIX *tkvox2ras = MRIxfmCRS2XYZtkreg(vol);
  int c;

  // max dist across a voxel
  voxR = sqrt(pow(vol->xsize,2.0)+pow(vol->ysize,2.0)+pow(vol->zsize,2.0));

  double SumPVF = 0.0;
#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(+ : SumPVF)
#endif
  for (c = 0; c < vol->width; c++) {
    if(debug>0) printf("c=%d\n",c);fflush(stdout);
    int r,s;
    double D,vpvf;
    MATRIX *crs, *ras=NULL;
    crs = MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
    for (r = 0; r < vol->height; r++) {
      for (s = 0; s < vol->depth; s++) {
	crs->rptr[1][1] = c;
	crs->rptr[2][1] = r;
	crs->rptr[3][1] = s;
	// Compute RAS in tkreg space
	ras = MatrixMultiplyD(tkvox2ras,crs,ras);
	// Get radius of the shape on the ray pointing through this voxel
	std::array<double,3> xyz = {ras->rptr[1][1],ras->rptr[2][1],ras->rptr[3][1]};
	double radius = RadiusFunc(xyz);
	// Distance to voxel
	D = sqrt(pow(ras->rptr[1][1],2.0)+pow(ras->rptr[2][1],2.0)+pow(ras->rptr[3][1],2.0));
	if(D < radius-voxR){
	  // fully inside, pvf=1
	  MRIsetVoxVal(vol,c,r,s,0,1.0);
	  SumPVF += 1.0;
	  continue;
	}
	if(D > radius+voxR){
	  // fully outside, pvf=0
	  MRIsetVoxVal(vol,c,r,s,0,0.0);
	  continue;
	}
	int count=0, ntot=0;
	double dc,dr,ds;
	for(dc = -0.5; dc < 0.5; dc += res){
	  for(dr = -0.5; dr < 0.5; dr += res){
	    for(ds = -0.5; ds < 0.5; ds += res){
	      ntot++;
	      crs->rptr[1][1] = c+dc;
	      crs->rptr[2][1] = r+dr;
	      crs->rptr[3][1] = s+ds;
	      ras = MatrixMultiplyD(tkvox2ras,crs,ras);
	      D = sqrt(pow(ras->rptr[1][1],2.0)+pow(ras->rptr[2][1],2.0)+pow(ras->rptr[3][1],2.0));
	      if(D < radius) count++; // subvoxel is fully inside radius
	    }
	  }
	}
	vpvf = (double)count/(double)ntot;
	MRIsetVoxVal(vol,c,r,s,0,vpvf);
	SumPVF += vpvf;
      }
    }
    MatrixFree(&crs);
    MatrixFree(&ras);
  }
  VolPVF = SumPVF;

  // This is just a check
  VolTrue = ShapeVolume();
  printf("VolTrue = %g, VolPVF = %g, err=%g %%\n",VolTrue,VolPVF,100*(VolTrue-VolPVF)/VolTrue);

  MatrixFree(&tkvox2ras);

  return(0);
}


