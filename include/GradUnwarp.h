#ifndef GRADUNWARP_H
#define GRADUNWARP_H

#include "gcamorph.h"
#include "mri.h"
//#include "vol_geom.h"

typedef struct
{
  int num;
  char A_or_B;
  int n;
  int m;
  float value;
  char xyz;
} COEFF;

class Siemens_B
{
public:
  Siemens_B(int coeffDim0, int nmax0, float R0, double **normfact0, float X, float Y, float Z);
  ~Siemens_B();

  float siemens_B_x(float **Alpha_x, float **Beta_x);
  float siemens_B_y(float **Alpha_y, float **Beta_y);
  float siemens_B_z(float **Alpha_z, float **Beta_z);
  float siemens_B(float **Alpha, float **Beta);
  void  siemens_legendre(int n, double x);

private:
  int coeffDim;
  int nmax;
  float R0_mm;

  double **normfact;

  double **P;
  float  *F;
  double *cosPhi, *sinPhi;
  float R, Theta, Phi;
};

class GradUnwarp
{
public:
  GradUnwarp(int nthreads = 1);
  ~GradUnwarp();

  void  read_siemens_coeff(const char *gradfilename);
  void  printCoeff();

  void  initSiemensLegendreNormfact();
  void  spharm_evaluate(float X, float Y, float Z, float *Dx, float *Dy, float *Dz);

  void  create_transtable(VOL_GEOM *vg, MATRIX *vox2ras, MATRIX *inv_vox2ras);
  void  load_transtable(const char* morphfile);
  void  save_transtable(const char* morphfile);

  MRI*  unwarp_volume_gradfile(MRI *warpedvol, MRI *unwarpedvol, MATRIX *vox2ras, MATRIX *inv_vox2ras, int interpcode, int sinchw);
  MRI*  unwarp_volume(MRI *warpedvol, MRI *unwarpedvol, int interpcode, int sinchw);
  MRIS* unwarp_surface_gradfile(MRIS *warpedsurf, MRIS *unwarpedsurf);
  MRIS* unwarp_surface(MRIS *warpedsurf, MRIS *unwarpedsurf);
  
private:
  int nthreads;

  FILE *fgrad;

  COEFF *coeff;

  int nmax;
  int mmax;

  int coeffCount;
  int coeffDim;

  bool Alpha_Beta_initialized;

  float R0;
  float **Alpha_x, **Alpha_y, **Alpha_z;
  float **Beta_x,  **Beta_y,  **Beta_z;

  // these variables are used to calculate siemens legendre
  // they are pre-calculated in initSiemensLegendreNormfact()
  double *minusonepow;
  double *factorials;
  double **normfact;

  GCAM *gcam;

private:
  void _skipCoeffComment();
  void _initCoeff();
  void _update_GCAMnode(int c, int r, int s, float fcs, float frs, float fss);
  void _assignUnWarpedVolumeValues(MRI* warpedvol, MRI* unwarpedvol, MRI_BSPLINE *bspline, int interpcode, int sinchw,
                                   int c, int r, int s, float fcs, float frs, float fss);
  void _printMatrix(MATRIX *matrix, const char *desc);
};

#endif
