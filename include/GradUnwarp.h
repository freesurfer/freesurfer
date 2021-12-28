#ifndef GRADUNWARP_H
#define GRADUNWARP_H

typedef struct
{
  int num;
  char A_or_B;
  int n;
  int m;
  float value;
  char xyz;
} COEFF;

class GradUnwarp
{
public:
  GradUnwarp();
  ~GradUnwarp();

  void  setup();
  void  list_coeff_files();

  void  read_siemens_coeff(const char *gradfilename);
  void  printCoeff();

  float siemens_B_x();
  float siemens_B_y();
  float siemens_B_z();
  float siemens_B(float **Alpha, float **Beta, float **F2);
  void  siemens_B0(float X, float Y, float Z);
  void  siemens_legendre(int n, double x, double **P);

  void  spharm_evaluate(float X, float Y, float Z, float *Dx, float *Dy, float *Dz);

  void  initSiemensLegendreNormfact();

  void  unwarp();
  void  unwap_volume();
  void  unwarp_surface();
  
private:
  FILE *fgrad;

  COEFF coeff[100];

  int coeffCount;
  int coeffDim;

  float R0;
  float **Alpha_x, **Alpha_y, **Alpha_z;
  float **Beta_x,  **Beta_y,  **Beta_z;

  // these variables are used to calculate siemens legendre
  double *minusonepow;
  double *factorials;
  double **normfact;
  
  double **P;
  float  *F;
  float **F2_x, **F2_y, **F2_z;
  double *cosPhi, *sinPhi;
  float R, Theta, Phi;

private:
  void _skipCoeffComment();
  void _initCoeff();
};

#endif
