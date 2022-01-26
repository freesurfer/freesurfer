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
  GradUnwarp();
  ~GradUnwarp();

  void  setup() {};
  void  list_coeff_files() {};

  void  read_siemens_coeff(const char *gradfilename);
  void  printCoeff();

  void initSiemensLegendreNormfact();
  void  spharm_evaluate(float X, float Y, float Z, float *Dx, float *Dy, float *Dz);

  void  unwarp() {};
  void  unwap_volume() {};
  void  unwarp_surface() {};
  
private:
  FILE *fgrad;

  COEFF coeff[100];

  int nmax;
  int mmax;

  int coeffCount;
  int coeffDim;

  float R0;
  float **Alpha_x, **Alpha_y, **Alpha_z;
  float **Beta_x,  **Beta_y,  **Beta_z;

  // these variables are used to calculate siemens legendre
  // they are pre-calculated in initSiemensLegendreNormfact()
  double *minusonepow;
  double *factorials;
  double **normfact;

private:
  void _skipCoeffComment();
  void _initCoeff();
};

#endif
