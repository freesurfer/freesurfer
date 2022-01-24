#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "GradUnwarp.h"
#include "legendre.h"

/*****************************************************************************/
/******************** Implementation of GradUnwarp class *********************/
/**********    3 environment variables to enable debug info:        **********/
/********** COEFF_READ_DEBUG, LEGENDRE_NORMFACT_DEBUG, PRN_LEGENDRE **********/
GradUnwarp::GradUnwarp()
{
  fgrad = NULL;

  nmax = 0;
  mmax = 0;

  coeffCount = 0;
  coeffDim = 0;

  R0 = 0;
  Alpha_x = NULL; Alpha_y = NULL; Alpha_z = NULL;
  Beta_x  = NULL; Beta_y  = NULL; Beta_z  = NULL;
}

GradUnwarp::~GradUnwarp()
{
  int i;
  for (i = 0; i < coeffDim; i++)
  {
    free(Alpha_x[i]);
    free(Alpha_y[i]);
    free(Alpha_z[i]);
    free(Beta_x[i]);
    free(Beta_y[i]);
    free(Beta_z[i]);
  }

  free(Alpha_x);
  free(Alpha_y);
  free(Alpha_z);
  free(Beta_x);
  free(Beta_y);
  free(Beta_z);  

  free(minusonepow);
  free(factorials);

  for (i = 0; i < coeffDim; i++)
    free(normfact[i]);
  free(normfact);
}

void GradUnwarp::read_siemens_coeff(const char *gradfilename)
{
  // check if gradfile has extension .grad
  // ...

  // open gradfilename
  fgrad = fopen(gradfilename, "r");
  if (fgrad == NULL)
  {
    printf("ERROR: could not read the coefficient file %s\n", gradfilename);
    return;
  }

  printf("==> reading coefficients from gradient coil file %s\n" , gradfilename) ;

  _skipCoeffComment();

  // hard-coded limits:
  // number of coeff entries - 100
  // length of each entry    - 1024
  char coeffline[1024];

  // skip the next line. (It contains an information about the system type.)
  fgets(coeffline, sizeof(coeffline), fgrad);

  printf("==> reading system type string from coeff.grad file...\n");
  printf("%s\n", coeffline);

  // check if first paramline contains "win_"; and, if so, parse
  fgets(coeffline, sizeof(coeffline), fgrad);
  if (strncmp(coeffline, " win_", 5) == 0)
  {

    // parse into the four parameters (these don't seem to be used anywhere...)
    //[dum, iThreshLo, dum, iThreshUp, dum, iAlgoTyp, dum, iWinDummy] ...
    //	= strread(paramline, " %10c%d%13c%d%13c%d%14c%d;");

    // read next line
    fgets(coeffline, sizeof(coeffline), fgrad);
  }

  // only extract radius and ignore rest
  sscanf(coeffline, "%f\n", &R0);
  printf("WARN: returning R0 = %f in units of METERS!\n", R0);
  R0 = R0 * 1000;  // R0 is now in mm
  
  // read next line, which contains gradient system mode "(0 = typ. tunnel magnet system; 1 = typ. open magnet system)"
  fgets(coeffline, sizeof(coeffline), fgrad);
  int CoSyMode;
  sscanf(coeffline, "%d = %*s\n", &CoSyMode);
  printf("CoSyMode = %d\n", CoSyMode);

  // skip the next 5 lines
  for (int ind = 0; ind < 5; ind++)
      fgets(coeffline, sizeof(coeffline), fgrad);


  /***********************************************************/
  /****** begin reading spherical harmonic coefficients ******/
  while (fgets(coeffline, sizeof(coeffline), fgrad) != NULL)
  {
    int len = strlen(coeffline);
    char* ptr = coeffline;
    char* endptr = ptr + len;

    // skip leading spaces and tabs, also empty lines
    while (ptr != endptr && (*ptr == ' ' || *ptr == '\t' || *ptr == '\n'))
      ptr++;

    if (*ptr == '\0')
      continue;

    // start filling coeff entries at index 1
    coeffCount++;

    if (getenv("COEFF_READ_DEBUG"))
      printf("entry #%d: %s\n", coeffCount, coeffline);

    sscanf(ptr, "%d %c(%d, %d) %f %c", 
                      &coeff[coeffCount].num, &coeff[coeffCount].A_or_B, &coeff[coeffCount].n, &coeff[coeffCount].m, 
                      &coeff[coeffCount].value, &coeff[coeffCount].xyz);
    nmax = (coeff[coeffCount].n > nmax) ? coeff[coeffCount].n : nmax;
    mmax = (coeff[coeffCount].m > mmax) ? coeff[coeffCount].m : mmax;
    //printf("%d %c (%d, %d) %f %c\n", coeff[coeffCount].num, 
    //       coeff[coeffCount].A_or_B, coeff[coeffCount].n, coeff[coeffCount].m, coeff[coeffCount].value, coeff[coeffCount].xyz);
  }

  fclose(fgrad);    

  coeffDim = (nmax > mmax) ? nmax+1 : mmax+1;
  _initCoeff();

  /**************************************************************************/
  /****** organize coefficient values ******/
  int n;
  for (n = 1; n <= coeffCount; n++)
  {
    float **arrPtr = NULL;
    if (coeff[n].A_or_B == 'A')
    {
      if (coeff[n].xyz == 'x')
        arrPtr = Alpha_x;
      else if (coeff[n].xyz == 'y')
        arrPtr = Alpha_y;
      else if (coeff[n].xyz == 'z')
        arrPtr = Alpha_z;
    }   
    else if (coeff[n].A_or_B == 'B')
    {
      if (coeff[n].xyz == 'x')
        arrPtr = Beta_x;
      else if (coeff[n].xyz == 'y')
        arrPtr = Beta_y;
      else if (coeff[n].xyz == 'z')
        arrPtr = Beta_z;
    }
    
    if (arrPtr == NULL)
    {
      printf("ERROR: unrecognized coefficient string: '%c%c'\n", coeff[n].A_or_B, coeff[n].xyz);
      continue;
    }

    int row = coeff[n].n;
    int col = coeff[n].m; 
    arrPtr[row+1][col+1] = coeff[n].value;
  }
}

void GradUnwarp::initSiemensLegendreNormfact()
{
  // initialize variables to pre-calculate normfact for siemens_legendre()
  minusonepow = new double[coeffDim];
  factorials  = new double[2*coeffDim];
  normfact = new double*[coeffDim];

  int n;
  for (n = 0; n < coeffDim; n++)
    normfact[n] = new double[coeffDim];

  // pre-calculate minusonepow, factorials, & normfact
  for (n = 0; n < coeffDim; n++)
    minusonepow[n] = pow((-1), n);

  for (n = 0; n < 2*coeffDim; n++)
    factorials[n] = factorial(n);

  if (getenv("LEGENDRE_NORMFACT_DEBUG"))
    printf("\n");

  for (n = 0; n < coeffDim; n++)
  {
    int m;
    for (m = 1; m <= n; m++)
    {
      normfact[n][m] = minusonepow[m] * sqrt((2*n+1)*factorials[n-m]/(2*factorials[n+m]));
   
      if (getenv("LEGENDRE_NORMFACT_DEBUG"))
        printf("normfact2[%2d][%2d] = %s%.6lf, pow((-1), %2d) * sqrt((%2d)*factorial(%2d)/(2*factorial(%2d)))\n",
               n, m, (normfact[n][m] > 0) ? " " : "", normfact[n][m], m, 2*n+1, n-m, n+m);
    }

    if (getenv("LEGENDRE_NORMFACT_DEBUG"))
      printf("\n");
  }
}

void GradUnwarp::spharm_evaluate(float X, float Y, float Z, float *Dx, float *Dy, float *Dz)
{
  Siemens_B *siemens_B = new Siemens_B(coeffDim, nmax, R0, normfact, X, Y, Z);

  float bx = siemens_B->siemens_B_x(Alpha_x, Beta_x);
  float by = siemens_B->siemens_B_y(Alpha_y, Beta_y);
  float bz = siemens_B->siemens_B_z(Alpha_z, Beta_z);

  //printf("bx=%lf, by=%lf, bz=%lf\n", bx, by, bz);

  *Dx = bx * R0;
  *Dy = by * R0;
  *Dz = bz * R0;

  delete siemens_B;
}

void GradUnwarp::printCoeff()
{
  const char *arrs = "AxAyAzBxByBz";

  float **arrPtr = NULL;

  int i = 0;
  for (; i < strlen(arrs); i++)
  {
    if (arrs[i] == 'A')
    {
      i++;
      if (arrs[i] == 'x')
        arrPtr = Alpha_x;
      else if (arrs[i] == 'y')
        arrPtr = Alpha_y;
      else if (arrs[i] == 'z')
        arrPtr = Alpha_z;
    }   
    else if (arrs[i] == 'B')
    {
      i++;
      if (arrs[i] == 'x')
        arrPtr = Beta_x;
      else if (arrs[i] == 'y')
        arrPtr = Beta_y;
      else if (arrs[i] == 'z')
        arrPtr = Beta_z;
    }
    
    if (arrPtr == NULL)
    {
      printf("ERROR: unrecognized coefficient string: '%c%c'\n", arrs[i-1], arrs[i]);
      continue;
    }

    printf("\n%s_%c = \n", (arrs[i-1] == 'A') ? "Alpha" : "Beta", arrs[i]);
    int row = 0;
    for (; row < coeffDim; row++)
    {
      int col = 0;
      for (; col < coeffDim; col++)
      {
        printf("\t%lf", arrPtr[row][col]); 
      }
      printf("\n");
    }
  }
}

void GradUnwarp::_skipCoeffComment()
{
  char line[1024];
  while (fgets(line, sizeof(line), fgrad) != NULL)
  {
    int len = strlen(line);
    char* ptr = line;
    char* endptr = ptr + len;

    // skip leading spaces and tabs, also empty lines
    while (ptr != endptr && (*ptr == ' ' || *ptr == '\t' || *ptr == '\n'))
      ptr++;

    // skip the comment lines. The comment section ends with line #*] END: 
    if (*ptr != '\0' && strncmp(ptr, "#*] END:", 8) == 0)
      return;
  }
}

void GradUnwarp::_initCoeff()
{
  printf("coeffDim = %d, coeffCount=%d\n", coeffDim, coeffCount);

  // allocate float *Alpha_x, *Alpha_y, *Alpha_z; float *Beta_x,  *Beta_y,  *Beta_z;
  // first row and first col of Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z are all zeros, they are not used.
  // this is to align with MATLAB matrix which has index starting at 1.
  int arrsize = coeffDim + 1;
  Alpha_x = new float*[arrsize];
  Alpha_y = new float*[arrsize];
  Alpha_z = new float*[arrsize];
  Beta_x  = new float*[arrsize];
  Beta_y  = new float*[arrsize];
  Beta_z  = new float*[arrsize];

  // initialize coefficient arrays
  int i;
  for (i = 0; i <= coeffDim; i++)
  {
    Alpha_x[i] = new float[arrsize];
    memset(Alpha_x[i], 0, sizeof(float)*arrsize);

    Alpha_y[i] = new float[arrsize];
    memset(Alpha_y[i], 0, sizeof(float)*arrsize);

    Alpha_z[i] = new float[arrsize];
    memset(Alpha_z[i], 0, sizeof(float)*arrsize);

    Beta_x[i]  = new float[arrsize];
    memset(Beta_x[i], 0, sizeof(float)*arrsize);

    Beta_y[i]  = new float[arrsize];
    memset(Beta_y[i], 0, sizeof(float)*arrsize);

    Beta_z[i]  = new float[arrsize];
    memset(Beta_z[i], 0, sizeof(float)*arrsize);
  }
}


/***********************************************************************************/
/************************ Implementation of Siemens_B class ************************/
Siemens_B::Siemens_B(int coeffDim0, int nmax0, float R0, double **normfact0, float X, float Y, float Z)
{
  coeffDim = coeffDim0;    // coeffDime = nmax + 1
  nmax = nmax0;
  R0_mm = R0;
  normfact = normfact0;

  P = new double*[coeffDim];

  int n;
  for (n = 0; n < coeffDim; n++)
  {
    // P[0][] is not used, MATLAB index starts at 1 
    P[n] = new double[coeffDim+1];
    memset(P[n], 0, (coeffDim+1)*sizeof(double));
  }

  F = new float[coeffDim];
  cosPhi = new double[coeffDim];
  sinPhi = new double[coeffDim];

  memset(F, 0, coeffDim*sizeof(float));
  memset(cosPhi, 0, coeffDim*sizeof(double));
  memset(sinPhi, 0, coeffDim*sizeof(double));

  // hack to avoid singularities at origin (R==0)
  X = X+0.0001;

  // convert to spherical coordinates
  R = sqrt(X*X + Y*Y + Z*Z);
  Theta = acos(Z/R);
  Phi = atan2(Y/R, X/R);

  // evaluate the Legendre polynomial (using Siemens's normalization)
  for (n = 0; n <= nmax; n++)
  {
    siemens_legendre(n, cos(Theta));

    F[n] = pow((R/R0_mm), n);
    cosPhi[n] = cos(n*Phi);
    sinPhi[n] = sin(n*Phi);
  }
}

Siemens_B::~Siemens_B()
{
  int n;
  for (n = 0; n < coeffDim; n++)
    free(P[n]);

  free(P);
  free(F);
  free(cosPhi);
  free(sinPhi);
}

float Siemens_B::siemens_B_x(float **Alpha_x, float **Beta_x)
{
  return Siemens_B::siemens_B(Alpha_x, Beta_x);
}

float Siemens_B::siemens_B_y(float **Alpha_y, float **Beta_y)
{
  return Siemens_B::siemens_B(Alpha_y, Beta_y);
}

float Siemens_B::siemens_B_z(float **Alpha_z, float **Beta_z)
{
  return Siemens_B::siemens_B(Alpha_z, Beta_z);
}

float Siemens_B::siemens_B(float **Alpha, float **Beta)
{
  float B = 0;
  int n;
  for (n = 0; n <= nmax; n++)
  {
    int m;
    for (m = 0; m <= n; m++)
    {
      float F2 = Alpha[n+1][m+1] * cosPhi[m] + Beta[n+1][m+1] * sinPhi[m];
      B = B + F[n]*P[n][m+1]*F2;
    }
  }

  return B;
}

void Siemens_B::siemens_legendre(int n, double x)
{
  int m;
  for (m = 0; m <= n; m++)
  {
    P[n][m+1] = gsl_sf_legendre_Plm_e(n, m, x);
  }

  if (getenv("PRN_LEGENDRE"))
  {
    printf("\nlegendre (n=%d, x=%lf) = \n", n, x);
    for (m = 0; m <= n+1; m++)
      printf("\tP[%d][%d] = %lf\n", n, m, P[n][m]);
  }

  for (m = 1; m <= n; m++)
  {
    P[n][m+1] *= normfact[n][m];
#if 0
    double normfact = pow((-1), m+1) * sqrt((2*n+1)*factorial(n-m-1)/(2*factorial(n+m+1)));
    P[m+1] *= normfact;

    printf("normfact[%d][%d] = %lf, pow((-1), %d) * sqrt((%d)*factorial(%d)/(2*factorial(%d)))\n",
           n, m, normfact, m+1, 2*n+1, n-m-1, n+m+1);
#endif
  }

  if (getenv("PRN_LEGENDRE"))
  {
    printf("\nsiemens_legendre (n=%d, x=%lf) = \n", n, x);
    for (m = 0; m <= n+1; m++)
      printf("\tP[%d][%d] = %lf\n", n, m, P[n][m]);
  }
}
 
