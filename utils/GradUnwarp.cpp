#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "GradUnwarp.h"
#include "legendre.h"

GradUnwarp::GradUnwarp()
{
  fgrad = NULL;

  coeffCount = 0;
  coeffDim = 0;

  R0 = 0;
  Alpha_x = NULL; Alpha_y = NULL; Alpha_z = NULL;
  Beta_x  = NULL; Beta_y  = NULL; Beta_z  = NULL;

  P = NULL;
}

GradUnwarp::~GradUnwarp()
{
  //if (P != NULL)
  //  free(P);

  int i = 0;
  for (; i < coeffDim; i++)
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
  {
    free(normfact[i]);
    free(P[i]);
  }
  free(normfact);
  free(P);
}

void GradUnwarp::setup()
{
}

void GradUnwarp::list_coeff_files()
{
}

void GradUnwarp::read_siemens_coeff(const char *gradfilename)
{
  // check if gradfile has extension .grad

  // open gradfilename
  fgrad = fopen(gradfilename, "r");
  if (fgrad == NULL)
  {
    printf("ERROR: could not read the coefficient file %s\n", gradfilename);
    return;
  }

  printf("==> reading coefficients from gradient coil file %s\n" , gradfilename) ;

  _skipCoeffComment();

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
  while (!feof(fgrad))
  {
    fgets(coeffline, sizeof(coeffline), fgrad);

    int len = strlen(coeffline);
    char* ptr = coeffline;
    char* endptr = ptr + len;

    // skip leading spaces and tabs, also empty lines
    while (ptr != endptr && (*ptr == ' ' || *ptr == '\t' || *ptr == '\n'))
      ptr++;

    if (*ptr == '\0')
      continue;

    sscanf(ptr, "%d %c(%d, %d) %f %c", 
                      &coeff[coeffCount].num, &coeff[coeffCount].A_or_B, &coeff[coeffCount].n, &coeff[coeffCount].m, 
                      &coeff[coeffCount].value, &coeff[coeffCount].xyz);
    int larger = (coeff[coeffCount].n > coeff[coeffCount].m) ? (coeff[coeffCount].n+1) : (coeff[coeffCount].m+1);
    coeffDim = (coeffDim < larger) ? larger : coeffDim; 
    //printf("%d %c (%d, %d) %f %c\n", coeff[coeffCount].num, 
    //       coeff[coeffCount].A_or_B, coeff[coeffCount].n, coeff[coeffCount].m, coeff[coeffCount].value, coeff[coeffCount].xyz);
    coeffCount++;
  }

  fclose(fgrad);    

  _initCoeff();

  /**************************************************************************/
  /****** organize coefficient values ******/
  int n = 0;
  for (; n < coeffCount; n++)
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
    arrPtr[row][col] = coeff[n].value;
  }
}

float GradUnwarp::siemens_B_x()
{
  return GradUnwarp::siemens_B(Alpha_x, Beta_x, F2_x);
}

float GradUnwarp::siemens_B_y()
{
  return GradUnwarp::siemens_B(Alpha_y, Beta_y, F2_y);
}

float GradUnwarp::siemens_B_z()
{
  return GradUnwarp::siemens_B(Alpha_z, Beta_z, F2_z);
}

void GradUnwarp::siemens_B0(float X, float Y, float Z)
{
  int n;
  for (n = 0; n < coeffDim; n++)
  {
    memset(P[n], 0, coeffDim*sizeof(double));
    memset(F2_x[n], 0, coeffDim*sizeof(float));
    memset(F2_y[n], 0, coeffDim*sizeof(float));
    memset(F2_z[n], 0, coeffDim*sizeof(float));
  }

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
  int nmax = coeffDim - 1;

  for (n = 0; n < nmax+1; n++)
  {
    siemens_legendre(n, cos(Theta), P);

    F[n] = pow((R/R0), n);
    cosPhi[n] = cos(n*Phi);
    sinPhi[n] = sin(n*Phi);
  }

  for (n = 0; n < nmax+1; n++)
  {
    int m = 0;
    for (; m < n+1; m++)
    {
      F2_x[n][m] = Alpha_x[n][m] * cosPhi[m] + Beta_x[n][m] * sinPhi[m];
      F2_y[n][m] = Alpha_y[n][m] * cosPhi[m] + Beta_y[n][m] * sinPhi[m];
      F2_z[n][m] = Alpha_z[n][m] * cosPhi[m] + Beta_z[n][m] * sinPhi[m];
    }
  }
}

float GradUnwarp::siemens_B(float **Alpha, float **Beta, float **F2)
{
  int nmax = coeffDim - 1;
  float B = 0;
  int n = 0;
  for (; n < nmax+1; n++)
  {
    int m = 0;
    for (; m < n+1; m++)
    {
      //float F2 = Alpha[n][m] * cos(m*Phi) + Beta[n][m] * sin(m*Phi);
      //float F2 = Alpha[n][m] * cosPhi[m] + Beta[n][m] * sinPhi[m];
      //B = B + F[n]*P[n][m]*F2;
      B = B + F[n]*P[n][m]*F2[n][m];
    }
  }

  return B;
}

void GradUnwarp::initSiemensLegendreNormfact()
{
  // initialize variables to pre-calculate normfact for siemens_legendre()
  minusonepow = new double[coeffDim];
  factorials  = new double[2*coeffDim];

  normfact = new double*[coeffDim];
  P = new double*[coeffDim];
  F2_x = new float*[coeffDim];
  F2_y = new float*[coeffDim];
  F2_z = new float*[coeffDim];
  int n;
  for (n = 0; n < coeffDim; n++)
  {
    normfact[n] = new double[coeffDim];
    P[n] = new double[coeffDim];

    F2_x[n] = new float[coeffDim];
    F2_y[n] = new float[coeffDim];
    F2_z[n] = new float[coeffDim];
  }

  F = new float[coeffDim];
  cosPhi = new double[coeffDim];
  sinPhi = new double[coeffDim];

  // pre-calculate minusonepow, factorials, & normfact
  for (n=0; n < coeffDim; n++)
  {
    minusonepow[n] = pow((-1), n);
  }

  for (n=0; n < 2*coeffDim; n++)
  {
    factorials[n] = factorial(n);
  }

#if PRN_LEGENDRE
  printf("\n");
#endif
  for (n=0; n < coeffDim; n++)
  {
    int m = 0;
    for (; m < n; m++)
    {
      normfact[n][m] = minusonepow[m+1] * sqrt((2*n+1)*factorials[n-m-1]/(2*factorials[n+m+1]));
   
#if PRN_LEGENDRE
      printf("normfact2[%2d][%2d] = %s%.6lf, pow((-1), %2d) * sqrt((%2d)*factorial(%2d)/(2*factorial(%2d)))\n",
             n, m, (normfact[n][m] > 0) ? " " : "", normfact[n][m], m+1, 2*n+1, n-m-1, n+m+1);
#endif
    }
#if PRN_LEGENDRE
    printf("\n");
#endif
  }
}

void GradUnwarp::siemens_legendre(int n, double x, double **P)
{
  int m = 0;
  for (; m < n+1; m++)
  {
    P[n][m] = gsl_sf_legendre_Plm_e(n, m, x);
  }

#if PRN_LEGENDRE
  printf("\nlegendre (n=%d) = \n", n);
  for (m = 0; m < n+1; m++)
  {
    printf("\tP[%d][%d] = %lf\n", n, m, P[n][m]);
  }
#endif

  m = 0;
  for (; m < n; m++)
  {
    P[n][m+1] *= normfact[n][m];
#if 0
    double normfact = pow((-1), m+1) * sqrt((2*n+1)*factorial(n-m-1)/(2*factorial(n+m+1)));
    P[m+1] *= normfact;

    printf("normfact[%d][%d] = %lf, pow((-1), %d) * sqrt((%d)*factorial(%d)/(2*factorial(%d)))\n",
           n, m, normfact, m+1, 2*n+1, n-m-1, n+m+1);
#endif
  }

#if PRN_LEGENDRE
  printf("\nsiemens_legendre (n=%d) = \n", n);
  for (m = 0; m < n+1; m++)
  {
    printf("\tP[%d][%d] = %lf\n", n, m, P[n][m]);
  }
#endif
}
 
void GradUnwarp::spharm_evaluate(float X, float Y, float Z, float *Dx, float *Dy, float *Dz)
{
  siemens_B0(X, Y, Z);

  float bx = siemens_B_x();
  float by = siemens_B_y();
  float bz = siemens_B_z();

  *Dx = bx * R0;
  *Dy = by * R0;
  *Dz = bz * R0;
}

void GradUnwarp::unwarp()
{
}

void GradUnwarp::unwap_volume()
{
}

void GradUnwarp::unwarp_surface()
{
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
  Alpha_x = new float*[coeffDim];
  Alpha_y = new float*[coeffDim];
  Alpha_z = new float*[coeffDim];
  Beta_x  = new float*[coeffDim];
  Beta_y  = new float*[coeffDim];
  Beta_z  = new float*[coeffDim];

  int i = 0;
  for (; i < coeffDim; i++)
  {
    Alpha_x[i] = new float[coeffDim];
    Alpha_y[i] = new float[coeffDim];
    Alpha_z[i] = new float[coeffDim];
    Beta_x[i]  = new float[coeffDim];
    Beta_y[i]  = new float[coeffDim];
    Beta_z[i]  = new float[coeffDim];
  }
}

