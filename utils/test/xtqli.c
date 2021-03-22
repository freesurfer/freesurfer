/*
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


/* Driver for routine tqli */
//
// On Linux the following fails with -O for the eigenvalues at 5, 9, 10
// Make sure that you don't optimize but with -ffloat-store
//
#include <stdio.h>
#include <math.h>

#include "numerics.h"

#define NP 10
#define TINY 1.0e-6

int main(void)
{
  int i,j,k;
  float *d,*e,*f,**a;
  static float c[NP][NP]=
    {
      {
        5.0, 4.3, 3.0, 2.0, 1.0, 0.0,-1.0,-2.0,-3.0,-4.0
      },
      { 4.3, 5.1, 4.0, 3.0, 2.0, 1.0, 0.0,-1.0,-2.0,-3.0},
      { 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0,-1.0,-2.0},
      { 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0,-1.0},
      { 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0},
      { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0},
      {-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0},
      {-2.0,-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0},
      {-3.0,-2.0,-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0},
      {-4.0,-3.0,-2.0,-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0}
    };

  d=vector(1,NP);
  e=vector(1,NP);
  f=vector(1,NP);
  a=matrix(1,NP,1,NP);
  for (i=1;i<=NP;i++)
    for (j=1;j<=NP;j++) a[i][j]=c[i-1][j-1];
  tred2(a,NP,d,e);
  tqli(d,e,NP,a);
  printf("\nEigenvectors for a real symmetric matrix\n");
  for (i=1;i<=NP;i++)
  {
    for (j=1;j<=NP;j++)
    {
      f[j]=0.0;
      for (k=1;k<=NP;k++)
        f[j] += (c[j-1][k-1]*a[k][i]);
    }
    printf("%s %3d %s %10.6f\n","eigenvalue",i," =",d[i]);
    printf("%11s %14s %9s\n","vector","mtrx*vect.","ratio");
    for (j=1;j<=NP;j++)
    {
      if (fabs(a[j][i]) < TINY)
        printf("%12.6f %12.6f %12s\n",
               a[j][i],f[j],"div. by 0");
      else
        printf("%12.6f %12.6f %12.6f\n",
               a[j][i],f[j],f[j]/a[j][i]);
    }
    printf("Verify the last column values are the same as the eigenvalue.\n");
    printf("Press ENTER to continue...\n");
    (void) getchar();
  }
  free_matrix(a,1,NP,1,NP);
  free_vector(f,1,NP);
  free_vector(e,1,NP);
  free_vector(d,1,NP);
  return 0;
}
