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
#include <math.h>
#include <stdio.h>

#include "numerics.h"

#define NP   10
#define TINY 1.0e-6

#define NRANSI
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define sqr(a)     (((a) == 0.0) ? 0.0 : (a) * (a))

double pythag(double a, double b) {
  double absa, absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb)
    return absa * sqrt(1.0 + sqr(absb / absa));
  else
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + sqr(absa / absb)));
}

void tqli(float *d, float *e, int n, float **z) {
  int    m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  for (i = 1; i < n; i++)
    e[i - 1] = e[i];
  e[n - 1] = 0.0;
  for (l = 0; l < n; l++) {
    iter = 0;
    do {
      for (m = l; m < n - 1; m++) {
        dd = fabs(d[m]) + fabs(d[m + 1]);
        if ((double)(fabs(e[m]) + dd) == dd)
          break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("ERROR: Too many iterations in tqli\n");
          exit(1);
        }
        g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
        s = c = 1.0;
        p     = 0.0;
        for (i = m - 1; i >= l; i--) {
          f        = s * e[i];
          b        = c * e[i];
          e[i + 1] = (r = pythag(f, g));
          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;
            break;
          }
          s        = f / r;
          c        = g / r;
          g        = d[i + 1] - p;
          r        = (d[i] - g) * s + 2.0 * c * b;
          d[i + 1] = g + (p = s * r);
          g        = c * r - b;
          for (k = 0; k < n; k++) {
            f           = z[k][i + 1];
            z[k][i + 1] = s * z[k][i] + c * f;
            z[k][i]     = c * z[k][i] - s * f;
          }
        }
        if (r == 0.0 && i >= l)
          continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }
}

void tred2(float **a, int n, float *d, float *e) {
  int    l, k, j, i;
  double scale, hh, h, g, f;

  for (i = n - 1; i > 0; i--) {
    l = i - 1;
    h = scale = 0.0;
    if (l > 0) {
      for (k = 0; k <= l; k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i] = a[i][l];
      else {
        for (k = 0; k <= l; k++) {
          a[i][k] /= scale;
          h += a[i][k] * a[i][k];
        }
        f    = a[i][l];
        g    = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale * g;
        h -= f * g;
        a[i][l] = f - g;
        f       = 0.0;
        for (j = 0; j <= l; j++) {
          a[j][i] = a[i][j] / h;
          g       = 0.0;
          for (k = 0; k <= j; k++)
            g += a[j][k] * a[i][k];
          for (k = j + 1; k <= l; k++)
            g += a[k][j] * a[i][k];
          e[j] = g / h;
          f += e[j] * a[i][j];
        }
        hh = f / (h + h);
        for (j = 0; j <= l; j++) {
          f    = a[i][j];
          e[j] = g = e[j] - hh * f;
          for (k = 0; k <= j; k++)
            a[j][k] -= (f * e[k] + g * a[i][k]);
        }
      }
    } else
      e[i] = a[i][l];
    d[i] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;
  /* Contents of this loop can be omitted if eigenvectors not
                  wanted except for statement d[i]=a[i][i]; */
  for (i = 0; i < n; i++) {
    l = i;
    if (d[i] != 0.0) {
      for (j = 0; j < l; j++) {
        g = 0.0;
        for (k = 0; k < l; k++)
          g += a[i][k] * a[k][j];
        for (k = 0; k < l; k++)
          a[k][j] -= g * a[k][i];
      }
    }
    d[i]    = a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < l; j++)
      a[j][i] = a[i][j] = 0.0;
  }
}

int main(void) {
  int          i, j, k;
  float *      d, *e, *f, **a;
  static float c[NP][NP] = {
      {5.0, 4.3, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -4.0},
      {4.3, 5.1, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0},
      {3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0},
      {2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0},
      {1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0},
      {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0},
      {-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0},
      {-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0},
      {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0},
      {-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0}};

  d = vector(1, NP);
  e = vector(1, NP);
  f = vector(1, NP);
  a = matrix(1, NP, 1, NP);
  for (i = 1; i <= NP; i++)
    for (j = 1; j <= NP; j++)
      a[i][j] = c[i - 1][j - 1];
  tred2(a, NP, d, e);
  tqli(d, e, NP, a);
  printf("\nEigenvectors for a real symmetric matrix\n");
  for (i = 1; i <= NP; i++) {
    for (j = 1; j <= NP; j++) {
      f[j] = 0.0;
      for (k = 1; k <= NP; k++)
        f[j] += (c[j - 1][k - 1] * a[k][i]);
    }
    printf("%s %3d %s %10.6f\n", "eigenvalue", i, " =", d[i]);
    printf("%11s %14s %9s\n", "vector", "mtrx*vect.", "ratio");
    for (j = 1; j <= NP; j++) {
      if (fabs(a[j][i]) < TINY)
        printf("%12.6f %12.6f %12s\n", a[j][i], f[j], "div. by 0");
      else
        printf("%12.6f %12.6f %12.6f\n", a[j][i], f[j], f[j] / a[j][i]);
    }
    printf("Verify the last column values are the same as the eigenvalue.\n");
    printf("Press ENTER to continue...\n");
    (void)getchar();
  }
  free_matrix(a, 1, NP, 1, NP);
  free_vector(f, 1, NP);
  free_vector(e, 1, NP);
  free_vector(d, 1, NP);
  return 0;
}
