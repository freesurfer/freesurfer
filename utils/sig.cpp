/**
 * @brief stats routines
 *
 * such things as FDR and beta and gamma incomplete functions
 */
/*
 * Original Author: Bruce Fischl and Doug Greve
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

#include <float.h>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "numerics.h"
#include "proto.h"
#include "sig.h"

/*----------------------------------------------------*/
/*!
  \fn double fdr2vwth(double *p, int np, double fdr)
  \brief returns the voxel-wise threshold given the
    False Discovery Rate. Complement to vwth2fdr().
  \param p - list of values between 0 and 1.
  \param np - number in p.
  \param fdr - False Discovery Rate, between 0 and 1.
  \return threshold between 0 and 1. If no values meet
    the FDR criterion, then a value slightly smaller
    than the smallest p is returned. This is assures
    that no voxels will be above/below threshold.

  Note: values in p will be changed (sorted ascending).
  Ref: http://www.sph.umich.edu/~nichols/FDR/FDR.m
  Thresholding of Statistical Maps in Functional Neuroimaging Using
  the False Discovery Rate.  Christopher R. Genovese, Nicole A. Lazar,
  Thomas E. Nichols (2002).  NeuroImage 15:870-878.
  ---------------------------------------------------------*/
double fdr2vwth(double *p, int np, double fdr)
{
  int n;
  double r;
  r = fdr / np;

  // Sort in ascending order
  qsort(p, np, sizeof(double), doublecompar);

  // Find the largest value of n for which p < r*(n+1)
  for (n = np - 1; n >= 0; n--)
    if (p[n] < r * (n + 1)) return (p[n]);

  // If p[n] never goes less than r*(n+1), then return a value
  // slightly smaller than the smallest value of p so that
  // nothing is above threshold
  return (0.9 * p[0]);
}

/*---------------------------------------------------------------
  vwth2fdr() - returns the False Discovery Rate given the voxel-wise
  threshold. Complement to fdr2vwth().
  p - list of values between 0 and 1.
  np - number in p.
  vwth - voxel-wise threshold between 0 and 1.
  Return: FDR between 0 and 1.
  Note: values in p will be changed (sorted ascending).
  ---------------------------------------------------------*/
double vwth2fdr(double *p, int np, double vwth)
{
  double fdr = 0, r = 0;
  int n;

  // Sort in ascending order
  qsort(p, np, sizeof(double), doublecompar);

  // Find the index n of the pvalue closest to vwth
  for (n = 1; n < np; n++)
    if (p[n - 1] <= vwth && vwth < p[n]) break;

  // Compute the corresponding FDR
  r = (double)n / np;  // Normalize rank for index n
  fdr = vwth / r;

  return (fdr);
}

/*---------------------------------------------------------------
  doublecompar() - qsort compare function to compare two doubles
  --------------------------------------------------------------*/
int doublecompar(const void *v1, const void *v2)
{
  double dv1, dv2;
  dv1 = *((double *)v1);
  dv2 = *((double *)v2);
  if (dv1 < dv2) return (-1);
  if (dv1 == dv2) return (0);
  return (+1);
}

#define MAXT 30.0
#define MAXDOF 200
/* note, everything here is Doug's fault */
double sigt(double t, int df)
{
  double sig = 0.0, sig1, sig2;

  if (df == 0)
    sig = 1;
  else if (df < MAXDOF && fabs(t) < MAXT)
    sig = OpenBetaIncomplete(0.5 * df, 0.5, df / (df + t * t));

  if (df > MAXDOF || sig < DBL_MIN) {
    sig1 = erfc(fabs(t) / sqrt(2.0));
    sig2 = 1.0 / sqrt(2.0 * M_PI) * exp(-0.5 * (t * t)); /* use Normal approximation */
    if (sig1 > sig2)
      sig = sig1;
    else
      sig = sig2;
  }

  if (!std::isfinite(sig)) printf("### Numerical error: sigt(%e,%d) = %e\n", t, df, sig);
  if (sig > 1.0) sig = 1.0;

  return sig;
}

float sigchisq(double chisq, int df)
{
  float p;

  if (FZERO(chisq) || df == 0) return (1.0f);
  p = OpenGammaIncomplete((float)df / 2, (float)chisq / 2.0);
  return (1 - p);
}
