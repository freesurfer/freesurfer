/**
 * @brief Compute B-Spline Coefficients and Interpolate MRI
 *
 * This Code is based on:
 *		P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
 *		IEEE Transactions on Medical Imaging,
 *		vol. 19, no. 7, pp. 739-758, July 2000.
 *
 * and code (for 2D images) obtained from
 *   http://bigwww.epfl.ch/thevenaz/interpolation/
 *   http://bigwww.epfl.ch/sage/pyramids/
 */
/*
 * Original Author: Martin Reuter
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
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "romp_support.h"

#include "error.h"
#include "macros.h"
#include "mriBSpline.h"

static double InitialCausalCoefficient(double c[],      /* coefficients */
                                       long DataLength, /* number of coefficients */
                                       double z,        /* actual pole */
                                       double Tolerance /* admissible relative error */
)

{ /* begin InitialCausalCoefficient */

  double Sum, zn, z2n, iz;
  long n, Horizon;

  /* this initialization corresponds to mirror boundaries */
  Horizon = DataLength;
  if (Tolerance > 0.0) {
    Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
  }
  if (Horizon < DataLength) {
    /* accelerated loop */
    zn = z;
    Sum = c[0];
    for (n = 1L; n < Horizon; n++) {
      Sum += zn * c[n];
      zn *= z;
    }
    return (Sum);
  }
  else {
    /* full loop */
    zn = z;
    iz = 1.0 / z;
    z2n = pow(z, (double)(DataLength - 1L));
    Sum = c[0] + z2n * c[DataLength - 1L];
    z2n *= z2n * iz;
    for (n = 1L; n <= DataLength - 2L; n++) {
      Sum += (zn + z2n) * c[n];
      zn *= z;
      z2n *= iz;
    }
    return (Sum / (1.0 - zn * zn));
  }
} /* end InitialCausalCoefficient */

static double InitialAntiCausalCoefficient(double c[],      /* coefficients */
                                           long DataLength, /* number of samples or coefficients */
                                           double z         /* actual pole */
)

{ /* begin InitialAntiCausalCoefficient */

  /* this initialization corresponds to mirror boundaries */
  return ((z / (z * z - 1.0)) * (z * c[DataLength - 2L] + c[DataLength - 1L]));
} /* end InitialAntiCausalCoefficient */

static void ConvertToInterpolationCoefficients(double c[],      /* input samples --> output coefficients */
                                               long DataLength, /* number of samples or coefficients */
                                               double z[],      /* poles */
                                               long NbPoles,    /* number of poles */
                                               double Tolerance /* admissible relative error */
)

{ /* begin ConvertToInterpolationCoefficients */

  double Lambda = 1.0;
  long n, k;

  /* special case required by mirror boundaries */
  if (DataLength == 1L) {
    return;
  }
  /* compute the overall gain */
  for (k = 0L; k < NbPoles; k++) {
    Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
  }
  /* apply the gain */
  for (n = 0L; n < DataLength; n++) {
    c[n] *= Lambda;
  }
  /* loop over all poles */
  for (k = 0L; k < NbPoles; k++) {
    /* causal initialization */
    c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
    /* causal recursion */
    for (n = 1L; n < DataLength; n++) {
      c[n] += z[k] * c[n - 1L];
    }
    /* anticausal initialization */
    c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
    /* anticausal recursion */
    for (n = DataLength - 2L; 0 <= n; n--) {
      c[n] = z[k] * (c[n + 1L] - c[n]);
    }
  }
} /* end ConvertToInterpolationCoefficients */

// get and set line routines need to be improved (speedup by grabbing a
// whole row in data block...)
static void getZLine(const MRI *mri, int x, int y, int f, double Line[])
{
  int z;
  for (z = 0; z < mri->depth; z++) Line[z] = MRIgetVoxVal(mri, x, y, z, f);
}
static void getZLineStop(const MRI *mri, int x, int y, int f, double Line[], int zstop)
{
  int z;
  for (z = 0; z < zstop; z++) Line[z] = MRIgetVoxVal(mri, x, y, z, f);
}
static void setZLine(MRI *mri, int x, int y, int f, double Line[])
{
  int z;
  if (x >= mri->width) {
    printf("X limit error\n");
    exit(1);
  }
  if (y >= mri->height) {
    printf("Y limit error\n");
    exit(1);
  }
  if (f >= mri->nframes) {
    printf("Frames limit error\n");
    exit(1);
  }
  for (z = 0; z < mri->depth; z++) {
    if (mri->type == MRI_UCHAR) {
      if (Line[z] > 255)
        MRIsetVoxVal(mri, x, y, z, f, 255);
      else if (Line[z] < 0)
        MRIsetVoxVal(mri, x, y, z, f, 0);
      else
        MRIsetVoxVal(mri, x, y, z, f, nint(Line[z]));
    }
    else
      MRIsetVoxVal(mri, x, y, z, f, (float)Line[z]);
  }
}
static void getYLine(const MRI *mri, int x, int z, int f, double Line[])
{
  int y;
  for (y = 0; y < mri->height; y++) Line[y] = MRIgetVoxVal(mri, x, y, z, f);
}
static void getYLineStop(const MRI *mri, int x, int z, int f, double Line[], int ystop)
{
  int y;
  for (y = 0; y < ystop; y++) Line[y] = MRIgetVoxVal(mri, x, y, z, f);
}
static void setYLine(MRI *mri, int x, int z, int f, double Line[])
{
  int y;
  if (x >= mri->width) {
    printf("X limit error\n");
    exit(1);
  }
  if (z >= mri->depth) {
    printf("Z limit error\n");
    exit(1);
  }
  if (f >= mri->nframes) {
    printf("Frames limit error\n");
    exit(1);
  }
  for (y = 0; y < mri->height; y++) {
    if (mri->type == MRI_UCHAR) {
      if (Line[y] > 255)
        MRIsetVoxVal(mri, x, y, z, f, 255);
      else if (Line[y] < 0)
        MRIsetVoxVal(mri, x, y, z, f, 0);
      else
        MRIsetVoxVal(mri, x, y, z, f, nint(Line[y]));
    }
    else
      MRIsetVoxVal(mri, x, y, z, f, (float)Line[y]);
  }
}
static void getXLine(const MRI *mri, int y, int z, int f, double Line[])
{
  int x;
  for (x = 0; x < mri->width; x++) Line[x] = MRIgetVoxVal(mri, x, y, z, f);
}
// static void	getXLineStop (const MRI	* mri, int	y, int z, int f, double Line[], int xstop)
// {
// 	int	x;
//   for (x=0;x<xstop;x++) Line[x] = MRIgetVoxVal(mri, x, y, z, f);
// }
static void setXLine(MRI *mri, int y, int z, int f, double Line[])
{
  int x;
  if (y >= mri->height) {
    printf("Y limit error\n");
    exit(1);
  }
  if (z >= mri->depth) {
    printf("Z limit error\n");
    exit(1);
  }
  if (f >= mri->nframes) {
    printf("Frames limit error\n");
    exit(1);
  }
  for (x = 0; x < mri->width; x++) {
    // printf("x %i\n",x);
    if (mri->type == MRI_UCHAR) {
      if (Line[x] > 255)
        MRIsetVoxVal(mri, x, y, z, f, 255);
      else if (Line[x] < 0)
        MRIsetVoxVal(mri, x, y, z, f, 0);
      else
        MRIsetVoxVal(mri, x, y, z, f, nint(Line[x]));
    }
    else
      MRIsetVoxVal(mri, x, y, z, f, (float)Line[x]);
  }
}

static int initPoles(double Pole[4], int SplineDegree)
{
  int NbPoles = 0;
  /* recover the poles from a lookup table */
  switch (SplineDegree) {
    case 2:
      NbPoles = 1;
      Pole[0] = sqrt(8.0) - 3.0;
      break;
    case 3:
      NbPoles = 1;
      Pole[0] = sqrt(3.0) - 2.0;
      break;
    case 4:
      NbPoles = 2;
      Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
      Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
      break;
    case 5:
      NbPoles = 2;
      Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0) - 13.0 / 2.0;
      Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0) - 13.0 / 2.0;
      break;
    case 6:
      NbPoles = 3;
      Pole[0] = -0.48829458930304475513011803888378906211227916123938;
      Pole[1] = -0.081679271076237512597937765737059080653379610398148;
      Pole[2] = -0.0014141518083258177510872439765585925278641690553467;
      break;
    case 7:
      NbPoles = 3;
      Pole[0] = -0.53528043079643816554240378168164607183392315234269;
      Pole[1] = -0.12255461519232669051527226435935734360548654942730;
      Pole[2] = -0.0091486948096082769285930216516478534156925639545994;
      break;
    case 8:
      NbPoles = 4;
      Pole[0] = -0.57468690924876543053013930412874542429066157804125;
      Pole[1] = -0.16303526929728093524055189686073705223476814550830;
      Pole[2] = -0.023632294694844850023403919296361320612665920854629;
      Pole[3] = -0.00015382131064169091173935253018402160762964054070043;
      break;
    case 9:
      NbPoles = 4;
      Pole[0] = -0.60799738916862577900772082395428976943963471853991;
      Pole[1] = -0.20175052019315323879606468505597043468089886575747;
      Pole[2] = -0.043222608540481752133321142979429688265852380231497;
      Pole[3] = -0.0021213069031808184203048965578486234220548560988624;
      break;
    default:
      printf("Invalid spline degree\n");
      exit(1);
  }
  return NbPoles;
}

static void BsplineWeights(int SplineDegree, const double *val, const int *index, double *weights)
{
  double w, w2, w4, t, t0, t1;
  switch (SplineDegree) {
    case 2:
      w = *val - (double)index[1];
      weights[1] = 3.0 / 4.0 - w * w;
      weights[2] = (1.0 / 2.0) * (w - weights[1] + 1.0);
      weights[0] = 1.0 - weights[1] - weights[2];
      break;
    case 3:
      w = *val - (double)index[1];
      weights[3] = (1.0 / 6.0) * w * w * w;
      weights[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - weights[3];
      weights[2] = w + weights[0] - 2.0 * weights[3];
      weights[1] = 1.0 - weights[0] - weights[2] - weights[3];
      break;
    case 4:
      w = *val - (double)index[2];
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      weights[0] = 1.0 / 2.0 - w;
      weights[0] *= weights[0];
      weights[0] *= (1.0 / 24.0) * weights[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      weights[1] = t1 + t0;
      weights[3] = t1 - t0;
      weights[4] = weights[0] + t0 + (1.0 / 2.0) * w;
      weights[2] = 1.0 - weights[0] - weights[1] - weights[3] - weights[4];
      break;
    case 5:
      w = *val - (double)index[2];
      w2 = w * w;
      weights[5] = (1.0 / 120.0) * w * w2 * w2;
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      weights[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - weights[5];
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      weights[2] = t0 + t1;
      weights[3] = t0 - t1;
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      weights[1] = t0 + t1;
      weights[4] = t0 - t1;
      break;
    case 6:
      w = *val - (double)index[3];
      weights[0] = 1.0 / 2.0 - w;
      weights[0] *= weights[0] * weights[0];
      weights[0] *= weights[0] / 720.0;
      weights[1] =
          (361.0 / 192.0 -
           w * (59.0 / 8.0 + w * (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w) * (1.0 / 2.0 + w))))) /
          120.0;
      weights[2] = (10543.0 / 960.0 +
                    w * (-289.0 / 16.0 + w * (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w * (-1.0 + w)))))) /
                   48.0;
      w2 = w * w;
      weights[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2 * (21.0 / 4.0 - w2))) / 36.0;
      weights[4] = (10543.0 / 960.0 +
                    w * (289.0 / 16.0 + w * (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w * (1.0 + w)))))) /
                   48.0;
      weights[6] = 1.0 / 2.0 + w;
      weights[6] *= weights[6] * weights[6];
      weights[6] *= weights[6] / 720.0;
      weights[5] = 1.0 - weights[0] - weights[1] - weights[2] - weights[3] - weights[4] - weights[6];
      break;
    case 7:
      w = *val - (double)index[3];
      weights[0] = 1.0 - w;
      weights[0] *= weights[0];
      weights[0] *= weights[0] * weights[0];
      weights[0] *= (1.0 - w) / 5040.0;
      w2 = w * w;
      weights[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w * (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
      weights[2] = (397.0 / 7.0 -
                    w * (245.0 / 3.0 + w * (-15.0 + w * (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w * (-5.0 + w))))))) /
                   240.0;
      weights[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2 * (-4.0 + w)))) / 144.0;
      weights[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w * (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
      weights[5] =
          (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w * (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
      weights[7] = w2;
      weights[7] *= weights[7] * weights[7];
      weights[7] *= w / 5040.0;
      weights[6] = 1.0 - weights[0] - weights[1] - weights[2] - weights[3] - weights[4] - weights[5] - weights[7];
      break;
    case 8:
      w = *val - (double)index[4];
      weights[0] = 1.0 / 2.0 - w;
      weights[0] *= weights[0];
      weights[0] *= weights[0];
      weights[0] *= weights[0] / 40320.0;
      w2 = w * w;
      weights[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2))) *
                   (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w * (-3.0 + w)))) / 5040.0;
      weights[2] =
          (82903.0 / 1792.0 +
           w * (-4177.0 / 32.0 +
                w * (2275.0 / 16.0 +
                     w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w * (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) /
          1440.0;
      weights[3] =
          (310661.0 / 1792.0 -
           w * (14219.0 / 64.0 +
                w * (-199.0 / 8.0 +
                     w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w * (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) /
          720.0;
      weights[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2 * (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
      weights[5] =
          (310661.0 / 1792.0 -
           w * (-14219.0 / 64.0 +
                w * (-199.0 / 8.0 +
                     w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w * (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) /
          720.0;
      weights[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2))) *
                   (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w * (3.0 + w)))) / 5040.0;
      weights[8] = 1.0 / 2.0 + w;
      weights[8] *= weights[8];
      weights[8] *= weights[8];
      weights[8] *= weights[8] / 40320.0;
      weights[6] =
          1.0 - weights[0] - weights[1] - weights[2] - weights[3] - weights[4] - weights[5] - weights[7] - weights[8];
      break;
    case 9:
      w = *val - (double)index[4];
      weights[0] = 1.0 - w;
      weights[0] *= weights[0];
      weights[0] *= weights[0];
      weights[0] *= weights[0] * (1.0 - w) / 362880.0;
      weights[1] =
          (502.0 / 9.0 +
           w * (-246.0 + w * (472.0 + w * (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 +
                                                                                   w * (24.0 + w * (-8.0 + w))))))))) /
          40320.0;
      weights[2] =
          (3652.0 / 9.0 -
           w * (2023.0 / 2.0 +
                w * (-952.0 + w * (938.0 / 3.0 +
                                   w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w * (14.0 + w * (-7.0 + w))))))))) /
          10080.0;
      weights[3] =
          (44117.0 / 42.0 +
           w * (-2427.0 / 2.0 +
                w * (66.0 + w * (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w * (-6.0 + w))))))))) /
          4320.0;
      w2 = w * w;
      weights[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2 * (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
      weights[5] =
          (44117.0 / 63.0 +
           w * (809.0 + w * (44.0 + w * (-868.0 / 3.0 +
                                         w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w * (-4.0 + w * (-4.0 + w))))))))) /
          2880.0;
      weights[6] =
          (3652.0 / 21.0 -
           w * (-867.0 / 2.0 +
                w * (-408.0 + w * (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w) * (2.0 + w))))))) /
          4320.0;
      weights[7] =
          (251.0 / 18.0 +
           w * (123.0 / 2.0 +
                w * (118.0 +
                     w * (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w * (-6.0 + w * (-2.0 + w))))))))) /
          10080.0;
      weights[9] = w2 * w2;
      weights[9] *= weights[9] * w / 362880.0;
      weights[8] = 1.0 - weights[0] - weights[1] - weights[2] - weights[3] - weights[4] - weights[5] - weights[6] -
                   weights[7] - weights[9];
      break;
    default:
      printf("Invalid spline degree\n");
      exit(1);
  }
}

/*--------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------
        Function: 	PyramidFilterSplinel2

        Purpose:  	Initializes down- and up-sampling filter arrays for
                                least squares splines of order 0 to 3.  (little l_2 norm)
                                 g : reduce filter
                                 h : expand filter

        Author:		Michael Unser, NIH, BEIP, May 1992

---------------------------------------------------------------------------- */
static void PyramidFilterSplinel2(double g[], long *ng, double *h, long *nh, long Order)
{
  switch (Order) {
    case 0L:
      *ng = 1L;
      *nh = 1L;
      break;

    case 1L:
      g[0] = 0.707107;
      g[1] = 0.292893;
      g[2] = -0.12132;
      g[3] = -0.0502525;
      g[4] = 0.0208153;
      g[5] = 0.00862197;
      g[6] = -0.00357134;
      g[7] = -0.0014793;
      g[8] = 0.000612745;
      *ng = 9L;
      h[0] = 1.;
      h[1] = 0.5;
      *nh = 2L;
      break;
    case 2L:
      g[0] = 0.617317;
      g[1] = 0.310754;
      g[2] = -0.0949641;
      g[3] = -0.0858654;
      g[4] = 0.0529153;
      g[5] = 0.0362437;
      g[6] = -0.0240408;
      g[7] = -0.0160987;
      g[8] = 0.0107498;
      g[9] = 0.00718418;
      g[10] = -0.00480004;
      g[11] = -0.00320734;
      g[12] = 0.00214306;
      g[13] = 0.00143195;
      g[14] = -0.0009568;
      g[15] = -0.000639312;
      *ng = 16L;
      h[0] = 1.;
      h[1] = 0.585786;
      h[2] = 0;
      h[3] = -0.100505;
      h[4] = 0;
      h[5] = 0.0172439;
      h[6] = 0;
      h[7] = -0.00295859;
      h[8] = 0;
      h[9] = 0.000507614;
      *nh = 10L;
      break;
    case 3L:
      g[0] = 0.596797;
      g[1] = 0.313287;
      g[2] = -0.0827691;
      g[3] = -0.0921993;
      g[4] = 0.0540288;
      g[5] = 0.0436996;
      g[6] = -0.0302508;
      g[7] = -0.0225552;
      g[8] = 0.0162251;
      g[9] = 0.0118738;
      g[10] = -0.00861788;
      g[11] = -0.00627964;
      g[12] = 0.00456713;
      g[13] = 0.00332464;
      g[14] = -0.00241916;
      g[15] = -0.00176059;
      g[16] = 0.00128128;
      g[17] = 0.000932349;
      g[18] = -0.000678643;
      g[19] = -0.000493682;
      *ng = 20L;
      h[0] = 1.;
      h[1] = 0.600481;
      h[2] = 0;
      h[3] = -0.127405;
      h[4] = 0;
      h[5] = 0.034138;
      h[6] = 0;
      h[7] = -0.00914725;
      h[8] = 0;
      h[9] = 0.002451;
      h[10] = 0;
      h[11] = -0.000656743;
      *nh = 12L;
      break;
    default:
      *ng = -1L;
      *nh = -1L;
      break;
  }
}

/* ----------------------------------------------------------------------------
        Function: 	PyramidFilterSplineL2

        Purpose:  	Initializes down- and up-sampling filter arrays for
                                L2 spline pyramid of order 0 to 5.
                                 g : reduce filter
                                 h : expand filter

        Author:		Michael Unser, NIH, BEIP, May 1992

---------------------------------------------------------------------------- */

static void PyramidFilterSplineL2(double g[], long *ng, double h[], long *nh, long Order)
{
  switch (Order) {
    case 0L:
      *ng = 1L;
      *nh = 1L;
      break;

    case 1L:
      g[0] = 0.683013;
      g[1] = 0.316987;
      g[2] = -0.116025;
      g[3] = -0.0849365;
      g[4] = 0.0310889;
      g[5] = 0.0227587;
      g[6] = -0.00833025;
      g[7] = -0.00609817;
      g[8] = 0.00223208;
      g[9] = 0.001634;
      g[10] = -0.000598085;
      g[11] = -0.000437829;
      g[12] = 0.000160256;
      g[13] = 0.000117316;
      *ng = 14L;
      h[0] = 1.;
      h[1] = 0.5;
      *nh = 2L;
      break;

    case 3L:
      g[0] = 0.594902;
      g[1] = 0.31431;
      g[2] = -0.0816632;
      g[3] = -0.0942586;
      g[4] = 0.0541374;
      g[5] = 0.0454105;
      g[6] = -0.0307778;
      g[7] = -0.0236728;
      g[8] = 0.0166858;
      g[9] = 0.0125975;
      g[10] = -0.00895838;
      g[11] = -0.00673388;
      g[12] = 0.00479847;
      g[13] = 0.00360339;
      g[14] = -0.00256892;
      g[15] = -0.00192868;
      g[16] = 0.00137514;
      g[17] = 0.00103237;
      g[18] = -0.000736093;
      g[19] = -0.000552606;
      g[20] = 0.000394017;
      g[21] = 0.000295799;
      g[22] = -0.00021091;
      g[23] = -0.000158335;
      g[24] = 0.000112896;
      *ng = 25L;
      h[0] = 1.;
      h[1] = 0.600481;
      h[2] = 0.0;
      h[3] = -0.127405;
      h[4] = 0;
      h[5] = 0.034138;
      h[6] = 0;
      h[7] = -0.00914725;
      h[8] = 0;
      h[9] = 0.002451;
      h[10] = 0;
      h[11] = -0.000656743;
      *nh = 12L;
      break;

    case 5L:
      g[0] = 0.564388;
      g[1] = 0.316168;
      g[2] = -0.0597634;
      g[3] = -0.0998708;
      g[4] = 0.0484525;
      g[5] = 0.0539099;
      g[6] = -0.0355614;
      g[7] = -0.033052;
      g[8] = 0.0246347;
      g[9] = 0.0212024;
      g[10] = -0.0166097;
      g[11] = -0.0138474;
      g[12] = 0.0110719;
      g[13] = 0.00911006;
      g[14] = -0.00734567;
      g[15] = -0.0060115;
      g[16] = 0.00486404;
      g[17] = 0.00397176;
      g[18] = -0.00321822;
      g[19] = -0.00262545;
      g[20] = 0.00212859;
      g[21] = 0.00173587;
      g[22] = -0.0014077;
      g[23] = -0.0011478;
      g[24] = 0.000930899;
      g[25] = 0.000758982;
      g[26] = -0.000615582;
      g[27] = -0.000501884;
      g[28] = 0.000407066;
      g[29] = 0.000331877;
      g[30] = -0.00026918;
      g[31] = -0.000219459;
      g[32] = 0.000178;
      g[33] = 0.00014512;
      g[34] = -0.000117706;
      *ng = 35L;
      h[0] = 1.;
      h[1] = 0.619879;
      h[2] = 0.0;
      h[3] = -0.167965;
      h[4] = 0;
      h[5] = 0.0686374;
      h[6] = 0;
      h[7] = -0.0293948;
      h[8] = 0.0;
      h[9] = 0.0126498;
      h[10] = 0;
      h[11] = -0.00544641;
      h[12] = 0.0;
      h[13] = 0.00234508;
      h[14] = 0;
      h[15] = -0.00100973;
      h[16] = 0.0;
      h[17] = 0.000434766;
      h[18] = 0;
      h[19] = -0.000187199;
      *nh = 20L;
      break;

    default:
      *ng = -1L;
      *nh = -1L;
      printf("Spline filters only defined for n=0,1,3,5\n");
      exit(1);

      break;
  }
}

/* ----------------------------------------------------------------------------
        Function:	PyramidFilterCentered

        Purpose:	Initializes down- and up-sampling filter arrays for
                                least squares CENTERED splines of order 0 to 4.  (little l_2 norm)
                                        g : reduce filter
                                        h : expand filter

        Note:		filter arrays should be defined as
                                        float g[20],h[20] filter arrays
                                        short *ng,*nh;	number of taps
                                        short Order;	order of the spline

        Author:		Patrick Brigger, NIH, BEIP	May 1996
                                Daniel Sage, EPFL, Biomedical Imaging Group, November 1999

---------------------------------------------------------------------------- */
static void PyramidFilterCentered(double g[], long *ng, double h[], long *nh, long Order)
{
  switch (Order) {
    case 0:
      g[0] = 1;
      *ng = 1;
      h[0] = 2;
      *nh = 1;
      break;

    case 1:
      g[0] = 1.;
      g[1] = 0.333333;
      g[2] = -0.333333;
      g[3] = -0.111111;
      g[4] = 0.111111;
      g[5] = 0.037037;
      g[6] = -0.037037;
      g[7] = -0.0123457;
      g[8] = 0.0123457;
      g[9] = 0.00411523;
      g[10] = -0.00411523;
      g[11] = -0.00137174;
      g[12] = 0.00137174;
      g[13] = 0.000457247;
      g[14] = -0.000457247;
      g[15] = -0.000152416;
      g[16] = 0.000152416;
      g[17] = 0.0000508053;
      g[18] = -0.0000508053;
      g[19] = -0.0000169351;
      g[20] = 0.0000169351;
      *ng = 21;
      h[0] = 1;
      h[1] = 0.5;
      *nh = 2;
      break;

    case 2:
      g[0] = 0.738417;
      g[1] = 0.307916;
      g[2] = -0.171064;
      g[3] = -0.0799199;
      g[4] = 0.0735791;
      g[5] = 0.03108;
      g[6] = -0.0307862;
      g[7] = -0.0128561;
      g[8] = 0.0128425;
      g[9] = 0.00535611;
      g[10] = -0.00535548;
      g[11] = -0.00223325;
      g[12] = 0.00223322;
      g[13] = 0.000931242;
      g[14] = -0.00093124;
      g[15] = -0.000388322;
      g[16] = 0.000388322;
      g[17] = 0.000161928;
      g[18] = -0.000161928;
      g[19] = -0.0000675233;
      g[20] = 0.0000675233;
      *ng = 21;
      h[0] = 1.20711;
      h[1] = 0.585786;
      h[2] = -0.12132;
      h[3] = -0.100505;
      h[4] = 0.0208153;
      h[5] = 0.0172439;
      h[6] = -0.00357134;
      h[7] = -0.00295859;
      h[8] = 0.000612745;
      h[9] = 0.000507614;
      h[10] = -0.00010513;
      *nh = 11;
      break;

    case 3:
      g[0] = 0.708792;
      g[1] = 0.328616;
      g[2] = -0.165157;
      g[3] = -0.114448;
      g[4] = 0.0944036;
      g[5] = 0.0543881;
      g[6] = -0.05193;
      g[7] = -0.0284868;
      g[8] = 0.0281854;
      g[9] = 0.0152877;
      g[10] = -0.0152508;
      g[11] = -0.00825077;
      g[12] = 0.00824629;
      g[13] = 0.00445865;
      g[14] = -0.0044582;
      g[15] = -0.00241009;
      g[16] = 0.00241022;
      g[17] = 0.00130278;
      g[18] = -0.00130313;
      g[19] = -0.000704109;
      g[20] = 0.000704784;
      *ng = 21;
      h[0] = 1.13726;
      h[1] = 0.625601;
      h[2] = -0.0870191;
      h[3] = -0.159256;
      h[4] = 0.0233167;
      h[5] = 0.0426725;
      h[6] = -0.00624769;
      h[7] = -0.0114341;
      h[8] = 0.00167406;
      h[9] = 0.00306375;
      h[10] = -0.000448564;
      h[11] = -0.000820929;
      h[12] = 0.000120192;
      h[13] = 0.000219967;
      h[14] = -0.0000322054;
      h[15] = -0.00005894;
      *nh = 16;
      break;

    case 4:
      g[0] = 0.673072;
      g[1] = 0.331218;
      g[2] = -0.139359;
      g[3] = -0.12051;
      g[4] = 0.086389;
      g[5] = 0.0611801;
      g[6] = -0.0542989;
      g[7] = -0.034777;
      g[8] = 0.033388;
      g[9] = 0.0206275;
      g[10] = -0.0203475;
      g[11] = -0.0124183;
      g[12] = 0.0123625;
      g[13] = 0.00751369;
      g[14] = -0.00750374;
      g[15] = -0.00455348;
      g[16] = 0.00455363;
      g[17] = 0.00276047;
      g[18] = -0.00276406;
      g[19] = -0.00167279;
      g[20] = 0.00167938;
      *ng = 21;
      h[0] = 1.14324;
      h[1] = 0.643609;
      h[2] = -0.0937888;
      h[3] = -0.194993;
      h[4] = 0.030127;
      h[5] = 0.0699433;
      h[6] = -0.0108345;
      h[7] = -0.0252663;
      h[8] = 0.00391424;
      h[9] = 0.00912967;
      h[10] = -0.00141437;
      h[11] = -0.00329892;
      h[12] = 0.000511068;
      h[13] = 0.00119204;
      h[14] = -0.00018467;
      h[15] = -0.000430732;
      h[16] = 0.0000667289;
      h[17] = 0.000155641;
      h[18] = -0.0000241119;
      h[19] = -0.0000562395;
      *nh = 20;
      break;

    default:
      g[0] = 1.;
      *ng = 1;
      h[0] = 2.;
      *nh = 1;
      printf("Spline filters only defined for n=0,1,2,3,4\n");
      exit(1);
  }
}

/* ----------------------------------------------------------------------------
        Function:	PyramidFilterCenteredL2

        Purpose:	Initializes the symmetric down- and up-sampling filter arrays for
                                L2 spline pyramid of order 0 to 5 when the downsampled grid is centered.
                                These filters have then to be followed by a Haar filter.
                                        g: reduce filter
                                        h: expand filter

        Note:		filter arrays should be defined as
                                        float g[35],h[35]	filter arrays
                                        short *ng,*nh	number of taps
                                        short Order	order of the spline

        Author:		Patrick Brigger, NIH, BEIP,	April 1996
                                Daniel Sage, EPFL, Biomedical Imaging Group, November 1999

---------------------------------------------------------------------------- */
static void PyramidFilterCenteredL2(double g[], long *ng, double h[], long *nh, long Order)
{
  switch (Order) {
    case 0:
      g[0] = 1.;
      *ng = 1;
      h[0] = 2.;
      *nh = 1;
      break;

    case 1:
      g[0] = 0.820272;
      g[1] = 0.316987;
      g[2] = -0.203044;
      g[3] = -0.0849365;
      g[4] = 0.0544056;
      g[5] = 0.0227587;
      g[6] = -0.0145779;
      g[7] = -0.00609817;
      g[8] = 0.00390615;
      g[9] = 0.001634;
      g[10] = -0.00104665;
      g[11] = -0.000437829;
      g[12] = 0.000280449;
      g[13] = 0.000117316;
      g[14] = -0.000075146;
      g[15] = -0.0000314347;
      g[16] = 0.0000201353;
      *ng = 17;
      h[0] = 1.20096;
      h[1] = 0.473076;
      h[2] = -0.0932667;
      h[3] = 0.0249907;
      h[4] = -0.00669625;
      h[5] = 0.00179425;
      h[6] = -0.000480769;
      h[7] = 0.000128822;
      h[8] = -0.0000345177;
      *nh = 9;
      break;

    case 2:
      g[0] = 0.727973;
      g[1] = 0.314545;
      g[2] = -0.167695;
      g[3] = -0.0893693;
      g[4] = 0.0768426;
      g[5] = 0.0354175;
      g[6] = -0.0331015;
      g[7] = -0.0151496;
      g[8] = 0.0142588;
      g[9] = 0.00651781;
      g[10] = -0.00613959;
      g[11] = -0.00280621;
      g[12] = 0.00264356;
      g[13] = 0.00120827;
      g[14] = -0.00113825;
      g[15] = -0.000520253;
      g[16] = 0.000490105;
      g[17] = 0.000224007;
      g[18] = -0.000211028;
      g[19] = -0.0000964507;
      g[20] = 0.0000908666;
      *ng = 21;
      h[0] = 1.20711;
      h[1] = 0.585786;
      h[2] = -0.12132;
      h[3] = -0.100505;
      h[4] = 0.0208153;
      h[5] = 0.0172439;
      h[6] = -0.00357134;
      h[7] = -0.00295859;
      h[8] = 0.000612745;
      h[9] = 0.000507614;
      h[10] = -0.00010513;
      *nh = 11;
      break;

    case 3:
      g[0] = 0.70222;
      g[1] = 0.328033;
      g[2] = -0.159368;
      g[3] = -0.113142;
      g[4] = 0.0902447;
      g[5] = 0.0530861;
      g[6] = -0.0492084;
      g[7] = -0.0274987;
      g[8] = 0.0264529;
      g[9] = 0.0146073;
      g[10] = -0.0141736;
      g[11] = -0.0078052;
      g[12] = 0.00758856;
      g[13] = 0.00417626;
      g[14] = -0.00406225;
      g[15] = -0.00223523;
      g[16] = 0.00217454;
      g[17] = 0.00119638;
      g[18] = -0.00116412;
      g[19] = -0.000640258;
      g[20] = 0.000623379;
      *ng = 21;
      h[0] = 1.15089;
      h[1] = 0.623278;
      h[2] = -0.0961988;
      h[3] = -0.155743;
      h[4] = 0.0259827;
      h[5] = 0.041346;
      h[6] = -0.0067263;
      h[7] = -0.0112084;
      h[8] = 0.00187221;
      h[9] = 0.00296581;
      h[10] = -0.000481593;
      h[11] = -0.000805427;
      h[12] = 0.000134792;
      h[13] = 0.000212736;
      h[14] = -0.00003447;
      *nh = 15;
      break;

    case 4:
      g[0] = 0.672101;
      g[1] = 0.331667;
      g[2] = -0.138779;
      g[3] = -0.121385;
      g[4] = 0.0864024;
      g[5] = 0.0618776;
      g[6] = -0.0545165;
      g[7] = -0.0352403;
      g[8] = 0.0335951;
      g[9] = 0.0209537;
      g[10] = -0.0205211;
      g[11] = -0.0126439;
      g[12] = 0.0124959;
      g[13] = 0.0076682;
      g[14] = -0.00760135;
      g[15] = -0.00465835;
      g[16] = 0.00462238;
      g[17] = 0.00283148;
      g[18] = -0.00281055;
      g[19] = -0.00172137;
      g[20] = 0.00170884;
      *ng = 21;
      h[0] = 1.14324;
      h[1] = 0.643609;
      h[2] = -0.0937888;
      h[3] = -0.194993;
      h[4] = 0.030127;
      h[5] = 0.0699433;
      h[6] = -0.0108345;
      h[7] = -0.0252663;
      h[8] = 0.00391424;
      h[9] = 0.00912967;
      h[10] = -0.00141437;
      h[11] = -0.00329892;
      h[12] = 0.000511068;
      h[13] = 0.00119204;
      h[14] = -0.00018467;
      h[15] = -0.000430732;
      h[16] = 0.0000667289;
      h[17] = 0.000155641;
      h[18] = -0.0000241119;
      h[19] = -0.0000562396;
      *nh = 20;
      break;

    default:
      g[0] = 1.;
      *ng = 1;
      h[0] = 2.;
      *nh = 1;
      printf("Spline filters only defined for n=0,1,2,3,4\n");
      exit(1);
  }
}

/* ----------------------------------------------------------------------------
        Function:	PyramidFilterCenteredL2Derivate

        Purpose:	Initializes the symmetric down- and up-sampling filter arrays for
                                L2 DERIVATIVE spline pyramid of order 0 to 5 when the downsampled
                                grid is centered.
                                These filters have then to be followed by a Derivative Haar filter.
                                        g : reduce filter
                                        h : expand filter
        Note:		filter arrays should be defined as
                                        float g[35],h[35]	filter arrays
                                        short *ng,*nh	number of taps
                                        short Order	order of the spline

        Author:		Patrick Brigger, NIH, BEIP,	April 1996
                                Daniel Sage, EPFL, Biomedical Imaging Group, November 1999

---------------------------------------------------------------------------- */
// static  void PyramidFilterCenteredL2Derivate(double g[],long *ng,double h[],long *nh,long Order)
// {
// 	switch (Order) {
// 	case 0 :
// 		g[0] = 1.;
// 		*ng=1;
// 		h[0] = 2.;
// 		*nh=1;
// 		break;
//
// 	case 1 :
// 		g[0]  =  0.820272;
// 		g[1]  =  0.316987;
// 		g[2]  = -0.203044;
// 		g[3]  = -0.0849365;
//    		g[4]  =  0.0544056;
//    		g[5]  =  0.0227587;
//    		g[6]  = -0.0145779;
//    		g[7]  = -0.00609817;
//    		g[8]  =  0.00390615;
//    		g[9]  =  0.001634;
//    		g[10] = -0.00104665;
//    		g[11] = -0.000437829;
//    		g[12] =  0.000280449;
//    		g[13] =  0.000117316;
//    		g[14] = -0.000075146;
//    		g[15] = -0.0000314347;
//    		g[16] =  0.0000201353;
//    		*ng = 17;
//    		h[0]  =  1.20096;
//    		h[1]  =  1.20096;
//    		h[2]  = -0.254809;
//    		h[3]  =  0.068276;
//    		h[4]  = -0.0182945;
//    		h[5]  =  0.004902;
//    		h[6]  = -0.00131349;
//    		h[7]  =  0.000351947;
//    		h[8]  = -0.000094304;
//    		h[9]  =  0.0000252687;
//    		*nh = 10;
// 		break;
//
// 	case 2 :
// 		g[0]  =  0.727973;
// 		g[1]  =  0.314545;
// 		g[2]  = -0.167695;
//    		g[3]  = -0.0893693;
//    		g[4]  =  0.0768426;
//    		g[5]  =  0.0354175;
//    		g[6]  = -0.0331015;
//    		g[7]  = -0.0151496;
//    		g[8]  =  0.0142588;
//    		g[9]  =  0.00651781;
//    		g[10] = -0.00613959;
//    		g[11] = -0.00280621;
//    		g[12] =  0.00264356;
//    		g[13] =  0.00120827;
//    		g[14] = -0.00113825;
//    		g[15] = -0.000520253;
//    		g[16] =  0.000490105;
//    		g[17] =  0.000224007;
//    		g[18] = -0.000211028;
//    		g[19] = -0.0000964507;
//   	    g[20] =  0.0000908666;
//   	    *ng = 21;
//   	    h[0]  =  1.20711;
//   	    h[1]  =  0.585786;
//   	    h[2]  = -0.12132;
//   	    h[3]  = -0.100505;
//    		h[4]  =  0.0208153;
//    		h[5]  =  0.0172439;
//    		h[6]  = -0.00357134;
//    		h[7]  = -0.00295859;
//    		h[8]  =  0.000612745;
//    		h[9]  =  0.000507614;
//    		h[10] = -0.00010513;
//    		*nh = 11;
//    		break;
//
// 	case 3 :
// 	    g[0]  =  0.70222;
// 	    g[1]  =  0.328033;
// 	    g[2]  = -0.159368;
// 	    g[3]  = -0.113142;
//    		g[4]  =  0.0902447;
//    		g[5]  =  0.0530861;
//    		g[6]  = -0.0492084;
//    		g[7]  = -0.0274987;
//    		g[8]  =  0.0264529;
//    		g[9]  =  0.0146073;
//    		g[10] = -0.0141736;
//    		g[11] = -0.0078052;
//    		g[12] =  0.00758856;
//    		g[13] =  0.00417626;
//    		g[14] = -0.00406225;
//    		g[15] = -0.00223523;
//    		g[16] =  0.00217454;
//    		g[17] =  0.00119638;
//    		g[18] = -0.00116412;
//    		g[19] = -0.000640258;
//    		g[20] =  0.000623379;
//    		*ng = 21;
// 		h[0]  =  1.15089;
// 		h[1]  =  0.623278;
// 		h[2]  = -0.0961988;
// 		h[3]  = -0.155743;
// 		h[4]  =  0.0259827;
// 		h[5]  =  0.041346;
// 		h[6]  = -0.0067263;
// 		h[7]  = -0.0112084;
// 		h[8]  =  0.00187221;
// 		h[9]  =  0.00296581;
// 		h[10] = -0.000481593;
// 		h[11] = -0.000805427;
// 		h[12] =  0.000134792;
// 		h[13] =  0.000212736;
// 		h[14] = -0.00003447;
//   		*nh = 15;
//    		break;
//
// 	case 4:
// 		g[0]  =  0.672101;
// 		g[1]  =  0.331667;
// 		g[2]  = -0.138779;
//    		g[3]  = -0.121385;
//    		g[4]  =  0.0864024;
//    		g[5]  =  0.0618776;
//    		g[6]  = -0.0545165;
//    		g[7]  = -0.0352403;
//    		g[8]  =  0.0335951;
//    		g[9]  =  0.0209537;
//    		g[10] = -0.0205211;
//    		g[11] = -0.0126439;
//    		g[12] =  0.0124959;
//    		g[13] =  0.0076682;
//    		g[14] = -0.00760135;
//    		g[15] = -0.00465835;
//    		g[16] =  0.00462238;
//    		g[17] =  0.00283148;
//    		g[18] = -0.00281055;
//    		g[19] = -0.00172137;
//    		g[20] =  0.00170884;
//    		*ng = 21;
//    		h[0]  =  1.14324;
//    		h[1]  =  0.643609;
//    		h[2]  = -0.0937888;
//    		h[3]  = -0.194993;
//    		h[4]  =  0.030127;
//    		h[5]  =  0.0699433;
//    		h[6]  = -0.0108345;
//    		h[7]  = -0.0252663;
//    		h[8]  =  0.00391424;
//    		h[9]  =  0.00912967;
//    		h[10] = -0.00141437;
//    		h[11] = -0.00329892;
//    		h[12] =  0.000511068;
//    		h[13] =  0.00119204;
//    		h[14] = -0.00018467;
//    		h[15] = -0.000430732;
//    		h[16] =  0.0000667289;
//    		h[17] =  0.000155641;
//    		h[18] = -0.0000241119;
//    		h[19] = -0.0000562396;
//    		*nh = 20;
//    		break;
//
// 	default :
// 		g[0]  = 1.;
// 		*ng=1;
// 		h[0]  = 2.;
// 		*nh=1;
// 		printf( "Spline filters only defined for n=0,1,2,3,4\n");
// 		exit(1);
// 	}
// }

/* ----------------------------------------------------------------------------

        Function:
                GetPyramidFilter

        Purpose:
                Get the coefficients of the filter (reduce and expand filter)
                Return the coefficients in g[ng] and in h[nh]

        Convention:
                g[ng] for the reduce filter
                h[nh] for the expansion filter

        Parameters:
                Filter is the name of the filter

                Order is the order for the filters based on splines
                        For the "Spline" filter, Order is 0, 1, 2 or 3
                        For the "Spline L2" filter, Order is 0, 1, 3 or 5
                        For the "Centered Spline" filter, Order is 0, 1, 2, 3 or 4
                        For the "Centered Spline L2" filter, Order is 0, 1, 2, 3 or 4

                IsCentered is a return value indicates if the filter is a centered filter
                        TRUE if it is a centered filter
                        FALSE if it is not a centered filter

---------------------------------------------------------------------------- */
static int GetPyramidFilter(const char *Filter, long Order, double g[], long *ng, double h[], long *nh, short *IsCentered)
{
  ng[0] = -1L;
  nh[0] = -1L;
  *IsCentered = FALSE;

  if (!strcmp(Filter, "Spline")) {
    PyramidFilterSplinel2(g, ng, h, nh, Order);
    *IsCentered = FALSE;
  }

  if (!strcmp(Filter, "Spline L2")) {
    PyramidFilterSplineL2(g, ng, h, nh, Order);
    *IsCentered = FALSE;
  }

  if (!strcmp(Filter, "Centered Spline")) {
    PyramidFilterCentered(g, ng, h, nh, Order);
    *IsCentered = TRUE;
  }

  if (!strcmp(Filter, "Centered Spline L2")) {
    PyramidFilterCenteredL2(g, ng, h, nh, Order);
    *IsCentered = TRUE;
  }

  if (ng[0] == -1L && nh[0] == -1L) {
    printf("This familly filters is unknown!\n");
    return (0);
  }
  return (1);
}

/* ----------------------------------------------------------------------------

        Function:
                ReduceStandard_1D

        Purpose:
                Basic function to reduce a 1D signal

        Parameters:
                In[NxIn] is the input signal	(NxIn should be greater than 2 and even)
                Out[NxIn/2] is the output signal
                g[ng] is an array that contains the coefficients of the filter

        Author:
                Michael Unser, NIH, BEIP, June 1994
                Daniel Sage, EPFL, Biomedical Imaging Group, April 1999

---------------------------------------------------------------------------- */
static void ReduceStandard_1D(double In[], long NxIn, double Out[], double g[], long ng)
{
  long k, i, i1, i2;
  long kk, kn, nred, n;

  nred = NxIn / 2L;
  n = nred * 2L;
  kn = n - 1L; /* kn=n-2; DS Modified */

  if (ng < 2L) { /* length filter < 2 */
    for (kk = 0L; kk < nred; kk++) {
      k = 2L * kk;
      i2 = k + 1L;
      if (i2 > n - 1L) i2 = kn - i2;
      Out[kk] = (In[k] + In[i2]) / 2.;
    }
  }

  else {
    for (kk = 0L; kk < nred; kk++) {
      k = 2L * kk;
      Out[kk] = In[k] * g[0];
      for (i = 1L; i < ng; i++) {
        i1 = k - i;
        i2 = k + i;
        if (i1 < 0L) {
          i1 = (-i1) % kn;
          if (i1 > n - 1) i1 = kn - i1;
        }
        if (i2 > n - 1L) {
          i2 = i2 % kn;
          if (i2 > n - 1L) i2 = kn - i2;
        }
        Out[kk] = Out[kk] + g[i] * (In[i1] + In[i2]);
      }
    }
  }
}

/* ----------------------------------------------------------------------------

        Function:
                ExpandStandard_1D

        Purpose:
                Basic function to expand a 1D signal

        Parameters:
                In[NxIn] is the input signal	(NxIn should be greater than 1)
                Out[NxIn*2] is the output signal
                w[nw] is an array that contains the coefficients of the filter

        Author:
                Michael Unser, NIH, BEIP, June 1994
                Daniel Sage, EPFL, Biomedical Imaging Group, April 1999

---------------------------------------------------------------------------- */
static void ExpandStandard_1D(double In[], long NxIn, double Out[], double h[], long nh)
{
  long k, j, i, i1, i2;
  long kn, nexp, n;

  nexp = NxIn * 2L;
  n = NxIn;
  kn = n - 1;

  if (nh < 2L) {
    for (i = 0L; i < NxIn; i++) {
      j = i * 2L;
      Out[j] = In[i];
      Out[j + 1] = In[i];
    }
  }

  else {
    for (i = 0L; i < nexp; i++) {
      Out[i] = 0.0;
      for (k = (i % 2L); k < nh; k += 2L) {
        i1 = (i - k) / 2L;
        if (i1 < 0L) {
          i1 = (-i1) % kn;
          if (i1 > kn) i1 = kn - i1;
        }
        Out[i] = Out[i] + h[k] * In[i1];
      }

      for (k = 2L - (i % 2L); k < nh; k += 2L) {
        i2 = (i + k) / 2L;
        if (i2 > kn) {
          i2 = i2 % kn;
          i2 = kn - i2;
          if (i2 > kn) i2 = kn - i2;
        }
        Out[i] = Out[i] + h[k] * In[i2];
      }
    }
  }
}

/* ----------------------------------------------------------------------------

        Function:	ReduceCentered_1D

        Purpose:	Reduces an image by a factor of two
                                The reduced image grid is between the finer grid

        Parameters:
                In[NxIn] is the input signal	(NxIn should be greater than 2 and even)
                Out[NxIn/2] is the output signal
                g[ng] is an array that contains the coefficients of the filter

        Author:
                Michael Unser, NIH, BEIP, June 1994
                Patrick Brigger, NIH, BEIP,	May 1996, modified
                Daniel Sage, EPFL, Biomedical Imaging Group, April 1999

---------------------------------------------------------------------------- */
static void ReduceCentered_1D(double In[], long NxIn, double Out[], double g[], long ng)
{
  double *y_tmp;
  long k, i, i1, i2;
  long kk, kn, nred, n;

  nred = NxIn / 2L;
  n = nred * 2L;
  kn = 2L * n;

  /* --- Allocated memory for a temporary buffer --- */
  y_tmp = (double *)malloc((size_t)(n * (long)sizeof(double)));
  if (y_tmp == (double *)NULL) {
    printf("Out of memory in reduce_centered!\n");
    exit(1);
  }

  /* --- Apply the symmetric filter to all the coefficients --- */
  for (k = 0L; k < n; k++) {
    y_tmp[k] = In[k] * g[0];
    for (i = 1L; i < ng; i++) {
      i1 = k - i;
      i2 = k + i;
      if (i1 < 0L) {
        i1 = (2L * n - 1 - i1) % kn;
        if (i1 >= n) i1 = kn - i1 - 1L;
      }
      if (i2 > (n - 1L)) {
        i2 = i2 % kn;
        if (i2 >= n) i2 = kn - i2 - 1L;
      }
      y_tmp[k] = y_tmp[k] + g[i] * (In[i1] + In[i2]);
    }
  }

  /* --- Now apply the Haar and perform downsampling --- */
  for (kk = 0L; kk < nred; kk++) {
    k = 2L * kk;
    Out[kk] = (y_tmp[k] + y_tmp[k + 1]) / 2.;
  }

  /* --- Free allocated memory --- */
  free(y_tmp);
}

/* ----------------------------------------------------------------------------

        Function:
                ExpandCentered_1D

        Purpose:
                Expands an 1D signal by a factor of two using the specified
                FIR interpolation filter. The expansion is based on
                a subsampled grid that is centered with respect to the
                finer grid

        Parameters:
                In[NxIn] is the input signal	(NxIn should be greater than 1)
                Out[NxIn"*2] is the output signal
                h[nh] is an array that contains the coefficients of the filter

        Author:
                Michael Unser, NIH, BEIP, June 1994
                Patrick Brigger, NIH, BEIP,	May 1996, modified
                Daniel Sage, EPFL, Biomedical Imaging Group, April 1999

---------------------------------------------------------------------------- */
static void ExpandCentered_1D(double In[], long NxIn, double Out[], double h[], long nh)
{
  long k, i, j, i1, k0, i2;
  long kk, kn, nexp, n;

  nexp = NxIn * 2L;
  n = NxIn;
  kn = 2L * n;
  k0 = (nh / 2L) * 2 - 1L;

  for (i = 0L; i < NxIn; i++) {
    j = i * 2L;
    Out[j] = In[i] * h[0];
    for (k = 2L; k < nh; k = k + 2L) {
      i1 = i - k / 2L;
      i2 = i + k / 2L;
      if (i1 < 0L) {                 /* Provide the correct border conditions. */
        i1 = (2L * n - 1 - i1) % kn; /* --> pseudo mirror image                */
        if (i1 >= n) i1 = kn - i1 - 1L;
      }
      if (i2 >= n) {
        i2 = (i2) % kn;
        if (i2 >= n) i2 = kn - i2 - 1L;
      }
      Out[j] = Out[j] + h[k] * (In[i1] + In[i2]);
    }
    Out[j + 1] = 0.;
    for (k = -k0; k < nh; k = k + 2L) {
      kk = abs(k); /* filter coeff. are shifted with respect to above. */
      i1 = i + (k + 1L) / 2L;
      if (i1 < 0L) {
        i1 = (2 * n - 1 - i1) % kn;
        if (i1 > n - 1) i1 = kn - i1 - 1L;
      }
      if (i1 >= n) {
        i1 = (i1) % kn;
        if (i1 >= n) i1 = kn - i1 - 1;
      }
      Out[j + 1L] = Out[j + 1L] + h[kk] * In[i1];
    }
  }

  /* Now apply the Haar[-x] and  */
  for (j = nexp - 1L; j > 0L; j--) Out[j] = (Out[j] + Out[j - 1]) / 2.0;
  Out[0] /= 2.0;
}

/* ----------------------------------------------------------------------------

        Function:
                Reduce_1D

        Purpose:
                Router function to call ReduceStandard_1D or ReduceCentered_1D

---------------------------------------------------------------------------- */
static void Reduce_1D(double In[], long NxIn, double Out[], double g[], long ng, short IsCentered)
{
  if (IsCentered)
    ReduceCentered_1D(In, NxIn, Out, g, ng);
  else
    ReduceStandard_1D(In, NxIn, Out, g, ng);
}
/* ----------------------------------------------------------------------------

        Function:
                Expand_1D

        Purpose:
                Router function to call ExpandStandard_1D or ExpandCentered_1D

---------------------------------------------------------------------------- */
static void Expand_1D(double In[], long NxIn, double Out[], double h[], long nh, short IsCentered)
{
  if (IsCentered)
    ExpandCentered_1D(In, NxIn, Out, h, nh);
  else
    ExpandStandard_1D(In, NxIn, Out, h, nh);
}

/*--------------------------------------------------------------------------*/

extern MRI_BSPLINE *MRIallocBSpline(int width, int height, int depth, int nframes)
{
  MRI_BSPLINE *bspline = (MRI_BSPLINE *)calloc(1, sizeof(MRI_BSPLINE));
  if (!bspline) ErrorExit(ERROR_NO_MEMORY, "MRIalloc: could not allocate MRI_BSPLINE\n");
  bspline->coeff = MRIallocSequence(width, height, depth, MRI_FLOAT, nframes);
  if (!bspline->coeff) ErrorExit(ERROR_NO_MEMORY, "MRIalloc: could not allocate MRI_BSPLINE\n");
  bspline->degree = -1;
  bspline->srctype = -1;
  return bspline;
}

extern int MRIfreeBSpline(MRI_BSPLINE **pbspline)
{
  MRI_BSPLINE *bspline = *pbspline;
  if (!bspline) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIfreeBSpline: null pointer\n"));
  MRIfree(&bspline->coeff);
  free(bspline);
  *pbspline = NULL;

  return (NO_ERROR);
}

extern MRI_BSPLINE *MRItoBSpline(const MRI *mri_src, MRI_BSPLINE *bspline, int degree)
{
  printf("MRItoBSpline degree %i\n", degree);

  double *Lines[_MAX_FS_THREADS];
  double Pole[4];
  int NbPoles;
  int x, y, z, f;
  int Width = mri_src->width;
  int Height = mri_src->height;
  int Depth = mri_src->depth;
  int Frames = mri_src->nframes;
  int nthreads = 1, tid = 0, i;

  //  double eps = 1e-5;
  //   if (fabs(mri_src->xsize - mri_src->ysize) > eps || fabs(mri_src->xsize = mri_src->zsize) > eps)
  //   {
  //     printf("ERROR MRItoBSpline: Input MRI not isotropic!\n");
  //     exit(1);
  //   }

  if (!bspline) {
    bspline = MRIallocBSpline(Width, Height, Depth, Frames);
    MRIcopyHeader(mri_src, bspline->coeff);
    bspline->coeff->type = MRI_FLOAT;
    bspline->coeff->outside_val = mri_src->outside_val;
  }

  if (!bspline->coeff) {
    bspline->coeff = MRIallocSequence(Width, Height, Depth, MRI_FLOAT, Frames);
    MRIcopyHeader(mri_src, bspline->coeff);
    bspline->coeff->type = MRI_FLOAT;
    bspline->coeff->outside_val = mri_src->outside_val;
  }

  bspline->degree = degree;
  bspline->srctype = mri_src->type;
  bspline->srcneg = FALSE;

  if (bspline->coeff->type != MRI_FLOAT) {
    printf("ERROR MRItoBSpline: BSpline Coeff MRI needs to be float!\n");
    exit(1);
  }

  if (bspline->coeff->width != Width || bspline->coeff->depth != Depth || bspline->coeff->height != Height ||
      bspline->coeff->nframes != Frames) {
    printf("ERROR MRItoBSpline: BSPline Coeff MRI dimensions must be same as Source MRI!\n");
    exit(1);
  }

  NbPoles = initPoles(Pole, degree);

/* convert the image samples into interpolation coefficients */
#ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif
  // printf("nthreads : %d\n",nthreads);
  for (i = 0; i < nthreads; i++) {
    Lines[i] = (double *)malloc((size_t)(Width * sizeof(double)));
    if (Lines[i] == (double *)NULL) {
      printf("ERROR MRItoBSpline: X Line allocation failed\n");
      exit(1);
    }
  }
  y = 0;
  f = 0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(tid, y, f, Lines) shared(Depth, Height, Width, Frames, Pole, NbPoles, bspline) \
    schedule(static, 1)
#endif
  for (z = 0; z < Depth; z++) {
    ROMP_PFLB_begin
    
#ifdef HAVE_OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    for (y = 0; y < Height; y++) {
      for (f = 0; f < Frames; f++) {
        // printf("f: %i  z: %i  y: %i\n",f,z,y);
        getXLine(mri_src, y, z, f, Lines[tid]);

        // check if we have a negative value
        // needs to be done only once and only if we did not find one already
        if (!bspline->srcneg) {
          for (x = 0; x < Width; x++) {
            if (Lines[tid][x] < 0.0) {
              bspline->srcneg = TRUE;  // write to shared, racing OK (only set to true)
              break;
            }
          }
        }

        if (Width > 1) {
          ConvertToInterpolationCoefficients(Lines[tid], Width, Pole, NbPoles, DBL_EPSILON);
        }
        setXLine(bspline->coeff, y, z, f, Lines[tid]);
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  for (i = 0; i < nthreads; i++) {
    free(Lines[i]);
  }

  if (Height > 1) {
    for (i = 0; i < nthreads; i++) {
      Lines[i] = (double *)malloc((size_t)(Height * sizeof(double)));
      if (Lines[i] == (double *)NULL) {
        printf("ERROR MRItoBSpline: Y Line allocation failed\n");
        exit(1);
      }
    }
    x = 0;
    f = 0;
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) firstprivate(tid, x, f, Lines) shared(Depth, Height, Width, Frames, Pole, NbPoles, bspline) \
    schedule(static, 1)
#endif
    for (z = 0; z < Depth; z++) {
      ROMP_PFLB_begin
      
#ifdef HAVE_OPENMP
      tid = omp_get_thread_num();
#else
      tid = 0;
#endif
      for (x = 0; x < Width; x++) {
        for (f = 0; f < Frames; f++) {
          // printf("f: %i  z: %i  x: %i\n",f,z,x);
          getYLine(bspline->coeff, x, z, f, Lines[tid]);
          ConvertToInterpolationCoefficients(Lines[tid], Height, Pole, NbPoles, DBL_EPSILON);
          setYLine(bspline->coeff, x, z, f, Lines[tid]);
        }
      }
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    for (i = 0; i < nthreads; i++) {
      free(Lines[i]);
    }
  }

  if (Depth > 1) {
    for (i = 0; i < nthreads; i++) {
      Lines[i] = (double *)malloc((size_t)(Depth * sizeof(double)));
      if (Lines[i] == (double *)NULL) {
        printf("ERROR MRItoBSpline: Z Line allocation failed\n");
        exit(1);
      }
    }
    x = 0;
    f = 0;
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) firstprivate(tid, x, f, Lines) shared(Depth, Height, Width, Frames, Pole, NbPoles, bspline) \
    schedule(static, 1)
#endif
    for (y = 0; y < Height; y++) {
      ROMP_PFLB_begin
      
#ifdef HAVE_OPENMP
      tid = omp_get_thread_num();
#else
      tid = 0;
#endif
      for (x = 0; x < Width; x++) {
        for (f = 0; f < Frames; f++) {
          // printf("f: %i  y: %i  x: %i\n",f,y,x);
          getZLine(bspline->coeff, x, y, f, Lines[tid]);
          ConvertToInterpolationCoefficients(Lines[tid], Depth, Pole, NbPoles, DBL_EPSILON);
          setZLine(bspline->coeff, x, y, f, Lines[tid]);
        }
      }
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    for (i = 0; i < nthreads; i++) {
      free(Lines[i]);
    }
  }

  return bspline;
}

extern int MRIsampleBSpline(const MRI_BSPLINE *bspline, /* input B-spline array of coefficients */
                            double x,                   /* x coordinate where to interpolate */
                            double y,                   /* y coordinate where to interpolate */
                            double z,                   /* y coordinate where to interpolate */
                            const int frame,
                            double *pval)

{
  int Width = bspline->coeff->width;
  int Height = bspline->coeff->height;
  int Depth = bspline->coeff->depth;
  int SplineDegree = bspline->degree;

  int OutOfBounds;
  double w, w2;
  double xWeight[10], yWeight[10], zWeight[10];
  double interpolated;
  int xIndex[10], yIndex[10], zIndex[10];
  int Width2 = 2 * Width - 2, Height2 = 2 * Height - 2, Depth2 = 2 * Depth - 2;
  long i, j, k, l;

  if (frame >= bspline->coeff->nframes || frame < 0) {
    *pval = bspline->coeff->outside_val;
    return (NO_ERROR);
  }

  //   if (frame >= 1)
  //   {
  //     ErrorReturn(ERROR_UNSUPPORTED,
  //                 (ERROR_UNSUPPORTED,
  //                  "MRIsampleBSpline: unsupported frame %d", frame)) ;
  //
  //   }

  OutOfBounds = MRIindexNotInVolume(bspline->coeff, x, y, z);
  if (OutOfBounds == 1) {
    /* unambiguoulsy out of bounds */
    *pval = bspline->coeff->outside_val;
    return (NO_ERROR);
  }

  /* compute the interpolation indexes */
  if (SplineDegree & 1) {
    i = floor(x) - SplineDegree / 2;
    j = floor(y) - SplineDegree / 2;
    k = floor(z) - SplineDegree / 2;
    for (l = 0; l <= SplineDegree; l++) {
      xIndex[l] = i++;
      yIndex[l] = j++;
      zIndex[l] = k++;
    }
  }
  else {
    i = (long)floor(x + 0.5) - SplineDegree / 2;
    j = (long)floor(y + 0.5) - SplineDegree / 2;
    k = (long)floor(z + 0.5) - SplineDegree / 2;
    for (l = 0; l <= SplineDegree; l++) {
      xIndex[l] = i++;
      yIndex[l] = j++;
      zIndex[l] = k++;
    }
  }

  /* compute the interpolation weights */
  if (Width == 1)
    xWeight[0] = 1.0;
  else
    BsplineWeights(SplineDegree, &x, xIndex, xWeight);
  if (Height == 1)
    yWeight[0] = 1.0;
  else
    BsplineWeights(SplineDegree, &y, yIndex, yWeight);
  if (Depth == 1)
    zWeight[0] = 1.0;
  else
    BsplineWeights(SplineDegree, &z, zIndex, zWeight);

  /* apply the mirror boundary conditions */
  for (k = 0; k <= SplineDegree; k++) {
    xIndex[k] = (Width == 1) ? (0)
                             : ((xIndex[k] < 0) ? (-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
                                                : (xIndex[k] - Width2 * (xIndex[k] / Width2)));
    if (Width <= xIndex[k]) {
      xIndex[k] = Width2 - xIndex[k];
    }
    yIndex[k] = (Height == 1) ? (0)
                              : ((yIndex[k] < 0) ? (-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
                                                 : (yIndex[k] - Height2 * (yIndex[k] / Height2)));
    if (Height <= yIndex[k]) {
      yIndex[k] = Height2 - yIndex[k];
    }
    zIndex[k] = (Depth == 1) ? (0)
                             : ((zIndex[k] < 0) ? (-zIndex[k] - Depth2 * ((-zIndex[k]) / Depth2))
                                                : (zIndex[k] - Depth2 * (zIndex[k] / Depth2)));
    if (Depth <= zIndex[k]) {
      zIndex[k] = Depth2 - zIndex[k];
    }
  }

  /* perform interpolation */
  interpolated = 0.0;
  int sdk = SplineDegree;
  int sdj = SplineDegree;
  int sdi = SplineDegree;
  if (Width == 1) sdi = 0;
  if (Height == 1) sdj = 0;
  if (Depth == 1) sdk = 0;
  for (k = 0; k <= sdk; k++) {
    w2 = 0.0;
    for (j = 0; j <= sdj; j++) {
      w = 0.0;
      for (i = 0; i <= sdi; i++) {
        w += xWeight[i] * MRIgetVoxVal(bspline->coeff, xIndex[i], yIndex[j], zIndex[k], frame);
      }
      w2 += yWeight[j] * w;
    }
    interpolated += zWeight[k] * w2;
  }

  // if input does not have negative values, clip to zero:
  if (!bspline->srcneg && interpolated < 0.0) interpolated = 0.0;

  /*  // clip interpolation to limits:
    // removed (let the caller decide)
    switch (bspline->srctype)
    {
    case MRI_UCHAR:
      if (interpolated < UCHAR_MIN) interpolated = UCHAR_MIN;
      if (interpolated > UCHAR_MAX) interpolated = UCHAR_MAX;
      break;
    case MRI_SHORT:
      if (interpolated < SHORT_MIN) interpolated = SHORT_MIN;
      if (interpolated > SHORT_MAX) interpolated = SHORT_MAX;
      break;
    case MRI_INT:
      if (interpolated < INT_MIN) interpolated = INT_MIN;
      if (interpolated > INT_MAX) interpolated = INT_MAX;
      break;
    case MRI_LONG:
      if (interpolated < LONG_MIN) interpolated = LONG_MIN;
      if (interpolated > LONG_MAX) interpolated = LONG_MAX;
      break;
    } */

  *pval = interpolated;
  return NO_ERROR;
}

extern int MRIsampleSeqBSpline(
    const MRI_BSPLINE *bspline, double x, double y, double z, float *valvect, int firstframe, int lastframe)
{
  int f;
  int OutOfBounds = MRIindexNotInVolume(bspline->coeff, x, y, z);
  if (OutOfBounds == 1) {
    /* unambiguously out of bounds */
    for (f = firstframe; f <= lastframe; f++) valvect[f] = bspline->coeff->outside_val;
    return (NO_ERROR);
  }

  int Width = bspline->coeff->width;
  int Height = bspline->coeff->height;
  int Depth = bspline->coeff->depth;
  int SplineDegree = bspline->degree;

  double w, w2;
  double xWeight[10], yWeight[10], zWeight[10];
  double interpolated;
  int xIndex[10], yIndex[10], zIndex[10];
  int Width2 = 2 * Width - 2, Height2 = 2 * Height - 2, Depth2 = 2 * Depth - 2;
  long i, j, k, l;

  /* compute the interpolation indexes */
  if (SplineDegree & 1) {
    i = floor(x) - SplineDegree / 2;
    j = floor(y) - SplineDegree / 2;
    k = floor(z) - SplineDegree / 2;
    for (l = 0; l <= SplineDegree; l++) {
      xIndex[l] = i++;
      yIndex[l] = j++;
      zIndex[l] = k++;
    }
  }
  else {
    i = (long)floor(x + 0.5) - SplineDegree / 2;
    j = (long)floor(y + 0.5) - SplineDegree / 2;
    k = (long)floor(z + 0.5) - SplineDegree / 2;
    for (l = 0; l <= SplineDegree; l++) {
      xIndex[l] = i++;
      yIndex[l] = j++;
      zIndex[l] = k++;
    }
  }

  /* compute the interpolation weights */
  if (Width == 1)
    xWeight[0] = 1.0;
  else
    BsplineWeights(SplineDegree, &x, xIndex, xWeight);
  if (Height == 1)
    yWeight[0] = 1.0;
  else
    BsplineWeights(SplineDegree, &y, yIndex, yWeight);
  if (Depth == 1)
    zWeight[0] = 1.0;
  else
    BsplineWeights(SplineDegree, &z, zIndex, zWeight);

  /* apply the mirror boundary conditions */
  for (k = 0; k <= SplineDegree; k++) {
    xIndex[k] = (Width == 1) ? (0)
                             : ((xIndex[k] < 0) ? (-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
                                                : (xIndex[k] - Width2 * (xIndex[k] / Width2)));
    if (Width <= xIndex[k]) {
      xIndex[k] = Width2 - xIndex[k];
    }
    yIndex[k] = (Height == 1) ? (0)
                              : ((yIndex[k] < 0) ? (-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
                                                 : (yIndex[k] - Height2 * (yIndex[k] / Height2)));
    if (Height <= yIndex[k]) {
      yIndex[k] = Height2 - yIndex[k];
    }
    zIndex[k] = (Depth == 1) ? (0)
                             : ((zIndex[k] < 0) ? (-zIndex[k] - Depth2 * ((-zIndex[k]) / Depth2))
                                                : (zIndex[k] - Depth2 * (zIndex[k] / Depth2)));
    if (Depth <= zIndex[k]) {
      zIndex[k] = Depth2 - zIndex[k];
    }
  }

  /* perform interpolation */
  int sdk = SplineDegree;
  int sdj = SplineDegree;
  int sdi = SplineDegree;
  if (Width == 1) sdi = 0;
  if (Height == 1) sdj = 0;
  if (Depth == 1) sdk = 0;

  for (f = firstframe; f <= lastframe; f++) {
    interpolated = 0.0;
    for (k = 0; k <= sdk; k++) {
      w2 = 0.0;
      for (j = 0; j <= sdj; j++) {
        w = 0.0;
        for (i = 0; i <= sdi; i++) {
          w += xWeight[i] * MRIgetVoxVal(bspline->coeff, xIndex[i], yIndex[j], zIndex[k], f);
        }
        w2 += yWeight[j] * w;
      }
      interpolated += zWeight[k] * w2;
    }

    // if input does not have negative values, clip to zero:
    if (!bspline->srcneg && interpolated < 0.0) interpolated = 0.0;

    valvect[f] = interpolated;
  }

  return NO_ERROR;
}

MRI *MRIlinearTransformBSpline(const MRI_BSPLINE *bspline, MRI *mri_dst, MATRIX *mA)
// mainly copied from mri.c
{
  int y1, y2, y3, width, height, depth, frame;
  // VECTOR *v_X, *v_Y ;   /* original and transformed coordinate systems */
  VECTOR *v_X[_MAX_FS_THREADS], *v_Y[_MAX_FS_THREADS]; /* original and transformed coordinate systems */
  MATRIX *mAinv;                                       /* inverse of mA */
  double val, x1, x2, x3;
  int nthreads = 1, tid = 0;
#ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();
#endif

  mAinv = MatrixInverse(mA, NULL); /* will sample from dst back to src */
  if (!mAinv) ErrorReturn(NULL, (ERROR_BADPARM, "MRIlinearTransformSpline: xform is singular"));

  if (!mri_dst) {
    // mri_dst = MRIclone(mri_src, NULL) ;
    mri_dst = MRIcloneDifferentType(bspline->coeff, bspline->srctype);
  }
  else
    MRIclear(mri_dst);

  width = mri_dst->width;
  height = mri_dst->height;
  depth = mri_dst->depth;

  for (tid = 0; tid < nthreads; tid++) {
    v_X[tid] = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
    v_Y[tid] = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
    v_Y[tid]->rptr[4][1] = 1.0f;
  }

  y2 = 0;
  y1 = 0;
  x1 = 0;
  x2 = 0;
  x3 = 0;
  frame = 0;
  val = 0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(tid, y2, y1, x1, x2, x3, frame, val, v_X, v_Y) \
    shared(depth, height, width, mAinv, bspline, mri_dst) schedule(static, 1)
#endif
  for (y3 = 0; y3 < depth; y3++) {
    ROMP_PFLB_begin
    
#ifdef HAVE_OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    V3_Z(v_Y[tid]) = y3;
    for (y2 = 0; y2 < height; y2++) {
      V3_Y(v_Y[tid]) = y2;
      for (y1 = 0; y1 < width; y1++) {
        V3_X(v_Y[tid]) = y1;
        MatrixMultiply(mAinv, v_Y[tid], v_X[tid]);

        x1 = V3_X(v_X[tid]);
        x2 = V3_Y(v_X[tid]);
        x3 = V3_Z(v_X[tid]);

        for (frame = 0; frame < bspline->coeff->nframes; frame++) {
          MRIsampleBSpline(bspline, x1, x2, x3, frame, &val);

          // will clip the val according to mri_dst type:
          MRIsetVoxVal(mri_dst, y1, y2, y3, frame, val);
        }
      }
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  for (tid = 0; tid < nthreads; tid++) {
    MatrixFree(&v_X[tid]);
    MatrixFree(&v_Y[tid]);
  }

  MatrixFree(&mAinv);

  mri_dst->ras_good_flag = 1;

  return (mri_dst);
}

MRI *MRIapplyRASlinearTransformBSpline(const MRI_BSPLINE *bspline, MRI *mri_dst, MATRIX *m_ras_xform)
{
  MATRIX *m_voxel_xform;

  m_voxel_xform = MRIrasXformToVoxelXform(bspline->coeff, mri_dst, m_ras_xform, NULL);
  // fprintf(stderr, "applying the vox to vox linear transform\n");
  // MatrixPrint(stderr, m_voxel_xform);
  mri_dst = MRIlinearTransformBSpline(bspline, mri_dst, m_voxel_xform);
  MatrixFree(&m_voxel_xform);
  return (mri_dst);
}

MRI *LTAtransformBSpline(const MRI_BSPLINE *bspline, MRI *mri_dst, LTA *lta)
// mainly copied from transform.c
{
  // int         y1, y2, y3, width, height, depth, xi, yi, zi, f ;
  // VECTOR      *v_X, *v_Y ;/* original and transformed coordinate systems */
  // double        x1, x2, x3 ;
  // MATRIX      *m_L, *m_L_inv ;
  LT *tran = &lta->xforms[0];
  MATRIX *r2i = 0;
  MATRIX *i2r = 0;
  MATRIX *tmp = 0;
  // MATRIX *v2v = 0;
  MRI *resMRI = 0;

  if (lta->num_xforms == 1) {
    /////////////////////////////////////////////////////////////////////////
    //  The transform was created using src and dst
    //           src vox ----> RAS
    //            |             |
    //            V             V
    //           dst vox ----> RAS
    //
    //  Note: RAS is the physical space.  vox is the way to
    // embed voxel in physical space
    //
    //  You can take arbitray embedding of the dst but keeping the RAS same.
    //
    //           src vox ----> RAS
    //            | v2v         | r2r
    //            V             V
    //           dst vox ----> RAS
    //            |             || identity
    //            V             ||
    //           dst' vox ---> RAS
    //
    //  where dst->dst' map is given by V2V' = r2i(dst')*i2r(dst)
    //  (here dst is from the lta, and dst' from mri_dst)
    //
    //  Note that in order to obtain src->dst' with r2r, you
    // "don't need any info from dst (only from dst')"
    //  since
    //          src->dst'  = r2i(dst')* r2r * i2r(src)
    //
    //  However, using v2v, you need the information from dst
    //  since
    //          src->dst'  = r2i(dst')* i2r(dst) * v2v
    //
    //  Similarly src geometry can be different, therefore we convert
    //  any passed v2v lta into r2r first and apply it on possible different
    //  image geometries src and dst.
    //
    ////////////////////////////////////////////////////////////////////////

    // when the mri_dst volume is not given
    if (!mri_dst) {
      if (tran->dst.valid == 1)  // transform dst is valid
      {
        // modify dst geometry using the transform dst value,
        // i.e. put the head in the same position in the dst'
        // volume as the transform dst was
        //         if (DIAG_VERBOSE_ON)
        //           fprintf(stderr, "INFO: Modifying dst geometry, "
        //                   "using the transform dst\n");
        // allocate dst space (take type from src and geometry from transform):
        mri_dst = MRIallocSequence(
            tran->dst.width, tran->dst.height, tran->dst.depth, bspline->srctype, bspline->coeff->nframes);
        // cp rest of header information from src:
        MRIcopyHeader(bspline->coeff, mri_dst);
        mri_dst->type = bspline->srctype;
        // make sure the geometry is taken from the transform, not from src:
        useVolGeomToMRI(&tran->dst, mri_dst);
      }
      else if (getenv("USE_AVERAGE305")) {
        fprintf(stderr,
                "INFO: Environmental variable "
                "USE_AVERAGE305 set\n");
        fprintf(stderr,
                "INFO: Modifying dst c_(r,a,s), "
                "using average_305 values\n");
        // use the same volume size as the src
        mri_dst = MRIcloneDifferentType(bspline->coeff, bspline->srctype);
        // reset talairach transform file name:
        mri_dst->transform_fname[0] = '\0';
        mri_dst->c_r = -0.0950;
        mri_dst->c_a = -16.5100;
        mri_dst->c_s = 9.7500;
        mri_dst->ras_good_flag = 1;
        // maye one should set also the other geometry entries
        // from the average ???
        //
        // now we cache transform and thus we have to
        // do the following whenever
        // we change direction cosines
        MRIreInitCache(mri_dst);
      }
      else {
        fprintf(stderr,
                "INFO: Transform dst volume "
                "info is not used (valid flag = 0).\n");
        // use the same volume size as the src
        mri_dst = MRIcloneDifferentType(bspline->coeff, bspline->srctype);
        // reset talairach transform file name:
        mri_dst->transform_fname[0] = '\0';
        // maybe also reset or concatenate the actual transform
        // if available in mri_dst->transform (not yet implemented) ...
      }
    }

    ////////////////////////////////////////////////////////////////////////
    if (lta->type == LINEAR_RAS_TO_RAS) {
      // don't need any info from dst(in lta) only from mri_dst:
      return (MRIapplyRASlinearTransformBSpline(bspline, mri_dst, lta->xforms[0].m_L));
    }
    else if (lta->type == LINEAR_VOX_TO_VOX)  // vox-to-vox
    {
      // convert to ras_to_ras using xforms from lta if available
      // this strips geometry information and allows to use the
      // lta on volumes with different geometries
      // if lta geomery information is missing (should not happen)
      // we use the geometry from the passed source and target image
      if (lta->xforms[0].dst.valid) {
        i2r = vg_i_to_r(&lta->xforms[0].dst);  // allocated
      }
      else {
        fprintf(stderr,
                "INFO: LTA dst geometry information missing!\n"
                "      We assume that the dst volume passed is the\n"
                "      same as the dst for the transform.\n");
        i2r = extract_i_to_r(mri_dst);
      }
      tmp = MatrixMultiply(i2r, lta->xforms[0].m_L, NULL);
      if (lta->xforms[0].src.valid) {
        r2i = vg_r_to_i(&lta->xforms[0].src);  // allocated
      }
      else {
        fprintf(stderr,
                "INFO: LTA src geometry information missing!\n"
                "      We assume that the src volume passed is the\n"
                "      same as the src for the transform.\n");
        r2i = extract_r_to_i(bspline->coeff);
      }
      tmp = MatrixMultiply(tmp, r2i, tmp);

      resMRI = MRIapplyRASlinearTransformBSpline(bspline, mri_dst, tmp);
      MatrixFree(&i2r);
      MatrixFree(&r2i);
      MatrixFree(&tmp);
      return resMRI;

      // old code only treated different dst not different src geometries
      /*if (lta->xforms[0].dst.valid)
      {
        i2r = vg_i_to_r(&lta->xforms[0].dst); // allocated
        r2i = extract_r_to_i(mri_dst);
        tmp = MatrixMultiply(i2r, lta->xforms[0].m_L, NULL);
        v2v = MatrixMultiply(r2i, tmp, NULL);
        resMRI = MRIlinearTransformBSpline(bspline, mri_dst, v2v);
        MatrixFree(&v2v);
        v2v = 0;
        MatrixFree(&i2r);
        i2r = 0;
        MatrixFree(&r2i);
        r2i = 0;
        MatrixFree(&tmp);
        tmp = 0;
        return resMRI;
      }
      else
      {
        fprintf(stderr, "INFO: assumes that the dst volume "
                "given is the same as the dst for the transform\n");
        return(MRIlinearTransformBSpline(bspline, mri_dst,
                                        lta->xforms[0].m_L)) ;
      }*/
    }
    else if (lta->type == LINEAR_PHYSVOX_TO_PHYSVOX) {
      // must have both transform src and dst geometry information
      LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      return (MRIapplyRASlinearTransformBSpline(bspline, mri_dst, lta->xforms[0].m_L));
    }
    else
      ErrorExit(ERROR_BADPARM, "LTAtransformBSpline: unknown linear transform\n");
  }

  ErrorExit(ERROR_BADPARM, "LTAtransformBSpline: numxforms != 1 not supported\n");

  // removed stuff here

  return (mri_dst);
}

// -------------------------------------------------------------------------------
//
//            UP DOWN SAMPLE
//
//--------------------------------------------------------------------------------

#define MAXF 200L                           /* Maximum size of the filter */
#define SPLINE "Spline"                     /* Spline filter (l2-norm) */
#define SPLINE_L2 "Spline L2"               /* Spline filter (L2-norm) */
#define SPLINE_CENT "Centered Spline"       /* Centered Spline filter (l2-norm) */
#define SPLINE_CENT_L2 "Centered Spline L2" /* Centered Spline filter (L2-norm) */

extern MRI *MRIdownsample2BSpline(const MRI *mri_src, MRI *mri_dst)
{
  double *InBuffer;  /* Input buffer to 1D process */
  double *OutBuffer; /* Output buffer to 1D process */
  double g[MAXF];    /* Coefficients of the reduce filter */
  long ng;           /* Number of coefficients of the reduce filter */
  double h[MAXF];    /* Coefficients of the expansion filter */
  long nh;           /* Number of coefficients of the expansion filter */
  short IsCentered;  /* Equal TRUE if the filter is a centered spline, FALSE otherwise */
  int kx, ky, kz, kf;

  /* Get the filter coefficients for the Spline (order = 3) filter*/
  if (!GetPyramidFilter(SPLINE_CENT, 3, g, &ng, h, &nh, &IsCentered)) {
    printf("MRIdownsample2BSpline ERROR: Unable to load the filter coefficients\n");
    exit(1);
  }

  int NxIn = mri_src->width;
  int NyIn = mri_src->height;
  int NzIn = mri_src->depth;
  int NfIn = mri_src->nframes;
  int NxOut = NxIn / 2;
  if (NxOut < 1) NxOut = 1;
  int NyOut = NyIn / 2;
  if (NyOut < 1) NyOut = 1;
  int NzOut = NzIn / 2;
  if (NzOut < 1) NzOut = 1;

  // MRIwrite(mri_src,"mrisrc.mgz");
  /* --- X processing --- */
  MRI *mri_tmp = MRIallocSequence(NxOut, NyIn, NzIn, MRI_FLOAT, NfIn);
  if (!mri_tmp) ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate tmp mri\n");
  InBuffer = (double *)malloc((size_t)(NxIn * (long)sizeof(double)));
  if (InBuffer == (double *)NULL) {
    MRIfree(&mri_tmp);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate inbuffer mem\n");
  }
  OutBuffer = (double *)malloc((size_t)(NxOut * (long)sizeof(double)));
  if (OutBuffer == (double *)NULL) {
    MRIfree(&mri_tmp);
    free(InBuffer);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate outbuffer mem\n");
  }
  for (kf = 0; kf < NfIn; kf++)
    for (kz = 0; kz < NzIn; kz++)
      for (ky = 0; ky < NyIn; ky++) {
        if (NxIn > 1) {
          getXLine(mri_src, ky, kz, kf, InBuffer);
          Reduce_1D(InBuffer, NxIn, OutBuffer, g, ng, IsCentered);
          // printf(" f %i  z %i  y %i  NxOut %i  width %i\n",kf,kz,ky,NxOut,mri_tmp->width);
          setXLine(mri_tmp, ky, kz, kf, OutBuffer);
        }
        else {
          getXLine(mri_src, ky, kz, kf, InBuffer);
          setXLine(mri_tmp, ky, kz, kf, InBuffer);
        }
      }
  free(InBuffer);
  free(OutBuffer);
  // MRIwrite(mri_tmp,"mri_tmp1.mgz");

  /* --- Y processing --- */
  MRI *mri_tmp2 = MRIallocSequence(NxOut, NyOut, NzIn, MRI_FLOAT, NfIn);
  if (!mri_tmp2) ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate tmp mri\n");
  InBuffer = (double *)malloc((size_t)(NyIn * (long)sizeof(double)));
  if (InBuffer == (double *)NULL) {
    MRIfree(&mri_tmp);
    MRIfree(&mri_tmp2);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate inbuffer mem\n");
  }
  OutBuffer = (double *)malloc((size_t)(NyOut * (long)sizeof(double)));
  if (OutBuffer == (double *)NULL) {
    MRIfree(&mri_tmp);
    MRIfree(&mri_tmp2);
    free(InBuffer);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate outbuffer mem\n");
  }
  for (kf = 0; kf < NfIn; kf++)
    for (kz = 0; kz < NzIn; kz++)
      for (kx = 0; kx < NxOut; kx++) {
        if (NyIn > 1) {
          getYLine(mri_tmp, kx, kz, kf, InBuffer);
          Reduce_1D(InBuffer, NyIn, OutBuffer, g, ng, IsCentered);
          setYLine(mri_tmp2, kx, kz, kf, OutBuffer);
        }
        else {
          getYLine(mri_tmp, kx, kz, kf, InBuffer);
          setYLine(mri_tmp2, kx, kz, kf, InBuffer);
        }
      }
  free(InBuffer);
  free(OutBuffer);
  MRIfree(&mri_tmp);
  // MRIwrite(mri_tmp2,"mri_tmp2.mgz");

  /* --- Z processing --- */
  int allocdst = 0;
  if (!mri_dst) {
    mri_dst = MRIallocSequence(NxOut, NyOut, NzOut, mri_src->type, NfIn);
    // mri_dst = MRIallocSequence(NxOut, NyOut, NzOut, MRI_FLOAT, NfIn) ;
    MRIcopyHeader(mri_src, mri_dst);
    allocdst = 1;
  }
  if (mri_dst->width != NxOut || mri_dst->height != NyOut || mri_dst->depth != NzOut || mri_dst->nframes != NfIn) {
    printf("ERROR MRIupsample2BSpline: MRI Dest dimensions not correct!\n");
    exit(1);
  }
  InBuffer = (double *)malloc((size_t)(NzIn * (long)sizeof(double)));
  if (InBuffer == (double *)NULL) {
    MRIfree(&mri_tmp2);
    if (allocdst == 1) MRIfree(&mri_dst);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate inbuffer mem\n");
  }
  OutBuffer = (double *)malloc((size_t)(NzOut * (long)sizeof(double)));
  if (OutBuffer == (double *)NULL) {
    MRIfree(&mri_tmp2);
    if (allocdst == 1) MRIfree(&mri_dst);
    free(InBuffer);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate outbuffer mem\n");
  }
  for (kf = 0; kf < NfIn; kf++)
    for (ky = 0; ky < NyOut; ky++)
      for (kx = 0; kx < NxOut; kx++) {
        if (NzIn > 1) {
          getZLine(mri_tmp2, kx, ky, kf, InBuffer);
          Reduce_1D(InBuffer, NzIn, OutBuffer, g, ng, IsCentered);
          setZLine(mri_dst, kx, ky, kf, OutBuffer);
        }
        else {
          getZLine(mri_tmp2, kx, ky, kf, InBuffer);
          setZLine(mri_dst, kx, ky, kf, InBuffer);
        }
      }
  free(InBuffer);
  free(OutBuffer);
  MRIfree(&mri_tmp2);

  mri_dst->imnr0 = mri_src->imnr0;
  mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1;
  mri_dst->xsize = mri_src->xsize * 2;
  mri_dst->ysize = mri_src->ysize * 2;
  mri_dst->zsize = mri_src->zsize * 2;
  mri_dst->thick = mri_src->thick * 2;
  mri_dst->ps = mri_src->ps * 2;

  // adjust cras
  // printf("COMPUTING new CRAS\n") ;
  VECTOR *C = VectorAlloc(4, MATRIX_REAL);
  // here divide by 2 (int) as odd images get clipped when downsampling
  VECTOR_ELT(C, 1) = mri_src->width / 2 + 0.5;
  VECTOR_ELT(C, 2) = mri_src->height / 2 + 0.5;
  VECTOR_ELT(C, 3) = mri_src->depth / 2 + 0.5;
  VECTOR_ELT(C, 4) = 1.0;
  MATRIX *V2R = extract_i_to_r(mri_src);
  MATRIX *P = MatrixMultiply(V2R, C, NULL);
  mri_dst->c_r = P->rptr[1][1];
  mri_dst->c_a = P->rptr[2][1];
  mri_dst->c_s = P->rptr[3][1];
  MatrixFree(&P);
  MatrixFree(&V2R);
  VectorFree(&C);

  MRIreInitCache(mri_dst);
  // printf("CRAS new: %2.3f %2.3f %2.3f\n",mri_dst->c_r,mri_dst->c_a,mri_dst->c_s);

  // MRIwrite(mri_dst,"mridst.mgz");
  // exit(1);

  return (mri_dst);
}

extern MRI *MRIupsample2BSpline(const MRI *mri_src, MRI *mri_dst)
{
  double *InBuffer;  /* Input buffer to 1D process */
  double *OutBuffer; /* Output buffer to 1D process */
  double g[MAXF];    /* Coefficients of the reduce filter */
  long ng;           /* Number of coefficients of the reduce filter */
  double h[MAXF];    /* Coefficients of the expansion filter */
  long nh;           /* Number of coefficients of the expansion filter */
  short IsCentered;  /* Equal TRUE if the filter is a centered spline, FALSE otherwise */
  int kx, ky, kz, kf;

  /* Get the filter coefficients for the Spline (order = 3) filter*/
  if (!GetPyramidFilter(SPLINE_CENT, 3, g, &ng, h, &nh, &IsCentered)) {
    printf("MRIupsample2BSpline ERROR: Unable to load the filter coeffiients\n");
    exit(1);
  }

  int NxIn = mri_src->width;
  int NyIn = mri_src->height;
  int NzIn = mri_src->depth;
  int NfIn = mri_src->nframes;
  int NxOut, NyOut, NzOut;
  if (NxIn <= 1)
    NxOut = 1L;
  else
    NxOut = NxIn * 2;
  if (NyIn <= 1)
    NyOut = 1;
  else
    NyOut = NyIn * 2;
  if (NzIn <= 1)
    NzOut = 1;
  else
    NzOut = NzIn * 2;

  int allocdst = 0;
  if (!mri_dst) {
    mri_dst = MRIallocSequence(NxOut, NyOut, NzOut, mri_src->type, NfIn);
    MRIcopyHeader(mri_src, mri_dst);
    allocdst = 1;
  }

  if (mri_dst->width != NxOut || mri_dst->height != NyOut || mri_dst->depth != NzOut || mri_dst->nframes != NfIn) {
    printf("ERROR MRIupsample2BSpline: MRI Dest dimensions not correct!\n");
    exit(1);
  }

  /* --- X processing --- */
  InBuffer = (double *)malloc((size_t)(NxIn * (long)sizeof(double)));
  if (InBuffer == (double *)NULL) {
    if (allocdst == 1) MRIfree(&mri_dst);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate inbuffer mem\n");
  }
  OutBuffer = (double *)malloc((size_t)(NxOut * (long)sizeof(double)));
  if (OutBuffer == (double *)NULL) {
    if (allocdst == 1) MRIfree(&mri_dst);
    free(InBuffer);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate outbuffer mem\n");
  }
  for (kf = 0; kf < NfIn; kf++)
    for (kz = 0; kz < NzIn; kz++)
      for (ky = 0; ky < NyIn; ky++) {
        if (NxIn > 1) {
          getXLine(mri_src, ky, kz, kf, InBuffer);
          Expand_1D(InBuffer, NxIn, OutBuffer, h, nh, IsCentered);
          setXLine(mri_dst, ky, kz, kf, OutBuffer);
        }
        else {
          getXLine(mri_src, ky, kz, kf, InBuffer);
          setXLine(mri_dst, ky, kz, kf, InBuffer);
        }
      }
  free(InBuffer);
  free(OutBuffer);

  /* --- Y processing --- */
  InBuffer = (double *)malloc((size_t)(NyIn * (long)sizeof(double)));
  if (InBuffer == (double *)NULL) {
    if (allocdst == 1) MRIfree(&mri_dst);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate inbuffer mem\n");
  }
  OutBuffer = (double *)malloc((size_t)(NyOut * (long)sizeof(double)));
  if (OutBuffer == (double *)NULL) {
    if (allocdst == 1) MRIfree(&mri_dst);
    free(InBuffer);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate outbuffer mem\n");
  }
  for (kf = 0; kf < NfIn; kf++)
    for (kz = 0; kz < NzIn; kz++)
      for (kx = 0; kx < NxOut; kx++) {
        if (NyIn > 1) {
          getYLineStop(mri_dst, kx, kz, kf, InBuffer, NyIn);
          Expand_1D(InBuffer, NyIn, OutBuffer, h, nh, IsCentered);
          setYLine(mri_dst, kx, kz, kf, OutBuffer);
        }
        else {
          getYLine(mri_src, kx, kz, kf, InBuffer);
          setYLine(mri_dst, kx, kz, kf, InBuffer);
        }
      }
  free(InBuffer);
  free(OutBuffer);

  /* --- Z processing --- */
  InBuffer = (double *)malloc((size_t)(NzIn * (long)sizeof(double)));
  if (InBuffer == (double *)NULL) {
    if (allocdst == 1) MRIfree(&mri_dst);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate inbuffer mem\n");
  }
  OutBuffer = (double *)malloc((size_t)(NzOut * (long)sizeof(double)));
  if (OutBuffer == (double *)NULL) {
    if (allocdst == 1) MRIfree(&mri_dst);
    free(InBuffer);
    ErrorExit(ERROR_NO_MEMORY, "MRIdownsample2BSpline: could not allocate outbuffer mem\n");
  }
  for (kf = 0; kf < NfIn; kf++)
    for (ky = 0; ky < NyOut; ky++)
      for (kx = 0; kx < NxOut; kx++) {
        if (NzIn > 1) {
          getZLineStop(mri_dst, kx, ky, kf, InBuffer, NzIn);
          Expand_1D(InBuffer, NzIn, OutBuffer, h, nh, IsCentered);
          setZLine(mri_dst, kx, ky, kf, OutBuffer);
        }
        else {
          getZLine(mri_src, kx, ky, kf, InBuffer);
          setZLine(mri_dst, kx, ky, kf, InBuffer);
        }
      }
  free(InBuffer);
  free(OutBuffer);

  mri_dst->imnr0 = mri_src->imnr0;
  mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1;

  mri_dst->xsize = mri_src->xsize / 2.0;
  mri_dst->ysize = mri_src->ysize / 2.0;
  mri_dst->zsize = mri_src->zsize / 2.0;

  // adjust cras
  // printf("COMPUTING new CRAS\n") ;
  VECTOR *C = VectorAlloc(4, MATRIX_REAL);
  // here divide by 2.0 (because also odd images get fully upsampled)
  VECTOR_ELT(C, 1) = mri_src->width / 2.0 - 0.25;
  VECTOR_ELT(C, 2) = mri_src->height / 2.0 - 0.25;
  VECTOR_ELT(C, 3) = mri_src->depth / 2.0 - 0.25;
  VECTOR_ELT(C, 4) = 1.0;
  MATRIX *V2R = extract_i_to_r(mri_src);
  MATRIX *P = MatrixMultiply(V2R, C, NULL);
  mri_dst->c_r = P->rptr[1][1];
  mri_dst->c_a = P->rptr[2][1];
  mri_dst->c_s = P->rptr[3][1];
  MatrixFree(&P);
  MatrixFree(&V2R);
  VectorFree(&C);

  MRIreInitCache(mri_dst);

  return (mri_dst);
}

// /* ----------------------------------------------------------------------------
//
// 	Function:
// 		Reduce_2D
//
// 	Purpose:
//  		Reduces an image by a factor of two in each dimension.
//
// 	Note:
//  		Expects the output array (Out) to be allocated.
//
// 	Parameters:
// 		Input image:  	In[NxIn*NyIn]
// 		Output image: 	Out[NxIn/2*NyIn/2]
// 		Filter:			g[ng] coefficients of the filter
//
// ---------------------------------------------------------------------------- */
// extern int Reduce_2D(
// 				float *In, long NxIn, long NyIn,
// 				float *Out,
// 				double g[], long ng,
// 				short IsCentered
// 				)
// {
// float	*Tmp;
// double 	*InBuffer;		/* Input buffer to 1D process */
// double	*OutBuffer;		/* Output buffer to 1D process */
// long	kx, ky;
// long 	NxOut;
// long 	NyOut;
//
// 	/* --- Define dimension of the output --- */
// 	NxOut = NxIn/2L;
// 	if (NxOut < 1L) NxOut = 1L;
//
// 	NyOut = NyIn/2L;
// 	if (NyOut < 1L) NyOut = 1L;
//
// 	/* --- Allocate a temporary image --- */
// 	Tmp = (float *)malloc((size_t)(NxOut*NyIn*(long)sizeof(float)));
// 	if (Tmp == (float *)NULL) {
// 		printf("Unable to allocate memory\n");
// 		exit(1);
// 	}
//
// 	/* --- X processing --- */
// 	if (NxIn > 1L) {
// 		InBuffer = (double *)malloc((size_t)(NxIn*(long)sizeof(double)));
// 		if (InBuffer == (double *)NULL) {
// 			free(Tmp);
// 			printf("Unable to allocate memory");
// 			return(ERROR);
// 		}
// 		OutBuffer = (double *)malloc((size_t)(NxOut*(long)sizeof(double)));
// 		if (OutBuffer == (double *)NULL) {
// 			free(Tmp);
// 			free(InBuffer);
// 			printf("Unable to allocate memory");
// 			return(ERROR);
// 		}
// 		for (ky=0L; ky<NyIn; ky++) {
// 			GetRow( In, NxIn, NyIn, ky, InBuffer, NxIn);
// 			Reduce_1D(InBuffer, NxIn, OutBuffer, g, ng, IsCentered);
// 	    	PutRow( Tmp, NxOut, NyIn, ky, OutBuffer, NxOut);
//     	}
// 		free(InBuffer);
// 		free(OutBuffer);
//     }
//     else
// 	   	memcpy( Tmp, In, (size_t)(NyIn*(long)sizeof(float)));
//
// 	/* --- Y processing --- */
// 	if (NyIn > 1L) {
// 		InBuffer = (double *)malloc((size_t)(NyIn*(long)sizeof(double)));
// 		if (InBuffer == (double *)NULL) {
// 			printf("Unable to allocate memory");
// 			return(ERROR);
// 		}
// 		OutBuffer = (double *)malloc((size_t)(NyOut*(long)sizeof(double)));
// 		if (OutBuffer == (double *)NULL) {
// 			free(InBuffer);
// 			printf("Unable to allocate memory");
// 			return(ERROR);
// 		}
// 	    for (kx=0L; kx<NxOut; kx++) {
// 			GetColumn( Tmp, NxOut, NyIn, kx, InBuffer, NyIn);
// 			Reduce_1D(InBuffer, NyIn, OutBuffer, g, ng,IsCentered);
// 	 	   	PutColumn( Out, NxOut, NyOut, kx, OutBuffer, NyOut);
//    		}
// 		free(InBuffer);
// 		free(OutBuffer);
// 	}
// 	else
// 	   	memcpy( Out, Tmp, (size_t)(NxOut*(long)sizeof(float)));
//
// 	/* --- Free the temporary image --- */
// 	free(Tmp);
//
// 	return (!ERROR);
// }
//
//
//
/* ----------------------------------------------------------------------------
//
// 	Function:
// 		Expand_2D
//
// 	Purpose:
//  		Expands an image by a factor of two in each dimension.
//
// 	Note:
//  		Expects the output array (Out) to be allocated.
//
// 	Parameters:
// 		Input volume:  	In[NxIn,NyIn]
// 		Output voulme: 	Out[NxIn*2,NyIn*2]
// 		Filter coef:	h[nh]
//
// ---------------------------------------------------------------------------- */
// extern int Expand_2D(
// 				float *In, long NxIn, long NyIn,
// 				float *Out,
// 				double h[], long nh,
// 				short IsCentered
// 				)
// {
// double  *InBuffer; 		/* Input buffer to 1D process */
// double  *OutBuffer;		/* Output buffer to 1D process */
// long 	kx, ky;
// long 	NxOut;
// long 	NyOut;
//
// 	if (NxIn <= 1L)
// 		NxOut = 1L;
// 	else
// 		NxOut = NxIn*2L;
//
// 	if (NyIn <= 1L)
// 		NyOut = 1L;
// 	else
// 		NyOut = NyIn*2L;
//
// 	/* --- X processing --- */
// 	if (NxIn > 1L) {
// 		InBuffer = (double *)malloc((size_t)(NxIn*(long)sizeof(double)));
// 		if (InBuffer == (double *)NULL) {
// 			printf("Unable to allocate memory");
// 			return(ERROR);
// 		}
// 		OutBuffer = (double *)malloc((size_t)(NxOut*(long)sizeof(double)));
// 		if (OutBuffer == (double *)NULL) {
// 			free(InBuffer);
// 			printf("Unable to allocate memory");
// 			return(ERROR);
// 		}
// 		for (ky=0L; ky<NyIn; ky++) {
// 			GetRow( In, NxIn, NyIn, ky, InBuffer, NxIn);
// 			Expand_1D(InBuffer, NxIn, OutBuffer , h, nh, IsCentered);
// 			PutRow( Out, NxOut, NyOut, ky, OutBuffer, NxOut);
// 		}
// 		free(InBuffer);
// 		free(OutBuffer);
// 	}
// 	else
// 	   	memcpy( Out, In, (size_t)(NxIn*NyIn*(long)sizeof(float)));
//
//
// 	/* --- Y processing --- */
// 	if (NyIn > 1L) {
// 		InBuffer = (double *)malloc((size_t)(NyIn*(long)sizeof(double)));
// 		if (InBuffer == (double *)NULL) {
// 			printf("Unable to allocate memory");
// 			return(ERROR);
// 		}
// 		OutBuffer = (double *)malloc((size_t)(NyOut*(long)sizeof(double)));
// 		if (OutBuffer == (double *)NULL) {
// 			free(InBuffer);
// 			printf("Unable to allocate memory");
// 			return(ERROR);
// 		}
// 		for (kx=0L; kx<NxOut; kx++) {
// 			GetColumn( Out, NxOut, NyOut, kx, InBuffer, NyIn);
// 			Expand_1D(InBuffer, NyIn, OutBuffer, h, nh, IsCentered);
//     		PutColumn( Out, NxOut, NyOut, kx, OutBuffer, NyOut);
// 		}
// 		free(InBuffer);
// 		free(OutBuffer);
// 	}
//
// 	return (ERROR);
//
// }
//
//
