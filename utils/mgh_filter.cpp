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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "typedefs.h"

#define SWAP(a, b) \
  tempr = (a);     \
  (a) = (b);       \
  (b) = tempr

/* Prototypes     */
static void four1(float data[], int nn, int isign);
static void bandpass(float data[], int nn, float lo, float hi);
void bpfilter(FLOATTYPE **data, int nchan, int nsamp, float lo, float hi);
/* End prototypes */

static void four1(float data[], int nn, int isign)
{
  int n, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta;
  float tempr, tempi;

  n = nn << 1;
  j = 1;
  for (i = 1; i < n; i += 2) {
    if (j > i) {
      SWAP(data[j], data[i]);
      SWAP(data[j + 1], data[i + 1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = 2 * mmax;
    theta = 6.28318530717959 / (isign * mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr = wr * data[j] - wi * data[j + 1];
        tempi = wr * data[j + 1] + wi * data[j];
        data[j] = data[i] - tempr;
        data[j + 1] = data[i + 1] - tempi;
        data[i] += tempr;
        data[i + 1] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}

#if 0
static void lowpass(float data[],int nn,float hi);
static void lowpass(float data[],int nn,float hi)
{
  float norm,f,w,wc;
  int i;

  norm = 1/(sqrt(2*3.1415926535)*hi);
  wc = hi;
  four1(data-1,nn,1);
  for (i=0;i<2*nn;i+=2)
  {
    w = nn/2-fabs((float)i/2-nn/2);
    /*    f = norm*exp(-(nn-fabs(i/2-nn))/(2*hi)); */
    f = 1/(1+pow(w/wc,3.0));
    data[i] *= f;
    data[i+1] *= f;
  }
  four1(data-1,nn,-1);
}
#endif

static void bandpass(float data[], int nn, float lo, float hi)
{
  // float norm,
  float f, fl, fh, w, wh, wl;
  int i;

  // norm = 1 / (sqrt(2 * M_PI) * hi);
  wh = hi;
  wl = lo;
  four1(data - 1, nn, 1);
  for (i = 0; i < 2 * nn; i += 2) {
    w = nn / 2 - fabs((float)i / 2 - nn / 2);
    /*    f = norm*exp(-(nn-fabs(i/2-nn))/(2*hi)); */
    /*    f = 1/(1+pow(wl/(w+0.001),3.0))*1/(1+pow(w/wh,3.0)); */
    fl = (wl > 0) ? 1 - 1 / (1 + pow(w / wl, 6.0)) : 1;
    fh = (wh > 0) ? 1 / (1 + pow(w / wh, 6.0)) : 0;
    f = sqrt(fl * fh);
    data[i] *= f;
    data[i + 1] *= f;
  }
  four1(data - 1, nn, -1);
}

void bpfilter(FLOATTYPE **data, int nchan, int nsamp, float lo, float hi)
{
  float *tmpvec;
  int i, j;

  tmpvec = (float *)malloc(nsamp * 2 * sizeof(float));
  for (i = 0; i < nchan; i++) {
    for (j = 0; j < nsamp; j++) {
      tmpvec[2 * j] = data[i][j];
      tmpvec[2 * j + 1] = 0;
    }
    bandpass(tmpvec, nsamp, lo, hi);
    for (j = 0; j < nsamp; j++) {
      data[i][j] = tmpvec[2 * j] / nsamp;
    }
  }
  free(tmpvec);
}
