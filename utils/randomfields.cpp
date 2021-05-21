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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "fsglm.h"
#include "mri.h"
#include "mri2.h"
#include "mrisurf.h"
#include "pdf.h"
#include "randomfields.h"
#include "utils.h"

#undef SIGN
#define SIGN(x) (((x) > 0) ? 1.0 : -1.0)

/*-------------------------------------------------------------------*/
int RFname2Code(RFS *rfs)
{
  int code = -1;
  if (!strcmp(rfs->name, "uniform")) code = RF_UNIFORM;
  if (!strcmp(rfs->name, "gaussian")) code = RF_GAUSSIAN;
  if (!strcmp(rfs->name, "z")) code = RF_Z;
  if (!strcmp(rfs->name, "t")) code = RF_T;
  if (!strcmp(rfs->name, "F")) code = RF_F;
  if (!strcmp(rfs->name, "chi2")) code = RF_CHI2;
  rfs->code = code;
  return (code);
}
/*-------------------------------------------------------------------*/
const char *RFcode2Name(RFS *rfs)
{
  switch (rfs->code) {
    case RF_UNIFORM:
      return ("uniform");
      break;
    case RF_GAUSSIAN:
      return ("gaussian");
      break;
    case RF_Z:
      return ("z");
      break;
    case RF_T:
      return ("t");
      break;
    case RF_F:
      return ("F");
      break;
    case RF_CHI2:
      return ("chi2");
      break;
  }
  return (NULL);
}

/*-------------------------------------------------------------------*/
RFS *RFspecInit(unsigned long int seed, sc_rng_type *rngtype)
{
  RFS *rfs;
  const sc_rng_type *sc_rng_intern_type = &intern_rng_type;

  rfs = (RFS *)calloc(sizeof(RFS), 1);

  /* Sc: only one type of rng is supported at this moment (ranlux389) */
  rfs->rngtype = sc_rng_intern_type;

  rfs->rng = sc_rng_alloc(rfs->rngtype);

  RFspecSetSeed(rfs, seed);
  return (rfs);
}

/*-------------------------------------------------------------------*/
int RFspecFree(RFS **prfs)
{
  sc_rng_free((*prfs)->rng);
  free((*prfs)->name);
  free(*prfs);
  *prfs = NULL;
  return (0);
}

/*-------------------------------------------------------------------*/
int RFspecSetSeed(RFS *rfs, unsigned long int seed)
{
  if (seed == 0)
    rfs->seed = PDFtodSeed();
  else
    rfs->seed = seed;
  sc_rng_set(rfs->rng, rfs->seed);
  return (0);
}

/*-------------------------------------------------------------------*/
int RFprint(FILE *fp, RFS *rfs)
{
  int n;
  fprintf(fp, "field %s\n", rfs->name);
  fprintf(fp, "code  %d\n", rfs->code);
  fprintf(fp, "nparams %d\n", rfs->nparams);
  for (n = 0; n < rfs->nparams; n++) fprintf(fp, " param %d  %lf\n", n, rfs->params[n]);
  fprintf(fp, "mean %lf\n", rfs->mean);
  fprintf(fp, "stddev %lf\n", rfs->stddev);
  fprintf(fp, "seed   %ld\n", rfs->seed);
  return (0);
}

/*-------------------------------------------------------------------*/
int RFnparams(RFS *rfs)
{
  int nparams = -1;
  switch (rfs->code) {
    case RF_UNIFORM:
      nparams = 2;
      break;
    case RF_GAUSSIAN:
      nparams = 2;
      break;
    case RF_Z:
      nparams = 0;
      break;
    case RF_T:
      nparams = 1;
      break;
    case RF_F:
      nparams = 2;
      break;
    case RF_CHI2:
      nparams = 1;
      break;
  }
  rfs->nparams = nparams;
  return (nparams);
}

/*-------------------------------------------------------------------*/
int RFexpectedMeanStddev(RFS *rfs)
{
  if (!strcmp(rfs->name, "uniform")) return (RFexpectedMeanStddevUniform(rfs));
  if (!strcmp(rfs->name, "gaussian")) return (RFexpectedMeanStddevGaussian(rfs));
  if (!strcmp(rfs->name, "t")) return (RFexpectedMeanStddevt(rfs));
  if (!strcmp(rfs->name, "F")) return (RFexpectedMeanStddevt(rfs));
  if (!strcmp(rfs->name, "chi2")) return (RFexpectedMeanStddevChi2(rfs));
  printf("ERROR: RFexpectedMeanStddev(): field type %s unknown\n", rfs->name);
  return (1);
}

/*-------------------------------------------------------------------*/
int RFsynth(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c, r, s, f;
  double v, m;

  if (RFname2Code(rfs) == -1) return (1);

  for (c = 0; c < rf->width; c++) {
    for (r = 0; r < rf->height; r++) {
      for (s = 0; s < rf->depth; s++) {
        if (binmask != NULL) {
          m = MRIgetVoxVal(binmask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        for (f = 0; f < rf->nframes; f++) {
          v = RFdrawVal(rfs);
          MRIsetVoxVal(rf, c, r, s, f, v);
        }
      }
    }
  }
  return (0);
}
/*!
  \fn MRI *RFstat2P(MRI *rf, RFS *rfs, MRI *binmask, int TwoSided, MRI *p)
  \brief Converts a stat to a p value. If TwoSided, then computes a p value
      based on an unsigned stat, but the sign is still passed to p.
*/
MRI *RFstat2P(MRI *rf, RFS *rfs, MRI *binmask, int TwoSided, MRI *p)
{
  int c, r, s, f = 0, m;
  double v, pval;

  if (RFname2Code(rfs) == -1) return (NULL);
  p = MRIclone(rf, p);

  for (c = 0; c < rf->width; c++) {
    for (r = 0; r < rf->height; r++) {
      for (s = 0; s < rf->depth; s++) {
        if (binmask != NULL) {
          m = (int)MRIgetVoxVal(binmask, c, r, s, 0);
          if (!m) continue;
        }
        for (f = 0; f < rf->nframes; f++) {
          v = MRIgetVoxVal(rf, c, r, s, f);
          if (TwoSided)
            pval = SIGN(v) * 2 * RFstat2PVal(rfs, fabs(v));
          else
            pval = RFstat2PVal(rfs, v);
          MRIsetVoxVal(p, c, r, s, f, pval);
        }
      }
    }
  }
  return (p);
}
/*!
  \fn MRI *RFz2p(MRI *z, MRI *mask, int TwoSided, MRI *p)
  \brief Converts z to p. If TwoSided, then computes the unsigned p,
         but p keeps the sign of z.
*/
MRI *RFz2p(MRI *z, MRI *mask, int TwoSided, MRI *p)
{
  RFS *rfs;
  rfs = RFspecInit(0, NULL);
  rfs->name = strcpyalloc("gaussian");
  rfs->params[0] = 0;
  rfs->params[1] = 1;
  p = RFstat2P(z, rfs, mask, TwoSided, p);
  return (p);
}

/*-------------------------------------------------------------------*/
MRI *RFp2Stat(MRI *p, RFS *rfs, MRI *binmask, MRI *rf)
{
  int c, r, s, f, m;
  double v, pval;

  if (RFname2Code(rfs) == -1) return (NULL);
  rf = MRIclone(p, rf);

  for (c = 0; c < rf->width; c++) {
    for (r = 0; r < rf->height; r++) {
      for (s = 0; s < rf->depth; s++) {
        if (binmask != NULL) {
          m = (int)MRIgetVoxVal(binmask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        for (f = 0; f < rf->nframes; f++) {
          pval = MRIgetVoxVal(p, c, r, s, f);
          // use fabs() in case it has been signed
          if (1 - fabs(pval) < FLT_EPSILON)
            v = 0;  // not sure what to do here
          else
            v = RFp2StatVal(rfs, fabs(pval));
          MRIsetVoxVal(rf, c, r, s, f, v);
        }
      }
    }
  }
  return (rf);
}

/*--------------------------------------------------------------------------*/
MRI *RFstat2Stat(MRI *rfin, RFS *rfsin, RFS *rfsout, MRI *binmask, MRI *rfout)
{
  MRI *p = NULL;

  if (RFname2Code(rfsin) == -1) return (NULL);
  if (RFname2Code(rfsout) == -1) return (NULL);

  p = RFstat2P(rfin, rfsin, binmask, 0, p);
  rfout = RFp2Stat(p, rfsout, binmask, rfout);
  MRIfree(&p);
  return (rfout);
}

/*-------------------------------------------------------------------*/
MRI *RFrescale(MRI *rf, RFS *rfs, MRI *binmask, MRI *rfout)
{
  int c, r, s, f, m;
  double v, gmean, gstddev, gmax;

  if (RFname2Code(rfs) == -1) return (NULL);
  RFexpectedMeanStddev(rfs);                            // expected
  RFglobalStats(rf, binmask, &gmean, &gstddev, &gmax);  // actual

  rfout = MRIclone(rf, rfout);

  for (c = 0; c < rf->width; c++) {
    for (r = 0; r < rf->height; r++) {
      for (s = 0; s < rf->depth; s++) {
        if (binmask != NULL) {
          m = (int)MRIgetVoxVal(binmask, c, r, s, 0);
          if (!m) continue;
        }
        for (f = 0; f < rf->nframes; f++) {
          v = MRIgetVoxVal(rf, c, r, s, f);
          v = (v - gmean) * (rfs->stddev / gstddev) + rfs->mean;
          MRIsetVoxVal(rfout, c, r, s, f, v);
        }
      }
    }
  }
  return (rfout);
}

/*-------------------------------------------------------------------*/
int RFglobalStats(MRI *rf, MRI *binmask, double *gmean, double *gstddev, double *max)
{
  int c, r, s, f, m;
  double v;
  double sum, sumsq;
  long nv;

  nv = 0;
  sum = 0;
  sumsq = 0;
  *max = -1000000;
  for (c = 0; c < rf->width; c++) {
    for (r = 0; r < rf->height; r++) {
      for (s = 0; s < rf->depth; s++) {
        if (binmask != NULL) {
          m = (int)MRIgetVoxVal(binmask, c, r, s, 0);
          if (!m) continue;
        }
        for (f = 0; f < rf->nframes; f++) {
          v = MRIgetVoxVal(rf, c, r, s, f);
          sum += v;
          sumsq += (v * v);
          nv++;
          if (*max < v) *max = v;
        }
      }
    }
  }
  *gmean = sum / nv;
  *gstddev = sqrt(sumsq / nv - (*gmean) * (*gmean));

  return (0);
}

/*-------------------------------------------------------------------*/
double RFdrawVal(RFS *rfs)
{
  if (!strcmp(rfs->name, "uniform")) {
    // params[0] = min
    // params[1] = max
    return (sc_ran_flat(rfs->rng, rfs->params[0], rfs->params[1]));
  }
  if (!strcmp(rfs->name, "gaussian")) {
    // params[0] = mean
    // params[1] = std
    return (sc_ran_gaussian(rfs->rng, rfs->params[1]) + rfs->params[0]);
  }
  if (!strcmp(rfs->name, "z")) {
    // nparams=0, gaussian with mean=0, std=1
    return (sc_ran_gaussian(rfs->rng, 1));
  }
  if (!strcmp(rfs->name, "t")) {
    // params[0] = dof
    return (sc_ran_tdist(rfs->rng, rfs->params[0]));
  }
  if (!strcmp(rfs->name, "F")) {
    // params[0] = numerator dof (rows in C)
    // params[1] = denominator dof
    return (sc_ran_fdist(rfs->rng, rfs->params[0], rfs->params[1]));
  }
  if (!strcmp(rfs->name, "chi2")) {
    // params[0] = dof
    return (sc_ran_chisq(rfs->rng, rfs->params[0]));
  }
  printf("ERROR: RFdrawVal(): field type %s unknown\n", rfs->name);
  return (10000000000.0);
}

/*-------------------------------------------------------------------
  RFstat2PVal() - returns the probability of getting a value of stat
  or greater from a random draw from the given distribution. Note:
  this is a one-sided test (where sidedness makes sense).
  -------------------------------------------------------------------*/
double RFstat2PVal(RFS *rfs, double stat)
{
  double p = -1;
  if (!strcmp(rfs->name, "uniform")) {
    // params[0] = min
    // params[1] = max
    p = sc_cdf_flat_Q(stat, rfs->params[0], rfs->params[1]);
  }
  if (!strcmp(rfs->name, "gaussian")) {
    // params[0] = mean
    // params[1] = std
    p = sc_cdf_gaussian_Q(stat - rfs->params[0], rfs->params[1]);
  }
  if (!strcmp(rfs->name, "z")) {
    // nparams=0, gaussian with mean=0, std=1
    p = sc_cdf_gaussian_Q(stat, 1);
  }
  if (!strcmp(rfs->name, "t")) {
    // params[0] = dof
    p = sc_cdf_tdist_Q(stat, rfs->params[0]);
  }
  if (!strcmp(rfs->name, "F")) {
    // params[0] = numerator dof (rows in C)
    // params[1] = denominator dof
    p = sc_cdf_fdist_Q(stat, rfs->params[0], rfs->params[1]);
  }
  if (!strcmp(rfs->name, "chi2")) {
    // params[0] = dof
    p = sc_cdf_chisq_Q(stat, rfs->params[0]);
  }
  if (p == -1) {
    printf("ERROR: RFstat2PVal(): field type %s unknown\n", rfs->name);
    return (10000000000.0);
  }
  if (std::isinf(p) || p < FLT_MIN) p = FLT_MIN;
  return (p);
}

/*-------------------------------------------------------------------
  RFp2StatVal() -
  -------------------------------------------------------------------*/
double RFp2StatVal(RFS *rfs, double p)
{
  if (!strcmp(rfs->name, "uniform")) {
    // params[0] = min
    // params[1] = max
    return (sc_cdf_flat_Qinv(p, rfs->params[0], rfs->params[1]));
  }
  if (!strcmp(rfs->name, "gaussian")) {
    // params[0] = mean
    // params[1] = std
    return (sc_cdf_gaussian_Qinv(p, rfs->params[1]) + rfs->params[0]);
  }
  if (!strcmp(rfs->name, "z")) {
    // nparams=0, gaussian with mean=0, std=1
    return (sc_cdf_gaussian_Qinv(p, 1));
  }
  if (!strcmp(rfs->name, "t")) {
    // params[0] = dof
    return (sc_cdf_tdist_Qinv(p, rfs->params[0]));
  }
  if (!strcmp(rfs->name, "F")) {
    // params[0] = numerator dof (rows in C)
    // params[1] = denominator dof
    // GSL does not have the Qinv for F :<.
    printf("ERROR: RFp2StatVal(): cannot invert RF field\n");
    return (10000000000.0);
  }
  if (!strcmp(rfs->name, "chi2")) {
    // params[0] = dof
    return (sc_cdf_chisq_Qinv(p, rfs->params[0]));
  }
  printf("ERROR: RFp2v(): field type %s unknown\n", rfs->name);
  return (10000000000.0);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevUniform(RFS *rfs)
{
  double min, max, d;
  min = rfs->params[0];
  max = rfs->params[1];
  d = max - min;
  rfs->mean = d / 2;
  rfs->stddev = sqrt((d * d) / 12.0);
  return (0);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevGaussian(RFS *rfs)
{
  rfs->mean = rfs->params[0];
  rfs->stddev = rfs->params[1];
  return (0);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevt(RFS *rfs)
{
  double dof;
  dof = rfs->params[0];
  rfs->mean = 0;
  rfs->stddev = sqrt(dof / (dof - 2));
  return (0);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevF(RFS *rfs)
{
  double ndof, ddof;
  ndof = rfs->params[0];  // numerator dof (rows in C)
  ddof = rfs->params[1];  // dof
  rfs->mean = ddof / (ddof - 2);
  rfs->stddev = 2 * (ddof * ddof) * (ndof + ddof - 2) / (ndof * ((ddof - 2) * (ddof - 2)) * (ddof - 4));
  return (0);
}

/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevChi2(RFS *rfs)
{
  double dof;
  dof = rfs->params[0];
  rfs->mean = dof;
  rfs->stddev = sqrt(2.0 * dof);
  return (0);
}

/*------------------------------------------------------------------
  RFar1ToGStd() - converts AR1 value to equivalent Gaussian standard
  devation.  Ie, if white Gaussian noise where filtered with a
  gaussian kernel with gstd, it would result in the given AR1.  The
  AR1 value is measured between two voxels a distance d appart. GStd
  will have the same units as d.
  -----------------------------------------------------------------*/
double RFar1ToGStd(double ar1, double d)
{
  double gstd;
  if (ar1 <= 0.0) return (0.0);
  gstd = d / sqrt(-4 * log(ar1));
  return (gstd);
}

/*---------------------------------------------------------------------
  RFar1ToFWHM() - converts AR1 value to equivalent FWHM.  Ie, if white
  Gaussian noise where filtered with a gaussian kernel with FWHM, it
  would result in the given AR1.  The AR1 value is measured between
  two voxels a distance d appart. FWHM will have the same units as d.
  ------------------------------------------------------*/
double RFar1ToFWHM(double ar1, double d)
{
  double gstd, fwhm;
  gstd = RFar1ToGStd(ar1, d);
  fwhm = gstd * sqrt(log(256.0));
  return (fwhm);
}

/*!
 \fn double RFprobZCluster(double clustersize, double vthresh,
  double fwhm, double searchsize, int dim)
 \brief Probability of a cluster in a z-field being >= clustersize
 \param clustersize - in physical units
 \param vzthresh - voxel-wise z threshold
 \param fwhm - of the z field in physical units
 \param searchsize - in physical units
 \param dim - dimension of the z field

  Friston, Worsley, Frackowiak, Mazziotta, Evans. Assessing the
  significance of focal activations using their spatial extent. HBM
  1994, 1:214-220.

  By design, this gives the same results as rft_zcluster_cdf.m in
  the FSFAST toolbox.

 */
double RFprobZCluster(double clustersize, double vzthresh, double fwhm, double searchsize, int dim)
{
  double u, phiu, k, S, W, Em, beta, Pnk, pcluster, D, pi;

  pi = M_PI;
  D = dim;
  u = vzthresh;

  // Equivalent p-value threshold. Note that the paper uses phi(-u),
  // but the results dont work out that way. Note that phi(-u) can be
  // used if (1-phiu) is used instead of phiu in the equation for beta
  // below. See Hayasaka and Nichols 2003, App A, equation 4.
  phiu = sc_cdf_gaussian_Q(u, 1);

  k = clustersize;  // cluster size to test (actual units, not resels)
  S = searchsize;   // search space (actual units, not resels)

  W = fwhm / sqrt(4 * log(2));

  // Expected number of clusters (Eq 2)
  // This form appears to go back to Hasofer 1978
  if (dim != 3) {
    Em = exp(-(u * u) / 2) * pow(u, (D - 1)) * pow(2 * pi, -(D + 1) / 2) * S / pow(W, D);
  }
  else {
    // This form with (pow(u,(D-1))-1) for D=3 found in
    // Worsley, et al, 1996, HBM 4:58-73, Table II. This is what
    // FSL and KJW's stat_threshold.m use. In simulations, it seems
    // more accurate.
    Em = exp(-(u * u) / 2) * (pow(u, (D - 1)) - 1) * pow(2 * pi, -(D + 1) / 2) * S / pow(W, D);
  }

  // Equation 3
  beta = pow(tgamma(D / 2 + 1) * Em / (S * phiu), 2 / D);

  // Prob that number of voxels in a cluster (n) exceeds k (Bet Eq 2 and 3)
  Pnk = exp(-beta * pow(k, 2.0 / D));

  // Prob of cluster of size k
  pcluster = 1 - exp(-Em * Pnk);

#if 0
  printf("--------------------------\n");
  printf("csize = %lf\n",clustersize);
  printf("vzthresh = %lf\n",vzthresh);
  printf("fwhm = %lf\n",fwhm);
  printf("searchsize = %lf\n",searchsize);
  printf("dim = %d\n",dim);
  printf("W = %lf\n",W);
  printf("Em = %lf\n",Em);
  printf("phiu = %lf\n",phiu);
  printf("beta = %lf\n",beta);
  printf("Pnk = %lf\n",Pnk);
  printf("pcluster = %lf\n",pcluster);
#endif

  return (pcluster);
}

/*!
 \fn double RFprobZClusterPThresh(double clustersize,
  double vpthresh, double fwhm, double searchsize, int dim)
 \brief Probability of a cluster in a z-field being >= clustersize
 \param clustersize - in physical units
 \param vpthresh - voxel-wise p-threshold (instead of z)
 \param fwhm - of the z field in physical units
 \param searchsize - in physical units
 \param dim - dimension of the z field
 Same as RFprobZCluster() but takes the p-value threshold used
 to threshold a p-field created from a z-field. This function
 just computes the equivalenet z-threshold and calls
 RFprobZCluster(). Note that the cluster are the set of voxels
 *below* the pvthresh.
 */
double RFprobZClusterPThresh(double clustersize, double vpthresh, double fwhm, double searchsize, int dim)
{
  double vzthresh, pcluster;
  vzthresh = sc_cdf_gaussian_Qinv(vpthresh, 1.0);
  pcluster = RFprobZCluster(clustersize, vzthresh, fwhm, searchsize, dim);
  return (pcluster);
}

/*!
 \fn double RFprobZClusterSigThresh(double clustersize,
  double vsigthresh, double fwhm, double searchsize, int dim)
 \brief Probability of a cluster in a z-field being >= clustersize
 \param clustersize - in physical units
 \param vpthresh - voxel-wise sig-threshold (-log10(p))
 \param fwhm - of the z field in physical units
 \param searchsize - in physical units
 \param dim - dimension of the z field
 Same as RFprobZCluster() but takes the significance (ie, -log10(p))
 value threshold used to threshold a p-field created from a
 z-field. This function just computes the equivalenet z-threshold and
 calls RFprobZCluster().
 */
double RFprobZClusterSigThresh(double clustersize, double vsigthresh, double fwhm, double searchsize, int dim)
{
  double vzthresh, pcluster, vpthresh;
  vpthresh = pow(10.0, -fabs(vsigthresh));
  vzthresh = sc_cdf_gaussian_Qinv(vpthresh, 1.0);
  pcluster = RFprobZCluster(clustersize, vzthresh, fwhm, searchsize, dim);
  return (pcluster);
}

/*!
 \fn MRI *RFp2z(MRI *p, MRI *mask, MRI *z)
 \brief Converts a p-value to a two-sided z
 */
MRI *RFp2z(MRI *p, MRI *mask, MRI *z)
{
  RFS *rfs;
  rfs = RFspecInit(0, NULL);
  rfs->name = strcpyalloc("gaussian");
  rfs->params[0] = 0;
  rfs->params[1] = 1;
  z = RFp2Stat(p, rfs, mask, z);
  return (z);
}

/*!
 \fn MRI *RFz1toz2(MRI *z1, MRI *mask, MRI *z2)
 \brief Converts a one-sided z-map to a two-sided z-map
 */
MRI *RFz1toz2(MRI *z1, MRI *mask, MRI *z2)
{
  MRI *p;
  int TwoSidedFlag = 1;

  // Convert to two-sided p-values
  p = RFz2p(z1, mask, TwoSidedFlag, NULL);
  // Now convert it back to z
  z2 = RFp2z(p, mask, z2);

  MRIfree(&p);
  return (z2);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
