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



#ifndef RANDOMFIELDS_H
#define RANDOMFIELDS_H

#include "numerics.h"
#include "mri.h"

#define RF_UNIFORM  1
#define RF_GAUSSIAN 2
#define RF_Z 3
#define RF_T 4
#define RF_F 5
#define RF_CHI2 6

typedef struct
{
  char   *name;
  int     code;
  int nparams;
  double params[20];
  double mean, stddev;
  sc_rng *rng;
  const sc_rng_type *rngtype;
  unsigned long int seed;
}
RANDOM_FIELD_SPEC, RFS;


RFS *RFspecInit(unsigned long int seed, sc_rng_type *rngtype);
int RFspecFree(RFS **prfs);
int RFname2Code(RFS *rfs);
const char *RFcode2Name(RFS *rfs);
int RFprint(FILE *fp, RFS *rfs);
int RFspecSetSeed(RFS *rfs,unsigned long int seed);
int RFnparams(RFS *rfs);
int RFexpectedMeanStddev(RFS *rfs);
int RFsynth(MRI *rf, RFS *rfs, MRI *binmask);
MRI *RFstat2P(MRI *rf, RFS *rfs, MRI *binmask, int TwoSided, MRI *p);
MRI *RFz2p(MRI *z, MRI *mask, int TwoSided, MRI *p);
MRI *RFp2Stat(MRI *rf, RFS *rfs, MRI *binmask, MRI *p);
MRI *RFstat2Stat(MRI *rfin, RFS *rfsin, RFS *rfsout, MRI *binmask, MRI *rfout);
int RFglobalStats(MRI *rf, MRI *binmask,
                  double *gmean, double *gstddev, double *max);
MRI *RFrescale(MRI *rf, RFS *rfs, MRI *binmask, MRI *rfout);

double RFdrawVal(RFS *rfs);
double RFstat2PVal(RFS *rfs, double v);
double RFp2StatVal(RFS *rfs, double p);

int RFexpectedMeanStddevUniform(RFS *rfs);
int RFexpectedMeanStddevGaussian(RFS *rfs);
int RFexpectedMeanStddevt(RFS *rfs);
int RFexpectedMeanStddevF(RFS *rfs);
int RFexpectedMeanStddevChi2(RFS *rfs);

double RFar1ToGStd(double ar1, double d);
double RFar1ToFWHM(double ar1, double d);

double RFprobZCluster(double clustersize, double vthresh,
                      double fwhm, double searchsize, int dim);
double RFprobZClusterPThresh(double clustersize, double vpthresh,
                             double fwhm, double searchsize, int dim);
double RFprobZClusterSigThresh(double clustersize, double vsigthresh,
                               double fwhm, double searchsize, int dim);

MRI *RFp2z(MRI *p, MRI *mask, MRI *z);
MRI *RFz1toz2(MRI *z1, MRI *mask, MRI *z2);
#endif



