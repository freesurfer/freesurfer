// $Id: randomfields.h,v 1.4 2006/08/30 20:54:39 czanner Exp $

#ifndef RANDOMFIELDS_H
#define RANDOMFIELDS_H

#if USE_SC_GSL_REPLACEMENT
  #include <gsl_wrapper.h>
#else
  #include <gsl/gsl_cdf.h>
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_randist.h>
#endif

#include "mri.h"

#define RF_UNIFORM  1
#define RF_GAUSSIAN 2
#define RF_Z 3
#define RF_T 4 
#define RF_F 5
#define RF_CHI2 6

typedef struct {
  char   *name;
  int     code;
  int nparams;
  double params[20];
  double mean, stddev;
#if USE_SC_GSL_REPLACEMENT
  sc_rng *rng;
  const sc_rng_type *rngtype;
#else
  gsl_rng *rng;
  const gsl_rng_type *rngtype;
#endif
  unsigned long int seed;
} RANDOM_FIELD_SPEC, RFS;


const char *RFSrcVersion(void);
#if USE_SC_GSL_REPLACEMENT
  RFS *RFspecInit(unsigned long int seed, sc_rng_type *rngtype);
#else
  RFS *RFspecInit(unsigned long int seed, gsl_rng_type *rngtype);
#endif
int RFspecFree(RFS **prfs);
int RFname2Code(RFS *rfs);
const char *RFcode2Name(RFS *rfs);
int RFprint(FILE *fp, RFS *rfs);

int RFspecSetSeed(RFS *rfs,unsigned long int seed);
int RFnparams(RFS *rfs);
int RFexpectedMeanStddev(RFS *rfs);
int RFsynth(MRI *rf, RFS *rfs, MRI *binmask);
MRI *RFstat2P(MRI *rf, RFS *rfs, MRI *binmask, MRI *p);
MRI *RFp2Stat(MRI *rf, RFS *rfs, MRI *binmask, MRI *p);
MRI *RFstat2Stat(MRI *rfin, RFS *rfsin, RFS *rfsout, MRI *binmask, MRI *rfout);
int RFglobalStats(MRI *rf, MRI *binmask, double *gmean, double *gstddev, double *max);
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

#endif



