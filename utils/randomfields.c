#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <float.h>
#if USE_SC_GSL_REPLACEMENT
#else
  #include <gsl/gsl_cdf.h>
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_randist.h>
#endif
#include "randomfields.h"
#include "utils.h"
#include "mri.h"
#include "mri2.h"
#include "mrisurf.h"
#include "fsglm.h"
#include "pdf.h"

/* --------------------------------------------- */
// Return the CVS version of this file.
const char *RFSrcVersion(void) { 
  return("$Id: randomfields.c,v 1.4 2006/08/30 20:55:08 czanner Exp $"); 
}

/*-------------------------------------------------------------------*/
int RFname2Code(RFS *rfs)
{
  int code = -1;
  if(!strcmp(rfs->name,"uniform"))  code = RF_UNIFORM;
  if(!strcmp(rfs->name,"gaussian")) code = RF_GAUSSIAN;
  if(!strcmp(rfs->name,"z"))        code = RF_Z;
  if(!strcmp(rfs->name,"t"))        code = RF_T;
  if(!strcmp(rfs->name,"F"))        code = RF_F;
  if(!strcmp(rfs->name,"chi2"))     code = RF_CHI2;
  rfs->code = code;
  return(code);
}
/*-------------------------------------------------------------------*/
const char *RFcode2Name(RFS *rfs)
{
  switch(rfs->code){
  case RF_UNIFORM:  return("uniform"); break;
  case RF_GAUSSIAN: return("gaussian"); break;
  case RF_Z:        return("z"); break;
  case RF_T:        return("t"); break;
  case RF_F:        return("F"); break;
  case RF_CHI2:     return("chi2"); break;
  }
  return(NULL);
}

/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
RFS *RFspecInit(unsigned long int seed, sc_rng_type *rngtype)
{
  RFS *rfs;  
  const sc_rng_type *sc_rng_intern_type = &intern_rng_type;
 
  rfs = (RFS*)calloc(sizeof(RFS),1);
 
  /* Sc: only one type of rng is supported at this moment (ranlux389) */
  rfs->rngtype = sc_rng_intern_type;
 
  rfs->rng = sc_rng_alloc(rfs->rngtype);
 
  RFspecSetSeed(rfs,seed);
  return(rfs);
}
#else
RFS *RFspecInit(unsigned long int seed, gsl_rng_type *rngtype)
{
  RFS *rfs;

  rfs = (RFS*)calloc(sizeof(RFS),1);

  if(rngtype == NULL) rfs->rngtype = gsl_rng_ranlux389;
  else                rfs->rngtype = rngtype;
  rfs->rng = gsl_rng_alloc(rfs->rngtype);

  RFspecSetSeed(rfs,seed);
  return(rfs);
}
#endif
/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
int RFspecFree(RFS **prfs)
{
  sc_rng_free((*prfs)->rng);
  free((*prfs)->name);
  free(*prfs);
  *prfs =NULL;
  return(0);
}
#else
int RFspecFree(RFS **prfs)
{
  gsl_rng_free((*prfs)->rng);
  free((*prfs)->name);
  free(*prfs);
  *prfs =NULL;
  return(0);
}
#endif
/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
int RFspecSetSeed(RFS *rfs,unsigned long int seed)
{
  if(seed == 0) rfs->seed = PDFtodSeed();
  else          rfs->seed = seed;
  sc_rng_set(rfs->rng, rfs->seed);
  return(0);
}
#else
int RFspecSetSeed(RFS *rfs,unsigned long int seed)
{
  if(seed == 0) rfs->seed = PDFtodSeed();
  else          rfs->seed = seed;
  gsl_rng_set(rfs->rng, rfs->seed);
  return(0);
}
#endif
/*-------------------------------------------------------------------*/
int RFprint(FILE *fp, RFS *rfs)
{
  int n;
  fprintf(fp,"field %s\n",rfs->name);
  fprintf(fp,"code  %d\n",rfs->code);
  fprintf(fp,"nparams %d\n",rfs->nparams);
  for(n=0; n < rfs->nparams; n++)
    fprintf(fp," param %d  %lf\n",n,rfs->params[n]);
  fprintf(fp,"mean %lf\n",rfs->mean);
  fprintf(fp,"stddev %lf\n",rfs->stddev);
  fprintf(fp,"seed   %ld\n",rfs->seed);
  return(0);
}
/*-------------------------------------------------------------------*/
int RFnparams(RFS *rfs)
{
  int nparams = -1;
  switch(rfs->code){
  case RF_UNIFORM:  nparams=2; break;
  case RF_GAUSSIAN: nparams=2; break;
  case RF_Z:        nparams=0; break;
  case RF_T:        nparams=1; break;
  case RF_F:        nparams=2; break;
  case RF_CHI2:     nparams=1; break;
  }
  rfs->nparams = nparams;
  return(nparams);
}
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddev(RFS *rfs)
{
  if(!strcmp(rfs->name,"uniform"))
    return(RFexpectedMeanStddevUniform(rfs));
  if(!strcmp(rfs->name,"gaussian"))
    return(RFexpectedMeanStddevGaussian(rfs));
  if(!strcmp(rfs->name,"t"))
    return(RFexpectedMeanStddevt(rfs));
  if(!strcmp(rfs->name,"F"))
    return(RFexpectedMeanStddevt(rfs));
  if(!strcmp(rfs->name,"chi2"))
    return(RFexpectedMeanStddevChi2(rfs));
  printf("ERROR: RFexpectedMeanStddev(): field type %s unknown\n",rfs->name);
  return(1);
}
/*-------------------------------------------------------------------*/
int RFsynth(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f;
  double v,m;

  if(RFname2Code(rfs) == -1) return(1);

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = MRIgetVoxVal(binmask,c,r,s,0);
	  if(m < 0.5) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = RFdrawVal(rfs);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
/*-------------------------------------------------------------------*/
MRI *RFstat2P(MRI *rf, RFS *rfs, MRI *binmask, MRI *p)
{
  int c,r,s,f=0,m;
  double v,pval;

  if(RFname2Code(rfs) == -1) return(NULL);
  p = MRIclone(rf,p);

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = MRIgetVoxVal(rf,c,r,s,f);
	  pval = RFstat2PVal(rfs,v);
	  MRIsetVoxVal(p,c,r,s,f,pval);
	}
      }
    }
  }
  return(p);
}
/*-------------------------------------------------------------------*/
MRI *RFp2Stat(MRI *p, RFS *rfs, MRI *binmask, MRI *rf)
{
  int c,r,s,f,m;
  double v,pval;

  if(RFname2Code(rfs) == -1) return(NULL);
  rf = MRIclone(p,rf);

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  pval = MRIgetVoxVal(p,c,r,s,f);
	  v = RFp2StatVal(rfs,pval);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(rf);
}
/*--------------------------------------------------------------------------*/
MRI *RFstat2Stat(MRI *rfin, RFS *rfsin, RFS *rfsout, MRI *binmask, MRI *rfout)
{
  MRI *p=NULL;

  if(RFname2Code(rfsin)  == -1) return(NULL);
  if(RFname2Code(rfsout) == -1) return(NULL);

  p     = RFstat2P(rfin, rfsin, binmask, p);
  rfout = RFp2Stat(p, rfsout, binmask, rfout);
  MRIfree(&p);
  return(rfout);
}

/*-------------------------------------------------------------------*/
MRI *RFrescale(MRI *rf, RFS *rfs, MRI *binmask, MRI *rfout)
{
  int c,r,s,f,m;
  double v, gmean, gstddev, gmax;

  if(RFname2Code(rfs) == -1) return(NULL);
  RFexpectedMeanStddev(rfs); // expected 
  RFglobalStats(rf, binmask, &gmean, &gstddev, &gmax); //actual

  rfout = MRIclone(rf,rfout);

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = MRIgetVoxVal(rf,c,r,s,f);
	  v = (v - gmean)*(rfs->stddev/gstddev) + rfs->mean;
	  MRIsetVoxVal(rfout,c,r,s,f,v);
	}
      }
    }
  }
  return(rfout);
}


/*-------------------------------------------------------------------*/
int RFglobalStats(MRI *rf, MRI *binmask, double *gmean, double *gstddev, double *max)
{
  int c,r,s,f,m;
  double v;
  double sum, sumsq;
  long nv;

  nv = 0;
  sum = 0;
  sumsq = 0;
  *max = -1000000;
  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = MRIgetVoxVal(rf,c,r,s,f);
	  sum += v;
	  sumsq += (v*v);
	  nv++;
	  if(*max < v) *max = v;
	}
      }
    }
  }
  *gmean = sum/nv;
  *gstddev = sqrt(sumsq/nv - (*gmean)*(*gmean));

  return(0);
}

/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
double RFdrawVal(RFS *rfs)
{
  if(!strcmp(rfs->name,"uniform")){
    //params[0] = min
    //params[1] = max
    return(sc_ran_flat(rfs->rng,rfs->params[0],rfs->params[1]));
  }
  if(!strcmp(rfs->name,"gaussian")){
    //params[0] = mean
    //params[1] = std
    return(sc_ran_gaussian(rfs->rng,rfs->params[1]) + rfs->params[0]);
  }
  if(!strcmp(rfs->name,"z")){
    //nparams=0, gaussian with mean=0, std=1
    return(sc_ran_gaussian(rfs->rng,1));
  }
  if(!strcmp(rfs->name,"t")){
    //params[0] = dof
    return(sc_ran_tdist(rfs->rng,rfs->params[0]));
  }
  if(!strcmp(rfs->name,"F")){
    //params[0] = numerator dof (rows in C)
    //params[1] = denominator dof
    return(sc_ran_fdist(rfs->rng,rfs->params[0],rfs->params[1]));
  }
  if(!strcmp(rfs->name,"chi2")){
    //params[0] = dof
    return(sc_ran_chisq(rfs->rng,rfs->params[0]));
  }
  printf("ERROR: RFdrawVal(): field type %s unknown\n",rfs->name);
  return(10000000000.0);
#else
double RFdrawVal(RFS *rfs)
{
  if(!strcmp(rfs->name,"uniform")){
    //params[0] = min
    //params[1] = max
    return(gsl_ran_flat(rfs->rng,rfs->params[0],rfs->params[1]));
  }
  if(!strcmp(rfs->name,"gaussian")){
    //params[0] = mean
    //params[1] = std
    return(gsl_ran_gaussian(rfs->rng,rfs->params[1]) + rfs->params[0]);
  }
  if(!strcmp(rfs->name,"z")){
    //nparams=0, gaussian with mean=0, std=1
    return(gsl_ran_gaussian(rfs->rng,1));
  }
  if(!strcmp(rfs->name,"t")){
    //params[0] = dof
    return(gsl_ran_tdist(rfs->rng,rfs->params[0]));
  }
  if(!strcmp(rfs->name,"F")){
    //params[0] = numerator dof (rows in C)
    //params[1] = denominator dof
    return(gsl_ran_fdist(rfs->rng,rfs->params[0],rfs->params[1]));
  }
  if(!strcmp(rfs->name,"chi2")){
    //params[0] = dof
    return(gsl_ran_chisq(rfs->rng,rfs->params[0]));
  }
  printf("ERROR: RFdrawVal(): field type %s unknown\n",rfs->name);
  return(10000000000.0);
}
#endif
/*-------------------------------------------------------------------
  RFtestVal() - returns the probability of getting a value of v or
  greater from a random draw from the given distribution. Note:
  this is a one-sided test (where sidedness makes sense).
  -------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
double RFstat2PVal(RFS *rfs, double stat)
{
  if(!strcmp(rfs->name,"uniform")){
    //params[0] = min
    //params[1] = max
    return(sc_cdf_flat_Q(stat, rfs->params[0], rfs->params[1]));
  }
  if(!strcmp(rfs->name,"gaussian")){
    //params[0] = mean
    //params[1] = std
    return(sc_cdf_gaussian_Q(stat-rfs->params[0],rfs->params[1]));
  }
  if(!strcmp(rfs->name,"z")){
    //nparams=0, gaussian with mean=0, std=1
    return(sc_cdf_gaussian_Q(stat,1));
  }
  if(!strcmp(rfs->name,"t")){
    //params[0] = dof
    return(sc_cdf_tdist_Q(stat,rfs->params[0]));
  }
  if(!strcmp(rfs->name,"F")){
    //params[0] = numerator dof (rows in C)
    //params[1] = denominator dof
    return(sc_cdf_fdist_Q(stat,rfs->params[0],rfs->params[1]));
  }
  if(!strcmp(rfs->name,"chi2")){
    //params[0] = dof
    return(sc_cdf_chisq_Q(stat,rfs->params[0]));
  }
  printf("ERROR: RFstat2PVal(): field type %s unknown\n",rfs->name);
  return(10000000000.0);
}
#else
double RFstat2PVal(RFS *rfs, double stat)
{
  if(!strcmp(rfs->name,"uniform")){
    //params[0] = min
    //params[1] = max
    return(gsl_cdf_flat_Q(stat, rfs->params[0], rfs->params[1]));
  }
  if(!strcmp(rfs->name,"gaussian")){
    //params[0] = mean
    //params[1] = std
    return(gsl_cdf_gaussian_Q(stat-rfs->params[0],rfs->params[1]));
  }
  if(!strcmp(rfs->name,"z")){
    //nparams=0, gaussian with mean=0, std=1
    return(gsl_cdf_gaussian_Q(stat,1));
  }
  if(!strcmp(rfs->name,"t")){
    //params[0] = dof
    return(gsl_cdf_tdist_Q(stat,rfs->params[0]));
  }
  if(!strcmp(rfs->name,"F")){
    //params[0] = numerator dof (rows in C)
    //params[1] = denominator dof
    return(gsl_cdf_fdist_Q(stat,rfs->params[0],rfs->params[1]));
  }
  if(!strcmp(rfs->name,"chi2")){
    //params[0] = dof
    return(gsl_cdf_chisq_Q(stat,rfs->params[0]));
  }
  printf("ERROR: RFstat2PVal(): field type %s unknown\n",rfs->name);
  return(10000000000.0);
}
#endif
/*-------------------------------------------------------------------
  RFp2StatVal() -
  -------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
double RFp2StatVal(RFS *rfs, double p)
{
  if(!strcmp(rfs->name,"uniform")){
    //params[0] = min
    //params[1] = max
    return(sc_cdf_flat_Qinv(p, rfs->params[0], rfs->params[1]));
  }
  if(!strcmp(rfs->name,"gaussian")){
    //params[0] = mean
    //params[1] = std
    return(sc_cdf_gaussian_Qinv(p,rfs->params[1])+rfs->params[0]);
  }
  if(!strcmp(rfs->name,"z")){
    //nparams=0, gaussian with mean=0, std=1
    return(sc_cdf_gaussian_Qinv(p,1));
  }
  if(!strcmp(rfs->name,"t")){
    //params[0] = dof
    return(sc_cdf_tdist_Qinv(p,rfs->params[0]));
  }
  if(!strcmp(rfs->name,"F")){
    //params[0] = numerator dof (rows in C)
    //params[1] = denominator dof
    // GSL does not have the Qinv for F :<.
    printf("ERROR: RFp2StatVal(): cannot invert RF field\n");
    return(10000000000.0);
  }
  if(!strcmp(rfs->name,"chi2")){
    //params[0] = dof
    return(sc_cdf_chisq_Qinv(p,rfs->params[0]));
  }
  printf("ERROR: RFp2v(): field type %s unknown\n",rfs->name);
  return(10000000000.0);
}
#else
double RFp2StatVal(RFS *rfs, double p)
{
  if(!strcmp(rfs->name,"uniform")){
    //params[0] = min
    //params[1] = max
    return(gsl_cdf_flat_Qinv(p, rfs->params[0], rfs->params[1]));
  }
  if(!strcmp(rfs->name,"gaussian")){
    //params[0] = mean
    //params[1] = std
    return(gsl_cdf_gaussian_Qinv(p,rfs->params[1])+rfs->params[0]);
  }
  if(!strcmp(rfs->name,"z")){
    //nparams=0, gaussian with mean=0, std=1
    return(gsl_cdf_gaussian_Qinv(p,1));
  }
  if(!strcmp(rfs->name,"t")){
    //params[0] = dof
    return(gsl_cdf_tdist_Qinv(p,rfs->params[0]));
  }
  if(!strcmp(rfs->name,"F")){
    //params[0] = numerator dof (rows in C)
    //params[1] = denominator dof
    // GSL does not have the Qinv for F :<.
    printf("ERROR: RFp2StatVal(): cannot invert RF field\n");
    return(10000000000.0);
  }
  if(!strcmp(rfs->name,"chi2")){
    //params[0] = dof
    return(gsl_cdf_chisq_Qinv(p,rfs->params[0]));
  }
  printf("ERROR: RFp2v(): field type %s unknown\n",rfs->name);
  return(10000000000.0);
}
#endif
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevUniform(RFS *rfs)
{
  double min,max,d;
  min = rfs->params[0];
  max  = rfs->params[1];
  d = max-min;
  rfs->mean = d/2;
  rfs->stddev = sqrt((d*d)/12.0);
  return(0);
}
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevGaussian(RFS *rfs)
{
  rfs->mean = rfs->params[0];
  rfs->stddev  = rfs->params[1];
  return(0);
}
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevt(RFS *rfs)
{
  double dof;
  dof = rfs->params[0];
  rfs->mean = 0;
  rfs->stddev  = sqrt(dof/(dof-2));
  return(0);
}
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevF(RFS *rfs)
{
  double ndof,ddof;
  ndof = rfs->params[0]; // numerator dof (rows in C)
  ddof = rfs->params[1]; // dof
  rfs->mean = ddof/(ddof-2);
  rfs->stddev  = 2*(ddof*ddof)*(ndof+ddof-2)/(ndof*((ddof-2)*(ddof-2))*(ddof-4));
  return(0);
}
/*-------------------------------------------------------------------*/
int RFexpectedMeanStddevChi2(RFS *rfs)
{
  double dof;
  dof = rfs->params[0];
  rfs->mean = dof;
  rfs->stddev = sqrt(2.0*dof);
  return(0);
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
  if(ar1 <= 0.0) return(0.0);
  gstd = d/sqrt(-4*log(ar1));
  return(gstd);
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
  gstd = RFar1ToGStd(ar1,d);
  fwhm = gstd*sqrt(log(256.0));
  return(fwhm);
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

#if 0

/*-------------------------------------------------------------------*/
int RFsynth(MRI *rf, RFS *rfs, MRI *binmask)
{
  if(!strcmp(rfs->name,"uniform"))
    return(RFsynthUniform(rf,rfs,binmask));
  if(!strcmp(rfs->name,"gaussian"))
    return(RFsynthGaussian(rf,rfs,binmask));
  if(!strcmp(rfs->name,"t"))
    return(RFsyntht(rf,rfs,binmask));
  if(!strcmp(rfs->name,"F"))
    return(RFsyntht(rf,rfs,binmask));
  printf("ERROR: RFsynth(): field type %s unknown\n",rfs->name);
  return(1);
}
/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
int RFsynthUniform(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f=0,m;
  double v,min,max;

  min = rfs->params[0];
  max = rfs->params[1];

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = sc_ran_flat(rfs->rng,min,max);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#else
int RFsynthUniform(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f=0,m;
  double v,min,max;

  min = rfs->params[0];
  max = rfs->params[1];

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = gsl_ran_flat(rfs->rng,min,max);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#endif

/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
int RFsynthGaussian(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f=0,m;
  double v,mean,std;

  mean = rfs->params[0];
  std  = rfs->params[1];

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = sc_ran_gaussian(rfs->rng,std) + mean;
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#else
int RFsynthGaussian(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f=0,m;
  double v,mean,std;

  mean = rfs->params[0];
  std  = rfs->params[1];

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = gsl_ran_gaussian(rfs->rng,std) + mean;
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#endif
/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
int RFsyntht(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f=0,m;
  double v,dof;

  dof = rfs->params[0];

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = sc_ran_tdist(rfs->rng,dof);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#else
int RFsyntht(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f=0,m;
  double v,dof;

  dof = rfs->params[0];

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = gsl_ran_tdist(rfs->rng,dof);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#endif
/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
int RFsynthF(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f,m;
  double v,ndof,ddof;

  ndof = rfs->params[0]; // numerator dof (rows in C)
  ddof = rfs->params[1]; // dof

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = sc_ran_fdist(rfs->rng,ndof,ddof);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#else
int RFsynthF(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f,m;
  double v,ndof,ddof;

  ndof = rfs->params[0]; // numerator dof (rows in C)
  ddof = rfs->params[1]; // dof

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = gsl_ran_fdist(rfs->rng,ndof,ddof);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#endif
/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
RFtestGaussian(MRI *rf, MRI *binmask, MRI *p)
{
  int c,r,s,f=0,m;
  double v,pval,mean,std;

  mean = rfs->params[0];
  std  = rfs->params[1];

  p = MRIclone(rf,p);

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = MRIgetVoxVal(rf,c,r,s,f);
	  pval = sc_cdf_gaussian_Q(v-mean,std);
	  MRIsetVoxVal(p,c,r,s,f,pval);
	}
      }
    }
  }
  return(p);
}
#else
RFtestGaussian(MRI *rf, MRI *binmask, MRI *p)
{
  int c,r,s,f=0,m;
  double v,pval,mean,std;

  mean = rfs->params[0];
  std  = rfs->params[1];

  p = MRIclone(rf,p);

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = MRIgetVoxVal(rf,c,r,s,f);
	  pval = gsl_cdf_gaussian_Q(v-mean,std);
	  MRIsetVoxVal(p,c,r,s,f,pval);
	}
      }
    }
  }
  return(p);
}
#endif
/*-------------------------------------------------------------------*/
#if USE_SC_GSL_REPLACEMENT
int RFtestt(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f,m;
  double v,pval,dof;

  dof = rfs->params[0];

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = sc_ran_tdist(rfs->rng,dof);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#else
int RFtestt(MRI *rf, RFS *rfs, MRI *binmask)
{
  int c,r,s,f,m;
  double v,pval,dof;

  dof = rfs->params[0];

  for(c=0; c < rf->width; c++){
    for(r=0; r < rf->height; r++){
      for(s=0; s < rf->depth; s++){
	if(binmask != NULL){
	  m = (int)MRIgetVoxVal(binmask,c,r,s,0);
	  if(!m) continue;
	}
	for(f=0; f < rf->nframes; f++){
	  v = gsl_ran_tdist(rfs->rng,dof);
	  MRIsetVoxVal(rf,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}
#endif

#endif
