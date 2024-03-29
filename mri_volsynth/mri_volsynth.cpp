/**
 * @brief Synthesizes data with various statistical properties.
 *
 */
/*
 * Original Author: Douglas N. Greve
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <ctype.h>

#include "version.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mri_identify.h"
#include "matrix.h"
#include "mri.h"
#include "mri2.h"
#include "fmriutils.h"
#include "MRIio_old.h"
#include "randomfields.h"
#include "mri_circulars.h"
#include "ctrpoints.h"
double round(double);

MRI *fMRIsqrt(MRI *mri, MRI *mrisqrt);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
static int  nth_is_arg(int nargc, char **argv, int nth);
static int checkfmt(char *fmt);
static int getfmtid(char *fname);
static int  isflag(char *flag);
//static int  stringmatch(char *str1, char *str2);

int main(int argc, char *argv[]) ;


const char *Progname = NULL;

int debug = 0;

char *volid = NULL;
char *vol_type;
int   volfmtid;
char *volfmt = NULL;

char *tempid = NULL; // Template
char *temp_type;
int   tempfmtid;
char *tempfmt = NULL;

int dim[4];
float res[4];
float cras[4];
float p0[4];
int usep0 = 0;
float cdircos[3], rdircos[3], sdircos[3];
const char *pdfname = "gaussian";
char *precision=NULL; 
int mritype; // precision
MRI *mri, *mrism, *mritemp, *mri2;
long seed = -1; /* < 0 for auto */
char *seedfile = NULL;
float fwhm = 0, gstd = 0, gmnnorm = 1;
int nframes = -1;
double TR = -1;
int delta_crsf[4];
int delta_crsf_speced = 0;
double delta_value = 1, delta_off_value = 0;
double gausmean=0, gausstd=1;
RFS *rfs;
int rescale = 0;
int numdof =  2;
int dendof = 20;
int AddOffset=0;
MRI *offset;
int OffsetFrame=0;
MRI *mask;
char *sum2file = NULL;
int NoOutput = 0;
MRI_REGION boundingbox;
double ValueA = 1;
double ValueB = 0;
double voxradius = -1;
double mmradius = -1;

int UseFFT = 0;
int SpikeTP = -1;
int DoCurv = 0;
char *subject=NULL, *hemi=NULL;
MRIS *surf;
int resSpeced=0,dimSpeced=0;
int NewVoxSizeSpeced=0;
int DoHSC=0; // force noise to be heteroscedastic
double HSCMin=0, HSCMax=0;
int DoTNorm=0;
int DoAbs=0;
MRI *fMRIhsynth(MRI *res, MRI *mask, int DoTNorm);
MPoint *ctrpoints=NULL, *crsctrpoints=NULL;
int nctrpoints=0, CPUseRealRAS;
int cgridspace=8, rgridspace=8, sgridspace=2;
int spherecenter[3], spherecenterset = 0;
double cubeedgemm = -1;
char *colortablefile = NULL;

/*---------------------------------------------------------------*/
int main(int argc, char **argv)
{
  int c,r,s,f,n,err=0;
  double val,rval;
  FILE *fp;
  MRI *mritmp;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  /* assign default geometry */
  cdircos[0] = 1.0;
  cdircos[1] = 0.0;
  cdircos[2] = 0.0;
  rdircos[0] = 0.0;
  rdircos[1] = 1.0;
  rdircos[2] = 0.0;
  sdircos[0] = 0.0;
  sdircos[1] = 0.0;
  sdircos[2] = 1.0;
  res[0] = 1.0;
  res[1] = 1.0;
  res[2] = 1.0;
  cras[0] = 0.0;
  cras[1] = 0.0;
  cras[2] = 0.0;
  res[3] = 2.0; /* TR */

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if(tempid != NULL) {
    printf("INFO: reading template header\n");
    if(! DoCurv) mritemp = MRIreadHeader(tempid,tempfmtid);
    else         mritemp = MRIread(tempid);
    if (mritemp == NULL) {
      printf("ERROR: reading %s header\n",tempid);
      exit(1);
    }
    if(NewVoxSizeSpeced){
      dim[0] = round(mritemp->width*mritemp->xsize/res[0]);
      dim[1] = round(mritemp->height*mritemp->ysize/res[1]);
      dim[2] = round(mritemp->depth*mritemp->zsize/res[2]);
      dim[3] = mritemp->nframes;
      res[3] = mritemp->tr;
      dimSpeced = 1;
    }
    if(dimSpeced){
      mritmp = MRIallocSequence(dim[0],dim[1],dim[2],MRI_FLOAT,dim[3]);
      MRIcopyHeader(mritemp,mritmp);
      MRIfree(&mritemp);
      mritemp = mritmp;
    }
    if(resSpeced){
      mritemp->xsize = res[0];
      mritemp->ysize = res[1];
      mritemp->zsize = res[2];
      mritemp->tr    = res[3];
    }
    else {
      res[0] = mritemp->xsize;
      res[1] = mritemp->ysize;
      res[2] = mritemp->zsize;
    }

    dim[0] = mritemp->width;
    dim[1] = mritemp->height;
    dim[2] = mritemp->depth;
    if (nframes > 0) dim[3] = nframes;
    else             dim[3] = mritemp->nframes;
    mritemp->nframes = dim[3];
  }

  if(mritemp) {
    if(SpikeTP >= mritemp->nframes){
      printf("ERROR: SpikeTP = %d >= mritemp->nframes = %d\n",
             SpikeTP,mritemp->nframes);
      exit(1);
    }
  }

  if(nctrpoints > 0){
    printf("Converting control points to voxel space\n");
    crsctrpoints = ControlPoints2Vox(ctrpoints, nctrpoints, CPUseRealRAS, mritemp);
  }

  printf("Synthesizing\n");
  srand48(seed);
  if (strcmp(pdfname,"gaussian")==0)
    mri = MRIrandn(dim[0], dim[1], dim[2], dim[3], gausmean, gausstd, NULL);
  else if (strcmp(pdfname,"uniform")==0)
    mri = MRIdrand48(dim[0], dim[1], dim[2], dim[3], 0, 1, NULL);
  else if (strcmp(pdfname,"const")==0)
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], ValueA, NULL);
  else if (strcmp(pdfname,"sphere")==0) {
    if(mmradius > 0)  voxradius = mmradius/((res[0]+res[1]+res[2])/3);
    if(voxradius < 0) voxradius = sqrt( pow(dim[0]/2.0,2)+pow(dim[1]/2.0,2)+pow(dim[2]/2.0,2) )/2.0;
    printf("voxradius = %lf\n",voxradius);
    if(!spherecenterset){
      spherecenter[0] = dim[0]/2.0;
      spherecenter[1] = dim[1]/2.0;
      spherecenter[2] = dim[2]/2.0;
    }
    mri = MRIsphereMask(dim[0], dim[1], dim[2], dim[3],
                        spherecenter[0],spherecenter[1],spherecenter[2],
                        voxradius, ValueA, NULL);
  } else if (strcmp(pdfname,"delta")==0) {
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], delta_off_value, NULL);
    if (delta_crsf_speced == 0) {
      delta_crsf[0] = dim[0]/2;
      delta_crsf[1] = dim[1]/2;
      delta_crsf[2] = dim[2]/2;
      delta_crsf[3] = dim[3]/2;
    }
    printf("delta set to %g at %d %d %d %d\n",delta_value,delta_crsf[0],
           delta_crsf[1],delta_crsf[2],delta_crsf[3]);
    MRIFseq_vox(mri,
                delta_crsf[0],
                delta_crsf[1],
                delta_crsf[2],
                delta_crsf[3]) = delta_value;
  } else if (strcmp(pdfname,"chi2")==0) {
    rfs = RFspecInit(seed,NULL);
    rfs->name = strcpyalloc("chi2");
    rfs->params[0] = dendof;
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    printf("Synthesizing chi2 with dof=%d\n",dendof);
    RFsynth(mri,rfs,NULL);
  } else if (strcmp(pdfname,"z")==0) {
    printf("Synthesizing z \n");
    rfs = RFspecInit(seed,NULL);
    rfs->name = strcpyalloc("gaussian");
    rfs->params[0] = 0; // mean
    rfs->params[1] = 1; // std
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    RFsynth(mri,rfs,NULL);
  } else if (strcmp(pdfname,"t")==0) {
    printf("Synthesizing t with dof=%d\n",dendof);
    rfs = RFspecInit(seed,NULL);
    rfs->name = strcpyalloc("t");
    rfs->params[0] = dendof;
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    RFsynth(mri,rfs,NULL);
  } else if (strcmp(pdfname,"tr")==0) {
    printf("Synthesizing t with dof=%d as ratio of z/sqrt(chi2)\n",dendof);
    rfs = RFspecInit(seed,NULL);
    // numerator
    rfs->name = strcpyalloc("gaussian");
    rfs->params[0] = 0; // mean
    rfs->params[1] = 1; // std
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    RFsynth(mri,rfs,NULL);
    // denominator
    rfs->name = strcpyalloc("chi2");
    rfs->params[0] = dendof;
    mri2 = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    RFsynth(mri2,rfs,NULL);
    fMRIsqrt(mri2,mri2); // sqrt of chi2
    mri = MRIdivide(mri,mri2,mri);
    MRIscalarMul(mri, mri, sqrt(dendof)) ;
    MRIfree(&mri2);
  } else if (strcmp(pdfname,"F")==0) {
    printf("Synthesizing F with num=%d den=%d\n",numdof,dendof);
    rfs = RFspecInit(seed,NULL);
    rfs->name = strcpyalloc("F");
    rfs->params[0] = numdof;
    rfs->params[1] = dendof;
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    RFsynth(mri,rfs,NULL);
  } else if (strcmp(pdfname,"Fr")==0) {
    printf("Synthesizing F with num=%d den=%d as ratio of two chi2\n",
           numdof,dendof);
    rfs = RFspecInit(seed,NULL);
    rfs->name = strcpyalloc("chi2");
    // numerator
    rfs->params[0] = numdof;
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    RFsynth(mri,rfs,NULL);
    // denominator
    rfs->params[0] = dendof;
    mri2 = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    RFsynth(mri2,rfs,NULL);
    mri = MRIdivide(mri,mri2,mri);
    MRIscalarMul(mri, mri, (double)dendof/numdof) ;
    MRIfree(&mri2);
  } else if (strcmp(pdfname,"voxcrs")==0) {
    // three frames. 1st=col, 2nd=row, 3rd=slice
    printf("Filling with vox CRS\n");
    mri = MRIconst(dim[0], dim[1], dim[2], 3, 0, NULL);
    for(c=0; c < mri->width; c ++){
      for(r=0; r < mri->height; r ++){
        for(s=0; s < mri->depth; s ++){
          MRIsetVoxVal(mri,c,r,s,0,c);
          MRIsetVoxVal(mri,c,r,s,1,r);
          MRIsetVoxVal(mri,c,r,s,2,s);
        }
      }
    }
  } else if (strcmp(pdfname,"boundingbox")==0) {
    printf("Setting bounding box \n");
    if(mritemp == NULL)
      mritemp = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    mri = MRIsetBoundingBox(mritemp,&boundingbox,ValueA,ValueB);
    if(!mri) exit(1);
  } 
  else if (strcmp(pdfname,"checker")==0) {
    printf("Checker \n");
    mri=MRIchecker(mritemp,NULL);
    if(!mri) exit(1);
  } 
  else if (strcmp(pdfname,"grid")==0) {
    printf("Grid %d %d %d\n",cgridspace,rgridspace,sgridspace);
    if(mritemp == NULL){
      mritemp = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
      mritemp->xsize = res[0];
      mritemp->ysize = res[1];
      mritemp->zsize = res[2];
      mritemp->tr = res[3];
    }
    mri=MRIgrid(mritemp,cgridspace,rgridspace,sgridspace,ValueA,NULL);
    if(!mri) exit(1);
  } 
  else if (strcmp(pdfname,"cube")==0) {
    int c0 = round(dim[0]/2.0);
    int r0 = round(dim[1]/2.0);
    int s0 = round(dim[2]/2.0);
    int c1 = round((dim[0] - cubeedgemm/res[0])/2.0);
    int c2 = round((dim[0] + cubeedgemm/res[0])/2.0);
    int r1 = round((dim[1] - cubeedgemm/res[1])/2.0);
    int r2 = round((dim[1] + cubeedgemm/res[1])/2.0);
    int s1 = round((dim[2] - cubeedgemm/res[2])/2.0);
    int s2 = round((dim[2] + cubeedgemm/res[2])/2.0);
    printf("Cube %lf (%d %d %d) %d %d %d %d %d %d\n",cubeedgemm,c0,r0,s0,c1,c2,r1,r2,s1,s2);
    printf("dim %d %d %d, res = %lf %lf %lf\n",dim[0],dim[1],dim[2],res[0],res[1],res[2]);
    if(c1 <= 0 || c1 >= dim[0]) exit(1);
    if(c2 <= 0 || c2 >= dim[0]) exit(1);
    if(r1 <= 0 || r1 >= dim[1]) exit(1);
    if(r2 <= 0 || r2 >= dim[1]) exit(1);
    if(s1 <= 0 || s1 >= dim[2]) exit(1);
    if(s2 <= 0 || s2 >= dim[2]) exit(1);
    if(mritemp == NULL){
      mritemp = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
      mritemp->xsize = res[0];
      mritemp->ysize = res[1];
      mritemp->zsize = res[2];
      mritemp->tr = res[3];
    }
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    if(!mri) exit(1);
    mri->xsize = res[0];
    mri->ysize = res[1];
    mri->zsize = res[2];
    mri->tr = res[3];
    MRIsetVoxVal(mri,c0,r0,s0,0,ValueA);
    MRIsetVoxVal(mri,c1,r1,s1,0,ValueA);
    MRIsetVoxVal(mri,c1,r1,s2,0,ValueA);
    MRIsetVoxVal(mri,c1,r2,s1,0,ValueA);
    MRIsetVoxVal(mri,c1,r2,s2,0,ValueA);
    MRIsetVoxVal(mri,c2,r1,s1,0,ValueA);
    MRIsetVoxVal(mri,c2,r1,s2,0,ValueA);
    MRIsetVoxVal(mri,c2,r2,s1,0,ValueA);
    MRIsetVoxVal(mri,c2,r2,s2,0,ValueA);
  } 
  else if (strcmp(pdfname,"sliceno")==0) {
    printf("SliceNo \n");
    if(mritemp == NULL){
      printf("ERROR: need --temp with sliceno\n");
      exit(1);
    }
    mri=MRIsliceNo(mritemp,NULL);
    if(!mri) exit(1);
  } 
  else if (strcmp(pdfname,"indexno")==0) {
    printf("IndexNo \n");
    if(mritemp == NULL){
      printf("ERROR: need --temp with indexno\n");
      exit(1);
    }
    mri=MRIindexNo(mritemp,NULL);
    if(!mri) exit(1);
  } 
  else if (strcmp(pdfname,"crs")==0) {
    printf("CRS \n");
    if(mritemp == NULL){
      printf("ERROR: need --temp with crs\n");
      exit(1);
    }
    mri=MRIcrs(mritemp,NULL);
    if(!mri) exit(1);
  } 
  else if (strcmp(pdfname,"cp")==0) {
    printf("Synthesizing control points volume \n");
    mri = MRIconst(dim[0], dim[1], dim[2], dim[3], 0, NULL);
    for(n=0; n < nctrpoints; n++){
      c = round(crsctrpoints[n].x);
      r = round(crsctrpoints[n].y);
      s = round(crsctrpoints[n].z);
      MRIFseq_vox(mri,c,r,s,0) += 1;
    }
  }
  else {
    printf("ERROR: pdf %s unrecognized, must be gaussian, uniform,\n"
	   "const, delta, checker\n", pdfname);
    exit(1);
  }
  if (tempid != NULL) {
    MRIcopyHeader(mritemp,mri);
    mri->type = MRI_FLOAT;
    // Override
    if(nframes > 0) mri->nframes = nframes;
    if(TR > 0) mri->tr = TR;
  } 
  else {
    if(mri == NULL) {
      usage_exit();
    }
    mri->xsize = res[0];
    mri->ysize = res[1];
    mri->zsize = res[2];
    mri->tr    = res[3];
    mri->x_r = cdircos[0];
    mri->x_a = cdircos[1];
    mri->x_s = cdircos[2];
    mri->y_r = rdircos[0];
    mri->y_a = rdircos[1];
    mri->y_s = rdircos[2];
    mri->z_r = sdircos[0];
    mri->z_a = sdircos[1];
    mri->z_s = sdircos[2];
    if(!usep0){
      mri->c_r = cras[0];
      mri->c_a = cras[1];
      mri->c_s = cras[2];
    } 
    else MRIp0ToCRAS(mri, p0[0], p0[1], p0[2]);
  }

  if (gstd > 0) {
    if(!UseFFT){
      printf("Smoothing\n");
      MRIgaussianSmooth(mri, gstd, gmnnorm, mri); /* gmnnorm = 1 = normalize */
    }
    else {
      printf("Smoothing with FFT \n");
      mri2 = MRIcopy(mri,NULL);
      mri = MRI_fft_gaussian(mri2, mri,
                             gstd, gmnnorm); /* gmnnorm = 1 = normalize */
    }
    if (rescale) {
      printf("Rescaling\n");
      if (strcmp(pdfname,"z")==0)     RFrescale(mri,rfs,NULL,mri);
      if (strcmp(pdfname,"chi2")==0)  RFrescale(mri,rfs,NULL,mri);
      if (strcmp(pdfname,"t")==0)     RFrescale(mri,rfs,NULL,mri);
      if (strcmp(pdfname,"tr")==0)    RFrescale(mri,rfs,NULL,mri);
      if (strcmp(pdfname,"F")==0)     RFrescale(mri,rfs,NULL,mri);
      if (strcmp(pdfname,"Fr")==0)    RFrescale(mri,rfs,NULL,mri);
    }
  }

  if(DoHSC){
    // This multiplies each frame by a random number
    // between HSCMin HSCMax to simulate heteroscedastisity
    printf("Applying HSC %lf %lf\n",HSCMin,HSCMax);
    for(f=0; f < mri->nframes; f++){
      rval = (HSCMax-HSCMin)*drand48() + HSCMin;
      if(debug) printf("%3d %lf\n",f,rval);
      for(c=0; c < mri->width; c ++){
	for(r=0; r < mri->height; r ++){
	  for(s=0; s < mri->depth; s ++){
	    val = MRIgetVoxVal(mri,c,r,s,f);
	    MRIsetVoxVal(mri,c,r,s,f,rval*val);
	  }
        }
      }
    }
  }

  if(AddOffset) {
    printf("Adding offset\n");
    offset = MRIread(tempid);
    if(offset == NULL) exit(1);
    if(OffsetFrame == -1) OffsetFrame = nint(offset->nframes/2);
    printf("Offset frame %d\n",OffsetFrame);
    mritmp = fMRIframe(offset, OffsetFrame, NULL);
    if(mritmp == NULL) exit(1);
    MRIfree(&offset);
    offset = mritmp;
    fMRIaddOffset(mri, offset, NULL, mri);
  }

  if(SpikeTP > 0){
    printf("Spiking time point %d\n",SpikeTP);
    for(c=0; c < mri->width; c ++){
      for(r=0; r < mri->height; r ++){
        for(s=0; s < mri->depth; s ++){
          MRIsetVoxVal(mri,c,r,s,SpikeTP,1e9);
        }
      }
    }
  }

  if(DoAbs){
    printf("Computing absolute value\n");
    MRIabs(mri,mri);
  }
  if(precision != NULL){
    printf("Changing precision to %s (no rescale)\n",precision);
    MRI *mri2 = MRISeqchangeType(mri, mritype, 0.0, 0.999, 1);
    if(mri2 == NULL){
      printf("ERROR: MRISeqchangeType\n");
      exit(1);
    }
    MRIfree(&mri);
    mri = mri2;
  }

  if(colortablefile){
    printf("Embedding ctab from %s\n",colortablefile);
    mri->ct = CTABreadASCII(colortablefile);
    if(mri->ct == NULL) exit(1);
  }

  if(!NoOutput){
    printf("Saving\n");
    if(!DoCurv) err = MRIwriteAnyFormat(mri,volid,volfmt,-1,NULL);
    else {
      printf("Saving in curv format\n");
      MRIScopyMRI(surf, mri, 0, "curv");
      err = MRISwriteCurvature(surf,volid);
    }
    if(err) exit(1);
  }

  if(sum2file){
    val = MRIsum2All(mri);
    fp = fopen(sum2file,"w");
    if(fp == NULL){
      printf("ERROR: opening %s\n",sum2file);
      exit(1);
    }
    printf("sum2all: %20.10lf\n",val);
    printf("vrf: %20.10lf\n",1/val);
    fprintf(fp,"%20.10lf\n",val);
  }

  return(0);
}
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  i, nargc , nargsused;
  char **pargv, *option ;
  char tmpstr[1000];

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))           print_help() ;
    else if (!strcasecmp(option, "--version"))   print_version() ;
    else if (!strcasecmp(option, "--debug"))     debug = 1;
    else if (!strcasecmp(option, "--abs"))       DoAbs = 1;
    else if (!strcasecmp(option, "--nogmnnorm")) gmnnorm = 0;
    else if (!strcasecmp(option, "--rescale"))   rescale=1;
    else if (!strcasecmp(option, "--norescale")) rescale=0;
    else if (!strcasecmp(option, "--no-output")) NoOutput = 1;
    else if (!strcasecmp(option, "--fft"))       UseFFT = 1;
    else if (!strcasecmp(option, "--offset"))    AddOffset=1;
    else if (!strcasecmp(option, "--offset-mid")){
      AddOffset = 1;
      OffsetFrame = -1;
    }
    else if (!strcmp(option, "--hsc")) {
      if (nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%lf",&HSCMin);
      sscanf(pargv[1],"%lf",&HSCMax);
      DoHSC = 1;
      nargsused = 2;
    }
    else if (!strcmp(option, "--sum2")) {
      if (nargc < 1) argnerr(option,1);
      sum2file = pargv[0];
      pdfname  = "delta";
      //NoOutput = 1;
      nframes  = 1;
      nargsused = 1;
    }
    else if (!strcmp(option, "--hsynth")) {
      // eres mask DoTnorm out
      if (nargc < 3) argnerr(option,3);
      mri = MRIread(pargv[0]);
      mask = MRIread(pargv[1]);
      sscanf(pargv[2],"%d",&DoTNorm);
      mri2 = fMRIhsynth(mri, mask, DoTNorm);
      MRIwrite(mri2,pargv[3]);
      exit(0);
      nargsused = 2;
    }
    else if (!strcmp(option, "--vol") || !strcmp(option, "--o")) {
      if (nargc < 1) argnerr(option,1);
      volid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        volfmt = pargv[1];
        nargsused ++;
        volfmtid = checkfmt(volfmt);
      } 
    }
    else if (!strcmp(option, "--temp") || !strcmp(option, "--template")) {
      if(DoCurv){
	printf("ERROR: cannot use --temp and --curv\n");
	exit(1);
      }
      if(nargc < 1) argnerr(option,1);
      tempid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        tempfmt = pargv[1];
        nargsused ++;
        tempfmtid = checkfmt(tempfmt);
      } else tempfmtid = getfmtid(tempid);
    } else if (!strcmp(option, "--curv")) {
      if(tempid != NULL){
	printf("ERROR: cannot use --temp and --curv\n");
	exit(1);
      }
      if(nargc < 2) argnerr(option,2);
      subject = pargv[0];
      hemi = pargv[1];
      DoCurv = 1;
      nargsused = 2;
      sprintf(tmpstr,"%s/%s/surf/%s.thickness",getenv("SUBJECTS_DIR"),subject,hemi);
      tempid = strcpyalloc(tmpstr);
    } 
    else if ( !strcmp(option, "--dim") ) {
      if (nargc < 4) argnerr(option,4);
      for (i=0;i<4;i++) sscanf(pargv[i],"%d",&dim[i]);
      nargsused = 4;
      dimSpeced = 1;
    } 
    else if ( !strcmp(option, "--dim-surf") ) {
      if(nargc < 1) argnerr(option,1);
      MRIS *surf = MRISread(pargv[0]);
      if(surf==NULL) exit(1);
      for (i=0;i<3;i++) dim[i] = 1;
      dim[0] = surf->nvertices;
      if(dim[3] == 0) dim[3] = 1;
      nargsused = 1;
      dimSpeced = 1;
    } 
    else if ( !strcmp(option, "--nframes") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nframes);
      dim[3] = nframes;
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--TR") || !strcmp(option, "--tr") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&TR);
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--res") ) {
      if (nargc < 4) argnerr(option,4);
      for (i=0;i<4;i++) sscanf(pargv[i],"%f",&res[i]);
      nargsused = 4;
      resSpeced = 1;
    } 
    else if ( !strcmp(option, "--vox-size") ) {
      if (nargc < 3) argnerr(option,3);
      for (i=0;i<3;i++) sscanf(pargv[i],"%f",&res[i]);
      nargsused = 4;
      resSpeced = 1;
      NewVoxSizeSpeced = 1;
    } 
    else if ( !strcmp(option, "--c_ras") ) {
      if (nargc < 3) argnerr(option,3);
      for (i=0;i<3;i++) sscanf(pargv[i],"%f",&cras[i]);
      nargsused = 3;
    } 
    else if ( !strcmp(option, "--p0") ) {
      if (nargc < 3) argnerr(option,3);
      for (i=0;i<3;i++) sscanf(pargv[i],"%f",&p0[i]);
      usep0 = 1;
      nargsused = 3;
    } 
    else if ( !strcmp(option, "--cdircos") ) {
      if (nargc < 3) argnerr(option,3);
      for (i=0;i<3;i++) sscanf(pargv[i],"%f",&cdircos[i]);
      nargsused = 3;
    } else if ( !strcmp(option, "--rdircos") ) {
      if (nargc < 3) argnerr(option,3);
      for (i=0;i<3;i++) sscanf(pargv[i],"%f",&rdircos[i]);
      nargsused = 3;
    } else if ( !strcmp(option, "--sdircos") ) {
      if (nargc < 3) argnerr(option,3);
      for (i=0;i<3;i++) sscanf(pargv[i],"%f",&sdircos[i]);
      nargsused = 3;
    } else if ( !strcmp(option, "--fwhm") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&fwhm);
      gstd = fwhm/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcmp(option, "--precision")) {
      if (nargc < 1) argnerr(option,1);
      precision = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--seed")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%ld",&seed);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--seedfile")) {
      if (nargc < 1) argnerr(option,1);
      seedfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--pdf")) {
      if (nargc < 1) argnerr(option,1);
      pdfname = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--bb")) {
      if (nargc < 6) argnerr(option,6);
      sscanf(pargv[0],"%d",&boundingbox.x);
      sscanf(pargv[1],"%d",&boundingbox.y);
      sscanf(pargv[2],"%d",&boundingbox.z);
      sscanf(pargv[3],"%d",&boundingbox.dx);
      sscanf(pargv[4],"%d",&boundingbox.dy);
      sscanf(pargv[5],"%d",&boundingbox.dz);
      pdfname = "boundingbox";
      nargsused = 6;
    }
    else if (!strcasecmp(option, "--checker"))
      pdfname = "checker";
    else if (!strcmp(option, "--dof-num")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&numdof);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--val-a")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&ValueA);
      nargsused = 1;
    } else if (!strcmp(option, "--val-b")) {
      if (nargc < 1) argnerr(option,1);
      if(*pargv[0] != '-' && !isdigit(*pargv[0])){
	printf("ERROR: --val-b must be a number\n");
	exit(1);
      }
      sscanf(pargv[0],"%lf",&ValueB);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--radius") || !strcmp(option, "--vox-radius")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&voxradius);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--mm-radius")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&mmradius);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--sphere-center")) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%d",&spherecenter[0]);
      sscanf(pargv[1],"%d",&spherecenter[1]);
      sscanf(pargv[2],"%d",&spherecenter[2]);
      spherecenterset = 1;
      nargsused = 3;
    } 
    else if (!strcmp(option, "--dof-den") || !strcmp(option, "--dof")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&dendof);
      nargsused = 1;
    } else if (!strcmp(option, "--gmean")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&gausmean);
      nargsused = 1;
    } else if (!strcmp(option, "--gstd")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&gausstd);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--grid")) {
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%d",&cgridspace);
      sscanf(pargv[1],"%d",&rgridspace);
      sscanf(pargv[2],"%d",&sgridspace);
      pdfname = "grid";
      nargsused = 3;
    } 
    else if (!strcmp(option, "--cube")) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&cubeedgemm);
      pdfname = "cube";
      nargsused = 1;
    } 
    else if (!strcmp(option, "--delta-crsf")) {
      if (nargc < 4) argnerr(option,4);
      sscanf(pargv[0],"%d",&delta_crsf[0]);
      sscanf(pargv[1],"%d",&delta_crsf[1]);
      sscanf(pargv[2],"%d",&delta_crsf[2]);
      sscanf(pargv[3],"%d",&delta_crsf[3]);
      delta_crsf_speced = 1;
      nargsused = 4;
    } else if (!strcmp(option, "--delta-val")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&delta_value);
      nargsused = 1;
    } else if (!strcmp(option, "--delta-val-off")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&delta_off_value);
      nargsused = 1;
    } else if (!strcmp(option, "--spike")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&SpikeTP);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--cp")) {
      if (nargc < 1) argnerr(option,1);
      printf("Reading control points from %s\n",pargv[0]);
      ctrpoints = MRIreadControlPoints(pargv[0], &nctrpoints, &CPUseRealRAS);
      printf("nctrpoints = %d, UseRealRAS=%d\n",nctrpoints, CPUseRealRAS);
      pdfname = "cp";
      nargsused = 1;
    } 
    else if (!strcmp(option, "--ctab")) {
      if (nargc < 1) argnerr(option,1);
      colortablefile = pargv[0];
      nargsused = 1;
    } 
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --o volid <fmt> : output volume path id and format\n");
  printf("\n");
  printf(" Get geometry from template\n");
  printf("   --template templateid : see also --curv\n");
  printf("   --nframes nframes : override template\n");
  printf("   --offset : use template as intensity offset\n");
  printf("   --offset-mid : use middle frame of template as intensity offset\n");
  printf("   --curv subject hemi : save output as curv (uses lh.thickness as template)\n");
  printf("\n");
  printf(" Specify geometry explicitly\n");
  printf("   --dim nc nr ns nf : required without template (overrides template)\n");
  printf("   --res dc dr ds df : voxel resolution (df is TR, use msec) (overrides template)\n");
  printf("   --vox-size dc dr ds : changes template voxel res AND dim\n");
  printf("   --tr TR : time between frames in msec \n");
  printf("   --cdircos x y z\n");
  printf("   --rdircos x y z\n");
  printf("   --sdircos x y z\n");
  printf("   --c_ras   c_r c_a c_s : 'center' voxel\n");
  printf("   --p0   p0r p0a p0s : first voxel\n");
  printf("   --precision precision : eg, float\n");
  printf("\n");
  printf(" Value distribution flags\n");
  printf("   --seed seed (default is time-based auto)\n");
  printf("   --seedfile fname : write seed value to this file\n");
  printf("   --pdf pdfname : <gaussian>, uniform, const, delta, \n");
  printf("      sphere, z, t, F, chi2, voxcrs, checker, sliceno, indexno\n");
  printf("      crs, grid\n");
  printf("   --bb c r s dc dr ds : bounding box (In=ValA, Out=ValB)\n");
  printf("   --gmean mean : use mean for gaussian (def is 0)\n");
  printf("   --gstd  std  : use std for gaussian standard dev (def is 1)\n");
  printf("   --delta-crsf col row slice frame : 0-based\n");
  printf("   --delta-val val : set delta value to val. Default is 1.\n");
  printf("   --delta-val-off offval : set delta background value to offval. "
         "Default is 0.\n");
  printf("   --grid dcol, drow, dslice\n");
  printf("   --dof dof : dof for t and chi2 \n");
  printf("   --dof-num numdof : numerator dof for F\n");
  printf("   --dof-den dendof : denomenator dof for F\n");
  printf("   --rescale : rescale z, t, F, or chi2 after smoothing\n");
  printf("   --val-a value : set ValA (default 1)\n");
  printf("   --val-b value : set ValB (default 0)\n");
  printf("   --vox-radius voxradius : radius (in voxels) for sphere\n");
  printf("   --mm-radius  radius : radius (in mm) for sphere, will be iso in voxels\n");
  printf("   --sphere-center col row slice\n");
  printf("   --hsc min max : multiply each frame by a random number bet min and max\n");
  printf("   --abs : compute absolute value\n");
  printf("   --cp control.dat : set control point voxels = 1 \n");
  printf("\n");
  printf(" Other arguments\n");
  printf("   --spike tp : set all values at time point tp to 1e9\n");
  printf("   --fwhm fwhm_mm : smooth by FWHM mm\n");
  printf("   --sum2 fname   : save sum vol^2 into fname (implies "
         "delta,nf=1,no-output)\n");
  printf("   --dim-surf : set dim to nvertices x 1 x 1 \n");
  printf("   --ctab colortable : embed ctab\n");
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("Synthesizes a volume with the given geometry and pdf. "
         "Default pdf \n");
  printf("is gaussian (mean 0, std 1). If uniform is chosen, then the min\n");
  printf("is 0 and the max is 1. If const is chosen, then all "
         "voxels are set\n");
  printf("to ValA. If delta, the middle voxel is set to 1, "
         "the rest to 0 unless\n");
  printf("the actual location is chosen with --delta-crsf. If sphere, then\n");
  printf("a spherical binary mask centered in the volume with radius half\n");
  printf("the field of view.\n");
  printf("Example: to create an empty volume (all voxels=0), named "
         "myvol.mgz,\nbased on parameters from an existing volume:\n");
  printf("  mri_volsynth --gstd 0 --vol myvol.mgz "
         "--temp some_existing_vol.mgz\n");
  printf("\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void) {
  struct timeval tv;
  FILE *fp;
  char tmpstr[1000];

  if(volid == NULL && !NoOutput) {
    printf("A volume path must be supplied unless --no-output\n");
    exit(1);
  }
  if (seed < 0) {
    gettimeofday(&tv, NULL);
    seed = tv.tv_sec + tv.tv_usec;
  }
  if (seedfile != NULL) {
    fp = fopen(seedfile,"w");
    if (fp == NULL) {
      printf("ERROR: cannot open seed file %s\n",seedfile);
      exit(1);
    }
    fprintf(fp,"%ld\n",seed);
    fclose(fp);
  }
  if(DoCurv && nframes > 1){
    printf("ERROR: cannot have more than 1 frame with curv output\n");
    exit(1);
  }
  if(DoCurv){
    sprintf(tmpstr,"%s/%s/surf/%s.orig",getenv("SUBJECTS_DIR"),subject,hemi);
    printf("Loading %s\n",tmpstr);
    surf = MRISread(tmpstr);
    if(!surf) exit(1);
  }
  if(!DoCurv) getfmtid(volid);

  if(precision != NULL){
    
    if(strcmp(StrLower(precision), "uchar") == 0)       mritype = MRI_UCHAR;
    else if(strcmp(StrLower(precision), "short") == 0)  mritype = MRI_SHORT;
    else if(strcmp(StrLower(precision), "int") == 0)    mritype = MRI_INT;
    else if(strcmp(StrLower(precision), "float") == 0)  mritype = MRI_FLOAT;
    else {
      printf("ERROR: precision %s not supported\n",precision);
      exit(1);
    }
  }
  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"volid    %s\n",volid);
  fprintf(fp,"volfmt   %s\n",volfmt);
  if (tempid == 0) {
    fprintf(fp,"dim    %3d %3d %3d %3d\n",
            dim[0],dim[1],dim[2],dim[3]);
    fprintf(fp,"res      %6.4f %6.4f %6.4f %6.4f\n",
            res[0],res[1],res[2],res[3]);
    fprintf(fp,"c_ras    %6.4f %6.4f %6.4f\n",
            cras[0], cras[1], cras[2]);
    fprintf(fp,"col   dircos  %6.4f %6.4f %6.4f\n",
            cdircos[0],cdircos[1],cdircos[2]);
    fprintf(fp,"row   dircos  %6.4f %6.4f %6.4f\n",
            rdircos[0],rdircos[1],rdircos[2]);
    fprintf(fp,"slice dircos  %6.4f %6.4f %6.4f\n",
            sdircos[0],sdircos[1],sdircos[2]);
  } else {
    fprintf(fp,"Using %s as template\n",tempid);
    if (nframes > 0) fprintf(fp,"nframes = %d\n",nframes);
  }
  fprintf(fp,"fwhm = %g, gstd  = %g\n",fwhm,gstd);
  //fprintf(fp,"precision %s\n",precision);
  fprintf(fp,"seed %ld\n",seed);
  fprintf(fp,"pdf   %s\n",pdfname);
  fprintf(fp,"SpikeTP %d\n",SpikeTP);
  fprintf(fp,"DoCurv %d\n",DoCurv);
  fprintf(fp,"DoAbs  %d\n",DoAbs);
  printf("Diagnostic Level %d\n",Gdiag_no);

  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);
  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth) {
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */
  /* check that there are enough args for nth to exist */
  if (nargc <= nth) return(0);
  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) return(0);
  return(1);
}
/*------------------------------------------------------------*/
static int getfmtid(char *fname) {
  int fmtid;
  fmtid = mri_identify(fname);
  if (fmtid == MRI_VOLUME_TYPE_UNKNOWN) {
    printf("ERROR: cannot determine format of %s\n",fname);
    exit(1);
  }
  return(fmtid);
}

/*------------------------------------------------------------*/
static int checkfmt(char *fmt) {
  int fmtid;

  if (fmt == NULL) return(MRI_VOLUME_TYPE_UNKNOWN);
  fmtid = string_to_type(fmt);
  if (fmtid == MRI_VOLUME_TYPE_UNKNOWN) {
    printf("ERROR: format string %s unrecognized\n",fmt);
    exit(1);
  }
  return(fmtid);
}
/*---------------------------------------------------------------*/
static int isflag(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}

/*---------------------------------------------------------------*/
MRI *fMRIsqrt(MRI *mri, MRI *mrisqrt) {
  int c,r,s,f;
  double val;

  if (mrisqrt == NULL) mrisqrt = MRIcopy(mri,NULL);
  if (mrisqrt == NULL) return(NULL);

  for (c=0; c < mri->width; c++) {
    for (r=0; r < mri->height; r++) {
      for (s=0; s < mri->depth; s++) {
        for (f=0; f < mri->nframes; f++) {
          val = MRIgetVoxVal(mri,c,r,s,f);
          if (val < 0.0) val = 0;
          else          val = sqrt(val);
          MRIsetVoxVal(mri,c,r,s,f,val);
        }
      }
    }
  }
  return(mrisqrt);
}

/*---------------------------------------------------------------*/
MRI *fMRIhsynth(MRI *res, MRI *mask, int DoTNorm)
{
  int c,r,s,f,nvox;
  double val;
  MRI *hsynth, *tvar=NULL;
  double *svar, svarsum, tstdvox;

  // Compute temporal variance at each voxel
  if(DoTNorm) tvar = fMRIcovariance(res, 0, -1, mask, NULL);

  // Compute spatial variance at each frame
  svar = (double *) calloc(res->nframes,sizeof(double));
  svarsum = 0;
  for (f=0; f < res->nframes; f++) {
    svar[f] = 0;
    nvox = 0;
    for (c=0; c < res->width; c++) {
      for (r=0; r < res->height; r++) {
	for (s=0; s < res->depth; s++) {
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) == 0) continue;
          val = MRIgetVoxVal(res,c,r,s,f);
	  if(DoTNorm){
	    tstdvox = sqrt(MRIgetVoxVal(tvar,c,r,s,0));
	    val /= tstdvox;
	  }
	  svar[f] += (val*val);
	  nvox ++;
        }
      }
    }
    svar[f] /= nvox;
    svarsum += svar[f];
  }
  for (f=0; f < res->nframes; f++) svar[f] /= (svarsum/res->nframes);
  for (f=0; f < res->nframes; f++) printf("%2d %g\n",f,svar[f]);

  // Synth noise that is both spatially and temporally white and gaussian
  hsynth = MRIrandn(res->width, res->height, res->depth, res->nframes, 0, 1, NULL);
  MRIcopyHeader(res,hsynth);

  // Scale by frame. Noise is still independent across 
  // space and time, but it is no longer homogeneous.
  for (f=0; f < res->nframes; f++) {
    for (c=0; c < res->width; c++) {
      for (r=0; r < res->height; r++) {
	for (s=0; s < res->depth; s++) {
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) == 0) {
	    MRIsetVoxVal(hsynth,c,r,s,f,0);
	    continue;
	  }
          val = MRIgetVoxVal(hsynth,c,r,s,f);
          val *= sqrt(svar[f]);
	  MRIsetVoxVal(hsynth,c,r,s,f,val);
        }
      }
    }
  }
  free(svar);
  if(DoTNorm) MRIfree(&tvar);
  return(hsynth);
}




