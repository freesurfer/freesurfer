/**
 * @file  mri_rbvpvc.c
 * @brief Peforms RBV Partial volume correction
 *
 * Implementation of Region-based Voxelwise (RBV) partial volume correction
 * as found in Thomas, et al, 2011, Eur J Nucl Med Mol Imaging, 38:1104-1119.
 * It also implements the Geometric Transfer Matrix (GTM) as it is needed by RBV.
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2013/11/10 21:06:30 $
 *    $Revision: 1.9 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

/*!
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/


// $Id: mri_rbvpvc.c,v 1.9 2013/11/10 21:06:30 greve Exp $

/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/

// Things to do:
// 1. Bounding box smoothing
// 2. FFT smoothing, sparse FFT
// 3. Subvoxel
// 4. Parallel
// 5. XtX symetric

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "timer.h"
#include "fmriutils.h"
#include "matfile.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_rbvpvc.c,v 1.9 2013/11/10 21:06:30 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *SrcVolFile=NULL,*SegVolFile=NULL,*MaskVolFile=NULL;
char *OutDir=NULL,*OutVolFile=NULL;
char *OutBetaFile=NULL,*OutXtXFile=NULL;
double psfFWHM=-1,cFWHM,rFWHM,sFWHM,cStd,rStd,sStd;
int GTMOnly=0;
int DoSegTest=0;
char tmpstr[5000];
int niterations = 0;
char *PVFFile=NULL, *SegTTypeFile=NULL;

MATRIX *MatrixGetDiag(MATRIX *M, VECTOR *d);
int VRFStats(MATRIX *iXtX, double *vrfmean, double *vrfmin, double *vrfmax);
MRI *IterativeRBV0PVC(MRI *src, MRI *seg, MRI *mask, int niters,
		      double cFWHM, double rFWHM, double sFWHM, 
		      MRI *dst);
MATRIX *BuildGTMPVF(MRI *seg, MATRIX *SegTType, MRI *pvf, MRI *mask, 
		    double cFWHM, double rFWHM, double sFWHM, 
		    MATRIX *X);
MATRIX *Seg2TissueType(MRI *seg);

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err,c,r,s,f,nmask;
  MRI *src, *seg, *mask=NULL,*rbv,*gtmsynth,*gtmsynthsm,*mritmp;
  int nsegs,msegs,nthseg,*segidlist0,*segidlist,segid;
  MATRIX *beta, *y, *X, *Xt, *XtX, *iXtX, *Xty;
  double val,XtXcond;
  double vrfmin,vrfmax,vrfmean;
  struct timeb  mytimer, timer;
  MRI *PVF=NULL;
  MATRIX *SegTType=NULL;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);

  // Convert FWHM to StdDev
  cStd = cFWHM/sqrt(log(256.0));
  rStd = rFWHM/sqrt(log(256.0));
  sStd = sFWHM/sqrt(log(256.0));
  printf("FWHM: %g %g %g\n",cFWHM,rFWHM,sFWHM);
  printf("Std:  %g %g %g\n",cStd,rStd,sStd);

  TimerStart(&timer);

  // Load data
  TimerStart(&mytimer);
  printf("Loading seg %s\n",SegVolFile);fflush(stdout);
  seg = MRIread(SegVolFile);
  if(seg==NULL) exit(1);

  if(DoSegTest == 0){
    printf("Loading input %s\n",SrcVolFile);fflush(stdout);
    src = MRIread(SrcVolFile);
    if(src==NULL) exit(1);
  } 
  else {
    printf("Using smoothed seg as input \n");
    src = MRIgaussianSmoothNI(seg, cStd, rStd, sStd, NULL);
  }
  err = MRIdimMismatch(seg, src, 0);
  if(err){
    printf("ERROR: seg and source dim mismatch %d\n",err);
    exit(1);
  }

  if(PVFFile){
    printf("Loading PVF %s\n",PVFFile);fflush(stdout);
    PVF = MRIread(PVFFile);
    if(PVF==NULL) exit(1);
    err = MRIdimMismatch(seg, PVF, 0);
    if(err){
      printf("ERROR: seg and PVF dim mismatch %d\n",err);
      exit(1);
    }
    printf("Computing SegTType\n");
    SegTType = Seg2TissueType(seg);
    //printf("Loading SegTType %s\n",SegTTypeFile);fflush(stdout);
    //SegTType = MatrixReadTxt(SegTTypeFile, NULL);
    if(SegTType == NULL){
      printf("ERROR: reading  %s\n",SegTTypeFile);
      exit(1);
    }
    //MatrixPrint(stdout,SegTType);
  }
  if(MaskVolFile){
    printf("Loading mask %s\n",MaskVolFile);fflush(stdout);
    mask = MRIread(MaskVolFile);
    if(mask==NULL) exit(1);
    err = MRIdimMismatch(seg, mask, 0);
    if(err){
      printf("ERROR: seg and mask dim mismatch %d\n",err);
      exit(1);
    }
  }
  printf("Data load time %4.1f sec\n",TimerStop(&mytimer)/1000.0);
  fflush(stdout);

  if(niterations > 0){
    printf("Starting iterative RBV\n");
    rbv = IterativeRBV0PVC(src, seg, mask, niterations,cFWHM, rFWHM, sFWHM, NULL);
    printf("Writing output to %s ...",OutVolFile);fflush(stdout); TimerStart(&mytimer) ;
    err = MRIwrite(rbv,OutVolFile);
    if(err){
      printf("ERROR: writing to %s\n",OutVolFile);
      exit(1);
    }
    printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
    printf("mri_rbvpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
    printf("mri_rbvpvc done\n");
    exit(0);
  }
  // get list of segmentation ids
  segidlist0 = MRIsegIdList(seg, &nsegs, 0);
  // remove 0 from the list
  segidlist = (int *)calloc(sizeof(int),nsegs);
  msegs = 0;
  for(nthseg = 0; nthseg < nsegs; nthseg++){
    if(segidlist0[nthseg] != 0){
      segidlist[msegs] = segidlist0[nthseg];
      msegs++;
    }
  }
  nsegs = msegs;
  free(segidlist0);
  printf("nsegs = %d, excluding segid=0\n",nsegs);

  // Create GTM matrix
  printf("Building GTM ... ");fflush(stdout); TimerStart(&mytimer) ;
  X = BuildGTMPVF(seg,SegTType,PVF,mask,cFWHM,rFWHM,sFWHM,NULL);
  if(X==NULL) exit(1);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  printf("nmask = %d\n",X->rows);
  //MatrixWriteTxt("X.mtx", X);
    
  // Fill y matrix - must be done in the same way that X was built!
  y = MatrixAlloc(X->rows,src->nframes,MATRIX_REAL);
  if(y==NULL){
    printf("ERROR: could not alloc y matrix %d\n",X->rows);
    exit(1);
  }
  nmask = 0;
  for(s=0; s < src->depth; s++){
    for(c=0; c < src->width; c++){
      for(r=0; r < src->height; r++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < src->nframes; f++)
	  y->rptr[nmask+1][f+1] = MRIgetVoxVal(src,c,r,s,f);
	nmask ++;
      }
    }
  }

  // Compute GTM means
  // beta = inv(X'*X)*X'*y
  // Should normalize X
  Xt = MatrixTranspose(X,NULL);
  printf("Computing  XtX ... ");fflush(stdout); TimerStart(&mytimer) ;
  XtX = MatrixMultiplyD(Xt,X,NULL);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  iXtX = MatrixInverse(XtX,NULL);
  printf("Computing  Xty ... ");fflush(stdout); TimerStart(&mytimer) ;
  Xty = MatrixMultiplyD(Xt,y,NULL);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);

  beta = MatrixMultiplyD(iXtX,Xty,NULL);

  XtXcond = MatrixConditionNumber(XtX);
  VRFStats(iXtX, &vrfmean, &vrfmin, &vrfmax);
  printf("XtX  Condition     %8.3f \n",XtXcond);
  printf("VRF  Mean/Min/Max  %8.3f %8.3f %8.3f \n",vrfmean,vrfmin,vrfmax);
  fflush(stdout);

  if(OutXtXFile){
    printf("Writing XtX to %s\n",OutXtXFile);
    if(err) exit(1);
    err=MatrixWriteTxt(OutXtXFile, XtX);
  }

  MatrixFree(&y);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&Xty);

  if(OutBetaFile){
    printf("Writing GTM estimates to %s\n",OutBetaFile);
    mritmp = MRIallocSequence(nsegs, 1, 1, 1, src->nframes);
    for(c=0; c < nsegs; c++){
      for(f=0; f < src->nframes; f++){
	MRIsetVoxVal(mritmp,c,0,0,f, beta->rptr[c+1][f+1]);
      }
    }
    err=MRIwrite(mritmp,OutBetaFile);
    if(err) exit(1);
    MRIfree(&mritmp);
    //err=MatrixWriteTxt(OutBetaFile, beta);
  }

  if(GTMOnly){
    printf("GTM-Only requested, so exiting now\n");
    printf("mri_rbvpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
    printf("mri_rbvpvc done\n");
    return(0);
    exit(0);
  }

  // Simulate Image by filling with GTM means
  printf("Synthesizing ... ");fflush(stdout); TimerStart(&mytimer) ;
  gtmsynth = MRIallocSequence(src->width, src->height, src->depth,
			      MRI_FLOAT, src->nframes);
  if(gtmsynth==NULL){
    printf("ERROR: could not alloc gtmsynth\n");
    exit(1);
  }
  MRIcopyHeader(src,gtmsynth);
  for(s=0; s < src->depth; s++){
    for(c=0; c < src->width; c++){
      for(r=0; r < src->height; r++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	segid = MRIgetVoxVal(seg,c,r,s,0);
	if(segid == 0) continue;
	for(nthseg=0; nthseg < nsegs; nthseg++)
	  if(segid == segidlist[nthseg]) break;
	for(f=0; f < src->nframes; f++)
	  MRIsetVoxVal(gtmsynth,c,r,s,f,beta->rptr[nthseg+1][f+1]);
      }
    }
  }
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);

  MatrixFree(&beta);

  // Smooth Simulated Images
  printf("Smoothing synthesized ... ");fflush(stdout); TimerStart(&mytimer) ;
  gtmsynthsm = MRIgaussianSmoothNI(gtmsynth, cStd, rStd, sStd, NULL);
  if(gtmsynthsm==NULL) exit(1);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);

  //MRIwrite(gtmsynth,"gtmsynth.nii");
  //MRIwrite(gtmsynthsm,"gtmsynthsm.nii");

  // Compute the final RBV
  printf("Computing RBV ... ");fflush(stdout); TimerStart(&mytimer) ;
  rbv = MRIallocSequence(src->width, src->height, src->depth,
			 MRI_FLOAT, src->nframes);
  if (rbv==NULL){
    printf("ERROR: could not alloc rbv\n");
    exit(1);
  }
  MRIcopyHeader(src,rbv);
  for(s=0; s < src->depth; s++){
    for(c=0; c < src->width; c++){
      for(r=0; r < src->height; r++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < src->nframes; f++){
	  val = (double)MRIgetVoxVal(src,c,r,s,f)*
	    MRIgetVoxVal(gtmsynth,c,r,s,f)/
	    (MRIgetVoxVal(gtmsynthsm,c,r,s,f)+FLT_EPSILON);
	  MRIsetVoxVal(rbv,c,r,s,f,val);
	}
      }
    }
  }
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);

  MRIfree(&gtmsynth);
  MRIfree(&gtmsynthsm);

  printf("Writing output to %s ...",OutVolFile);fflush(stdout); TimerStart(&mytimer) ;
  err = MRIwrite(rbv,OutVolFile);
  if(err){
    printf("ERROR: writing to %s\n",OutVolFile);
    exit(1);
  }
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);

  printf("mri_rbvpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
  printf("mri_rbvpvc done\n");
  return(0);
  exit(0);
}
/*--------------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--gtm-only")) GTMOnly = 1;
    else if (!strcasecmp(option, "--seg-test")) DoSegTest = 1;
    else if (!strcasecmp(option, "--old-dil"))setenv("GTMDILATEOLD","1",1);

    else if(!strcasecmp(option, "--src") || !strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      SrcVolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--seg")) {
      if (nargc < 1) CMDargNErr(option,1);
      SegVolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--pvf")) {
      if (nargc < 2) CMDargNErr(option,1);
      PVFFile = pargv[0];
      SegTTypeFile = pargv[1];
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      MaskVolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--psf")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&psfFWHM);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--niters")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&niterations);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--od")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutDir = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutVolFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--beta") || !strcasecmp(option, "--gtm-means")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutBetaFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--xtx")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutXtXFile = pargv[0];
      nargsused = 1;
    } 
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/*---------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --src volfile : source data to PVC\n");
  printf("   --seg volfile : segmentation to define regions for RBV\n");
  printf("   --mask volfile : ignore areas outside of the mask\n");
  printf("   --pvf pvffile segttype: Non-binary voxelwise PVF\n");
  printf("   --o  outvolfile : PVC'ed input\n");
  printf("   --gtm-means volfile : save ROI means in volume format\n");
  printf("   --gtm-only : only perform GTM\n");
  printf("   --xtx xtx.mtx : save X'*X into xtx.mtx\n");
  printf("   --niters N : use iterative method instead of GTM\n");
  printf("   --seg-test : replace input with seg smoothed by psf\n");
  //printf("   --od outdir     : output directory\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*---------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void check_options(void) 
{
  if(SrcVolFile == NULL && ! DoSegTest){
    printf("ERROR: must spec source volume\n");
    exit(1);
  }
  if(SegVolFile == NULL){
    printf("ERROR: must spec segmentation volume\n");
    exit(1);
  }
  if(psfFWHM < 0){
    printf("ERROR: must spec psf FWHM\n");
    exit(1);
  }
  if(!GTMOnly){
    if(OutVolFile == NULL && OutDir == NULL){
      printf("ERROR: must spec output with --od or --o\n");
      exit(1);
    }
  }
  else {
    if(OutBetaFile == NULL){
      printf("ERROR: must spec output with --gtm-means when using --gtm-only \n");
      exit(1);
    }
  }

  cFWHM = psfFWHM;
  rFWHM = psfFWHM;
  sFWHM = psfFWHM;

  return;
}
/*---------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  return;
}

/*-----------------------------------------------*/
MRI *IterativeRBV0PVC(MRI *src, MRI *seg, MRI *mask, int niters,
		      double cFWHM, double rFWHM, double sFWHM, 
		      MRI *dst)
{
  int *segidlist,segid,nsegs,nthseg,*nperseg;
  double **segmeans;
  int err, c,r,s, nthframe, nthiter;
  MRI *srcsynth=NULL,*srcsynthsm=NULL,*isrc;
  double cStd,rStd,sStd,val;

  err = MRIdimMismatch(seg, src, 0);
  if(err){
    printf("ERROR: BuildGTM(): src seg dim mismatch %d\n",err);
    return(NULL);
  }

  cStd = cFWHM/sqrt(log(256.0));
  rStd = rFWHM/sqrt(log(256.0));
  sStd = sFWHM/sqrt(log(256.0));
  
  segidlist = MRIsegIdList(seg, &nsegs, 0);

  segmeans = (double **) calloc(sizeof(double*),nsegs);
  for(nthseg=0; nthseg < nsegs; nthseg++)
    segmeans[nthseg] = (double *) calloc(sizeof(double),src->nframes);
  nperseg = (int *) calloc(sizeof(int),nsegs);

  if(dst == NULL){
    dst = MRIallocSequence(src->width,src->height,src->depth,
			   MRI_FLOAT,src->nframes);
    MRIcopyHeader(src,dst);
  }
  err = MRIdimMismatch(dst, src, 0);
  if(err){
    printf("ERROR: BuildGTM(): src dst dim mismatch %d\n",err);
    return(NULL);
  }

  srcsynth = MRIallocSequence(src->width,src->height,src->depth,
		      MRI_FLOAT,src->nframes);
  MRIcopyHeader(src,srcsynth);

  isrc = src;
  for(nthiter=0; nthiter < niters; nthiter++){
    printf("#@# iter %3d/%d ---------\n",nthiter,niters);
    fflush(stdout);

    // Compute seg means
    printf("  computing seg means\n");    fflush(stdout);
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	for(s=0; s < seg->depth; s++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	  segid = MRIgetVoxVal(seg,c,r,s,0);
	  for(nthseg=0; nthseg < nsegs; nthseg++)
	    if(segid == segidlist[nthseg]) break;
	  nperseg[nthseg]++;
	  for(nthframe=0; nthframe < src->nframes; nthframe++)
	    segmeans[nthseg][nthframe] += 
	      MRIgetVoxVal(isrc,c,r,s,nthframe);
	}
      }
    }

    for(nthseg=0; nthseg < nsegs; nthseg++)
      for(nthframe=0; nthframe < src->nframes; nthframe++)
	segmeans[nthseg][nthframe] /= nperseg[nthseg];

    // Create an image with segmeans
    printf("  synthesize\n");    fflush(stdout);
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	for(s=0; s < seg->depth; s++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	  segid = MRIgetVoxVal(seg,c,r,s,0);
	  for(nthseg=0; nthseg < nsegs; nthseg++)
	    if(segid == segidlist[nthseg]) break;
	  for(nthframe=0; nthframe < src->nframes; nthframe++)
	    MRIsetVoxVal(srcsynth,c,r,s,nthframe,segmeans[nthseg][nthframe]);
	}
      }
    }

    // Smooth the syntesized image
    printf("  smooth\n");    fflush(stdout);
    srcsynthsm = MRIgaussianSmoothNI(srcsynth, cStd, rStd, sStd, srcsynthsm);
    if(srcsynthsm == NULL) return(NULL);

    // Apply
    printf("  apply\n");    fflush(stdout);
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	for(s=0; s < seg->depth; s++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	  for(nthframe=0; nthframe < src->nframes; nthframe++){
	    val = MRIgetVoxVal(src,c,r,s,nthframe)*
	      MRIgetVoxVal(srcsynth,c,r,s,nthframe)/
	      MRIgetVoxVal(srcsynthsm,c,r,s,nthframe);
	    MRIsetVoxVal(dst,c,r,s,nthframe,val);
	  }
	}
      }
    }
    isrc = dst;

  }// Iteration

  MRIfree(&srcsynth);
  MRIfree(&srcsynthsm);
  free(nperseg);
  for(nthseg=0; nthseg < nsegs; nthseg++)
    free(segmeans[nthseg]);
  free(segmeans);

  return(dst);
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
MRI *MRIdownSmoothUp(MRI *src, int Fc, int Fr, int Fs, 
		     double cFWHM, double rFWHM, double sFWHM, 
		     MRI *dst)
{
  double cStd,rStd,sStd;
  MRI *dvol;

  if(Fr != Fc || Fs != Fc){
    printf("ERROR: MRIdownSmoothUp(): sampling factor must be iso\n");
    return(NULL);
  }

  if(Fc == 1 && Fr == 1 && Fs == 1)
    dst = MRIupsampleN(src, dst, Fc);
  else{
    printf("    downsample\n");fflush(stdout);
    dvol = MRIdownsampleN(src, NULL, Fc, Fr, Fs, 0);
    if(dvol == NULL) return(NULL);
    
    // Smooth unmasked
    printf("    smooth\n");fflush(stdout);
    cStd = cFWHM/sqrt(log(256.0));
    rStd = rFWHM/sqrt(log(256.0));
    sStd = sFWHM/sqrt(log(256.0));
    dvol = MRIgaussianSmoothNI(dvol, cStd, rStd, sStd, dvol);
    if(dvol == NULL) return(NULL);
    printf("    upsample\n");fflush(stdout);
    dst = MRIupsampleN(dvol, dst, Fc);
    MRIfree(&dvol);
  }

  return(dst);
}
/*---------------------------------------------------------------------*/
MATRIX *BuildGTM(MRI *seg, MATRIX *SegTType, MRI *pvf, MRI *mask,
		 int Fc, int Fr, int Fs, 
		 double cFWHM, double rFWHM, double sFWHM, 
		 MATRIX *X)
{
  int err, c,r,s,nmask, nsegs,nthseg, n;
  int segid, segtt;
  MRI *pvfmask=NULL,*pvfmasksm=NULL;
  MRI *pvflist[100];

  err = MRIdimMismatch(seg, pvf, 0);
  if(err){
    printf("ERROR: BuildGTM(): dim mismatch %d\n",err);
    return(NULL);
  }
  nsegs = SegTType->rows;

  // Count number of voxels in the mask
  nmask = 0;
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	nmask ++;
      }
    }
  }
  if(X==NULL) X = MatrixAlloc(nmask,nsegs,MATRIX_REAL);
  if(X->rows != nmask || X->cols != nsegs){
    printf("ERROR: BuildGTM(): X dim mismatch\n");
    return(NULL);
  }
  printf("nmask = %d, nsegs = %d\n",nmask,nsegs);

  for(n=0; n < pvf->nframes; n++)
    pvflist[n] = fMRIframe(pvf, n, NULL);

  // indWM = [2 41 7 46 16 28 60];
  // indCSF = 24;
  for(nthseg=0; nthseg < nsegs; nthseg++){
    segid = (int)SegTType->rptr[nthseg+1][1];
    segtt = (int)SegTType->rptr[nthseg+1][2];
    printf("#@# %3d/%d %3d %d ---------\n",nthseg,nsegs,segid,segtt);
    fflush(stdout);
    printf("  masking \n"); fflush(stdout);
    pvfmask = MRImask(pvflist[segtt], seg, pvfmask, segid, 0.0);
    printf("  down-smooth-up \n"); fflush(stdout);
    pvfmasksm = MRIdownSmoothUp(pvfmask, Fc, Fr, Fs, 
				cFWHM, rFWHM, sFWHM, pvfmasksm);
    printf("  filling \n"); fflush(stdout);
    nmask = 0;
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	for(s=0; s < seg->depth; s++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	  X->rptr[nmask+1][nthseg+1] = MRIgetVoxVal(pvfmasksm,c,r,s,0);
	  nmask ++;
	}
      }
    }
  } // seg

  MRIfree(&pvfmasksm);
  MRIfree(&pvfmask);
  for(n=0; n < pvf->nframes; n++) MRIfree(&pvflist[n]);

  return(X);
}

/*----------------------------------------------------------*/
MATRIX *Seg2TissueType(MRI *seg)
{
  int *segidlist, nsegs, nthseg, has0, segid, segtt,n;
  MATRIX *SegTType;

  segidlist = MRIsegIdList(seg, &nsegs, 0);

  //for(nthseg=0; nthseg < nsegs; nthseg++)
  //printf("%d\n",segidlist[nthseg]);

  has0 = 0;
  for(nthseg=0; nthseg < nsegs; nthseg++)
    if(segidlist[nthseg] == 0) has0 = 1;

  SegTType = MatrixAlloc(nsegs-has0,2,MATRIX_REAL);

  n = 0;
  for(nthseg=0; nthseg < nsegs; nthseg++){
    segid = segidlist[nthseg];
    if(segid == 0) continue;
    
    segtt = -1;
    if(segid == 2 || segid == 41 || segid == 7 ||
       segid == 46 || segid == 16 || segid == 28 ||
       segid == 60) segtt = 1; // WM
    else if(segid == 24) segtt = 2; // CSF
    else segtt = 0; // GM

    if(segtt == -1){
      printf("ERROR: Seg2TissueType(): segid %d unknown\n",segid);
      MatrixFree(&SegTType);
      return(NULL);
    }

    SegTType->rptr[n+1][1] = segid;
    SegTType->rptr[n+1][2] = segtt;
    n++;
  }
  return(SegTType);
}


/*---------------------------------------------------------------*/
MATRIX *GTM0PVC(MRI *src, MRI *seg, MRI *mask, 
		double cFWHM, double rFWHM, double sFWHM, 
		MATRIX **pX)
{
  int err, c,r,s,f,nmask;
  MATRIX *beta, *y, *X, *Xt, *XtX, *iXtX, *Xty;
  struct timeb  mytimer;

  err = MRIdimMismatch(seg, src, 0);
  if(err){
    printf("ERROR: SolveGTM(): src seg dim mismatch %d\n",err);
    return(NULL);
  }

  printf(" building X ... ");fflush(stdout);
  TimerStart(&mytimer) ;
  X = BuildGTM0(seg,mask,cFWHM,rFWHM,sFWHM,NULL);
  if(X==NULL) return(NULL);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  *pX = X;

  //MatrixWriteTxt("X.mtx", X);

  TimerStart(&mytimer) ;
  y = MatrixAlloc(X->rows,src->nframes,MATRIX_REAL);
  if(y==NULL){
    printf("ERROR: SolveGTM(): could not alloc y %d\n",X->rows);
    return(NULL);
  }
    
  // Fill y - must be done in the same way that X was built!
  nmask = 0;
  for(s=0; s < seg->depth; s++){
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < src->nframes; f++)
	  y->rptr[nmask+1][f+1] = MRIgetVoxVal(src,c,r,s,f);
	nmask ++;
      }
    }
  }

  // beta = inv(X'*X)*X'*y
  Xt = MatrixTranspose(X,NULL);
  printf("  XtX ... ");fflush(stdout);
  TimerStart(&mytimer) ;
  XtX = MatrixMultiplyD(Xt,X,NULL);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  iXtX = MatrixInverse(XtX,NULL);
  Xty = MatrixMultiplyD(Xt,y,NULL);
  beta = MatrixMultiplyD(iXtX,Xty,NULL);

  MatrixFree(&y);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&Xty);
  // Dont free X!

  return(beta);
}

/*---------------------------------------------------------------*/
MRI *GTM0PVCSynth(MRI *src, MRI *seg, MRI *mask, 
		double cFWHM, double rFWHM, double sFWHM, 
		  MRI *dst)
{
  MATRIX *beta, *X;
  struct timeb  mytimer;
  int c,r,s,f, *segidlist,segid,nthseg,nsegs,has0;

  if(dst == NULL){
    dst = MRIallocSequence(src->width, src->height, src->depth,
			   MRI_FLOAT, src->nframes);
    if (dst==NULL){
      printf("ERROR: GTM0PVCSynth(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(src,dst);
  }
  if(MRIdimMismatch(src,dst,1)){
    printf("ERROR: GTM0PVCSynth(): dim mismatch\n");
    return(NULL);
  }

  segidlist = MRIsegIdList(seg, &nsegs, 0);
  has0 = 0;
  for(nthseg=0; nthseg < nsegs; nthseg++)
    if(segidlist[nthseg]==0) has0 = 1;

  printf(" computing beta ... \n");fflush(stdout);
  TimerStart(&mytimer) ;
  beta = GTM0PVC(src, seg, mask, cFWHM, rFWHM, sFWHM, &X);
  if(beta==NULL) return(NULL);
  printf(" ... betatime = %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  MatrixFree(&X);

  //MatrixWriteTxt("beta.mtx", beta);

  // Fill y - must be done in the same way that X was built!
  for(s=0; s < seg->depth; s++){
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	segid = MRIgetVoxVal(seg,c,r,s,0);
	if(segid == 0) continue;
	for(nthseg=0; nthseg < nsegs; nthseg++)
	  if(segid == segidlist[nthseg]) break;
	for(f=0; f < src->nframes; f++)
	  MRIsetVoxVal(dst,c,r,s,f,beta->rptr[nthseg+1-has0][f+1]);
      }
    }
  }

  MatrixFree(&beta);

  return(dst);
}


/*-------------------------------------------------------------------*/
MRI *RBV0PVC(MRI *src, MRI *seg, MRI *mask, 
	     double cFWHM, double rFWHM, double sFWHM, 
	     MRI *rbv)
{
  MRI *srcsynth,*srcsynthsm;
  double cStd,rStd,sStd,val;
  int c,r,s,f;

  if(rbv == NULL){
    rbv = MRIallocSequence(src->width, src->height, src->depth,
			   MRI_FLOAT, src->nframes);
    if (rbv==NULL){
      printf("ERROR: RBV0PVC(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(src,rbv);
  }
  if(MRIdimMismatch(src,rbv,1)){
    printf("ERROR: RBV0PVC(): dim mismatch\n");
    return(NULL);
  }

  srcsynth = GTM0PVCSynth(src, seg, mask, cFWHM, rFWHM, sFWHM, NULL);
  if(srcsynth==NULL) return(NULL);

  cStd = cFWHM/sqrt(log(256.0));
  rStd = rFWHM/sqrt(log(256.0));
  sStd = sFWHM/sqrt(log(256.0));
  srcsynthsm = MRIgaussianSmoothNI(srcsynth, cStd, rStd, sStd, NULL);
  if(srcsynthsm==NULL) return(NULL);

  //MRIwrite(srcsynth,"srcsynth.nii");
  //MRIwrite(srcsynthsm,"srcsynthsm.nii");

  for(s=0; s < seg->depth; s++){
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < src->nframes; f++){
	  val = (double)MRIgetVoxVal(src,c,r,s,f)*
	    MRIgetVoxVal(srcsynth,c,r,s,f)/
	    (MRIgetVoxVal(srcsynthsm,c,r,s,f)+FLT_EPSILON);
	  MRIsetVoxVal(rbv,c,r,s,f,val);
	}
      }
    }
  }

  MRIfree(&srcsynth);
  MRIfree(&srcsynthsm);

  return(rbv);
}
/***************------------------------******************/
int VRFStats(MATRIX *iXtX, double *vrfmean, double *vrfmin, double *vrfmax)
{
  int n;
  double vrf;

  *vrfmean = 0;
  *vrfmax = 0;
  *vrfmin = 0;
  for(n=0; n < iXtX->rows; n++){
    vrf = (double) 1.0/iXtX->rptr[n+1][n+1];
    if(n==0){
      *vrfmax = vrf;
      *vrfmin = vrf;
    }
    if(*vrfmax < vrf) *vrfmax = vrf;
    if(*vrfmin > vrf) *vrfmin = vrf;
    *vrfmean += vrf;
  }
  *vrfmean /= iXtX->rows;
  return(0);
}

/***************------------------------******************/
MATRIX *MatrixGetDiag(MATRIX *M, VECTOR *d)
{
  int n;
  if(d==NULL) d = MatrixAlloc(M->rows,1,MATRIX_REAL);
  for(n=0; n < M->rows; n++)
    d->rptr[n+1][1] = M->rptr[n+1][n+1];
  return(d);
}

/***************------------------------******************/
MATRIX *BuildGTMPVF(MRI *seg, MATRIX *SegTType, MRI *pvf, MRI *mask, 
		    double cFWHM, double rFWHM, double sFWHM, 
		    MATRIX *X)
{
  int c,r,s,nmask, nsegs,nthseg, mthseg, kthseg, segid, *segidlist,has0,segtt;
  double cStd,rStd,sStd,val;
  MRI *roimask=NULL,*roimasksm=NULL,*segttvol[10],*mritmp,*pvftt[10];
  int nthtt,ndil=3,nchanges;
  int UseOld = 0;

  if(pvf == NULL) return(BuildGTM0(seg,mask,cFWHM,rFWHM,sFWHM,X));

  if(getenv("GTMDILATEOLD")) sscanf(getenv("GTMDILATEOLD"),"%d",&UseOld);
  if(UseOld) {printf("\nUsing old dilation method\n");fflush(stdout);}

  cStd = cFWHM/sqrt(log(256.0));
  rStd = rFWHM/sqrt(log(256.0));
  sStd = sFWHM/sqrt(log(256.0));

  // Count number of voxels in the mask
  nmask = 0;
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	nmask ++;
      }
    }
  }
  if(Gdiag_no > 0) printf("BuildGTM0(): nmask = %d, nsegs = %d\n",nmask,nsegs);

  // Get list of segids from the segmentation
  segidlist = MRIsegIdList(seg, &nsegs, 0);
  has0 = 0;
  for(nthseg=0; nthseg < nsegs; nthseg++)
    if(segidlist[nthseg] == 0) has0 = 1;

  printf("nsegs = %d, has0 = %d, SegTType %d\n",nsegs,has0,SegTType->rows);
  if(nsegs-has0 != SegTType->rows){
    printf("ERROR:Seg dim mismatch\n");
    return(NULL);
  }
  // Check that every segid in SegTType is in segidlist
  for(nthseg=0; nthseg < nsegs; nthseg++){
    segid = segidlist[nthseg];
    if(segid == 0) continue;
    segtt = -1;
    for(kthseg=0; kthseg < nsegs; kthseg++){
      if((int)SegTType->rptr[kthseg+1][1] == segid){
	segtt = (int)SegTType->rptr[kthseg+1][2];
	break;
      }
    }
    if(segtt == -1){
      printf("ERROR: cannot find a match for seg %d in SegTType\n",segid);
      return(NULL);
    }
  }

  // Create TT spec PVF MRI binarized/thresholded at 0.25
  for(nthtt = 0; nthtt < pvf->nframes; nthtt++){
    pvftt[nthtt] = fMRIframe(pvf, nthtt, NULL);
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	for(s=0; s < seg->depth; s++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	  val = MRIgetVoxVal(pvftt[nthtt],c,r,s,0);
	  if(val < 0.25) MRIsetVoxVal(pvftt[nthtt],c,r,s,0,0.0);
	  else           MRIsetVoxVal(pvftt[nthtt],c,r,s,0,1.0);
	}
      }
    }
  }

  printf("Creating tissue type specific seg\n");
  // Alloc tissue type specific seg
  for(nthtt = 0; nthtt < pvf->nframes; nthtt++){
    segttvol[nthtt] = MRIalloc(seg->width,seg->height,seg->depth,MRI_INT);
    MRIcopyHeader(seg,segttvol[nthtt]);
  }
  // Fill tissue type specific seg
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	segid = MRIgetVoxVal(seg,c,r,s,0);
	if(segid == 0) continue;
        segtt = -1;
	for(kthseg=0; kthseg < nsegs; kthseg++){
	  if((int)SegTType->rptr[kthseg+1][1] == segid){
	    segtt = (int)SegTType->rptr[kthseg+1][2];
	    break;
	  }
	}
	MRIsetVoxVal(segttvol[segtt],c,r,s,0,segid);
      }
    }
  }
  printf("Dilating tissue type specific seg by %d (%d)\n",ndil,UseOld);
  if(ndil != 0){
    for(nthtt = 0; nthtt < pvf->nframes; nthtt++){
      printf("  TType %d -----------\n",nthtt);
      if(UseOld){
	if(nthtt != 2) // gm and wm
	  mritmp = MRIdilateSegmentation(segttvol[nthtt], NULL, ndil, pvftt[nthtt], &nchanges);
	else // csf
	  mritmp = MRIdilateSegmentation(segttvol[nthtt], NULL, 3, NULL, &nchanges);
      }
      else mritmp = MRIdilateSegmentation(segttvol[nthtt], NULL, ndil, mask, &nchanges);
      MRIfree(&segttvol[nthtt]);
      segttvol[nthtt] = mritmp;
      printf("  TType %d had  %d changes\n",nthtt,nchanges);
    }
  }

  // Not sure why I threshold at 0.25 and do a fixed ndil=3 for CSF
  // Maybe there were so many voxels with very small PVF that they
  // overwhelmed the high PVF voxels. These would be the ones that
  // would be the most in error if the PSF is not accurate.
  // This matches what I did int matlab, but I don't know why I did 
  // it in the first place. CSF PVF will include sulcal but the seg will not.
  // Seems like thresholding should be done after smoothing, if at all

  // Alloc the ROI mask
  roimask = MRIconst(seg->width,seg->height,seg->depth,1,0.0,NULL);
  MRIcopyHeader(seg,roimask);

  // Alloc the GTM design matrix
  if(X==NULL) X = MatrixAlloc(nmask,nsegs-has0,MATRIX_REAL);
  if(X->rows != nmask || X->cols != nsegs-has0){
    printf("ERROR: BuildGTM0(): X dim mismatch\n");
    return(NULL);
  }

  // Create a regressor for each segid
  printf("Computing regressors\n");
  mthseg = 0;
  for(nthseg=0; nthseg < nsegs; nthseg++){
    segid = segidlist[nthseg];
    if(segid == 0) continue;

    // Get the tissue type for this seg
    segtt = -1;
    for(kthseg=0; kthseg < nsegs; kthseg++){
      if((int)SegTType->rptr[kthseg+1][1] == segid){
	segtt = (int)SegTType->rptr[kthseg+1][2];
	break;
      }
    }
    //printf("#@# %3d/%d %3d %d ---------\n",nthseg,nsegs,segid,segtt);
    if(Gdiag_no > 0) {
      printf("BuildGTM0(): #@# %3d/%d %3d ---\n",mthseg,nsegs-has0,segid); 
      fflush(stdout);
    }

    // Create a mask of the seg
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	for(s=0; s < seg->depth; s++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	  if(MRIgetVoxVal(segttvol[segtt],c,r,s,0) == segid) val = MRIgetVoxVal(pvf,c,r,s,segtt);
	  else                                               val = 0;
	  MRIsetVoxVal(roimask,c,r,s,0,val);
	}
      }
    }
    // Smooth the mask
    roimasksm = MRIgaussianSmoothNI(roimask, cStd, rStd, sStd, roimasksm);
    // Fill X
    // Creating X in this order makes it consistent with matlab
    // Note: y must be ordered in the same way.
    nmask = 0;
    for(s=0; s < seg->depth; s++){
      for(c=0; c < seg->width; c++){
	for(r=0; r < seg->height; r++){
	  if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	  X->rptr[nmask+1][mthseg+1] = MRIgetVoxVal(roimasksm,c,r,s,0);
	  nmask ++;
	}
      }
    }
    mthseg ++;
  } // seg

  for(nthtt = 0; pvf->nframes < nthtt; nthtt++) {
    MRIfree(&segttvol[nthtt]);
    MRIfree(&pvftt[nthtt]);
  }
  MRIfree(&roimask);
  MRIfree(&roimasksm);
  free(segidlist);

  return(X);
}
