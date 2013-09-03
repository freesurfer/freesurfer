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
 *    $Date: 2013/09/03 21:32:33 $
 *    $Revision: 1.2 $
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


// $Id: mri_rbvpvc.c,v 1.2 2013/09/03 21:32:33 greve Exp $

/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


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

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_rbvpvc.c,v 1.2 2013/09/03 21:32:33 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *SrcVolFile=NULL,*SegVolFile=NULL,*MaskVolFile=NULL;
char *OutDir=NULL,*OutVolFile=NULL;
char *OutBetaFile=NULL;
double psfFWHM=-1,cFWHM,rFWHM,sFWHM,cStd,rStd,sStd;
int GTMOnly=0;
char tmpstr[5000];

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err,c,r,s,f,nmask;
  MRI *src, *seg, *mask=NULL,*rbv,*gtmsynth,*gtmsynthsm,*mritmp;
  int nsegs,msegs,nthseg,*segidlist0,*segidlist,segid;
  MATRIX *beta, *y, *X, *Xt, *XtX, *iXtX, *Xty;
  double val,XtXcond;
  struct timeb  mytimer, timer;

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

  // Load the data
  printf("Loading seg %s\n",SegVolFile);fflush(stdout);
  seg = MRIread(SegVolFile);
  if(seg==NULL) exit(1);
  printf("Loading input %s\n",SrcVolFile);fflush(stdout);
  src = MRIread(SrcVolFile);
  if(src==NULL) exit(1);
  err = MRIdimMismatch(seg, src, 0);
  if(err){
    printf("ERROR: seg and source dim mismatch %d\n",err);
    exit(1);
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
  TimerStart(&timer);

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
  X = BuildGTM0(seg,mask,cFWHM,rFWHM,sFWHM,NULL);
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
  XtXcond = MatrixConditionNumber(XtX);
  printf("XtX Condition  %8.3f \n",XtXcond);fflush(stdout);

  iXtX = MatrixInverse(XtX,NULL);
  printf("Computing  Xty ... ");fflush(stdout); TimerStart(&mytimer) ;
  Xty = MatrixMultiplyD(Xt,y,NULL);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);

  beta = MatrixMultiplyD(iXtX,Xty,NULL);

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
  MatrixFree(&beta);

  // Smooth Simulated Images
  gtmsynthsm = MRIgaussianSmoothNI(gtmsynth, cStd, rStd, sStd, NULL);
  if(gtmsynthsm==NULL) exit(1);

  //MRIwrite(gtmsynth,"gtmsynth.nii");
  //MRIwrite(gtmsynthsm,"gtmsynthsm.nii");

  // Compute the final RBV
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

  MRIfree(&gtmsynth);
  MRIfree(&gtmsynthsm);

  printf("Writing output to %s\n",OutVolFile);
  err = MRIwrite(rbv,OutVolFile);
  if(err){
    printf("ERROR: writing to %s\n",OutVolFile);
    exit(1);
  }

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
  printf("   --o  outvolfile : PVC'ed input\n");
  printf("   --gtm-means volfile : save ROI means in volume format\n");
  printf("   --gtm-only : only perform GTM\n");
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
  if(SrcVolFile == NULL){
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
