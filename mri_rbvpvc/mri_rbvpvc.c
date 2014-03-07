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
 *    $Date: 2014/03/07 23:12:57 $
 *    $Revision: 1.27 $
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


// $Id: mri_rbvpvc.c,v 1.27 2014/03/07 23:12:57 greve Exp $

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
#include <sys/stat.h>
#include <sys/types.h>

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
#include "cma.h"
#include "mrimorph.h"
#include "region.h"
#include "resample.h"
#include "numerics.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_rbvpvc.c,v 1.27 2014/03/07 23:12:57 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

typedef struct 
{
  MRI *yvol;
  MRI *anatseg;
  LTA *anatseg2pet;
  MRI *mask;
  double cFWHM, rFWHM, sFWHM;
  double PadThresh;

  MRI *segpvf;
  MRI *gtmseg;
  int nPad;
  int nmask;
  int nsegs,*segidlist;
  int dof;
  COLOR_TABLE *ctGTMSeg;
  MATRIX *X,*X0;
  MATRIX *y,*Xt, *XtX, *iXtX, *Xty, *beta, *res, *yhat;
  double XtXcond;
  MATRIX *rvar;
  double cStd, rStd, sStd;
  MRI *ysynth,*ysynthsm;
  int *nperseg;
  MATRIX *segstats;
  MRI *rbv;
  MRI *mg;
  double mg_gmthresh;
  int *mg_refids;
  int n_mg_refids;
  MATRIX *mg_reftac;
} GTM;

GTM *GTMalloc();
int GTMfree(GTM **pGTM);
int GTMmatrixY(GTM *gtm);
int GTMsetNMask(GTM *gtm);
int GTMpsfStd(GTM *gtm);
int GTMsegidlist(GTM *gtm);
int GTMnPad(GTM *gtm);
int GTMbuildX(GTM *gtm);
int GTMsolve(GTM *gtm);
int GTMsegrvar(GTM *gtm);
int GTMsmoothSynth(GTM *gtm);
int GTMsynth(GTM *gtm);
int GTMsynth2(GTM *gtm);
int GTMrbv(GTM *gtm);
int GTMmgpvc(GTM *gtm);
MATRIX *GTMvol2mat(GTM *gtm, MRI *vol, MATRIX *m);
MRI *GTMmat2vol(GTM *gtm, MATRIX *m, MRI *vol);

typedef struct 
{
  GTM  *gtm;
  LTA  *anat2pet0,*anat2pet; // has subject name
  char *gtmsegfile;
  char *pvfsegfile;
  char *wsurf;
  char *psurf;
  int USF;
  COLOR_TABLE *ctPVFSeg;
  MRI *gtmseganat, *gtmsegmu, *pvfseganat, *pvfseg;
  MRIS *lhw, *lhp, *rhw, *rhp;
  LTA *anat2gtmsegmu, *gtmsegmu2anat,*gtmsegmu2pet;
  LTA *pvfseg2anat,*pvfseg2pet,*anat2pvfseg;
  float params[100];
  int nparams;
  double ftol;
  double linmintol;
  float fret;
  int niters,nitersmax;
  int nCostEvaluations;
  double tLastEval;
} GTMOPT;
int GTMOPTsetup(GTMOPT *gtmopt);
double GTMOPTcost(GTMOPT *gtmopt);

float compute_powell_cost(float *p);
int MinPowell();

char *SrcVolFile=NULL,*SegVolFile=NULL,*MaskVolFile=NULL;
char *OutDir=NULL,*RBVVolFile=NULL;
char *OutBetaFile=NULL,*OutXtXFile=NULL;
double psfFWHM=-1;
char tmpstr[5000];
char *PVFFile=NULL, *SegTTypeFile=NULL;
double ApplyFWHM=0;
char *Xfile=NULL;
char *VRFStatsFile=NULL;
char *eresFile=NULL, *yhatFile=NULL, *yhatFile0=NULL;
char *SynthFile=NULL;
MRI *mritmp;
char *RVarFile=NULL;
int RVarOnly=0;
int nthreads=1;
int nReplace = 0, SrcReplace[1000], TrgReplace[1000];
int ttReduce = 0;
char *MGPVCFile=NULL;

int VRFStats(MATRIX *iXtX, double *vrfmean, double *vrfmin, double *vrfmax);
int WriteVRFStats(char *fname, GTM *gtm);


MATRIX *GTMttest(MATRIX *beta,MATRIX *iXtX, double *rvar, int nthseg);

GTM *gtm;
GTMOPT *gtmopt;
GTMOPT *gtmopt_powell;

LTA *LTAapplyAffineParametersTKR(LTA *inlta, const float *p, const int np, LTA *outlta);
int DoOpt=0;
char *SUBJECTS_DIR;
int CheckX(MATRIX *X);
int dngtest(LTA *aseg2vol);

MRI *MRIaseg2volFrame(MRI *aseg, LTA *aseg2vol, int *segidlist, int nsegs, MRI *out);

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err,c,f,nTT,tt;
  MRI *pvfstack;
  double vrfmin,vrfmax,vrfmean;
  struct timeb  mytimer, timer;
  MRI *gtmres;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  gtm = GTMalloc();
  gtm->ctGTMSeg = TissueTypeSchema(NULL,"default-jan-2014");
  gtm->nDils = 4;
  gtm->n_mg_refids = 2;
  gtm->mg_refids = (int*) calloc(sizeof(int),gtm->n_mg_refids);
  gtm->mg_refids[0] =  2;
  gtm->mg_refids[1] = 41;

  //ctSeg = TissueTypeSchema(NULL,"default-jan-2014");

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
  nTT = gtm->ctGTMSeg->ctabTissueType->nentries-1;

#ifdef _OPENMP
  printf("%d avail.processors, using %d\n",omp_get_num_procs(),omp_get_max_threads());
#endif

  if(OutDir != NULL) {
    printf("Creating output directory %s\n",OutDir);
    err = mkdir(OutDir,0777);
    if (err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n",OutDir);
      perror(NULL);
      return(1);
    }
  }

  TimerStart(&timer);

  // Load data
  TimerStart(&mytimer);
  printf("Loading seg for gtm %s\n",SegVolFile);fflush(stdout);
  gtm->gtmseg = MRIread(SegVolFile);
  if(gtm->gtmseg==NULL) exit(1);
  if(nReplace > 0) {
    printf("Replacing %d\n",nReplace);
    for(f=0; f < nReplace; f++) printf("%2d:  %4d %4d\n",f+1,SrcReplace[f],TrgReplace[f]);
    mritmp = MRIreplaceList(gtm->gtmseg, SrcReplace, TrgReplace, nReplace, NULL);
    MRIfree(&gtm->gtmseg);
    gtm->gtmseg = mritmp;
  }
  if(ttReduce > 0) {
    MRI *ttseg;
    COLOR_TABLE *ctTT;
    printf("Reducing seg to tissue type seg\n");
    ttseg = MRIseg2TissueType(gtm->gtmseg, gtm->ctGTMSeg, NULL);
    if(ttseg == NULL) exit(1);
    MRIfree(&gtm->gtmseg);
    gtm->gtmseg = ttseg;
    ctTT = CTABdeepCopy(gtm->ctGTMSeg->ctabTissueType);
    for(f=0; f < ctTT->nentries; f++) ctTT->entries[f]->TissueType=f;
    ctTT->ctabTissueType = CTABdeepCopy(gtm->ctGTMSeg->ctabTissueType);
    CTABfree(&gtm->ctGTMSeg);
    gtm->ctGTMSeg = ctTT;
  }

  printf("Loading input %s\n",SrcVolFile);
  gtm->yvol = MRIread(SrcVolFile);
  if(gtm->yvol==NULL) exit(1);
  if(ApplyFWHM > 0){
    double cStdApply, rStdApply, sStdApply;
    MRI *mritmp;
    printf("Smoothing input by %g mm FWHM \n",ApplyFWHM);
    cStdApply = ApplyFWHM/sqrt(log(256.0));
    rStdApply = ApplyFWHM/sqrt(log(256.0));
    sStdApply = ApplyFWHM/sqrt(log(256.0));
    mritmp = MRIgaussianSmoothNI(gtm->yvol, cStdApply, rStdApply, sStdApply, NULL);
    MRIfree(&gtm->yvol);
    gtm->yvol = mritmp;
  }
  err = MRIdimMismatch(gtm->gtmseg, gtm->yvol, 0);
  if(err){
    printf("ERROR: seg and source dim mismatch %d\n",err);
    exit(1);
  }
  if(MaskVolFile){
    printf("Loading mask %s\n",MaskVolFile);fflush(stdout);
    gtm->mask = MRIread(MaskVolFile);
    if(gtm->mask==NULL) exit(1);
    err = MRIdimMismatch(gtm->gtmseg, gtm->mask, 0);
    if(err){
      printf("ERROR: seg and mask dim mismatch %d\n",err);
      exit(1);
    }
  }
  else gtm->mask=NULL; // should already be NULL

  printf("Data load time %4.1f sec\n",TimerStop(&mytimer)/1000.0);

  GTMsetNMask(gtm);
  GTMsegidlist(gtm);
  GTMmatrixY(gtm);
  GTMpsfStd(gtm); 
  GTMnPad(gtm);   
  //GTMdilateSeg(gtm);
  printf("nmask = %d, nsegs = %d, excluding segid=0\n",gtm->nmask,gtm->nsegs);
  printf("FWHM: %g %g %g\n",gtm->cFWHM,gtm->rFWHM,gtm->sFWHM);
  printf("Std:  %g %g %g\n",gtm->cStd,gtm->rStd,gtm->sStd);
  printf("nPad %d, PadThresh %g\n",gtm->nPad,gtm->PadThresh);

  if(DoOpt){
  printf("\nRunning optimization\n");
  gtmopt = (GTMOPT *) calloc(sizeof(GTMOPT),1);
  gtmopt->anat2pet0 = LTAread("register.lta");
  gtmopt->anat2pet  = LTAcopy(gtmopt->anat2pet0,NULL);
  gtmopt->gtm = gtm;
  GTMOPTsetup(gtmopt);
  gtmopt_powell = gtmopt;
  MinPowell();
  LTAwrite(gtmopt->anat2pet,"gtmopt.reg.lta");
  //for(f=0; f < 100; f++){
  //printf("#@# %d ------------------------------------------\n",f);
  //GTMOPTcost(gtmopt);
  //}
  printf("Done optimization\n\n");
  //exit(1);
  }

  // Create GTM matrix
  printf("Building GTM ... ");fflush(stdout); 
  TimerStart(&mytimer) ;
  if(1){
    printf("Reading reg\n");fflush(stdout);
    gtm->anatseg2pet = LTAread("register.lta");
    if(gtm->anatseg2pet==NULL) exit(1);
    printf("Reading anat\n");fflush(stdout);
    sprintf(tmpstr,"%s/%s/mri/petseg+head.mgz",SUBJECTS_DIR,gtm->anatseg2pet->subject);
    gtm->anatseg = MRIread(tmpstr);
    if(gtm->anatseg==NULL) exit(1);
    GTMbuildX(gtm);
    printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
    CheckX(gtm->X0);
  }

  if(gtm->X==NULL) exit(1);
  if(Xfile) {
    printf("Writing X to %s\n",Xfile);
    MatrixWriteTxt(Xfile, gtm->X);
  }
  printf("Solving ...\n");
  TimerStart(&mytimer) ; 
  GTMsolve(gtm);
  printf("Time to solve %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);

  if(OutBetaFile){
    printf("Writing GTM estimates to %s\n",OutBetaFile);
    mritmp = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, gtm->yvol->nframes);
    for(c=0; c < gtm->nsegs; c++){
      for(f=0; f < gtm->yvol->nframes; f++){
	MRIsetVoxVal(mritmp,c,0,0,f, gtm->beta->rptr[c+1][f+1]);
      }
    }
    err=MRIwrite(mritmp,OutBetaFile);
    if(err) exit(1);
    MRIfree(&mritmp);
    //err=MatrixWriteTxt(OutBetaFile, beta);
  }

  printf("rvar = %g\n",gtm->rvar->rptr[1][1]);
  gtm->XtXcond = MatrixConditionNumber(gtm->XtX);
  printf("XtX  Condition     %8.3f \n",gtm->XtXcond);

  printf("Synthesizing ... ");fflush(stdout); TimerStart(&mytimer) ;
  GTMsynth(gtm);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  if(yhatFile0) MRIwrite(gtm->ysynth,yhatFile0);

  printf("Smoothing synthesized ... ");fflush(stdout); TimerStart(&mytimer) ;
  GTMsmoothSynth(gtm);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  if(yhatFile) MRIwrite(gtm->ysynthsm,yhatFile);

  gtmres = GTMmat2vol(gtm,gtm->res,NULL);
  if(eresFile) MRIwrite(gtmres,eresFile);

  if(RVarFile) {
    FILE *fp;
    fp = fopen(RVarFile,"w");
    for(f=0; f < gtm->yvol->nframes; f++)
      fprintf(fp,"%30.20f\n",gtm->rvar->rptr[1][f+1]);
    fclose(fp);
    if(RVarOnly){
      printf("rvar-only requested so exiting now\n");
      printf("mri_rbvpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
      exit(0);
    }
  }
  GTMsegrvar(gtm);
  if(VRFStatsFile) WriteVRFStats(VRFStatsFile, gtm);

  VRFStats(gtm->iXtX, &vrfmean, &vrfmin, &vrfmax);
  printf("XtX  Condition     %8.3f \n",gtm->XtXcond);
  printf("VRF  Mean/Min/Max  %8.3f %8.3f %8.3f \n",vrfmean,vrfmin,vrfmax);
  fflush(stdout);


  if(MGPVCFile != NULL){
    printf("MG PVC\n");
    GTMmgpvc(gtm);
    err = MRIwrite(gtm->mg,MGPVCFile);
    if(err) exit(1);
    printf("done with mgpvc\n");
  }
    
  if(RBVVolFile){
    printf("Computing RBV\n");
    GTMrbv(gtm);
    printf("Writing output to %s ...",RBVVolFile);fflush(stdout); TimerStart(&mytimer) ;
    err = MRIwrite(gtm->rbv,RBVVolFile);
    if(err){
      printf("ERROR: writing to %s\n",RBVVolFile);
      exit(1);
    }
    printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
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
    //else if (!strcasecmp(option, "--seg-test")) DoSegTest = 1;
    else if (!strcasecmp(option, "--old-dil"))setenv("GTMDILATEOLD","1",1);
    else if (!strcasecmp(option, "--gtm-only")) ; // not used anymore
    else if (!strcasecmp(option, "--opt")) DoOpt=1;
    else if (!strcasecmp(option, "--ttype+head"))
      gtm->ctGTMSeg = TissueTypeSchema(NULL,"default-jan-2014+head");
    else if (!strcasecmp(option, "--tt-reduce")) ttReduce = 1;

    else if (!strcmp(option, "--sd") || !strcmp(option, "-SDIR")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      nargsused = 1;
    } 
    else if (!strcmp(option, "--reg")){
      if(nargc < 1) CMDargNErr(option,1);
      gtm->anatseg2pet = LTAread(pargv[0]);
      if(gtm->anatseg2pet==NULL) exit(1);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--anat-seg")){
      if(nargc < 1) CMDargNErr(option,1);
      gtm->anatseg = MRIread(pargv[0]);
      if(gtm->anatseg==NULL) exit(1);
      nargsused = 1;
    } 
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
      if (nargc < 1) CMDargNErr(option,1);
      PVFFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      MaskVolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--gdiag")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--psf")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&psfFWHM);
      gtm->cFWHM = psfFWHM;
      gtm->rFWHM = psfFWHM;
      gtm->sFWHM = psfFWHM;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--psf-col")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&gtm->cFWHM);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--psf-row")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&gtm->rFWHM);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--psf-slice")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&gtm->sFWHM);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--threads")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
#ifdef _OPENMP
      omp_set_num_threads(nthreads);
#endif
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--apply-fwhm")){
      // apply to input for testing
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&ApplyFWHM);
      nargsused = 1;
    } 
    //else if (!strcasecmp(option, "--niters")){
    //if(nargc < 1) CMDargNErr(option,1);
    //sscanf(pargv[0],"%d",&niterations);
    //nargsused = 1;
    //} 
    else if (!strcasecmp(option, "--ndil")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&gtm->nDils);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o") || !strcasecmp(option, "--rbv")) {
      if (nargc < 1) CMDargNErr(option,1);
      RBVVolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--eres")) {
      if (nargc < 1) CMDargNErr(option,1);
      eresFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--yhat")) {
      if (nargc < 1) CMDargNErr(option,1);
      yhatFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--yhat0")) {
      if (nargc < 1) CMDargNErr(option,1);
      yhatFile0 = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mgpvc")) {
      if (nargc < 2) CMDargNErr(option,2);
      MGPVCFile = pargv[0];
      sscanf(pargv[1],"%lf",&gtm->mg_gmthresh);
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--beta") || !strcasecmp(option, "--gtm-means")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutBetaFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--X")) {
      if (nargc < 1) CMDargNErr(option,1);
      Xfile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--xtx")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutXtXFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--vrf")) {
      if (nargc < 1) CMDargNErr(option,1);
      VRFStatsFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--rvar")) {
      if (nargc < 1) CMDargNErr(option,1);
      RVarFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--rvar-only"))
      RVarOnly = 1;
    else if(!strcasecmp(option, "--odir")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutDir = pargv[0];
      sprintf(tmpstr,"%s/rvar.dat",OutDir);
      RVarFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/vrf.dat",OutDir);
      VRFStatsFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/eres.nii.gz",OutDir);
      eresFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/yhat.nii.gz",OutDir);
      yhatFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/gtm.nii.gz",OutDir);
      OutBetaFile = strcpyalloc(tmpstr);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--replace")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&SrcReplace[nReplace]);
      sscanf(pargv[1],"%d",&TrgReplace[nReplace]);
      nReplace++;
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--synth")) {
      MRI *pvfstack;
      int tt,nTT;
      if (nargc < 5) CMDargNErr(option,5);
      OutBetaFile = pargv[0];
      SegVolFile = pargv[1];
      PVFFile = pargv[2];
      MaskVolFile = pargv[3];
      SynthFile = pargv[4];
      mritmp = MRIread(OutBetaFile);
      gtm->beta = fMRItoMatrix(mritmp,NULL);
      gtm->beta = MatrixTranspose(gtm->beta,NULL);
      gtm->gtmseg = MRIread(SegVolFile);
      pvfstack = MRIread(PVFFile);
      nTT = gtm->ctGTMSeg->ctabTissueType->nentries-1;
      gtm->pvf = (MRI **) calloc(sizeof(MRI*),nTT);
      for(tt=0; tt<nTT; tt++) gtm->pvf[tt] = fMRIframe(pvfstack,tt,NULL);
      gtm->mask = MRIread(MaskVolFile);
      GTMsetNMask(gtm);
      GTMsegidlist(gtm);
      GTMdilateSeg(gtm);
      GTMsynth(gtm);
      MRIwrite(gtm->ysynth,SynthFile);
      exit(0);
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
  printf("   --seg segfile : segmentation to define regions for RBV\n");
  printf("   --mask volfile : ignore areas outside of the mask\n");
  printf("   --psf psfmm : scanner PSF FWHM in mm\n");
  printf("   --pvf pvffile : Non-binary voxelwise PVF\n");
  printf("   --gtm-means volfile : save ROI means in volume format\n");
  printf("   --rbv rbvfile : PVC'ed input\n");
  printf("   --mgpvc mgpvcfile gmthresh\n");
  printf("   --vrf vrfstatsfile\n");
  printf("   --yhat gtm yhat file\n");
  printf("   --eres  gtm residual file\n");
  printf("   --seg-test : replace input with seg smoothed by psf\n");
  printf("   --xtx xtx.mtx : save X'*X into xtx.mtx\n");
  printf("   --X Xfile : save X matrix (it will be big)\n");
  printf("   --niters N : use iterative method instead of GTM\n");
  printf("   --odir outdir     : output directory\n");
  printf("   --synth gtmbeta seg pvf mask out\n");
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
  //if(SrcVolFile == NULL && ! DoSegTest){
  if(SrcVolFile == NULL){
    printf("ERROR: must spec source volume\n");
    exit(1);
  }
  if(SegVolFile == NULL){
    printf("ERROR: must spec segmentation volume\n");
    exit(1);
  }
  if(gtm->cFWHM < 0 || gtm->rFWHM < 0 || gtm->sFWHM < 0){
    printf("ERROR: must spec psf FWHM\n");
    exit(1);
  }
  if(RBVVolFile == NULL && OutDir == NULL && MGPVCFile == NULL 
     && OutBetaFile == NULL && RVarFile == NULL && yhatFile == NULL){
    printf("ERROR: must spec an output with --o, --gtm-means,  or --mgpvc\n");
    exit(1);
  }
  if(nthreads != 1){
#ifdef _OPENMP
    printf("Setting maximum number of threads to %d\n",nthreads);
    omp_set_num_threads(nthreads);
#endif
  }
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
/*--------------------------------------------------------------------------*/
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
/*--------------------------------------------------------------------------*/
int WriteVRFStats(char *fname, GTM *gtm)
{
  int n, segid, nvox;
  double vrf;
  FILE *fp;
  CTE *cte;
  CT *ttctab;

  ttctab = gtm->ctGTMSeg->ctabTissueType;

  fp = fopen(fname,"w");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",fname);
    return(1);
  }

  for(n=0; n < gtm->iXtX->rows; n++){
    segid = gtm->segidlist[n];
    nvox = MRIcountMatches(gtm->gtmseg, segid, 0, gtm->mask);
    vrf = (double) 1.0/gtm->iXtX->rptr[n+1][n+1];
    cte = gtm->ctGTMSeg->entries[segid];
    fprintf(fp,"%3d %4d %-31s %-13s %6d %8.3f",n+1,segid,cte->name,
	    ttctab->entries[cte->TissueType]->name,nvox,vrf);
    //printf("%3d %4d %-31s %-13s %6d %8.3f",n+1,segid,cte->name,
    //ttctab->entries[cte->TissueType]->name,nvox,vrf);
    if(gtm->beta)
      fprintf(fp,"   %10.3f",gtm->beta->rptr[n+1][1]);
    if(gtm->segstats){
      fprintf(fp,"   %10.4f",sqrt(gtm->segstats->rptr[n+1][1]));
      //printf("   %10.3f",gtm->rptr[n+1][1]);
    }
    fprintf(fp,"\n");
    //printf(" \n");
  }
  fclose(fp);
  fflush(stdout);

  return(0);
}
/*-------------------------------------------------------------------------------*/
MATRIX *GTMttest(MATRIX *beta,MATRIX *iXtX, double *rvar, int nthseg)
{
  MATRIX *C,*Ct,*CiXtX,*CiXtXCt,*gamma,*t;
  int f,nframes;

  nframes = beta->cols;
  t = MatrixAlloc(1,nframes,MATRIX_REAL);

  C = MatrixAlloc(1,beta->rows,MATRIX_REAL);
  C->rptr[1][nthseg+1] = 1;
  gamma = MatrixMultiply(C,beta,NULL);
  Ct = MatrixTranspose(C,NULL);
  CiXtX = MatrixMultiply(C,iXtX,NULL);
  CiXtXCt = MatrixMultiply(CiXtX,Ct,NULL);
  for(f=0; f < nframes; f++)
    t->rptr[1][f+1] = gamma->rptr[1][f+1]/sqrt(CiXtXCt->rptr[1][1]*rvar[f]);
  return(t);
}
/*-------------------------------------------------------------------------------*/



/*-------------------------------------------------------------------------------*/
int GTMOPTsetup(GTMOPT *gtmopt)
{
  char *subject, *SUBJECTS_DIR;
  //gtmopt->gtmsegfile="petseg+head.mgz";
  //gtmopt->pvfsegfile="petseg+head.mgz";
  gtmopt->gtmsegfile="petseg.mgz";
  gtmopt->pvfsegfile="petseg.mgz";
  gtmopt->wsurf = "white";
  gtmopt->psurf = "pial";
  gtmopt->USF = 1;
  gtmopt->ftol= 1e-8;
  gtmopt->linmintol= 5e-3;
  gtmopt->nitersmax = 4;
  gtmopt->nparams = 6;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  subject = gtmopt->anat2pet->subject;

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,gtmopt->pvfsegfile);
  printf("Loading %s\n",tmpstr);
  gtmopt->pvfseganat = MRIread(tmpstr);
  if(gtmopt->pvfseganat==NULL) exit(1);
  
  sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,gtmopt->wsurf);
  gtmopt->lhw = MRISread(tmpstr);
  if(gtmopt->lhw==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,gtmopt->psurf);
  gtmopt->lhp = MRISread(tmpstr);
  if(gtmopt->lhp==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,gtmopt->wsurf);
  gtmopt->rhw = MRISread(tmpstr);
  if(gtmopt->rhw==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,gtmopt->psurf);
  gtmopt->rhp = MRISread(tmpstr);
  if(gtmopt->rhp==NULL) exit(1);

  // Create a high resolution segmentation used to create PVF
  gtmopt->pvfseg = MRIhiresSeg(gtmopt->pvfseganat,
			       gtmopt->lhw,gtmopt->lhp,gtmopt->rhw,gtmopt->rhp, 
			       0, gtmopt->USF, &gtmopt->anat2pvfseg);
  if(gtmopt->pvfseg == NULL) return(1);
  gtmopt->pvfseg2anat = LTAinvert(gtmopt->anat2pvfseg,NULL);

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,gtmopt->gtmsegfile);
  gtmopt->gtmseganat = MRIread(tmpstr);
  if(gtmopt->gtmseganat==NULL) exit(1);
  gtmopt->gtmsegmu = MRImaskAndUpsample(gtmopt->gtmseganat, NULL, gtmopt->USF, 2, 0, &gtmopt->anat2gtmsegmu);
  gtmopt->gtmsegmu2anat = LTAinvert(gtmopt->anat2gtmsegmu,NULL);

  GTMsetNMask(gtmopt->gtm);
  GTMmatrixY(gtmopt->gtm);
  GTMpsfStd(gtmopt->gtm); 
  GTMnPad(gtmopt->gtm);   

  return(0);
}
/*--------------------------------------------------------------------------*/
double GTMOPTcost(GTMOPT *gtmopt)
{
  MRI *hitvol;
  int nTT,n;
  LTA *ltaArray[2];
  struct timeb timer;
  TimerStart(&timer);

  //printf("GTMOPTcost() USF=%d\n",gtmopt->USF);
  nTT = gtmopt->gtm->ctGTMSeg->ctabTissueType->nentries-1;

  ltaArray[0] = gtmopt->pvfseg2anat;
  ltaArray[1] = gtmopt->anat2pet;
  gtmopt->pvfseg2pet = LTAconcat(ltaArray,2,1);
  //printf("Computing PVF\n");
  if(gtmopt->gtm->pvf != NULL)
    for(n=0; n < nTT; n++) 
      if(gtmopt->gtm->pvf[n]) MRIfree(&(gtmopt->gtm->pvf[n]));
  //  gtmopt->gtm->pvf = MRIpartialVolumeFraction(gtmopt->pvfseg2pet,gtmopt->pvfseg, 
  //					      -1, gtmopt->gtm->ctGTMSeg,gtmopt->gtm->pvf);
  if(gtmopt->gtm->pvf == NULL) return(-1);
  //printf("  PVF t = %g\n",TimerStop(&timer)/1000.0);fflush(stdout);
  LTAfree(&gtmopt->pvfseg2pet);

  ltaArray[0] = gtmopt->gtmsegmu2anat;
  ltaArray[1] = gtmopt->anat2pet;
  gtmopt->gtmsegmu2pet = LTAconcat(ltaArray,2,1);
  if(gtmopt->gtm->gtmseg) MRIfree(&gtmopt->gtm->gtmseg);
  gtmopt->gtm->gtmseg = 
    MRIaseg2volMU(gtmopt->gtmsegmu, gtmopt->gtmsegmu2pet, 0.0, &hitvol, -1, gtmopt->gtm->ctGTMSeg);
  if(gtmopt->gtm->gtmseg == NULL) return(-1);
  //printf("  aseg2vol t = %g\n",TimerStop(&timer)/1000.0);fflush(stdout);
  LTAfree(&gtmopt->gtmsegmu2pet);

  GTMsegidlist(gtmopt->gtm);
  GTMdilateSeg(gtm);
  GTMbuildX(gtm);
  //printf("  buildX t = %g\n",TimerStop(&timer)/1000.0);fflush(stdout);
  if(gtm->X==NULL) return(-1);
  GTMsolve(gtm);

  gtmopt->tLastEval = TimerStop(&timer)/1000.0;
  return(0);
}
/*------------------------------------------------------------------*/
GTM *GTMalloc()
{
  GTM *gtm;
  gtm = (GTM *) calloc(sizeof(GTM),1);
  gtm->PadThresh = .0001;
  return(gtm);
}
/*------------------------------------------------------------------*/
int GTMfree(GTM **pGTM)
{
  int tt,nTT;
  GTM *gtm = *pGTM;

  nTT = gtm->ctGTMSeg->ctabTissueType->nentries-1;
  for(tt=0; tt < nTT; tt++){
    MRIfree(&gtm->pvf[tt]);
    MRIfree(&gtm->segttdil[tt]);
  }
  free(gtm->pvf);
  MRIfree(&gtm->yvol);
  MRIfree(&gtm->gtmseg);
  MRIfree(&gtm->mask);
  MatrixFree(&gtm->X);
  MatrixFree(&gtm->y);
  MatrixFree(&gtm->Xt);
  MatrixFree(&gtm->XtX);
  MatrixFree(&gtm->iXtX);
  MatrixFree(&gtm->Xty);
  MatrixFree(&gtm->beta);
  MatrixFree(&gtm->res);
  MatrixFree(&gtm->yhat);
  MRIfree(&gtm->ysynth);
  free(gtm);
  *pGTM=NULL;
  return(0);
}
/*------------------------------------------------------------------*/
int GTMmatrixY(GTM *gtm)
{
  gtm->y = GTMvol2mat(gtm, gtm->yvol, NULL);
  return(0);
}
/*------------------------------------------------------------------*/
int GTMsetNMask(GTM *gtm)
{
  if(gtm->mask) gtm->nmask = MRIcountAboveThreshold(gtm->mask, 0.5);
  else gtm->nmask = gtm->yvol->width*gtm->yvol->height*gtm->yvol->depth;
  return(0);
}
/*------------------------------------------------------------------*/
int GTMpsfStd(GTM *gtm)
{
  gtm->cStd = gtm->cFWHM/sqrt(log(256.0));
  gtm->rStd = gtm->rFWHM/sqrt(log(256.0));
  gtm->sStd = gtm->sFWHM/sqrt(log(256.0));
  return(0);
}
/*------------------------------------------------------------------*/
int GTMsegidlist(GTM *gtm)
{
  int *segidlist0,msegs,nthseg;
  segidlist0 = MRIsegIdList(gtm->anatseg, &gtm->nsegs, 0);
  // remove 0 from the list
  gtm->segidlist = (int *)calloc(sizeof(int),gtm->nsegs);
  msegs = 0;
  for(nthseg = 0; nthseg < gtm->nsegs; nthseg++){
    if(segidlist0[nthseg] != 0){
      gtm->segidlist[msegs] = segidlist0[nthseg];
      msegs++;
    }
  }
  gtm->nsegs = msegs;
  free(segidlist0);
  return(0);
}
/*------------------------------------------------------------------*/
int GTMnPad(GTM *gtm)
{
  double maxFWHM, maxStd;
  maxFWHM = MAX(gtm->cFWHM/gtm->gtmseg->xsize,
		MAX(gtm->rFWHM/gtm->gtmseg->ysize,gtm->sFWHM/gtm->gtmseg->zsize));
  maxStd = maxFWHM*sqrt(log(256.0));
  gtm->nPad = ceil(sqrt(-log(gtm->PadThresh*maxStd*sqrt(2*M_PI))*2*maxStd));
  if(Gdiag_no > 0) printf("maxFWHM = %g (voxels), PadThresh=%g, nPad=%d\n",maxFWHM,gtm->PadThresh,gtm->nPad);
  return(0);
}
/*------------------------------------------------------------------*/
int GTMsolve(GTM *gtm)
{
  struct timeb timer;
  int n,f;
  double sum;

  gtm->Xt = MatrixTranspose(gtm->X,gtm->Xt);
  printf("Computing  XtX ... ");fflush(stdout);
  TimerStart(&timer);
  gtm->XtX = MatrixMtM(gtm->X,gtm->XtX);
  printf(" %4.1f sec\n",TimerStop(&timer)/1000.0);fflush(stdout);
  gtm->iXtX = MatrixInverse(gtm->XtX,gtm->iXtX);
  if(gtm->iXtX==NULL){
    gtm->XtXcond = MatrixConditionNumber(gtm->XtX);
    printf("ERROR: matrix cannot be inverted, cond=%g\n",gtm->XtXcond);
    return(1);
  }
  gtm->Xty = MatrixMultiplyD(gtm->Xt,gtm->y,gtm->Xty);
  gtm->beta = MatrixMultiplyD(gtm->iXtX,gtm->Xty,gtm->beta);
  gtm->yhat = MatrixMultiplyD(gtm->X,gtm->beta,gtm->yhat);
  gtm->res  = MatrixSubtract(gtm->y,gtm->yhat,gtm->res);
  gtm->dof = gtm->X->rows - gtm->X->cols;
  if(gtm->rvar==NULL) gtm->rvar = MatrixAlloc(1,gtm->res->cols,MATRIX_REAL);
  for(f=0; f < gtm->res->cols; f++){
    sum = 0;
    for(n=0; n < gtm->res->rows; n++) sum += (gtm->res->rptr[n+1][f+1]*gtm->res->rptr[n+1][f+1]);
    gtm->rvar->rptr[1][f+1] = sum/gtm->dof;
  }
  return(0);
}
/*-----------------------------------------------------------------*/
MRI *GTMmat2vol(GTM *gtm, MATRIX *m, MRI *vol)
{
  int k,c,r,s,f;

  if(vol == NULL){
    vol = MRIallocSequence(gtm->yvol->width,gtm->yvol->height,gtm->yvol->depth,MRI_FLOAT,m->cols);
    if(vol == NULL) return(NULL);
    MRIcopyHeader(gtm->yvol,vol);
    MRIcopyPulseParameters(gtm->yvol,vol);
  }
  if(MRIdimMismatch(gtm->yvol,vol,0)){
    printf("ERROR: GTMmat2vol() dim mismatch\n");
    return(NULL);
  }

  k = 0;
  for(s=0; s < gtm->yvol->depth; s++){
    for(c=0; c < gtm->yvol->width; c++){
      for(r=0; r < gtm->yvol->height; r++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < gtm->yvol->nframes; f++)
	  MRIsetVoxVal(vol,c,r,s,f,m->rptr[k+1][f+1]);
	k++;
      }
    }
  }
  return(vol);
}
/*-----------------------------------------------------------------*/
MATRIX *GTMvol2mat(GTM *gtm, MRI *vol, MATRIX *m)
{
  int k,c,r,s,f;
  
  if(m==NULL){
    m = MatrixAlloc(gtm->nmask,vol->nframes,MATRIX_REAL);
    if(m==NULL){
      printf("ERROR: GTMvol2mat(): could not alloc matrix %d %d\n",gtm->nmask,vol->nframes);
      return(NULL);
    }
  }
  if(m->rows != gtm->nmask){
    printf("ERROR: GTMvol2mat(): row mismatch %d %d\n",m->rows,gtm->nmask);
    return(NULL);
  }
  if(m->cols != vol->nframes){
    printf("ERROR: GTMvol2mat(): col mismatch %d %d\n",m->cols,vol->nframes);
    return(NULL);
  }

  k = 0;
  for(s=0; s < vol->depth; s++){
    for(c=0; c < vol->width; c++){
      for(r=0; r < vol->height; r++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < vol->nframes; f++)
	  m->rptr[k+1][f+1] = MRIgetVoxVal(vol,c,r,s,f);
	k++;
      }
    }
  }
  return(m);
}
/*-----------------------------------------------------*/
int GTMsegrvar(GTM *gtm)
{
  int c,r,s,f,k;
  int nthseg=0,segid;
  double v;

  gtm->segstats = MatrixAlloc(gtm->nsegs,gtm->beta->cols,MATRIX_REAL);
  gtm->nperseg = (int *)calloc(sizeof(int),gtm->nsegs);

  k=0;
  for(s=0; s < gtm->yvol->depth; s++){
    for(c=0; c < gtm->yvol->width; c++){
      for(r=0; r < gtm->yvol->height; r++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	segid = MRIgetVoxVal(gtm->gtmseg,c,r,s,0);
	if(segid != 0){
	  for(nthseg=0; nthseg < gtm->nsegs; nthseg++) if(gtm->segidlist[nthseg]==segid) break;
	  gtm->nperseg[nthseg] ++;
	}
	for(f=0; f < gtm->beta->cols; f++){
	  v = gtm->res->rptr[k+1][f+1];
	  if(segid != 0) gtm->segstats->rptr[nthseg+1][f+1] += v*v;
	}
	k++;
      }// r 
    }// c
  }// s

  for(f=0; f < gtm->beta->cols; f++) {
    for(nthseg=0; nthseg < gtm->nsegs; nthseg++){
      v = gtm->segstats->rptr[nthseg+1][f+1];
      gtm->segstats->rptr[nthseg+1][f+1] = v/gtm->nperseg[nthseg];
    }
  }
  return(0);
}

/*--------------------------------------------------------------------------*/
int GTMrbv(GTM *gtm)
{
  int c,r,s,f;
  double val;

  if(gtm->rbv) MRIfree(&gtm->rbv);
  gtm->rbv = MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth,
			      MRI_FLOAT, gtm->yvol->nframes);
  MRIcopyHeader(gtm->yvol,gtm->rbv);
  for(s=0; s < gtm->yvol->depth; s++){
    for(c=0; c < gtm->yvol->width; c++){
      for(r=0; r < gtm->yvol->height; r++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < gtm->yvol->nframes; f++){
	  val = (double)MRIgetVoxVal(gtm->yvol,c,r,s,f)*
	    MRIgetVoxVal(gtm->ysynth,c,r,s,f)/
	      (MRIgetVoxVal(gtm->ysynthsm,c,r,s,f)+FLT_EPSILON);
	    MRIsetVoxVal(gtm->rbv,c,r,s,f,val);
	  }
	}
      }
    }
  return(0);
}

/*--------------------------------------------------------------------------*/
int GTMsmoothSynth(GTM *gtm)
{
  if(gtm->ysynth == NULL) GTMsynth(gtm);
  gtm->ysynthsm = MRIgaussianSmoothNI(gtm->ysynth, gtm->cStd, gtm->rStd, gtm->sStd, gtm->ysynthsm);
  return(0);
}

/*--------------------------------------------------------------------------*/
int GTMmgpvc(GTM *gtm)
{
  int c,r,s,f,n,nthseg,segid;
  double sum,vgmpsf,vwmpsf,vwmtac,vtac,vmgtac;
  MRI *gmpvf,*gmpvfpsf,*wmpvfpsf;

  gtm->mg_reftac = MatrixAlloc(gtm->yvol->nframes,1,MATRIX_REAL);

  for(f=0; f < gtm->yvol->nframes; f++){
    sum = 0;
    for(n=0; n < gtm->n_mg_refids; n++){
      segid = -1;
      for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
	segid = gtm->segidlist[nthseg];
	if(segid == gtm->mg_refids[n]) break;
      }
      if(segid == -1){
	printf("ERROR: could not find a match for n=%d %d\n",n,gtm->mg_refids[n]);
	return(1);
      }
      sum += gtm->beta->rptr[nthseg+1][f+1];
      printf("n=%d, nthseg=%d %g\n",n,nthseg,gtm->beta->rptr[nthseg+1][f+1]);
    }
    gtm->mg_reftac->rptr[f+1][1] = sum/gtm->n_mg_refids;
    printf("wm tac %2d %g\n",f,gtm->mg_reftac->rptr[f+1][1]);
  }

  if(gtm->mg) MRIfree(&gtm->mg);
  gtm->mg = MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth,
			   MRI_FLOAT, gtm->yvol->nframes);
  if(gtm->mg == NULL) return(1);
  MRIcopyHeader(gtm->yvol,gtm->mg);

  gmpvf = MRIadd(gtm->pvf[0],gtm->pvf[1],NULL);

  gmpvfpsf = MRIgaussianSmoothNI(gmpvf,gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  wmpvfpsf = MRIgaussianSmoothNI(gtm->pvf[2],gtm->cStd, gtm->rStd, gtm->sStd, NULL);

  for(c=0; c < gtm->yvol->width; c++){
    for(r=0; r < gtm->yvol->height; r++){
      for(s=0; s < gtm->yvol->depth; s++){
	if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue; 
	vgmpsf = MRIgetVoxVal(gmpvfpsf,c,r,s,0);
	if(vgmpsf < gtm->mg_gmthresh) continue; 
	vwmpsf = MRIgetVoxVal(wmpvfpsf,c,r,s,0);
	for(f=0; f < gtm->yvol->nframes; f++){
	  vwmtac = gtm->mg_reftac->rptr[f+1][1];
	  vtac = MRIgetVoxVal(gtm->yvol,c,r,s,f);
	  vmgtac = (vtac - vwmpsf*vwmtac)/vgmpsf;
	  MRIsetVoxVal(gtm->mg,c,r,s,f, vmgtac);
	}
      }
    }
  }
  MRIfree(&gmpvf);
  MRIfree(&gmpvfpsf);
  MRIfree(&wmpvfpsf);
  return(0);
}
/*---------------------------------------------------------*/
int MinPowell()
{
  extern GTMOPT *gtmopt_powell;
  GTMOPT *gtmopt = gtmopt_powell;
  float *pPowel, **xi;
  int    r, c, n,dof;
  struct timeb timer;

  TimerStart(&timer);
  dof = gtmopt->nparams;

  printf("Init Powel Params dof = %d\n",dof);
  pPowel = vector(1, dof) ;
  for(n=0; n < dof; n++) pPowel[n+1] = gtmopt->params[n];
  //pPowel[1] =  .0897;
  //pPowel[2] =  .2102;
  //pPowel[3] = -.2905;
  //pPowel[4] = 0;
  //pPowel[5] = 0;
  //pPowel[6] = 0;

  xi = matrix(1, dof, 1, dof) ;
  for (r = 1 ; r <= dof ; r++) {
    for (c = 1 ; c <= dof ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  OpenPowell2(pPowel, xi, dof, gtmopt->ftol, gtmopt->linmintol, gtmopt->nitersmax, 
	      &gtmopt->niters, &gtmopt->fret, compute_powell_cost);
  printf("Powell done niters = %d\n",gtmopt->niters);
  printf("OptTimeSec %4.1f sec\n",TimerStop(&timer)/1000.0);
  printf("OptTimeMin %5.2f sec\n",(TimerStop(&timer)/1000.0)/60);
  printf("nEvals %d\n",gtmopt->nCostEvaluations);
  printf("EvalTimeSec %4.1f sec\n",(TimerStop(&timer)/1000.0)/gtmopt->nCostEvaluations);
  fflush(stdout);
  for(n=0; n < dof; n++) gtmopt->params[n] = pPowel[n+1];
  gtmopt->anat2pet = LTAcopy(gtmopt->anat2pet0,gtmopt->anat2pet);
  LTAapplyAffineParametersTKR(gtmopt->anat2pet, gtmopt->params, gtmopt->nparams, gtmopt->anat2pet);
  printf("#@# %4d  ",gtmopt->nCostEvaluations);
  for(n=0; n<gtmopt->nparams; n++) printf("%10.6f ",gtmopt->params[n]);
  printf("\n");

  free_matrix(xi, 1, dof, 1, dof);
  free_vector(pPowel, 1, dof);
  return(NO_ERROR) ;
}
/*--------------------------------------------------------------------------*/
float compute_powell_cost(float *pPowel) 
{
  extern GTMOPT *gtmopt_powell;
  GTMOPT *gtmopt = gtmopt_powell;
  int n,newmin;
  float pp[100],curcost;
  static float initcost=-1,mincost=-1,ppmin[100];
  FILE *fp=NULL;
  
  for(n=0; n<gtmopt->nparams; n++) pp[n] = pPowel[n+1];

  gtmopt->anat2pet = LTAcopy(gtmopt->anat2pet0,gtmopt->anat2pet);
  LTAapplyAffineParametersTKR(gtmopt->anat2pet, pp, gtmopt->nparams, gtmopt->anat2pet);
  GTMOPTcost(gtmopt);
  curcost = gtmopt->gtm->rvar->rptr[1][1];
  newmin = 0;
  if(initcost<0) {
    newmin = 1;
    initcost = curcost;
    mincost = curcost;
    for(n=0; n<gtmopt->nparams; n++) ppmin[n] = pp[n];
    printf("InitialCost %f\n",initcost);
    fp=fopen("costs.dat","w");
  }
  else fp=fopen("costs.dat","a");

  if(mincost > curcost) {
    newmin = 1;
    mincost = curcost;
    for(n=0; n<gtmopt->nparams; n++) ppmin[n] = pp[n];
  }

  fprintf(fp,"%4d  ",gtmopt->nCostEvaluations);
  for(n=0; n<gtmopt->nparams; n++) fprintf(fp,"%12.8f ",pp[n]);
  fprintf(fp,"  %12.10f \n",curcost/initcost);
  fflush(fp);
  fclose(fp);

  if(newmin){
    printf("#@# %4d  ",gtmopt->nCostEvaluations);
    for(n=0; n<gtmopt->nparams; n++) printf("%7.3f ",ppmin[n]);
    printf("  %5.1f %8.6f\n",gtmopt->tLastEval,mincost/initcost);
    fflush(stdout);
  }

  gtmopt->nCostEvaluations++;
  return((float)gtmopt->gtm->rvar->rptr[1][1]);
}

/*--------------------------------------------------------------------------*/
LTA *LTAapplyAffineParametersTKR(LTA *inlta, const float *p, const int np, LTA *outlta)
{
  static MATRIX *M;
  static double angles[3];
  MATRIX *R;
  int intype;

  intype = inlta->type;
  outlta = LTAcopy(inlta,outlta);
  outlta = LTAchangeType(outlta,REGISTER_DAT);
  R = outlta->xforms[0].m_L;

  // New R = Mshear*Mscale*Mtrans*Mrot*R (consistent with bbr)
  if(np >= 6){
    angles[0] = p[3]*(M_PI/180);
    angles[1] = p[4]*(M_PI/180);
    angles[2] = p[5]*(M_PI/180);
    M = MRIangles2RotMat(angles);
    MatrixMultiply(M,R,R);    
    MatrixFree(&M);
  }
  if(np >= 3){
    // shift
    M = MatrixIdentity(4,M);
    M->rptr[1][4] = p[0];
    M->rptr[2][4] = p[1];
    M->rptr[3][4] = p[2];
    MatrixMultiply(M,R,R);    
  }
  if(np >= 9){
    // scale
    if(p[6]==0 || p[7]==0 || p[8]==0){
      printf("ERROR: LTAapplyAffineParametersTKR(): 0 scale %f %f %f\n",p[6],p[7],p[8]);
      outlta->xforms[0].m_L = NULL;
      return(NULL);
    }
    M = MatrixIdentity(4,M);
    M->rptr[1][1] = p[6];
    M->rptr[2][2] = p[7];
    M->rptr[3][3] = p[8];
    MatrixMultiply(M,R,R);    
  }
  if(np >= 12){
    // shear
    M = MatrixIdentity(4,M);
    M->rptr[1][2] = p[9];
    M->rptr[1][3] = p[10];
    M->rptr[2][3] = p[11];
    MatrixMultiply(M,R,R);    
  }

  if(0){
    int n;
    printf("---------------------------\n");
    printf("Parameters: ");
    for(n=0; n<np; n++) printf("%6.4f ",p[n]);
    printf("\n");fflush(stdout);
    MatrixPrint(stdout,R);fflush(stdout);
    printf("---------------------------\n");
    fflush(stdout);
  }

  outlta = LTAchangeType(outlta,intype);
  return(outlta);
}
/*------------------------------------------------------------------*/
int GTMsynth(GTM *gtm)
{
  MATRIX *yhat;

  if(gtm->ysynth) MRIfree(&gtm->ysynth);
  gtm->ysynth = MRIallocSequence(gtm->gtmseg->width,gtm->gtmseg->height,gtm->gtmseg->depth,MRI_FLOAT,
				 gtm->beta->cols);
  if(gtm->yvol) {
    MRIcopyHeader(gtm->yvol,gtm->ysynth);
    MRIcopyPulseParameters(gtm->yvol,gtm->ysynth);
  }
  yhat = MatrixMultiply(gtm->X0,gtm->beta,NULL);
  GTMmat2vol(gtm, yhat, gtm->ysynth);
  MatrixFree(&yhat);
    
  return(0);
}
/*------------------------------------------------------------*/
int CheckX(MATRIX *X)
{
  int r,c,count;
  double sum,d,dmax;

  count=0;
  dmax = -1;
  for(r=0; r < X->rows; r++){
    sum = 0;
    for(c=0; c < X->cols; c++){
      sum += X->rptr[r+1][c+1];
    }
    d = abs(sum-1);
    if(d>.00001) count++;
    if(dmax < d) dmax = d;
  }

  printf("CheckX: count=%d, dmax=%g\n",count,dmax);
  return(count);
}



int dngtest(LTA *aseg2vol)
{
  MRI *seg,*segf=NULL,*aseg,*mask,*pvf,*seg2,*newseg;
  char *subject, *SUBJECTS_DIR;
  char tmpstr[5000];
  int f,nsegs, *segidlist; //nPad=2,USF=1;
  LTA *seg2vol,*lta,*aseg2hrseg,*ltaArray[2];
  //LTA *anat2seg, *seg2anat,*ltaArray[2];
  COLOR_TABLE *ct;
  char *wsurf="white", *psurf="pial";
  MRIS *lhw, *lhp, *rhw, *rhp;

#ifdef _OPENMP
  omp_set_num_threads(8);
#endif

  ct = TissueTypeSchema(NULL,"default-jan-2014+head");

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  SUBJECTS_DIR = "/autofs/cluster/con_009/users/greve/fdg-pvc/FSMR";
  subject = aseg2vol->subject;

  mask = MRIread("mask.nii.gz");

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,"aparc+aseg.mgz");
  //sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,"petseg+head.mgz");
  //sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,"petseg.mgz");
  //sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,"aseg.mgz");
  printf("Loading %s\n",tmpstr);
  aseg = MRIread(tmpstr);
  if(aseg==NULL) exit(1);
  if(aseg->type != MRI_INT) aseg = MRIchangeType(aseg, MRI_INT, 0,0,1);
  aseg = MRIaddExtraCerebralCSF(aseg, 3, aseg);

  sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,wsurf);
  lhw = MRISread(tmpstr);
  if(lhw==NULL) exit(1);
  sprintf(tmpstr,"%s/%s/label/lh.aparc.annot",SUBJECTS_DIR,subject);
  MRISreadAnnotation(lhw, tmpstr);

  sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,psurf);
  lhp = MRISread(tmpstr);
  if(lhp==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,wsurf);
  rhw = MRISread(tmpstr);
  if(rhw==NULL) exit(1);
  sprintf(tmpstr,"%s/%s/label/rh.aparc.annot",SUBJECTS_DIR,subject);
  MRISreadAnnotation(rhw, tmpstr);

  sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,psurf);
  rhp = MRISread(tmpstr);
  if(rhp==NULL) exit(1);

  printf("computing hires seg\n");
  newseg = MRIhiresSeg(aseg,lhw,lhp,rhw,rhp,1,1,&aseg2hrseg);
  ltaArray[0] = LTAinvert(aseg2hrseg,NULL);
  ltaArray[1] = aseg2vol;
  lta = LTAconcat(ltaArray,2,1);

  seg = newseg;
  seg2vol = lta;
  segidlist = MRIsegIdListNot0(seg, &nsegs, 0);

  printf("computing seg with MRIseg2SegPVF\n");fflush(stdout);
  seg2 = MRIseg2SegPVF(seg, seg2vol, 0.5, segidlist, nsegs, mask, 1, ct, NULL);
  MRIwrite(seg2,"dng.apas.nii");
  exit(0);

  printf("aseg2volframe %d\n",nsegs);
  //for(f=0; f < nsegs; f++) printf("%2d %4d\n",f,segidlist[f]);
  for(f=0; f < 2; f++){
    printf("f = %d -------------------------------------------\n",f);
    segf = MRIseg2SegPVF(seg, seg2vol, 0.5, segidlist, nsegs, mask, 0, NULL, segf);
    //sprintf(tmpstr,"segf.%02d.mgh",f); MRIwrite(segf,tmpstr);
    printf("\n");
  }

  printf("Computing seg with MRIsegPVF2Seg\n");fflush(stdout);
  seg2 = MRIsegPVF2Seg(segf, segidlist, nsegs, ct, mask, NULL);
  MRIwrite(seg2,"dng.seg.nii");


  pvf=NULL;
  printf("saving pvf\n");fflush(stdout);
  pvf = MRIsegPVF2TissueTypePVF(segf, segidlist, nsegs, ct, mask, NULL);
  MRIwrite(pvf,"dng.pvf.nii");

  //printf("saving\n");
  //MRIwrite(segf,"segf.nii");
  printf("done\n");
  return(0);

  sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,wsurf);
  lhw = MRISread(tmpstr);
  if(lhw==NULL) exit(1);
  sprintf(tmpstr,"%s/%s/label/lh.aparc.annot",SUBJECTS_DIR,subject);
  MRISreadAnnotation(lhw, tmpstr);

  sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,psurf);
  lhp = MRISread(tmpstr);
  if(lhp==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,wsurf);
  rhw = MRISread(tmpstr);
  if(rhw==NULL) exit(1);
  sprintf(tmpstr,"%s/%s/label/rh.aparc.annot",SUBJECTS_DIR,subject);
  MRISreadAnnotation(rhw, tmpstr);

  sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,psurf);
  rhp = MRISread(tmpstr);
  if(rhp==NULL) exit(1);

  printf("computing hires seg\n");
  newseg = MRIhiresSeg(aseg,lhw,lhp,rhw,rhp,1,1,&lta);
  MRIwrite(newseg,"newseg.annot.usf1.mgh");

  printf("computing hires seg\n");
  newseg = MRIhiresSeg(aseg,lhw,lhp,rhw,rhp,1,2,&lta);
  MRIwrite(newseg,"newseg.annot.usf2.mgh");

  exit(1);


}

int GTMbuildX(GTM *gtm)
{
  int nthseg;
  struct timeb timer;

  if(gtm->X==NULL || gtm->X->rows != gtm->nmask || gtm->X->cols != gtm->nsegs){
    if(gtm->X) MatrixFree(&gtm->X);
    gtm->X = MatrixAlloc(gtm->nmask,gtm->nsegs,MATRIX_REAL);
    if(gtm->X == NULL){
      printf("ERROR: GTMbuildX(): could not alloc X %d %d\n",gtm->nmask,gtm->nsegs);
      return(1);
    }
  }
  if(gtm->X0==NULL || gtm->X0->rows != gtm->nmask || gtm->X0->cols != gtm->nsegs){
    if(gtm->X0) MatrixFree(&gtm->X0);
    gtm->X0 = MatrixAlloc(gtm->nmask,gtm->nsegs,MATRIX_REAL);
    if(gtm->X0 == NULL){
      printf("ERROR: GTMbuildX(): could not alloc X0 %d %d\n",gtm->nmask,gtm->nsegs);
      return(1);
    }
  }
  gtm->dof = gtm->X->rows - gtm->X->cols;

  TimerStart(&timer);
  printf("computing seg pvf \n");fflush(stdout);
  gtm->segpvf = MRIseg2SegPVF(gtm->anatseg, gtm->anatseg2pet, 0.5, gtm->segidlist, 
			       gtm->nsegs, gtm->mask, 0, NULL, gtm->segpvf);
  if(gtm->segpvf==NULL) return(1);


  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
  for(nthseg = 0; nthseg < gtm->nsegs; nthseg++){
    int segid,k,c,r,s;
    MRI *nthsegpvf=NULL,*nthsegpvfbb=NULL,*nthsegpvfbbsm=NULL;
    MRI_REGION *region;
    segid = gtm->segidlist[nthseg];
    //printf("nthseg = %d, %d %6.4f\n",nthseg,segid,TimerStop(&timer)/1000.0);fflush(stdout);
    nthsegpvf = fMRIframe(gtm->segpvf,nthseg,NULL);
    region      = REGIONgetBoundingBox(nthsegpvf,gtm->nPad);
    nthsegpvfbb = MRIextractRegion(nthsegpvf, NULL, region) ;
    nthsegpvfbbsm = MRIgaussianSmoothNI(nthsegpvfbb, gtm->cStd, gtm->rStd, gtm->sStd, nthsegpvfbbsm);
    // Fill X, creating X in this order makes it consistent with matlab
    // Note: y must be ordered in the same way.
    k = 0;
    for(s=0; s < gtm->gtmseg->depth; s++){
      for(c=0; c < gtm->gtmseg->width; c++){
	for(r=0; r < gtm->gtmseg->height; r++){
	  if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	  k ++;
	  if(c < region->x || c >= region->x+region->dx)  continue;
	  if(r < region->y || r >= region->y+region->dy)  continue;
	  if(s < region->z || s >= region->z+region->dz)  continue;
	  // do not use k+1 here because it has already been incr above
	  gtm->X->rptr[k][nthseg+1] = 
	    MRIgetVoxVal(nthsegpvfbbsm,c-region->x,r-region->y,s-region->z,0);
	  gtm->X0->rptr[k][nthseg+1] = 
	    MRIgetVoxVal(nthsegpvfbb,c-region->x,r-region->y,s-region->z,0);
	}
      }
    }
    MRIfree(&nthsegpvf);
    MRIfree(&nthsegpvfbb);
    MRIfree(&nthsegpvfbbsm);
  }
  printf(" Build time %6.4f\n",TimerStop(&timer)/1000.0);fflush(stdout);

  return(0);

}

