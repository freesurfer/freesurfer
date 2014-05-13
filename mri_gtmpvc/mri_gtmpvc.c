/**
 * @file  mri_gtmpvc.c
 * @brief Peforms Partial volume correction based on the geometric transfer matrix GTM
 *
 *
 * Implementation of geometric transfer matrix (GTM). Also includes Muller-Gartner (MG) and
 * Region-based Voxelwise (RBV) partial volume correction.
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2014/05/13 21:19:28 $
 *    $Revision: 1.11 $
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


// $Id: mri_gtmpvc.c,v 1.11 2014/05/13 21:19:28 greve Exp $

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
#include "registerio.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "cma.h"
#include "mri_identify.h"
#include "gtm.h"

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

static char vcid[] = "$Id: mri_gtmpvc.c,v 1.11 2014/05/13 21:19:28 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

typedef struct 
{
  char *name; // contrast name
  MATRIX *C;
  MATRIX *gamma;
  MATRIX *gammavar;
  MATRIX *t;
  MATRIX *p;
} GTM_CONTRAST, GTMCON;


typedef struct 
{
  MRI *yvol; // source (PET) data
  MRI *anatseg; // segmentation in anatomical space, may be hires (ie, not conformed)
  LTA *anat2pet; // registration from anat to pet space
  LTA *anat2seg; // registration from FS conformed space to seg
  LTA *seg2anat; // inverse registration
  LTA *seg2pet;  // registration from segmentation to source volume
  MRI *mask;     // binary mask in PET space
  double cFWHM, rFWHM, sFWHM; // assumed FWHM of PSF
  double cStd, rStd, sStd; // PSF FWHM converted to standard dev

  double automask_fwhm,automask_thresh; // Use this FWHM instead of PSF when computing mask
  int automask_reduce_fov; // Flag: reduce PET FoV to be tight to mask.
  MRI_REGION *automaskRegion; // Keep reduction region in case reinstate at original FoV
  MRI *yvol_full_fov; // Used to keep header of original FoV source data
  double PadThresh; // Used to dilate mask based on automask_fwhm
  int nPad; // Actual amount of dilation

  int rescale,n_scale_refids,scale_refids[100];
  double scale;

  MRI *segpvf; // PVF of each seg (one in each frame)
  MRI *ttpvf; // PVF of each tissue type (one in each frame), used by MG
  MRI *gtmseg; // Segmentation in PET space
  int nmask; // number of voxels in the mask

  int nsegs,*segidlist; // number of segments in segmentation, list of segids
  int *nperseg; // number of voxels per seg
  COLOR_TABLE *ctGTMSeg; // color table of segments
  int nReplace , SrcReplace[1000], TrgReplace[1000]; // for replacing segs
  MATRIX *nvox; // same as nperseg, doh
  MATRIX *vrf;  // variance reduction factor for each seg
  MATRIX *segrvar; // residual variance in each seg

  // GLM stuff for GTM
  MATRIX *X,*X0;
  MATRIX *y,*Xt, *XtX, *iXtX, *Xty, *beta, *res, *yhat,*betavar;
  MATRIX *rvar,*rvargm; // residual variance, all vox and only GM
  int dof;
  double XtXcond;
  MRI *ysynth,*ysynthsm; // synthesized vs yhat
  MATRIX *skew,*kurtosis; // for QAe
  int DoVoxFracCor; // Flag to correct for volume fraction effect

  int nContrasts;
  GTMCON *contrasts[100];

  MRI *rbv; // RBV computed volume
  int mask_rbv_to_brain; // Reduce FoV of RBV to be tight to brain
  MRI *yseg; // source volume trilin resampled to seg space (used with RBV)
  MRI *yhat0seg; // unsmoothed yhat created in seg space (used with RBV)
  MRI *yhatseg;  // smoothed yhat in seg space (used with RBV)
  MRI *rbvsegmean; // seg mean in RBV, used for QA

  int DoMGPVC; // Muller-Gartner
  MRI *mg; // MG output volume
  double mg_gmthresh; // GM PVF threshold 
  int n_mg_refids,mg_refids[100]; // WM reference seg IDs
  MATRIX *mg_reftac; // WM reference TAC

  int DoKMRef; // Kinetic Modeling Reference TAC
  int n_km_refids,km_refids[100]; // KM reference seg IDs
  MATRIX *km_reftac; // KM reference TAC

  int DoKMHB; // Kinetic Modeling HiBinding TAC
  int n_km_hbids,km_hbids[100];// KM HiBinding seg IDs
  MATRIX *km_hbtac; // KM HiBinding TAC

  char *OutDir; // output folder
  FILE *logfp;  // log file pointer
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
MRI *GTMsegSynth(GTM *gtm);
int GTMrbv(GTM *gtm);
int GTMmgpvc(GTM *gtm);
MATRIX *GTMvol2mat(GTM *gtm, MRI *vol, MATRIX *m);
MRI *GTMmat2vol(GTM *gtm, MATRIX *m, MRI *vol);
int GTMprintReplaceList(FILE *fp, const int nReplace, const int *ReplaceThis, const int *WithThat);
int GTMcheckReplaceList(const int nReplace, const int *ReplaceThis, const int *WithThat);
int GTMloadReplacmentList(const char *fname, int *nReplace, int *ReplaceThis, int *WithThat);
int GTMcheckX(MATRIX *X);
int GTMautoMask(GTM *gtm);
int GTMrvarGM(GTM *gtm);

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
char *OutBetaFile=NULL,*OutBetaVarFile=NULL,*OutXtXFile=NULL;
char *SrcBetaFile=NULL,*SynthFile=NULL;
int SynthOnly=0;
MATRIX *srcbeta;
double psfFWHM=-1;
char tmpstr[5000],logfile[5000];
char *PVFFile=NULL, *SegTTypeFile=NULL;
double ApplyFWHM=0;
char *Xfile=NULL,*X0file=NULL,*ymatfile=NULL,*betamatfile=NULL;
int SaveX=0, SaveX0=0, SaveYMat=0, SaveBetaMat=0;
char *VRFStatsFile=NULL;
char *eresFile=NULL, *yhatFile=NULL, *yhat0File=NULL,*yhatFullFoVFile=NULL;
char *OutSegFile=NULL;
MRI *mritmp;
char *RVarFile=NULL,*SkewFile=NULL,*KurtosisFile=NULL;
int RVarOnly=0;
int nthreads=1;
int ttReduce = 0;
char *MGPVCFile=NULL;
char *regfile;
int regidentity = 0;
int regtype;
int SaveEres=0, SaveYhat=0,SaveYhat0=0, SaveInput=0,SaveYhatFullFoV=0;

int VRFStats(GTM *gtm, double *vrfmean, double *vrfmin, double *vrfmax);
int WriteVRFStats(char *fname, GTM *gtm);


int GTMttest(GTM *gtm);
int GTMmgRefIds(GTM *gtm);
int GTMprintMGRefTAC(GTM *gtm, FILE *fp);
int GTMwriteMGRefTAC(GTM *gtm, char *filename);
int GTMrescale(GTM *gtm);
int GTMwriteContrasts(GTM *GTM);
int GTMprintRefIds(GTM *gtm, FILE *fp);
int GTMcheckRefIds(GTM *gtm);
int GTMrefTAC(GTM *gtm);

GTM *gtm;
GTMOPT *gtmopt;
GTMOPT *gtmopt_powell;

LTA *LTAapplyAffineParametersTKR(LTA *inlta, const float *p, const int np, LTA *outlta);
int DoOpt=0;
char *SUBJECTS_DIR;

int DoRBV=0;
int AutoMask=0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err,c,f,n;
  double vrfmin,vrfmax,vrfmean;
  struct timeb  mytimer, timer;
  MRI *gtmres;
  FILE *logfp,*fp;
  char *stem;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  vg_isEqual_Threshold = 10e-4;

  gtm = GTMalloc();
  //gtm->ctGTMSeg = TissueTypeSchema(NULL,"default-jan-2014+head");
  //gtm->ctGTMSeg = TissueTypeSchema(NULL,"default-jan-2014");

  // by default, rescale to Cerebellum WM
  gtm->rescale = 1; 
  gtm->n_scale_refids = 2;
  gtm->scale_refids[0] = 7;
  gtm->scale_refids[1] = 46;
  gtm->mask_rbv_to_brain = 1;
  gtm->nReplace = 0;
  gtm->nContrasts = 0;
  gtm->automask_reduce_fov = 1;
  gtm->DoVoxFracCor = 1;

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
    sprintf(logfile,"%s/mri_gtmpvc.log",OutDir);
    logfp = fopen(logfile,"w");
  }
  else logfp = stdout;
  gtm->logfp = logfp;
  dump_options(logfp);

  TimerStart(&timer);

  // Load seg
  TimerStart(&mytimer);
  printf("Loading seg for gtm %s\n",SegVolFile);fflush(stdout);
  fprintf(logfp,"Loading seg for gtm %s\n",SegVolFile);fflush(logfp);
  gtm->anatseg = MRIread(SegVolFile);
  if(gtm->anatseg==NULL) exit(1);

  stem = IDstemFromName(SegVolFile);
  sprintf(tmpstr,"%s.ctab",stem);
  printf("Loading seg ctab %s\n",tmpstr);fflush(stdout);
  fprintf(logfp,"Loading ctab %s\n",tmpstr);fflush(logfp);
  gtm->ctGTMSeg = CTABreadASCII(tmpstr);
  if(gtm->ctGTMSeg == NULL) exit(1);

  sprintf(tmpstr,"%s.lta",stem);
  if(fio_FileExistsReadable(tmpstr)){
    printf("Reading %s\n",tmpstr);
    fprintf(logfp,"Reading %s\n",tmpstr);
    gtm->anat2seg = LTAread(tmpstr);
    if(gtm->anat2seg == NULL) exit(1);
  }
  else{
    printf("Assuming seg is in conformed space\n");
    fprintf(logfp,"Assuming seg is in conformed space\n");
    gtm->anat2seg = TransformRegDat2LTA(gtm->anatseg, gtm->anatseg, NULL);
  }
  gtm->seg2anat = LTAinvert(gtm->anat2seg,NULL);
  gtm->seg2pet = LTAconcat2(gtm->seg2anat,gtm->anat2pet,1);
  if(gtm->seg2pet == NULL) {
    printf("ERROR: LTAconcat()\n");
    printf("mri_gtmpvc exited with errors\n");
    exit(1);
  }

  if(gtm->nReplace > 0) {
    printf("Replacing %d\n",gtm->nReplace);
    mritmp = MRIreplaceList(gtm->anatseg, gtm->SrcReplace, gtm->TrgReplace, gtm->nReplace, NULL);
    MRIfree(&gtm->anatseg);
    gtm->anatseg = mritmp;
    sprintf(tmpstr,"%s/seg.replace.list",OutDir);
    fp = fopen(tmpstr,"w");
    GTMprintReplaceList(fp, gtm->nReplace, gtm->SrcReplace, gtm->TrgReplace);
    fclose(fp);
  }
  if(CheckSegTissueType(gtm->anatseg, gtm->ctGTMSeg)){
    printf("Failed tissue type check\n");
    fprintf(logfp,"Failed tissue type check\n");
    exit(1);
  }
  gtm->ctGTMSeg = CTABpruneCTab(gtm->ctGTMSeg, gtm->anatseg);

  if(ttReduce > 0) {
    MRI *ttseg;
    COLOR_TABLE *ctTT;
    printf("Reducing seg to tissue type seg\n");
    fprintf(logfp,"Reducing seg to tissue type seg\n");
    ttseg = MRIseg2TissueType(gtm->anatseg, gtm->ctGTMSeg, NULL);
    if(ttseg == NULL) exit(1);
    MRIfree(&gtm->anatseg);
    gtm->anatseg = ttseg;
    ctTT = CTABdeepCopy(gtm->ctGTMSeg->ctabTissueType);
    for(f=0; f < ctTT->nentries; f++) ctTT->entries[f]->TissueType=f;
    ctTT->ctabTissueType = CTABdeepCopy(gtm->ctGTMSeg->ctabTissueType);
    CTABfree(&gtm->ctGTMSeg);
    gtm->ctGTMSeg = ctTT;
  }

  if(ApplyFWHM > 0){
    double cStdApply, rStdApply, sStdApply;
    MRI *mritmp;
    printf("Smoothing input by %g mm FWHM \n",ApplyFWHM);
    fprintf(logfp,"Smoothing input by %g mm FWHM \n",ApplyFWHM);
    cStdApply = ApplyFWHM/sqrt(log(256.0));
    rStdApply = ApplyFWHM/sqrt(log(256.0));
    sStdApply = ApplyFWHM/sqrt(log(256.0));
    mritmp = MRIgaussianSmoothNI(gtm->yvol, cStdApply, rStdApply, sStdApply, NULL);
    MRIfree(&gtm->yvol);
    gtm->yvol = mritmp;
  }

  // These have to be done before automask
  GTMpsfStd(gtm); 
  GTMnPad(gtm);   

  if(MaskVolFile){
    printf("Loading mask %s\n",MaskVolFile);fflush(stdout);
    fprintf(logfp,"Loading mask %s\n",MaskVolFile);
    gtm->mask = MRIread(MaskVolFile);
    if(gtm->mask==NULL) exit(1);
    err = MRIdimMismatch(gtm->yvol, gtm->mask, 0);
    if(err){
      printf("ERROR: seg and mask dim mismatch %d\n",err);
      exit(1);
    }
  }
  else if(AutoMask){
    printf("Computing auto mask \n");
    GTMautoMask(gtm);
    sprintf(tmpstr,"%s/mask.nii.gz",OutDir);
    MRIwrite(gtm->mask,tmpstr);
  }
  else  gtm->mask=NULL; // should already be NULL

  printf("Data load time %4.1f sec\n",TimerStop(&mytimer)/1000.0);
  fprintf(logfp,"Data load time %4.1f sec\n",TimerStop(&mytimer)/1000.0);
  fflush(stdout);fflush(logfp);

  GTMsetNMask(gtm);
  GTMsegidlist(gtm);
  GTMmatrixY(gtm);
  printf("nmask = %d, nsegs = %d, excluding segid=0\n",gtm->nmask,gtm->nsegs);
  printf("FWHM: %g %g %g\n",gtm->cFWHM,gtm->rFWHM,gtm->sFWHM);
  printf("Std:  %g %g %g\n",gtm->cStd,gtm->rStd,gtm->sStd);
  printf("nPad %d, PadThresh %g\n",gtm->nPad,gtm->PadThresh);
  fprintf(logfp,"nmask = %d, nsegs = %d, excluding segid=0\n",gtm->nmask,gtm->nsegs);
  fprintf(logfp,"FWHM: %g %g %g\n",gtm->cFWHM,gtm->rFWHM,gtm->sFWHM);
  fprintf(logfp,"Std:  %g %g %g\n",gtm->cStd,gtm->rStd,gtm->sStd);
  fprintf(logfp,"nPad %d, PadThresh %g\n",gtm->nPad,gtm->PadThresh);

  for(n=0; n < gtm->nContrasts; n++){
    if(gtm->contrasts[n]->C->cols != gtm->nsegs){
      printf("ERROR: contrast %d %s has %d columns, expecting %d\n",
	     n, gtm->contrasts[n]->name, gtm->contrasts[n]->C->cols,gtm->nsegs);
      fprintf(logfp,"ERROR: contrast %d %s has %d columns, expecting %d\n",
	     n, gtm->contrasts[n]->name, gtm->contrasts[n]->C->cols,gtm->nsegs);
      exit(1);
    }
  }
  sprintf(tmpstr,"%s/seg.ctab",OutDir);
  CTABwriteFileASCIItt(gtm->ctGTMSeg,tmpstr);

  printf("Checking Ref Ids\n");
  err=GTMcheckRefIds(gtm);
  if(err) exit(1);
  GTMprintRefIds(gtm, stdout);
  GTMprintRefIds(gtm, logfp);

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

  /* Create the "SegPVF". This is a multi-frame volume, each frame corresponds to a
     different Seg ID. The value is the PVF of that SegID. This is independent of
     the PSF (but accounts for the volume fraction effect). */
  printf("Computing Seg PVF \n");fflush(stdout);
  gtm->segpvf = MRIseg2SegPVF(gtm->anatseg, gtm->seg2pet, 0.5, gtm->segidlist, 
			      gtm->nsegs, gtm->mask, 0, NULL, gtm->segpvf);
  if(gtm->segpvf==NULL) exit(1);

  /* This creates a segmentation in the PET space based upon which Seg
     has the greated PVF (independent of PSF). */
  printf("Computing Seg in input space \n");fflush(stdout);
  gtm->gtmseg = MRIsegPVF2Seg(gtm->segpvf, gtm->segidlist, gtm->nsegs, 
			      gtm->ctGTMSeg, gtm->mask, gtm->gtmseg);

  // Create GTM matrix
  printf("Building GTM DoVoxFracCor=%d... ",gtm->DoVoxFracCor);fflush(stdout); 
  TimerStart(&mytimer) ;
  GTMbuildX(gtm);
  if(gtm->X==NULL) exit(1);
  printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  fprintf(logfp,"GTM-Build-time %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(logfp);

  // Create GTM seg in pet space (used by GTMsynth)
  gtm->gtmseg = MRIsegPVF2Seg(gtm->segpvf, gtm->segidlist, gtm->nsegs, 
			      gtm->ctGTMSeg, gtm->mask, gtm->gtmseg);
  if(OutSegFile){
    err=MRIwrite(gtm->gtmseg,OutSegFile);
    if(err) exit(1);
  }

  if(SrcBetaFile){
    printf("Synthsizing using supplied beta %s\n",SrcBetaFile);
    gtm->beta = srcbeta;
    GTMsynth(gtm);
    MRIwrite(gtm->ysynth,SynthFile);
    if(SynthOnly){
      printf("SynthOnly requested so exiting now\n");
      printf("mri_gtmpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
      exit(0);
    }
    MRIfree(&gtm->yvol);
    GTMsmoothSynth(gtm);
    gtm->yvol = MRIcopy(gtm->ysynthsm,NULL);
  }

  if(SaveX0) {
    printf("Writing X0 to %s\n",Xfile);
    MatlabWrite(gtm->X0, X0file,"X0");
  }
  if(SaveX) {
    printf("Writing X to %s\n",Xfile);
    MatlabWrite(gtm->X, Xfile,"X");
    //MatrixWriteTxt(Xfile, gtm->X);
  }
  printf("Solving ...\n");
  TimerStart(&mytimer) ; 
  err=GTMsolve(gtm); // also rescales everything if desired
  if(err) exit(1);
  printf("Time to solve %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  fprintf(logfp,"GTM-Solve-Time %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(logfp);

  if(SaveInput){
    sprintf(tmpstr,"%s/input.rescaled.nii.gz",OutDir);
    printf("Writing input to %s\n",tmpstr);
    MRIwrite(gtm->yvol,tmpstr);
  }

  if(ymatfile) {
    printf("Writing y to %s\n",ymatfile);
    MatlabWrite(gtm->y, ymatfile, "y");
    MatlabWrite(gtm->beta, "beta.mat", "beta");
  }
  if(betamatfile) {
    printf("Writing beta to %s\n",betamatfile);
    MatlabWrite(gtm->beta, betamatfile, "beta");
  }

  GTMrefTAC(gtm);

  sprintf(tmpstr,"%s/hrseg2pet.lta",OutDir);
  LTAwrite(gtm->seg2pet,tmpstr);
  sprintf(tmpstr,"%s/anat2pet.lta",OutDir);
  LTAwrite(gtm->anat2pet,tmpstr);
  //sprintf(tmpstr,"%s/anat2hrseg.lta",OutDir);
  //LTAwrite(gtm->anat2seg,tmpstr);


  if(gtm->rescale){
    // Rescaling is done during GTMsolve()
    printf("rescale factor %20.15lf\n",gtm->scale);
    fprintf(logfp,"rescale factor %20.15lf\n",gtm->scale);
    sprintf(tmpstr,"%s/scale.dat",OutDir);
    fp = fopen(tmpstr,"w");
    fprintf(fp,"%20.15lf\n",gtm->scale);
    fclose(fp);
  }

  GTMsegrvar(gtm);
  VRFStats(gtm, &vrfmean, &vrfmin, &vrfmax);

  // Create GTM pvf in pet space (why?)
  gtm->ttpvf = MRIsegPVF2TissueTypePVF(gtm->segpvf, gtm->segidlist, gtm->nsegs, 
				       gtm->ctGTMSeg, gtm->mask, gtm->ttpvf);
  sprintf(tmpstr,"%s/pvf.nii.gz",OutDir);
  MRIwrite(gtm->ttpvf,tmpstr);

  printf("Writing GTM beta estimates to %s\n",OutBetaFile);
  mritmp = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, gtm->yvol->nframes);
  for(c=0; c < gtm->nsegs; c++){
    for(f=0; f < gtm->yvol->nframes; f++)
      MRIsetVoxVal(mritmp,c,0,0,f, gtm->beta->rptr[c+1][f+1]);
  }
  err=MRIwrite(mritmp,OutBetaFile);
  if(err) exit(1);
  MRIfree(&mritmp);

  printf("Writing var of GTM estimates to %s\n",OutBetaVarFile);
  mritmp = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, gtm->yvol->nframes);
  for(c=0; c < gtm->nsegs; c++){
    for(f=0; f < gtm->yvol->nframes; f++)
      MRIsetVoxVal(mritmp,c,0,0,f, gtm->betavar->rptr[c+1][f+1]);
  }
  err=MRIwrite(mritmp,OutBetaVarFile);
  if(err) exit(1);
  MRIfree(&mritmp);

  printf("rvar = %g\n",gtm->rvar->rptr[1][1]);
  fprintf(logfp,"rvar = %g\n",gtm->rvar->rptr[1][1]);

  // Note: not going to rescale matrix because scaling should not 
  // be an issue here
  gtm->XtXcond = MatrixConditionNumber(gtm->XtX);
  printf("XtX  Condition     %8.3f \n",gtm->XtXcond);
  fprintf(logfp,"XtX  Condition     %8.3f \n",gtm->XtXcond);

  if(yhat0File || yhatFile || yhatFullFoVFile){
    printf("Synthesizing ... ");fflush(stdout); TimerStart(&mytimer) ;
    GTMsynth(gtm);
    printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  }
  if(yhat0File) MRIwrite(gtm->ysynth,yhat0File);

  if(yhatFile|| yhatFullFoVFile){
    printf("Smoothing synthesized ... ");fflush(stdout); TimerStart(&mytimer) ;
    GTMsmoothSynth(gtm);
    printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
    if(yhatFile) MRIwrite(gtm->ysynthsm,yhatFile);
  }

  if(yhatFullFoVFile){
    printf("Restoring FoV to yhat, saving in %s\n",yhatFullFoVFile);
    mritmp = MRIinsertRegion(gtm->ysynthsm, gtm->automaskRegion, gtm->yvol_full_fov, NULL);
    MRIwrite(mritmp,yhatFullFoVFile);
    MRIfree(&mritmp);
  }

  if(eresFile){
    gtmres = GTMmat2vol(gtm,gtm->res,NULL);
    MRIwrite(gtmres,eresFile);
    MRIfree(&gtmres);
  }

  fp = fopen(RVarFile,"w");
  for(f=0; f < gtm->yvol->nframes; f++)
    fprintf(fp,"%30.20f\n",gtm->rvar->rptr[1][f+1]);
  fclose(fp);
  if(RVarOnly){
    printf("rvar-only requested so exiting now\n");
    printf("mri_gtmpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
    exit(0);
  }
  GTMrvarGM(gtm);

  // Write the number of voxels in the mask
  sprintf(tmpstr,"%s/nmask.dat",OutDir);
  fp = fopen(tmpstr,"w");
  fprintf(fp,"%d\n",gtm->nmask);
  fclose(fp);

  if(KurtosisFile) {
    fp = fopen(KurtosisFile,"w");
    for(f=0; f < gtm->yvol->nframes; f++)
      fprintf(fp,"%30.20f\n",gtm->kurtosis->rptr[1][f+1]);
    fclose(fp);
  }
  if(SkewFile) {
    fp = fopen(SkewFile,"w");
    for(f=0; f < gtm->yvol->nframes; f++)
      fprintf(fp,"%30.20f\n",gtm->skew->rptr[1][f+1]);
    fclose(fp);
  }

  MatlabWrite(gtm->XtX, OutXtXFile,"XtX");

  printf("XtX  Condition     %8.3f \n",gtm->XtXcond);
  printf("VRF  Mean/Min/Max  %8.3f %8.3f %8.3f \n",vrfmean,vrfmin,vrfmax);
  fflush(stdout);
  fprintf(logfp,"XtX  Condition     %8.3f \n",gtm->XtXcond);
  fprintf(logfp,"VRF  Mean/Min/Max  %8.3f %8.3f %8.3f \n",vrfmean,vrfmin,vrfmax);
  fflush(logfp);
  WriteVRFStats(VRFStatsFile, gtm);

  // Write condition number to a dat file
  sprintf(tmpstr,"%s/cond.dat",OutDir);
  fp = fopen(tmpstr,"w");
  fprintf(fp,"%20.15lf\n",gtm->XtXcond);
  fclose(fp);

  // Res variance in each seg
  sprintf(tmpstr,"%s/seg.rvar.nii.gz",OutDir);
  printf("Writing seg rvar estimates to %s\n",tmpstr);
  mritmp = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, gtm->yvol->nframes);
  for(c=0; c < gtm->nsegs; c++){
    for(f=0; f < gtm->yvol->nframes; f++){
      MRIsetVoxVal(mritmp,c,0,0,f, gtm->segrvar->rptr[c+1][f+1]);
    }
  }
  err=MRIwrite(mritmp,tmpstr);
  if(err) exit(1);
  MRIfree(&mritmp);

  // VRF in each seg
  sprintf(tmpstr,"%s/seg.vrf.nii.gz",OutDir);
  printf("Writing seg vrf to %s\n",tmpstr);
  mritmp = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, 1);
  for(c=0; c < gtm->nsegs; c++) MRIsetVoxVal(mritmp,c,0,0,0, gtm->vrf->rptr[c+1][1]);
  err=MRIwrite(mritmp,tmpstr);
  if(err) exit(1);
  MRIfree(&mritmp);

  // NVox in each seg
  sprintf(tmpstr,"%s/seg.nvox.nii.gz",OutDir);
  printf("Writing seg nvox to %s\n",tmpstr);
  mritmp = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, 1);
  for(c=0; c < gtm->nsegs; c++) MRIsetVoxVal(mritmp,c,0,0,0, gtm->nvox->rptr[c+1][1]);
  err=MRIwrite(mritmp,tmpstr);
  if(err) exit(1);
  MRIfree(&mritmp);

  if(gtm->DoMGPVC){
    printf("Performing MG PVC\n");
    sprintf(tmpstr,"%s/mg.nii.gz",OutDir);
    MGPVCFile = strcpyalloc(tmpstr);
    fprintf(logfp,"MG PVC\n");
    GTMmgpvc(gtm);
    err = MRIwrite(gtm->mg,MGPVCFile);
    if(err) exit(1);
    printf("done with mgpvc\n");
    fprintf(logfp,"done with mgpvc\n");
    stem = IDstemFromName(MGPVCFile);
    sprintf(tmpstr,"%s.reftac.dat",stem);
    GTMwriteMGRefTAC(gtm, tmpstr);
  }
    
  if(DoRBV){
    sprintf(tmpstr,"%s/rbv.nii.gz",OutDir);
    RBVVolFile = strcpyalloc(tmpstr);
    printf("Computing RBV\n");
    GTMrbv(gtm);
    printf("Writing output to %s ...",RBVVolFile);fflush(stdout); TimerStart(&mytimer) ;
    err = MRIwrite(gtm->rbv,RBVVolFile);
    if(err){
      printf("ERROR: writing to %s\n",RBVVolFile);
      exit(1);
    }
    printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
    sprintf(tmpstr,"%s/rbv.segmean.nii.gz",OutDir);
    MRIwrite(gtm->rbvsegmean,tmpstr);
  }
  
  if(gtm->nContrasts > 0){
    printf("Testing %d contrasts\n",gtm->nContrasts);
    GTMttest(gtm);
    printf("Writing contrasts\n");
    GTMwriteContrasts(gtm);
  }

  fprintf(logfp,"mri_gtmpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
  fprintf(logfp,"mri_gtmpvc done\n");
  fclose(logfp);
  printf("mri_gtmpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
  printf("mri_gtmpvc done\n");
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
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if(!strcasecmp(option, "--help"))  print_help() ;
    else if(!strcasecmp(option, "--version")) print_version() ;
    else if(!strcasecmp(option, "--debug"))   debug = 1;
    else if(!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if(!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if(!strcasecmp(option, "--no-vox-frac-cor")) gtm->DoVoxFracCor=0;
    else if(!strcasecmp(option, "--no-vox-frac")) gtm->DoVoxFracCor=0;
    else if(!strcasecmp(option, "--no-vfc"))      gtm->DoVoxFracCor=0;
    else if(!strcasecmp(option, "--auto-mask")){
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%lf",&gtm->automask_fwhm);
      sscanf(pargv[1],"%lf",&gtm->automask_thresh);
      AutoMask=1;
      nargsused = 2;
    }
    else if(!strcasecmp(option, "--no-reduce-fov")) gtm->automask_reduce_fov = 0;
    else if(!strcasecmp(option, "--opt")) DoOpt=1;
    else if(!strcmp(option, "--sd") || !strcmp(option, "-SDIR")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      nargsused = 1;
    } 
    else if(!strcmp(option, "--ctab")){
      if(nargc < 1) CMDargNErr(option,1);
      gtm->ctGTMSeg = CTABreadASCII(pargv[0]);
      if(gtm->ctGTMSeg==NULL) exit(1);
      if(gtm->ctGTMSeg->ctabTissueType == NULL){
	printf("ERROR: ctab %s does not have a tissue type designation\n",pargv[0]);
	exit(1);
      }
      nargsused = 1;
    }

    else if(!strcasecmp(option, "--tt-reduce")) ttReduce = 1;
    else if(!strcasecmp(option, "--no-mask_rbv_to_brain")) gtm->mask_rbv_to_brain = 0;
    else if(!strcasecmp(option, "--default-seg-merge"))
      GTMdefaultSegReplacmentList(&gtm->nReplace,&(gtm->SrcReplace[0]),&(gtm->TrgReplace[0]));

    else if(!strcmp(option, "--replace-file")){
      if(nargc < 1) CMDargNErr(option,1);
      int err=GTMloadReplacmentList(pargv[0],&gtm->nReplace,&(gtm->SrcReplace[0]),&(gtm->TrgReplace[0]));
      if(err) exit(1);
      nargsused = 1;
    }

    else if(!strcmp(option, "--reg")){
      if(nargc < 1) CMDargNErr(option,1);
      regfile = pargv[0];
      regtype = TransformFileNameType(regfile);
      if(regtype == REGISTER_DAT){
	printf("ERROR: cannot use register.dat-style registration matrix, must be LTA\n");
	exit(1);
	//gtm->anatseg2pet = ltaReadRegisterDat(regfile, SrcVolFile, SegVolFile);
      }
      gtm->anat2pet = LTAread(regfile);
      if(gtm->anat2pet == NULL){
	printf("ERROR: reading %s\n",regfile);
	exit(1);
      }
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--reg-identity")) regidentity = 1;
    else if(!strcasecmp(option, "--identity")) regidentity = 1;
    else if(!strcasecmp(option, "--src") || !strcasecmp(option, "--i")) {
      if(nargc < 1) CMDargNErr(option,1);
      SrcVolFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--seg")) {
      if(nargc < 1) CMDargNErr(option,1);
      SegVolFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--mask")) {
      if(nargc < 1) CMDargNErr(option,1);
      MaskVolFile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--gdiag")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--psf")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&psfFWHM);
      gtm->cFWHM = psfFWHM;
      gtm->rFWHM = psfFWHM;
      gtm->sFWHM = psfFWHM;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--psf-col")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&gtm->cFWHM);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--psf-row")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&gtm->rFWHM);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--psf-slice")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&gtm->sFWHM);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--apply-fwhm")){
      // apply to input for testing
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&ApplyFWHM);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--rbv")) DoRBV = 1;
    else if(!strcasecmp(option, "--no-rescale")) gtm->rescale = 0;
    else if(!strcasecmp(option, "--rescale")) {
      if(nargc < 1) CMDargNErr(option,1);
      gtm->rescale = 1;
      gtm->n_scale_refids = 0;
      int nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) ){
	sscanf(pargv[nth],"%d",&gtm->scale_refids[gtm->n_scale_refids]);
	gtm->n_scale_refids++;
	nth ++;
      }
      nargsused = nth;
    }
    else if(!strcasecmp(option, "--mgpvc") || !strcasecmp(option, "--mg")){
      if(nargc < 1) CMDargNErr(option,1);
      gtm->DoMGPVC = 1;
      sscanf(pargv[0],"%lf",&gtm->mg_gmthresh);
      nargsused = 1;
      if(CMDnthIsArg(nargc, pargv, 1)){
	gtm->n_mg_refids = 0;
	int nth = 1;
	while(CMDnthIsArg(nargc, pargv, nth) ){
	  sscanf(pargv[nth],"%d",&gtm->mg_refids[gtm->n_mg_refids]);
	  gtm->n_mg_refids++;
	  nth ++;
	}
	nargsused += gtm->n_mg_refids;
      }
    } 
    else if(!strcasecmp(option, "--mg-ref-cerebral-wm")) {
      int m=0;
      gtm->mg_refids[m] =  2; m++;
      gtm->mg_refids[m] = 41; m++;
      gtm->n_mg_refids = m;
    }
    else if(!strcasecmp(option, "--mg-ref-lobes-wm")) {
      int m=0;
      gtm->mg_refids[m] = 3201; m++;
      gtm->mg_refids[m] = 3203; m++;
      gtm->mg_refids[m] = 3204; m++;
      gtm->mg_refids[m] = 3205; m++;
      gtm->mg_refids[m] = 3206; m++;
      gtm->mg_refids[m] = 3207; m++;
      gtm->mg_refids[m] = 4201; m++;
      gtm->mg_refids[m] = 4203; m++;
      gtm->mg_refids[m] = 4204; m++;
      gtm->mg_refids[m] = 4205; m++;
      gtm->mg_refids[m] = 4206; m++;
      gtm->mg_refids[m] = 4207; m++;
      gtm->mg_refids[m] = 5001; m++;
      gtm->mg_refids[m] = 5002; m++;
      gtm->n_mg_refids = m;
    }
    else if(!strcasecmp(option, "--km-ref")) {
      if(nargc < 1) CMDargNErr(option,1);
      gtm->DoKMRef = 1;
      gtm->n_km_refids = 0;
      int nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) ){
	sscanf(pargv[nth],"%d",&gtm->km_refids[gtm->n_km_refids]);
	gtm->n_km_refids++;
	nth ++;
      }
      nargsused = gtm->n_km_refids;
    } 
    else if(!strcasecmp(option, "--km-hb")) {
      if(nargc < 1) CMDargNErr(option,1);
      gtm->DoKMHB = 1;
      gtm->n_km_hbids = 0;
      int nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) ){
	sscanf(pargv[nth],"%d",&gtm->km_hbids[gtm->n_km_hbids]);
	gtm->n_km_hbids++;
	nth ++;
      }
      nargsused = gtm->n_km_hbids;
    } 
    else if(!strcasecmp(option, "--save-input")) SaveInput = 1;
    else if(!strcasecmp(option, "--X0"))   SaveX0 = 1;
    else if(!strcasecmp(option, "--X"))    SaveX = 1;
    else if(!strcasecmp(option, "--y"))    SaveYMat = 1;
    else if(!strcasecmp(option, "--beta")) SaveBetaMat = 1;

    else if(!strcasecmp(option, "--rvar-only")) RVarOnly = 1;
    else if(!strcasecmp(option, "--o")) {
      if(nargc < 1) CMDargNErr(option,1);
      OutDir = pargv[0];
      gtm->OutDir = OutDir;
      sprintf(tmpstr,"%s/rvar.dat",OutDir);
      RVarFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/skew.dat",OutDir);
      SkewFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/kurtosis.dat",OutDir);
      KurtosisFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/vrf.dat",OutDir);
      VRFStatsFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/beta.nii.gz",OutDir);
      OutBetaFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/beta.var.nii.gz",OutDir);
      OutBetaVarFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/seg.nii.gz",OutDir);
      OutSegFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/xtx.mat",OutDir);
      OutXtXFile = strcpyalloc(tmpstr);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--save-eres")) SaveEres=1;
    else if(!strcasecmp(option, "--save-yhat")) SaveYhat=1;
    else if(!strcasecmp(option, "--save-yhat0")) SaveYhat0=1;
    else if(!strcasecmp(option, "--save-yhat-full-fov")) SaveYhatFullFoV=1;
    else if(!strcasecmp(option, "--replace")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&gtm->SrcReplace[gtm->nReplace]);
      sscanf(pargv[1],"%d",&gtm->TrgReplace[gtm->nReplace]);
      gtm->nReplace++;
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--merge-hypos")) {
      gtm->SrcReplace[gtm->nReplace]=78; gtm->TrgReplace[gtm->nReplace]=77; gtm->nReplace++;
      gtm->SrcReplace[gtm->nReplace]=79; gtm->TrgReplace[gtm->nReplace]=77; gtm->nReplace++;
    } 
    else if(!strcasecmp(option, "--vg-thresh")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&vg_isEqual_Threshold);
      nargsused = 1;
    }
    else if (!strcmp(option, "--C")) {
      if(nargc < 1) CMDargNErr(option,1);
      int nc = gtm->nContrasts;
      gtm->contrasts[nc] = (GTMCON *) calloc(sizeof(GTMCON),1);
      gtm->contrasts[nc]->C = MatrixReadTxt(pargv[0], NULL);
      if(gtm->contrasts[nc]->C == NULL) exit(1);
      gtm->contrasts[nc]->name = fio_basename(pargv[0],".mtx"); //strip .mtx
      gtm->nContrasts++;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--synth")) {
      if(nargc < 2) CMDargNErr(option,2);
      SrcBetaFile = pargv[0];
      SynthFile   = pargv[1];
      mritmp = MRIread(SrcBetaFile);
      if(mritmp == NULL) exit(1);
      MATRIX *srcbetaT = fMRItoMatrix(mritmp,NULL);
      srcbeta = MatrixTranspose(srcbetaT,NULL);
      MatrixFree(&srcbetaT);
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--threads")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--max-threads")){
      nthreads = 1;
      #ifdef _OPENMP
      nthreads = omp_get_max_threads();
      omp_set_num_threads(nthreads);
      #endif
    } 
    else if(!strcasecmp(option, "--max-threads-1")){
      nthreads = 1;
      #ifdef _OPENMP
      nthreads = omp_get_max_threads()-1;
      if(nthreads < 0) nthreads = 1;
      omp_set_num_threads(nthreads);
      #endif
    } 
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(CMDsingleDash(option))
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
  printf("   --i   inputvol : source data to PVC\n");
  printf("   --psf psfmm : scanner PSF FWHM in mm\n");
  printf("   --seg segfile : anatomical segmentation to define regions for GTM\n");
  printf("   --reg reg.lta : LTA registration file that maps anatomical to PET\n");
  printf("   --o   outdir    : output directory\n");
  printf("\n");
  printf("   --mask volfile : ignore areas outside of the mask (in input vol space)\n");
  printf("   --auto-mask FWHM thresh : automatically compute mask\n");
  printf("   --no-reduce-fov : do not reduce FoV to encompass automask\n");
  printf("   --C contrast.mtx : univariate contrast to test (ascii text file)\n");
  printf("\n");
  printf("   --default-seg-merge : default schema for merging ROIs\n");
  printf("   --merge-hypos : merge left and right hypointensites into to ROI\n");
  printf("   --tt-reduce : reduce segmentation to that of a tissue type\n");
  printf("   --replace Id1 Id2 : replace seg Id1 with seg Id2\n");
  printf("   --replace-file : file with a list of Ids to replace\n");
  printf("   --reg-identity : assume that input is in anatomical space \n");
  printf("   --no-rescale   : do not global rescale such that mean of cerebellum WM is 100\n");
  printf("\n");
  printf("   --no-vox-frac-cor : do not use voxel fraction correction (with --psf 0 turns off PVC entirely)\n");
  printf("   --rbv            : perform RBV PVC\n");
  printf("   --mg gmthresh RefId1 RefId2 ...: perform Mueller-Gaertner PVC, gmthresh is min gm pvf bet 0 and 1\n");
  printf("   --mg-ref-cerebral-wm : set MG RefIds to 2 and 41\n");
  printf("   --mg-ref-lobes-wm : set MG RefIds to those for lobes when using wm subseg\n");
  printf("   --km-ref RefId1 RefId2 ... : compute reference TAC for KM as mean of given RefIds\n");
  printf("   --km-hb  RefId1 RefId2 ... : compute HiBinding TAC for KM as mean of given RefIds\n");
  printf("\n");
  printf("   --X : save X matrix in matlab4 format as X.mat (it will be big)\n");
  printf("   --y : save y matrix in matlab4 format as y.mat\n");
  printf("   --beta : save beta matrix in matlab4 format as beta.mat\n");
  printf("   --X0 : save X0 matrix in matlab4 format as X0.mat (it will be big)\n");
  printf("   --save-input : saves rescaled input as input.rescaled.nii.gz\n");
  printf("   --save-eres : saves residual error\n");
  printf("   --save-yhat : saves yhat\n");
  printf("   --save-yhat-full-fov : saves yhat in full FoV (if FoV was reduced)\n");
  printf("   --save-yhat0 : saves yhat prior to smoothing\n");
  printf("\n");
  printf("   --synth gtmbeta synthvolume : synthesize unsmoothed volume with gtmbeta as input\n");
  printf("\n");
  #ifdef _OPENMP
  printf("   --threads N : use N threads (with Open MP)\n");
  printf("   --max-threads : use the maximum allowable number of threads for this computer\n");
  printf("   --max-threads-1 : use one less than the maximum allowable number of threads for this computer\n");
  #endif
  printf("   --sd SUBJECTS_DIR\n");
  printf("   --vg-thresh thrshold : threshold for  'ERROR: LTAconcat(): LTAs 0 and 1 do not match'\n");
  printf("   --gdiag diagno : set diagnostic level\n");
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
  if(! fio_FileExistsReadable(SegVolFile)){
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,gtm->anat2pet->subject,SegVolFile);
    if(! fio_FileExistsReadable(tmpstr)){
      printf("ERROR: cannot find seg vol %s or %s\n",SegVolFile,tmpstr);
      exit(1);
    }
    SegVolFile = strcpyalloc(tmpstr);
  }

  if(gtm->cFWHM < 0 || gtm->rFWHM < 0 || gtm->sFWHM < 0){
    printf("ERROR: must spec psf FWHM\n");
    exit(1);
  }

  if(OutDir == NULL && SynthFile != NULL) SynthOnly=1;
  if(OutDir == NULL && SynthFile == NULL){
    printf("ERROR: must spec an output with --o or --synth\n");
    exit(1);
  }

  if(OutDir){
    if(SaveEres){
      sprintf(tmpstr,"%s/eres.nii.gz",OutDir);
      eresFile = strcpyalloc(tmpstr);
    }
    if(SaveYhat){
      sprintf(tmpstr,"%s/yhat.nii.gz",OutDir);
      yhatFile = strcpyalloc(tmpstr);
    }
    if(SaveYhatFullFoV){
      sprintf(tmpstr,"%s/yhat.fullfov.nii.gz",OutDir);
      yhatFullFoVFile = strcpyalloc(tmpstr);
    }
    if(SaveYhat0){
      sprintf(tmpstr,"%s/yhat0.nii.gz",OutDir);
      yhat0File = strcpyalloc(tmpstr);
    }
  }

  if(nthreads != 1){
    #ifdef _OPENMP
    printf("Setting maximum number of threads to %d\n",nthreads);
    omp_set_num_threads(nthreads);
    #endif
  }
  printf("Loading input %s\n",SrcVolFile);
  gtm->yvol = MRIread(SrcVolFile);
  if(gtm->yvol==NULL) exit(1);

  if(regidentity){
    printf("Using identity registration\n");
    gtm->anat2pet = TransformRegDat2LTA(gtm->yvol, gtm->yvol, NULL);
  }
  else {
    if(regfile == NULL){
      printf("ERROR: must spec regfile \n");
      exit(1);
    }
  }

  if(SaveX0) {
    sprintf(tmpstr,"%s/X0.mat",OutDir);
    X0file = strcpyalloc(tmpstr);
  }
  if(SaveX) {
    sprintf(tmpstr,"%s/X.mat",OutDir);
    Xfile = strcpyalloc(tmpstr);
  }
  if(SaveYMat) {
    sprintf(tmpstr,"%s/y.mat",OutDir);
    ymatfile = strcpyalloc(tmpstr);
  }
  if(SaveBetaMat) {
    sprintf(tmpstr,"%s/beta.mat",OutDir);
    betamatfile = strcpyalloc(tmpstr);
  }

  GTMcheckReplaceList(gtm->nReplace, gtm->SrcReplace, gtm->TrgReplace);

  if(gtm->DoMGPVC && gtm->n_mg_refids == 0){
    printf("ERROR: you must specify the ref segids for MG PVC either as extra \n"
	   "options to --mg or with --mg-ref-cerebral-wm\n");
    exit(1);
  }

  return;
}
/*---------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"setenv SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"cd %s\n",cwd);
  fprintf(fp,"%s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"vgthresh   %lf\n",vg_isEqual_Threshold);
  fprintf(fp,"nReplace   %d\n",gtm->nReplace);
  if(0 && gtm->nReplace > 0){
    fprintf(fp,"SegId replacement list\n");
    GTMprintReplaceList(fp, gtm->nReplace, gtm->SrcReplace, gtm->TrgReplace);
  }
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
/*
  \fn int VRFStats(GTM *gtm, double *vrfmean, double *vrfmin, double *vrfmax)
  \brief Computes the variance reduction factor (VRF) for each seg, the number
  of voxels for the seg, and betavar for each seg and frame. Also computes
  the max and min VRF.
*/
int VRFStats(GTM *gtm, double *vrfmean, double *vrfmin, double *vrfmax)
{
  int n,nvox,segid,f;
  double vrf;

  *vrfmean = 0;
  *vrfmax = 0;
  *vrfmin = 0;
  for(n=0; n < gtm->iXtX->rows; n++){
    vrf = (double) 1.0/gtm->iXtX->rptr[n+1][n+1];
    if(n==0){
      *vrfmax = vrf;
      *vrfmin = vrf;
    }
    if(*vrfmax < vrf) *vrfmax = vrf;
    if(*vrfmin > vrf) *vrfmin = vrf;
    *vrfmean += vrf;
  }
  *vrfmean /= gtm->iXtX->rows;

  if(gtm->betavar) MatrixFree(&gtm->betavar);
  gtm->betavar = MatrixAlloc(gtm->iXtX->rows,gtm->beta->cols,MATRIX_REAL);

  if(gtm->vrf) MatrixFree(&gtm->vrf);
  gtm->vrf = MatrixAlloc(gtm->iXtX->rows,1,MATRIX_REAL);

  if(gtm->nvox) MatrixFree(&gtm->nvox);
  gtm->nvox = MatrixAlloc(gtm->iXtX->rows,1,MATRIX_REAL);

  for(n=0; n < gtm->iXtX->rows; n++){
    segid = gtm->segidlist[n];
    nvox = MRIcountMatches(gtm->gtmseg, segid, 0, gtm->mask);
    vrf = (double) 1.0/gtm->iXtX->rptr[n+1][n+1];
    gtm->vrf->rptr[n+1][1]  = vrf;
    gtm->nvox->rptr[n+1][1] = nvox;
    for(f=0; f < gtm->beta->cols; f++)
      gtm->betavar->rptr[n+1][f+1] = gtm->rvar->rptr[1][f+1]/vrf;
  }

  return(0);
}
/*--------------------------------------------------------------------------*/
/*
  \fn int WriteVRFStats(char *fname, GTM *gtm)
  \brief Creates the vrf.dat file in the output folder. This is a text file
  that reports several statistics including the variance reduction factor (VRF).
*/
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
    vrf = gtm->vrf->rptr[n+1][1];
    nvox = gtm->nvox->rptr[n+1][1];
    cte = gtm->ctGTMSeg->entries[segid];
    fprintf(fp,"%3d %4d %-31s %-13s %6d %8.3f",n+1,segid,cte->name,
	    ttctab->entries[cte->TissueType]->name,nvox,vrf);
    //printf("%3d %4d %-31s %-13s %6d %8.3f",n+1,segid,cte->name,
    //ttctab->entries[cte->TissueType]->name,nvox,vrf);
    if(gtm->beta)    fprintf(fp,"   %10.3f",gtm->beta->rptr[n+1][1]);
    if(gtm->segrvar) fprintf(fp,"   %10.4f",sqrt(gtm->segrvar->rptr[n+1][1]));
    fprintf(fp,"\n");
  }
  fclose(fp);
  fflush(stdout);

  return(0);
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
			       gtmopt->USF, &gtmopt->anat2pvfseg);
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
  int nTT;
  LTA *ltaArray[2];
  struct timeb timer;
  TimerStart(&timer);

  //printf("GTMOPTcost() USF=%d\n",gtmopt->USF);
  nTT = gtmopt->gtm->ctGTMSeg->ctabTissueType->nentries-1;

  ltaArray[0] = gtmopt->pvfseg2anat;
  ltaArray[1] = gtmopt->anat2pet;
  gtmopt->pvfseg2pet = LTAconcat(ltaArray,2,1);
  //printf("Computing PVF\n");
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
  GTMbuildX(gtm);
  //printf("  buildX t = %g\n",TimerStop(&timer)/1000.0);fflush(stdout);
  if(gtm->X==NULL) return(-1);
  GTMsolve(gtm);

  gtmopt->tLastEval = TimerStop(&timer)/1000.0;
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn GTM *GTMalloc()
  \brief Allocates the GTM structure but nothing in the structure.
   sets PadThresh = .0001;
*/
GTM *GTMalloc()
{
  GTM *gtm;
  gtm = (GTM *) calloc(sizeof(GTM),1);
  gtm->PadThresh = .0001;
  return(gtm);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMfree(GTM **pGTM)
  \brief Frees a lot of stuff, but not everything.
*/
int GTMfree(GTM **pGTM)
{
  GTM *gtm = *pGTM;

  MRIfree(&gtm->yvol);
  //MRIfree(&gtm->gtmseg);
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
/*
  \fn int GTMmatrixY(GTM *gtm)
  \brief Converts the input gtm->yvol to a matrix using GTMvol2mat().
  It is important that this be done conistently with X, etc.
*/
int GTMmatrixY(GTM *gtm)
{
  gtm->y = GTMvol2mat(gtm, gtm->yvol, NULL);
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMsetNMask(GTM *gtm)
  \brief Computes the number of voxels in the mask. If the mask is
  NULL, then just computes the number of voxels in the input.
*/
int GTMsetNMask(GTM *gtm)
{
  if(gtm->mask) gtm->nmask = MRIcountAboveThreshold(gtm->mask, 0.5);
  else gtm->nmask = gtm->yvol->width*gtm->yvol->height*gtm->yvol->depth;
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMpsfStd(GTM *gtm)
  \brief Convert the PSF {crs}FWHM to a standard deviation.
*/
int GTMpsfStd(GTM *gtm)
{
  gtm->cStd = gtm->cFWHM/sqrt(log(256.0));
  gtm->rStd = gtm->rFWHM/sqrt(log(256.0));
  gtm->sStd = gtm->sFWHM/sqrt(log(256.0));
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMsegidlist(GTM *gtm)
  \brief Compute a sorted list of segmentation IDs from the segmentation
  itself (excludes 0). Just runs MRIsegIdListNot0().
*/
int GTMsegidlist(GTM *gtm)
{
  gtm->segidlist = MRIsegIdListNot0(gtm->anatseg, &gtm->nsegs, 0);
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMnPad(GTM *gtm)
  \brief Computes the number of voxels used to pad the tightest
  fitting bounding box around the nonzero voxels of the seg. There
  must be enough padding to account for spill-over created by
  smoothing the input by the PSF. It works by determining the distance
  from the center of a Gaussian such that the kernel equals PadThresh
  (a value between 0 and 1). The FWHM of the Guassian is the maximum
  of the {crs}FWHM. The way this should be interpreted is that any
  spill-over signal outside of the brain that is excluded by the
  bounding box will be no larger than PadThresh times the unsmoothed
  signal. Eg, if PadThresh is .001, then nPad will be computed such
  that the spill-over signal will be 0.1% less than the original
  signal. Using the bounding box can greatly reduce the size of 
  the input (esp for PET).
*/
int GTMnPad(GTM *gtm)
{
  double maxFWHM, maxStd;
  maxFWHM = MAX(gtm->cFWHM/gtm->yvol->xsize,
		MAX(gtm->rFWHM/gtm->yvol->ysize,gtm->sFWHM/gtm->yvol->zsize));
  if(maxFWHM > 0){
    // The case where psf=0, just correcting for volume fraction
    maxStd = maxFWHM*sqrt(log(256.0));
    gtm->nPad = ceil(sqrt(-log(gtm->PadThresh*maxStd*sqrt(2*M_PI))*2*maxStd));
    printf("maxFWHM = %g (voxels), PadThresh=%g, nPad=%d\n",maxFWHM,gtm->PadThresh,gtm->nPad);
  }
  else gtm->nPad = 1;
  return(0);
}
/*------------------------------------------------------------------*/
/*
  \fn int GTMsolve(GTM *gtm)
  \brief Solves the GTM using a GLM. X must already have been created.
  Computes Xt, XtX, iXtX, beta, yhat, res, dof, rvar, kurtosis, and skew.
  Also will rescale if rescaling. Returns 1 and computes condition
  number if matrix cannot be inverted. Otherwise returns 0.
*/
int GTMsolve(GTM *gtm)
{
  struct timeb timer;
  int n,f;
  double sum;

  if(gtm->X == NULL){
    printf("ERROR: GTMsolve(): must build design matrix first\n");
    exit(1);
  }

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
  if(gtm->rescale) GTMrescale(gtm);
  gtm->yhat = MatrixMultiplyD(gtm->X,gtm->beta,gtm->yhat);
  gtm->res  = MatrixSubtract(gtm->y,gtm->yhat,gtm->res);
  gtm->dof = gtm->X->rows - gtm->X->cols;
  if(gtm->rvar==NULL) gtm->rvar = MatrixAlloc(1,gtm->res->cols,MATRIX_REAL);
  for(f=0; f < gtm->res->cols; f++){
    sum = 0;
    for(n=0; n < gtm->res->rows; n++) sum += ((double)gtm->res->rptr[n+1][f+1]*gtm->res->rptr[n+1][f+1]);
    gtm->rvar->rptr[1][f+1] = sum/gtm->dof;
  }
  gtm->kurtosis = MatrixKurtosis(gtm->res,gtm->kurtosis);
  gtm->skew     = MatrixSkew(gtm->res,gtm->skew);
  return(0);
}
/*-----------------------------------------------------------------*/
/*
  \fn MRI *GTMmat2vol(GTM *gtm, MATRIX *m, MRI *vol)
  \brief Converts a matrix an MRI volume with rows=frames
  and columns=spatial dims. It is done row fastest, then col, then
  slice which makes it consistent with matlab. Any place that operates
  on the matrix data must be consistent when going from vol to matrix
  and back. See also GTMbuildX() and GTMvol2mat().
*/
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
/*
  \fn MATRIX *GTMvol2mat(GTM *gtm, MRI *vol, MATRIX *m)
  \brief Converts an MRI volume into a matrix with rows=frames
  and columns=spatial dims. It is done row fastest, then col, then
  slice which makes it consistent with matlab. Any place that operates
  on the matrix data must be consistent when going from vol to matrix
  and back. See also GTMbuildX() and GTMmat2vol(). 
*/
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
  for(s=0; s < vol->depth; s++){ // crs order is important here!
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
/*--------------------------------------------------------------------------*/
/*
  \fn int GTMsegrvar(GTM *gtm)
  \brief Computes the residual variance in each segmentation. Not perfect
  because the resdiual has spill-over. Hopefully it is meaningful for something.
*/
int GTMsegrvar(GTM *gtm)
{
  int c,r,s,f,k;
  int nthseg=0,segid;
  double v;

  gtm->segrvar = MatrixAlloc(gtm->nsegs,gtm->beta->cols,MATRIX_REAL);
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
	  if(segid != 0) gtm->segrvar->rptr[nthseg+1][f+1] += v*v;
	}
	k++;
      }// r 
    }// c
  }// s

  for(f=0; f < gtm->beta->cols; f++) {
    for(nthseg=0; nthseg < gtm->nsegs; nthseg++){
      v = gtm->segrvar->rptr[nthseg+1][f+1];
      gtm->segrvar->rptr[nthseg+1][f+1] = v/gtm->nperseg[nthseg];
    }
  }
  return(0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMrbv(GTM *gtm)
  \brief Performs Region-based Voxel-wise PVF correction. Can reduce the FoV
  of the output to a bounding box that tightly fits around the brain to
  reduce memory and disk space. The output shares a RAS space with the 
  anatseg (which shares a RAS space with the conformed anat) so a new
  registration is not needed.
 */
int GTMrbv(GTM *gtm)
{
  int c,r,s,f,nthseg,segid;
  double val,v,vhat0,vhat,v2;
  LTA *lta;
  struct timeb mytimer;
  MATRIX *nhits;

  if(gtm->rbv)      MRIfree(&gtm->rbv);
  if(gtm->yhat0seg) MRIfree(&gtm->yhat0seg);
  if(gtm->yhatseg)  MRIfree(&gtm->yhatseg);
  if(gtm->yseg)     MRIfree(&gtm->yseg);

  TimerStart(&mytimer) ; 

  printf("   Synthesizing unsmoothed input in seg space... ");fflush(stdout);
  gtm->yhat0seg = GTMsegSynth(gtm);
  if(gtm->yhat0seg == NULL){
    printf("ERROR: GTMrbv() could not synthesize yhat0seg\n");
    return(1);
  }
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);


  printf("   Smoothing synthesized in seg space... ");fflush(stdout);
  gtm->yhatseg = MRIgaussianSmoothNI(gtm->yhat0seg, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  if(gtm->yhatseg == NULL){
    printf("ERROR: GTMrbv() could not smooth yhatseg\n");
    return(1);
  }
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);

  printf("   Sampling input to seg space with trilin... ");fflush(stdout);
  gtm->yseg = MRIallocSequence(gtm->anatseg->width, gtm->anatseg->height, gtm->anatseg->depth,
			      MRI_FLOAT, gtm->yvol->nframes);
  if(gtm->yseg == NULL){
    printf("ERROR: GTMrbv() could not alloc yseg\n");
    return(1);
  }
  MRIcopyHeader(gtm->anatseg,gtm->yseg);
  MRIcopyPulseParameters(gtm->yvol,gtm->yseg);
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);

  lta = LTAcopy(gtm->seg2pet,NULL);
  LTAchangeType(lta,LINEAR_VOX_TO_VOX);
  MRIvol2Vol(gtm->yvol,gtm->yseg,(lta->xforms[0].m_L),SAMPLE_TRILINEAR, 0.0);
  LTAfree(&lta);

  printf("   Computing RBV ... ");fflush(stdout);
  gtm->rbv = MRIallocSequence(gtm->anatseg->width, gtm->anatseg->height, gtm->anatseg->depth,
			      MRI_FLOAT, gtm->yvol->nframes);
  if(gtm->rbv == NULL){
    printf("ERROR: GTMrbv() could not alloc rbv\n");
    return(1);
  }
  MRIcopyHeader(gtm->anatseg,gtm->rbv);
  MRIcopyPulseParameters(gtm->yvol,gtm->rbv);

  gtm->rbvsegmean = MRIallocSequence(gtm->nsegs,1,1,MRI_FLOAT,gtm->yvol->nframes);
  nhits = MatrixAlloc(gtm->beta->rows,1,MATRIX_REAL);
  for(s=0; s < gtm->anatseg->depth; s++){ // crs order not important
    for(c=0; c < gtm->anatseg->width; c++){
      for(r=0; r < gtm->anatseg->height; r++){
	segid = MRIgetVoxVal(gtm->anatseg,c,r,s,0);
	if(segid < 0.5) continue;
	for(nthseg=0; nthseg < gtm->nsegs; nthseg++) if(gtm->segidlist[nthseg] == segid) break;
	nhits->rptr[nthseg+1][1] ++;
	for(f=0; f < gtm->yvol->nframes; f++){
	  v     = MRIgetVoxVal(gtm->yseg,c,r,s,f);
	  vhat0 = MRIgetVoxVal(gtm->yhat0seg,c,r,s,f);
	  vhat  = MRIgetVoxVal(gtm->yhatseg,c,r,s,f);
	  val = v*vhat0/(vhat+FLT_EPSILON); // RBV equation
	  MRIsetVoxVal(gtm->rbv,c,r,s,f,val);
	  v2 = MRIgetVoxVal(gtm->rbvsegmean,nthseg,0,0,f); // track seg means for QA
	  MRIsetVoxVal(gtm->rbvsegmean,nthseg,0,0,f,v2+val);
	}
      }
    }
  }
  MRIfree(&gtm->yseg);
  MRIfree(&gtm->yhat0seg);
  MRIfree(&gtm->yhatseg);
  printf("  t = %4.2f min\n",TimerStop(&mytimer)/60000.0);fflush(stdout);

  // track seg means for QA
  for(nthseg=0; nthseg < gtm->nsegs; nthseg++){
    for(f=0; f < gtm->yvol->nframes; f++){
      val = MRIgetVoxVal(gtm->rbvsegmean,nthseg,0,0,f)/nhits->rptr[nthseg+1][1];
      MRIsetVoxVal(gtm->rbvsegmean,nthseg,0,0,f,val);
    }
  }
  MatrixFree(&nhits);

  if(gtm->mask_rbv_to_brain){
    // Reduce RBV to a tight mask around the brain. This can greatly reduce 
    // memory requirements. The RAS space is still that of the anatseg
    // (and so also that of the conformed anat) so no new registration is necessary
    printf("   masking RBV to brain\n");
    int n, nReplace, ReplaceThis[1000], WithThat[1000];
    MRI *segtmp,*rbvtmp;
    MRI_REGION *region;
    nReplace = 0;
    for(n=0; n < gtm->ctGTMSeg->nentries; n++){
      if(gtm->ctGTMSeg->entries[n] == NULL)  continue;
      if(gtm->ctGTMSeg->entries[n]->TissueType != 5) continue; // should not hard-code
      ReplaceThis[nReplace] = n;
      WithThat[nReplace] = 0;
      nReplace++;
    }
    printf("  replacing head voxels with 0\n");
    segtmp = MRIreplaceList(gtm->anatseg, ReplaceThis, WithThat, nReplace, NULL);
    printf("  computing bounding box  ");
    region = REGIONgetBoundingBox(segtmp,10);
    REGIONprint(stdout, region);
    printf("  extracting bounding box\n");
    rbvtmp = MRIextractRegion(gtm->rbv, NULL, region);
    free(region); region=NULL;
    MRIfree(&segtmp);
    MRIfree(&gtm->rbv);
    gtm->rbv = rbvtmp;
  }
  printf("  RBV took %4.2f min\n",TimerStop(&mytimer)/60000.0);

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
/*
  \fn int GTMmgpvc(GTM *gtm)
  \brief Performs Muller-Gartner PVC. Hardcodes tissue type IDs to 
  be 0=cortex, 1=subcortexgm, 2=WM.
 */
int GTMmgpvc(GTM *gtm)
{
  int c,r,s,f,n,nthseg,segid,found,nhits;
  double sum,vgmpsf,vwmpsf,vwmtac,vtac,vmgtac;
  MRI *ctxpvf, *subctxpvf, *wmpvf, *gmpvf,*gmpvfpsf,*wmpvfpsf;

  // Compute the MG reference TAC
  gtm->mg_reftac = MatrixAlloc(gtm->yvol->nframes,1,MATRIX_REAL);
  for(f=0; f < gtm->yvol->nframes; f++){
    nhits = 0;
    sum = 0;
    for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
      segid = gtm->segidlist[nthseg];
      found = 0;
      for(n=0; n < gtm->n_mg_refids; n++) {
	if(segid == gtm->mg_refids[n]) {
	  found = 1;
	  nhits++;
	  break;
	}
      }
      if(!found) continue;
      sum += gtm->beta->rptr[nthseg+1][f+1];
      printf("   n=%d, nthseg=%d %g\n",n,nthseg,gtm->beta->rptr[nthseg+1][f+1]);
    }
    gtm->mg_reftac->rptr[f+1][1] = sum/nhits;
    printf("   wm tac %2d %2d %g\n",f,nhits,gtm->mg_reftac->rptr[f+1][1]);
  }

  if(gtm->mg) MRIfree(&gtm->mg);
  gtm->mg = MRIallocSequence(gtm->yvol->width, gtm->yvol->height, gtm->yvol->depth,
			   MRI_FLOAT, gtm->yvol->nframes);
  if(gtm->mg == NULL) return(1);
  MRIcopyHeader(gtm->yvol,gtm->mg);

  // Compute gray matter PVF with smoothing
  ctxpvf    = fMRIframe(gtm->ttpvf,0,NULL); // cortex PVF
  subctxpvf = fMRIframe(gtm->ttpvf,1,NULL); // subcortex GM PVF
  gmpvf = MRIadd(ctxpvf,subctxpvf,NULL); // All GM PVF
  // Smooth GM PVF by PSF
  gmpvfpsf = MRIgaussianSmoothNI(gmpvf,gtm->cStd, gtm->rStd, gtm->sStd, NULL);

  // WM PVF
  wmpvf = fMRIframe(gtm->ttpvf,2,NULL); 
  // Smooth GM PVF by PSF
  wmpvfpsf = MRIgaussianSmoothNI(wmpvf,gtm->cStd, gtm->rStd, gtm->sStd, NULL);

  // Finally, do the actual MG correction
  for(c=0; c < gtm->yvol->width; c++){ // crs order not important
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
  MRIfree(&ctxpvf);
  MRIfree(&subctxpvf);
  MRIfree(&wmpvf);
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
/*
  \fn int GTMsynth(GTM *gtm)
  \brief Synthesizes the unsmoothed PET image by computing 
   ysynth = X0*beta and then re-packing the result into a volume.
   This is mainly useful for debugging.
 */
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
/*------------------------------------------------------------------*/
/*
  \fn int GTMcheckX(MATRIX *X)
  \brief Checks that all colums sum to 1
 */
int GTMcheckX(MATRIX *X)
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

  printf("GTMcheckX: count=%d, dmax=%g\n",count,dmax);
  return(count);
}
/*------------------------------------------------------------------------------*/
/*
  \fn int GTMbuildX(GTM *gtm)
  \brief Builds the GTM design matrix both with (X) and without (X0) PSF.  If 
  gtm->DoVoxFracCor=1 then corrects for volume fraction effect.
*/
int GTMbuildX(GTM *gtm)
{
  int nthseg;
  struct timeb timer;

  if(gtm->X==NULL || gtm->X->rows != gtm->nmask || gtm->X->cols != gtm->nsegs){
    // Alloc or realloc X
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

  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
  for(nthseg = 0; nthseg < gtm->nsegs; nthseg++){
    int segid,k,c,r,s;
    MRI *nthsegpvf=NULL,*nthsegpvfbb=NULL,*nthsegpvfbbsm=NULL;
    MRI_REGION *region;
    segid = gtm->segidlist[nthseg];
    //printf("nthseg = %d, %d %6.4f\n",nthseg,segid,TimerStop(&timer)/1000.0);fflush(stdout);
    if(gtm->DoVoxFracCor)
      nthsegpvf = fMRIframe(gtm->segpvf,nthseg,NULL); // extract PVF for this seg
    else
      nthsegpvf = MRIbinarizeMatch(gtm->gtmseg,&segid,1,0,NULL); // or get binary mask
    // Extract a region for the seg. This speeds up smoothing.
    region      = REGIONgetBoundingBox(nthsegpvf,gtm->nPad); // tight+pad bounding box 
    nthsegpvfbb = MRIextractRegion(nthsegpvf, NULL, region) ; // extract BB
    nthsegpvfbbsm = MRIgaussianSmoothNI(nthsegpvfbb, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
    // Fill X, creating X in this order makes it consistent with matlab
    // Note: y must be ordered in the same way. See GTMvol2mat()
    k = 0;
    for(s=0; s < gtm->yvol->depth; s++){
      for(c=0; c < gtm->yvol->width; c++){
	for(r=0; r < gtm->yvol->height; r++){
	  if(gtm->mask && MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	  k ++; // have to incr here in case continue below
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

/*--------------------------------------------------------------------------*/
/*
  \fn MRI *GTMsegSynth(GTM *gtm)
  \brief Creates a volume that is the same size as the anatseg in which the value
  at a voxel equals the beta of the seg ID at that voxel. This is used for RBV.
*/
MRI *GTMsegSynth(GTM *gtm)
{
  int c,r,s,f,segid,segno;
  MRI *synth;

  synth = MRIallocSequence(gtm->anatseg->width,gtm->anatseg->height,gtm->anatseg->depth,MRI_FLOAT,gtm->beta->cols);
  if(synth == NULL){
    printf("ERROR: GTMsegSynth(): could not alloc\n");
    return(NULL);
  }
  MRIcopyHeader(gtm->anatseg,synth);
  MRIcopyPulseParameters(gtm->yvol,synth);

  for(c=0; c < gtm->anatseg->width; c++){ // crs order does not matter here
    for(r=0; r < gtm->anatseg->height; r++){
      for(s=0; s < gtm->anatseg->depth; s++){
	segid = MRIgetVoxVal(gtm->anatseg,c,r,s,0);
	if(segid == 0) continue;
	for(segno=0; segno < gtm->nsegs; segno++)
	  if(segid == gtm->segidlist[segno]) break;
	if(segno == gtm->nsegs){
	  printf("ERROR: GTMsegSynth(): could not find a match for segid=%d\n",segid);
	  for(segno=0; segno < gtm->nsegs; segno++) printf("%3d %5d\n",segno,gtm->segidlist[segno]);
	  return(NULL);
	}
	for(f=0; f < gtm->beta->cols; f++)
	  MRIsetVoxVal(synth,c,r,s,f,gtm->beta->rptr[segno+1][f+1]);
      }
    }
  }

  return(synth);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMprintMGRefTAC(GTM *gtm, FILE *fp)
  \brief Prints the Muller-Gartner WM reference TAC to the given file pointer
*/
int GTMprintMGRefTAC(GTM *gtm, FILE *fp)
{
  int f;
  for(f=0; f < gtm->yvol->nframes; f++)
    fprintf(fp,"%3d %10.5f\n",f,gtm->mg_reftac->rptr[f+1][1]);

  return(0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMwriteMGRefTAC(GTM *gtm, char *filename)
  \brief Writes the Muller-Gartner WM reference TAC to the given file
*/
int GTMwriteMGRefTAC(GTM *gtm, char *filename)
{
  FILE *fp;
  fp = fopen(filename,"w");
  GTMprintMGRefTAC(gtm, fp);
  fclose(fp);
  return(0);
}

/*--------------------------------------------------------------------------*/
/*
  \fn int GTMrescale(GTM *gtm)
  \brief Computes global rescaling factor and applies it to yvol, beta, and y.
  The factor = 100/mean(beta_i) where i is the list of scale seg IDs (scale_ref_ids)
*/
int GTMrescale(GTM *gtm)
{
  int f,n,nthseg,segid,found,nhits;
  double sum;

  // global rescaling
  nhits = 0;
  sum = 0;
  for(f=0; f < gtm->yvol->nframes; f++){
    for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
      segid = gtm->segidlist[nthseg];
      found = 0;
      for(n=0; n < gtm->n_scale_refids; n++) {
	if(segid == gtm->scale_refids[n]) {
	  found = 1;
	  nhits++;
	  break;
	}
      }
      if(!found) continue;
      sum += gtm->beta->rptr[nthseg+1][f+1];
    }
  }
  gtm->scale = 100/(sum/nhits);
  printf("gtm multiplicative scale %g\n",gtm->scale);

  MRImultiplyConst(gtm->yvol,gtm->scale,gtm->yvol);
  MatrixScalarMul(gtm->beta, gtm->scale,gtm->beta);
  MatrixScalarMul(gtm->y,    gtm->scale,gtm->y);

  return(0);
}

/*-------------------------------------------------------------------------------*/
/*
  \fn int GTMttest(GTM *gtm)
  \brief Computes a t-test for each contrast. This includes computing gamma, 
  gammavar, t, and p.
*/
int GTMttest(GTM *gtm)
{
  MATRIX *Ct,*CiXtX,*CiXtXCt;
  int n,f,nframes;
  GTMCON *gtmc;
  double F;

  nframes = gtm->beta->cols;

  for(n=0; n < gtm->nContrasts; n++){
    gtmc = gtm->contrasts[n];
    gtmc->gamma = MatrixMultiply(gtmc->C,gtm->beta,NULL);
    gtmc->gammavar = MatrixAlloc(1,nframes,MATRIX_REAL);
    gtmc->t = MatrixAlloc(1,nframes,MATRIX_REAL);
    gtmc->p = MatrixAlloc(1,nframes,MATRIX_REAL);
    Ct = MatrixTranspose(gtmc->C,NULL);
    CiXtX = MatrixMultiply(gtmc->C,gtm->iXtX,NULL);
    CiXtXCt = MatrixMultiply(CiXtX,Ct,NULL);
    for(f=0; f < nframes; f++){
      gtmc->gammavar->rptr[1][f+1] = CiXtXCt->rptr[1][1]*gtm->rvar->rptr[1][f+1];
      gtmc->t->rptr[1][f+1] = gtmc->gamma->rptr[1][f+1]/sqrt(gtmc->gammavar->rptr[1][f+1]);
      F = gtmc->t->rptr[1][f+1]*gtmc->t->rptr[1][f+1];
      gtmc->p->rptr[1][f+1] = sc_cdf_fdist_Q(F,gtmc->C->rows,gtm->dof);
    }
    MatrixFree(&Ct);
    MatrixFree(&CiXtX);
    MatrixFree(&CiXtXCt);
  }
  return(0);
}
/*------------------------------------------------------------------------------*/
/*
  \fn int GTMwriteContrasts(GTM *GTM)
  \brief Creates an ASCII file in the output folder for each contrast with gamma, 
  gammavar, t, and p
 */
int GTMwriteContrasts(GTM *GTM)
{
  int n,nframes,f;
  GTMCON *gtmc;
  char tmpstr[5000];
  FILE *fp;

  nframes = gtm->beta->cols;

  for(n=0; n < gtm->nContrasts; n++){
    gtmc = gtm->contrasts[n];
    sprintf(tmpstr,"%s/%s.mat",gtm->OutDir,gtmc->name);
    MatlabWrite(gtmc->C, tmpstr,"C");
    sprintf(tmpstr,"%s/%s.dat",gtm->OutDir,gtmc->name);
    fp = fopen(tmpstr,"w");
    for(f=0; f < nframes; f++){
      fprintf(fp,"%2d %20.15f %20.15f %20.15f %20.15f\n",f,
	      gtmc->gamma->rptr[1][1],gtmc->gammavar->rptr[1][1],
	      gtmc->t->rptr[1][1],gtmc->p->rptr[1][1]);
    }
    fclose(fp);
  }
  return(0);
}

/*-------------------------------------------------------------------------*/
/*
  \fn int GTMcheckRefIds(GTM *gtm)
  \brief Checks the segmentation IDs used for rescaling, MG, KM Ref, and KM HB
  to make sure that they are in the segmentation.
 */
int GTMcheckRefIds(GTM *gtm)
{
  int n,m,ok;

  if(gtm->rescale){
    for(n=0; n < gtm->n_scale_refids; n++){
      ok = 0;
      for(m=0; m < gtm->nsegs; m++){
	if(gtm->segidlist[m] == gtm->scale_refids[n]){
	  ok = 1;
	  break;
	}
      }
      if(! ok) {
	printf("ERROR: scale reference id %d cannot be found in seg id list\n",gtm->scale_refids[n]);
	fprintf(gtm->logfp,"ERROR: scale reference id %d cannot be found in seg id list\n",gtm->scale_refids[n]);
	return(1);
      }
    }
  }

  if(gtm->DoMGPVC){
    for(n=0; n < gtm->n_mg_refids; n++){
      ok = 0;
      for(m=0; m < gtm->nsegs; m++){
	if(gtm->segidlist[m] == gtm->mg_refids[n]){
	  ok = 1;
	  break;
	}
      }
      if(! ok) {
	printf("ERROR: MG reference id %d cannot be found in seg id list\n",gtm->mg_refids[n]);
	fprintf(gtm->logfp,"ERROR: scale reference id %d cannot be found in seg id list\n",gtm->mg_refids[n]);
	return(1);
      }
    }
  }

  if(gtm->DoKMRef){
    for(n=0; n < gtm->n_km_refids; n++){
      ok = 0;
      for(m=0; m < gtm->nsegs; m++){
	if(gtm->segidlist[m] == gtm->km_refids[n]){
	  ok = 1;
	  break;
	}
      }
      if(! ok) {
	printf("ERROR: KM reference id %d cannot be found in seg id list\n",gtm->km_refids[n]);
	fprintf(gtm->logfp,"ERROR: scale reference id %d cannot be found in seg id list\n",gtm->km_refids[n]);
	return(1);
      }
    }
  }

  if(gtm->DoKMHB){
    for(n=0; n < gtm->n_km_hbids; n++){
      ok = 0;
      for(m=0; m < gtm->nsegs; m++){
	if(gtm->segidlist[m] == gtm->km_hbids[n]){
	  ok = 1;
	  break;
	}
      }
      if(! ok) {
	printf("ERROR: KM high binding id %d cannot be found in seg id list\n",gtm->km_hbids[n]);
	fprintf(gtm->logfp,"ERROR: scale high binding id %d cannot be found in seg id list\n",gtm->km_hbids[n]);
	return(1);
      }
    }
  }

  return(0);
}

/*
  \fn int GTMprintRefIds(GTM *gtm, FILE *fp)
  \brief Prints the segmentation IDs used for rescaling, MG, KM Ref, and KM HB
  to the given FILE pointer.
 */
int GTMprintRefIds(GTM *gtm, FILE *fp)
{
  int n;

  if(gtm->rescale){
    fprintf(fp,"Segmentations used for rescaling\n");
    for(n=0; n < gtm->n_scale_refids; n++) 
      fprintf(fp,"%4d %s\n",gtm->scale_refids[n],gtm->ctGTMSeg->entries[gtm->scale_refids[n]]->name);
  }

  if(gtm->DoMGPVC){
    fprintf(fp,"Segmentations used for MG PVC WM reference\n");
    for(n=0; n < gtm->n_mg_refids; n++) 
      fprintf(fp,"%4d %s\n",gtm->mg_refids[n],gtm->ctGTMSeg->entries[gtm->mg_refids[n]]->name);
  }

  if(gtm->DoKMRef){
    fprintf(fp,"Segmentations used for KM reference TAC\n");
    for(n=0; n < gtm->n_km_refids; n++)
      fprintf(fp,"%4d %s\n",gtm->km_refids[n],gtm->ctGTMSeg->entries[gtm->km_refids[n]]->name);
  }

  if(gtm->DoKMHB){
    fprintf(fp,"Segmentations used for KM High Binding TAC\n");
    for(n=0; n < gtm->n_km_hbids; n++)
      fprintf(fp,"%4d %s\n",gtm->km_hbids[n],gtm->ctGTMSeg->entries[gtm->km_hbids[n]]->name);
  }
  fflush(fp);
  return(0);
}
/*
  \fn int GTMrefTAC(GTM *gtm)
  \brief Computes KM reference and hibinding TACs. Also writes them
  to the output folder. The HB TAC is written as both an ascii file
  and a nii.gz (the later needed for KM analysis with mri_glmfit)
 */
int GTMrefTAC(GTM *gtm)
{
  int f,n,nthseg,segid;
  double sum;
  char tmpstr[5000];
  FILE *fp;

  if(gtm->DoKMRef){
    sprintf(tmpstr,"%s/km.ref.tac.dat",gtm->OutDir);
    fp = fopen(tmpstr,"w");
    if(fp==NULL) return(1);
    gtm->km_reftac = MatrixAlloc(gtm->yvol->nframes,1,MATRIX_REAL);
    for(f=0; f < gtm->yvol->nframes; f++){
      sum = 0;
      for(n=0; n < gtm->n_km_refids; n++) {
	for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
	  segid = gtm->segidlist[nthseg];
	  if(segid == gtm->km_refids[n]) break;
	}
	sum += gtm->beta->rptr[nthseg+1][f+1];
      }
      gtm->km_reftac->rptr[f+1][1] = sum/gtm->n_km_refids;
      fprintf(fp,"%20.15lf\n",sum/gtm->n_km_refids);
    }
    fclose(fp);
  }

  if(gtm->DoKMHB){
    sprintf(tmpstr,"%s/km.hb.tac.dat",gtm->OutDir);
    fp = fopen(tmpstr,"w");
    if(fp==NULL) return(1);
    gtm->km_hbtac = MatrixAlloc(gtm->yvol->nframes,1,MATRIX_REAL);
    for(f=0; f < gtm->yvol->nframes; f++){
      sum = 0;
      for(n=0; n < gtm->n_km_hbids; n++) {
	for(nthseg = 0; nthseg < gtm->nsegs; nthseg++) {
	  segid = gtm->segidlist[nthseg];
	  if(segid == gtm->km_hbids[n]) break;
	}
	sum += gtm->beta->rptr[nthseg+1][f+1];
      }
      gtm->km_hbtac->rptr[f+1][1] = sum/gtm->n_km_hbids;
      fprintf(fp,"%20.15lf\n",gtm->km_hbtac->rptr[f+1][1]);
    }
    fclose(fp);

    MRI *mritmp = MRIallocSequence(1, 1, 1, MRI_FLOAT, gtm->yvol->nframes);
    for(f=0; f < gtm->yvol->nframes; f++)
      MRIsetVoxVal(mritmp,0,0,0,f, gtm->km_hbtac->rptr[f+1][1]);
    sprintf(tmpstr,"%s/km.hb.tac.nii.gz",gtm->OutDir);
    MRIwrite(mritmp,tmpstr);
    MRIfree(&mritmp);
  }

  return(0);
}

/*
  \fn int GTMprintReplaceList(FILE *fp, const int nReplace, const int *ReplaceThis, const int *WithThat)
  \brief Prints the replacement list to the FILE pointer.
 */
int GTMprintReplaceList(FILE *fp, const int nReplace, const int *ReplaceThis, const int *WithThat)
{
  int n;
  for(n=0; n < nReplace; n++)
    fprintf(fp,"%5d %5d\n",ReplaceThis[n],WithThat[n]);
  return(0);
}

/*
  \fn int GTMcheckReplaceList(const int nReplace, const int *ReplaceThis, const int *WithThat)
  \brief Checks replacement list to make sure that no item in ReplaceThis list appears in
  the WithThat list.
 */
int GTMcheckReplaceList(const int nReplace, const int *ReplaceThis, const int *WithThat)
{
  int n,m;
  for(n=0; n < nReplace; n++){
    for(m=0; m < nReplace; m++){
      if(ReplaceThis[n] == WithThat[m]){
	printf("ERROR: item %d appears as both source and target seg id in replacement list\n",ReplaceThis[n]);
	return(1);
      }
    }
  }
  return(0);
}

/*
  \fn int GTMloadReplacmentList(const char *fname, int *nReplace, int *ReplaceThis, int *WithThat)
  \brief Loads in data from a file. The file should have two columns of numbers. The first
  is the segmentation ID to replace the second is the segmentation ID to replace it with.
  It will ignore any line that begins with a #. 
 */
int GTMloadReplacmentList(const char *fname, int *nReplace, int *ReplaceThis, int *WithThat)
{
  FILE *fp;
  int nlist,r,nth;
  char tmpstr[1001], *s;

  fp = fopen(fname,"r");
  if(fp==NULL){
    int e=errno;
    printf("ERROR: could not open %s %d\n",fname,e);
    printf("%s\n",strerror(e));
    return(1);
  }
  nlist = *nReplace;

  nth = -1;
  while(1){
    nth++;
    s = fgets(tmpstr,1000,fp);
    if(s==NULL) break;
    if(tmpstr[0] == '#') continue;
    r = CountItemsInString(tmpstr);
    if(r != 2){
      printf("ERROR: line %d in %s has the wrong format (has %d elements, expecting 2)\n",nth,fname,r);
      return(1);
    }
    sscanf(tmpstr,"%d %d",&ReplaceThis[nlist],&WithThat[nlist]); 
    nlist++;
  }
  if(nlist == 0){
    printf("ERROR: could not find any replacement items in %s\n",fname);
    return(1);
  }
  printf("Read in %d replacement items from %s\n",nlist,fname);
  *nReplace += nlist;
  return(0);
}
/*------------------------------------------------------------*/
/*
  \fn int GTMautoMask(GTM *gtm)
  \brief Computes a mask in PET space based on the segmentation and
  smoothing level. Optionally, it reduces the FoV of PET to be tight
  with the mask. If this is done, it recomputes the registration
  matrices.  The header of the full FoV PET is kept in
  gtm->yvol_full_fov but the pixel data are freed to reduce storage.
*/
int GTMautoMask(GTM *gtm)
{
  LTA *lta,*pet2bbpet,*seg2bbpet,*anat2bbpet;
  double std;
  MRI *masktmp,*yvoltmp;

  gtm->mask = MRIalloc(gtm->yvol->width,gtm->yvol->height,gtm->yvol->depth,MRI_FLOAT);
  MRIcopyHeader(gtm->yvol,gtm->mask);
  if(LTAmriIsSource(gtm->seg2pet,gtm->yvol)) lta = LTAcopy(gtm->seg2pet,NULL);
  else                                       lta = LTAinvert(gtm->seg2pet,NULL);        
  LTAchangeType(lta,LINEAR_VOX_TO_VOX);

  // Sample the seg into the pet space
  MRIvol2Vol(gtm->anatseg,gtm->mask,lta->xforms[0].m_L,SAMPLE_NEAREST,0);
  LTAfree(&lta);
  // Threshold it at 0.5 to give the mask
  MRIbinarize(gtm->mask,gtm->mask,0.5,0,1);
  // Smooth binary mask
  std = gtm->automask_fwhm/sqrt(log(256.0)); // convert fwhm to std
  MRIgaussianSmoothNI(gtm->mask, std,std,std, gtm->mask);
  // Binarize again to get final mask
  MRIbinarize(gtm->mask,gtm->mask,gtm->automask_thresh,0,1);

  if(gtm->automask_reduce_fov){
    printf("Automask, reducing FOV\n");
    gtm->automaskRegion = REGIONgetBoundingBox(gtm->mask,1);
    printf("region %d %d %d reduced to ",gtm->yvol->width,gtm->yvol->height,gtm->yvol->depth);
    REGIONprint(stdout, gtm->automaskRegion);
    fflush(stdout);
    masktmp = MRIextractRegion(gtm->mask, NULL, gtm->automaskRegion);
    yvoltmp = MRIextractRegion(gtm->yvol, NULL, gtm->automaskRegion);
    pet2bbpet = TransformRegDat2LTA(gtm->yvol, yvoltmp, NULL);
    seg2bbpet = LTAconcat2(gtm->seg2pet,pet2bbpet, 1);
    if(LTAmriIsSource(gtm->anat2pet,gtm->yvol)) lta = LTAinvert(gtm->anat2pet,NULL);        
    else                                        lta = LTAcopy(gtm->anat2pet,NULL);
    anat2bbpet = LTAconcat2(lta,pet2bbpet,1);
    LTAfree(&lta);
    gtm->yvol_full_fov = MRIallocHeader(gtm->yvol->width,gtm->yvol->height,gtm->yvol->depth,MRI_FLOAT,1);
    MRIcopyHeader(gtm->yvol,gtm->yvol_full_fov);
    MRIcopyPulseParameters(gtm->yvol,gtm->yvol_full_fov);
    MRIfree(&gtm->yvol);
    gtm->yvol = yvoltmp;
    MRIfree(&gtm->mask);
    gtm->mask = masktmp;
    LTAfree(&gtm->seg2pet);
    gtm->seg2pet = seg2bbpet;
    LTAfree(&gtm->anat2pet);
    gtm->anat2pet = anat2bbpet;
  }

  return(0);
}
/*
  \fn int GTMrvarGM(GTM *gtm)
  \brief Computes residual variance only in GM structures. The measure
  is not perfect because there will be spill-over from adjacent
  structures. But it will not include much from outside of the head (as
  the standard rvar measure will). The result is stored in
  gtm->rvargm. This hard-codes GM as 1 and 2 in the tissue type CTAB.
  Also writes rvar.gm.dat in the output folder. The value is computed
  as the sum of res.^2 in GM divided by the number of GM voxels.x  
*/
int GTMrvarGM(GTM *gtm)
{
  COLOR_TABLE *ct;
  int f,s,c,r,n,nhits,segid,tt;
  double sum;
  FILE *fp;
  char tmpstr[2000];

  if(gtm->rvargm==NULL) gtm->rvargm = MatrixAlloc(1,gtm->res->cols,MATRIX_REAL);
  ct = gtm->ctGTMSeg;

  for(f=0; f < gtm->res->cols; f++){
    sum = 0;
    n = -1;
    nhits = 0;
    // slice, col, row order is important here as is skipping the mask
    for(s=0; s < gtm->yvol->depth; s++){
      for(c=0; c < gtm->yvol->width; c++){
	for(r=0; r < gtm->yvol->height; r++){
	  if(MRIgetVoxVal(gtm->mask,c,r,s,0) < 0.5) continue;
	  n++;
	  segid = MRIgetVoxVal(gtm->gtmseg,c,r,s,0);
	  tt = ct->entries[segid]->TissueType;
	  if(tt != 1 && tt != 2) continue; // not GM (hardcoded)
	  sum += ((double)gtm->res->rptr[n+1][f+1]*gtm->res->rptr[n+1][f+1]);
	  nhits ++;
	}
      }
    }
    gtm->rvargm->rptr[1][f+1] = sum/nhits;
    printf("rvargm %2d %6.4f\n",f,gtm->rvargm->rptr[1][f+1]);
  }

  sprintf(tmpstr,"%s/rvar.gm.dat",gtm->OutDir);
  fp = fopen(tmpstr,"w");
  for(f=0; f < gtm->res->cols; f++) fprintf(fp,"%30.20f\n",gtm->rvargm->rptr[1][f+1]);
  fclose(fp);

  return(0);
}






