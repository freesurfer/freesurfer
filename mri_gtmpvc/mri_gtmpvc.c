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
 *    $Date: 2014/07/11 19:18:13 $
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


// $Id: mri_gtmpvc.c,v 1.27 2014/07/11 19:18:13 greve Exp $

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
#include "pdf.h"

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

static char vcid[] = "$Id: mri_gtmpvc.c,v 1.27 2014/07/11 19:18:13 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

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
char *OutDir=NULL,*AuxDir,*RBVVolFile=NULL;
char *OutBetaFile=NULL,*OutBetaVarFile=NULL,*OutXtXFile=NULL;
char *SrcBetaFile=NULL;
int SynthOnly=0,SaveSynth=0;
int GTMSynthSeed=0,GTMSynthReps=1;
double SynthPSFFWHMCol=0,SynthPSFFWHMRow=0,SynthPSFFWHMSlice=0;
double SynthPSFStdCol=0,SynthPSFStdRow=0,SynthPSFStdSlice=0;
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
int SaveRBVSeg=0;

GTM *gtm;
GTMOPT *gtmopt;
GTMOPT *gtmopt_powell;

LTA *LTAapplyAffineParametersTKR(LTA *inlta, const float *p, const int np, LTA *outlta);
int DoOpt=0;
char *SUBJECTS_DIR;

int DoRBV=0;
int AutoMask=0;
float pxfm[6]={0,0,0,0,0,0};
int ApplyXFM=0;
int DoSimulation = 0;
int DoVoxFracCorTmp;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err,c,f,n;
  double vrfmin,vrfmax,vrfmean;
  struct timeb  mytimer, timer;
  MRI *gtmres;
  FILE *logfp,*fp;
  char *stem;
  LTA *ltatmp;

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
  gtm->rbvsegres = 0;

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

  DoVoxFracCorTmp = gtm->DoVoxFracCor;

#ifdef _OPENMP
  printf("%d avail.processors, using %d\n",omp_get_num_procs(),omp_get_max_threads());
#endif

  printf("Creating output directory %s\n",OutDir);
  err = mkdir(OutDir,0777);
  if(err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",OutDir);
    perror(NULL);
    return(1);
  }
  sprintf(logfile,"%s/mri_gtmpvc.log",OutDir);
  logfp = fopen(logfile,"w");
  gtm->logfp = logfp;
  dump_options(logfp);

  err = mkdir(AuxDir,0777);
  if(err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",AuxDir);
    perror(NULL);
    return(1);
  }
  gtm->AuxDir = AuxDir;

  TimerStart(&timer);

  // Load seg
  TimerStart(&mytimer);
  printf("Loading seg for gtm %s\n",SegVolFile);fflush(stdout);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  fprintf(logfp,"Loading seg for gtm %s\n",SegVolFile);fflush(logfp);
  gtm->anatseg = MRIread(SegVolFile);
  if(gtm->anatseg==NULL) exit(1);
  if(Gdiag_no > 0)printf("  done loading seg\n");fflush(stdout);
  if(Gdiag_no > 0) PrintMemUsage(stdout);

  stem = IDstemFromName(SegVolFile);
  sprintf(tmpstr,"%s.ctab",stem);
  printf("Loading seg ctab %s\n",tmpstr);fflush(stdout);
  fprintf(logfp,"Loading ctab %s\n",tmpstr);fflush(logfp);
  gtm->ctGTMSeg = CTABreadASCII(tmpstr);
  if(gtm->ctGTMSeg == NULL) exit(1);
  if(Gdiag_no > 0) printf("  done loading ctab\n");fflush(stdout);

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
    printf("Replacing %d\n",gtm->nReplace);fflush(stdout);
    mritmp = MRIreplaceList(gtm->anatseg, gtm->SrcReplace, gtm->TrgReplace, gtm->nReplace, NULL);
    if(Gdiag_no > 0) printf("  done replacing\n");fflush(stdout);
    MRIfree(&gtm->anatseg);
    gtm->anatseg = mritmp;
    sprintf(tmpstr,"%s/seg.replace.list",AuxDir);
    fp = fopen(tmpstr,"w");
    GTMprintReplaceList(fp, gtm->nReplace, gtm->SrcReplace, gtm->TrgReplace);
    fclose(fp);
  }
  if(CheckSegTissueType(gtm->anatseg, gtm->ctGTMSeg)){
    printf("Failed tissue type check\n");
    fprintf(logfp,"Failed tissue type check\n");
    exit(1);
  }
  printf("Pruning ctab\n"); fflush(stdout);
  gtm->ctGTMSeg = CTABpruneCTab(gtm->ctGTMSeg, gtm->anatseg);
  if(Gdiag_no > 0) printf("  done pruning ctab\n"); fflush(stdout);

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
    printf("Computing auto mask \n");fflush(stdout);
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    PrintMemUsage(logfp);
    GTMautoMask(gtm);
    if(Gdiag_no > 0) printf("  done auto mask \n");fflush(stdout);
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    PrintMemUsage(logfp);
    sprintf(tmpstr,"%s/mask.nii.gz",AuxDir);
    MRIwrite(gtm->mask,tmpstr);
  }
  else  gtm->mask=NULL; // should already be NULL

  printf("Data load time %4.1f sec\n",TimerStop(&mytimer)/1000.0);
  fprintf(logfp,"Data load time %4.1f sec\n",TimerStop(&mytimer)/1000.0);
  fflush(stdout);fflush(logfp);

  GTMsetNMask(gtm);
  GTMsegidlist(gtm);

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
  sprintf(tmpstr,"%s/seg.ctab",AuxDir);
  CTABwriteFileASCIItt(gtm->ctGTMSeg,tmpstr);

  printf("Checking Ref Ids\n");
  err=GTMcheckRefIds(gtm);
  if(err) exit(1);
  GTMprintRefIds(gtm, stdout);
  GTMprintRefIds(gtm, logfp);

  /* Create the "SegPVF". This is a multi-frame volume, each frame corresponds to a
     different Seg ID. The value is the PVF of that SegID. This is independent of
     the PSF (but accounts for the volume fraction effect). */
  printf("Computing Seg PVF \n");fflush(stdout);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  PrintMemUsage(logfp);
  gtm->segpvf = MRIseg2SegPVF(gtm->anatseg, gtm->seg2pet, 0.5, gtm->segidlist, 
			      gtm->nsegs, gtm->mask, 0, NULL, gtm->segpvf);
  if(gtm->segpvf==NULL) exit(1);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  PrintMemUsage(logfp);

  /* This creates a segmentation in the PET space based upon which Seg
     has the greated PVF (independent of PSF). (used by GTMsynth) */
  printf("Computing Seg in input space \n");fflush(stdout);
  gtm->gtmseg = MRIsegPVF2Seg(gtm->segpvf, gtm->segidlist, gtm->nsegs, 
			      gtm->ctGTMSeg, gtm->mask, gtm->gtmseg);
  if(OutSegFile){
    err=MRIwrite(gtm->gtmseg,OutSegFile);
    if(err) exit(1);
  }

  // Create GTM matrix
  if(DoSimulation && !gtm->DoVoxFracCor) {
    // DoVoxFracCorTmp keeps track of the old value of gtm->DoVoxFracCor
    printf("Turning VoxFracCor ON for simulation. It will be turned OFF for anlaysis.\n");
    fprintf(logfp,"Turning VoxFracCor ON for simulation. It will be turned OFF for anlaysis.\n");
    gtm->DoVoxFracCor = 1;
  }
  printf("Building GTM DoVoxFracCor=%d\n",gtm->DoVoxFracCor);fflush(stdout); 
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  PrintMemUsage(logfp);
  TimerStart(&mytimer) ;
  GTMbuildX(gtm);
  if(gtm->X==NULL) exit(1);
  printf(" gtm build time %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  fprintf(logfp,"GTM-Build-time %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(logfp);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  PrintMemUsage(logfp);

  // Create GTM pvf in pet space (why?)
  gtm->ttpvf = MRIsegPVF2TissueTypePVF(gtm->segpvf, gtm->segidlist, gtm->nsegs, 
				       gtm->ctGTMSeg, gtm->mask, gtm->ttpvf);
  sprintf(tmpstr,"%s/pvf.nii.gz",AuxDir);
  MRIwrite(gtm->ttpvf,tmpstr);

  sprintf(tmpstr,"%s/hrseg2pet.lta",AuxDir);
  LTAwrite(gtm->seg2pet,tmpstr);
  sprintf(tmpstr,"%s/anat2pet.lta",AuxDir);
  LTAwrite(gtm->anat2pet,tmpstr);
  ltatmp = LTAinvert(gtm->anat2pet,NULL);
  sprintf(tmpstr,"%s/pet2anat.lta",AuxDir);
  LTAwrite(ltatmp,tmpstr);
  LTAfree(&ltatmp);
  //sprintf(tmpstr,"%s/anat2hrseg.lta",OutDir);
  //LTAwrite(gtm->anat2seg,tmpstr);

  if(DoSimulation){
    printf("Synthsizing using supplied beta %s\n",SrcBetaFile);
    gtm->beta = srcbeta;
    GTMsynth(gtm,GTMSynthSeed,GTMSynthReps);
    printf("Smoothing synthsized %g %g %g\n",SynthPSFStdCol,SynthPSFStdRow,SynthPSFStdSlice);
    gtm->ysynthsm = MRIgaussianSmoothNI(gtm->ysynth,SynthPSFStdCol,SynthPSFStdRow,SynthPSFStdSlice,NULL);
    if(SaveSynth){
      sprintf(tmpstr,"%s/synth.nii.gz",OutDir);
      MRIwrite(gtm->ysynthsm,tmpstr);
    }
    if(SynthOnly){
      printf("SynthOnly requested so exiting now\n");
      printf("mri_gtmpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
      fprintf(logfp,"SynthOnly requested so exiting now\n");
      fprintf(logfp,"mri_gtmpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
      exit(0);
    }
    MRIfree(&gtm->yvol);
    gtm->yvol = MRIcopy(gtm->ysynthsm,NULL);
    MRIfree(&gtm->ysynth);
    MRIfree(&gtm->ysynthsm);
    MatrixFree(&gtm->beta);
    if(!DoVoxFracCorTmp) {
      // DoVoxFracCorTmp keeps track of the old value of gtm->DoVoxFracCor
      printf("Turning VoxFracCor back OFF for anlaysis.\n");
      fprintf(logfp,"Turning VoxFracCor back OFF for anlaysis.\n");
      gtm->DoVoxFracCor = 0;
      printf("Rebuilding GTM DoVoxFracCor=%d\n",gtm->DoVoxFracCor);
      fprintf(logfp,"Rebuilding GTM DoVoxFracCor=%d\n",gtm->DoVoxFracCor);
      if(Gdiag_no > 0) PrintMemUsage(stdout);
      PrintMemUsage(logfp);
      TimerStart(&mytimer);
      MatrixFree(&gtm->X);
      MatrixFree(&gtm->X0);
      GTMbuildX(gtm);
      if(gtm->X==NULL) exit(1);
      printf(" gtm build time %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
      fprintf(logfp,"GTM-rebuild-time %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(logfp);
      if(Gdiag_no > 0) PrintMemUsage(stdout);
      PrintMemUsage(logfp);
      // Create GTM pvf in pet space (why?)
      gtm->ttpvf = MRIsegPVF2TissueTypePVF(gtm->segpvf, gtm->segidlist, gtm->nsegs, 
					   gtm->ctGTMSeg, gtm->mask, gtm->ttpvf);
      sprintf(tmpstr,"%s/pvf.nii.gz",AuxDir);
      MRIwrite(gtm->ttpvf,tmpstr);
    }
  }

  //printf("Freeing segpvf\n"); fflush(stdout);
  //MRIfree(&gtm->segpvf);
  if(SaveX0) {
    printf("Writing X0 to %s\n",Xfile);
    MatlabWrite(gtm->X0, X0file,"X0");
  }
  if(SaveX) {
    printf("Writing X to %s\n",Xfile);
    MatlabWrite(gtm->X, Xfile,"X");
  }

  printf("Solving ...\n");
  GTMmatrixY(gtm);
  TimerStart(&mytimer) ; 
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  PrintMemUsage(logfp);
  err=GTMsolve(gtm); // also rescales everything if desired
  if(err) exit(1);
  printf("Time to solve %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
  fprintf(logfp,"GTM-Solve-Time %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(logfp);
  if(Gdiag_no > 0) PrintMemUsage(stdout);
  PrintMemUsage(logfp);

  printf("Freeing X\n");
  MatrixFree(&gtm->X);

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

  if(gtm->rescale){
    // Rescaling is done during GTMsolve()
    printf("rescale factor %20.15lf\n",gtm->scale);
    fprintf(logfp,"rescale factor %20.15lf\n",gtm->scale);
    sprintf(tmpstr,"%s/scale.dat",AuxDir);
    fp = fopen(tmpstr,"w");
    fprintf(fp,"%20.15lf\n",gtm->scale);
    fclose(fp);
  }

  GTMsegrvar(gtm);
  VRFStats(gtm, &vrfmean, &vrfmin, &vrfmax);

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
    GTMsynth(gtm,GTMSynthSeed,GTMSynthReps);
    printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
    fprintf(logfp,"GTMSynth: Seed=%d, Reps=%d\n",GTMSynthSeed,GTMSynthReps);
  }
  if(yhat0File) MRIwrite(gtm->ysynth,yhat0File);
  
  printf("Freeing X0\n");
  MatrixFree(&gtm->X0);

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
  sprintf(tmpstr,"%s/nmask.dat",AuxDir);
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
  sprintf(tmpstr,"%s/cond.dat",AuxDir);
  fp = fopen(tmpstr,"w");
  fprintf(fp,"%20.15lf\n",gtm->XtXcond);
  fclose(fp);

  // Res variance in each seg
  sprintf(tmpstr,"%s/seg.rvar.nii.gz",AuxDir);
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
  sprintf(tmpstr,"%s/seg.vrf.nii.gz",AuxDir);
  printf("Writing seg vrf to %s\n",tmpstr);
  mritmp = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, 1);
  for(c=0; c < gtm->nsegs; c++) MRIsetVoxVal(mritmp,c,0,0,0, gtm->vrf->rptr[c+1][1]);
  err=MRIwrite(mritmp,tmpstr);
  if(err) exit(1);
  MRIfree(&mritmp);

  // NVox in each seg
  sprintf(tmpstr,"%s/seg.nvox.nii.gz",AuxDir);
  printf("Writing seg nvox to %s\n",tmpstr);
  mritmp = MRIallocSequence(gtm->nsegs, 1, 1, MRI_FLOAT, 1);
  for(c=0; c < gtm->nsegs; c++) MRIsetVoxVal(mritmp,c,0,0,0, gtm->nvox->rptr[c+1][1]);
  err=MRIwrite(mritmp,tmpstr);
  if(err) exit(1);
  MRIfree(&mritmp);

  if(gtm->DoMGPVC){
    printf("Performing MG PVC\n");
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    sprintf(tmpstr,"%s/mg.nii.gz",OutDir);
    MGPVCFile = strcpyalloc(tmpstr);
    fprintf(logfp,"MG PVC\n");
    GTMmgpvc(gtm);
    err = MRIwrite(gtm->mg,MGPVCFile);
    if(err) exit(1);
    if(Gdiag_no > 0) printf("done with mgpvc\n");
    fprintf(logfp,"done with mgpvc\n");
    sprintf(tmpstr,"%s/mg.reftac.dat",AuxDir);
    GTMwriteMGRefTAC(gtm, tmpstr);
    sprintf(tmpstr,"%s/gm.pvf.psf.nii.gz",AuxDir);
    MRIwrite(gtm->gmpvfpsf, tmpstr);
  }
    
  if(DoRBV){
    printf("Computing RBV Seg\n");
    GTMrbvseg(gtm);
    sprintf(tmpstr,"%s/rbv.nii.gz",OutDir);
    RBVVolFile = strcpyalloc(tmpstr);
    printf("Computing RBV\n");
    if(Gdiag_no > 0) PrintMemUsage(stdout);
    GTMrbv(gtm);
    printf("Writing output to %s ...",RBVVolFile);fflush(stdout); TimerStart(&mytimer) ;
    err = MRIwrite(gtm->rbv,RBVVolFile);
    if(err){
      printf("ERROR: writing to %s\n",RBVVolFile);
      exit(1);
    }
    printf(" %4.1f sec\n",TimerStop(&mytimer)/1000.0);fflush(stdout);
    sprintf(tmpstr,"%s/aux/anat2rbv.lta",OutDir);
    LTAwrite(gtm->anat2rbv,tmpstr);
    sprintf(tmpstr,"%s/aux/anat2rbv.lta",OutDir);
    LTAwrite(gtm->anat2rbv,tmpstr);
    ltatmp = LTAinvert(gtm->anat2rbv,NULL);
    sprintf(tmpstr,"%s/rbv2anat.lta",AuxDir);
    LTAwrite(ltatmp,tmpstr);
    LTAfree(&ltatmp);

    sprintf(tmpstr,"%s/rbv.segmean.nii.gz",AuxDir);
    MRIwrite(gtm->rbvsegmean,tmpstr);
    if(SaveRBVSeg){
      sprintf(tmpstr,"%s/seg.rbv.nii.gz",AuxDir);
      MRIwrite(gtm->rbvsegmasked,tmpstr);
    }
  }

  // Free the data from gtm->yvol, keep header
  // Not very useful here, but have to wait until after RBV and MG
  printf("Freeing y\n"); fflush(stdout);
  mritmp = MRIcopyHeader(gtm->yvol,NULL);
  MRIcopyPulseParameters(gtm->yvol,mritmp);
  MRIfree(&gtm->yvol);
  gtm->yvol = mritmp;
  
  if(gtm->nContrasts > 0){
    printf("Testing %d contrasts\n",gtm->nContrasts);
    GTMttest(gtm);
    printf("Writing contrasts\n");
    GTMwriteContrasts(gtm);
  }
  PrintMemUsage(stdout);
  PrintMemUsage(logfp);

  fprintf(logfp,"mri_gtmpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
  fprintf(logfp,"mri_gtmpvc done\n");
  fclose(logfp);
  printf("mri_gtmpvc-runtime %5.2f min\n",TimerStop(&timer)/60000.0);
  printf("mri_gtmpvc done\n");
  return(0);
  exit(0);
} // end of main
/*--------------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused, n;
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
    else if(!strcasecmp(option, "--default-seg-merge-choroid")){
      GTMdefaultSegReplacmentList(&gtm->nReplace,&(gtm->SrcReplace[0]),&(gtm->TrgReplace[0]));
      // Last two itmes are choroid.
      gtm->nReplace -= 2;
    }

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
    else if(!strcmp(option, "--xfm")){
      if(nargc < 6) CMDargNErr(option,6);
      for(n=0; n < 6; n++) sscanf(pargv[n],"%f",&pxfm[n]);
      ApplyXFM=1;
      nargsused = 6;
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
    else if(!strcasecmp(option, "--rbv-res")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&gtm->rbvsegres);
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
      sprintf(tmpstr,"%s/aux",OutDir);
      AuxDir = strcpyalloc(tmpstr);
      gtm->OutDir = OutDir;
      sprintf(tmpstr,"%s/rvar.dat",AuxDir);
      RVarFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/skew.dat",AuxDir);
      SkewFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/kurtosis.dat",AuxDir);
      KurtosisFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/gtm.stats.dat",OutDir);
      VRFStatsFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/beta.nii.gz",OutDir);
      OutBetaFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/beta.var.nii.gz",OutDir);
      OutBetaVarFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/seg.nii.gz",AuxDir);
      OutSegFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/xtx.mat",AuxDir);
      OutXtXFile = strcpyalloc(tmpstr);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--save-eres")) SaveEres=1;
    else if(!strcasecmp(option, "--save-yhat")) SaveYhat=1;
    else if(!strcasecmp(option, "--save-yhat-with-noise")){
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&GTMSynthSeed);
      if(GTMSynthSeed <= 0) GTMSynthSeed = PDFtodSeed();
      sscanf(pargv[1],"%d",&GTMSynthReps);
      SaveYhat=1;
      nargsused = 2;
    }
    else if(!strcasecmp(option, "--save-yhat0")) SaveYhat0=1;
    else if(!strcasecmp(option, "--save-yhat-full-fov")) SaveYhatFullFoV=1;
    else if(!strcasecmp(option, "--save-rbv-seg")) SaveRBVSeg = 1;
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
      if(nargc < 6) CMDargNErr(option,6);
      SrcBetaFile = pargv[0];
      sscanf(pargv[1],"%lf",&SynthPSFFWHMCol);
      sscanf(pargv[2],"%lf",&SynthPSFFWHMRow);
      sscanf(pargv[3],"%lf",&SynthPSFFWHMSlice);
      sscanf(pargv[4],"%d",&GTMSynthSeed);
      if(GTMSynthSeed < 0) GTMSynthSeed = PDFtodSeed();
      sscanf(pargv[5],"%d",&GTMSynthReps);
      mritmp = MRIread(SrcBetaFile);
      if(mritmp == NULL) exit(1);
      MATRIX *srcbetaT = fMRItoMatrix(mritmp,NULL);
      srcbeta = MatrixTranspose(srcbetaT,NULL);
      MatrixFree(&srcbetaT);
      SynthPSFStdCol   = SynthPSFFWHMCol/sqrt(log(256.0));
      SynthPSFStdRow   = SynthPSFFWHMRow/sqrt(log(256.0));
      SynthPSFStdSlice = SynthPSFFWHMSlice/sqrt(log(256.0));
      DoSimulation = 1;
      nargsused = 6;
    } 
    else if(!strcasecmp(option, "--synth-only")) {SynthOnly = 1;SaveSynth = 1;}
    else if(!strcasecmp(option, "--save-synth")) SaveSynth = 1;
    else if(!strcasecmp(option, "--synth-save")) SaveSynth = 1;
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
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
    else if(!strcasecmp(option, "--max-threads-1") || !strcasecmp(option, "--max-threads-minus-1")){
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
  printf("   --default-seg-merge-choroid : default schema for merging ROIs but keeps choroid\n");
  printf("\n");
  printf("   --no-vox-frac-cor : do not use voxel fraction correction (with --psf 0 turns off PVC entirely)\n");
  printf("   --rbv             : perform RBV PVC\n");
  printf("   --rbv-res voxsize : set RBV voxel resolution (good for when standard res takes too much memory)\n");
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
  printf("   --save-yhat-with-noise seed nreps : saves yhat with noise, seed < 0 for TOD\n");
  printf("   --save-yhat-full-fov : saves yhat in full FoV (if FoV was reduced)\n");
  printf("   --save-yhat0 : saves yhat prior to smoothing\n");
  printf("\n");
  printf("   --synth gtmbeta C R S seed nreps : synthesize volume with gtmbeta as input\n");
  printf("       spec all other inputs the same; CRS are PSF for col, row, slice\n");
  printf("       seed=0 for no noise, -1 for TOD seed. Turns on VFC for synth\n");
  printf("       but returns it to its specified value for analysis\n");
  printf("   --synth-only : exit after doing synthesis (implies --synth-save)\n");
  printf("   --synth-save : with --synth saves synthesized volume to outdir/synth.nii.gz\n");
  printf("\n");
  #ifdef _OPENMP
  printf("   --threads N : use N threads (with Open MP)\n");
  printf("   --max-threads : use the maximum allowable number of threads for this computer\n");
  printf("   --max-threads-minus-1 : use one less than the maximum allowable number of threads for this computer\n");
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

  sprintf(tmpstr,"%s/%s/mri/orig.mgz",SUBJECTS_DIR,gtm->anat2pet->subject);
  gtm->anatconf = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
  if(gtm->anatconf == NULL) exit(1);

  if(gtm->cFWHM < 0 || gtm->rFWHM < 0 || gtm->sFWHM < 0){
    printf("ERROR: must spec psf FWHM\n");
    exit(1);
  }

  if(OutDir == NULL){
    printf("ERROR: must spec an output folder with --o\n");
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
    omp_set_num_threads(nthreads);fflush(stdout);
    #endif
  }
  printf("Loading input %s\n",SrcVolFile);fflush(stdout);
  gtm->yvol = MRIread(SrcVolFile);
  if(gtm->yvol==NULL) exit(1);
  gtm->nframes = gtm->yvol->nframes ;
  printf("  done loading input %d frames\n",gtm->nframes);fflush(stdout);

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
  if(ApplyXFM){
    printf("Applying xfm parameters\n");
    LTAapplyAffineParametersTKR(gtm->anat2pet, pxfm, 6, gtm->anat2pet);
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
  int n;
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
  fprintf(fp,"ApplyXFM %d: ",ApplyXFM);
  for(n=0; n < 6; n++) fprintf(fp,"%6.4f ",pxfm[n]);
  fprintf(fp,"\n");


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


