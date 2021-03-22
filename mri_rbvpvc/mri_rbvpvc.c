/**
 * @brief Peforms RBV Partial volume correction
 *
 * Implementation of Region-based Voxelwise (RBV) partial volume correction
 * as found in Thomas, et al, 2011, Eur J Nucl Med Mol Imaging, 38:1104-1119.
 * It also implements the Geometric Transfer Matrix (GTM) as it is needed by RBV.
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

/*!
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/



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
#include "romp_support.h"
#endif

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

typedef struct 
{
  MRI *yvol;
  MRI *anatseg;
  LTA *anat2pet;
  LTA *anat2seg;
  LTA *seg2anat;
  LTA *seg2pet;
  MRI *mask;
  double cFWHM, rFWHM, sFWHM;
  double PadThresh;
  MRI *yseg,*yhat0seg,*yhatseg;
  int rescale,n_scale_refids,scale_refids[100];
  double scale;
  MRI *segpvf;
  MRI *ttpvf;
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
  MATRIX *skew,*kurtosis;
  double cStd, rStd, sStd;
  MRI *ysynth,*ysynthsm;
  int *nperseg;
  MATRIX *nvox,*vrf,*segrvar;
  MRI *rbv;
  int mask_rbv_to_brain;
  MRI *mg;
  double mg_gmthresh;
  char *mg_ref_schema;
  int n_mg_refids,mg_refids[100];
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
MRI *GTMsegSynth(GTM *gtm);
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
char *SrcBetaFile=NULL,*SynthFile=NULL;
int SynthOnly=0;
MATRIX *srcbeta;
double psfFWHM=-1;
char tmpstr[5000],logfile[5000];
char *PVFFile=NULL, *SegTTypeFile=NULL;
double ApplyFWHM=0;
char *Xfile=NULL;
char *VRFStatsFile=NULL;
char *eresFile=NULL, *yhatFile=NULL, *yhat0File=NULL;
char *OutSegFile=NULL;
MRI *mritmp;
char *RVarFile=NULL,*SkewFile=NULL,*KurtosisFile=NULL;
int RVarOnly=0;
int nthreads=1;
int nReplace = 0, SrcReplace[1000], TrgReplace[1000];
int ttReduce = 0;
char *MGPVCFile=NULL;
char *regfile;
int regidentity = 0;
int regtype;
int SaveEres=0, SaveYhat=0,SaveYhat0=0;

int VRFStats(GTM *gtm, double *vrfmean, double *vrfmin, double *vrfmax);
int WriteVRFStats(char *fname, GTM *gtm);


MATRIX *GTMttest(MATRIX *beta,MATRIX *iXtX, double *rvar, int nthseg);
int GTMmgRefIds(GTM *gtm);
int GTMprintMGRefTAC(GTM *gtm, FILE *fp);
int GTMwriteMGRefTAC(GTM *gtm, char *filename);
int GTMrescale(GTM *gtm);

GTM *gtm;
GTMOPT *gtmopt;
GTMOPT *gtmopt_powell;

LTA *LTAapplyAffineParametersTKR(LTA *inlta, const float *p, const int np, LTA *outlta);
int DoOpt=0;
char *SUBJECTS_DIR;
int CheckX(MATRIX *X);
int dngtest(LTA *aseg2vol);

int DoMGPVC=0, DoRBV=0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err,c,f,nTT,n;
  double vrfmin,vrfmax,vrfmean;
  Timer mytimer, timer;
  MRI *gtmres;
  FILE *logfp,*fp;
  char *stem;

  nargs = handleVersionOption(argc, argv, "mri_rbvpvc");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  vg_isEqual_Threshold = 10e-4;

  gtm = GTMalloc();
  gtm->ctGTMSeg = TissueTypeSchema(NULL,"default-jan-2014");
  gtm->mg_ref_schema = "basic-or-lobes";
  gtm->rescale = 1; 
  gtm->n_scale_refids = 2;
  gtm->scale_refids[0] = 7;
  gtm->scale_refids[1] = 46;
  gtm->mask_rbv_to_brain = 1;

  GTMmgRefIds(gtm);

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
    sprintf(logfile,"%s/mri_rbvpvc.log",OutDir);
    logfp = fopen(logfile,"w");
  }
  else logfp = stdout;
  dump_options(logfp);

  timer.reset();

  stem = IDstemFromName(SegVolFile);
  sprintf(tmpstr,"%s.lta",stem);
  gtm->anat2seg = LTAread(tmpstr);
  if(gtm->anat2seg == NULL) exit(1);
  gtm->seg2anat = LTAinvert(gtm->anat2seg,NULL);
  gtm->seg2pet = LTAconcat2(gtm->seg2anat,gtm->anat2pet,1);
  if(gtm->seg2pet == NULL) {
    printf("ERROR: LTAconcat()\n");
    printf("mri_rbvpvc exited with errors\n");
    exit(1);
  }

  // Load seg
  mytimer.reset();
  printf("Loading seg for gtm %s\n",SegVolFile);fflush(stdout);
  fprintf(logfp,"Loading seg for gtm %s\n",SegVolFile);fflush(logfp);
  gtm->anatseg = MRIread(SegVolFile);
  if(gtm->anatseg==NULL) exit(1);
  if(nReplace > 0) {
    printf("Replacing %d\n",nReplace);
    for(f=0; f < nReplace; f++) printf("%2d:  %4d %4d\n",f+1,SrcReplace[f],TrgReplace[f]);
    for(f=0; f < nReplace; f++) fprintf(logfp,"%2d:  %4d %4d\n",f+1,SrcReplace[f],TrgReplace[f]);
    mritmp = MRIreplaceList(gtm->anatseg, SrcReplace, TrgReplace, nReplace, NULL, NULL);
    MRIfree(&gtm->anatseg);
    gtm->anatseg = mritmp;
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

  // could check dims against LTA

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
  else gtm->mask=NULL; // should already be NULL

  printf("Data load time %4.1f sec\n",mytimer.seconds());
  fprintf(logfp,"Data load time %4.1f sec\n",mytimer.seconds());
  fflush(stdout);fflush(logfp);

  GTMsetNMask(gtm);
  GTMsegidlist(gtm);
  GTMmatrixY(gtm);
  GTMpsfStd(gtm); 
  GTMnPad(gtm);   
  printf("nmask = %d, nsegs = %d, excluding segid=0\n",gtm->nmask,gtm->nsegs);
  printf("FWHM: %g %g %g\n",gtm->cFWHM,gtm->rFWHM,gtm->sFWHM);
  printf("Std:  %g %g %g\n",gtm->cStd,gtm->rStd,gtm->sStd);
  printf("nPad %d, PadThresh %g\n",gtm->nPad,gtm->PadThresh);
  fprintf(logfp,"nmask = %d, nsegs = %d, excluding segid=0\n",gtm->nmask,gtm->nsegs);
  fprintf(logfp,"FWHM: %g %g %g\n",gtm->cFWHM,gtm->rFWHM,gtm->sFWHM);
  fprintf(logfp,"Std:  %g %g %g\n",gtm->cStd,gtm->rStd,gtm->sStd);
  fprintf(logfp,"nPad %d, PadThresh %g\n",gtm->nPad,gtm->PadThresh);

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
  mytimer.reset() ;
  GTMbuildX(gtm);
  printf(" %4.1f sec\n",mytimer.seconds());fflush(stdout);
  if(gtm->X==NULL) exit(1);
  fprintf(logfp,"GTM-Build-time %4.1f sec\n",mytimer.seconds());fflush(logfp);

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
      printf("mri_rbvpvc-runtime %5.2f min\n",timer.minutes());
      exit(0);
    }
    MRIfree(&gtm->yvol);
    GTMsmoothSynth(gtm);
    gtm->yvol = MRIcopy(gtm->ysynthsm,NULL);
  }

  if(Xfile) {
    printf("Writing X to %s\n",Xfile);
    MatrixWriteTxt(Xfile, gtm->X);
  }
  printf("Solving ...\n");
  mytimer.reset() ; 
  GTMsolve(gtm);
  printf("Time to solve %4.1f sec\n",mytimer.seconds());fflush(stdout);
  fprintf(logfp,"GTM-Solve-Time %4.1f sec\n",mytimer.seconds());fflush(logfp);

  sprintf(tmpstr,"%s/hrseg2pet.lta",OutDir);
  LTAwrite(gtm->seg2pet,tmpstr);
  sprintf(tmpstr,"%s/anat2pet.lta",OutDir);
  LTAwrite(gtm->anat2pet,tmpstr);
  //sprintf(tmpstr,"%s/anat2hrseg.lta",OutDir);
  //LTAwrite(gtm->anat2seg,tmpstr);

  sprintf(tmpstr,"%s/seg.ctab",OutDir);
  CTABwriteFileASCIItt(gtm->ctGTMSeg,tmpstr);

  if(gtm->rescale){
    // Rescaling is done during GTMsolve()
    printf("Rescaling using segids ");
    for(n=0; n < gtm->n_scale_refids; n++) printf("%4d ",gtm->scale_refids[n]);
    printf("\n");
    fprintf(logfp,"Rescaling using segids ");
    for(n=0; n < gtm->n_scale_refids; n++) fprintf(logfp,"%4d ",gtm->scale_refids[n]);
    fprintf(logfp,"\n");
    fprintf(logfp,"rescale factor %20.15lf\n",gtm->scale);
    sprintf(tmpstr,"%s/scale.dat",OutDir);
    fp = fopen(tmpstr,"w");
    fprintf(fp,"%20.15lf\n",gtm->scale);
    fclose(fp);
  }

  // Create GTM pvf in pet space (why?)
  gtm->ttpvf = MRIsegPVF2TissueTypePVF(gtm->segpvf, gtm->segidlist, gtm->nsegs, 
				       gtm->ctGTMSeg, gtm->mask, gtm->ttpvf);
  sprintf(tmpstr,"%s/pvf.nii.gz",OutDir);
  MRIwrite(gtm->ttpvf,tmpstr);

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

  printf("rvar = %g\n",gtm->rvar->rptr[1][1]);
  fprintf(logfp,"rvar = %g\n",gtm->rvar->rptr[1][1]);
  gtm->XtXcond = MatrixConditionNumber(gtm->XtX);
  printf("XtX  Condition     %8.3f \n",gtm->XtXcond);
  fprintf(logfp,"XtX  Condition     %8.3f \n",gtm->XtXcond);

  if(yhat0File || yhatFile){
    printf("Synthesizing ... ");fflush(stdout); mytimer.reset() ;
    GTMsynth(gtm);
    printf(" %4.1f sec\n",mytimer.seconds());fflush(stdout);
  }
  if(yhat0File) MRIwrite(gtm->ysynth,yhat0File);

  if(yhatFile){
    printf("Smoothing synthesized ... ");fflush(stdout); mytimer.reset() ;
    GTMsmoothSynth(gtm);
    printf(" %4.1f sec\n",mytimer.seconds());fflush(stdout);
    MRIwrite(gtm->ysynthsm,yhatFile);
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
    printf("mri_rbvpvc-runtime %5.2f min\n",timer.minutes());
    exit(0);
  }

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

  GTMsegrvar(gtm);
  VRFStats(gtm, &vrfmean, &vrfmin, &vrfmax);

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

  if(DoMGPVC){
    printf("Performing MG PVC, refidschema %s\n",gtm->mg_ref_schema);
    sprintf(tmpstr,"%s/mg.nii.gz",OutDir);
    MGPVCFile = strcpyalloc(tmpstr);
    fprintf(logfp,"MG PVC\n");
    GTMmgpvc(gtm);
    err = MRIwrite(gtm->mg,MGPVCFile);
    if(err) exit(1);
    printf("done with mgpvc\n");
    fprintf(logfp,"done with mgpvc\n");
    fprintf(logfp,"MG Ref Schema %s\n",gtm->mg_ref_schema);
    for(c=0; c < gtm->n_mg_refids; c++) fprintf(logfp,"%4d ",gtm->mg_refids[c]);
    fprintf(logfp,"\n");
    stem = IDstemFromName(MGPVCFile);
    sprintf(tmpstr,"%s.reftac.dat",stem);
    GTMwriteMGRefTAC(gtm, tmpstr);
  }
    
  if(DoRBV){
    sprintf(tmpstr,"%s/rbv.nii.gz",OutDir);
    RBVVolFile = strcpyalloc(tmpstr);
    printf("Computing RBV\n");
    GTMrbv(gtm);
    printf("Writing output to %s ...",RBVVolFile);fflush(stdout); mytimer.reset() ;
    err = MRIwrite(gtm->rbv,RBVVolFile);
    if(err){
      printf("ERROR: writing to %s\n",RBVVolFile);
      exit(1);
    }
    printf(" %4.1f sec\n",mytimer.seconds());fflush(stdout);
  }


  fprintf(logfp,"mri_rbvpvc-runtime %5.2f min\n",timer.minutes());
  fprintf(logfp,"mri_rbvpvc done\n");
  fclose(logfp);
  printf("mri_rbvpvc-runtime %5.2f min\n",timer.minutes());
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
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if(!strcasecmp(option, "--help"))  print_help() ;
    else if(!strcasecmp(option, "--version")) print_version() ;
    else if(!strcasecmp(option, "--debug"))   debug = 1;
    else if(!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if(!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if(!strcasecmp(option, "--opt")) DoOpt=1;
    else if(!strcasecmp(option, "--ttype+head"))
      gtm->ctGTMSeg = TissueTypeSchema(NULL,"default-jan-2014+head");
    else if(!strcasecmp(option, "--tt-reduce")) ttReduce = 1;
    else if(!strcasecmp(option, "--no-mask_rbv_to_brain")) gtm->mask_rbv_to_brain = 0;

    else if(!strcmp(option, "--sd") || !strcmp(option, "-SDIR")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
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
    else if(!strcasecmp(option, "--rescale")) gtm->rescale = 1;
    else if(!strcasecmp(option, "--no-rescale")) gtm->rescale = 0;

    else if(!strcasecmp(option, "--mgpvc")) {
      if(nargc < 1) CMDargNErr(option,1);
      DoMGPVC = 1;
      sscanf(pargv[0],"%lf",&gtm->mg_gmthresh);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--X")) {
      if(nargc < 1) CMDargNErr(option,1);
      Xfile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--rvar-only")) RVarOnly = 1;
    else if(!strcasecmp(option, "--o")) {
      if(nargc < 1) CMDargNErr(option,1);
      OutDir = pargv[0];
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
      sprintf(tmpstr,"%s/seg.nii.gz",OutDir);
      OutSegFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/xtx.mat",OutDir);
      OutXtXFile = strcpyalloc(tmpstr);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--save-eres")) SaveEres=1;
    else if(!strcasecmp(option, "--save-yhat")) SaveYhat=1;
    else if(!strcasecmp(option, "--save-yhat0")) SaveYhat0=1;
    else if(!strcasecmp(option, "--replace")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&SrcReplace[nReplace]);
      sscanf(pargv[1],"%d",&TrgReplace[nReplace]);
      nReplace++;
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--merge-hypos")) {
      SrcReplace[nReplace]=78; TrgReplace[nReplace]=77; nReplace++;
      SrcReplace[nReplace]=79; TrgReplace[nReplace]=77; nReplace++;
    } 
    else if(!strcasecmp(option, "--vg-thresh")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&vg_isEqual_Threshold);
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
  printf("   --no-rescale   : do not global rescale such that mean of cerebellum WM is 100\n");
  printf("   --mask volfile : ignore areas outside of the mask (in input vol space)\n");
  printf("   --ttype+head : use tissue type def that includes head segmentation\n");
  printf("   --merge-hypos : merge left and right hypointensites into to ROI\n");
  printf("   --tt-reduce : reduce segmentation to that of a tissue type\n");
  printf("   --replace Id1 Id2 : replace seg Id1 with seg Id2\n");
  printf("\n");
  printf("   --rbv            : perform RBV PVC\n");
  printf("   --mgpvc gmthresh : perform Mueller-Gaertner PVC, gmthresh is min gm pvf bet 0 and 1\n");
  printf("   --X Xfile : save X matrix (it will be big)\n");
  printf("\n");
  printf("   --synth gtmbeta synthvolume : synthesize unsmoothed volume with gtmbeta as input\n");
  printf("\n");
  #ifdef _OPENMP
  printf("   --threads N : use N threads (with Open MP)\n");
  printf("   --threads-max : use the maximum allowable number of threads for this computer\n");
  printf("   --threads-max-1 : use one less than the maximum allowable number of threads for this computer\n");
  #endif
  printf("   --sd SUBJECTS_DIR\n");
  printf("   --vg-thresh thrshold : threshold for  'ERROR: LTAconcat(): LTAs 0 and 1 do not match'\n");
  printf("   --gdiag diagno : set diagnostic level\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
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
  std::cout << getVersion() << std::endl;
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

  return;
}
/*---------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"setenv SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"cd %s\n",cwd);
  fprintf(fp,"%s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"vgthresh   %lf\n",vg_isEqual_Threshold);
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
int VRFStats(GTM *gtm, double *vrfmean, double *vrfmin, double *vrfmax)
{
  int n,nvox,segid;
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
  }

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

  gtm->vrf = MatrixAlloc(gtm->iXtX->rows,1,MATRIX_REAL);

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
  Timer timer;
  timer.reset();

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
  LTAfree(&gtmopt->gtmsegmu2pet);

  GTMsegidlist(gtmopt->gtm);
  GTMbuildX(gtm);
  if(gtm->X==NULL) return(-1);
  GTMsolve(gtm);

  gtmopt->tLastEval = timer.seconds();
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
  gtm->segidlist = MRIsegIdListNot0(gtm->anatseg, &gtm->nsegs, 0);
  return(0);
}
/*------------------------------------------------------------------*/
int GTMnPad(GTM *gtm)
{
  double maxFWHM, maxStd;
  maxFWHM = MAX(gtm->cFWHM/gtm->yvol->xsize,
		MAX(gtm->rFWHM/gtm->yvol->ysize,gtm->sFWHM/gtm->yvol->zsize));
  maxStd = maxFWHM*sqrt(log(256.0));
  gtm->nPad = ceil(sqrt(-log(gtm->PadThresh*maxStd*sqrt(2*M_PI))*2*maxStd));
  if(Gdiag_no > 0) printf("maxFWHM = %g (voxels), PadThresh=%g, nPad=%d\n",maxFWHM,gtm->PadThresh,gtm->nPad);
  return(0);
}
/*------------------------------------------------------------------*/
int GTMsolve(GTM *gtm)
{
  Timer timer;
  int n,f;
  double sum;

  gtm->Xt = MatrixTranspose(gtm->X,gtm->Xt);
  printf("Computing  XtX ... ");fflush(stdout);
  timer.reset();
  gtm->XtX = MatrixMtM(gtm->X,gtm->XtX);
  printf(" %4.1f sec\n", timer.seconds());fflush(stdout);

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
int GTMrbv(GTM *gtm)
{
  int c,r,s,f;
  double val,v,vhat0,vhat;
  LTA *lta;
  Timer mytimer;


  if(gtm->rbv)      MRIfree(&gtm->rbv);
  if(gtm->yhat0seg) MRIfree(&gtm->yhat0seg);
  if(gtm->yhatseg)  MRIfree(&gtm->yhatseg);
  if(gtm->yseg)     MRIfree(&gtm->yseg);

  mytimer.reset() ; 

  printf("   Synthesizing in seg space... ");fflush(stdout);
  gtm->yhat0seg = GTMsegSynth(gtm);
  if(gtm->yhat0seg == NULL){
    printf("ERROR: GTMrbv() could not synthesize yhat0seg\n");
    return(1);
  }
  printf("  t = %4.2f min\n",mytimer.minutes());fflush(stdout);


  printf("   Smoothing in seg space... ");fflush(stdout);
  gtm->yhatseg = MRIgaussianSmoothNI(gtm->yhat0seg, gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  if(gtm->yhatseg == NULL){
    printf("ERROR: GTMrbv() could not smooth yhatseg\n");
    return(1);
  }
  printf("  t = %4.2f min\n",mytimer.minutes());fflush(stdout);

  printf("   Sampling input to seg space... ");fflush(stdout);
  gtm->yseg = MRIallocSequence(gtm->anatseg->width, gtm->anatseg->height, gtm->anatseg->depth,
			      MRI_FLOAT, gtm->yvol->nframes);
  if(gtm->yseg == NULL){
    printf("ERROR: GTMrbv() could not alloc yseg\n");
    return(1);
  }
  MRIcopyHeader(gtm->anatseg,gtm->yseg);
  MRIcopyPulseParameters(gtm->yvol,gtm->yseg);
  printf("  t = %4.2f min\n",mytimer.minutes());fflush(stdout);

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
  for(s=0; s < gtm->anatseg->depth; s++){
    for(c=0; c < gtm->anatseg->width; c++){
      for(r=0; r < gtm->anatseg->height; r++){
	if(MRIgetVoxVal(gtm->anatseg,c,r,s,0) < 0.5) continue;
	for(f=0; f < gtm->yvol->nframes; f++){
	  v     = MRIgetVoxVal(gtm->yseg,c,r,s,f);
	  vhat0 = MRIgetVoxVal(gtm->yhat0seg,c,r,s,f);
	  vhat  = MRIgetVoxVal(gtm->yhatseg,c,r,s,f);
	  val = v*vhat0/(vhat+FLT_EPSILON);
	  MRIsetVoxVal(gtm->rbv,c,r,s,f,val);
	}
      }
    }
  }
  MRIfree(&gtm->yseg);
  MRIfree(&gtm->yhat0seg);
  MRIfree(&gtm->yhatseg);
  printf("  t = %4.2f min\n",mytimer.minutes());fflush(stdout);

  if(gtm->mask_rbv_to_brain){
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
    segtmp = MRIreplaceList(gtm->anatseg, ReplaceThis, WithThat, nReplace, NULL, NULL);
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
  printf("  RBV took %4.2f min\n",mytimer.minutes());

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
  int c,r,s,f,n,nthseg,segid,found,nhits;
  double sum,vgmpsf,vwmpsf,vwmtac,vtac,vmgtac;
  MRI *ctxpvf, *subctxpvf, *wmpvf, *gmpvf,*gmpvfpsf,*wmpvfpsf;

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

  ctxpvf    = fMRIframe(gtm->ttpvf,0,NULL);
  subctxpvf = fMRIframe(gtm->ttpvf,1,NULL);
  wmpvf     = fMRIframe(gtm->ttpvf,2,NULL);

  gmpvf = MRIadd(ctxpvf,subctxpvf,NULL);

  gmpvfpsf = MRIgaussianSmoothNI(gmpvf,gtm->cStd, gtm->rStd, gtm->sStd, NULL);
  wmpvfpsf = MRIgaussianSmoothNI(wmpvf,gtm->cStd, gtm->rStd, gtm->sStd, NULL);

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
  Timer timer;

  timer.reset();
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
  printf("OptTimeSec %4.1f sec\n",timer.seconds());
  printf("OptTimeMin %5.2f min\n",timer.minutes());
  printf("nEvals %d\n",gtmopt->nCostEvaluations);
  printf("EvalTimeSec %4.1f sec\n", timer.seconds() / gtmopt->nCostEvaluations);
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
  newseg = MRIhiresSeg(aseg,lhw,lhp,rhw,rhp,1,&aseg2hrseg);
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
  newseg = MRIhiresSeg(aseg,lhw,lhp,rhw,rhp,1,&lta);
  MRIwrite(newseg,"newseg.annot.usf1.mgh");

  printf("computing hires seg\n");
  newseg = MRIhiresSeg(aseg,lhw,lhp,rhw,rhp,2,&lta);
  MRIwrite(newseg,"newseg.annot.usf2.mgh");

  exit(1);


}

int GTMbuildX(GTM *gtm)
{
  int nthseg;
  Timer timer;

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

  timer.reset();
  printf("computing seg pvf \n");fflush(stdout);
  gtm->segpvf = MRIseg2SegPVF(gtm->anatseg, gtm->seg2pet, 0.5, gtm->segidlist, 
			       gtm->nsegs, gtm->mask, 0, NULL, gtm->segpvf);
  if(gtm->segpvf==NULL) return(1);


  #ifdef _OPENMP
  #pragma omp parallel for if_ROMP(experimental) 
  #endif
  for(nthseg = 0; nthseg < gtm->nsegs; nthseg++){
    int segid,k,c,r,s;
    MRI *nthsegpvf=NULL,*nthsegpvfbb=NULL,*nthsegpvfbbsm=NULL;
    MRI_REGION *region;
    segid = gtm->segidlist[nthseg];
    nthsegpvf = fMRIframe(gtm->segpvf,nthseg,NULL);
    region      = REGIONgetBoundingBox(nthsegpvf,gtm->nPad);
    nthsegpvfbb = MRIextractRegion(nthsegpvf, NULL, region) ;
    nthsegpvfbbsm = MRIgaussianSmoothNI(nthsegpvfbb, gtm->cStd, gtm->rStd, gtm->sStd, nthsegpvfbbsm);
    // Fill X, creating X in this order makes it consistent with matlab
    // Note: y must be ordered in the same way.
    k = 0;
    for(s=0; s < gtm->yvol->depth; s++){
      for(c=0; c < gtm->yvol->width; c++){
	for(r=0; r < gtm->yvol->height; r++){
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
  printf(" Build time %6.4f\n", timer.seconds());fflush(stdout);

  return(0);

}

/*--------------------------------------------------------------------------*/
MRI *GTMsegSynth(GTM *gtm)
{
  int c,r,s,f,segid,segno;
  MRI *synth;

  synth = MRIallocSequence(gtm->anatseg->width,gtm->anatseg->height,gtm->anatseg->depth,MRI_FLOAT,gtm->beta->cols);
  MRIcopyHeader(gtm->anatseg,synth);
  MRIcopyPulseParameters(gtm->yvol,synth);

  for(c=0; c < gtm->anatseg->width; c++){
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
int GTMmgRefIds(GTM *gtm)
{
  int m;
  if(strcmp(gtm->mg_ref_schema,"basic")==0){
    m=0;
    gtm->mg_refids[m] =  2; m++;
    gtm->mg_refids[m] = 41; m++;
    gtm->n_mg_refids = m;
    return(gtm->n_mg_refids);
  }
  if(strcmp(gtm->mg_ref_schema,"lobes")==0){
    m=0;
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
    return(gtm->n_mg_refids);
  }

  if(strcmp(gtm->mg_ref_schema,"basic-or-lobes")==0){
    m=0;
    gtm->mg_refids[m] =  2; m++;
    gtm->mg_refids[m] = 41; m++;
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
    return(gtm->n_mg_refids);
  }

  printf("MG Reference ID schema %s unrecognized\n",gtm->mg_ref_schema);
  return(0);
}


/*--------------------------------------------------------------------------*/
int GTMprintMGRefTAC(GTM *gtm, FILE *fp)
{
  int f;
  for(f=0; f < gtm->yvol->nframes; f++)
    fprintf(fp,"%3d %10.5f\n",f,gtm->mg_reftac->rptr[f+1][1]);

  return(0);
}

/*--------------------------------------------------------------------------*/
int GTMwriteMGRefTAC(GTM *gtm, char *filename)
{
  FILE *fp;
  fp = fopen(filename,"w");
  GTMprintMGRefTAC(gtm, fp);
  fclose(fp);
  return(0);
}

/*--------------------------------------------------------------------------*/
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
  printf("gtm multiplicative scale %g (nhits=%2d) \n",gtm->scale,nhits);

  MRImultiplyConst(gtm->yvol,gtm->scale,gtm->yvol);
  MatrixScalarMul(gtm->beta, gtm->scale,gtm->beta);
  MatrixScalarMul(gtm->y,    gtm->scale,gtm->y);

  return(0);
}
