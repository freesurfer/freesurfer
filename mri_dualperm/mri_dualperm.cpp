/**
 * @brief dual permutation test and other spatial correlation functions
 *
 *
 */
/*
 * Original Author: Douglas N. Greve
 *
 * Copyright © 2025 The General Hospital Corporation (Boston, MA) "MGH"
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
/* 
To do: 
1. Add --X arg
2. Ability to test t or z or sig
3. Add a continuous default design
4. Add to log file
5. Rand - make sure no overlap between modes

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

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
#include "region.h"
#include "resample.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "randomfields.h"
#include "fsgdf.h"
#include "fmriutils.h"
#include "pdf.h"


#ifdef _OPENMP
#include "romp_support.h"
#endif

class DualPerm {
public:
  MRIGLM *mode[2]={NULL,NULL};
  int permtype[2]={0,0}; // 0=none, 1=sign, 2=shuffle, 3=signshuff, pstack
  MRI *pstack[2]={NULL,NULL};
  int pstacksave[2]={0,0};
  int perm1=1, perm2=1, perm12=1;
  int nperm = 0;
  unsigned long int seed = -1;
  int debug = 0;
  FILE *logfp=NULL;
  int mergedirs(const char *outdir, const char *srcdirs[], int nscrdirs);
};

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;
std::vector<int> ReadFrameList(char *fname, int ht, int nmax);
std::vector<int> randperm2(int ntot, int nlist, unsigned long int seed=0);
MRI *GetSubSet(MRI *mri, int nsubset, char *subsettype, char *framelistfile, unsigned long int seed,char *outfile=NULL);
std::vector<std::vector<int>> GetPStackNos(int nperm, int n1, int n2, int seed);

class ModeArg {
public:
  char *modefile=NULL;
  char *maskfile=NULL;
  char *pstackfile=NULL;
  int pstacksave=0;
  char *fsgdfile=NULL;
  int osgm=0;
  int tsgd=0;
  int nsubset=0;
  char *subsettype=NULL;
  int prune = 1;
  char *framelistfile=NULL;
  int residualize=1;
  std::vector<int> framelist;
};

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;
char *outdir=NULL;
int threads = 1;
char *SUBJECTS_DIR=NULL;
int modeno = 0;
unsigned long int seed = -1;
int nperm = 0;
ModeArg marg[2];
DualPerm dp;
float prune_thr = FLT_MIN;
MATRIX *PermutationSquences(int nrows, int nperm, int ptype, unsigned long int seed);
MATRIX *ApplyPermutation(MATRIX *X, int permno, MATRIX *shuffle, MATRIX *fsign);
double GetPVal(double val, std::vector<double> vallist, int psign);
int writepvals(const char *outdir, double cc0, std::vector<double> cc1,std::vector<double> cc2,std::vector<double> cc12);
int nmodes = 2;
int SavePX = 0;
int SaveInput = 0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err;
  char logfile[1000];
  char fname[1000];
  char glmdir[1000],cdir[1000];
  FILE *fp=NULL;

  nargs = handleVersionOption(argc, argv, "mri_gtmpvc");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

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

  Timer timer, mytimer;

  printf("Creating output directory %s\n",outdir);
  err = mkdir(outdir,0777);
  if(err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",outdir);
    perror(NULL);
    return(1);
  }
  sprintf(logfile,"%s/mri_dualperm.log",outdir);
  FILE *logfp = fopen(logfile,"w");
  dump_options(logfp);
  sprintf(fname,"%s/seed.txt",outdir);
  fp = fopen(fname,"w");
  fprintf(fp,"%lu\n",seed);
  fclose(fp);

  for(int modeno=0; modeno < nmodes; modeno++){
    ModeArg *ma = &marg[modeno];
    dp.mode[modeno] = (MRIGLM *) calloc(sizeof(MRIGLM),1);
    MRIGLM *dpm = dp.mode[modeno];

    dpm->y = MRIread(ma->modefile);
    if(dpm->y==NULL) exit(1);

    if(ma->maskfile){
      dpm->mask = MRIread(ma->maskfile);
      if(dpm->mask==NULL) exit(1);
    }
    // Prune mask based on all data (before subtype) so that mask is consistent
    if(ma->prune) dpm->mask = MRIframeBinarize(dpm->y,prune_thr,dpm->mask);
    if(dpm->mask == NULL) exit(1);
    sprintf(fname,"%s/residualize%d.dat",outdir,modeno+1);
    fp = fopen(fname,"w");
    fprintf(fp,"%d\n",ma->residualize);
    fclose(fp);
  }
  // Create a single mask using masks from both modes
  if(nmodes == 2){
    dp.mode[0]->mask = MRIframeBinarize(dp.mode[1]->mask,prune_thr,dp.mode[0]->mask);
    dp.mode[1]->mask = MRIframeBinarize(dp.mode[0]->mask,prune_thr,dp.mode[1]->mask);
  }

  MRIGLM *pmode[2];
  MATRIX *fsign[2]={NULL,NULL};
  MATRIX *shuffle[2]={NULL,NULL};
  MRI *modemap[2]={NULL,NULL};
  for(int modeno=0; modeno < nmodes; modeno++){
    ModeArg *ma = &marg[modeno];
    MRIGLM *dpm = dp.mode[modeno];
    if(ma->pstackfile){
      if(dpm->y->nframes != 1){
	printf("ERROR: mode%d nframes=%d, must be 1 when using pstack\n",modeno,dpm->y->nframes);
	exit(1);
      }
      dp.pstack[modeno] = MRIread(ma->pstackfile);
      if(dp.pstack[modeno]==NULL) exit(1);
      modemap[modeno] = dpm->y;
      continue; // nothing more to do when pstack is passed
    }

    if(ma->subsettype){
      sprintf(fname,"%s/framelist%d.txt",outdir,modeno+1);
      MRI *mritmp = GetSubSet(dpm->y, ma->nsubset, ma->subsettype, ma->framelistfile, seed+modeno+1,fname);
      if(!mritmp) exit(1);
      dpm->y = mritmp;
    }

    dpm->glm = GLMalloc();
    if(ma->osgm){
      dpm->Xg = MatrixConstVal(1.0,dpm->y->nframes,1,NULL);
      dpm->glm->ncontrasts = 1;
      dpm->glm->Cname[0] = strcpyalloc("osgm");
      dpm->glm->C[0] = MatrixConstVal(1.0,1,1,NULL);
    }
    if(ma->tsgd){
      dpm->Xg = MatrixAlloc(dpm->y->nframes,2,MATRIX_REAL);
      for(int n=0; n < dpm->y->nframes; n++){
	if(n < dpm->y->nframes/2) dpm->Xg->rptr[n+1][1] = 1;
	else                      dpm->Xg->rptr[n+1][2] = 1;
      }
      dpm->glm->ncontrasts = 1;
      dpm->glm->Cname[0] = strcpyalloc("tsgd");
      dpm->glm->C[0] = MatrixAlloc(1,2,MATRIX_REAL);
      dpm->glm->C[0]->rptr[1][1] = +1;
      dpm->glm->C[0]->rptr[1][2] = -1;
    }
    if(ma->fsgdfile){
      char gd2mtx_method[1000];
      sprintf(gd2mtx_method,"none");
      FSGD *fsgd = gdfRead(ma->fsgdfile,gd2mtx_method,0);
      if(fsgd==NULL) exit(1);
      dpm->Xg = gdfMatrix(fsgd,fsgd->gd2mtx_method,NULL);
      if(dpm->Xg==NULL) exit(1);
      if(dpm->Xg->rows != dpm->y->nframes){
	printf("ERROR: fsgd %s has %d rows, expecting %d\n",ma->fsgdfile,dpm->Xg->rows,dpm->y->nframes);
	exit(1);
      }
      if(fsgd->nContrasts == 0){
	printf("ERROR: fsgd %s has no contrasts\n",ma->fsgdfile);
	exit(1);
      }
      if(fsgd->nContrasts > 1){
	printf("WARNING: fsgd %s has %d contrasts, only using first one\n",ma->fsgdfile,fsgd->nContrasts);
	exit(1);
      }
      dpm->glm->C[0] = MatrixCopy(fsgd->C[0],NULL);
      dpm->glm->Cname[0] = strcpyalloc(fsgd->ContrastName[0]);
      dpm->glm->ncontrasts = 1;
    }
    if(dpm->glm->C[0]->rows != 1){
      printf("ERROR: contrast has %d rows, must be 1\n",dpm->glm->C[0]->rows);
      exit(1);
    }
    if(dpm->Xg->cols == 1 && (dp.permtype[modeno] == 2 || dp.permtype[modeno] == 3)){
      printf("ERROR: design matrix has only 1 column, but shuffling requested\n");
      exit(1);
    }

    double Xcond = MatrixNSConditionNumber(dpm->Xg);
    printf("Matrix condition is %g\n",Xcond);
    sprintf(glmdir,"%s/glm%d",outdir,modeno+1);
    sprintf(cdir,"%s/%s",glmdir,dpm->glm->Cname[0]);
    err = mkdir(glmdir,0777);
    if(err != 0 && errno != EEXIST) exit(err);
    err = mkdir(cdir,0777);
    if(err != 0 && errno != EEXIST) exit(err);
    sprintf(fname,"%s/Xg.dat",glmdir);
    MatrixWriteTxt(fname, dpm->Xg);
    sprintf(fname,"%s/C.dat",glmdir);
    MatrixWriteTxt(fname, dpm->glm->C[0]);
    sprintf(fname,"%s/mask.nii.gz",glmdir);
    MRIwrite(dpm->mask,fname);
    MRIglmFitAndTest(dpm);
    sprintf(fname,"%s/beta.nii.gz",glmdir);
    MRIwrite(dpm->beta,fname);
    sprintf(fname,"%s/rvar.nii.gz",glmdir);
    MRIwrite(dpm->rvar,fname);
    sprintf(fname,"%s/gamma.nii.gz",cdir);
    MRIwrite(dpm->gamma[0],fname);
    if(dpm->mask) MRImask(dpm->sig[0],dpm->mask,dpm->sig[0],0.0,0.0);
    MRIsetSign(dpm->sig[0],dpm->gamma[0],0);
    sprintf(fname,"%s/sig.nii.gz",cdir);
    MRIwrite(dpm->sig[0],fname);
    int nmask = MRIcountNonzero(dpm->mask);
    sprintf(fname,"%s/nmask.dat",glmdir);
    fp = fopen(fname,"w");
    fprintf(fp,"%d\n",nmask);
    fclose(fp);
    printf("nmask %d\n",nmask);

    modemap[modeno] = dpm->gamma[0];
    pmode[modeno] = MRIglmDeepCopy(dpm);
    if(ma->residualize) pmode[modeno]->y = MRIcopy(dpm->eres,NULL);
    if(SaveInput){
      // Do this after any residualization
      char inputoutputname[1000];
      sprintf(inputoutputname,"%s/stack%d.nii.gz",outdir,modeno+1);
      MRIwrite(dpm->y,inputoutputname);
      sprintf(inputoutputname,"%s/res%d.nii.gz",outdir,modeno+1);
      MRIwrite(dpm->eres,inputoutputname);
      printf("%d resid %d\n",modeno,ma->residualize);
    }

    if(nperm > 0) {
      if(dp.permtype[modeno] == 1 || dp.permtype[modeno] == 3){
	fsign[modeno] = PermutationSquences(dpm->y->nframes, nperm, 2, seed+modeno+1);
	sprintf(fname,"%s/fsign%d.mat",outdir,modeno+1);
	MatrixWrite(fsign[modeno], fname, "fsign");
      }
      if(dp.permtype[modeno] == 2 || dp.permtype[modeno] == 3){
	shuffle[modeno] = PermutationSquences(dpm->y->nframes, nperm, 1, seed+modeno+1);
	sprintf(fname,"%s/shuffle%d.mat",outdir,modeno+1);
	MatrixWrite(shuffle[modeno], fname, "shuffle");
      }
    }
  } // mode

  if(marg[0].pstackfile || marg[1].pstackfile){
    int n1=0, n2=0;
    if(marg[0].pstackfile) n1 = dp.pstack[0]->nframes;
    if(marg[1].pstackfile) n2 = dp.pstack[1]->nframes;
    std::vector<std::vector<int>> psno = GetPStackNos(nperm, n1, n2, seed);
    sprintf(fname,"%s/pstack.dat",outdir);
    FILE *fpstack = fopen(fname,"w");
    for(int n=0; n < psno.size(); n++){
      fprintf(fpstack,"%d %d\n",psno[n][0],psno[n][1]);
    }
    fclose(fpstack);
  }

  // Write out the basic data and set up a file pointer for loop
  FILE *fp1=NULL, *fp2=NULL, *fp12=NULL;
  std::vector<double> cc0;
  double cc0val = 0;
  if(nmodes == 2){
    SpatialCor sc(modemap[0],modemap[1],dp.mode[0]->mask);
    sc.pearsoncor();
    sprintf(fname,"%s/cc.dat",outdir);
    cc0.push_back(sc.pcc);
    printf("cc0 = %12.9lf\n",sc.pcc);
    cc0val = sc.pcc;
    fp = fopen(fname,"w");
    fprintf(fp,"%12.9lf\n",sc.pcc);
    fclose(fp);
    sprintf(fname,"%s/cc.info.dat",outdir);
    sc.print(fname);
    if(dp.permtype[0] != 0){
      sprintf(fname,"%s/cc1.perm.dat",outdir);
      fp1 = fopen(fname,"w");
    }
    if(dp.permtype[1] != 0){
      sprintf(fname,"%s/cc2.perm.dat",outdir);
      fp2 = fopen(fname,"w");
    }
    if(dp.permtype[0] != 0 && dp.permtype[1] != 0){
      sprintf(fname,"%s/cc12.perm.dat",outdir);
      fp12 = fopen(fname,"w");
    }
  }

  // permutation loop
  MRI *pmodemap[2] = {NULL,NULL};
  MRI *pstacksave[2] = {NULL,NULL};
  std::vector<double> cc, cc1, cc2, cc12;
  for(int n=0; n < nperm; n++){
    for(int modeno=0; modeno < nmodes; modeno++){
      if(dp.permtype[modeno] == 0) continue;
      MRIGLM *dpm = dp.mode[modeno];
      MRIGLM *dppm = pmode[modeno];
      ModeArg *ma  = &marg[modeno];
      if(marg[modeno].pstackfile==NULL){
	if(dppm->Xg) MatrixFree(&dppm->Xg);
	dppm->Xg = ApplyPermutation(dpm->Xg, n, shuffle[modeno], fsign[modeno]);
	if(SavePX){
	  char pxfile[1000];
	  sprintf(pxfile,"%s/pX.m%d.%04d.mtx",outdir,modeno+1,n);
	  MatrixWriteTxt(pxfile, dppm->Xg);
	}
	MRIglmFitAndTest(dppm,1);
	// readd: add gamma from dpm to the gamma of dppm
	if(ma->residualize == 2) MRIadd(dpm->gamma[0],dppm->gamma[0],dppm->gamma[0]);
	if(ma->pstacksave){
	  if(pstacksave[modeno]==NULL) {
	    pstacksave[modeno] = MRIcloneBySpace(dpm->mask, MRI_FLOAT, nperm);
	    if(!pstacksave[modeno]) exit(1);
	  }
	  fMRIinsertFrame(dppm->gamma[0], 0, pstacksave[modeno], n);
	} 
	pmodemap[modeno] = dppm->gamma[0] ;
      }
      else {
	pmodemap[modeno] = fMRIframe(dp.pstack[modeno], n, pmodemap[modeno]);
      }
    }

    printf("%4d ",n);
    if(nmodes == 1) {
      printf("\n"); fflush(stdout);
      continue;
    }
    if(dp.permtype[0] != 0){
      SpatialCor sc(pmodemap[0],modemap[1],dp.mode[0]->mask);
      sc.pearsoncor();
      cc1.push_back(sc.pcc);
      fprintf(fp1,"%21.18lf\n",sc.pcc);
      printf("%12.9lf ",sc.pcc);fflush(fp1);
      sprintf(fname,"%s/cc1.info.dat",outdir);
      sc.printline(fname);
      // add opposite sign?
    }
    if(dp.permtype[1] != 0){
      SpatialCor sc(modemap[0],pmodemap[1],dp.mode[0]->mask);
      sc.pearsoncor();
      cc2.push_back(sc.pcc);
      fprintf(fp2,"%21.18lf\n",sc.pcc);
      printf("%12.9lf ",sc.pcc);fflush(fp2);
      sprintf(fname,"%s/cc2.info.dat",outdir);
      sc.printline(fname);
      // add opposite sign?
    }
    if(dp.permtype[0] != 0 && dp.permtype[1] != 0){
      SpatialCor sc(pmodemap[0],pmodemap[1],dp.mode[0]->mask);
      sc.pearsoncor();
      cc12.push_back(sc.pcc);
      fprintf(fp12,"%21.18lf\n",sc.pcc); fflush(fp12);
      printf("%12.9lf ",sc.pcc);
      sprintf(fname,"%s/cc12.info.dat",outdir);
      sc.printline(fname);
      // add opposite sign?
    }
    printf("\n"); fflush(stdout);
  }

  if(nmodes == 2){
    if(fp1) fclose(fp1);
    if(fp2) fclose(fp2);
    if(fp12) fclose(fp12);

    double p1pos=0, p1neg=0, p1abs=0, p2pos=0, p2neg=0, p2abs=0, p12pos=0, p12neg=0, p12abs=0;
    double cc1std=0,cc2std=0,cc12std=0;
    if(dp.permtype[0] != 0){
      sprintf(fname,"%s/p1.dat",outdir);
      fp = fopen(fname,"w");
      p1pos = GetPVal(cc0[0], cc1, +1); fprintf(fp,"%31.28lf ",p1pos);
      p1neg = GetPVal(cc0[0], cc1, -1); fprintf(fp,"%31.28lf ",p1neg);
      p1abs = GetPVal(cc0[0], cc1,  0); fprintf(fp,"%31.28lf ",p1abs);
      fprintf(fp,"\n");
      fclose(fp);
      cc1std = stddevvector(cc1);
      sprintf(fname,"%s/cc1.std.dat",outdir);
      fp = fopen(fname,"w");
      fprintf(fp,"%31.28lf\n",cc1std);
      fclose(fp);
    }
    if(dp.permtype[1] != 0){
      sprintf(fname,"%s/p2.dat",outdir);
      fp = fopen(fname,"w");
      p2pos = GetPVal(cc0[0], cc2, +1); fprintf(fp,"%31.28lf ",p2pos);
      p2neg = GetPVal(cc0[0], cc2, -1); fprintf(fp,"%31.28lf ",p2neg);
      p2abs = GetPVal(cc0[0], cc2,  0); fprintf(fp,"%31.28lf ",p2abs);
      fprintf(fp,"\n");
      fclose(fp);
      cc2std = stddevvector(cc2);
      sprintf(fname,"%s/cc2.std.dat",outdir);
      fp = fopen(fname,"w");
      fprintf(fp,"%31.28lf\n",cc2std);
      fclose(fp);
    }
    if(dp.permtype[0] != 0 && dp.permtype[1] != 0){
      sprintf(fname,"%s/p12.dat",outdir);
      fp = fopen(fname,"w");
      p12pos = GetPVal(cc0[0], cc12, +1); fprintf(fp,"%31.28lf ",p12pos);
      p12neg = GetPVal(cc0[0], cc12, -1); fprintf(fp,"%31.28lf ",p12neg);
      p12abs = GetPVal(cc0[0], cc12,  0); fprintf(fp,"%31.28lf ",p12abs);
      fprintf(fp,"\n");
      fclose(fp);

      sprintf(fname,"%s/p12.max.dat",outdir);
      fp = fopen(fname,"w");
      if(p1pos > p2pos) fprintf(fp,"%31.28lf ",p1pos);
      else              fprintf(fp,"%31.28lf ",p2pos);
      if(p1neg > p2neg) fprintf(fp,"%31.28lf ",p1neg);
      else              fprintf(fp,"%31.28lf ",p2neg);
      if(p1abs > p2abs) fprintf(fp,"%31.28lf ",p1abs);
      else              fprintf(fp,"%31.28lf ",p2abs);
      fprintf(fp,"\n");
      fclose(fp);

      cc12std = stddevvector(cc12);
      sprintf(fname,"%s/cc12.std.dat",outdir);
      fp = fopen(fname,"w");
      fprintf(fp,"%31.28lf\n",cc12std);
      fclose(fp);

      sprintf(fname,"%s/p.maxdv.dat",outdir);
      fp = fopen(fname,"w");
      int pick = 0;
      if(fabs(cc12std-cc1std) > fabs(cc12std-cc2std)) {
	pick = 1;
	fprintf(fp,"%31.28lf %31.28lf %31.28lf\n",p1pos,p1neg,p1abs);
	printf("MaxDV %31.28lf %31.28lf %31.28lf\n",p1pos,p1neg,p1abs);
      }
      else {
	pick = 2;
	fprintf(fp,"%31.28lf %31.28lf %31.28lf\n",p2pos,p2neg,p2abs);
	printf("MaxDV %31.28lf %31.28lf %31.28lf\n",p2pos,p2neg,p2abs);
      }
      fclose(fp);

      sprintf(fname,"%s/cc.perm.std.dat",outdir);
      fp = fopen(fname,"w");
      fprintf(fp,"%31.28lf %31.28lf %31.28lf %d\n",cc1std,cc2std,cc12std,pick);
      fclose(fp);
      
      printf("cc0 %31.28lf\n",cc0val);
      printf("pick %d\n",pick);
      printf("CCSTD: %31.28lf %31.28lf %31.28lf %d\n",cc1std,cc2std,cc12std,pick);
      
    }
  }

  for(int modeno=0; modeno < nmodes; modeno++){
    ModeArg *ma  = &marg[modeno];
    if(ma->pstacksave){
      sprintf(fname,"%s/surrogates%d.nii.gz",outdir,modeno+1);
      printf("Writing mode %d surrogates to %s\n",modeno+1,fname);
      err = MRIwrite(pstacksave[modeno],fname);
      if(err) exit(1);
    }
  }

  fprintf(logfp,"#VMPC# mris_dualperm VmPeak  %d\n",GetVmPeak());
  fprintf(logfp,"mris_dualperm-runtime %5.2f min\n",timer.minutes());
  fprintf(logfp,"mris_dualperm done\n");
  fclose(logfp);
  printf("#VMPC# mris_dualperm VmPeak  %d\n",GetVmPeak());
  printf("mris_dualperm-runtime %5.2f min\n",timer.minutes());
  printf("mris_dualperm done\n");
  return(0);
  exit(0);

} // end of main

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
    else if(!strcasecmp(option, "--mode1-only")) nmodes=1;
    else if(!strcasecmp(option, "--map1-only")) nmodes=1;
    else if(!strcasecmp(option, "--o")) {
      if(nargc < 1) CMDargNErr(option,1);
      outdir = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--mode1") || !strcasecmp(option, "--mode2") ||
	    !strcasecmp(option, "--map1") || !strcasecmp(option, "--map2")) {
      if(!strcasecmp(option, "--mode1") || !strcasecmp(option, "--map1")) modeno = 1;
      if(!strcasecmp(option, "--mode2") || !strcasecmp(option, "--map2")) modeno = 2;
      if(nargc < 1) CMDargNErr(option,1);
      marg[modeno-1].modefile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--mask1") || !strcasecmp(option, "--mask2")) {
      if(!strcasecmp(option, "--mask1")) modeno = 1;
      if(!strcasecmp(option, "--mask2")) modeno = 2;
      if(nargc < 1) CMDargNErr(option,1);
      marg[modeno-1].maskfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--fsgd1") || !strcasecmp(option, "--fsgd2")) {
      if(!strcasecmp(option, "--fsgd1")) modeno = 1;
      if(!strcasecmp(option, "--fsgd2")) modeno = 2;
      if(nargc < 1) CMDargNErr(option,1);
      marg[modeno-1].fsgdfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--osgm1") || !strcasecmp(option, "--osgm2")) {
      if(!strcasecmp(option, "--osgm1")) modeno = 1;
      if(!strcasecmp(option, "--osgm2")) modeno = 2;
      marg[modeno-1].osgm = 1;
    }
    else if(!strcasecmp(option, "--tsgd1") || !strcasecmp(option, "--tsgd2")) {
      if(!strcasecmp(option, "--tsgd1")) modeno = 1;
      if(!strcasecmp(option, "--tsgd2")) modeno = 2;
      marg[modeno-1].tsgd = 1;
    }
    else if(!strcasecmp(option, "--subset1") || !strcasecmp(option, "--subset2")) {
      if(!strcasecmp(option, "--subset1")) modeno = 1;
      if(!strcasecmp(option, "--subset2")) modeno = 2;
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&marg[modeno-1].nsubset);
      marg[modeno-1].subsettype = pargv[1];
      nargsused = 2;
    }
    else if(!strcasecmp(option, "--no-residualize1") || !strcasecmp(option, "--no-residualize2")) {
      if(!strcasecmp(option, "--no-residualize1")) modeno = 1;
      if(!strcasecmp(option, "--no-residualize2")) modeno = 2;
      marg[modeno-1].residualize = 0;
    }
    else if(!strcasecmp(option, "--readd1") || !strcasecmp(option, "--readd2")) {
      if(!strcasecmp(option, "--readd1")) modeno = 1;
      if(!strcasecmp(option, "--readd2")) modeno = 2;
      marg[modeno-1].residualize = 2;
    }
    else if(!strcasecmp(option, "--surrogates1") || !strcasecmp(option, "--surrogates2")) {
      if(!strcasecmp(option, "--surrogates1")) modeno = 1;
      if(!strcasecmp(option, "--surrogates2")) modeno = 2;
      if(nargc < 1) CMDargNErr(option,1);
      marg[modeno-1].pstackfile = pargv[0];
      dp.permtype[modeno-1] = 4;
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--save-surrogates1") || !strcasecmp(option, "--save-surrogates2")) {
      if(!strcasecmp(option, "--save-surrogates1")) modeno = 1;
      if(!strcasecmp(option, "--save-surrogates2")) modeno = 2;
      marg[modeno-1].pstacksave = 1;
    }
    else if(!strcasecmp(option, "--list1") || !strcasecmp(option, "--list2")) {
      if(!strcasecmp(option, "--list1")) modeno = 1;
      if(!strcasecmp(option, "--list2")) modeno = 2;
      if(nargc < 1) CMDargNErr(option,1);
      marg[modeno-1].framelistfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--ptype1") || !strcasecmp(option, "--ptype2")) {
      if(!strcasecmp(option, "--ptype1")) modeno = 1;
      if(!strcasecmp(option, "--ptype2")) modeno = 2;
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&dp.permtype[modeno-1]); //0=none,1=sign,2=shuffle,3=sign+shuffle
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--pX-save")){
      SavePX = 1;
    }
    else if(!strcasecmp(option, "--save-input") || !strcasecmp(option, "--input-save")){
      SaveInput = 1;
    }
    else if(!strcasecmp(option, "--nperm")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nperm);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--threads")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&threads);
      #ifdef _OPENMP
      omp_set_num_threads(threads);
      #endif
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--seed")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lu",&seed);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--merge")){
      //outdir srcdir1 srcdir2 ...
      if(nargc < 3) CMDargNErr(option,3);
      char const *srcdirs[50];
      int nsrcdirs=nargc-1;
      for(int n=1; n < nargc; n++) srcdirs[n-1] = pargv[n];
      int err = dp.mergedirs(pargv[0],srcdirs,nsrcdirs);
      exit(err);
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
  printf("USAGE: %s   See: surfer.nmr.mgh.harvard.edu/fswiki/DualPerm\n",Progname) ;
  printf("\n");
  printf("   --o outdir : output folder\n");
  printf("   --nperm npermutations : number of permutations\n");
  printf("   --mapN mapNfile : input stack for modality N=1,2\n");
  printf("   --maskN maskNfile : input mask for modality N=1,2\n");
  printf("   --fsgdN fsgdNfile : FreeSurfer Group Descriptor File for modality N=1,2\n"); 
  printf("     see surfer.nmr.mgh.harvard.edu/fswiki/FsgdExamples\n");
  printf("   --osgmN : analyze modality N as a one-sample-group-mean (OSGM), instead of --fsgd \n");
  printf("   --ptypeN permtype (0=noperm, 1=sign, 2=shuffle, 3=sign+shufflen\n");
  printf("   --no-residualizeN : do not residualize modality N before permuting (not recommended)\n");
  printf("   --save-surrogatesN : save all surrogates in one file (a stack) for modality N=1,2\n");
  printf("   --save-input : save input maps\n");
  printf("   --map1-only only analyze modality 1 (for creating atlases and surrogatess); modality 2 not needed\n");
  printf("   --surrogatesN surrogatesNfile : input surrogate stack for modality N\n");
  printf("     as output from another run of mri_dualperm, possibly using --map1-only.\n");
  printf("     If using this option, then do not supply --mapN, --maskN, or --fsgdN\n");
  #ifdef _OPENMP
  printf("   --threads nthreads : use nthreads threads (with Open MP, speeds up a little but not a lot)\n");
  #endif
  printf("   \n");
  printf("   These options are helpful for testing\n");
  printf("   --seed seed : random seed. If not specified, then uses time-of-day\n");
  printf("   --tsgdN : analyze modality N as a two-sample-group-diff (TSGD) \n");
  printf("   --subsetN nsubset type (first, last, rand)\n");
  printf("   --listN framelistfile\n");
  printf("   --pX-save : save permuted design matrices\n");
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
  if(outdir == NULL){
    printf("ERROR: must spec outdir\n");
    exit(1);
  }
  for(int modeno=0; modeno < nmodes; modeno++){
    ModeArg *ma = &marg[modeno];
    if(ma->modefile == NULL){
      printf("ERROR: must spec --mode%d\n",modeno+1);
      exit(1);
    }
    if(ma->pstackfile == NULL){
      if(ma->osgm + ma->tsgd + (ma->fsgdfile!=NULL) == 0){
	printf("ERROR: mode %d must spec a design with one of --osgm --tsgd or --fsgd\n",modeno+1);
	exit(1);
      }
      if(ma->osgm + ma->tsgd + (ma->fsgdfile!=NULL) > 1){
	printf("ERROR: mode %d can only spec one of --osgm --tsgd or --fsgd\n",modeno+1);
	exit(1);
      }
    }
    else{
      if(ma->osgm || ma->tsgd || ma->fsgdfile || ma->nsubset != 0 || ma->framelistfile){
	printf("ERROR: mode %d cannot spec  --pstack and --osgm or --tsgd or --fsgd or --subset or --list\n",modeno+1);
	exit(1);
      }
    }
  }
  printf("permtype %d %d\n",dp.permtype[0],dp.permtype[1]);
  if(dp.permtype[0] == 0 && dp.permtype[1] == 0 && nperm > 0){
    printf("ERROR: nperm > 0 but neither mode is being permuted\n");
    exit(1);
  }
  if(seed == 0) seed = PDFtodSeed();

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
  fprintf(fp,"nperm   %d\n",nperm);
  fprintf(fp,"seed   %lu\n",seed);
  return;
}

std::vector<int> ReadFrameList(char *fname)
{
  std::vector<int> framelist;
  FILE *fp = fopen(fname,"r");
  if(!fp){
    printf("ERROR: ReadFrameList: could not open %s\n",fname);
    return(framelist);
  }
  char *line;
  size_t len = 1000;
  line = (char*)calloc(sizeof(char),len);
  while(1){
    size_t nread = getline(&line,&len,fp);
    if(nread == -1) break;
    int f;
    // read the first item on the line and ignnore everything else
    sscanf(line,"%d",&f);
    framelist.push_back(f);
  }
  fclose(fp);
  free(line);
  
  printf("nlist %d\n",(int)framelist.size());
  return(framelist);
}


std::vector<int> randperm2(int ntot, int nlist, unsigned long int seed)
{
  std::vector<int> rp;
  if(nlist > ntot) {
    printf("ERROR: randperm2(): ntot=%d < nlist=%d\n",ntot,nlist);
    return(rp);
  }

  RFS *rfs = RFspecInit(0, NULL);// could set seed here too
  rfs->name = strcpyalloc("uniform");
  rfs->params[0] = 0;
  rfs->params[1] = 1;
  RFspecSetSeed(rfs, seed);

  std::vector<int> v;
  for(int n=0; n < ntot; n++) v.push_back(n);
  for(int n=0; n < ntot; n++){
    int n2 = (int)floor(RFdrawVal(rfs)*ntot);
    int tmp = v[n];
    v[n] = v[n2];
    v[n2] = tmp;
  }

  for(int n=0; n < nlist; n++) rp.push_back(v[n]);
  //for(int n=0; n < nlist; n++) printf("%3d %3d\n",n,v[n]);
  //fflush(stdout);

  RFspecFree(&rfs);
  return(rp);
}

MRI *GetSubSet(MRI *mri, int nsubset, char *subsettype, char *framelistfile, unsigned long int seed, char *outfile)
{
  if(mri->nframes < nsubset){
    printf("ERROR: nframes = %d < nsubset = %d\n",mri->nframes,nsubset);
    return(NULL);
  }
  std::vector<int> framelist0;
  if(framelistfile) {
    framelist0 = ReadFrameList(framelistfile);
    if(framelist0.size()==0) return(NULL);
  }
  else for(int n=0; n < mri->nframes; n++) framelist0.push_back(n); 

  std::vector<int> framelist;
  if(strcmp(subsettype,"first")==0) for(int n=0; n < nsubset; n++) framelist.push_back(framelist0[n]);
  if(strcmp(subsettype,"last")==0)  for(int n=0; n < nsubset; n++) framelist.push_back(framelist0[mri->nframes-n-1]);
  if(strcmp(subsettype,"rand")==0) framelist = randperm2(mri->nframes, nsubset, seed);

  if(outfile){
    FILE *fp = fopen(outfile,"w");
    for(int f=0; f < nsubset; f++) fprintf(fp,"%4d\n",framelist[f]);
    fclose(fp);
  }

  MRI *out = MRIcloneBySpace(mri,mri->type,nsubset);
  if(!out) return(NULL);

  MRIcopyPulseParameters(mri,out);

  for(int c=0; c < mri->width; c++){
    for(int r=0; r < mri->height; r++){
      for(int s=0; s < mri->depth; s++){
	for(int f=0; f < nsubset; f++){
	  double val = MRIgetVoxVal(mri,c,r,s,framelist[f]);
	  MRIsetVoxVal(out,c,r,s,f,val);
	}
      }
    }
  }

  return(out);
}

MATRIX *PermutationSquences(int nrows, int nperm, int ptype, unsigned long int seed)
{
  // Note: seed = 0 means to pick randomly
  MATRIX *ps = MatrixAlloc(nrows,nperm,MATRIX_REAL);
  if(ptype == 1){ // shuffle
    for(int n=0; n < nperm; n++){
      std::vector<int> rp = randperm2(nrows,nrows,seed+n);
      for(int r=0; r < nrows; r++) ps->rptr[r+1][n+1] = rp[r];
    }
  }
  else{  // ptype == 2
    RFS *rfs = RFspecInit(0, NULL);
    rfs->name = strcpyalloc("uniform");
    rfs->params[0] = 0;
    rfs->params[1] = 1;
    RFspecSetSeed(rfs, seed);
    for(int n=0; n < nperm; n++){
      for(int r=0; r < nrows; r++) {
	double v = RFdrawVal(rfs);
	if(v < 0.5) ps->rptr[r+1][n+1] = -1;
	else        ps->rptr[r+1][n+1] = +1;
      }
    }
    RFspecFree(&rfs);
  }
  return(ps);
}

MATRIX *ApplyPermutation(MATRIX *X, int permno, MATRIX *shuffle, MATRIX *fsign)
{
  if(!shuffle && !fsign){
    printf("ERROR: ApplyPermutation(): neither shuffle nor fsign are non-null\n");
    return(NULL);
  }
  if(shuffle && fsign){
    if(shuffle->rows != fsign->rows || shuffle->cols != fsign->cols){
      printf("ERROR: ApplyPermutation(): dim mismatch between shuffle and fsign %d %d  %d %d\n",
	     shuffle->rows,fsign->rows,shuffle->cols,fsign->cols);
      return(NULL);
    }
  }
  int nrows = 0, ncols = 0;
  if(shuffle) {nrows = shuffle->rows; ncols = shuffle->cols;}
  if(fsign) {nrows = fsign->rows; ncols = fsign->cols;}
  if(permno >= ncols){
    printf("ERROR: ApplyPermutation(): permno =%d >= ncols = %d\n",permno,ncols);
    return(NULL);
  }
  if(X->rows != nrows){
    printf("ERROR: ApplyPermutation(): X has %d rows, but shuffle/fsign has  %d\n",X->rows,nrows);
    return(NULL);
  }

  MATRIX *Xp = MatrixAlloc(X->rows,X->cols,MATRIX_REAL);
  for(int r=0; r < nrows; r++){
    int rp = r;
    if(shuffle) rp = shuffle->rptr[r+1][permno+1];
    double s = 1;
    if(fsign) s = fsign->rptr[r+1][permno+1];
    for(int c=0; c < X->cols; c++) {
      Xp->rptr[r+1][c+1] = s*X->rptr[rp+1][c+1];
    }
  }
  return(Xp);
}


double GetPVal(double val, std::vector<double> vallist, int psign)
{
  int count=0;
  for(int n=0; n < vallist.size(); n++){
    if(psign > 0  && vallist[n] >= val) count++;
    if(psign < 0  && vallist[n] <= val) count++;
    if(psign == 0 && fabs(vallist[n]) >= fabs(val)) count++;
  } 
  double p = double(count)/vallist.size();
  return(p);
}

// When a pstack is used, this gets the frames to use at each permutation step
std::vector<std::vector<int>> GetPStackNos(int nperm, int n1, int n2, int seed)
{
  // n? is the number of frames in pstack?
  if(n2 == 0 && n1 < nperm){
    printf("ERROR: GetPStackNos(): n2=0 and n1=%d < nperm=%d\n",n1,nperm);
    exit(1);// should probably return empty vector
  }
  if(n1 == 0 && n2 < nperm){
    printf("ERROR: GetPStackNos(): n1=0 n2=%d < nperm=%d\n",n2,nperm);
    exit(1);
  }
  if(n1>0 && n2>0 && n1*n2 < nperm){
    printf("ERROR: GetPStackNos(): n1=%d n2=%d n1*n2=%d < nperm=%d\n",n1,n2,n1*n2,nperm);
    exit(1);
  }
  std::vector<std::vector<int>> psno;
  if(n2 == 0){ // only n1, no repeats, just use 1 to nperm
    for(int n=0; n < nperm; n++){
      std::vector<int> nn = {n,0};
      psno.push_back(nn);
    }
    return(psno);
  }
  if(n1 == 0){ // only n2, no repeats, just use 1 to nperm
    for(int n=0; n < nperm; n++){
      std::vector<int> nn = {0,n};
      psno.push_back(nn);
    }
    return(psno);
  }
  if(n1 >= nperm && n2 >= nperm ){
    // n1 and n2, no repeats, just use 1-nperm
    for(int n=0; n < nperm; n++){
      std::vector<int> nn = {n,n};
      psno.push_back(nn);
    }
    return(psno);
  }
  // to get here n1>0 and n2>0 and either n1 or n2 > nperm so have to repeat

  RFS *rfs = RFspecInit(0, NULL);
  rfs->name = strcpyalloc("uniform");
  rfs->params[0] = 0;
  rfs->params[1] = 1;
  RFspecSetSeed(rfs, seed);
  
  int ntries=0;
  while(ntries < nperm*2 && psno.size() < nperm){
    int i1 = (int)floor(RFdrawVal(rfs)*n1);
    int i2 = (int)floor(RFdrawVal(rfs)*n2);
    // this makes sure that all entries are tried before going random
    // nothing in random to balance entries
    if(ntries < n1) i1 = ntries;
    if(ntries < n2) i2 = ntries;
    int hit = 0;
    for(int k=0; k < psno.size(); k++){
      if(psno[k][0]==i1 && psno[k][1]==i2) {
	hit =1;
	break;
      }
    }
    if(!hit){
      std::vector<int> nn = {i1,i2};
      psno.push_back(nn);
    }
    ntries++;
  }

  return(psno);
}


int DualPerm::mergedirs(const char *outdir, const char *srcdirs[], int nsrcdirs)
{
  char fname[2000];
  FILE *fp=NULL;

  long unsigned int seedlist[30];
  for(int n = 0; n < nsrcdirs; n++){
    sprintf(fname,"%s/seed.txt",srcdirs[n]);
    fp = fopen(fname,"r");
    fscanf(fp,"%lu\n",&seedlist[n]);
    fclose(fp);
    if(n>0 && seedlist[0] != seedlist[n]){
      printf("ERROR: merge(): seeds do not match %lu %lu\n",seedlist[0],seedlist[n]);
      return(1);
    }
  }

  int err = mkdir(outdir,0777);
  if(err != 0 && errno != EEXIST) {
    printf("ERROR: merge(): creating %s\n",outdir);
    return(err);
  }

  sprintf(fname,"%s/merge.txt",outdir);
  fp = fopen(fname,"w");
  for(int n = 0; n < nsrcdirs; n++) fprintf(fp,"%s\n",srcdirs[n]);
  fclose(fp);

  sprintf(fname,"%s/seed.txt",outdir);
  fp = fopen(fname,"w");
  fprintf(fp,"%lu\n",seedlist[0]);
  fclose(fp);

  SpatialCor sc;
  for(int n = 0; n < nsrcdirs; n++){
    SpatialCor scn;
    sprintf(fname,"%s/cc.info.dat",srcdirs[n]);
    scn.read(fname);
    sc.merge(scn);
  }
  sc.complete();
  sprintf(fname,"%s/cc.info.dat",outdir);
  sc.print(fname);

  std::vector<double> cc1, cc2, cc12;
  char kstr[3];
  for(int k=0; k<3; k++){
    if(k==0 || k==1) sprintf(kstr,"%d",k+1);
    else sprintf(kstr,"12");
    FILE *nfp[nsrcdirs];
    for(int n = 0; n < nsrcdirs; n++){
      sprintf(fname,"%s/cc%s.info.dat",srcdirs[n],kstr);
      nfp[n] = fopen(fname,"r");
      if(nfp[n] == NULL) {
	printf("ERROR: opening %s\n",fname);
	return(1);
      }
    }
    sprintf(fname,"%s/cc%s.info.dat",outdir,kstr);
    if(fio_FileExistsReadable(fname)) {
      int r = unlink(fname);
      if(r != 0 && r != ENOENT){
	printf("ERROR: %d deleting %s\n",r,fname);
	return(1);
      }
    }
    int stop = 0;
    while(!stop){
      SpatialCor scp;
      SpatialCor scn;
      for(int n = 0; n < nsrcdirs; n++){
	int nread = scn.readline(nfp[n]);
	if(nread == EOF){
	  stop = 1;
	  break;
	}
	scp.merge(scn);
      }
      if(stop) break;
      scp.nmask = sc.nmask;
      scp.complete();
      scp.printline(fname);
      if(k==0) cc1.push_back(scp.pcc);
      if(k==1) cc2.push_back(scp.pcc);
      if(k==2) cc12.push_back(scp.pcc);
    }
    for(int n = 0; n < nsrcdirs; n++) fclose(nfp[n]);
  }
  writepvals(outdir, sc.pcc, cc1, cc2, cc12);


  return(0);
}



int writepvals(const char *outdir, double cc0, std::vector<double> cc1,std::vector<double> cc2,std::vector<double> cc12)
{
  double p1pos=0, p1neg=0, p1abs=0, p2pos=0, p2neg=0, p2abs=0, p12pos=0, p12neg=0, p12abs=0;
  FILE *fp;
  char fname[2000];
  if(cc1.size() != 0){
    sprintf(fname,"%s/p1.dat",outdir);
    fp = fopen(fname,"w");
    p1pos = GetPVal(cc0, cc1, +1); fprintf(fp,"%31.28lf ",p1pos);
    p1neg = GetPVal(cc0, cc1, -1); fprintf(fp,"%31.28lf ",p1neg);
    p1abs = GetPVal(cc0, cc1,  0); fprintf(fp,"%31.28lf ",p1abs);
    fprintf(fp,"\n");
    fclose(fp);
  }
  if(cc2.size() != 0){
    sprintf(fname,"%s/p2.dat",outdir);
    fp = fopen(fname,"w");
    p2pos = GetPVal(cc0, cc2, +1); fprintf(fp,"%31.28lf ",p2pos);
    p2neg = GetPVal(cc0, cc2, -1); fprintf(fp,"%31.28lf ",p2neg);
    p2abs = GetPVal(cc0, cc2,  0); fprintf(fp,"%31.28lf ",p2abs);
    fprintf(fp,"\n");
    fclose(fp);
  }
  if(cc12.size() != 0){
    sprintf(fname,"%s/p12.dat",outdir);
    fp = fopen(fname,"w");
    p12pos = GetPVal(cc0, cc12, +1); fprintf(fp,"%31.28lf ",p12pos);
    p12neg = GetPVal(cc0, cc12, -1); fprintf(fp,"%31.28lf ",p12neg);
    p12abs = GetPVal(cc0, cc12,  0); fprintf(fp,"%31.28lf ",p12abs);
    fprintf(fp,"\n");
    fclose(fp);
    
    sprintf(fname,"%s/p12.max.dat",outdir);
    fp = fopen(fname,"w");
    if(p1pos > p2pos) fprintf(fp,"%31.28lf ",p1pos);
    else              fprintf(fp,"%31.28lf ",p2pos);
    if(p1neg > p2neg) fprintf(fp,"%31.28lf ",p1neg);
    else              fprintf(fp,"%31.28lf ",p2neg);
    if(p1abs > p2abs) fprintf(fp,"%31.28lf ",p1abs);
    else              fprintf(fp,"%31.28lf ",p2abs);
    fprintf(fp,"\n");
    fclose(fp);
  }
  return(0);
}




