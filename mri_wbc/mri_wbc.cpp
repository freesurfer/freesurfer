/**
 * @brief Computes whole-brain connectivity
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


/*!
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve
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
#include "mrisurf.h"
#include "romp_support.h"
#include "timer.h"
#include "mrimorph.h"
#include "fmriutils.h"
#include "fsenv.h"
#include "matfile.h"
#include "icosahedron.h"

double round(double x);

typedef struct {
  double distthresh;
  int nframes,voldim;
  double volres,*wf;
  int nshorttarg[2],nlongtarg[2];
  int nshort[2],nlong[2];
  int nshortvol,nlongvol;
  int v0,c0,r0,s0,c2;
  int ForceFail;
} WBCSYNTH;

typedef struct {
  MRI *fvol, *flh, *frh;
  int nframes;
  MRIS *lh, *rh, *lh2, *rh2;
  LABEL *lhlabel, *rhlabel;
  MRI *volmask, *lhmask, *rhmask, *lhvtxvol, *rhvtxvol;
  int nvolmask, nlhmask, nrhmask, ntot;
  MRI *volcon, *lhcon, *rhcon;
  MRI *volconS, *lhconS, *rhconS;
  MRI *volconL, *lhconL, *rhconL;
  MRI *volrhomean, *lhrhomean, *rhrhomean;
  MRI *coordtype, *vertexno, *xyz, *xyz2, *vvol;
  MRI *f, *fnorm;
  MRI *rhomean; 
  double rholist[100];
  int nrholist;
  int DoDist;
  double distthresh;
  MRI *con,*conS,*conL;
  int DoMat,DoTest;
  MATRIX *M;
  WBCSYNTH *wbcsynth;
} WBC;

MRI *WholeBrainCon(WBC *wbc);
int WBCfinish(WBC *wbc);
int WBCprep(WBC *wbc);
int WBCnframes(WBC *wbc);
int WBCsynth(WBC *wbc);
WBC *WBCtestSynth(WBC *wbc);
int WBCtestCheck(WBC *wbc);
int Index2UpperSubscript(int N, long i, int *r, int *c);
int Index2UpperSubscriptTest(int seed);

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

typedef struct {
  char *fvol;
  char *volmask;
  char *flh, *lhsurface, *lhsurface2, *lhlabel, *lhmask;
  char *frh, *rhsurface, *rhsurface2, *rhlabel, *rhmask;
  double rholist[100];
  int nrholist;
  double distthresh;
  int DoDist;
  char *outdir;
  char *volcon, *lhcon, *rhcon;
  char *volconS, *lhconS, *rhconS;
  char *volconL, *lhconL, *rhconL;
  char *volrhomean,*lhrhomean,*rhrhomean;
  int DoMat;
  char *matfile;
  int DoTest,ForceFail,SaveTest;
} CMDARGS;

CMDARGS *cmdargs;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err,n;
  WBC *wbc;

  cmdargs = (CMDARGS *)calloc(sizeof(CMDARGS),1);
  cmdargs->distthresh = 10;
  cmdargs->DoDist = 0;
  cmdargs->DoMat = 0;
  cmdargs->nrholist = 0;
  cmdargs->ForceFail = 0;

  nargs = handleVersionOption(argc, argv, "mri_wbc");
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

  if(cmdargs->outdir){
    err = mkdir(cmdargs->outdir,0777);
    if(err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n",cmdargs->outdir);
      perror(NULL);
      exit(1);
    }
  }

  wbc = (WBC *)calloc(sizeof(WBC),1);
  wbc->distthresh = cmdargs->distthresh;
  wbc->DoDist = cmdargs->DoDist;
  wbc->DoMat = cmdargs->DoMat;
  for(n=0; n < cmdargs->nrholist; n++)
    wbc->rholist[n] = cmdargs->rholist[n];
  wbc->nrholist = cmdargs->nrholist;
  wbc->DoTest = cmdargs->DoTest;

  if(wbc->DoTest == 0){
    if(cmdargs->fvol){
      wbc->fvol = MRIread(cmdargs->fvol);
      if(wbc->fvol == NULL) exit(1);
    }
    if(cmdargs->volmask){
      wbc->volmask = MRIread(cmdargs->volmask);
      if(wbc->volmask == NULL) exit(1);
    }
    if(cmdargs->flh){
      wbc->flh = MRIread(cmdargs->flh);
      if(wbc->flh == NULL) exit(1);
      wbc->lh  = MRISread(cmdargs->lhsurface);
      if(wbc->lh == NULL) exit(1);
      if(cmdargs->lhlabel){
	wbc->lhlabel = LabelRead(NULL, cmdargs->lhlabel);
	if(wbc->lhlabel == NULL) exit(1);
      }
      if(cmdargs->lhmask){
	wbc->lhmask = MRIread(cmdargs->lhmask);
	if(wbc->lhmask == NULL) exit(1);
      }
      if(cmdargs->lhsurface2) {
	wbc->lh2 = MRISread(cmdargs->lhsurface2);
	if(wbc->lh2 == NULL) exit(1);
      }
    }
    if(cmdargs->frh){
      wbc->frh = MRIread(cmdargs->frh);
      if(wbc->frh == NULL) exit(1);
      wbc->rh  = MRISread(cmdargs->rhsurface);
      if(wbc->rh == NULL) exit(1);
      if(cmdargs->rhlabel){
	wbc->rhlabel = LabelRead(NULL, cmdargs->rhlabel);
	if(wbc->rhlabel == NULL) exit(1);
      }
      if(cmdargs->rhmask){
	wbc->rhmask = MRIread(cmdargs->rhmask);
	if(wbc->rhmask == NULL) exit(1);
      }
      if(cmdargs->rhsurface2) {
	wbc->rh2 = MRISread(cmdargs->rhsurface2);
	if(wbc->rh2 == NULL) exit(1);
      }
    }
  }
  else {
    wbc->wbcsynth = (WBCSYNTH *)calloc(sizeof(WBCSYNTH),1);
    wbc->wbcsynth->volres = 3.5;
    wbc->wbcsynth->voldim = 32;
    wbc->wbcsynth->nframes = 5;
    wbc->wbcsynth->nshorttarg[0] = 10;
    wbc->wbcsynth->nlongtarg[0] = 20;
    wbc->wbcsynth->nshorttarg[1] = 15;
    wbc->wbcsynth->nlongtarg[1] = 25;
    wbc->wbcsynth->ForceFail = cmdargs->ForceFail;
    WBCtestSynth(wbc);
  }

  err = WBCnframes(wbc);
  if(err) exit(1);

  WBCprep(wbc);
  WholeBrainCon(wbc);
  WBCfinish(wbc);

  if(wbc->DoMat){
    printf("Writing matfile %s\n",cmdargs->matfile);
    err = MatlabWrite(wbc->M,cmdargs->matfile,"M");
    if(err) exit(1);    
  }
  if(wbc->DoTest){
    err = WBCtestCheck(wbc);
  }

  if(cmdargs->volcon){
    err = MRIwrite(wbc->volcon,cmdargs->volcon);
    if(err) exit(1);
    if(wbc->DoDist){
      err = MRIwrite(wbc->volconS,cmdargs->volconS);
      if(err) exit(1);
      err = MRIwrite(wbc->volconL,cmdargs->volconL);
      if(err) exit(1);
    }
  }
  if(cmdargs->lhcon){
    err = MRIwrite(wbc->lhcon,cmdargs->lhcon);
    if(err) exit(1);
    if(wbc->DoDist){
      err = MRIwrite(wbc->lhconS,cmdargs->lhconS);
      if(err) exit(1);
      err = MRIwrite(wbc->lhconL,cmdargs->lhconL);
      if(err) exit(1);
    }
  }
  if(cmdargs->rhcon){
    err = MRIwrite(wbc->rhcon,cmdargs->rhcon);
    if(err) exit(1);
    if(wbc->DoDist){
      err = MRIwrite(wbc->rhconS,cmdargs->rhconS);
      if(err) exit(1);
      err = MRIwrite(wbc->rhconL,cmdargs->rhconL);
      if(err) exit(1);
    }
  }

  if(cmdargs->volrhomean){
    err = MRIwrite(wbc->volrhomean,cmdargs->volrhomean);
    if(err) exit(1);
  }
  if(cmdargs->lhrhomean){
    err = MRIwrite(wbc->lhrhomean,cmdargs->lhrhomean);
    if(err) exit(1);
  }
  if(cmdargs->rhrhomean){
    err = MRIwrite(wbc->rhrhomean,cmdargs->rhrhomean);
    if(err) exit(1);
  }

  exit(0);
}

/* -------------------------------------------------------- */
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

    else if (!strcasecmp(option, "--fvol")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->fvol = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->outdir = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--volmask")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volmask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lh")) {
      if(nargc < 2) CMDargNErr(option,2);
      cmdargs->flh = pargv[0];
      cmdargs->lhsurface = pargv[1];
      nargsused = 2;
      if(CMDnthIsArg(nargc, pargv, 2)) {
	cmdargs->lhsurface = pargv[2];
        nargsused ++;
      } 
    } 
    else if (!strcasecmp(option, "--lhmask")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhmask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhlabel")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhlabel = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rh")) {
      if(nargc < 2) CMDargNErr(option,2);
      cmdargs->frh = pargv[0];
      cmdargs->rhsurface = pargv[1];
      nargsused = 2;
      if(CMDnthIsArg(nargc, pargv, 2)) {
	cmdargs->rhsurface = pargv[2];
        nargsused ++;
      } 
    } 
    else if (!strcasecmp(option, "--rhmask")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhmask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhlabel")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhlabel = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rho")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->rholist[cmdargs->nrholist]);
      cmdargs->nrholist++;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--dist")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->distthresh);
      cmdargs->DoDist = 1;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--volrhomean")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volrhomean = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhrhomean")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhrhomean = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhrhomean")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhrhomean = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--volcon")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volcon = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhcon")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhcon = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhcon")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhcon = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--volconS")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volconS = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhconS")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhconS = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhconS")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhconS = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--volconL")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volconL = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhconL")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhconL = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhconL")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhconL = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mat")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->matfile = pargv[0];
      cmdargs->DoMat = 1;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--uppersub-test")) {
      if(nargc < 1) CMDargNErr(option,1);
      int N, err;
      sscanf(pargv[0],"%d",&N);
      err = Index2UpperSubscriptTest(N);
      exit(err);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--test")) {
      cmdargs->DoTest = 1;
      nargsused = 0;
    }
    else if (!strcasecmp(option, "--test-fail")) {
      cmdargs->DoTest = 1;
      cmdargs->ForceFail = 1;
      nargsused = 0;
    }
    else if (!strcasecmp(option, "--save-test")) {
      cmdargs->DoTest = 1;
      cmdargs->SaveTest = 1;
      nargsused = 0;
    }
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      int nthreads;
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef HAVE_OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag-show"))    Gdiag = (Gdiag & DIAG_SHOW);
    else if (!strcasecmp(option, "--diag-verbose")) Gdiag = (Gdiag & DIAG_VERBOSE);
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
/* -------------------------------------------------------- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --fvol volfile : functional volume\n");
  printf("   --volmask volfile : mask for functional volume \n");
  printf("\n");
  printf("   --lh lhffile lh.surface <lh.inflated>: functional surface lh\n");
  printf("   --lhmask lhfile : mask for lh functional surface\n");
  printf("   --lhlabel lhlabel : label mask for lh functional surface\n");
  printf("\n");
  printf("   --rh rhffile rh.surface <rh.inflated>: functional surface rh\n");
  printf("   --rhmask rhfile : mask for rh functional surface\n");
  printf("   --rhlabel rhlabel : label mask for rh functional surface\n");
  printf("\n");
  printf("   --rho rhothresh\n");
  printf("   --dist distthresh\n");
  printf("\n");
  printf("   --threads nthreads\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* -------------------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void check_options(void) {
  char tmpstr[2000];
  if((cmdargs->outdir == NULL && !cmdargs->DoTest) || 
     (cmdargs->outdir == NULL && cmdargs->SaveTest)){
    printf("ERROR: no output\n");
    exit(1);
  }
  if(cmdargs->fvol == NULL && cmdargs->flh == NULL && 
     cmdargs->frh == NULL && !cmdargs->DoTest){
    printf("ERROR: no input\n");
    exit(1);
  }
  if(cmdargs->fvol != NULL || cmdargs->SaveTest){
    sprintf(tmpstr,"%s/vol.con.nii.gz",cmdargs->outdir);
    cmdargs->volcon = strcpyalloc(tmpstr);
    if(cmdargs->DoDist){
      sprintf(tmpstr,"%s/vol.conS.nii.gz",cmdargs->outdir);
      cmdargs->volconS = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/vol.conL.nii.gz",cmdargs->outdir);
      cmdargs->volconL = strcpyalloc(tmpstr);
    }
  }
  if(cmdargs->flh != NULL || cmdargs->SaveTest){
    sprintf(tmpstr,"%s/lh.con.nii.gz",cmdargs->outdir);
    cmdargs->lhcon = strcpyalloc(tmpstr);
    if(cmdargs->DoDist){
      sprintf(tmpstr,"%s/lh.conS.nii.gz",cmdargs->outdir);
      cmdargs->lhconS = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/lh.conL.nii.gz",cmdargs->outdir);
      cmdargs->lhconL = strcpyalloc(tmpstr);
    }
  }
  if(cmdargs->frh != NULL || cmdargs->SaveTest){
    sprintf(tmpstr,"%s/rh.con.nii.gz",cmdargs->outdir);
    cmdargs->rhcon = strcpyalloc(tmpstr);
    if(cmdargs->DoDist){
      sprintf(tmpstr,"%s/rh.conS.nii.gz",cmdargs->outdir);
      cmdargs->rhconS = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/rh.conL.nii.gz",cmdargs->outdir);
      cmdargs->rhconL = strcpyalloc(tmpstr);
    }
  }
  if(cmdargs->fvol != NULL || cmdargs->SaveTest){
    sprintf(tmpstr,"%s/vol.rho.mean.nii.gz",cmdargs->outdir);
    cmdargs->volrhomean = strcpyalloc(tmpstr);
  }
  if(cmdargs->flh != NULL || cmdargs->SaveTest){
    sprintf(tmpstr,"%s/lh.rho.mean.nii.gz",cmdargs->outdir);
    cmdargs->lhrhomean = strcpyalloc(tmpstr);
  }
  if(cmdargs->frh != NULL || cmdargs->SaveTest){
    sprintf(tmpstr,"%s/rh.rho.mean.nii.gz",cmdargs->outdir);
    cmdargs->rhrhomean = strcpyalloc(tmpstr);
  }
  if(cmdargs->nrholist == 0){
    cmdargs->rholist[0] = 0.2;
    cmdargs->nrholist = 1;
  }

  return;
}
/* -------------------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  return;
}

/*******************************************************************************/
MRI *WholeBrainCon(WBC *wbc)
{
  MRI **conth, **conSth=NULL, **conLth=NULL, **rhomean;
  int nthreads,nthrho;
  long nperthread0,*nperthread,ntot;
  int threadno;
  long npairs,ia,ib;
  int k1; // note: these are redefined in loop
  Timer timer;
  double **pf, val;
  int *k1a,*k1b,*k2a,*k2b;

  // Load time courses into a pointer structure
  pf = (double **) calloc(sizeof(double *),wbc->ntot);
  for(k1 = 0; k1 < wbc->ntot; k1++){
    int t;
    pf[k1] = (double *) calloc(sizeof(double),wbc->f->nframes);
    for(t=0; t < wbc->f->nframes; t++) 
      pf[k1][t] = MRIgetVoxVal(wbc->fnorm,k1,0,0,t);    
  }

  // This is for testing. In most apps, M will be huge
  if(wbc->DoMat) wbc->M = MatrixAlloc(wbc->ntot,wbc->ntot,MATRIX_REAL);

  nthreads = 1;
  #ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();
  #endif
  npairs = (long)wbc->ntot*(wbc->ntot-1)/2;
  nperthread0 = (long)round((double)npairs/nthreads - 1); // don't use nint(), need long
  if(nperthread0 <= 0) nperthread0 = 1;
  printf("ntot = %d, nthreads = %d, npairs = %ld, nperthread0 = %ld\n",
	 wbc->ntot,nthreads,npairs,nperthread0);

  // Compute number of pairs per thread
  nperthread = (long *) calloc(sizeof(long),nthreads);
  ntot = 0;
  for(threadno=0; threadno < nthreads; threadno++){
    nperthread[threadno] = nperthread0; // number of pairs per thread
    ntot += nperthread0;
  }
  if(ntot != npairs) nperthread[nthreads-1] += npairs-ntot;

  // Allocate con volumes for each thread
  rhomean = (MRI **) calloc(sizeof(MRI*),nthreads);
  conth   = (MRI **) calloc(sizeof(MRI*),nthreads);
  if(wbc->DoDist){
    conSth = (MRI **) calloc(sizeof(MRI*),nthreads);
    conLth = (MRI **) calloc(sizeof(MRI*),nthreads);
  }
  for(threadno=0; threadno < nthreads; threadno++){
    rhomean[threadno] = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, 1);
    conth[threadno] = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
    if(wbc->DoDist){
      conSth[threadno] = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
      conLth[threadno] = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
    }
  }

  // This is just a test
  ntot = 0;
  for(threadno=0; threadno < nthreads; threadno++){
    printf("thread %d %ld\n",threadno,nperthread[threadno]);
    ntot += nperthread[threadno];
  }
  printf("ntotpairs = %ld vs npairs %ld, diff %ld\n",ntot,npairs,ntot-npairs); // should be equal

  k1a = (int *)calloc(sizeof(int),nthreads);
  k1b = (int *)calloc(sizeof(int),nthreads);
  k2a = (int *)calloc(sizeof(int),nthreads);
  k2b = (int *)calloc(sizeof(int),nthreads);
  ia = 0;
  for(threadno=0; threadno < nthreads; threadno++){
    ib = ia + nperthread[threadno];
    Index2UpperSubscript(wbc->ntot, ia, &k1a[threadno], &k2a[threadno]);
    Index2UpperSubscript(wbc->ntot, ib-1, &k1b[threadno], &k2b[threadno]);
    printf("thread %d %12ld %12ld   %7d %7d   %7d %7d\n",threadno,ia,ib,
	   k1a[threadno],k2a[threadno], k1b[threadno],k2b[threadno]);
    ia = ib;
  }

  printf("rho thresholds (%d): ",wbc->nrholist);
  for(nthrho = 0; nthrho < wbc->nrholist; nthrho++)
    printf("%lf ",wbc->rholist[nthrho]);
  printf("\n");
  
  printf("Starting WBC loop\n"); fflush(stdout);
  timer.reset();
  ROMP_PF_begin
  #ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) 
  #endif
  for(threadno = 0; threadno < nthreads; threadno ++){
    ROMP_PFLB_begin
    int  k1, k2, t, n, thno, q, ct1, ct2;
    int k1start,k1stop,k2start,k2stop,k2min,k2max,nthrho;
    long nthpair;
    double rho,dx,dy,dz,dist,*pf1,*pf2,x1,y1,z1,x2,y2,z2;
    double rhothresh;
    MRI *conDth;
    thno = threadno;
    #ifdef HAVE_OPENMP
    thno = omp_get_thread_num(); // actual thread number
    #endif

    k1start = k1a[thno];
    k1stop  = k1b[thno];
    k2start = k2a[thno];
    k2stop  = k2b[thno];

    nthpair = 0;
    for(k1=k1start; k1<=k1stop; k1++){
      if(k1 == k1start) k2min = k2start;
      else              k2min = k1+1;
      if(k1 == k1stop)  k2max = k2stop;
      else              k2max = wbc->ntot-1;
      for(k2=k2min; k2<=k2max; k2++){
	nthpair ++;
	if(thno == 0 && (nthpair == 1 || (nthpair % (long)1e8 == 0))){
	  printf("%4.1f%%, t=%5.1f min\n",nthreads*100.0*nthpair/npairs, timer.minutes());
	  fflush(stdout);
	}

	// compute rho
	rho = 0;
	pf1 = pf[k1];
	pf2 = pf[k2];
	for(t=0; t < wbc->f->nframes; t++){
	  rho += (*pf1) * (*pf2);
	  pf1++;
	  pf2++;
	}
	if(wbc->M != NULL){
	  wbc->M->rptr[k1+1][k2+1] = rho;
	  wbc->M->rptr[k2+1][k1+1] = rho;
	}

	val = MRIgetVoxVal(rhomean[thno],k1,0,0,0);
	MRIsetVoxVal(rhomean[thno],k1,0,0,0,val+rho);

	for(nthrho = 0; nthrho < wbc->nrholist; nthrho++){
	  rhothresh = wbc->rholist[nthrho];
	  if(fabs(rho) < rhothresh) continue;
	
	  n = MRIgetVoxVal(conth[thno],k1,0,0,nthrho);
	  MRIsetVoxVal(conth[thno],k1,0,0,nthrho,n+1);
	  n = MRIgetVoxVal(conth[thno],k2,0,0,nthrho);
	  MRIsetVoxVal(conth[thno],k2,0,0,nthrho,n+1);
	} // rho thresh
	
	if(!wbc->DoDist) continue;
	  
	ct1 = MRIgetVoxVal(wbc->coordtype,k1,0,0,0);
	ct2 = MRIgetVoxVal(wbc->coordtype,k2,0,0,0);
	if(ct1 == ct2){
	  x1 = MRIgetVoxVal(wbc->xyz,k1,0,0,0);
	  y1 = MRIgetVoxVal(wbc->xyz,k1,0,0,1);
	  z1 = MRIgetVoxVal(wbc->xyz,k1,0,0,2);
	  x2 = MRIgetVoxVal(wbc->xyz,k2,0,0,0);
	  y2 = MRIgetVoxVal(wbc->xyz,k2,0,0,1);
	  z2 = MRIgetVoxVal(wbc->xyz,k2,0,0,2);
	  dx = x1-x2;
	  dy = y1-y2;
	  dz = z1-z2;
	  dist = sqrt(dx*dx + dy*dy + dz*dz);
	  if(dist < wbc->distthresh) q = 1; // short dist
	  else                       q = 2; // long dist
	  if(q==1 && ((wbc->lh2 && ct1 == 1) || (wbc->rh2 && ct1 == 2))){
	    /* Surface point that has been declared short using the
	    first surface (prob ?h.white). Now check with the second
	    surface (prob ?h.inflated). This can only cause a short to
	    be recoded to long. The problem this is attempting to
	    solve is that it is very expensive to compute the actual
	    distance between two surface points. The actual distance is
	    always >= the euclidean dist between the ?h.white points, so,
	    if that distance is > thresh, it must be a long connection.
	    If it is shorter, then it could be a short connection or 
	    it could be a long with a short euclidean (eg, jumping across
	    a gyrus or sulcus) so this checks the inflated. The problem
	    with the inflated is that it is not distortion-free.
	    */
	    x1 = MRIgetVoxVal(wbc->xyz2,k1,0,0,0);
	    y1 = MRIgetVoxVal(wbc->xyz2,k1,0,0,1);
	    z1 = MRIgetVoxVal(wbc->xyz2,k1,0,0,2);
	    x2 = MRIgetVoxVal(wbc->xyz2,k2,0,0,0);
	    y2 = MRIgetVoxVal(wbc->xyz2,k2,0,0,1);
	    z2 = MRIgetVoxVal(wbc->xyz2,k2,0,0,2);
	    dx = x1-x2;
	    dy = y1-y2;
	    dz = z1-z2;
	    dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist >= wbc->distthresh) q = 2; // long dist
	  }
	}
	else q = 2; // long dist
	if(q == 1) conDth = conSth[thno];
	if(q == 2) conDth = conLth[thno];

	for(nthrho = 0; nthrho < wbc->nrholist; nthrho++){
	  rhothresh = wbc->rholist[nthrho];
	  if(fabs(rho) < rhothresh) continue;
	  n = MRIgetVoxVal(conDth,k1,0,0,nthrho);
	  MRIsetVoxVal(conDth,k1,0,0,nthrho,n+1);
	  n = MRIgetVoxVal(conDth,k2,0,0,nthrho);
	  MRIsetVoxVal(conDth,k2,0,0,nthrho,n+1);
	} // rho thresh

      } // k2
    } // k1
    ROMP_PFLB_end
  } // thread
  ROMP_PF_end

  // Sum up the threads
  for(threadno=0; threadno < nthreads; threadno++){
    MRIadd(wbc->rhomean,rhomean[threadno],wbc->rhomean);
    MRIadd(wbc->con,conth[threadno],wbc->con);
    if(wbc->DoDist){
      MRIadd(wbc->conS,conSth[threadno],wbc->conS);
      MRIadd(wbc->conL,conLth[threadno],wbc->conL);
    }
  }

  // Divide number of connections by total possible
  printf("Scaling by ntot-1 %d\n",wbc->ntot-1);
  MRImultiplyConst(wbc->rhomean,1.0/(wbc->ntot-1),wbc->rhomean);
  MRImultiplyConst(wbc->con,1.0/(wbc->ntot-1),wbc->con);
  if(wbc->DoDist){
    MRImultiplyConst(wbc->conS,1.0/(wbc->ntot-1),wbc->conS);
    MRImultiplyConst(wbc->conL,1.0/(wbc->ntot-1),wbc->conL);
  }

  // Clean up
  for(threadno=0; threadno < nthreads; threadno++){
    MRIfree(&rhomean[threadno]);
    MRIfree(&conth[threadno]);
    if(wbc->DoDist){
      MRIfree(&conSth[threadno]);
      MRIfree(&conLth[threadno]);
    }
  }
  free(k1a);
  free(k1b);
  free(k2a);
  free(k2b);
  free(conth);
  for(k1 = 0; k1 < wbc->ntot; k1++) free(pf[k1]);
  free(pf);

  printf(" t = %6.4f, ntot = %d\n", timer.seconds(), wbc->ntot);fflush(stdout);
  return(wbc->con);
}

int WBCprep(WBC *wbc)
{
  int nthvox, c, r, s, t, k, vtxno;
  MATRIX *V, *crs, *xyz=NULL;
  double val,voxvolume;

  if(wbc->fvol)     wbc->nframes = wbc->fvol->nframes;
  else if(wbc->flh) wbc->nframes = wbc->flh->nframes;
  else if(wbc->frh) wbc->nframes = wbc->frh->nframes;

  if(wbc->lhlabel) wbc->lhmask = MRISlabel2Mask(wbc->lh,wbc->lhlabel,NULL);
  if(wbc->rhlabel) wbc->rhmask = MRISlabel2Mask(wbc->rh,wbc->rhlabel,NULL);

  wbc->nvolmask = 0;
  wbc->nlhmask = 0;
  wbc->nrhmask = 0;
  if(wbc->fvol){
    if(wbc->volmask) wbc->nvolmask = MRIcountAboveThreshold(wbc->volmask, 0.5);
    else wbc->nvolmask = wbc->fvol->width * wbc->fvol->height * wbc->fvol->depth;
  }
  if(wbc->flh){
    if(wbc->lhmask) wbc->nlhmask  = MRIcountAboveThreshold(wbc->lhmask, 0.5);
    else            wbc->nlhmask  = wbc->flh->width;
  }
  if(wbc->frh){
    if(wbc->rhmask) wbc->nrhmask  = MRIcountAboveThreshold(wbc->rhmask, 0.5);
    else            wbc->nrhmask  = wbc->frh->width;
  }

  wbc->ntot = wbc->nvolmask + wbc->nlhmask + wbc->nrhmask;
  printf("nmask %d %d %d   %d\n",wbc->nvolmask,wbc->nlhmask,wbc->nrhmask,wbc->ntot);

  wbc->coordtype = MRIalloc(wbc->ntot,1,1,MRI_INT);
  wbc->vvol      = MRIallocSequence(wbc->ntot,1,1,MRI_FLOAT,1); // vertex/voxel volume
  wbc->vertexno  = MRIalloc(wbc->ntot,1,1,MRI_INT);
  wbc->xyz       = MRIallocSequence(wbc->ntot,1,1,MRI_FLOAT,3);
  if(wbc->lh2 || wbc->rh2) wbc->xyz2 = MRIallocSequence(wbc->ntot,1,1,MRI_FLOAT,3);
  wbc->f = MRIallocSequence(wbc->ntot,1,1,MRI_FLOAT,wbc->nframes);

  nthvox = 0;
  if(wbc->fvol){
    // Not sure it is necessary/useful to map into a particular mm space as
    // long as it is a mm space so the distance threshold is right. This tkreg
    // space does not match that of the surface, but this should be ok becase
    // we never compute the distance from volume to surface, it is always
    // long distance. If this turns out not to be the case, then will need
    // a registraion matrix
    voxvolume = wbc->fvol->xsize * wbc->fvol->ysize * wbc->fvol->zsize;
    V = MRIxfmCRS2XYZtkreg(wbc->fvol);
    crs = MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
    for(s=0; s < wbc->fvol->depth; s++){
      for(c=0; c < wbc->fvol->width; c++){
	for(r=0; r < wbc->fvol->height; r++){
	  if(wbc->volmask && MRIgetVoxVal(wbc->volmask,c,r,s,0) < 0.5) continue;
	  MRIsetVoxVal(wbc->coordtype,nthvox,0,0,0, 0);
	  MRIsetVoxVal(wbc->vvol,nthvox,0,0,0, voxvolume);
	  crs->rptr[1][1] = c;
	  crs->rptr[2][1] = r;
	  crs->rptr[3][1] = s;
	  xyz = MatrixMultiplyD(V,crs,xyz);
	  for(k=0; k < 3; k++) MRIsetVoxVal(wbc->xyz,nthvox,0,0,k, xyz->rptr[k+1][1]);
	  for(t=0; t < wbc->nframes; t++){
	    val = MRIgetVoxVal(wbc->fvol,c,r,s,t);
	    MRIsetVoxVal(wbc->f,nthvox,0,0,t,val);
	  } // time
	  nthvox ++;
	} // row
      } // col
    } // slice
    MatrixFree(&V);
    MatrixFree(&xyz);
    MatrixFree(&crs);    
  }
  if(wbc->flh){
    for(vtxno = 0; vtxno < wbc->flh->width; vtxno++){
      if(wbc->lhmask && MRIgetVoxVal(wbc->lhmask,vtxno,0,0,0) < 0.5) continue;
      MRIsetVoxVal(wbc->vertexno,nthvox,0,0,0, vtxno);
      MRIsetVoxVal(wbc->coordtype,nthvox,0,0,0, 1);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,0, wbc->lh->vertices[vtxno].x);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,1, wbc->lh->vertices[vtxno].y);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,2, wbc->lh->vertices[vtxno].z);
      if(wbc->lh2){
	MRIsetVoxVal(wbc->xyz2,nthvox,0,0,0, wbc->lh->vertices[vtxno].x);
	MRIsetVoxVal(wbc->xyz2,nthvox,0,0,1, wbc->lh->vertices[vtxno].y);
	MRIsetVoxVal(wbc->xyz2,nthvox,0,0,2, wbc->lh->vertices[vtxno].z);
      }
      for(t=0; t < wbc->flh->nframes; t++){
	val = MRIgetVoxVal(wbc->flh,vtxno,0,0,t);
	MRIsetVoxVal(wbc->f,nthvox,0,0,t,val);
      } // time
      nthvox ++;
    } // vertex
  }
  if(wbc->frh){
    for(vtxno = 0; vtxno < wbc->frh->width; vtxno++){
      if(wbc->rhmask && MRIgetVoxVal(wbc->rhmask,vtxno,0,0,0) < 0.5) continue;
      MRIsetVoxVal(wbc->vertexno,nthvox,0,0,0, vtxno);
      MRIsetVoxVal(wbc->coordtype,nthvox,0,0,0, 2);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,0, wbc->rh->vertices[vtxno].x);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,1, wbc->rh->vertices[vtxno].y);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,2, wbc->rh->vertices[vtxno].z);
      if(wbc->rh2){
	MRIsetVoxVal(wbc->xyz2,nthvox,0,0,0, wbc->rh->vertices[vtxno].x);
	MRIsetVoxVal(wbc->xyz2,nthvox,0,0,1, wbc->rh->vertices[vtxno].y);
	MRIsetVoxVal(wbc->xyz2,nthvox,0,0,2, wbc->rh->vertices[vtxno].z);
      }
      for(t=0; t < wbc->frh->nframes; t++){
	val = MRIgetVoxVal(wbc->frh,vtxno,0,0,t);
	MRIsetVoxVal(wbc->f,nthvox,0,0,t,val);
      } // time
      nthvox ++;
    } // vertex
  }

  // normalize by stddev
  wbc->fnorm = MRIframeNorm(wbc->f, NULL, NULL);  

  wbc->rhomean = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, 1);
  wbc->con = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
  if(wbc->DoDist){
    wbc->conS = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
    wbc->conL = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
  }

  return(0);
}

int WBCfinish(WBC *wbc)
{
  int nthvox, c,r,s,vtxno, nthrho;
  double val;

  nthvox = 0;
  if(wbc->fvol){
    wbc->volrhomean = MRIallocSequence(wbc->fvol->width,wbc->fvol->height,wbc->fvol->depth,MRI_FLOAT,1);
    MRIcopyHeader(wbc->fvol,wbc->volrhomean);
    MRIcopyPulseParameters(wbc->fvol,wbc->volrhomean);
    wbc->volcon = MRIallocSequence(wbc->fvol->width,wbc->fvol->height,wbc->fvol->depth,MRI_FLOAT,wbc->nrholist);
    MRIcopyHeader(wbc->fvol,wbc->volcon);
    MRIcopyPulseParameters(wbc->fvol,wbc->volcon);
    if(wbc->DoDist){
      wbc->volconS = MRIallocSequence(wbc->fvol->width,wbc->fvol->height,wbc->fvol->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->fvol,wbc->volconS);
      MRIcopyPulseParameters(wbc->fvol,wbc->volconS);
      wbc->volconL = MRIallocSequence(wbc->fvol->width,wbc->fvol->height,wbc->fvol->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->fvol,wbc->volconL);
      MRIcopyPulseParameters(wbc->fvol,wbc->volconL);
    }
    for(s=0; s < wbc->fvol->depth; s++){
      for(c=0; c < wbc->fvol->width; c++){
	for(r=0; r < wbc->fvol->height; r++){
	  if(wbc->volmask && MRIgetVoxVal(wbc->volmask,c,r,s,0) < 0.5) continue;
	  val = MRIgetVoxVal(wbc->rhomean,nthvox,0,0,0);
	  MRIsetVoxVal(wbc->volrhomean,c,r,s,0,val);
	  for(nthrho=0; nthrho < wbc->nrholist; nthrho++){
	    val = MRIgetVoxVal(wbc->con,nthvox,0,0,nthrho);
	    MRIsetVoxVal(wbc->volcon,c,r,s,nthrho,val);
	    if(wbc->DoDist){
	      val = MRIgetVoxVal(wbc->conS,nthvox,0,0,nthrho);
	      MRIsetVoxVal(wbc->volconS,c,r,s,nthrho,val);
	      val = MRIgetVoxVal(wbc->conL,nthvox,0,0,nthrho);
	      MRIsetVoxVal(wbc->volconL,c,r,s,nthrho,val);
	    }
	  } // rho
	  nthvox ++;
	} // row
      } // col
    } // slice
  }
  if(wbc->flh){
    wbc->lhrhomean = MRIallocSequence(wbc->flh->width,wbc->flh->height,wbc->flh->depth,MRI_FLOAT,1);
    MRIcopyHeader(wbc->flh,wbc->lhrhomean);
    MRIcopyPulseParameters(wbc->flh,wbc->lhrhomean);
    wbc->lhcon = MRIallocSequence(wbc->flh->width,wbc->flh->height,wbc->flh->depth,MRI_FLOAT,wbc->nrholist);
    MRIcopyHeader(wbc->flh,wbc->lhcon);
    MRIcopyPulseParameters(wbc->flh,wbc->lhcon);
    if(wbc->DoDist){
      wbc->lhconS = MRIallocSequence(wbc->flh->width,wbc->flh->height,wbc->flh->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->flh,wbc->lhconS);
      MRIcopyPulseParameters(wbc->flh,wbc->lhconS);
      wbc->lhconL = MRIallocSequence(wbc->flh->width,wbc->flh->height,wbc->flh->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->flh,wbc->lhconL);
      MRIcopyPulseParameters(wbc->flh,wbc->lhconL);
    }
    for(vtxno = 0; vtxno < wbc->flh->width; vtxno++){
      if(wbc->lhmask && MRIgetVoxVal(wbc->lhmask,vtxno,0,0,0) < 0.5) continue;
      val = MRIgetVoxVal(wbc->rhomean,nthvox,0,0,0);
      MRIsetVoxVal(wbc->lhrhomean,vtxno,0,0,0,val);
      for(nthrho=0; nthrho < wbc->nrholist; nthrho++){
	val = MRIgetVoxVal(wbc->con,nthvox,0,0,nthrho);
	MRIsetVoxVal(wbc->lhcon,vtxno,0,0,nthrho,val);
	if(wbc->DoDist){
	  val = MRIgetVoxVal(wbc->conS,nthvox,0,0,nthrho);
	  MRIsetVoxVal(wbc->lhconS,vtxno,0,0,nthrho,val);
	  val = MRIgetVoxVal(wbc->conL,nthvox,0,0,nthrho);
	  MRIsetVoxVal(wbc->lhconL,vtxno,0,0,nthrho,val);
	}
      } // rho
      nthvox ++;
    } // vertex
  }
  if(wbc->frh){
    wbc->rhrhomean = MRIallocSequence(wbc->frh->width,wbc->frh->height,wbc->frh->depth,MRI_FLOAT,1);
    MRIcopyHeader(wbc->frh,wbc->rhrhomean);
    MRIcopyPulseParameters(wbc->frh,wbc->rhrhomean);
    wbc->rhcon = MRIallocSequence(wbc->frh->width,wbc->frh->height,wbc->frh->depth,MRI_FLOAT,wbc->nrholist);
    MRIcopyHeader(wbc->frh,wbc->rhcon);
    MRIcopyPulseParameters(wbc->frh,wbc->rhcon);
    if(wbc->DoDist){
      wbc->rhconS = MRIallocSequence(wbc->frh->width,wbc->frh->height,wbc->frh->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->frh,wbc->rhconS);
      MRIcopyPulseParameters(wbc->frh,wbc->rhconS);
      wbc->rhconL = MRIallocSequence(wbc->frh->width,wbc->frh->height,wbc->frh->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->frh,wbc->rhconL);
      MRIcopyPulseParameters(wbc->frh,wbc->rhconL);
    }
    for(vtxno = 0; vtxno < wbc->frh->width; vtxno++){
      if(wbc->rhmask && MRIgetVoxVal(wbc->rhmask,vtxno,0,0,0) < 0.5) continue;
      val = MRIgetVoxVal(wbc->rhomean,nthvox,0,0,0);
      MRIsetVoxVal(wbc->rhrhomean,vtxno,0,0,0,val);
      for(nthrho=0; nthrho < wbc->nrholist; nthrho++){
	val = MRIgetVoxVal(wbc->con,nthvox,0,0,nthrho);
	MRIsetVoxVal(wbc->rhcon,vtxno,0,0,nthrho,val);
	if(wbc->DoDist){
	  val = MRIgetVoxVal(wbc->conS,nthvox,0,0,nthrho);
	  MRIsetVoxVal(wbc->rhconS,vtxno,0,0,nthrho,val);
	  val = MRIgetVoxVal(wbc->conL,nthvox,0,0,nthrho);
	  MRIsetVoxVal(wbc->rhconL,vtxno,0,0,nthrho,val);
	}
      } // rho
      nthvox ++;
    } // vertex
  }
  return(0);
}

/*!
  \fn int Index2UpperSubscript(int N, long i, int *r, int *c)
  \brief Computes the row and col of the ith index in an upper
  triangular N-x-N matrix ignoring the diagonal and lower triangular
  components. The application for this involves evaluating pairs
  of items in a list excluding self-comparisons (diag) and 
  reverse comparisons (lower). See also Index2UpperSubscriptTest().
  Row, col, and index are all 0-based.
 */
int Index2UpperSubscript(int N, long i, int *r, int *c)
{
  long i1,r1,c1,ir1;
  i1 = i + 1;
  r1 = ceil((-(1-2.0*N) - sqrt( pow((1-2.0*N),2) - 8.0*i1))/2.0);
  ir1 = N*r1 - (r1*(r1+1))/2;
  c1 = N - (ir1-i1);
  *r = r1 - 1;
  *c = c1 - 1;
  return(0);
}

/*!
  \fn int Index2UpperSubscriptTest(int N)
  \brief Test for Index2UpperSubscript()
 */
int Index2UpperSubscriptTest(int N)
{
  int r, c, rt, ct;
  long i,err;

  printf("Index2UpperSubscriptTest(): N = %d\n",N);
  err = 0;
  i = 0;
  for(r=0; r < N-1; r++){
    for(c=r+1; c < N; c++){
      Index2UpperSubscript(N, i, &rt, &ct);
      if(r != rt || c != ct) {
	err++;
	printf("ERROR: %4d %4d %4ld  %4d %4d\n",r,c,i,rt,ct);
      }
      i = i + 1;
    }
  }
  printf("Found %ld errors\n",err);
  if(err) return(1);

  return(0);
}

/*!
  \fn int WBCnframes(WBC *wbc)
  \brief get the number of frames and checks that all are
  consistent. Returns 0 if ok, 1 if an error.
 */
int WBCnframes(WBC *wbc)
{
  if(wbc->fvol){
    wbc->nframes = wbc->fvol->nframes;
    if(wbc->flh){
      if(wbc->fvol->nframes != wbc->flh->nframes){
	printf("ERROR: nframes mismatch fvol and flh %d %d\n",
	       wbc->fvol->nframes,wbc->flh->nframes);
	return(1);
      }
    }
    if(wbc->frh){
      if(wbc->fvol->nframes != wbc->frh->nframes){
	printf("ERROR: nframes mismatch fvol and frh %d %d\n",
	       wbc->fvol->nframes,wbc->frh->nframes);
	return(1);
      }
    }
  }
  if(wbc->flh && wbc->frh){
    wbc->nframes = wbc->flh->nframes;
    if(wbc->flh->nframes != wbc->frh->nframes){
      printf("ERROR: nframes mismatch flh and frh %d %d\n",
	     wbc->flh->nframes,wbc->frh->nframes);
      return(1);
    }
  }
  printf("nframes = %d\n",wbc->nframes);
  return(0);
}

WBC *WBCtestSynth(WBC *wbc)
{
  FSENV *fsenv;
  double x0,y0,z0,x,y,z,dx,dy,dz,d,volres,*wf,dmax;
  int nshort,nlong,t,vno,hemi,voldim,nframes,c0,r0,s0,c2;
  char tmpstr[2000];
  const char *hemistr;
  MRIS *surf,*surf2;
  MRI *func,*mask;

  fsenv = FSENVgetenv();
  volres = wbc->wbcsynth->volres;
  voldim = wbc->wbcsynth->voldim;
  nframes = wbc->wbcsynth->nframes;
  wbc->nframes = nframes;

  if(wbc->wbcsynth->ForceFail) wbc->wbcsynth->distthresh = 2*wbc->distthresh;
  else                         wbc->wbcsynth->distthresh = wbc->distthresh;

  dmax = 0;
  //if(fabs(wbc->distthresh - wbc->wbcsynth->distthresh) < .00001) dmax = 10e10;

  // create synthetic time seriesa
  wf = (double *) calloc(sizeof(double),wbc->nframes);
  for(t=0; t < wbc->nframes; t++) wf[t] = drand48()-0.5;
  wbc->wbcsynth->wf = wf;

  wbc->fvol = MRIallocSequence(voldim,voldim,voldim,MRI_FLOAT,wbc->nframes);
  wbc->fvol->xsize = volres;
  wbc->fvol->ysize = volres;
  wbc->fvol->zsize = volres;
  wbc->fvol->tr    = 2000;

  wbc->volmask = MRIallocSequence(voldim,voldim,voldim,MRI_FLOAT,1);
  wbc->volmask->xsize = volres;
  wbc->volmask->ysize = volres;
  wbc->volmask->zsize = volres;
  wbc->volmask->tr    = 2000;

  c0 = nint(wbc->fvol->width/2.0);
  r0 = nint(wbc->fvol->height/2.0);
  s0 = nint(wbc->fvol->depth/2.0);
  // Column for "long" connections
  c2 = c0 + ceil(wbc->wbcsynth->distthresh/volres)+1;
  if(wbc->wbcsynth->ForceFail) c2 = c0 + floor(wbc->distthresh/volres)-1;
  printf("synth vol: %d %d %d, c2=%d, res=%lf, dim=%d\n",
	 c0,r0,s0,c2,volres,voldim);
  wbc->wbcsynth->c0 = c0;
  wbc->wbcsynth->r0 = r0;
  wbc->wbcsynth->s0 = s0;
  wbc->wbcsynth->c2 = c2;
  for(t=0; t < wbc->nframes; t++){
    // Set crs0 to have 2 short connections 
    MRIsetVoxVal(wbc->fvol,c0,r0,s0,t,wf[t]);
    MRIsetVoxVal(wbc->fvol,c0,r0+1,s0,t,wf[t]);
    MRIsetVoxVal(wbc->fvol,c0,r0-1,s0,t,wf[t]);
    // and 2 longs plus the longs from the surfaces
    MRIsetVoxVal(wbc->fvol,c2,r0,s0,t,wf[t]);
    MRIsetVoxVal(wbc->fvol,c2,r0+1,s0,t,wf[t]);
    if(wbc->wbcsynth->ForceFail) MRIsetVoxVal(wbc->fvol,c2,r0+1,s0+1,t,wf[t]);
  }
  MRIsetVoxVal(wbc->volmask,c0,r0,s0,0,1);
  MRIsetVoxVal(wbc->volmask,c0,r0+1,s0,0,1);
  MRIsetVoxVal(wbc->volmask,c0,r0-1,s0,0,1);
  MRIsetVoxVal(wbc->volmask,c2,r0,s0,0,1);
  MRIsetVoxVal(wbc->volmask,c2,r0+1,s0,0,1);
  MRIsetVoxVal(wbc->volmask,c2,r0+1,s0+1,0,1);   // mask but no signal, unless ForceFail
  MRIsetVoxVal(wbc->volmask,c2,r0+1,s0-1,0,1); // mask but no signal
  wbc->wbcsynth->nshortvol = 2; // number of short con to/from crs0
  wbc->wbcsynth->nlongvol = 2;  // number of long con to/from crs0

  wbc->wbcsynth->v0 = 0;
  for(hemi = 0; hemi < 2; hemi++){
    if(hemi == 0) hemistr = "lh";
    else          hemistr = "rh";

    sprintf(tmpstr,"%s/subjects/fsaverage/surf/%s.white",fsenv->FREESURFER_HOME,hemistr);
    surf = MRISread(tmpstr);
    sprintf(tmpstr,"%s/subjects/fsaverage/surf/%s.inflated",fsenv->FREESURFER_HOME,hemistr);
    surf2 = MRISread(tmpstr);
    //sprintf(tmpstr,"%s/subjects/fsaverage/label/%s.aparc.annot",fsenv->FREESURFER_HOME,hemistr);
    //err = MRISreadAnnotation(surf, tmpstr);

    func = MRIallocSequence(surf->nvertices,1,1,MRI_FLOAT,wbc->nframes);
    func->tr = 2000;
    mask = MRIallocSequence(surf->nvertices,1,1,MRI_FLOAT,1);

    // set the first vertex waveform
    MRIsetVoxVal(mask,0,0,0,0,1);
    for(t=0; t < wbc->nframes; t++) MRIsetVoxVal(func,0,0,0,t,wf[t]);
    
    x0 = surf->vertices[0].x;
    y0 = surf->vertices[0].y;
    z0 = surf->vertices[0].z;
    nshort = 0;
    nlong  = 0;
    for(vno = 0; vno < surf->nvertices; vno++){
      x = surf->vertices[vno].x;
      y = surf->vertices[vno].y;
      z = surf->vertices[vno].z;
      dx = (x-x0);
      dy = (y-y0);
      dz = (z-z0);
      d = sqrt(dx*dx + dy*dy + dz*dz);
      if(d < .9*wbc->wbcsynth->distthresh && nshort < wbc->wbcsynth->nshorttarg[hemi]){
	for(t=0; t < wbc->nframes; t++) MRIsetVoxVal(func,vno,0,0,t,wf[t]);
	MRIsetVoxVal(mask,vno,0,0,0,1);
	//printf("%s short %d %d\n",hemistr,nshort,vno);
	nshort ++;
      }
      if(d > 1.1*wbc->wbcsynth->distthresh && nlong < wbc->wbcsynth->nlongtarg[hemi]){
	for(t=0; t < wbc->nframes; t++) MRIsetVoxVal(func,vno,0,0,t,wf[t]);
	MRIsetVoxVal(mask,vno,0,0,0,1);
	//printf("%s long %d %d\n",hemistr,nlong,vno);
	nlong ++;
      }
    }
    wbc->wbcsynth->nshort[hemi] = nshort - 1; // remove self
    wbc->wbcsynth->nlong[hemi] = nlong;
    printf("%s short %d long %d\n",hemistr,nshort,nlong);
    if(hemi == 0) {
      wbc->lh = surf;
      wbc->lh2 = surf2;
      wbc->flh = func;
      wbc->lhmask = mask;
    }
    else {
      wbc->rh = surf;
      wbc->rh2 = surf2;
      wbc->frh = func;
      wbc->rhmask = mask;
    }
  } // hemi


  return(wbc);
}

int WBCtestCheck(WBC *wbc)
{
  int n,nexp,v0,c0,r0,s0,err;

  v0 = wbc->wbcsynth->v0;
  c0 = wbc->wbcsynth->c0;
  r0 = wbc->wbcsynth->r0;
  s0 = wbc->wbcsynth->s0;

  printf("WBCtestCheck(): ------------\n");
  printf("Vol short=%d, long=%d\n",wbc->wbcsynth->nshortvol,wbc->wbcsynth->nlongvol);
  printf("LH short=%d, long=%d\n",wbc->wbcsynth->nshort[0],wbc->wbcsynth->nlong[0]);
  printf("RH short=%d, long=%d\n",wbc->wbcsynth->nshort[1],wbc->wbcsynth->nlong[1]);
  printf("DistThresh: wbc = %lf, synth = %lf\n",wbc->distthresh,wbc->wbcsynth->distthresh);
  printf("ForceFail: %d\n",wbc->wbcsynth->ForceFail);

  err = 0;

  // expected number of cons from lhv0, rhv0, or crs0 to/from everyone else
  nexp = wbc->wbcsynth->nshortvol + wbc->wbcsynth->nlongvol +
    wbc->wbcsynth->nshort[0] + wbc->wbcsynth->nlong[0] +
    wbc->wbcsynth->nshort[1] + wbc->wbcsynth->nlong[1];
  // The above numbers exclude self from the shorts. There are 3 "selves"
  // so add 3, but then subtract one to account for the self for which
  // the comutation is being done, so +2.
  nexp += 2;
  printf("Expecting %d total connections\n",nexp);
  n = nint(MRIgetVoxVal(wbc->volcon,c0,r0,s0,0)*(wbc->ntot-1));
  printf("  Vol %d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");
  n = nint(MRIgetVoxVal(wbc->lhcon,v0,0,0,0)*(wbc->ntot-1));
  printf("  lh %d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");
  n = nint(MRIgetVoxVal(wbc->rhcon,v0,0,0,0)*(wbc->ntot-1));
  printf("  rh %d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");

  nexp = wbc->wbcsynth->nshortvol;
  printf("Volume short: expecting %3d  ",nexp);
  n = nint(MRIgetVoxVal(wbc->volconS,c0,r0,s0,0)*(wbc->ntot-1));
  printf("found %4d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");

  nexp = wbc->wbcsynth->nlongvol +
    wbc->wbcsynth->nshort[0] + wbc->wbcsynth->nlong[0] +
    wbc->wbcsynth->nshort[1] + wbc->wbcsynth->nlong[1];
  nexp += 2;
  printf("Volume long:  expecting %3d  ",nexp);
  n = nint(MRIgetVoxVal(wbc->volconL,c0,r0,s0,0)*(wbc->ntot-1));
  printf("found %4d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");

  nexp = wbc->wbcsynth->nshort[0];
  printf("lh short:       expecting %3d  ",nexp);
  n = nint(MRIgetVoxVal(wbc->lhconS,v0,0,0,0)*(wbc->ntot-1));
  printf("found %4d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");

  nexp = wbc->wbcsynth->nshortvol + wbc->wbcsynth->nlongvol +
    wbc->wbcsynth->nlong[0] +
    wbc->wbcsynth->nshort[1] + wbc->wbcsynth->nlong[1];
  nexp += 2;
  printf("lh long:        expecting %3d  ",nexp);
  n = nint(MRIgetVoxVal(wbc->lhconL,v0,0,0,0)*(wbc->ntot-1));
  printf("found %4d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");

  nexp = wbc->wbcsynth->nshort[1];
  printf("rh short:       expecting %3d  ",nexp);
  n = nint(MRIgetVoxVal(wbc->rhconS,v0,0,0,0)*(wbc->ntot-1));
  printf("found %4d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");

  nexp = wbc->wbcsynth->nshortvol + wbc->wbcsynth->nlongvol +
    wbc->wbcsynth->nshort[0] + wbc->wbcsynth->nlong[0] +
    wbc->wbcsynth->nlong[1];
  nexp += 2;
  printf("rh long:        expecting %3d  ",nexp);
  n = nint(MRIgetVoxVal(wbc->rhconL,v0,0,0,0)*(wbc->ntot-1));
  printf("found %4d ",n);
  if(n != nexp){
    printf(" ERROR\n");
    err++;
  } else  printf(" PASS\n");

  if(err){
    printf("Overall FAILED with %d errors\n",err);
    if(wbc->wbcsynth->ForceFail){
      printf(" ... but, a failure was forced, so this is a good thing\n");
      if(err != 9){
	printf(" ... however, expected 9 errors but got %d\n",err);
	return(err);
      }
      return(0);
    }
    return(err);
  }

  printf("Overall PASS\n");
  return(0);
}
