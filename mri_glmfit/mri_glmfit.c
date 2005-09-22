// mri_glmfit.c

// Incorporate glm into glmmri
// Add support for fsgdf
// Check to make sure no two contrast names are the same
// Check to make sure no two contrast mtxs are the same
// Save config in output dir, and Xg and Cs
// Links to source data
// PCA
// Permute X
// Leave out one
// cpmf
// Cleanup

// p-to-z
// Rewrite MatrixReadTxt to ignore # and % and empty lines
// Auto-det/read matlab4 matrices
// Profiling

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "matrix.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"
#include "gsl/gsl_cdf.h"
#include "pdf.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_glmfit.c,v 1.17 2005/09/22 23:12:38 greve Exp $";
char *Progname = NULL;

char *yFile = NULL, *XFile=NULL, *betaFile=NULL, *rvarFile=NULL;
char *yhatFile=NULL, *eresFile=NULL, *wFile=NULL, *maskFile;
char *condFile=NULL;
char *GLMDir=NULL;
char *pvrFiles[50];
int synth = 0;
int yhatSave=0;
int condSave=0;

MRI *mritmp=NULL, *sig=NULL;

int debug = 0, checkoptsonly = 0;
char tmpstr[2000];

int nContrasts=0;
char *CFile[100];
int err,c,r,s;
int SynthSeed = -1;
double Xcond;

int npvr=0;
MRIGLM *mriglm;

char *voxdumpdir=NULL;
int voxdump[3];

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs,n;
  struct utsname uts;
  char *cmdline, cwd[2000];

  /* rkt: check for and handle version tag */
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
  mriglm = (MRIGLM *) calloc(sizeof(MRIGLM),1);

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  if(checkoptsonly) return(0);

  printf("\n");
  printf("%s\n",vcid);
  printf("cwd %s\n",cwd);
  printf("cmdline %s\n",cmdline);
  printf("sysname  %s\n",uts.sysname);
  printf("hostname %s\n",uts.nodename);
  printf("machine  %s\n",uts.machine);

  if(SynthSeed < 0) SynthSeed = PDFtodSeed();
  printf("SynthSeed = %d\n",SynthSeed);

  dump_options(stdout);

  if(GLMDir != NULL){
    printf("Creating output directory %s\n",GLMDir);
    err = mkdir(GLMDir,(mode_t)-1);
    if(err != 0 && errno != EEXIST){
      printf("ERROR: creating directory %s\n",GLMDir);
      perror(NULL);    
      return(1);
    }
  }

  mriglm->npvr = npvr;
  mriglm->ncontrasts = nContrasts;
  mriglm->yhatsave = yhatSave;
  mriglm->condsave = condSave;

  //-----------------------------------------------------
  mriglm->Xg = MatrixReadTxt(XFile, NULL);
  if(mriglm->Xg==NULL){
    printf("ERROR: loading X %s\n",XFile);
    exit(1);
  }

  // Total number of columns in X, including PVRs
  mriglm->Xcols = mriglm->Xg->cols + npvr;

  // Compute and check DOF -------------------------------------
  mriglm->DOF = mriglm->Xg->rows - mriglm->Xcols;
  if(mriglm->DOF < 1){
    printf("ERROR: DOF = %g\n",mriglm->DOF);
    exit(1);
  }

  // Check the condition of the global matrix -----------------
  Xcond = MatrixNSConditionNumber(mriglm->Xg);
  printf("Matrix condition is %g\n",Xcond);
  if(Xcond > 10000){
    printf("ERROR: matrix is ill-conditioned or badly scaled, condno = %g\n",Xcond);
    exit(1);
  }

  // Load the contrast matrices ---------------------------------
  if(nContrasts > 0){
    for(n=0; n < nContrasts; n++){
      mriglm->C[n] = MatrixReadTxt(CFile[n], NULL);
      if(mriglm->C[n] == NULL){
	printf("ERROR: loading C %s\n",CFile[n]);
	exit(1);
      }
      if(mriglm->C[n]->cols != mriglm->Xcols){
	printf("ERROR: dimension mismatch between X and contrast %s",CFile[n]);
	printf("       X has %d cols, C has %d cols\n",
	       mriglm->Xcols,mriglm->C[n]->cols);
	exit(1);
      }
      // Get its name
      mriglm->cname[n] = fio_basename(CFile[n],".mat");
      // Should check to make sure no two are the same
    }
  }

  // Load or synthesize input  -----------------------------------------------
  if(synth == 0){
    printf("Loading y from %s\n",yFile);
    mriglm->y = MRIread(yFile);
    if(mriglm->y == NULL){
      printf("ERROR: loading y %s\n",yFile);
      exit(1);
    }
  }
  else{
    mritmp = MRIreadHeader(yFile,MRI_VOLUME_TYPE_UNKNOWN);
    if(mritmp == NULL){
      printf("ERROR: reading header for y %s\n",yFile);
      exit(1);
    }
    printf("Synthesizing y with white noise\n");
    mriglm->y =  MRIrandn(mritmp->width, mritmp->height, mritmp->depth, 
			 mritmp->nframes,0,1,NULL);
    MRIfree(&mritmp);
  }

  // Check some dimensions -----------------------------------
  if(mriglm->y->nframes != mriglm->Xg->rows){
    printf("ERROR: dimension mismatch between y and X.\n");
    printf("  y has %d inputs, X has %d rows.\n",
	   mriglm->y->nframes,mriglm->Xg->rows);
    exit(1);
  }

  // Load the mask file ----------------------------------
  if(maskFile != NULL){
    mriglm->mask = MRIread(maskFile);
    if(mriglm->mask  == NULL){
      printf("ERROR: reading mask file %s\n",maskFile);
      exit(1);
    }
  }
  else mriglm->mask = NULL;

  // Load the weight file ----------------------------------
  if(wFile != NULL){
    mriglm->w = MRIread(wFile);
    if(mriglm->w  == NULL){
      printf("ERROR: reading weight file %s\n",wFile);
      exit(1);
    }
    MRInormWeights(mriglm->w, 1, 1, mriglm->mask, mriglm->w);
    MRIwrite(mriglm->w,"wn.mgh");
  }
  else mriglm->w = NULL;

  // Load PVRs -----------------------------------
  if(mriglm->npvr > 0){
    for(n=0; n < mriglm->npvr; n++){
      mriglm->pvr[n] = MRIread(pvrFiles[n]);
      if(mriglm->pvr[n] == NULL) exit(1);
    }
  }

  if(voxdumpdir != NULL){
    printf("Dumping voxel %d %d %d to %s\n",
	   voxdump[0],voxdump[1],voxdump[2],voxdumpdir);
    printf("But not really\n");
  }

  printf("Starting fit\n");
  MRIglmFit(mriglm);

  printf("Writing results\n");
  MRIwrite(mriglm->beta,betaFile);
  MRIwrite(mriglm->rvar,rvarFile);
  if(mriglm->yhatsave) MRIwrite(mriglm->yhat,yhatFile);
  if(mriglm->condsave) MRIwrite(mriglm->cond,condFile);
  if(eresFile) MRIwrite(mriglm->eres,eresFile);    
  
  for(n=0; n < mriglm->ncontrasts; n++){
    printf("  %s\n",mriglm->cname[n]);
    
    // Create output directory for contrast
    sprintf(tmpstr,"%s/%s",GLMDir,mriglm->cname[n]);
    mkdir(tmpstr,(mode_t)-1);
    
    // Save gamma and F
    sprintf(tmpstr,"%s/%s/gamma.mgh",GLMDir,mriglm->cname[n]);
    MRIwrite(mriglm->gamma[n],tmpstr);
    sprintf(tmpstr,"%s/%s/F.mgh",GLMDir,mriglm->cname[n]);
    MRIwrite(mriglm->F[n],tmpstr);
    
    // Compute and Save p-values
    sig = fMRIsigF(mriglm->F[n], mriglm->DOF, mriglm->C[n]->rows, NULL);
    sprintf(tmpstr,"%s/%s/p.mgh",GLMDir,mriglm->cname[n]);
    MRIwrite(sig,tmpstr);
    
    // Compute and Save -log10 p-values
    MRIlog10(sig,sig,1);
    sprintf(tmpstr,"%s/%s/sig.mgh",GLMDir,mriglm->cname[n]);
    MRIwrite(sig,tmpstr);
    
    MRIfree(&sig);
  }
  
  printf("mri_glmfit done\n");
  return(0); 
  exit(0);

}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused, msec, niters;
  char **pargv, *option ;
  double rvartmp;

  if(argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--synth"))   synth = 1;
    else if (!strcasecmp(option, "--yhatsave")) yhatSave = 1;
    else if (!strcasecmp(option, "--condsave")) condSave = 1;

    else if (!strcasecmp(option, "--seed")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&SynthSeed);
      srand48(SynthSeed);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--voxdump")){
      if(nargc < 4) CMDargNErr(option,4);
      voxdumpdir = pargv[0];
      sscanf(pargv[1],"%d",&voxdump[0]);
      sscanf(pargv[2],"%d",&voxdump[1]);
      sscanf(pargv[3],"%d",&voxdump[2]);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--profile")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&niters);
      if(SynthSeed < 0) SynthSeed = PDFtodSeed();
      srand48(SynthSeed);
      printf("Starting GLM profile over %d iterations. Seed=%d\n",niters,SynthSeed);
      msec = GLMprofile(200, 20, 5, niters);
      nargsused = 1;
      exit(0);
    }
    else if (!strcasecmp(option, "--resynthtest")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&niters);
      if(SynthSeed < 0) SynthSeed = PDFtodSeed();
      srand48(SynthSeed);
      printf("Starting GLM resynth test over %d iterations. Seed=%d\n",niters,SynthSeed);
      err = GLMresynthTest(niters, &rvartmp);
      if(err){
	printf("Failed. rvar = %g\n",rvartmp);
	exit(1);
      }
      printf("Passed. rvarmax = %g\n",rvartmp);
      exit(0);
      nargsused = 1;
    }
    else if (!strcmp(option, "--y")){
      if(nargc < 1) CMDargNErr(option,1);
      yFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--mask")){
      if(nargc < 1) CMDargNErr(option,1);
      maskFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--w")){
      if(nargc < 1) CMDargNErr(option,1);
      wFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--X")){
      if(nargc < 1) CMDargNErr(option,1);
      XFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--pvr")){
      if(nargc < 1) CMDargNErr(option,1);
      pvrFiles[npvr] = pargv[0];
      npvr++;
      nargsused = 1;
    }
    else if (!strcmp(option, "--glmdir")){
      if(nargc < 1) CMDargNErr(option,1);
      GLMDir = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--beta")){
      if(nargc < 1) CMDargNErr(option,1);
      betaFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--rvar")){
      if(nargc < 1) CMDargNErr(option,1);
      rvarFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--yhat")){
      if(nargc < 1) CMDargNErr(option,1);
      yhatFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--eres")){
      if(nargc < 1) CMDargNErr(option,1);
      eresFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--C")){
      if(nargc < 1) CMDargNErr(option,1);
      CFile[nContrasts] = pargv[0];
      nContrasts++;
      nargsused = 1;
    }
    else{
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
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --y input volume \n");
  printf("   --X design matrix file\n");
  printf("   --pvr pvr1 <--prv pvr2 ...>\n");
  printf("   --w weight volume\n");
  printf("   --mask mask volume\n");
  printf("   --C contrast1.mat <--C contrast2.mat ...>\n");
  printf("\n");
  printf("   --glmdir dir : save outputs to dir\n");
  printf("\n");
  printf("   --saveyhat \n");
  printf("   --savecond \n");
  printf("\n");
  printf("   --synth : replace input with gaussian \n");
  printf("   --seed seed \n");
  printf("\n");
  printf("   --resynthtest niters \n");
  printf("   --profile     niters \n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(yFile == NULL){
    printf("ERROR: must specify input y file\n");
    exit(1);
  }
  if(XFile == NULL){
    printf("ERROR: must specify an input X file\n");
    exit(1);
  }
  if(GLMDir != NULL){
    sprintf(tmpstr,"%s/beta.mgh",GLMDir);
    betaFile = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/rvar.mgh",GLMDir);
    rvarFile = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/eres.mgh",GLMDir);
    eresFile = strcpyalloc(tmpstr);
    if(yhatSave){
      sprintf(tmpstr,"%s/yhat.mgh",GLMDir);
      yhatFile = strcpyalloc(tmpstr);
    }
    if(condSave){
      sprintf(tmpstr,"%s/cond.mgh",GLMDir);
      condFile = strcpyalloc(tmpstr);
    }
  }
  else{
    if(betaFile == NULL){
      printf("ERROR: must specify an output dir or beta file\n");
      exit(1);
    }
    if(rvarFile == NULL){
      printf("ERROR: must specify an output dir or rvar file\n");
      exit(1);
    }
  }

  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"y    %s\n",yFile);
  fprintf(fp,"X    %s\n",XFile);
  fprintf(fp,"beta %s\n",betaFile);
  fprintf(fp,"rvar %s\n",rvarFile);
  if(maskFile) fprintf(fp,"mask %s\n",maskFile);
  if(yhatFile) fprintf(fp,"yhat %s\n",yhatFile);
  if(eresFile) fprintf(fp,"eres %s\n",eresFile);
  if(condFile) fprintf(fp,"cond %s\n",condFile);

  return;
}


