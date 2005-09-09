// mri_glmfit.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

#undef X

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_glmfit.c,v 1.1 2005/09/09 22:37:15 greve Exp $";
char *Progname = NULL;

char *yFile = NULL, *XFile=NULL, *betaFile=NULL, *rvarFile=NULL;
char *yhatFile=NULL, *eresFile=NULL;
char *OutDir=NULL;
int synth = 0;
int yhatSave=0;

MRI *y, *beta, *rvar, *yhat, *eres, *mritmp;

int debug = 0, checkoptsonly = 0;
char tmpstr[2000];

MATRIX *X; /* design matrix */
MATRIX *H=NULL, *Xt=NULL, *XtX=NULL, *iXtX=NULL, *Q=NULL, *R=NULL;
float DOF;
int err;

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs;
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

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  if(checkoptsonly){
    return(0);
    exit(0);
  }

  printf("\n");
  printf("%s\n",vcid);
  printf("cwd %s\n",cwd);
  printf("cmdline %s\n",cmdline);
  printf("sysname  %s\n",uts.sysname);
  printf("hostname %s\n",uts.nodename);
  printf("machine  %s\n",uts.machine);
  dump_options(stdout);

  if(OutDir != NULL){
    printf("Creating output directory %s\n",OutDir);
    err = mkdir(OutDir,(mode_t)-1);
    if(err != 0 && errno != EEXIST){
      printf("ERROR: creating directory %s\n",OutDir);
      perror(NULL);    
      return(1);
    }
  }

  //-----------------------------------------------------
  printf("Loading design matrix from %s\n",XFile);
  X = MatrixReadTxt(XFile, NULL);
  if(X==NULL){
    printf("ERROR: loading X %s\n",XFile);
    exit(1);
  }

  //-----------------------------------------------------
  if(synth == 0){
    printf("Loading y from %s\n",yFile);
    y = MRIread(yFile);
    if(y == NULL){
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
    y =  MRIrandn(mritmp->width, mritmp->depth, mritmp->height, 
		  mritmp->nframes,0,1,NULL);
    MRIfree(&mritmp);
  }

  //-----------------------------------------------------
  if(y->nframes != X->rows){
    printf("ERROR: dimension mismatch between y and X.\n");
    printf("  y has %d inputs, X has %d rows.\n",y->nframes,X->rows);
    exit(1);
  }

  if(y->nframes <= X->cols){
    printf("ERROR: number of columns in X exceeds number of inputs\n");
    exit(1);
  }

  Xt = MatrixTranspose(X,NULL);
  XtX = MatrixMultiply(Xt,X,NULL);
  iXtX = MatrixInverse(XtX,NULL);
  if(iXtX==NULL){
    printf("ERROR: could not compute psuedo inverse of X\n");
    exit(1);
  }
  /* Q is the matrix that when multiplied by y gives beta */
  Q = MatrixMultiply(iXtX,Xt,NULL);
  /* H is the matrix that when multiplied by y gives the signal estimate */
  H = MatrixMultiply(X,Q,NULL);
  /* R is the matrix that when multiplied by y gives the residual error */
  R = MatrixSubtract(MatrixIdentity(y->nframes,NULL),H,NULL);
  DOF = MatrixTrace(R);
  printf("DOF = %f\n",DOF);

  //-----------------------------------------------------
  printf("Computing beta\n");
  beta = fMRImatrixMultiply(y, Q, NULL);
  MRIwrite(beta,betaFile);

  printf("Computing residuals\n");
  eres = fMRImatrixMultiply(y, R, NULL);
  if(eresFile) MRIwrite(eres,eresFile);

  printf("Computing residual variance\n");
  rvar = fMRIvariance(eres,DOF,0,NULL);
  MRIwrite(rvar,rvarFile);

  if(yhatFile){ 
    printf("Computing signal estimate\n");
    yhat = fMRImatrixMultiply(y, H, NULL);
    MRIwrite(yhat,yhatFile);
  }

  printf("mri_glmfit done\n");
  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

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

    else if (!strcmp(option, "--y")){
      if(nargc < 1) argnerr(option,1);
      yFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--X")){
      if(nargc < 1) argnerr(option,1);
      XFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--outdir")){
      if(nargc < 1) argnerr(option,1);
      OutDir = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--beta")){
      if(nargc < 1) argnerr(option,1);
      betaFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--rvar")){
      if(nargc < 1) argnerr(option,1);
      rvarFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--yhat")){
      if(nargc < 1) argnerr(option,1);
      yhatFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--eres")){
      if(nargc < 1) argnerr(option,1);
      eresFile = pargv[0];
      nargsused = 1;
    }
    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
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
  printf("   --outdir dir : save outputs to dir\n");
  printf("\n");
  printf("   --beta regression coeffient volume\n");
  printf("   --rvar residual variance\n");
  printf("   --yhat input estimate\n");
  printf("   --eres residual\n");
  printf("\n");
  printf("   --synth : replace input with gaussian \n");
  printf("\n");
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
static void argnerr(char *option, int n)
{
  if(n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
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
  if(OutDir != NULL){
    sprintf(tmpstr,"%s/beta.mgh",OutDir);
    betaFile = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/rvar.mgh",OutDir);
    rvarFile = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/eres.mgh",OutDir);
    eresFile = strcpyalloc(tmpstr);
    if(yhatSave){
      sprintf(tmpstr,"%s/yhat.mgh",OutDir);
      yhatFile = strcpyalloc(tmpstr);
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
  if(yhatFile) fprintf(fp,"yhat %s\n",yhatFile);
  if(eresFile) fprintf(fp,"eres %s\n",eresFile);

  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}


