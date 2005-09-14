// mri_glmfit.c

// Cleanup
// Weighting
// Save config in output dir
// Links to source data

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

#undef X

typedef struct{
  MRI *y;            // Input data
  MATRIX *X;         // Global regressors
  int npvr;          // Number of per-voxel regressors
  MRI *pvr[50];      // Per-voxel regressors (local)
  float DOF;         // DOF

  MRI *iXtX;         // inv(X'X), where X is global and pv
  MRI *beta;         // beta = inv(X'X)*X'*y
  MRI *yhat;         // yhat = X*beta
  MRI *eres;         // eres = y - yhat
  MRI *rvar;         // rvar = sum(eres.^2)/DOF;

  int ncontrasts;    // Number of contrasts
  char *cname[100];  // Contrast names
  MATRIX *C[100];    // Contrast matrices
  MATRIX *Ct[100];   // transposes Contrast matrices
  MRI *gamma[100];   // gamma = C*beta
  MRI *F[100];       // F = gamma'*iXtX*gamma/(rvar*J)
  MRI *sig[100];     // sig = significance of the F
} GLMPV;

int MRIglmpvFit(GLMPV *glmpv);
MATRIX *MRItoMatrix(MRI *mri, int c, int r, int s, 
		    int Mrows, int Mcols, MATRIX *M);
MATRIX *MRItoSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M);
int MRIfromMatrix(MRI *mri, int c, int r, int s, MATRIX *M);
int MRIfromSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_glmfit.c,v 1.9 2005/09/14 18:25:08 greve Exp $";
char *Progname = NULL;

char *yFile = NULL, *XFile=NULL, *betaFile=NULL, *rvarFile=NULL;
char *yhatFile=NULL, *eresFile=NULL;
char *GLMDir=NULL;
char *pvrFiles[50];
int synth = 0;
int yhatSave=0;

MRI *y, *beta, *rvar, *yhat, *eres, *mritmp;
MRI *gam, *rstat, *sig;

int debug = 0, checkoptsonly = 0;
char tmpstr[2000];

MATRIX *X; /* design matrix */
MATRIX *H=NULL, *Xt=NULL, *XtX=NULL, *iXtX=NULL, *Q=NULL, *R=NULL;
MATRIX *C, *M;
int nContrasts=0;
char *CFile[100];
char *ContrastName;
float DOF;
int err,c,r,s;

int npvr=0;
GLMPV *glmpv;

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

  //-----------------------------------------------------
  printf("Loading design matrix from %s\n",XFile);
  X = MatrixReadTxt(XFile, NULL);
  if(X==NULL){
    printf("ERROR: loading X %s\n",XFile);
    exit(1);
  }

  //-----------------------------------------------------
  if(nContrasts > 0){
    for(n=0; n < nContrasts; n++){
      C = MatrixReadTxt(CFile[n], NULL);
      if(C==NULL){
	printf("ERROR: loading C %s\n",CFile[n]);
	exit(1);
      }
      if(C->cols != X->cols + npvr){
	printf("ERROR: dimension mismatch between X and contrast %s",CFile[n]);
	printf("       X has %d cols, C has %d cols\n",X->cols,C->cols);
	exit(1);
      }
      MatrixFree(&C);
    }
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

  if(npvr){
    glmpv = (GLMPV *) calloc(sizeof(GLMPV),1);
    glmpv->npvr = npvr;
    glmpv->ncontrasts = nContrasts;
    glmpv->y = y;
    glmpv->X = X;
    for(n=0; n < npvr; n++){
      glmpv->pvr[n] = MRIread(pvrFiles[n]);
      if(glmpv->pvr[n] == NULL) exit(1);
    }
    for(n=0; n < nContrasts; n++){
      glmpv->C[n] = MatrixReadTxt(CFile[n], NULL);
      if(glmpv->C[n]==NULL){
	printf("ERROR: loading C %s\n",CFile[n]);
	exit(1);
      }
    }
    printf("Starting Fit\n");
    MRIglmpvFit(glmpv);
    MRIwrite(glmpv->beta,betaFile);
    MRIwrite(glmpv->rvar,rvarFile);
    if(yhatFile) MRIwrite(glmpv->yhat,yhatFile);
    if(eresFile) MRIwrite(glmpv->eres,eresFile);    

    for(n=0; n < nContrasts; n++){
      ContrastName = fio_basename(CFile[n],".mat");
      sprintf(tmpstr,"%s/%s",GLMDir,ContrastName);
      mkdir(tmpstr,(mode_t)-1);
      printf("%s\n",ContrastName);

      sprintf(tmpstr,"%s/%s/gamma.mgh",GLMDir,ContrastName);
      MRIwrite(glmpv->gamma[n],tmpstr);

      sprintf(tmpstr,"%s/%s/stat.mgh",GLMDir,ContrastName);
      MRIwrite(glmpv->F[n],tmpstr);

      //sig = fMRIsigT(rstat, DOF, NULL);

      //sig = fMRIsigF(glmpv->F[n], glmpv->DOF, glmpv->C[n]->rows, NULL);
      sig = fMRIsigF(glmpv->F[n], glmpv->DOF, glmpv->C[n]->rows, NULL);
      printf("DOF = %g, J = %d, F = %g, p = %g\n",
	     glmpv->DOF, glmpv->C[n]->rows,
	     MRIgetVoxVal(glmpv->F[n],0,0,0,0),
	     MRIgetVoxVal(sig,0,0,0,0));


      sprintf(tmpstr,"%s/%s/p.mgh",GLMDir,ContrastName);
      MRIwrite(sig,tmpstr);

      MRIlog10(sig,sig,1);
      sprintf(tmpstr,"%s/%s/sig.mgh",GLMDir,ContrastName);
      MRIwrite(sig,tmpstr);

      MRIfree(&sig);
    }

    printf("mri_glmfit done\n");
    return(0); exit(0);

  }

  //------------------------------------------------------
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

  if(nContrasts == 0){
    printf("mri_glmfit done\n");
    return(0);
  }

  //---------------------------------------------------
  for(n=0; n < nContrasts; n++){
    C = MatrixReadTxt(CFile[n], NULL);
    if(C==NULL){
      printf("ERROR: loading C %s\n",CFile[n]);
      exit(1);
    }
    ContrastName = fio_basename(CFile[n], ".mat");
    sprintf(tmpstr,"%s/%s",GLMDir,ContrastName);
    mkdir(tmpstr,(mode_t)-1);
    
    printf("%s\n",ContrastName);

    gam = fMRImatrixMultiply(beta,C,NULL);
    sprintf(tmpstr,"%s/%s/gamma.mgh",GLMDir,ContrastName);
    MRIwrite(gam,tmpstr);

    if(C->rows == 1) rstat = fMRIcomputeT(gam, X, C, rvar, NULL);
    else             rstat = fMRIcomputeF(gam, X, C, rvar, NULL);
    sprintf(tmpstr,"%s/%s/stat.mgh",GLMDir,ContrastName);
    MRIwrite(rstat,tmpstr);

    if(C->rows == 1) sig = fMRIsigT(rstat, DOF, NULL);
    else             sig = fMRIsigF(rstat, DOF, C->rows, NULL);
    MRIlog10(sig,sig,1);
    sprintf(tmpstr,"%s/%s/sig.mgh",GLMDir,ContrastName);
    MRIwrite(sig,tmpstr);

    //-- Should do a PMF
    
    MatrixFree(&C);
    MRIfree(&gam);
    MRIfree(&rstat);
    MRIfree(&sig);

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
      if(nargc < 1) CMDargNErr(option,1);
      yFile = pargv[0];
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
  printf("   --glmdir dir : save outputs to dir\n");
  printf("   --C contrast1.mat <--C contrast2.mat ...>\n");
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
MATRIX *MRItoMatrix(MRI *mri, int c, int r, int s, 
		    int Mrows, int Mcols, MATRIX *M)
{
  int mr, mc, f;

  if(M==NULL) M = MatrixAlloc(Mrows,Mcols,MATRIX_REAL);
  else{
    if(M->rows != Mrows || M->cols != Mcols){
      printf("ERROR: Matrix dim mismatch\n");
    }
  }

  if(mri->nframes != Mrows*Mcols){
    printf("ERROR: MRItoMatrix: MRI frames = %d, does not equal\n",
	   mri->nframes);
    printf("       matrix dim = %dx%d = %d",Mrows,Mcols,Mrows*Mcols);
    return(NULL);
  }

  f = 0;
  for(mr=1; mr <= Mrows; mr++){
    for(mc=1; mc <= Mcols; mc++){
      M->rptr[mr][mc] = MRIgetVoxVal(mri,c,r,s,f);
      f++;
    }
  }
  return(M);
}

/*---------------------------------------------------------------*/
MATRIX *MRItoSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M)
{
  int mr, mc, f, Msize, nframesexp;

  if(M==NULL){
    Msize = (int)(round( (sqrt(8.0*mri->nframes + 1.0) - 1.0 )/2.0 ));
    printf("Msize = %d\n",Msize);
    M = MatrixAlloc(Msize,Msize,MATRIX_REAL);
  }

  nframesexp = M->rows*(M->rows+1)/2;
  if(mri->nframes != nframesexp){
    printf("ERROR: MRItoSymMatrix: MRI frames = %d, does not support sym\n",
	   mri->nframes);
    return(NULL);
  }

  f = 0;
  for(mr=1; mr <= M->rows; mr++){
    for(mc=mr; mc <= M->cols; mc++){
      M->rptr[mr][mc] = MRIgetVoxVal(mri,c,r,s,f);
      M->rptr[mc][mr] = MRIgetVoxVal(mri,c,r,s,f);
      f++;
    }
  }
  return(M);
}

/*---------------------------------------------------------------*/
int MRIfromMatrix(MRI *mri, int c, int r, int s, MATRIX *M)
{
  int mr,mc,f;

  if(mri->nframes != M->rows*M->cols){
    printf("ERROR: MRIfromMatrix: MRI frames = %d, does not equal\n",
	   mri->nframes);
    printf("       matrix dim = %dx%d = %d",M->rows,M->cols,M->rows*M->cols);
    return(1);
  }

  f = 0;
  for(mr=1; mr <= M->rows; mr++){
    for(mc=1; mc <= M->cols; mc++){
      MRIsetVoxVal(mri,c,r,s,f,M->rptr[mr][mc]);
      f++;
    }
  }
  return(0);
}

/*---------------------------------------------------------------*/
int MRIfromSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M)
{
  int mr,mc,f, nframesexp;

  nframesexp = M->rows*(M->rows+1)/2;
  if(mri->nframes != nframesexp){
    printf("ERROR: MRIfromSumMatrix: MRI frames = %d, does not equal\n",
	   mri->nframes);
    printf("       matrix dim = %dx%d = %d",M->rows,M->cols,M->rows*M->cols);
    return(1);
  }

  f = 0;
  for(mr=1; mr <= M->rows; mr++){
    for(mc=mr; mc <= M->cols; mc++){
      MRIsetVoxVal(mri,c,r,s,f,M->rptr[mr][mc]);
      f++;
    }
  }
  return(0);
}

/*---------------------------------------------------------------*/
int MRIglmpvFit(GLMPV *glmpv)
{
  int c,r,s,n,f,nc,nr,ns,n_ill_cond;
  MATRIX *X, *Xt=NULL, *XtX=NULL, *iXtX=NULL, *y, *Xty=NULL;
  MATRIX *beta=NULL,*yhat=NULL,*eres=NULL, *Mtmp;
  MATRIX *C[50],*Ct[50];
  MATRIX *gam=NULL, *gamt=NULL, *CiXtX=NULL,*CiXtXCt=NULL;
  MATRIX *iCiXtXCt=NULL,*gtiCiXtXCt=NULL, *gtiCiXtXCtg=NULL;
  double rvar,F;

  glmpv->DOF = glmpv->X->rows - (glmpv->X->cols + glmpv->npvr);

  X = MatrixAlloc(glmpv->X->rows, glmpv->X->cols + glmpv->npvr, MATRIX_REAL);
  y = MatrixAlloc(glmpv->X->rows, 1, MATRIX_REAL);

  nc = glmpv->y->width;
  nr = glmpv->y->height;
  ns = glmpv->y->depth;
  glmpv->beta = MRIallocSequence(nc, nr, ns, MRI_FLOAT, X->cols) ;
  glmpv->yhat = MRIallocSequence(nc, nr, ns, MRI_FLOAT, glmpv->y->nframes) ;
  glmpv->eres = MRIallocSequence(nc, nr, ns, MRI_FLOAT, glmpv->y->nframes) ;
  glmpv->rvar = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);

  for(n = 0; n < glmpv->ncontrasts; n++){
    glmpv->gamma[n] = MRIallocSequence(nc,nr,ns,MRI_FLOAT, glmpv->C[n]->rows);
    glmpv->F[n] = MRIallocSequence(nc, nr, ns,MRI_FLOAT, 1);
    C[n] = glmpv->C[n];
    Ct[n] = MatrixTranspose(C[n],NULL);
  }

  // pre-load X
  for(f = 1; f <= X->rows; f++){
    for(n = 1; n <= glmpv->X->cols; n++){
      X->rptr[f][n] = glmpv->X->rptr[f][n];
    }
  }

  //--------------------------------------------
  n_ill_cond = 0;
  for(c=0; c < nc; c++){
    printf("%d \n",c);
    for(r=0; r < nr; r++){
      for(s=0; s< ns; s++){

	// Load y and the per-vox reg --------------------------
	for(f = 1; f <= X->rows; f++){
	  y->rptr[f][1] = MRIgetVoxVal(glmpv->y,c,r,s,f-1);
	  for(n = 1; n <= glmpv->npvr; n++){
	    X->rptr[f][n+glmpv->X->cols] = 
	      MRIgetVoxVal(glmpv->pvr[n-1],c,r,s,f-1);
	  }
	}

	Xt   = MatrixTranspose(X,Xt);
	XtX  = MatrixMultiply(Xt,X,XtX);
	Mtmp = MatrixInverse(XtX,iXtX);
	if(Mtmp == NULL){
	  MatrixPrint(stdout,X);
	  exit(1);
	  n_ill_cond++;
	  continue;
	}
	iXtX = Mtmp;

	Xty  = MatrixMultiply(Xt,y,Xty);
	beta = MatrixMultiply(iXtX,Xty,beta);
	yhat = MatrixMultiply(X,beta,yhat);
	eres = MatrixSubtract(y, yhat, eres);
	
	rvar = 0;
	for(f = 1; f <= eres->rows; f++)
	  rvar += (eres->rptr[f][1] * eres->rptr[f][1]);
	rvar /= glmpv->DOF;
	MRIsetVoxVal(glmpv->rvar,c,r,s,0,rvar);

	MRIfromMatrix(glmpv->beta, c, r, s, beta);
	MRIfromMatrix(glmpv->yhat, c, r, s, yhat);
	MRIfromMatrix(glmpv->eres, c, r, s, eres);

	for(n = 0; n < glmpv->ncontrasts; n++){
	  gam         = MatrixMultiply(C[n],beta,gam);
	  gamt        = MatrixTranspose(gam,gamt);
	  CiXtX       = MatrixMultiply(C[n],iXtX,CiXtX);
	  CiXtXCt     = MatrixMultiply(CiXtX,Ct[n],CiXtXCt);
	  iCiXtXCt    = MatrixInverse(CiXtXCt,iCiXtXCt);
	  gtiCiXtXCt  = MatrixMultiply(gamt,iCiXtXCt,gtiCiXtXCt);
	  gtiCiXtXCtg = MatrixMultiply(gtiCiXtXCt,gam,gtiCiXtXCtg);
	  F           = gtiCiXtXCtg->rptr[1][1]/(rvar/C[n]->rows);
	  MRIfromMatrix(glmpv->gamma[n], c, r, s, gam);
	  MRIsetVoxVal(glmpv->F[n],c,r,s,0,F);
	}

      }
    }
  }

  printf("n_ill_cond = %d\n",n_ill_cond);
  return(0);
}




#if 0
  printf("Packing and unpacking\n");
  for(c=0; c<y->width; c++){
    for(r=0; r<y->height; r++){
      for(s=0; s<y->depth; s++){
	M = MRItoSymMatrix(y,c,r,s,M);
	MRIfromSymMatrix(y,c,r,s,M);
      }
    }
  }
  printf("Done Packing and unpacking\n");    
  MRIwrite(y,"y2.mgh");
#endif
