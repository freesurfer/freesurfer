/*
  Name:    mris_glm.c
  Author:  Douglas N. Greve 
  email:   analysis-bugs@nmr.mgh.harvard.edu
  Date:    2/27/02
  Purpose: Computes glm inferences on the surface.
  $Id: mris_glm.c,v 1.7 2002/10/29 17:46:57 greve Exp $

Things to do:
  0. Documentation.
  1. More sophisticated derived variable.
  2. Input volume directly.
  3. Add ability to load in y directly
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "error.h"
#include "diag.h"
#include "proto.h"

#include "icosahedron.h"
#include "matrix.h"
#include "matfile.h"
#include "mri.h"
#include "MRIio_old.h"
#include "mri_identify.h"
#include "sig.h"
#include "fmriutils.h"

#ifdef X
#undef X
#endif

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  isflag(char *flag);
static int  nth_is_arg(int nargc, char **argv, int nth);
static int  singledash(char *flag);
static int  stringmatch(char *str1, char *str2);
static int checkfmt(char *fmt);
static int getfmtid(char *fname);
static int IsSurfFmt(char *fmt);
MRIS *MRISloadSurfSubject(char *subj, char *hemi, char *surfid, 
			  char *SUBJECTS_DIR);

int ReadAsciiMatrixNRows(char *desmtxfname);
int ReadAsciiMatrixSize(char *desmtxfname, int *pnrows, int *pncols);
int ReadDesignMatrix(char *desmtxfname);
MATRIX *ReadAsciiMatrix(char *asciimtxfname);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mris_glm.c,v 1.7 2002/10/29 17:46:57 greve Exp $";
char *Progname = NULL;

char *hemi        = NULL;
char *desmtxfname = NULL;
char *xmatfile = NULL;
int  nsmooth   = 0;
int  frame     = 0;    

char *surfmeasure    = NULL;

char *surfregid   = "sphere.reg";
int  ninputs = 0;
char *inputlist[1000];
char *inputfmt = NULL;
int   inputfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *subjectlistfile;
int  nsubjects = 0;
char *subjectlist[1000];
int  nregressors = 0;
MATRIX *X; /* design matrix */

char   *conmtxfname;
MATRIX *C; /* contrast vector */

char *betaid  = NULL;
char *betafmt = NULL;
int  betafmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *beta_in_id  = NULL;
char *beta_in_fmt = NULL;
int  beta_in_fmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *cesid  = NULL;
char *cesfmt = NULL;
int  cesfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *eresid  = NULL;
char *eresfmt = NULL;
int   eresfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *yid  = NULL;
char *yfmt = NULL;
int   yfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *yhatid  = NULL;
char *yhatfmt = NULL;
int   yhatfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *eresvarid  = NULL;
char *eresvarfmt = NULL;
int   eresvarfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *eresvar_in_id  = NULL;
char *eresvar_in_fmt = NULL;
int   eresvar_in_fmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *tid  = NULL;
char *tfmt = NULL;
int   tfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *sigid  = NULL;
char *sigfmt = NULL;
int   sigfmtid = MRI_VOLUME_TYPE_UNKNOWN;

int  IcoOrder     = 7;
float IcoRadius = 100.0;
char *regsurf     = "sphere.reg";
char *trgsubject = NULL;

int dof;
int debug = 0;
int synthseed = 0;

MRI_SURFACE *IcoSurf, *SurfReg;
MRI *SrcVals, *beta, *yhat, *eres, *eresvar, *ces, *t, *sig;
MATRIX *T, *Xt, *XtX, *iXtX, *Q, *R;
MRI *tmpmri, *tmpmri2, *SrcHits, *SrcDist, *TrgHits, *TrgDist;

float DOF;
char *SUBJECTS_DIR;

/*---------------------------------------------------------------*/
int main(int argc, char **argv)
{
  int vtx,nthsubj;
  char *subject;
  char *inputfname;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  dump_options(stdout);

  printf("%s\n",vcid);

  if(xmatfile != NULL) MatlabWrite(X,xmatfile,"X");

  /* X is the design matrix */
  Xt = MatrixTranspose(X,NULL);
  XtX = MatrixMultiply(Xt,X,NULL);
  iXtX = MatrixInverse(XtX,NULL);

  /* Q is the matrix that when multiplied by y gives beta */
  Q = MatrixMultiply(iXtX,Xt,NULL);

  /* T is the matrix that when multiplied by y gives the signal estimate */
  T = MatrixMultiply(X,Q,NULL);

  /* R is the matrix that when multiplied by y gives the residual error */
  R = MatrixSubtract(MatrixIdentity(nsubjects,NULL),T,NULL);
  DOF = MatrixTrace(R);

  printf("Design Matrix ------------------------------------\n");
  MatrixPrint(stdout,X);
  //printf("Q ------------------------------------\n");
  //MatrixPrint(stdout,Q);
  //printf("T ------------------------------------\n");
  //MatrixPrint(stdout,T);
  //printf("R ------------------------------------\n");
  //MatrixPrint(stdout,R);
  printf("Design Covariance Matrix ------------------------------------\n");
  MatrixPrint(stdout,XtX);
  if(C != NULL){
    printf("Contrast Matrix: -----------------------------\n");
    MatrixPrint(stdout,C);
  }
  printf("DOF = %g\n",DOF);
  if(DOF < 1){
    printf("ERROR: zero (or fewer) degrees of freedom\n");
    exit(1);
  }
  fflush(stdout);

  if(beta_in_id == NULL){
    /* Compute beta and reserr from data (real or synth) */

    if(stringmatch(trgsubject,"ico")){
      /* Use Icosahedron as target surface */
      printf("INFO: loading ico (order = %d)\n",IcoOrder);
      fflush(stdout);
      IcoSurf = ReadIcoByOrder(IcoOrder, IcoRadius);
    }
    else{
      /* Use target subject (still called IcoSurf) */
      IcoSurf = MRISloadSurfSubject(trgsubject,hemi,surfregid,SUBJECTS_DIR);
      if(IcoSurf == NULL){
	printf("ERROR: could not load registration surface\n");
	exit(1);
      }
    }

    /* Alloc enough data for all the subjects */
    SrcVals = MRIallocSequence(IcoSurf->nvertices, 1, 1,MRI_FLOAT,nsubjects);
    
    if(synthseed == 0){

      /* Load in all the subjects' data, one at a time */      
      for(nthsubj = 0; nthsubj < nsubjects; nthsubj++){
	
	subject = subjectlist[nthsubj];
	printf("nthsubj = %d/%d, %s ---------------------------------------\n",
	       nthsubj+1,nsubjects,subject);
	fflush(stdout);
	
	/* Read in the spherical registration surface */
	SurfReg = MRISloadSurfSubject(subject,hemi,surfregid,SUBJECTS_DIR);
	if(SurfReg == NULL){
	  printf("ERROR: could not load registration surface\n");
	  exit(1);
	}
	
	if(surfmeasure != NULL) inputfname = surfmeasure;
	else                    inputfname = inputlist[nthsubj];

	/* Read in the input for this subject. */
	/* If the input is a volume, then read in the registration file,
	   check the subject name, load white surface, load thickness
	   for projection, resample to the white surface, save in tmpmri */

	printf("  INFO: loading input %s \n",inputfname);fflush(stdout);
	if(stringmatch(inputfmt,"curv"))
	  tmpmri = MRISloadSurfVals(inputfname,"curv",NULL,subject,hemi,NULL);
	else if(stringmatch(inputfmt,"paint") || stringmatch(inputfmt,"w") ||
		stringmatch(inputfmt,"wfile")){
	  MRISreadValues(SurfReg,inputfname);
	  tmpmri = MRIcopyMRIS(NULL, SurfReg, 0, "val");
	}
	else{
	  if(inputfmt != NULL) tmpmri = MRIreadType(inputfname,inputfmtid);
	  else                 tmpmri = MRIread(inputfname);
	}
	if(tmpmri == NULL){
	  printf("ERROR: could not load %s\n",inputfname);
	  exit(1);
	}

	/* Extract frame (Future: compute derived variable) */
	if(tmpmri->nframes <= frame){
	  printf("ERROR: nframes (%d) <= frame (%d)\n",
		 tmpmri->nframes,frame);
	  exit(1);
	}
	if(tmpmri->nframes > 0){
	  /* extract frame */
	  tmpmri2 = MRIallocSequence(SurfReg->nvertices,1,1,MRI_FLOAT,1);
	  for(vtx=0; vtx < SurfReg->nvertices; vtx++)
	    MRIFseq_vox(tmpmri2,vtx,0,0,0) = MRIFseq_vox(tmpmri,vtx,0,0,frame);
	  MRIfree(&tmpmri);
	  tmpmri = tmpmri2;
	}
	
	/* Smooth on the native surface */
	if(nsmooth > 0)
	  MRISsmoothMRI(SurfReg, tmpmri, nsmooth, tmpmri);

	/*------- Resample to target subject -------------------*/
	if(!stringmatch(trgsubject,subject)){
	  printf("  INFO: resampling to %s\n",trgsubject);fflush(stdout);
	  tmpmri2 = surf2surf_nnfr(tmpmri, SurfReg, IcoSurf,
				   &SrcHits,&SrcDist,&TrgHits,&TrgDist,1,1);
	  if(tmpmri2 == NULL){
	    printf("ERROR: could not resample to %s.\n",trgsubject);
	    exit(1);
	  }
	  MRIfree(&SrcHits);
	  MRIfree(&SrcDist);
	  MRIfree(&TrgHits);
	  MRIfree(&TrgDist);
	  MRIfree(&tmpmri);
	  tmpmri = tmpmri2;
	}
	
	/* Copy into SrcVals structure */
	for(vtx = 0; vtx < IcoSurf->nvertices; vtx++)
	  MRIFseq_vox(SrcVals,vtx,0,0,nthsubj)=MRIFseq_vox(tmpmri2,vtx,0,0,0);
	
	MRIfree(&tmpmri);
	MRISfree(&SurfReg);
	
	printf("\n\n");
	fflush(stdout);
      }
    }
    else{
      printf("INFO: synthesizing randn, seed = %d\n",synthseed);fflush(stdout);
      srand48(synthseed);
      MRIrandn(IcoSurf->nvertices, 1, 1,nsubjects, 0.0, 1, SrcVals);
      if(nsmooth > 0) MRISsmoothMRI(IcoSurf, SrcVals, nsmooth, SrcVals);
    }

    if(yid != NULL)
      if(MRIwriteAnyFormat(SrcVals,yid,yfmt,-1,NULL)) exit(1);

    /* Future: run permutation loop here to get null dist of t */
    /* Must be done in one shot with estimation */

    printf("INFO: computing beta \n");fflush(stdout);
    beta = fMRImatrixMultiply(SrcVals, Q, NULL);
    if(betaid != NULL) 
      if(MRIwriteAnyFormat(beta,betaid,betafmt,-1,NULL)) exit(1);

    printf("INFO: computing eres \n");fflush(stdout);
    eres = fMRImatrixMultiply(SrcVals, R, NULL);
    if(eresid != NULL)
      if(MRIwriteAnyFormat(eres,eresid,eresfmt,-1,NULL)) exit(1);

    if(yhatid != NULL){
      printf("INFO: computing yhat \n");fflush(stdout);
      yhat = fMRImatrixMultiply(SrcVals, T, NULL);
      if(MRIwriteAnyFormat(yhat,yhatid,yhatfmt,-1,NULL)) exit(1);
      MRIfree(&yhat);
    }
    MRIfree(&SrcVals);

    printf("INFO: computing var \n");fflush(stdout);
    eresvar = fMRIvariance(eres,DOF,0,NULL);
    if(eresvarid != NULL)
      if(MRIwriteAnyFormat(eresvar,eresvarid,eresvarfmt,0,IcoSurf)) exit(1);
    MRIfree(&eres);
  }
  else{
    /*--- Load previously computed beta and res err var ----*/
    printf("INFO: loading beta_in %s\n",beta_in_id);fflush(stdout);
    if(beta_in_fmt != NULL) beta = MRIreadType(beta_in_id,beta_in_fmtid);
    else                    beta = MRIread(beta_in_id);
    if(beta == NULL){
      printf("ERROR: loading %s\n",beta_in_id);
      exit(1);
    }
    if(beta->nframes != nregressors){
      printf("ERROR: number of components in input beta does not match\n");
      printf("       the number of regressors in the design matrix\n");
      exit(1);
    }

    printf("INFO: loading var_in %s\n",eresvar_in_id);
    if(eresvar_in_fmt != NULL) eresvar = MRIreadType(eresvar_in_id,eresvar_in_fmtid);
    else                       eresvar = MRIread(eresvar_in_id);
    if(eresvar == NULL){
      printf("ERROR: loading %s\n",eresvar_in_id);
      exit(1);
    }
    if(eresvar->nframes != 1){
      printf("ERROR: number of frames in input var does not equal 1.\n");
      exit(1);
    }
    
    if(beta->width  != eresvar->width ||
       beta->height != eresvar->height ||
       beta->depth  != eresvar->depth){
      printf("ERROR: dimension mismatch between input beta and var.\n");
      exit(1);
    }
  }


  if(tid != NULL || sigid != NULL || cesid != NULL){
    printf("INFO: computing contrast effect size \n");fflush(stdout);
    ces = fMRImatrixMultiply(beta,C,NULL);
    if(cesid != NULL){
      if(IsSurfFmt(cesfmt) && IcoSurf == NULL)
	IcoSurf = MRISloadSurfSubject(trgsubject,hemi,surfregid,SUBJECTS_DIR);
      if(MRIwriteAnyFormat(ces,cesid,cesfmt,0,IcoSurf)) exit(1);
    }
  }

  if(tid != NULL || sigid != NULL){
    printf("INFO: computing t \n"); fflush(stdout);
    t = fMRIcomputeT(ces, X, C, eresvar, NULL);
    if(tid != NULL) {
      if(IsSurfFmt(tfmt) && IcoSurf == NULL)
	IcoSurf = MRISloadSurfSubject(trgsubject,hemi,surfregid,SUBJECTS_DIR);
      if(MRIwriteAnyFormat(t,tid,tfmt,0,IcoSurf)) exit(1);
    }
  }

  if(sigid != NULL){
    printf("INFO: computing t significance \n");fflush(stdout);
    sig = fMRIsigT(t, DOF, NULL);
    MRIlog10(sig,sig,1);
    if(sigfmt != NULL){
      if(IsSurfFmt(sigfmt) && IcoSurf == NULL)
      IcoSurf = MRISloadSurfSubject(trgsubject,hemi,surfregid,SUBJECTS_DIR);
      if(MRIwriteAnyFormat(sig,sigid,sigfmt,0,IcoSurf)) exit(1);
    }
  }

  printf("done \n");
  return(0);
}
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv)
{
  extern MATRIX *X;
  int  nargc , nargsused;
  char **pargv, *option ;
  int m;
  float fvtmp[1000];
  float Xcondition;

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

    else if (!strcmp(option, "--synth")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&synthseed);
      nargsused = 1;
    }
    else if (!strcmp(option, "--icoorder")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&IcoOrder);
      nargsused = 1;
    }
    else if (!strcmp(option, "--hemi")){
      if(nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--sd")){
      if(nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      nargsused = 1;
    }
    else if( !strcmp(option, "--trgsubj") || !strcmp(option, "--ts") ){
      if(nargc < 1) argnerr(option,1);
      trgsubject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--surfmeas")){
      if(nargc < 1) argnerr(option,1);
      surfmeasure = pargv[0];
      inputfmt = "curv";
      inputfmtid = checkfmt(inputfmt);
      nargsused = 1;
    }
    else if (!strcmp(option, "--i")){
      if(nargc < 2) argnerr(option,2);
      nargsused = 0;
      while( nth_is_arg(nargc, pargv, nargsused) ){
	inputlist[nargsused] = pargv[nargsused];
	nargsused ++;
	ninputs++;
      }
      printf("INFO: found %d input files on cmdline \n",ninputs);
    }
    else if (!strcmp(option, "--ifmt")){
      if(nargc < 1) argnerr(option,1);
      inputfmt = pargv[0];
      inputfmtid = checkfmt(inputfmt);
      nargsused = 1;
    }
    else if ( !strcmp(option, "--design") ){
      if(nargc < 1) argnerr(option,1);
      desmtxfname = pargv[0];
      nargsused = 1;
      ReadDesignMatrix(desmtxfname);
      Xcondition = MatrixNSConditionNumber(X);
      printf("INFO: Design Matrix Condition Number is %g\n",Xcondition);
      if(xmatfile != NULL && debug) {
	printf("INFO: Writing mat file to %s\n",xmatfile);
	if(MatlabWrite(X,xmatfile,"X")){
	  printf("ERROR: Writing mat file to %s\n",xmatfile);
	  exit(1);
	}
      }
      if(Xcondition > 100000){
	printf("ERROR: Design matrix is badly conditioned, check for linear\n"
	       "dependency  between columns (ie, two or more columns \n"
	       "that add up to another column).\n\n");
	exit(1);
      }
    }
    else if ( !strcmp(option, "--xmat") ){
      if(nargc < 1) argnerr(option,1);
      xmatfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--nsmooth")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nsmooth);
      nargsused = 1;
    }
    else if (!strcmp(option, "--frame")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
    }
    else if (!strcmp(option, "--contrast")){
      if(nargc < 1) argnerr(option,1);
      conmtxfname = pargv[0];
      C = ReadAsciiMatrix(conmtxfname);
      if(C==NULL) exit(1);
      nargsused = 1;
    }
    else if (!strcmp(option, "--gcv")){
      if(nargc < 1) argnerr(option,1);
      nargsused = 0;
      while( nth_is_arg(nargc, pargv, nargsused) ){
	sscanf(pargv[nargsused],"%f",&fvtmp[nargsused]); 
	nargsused ++;
      }
      printf("INFO: found %d elements in group contrast vector\n",
	     nargsused);
      C = MatrixAlloc(1,nargsused,MATRIX_REAL);
      for(m=0; m < nargsused; m++) C->rptr[1][m+1] = fvtmp[m];
      MatrixPrint(stdout,C);
    }
    else if (!strcmp(option, "--beta")){
      if(nargc < 1) argnerr(option,1);
      betaid = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	betafmt = pargv[1]; nargsused ++;
	betafmtid = checkfmt(betafmt);
      }
      else betafmtid = getfmtid(betaid);
      if(IsSurfFmt(betafmt)){
	printf("ERROR: cannot use curv or paint as output format for beta\n");
	exit(1);
      }
    }
    else if (!strcmp(option, "--beta_in")){
      if(nargc < 1) argnerr(option,1);
      beta_in_id = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	beta_in_fmt = pargv[1]; nargsused ++;
	beta_in_fmtid = checkfmt(beta_in_fmt);
      }
    }
    else if (!strcmp(option, "--ces")){
      if(nargc < 1) argnerr(option,1);
      cesid = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	cesfmt = pargv[1]; nargsused ++;
	cesfmtid = checkfmt(cesfmt);
      }
      else cesfmtid = getfmtid(cesid);
    }
    else if (!strcmp(option, "--eres")){
      if(nargc < 1) argnerr(option,1);
      eresid = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	eresfmt = pargv[1]; nargsused ++;
	eresfmtid = checkfmt(eresfmt);
      }
      else eresfmtid = getfmtid(eresid);
    }
    else if (!strcmp(option, "--y")){
      if(nargc < 1) argnerr(option,1);
      yid = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	yfmt = pargv[1]; nargsused ++;
	yfmtid = checkfmt(yfmt);
      }
      else yfmtid = getfmtid(yid);
    }
    else if (!strcmp(option, "--yhat")){
      if(nargc < 1) argnerr(option,1);
      yhatid = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	yhatfmt = pargv[1]; nargsused ++;
	yhatfmtid = checkfmt(yhatfmt);
      }
      else yhatfmtid = getfmtid(yhatid);
    }
    else if (!strcmp(option, "--var")){
      if(nargc < 1) argnerr(option,1);
      eresvarid = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	eresvarfmt = pargv[1]; nargsused ++;
	eresvarfmtid = checkfmt(eresvarfmt);
      }
      else eresvarfmtid = getfmtid(eresvarid);
      if(stringmatch(eresvarfmt,"curv")){
	printf("ERROR: cannot use curv as output format for var\n");
	exit(1);
      }
    }
    else if (!strcmp(option, "--var_in")){
      if(nargc < 1) argnerr(option,1);
      eresvar_in_id = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	eresvar_in_fmt = pargv[1]; nargsused ++;
	eresvar_in_fmtid = checkfmt(eresvar_in_fmt);
      }
      else eresvar_in_fmtid = getfmtid(eresvar_in_id);
    }
    else if (!strcmp(option, "--t")){
      if(nargc < 1) argnerr(option,1);
      tid = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	tfmt = pargv[1]; nargsused ++;
	tfmtid = checkfmt(tfmt);
      }
      else tfmtid = getfmtid(tid);
      if(stringmatch(tfmt,"curv")){
	printf("ERROR: cannot use curv as output format for t\n");
	exit(1);
      }
    }
    else if (!strcmp(option, "--sigt")){
      if(nargc < 1) argnerr(option,1);
      sigid = pargv[0];
      nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	sigfmt = pargv[1]; nargsused ++;
	sigfmtid = checkfmt(sigfmt);
      }
      else sigfmtid = getfmtid(sigid);
      if(stringmatch(sigfmt,"curv")){
	printf("ERROR: cannot use curv as output format for sigt\n");
	exit(1);
      }
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
  printf("USAGE: mris_glm (formerly mri_surfglm)\n") ;
  printf("\n");
  printf("Raw Data Input Options\n");
  printf("   --design     fname : name of design matrix (ascii)\n");
  printf("   --surfmeas name  : input file or name of surface measure\n");
  printf("   --i input1 input2 ...> : input file list\n");
  printf("   --frame      M     : use 0-based Mth frame (default is 0)\n");
  printf("   --hemi       hemi  : hemisphere (lh or rh) \n");
  printf("   --trgsubj    subject : target subject \n");
  printf("   --icoorder   order : order of icosahedral tesselation (default 7)\n");
  printf("   --nsmooth    N     : number of smoothing iterations\n");
  printf("\n");
  printf("Processed Data Input Options\n");
  printf("   --beta_in    name <fmt> : parameter estimates from previous \n");
  printf("   --var_in     name <fmt> : reserr var from previous \n");
  printf("\n");
  printf("Estimation Output Options\n");
  printf("   --y       name <fmt> : input data after resampling and smoothing\n");
  printf("   --beta    name <fmt> : parameter estimates \n");
  printf("   --var     name <fmt> : variance of residual error \n");
  printf("   --yhat    name <fmt> : signal estimate\n");
  printf("   --eres    name <fmt> : residual error \n");
  printf("   --xmat    matfile    : save design matrix in matlab4 format\n");
  printf("\n");
  printf("Group Contrast Options\n");
  printf("   --contrast   fname : file containing group contrast matrix (ascii)\n");
  printf("   --gcv c1 ... cN : group contrast vector specified on cmdline\n");
  printf("   --ces     name <fmt> : contrast effect size  \n");
  printf("   --t       name <fmt> : t-ratio of contrast \n");
  printf("   --sigt    name <fmt> : signficance of t-ratio (ie, t-Test) \n");
  printf("\n");
  printf("   --synth seed : substitute white gaussian noise for data\n");
  printf("   --sd    subjectsdir : default is env SUBJECTS_DIR\n");
  printf("\n");
  printf("   --version : print version and exit\n");
  printf("   --help : be VERY careful with this, you MIGHT learn something.\n");
  printf("\n");
  printf("%s\n",vcid);
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(

"
SUMMARY

NOTE: this program replaces mri_surfglm (this is just a name change).

This program performs inter-subject/group averaging and inference on
the surface by fitting a GLM model at each vertex. The model consists
of subject parameters (eg, age, gender, etc). The model is the same
across all vertices, though the fit may be (will be) different. The
user must supply a matrix that represents the GLM. While estimation
and inference can be performed in a single call to mris_glm, the tasks
can also be separated, which can be much more convenient. Inferences
are not corrected for multiple comparisons.

MATHEMATICAL BACKGROUND

The forward model is given by:

    y = XB + n

where X is the Ns-by-Nb design matrix, y is the Ns-by-Nv raw data set,
B is the Nb-by-Nv regression parameters, and n is noise. Ns is the
number of subjects, Nb is the number of regressors, and Nv is the 
number of vertices. y will have been preprocessed in possibly two
ways: (1) it will be sampled on the surface of the target subject, 
and (2) it may be spatially smoothed (prior to resampling) (see
--nsmooth). 

During the estimation stage, the forward model is inverted to
solve for B:

    B = inv(X'*X)*X'y

This is performed at each vertex, the result of which can be saved
with the --beta flag.

The signal estimate (which can be saved with --yhat) is computed as 

    yhat = B*X

The residual error (which can be saved with --eres) is computed as 

    eres = y - yhat

The noise variance estimate (computed for each vertex) is computed
as the sum of the squares of the residual error divided by the DOF.
The DOF equals the number of rows of X minus the number of columns.
The noise variance can be saved with --var.

A contrast vector C has as many elements as columns of X. The 
contrast effect size (--ces) is then computed as:

   G = C*B

The t-ratio (--t) for the contrast is then given by:

   t = G/sqrt(var * C*inv(X'X)*C')

The signifiance of the t-ratio (based on a double-sided t-test
and uncorrected for multiple comparisons) can be saved with 
the --sigt flag.


COMMAND-LINE ARGUMENTS

--design fname

File name for design matrix. The design matrix must be in an ASCII
file. The first column must be the id assigned to the subject 
during FreeSurfer reconstruction. The following columns are the
the design matrix. It is possible for the design matrix to be
ill-conditioned. This means that two columns are identical or that
one column is equal to a weighted sum of any of the other columns.
In this case, the matrix cannot be inverted and so the analysis
must stop. The matrix is judged to be ill-conditioned if its 
condition number is greater than 100000. It is possible to test
the condition of a matrix without having to include all the 
command-line options needed for a full analysis. Just run:
    mris_glm --design fname
It will print out an INFO line with the condition number.

--surfmeas name 

This is one of the two ways the user can specify the raw data (ie, the
input to the estimation stage). This method requires that the data
file for each subject reside in the FreeSurfer anatomical directory
under the surf subdirectory. This frees the user from having to list
all the inputs on the command-line. The name can be one of two things,
depending upon the format designation (--ifmt). If the format is curv
(or unspecified), then mris_glm will construct the name of the
input file as hemi.name (where hemi is either lh or rh as specified 
by --hemi). If the format is anything else, then it looks for a file
called name in the surf subdirectory. Only specify a raw data input
when performing estimation.

--i input1 input2 ...

This is second method that the user can specify the raw data (ie, the
input to the estimation stage). This method allows the user to specify
all of the input files on the command line. There must be as many
files listed as there are rows in the design matrix. The format (same
for all inputs) is as designated with the --ifmt flag. If the format
is unspecified, mris_glm will attempt to determine the format.
curv format cannot be used with explicit inputs. Only specify a raw 
data input when performing estimation.

--ifmt format

This flag is used to specify the format of the input raw data files.
When a surface measure is specified, the default format becomes curv,
otherwise the mris_glm will attempt to determine the format from
the file name if a format is not explicitly give. It is not possible
to use curv format with explicit inputs (ie, --i). Valid formats are:
bshort, bfloat, COR, analyze, spm, paint (or w or wfile), and curv.

--frame M

This allows the user to specify that the Mth frame of the raw data (if 
the input format supports multiple frames) should be used as input. 
M is zero-based. Default is 0.

--hemi hemisphere

Specify that the input data either the left (lh) or right (rh) hemisphere.

--trgsubject subject

Resample each input data set from the individual's surface to that of the
target subject. 

--icoorder order

When the target subject is ico, this specifies the order of the 
icosahedron. Default is 7 (163842 vertices). 

--nsmooth N

Perform N iterations of nearest-neighbor spatial smoothing. Note: 
smoothing is  performed on the surface of each source subject before 
resampling to the target subject.

--beta_in betaname <fmt>

This flag (with --var_in) allows the user to use a previous estimation 
as input to the contrast/inference computation. This arguments should 
be identical to that specified with --beta when doing the estimation.
Note: the user must also specify a pointer to the residual error
variance using --var_in.

--var_in betaname <fmt>

This flag (with --beta_in) allows the user to use a previous estimation 
as input to the contrast/inference computation. This arguments should be
identical to that specified with --var when doing the estimation.

--y name <fmt>

Save the raw data (after resampling and smoothing) into a single volume.
This is mainly good for debugging. fmt is the format (see OUTPUT FORMATS).

--beta name <fmt>

Save the result of estimation (ie, the map regression coefficients). 
This can be used in subsequent calls to mris_glm to perform.
inference. fmt is the format (see OUTPUT FORMATS).

--var name <fmt>

Save the estimate of the noise variances (estimated as the variance
of the residual error). This can be used in subsequent calls to 
mris_glm to perform inference.

--yhat name <fmt>

Save the estimate of the signal. This is only good for debugging.
fmt is the format (see OUTPUT FORMATS).

--eres name <fmt>

Save the residual error. This is only good for debugging.
fmt is the format (see OUTPUT FORMATS).

--xmat name 

Save the design matrix in matlab4 format. This is only good for 
debugging. fmt is the format (see OUTPUT FORMATS).

--contrast fname

Load the group contrast vector from the ascii file fname. The 
contrast vector should have as many entries as there are columns
in the design matrix. The contrast vector can also be specified 
on the command-line with --gcv.

--gcv c1 ... cN

Specify the group contrast vector (gcv) on the command-line as the
values c1 ... cN. The contrast vector should have as many entries as 
there are columns in the design matrix. The contrast vector can also 
be specified in a file with --contrast.

--ces name <fmt>

Save the contrast effect size from the contrast vector specified
with --contrast or --gcv. fmt is the format (see OUTPUT FORMATS).

--t name <fmt>

Save the t-ratio of the contrast.

--sigt name <fmt>

Save the signficance (p-value) of the t-ratio of the contrast. The
value is actually the -log10 of the significance with the same
sign as the t-ratio from which it was computed. The significance
is computed from a double-sided t-test and is NOT corrected for
multiple comparisons across space. fmt is the format (see OUTPUT 
FORMATS).

--synth seed

Substitute white gaussian noise for the data. Good for debugging.

--sd subjectsdir

Look for FreeSurfer reconstructions in subjectsdir. If unspecified,
the SUBJECTS_DIR envionment variable will be used.


OUTPUT FORMATS:

Output formats can be designated by specifying a format string
following the the output name. Valid strings include bfloat,
bshort, spm, analyze, analyze4d, COR, paint, w, and wfile. Paint,
w, and wfile cannot be used for beta, y, yhat, or eres.


EXAMPLES:

1. Analyze thickness maps based on gender and age for 5 hypothetical
subjects: subj1 (m, 22) , subj2 (m, 57), subj3 (f, 33), subj4 (f, 65), 
subj5 (m, 27). The design matrix would look something like:

   subj1  1  0  22
   subj2  1  0  57
   subj3  0  1  33
   subj4  0  1  65
   subj5  1  0  27

The first column is the name of the subject as it appears in the
FreeSurfer SubjectsDir. The second and third columns categorically
code gender (first column male, second column female). The last
column codes age. Assume this matrix is stored in a file called
genage.mtx

  mris_glm --design genage.mtx --hemi lh --surfmeas thickness 
    --trgsubj average7 --nsmooth 50 --beta beta bfloat 
    --var var bfloat 

This will read the thickness maps for each of the subjects, smooth
it with 50 iterations of nearest-neighbor smoothing, resample to
the average7 subject and save the regression coeffients and
noise variance, both in bfloat format No inference was performed.

2. Test the data in Example 1 for an effect of age:

  mris_glm --design genage.mtx --hemi lh --trgsubj average7 
      --beta_in beta bfloat --var_in var bfloat 
      --gcv 0 0 1 --sigt ./age-sigt-lh.w paint

3. Test the data in Example 1 for a difference between males 
and females with age regressed out:

  mris_glm --design genage.mtx --hemi lh --trgsubj average7 
      --beta_in beta bfloat --var_in var bfloat 
      --gcv 1 -1 0 --sigt ./gender-sigt-lh.w paint

4. Perform the same analysis as done in Example 1, but use 
values that have been painted onto the surface of each subject
(this could have come from an fMRI analysis):

  mris_glm --design genage.mtx --hemi lh
    --ifmt paint --i ./subj1-data-lh.w ./subj2-data-lh.w 
    ./subj3-data-lh.w ./subj4-data-lh.w ./subj5-data-lh.w
    --trgsubj average7 --nsmooth 50 --beta beta bfloat 
    --var var bfloat 

BUGS

No correction for multiple comparisons.

BUG REPORTING

If you want your bug report or question to have a prayer of being
answered, make sure to include the following information (send to
freesurfer@surfer.nmr.mgh.harvard.edu): (1) version of mris_glm
(run with --version), (2) full command-line used, (3) terminal
output, (4) description of the problem. 



"


);


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
  if(desmtxfname == NULL){
    printf("ERROR: must specify a design matrix\n");
    exit(1);
  }

  if( (beta_in_id != NULL && eresvar_in_id == NULL) ||
      (beta_in_id == NULL && eresvar_in_id != NULL) ){
    printf("ERROR: if using a previous estimate, must specify both\n");
    printf("       beta and var inputs\n");
    exit(1);
  }

  if(ninputs == 0 && surfmeasure == NULL && beta_in_id == NULL){
    printf("ERROR: must specify an input, either surface measure or \n");
    printf("       a list of files or previous estimate (beta_in, var_in)\n");
    exit(1);
  }
  if(surfmeasure != NULL && beta_in_id != NULL){
    printf("ERROR: cannot specify both surface measure and beta\n");
    exit(1);
  }
  if(ninputs != 0 && beta_in_id != NULL){
    printf("ERROR: cannot specify both --i input and beta\n");
    exit(1);
  }
  if(surfmeasure != NULL && ninputs != 0){
    printf("ERROR: cannot specify both surface measure and --i input\n");
    exit(1);
  }

  if(eresid != NULL && beta_in_id != NULL){
    printf("ERROR: cannot specify an eres output with a previous estimate\n");
    printf("       as input.\n");
    exit(1);
  }

  if(yid != NULL && beta_in_id != NULL){
    printf("ERROR: cannot specify a y output with a previous estimate.\n");
    printf("       as input.\n");
    exit(1);
  }

  if(ces != NULL || tid != NULL || sigid != NULL){
    if(C == NULL){
      printf("ERROR: must specify a contrast matrix with ces, t, or sig output\n");
      exit(1);
    }
    if(C != NULL){
      if(C->cols != nregressors){
	printf("ERROR: dimension mismatch: ncols in contrast matrix\n");
	printf("       is %d but the number of regressors is %d\n",
	       C->cols,nregressors);
	exit(1);
      }      
      if(C->rows != 1){
	printf("ERROR: the contrast matrix can only have one row.\n");
	exit(1);
      }
    }
  }

  if(beta_in_id != NULL && cesid == NULL && tid == NULL && sigid == NULL){
    printf("ERROR: nothing to do. You have specified a previous\n");
    printf("       estimate as input but have not specified an output.\n");
    printf("       The output should be ces, t, or sigt\n");
    exit(1);
  }

  if(trgsubject == NULL){
    printf("ERROR: no target subject specified.\n");
    exit(1);
  }
  if(hemi == NULL){
    printf("ERROR: no hemisphere specified.\n");
    exit(1);
  }


  if(SUBJECTS_DIR == NULL){
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if(SUBJECTS_DIR==NULL){
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  if(beta_in_id == NULL){
    if(hemi == NULL){
      printf("ERROR: must specify a hemisphere.\n");
      exit(1);
    }
    if(strcmp(hemi,"lh") && strcmp(hemi,"rh")){
      printf("ERROR: hemi = %s, must be either lh or rh\n",hemi);
      exit(1);
    }
  }

  if(ninputs != 0 && ninputs != X->rows){
    printf("ERROR: number of inputs does not equal the number of rows.\n");
    printf("       in the design matrix\n");
    exit(1);
  }



  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
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
/*---------------------------------------------------------------*/
static int isflag(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth)
{
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if(nargc <= nth) return(0); 

  /* check whether the nth arg is a flag */
  if(isflag(argv[nth])) return(0);

  return(1);
}
/*------------------------------------------------------------*/
static int stringmatch(char *str1, char *str2)
{
  if(! strcmp(str1,str2)) return(1);
  return(0);
}
/*------------------------------------------------------------*/
static int getfmtid(char *fname)
{
  int fmtid;
  fmtid = mri_identify(fname);
  if(fmtid == MRI_VOLUME_TYPE_UNKNOWN){
    printf("ERROR: cannot determine format of %s\n",fname);
    exit(1);
  }
  return(fmtid);
}

/*------------------------------------------------------------*/
static int checkfmt(char *fmt)
{
  int fmtid;

  if(fmt == NULL) return(MRI_VOLUME_TYPE_UNKNOWN);
  if(stringmatch(fmt,"curv") ||
     stringmatch(fmt,"paint") ||
     stringmatch(fmt,"wfile") ||
     stringmatch(fmt,"w")) return(MRI_VOLUME_TYPE_UNKNOWN);

  fmtid = string_to_type(fmt);
  if(fmtid == MRI_VOLUME_TYPE_UNKNOWN){
    printf("ERROR: format string %s unrecognized\n",fmt);
    exit(1);
  }
  return(fmtid);
}

/*------------------------------------------------------------*/
static int IsSurfFmt(char *fmt)
{
  if(fmt == NULL) return(0);
  if(stringmatch(fmt,"curv") ||
     stringmatch(fmt,"paint") ||
     stringmatch(fmt,"wfile") ||
     stringmatch(fmt,"w")) return(1);
  return(0);
}

/*------------------------------------------------------------*/
int ReadDesignMatrix(char *desmtxfname)
{
  extern MATRIX *X;
  extern char *subjectlist[1000];
  extern int nsubjects, nregressors;
  int nrows,ncols, r,c;
  char tmpstring[1001];
  FILE *fp;

  ReadAsciiMatrixSize(desmtxfname, &nrows, &ncols);
  nsubjects = nrows;
  nregressors = ncols-1;
  X = MatrixAlloc(nrows,ncols-1,MATRIX_REAL);

  fp = fopen(desmtxfname,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",desmtxfname);
    exit(1);
  }

  for(r=0; r < nrows; r++){

    fscanf(fp,"%s",tmpstring);
    subjectlist[r] = (char *)calloc(strlen(tmpstring)+1,sizeof(char));
    memcpy(subjectlist[r],tmpstring,strlen(tmpstring)+1);
    //printf("%2d %s\n",r+1,subjectlist[r]);

    for(c=0; c < ncols-1; c++)
      fscanf(fp,"%f",&(X->rptr[r+1][c+1]));

  }
  //MatrixPrint(stdout,X);

  fclose(fp);
  return(0);
}
/*------------------------------------------------------------*/
int ReadAsciiMatrixNRows(char *desmtxfname)
{
  FILE *fp;
  int nrows;
  char tmpstring[2001];

  fp = fopen(desmtxfname,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",desmtxfname);
    return(-1);
  }

  nrows = 0;
  while(fgets(tmpstring,2000,fp) != NULL)  nrows ++;
  fclose(fp);

  if(nrows == 0){
    printf("ERROR: no data found in %s\n",desmtxfname);
    return(-1);
  }

  if(nrows < 0){
    printf("ERROR: reading number of rows from %s\n",desmtxfname);
    return(-1);
  }

  return(nrows);
}
/*-----------------------------------------------------------------*/
int ReadAsciiMatrixSize(char *asciimtxfname, int *pnrows, int *pncols)
{
  FILE *fp;
  int nrows, nitems, ncols;
  char tmpstring[2001];

  nrows = ReadAsciiMatrixNRows(asciimtxfname);
  if(nrows < 0) return(1);

  fp = fopen(asciimtxfname,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",asciimtxfname);
    return(1);
  }

  nitems = 0;
  while(fscanf(fp,"%s",tmpstring) != EOF)  nitems ++;
  fclose(fp);

  if(nitems == 0){
    printf("ERROR: no items found in %s\n",asciimtxfname);
    return(1);
  }

  ncols = nitems/nrows;
  if(ncols*nrows != nitems){
    printf("ERROR: number of items (%d) not divisible by nrows (%d)\n",
	   nitems,nrows);
    return(1);
  }

  if(ncols < -1){
    printf("ERROR: reading number of cols from %s\n",desmtxfname);
    return(-1);
  }

  *pnrows = nrows;
  *pncols = ncols;

  return(0);
}
/*------------------------------------------------------------*/
MATRIX *ReadAsciiMatrix(char *asciimtxfname)
{
  int err, nrows,ncols, r,c, nread;
  FILE *fp;

  err = ReadAsciiMatrixSize(asciimtxfname, &nrows, &ncols);
  if(err) return(NULL);

  C = MatrixAlloc(nrows,ncols,MATRIX_REAL);

  fp = fopen(asciimtxfname,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",asciimtxfname);
    exit(1);
  }

  for(r=0; r < nrows; r++){
    for(c=0; c < ncols; c++){
      nread = fscanf(fp,"%f",&(C->rptr[r+1][c+1]));
      if(nread != 1){
	printf("ERROR: ReadAsciiMatrix: could not read item %d,%d\n",r,c);
	MatrixFree(&C);
	return(NULL);
      }
    }
  }

  fclose(fp);
  return(C);
}
/*--------------------------------------------------------------------*/
MRIS *MRISloadSurfSubject(char *subj, char *hemi, char *surfid, 
			  char *SUBJECTS_DIR)
{
  MRIS *Surf;
  char fname[2000];

  if(SUBJECTS_DIR == NULL){
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if(SUBJECTS_DIR==NULL){
      printf("ERROR: SUBJECTS_DIR not defined in environment\n");
      return(NULL);
    }
  }

  sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subj,hemi,surfid);
  printf("  INFO: loading surface  %s\n",fname);
  fflush(stdout);
  Surf = MRISread(fname) ;
  if(Surf == NULL){
    printf("ERROR: could not load registration surface\n");
    exit(1);
  }
  printf("nvertices = %d\n",Surf->nvertices);fflush(stdout);

  return(Surf);
}
