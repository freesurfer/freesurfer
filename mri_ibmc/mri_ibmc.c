/**
 * @file  mri_ibmc.c
 * @brief Intersection-based motion correction
 *
 * Intersection-based motion correction based on Kim, et al, TMI,
 * 2010. Registers three volumes.
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2011/07/05 23:00:37 $
 *    $Revision: 1.9 $
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
#include <float.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "fsenv.h"
#include "registerio.h"
#include "matrix.h"
#include "mri2.h"
#include "transform.h"
#include "timer.h"
#include "numerics.h"
#include "resample.h"


double round(double);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_ibmc.c,v 1.9 2011/07/05 23:00:37 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;
MATRIX *MRIangles2RotMatB(double *angles, MATRIX *R);

#define IBMC_A2B_EQUAL 0
#define IBMC_A2B_RIGID 1

/*---------------------------------------------------------------*/
typedef struct 
{
  MRI *mri;       // pixel data for slice
  MATRIX *Rbeta;  // tkreg matrix for slice beyond Rv
  MATRIX *Rs;     // tkreg matrix for slice 
  double beta[6]; // motion params
  int Vno, Sno;   // volume stack and slice numbers
  MATRIX *D, *iD; // index stack-to-slice (and inv)
  MATRIX *Ts;     // tkr vox2ras
  MATRIX *iRs;    // inv(Rs)
  MATRIX *TsDiTv; // cache
  MATRIX *Rs0;    // tkreg matrix for slice if Rbeta=I
} IBMC_SLICE;
typedef struct
{
  int a2bmethod;    // method to convert alpha to beta
  int nalpha;       // number of parameters for this stack
  double *alpha;    // alpha parameters for this stack
} IBMC_PARAMS;
/*---------------------------------------------------------------*/
typedef struct 
{
  MRI *mri;
  MATRIX *Rv,*iRv;  // Stack tkreg (for init)
  MATRIX *Tv,*iTv;  // Stack tkr vox2ras
  IBMC_SLICE **slice; // slices in the stack
  int Vno;            // volume stack number 
  IBMC_PARAMS *params;
} IBMC_STACK;
/*---------------------------------------------------------------*/
typedef struct 
{
  int nstackss; // Number of input volume stacks
  IBMC_STACK *stacks; // Input stacks
} IBMC;
/*---------------------------------------------------------------*/
IBMC_SLICE *IBMCinitSlice(IBMC_STACK *stack, int Sno)
{
  IBMC_SLICE *slice;

  slice = (IBMC_SLICE *) calloc(1,sizeof(IBMC_SLICE));
  slice->Sno = Sno;
  slice->mri = MRIcrop(stack->mri, 0,0,Sno, stack->mri->width-1,stack->mri->height-1,Sno);
  slice->Ts  = MRIxfmCRS2XYZtkreg(slice->mri);
  slice->D = MatrixIdentity(4,NULL);
  slice->D->rptr[3][4] = -Sno;
  slice->iD   = MatrixInverse(slice->D,NULL);
  // Cache TsDiTv = Ts*D*inv(Tv)
  slice->TsDiTv  = MatrixMultiply(slice->Ts,slice->D,NULL);
  slice->TsDiTv  = MatrixMultiply(slice->TsDiTv,stack->iTv,slice->TsDiTv);
  // Rs0 = Ts*D*inv(Tv)*Rv = Rs when Rbeta=I
  slice->Rs0  = MatrixMultiply(slice->TsDiTv,stack->Rv,NULL);
  slice->Rbeta = MatrixAlloc(4,4,MATRIX_REAL);
  slice->Rs    = MatrixAlloc(4,4,MATRIX_REAL);
  slice->iRs   = MatrixAlloc(4,4,MATRIX_REAL);
  return(slice);
}
/*---------------------------------------------------------------*/
int IBMCnAlpha(IBMC_STACK *stack)
{
  stack->params->nalpha = -1;
  switch(stack->params->a2bmethod){
  case IBMC_A2B_EQUAL:
    stack->params->nalpha = 6*stack->mri->depth;
    break;
  case IBMC_A2B_RIGID:
    stack->params->nalpha = 6;
    break;
  default:
    printf("ERROR: A2B method %d not regocnized\n",stack->params->a2bmethod);
  }
  return(stack->params->nalpha);
}

/*---------------------------------------------------------------*/
IBMC_PARAMS *IBMCallocParams(int nalpha)
{
  IBMC_PARAMS *p;
  p = (IBMC_PARAMS *) calloc(sizeof(IBMC_PARAMS),1);
  p->nalpha = nalpha;
  p->alpha = (double *) calloc(sizeof(double),p->nalpha);
  return(p);
}

/*---------------------------------------------------------------*/
int IBMCrandParams(IBMC_PARAMS *p)
{
  int k;
  for(k=0; k < p->nalpha; k++) p->alpha[k] = 2*(drand48()-0.5);
  return(0);
}

/*---------------------------------------------------------------*/
IBMC_PARAMS *IBMCreadParams(char *fname)
{
  IBMC_PARAMS *p=NULL;
  FILE *fp;
  int k,magic, a2bmethod, nalpha, IsASCII, version;
  size_t nitems;
  char tmpstr[101];

  fp = fopen(fname,"r");
  if(fp == NULL){
    printf("ERROR: could not open %s for reading\n",fname);
    return(NULL);
  }

  fscanf(fp,"%*s %*s %d",&version);
  fscanf(fp,"%*s %d",&a2bmethod);
  fscanf(fp,"%*s %d",&nalpha);
  fscanf(fp,"%*s %d",&IsASCII);
  fprintf(stderr,"IBMCreadParams: %s %d %d %d\n",fname,a2bmethod,nalpha,IsASCII);
  p = IBMCallocParams(nalpha);
  p->a2bmethod = a2bmethod;

  if(IsASCII){
    // As written by IBMCprintParams()
    for(k=0; k < nalpha; k++) fscanf(fp,"%*d %lf",&(p->alpha[k]));
  }
  else{
    // As written by IBMCwriteParams()
    fgets(tmpstr,100,fp); // swallow the new line
    fread(&magic,sizeof(int),1,fp);
    nitems = fread(p->alpha,sizeof(double),nalpha,fp);
    if(nitems != nalpha){
      printf("ERROR: %s tried to read %d, actually read %ld\n",fname,nalpha,nitems);
      return(NULL);
    }
  }
  fclose(fp);
  
  return(p);
}
/*---------------------------------------------------------------*/
int IBMCwriteParams(char *fname, IBMC_PARAMS *p)
{
  FILE *fp;
  int magic=1;

  fp = fopen(fname,"w");
  if(fp == NULL){
    printf("ERROR: could not open %s\n",fname);
    return(1);
  }
  fprintf(fp,"IBMC_PARAMS Version 1\n");
  fprintf(fp,"a2bmethod %d\n",p->a2bmethod);
  fprintf(fp,"nalpha    %d\n",p->nalpha);
  fprintf(fp,"IsASCII   %d\n",0);
  fwrite(&magic,sizeof(int),1,fp);
  fwrite(p->alpha,sizeof(double),p->nalpha,fp);

  fclose(fp);
  return(0);
}
/*---------------------------------------------------------------*/
int IBMCprintParams(FILE *fp, IBMC_PARAMS *p)
{
  int k;
  fprintf(fp,"IBMC_PARAMS Version 1\n");
  fprintf(fp,"a2bmethod %d\n",p->a2bmethod);
  fprintf(fp,"nalpha    %d\n",p->nalpha);
  fprintf(fp,"IsASCII   %d\n",1);
  for(k=0; k < p->nalpha; k++)
    fprintf(fp,"%3d %10.5lf\n",k,p->alpha[k]);
  return(0);
}


/*---------------------------------------------------------------*/
IBMC_STACK *IBMCinitStack(MRI *vol, MATRIX *Rv)
{
  IBMC_STACK *stack;
  int nthslice;

  stack = (IBMC_STACK *)calloc(1,sizeof(IBMC_STACK));
  stack->mri = vol;
  stack->Rv = Rv;
  stack->iRv = MatrixInverse(stack->Rv,NULL);
  stack->slice = (IBMC_SLICE **)calloc(vol->depth,sizeof(IBMC_SLICE *));
  stack->Tv  = MRIxfmCRS2XYZtkreg(stack->mri);
  stack->iTv = MatrixInverse(stack->Tv,NULL);

  for(nthslice = 0; nthslice < stack->mri->depth; nthslice++){
    stack->slice[nthslice] = IBMCinitSlice(stack, nthslice);
    //printf("%2d %3d\n",nthslice,stack->slice[0]->vol->depth);
  }
  return(stack);
}

/*---------------------------------------------------------------*/
int IBMCwriteSlices(IBMC_STACK *stack, char *base)
{
  char tmpstr[2000];
  int nthslice;

  for(nthslice = 0; nthslice < stack->mri->depth; nthslice++){
    sprintf(tmpstr,"%s.%03d.nii",base,nthslice);
    //printf("%2d %3d %s\n",nthslice,stack->slice[nthslice]->mri->depth,tmpstr);
    MRIwrite(stack->slice[nthslice]->mri,tmpstr);
    sprintf(tmpstr,"%s.%03d.Rs0.dat",base,nthslice);
    regio_write_register(tmpstr, "fsf01anat", stack->slice[nthslice]->mri->xsize,
			 stack->slice[nthslice]->mri->zsize, .15, 
			 stack->slice[nthslice]->Rs0,FLT2INT_ROUND);
    sprintf(tmpstr,"%s.%03d.Rs.dat",base,nthslice);
    regio_write_register(tmpstr, "fsf01anat", stack->slice[nthslice]->mri->xsize,
			 stack->slice[nthslice]->mri->zsize, .15, 
			 stack->slice[nthslice]->Rs,FLT2INT_ROUND);
  }
  return(0);
}
/*---------------------------------------------------------------*/
// Copies stack pixels into an mri volume
MRI *IBMCcopyStack2MRI(IBMC_STACK *stack, MRI *mri)
{
  int c, r, s, f;
  double v;

  if(mri == NULL){
    mri = MRIallocSequence(stack->mri->width,stack->mri->height,
			   stack->mri->depth,MRI_FLOAT,stack->mri->nframes);
    MRIcopyHeader(stack->mri,mri);
    MRIcopyPulseParameters(stack->mri,mri);
  }

  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
	for(f=0; f < mri->nframes; f++){
	  v = MRIgetVoxVal(stack->slice[s]->mri,c,r,0,f);
	  MRIsetVoxVal(mri,c,r,s,f,nint(v));
	}
      }
    }
  }

  return(mri);
}
/*---------------------------------------------------------------*/
int IBMCalpha2Beta(IBMC_STACK *stack)
{
  int s,k,nthalpha;
  static double angles[3];
  IBMC_SLICE *slice;

  switch(stack->params->a2bmethod){
  case IBMC_A2B_EQUAL:
    // Just copy alphas to betas
    nthalpha = 0;
    for(s=0; s < stack->mri->depth; s++){
      for(k=0; k < 6; k++){
	stack->slice[s]->beta[k] = stack->params->alpha[nthalpha];
	nthalpha ++;
      }
    }
    break;
  case IBMC_A2B_RIGID:
    // The 6 alphas become the 6 betas for all slices
    for(s=0; s < stack->mri->depth; s++){
      for(k=0; k < 6; k++)
	stack->slice[s]->beta[k] = stack->params->alpha[k];
    }
    break;
  }

  // Set the matrix transform
  for(s=0; s < stack->mri->depth; s++){
    slice = stack->slice[s];
    for(k=0; k < 3; k++) angles[k] = slice->beta[k+3]*M_PI/180;
    slice->Rbeta = MRIangles2RotMatB(angles, slice->Rbeta);
    for(k=0; k < 3; k++) slice->Rbeta->rptr[k+1][4] = slice->beta[k];
    // For consistency, it must be Rs=Ts*D*inv(Tv)*Rbeta*Rv
    slice->Rs = MatrixMultiply(slice->TsDiTv,slice->Rbeta,slice->Rs);
    slice->Rs = MatrixMultiply(slice->Rs,stack->Rv,slice->Rs);
    slice->iRs = MatrixInverse(slice->Rs,slice->iRs);
  }

  return(0);
}
/*---------------------------------------------------------------*/
// Convert a vol (eg, anat) into the stack/slice space. For synthesis
int IBMCsynthStack(MRI *src, IBMC_STACK *dst)
{
  int nthslice;
  MRI *slcmri;
  for(nthslice=0; nthslice < dst->mri->depth; nthslice++){
    slcmri = dst->slice[nthslice]->mri;
    //printf("Rs %d --------------------- \n",nthslice);
    //MatrixPrint(stdout,dst->slice[nthslice]->Rs);
    //printf("Rbeta %d --------------------- \n",nthslice);
    //MatrixPrint(stdout,dst->slice[nthslice]->Rbeta);
    MRIconst(slcmri->width,slcmri->height,slcmri->depth,slcmri->nframes,0,slcmri);
    MRIvol2VolTkRegVSM(src, slcmri, dst->slice[nthslice]->iRs, SAMPLE_TRILINEAR, 0, NULL);
  }
  return(0);
}

/*---------------------------------------------------------------*/
// probably should not use
int IBMCsetSliceBeta(IBMC_SLICE *slice, const double beta[6])
{
  int k;
  static double angles[3];

  for(k=0; k < 6; k++) slice->beta[k] = beta[k];
  for(k=0; k < 3; k++) angles[k] = beta[k+3]*M_PI/180;
  slice->Rs = MRIangles2RotMatB(angles, slice->Rs);
  for(k=0; k < 3; k++) slice->Rs->rptr[k+1][4] = beta[k];
  slice->Rs = MatrixMultiply(slice->Rs0,slice->Rs,slice->Rs);
  slice->iRs = MatrixInverse(slice->Rs,slice->iRs);
  return(0);
}

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
MRI *targ=NULL, *mov=NULL;
MATRIX *Rtarg=NULL,*Rmov=NULL;
char *AlphaFile=NULL;
char tmpstr[2000];

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,k,nthslice;
  IBMC_STACK *stack;
  double beta[6];
  MRI *vol;

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
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);

  stack = IBMCinitStack(mov, Rmov);
  IBMCwriteSlices(stack, "base0");

  for(k=0; k<6; k++) beta[k] = 0;
  for(nthslice=0; nthslice < stack->mri->depth; nthslice++)
    IBMCsetSliceBeta(stack->slice[nthslice], beta);
  IBMCsynthStack(targ, stack);
  IBMCwriteSlices(stack, "base");
  vol = IBMCcopyStack2MRI(stack, NULL);
  MRIwrite(vol,"vol.nii");

  beta[3] = 2;
  for(nthslice=0; nthslice < stack->mri->depth; nthslice++)
    IBMCsetSliceBeta(stack->slice[nthslice], beta);
  IBMCsynthStack(targ, stack);
  IBMCwriteSlices(stack, "base2");

  vol = IBMCcopyStack2MRI(stack, NULL);
  MRIwrite(vol,"vol2.nii");

  return 0;
}

/*------------------------------------------------------------------*/
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
    else if (!strcasecmp(option, "--debug"))  debug = 1;
    else if (!strcasecmp(option, "--diag")) Gdiag_no = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--targ")) {
      if(nargc < 1) CMDargNErr(option,1);
      printf("Reading targ\n");
      targ = MRIread(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mov")) {
      if(nargc < 2) CMDargNErr(option,2);
      printf("Reading mov\n");
      mov = MRIread(pargv[0]);
      printf("Reading Rmov\n");
      Rmov = regio_read_registermat(pargv[1]);
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--alpha")) {
      if(nargc < 1) CMDargNErr(option,1);
      AlphaFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--print-params")) {
      IBMC_PARAMS *p;
      if(nargc < 1) CMDargNErr(option,1);
      p = IBMCreadParams(pargv[0]);
      IBMCprintParams(stdout,p);
      nargsused = 1;
      exit(0);
    } 
    else if (!strcasecmp(option, "--synth")) {
      MRI *SynthSrc, *SynthTemp, *SynthVol;
      MATRIX *SynthReg;
      IBMC_STACK *SynthStack;
      if(nargc < 5) CMDargNErr(option,5);
      SynthSrc = MRIread(pargv[0]);
      SynthTemp = MRIread(pargv[1]);
      SynthReg  = regio_read_registermat(pargv[2]);      
      SynthStack = IBMCinitStack(SynthTemp, SynthReg);
      SynthStack->params = IBMCreadParams(pargv[3]);
      IBMCalpha2Beta(SynthStack);
      IBMCsynthStack(SynthSrc,SynthStack);
      IBMCwriteSlices(SynthStack, pargv[4]);
      SynthVol = IBMCcopyStack2MRI(SynthStack, NULL);
      sprintf(tmpstr,"%s.nii",pargv[4]);
      MRIwrite(SynthVol,tmpstr);
      //sprintf(tmpstr,"%s.ibmc",pargv[4]);
      //IBMCwriteParams(tmpstr,SynthStack->params);
      nargsused = 5;
      exit(0);
    } 
    else if (!strcasecmp(option, "--synth-params")) {
      if(nargc < 1) CMDargNErr(option,1);
      IBMC_PARAMS *p;
      int k;
      p= IBMCallocParams(30*6);
      p->a2bmethod = IBMC_A2B_EQUAL;
      for(k=0;k<180;k+=6) p->alpha[k+3] = 30.0*(k-90)/180;
      //for(s=0;s<30;s++) for(k=0;k<6;k++) p->alpha[k+s*6] = k+1;
      IBMCwriteParams(pargv[0],p);
      nargsused = 1;
      exit(0);
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
/*------------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*------------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --v1 vol1 <regfile>: template volume \n");
  printf("   --v2 vol2 <regfile>: template volume \n");
  printf("   --v3 vol3 <regfile>: template volume \n");
  printf("\n");
  printf("  --nmax nmax   : max number of powell iterations (def 36)\n");
  printf("  --tol   tol   : powell inter-iteration tolerance on cost\n");
  printf("       This is the fraction of the cost that the difference in \n");
  printf("       successive costs must drop below to stop the optimization.  \n");
  printf("  --tol1d tol1d : tolerance on powell 1d minimizations\n");
  printf("\n");
  printf("  --synth source template src2temp.tkreg param output\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*------------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/*------------------------------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*------------------------------------------------------------------*/
static void check_options(void) {
  return;
}
/*------------------------------------------------------------------*/
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
/*-----------------------------------------------------*/
MATRIX *MRIangles2RotMatB(double *angles, MATRIX *R)
{
  double gamma, beta, alpha;
  int r,c;
  MATRIX *R3, *Rx, *Ry, *Rz;

  gamma = angles[0];
  beta  = angles[1];
  alpha = angles[2];

  //printf("angles %g %g %g\n",angles[0],angles[1],angles[2]);

  Rx = MatrixZero(3,3,NULL);
  Rx->rptr[1][1] = +1;
  Rx->rptr[2][2] = +cos(gamma);
  Rx->rptr[2][3] = -sin(gamma);
  Rx->rptr[3][2] = +sin(gamma);
  Rx->rptr[3][3] = +cos(gamma);
  //printf("Rx ----------------\n");
  //MatrixPrint(stdout,Rx);

  Ry = MatrixZero(3,3,NULL);
  Ry->rptr[1][1] = +cos(beta);
  Ry->rptr[1][3] = +sin(beta);
  Ry->rptr[2][2] = 1;
  Ry->rptr[3][1] = -sin(beta);
  Ry->rptr[3][3] = +cos(beta);
  //printf("Ry ----------------\n");
  //MatrixPrint(stdout,Ry);

  Rz = MatrixZero(3,3,NULL);
  Rz->rptr[1][1] = +cos(alpha);
  Rz->rptr[1][2] = -sin(alpha);
  Rz->rptr[2][1] = +sin(alpha);
  Rz->rptr[2][2] = +cos(alpha);
  Rz->rptr[3][3] = +1;
  //printf("Rz ----------------\n");
  //MatrixPrint(stdout,Rz);

  // This will be a 3x3 matrix
  R3 = MatrixMultiply(Rz,Ry,NULL);
  R3 = MatrixMultiply(R3,Rx,R3);

  // Stuff 3x3 into a 4x4 matrix, with (4,4) = 1
  R = MatrixZero(4,4,R);
  for(c=1; c <= 3; c++){
    for(r=1; r <= 3; r++){
      R->rptr[r][c] = R3->rptr[r][c];
    }
  }
  R->rptr[4][4] = 1;

  MatrixFree(&Rx);
  MatrixFree(&Ry);
  MatrixFree(&Rz);
  MatrixFree(&R3);

  //printf("R ----------------\n");
  //MatrixPrint(stdout,R);

  return(R);
}

/*--------------------------------------------------------*/
// Not used
MRI *MRIextractSlice(MRI *vol, MRI *slice, int Sno)
{
  int width, height, w, h, f;
  MATRIX *Sv,*P0, *crs;
  double v;

  if(Sno >= vol->depth){
    printf("ERROR: MRIextractSlice: slice no %d too large (%d)\n",
	   Sno,vol->depth);
    return(NULL);
  }

  width = vol->width;
  height = vol->height;
  if(slice == NULL){
    slice = MRIallocSequence(width,height, 1, MRI_FLOAT, vol->nframes);
    MRIcopyHeader(vol,slice);
  }
  else {
    if(slice->width != width || slice->height != height) {
      printf("ERROR: MRIextractSlice: size mismatch\n");
      return(NULL);
    }
  }

  for(f=0; f < vol->nframes; f++){
    for(w=0; w < width; w++){
      for(h=0; h < height; h++){
	v = MRIgetVoxVal(vol,w,h,Sno,f);
	MRIsetVoxVal(slice,w,h,0,f,v);
      }
    }
  }

  crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[1][3] = Sno;
  crs->rptr[1][4] = 1;
  Sv = MRIxfmCRS2XYZ(vol, 0);
  P0 = MatrixMultiply(Sv,crs,NULL);
  MRIp0ToCRAS(slice, P0->rptr[1][1], P0->rptr[1][2], P0->rptr[1][3]);
  MatrixFree(&crs);
  MatrixFree(&Sv);
  MatrixFree(&P0);
  
  return(slice);
}

