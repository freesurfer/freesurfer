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
 *    $Date: 2011/08/16 22:24:02 $
 *    $Revision: 1.13 $
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
#include <sys/stat.h>
#include <sys/types.h>

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

static char vcid[] = "$Id: mri_ibmc.c,v 1.13 2011/08/16 22:24:02 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;
MATRIX *MRIangles2RotMatB(double *angles, MATRIX *R);

#define IBMC_NL_MAX 500

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
  MATRIX *Ts, *iTs;  // tkr vox2ras (and inv)
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
  char *subjname;
} IBMC_STACK;
/*---------------------------------------------------------------*/
typedef struct 
{
  IBMC_SLICE *sliceA;
  IBMC_SLICE *sliceB;
  MATRIX *M; // M = inv(TsB)*RsB*inv(RsA)*TsA
  int nL;
  double colA[IBMC_NL_MAX], rowA[IBMC_NL_MAX];
  double colB[IBMC_NL_MAX], rowB[IBMC_NL_MAX];
  double IA[IBMC_NL_MAX], IB[IBMC_NL_MAX];
  double cost;
  double betaAprev[6],betaBprev[6]; // for caching values at previous step
  int InitNeeded;
} IBMC_PAIR;
/*---------------------------------------------------------------*/
typedef struct 
{
  int nstacks; // Number of input volume stacks
  IBMC_STACK *stack[20]; // Input stacks
  int dof; // Total number of parameters to opt = sum(nalpha)-6
  int npairs;
  IBMC_PAIR **pair;
  double cost;
  int DoSmooth;
} IBMC;

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
IBMC_PAIR *IBMCallocPair(void)
{
  IBMC_PAIR *pair;
  pair = (IBMC_PAIR *) calloc(sizeof(IBMC_PAIR),1);
  pair->M = MatrixIdentity(4,NULL);
  pair->InitNeeded = 1;
  return(pair);
}
/*---------------------------------------------------------------*/
IBMC_SLICE *IBMCinitSlice(IBMC_STACK *stack, int Sno)
{
  IBMC_SLICE *slice;

  slice = (IBMC_SLICE *) calloc(1,sizeof(IBMC_SLICE));
  slice->Sno = Sno;
  slice->Vno = stack->Vno;
  slice->mri = MRIcrop(stack->mri, 0,0,Sno, stack->mri->width-1,stack->mri->height-1,Sno);
  slice->Ts  = MRIxfmCRS2XYZtkreg(slice->mri);
  slice->iTs  = MatrixInverse(slice->Ts,NULL);
  slice->D = MatrixIdentity(4,NULL);
  slice->D->rptr[3][4] = -Sno;
  slice->iD   = MatrixInverse(slice->D,NULL);
  // Cache TsDiTv = Ts*D*inv(Tv)
  slice->TsDiTv  = MatrixMultiply(slice->Ts,slice->D,NULL);
  slice->TsDiTv  = MatrixMultiply(slice->TsDiTv,stack->iTv,slice->TsDiTv);
  // Rs0 = Ts*D*inv(Tv)*Rv = Rs when Rbeta=I
  slice->Rs0  = MatrixMultiply(slice->TsDiTv,stack->Rv,NULL);
  slice->Rbeta = MatrixIdentity(4,NULL);
  slice->Rs    = MatrixAlloc(4,4,MATRIX_REAL);
  slice->iRs   = MatrixAlloc(4,4,MATRIX_REAL);
  return(slice);
}
/*---------------------------------------------------------------*/
int IBMCinitParams(IBMC_STACK *stack, int a2bmethod)
{
  int AllocAlpha = 0;

  if(stack->params==NULL){
    stack->params = (IBMC_PARAMS *) calloc(sizeof(IBMC_PARAMS),1);
    AllocAlpha = 1;
  }
 
  stack->params->a2bmethod = a2bmethod;
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
  if(AllocAlpha) 
    stack->params->alpha = (double *) calloc(sizeof(double),stack->params->nalpha);

 return(0);
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
  //fprintf(stderr,"IBMCreadParams: %s %d %d %d\n",fname,a2bmethod,nalpha,IsASCII);
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
      printf("ERROR: %s tried to read %d, actually read %d\n",fname,nalpha,(int)nitems);
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
int IBMCwriteParamsAscii(char *fname, IBMC_PARAMS *p)
{
  FILE *fp;

  fp = fopen(fname,"w");
  if(fp == NULL){
    printf("ERROR: could not open %s\n",fname);
    return(1);
  }
  IBMCprintParams(fp, p);
  fclose(fp);
  return(0);
}


/*---------------------------------------------------------------*/
IBMC_STACK *IBMCinitStack(MRI *vol, MATRIX *Rv, int Vno)
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
  stack->Vno = Vno;
  stack->params = NULL;

  for(nthslice = 0; nthslice < stack->mri->depth; nthslice++){
    stack->slice[nthslice] = IBMCinitSlice(stack, nthslice);
    //printf("%2d %3d\n",nthslice,stack->slice[0]->vol->depth);
  }
  return(stack);
}

/*---------------------------------------------------------------*/
int IBMCwriteSlices(IBMC_STACK *stack, char *base, char *subject)
{
  char tmpstr[2000];
  int nthslice;

  for(nthslice = 0; nthslice < stack->mri->depth; nthslice++){
    sprintf(tmpstr,"%s.%03d.nii",base,nthslice);
    //printf("%2d %3d %s\n",nthslice,stack->slice[nthslice]->mri->depth,tmpstr);
    MRIwrite(stack->slice[nthslice]->mri,tmpstr);
    sprintf(tmpstr,"%s.%03d.Rs0.dat",base,nthslice);
    regio_write_register(tmpstr, subject, stack->slice[nthslice]->mri->xsize,
			 stack->slice[nthslice]->mri->zsize, .15, 
			 stack->slice[nthslice]->Rs0,FLT2INT_ROUND);
    sprintf(tmpstr,"%s.%03d.Rs.dat",base,nthslice);
    regio_write_register(tmpstr, subject, stack->slice[nthslice]->mri->xsize,
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
	  MRIsetVoxVal(mri,c,r,s,f,v);
	}
      }
    }
  }

  return(mri);
}
/*---------------------------------------------------------------*/
int IBMCalpha2Beta(IBMC *ibmc, int nthstack)
{
  int s,k,nthalpha;
  static double angles[3];
  IBMC_STACK *stack;
  IBMC_SLICE *slice;

  stack = ibmc->stack[nthstack];

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
/*------------------------------------------------------------*/
int IBMCpairCost(IBMC *ibmc, int nthpair)
{
  extern int ForceUpdate;
  IBMC_SLICE *sA, *sB;
  int nth,k,UpdateNeeded;
  double M31, M32, M34, M11, M12, M14, M21, M22, M24, c;
  float valvect=0;
  static double colA[IBMC_NL_MAX], rowA[IBMC_NL_MAX];
  static double colB[IBMC_NL_MAX], rowB[IBMC_NL_MAX];
  static double diff[IBMC_NL_MAX],diffsm[IBMC_NL_MAX];
  int nL,mth;
  IBMC_PAIR *pair;

  pair = ibmc->pair[nthpair];

  sA = pair->sliceA;
  sB = pair->sliceB;

  if(pair->InitNeeded){
    for(k=0; k<6; k++) {
      pair->betaAprev[k] = sA->beta[k];
      pair->betaBprev[k] = sB->beta[k];
    }
  }

  UpdateNeeded = 0;
  for(k=0; k<6; k++){
    if(pair->betaAprev[k] != sA->beta[k]) UpdateNeeded = 1;
    if(pair->betaBprev[k] != sB->beta[k]) UpdateNeeded = 1;
  }
  if(!UpdateNeeded  && !ForceUpdate && !pair->InitNeeded) return(0);
  pair->InitNeeded=0;

  pair->M = MatrixMultiply(sB->iTs,sB->Rs, pair->M);
  MatrixMultiply(pair->M, sA->iRs, pair->M);
  MatrixMultiply(pair->M, sA->Ts, pair->M);
  
  #if 0
  printf("iTsB \n");
  MatrixPrint(stdout,sB->iTs);
  printf("RsB \n");
  MatrixPrint(stdout,sB->Rs);
  printf("iRsA \n");
  MatrixPrint(stdout,sA->iRs);
  printf("TsA \n");
  MatrixPrint(stdout,sA->Ts);
  printf("M \n");
  MatrixPrint(stdout,pair->M);
  #endif

  M31 = pair->M->rptr[3][1];
  M32 = pair->M->rptr[3][2];
  M34 = pair->M->rptr[3][4];

  if(fabs(M32) > fabs(M31)){
    nL = sA->mri->width;
    for(nth=0; nth < nL; nth++) {
      colA[nth] = nth; //+0.5; // +0.5 force interp
      rowA[nth] = -(M34+M31*colA[nth])/M32;
    }
  }
  else {
    nL = sA->mri->height;
    for(nth=0; nth < nL; nth++) {
      rowA[nth] = nth; //+0.5; // +0.5 force interp
      colA[nth] = -(M34+M32*rowA[nth])/M31;
    }
  }

  M11 = pair->M->rptr[1][1];
  M12 = pair->M->rptr[1][2];
  M14 = pair->M->rptr[1][4];

  M21 = pair->M->rptr[2][1];
  M22 = pair->M->rptr[2][2];
  M24 = pair->M->rptr[2][4];

  mth = 0;
  for(nth=0; nth < nL; nth++) {
    colB[nth] = M11*colA[nth] + M12*rowA[nth] + M14;
    rowB[nth] = M21*colA[nth] + M22*rowA[nth] + M24;

    if(colA[nth] < 0 || colA[nth] > sA->mri->width-1) continue;
    if(rowA[nth] < 0 || rowA[nth] > sA->mri->height-1) continue;
    if(colB[nth] < 0 || colB[nth] > sB->mri->width-1) continue;
    if(rowB[nth] < 0 || rowB[nth] > sB->mri->height-1) continue;
    pair->colA[mth] = colA[nth];
    pair->rowA[mth] = rowA[nth];
    pair->colB[mth] = colB[nth];
    pair->rowB[mth] = rowB[nth];

    MRIsampleSeqVolume(sA->mri, colA[nth], rowA[nth], 0, &valvect, 0, 0);
    pair->IA[mth] = valvect;
    //pair->IA[nth] = MRIgetVoxVal(sA->mri,,irsA,p->SnoA,0);

    MRIsampleSeqVolume(sB->mri, colB[nth], rowB[nth], 0, &valvect, 0, 0);
    pair->IB[mth] = valvect;

    diff[mth] = pair->IA[mth] - pair->IB[mth];

    //printf("   %2d   %2d %2d   %2d %2d   %6.4f %6.4f   %6.4f\n",
    //   mth, (int)round(pair->colA[mth]), (int) round(pair->rowA[mth]), 
    //   (int)round(pair->colB[mth]),(int)round(pair->rowB[mth]),
    //   pair->IA[mth],pair->IB[mth],c);
    mth ++;
  }
  pair->nL = mth;

  if(ibmc->DoSmooth){
    for(mth=0; mth < pair->nL; mth++) {
      if(mth==0 || mth == pair->nL-1) {
	diffsm[mth] = diff[mth];
	continue;
      }
      diffsm[mth] = (diff[mth-1] + diff[mth] + diff[mth+1])/3.0;
    }
  }
  else {
    for(mth=0; mth < pair->nL; mth++) diffsm[mth] = diff[mth];
  }

  pair->cost = 0;
  for(mth=0; mth < pair->nL; mth++) {
    c = diffsm[mth] * diffsm[mth];
    pair->cost += c;
  }

  for(k=0; k<6; k++) {
    pair->betaAprev[k] = sA->beta[k];
    pair->betaBprev[k] = sB->beta[k];
  }

  return(0);
}
/*---------------------------------------------------------------*/
double IBMCcost(IBMC *ibmc)
{
  int k, nthpair;
  double cost;

  for(k=0; k < ibmc->nstacks; k++)  IBMCalpha2Beta(ibmc,k);

  cost=0;
  for(nthpair = 0; nthpair < ibmc->npairs; nthpair++){
    IBMCpairCost(ibmc,nthpair);
    cost += ibmc->pair[nthpair]->cost;
  }
  cost /= ibmc->npairs;
  ibmc->cost = cost;
  return(cost);
}

/*---------------------------------------------------------------*/
int IBMCinit(IBMC *ibmc, int a2bmethod)
{
  int k, nthpair, s0, s1, s2;

  ibmc->dof = 0;
  for(k=0; k < ibmc->nstacks; k++){
    IBMCinitParams(ibmc->stack[k],a2bmethod);
    ibmc->dof += ibmc->stack[k]->params->nalpha;
  }
  ibmc->dof -= 6; 
  printf("dof = %d\n",ibmc->dof); 

  ibmc->npairs = 
    ibmc->stack[0]->mri->depth * ibmc->stack[1]->mri->depth +
    ibmc->stack[0]->mri->depth * ibmc->stack[2]->mri->depth +
    ibmc->stack[1]->mri->depth * ibmc->stack[2]->mri->depth;
  printf("npairs = %d\n",ibmc->npairs);

  ibmc->pair = (IBMC_PAIR **) calloc(ibmc->npairs,sizeof(IBMC_PAIR *));
  if(ibmc->pair == NULL){
    printf("ERROR: could not alloc %d pairs\n",ibmc->npairs);
    return(1);
  }

  printf("Initializing %d pairs\n",ibmc->npairs);
  nthpair = 0;
  for(s0=0; s0 < ibmc->stack[0]->mri->depth; s0++){
    for(s1=0; s1 < ibmc->stack[1]->mri->depth; s1++){
      ibmc->pair[nthpair] = IBMCallocPair();
      ibmc->pair[nthpair]->sliceA = ibmc->stack[0]->slice[s0];
      ibmc->pair[nthpair]->sliceB = ibmc->stack[1]->slice[s1];
      nthpair++;
    }
  }
  for(s0=0; s0 < ibmc->stack[0]->mri->depth; s0++){
    for(s2=0; s2 < ibmc->stack[2]->mri->depth; s2++){
      ibmc->pair[nthpair] = IBMCallocPair();
      ibmc->pair[nthpair]->sliceA = ibmc->stack[0]->slice[s0];
      ibmc->pair[nthpair]->sliceB = ibmc->stack[2]->slice[s2];
      nthpair++;
    }
  }
  for(s1=0; s1 < ibmc->stack[1]->mri->depth; s1++){
    for(s2=0; s2 < ibmc->stack[2]->mri->depth; s2++){
      ibmc->pair[nthpair] = IBMCallocPair();
      ibmc->pair[nthpair]->sliceA = ibmc->stack[1]->slice[s1];
      ibmc->pair[nthpair]->sliceB = ibmc->stack[2]->slice[s2];
      nthpair++;
    }
  }
  printf("Done initializing pairs %d\n",nthpair);
  return(0);
}
/*---------------------------------------------------------------*/
int IBMCsetParam(IBMC *ibmc, int nthparam, double paramval);
int IBMCzeroParams(IBMC *ibmc);
int IBMCprofile(IBMC *ibmc, char *ProfileFile);
double *IBMCgetParams(IBMC *ibmc, double *params);
int IBMCsetParams(IBMC *ibmc, double *params);
int IBMClineMin(IBMC *ibmc);
int IBMCwriteStackReg(IBMC_STACK *stack, char *fname, char *subject);
float compute_powell_cost(float *params);
int MinPowell(double ftol, double linmintol, int nmaxiters);
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
MRI *targ=NULL, *mov=NULL;
MATRIX *Rtarg=NULL,*Rmov=NULL;
char *AlphaFile=NULL;
char tmpstr[2000];
IBMC *ibmc;
int a2bmethod = IBMC_A2B_RIGID;
int nCostEvaluations;
int nMaxItersPowell = 36;
double TolPowell = 1e-8; //1e-8;
double LinMinTolPowell = 1e-8; //1e-8;
int DoProfile=0;
char *ProfileFile=NULL;
int ProfileParam=0;
char *outdir = NULL;
char *ParBase=NULL;
int ForceUpdate=0;
int DoLineMin = 0;
char *subject = "subject-unkown";

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,k,err;
  double p;
  FILE *fp;

  ibmc = (IBMC *) calloc(sizeof(IBMC),1);
  ibmc->DoSmooth = 0;

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

  IBMCinit(ibmc,a2bmethod);

  if(DoProfile){
    printf("Profiling %d %s\n",ProfileParam,ProfileFile);
    IBMCprofile(ibmc, ProfileFile);
    exit(0);
    IBMCzeroParams(ibmc);
    fp = fopen(ProfileFile,"w");
    for(p = -5; p < 5; p+= .01){
      IBMCsetParam(ibmc, ProfileParam, p);
      IBMCcost(ibmc);
      printf("%6.4lf %10.5lf\n",p,ibmc->cost);
      fprintf(fp,"%6.4lf %10.5lf\n",p,ibmc->cost);
    }
    fclose(fp);
    exit(0);
  }

  printf("Creating output directory %s\n",outdir);
  err = mkdir(outdir,0777);
  if (err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",outdir);
    perror(NULL);
    return(1);
  }

  IBMCcost(ibmc);
  printf("Initial cost %g\n",ibmc->cost);

  for(k=0; k < ibmc->nstacks; k++){
    sprintf(tmpstr,"%s/init.%s%d.ibmc",outdir,ParBase,k);
    IBMCwriteParamsAscii(tmpstr,ibmc->stack[k]->params);
    printf("Init Stack %d ------\n",k);
    IBMCprintParams(stdout,ibmc->stack[k]->params);
  }

  if(DoLineMin){
    IBMClineMin(ibmc);
    printf("Params at min after line search\n");
    for(k=0; k < ibmc->nstacks; k++){
      sprintf(tmpstr,"%s/linmin.%s%d.ibmc",outdir,ParBase,k);
      IBMCwriteParamsAscii(tmpstr,ibmc->stack[k]->params);
      printf("Stack %d ------\n",k);
      IBMCprintParams(stdout,ibmc->stack[k]->params);
    }
    printf("\n");
  }


  MinPowell(TolPowell, LinMinTolPowell, nMaxItersPowell);

  for(k=0; k < ibmc->nstacks; k++){
    sprintf(tmpstr,"%s/%s%d.ibmc",outdir,ParBase,k);
    IBMCwriteParamsAscii(tmpstr,ibmc->stack[k]->params);
    printf("Stack %d ------\n",k);
    IBMCprintParams(stdout,ibmc->stack[k]->params);
    sprintf(tmpstr,"%s/%s%d.reg.dat",outdir,ParBase,k);
    IBMCwriteStackReg(ibmc->stack[k], tmpstr, subject);

    if(a2bmethod == IBMC_A2B_EQUAL){
      sprintf(tmpstr,"%s/stack%d",outdir,k);      
      err = mkdir(tmpstr,0777);
      sprintf(tmpstr,"%s/stack%d/v",outdir,k);      
      IBMCwriteSlices(ibmc->stack[k], tmpstr, subject);
    }

  }

  printf("mri_ibmc done\n");

  exit(0);
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
    else if (!strcasecmp(option, "--force-update"))  ForceUpdate = 1;
    else if (!strcasecmp(option, "--diag")) Gdiag_no = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--a2b-rigid")) a2bmethod = IBMC_A2B_RIGID;
    else if (!strcasecmp(option, "--a2b-equal")) a2bmethod = IBMC_A2B_EQUAL;
    else if (!strcasecmp(option, "--line-min"))  DoLineMin = 1;
    else if (!strcasecmp(option, "--smooth"))  ibmc->DoSmooth = 1;
    else if (!strcasecmp(option, "--low-tol")){
      TolPowell = 1e-1;
      LinMinTolPowell = 1e-1;
    }

    else if (!strcasecmp(option, "--o")) {
      if(nargc < 1) CMDargNErr(option,1);
      outdir = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--par")) {
      if(nargc < 1) CMDargNErr(option,1);
      ParBase = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--s")) {
      if(nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--profile")) {
      if(nargc < 2) CMDargNErr(option,2);
      DoProfile = 1;
      ProfileFile = pargv[0];
      sscanf(pargv[1],"%d",&ProfileParam);
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--stack")) {
      if(nargc < 2) CMDargNErr(option,2);
      MRI *mri;    // dont free
      MATRIX *reg; // dont free
      printf("Reading stack %d  %s\n",ibmc->nstacks,pargv[0]);
      mri = MRIread(pargv[0]);
      reg = regio_read_registermat(pargv[1]);
      ibmc->stack[ibmc->nstacks] = IBMCinitStack(mri, reg,ibmc->nstacks);
      nargsused = 2;
      if(CMDnthIsArg(nargc, pargv, 2)) {
        ibmc->stack[ibmc->nstacks]->params = IBMCreadParams(pargv[2]);
	printf("Loaded Stack %d ------\n",ibmc->nstacks);
	IBMCprintParams(stdout,ibmc->stack[ibmc->nstacks]->params);
        nargsused ++;
      }
      ibmc->nstacks++;
    } 
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
      SynthStack = IBMCinitStack(SynthTemp, SynthReg, 0);
      SynthStack->params = IBMCreadParams(pargv[3]);
      ibmc->stack[0] = SynthStack;
      IBMCalpha2Beta(ibmc,0);
      IBMCsynthStack(SynthSrc,SynthStack);
      IBMCwriteSlices(SynthStack, pargv[4],subject);
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
  if(!DoProfile){
    if(outdir==NULL){
      printf("ERROR: must spec outdir\n");
      exit(1);
    }
    if(ParBase==NULL){
      printf("ERROR: must spec par base\n");
      exit(1);
    }
  }
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

/*--------------------------------------------------------*/
float compute_powell_cost(float *params)
{
  extern IBMC *ibmc;
  extern int nCostEvaluations;
  extern char *ParBase;
  static double copt = 10e10;
  static int first=1;
  static float *paramsprev, *beta, cost0;
  static struct timeb timer;
  double secCostTime=0;
  float cost;
  int newopt,k,kbeta=0,n,r,c,r0;

  if(first){
    paramsprev = vector(1,ibmc->dof);
    beta = (float *) calloc(ibmc->dof,sizeof(float));
    TimerStart(&timer);
    kbeta=1;
  }
  else {
    kbeta=0;
    for(k=0; k < ibmc->dof; k++) {
      if(paramsprev[k+1] != params[k+1]){
	kbeta = k+1;
	break;
      }
    }
  }

  n = 0;
  for(c=0; c < ibmc->nstacks; c++) {
    r0 = 0;
    if(c==0) r0 = 6;
    for(r=r0; r < ibmc->stack[c]->params->nalpha; r++) {
      ibmc->stack[c]->params->alpha[r] = params[n+1];
      n++;
    }
  }

  cost = IBMCcost(ibmc);
  nCostEvaluations ++;
  newopt = 0;
  if(copt >= cost){
    copt = cost;
    newopt = 1;
  }
  if(first) {
    cost0 = cost;
    if(cost0 < FLT_MIN) cost0 = 1;
    printf("Initial cost %g\n",cost);
  }

  if(newopt && kbeta != 0){
    secCostTime = TimerStop(&timer)/1000.0;
    printf("%4d %3d  %9.6f    %9.6f %9.6f   t=%7.3f\n",
	   nCostEvaluations,kbeta-1,params[kbeta],cost,copt,
	   secCostTime/60.0);
    fflush(stdout);
    for(k=0; k < ibmc->nstacks; k++){
      sprintf(tmpstr,"%s/curopt.%s%d.ibmc",outdir,ParBase,k);
      IBMCwriteParamsAscii(tmpstr,ibmc->stack[k]->params);
    }
  }

  for(k=0; k < ibmc->dof; k++) paramsprev[k+1] = params[k+1];

  first=0;
  return(cost);
}
/*---------------------------------------------------------*/
int MinPowell(double ftol, double linmintol, int nmaxiters)
{
  float **xi, fret;
  int    r, r0,c, n;
  int niters,err;
  float *pPowel;

  xi = matrix(1, ibmc->dof, 1, ibmc->dof) ;
  for (r = 1 ; r <= ibmc->dof ; r++) {
    for (c = 1 ; c <= ibmc->dof ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  pPowel = vector(1, ibmc->dof) ;
  n = 0;
  for(c=0; c < ibmc->nstacks; c++) {
    r0 = 0;
    if(c==0) r0 = 6;
    for(r=r0; r < ibmc->stack[c]->params->nalpha; r++) {
      pPowel[n+1] = ibmc->stack[c]->params->alpha[r];
      n++;
    }
  }

  printf("Starting Powell nitersmax = %d, ftol=%g  linmintol %g, \n",
	 nmaxiters,ftol,linmintol);
  err=OpenPowell2(pPowel, xi, ibmc->dof, ftol, linmintol, nmaxiters, 
	      &niters, &fret, compute_powell_cost);
  printf("Powell done niters = %d, err=%d, fret = %f\n",niters,err,fret);

  // Stuff params back into IBMC
  n = 0;
  for(c=0; c < ibmc->nstacks; c++) {
    r0 = 0;
    if(c==0) r0 = 6;
    for(r=r0; r < ibmc->stack[c]->params->nalpha; r++) {
      ibmc->stack[c]->params->alpha[r] = pPowel[n+1];
      n++;
    }
  }
  // Update IBMC and compute final cost
  IBMCcost(ibmc);

  free_matrix(xi, 1, ibmc->dof, 1, ibmc->dof);
  return(niters);
}
/*---------------------------------------------------------*/
int IBMCsetParam(IBMC *ibmc, int nthparam, double paramval)
{
  int  r, c, n;
  n = 0;
  for(c=0; c < ibmc->nstacks; c++) {
    for(r=0; r < ibmc->stack[c]->params->nalpha; r++) {
      if(n == nthparam) ibmc->stack[c]->params->alpha[r] = paramval;
      n++;
    }
  }
  return(0);
}
/*---------------------------------------------------------*/
int IBMCzeroParams(IBMC *ibmc)
{
  int  s, p;
  for(s=0; s < ibmc->nstacks; s++) {
    for(p=0; p < ibmc->stack[s]->params->nalpha; p++)
      ibmc->stack[s]->params->alpha[p] = 0;
  }
  return(0);
}

/*---------------------------------------------------------*/
int IBMCprofile(IBMC *ibmc, char *ProfileFile)
{
  int  nthp, nv, nthv;
  double v, vmin, vmax, dv;
  MATRIX *C;
  FILE *fp;

  vmin = -3;
  vmax = +3;
  dv = .05;
  nv = round((vmax-vmin)/dv);

  C = MatrixAlloc(nv,ibmc->dof,MATRIX_REAL);

  printf("Starting Profile\n");
  nthp = 0;
  for(nthp = 0; nthp < ibmc->dof; nthp++){
    printf("Param %3d/%d\n",nthp+1,ibmc->dof);
    IBMCzeroParams(ibmc);
    for(nthv=0; nthv < nv; nthv++){
      v = vmin + nthv*dv;
      IBMCsetParam(ibmc, nthp, v);
      IBMCcost(ibmc);
      C->rptr[nthv+1][nthp+1] = ibmc->cost;
    }
  }

  fp = fopen(ProfileFile,"w");
  for(nthv=0; nthv < nv; nthv++){
    v = vmin + nthv*dv;
    fprintf(fp,"%6.3lf ",v);
    for(nthp = 0; nthp < ibmc->dof; nthp++)
      fprintf(fp,"%6.2f ",C->rptr[nthv+1][nthp+1]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  MatrixFree(&C);
  return(0);
}
/*---------------------------------------------------------------------*/
double *IBMCgetParams(IBMC *ibmc, double *params)
{
  int  nths, ntha, nthp;

  if(params == NULL) params = (double *) calloc(sizeof(double),ibmc->dof);

  nthp = 0;
  for(nths=0; nths < ibmc->nstacks; nths++) {
    for(ntha=0; ntha < ibmc->stack[nths]->params->nalpha; ntha++){
      params[nthp] = ibmc->stack[nths]->params->alpha[ntha];
      nthp ++;
    }
  }
  return(params);
}
/*---------------------------------------------------------------------*/
int IBMCsetParams(IBMC *ibmc, double *params)
{
  int  nths, ntha, nthp;

  nthp = 0;
  for(nths=0; nths < ibmc->nstacks; nths++) {
    for(ntha=0; ntha < ibmc->stack[nths]->params->nalpha; ntha++){
      ibmc->stack[nths]->params->alpha[ntha] = params[nthp];
      nthp ++;
    }
  }
  return(0);
}
/*---------------------------------------------------------------------*/
int IBMClineMin(IBMC *ibmc)
{
  int  nthp, nv, nthv;
  double v, vmin, vmax, dv, cmin;
  double *optparams=NULL;

  vmin = -5;
  vmax = +5;
  dv = .5;
  nv = round((vmax-vmin)/dv);

  printf("Starting Line Min init cost = %g\n",ibmc->cost);
  nthp = 0;
  cmin=ibmc->cost;
  optparams = IBMCgetParams(ibmc,optparams);
  for(nthp = 0; nthp < ibmc->dof; nthp++){
    printf("  Param %3d/%d\n",nthp+1,ibmc->dof);
    fflush(stdout);
    IBMCsetParams(ibmc,optparams);    
    for(nthv=0; nthv < nv; nthv++){
      v = vmin + nthv*dv;
      IBMCsetParam(ibmc, nthp+6, v);
      IBMCcost(ibmc);
      if(cmin > ibmc->cost){
	cmin = ibmc->cost;
	optparams = IBMCgetParams(ibmc,optparams);
	printf("     %g\n",ibmc->cost);
	fflush(stdout);
      }	
    }
  }
  IBMCsetParams(ibmc,optparams);
  IBMCcost(ibmc);
  printf("Ending Line Min init cost = %g\n",ibmc->cost);

  return(0);
}
/*---------------------------------------------------------------------*/
int IBMCwriteStackReg(IBMC_STACK *stack, char *fname, char *subject)
{
  MATRIX *R;
  // Really only useful for RIGID, slice number does not matter with RIGID
  R = MatrixMultiply(stack->slice[0]->Rbeta,stack->Rv,NULL);
  regio_write_register(fname, subject, stack->slice[0]->mri->xsize,
		       stack->slice[0]->mri->zsize, .15, R,FLT2INT_ROUND);
  MatrixFree(&R);
  return(0);
}





