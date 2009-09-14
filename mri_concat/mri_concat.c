/**
 * @file  mri_concat.c
 * @brief Concatenates input data sets.
 *
 * EXAMPLES:
 *   mri_concat --i f1.mgh --i f2.mgh --o cout.mgh
 *   mri_concat f1.mgh f2.mgh --o cout.mgh
 *   mri_concat f*.mgh --o cout.mgh
 *   mri_concat f*.mgh --o coutmn.mgh --mean
 *   mri_concat f*.mgh --o coutdiff.mgh --paired-diff
 *   mri_concat f*.mgh --o coutdiff.mgh --paired-diff-norm
 *   mri_concat f*.mgh --o coutdiff.mgh --paired-diff-norm1
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2009/09/14 16:56:29 $
 *    $Revision: 1.37 $
 *
 * Copyright (C) 2002-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "fmriutils.h"
#include "version.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
//static int  singledash(char *flag);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_concat.c,v 1.37 2009/09/14 16:56:29 greve Exp $";
char *Progname = NULL;
int debug = 0;
char *inlist[5000];
int ninputs = 0;
char *out = NULL;
MRI *mritmp, *mritmp0, *mriout, *mask=NULL;
char *maskfile = NULL;
int DoMean=0;
int DoMeanDivN=0;
int DoSum=0;
int DoVar=0;
int DoStd=0;
int DoMax=0;
int DoMaxIndex=0;
int DoMin=0;
int DoPaired=0;
int DoPairedAvg=0;
int DoPairedSum=0;
int DoPairedDiff=0;
int DoPairedDiffNorm=0;
int DoPairedDiffNorm1=0;
int DoPairedDiffNorm2=0;
int DoASL = 0;
int DoVote=0;
int DoSort=0;
int DoCombine=0;
int DoKeepDatatype=0;

int DoMultiply=0;
double MultiplyVal=0;

int DoAdd=0;
double AddVal=0;
int DoBonfCor=0;
int DoAbs = 0;
int DoPos = 0;
int DoNeg = 0;

char *matfile = NULL;
MATRIX *M = NULL;
int ngroups = 0;
MATRIX *GroupedMeanMatrix(int ngroups, int ntotal);

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, nthin, nframestot=0, nr=0,nc=0,ns=0, fout;
  int r,c,s,f,nframes,err;
  double v, v1, v2, vavg;
  int inputDatatype=MRI_UCHAR;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if(maskfile){
    printf("Loading mask %s\n",maskfile);
    mask = MRIread(maskfile);
    if(mask == NULL) exit(1);
  }

  printf("ninputs = %d\n",ninputs);
  printf("Checking inputs\n");
  for(nthin = 0; nthin < ninputs; nthin++) {
    if(Gdiag_no > 0 || debug) {
      printf("Checking %2d %s\n",nthin,inlist[nthin]);
      fflush(stdout);
    }
    mritmp = MRIreadHeader(inlist[nthin],MRI_VOLUME_TYPE_UNKNOWN);
    if (mritmp == NULL) {
      printf("ERROR: reading %s\n",inlist[nthin]);
      exit(1);
    }
    if (nthin == 0) {
      nc = mritmp->width;
      nr = mritmp->height;
      ns = mritmp->depth;
    }
    if (mritmp->width != nc ||
        mritmp->height != nr ||
        mritmp->depth != ns) {
      printf("ERROR: dimension mismatch between %s and %s\n",
             inlist[0],inlist[nthin]);
      exit(1);
    }

    nframestot += mritmp->nframes;
    inputDatatype = mritmp->type; // used by DoKeepDatatype option
    MRIfree(&mritmp);
  }

  if(DoCombine) nframestot = 1; // combine creates a single-frame volume
  printf("nframestot = %d\n",nframestot);

  if(ngroups != 0){
    printf("Creating grouped mean matrix ngroups=%d, nper=%d\n",ngroups,nframestot/ngroups);
    M = GroupedMeanMatrix(ngroups,nframestot);
    if(M==NULL) exit(1);
    if(debug) MatrixPrint(stdout,M);
  }

  if(M != NULL){
    if(nframestot != M->cols){
      printf("ERROR: dimension mismatch between inputs (%d) and matrix (%d)\n",
	     nframestot,M->rows);
      exit(1);
    }
  }

  if (DoPaired) {
    if (remainder(nframestot,2) != 0) {
      printf("ERROR: --paired-xxx specified but there are an "
             "odd number of frames\n");
      exit(1);
    }
  }

  printf("Allocing output\n");
  fflush(stdout);
  int datatype=MRI_FLOAT;
  if (DoKeepDatatype) {
    datatype = inputDatatype;
  }
  mriout = MRIallocSequence(nc,nr,ns,datatype,nframestot);
  if (mriout == NULL) exit(1);

  fout = 0;
  for (nthin = 0; nthin < ninputs; nthin++) {
    if(Gdiag_no > 0 || debug) {
      printf("Loading %dth input %s\n",nthin+1,fio_basename(inlist[nthin],NULL));
      fflush(stdout);
    }
    mritmp = MRIread(inlist[nthin]);
    if(mritmp == NULL){
      printf("ERROR: loading %s\n",inlist[nthin]);
      exit(1);
    }
    if(nthin == 0) {
      MRIcopyHeader(mritmp, mriout);
      //mriout->nframes = nframestot;
    }
    if(DoAbs) {
      if(Gdiag_no > 0 || debug) printf("Removing sign from input\n");
      MRIabs(mritmp,mritmp);
    }
    if(DoPos) {
      if(Gdiag_no > 0 || debug) printf("Setting input negatives to 0.\n");
      MRIpos(mritmp,mritmp);
    }
    if(DoNeg) {
      if(Gdiag_no > 0 || debug) printf("Setting input positives to 0.\n");
      MRIneg(mritmp,mritmp);
    }
    for(f=0; f < mritmp->nframes; f++) {
      for(c=0; c < nc; c++) {
        for(r=0; r < nr; r++) {
          for(s=0; s < ns; s++) {
            v = MRIgetVoxVal(mritmp,c,r,s,f);
            if (DoCombine) {
              if (v > 0) MRIsetVoxVal(mriout,c,r,s,0,v);
            }
            else MRIsetVoxVal(mriout,c,r,s,fout,v);
          }
        }
      }
      fout++;
    }
    MRIfree(&mritmp);
  }

  if(DoASL){
    printf("Computing ASL matrix matrix\n");
    M = ASLinterpMatrix(mriout->nframes);
  }

  if(M != NULL){
    printf("Multiplying by matrix\n");
    mritmp = fMRImatrixMultiply(mriout, M, NULL);
    if(mritmp == NULL) exit(1);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoPaired) {
    printf("Combining pairs\n");
    mritmp = MRIcloneBySpace(mriout,-1,mriout->nframes/2);
    for (c=0; c < nc; c++) {
      for (r=0; r < nr; r++) {
        for (s=0; s < ns; s++) {
          fout = 0;
          for (f=0; f < mriout->nframes; f+=2) {
            v1 = MRIgetVoxVal(mriout,c,r,s,f);
            v2 = MRIgetVoxVal(mriout,c,r,s,f+1);
            v = 0;
            if(DoPairedAvg)  v = (v1+v2)/2.0;
            if(DoPairedSum)  v = (v1+v2);
            if(DoPairedDiff) v = v1-v2; // difference
            if(DoPairedDiffNorm) {
              v = v1-v2; // difference
              vavg = (v1+v2)/2.0;
              if (vavg != 0.0) v = v/vavg;
            }
            if(DoPairedDiffNorm1) {
              v = v1-v2; // difference
              if (v1 != 0.0) v = v/v1;
              else v = 0;
            }
            if(DoPairedDiffNorm2) {
              v = v1-v2; // difference
              if (v2 != 0.0) v = v/v2;
              else v = 0;
            }
            MRIsetVoxVal(mritmp,c,r,s,fout,v);
            fout++;
          }
        }
      }
    }
    MRIfree(&mriout);
    mriout = mritmp;
  }
  nframes = mriout->nframes;
  printf("nframes = %d\n",nframes);

  if(DoBonfCor){
    DoAdd = 1;
    AddVal = -log10(mriout->nframes);
  }

  if(DoMean) {
    printf("Computing mean across frames\n");
    mritmp = MRIframeMean(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }
  if(DoMeanDivN) {
    printf("Computing mean2 = sum/(nframes^2)\n");
    mritmp = MRIframeSum(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
    MRImultiplyConst(mriout, 1.0/(nframes*nframes), mriout);
  }
  if(DoSum) {
    printf("Computing sum across frames\n");
    mritmp = MRIframeSum(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoStd || DoVar) {
    printf("Computing std/var across frames\n");
    if(mriout->nframes < 2){
      printf("ERROR: cannot compute std from one frame\n");
      exit(1);
    }
    //mritmp = fMRIvariance(mriout, -1, 1, NULL);
    mritmp = fMRIcovariance(mriout, 0, -1, NULL, NULL);
    if(DoStd) MRIsqrt(mritmp, mritmp);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoMax) {
    printf("Computing max across all frames \n");
    mritmp = MRIvolMax(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoMaxIndex) {
    printf("Computing max index across all frames \n");
    mritmp = MRIvolMaxIndex(mriout,1,NULL,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoMin) {
    printf("Computing min across all frames \n");
    mritmp = MRIvolMin(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoSort) {
    printf("Sorting \n");
    mritmp = MRIsort(mriout,mask,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoVote) {
    printf("Voting \n");
    mritmp = MRIvote(mriout,mask,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoMultiply){
    printf("Multiplying by %lf\n",MultiplyVal);
    MRImultiplyConst(mriout, MultiplyVal, mriout);
  }

  if(DoAdd){
    printf("Adding %lf\n",AddVal);
    MRIaddConst(mriout, AddVal, mriout);
  }

  printf("Writing to %s\n",out);
  err = MRIwrite(mriout,out);
  if(err) exit(err);

  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
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
    else if (!strcasecmp(option, "--mean"))   DoMean = 1;
    else if (!strcasecmp(option, "--mean-div-n")) DoMeanDivN = 1;
    else if (!strcasecmp(option, "--mean2"))      DoMeanDivN = 1;
    else if (!strcasecmp(option, "--sum"))    DoSum = 1;
    else if (!strcasecmp(option, "--std"))    DoStd = 1;
    else if (!strcasecmp(option, "--var"))    DoVar = 1;
    else if (!strcasecmp(option, "--abs"))    DoAbs = 1;
    else if (!strcasecmp(option, "--pos"))    DoPos = 1;
    else if (!strcasecmp(option, "--neg"))    DoNeg = 1;
    else if (!strcasecmp(option, "--max"))    DoMax = 1;
    else if (!strcasecmp(option, "--max-index")) DoMaxIndex = 1;
    else if (!strcasecmp(option, "--min"))    DoMin = 1;
    else if (!strcasecmp(option, "--vote"))   DoVote = 1;
    else if (!strcasecmp(option, "--sort"))   DoSort = 1;
    else if (!strcasecmp(option, "--max-bonfcor")){
      DoMax = 1;
      DoBonfCor = 1;
    }
    else if (!strcasecmp(option, "--asl")){
      DoASL = 1;
    }
    else if (!strcasecmp(option, "--paired-avg")){
      DoPaired = 1;
      DoPairedAvg = 1;
    }
    else if (!strcasecmp(option, "--paired-sum")){
      DoPaired = 1;
      DoPairedSum = 1;
    }
    else if (!strcasecmp(option, "--paired-diff")){
      DoPaired = 1;
      DoPairedDiff = 1;
    }
    else if (!strcasecmp(option, "--paired-diff-norm")) {
      DoPairedDiff = 1;
      DoPairedDiffNorm = 1;
      DoPaired = 1;
    } else if (!strcasecmp(option, "--paired-diff-norm1")) {
      DoPairedDiff = 1;
      DoPairedDiffNorm1 = 1;
      DoPaired = 1;
    } else if (!strcasecmp(option, "--paired-diff-norm2")) {
      DoPairedDiff = 1;
      DoPairedDiffNorm2 = 1;
      DoPaired = 1;
    } else if (!strcasecmp(option, "--combine")) {
      DoCombine = 1;
    } else if (!strcasecmp(option, "--keep-datatype")) {
      DoKeepDatatype = 1;
    } else if ( !strcmp(option, "--i") ) {
      if (nargc < 1) argnerr(option,1);
      inlist[ninputs] = pargv[0];
      ninputs ++;
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--mtx") ) {
      if (nargc < 1) argnerr(option,1);
      matfile = pargv[0];
      M = MatrixReadTxt(matfile, NULL);
      if(M==NULL){
	printf("ERROR: reading %s\n",matfile);
	exit(1);
      }
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--gmean") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&ngroups);
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--o") ) {
      if (nargc < 1) argnerr(option,1);
      out = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--mask") ) {
      if (nargc < 1) argnerr(option,1);
      maskfile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--mul") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&MultiplyVal);
      DoMultiply = 1;
      nargsused = 1;
    } else if ( !strcmp(option, "--add") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&AddVal);
      DoAdd = 1;
      nargsused = 1;
    } else {
      inlist[ninputs] = option;
      ninputs ++;
      //fprintf(stderr,"ERROR: Option %s unknown\n",option);
      //if(singledash(option))
      //fprintf(stderr,"       Did you really mean -%s ?\n",option);
      //exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --i invol <--i invol ...> (don't need --i) \n");
  printf("   --o out \n");
  printf("\n");
  printf("   --paired-sum  : compute paired sum (1+2, 3d+4, etc) \n");
  printf("   --paired-avg  : compute paired avg (1+2, 3d+4, etc) \n");
  printf("   --paired-diff : compute paired diff (1-2, 3-4, etc) \n");
  printf("   --paired-diff-norm : same as paired-diff but scale by TP1,2 average \n");
  printf("   --paired-diff-norm1 : same as paired-diff but scale by TP1 \n");
  printf("   --paired-diff-norm2 : same as paired-diff but scale by TP2 \n");
  printf("   --mtx matrix.asc    : multiply by matrix in ascii file\n");
  printf("   --gmean Ng          : create matrix to average Ng groups, Nper=Ntot/Ng\n");
  printf("\n");
  printf("   --combine : combine non-zero values into a one-frame volume\n");
  printf("             (useful to combine lh.ribbon.mgz and rh.ribbon.mgz)\n");
  printf("   --keep-datatype : write output in same datatype as input\n");
  printf("                    (default is to write output in Float format)\n");
  printf("\n");
  printf("   --abs  : take abs of input\n");
  printf("   --pos  : set input negatives to 0\n");
  printf("   --neg  : set input postives to 0\n");
  printf("   --mean : compute mean of concatenated volumes\n");
  printf("   --mean-div-n : compute mean/nframes (good for var) \n");
  printf("   --sum  : compute sum of concatenated volumes\n");
  printf("   --var  : compute var  of concatenated volumes\n");
  printf("   --std  : compute std  of concatenated volumes\n");
  printf("   --max  : compute max  of concatenated volumes\n");
  printf("   --max-index  : compute index of max of concatenated volumes (1-based)\n");
  printf("   --min  : compute min of concatenated volumes\n");
  printf("   --vote : most frequent value at each voxel and fraction of occurances\n");
  printf("   --sort : sort each voxel by ascending frame value\n");
  printf("\n");
  printf("   --max-bonfcor  : compute max and bonferroni correct (assumes -log10(p))\n");
  printf("   --mul mulval   : multiply by mulval\n");
  printf("   --add addval   : add addval\n");
  printf("\n");
  printf("   --mask maskfile : mask used with --vote or --sort\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf("Concatenates input data sets.\n");
  printf("EXAMPLES:\n");
  printf("  mri_concat --i f1.mgh --i f2.mgh --o cout.mgh\n");
  printf("  mri_concat f1.mgh f2.mgh --o cout.mgh\n");
  printf("  mri_concat f*.mgh --o cout.mgh\n");
  printf("  mri_concat f*.mgh --o coutmn.mgh --mean\n");
  printf("  mri_concat f*.mgh --o coutdiff.mgh --paired-diff\n");
  printf("  mri_concat f*.mgh --o coutdiff.mgh --paired-diff-norm\n");
  printf("  mri_concat f*.mgh --o coutdiff.mgh --paired-diff-norm1\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void) {
  if (ninputs == 0) {
    printf("ERROR: no inputs specified\n");
    exit(1);
  }
  if (out == NULL) {
    printf("ERROR: no output specified\n");
    exit(1);
  }
  if(DoPairedDiff && DoPairedAvg) {
    printf("ERROR: cannot specify both --paried-diff-xxx and --paried-avg \n");
    exit(1);
  }
  if (DoPairedDiffNorm1 && DoPairedDiffNorm2) {
    printf("ERROR: cannot specify both --paried-diff-norm1 and --paried-diff-norm2 \n");
    exit(1);
  }
  if (DoPairedDiffNorm && DoPairedDiffNorm1) {
    printf("ERROR: cannot specify both --paried-diff-norm and --paried-diff-norm1 \n");
    exit(1);
  }
  if (DoPairedDiffNorm && DoPairedDiffNorm2) {
    printf("ERROR: cannot specify both --paried-diff-norm and --paried-diff-norm2 \n");
    exit(1);
  }
  if(DoMean && DoStd){
    printf("ERROR: cannot --mean and --std\n");
    exit(1);
  }
  if(mask && !DoVote && !DoSort){
    printf("ERROR: --mask only valid with --vote or --sort\n");
    exit(1);
  }
  if(DoStd && DoVar){
    printf("ERROR: cannot compute std and var, you bonehead.\n");
    exit(1);
  }
  if(DoAbs + DoPos + DoNeg > 1){
    printf("ERROR: do not use more than one of --abs, --pos, --neg\n");
    exit(1);
  }


  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  return;
}

MATRIX *GroupedMeanMatrix(int ngroups, int ntotal)
{
  int nper,r,c;
  MATRIX *M;

  nper = ntotal/ngroups;
  if(nper*ngroups != ntotal){
    printf("ERROR: --gmean, Ng must be integer divisor of Ntotal\n");
    return(NULL);
  }

  M = MatrixAlloc(nper,ntotal,MATRIX_REAL);
  for(r=1; r <= nper; r++){
    for(c=r; c <= ntotal; c += nper)
      M->rptr[r][c] = 1.0/ngroups;
  }
  return(M);
}
