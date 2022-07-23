/**
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>
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
#include "mri_identify.h"
#include "cmdargs.h"

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

const char *Progname = NULL;
int debug = 0;
#define NInMAX 400000 // such a large number may break valgrind
char *inlist[NInMAX];
int ninputs = 0;
char flist[NInMAX][2000];
int nthf=0;
char *out = NULL;
MRI *mritmp, *mritmp0, *mriout, *mask=NULL;
char *maskfile = NULL;
int DoMean=0;
int DoNormMean=0;
int DoNorm1=0;
int DoMeanDivN=0;
int DoMedian=0;
int DoSum=0;
int DoVar=0;
int DoStd=0;
int DoMax=0;
int DoMaxIndex=0;
int DoMaxIndexPrune = 0;
int DoMaxIndexAdd = 0 ;
int MaxIndexAdd = 0 ;
int DoMin=0;
int DoConjunction=0;
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
MATRIX *FrameWeight = NULL;
int ngroups = 0;
MATRIX *GroupedMeanMatrix(int ngroups, int ntotal);
char tmpstr[2000];

int DoPCA = 0;
MRI *PCAMask = NULL;
char *PCAMaskFile = NULL;
int DoSCM = 0; // spat cor matrix
int DoCheck = 1;
int DoTAR1 = 0, TAR1DOFAdjust = 1;
int NReplications = 0;

int DoPrune = 0;
MRI *PruneMask = NULL;

int DoRMS = 0; // compute root-mean-square on multi-frame input
int DoCumSum = 0;
int DoFNorm = 0;
char *rusage_file=NULL;
//MRI *MRIzconcat(MRI *mri1, MRI *mri2, int nskip, MRI *out);


/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs, nthin, nframestot=0, nr=0,nc=0,ns=0, fout;
  int r,c,s,f,outf,nframes,err,nthrep,AllZero;
  double v, v1, v2, vavg, vsum;
  int inputDatatype=MRI_UCHAR;
  MATRIX *Upca=NULL,*Spca=NULL;
  MRI *Vpca=NULL;
  char *stem;

  nargs = handleVersionOption(argc, argv, "mri_concat");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0)
  {
    usage_exit();
  }

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if(maskfile)
  {
    printf("Loading mask %s\n",maskfile);
    mask = MRIread(maskfile);
    if(mask == NULL)
    {
      exit(1);
    }
  }

  printf("ninputs = %d\n",ninputs);
  if(DoCheck)
  {
    printf("Checking inputs\n");
    for(nthin = 0; nthin < ninputs; nthin++)
    {
      if(Gdiag_no > 0 || debug)
      {
        printf("Checking %2d %s\n",nthin,inlist[nthin]);
        fflush(stdout);
      }
      mritmp = MRIreadHeader(inlist[nthin],MRI_VOLUME_TYPE_UNKNOWN);
      if (mritmp == NULL)
      {
        printf("ERROR: reading %s\n",inlist[nthin]);
        exit(1);
      }
      if (nthin == 0)
      {
        nc = mritmp->width;
        nr = mritmp->height;
        ns = mritmp->depth;
      }
      if (mritmp->width != nc ||
          mritmp->height != nr ||
          mritmp->depth != ns)
      {
        printf("ERROR: dimension mismatch between %s and %s\n",
               inlist[0],inlist[nthin]);
        exit(1);
      }

      nframestot += mritmp->nframes;
      inputDatatype = mritmp->type; // used by DoKeepDatatype option
      MRIfree(&mritmp);
    }
  }
  else
  {
    printf("NOT Checking inputs, assuming nframestot = ninputs\n");
    nframestot = ninputs;
    mritmp = MRIreadHeader(inlist[0],MRI_VOLUME_TYPE_UNKNOWN);
    if (mritmp == NULL)
    {
      printf("ERROR: reading %s\n",inlist[0]);
      exit(1);
    }
    nc = mritmp->width;
    nr = mritmp->height;
    ns = mritmp->depth;
    MRIfree(&mritmp);
  }

  printf("nframestot = %d\n",nframestot);

  if (DoRMS)
  {
    if (ninputs != 1)
    {
      printf("ERROR: --rms supports only single input w/ multiple frames\n");
      exit (1);
    }
    if (nframestot == 1)
    {
      printf("ERROR: --rms input must have multiple frames\n");
      exit (1);
    }
  }

  if(ngroups != 0)
  {
    printf("Creating grouped mean matrix ngroups=%d, nper=%d\n",
           ngroups,nframestot/ngroups);
    M = GroupedMeanMatrix(ngroups,nframestot);
    if(M==NULL)
    {
      exit(1);
    }
    if(debug)
    {
      MatrixPrint(stdout,M);
    }
  }

  if(M != NULL) {
    if(nframestot != M->cols){
      printf("ERROR: dimension mismatch between inputs (%d) and matrix (%d)\n",nframestot,M->rows);
      exit(1);
    }
  }

  if (DoPaired)
  {
    if (remainder(nframestot,2) != 0)
    {
      printf("ERROR: --paired-xxx specified but there are an "
             "odd number of frames\n");
      exit(1);
    }
  }

  if(FrameWeight != NULL) {
    if(nframestot != FrameWeight->rows){
      printf("ERROR: dimension mismatch between inputs (%d) and FrameWeight (%d)\n",nframestot,FrameWeight->rows);
      exit(1);
    }
  }

  printf("Allocing output\n");
  fflush(stdout);
  int datatype=MRI_FLOAT;
  if (DoKeepDatatype)
  {
    datatype = inputDatatype;
  }
  if (DoRMS)
  {
    // RMS always has single frame output
    mriout = MRIallocSequence(nc,nr,ns,datatype,1);
  }
  else
  {
    mriout = MRIallocSequence(nc,nr,ns,datatype,nframestot);
  }
  if (mriout == NULL)
  {
    exit(1);
  }
  printf("Done allocing\n");

  fout = 0;
  for (nthin = 0; nthin < ninputs; nthin++)
  {
    if (DoRMS) break; // MRIrms reads the input frames
    if(Gdiag_no > 0 || debug)
    {
      printf("Loading %dth input %s\n",
             nthin+1,fio_basename(inlist[nthin],NULL));
      fflush(stdout);
    }
    mritmp = MRIread(inlist[nthin]);
    if(mritmp == NULL)
    {
      printf("ERROR: loading %s\n",inlist[nthin]);
      exit(1);
    }
    if(nthin == 0)
    {
      MRIcopyHeader(mritmp, mriout);
      //mriout->nframes = nframestot;
    }
    if(DoAbs)
    {
      if(Gdiag_no > 0 || debug)
      {
        printf("Removing sign from input\n");
      }
      MRIabs(mritmp,mritmp);
    }
    if(DoPos)
    {
      if(Gdiag_no > 0 || debug)
      {
        printf("Setting input negatives to 0.\n");
      }
      MRIpos(mritmp,mritmp);
    }
    if(DoNeg)
    {
      if(Gdiag_no > 0 || debug)
      {
        printf("Setting input positives to 0.\n");
      }
      MRIneg(mritmp,mritmp);
    }
    for(f=0; f < mritmp->nframes; f++) {
      for(c=0; c < nc; c++)      {
        for(r=0; r < nr; r++)        {
          for(s=0; s < ns; s++)          {
            v = MRIgetVoxVal(mritmp,c,r,s,f);
	    if(FrameWeight != NULL) v *= FrameWeight->rptr[fout+1][1];
            MRIsetVoxVal(mriout,c,r,s,fout,v);
          }
        }
      }
      fout++;
    }
    MRIfree(&mritmp);
  }

  if(DoCombine)
  {
    // Average frames from non-zero voxels
    int nhits;
    mritmp = MRIallocSequence(nc,nr,ns,MRI_FLOAT,1);
    MRIcopyHeader(mriout,mritmp);
    for(c=0; c < nc; c++)
    {
      for(r=0; r < nr; r++)
      {
        for(s=0; s < ns; s++)
        {
          nhits = 0;
          vsum = 0;
          for(f=0; f < mriout->nframes; f++)
          {
            v = MRIgetVoxVal(mriout,c,r,s,f);
            if (v > 0)
            {
              vsum += v;
              nhits ++;
            }
          }
          if(nhits > 0 )
          {
            MRIsetVoxVal(mritmp,c,r,s,0,vsum/nhits);
          }
        } // for s
      }// for r
    } // for c
    MRIfree(&mriout);
    mriout = mritmp;
  } // do combine

  if(DoPrune)
  {
    // This computes the prune mask, applied below
    printf("Computing prune mask \n");
    PruneMask = MRIframeBinarize(mriout,FLT_MIN,NULL);
    printf("Found %d voxels in prune mask\n",MRInMask(PruneMask));
  }

  if(DoNormMean)
  {
    printf("Normalizing by mean across frames\n");
    MRInormalizeFramesMean(mriout);
  }
  if(DoNorm1)
  {
    printf("Normalizing by first across frames\n");
    MRInormalizeFramesFirst(mriout);
  }

  if(DoASL)
  {
    printf("Computing ASL matrix matrix\n");
    M = ASLinterpMatrix(mriout->nframes);
  }

  if(M != NULL) {
    printf("Multiplying by matrix\n");
    mritmp = fMRImatrixMultiply(mriout, M, NULL);
    if(mritmp == NULL) exit(1);
    MRIfree(&mriout);
    mriout = mritmp;
  }


  if(DoPaired)
  {
    printf("Combining pairs\n");
    mritmp = MRIcloneBySpace(mriout,-1,mriout->nframes/2);
    for (c=0; c < nc; c++)
    {
      for (r=0; r < nr; r++)
      {
        for (s=0; s < ns; s++)
        {
          fout = 0;
          for (f=0; f < mriout->nframes; f+=2)
          {
            v1 = MRIgetVoxVal(mriout,c,r,s,f);
            v2 = MRIgetVoxVal(mriout,c,r,s,f+1);
            v = 0;
            if(DoPairedAvg) v = (v1+v2)/2.0;
            if(DoPairedSum) v = (v1+v2);
            if(DoPairedDiff) v = v1-v2;  // difference
            if(DoPairedDiffNorm){
              v = v1-v2; // difference
              vavg = (v1+v2)/2.0;
              if (vavg != 0.0) v = v/vavg;
	      else             v = 0;
            }
            if(DoPairedDiffNorm1)
            {
              v = v1-v2; // difference
              if (v1 != 0.0) v = v/v1;
              else           v = 0;
            }
            if(DoPairedDiffNorm2)
            {
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

  if(DoBonfCor)
  {
    DoAdd = 1;
    AddVal = -log10(mriout->nframes);
  }

  if(DoMean)
  {
    printf("Computing mean across frames\n");
    mritmp = MRIframeMean(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }
  if(DoMedian)
  {
    printf("Computing median across frames\n");
    mritmp = MRIframeMedian(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }
  if(DoMeanDivN)
  {
    printf("Computing mean2 = sum/(nframes^2)\n");
    mritmp = MRIframeSum(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
    MRImultiplyConst(mriout, 1.0/(nframes*nframes), mriout);
  }
  if(DoSum)
  {
    printf("Computing sum across frames\n");
    mritmp = MRIframeSum(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }
  if(DoFNorm)
  {
    printf("Normalizing across frames\n");
    mritmp = MRIframeNorm(mriout,NULL,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }
  if(DoTAR1)
  {
    printf("Computing temoral AR1 %d\n",mriout->nframes-TAR1DOFAdjust);
    mritmp = fMRItemporalAR1(mriout,TAR1DOFAdjust,NULL,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoStd || DoVar)
  {
    printf("Computing std/var across frames\n");
    if(mriout->nframes < 2)
    {
      printf("ERROR: cannot compute std from one frame\n");
      exit(1);
    }
    //mritmp = fMRIvariance(mriout, -1, 1, NULL);
    mritmp = fMRIcovariance(mriout, 0, -1, NULL, NULL);
    if(DoStd)
    {
      MRIsqrt(mritmp, mritmp);
    }
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoMax)
  {
    printf("Computing max across all frames \n");
    mritmp = MRIvolMax(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoMaxIndex){
    printf("Computing max index across all frames \n");
    mritmp = MRIvolMaxIndex(mriout,1,NULL,NULL);
    if(DoMaxIndexPrune){
      // Set to 0 any voxels that are all 0 in each input frame
      // Note: not the same as --prune (which sets to 0 if ANY frame=0)
      printf("Pruning max index\n");
      for (c=0; c < nc; c++){
	for (r=0; r < nr; r++) {
	  for (s=0; s < ns; s++){
	    AllZero = 1;
	    for (f=0; f < mriout->nframes; f++){
	      if(fabs(MRIgetVoxVal(mriout,c,r,s,f))>0){
		AllZero = 0;
		break;
	      }
	    }
	    if(AllZero) MRIsetVoxVal(mritmp,c,r,s,0, 0);
	  }
	}
      }
    }
    MRIfree(&mriout);
    mriout = mritmp;
    if(DoMaxIndexAdd){
      printf("Adding %d to index\n",MaxIndexAdd);
      // This adds a value only to the non-zero voxels
      for (c=0; c < nc; c++){
	for (r=0; r < nr; r++) {
	  for (s=0; s < ns; s++){
	    f = MRIgetVoxVal(mriout,c,r,s,0);
	    if(f == 0) continue;
	    MRIsetVoxVal(mriout,c,r,s,0,f+MaxIndexAdd);
	  }
	}
      }
    }
  }

  if(DoConjunction)
  {
    printf("Computing conjunction across all frames \n");
    mritmp = MRIconjunct(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoMin)
  {
    printf("Computing min across all frames \n");
    mritmp = MRIvolMin(mriout,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoSort)
  {
    printf("Sorting \n");
    mritmp = MRIsort(mriout,mask,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoVote)
  {
    printf("Voting \n");
    mritmp = MRIvote(mriout,mask,NULL);
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoMultiply)
  {
    printf("Multiplying by %lf\n",MultiplyVal);
    MRImultiplyConst(mriout, MultiplyVal, mriout);
  }

  if(DoAdd)
  {
    printf("Adding %lf\n",AddVal);
    MRIaddConst(mriout, AddVal, mriout);
  }

  if(DoCumSum){
    printf("Computing cumulative sum\n");
    fMRIcumSum(mriout, mask, mriout);
  }

  if(DoSCM)
  {
    printf("Computing spatial correlation matrix (%d)\n",mriout->nframes);
    mritmp = fMRIspatialCorMatrix(mriout);
    if(mritmp == NULL)
    {
      exit(1);
    }
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoPCA)
  {
    // Saves only non-zero components
    printf("Computing PCA\n");
    if(PCAMaskFile)
    {
      printf("  PCA Mask %s\n",PCAMaskFile);
      PCAMask = MRIread(PCAMaskFile);
      if(PCAMask == NULL)
      {
        exit(1);
      }
    }
    err=MRIpca(mriout, &Upca, &Spca, &Vpca, PCAMask);
    if(err)
    {
      exit(1);
    }
    stem = IDstemFromName(out);
    sprintf(tmpstr,"%s.u.mtx",stem);
    MatrixWriteTxt(tmpstr, Upca);
    sprintf(tmpstr,"%s.stats.dat",stem);
    WritePCAStats(tmpstr,Spca);
    MRIfree(&mriout);
    mriout = Vpca;
  }

  if(NReplications > 0)
  {
    printf("NReplications %d\n",NReplications);
    mritmp = MRIallocSequence(mriout->width,
                              mriout->height,
                              mriout->depth,
                              mriout->type,
                              mriout->nframes*NReplications);
    if(mritmp == NULL)
    {
      exit(1);
    }
    printf("Done allocing\n");
    MRIcopyHeader(mriout,mritmp);
    for(c=0; c < mriout->width; c++)
    {
      for(r=0; r < mriout->height; r++)
      {
        for(s=0; s < mriout->depth; s++)
        {
          outf = 0;
          for(nthrep = 0; nthrep < NReplications; nthrep++)
          {
            for(f=0; f < mriout->nframes; f++)
            {
              v = MRIgetVoxVal(mriout,c,r,s,f);
              MRIsetVoxVal(mritmp,c,r,s,outf,v);
              outf ++;
            }
          }
        }
      }
    }
    MRIfree(&mriout);
    mriout = mritmp;
  }

  if(DoPrune)
  {
    // Apply prune mask that was computed above
    printf("Applying prune mask \n");
    MRImask(mriout, PruneMask, mriout, 0, 0);
  }

  if(DoRMS)
  {
    printf("Computing RMS across input frames\n");
    mritmp = MRIread(inlist[0]);
    MRIcopyHeader(mritmp, mriout);
    MRIrms(mritmp,mriout);
  }

  printf("Writing to %s\n",out);
  err = MRIwrite(mriout,out);
  if(err) exit(err);

  if(debug) PrintRUsage(RUSAGE_SELF, "mri_ca_label ", stdout);
  if(rusage_file) WriteRUsage(RUSAGE_SELF, "", rusage_file);

  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused, rt;
  char **pargv, *option, listfile[2000] ;
  FILE *fp0;

  if (argc < 1)
  {
    usage_exit();
  }

  nargc   = argc;
  pargv = argv;
  while (nargc > 0)
  {

    option = pargv[0];
    if (debug)
    {
      printf("%d %s\n",nargc,option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))
    {
      print_help() ;
    }
    else if (!strcasecmp(option, "--version"))
    {
      print_version() ;
    }
    else if (!strcasecmp(option, "--debug"))
    {
      debug = 1;
    }
    else if (!strcasecmp(option, "--check"))
    {
      DoCheck = 1;
    }
    else if (!strcasecmp(option, "--no-check"))
    {
      DoCheck = 0;
    }
    else if (!strcasecmp(option, "--mean"))
    {
      DoMean = 1;
    }
    else if (!strcasecmp(option, "--median"))
    {
      DoMedian = 1;
    }
    else if (!strcasecmp(option, "--mean-div-n"))
    {
      DoMeanDivN = 1;
    }
    else if (!strcasecmp(option, "--mean2"))
    {
      DoMeanDivN = 1;
    }
    else if (!strcasecmp(option, "--sum"))
    {
      DoSum = 1;
    }
    else if (!strcasecmp(option, "--std"))
    {
      DoStd = 1;
    }
    else if (!strcasecmp(option, "--var"))
    {
      DoVar = 1;
    }
    else if (!strcasecmp(option, "--abs"))
    {
      DoAbs = 1;
    }
    else if (!strcasecmp(option, "--fnorm"))
    {
      DoFNorm = 1;
    }
    else if (!strcasecmp(option, "--pos"))
    {
      DoPos = 1;
    }
    else if (!strcasecmp(option, "--neg"))
    {
      DoNeg = 1;
    }
    else if (!strcasecmp(option, "--max"))
    {
      DoMax = 1;
    }
    else if (!strcasecmp(option, "--max-index"))
    {
      DoMaxIndex = 1;
    }
    else if (!strcasecmp(option, "--max-index-prune"))
    {
      DoMaxIndex = 1;
      DoMaxIndexPrune = 1;
    }
    else if (!strcasecmp(option, "--max-index-add"))
    {
      if (nargc < 1) argnerr(option,1);
      if(! isdigit(pargv[0][0]) && pargv[0][0] != '-' && pargv[0][0] != '+'){
        printf("ERROR: value passed to the --max-index-add flag must be a number\n");
        printf("       If you want to add two images, use --sum or fscalc\n");
        exit(1);
      }
      sscanf(pargv[0],"%d",&MaxIndexAdd);
      DoMaxIndexAdd = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--min"))
    {
      DoMin = 1;
    }
    else if (!strcasecmp(option, "--conjunct"))
    {
      DoConjunction = 1;
    }
    else if (!strcasecmp(option, "--vote"))
    {
      DoVote = 1;
    }
    else if (!strcasecmp(option, "--sort"))
    {
      DoSort = 1;
    }
    else if (!strcasecmp(option, "--norm-mean"))
    {
      DoNormMean = 1;
    }
    else if (!strcasecmp(option, "--norm1"))
    {
      DoNorm1 = 1;
    }
    else if (!strcasecmp(option, "--prune")) DoPrune = 1;
    else if (!strcasecmp(option, "--no-prune")) DoPrune = 0;

    else if (!strcasecmp(option, "--max-bonfcor"))
    {
      DoMax = 1;
      DoBonfCor = 1;
    }
    else if (!strcasecmp(option, "--asl"))
    {
      DoASL = 1;
    }
    else if (!strcasecmp(option, "--paired-avg"))
    {
      DoPaired = 1;
      DoPairedAvg = 1;
    }
    else if (!strcasecmp(option, "--paired-sum"))
    {
      DoPaired = 1;
      DoPairedSum = 1;
    }
    else if (!strcasecmp(option, "--paired-diff"))
    {
      DoPaired = 1;
      DoPairedDiff = 1;
    }
    else if (!strcasecmp(option, "--paired-diff-norm"))
    {
      DoPairedDiff = 1;
      DoPairedDiffNorm = 1;
      DoPaired = 1;
    }
    else if (!strcasecmp(option, "--paired-diff-norm1"))
    {
      DoPairedDiff = 1;
      DoPairedDiffNorm1 = 1;
      DoPaired = 1;
    }
    else if (!strcasecmp(option, "--paired-diff-norm2"))
    {
      DoPairedDiff = 1;
      DoPairedDiffNorm2 = 1;
      DoPaired = 1;
    }
    else if (!strcasecmp(option, "--combine"))
    {
      DoCombine = 1;
    }
    else if (!strcasecmp(option, "--cumsum")) DoCumSum = 1;
    else if (!strcasecmp(option, "--rms"))
    {
      DoRMS = 1;
    }
    else if (!strcasecmp(option, "--keep-datatype"))
    {
      DoKeepDatatype = 1;
    }
    else if (!strcasecmp(option, "--pca"))
    {
      DoPCA = 1;
    }
    else if (!strcasecmp(option, "--chunk")) setenv("FS_USE_MRI_CHUNK","1",1);
    else if (!strcasecmp(option, "--no-chunk") ) unsetenv("FS_USE_MRI_CHUNK");
      
    else if (!strcasecmp(option, "--scm"))
    {
      DoSCM = 1;
    }
    else if ( !strcmp(option, "--pca-mask") )
    {
      if(nargc < 1)
      {
        argnerr(option,1);
      }
      PCAMaskFile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--i") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      inlist[ninputs] = pargv[0];
      ninputs ++;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--mtx") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      matfile = pargv[0];
      M = MatrixReadTxt(matfile, NULL);
      if(M==NULL)
      {
        printf("ERROR: reading %s\n",matfile);
        exit(1);
      }
      nargsused = 1;
    }
    else if( !strcmp(option, "--w") || !strcmp(option, "--wn") )
    {
      // Frame weight
      if(nargc < 1) argnerr(option,1);
      matfile = pargv[0];
      FrameWeight = MatrixReadTxt(matfile, NULL);
      if(FrameWeight==NULL){
        printf("ERROR: reading %s\n",matfile);
        exit(1);
      }
      if(!strcmp(option, "--wn")){
	printf("Normalizing frame weights to sum to 1\n");
	double sum=0;
	for(int n=1; n <= FrameWeight->rows; n++) sum += FrameWeight->rptr[n][1];
	for(int n=1; n <= FrameWeight->rows; n++) FrameWeight->rptr[n][1] /= sum;
      }
      nargsused = 1;
    }
    else if ( !strcmp(option, "--gmean") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&ngroups);
      nargsused = 1;
    }
    else if ( !strcmp(option, "--rep") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&NReplications);
      nargsused = 1;
    }
    else if ( !strcmp(option, "--o") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      out = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--rusage") )
    {
      if (nargc < 1) argnerr(option,1);
      rusage_file = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--mask") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      maskfile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--mul") ){
      if (nargc < 1)argnerr(option,1);
      if(! isdigit(pargv[0][0]) && pargv[0][0] != '-' && 
	 pargv[0][0] != '+' && pargv[0][0] != '.'){
        printf("ERROR: value passed to the --mul flag must be a number\n");
        printf("       If you want to multiply two images, use fscalc\n");
        exit(1);
      }
      sscanf(pargv[0],"%lf",&MultiplyVal);
      DoMultiply = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--zconcat") ) {
      // --zconcat mri1 mri2 nskip out
      if(nargc < 4) argnerr(option,4);
      MRI *mri1 = MRIread(pargv[0]);
      if(mri1 == NULL) exit(1);
      MRI *mri2 = MRIread(pargv[1]);
      if(mri2 == NULL) exit(1);
      int nskip;
      sscanf(pargv[2],"%d",&nskip);
      MRI *out = MRIzconcat(mri1,mri2,nskip,NULL);
      if(out == NULL) exit(1);
      int err = MRIwrite(out,pargv[3]);
      exit(err);
    }
    else if ( !strcmp(option, "--add") ) {
      if (nargc < 1) argnerr(option,1);
      if(! isdigit(pargv[0][0]) && pargv[0][0] != '-' && pargv[0][0] != '+'){
        printf("ERROR: value passed to the --add flag must be a number\n");
        printf("       If you want to add two images, use --sum or fscalc\n");
        exit(1);
      }
      sscanf(pargv[0],"%lf",&AddVal);
      DoAdd = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--tar1"))
    {
      if(nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&TAR1DOFAdjust);
      DoTAR1 = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--f"))
    {
      if(nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%s",listfile);
      fp0 = fopen(listfile,"r");
      if(fp0==NULL)
      {
        printf("ERROR: opening %s\n",listfile);
        exit(1);
      }
      while(1)
      {
        rt = fscanf(fp0,"%s",flist[nthf]);
        if(rt == EOF)
        {
          fclose(fp0);
          break;
        }
        inlist[ninputs] = flist[nthf];
        nthf++;
        ninputs ++;
      }
      nargsused = 1;
    }
    else
    {
      inlist[ninputs] = option;
      ninputs ++;
      if(ninputs > NInMAX)
      {
        printf("ERROR: ninputs > %d\n",NInMAX);
        exit(1);
      }
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
  printf("   --o out \n");
  printf("   --i invol <--i invol ...> (don't need --i) \n");
  printf("   --f listfile : list file has a text list of files (up to %d)\n",
         NInMAX);
  printf("\n");
  printf("   --paired-sum  : compute paired sum (1+2, 3d+4, etc) \n");
  printf("   --paired-avg  : compute paired avg (1+2, 3d+4, etc) \n");
  printf("   --paired-diff : compute paired diff (1-2, 3-4, etc) \n");
  printf("   --paired-diff-norm : same as paired-diff but scale by TP1,2 average \n");
  printf("   --paired-diff-norm1 : same as paired-diff but scale by TP1 \n");
  printf("   --paired-diff-norm2 : same as paired-diff but scale by TP2 \n");
  printf("   --norm-mean         : normalize frames by mean of all TP\n");
  printf("   --norm1             : normalize frames by TP1 \n");
  printf("   --mtx matrix.asc    : multiply by matrix in ascii file\n");
  printf("   --w frameweight.asc : weight each frame by values in ascii file (one val for each frame)\n");
  printf("   --wn frameweight.asc : same as --w but normalizes the frames to sum to 1\n");
  printf("   --gmean Ng          : create matrix to average Ng groups, Nper=Ntot/Ng\n");
  printf("\n");
  printf("   --combine : average frames from non-zero voxels\n");
  printf("             (useful to combine lh.ribbon.mgz and rh.ribbon.mgz)\n");
  printf("   --keep-datatype : write output in same datatype as input\n");
  printf("                    (default is to write output in Float format)\n");
  printf("\n");
  printf("   --abs  : take abs of input\n");
  printf("   --pos  : set input negatives to 0\n");
  printf("   --neg  : set input postives to 0\n");
  printf("   --mean : compute mean of concatenated volumes\n");
  printf("   --median : compute median of concatenated volumes\n");
  printf("   --mean-div-n : compute mean/nframes (good for var) \n");
  printf("   --sum  : compute sum of concatenated volumes\n");
  printf("   --var  : compute var  of concatenated volumes\n");
  printf("   --std  : compute std  of concatenated volumes\n");
  printf("   --max  : compute max  of concatenated volumes\n");
  printf("   --max-index  : compute index of max of concatenated volumes (1-based)\n");
  printf("   --max-index-prune  : max index setting to 0 any voxel where all frames are 0 (not the same as --prune)\n");
  printf("   --max-index-add val  : add val to non-zero max indices)\n");
  printf("   --min  : compute min of concatenated volumes\n");
  printf("   --rep N : replicate N times (over frame)\n");
  printf("   --fnorm : normalize time series at each voxel (remove mean, divide by SSS)\n");
  printf("   --conjunct  : compute voxel-wise conjunction concatenated volumes\n");
  printf("   --vote : most frequent value at each voxel and fraction of occurances\n");
  printf("   --sort : sort each voxel by ascending frame value\n");
  printf("   --tar1 dofadjust : compute temporal ar1\n");
  printf("   --prune : set vox to 0 unless all frames are non-zero\n");
  printf("   --pca  : output is pca. U is output.u.mtx and S is output.stats.dat\n");
  printf("   --pca-mask mask  : Only use voxels whose mask > 0.5\n");
  printf("   --scm  : compute spatial covariance matrix (can be huge!)\n");
  printf("   --zconcat bot top nskip out : concat in the slice direction skipping nskip slices of the top\n");
  printf("\n");
  printf("   --max-bonfcor  : compute max and bonferroni correct (assumes -log10(p))\n");
  printf("   --mul mulval   : multiply by mulval\n");
  printf("   --add addval   : add addval\n");
  printf("\n");
  printf("   --mask maskfile : mask used with --vote or --sort\n");
  printf("   --rms : root mean square (eg. combine memprage)\n");
  printf("           (square, sum, div-by-nframes, square root)\n");
  printf("   --no-check : do not check inputs (faster)\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
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
  printf("\n");
  printf("Conjunction takes the min of the abs across frames\n");
  printf("at each voxel. The output value at the voxel is the min, \n");
  printf("including the true sign of the min. Eg, if the two frames are:\n");
  printf("   +2.1 and +3.4 --> +2.1\n");
  printf("   -2.1 and -3.4 --> -2.1\n");
  printf("   +2.1 and -3.4 --> +2.1\n");
  printf("   -2.1 and +3.4 --> -2.1\n");
  printf("See: Nichols, Brett,Andersson, Wager, and Poline\n");
  printf("NeuroImage, Volume 25, Issue 3, 15 April 2005, 653-660.\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1)
  {
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  }
  else
  {
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  }
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if (ninputs == 0)
  {
    printf("ERROR: no inputs specified\n");
    exit(1);
  }
  if (out == NULL)
  {
    printf("ERROR: no output specified\n");
    exit(1);
  }
  if(DoPairedDiff && DoPairedAvg)
  {
    printf("ERROR: cannot specify both --paried-diff-xxx and --paried-avg \n");
    exit(1);
  }
  if (DoPairedDiffNorm1 && DoPairedDiffNorm2)
  {
    printf("ERROR: cannot specify both --paried-diff-norm1 and --paried-diff-norm2 \n");
    exit(1);
  }
  if (DoPairedDiffNorm && DoPairedDiffNorm1)
  {
    printf("ERROR: cannot specify both --paried-diff-norm and --paried-diff-norm1 \n");
    exit(1);
  }
  if (DoPairedDiffNorm && DoPairedDiffNorm2)
  {
    printf("ERROR: cannot specify both --paried-diff-norm and --paried-diff-norm2 \n");
    exit(1);
  }
  if(DoMean && DoStd)
  {
    printf("ERROR: cannot --mean and --std\n");
    exit(1);
  }
  if(mask && !DoVote && !DoSort)
  {
    printf("ERROR: --mask only valid with --vote or --sort\n");
    exit(1);
  }
  if(DoStd && DoVar)
  {
    printf("ERROR: cannot compute std and var, you bonehead.\n");
    exit(1);
  }
  if(DoAbs + DoPos + DoNeg > 1)
  {
    printf("ERROR: do not use more than one of --abs, --pos, --neg\n");
    exit(1);
  }


  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  return;
}

MATRIX *GroupedMeanMatrix(int ngroups, int ntotal)
{
  int nper,r,c;
  MATRIX *M;

  nper = ntotal/ngroups;
  if(nper*ngroups != ntotal)
  {
    printf("ERROR: --gmean, Ng must be integer divisor of Ntotal\n");
    return(NULL);
  }

  M = MatrixAlloc(nper,ntotal,MATRIX_REAL);
  for(r=1; r <= nper; r++)
  {
    for(c=r; c <= ntotal; c += nper)
    {
      M->rptr[r][c] = 1.0/ngroups;
    }
  }
  return(M);
}

#if 0
/*!
  \fn MRI *MRIzconcat(MRI *mri1, MRI *mri2, int nskip, MRI *out)
  \brief Concatenates mri2 onto mri1 in the z (slice) direction. The
  first nskip slices are removed from mri2 before concatenation. The
  original app for this was to combine two hires suscept slabs (JonP)
  where the top slice of the bottom slab overlapped with the first
  slice of the top slab. The geometry is such that it agrees with the
  bottom slab (mri1). 
*/
MRI *MRIzconcat(MRI *mri1, MRI *mri2, int nskip, MRI *out)
{
  int nslices = mri1->depth + mri2->depth - nskip;
  if(out == NULL) {
    out = MRIallocSequence(mri1->width, mri1->height, nslices, mri1->type, mri1->nframes);
    MRIcopyHeader(mri1, out);
  }
  if(mri1->width != mri2->width){
    printf("ERROR: MRIzconcat(): mri1 and mri2 mismatch width %d %d\n",mri1->width,mri2->width);
    return (NULL);
  }
  if(mri1->height != mri2->height){
    printf("ERROR: MRIzconcat(): mri1 and mri2 mismatch height %d %d\n",mri1->height,mri2->height);
    return (NULL);
  }
  if(mri1->nframes != out->nframes) {
    printf("ERROR: MRIzconcat(): out nframes mismatch %d %d\n",mri1->nframes, out->nframes);
    return (NULL);
  }
  if(mri1->width != out->width) {
    printf("ERROR: MRIzconcat(): out width mismatch %d %d\n",mri1->width, out->width);
    return (NULL);
  }
  if(mri1->height != out->height) {
    printf("ERROR: MRIzconcat(): out height mismatch %d %d\n",mri1->height, out->height);
    return (NULL);
  }
  if(nslices != out->depth) {
    printf("ERROR: MRIzconcat(): out depth mismatch %d %d\n",nslices, out->depth);
    return (NULL);
  }


  MRIcopyPulseParameters(mri1, out);

  int c;
  for(c=0; c < mri1->width; c++){
    int r,s,f;
    for(r=0; r < mri1->width; r++){
      for(f=0; f < mri1->nframes; f++){
	int sout = 0;
	for(s=0; s < mri1->depth; s++){
	  double v = MRIgetVoxVal(mri1,c,r,s,f);
	  MRIsetVoxVal(out,c,r,sout,f,v);
	  sout++;
	}
	for(s=nskip; s < mri2->depth; s++){
	  double v = MRIgetVoxVal(mri2,c,r,s,f);
	  MRIsetVoxVal(out,c,r,sout,f,v);
	  sout++;
	}
      } //f
    } //r
  } //c

  // Need to fix geometry because we want to simply extend mri1
  MATRIX *M = MRIxfmCRS2XYZ(mri1, 0);
  MRIsetVox2RASFromMatrix(out, M);
  MatrixFree(&M);

  return(out);
}
#endif


