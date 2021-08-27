/**
 * @brief a tool for automatically scheduling events for fMRI experiments
 *
 * optseq2 is a tool for automatically scheduling events for
 * rapid-presentation event-related (RPER) fMRI experiments (the schedule
 * is the order and timing of events).
 */
/*
 * Original Author: Doug Greve
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


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <sys/time.h>

#include "error.h"
#include "diag.h"
#include "proto.h"

#include "matrix.h"
#include "mri.h"
#include "matfile.h"
#include "evschutils.h"
#include "version.h"
#include "numerics.h"

/* Things to do:
   1. Automatically compute Ntp such that Null has as much time
      as the average of the non-null stimuli, or, automatically
      compute the number of stimuli based on the nomial values
      such that, given the Ntp, the Null has as much time at
      the non-null.
   1. Assume HRF
   2. Temporal Filtering
   3. Higher order PolyFit
   4. Compensation VRF for different Nper.
   5. Write function to print an update
   8. Add fourier
   9. Compare to optseq

How much of an effect will TPreScan have?
Does CB1 Pre-opt limit the max AVRF?
Does VRF StdDev optimization limit the max AVRF?
How do the various costs perform under stress?
  Stress = different numbers of conditions
           more conditions
           fewer nulls
How different are the schedules depending upon the model?

Can it be shown that optimizing with an FIR will reduce
the impact of model error when an HRF is assumed?

Can something be done to affect the off-diagonals?

*/

#ifdef X
#undef X
#endif

const char *Progname = NULL;

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
static int  stringmatch(const char *str1, const char *str2);
static int PrintUpdate(FILE *fp, int n);
static int CheckIntMult(float val, float res, float tol);
static MATRIX * ContrastMatrix(float *EVContrast,
                               int nEVs, int nPer, int nNuis, int SumDelays);
static MATRIX * AR1WhitenMatrix(double rho, int N);
int debug = 0;

int   Ntp = -1;
float TR = -1.0;
float TPreScan = 0.0;

int   PSDSpeced = 0;
float PSDMin;
float PSDMax;
float dPSD = -1;
float PSDWindow;
int   nPSDWindow;

int   nEvTypes = 0;
float  EvDuration[500];
int    EvReps[500];
int    EvRepsNom[500];
float  PctVarEvReps = 0.0;
int    VarEvRepsPerCond = 0;
char  *EvLabel[500];
float TStimTot;
float TScanTot;

int  PolyOrder = -1;
char  *outstem = NULL;
char  *mtxstem = NULL;
char  *infilelist[1000];
int   nInFiles=0;

int nSearch   = -1; /* target number */
int nSearched;      /* actual number */
int NoSearch  = 0;
float tSearch   = -1;
float tSearched;
float PctUpdate = 10;
int Update = 1;
long seed = -1;
int nKeep = -1;
int nCB1Opt = 0;
float tNullMin = 0.0;
float tNullMax = -1;

EVSCH **EvSchList;
char *SvAllFile=NULL;
FILE *fpSvAll;

int nTaskAvgs;
int   CostId = EVS_COST_EFF;
const char *CostString;
float CostSum, CostSum2, CostAvg, CostStd, SumCorrect, Sum2Correct;
float EffMax, VRFAvgMax;
float VRFAvgStd_Cost_Ratio;
int nSince; /* niterations since one of the kept schedules  has changed */

char *SumFile = NULL;
char *LogFile = NULL;
float PctDone, PctDoneLast, PctDoneSince;
int UpdateNow;
MATRIX *C=NULL;
float EVContrast[1000];
int ContrastSumDelays = 0, nEVContrast=0;
char *CMtxFile=NULL;
double ar1rho = 0;

int penalize = 0;
double penalpha = 0, penT = 0, pendtmin = 0;

/*-------------------------------------------------------------*/
int main(int argc, char **argv) {
  EVSCH *EvSch;
  MATRIX *Xfir=NULL, *Xpoly=NULL, *X=NULL, *Xt=NULL,
                                     *XtX=NULL, *XtXIdeal=NULL, *W=NULL;
  int m,n, nthhit=0;
  //float eff, cb1err, vrfavg, vrfstd, vrfmin, vrfmax, vrfrange;
  char fname[2000];
  FILE *fpsum, *fplog;
  struct timeval tod;
  long tNow, tStart;
  float ftmp=0, effxtxideal=0;
  int Singular;
  int nargs;

  nargs = handleVersionOption(argc, argv, "optseq2");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  infilelist[0] = NULL;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  dump_options(stdout);

  // Create whitening matrix
  if (ar1rho != 0) {
    W = AR1WhitenMatrix(ar1rho,Ntp);
    //MatlabWrite(W,"W.mat","W");
  }

  /* Create PolyFit Matrix */
  if (PolyOrder >= 0) {
    Xpoly = MatrixAlloc(Ntp,PolyOrder+1,MATRIX_REAL);
    for (n=0;n<Ntp;n++) Xpoly->rptr[n+1][1] = 1.0;
    ftmp = Ntp/2.0;
    if (PolyOrder >= 1)
      for (n=0;n<Ntp;n++) Xpoly->rptr[n+1][2] = (n-ftmp)/ftmp;
    if (PolyOrder >= 2)
      for (n=0;n<Ntp;n++) Xpoly->rptr[n+1][3] = pow((n-ftmp),2.0)/(ftmp*ftmp);
  } else Xpoly = NULL;

  /* To do: Create Temporal Filter Matrix */

  /* Check that PSDMin, PSDMax, and EvDur are integer mults of dPSD */
  if (! CheckIntMult(PSDMin,dPSD,dPSD/100)) {
    printf("ERROR: PSDMin is not an integer multiple of dPSD\n");
    exit(1);
  }
  if (! CheckIntMult(PSDMax,dPSD,dPSD/100)) {
    printf("ERROR: PSDMax is not an integer multiple of dPSD\n");
    exit(1);
  }
  for (n=0; n < nEvTypes; n++) {
    if (! CheckIntMult(EvDuration[n]+tNullMin,dPSD,dPSD/100)) {
      if (tNullMin > 0)
        printf("ERROR: Duration of EventType %d (%g sec) + "
               "tNullMin (%g sec) is %g sec, which is not an integer "
               "multiple of the dPSD (%g sec).\n",n+1,
               EvDuration[n],tNullMin,EvDuration[n]+tNullMin,dPSD);
      else
        printf("ERROR: Duration of EventType %d (%g sec) is not an integer multiple "
               "of the dPSD (%g sec)\n",n+1,EvDuration[n],dPSD);
      exit(1);
    }
  }

  /* Check DOF constraint */
  PSDWindow = PSDMax-PSDMin;
  nPSDWindow = rint(PSDWindow/dPSD);
  nTaskAvgs = nPSDWindow*nEvTypes; /* for FIR only */
  printf("nTaskAvgs = %d\n",nTaskAvgs);
  if (nTaskAvgs >= Ntp) {
    printf("\nERROR: DOF Constraint Violation: number of estimates (%d)\n"
           "  is greater than or equal to the number of time points (%d).\n",
           nTaskAvgs,Ntp);
    printf("  To fix: (1) increase the number of time points to %d, or\n"
           "  (2) decrease the number of event types, or, (3) reduce \n"
           "  the size of the PSD window, or (4) increase the dPSD.\n\n",
           nTaskAvgs+1);
    exit(1);
  }

  /* Create the contrast matrix */
  if(nEVContrast > 0) {
    C = ContrastMatrix(EVContrast, nEvTypes,nPSDWindow,
                       PolyOrder+1,ContrastSumDelays);
  }
  if(C){
    if(C->cols != (nEvTypes*nPSDWindow+PolyOrder+1)){
      printf("ERROR: number of cols in C (%d) does not equal "
	     "the number of cols in X (%d)\n",C->cols,nEvTypes*nPSDWindow);
      exit(1);
    }
    printf("C Contrast Matrix ----------------------\n");
    MatrixPrint(stdout,C);
    printf("================== ----------------------\n");
    if(CMtxFile != NULL)  MatrixWriteTxt(CMtxFile,C);
  }

  /* Alloc the event list and load inputs */
  EvSchList = (EVSCH **) calloc(sizeof(EVSCH*),nKeep);
  if (nInFiles > 0) {
    for (n=0;n<nInFiles;n++) {
      //printf("INFO: reading %s\n",infilelist[n]);
      EvSchList[n] = EVSreadPar(infilelist[n]);
      Xfir = EVSfirMtxAll(EvSchList[n], 0, TR, Ntp, PSDMin, PSDMax, dPSD);
      Singular = EVSdesignMtxStats(Xfir, Xpoly, EvSchList[n],C,W);
      MatrixFree(&Xfir);
      if (Singular) continue;
      EVScb1Error(EvSchList[n]);
      EVScost(EvSchList[n], CostId, &VRFAvgStd_Cost_Ratio);
      PrintUpdate(stdout,n);
      //printf("%g %g %g %g %g %g %g\n",stats[0],EvSchList[n]->cb1err,
      //     stats[1],stats[2], stats[3],stats[4],stats[5]);
    }
    if(! NoSearch) EVSsort(EvSchList,nInFiles);
  }

  if(NoSearch) goto PastSearch;

  /* Check that the scan time is sufficient based on max possible */
  TScanTot = Ntp*TR + TPreScan;
  TStimTot = 0.0;
  for (n=0; n < nEvTypes; n++)
    TStimTot += ( (1.0+PctVarEvReps/100)*EvRepsNom[n] * EvDuration[n]);
  if (TStimTot >= TScanTot) {
    printf("ERROR: Time Constraint Violation: the total amount of  \n"
           " stimulation time (%g sec) equals or exceeds the total amount \n"
           " of scanning time (%g sec).\n",
           TStimTot,TScanTot);
    exit(1);
  }
  for (n=0; n < nEvTypes; n++) EvReps[n] = EvRepsNom[n];

  /* Need to warn if scan time is close to insufficient */

  /* Compute the Ideal XtX for the Nominal Ev Reps */
  XtXIdeal = EVSfirXtXIdeal(nEvTypes, EvRepsNom, EvDuration, TR, Ntp,
                            PSDMin, PSDMax, dPSD);
  if (outstem != NULL) {
    sprintf(fname,"%s.xtxideal.mat",outstem);
    MatlabWrite(XtXIdeal,fname,"xtxideal");
  }
  effxtxideal = 1.0/MatrixTrace(MatrixInverse(XtXIdeal,NULL)); /* need dealloc */

  /* --------- Prep for Search --------------*/
  nthhit = 0;
  nSearched = 0;
  gettimeofday(&tod,NULL);
  tStart = tod.tv_sec;
  nSince = 0;
  CostSum = 0.0;
  CostSum2 = 0.0;
  SumCorrect = 0.0;
  Sum2Correct = 0.0;
  EffMax = 0.0;
  VRFAvgMax = 0.0;
  PctDone = 0.0;
  PctDoneLast = 0.0;
  UpdateNow = 1;
  if (SvAllFile != NULL) {
    fpSvAll = fopen(SvAllFile,"w");
    if (fpSvAll == NULL) {
      printf("ERROR: cannot open %s for writing\n",SvAllFile);
      exit(1);
    }
  }
  printf("INFO: LogFile is %s\n",LogFile);
  fplog = fopen(LogFile,"w");
  if (fplog == NULL) {
    printf("ERROR: cannot open %s\n",LogFile);
    exit(1);
  }
  dump_options(fplog);
  fprintf(fplog,"\n");
  fprintf(fplog,"Pct nSrch MinSrch Cost Eff Cb1Err VRFAvg VRFStd VRFMin VRFMax VRFRange nSince\n");
  fprintf(fplog,"\nBeginUpdateLog\n");

  /* ------------->>>>>>>----- Search -----<<<<<<<<<---------------------*/
  while (1) {

    /* Termination Condition */
    gettimeofday(&tod,NULL);
    tNow = tod.tv_sec;
    tSearched = (tNow-tStart)/3600.0;
    if ( (tSearch > 0)  && (tSearched >= tSearch) )break;
    if ( (nSearch > 0)  && (nSearched >= nSearch) ) break;
    nSearched++;

    /* Randomly select Number of Event Repetitions */
    if (PctVarEvReps > 0.0) {
      if (!VarEvRepsPerCond) ftmp = 1.0+2*(drand48()-0.5)*PctVarEvReps/100;
      for (m=0; m < nEvTypes; m++) {
        if (VarEvRepsPerCond) ftmp = 1.0+2*(drand48()-0.5)*PctVarEvReps/100;
        EvReps[m] = (int)nint(ftmp*EvRepsNom[m]);
      }
      MatrixFree(&XtXIdeal);
      XtXIdeal = EVSfirXtXIdeal(nEvTypes, EvReps, EvDuration,
                                TR, Ntp, PSDMin, PSDMax, dPSD);
    }

    /* Synthesize a Sequence and Schedule */
    EvSch = EVSsynth(nEvTypes, EvReps, EvDuration, dPSD,
                     TR*Ntp, TPreScan, nCB1Opt, tNullMin, tNullMax);
    if (EvSch==NULL) {
      printf("ERROR: syntheszing schedule\n");
      exit(1);
    }
    EvSch->nthsearched = nSearched;
    if (penalize) EVSrefractory(EvSch, penalpha, penT, pendtmin);

    /* Construct the FIR Design Matrix */
    Xfir = EVSfirMtxAll(EvSch, 0, TR, Ntp, PSDMin, PSDMax, dPSD);

    /* Compute XtXIdeal Error */
    Xt = MatrixTranspose(Xfir,Xt);
    XtX = MatrixMultiply(Xt,Xfir,XtX);

    EvSch->idealxtxerr = 0;
    for (m=1; m <= Xfir->cols; m++) {
      for (n=1; n <= Xfir->cols; n++) {
        EvSch->idealxtxerr += fabs(XtX->rptr[m][n]-XtXIdeal->rptr[m][n]);
      }
    }
    MatrixFree(&Xt);
    MatrixFree(&XtX);

    Singular = EVSdesignMtxStats(Xfir, Xpoly, EvSch, C, W);
    MatrixFree(&Xfir);

    if (Singular) continue;

    /* Compute the Cost (to be maximized) */
    EVScost(EvSch, CostId, &VRFAvgStd_Cost_Ratio);
//  CostSum += EvSch->cost;
    { // Kahan summation algorithm for correction of sum error accumulation:
      // http://en.wikipedia.org/wiki/Kahan_summation_algorithm
      float y = EvSch->cost - SumCorrect;
      float t = CostSum + y;
      SumCorrect = (t - CostSum) - y;
      CostSum = t;
    }
//  CostSum2 += (EvSch->cost * EvSch->cost);
    { // Kahan summation algorithm for correction of sum error accumulation:
      float y = (EvSch->cost * EvSch->cost) - Sum2Correct;
      float t = CostSum2 + y;
      Sum2Correct = (t - CostSum2) - y;
      CostSum2 = t;
    }
    if (EffMax < EvSch->eff)       EffMax    = EvSch->eff;
    if (VRFAvgMax < EvSch->vrfavg) VRFAvgMax = EvSch->vrfavg;

    /* Save data on each iteration to a file */
    if (SvAllFile != NULL) {
      fprintf(fpSvAll,"%g  %g  %g  %g  %g  %g  %g %g",
              EvSch->cost,EvSch->eff,EvSch->cb1err,EvSch->vrfavg,
              EvSch->vrfstd,EvSch->vrfmin,EvSch->vrfmax,EvSch->idealxtxerr);
      if (PctVarEvReps > 0.0)
        for (m=0; m < nEvTypes; m++) fprintf(fpSvAll,"%d ",EvReps[m]);
      fprintf(fpSvAll,"\n");
    }

    if (nthhit < nKeep && nInFiles == 0) {
      EvSchList[nthhit] = EvSch;
      if (nthhit == nKeep-1) EVSsort(EvSchList,nKeep);
    } else {
      if (EvSch->cost > EvSchList[nKeep-1]->cost) {
        /* Print update before and after the list changes */
        PrintUpdate(fplog,0);
        PrintUpdate(stdout,0);

        EVSfree(&EvSchList[nKeep-1]);
        EvSchList[nKeep-1] = EvSch;
        EVSsort(EvSchList,nKeep);
        nSince = 0;

        PrintUpdate(fplog,0);
        PrintUpdate(stdout,0);
      } else {
        EVSfree(&EvSch);
        nSince++;
      }
    }

    /* Print an update to the terminal */
    if (nSearch > 0) PctDone = 100*nSearched/nSearch;
    else            PctDone = 100*tSearched/tSearch;
    PctDoneSince = PctDone - PctDoneLast;

    if (Update && (PctDoneSince > PctUpdate || UpdateNow ) ) {
      PrintUpdate(fplog,0);
      PrintUpdate(stdout,0);
      PctDoneLast = PctDone;
      UpdateNow = 0;
    }

    nthhit ++;

  }/*----------- Done Search Loop ----------------------------*/
  /*-----------------------------------------------------------*/

  if (Update) {
    PrintUpdate(fplog,0);
    PrintUpdate(stdout,0);
  }
  fprintf(fplog,"EndUpdateLog\n");
  fprintf(fplog,"\nDone\n");
  fclose(fplog);

  /*---------------- Clean-up after loop ------------------------*/
  if (SvAllFile != NULL) fclose(fpSvAll);

  printf("INFO: searched %d iterations for %f hours\n",
         nSearched,tSearched);
  printf("INFO: %g iterations per second\n",nSearched/(tSearched*3600.0));

  if ( (nSearch-nthhit) > 0) {
    printf("INFO: %d/%d schedules were ill-conditioned \n",
           nSearch-nthhit,nSearch);
  }
  CostAvg = CostSum/nthhit;
  CostStd = sqrt((double)CostSum2/(double)nthhit - 
                 (double)CostAvg*(double)CostAvg);

  /*-------- Check for ill-conditioned schedules ----------------------*/
  if (nthhit == 0) {
    printf("ERROR: all schedules found were ill-conditioned. This \n"
           "probably means that you need more scan time (ie, a \n"
           "greater number of time points) or fewer repetitions.\n");
    exit(1);
  }

  if (nthhit < nKeep) {
    printf("WARNING: optseq could only find %d well-conditioned schedules.\n"
           "Try: (1) increasing the number of search iterations, or \n"
           "(2) increasing the number of time points, or (3) reducing \n"
           "the number of repetitions (or any combination of those).\n"
           "I'll proceed keeping only the well-conditioned schedules.\n",
           nthhit);
    nKeep = nthhit;
  }

PastSearch:

  /* Summarize and save the results */
  fpsum = fopen(SumFile,"w");
  if (fpsum == NULL) {
    printf("ERROR: cannot open %s for writing\n",SumFile);
    printf("Printing summary to stdout\n");
    fpsum = stdout;
  }

  dump_options(fpsum);
  fprintf(fpsum,"\n");
  fprintf(fpsum,"\n");
  if (C != NULL) {
    fprintf(fpsum,"EV Contrast: ");
    for (m=0; m < nEvTypes; m++) fprintf(fpsum,"%6.3f ",EVContrast[m]);
    fprintf(fpsum,"\n");
    fprintf(fpsum,"SumDelays = %d\n",ContrastSumDelays);
    fprintf(fpsum,"Contrast Matrix: ----------------- \n");
    MatrixPrint(fpsum,C);
    fprintf(fpsum,"---------------------------------- \n");
  }


  fprintf(fpsum,"\n");
  fprintf(fpsum,"\n");
  if (!NoSearch) {
    fprintf(fpsum,"Searched %d iterations for %f hours\n",
            nSearched,tSearched);
    fprintf(fpsum,"Rate: %g iterations per second\n",
            nSearched/(tSearched*3600.0));

    if ( (nSearch-nthhit) > 0) {
      fprintf(fpsum,"INFO: %d/%d schedules were ill-conditioned \n",
              nSearch-nthhit,nSearch);
    }
    fprintf(fpsum,"Number of iterations since last substitution %d\n",nSince);
    fprintf(fpsum,"Cost Avg/Std: %.7f %.7f\n",CostAvg,CostStd);
    fprintf(fpsum,"Max Eff Encountered:    %g\n",EffMax);
    fprintf(fpsum,"Max VRFAvg Encountered: %g\n",VRFAvgMax);
  } else {
    fprintf(fpsum,"No search was performed.\n");
    CostAvg=0;
    CostStd=-1000;
  }
  fprintf(fpsum,"Eff of Nominal Ideal XtX: %g\n",effxtxideal);
  fprintf(fpsum,"\n");
  fprintf(fpsum,"Rank     Cost   ZCost NthIter  Eff   CB1Err   VRFAvg VRFStd VRFMin VRFMax VRFRng IXtXErr  ");
  if (PctVarEvReps > 0.0) fprintf(fpsum,"NReps");
  fprintf(fpsum,"\n");

  for (n=0;n<nKeep;n++) {
    fprintf(fpsum,
            "%3d %9.4f   %5.3f %5d %7.3f %7.6f %6.2f %6.2f %6.2f %6.2f %6.2f   %6.2f",
            n+1,EvSchList[n]->cost,(EvSchList[n]->cost-CostAvg)/CostStd,
            EvSchList[n]->nthsearched,
            EvSchList[n]->eff,EvSchList[n]->cb1err,
            EvSchList[n]->vrfavg,EvSchList[n]->vrfstd,
            EvSchList[n]->vrfmin,EvSchList[n]->vrfmax,EvSchList[n]->vrfrange,
            EvSchList[n]->idealxtxerr);
    if (PctVarEvReps > 0.0)
      for (m=0; m < nEvTypes; m++)
        fprintf(fpsum,"%3d ",EvSchList[n]->nEvReps[m]);
    fprintf(fpsum,"\n");

    if(outstem != NULL){
      sprintf(fname,"%s-%03d.par",outstem,n+1);
      EVSwritePar(fname,EvSchList[n],EvLabel,TPreScan,TR*Ntp);
    }

    if (mtxstem != NULL) {
      sprintf(fname,"%s_%03d.mat",mtxstem,n+1);
      Xfir = EVSfirMtxAll(EvSchList[n], 0, TR, Ntp, PSDMin, PSDMax, dPSD);
      X = MatrixHorCat(Xfir,Xpoly,NULL);
      MatlabWrite(X,fname,"X");
      MatrixFree(&Xfir);
      MatrixFree(&X);
    }
  }

  fprintf(fpsum,"\n");
  fprintf(fpsum,"Ideal CB1 PMtx--------------------------\n");
  MatrixPrint(fpsum,EVScb1IdealProbMatrix(EvSchList[0]));
  for (n=0;n<nKeep;n++) {
    fprintf(fpsum,"Schedule %d CB1 PMtx--------------------------\n",n+1);
    MatrixPrint(fpsum,EVScb1ProbMatrix(EvSchList[n]));
  }


  if (fpsum != stdout) fclose(fpsum);

  return(0);
}
/* ------------------------------------------------------------------ */
/*-------------------------------------------------------------*/
/*---------------------------------------------------------*/
/*-------------------------------------------------------------*/
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  char fname[2000], *instem;
  FILE *fp;

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
    else if (!strcasecmp(option, "--update"))    Update = 1;
    else if (!strcasecmp(option, "--noupdate"))  Update = 0;
    else if (!strcasecmp(option, "--nosearch"))  NoSearch = 1;
    else if (!strcasecmp(option, "--sumdelays")) ContrastSumDelays = 1;

    else if (stringmatch(option, "--nsearch")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nSearch);
      nargsused = 1;
    } else if (stringmatch(option, "--tsearch")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&tSearch);
      nargsused = 1;
    } else if (stringmatch(option, "--pctupdate")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&PctUpdate);
      nargsused = 1;
    } else if (stringmatch(option, "--focb")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nCB1Opt);
      nargsused = 1;
    } else if (stringmatch(option, "--nkeep")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nKeep);
      nargsused = 1;
    } else if (stringmatch(option, "--seed")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%ld",&seed);
      nargsused = 1;
    } else if (stringmatch(option, "--ntp")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&Ntp);
      nargsused = 1;
    } else if (stringmatch(option, "--tr")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&TR);
      nargsused = 1;
    } else if (stringmatch(option, "--tprescan")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&TPreScan);
      nargsused = 1;
    } else if (stringmatch(option, "--psdwin")) {
      if (nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%f",&PSDMin);
      sscanf(pargv[1],"%f",&PSDMax);
      nargsused = 2;
      if (nth_is_arg(nargc,pargv,2)) {
        sscanf(pargv[2],"%f",&dPSD);
        nargsused++;
      }
      PSDSpeced = 1;
    } else if (stringmatch(option, "--ev")) {
      if (nargc < 3) argnerr(option,3);
      EvLabel[nEvTypes] = pargv[0];
      sscanf(pargv[1],"%f",&EvDuration[nEvTypes]);
      sscanf(pargv[2],"%d",&EvRepsNom[nEvTypes]);
      nargsused = 3;
      nEvTypes++;
    } else if (stringmatch(option, "--polyfit")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&PolyOrder);
      nargsused = 1;
    } else if (stringmatch(option, "--tnullmax")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&tNullMax);
      nargsused = 1;
    } else if (stringmatch(option, "--tnullmin")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&tNullMin);
      nargsused = 1;
    } else if (stringmatch(option, "--ar1")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&ar1rho);
      if (abs(ar1rho) > 1) {
        printf("ERROR: ar1 = %g, must be between -1 and 1\n",ar1rho);
        exit(1);
      }
      nargsused = 1;
    } else if (stringmatch(option, "--pen")) {
      if (nargc < 3) argnerr(option,3);
      penalize = 1;
      sscanf(pargv[0],"%lf",&penalpha);
      sscanf(pargv[1],"%lf",&penT);
      sscanf(pargv[2],"%lf",&pendtmin);
      if (penalpha < 0 || penalpha > 1) {
        printf("ERROR: alpha = %g, must be between 0 and 1\n",penalpha);
        exit(1);
      }
      nargsused = 3;
    } else if (stringmatch(option, "--repvar")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&PctVarEvReps);
      nargsused = 1;
      if (nth_is_arg(nargc,pargv,1)) {
        VarEvRepsPerCond = 1;
        nargsused++;
      } else VarEvRepsPerCond = 0;
    } 
    else if (stringmatch(option, "--evc")) {
      if (nargc < 1) argnerr(option,1);
      if(C){
	printf("ERROR: cannot --evc and --C\n");
	exit(1);
      }
      nEVContrast = 0;
      while (nth_is_arg(nargc,pargv,nEVContrast)) {
        sscanf(pargv[nEVContrast],"%f",&EVContrast[nEVContrast]);
        nEVContrast++;
      }
      if (nEVContrast==0) {
        printf("ERROR: no elements to contrast vector\n");
        exit(1);
      }
      nargsused = nEVContrast;
    } 
    else if (stringmatch(option, "--C")) {
      if(nargc < 1) argnerr(option,1);
      if(nEVContrast != 0){
	printf("ERROR: cannot --evc and --C\n");
	exit(1);
      }
      C = MatrixReadTxt(pargv[0], NULL);
      if (C == NULL) {
	printf("ERROR: loading C %s\n",pargv[0]);
	exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcmp(option, "--o")) {
      if (nargc < 1) argnerr(option,1);
      outstem = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--cmtx")) {
      if (nargc < 1) argnerr(option,1);
      CMtxFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--mtx")) {
      if (nargc < 1) argnerr(option,1);
      mtxstem = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--sum")) {
      if (nargc < 1) argnerr(option,1);
      SumFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--log")) {
      if (nargc < 1) argnerr(option,1);
      LogFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--sviter")) {
      if (nargc < 1) argnerr(option,1);
      SvAllFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--in")) {
      if (nargc < 1) argnerr(option,1);
      infilelist[nInFiles] = pargv[0];
      nInFiles++;
      nargsused = 1;
    } else if (!strcmp(option, "--i")) {
      if (nargc < 1) argnerr(option,1);
      instem = pargv[0];
      nargsused = 1;

      nInFiles = 0;
      sprintf(fname,"%s-%03d.par",instem,nInFiles+1);
      fp = fopen(fname,"r");
      while (fp != NULL) {
        fclose(fp);
        infilelist[nInFiles] = (char*)calloc(sizeof(char),strlen(fname)+1);
        memmove(infilelist[nInFiles],fname,strlen(fname));
        nInFiles++;
        sprintf(fname,"%s-%03d.par",instem,nInFiles+1);
        fp = fopen(fname,"r");
      }
    } else if (!strcmp(option, "--cost")) {
      if (nargc < 1) argnerr(option,1);
      CostString = pargv[0];
      nargsused = 1;
      CostId = EVScostId(CostString);
      if (CostId == EVS_COST_UNKNOWN) {
        printf("ERROR: Cost %s unrecognized\n",CostString);
        exit(1);
      }
      if (CostId == EVS_COST_VRFAVGSTD) {
        if (nargc < 2 || isflag(pargv[1]) ) {
          printf("ERROR: --cost vrfavgstd requires one parameter\n");
          exit(1);
        }
        sscanf(pargv[1],"%f",&VRFAvgStd_Cost_Ratio);
        nargsused ++;
        printf("VRFAvgStd_Cost_Ratio = %g\n",VRFAvgStd_Cost_Ratio);
      }

    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
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
  printf("USAGE: optseq2 \n") ;
  printf("\n");
  printf("Data Acquistion Parameters\n");
  printf("\n");
  printf("  --ntp Ntp : number of time points\n");
  printf("  --tr TR : temporal resolution of acquisition (in sec)\n");
  printf("  --tprescan t : start events t sec before first acquisition\n");

  printf("\n");
  printf("Event Response and Nuisance Descriptors\n");
  printf("\n");
  printf("  --psdwin psdmin psdmax <dPSD> : PSD window specifications\n");
  printf("  --ev label duration nrepetitions\n");
  printf("  --repvar pct <per-evt>: allow nrepetitions to vary by +/- percent\n");
  printf("  --polyfit order  \n");
  printf("  --tnullmin tnullmin : limit min null duration to tnullmin sec  \n");
  printf("  --tnullmax tnullmax : limit max null duration to tnullmax sec  \n");

  printf("\n");
  printf("Searching and Cost Parameters\n");
  printf("\n");
  printf("  --nsearch n : search over n schedules\n");
  printf("  --tsearch t : search for t hours\n");
  printf("  --focb    n : pre-optimize first order counter-balancing\n");
  printf("  --ar1 rho : optimize assuming whitening with AR1\n");
  printf("  --pen alpha T dtmin: penalize for presentations being too close\n");
  printf("  --evc c1 c2 ... cN : event contrast\n");
  printf("  --C cmtx : load contrast from ascii cmtx\n");
  printf("  --cost name <params>: eff, vrfavg, vrfavgstd, effinv\n");
  printf("\n");
  printf("  --sumdelays : sum delays when forming contrast matrix\n");
  printf("  --seed seedval : initialize random number generator to seedval\n");

  printf("\n");
  printf("Output Options\n");
  printf("\n");
  printf("  --nkeep   n : keep n schedules\n");
  printf("  --o outstem  : save schedules in outstem-RRR.par\n");
  printf("  --mtx mtxstem  : save design matrices in mtxstem_RRR.mat\n");
  printf("  --cmtx cmtxfile  : save contrast matrix in ascii cmtxfile\n");
  printf("  --sum file : save summary in file (outstem.sum)\n");
  printf("  --log file : save log in file (outstem.log)\n");
  printf("  --pctupdate pct : print an update after each pct done\n");
  printf("  --sviter file : save info from each iteration in file  \n");

  printf("\n");
  printf("Input/Initialization Options\n");
  printf("\n");
  printf("  --i instem : initialize with instem-RRR.par  \n");
  printf("  --in input-schedule <--in input-schedule >  \n");
  printf("  --nosearch  : just print output for input files\n");

  printf("\n");
  printf("Help, Documentation, and Bug Reporting\n");
  printf("  --help : print help page\n");
  printf("  --version : print version string \n");

  printf("\n");
  printf("%s\n",getVersion().c_str());
  printf("\n");
  printf("\n");
  printf("Optseq Home Page: \n"
         "   http://surfer.nmr.mgh.harvard.edu/optseq\n");
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n"
         "\n"
         "SUMMARY\n"
         "\n"
         "optseq2 is a tool for automatically scheduling events for\n"
         "rapid-presentation event-related (RPER) fMRI experiments (the schedule\n"
         "is the order and timing of events). Events in RPER are presented\n"
         "closely enough in time that their hemodynamic responses will\n"
         "overlap. This requires that the onset times of the events be jittered\n"
         "in order to remove the overlap from the estimate of the hemodynamic \n"
         "response. RPER is highly resistant to habituation, expectation, and set\n"
         "because the subject does not know when the next stimulus will appear\n"
         "or which stimulus type it will be. RPER is also more efficient than\n"
         "fixed-interval event related (FIER) because more stimuli can be\n"
         "presented within a given scanning interval at the cost of assuming \n"
         "that the overlap in the hemodynamic responses will be linear. In SPM\n"
         "parlance, RPER is referred to as 'stochastic design'.\n"
         "\n"
         "The flexibility of RPER means that there are a huge number of possible\n"
         "schedules, and they are not equal. optseq2 randomly samples the space\n"
         "of possible schedules and returns the 'best' one, where the user can\n"
         "control the definition of 'best'. Cost functions include: average\n"
         "efficiency, average variance reduction factor (VRF), and a weighted\n"
         "combination of average and stddev of the VRF. The user can also\n"
         "specify that the first order counter-balancing of the sequence of \n"
         "event-types be pre-optimized.\n"
         "\n"
         "Visit the Optseq Home Page at: \n"
         "   http://surfer.nmr.mgh.harvard.edu/optseq\n"
         "\n"
         "COMMAND-LINE ARGUMENTS\n"
         "\n"
         "--ntp Ntp\n"
         "\n"
         "Number of time points to be acquired during the scan.  This should be\n"
         "for one 'run' not for the entire session. The Total Scanning Time\n"
         "is the number of time points times the TR plus the prescan period,\n"
         "ie, tScanTot = Ntp*TR+tPreScan. \n"
         "\n"
         "--tr TR\n"
         "\n"
         "Time between functional volumes (in seconds).\n"
         "\n"
         "--tprescan tPreScan\n"
         "\n"
         "Time before the acquisition of the first volume to be processed to\n"
         "begin stimulation.\n"
         "\n"
         "--psdwin PSDMin PSDMax <dPSD>\n"
         "\n"
         "Specifications for the FIR event response window. It will be assumed that \n"
         "the entire response can be captured within this window. PSDMin is the \n"
         "minimum PostStimulus Delay (PSD), PSDMax is the maximum PSD. dPSD \n"
         "is the sampling interval within the window. dPSD is optional; if \n"
         "left unset, it will default to the TR. dPSD controls how finely spaced  \n"
         "the event onsets can be scheduled (ie, the onsets will only appear at  \n"
         "integer multiples of the dPSD). \n"
         "\n"
         "--ev label duration nrepetitions \n"
         "\n"
         "Event Type specification. The label is just a text label (which may be \n"
         "more informative than a numeric id). Duration is the number of seconds \n"
         "that the stimulus will be presented; it should be an integer multiple \n"
         "of the dPSD (see --psdwin). Nrepetitions is the number of times that \n"
         "this event type will be presented during the course of the run. The \n"
         "number of repetitions can be optimized using the --repvar option. Use \n"
         "a different --ev flag for each event type. NOTE: DO NOT INCLUDE THE \n"
         "NULL STIMULUS AS AN EVENT TYPE.  The total stimulation time, tStimTot, \n"
         "equals the product of the duration and the number of repetitions \n"
         "summed over all the events. It should be obvious that the total \n"
         "stimulation time must be less than the total scanning time. \n"
         "\n"
         "--repvar pct <per-evt> \n"
         "\n"
         "Allow the number of repetitions of each event type to randomly vary by \n"
         "+/- pct percent from the number specified with --ev. This allows the \n"
         "user to optimize over the number of repetitions. The total stimulation \n"
         "time is computed from the maximum possible number of repetitions. If \n"
         "only the percentage is given, then the relative number of repetitions \n"
         "of each event type will stay constant. If the string 'per-evt' is  \n"
         "appended, then the number of reps for each event type can change \n"
         "independently to each other. \n"
         "\n"
         "--polyfit order \n"
         "\n"
         "Add polynomial regressors as nuisance variables. Order N includes the\n"
         "Nth order polynomial was well as all lower orders. Max order is currently\n"
         "2. Order 0 is a baseline offset; Order 1 is a linear trend; Order 2\n"
         "is a quadradic trend. Cost functions will not explicitly include the \n"
         "nuisance variables. \n"
         "\n"
         "--tnullmin tNullMin \n"
         "\n"
         "Force the NULL stimulus to be at least tNullMin sec between stimuli.\n"
         "Note that this means that the stimulus duration + tNullMin must be\n"
         "an integer multiple of the dPSD.\n"
         "\n"
         "--tnullmax tNullMax \n"
         "\n"
         "Limit the maximum duration of the NULL stimulus to be tNullMax sec.\n"
         " Note: it may not be possible for a given parameter set to keep the NULL \n"
         "stimulus below a certain amount. In this case, the following error \n"
         "message will be printed out 'ERROR: could not enforce tNullMax'. By\n"
         "default, tNullMax is infinite. \n"
         " \n"
         "--nsearch Nsearch \n"
         " \n"
         "Search over Nsearch iterations. optseq will randomly construct Nsearch \n"
         "schedules, compute the cost of each one, and keep the ones with the \n"
         "highest cost. It is not permitted to specify both Nsearch and Tsearch. \n"
         " \n"
         "--tsearch Tsearch \n"
         " \n"
         "Search for Tsearch hours. optseq will randomly construct as many \n"
         "schedules schedules as it can in the given time, compute the cost of \n"
         "each one, and keep the ones with the highest cost.  It is not \n"
         "permitted to specify both Nsearch and Tsearch. \n"
         " \n"
         "--focb nCB1Opt \n"
         " \n"
         "Pre-optimize the first order counter-balancing (FOCB) of the event \n"
         "sequence. This will cause optseq2 to construct nCB1Opt random \n"
         "sequences and keep the one with the best FOCB properties. This will be \n"
         "done for each iteration. Counter balance optimization is not allowed \n"
         "when there is only one event type. \n"
         " \n"
         "--ar1 rho \n"
         " \n"
         "Optimize while whitening with an AR(1) model with parameter rho. rho must\n"
         "be between -1 and +1.\n"
         " \n"
         "--pen alpha T dtmin\n"
         " \n"
         "Penalize for one presentation starting too soon after the previous \n"
         "presentation. The weight is computed as 1 - alpha*exp(-(dt+dtmin)/T), where \n"
         "dt is the time from the offset of the previous stimulus to the onset\n"
         "of the next stimulus. The basic idea here is that the second stimulus \n"
         "will be reduced in amplitude by the weight factor. alpha and T were fit\n"
         "from data presented in Huettel and McCarthy (NI, 2000) to be alpha=0.8 \n"
         "and T = 2.2 sec. \n"
         " \n"
         "--evc C1 C2 ... CN \n"
         " \n"
         "Optimize based on a contrast of the event types. Ci is the contrast \n"
         "weight for event type i. There must be as many weights as event types. \n"
         "Weights are NOT renormalized such that the sum to 1. \n"
         " \n"
         "--cost costname <params> \n"
         " \n"
         "Specify cost function. Legal values are eff, vrfavg, \n"
         "vrfavgstd, effinv. Default is eff. params as any parameters which accompany \n"
         "the given cost function. eff is the cost function which maximizes \n"
         "efficiency (no parameters). vrfavg is the cost function which \n"
         "maximizes the average Variance Reduction Factor (VRF) (no \n"
         "parameters). vrfavgstd maximizes a weighted combination of the average \n"
         "and stddev VRF; there is one parameter, the weight give to the stddev \n"
         "component. effinv optimizes 1/eff to find the worst schedules\n"
         " \n"
         "--sumdelays \n"
         " \n"
         "Sum the delay regression parameters when computing contrast matrix. \n"
         "The event contrast (--evc) specifies how to weight the events when \n"
         "forming the contrast vector. However, there are multiple coefficients \n"
         "per event type corresponding to the delay in the FIR window. By default, \n"
         "a separate row in the contrast matrix is provided for each delay. To \n"
         "sum across the delays instead, use --sumdelays. The contrast matrix\n"
         "will have only one row in this case. \n"
         " \n"
         "--seed seedval \n"
         " \n"
         "Initialize the random number generator to seedval. If no seedval is \n"
         "specified, then one will be picked based on the time of day. optseq2 \n"
         "uses drand48(). \n"
         " \n"
         "--pctupdate pct \n"
         " \n"
         "Print an update line to stdout and the log file after completing each \n"
         "pct percent of the search.  \n"
         " \n"
         "--nkeep nKeep \n"
         " \n"
         "Save nKeep of the best schedules. Increasing this number does not \n"
         "substantially increase the search time, so it is a good idea to  \n"
         "specify more than you think you will need. \n"
         " \n"
         "--o outstem \n"
         " \n"
         "Save schedules in outstem-RRR.par, where RRR is the 3-digit \n"
         "zero-padded schedule rank number (there will be nKeep of them). \n"
         "The schedules will be saved in the Paradigm File Format (see below).\n"
         " \n"
         "--mtx mtxstem \n"
         " \n"
         "Save the FIR design matrices to mtxstem_RRR.mat in Matlab 4 binary \n"
         "format. \n"
         " \n"
         "--cmtx cmtxfile \n"
         " \n"
         "Save the contrast matrix in Matlab 4 binary format. \n"
         " \n"
         "--sum summaryfile \n"
         " \n"
         "optseq2 will create a file which summarizes the search, including \n"
         "all the input parameters as well as characteristics of each of \n"
         "the schedules kept. By default, the summary file will be outstem.sum, \n"
         "but it can be specified explicitly using this flag. See THE SUMMARY  \n"
         "FILE below. \n"
         " \n"
         "--log logfile \n"
         " \n"
         "During the course of the search, optseq2 will print information about \n"
         "the current search status to stdio and to the log file. By default \n"
         "the log file will be outstem.log. The log file will contain a summary \n"
         "of input arguments as well as a series of status lines. A status line \n"
         "will be printed each time there is a change in the list of nKeep best \n"
         "schedules as well as at prespecified regular intervals. By default, \n"
         "the interval is 10%% of the search time, but this can be changed \n"
         "with --pctupdate. Each status line has 12 columns: (1) percent complete, \n"
         "(2) iteration number, (3) minutes since start, (4) best cost, \n"
         "(5) efficiency, (6) CB1Error, (7) vrfavg, (8) vrfstd, (9) vrfmin,  \n"
         "(10) vrfmax, (11) vrfrange, and (12) number of iterations since  \n"
         "last substitution. \n"
         " \n"
         "--pctupdate pct \n"
         " \n"
         "Print a search status to stdio and the log file at regular intervals \n"
         "corresponding to pct percent of the search time. Default is 10%%. \n"
         " \n"
         "--sviter SvIterFile  \n"
         " \n"
         "Save information summary about all the schedules to SvIterFile in \n"
         "ASCII format. Each line will have 7 columns corresponding to: \n"
         "(1) cost, (2) efficiency, (3) cb1err, (4) vrfavg, (5) vrfstd,  \n"
         "(6) vrfmin, (7) vrfmax. This is mainly for exploring the distribution \n"
         "of the various costs. WARNING: this file can grow to be very large. \n"
         " \n"
         "--i instem \n"
         " \n"
         "Load all input schedules that match instem-RRR.par. These can be used \n"
         "to initialize the search (for example, if you want to continue a \n"
         "previous optimization). It is also possible to only generate a summary \n"
         "and/or design matrices of the given input schedules by include the  \n"
         "--nosearch flag. This can be useful for testing schedules that were  \n"
         "optimized under one cost function against another cost function or \n"
         "for testing independently generated schedules. See also --in. \n"
         " \n"
         "--in input-schedule <--in input-schedule > \n"
         " \n"
         "This does the same thing as --i except that each file is specified \n"
         "separately.  \n"
         " \n"
         "--nosearch \n"
         " \n"
         "Do not search for optimal schedules. This can only be used when \n"
         "reading schedules in using --i or --in. See --i for more information. \n"
         "\n"
         "ALGORITHM OVERVIEW\n"
         "\n"
         "optseq2 randomly searches the space of schedules given the constraints \n"
         "on the command-line and keeps the ones that maximize the given cost \n"
         "function.  Each search iteration begins by creating a random order of \n"
         "events with the appropriate number of repetitions for each event \n"
         "type. First order counter-balancing optimization, if done, is \n"
         "performed here. Next, the timing is generated by inserting random \n"
         "amounts of NULL stimulus so that the total stimulation time plus null \n"
         "time is equal to the total scan time.  Event onset times are \n"
         "constrained to be integer multiples of dPSD. An FIR design matrix is \n"
         "created from this schedule. The FIR peristimulus window begins at \n"
         "PSDMin and ends at PSDMax and is incremented by dPSD. If polynomial \n"
         "regressors are specified, they are appended to the FIR matrix to give \n"
         "the final design matrix, hereafter referred as X. The various costs \n"
         "are computed from X. The forward model is then y = XB+n, which has the \n"
         "solution Bhat = inv(XtX)Xy. A contrast is Ghat = C*Bhat, where C is the\n"
         "the contrast matrix.\n"
         " \n"
         "CONTRAST MATRIX \n"
         " \n"
         "By default, the contrast matrix is the identity over all task-related \n"
         "components. The contrast matrix can be changed by specifying --evc  \n"
         "(and possibly --sumdelays). \n"
         " \n"
         "COST FUNCTIONS \n"
         " \n"
         "First-Order Counter-Balancing (FOCB). The FOCB matrix is the \n"
         "Nevt-by-Nevt matrix of probabilities that one event type follows \n"
         "another, where Nevt is the number of event types (excluding the NULL \n"
         "condition). This is computed only from the sequence of events and is \n"
         "independent of the timing (this is why it is referred to as \n"
         "'pre-optimization'). The ideal FOCB matrix can be computed from the \n"
         "number of repetitions for each event type.  The FOCB cost matrix is \n"
         "computed by subtracting the actual probability from the ideal and then \n"
         "dividing by the ideal. The final cost is computed by averaging the \n"
         "absolute values of all elements in the cost matrix. This cost is \n"
         "minimized during pre-optimization. FOCB optimization can be combined \n"
         "with any other cost function. Note: FOCB requires that there be at \n"
         "least 2 event types. \n"
         " \n"
         "Efficiency (eff). Efficiency is defined as eff = 1/trace(C*inv(Xt*X)*Ct) \n"
         "(note: any nuisance regressors are not included in the computation of \n"
         "the trace but are included in the computation of the inverse). The \n"
         "quantity trace(C*inv(XtX)*Ct) is a measure of the sum square error in Ghat \n"
         "(ie, G-Bhat) relative to the noise inherent in the experiment. Therefore,  \n"
         "maximizing eff is a way of finding a schedule that will result in, on \n"
         "average, the least error in Ghat. \n"
         " \n"
         "Average Variance Reduction Factor (vrfavg). The Variance Reduction Factor \n"
         "(VRF) is the amount by which the variance of an individual estimator (ie,  \n"
         "a component of Ghat) is reduced relative to the noise inherent in the \n"
         "experiment. The VRF for a estimator is the inverse of the corresponding \n"
         "component on the diagonal of C*inv(XtX)*Ct. The average VRF is this value \n"
         "averaged across all estimators.  This will yield similar results as when \n"
         "the efficiency is optimized. \n"
         " \n"
         "Average/StdDev Variance Reduction Factor (vrfavgstd). The cost is defined \n"
         "as cost = vrfavg - W*vrfstd, where vrfstd is the standard deviation of  \n"
         "the VRFs and W is a weighting factor specified as a parameter on the  \n"
         "command-line.  This penalizes schedules that result in large variations \n"
         "in the individual VRFs of the estimators. There is currently a bug in \n"
         "the implementation that causes it to mis-state the cost when the number \n"
         "of repetitions are different for different event types. Also, only use \n"
         "this cost when using a prescan window equal to or greater than the  \n"
         "PSD window (otherwise there will be a tendency not to schedule events \n"
         "near the end of the run). \n"
         " \n"
         "THE SUMMARY FILE \n"
         " \n"
         "The summary file summarizes the conditions under which the search was  \n"
         "performed as well as the properties of each schedule found. It also \n"
         "includes the number of iterations searched and the time it took to \n"
         "search them as well as the average and standard deviation of the cost \n"
         "measured over all schedules. It also includes the maximum efficiency \n"
         "and average VRF over all schedules (these will be the same as the  \n"
         "best schedule if the eff or vrfavg cost functions were chosen). \n"
         " \n"
         "Each schedule is summarized in a table with the following columns: \n"
         "(1) Rank, (2) Cost, (3) ZCost, (4) Iteration Number (NthIter), (5)  \n"
         "Efficiency (Eff), (6) FOBC Error (CB1Err), (7) Average VRF (VRFAvg),  \n"
         "(7) StdDev VRF (VRFStd), \n"
         "(8) Minimum VRF (VRFMin), (9) Maximum VRF (VRFMax), and (10) VRF \n"
         "Range (VRFRng). Many of these measures have been described above. \n"
         "ZCost is the number of standard deviations from the average cost \n"
         "(over all schedules). The Iteration Number is the search iteration that \n"
         "that schedule was found on. The first-order counter-balancing  \n"
         "measures come after this table. First, the ideal FOCB probability \n"
         "matrix is printed followed by the actual matrix for each of the \n"
         "schedules. Note: the printed ideal matrix is based on the nominal number \n"
         "of repetitions. See BUGS.\n"
         " \n"
         "CHOOSING PARAMETERS SETS \n"
         " \n"
         "There are several parameters that must be chosen as a group because \n"
         "they rely and/or relate to each other. These parameters are: (1) the \n"
         "number of time points (Ntp), (2) the TR, (3) the prescan window (tPreScan), \n"
         "(4) the duration of each event type (tEv), and (5) the number of repetitions \n"
         "of each event type (nReps). The most basic relationship requires that \n"
         "the total amount of stimulation time (tStimTot) be less than or equal to \n"
         "the total  amount of scan time (tScanTot), where tStimTot = sum(tEv*nReps) \n"
         "(summed over all conditions), and tScanTot = Ntp*TR+tPreScan, so \n"
         " \n"
         "                sum(tEv*nReps) <= Ntp*TR+tPreScan                   (1) \n"
         " \n"
         "If this constraint is not met, you will receive a 'Time Constraint \n"
         "Violation' Error. The total amount of time dedicated to the Null stimulus \n"
         "(tNullTot) is equal to the difference between the total scan time and \n"
         "the total stimulation time: \n"
         " \n"
         "          tNullTot = Ntp*TR+tPreScan - sum(tEv*nReps)               (2) \n"
         " \n"
         "If the parameters are chosen such that equality results in equation (1), \n"
         "then there will not be any time for the Null stimulus, which is generally \n"
         "not good because inter-stimulus jitter is dependent upon inserting  \n"
         "random amounts of Null between non-Null stimuli.  \n"
         " \n"
         "A rule of thumb is to allocate as much time for the Null as one would for  \n"
         "any other stimulus. This can be done by choosing parameters such that \n"
         " \n"
         "        sum(tEv*nReps)(nEv+1)/nEv = Ntp*TR+tPreScan                 (3) \n"
         " \n"
         "where nEv is the number of event types. The schedule can be optimized \n"
         "around this point by allowing the number of repetitions to vary around \n"
         "this nominal value. \n"
         " \n"
         "There is also a DOF constraint which requires that the number of parameters \n"
         "estimated be less than the number of time points, ie \n"
         " \n"
         "       Nbeta = nPSD*nEv+(PolyOrder+1) < Ntp                         (4)\n"
         " \n"
         "where Nbeta is the number of parameters, nPSD is the number of elements \n"
         "in the post-stimulus time window (ie, (PSDMax-PSDMin)/dPSD), and PolyOrder \n"
         "is the order of the nuisance polynomial specified with -polyfit. If this \n"
         "constraint is not met, you will receive a 'DOF Constraint Violation' Error. \n"
         " \n"
         "CHOOSING THE SEARCH TERMINATION CRITERIA \n"
         " \n"
         "The search is terminated when either the maximum number of iterations \n"
         "has been reached or the maximum search time has been reached. It is  \n"
         "impossible  to determine how many iterations to search over because \n"
         "it is not possible to globally determine what the best schedule is nor \n"
         "would it be possible to determine how long it would take a random search \n"
         "to get there. That said, there are some rules of thumb that can be followed, \n"
         "the most basic being that if a 'large' number of schedules have been  \n"
         "searched and the best cost has not changed 'much' in a 'long time', then \n"
         "you are done. Of course, you still have to define 'large', 'much', and \n"
         "'long time'. The summary file can help with this. In particular, there is a \n"
         "line with the number of iterations since the last substitution (ie,  \n"
         "the number of iterations since one of the best nKeep schedules changed). \n"
         "This can be used to judge how long a 'long time' is. The same information \n"
         "can be extracted from the NthIter column of the summary table. At a minimum \n"
         "let it run for 10000 iterations. \n"
         " \n"
         "PARADIGM FILE FORMAT\n"
         "The schedules will be saved in 'paradigm file' format. This format has four \n"
         "columns: time, numeric event id, event duration, and event label. A numeric \n"
         "id of 0 indicates the Null Stimulus. \n"
         " \n"
         "BUGS \n"
         " \n"
         " Also see the Optseq Home page at http://surfer.nmr.mgh.harvard.edu/optseq\n"
         " \n"
         " The vrfavgstd cost function does not work properly if the number of reps\n"
         " is different for different event types. A prescan window should also be\n"
         " specified\n"
         " \n"
         " The ideal counter-balance matrix reported in the summary file will be\n"
         " for the nominal number of reps when the user has selected to optimize\n"
         " over the number of reps making comparisons between the actual and \n"
         " ideal inappropriate (the FOCB error reported for each will be correct\n"
         " however).\n"
         " \n"
         "BUG REPORTING \n"
         " \n"
         " optseq2 is offered with no guarantees and with no obligation or implication\n"
         " that it will be supported in any way. Having said that, you can send\n"
         " bug reports/questions to: analysis-bugs@nmr.mgh.harvard.edu. You must\n"
         " include the following information: (1) optseq2 version, (2) computer \n"
         " operating system, (3) command-line, (4) description of the problem, and\n"
         " (5) log and/or summary files.\n"
         " \n"
         "AUTHOR \n"
         " \n"
         " optseq2 was written by Douglas N. Greve in the Summber of '02.\n"
         " \n"
         "REFERENCES \n"
         " \n"
         "Dale, A.M., Greve, D.N., and Burock, M.A. (1999) Optimal Stimulus \n"
         "equences for Event-Related fMRI.  5th International Conference on \n"
         "Functional Mapping of the Human Brain. Duesseldorf, Germany. June \n"
         "11-16. \n"
         " \n"
         "Dale, A.M. (1999). Optimal experimental design for event-related \n"
         "fMRI. Human Brain Mapping 8:109-114. \n"
         "\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
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
  struct timeval tod;
  FILE *fp;
  char ctmp[2000];

  if (NoSearch && nInFiles == 0) {
    printf("ERROR: must specify input files with --nosearch\n");
    exit(1);
  }

  if (!NoSearch) {
    if (nSearch < 0 && tSearch < 0) {
      printf("ERROR: must specify search time or number of search iterations\n");
      exit(1);
    }
    if (nSearch > 0 && tSearch > 0) {
      printf("ERROR: cannot specify both search time and number of search iterations\n");
      exit(1);
    }
    if (nEvTypes == 0) {
      printf("ERROR: must specify at least one event type\n");
      exit(1);
    }

    if (nCB1Opt > 0 && nEvTypes == 1) {
      printf("ERROR: cannot counter-balance sequence with one event type\n");
      exit(1);
    }
    if (seed < 0) {  /* Set seed */
      gettimeofday(&tod,NULL);
      seed = tod.tv_usec;
      srand48(seed);
      printf("INFO: Setting srand48() seed to %ld\n",seed);
    } else srand48(seed);
  }

  if (nInFiles > 0 && nKeep > 0) {
    printf("ERROR: cannot spec input file and nkeep\n");
    exit(1);
  }
  if (nInFiles > 0) {
    nKeep = nInFiles;
  }

  if (nKeep < 0) {
    printf("ERROR: must specify number of schedules to keep\n");
    exit(1);
  }
  if (Ntp < 0) {
    printf("ERROR: must specify number of time points (--ntp)\n");
    exit(1);
  }
  if (TR < 0) {
    printf("ERROR: must specify TR (--tr)\n");
    exit(1);
  }
  if (!PSDSpeced) {
    printf("ERROR: must specify the post-stimulus delay window (--psdwin)\n");
    exit(1);
  }
  if (dPSD < 0.0) dPSD = TR;
  if (! CheckIntMult(TR,dPSD,dPSD/100)) {
    printf("ERROR: TR (%g) must be a multiple of dPSD (%g)\n",TR,dPSD);
    exit(1);
  }

  if (outstem == NULL && infilelist[0] == NULL) {
    printf("ERROR: must specify either an input or an output\n");
    exit(1);
  }
  if (outstem != NULL) {
    sprintf(ctmp,"%s-001.par",outstem);
    fp = fopen(ctmp,"w");
    if (fp == NULL) {
      printf("ERROR: cannot open %s for writing\n",ctmp);
      exit(1);
    }
    fclose(fp);

    if (SumFile == NULL) {
      SumFile = (char *) calloc(sizeof(char), strlen(outstem)+100);
      sprintf(SumFile,"%s.sum",outstem);
    }
    if (LogFile == NULL) {
      LogFile = (char *) calloc(sizeof(char), strlen(outstem)+100);
      sprintf(LogFile,"%s.log",outstem);
    }
  }

  if (PolyOrder > 2) {
    printf("ERROR: PolyOrder=%d, Max is 2\n",PolyOrder);
    exit(1);
  }

  CostString = EVScostString(CostId);

  if (SumFile == NULL && outstem == NULL) {
    printf("ERROR: must specify a summary file explicity when not "
           "specifying an output\n");
    exit(1);
  }

  if (LogFile == NULL && outstem == NULL) {
    printf("ERROR: must specify a log file explicity when not "
           "specifying an output\n");
    exit(1);
  }

  /* Check that the summary file is writable */
  fp = fopen(SumFile,"w");
  if (fp == NULL) {
    printf("ERROR: cannot open %s for writing\n",SumFile);
    exit(1);
  }
  fclose(fp);

  if (! strcmp(SumFile,LogFile) ) {
    printf("ERROR: do not set the summary and log files to the same file\n");
    exit(1);
  }

  if (nEVContrast > 0 && nEVContrast != nEvTypes) {
    printf("ERROR: contrast vector length != number of event types\n");
    exit(1);
  }


  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  int n;

  fprintf(fp,"optseq2\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"NoSearch  = %d\n",NoSearch);
  if (nSearch > 0) fprintf(fp,"nSearch  = %d\n",nSearch);
  if (tSearch > 0) fprintf(fp,"tSearch  = %g hours\n",tSearch);
  fprintf(fp,"nKeep    = %d\n",nKeep);
  fprintf(fp,"PctUpdate  = %f\n",PctUpdate);
  fprintf(fp,"nCB1Opt  = %d\n",nCB1Opt);
  fprintf(fp,"seed     = %ld\n",seed);
  fprintf(fp,"Ntp  = %d\n",Ntp);
  fprintf(fp,"TR   = %g\n",TR);
  fprintf(fp,"TPreScan   = %g\n",TPreScan);
  fprintf(fp,"PSD Window   = %g %g %g\n",PSDMin,PSDMax,dPSD);
  fprintf(fp,"nEvTypes = %d\n",nEvTypes);
  fprintf(fp,"EvNo    Label Duration nRepsNom\n");
  for (n=0;n<nEvTypes;n++)
    fprintf(fp,"%2d %10s %6.3f  %3d\n",n+1,EvLabel[n],
            EvDuration[n],EvRepsNom[n]);
  fprintf(fp,"PctVarEvReps = %g\n",PctVarEvReps);
  fprintf(fp,"VarEvRepsPerCond = %d\n",VarEvRepsPerCond);
  fprintf(fp,"PolyOrder = %d\n",PolyOrder);
  fprintf(fp,"tNullMax = %g\n",tNullMax);
  fprintf(fp,"tNullMin = %g\n",tNullMin);
  if (outstem != NULL) printf("outstem = %s\n",outstem);
  if (SvAllFile != NULL) printf("SvAllFile = %s\n",SvAllFile);
  if (nInFiles != 0) {
    fprintf(fp,"nInFiles: %d\n",nInFiles);
    for (n=0; n < nInFiles; n++)
      fprintf(fp,"%2d infile = %s\n",n,infilelist[n]);
  }
  fprintf(fp,"AR1 = %g\n",ar1rho);
  if (!penalize) fprintf(fp,"No refractory penalty\n");
  else fprintf(fp,"Refractory penalty: alpha = %g, T = %g, dtmin = %g\n",
                 penalpha,penT,pendtmin);
  fprintf(fp,"Cost = %s\n",CostString);
  fprintf(fp,"OutStem = %s\n",outstem);
  fprintf(fp,"Summary File = %s\n",SumFile);

  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int isflag(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth) {
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if (nargc <= nth) {
    //printf("nthisarg: not enough %d %d\n",nargc,nth);
    return(0);
  }

  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) {
    //printf("nthisarg: %s is a flag\n",argv[nth]);
    return(0);
  }

  return(1);
}
/*------------------------------------------------------------*/
static int stringmatch(const char *str1, const char *str2) {
  if (! strcmp(str1,str2)) return(1);
  return(0);
}
/*------------------------------------------------------------*/
static int PrintUpdate(FILE *fp, int n) {
  if (EvSchList[n] == NULL) return(1);

  fprintf(fp,
          "%2.0f %7d  %6.1f   %g  %g  %6.4f  %g  %g  %g  %g  %g %g %d\n",
          PctDone, nSearched, tSearched*60,
          EvSchList[n]->cost,EvSchList[n]->eff,
          EvSchList[n]->cb1err,EvSchList[n]->vrfavg,
          EvSchList[n]->vrfstd,EvSchList[n]->vrfmin,
          EvSchList[n]->vrfmax,EvSchList[n]->vrfrange,
          EvSchList[n]->idealxtxerr,
          nSince);
  fflush(fp);
  return(0);
}

/*------------------------------------------------------------
  CheckIntMult() - checks that val is an integer multiple
  of res (to within tolerance tol). Returns 1 if it is an
  integer multiple and 0 otherwise.
  ------------------------------------------------------------*/
static int CheckIntMult(float val, float res, float tol) {
  if ( fabs( val/res - rint(val/res) ) > tol ) return(0);
  return(1);
}
/*------------------------------------------------------------*/
static MATRIX * ContrastMatrix(float *EVContrast,
                               int nEVs, int nPer, int nNuis, int SumDelays) {
  int n, d, m, nbeta;

  nbeta = nPer*nEVs + nNuis;

  if (SumDelays) {
    C = MatrixZero(1,nbeta,NULL);
    m = 0;
    for (n=0; n<nEVs; n++) {
      for (d=0; d<nPer; d++) {
        C->rptr[1][m+1] = EVContrast[n];
        m = m + 1;
      }
    }
    return(C);
  }


  C = MatrixZero(nPer,nbeta,NULL);

  for (n=0; n<nEVs; n++) {
    for (d=0; d<nPer; d++) {
      m = d + n*nPer;
      C->rptr[d+1][m+1] = EVContrast[n];
    }
  }
  return(C);
}
/*------------------------------------------------------------*/
static MATRIX * AR1WhitenMatrix(double rho, int N) {
  int m,n,d;
  double v;
  MATRIX *W,*M;
  MATRIX *G=NULL;

  // First create AR1 covariance matrix M
  M = MatrixAlloc(N,N,MATRIX_REAL);
  for (m=0;m<N;m++) {
    for (n=0;n<N;n++) {
      d = abs(m-n);
      v = pow(rho,(double)d);
      M->rptr[m+1][n+1] = v;
      //gsl_matrix_set(M,m,n,v);
    }
  }

  // Invert
  M = MatrixInverse(M,M);

  // Copy to matrix
  G = MatrixCopy(M,G);

  // Perform cholesky decomposition
  // M = G*G'; Note: in matlab: M = G'*G;
  sc_linalg_cholesky_decomp(G);

  // Keep the upper triangular part
  W = MatrixAlloc(N,N,MATRIX_REAL);
  for (m=1;m<=N;m++) {
    for (n=1;n<=N;n++) {
      if (m<=n) W->rptr[m][n] = G->rptr[m][n];
    }
  }

  MatrixFree(&G);
  MatrixFree(&M);

  return(W);
}
