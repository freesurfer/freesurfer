/**
 * @brief compute/optimize cost function of segmentation-based registration
 *
 */
/*
 * Original Author: Greg Grev
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


/*
  BEGINUSAGE --------------------------------------------------------------

  mri_segreg

  --init-reg regfile
  --regheader subject
  --mov fvol
  --out-reg outreg : reg at lowest cost

  --cost costfile
  --sum sumfile : default is outreg.sum
  --cur-reg curregfile : reg at current optimum

  --o out : save final output

  --label labelfile : limit calculation to specified label

  --brute_trans min max delta : brute force translation in all directions
  --brute min max delta : brute force in all directions
  --1dpreopt min max delta : brute force in PE direction

  --fwhm fwhm : smooth input by fwhm mm
  --abs       : compute abs of mov
  --subsamp nsub : only sample every nsub vertices

  --preopt-file file : save preopt results in file
  --preopt-dim dim : 0-5 (def 2) (0=TrLR,1=TrSI,2=TrAP,3=RotLR,4=RotSI,5=RotAP)
  --preopt-only : only preopt, so not optimize

  --slope slope    : set cost slope
  --offset offset  : cost offset (pct)
  --T1, --t1         : assume T1 gray/white contrast
  --T2, --t2, --bold : assume T2/BOLD gray/white contrast (default)
  --gm-gt-wm slope : set cost slope, spec gray matter brighter than WM
  --wm-gt-gm slope : set cost slope, spec WM brighter than gray matter
  --penalty-abs : remove sign from contrast
  --cf cfile  : save cost function values (pct,cost)

  --mask : mask out expected B0 regions
  --no-mask : do not mask out (default)
  --no-cortex-label : do not use cortex label

  --lh-only : only use left hemisphere
  --rh-only : only use right hemisphere

  --gm-proj-frac frac : fraction of cortical thickness (default, 0.5)
  --gm-proj-abs  dist : absolute distance into cortex

  --wm-proj-frac frac : fraction of cortical thickness
  --wm-proj-abs  dist : absolute distance into WM (default 2mm)

  --projabs : project -1mm and + 2mm
  --proj-frac frac : projection fraction

  --frame nthframe : use given frame in input (default = 0)
  --mid-frame : use use middle frame

  --trans Tx Ty Tz : translate input reg (mm)
  --rot   Ax Ay Zz : rotate input reg (degrees)

  --trans-rand Tmax : uniformly dist -Tmax to +Tmax (Tx,Ty,Tz)
  --rot-rand   Amax : uniformly dist -Amax to +Amax (Ax,Ay,Az)

  --interp interptype : interpolation trilinear or nearest (def is trilin)
  --no-crop: do not crop anat (crops by default)
  --profile : print out info about exec time

  --noise stddev : add noise with stddev to input for testing sensitivity
  --seed randseed : for use with --noise

  --nmax nmax   : max number of powell iterations (def 36)
  --tol   tol   : powell inter-iteration tolerance on cost. 
     This is the fraction of the cost that the difference in 
     successive costs must drop below to stop the optimization. 

  --tol1d tol1d : tolerance on powell 1d minimizations

  --1dmin : use brute force 1D minimizations instead of powell
  --n1dmin n1dmin : number of 1d minimization (default = 3)

  --mincost MinCostFile
  --param   ParamFile
  --rms     RMSDiffFile : saves Tx Ty Tz Ax Ay Az RMSDiff MinCost 
              WMMean CtxMean PctContrast C0 Slope NSubSamp UseMask
  --surf surfname : use ?h.surfname instead of lh.white
  --surf-cost basename : saves as basename.?h.mgh
  --init-surf-cost basename0 : saves init cost as basename0.?h.mgh
  --surf-cost-diff diffbase : saves final-init cost as diffbase.?h.mgh
  --surf-con basename : saves final contrast as basename.?h.mgh
  --cost-eval file.dat : save costs at each iteration

  ENDUSAGE ---------------------------------------------------------------
*/

/*
  BEGINHELP --------------------------------------------------------------

  FORMATS

  Data file format can be specified implicitly (through the path name)
  or explicitly. All formats accepted by mri_convert can be used.

  BUGS

  sinc interpolation is broken except for maybe COR to COR.


  BUG REPORTING

  Report bugs to analysis-bugs@nmr.mgh.harvard.edu. Include the following
  formatted as a list as follows: (1) command-line, (2) directory where
  the program was run (for those in the MGH-NMR Center), (3) version,
  (4) text output, (5) description of the problem.

  SEE ALSO

  mri_vol2vol mri_convert, tkregister2


  ENDHELP --------------------------------------------------------------

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <sys/utsname.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"

#include "matrix.h"
#include "mri.h"
#include "version.h"
#include "mri2.h"
#include "mri_identify.h"
#include "MRIio_old.h"
#include "registerio.h"
#include "resample.h"
#include "gca.h"
#include "gcamorph.h"
#include "fio.h"
#include "cmdargs.h"
#include "pdf.h"
#include "timer.h"
#include "fmriutils.h"
#include "numerics.h"
#include "annotation.h"
#include "transform.h"
#include "label.h"

#ifdef X
#undef X
#endif

int MRISbbrSurfs(char *subject);

double *GetSurfCosts(MRI *mov, MRI *notused, MATRIX *R0, MATRIX *R,
		     double *p, int dof, double *costs);
int MinPowell(MRI *mov, MRI *notused, MATRIX *R, double *params,
	      int dof, double ftol, double linmintol, int nmaxiters,
	      char *costfile, double *costs, int *niters);
float compute_powell_cost(float *p) ;
double RelativeSurfCost(MRI *mov, MATRIX *R0);

char *costfile_powell = NULL;

// For some reason, this does not seemed to be defined in math.h
double round(double x);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
#include "tags.h"
static int istringnmatch(const char *str1, const char *str2, int n);
double VertexCost(double vctx, double vwm, double slope, 
		  double center, double sign, double *pct);


int main(int argc, char *argv[]) ;

const char *Progname = NULL;

int debug = 0, gdiagno = -1;

char *movvolfile=NULL;
char *regfile=NULL;
char *outregfile=NULL;
char *sumfile=NULL;
char *curregfile=NULL;

const char *interpmethod = "trilinear";
int   interpcode = SAMPLE_TRILINEAR;
int   sinchw;

MRI *mov, *out;

MATRIX *R0;

char *SUBJECTS_DIR=NULL;
char *subject = NULL;

float ipr, bpr, intensity=1.0;
int float2int,err, nargs;

char tmpstr[2000];

char *SegRegCostFile = NULL;
char  *fspec;
MRI *anat;
MRI *noise=NULL;
MRI *mritmp;
char *outfile=NULL;

int frame = 0;

int SynthSeed = -1;
int AddNoise = 0;
double NoiseStd;

int DoProfile = 0;
int n1dmin = 3;

int nMaxItersPowell = 36;
double TolPowell = 1e-8;
double LinMinTolPowell = 1e-8;

#define NMAX 100
int ntx=0, nty=0, ntz=0, nax=0, nay=0, naz=0;
double txlist[NMAX],tylist[NMAX],tzlist[NMAX];
double axlist[NMAX],aylist[NMAX],azlist[NMAX];

MRIS *lhwm, *rhwm, *lhctx, *rhctx;
MRI *lhsegmask=NULL, *rhsegmask=NULL;
MRI *lhlabel, *rhlabel;
int UseMask = 0;
char *cmdline2;
struct utsname uts;

static int UseLabel = 0 ;
static LABEL *mask_label = NULL ;

int PenaltySign  = -1;
double PenaltySlope = .5;
double PenaltyCenter = 0;
int DoMidFrame = 0;

int Do1DPreOpt = 0;
double PreOptMin, PreOptMax, PreOptDelta, PreOpt, PreOptAtMin;
double PreOptMinTrans, PreOptMaxTrans, PreOptDeltaTrans=-1 ;
int PreOptDim = 2;
char *PreOptFile = NULL;
int PreOptOnly = 0;

char *surfcostbase=NULL, *lhcostfile=NULL, *rhcostfile=NULL;
char *surfconbase=NULL, *lhconfile=NULL, *rhconfile=NULL;
char *surfcost0base=NULL, *lhcost0file=NULL, *rhcost0file=NULL;
char *surfcostdiffbase=NULL;

char *TargConLHFile=NULL,*TargConRHFile=NULL;
MRI *TargConLH=NULL,*TargConRH=NULL;

int UseLH = 1;
int UseRH = 1;

MATRIX *MrotPre=NULL,*MtransPre=NULL,*MscalePre=NULL,*MshearPre=NULL;
double TransRandMax = 0;
double RotRandMax = 0;
char *MinCostFile=NULL;
char *InitCostFile=NULL;

int DoRMSDiff = 1; // this is fast, so ok to do automatically
char *RMSDiffFile = NULL;

int DoSmooth = 0;
double fwhm = 0, gstd = 0;
int nsubsamp = 1;
int nsubsampbrute = 100;

int  DoGMProjFrac = 1;   // default
double GMProjFrac = 0.5; // default
int  DoWMProjFrac = 0;
double WMProjFrac = 0.5;

int  DoGMProjAbs = 0;
double GMProjAbs = +2.0;
int  DoWMProjAbs = 1;    // default
double WMProjAbs =  2.0; // default

int BruteForce = 0;
int   regheader=0;

int DoAbs = 0;
MRI *lhcost=NULL, *lhcost0=NULL, *rhcost=NULL, *rhcost0=NULL;
MRI *lhcostdiff=NULL, *rhcostdiff=NULL;
MRI *lhcon=NULL, *rhcon=NULL;

MRI *lhCortexLabel=NULL, *rhCortexLabel=NULL;
int UseCortexLabel = 1;
int nCostEvaluations=0;
MRI *vsm=NULL;
char *vsmfile = NULL;
double angles[3],xyztrans[3],scale[3],shear[3];
const char *surfname = "white";
int dof = 6; 
char *RelCostFile = NULL;
char *ParamFile = NULL;
int InitSurfCostOnly=0;

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  double costs[8], mincost, p[12], pmin[6];
  double tx, ty, tz, ax, ay, az;
  int nth, n, vno;
  MATRIX *R=NULL, *R00=NULL, *Rdiff=NULL;
  Timer mytimer;
  double secCostTime;
  FILE *fp, *fpMinCost, *fpInitCost, *fpRMSDiff, *fpPreOpt=NULL, *fpRelCost, *fpParam;
  double rmsDiffSum, rmsDiffSum2, rmsDiffMean=0, rmsDiffMax=0, d, rmsDiffStd=0;
  double rcost0, rcost;
  VERTEX *v;
  int nsubsampsave = 0;
  LABEL *label;
  int PrintT1Warning = 0, PrintT2Warning = 0;

  memset(pmin,0,sizeof(pmin));

  std::string cmdline = getAllInfo(argc, argv, "mri_segreg");

  nargs = handleVersionOption(argc, argv, "mri_segreg");
  if(nargs && argc - nargs == 1) exit (0);

  cmdline2 = argv2cmdline(argc,argv);
  uname(&uts);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  // do this before parse
  if(SynthSeed < 0) SynthSeed = PDFtodSeed(); 
  srand48(SynthSeed);

  parse_commandline(argc, argv);
  if(gdiagno > -1) Gdiag_no = gdiagno;
  check_options();
  dump_options(stdout);

  fp = fopen(sumfile,"w");
  dump_options(fp);
  fflush(fp);

  // voxel shift map
  if(vsmfile){
    printf("Loading vsm\n");
    vsm = MRIread(vsmfile);
    if(vsm == NULL) exit(1);
  }

  printf("Loading mov\n");
  mov = MRIread(movvolfile);
  if (mov == NULL) exit(1);

  if(DoAbs){
    printf("Computing abs\n");
    MRIabs(mov,mov);
  }

  // Load an anatomical for reference
  sprintf(tmpstr,"%s/%s/mri/orig.mgz",SUBJECTS_DIR,subject);
  anat = MRIread(tmpstr); // Just need a template
  if(anat == NULL) exit(1);
  mritmp = MRIchangeType(anat,MRI_FLOAT,0,0,0);
  MRIfree(&anat);
  anat = mritmp;

  if(regheader) {
    printf("Computing registration based on scanner-to-scanner\n");
    R0 = MRItkRegMtx(anat,mov,NULL);
  }

  R00 = MatrixCopy(R0,NULL);

  if(MrotPre || MtransPre || MscalePre || MshearPre){
    printf("Applying Pre Transform to input reg\n");
    if(MrotPre){
      printf("Rot:\n");
      MatrixPrint(stdout,MrotPre);
      fprintf(fp,"Rot:\n");
      MatrixPrint(fp,MrotPre);
      R0 = MatrixMultiply(MrotPre,R0,R0);
    }
    if(MtransPre){
      printf("Trans Mtx:\n");
      MatrixPrint(stdout,MtransPre);
      fprintf(fp,"After Trans:\n");
      MatrixPrint(fp,MtransPre);
      R0 = MatrixMultiply(MtransPre,R0,R0);
    }
    if(MscalePre){
      printf("Scale Mtx:\n");
      MatrixPrint(stdout,MscalePre);
      fprintf(fp,"After Scale:\n");
      MatrixPrint(fp,MscalePre);
      R0 = MatrixMultiply(MscalePre,R0,R0);
    }
    if(MshearPre){
      printf("Shear Mtx:\n");
      MatrixPrint(stdout,MshearPre);
      fprintf(fp,"After Shear:\n");
      MatrixPrint(fp,MshearPre);
      R0 = MatrixMultiply(MshearPre,R0,R0);
    }
    printf("New input reg:\n");
    MatrixPrint(stdout,R0);
  }

  if(DoMidFrame) frame = nint(mov->nframes/2);

  if(mov->nframes > 1){
    printf("Extracting frame %d\n",frame);
    mritmp = fMRIframe(mov, frame, NULL);
    MRIfree(&mov);
    mov = mritmp;
  }

  if(DoSmooth){
    printf("Smoothing input by fwhm=%lf, gstd=%lf\n",fwhm,gstd);
    MRImaskedGaussianSmooth(mov, NULL, gstd, mov);
  }

  // Loads surfaces, creates masks, does not actually reg.
  MRISbbrSurfs(subject);

  if(UseLH && UseCortexLabel){
    // Load the LH cortex label
    label = LabelRead(subject, "lh.cortex.label");
    if(label == NULL){
      printf("Warning cannot find lh.cortex.label ... continuing\n");
      lhCortexLabel = NULL;
    }
    else{
      printf("Using lh.cortex.label\n");
      lhCortexLabel = MRISlabel2Mask(lhwm, label, NULL);
      if(lhCortexLabel == NULL) exit(1);
      LabelFree(&label);
    }
  }
  if(UseRH && UseCortexLabel){
    // Load the RH cortex label
    label = LabelRead(subject, "rh.cortex.label");
    if(label == NULL){
      printf("Warning cannot find rh.cortex.label ... continuing\n");
      rhCortexLabel = NULL;
    }
    else{
      printf("Using rh.cortex.label\n");
      rhCortexLabel = MRISlabel2Mask(rhwm, label, NULL);
      if(rhCortexLabel == NULL) exit(1);
      LabelFree(&label);
    }
  }



  if(AddNoise){
    // Seed the random number generator just in case
    printf("Adding noise, Seed = %d, stddev = %lf\n",SynthSeed,NoiseStd);
    fprintf(fp,"Adding noise, Seed = %d, stddev = %lf\n",SynthSeed,NoiseStd);
    noise = MRIrandn(mov->width,mov->height,mov->depth,mov->nframes,
                     0,NoiseStd,NULL);
    mov = MRIadd(mov, noise, mov);
  }

  R = MatrixCopy(R0,NULL);

  // Init parameters
  for(nth=0; nth < 12; nth++) p[nth] = 0.0;
  p[6] = 1.0; p[7] = 1.0; p[8] = 1.0;

  // Compute cost at initial
  if(surfcost0base){
    if(UseLH){
      sprintf(tmpstr,"%s.lh.mgh",surfcost0base);
      lhcost0file = strcpyalloc(tmpstr);
    }
    if(UseRH){
      sprintf(tmpstr,"%s.rh.mgh",surfcost0base);
      rhcost0file = strcpyalloc(tmpstr);
    }
    nsubsampsave=nsubsamp;
    nsubsamp=1;
    GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
    nsubsamp=nsubsampsave;
    if(UseLH) lhcost0 = MRIcopy(lhcost,NULL);
    if(UseRH) rhcost0 = MRIcopy(rhcost,NULL);
    if(InitSurfCostOnly){
      if(lhcost0file)  MRIwrite(lhcost0,lhcost0file);
      if(rhcost0file)  MRIwrite(rhcost0,rhcost0file);
      printf("Init surf cost only requested, so exiting now\n");
      printf("mri_segreg done\n");
      exit(0);
    }
    free(lhcost0file);lhcost0file=NULL;
    free(rhcost0file);rhcost0file=NULL;
    /*If the lhcost0file is not set to NULL, then it adds overhead to
     GetSurfCosts(). Not sure why I don't just save the cost0 file
     and don't mess with it again.*/
  }

  if(DoProfile){
    printf("Profiling over 100 iterations, nsubsamp = %d\n",nsubsamp);
    mytimer.reset() ;
    for(n=0; n < 100; n++) GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
    secCostTime = mytimer.seconds();
    printf("ttot = %g\n",secCostTime);
    printf("tper = %g\n",secCostTime/100);
    exit(0);
  }

  // Compute relative initial cost 
  rcost0 = RelativeSurfCost(mov, R0);
  if(RelCostFile){
    fpRelCost = fopen(RelCostFile,"w");
    fprintf(fpRelCost,"%lf\n",rcost0);
    exit(0);
  }

  GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
  fprintf(fp,"Initial costs ----------------\n");
  fprintf(fp,"Number of surface hits %d\n",(int)costs[0]);  
  fprintf(fp,"WM  Intensity0 %10.4lf +/- %8.4lf\n",costs[1],costs[2]); 
  fprintf(fp,"Ctx Intensity0 %10.4lf +/- %8.4lf\n",costs[4],costs[5]); 
  fprintf(fp,"Pct Contrast0  %10.4lf +/- %8.4lf\n",costs[6],costs[3]); 
  fprintf(fp,"Cost %8.4lf\n",costs[7]); 
  fprintf(fp,"RelCost %8.4lf\n",rcost0);
  fflush(fp);

  printf("Initial costs ----------------\n");
  printf("Number of surface hits %d\n",(int)costs[0]);  
  printf("WM  Intensity %10.4lf +/- %8.4lf\n",costs[1],costs[2]); 
  printf("Ctx Intensity %10.4lf +/- %8.4lf\n",costs[4],costs[5]); 
  printf("Pct Contrast  %10.4lf +/- %8.4lf\n",costs[6],costs[3]); 
  printf("Cost %8.4lf\n",costs[7]); 
  printf("RelCost %8.4lf\n",rcost0);
  if(InitCostFile){
    fpInitCost = fopen(InitCostFile,"w");
    //MinCost, WMMean, CtxMean, PctContrast
    fprintf(fpInitCost,"%lf %lf %lf %lf \n",costs[7],costs[1],costs[4],costs[6]);
    fclose(fpInitCost);
  }

  if(costs[6] < 0 && PenaltySign < 0) PrintT1Warning = 1;
  if(PrintT1Warning){
    fprintf(fp,"\n\n");
    fprintf(fp,"WARNING: initial G-W contrast is negative, but expecting positive.\n");
    fprintf(fp,"If the mov data has a T1 contrast, re-run with --T1\n");
    fprintf(fp,"\n\n");
    printf("\n\n");
    printf("WARNING: initial G-W contrast is negative, but expecting positive.\n");
    printf("If the mov data has a T1 contrast, re-run with --T1\n");
    printf("\n\n");
  }

  if(costs[6] > 0 && PenaltySign > 0) PrintT2Warning = 1;
  if(PrintT2Warning){
    fprintf(fp,"\n\n");
    fprintf(fp,"WARNING: initial G-W contrast is positive, but expecting negative.\n");
    fprintf(fp,"If the mov data has a T2 contrast, re-run with --T2\n");
    fprintf(fp,"\n\n");
    printf("\n\n");
    printf("WARNING: initial G-W contrast is positive, but expecting negative.\n");
    printf("If the mov data has a T2 contrast, re-run with --T2\n");
    printf("\n\n");
  }

  if(Do1DPreOpt){
    printf("\n");
    printf("Performing 1D preopt %g %g %g\n",PreOptMin,PreOptMax,PreOptDelta);
    fprintf(fp,"\n");
    fprintf(fp,"Performing 1D preopt %g %g %g\n",PreOptMin,PreOptMax,PreOptDelta);
    mincost = 10e10;
    PreOptAtMin = 0;
    nth = 0;
    if(PreOptFile) fpPreOpt = fopen(PreOptFile,"w");
    for(PreOpt = PreOptMin; PreOpt <= PreOptMax; PreOpt += PreOptDelta){
      p[PreOptDim] = PreOpt;
      GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
      if(costs[7] < mincost) {
	mincost = costs[7];
	PreOptAtMin = PreOpt;
      }
      printf("%8.4lf %8.4lf %8.4lf \n",PreOpt,costs[7],mincost);
      fprintf(fp,"%8.4lf %8.4lf %8.4lf \n",PreOpt,costs[7],mincost);
      if(PreOptFile) {
	fprintf(fpPreOpt,"%8.8lf %8.8lf \n",PreOpt,costs[7]);
      }
      nth ++;
    }
    if(PreOptFile) fclose(fpPreOpt);
    if(PreOptOnly) {
      printf("PreOptOnly specified, so exiting now\n");
      exit(0);
    }
    p[PreOptDim] = PreOptAtMin; // phase encode direction
    GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
    MatrixCopy(R0,R);
    printf("\n");
    fprintf(fp,"\n");
  }

  if(BruteForce){
    n = pow((PreOptMax-PreOptMin)/PreOptDelta+1.0,6.0);
    printf("\n");
    printf("------------------------------------\n");
    nsubsampsave = nsubsamp;
    nsubsamp = nsubsampbrute;
    printf("Brute force preopt %g %g %g, n = %d\n",PreOptMin,PreOptMax,PreOptDelta,n);
    fprintf(fp,"\n");
    fprintf(fp,"Brute force preopt %g %g %g, n = %d\n",PreOptMin,PreOptMax,PreOptDelta,n);
    fflush(stdout); fflush(fp);
    mincost = 10e10;
    PreOptAtMin = 0;
    mytimer.reset() ;
    nth = 0;
    if (PreOptDeltaTrans < 0)  // no separate range specified for translations
    {
      PreOptMinTrans = PreOptMin ; PreOptMaxTrans = PreOptMax ; PreOptDeltaTrans = PreOptDelta ;
    }
    if (UseLH && mask_label)
      LabelRipRestOfSurface(mask_label, lhwm) ;
    else if (UseRH && mask_label)
      LabelRipRestOfSurface(mask_label, rhwm) ;
    if(PreOptFile) fpPreOpt = fopen(PreOptFile,"w");
    for(tx = PreOptMinTrans; tx <= PreOptMaxTrans; tx += PreOptDeltaTrans){
      for(ty = PreOptMinTrans; ty <= PreOptMaxTrans; ty += PreOptDeltaTrans){
        for(tz = PreOptMinTrans; tz <= PreOptMaxTrans; tz += PreOptDeltaTrans){
          for(ax = PreOptMin; ax <= PreOptMax; ax += PreOptDelta){
            for(ay = PreOptMin; ay <= PreOptMax; ay += PreOptDelta){
              for(az = PreOptMin; az <= PreOptMax; az += PreOptDelta){
                p[0] = tx;
                p[1] = ty;
                p[2] = tz;
                p[3] = ax;
                p[4] = ay;
                p[5] = az;
                GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
                if(costs[7] < mincost) {
                  mincost = costs[7];
                  for(n=0; n < 6; n++) pmin[n] = p[n];
                  secCostTime = mytimer.seconds() ;
                  printf("%6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf %8.4lf %4.1f\n",
                         nth,tx,ty,tz,ax,ay,az,costs[7],mincost,secCostTime/60);
                  fprintf(fp,"%6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf %8.4lf %4.1f\n",
                          nth,tx,ty,tz,ax,ay,az,costs[7],mincost,secCostTime/60);
                  fflush(stdout); fflush(fp);
                  
                } else {
                  if(nth == 0 || nth%1000 == 0 || debug){
                    secCostTime = mytimer.seconds() ;
                    printf("%6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf %8.4lf %4.1f\n",
                           nth,tx,ty,tz,ax,ay,az,costs[7],mincost,secCostTime/60);
                    fprintf(fp,"%6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf %8.4lf %4.1f\n",
                            nth,tx,ty,tz,ax,ay,az,costs[7],mincost,secCostTime/60);
                    fflush(stdout); fflush(fp);
                  }
                }
                if(PreOptFile) 
                  fprintf(fpPreOpt,"%8.8lf %8.8lf %8.8lf %8.8lf %8.8lf %8.8lf    %8.8lf\n",
                          tx,ty,tz,ax,ay,az,costs[7]);
                nth ++;
              }
            }
          }
        }
      }
    }
    // Assign min found above to p vector
    for(n=0; n < 6; n++) p[n] = pmin[n];
    GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
    MatrixCopy(R0,R);
    printf("Brute Force --------------------------\n");
    printf("Min cost was %lf\n",mincost);
    printf("Number of iterations %5d\n",nth);
    secCostTime = mytimer.seconds() ;
    printf("Search time %lf sec\n",secCostTime);
    printf("Parameters at best (transmm, rotdeg)\n");
    printf("%7.3lf %7.3lf %7.3lf %6.3lf %6.3lf %6.3lf \n",
	   p[0],p[1],p[2],p[3],p[4],p[5]);
    printf("--------------------------------------------\n");
    
    if(PreOptFile) fclose(fpPreOpt);
    if(PreOptOnly) {
      printf("PreOptOnly specified, so exiting now\n");
      exit(0);
    }
    printf("\n");
    fprintf(fp,"\n");
    // Restore
    interpmethod = "trilinear";
    interpcode = SAMPLE_TRILINEAR;
    nsubsamp = nsubsampsave;
  }

  mytimer.reset() ;
  printf("Starting Powell Minimization\n");
  MinPowell(mov, NULL, R, p, dof, TolPowell, LinMinTolPowell,
	    nMaxItersPowell,SegRegCostFile, costs, &nth);
  secCostTime = mytimer.seconds() ;

  // Compute relative final cost 
  rcost = RelativeSurfCost(mov, R);

  // Recompute at optimal. This forces MRI *out to be the output at best reg
  if(surfcostbase){
    sprintf(tmpstr,"%s.lh.mgh",surfcostbase);
    lhcostfile = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s.rh.mgh",surfcostbase);
    rhcostfile = strcpyalloc(tmpstr);
  }
  if(surfconbase){
    sprintf(tmpstr,"%s.lh.mgh",surfconbase);
    lhconfile = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s.rh.mgh",surfconbase);
    rhconfile = strcpyalloc(tmpstr);
  }
  GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
  if(surfcostdiffbase){
    if(UseLH){
      sprintf(tmpstr,"%s.lh.mgh",surfcostdiffbase);
      printf("Writing lh diff to %s\n",tmpstr);
      mritmp = MRIsubtract(lhcost,lhcost0,NULL);
      MRIwrite(mritmp,tmpstr);
      MRIfree(&mritmp);
    }
    if(UseRH){
      sprintf(tmpstr,"%s.rh.mgh",surfcostdiffbase);
      printf("Writing rh diff to %s\n",tmpstr);
      mritmp = MRIsubtract(rhcost,rhcost0,NULL);
      MRIwrite(mritmp,tmpstr);
      MRIfree(&mritmp);
    }
  }

  if(MinCostFile){
    fpMinCost = fopen(MinCostFile,"w");
    //MinCost, WMMean, CtxMean, PctContrast
    fprintf(fpMinCost,"%lf %lf %lf %lf \n",costs[7],costs[1],costs[4],costs[6]);
    fclose(fpMinCost);
  }

  printf("Number of iterations %5d\n",nth);
  printf("Min cost was %lf\n",costs[7]);
  printf("Number of FunctionCalls %5d\n",nCostEvaluations);
  printf("TolPowell %lf\n",TolPowell);
  printf("nMaxItersPowell %d\n",nMaxItersPowell);
  printf("OptimizationTime %lf sec\n",secCostTime);
  printf("Parameters at optimum (transmm) %8.5lf %8.5lf %8.5lf\n",
	 p[0],p[1],p[2]);
  printf("Parameters at optimum (rotdeg) %8.5lf %8.5lf %8.5lf \n",
	 p[3],p[4],p[5]);
  if(dof > 6){
    printf("Parameters at optimum (scale) ");
    printf("%8.5lf %8.5lf %8.5lf\n",p[6],p[7],p[8]);
  }
  if(dof > 9){
    printf("Parameters at optimum (shear) ");
    printf("%8.5lf %8.5lf %8.5lf\n",p[9],p[10],p[11]);
  }
  if(ParamFile){
    // Write out transmm rotdeg scale shear (each with 3)
    fpParam = fopen(ParamFile,"w");
    fprintf(fpParam,"%lf %lf %lf %lf %lf %lf  ",p[0],p[1],p[2],p[3],p[4],p[5]);
    if(dof > 6) {
      fprintf(fpParam,"%lf %lf %lf ",p[6],p[7],p[8]);
      if(dof > 9) fprintf(fpParam,"%lf %lf %lf ",p[9],p[10],p[11]);
      else        fprintf(fpParam,"%lf %lf %lf ",0.0,0.0,0.0);
    }
    else fprintf(fpParam,"%lf %lf %lf %lf %lf %lf ",1.0,1.0,1.0,0.0,0.0,0.0);
    fprintf(fpParam,"\n");
    fclose(fpParam);
  }

  GetSurfCosts(mov, NULL, R0, R, p, dof, costs);
  fprintf(fp,"Final costs ----------------\n");
  fprintf(fp,"Number of surface hits %d\n",(int)costs[0]);  
  fprintf(fp,"WM  Intensity0 %10.4lf +/- %8.4lf\n",costs[1],costs[2]); 
  fprintf(fp,"Ctx Intensity0 %10.4lf +/- %8.4lf\n",costs[4],costs[5]); 
  fprintf(fp,"Pct Contrast0  %10.4lf +/- %8.4lf\n",costs[6],costs[3]); 
  fprintf(fp,"Cost %8.4lf\n",costs[7]); 
  fprintf(fp,"RelCost %8.4lf\n",rcost0);
  fflush(fp);

  printf("Final costs ----------------\n");
  printf("Number of surface hits %d\n",(int)costs[0]);  
  printf("WM  Intensity %10.4lf +/- %8.4lf\n",costs[1],costs[2]); 
  printf("Ctx Intensity %10.4lf +/- %8.4lf\n",costs[4],costs[5]); 
  printf("Pct Contrast  %10.4lf +/- %8.4lf\n",costs[6],costs[3]); 
  printf("Cost %8.4lf\n",costs[7]); 
  printf("RelCost %8.4lf\n",rcost0);

  printf("Reg at min cost was \n");
  MatrixPrint(stdout,R);
  printf("\n");
  
  fprintf(fp,"Number of iterations %5d\n",nth);
  fprintf(fp,"Min cost was %lf\n",costs[7]);
  fprintf(fp,"RelMinCost %8.4lf\n",rcost);
  fprintf(fp,"Number of FunctionCalls %5d\n",nCostEvaluations);
  fprintf(fp,"OptimizationTime %lf sec\n",secCostTime);
  fprintf(fp,"Parameters at optimum (transmm, rotdeg)\n");
  fprintf(fp,"%8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf \n",
	 p[0],p[1],p[2],p[3],p[4],p[5]);
  if(dof > 6){
    fprintf(fp,"Parameters at optimum (scale)  ");
    fprintf(fp,"%8.5lf %8.5lf %8.5lf\n",p[6],p[7],p[8]);
  }
  if(dof > 9){
    fprintf(fp,"Parameters at optimum (shear)  ");
    fprintf(fp,"%8.5lf %8.5lf %8.5lf\n",p[9],p[10],p[11]);
  }
  fprintf(fp,"Number of surface hits %d\n",(int)costs[0]);  
  fprintf(fp,"WM  Intensity %10.4lf +/- %8.4lf\n",costs[1],costs[2]); 
  fprintf(fp,"Ctx Intensity %10.4lf +/- %8.4lf\n",costs[4],costs[5]); 
  fprintf(fp,"Pct Contrast  %10.4lf +/- %8.4lf\n",costs[6],costs[3]); 
  fprintf(fp,"Cost at optimum %8.4lf\n",costs[7]); 
  fprintf(fp,"RelCostOptimum %8.4lf\n",rcost);

  fprintf(fp,"nhits=%7d wmmn=%10.4lf wmstd=%8.4lf ",
	 (int)costs[0],costs[1],costs[2]); 
  fprintf(fp,"constd=%10.4lf ctxmn=%10.4lf ctxstd=%8.4lf ",
	 costs[3],costs[4],costs[5]); // constd ctxmn ctxstd
  fprintf(fp,"conmn = %8.4lf costmn = %8.4lf ",costs[6],costs[7]); // conmn costmn
  fprintf(fp,"\n");
  fprintf(fp,
	  "#C# %7.3lf %7.3lf %7.3lf %6.3lf %6.3lf %6.3lf %lf %lf %lf %lf %lf %lf %lf %lf %5d %d\n",
	  p[0],p[1],p[2],p[3],p[4],p[5],rmsDiffMean,rmsDiffMax,costs[7],costs[1],costs[4],costs[6],
	  PenaltyCenter,PenaltySlope,nsubsamp,UseMask);
  
  fprintf(fp,"Reg at min cost was \n");
  MatrixPrint(fp,R);
  fprintf(fp,"\n");
  
  if(outregfile){
    int type = TransformFileNameType(outregfile);
    printf("Writing optimal reg to %s, type = %d \n",outregfile,type);
    fflush(stdout);
    if (type == TRANSFORM_ARRAY_TYPE) {
      LTA *lta = LTAalloc(1, NULL) ;
      LT  *lt;
      
      printf("saving transform to LTA file\n") ;
      strcpy(lta->subject, subject) ;
      lta->fscale = intensity ;
      lt = &lta->xforms[0] ;
      lt->m_L = MatrixCopy(R,NULL) ;
      getVolGeom(anat, &lt->src) ;
      getVolGeom(mov, &lt->dst) ;
      strcpy(lt->dst.fname, movvolfile) ;
      lta->type = REGISTER_DAT ;
      LTAwriteEx(lta, outregfile) ;
      LTAfree(&lta) ;
    }
    else // this is disabled
      regio_write_register(outregfile,subject,mov->xsize,
                           mov->zsize,intensity,R,FLT2INT_ROUND);

  }
  
  if(outfile) {
    // Allocate the output
    out = MRIallocSequence(anat->width,  anat->height,
			   anat->depth,  MRI_FLOAT, 1);
    MRIcopyHeader(anat,out);
    MRIvol2VolTkRegVSM(mov, out, R, SAMPLE_TRILINEAR, sinchw, vsm);
    printf("Writing output volume to %s \n",outfile);
    MRIwrite(out,outfile);
  }

  printf("Original Reg \n");
  fflush(stdout);
  MatrixPrint(stdout,R00);
  printf("\n");
  
  Rdiff = MatrixSubtract(R00,R,NULL);
  printf("Original Reg - Optimal Reg\n");
  MatrixPrint(stdout,Rdiff);
  printf("\n");

  fprintf(fp,"Original Reg - Optimal Reg\n");
  MatrixPrint(fp,Rdiff);
  fprintf(fp,"\n");

  if(DoRMSDiff){
    rmsDiffSum = 0;
    rmsDiffSum2 = 0;
    rmsDiffMax = 0;
    if(lhwm){
      printf("Computing change in lh position\n");
      MRISmatrixMultiply(lhwm,Rdiff);
      for(vno = 0; vno < lhwm->nvertices; vno++){
	v = &(lhwm->vertices[vno]);
	d = sqrt(v->x*v->x + v->y*v->y + v->z*v->z );
	rmsDiffSum += d;
	rmsDiffSum2 += (d*d);
	if(rmsDiffMax < d) rmsDiffMax = d;
      }
      rmsDiffMean = rmsDiffSum/(lhwm->nvertices);
      printf("LH rmsDiffMean %lf\n",rmsDiffMean);
    }
    if(rhwm){
      printf("Computing change in rh position\n");
      MRISmatrixMultiply(rhwm,Rdiff);
      for(vno = 0; vno < rhwm->nvertices; vno++){
	v = &(rhwm->vertices[vno]);
	d = sqrt(v->x*v->x + v->y*v->y + v->z*v->z );
	rmsDiffSum += d;
	rmsDiffSum2 += (d*d);
	if(rmsDiffMax < d) rmsDiffMax = d;
      }
      if(lhwm) n = (lhwm->nvertices + rhwm->nvertices);
      else     n = rhwm->nvertices;
      rmsDiffMean = rmsDiffSum/n;
      rmsDiffStd = sqrt((rmsDiffSum2-n*rmsDiffMean*rmsDiffMean)/(n-1));
    }
    printf("Surface-RMS-Diff-mm %lf %lf %lf\n",rmsDiffMean,rmsDiffStd,rmsDiffMax);
    fprintf(fp,"Surface-RMS-Diff-mm %lf %lf %lf\n",rmsDiffMean,rmsDiffStd,rmsDiffMax);
    if(RMSDiffFile){
      fpRMSDiff = fopen(RMSDiffFile,"w");
      //rmsDiffMean rmsDiffStd rmsDiffMax MinCost Nevals Tx Ty Tz Rx Ry Rz  WMMean CtxMean PctContrast 
      fprintf(fpRMSDiff,
	      "%12.9lf %12.9lf %12.9lf   %12.9lf %4d  %7.3lf %7.3lf %7.3lf %6.3lf %6.3lf %6.3lf  %lf %lf %lf \n",
	      rmsDiffMean,rmsDiffStd,rmsDiffMax,costs[7],nCostEvaluations,
	      p[0],p[1],p[2],p[3],p[4],p[5],
	      costs[1],costs[4],costs[6]);
	      
      fclose(fpRMSDiff);
    }    
  }

  if(PrintT1Warning){
    fprintf(fp,"\n\n");
    fprintf(fp,"WARNING: initial G-W contrast was negative, but expected positive.\n");
    fprintf(fp,"If the mov data has a T1 contrast, re-run with --T1\n");
    fprintf(fp,"\n\n");
    printf("\n\n");
    printf("WARNING: initial G-W contrast was negative, but expected positive.\n");
    printf("If the mov data has a T1 contrast, re-run with --T1\n");
    printf("\n\n");
  }
  if(PrintT2Warning){
    fprintf(fp,"\n\n");
    fprintf(fp,"WARNING: initial G-W contrast was positive, but expected negative.\n");
    fprintf(fp,"If the mov data has a T2 contrast, re-run with --T2\n");
    fprintf(fp,"\n\n");
    printf("\n\n");
    printf("WARNING: initial G-W contrast was positive, but expected negative.\n");
    printf("If the mov data has a T2 contrast, re-run with --T2\n");
    printf("\n\n");
  }


  fclose(fp);
  printf("mri_segreg done\n");

  exit(0);
  return(0);
}


/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option;
  int err,nv,n;
  double vmin, vmax, vdelta, v, c, d;
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

    if (!strcasecmp(option,      "--help"))     print_help() ;
    else if (!strcasecmp(option, "--version"))  print_version() ;
    else if (!strcasecmp(option, "--debug"))    debug = 1;
    else if (!strcasecmp(option, "--profile"))  DoProfile = 1;
    else if (!strcasecmp(option, "--abs"))       DoAbs = 1;
    else if (!strcasecmp(option, "--no-abs"))    DoAbs = 0;
    else if (!strcasecmp(option, "--no-cortex-label")) UseCortexLabel = 0;
    else if (!strcasecmp(option, "--brute_trans")){
      if (nargc < 3) argnerr(option,3);
      nargsused = 3;
      sscanf(pargv[0],"%lf",&PreOptMinTrans);
      sscanf(pargv[1],"%lf",&PreOptMaxTrans);
      sscanf(pargv[2],"%lf",&PreOptDeltaTrans);
      BruteForce = 1;
    }
    else if (!strcasecmp(option, "--brute")){
      if (nargc < 3) argnerr(option,3);
      nargsused = 3;
      sscanf(pargv[0],"%lf",&PreOptMin);
      sscanf(pargv[1],"%lf",&PreOptMax);
      sscanf(pargv[2],"%lf",&PreOptDelta);
      BruteForce = 1;
    }
    else if (!strcasecmp(option, "--gm-proj-frac")) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&GMProjFrac);
      DoGMProjFrac = 1;
      DoGMProjAbs = 0;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--wm-proj-frac")) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&WMProjFrac);
      DoWMProjFrac = 1;
      DoWMProjAbs = 0;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--gm-proj-abs")) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&GMProjAbs);
      DoGMProjAbs = 1;
      DoGMProjFrac = 0;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--wm-proj-abs")) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&WMProjAbs);
      DoWMProjAbs = 1;
      DoWMProjFrac = 0;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mid-frame")) DoMidFrame = 1;
    else if (!strcasecmp(option, "--no-mask")) UseMask = 0;
    else if (!strcasecmp(option, "--mask"))    UseMask = 1;
    else if (!strcasecmp(option, "--lh-mask")){
      if(nargc < 1) argnerr(option,1);
      lhsegmask = MRIread(pargv[0]);
      if(lhsegmask==NULL) exit(1);
      UseMask = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--rh-mask")){
      if(nargc < 1) argnerr(option,1);
      rhsegmask = MRIread(pargv[0]);
      if(rhsegmask==NULL) exit(1);
      UseMask = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--sd")) {
      if(nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0] ;
      printf("using %s as SUBJECTS_DIR\n", SUBJECTS_DIR) ;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--label")) {
      mask_label = LabelRead(NULL, pargv[0]) ;
      if (mask_label == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname, pargv[0]) ;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--surf")) {
      surfname = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--lh-only")){ UseLH = 1; UseRH = 0;}
    else if (!strcasecmp(option, "--rh-only")){ UseLH = 0; UseRH = 1;}
    else if (istringnmatch(option, "--surf",0)) {
      if(nargc < 1) argnerr(option,1);
      surfname = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--surf-cost",0)) {
      if(nargc < 1) argnerr(option,1);
      surfcostbase = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--surf-con",0)) {
      if(nargc < 1) argnerr(option,1);
      surfconbase = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--init-surf-cost",0)) {
      if(nargc < 1) argnerr(option,1);
      surfcost0base = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--init-surf-cost-only",0)) {
      InitSurfCostOnly=1;
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--surf-cost-diff",0)) {
      if(nargc < 1) argnerr(option,1);
      surfcostdiffbase = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--mov",0)) {
      if (nargc < 1) argnerr(option,1);
      movvolfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--vsm",0)) {
      if (nargc < 1) argnerr(option,1);
      vsmfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--frame",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--dof",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&dof);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--9",0)) {
      dof = 9;
    } 
    else if (istringnmatch(option, "--n1dmin",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&n1dmin);
      nargsused = 1;
    } else if (istringnmatch(option, "--1dpreopt",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&PreOptMin);
      sscanf(pargv[1],"%lf",&PreOptMax);
      sscanf(pargv[2],"%lf",&PreOptDelta);
      Do1DPreOpt = 1;
      nargsused = 3;
    } 
    else if (istringnmatch(option, "--preopt-file",0)) {
      if(nargc < 1) argnerr(option,1);
      PreOptFile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--preopt-only",0)) PreOptOnly = 1;
    else if (istringnmatch(option, "--preopt-dim",0)) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&PreOptDim);
      nargsused = 1;
    } 
    else if(istringnmatch(option, "--nsub",0) ||
	    istringnmatch(option, "--subsamp",0) ||
	    istringnmatch(option, "--skip",0)) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nsubsamp);
      nargsused = 1;
    } 
    else if(istringnmatch(option, "--subsamp-brute",0)){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nsubsampbrute);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--nmax",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nMaxItersPowell);
      nargsused = 1;
    } else if (istringnmatch(option, "--tol",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&TolPowell);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--tol1d",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&LinMinTolPowell);
      nargsused = 1;
    }
    else if (istringnmatch(option, "--o",0)) {
      if (nargc < 1) argnerr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--cost-eval",0)) {
      if (nargc < 1) argnerr(option,1);
      costfile_powell = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--relcost",0)) {
      if (nargc < 1) argnerr(option,1);
      RelCostFile = pargv[0];
      nargsused = 1;
    } 
    else if(istringnmatch(option, "--init-reg",0) ||
	    istringnmatch(option, "--reg",0)) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      err = regio_read_register(regfile, &subject, &ipr, &bpr,
                                &intensity, &R0, &float2int);
      if (err) exit(1);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--regheader")){
      if(nargc < 1) argnerr(option,1);
      regheader = 1;
      subject = pargv[0];
      nargsused = 1;
    }
    else if (istringnmatch(option, "--rot",0)) {
      if (nargc < 3) argnerr(option,3);
      // Angles are in degrees
      sscanf(pargv[0],"%lf",&angles[0]);
      sscanf(pargv[1],"%lf",&angles[1]);
      sscanf(pargv[2],"%lf",&angles[2]);
      angles[0] *= (M_PI/180);
      angles[1] *= (M_PI/180);
      angles[2] *= (M_PI/180);
      MrotPre = MRIangles2RotMat(angles);
      nargsused = 3;
    } 
    else if (istringnmatch(option, "--rot-rand",0)) {
      if(nargc < 1) argnerr(option,1);
      // Rotation in deg
      sscanf(pargv[0],"%lf",&RotRandMax);
      angles[0] = 2.0*(drand48()-0.5)*RotRandMax*(M_PI/180);
      angles[1] = 2.0*(drand48()-0.5)*RotRandMax*(M_PI/180);
      angles[2] = 2.0*(drand48()-0.5)*RotRandMax*(M_PI/180);
      MrotPre = MRIangles2RotMat(angles);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--trans",0)) {
      if (nargc < 3) argnerr(option,3);
      // Translation in mm
      sscanf(pargv[0],"%lf",&xyztrans[0]);
      sscanf(pargv[1],"%lf",&xyztrans[1]);
      sscanf(pargv[2],"%lf",&xyztrans[2]);
      MtransPre = MatrixIdentity(4,NULL);
      MtransPre->rptr[1][4] = xyztrans[0];
      MtransPre->rptr[2][4] = xyztrans[1];
      MtransPre->rptr[3][4] = xyztrans[2];
      nargsused = 3;
    } 
    else if (istringnmatch(option, "--scale",0)) {
      if (nargc < 3) argnerr(option,3);
      // Scale
      sscanf(pargv[0],"%lf",&scale[0]);
      sscanf(pargv[1],"%lf",&scale[1]);
      sscanf(pargv[2],"%lf",&scale[2]);
      MscalePre = MatrixIdentity(4,NULL);
      MscalePre->rptr[1][1] = scale[0];
      MscalePre->rptr[2][2] = scale[1];
      MscalePre->rptr[3][3] = scale[2];
      nargsused = 3;
    } 
    else if (istringnmatch(option, "--shear",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&shear[0]);
      sscanf(pargv[1],"%lf",&shear[1]);
      sscanf(pargv[2],"%lf",&shear[2]);
      MshearPre = MatrixIdentity(4,NULL);
      MshearPre->rptr[1][2] = shear[0];
      MshearPre->rptr[1][3] = shear[1];
      MshearPre->rptr[2][3] = shear[2];
      nargsused = 3;
    } 
    else if (istringnmatch(option, "--trans-rand",0)) {
      if(nargc < 1) argnerr(option,1);
      // Translation in mm
      sscanf(pargv[0],"%lf",&TransRandMax);
      xyztrans[0] = 2.0*(drand48()-0.5)*TransRandMax;
      xyztrans[1] = 2.0*(drand48()-0.5)*TransRandMax;
      xyztrans[2] = 2.0*(drand48()-0.5)*TransRandMax;
      MtransPre = MatrixIdentity(4,NULL);
      MtransPre->rptr[1][4] = xyztrans[0];
      MtransPre->rptr[2][4] = xyztrans[1];
      MtransPre->rptr[3][4] = xyztrans[2];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--out-reg",0)) {
      if (nargc < 1) argnerr(option,1);
      outregfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--cur-reg",0)) {
      if (nargc < 1) argnerr(option,1);
      curregfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--sum",0)) {
      if (nargc < 1) argnerr(option,1);
      sumfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--gm-gt-wm",0)) {
      if (nargc < 1) argnerr(option,1);
      PenaltySign = -1;
      sscanf(pargv[0],"%lf",&PenaltySlope);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--wm-gt-gm",0)) {
      if (nargc < 1) argnerr(option,1);
      PenaltySign = +1;
      sscanf(pargv[0],"%lf",&PenaltySlope);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--slope",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&PenaltySlope);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--T1",0) ||
	     istringnmatch(option, "--t1",0))   PenaltySign = +1;
    else if (istringnmatch(option, "--T2",0) ||
	     istringnmatch(option, "--t2",0) ||
	     istringnmatch(option, "--bold",0)) PenaltySign = -1;
    else if (istringnmatch(option, "--c0",0) ||
	     istringnmatch(option, "--offset",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&PenaltyCenter);
      nargsused = 1;
    } else if (istringnmatch(option, "--wm-gt-gm",0)) {
      if (nargc < 1) argnerr(option,1);
      PenaltySign = +1;
      sscanf(pargv[0],"%lf",&PenaltySlope);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--penalty-abs",0)) {
      // no direction of contrast expected
      PenaltySign = 0;
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--ignore-neg",0)) {
      PenaltySign = -2;
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--fwhm",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&fwhm);
      gstd = fwhm/sqrt(log(256.0));
      DoSmooth = 1;
      nargsused = 1;
    } else if (istringnmatch(option, "--noise",0)) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&NoiseStd);
      AddNoise = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--seed")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&SynthSeed);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--tx-mmd",0)) {
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("ntx = %d  %lf %lf %lf\n",nv,vmin,vmax,vdelta);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
        txlist[ntx] = vmin + vdelta*n;
        ntx++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--ty-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("nty = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
        tylist[nty] = vmin + vdelta*n;
        nty++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--tz-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("ntz = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
        tzlist[ntz] = vmin + vdelta*n;
        ntz++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--ax-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("nax = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
        axlist[nax] = vmin + vdelta*n;
        nax++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--ay-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("nay = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
        aylist[nay] = vmin + vdelta*n;
        nay++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--az-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("naz = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
        azlist[naz] = vmin + vdelta*n;
        naz++;
      }
      nargsused = 3;
    } else if ( !strcmp(option, "--gdiagno") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&gdiagno);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--cost",0)) {
      if (nargc < 1) argnerr(option,1);
      SegRegCostFile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--targ-con",0)) {
      if(nargc < 2) argnerr(option,2);
      TargConLHFile = pargv[0];
      TargConRHFile = pargv[1];
      TargConLH = MRIread(TargConLHFile);
      TargConRH = MRIread(TargConRHFile);
      nargsused = 2;
    } 
    else if (istringnmatch(option, "--mincost",0)) {
      if (nargc < 1) argnerr(option,1);
      MinCostFile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--initcost",0)) {
      if (nargc < 1) argnerr(option,1);
      InitCostFile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--param",0)) {
      if (nargc < 1) argnerr(option,1);
      ParamFile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--rms",0)) {
      if (nargc < 1) argnerr(option,1);
      RMSDiffFile = pargv[0];
      DoRMSDiff = 1;
      nargsused = 1;
    } else if (istringnmatch(option, "--cf",0)) {
      if (nargc < 1) argnerr(option,1);
      nargsused = 1;
      fp = fopen(pargv[0],"w");
      for(v = -30; v < +30; v++){
	c = VertexCost(10*(1+v/200)/(1-v/200), 10, PenaltySlope, PenaltyCenter, PenaltySign, &d);
	fprintf(fp,"%lf %lf\n",d,c);
      }
      fclose(fp);
      exit(0);
    } else if (istringnmatch(option, "--interp",8)) {
      if (nargc < 1) argnerr(option,1);
      interpmethod = pargv[0];
      nargsused = 1;
      if (!strcmp(interpmethod,"sinc") && CMDnthIsArg(nargc, pargv, 1)) {
        sscanf(pargv[1],"%d",&sinchw);
        nargsused ++;
      }
    } 
    else if (istringnmatch(option, "--trilinear",6)) {
      interpmethod = "trilinear";
      interpcode = SAMPLE_TRILINEAR;
    }
    else if (istringnmatch(option, "--nearest",7)) {
      interpmethod = "nearest";
      interpcode = SAMPLE_NEAREST;
    }
    else {
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
static void print_usage(void) 
{
printf("\n");
printf("  mri_segreg\n");
printf("\n");
printf("  --init-reg regfile\n");
printf("  --regheader subject\n");
printf("  --mov fvol\n");
printf("  --out-reg outreg : reg at lowest cost\n");
printf("\n");
printf("  --cost costfile\n");
printf("  --sum sumfile : default is outreg.sum\n");
printf("    --cur-reg curregfile : reg at current optimum\n");
printf("  --o out : save final output\n");
printf("\n");
printf("  --brute_trans min max delta : brute force translation in all directions\n");
printf("  --brute min max delta : brute force in all directions\n");
printf("  --1dpreopt min max delta : brute force in PE direction\n");
printf("\n");
printf("  --fwhm fwhm : smooth input by fwhm mm\n");
printf("  --abs       : compute abs of mov\n");
printf("  --subsamp nsub : only sample every nsub vertices\n");
printf("  --subsamp-brute nsub : only sample every nsub vertices during brute-force search (%d)\n",nsubsampbrute);
printf("\n");
printf("  --preopt-file file : save preopt results in file\n");
printf("  --preopt-dim dim : 0-5 (def 2) (0=TrLR,1=TrSI,2=TrAP,3=RotLR,4=RotSI,5=RotAP)\n");
printf("  --preopt-only : only preopt, so not optimize\n");
printf("\n");
printf("  --T1, --t1         : assume T1 gray/white contrast\n");
printf("  --T2, --t2, --bold : assume T2/BOLD gray/white contrast (default)\n");
printf("  --slope slope    : set cost slope\n");
printf("  --gm-gt-wm slope : set cost slope, spec gray matter brighter than WM\n");
printf("  --wm-gt-gm slope : set cost slope, spec WM brighter than gray matter\n");
printf("  --penalty-abs : remove sign from contrast\n");
printf("  --c0 offset : cost offset (pct)\n");
printf("  --cf cfile  : save cost function values (pct,cost)\n");
printf("\n");
printf("  --label <label file> : only use the portion of the surface in the label (only for --?h-only\n");
printf("\n");
printf("  --mask : mask out expected B0 regions\n");
printf("  --no-mask : do not mask out (default)\n");
printf("  --no-cortex-label : do not use cortex label\n");
printf("\n");
printf("  --lh-only : only use left hemisphere\n");
printf("  --rh-only : only use right hemisphere\n");
printf("\n");
printf("  --gm-proj-frac frac : fraction of cortical thickness (default, 0.5)\n");
printf("  --gm-proj-abs  dist : absolute distance into cortex\n");
printf("\n");
printf("  --wm-proj-frac frac : fraction of cortical thickness\n");
printf("  --wm-proj-abs  dist : absolute distance into WM (default 2mm)\n");
printf("\n");
printf("  --projabs : project -1mm and + 2mm\n");
printf("  --proj-frac frac : projection fraction\n");
printf("\n");
printf("  --frame nthframe : use given frame in input (default = 0)\n");
printf("  --mid-frame : use use middle frame\n");
printf("\n");
printf("  --trans Tx Ty Tz : translate input reg (mm)\n");
printf("  --rot   Ax Ay Zz : rotate input reg (degrees)\n");
printf("\n");
printf("  --trans-rand Tmax : uniformly dist -Tmax to +Tmax (Tx,Ty,Tz)\n");
printf("  --rot-rand   Amax : uniformly dist -Amax to +Amax (Ax,Ay,Az)\n");
printf("\n");
printf("  --interp interptype : interpolation trilinear or nearest (def is trilin)\n");
printf("  --profile : print out info about exec time\n");
printf("\n");
printf("  --noise stddev : add noise with stddev to input for testing sensitivity\n");
printf("  --seed randseed : for use with --noise\n");
printf("\n");
printf("  --aseg : use aseg instead of segreg.mgz\n");
printf("\n");
printf("  --nmax nmax   : max number of powell iterations (def 36)\n");
printf("  --tol   tol   : powell inter-iteration tolerance on cost\n");
printf("       This is the fraction of the cost that the difference in \n");
printf("       successive costs must drop below to stop the optimization.  \n");
printf("  --tol1d tol1d : tolerance on powell 1d minimizations\n");
printf("\n");
printf("  --1dmin : use brute force 1D minimizations instead of powell\n");
printf("  --n1dmin n1dmin : number of 1d minimization (default = 3)\n");
printf("\n");
printf("  --mincost MinCostFile\n");
printf("  --initcost InitCostFile\n");
printf("  --param   ParamFile\n");
printf("  --rms     RMSDiffFile : saves Tx Ty Tz Ax Ay Az RMSDiff MinCost \n");
printf("              WMMean CtxMean PctContrast C0 Slope NSubSamp UseMask\n");
printf("  --surf surfname : use ?h.surfname instead of lh.white\n");
printf("  --surf-cost basename : saves final cost as basename.?h.mgh\n");
printf("  --init-surf-cost basename0 : saves init cost as basename0.?h.mgh\n");
printf("  --surf-cost-diff diffbase : saves final-init cost as diffbase.?h.mgh\n");
printf("  --surf-con basename : saves final contrast as basename.?h.mgh\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n%s\n\n",getVersion().c_str());
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) 
{
  if (SUBJECTS_DIR == NULL)
  {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL) {
      printf("ERROR: SUBJECTS_DIR undefined.\n");
      exit(1);
    }
  }

  if (movvolfile == NULL) {
    printf("ERROR: No mov volume supplied.\n");
    exit(1);
  }
  if(regfile == NULL && ! regheader) {
    printf("ERROR: need --reg or --regheader.\n");
    exit(1);
  }
  if(outregfile == NULL) {
    printf("ERROR: need --out-reg.\n");
    exit(1);
  }
  if(SegRegCostFile == NULL) {
    sprintf(tmpstr,"%s.cost",outregfile);
    SegRegCostFile = strcpyalloc(tmpstr);
  }

  interpcode = MRIinterpCode(interpmethod);
  if (interpcode < 0) {
    printf("ERROR: interpolation method %s unrecognized\n",interpmethod);
    printf("       legal values are nearest, trilin, and sinc\n");
    exit(1);
  }

  if(sumfile == NULL) {
    sprintf(tmpstr,"%s.sum",outregfile);
    sumfile = strcpyalloc(tmpstr);
  }

  if(ntx == 0) {ntx=1; txlist[0] = 0;}
  if(nty == 0) {nty=1; tylist[0] = 0;}
  if(ntz == 0) {ntz=1; tzlist[0] = 0;}
  if(nax == 0) {nax=1; axlist[0] = 0;}
  if(nay == 0) {nay=1; aylist[0] = 0;}
  if(naz == 0) {naz=1; azlist[0] = 0;}

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  int n;
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"setenv SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"cd %s\n",getcwd(tmpstr,sizeof(tmpstr)));
  fprintf(fp,"%s\n",cmdline2);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"movvol %s\n",movvolfile);
  fprintf(fp,"regfile %s\n",regfile);
  fprintf(fp,"subject %s\n",subject);
  fprintf(fp,"dof %d\n",dof);
  if(outregfile) fprintf(fp,"outregfile %s\n",outregfile);
  if(outfile) fprintf(fp,"outfile %s\n",outfile);
  fprintf(fp,"UseMask %d\n",UseMask);
  fprintf(fp,"UseLH %d\n",UseLH);
  fprintf(fp,"UseRH %d\n",UseRH);
  fprintf(fp,"nsubsamp %d\n",nsubsamp);
  fprintf(fp,"PenaltySign  %d\n",PenaltySign);
  fprintf(fp,"PenaltySlope %lf\n",PenaltySlope);
  fprintf(fp,"PenaltyCenter %lf\n",PenaltyCenter);
  fprintf(fp,"surfname %s\n",surfname);
  if(DoGMProjFrac) fprintf(fp,"GMProjFrac %lf\n",GMProjFrac);
  if(DoGMProjAbs)  fprintf(fp,"GMProjAbs %lf\n",GMProjAbs);
  if(DoWMProjFrac) fprintf(fp,"WMProjFrac %lf\n",WMProjFrac);
  if(DoWMProjAbs)  fprintf(fp,"WMProjAbs %lf\n",WMProjAbs);
  fprintf(fp,"lhcostfile %s\n",lhcostfile);
  fprintf(fp,"rhcostfile %s\n",rhcostfile);
  fprintf(fp,"interp  %s (%d)\n",interpmethod,interpcode);
  fprintf(fp,"frame  %d\n",frame);
  fprintf(fp,"TolPowell %lf\n",TolPowell);
  fprintf(fp,"nMaxItersPowell %d\n",nMaxItersPowell);
  fprintf(fp,"n1dmin  %d\n",n1dmin);
  if(interpcode == SAMPLE_SINC) fprintf(fp,"sinc hw  %d\n",sinchw);
  fprintf(fp,"Profile   %d\n",DoProfile);
  fprintf(fp,"Gdiag_no  %d\n",Gdiag_no);
  fprintf(fp,"AddNoise  %d (%g)\n",AddNoise,NoiseStd);
  fprintf(fp,"SynthSeed %d\n",SynthSeed);
  fprintf(fp,"TransRandMax %lf\n",TransRandMax);
  fprintf(fp,"RotRandMax %lf\n",RotRandMax);
  fprintf(fp,"Translations %lf %lf %lf\n", xyztrans[0],xyztrans[1],xyztrans[2]);
  fprintf(fp,"Rotations   %lf %lf %lf\n",angles[0],angles[1],angles[2]);

  fprintf(fp,"Input reg\n");
  MatrixPrint(fp,R0);
  fprintf(fp,"\n");

  if(0){
    fprintf(fp,"ntx %d\n",ntx);
    if(ntx > 0){
      fprintf(fp," tx values\n");
      for(n=0; n < ntx; n++) printf("    %2d %g\n",n+1,txlist[n]);
    }
    fprintf(fp,"nty %d\n",nty);
    if(nty > 0){
      fprintf(fp," ty values\n");
      for(n=0; n < nty; n++) printf("    %2d %g\n",n+1,tylist[n]);
    }
    fprintf(fp,"ntz %d\n",ntz);
    if(ntz > 0){
      fprintf(fp," tz values\n");
      for(n=0; n < ntz; n++) printf("    %2d %g\n",n+1,tzlist[n]);
    }
    fprintf(fp,"nax %d\n",nax);
    if(nax > 0){
      fprintf(fp," ax values\n");
      for(n=0; n < nax; n++) printf("    %2d %g\n",n+1,axlist[n]);
    }
    fprintf(fp,"nay %d\n",nay);
    if(nay > 0){
      fprintf(fp," ay values\n");
      for(n=0; n < nay; n++) printf("    %2d %g\n",n+1,aylist[n]);
    }
    fprintf(fp,"naz %d\n",naz);
    if(naz > 0){
      fprintf(fp," az values\n");
      for(n=0; n < naz; n++) printf("    %2d %g\n",n+1,azlist[n]);
    }
  }
  return;
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
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);
  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*------------------------------------------------------------
  istringnmatch() - compare the first n characters of two strings,
  return a 1 if they match (ignoring case), a zero otherwise. If
  n=0, then do a full comparison.
  ------------------------------------------------------------*/
static int istringnmatch(const char *str1, const char *str2, int n) {
  if (n > 0  && ! strncasecmp(str1,str2,n)) return(1);
  if (n <= 0 && ! strcasecmp(str1,str2)) return(1);
  return(0);
}

/*---------------------------------------------------------*/
float compute_powell_cost(float *p) 
{
  extern MRI *mov;
  extern char *costfile_powell;
  extern int nCostEvaluations;
  extern int dof;
  static MATRIX *R = NULL;
  static double copt = -1;
  static double cprev = -1;
  double costs[8], pp[12], cdelta;
  int n, newopt;
  FILE *fp;

  if(R==NULL) R = MatrixAlloc(4,4,MATRIX_REAL);
  for(n=0; n < dof; n++) pp[n] = p[n+1];
  
  GetSurfCosts(mov, NULL, R0, R, pp, dof, costs);

  // This is for a fast check on convergence
  //costs[7] = 0;
  //for(n=0; n < 6; n++) costs[7] += ((pp[n]-n)*(pp[n]-n)+1);

  if(copt == -1) copt = costs[7];
  newopt = 0;
  if(copt > costs[7]){
    copt = costs[7];
    newopt = 1;
  }

  cdelta = 0;
  if(cprev < 0) cprev = costs[7];
  cdelta = 0.5*fabs(costs[7]-cprev)/(costs[7]+cprev);

  if(costfile_powell != NULL){
    // write costs to file
    if(nCostEvaluations == 0) fp = fopen(costfile_powell,"w");
    else                      fp = fopen(costfile_powell,"a");
    fprintf(fp,"%4d ",nCostEvaluations);
    fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[0],pp[1],pp[2]);
    fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
    if(dof > 6) fprintf(fp,"sc: %4.3lf %4.3lf %4.3lf ",pp[6],pp[7],pp[8]);
    if(dof > 9) fprintf(fp,"sh: %6.3lf %6.3lf %6.3lf ",pp[9],pp[10],pp[11]);
    fprintf(fp,"  %8.5lf %8.5lf %7d\n",costs[7],copt,(int)costs[0]);
    //fprintf(fp,"%7d %10.4lf %8.4lf ",
	//    (int)costs[0],costs[1],costs[2]); // WM  n mean std
    //fprintf(fp,"%10.4lf %10.4lf %8.4lf ",
	//    costs[3],costs[4],costs[5]); // CTX n mean std
    //fprintf(fp,"%8.4lf %12.10lf   %12.10lf ",costs[6],costs[7],copt); // t, cost=1/t
    //fprintf(fp,"\n");
    fclose(fp);
  }

  if(newopt){
    // If there is a new optimum, print it out
    fp = stdout;
    fprintf(fp,"%4d ",nCostEvaluations);
    fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[0],pp[1],pp[2]);
    fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
    if(dof > 6) fprintf(fp,"sc: %4.3lf %4.3lf %4.3lf ",pp[6],pp[7],pp[8]);
    if(dof > 9) fprintf(fp,"sh: %6.3lf %6.3lf %6.3lf ",pp[9],pp[10],pp[11]);
    fprintf(fp,"  %12.10lf\n",costs[7]);
    fflush(stdout); 
    // write out the reg current from the current opt (because I'm impatient)
    if(curregfile)
      regio_write_register(curregfile,subject,mov->xsize,
			   mov->zsize,intensity,R,FLT2INT_ROUND);
  }

  nCostEvaluations++;
  cprev = costs[7];
  return((float)costs[7]);
}
/*---------------------------------------------------------------
  MRISbbrSurfs() - creates surfaces used with BBR by projecting
  white in and out. Can also create B0 mask.
  ---------------------------------------------------------------*/
int MRISbbrSurfs(char *subject)
{
  extern MRI *lhsegmask, *rhsegmask;
  extern MRIS *lhwm, *rhwm, *lhctx, *rhctx;
  extern int UseLH, UseRH, UseMask;

  extern int DoWMProjFrac;
  extern double WMProjFrac;
  extern int DoGMProjFrac;
  extern double GMProjFrac;
  extern int DoWMProjAbs;
  extern double WMProjAbs;
  extern int DoGMProjAbs;
  extern double GMProjAbs;

  //  char *SUBJECTS_DIR;
  char tmpstr[2000];
  int c,n,v;
  float  fx, fy, fz;
  int annot,B0Annots[10], nB0Annots=10;

  //  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  if(UseLH){
    printf("Projecting LH Surfs\n");
    // Load the LH white surface, project it into WM and Ctx
    printf("Loading lh.%s surf\n",surfname);
    sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,surfname);
    lhwm = MRISread(tmpstr); // starts as white, projected in
    if(lhwm == NULL) exit(1);
    lhctx = MRISread(tmpstr); // starts as white, projected out

    if(DoWMProjFrac) {
      printf("Loading lh.thickness for WM\n");
      sprintf(tmpstr,"%s/%s/surf/lh.thickness",SUBJECTS_DIR,subject);
      err = MRISreadCurvatureFile(lhwm, tmpstr);
      if(err) exit(1);
    }
    if(DoGMProjFrac) {
      printf("Loading lh.thickness for GM\n");
      sprintf(tmpstr,"%s/%s/surf/lh.thickness",SUBJECTS_DIR,subject);
      err = MRISreadCurvatureFile(lhctx, tmpstr);
      if(err) exit(1);
    }

    printf("GM Proj: %d %lf %lf\n",DoGMProjFrac,GMProjFrac,GMProjAbs);
    printf("WM Proj: %d %lf %lf\n",DoWMProjFrac,WMProjFrac,WMProjAbs);

    MRISfreeDistsButNotOrig(lhwm);
    MRISfreeDistsButNotOrig(lhctx);
      // MRISsetXYZ will invalidate all of these,
      // so make sure they are recomputed before being used again!

    for(n = 0; n < lhwm->nvertices; n++){
      if(DoWMProjAbs)  ProjNormDist(&fx, &fy, &fz, lhwm,  n, -WMProjAbs);
      if(DoWMProjFrac) ProjNormFracThick(&fx, &fy, &fz, lhwm,  n, -WMProjFrac);
      MRISsetXYZ(lhwm,n,fx,fy,fz);
      if(DoGMProjAbs)  ProjNormDist(&fx, &fy, &fz, lhctx,  n, +GMProjAbs);
      if(DoGMProjFrac) ProjNormFracThick(&fx, &fy, &fz, lhctx,  n, +GMProjFrac);
      MRISsetXYZ(lhctx,n,fx,fy,fz);
    }
    if (UseLabel)
    {
      int         vno ;
      MRI_SURFACE *mris ;
      VERTEX      *v ;
      MRI         *mri ;

      if (UseRH)
      {
        mris = rhwm ;
        mri = rhlabel = MRIalloc(rhwm->nvertices,1,1,MRI_INT);
      }
      else
      {
        mris = lhwm ;
        mri = rhlabel = MRIalloc(rhwm->nvertices,1,1,MRI_INT);
      }
      
      MRISclearMarks(mris) ;
      LabelMark(mask_label, mris) ;
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        MRIsetVoxVal(mri,n,0,0,0, v->marked);
      }
    }

    if(UseMask && lhsegmask==NULL){
      sprintf(tmpstr,"%s/%s/label/lh.aparc.annot",SUBJECTS_DIR,subject);
      printf("Reading %s\n",tmpstr);
      err = MRISreadAnnotation(lhwm, tmpstr);
      if(err) exit(1);
      B0Annots[0] = CTABentryNameToAnnotation("middletemporal",lhwm->ct);
      B0Annots[1] = CTABentryNameToAnnotation("inferiortemporal",lhwm->ct);
      B0Annots[2] = CTABentryNameToAnnotation("temporalpole",lhwm->ct);
      B0Annots[3] = CTABentryNameToAnnotation("fusiform",lhwm->ct);
      B0Annots[4] = CTABentryNameToAnnotation("entorhinal",lhwm->ct);
      B0Annots[5] = CTABentryNameToAnnotation("medialorbitofrontal",lhwm->ct);
      B0Annots[6] = CTABentryNameToAnnotation("caudalanteriorcingulate",lhwm->ct);
      B0Annots[7] = CTABentryNameToAnnotation("rostralanteriorcingulate",lhwm->ct);
      B0Annots[8] = CTABentryNameToAnnotation("unknown",lhwm->ct);
      B0Annots[9] = CTABentryNameToAnnotation("corpuscallosum",lhwm->ct);
      lhsegmask = MRIalloc(lhwm->nvertices,1,1,MRI_INT);
      for(n = 0; n < lhwm->nvertices; n++){
        annot = lhwm->vertices[n].annotation;
        v = 1;
        for(c=0; c < nB0Annots; c++){
          if(annot == B0Annots[c]){
            v = 0; 
            break;
          }
        }
        MRIsetVoxVal(lhsegmask,n,0,0,0, v);
      }
      //MRIwrite(lhsegmask,"lh.segmask.mgh");
    }
  }

  if(UseRH){
    printf("Projecting RH Surfs\n");
    // Load the RH white surface, project it into WM and Ctx
    printf("Loading rh.%s surf\n",surfname);
    sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,surfname);
    rhwm = MRISread(tmpstr);
    if(rhwm == NULL) exit(1);
    rhctx = MRISread(tmpstr);
    
    if(DoWMProjFrac) {
      printf("Loading rh.thickness for WM\n");
      sprintf(tmpstr,"%s/%s/surf/rh.thickness",SUBJECTS_DIR,subject);
      err = MRISreadCurvatureFile(rhwm, tmpstr);
      if(err) exit(1);
    }
    if(DoGMProjFrac){
      printf("Loading rh.thickness for GM\n");
      sprintf(tmpstr,"%s/%s/surf/rh.thickness",SUBJECTS_DIR,subject);
      err = MRISreadCurvatureFile(rhctx, tmpstr);
      if(err) exit(1);
    }
    
    printf("Projecting RH Surfs\n");
    MRISfreeDistsButNotOrig(rhwm);
    MRISfreeDistsButNotOrig(rhctx);
      // MRISsetXYZ will invalidate all of these,
      // so make sure they are recomputed before being used again!
    
    for(n = 0; n < rhwm->nvertices; n++){
      if(DoWMProjAbs)  ProjNormDist(&fx, &fy, &fz, rhwm,  n, -WMProjAbs);
      if(DoWMProjFrac) ProjNormFracThick(&fx, &fy, &fz, rhwm,  n, -WMProjFrac);
      MRISsetXYZ(rhwm,n,fx,fy,fz);
      if(DoGMProjAbs)  ProjNormDist(&fx, &fy, &fz, rhctx,  n, +GMProjAbs);
      if(DoGMProjFrac) ProjNormFracThick(&fx, &fy, &fz, rhctx,  n, +GMProjFrac);
      MRISsetXYZ(rhctx,n,fx,fy,fz);
    }
    if(UseMask && rhsegmask==NULL){
      sprintf(tmpstr,"%s/%s/label/rh.aparc.annot",SUBJECTS_DIR,subject);
      printf("Reading %s\n",tmpstr);
      err = MRISreadAnnotation(rhwm, tmpstr);
      if(err) exit(1);
      B0Annots[0] = CTABentryNameToAnnotation("middletemporal",lhwm->ct);
      B0Annots[1] = CTABentryNameToAnnotation("inferiortemporal",lhwm->ct);
      B0Annots[2] = CTABentryNameToAnnotation("temporalpole",lhwm->ct);
      B0Annots[3] = CTABentryNameToAnnotation("fusiform",lhwm->ct);
      B0Annots[4] = CTABentryNameToAnnotation("entorhinal",lhwm->ct);
      B0Annots[5] = CTABentryNameToAnnotation("medialorbitofrontal",lhwm->ct);
      B0Annots[6] = CTABentryNameToAnnotation("caudalanteriorcingulate",lhwm->ct);
      B0Annots[7] = CTABentryNameToAnnotation("rostralanteriorcingulate",lhwm->ct);
      B0Annots[8] = CTABentryNameToAnnotation("unknown",lhwm->ct);
      B0Annots[9] = CTABentryNameToAnnotation("corpuscallosum",lhwm->ct);
      rhsegmask = MRIalloc(rhwm->nvertices,1,1,MRI_INT);
      for(n = 0; n < rhwm->nvertices; n++){
	annot = rhwm->vertices[n].annotation;
	v = 1;
	for(c=0; c < nB0Annots; c++){
	  if(annot == B0Annots[c]){
	    v = 0; 
	    break;
	  }
	}
	MRIsetVoxVal(rhsegmask,n,0,0,0, v);
      }
      //MRIwrite(rhsegmask,"rh.segmask.mgh");
    }
  }

  return(0);
}

/*------------------------------------------------------
  vctx - intensity in cortex
  vwm - intensity in wm
  slope - cost function slope at center
  center - cost function center (percent contrast)
  sign - 0 (abs), +1 (gm>wm), -1 (gm<wm)
  --------------------------------------------------------*/
double VertexCost(double vctx, double vwm, double slope, 
		  double center, double sign, double *pct)
{
  double d,a=0,c;
  d = 100*(vctx-vwm)/((vctx+vwm)/2.0); // percent contrast
  if(sign ==  0) a = -fabs(slope*(d-center)); // not sure this is useful
  if(sign == -1) a = -(slope*(d-center));
  if(sign == +1) a = +(slope*(d-center));
  if(sign == -2){
    if(d >= 0) a = -(slope*(d-center));
    else       a = 0;
  }
  c = 1+tanh(a);
  *pct = d;
  return(c);
}

/*-------------------------------------------------------*/
double *GetSurfCosts(MRI *mov, MRI *notused, MATRIX *R0, MATRIX *R,
		     double *p, int dof, double *costs)
{
  static MRI *vlhwm=NULL, *vlhctx=NULL, *vrhwm=NULL, *vrhctx=NULL;
  extern MRI *lhcost, *rhcost;
  extern MRI *lhcon, *rhcon;
  extern char *lhcostfile, *rhcostfile;
  extern char *lhconfile, *rhconfile;
  extern int UseMask, UseLH, UseRH;
  extern MRI *lhsegmask, *rhsegmask;
  extern MRI *lhCortexLabel, *rhCortexLabel;
  extern MRIS *lhwm, *rhwm, *lhctx, *rhctx;
  extern int PenaltySign;
  extern double PenaltySlope;
  extern int nsubsamp;
  extern int interpcode;
  double angles[3],d,dsum,dsum2,dstd,dmean,vwm,vctx,c,csum,csum2,cstd,cmean,val;
  MATRIX *Mrot=NULL, *Mtrans=NULL, *Mscale=NULL, *Mshear=NULL;
  int nhits,n;
  //FILE *fp;

  if(R==NULL){
    printf("ERROR: GetSurfCosts(): R cannot be NULL\n");
    return(NULL);
  }

  Mtrans = MatrixIdentity(4,NULL);
  if(dof > 0){
    Mtrans->rptr[1][4] = p[0];
    Mtrans->rptr[2][4] = p[1];
    Mtrans->rptr[3][4] = p[2];
  }

  if(dof > 3){
    angles[0] = p[3]*(M_PI/180);
    angles[1] = p[4]*(M_PI/180);
    angles[2] = p[5]*(M_PI/180);
    Mrot = MRIangles2RotMat(angles);
  } else Mrot = MatrixIdentity(4,NULL);

  Mscale = MatrixIdentity(4,NULL);
  if(dof > 6){
    Mscale->rptr[1][1] = p[6];
    Mscale->rptr[2][2] = p[7];
    Mscale->rptr[3][3] = p[8];
  }

  Mshear = MatrixIdentity(4,NULL);
  if(dof > 9){
    Mshear->rptr[1][2] = p[9];
    Mshear->rptr[1][3] = p[10];
    Mshear->rptr[2][3] = p[11];
  }

  // R = Mshear*Mscale*Mtrans*Mrot*R0
  R = MatrixMultiply(Mrot,R0,R);
  R = MatrixMultiply(Mtrans,R,R);
  R = MatrixMultiply(Mscale,R,R);
  R = MatrixMultiply(Mshear,R,R);
  MatrixFree(&Mrot);
  MatrixFree(&Mtrans);
  MatrixFree(&Mscale);
  MatrixFree(&Mshear);

  //printf("Trans: %g %g %g\n",p[0],p[1],p[2]);
  //printf("Rot:   %g %g %g\n",p[3],p[4],p[5]);
  //printf("Scale: %g %g %g\n",p[6],p[7],p[8]);

  if(UseLH){
    vlhwm  = MRIvol2surfVSM(mov,R,lhwm,  vsm, interpcode, NULL, 0, 0, nsubsamp, vlhwm);
    vlhctx = MRIvol2surfVSM(mov,R,lhctx, vsm, interpcode, NULL, 0, 0, nsubsamp, vlhctx);
    if(lhcost == NULL) lhcost = MRIclone(vlhctx,NULL);
    if(lhcon == NULL)  lhcon  = MRIclone(vlhctx,NULL);
  }
  if(UseRH){
    vrhwm  = MRIvol2surfVSM(mov,R,rhwm,  vsm, interpcode, NULL, 0, 0, nsubsamp, vrhwm);
    vrhctx = MRIvol2surfVSM(mov,R,rhctx, vsm, interpcode, NULL, 0, 0, nsubsamp, vrhctx);
    if(rhcost == NULL) rhcost = MRIclone(vrhctx,NULL);
    if(rhcon == NULL)  rhcon  = MRIclone(vrhctx,NULL);
  }

  for(n = 0; n < 8; n++) costs[n] = 0;

  dsum = 0.0;
  dsum2 = 0.0;
  csum = 0.0;
  csum2 = 0.0;
  nhits = 0;

  //fp = fopen("tmp.dat","w");
  if(UseLH){
    for(n = 0; n < lhwm->nvertices; n += nsubsamp){
      if (lhwm->vertices[n].ripflag != 0) continue ;
      if(lhcostfile || lhcost0file) MRIsetVoxVal(lhcost,n,0,0,0,0.0);
      if(lhconfile)  MRIsetVoxVal(lhcon,n,0,0,0,0.0);
      if(lhCortexLabel && MRIgetVoxVal(lhCortexLabel,n,0,0,0) < 0.5) continue;
      if(UseMask && MRIgetVoxVal(lhsegmask,n,0,0,0) < 0.5) continue;
      if(UseLabel && MRIgetVoxVal(lhlabel,n,0,0,0) < 0.5) continue;
      vwm = MRIgetVoxVal(vlhwm,n,0,0,0);
      if(vwm == 0.0) continue;
      vctx = MRIgetVoxVal(vlhctx,n,0,0,0);
      if(vctx == 0.0) continue;
      nhits++;
      costs[1] += vwm;
      costs[2] += (vwm*vwm);
      costs[4] += vctx;
      costs[5] += (vctx*vctx);
      c = VertexCost(vctx, vwm, PenaltySlope, PenaltyCenter, PenaltySign, &d);
      if(TargConLH){
	val = MRIgetVoxVal(TargConLH,n,0,0,0);
	c = (d-val)*(d-val);
      }
      dsum += d;
      dsum2 += (d*d);
      csum += c;
      csum2 += (c*c);
      if(lhcostfile || lhcost0file) MRIsetVoxVal(lhcost,n,0,0,0,c);
      if(lhconfile)  MRIsetVoxVal(lhcon,n,0,0,0,d);
      //fprintf(fp,"%6d %lf %lf %lf %lf %lf %lf %lf\n",nhits,vwm,vctx,d,dsum,a,c,csum);
    }
  }
  if(UseRH){
    for(n = 0; n < rhwm->nvertices; n += nsubsamp){
      if (rhwm->vertices[n].ripflag != 0)
        continue ;
      if(rhcostfile || rhcost0file) MRIsetVoxVal(rhcost,n,0,0,0,0.0);
      if(rhconfile)  MRIsetVoxVal(rhcon,n,0,0,0,0.0);
      if(rhCortexLabel && MRIgetVoxVal(rhCortexLabel,n,0,0,0) < 0.5) continue;
      if(UseMask && MRIgetVoxVal(rhsegmask,n,0,0,0) < 0.5) continue;
      if(UseLabel && MRIgetVoxVal(rhlabel,n,0,0,0) < 0.5) continue;
      vwm = MRIgetVoxVal(vrhwm,n,0,0,0);
      if(vwm == 0.0) continue;
      vctx = MRIgetVoxVal(vrhctx,n,0,0,0);
      if(vctx == 0.0) continue;
      nhits++;
      costs[1] += vwm;
      costs[2] += (vwm*vwm);
      costs[4] += vctx;
      costs[5] += (vctx*vctx);
      c = VertexCost(vctx, vwm, PenaltySlope, PenaltyCenter, PenaltySign, &d);
      if(TargConRH){
	val = MRIgetVoxVal(TargConRH,n,0,0,0);
	c = (d-val)*(d-val);
      }
      dsum += d;
      dsum2 += (d*d);
      csum += c;
      csum2 += (c*c);
      if(rhcostfile || rhcost0file) MRIsetVoxVal(rhcost,n,0,0,0,c);
      if(rhconfile)  MRIsetVoxVal(rhcon,n,0,0,0,d);
      //fprintf(fp,"%6d %lf %lf %lf %lf %lf %lf %lf\n",nhits,vwm,vctx,d,dsum,a,c,csum);
    }
    //fclose(fp);
  }

  dmean = dsum/nhits;
  dstd  = sum2stddev(dsum,dsum2,nhits);
  //printf("dsum=%g, dsum2=%g, dstd=%g\n",dsum,dsum2,dstd);
  cmean = csum/nhits;
  cstd  = sum2stddev(csum,csum2,nhits);
  //printf("csum=%g, csum2=%g, cstd=%g\n",csum,csum2,cstd);

  costs[0] = nhits;
  costs[2] = sum2stddev(costs[1],costs[2],nhits); // wm std
  costs[1] = costs[1]/nhits; // wm mean
  costs[3] = dstd; // std in percent contrast
  costs[5] = sum2stddev(costs[4],costs[5],nhits); // ctx std
  costs[4] = costs[4]/nhits; // ctx mean
  costs[6] = dmean; // percent contrast
  costs[7] = cmean;
  if(nhits == 0) costs[7] = 10.0; 

  if(UseLH){
    //MRIfree(&vlhwm);
    //MRIfree(&vlhctx);
    if(lhcost0file){
      printf("Writing lh init cost to %s\n",lhcost0file);
      MRIwrite(lhcost,lhcost0file);
    }
    if(lhcostfile){
      printf("Writing lh cost to %s\n",lhcostfile);
      MRIwrite(lhcost,lhcostfile);
    }
    if(lhconfile){
      printf("Writing lh con to %s\n",lhconfile);
      MRIwrite(lhcon,lhconfile);
    }
  }

  if(UseRH){
    //MRIfree(&vrhwm);
    //MRIfree(&vrhctx);
    if(rhcost0file){
      printf("Writing rh init cost to %s\n",rhcost0file);
      MRIwrite(rhcost,rhcost0file);
    }
    if(rhcostfile){
      printf("Writing rh cost to %s\n",rhcostfile);
      MRIwrite(rhcost,rhcostfile);
    }
    if(rhconfile){
      printf("Writing rh con to %s\n",rhconfile);
      MRIwrite(rhcon,rhconfile);
    }
  }

  return(costs);
}

/*---------------------------------------------------------*/
int MinPowell(MRI *mov, MRI *notused, MATRIX *R, double *params,
	      int dof, double ftol, double linmintol, int nmaxiters,
	      char *costfile, double *costs, int *niters)
{
  MATRIX *R0;
  float *pPowel, **xi, fret;
  int    r, c, n;

  printf("Init Powel Params dof = %d\n",dof);
  pPowel = vector(1, dof) ;
  for(n=0; n < dof; n++) {
    pPowel[n+1] = params[n];
    printf("%d %g\n",n,params[n]);
  }

  xi = matrix(1, dof, 1, dof) ;
  for (r = 1 ; r <= dof ; r++) {
    for (c = 1 ; c <= dof ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  R0 = MatrixCopy(R,NULL);

  OpenPowell2(pPowel, xi, dof, ftol, linmintol, nmaxiters, 
	      niters, &fret, compute_powell_cost);
  printf("Powell done niters = %d\n",*niters);

  for(n=0; n < dof; n++) params[n] = pPowel[n+1];
  GetSurfCosts(mov, NULL, R0, R, params, dof, costs);

  free_matrix(xi, 1, dof, 1, dof);
  free_vector(pPowel, 1, dof);
  return(NO_ERROR) ;
}

/*-------------------------------------------------------*/
double RelativeSurfCost(MRI *mov, MATRIX *R0)
{
  double params[6], costs[8], cost0, costsum, costavg, rX, rY, rZ, rcost;
  MATRIX *R;
  int n;

  printf("Computing relative cost\n");

  // Get costs at R0 (dof=6)
  for(n=0; n<6; n++) params[n] = 0.0;
  R = MatrixIdentity(4,NULL);
  GetSurfCosts(mov, NULL, R0, R, params, 6, costs);
  cost0 = costs[7];

  // Now rotate in each dimension by 25 degrees to get 
  // an idea of what the cost is away from R0
  costsum = 0.0;
  n = 0;
  for(rX = -25; rX <= 25; rX += 50){
    for(rY = -25; rY <= 25; rY += 50){
      for(rZ = -25; rZ <= 25; rZ += 50){
	params[3] = rX;
	params[4] = rY;
	params[5] = rZ;
	GetSurfCosts(mov, NULL, R0, R, params, 6, costs);
	costsum += costs[7];
	printf("%2d  %5.1f %5.1f %5.1f   %8.6f\n",n,rX,rY,rZ,costs[7]);
	n++;
      }
    }
  }

  costavg = costsum/n;
  rcost = cost0/costavg;

  printf("REL: %2d  %8.6f %11.6f  %8.6f rel = %g \n",n,cost0,costsum,costavg,rcost);

  MatrixFree(&R);
  return(rcost);
}
