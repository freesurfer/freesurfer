/**
 * @file  mri_segreg.c
 * @brief compute/optimize cost function of segmentation-based registration
 *
 */
/*
 * Original Author: Greg Grev
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2008/08/20 16:11:13 $
 *    $Revision: 1.62 $
 *
 * Copyright (C) 2007,
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


/*
  BEGINUSAGE --------------------------------------------------------------

  mri_segreg

  --init-reg regfile
  --regheader subject
  --mov fvol
  --out-reg outreg : reg at lowest cost

  --cost costfile
  --sum sumfile : default is outreg.sum
  --o out : save final output

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
  --cf cfile  : save cost function values (pct,cost)

  --mask : mask out regions edge and B0 regions
  --no-mask : do not mask out (default)

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

  --aseg : use aseg instead of segreg.mgz

  --nmax nmax   : max number of powell iterations (def 36)
  --tol   tol   : powell inter-iteration tolerance on cost. 
     This is the fraction of the cost that the difference in 
     successive costs must drop below to stop the optimization. 

  --tol1d tol1d : tolerance on powell 1d minimizations

  --1dmin : use brute force 1D minimizations instead of powell
  --n1dmin n1dmin : number of 1d minimization (default = 3)

  --mincost MinCostFile
  --rms     RMSDiffFile : saves Tx Ty Tz Ax Ay Az RMSDiff MinCost 
              WMMean CtxMean PctContrast C0 Slope NSubSamp UseMask
  --surf-cost basename : saves as basename.?h.mgh
  --init-surf-cost basename0 : saves init cost as basename0.?h.mgh
  --surf-cost-diff diffbase : saves final-init cost as diffbase.?h.mgh
  --surf-con basename : saves final contrast as basename.?h.mgh

  --mksegreg subject : create segreg.mgz and exit

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

#ifdef X
#undef X
#endif

int MRIsegReg(char *subject);
int MRISsegReg(char *subject, int ForceReRun);

double *GetCosts(MRI *mov, MRI *seg,
                 MATRIX *R0, MATRIX *R,
                 double *p, double *costs);
double *GetSurfCosts(MRI *mov, MRI *notused, MATRIX *R0, MATRIX *R,
		     double *p, double *costs);
int Min1D(MRI *mov, MRI *seg,
          MATRIX *R, double *p, char *costfile, double *costs);
float compute_powell_cost(float *p) ;
int MinPowell(MRI *mov, MRI *seg, MATRIX *R, double *params,
	      double ftol, double linmintol, int nmaxiters,
	      char *costfile, double *costs, int *niters);
MRI *mov_powell = NULL;
MRI *seg_powell = NULL;
MATRIX *R0_powell = NULL;
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
static int istringnmatch(char *str1, char *str2, int n);
double VertexCost(double vctx, double vwm, double slope, 
		  double center, double sign, double *pct);


int main(int argc, char *argv[]) ;

static char vcid[] =
"$Id: mri_segreg.c,v 1.62 2008/08/20 16:11:13 greve Exp $";
char *Progname = NULL;

int debug = 0, gdiagno = -1;

char *movvolfile=NULL;
char *regfile=NULL;
char *outregfile=NULL;
char *sumfile=NULL;

char *interpmethod = "nearest";
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
MRI *segreg, *segreg0;
MRI *noise=NULL;
MRI *mritmp;
MRI *inorm=NULL;
char *inormfile=NULL;
char *outfile=NULL;

int frame = 0;

int SynthSeed = -1;
int AddNoise = 0;
double NoiseStd;

int UseASeg = 0;
int DoCrop = 1;
int DoProfile = 0;
int n1dmin = 3;

int DoPowell = 1;
int nMaxItersPowell = 36;
double TolPowell = 1e-8;
double LinMinTolPowell = 1e-8;

int MkSegReg = 0;

#define NMAX 100
int ntx=0, nty=0, ntz=0, nax=0, nay=0, naz=0;
double txlist[NMAX],tylist[NMAX],tzlist[NMAX];
double axlist[NMAX],aylist[NMAX],azlist[NMAX];

int UseSurf = 1;
MRIS *lhwm, *rhwm, *lhctx, *rhctx;
MRI *lhsegmask, *rhsegmask;
int UseMask = 0;
char *cmdline2;
struct utsname uts;

int PenaltySign  = -1;
double PenaltySlope = .5;
double PenaltyCenter = 0;
int DoMidFrame = 0;

int Do1DPreOpt = 0;
double PreOptMin, PreOptMax, PreOptDelta, PreOpt, PreOptAtMin;
int PreOptDim = 2;
char *PreOptFile = NULL;
int PreOptOnly = 0;

char *surfcostbase=NULL, *lhcostfile=NULL, *rhcostfile=NULL;
char *surfconbase=NULL, *lhconfile=NULL, *rhconfile=NULL;
char *surfcost0base=NULL, *lhcost0file=NULL, *rhcost0file=NULL;
char *surfcostdiffbase=NULL;

int UseLH = 1;
int UseRH = 1;

MATRIX *MrotPre=NULL,*MtransPre=NULL;
double TransRandMax = 0;
double RotRandMax = 0;
char *MinCostFile=NULL;

int DoRMSDiff = 1; // this is fast, so ok to do automatically
char *RMSDiffFile = NULL;

int DoSmooth = 0;
double fwhm = 0, gstd = 0;
int nsubsamp = 1;

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
int nCostEvaluations=0;

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  char cmdline[CMD_LINE_LEN] ;
  double costs[8], mincost, p[6], pmin[6];
  double tx, ty, tz, ax, ay, az;
  int nth,nthtx, nthty, nthtz, nthax, nthay, nthaz, ntot, n, err,vno;
  MATRIX *Ttemp=NULL, *invTtemp=NULL, *Stemp=NULL, *invStemp=NULL;
  MATRIX *R=NULL, *Rcrop=NULL, *R00=NULL, *Rdiff=NULL;
  MATRIX *Rmin=NULL, *Scrop=NULL, *invScrop=NULL, *invTcrop=NULL, *Tcrop=NULL;
  struct timeb  mytimer;
  double secCostTime;
  MRI_REGION box;
  FILE *fp, *fpMinCost, *fpRMSDiff, *fpPreOpt=NULL;
  double rmsDiffSum, rmsDiffMean=0, rmsDiffMax=0, d;
  VERTEX *v;
  int nsubsampsave = 0;
  LABEL *label;
  int PrintT1Warning = 0, PrintT2Warning = 0;

  make_cmd_version_string
    (argc, argv,
     "$Id: mri_segreg.c,v 1.62 2008/08/20 16:11:13 greve Exp $",
     "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mri_segreg.c,v 1.62 2008/08/20 16:11:13 greve Exp $",
     "$Name:  $");
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

  printf("Loading mov\n");
  mov = MRIread(movvolfile);
  if (mov == NULL) exit(1);

  if(DoAbs){
    printf("Computing abs\n");
    MRIabs(mov,mov);
  }

  if(regheader) {
    printf("Computing registration based on scanner-to-scanner\n");
    sprintf(tmpstr,"%s/%s/mri/orig.mgz",SUBJECTS_DIR,subject);
    segreg = MRIread(tmpstr); // Just need a template
    if(segreg == NULL) exit(1);
    R0 = MRItkRegMtx(segreg,mov,NULL);
  }

  R00 = MatrixCopy(R0,NULL);

  if(MrotPre || MtransPre){
    printf("Applying Pre Transform to input reg\n");
    if(MrotPre){
      printf("Rot:\n");
      MatrixPrint(stdout,MrotPre);
      fprintf(fp,"Rot:\n");
      MatrixPrint(fp,MrotPre);
      R0 = MatrixMultiply(MrotPre,R0,R0);
    }
    if(MtransPre){
      printf("Trans:\n");
      MatrixPrint(stdout,MtransPre);
      fprintf(fp,"Trans:\n");
      MatrixPrint(fp,MtransPre);
      R0 = MatrixMultiply(MtransPre,R0,R0);
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

  if(inormfile){
    printf("Loading inorm\n");
    inorm = MRIread(inormfile);
    if(inorm == NULL) exit(1);
  }

  if(UseSurf){
    MRISsegReg(subject,0);
    sprintf(tmpstr,"%s/%s/mri/orig.mgz",SUBJECTS_DIR,subject);
    segreg = MRIread(tmpstr); // Just need a template
    if(segreg == NULL) exit(1);
    mritmp = MRIchangeType(segreg,MRI_FLOAT,0,0,0);
    MRIfree(&segreg);
    segreg = mritmp;
    if(UseLH){
      // Load the LH cortex label
      label = LabelRead(subject, "lh.cortex.label");
      if(label == NULL) exit(1);
      lhCortexLabel = MRISlabel2Mask(lhwm, label, NULL);
      LabelFree(&label);
    }
    if(UseRH){
      // Load the RH cortex label
      label = LabelRead(subject, "rh.cortex.label");
      if(label == NULL) exit(1);
      rhCortexLabel = MRISlabel2Mask(rhwm, label, NULL);
      LabelFree(&label);
    }

  } else {
  if(!UseASeg){
    if(!MkSegReg){
      printf("Loading segreg\n");
      sprintf(tmpstr,"%s/%s/mri/segreg.mgz",SUBJECTS_DIR,subject);
      if(!fio_FileExistsReadable(tmpstr)){
	printf("\n");
	printf("WARNING: cannot find %s\n",tmpstr);
	printf("So, I'm going to create it on-the-fly.\n");
	MkSegReg = 1;
      }
    }
    if(MkSegReg){
      printf("Creating segreg\n");
      err = MRIsegReg(subject);
      if(err) exit(err);
    }
    segreg = MRIread(tmpstr);
    if(segreg == NULL) exit(1);
    free(fspec);
  } else {
    printf("Loading aseg\n");
    sprintf(tmpstr,"%s/%s/mri/aseg",SUBJECTS_DIR,subject);
    fspec = IDnameFromStem(tmpstr);
    segreg = MRIread(fspec);
    if(segreg == NULL) exit(1);
    free(fspec);
  }
  }
  segreg0 = segreg;


  // Cropping reduces the size of the target volume down to the
  // voxels that really matter. This can greatly increase the speed
  // (factor of 2-3). Requires that the registration matrix be
  // recomputed for the smaller volume.
  if(DoCrop){
    printf("Cropping\n");

    // Prepare to adjust the input reg matrix
    Ttemp    = MRIxfmCRS2XYZtkreg(segreg); // Vox-to-tkRAS Matrices
    invTtemp = MatrixInverse(Ttemp,NULL);
    Stemp    = MRIxfmCRS2XYZ(segreg,0); // Vox-to-ScannerRAS Matrices
    invStemp = MatrixInverse(Stemp,NULL);

    err = MRIboundingBox(segreg, 0.5, &box);
    if(err) exit(1);
    printf("BBbox start: %d %d %d, delta = %d %d %d\n",
           box.x,box.y,box.z,box.dx,box.dy,box.dz);
    // Now crop
    segreg  = MRIcrop(segreg0,
                      box.x, box.y, box.z,
                      box.x+box.dx, box.y+box.dy, box.z+box.dz);
    if(segreg == NULL) exit(1);
    if(inorm){
      // Crop inorm too
      mritmp = MRIcrop(inorm,
                       box.x, box.y, box.z,
                       box.x+box.dx, box.y+box.dy, box.z+box.dz);
      if(mritmp == NULL) exit(1);
      MRIfree(&inorm);
      inorm = mritmp;
    }

    Tcrop  = MRIxfmCRS2XYZtkreg(segreg); // Vox-to-tkRAS Matrices
    invTcrop = MatrixInverse(Tcrop,NULL);
    Scrop  = MRIxfmCRS2XYZ(segreg,0); // Vox-to-ScannerRAS Matrices
    invScrop = MatrixInverse(Scrop,NULL);

    // Now adjust input reg
    // Rc = R*T*inv(S)*Sc*inv(Tc)
    R0 = MatrixMultiply(R0,Ttemp,R0);
    R0 = MatrixMultiply(R0,invStemp,R0);
    R0 = MatrixMultiply(R0,Scrop,R0);
    R0 = MatrixMultiply(R0,invTcrop,R0);

    printf("Input reg after cropping\n");
    MatrixPrint(stdout,R0);
    printf("\n");

    if(0){
      MatrixFree(&Ttemp);
      MatrixFree(&Stemp);
      MatrixFree(&invStemp);
      MatrixFree(&Tcrop);
      MatrixFree(&invTcrop);
      MatrixFree(&Scrop);
    }

    //MRIwrite(segreg,"segregcrop.mgz");
    //regio_write_register("crop.reg",subject,mov->xsize,
    //     mov->zsize,intensity,R0,FLT2INT_ROUND);
    //exit(1);
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

  // Allocate the output
  out = MRIallocSequence(segreg->width,
                         segreg->height,
                         segreg->depth,
                         MRI_FLOAT, 1);
  MRIcopyHeader(segreg,out);

  // Compute cost at initial
  for(nth=0; nth < 6; nth++) p[nth] = 0.0;
  if(!UseSurf)  GetCosts(mov,     segreg, R0, R, p, costs);
  else{
    if(surfcost0base){
      if(UseLH){
	sprintf(tmpstr,"%s.lh.mgh",surfcost0base);
	lhcost0file = strcpyalloc(tmpstr);
      }
      if(UseRH){
	sprintf(tmpstr,"%s.rh.mgh",surfcost0base);
	rhcost0file = strcpyalloc(tmpstr);
      }
    }
    GetSurfCosts(mov, segreg, R0, R, p, costs);
    if(UseLH) lhcost0 = MRIcopy(lhcost,NULL);
    if(UseRH) rhcost0 = MRIcopy(rhcost,NULL);
    //if(lhcost0file)  MRIwrite(lhcost0,lhcost0file);
    //if(rhcost0file)  MRIwrite(rhcost0,rhcost0file);
    free(lhcost0file);lhcost0file=NULL;
    free(rhcost0file);rhcost0file=NULL;
  }

  if(DoProfile){
    printf("Profiling over 100 iterations\n");
    TimerStart(&mytimer) ;
    for(n=0; n < 100; n++) GetSurfCosts(mov, segreg, R0, R, p, costs);
    secCostTime = TimerStop(&mytimer)/1000.0;
    printf("ttot = %g\n",secCostTime);
    printf("tper = %g\n",secCostTime/100);
    exit(0);
  }

  fprintf(fp,"Initial costs ----------------\n");
  fprintf(fp,"Number of surface hits %d\n",(int)costs[0]);  
  fprintf(fp,"WM  Intensity0 %10.4lf +/- %8.4lf\n",costs[1],costs[2]); 
  fprintf(fp,"Ctx Intensity0 %10.4lf +/- %8.4lf\n",costs[4],costs[5]); 
  fprintf(fp,"Pct Contrast0  %10.4lf +/- %8.4lf\n",costs[6],costs[3]); 
  fprintf(fp,"Cost %8.4lf\n",costs[7]); 
  fflush(fp);

  printf("Initial costs ----------------\n");
  printf("Number of surface hits %d\n",(int)costs[0]);  
  printf("WM  Intensity %10.4lf +/- %8.4lf\n",costs[1],costs[2]); 
  printf("Ctx Intensity %10.4lf +/- %8.4lf\n",costs[4],costs[5]); 
  printf("Pct Contrast  %10.4lf +/- %8.4lf\n",costs[6],costs[3]); 
  printf("Cost %8.4lf\n",costs[7]); 

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
    for(nth=0; nth < 6; nth++) p[nth] = 0.0;
    nth = 0;
    if(PreOptFile) fpPreOpt = fopen(PreOptFile,"w");
    for(PreOpt = PreOptMin; PreOpt <= PreOptMax; PreOpt += PreOptDelta){
      p[PreOptDim] = PreOpt;
      GetSurfCosts(mov, segreg, R0, R, p, costs);
      if(costs[7] < mincost) {
	mincost = costs[7];
	PreOptAtMin = PreOpt;
      }
      printf("%8.4lf %8.4lf %8.4lf \n",PreOpt,costs[7],mincost);
      fprintf(fp,"%8.4lf %8.4lf %8.4lf \n",PreOpt,costs[7],mincost);
      if(PreOptFile) fprintf(fpPreOpt,"%8.4lf %8.4lf \n",PreOpt,costs[7]);
      nth ++;
    }
    if(PreOptFile) fclose(fpPreOpt);
    if(PreOptOnly) {
      printf("PreOptOnly specified, so exiting now\n");
      exit(0);
    }
    p[PreOptDim] = PreOptAtMin; // phase encode direction
    GetSurfCosts(mov, segreg, R0, R, p, costs);
    MatrixCopy(R0,R);
    printf("\n");
    fprintf(fp,"\n");
  }

  if(BruteForce){
    n = pow((PreOptMax-PreOptMin)/PreOptDelta+1.0,6.0);
    printf("\n");
    printf("------------------------------------\n");
    nsubsampsave = nsubsamp;
    nsubsamp = 50;
    printf("Bruce force preopt %g %g %g, n = %d\n",PreOptMin,PreOptMax,PreOptDelta,n);
    fprintf(fp,"\n");
    fprintf(fp,"Bruce force preopt %g %g %g, n = %d\n",PreOptMin,PreOptMax,PreOptDelta,n);
    fflush(stdout); fflush(fp);
    mincost = 10e10;
    PreOptAtMin = 0;
    for(n=0; n < 6; n++) p[n] = 0.0;
    TimerStart(&mytimer) ;
    nth = 0;
    if(PreOptFile) fpPreOpt = fopen(PreOptFile,"w");
    for(tx = PreOptMin; tx <= PreOptMax; tx += PreOptDelta){
      for(ty = PreOptMin; ty <= PreOptMax; ty += PreOptDelta){
	for(tz = PreOptMin; tz <= PreOptMax; tz += PreOptDelta){
	  for(ax = PreOptMin; ax <= PreOptMax; ax += PreOptDelta){
	    for(ay = PreOptMin; ay <= PreOptMax; ay += PreOptDelta){
	      for(az = PreOptMin; az <= PreOptMax; az += PreOptDelta){
		p[0] = tx;
		p[1] = ty;
		p[2] = tz;
		p[3] = ax;
		p[4] = ay;
		p[5] = az;
		GetSurfCosts(mov, segreg, R0, R, p, costs);
		if(costs[7] < mincost) {
		  mincost = costs[7];
		  for(n=0; n < 6; n++) pmin[n] = p[n];
		  secCostTime = TimerStop(&mytimer)/1000.0 ;
		  printf("%6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf %8.4lf %4.1f\n",
			 nth,tx,ty,tz,ax,ay,az,costs[7],mincost,secCostTime/60);
		  fprintf(fp,"%6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf %8.4lf %4.1f\n",
			  nth,tx,ty,tz,ax,ay,az,costs[7],mincost,secCostTime/60);
		  fflush(stdout); fflush(fp);

		} else {
		  if(nth == 0 || nth%1000 == 0){
		    secCostTime = TimerStop(&mytimer)/1000.0 ;
		    printf("%6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf %8.4lf %4.1f\n",
			   nth,tx,ty,tz,ax,ay,az,costs[7],mincost,secCostTime/60);
		    fprintf(fp,"%6d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf %8.4lf %4.1f\n",
			    nth,tx,ty,tz,ax,ay,az,costs[7],mincost,secCostTime/60);
		    fflush(stdout); fflush(fp);
		  }
		}
		if(PreOptFile) 
		  fprintf(fpPreOpt,"%8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf    %8.4lf\n",
			  tx,ty,tz,ax,ay,az,costs[7]);
		nth ++;
	      }
	    }
	  }
	}
      }
    }
    for(n=0; n < 6; n++) p[n] = pmin[n];
    GetSurfCosts(mov, segreg, R0, R, p, costs);
    MatrixCopy(R0,R);
    printf("Brute Force --------------------------\n");
    printf("Min cost was %lf\n",mincost);
    printf("Number of iterations %5d\n",nth);
    secCostTime = TimerStop(&mytimer)/1000.0 ;
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

  TimerStart(&mytimer) ;
  if(DoPowell) {
    printf("Staring Powell Minimization\n");
    MinPowell(mov, segreg, R, p, TolPowell, LinMinTolPowell,
	      nMaxItersPowell,SegRegCostFile, costs, &nth);
  }
  else{
    printf("Staring 1D Minimization\n");
    // 1D minimization
    nth = 0;
    for(n=0; n < n1dmin; n++){
      printf("n = %d --------------------------------------\n",n);
      R = MatrixCopy(R0,NULL);
      nth += Min1D(mov, segreg, R, p, SegRegCostFile, costs);
      printf("\n");
    }
  }
  secCostTime = TimerStop(&mytimer)/1000.0 ;

  // Recompute at optimal. This forces MRI *out to be the output at best reg
  if(!UseSurf) GetCosts(mov, segreg, R0, R, p, costs);
  else{
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
    GetSurfCosts(mov, segreg, R0, R, p, costs);
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
  printf("Parameters at optimum (transmm, rotdeg)\n");
  printf("%8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf \n",
	 p[0],p[1],p[2],p[3],p[4],p[5]);
  printf("Costs at optimum\n");
  printf("%7d %10.4lf %8.4lf ",
	 (int)costs[0],costs[1],costs[2]); // WM  n mean std
  printf("%10.4lf %10.4lf %8.4lf ",
	 costs[3],costs[4],costs[5]); // CTX n mean std
  printf("%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
  printf("\n");
  
  printf("Reg at min cost was \n");
  MatrixPrint(stdout,R);
  printf("\n");
  
  fprintf(fp,"Number of iterations %5d\n",nth);
  fprintf(fp,"Min cost was %lf\n",costs[7]);
  fprintf(fp,"Number of FunctionCalls %5d\n",nCostEvaluations);
  fprintf(fp,"OptimizationTime %lf sec\n",secCostTime);
  fprintf(fp,"Parameters at optimum (transmm, rotdeg)\n");
  fprintf(fp,"%8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf \n",
	 p[0],p[1],p[2],p[3],p[4],p[5]);
  fprintf(fp,"Number of surface hits %d\n",(int)costs[0]);  
  fprintf(fp,"WM  Intensity %10.4lf +/- %8.4lf\n",costs[1],costs[2]); 
  fprintf(fp,"Ctx Intensity %10.4lf +/- %8.4lf\n",costs[4],costs[5]); 
  fprintf(fp,"Pct Contrast  %10.4lf +/- %8.4lf\n",costs[6],costs[3]); 
  fprintf(fp,"Cost at optimum %8.4lf\n",costs[7]); 
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
  
  if(DoCrop){
    printf("Uncropping reg\n");
    // R = Rc*Tc*inv(Sc)*S*inv(T)
    Rcrop = MatrixCopy(R,NULL);
    R = MatrixMultiply(Rcrop,Tcrop,R);
    R = MatrixMultiply(R,invScrop,R);
    R = MatrixMultiply(R,Stemp,R);
    R = MatrixMultiply(R,invTtemp,R);
    printf("Uncropped reg at min cost was \n");
    MatrixPrint(stdout,R);
    printf("\n");
  }
  
  if(outregfile){
    printf("Writing optimal reg to %s \n",outregfile);
    fflush(stdout);
    regio_write_register(outregfile,subject,mov->xsize,
			 mov->zsize,intensity,R,FLT2INT_ROUND);
  }
  
  if(outfile) {
    // This changes values in segreg0
    printf("Writing output volume to %s \n",outfile);
    fflush(stdout);
    MRIsetValues(segreg0,0.0);
    MRIvol2VolTkReg(mov, segreg0, R, SAMPLE_TRILINEAR, sinchw);
    MRIwrite(segreg0,outfile);
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
    rmsDiffMax = 0;
    if(lhwm){
      printf("Computing change in lh position\n");
      MRISmatrixMultiply(lhwm,Rdiff);
      for(vno = 0; vno < lhwm->nvertices; vno++){
	v = &(lhwm->vertices[vno]);
	d = sqrt(v->x*v->x + v->y*v->y + v->z*v->z );
	rmsDiffSum += d;
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
	if(rmsDiffMax < d) rmsDiffMax = d;
      }
      if(lhwm)
	rmsDiffMean = rmsDiffSum/(lhwm->nvertices + rhwm->nvertices);
      else
	rmsDiffMean = rmsDiffSum/(rhwm->nvertices);
    }
    printf("Surface RMS Diff (mm) %lf %lf\n",rmsDiffMean,rmsDiffMax);
    fprintf(fp,"Surface RMS Diff (mm) %lf %lf\n",rmsDiffMean,rmsDiffMax);
    if(RMSDiffFile){
      fpRMSDiff = fopen(RMSDiffFile,"w");
      //rmsDiffMean rmsDiffMax MinCost Nevals Tx Ty Tz Rx Ry Rz  WMMean CtxMean PctContrast 
      fprintf(fpRMSDiff,
	      "%12.9lf %12.9lf %12.9lf %4d  %7.3lf %7.3lf %7.3lf %6.3lf %6.3lf %6.3lf  %lf %lf %lf \n",
	      rmsDiffMean,rmsDiffMax,costs[7],nCostEvaluations,
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
  /*- end main ----------------------------------------------------*/

  ntot = ntx*nty*ntz*nax*nay*naz;
  printf("Staring loop over %d\n",ntot);
  nth = 0;
  mincost = 10e10;
  for(nthtx = 0; nthtx < ntx; nthtx++){
    tx = txlist[nthtx];
    for(nthty = 0; nthty < nty; nthty++){
      ty = tylist[nthty];
      for(nthtz = 0; nthtz < ntz; nthtz++){
        tz = tzlist[nthtz];
        for(nthax = 0; nthax < nax; nthax++){
          ax = axlist[nthax];
          for(nthay = 0; nthay < nay; nthay++){
            ay = aylist[nthay];
            for(nthaz = 0; nthaz < naz; nthaz++){
              az = azlist[nthaz];
              nth ++;

              p[0] = tx;
              p[1] = ty;
              p[2] = tz;
              p[3] = ax;
              p[4] = ay;
              p[5] = az;
              TimerStart(&mytimer) ;
              GetCosts(mov, segreg, R0, R, p, costs);
              secCostTime = TimerStop(&mytimer)/1000.0 ;

              // write costs to file
              fp = fopen(SegRegCostFile,"a");
              fprintf(fp,"%7.3lf %7.3lf %7.3lf ",tx,ty,tz);
              fprintf(fp,"%6.3lf %6.3lf %6.3lf ",ax,ay,az);
              fprintf(fp,"%7d %10.4lf %8.4lf ",
                      (int)costs[0],costs[1],costs[2]); // WM  n mean std
              fprintf(fp,"%7d %10.4lf %8.4lf ",
                      (int)costs[3],costs[4],costs[5]); // CTX n mean std
              fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
              fprintf(fp,"\n");
              fclose(fp);

              fp = stdout;
              fprintf(fp,"%5d ",nth);
              fprintf(fp,"%7.3lf %7.3lf %7.3lf ",tx,ty,tz);
              fprintf(fp,"%6.3lf %6.3lf %6.3lf ",ax,ay,az);
              fprintf(fp,"%7d %10.4lf %8.4lf ",
                      (int)costs[0],costs[1],costs[2]); // WM  n mean std
              fprintf(fp,"%7d %10.4lf %8.4lf ",
                      (int)costs[3],costs[4],costs[5]); // CTX n mean std
              fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
              if(DoProfile) fprintf(fp,"%4.2lf ",secCostTime);
              printf("\n");
              fflush(stdout);

              if(mincost > costs[7]){
                mincost = costs[7];
                Rmin = MatrixCopy(R,Rmin);
                if(outregfile){
                  regio_write_register(outregfile,subject,mov->xsize,
                                       mov->zsize,intensity,Rmin,FLT2INT_ROUND);
                }
              }

              // clean up
            }
          }
        }
      }
    }
  }

  printf("min cost was %lf\n",mincost);
  printf("Reg at min cost was \n");
  MatrixPrint(stdout,Rmin);
  printf("\n");

  if(outregfile){
    printf("Writing optimal reg to %s \n",outregfile);
    regio_write_register(outregfile,subject,mov->xsize,
                         mov->zsize,intensity,Rmin,FLT2INT_ROUND);
  }

  printf("\n");
  printf("mri_segreg done\n");

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
  double angles[3],xyztrans[3];
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
    else if (!strcasecmp(option, "--aseg"))     UseASeg = 1;
    else if (!strcasecmp(option, "--no-crop"))  DoCrop = 0;
    else if (!strcasecmp(option, "--crop"))     DoCrop = 1;
    else if (!strcasecmp(option, "--profile"))  DoProfile = 1;
    else if (!strcasecmp(option, "--powell"))   DoPowell = 1;
    else if (!strcasecmp(option, "--no-powell")) DoPowell = 0;
    else if (!strcasecmp(option, "--abs"))       DoAbs = 1;
    else if (!strcasecmp(option, "--no-abs"))    DoAbs = 0;
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
    else if (!strcasecmp(option, "--lh-only")){ UseLH = 1; UseRH = 0;}
    else if (!strcasecmp(option, "--rh-only")){ UseLH = 0; UseRH = 1;}
    else if (!strcasecmp(option, "--surf")){
      UseSurf = 1;
      DoCrop = 0;
    }
    else if (!strcasecmp(option, "--no-surf"))   UseSurf = 0;
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
    else if (istringnmatch(option, "--surf-cost-diff",0)) {
      if(nargc < 1) argnerr(option,1);
      surfcostdiffbase = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--mov",0)) {
      if (nargc < 1) argnerr(option,1);
      movvolfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--frame",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
    } else if (istringnmatch(option, "--n1dmin",0)) {
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
    else if (istringnmatch(option, "--nmax",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nMaxItersPowell);
      nargsused = 1;
    } else if (istringnmatch(option, "--tol",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&TolPowell);
      nargsused = 1;
    } else if (istringnmatch(option, "--tol1d",0)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&LinMinTolPowell);
      nargsused = 1;
    } else if (istringnmatch(option, "--inorm",0)) {
      if (nargc < 1) argnerr(option,1);
      inormfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--o",0)) {
      if (nargc < 1) argnerr(option,1);
      outfile = pargv[0];
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
    } else if (istringnmatch(option, "--trans",0)) {
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
    else if (istringnmatch(option, "--trans-rand",0)) {
      if(nargc < 1) argnerr(option,1);
      // Translation in mm
      sscanf(pargv[0],"%lf",&TransRandMax);
      MtransPre = MatrixIdentity(4,NULL);
      MtransPre->rptr[1][4] = 2.0*(drand48()-0.5)*TransRandMax;
      MtransPre->rptr[2][4] = 2.0*(drand48()-0.5)*TransRandMax;
      MtransPre->rptr[3][4] = 2.0*(drand48()-0.5)*TransRandMax;
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--out-reg",0)) {
      if (nargc < 1) argnerr(option,1);
      outregfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--sum",0)) {
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
    } else if (istringnmatch(option, "--abs",0)) {
      // no direction of contrast expected
      PenaltySign = 0;
      nargsused = 1;
    } else if (istringnmatch(option, "--fwhm",0)) {
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
    } else if (istringnmatch(option, "--mksegreg",0)) {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      MkSegReg = 1;
      nargsused = 1;
    } else if (istringnmatch(option, "--tx-mmd",0)) {
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
    } else if (istringnmatch(option, "--cost",0)) {
      if (nargc < 1) argnerr(option,1);
      SegRegCostFile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--mincost",0)) {
      if (nargc < 1) argnerr(option,1);
      MinCostFile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--rms",0)) {
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
printf("  --o out : save final output\n");
printf("\n");
printf("  --brute min max delta : brute force in all directions\n");
printf("  --1dpreopt min max delta : brute force in PE direction\n");
printf("\n");
printf("  --fwhm fwhm : smooth input by fwhm mm\n");
printf("  --abs       : compute abs of mov\n");
printf("  --subsamp nsub : only sample every nsub vertices\n");
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
printf("  --c0 offset : cost offset (pct)\n");
printf("  --cf cfile  : save cost function values (pct,cost)\n");
printf("\n");
printf("  --mask : mask out regions edge and B0 regions\n");
printf("  --no-mask : do not mask out (default)\n");
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
printf("  --no-crop: do not crop anat (crops by default)\n");
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
printf("  --rms     RMSDiffFile : saves Tx Ty Tz Ax Ay Az RMSDiff MinCost \n");
printf("              WMMean CtxMean PctContrast C0 Slope NSubSamp UseMask\n");
printf("  --surf-cost basename : saves final cost as basename.?h.mgh\n");
printf("  --init-surf-cost basename0 : saves init cost as basename0.?h.mgh\n");
printf("  --surf-cost-diff diffbase : saves final-init cost as diffbase.?h.mgh\n");
printf("  --surf-con basename : saves final contrast as basename.?h.mgh\n");

printf("\n");
printf("  --mksegreg subject : create segreg.mgz and exit\n");
printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n%s\n\n",vcid);
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) 
{
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL) {
    printf("ERROR: SUBJECTS_DIR undefined.\n");
    exit(1);
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

  if(UseSurf) DoCrop = 0;

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
  if(inormfile) fprintf(fp,"inorm %s\n",inormfile);
  if(outregfile) fprintf(fp,"outregfile %s\n",outregfile);
  if(outfile) fprintf(fp,"outfile %s\n",outfile);
  fprintf(fp,"UseSurf %d\n",UseSurf);
  fprintf(fp,"UseMask %d\n",UseMask);
  fprintf(fp,"UseLH %d\n",UseLH);
  fprintf(fp,"UseRH %d\n",UseRH);
  fprintf(fp,"nsubsamp %d\n",nsubsamp);
  fprintf(fp,"PenaltySign  %d\n",PenaltySign);
  fprintf(fp,"PenaltySlope %lf\n",PenaltySlope);
  fprintf(fp,"PenaltyCenter %lf\n",PenaltyCenter);
  if(DoGMProjFrac) fprintf(fp,"GMProjFrac %lf\n",GMProjFrac);
  if(DoGMProjAbs)  fprintf(fp,"GMProjAbs %lf\n",GMProjAbs);
  if(DoWMProjFrac) fprintf(fp,"WMProjFrac %lf\n",WMProjFrac);
  if(DoWMProjAbs)  fprintf(fp,"WMProjAbs %lf\n",WMProjAbs);
  fprintf(fp,"lhcostfile %s\n",lhcostfile);
  fprintf(fp,"rhcostfile %s\n",rhcostfile);
  fprintf(fp,"interp  %s (%d)\n",interpmethod,interpcode);
  fprintf(fp,"frame  %d\n",frame);
  fprintf(fp,"DoPowell  %d\n",DoPowell);
  fprintf(fp,"TolPowell %lf\n",TolPowell);
  fprintf(fp,"nMaxItersPowell %d\n",nMaxItersPowell);
  fprintf(fp,"n1dmin  %d\n",n1dmin);
  if(interpcode == SAMPLE_SINC) fprintf(fp,"sinc hw  %d\n",sinchw);
  fprintf(fp,"Crop      %d\n",DoCrop);
  fprintf(fp,"Profile   %d\n",DoProfile);
  fprintf(fp,"Gdiag_no  %d\n",Gdiag_no);
  fprintf(fp,"AddNoise  %d (%g)\n",AddNoise,NoiseStd);
  fprintf(fp,"SynthSeed %d\n",SynthSeed);
  fprintf(fp,"TransRandMax %lf\n",TransRandMax);
  fprintf(fp,"RotRandMax %lf\n",RotRandMax);
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
static int istringnmatch(char *str1, char *str2, int n) {
  if (n > 0  && ! strncasecmp(str1,str2,n)) return(1);
  if (n <= 0 && ! strcasecmp(str1,str2)) return(1);
  return(0);
}

/*-------------------------------------------------------*/
double *GetCosts(MRI *mov, MRI *seg, MATRIX *R0, MATRIX *R,
                 double *p, double *costs)
{
  double angles[3];
  MATRIX *Mrot=NULL, *Mtrans=NULL, *vox2vox = NULL;
  MATRIX *Tin, *invTin, *Ttemp;
  extern MRI *out, *inorm;
  extern int interpcode, sinchw;

  if(R==NULL){
    printf("ERROR: GetCosts(): R cannot be NULL\n");
    return(NULL);
  }

  Tin      = MRIxfmCRS2XYZtkreg(mov);
  invTin   = MatrixInverse(Tin,NULL);
  Ttemp    = MRIxfmCRS2XYZtkreg(seg);

  Mtrans = MatrixIdentity(4,NULL);
  Mtrans->rptr[1][4] = p[0];
  Mtrans->rptr[2][4] = p[1];
  Mtrans->rptr[3][4] = p[2];

  angles[0] = p[3]*(M_PI/180);
  angles[1] = p[4]*(M_PI/180);
  angles[2] = p[5]*(M_PI/180);
  Mrot = MRIangles2RotMat(angles);

  // R = Mtrans*Mrot*R0
  R = MatrixMultiply(Mrot,R0,R);
  R = MatrixMultiply(Mtrans,R,R);

  // vox2vox = invTin*R*Ttemp
  vox2vox = MatrixMultiply(invTin,R,vox2vox);
  MatrixMultiply(vox2vox,Ttemp,vox2vox);

  // Zero output
  MRIsetValues(out,0.0);

  // resample
  MRIvol2Vol(mov,out,vox2vox,interpcode,sinchw);

  if(inorm)  MRImultiply(out,inorm,out);

  // compute costs
  costs = SegRegCost(segreg,out,costs);

  MatrixFree(&Mrot);
  MatrixFree(&Mtrans);
  MatrixFree(&vox2vox);
  MatrixFree(&Tin);
  MatrixFree(&invTin);
  MatrixFree(&Ttemp);

  return(costs);
}

/*---------------------------------------------------------------------*/
int Min1D(MRI *mov, MRI *seg, MATRIX *R, double *p,
          char *costfile, double *costs)
{
  extern int UseSurf;
  double q, q0, pp[6], c, copt=0, qopt=0, costsopt[8], qdelta;
  int nthp, nth, n, hit, nthq;
  MATRIX *R0, *Rtmp;
  FILE *fp;

  qdelta = .2;

  if(R==NULL) exit(1);
  if(p==NULL) exit(1);
  if(costs==NULL) exit(1);

  for(nthp = 0; nthp < 6; nthp++) pp[nthp] = p[nthp];
  R0 = MatrixCopy(R,NULL);
  Rtmp = MatrixAlloc(4,4,MATRIX_REAL);

  if(!UseSurf)  GetCosts(    mov, seg, R0, Rtmp, pp, costs);
  else          GetSurfCosts(mov, seg, R0, Rtmp, pp, costs);
  
  copt = costs[7];

  fp = stdout;
  fprintf(fp,"init1d ");
  fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
  fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
  fprintf(fp,"%7d %10.4lf %8.4lf ",
          (int)costs[0],costs[1],costs[2]); // WM  n mean std
  fprintf(fp,"%7d %10.4lf %8.4lf ",
          (int)costs[3],costs[4],costs[5]); // CTX n mean std
  fprintf(fp,"%8.4lf %8.4lf   %8.4lf ",costs[6],costs[7],copt); // t, cost=1/t
  printf("\n");
  printf("\n");
  fflush(stdout);

  nth = 0;
  for(nthp = 0; nthp < 6; nthp++){
    qopt = 0;
    hit = 0;
    q0 = pp[nthp];

    nthq = 0;
    for(q = -2; q <= 2; q += qdelta){
      nth ++;
      nthq ++;
      pp[nthp] = q;

      GetCosts(mov, seg, R0, Rtmp, pp, costs);

      if(costfile != NULL){
        // write costs to file
        fp = fopen(costfile,"a");
        fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
        fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
        fprintf(fp,"%7d %10.4lf %8.4lf ",
                (int)costs[0],costs[1],costs[2]); // WM  n mean std
        fprintf(fp,"%7d %10.4lf %8.4lf ",
                (int)costs[3],costs[4],costs[5]); // CTX n mean std
        fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
        fprintf(fp,"\n");
        fclose(fp);
      }

      fp = stdout;
      fprintf(fp,"%4d %2d ",nth,nthq);
      fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
      fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
      fprintf(fp,"%7d %10.4lf %8.4lf ",
              (int)costs[0],costs[1],costs[2]); // WM  n mean std
      fprintf(fp,"%7d %10.4lf %8.4lf ",
              (int)costs[3],costs[4],costs[5]); // CTX n mean std
      fprintf(fp,"%8.4lf %8.4lf   %8.4lf ",
              costs[6],costs[7],copt); // t, cost=1/t

      c = costs[7];
      if(c < copt){
        copt = c;
        qopt = q;
        MatrixCopy(Rtmp,R);
        for(n=0; n<8; n++) costsopt[n] = costs[n];
        hit = 1;
      }
      printf("%d\n",hit);
      fflush(stdout);

    }
    if(hit) pp[nthp] = qopt;
    else    pp[nthp] = q0;
    printf("\n");
  } // loop over params

  for(nthp = 0; nthp < 6; nthp++) p[nthp] = pp[nthp];
  for(n=0; n<8; n++) costs[n] = costsopt[n];

  fp = stdout;
  printf("\n");
  fprintf(fp,"final1d ");
  fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
  fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
  fprintf(fp,"%7d %10.4lf %8.4lf ",
          (int)costs[0],costs[1],costs[2]); // WM  n mean std
  fprintf(fp,"%7d %10.4lf %8.4lf ",
          (int)costs[3],costs[4],costs[5]); // CTX n mean std
  fprintf(fp,"%8.4lf %8.4lf   %8.4lf ",costs[6],costs[7],copt); // t, cost=1/t
  printf("\n");
  printf("\n");
  fflush(stdout);

  return(nth);
}
/*---------------------------------------------------------*/
int MinPowell(MRI *mov, MRI *seg, MATRIX *R, double *params,
	      double ftol, double linmintol, int nmaxiters,
	      char *costfile, double *costs, int *niters)
{
  extern MRI *mov_powell;
  extern MRI *seg_powell;
  extern MATRIX *R0_powell;
  extern char *costfile_powell;
  MATRIX *R0;
  float *p, **xi, fret;
  int    r, c, n,  nparams;

  nparams = 6;

  p = vector(1, nparams) ;
  for(n=0; n < nparams; n++) p[n+1] = params[n];

  xi = matrix(1, nparams, 1, nparams) ;
  for (r = 1 ; r <= nparams ; r++) {
    for (c = 1 ; c <= nparams ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  mov_powell = mov;
  seg_powell = seg;
  R0 = MatrixCopy(R,NULL);
  R0_powell = R0;
  costfile_powell = costfile;

  OpenPowell2(p, xi, nparams, ftol, linmintol, nmaxiters, 
	      niters, &fret, compute_powell_cost);
  printf("Powell done niters = %d\n",*niters);

  for(n=0; n < nparams; n++) params[n] = p[n+1];
  GetCosts(mov_powell, seg_powell, R0_powell, R, params, costs);

  free_matrix(xi, 1, nparams, 1, nparams);
  free_vector(p, 1, nparams);
  return(NO_ERROR) ;
}
/*---------------------------------------------------------*/
float compute_powell_cost(float *p) 
{
  extern MRI *mov_powell;
  extern MRI *seg_powell;
  extern MATRIX *R0_powell;
  extern char *costfile_powell;
  extern int nCostEvaluations;
  static MATRIX *R = NULL;
  static double copt = -1;
  double costs[8], pp[6];
  int n, newopt;
  FILE *fp;

  if(R==NULL) R = MatrixAlloc(4,4,MATRIX_REAL);
  for(n=0; n < 6; n++) pp[n] = p[n+1];

  
  if(!UseSurf)  GetCosts(mov_powell,     seg_powell, R0_powell, R, pp, costs);
  else          GetSurfCosts(mov_powell, seg_powell, R0_powell, R, pp, costs);

  // This is for a fast check on convergence
  //costs[7] = 0;
  //for(n=0; n < 6; n++) costs[7] += ((pp[n]-n)*(pp[n]-n)+1);

  if(copt == -1) copt = costs[7];
  newopt = 0;
  if(copt > costs[7]){
    copt = costs[7];
    newopt = 1;
  }


  if(costfile_powell != NULL){
    // write costs to file
    fp = fopen(costfile_powell,"a");
    fprintf(fp,"%4d ",nCostEvaluations);
    fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
    fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
    fprintf(fp,"%7d %10.4lf %8.4lf ",
	    (int)costs[0],costs[1],costs[2]); // WM  n mean std
    fprintf(fp,"%10.4lf %10.4lf %8.4lf ",
	    costs[3],costs[4],costs[5]); // CTX n mean std
    fprintf(fp,"%8.4lf %12.10lf   %12.10lf ",costs[6],costs[7],copt); // t, cost=1/t
    fprintf(fp,"\n");
    fclose(fp);
  }

  if(newopt){
    fp = stdout;
    fprintf(fp,"%4d ",nCostEvaluations);
    fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
    fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
    fprintf(fp,"%7d %10.4lf %8.4lf ",
	    (int)costs[0],costs[1],costs[2]); // WM  n mean std
    fprintf(fp,"%10.4lf %10.4lf %8.4lf ",
	    costs[3],costs[4],costs[5]); // CTX n mean std
    fprintf(fp,"%8.4lf %8.5lf   %8.5lf ",costs[6],costs[7],copt); // t, cost=1/t
    fprintf(fp,"\n");
    fflush(stdout); 
  }

  nCostEvaluations++;
  return((float)costs[7]);
}

/*---------------------------------------------------------------*/
int MRIsegReg(char *subject)
{
  char *SUBJECTS_DIR;
  char tmpstr[2000];
  int c,r,s,n,nwm,nctx,segid,v,z=0,nErode3d;
  MRI *apas, *brain, *wm, *ctx, *segreg;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  printf("Reading aparc+aseg\n");
  sprintf(tmpstr,"%s/%s/mri/aparc+aseg.mgz",SUBJECTS_DIR,subject);
  apas = MRIread(tmpstr);
  if(apas==NULL) exit(1);

  printf("Reading brainmask\n");
  sprintf(tmpstr,"%s/%s/mri/brainmask.mgz",SUBJECTS_DIR,subject);
  brain = MRIread(tmpstr);
  if(brain==NULL) exit(1);

  printf("Zeroing B0 regions\n");
  for(c=0; c < brain->width; c++){
    for(r=0; r < brain->height; r++){
      for(s=0; s < brain->depth; s++){
	v = MRIgetVoxVal(apas,c,r,s,0);
	z = 0;
	if(v == 1006 || v == 2006) z = 1; // entorhinal
	if(v == 1009 || v == 2009) z = 1; // inferiortemp
	if(v == 1012 || v == 2012) z = 1; // lat orb front
	if(v == 1014 || v == 2014) z = 1; // med orb front
	if(v == 1016 || v == 2016) z = 1; // parahipp
	if(v == 1032 || v == 2032) z = 1; // frontal pole
	if(v == 1033 || v == 2033) z = 1; // temporal pole
	if(z) MRIsetVoxVal(brain,c,r,s,0, 0);
      }
    }
  }

  nErode3d = 5;
  printf("Eroding brain mask %d\n",nErode3d);
  for(n=0; n<nErode3d; n++) MRIerode(brain,brain);

  printf("Creating WM mask\n");
  wm = MRIclone(apas,NULL);
  nwm = 0;
  for(c=0; c < brain->width; c++){
    for(r=0; r < brain->height; r++){
      for(s=0; s < brain->depth; s++){
	v = MRIgetVoxVal(brain,c,r,s,0);
	z = 0;
	if(v != 0) {
	  segid = MRIgetVoxVal(apas,c,r,s,0);
	  if(segid == 2 || segid == 41) {
	    z = 41;
	    nwm ++;
	  }
	}
	MRIsetVoxVal(wm,c,r,s,0, z);
      }
    }
  }
  printf("Found %d WM voxels\n",nwm);

  nErode3d = 2;
  printf("Eroding WM mask %d\n",nErode3d);
  for(n=0; n<nErode3d; n++) MRIerode(wm,wm);

  printf("Creating Cortex mask\n");
  ctx = MRIclone(apas,NULL);
  nctx = 0;
  for(c=0; c < brain->width; c++){
    for(r=0; r < brain->height; r++){
      for(s=0; s < brain->depth; s++){
	v = MRIgetVoxVal(brain,c,r,s,0);
	z = 0;
	if(v != 0){
	  segid = MRIgetVoxVal(apas,c,r,s,0);
	  if(segid == 3 || segid == 42 ||
	     (segid >= 1000 && segid <= 1034 ) ||
	     (segid >= 2000 && segid <= 2034 ) ){
	    z = 3;
	    nctx++;
	  }
	}
	MRIsetVoxVal(ctx,c,r,s,0, z);
      }
    }
  }
  printf("Found %d Ctx voxels\n",nctx);

  segreg = MRIadd(wm,ctx,NULL);

  printf("Writing segreg.mgz\n");
  sprintf(tmpstr,"%s/%s/mri/segreg.mgz",SUBJECTS_DIR,subject);
  MRIwrite(segreg,tmpstr);

  MRIfree(&brain);
  MRIfree(&apas);
  MRIfree(&wm);
  MRIfree(&ctx);
  MRIfree(&segreg);

  return(0);
}

/*---------------------------------------------------------------*/
int MRISsegReg(char *subject, int ForceReRun)
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

  char *SUBJECTS_DIR;
  char tmpstr[2000];
  int c,r,s,n,v,z=0;
  MRI *apas, *brain;
  extern MRI *lhCortexLabel, *rhCortexLabel;
  MRI *lhwmmask, *rhwmmask, *lhctxmask, *rhctxmask;
  float  fx, fy, fz;
  int ReRun=0;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  if(UseLH){
    // Load the LH white surface, project it into WM and Ctx
    printf("Loading lh.white surf\n");
    sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
    lhwm = MRISread(tmpstr); // starts as white, projected in
    if(lhwm == NULL) exit(1);
    lhctx = MRISread(tmpstr); // starts as white, projected out

    printf("Loading lh.thickness\n");
    sprintf(tmpstr,"%s/%s/surf/lh.thickness",SUBJECTS_DIR,subject);
    err = MRISreadCurvatureFile(lhwm, tmpstr);
    if(err) exit(1);
    err = MRISreadCurvatureFile(lhctx, tmpstr);

    printf("GM Proj: %d %lf %lf\n",DoGMProjFrac,GMProjFrac,GMProjAbs);
    printf("WM Proj: %d %lf %lf\n",DoWMProjFrac,WMProjFrac,WMProjAbs);
    printf("Projecting LH Surfs\n");
    for(n = 0; n < lhwm->nvertices; n++){
      if(DoWMProjAbs)  ProjNormDist(&fx, &fy, &fz, lhwm,  n, -WMProjAbs);
      if(DoWMProjFrac) ProjNormFracThick(&fx, &fy, &fz, lhwm,  n, -WMProjFrac);
      lhwm->vertices[n].x = fx;
      lhwm->vertices[n].y = fy;
      lhwm->vertices[n].z = fz;
      if(DoGMProjAbs)  ProjNormDist(&fx, &fy, &fz, lhctx,  n, +GMProjAbs);
      if(DoGMProjFrac) ProjNormFracThick(&fx, &fy, &fz, lhctx,  n, +GMProjFrac);
      lhctx->vertices[n].x = fx;
      lhctx->vertices[n].y = fy;
      lhctx->vertices[n].z = fz;
    }
  }

  printf("Projecting RH Surfs\n");
  if(UseRH){
    // Load the RH white surface, project it into WM and Ctx
    printf("Loading rh.white surf\n");
    sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
    rhwm = MRISread(tmpstr);
    if(rhwm == NULL) exit(1);
    rhctx = MRISread(tmpstr);
    
    printf("Loading rh.thickness\n");
    sprintf(tmpstr,"%s/%s/surf/rh.thickness",SUBJECTS_DIR,subject);
    err = MRISreadCurvatureFile(rhwm, tmpstr);
    if(err) exit(1);
    err = MRISreadCurvatureFile(rhctx, tmpstr);
    
    printf("Projecting RH Surfs\n");
    for(n = 0; n < rhwm->nvertices; n++){
      if(DoWMProjAbs)  ProjNormDist(&fx, &fy, &fz, rhwm,  n, -WMProjAbs);
      if(DoWMProjFrac) ProjNormFracThick(&fx, &fy, &fz, rhwm,  n, -WMProjFrac);
      rhwm->vertices[n].x = fx;
      rhwm->vertices[n].y = fy;
      rhwm->vertices[n].z = fz;
      if(DoGMProjAbs)  ProjNormDist(&fx, &fy, &fz, rhctx,  n, +GMProjAbs);
      if(DoGMProjFrac) ProjNormFracThick(&fx, &fy, &fz, rhctx,  n, +GMProjFrac);
      rhctx->vertices[n].x = fx;
      rhctx->vertices[n].y = fy;
      rhctx->vertices[n].z = fz;
    }
  }

  if(!UseMask) return(0);

  // Determine whether segreg needs to be recomputed or not
  ReRun = 0;
  sprintf(tmpstr,"%s/%s/surf/lh.segreg.mgh",SUBJECTS_DIR,subject);
  if(!fio_FileExistsReadable(tmpstr)) ReRun = 1;
  sprintf(tmpstr,"%s/%s/surf/rh.segreg.mgh",SUBJECTS_DIR,subject);
  if(!fio_FileExistsReadable(tmpstr)) ReRun = 1;

  if(!ReRun && !ForceReRun){
    printf("Surf segreg already exists, reading in\n");
    sprintf(tmpstr,"%s/%s/surf/lh.segreg.mgh",SUBJECTS_DIR,subject);
    lhsegmask = MRIread(tmpstr);
    if(lhsegmask == NULL) exit(1);
    sprintf(tmpstr,"%s/%s/surf/rh.segreg.mgh",SUBJECTS_DIR,subject);
    rhsegmask = MRIread(tmpstr);
    if(rhsegmask == NULL) exit(1);
    return(0); // Return here
  }
  if(!ReRun && ForceReRun)
    printf("Surf segreg already exists, but rerun forced\n");
  else
    printf("Surf segreg does not exist ... computing \n");


  /* First, create a mask by zeroing the brainmask in dropout regions
     then eroding that by 5 voxels. This gets rid of both the edges of
     the brain and areas near the dropout.
   */

  printf("Reading aparc+aseg\n");
  sprintf(tmpstr,"%s/%s/mri/aparc+aseg.mgz",SUBJECTS_DIR,subject);
  apas = MRIread(tmpstr);
  if(apas==NULL) exit(1);
    
  printf("Reading brainmask\n");
  sprintf(tmpstr,"%s/%s/mri/brainmask.mgz",SUBJECTS_DIR,subject);
  brain = MRIread(tmpstr);
  if(brain==NULL) exit(1);
    
  printf("Zeroing B0 regions\n");
  for(c=0; c < brain->width; c++){
    for(r=0; r < brain->height; r++){
      for(s=0; s < brain->depth; s++){
	if(MRIgetVoxVal(brain,c,r,s,0) < 0.5) continue;
	v = MRIgetVoxVal(apas,c,r,s,0);
	z = 0;
	if(v == 1006 || v == 2006) z = 1; // entorhinal
	if(v == 1009 || v == 2009) z = 1; // inferiortemp
	if(v == 1012 || v == 2012) z = 1; // lat orb front
	if(v == 1014 || v == 2014) z = 1; // med orb front
	if(v == 1016 || v == 2016) z = 1; // parahipp
	if(v == 1032 || v == 2032) z = 1; // frontal pole
	if(v == 1033 || v == 2033) z = 1; // temporal pole
	if(z) MRIsetVoxVal(brain,c,r,s,0, 0);
	else  MRIsetVoxVal(brain,c,r,s,0, 1);
      }
    }
  }

  printf("Eroding mask %d\n",5);
  for(n=0; n < 5; n++) MRIerode(brain,brain);

  // Now sample the new mask on the both the wm and ctx LH surfaces
  lhwmmask = vol2surf_linear(brain,NULL,NULL,NULL,NULL,
			     lhwm, 0,SAMPLE_NEAREST, FLT2INT_ROUND, NULL, 0, 1);
  lhctxmask = vol2surf_linear(brain,NULL,NULL,NULL,NULL,
			      lhctx, 0,SAMPLE_NEAREST, FLT2INT_ROUND, NULL, 0, 1);

  // Create the final mask by forcing each vertex to be in both
  // the wm and ctx masks as well as in the cortex label
  lhsegmask = MRIalloc(lhwm->nvertices,1,1,MRI_INT);
  for(n = 0; n < lhwm->nvertices; n++){
    if(MRIgetVoxVal(lhCortexLabel,n,0,0,0) < 0.5 ||
       MRIgetVoxVal(lhwmmask,n,0,0,0) < 0.5 ||
       MRIgetVoxVal(lhctxmask,n,0,0,0) < 0.5)
      MRIsetVoxVal(lhsegmask,n,0,0,0,  0);
    else
      MRIsetVoxVal(lhsegmask,n,0,0,0,  1);
  }

  // Finally, Write out the mask 
  sprintf(tmpstr,"%s/%s/surf/lh.segreg.mgh",SUBJECTS_DIR,subject);
  MRIwrite(lhsegmask,tmpstr);

  // -----------------------------------------------------
  // Do the same thing for the RH 
  rhwmmask = vol2surf_linear(brain,NULL,NULL,NULL,NULL,
			     rhwm, 0,SAMPLE_NEAREST, FLT2INT_ROUND, NULL, 0, 1);

  rhctxmask = vol2surf_linear(brain,NULL,NULL,NULL,NULL,
			      rhctx, 0,SAMPLE_NEAREST, FLT2INT_ROUND, NULL, 0, 1);

  rhsegmask = MRIalloc(rhwm->nvertices,1,1,MRI_FLOAT);
  for(n = 0; n < rhwm->nvertices; n++){
    if(MRIgetVoxVal(rhCortexLabel,n,0,0,0) < 0.5 ||
       MRIgetVoxVal(rhwmmask,n,0,0,0) < 0.5 ||
       MRIgetVoxVal(rhctxmask,n,0,0,0) < 0.5)
      MRIsetVoxVal(rhsegmask,n,0,0,0,  0);
    else
      MRIsetVoxVal(rhsegmask,n,0,0,0,  1);
  }
  sprintf(tmpstr,"%s/%s/surf/rh.segreg.mgh",SUBJECTS_DIR,subject);
  MRIwrite(rhsegmask,tmpstr);

  //-----------------------------------------------------------------

  if(0) {
    MRIwrite(brain,"be5.mgh");
    MRIwrite(lhwmmask,"lh.wm.b5.mgh");
    MRIwrite(lhctxmask,"lh.ctx.b5.mgh");
    MRIwrite(rhctxmask,"rh.ctx.b5.mgh");
    MRIwrite(rhwmmask,"rh.wm.b5.mgh");
    MRISwrite(rhwm,"./rh.wm");
    MRISwrite(rhctx,"./rh.ctx");
    MRISwrite(lhwm,"./lh.wm");
    MRISwrite(lhctx,"./lh.ctx");
  }
  MRIfree(&lhwmmask);
  MRIfree(&lhctxmask);
  MRIfree(&rhctxmask);
  MRIfree(&rhwmmask);
  MRIfree(&brain);
  MRIfree(&apas);

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
  if(sign ==  0) a = -abs(slope*(d-center)); // not sure this is useful
  if(sign == -1) a = -(slope*(d-center));
  if(sign == +1) a = +(slope*(d-center));
  c = 1+tanh(a);
  *pct = d;
  return(c);
}

/*-------------------------------------------------------*/
double *GetSurfCosts(MRI *mov, MRI *notused, MATRIX *R0, MATRIX *R,
		     double *p, double *costs)
{
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
  double angles[3],d,dsum,dsum2,dstd,dmean,vwm,vctx,c,csum,csum2,cstd,cmean;
  MATRIX *Mrot=NULL, *Mtrans=NULL;
  MRI *vlhwm, *vlhctx, *vrhwm, *vrhctx;
  int nhits,n;
  //FILE *fp;

  if(R==NULL){
    printf("ERROR: GetSurfCosts(): R cannot be NULL\n");
    return(NULL);
  }

  Mtrans = MatrixIdentity(4,NULL);
  Mtrans->rptr[1][4] = p[0];
  Mtrans->rptr[2][4] = p[1];
  Mtrans->rptr[3][4] = p[2];

  angles[0] = p[3]*(M_PI/180);
  angles[1] = p[4]*(M_PI/180);
  angles[2] = p[5]*(M_PI/180);
  Mrot = MRIangles2RotMat(angles);

  // R = Mtrans*Mrot*R0
  R = MatrixMultiply(Mrot,R0,R);
  R = MatrixMultiply(Mtrans,R,R);

  //printf("Trans: %g %g %g\n",p[0],p[1],p[2]);
  //printf("Rot:   %g %g %g\n",p[3],p[4],p[5]);

  if(UseLH){
    vlhwm = vol2surf_linear(mov,NULL,NULL,NULL,R,
			    lhwm, 0,interpcode, 
			    FLT2INT_ROUND, NULL, 0, nsubsamp);
    vlhctx = vol2surf_linear(mov,NULL,NULL,NULL,R,
			     lhctx, 0,interpcode, 
			     FLT2INT_ROUND, NULL, 0, nsubsamp);
    if(lhcost == NULL) lhcost = MRIclone(vlhctx,NULL);
    if(lhcon == NULL)  lhcon  = MRIclone(vlhctx,NULL);
  }
  if(UseRH){
    vrhwm = vol2surf_linear(mov,NULL,NULL,NULL,R,
			    rhwm, 0,interpcode, 
			    FLT2INT_ROUND, NULL, 0, nsubsamp);
    vrhctx = vol2surf_linear(mov,NULL,NULL,NULL,R,
			     rhctx, 0,interpcode, 
			     FLT2INT_ROUND, NULL, 0, nsubsamp);
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
      if(lhcostfile || lhcost0file) MRIsetVoxVal(lhcost,n,0,0,0,0.0);
      if(lhconfile)  MRIsetVoxVal(lhcon,n,0,0,0,0.0);
      if(lhCortexLabel && MRIgetVoxVal(lhCortexLabel,n,0,0,0) < 0.5) continue;
      if(UseMask && MRIgetVoxVal(lhsegmask,n,0,0,0) < 0.5) continue;
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
      if(rhcostfile || rhcost0file) MRIsetVoxVal(rhcost,n,0,0,0,0.0);
      if(rhconfile)  MRIsetVoxVal(rhcon,n,0,0,0,0.0);
      if(rhCortexLabel && MRIgetVoxVal(rhCortexLabel,n,0,0,0) < 0.5) continue;
      if(UseMask && MRIgetVoxVal(rhsegmask,n,0,0,0) < 0.5) continue;
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

  MatrixFree(&Mrot);
  MatrixFree(&Mtrans);

  if(UseLH){
    MRIfree(&vlhwm);
    MRIfree(&vlhctx);
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
    MRIfree(&vrhwm);
    MRIfree(&vrhctx);
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

