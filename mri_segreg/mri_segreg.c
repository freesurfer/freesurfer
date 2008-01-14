/**
 * @file  mri_segreg.c
 * @brief compute/optimize cost function of segmentation-based registration
 *
 */
/*
 * Original Author: Greg Grev
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2008/01/14 07:17:10 $
 *    $Revision: 1.32 $
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

  --reg regfile
  --mov fvol
  --inorm inorm : intensity profile in reg with anat

  --o out : save final output
  --out-reg outreg : reg at lowest cost
  --cost costfile
  --sum sumfile : def is outreg.sum
  --no-surf : do not use surface-based method
  --1dpreopt min max delta : brute force in PE direction

  --frame nthframe : use given frame in input (default = 0)
  --mid-frame : use use middle frame

  --interp interptype : interpolation trilinear or nearest (def is trilin)
  --no-crop: do not crop anat (crops by default)
  --profile : print out info about exec time

  --noise stddev : add noise with stddev to input for testing sensitivity
  --seed randseed : for use with --noise

  --aseg : use aseg instead of segreg.mgz

  --nmax nmax   : max number of powell iterations (def 36)
  --tol   tol   : powell inter-iteration tolerance on cost
  --tol1d tol1d : tolerance on powell 1d minimizations

  --1dmin : use brute force 1D minimizations instead of powell
  --n1dmin n1dmin : number of 1d minimization (default = 3)

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


int main(int argc, char *argv[]) ;

static char vcid[] =
"$Id: mri_segreg.c,v 1.32 2008/01/14 07:17:10 greve Exp $";
char *Progname = NULL;

int debug = 0, gdiagno = -1;

char *movvolfile=NULL;
char *regfile=NULL;
char *outregfile=NULL;
char *sumfile=NULL;

char *interpmethod = "trilinear";
int   interpcode = 0;
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
double TolPowell = .0001;
double LinMinTolPowell = .1;

int MkSegReg = 0;

#define NMAX 100
int ntx=0, nty=0, ntz=0, nax=0, nay=0, naz=0;
double txlist[NMAX],tylist[NMAX],tzlist[NMAX];
double axlist[NMAX],aylist[NMAX],azlist[NMAX];

int UseSurf = 1;
MRIS *lhwm, *rhwm, *lhctx, *rhctx;
double SurfProj = 0.5;
MRI *lhsegmask, *rhsegmask;
char *cmdline2;
struct utsname uts;

int PenaltySign  = -1;
double PenaltySlope = .5;
int DoMidFrame = 0;

int Do1DPreOpt = 0;
double PreOptMin, PreOptMax, PreOptDelta, PreOpt, PreOptAtMin;

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  char cmdline[CMD_LINE_LEN] ;
  double costs[8], mincost, p[6];
  double tx, ty, tz, ax, ay, az;
  int nth,nthtx, nthty, nthtz, nthax, nthay, nthaz, ntot, n, err;
  MATRIX *Ttemp=NULL, *invTtemp=NULL, *Stemp=NULL, *invStemp=NULL;
  MATRIX *R=NULL, *Rcrop=NULL, *R00=NULL, *Rdiff=NULL;
  MATRIX *Rmin=NULL, *Scrop=NULL, *invScrop=NULL, *invTcrop=NULL, *Tcrop=NULL;
  struct timeb  mytimer;
  double secCostTime;
  MRI_REGION box;
  FILE *fp;

  make_cmd_version_string
    (argc, argv,
     "$Id: mri_segreg.c,v 1.32 2008/01/14 07:17:10 greve Exp $",
     "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mri_segreg.c,v 1.32 2008/01/14 07:17:10 greve Exp $",
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

  parse_commandline(argc, argv);
  if(gdiagno > -1) Gdiag_no = gdiagno;
  check_options();
  dump_options(stdout);

  fp = fopen(sumfile,"w");
  dump_options(fp);
  fflush(fp);

  R00 = MatrixCopy(R0,NULL);

  printf("Loading mov\n");
  mov = MRIread(movvolfile);
  if (mov == NULL) exit(1);

  if(DoMidFrame) frame = nint(mov->nframes/2);

  if(mov->nframes > 1){
    printf("Extracting frame %d\n",frame);
    mritmp = fMRIframe(mov, frame, NULL);
    MRIfree(&mov);
    mov = mritmp;
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
    if (SynthSeed < 0) SynthSeed = PDFtodSeed();
    srand48(SynthSeed);
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
  else          GetSurfCosts(mov, segreg, R0, R, p, costs);

  printf("Initial cost is %lf\n",costs[7]);
  printf("Initial Costs\n");
  printf("%7d %10.4lf %8.4lf ",
         (int)costs[0],costs[1],costs[2]); // WM  n mean std
  printf("%10.4lf %10.4lf %8.4lf ",
         costs[3],costs[4],costs[5]); // CTX n mean std
  printf("%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
  printf("\n");

  fprintf(fp,"Initial cost is %lf\n",costs[7]);
  fprintf(fp,"Initial Costs\n");
  fprintf(fp,"%7d %10.4lf %8.4lf ",
         (int)costs[0],costs[1],costs[2]); // WM  n mean std
  printf("%10.4lf %10.4lf %8.4lf ",
         costs[3],costs[4],costs[5]); // CTX n mean std
  fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
  fprintf(fp,"\n");
  fflush(fp);

  //exit(1);

  if(Do1DPreOpt){
    printf("\n");
    printf("Performing 1D preopt %g %g %g\n",PreOptMin,PreOptMax,PreOptDelta);
    fprintf(fp,"\n");
    fprintf(fp,"Performing 1D preopt %g %g %g\n",PreOptMin,PreOptMax,PreOptDelta);
    mincost = 10e10;
    PreOptAtMin = 0;
    for(nth=0; nth < 6; nth++) p[nth] = 0.0;
    nth = 0;
    for(PreOpt = PreOptMin; PreOpt <= PreOptMax; PreOpt += PreOptDelta){
      p[2] = PreOpt;
      GetSurfCosts(mov, segreg, R0, R, p, costs);
      if(costs[7] < mincost) {
	mincost = costs[7];
	PreOptAtMin = PreOpt;
      }
      printf("%8.4lf %8.4lf %8.4lf \n",PreOpt,costs[7],mincost);
      fprintf(fp,"%8.4lf %8.4lf %8.4lf \n",PreOpt,costs[7],mincost);
      //sprintf(tmpstr,"myreg.%02d.dat",nth);
      //regio_write_register(tmpstr,subject,mov->xsize,
      //mov->zsize,intensity,R,FLT2INT_ROUND);
      nth ++;
    }
    p[2] = PreOptAtMin; // phase encode direction
    GetSurfCosts(mov, segreg, R0, R, p, costs);
    //regio_write_register("myreg.reg",subject,mov->xsize,
    //		 mov->zsize,intensity,R,FLT2INT_ROUND);
    MatrixCopy(R0,R);
    printf("\n");
    fprintf(fp,"\n");
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
  if(!UseSurf)  GetCosts(mov,     segreg, R0, R, p, costs);
  else          GetSurfCosts(mov, segreg, R0, R, p, costs);
  

  printf("Min cost was %lf\n",costs[7]);
  printf("Number of iterations %5d\n",nth);
  printf("Optmization time %lf sec\n",secCostTime);
  printf("Parameters at optimum (transmm, rotdeg)\n");
  printf("%7.3lf %7.3lf %7.3lf %6.3lf %6.3lf %6.3lf \n",
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
  
  fprintf(fp,"Min cost was %lf\n",costs[7]);
  fprintf(fp,"Number of iterations %5d\n",nth);
  fprintf(fp,"Optmization time %lf sec\n",secCostTime);
  fprintf(fp,"Parameters at optimum (transmm, rotdeg)\n");
  fprintf(fp,"%7.3lf %7.3lf %7.3lf %6.3lf %6.3lf %6.3lf \n",
	 p[0],p[1],p[2],p[3],p[4],p[5]);
  fprintf(fp,"Costs at optimum\n");
  fprintf(fp,"%7d %10.4lf %8.4lf ",
	 (int)costs[0],costs[1],costs[2]); // WM  n mean std
  fprintf(fp,"%10.4lf %10.4lf %8.4lf ",
	 costs[3],costs[4],costs[5]); // CTX n mean std
  fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
  fprintf(fp,"\n");
  
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
    MRIvol2VolTkReg(mov, segreg0, R, interpcode, sinchw);
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
  char **pargv, *option ;
  int err,nv,n;
  double vmin, vmax, vdelta;

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
    else if (!strcasecmp(option, "--1dmin"))     DoPowell = 0;
    else if (!strcasecmp(option, "--mid-frame")) DoMidFrame = 1;
    else if (!strcasecmp(option, "--surf")){
      UseSurf = 1;
      DoCrop = 0;
    }
    else if (!strcasecmp(option, "--no-surf"))   UseSurf = 0;
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
    } else if (istringnmatch(option, "--nmax",0)) {
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
    } else if (istringnmatch(option, "--o",0)) {
      if (nargc < 1) argnerr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--reg",0)) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      err = regio_read_register(regfile, &subject, &ipr, &bpr,
                                &intensity, &R0, &float2int);
      if (err) exit(1);
      nargsused = 1;
    } else if (istringnmatch(option, "--out-reg",0)) {
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
    } else if (istringnmatch(option, "--wm-gt-gm",0)) {
      if (nargc < 1) argnerr(option,1);
      PenaltySign = +1;
      sscanf(pargv[0],"%lf",&PenaltySlope);
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
    } else if (istringnmatch(option, "--interp",8)) {
      if (nargc < 1) argnerr(option,1);
      interpmethod = pargv[0];
      nargsused = 1;
      if (!strcmp(interpmethod,"sinc") && CMDnthIsArg(nargc, pargv, 1)) {
        sscanf(pargv[1],"%d",&sinchw);
        nargsused ++;
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
printf("\n");
printf("  mri_segreg\n");
printf("\n");
printf("  --reg regfile\n");
printf("  --mov fvol\n");
printf("  --inorm inorm : intensity profile in reg with anat\n");
printf("\n");
printf("  --o out : save final output\n");
printf("  --out-reg outreg : reg at lowest cost\n");
printf("  --cost costfile\n");
printf("  --sum sumfile : def is outreg.sum\n");
printf("  --no-surf : do not use surface-based method\n");
printf("  --1dpreopt min max delta : brute force in PE direction\n");
printf("\n");
printf("  --frame nthframe : use given frame in input (default = 0)\n");
printf("  --mid-frame : use middle frame\n");
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
printf("  --tol1d tol1d : tolerance on powell 1d minimizations\n");
printf("\n");
printf("  --1dmin : use brute force 1D minimizations instead of powell\n");
printf("  --n1dmin n1dmin : number of 1d minimization (default = 3)\n");
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
  if(regfile == NULL) {
    printf("ERROR: need --reg.\n");
    exit(1);
  }
  if(outregfile == NULL) {
    printf("ERROR: need --out-reg.\n");
    exit(1);
  }
  if (SegRegCostFile == NULL) {
    printf("ERROR: need --cost.\n");
    exit(1);
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
  fprintf(fp,"PenaltySign  %d\n",PenaltySign);
  fprintf(fp,"PenaltySlope %lf\n",PenaltySlope);
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
  static MATRIX *R = NULL;
  static int nth=0;
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
    fprintf(fp,"%4d ",nth);
    fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
    fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
    fprintf(fp,"%7d %10.4lf %8.4lf ",
	    (int)costs[0],costs[1],costs[2]); // WM  n mean std
    fprintf(fp,"%10.4lf %10.4lf %8.4lf ",
	    costs[3],costs[4],costs[5]); // CTX n mean std
    fprintf(fp,"%8.4lf %8.5lf   %8.5lf ",costs[6],costs[7],copt); // t, cost=1/t
    fprintf(fp,"\n");
    fclose(fp);
  }

  if(newopt){
    fp = stdout;
    fprintf(fp,"%4d ",nth);
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

  nth++;
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
  extern double SurfProj;
  char *SUBJECTS_DIR;
  char tmpstr[2000];
  int c,r,s,n,v,z=0;
  LABEL *label;
  MRI *apas, *brain;
  MRI *lhCortexLabel, *rhCortexLabel, *lhwmmask, *rhwmmask, *lhctxmask, *rhctxmask;
  float  fx, fy, fz;
  int ReRun=0;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  // Load the LH white surface, project it into WM and Ctx
  printf("Loading lh.white surf\n");
  sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
  lhwm = MRISread(tmpstr);
  if(lhwm == NULL) exit(1);
  lhctx = MRISread(tmpstr);

  printf("Loading lh.thickness\n");
  sprintf(tmpstr,"%s/%s/surf/lh.thickness",SUBJECTS_DIR,subject);
  err = MRISreadCurvatureFile(lhwm, tmpstr);
  if(err) exit(1);
  err = MRISreadCurvatureFile(lhctx, tmpstr);

  printf("Projecting LH Surfs by %g\n",SurfProj);
  for(n = 0; n < lhwm->nvertices; n++){
    ProjNormFracThick(&fx, &fy, &fz, lhwm,  n, -SurfProj);
    lhwm->vertices[n].x = fx;
    lhwm->vertices[n].y = fy;
    lhwm->vertices[n].z = fz;
    ProjNormFracThick(&fx, &fy, &fz, lhctx, n, +SurfProj);
    lhctx->vertices[n].x = fx;
    lhctx->vertices[n].y = fy;
    lhctx->vertices[n].z = fz;
  }

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

  printf("Projecting RH Surfs by %g\n",SurfProj);
  for(n = 0; n < rhwm->nvertices; n++){
    ProjNormFracThick(&fx, &fy, &fz, rhwm,n,-SurfProj);
    rhwm->vertices[n].x = fx;
    rhwm->vertices[n].y = fy;
    rhwm->vertices[n].z = fz;
    ProjNormFracThick(&fx, &fy, &fz, rhctx,n,+SurfProj);
    rhctx->vertices[n].x = fx;
    rhctx->vertices[n].y = fy;
    rhctx->vertices[n].z = fz;
  }

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
			     lhwm, 0,SAMPLE_NEAREST, FLT2INT_ROUND, NULL, 0);
  lhctxmask = vol2surf_linear(brain,NULL,NULL,NULL,NULL,
			      lhctx, 0,SAMPLE_NEAREST, FLT2INT_ROUND, NULL, 0);

  // Load the LH cortex label
  label = LabelRead(subject, "lh.cortex.label");
  if(label == NULL) exit(1);
  lhCortexLabel = MRISlabel2Mask(lhwm, label, NULL);
  LabelFree(&label);

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
  // Load the RH cortex label
  label = LabelRead(subject, "rh.cortex.label");
  if(label == NULL) exit(1);
  rhCortexLabel = MRISlabel2Mask(rhwm, label, NULL);
  LabelFree(&label);

  rhwmmask = vol2surf_linear(brain,NULL,NULL,NULL,NULL,
			     rhwm, 0,SAMPLE_NEAREST, FLT2INT_ROUND, NULL, 0);

  rhctxmask = vol2surf_linear(brain,NULL,NULL,NULL,NULL,
			      rhctx, 0,SAMPLE_NEAREST, FLT2INT_ROUND, NULL, 0);

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
  MRIfree(&lhCortexLabel);
  MRIfree(&rhctxmask);
  MRIfree(&rhwmmask);
  MRIfree(&rhCortexLabel);
  MRIfree(&brain);
  MRIfree(&apas);

  return(0);
}

/*-------------------------------------------------------*/
double *GetSurfCosts(MRI *mov, MRI *notused, MATRIX *R0, MATRIX *R,
		     double *p, double *costs)
{
  extern MRI *lhsegmask, *rhsegmask;
  extern MRIS *lhwm, *rhwm, *lhctx, *rhctx;
  extern int PenaltySign;
  extern double PenaltySlope;
  double angles[3],d,dsum,dsum2,dstd,dmean,vwm,vctx,c,csum,csum2,cstd,cmean,a=0;
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

  vlhwm = vol2surf_linear(mov,NULL,NULL,NULL,R,
			  lhwm, 0,SAMPLE_TRILINEAR, 
			  FLT2INT_ROUND, NULL, 0);
  vlhctx = vol2surf_linear(mov,NULL,NULL,NULL,R,
			  lhctx, 0,SAMPLE_TRILINEAR, 
			  FLT2INT_ROUND, NULL, 0);
  vrhwm = vol2surf_linear(mov,NULL,NULL,NULL,R,
			  rhwm, 0,SAMPLE_TRILINEAR, 
			  FLT2INT_ROUND, NULL, 0);
  vrhctx = vol2surf_linear(mov,NULL,NULL,NULL,R,
			  rhctx, 0,SAMPLE_TRILINEAR, 
			  FLT2INT_ROUND, NULL, 0);

  for(n = 0; n < 8; n++) costs[n] = 0;

  dsum = 0.0;
  dsum2 = 0.0;
  csum = 0.0;
  csum2 = 0.0;
  nhits = 0;

  //fp = fopen("tmp.dat","w");
  for(n = 0; n < lhwm->nvertices; n++){
    if(MRIgetVoxVal(lhsegmask,n,0,0,0) < 0.5) continue;
    vwm = MRIgetVoxVal(vlhwm,n,0,0,0);
    if(vwm == 0.0) continue;
    vctx = MRIgetVoxVal(vlhctx,n,0,0,0);
    if(vctx == 0.0) continue;
    nhits++;
    costs[1] += vwm;
    costs[2] += (vwm*vwm);
    costs[4] += vctx;
    costs[5] += (vctx*vctx);
    d = 100*(vctx-vwm)/((vctx+vwm)/2.0);
    dsum += d;
    dsum2 += (d*d);
    if(PenaltySign ==  0) a = -abs(PenaltySlope*d);
    if(PenaltySign == -1) a =    -(PenaltySlope*d);
    if(PenaltySign == +1) a =    +(PenaltySlope*d);
    c = 1+tanh(a);
    csum += c;
    csum2 += (c*c);
    //fprintf(fp,"%6d %lf %lf %lf %lf %lf %lf %lf\n",nhits,vwm,vctx,d,dsum,a,c,csum);
  }
  for(n = 0; n < rhwm->nvertices; n++){
    if(MRIgetVoxVal(rhsegmask,n,0,0,0) < 0.5) continue;
    vwm = MRIgetVoxVal(vrhwm,n,0,0,0);
    if(vwm == 0.0) continue;
    vctx = MRIgetVoxVal(vrhctx,n,0,0,0);
    if(vctx == 0.0) continue;
    nhits++;
    costs[1] += vwm;
    costs[2] += (vwm*vwm);
    costs[4] += vctx;
    costs[5] += (vctx*vctx);
    d = 100*(vctx-vwm)/((vctx+vwm)/2.0);
    dsum += d;
    dsum2 += (d*d);
    if(PenaltySign ==  0) a = -abs(PenaltySlope*d); // not sure this is useful
    if(PenaltySign == -1) a =    -(PenaltySlope*d);
    if(PenaltySign == +1) a =    +(PenaltySlope*d);
    c = 1+tanh(a);
    csum += c;
    csum2 += (c*c);
    //fprintf(fp,"%6d %lf %lf %lf %lf %lf %lf %lf\n",nhits,vwm,vctx,d,dsum,a,c,csum);
  }
  //fclose(fp);

  dmean = dsum/nhits;
  dstd  = sum2stddev(dsum,dsum2,nhits);
  //printf("dsum=%g, dsum2=%g, dstd=%g\n",dsum,dsum2,dstd);
  cmean = csum/nhits;
  cstd  = sum2stddev(csum,csum2,nhits);
  //printf("csum=%g, csum2=%g, cstd=%g\n",csum,csum2,cstd);

  costs[0] = nhits;
  costs[2] = sum2stddev(costs[1],costs[2],nhits);
  costs[1] = costs[1]/nhits;
  costs[3] = dstd;
  costs[5] = sum2stddev(costs[4],costs[5],nhits);
  costs[4] = costs[4]/nhits;
  costs[6] = dmean;
  costs[7] = cmean;

  MatrixFree(&Mrot);
  MatrixFree(&Mtrans);
  MRIfree(&vlhwm);
  MRIfree(&vlhctx);
  MRIfree(&vrhwm);
  MRIfree(&vrhctx);

  return(costs);
}

