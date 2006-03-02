/*
BEGINHELP

Estimates the smoothness of a volume-based data set. The smoothness is
measured as the Full-Width-Half-Max (FWHM) in mm.

--i inputvol

Input data. Format must be something readable by mri_convert
(eg, mgh, mgz, img, nii). Alternately, one can synthesize
white gaussian noise with --synth and --synth-frames.

--mask maskfile

Compute FWHM only over voxels in the given mask. Format can be anything
accepted by mri_convert.

--mask-inv

Invert mask, ie, compute FWHM only over voxels outside the given mask.

--X x.mat

Detrend data with the matrix in x.mat. Ie, y = (I-inv(X'*X)*X')*y, where
y is the input. x.mat must be a matlab4 matrix.

--sum sumfile

Prints ascii summary to sumfile.

--fwhm fwhm

Smooth by fwhm mm before estimating the fwhm. This is mainly good for 
debuggging. But with --out can also be used to smooth data.

--o outvol

Save (possibly synthesized and/or smoothed) data to outfile. Automatically
detects format. Format must be one accepted as by mri_convert. 

--synth 

Synthesize input with white gaussian noise. Ten frames are used by default,
but this can be changed with --synth-frames.

--synth-frames nframes

Synthesize input with white gaussian noise with the given number of frames.
Implies --synth.


ENDHELP
*/

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
#include "fmriutils.h"
#include "mrisurf.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "annotation.h"
#include "cmdargs.h"
#include "timer.h"
#include "matfile.h"
#include "randomfields.h"
#include "icosahedron.h"
#include "pdf.h"
#include "matfile.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_fwhm.c,v 1.2 2006/03/02 00:40:32 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *inpath=NULL;
char *outpath=NULL;
char *sumfile=NULL;
MRI *InVals=NULL;

char *maskpath=NULL;
MRI *mask=NULL;
int maskinv = 0;
double maskthresh = 0.5;
char *outmaskpath=NULL;

MRI *mritmp=NULL;
char tmpstr[2000];

double infwhm = 0, ingstd = 0;
int synth = 0, nframes = 10;
int SynthSeed = -1;

char *Xfile=NULL;
MATRIX *X=NULL;
int DetrendOrder = -1;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int nargs, n, Ntp, nsearch;
  double fwhm = 0, nresels, voxelvolume, nvoxperresel;
  double car1mn, rar1mn,sar1mn,cfwhm,rfwhm,sfwhm, ftmp; 
  FILE *fp;

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

  if(SynthSeed < 0) SynthSeed = PDFtodSeed();
  if(debug) dump_options(stdout);

  if(!synth){
    InVals = MRIread(inpath);
    if(InVals == NULL) exit(1);
    if(InVals->type != MRI_FLOAT){
      mritmp = MRISeqchangeType(InVals, MRI_FLOAT, 0, 0, 0);
      MRIfree(&InVals);
      InVals = mritmp;
    }
  }
  else{
    InVals = MRIreadHeader(inpath,MRI_VOLUME_TYPE_UNKNOWN);
    if(InVals == NULL) exit(1);
    if(InVals->type != MRI_FLOAT) InVals->type = MRI_FLOAT;
    printf("Synthesizing %d frames, Seed = %d\n",nframes,SynthSeed);
    InVals = MRIrandn(InVals->width, InVals->height, InVals->depth, nframes,0, 1, NULL);
  }
  voxelvolume = InVals->xsize * InVals->ysize * InVals->zsize ;
  printf("voxelvolume %g mm3\n",voxelvolume);

  if(maskpath){
    printf("Loading mask %s\n",maskpath);
    mask = MRIread(maskpath);
    if(mask==NULL) exit(1);
  }
  if(mask){
    MRIbinarize(mask, mask, maskthresh, 0, 1);
    if(maskinv) MRImaskInvert(mask,mask);
    nsearch = MRInMask(mask);
  }
  else nsearch = InVals->width * InVals->height * InVals->depth;
  printf("Search region is %d voxels = %lf mm3\n",nsearch,nsearch*voxelvolume);

  if(outmaskpath) MRIwrite(mask,outmaskpath);

  if(DetrendOrder >= 0){
    Ntp = InVals->nframes;
    printf("Polynomial detrending, order = %d\n",DetrendOrder);
    X = MatrixAlloc(Ntp,DetrendOrder+1,MATRIX_REAL);
    for(n=0;n<Ntp;n++) X->rptr[n+1][1] = 1.0;
    ftmp = Ntp/2.0;
    if(DetrendOrder >= 1)
      for(n=0;n<Ntp;n++) X->rptr[n+1][2] = (n-ftmp)/ftmp;
    if(DetrendOrder >= 2)
      for(n=0;n<Ntp;n++) X->rptr[n+1][3] = pow((n-ftmp),2.0)/(ftmp*ftmp);
  }

  if(X){
    printf("Detrending\n");
    if(X->rows != InVals->nframes){
      printf("ERROR: dimension mismatch between X and input\n");
      exit(1);
    }
    mritmp = fMRIdetrend(InVals,X);
    if(mritmp == NULL) exit(1);
    MRIfree(&InVals);
    InVals = mritmp;
  }

  if(infwhm > 0){
    printf("Smoothing input by fwhm=%lf, gstd=%lf\n",infwhm,ingstd);
    MRImaskedGaussianSmooth(InVals, mask, ingstd, InVals);
  }

  printf("Computing spatial AR1 in volume.\n");
  fMRIspatialAR1Mean(InVals, mask, &car1mn, &rar1mn, &sar1mn);
  cfwhm = RFar1ToFWHM(car1mn, InVals->xsize);
  rfwhm = RFar1ToFWHM(rar1mn, InVals->ysize);
  sfwhm = RFar1ToFWHM(sar1mn, InVals->zsize);
  fwhm = sqrt((cfwhm*cfwhm + rfwhm*rfwhm + sfwhm*sfwhm)/3.0);
  printf("ar1mn = (%lf,%lf,%lf)\n",car1mn,rar1mn,sar1mn);
  printf("colfwhm   = %lf\n",cfwhm);
  printf("rowfwhm   = %lf\n",rfwhm);
  printf("slicefwhm = %lf\n",sfwhm);
  printf("fwhm = %lf\n",fwhm);

  nvoxperresel = pow(fwhm,3.0)/voxelvolume;
  nresels = voxelvolume*nsearch/(pow(fwhm,3.0));
  printf("nresels %lf\n",nresels);
  printf("nvoxperresel %lf\n",nvoxperresel);

  fflush(stdout);

  if(sumfile){
    fp = fopen(sumfile,"w");
    if(fp == NULL){
      printf("ERROR: opening %s\n",sumfile);
      exit(1);
    }
    dump_options(fp);
    fprintf(fp,"fwhm %lf\n",fwhm);
    fclose(fp);
  }

  if(outpath) MRIwrite(InVals,outpath);

  return 0;
}
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
    else if (!strcasecmp(option, "--synth")) synth = 1;
    else if (!strcasecmp(option, "--nosynth")) synth = 0;
    else if (!strcasecmp(option, "--mask-inv")) maskinv = 1;

    else if (!strcasecmp(option, "--i")){
      if(nargc < 1) CMDargNErr(option,1);
      inpath = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mask")){
      if(nargc < 1) CMDargNErr(option,1);
      maskpath = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mask-thresh")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&maskthresh);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--out-mask")){
      if(nargc < 1) CMDargNErr(option,1);
      outmaskpath = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--sum")){
      if(nargc < 1) CMDargNErr(option,1);
      sumfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--fwhm")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&infwhm);
      ingstd = infwhm/sqrt(log(256.0));
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--synth-frames")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nframes);
      synth = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      outpath = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--X")){
      if(nargc < 1) CMDargNErr(option,1);
      Xfile = pargv[0];
      //X = MatrixReadTxt(Xfile, NULL);
      X = MatlabRead(Xfile);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--detrend")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&DetrendOrder);
      if(DetrendOrder > 2){
	printf("ERROR: cannot have detrending order > 2\n");
	exit(1);
      }
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
  printf("USAGE: %s\n",Progname) ;
  printf("\n");
  printf("   --i inputvol\n");
  printf("   --mask maskvol : binary mask\n");
  printf("   --mask-inv : invert mask\n");
  printf("   --X x.mat : matlab4 detrending matrix\n");
  printf("   --detrend order : polynomial detrending\n");
  printf("   --sum sumfile\n");
  printf("   \n");
  printf("   --fwhm fwhm : apply before measuring\n");
  printf("   --o output\n");
  printf("\n");
  printf("   --synth \n");
  printf("   --synth-frames nframes : default is 10 \n");
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
  if(inpath == NULL && !synth){
    printf("ERROR: need to specify --in or --synth\n");
    exit(1);
  }
  if(Xfile && DetrendOrder >= 0){
    printf("ERROR: cannot use --X and --detrend\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"%s\n",Progname);
  fprintf(fp,"FREESURFER_HOME %s\n",getenv("FREESURFER_HOME"));
  fprintf(fp,"cwd       %s\n",cwd);
  fprintf(fp,"cmdline   %s\n",cmdline);
  fprintf(fp,"timestamp %s\n",VERcurTimeStamp());
  fprintf(fp,"sysname   %s\n",uts.sysname);
  fprintf(fp,"hostname  %s\n",uts.nodename);
  fprintf(fp,"machine   %s\n",uts.machine);
  fprintf(fp,"user      %s\n",VERuser());
  if(inpath)   fprintf(fp,"input  %s\n",inpath);
  if(maskpath) fprintf(fp,"mask  %s\n",maskpath);
  fprintf(fp,"infwhm %lf\n",infwhm);
  fprintf(fp,"ingstd %lf\n",ingstd);
  fprintf(fp,"synth %d\n",synth);
  if(synth){
    fprintf(fp,"synth-frames %d\n",nframes);
    fprintf(fp,"seed  %d\n",SynthSeed);
  }
  if(Xfile) fprintf(fp,"xfile %s\n",Xfile);
  if(DetrendOrder >= 0) fprintf(fp,"detrend %d\n",DetrendOrder);
  if(sumfile) fprintf(fp,"sumfile  %s\n",sumfile);
  if(outpath) fprintf(fp,"out  %s\n",outpath);

  return;
}
