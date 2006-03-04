/*
BEGINHELP

Estimates the smoothness of a volume-based data set. The smoothness is
measured as the Full-Width-Half-Max (FWHM) in mm and resel size.

--i inputvol

Input data. Format must be something readable by mri_convert
(eg, mgh, mgz, img, nii). Alternately, one can synthesize
white gaussian noise with --synth and --synth-frames.

--mask maskfile

Compute FWHM only over voxels in the given mask. Format can be anything
accepted by mri_convert.

--mask-thresh thresh

Threshold mask at thresh. Default is 0.5 (ie, it expects a binary mask).

--auto-mask rthresh

Compute a mask based on a fraction (rthresh) of the global mean. If 
rthresh is 0.1, then all voxels in the mean image above 0.1*globalmean
are in the mask.

--mask-inv

Invert mask, ie, compute FWHM only over voxels outside the given mask.

--out-mask outmaskfile

Save final mask to outmaskfile. This mask does not include any eroding.

--nerode n

Erode mask n times prior to computing the fwhm. This applies only
to the fwhm computation and not the smoothing. This assures that
the fwhm is not computed near the edges of the mask.

--X x.mat

Detrend data with the matrix in x.mat. Ie, y = (I-inv(X'*X)*X')*y, where
y is the input. x.mat must be a matlab4 matrix. Not with --detrend.

--detrend order

Detrend data with polynomial of given order. Not with --X.

--sum sumfile

Prints ascii summary to sumfile.

--fwhm fwhm

Smooth by fwhm mm before estimating the fwhm. This is mainly good for 
debuggging. But with --out can also be used to smooth data.

--o outvol

Save (possibly synthesized and/or smoothed and/or detrended and/or masked) data 
to outfile. Automatically detects format. Format must be one accepted as by 
mri_convert. 

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

MRI *MRImaskedGaussianSmoothTo(MRI *invol, MRI *mask, double ToFWHM, 
			       double *pByFWHM, double tol, double *pToFWHMActual, 
			       MRI *outvol);

int getybest(double xa, double ya, double xb, double yb, double xc, double yc,
	     double *xbest, double *ybest, double ytarg);
double EvalFWHM(MRI *vol, MRI *mask);

MRI * MRIbinarize2(MRI *mri_src, MRI *mri_dst, 
		   double threshold, double low_val, double hi_val);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_fwhm.c,v 1.5 2006/03/04 19:31:16 greve Exp $";
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
double byfwhm, tofwhm, tofwhmact;
int synth = 0, nframes = 10;
int SynthSeed = -1;

char *Xfile=NULL;
MATRIX *X=NULL;
int DetrendOrder = -1;
int automask = 0;
double automaskthresh = .1;
int nerode = 0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int nargs, n, Ntp, nsearch, nsearch2=0;
  double fwhm = 0, nresels, voxelvolume, nvoxperresel, reselvolume;
  double car1mn, rar1mn,sar1mn,cfwhm,rfwhm,sfwhm, ftmp; 
  double gmean, gstd, gmax;
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

  // ------------- load or synthesize input ---------------------
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

  // -------------------- handle masking ------------------------
  if(maskpath){
    printf("Loading mask %s\n",maskpath);
    mask = MRIread(maskpath);
    if(mask==NULL) exit(1);
    MRIbinarize2(mask, mask, maskthresh, 0, 1);
  }
  if(automask){
    RFglobalStats(InVals, NULL, &gmean, &gstd, &gmax);
    maskthresh = gmean * automaskthresh;
    printf("Computing mask, relative threshold = %g, gmean = %g, absthresh = %g\n",
	   automaskthresh,gmean,maskthresh);
    mritmp = MRIframeMean(InVals,NULL);
    //MRIwrite(mritmp,"fmean.mgh");
    mask = MRIbinarize2(mritmp, NULL, maskthresh, 0, 1);
    MRIfree(&mritmp);
  }
  if(mask){
    if(maskinv){
      printf("Inverting mask\n");
      MRImaskInvert(mask,mask);
    }
    nsearch = MRInMask(mask);
    if(nsearch == 0){
      printf("ERROR: no voxels found in mask\n");
      exit(1);
    }
    //---- Save mask -----
    if(outmaskpath) MRIwrite(mask,outmaskpath);
  }
  else nsearch = InVals->width * InVals->height * InVals->depth;
  printf("Search region is %d voxels = %lf mm3\n",nsearch,nsearch*voxelvolume);


  //------------------------ Detrend ------------------
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

  // ------------ Smooth Input ---------------------------------
  if(infwhm > 0){
    printf("Smoothing input by fwhm=%lf, gstd=%lf\n",infwhm,ingstd);
    MRImaskedGaussianSmooth(InVals, mask, ingstd, InVals);
  }
  if(tofwhm > 0){
    printf("Attempting to smooth to %g fwhm\n",tofwhm);
    mritmp = MRImaskedGaussianSmoothTo(InVals, mask, tofwhm, &byfwhm, .5, &tofwhmact, InVals);
    if(mritmp == NULL) exit(1);
    printf("Smoothed by %g to %g \n",byfwhm,tofwhmact);
  }

  if(nerode > 0){
    printf("Eroding mask %d times\n",nerode);
    for(n=0; n<nerode; n++) MRIerode(mask,mask);
    nsearch2 = MRInMask(mask);
    if(nsearch2 == 0){
      printf("ERROR: no voxels found in mask after eroding\n");
      exit(1);
    }
    printf("%d voxels in mask after eroding\n",nsearch2);
  }

  // ----------- Compute smoothness -----------------------------
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
  printf("outfwhm = %lf\n",fwhm);

  reselvolume = cfwhm*rfwhm*sfwhm;
  nvoxperresel = reselvolume/voxelvolume;
  nresels = voxelvolume*nsearch/reselvolume;
  printf("reselvolume %lf\n",reselvolume);
  printf("nresels %lf\n",nresels);
  printf("nvoxperresel %lf\n",nvoxperresel);

  fflush(stdout);

  // ---------- Save summary file ---------------------
  if(sumfile){
    fp = fopen(sumfile,"w");
    if(fp == NULL){
      printf("ERROR: opening %s\n",sumfile);
      exit(1);
    }
    dump_options(fp);
    fprintf(fp,"nsearch2        %d\n",nsearch2);
    fprintf(fp,"searchspace_vox %d\n",nsearch);
    fprintf(fp,"searchspace_mm3 %lf\n",nsearch*voxelvolume);
    fprintf(fp,"voxelvolume_mm3 %g\n",voxelvolume);
    fprintf(fp,"voxelsize_mm %g %g %g\n",InVals->xsize,InVals->ysize,InVals->zsize);
    fprintf(fp,"ar1mn  %lf %lf %lf\n",car1mn,rar1mn,sar1mn);
    fprintf(fp,"colfwhm_mm      %lf\n",cfwhm);
    fprintf(fp,"rowfwhm_mm      %lf\n",rfwhm);
    fprintf(fp,"slicefwhm_mm    %lf\n",sfwhm);
    fprintf(fp,"outfwhm_mm      %lf\n",fwhm);
    fprintf(fp,"reselvolume_mm3 %lf\n",reselvolume);
    fprintf(fp,"nresels         %lf\n",nresels);
    fprintf(fp,"nvox_per_resel  %lf\n",nvoxperresel);
    fclose(fp);
  }

  if(outpath) MRIwrite(InVals,outpath);

  return 0;
}
/* ------------------=================--------------------------- */
/* ------------------=================--------------------------- */
/* ------------------=================--------------------------- */


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
    else if (!strcasecmp(option, "--auto-mask")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&automaskthresh);
      automask = 1;
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
    else if (!strcasecmp(option, "--tofwhm")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&tofwhm);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--nerode")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nerode);
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
  printf("\n");
  printf("   --mask maskvol : binary mask\n");
  printf("   --mask-thresh absthresh : threshold for mask (default is .5)\n");
  printf("   --auto-mask rthresh : compute mask\n");
  printf("   --mask-inv : invert mask\n");
  printf("   --nerode n : erode mask n times prior to computing fwhm\n");
  printf("   \n");
  printf("   --X x.mat : matlab4 detrending matrix\n");
  printf("   --detrend order : polynomial detrending\n");
  printf("   \n");
  printf("   --sum sumfile\n");
  printf("   \n");
  printf("   --fwhm fwhm : apply before measuring\n");
  printf("   \n");
  printf("   --o output : save input after detrending, masking, smoothing\n");
  printf("   --out-mask outmask : save final mask\n");
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
printf("\n");
printf("Estimates the smoothness of a volume-based data set. The smoothness is\n");
printf("measured as the Full-Width-Half-Max (FWHM) in mm and resel size.\n");
printf("\n");
printf("--i inputvol\n");
printf("\n");
printf("Input data. Format must be something readable by mri_convert\n");
printf("(eg, mgh, mgz, img, nii). Alternately, one can synthesize\n");
printf("white gaussian noise with --synth and --synth-frames.\n");
printf("\n");
printf("--mask maskfile\n");
printf("\n");
printf("Compute FWHM only over voxels in the given mask. Format can be anything\n");
printf("accepted by mri_convert.\n");
printf("\n");
printf("--mask-thresh thresh\n");
printf("\n");
printf("Threshold mask at thresh. Default is 0.5 (ie, it expects a binary mask).\n");
printf("\n");
printf("--auto-mask rthresh\n");
printf("\n");
printf("Compute a mask based on a fraction (rthresh) of the global mean. If \n");
printf("rthresh is 0.1, then all voxels in the mean image above 0.1*globalmean\n");
printf("are in the mask.\n");
printf("\n");
printf("--mask-inv\n");
printf("\n");
printf("Invert mask, ie, compute FWHM only over voxels outside the given mask.\n");
printf("\n");
printf("--out-mask outmaskfile\n");
printf("\n");
printf("Save final mask to outmaskfile. This mask does not include any eroding.\n");
printf("\n");
printf("--nerode n\n");
printf("\n");
printf("Erode mask n times prior to computing the fwhm. This applies only\n");
printf("to the fwhm computation and not the smoothing. This assures that\n");
printf("the fwhm is not computed near the edges of the mask.\n");
printf("\n");
printf("--X x.mat\n");
printf("\n");
printf("Detrend data with the matrix in x.mat. Ie, y = (I-inv(X'*X)*X')*y, where\n");
printf("y is the input. x.mat must be a matlab4 matrix. Not with --detrend.\n");
printf("\n");
printf("--detrend order\n");
printf("\n");
printf("Detrend data with polynomial of given order. Not with --X.\n");
printf("\n");
printf("--sum sumfile\n");
printf("\n");
printf("Prints ascii summary to sumfile.\n");
printf("\n");
printf("--fwhm fwhm\n");
printf("\n");
printf("Smooth by fwhm mm before estimating the fwhm. This is mainly good for \n");
printf("debuggging. But with --out can also be used to smooth data.\n");
printf("\n");
printf("--o outvol\n");
printf("\n");
printf("Save (possibly synthesized and/or smoothed and/or detrended and/or masked) data \n");
printf("to outfile. Automatically detects format. Format must be one accepted as by \n");
printf("mri_convert. \n");
printf("\n");
printf("--synth \n");
printf("\n");
printf("Synthesize input with white gaussian noise. Ten frames are used by default,\n");
printf("but this can be changed with --synth-frames.\n");
printf("\n");
printf("--synth-frames nframes\n");
printf("\n");
printf("Synthesize input with white gaussian noise with the given number of frames.\n");
printf("Implies --synth.\n");
printf("\n");
printf("\n");
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
  if(maskpath && automask){
    printf("ERROR: cannot use --mask and --auto-mask\n");
    exit(1);
  }
  if(outmaskpath && maskpath==NULL && !automask){
    printf("ERROR: cannot use --outmask without --mask or --auto-mask\n");
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
  fprintf(fp,"input  %s\n",inpath);
  fprintf(fp,"synth %d\n",synth);
  if(Xfile) fprintf(fp,"xfile %s\n",Xfile);
  if(DetrendOrder >= 0) fprintf(fp,"detrend %d\n",DetrendOrder);
  fprintf(fp,"infwhm %lf\n",infwhm);
  fprintf(fp,"ingstd %lf\n",ingstd);
  if(maskpath) fprintf(fp,"mask  %s\n",maskpath);
  if(automask) fprintf(fp,"automaskthresh %g\n",automaskthresh);
  if(mask) {
    fprintf(fp,"maskthresh %g\n",maskthresh);
    fprintf(fp,"maskinv %d\n",maskinv);
    fprintf(fp,"nerode  %d\n",nerode);
  }
  if(outmaskpath) fprintf(fp,"outmask  %s\n",outmaskpath);
  if(synth){
    fprintf(fp,"synth-frames %d\n",nframes);
    fprintf(fp,"seed  %d\n",SynthSeed);
  }
  if(outpath) fprintf(fp,"out  %s\n",outpath);

  return;
}

/*-------------------------------------------------------------------------------
  MRIbinarize2() - same as MRIbinarize() but passes theshold, low, and hi as
  doubles instead of UCHARs.
  -------------------------------------------------------------------------------*/
MRI * MRIbinarize2(MRI *mri_src, MRI *mri_dst, 
		   double threshold, double low_val, double hi_val)
{
  int     width, height, depth, x, y, z, f ;
  double  val;

  if(!mri_dst)  mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (f = 0 ; f < mri_src->nframes ; f++){
    for (z = 0 ; z < depth ; z++){
      for (y = 0 ; y < height ; y++) {
	for (x = 0 ; x < width ; x++){
	  val = MRIgetVoxVal(mri_src, x, y, z, f);
	  if(val > threshold) val = hi_val ;
	  else                val = low_val ;
	  MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
	}
      }
    }
  }
  
  return(mri_dst) ;
}


/*------------------------------------------------------------------------*/
double EvalFWHM(MRI *vol, MRI *mask)
{
  double car1mn, rar1mn, sar1mn;
  double cfwhm,rfwhm,sfwhm,fwhm;
  fMRIspatialAR1Mean(vol, mask, &car1mn, &rar1mn, &sar1mn);
  cfwhm = RFar1ToFWHM(car1mn, vol->xsize);
  rfwhm = RFar1ToFWHM(rar1mn, vol->ysize);
  sfwhm = RFar1ToFWHM(sar1mn, vol->zsize);
  fwhm = sqrt((cfwhm*cfwhm + rfwhm*rfwhm + sfwhm*sfwhm)/3.0);
  return(fwhm);
}
/*------------------------------------------------------------------------*/
MRI *MRImaskedGaussianSmoothTo(MRI *invol, MRI *mask, double ToFWHM, 
			       double *pByFWHM, double tol, double *pToFWHMActual, 
			       MRI *outvol)
{
  double SrcFWHM, ByGStd;
  MRI *volsm;
  double ya, yb, yc, xa, xb, xc, s1, s2;
  double xn,yn;
  double C, R, err;
  int nth;

  C = (3.0-sqrt(5.0))/2.0; // golden mean
  R = 1-C;

  SrcFWHM = EvalFWHM(invol, mask);
  if(SrcFWHM > ToFWHM){
    printf("ERROR: MRImaskedGaussianSmoothTo(): source fwhm = %g > to fwhm = %g\n",
	   SrcFWHM,ToFWHM);
    return(NULL);
  }
  xa = 0;
  ya = SrcFWHM;
  printf("Unsmoothed actual fwhm of %g\n",ya);

  // Check whether we are close enough already
  if(fabs(SrcFWHM-ToFWHM) < tol){
    *pByFWHM = 0.0;
    *pToFWHMActual = SrcFWHM;
    outvol = MRIcopy(invol,outvol);
    return(outvol);
  }

  volsm = MRIcopy(invol,NULL); // allocate

  // First point in the bracket

  // Second point in the bracket
  xb = sqrt(ToFWHM*ToFWHM - SrcFWHM*SrcFWHM); // power law
  ByGStd = xb/sqrt(log(256.0));
  printf("Trying smoothing by fwhm %g  ",xb); fflush(stdout);
  MRImaskedGaussianSmooth(invol, mask, ByGStd, volsm);
  yb = EvalFWHM(volsm, mask);
  printf("results in actual fwhm of %g\n",yb);
  // Check whether we are close enough now
  if(fabs(yb-ToFWHM) < tol){
    *pByFWHM = xb;
    *pToFWHMActual = yb;
    outvol = MRIcopy(volsm,outvol);
    MRIfree(&volsm);
    return(outvol);
  }

  // Third point in the bracket. Not sure how to choose this point, 
  // it needs to be far enough to bracket the min, but too far
  // and we end up doing to many evaluations.
  xc = xb + (xb - xa); // 
  ByGStd = xc/sqrt(log(256.0));
  printf("Trying smoothing by fwhm %g  ",xc); fflush(stdout);
  MRImaskedGaussianSmooth(invol, mask, ByGStd, volsm);
  yc = EvalFWHM(volsm, mask);
  printf("results in actual fwhm of %g\n",yc);
  // Check whether we are close enough now
  if(fabs(yc-ToFWHM) < tol){
    *pByFWHM = xc;
    *pToFWHMActual = yc;
    outvol = MRIcopy(volsm,outvol);
    MRIfree(&volsm);
    return(outvol);
  }
  if(yc < ToFWHM){
    // Did not step far enough out
    printf("ERROR: did not step far enough out\n");
    return(NULL);
  }

  // ok, we've brackated the min, now chase it down like a scared rabbit
  nth = 0;
  getybest(xa, ya, xb, yb, xc, yc, pByFWHM, pToFWHMActual, ToFWHM);
  err = fabs(*pToFWHMActual-ToFWHM);
  while( err > tol){
    if(nth > 20){
      printf("ERROR: timed out at ntry=%d\n",nth);
      MRIfree(&volsm);
      return(NULL);
    }
    printf("ntry=%3d  by=%6.4lf  to=%6.4lf (%6.4lf,%6.4lf)  err=%g\n",nth, 
	   *pByFWHM, *pToFWHMActual, xa, xc, err);
    s1 = xb - xa;
    s2 = xc - xb;
    if(s1 > s2){ // s1 is bigger
      xn = xa + R*s1;
      printf("   Trying smoothing by fwhm %g  ",xn); fflush(stdout);
      ByGStd = xn/sqrt(log(256.0));
      MRImaskedGaussianSmooth(invol, mask, ByGStd, volsm);
      yn = EvalFWHM(volsm, mask);
      printf("results in actual fwhm of %g\n",yn);
      if(yn < yb){
	xc = xb;
	yc = yb;
	xb = xn;
	yb = yn;
      } else {
	xa = xn; //a replaced by new
	ya = yn;
	// b and c stay the same
      }
    } else { // s2 is bigger
      xn = xb + C*s2;
      ByGStd = xn/sqrt(log(256.0));
      printf("   Trying smoothing by fwhm %g  ",xn); fflush(stdout);
      MRImaskedGaussianSmooth(invol, mask, ByGStd, volsm);
      yn = EvalFWHM(volsm, mask);
      printf("results in actual fwhm of %g\n",yn);
      if(yn < yb){
	xa = xb; //a replaced by  b
	ya = yb;
	xb = xn; //b replace by new
	yb = yn;
	// c stays the same
      } else {
	xc = xn; //c replaced by new
	yc = yn;
	// a and b stay the same
      }
    }
    getybest(xa, ya, xb, yb, xc, yc, pByFWHM, pToFWHMActual, ToFWHM);
    err = fabs(*pToFWHMActual-ToFWHM);
    nth++;
  }
  printf("ntry=%3d  by=%6.4lf  to=%6.4lf targ=%6.4lf  err=%g\n",nth, *pByFWHM, *pToFWHMActual, ToFWHM, err);
  outvol = MRIcopy(volsm,outvol);
  MRIfree(&volsm);
  return(outvol);
}


int getybest(double xa, double ya, double xb, double yb, double xc, double yc,
	     double *xbest, double *ybest, double ytarg)
{
  double ea, eb, ec;

  ea = fabs(ya-ytarg);
  eb = fabs(yb-ytarg);
  ec = fabs(yc-ytarg);

  if(ea < eb && ea < ec){
    *xbest = xa;
    *ybest = ya;
    return(0);
  }
  if(eb < ec){
    *xbest = xb;
    *ybest = yb;
    return(0);
  }
  *xbest = xc;
  *ybest = yc;
  return(0);
}
