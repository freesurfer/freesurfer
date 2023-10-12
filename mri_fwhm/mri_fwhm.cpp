/*
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
BEGINHELP

FreeSurfer program to estimate the global Gaussian smoothness of a
multi-frame, volume-based data set. The smoothness is measured as the
Full-Width-Half-Max (FWHM) in mm.  Gaussian stddev =
fwhm/sqrt(log(256.0)). The voxels used in the fwhm estimation can be
constrained to be inside of a mask. It is STRONGLY recommended that a
masked be used. The mask can be specified explictly or computed
automatically. By default, the time course will be detrended by
removing the mean. Higher order polynomial detrending is
possible. Alternatively, the user can specify a detrending matrix. The
data can be smoothed BY a given fwhm or TO a given fwhm prior to
estimating the fwhm.  The resulting data can then be saved (thus
turning this program into a smoother). If smoothing is to be done, it
will only be done inside the mask (except see --save-unmasked).

For more info, see surfer.nmr.mgh.harvard.edu. For help, send email
and summary file to  freesurfer@nmr.mgh.harvard.edu.

--i inputvol

Input data. Format must be something readable by mri_convert
(eg, mgh, mgz, img, nii, nii.gz). Alternately, one can synthesize
white gaussian noise with --synth and --synth-frames in which
case inputvol is used as a dimension template.

--o outputvol

Save input after smoothing. See also --save-detended and --save-unmasked.

--smooth-only

Does not attempt to compute FWHM. Smooths the input, saves to outputvol,
and exists. Respects --save-unmasked, but not --save-detended. This allows
for data sets with fewer than 10 frames to be smoothed.

--save-detended

Save input after smoothing, masking, and detrending.

--save-unmasked

Save input after smoothing and detrending, but do not mask while
smoothing. Note that the output will be different even inside
the mask because the smoother handles voxels on the boundary
of the mask differently than those at the center.

--mask maskvol

Compute FWHM only over voxels in the given mask. Format can be anything
accepted by mri_convert. If smoothing is to be done, it will only be
done inside the mask. It is strongly recommended that a masked be used
(see also --auto-mask and --save-unmasked).

--mask-thresh thresh

Threshold mask at thresh. Default is 0.5 (ie, it expects a binary mask).

--auto-mask rthresh

Compute a mask based on a fraction (rthresh) of the global mean. If
rthresh is 0.1, then all voxels in the mean image above 0.1*globalmean
are in the mask.

--mask-inv

Invert mask, ie, compute FWHM only over voxels outside the given mask.

--nerode n

Erode mask n times (ie, make it smaller). Occurs after any mask inversion.

--out-mask outmaskvol

Save final mask to outmaskvol.

--ar1 ar1path

Save spatial AR1 volume (6 frames)

--ar1red ar1redpath

Save spatial AR1 volume (3 frames), averaging the 6

--X x.mat

Detrend/residualize data with the matrix in x.mat. Ie,
y = (I-inv(X'*X)*X')*y, where y is the input. x.mat must
be a matlab4 matrix. Not with --detrend.

--detrend order

Detrend data with polynomial of given order. Not with --X. Note:
if neither --X nor --detrend are specified, then detrending
order of 0 is used (ie, the mean is removed).

--fwhm fwhm

Smooth BY fwhm mm before estimating the fwhm. This is mainly good for
debuggging. But with --out can also be used to smooth data.

--to-fwhm tofwhm

Smooth TO tofwhm mm. This is idea proposed by Lee Friedman
(lfriedman10@comcast.net) and fBIRN (www.birn.net) as a way to
reduce variation across data sets, particularly when the data
may have come from different scanners. The method implemented
here uses an iterative approach in which the data are smoothed
BY varying amounmts until the resulting fwhm is within a certain
tolerance of tofwhm. By default, the tolerance is 0.5 mm, but
can be changed with --to-fwhm-tol. It will iterate at most 20
times (can be changed with --to-fwhm-nmax). An error will be
returned if the tofwhm is less than the inherent fwhm or the max
number of iterations is exceeded. The minimization is done with
a Golden Section Search (see Numerical Recipes in C).

--to-fwhm-tol tol

Keep iterating the tofwhm search until the result is within tol
of the desired fwhm (or until the maximum number of iterations
is reached).

--to-fwhm-nmax niterationsmax

Maximum number of iterations. Default is 20.

--to-fwhm-file file

Save some results of the tofwhm minimization to file. Good for
debugging. Results also saved in summary file.

--sum sumfile

Prints summary to ascii sumfile. Send this file when requesting
help or more information.

--dat datfile

Prints only the final fwhm estimate into this file.

--synth

Synthesize input with white gaussian noise. Ten frames are used by default,
but this can be changed with --synth-frames. Uses input volume as template.
This functionality is useful for degugging. Eg, when using --synth and --fwhm,
it should measure the resulting fwhm to be that passed by --fwhm. Can do
the same with --to-fwhm.

--synth-frames nframes

Synthesize input with white gaussian noise with the given number of frames.
Implies --synth.

--tr TRms

Set TR in msec (generally not too useful)

--nthreads nthreads

Set OPEN MP threads

--inorm

Spatial intensity normalization. Subtract the in-mask mean and divide by the in-mask 
stddev. 

EXAMPLES:

1. Measure the fwhm of an input data set, compute mask automatically by
   thresholding input at 20%% of global mean. The time series will be
   have its mean removed prior to computing the fwhm. Save result in
   a summary file (one example uses mgh format, the other gzipped NIFTI):

      mri_fwhm --i f.mgh    --auto-mask .2 --sum f.fwhm.sum
      mri_fwhm --i f.nii.gz --auto-mask .2 --sum f.fwhm.sum

2. Same as above, but smooth input BY 5mm fwhm first. Save the
   smoothed output in fsm5.mgh. Save the mask to automask.nii.
   Note: mask is computed on unsmoothed data.

      mri_fwhm --i f.mgh --auto-mask .2 --sum f.fwhm5.sum
        --fwhm 5 --o fsm5.mgh --out-mask automask.nii

3. Same as above, but smooth input TO 5 +/- .1mm fwhm first.
   Save the smoothed output in fto5.mgh.

      mri_fwhm --i f.mgh --auto-mask .2 --sum f.fwhm5.sum
        --to-fwhm-tol .1 --to-fwhm 5 --o fto5.mgh

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
#include "mrinorm.h"

#include "romp_support.h"


MRI *MRImaskedGaussianSmoothTo(MRI *invol, MRI *mask, double ToFWHM,
                               double tol, int nitersmax,
                               double *pByFWHM, double *pToFWHMActual, int *niters,
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

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *inpath=NULL;
char *outpath=NULL;
char *sumfile=NULL;
char *datfile=NULL;
MRI *InVals=NULL;
MRI *InValsCopy=NULL;
int InValsType = MRI_VOLUME_TYPE_UNKNOWN;

char *maskpath=NULL;
MRI *mask=NULL;
int maskinv = 0;
double maskthresh = 0.5;
char *outmaskpath=NULL;

MRI *mritmp=NULL;
char tmpstr[2000];

MRI *ar1;
char *ar1path=NULL;

double infwhm = 0,ingstd =0;
double infwhmc=0, infwhmr=0, infwhms=0;
double ingstdc=0, ingstdr=0, ingstds=0;
double byfwhm;
double bygstd;
char *tofwhmfile = NULL;
double  tofwhm, togstd, tofwhmact, tofwhmtol=0.5;
int tofwhmnitersmax = 20;
int tofwhmniters;
int synth = 0, nframes = -1;
int SynthSeed = -1;

char *Xfile=NULL;
MATRIX *X=NULL;
int DetrendOrder = -1;
int SaveDetrended = 0;
int SaveUnmasked = 0;
int automask = 0;
double automaskthresh = .1;
int nerode = 0;
int SmoothOnly = 0;
int nframesmin = 10;
int DoSqr = 0; // take square of input before smoothing
int DoMedian = 0, MedianWidth=0;

char *sum2file = NULL;
char *arNfname = NULL;
int arNlags;

int DoAR2;

double TR=0.0;
int SetTR=0;

MB2D *mb2drad=NULL,*mb2dtan=NULL;
int DoSpatialINorm = 0;
char *ar1redfile=NULL, *fwhmvolfile=NULL;
char *fwhmvolmnfile=NULL,*fwhmdatfile=NULL;

int fMRIspatialFWHMMean(MRI *fhwmvol, MRI *mask, double *cfwhmmn, double *rfwhmmn, double *sfwhmmn);

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, n, Ntp, nsearch, nsearch2=0;
  double fwhm = 0, nresels, voxelvolume, nvoxperresel, reselvolume;
  double car1mn, rar1mn,sar1mn,cfwhm,rfwhm,sfwhm, ftmp;
  double car2mn, rar2mn,sar2mn;
  double gmean, gstd, gmax;
  FILE *fp;

  sprintf(tmpstr, "S%sER%sRONT%sOR", "URF", "_F", "DO") ;
  setenv(tmpstr,"1",0);

  nargs = handleVersionOption(argc, argv, "mri_fwhm");
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

  if (SynthSeed < 0) SynthSeed = PDFtodSeed();
  if (debug) dump_options(stdout);

  // ------------- load or synthesize input ---------------------
  InVals = MRIreadType(inpath,InValsType);
  if(InVals == NULL) exit(1);
  if(SetTR){
    printf("Setting TR to %g ms\n",TR);
    InVals->tr = TR;
  }
  if((nframes < 0 && synth) || !synth) nframes = InVals->nframes;
  if(nframes < nframesmin && !SmoothOnly && !sum2file) {
    printf("ERROR: nframes = %d, need at least %d\n",
	   nframes,nframesmin);
    exit(1);
  }
  if (InVals->type != MRI_FLOAT) {
    mritmp = MRISeqchangeType(InVals, MRI_FLOAT, 0, 0, 0);
    MRIfree(&InVals);
    InVals = mritmp;
  }
  if(synth) {
    printf("Synthesizing %d frames, Seed = %d\n",nframes,SynthSeed);
    mritmp = MRIcloneBySpace(InVals,MRI_FLOAT,nframes);
    MRIfree(&InVals);
    MRIrandn(mritmp->width, mritmp->height, mritmp->depth, 
	     nframes, 0, 1, mritmp);
    InVals = mritmp;
  }
  voxelvolume = InVals->xsize * InVals->ysize * InVals->zsize ;
  printf("voxelvolume %g mm3\n",voxelvolume);

  if(DoSqr){
    printf("Computing square of input\n");
    MRIsquare(InVals,NULL,InVals);
  }

  // -------------------- handle masking ------------------------
  if (maskpath) {
    printf("Loading mask %s\n",maskpath);
    mask = MRIread(maskpath);
    if(mask==NULL) exit(1);
    if(MRIdimMismatch(mask,InVals,0)){
      printf("ERROR: dimension mismatch between mask and input\n");
      exit(1);
    }
    MRIbinarize2(mask, mask, maskthresh, 0, 1);
  }
  if (automask) {
    RFglobalStats(InVals, NULL, &gmean, &gstd, &gmax);
    maskthresh = gmean * automaskthresh;
    printf("Computing mask, relative threshold = %g, gmean = %g, absthresh = %g\n",
           automaskthresh,gmean,maskthresh);
    mritmp = MRIframeMean(InVals,NULL);
    //MRIwrite(mritmp,"fmean.mgh");
    mask = MRIbinarize2(mritmp, NULL, maskthresh, 0, 1);
    MRIfree(&mritmp);
  }
  if (mask) {
    if (maskinv) {
      printf("Inverting mask\n");
      MRImaskInvert(mask,mask);
    }
    nsearch = MRInMask(mask);
    if (nsearch == 0) {
      printf("ERROR: no voxels found in mask\n");
      exit(1);
    }
    // Erode the mask -----------------------------------------------
    if (nerode > 0) {
      printf("Eroding mask %d times\n",nerode);
      for (n=0; n<nerode; n++) MRIerode(mask,mask);
      nsearch2 = MRInMask(mask);
      if (nsearch2 == 0) {
        printf("ERROR: no voxels found in mask after eroding\n");
        exit(1);
      }
      printf("%d voxels in mask after eroding\n",nsearch2);
    }
    //---- Save mask -----
    if (outmaskpath) MRIwrite(mask,outmaskpath);
  } else nsearch = InVals->width * InVals->height * InVals->depth;
  printf("Search region is %d voxels = %lf mm3\n",nsearch,nsearch*voxelvolume);

  if(DoMedian == 0){
    if( (infwhm > 0 || infwhmc > 0 || infwhmr > 0 || infwhms > 0 || mb2drad || mb2dtan) && SmoothOnly) {
      if(SaveUnmasked) mritmp = NULL;
      else             mritmp = mask;
      if(infwhm > 0) {
	printf("Smoothing input by fwhm=%lf, gstd=%lf\n",infwhm,ingstd);
	MRImaskedGaussianSmooth(InVals, mritmp, ingstd, InVals);
      }
      if(infwhmc > 0 || infwhmr > 0 || infwhms > 0) {
	printf("Smoothing input by fwhm=(%lf,%lf,%lf) gstd=(%lf,%lf,%lf)\n",
	       infwhmc,infwhmr,infwhms,ingstdc,ingstdr,ingstds);
	MRIgaussianSmoothNI(InVals, ingstdc, ingstdr, ingstds, InVals);
      }
      if(mb2drad){
	printf("Applying radial motion blur slope=%lf\n",mb2drad->slope);
	mb2drad->DeltaD = InVals->xsize/2.0;
	mb2drad->c0 = InVals->width/2.0; // center of volume
	mb2drad->r0 = InVals->height/2.0; // center of volume
	mritmp = MRImotionBlur2D(InVals, mb2drad, NULL);
	MRIfree(&InVals);
	InVals = mritmp;
      }
      if(mb2dtan){
	printf("Applying tangential motion blur slope=%lf\n",mb2dtan->slope);
	mb2dtan->DeltaD = InVals->xsize/2.0;
	mb2dtan->c0 = InVals->width/2.0; // center of volume
	mb2dtan->r0 = InVals->height/2.0; // center of volume
	mritmp = MRImotionBlur2D(InVals, mb2dtan, NULL);
	MRIfree(&InVals);
	InVals = mritmp;
      }
      if(DoSpatialINorm){
	mritmp = SpatialINorm(InVals, mask, NULL);
	MRIfree(&InVals);
	InVals = mritmp;
      }
      printf("Saving to %s\n",outpath);
      MRIwrite(InVals,outpath);
      printf("SmoothOnly requested, so exiting now\n");
      exit(0);
    }
  } else {
    printf("Running median filter %d\n",MedianWidth);
    mritmp = MRImedian(InVals, NULL, MedianWidth, NULL);
    MRIfree(&InVals);
    InVals = mritmp;
    if(SmoothOnly){
      printf("Saving to %s\n",outpath);
      MRIwrite(InVals,outpath);
      printf("SmoothOnly requested, so exiting now\n");
      exit(0);
    }
  }

  // Make a copy, if needed, prior to doing anything to data
  if(outpath) InValsCopy = MRIcopy(InVals,NULL);

  // Compute variance reduction factor -------------------
  if(sum2file){
    ftmp = MRIsum2All(InVals);
    fp = fopen(sum2file,"w");
    if(fp == NULL){
      printf("ERROR: opening %s\n",sum2file);
      exit(1);
    }
    printf("sum2all: %20.10lf\n",ftmp);
    printf("vrf: %20.10lf\n",1/ftmp);
    fprintf(fp,"%20.10lf\n",ftmp);
    exit(0);
  }

  //------------------------ Detrend ------------------
  if(DetrendOrder >= 0) {
    Ntp = InVals->nframes;
    printf("Polynomial detrending, order = %d\n",DetrendOrder);
    X = MatrixAlloc(Ntp,DetrendOrder+1,MATRIX_REAL);
    for (n=0;n<Ntp;n++) X->rptr[n+1][1] = 1.0;
    ftmp = Ntp/2.0;
    if (DetrendOrder >= 1)
      for (n=0;n<Ntp;n++) X->rptr[n+1][2] = (n-ftmp)/ftmp;
    if (DetrendOrder >= 2)
      for (n=0;n<Ntp;n++) X->rptr[n+1][3] = pow((n-ftmp),2.0)/(ftmp*ftmp);
  }
  if(X){
    printf("Detrending\n");
    if (X->rows != InVals->nframes) {
      printf("ERROR: dimension mismatch between X and input\n");
      exit(1);
    }
    mritmp = fMRIdetrend(InVals,X);
    if (mritmp == NULL) exit(1);
    MRIfree(&InVals);
    InVals = mritmp;
  }

  // ------------ Smooth Input isotropically BY infwhm -------------------------
  if(infwhm > 0) {
    printf("Smoothing input by fwhm=%lf, gstd=%lf\n",infwhm,ingstd);
    MRImaskedGaussianSmooth(InVals, mask, ingstd, InVals);
  }
  // ------------ Smooth Input nonisotropically BY infwhm -------------------------
  if(infwhmc > 0 || infwhmr > 0 || infwhms > 0) {
    printf("Smoothing input by fwhm=(%lf,%lf,%lf) gstd=(%lf,%lf,%lf)\n",
	   infwhmc,infwhmr,infwhms,ingstdc,ingstdr,ingstds);
    MRIgaussianSmoothNI(InVals, ingstdc, ingstdr, ingstds, InVals);
  }
  if(mb2drad){
    printf("Applying radial motion blur\n");
    mritmp = MRImotionBlur2D(InVals, mb2drad, NULL);
    MRIfree(&InVals);
    InVals = mritmp;
  }
  if(mb2dtan){
    printf("Applying tangential motion blur\n");
    mritmp = MRImotionBlur2D(InVals, mb2dtan, NULL);
    MRIfree(&InVals);
    InVals = mritmp;
  }

  // ------------ Smooth Input TO fwhm -------------------------
  if (tofwhm > 0) {
    printf("Attempting to smooth to %g +/- %g mm fwhm (nitersmax=%d)\n",
           tofwhm,tofwhmtol,tofwhmnitersmax);
    mritmp = MRImaskedGaussianSmoothTo(InVals, mask, tofwhm,
                                       tofwhmtol, tofwhmnitersmax,
                                       &byfwhm, &tofwhmact, &tofwhmniters,
                                       InVals);
    if (mritmp == NULL) exit(1);
    printf("Smoothed by %g to %g in %d iterations\n",
           byfwhm,tofwhmact,tofwhmniters);
    if (tofwhmfile) {
      fp = fopen(tofwhmfile,"w");
      if (!fp) {
        printf("ERROR: opening %s\n",tofwhmfile);
        exit(1);
      }
      fprintf(fp,"tofwhm    %lf\n",tofwhm);
      fprintf(fp,"tofwhmtol %lf\n",tofwhmtol);
      fprintf(fp,"tofwhmact %lf\n",tofwhmact);
      fprintf(fp,"byfwhm    %lf\n",byfwhm);
      fprintf(fp,"niters    %d\n",tofwhmniters);
      fprintf(fp,"nitersmax %d\n",tofwhmnitersmax);
      fclose(fp);
    }
  }

  // ------ Save smoothed/detrended ------------------------------
  if(outpath) {
    // This is a bit of a hack in order to be able to save undetrended
    // Operates on InValsCopy, which has not been modified (requires
    // smoothing twice, which is silly:).
    printf("Saving to %s\n",outpath);
    // Smoothed output will not be masked
    if (SaveDetrended && X) {
      mritmp = fMRIdetrend(InValsCopy,X);
      if (mritmp == NULL) exit(1);
      MRIfree(&InValsCopy);
      InValsCopy = mritmp;
    }
    if (SaveUnmasked) mritmp = NULL;
    else             mritmp = mask;
    if(infwhm > 0)
      MRImaskedGaussianSmooth(InValsCopy, mritmp, ingstd, InValsCopy);
    if(infwhmc > 0 || infwhmr > 0 || infwhms > 0) 
      MRIgaussianSmoothNI(InValsCopy, ingstdc, ingstdr, ingstds, InValsCopy);
    if(tofwhm > 0) {
      bygstd = byfwhm/sqrt(log(256.0));
      MRImaskedGaussianSmooth(InValsCopy, mritmp, bygstd, InValsCopy);
    }
    if(mb2drad){
      printf("Applying radial motion blur\n");
      mritmp = MRImotionBlur2D(InValsCopy, mb2drad, NULL);
      MRIfree(&InValsCopy);
      InValsCopy = mritmp;
    }
    if(mb2dtan){
      printf("Applying tangential motion blur\n");
      mritmp = MRImotionBlur2D(InValsCopy, mb2dtan, NULL);
      MRIfree(&InValsCopy);
      InValsCopy = mritmp;
    }

    MRIwrite(InValsCopy,outpath);
    MRIfree(&InValsCopy);
  }

  if(arNfname){
    printf("Computing spatial ARN %d in volume.\n",arNlags);
    ar1 = fMRIspatialARN(InVals, mask, arNlags, NULL);
    if(ar1 == NULL) exit(1);
    MRIwrite(ar1,arNfname);
  }

  // ----------- Compute smoothness -----------------------------
  printf("Computing spatial AR1 in volume.\n");
  ar1 = fMRIspatialAR1(InVals, mask, NULL);
  if (ar1 == NULL) exit(1);
  fMRIspatialAR1Mean(ar1, mask, &car1mn, &rar1mn, &sar1mn);

  cfwhm = RFar1ToFWHM(car1mn, InVals->xsize);
  rfwhm = RFar1ToFWHM(rar1mn, InVals->ysize);
  sfwhm = RFar1ToFWHM(sar1mn, InVals->zsize);
  fwhm = sqrt((cfwhm*cfwhm + rfwhm*rfwhm + sfwhm*sfwhm)/3.0);
  printf("ar1mn = (%lf,%lf,%lf)\n",car1mn,rar1mn,sar1mn);
  printf("colfwhm   = %lf\n",cfwhm);
  printf("rowfwhm   = %lf\n",rfwhm);
  printf("slicefwhm = %lf\n",sfwhm);
  printf("outfwhm = %lf\n",fwhm);
  if(fwhmdatfile){
    FILE *fp = fopen(fwhmdatfile,"w");
    fprintf(fp,"%8.4lf %8.4lf %8.4lf\n",cfwhm,rfwhm,sfwhm);
    fclose(fp);
  }

  reselvolume = cfwhm*rfwhm*sfwhm;
  nvoxperresel = reselvolume/voxelvolume;
  nresels = voxelvolume*nsearch/reselvolume;
  printf("reselvolume %lf\n",reselvolume);
  printf("nresels %lf\n",nresels);
  printf("nvoxperresel %lf\n",nvoxperresel);

  if(ar1redfile || fwhmvolfile || fwhmvolmnfile){
    // Reduce the 6 measures down to 3 by averaging
    MRI *ar1red = fMRIspatialARreduce(ar1,mask,NULL);
    if(ar1redfile){
      MRIwrite(ar1red,ar1redfile);
    }
    if(fwhmvolfile || fwhmvolmnfile ){
      // Compute fwhm at each voxel
      double fvoxsize[] = {ar1->xsize,ar1->ysize,ar1->zsize};
      MRI *fwhmvol = fMRIarToFWHM(ar1red, 1, fvoxsize, mask, NULL);
      if(fwhmvolfile) MRIwrite(fwhmvol,fwhmvolfile);
      double cfwhmmn, rfwhmmn, sfwhmmn;
      fMRIspatialFWHMMean(fwhmvol, mask, &cfwhmmn, &rfwhmmn, &sfwhmmn);
      printf("vcolfwhm   = %lf\n",cfwhmmn);
      printf("vrowfwhm   = %lf\n",rfwhmmn);
      printf("vslicefwhm = %lf\n",sfwhmmn);
      if(fwhmvolmnfile){
	FILE *fp = fopen(fwhmvolmnfile,"w");
	fprintf(fp,"%8.4lf %8.4lf %8.4lf\n",cfwhmmn,rfwhmmn,sfwhmmn);
	fclose(fp);
      }
    }
  }

  if(DoAR2){
    printf("Computing spatial AR2 in volume.\n");
    fMRIspatialAR2Mean(InVals, mask, &car2mn, &rar2mn, &sar2mn);
    printf("ar2mn = (%lf,%lf,%lf)\n",car2mn,rar2mn,sar2mn);
  }

  if(ar1path) MRIwrite(ar1,ar1path);

  fflush(stdout);

  // ---------- Save summary file ---------------------
  if(sumfile) {
    fp = fopen(sumfile,"w");
    if (fp == NULL) {
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

  if(datfile) {
    fp = fopen(datfile,"w");
    if(fp == NULL) {
      printf("ERROR: opening %s\n",datfile);
      exit(1);
    }
    fprintf(fp,"%lf\n",fwhm);
    fclose(fp);
  }


  printf("mri_fwhm done\n");

  return 0;
}
/* ------------------=================--------------------------- */
/* ------------------=================--------------------------- */
/* ------------------=================--------------------------- */


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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--synth")) synth = 1;
    else if (!strcasecmp(option, "--nosynth")) synth = 0;
    else if (!strcasecmp(option, "--mask-inv")) maskinv = 1;
    else if (!strcasecmp(option, "--save-detrended")) SaveDetrended = 1;
    else if (!strcasecmp(option, "--save-unmasked")) SaveUnmasked = 1;
    else if (!strcasecmp(option, "--smooth-only")) SmoothOnly = 1;
    else if (!strcasecmp(option, "--so")) SmoothOnly = 1;
    else if (!strcasecmp(option, "--sqr")) DoSqr = 1;
    else if (!strcasecmp(option, "--ispm")) InValsType = MRI_ANALYZE_FILE;
    else if (!strcasecmp(option, "--ar2")) DoAR2 = 1;
    else if (!strcasecmp(option, "--gdiag")) Gdiag_no = 1;
    else if (!strcasecmp(option, "--inorm")) DoSpatialINorm = 1;
    else if (!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      inpath = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskpath = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--mask-thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&maskthresh);
      nargsused = 1;
    } else if (!strcasecmp(option, "--auto-mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&automaskthresh);
      automask = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--out-mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      outmaskpath = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--sum")) {
      if (nargc < 1) CMDargNErr(option,1);
      sumfile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--dat")) {
      if (nargc < 1) CMDargNErr(option,1);
      datfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--fwhm")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&infwhm);
      ingstd = infwhm/sqrt(log(256.0));
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mb-rad")) {
      if (nargc < 2) CMDargNErr(option,2);
      mb2drad = (MB2D *) calloc(sizeof(MB2D),1);
      mb2drad->type = MB_RADIAL;
      mb2drad->cutoff = 4;// number of stddevs to cut off kernel
      mb2drad->Interp = SAMPLE_NEAREST;
      sscanf(pargv[0],"%lf",&mb2drad->offset);
      sscanf(pargv[1],"%lf",&mb2drad->slope);
      nargsused = 2;
    }
    else if (!strcasecmp(option, "--mb-tan")) {
      if (nargc < 2) CMDargNErr(option,2);
      mb2dtan = (MB2D *) calloc(sizeof(MB2D),1);
      mb2dtan->type = MB_TANGENTIAL;
      mb2dtan->cutoff = 4;
      mb2dtan->Interp = SAMPLE_NEAREST;
      sscanf(pargv[0],"%lf",&mb2dtan->offset);
      sscanf(pargv[1],"%lf",&mb2dtan->slope);
      nargsused = 2;
    }
    else if (!strcasecmp(option, "--median")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&MedianWidth);
      DoMedian = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--fwhmc")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&infwhmc);
      ingstdc = infwhmc/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcasecmp(option, "--fwhmr")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&infwhmr);
      ingstdr = infwhmr/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcasecmp(option, "--fwhms")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&infwhms);
      ingstds = infwhms/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcasecmp(option, "--gstd")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&ingstd);
      infwhm = ingstd*sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcasecmp(option, "--to-fwhm")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&tofwhm);
      nargsused = 1;
    } else if (!strcasecmp(option, "--to-gstd")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&togstd);
      tofwhm = togstd*sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcasecmp(option, "--to-fwhm-tol")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&tofwhmtol);
      nargsused = 1;
    } else if (!strcasecmp(option, "--to-fwhm-nmax")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&tofwhmnitersmax);
      nargsused = 1;
    } else if (!strcasecmp(option, "--to-fwhm-file")) {
      if (nargc < 1) CMDargNErr(option,1);
      tofwhmfile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--nerode")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nerode);
      nargsused = 1;
    } else if (!strcasecmp(option, "--nframesmin")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nframesmin);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--synth-frames")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nframes);
      synth = 1;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--g2")) {
      if (nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%f",&smni_cw1);
      if(smni_cw1 > 1){
	printf("ERROR: g2 w1 = %g > 1\n",smni_cw1);
	exit(1);
      }
      sscanf(pargv[1],"%f",&smni_cstd2); // Actually fwhm
      smni_cstd2 = smni_cstd2/sqrt(log(256.0));
      smni_rw1 = smni_cw1;
      smni_sw1 = smni_cw1;
      smni_rstd2 = smni_cstd2;
      smni_sstd2 = smni_cstd2;
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--g2c")) {
      if (nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%f",&smni_cw1);
      if(smni_cw1 > 1){
	printf("ERROR: g2 w1 = %g > 1\n",smni_cw1);
	exit(1);
      }
      sscanf(pargv[1],"%f",&smni_cstd2); // Actually fwhm
      smni_cstd2 = smni_cstd2/sqrt(log(256.0));
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--g2r")) {
      if (nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%f",&smni_rw1);
      if(smni_rw1 > 1){
	printf("ERROR: g2 w1 = %g > 1\n",smni_rw1);
	exit(1);
      }
      sscanf(pargv[1],"%f",&smni_rstd2); // Actually fwhm
      smni_rstd2 = smni_rstd2/sqrt(log(256.0));
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--g2s")) {
      if (nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%f",&smni_sw1);
      if(smni_sw1 > 1){
	printf("ERROR: g2 w1 = %g > 1\n",smni_sw1);
	exit(1);
      }
      sscanf(pargv[1],"%f",&smni_sstd2); // Actually fwhm
      smni_sstd2 = smni_sstd2/sqrt(log(256.0));
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--seed")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&SynthSeed);
      synth = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--tr")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&TR);
      SetTR = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      outpath = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ar1")) {
      if (nargc < 1) CMDargNErr(option,1);
      ar1path = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ar1red")) {
      if (nargc < 1) CMDargNErr(option,1);
      ar1redfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--fwhmvol")) {
      if (nargc < 1) CMDargNErr(option,1);
      fwhmvolfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--fwhmvolmn")) {
      if (nargc < 1) CMDargNErr(option,1);
      fwhmvolmnfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--fwhmdat")) {
      if (nargc < 1) CMDargNErr(option,1);
      fwhmdatfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--arN")) {
      if (nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&arNlags);
      arNfname = pargv[1];
      nargsused = 2;
    }
    else if (!strcmp(option, "--sum2")) {
      if (nargc < 1) CMDargNErr(option,1);
      sum2file = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--X")) {
      if (nargc < 1) CMDargNErr(option,1);
      Xfile = pargv[0];
      //X = MatrixReadTxt(Xfile, NULL);
      X = MatlabRead(Xfile);
      if (X == NULL) {
        printf("ERROR: reading %s\n",Xfile);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcasecmp(option, "--detrend")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&DetrendOrder);
      if (DetrendOrder > 2) {
        printf("ERROR: cannot have detrending order > 2\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcasecmp(option, "--in_nspmzeropad")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&N_Zero_Pad_Input);
      InValsType = MRI_ANALYZE_FILE;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      int nthreads;
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--soap")){
      // --soap src mask niters out
      if(nargc < 4) CMDargNErr(option,4);
      MRI *src  = MRIread(pargv[0]);
      MRI *mask = MRIread(pargv[1]);
      int niter;
      sscanf(pargv[2],"%d",&niter);
      MRI *out = MRIsoapBubble(src, mask, NULL, niter, .01);
      MRIwrite(out,pargv[3]);
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
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s\n",Progname) ;
  printf("\n");
  printf("   --i inputvol  : input volume\n");
  printf("   --o outputvol : save input after smoothing\n");
  printf("   --save-detrended : detrend output when saving\n");
  printf("   --save-unmasked  : do not mask outputvol\n");
  printf("   --smooth-only    : smooth and save, do not compute fwhm (--so)\n");
  printf("\n");
  printf("   --mask maskvol : binary mask\n");
  printf("   --mask-thresh absthresh : threshold for mask (default is .5)\n");
  printf("   --auto-mask rthresh : compute mask\n");
  printf("   --nerode n : erode mask n times prior to computing fwhm\n");
  printf("   --mask-inv : invert mask\n");
  printf("   --out-mask outmask : save final mask\n");
  printf("\n");
  printf("   --X x.mat : matlab4 detrending matrix\n");
  printf("   --detrend order : polynomial detrending (default 0)\n");
  printf("   --sqr : compute square of input before smoothing\n");
  printf("\n");
  printf("   --fwhm fwhm : smooth BY fwhm before measuring\n");
  printf("   --gstd gstd : same as --fwhm but specified as the stddev\n");
  printf("   --median width : perform median filtering instead of gaussian\n");
  printf("   --fwhmc, --fwhmr, --fwhms to control each axis separately\n");
  printf("   --g2 w1 fwhm2 : gaussian mixture (w1*g1 + (1-w1)*g2)\n");
  printf("   --g2{crs} w1{crs} fwhm2{crs} : assign a mixture model to given axis\n");
  printf("\n");
  printf("   --to-fwhm tofwhm : smooth TO fwhm\n");
  printf("   --to-fwhm-tol tolerance : smooth to fwhm +/- tol (def .5mm)\n");
  printf("   --to-fwhm-nmax nitersmax : maximum number of iterations (def 20)\n");
  printf("   --to-fwhm-file file : save to-fwhm params in file\n");
  printf("\n");
  printf("   --sum sumfile : summary/log\n");
  printf("   --dat datfile : only the final fwhm estimate (base on mean ar1)\n");
  printf("   --fwhmdat fwhmdatfile : compute and save the fwhm of each dim (based on mean ar1)\n");
  printf("   --fwhmvolmn fwhmvolmnfile : compute and save the fwhm of each dim (based on fwhmvol)\n");
  printf("   --fwhmvol fwhmvol : Save 3 frame map of the FWHM at each voxel.\n");
  printf("\n");
  printf("   --synth \n");
  printf("   --synth-frames nframes : default is 10 \n");
  printf("\n");
  printf("   --nframesmin n : require at least this many frames\n");
  printf("   --ispm : input is spm-analyze. Set --i to stem.\n");
  printf("   --in_nspmzeropad nz : zero-padding for spm-analyze\n");
  printf("   --nthreads nthreads : Set OPEN MP threads\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
printf("\n");
printf("FreeSurfer program to estimate the global Gaussian smoothness of a\n");
printf("multi-frame, volume-based data set. The smoothness is measured as the\n");
printf("Full-Width-Half-Max (FWHM) in mm.  Gaussian stddev =\n");
printf("fwhm/sqrt(log(256.0)). The voxels used in the fwhm estimation can be\n");
printf("constrained to be inside of a mask. It is STRONGLY recommended that a\n");
printf("masked be used. The mask can be specified explictly or computed\n");
printf("automatically. By default, the time course will be detrended by\n");
printf("removing the mean. Higher order polynomial detrending is\n");
printf("possible. Alternatively, the user can specify a detrending matrix. The\n");
printf("data can be smoothed BY a given fwhm or TO a given fwhm prior to\n");
printf("estimating the fwhm.  The resulting data can then be saved (thus\n");
printf("turning this program into a smoother). If smoothing is to be done, it\n");
printf("will only be done inside the mask (except see --save-unmasked).\n");
printf("\n");
printf("For more info, see surfer.nmr.mgh.harvard.edu. For help, send email\n");
printf("and summary file to  freesurfer@nmr.mgh.harvard.edu.\n");
printf("\n");
printf("--i inputvol\n");
printf("\n");
printf("Input data. Format must be something readable by mri_convert\n");
printf("(eg, mgh, mgz, img, nii, nii.gz). Alternately, one can synthesize\n");
printf("white gaussian noise with --synth and --synth-frames in which\n");
printf("case inputvol is used as a dimension template.\n");
printf("\n");
printf("--o outputvol\n");
printf("\n");
printf("Save input after smoothing. See also --save-detended and --save-unmasked.\n");
printf("\n");
printf("--smooth-only\n");
printf("\n");
printf("Does not attempt to compute FWHM. Smooths the input, saves to outputvol,\n");
printf("and exists. Respects --save-unmasked, but not --save-detended. This allows\n");
printf("for data sets with fewer than 10 frames to be smoothed.\n");
printf("\n");
printf("--save-detended\n");
printf("\n");
printf("Save input after smoothing, masking, and detrending.\n");
printf("\n");
printf("--save-unmasked\n");
printf("\n");
printf("Save input after smoothing and detrending, but do not mask while\n");
printf("smoothing. Note that the output will be different even inside\n");
printf("the mask because the smoother handles voxels on the boundary\n");
printf("of the mask differently than those at the center.\n");
printf("\n");
printf("--mask maskvol\n");
printf("\n");
printf("Compute FWHM only over voxels in the given mask. Format can be anything\n");
printf("accepted by mri_convert. If smoothing is to be done, it will only be\n");
printf("done inside the mask. It is strongly recommended that a masked be used\n");
printf("(see also --auto-mask and --save-unmasked).\n");
printf("\n");
printf("--mask-thresh thresh\n");
printf("\n");
printf("Threshold mask at thresh. Default is 0.5 (ie, it expects a binary mask).\n");
printf("\n");
printf("--auto-mask rthresh\n");
printf("\n");
printf("Compute a mask based on a fraction (rthresh) of the global mean. If\n");
printf("rthresh is 0.1, then all voxels in the mean image above 0.1*globalmean\n");
printf("are in the mask.\n");
printf("\n");
printf("--mask-inv\n");
printf("\n");
printf("Invert mask, ie, compute FWHM only over voxels outside the given mask.\n");
printf("\n");
printf("--nerode n\n");
printf("\n");
printf("Erode mask n times (ie, make it smaller). Occurs after any mask inversion.\n");
printf("\n");
printf("--out-mask outmaskvol\n");
printf("\n");
printf("Save final mask to outmaskvol.\n");
printf("\n");
printf("--ar1 ar1path\n");
printf("\n");
printf("Save spatial AR1 volume. (6 frames)\n");
printf("\n");
printf("--ar1red ar1redpath\n");
printf("\n");
printf("Save spatial AR1 volume, 3 frames, average of the 6 frame.\n");
printf("\n");
printf("--fwhmvol fwhmvol\n");
printf("\n");
printf("Save 3 frame map of the FWHM at each voxel.\n");
printf("\n");
printf("--fwhmvolmn fwhmvolmn.dat\n");
printf("\n");
printf("Compute mean fwhm from fwhmvol and write into dat file\n");
printf("\n");
printf("--fwhmdat fwhm.dat\n");
printf("\n");
printf("Compute mean fwhm from ar1mean and write into dat file\n");
printf("\n");
printf("--X x.mat\n");
printf("\n");
printf("Detrend/residualize data with the matrix in x.mat. Ie,\n");
printf("y = (I-inv(X'*X)*X')*y, where y is the input. x.mat must\n");
printf("be a matlab4 matrix. Not with --detrend.\n");
printf("\n");
printf("--detrend order\n");
printf("\n");
printf("Detrend data with polynomial of given order. Not with --X. Note:\n");
printf("if neither --X nor --detrend are specified, then detrending\n");
printf("order of 0 is used (ie, the mean is removed).\n");
printf("\n");
printf("--fwhm fwhm\n");
printf("\n");
printf("Smooth BY fwhm mm before estimating the fwhm. This is mainly good for\n");
printf("debuggging. But with --out can also be used to smooth data.\n");
printf("\n");
printf("--to-fwhm tofwhm\n");
printf("\n");
printf("Smooth TO tofwhm mm. This is idea proposed by Lee Friedman\n");
printf("(lfriedman10@comcast.net) and fBIRN (www.birn.net) as a way to\n");
printf("reduce variation across data sets, particularly when the data\n");
printf("may have come from different scanners. The method implemented\n");
printf("here uses an iterative approach in which the data are smoothed\n");
printf("BY varying amounmts until the resulting fwhm is within a certain\n");
printf("tolerance of tofwhm. By default, the tolerance is 0.5 mm, but\n");
printf("can be changed with --to-fwhm-tol. It will iterate at most 20\n");
printf("times (can be changed with --to-fwhm-nmax). An error will be\n");
printf("returned if the tofwhm is less than the inherent fwhm or the max\n");
printf("number of iterations is exceeded. The minimization is done with\n");
printf("a Golden Section Search (see Numerical Recipes in C).\n");
printf("\n");
printf("--to-fwhm-tol tol\n");
printf("\n");
printf("Keep iterating the tofwhm search until the result is within tol\n");
printf("of the desired fwhm (or until the maximum number of iterations\n");
printf("is reached).\n");
printf("\n");
printf("--to-fwhm-nmax niterationsmax\n");
printf("\n");
printf("Maximum number of iterations. Default is 20.\n");
printf("\n");
printf("--to-fwhm-file file\n");
printf("\n");
printf("Save some results of the tofwhm minimization to file. Good for\n");
printf("debugging. Results also saved in summary file.\n");
printf("\n");
printf("--sum sumfile\n");
printf("\n");
printf("Prints summary to ascii sumfile. Send this file when requesting\n");
printf("help or more information.\n");
printf("\n");
printf("--dat datfile\n");
printf("\n");
printf("Prints only the final fwhm estimate into this file.\n");
printf("\n");
printf("--synth\n");
printf("\n");
printf("Synthesize input with white gaussian noise. Ten frames are used by default,\n");
printf("but this can be changed with --synth-frames. Uses input volume as template.\n");
printf("This functionality is useful for degugging. Eg, when using --synth and --fwhm,\n");
printf("it should measure the resulting fwhm to be that passed by --fwhm. Can do\n");
printf("the same with --to-fwhm.\n");
printf("\n");
printf("--synth-frames nframes\n");
printf("\n");
printf("Synthesize input with white gaussian noise with the given number of frames.\n");
printf("Implies --synth.\n");
printf("\n");
printf("--tr TRms\n");
printf("\n");
printf("Set TR in msec (generally not too useful)\n");
printf("\n");
printf("--nthreads nthreads\n");
printf("\n");
printf("Set OPEN MP threads\n");
printf("\n");
printf("--inorm\n");
printf("\n");
printf("Spatial intensity normalization. Subtract the in-mask mean and divide by the in-mask \n");
printf("stddev. \n");
printf("\n");
printf("EXAMPLES:\n");
printf("\n");
printf("1. Measure the fwhm of an input data set, compute mask automatically by\n");
printf("   thresholding input at 20%% of global mean. The time series will be\n");
printf("   have its mean removed prior to computing the fwhm. Save result in\n");
printf("   a summary file (one example uses mgh format, the other gzipped NIFTI):\n");
printf("\n");
printf("      mri_fwhm --i f.mgh    --auto-mask .2 --sum f.fwhm.sum\n");
printf("      mri_fwhm --i f.nii.gz --auto-mask .2 --sum f.fwhm.sum\n");
printf("\n");
printf("2. Same as above, but smooth input BY 5mm fwhm first. Save the\n");
printf("   smoothed output in fsm5.mgh. Save the mask to automask.nii.\n");
printf("   Note: mask is computed on unsmoothed data.\n");
printf("\n");
printf("      mri_fwhm --i f.mgh --auto-mask .2 --sum f.fwhm5.sum\n");
printf("        --fwhm 5 --o fsm5.mgh --out-mask automask.nii\n");
printf("\n");
printf("3. Same as above, but smooth input TO 5 +/- .1mm fwhm first.\n");
printf("   Save the smoothed output in fto5.mgh.\n");
printf("\n");
printf("      mri_fwhm --i f.mgh --auto-mask .2 --sum f.fwhm5.sum\n");
printf("        --to-fwhm-tol .1 --to-fwhm 5 --o fto5.mgh\n");
printf("\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if (inpath == NULL && !synth) {
    printf("ERROR: need to specify --in or --synth\n");
    exit(1);
  }
  if (Xfile && DetrendOrder >= 0) {
    printf("ERROR: cannot use --X and --detrend\n");
    exit(1);
  }
  // At least remove the mean
  if (Xfile == NULL && DetrendOrder < 0) DetrendOrder = 0;
  if (maskpath && automask) {
    printf("ERROR: cannot use --mask and --auto-mask\n");
    exit(1);
  }
  if (outmaskpath && maskpath==NULL && !automask) {
    printf("ERROR: cannot use --outmask without --mask or --auto-mask\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
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
  if (Xfile) fprintf(fp,"xfile %s\n",Xfile);
  if (DetrendOrder >= 0) fprintf(fp,"detrend %d\n",DetrendOrder);
  fprintf(fp,"infwhm %lf\n",infwhm);
  fprintf(fp,"ingstd %lf\n",ingstd);
  fprintf(fp,"tofwhm %lf\n",tofwhm);
  if (tofwhm > 0) {
    fprintf(fp,"byfwhm-tofwhm     %lf\n",byfwhm);
    fprintf(fp,"tofwhm-tol        %lf\n",tofwhmtol);
    fprintf(fp,"tofwhm-niters     %d\n",tofwhmniters);
    fprintf(fp,"tofwhm-niters-max %d\n",tofwhmnitersmax);
  }
  if (maskpath) fprintf(fp,"mask  %s\n",maskpath);
  if (automask) fprintf(fp,"automaskthresh %g\n",automaskthresh);
  if (mask) {
    fprintf(fp,"maskthresh %g\n",maskthresh);
    fprintf(fp,"maskinv %d\n",maskinv);
    fprintf(fp,"nerode  %d\n",nerode);
  }
  if (outmaskpath) {
    fprintf(fp,"outmask  %s\n",outmaskpath);
    fprintf(fp,"SaveDetrended  %d\n",SaveDetrended);
    fprintf(fp,"SaveUnmasked   %d\n",SaveUnmasked);
  }
  if (synth) {
    fprintf(fp,"synth-frames %d\n",nframes);
    fprintf(fp,"seed  %d\n",SynthSeed);
  }
  if (outpath) fprintf(fp,"out  %s\n",outpath);

  return;
}

/*-------------------------------------------------------------------------------
  MRIbinarize2() - same as MRIbinarize() but passes theshold, low, and hi as
  doubles instead of UCHARs.
  -------------------------------------------------------------------------------*/
MRI * MRIbinarize2(MRI *mri_src, MRI *mri_dst,
                   double threshold, double low_val, double hi_val) {
  int     width, height, depth, x, y, z, f ;
  double  val;

  if (!mri_dst)  mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (f = 0 ; f < mri_src->nframes ; f++) {
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          val = MRIgetVoxVal(mri_src, x, y, z, f);
          if (val > threshold) val = hi_val ;
          else                val = low_val ;
          MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
        }
      }
    }
  }

  return(mri_dst) ;
}


/*------------------------------------------------------------------------*/
double EvalFWHM(MRI *vol, MRI *mask) {
  double car1mn, rar1mn, sar1mn;
  double cfwhm,rfwhm,sfwhm,fwhm;
  static MRI *ar1 = NULL;
  ar1 = fMRIspatialAR1(vol, mask, ar1);
  fMRIspatialAR1Mean(ar1, mask, &car1mn, &rar1mn, &sar1mn);
  cfwhm = RFar1ToFWHM(car1mn, vol->xsize);
  rfwhm = RFar1ToFWHM(rar1mn, vol->ysize);
  sfwhm = RFar1ToFWHM(sar1mn, vol->zsize);
  fwhm = sqrt((cfwhm*cfwhm + rfwhm*rfwhm + sfwhm*sfwhm)/3.0);
  return(fwhm);
}
/*------------------------------------------------------------------------*/
MRI *MRImaskedGaussianSmoothTo(MRI *invol, MRI *mask, double ToFWHM,
                               double tol, int nitersmax,
                               double *pByFWHM, double *pToFWHMActual, int *niters,
                               MRI *outvol) {
  double SrcFWHM, ByGStd;
  MRI *volsm;
  double ya, yb, yc, xa, xb, xc, s1, s2;
  double xn,yn;
  double C, R, err;
  int nth;

  C = (3.0-sqrt(5.0))/2.0; // golden mean
  R = 1-C;

  *niters = 0;
  SrcFWHM = EvalFWHM(invol, mask);
  if (SrcFWHM > ToFWHM) {
    printf("ERROR: MRImaskedGaussianSmoothTo(): the inherent smoothness\n");
    printf("of the data set is about fwhm=%gmm, which is more than the\n",SrcFWHM);
    printf("amount that you want to smooth it to (%gmm). It is impossible\n",ToFWHM);
    printf("to 'unsmooth' the data.\n");
    return(NULL);
  }
  xa = 0;
  ya = SrcFWHM;
  printf("Unsmoothed actual fwhm is %g\n",ya);

  // Check whether we are close enough already
  if (fabs(SrcFWHM-ToFWHM) < tol) {
    *pByFWHM = 0.0;
    *pToFWHMActual = SrcFWHM;
    outvol = MRIcopy(invol,outvol);
    return(outvol);
  }

  volsm = MRIcopy(invol,NULL); // allocate

  // First point in the bracket

  // Second point in the bracket
  (*niters)++;
  xb = sqrt(ToFWHM*ToFWHM - SrcFWHM*SrcFWHM); // power law
  ByGStd = xb/sqrt(log(256.0));
  printf("Trying smoothing by fwhm %g  ",xb);
  fflush(stdout);
  MRImaskedGaussianSmooth(invol, mask, ByGStd, volsm);
  yb = EvalFWHM(volsm, mask);
  printf("results in actual fwhm of %g\n",yb);
  // Check whether we are close enough now
  if (fabs(yb-ToFWHM) < tol) {
    *pByFWHM = xb;
    *pToFWHMActual = yb;
    outvol = MRIcopy(volsm,outvol);
    MRIfree(&volsm);
    return(outvol);
  }

  // Third point in the bracket. Not sure how to choose this point,
  // it needs to be far enough to bracket the min, but too far
  // and we end up doing to many evaluations.
  (*niters)++;
  xc = xb + (xb - xa); //
  ByGStd = xc/sqrt(log(256.0));
  printf("Trying smoothing by fwhm %g  ",xc);
  fflush(stdout);
  MRImaskedGaussianSmooth(invol, mask, ByGStd, volsm);
  yc = EvalFWHM(volsm, mask);
  printf("results in actual fwhm of %g\n",yc);
  // Check whether we are close enough now
  if (fabs(yc-ToFWHM) < tol) {
    *pByFWHM = xc;
    *pToFWHMActual = yc;
    outvol = MRIcopy(volsm,outvol);
    MRIfree(&volsm);
    return(outvol);
  }
  if (yc < ToFWHM) {
    // Did not step far enough out
    printf("ERROR: did not step far enough out\n");
    return(NULL);
  }

  // ok, we've brackated the min, now chase it down like a scared rabbit
  // using a golden section search (see numerical recipes in C).
  printf("Beginning golden section search.\n");
  printf("Expecting roughly %d iterations will be needed.\n",
         (int)ceil(log(tol/(xc-xa))/log(R)));
  nth = 0;
  getybest(xa, ya, xb, yb, xc, yc, pByFWHM, pToFWHMActual, ToFWHM);
  err = fabs(*pToFWHMActual-ToFWHM);
  while ( err > tol) {
    if (nth > nitersmax) {
      printf("ERROR: searched timed out at niters=%d\n",nth);
      MRIfree(&volsm);
      return(NULL);
    }
    printf("n=%d by=(%4.2lf,%4.2lf,%4.2lf) to=(%4.2lf,%4.2lf,%4.2lf) err=%g\n",
           nth, xa, xb, xc, ya, yb, yc, err);
    fflush(stdout);
    s1 = xb - xa;
    s2 = xc - xb;
    if (s1 > s2) { // s1 is bigger
      xn = xa + R*s1;
      printf("   Trying smoothing by fwhm %g  ",xn);
      fflush(stdout);
      ByGStd = xn/sqrt(log(256.0));
      MRImaskedGaussianSmooth(invol, mask, ByGStd, volsm);
      yn = EvalFWHM(volsm, mask);
      printf("results in actual fwhm of %g\n",yn);
      if (fabs(yn-ToFWHM) < fabs(ToFWHM-yb)) {
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
      printf("   Trying smoothing by fwhm %g  ",xn);
      fflush(stdout);
      MRImaskedGaussianSmooth(invol, mask, ByGStd, volsm);
      yn = EvalFWHM(volsm, mask);
      printf("results in actual fwhm of %g\n",yn);
      if (fabs(yn-ToFWHM) < fabs(ToFWHM-yb)) {
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
    (*niters)++;
  }
  printf("niters=%3d  by=%6.4lf  to=%6.4lf targ=%6.4lf  err=%g\n",
         nth, *pByFWHM, *pToFWHMActual, ToFWHM, err);
  outvol = MRIcopy(volsm,outvol);
  MRIfree(&volsm);
  return(outvol);
}


int getybest(double xa, double ya, double xb, double yb, double xc, double yc,
             double *xbest, double *ybest, double ytarg) {
  double ea, eb, ec;

  ea = fabs(ya-ytarg);
  eb = fabs(yb-ytarg);
  ec = fabs(yc-ytarg);

  if (ea < eb && ea < ec) {
    *xbest = xa;
    *ybest = ya;
    return(0);
  }
  if (eb < ec) {
    *xbest = xb;
    *ybest = yb;
    return(0);
  }
  *xbest = xc;
  *ybest = yc;
  return(0);
}


int fMRIspatialFWHMMean(MRI *fwhmvol, MRI *mask, double *cfwhmmn, double *rfwhmmn, double *sfwhmmn)
{
  int c, r, s;
  long nhits;
  double m, cfwhmsum, rfwhmsum, sfwhmsum;

  cfwhmsum = 0.0;
  rfwhmsum = 0.0;
  sfwhmsum = 0.0;
  nhits = 0;
  for (c = 1; c < fwhmvol->width - 1; c++) {
    for (r = 1; r < fwhmvol->height - 1; r++) {
      for (s = 1; s < fwhmvol->depth - 1; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        if(MRIgetVoxVal(fwhmvol, c, r, s, 0) == 0) continue;
        cfwhmsum += MRIgetVoxVal(fwhmvol, c, r, s, 0);
        rfwhmsum += MRIgetVoxVal(fwhmvol, c, r, s, 1);
        sfwhmsum += MRIgetVoxVal(fwhmvol, c, r, s, 2);
        nhits++;
      }
    }
  }

  *cfwhmmn = (cfwhmsum / nhits);
  *rfwhmmn = (rfwhmsum / nhits);
  *sfwhmmn = (sfwhmsum / nhits);

  return(0);
}
