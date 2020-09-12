/**
 * @brief GLM analysis with or without FSGD files
 *
 * Performs general linear model (GLM) analysis in the volume or the
 * surface.  Options include simulation for correction for multiple
 * comparisons, weighted LMS, variance smoothing, PCA/SVD analysis of
 * residuals, per-voxel design matrices, and 'self' regressors. This
 * program performs both the estimation and inference. This program
 * is meant to replace mris_glm (which only operated on surfaces).
 * This program can be run in conjunction with mris_preproc.
 */
/*
 * Original Author: Douglas N Greve
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

BEGINUSAGE --------------------------------------------------------------

USAGE: ./mri_glmfit

   --glmdir dir : save outputs to dir

   --y inputfile
   --table stats-table : as output by asegstats2table or aparcstats2table 
   --fsgd FSGDF <gd2mtx> : freesurfer descriptor file
   --X design matrix file
   --C contrast1.mtx <--C contrast2.mtx ...>
   --osgm : construct X and C as a one-sample group mean
   --no-contrasts-ok : do not fail if no contrasts specified
   --dti bvals bvecs : do DTI analysis using bvals and bvecs
   --dti siemensdicom : do DTI analysis extracting bvals and bvecs from dicom
   --dti-X X.txt : do DTI analysis using provided matrix

   --pvr pvr1 <--prv pvr2 ...> : per-voxel regressors
   --selfreg col row slice   : self-regressor from index col row slice

   --wls yffxvar : weighted least squares
   --yffxvar yffxvar : for fixed effects analysis
   --ffxdof DOF : dof for fixed effects analysis
   --ffxdofdat ffxdof.dat : text file with dof for fixed effects analysis

   --w weightfile : weight for each input at each voxel
   --w-inv : invert weights
   --w-sqrt : sqrt of (inverted) weights

   --fwhm fwhm : smooth input by fwhm
   --var-fwhm fwhm : smooth variance by fwhm
   --no-mask-smooth : do not mask when smoothing
   --no-est-fwhm : turn off FWHM output estimation

   --mask maskfile : binary mask
   --label labelfile : use label as mask, surfaces only
   --no-mask : do NOT use a mask (same as --no-cortex)
   --no-cortex : do NOT use subjects ?h.cortex.label as --label
   --mask-inv : invert mask
   --prune : remove voxels that do not have a non-zero value at each frame (def)
   --no-prune : do not prune
   --logy : compute natural log of y prior to analysis
   --no-logy : do not compute natural log of y prior to analysis
   --rm-spatial-mean : subtract the (masked) mean from each frame
   --yhat-save : save signal estimate (yhat)
   --eres-save : save residual error (eres)
   --eres-scm : save residual error spatial correlation matrix (eres.scm). Big!
   --y-out y.out.mgh : save input after any pre-processing

   --surf subject hemi <surfname> : needed for some flags (uses white by default)

   --skew : compute skew and p-value for skew
   --kurtosis : compute kurtosis and p-value for kurtosis

   --sim nulltype nsim thresh csdbasename : simulation perm, mc-full, mc-z
   --sim-sign signstring : abs, pos, or neg. Default is abs.
   --uniform min max : use uniform distribution instead of gaussian

   --pca : perform pca/svd analysis on residual
   --tar1 : compute and save temporal AR1 of residual
   --save-yhat : flag to save signal estimate
   --save-cond  : flag to save design matrix condition at each voxel
   --voxdump col row slice  : dump voxel GLM and exit

   --seed seed : used for synthesizing noise
   --synth : replace input with gaussian

   --resynthtest niters : test GLM by resynthsis
   --profile     niters : test speed

   --mrtm1 RefTac TimeSec : perform MRTM1 kinetic modeling
   --mrtm2 RefTac TimeSec k2prime : perform MRTM2 kinetic modeling

   --perm-force : force perumtation test, even when design matrix is not orthog
   --diag Gdiag_no : set diagnositc level
   --diag-cluster : save sig volume and exit from first sim loop
   --debug     turn on debugging
   --checkopts don't run anything, just check options and exit
   --help      print out information on how to use this program
   --version   print out version and exit
   --no-fix-vertex-area : turn off fixing of vertex area (for back comapt only)
   --allowsubjrep allow subject names to repeat in the fsgd file (must appear
                  before --fsgd)
   --allow-zero-dof : mostly for very special purposes
   --illcond : allow ill-conditioned design matrices
   --sim-done SimDoneFile : create DoneFile when simulation finished 

ENDUSAGE --------------------------------------------------------------

BEGINHELP --------------------------------------------------------------

OUTLINE:
  SUMMARY
  MATHEMATICAL BACKGROUND
  COMMAND-LINE ARGUMENTS
  MONTE CARLO SIMULATION AND CORRECTION FOR MULTIPLE COMPARISONS

SUMMARY

Performs general linear model (GLM) analysis in the volume or the
surface.  Options include simulation for correction for multiple
comparisons, weighted LMS, variance smoothing, PCA/SVD analysis of
residuals, per-voxel design matrices, and 'self' regressors. This
program performs both the estimation and inference. This program
is meant to replace mris_glm (which only operated on surfaces).
This program can be run in conjunction with mris_preproc.

MATHEMATICAL BACKGROUND

This brief intoduction to GLM theory is provided to help the user
understand what the inputs and outputs are and how to set the
various parameters. These operations are performed at each voxel
or vertex separately (except with --var-fwhm).

The forward model is given by:

    y = W*X*B + n

where X is the Ns-by-Nb design matrix, y is the Ns-by-Nv input data
set, B is the Nb-by-Nv regression parameters, and n is noise. Ns is
the number of inputs, Nb is the number of regressors, and Nv is the
number of voxels/vertices (all cols/rows/slices). y may be surface
or volume data and may or may not have been spatially smoothed. W
is a diagonal weighing matrix.

During the estimation stage, the forward model is inverted to
solve for B:

    B = inv(X'W'*W*X)*X'W'y

The signal estimate is computed as

    yhat = B*X

The residual error is computed as

    eres = y - yhat

For random effects analysis, the noise variance estimate (rvar) is
computed as the sum of the squares of the residual error divided by
the DOF.  The DOF equals the number of rows of X minus the number of
columns. For fixed effects analysis, the noise variance is estimated
from the lower-level variances passed with --yffxvar, and the DOF
is the sum of the DOFs from the lower level.

A contrast matrix C has J rows and as many columns as columns of
X. The contrast is then computed as:

   G = C*B

The F-ratio for the contrast is then given by:

   F = G'*inv(C*inv(X'W'*W*X))*C')*G/(J*rvar)

The F is then used to compute a p-value.  Note that when J=1, this
reduces to a two-tailed t-test.


COMMAND-LINE ARGUMENTS

--glmdir dir

Directory where output will be saved. Not needed with --sim.

The outputs will be saved in mgh format as:
  mri_glmfit.log - execution parameters (send with bug reports)
  beta.mgh - all regression coefficients (B above)
  eres.mgh - residual error
  rvar.mgh - residual error variance
  rstd.mgh - residual error stddev (just sqrt of rvar)
  y.fsgd - fsgd file (if one was input)
  wn.mgh - normalized weights (with --w or --wls)
  yhat.mgh - signal estimate (with --save-yhat)
  mask.mgh - final mask (when a mask is used)
  cond.mgh - design matrix condition at each voxel (with --save-cond)
  contrast1name/ - directory for each contrast (see --C)
    C.dat - copy of contrast matrix
    gamma.mgh - contrast (G above)
    F.mgh - F-ratio (t = sign(gamma)*sqrt(F) for t-test contrasts)
    sig.mgh - significance from F-test (actually -log10(p))
    z.mgh - z map computed from the p-value
    pcc.mgh - partial correlation coefficient (for t-tests)

--y inputfile

Path to input file with each frame being a separate input. This can be
volume or surface-based, but the file must be one of the 'volume'
formats (eg, mgh, img, nii, etc) accepted by mri_convert. See
mris_preproc for an easy way to generate this file for surface data.
Not with --table.

--table stats-table

Use text table as input instead of --y. The stats-table is that of
the form produced by asegstats2table or aparcstats2table.

--fsgd fname <gd2mtx>

Specify the global design matrix with a FreeSurfer Group Descriptor
File (FSGDF).  See http://surfer.nmr.mgh.harvard.edu/docs/fsgdf.txt
for more info.  The gd2mtx is the method by which the group
description is converted into a design matrix. Legal values are doss
(Different Offset, Same Slope) and dods (Different Offset, Different
Slope). doss will create a design matrix in which each class has it's
own offset but forces all classes to have the same slope. dods models
each class with it's own offset and slope. In either case, you'll need
to know the order of the regressors in order to correctly specify the
contrast vector. For doss, the first NClass columns refer to the
offset for each class.  The remaining columns are for the continuous
variables. In dods, the first NClass columns again refer to the offset
for each class.  However, there will be NClass*NVar more columns (ie,
one column for each variable for each class). The first NClass columns
are for the first variable, etc. If neither of these models works for
you, you will have to specify the design matrix manually (with --X).

--no-rescale-x

By default the inverse of the covariance of the desgin matrix is
computed by rescaling each column of the design matrix prior to the 
inverse computation, then rescaling back afterwards. This helps
with designs that are badly scaled. This is completely transparent
to the user. This flag turns this feature off so that the inverse
is computed directly from the design matrix.

--X design matrix file

Explicitly specify the design matrix. Can be in simple text or in matlab4
format. If matlab4, you can save a matrix with save('X.mat','X','-v4');

--C contrast1.mtx <--C contrast2.mtx ...>

Specify one or more contrasts to test. The contrast.mtx file is an
ASCII text file with the contrast matrix in it (make sure the last
line is blank). The name can be (almost) anything. If the extension is
.mtx, .mat, .dat, or .con, the extension will be stripped of to form
the directory output name.  The output will be saved in
glmdir/contrast1, glmdir/contrast2, etc. Eg, if --C norm-v-cont.mtx,
then the ouput will be glmdir/norm-v-cont.

--osgm

Construct X and C as a one-sample group mean. X is then a one-column
matrix filled with all 1s, and C is a 1-by-1 matrix with value 1.
You cannot specify both --X and --osgm. A contrast cannot be specified
either. The contrast name will be osgm.

--pvr pvr1 <--prv pvr2 ...>

Per-voxel (or vertex) regressors (PVR). Normally, the design matrix is
'global', ie, the same matrix is used at each voxel. This option allows the
user to specify voxel-specific regressors to append to the design
matrix. Note: the contrast matrices must include columns for these
components.

--selfreg col row slice

Create a 'self-regressor' from the input data based on the waveform at
index col row slice. This waveform is residualized and then added as a
column to the design matrix. Note: the contrast matrices must include
columns for this component.

--wls yffxvar : weighted least squares

Perform weighted least squares (WLS) random effects analysis instead
of ordinary least squares (OLS).  This requires that the lower-level
variances be available.  This is often the case with fMRI analysis but
not with an anatomical analysis. Note: this should not be confused
with fixed effects analysis. The weights will be inverted,
square-rooted, and normalized to sum to the number of inputs for each
voxel. Same as --w yffxvar --w-inv --w-sqrt (see --w below).

--yffxvar yffxvar      : for fixed effects analysis
--ffxdof DOF           : DOF for fixed effects analysis
--ffxdofdat ffxdof.dat : text file with DOF for fixed effects analysis

Perform fixed-effect analysis. This requires that the lower-level variances
be available.  This is often the case with fMRI analysis but not with
an anatomical analysis. Note: this should not be confused with weighted
random effects analysis (wls). The dof is the sum of the DOFs from the
lower levels.

--w weightfile
--w-inv
--w-sqrt

Perform weighted LMS using per-voxel weights from the weightfile. The
data in weightfile must have the same dimensions as the input y
file. If --w-inv is flagged, then the inverse of each weight is used
as the weight.  If --w-sqrt is flagged, then the square root of each
weight is used as the weight.  If both are flagged, the inverse is
done first. The final weights are normalized so that the sum at each
voxel equals the number of inputs. The normalized weights are then
saved in glmdir/wn.mgh.  The --w-inv and --w-sqrt flags are useful
when passing contrast variances from a lower level analysis to a
higher level analysis (as is often done in fMRI).

--fwhm fwhm

Smooth input with a Gaussian kernel with the given full-width/half-maximum
(fwhm) specified in mm. If the data are surface-based, then you must
specify --surf, otherwise mri_glmfit assumes that the input is a volume
and will perform volume smoothing.

--var-fwhm fwhm

Smooth residual variance map with a Gaussian kernel with the given
full-width/half-maximum (fwhm) specified in mm. If the data are
surface-based, then you must specify --surf, otherwise mri_glmfit
assumes that the input is a volume and will perform volume smoothing.

--mask maskfile
--label labelfile
--mask-inv
--cortex 

Only perform analysis where mask=1. All other voxels will be set to 0.
If using surface, then labelfile will be converted to a binary mask
(requires --surf). By default, the label file for surfaces is 
?h.cortex.label. To force a no-mask with surfaces, use --no-mask or 
--no-cortex. If --mask-inv is flagged, then performs analysis
only where mask=0. If performing a simulation (--sim), map maximums
and clusters will only be searched for in the mask. The final binary
mask will automatically be saved in glmdir/mask.mgh

--prune
--no-prune

This happens by default. Use --no-prune to turn it off. Remove voxels
from the analysis if the ALL the frames at that voxel
do not have an absolute value that exceeds zero (actually FLT_MIN, 
or whatever is set by --prune_thr). This helps to prevent the situation 
where some frames are 0 and others are not. If no mask is supplied, 
a mask is created and saved. If a mask is supplied, it is pruned, and 
the final mask is saved. Do not use with --sim. Rather, run the non-sim 
analysis with --prune, then pass the created mask when running simulation. 
It is generally a good idea to prune. --no-prune will turn off pruning 
if it had been turned on. For DTI, only the first frame is used to 
create the mask.

--prune_thr threshold

Use threshold to create the mask using pruning. Default is FLT_MIN

--surf subject hemi <surfname>

Specify that the input has a surface geometry from the hemisphere of the
given FreeSurfer subject. This is necessary for smoothing surface data
(--fwhm or --var-fwhm), specifying a label as a mask (--label), or
running a simulation (--sim) on surface data. If --surf is not specified,
then mri_glmfit will assume that the data are volume-based and use
the geometry as specified in the header to make spatial calculations.
By default, the white surface is used, but this can be overridden by
specifying surfname.

--pca

Flag to perform PCA/SVD analysis on the residual. The result is stored
in glmdir/pca-eres as v.mgh (spatial eigenvectors), u.mtx (frame
eigenvectors), sdiag.mat (singular values). eres = u*s*v'. The matfiles
are just ASCII text. The spatial EVs can be loaded as overlays in
tkmedit or tksurfer. In addition, there is stats.dat with 5 columns:
  (1) component number
  (2) variance spanned by that component
  (3) cumulative variance spanned up to that component
  (4) percent variance spanned by that component
  (5) cumulative percent variance spanned up to that component

--save-yhat

Flag to save the signal estimate (yhat) as glmdir/yhat.mgh. Normally, this
pis not very useful except for debugging.

--save-cond

Flag to save the condition number of the design matrix at eaach voxel.
Normally, this is not very useful except for debugging. It is totally
useless if not using weights or PVRs.

--nii, --nii.gz

Use nifti (or compressed nifti) as output format instead of mgh. This
will work with surfaces, but you will not be able to open the output
nifti files with non-freesurfer software.

--seed seed

Use seed as the seed for the random number generator. By default, mri_glmfit
will select a seed based on time-of-day. This only has an effect with
--sim or --synth.

--synth

Replace input data with whise gaussian noise. This is good for testing.

--voxdump col row slice

Save GLM data for a single voxel in directory glmdir/voxdump-col-row-slice.
Exits immediately. Good for debugging.


MONTE CARLO SIMULATION AND CORRECTION FOR MULTIPLE COMPARISONS

One method for correcting for multiple comparisons is to perform simulations
under the null hypothesis and see how often the value of a statistic
from the 'true' analysis is exceeded. This frequency is then interpreted
as a p-value which has been corrected for multiple comparisons. This
is especially useful with surface-based data as traditional random
field theory is harder to implement. This simulator is roughly based
on FSLs permuation simulator (randomise) and AFNIs null-z simulator
(AlphaSim). Note that FreeSurfer also offers False Discovery Rate (FDR)
correction in tkmedit and tksurfer.

The estimation, simulation, and correction are done in three distinct
phases:
  1. Estimation: run the analysis on your data without simulation.
     At this point you can view your results (see if FDR is
     sufficient:).
  2. Simulation: run the simulator with the same parameters
     as the estimation to get the Cluster Simulation Data (CSD).
  3. Clustering: run mri_surfcluster or mri_volcluster with the CSD
     from the simulator and the output of the estimation. These
     programs will print out clusters along with their p-values.

The Estimation step is described in detail above. The simulation
is invoked by calling mri_glmfit with the following arguments:

--sim nulltype nsim thresh csdbasename
--sim-sign sign

It is not necessary to specify --glmdir (it will be ignored). If
you are analyzing surface data, then include --surf.

nulltype is the method of generating the null data. Legal values are:
  (1) perm - perumation, randomly permute rows of X (cf FSL randomise)
  (2) mc-full - replace input with white gaussian noise
  (3) mc-z - do not actually do analysis, just assume the output
      is z-distributed (cf ANFI AlphaSim)
nsim - number of simulation iterations to run (see below)
thresh - threshold, specified as -log10(pvalue) to use for clustering
csdbasename - base name of the file to store the CSD data in. Each
  contrast will get its own file (created by appending the contrast
  name to the base name). A '.csd' is appended to each file name.

Multiple simulations can be run in parallel by specifying different
csdbasenames. Then pass the multiple CSD files to mri_surfcluster
and mri_volcluster. The Full CSD file is written on each iteration,
which means that the CSD file will be valid if the simulation
is aborted or crashes.

In the cases where the design matrix is a single columns of ones
(ie, one-sample group mean), it makes no sense to permute the
rows of the design matrix. mri_glmfit automatically checks
for this case. If found, the design matrix is rebuilt on each
permutation with randomly selected +1 and -1s. Same as the -1
option to FSLs randomise.

--sim-sign sign

sign is either abs (default), pos, or neg. pos/neg tell mri_glmfit to
perform a one-tailed test. In this case, the contrast matrix can
only have one row.

--uniform min max

For mc-full, synthesize input as a uniform distribution between min
and max. 

ENDHELP --------------------------------------------------------------

*/

// Things to do:

// Save some sort of config in output dir.

// Leave-one-out
// Copies or Links to source data?

// Check to make sure no two contrast names are the same
// Check to make sure no two contrast mtxs are the same
// p-to-z
// Rewrite MatrixReadTxt to ignore # and % and empty lines
// Auto-det/read matlab4 matrices

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <float.h>
#include <errno.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "matrix.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "randomfields.h"
#include "dti.h"
#include "image.h"
#include "stats.h"

int MRISmaskByLabel(MRI *y, MRIS *surf, LABEL *lb, int invflag);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
static int SmoothSurfOrVol(MRIS *surf, MRI *mri, MRI *mask, double SmthLevel);

int main(int argc, char *argv[]) ;

const char *Progname = "mri_glmfit";

int SynthSeed = -1;

char *XFile=NULL, *betaFile=NULL, *rvarFile=NULL;
char *yhatFile=NULL, *eresFile=NULL, *maskFile=NULL;
char *wgFile=NULL,*wFile=NULL;
char *eresSCMFile=NULL;
char *condFile=NULL;
char *GLMDir=NULL;
char *pvrFiles[50];
int yhatSave=0;
int eresSave=0;
int SaveFWHMMap=0;
int eresSCMSave=0;
int condSave=0;
std::string yFile, yOutFile, yffxvarFile;

char *labelFile=NULL;
LABEL *clabel=NULL;
int   maskinv = 0;
int   nmask, nvoxels;
float maskfraction, voxelsize;
int   prunemask = 1;

MRI *mritmp=NULL, *mritmp2=NULL, *sig=NULL, *rstd, *fsnr;

int debug = 0, checkoptsonly = 0;
char tmpstr[2000];

int nContrasts=0;
char *CFile[100];
int err,c,r,s;
double Xcond;

int npvr=0;
MRIGLM *mriglm=NULL, *mriglmtmp=NULL;

double FWHM=0;
double SmoothLevel=0;
double VarFWHM=0;
double VarSmoothLevel=0;
int UseMaskWithSmoothing = 1;
double ResFWHM;

char voxdumpdir[1000];
int voxdump[3];
int voxdumpflag = 0;

char *fsgdfile = NULL;
FSGD *fsgd=NULL;
const char  *gd2mtx_method = "none";
int fsgdReScale = 0; 
int ReScaleX = 1; 

int nSelfReg = 0;
int crsSelfReg[100][3];
char *SUBJECTS_DIR;
int cmax, rmax, smax;
double Fmax, sigmax;

int pcaSave=0;
int npca = -1;
MATRIX *Upca=NULL,*Spca=NULL;
MRI *Vpca=NULL;

struct utsname uts;
char *cmdline, cwd[2000];

char *MaxVoxBase = NULL;
int DontSave = 0;
int DontSaveWn = 0;

int DoSim=0;
int synth = 0;
int PermForce = 0;
int UseUniform = 0;
double UniformMin = 0;
double UniformMax = 0;

SURFCLUSTERSUM *SurfClustList;
int nClusters;
char *subject=NULL, *hemi=NULL, *simbase=NULL;
MRI_SURFACE *surf=NULL;
int nsim,nthsim;
double csize;
MRI *fwhmmap = NULL;

VOLCLUSTER **VolClustList;

int DiagCluster=0;

double InterVertexDistAvg, InterVertexDistStdDev, avgvtxarea;
double ar1mn, ar1std, ar1max;
double eresgstd, eresfwhm, searchspace;
double car1mn, rar1mn,sar1mn,cfwhm,rfwhm,sfwhm;
MRI *ar1=NULL, *tar1=NULL, *z=NULL, *zabs=NULL, *cnr=NULL;

CSD *csd;
RFS *rfs;
int weightinv=0, weightsqrt=0;

int OneSamplePerm=0;
int OneSampleGroupMean=0;
int PermNonStatCor = 0;
Timer mytimer;
int ReallyUseAverage7 = 0;
int logflag = 0; // natural log
float prune_thr = FLT_MIN;

DTI *dti;
int usedti = 0;
int usepruning = 0;
MRI *lowb, *tensor, *evals, *evec1, *evec2, *evec3;
MRI  *fa, *ra, *vr, *adc, *dwi, *dwisynth,*dwires,*dwirvar;
MRI  *ivc, *k, *pk;
char *bvalfile=NULL, *bvecfile=NULL;

int useasl = 0;
double asl1val = 1, asl2val = 0;

int useqa = 0;

const char *format = "mgh";
const char *surfname = "white";

int SubSample = 0;
int SubSampStart = 0;
int SubSampDelta = 0;

int DoDistance = 0;

int DoTemporalAR1 = 0;
int DoFFx = 0;
IMAGE *I;

int IllCondOK = 0;
int NoContrastsOK = 0;
int ComputeFWHM = 1;

int UseStatTable = 0;
STAT_TABLE *StatTable=NULL, *OutStatTable=NULL, *GammaStatTable=NULL;
int  UseCortexLabel = 1;

char *SimDoneFile = NULL;
int tSimSign = 0;
int FWHMSet = 0;
int DoKurtosis = 0;
int DoSkew = 0;

char *Gamma0File[GLMMAT_NCONTRASTS_MAX];
MATRIX *Xtmp=NULL, *Xnorm=NULL;
char *XOnlyFile = NULL;

char *frameMaskFile = NULL;

int DoSimThreshLoop = 0;
int  nThreshList = 4, nthThresh;
float ThreshList[4] = {1.3,  2.0,  2.3,  3.0};
int  nSignList = 3, nthSign;
int SignList[3] = {-1,0,1};
CSD *csdList[5][3][20];

MATRIX *RTM_Cr, *RTM_intCr, *RTM_TimeSec, *RTM_TimeMin;
int DoMRTM1=0;
int DoMRTM2=0;
int DoLogan=0;
double MRTM2_k2p=0, Logan_Tstar=0;
MATRIX *MRTM2_x1;

int nRandExclude=0,  *ExcludeFrames=NULL, nExclude=0;
char *ExcludeFrameFile=NULL;
MATRIX *MatrixExcludeFrames(MATRIX *Src, int *ExcludeFrames, int nExclude);
MRI *fMRIexcludeFrames(MRI *f, int *ExcludeFrames, int nExclude, MRI *fex);
int AllowZeroDOF=0;
MRI *BindingPotential(MRI *k2, MRI *k2a, MRI *mask, MRI *bp);
int DoReshape = 0;
MRI *MRIconjunct3(MRI *sig1, MRI *sig2, MRI *sig3, MRI *mask, MRI *c123);

int NSplits=0, SplitNo=0;
int SplitMin, SplitMax, nPerSplit, RandSplit;
int DoFisher = 0; 
int DoPCC=1;
int RmSpatialMean = 0;

double GLMEfficiency(MATRIX *X, MATRIX *C);
int GLMdiagnoseDesignMatrix(MATRIX *X);
MRI *fMRIskew(MRI *y, MRI *mask);
MRI *MRIpskew(MRI *kvals, int dof, MRI *mask, int nsamples);
MRI *MRIremoveSpatialMean(MRI *vol, MRI *mask, MRI *out);
int MRIloganize(MATRIX **X, MRI **Ct, MRI **intCt, const MATRIX *t, const double tstar);

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, n,m;
  int msecFitTime;
  MATRIX *wvect=NULL, *Mtmp=NULL, *Xselfreg=NULL, *Ex=NULL, *XgNew=NULL;
  MATRIX *Ct, *CCt;
  FILE *fp;
  double Ccond, dtmp, threshadj, eff;
  const char *tmpstr2=NULL;

  eresfwhm = -1;
  csd = CSDalloc();
  csd->threshsign = 0; //0=abs,+1,-1

  nargs = handleVersionOption(argc, argv, "mri_glmfit");
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
  mriglm = (MRIGLM *) calloc(sizeof(MRIGLM),1);
  mriglm->glm = GLMalloc();
  mriglm->ffxdof = 0;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);

  mriglm->glm->DoPCC = DoPCC;
  mriglm->glm->ReScaleX = ReScaleX;

  // Seed the random number generator just in case
  if (SynthSeed < 0) SynthSeed = PDFtodSeed();
  srand48(SynthSeed);

  if (surf != NULL) {
    MRIScomputeMetricProperties(surf);
    InterVertexDistAvg    = surf->avg_vertex_dist;
    InterVertexDistStdDev = surf->std_vertex_dist;
    avgvtxarea = surf->avg_vertex_area;
    printf("Number of vertices %d\n",surf->nvertices);
    printf("Number of faces    %d\n",surf->nfaces);
    printf("Total area         %lf\n",surf->total_area);
    printf("AvgVtxArea       %lf\n",avgvtxarea);
    printf("AvgVtxDist       %lf\n",InterVertexDistAvg);
    printf("StdVtxDist       %lf\n",InterVertexDistStdDev);
  }

  // Compute number of iterations for surface smoothing
  if (FWHM > 0 && surf != NULL) {
    SmoothLevel = MRISfwhm2nitersSubj(FWHM, subject, hemi, surfname);
    printf("Surface smoothing by fwhm=%lf, niters=%lf\n",FWHM,SmoothLevel);
  } else SmoothLevel = FWHM;

  if (VarFWHM > 0 && surf != NULL) {
    VarSmoothLevel = MRISfwhm2nitersSubj(VarFWHM, subject, hemi, "white");
    printf("Variance surface smoothing by fwhm=%lf, niters=%lf\n",
           VarFWHM,VarSmoothLevel);
  } else VarSmoothLevel = VarFWHM;

  dump_options(stdout);

  // Create the output directory
  if (! DontSave) {
    if (GLMDir != NULL) {
      printf("Creating output directory %s\n",GLMDir);
      err = mkdir(GLMDir,0777);
      if (err != 0 && errno != EEXIST) {
        printf("ERROR: creating directory %s\n",GLMDir);
        perror(NULL);
        return(1);
      }
    }
    sprintf(tmpstr,"%s/mri_glmfit.log",GLMDir);
    fp = fopen(tmpstr,"w");
    dump_options(fp);
    fclose(fp);
    if(subject){
      sprintf(tmpstr,"%s/surface",GLMDir);
      fp = fopen(tmpstr,"w");
      fprintf(fp,"%s %s\n",subject,hemi);
      fclose(fp);
    }
  }

  mriglm->npvr     = npvr;
  mriglm->yhatsave = yhatSave;
  mriglm->condsave = condSave;

  // Load input--------------------------------------
  printf("Loading y from %s\n",yFile.c_str());
  fflush(stdout);
  if(! UseStatTable){
    mriglm->y = MRIread(yFile.c_str());
    printf("   ... done reading.\n");
    fflush(stdout);
    if (mriglm->y == NULL) {
      printf("ERROR: loading y %s\n",yFile.c_str());
      exit(1);
    }
    nvoxels = mriglm->y->width * mriglm->y->height * mriglm->y->depth;
    if(nvoxels == 163842 && surf == NULL){
      printf("ERROR: you must use '--surface subject hemi' with surface data\n");
      exit(1);
    }
    if(DoReshape && surf != NULL){
      printf("Forcing reshape to 1d\n");
      mritmp = mri_reshape(mriglm->y,nvoxels,1, 1, mriglm->y->nframes);
      MRIfree(&mriglm->y);
      mriglm->y = mritmp;
    }
    if(DoFisher){
      printf("Computing fisher transform\n");      
      mritmp = MRIfisherTransform(mriglm->y,NULL,NULL);
      if(mritmp == NULL) exit(1);
      MRIfree(&mriglm->y);
      mriglm->y = mritmp;
    }
  }
  else {
    StatTable = LoadStatTable(yFile.c_str());
    if(StatTable == NULL) {
      printf("ERROR: loading y %s as a stat table\n",yFile.c_str());
      exit(1);
    }
    mriglm->y = StatTable->mri;
  }
  if (mriglm->y->type != MRI_FLOAT) {
    printf("INFO: changing y type to float\n");
    mritmp = MRISeqchangeType(mriglm->y,MRI_FLOAT,0,0,0);
    if (mritmp == NULL) {
      printf("ERROR: could change type\n");
      exit(1);
    }
    MRIfree(&mriglm->y);
    mriglm->y = mritmp;
  }

  if(DoFFx){
    // Load yffx var--------------------------------------
    printf("Loading yffxvar from %s\n",yffxvarFile.c_str());
    mriglm->yffxvar = MRIread(yffxvarFile.c_str());
    if(mriglm->yffxvar == NULL) {
      printf("ERROR: loading yffxvar %s\n",yffxvarFile.c_str());
      exit(1);
    }
  }

  if(SubSample){
    printf("Subsampling start=%d delta = %d, nframes = %d \n",
           SubSampStart, SubSampDelta,
           mriglm->y->nframes);
    if( (mriglm->y->nframes % SubSampDelta) != 0){
      printf("ERROR: delta is not an interger divisor of the frames\n");
      exit(1);
    }
    if( SubSampStart > SubSampDelta ){
      printf("ERROR: subsample start > delta\n");
      exit(1);
    }
    mritmp = fMRIsubSample(mriglm->y, SubSampStart, SubSampDelta, -1, NULL);
    if(mritmp == NULL) exit(1);
    MRIfree(&mriglm->y);
    mriglm->y = mritmp;
  }

  if(DoDistance){
    mritmp = fMRIdistance(mriglm->y, NULL);
    MRIfree(&mriglm->y);
    mriglm->y = mritmp;
  }

  nvoxels = mriglm->y->width * mriglm->y->height * mriglm->y->depth;

  // X ---------------------------------------------------------
  //Load global X------------------------------------------------
  if((XFile != NULL) | usedti) {
    if(usedti == 0 || usedti == 3) {
      mriglm->Xg = MatrixReadTxt(XFile, NULL);
      if (mriglm->Xg==NULL) mriglm->Xg = MatlabRead(XFile);
      if (mriglm->Xg==NULL) {
        printf("ERROR: loading X %s\n",XFile);
        printf("Could not load as text or matlab4");
        exit(1);
      }
      if(usedti == 3){
	dti = (DTI*) calloc(sizeof(DTI),1);
	dti->B = mriglm->Xg;
      }
    } 
    else {
      printf("Using DTI\n");
      if(XFile != NULL) dti = DTIstructFromSiemensAscii(XFile);
      else              dti = DTIstructFromBFiles(bvalfile,bvecfile);
      if(dti==NULL) exit(1);
      sprintf(tmpstr,"%s/bvals.dat",GLMDir);
      DTIwriteBValues(dti->bValue, tmpstr);
      //DTIfslBValFile(dti,tmpstr);
      sprintf(tmpstr,"%s/bvecs.dat",GLMDir);
      DTIwriteBVectors(dti->GradDir,tmpstr);
      //DTIfslBVecFile(dti,tmpstr);
      mriglm->Xg = MatrixCopy(dti->B,NULL);
    }
  }

  if (useasl) {
    mriglm->Xg = MatrixConstVal(1.0, mriglm->y->nframes, 3, NULL);
    for (n=0; n < mriglm->y->nframes; n += 2) {
      mriglm->Xg->rptr[n+1][2] = asl1val;
      if(n+2 >= mriglm->y->nframes) break;
      mriglm->Xg->rptr[n+2][2] = asl2val;
    }
    for(n=0; n < mriglm->y->nframes; n ++)
      mriglm->Xg->rptr[n+1][3] = n - mriglm->y->nframes/2.0;
    nContrasts = 3;
    mriglm->glm->ncontrasts = nContrasts;
    mriglm->glm->Cname[0] = "perfusion";
    mriglm->glm->C[0] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[0]->rptr[1][1] = 0;
    mriglm->glm->C[0]->rptr[1][2] = 1;
    mriglm->glm->Cname[1] = "control";
    mriglm->glm->C[1] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[1]->rptr[1][1] = 1;
    mriglm->glm->C[1]->rptr[1][2] = 0;
    mriglm->glm->Cname[2] = "label";
    mriglm->glm->C[2] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[2]->rptr[1][1] = 1;
    mriglm->glm->C[2]->rptr[1][2] = 1;
  }
  if(useqa) {
    // Set up a model with const, linear, and quad
    mriglm->Xg = MatrixConstVal(1.0, mriglm->y->nframes, 3, NULL);
    dtmp = 0;
    for (n=0; n < mriglm->y->nframes; n += 1) {
      mriglm->Xg->rptr[n+1][2] = n - (mriglm->y->nframes-1.0)/2.0;
      mriglm->Xg->rptr[n+1][3] = n*n;
      dtmp += (n*n);
    }
    for (n=0; n < mriglm->y->nframes; n += 1) {
      mriglm->Xg->rptr[n+1][3] -= dtmp/mriglm->y->nframes;
    }
    // Test mean, linear and quad
    // mean is not an important test, but snr = cnr
    nContrasts = 3;
    mriglm->glm->ncontrasts = nContrasts;
    mriglm->glm->Cname[0] = "mean";
    mriglm->glm->C[0] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[0]->rptr[1][1] = 1;
    mriglm->glm->C[0]->rptr[1][2] = 0;
    mriglm->glm->C[0]->rptr[1][3] = 0;
    mriglm->glm->Cname[1] = "linear";
    mriglm->glm->C[1] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[1]->rptr[1][1] = 0;
    mriglm->glm->C[1]->rptr[1][2] = 1;
    mriglm->glm->C[1]->rptr[1][3] = 0;
    mriglm->glm->Cname[2] = "quad";
    mriglm->glm->C[2] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[2]->rptr[1][1] = 0;
    mriglm->glm->C[2]->rptr[1][2] = 0;
    mriglm->glm->C[2]->rptr[1][3] = 1;
  }
  if (fsgd != NULL) {
    printf("INFO: gd2mtx_method is %s\n",gd2mtx_method);
    mriglm->Xg = gdfMatrix(fsgd,gd2mtx_method,NULL);
    if (mriglm->Xg==NULL) exit(1);
  }
  if (OneSampleGroupMean) {
    mriglm->Xg = MatrixConstVal(1.0,mriglm->y->nframes,1,NULL);
  }

  if(NSplits > 0){
    nPerSplit = floor((double)mriglm->y->nframes/NSplits);
    SplitMin = SplitNo*nPerSplit;
    SplitMax = SplitMin + nPerSplit;
    if((SplitNo==NSplits-1) && (SplitMax<mriglm->y->nframes-1)) SplitMax = mriglm->y->nframes-1;
    nExclude = mriglm->y->nframes-(SplitMax-SplitMin);
    printf("nframes=%d, NSplits=%d, SplitNo=%d, nPerSplit=%d, SplitMin=%d, SplitMax=%d, nEx=%d, RandSplit=%d\n",
	   mriglm->y->nframes,NSplits,SplitNo,nPerSplit,SplitMin,SplitMax,nExclude,RandSplit);

    Ex = MatrixConstVal(0,mriglm->y->nframes,1,NULL);
    for(n=0; n<mriglm->y->nframes; n++) Ex->rptr[n+1][1] = n;

    if(RandSplit) {
      /* Make sure to use the same seed for each split group to assure
         that each group is unqiue. Or run it once and get
         glmdir/synthseed.dat, and use that for future splits.*/
      // MatrixRandPermRows() creates a random order of frames to exclude, but 
      // this is not important because it keeps the right order with the design matrix
      MatrixRandPermRows(Ex);
    }

    ExcludeFrames = (int *) calloc(sizeof(int),nExclude);
    m = 0;
    for(n=0; n<mriglm->y->nframes; n++){
      if(n < SplitMin || n >= SplitMax){
	ExcludeFrames[m] = Ex->rptr[n+1][1];
	m++;
      }
    }
  }

  // Randomly create frames to exclude
  if(nRandExclude > 0){
    ExcludeFrames = (int *) calloc(sizeof(int),nRandExclude);
    Ex = MatrixConstVal(0,mriglm->y->nframes,1,NULL);
    for(n=0; n<nRandExclude; n++) Ex->rptr[n+1][1] = 1;
    MatrixRandPermRows(Ex);
    nExclude = 0;
    for(n=0; n<mriglm->y->nframes; n++){
      if(Ex->rptr[n+1][1]){
	ExcludeFrames[nExclude] = n;
	nExclude++;
      }
    }
    MatrixFree(&Ex);
  }
  // Exclude frames from both design matrix and data
  if(ExcludeFrames){
    sprintf(tmpstr,"%s/exclude-frames.dat",GLMDir);
    fp = fopen(tmpstr,"w");
    for(m=0; m<nExclude; m++) fprintf(fp,"%d\n",ExcludeFrames[m]);
    fclose(fp);
    XgNew = MatrixExcludeFrames(mriglm->Xg, ExcludeFrames,nExclude);
    MatrixFree(&mriglm->Xg);
    mriglm->Xg = XgNew; 
    mritmp = fMRIexcludeFrames(mriglm->y, ExcludeFrames,nExclude,NULL);
    MRIfree(&mriglm->y);
    mriglm->y = mritmp;
    if(mriglm->w){
      mritmp = fMRIexcludeFrames(mriglm->w, ExcludeFrames,nExclude,NULL);
      MRIfree(&mriglm->w);
      mriglm->w = mritmp;
    }
  }

  if(frameMaskFile){
    mriglm->FrameMask = MRIread(frameMaskFile);
    if(mriglm->FrameMask == NULL) exit(1);
  }

  // MRTM1 ------------------------------------
  if(DoMRTM1) {
    printf("Performing MRTM1\n"); fflush(stdout);
    mriglm->Xg = MatrixHorCat(RTM_Cr,RTM_intCr,NULL);
    mriglm->npvr = 1;
    printf("Computing integral of input ..."); fflush(stdout);
    mriglm->pvr[0] = fMRIcumTrapZ(mriglm->y,RTM_TimeMin,NULL,NULL);
    printf("done.\n"); fflush(stdout);
    nContrasts = 4;
    mriglm->glm->ncontrasts = nContrasts;
    //------------------------------------------
    mriglm->glm->Cname[0] = "R1";
    mriglm->glm->C[0] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[0]->rptr[1][1] = 1;
    //------------------------------------------
    mriglm->glm->Cname[1] = "k2";
    mriglm->glm->C[1] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[1]->rptr[1][2] = 1;
    //------------------------------------------
    mriglm->glm->Cname[2] = "k2a";
    mriglm->glm->C[2] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[2]->rptr[1][3] = -1;
    //------------------------------------------
    mriglm->glm->Cname[3] = "k2-k2a";
    mriglm->glm->C[3] = MatrixConstVal(0.0, 1, 3, NULL);
    mriglm->glm->C[3]->rptr[1][2] = +1;
    mriglm->glm->C[3]->rptr[1][3] = +1;
    //------------------------------------------
    sprintf(tmpstr,"%s/time.min.dat",GLMDir);
    MatrixWriteTxt(tmpstr, RTM_TimeMin);
  }
  // MRTM2 ------------------------------------
  if(DoMRTM2) {
    printf("Performing MRTM2\n"); fflush(stdout);
    mriglm->Xg = MRTM2_x1;
    mriglm->wg = NULL;
    mriglm->npvr = 1;
    printf("Computing integral of input ..."); fflush(stdout);
    mriglm->pvr[0] = fMRIcumTrapZ(mriglm->y,RTM_TimeMin,NULL,NULL);
    printf("done.\n"); fflush(stdout);
    mriglm->pvr[0] = MRImultiplyConst(mriglm->pvr[0], -1, mriglm->pvr[0]);
    nContrasts = 3;
    mriglm->glm->ncontrasts = nContrasts;
    //------------------------------------------
    mriglm->glm->Cname[0] = "k2";
    mriglm->glm->C[0] = MatrixConstVal(0.0, 1, 2, NULL);
    mriglm->glm->C[0]->rptr[1][1] = 1;
    //------------------------------------------
    mriglm->glm->Cname[1] = "k2a";
    mriglm->glm->C[1] = MatrixConstVal(0.0, 1, 2, NULL);
    mriglm->glm->C[1]->rptr[1][2] = 1;
    //------------------------------------------
    mriglm->glm->Cname[2] = "k2-k2a";
    mriglm->glm->C[2] = MatrixConstVal(0.0, 1, 2, NULL);
    mriglm->glm->C[2]->rptr[1][1] = +1;
    mriglm->glm->C[2]->rptr[1][2] = -1;
    //------------------------------------------
    sprintf(tmpstr,"%s/time.min.dat",GLMDir);
    MatrixWriteTxt(tmpstr, RTM_TimeMin);
  }
  // Logan ------------------------------------
  if(DoLogan) {
    // This is not the most beautiful code:). 
    // Model integralYi = [integralRef Yi]*beta
    //  BPnd = beta[1]-1
    //  The equation is solved only for time points > tstar
    printf("Performing Logan\n"); fflush(stdout);
    mriglm->Xg = RTM_intCr;
    mriglm->npvr = 1;
    printf("Computing integral of input ..."); fflush(stdout);
    mriglm->pvr[0] = mriglm->y;
    mriglm->y = fMRIcumTrapZ(mriglm->y,RTM_TimeMin,NULL,NULL);
    printf("done.\n"); fflush(stdout);
    printf("Loganizing\n"); fflush(stdout);
    MRIloganize(&(mriglm->Xg), &(mriglm->y), &(mriglm->pvr[0]),RTM_TimeMin,Logan_Tstar/60);
    MRIwrite(mriglm->y,"tmp.y.mgh");
    MRIwrite(mriglm->pvr[0],"tmp.pvr.mgh");
    nContrasts = 1;
    mriglm->glm->ncontrasts = nContrasts;
    //------------------------------------------
    mriglm->glm->Cname[0] = "dvr";
    mriglm->glm->C[0] = MatrixConstVal(0.0, 1, 2, NULL);
    mriglm->glm->C[0]->rptr[1][1] = 1;
    //------------------------------------------
    sprintf(tmpstr,"%s/time.min.dat",GLMDir);
    MatrixWriteTxt(tmpstr, RTM_TimeMin);
  }

  if(! DontSave) {
    if(GLMDir != NULL) {
      sprintf(tmpstr,"%s/Xg.dat",GLMDir);
      printf("Saving design matrix to %s\n",tmpstr);
      MatrixWriteTxt(tmpstr, mriglm->Xg);
    }
  }

  // Check the condition of the global matrix -----------------
  printf("Computing normalized matrix\n"); fflush(stdout);
  Xnorm = MatrixNormalizeCol(mriglm->Xg,NULL,NULL);
  Xcond = MatrixNSConditionNumber(Xnorm);
  printf("Normalized matrix condition is %g\n",Xcond);
  if(Xcond > 10000 && ! IllCondOK) {
    printf("Design matrix ------------------\n");
    MatrixPrint(stdout,mriglm->Xg);
    printf("--------------------------------\n");
    printf("ERROR: matrix is ill-conditioned or badly scaled, condno = %g\n",Xcond);
    printf("\n\n--------------------------------\n");
    printf("-------- ERROR: READ THIS -----------------\n");
    printf("--vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv-----\n");
    printf("Possible problem with experimental design:\n");
    printf("Check for duplicate entries and/or lack of range of\n"
           "continuous variables within a class.\n");
    printf("If you seek help with this problem, make sure to send:\n");
    printf("  1. Your command line:\n");
    printf("    %s\n",cmdline);
    printf("  2. The terminal output of this program (ie, everything it prints to the screen)\n");
    if(fsgdfile){
      printf("  3. The FSGD file (%s)\n",fsgdfile);
      printf("  4. The design matrix %s/Xg.dat\n",GLMDir);
    }
    else  printf("  3. The design matrix %s/Xg.dat\n",GLMDir);
    printf("Attempting to diagnose further \n");
    if(fsgd) gdfCheckNPerClass(fsgd);
    err = GLMdiagnoseDesignMatrix(mriglm->Xg);
    // Could check to make sure that there are enough members of the class given 
    // the number of variables.
    if(err == 0) printf(" ... could not determine the cause of the problem\n");
    printf("--^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^----\n\n\n");
    exit(1);
  }
  Xcond = MatrixNSConditionNumber(mriglm->Xg);
  printf("Matrix condition is %g\n",Xcond);
  if(Xcond > 10000 && !ReScaleX){
    printf("\n");
    printf("WARNING: matrix may be badly scaled!\n");
    printf("You might want to re-run with --rescale-x\n");
    printf("\n");
  }
  fflush(stdout);

  // Load Per-Voxel Regressors -----------------------------------
  if(mriglm->npvr > 0 && !DoMRTM1 && !DoMRTM2 && !DoLogan) {
    for (n=0; n < mriglm->npvr; n++) {
      mriglm->pvr[n] = MRIread(pvrFiles[n]);
      if (mriglm->pvr[n] == NULL) exit(1);
      if (mriglm->pvr[n]->nframes != mriglm->Xg->rows) {
        printf("ERROR: dimension mismatch between pvr and X. ");
        printf("pvr has %d frames, X has %d rows.\n",
               mriglm->pvr[n]->nframes,mriglm->Xg->rows);
        printf("PVR %d %s\n",n,pvrFiles[n]);
        exit(1);
      }
    }
  }

  mriglm->mask = NULL;
  // Load the mask file ----------------------------------
  if(maskFile != NULL) {
    mriglm->mask = MRIread(maskFile);
    if(mriglm->mask  == NULL) {
      printf("ERROR: reading mask file %s\n",maskFile);
      exit(1);
    }
    err = MRIdimMismatch(mriglm->mask,mriglm->y,0);
    if(err){
      printf("ERROR: dimension mismatch %d between y and mask\n",err);
      exit(1);
    }
  }
  // Load the label mask file ----------------------------------
  if (labelFile != NULL) {
    clabel = LabelRead(NULL, labelFile);
    if (clabel == NULL) {
      printf("ERROR reading %s\n",labelFile);
      exit(1);
    }
    printf("Found %d points in label.\n",clabel->n_points);
    mriglm->mask = MRISlabel2Mask(surf, clabel, NULL);
    MRIcopyHeader(mriglm->y,mriglm->mask);
    MRIcopyPulseParameters(mriglm->y,mriglm->mask);
    mritmp = mri_reshape(mriglm->mask, mriglm->y->width,
                         mriglm->y->height, mriglm->y->depth, 1);
    MRIfree(&mriglm->mask);
    mriglm->mask = mritmp;
  }
  if(prunemask) {
    printf("Pruning voxels by thr: %e\n", prune_thr);
    if(usedti){
      // NOTE: for DWI volumes
      MRI* firstFrameVol;
      if (!usepruning) 
	prune_thr = 50; // needs to be larger than 0 to get meaningful mask!
      firstFrameVol = MRIcopyFrame(mriglm->y,NULL, 0, 0);
      mriglm->mask = MRIframeBinarize(firstFrameVol,prune_thr,mriglm->mask);
      MRIfree(&firstFrameVol);
    }
    else{
      mriglm->mask = MRIframeBinarize(mriglm->y,prune_thr,mriglm->mask);
    }
  }

  if (mriglm->mask && maskinv)
    MRImaskInvert(mriglm->mask,mriglm->mask);
  if (surf && mriglm->mask)
    MRISremoveRippedFromMask(surf, mriglm->mask, mriglm->mask);

  if (mriglm->mask) {
    nmask = MRInMask(mriglm->mask);
    printf("Found %d voxels in mask\n",nmask);
    if(nmask == 0){
      printf("ERROR: no voxels found in the mask\n");
      if(prunemask)
	printf("  make sure at least one voxel has a non-zero value for each input\n");
      exit(1);
    }
    if (!DontSave) {
      sprintf(tmpstr,"%s/mask.%s",GLMDir,format);
      printf("Saving mask to %s\n",tmpstr);
      MRIwrite(mriglm->mask,tmpstr); 
    }
  } else nmask = nvoxels;
  maskfraction = (double)nmask/nvoxels;

  if(RmSpatialMean){
    printf("Removing spatial mean from each frame\n");
    mritmp = MRIremoveSpatialMean(mriglm->y, mriglm->mask, NULL);
    MRIfree(&mriglm->y);
    mriglm->y = mritmp;
  }

  if(surf != NULL)  {
    searchspace = 0;
    MRI *mritmp = NULL;
    if (mriglm->mask && mriglm->mask->height != surf->nvertices) {
        printf("Reshaping mriglm->mask...\n");
        mritmp = mri_reshape(mriglm->mask, 
                             surf->nvertices,
                             1, 1, mriglm->mask->nframes);
	if(mritmp == NULL) exit(1);
    }
    for(n=0; n < surf->nvertices; n++){
      if(mritmp && MRIgetVoxVal(mritmp,n,0,0,0) < 0.5) continue;
      searchspace += surf->vertices[n].area;
    }
    if (mritmp) MRIfree(&mritmp);
    if (surf->group_avg_surface_area > 0)
      searchspace *= (surf->group_avg_surface_area/surf->total_area);
  } 
  else{
    voxelsize = mriglm->y->xsize * mriglm->y->ysize * mriglm->y->zsize;
    searchspace = nmask * voxelsize;
  }
  printf("search space = %lf\n",searchspace);

  // Check number of frames ----------------------------------
  if (mriglm->y->nframes != mriglm->Xg->rows) {
    printf("ERROR: dimension mismatch between y and X.\n");
    printf("  y has %d inputs, X has %d rows.\n",
           mriglm->y->nframes,mriglm->Xg->rows);
    exit(1);
  }
  // Load the weight file ----------------------------------
  if (wFile != NULL) {
    mriglm->w = MRIread(wFile);
    if (mriglm->w  == NULL) {
      printf("ERROR: reading weight file %s\n",wFile);
      exit(1);
    }
    // Check number of frames
    if (mriglm->y->nframes != mriglm->w->nframes) {
      printf("ERROR: dimension mismatch between y and w.\n");
      printf("  y has %d frames, w has %d frames.\n",
             mriglm->y->nframes,mriglm->w->nframes);
      exit(1);
    }
    // Invert, Sqrt, and Normalize the weights
    mritmp = MRInormWeights(mriglm->w, weightsqrt, weightinv,
                            mriglm->mask, mriglm->w);
    if (mritmp==NULL) exit(1);
    if (weightsqrt || weightinv) {
      sprintf(tmpstr,"%s/wn.%s",GLMDir,format);
      if(!DontSave && !DontSaveWn) MRIwrite(mriglm->w,tmpstr);
    }
  } 
  else if(wgFile != NULL){
    // Global weight to use at every voxel.
    mriglm->wg = MatrixReadTxt(wgFile,NULL);
    if(mriglm->wg==NULL) exit(1);
    if (mriglm->y->nframes != mriglm->wg->rows) {
      printf("ERROR: dimension mismatch between y and wg.\n");
      printf("  y has %d frames, w has %d frames.\n",
             mriglm->y->nframes,mriglm->wg->rows);
      exit(1);
    }
  }
  else if(!DoMRTM1 && !DoMRTM2 && !DoLogan) {
    mriglm->w = NULL;
    mriglm->wg = NULL;
  }

  if(mriglm->wg != NULL){
    sprintf(tmpstr,"%s/wg.mtx",GLMDir);
    MatrixWriteTxt(tmpstr,mriglm->wg);
  }

  if(synth) {
    if(! UseUniform){
      printf("Replacing input data with synthetic white gaussian noise\n");
      MRIrandn(mriglm->y->width,mriglm->y->height,mriglm->y->depth,
	       mriglm->y->nframes,0,1,mriglm->y);
    }
    else {
      printf("Replacing input data with synthetic white uniform noise\n");
      MRIdrand48(mriglm->y->width,mriglm->y->height,mriglm->y->depth,
	       mriglm->y->nframes,UniformMin,UniformMax,mriglm->y);
    }
  }
  if (logflag) {
    printf("Computing natural log of input\n");
    if (usedti) dwi = MRIcopy(mriglm->y,NULL);
    MRIlog(mriglm->y,mriglm->mask,-1,1,mriglm->y);
  }
  if (FWHM > 0 && (!DoSim || !strcmp(csd->simtype,"perm")) ) {
    printf("Smoothing input by fwhm %lf \n",FWHM);
    SmoothSurfOrVol(surf, mriglm->y, mriglm->mask, SmoothLevel);
    printf("   ... done\n");
  }

  // Handle self-regressors
  if (nSelfReg > 0) {
    if (DoSim && !strcmp(csd->simtype,"perm")) {
      printf("ERROR: cannot use --selfreg with perm simulation\n");
      exit(1);
    }

    for (n=0; n<nSelfReg; n++) {
      c = crsSelfReg[n][0];
      r = crsSelfReg[n][1];
      s = crsSelfReg[n][2];
      printf("Self regressor %d     %d %d %d\n",n,c,r,s);
      if (c < 0 || c >= mriglm->y->width ||
          r < 0 || r >= mriglm->y->height ||
          s < 0 || s >= mriglm->y->depth) {
        printf("ERROR: %d self regressor is out of the volume (%d,%d,%d)\n",
               n,c,r,s);
        exit(1);
      }
      MRIglmLoadVox(mriglm,c,r,s,0,NULL);
      GLMxMatrices(mriglm->glm);
      GLMfit(mriglm->glm);
      GLMdump("selfreg",mriglm->glm);
      Mtmp = MatrixHorCat(Xselfreg,mriglm->glm->eres,NULL);
      if (n > 0) MatrixFree(&Xselfreg);
      Xselfreg = Mtmp;
    }

    // Need new mriglm, so alloc and copy the old one
    mriglmtmp = (MRIGLM *) calloc(sizeof(MRIGLM),1);
    mriglmtmp->glm = GLMalloc();
    mriglmtmp->y = mriglm->y;
    mriglmtmp->w = mriglm->w;
    mriglmtmp->Xg = MatrixHorCat(mriglm->Xg,Xselfreg,NULL);
    for (n=0; n < mriglm->npvr; n++) mriglmtmp->pvr[n] = mriglmtmp->pvr[n];
    //MRIglmFree(&mriglm);
    mriglm = mriglmtmp;
  }
  MRIglmNRegTot(mriglm);

  if (DoSim && !strcmp(csd->simtype,"perm")) {
    if (MatrixColsAreNotOrthog(mriglm->Xg)) {
      if (PermForce)
        printf("INFO: design matrix is not orthogonal, but perm forced\n");
      else {
        printf("ERROR: design matrix is not orthogonal, "
               "cannot be used with permutation.\n");
        printf("If this something you really want to do, "
               "run with --perm-force\n");
        exit(1);
      }
    }
  }

  // Load the contrast matrices ---------------------------------
  mriglm->glm->ncontrasts = nContrasts;
  if(nContrasts > 0) {
    for(n=0; n < nContrasts; n++) {
      if (! useasl && ! useqa  && !(fsgd != NULL && fsgd->nContrasts != 0) && !DoMRTM1 && !DoMRTM2 && !DoLogan) {
        // Get its name
        mriglm->glm->Cname[n] =
          fio_basename(CFile[n],".mat"); //strip .mat
        mriglm->glm->Cname[n] =
          fio_basename(mriglm->glm->Cname[n],".mtx"); //strip .mtx
        mriglm->glm->Cname[n] =
          fio_basename(mriglm->glm->Cname[n],".dat"); //strip .dat
        mriglm->glm->Cname[n] =
          fio_basename(mriglm->glm->Cname[n],".con"); //strip .con
        // Read it in
        mriglm->glm->C[n] = MatrixReadTxt(CFile[n], NULL);
        if (mriglm->glm->C[n] == NULL) {
          printf("ERROR: loading C %s\n",CFile[n]);
          exit(1);
        }
	if(Gamma0File[n]){
	  printf("Loading gamma0 file %s\n",Gamma0File[n]);
	  mriglm->glm->gamma0[n] = MatrixReadTxt(Gamma0File[n], NULL);
	  if(mriglm->glm->gamma0[n] == NULL) {
	    printf("ERROR: loading gamma0 %s\n",Gamma0File[n]);
	    exit(1);
	  }
	  mriglm->glm->UseGamma0[n] = 1;
        }
      }
      if(fsgd && fsgd->nContrasts != 0) {
        mriglm->glm->C[n] = MatrixCopy(fsgd->C[n],NULL);
        mriglm->glm->Cname[n] = strcpyalloc(fsgd->ContrastName[n]);
      }
      // Check it's dimension
      if (mriglm->glm->C[n]->cols != mriglm->nregtot) {
        printf("ERROR: dimension mismatch between X and contrast %s",CFile[n]);
        printf("       X has %d cols, C has %d cols\n",
               mriglm->nregtot,mriglm->glm->C[n]->cols);
        exit(1);
      }
      // Check it's condition
      Ct  = MatrixTranspose(mriglm->glm->C[n],NULL);
      CCt = MatrixMultiplyD(mriglm->glm->C[n],Ct,NULL);
      Ccond = MatrixConditionNumber(CCt);
      if (Ccond > 1000) {
        printf("ERROR: contrast %s is ill-conditioned (%g)\n",
               mriglm->glm->Cname[n],Ccond);
        MatrixPrint(stdout,mriglm->glm->C[n]);
        exit(1);
      }
      MatrixFree(&Ct);
      MatrixFree(&CCt);
    }
  }

  if(OneSampleGroupMean) {
    nContrasts = 1;
    mriglm->glm->ncontrasts = nContrasts;
    mriglm->glm->Cname[0] = strcpyalloc("osgm");
    mriglm->glm->C[0] = MatrixConstVal(1.0,1,1,NULL);
  }

  // Check for one-sample group mean with permutation simulation
  if (!strcmp(csd->simtype,"perm")) {
    OneSamplePerm=1;
    if (mriglm->nregtot == 1) {
      for (n=0; n < mriglm->y->nframes; n++) {
        if (mriglm->Xg->rptr[n+1][1] != 1) {
          OneSamplePerm=0;
          break;
        }
      }
    } else OneSamplePerm=0;
    if(OneSamplePerm)
      printf("Design detected as one-sample group mean, "
             "adjusting permutation simulation\n");
  }

  // At this point, y, X, and C have been loaded, now pre-alloc
  GLMallocX(mriglm->glm,mriglm->y->nframes,mriglm->nregtot);
  GLMallocY(mriglm->glm);

  // Check DOF
  if(! DoFFx){
    GLMdof(mriglm->glm);
    printf("DOF = %g\n",mriglm->glm->dof);
    if(mriglm->glm->dof < 1) {
      if(!usedti || mriglm->glm->dof < 0){
	if(! AllowZeroDOF){
	  printf("ERROR: DOF = %g\n",mriglm->glm->dof);
	  exit(1);
	} 
	else {
	  mriglm->glm->AllowZeroDOF = 1;
	  mriglm->glm->dof = 1;
	}
      } else
	printf("WARNING: DOF = %g\n",mriglm->glm->dof);
    }
  }

  // Compute Contrast-related matrices
  if(DoPCC) mriglm->glm->X = mriglm->Xg;
  GLMcMatrices(mriglm->glm);

  if (pcaSave) {
    if (npca < 0) npca = mriglm->y->nframes;
    if (npca > mriglm->y->nframes) {
      printf("ERROR: npca = %d, max can be %d\n",npca,mriglm->y->nframes);
      exit(1);
    }
  }

  // Dump a voxel
  if (voxdumpflag) {
    sprintf(voxdumpdir,"%s/voxdump-%d-%d-%d",GLMDir,
            voxdump[0],voxdump[1],voxdump[2]);
    printf("Dumping voxel %d %d %d to %s\n",
           voxdump[0],voxdump[1],voxdump[2],voxdumpdir);
    MRIglmLoadVox(mriglm,voxdump[0],voxdump[1],voxdump[2],0,NULL);
    GLMxMatrices(mriglm->glm);
    GLMfit(mriglm->glm);
    GLMtest(mriglm->glm);
    GLMdump(voxdumpdir,mriglm->glm);
    if (mriglm->w) {
      wvect = MRItoVector(mriglm->w,voxdump[0],voxdump[1],voxdump[2],NULL);
      sprintf(tmpstr,"%s/w.dat",voxdumpdir);
      MatrixWriteTxt(tmpstr,wvect);
    }
    exit(0);
  }

  if(UseStatTable){
    OutStatTable = AllocStatTable(StatTable->ncols,nContrasts);
    OutStatTable->measure = strcpyalloc(StatTable->measure);
    for(n=0; n < StatTable->ncols; n++)
      OutStatTable->rownames[n] = strcpyalloc(StatTable->colnames[n]);
    for(n=0; n < nContrasts; n++)
      OutStatTable->colnames[n] = strcpyalloc(mriglm->glm->Cname[n]);
    GammaStatTable = AllocStatTable(StatTable->ncols,nContrasts);
    GammaStatTable->measure = strcpyalloc(StatTable->measure);
    for(n=0; n < StatTable->ncols; n++)
      GammaStatTable->rownames[n] = strcpyalloc(StatTable->colnames[n]);
    for(n=0; n < nContrasts; n++)
      GammaStatTable->colnames[n] = strcpyalloc(mriglm->glm->Cname[n]);
  }

  // Don't do sim --------------------------------------------------------
  if (!DoSim) {
    // Now do the estimation and testing
    mytimer.reset() ;

    if (VarFWHM > 0) {
      printf("Starting fit\n");  fflush(stdout);
      MRIglmFit(mriglm);
      printf("Variance smoothing\n");
      SmoothSurfOrVol(surf, mriglm->rvar, mriglm->mask, VarSmoothLevel);
      printf("Starting test\n");   fflush(stdout);
      MRIglmTest(mriglm);
    } else {
      printf("Starting fit and test\n");   fflush(stdout);
      MRIglmFitAndTest(mriglm);
    }
    msecFitTime = mytimer.milliseconds();
    printf("Fit completed in %g minutes\n",msecFitTime/(1000*60.0));  fflush(stdout);
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  if (DoSim) {
    csd->seed = SynthSeed;
    if (surf != NULL) {
      strcpy(csd->anattype,"surface");
      strcpy(csd->subject,subject);
      strcpy(csd->hemi,hemi);
    } 
    else  strcpy(csd->anattype,"volume");
    csd->searchspace = searchspace;
    csd->nreps = nsim;
    CSDallocData(csd);
    if (!strcmp(csd->simtype,"mc-z")) {
      rfs = RFspecInit(SynthSeed,NULL);
      rfs->name = strcpyalloc("gaussian");
      rfs->params[0] = 0;
      rfs->params[1] = 1;
      //z = MRIalloc(mriglm->y->width,mriglm->y->height,mriglm->y->depth,MRI_FLOAT);
      z = MRIcloneBySpace(mriglm->y,MRI_FLOAT,1);
      zabs = MRIcloneBySpace(mriglm->y,MRI_FLOAT,1);
    }
    if (!strcmp(csd->simtype,"mc-t")) {
      rfs = RFspecInit(SynthSeed,NULL);
      rfs->name = strcpyalloc("t");
      rfs->params[0] = mriglm->glm->dof;
      z = MRIcloneBySpace(mriglm->y,MRI_FLOAT,1);
      zabs = MRIcloneBySpace(mriglm->y,MRI_FLOAT,1);
    }
    printf("thresh = %g, threshadj = %g \n",csd->thresh,csd->thresh-log10(2.0));

    if(!DoSimThreshLoop){
      nThreshList = 1;
      ThreshList[0] = csd->thresh;
      nSignList = 1;
      SignList[0] = tSimSign;
      DoSimThreshLoop = 1;
    }

    if(DoSimThreshLoop){
      for(nthThresh = 0; nthThresh < nThreshList; nthThresh++){
	for(nthSign = 0; nthSign < nSignList; nthSign++){
	  for (n=0; n < mriglm->glm->ncontrasts; n++) {
	    csdList[nthThresh][nthSign][n] = CSDcopy(csd,NULL);
	    csdList[nthThresh][nthSign][n]->thresh = ThreshList[nthThresh];
	    csdList[nthThresh][nthSign][n]->threshsign = SignList[nthSign];
	    csdList[nthThresh][nthSign][n]->seed = csd->seed;
	  }
	}
      }
    }

    printf("\n\nStarting simulation sim over %d trials\n",nsim);
    mytimer.reset() ;
    for (nthsim=0; nthsim < nsim; nthsim++) {
      msecFitTime = mytimer.milliseconds();
      if(debug) printf("%d/%d t=%g ---------------------------------\n",
             nthsim+1,nsim,msecFitTime/(1000*60.0));

      if (!strcmp(csd->simtype,"mc-full")) {
	if(! UseUniform)
	  MRIrandn(mriglm->y->width,mriglm->y->height,mriglm->y->depth,
		   mriglm->y->nframes,0,1,mriglm->y);
	else
	  MRIdrand48(mriglm->y->width,mriglm->y->height,mriglm->y->depth,
		  mriglm->y->nframes,UniformMin,UniformMax,mriglm->y);
        if(logflag) MRIlog(mriglm->y,mriglm->mask,-1,1,mriglm->y);
        if(FWHM > 0)
          SmoothSurfOrVol(surf, mriglm->y, mriglm->mask, SmoothLevel);
      }
      if (!strcmp(csd->simtype,"perm")) {
        if (!OneSamplePerm) MatrixRandPermRows(mriglm->Xg);
        else {
          for (n=0; n < mriglm->y->nframes; n++) {
            if (drand48() > 0.5) m = +1;
            else                m = -1;
            mriglm->Xg->rptr[n+1][1] = m;
          }
          //MatrixPrint(stdout,mriglm->Xg);
        }
      }

      // Variance smoothing
      if (!strcmp(csd->simtype,"mc-full") || !strcmp(csd->simtype,"perm")) {
        // If variance smoothing, then need to test and fit separately
        if (VarFWHM > 0) {
          if(!DoSim) printf("Starting fit\n");
          MRIglmFit(mriglm);
          if(!DoSim) printf("Variance smoothing\n");
          SmoothSurfOrVol(surf, mriglm->rvar, mriglm->mask, VarSmoothLevel);
          if(!DoSim) printf("Starting test\n");
          MRIglmTest(mriglm);
        }
	else {
          if(!DoSim) printf("Starting fit and test\n");
          MRIglmFitAndTest(mriglm);
	  // If using permutation with non-stationary correction, compute fwhmmap here
	  if(!strcmp(csd->simtype,"perm") && PermNonStatCor) {
	    if(ar1)     MRIfree(&ar1);
	    if(fwhmmap) MRIfree(&fwhmmap);
	    ar1 = MRISar1(surf, mriglm->eres, mriglm->mask, NULL);
	    fwhmmap = MRISfwhmFromAR1Map(surf, mriglm->mask, ar1);
	    // MRISsmoothMRI(surf, fwhmmap, SmthLevel, mriglm->mask, fwhmmap)
	  }
        }
      }

      for(nthThresh = 0; nthThresh < nThreshList; nthThresh++){
	for(nthSign = 0; nthSign < nSignList; nthSign++){
	  // Go through each contrast.
	  for (n=0; n < mriglm->glm->ncontrasts; n++) {
	    if(DoSimThreshLoop) {
	      csd = csdList[nthThresh][nthSign][n];
	      tSimSign = SignList[nthSign];
	    }
	    if(debug) printf("%2d %d %5.1f  %d %2d %5.1f\n",nthsim,nthThresh,
			     csd->thresh,nthSign,tSimSign,mytimer.seconds());

	    // Change sign to abs for F-tests
	    csd->threshsign = tSimSign;
	    if(mriglm->glm->C[n]->rows > 1) csd->threshsign = 0;
	    
	    // Adjust threshold for one- or two-sided
	    if(csd->threshsign == 0) threshadj = csd->thresh;
	    else threshadj = csd->thresh - log10(2.0); // one-sided test

	    if (!strcmp(csd->simtype,"mc-full") || !strcmp(csd->simtype,"perm")) {
	      sig = MRIlog10(mriglm->p[n],NULL,sig,1);
	      // If test is not ABS then apply the sign
	      if(csd->threshsign != 0) MRIsetSign(sig,mriglm->gamma[n],0);
	      sigmax = MRIframeMax(sig,0,mriglm->mask,csd->threshsign,
				   &cmax,&rmax,&smax);
	      // Get Fmax at sig max 
	      Fmax = MRIgetVoxVal(mriglm->F[n],cmax,rmax,smax,0);
	      if(csd->threshsign != 0) Fmax = Fmax*SIGN(sigmax);
	    } 
	    else {
	      // mc-z or mc-t: synth z-field, smooth, rescale,
	      // compute p, compute sig
	      // This should do the same thing as AFNI's AlphaSim
	      // Synth and rescale without the mask, otherwise smoothing
	      // smears the 0s into the mask area. Also, the stuff outisde
	      // the mask area wont get zeroed.
	      if(nthThresh == 0 && nthSign == 0) {
		RFsynth(z,rfs,mriglm->mask); // z or t, as needed
		if (SmoothLevel > 0) {
		  SmoothSurfOrVol(surf, z, mriglm->mask, SmoothLevel);
		  if(DiagCluster) {
		    sprintf(tmpstr,"./%s-zsm0.%s",mriglm->glm->Cname[n],format);
		    printf("Saving z into %s\n",tmpstr);
		    MRIwrite(z,tmpstr);
		    // Exits below
		  }
		  RFrescale(z,rfs,mriglm->mask,z);
		}
	      }
	      if(DiagCluster) {
		sprintf(tmpstr,"./%s-zsm1.%s",mriglm->glm->Cname[n],format);
		printf("Saving z into %s\n",tmpstr);
		MRIwrite(z,tmpstr);
		// Exits below
	      }
	      // Slightly tortured way to get the right p-values because
	      //   RFstat2P() computes one-sided, but I handle sidedness
	      //   during thresholding.
	      // First, use zabs to get a two-sided pval bet 0 and 0.5
	      zabs = MRIabs(z,zabs);
	      mriglm->p[n] = RFstat2P(zabs,rfs,mriglm->mask,0,mriglm->p[n]);
	      // Next, mult pvals by 2 to get two-sided bet 0 and 1
	      MRIscalarMul(mriglm->p[n],mriglm->p[n],2);
	      // sig = -log10(p)
	      sig = MRIlog10(mriglm->p[n],NULL,sig,1);
	      // If test is not ABS then apply the sign
	      if(csd->threshsign != 0) MRIsetSign(sig,z,0);

	      sigmax = MRIframeMax(sig,0,mriglm->mask,csd->threshsign,
				   &cmax,&rmax,&smax);
	      Fmax = MRIgetVoxVal(z,cmax,rmax,smax,0);
	      if(csd->threshsign == 0) Fmax = fabs(Fmax);
	    }
	    if(mriglm->mask) MRImask(sig,mriglm->mask,sig,0.0,0.0);

	    if(surf) {
	      // surface clustering -------------
	      MRIScopyMRI(surf, sig, 0, "val");
	      if(debug || Gdiag_no > 0) printf("Clustering on surface %lf\n",
					       mytimer.seconds());
	      SurfClustList = sclustMapSurfClusters(surf,threshadj,-1,csd->threshsign,
						    0,&nClusters,NULL,fwhmmap);
	      csize = sclustMaxClusterArea(SurfClustList, nClusters);
	    } 
	    else {
	      // volume clustering -------------
	      if (debug) printf("Clustering on volume\n");
	      VolClustList = clustGetClusters(sig, 0, threshadj,-1,csd->threshsign,0,
					      mriglm->mask, &nClusters, NULL);
	      csize = voxelsize*clustMaxClusterCount(VolClustList,nClusters);
	      if (Gdiag_no > 0) clustDumpSummary(stdout,VolClustList,nClusters);
	      clustFreeClusterList(&VolClustList,nClusters);
	    }
	    if(debug) printf("%s %d nc=%d  maxcsize=%g  sigmax=%g  Fmax=%g\n",
			     mriglm->glm->Cname[n],nthsim,nClusters,csize,sigmax,Fmax);

	    // Re-write the full CSD file each time. Should not take that
	    // long and assures output can be used immediately regardless
	    // of whether the job terminated properly or not
	    strcpy(csd->contrast,mriglm->glm->Cname[n]);
	    if(DoSimThreshLoop && (nThreshList > 1 || nSignList > 1) ){
	      if(round(csd->threshsign) ==  0) tmpstr2 = "abs"; 
	      if(round(csd->threshsign) == +1) tmpstr2 = "pos"; 
	      if(round(csd->threshsign) == -1) tmpstr2 = "neg"; 
	      //sprintf(tmpstr,"%s-%s.th%04d.%s.csd",simbase,mriglm->glm->Cname[n],
	      //      (int)round(csd->thresh*100),tmpstr2);
	      sprintf(tmpstr,"%s.th%02d.%s.j001-%s.csd",simbase,
		      (int)round(csd->thresh*10),tmpstr2,mriglm->glm->Cname[n]);
	    }
	    else
	      sprintf(tmpstr,"%s-%s.csd",simbase,mriglm->glm->Cname[n]);
	    if(debug) printf("csd %s \n",tmpstr);
	    fflush(stdout);
	    fp = fopen(tmpstr,"w");
	    if (fp == NULL) {
	      printf("ERROR: opening %s\n",tmpstr);
	      exit(1);
	    }
	    fprintf(fp,"# ClusterSimulationData 2\n");
	    fprintf(fp,"# mri_glmfit simulation sim\n");
	    fprintf(fp,"# hostname %s\n",uts.nodename);
	    fprintf(fp,"# machine  %s\n",uts.machine);
	    fprintf(fp,"# runtime_min %g\n",msecFitTime/(1000*60.0));
	    fprintf(fp,"# FixVertexAreaFlag %d\n",MRISgetFixVertexAreaValue());
	    if (mriglm->mask) fprintf(fp,"# masking 1\n");
	    else             fprintf(fp,"# masking 0\n");
	    fprintf(fp,"# num_dof %d\n",mriglm->glm->C[n]->rows);
	    fprintf(fp,"# den_dof %g\n",mriglm->glm->dof);
	    fprintf(fp,"# SmoothLevel %g\n",SmoothLevel);
	    csd->nreps = nthsim+1;
	    csd->nClusters[nthsim] = nClusters;
	    csd->MaxClusterSize[nthsim] = csize;
	    csd->MaxSig[nthsim] = sigmax;
	    csd->MaxStat[nthsim] = Fmax;
	    CSDprint(fp, csd);
	    fclose(fp);
	    if(debug) CSDprint(stdout, csd);

	    if(DiagCluster) {
	      sprintf(tmpstr,"./%s-sig.%s",mriglm->glm->Cname[n],format);
	      printf("Saving sig into %s and exiting ... \n",tmpstr);
	      MRIwrite(sig,tmpstr);
	      exit(1);
	    }
	    free(SurfClustList);
	  } // contrasts
	} // sign list
      } // thresh list
      //MRIfree(&sig);

    }// simulation loop
    if(SimDoneFile){
      fp = fopen(SimDoneFile,"w");
      fclose(fp);
    }
    msecFitTime = mytimer.milliseconds();
    printf("mri_glmfit simulation done %g\n\n\n",msecFitTime/(1000*60.0));
    exit(0);
  }
  //--------------------------------------------------------------------------

  if (MaxVoxBase != NULL) {
    for (n=0; n < mriglm->glm->ncontrasts; n++) {
      sig    = MRIlog10(mriglm->p[n],NULL,sig,1);
      sigmax = MRIframeMax(sig,0,mriglm->mask,0,&cmax,&rmax,&smax);
      Fmax = MRIgetVoxVal(mriglm->F[n],cmax,rmax,smax,0);
      sprintf(tmpstr,"%s-%s",MaxVoxBase,mriglm->glm->Cname[n]);
      fp = fopen(tmpstr,"a");
      fprintf(fp,"%e  %e    %d %d %d     %d\n",
              sigmax,Fmax,cmax,rmax,smax,SynthSeed);
      fclose(fp);
      MRIfree(&sig);
    }
  }

  if (DontSave) exit(0);

  if(ComputeFWHM) {
    // Compute fwhm of residual
    if (surf != NULL) {
      printf("Computing spatial AR1 on surface\n");
      ar1 = MRISar1(surf, mriglm->eres, mriglm->mask, NULL);
      sprintf(tmpstr,"%s/sar1.%s",GLMDir,format);
      MRIwrite(ar1,tmpstr);
      RFglobalStats(ar1, mriglm->mask, &ar1mn, &ar1std, &ar1max);
      eresfwhm = MRISfwhmFromAR1(surf, ar1mn);
      eresgstd = eresfwhm/sqrt(log(256.0));
      printf("Residual: ar1mn=%lf, ar1std=%lf, gstd=%lf, fwhm=%lf\n",
             ar1mn,ar1std,eresgstd,eresfwhm);
      if(SaveFWHMMap){
	printf("Computing map of FWHM\n");
	fwhmmap = MRISfwhmFromAR1Map(surf, mriglm->mask, ar1);
	sprintf(tmpstr,"%s/fwhm.%s",GLMDir,format);
	MRIwrite(fwhmmap,tmpstr);
	MRIfree(&fwhmmap);
      }
      MRIfree(&ar1);
    } else {
      printf("Computing spatial AR1 in volume.\n");
      ar1 = fMRIspatialAR1(mriglm->eres, mriglm->mask, NULL);
      if (ar1 == NULL) exit(1);
      sprintf(tmpstr,"%s/sar1.%s",GLMDir,format);
      MRIwrite(ar1,tmpstr);
      fMRIspatialAR1Mean(ar1, mriglm->mask, &car1mn, &rar1mn, &sar1mn);
      cfwhm = RFar1ToFWHM(car1mn, mriglm->eres->xsize);
      rfwhm = RFar1ToFWHM(rar1mn, mriglm->eres->ysize);
      sfwhm = RFar1ToFWHM(sar1mn, mriglm->eres->zsize);
      eresfwhm = sqrt((cfwhm*cfwhm + rfwhm*rfwhm + sfwhm*sfwhm)/3.0);
      printf("Residual: ar1mn = (%lf,%lf,%lf) fwhm = (%lf,%lf,%lf) %lf\n",
             car1mn,rar1mn,sar1mn,cfwhm,rfwhm,sfwhm,eresfwhm);
      MRIfree(&ar1);
    }
    sprintf(tmpstr,"%s/fwhm.dat",GLMDir);
    fp = fopen(tmpstr,"w");
    fprintf(fp,"%f\n",eresfwhm);
    fclose(fp);
  }

  if(DoTemporalAR1){
    printf("Computing temporal AR1\n");
    tar1 = fMRItemporalAR1(mriglm->eres,
                           mriglm->beta->nframes,
                           mriglm->mask,
                           NULL);
    sprintf(tmpstr,"%s/tar1.%s",GLMDir,format);
    MRIwrite(tar1,tmpstr);
  }

  // Save estimation results
  printf("Writing results\n");
  MRIwrite(mriglm->beta,betaFile);
  MRIwrite(mriglm->rvar,rvarFile);

  rstd = MRIsqrt(mriglm->rvar,NULL);
  sprintf(tmpstr,"%s/rstd.%s",GLMDir,format);
  MRIwrite(rstd,tmpstr);

  if(mriglm->yhatsave) MRIwrite(mriglm->yhat,yhatFile);
  if(mriglm->condsave) MRIwrite(mriglm->cond,condFile);
  if(eresFile) MRIwrite(mriglm->eres,eresFile);
  if(eresSCMFile){
    printf("Computing residual spatial correlation matrix\n");
    mritmp = fMRIspatialCorMatrix(mriglm->eres);
    if(mritmp == NULL) exit(1);
    MRIwrite(mritmp,eresSCMFile);
    MRIfree(&mritmp);
  }
  if(yOutFile.size() != 0){
    err = MRIwrite(mriglm->y,yOutFile.c_str());
    if(err) {
      exit(1);
    }
  }

  sprintf(tmpstr,"%s/Xg.dat",GLMDir);
  MatrixWriteTxt(tmpstr, mriglm->Xg);

  sprintf(tmpstr,"%s/dof.dat",GLMDir);
  fp = fopen(tmpstr,"w");  
  if(DoFFx) fprintf(fp,"%d\n",(int)mriglm->ffxdof);
  else      fprintf(fp,"%d\n",(int)mriglm->glm->dof);
  fclose(fp);

  if(useqa){
    // Compute FSNR
    float fsnrmin, fsnrmax, fsnrrange, fsnrmean, fsnrstd;
    printf("Computing FSNR\n");
    fsnr = MRIdivide(mriglm->gamma[0],rstd,NULL);
    sprintf(tmpstr,"%s/fsnr.%s",GLMDir,format);
    MRIwrite(fsnr,tmpstr);
    MRIsegStats(mriglm->mask, 1, fsnr, 0,
		&fsnrmin, &fsnrmax, &fsnrrange, &fsnrmean, &fsnrstd);
    sprintf(tmpstr,"%s/fsnr.dat",GLMDir);
    fp = fopen(tmpstr,"w");  
    fprintf(fp,"%f %f\n",fsnrmean,fsnrstd);
    fclose(fp);
    // Write out mean
    sprintf(tmpstr,"%s/mean.%s",GLMDir,format);
    MRIwrite(mriglm->gamma[0],tmpstr);
  }

  // Save the contrast results
  for (n=0; n < mriglm->glm->ncontrasts; n++) {
    printf("  %s\n",mriglm->glm->Cname[n]);

    // Create output directory for contrast
    sprintf(tmpstr,"%s/%s",GLMDir,mriglm->glm->Cname[n]);
    err = mkdir(tmpstr,0777);
    if (err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n",GLMDir);
      perror(NULL);
      return(1);
    }

    // Dump contrast matrix
    sprintf(tmpstr,"%s/%s/C.dat",GLMDir,mriglm->glm->Cname[n]);
    MatrixWriteTxt(tmpstr, mriglm->glm->C[n]);

    // Save gamma
    sprintf(tmpstr,"%s/%s/gamma.%s",GLMDir,mriglm->glm->Cname[n],format);
    MRIwrite(mriglm->gamma[n],tmpstr);

    if(mriglm->glm->C[n]->rows == 1){
      // Save gammavar
      sprintf(tmpstr,"%s/%s/gammavar.%s",GLMDir,mriglm->glm->Cname[n],format);
      MRIwrite(mriglm->gammaVar[n],tmpstr);
    }

    // Save F
    sprintf(tmpstr,"%s/%s/F.%s",GLMDir,mriglm->glm->Cname[n],format);
    MRIwrite(mriglm->F[n],tmpstr);

    // Compute and Save -log10 p-values
    sig=MRIlog10(mriglm->p[n],NULL,sig,1);
    if (mriglm->mask) MRImask(sig,mriglm->mask,sig,0.0,0.0);

    if(UseStatTable){
      for(m=0; m < OutStatTable->nrows; m++){
	//OutStatTable->data[m][n] = MRIgetVoxVal(mriglm->p[n],m,0,0,0);
	OutStatTable->data[m][n] = MRIgetVoxVal(sig,m,0,0,0);
	if(mriglm->glm->C[n]->rows == 1)
	  OutStatTable->data[m][n] *= SIGN(MRIgetVoxVal(mriglm->gamma[n],m,0,0,0));
	GammaStatTable->data[m][n] = MRIgetVoxVal(mriglm->gamma[n],m,0,0,0);
      }
    }

    // Compute and Save CNR and PCC
    if(mriglm->glm->C[n]->rows == 1){
      cnr = MRIdivide(mriglm->gamma[n], rstd, NULL) ;
      sprintf(tmpstr,"%s/%s/cnr.%s",GLMDir,mriglm->glm->Cname[n],format);
      MRIwrite(cnr,tmpstr);
      MRIfree(&cnr);
      if(mriglm->glm->Dt[n] != NULL){
	// Write out the pcc
	sprintf(tmpstr,"%s/%s/pcc.%s",GLMDir,mriglm->glm->Cname[n],format);
	MRIwrite(mriglm->pcc[n],tmpstr);
      }
    }

    // If it is t-test (ie, one row) then apply the sign
    if (mriglm->glm->C[n]->rows == 1) MRIsetSign(sig,mriglm->gamma[n],0);

    // Write out the sig
    sprintf(tmpstr,"%s/%s/sig.%s",GLMDir,mriglm->glm->Cname[n],format);
    MRIwrite(sig,tmpstr);

    // Write out the z
    sprintf(tmpstr,"%s/%s/z.%s",GLMDir,mriglm->glm->Cname[n],format);
    MRIwrite(mriglm->z[n],tmpstr);

    // Find and save the max sig
    sigmax = MRIframeMax(sig,0,mriglm->mask,0,&cmax,&rmax,&smax);
    Fmax = MRIgetVoxVal(mriglm->F[n],cmax,rmax,smax,0);
    printf("    maxvox sig=%g  F=%g  at  index %d %d %d    seed=%d\n",
           sigmax,Fmax,cmax,rmax,smax,SynthSeed);

    sprintf(tmpstr,"%s/%s/maxvox.dat",GLMDir,mriglm->glm->Cname[n]);
    fp = fopen(tmpstr,"w");
    fprintf(fp,"%e  %e    %d %d %d     %d\n",
            sigmax,Fmax,cmax,rmax,smax,SynthSeed);
    fclose(fp);

    MRIfree(&sig);

    if(mriglm->npvr == 0){
      eff = GLMEfficiency(mriglm->Xg,mriglm->glm->C[n]);
      sprintf(tmpstr,"%s/%s/efficiency.dat",GLMDir,mriglm->glm->Cname[n]);
      fp = fopen(tmpstr,"w");
      fprintf(fp,"%g\n",eff);
      fclose(fp);
    }
  } // contrasts
  
  if(UseStatTable){
    PrintStatTable(stdout, OutStatTable);
    sprintf(tmpstr,"%s/sig.table.dat",GLMDir);
    WriteStatTable(tmpstr, OutStatTable);
    sprintf(tmpstr,"%s/input.table.dat",GLMDir);
    WriteStatTable(tmpstr, StatTable);
    sprintf(tmpstr,"%s/gamma.table.dat",GLMDir);
    WriteStatTable(tmpstr, GammaStatTable);
  }

  if (usedti) {
    printf("Saving DTI Analysis\n");
    lowb = DTIbeta2LowB(mriglm->beta, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/lowb.%s",GLMDir,format);
    MRIwrite(lowb,tmpstr);

    tensor = DTIbeta2Tensor(mriglm->beta, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/tensor.%s",GLMDir,format);
    MRIwrite(tensor,tmpstr);

    evals=NULL;
    evec1=NULL;
    evec2=NULL;
    evec3=NULL;
    DTItensor2Eig(tensor, mriglm->mask, &evals, &evec1, &evec2, &evec3);
    sprintf(tmpstr,"%s/eigvals.%s",GLMDir,format);
    MRIwrite(evals,tmpstr);
    sprintf(tmpstr,"%s/eigvec1.%s",GLMDir,format);
    MRIwrite(evec1,tmpstr);
    sprintf(tmpstr,"%s/eigvec2.%s",GLMDir,format);
    MRIwrite(evec2,tmpstr);
    sprintf(tmpstr,"%s/eigvec3.%s",GLMDir,format);
    MRIwrite(evec3,tmpstr);

    printf("Computing fa\n");
    fa = DTIeigvals2FA(evals, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/fa.%s",GLMDir,format);
    MRIwrite(fa,tmpstr);

    printf("Computing ra\n");
    ra = DTIeigvals2RA(evals, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/ra.%s",GLMDir,format);
    MRIwrite(ra,tmpstr);

    printf("Computing vr\n");
    vr = DTIeigvals2VR(evals, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/vr.%s",GLMDir,format);
    MRIwrite(vr,tmpstr);

    printf("Computing radial diffusivity\n");
    vr = DTIradialDiffusivity(evals, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/radialdiff.%s",GLMDir,format);
    MRIwrite(vr,tmpstr);

    printf("Computing adc\n");
    adc = DTItensor2ADC(tensor, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/adc.%s",GLMDir,format);
    MRIwrite(adc,tmpstr);

    printf("Computing ivc\n");
    ivc = DTIivc(evec1, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/ivc.%s",GLMDir,format);
    MRIwrite(ivc,tmpstr);

    if(mriglm->glm->dof > 0){
      printf("Computing dwisynth\n");
      dwisynth = DTIsynthDWI(dti->B, mriglm->beta, mriglm->mask, NULL);
      //sprintf(tmpstr,"%s/dwisynth.%s",GLMDir,format);
      //MRIwrite(dwisynth,tmpstr);

      printf("Computing dwires\n");
      dwires = MRIsum(dwi, dwisynth, 1, -1, mriglm->mask, NULL);
      if(eresSave){
	sprintf(tmpstr,"%s/dwires.%s",GLMDir,format);
	MRIwrite(dwires,tmpstr);
      }

      printf("Computing dwi rvar\n");
      dwirvar = fMRIcovariance(dwires, 0, mriglm->beta->nframes, 0, NULL);
      sprintf(tmpstr,"%s/dwirvar.%s",GLMDir,format);
      MRIwrite(dwirvar,tmpstr);

      MRIfree(&dwisynth);
    }

    MRIfree(&lowb);
    MRIfree(&tensor);
    MRIfree(&evals);
    MRIfree(&evec1);
    MRIfree(&evec2);
    MRIfree(&evec3);
    MRIfree(&fa);
    MRIfree(&ra);
    MRIfree(&vr);
    MRIfree(&adc);
  }

  if(DoMRTM1){
    MRI *sig1, *sig2, *sig3, *c123;
    printf("Computing binding potentials\n");
    mritmp = BindingPotential(mriglm->gamma[1],mriglm->gamma[2], mriglm->mask, NULL);
    sprintf(tmpstr,"%s/bp.%s",GLMDir,format);
    err = MRIwrite(mritmp,tmpstr);
    if(err) exit(1);
    MRIfree(&mritmp);
    printf("Computing conjunction of k2, k2a, and k2-k2a\n");
    sig1 = MRIlog10(mriglm->p[1],NULL,NULL,1); // k2
    MRIsetSign(sig1,mriglm->gamma[1],0);
    //if(mriglm->mask) MRImask(sig1,mriglm->mask,sig1,0.0,0.0);
    sig2 = MRIlog10(mriglm->p[2],NULL,NULL,1); // k2a
    MRIsetSign(sig2,mriglm->gamma[2],0);
    //if(mriglm->mask) MRImask(sig2,mriglm->mask,sig3,0.0,0.0);
    sig3 = MRIlog10(mriglm->p[3],NULL,NULL,1); // k2-k2a
    MRIsetSign(sig3,mriglm->gamma[3],0);
    //if(mriglm->mask) MRImask(sig3,mriglm->mask,sig3,0.0,0.0);
    c123 = MRIconjunct3(sig1, sig2, sig3, mriglm->mask, NULL);
    sprintf(tmpstr,"%s/bp.kconjunction.%s",GLMDir,format);
    err = MRIwrite(c123,tmpstr);
    if(err) exit(1);
    MRIfree(&sig1);MRIfree(&sig2);MRIfree(&sig3);
    MRIfree(&c123);

    MRI *k2prime = MRIdivide(mriglm->gamma[1],mriglm->gamma[0],NULL);
    sprintf(tmpstr,"%s/k2prime.%s",GLMDir,format);
    err = MRIwrite(k2prime,tmpstr);
    if(err) exit(1);
    if(nmask == 1){
      sprintf(tmpstr,"%s/k2prime.dat",GLMDir);
      FILE *fp = fopen(tmpstr,"w");
      if(fp == NULL) exit(1);
      fprintf(fp,"%20.15f\n",MRIgetVoxVal(k2prime,0,0,0,0));
      fclose(fp);
    }
    MRIfree(&k2prime);
  }

  if(DoMRTM2){
    printf("Computing binding potentials\n");
    mritmp = BindingPotential(mriglm->gamma[0],mriglm->gamma[1], mriglm->mask, NULL);
    sprintf(tmpstr,"%s/bp.%s",GLMDir,format);
    err = MRIwrite(mritmp,tmpstr);
    if(err) exit(1);
    MRIfree(&mritmp);
  }
  if(DoLogan){
    printf("Computing binding potentials\n");
    mritmp = MRIcloneBySpace(mriglm->y,MRI_FLOAT,1);
    for(c=0; c < mriglm->y->width; c++) {
      for(r=0; r < mriglm->y->height; r++) {
	for(s=0; s < mriglm->y->depth; s++) {
	  if(mriglm->mask && MRIgetVoxVal(mriglm->mask,c,r,s,0)<0.5) continue;
          double v = MRIgetVoxVal(mriglm->beta,c,r,s,0);
	  MRIsetVoxVal(mritmp,c,r,s,0,v-1);
	}//s
      }//r
    }//s
    sprintf(tmpstr,"%s/bp.%s",GLMDir,format);
    err = MRIwrite(mritmp,tmpstr);
    if(err) exit(1);
    MRIfree(&mritmp);
  }

  sprintf(tmpstr,"%s/X.mat",GLMDir);
  MatlabWrite(mriglm->Xg,tmpstr,"X");
  if(0){
    // Does not work properly. Image values are not right. Rescale?
    I = ImageFromMatrix(mriglm->Xg, NULL);
    sprintf(tmpstr,"%s/X.rgb",GLMDir);
    printf("Writing rgb of design matrix to %s\n",tmpstr);
    ImageWrite(I, tmpstr);
    ImageFree(&I);
  }

  // --------- Save FSGDF stuff --------------------------------
  if (fsgd != NULL) {
    if (strlen(fsgd->measname) == 0) {
      strcpy(fsgd->measname,"external");
    }

    const size_t FSGD_datafile_size = 1000;
    if(!yOutFile.empty()) {
      std::string truncName(yOutFile);
      if( yOutFile.size() > (FSGD_datafile_size-1) ) {
	truncName.erase(FSGD_datafile_size);
      }
      snprintf(fsgd->datafile,FSGD_datafile_size,"%s",truncName.c_str());
    }  else {
      std::string truncName(yFile);
      if( yFile.size() > (FSGD_datafile_size-1) ) {
	truncName.erase(FSGD_datafile_size);
      }
      snprintf(fsgd->datafile,FSGD_datafile_size,"%s",truncName.c_str());
    }

    if (surf) {
      strcpy(fsgd->tessellation,"surface");
    } else {
      strcpy(fsgd->tessellation,"volume");
    }

    sprintf(fsgd->DesignMatFile,"X.mat");
    sprintf(tmpstr,"%s/y.fsgd",GLMDir);
    fsgd->ResFWHM = eresfwhm;
    fsgd->LogY = logflag;

    fp = fopen(tmpstr,"w");
    gdfPrintHeader(fp,fsgd);
    fprintf(fp,"Creator          %s\n",Progname);
    fprintf(fp,"SUBJECTS_DIR     %s\n",SUBJECTS_DIR);
    fprintf(fp,"SynthSeed        %d\n",SynthSeed);
    fclose(fp);
  }

  if(DoKurtosis){
    // Compute and save kurtosis
    printf("Computing kurtosis of residuals\n");
    k = fMRIkurtosis(mriglm->eres,mriglm->mask);
    sprintf(tmpstr,"%s/kurtosis.%s",GLMDir,format);
    MRIwrite(k,tmpstr);
    pk = MRIpkurtosis(k, mriglm->glm->dof, mriglm->mask, 10000);
    sprintf(tmpstr,"%s/kurtosis.sig.%s",GLMDir,format);
    MRIwrite(pk,tmpstr);
    MRIfree(&k);
    MRIfree(&pk);
  }
  if(DoSkew){
    // Compute and save skew
    printf("Computing skew of residuals\n");
    k = fMRIskew(mriglm->eres,mriglm->mask);
    sprintf(tmpstr,"%s/skew.%s",GLMDir,format);
    MRIwrite(k,tmpstr);
    pk = MRIpskew(k, mriglm->glm->dof, mriglm->mask, 10000);
    sprintf(tmpstr,"%s/skew.sig.%s",GLMDir,format);
    MRIwrite(pk,tmpstr);
    MRIfree(&pk);
    MRIfree(&k);
  }

  // Compute and save PCA
  if (pcaSave) {
    printf("Computing PCA (%d)\n",npca);
    sprintf(tmpstr,"%s/pca-eres",GLMDir);
    mkdir(tmpstr,0777);
    err=MRIpca(mriglm->eres, &Upca, &Spca, &Vpca, mriglm->mask);
    if (err) exit(1);
    sprintf(tmpstr,"%s/pca-eres/v.%s",GLMDir,format);
    MRIwrite(Vpca,tmpstr);
    sprintf(tmpstr,"%s/pca-eres/u.mtx",GLMDir);
    MatrixWriteTxt(tmpstr, Upca);
    sprintf(tmpstr,"%s/pca-eres/sdiag.mat",GLMDir);
    MatrixWriteTxt(tmpstr, Spca);
    sprintf(tmpstr,"%s/pca-eres/stats.dat",GLMDir);
    WritePCAStats(tmpstr,Spca);
  }

  if(RandSplit){
    sprintf(tmpstr,"%s/synthseed.dat",GLMDir);
    fp = fopen(tmpstr,"w");
    fprintf(fp,"%d",SynthSeed);
    fclose(fp);
  }

  // re-write the log file, adding a few things
  sprintf(tmpstr,"%s/mri_glmfit.log",GLMDir);
  fp = fopen(tmpstr,"w");
  dump_options(fp);
  fprintf(fp,"ResidualFWHM %lf\n",eresfwhm);
  fprintf(fp,"SearchSpace %lf\n",searchspace);
  if (surf != NULL)  fprintf(fp,"anattype surface\n");
  else              fprintf(fp,"anattype volume\n");
  fclose(fp);

  printf("mri_glmfit done\n");
  return(0);
  exit(0);

}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused, msec, niters, frameno,k;
  char **pargv, *option ;
  double rvartmp;
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
    else if (!strcasecmp(option, "--reshape"))   DoReshape = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--save-yhat")) yhatSave = 1;
    else if (!strcasecmp(option, "--yhat-save")) yhatSave = 1;
    else if (!strcasecmp(option, "--save-eres")) eresSave = 1;
    else if (!strcasecmp(option, "--eres-save")) eresSave = 1;
    else if (!strcasecmp(option, "--save-fwhm-map")) {SaveFWHMMap=1;ComputeFWHM = 1;}
    else if (!strcasecmp(option, "--eres-scm")) eresSCMSave = 1;
    else if (!strcasecmp(option, "--save-cond")) condSave = 1;
    else if (!strcasecmp(option, "--dontsave")) DontSave = 1;
    else if (!strcasecmp(option, "--dontsavewn")) DontSaveWn = 1;
    else if (!strcasecmp(option, "--synth"))   synth = 1;
    else if (!strcasecmp(option, "--mask-inv"))  maskinv = 1;
    else if (!strcasecmp(option, "--prune"))    prunemask = 1;
    else if (!strcasecmp(option, "--no-prune")) prunemask = 0;
    else if (!strcasecmp(option, "--w-inv"))  weightinv = 1;
    else if (!strcasecmp(option, "--w-sqrt")) weightsqrt = 1;
    else if (!strcasecmp(option, "--perm-1")) OneSamplePerm = 1;
    else if (!strcasecmp(option, "--osgm"))   {OneSampleGroupMean = 1; DoPCC = 0;}
    else if (!strcasecmp(option, "--diag-cluster")) DiagCluster = 1;
    else if (!strcasecmp(option, "--perm-force")) PermForce = 1;
    else if (!strcasecmp(option, "--perm-nonstatcor")) PermNonStatCor = 1;
    else if (!strcasecmp(option, "--logy")) logflag = 1;
    else if (!strcasecmp(option, "--no-logy")) logflag = 0;
    else if (!strcasecmp(option, "--kurtosis")) DoKurtosis = 1;
    else if (!strcasecmp(option, "--skew")) DoSkew = 1;
    else if (!strcasecmp(option, "--rm-spatial-mean")) RmSpatialMean = 1;
    else if (!strcasecmp(option, "--allow-zero-dof")) AllowZeroDOF = 1;
    else if (!strcasecmp(option, "--prune_thr")){
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&prune_thr); 
      usepruning = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--nii")) format = "nii";
    else if (!strcasecmp(option, "--nii.gz")) format = "nii.gz";
    else if (!strcasecmp(option, "--mgh")) format = "mgh";
    else if (!strcasecmp(option, "--mgz")) format = "mgz";
    else if (!strcasecmp(option, "--allowsubjrep"))
      fsgdf_AllowSubjRep = 1; /* external, see fsgdf.h */
    else if (!strcasecmp(option, "--fsgd-rescale")) fsgdReScale = 1; 
    else if (!strcasecmp(option, "--fisher"))      DoFisher = 1; 
    else if (!strcasecmp(option, "--pcc"))         DoPCC = 1; 
    else if (!strcasecmp(option, "--no-pcc"))      DoPCC = 0; 
    else if (!strcasecmp(option, "--rescale-x"))   ReScaleX = 1; 
    else if (!strcasecmp(option, "--no-rescale-x")) ReScaleX = 0; 
    else if (!strcasecmp(option, "--tar1")) DoTemporalAR1 = 1;
    else if (!strcasecmp(option, "--no-tar1")) DoTemporalAR1 = 0;
    else if (!strcasecmp(option, "--qa")) {
      useqa = 1;
      DoTemporalAR1 = 1;
    }
    else if (!strcasecmp(option, "--no-mask-smooth")) UseMaskWithSmoothing = 0;
    else if (!strcasecmp(option, "--distance")) DoDistance = 1;
    else if (!strcasecmp(option, "--illcond")) IllCondOK = 1;
    else if (!strcasecmp(option, "--no-illcond")) IllCondOK = 0;
    else if (!strcasecmp(option, "--asl")) useasl = 1;
    else if (!strcasecmp(option, "--asl-rev")){
      useasl = 1;
      asl1val = 0;
      asl2val = 1;
    }
    else if (!strcasecmp(option, "--no-contrasts-ok")) NoContrastsOK = 1;
    else if (!strcmp(option, "--no-fix-vertex-area")) {
      printf("Turning off fixing of vertex area\n");
      MRISsetFixVertexAreaValue(0);
    } 
    else if (!strcasecmp(option, "--diag")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag-show")) {
      Gdiag = (Gdiag & DIAG_SHOW);
    } 
    else if (!strcasecmp(option, "--diag-verbose")) {
      Gdiag = (Gdiag & DIAG_VERBOSE);
    } 
    else if (!strcasecmp(option, "--sim")) {
      if (nargc < 4) CMDargNErr(option,4);
      if (CSDcheckSimType(pargv[0])) {
        printf("ERROR: simulation type %s unrecognized, supported values are\n"
               "  perm, mc-full, mc-z\n", pargv[0]);
        exit(1);
      }
      strcpy(csd->simtype,pargv[0]);
      sscanf(pargv[1],"%d",&nsim);
      sscanf(pargv[2],"%lf",&csd->thresh);
      simbase = pargv[3]; // basename
      printf("simbase %s\n",simbase);
      DoSim = 1;
      DontSave = 1;
      prunemask = 0;
      DoPCC = 0;
      nargsused = 4;
    } 
    else if(!strcasecmp(option, "--sim-thresh-loop")) DoSimThreshLoop = 1;
    else if(!strcasecmp(option, "--sim-thresh-loop-pos")){
      DoSimThreshLoop = 1;
      nSignList = 1;
      SignList[0] = +1; // pos
    }
    else if (!strcasecmp(option, "--uniform")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%lf",&UniformMin);
      sscanf(pargv[1],"%lf",&UniformMax);
      UseUniform = 1;
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--sim-sign")) {
      // this applies only to t-tests
      if (nargc < 1) CMDargNErr(option,1);
      if (!strcmp(pargv[0],"abs"))      tSimSign = 0;
      else if (!strcmp(pargv[0],"pos")) tSimSign = +1;
      else if (!strcmp(pargv[0],"neg")) tSimSign = -1;
      else {
        printf("ERROR: --sim-sign argument %s unrecognized\n",pargv[0]);
        exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rand-exclude")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nRandExclude);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--exclude-frame")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&frameno);
      nExclude = 1;
      ExcludeFrames = (int *)calloc(1,sizeof(int));
      ExcludeFrames[0] = frameno;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--exclude-frame-file")) {
      if(nargc < 1) CMDargNErr(option,1);
      ExcludeFrameFile = pargv[0];
      ExcludeFrames = (int *)calloc(1000,sizeof(int));
      fp = fopen(ExcludeFrameFile,"r");
      if(fp == NULL) exit(1);
      nExclude = 0;
      while(1){
	k=fscanf(fp,"%d",&ExcludeFrames[nExclude]);
	if(k==EOF) break;
	nExclude ++;
      }
      fclose(fp);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rand-split")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&NSplits);
      sscanf(pargv[1],"%d",&SplitNo);
      if(SplitNo > NSplits-1){
	printf("ERROR: SplitNo = %d > NSplits-1 = %d\n",SplitNo,NSplits-1);
	exit(1);
      }
      RandSplit = 1;
      nargsused = 2;
    } 
    else if (!strcmp(option, "--really-use-average7")) ReallyUseAverage7 = 1;
    else if (!strcasecmp(option, "--surf") || !strcasecmp(option, "--surface")) {
      if (nargc < 2) CMDargNErr(option,1);
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      if (SUBJECTS_DIR == NULL) {
        printf("ERROR: SUBJECTS_DIR not defined in environment\n");
        exit(1);
      }
      subject = pargv[0];
      if (!strcmp(subject,"average7")) {
        if (!ReallyUseAverage7) {
          printf("\n");
          printf("ERROR: you have selected subject average7. It is recommended that\n");
          printf("you use the fsaverage subject in $FREESURFER_HOME/subjects.\n");
          printf("If you really want to use average7, re-run this program with\n");
          printf("--really-use-average7 as the first argument.\n");
          printf("\n");
          exit(1);
        } else {
          printf("\n");
          printf("INFO: you have selected subject average7 (and REALLY want to use it)\n");
          printf("instead of fsaverage. So I'm going to turn off fixing of vertex area\n");
          printf("to maintain compatibility with the pre-stable3 release.\n");
          printf("\n");
          MRISsetFixVertexAreaValue(0);
        }
      }
      hemi = pargv[1];
      nargsused = 2;
      if(nargc > 2 && !CMDisFlag(pargv[2])){
        surfname = pargv[2];
        nargsused++;
      }
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
      printf("Reading source surface %s\n",tmpstr);
      surf = MRISread(tmpstr) ;
      if (!surf)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read surface %s", Progname, tmpstr) ;
    } else if (!strcasecmp(option, "--seed")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&SynthSeed);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--smooth") ||
             !strcasecmp(option, "--fwhm")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&FWHM);
      csd->nullfwhm = FWHM;
      printf("FWHM = %f\n",FWHM);
      if(std::isnan(FWHM)){
	printf("ERROR: input FWHM is NaN (not a number).\n");
	printf("  Check the mask in the glm directory.\n");
	exit(1);
      }
      if(FWHM < 0){
	printf("ERROR: input FWHM = %f < 0.\n",FWHM);
	exit(1);
      }
      FWHMSet = 1;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--no-fwhm-est")) ComputeFWHM = 0;
    else if (!strcasecmp(option, "--no-est-fwhm")) ComputeFWHM = 0;
    else if (!strcasecmp(option, "--var-smooth") ||
               !strcasecmp(option, "--var-fwhm")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&VarFWHM);
      csd->varfwhm = VarFWHM;
      nargsused = 1;
    } else if (!strcasecmp(option, "--voxdump")) {
      if (nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%d",&voxdump[0]);
      sscanf(pargv[1],"%d",&voxdump[1]);
      sscanf(pargv[2],"%d",&voxdump[2]);
      voxdumpflag = 1;
      nargsused = 3;
    } else if (!strcasecmp(option, "--selfreg")) {
      if (nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%d",&crsSelfReg[nSelfReg][0]);
      sscanf(pargv[1],"%d",&crsSelfReg[nSelfReg][1]);
      sscanf(pargv[2],"%d",&crsSelfReg[nSelfReg][2]);
      nSelfReg++;
      nargsused = 3;
    } else if (!strcasecmp(option, "--profile")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&niters);
      if (SynthSeed < 0) SynthSeed = PDFtodSeed();
      srand48(SynthSeed);
      printf("Starting GLM profile over %d iterations. Seed=%d\n",
             niters,SynthSeed);
      msec = GLMprofile(200, 20, 5, niters);
      printf(" ... msec = %d\n",msec);
      nargsused = 1;
      exit(0);
    } else if (!strcasecmp(option, "--resynthtest")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&niters);
      if (SynthSeed < 0) SynthSeed = PDFtodSeed();
      srand48(SynthSeed);
      printf("Starting GLM resynth test over %d iterations. Seed=%d\n",
             niters,SynthSeed);
      err = GLMresynthTest(niters, &rvartmp);
      if (err) {
        printf("Failed. rvar = %g\n",rvartmp);
        exit(1);
      }
      printf("Passed. rvarmax = %g\n",rvartmp);
      exit(0);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--y")) {
      if (nargc < 1) CMDargNErr(option,1);
      yFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--y-out")) {
      if (nargc < 1) CMDargNErr(option,1);
      yOutFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--table")) {
      if (nargc < 1) CMDargNErr(option,1);
      yFile = fio_fullpath(pargv[0]);
      UseStatTable = 1;
      ComputeFWHM = 0;
      prunemask = 0;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--yffxvar")) {
      if(nargc < 1) CMDargNErr(option,1);
      yffxvarFile = fio_fullpath(pargv[0]);
      DoFFx = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--ffxdof")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&mriglm->ffxdof);
      DoFFx = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--ffxdofdat")) {
      if(nargc < 1) CMDargNErr(option,1);
      fp = fopen(pargv[0],"r");
      if(fp == NULL){
        printf("ERROR: opening %s\n",pargv[0]);
        exit(1);
      }
      fscanf(fp,"%d",&mriglm->ffxdof);
      fclose(fp);
      DoFFx = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--frame-mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      frameMaskFile = pargv[0];
      DoPCC = 0;
      nargsused = 1;
    } else if (!strcmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskFile = pargv[0];
      labelFile = NULL;
      UseCortexLabel = 0;      
      nargsused = 1;
    } 
    else if (!strcmp(option, "--label")) {
      if (nargc < 1) CMDargNErr(option,1);
      labelFile = pargv[0];
      maskFile = NULL;
      UseCortexLabel = 0;      
      nargsused = 1;
    } 
    else if (!strcmp(option, "--cortex"))   UseCortexLabel = 1;
    else if (!strcmp(option, "--no-mask") || !strcmp(option, "--no-cortex")) {
      labelFile = NULL;
      maskFile = NULL;
      UseCortexLabel = 0;
    }
    else if (!strcmp(option, "--w")) {
      if (nargc < 1) CMDargNErr(option,1);
      wFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--wg")) {
      if (nargc < 1) CMDargNErr(option,1);
      wgFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--wls")) {
      if (nargc < 1) CMDargNErr(option,1);
      wFile = pargv[0];
      weightinv = 1;
      weightsqrt = 1;
      nargsused = 1;
      DoPCC = 0; 
    } 
    else if (!strcmp(option, "--X")) {
      if (nargc < 1) CMDargNErr(option,1);
      XFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--dti-X")) {
      if (nargc < 1) CMDargNErr(option,1);
      XFile = pargv[0];
      usedti=3;
      DoPCC = 0;
      logflag = 1;
      format = "nii.gz";
      ComputeFWHM = 0;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--dti")) {
      if(nargc < 1) CMDargNErr(option,1);
      if(CMDnthIsArg(nargc, pargv, 1)){
        bvalfile = pargv[0];
        bvecfile = pargv[1];
	DoPCC = 0;
	usedti=1;
        nargsused = 2;
      }
      else {
	// file with siemens ascii header
        XFile = pargv[0];
	usedti=2;
        nargsused = 1;
      }
      DoPCC = 0;
      logflag = 1;
      format = "nii.gz";
      ComputeFWHM = 0;
    } 
    else if (!strcmp(option, "--mrtm1")) {
      // --mrtm1 cr.dat time.sec.dat 
      // PET Kinetic Modeling, multilinear reference tissue model 1
      // k2 and k2a are per-min
      if(nargc < 2) CMDargNErr(option,1);
      DoMRTM1=1;
      RTM_Cr = MatrixReadTxt(pargv[0], NULL);
      if(RTM_Cr == NULL) exit(1);
      RTM_TimeSec = MatrixReadTxt(pargv[1], NULL);
      if(RTM_TimeSec == NULL) exit(1);
      RTM_TimeMin = MatrixAlloc(RTM_TimeSec->rows,1,MATRIX_REAL);
      for(k=0; k < RTM_TimeSec->rows; k++)
	RTM_TimeMin->rptr[k+1][1] = RTM_TimeSec->rptr[k+1][1]/60;
      RTM_intCr = MatrixCumTrapZ(RTM_Cr, RTM_TimeMin, NULL);
      prunemask = 0;
      NoContrastsOK = 1;
      DoPCC = 0;
      nargsused = 2;
    } 
    else if (!strcmp(option, "--mrtm2")) {
      // --mrtm2 cr.dat time.sec.dat k2pmin 
      // PET Kinetic Modeling, multilinear reference tissue model 2
      // k2 and k2a are per-min
      // R1 = k2/k2p --> not saved
      if(nargc < 3) CMDargNErr(option,1);
      DoMRTM2=1;
      RTM_Cr = MatrixReadTxt(pargv[0], NULL);
      if(RTM_Cr == NULL) exit(1);
      RTM_TimeSec = MatrixReadTxt(pargv[1], NULL);
      if(RTM_TimeSec == NULL) exit(1);
      RTM_TimeMin = MatrixAlloc(RTM_TimeSec->rows,1,MATRIX_REAL);
      for(k=0; k < RTM_TimeSec->rows; k++)
	RTM_TimeMin->rptr[k+1][1] = RTM_TimeSec->rptr[k+1][1]/60;
      sscanf(pargv[2],"%lf",&MRTM2_k2p);
      printf("MRTM2 k2p %g\n",MRTM2_k2p);
      RTM_intCr = MatrixCumTrapZ(RTM_Cr, RTM_TimeMin, NULL);
      MRTM2_x1 = MatrixAlloc(RTM_Cr->rows,1,MATRIX_REAL);
      for(k=0; k < RTM_Cr->rows; k++)
	MRTM2_x1->rptr[k+1][1] = (RTM_Cr->rptr[k+1][1]/MRTM2_k2p + RTM_intCr->rptr[k+1][1]);
      prunemask = 0;
      NoContrastsOK = 1;
      DoPCC = 0;
      nargsused = 3;
    } 
    else if (!strcmp(option, "--logan")) {
      // --logan cr.dat time.sec.dat tstar
      // Logan non-invasive
      if(nargc < 3) CMDargNErr(option,1);
      DoLogan=1;
      RTM_Cr = MatrixReadTxt(pargv[0], NULL);
      if(RTM_Cr == NULL) exit(1);
      RTM_TimeSec = MatrixReadTxt(pargv[1], NULL);
      if(RTM_TimeSec == NULL) exit(1);
      RTM_TimeMin = MatrixAlloc(RTM_TimeSec->rows,1,MATRIX_REAL);
      for(k=0; k < RTM_TimeSec->rows; k++)
	RTM_TimeMin->rptr[k+1][1] = RTM_TimeSec->rptr[k+1][1]/60;
      sscanf(pargv[2],"%lf",&Logan_Tstar);
      printf("Logan tstar %g\n",Logan_Tstar);
      RTM_intCr = MatrixCumTrapZ(RTM_Cr, RTM_TimeMin, NULL);
      prunemask = 0;
      NoContrastsOK = 1;
      DoPCC = 0;
      nargsused = 3;
    } 
    else if (!strcmp(option, "--pvr")) {
      if (nargc < 1) CMDargNErr(option,1);
      pvrFiles[npvr] = pargv[0];
      DoPCC = 0;
      npvr++;
      nargsused = 1;
    } else if (!strcmp(option, "--glmdir") || !strcmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      GLMDir = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--beta")) {
      if (nargc < 1) CMDargNErr(option,1);
      betaFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--rvar")) {
      if (nargc < 1) CMDargNErr(option,1);
      rvarFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--yhat")) {
      if (nargc < 1) CMDargNErr(option,1);
      yhatFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--eres")) {
      if (nargc < 1) CMDargNErr(option,1);
      eresFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--C")) {
      if (nargc < 1) CMDargNErr(option,1);
      CFile[nContrasts] = pargv[0];
      Gamma0File[nContrasts] = NULL;
      nargsused = 1;
      if(CMDnthIsArg(nargc, pargv, 1)) {
	Gamma0File[nContrasts] = pargv[1];
        nargsused++;
      }
      nContrasts++;
    } 
    else if (!strcmp(option, "--pca")) {
      if (CMDnthIsArg(nargc, pargv, 0)) {
        sscanf(pargv[0],"%d",&niters);
        nargsused = 1;
      }
      pcaSave = 1;
    } 
    else if ( !strcmp(option, "--fsgd") ) {
      if (nargc < 1) CMDargNErr(option,1);
      fsgdfile = pargv[0];
      nargsused = 1;
      fsgd = gdfRead(fsgdfile,0);
      if (fsgd==NULL) exit(1);
      if(CMDnthIsArg(nargc, pargv, 1)) {
        gd2mtx_method = pargv[1];
        nargsused ++;
        if (gdfCheckMatrixMethod(gd2mtx_method)) exit(1);
      } 
      else {
	if(strcmp(fsgd->DesignMatMethod,"DOSS") == 0 ||
	   strcmp(fsgd->DesignMatMethod,"doss") == 0)
	  gd2mtx_method = "doss";
	else gd2mtx_method = "dods";
      }
      printf("INFO: gd2mtx_method is %s\n",gd2mtx_method);
      strcpy(fsgd->DesignMatMethod,gd2mtx_method);
    } 
    else if ( !strcmp(option, "--xonly") ) {
      if(nargc < 1) CMDargNErr(option,1);
      XOnlyFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--maxvox")) {
      if (nargc < 1) CMDargNErr(option,1);
      MaxVoxBase = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--subsample")) {
      if (nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&SubSampStart);
      sscanf(pargv[1],"%d",&SubSampDelta);
      SubSample = 1;
      nargsused = 2;
    } 
    else if (!strcmp(option, "--sim-done")) {
      if(nargc < 1) CMDargNErr(option,1);
      SimDoneFile = pargv[0];
      nargsused = 1;
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


/* --------------------------------------------- */
static void print_usage(void) 
{
printf("\n");
printf("USAGE: ./mri_glmfit\n");
printf("\n");
printf("   --glmdir dir : save outputs to dir\n");
printf("\n");
printf("   --y inputfile\n");
printf("   --table stats-table : as output by asegstats2table or aparcstats2table \n");
printf("   --fsgd FSGDF <gd2mtx> : freesurfer descriptor file\n");
printf("   --X design matrix file\n");
printf("   --C contrast1.mtx <--C contrast2.mtx ...>\n");
printf("   --osgm : construct X and C as a one-sample group mean\n");
printf("   --no-contrasts-ok : do not fail if no contrasts specified\n");
printf("   --dti bvals bvecs : do DTI analysis using bvals and bvecs\n");
printf("   --dti siemensdicom : do DTI analysis extracting bvals and bvecs from dicom\n");
printf("   --dti-X X.txt : do DTI analysis using provided matrix\n");
printf("\n");
printf("   --pvr pvr1 <--prv pvr2 ...> : per-voxel regressors\n");
printf("   --selfreg col row slice   : self-regressor from index col row slice\n");
printf("\n");
printf("   --wls yffxvar : weighted least squares\n");
printf("   --yffxvar yffxvar : for fixed effects analysis\n");
printf("   --ffxdof DOF : dof for fixed effects analysis\n");
printf("   --ffxdofdat ffxdof.dat : text file with dof for fixed effects analysis\n");
printf("\n");
printf("   --w weightfile : weight for each input at each voxel\n");
printf("   --w-inv : invert weights\n");
printf("   --w-sqrt : sqrt of (inverted) weights\n");
printf("\n");
printf("   --fwhm fwhm : smooth input by fwhm\n");
printf("   --var-fwhm fwhm : smooth variance by fwhm\n");
printf("   --no-mask-smooth : do not mask when smoothing\n");
printf("   --no-est-fwhm : turn off FWHM output estimation\n");
printf("\n");
printf("   --mask maskfile : binary mask\n");
printf("   --label labelfile : use label as mask, surfaces only\n");
printf("   --no-mask : do NOT use a mask (same as --no-cortex)\n");
printf("   --no-cortex : do NOT use subjects ?h.cortex.label as --label\n");
printf("   --mask-inv : invert mask\n");
printf("   --prune : remove voxels that do not have a non-zero value at each frame (def)\n");
printf("   --no-prune : do not prune\n");
printf("   --logy : compute natural log of y prior to analysis\n");
printf("   --no-logy : do not compute natural log of y prior to analysis\n");
printf("   --rm-spatial-mean : subtract the (masked) mean from each frame\n");
printf("   --yhat-save : save signal estimate (yhat)\n");
printf("   --eres-save : save residual error (eres)\n");
printf("   --eres-scm : save residual error spatial correlation matrix (eres.scm). Big!\n");
printf("   --save-fwhm-map : save voxel-wise map of FWHM estimates\n");
printf("   --y-out y.out.mgh : save input after any pre-processing\n");
printf("\n");
printf("   --surf subject hemi <surfname> : needed for some flags (uses white by default)\n");
printf("\n");
printf("   --skew : compute skew and p-value for skew\n");
printf("   --kurtosis : compute kurtosis and p-value for kurtosis\n");
printf("\n");
printf("   --sim nulltype nsim thresh csdbasename : simulation perm, mc-full, mc-z\n");
printf("   --sim-sign signstring : abs, pos, or neg. Default is abs.\n");
printf("   --uniform min max : use uniform distribution instead of gaussian\n");
printf("\n");
printf("   --pca : perform pca/svd analysis on residual\n");
printf("   --tar1 : compute and save temporal AR1 of residual\n");
printf("   --save-yhat : flag to save signal estimate\n");
printf("   --save-cond  : flag to save design matrix condition at each voxel\n");
printf("   --voxdump col row slice  : dump voxel GLM and exit\n");
printf("\n");
printf("   --seed seed : used for synthesizing noise\n");
printf("   --synth : replace input with gaussian\n");
printf("\n");
printf("   --resynthtest niters : test GLM by resynthsis\n");
printf("   --profile     niters : test speed\n");
printf("\n");
printf("   --mrtm1 RefTac TimeSec : perform MRTM1 kinetic modeling\n");
printf("   --mrtm2 RefTac TimeSec k2prime : perform MRTM2 kinetic modeling\n");
printf("   --logan RefTac TimeSec tstar   : perform Logan kinetic modeling\n");
printf("\n");
printf("   --perm-force : force perumtation test, even when design matrix is not orthog\n");
printf("   --diag Gdiag_no : set diagnositc level\n");
printf("   --diag-cluster : save sig volume and exit from first sim loop\n");
printf("   --debug     turn on debugging\n");
printf("   --checkopts don't run anything, just check options and exit\n");
printf("   --help      print out information on how to use this program\n");
printf("   --version   print out version and exit\n");
printf("   --no-fix-vertex-area : turn off fixing of vertex area (for back comapt only)\n");
printf("   --allowsubjrep allow subject names to repeat in the fsgd file (must appear\n");
printf("                  before --fsgd)\n");
printf("   --allow-zero-dof : mostly for very special purposes\n");
printf("   --illcond : allow ill-conditioned design matrices\n");
printf("   --sim-done SimDoneFile : create DoneFile when simulation finished \n");
printf("\n");
printf("\n");
}


/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
printf("\n");
printf("OUTLINE:\n");
printf("  SUMMARY\n");
printf("  MATHEMATICAL BACKGROUND\n");
printf("  COMMAND-LINE ARGUMENTS\n");
printf("  MONTE CARLO SIMULATION AND CORRECTION FOR MULTIPLE COMPARISONS\n");
printf("\n");
printf("SUMMARY\n");
printf("\n");
printf("Performs general linear model (GLM) analysis in the volume or the\n");
printf("surface.  Options include simulation for correction for multiple\n");
printf("comparisons, weighted LMS, variance smoothing, PCA/SVD analysis of\n");
printf("residuals, per-voxel design matrices, and 'self' regressors. This\n");
printf("program performs both the estimation and inference. This program\n");
printf("is meant to replace mris_glm (which only operated on surfaces).\n");
printf("This program can be run in conjunction with mris_preproc.\n");
printf("\n");
printf("MATHEMATICAL BACKGROUND\n");
printf("\n");
printf("This brief intoduction to GLM theory is provided to help the user\n");
printf("understand what the inputs and outputs are and how to set the\n");
printf("various parameters. These operations are performed at each voxel\n");
printf("or vertex separately (except with --var-fwhm).\n");
printf("\n");
printf("The forward model is given by:\n");
printf("\n");
printf("    y = W*X*B + n\n");
printf("\n");
printf("where X is the Ns-by-Nb design matrix, y is the Ns-by-Nv input data\n");
printf("set, B is the Nb-by-Nv regression parameters, and n is noise. Ns is\n");
printf("the number of inputs, Nb is the number of regressors, and Nv is the\n");
printf("number of voxels/vertices (all cols/rows/slices). y may be surface\n");
printf("or volume data and may or may not have been spatially smoothed. W\n");
printf("is a diagonal weighing matrix.\n");
printf("\n");
printf("During the estimation stage, the forward model is inverted to\n");
printf("solve for B:\n");
printf("\n");
printf("    B = inv(X'W'*W*X)*X'W'y\n");
printf("\n");
printf("The signal estimate is computed as\n");
printf("\n");
printf("    yhat = B*X\n");
printf("\n");
printf("The residual error is computed as\n");
printf("\n");
printf("    eres = y - yhat\n");
printf("\n");
printf("For random effects analysis, the noise variance estimate (rvar) is\n");
printf("computed as the sum of the squares of the residual error divided by\n");
printf("the DOF.  The DOF equals the number of rows of X minus the number of\n");
printf("columns. For fixed effects analysis, the noise variance is estimated\n");
printf("from the lower-level variances passed with --yffxvar, and the DOF\n");
printf("is the sum of the DOFs from the lower level.\n");
printf("\n");
printf("A contrast matrix C has J rows and as many columns as columns of\n");
printf("X. The contrast is then computed as:\n");
printf("\n");
printf("   G = C*B\n");
printf("\n");
printf("The F-ratio for the contrast is then given by:\n");
printf("\n");
printf("   F = G'*inv(C*inv(X'W'*W*X))*C')*G/(J*rvar)\n");
printf("\n");
printf("The F is then used to compute a p-value.  Note that when J=1, this\n");
printf("reduces to a two-tailed t-test.\n");
printf("\n");
printf("\n");
printf("COMMAND-LINE ARGUMENTS\n");
printf("\n");
printf("--glmdir dir\n");
printf("\n");
printf("Directory where output will be saved. Not needed with --sim.\n");
printf("\n");
printf("The outputs will be saved in mgh format as:\n");
printf("  mri_glmfit.log - execution parameters (send with bug reports)\n");
printf("  beta.mgh - all regression coefficients (B above)\n");
printf("  eres.mgh - residual error\n");
printf("  rvar.mgh - residual error variance\n");
printf("  rstd.mgh - residual error stddev (just sqrt of rvar)\n");
printf("  y.fsgd - fsgd file (if one was input)\n");
printf("  wn.mgh - normalized weights (with --w or --wls)\n");
printf("  yhat.mgh - signal estimate (with --save-yhat)\n");
printf("  mask.mgh - final mask (when a mask is used)\n");
printf("  cond.mgh - design matrix condition at each voxel (with --save-cond)\n");
printf("  contrast1name/ - directory for each contrast (see --C)\n");
printf("    C.dat - copy of contrast matrix\n");
printf("    gamma.mgh - contrast (G above)\n");
printf("    F.mgh - F-ratio (t = sign(gamma)*sqrt(F) for t-test contrasts)\n");
printf("    sig.mgh - significance from F-test (actually -log10(p))\n");
printf("    z.mgh - z map computed from the p-value\n");
printf("    pcc.mgh - partial correlation coefficient (for t-tests)\n");
printf("\n");
printf("--y inputfile\n");
printf("\n");
printf("Path to input file with each frame being a separate input. This can be\n");
printf("volume or surface-based, but the file must be one of the 'volume'\n");
printf("formats (eg, mgh, img, nii, etc) accepted by mri_convert. See\n");
printf("mris_preproc for an easy way to generate this file for surface data.\n");
printf("Not with --table.\n");
printf("\n");
printf("--table stats-table\n");
printf("\n");
printf("Use text table as input instead of --y. The stats-table is that of\n");
printf("the form produced by asegstats2table or aparcstats2table.\n");
printf("\n");
printf("--fsgd fname <gd2mtx>\n");
printf("\n");
printf("Specify the global design matrix with a FreeSurfer Group Descriptor\n");
printf("File (FSGDF).  See http://surfer.nmr.mgh.harvard.edu/docs/fsgdf.txt\n");
printf("for more info.  The gd2mtx is the method by which the group\n");
printf("description is converted into a design matrix. Legal values are doss\n");
printf("(Different Offset, Same Slope) and dods (Different Offset, Different\n");
printf("Slope). doss will create a design matrix in which each class has it's\n");
printf("own offset but forces all classes to have the same slope. dods models\n");
printf("each class with it's own offset and slope. In either case, you'll need\n");
printf("to know the order of the regressors in order to correctly specify the\n");
printf("contrast vector. For doss, the first NClass columns refer to the\n");
printf("offset for each class.  The remaining columns are for the continuous\n");
printf("variables. In dods, the first NClass columns again refer to the offset\n");
printf("for each class.  However, there will be NClass*NVar more columns (ie,\n");
printf("one column for each variable for each class). The first NClass columns\n");
printf("are for the first variable, etc. If neither of these models works for\n");
printf("you, you will have to specify the design matrix manually (with --X).\n");
printf("\n");
printf("--no-rescale-x\n");
printf("\n");
printf("By default the inverse of the covariance of the desgin matrix is\n");
printf("computed by rescaling each column of the design matrix prior to the \n");
printf("inverse computation, then rescaling back afterwards. This helps\n");
printf("with designs that are badly scaled. This is completely transparent\n");
printf("to the user. This flag turns this feature off so that the inverse\n");
printf("is computed directly from the design matrix.\n");
printf("\n");
printf("--X design matrix file\n");
printf("\n");
printf("Explicitly specify the design matrix. Can be in simple text or in matlab4\n");
printf("format. If matlab4, you can save a matrix with save('X.mat','X','-v4');\n");
printf("\n");
printf("--C contrast1.mtx <--C contrast2.mtx ...>\n");
printf("\n");
printf("Specify one or more contrasts to test. The contrast.mtx file is an\n");
printf("ASCII text file with the contrast matrix in it (make sure the last\n");
printf("line is blank). The name can be (almost) anything. If the extension is\n");
printf(".mtx, .mat, .dat, or .con, the extension will be stripped of to form\n");
printf("the directory output name.  The output will be saved in\n");
printf("glmdir/contrast1, glmdir/contrast2, etc. Eg, if --C norm-v-cont.mtx,\n");
printf("then the ouput will be glmdir/norm-v-cont.\n");
printf("\n");
printf("--osgm\n");
printf("\n");
printf("Construct X and C as a one-sample group mean. X is then a one-column\n");
printf("matrix filled with all 1s, and C is a 1-by-1 matrix with value 1.\n");
printf("You cannot specify both --X and --osgm. A contrast cannot be specified\n");
printf("either. The contrast name will be osgm.\n");
printf("\n");
printf("--pvr pvr1 <--prv pvr2 ...>\n");
printf("\n");
printf("Per-voxel (or vertex) regressors (PVR). Normally, the design matrix is\n");
printf("'global', ie, the same matrix is used at each voxel. This option allows the\n");
printf("user to specify voxel-specific regressors to append to the design\n");
printf("matrix. Note: the contrast matrices must include columns for these\n");
printf("components.\n");
printf("\n");
printf("--selfreg col row slice\n");
printf("\n");
printf("Create a 'self-regressor' from the input data based on the waveform at\n");
printf("index col row slice. This waveform is residualized and then added as a\n");
printf("column to the design matrix. Note: the contrast matrices must include\n");
printf("columns for this component.\n");
printf("\n");
printf("--wls yffxvar : weighted least squares\n");
printf("\n");
printf("Perform weighted least squares (WLS) random effects analysis instead\n");
printf("of ordinary least squares (OLS).  This requires that the lower-level\n");
printf("variances be available.  This is often the case with fMRI analysis but\n");
printf("not with an anatomical analysis. Note: this should not be confused\n");
printf("with fixed effects analysis. The weights will be inverted,\n");
printf("square-rooted, and normalized to sum to the number of inputs for each\n");
printf("voxel. Same as --w yffxvar --w-inv --w-sqrt (see --w below).\n");
printf("\n");
printf("--yffxvar yffxvar      : for fixed effects analysis\n");
printf("--ffxdof DOF           : DOF for fixed effects analysis\n");
printf("--ffxdofdat ffxdof.dat : text file with DOF for fixed effects analysis\n");
printf("\n");
printf("Perform fixed-effect analysis. This requires that the lower-level variances\n");
printf("be available.  This is often the case with fMRI analysis but not with\n");
printf("an anatomical analysis. Note: this should not be confused with weighted\n");
printf("random effects analysis (wls). The dof is the sum of the DOFs from the\n");
printf("lower levels.\n");
printf("\n");
printf("--w weightfile\n");
printf("--w-inv\n");
printf("--w-sqrt\n");
printf("\n");
printf("Perform weighted LMS using per-voxel weights from the weightfile. The\n");
printf("data in weightfile must have the same dimensions as the input y\n");
printf("file. If --w-inv is flagged, then the inverse of each weight is used\n");
printf("as the weight.  If --w-sqrt is flagged, then the square root of each\n");
printf("weight is used as the weight.  If both are flagged, the inverse is\n");
printf("done first. The final weights are normalized so that the sum at each\n");
printf("voxel equals the number of inputs. The normalized weights are then\n");
printf("saved in glmdir/wn.mgh.  The --w-inv and --w-sqrt flags are useful\n");
printf("when passing contrast variances from a lower level analysis to a\n");
printf("higher level analysis (as is often done in fMRI).\n");
printf("\n");
printf("--fwhm fwhm\n");
printf("\n");
printf("Smooth input with a Gaussian kernel with the given full-width/half-maximum\n");
printf("(fwhm) specified in mm. If the data are surface-based, then you must\n");
printf("specify --surf, otherwise mri_glmfit assumes that the input is a volume\n");
printf("and will perform volume smoothing.\n");
printf("\n");
printf("--var-fwhm fwhm\n");
printf("\n");
printf("Smooth residual variance map with a Gaussian kernel with the given\n");
printf("full-width/half-maximum (fwhm) specified in mm. If the data are\n");
printf("surface-based, then you must specify --surf, otherwise mri_glmfit\n");
printf("assumes that the input is a volume and will perform volume smoothing.\n");
printf("\n");
printf("--mask maskfile\n");
printf("--label labelfile\n");
printf("--mask-inv\n");
printf("--cortex \n");
printf("\n");
printf("Only perform analysis where mask=1. All other voxels will be set to 0.\n");
printf("If using surface, then labelfile will be converted to a binary mask\n");
printf("(requires --surf). By default, the label file for surfaces is \n");
printf("?h.cortex.label. To force a no-mask with surfaces, use --no-mask or \n");
printf("--no-cortex. If --mask-inv is flagged, then performs analysis\n");
printf("only where mask=0. If performing a simulation (--sim), map maximums\n");
printf("and clusters will only be searched for in the mask. The final binary\n");
printf("mask will automatically be saved in glmdir/mask.mgh\n");
printf("\n");
printf("--prune\n");
printf("--no-prune\n");
printf("\n");
printf("This happens by default. Use --no-prune to turn it off. Remove voxels\n");
printf("from the analysis if the ALL the frames at that voxel\n");
printf("do not have an absolute value that exceeds zero (actually FLT_MIN, \n");
printf("or whatever is set by --prune_thr). This helps to prevent the situation \n");
printf("where some frames are 0 and others are not. If no mask is supplied, \n");
printf("a mask is created and saved. If a mask is supplied, it is pruned, and \n");
printf("the final mask is saved. Do not use with --sim. Rather, run the non-sim \n");
printf("analysis with --prune, then pass the created mask when running simulation. \n");
printf("It is generally a good idea to prune. --no-prune will turn off pruning \n");
printf("if it had been turned on. For DTI, only the first frame is used to \n");
printf("create the mask.\n");
printf("\n");
printf("--prune_thr threshold\n");
printf("\n");
printf("Use threshold to create the mask using pruning. Default is FLT_MIN\n");
printf("\n");
printf("--surf subject hemi <surfname>\n");
printf("\n");
printf("Specify that the input has a surface geometry from the hemisphere of the\n");
printf("given FreeSurfer subject. This is necessary for smoothing surface data\n");
printf("(--fwhm or --var-fwhm), specifying a label as a mask (--label), or\n");
printf("running a simulation (--sim) on surface data. If --surf is not specified,\n");
printf("then mri_glmfit will assume that the data are volume-based and use\n");
printf("the geometry as specified in the header to make spatial calculations.\n");
printf("By default, the white surface is used, but this can be overridden by\n");
printf("specifying surfname.\n");
printf("\n");
printf("--pca\n");
printf("\n");
printf("Flag to perform PCA/SVD analysis on the residual. The result is stored\n");
printf("in glmdir/pca-eres as v.mgh (spatial eigenvectors), u.mtx (frame\n");
printf("eigenvectors), sdiag.mat (singular values). eres = u*s*v'. The matfiles\n");
printf("are just ASCII text. The spatial EVs can be loaded as overlays in\n");
printf("tkmedit or tksurfer. In addition, there is stats.dat with 5 columns:\n");
printf("  (1) component number\n");
printf("  (2) variance spanned by that component\n");
printf("  (3) cumulative variance spanned up to that component\n");
printf("  (4) percent variance spanned by that component\n");
printf("  (5) cumulative percent variance spanned up to that component\n");
printf("\n");
printf("--save-yhat\n");
printf("\n");
printf("Flag to save the signal estimate (yhat) as glmdir/yhat.mgh. Normally, this\n");
printf("pis not very useful except for debugging.\n");
printf("\n");
printf("--save-cond\n");
printf("\n");
printf("Flag to save the condition number of the design matrix at eaach voxel.\n");
printf("Normally, this is not very useful except for debugging. It is totally\n");
printf("useless if not using weights or PVRs.\n");
printf("\n");
printf("--nii, --nii.gz\n");
printf("\n");
printf("Use nifti (or compressed nifti) as output format instead of mgh. This\n");
printf("will work with surfaces, but you will not be able to open the output\n");
printf("nifti files with non-freesurfer software.\n");
printf("\n");
printf("--seed seed\n");
printf("\n");
printf("Use seed as the seed for the random number generator. By default, mri_glmfit\n");
printf("will select a seed based on time-of-day. This only has an effect with\n");
printf("--sim or --synth.\n");
printf("\n");
printf("--synth\n");
printf("\n");
printf("Replace input data with whise gaussian noise. This is good for testing.\n");
printf("\n");
printf("--voxdump col row slice\n");
printf("\n");
printf("Save GLM data for a single voxel in directory glmdir/voxdump-col-row-slice.\n");
printf("Exits immediately. Good for debugging.\n");
printf("\n");
printf("\n");
printf("MONTE CARLO SIMULATION AND CORRECTION FOR MULTIPLE COMPARISONS\n");
printf("\n");
printf("One method for correcting for multiple comparisons is to perform simulations\n");
printf("under the null hypothesis and see how often the value of a statistic\n");
printf("from the 'true' analysis is exceeded. This frequency is then interpreted\n");
printf("as a p-value which has been corrected for multiple comparisons. This\n");
printf("is especially useful with surface-based data as traditional random\n");
printf("field theory is harder to implement. This simulator is roughly based\n");
printf("on FSLs permuation simulator (randomise) and AFNIs null-z simulator\n");
printf("(AlphaSim). Note that FreeSurfer also offers False Discovery Rate (FDR)\n");
printf("correction in tkmedit and tksurfer.\n");
printf("\n");
printf("The estimation, simulation, and correction are done in three distinct\n");
printf("phases:\n");
printf("  1. Estimation: run the analysis on your data without simulation.\n");
printf("     At this point you can view your results (see if FDR is\n");
printf("     sufficient:).\n");
printf("  2. Simulation: run the simulator with the same parameters\n");
printf("     as the estimation to get the Cluster Simulation Data (CSD).\n");
printf("  3. Clustering: run mri_surfcluster or mri_volcluster with the CSD\n");
printf("     from the simulator and the output of the estimation. These\n");
printf("     programs will print out clusters along with their p-values.\n");
printf("\n");
printf("The Estimation step is described in detail above. The simulation\n");
printf("is invoked by calling mri_glmfit with the following arguments:\n");
printf("\n");
printf("--sim nulltype nsim thresh csdbasename\n");
printf("--sim-sign sign\n");
printf("\n");
printf("It is not necessary to specify --glmdir (it will be ignored). If\n");
printf("you are analyzing surface data, then include --surf.\n");
printf("\n");
printf("nulltype is the method of generating the null data. Legal values are:\n");
printf("  (1) perm - perumation, randomly permute rows of X (cf FSL randomise)\n");
printf("  (2) mc-full - replace input with white gaussian noise\n");
printf("  (3) mc-z - do not actually do analysis, just assume the output\n");
printf("      is z-distributed (cf ANFI AlphaSim)\n");
printf("nsim - number of simulation iterations to run (see below)\n");
printf("thresh - threshold, specified as -log10(pvalue) to use for clustering\n");
printf("csdbasename - base name of the file to store the CSD data in. Each\n");
printf("  contrast will get its own file (created by appending the contrast\n");
printf("  name to the base name). A '.csd' is appended to each file name.\n");
printf("\n");
printf("Multiple simulations can be run in parallel by specifying different\n");
printf("csdbasenames. Then pass the multiple CSD files to mri_surfcluster\n");
printf("and mri_volcluster. The Full CSD file is written on each iteration,\n");
printf("which means that the CSD file will be valid if the simulation\n");
printf("is aborted or crashes.\n");
printf("\n");
printf("In the cases where the design matrix is a single columns of ones\n");
printf("(ie, one-sample group mean), it makes no sense to permute the\n");
printf("rows of the design matrix. mri_glmfit automatically checks\n");
printf("for this case. If found, the design matrix is rebuilt on each\n");
printf("permutation with randomly selected +1 and -1s. Same as the -1\n");
printf("option to FSLs randomise.\n");
printf("\n");
printf("--sim-sign sign\n");
printf("\n");
printf("sign is either abs (default), pos, or neg. pos/neg tell mri_glmfit to\n");
printf("perform a one-tailed test. In this case, the contrast matrix can\n");
printf("only have one row.\n");
printf("\n");
printf("--uniform min max\n");
printf("\n");
printf("For mc-full, synthesize input as a uniform distribution between min\n");
printf("and max. \n");
printf("\n");
  exit(1) ;
}

/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if(XFile == NULL && bvalfile == NULL && fsgdfile == NULL &&
     ! OneSampleGroupMean && ! useasl && !useqa && !DoMRTM1 && !DoMRTM2 && !DoLogan) {
    printf("ERROR: must specify an input X file or fsgd file or --osgm\n");
    exit(1);
  }
  if (XFile && fsgdfile ) {
    printf("ERROR: cannot specify both X file and fsgd file\n");
    exit(1);
  }
  if (XFile && OneSampleGroupMean) {
    printf("ERROR: cannot specify both X file and --osgm\n");
    exit(1);
  }
  if(XOnlyFile != NULL){
    if(fsgdfile == NULL) {
      printf("ERROR: you must spec --fsgd with --xonly\n");
      exit(1);
    }
    printf("INFO: gd2mtx_method is %s\n",gd2mtx_method);
    Xtmp = gdfMatrix(fsgd,gd2mtx_method,NULL);
    if(Xtmp==NULL) exit(1);
    Xnorm = MatrixNormalizeCol(Xtmp,NULL,NULL);
    Xcond = MatrixNSConditionNumber(Xnorm);
    printf("Matrix condition is %g\n",Xcond);
    MatrixWriteTxt(XOnlyFile, Xtmp);
    exit(0);
  }

  if(yFile.size() == 0) {
    printf("ERROR: must specify input y file\n");
    exit(1);
  }
  if (nContrasts > 0 && OneSampleGroupMean) {
    printf("ERROR: cannot specify --C with --osgm\n");
    exit(1);
  }

  if(OneSampleGroupMean || usedti || useasl || useqa) NoContrastsOK = 1;
  if(fsgdfile && fsgd->nContrasts != 0 && nContrasts != 0){
    printf("ERROR: cannot have contrasts in FSGD and on the command-line\n");
    exit(1);
  }
  if(fsgdfile && fsgd->nContrasts != 0) nContrasts = fsgd->nContrasts;

  if(nContrasts == 0 && ! NoContrastsOK) {
    printf("ERROR: no contrasts specified.\n");
    exit(1);
  }
  if (GLMDir == NULL && !DontSave) {
    printf("ERROR: must specify GLM output dir\n");
    exit(1);
  }
  if (GLMDir != NULL) {
    sprintf(tmpstr,"%s/beta.%s",GLMDir,format);
    betaFile = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/rvar.%s",GLMDir,format);
    rvarFile = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/eres.%s",GLMDir,format);
    if(eresSave) eresFile = strcpyalloc(tmpstr);
    if(eresSCMSave){
      sprintf(tmpstr,"%s/eres.scm.%s",GLMDir,format);
      eresSCMFile = strcpyalloc(tmpstr);
    }
    if (yhatSave) {
      sprintf(tmpstr,"%s/yhat.%s",GLMDir,format);
      yhatFile = strcpyalloc(tmpstr);
    }
    if (condSave) {
      sprintf(tmpstr,"%s/cond.%s",GLMDir,format);
      condFile = strcpyalloc(tmpstr);
    }
  }
  if (SUBJECTS_DIR == NULL) {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL) {
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }
  if(UseCortexLabel && maskFile){
    printf("ERROR: cannot specify both --cortex and --mask\n");
    exit(1);
  }
  if(UseCortexLabel && surf){
    sprintf(tmpstr,"%s/%s/label/%s.cortex.label",SUBJECTS_DIR,subject,hemi);
    labelFile = strcpyalloc(tmpstr);
  }
  if(labelFile != NULL && surf==NULL) {
    printf("ERROR: need --surf with --label\n");
    exit(1);
  }
  if(prunemask && DoSim) {
    printf("ERROR: do not use --prune with --sim\n");
    exit(1);
  }
  if(DoSim && VarFWHM > 0 &&
      (!strcmp(csd->simtype,"mc-z") || !strcmp(csd->simtype,"mc-t"))) {
    printf("ERROR: cannot use variance smoothing with mc-z or "
           "mc-t simulation\n");
    exit(1);
  }
  if(DoFFx){
    if(mriglm->ffxdof == 0){
      printf("ERROR: you need to specify the dof with FFx\n");
      exit(1);
    }
    if(yffxvarFile.size() == 0){
      printf("ERROR: you need to specify the yffxvar file with FFx\n");
      exit(1);
    }
  }
  if(UseUniform){
    if(DoSim == 0 && synth == 0){
      printf("ERROR: need --sim or --synth with --uniform\n");
      exit(1);
    }
    if(DoSim && strcmp(csd->simtype,"mc-full") != 0){
      printf("ERROR: must use mc-full with --uniform\n");
      exit(1);
    }
    if(UniformMax <= UniformMin){
      printf("ERROR: uniform max (%lf) <= min (%lf)\n",UniformMax,UniformMin);
      exit(1);
    }
  }
  if(DoSim && FWHMSet == 0){
    printf("ERROR: you must supply --fwhm with --sim, even if it is 0\n");
    exit(1);
  }
  // Not sure why this was here, but it does not seem to apply anymore
  //if(DoSimThreshLoop && strcmp(csd->simtype,"mc-z")){
  //printf("ERROR: you can only use --sim-thresh-loop with mc-z\n");
    //exit(1);
  //}

  if(fsgdReScale){
    if(fsgdfile == NULL){
      printf("ERROR: need an fsgd file to use --fsgd-rescale\n");
      exit(1);
    }
    fsgd->ReScale = 1;
  }
  if(wFile && wgFile){
    printf("ERROR: cannot have --w and --wg\n");
    exit(1);
  }
  if(DoPCC && npvr != 0){
    printf("ERROR: cannot have compute pcc with pvr\n");
    exit(1);
  }


  return;
}


/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  int n;

  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"FixVertexAreaFlag = %d\n",MRISgetFixVertexAreaValue());

  fprintf(fp,"UseMaskWithSmoothing     %d\n",UseMaskWithSmoothing);
  if (FWHM > 0) {
    fprintf(fp,"fwhm     %lf\n",FWHM);
    if (surf != NULL) fprintf(fp,"niters    %lf\n",SmoothLevel);
  }
  if (VarFWHM > 0) {
    fprintf(fp,"varfwhm  %lf\n",VarFWHM);
    if (surf != NULL) fprintf(fp,"varniters %lf\n",VarSmoothLevel);
  }

  if (synth) fprintf(fp,"SynthSeed = %d\n",SynthSeed);
  fprintf(fp,"OneSampleGroupMean %d\n",OneSampleGroupMean);

  fprintf(fp,"y    %s\n",yFile.c_str());
  fprintf(fp,"logyflag %d\n",logflag);
  if (XFile)     fprintf(fp,"X    %s\n",XFile);
  fprintf(fp,"usedti  %d\n",usedti);
  if (fsgdfile)  fprintf(fp,"FSGD %s\n",fsgdfile);
  if (labelFile) fprintf(fp,"labelmask  %s\n",labelFile);
  if (maskFile)  fprintf(fp,"mask %s\n",maskFile);
  if (labelFile || maskFile) fprintf(fp,"maskinv %d\n",maskinv);

  fprintf(fp,"glmdir %s\n",GLMDir);
  fprintf(fp,"IllCondOK %d\n",IllCondOK);
  fprintf(fp,"ReScaleX %d\n",ReScaleX);

  for (n=0; n < nSelfReg; n++) {
    fprintf(fp,"SelfRegressor %d  %4d %4d %4d\n",n+1,
            crsSelfReg[n][0],crsSelfReg[n][1],crsSelfReg[n][2]);
  }
  if(SubSample){
    fprintf(fp,"SubSampStart %d\n",SubSampStart);
    fprintf(fp,"SubSampDelta %d\n",SubSampDelta);
  }
  if(DoDistance)  fprintf(fp,"DoDistance %d\n",DoDistance);
  fprintf(fp,"DoFFx %d\n",DoFFx);
  if(DoFFx){
    fprintf(fp,"FFxDOF %d\n",mriglm->ffxdof);
    fprintf(fp,"yFFxVar %s\n",yffxvarFile.c_str());
  }
  if(wgFile) fprintf(fp,"wgFile %s\n",wgFile);
  if(wFile){
    fprintf(fp,"wFile %s\n",wFile);
    fprintf(fp,"weightinv  %d\n",weightinv);
    fprintf(fp,"weightsqrt %d\n",weightsqrt);
  }
  if(UseUniform)
    fprintf(fp,"Uniform %lf %lf\n",UniformMin,UniformMax);

  return;
}


/*--------------------------------------------------------------------*/
static int SmoothSurfOrVol(MRIS *surf, MRI *mri, MRI *mask, double SmthLevel) {
  extern int DoSim;
  extern Timer mytimer;
  extern int UseMaskWithSmoothing;
  double gstd;

  if (surf == NULL) {
    gstd = SmthLevel/sqrt(log(256.0));
    if (!DoSim || debug || Gdiag_no > 0)
      printf("  Volume Smoothing by FWHM=%lf, Gstd=%lf, t=%lf\n",
             SmthLevel,gstd,mytimer.seconds());
    if(UseMaskWithSmoothing)
      MRImaskedGaussianSmooth(mri, mask, gstd, mri);
    else
      MRImaskedGaussianSmooth(mri, NULL, gstd, mri);
    if (!DoSim || debug || Gdiag_no > 0)
      printf("  Done Volume Smoothing t=%lf\n",mytimer.seconds());
  }
  else {
    if (!DoSim || debug || Gdiag_no > 0)
      printf("  Surface Smoothing by %d iterations, t=%lf\n",
             (int)SmthLevel,mytimer.seconds());
    if(UseMaskWithSmoothing)
      MRISsmoothMRI(surf, mri, SmthLevel, mask, mri);
    else
      MRISsmoothMRI(surf, mri, SmthLevel, NULL, mri);
    if (!DoSim || debug || Gdiag_no > 0)
      printf("  Done Surface Smoothing t=%lf\n",mytimer.seconds());
  }
  return(0);
}


/*--------------------------------------------------------------------*/
int MRISmaskByLabel(MRI *y, MRIS *surf, LABEL *lb, int invflag) {
  int **crslut, *lbmask, vtxno, n, c, r, s, f;

  lbmask = (int*) calloc(surf->nvertices,sizeof(int));

  // Set each label vertex in lbmask to 1
  for (n=0; n<lb->n_points; n++) {
    vtxno = lb->lv[n].vno;
    lbmask[vtxno] = 1;
  }

  crslut = MRIScrsLUT(surf, y);
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    if (lbmask[vtxno] && !invflag) continue; // in-label and not inverting
    if (!lbmask[vtxno] && invflag) continue; // out-label but inverting
    c = crslut[0][vtxno];
    r = crslut[1][vtxno];
    s = crslut[2][vtxno];
    for (f=0; f < y->nframes; f++) MRIsetVoxVal(y,c,r,s,f,0);
  }

  free(lbmask);
  MRIScrsLUTFree(crslut);
  return(0);
}

/*--------------------------------------------------------------*/
MRI *BindingPotential(MRI *k2, MRI *k2a, MRI *mask, MRI *bp)
{
  int c, r, s;
  double k2v, k2av, bpv;

  if (bp==NULL){
    bp = MRIallocSequence(k2->width,k2->height,k2->depth,MRI_FLOAT,1);
    if(bp==NULL){
      printf("ERROR: BindingPotential(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(k2,bp);
  }

  for(c=0; c < k2->width; c++)  {
    for(r=0; r < k2->height; r++)    {
      for(s=0; s < k2->depth; s++)   {
	if(mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5){
	  MRIsetVoxVal(bp,c,r,s,0,0.0);
	  continue;
	}
	k2v  = MRIgetVoxVal(k2,c,r,s,0);
	k2av = MRIgetVoxVal(k2a,c,r,s,0);
	bpv = k2v/(k2av+DBL_EPSILON) - 1.0;
	MRIsetVoxVal(bp,c,r,s,0, bpv);
      }
    }
  }

  return(bp);
}


/*--------------------------------------------------------------*/
MRI *MRIconjunct3(MRI *sig1, MRI *sig2, MRI *sig3, MRI *mask, MRI *c123)
{
  int c, r, s;
  MRI *f3;
  double sigv;

  f3 = MRIallocSequence(sig1->width,sig1->height,sig1->depth,MRI_FLOAT,3);

  for(c=0; c < sig1->width; c++)  {
    for(r=0; r < sig1->height; r++)    {
      for(s=0; s < sig1->depth; s++)   {
	if(mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
	sigv = MRIgetVoxVal(sig1,c,r,s,0);
	MRIsetVoxVal(f3,c,r,s,0, sigv);
	sigv = MRIgetVoxVal(sig2,c,r,s,0);
	MRIsetVoxVal(f3,c,r,s,1, sigv);
	sigv = MRIgetVoxVal(sig3,c,r,s,0);
	MRIsetVoxVal(f3,c,r,s,2, sigv);
      }
    }
  }

  c123 = MRIconjunct(f3, c123);
  MRIfree(&f3);
  return(c123);
}

/*--------------------------------------------------------------*/
double GLMEfficiency(MATRIX *X, MATRIX *C)
{
  double efficiency;
  MATRIX *Xt, *Ct, *XtX, *iXtX, *A, *M;

  Xt = MatrixTranspose(X,NULL);
  Ct = MatrixTranspose(C,NULL);

  XtX = MatrixMultiplyD(Xt,X,NULL);
  iXtX = MatrixInverse(XtX,NULL);
  // M = C*inv(X'*X)*C'
  A = MatrixMultiplyD(C,iXtX,NULL);
  M = MatrixMultiplyD(A,Ct,NULL);

  efficiency = 1/MatrixTrace(M);

  MatrixFree(&Xt);
  MatrixFree(&Ct);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&A);
  MatrixFree(&M);

  return(efficiency);
}


/*!
  \fn int GLMdiagnoseDesignMatrix(MATRIX *X)
  \brief Simple routines to determine whether a design matrix is ill-conditioned
  because of some common mistakes.
 */
int GLMdiagnoseDesignMatrix(MATRIX *X)
{
  int ret,r,c,c2,all0;
  float *Xsumsq;
  float XsumsqMin,XsumsqMax;
  int cXsumsqMin, cXsumsqMax;

  ret = 0;
  
  // Check the scale
  Xsumsq = (float *)calloc(X->cols,sizeof(float));
  XsumsqMin = 10e10;
  XsumsqMax = 0;
  cXsumsqMin = 0;
  cXsumsqMax = 0;
  for(c=1; c <= X->cols; c++){
    for(r=1; r <= X->rows; r++) 
      Xsumsq[c-1] += (X->rptr[r][c]*X->rptr[r][c]);
    Xsumsq[c-1] = sqrt(Xsumsq[c-1]);
    if(XsumsqMax < Xsumsq[c-1]){
      XsumsqMax = Xsumsq[c-1];
      cXsumsqMax = c;
    }
    if(XsumsqMin > Xsumsq[c-1]){
      XsumsqMin = Xsumsq[c-1];
      cXsumsqMin = c;
    }
  }
  printf("SumSq: Min=%f (col %d), Max=%f (col %d)\n",XsumsqMin,cXsumsqMin,XsumsqMax,cXsumsqMax);
  if(XsumsqMax/XsumsqMin > 10){
    printf(" The scale is much different between columns %d and %d, you may want to \n",cXsumsqMin,cXsumsqMax);
    printf(" normalize by subtracting the mean and dividing by the standard deviation.\n");
    ret=1;
  }
  free(Xsumsq);

  // Check whether any column is all 0s
  for(c=1; c <= X->cols; c++){
    all0 = 1;
    for(r=1; r <= X->rows; r++) {
      if(fabs(X->rptr[r][c])>FLT_EPSILON) {
	all0 = 0;
	break;
      }
    }
    if(all0) {
      printf("Column %d,  all values are 0\n",c);
      ret = 1;
    }
  }

  // Check whether any column equals any other column
  for(c=1; c <= X->cols-1; c++){
    for(c2=c+1; c2 <= X->cols; c2++){
      all0 = 1;
      for(r=1; r <= X->rows; r++) {
	if(fabs(X->rptr[r][c]-X->rptr[r][c2])>FLT_EPSILON) {
	  all0 = 0;
	  break;
	}
      }
      if(all0) {
	printf("Columns %d and %d are the same\n",c,c2);
	ret = 1;
      }
    }
  }

  return(ret);
}



/*!
  \fn MRI *fMRIskew(MRI *y, MRI *mask)
  \brief Computes a voxel-wise map of skew.
  See also MRIkurtosis()
*/
MRI *fMRIskew(MRI *y, MRI *mask)
{
  MRI *k;
  int c, r, s, f;
  double v,mn,m3=0,m2=0,delta,skew;
  k = MRIallocSequence(y->width,y->height,y->depth,MRI_FLOAT,1);
  MRIcopyHeader(y,k);

  for(c=0; c < y->width; c++){
    for(r=0; r < y->height; r++){
      for(s=0; s < y->depth; s++){
	if(mask){
	  v = MRIgetVoxVal(mask,c,r,s,0);
	  if(v < 0.0001) {
	    MRIsetVoxVal(k,c,r,s,0,0.0);
	    continue;
	  }
	}
	mn = 0;
	for(f=0; f < y->nframes; f++) mn += MRIgetVoxVal(y,c,r,s,f);
	mn /= y->nframes;
	m2 = 0;
	m3 = 0;
	for(f=0; f < y->nframes; f++){
	  delta = MRIgetVoxVal(y,c,r,s,f)-mn;
	  m2 += pow(delta,2.0); // sum of squares
	  m3 += pow(delta,3.0); // sum of cubics
	}
	if(m2 < FLT_EPSILON) continue;
	m2 /= y->nframes;
	m3 /= (y->nframes-1);
	skew = m3/pow(m2,1.5);
	MRIsetVoxVal(k,c,r,s,0,skew);
      }
    }
  }
  return(k);
}

/*!
  \fn MRI *MRIpskew(MRI *kvals, int dof, MRI *mask, int nsamples)
  \brief Computes the p-value for a skew map. Uses simulation.
  \param kvals - source volume with skew values.
  \param dof - dof that the skew was computed from
  \param mask - skip voxels where mask < 0.0001
  \param nsamples - samples to use in the simulation (eg, 10000)
  See also MRIpkurtosis()
*/
MRI *MRIpskew(MRI *kvals, int dof, MRI *mask, int nsamples)
{
  MRI *nmri, *kmri, *pkmri;
  double *ksynth,pk,kvox,v;
  int m,c,r,s,f,ind;

  nmri = MRIrandn(nsamples,1,1,dof,0,1,NULL);
  kmri = fMRIskew(nmri,NULL);

  ksynth = (double*) calloc(sizeof(double),nsamples);
  for(m=0; m < nsamples; m++) ksynth[m] = MRIgetVoxVal(kmri,m,0,0,0);

  qsort(ksynth,nsamples,sizeof(double),CompareDoubles);

  pkmri = MRIclone(kvals,NULL);
  for(c=0; c < pkmri->width; c++){
    for(r=0; r < pkmri->height; r++){
      for(s=0; s < pkmri->depth; s++){
	if(mask){
	  v = MRIgetVoxVal(mask,c,r,s,0);
	  if(v < 0.0001) {
	    for(f=0; f < pkmri->nframes; f++) MRIsetVoxVal(pkmri,c,r,s,f,0.0);
	    continue;
	  }
	}
	for(f=0; f < pkmri->nframes; f++){
	  kvox = MRIgetVoxVal(kvals,c,r,s,f);
	  ind = PDFsearchOrderedTable(kvox,ksynth,nsamples);
	  pk = 1.0 - (double)ind/nsamples;
	  MRIsetVoxVal(pkmri,c,r,s,f,-log10(pk));
	}
      }
    }
  }
  MRIfree(&nmri);
  MRIfree(&kmri);
  free(ksynth);
  return(pkmri);
}

MRI *MRIremoveSpatialMean(MRI *vol, MRI *mask, MRI *out)
{
  int c, r, s,f;
  long nhits;
  double v, vmean;

  if(out == NULL) {
    out = MRIallocSequence(vol->width,vol->height,vol->depth,MRI_FLOAT,vol->nframes);
    MRIcopyHeader(vol,out);
    MRIcopyPulseParameters(vol,out);
  }

  for (f=0; f < vol->nframes; f++){
    // Sum the values for this frame
    nhits = 0;
    v = 0;
    for(c=0; c < vol->width; c++) {
      for(r=0; r < vol->height; r++) {
	for(s=0; s < vol->depth; s++) {
	  if(mask && MRIgetVoxVal(mask,c,r,s,0)<0.5) continue;
	  nhits++;
          v += MRIgetVoxVal(vol,c,r,s,f);
	}//s
      }//r
    }//s
    // Compute the mean for this frame
    vmean = v/nhits;
    if(Gdiag_no > 0) printf("%2d %g\n",f,vmean);
    // Subtract the spatial mean of this frame
    for(c=0; c < vol->width; c++) {
      for(r=0; r < vol->height; r++) {
	for(s=0; s < vol->depth; s++) {
	  if(mask && MRIgetVoxVal(mask,c,r,s,0)<0.5) continue;
          v = MRIgetVoxVal(vol,c,r,s,f) - vmean;
	  MRIsetVoxVal(out,c,r,s,f,v);
	}//s
      }//r
    }//s
  }// f
  return(out);
}

/*!
 \fn int MRIloganize(MATRIX **X, MRI **Ct, MRI **intCt, const MATRIX *t, const double tstar)
 \brief Removes time points less than tstar. X, Ct, and intCt are changed
*/
int MRIloganize(MATRIX **X, MRI **Ct, MRI **intCt, const MATRIX *t, const double tstar)
{
  int n,n0,nkeep=0;
  for(n = (*Ct)->nframes-1; n >= 0; n--){
    if(t->rptr[n+1][1] < tstar) break;
    nkeep++;
  }
  n0 = n + 1;
  printf("Loganize tstar = %g, n0 = %d, nkeep=%d\n",tstar,n0,nkeep);

  MATRIX *XX = MatrixAlloc(nkeep,1,MATRIX_REAL);
  int k=0;
  for(n=n0; n < (*Ct)->nframes; n++){
    XX->rptr[k+1][1] = (*X)->rptr[n+1][1];
    k++;
  }

  MRI *Ct2 = MRIallocSequence((*Ct)->width,(*Ct)->height,(*Ct)->depth,MRI_FLOAT,nkeep);
  MRIcopyHeader((*Ct),Ct2);
  MRIcopyPulseParameters((*Ct),Ct2);
  MRI *intCt2 = MRIallocSequence((*Ct)->width,(*Ct)->height,(*Ct)->depth,MRI_FLOAT,nkeep);
  MRIcopyHeader((*Ct),intCt2);
  MRIcopyPulseParameters((*Ct),intCt2);
  int c,r,s;
  double v;
  for(c=0; c < (*Ct)->width; c++){
    for(r=0; r < (*Ct)->height; r++){
      for(s=0; s < (*Ct)->depth; s++){
	k=0;
	for(n=n0; n < (*Ct)->nframes; n++){
	  v = MRIgetVoxVal((*Ct),c,r,s,n);
	  MRIsetVoxVal(Ct2,c,r,s,k,v);
	  v = MRIgetVoxVal((*intCt),c,r,s,n);
	  MRIsetVoxVal(intCt2,c,r,s,k,v);
	  k++;
	}
      }
    }
  }
  MRIfree(Ct);
  MRIfree(intCt);
  MatrixFree(X);
  *Ct = Ct2;
  *intCt = intCt2;
  *X = XX;
  printf("nframes = %d\n",(*Ct)->nframes);

  return(0);
}




