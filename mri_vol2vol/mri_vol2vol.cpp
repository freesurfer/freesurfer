/**
 * @brief converts values in one volume to another volume
 *
 * Resamples a volume into another field-of-view using various types of 
 * matrices (FreeSurfer, FSL, SPM, and MNI). This is meant to be used
 * in conjunction with tkregister2.
 *
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

/*
BEGINUSAGE --------------------------------------------------------------

mri_vol2vol

  --mov  movvol       : input (or output template with --inv)
  --targ targvol      : output template (or input with --inv)
  --o    outvol       : output volume
  --disp dispvol      : displacement volume
  --downsample N1 N2 N3 : downsample input volume (do not include a targ or regsitration)
         sets --fill-average, --fill-upsample 2, and --regheader

  --reg  register.dat : tkRAS-to-tkRAS matrix   (tkregister2 format)
  --lta  register.lta : Linear Transform Array (usually only 1 transform)
  --lta-inv  register.lta : LTA, invert (may not be the same as --lta --inv with --fstal)
  --fsl  register.fsl : fslRAS-to-fslRAS matrix (FSL format)
  --xfm  register.xfm : ScannerRAS-to-ScannerRAS matrix (MNI format)
  --regheader         : ScannerRAS-to-ScannerRAS matrix = identity
  --mni152reg         : target MNI152 space (need FSL installed)
  --s subject         : set matrix = identity and use subject for any templates

  --inv               : sample from targ to mov

  --tal               : map to a sub FOV of MNI305 (with --reg only)
  --talres resolution : set voxel size 1mm or 2mm (def is 1)
  --talxfm xfmfile    : default is talairach.xfm (looks in mri/transforms)

  --m3z morph    : non-linear morph encoded in the m3z format
  --noDefM3zPath : flag indicating that the code should not be looking for 
       the non-linear m3z morph in the default location (subj/mri/transforms), but should use 
       the morph name as is
  --inv-morph    : compute and use the inverse of the m3z morph

  --fstarg <vol>      : optionally use vol from subject in --reg as target. default is orig.mgz 
  --crop scale        : crop and change voxel size
  --slice-crop sS sE  : crop output slices to be within sS and sE
  --slice-reverse     : reverse order of slices, update vox2ras
  --slice-bias alpha  : apply half-cosine bias field

  --trilin            : trilinear interpolation (default)
  --nearest           : nearest neighbor interpolation
  --cubic             : cubic B-Spline interpolation
  --interp interptype : interpolation cubic, trilin, nearest (def is trilin)
  --fill-average      : compute mean of all source voxels in a given target voxel
  --fill-conserve     : compute sum  of all source voxels in a given target voxel
  --fill-upsample USF : source upsampling factor for --fill-xxx (default is 2)

  --mul mulval   : multiply output by mulval

  --vsm vsmvol <pedir> : Apply a voxel shift map. pedir: +/-1=+/-x, +/-2=+/-y, +/-3=+/-z (default +2)
  --vsm-pedir pedir : phase encode direction for vsm

  --precision precisionid : output precision (def is float)
  --keep-precision  : set output precision to that of input
  --kernel            : save the trilinear interpolation kernel instead

  --gcam mov srclta gcam dstlta vsm interp out
     srclta, gcam, or vsm can be set to 0 to indicate identity
     direction is automatically determined from srclta and dstlta
     interp 0=nearest, 1=trilin, 5=cubicbspline

  --spm-warp mov movlta warp interp output
     mov is the input to be mapped 
     movlta maps mov to the vbm input space (use 0 to ignore)
       if movlta=0, then input is anything that shares a RAS space with the VBM input
     warp is typically y_rinput.nii
     interp 0=nearest, 1=trilin

  --no-resample : do not resample, just change vox2ras matrix
  --no-resample-scale : do not resample, just change vox2ras matrix, using scale=voxsize

  --rot   Ax Ay Az : rotation angles (deg) to apply to reg matrix
  --trans Tx Ty Tz : translation (mm) to apply to reg matrix
  --shear Sxy Sxz Syz : xz is in-plane
  --reg-final regfinal.dat : final reg after rot and trans (but not inv)

  --synth : replace input with white gaussian noise
  --seed seed : seed for synth (def is to set from time of day)

  --save-reg : write out output volume registration matrix

  --help : go ahead, make my day
  --debug
  --version



ENDUSAGE ---------------------------------------------------------------
*/

/*
BEGINHELP --------------------------------------------------------------

Resamples a volume into another field-of-view using various types of
matrices (FreeSurfer, FSL, SPM, and MNI). This is meant to be used
in conjunction with tkregister2.

FLAGS AND ARGUMENTS

--mov movvol

This volume must have the same geometry as the --mov volume passed to
tkregister2 when creating/checking the registration file. By default,
this will be the input volume that will be resampled. If --inv is
specified, then this will become the geometry template for the output
instead.

--targ targvol

This volume must have the same geometry as the --targ volume passed to
tkregister2 when creating/checking the registration file. By default,
this will be the volume will be the geometry template for the output.
If --inv is specified, then this becomes the input volume that will be
resampled instead. The target volume can be implicitly specified with
--tal or --fstarg.

--reg register.dat

This simple text file contains the freesurfer registration matrix. It
is the same as the file passed to and generated by tkregister2 with
the --reg flag. If --tal or --fstarg is specified, then the subject
is obtained from the regfile.

--fsl register.fsl

Registration matrix created with the FSL flirt program using targ as
the reference and mov as input. Note: you cannot use any of the files
from $FSLDIR/etc/standard as mov or targ. These volumes do not have
geometry information in them, and FreeSurfer and FSL will default to
different things. Same as in tkregister2.

--xfm register.xfm

MNI-style registration matrix (eg, like one created with mritotal).
This matrix maps from mov Scanner-RAS to targ Scanner-RAS, where
'Scanner-RAS' is the vox2ras matrix as found in each file.
Same as in tkregister2.

--regheader

Create a registration matrix assuuming that the mov Scanner-RAS and
targ Scanner-RAS are the same. This is the same as using a register.xfm
with the identity matrix in it. This can be used with some SPM
registrations (which change only the matrix in the .mat file).
Same as in tkregister2.

--mni152reg 

Target MNI152 space. If the mov volume is in the native space of an
individual, then also supply a registration (--reg). This registration
is concatenated with that in subject/mri/transforms/reg.mni152.2mm.dat
(created with mni152reg) to produce a registration from the mov vol
to MNI152 (defined by $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz).
If the data are in fsaverage 2mm space, then do not supply a --reg.
Instead, $FREESURFER_HOME/average/mni152.register.dat is used
as the registration. Do not supply a target volume with --mni152reg.

--inv

Invert the transform. The movvol becomes the geometry template for the
output, and the targvol becomes the input that will be resampled.

--o outvol

Output volume.  By default, this will be the movvol resmapled into the
targvol space (and so will have the same geometry as the targvol). If
--inv is specified, then this will be the targvol resmapled into the
movvol space (and so will have the same geometry as the movvol). By
default, the output volume will be float, but this can be changed
with --precision. By default, the interpolation will be done with
trilinear, but this can be changed with --interp.

--tal

Resample the movvol to talairach (ie, MNI305) space. The talairach
matrix is obtained from talairach.xfm from
SUBJECTS_DIR/subjid/transforms. subjid is read from the register.dat
file. Requires --reg. Do not specify --targ as the target volume is
implicitly set to $FREESURFER_HOME/average/mni305.cor.subfovV.mgz,
where V is either 1 (for 1mm) or 2 (for 2mm). 2mm is used by default,
but this can be changed with --talres.  mni305.cor.subfovV.mgz the
MNI305 (1mm or 2mm isotropic) volume in a reduced FOV that covers only
the brain. Reducing the FOV saves space relative to the 256^3 COR FOV.
The transformation matrix is computed as R*inv(Xtal)*inv(Rtal), where
Xtal is talairach.xfm matrix, R is the matrix in the regfile, and Rtal
maps from the talairach COR FOV to the SubFOV (mni305.cor.subfovV.reg).
If you want to sample the targvol from talairach space into the movvol
space, then specify --inv. SUBJECTS_DIR is read from the environment
or can be specified with --sd.

--fstalres resmm

Set the resolution of the output when using --fstal. By default, it
is 2 mm, but can be changed to 1.0 mm with --fstalres 1

--fstarg <vol>

Set target to vol from the subject found in register.dat
file. If vol is not specified, uses orig.mgz. Requires --reg.  
Same as tkregister2.

--crop scale

Crop mov volume down to minimum size to fit non-zero voxels. The size of
the voxels is reduced by scale (ie, --crop 2 would crop and reduce the
voxel size by a factor of 2, eg 1.0 mm becomes 0.5 mm).

--slice-crop start end

Crop output volume to be within slices start and end. The geometry is 
updated to reflect the new limits.

--interp method

Interpolate the output based on the given method. Legal values are:
cubic, trilin and nearest. trilin is the default. Can also use
--cubic, --trilin or --nearest.

--soap soap_ctl_point_fname num_iters

perform soap bubble smoothing on the input volume, using the volume specified by soap_ctl_point_fname
as the fixed (control) points for num_iters iterations

--precision precisionid

Set output precision to precisionid. Legal values are uchar, short,
int, long, and float. Default is float.

--keep-precision  : set output precision to that of input

--kernel

Save the trilinear interpolation kernel at each voxel instead of the
interpolated image.

--nomr

Don't copy the template MR parameters, but instead preserve the input volume 
ones

--help

Prints out all this information.

--gdiagno diagnostic level

Sets the diagnostic level (only good for debuggin').

--version

Print out version string and exit.


EXAMPLES:

Below are some exampls of how one might use mri_vol2vol. They are not
exhaustive of all the possible combinations of options. Typically, one
uses a template to establish the registration, then resamples data
that are in correspondence with the template.

1. If a functional volume is f.bhdr (or f.nii.gz, or f.mgh, etc), and the
subject is bert, and the registration file is register.dat, then
running the following command should show that they are in
registration:

tkregister2 --reg register.dat --mov f.nii.gz

If they are not, then fix it because nothing below is going to work. You
can also check the registration with:

tkmedit bert orig.mgz  -overlay f.nii.gz -overlay-reg register.dat

The register.dat will look something like this
----------------- register.dat --------------------------
bert
3.125
5.000
0.150000
1.000000e+00 0.000000e+00 0.000000e+00 -2.252487e+00
0.000000e+00 -8.902127e-01 4.555448e-01 2.342102e+00
0.000000e+00 4.555449e-01 8.902128e-01 -2.159538e-01
0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00
round
----------------- register.dat --------------------------

1.A. To resample the functional into anatomical space:

mri_vol2vol --reg register.dat --mov f.nii.gz --fstarg \\
   --o f-in-anat.mgh

f-in-anat.mgh will have the same size and geometry as
bert/mri/orig.mgz.  You can test the result in two ways:

  # This will show the low-res functional alighned with its resampled self
  tkregister2 --reg register.dat --mov f.nii.gz --targ f-in-anat.mgh

  # This will show the resampled functional aligned with the anatomical
  tkregister2 --mov f-in-anat.mgh --targ $SUBJECTS_DIR/bert/mri/orig.mgz \\
     --regheader --reg /tmp/reg

1.B. To resample the anatomical into the functional space. This is
basically the same command line as 1.A, but --inv has been included
and the name of the output is changed.

mri_vol2vol --reg register.dat --mov f.nii.gz --fstarg \\
  --o anat-in-func.mgh --inv

anat-in-func.mgh will be the same size and geometry as f.nii.gz.
You can test the result in two ways:

  # This will show the low-res anat aligned with its hires self
  tkregister2 --reg register.dat --mov anat-in-func.mgh

  # This will show the resampled anat aligned with the functional
  tkregister2 --mov anat-in-func.mgh --targ f.nii.gz \\
     --regheader --reg /tmp/reg

1.C Map functional to anatomical without resampling. Rather, change
the vox2ras (sform/qform) matrix. This is the same cmd line as 1.A,
but --no-resample as been added.

mri_vol2vol --reg register.dat --mov f.nii.gz --fstarg \\
   --o f.new.vox2ras.nii.gz --no-resample

f.new.vox2ras.nii.gz will have the same dimension and voxel size
as f.nii.gz, but its vox2ras (sform/qform) matrix will have changed.
You can check the registration in two ways:

  # The registration is created implicitly from the vox2ras matrix
  # (that is what --regheader does). There's no need to specify
  # and input registration
  tkregister2 --mov f.new.vox2ras.nii.gz --s bert --regheader --reg /tmp/reg

  # Display the functional as an overlay in tkmedit (no registration
  # needed).
  tkmedit bert orig.mgz -overlay f.new.vox2ras.nii.gz

1.D Map a binary mask in functional space to anatomical space. This is
basically the same cmd line as 1.A, but --interp nearest has been
added so that it does not try to interpolate the mask (ie, it will
still be binary after resampling):

mri_vol2vol --reg register.dat --mov mask.nii.gz --fstarg \\
   --o mask-in-anat.mgh --interp nearest

2. Map functional to/from talairach (MNI305) space. This uses a
two-stage registration: func-to-anat (register.dat) and
anat-to-talairach (talairach.xfm).

Make sure that sure the func-to-anat reg is correct as was done in
Example 1. Next, make sure that the anat-to-tal is correct with:

tkregister2 --s bert --fstal

2.A Map functional to talairach (MNI305) space with 2mm isotropic
resolution. This is very similar to 1.A with the addition of --tal
and --talres 2.

mri_vol2vol --mov f.nii.gz --reg register.dat \\
     --o f-in-tal.2mm.mgh --tal --talres 2

f-in-tal.2mm.mgh will be 2mm isotropic with the same geometry as
$FREESURFER_HOME/average/mni305.cor.subfov2.mgz. This command will
also create f-in-tal.2mm.mgh.reg, which will register the volume with
any average MNI305 FreeSurfer subject (fsaverage is used by default).
The resampling can be checked with:

  # This will show the functional with the fsaverage anatomical
  tkregister2 --mov f-in-tal.2mm.mgh --reg f-in-tal.2mm.mgh.reg

2.B Map functional to talairach (MNI305) space with 1mm isotropic
resolution. Same as 2.A but use --talres 1.

mri_vol2vol --mov f.nii.gz --reg register.dat \\
     --o f-in-tal.1mm.mgh --tal --talres 1

f-in-tal.1mm.mgh will take up 8 times as much space as f-in-tal.2mm.mgh

3. Apply an MNI transform to data by resampling the anatomical orig.mgz
into talairach space using bert/mri/transforms/talairach.xfm:

First, check that the talairach.xfm is correct (this is basically the same
thing as 'tkregister2 --s bert --fstal' in Example 2:

 cd bert/mri
 tkregister2 --targ orig.mgz \\
     --mov $FREESURFER_HOME/average/mni305.cor.mgz \\
     --xfm transforms/talairach.xfm --reg /tmp/reg

 Now resample
 mri_vol2vol --mov orig.mgz \\
     --targ $FREESURFER_HOME/average/mni305.cor.mgz \\
     --xfm transforms/talairach.xfm  \\
     --o orig-in-mni305.mgz

 Now test the resampling:
 tkregister2 --mov orig-in-mni305.mgz \\
    --targ $FREESURFER_HOME/average/mni305.cor.mgz \\
    --reg /tmp/reg --regheader


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

mri_convert, tkregister2


ENDHELP --------------------------------------------------------------

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h>
#include <errno.h>

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
#include "transform.h"
#include "gca.h"
#include "gcamorph.h"
#include "fio.h"
#include "pdf.h"
#include "cmdargs.h"
#include "mri_circulars.h"
#include "mriBSpline.h"
#include "timer.h"
#include "mrinorm.h"

#ifdef X
#undef X
#endif

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
static int isflag(char *flag);
static int nth_is_arg(int nargc, char **argv, int nth);
#include "tags.h"
static int istringnmatch(const char *str1, const char *str2, int n);
static MATRIX *LoadRtal(int talres);
MATRIX *LoadRfsl(char *fname);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

int debug = 0, gdiagno = -1;

static int soap_bubble_iters = 0 ;
static MRI *mri_soap_ctrl = NULL ;

char *movvolfile=NULL;
char *targvolfile=NULL;
int  fstarg = 0;
const char *fstargfile = "orig.mgz";
char *outvolfile=NULL, *outdir=NULL;
char *regfile=NULL;
char *xfmfile=NULL;
char *fslregfile=NULL;
char *tempvolfile=NULL;
int  invert=0;
int  fstal=0;
LTA  *lta=NULL;
int  usedltageom=0;
int  fstalres = 2; // Can only be 1 or 2
const char *precision = "float";
int   precisioncode = MRI_FLOAT;
const char *interpmethod = "trilinear";
int   interpcode = 0;
int   sinchw;
int   regheader=0;
int   noresample=0;

MRI *mov, *targ, *out;
MRI *in, *mri_template;
MRI *tmpmri;

MATRIX *R=NULL, *R2=NULL, *invR, *XFM;
MATRIX *vox2vox, *vox2ras;
MATRIX *Xtal,*invXtal,*Rtal,*invRtal;
MATRIX *Tin, *invTin, *Sin, *invSin;
MATRIX *Ttemp, *invTtemp, *Stemp, *invStemp;
MATRIX *Rfsl, *Rfsl2;

char *FSH=NULL;
char *SUBJECTS_DIR=NULL;
const char *talxfmfile = "talairach.xfm";

char *talsubject = NULL;
char *subject = NULL;
const char *MNIsubject = "21_vc716";
const char *subject_outreg = NULL;

int dont_irescale = 1;
float minrescale = 0.0, maxrescale = 255.0;

float ipr, bpr, intensity;
int float2int,err, nargs;
int SaveReg=0;

int DoKernel = 0;
int DoSaveInputMR = 1 ; // this is now the default behavior
int DoDelta  = 0;

char tmpstr[2000];

int DoMorph = 0;
int InvertMorph = 0;
TRANSFORM *Rtransform;  //types : M3D, M3Z, LTA, FSLMAT, DAT, OCT(TA), XFM
GCAM      *gcam;
GCAM      *MNIgcam;
char gcamfile[1000];
char MNIgcamfile[1000];
MRI_REGION region;
const char *m3zfile = "talairach.m3z";

double angles[3] = {0,0,0};
MATRIX *Mrot = NULL;
double xyztrans[3] = {0,0,0};
MATRIX *Mtrans = NULL;
double shear[3] = {0,0,0};
MATRIX *Mshear = NULL;

char *SegRegCostFile = NULL;
char  *fspec;
MRI *regseg;
int CostOnly = 0;
char *RegFileFinal=NULL;

int SynthSeed = -1;
int synth = 0;
int DoCrop = 0;
double CropScale = 0;

char *DispFile = NULL;
MRI *DispMap = NULL;

int slice_crop_flag = 0;
int slice_crop_start, slice_crop_stop;
int SliceReverse = 0;
int SliceBias  = 0;
double SliceBiasAlpha = 1.0;

int useold = 1;
MRI *vsm = NULL;
char *vsmvolfile=NULL;

int defM3zPath = 1; // use default path to the m3z file
int TargMNI152 = 0;
int keepprecision = 0;
int DoFill=0;
int DoFillConserve=0;
int FillUpsample=2;
MRI *MRIvol2volGCAM(MRI *src, LTA *srclta, GCA_MORPH *gcam, LTA *dstlta, MRI *vsm, int sample_type, MRI *dst, int pedir=2);
MRI *MRIvol2volGCAM0(MRI *src, LTA *srclta, GCA_MORPH *gcam, LTA *dstlta, MRI *vsm, int sample_type, MRI *dst);
int DoMultiply=0;
double MultiplyVal=0;
int DownSample[3] = {0,0,0}; // downsample source
int pedir = 2; // for VSM 1=x, 2=y, 3=z

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  char regfile0[1000];
  double costs[8];
  FILE *fp;
  int n,err;
  MRI *crop, *cropnew, *mri;
  MRI_REGION box;
  LTA *ltareg;

  vg_isEqual_Threshold = 10e-4;

  nargs = handleVersionOption(argc, argv, "mri_vol2vol");
  if(nargs && argc - nargs == 1) exit (0);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  FSH = getenv("FREESURFER_HOME");
  if(FSH==NULL) {
    printf("ERROR: FREESURFER_HOME undefined.\n");
    exit(1);
  }
  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  if(gdiagno > -1) Gdiag_no = gdiagno;
  check_options();

  // Seed the random number generator just in case
  if (SynthSeed < 0) SynthSeed = PDFtodSeed();
  srand48(SynthSeed);

  dump_options(stdout);

  if(DoCrop){
    //printf("\n"); 
    //printf("DoCrop \n"); 
    printf("Crop %lf\n",CropScale);
    mov = MRIread(movvolfile);
    if(mov == NULL) exit(1);
    err = MRIboundingBox(mov, 0, &box);
    if(err) exit(1);
    crop  = MRIcrop(mov,box.x, box.y, box.z,
		    box.x+box.dx, box.y+box.dy, box.z+box.dz);
    //MRIwrite(crop,"crop.mgh");

    cropnew = MRIalloc(nint(crop->width*CropScale),
		       nint(crop->height*CropScale),
		       nint(crop->depth*CropScale),
		       mov->type);
    cropnew->x_r = crop->x_r;
    cropnew->x_a = crop->x_a;
    cropnew->x_s = crop->x_s;

    cropnew->y_r = crop->y_r;
    cropnew->y_a = crop->y_a;
    cropnew->y_s = crop->y_s;

    cropnew->z_r = crop->z_r;
    cropnew->z_a = crop->z_a;
    cropnew->z_s = crop->z_s;

    cropnew->c_r = crop->c_r;
    cropnew->c_a = crop->c_a;
    cropnew->c_s = crop->c_s;

    cropnew->xsize = crop->xsize/CropScale;
    cropnew->ysize = crop->ysize/CropScale;
    cropnew->zsize = crop->zsize/CropScale;

    printf("vol2vol\n");
    err = MRIvol2Vol(crop,cropnew,NULL,SAMPLE_NEAREST,0);
    if(err) exit(1);
    printf("Saving\n");
    err = MRIwrite(cropnew,outvolfile);
    printf("#VMPC# mri_vol2vol VmPeak %d\n",GetVmPeak());
    printf("mri_vol2vol done\n");
    exit(err);
  }

  /*-----------------------------------------------------*/
  if(fstal) {
    // Recompute R for converting to/from talairach space
    // and set the target volume file
    printf("\n"); 
    printf("Compute R for talairach space\n");
    Xtal = DevolveXFM(subject, NULL, talxfmfile);
    invXtal = MatrixInverse(Xtal,NULL);
    if(Xtal == NULL) exit(1);
    if(fstalres > 0) Rtal = LoadRtal(fstalres);
    else             Rtal = MatrixIdentity(4,NULL);
    invRtal = MatrixInverse(Rtal,NULL);
    if(1 || Gdiag_no > 0) {
      printf("matrix from regfile ---------------------- \n");
      MatrixPrint(stdout,R);
      printf("Xtal ---------------------- \n");
      MatrixPrint(stdout,Xtal);
      printf("Rtal ---------------------- \n");
      MatrixPrint(stdout,Rtal);
    }
    // Recompute: R = R*inv(Xtal)*inv(Rtal)
    R = MatrixMultiply(R,invXtal,R);
    R = MatrixMultiply(R,invRtal,R);
  }

  if(!invert) {
    // dont invert
    //printf("\n"); 
    //printf("Don't invert!\n"); 
    mov = MRIread(movvolfile);
    if (mov == NULL) exit(1);
    if (targvolfile != NULL ) targ = MRIreadHeader(targvolfile,MRI_VOLUME_TYPE_UNKNOWN);
    else if (lta != NULL && !fstal)
    {
       targ = MRIclone(mov,targ);
       MRIcopyVolGeomToMRI(targ,&lta->xforms[0].dst); 
       usedltageom = 1;
    }
    else {
      // Downsample
      targ = MRIallocHeader(mov->width,mov->height,mov->depth,MRI_FLOAT,mov->nframes);
      MRIcopyHeader(mov, targ);
      MRIcopyPulseParameters(mov, targ);
      targ->width  = ceil( mov->width  / DownSample[0]);
      targ->height = ceil( mov->height / DownSample[1]);
      targ->depth  = ceil( mov->depth  / DownSample[2]);
      targ->xsize *= DownSample[0];
      targ->ysize *= DownSample[1];
      targ->zsize *= DownSample[2];
      targ->xstart = -targ->xsize*targ->width / 2;
      targ->xend   =  targ->xsize*targ->width / 2;
      targ->ystart = -targ->ysize*targ->height/ 2;
      targ->yend   =  targ->ysize*targ->height/ 2;
      targ->zstart = -targ->zsize*targ->depth / 2;
      targ->zend   =  targ->zsize*targ->depth / 2;
    }
    if (targ == NULL) exit(1);
    in = mov;
    mri_template = targ;
    tempvolfile = targvolfile;
  }
  else{
    //invert
    printf("\n"); 
    printf("Invert!\n"); 
    if (targvolfile != NULL ) targ = MRIread(targvolfile);
    if(targ == NULL) exit(1);
    if (movvolfile != NULL) mov = MRIreadHeader(movvolfile,MRI_VOLUME_TYPE_UNKNOWN);
    else if (lta != NULL && !fstal)
    {
       mov = MRIclone(targ,mov);
       MRIcopyVolGeomToMRI(mov,&lta->xforms[0].src); 
       usedltageom = 1;
    }    
    if(mov == NULL) exit(1);
    in = targ;
    mri_template = mov;
    tempvolfile = movvolfile;
  }
  if (mri_soap_ctrl)
  {
    MRI *mri_out ;

    MRIbinarize(mri_soap_ctrl, mri_soap_ctrl, 1, 0, CONTROL_MARKED) ;
    mri_out = MRIsoapBubble(mov, mri_soap_ctrl, NULL, soap_bubble_iters, 0);
    MRIwrite(mri_out, outvolfile) ;
    printf("#VMPC# mri_vol2vol VmPeak %d\n",GetVmPeak());
    exit(0) ;
  }
  if(synth) {
    printf("\n"); 
    printf("Replacing input data with synthetic white noise\n");
    MRIrandn(in->width,in->height,in->depth,in->nframes,0,1,in);
  }

  if(regheader) {
    printf("\n"); 
    printf("Computing registration based on scanner-to-scanner\n");
    R = MRItkRegMtx(targ,mov,XFM);
  }

  if(fslregfile) {
    printf("\n"); 
    printf("Computing registration based on fsl registration\n");
    R = MRIfsl2TkReg(targ, mov, Rfsl);
    if(Rfsl2) R2 = MRIfsl2TkReg(targ, mov, Rfsl2);
  }

  if(R == NULL)
    ErrorExit(ERROR_BADPARM, "ERROR: no registration (R) is specified\n") ;

  //printf("\n"); 
  //printf("Registration:\n");
  //MatrixPrint(stdout,R);
  //printf("\n");

  if(R2){
    R = MatrixSubtract(R,R2,R);
    for(n=1; n<4; n++) R->rptr[n][n] += 1.0;
  }

  if(Mrot){
    printf("Applying rotation matrix (R=M*R)\n");
    printf("Current Reg Matrix is:\n");
    MatrixPrint(stdout,R);
    printf("  Angles (deg): %lf %lf %lf\n",angles[0]*180/M_PI,angles[1]*180/M_PI,angles[2]*180/M_PI);
    printf("  Angles (rad): %lf %lf %lf\n",angles[0],angles[1],angles[2]);
    printf("  Rotation matrix:\n");
    MatrixPrint(stdout,Mrot);
    R = MatrixMultiply(Mrot,R,R);
  }

  if(Mtrans){
    printf("Applying translation matrix (R=M*R)\n");
    printf("Current Reg Matrix is:\n");
    MatrixPrint(stdout,R);
    printf("  Trans (mm): %lf %lf %lf\n",xyztrans[0],xyztrans[1],xyztrans[2]);
    printf("  Translation matrix:\n");
    MatrixPrint(stdout,Mtrans);
    R = MatrixMultiply(Mtrans,R,R);
  }
  if(Mshear){
    printf("Applying shear matrix (R=M*R)\n");
    printf("Current Reg Matrix is:\n");
    MatrixPrint(stdout,R);
    printf("  Shear: %lf %lf %lf\n",shear[0],shear[1],shear[2]);
    printf("  Shear matrix:\n");
    MatrixPrint(stdout,Mshear);
    R = MatrixMultiply(Mshear,R,R);
  }

  if(RegFileFinal){
    printf("Writing final tkRAS-to-tkRAS Matrix to %s\n",RegFileFinal);
    regio_write_register(RegFileFinal,subject,in->xsize,
			 in->zsize,1,R,FLT2INT_ROUND);
  }

  if(invert) {
    printf("Inverting registration\n");
    R = MatrixInverse(R,NULL);
  }

  printf("\n");
  printf("Final tkRAS-to-tkRAS Matrix is:\n");
  MatrixPrint(stdout,R);
  printf("\n");

  if(DispFile){
    printf("Computing affine displacment\n");
    DispMap = MRIaffineDisplacment(in, R);
    MRIwrite(DispMap,DispFile);
    if(outvolfile == NULL) {
      printf("No other output specified, so exiting now\n");
      printf("#VMPC# mri_vol2vol VmPeak %d\n",GetVmPeak());
      exit(0);
    }
  }

  invR = MatrixInverse(R,NULL);

  // Vox-to-tkRAS Matrices
  Tin      = MRIxfmCRS2XYZtkreg(in);
  invTin   = MatrixInverse(Tin,NULL);
  Ttemp    = MRIxfmCRS2XYZtkreg(mri_template);
  invTtemp = MatrixInverse(Ttemp,NULL);

  // Vox-to-ScannerRAS Matrices
  Sin      = MRIxfmCRS2XYZ(in,0);
  invSin   = MatrixInverse(Sin,NULL);
  Stemp    = MRIxfmCRS2XYZ(mri_template,0);
  invStemp = MatrixInverse(Stemp,NULL);

  if(noresample) {
    printf("Not resampling, only changing vox2ras matrix\n");
    // Compte new vox2ras instead of resampling
    // vox2ras = Stemp * invTtemp * invR * Tin
    vox2ras = MatrixMultiply(Stemp,invTtemp,NULL);
    MatrixMultiply(vox2ras,invR,vox2ras);
    MatrixMultiply(vox2ras,Tin,vox2ras);
    MRIsetVoxelToRasXform(in,vox2ras);
    if(slice_crop_flag){
      printf("Cropping slices from %d to %d\n",slice_crop_start,slice_crop_stop);
      crop  = MRIcrop(in, 0, 0, slice_crop_start,
		      in->width-1, in->height-1,slice_crop_stop);
      if(crop == NULL) exit(1);
      MRIfree(&in);
      in = crop;
    }
    if(SliceReverse){
      printf("Reversing slices, updating vox2ras\n");
      mri = MRIreverseSlices(in, NULL);
      if(mri == NULL) exit(1);
      MRIfree(&in);
      in = mri;
    }
    if(SliceBias){
      printf("Applying Half-Cosine Slice Bias, Alpha = %g\n",SliceBiasAlpha);
      MRIhalfCosBias(in, SliceBiasAlpha, in);
    }
    err = MRIwrite(in,outvolfile);
    if(err){
      printf("ERROR: writing %s\n",outvolfile);
      exit(1);
    }
    printf("To check registration, run:\n");
    printf("\n");
    printf("  tkregister2 --mov %s --targ %s --regheader --reg /tmp/reg \n",
           outvolfile,tempvolfile);
    printf("\n");
    printf("mri_vol2vol done\n");
    return(0);
  }

  // Only gets here if resampling
  // vox2vox converts a template vox to input vox
  // vox2vox = invTin * R * Ttemp
  vox2vox = MatrixMultiply(invTin,R,NULL);
  MatrixMultiply(vox2vox,Ttemp,vox2vox);

  printf("\n");
  printf("Vox2Vox Matrix is:\n");
  MatrixPrint(stdout,vox2vox);
  printf("\n");

  // Allocate the output
  mri_template->type = precisioncode;
  if (DoSaveInputMR)  // it is now on by default
  {
    mri_template->tr = in->tr ;
    mri_template->ti = in->ti ;
    mri_template->flip_angle = in->flip_angle ;
    mri_template->te = in->te ;

  }
  if (!DoMorph) {
    if(DoKernel) {
      out = MRIcloneBySpace(mri_template,MRI_FLOAT,8);
      printf("Computing Trilinear Kernel\n");
      MRIvol2VolTLKernel(in,out,vox2vox);
    } else if(DoDelta) {
      printf("Computing Delta\n");
      out = MRIvol2VolDelta(in,mri_template,R);
    } 
    else if(DoFill){
      printf("Running MRIvol2VolFill(), DoConserve=%d, US=%d\n",DoFillConserve,FillUpsample);
      if(! invert) ltareg = TransformRegDat2LTA(targ, mov, R);
      else   	   ltareg = TransformRegDat2LTA(mov, targ, R);
      out = MRIvol2VolFill(in, NULL, ltareg, FillUpsample, DoFillConserve, out);
    }
    else {
      out = MRIcloneBySpace(mri_template,-1,in->nframes);
      printf("Resampling\n");
      if(useold) MRIvol2Vol(in,out,vox2vox,interpcode,sinchw);
      if(!useold){
	if(vsmvolfile){
	  printf("Reading %s\n",vsmvolfile);
	  vsm = MRIread(vsmvolfile);
	  if(vsm == NULL) exit(1);
	}
	MRIvol2VolVSM(in,out,vox2vox,interpcode,sinchw,vsm,pedir);
      }
    }
  }
  else {
    Rtransform = (TRANSFORM *)calloc(sizeof(TRANSFORM),1);
    Rtransform->xform = (void *)TransformRegDat2LTA(mri_template, mov, R); // LZ: this is where the morphing goes wrong.....

    printf("Reading gcam\n");
    if (defM3zPath)
      sprintf(gcamfile,"%s/%s/mri/transforms/%s",
	      SUBJECTS_DIR,subject,m3zfile);
    else
      sprintf(gcamfile,"%s", m3zfile);

    if(! InvertMorph){
      //mri_vol2vol --mov orig.mgz --morph --s subject --o orig.morphed.mgz
      gcam = GCAMread(gcamfile);
      if(gcam == NULL) exit(1);
      printf("Applying reg to gcam\n");
      GCAMapplyTransform(gcam, Rtransform);  //voxel2voxel
      printf("Applying morph to input\n");
      out = GCAMmorphToAtlas(in, gcam, NULL, -1, interpcode);

      //sprintf(MNIgcamfile,"%s/transforms/talairach.m3z", fio_dirname(gcam->atlas.fname));
      //printf("The MNI gcam fname is: %s\n", MNIgcamfile);      
      //MNIgcam = GCAMread(MNIgcamfile);
      //out = GCAMmorphToAtlasToMNI(in, gcam, MNIgcam, NULL, -1, interpcode);
    }
    else{
      //mri_vol2vol --mov orig.morphed.mgz --inv-morph --s subject --o origB.mgz
      gcam = GCAMread(gcamfile);
      if(gcam == NULL) exit(1);
      {
	MRI *mri_tmp ;
	mri_tmp = MRIalloc(gcam->image.width, gcam->image.height, gcam->image.depth, MRI_FLOAT) ;
	useVolGeomToMRI(&gcam->image, mri_tmp);
	
	GCAMinvert(gcam, mri_tmp) ;
	MRIfree(&mri_tmp) ;
      }
      printf("Applying reg to gcam\n");
      GCAMapplyTransform(gcam, Rtransform);  //voxel2voxel
      if (0) { //(in->type != MRI_UCHAR){
	printf("Changing type to uchar\n");
	tmpmri = MRISeqchangeType(in, MRI_UCHAR, 0 , 255, 1);
	MRIfree(&in);
	in = tmpmri;
      }
      printf("Applying inverse morph to input\n");
      gcam->gca = gcaAllocMax(1, 1, 1,
			      in->width, in->height,
			      in->depth,
			      0, 0) ;
      out = GCAMmorphFromAtlas(in, gcam, NULL, interpcode);
    }
    if(out == NULL) exit(1);
    
    if(0){
    printf("Extracting region\n");
    region.x = 51;
    region.y = 0;
    region.z = 11;
    region.dx = 156;
    region.dy = 216;
    region.dz = 240;
    tmpmri = MRIextractRegion(out, NULL, &region) ;
    MRIfree(&out);
    out = tmpmri;
    }
  }

  if(SegRegCostFile){
    sprintf(tmpstr,"%s/%s/mri/regseg",SUBJECTS_DIR,subject);
    fspec = IDnameFromStem(tmpstr);
    regseg = MRIread(fspec);
    if(regseg == NULL) exit(1);
    free(fspec);
    SegRegCost(regseg,out,costs);
    fp = fopen(SegRegCostFile,"a");
    for(n=0; n<3; n++) fprintf(fp,"%7.3lf ",xyztrans[n]);
    for(n=0; n<3; n++) fprintf(fp,"%5.1lf ",angles[n]*180/M_PI);
    fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[0],costs[1],costs[2]); // WM  n mean std
    fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[3],costs[4],costs[5]); // CTX n mean std
    fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
    fprintf(fp,"\n");
    fclose(fp);
    if(CostOnly) exit(0);
  }

  if(slice_crop_flag){
    printf("Cropping slices from %d to %d\n",slice_crop_start,slice_crop_stop);
    crop  = MRIcrop(out, 0, 0, slice_crop_start,
		    out->width-1, out->height-1,slice_crop_stop);
    if(crop == NULL) exit(1);
    MRIfree(&out);
    out = crop;
  }
  if(DoMultiply) {
    printf("Multiplying by %lf\n",MultiplyVal);
    MRImultiplyConst(out, MultiplyVal, out);
  }

  if(mov->ct) out->ct = CTABdeepCopy(mov->ct);

  err = MRIwrite(out,outvolfile);
  if(err){
    printf("ERROR: writing %s\n",outvolfile);
    exit(1);
  }

  if(fstal) {
    R = Rtal;
    subject_outreg = "fsaverage";
  }
  else{
    R = MatrixIdentity(4,NULL);
    if(subject != NULL) subject_outreg = subject;
    else                subject_outreg = "subject-unknown";
    printf("Output registration matrix is identity\n");
  }

  if(SaveReg) {
    sprintf(regfile0,"%s.lta",outvolfile);
    printf("INFO: writing registration to %s\n",regfile0);
    if (lta) LTAfree(&lta);
    lta = LTAalloc(1, NULL) ;
    lta->xforms[0].sigma = 10000.000f ;
    lta->xforms[0].x0 = lta->xforms[0].y0 = lta->xforms[0].z0 = 0 ;
    strcpy(lta->subject, subject_outreg); 
    MatrixCopy(R, lta->xforms[0].m_L) ;
    lta->type = REGISTER_DAT;
    lta->fscale  = 1;
    LTAmodifySrcDstGeom(lta,in,mri_template);
    LTAchangeType(lta,LINEAR_VOX_TO_VOX);
    LTAwrite(lta,regfile0);    
    LTAfree(&lta);
    
     
    sprintf(regfile0,"%s.reg",outvolfile);
    printf("INFO: writing registration matrix to %s\n",regfile0);
    regio_write_register(regfile0,subject_outreg,out->xsize,
                         out->zsize,1,R,FLT2INT_ROUND);
    if (!fstal) {
      if (! usedltageom ){
        printf("To check registration, run:\n");
        printf("\n");
        printf("  tkregister2 --mov %s --targ %s --reg %s \n",
               outvolfile,tempvolfile,regfile0);
      }
    }  else {
      printf("To check registration, run:\n");
      printf("\n");
      printf("  tkregister2 --s %s --surf white --reg %s --mov %s \n",
             subject_outreg,regfile0,outvolfile);
    }
  }

  printf("\n");
  printf("mri_vol2vol done\n");

  return(0);
}


/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  int err;
  char tmp[1000];

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
    else if (!strcasecmp(option, "--tal"))      fstal = 1;
    else if (!strcasecmp(option, "--inv"))      invert = 1;
    else if (!strcasecmp(option, "--no-resample")) noresample = 1;
    else if (!strcasecmp(option, "--no-resample-scale")){
      if(nargc < 1) argnerr(option,1);
      setenv("FS_SetVoxToRasXform_Change_VoxSize",pargv[0],1);
      noresample = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--regheader")) regheader = 1;
    else if (!strcasecmp(option, "--kernel"))    DoKernel = 1;
    else if (!strcasecmp(option, "--nomr"))      DoSaveInputMR = 1;
    else if (!strcasecmp(option, "--mr"))        DoSaveInputMR = 0;
    else if (!strcasecmp(option, "--delta"))     DoDelta = 1;
    else if (!strcasecmp(option, "--no-save-reg"))  SaveReg = 0;
    else if (!strcasecmp(option, "--save-reg"))  SaveReg = 1;
    else if (!strcasecmp(option, "--cost-only"))  CostOnly = 1;
    else if (!strcasecmp(option, "--synth"))   synth = 1;
    else if (!strcasecmp(option, "--soap"))   
    {
      if (nargc < 2) argnerr(option,2);
      nargsused = 2;
      mri_soap_ctrl  = MRIread(pargv[0]) ;
      if (mri_soap_ctrl == NULL)
	ErrorExit(ERROR_NOFILE, "") ;
      sscanf(pargv[1],"%d",&soap_bubble_iters);
      printf("performing soap bubble smoothing using %s for %d iterations\n", pargv[0], soap_bubble_iters) ;
    }
    else if (!strcasecmp(option, "--new"))   useold = 0;
    else if (!strcasecmp(option, "--fill-average")){
      DoFill = 1;
      DoFillConserve = 0;
    }
    else if (!strcasecmp(option, "--fill-conserve")){
      DoFill = 1;
      DoFillConserve = 1;
    }
    else if (!strcasecmp(option, "--fill-upsample")){
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&FillUpsample);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--downsample")){
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%d",&DownSample[0]);
      sscanf(pargv[1],"%d",&DownSample[1]);
      sscanf(pargv[2],"%d",&DownSample[2]);
      DoFill = 1;
      DoFillConserve = 0;
      FillUpsample = 2;
      regheader = 1;
      nargsused = 3;
    }
    else if (!strcasecmp(option, "--morph")) {
      DoMorph = 1;
      fstarg = 1;
    } 
    else if (!strcasecmp(option, "--fstarg")){
      fstarg = 1;
      if(CMDnthIsArg(nargc, pargv, 0)){
        fstargfile = pargv[0];
        nargsused = 1;
      } 
      printf("fstargfile %s\n",fstargfile);
    }
    else if (!strcasecmp(option, "--inv-morph")) {
      DoMorph = 1;
      InvertMorph = 1;
      invert = 1;
    } else if (istringnmatch(option, "--m3z",0)) {
      if (nargc < 1) argnerr(option,1);
      m3zfile = pargv[0]; DoMorph = 1;
      nargsused = 1;
    } else if (istringnmatch(option, "--noDefM3zPath",0)) {
      defM3zPath = 0; // use the m3z file as it is; no assumed location
      if(R == NULL) R = MatrixIdentity(4,NULL); // as subjid is not neccesary any more
      printf("Using the m3z file as it is; no assumed location.\n");
    } else if (istringnmatch(option, "--mov",0)) {
      if (nargc < 1) argnerr(option,1);
      movvolfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--vsm",0)) {
      if(nargc < 1) argnerr(option,1);
      vsmvolfile = pargv[0];
      useold = 0;
      if(CMDnthIsArg(nargc, pargv, 1)) sscanf(pargv[1],"%d",&pedir);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--vsm-pedir")) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&pedir);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--targ",0)) {
      if (nargc < 1) argnerr(option,1);
      targvolfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--reg",0)) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      err = regio_read_register(regfile, &subject, &ipr, &bpr,
                                &intensity, &R, &float2int);
      printf("\n");
      printf("Matrix from regfile:\n");
      MatrixPrint(stdout,R);
      printf("\n");
      if (err) exit(1);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--lta",0) || istringnmatch(option, "--lta-inv",0)) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      if(stricmp(FileNameExtension(regfile, tmp), "LTA")){
	if (stricmp(FileNameExtension(regfile, tmp), "identity.nofile")){
	  printf("LTA registration file needs to have .lta extension!\n");
	  exit(1);        
	}
      }
      lta = LTAread(regfile) ;
      if(lta == NULL){
	printf("ERROR reading LTA %s !\n",regfile);        
	exit(1) ;
      }
      if (!lta->xforms[0].src.valid){
	printf("ERROR LTA %s has no valid src geometry!\n",regfile);        
	exit(1) ;       
      }
      if (!lta->xforms[0].dst.valid){
	printf("ERROR LTA %s has no valid dst geometry!\n",regfile);        
	exit(1) ; 
      }
      if(lta->subject[0]==0) strcpy(lta->subject, "subject-unknown"); 
      subject = (char *) calloc(strlen(lta->subject)+2,sizeof(char));
      strcpy(subject, lta->subject) ;
      if(istringnmatch(option, "--lta-inv",0)){
	printf("Inverting LTA\n");
	LTAinvert(lta,lta);
      }
      intensity = lta->fscale ;
      float2int = FLT2INT_ROUND ;
      LTAchangeType(lta, REGISTER_DAT);
      R = lta->xforms[0].m_L;
      ipr = lta->xforms[0].src.xsize ;
      bpr = lta->xforms[0].src.zsize ;
      printf("\n");
      printf("Matrix from LTA:\n");
      MatrixPrint(stdout,R);
      printf("\n");
      //err = regio_read_register(regfile, &subject, &ipr, &bpr,
      //                          &intensity, &R, &float2int);
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--mni152reg",0)) {
      TargMNI152 = 1;
      sprintf(tmpstr,"%s/data/standard/MNI152_T1_2mm.nii.gz",getenv("FSLDIR"));
      targvolfile = strcpyalloc(tmpstr);
    } 
    else if (istringnmatch(option, "--reg-final",0)) {
      if (nargc < 1) argnerr(option,1);
      RegFileFinal = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--s",0)) {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      R = MatrixIdentity(4,NULL);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--copy-ctab")) {
      setenv("FS_COPY_HEADER_CTAB","1",1);
    } 
    else if (istringnmatch(option, "--fsl",0) ||
               istringnmatch(option, "--fslreg",0)) {
      if(nargc < 1) argnerr(option,1);
      fslregfile = pargv[0];
      Rfsl = LoadRfsl(fslregfile);
      if(Rfsl == NULL) exit(1);
      nargsused = 1;
      if(nargc > 1 && !CMDisFlag(pargv[1])){
	Rfsl2 = LoadRfsl(pargv[1]);
	if(Rfsl2 == NULL) exit(1);
	nargsused++;
      }
    } else if (istringnmatch(option, "--xfm",0)) {
      if (nargc < 1) argnerr(option,1);
      xfmfile = pargv[0];
      err = regio_read_mincxfm(xfmfile, &XFM, NULL);
      if (err) exit(1);
      regheader = 1;
      nargsused = 1;
    } else if (istringnmatch(option, "--talxfm",0)) {
      if (nargc < 1) argnerr(option,1);
      talxfmfile = pargv[0];
      fstal = 1;
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--out",0) ||
               istringnmatch(option, "--o",0)) {
      if (nargc < 1) argnerr(option,1);
      outvolfile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--disp",0)){
      if(nargc < 1) argnerr(option,1);
      DispFile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--talres",8)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&fstalres);
      if (fstalres != 1 && fstalres != 2) {
        printf("ERROR: tal res %d invalid. Only use 1 or 2\n",fstalres);
        exit(1);
      }
      nargsused = 1;
    } else if (istringnmatch(option, "--crop",6)) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&CropScale);
      DoCrop = 1;
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
    else if(istringnmatch(option, "--slice-crop",12)) {
      if(nargc < 2) argnerr(option,2);
      slice_crop_flag = 1;
      sscanf(pargv[0],"%d",&slice_crop_start);
      sscanf(pargv[1],"%d",&slice_crop_stop);
      nargsused = 2;
    } 
    else if(istringnmatch(option, "--slice-reverse",0)) SliceReverse = 1;
    else if(istringnmatch(option, "--slice-bias",0)){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&SliceBiasAlpha);
      SliceBias = 1;
      nargsused = 1;
    }
    else if (istringnmatch(option, "--interp",8)) {
      if (nargc < 1) argnerr(option,1);
      interpmethod = pargv[0];
      nargsused = 1;
      if (!strcmp(interpmethod,"sinc") && nth_is_arg(nargc, pargv, 1)) {
        sscanf(pargv[1],"%d",&sinchw);
        nargsused ++;
      }
    } 
    else if (istringnmatch(option, "--trilin",6)) {
      interpmethod = "trilinear";
    } 
    else if (istringnmatch(option, "--nearest",7)) {
      interpmethod = "nearest";
    } 
    else if (istringnmatch(option, "--cubic",0)) {
      interpmethod = "cubic";
    } 
    else if (istringnmatch(option, "--precision",0)) {
      if (nargc < 1) argnerr(option,1);
      precision = pargv[0];
      precisioncode = MRIprecisionCode(precision);
      if (precisioncode < 0) {
        printf("ERROR: precision %s unrecognized\n",precision);
        printf("       legal values are uchar, short, int, long, and float\n");
        exit(1);
      }
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--keep-precision",0)) {
      keepprecision = 1;
    } 
    else if (!strcasecmp(option, "--seed")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&SynthSeed);
      synth = 1;
      nargsused = 1;
    } else if (istringnmatch(option, "--sd",4)) {
      if (nargc < 1) argnerr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    } else if (istringnmatch(option, "--rot",0)) {
      if (nargc < 3) argnerr(option,3);
      // Angles are in degrees
      sscanf(pargv[0],"%lf",&angles[0]);
      sscanf(pargv[1],"%lf",&angles[1]);
      sscanf(pargv[2],"%lf",&angles[2]);
      angles[0] *= (M_PI/180);
      angles[1] *= (M_PI/180);
      angles[2] *= (M_PI/180);
      Mrot = MRIangles2RotMat(angles);
      nargsused = 3;
    } else if (istringnmatch(option, "--trans",0)) {
      if (nargc < 3) argnerr(option,3);
      // Translation in mm
      sscanf(pargv[0],"%lf",&xyztrans[0]);
      sscanf(pargv[1],"%lf",&xyztrans[1]);
      sscanf(pargv[2],"%lf",&xyztrans[2]);
      Mtrans = MatrixIdentity(4,NULL);
      Mtrans->rptr[1][4] = xyztrans[0];
      Mtrans->rptr[2][4] = xyztrans[1];
      Mtrans->rptr[3][4] = xyztrans[2];
      nargsused = 3;
    } else if (istringnmatch(option, "--shear",0)) {
      if (nargc < 3) argnerr(option,3);
      // Shear
      sscanf(pargv[0],"%lf",&shear[0]);
      sscanf(pargv[1],"%lf",&shear[1]);
      sscanf(pargv[2],"%lf",&shear[2]);
      Mshear = MatrixIdentity(4,NULL);
      Mshear->rptr[1][2] = shear[0]; // xy/col-slice
      Mshear->rptr[1][3] = shear[1]; // xz/col-row - in-plane
      Mshear->rptr[2][3] = shear[2]; // yz/row-slice
      nargsused = 3;
    } else if ( !strcmp(option, "--gdiagno") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&gdiagno);
      nargsused = 1;
    } else if (istringnmatch(option, "--subject",0)) {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--cost",0)) {
      if (nargc < 1) argnerr(option,1);
      SegRegCostFile = pargv[0];
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--spm-warp",0)) {
      LTA *srclta;
      if(nargc < 5){
	printf("  --spm-warp mov movlta warp interp out\n");
	argnerr(option,6);
      }
      mov = MRIread(pargv[0]);
      if(mov == NULL) exit(1);
      if(strcmp(pargv[1],"0")!=0){
	printf("Loading source LTA %s\n",pargv[1]);
	srclta = LTAread(pargv[1]);
	if(srclta == NULL) exit(1);
      } else srclta = NULL;
      MRI *warp = MRIread(pargv[2]);
      if(warp == NULL) exit(1);
      sscanf(pargv[3],"%d",&interpcode);
      printf("Running MRIapplySpmVbmWarp() interp = %d\n",interpcode);
      out = MRIapplySpmWarp(mov, srclta, warp, 0, interpcode, NULL);
      if(out == NULL) exit(1);
      printf("Writing to %s\n",pargv[4]);
      err = MRIwrite(out,pargv[4]);
      if(err) exit(1);
      MRIfree(&warp);
      printf("mri_vol2vol spm-vbm done\n");
      printf("#VMPC# mri_vol2vol VmPeak %d\n",GetVmPeak());
      exit(0);
    }
    else if(istringnmatch(option, "--gcam",0) || istringnmatch(option, "--gcam0",0)) {
      LTA *srclta, *dstlta=NULL;
      if(nargc < 7){
	printf("  --gcam mov srclta gcam dstlta vsm interp out\n");
	argnerr(option,7);
      }
      printf("Loading mov %s\n",pargv[0]);
      mov = MRIread(pargv[0]);
      if(mov == NULL) exit(1);
      if(strcmp(pargv[1],"0")!=0){
	printf("Loading source LTA %s\n",pargv[1]);
	srclta = LTAread(pargv[1]);
	if(srclta == NULL) exit(1);
      } else srclta = NULL;
      if(strcmp(pargv[2],"0")!=0){
	printf("Loading GCAM %s\n",pargv[2]);
	gcam = GCAMread(pargv[2]);
	if(gcam == NULL) exit(1);
      } else gcam = NULL;
      if(strcmp(pargv[3],"0")!=0){
	printf("Loading destination LTA %s\n",pargv[3]);
	dstlta = LTAread(pargv[3]);
      }
      if(strcmp(pargv[4],"0")!=0){
	vsm = MRIread(pargv[4]);
	if(vsm == NULL) exit(1);
      } else vsm = NULL;
      sscanf(pargv[5],"%d",&interpcode);
      targvolfile = pargv[6];
      if(getenv("MY_MORPHS_DO_NOT_CONFORM_DEAL_WITH_IT") != NULL || istringnmatch(option, "--gcam0",0))
	out = MRIvol2volGCAM0(mov, srclta, gcam, dstlta, vsm, interpcode, NULL);
      else
	out = MRIvol2volGCAM(mov, srclta, gcam, dstlta, vsm, interpcode, NULL,pedir);
      if(out == NULL) exit(1);
      printf("Writing to %s\n",targvolfile);
      err = MRIwrite(out,targvolfile);
      if(err) exit(1);
      printf("mri_vol2vol gcam done\n");
      printf("#VMPC# mri_vol2vol VmPeak %d\n",GetVmPeak());
      exit(0);
      nargsused = 7;
    } 
    else if(istringnmatch(option, "--map-point",0) || istringnmatch(option, "--map-point-inv-lta",0)) {
      // --map-point a b c incoords lta outcoords outfile
      if(nargc < 7){
	printf("  --map-point a b c incoords lta outcoords outfile\n");
	printf("  coords: 1=tkras, 2=scannerras, 3=vox\n");
	argnerr(option,7);
      }
      double a,b,c;
      int incoords, outcoords;
      LTA *lta;
      sscanf(pargv[0],"%lf",&a);
      sscanf(pargv[1],"%lf",&b);
      sscanf(pargv[2],"%lf",&c);
      sscanf(pargv[3],"%d",&incoords);
      lta = LTAread(pargv[4]);
      if(!lta) exit(1);
      if(istringnmatch(option, "--map-point-inv-lta",0)) LTAinvert(lta,lta);
      sscanf(pargv[5],"%d",&outcoords);
      MATRIX *M = lta->get_matrix(incoords,outcoords,NULL);
      if(!M) exit(1);
      MATRIX *src = MatrixAlloc(4,1,MATRIX_REAL);
      src->rptr[1][1] = a;
      src->rptr[2][1] = b;
      src->rptr[3][1] = c;
      src->rptr[4][1] = 1;
      MATRIX *trg = MatrixMultiplyD(M,src,NULL);
      printf("%8.4f %8.4f %8.4f \n",trg->rptr[1][1],trg->rptr[2][1],trg->rptr[3][1]);
      if(strcmp(pargv[6],"nofile")!=0) {
	FILE *fp = fopen(pargv[6],"w");
	if(!fp){
	  printf("ERROR: could not open %s\n",pargv[6]);
	  exit(1);
	}
	printf("writing to %s\n",pargv[6]);
	fprintf(fp,"%8.4f %8.4f %8.4f \n",trg->rptr[1][1],trg->rptr[2][1],trg->rptr[3][1]);
	fclose(fp);
      }
      exit(0);
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
static void print_usage(void) {
printf("\n");
printf("mri_vol2vol\n");
printf("\n");
printf("  --mov  movvol       : input (or output template with --inv)\n");
printf("  --targ targvol      : output template (or input with --inv)\n");
printf("  --o    outvol       : output volume\n");
printf("  --disp dispvol      : displacement volume\n");
printf("  --downsample N1 N2 N3 : downsample factor (eg, 2) (do not include a targ or regsitration)\n");
printf("         sets --fill-average, --fill-upsample 2, and --regheader\n");
printf("\n");
printf("  --reg  register.dat : tkRAS-to-tkRAS matrix   (tkregister2 format)\n");
printf("  --lta  register.lta : Linear Transform Array (usually only 1 transform)\n");
printf("  --lta-inv  register.lta : LTA, invert (may not be the same as --lta --inv with --fstal)\n");
printf("  --fsl  register.fsl : fslRAS-to-fslRAS matrix (FSL format)\n");
printf("  --xfm  register.xfm : ScannerRAS-to-ScannerRAS matrix (MNI format)\n");
printf("  --regheader         : ScannerRAS-to-ScannerRAS matrix = identity\n");
printf("  --mni152reg         : target MNI152 space (need FSL installed)\n");
printf("  --s subject         : set matrix = identity and use subject for any templates\n");
printf("\n");
printf("  --inv               : sample from targ to mov\n");
printf("\n");
printf("  --tal               : map to a sub FOV of MNI305 (with --reg only)\n");
printf("  --talres resolution : set voxel size 1mm or 2mm (def is 1)\n");
printf("  --talxfm xfmfile    : default is talairach.xfm (looks in mri/transforms)\n");
printf("\n");
printf("  --m3z morph    : non-linear morph encoded in the m3z format\n");
printf("  --noDefM3zPath : flag indicating that the code should not be looking for \n");
printf("       the non-linear m3z morph in the default location (subj/mri/transforms), but should use \n");
printf("       the morph name as is\n");
printf("  --inv-morph    : compute and use the inverse of the m3z morph\n");
printf("\n");
printf("  --fstarg <vol>      : optionally use vol from subject in --reg as target. default is orig.mgz \n");
printf("  --crop scale        : crop and change voxel size\n");
printf("  --slice-crop sS sE  : crop output slices to be within sS and sE\n");
printf("  --slice-reverse     : reverse order of slices, update vox2ras\n");
printf("  --slice-bias alpha  : apply half-cosine bias field\n");
printf("\n");
printf("  --trilin            : trilinear interpolation (default)\n");
printf("  --nearest           : nearest neighbor interpolation\n");
printf("  --cubic             : cubic B-Spline interpolation\n");
printf("  --interp interptype : interpolation cubic, trilin, nearest (def is trilin)\n");
printf("  --fill-average      : compute mean of all source voxels in a given target voxel\n");
printf("  --fill-conserve     : compute sum  of all source voxels in a given target voxel\n");
printf("  --fill-upsample USF : source upsampling factor for --fill-{avg,cons} (default is 2)\n");
printf("\n");
printf("  --mul mulval   : multiply output by mulval\n");
printf("\n");
printf("  --vsm vsmvol <pedir> : Apply a voxel shift map. pedir: +/-1=+/-x, +/-2=+/-y, +/-3=+/-z (default +2)\n");
printf("  --vsm-pedir pedir : set pedir +/-1=+/-x, +/-2=+/-y, +/-3=+/-z (default +2)\n");
printf("\n");
printf("  --precision precisionid : output precision (def is float)\n");
printf("  --keep-precision  : set output precision to that of input\n");
printf("  --kernel            : save the trilinear interpolation kernel instead\n");
printf("   --copy-ctab : setenv FS_COPY_HEADER_CTAB to copy any ctab in the mov header\n");
printf("\n");
printf("  --gcam mov srclta gcam dstlta vsm interp out\n");
printf("     srclta, gcam, or vsm can be set to 0 to indicate identity (not regheader)\n");
printf("     if dstlta is 0, then uses gcam atlas geometry as output target\n");
printf("     direction is automatically determined from srclta and dstlta\n");
printf("     interp 0=nearest, 1=trilin, 5=cubicbspline\n");
printf("     vsm pedir can be set with --vsm-pedir\n");
printf("     DestVol -> dstLTA -> CVSVol -> gcam -> AnatVol -> srcLTA -> B0UnwarpedVol -> VSM -> MovVol (b0Warped)\n");
printf("\n");
printf("  --spm-warp mov movlta warp interp output\n");
printf("     mov is the input to be mapped \n");
printf("     movlta maps mov to the vbm input space (use 0 to ignore)\n");
printf("       if movlta=0, then input is anything that shares a RAS space with the VBM input\n");
printf("     warp is typically y_rinput.nii\n");
printf("     interp 0=nearest, 1=trilin\n");
printf("\n");
printf("  --map-point a b c incoords lta outcoords outfile : stand-alone option to map a point to another space\n");
printf("     coords: 1=tkras, 2=scannerras, 3=vox; outfile can be nofile\n");
printf("  --map-point-inv-lta a b c incoords lta outcoords outfile\n");
printf("      same as --map-point but inverts the lta\n");
printf("\n");
printf("\n");
printf("  --no-resample : do not resample, just change vox2ras matrix\n");
printf("  --no-resample-scale : do not resample, just change vox2ras matrix, using scale=voxsize\n");
printf("\n");
printf("  --rot   Ax Ay Az : rotation angles (deg) to apply to reg matrix\n");
printf("  --trans Tx Ty Tz : translation (mm) to apply to reg matrix\n");
printf("  --shear Sxy Sxz Syz : xz is in-plane\n");
printf("  --reg-final regfinal.dat : final reg after rot and trans (but not inv)\n");
printf("\n");
printf("  --synth : replace input with white gaussian noise\n");
printf("  --seed seed : seed for synth (def is to set from time of day)\n");
printf("\n");
printf("  --save-reg : write out output volume registration matrix\n");
printf("\n");
printf("  --help : go ahead, make my day\n");
printf("  --debug\n");
printf("  --version\n");
printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n%s\n\n",getVersion().c_str());
printf("\n");
printf("Resamples a volume into another field-of-view using various types of\n");
printf("matrices (FreeSurfer, FSL, SPM, and MNI). This is meant to be used\n");
printf("in conjunction with tkregister2.\n");
printf("\n");
printf("FLAGS AND ARGUMENTS\n");
printf("\n");
printf("--mov movvol\n");
printf("\n");
printf("This volume must have the same geometry as the --mov volume passed to\n");
printf("tkregister2 when creating/checking the registration file. By default,\n");
printf("this will be the input volume that will be resampled. If --inv is\n");
printf("specified, then this will become the geometry template for the output\n");
printf("instead.\n");
printf("\n");
printf("--targ targvol\n");
printf("\n");
printf("This volume must have the same geometry as the --targ volume passed to\n");
printf("tkregister2 when creating/checking the registration file. By default,\n");
printf("this will be the volume will be the geometry template for the output.\n");
printf("If --inv is specified, then this becomes the input volume that will be\n");
printf("resampled instead. The target volume can be implicitly specified with\n");
printf("--tal or --fstarg.\n");
printf("\n");
printf("--reg register.dat\n");
printf("\n");
printf("This simple text file contains the freesurfer registration matrix. It\n");
printf("is the same as the file passed to and generated by tkregister2 with\n");
printf("the --reg flag. If --tal or --fstarg is specified, then the subject\n");
printf("is obtained from the regfile.\n");
printf("\n");
printf("--fsl register.fsl\n");
printf("\n");
printf("Registration matrix created with the FSL flirt program using targ as\n");
printf("the reference and mov as input. Note: you cannot use any of the files\n");
printf("from $FSLDIR/etc/standard as mov or targ. These volumes do not have\n");
printf("geometry information in them, and FreeSurfer and FSL will default to\n");
printf("different things. Same as in tkregister2.\n");
printf("\n");
printf("--xfm register.xfm\n");
printf("\n");
printf("MNI-style registration matrix (eg, like one created with mritotal).\n");
printf("This matrix maps from mov Scanner-RAS to targ Scanner-RAS, where\n");
printf("'Scanner-RAS' is the vox2ras matrix as found in each file.\n");
printf("Same as in tkregister2.\n");
printf("\n");
printf("--regheader\n");
printf("\n");
printf("Create a registration matrix assuuming that the mov Scanner-RAS and\n");
printf("targ Scanner-RAS are the same. This is the same as using a register.xfm\n");
printf("with the identity matrix in it. This can be used with some SPM\n");
printf("registrations (which change only the matrix in the .mat file).\n");
printf("Same as in tkregister2.\n");
printf("\n");
printf("--mni152reg \n");
printf("\n");
printf("Target MNI152 space. If the mov volume is in the native space of an\n");
printf("individual, then also supply a registration (--reg). This registration\n");
printf("is concatenated with that in subject/mri/transforms/reg.mni152.2mm.dat\n");
printf("(created with mni152reg) to produce a registration from the mov vol\n");
printf("to MNI152 (defined by $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz).\n");
printf("If the data are in fsaverage 2mm space, then do not supply a --reg.\n");
printf("Instead, $FREESURFER_HOME/average/mni152.register.dat is used\n");
printf("as the registration. Do not supply a target volume with --mni152reg.\n");
printf("\n");
printf("--inv\n");
printf("\n");
printf("Invert the transform. The movvol becomes the geometry template for the\n");
printf("output, and the targvol becomes the input that will be resampled.\n");
printf("\n");
printf("--o outvol\n");
printf("\n");
printf("Output volume.  By default, this will be the movvol resmapled into the\n");
printf("targvol space (and so will have the same geometry as the targvol). If\n");
printf("--inv is specified, then this will be the targvol resmapled into the\n");
printf("movvol space (and so will have the same geometry as the movvol). By\n");
printf("default, the output volume will be float, but this can be changed\n");
printf("with --precision. By default, the interpolation will be done with\n");
printf("trilinear, but this can be changed with --interp.\n");
printf("\n");
printf("  --keep-precision  : set output precision to that of input\n");
printf("\n");
printf("--tal\n");
printf("\n");
printf("Resample the movvol to talairach (ie, MNI305) space. The talairach\n");
printf("matrix is obtained from talairach.xfm from\n");
printf("SUBJECTS_DIR/subjid/transforms. subjid is read from the register.dat\n");
printf("file. Requires --reg. Do not specify --targ as the target volume is\n");
printf("implicitly set to $FREESURFER_HOME/average/mni305.cor.subfovV.mgz,\n");
printf("where V is either 1 (for 1mm) or 2 (for 2mm). 2mm is used by default,\n");
printf("but this can be changed with --talres.  mni305.cor.subfovV.mgz the\n");
printf("MNI305 (1mm or 2mm isotropic) volume in a reduced FOV that covers only\n");
printf("the brain. Reducing the FOV saves space relative to the 256^3 COR FOV.\n");
printf("The transformation matrix is computed as R*inv(Xtal)*inv(Rtal), where\n");
printf("Xtal is talairach.xfm matrix, R is the matrix in the regfile, and Rtal\n");
printf("maps from the talairach COR FOV to the SubFOV (mni305.cor.subfovV.reg).\n");
printf("If you want to sample the targvol from talairach space into the movvol\n");
printf("space, then specify --inv. SUBJECTS_DIR is read from the environment\n");
printf("or can be specified with --sd.\n");
printf("\n");
printf("--fstalres resmm\n");
printf("\n");
printf("Set the resolution of the output when using --fstal. By default, it\n");
printf("is 2 mm, but can be changed to 1.0 mm with --fstalres 1\n");
printf("\n");
printf("--fstarg <vol>\n");
printf("\n");
printf("Set target to vol from the subject found in register.dat\n");
printf("file. If vol is not specified, uses orig.mgz. Requires --reg.  \n");
printf("Same as tkregister2.\n");
printf("\n");
printf("--crop scale\n");
printf("\n");
printf("Crop mov volume down to minimum size to fit non-zero voxels. The size of\n");
printf("the voxels is reduced by scale (ie, --crop 2 would crop and reduce the\n");
printf("voxel size by a factor of 2, eg 1.0 mm becomes 0.5 mm).\n");
printf("\n");
printf("--slice-crop start end\n");
printf("\n");
printf("Crop output volume to be within slices start and end. The geometry is \n");
printf("updated to reflect the new limits.\n");
printf("\n");
printf("--interp method\n");
printf("\n");
printf("Interpolate the output based on the given method. Legal values are:\n");
printf("cubic, trilin and nearest. trilin is the default. Can also use\n");
printf("--cubic, --trilin or --nearest.\n");
printf("\n");
printf("--precision precisionid\n");
printf("\n");
printf("Set output precision to precisionid. Legal values are uchar, short,\n");
printf("int, long, and float. Default is float.\n");
printf("\n");
printf("--kernel\n");
printf("\n");
printf("Save the trilinear interpolation kernel at each voxel instead of the\n");
printf("interpolated image.\n");
printf("\n");
printf("--nomr\n");
printf("\n");
printf("Don't copy the template MR parameters, but instead preserve the input volume \n");
printf("ones\n");
printf("\n");
printf("--help\n");
printf("\n");
printf("Prints out all this information.\n");
printf("\n");
printf("--gdiagno diagnostic level\n");
printf("\n");
printf("Sets the diagnostic level (only good for debuggin').\n");
printf("\n");
printf("--version\n");
printf("\n");
printf("Print out version string and exit.\n");
printf("\n");
printf("\n");
printf("EXAMPLES:\n");
printf("\n");
printf("Below are some exampls of how one might use mri_vol2vol. They are not\n");
printf("exhaustive of all the possible combinations of options. Typically, one\n");
printf("uses a template to establish the registration, then resamples data\n");
printf("that are in correspondence with the template.\n");
printf("\n");
printf("1. If a functional volume is f.bhdr (or f.nii.gz, or f.mgh, etc), and the\n");
printf("subject is bert, and the registration file is register.dat, then\n");
printf("running the following command should show that they are in\n");
printf("registration:\n");
printf("\n");
printf("tkregister2 --reg register.dat --mov f.nii.gz\n");
printf("\n");
printf("If they are not, then fix it because nothing below is going to work. You\n");
printf("can also check the registration with:\n");
printf("\n");
printf("tkmedit bert orig.mgz  -overlay f.nii.gz -overlay-reg register.dat\n");
printf("\n");
printf("The register.dat will look something like this\n");
printf("----------------- register.dat --------------------------\n");
printf("bert\n");
printf("3.125\n");
printf("5.000\n");
printf("0.150000\n");
printf("1.000000e+00 0.000000e+00 0.000000e+00 -2.252487e+00\n");
printf("0.000000e+00 -8.902127e-01 4.555448e-01 2.342102e+00\n");
printf("0.000000e+00 4.555449e-01 8.902128e-01 -2.159538e-01\n");
printf("0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00\n");
printf("round\n");
printf("----------------- register.dat --------------------------\n");
printf("\n");
printf("1.A. To resample the functional into anatomical space:\n");
printf("\n");
printf("mri_vol2vol --reg register.dat --mov f.nii.gz --fstarg \\\n");
printf("   --o f-in-anat.mgh\n");
printf("\n");
printf("f-in-anat.mgh will have the same size and geometry as\n");
printf("bert/mri/orig.mgz.  You can test the result in two ways:\n");
printf("\n");
printf("  # This will show the low-res functional alighned with its resampled self\n");
printf("  tkregister2 --reg register.dat --mov f.nii.gz --targ f-in-anat.mgh\n");
printf("\n");
printf("  # This will show the resampled functional aligned with the anatomical\n");
printf("  tkregister2 --mov f-in-anat.mgh --targ $SUBJECTS_DIR/bert/mri/orig.mgz \\\n");
printf("     --regheader --reg /tmp/reg\n");
printf("\n");
printf("1.B. To resample the anatomical into the functional space. This is\n");
printf("basically the same command line as 1.A, but --inv has been included\n");
printf("and the name of the output is changed.\n");
printf("\n");
printf("mri_vol2vol --reg register.dat --mov f.nii.gz --fstarg \\\n");
printf("  --o anat-in-func.mgh --inv\n");
printf("\n");
printf("anat-in-func.mgh will be the same size and geometry as f.nii.gz.\n");
printf("You can test the result in two ways:\n");
printf("\n");
printf("  # This will show the low-res anat aligned with its hires self\n");
printf("  tkregister2 --reg register.dat --mov anat-in-func.mgh\n");
printf("\n");
printf("  # This will show the resampled anat aligned with the functional\n");
printf("  tkregister2 --mov anat-in-func.mgh --targ f.nii.gz \\\n");
printf("     --regheader --reg /tmp/reg\n");
printf("\n");
printf("1.C Map functional to anatomical without resampling. Rather, change\n");
printf("the vox2ras (sform/qform) matrix. This is the same cmd line as 1.A,\n");
printf("but --no-resample as been added.\n");
printf("\n");
printf("mri_vol2vol --reg register.dat --mov f.nii.gz --fstarg \\\n");
printf("   --o f.new.vox2ras.nii.gz --no-resample\n");
printf("\n");
printf("f.new.vox2ras.nii.gz will have the same dimension and voxel size\n");
printf("as f.nii.gz, but its vox2ras (sform/qform) matrix will have changed.\n");
printf("You can check the registration in two ways:\n");
printf("\n");
printf("  # The registration is created implicitly from the vox2ras matrix\n");
printf("  # (that is what --regheader does). There's no need to specify\n");
printf("  # and input registration\n");
printf("  tkregister2 --mov f.new.vox2ras.nii.gz --s bert --regheader --reg /tmp/reg\n");
printf("\n");
printf("  # Display the functional as an overlay in tkmedit (no registration\n");
printf("  # needed).\n");
printf("  tkmedit bert orig.mgz -overlay f.new.vox2ras.nii.gz\n");
printf("\n");
printf("1.D Map a binary mask in functional space to anatomical space. This is\n");
printf("basically the same cmd line as 1.A, but --interp nearest has been\n");
printf("added so that it does not try to interpolate the mask (ie, it will\n");
printf("still be binary after resampling):\n");
printf("\n");
printf("mri_vol2vol --reg register.dat --mov mask.nii.gz --fstarg \\\n");
printf("   --o mask-in-anat.mgh --interp nearest\n");
printf("\n");
printf("2. Map functional to/from talairach (MNI305) space. This uses a\n");
printf("two-stage registration: func-to-anat (register.dat) and\n");
printf("anat-to-talairach (talairach.xfm).\n");
printf("\n");
printf("Make sure that sure the func-to-anat reg is correct as was done in\n");
printf("Example 1. Next, make sure that the anat-to-tal is correct with:\n");
printf("\n");
printf("tkregister2 --s bert --fstal\n");
printf("\n");
printf("2.A Map functional to talairach (MNI305) space with 2mm isotropic\n");
printf("resolution. This is very similar to 1.A with the addition of --tal\n");
printf("and --talres 2.\n");
printf("\n");
printf("mri_vol2vol --mov f.nii.gz --reg register.dat \\\n");
printf("     --o f-in-tal.2mm.mgh --tal --talres 2\n");
printf("\n");
printf("f-in-tal.2mm.mgh will be 2mm isotropic with the same geometry as\n");
printf("$FREESURFER_HOME/average/mni305.cor.subfov2.mgz. This command will\n");
printf("also create f-in-tal.2mm.mgh.reg, which will register the volume with\n");
printf("any average MNI305 FreeSurfer subject (fsaverage is used by default).\n");
printf("The resampling can be checked with:\n");
printf("\n");
printf("  # This will show the functional with the fsaverage anatomical\n");
printf("  tkregister2 --mov f-in-tal.2mm.mgh --reg f-in-tal.2mm.mgh.reg\n");
printf("\n");
printf("2.B Map functional to talairach (MNI305) space with 1mm isotropic\n");
printf("resolution. Same as 2.A but use --talres 1.\n");
printf("\n");
printf("mri_vol2vol --mov f.nii.gz --reg register.dat \\\n");
printf("     --o f-in-tal.1mm.mgh --tal --talres 1\n");
printf("\n");
printf("f-in-tal.1mm.mgh will take up 8 times as much space as f-in-tal.2mm.mgh\n");
printf("\n");
printf("3. Apply an MNI transform to data by resampling the anatomical orig.mgz\n");
printf("into talairach space using bert/mri/transforms/talairach.xfm:\n");
printf("\n");
printf("First, check that the talairach.xfm is correct (this is basically the same\n");
printf("thing as 'tkregister2 --s bert --fstal' in Example 2:\n");
printf("\n");
printf(" cd bert/mri\n");
printf(" tkregister2 --targ orig.mgz \\\n");
printf("     --mov $FREESURFER_HOME/average/mni305.cor.mgz \\\n");
printf("     --xfm transforms/talairach.xfm --reg /tmp/reg\n");
printf("\n");
printf(" Now resample\n");
printf(" mri_vol2vol --mov orig.mgz \\\n");
printf("     --targ $FREESURFER_HOME/average/mni305.cor.mgz \\\n");
printf("     --xfm transforms/talairach.xfm  \\\n");
printf("     --o orig-in-mni305.mgz\n");
printf("\n");
printf(" Now test the resampling:\n");
printf(" tkregister2 --mov orig-in-mni305.mgz \\\n");
printf("    --targ $FREESURFER_HOME/average/mni305.cor.mgz \\\n");
printf("    --reg /tmp/reg --regheader\n");
printf("\n");
printf("\n");
printf("FORMATS\n");
printf("\n");
printf("Data file format can be specified implicitly (through the path name)\n");
printf("or explicitly. All formats accepted by mri_convert can be used.\n");
printf("\n");
printf("BUGS\n");
printf("\n");
printf("sinc interpolation is broken except for maybe COR to COR.\n");
printf("\n");
printf("\n");
printf("BUG REPORTING\n");
printf("\n");
printf("Report bugs to analysis-bugs@nmr.mgh.harvard.edu. Include the following\n");
printf("formatted as a list as follows: (1) command-line, (2) directory where\n");
printf("the program was run (for those in the MGH-NMR Center), (3) version,\n");
printf("(4) text output, (5) description of the problem.\n");
printf("\n");
printf("SEE ALSO\n");
printf("\n");
printf("mri_convert, tkregister2\n");
printf("\n");
printf("\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR==NULL) {
    printf("ERROR: SUBJECTS_DIR undefined.\n");
    exit(1);
  }
  if(movvolfile == NULL && ( lta == NULL || ! invert) ) {
    printf("ERROR: No mov volume supplied.\n");
    exit(1);
  }
  if(fstarg && targvolfile != NULL) {
    printf("ERROR: Do not specify a targ volume with --fstarg.\n");
    exit(1);
  }
  if(fstal && targvolfile != NULL) {
    printf("ERROR: Do not specify a targ volume with --tal.\n");
    exit(1);
  }
  if(fstal && fstarg) {
    printf("ERROR: cannot specify a --tal and --fstarg.\n");
    exit(1);
  }

  if(fstarg) {
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,fstargfile);
    if (!fio_FileExistsReadable(tmpstr))
      sprintf(tmpstr,"%s/%s/mri/orig",SUBJECTS_DIR,subject);
    targvolfile = strcpyalloc(tmpstr);
    printf("Using %s as targ volume\n",targvolfile);
  }

  if(fstal){
    if (fstalres > 0) {
      sprintf(tmpstr,"%s/average/mni305.cor.subfov%d.mgz",FSH,fstalres);
      targvolfile = strcpyalloc(tmpstr);
    } else {
      sprintf(tmpstr,"%s/average/mni305.cor.mgz",FSH);
      targvolfile = strcpyalloc(tmpstr);
    }
  }

  if (targvolfile != NULL && DownSample[0]){
    printf("ERROR: cannot spec target file and downsample.\n");
    exit(1);
  }

  if (targvolfile == NULL && DownSample[0]==0){
    printf("ERROR: No target volume supplied.\n");
    exit(1);
  }

  if(outvolfile == NULL && DispFile == NULL) {
    printf("ERROR: No output volume supplied.\n");
    exit(1);
  }
  if(outvolfile){
    outdir = fio_dirname(outvolfile);
    err = mkdir(outdir,0777);
    if (err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n",outdir);
      exit(1);
    }
  }
  if(fstal && regheader){
    printf("ERROR: cannot use --tal and --regheader\n");
    exit(1);
  }
  if(keepprecision){
    mov = MRIreadHeader(movvolfile,MRI_VOLUME_TYPE_UNKNOWN);
    if(mov==NULL) exit(1);
    precisioncode = mov->type;
    precision = MRIprecisionString(precisioncode);
    MRIfree(&mov);
  }
  if(TargMNI152){
    if(regfile == NULL){
      sprintf(tmpstr,"%s/average/mni152.register.dat",getenv("FREESURFER_HOME"));
      regfile = strcpyalloc(tmpstr);
      err = regio_read_register(regfile, &subject, &ipr, &bpr,
                                &intensity, &R, &float2int);
      if (err) exit(1);
    } else {
      MATRIX *Q, *invQ, *M;
      sprintf(tmpstr,"%s/%s/mri/transforms/reg.mni152.2mm.dat",SUBJECTS_DIR,subject);
      Q = regio_read_registermat(tmpstr);
      if(Q == NULL){
	printf("Run mni152reg --s %s to create reg.mni152.2mm.dat\n",subject);
	exit(1);
      }
      invQ = MatrixInverse(Q,NULL);
      M = MatrixMultiply(R,invQ,NULL);
      MatrixFree(&R);
      R = MatrixCopy(M,NULL);
      MatrixFree(&Q);
      MatrixFree(&invQ);
      MatrixFree(&M);
      printf("New MNI152 R\n");
      MatrixPrint(stdout,R);
    }
  }
  
  if(lta != NULL && !fstal){
    MRI *mrimovtmp,*mritrgtmp;
    printf("%s %s\n",movvolfile,targvolfile);
    mrimovtmp = MRIreadHeader(movvolfile,MRI_VOLUME_TYPE_UNKNOWN);
    if(mrimovtmp == NULL) exit(1);
    mritrgtmp = MRIreadHeader(targvolfile,MRI_VOLUME_TYPE_UNKNOWN);
    if(mritrgtmp == NULL) exit(1);
    lta = LTAchangeType(lta,LINEAR_RAS_TO_RAS);
    LTAmodifySrcDstGeom(lta, mrimovtmp, mritrgtmp);
    R = TransformLTA2RegDat(lta);
    ipr = lta->xforms[0].src.xsize ;
    bpr = lta->xforms[0].src.zsize ;
    MRIfree(&mrimovtmp);
    MRIfree(&mritrgtmp);
  }

  if(!fstal && !DoCrop && !fstarg && targvolfile == NULL &&  ( lta == NULL || invert) && !DownSample[0]) {
    printf("ERROR: No targ volume supplied.\n");
    exit(1);
  }
  if(DoCrop && targvolfile != NULL) {
    printf("ERROR: Do not specify a targ volume with --crop.\n");
    exit(1);
  }
  if(xfmfile != NULL && regfile != NULL) {
    printf("ERROR: cannot specify both --xfm and --reg.\n");
    exit(1);
  }
  if (regheader && regfile != NULL) {
    printf("ERROR: cannot specify both --regheader and --reg.\n");
    exit(1);
  }
  if(fstarg && regfile == NULL && subject == NULL) {
    printf("ERROR: Need --reg with --fstarg.\n");
    exit(1);
  }

  interpcode = MRIinterpCode(interpmethod);
  if (interpcode < 0) {
    printf("ERROR: interpolation method %s unrecognized\n",interpmethod);
    printf("       legal values are nearest, trilin, and sinc\n");
    exit(1);
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"movvol %s\n",movvolfile);
  if (targvolfile)
    fprintf(fp,"targvol %s\n",targvolfile);
  fprintf(fp,"outvol %s\n",outvolfile);
  if (regfile) fprintf(fp,"regfile %s\n",regfile);
  if (xfmfile) fprintf(fp,"xfmfile %s\n",xfmfile);
  fprintf(fp,"invert %d\n",invert);
  fprintf(fp,"tal    %d\n",fstal);
  fprintf(fp,"talres %d\n",fstalres);
  fprintf(fp,"regheader %d\n",regheader);
  fprintf(fp,"noresample %d\n",noresample);
  fprintf(fp,"interp  %s (%d)\n",interpmethod,interpcode);
  if (interpcode == SAMPLE_SINC) fprintf(fp,"sinc hw  %d\n",sinchw);
  fprintf(fp,"precision  %s (%d)\n",precision,precisioncode);
  fprintf(fp,"Gdiag_no  %d\n",Gdiag_no);

  if(DoMorph){
    fprintf(fp,"Morphing\n");
    fprintf(fp,"InvertMorph %d\n",InvertMorph);
  }

  fprintf(fp,"Synth      %d\n",synth);
  fprintf(fp,"SynthSeed  %d\n",SynthSeed);

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
  if (nargc <= nth) return(0);

  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) return(0);

  return(1);
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
static MATRIX *LoadRtal(int talres) {
  char *FSH;
  char rtalfile[2000];
  float ipr, bpr, intensity;
  int float2int, err;
  MATRIX *Rtal;

  FSH = getenv("FREESURFER_HOME");
  if (FSH==NULL) {
    printf("ERROR: FREESURFER_HOME undefined.\n");
    exit(1);
  }
  sprintf(rtalfile,"%s/average/mni305.cor.subfov%d.reg",FSH,talres);
  err = regio_read_register(rtalfile, &subject, &ipr, &bpr,
                            &intensity, &Rtal, &float2int);
  if (err) exit(1);
  return(Rtal);
}
/*-----------------------------------------------------*/
MATRIX *LoadRfsl(char *fname) {
  MATRIX *FSLRegMat;
  FILE *fp;
  int i,j,n;

  fp = fopen(fname,"r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n",fname);
    exit(1);
  }
  FSLRegMat = MatrixAlloc(4,4,MATRIX_REAL);
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      n = fscanf(fp,"%f",&(FSLRegMat->rptr[i+1][j+1]));
      if (n != 1) {
        printf("ERROR: reading %s, row %d, col %d\n",fname,i,j);
        return(NULL);
      }
    }
  }
  return(FSLRegMat);
}
/*
  \fn MRI *MRIvol2volGCAM(MRI *src, LTA *srclta, GCA_MORPH *gcam, LTA *dstlta, MRI *vsm, int sample_type, MRI *dst)
  \brief Converts one volume into another using as many as four transforms: VSM, src linear, gcam/m3z, dst linear.
  Any may be NULL except dstlta in which case they are assumed to be the identity. This function allows for
  transforming from the functional space to a sub FoV of CVS space including B0 distortion correction.
  This is kind of confusing because the transforms all go backwards
  DestVol --> dstLTA --> CVSVol --> gcam --> AnatVol --> srcLTA --> B0UnwarpedVol --> VSM --> MovVol (b0Warped)
  DestVol --> dstLTA --> srcLTA --> B0UnwarpedVol --> VSM --> MovVol (b0Warped)
  DestVol --> dstLTA --> B0UnwarpedVol --> VSM --> MovVol (b0Warped)
  Three more possibilites with removing VSM
 */
MRI *MRIvol2volGCAM(MRI *src, LTA *srclta, GCA_MORPH *gcam, LTA *dstlta, MRI *vsm, int sample_type, MRI *dst, int pedir)
{
  int c,r,s,f,out_of_gcam,cvsm,rvsm,iss;
  //VOL_GEOM *vgdst_src,*vgdst_dst;
  MATRIX *crsDst, *crsGCAM=NULL, *crsAnat=NULL, *crsSrc=NULL, *Vdst, *Vsrc;
  double val,v;
  MRI_BSPLINE * bspline = NULL;
  float dvsm, *valvect;
  Timer timer;

  printf("MRIvol2volGCAM(): ---------+++++++++++++++++++----------------------\n");

  if(!vsm) printf("MRIvol2volGCAM(): VSM not used\n");
  if(!gcam) printf("MRIvol2volGCAM(): GCAM not used\n");
  if(!srclta) printf("MRIvol2volGCAM(): Source LTA not used\n");
  printf("MRIvol2volGCAM(): interpolation type is %d %s\n",sample_type,MRIinterpString(sample_type));

  if(vsm){
    if(MRIdimMismatch(src,vsm,0)){
      printf("ERROR: MRIvol2volGCAM(): src-vsm dim mismatch\n");
      return(NULL);
    }
  }
  int free_dstlta = 0;
  if(!dstlta){
    if(gcam==NULL){
      printf("ERROR: MRIvol2volGCAM(): both dstlta and gcam are NULL\n");
      return(NULL);
    }
    printf("Dst LTA is null, so just using gcam atas as output target\n");
    dstlta = TransformRegDat2LTA(&gcam->atlas, &gcam->atlas, NULL); // nothing to do with reg.dat
    free_dstlta = 1;
  }

  // Make sure that the source lta points in the right direction
  LTA *srcltacopy = srclta;
  if(srclta){
    srcltacopy = LTAcopy(srclta,NULL); // don't modify the passed lta
    if(LTAmriIsSource(srcltacopy, src)){
      printf("MRIvol2volGCAM(): Inverting Source LTA\n");
      LTAinvert(srclta, srcltacopy); // have to actually invert it for below
    }
    else if(!LTAmriIsTarget(srcltacopy, src)){
      printf("ERROR: MRIvol2volGCAM(): neither src nor dst of Source LTA matches the mov vol geom\n");
      printf("mov vol geom ======================\n");
      src->vgprint();
      printf("Source LTA ======================\n");
      LTAprint(stdout, srcltacopy);
      printf("ERROR: MRIvol2volGCAM(): neither src nor dst of Source LTA matches the mov vol geom\n");
      return(NULL);
    }
    if(gcam && !vg_isEqual(&srcltacopy->xforms[0].src,&(gcam->image))){
      // make sure that the gcam->image matches the src of the Source LTA
      printf("ERROR: MRIvol2volGCAM(): gcam image vol geom does not match Source LTA src vol geom\n");
      printf("gcam image vol geom ======================\n");
      gcam->image.vgprint();
      printf("Source LTA ======================\n");
      LTAprint(stdout, srcltacopy);
      printf("ERROR: MRIvol2volGCAM(): gcam image vol geom does not match Source LTA src vol geom\n");
      return(NULL);
    }
    // Finally, make sure that it is vox2vox
    if(srcltacopy->type != LINEAR_VOX_TO_VOX){
      printf("MRIvol2volGCAM(): Chaning Source LTA type to vox2vox\n");
      LTAchangeType(srcltacopy, LINEAR_VOX_TO_VOX) ;
    }
    Vsrc = MatrixCopy(srcltacopy->xforms[0].m_L,NULL);
  }
  else Vsrc = MatrixIdentity(4,NULL);

  // Make sure the dest lta matches geom and points in the right direction
  int InvertDstLTA=0;
  if(gcam){
    // given CRS in the "atlas" space, the GCAM returns the CRS in the image space
    // gcam->{image,atlas}
    //VOL_GEOM   image;   /* destination/target/output of the transform  */
    //VOL_GEOM   atlas ;  /* source/move/input of the transform       */
    if(vg_isEqual(&dstlta->xforms[0].dst,&(gcam->atlas))) InvertDstLTA=0;
    else if(vg_isEqual(&dstlta->xforms[0].src,&(gcam->atlas))) InvertDstLTA=1;
    else {
      printf("ERROR: MRIvol2volGCAM(): neither src nor dst of Dest LTA matches the atlas vol geom in the GCAM\n");
      printf("gcam atlas vol geom ======================\n");
      gcam->atlas.vgprint();
      printf("Dest LTA ======================\n");
      LTAprint(stdout, dstlta);
      printf("ERROR: MRIvol2volGCAM(): neither src nor dst of Dest LTA matches the atlas vol geom in the GCAM\n");
      return(NULL);
    }
  }
  else if(srcltacopy){ // gcam not specified, so check against the srclta; srclta must have been inverted if needed
    if(vg_isEqual(&dstlta->xforms[0].dst,&(srcltacopy->xforms[0].src))) InvertDstLTA=0;
    if(vg_isEqual(&dstlta->xforms[0].src,&(srcltacopy->xforms[0].src))) InvertDstLTA=1;
    else {
      printf("ERROR: MRIvol2volGCAM(): neither src nor dst of Dest LTA matches the src vol geom of the Source LTA\n");
      printf("Source LTA src vol geom ======================\n");
      srcltacopy->xforms[0].src.vgprint();
      printf("Dest LTA ======================\n");
      LTAprint(stdout, dstlta);
      printf("ERROR: MRIvol2volGCAM(): neither src nor dst of Dest LTA matches the src vol geom of the Source LTA\n");
      return(NULL);
    }
  }
  else { // must match mov/src volume
    if(LTAmriIsTarget(dstlta, src))      InvertDstLTA = 0;
    else if(LTAmriIsSource(dstlta, src)) InvertDstLTA = 1;
    else {
      printf("ERROR: MRIvol2volGCAM(): neither src nor dst of Dest LTA matches the mov/src vol geom\n");
      printf("mov vol geom ======================\n");
      src->vgprint();
      printf("Dest LTA ======================\n");
      LTAprint(stdout, dstlta);
      printf("ERROR: MRIvol2volGCAM(): neither src nor dst of Dest LTA matches the mov/src vol geom\n");
      return(NULL);
    }
  }

  LTA *dstltacopy;
  if(InvertDstLTA){
    printf("MRIvol2volGCAM(): Inverting Dest LTA\n");
    dstltacopy = LTAinvert(dstlta,NULL);
  }
  else dstltacopy = LTAcopy(dstlta,NULL); // don't modify the passed lta
  if(dstltacopy->type != LINEAR_VOX_TO_VOX){
    printf("MRIvol2volGCAM(): Chaning Destination LTA type to vox2vox\n");
    LTAchangeType(dstltacopy, LINEAR_VOX_TO_VOX) ;
  }
  Vdst = MatrixCopy(dstltacopy->xforms[0].m_L,NULL);
  // Vdst = MatrixInverse(dstlta->xforms[0].m_L,NULL); // was in original

  if(dst == NULL){
    dst = MRIallocFromVolGeom(&(dstltacopy->xforms[0].src), MRI_FLOAT, src->nframes, 0);
    if(dst==NULL) return(NULL);
    MRIcopyPulseParameters(src,dst);
  }
  // else should check

  if(sample_type == SAMPLE_CUBIC_BSPLINE) bspline = MRItoBSpline(src,NULL,3);
  valvect = (float *) calloc(sizeof(float),src->nframes);

  crsDst = MatrixAlloc(4,1,MATRIX_REAL);
  crsDst->rptr[4][1] = 1;
  crsAnat = MatrixAlloc(4,1,MATRIX_REAL);
  crsAnat->rptr[4][1] = 1;
  // scroll thru the CRS in the output/dest volume
  timer.reset();
  for(c=0; c < dst->width; c++){
    for(r=0; r < dst->height; r++){
      for(s=0; s < dst->depth; s++){
	// CRS in destination volume
	crsDst->rptr[1][1] = c;
	crsDst->rptr[2][1] = r;
	crsDst->rptr[3][1] = s;

	// Compute the CRS in the GCAM "atlas" CRS space
	crsGCAM = MatrixMultiplyD(Vdst,crsDst,crsGCAM);

	if(gcam){
	  // Compute the CRS in the anatomical "image" space
	  out_of_gcam = GCAMsampleMorph(gcam, 
					crsGCAM->rptr[1][1],crsGCAM->rptr[2][1],crsGCAM->rptr[3][1],
					&crsAnat->rptr[1][1],&crsAnat->rptr[2][1],&crsAnat->rptr[3][1]);
	  if(out_of_gcam) continue;
	}
	else crsAnat = MatrixCopy(crsGCAM,crsAnat);

	// Compute the CRS in the Source Space (eg, functional B0warped space )
	crsSrc = MatrixMultiply(Vsrc,crsAnat,crsSrc);

        if(vsm){
          /* crsSrc is the CRS in the undistored source space. This
	     code computes the crsSrc in the space distorted by B0
	     inhomogeneity (ie, the space of MRI *src).  This is just
	     a change in the row value as given by the voxel shift map
	     (VSM). The VSM must have the same dimensions as src. */
	  cvsm = floor(crsSrc->rptr[1][1]);
	  rvsm = floor(crsSrc->rptr[2][1]);
	  iss = nint(crsSrc->rptr[3][1]);

	  if(cvsm < 0 || cvsm+1 >= src->width)  continue;
	  if(rvsm < 0 || rvsm+1 >= src->height) continue;
	  if(iss < 0  || iss+1  >= src->depth) continue;
	  // Dont sample outside the BO mask indicated by vsm=0
	  v = MRIgetVoxVal(vsm,cvsm,rvsm,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm+1,rvsm,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm,rvsm+1,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm+1,rvsm+1,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  /* Performs 3D interpolation. May want to use iss instead of crsSrc->rptr[3][1]
	     to make it a 2D interpolation. Not sure.*/
          MRIsampleSeqVolume(vsm, crsSrc->rptr[1][1],crsSrc->rptr[2][1],crsSrc->rptr[3][1], &dvsm, 0, 0);
	  if(dvsm == 0) continue;
	  if(pedir<0) dvsm *= dvsm;
          crsSrc->rptr[pedir][1] += dvsm;
        }

	// Check for out of the source FoV
	if(crsSrc->rptr[1][1] < 0 || crsSrc->rptr[1][1] >= src->width)  continue;
	if(crsSrc->rptr[2][1] < 0 || crsSrc->rptr[2][1] >= src->height) continue;
	if(crsSrc->rptr[3][1] < 0 || crsSrc->rptr[3][1] >= src->depth)  continue;

        if(sample_type != SAMPLE_TRILINEAR)
	  for(f=0; f < src->nframes; f++){
	    if(sample_type == SAMPLE_CUBIC_BSPLINE)
	      MRIsampleBSpline(bspline, crsSrc->rptr[1][1],crsSrc->rptr[2][1],crsSrc->rptr[3][1], f, &val);
	    else
	      MRIsampleVolumeFrameType(src,
				       crsSrc->rptr[1][1],crsSrc->rptr[2][1],crsSrc->rptr[3][1],
				       f, sample_type, &val) ;
	    MRIsetVoxVal(dst,c,r,s,f, val);
	  }
	else {
	  // This will do the same as above, it is just faster with multiple frames
          MRIsampleSeqVolume(src, crsSrc->rptr[1][1],crsSrc->rptr[2][1],crsSrc->rptr[3][1],
			     valvect,0, src->nframes-1) ;
	  for(f=0; f < src->nframes; f++)  MRIsetVoxVal(dst,c,r,s,f, valvect[f]);
	}

      } // s
    } // r
  } // c

  MatrixFree(&crsDst);
  MatrixFree(&crsGCAM);
  MatrixFree(&crsAnat);
  MatrixFree(&crsSrc);
  MatrixFree(&Vdst);
  MatrixFree(&Vsrc);
  if(free_dstlta) LTAfree(&dstlta);
  if(bspline) MRIfreeBSpline(&bspline);
  free(valvect);
  if(srcltacopy) LTAfree(&srcltacopy);
  if(dstltacopy) LTAfree(&dstltacopy);

  printf("MRIvol2volGCAM: t=%6.4f\n", timer.seconds());
  fflush(stdout);

  return(dst);
}


// This is the "old" version without all the checks and flexibility of
// the new one. This can be called from the command line because the
// new one has a lot of changes and I'm afraid I might have broken
// something.
MRI *MRIvol2volGCAM0(MRI *src, LTA *srclta, GCA_MORPH *gcam, LTA *dstlta, MRI *vsm, int sample_type, MRI *dst)
{
  int c,r,s,f,out_of_gcam,cvsm,rvsm,iss;
  VOL_GEOM *vgdst_src,*vgdst_dst;
  MATRIX *crsDst, *crsGCAM=NULL, *crsAnat=NULL, *crsSrc=NULL, *Vdst, *Vsrc;
  double val,v;
  MRI_BSPLINE * bspline = NULL;
  float drvsm, *valvect;
  Timer timer;

  printf("MRIvol2volGCAM0(): ===========================================================\n");
  if(getenv("MY_MORPHS_DO_NOT_CONFORM_DEAL_WITH_IT") != NULL) printf("MY_MORPHS_DO_NOT_CONFORM_DEAL_WITH_IT is set\n");
  else printf("MY_MORPHS_DO_NOT_CONFORM_DEAL_WITH_IT is NOT set\n");

  if(!vsm) printf("MRIvol2volGCAM(): VSM not used\n");
  if(!gcam) printf("MRIvol2volGCAM(): GCAM not used\n");
  if(!srclta) printf("MRIvol2volGCAM(): Source LTA not used\n");
  printf("MRIvol2volGCAM(): interpolation type is %d %s\n",sample_type,MRIinterpString(sample_type));

  if(vsm){
    if(MRIdimMismatch(src,vsm,0)){
      printf("ERROR: MRIvol2volGCAM(): src-vsm dim mismatch\n");
      return(NULL);
    }
  }

  if(srclta){
    if(srclta->type != LINEAR_VOX_TO_VOX){
      printf("MRIvol2volGCAM(): Chaning Source LTA type to vox2vox\n");
      LTAchangeType(srclta, LINEAR_VOX_TO_VOX) ;
    }
    if(LTAmriIsSource(srclta, src)){
      printf("MRIvol2volGCAM(): Inverting Source LTA\n");
      Vsrc = MatrixInverse(srclta->xforms[0].m_L,NULL);
    }
    else Vsrc = MatrixCopy(srclta->xforms[0].m_L,NULL);
  }
  else Vsrc = MatrixIdentity(4,NULL);

  vgdst_src = &(dstlta->xforms[0].src);
  vgdst_dst = &(dstlta->xforms[0].dst);
  // check that vgdst_src is 256^3, 1mm
  if(vgdst_src->width != 256 || vgdst_src->height != 256 ||
     vgdst_src->depth != 256 || vgdst_src->xsize != 1 ||
     vgdst_src->ysize != 1   || vgdst_src->zsize != 1){
    if(vgdst_dst->width != 256 || vgdst_dst->height != 256 ||
       vgdst_dst->depth != 256 || vgdst_dst->xsize != 1 ||
       vgdst_dst->ysize != 1   || vgdst_dst->zsize != 1)
      if(getenv("MY_MORPHS_DO_NOT_CONFORM_DEAL_WITH_IT") == NULL) {
        printf("ERROR: MRIvol2volGCAM(): neither src nor dst VG of Dest LTA is conformed\n");
        return(NULL);
      }
      else
        printf("WARN: MRIvol2volGCAM(): neither src nor dst VG of Dest LTA is conformed\n");
    else {
      printf("MRIvol2volGCAM(): Inverting Destination LTA\n");
      LTAinvert(dstlta,dstlta);
    }
  }
  else printf("MRIvol2volGCAM(): NOT inverting Destination LTA\n");
  vgdst_src = &(dstlta->xforms[0].src);
  vgdst_dst = &(dstlta->xforms[0].dst);
  if(dstlta->type != LINEAR_VOX_TO_VOX){
    printf("MRIvol2volGCAM(): Chaning Destination LTA type to vox2vox\n");
    LTAchangeType(dstlta, LINEAR_VOX_TO_VOX) ;
  }
  Vdst = MatrixInverse(dstlta->xforms[0].m_L,NULL);

  if(dst == NULL){
    dst = MRIallocFromVolGeom(vgdst_dst, MRI_FLOAT, src->nframes, 0);
    if(dst==NULL) return(NULL);
    MRIcopyPulseParameters(src,dst);
  }

  if(sample_type == SAMPLE_CUBIC_BSPLINE) bspline = MRItoBSpline(src,NULL,3);
  valvect = (float *) calloc(sizeof(float),src->nframes);

  crsDst = MatrixAlloc(4,1,MATRIX_REAL);
  crsDst->rptr[4][1] = 1;
  crsAnat = MatrixAlloc(4,1,MATRIX_REAL);
  crsAnat->rptr[4][1] = 1;
  // scroll thru the CRS in the output/dest volume
  timer.reset();
  for(c=0; c < dst->width; c++){
    for(r=0; r < dst->height; r++){
      for(s=0; s < dst->depth; s++){
	// CRS in destination volume
	crsDst->rptr[1][1] = c;
	crsDst->rptr[2][1] = r;
	crsDst->rptr[3][1] = s;

	// Compute the CRS in the GCAM
	crsGCAM = MatrixMultiplyD(Vdst,crsDst,crsGCAM);

	if(gcam){
	  // Compute the CRS in the anatomical
	  out_of_gcam = GCAMsampleMorph(gcam, 
					crsGCAM->rptr[1][1],crsGCAM->rptr[2][1],crsGCAM->rptr[3][1],
					&crsAnat->rptr[1][1],&crsAnat->rptr[2][1],&crsAnat->rptr[3][1]);
	  if(out_of_gcam) continue;
	}
	else crsAnat = MatrixCopy(crsGCAM,crsAnat);

	// Compute the CRS in the Source Space
	crsSrc = MatrixMultiply(Vsrc,crsAnat,crsSrc);

        if(vsm){
          /* crsSrc is the CRS in the undistored source space. This
	     code computes the crsSrc in the space distorted by B0
	     inhomogeneity (ie, the space of MRI *src).  This is just
	     a change in the row value as given by the voxel shift map
	     (VSM). The VSM must have the same dimensions as src. */
	  cvsm = floor(crsSrc->rptr[1][1]);
	  rvsm = floor(crsSrc->rptr[2][1]);
	  iss = nint(crsSrc->rptr[3][1]);

	  if(cvsm < 0 || cvsm+1 >= src->width)  continue;
	  if(rvsm < 0 || rvsm+1 >= src->height) continue;
	  if(iss < 0  || iss+1  >= src->depth) continue;
	  // Dont sample outside the BO mask indicated by vsm=0
	  v = MRIgetVoxVal(vsm,cvsm,rvsm,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm+1,rvsm,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm,rvsm+1,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm+1,rvsm+1,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  /* Performs 3D interpolation. May want to use iss instead of crsSrc->rptr[3][1]
	     to make it a 2D interpolation. Not sure.*/
          MRIsampleSeqVolume(vsm, crsSrc->rptr[1][1],crsSrc->rptr[2][1],crsSrc->rptr[3][1], &drvsm, 0, 0);
	  if(drvsm == 0) continue;
          crsSrc->rptr[2][1] += drvsm;
        }

	// Check for out of the source FoV
	if(crsSrc->rptr[1][1] < 0 || crsSrc->rptr[1][1] >= src->width)  continue;
	if(crsSrc->rptr[2][1] < 0 || crsSrc->rptr[2][1] >= src->height) continue;
	if(crsSrc->rptr[3][1] < 0 || crsSrc->rptr[3][1] >= src->depth)  continue;

        if(sample_type != SAMPLE_TRILINEAR)
	  for(f=0; f < src->nframes; f++){
	    if(sample_type == SAMPLE_CUBIC_BSPLINE)
	      MRIsampleBSpline(bspline, crsSrc->rptr[1][1],crsSrc->rptr[2][1],crsSrc->rptr[3][1], f, &val);
	    else
	      MRIsampleVolumeFrameType(src,
				       crsSrc->rptr[1][1],crsSrc->rptr[2][1],crsSrc->rptr[3][1],
				       f, sample_type, &val) ;
	    MRIsetVoxVal(dst,c,r,s,f, val);
	  }
	else {
	  // This will do the same as above, it is just faster with multiple frames
          MRIsampleSeqVolume(src, crsSrc->rptr[1][1],crsSrc->rptr[2][1],crsSrc->rptr[3][1],
			     valvect,0, src->nframes-1) ;
	  for(f=0; f < src->nframes; f++)  MRIsetVoxVal(dst,c,r,s,f, valvect[f]);
	}

      } // s
    } // r
  } // c

  MatrixFree(&crsDst);
  MatrixFree(&crsGCAM);
  MatrixFree(&crsAnat);
  MatrixFree(&crsSrc);
  MatrixFree(&Vdst);
  MatrixFree(&Vsrc);
  if(bspline) MRIfreeBSpline(&bspline);
  free(valvect);

  printf("MRIvol2volGCAM: t=%6.4f\n", timer.seconds());
  fflush(stdout);

  return(dst);
}

#if 0
  if(0){
    vgdst_src = &(dstlta->xforms[0].src);
    vgdst_dst = &(dstlta->xforms[0].dst);
    // check that vgdst_src is 256^3, 1mm
    if(vgdst_src->width != 256 || vgdst_src->height != 256 ||
       vgdst_src->depth != 256 || vgdst_src->xsize != 1 ||
       vgdst_src->ysize != 1   || vgdst_src->zsize != 1){
      if(vgdst_dst->width != 256 || vgdst_dst->height != 256 ||
	 vgdst_dst->depth != 256 || vgdst_dst->xsize != 1 ||
	 vgdst_dst->ysize != 1   || vgdst_dst->zsize != 1)
	if(getenv("MY_MORPHS_DO_NOT_CONFORM_DEAL_WITH_IT") == NULL) {
	  printf("ERROR: MRIvol2volGCAM(): neither src nor dst VG of Dest LTA is conformed\n");
	  return(NULL);
	}
	else
	  printf("WARN: MRIvol2volGCAM(): neither src nor dst VG of Dest LTA is conformed\n");
      else {
	printf("MRIvol2volGCAM(): Inverting Destination LTA\n");
	LTAinvert(dstlta,dstlta);
      }
    }
    else printf("MRIvol2volGCAM(): NOT inverting Destination LTA\n");
  }


#endif
