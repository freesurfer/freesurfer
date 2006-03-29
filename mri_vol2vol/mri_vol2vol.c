/*
  Name:    mri_vol2vol
  Author:  Douglas N. Greve 
  email:   analysis-bugs@nmr.mgh.harvard.edu
  Date:    2/27/02
  Purpose: converts values in one volume to another volume
  $Id: mri_vol2vol.c,v 1.17 2006/03/29 07:15:35 greve Exp $

*/

/*
BEGINHELP --------------------------------------------------------------

Resamples a volume into another field-of-view. 

FLAGS AND ARGUMENTS

--i invol <fmt>

Input volume, and, optionally, the intput format.  If the format is
not included, the format will be inferred from the path name. See
FORMATS below.

--o outvol <fmt>

Path name of the output volume. If the format is not included, the
format will be inferred from the path name. See FORMATS below.
A register.dat-style file will also be produced. It will be called
outvol.reg.

--t templatevol <fmt>

This is the volume that will be used as a geometry template for the
output volume in terms of the field-of-view, geometry, and
precision. Some of the template parameters can be overridden as
explained below. If the format is not included, the format will be
inferred from the path name. See FORMATS below. If no template volume
is specified, the input is used as the template.

--xfm in-to-temp-xfmfile
--reg temp-to-in-xfmfile

The xfmfile contains the matrix that maps the XYZ of one volume to
that of another. There can be some confusion as to which direction
this goes. For historical reasons, matrices that are compatible with
tkregister2 actually map from the targ/template to the mov/input. But
one can also use matrices that were created with the MNI mritotal or
Register tools which more reasonably map from input to template. When
using tkregister2 xfmfiles, specify them with --reg and mri_vol2vol
will automatically invert them to go in the right direction. When
using an xfmfile that already goes in the right direction, pass it
with --xfm. 

--invxfm

Invert the xfm matrix before applying transform. Will not have an 
effect with --reg.

--noinvxfm

Do not invert the xfm matrix before applying transform. 

--fstal

Resample the input volume to talairach space. The xfm file must be a
register.dat-format matrix. The talairach matrix is obtained from
talairach.xfm from SUBJECTS_DIR/subjid/transforms. SUBJECTS_DIR is
read from the environment or can be specified with --sd. subjid is
read from the xfm file. The transformation matrix is then computed as
inv(R*inv(Xtal)*inv(Rtal)), where Xtal is talairach.xfm matrix, R is
the matrix in the xfm file, and Rtal maps from the talairach COR FOV
to a reduced FOV that covers only the brain. Reducing the FOV saves
space relative to the 256^3 COR FOV. By default, the output will be
2mm isotropic, but this can be changed with --fstalres. Specify the 
xfm with --reg. If you want to go from talairach space back to the input,
then specify --noinvxfm. The talairach subject is assumed to be fsaverage.

--fstalres resmm

Set the resolution of the output when using --fstal. By default, it
is 2 mm, but can be changed to 1.0 mm.

--interp method

Interpolate the output based on the given method. Legal values are:
nearest, trilin, and sinc. trilin is the default. sinc requires one
parameter (hw). sinc probably does not work.

--s subjectname

Subject name for output registration file. This has not effect on the
actual reslicing but can have an effect when the registration file
is used in subsequent processing.

--precision precisionid 

Set output precision to precisionid. Legal values are uchar, short,
int, long, and float. Default is float.

--precision-temp 

Set output precision to be that of the template. Overrides default of float.

--voxres colres rowres sliceres 

Set output voxel resolution (in mm). Overrides voxres in template.
Keeps the same field-of-view, so the output dimension is automatically
changed accordingly (so it cannot be used with --voxdim).

--voxres-in-plane colres rowres 

Same as --voxres but only on the in-plane components.

--voxdim ncols nrows nslices 

Set output voxel dimension. Overrides voxres in template. Keeps
the same field-of-view, so the output resolution is automatically
changed accordingly (so it cannot be used with --voxres).

--voxdim-in-plane ncols nrows 

Same as --voxdim but only on the in-plane components.

--irescale min max

Rescale intensity to be between min and max. This can be useful when
the output precision is less than that of the input. When the output
is COR, the values are automatically rescaled to 0 to 255.



--gdiagno diagnostic level

Sets the diagnostic level (only good for debuggin').

--version

Print out version string and exit.

--help 

Prints out all this information.

ALGORITH/TKREGISTER MATRIX CONENTION 

To convert a volume from one space/FOV to another, one needs to know
how to convert CRS (ie, col, row, slice) in the target FOV to that
in the source. This is referred to as the voxel-to-voxel transform,
and its matrix is called V.

CRSin = V * CRSout
V = inv(Tin*X)*Tout

where T is a matrix that converts CRS to XYZ. X is the matrix
specified with -xfm and maps from XYZin to XYZout. The X matrix is
only meaningful in terms of what Tin and Tout are.  The TkRegister
Convention defines T to be:

T = [-dc  0   0  Nc/2
      0   0  ds -Ns/2
      0 -dr   0  Nr/2
      0   0   0  1];

where dc, dr, and ds are the resolutions of the columns, rows, and 
slices, respectively, and Nc, Nr, and Ns are the number of columns,
rows, and slices, respectively. Column is the fastest dimension,
Row is the next fastest, and Slice is the slowest.

EXAMPLES:

If a functional volume is f.bhdr (or f.nii.gz, or f.mgh, etc), and the
subject is bert, and the registration file is register.dat, then
running the following command should show that they are in
registration:

tkregister2 --reg register.dat --mov f.bhdr

If they are not, then fix it because nothing below is going to work.

1. Resample a functional volume to talairach space at 1 mm res

   mri_vol2vol --i f.bhdr --o ftal.mgh \
     --reg register.dat --fstal --fstalres 1

   register.dat registers the subject anatomical and functional. Note
   that a template is not needed (it is generated internally). The
   registration of the ftal volume with the talairach subject can then
   be checked with: tkregister2 --mov ftal.mgh --reg ftal.mgh.reg. This
   is specially designed to be in registration with a FreeSurfer average
   subject (ie, one created by make_average_subject) such as the default
   fsaverage. Accordingly, you do not need a registration file when
   displaying ftal.mgh as an overlay on the average subject, eg:
     tkmedit fsaverage orig.mgz -overlay ftal.mgh
   will overlay ftal.mgh correctly.

2. Resample an anatomical volume into the functional space with a
   1 mm in-plane resolution:

  mri_vol2vol --i $SUBJECTS_DIR/mysubj/mri/orig.mgz
    --t mysubjsess/bold/001/f.bhdr --s mysubj 
    --reg  mysubjsess/bold/register.dat --noinvxfm 
    --o  mysubjsess/anatfunc/999/f.mgh
    --voxres-in-plane 1 1

3. Resample a subcortical segmentation to functional space. NOTE: THIS
   IS NOT THE RECOMMENDED METHOD FOR THIS. TRY USING mri_label2vol
   INSTEAD.  This uses nearest-neighbor interp because the
   segmentation values are categorical, not continuous.

   mri_vol2vol --i $SUBJECTS_DIR/subjid/mri/aseg.mgz
               --o aseg.mgh
               --t func.bhdr
               --reg  register.dat --noinvxfm
               --interp nearest

4. Resample a structural volume to talairach space.
   NOTE: BUG!! THIS WILL LIKELY NOT WORK!!!!
   cd  $SUBJECTS_DIR/subjid/mri/
   mri_vol2vol --i orig.mgz --o orig-tal.mgz
               --t $SUBJECTS_DIR/fsaverage/mri/orig.mgz
               --xfm transforms/talairach.xfm 

   NOTE: this should give the same result as:
   mri_convert orig.mgz orig-tal.mgz --apply_transform transforms/talairach.xfm 

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
static int istringnmatch(char *str1, char *str2, int n);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_vol2vol.c,v 1.17 2006/03/29 07:15:35 greve Exp $";
char *Progname = NULL;

int debug = 0, gdiagno = -1;

char *tempvolpath;
char *tempvolfmt = NULL;
int   tempvolfmtid = 0;

char *outvolpath;
char *outvolfmt = NULL;
int   outvolfmtid = 0;
char *outprecision = "float";
int   outprecisioncode;

float outvoxres[3];
int   force_outvoxres = 0;
int   force_outvoxres_in_plane = 0;

int   outvoxdim[3];
int   force_outvoxdim = 0;
int   force_outvoxdim_in_plane = 0;

float outcenter[3];
int   force_outcenter = 0;
float shiftcenter[3] = {0,0.0};

char *involpath;
char *involfmt = NULL;
int   involfmtid = 0;

char *xfmfile = NULL;
int  invertxfm = 0;
char *interpmethod = "trilinear";
int   interpcode = 0;
int   sinchw;

MRI *InVol;
MRI *TempVol,*OutVol;
MRI *tmpmri;

MATRIX *Vt2s, *X, *invX, *Xtmp, *invXtal;
MATRIX *Tin, *Tout, *invTin;
MATRIX *Xtal, *R, *Rtal, *invRtal;
MATRIX *Ttemp, *Mtemp, *D;
MATRIX *invTtemp;

char *FSH=NULL;
char *subjectsdir = NULL;   /* SUBJECTS_DIR */
char *talxfmfile = "talairach.xfm";
char talxfmpath[1000];
int fstalairach = 0;
char *talsubject = NULL;
int talres = 2; // Can only be 1 or 2
char *subject = NULL;

int dont_irescale = 1;
float minrescale = 0.0, maxrescale = 255.0;
char fname[1000], dname[1000];
int ModInput = 0;
float ipr, bpr, intensity, xfov, yfov, zfov;
int float2int,err, nargs;


/*---------------------------------------------------------------*/
int main(int argc, char **argv)
{
  char *trgsubject;
  char regfile[1000];
  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string (argc, argv, "$Id: mri_vol2vol.c,v 1.17 2006/03/29 07:15:35 greve Exp $", "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_vol2vol.c,v 1.17 2006/03/29 07:15:35 greve Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  FSH = getenv("FREESURFER_HOME");
  if(FSH==NULL){
    printf("ERROR: FREESURFER_HOME undefined.\n");
    exit(1);
  }

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  if(gdiagno > -1) Gdiag_no = gdiagno;

  check_options();

  /* Read in the template volume header */
  TempVol = MRIreadHeader(tempvolpath,tempvolfmtid);
  if(TempVol == NULL){
    printf("ERROR: reading %s header\n",tempvolpath);
    exit(1);
  }

  if(outprecision != NULL){
    outprecisioncode = MRIprecisionCode(outprecision);
    if(outprecisioncode < 0){
      printf("ERROR: precision %s unrecognized\n",outprecision);
      printf("       legal values are uchar, short, int, long, and float\n");
      exit(1);
    }
  }
  else{
    outprecisioncode = TempVol->type;
    outprecision = MRIprecisionString(TempVol->type);
  }
  if(outvolfmtid == MRI_CORONAL_SLICE_DIRECTORY){
    printf("INFO: forcing output to be uchar for COR\n");
    outprecisioncode = MRI_UCHAR;
    outprecision = MRIprecisionString(MRI_UCHAR);
  }


  xfov = TempVol->xsize * TempVol->width;
  yfov = TempVol->ysize * TempVol->height;
  zfov = TempVol->zsize * TempVol->depth;
  if(force_outvoxres || force_outvoxres_in_plane){
    // Change pixel size while keeping FOV
    TempVol->xsize = outvoxres[0];
    TempVol->ysize = outvoxres[1];
    TempVol->width  = (int)(round(xfov/TempVol->xsize));
    TempVol->height = (int)(round(yfov/TempVol->ysize));
    if(! force_outvoxres_in_plane){
      TempVol->zsize = outvoxres[2];
      TempVol->depth  = (int)(round(zfov/TempVol->zsize));
    }
  }
  if(force_outvoxdim || force_outvoxdim_in_plane){
    // Change dimension while keeping FOV
    TempVol->width  = outvoxdim[0];
    TempVol->height = outvoxdim[1];
    TempVol->xsize = xfov/TempVol->width;
    TempVol->ysize = yfov/TempVol->height;
    if(! force_outvoxdim_in_plane){
      TempVol->depth  = outvoxdim[2];
      TempVol->zsize = zfov/TempVol->depth;
    }
  }
  if(force_outcenter){
    /* This does nothing */
    TempVol->c_r = outcenter[0];
    TempVol->c_a = outcenter[1];
    TempVol->c_s = outcenter[2];
  }

  /* Fix the tkregister matrix */
  //Ma2v = MRIfixTkReg(RefAnat, TempVol, Ma2vTKR);

  /* Dump some info before staring the main program */
  dump_options(stdout);

  printf("INFO: reading  %s as %s\n",involpath,involfmt);
  InVol =  MRIreadType(involpath,involfmtid);
  if(InVol == NULL){
    printf("ERROR: could not read %s as %s\n",involpath,involfmt);
    exit(1);
  }
  Tin  = MRIxfmCRS2XYZtkreg(InVol); 
  MRIaddCommandLine(InVol, cmdline) ;
  printf("Done loading source volume \n");

  /* --------- read in transform ------------------*/
  if(xfmfile != NULL){
    printf("INFO: reading xfm file %s, trying as reg.dat \n",xfmfile);
    err = regio_read_register(xfmfile, &trgsubject, &ipr, &bpr, 
			      &intensity, &X, &float2int);
    printf("Input XFM: ----------------\n");
    MatrixPrint(stdout,X);
    if(!err){
      if(float2int == FLT2INT_TKREG){
	printf("INFO: Fixing tkregister matrix\n");
	printf("Original Reg Matrix: ----------------\n");
	MatrixPrint(stdout,X);
	Xtmp = MRIfixTkReg(InVol,X);
	MatrixFree(&X);
	X = Xtmp;
	printf("New Reg Matrix: ----------------\n");
	MatrixPrint(stdout,X);
      }
      if(fstalairach){
	printf("INFO: recomputing xfm for talairach transform\n");
	Xtal = DevolveXFM(trgsubject, NULL, NULL);
	invXtal = MatrixInverse(Xtal,NULL);
	printf("Xtal: ------------------------------\n");
	MatrixPrint(stdout,Xtal);
	printf("inv(Xtal): ------------------------------\n");
	MatrixPrint(stdout,invXtal);
	// X = X*inv(Xtal)*inv(Rtal)
	X = MatrixMultiply(X,invXtal,X); 
	invRtal = MatrixInverse(Rtal,NULL);
	X = MatrixMultiply(X,invRtal,X); 
      }
    }
    else{
      if(fstalairach){
	printf("ERROR: xfmfile must be a register.dat-format with "
	       "--fstalairach\n");
	exit(1);
      }
      printf("INFO: reading xfm file %s, trying as MINC xfm \n",xfmfile);
      err = regio_read_mincxfm(xfmfile, &X,NULL);
      if(err) exit(1);
    }
  }
  else X = MatrixIdentity(4,NULL);
  if(invertxfm) MatrixInverse(X,X);

  printf("Final XFM (including any inversion) ----------------\n");
  MatrixPrint(stdout,X);

  if(ModInput){
    printf("Modifying input vox2ras instead of resampling\n");
    Ttemp = MRIxfmCRS2XYZtkreg(TempVol); 
    invTtemp = MatrixInverse(Ttemp,NULL);
    Mtemp = MRIxfmCRS2XYZ(TempVol,0); 
    D = MatrixMultiply(Mtemp,invTtemp,NULL);
    D = MatrixMultiply(D,X,D);
    D = MatrixMultiply(D,Tin,D);
    printf("D ---------------------\n");
    MatrixPrint(stdout,D);
    MRIsetVoxelToRasXform(InVol,D);
    printf("INFO: writing output volume to %s (%s)\n",
	   outvolpath,outvolfmt);
    MRIwriteType(InVol,outvolpath,outvolfmtid);
    exit(1);
  }

  /*---------- Allocate the output volume ------------------*/
  OutVol = MRIallocSequence(TempVol->width, TempVol->height, 
          TempVol->depth, MRI_FLOAT, InVol->nframes );
  if(OutVol == NULL){
    printf("ERROR: could not alloc output volume MRI\n");
    exit(1);
  }
  MRIcopyHeader(TempVol,OutVol);
  OutVol->nframes = InVol->nframes;
  OutVol->type = MRI_FLOAT; /*Keep float until the end*/


  /* Construct the matrix to map from out CRS to in CRS */
  /* Vt2s = inv(Tin)*inv(X)*Tout */
  Tin  = MRIxfmCRS2XYZtkreg(InVol); 
  Tout = MRIxfmCRS2XYZtkreg(OutVol);
  /* Apply shift */
  //Tout->rptr[1][4] += shiftcenter[0];
  //Tout->rptr[2][4] += shiftcenter[1];
  //Tout->rptr[3][4] += shiftcenter[2];
  X->rptr[1][4] += shiftcenter[0];
  X->rptr[2][4] += shiftcenter[1];
  X->rptr[3][4] += shiftcenter[2];

  invTin = MatrixInverse(Tin,NULL);
  //R = MRItkRegMtx(OutVol, InVol, NULL);
  invX = MatrixInverse(X,NULL);
  Vt2s = MatrixMultiply(invTin,invX,NULL); 
  MatrixMultiply(Vt2s,Tout,Vt2s); 

  printf("Tin: ------------------------------\n");
  MatrixPrint(stdout,Tin);
  printf("Tout: ------------------------------\n");
  MatrixPrint(stdout,Tout);
  printf("X: ------------------------------\n");
  MatrixPrint(stdout,X);
  printf("invX: ------------------------------\n");
  MatrixPrint(stdout,MatrixInverse(X,NULL));
  printf("OutVox to InVox XFM: ------------------------------\n");
  MatrixPrint(stdout,Vt2s);
  printf("--------------------------------------------------\n");

  printf("INFO: resampling volume to volume\n");  
  //MRIvol2Vol(InVol,OutVol,Vt2s,interpcode,-1,-1,sinchw);
  MRIvol2Vol(InVol,OutVol,Vt2s,interpcode,sinchw);
  OutVol->imnr0 = 1;
  OutVol->imnr1 = OutVol->depth;

  if(outvolpath != NULL){
    if(OutVol->type != outprecisioncode){
      printf("INFO: changing type to %s\n",
	     MRIprecisionString(outprecisioncode));
      printf("outtype = %d, intype = %d (uchar = %d)\n",
	     OutVol->type,InVol->type,MRI_UCHAR);
      if(outprecisioncode == MRI_UCHAR && InVol->type != MRI_UCHAR &&
	 dont_irescale == 1){
	printf("INFO: forcing rescale to 0 to 255 for uchar ouptput\n");
	dont_irescale = 0;
	minrescale = 0;
	maxrescale = 255;
      }
      printf("dont_irescale = %d\n",dont_irescale);
      tmpmri = MRISeqchangeType(OutVol, outprecisioncode, 
				minrescale, maxrescale, dont_irescale);
      if(tmpmri == NULL){
	printf("ERROR: changing type\n");
	exit(1);
      }
      MRIfree(&OutVol);
      OutVol = tmpmri;
    }

    printf("INFO: writing output volume to %s (%s)\n",
	   outvolpath,outvolfmt);
    MRIwriteType(OutVol,outvolpath,outvolfmtid);

    if(outvolfmtid == MRI_CORONAL_SLICE_DIRECTORY)
      sprintf(regfile,"%s/COR.reg",outvolpath);
    if(outvolfmtid == BSHORT_FILE || outvolfmtid == BFLOAT_FILE){
      decompose_b_fname(outvolpath, dname, fname);
      sprintf(regfile,"%s/%s.reg",dname,fname);
    }
    else  sprintf(regfile,"%s.reg",outvolpath);
    printf("INFO: writing registration matrix to %s\n",regfile);
    if(fstalairach){
      R = Rtal;
      subject = talsubject;
    }
    else  {
      R = MatrixIdentity(4,NULL);
      // On 6/3/28, it used to be this:
      //    R = MRItkRegMtx(TempVol,InVol,NULL);
      // I have no idea why it was that, seems clearly wrong
    }
    regio_write_register(regfile,subject,OutVol->xsize,
			 OutVol->zsize,1,R,FLT2INT_ROUND);

    printf("To check registration, run:\n");
    printf("\n");
    printf("  tkregister2 --mov %s --reg %s ",
	   outvolpath,regfile);

    if(!fstalairach) printf("--targ %s\n",tempvolpath);
    else printf("\n");

    printf("\n");
    printf("\n");

  }

  printf("done\n");
  return(0);
}


/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
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
    else if (!strcasecmp(option, "--fstal"))       fstalairach = 1;
    else if (!strcasecmp(option, "--fstalairach")) fstalairach = 1;
    else if (!strcasecmp(option, "--invxfm"))   invertxfm = 1;
    else if (!strcasecmp(option, "--noinvxfm"))   invertxfm = 0;
    else if (!strcasecmp(option, "--modinput"))  ModInput = 1;
    else if (istringnmatch(option, "--precision-temp",0)) outprecision=NULL;

    else if ( !strcmp(option, "--gdiagno") ) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&gdiagno);
      nargsused = 1;
    }

    else if(istringnmatch(option, "--in",0) || istringnmatch(option, "--i",0)){
      if(nargc < 1) argnerr(option,1);
      involpath = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	involfmt = pargv[1]; nargsused ++;
	involfmtid = string_to_type(involfmt);
      }
    }
    else if(istringnmatch(option, "--out",0) || istringnmatch(option, "--o",0)){
      if(nargc < 1) argnerr(option,1);
      outvolpath = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	outvolfmt = pargv[1]; nargsused ++;
	outvolfmtid = string_to_type(outvolfmt);
      }
    }
    else if (istringnmatch(option, "--precision",0)){
      if(nargc < 1) argnerr(option,1);
      outprecision = pargv[0]; nargsused = 1;
    }
    else if (istringnmatch(option, "--voxres",0)){
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%f",&outvoxres[0]);
      sscanf(pargv[1],"%f",&outvoxres[1]);
      sscanf(pargv[2],"%f",&outvoxres[2]);
      force_outvoxres = 1;
      nargsused = 3;
    }
    else if (istringnmatch(option, "--voxres-in-plane",0)){
      if(nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%f",&outvoxres[0]);
      sscanf(pargv[1],"%f",&outvoxres[1]);
      force_outvoxres_in_plane = 1;
      nargsused = 2;
    }
    else if (istringnmatch(option, "--voxdim",0)){
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%d",&outvoxdim[0]);
      sscanf(pargv[1],"%d",&outvoxdim[1]);
      sscanf(pargv[2],"%d",&outvoxdim[2]);
      force_outvoxdim = 1;
      nargsused = 3;
    }
    else if (istringnmatch(option, "--voxdim-in-plane",0)){
      if(nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%d",&outvoxdim[0]);
      sscanf(pargv[1],"%d",&outvoxdim[1]);
      force_outvoxdim_in_plane = 1;
      nargsused = 2;
    }
    else if (istringnmatch(option, "--center",0)){
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%f",&outcenter[0]);
      sscanf(pargv[1],"%f",&outcenter[1]);
      sscanf(pargv[2],"%f",&outcenter[2]);
      force_outcenter = 1;
      nargsused = 3;
    }
    else if (istringnmatch(option, "--shift",0)){
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%f",&shiftcenter[0]);
      sscanf(pargv[1],"%f",&shiftcenter[1]);
      sscanf(pargv[2],"%f",&shiftcenter[2]);
      nargsused = 3;
    }
    else if(istringnmatch(option, "--template",6) || istringnmatch(option, "--t",6)){
      if(nargc < 1) argnerr(option,1);
      tempvolpath = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	tempvolfmt = pargv[1]; nargsused ++;
	tempvolfmtid = string_to_type(tempvolfmt);
      }
    }

    else if (istringnmatch(option, "--sd",4)){
      if(nargc < 1) argnerr(option,1);
      subjectsdir = pargv[0]; nargsused = 1;
    }

    else if (istringnmatch(option, "--irescale",10)){
      if(nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%f",&minrescale); 
      sscanf(pargv[1],"%f",&maxrescale); 
      dont_irescale = 0;
      if(minrescale > maxrescale){
	printf("ERROR: rescale min = %g > max = %g\n",
	       minrescale, maxrescale);
	exit(1);
      }
      nargsused = 2;
    }

    else if (istringnmatch(option, "--xfm",8)){
      if(nargc < 1) argnerr(option,1);
      xfmfile = pargv[0]; nargsused = 1;
    }
    else if (istringnmatch(option, "--reg",8)){
      if(nargc < 1) argnerr(option,1);
      xfmfile = pargv[0]; nargsused = 1;
      invertxfm = 1;
    }

    else if (istringnmatch(option, "--talres",8)){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&talres); 
      if(talres != 1 && talres != 2){
	printf("ERROR: tal res %d invalid. Only use 1 or 2\n",talres);
	exit(1);
      }
      nargsused = 1;
    }

    else if (istringnmatch(option, "--interp",8)){
      if(nargc < 1) argnerr(option,1);
      interpmethod = pargv[0]; nargsused = 1;
      if(!strcmp(interpmethod,"sinc") && nth_is_arg(nargc, pargv, 1)){
	sscanf(pargv[1],"%d",&sinchw); 
	nargsused ++;
      }
    }

    else if (istringnmatch(option, "--subject",3)){
      if(nargc < 1) argnerr(option,1);
      subject = pargv[0]; nargsused = 1;
    }

    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("  --in   input  volume path <fmt>\n");
  printf("  --out  output volume path <fmt>\n");
  printf("  \n");
  printf("  --temp output template volume <fmt>\n");
  //printf("  --conform : use COR as output template <fmt>\n");
  printf("  --precision precision : overrides default (float)\n");
  printf("  --precision-temp : uses precision of template\n");
  printf("  --voxres colres rowres sliceres : override template\n");
  printf("  --voxres-in-plane colres rowres : override template\n");
  printf("  --voxdim ncols  nrows  nslices : override template\n");
  printf("  --voxdim-in-plane ncols  nrows : override template\n");
  printf("  --center x y z : override template\n");
  printf("\n");
  printf("  --irescale min max : rescale intensities to min/max\n");
  printf("  \n");
  printf("  --xfm xfmfile : apply transform (assumes identity if not give)\n");
  printf("  --reg regfile : same as --xfm with --invxfm for register.dat\n");
  printf("  --fstal : resample volume into talairach space (needs xfm).\n");
  printf("  --talres : Output talairach resolution (1, 1.5, or <2>)\n");
  printf("  --invxfm : invert transform before applying\n");
  printf("  \n");
  printf("  --interp method : nearest, <trilin>, sinc \n");
  printf("  --s subjectname : subject name for output reg file\n");
  printf("  --modinput : change vox2ras, do not resample input \n");
  printf("  \n");
  printf("  --help    : hidden secrets of success\n");
  printf("  --gdiagno number : set diag level\n");
  printf("  --version : print version and exit\n");
  printf("  \n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;

  printf("\n%s\n\n",vcid);
printf("\n");
printf("Resamples a volume into another field-of-view. \n");
printf("\n");
printf("FLAGS AND ARGUMENTS\n");
printf("\n");
printf("--i invol <fmt>\n");
printf("\n");
printf("Input volume, and, optionally, the intput format.  If the format is\n");
printf("not included, the format will be inferred from the path name. See\n");
printf("FORMATS below.\n");
printf("\n");
printf("--o outvol <fmt>\n");
printf("\n");
printf("Path name of the output volume. If the format is not included, the\n");
printf("format will be inferred from the path name. See FORMATS below.\n");
printf("A register.dat-style file will also be produced. It will be called\n");
printf("outvol.reg.\n");
printf("\n");
printf("--t templatevol <fmt>\n");
printf("\n");
printf("This is the volume that will be used as a geometry template for the\n");
printf("output volume in terms of the field-of-view, geometry, and\n");
printf("precision. Some of the template parameters can be overridden as\n");
printf("explained below. If the format is not included, the format will be\n");
printf("inferred from the path name. See FORMATS below. If no template volume\n");
printf("is specified, the input is used as the template.\n");
printf("\n");
printf("--xfm in-to-temp-xfmfile\n");
printf("--reg temp-to-in-xfmfile\n");
printf("\n");
printf("The xfmfile contains the matrix that maps the XYZ of one volume to\n");
printf("that of another. There can be some confusion as to which direction\n");
printf("this goes. For historical reasons, matrices that are compatible with\n");
printf("tkregister2 actually map from the targ/template to the mov/input. But\n");
printf("one can also use matrices that were created with the MNI mritotal or\n");
printf("Register tools which more reasonably map from input to template. When\n");
printf("using tkregister2 xfmfiles, specify them with --reg and mri_vol2vol\n");
printf("will automatically invert them to go in the right direction. When\n");
printf("using an xfmfile that already goes in the right direction, pass it\n");
printf("with --xfm. \n");
printf("\n");
printf("--invxfm\n");
printf("\n");
printf("Invert the xfm matrix before applying transform. Will not have an \n");
printf("effect with --reg.\n");
printf("\n");
printf("--noinvxfm\n");
printf("\n");
printf("Do not invert the xfm matrix before applying transform. \n");
printf("\n");
printf("--fstal\n");
printf("\n");
printf("Resample the input volume to talairach space. The xfm file must be a\n");
printf("register.dat-format matrix. The talairach matrix is obtained from\n");
printf("talairach.xfm from SUBJECTS_DIR/subjid/transforms. SUBJECTS_DIR is\n");
printf("read from the environment or can be specified with --sd. subjid is\n");
printf("read from the xfm file. The transformation matrix is then computed as\n");
printf("inv(R*inv(Xtal)*inv(Rtal)), where Xtal is talairach.xfm matrix, R is\n");
printf("the matrix in the xfm file, and Rtal maps from the talairach COR FOV\n");
printf("to a reduced FOV that covers only the brain. Reducing the FOV saves\n");
printf("space relative to the 256^3 COR FOV. By default, the output will be\n");
printf("2mm isotropic, but this can be changed with --fstalres. Specify the \n");
printf("xfm with --reg. If you want to go from talairach space back to the input,\n");
printf("then specify --noinvxfm. The talairach subject is assumed to be fsaverage.\n");
printf("\n");
printf("--fstalres resmm\n");
printf("\n");
printf("Set the resolution of the output when using --fstal. By default, it\n");
printf("is 2 mm, but can be changed to 1.0 mm.\n");
printf("\n");
printf("--interp method\n");
printf("\n");
printf("Interpolate the output based on the given method. Legal values are:\n");
printf("nearest, trilin, and sinc. trilin is the default. sinc requires one\n");
printf("parameter (hw). sinc probably does not work.\n");
printf("\n");
printf("--s subjectname\n");
printf("\n");
printf("Subject name for output registration file. This has not effect on the\n");
printf("actual reslicing but can have an effect when the registration file\n");
printf("is used in subsequent processing.\n");
printf("\n");
printf("--precision precisionid \n");
printf("\n");
printf("Set output precision to precisionid. Legal values are uchar, short,\n");
printf("int, long, and float. Default is float.\n");
printf("\n");
printf("--precision-temp \n");
printf("\n");
printf("Set output precision to be that of the template. Overrides default of float.\n");
printf("\n");
printf("--voxres colres rowres sliceres \n");
printf("\n");
printf("Set output voxel resolution (in mm). Overrides voxres in template.\n");
printf("Keeps the same field-of-view, so the output dimension is automatically\n");
printf("changed accordingly (so it cannot be used with --voxdim).\n");
printf("\n");
printf("--voxres-in-plane colres rowres \n");
printf("\n");
printf("Same as --voxres but only on the in-plane components.\n");
printf("\n");
printf("--voxdim ncols nrows nslices \n");
printf("\n");
printf("Set output voxel dimension. Overrides voxres in template. Keeps\n");
printf("the same field-of-view, so the output resolution is automatically\n");
printf("changed accordingly (so it cannot be used with --voxres).\n");
printf("\n");
printf("--voxdim-in-plane ncols nrows \n");
printf("\n");
printf("Same as --voxdim but only on the in-plane components.\n");
printf("\n");
printf("--irescale min max\n");
printf("\n");
printf("Rescale intensity to be between min and max. This can be useful when\n");
printf("the output precision is less than that of the input. When the output\n");
printf("is COR, the values are automatically rescaled to 0 to 255.\n");
printf("\n");
printf("\n");
printf("\n");
printf("--gdiagno diagnostic level\n");
printf("\n");
printf("Sets the diagnostic level (only good for debuggin').\n");
printf("\n");
printf("--version\n");
printf("\n");
printf("Print out version string and exit.\n");
printf("\n");
printf("--help \n");
printf("\n");
printf("Prints out all this information.\n");
printf("\n");
printf("ALGORITH/TKREGISTER MATRIX CONENTION \n");
printf("\n");
printf("To convert a volume from one space/FOV to another, one needs to know\n");
printf("how to convert CRS (ie, col, row, slice) in the target FOV to that\n");
printf("in the source. This is referred to as the voxel-to-voxel transform,\n");
printf("and its matrix is called V.\n");
printf("\n");
printf("CRSin = V * CRSout\n");
printf("V = inv(Tin*X)*Tout\n");
printf("\n");
printf("where T is a matrix that converts CRS to XYZ. X is the matrix\n");
printf("specified with -xfm and maps from XYZin to XYZout. The X matrix is\n");
printf("only meaningful in terms of what Tin and Tout are.  The TkRegister\n");
printf("Convention defines T to be:\n");
printf("\n");
printf("T = [-dc  0   0  Nc/2\n");
printf("      0   0  ds -Ns/2\n");
printf("      0 -dr   0  Nr/2\n");
printf("      0   0   0  1];\n");
printf("\n");
printf("where dc, dr, and ds are the resolutions of the columns, rows, and \n");
printf("slices, respectively, and Nc, Nr, and Ns are the number of columns,\n");
printf("rows, and slices, respectively. Column is the fastest dimension,\n");
printf("Row is the next fastest, and Slice is the slowest.\n");
printf("\n");
printf("EXAMPLES:\n");
printf("\n");
printf("If a functional volume is f.bhdr (or f.nii.gz, or f.mgh, etc), and the\n");
printf("subject is bert, and the registration file is register.dat, then\n");
printf("running the following command should show that they are in\n");
printf("registration:\n");
printf("\n");
printf("tkregister2 --reg register.dat --mov f.bhdr\n");
printf("\n");
printf("If they are not, then fix it because nothing below is going to work.\n");
printf("\n");
printf("1. Resample a functional volume to talairach space at 1 mm res\n");
printf("\n");
printf("   mri_vol2vol --i f.bhdr --o ftal.mgh \\n");
printf("     --reg register.dat --fstal --fstalres 1\n");
printf("\n");
printf("   register.dat registers the subject anatomical and functional. Note\n");
printf("   that a template is not needed (it is generated internally). The\n");
printf("   registration of the ftal volume with the talairach subject can then\n");
printf("   be checked with: tkregister2 --mov ftal.mgh --reg ftal.mgh.reg. This\n");
printf("   is specially designed to be in registration with a FreeSurfer average\n");
printf("   subject (ie, one created by make_average_subject) such as the default\n");
printf("   fsaverage. Accordingly, you do not need a registration file when\n");
printf("   displaying ftal.mgh as an overlay on the average subject, eg:\n");
printf("     tkmedit fsaverage orig.mgz -overlay ftal.mgh\n");
printf("   will overlay ftal.mgh correctly.\n");
printf("\n");
printf("2. Resample an anatomical volume into the functional space with a\n");
printf("   1 mm in-plane resolution:\n");
printf("\n");
printf("  mri_vol2vol --i $SUBJECTS_DIR/mysubj/mri/orig.mgz\n");
printf("    --t mysubjsess/bold/001/f.bhdr --s mysubj \n");
printf("    --reg  mysubjsess/bold/register.dat --noinvxfm \n");
printf("    --o  mysubjsess/anatfunc/999/f.mgh\n");
printf("    --voxres-in-plane 1 1\n");
printf("\n");
printf("3. Resample a subcortical segmentation to functional space. NOTE: THIS\n");
printf("   IS NOT THE RECOMMENDED METHOD FOR THIS. TRY USING mri_label2vol\n");
printf("   INSTEAD.  This uses nearest-neighbor interp because the\n");
printf("   segmentation values are categorical, not continuous.\n");
printf("\n");
printf("   mri_vol2vol --i $SUBJECTS_DIR/subjid/mri/aseg.mgz\n");
printf("               --o aseg.mgh\n");
printf("               --t func.bhdr\n");
printf("               --reg  register.dat --noinvxfm\n");
printf("               --interp nearest\n");
printf("\n");
printf("4. Resample a structural volume to talairach space.\n");
printf("   NOTE: BUG!! THIS WILL LIKELY NOT WORK!!!!\n");
printf("   cd  $SUBJECTS_DIR/subjid/mri/\n");
printf("   mri_vol2vol --i orig.mgz --o orig-tal.mgz\n");
printf("               --t $SUBJECTS_DIR/fsaverage/mri/orig.mgz\n");
printf("               --xfm transforms/talairach.xfm \n");
printf("\n");
printf("   NOTE: this should give the same result as:\n");
printf("   mri_convert orig.mgz orig-tal.mgz --apply_transform transforms/talairach.xfm \n");
printf("\n");
printf("FORMATS\n");
printf("\n");
printf("Data file format can be specified implicitly (through the path name)\n");
printf("or explicitly. All formats accepted by mri_convert can be used. \n");
printf("\n");
printf("BUGS\n");
printf("\n");
printf("sinc interpolation is broken except for maybe COR to COR.\n");
printf("\n");
printf("\n");
printf("BUG REPORTING\n");
printf("\n");
printf("Report bugs to analysis-bugs@nmr.mgh.harvard.edu. Include the following \n");
printf("formatted as a list as follows: (1) command-line, (2) directory where\n");
printf("the program was run (for those in the MGH-NMR Center), (3) version, \n");
printf("(4) text output, (5) description of the problem.\n");
printf("\n");
printf("SEE ALSO \n");
printf("\n");
printf("mri_convert, tkregister2\n");
printf("\n");
printf("\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(involpath == NULL){
    printf("ERROR: No input supplied.\n");
    exit(1);
  }
  if(!ModInput && outvolpath == NULL){
    printf("ERROR: No output supplied.\n");
    exit(1);
  }

  if(involfmt == NULL){
    involfmtid = mri_identify(involpath);
    if(involfmtid == MRI_VOLUME_TYPE_UNKNOWN){
      printf("ERROR: cannot recognize the type of %s\n",involpath);
      exit(1);
    }
    involfmt = type_to_string(involfmtid);
  }
  else{
    involfmtid = string_to_type(involfmt);
    if(involfmtid == MRI_VOLUME_TYPE_UNKNOWN){
      printf("ERROR: cannot recognize format %s\n",involfmt);
      exit(1);
    }
  }

  if(outvolfmt == NULL){
    outvolfmtid = mri_identify(outvolpath);
    if(outvolfmtid == MRI_VOLUME_TYPE_UNKNOWN){
      printf("ERROR: cannot recognize the type of %s\n",outvolpath);
      exit(1);
    }
    outvolfmt = type_to_string(outvolfmtid);
  }
  else{
    outvolfmtid = string_to_type(outvolfmt);
    if(outvolfmtid == MRI_VOLUME_TYPE_UNKNOWN){
      printf("ERROR: cannot recognize format %s\n",outvolfmt);
      exit(1);
    }
  }

  if(fstalairach){
    if(tempvolpath != NULL){
      printf("ERROR: cannot specify --temp and --fstal\n");
      exit(1);
    }
    sprintf(fname,"%s/average/mni305.cor.subfov%d.mgz",FSH,(int)talres);
    tempvolpath = strcpyalloc(fname);
    sprintf(fname,"%s/average/mni305.cor.subfov%d.reg",FSH,(int)talres);
    err = regio_read_register(fname, &talsubject, &ipr, &bpr, 
			      &intensity, &Rtal, &float2int);
  }

  if(tempvolpath == NULL){
    printf("INFO: no template specified, using input as template\n");
    tempvolpath  = involpath;
    tempvolfmt   = involfmt;
    tempvolfmtid = involfmtid;
  }
  else{
    if(tempvolfmt == NULL){
      tempvolfmtid = mri_identify(tempvolpath);
      if(tempvolfmtid == MRI_VOLUME_TYPE_UNKNOWN){
	printf("ERROR: cannot recognize the type of %s\n",tempvolpath);
	exit(1);
      }
      tempvolfmt = type_to_string(tempvolfmtid);
    }
    else{
      tempvolfmtid = string_to_type(tempvolfmt);
      if(tempvolfmtid == MRI_VOLUME_TYPE_UNKNOWN){
	printf("ERROR: cannot recognize format %s\n",tempvolfmt);
	exit(1);
      }
    }
  }

  interpcode = MRIinterpCode(interpmethod);
  if(interpcode < 0){
    printf("ERROR: interpolation method %s unrecognized\n",interpmethod);
    printf("       legal values are nearest, trilin, and sinc\n");
    exit(1);
  }

  if(fstalairach && xfmfile == NULL){
    printf("ERROR: an xfmfile (register.dat) must be supplied when\n"
	   "       using --fstal\n");
    exit(1);
  }
  if(force_outvoxres && force_outvoxdim){
    printf("ERROR: cannot change ouput voxel resolution and output\n");
    printf("       voxel dimension\n");
    exit(1);
  }


  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  MRI *mri;

  fprintf(fp,"invol   path  %s\n",involpath);
  fprintf(fp,"outvol  path  %s\n",outvolpath);
  if(xfmfile) fprintf(fp,"xfm file    %s\n",xfmfile);
  fprintf(fp,"interp  %s (%d)\n",interpmethod,interpcode);
  fprintf(fp,"precision  %s (%d)\n",outprecision,outprecisioncode);
  if(interpcode == SAMPLE_SINC) fprintf(fp,"sinc hw  %d\n",sinchw);

  if(!tempvolpath) return;

  fprintf(fp,"template path  %s\n",tempvolpath);
  mri = TempVol;
  fprintf(fp, "%6.6s = %d\n", "height", mri->height);
  fprintf(fp, "%6.6s = %d\n", "width", mri->width);
  fprintf(fp, "%6.6s = %d\n", "depth", mri->depth);
  fprintf(fp, "%6.6s = %f\n", "xsize", mri->xsize);
  fprintf(fp, "%6.6s = %f\n", "ysize", mri->ysize);
  fprintf(fp, "%6.6s = %f\n", "zsize", mri->zsize);
  fprintf(fp, "%6.6s = %f %f %f\n", "cdc ", mri->x_r, mri->x_a, mri->x_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "rdc ", mri->y_r, mri->y_a, mri->y_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "sdc ", mri->z_r, mri->z_a, mri->z_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "xyz0", mri->c_r, mri->c_a, mri->c_s);
  
  fprintf(fp,"Gdiag_no  %d\n",Gdiag_no);

  return;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if(n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);
  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int isflag(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth)
{
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if(nargc <= nth) return(0); 

  /* check whether the nth arg is a flag */
  if(isflag(argv[nth])) return(0);

  return(1);
}

/*------------------------------------------------------------
  istringnmatch() - compare the first n characters of two strings,
  return a 1 if they match (ignoring case), a zero otherwise. If
  n=0, then do a full comparison.
  ------------------------------------------------------------*/
static int istringnmatch(char *str1, char *str2, int n)
{
  if(n > 0  && ! strncasecmp(str1,str2,n)) return(1);
  if(n <= 0 && ! strcasecmp(str1,str2)) return(1);
  return(0);
}
