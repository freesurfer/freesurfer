/*
  Name:    mri_vol2vol
  Author:  Douglas N. Greve 
  email:   analysis-bugs@nmr.mgh.harvard.edu
  Date:    2/27/02
  Purpose: converts values in one volume to another volume
  $Id: mri_vol2vol.c,v 1.1 2004/05/19 17:30:41 greve Exp $

  Things to do:
    1. Add ability to spec output center XYZ.
    2. Handle selxavg files (frame power).
    3. Test with no template.
    4. Add ability to spec another tal.xfm
    5. Unwarping
    6. Synthesizing
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
#include "mri2.h"
#include "mri_identify.h"
#include "MRIio_old.h"
#include "registerio.h"
#include "resample.h"

#ifdef X
#undef X
#endif

int MRIvol2Vol(MRI *src, MRI *targ, MATRIX *Vt2s, int InterpCode, 
	       int firstframe, int lastframe, float param);


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
static int istringnmatch(char *str1, char *str2, int n);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_vol2vol.c,v 1.1 2004/05/19 17:30:41 greve Exp $";
char *Progname = NULL;

int debug = 0, gdiagno = -1;

char *tempvolpath;
char *tempvolfmt = NULL;
int   tempvolfmtid = 0;

char *outvolpath;
char *outvolfmt = NULL;
int   outvolfmtid = 0;
char *outprecision = NULL;
int   outprecisioncode;
float outvoxres[3];
int   force_outvoxres = 0;
int    outvoxdim[3];
int   force_outvoxdim = 0;
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

MATRIX *Vt2s, *X, *invX, *Xtmp;
MATRIX *Tin, *Tout, *invTin;
MATRIX *Xtal, *R;

int fstalairach = 0;
char *subjectsdir = NULL;   /* SUBJECTS_DIR */
char *talxfmfile = "talairach.xfm";
char talxfmpath[1000];

int dont_irescale = 1;
float minrescale = 0.0, maxrescale = 255.0;

/*---------------------------------------------------------------*/
int main(int argc, char **argv)
{
  float ipr, bpr, intensity;
  int float2int,err;
  char *trgsubject;
  char regfile[1000];

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

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

  if(force_outvoxres){
    TempVol->xsize = outvoxres[0];
    TempVol->ysize = outvoxres[1];
    TempVol->zsize = outvoxres[2];
  }
  if(force_outvoxdim){
    TempVol->width = outvoxdim[0];
    TempVol->height = outvoxdim[1];
    TempVol->depth = outvoxdim[2];
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
  printf("Done loading source volume \n");

  /* --------- read in transform ------------------*/
  if(xfmfile != NULL){
    printf("INFO: reading xfm file %s, trying as reg.dat \n",xfmfile);
    err = regio_read_register(xfmfile, &trgsubject, &ipr, &bpr, 
			      &intensity, &X, &float2int);
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
	if(subjectsdir == NULL) subjectsdir = getenv("SUBJECTS_DIR");
	if (subjectsdir==NULL) {
	  printf("ERROR: SUBJECTS_DIR undefined. Use setenv or --sd\n");
	  exit(1);
	}
	sprintf(talxfmpath,"%s/%s/mri/transforms/%s",
		subjectsdir,trgsubject,talxfmfile);
	err = regio_read_mincxfm(talxfmpath, &Xtal);
	if(err) exit(1);
	invX = MatrixInverse(X,NULL);
        MatrixMultiply(Xtal,invX,X); 
	printf("Xtal: ------------------------------\n");
	MatrixPrint(stdout,Xtal);
	printf("inv(reg): ------------------------------\n");
	MatrixPrint(stdout,invX);
      }
    }
    else{
      if(fstalairach){
	printf("ERROR: xfmfile must be a register.dat-format with "
	       "--fstalairach\n");
	exit(1);
      }
      printf("INFO: reading xfm file %s, trying as MINC xfm \n",xfmfile);
      err = regio_read_mincxfm(xfmfile, &X);
      if(err) exit(1);
    }
  }
  else X = MatrixIdentity(4,NULL);


  /*---------- Allocate the output volume ------------------*/
  OutVol = MRIallocSequence(TempVol->width, TempVol->height, 
          TempVol->depth, MRI_FLOAT, InVol->nframes );
  if(OutVol == NULL){
    printf("ERROR: could not alloc output volume MRI\n");
    exit(1);
  }
  MRIcopyHeader(TempVol,OutVol);
  OutVol->nframes = InVol->nframes;

  if(invertxfm) MatrixInverse(X,X);

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
  MRIvol2Vol(InVol,OutVol,Vt2s,interpcode,-1,-1,sinchw);
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
      tmpmri = MRISeqchangeType(OutVol, TempVol->type, 
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
    if(fstalairach){
      R = MatrixIdentity(4,NULL);
      if(outvolfmtid == MRI_CORONAL_SLICE_DIRECTORY)
	sprintf(regfile,"%s/COR.reg",outvolpath);
      else
	sprintf(regfile,"%s.reg",outvolpath);
      printf("INFO: writing registration matrix to %s\n",regfile);
      regio_write_register(regfile,"talairach",OutVol->xsize,
			   OutVol->zsize,1,R,FLT2INT_ROUND);
    }

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
    else if (!strcasecmp(option, "--fstal"))   fstalairach = 1;
    else if (!strcasecmp(option, "--fstaliarach"))   fstalairach = 1;
    else if (!strcasecmp(option, "--invxfm"))   invertxfm = 1;

    else if ( !strcmp(option, "--gdiagno") ) {
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&gdiagno);
      nargsused = 1;
    }

    else if (istringnmatch(option, "--in",0)){
      if(nargc < 1) argnerr(option,1);
      involpath = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	involfmt = pargv[1]; nargsused ++;
	involfmtid = string_to_type(involfmt);
      }
    }
    else if (istringnmatch(option, "--out",0)){
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
    else if (istringnmatch(option, "--voxdim",0)){
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%d",&outvoxdim[0]);
      sscanf(pargv[1],"%d",&outvoxdim[1]);
      sscanf(pargv[2],"%d",&outvoxdim[2]);
      force_outvoxdim = 1;
      nargsused = 3;
    }

    else if (istringnmatch(option, "--template",6)){
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

    else if (istringnmatch(option, "--interp",8)){
      if(nargc < 1) argnerr(option,1);
      interpmethod = pargv[0]; nargsused = 1;
      if(!strcmp(interpmethod,"sinc") && nth_is_arg(nargc, pargv, 1)){
	sscanf(pargv[1],"%d",&sinchw); 
	nargsused ++;
      }
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
  printf("  --precision precision : overrides template precision\n");
  printf("  --voxres colres rowres sliceres : override template\n");
  printf("  --voxdim ncols  nrows  nslices : override template\n");
  printf("  --center x y z : override template\n");
  printf("\n");
  printf("  --irescale min max : rescale intensities to min/max\n");
  printf("  \n");
  printf("  --xfm xfmfile : apply transform \n");
  printf("  --fstal : resample volume into talairach space (needs xfm).\n");
  printf("  --invxfm : invert transform before applying\n");
  printf("  \n");
  printf("  --interp method : nearest, <trilin>, sinc \n");
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

  printf(

"Resamples a volume into another field-of-view. \n"
"\n"
"FLAGS AND ARGUMENTS\n"
"\n"
"--in input volume path <fmt>\n"
"\n"
"Input volume, and, optionally, the intput format.\n"
"If the format is not included, the format will be inferred from the\n"
"path name. See FORMATS below.\n"
"\n"
"--out output volume <fmt>\n"
"\n"
"Path name of the output volume. If the format is not included, the\n"
"format will be inferred from the path name. See FORMATS below.\n"
"\n"
"--temp template volume <fmt>\n"
"\n"
"This is the volume that will be used as a template for the output\n"
"volume in terms of the field-of-view, geometry, and precision. Some\n"
"of the template parameters can be overridden as explained below. If\n"
"the format is not included, the format will be inferred from the path\n"
"name. See FORMATS below. If no template volume is specified, the \n"
"input is used as the template.\n"
"\n"
"--precision precisionid \n"
"\n"
"Set output precision to precisionid. Legal values are uchar, short,\n"
"int, long, and float. Overrides precision in template.\n"
"\n"
"--voxres colres rowres sliceres \n"
"\n"
"Set output voxel resolution (in mm). Overrides voxres in template.\n"
"\n"
"--voxdim ncols nrows nslices \n"
"\n"
"Set output voxel dimension. Overrides voxres in template.\n"
"\n"
"--xfm xfmfile\n"
"\n"
"Apply matrix to input XYZ to get output XYZ. Matrix is assumed to be \n"
"tkregister-style. Note: the output of tkregister is a matrix that \n"
"maps XYZCor to XYZFunc, so, to resample a functional into the space\n"
"of the structural, make sure to invert the transform with -invxfm.\n"
"If no xfm is specified, the identity is assumed.\n"
"\n"
"--invxfm\n"
"\n"
"Invert the xfm matrix before applying transform.\n"
"\n"
"--fstal\n"
"\n"
"Resample the input volume to talairach space. The xfm file must be\n"
"a register.dat-format matrix. The talairach matrix is obtained from\n"
"talairach.xfm from SUBJECTS_DIR/subjid/transforms. SUBJECTS_DIR is\n"
"read from the environment or can be specified with --sd. subjid is\n"
"read from the xfm file. The transformation matrix is then computed\n"
"as Xtal*inv(R), where Xtal is talairach.xfm matrix and R is the matrix\n"
"in the xfm file. Don't use --invxfm unless you want to go from \n"
"talairach space back to the input.\n"
"\n"
"--irescale min max\n"
"\n"
"Rescale intensity to be between min and max. This can be useful when\n"
"the output precision is less than that of the input. When the output\n"
"is COR, the values are automatically rescaled to 0 to 255.\n"
"\n"
"--interp method\n"
"\n"
"Interpolate the output based on the given method. Legal values\n"
"are: nearest, trilin, and sinc. trilin is the default. sinc requires\n"
"one parameter (hw). sinc probably does not work.\n"
"\n"
"--gdiagno diagnostic level\n"
"\n"
"Sets the diagnostic level (only good for debuggin').\n"
"\n"
"--version\n"
"\n"
"Print out version string and exit.\n"
"\n"
"--help \n"
"\n"
"Prints out all this information.\n"
"\n"
"ALGORITH/TKREGISTER MATRIX CONENTION 

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

1. Resample a functional volume to talairach space.

   mri_vol2vol --in f_000.bfloat --out ftal_000.bfloat
     --temp $SUBJECTS_DIR/talairach/mri/orig 
     --xfm register.dat --fstal

2. Resample a structural volume to talairach space.

   cd  $SUBJECTS_DIR/subjid/mri/
   mkdir -p orig-tal
   mri_vol2vol --in orig --out orig-tal
               --temp $SUBJECTS_DIR/talairach 
               --xfm transforms/talairach.xfm 

   NOTE: this should give the same result as:
   mri_convert orig orig-tal --apply_transform transforms/talairach.xfm 

3. Resample a subcortical segmentation to functional space. It uses
   nearest-neighbor interp because the segmentation values are
   categorical, not continuous. 

   mri_vol2vol --in   $SUBJECTS_DIR/subjid/mri/aseg 
               --out  aseg_000.bshort
               --temp func_000.bshort
               --xfm  register.dat
               --interp nearest
"
"\n"
"FORMATS\n"
"\n"
"Data file format can be specified implicitly (through the path name)\n"
"or explicitly. All formats accepted by mri_convert can be used. \n"
"\n"
"BUGS\n"
"\n"
"sinc interpolation is broken except for maybe COR to COR.\n"
"\n"
"BUG REPORTING\n"
"\n"
"Report bugs to analysis-bugs@nmr.mgh.harvard.edu. Include the following \n"
"formatted as a list as follows: (1) command-line, (2) directory where\n"
"the program was run (for those in the MGH-NMR Center), (3) version, \n"
"(4) text output, (5) description of the problem.\n"
"\n"
"SEE ALSO \n"
"\n"
"mri_convert\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(involpath == NULL){
    printf("ERROR: No input supplied.\n");
    exit(1);
  }
  if(outvolpath == NULL){
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

/*---------------------------------------------------------------
  MRIvol2Vol() - samples the values of one volume into that of
  another. Handles multiple frames. Can do nearest-neighbor,
  trilinear, and sinc interpolation (sinc may be flaky).

  Vt2s is the 4x4 matrix which converts CRS in the target to
  CRS in the source.

  InterpCode is either: SAMPLE_NEAREST, SAMPLE_TRILINEAR, or
  SAMPLE_SINC.

  If either firstframe or lastframe are less than zero, then all
  frames are resampled. Regarless of the values of first and last
  frame, the target volume must have the same number of frames as the
  source volume. If a subset of frames are resampled, the only the
  corresponding output frames are affected. Ie, if only the 3rd frame
  is resampled, the output will go into the 3rd frame of the target
  volume (not the first frame). The first and last frames are
  specifiable for applications that maybe time sensative.

  param is a generic parameter. For sinc, param is the hw parameter,
  otherwise, it currently has no meaning.
  ---------------------------------------------------------------*/
int MRIvol2Vol(MRI *src, MRI *targ, MATRIX *Vt2s, int InterpCode, 
	       int firstframe, int lastframe, float param)
{
  int   ct,  rt,  st,  f;
  int   ics, irs, iss;
  float fcs, frs, fss;
  float *valvect;
  int sinchw;
  Real rval;

  if(src->nframes != targ->nframes){
    printf("ERROR: MRIvol2vol: source and target have different number "
	   "of frames\n");
    return(1);
  }

  if(firstframe < 0 || lastframe < 0){
    firstframe = 0;
    lastframe = src->nframes-1;
  }
  if(firstframe >= src->nframes){
    printf("ERROR: MRIvol2vol: fristframe = %d exceeds number of frames\n",
	   firstframe);
    return(1);
  }
  if(lastframe >= src->nframes){
    printf("ERROR: MRIvol2vol: lastframe = %d exceeds number of frames\n",
	   lastframe);
    return(1);
  }

  sinchw = nint(param);

  valvect = (float *) calloc(sizeof(float),src->nframes);

  for(ct=0; ct < targ->width; ct++){
    for(rt=0; rt < targ->height; rt++){
      for(st=0; st < targ->depth; st++){
	
	/* Column in source corresponding to CRS in Target */
	fcs = Vt2s->rptr[1][1] * ct + Vt2s->rptr[1][2] * rt +
	      Vt2s->rptr[1][3] * st + Vt2s->rptr[1][4] ;
	ics = nint(fcs);
	if(ics < 0 || ics >= src->width) continue;

	/* Row in source corresponding to CRS in Target */
	frs = Vt2s->rptr[2][1] * ct + Vt2s->rptr[2][2] * rt +
	      Vt2s->rptr[2][3] * st + Vt2s->rptr[2][4] ;
	irs = nint(frs);
	if(irs < 0 || irs >= src->height) continue;

	/* Slice in source corresponding to CRS in Target */
	fss = Vt2s->rptr[3][1] * ct + Vt2s->rptr[3][2] * rt +
	      Vt2s->rptr[3][3] * st + Vt2s->rptr[3][4] ;
	iss = nint(fss);
	if(iss < 0 || iss >= src->depth) continue;

	/* Assign output volume values */
	if(InterpCode == SAMPLE_TRILINEAR)
	  MRIsampleSeqVolume(src, fcs, frs, fss, valvect, 
			     firstframe, lastframe) ;
	else{
	  for(f=firstframe; f <= lastframe ; f++){
	    switch(InterpCode){
	    case SAMPLE_NEAREST:
	      valvect[f] = MRIgetVoxVal(src,ics,irs,iss,f);
	      break ;
	    case SAMPLE_SINC:      /* no multi-frame */
	      MRIsincSampleVolume(src, fcs, frs, fss, sinchw, &rval) ;
	      valvect[f] = rval;
	      break ;
	    }
	  }
	}

	for(f=0; f < src->nframes; f++)
	  MRIsetVoxVal(targ,ct,rt,st,f,valvect[f]);

      } /* target col */
    } /* target row */
  } /* target slice */

  free(valvect);

  return(0);
}


