/**
 * @file  mris_calc.c
 * @brief A *very* simple calculator for FreeSurfer curvature files.
 *
 * 'mris_calc' performs some simple mathematical operations on
 * FreeSurfer 'curv' format files. These files typically contain
 * curvature information, but are also used as more generic
 * placeholders for other data, such as thickness information.
 *
 */
/*
 * Original Author: Rudolph Pienaar
 * CVS Revision Info:
 *    $Author: rudolph $
 *    $Date: 2008/11/13 21:29:16 $
 *    $Revision: 1.9 $
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
#define _ISOC99_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <assert.h>
#include <errno.h>

#include <getopt.h>
#include <stdarg.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "fio.h"
#include "xDebug.h"
#include "mri_identify.h"

#define  STRBUF   65536
#define  MAX_FILES    1000
#define  CO( x )    fprintf(stdout, ( x ))
#define  CE( x )    fprintf(stderr, ( x ))
#define  START_i    3

static const char vcid[] =
"$Id: mris_calc.c,v 1.9 2008/11/13 21:29:16 rudolph Exp $";

// ----------------------------------------------------------------------------
// DECLARATION
// ----------------------------------------------------------------------------

// ------------------
// typedefs and enums
// ------------------

typedef enum _FILETYPE {
	e_Unknown, e_VolumeFile, e_CurvatureFile, e_SurfaceFile, e_FloatArg
} e_FILETYPE;

char* 	Gppch_filetype[] = {
	"Unknown",
	"Volume",
	"Curvature",
	"Surface",
	"FloatArg"
};

char*	Gppch_fileExt[] = {
	".null",
	".mgz",
	".crv",
	".surf",
	".dummy"
};

typedef enum _FILEACCESS {
	e_UNSPECIFIED		= -10,
	e_WRONGMAGICNUMBER	= -1, 
	e_OK			=  0, 
	e_READACCESSERROR	=  1, 
	e_WRITEACCESSERROR	=  2
} e_FILEACCESS;

typedef enum _operation {
  e_mul,
  e_div,
  e_add,
  e_sub,
  e_set,
  e_abs
} e_operation;

const char* Gppch_operation[] = {
  "multiply",
  "divide",
  "add",
  "subtract",
  "set",
  "abs",
  "min",
  "max",
  "size",
  "mini",
  "maxi",
  "mean",
  "std",
  "stats",
  "ascii"
};


// -------------------------------
// Global "class" member variables
//--------------------------------

char*           	G_pch_progname ;
char*           	Progname ;

static int      	G_verbosity             = 0;
static FILE*    	G_FP                    = NULL;

// Input 1
static int      	G_sizeCurv1             = 0;
static char*    	G_pch_curvFile1         = NULL;
static float*   	G_pf_arrayCurv1         = NULL;
static int      	G_nfaces                = 0;
static int      	G_valsPerVertex         = 0;
static e_FILETYPE	G_eFILETYPE1		= e_Unknown;

// Input 2
static int      	Gb_curvFile2            = 0;	//  The second input 
static char*    	G_pch_curvFile2         = NULL; //+ file is optional.
static int      	G_sizeCurv2             = 0;
static float*   	G_pf_arrayCurv2         = NULL;
static e_FILETYPE	G_eFILETYPE2		= e_Unknown;

// "Helper" pointers
static MRI*		Gp_MRI			= NULL; // Pointer to most
							//+ recently read 
							//+ MRI_VOLUME struct
							//+ and used in volume
							//+ handling... mostly
							//+ for volume size.

// Operation to perform on input1 and input2
static char*    	G_pch_operator          = NULL;

// Output file
static int      	G_sizeCurv3             = 0;
static char    		G_pch_curvFile3[STRBUF];
static float*   	G_pf_arrayCurv3         = NULL;
static short		Gb_file3		= 0;

//----------------
// "class" methods
//----------------

// Housekeeping functions
static void synopsis_show(void);
static void version_print(void);

// Setup functions
static int  options_parse(
  int argc,
  char* apchv[]
  );
static int  options_print(void);

// Simple functions on two float arguments
double fn_mul(float af_A, float af_B)   {return (af_A * af_B);}
double fn_div(float af_A, float af_B)   {return (af_B != 0 ? (af_A / af_B) : 0.0);}
double fn_add(float af_A, float af_B)   {return (af_A + af_B);}
double fn_sub(float af_A, float af_B)   {return (af_A - af_B);}
double fn_set(float af_A, float af_B)   {return (af_B);}

// Simple functions on one argument
double fn_abs(float af_A)		{return fabs(af_A);}

double fn_min(float af_A) {
    static float        f_min  = 0;
    static int          count  = 0;
    if(!count) f_min = af_A;
    if(af_A <= f_min)
        f_min = af_A;
    count++;
    return f_min;
}
double fn_mini(float af_A) {
    static float        f_min  = 0;
    static int          mini   = -1;
    static int          count  = 0;
    if(!count) f_min = af_A;
    if(af_A <= f_min) {
        f_min   = af_A;
        mini    = count;
    }
    count++;
    return (float) mini;
}
double fn_max(float af_A) {
    static float        f_max  = 0;
    static int          count  = 0;
    if(!count) f_max = af_A;
    if(af_A >= f_max)
        f_max = af_A;
    count++;
    return f_max;
}
double fn_maxi(float af_A) {
    static float        f_max  = 0;
    static int          maxi   = -1;
    static int          count  = 0;
    if(!count) f_max = af_A;
    if(af_A >= f_max) {
        f_max   = af_A;
        maxi    = count;
    }
    count++;
    return maxi;
}

//  The following two functions are identical. When the "stats" operator
//+ is called, the fn_sum() cannot be re-used since the static f_sum
//+ will corrupt the recalculation of the sum as determined by a call
//+ to fn_mean() 
double fn_sum(float af_A) {
    static double f_sum  = 0.;
    f_sum += af_A;
    return f_sum;
}

double fn2_sum(float af_A) {
    static double f_sum  = 0.;
    f_sum += af_A;
    return f_sum;
}

double fn_sum2(float af_A) {
    static double f_sum2 = 0.;
    f_sum2  += af_A * af_A;
    return f_sum2;
}

double fn_mean(float af_A) {
    static int  count   = 0;
    float       f_sum   = 0.;
    float       f_mean  = 0.;
    f_sum       = fn_sum(af_A);
    f_mean      = f_sum / ++count;
    return(f_mean);
}

double fn_dev(float af_A) {
    double      f_sum   = 0.0;
    double      f_sum2  = 0.0;
    static int  count   = 1;
    double      f_dev   = 0.;
    f_sum       = fn2_sum(af_A);
    f_sum2      = fn_sum2(af_A);
    f_dev       = (count*f_sum2 - f_sum*f_sum)/count;
    count++;
    return f_dev;
}

// Info functions
e_FILETYPE
fileType_find(
  char*		apch_inputFile
  );

// I/O functions
short CURV_arrayProgress_print(
  int   	asize,
  int   	acurrent,
  char* 	apch_message
  );

e_FILEACCESS
VOL_fileRead(
  char* 	apch_VolFileName,
  int*  	ap_vectorSize,
  float*  	apf_volData[]
  );

e_FILEACCESS
VOL_fileWrite(
  char* 	apch_VolFileName,
  int  		a_vectorSize,
  float*  	apf_volData
  );

e_FILEACCESS
CURV_fileRead(
  char* 	apch_curvFileName,
  int*  	ap_vectorSize,
  float*  	apf_curv[]
  );

e_FILEACCESS
CURV_fileWrite(
  char* 	apch_curvFileName,
  int  		a_vectorSize,
  float*  	apf_curv
  );

short   b_outCurvFile_write(char* apch_operator);
short   CURV_process(void);
short   CURV_functionRunABC( double (*F)(float f_A, float f_B) );
double  CURV_functionRunAC( double (*F)(float f_A) );

int main(int argc, char *argv[]) ;

// ----------------------------------------------------------------------------
// IMPLEMENTATION
// ----------------------------------------------------------------------------

static void
synopsis_show(void) {
  char  pch_synopsis[STRBUF];

  sprintf(pch_synopsis, "    \n\
    NAME \n\
 \n\
          mris_calc \n\
 \n\
    SYNOPSIS \n\
 \n\
          mris_calc [OPTIONS] <file1> <ACTION> [<file2> | <floatNumber>] \n\
 \n\
    DESCRIPTION \n\
 \n\
	'mris_calc' is a simple calculator that operates on FreeSurfer \n\
	curvatures and volumes. \n\
 \n\
	In most cases, the calculator functions with three arguments: \n\
        two inputs and an <ACTION> linking them. Some actions, however, \n\
        operate with only one input <file1>. \n\
 \n\
        In all cases, the first input <file1> is the name of a FreeSurfer \n\
        curvature overlay (e.g. rh.curv) or volume file (e.g. orig.mgz). \n\
 \n\
        For two inputs, the calculator first assumes that the second input \n\
        is a file. If, however, this second input file doesn't exist, the \n\
        calculator assumes it refers to a float number, which is then \n\
        processed according to <ACTION>. \n\
 \n\
    OPTIONS \n\
 \n\
    	--output <outputCurvFile> \n\
   	 -o <outputCurvFile> \n\
 \n\
    	By default, 'mris_calc' will save the output of the calculation to a \n\
    	file in the current working directory with filestem 'out'. The file \n\
    	extension is automatically set to the appropriate filetype based on \n\
    	the input. For any volume type, the output defaults to '.mgz' and for \n\
    	curvature inputs, the output defaults to '.crv'. \n\
 \n\
    	--version \n\
    	-v \n\
 \n\
	Print out version number. \n\
 \n\
    	--verbosity <value> \n\
 \n\
	Set the verbosity of the program. Any positive value will trigger \n\
	verbose output, displaying intermediate results. The <value> can be \n\
	set arbitrarily. Useful mostly for debugging. \n\
 \n\
    ACTION \n\
 \n\
	The action to be perfomed on the two curvature files. This is a \n\
	text string that defines the mathematical operation to execute. For two \n\
	inputs, this action is applied in an indexed element-by-element fashion, \n\
	i.e. <file1>[n] <ACTION> <file2>[n] where 'n' is an index counter into \n\
	the data space. \n\
 \n\
	ACTION  INPUTS OUTPUTS	                EFFECT \n\
	  mul	   2      1     <outputFile> = <file1> * <file2> \n\
	  div	   2	  1     <outputFile> = <file1> / <file2> \n\
	  add      2	  1     <outputFile> = <file1> + <file2> \n\
	  sub	   2	  1     <outputFile> = <file1> - <file2> \n\
	  set	   2	  1     <outputFile> = <file1> \n\
 \n\
          ascii    1      1     <outputFile> = ascii <file1> \n\
	  abs	   1      1     <outputFile> = abs(<file1>) \n\
 \n\
          size     1      0     print the size (number of elements) of <file1> \n\
          min      1      0     print the min value (and index) of <file1> \n\
          max      1      0     print the max value (and index) of <file1> \n\
          mean     1      0     print the mean value of <file1> \n\
          std      1      0     print the standard deviation of <file1> \n\
          stats    1      0     process 'size', 'min', 'max', 'mean', 'std' \n\
 \n\
	The 'set' command is somewhat different in that for practical purposes \n\
	the contents of <file1> are ignored. It is still important to \n\
	specifiy a valid <file1> since it is parsed by 'mris_calc' \n\
	to determine the size and output filetype file to create. In most \n\
	instances, <file2> will denote a float value, and not an actual \n\
	curvature file, i.e. 'mris_calc rh.area set 0.005' will create \n\
	an output curvature, 'out.crv' of the same size as rh.area, and with \n\
	each element set to 0.005. Similarly for volumes, \n\
	'mris_calc orig.mgz set 127' will create a volume file, out.mgz with \n\
	all voxel intensities set to 127. \n\
 \n\
        The 'ascii' command converts <file1> to a text format file, \n\
        suitable for reading into MatLAB, for example. Note that volumes \n\
        are written out as a 1D linear array. \n\
 \n\
        Note also that the standard deviation can suffer from float rounding \n\
        errors and is only accurate to 4 digits of precision. \n\
 \n\
    ARBITRARY FLOATS AS SECOND INPUT ARGUMENT \n\
 \n\
	If a second input argument is specified, 'mris_calc' will attempt to \n\
        open the argument following <ACTION> as if it were a curvature file. \n\
        Should this file not exist, 'mris_calc' will attempt to parse the \n\
        argument as if it were a float value. \n\
 \n\
	In such a case, 'mris_calc' will create a dummy internal \n\
	array structure and set all its elements to this float value. \n\
 \n\
    NOTES \n\
 \n\
	<file1> and <file2> should typically be generated on the \n\
	same subject. \n\
 \n\
    EXAMPLES \n\
 \n\
    	$>mris_calc rh.area mul rh.thickness \n\
 \n\
	Multiply each value in <rh.area> with the corresponding value \n\
	in <rh.thickness>, creating a new file called 'out.crv' that \n\
	contains the result. \n\
 \n\
    	$>mris_calc --output rh.weightedCortex rh.area mul rh.thickness \n\
 \n\
	Same as above, but give the ouput file the more meaningful name \n\
	of 'rh.weightedCortex'. \n\
 \n\
        $>mris_calc rh.area max \n\
 \n\
        Determine the maximum value in 'rh.area' and print to stdout. In \n\
        addition to the max value, the index offset in 'rh.area' containing \n\
        this value is also printed. \n\
 \n\
        $>mris_calc rh.area stats \n\
 \n\
        Determine the size, min, max, mean, and std of 'rh.area'. \n\
 \n\
	$>mris_calc orig.mgz sub brainmask.mgz \n\
 \n\
	Subtract the brainmask.mgz volume from the orig.mgz volume. Result is \n\
	saved by default to out.mgz. \n\
 \n\
    ADVANCED EXAMPLES \n\
 \n\
	Consider the case when calculating the right hemisphere pseudo volume \n\
	formed by the FreeSurfer generated white matter 'rh.area' curvature \n\
	file, and the cortical thickness, 'rh.thickness'. Imagine this is to \n\
	be expressed as a percentage of intercranial volume. \n\
 \n\
	$>mris_calc -o rh.cortexVol rh.area mul rh.thickness \n\
	Calculate the volume and store in a curvature format: \n\
 \n\
	Now, find the intercranial volume (ICV) in the corresponding output \n\
	file generated by FreeSurfer for this subject. Assume ICV = 100000. \n\
 \n\
	$>mris_calc -o rh.cortexVolICV rh.cortexVol div 100000 \n\
	Here the second <ACTION> argument is a number and not a curvature file. \n\
 \n\
	We could have achieved the same effect by first creating an \n\
	intermediate curvature file, 'rh.ICV' with each element set to \n\
	the ICV, and then divided by this curvature: \n\
 \n\
	$>mris_calc -o rh.ICV rh.area set 100000 \n\
	$>mris_calc -o rh.cortexVolICV rh.cortexVol div rh.ICV \n\
 \n\
 \n\
\n");

  fprintf(stdout,pch_synopsis);
  exit(1) ;
}

void
simpleSynopsis_show(void) {
  char  pch_errorMessage[STRBUF];

  sprintf(pch_errorMessage, "Insufficient number of arguments.");
  sprintf(pch_errorMessage,
          "%s\nYou should specify '<input1> <ACTION> [<input2> | <floatNumber>]'",
          pch_errorMessage);
  sprintf(pch_errorMessage,
          "%s\nUse a '-u' for full usage instructions.",
          pch_errorMessage);
  ErrorExit(10, "%s: %s", G_pch_progname, pch_errorMessage);
}

void *
xmalloc (size_t size)
{
  register void *value = malloc (size);
  if (value == 0)
    ErrorExit(10, "%s: virtual memory exhausted.", G_pch_progname);
  return value;
}

void
verbosity_set(void) {
  if(G_verbosity) {
    G_FP  = stdout;
    options_print();
  } else
    if((G_FP = fopen("/dev/null", "w")) == NULL)
      ErrorExit(  ERROR_NOFILE,
                  "%s: Could not open /dev/null for console sink!",
                  G_pch_progname);
}

void
init(void) {
  InitDebugging( "mris_calc" );
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
}

void 
error_exit(
	char* 	apch_action, 
	char* 	apch_error, 
	int 	exitCode) {

  char	pch_errorMessage[STRBUF];
  strcpy(pch_errorMessage, "");

  sprintf(pch_errorMessage, "\n%s:", G_pch_progname);
  sprintf(pch_errorMessage, 
	"%s\n\tSorry, but I seem to have encountered an error.",
	pch_errorMessage);
  sprintf(pch_errorMessage, 
	"%s\n\tWhile %s,", pch_errorMessage, apch_action);
  sprintf(pch_errorMessage,
	"%s\n\t%s\n", pch_errorMessage, apch_error);

  fprintf(stderr, pch_errorMessage);
  exit(exitCode);

}

void
error_unmatchedSizes(void) {
  error_exit(	"checking on input filetype sizes",
		"I found a size mismatch, i.e. len(input1)!=len(input2)",
		10);
}

void
error_noVolumeStruct(void) {
  error_exit(	"checking on output volume",
		"it seems that internal volume structure is invalid.",
		10);
}

void
error_volumeWriteSizeMismatch(void) {
  error_exit(	"checking on output volume",
		"I found a size mismatch, internal data does not fit into vol.",
		10);
}

void
error_incompatibleFileTypes(void) {
  error_exit(	"checking on input filetypes",
		"it seems that you specified incompatible types.",
		10);
}

void
error_zeroLengthInput(void) {
  error_exit(	"checking input 1",
		"it seems that the input file has no contents.",
		10);
}

void
error_surfacesNotHandled(void) {
  error_exit(	"checking inputs",
		"FreeSurfer surface files are not handled yet.",
		10);
}

void
output_init(void) {
  int i;
  G_sizeCurv3   = G_sizeCurv1;
  G_pf_arrayCurv3 = (float*) xmalloc(G_sizeCurv1 * sizeof(float));
  for(i=0; i<G_sizeCurv3; i++)
    G_pf_arrayCurv3[i] = 0.0;
  if(!Gb_file3) strcat(G_pch_curvFile3, Gppch_fileExt[G_eFILETYPE1]);
}

e_FILETYPE
fileType_find(
  char*		apch_inputFile) {

  int		type, len;
  float		f			= -1.0;
  char**	ppch_end		= &apch_inputFile;

  // First, check if we have a valid float arg conversion
  errno = 0;
  f 	= strtof(apch_inputFile, ppch_end);
  len	= (int) strlen(apch_inputFile);
  if(!len) {
	return(e_FloatArg);
  }

  // Check if input is a volume file...
  type = mri_identify(apch_inputFile);
  if(type != MRI_VOLUME_TYPE_UNKNOWN && type != MRI_CURV_FILE) 
	return(e_VolumeFile);
  
  // Check if input is a curvature file...
  if(type == MRI_CURV_FILE) return(e_CurvatureFile);

  // Assume that the input is a surface file...
  return(e_SurfaceFile);
}

void
fileIO_errorHander(
  char* 	pch_filename, 
  e_FILEACCESS eERROR) {
  switch(eERROR) {
	case e_UNSPECIFIED:
	break;
	case e_OK:
	break;
	case e_READACCESSERROR:
    	  ErrorExit(ERROR_BADPARM,
              "\n%s: could not establish read access to '%s'.\n",
              G_pch_progname, pch_filename);  
	break;
	case e_WRONGMAGICNUMBER:
    	  ErrorExit(ERROR_BADPARM,
              "\n%s: curvature file '%s' has wrong magic number.\n",
              G_pch_progname, pch_filename);  
	break;
	case e_WRITEACCESSERROR:
    	  ErrorExit(ERROR_BADPARM,
              "\n%s: could not establish write access to '%s'.\n",
              G_pch_progname, pch_filename);  
	break;
  }
}

e_FILEACCESS
fileRead(
  char* 	apch_fileName,
  int*  	ap_vectorSize,
  float*  	apf_curv[],
  e_FILETYPE*	aeFILETYPE) {

  e_FILEACCESS	eACCESS		= e_UNSPECIFIED;

  *aeFILETYPE	= fileType_find(apch_fileName);
  switch(*aeFILETYPE) {
	case e_Unknown:
	break;
	case e_FloatArg:
	break;
	case e_VolumeFile:
	eACCESS	= VOL_fileRead( apch_fileName, ap_vectorSize, apf_curv);
	break;
	case e_CurvatureFile:
	eACCESS	= CURV_fileRead(apch_fileName, ap_vectorSize, apf_curv);
	break;
	case e_SurfaceFile:
	error_surfacesNotHandled();
	break;
  }
  return eACCESS;
}

e_FILEACCESS
fileWrite(
  char* 	apch_fileName,
  int  		a_vectorSize,
  float*  	apf_curv
) {

  e_FILEACCESS	eACCESS		= e_UNSPECIFIED;

  switch(G_eFILETYPE1) {
	case e_Unknown:
	break;
	case e_FloatArg:
	break;
	case e_VolumeFile:
	eACCESS	= VOL_fileWrite( apch_fileName, a_vectorSize, apf_curv);
	break;
	case e_CurvatureFile:
	eACCESS	= CURV_fileWrite(apch_fileName, a_vectorSize, apf_curv);
	break;
	case e_SurfaceFile:
	error_surfacesNotHandled();
	break;
  }
  return eACCESS;
}

void
debuggingInfo_display() {
  cprintd("Size of input1: ", 	G_sizeCurv1);
  cprints("Type of input1: ", 	Gppch_filetype[G_eFILETYPE1]);
  cprints("ACTION:", 		G_pch_operator);
  cprintd("Size of input2: ", 	G_sizeCurv2);
  cprints("Type of input2: ", 	Gppch_filetype[G_eFILETYPE2]);
}

int
main(
  int   argc,
  char  *argv[]
  ) {

  int   	nargs, i;
  short 	ret		= 0;
  float		f_curv2;
  e_FILEACCESS	eACCESS;

  init();
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_calc.c,v 1.9 2008/11/13 21:29:16 rudolph Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  G_pch_progname = argv[0] ;
  strcpy(G_pch_curvFile3, "out");

  for (; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = options_parse(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  // Process command line surface options and operator
  if((argc != 4) && (argc != 3)) simpleSynopsis_show();
  if(argc==4) Gb_curvFile2              = 1;
  G_pch_curvFile1                       = argv[1];
  G_pch_operator                        = argv[2];
  if(Gb_curvFile2) G_pch_curvFile2      = argv[3];
  verbosity_set();

  // Read in relevant input files -
  eACCESS	= fileRead(G_pch_curvFile1, 
				&G_sizeCurv1, 
				&G_pf_arrayCurv1, 
				&G_eFILETYPE1);  
  if(eACCESS != e_OK) fileIO_errorHander(G_pch_curvFile1, eACCESS);

  // Second input file is optional, and, if specified could in fact
  //+ denote a float value and not an actual file per se.
  if(Gb_curvFile2) {
      eACCESS 	= fileRead(G_pch_curvFile2, 
				&G_sizeCurv2, 
				&G_pf_arrayCurv2,
				&G_eFILETYPE2);
      if(G_eFILETYPE2 == e_FloatArg) {
    	f_curv2   	= atof(G_pch_curvFile2);
    	G_sizeCurv2 	= G_sizeCurv1;
    	G_pf_arrayCurv2 = (float*) malloc(G_sizeCurv1 * sizeof(float));
    	for(i=0; i<G_sizeCurv2; i++)
      	    G_pf_arrayCurv2[i] = f_curv2;
  	}
  }
  if(G_verbosity) debuggingInfo_display();
  if(!G_sizeCurv1)
	error_zeroLengthInput();
  if(G_sizeCurv1 != G_sizeCurv2 && argc==4) 
	error_unmatchedSizes();
  if(G_eFILETYPE1 != G_eFILETYPE2 && G_eFILETYPE2 != e_FloatArg && argc==4)
	error_incompatibleFileTypes();

  output_init();
  ret     = CURV_process();
  if(G_verbosity) printf("\n");
  return(0);
}

static int
options_print(void) {
  cprints("Input file 1:", G_pch_curvFile1);
  cprints("ACTION:", G_pch_operator);
  if(G_pch_curvFile2) {
	cprints("Input file 2:", G_pch_curvFile2);
  }
  cprints("Output file:", G_pch_curvFile3);

  return 1;
}

static int
options_parse(int argc, char *argv[]) {
  int    nargs    = 0;
  char*  option;
  char*  pch_text;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-output") || (toupper(*option) == 'O')) {
    strcpy(G_pch_curvFile3, (argv[2]));
    Gb_file3		= 1;
    nargs     		= 1;
  } else if (!stricmp(option, "-version") || (toupper(*option) == 'V')) {
    version_print();
  } else if (!stricmp(option, "-verbosity")) {
    G_verbosity		= atoi(argv[2]);
    nargs     		= 1;
  } else switch (toupper(*option)) {
  case 'T':
    pch_text  = "void";
    break ;
  case '?':
  case 'U':
    synopsis_show() ;
    exit(1) ;
    break ;
  case '-':
    break;
  default:
    fprintf
      (stderr, 
       "Unknown option '%s'. Looking for help? Use '-u' instead.\n", 
       argv[1]) ;
    exit(1) ;
    break ;
  }
  return(nargs) ;
}

static void
version_print(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


/*!
  \fn VOL_fileRead(char* apch_volFileName, int* ap_vectorSize, float* apf_data)
  \brief Read a FreeSurfer volume file into a float array
  \param apch_curvFileName The name of a FreeSurfer volume file.
  \param ap_vectorSize Pointer to the size (i.e. number of elements) in volume.
  \param apf_curv 1D Array containing the volume values.
  \return If volume file is successfully read, return e_OK. If file could not be opened, return e_READACCESSERROR.
*/
e_FILEACCESS
VOL_fileRead(
  char* 	apch_volFileName,
  int*  	ap_vectorSize,
  float*  	apf_data[]
  ) {

  char  	pch_readMessage[STRBUF];
  int		i, j, k;
  int		I 				= 0; 
  MRI*		pMRI				= NULL;
  float*  	pf_data				= NULL;

  if(G_verbosity) cprints("Reading...", "");
  if( (pMRI = MRIread(apch_volFileName)) == NULL)
    return e_READACCESSERROR;
  if(G_verbosity) cprints("", "ok");
  Gp_MRI		= pMRI;		// Global pointer.
  *ap_vectorSize	= pMRI->width * pMRI->height * pMRI->depth;
  pf_data   		= (float*) xmalloc(*ap_vectorSize * sizeof(float));
  sprintf(pch_readMessage, "Packing %s", apch_volFileName);
  for(i=0; i<pMRI->width; i++)		// 'x', i.e. columns in slice
    for(j=0; j<pMRI->height; j++)	// 'y', i.e. rows in slice
      for(k=0; k<pMRI->depth; k++) {	// 'z', i.e. # of slices
	CURV_arrayProgress_print(*ap_vectorSize, I, pch_readMessage);
	pf_data[I++]	= MRIgetVoxVal(pMRI, i, j, k, 0);
      }
  *apf_data = pf_data;
  return(e_OK);
}

/*!
  \fn VOL_fileWrite(char* apch_volFileName, int a_vectorSize, float* apf_data)
  \brief Write the 1D array data to a FreeSurfer volume
  \param apch_volFileName The name of a FreeSurfer volume file.
  \param a_vectorSize Size (i.e. number of elements) of array.
  \param apf_data Array containing the data values.
  \return If volume file is successfully written, return e_OK, else return e_WRITEACCESSERROR.
*/
e_FILEACCESS
VOL_fileWrite(
  char* 	apch_volFileName,
  int  		a_vectorSize,
  float*  	apf_data
) {
  //
  // PRECONDITIONS
  // o The Gp_MRI *must* have been set with a previous call to VOL_fileRead.
  //
  // POSTCONDITIONS
  // o Gp_MRI will have its voxel data set element-by-element to apf_data.
  // o Gp_MRI saved to <apch_volFileName>.
  //
  
  int		volSize;
  int		i, j, k;
  int		I				= 0;
  char  	pch_readMessage[STRBUF];

  if(!Gp_MRI) 			error_noVolumeStruct();
  volSize	= Gp_MRI->width * Gp_MRI->height * Gp_MRI->depth;
  if(volSize != a_vectorSize)	error_volumeWriteSizeMismatch();
  sprintf(pch_readMessage, "Packing %s", apch_volFileName);
  for(i=0; i<Gp_MRI->width; i++)		// 'x', i.e. columns in slice
    for(j=0; j<Gp_MRI->height; j++)		// 'y', i.e. rows in slice
      for(k=0; k<Gp_MRI->depth; k++) {		// 'z', i.e. # of slices
	CURV_arrayProgress_print(a_vectorSize, I, pch_readMessage);
	MRIsetVoxVal(Gp_MRI, i, j, k, 0, apf_data[I++]);
      }
  if(G_verbosity) cprints("Saving", "");
  return(MRIwrite(Gp_MRI, apch_volFileName));
  if(G_verbosity) cprints("", "ok");
}

/*!
  \fn CURV_fileRead(char* apch_curvFileName, int* ap_vectorSize, float* apf_curv)
  \brief Read a FreeSurfer curvature file into a float array
  \param apch_curvFileName The name of a FreeSurfer curvature file.
  \param ap_vectorSize Pointer to the size (i.e. number of elements) in file.
  \param apf_curv Array containing the curvature values.
  \return If curvature file is successfully read, return e_OK. If file could not be opened, return e_READACCESSERROR. If file not a curvature format file, return e_WRONGMAGICNUMBER.
*/
#define   NEW_VERSION_MAGIC_NUMBER  16777215
e_FILEACCESS
CURV_fileRead(
  char* apch_curvFileName,
  int*  ap_vectorSize,
  float*  apf_curv[]
  ) {

  FILE* 	FP_curv;
  int   	vnum;
  int   	nvertices;
  int   	i;
  char  	pch_readMessage[STRBUF];
  float*  	pf_data				= NULL;

  if((FP_curv = fopen(apch_curvFileName, "r")) == NULL) return(e_READACCESSERROR);
  fread3(&vnum, FP_curv);
  if(vnum == NEW_VERSION_MAGIC_NUMBER) {
    nvertices = freadInt(FP_curv);
    G_nfaces  = freadInt(FP_curv);
    G_valsPerVertex = freadInt(FP_curv);
    *ap_vectorSize  = nvertices;
    pf_data   = (float*) xmalloc(nvertices * sizeof(float));
    sprintf(pch_readMessage, "Reading %s", apch_curvFileName);
    for(i=0; i<nvertices; i++) {
      CURV_arrayProgress_print(nvertices, i, pch_readMessage);
      pf_data[i]  = freadFloat(FP_curv);
    }
  } else return(e_WRONGMAGICNUMBER);
/*    ErrorExit(ERROR_BADPARM,
              "\n%s: curvature file '%s' has wrong magic number.\n",
              G_pch_progname, apch_curvFileName);*/
  *apf_curv   = pf_data;
  fclose(FP_curv);
  return(e_OK);
}

/*!
  \fn ascii_fileWrite(char* apch_fileName, float* apf_data)
  \brief Write internal data array to an ascii text file
  \param apch_fileName Output filename.
  \param apf_data Array data to output.
  \return If output file is successfully written, return e_OK, else return e_WRITEACCESSERROR.
*/
e_FILEACCESS
ascii_fileWrite(
  char* 	apch_fileName,
  float*  	apf_data
  ) {
  FILE* FP_curv;
  int   i;
  char  pch_readMessage[STRBUF];

  if((FP_curv = fopen(apch_fileName, "w")) == NULL)
    return(e_WRITEACCESSERROR);
  sprintf(pch_readMessage, "Writing %s", apch_fileName);
  for(i=0; i<G_sizeCurv1; i++) {
    CURV_arrayProgress_print(G_sizeCurv1, i, pch_readMessage);
    fprintf(FP_curv, "%f\n", apf_data[i]);
  }
  fclose(FP_curv);
  return(e_OK);
}


/*!
  \fn CURV_fileWrite(char* apch_curvFileName, int* ap_vectorSize, float* apf_curv)
  \brief Write a FreeSurfer curvature array to a file
  \param apch_curvFileName The name of a FreeSurfer curvature file.
  \param a_vectorSize Size (i.e. number of elements) of data array.
  \param apf_curv Array containing the curvature values.
  \return If curvature file is successfully written, return e_OK, else return e_WRITEACCESSERROR.
*/
e_FILEACCESS
CURV_fileWrite(
  char* 	apch_curvFileName,
  int  		a_vectorSize,
  float*  	apf_curv
  ) {
  FILE* FP_curv;
  int   i;
  char  pch_readMessage[STRBUF];

  if((FP_curv = fopen(apch_curvFileName, "w")) == NULL)
    return(e_WRITEACCESSERROR);
  fwrite3(NEW_VERSION_MAGIC_NUMBER, FP_curv);
  fwriteInt(a_vectorSize, FP_curv);
  fwriteInt(G_nfaces, FP_curv);
  fwriteInt(G_valsPerVertex, FP_curv);
  sprintf(pch_readMessage, "Writing %s", apch_curvFileName);
  for(i=0; i<a_vectorSize; i++) {
    CURV_arrayProgress_print(a_vectorSize, i, pch_readMessage);
    fwriteFloat(apf_curv[i], FP_curv);
  }
  fclose(FP_curv);
  return(e_OK);
}

short
CURV_arrayProgress_print(
  int   asize,
  int   acurrent,
  char* apch_message
  ) {
  //
  // PRECONDITIONS
  //  o <acurrent> is the current index being processed in a stream.
  //  o If <apch_message> is non-NULL, then prefix the progress bar
  //    with <apch_message> (and terminate progress bar with [ ok ]).
  //
  // POSTCONDITIONS
  //  o For every 5% of processed asize a "#" is written to G_FP
  //

  static int    fivePerc  = 0;
  fivePerc        = 0.05 * asize;

  if(!acurrent) {
    if(apch_message != NULL)
      fprintf(G_FP, "%*s", G_LC, apch_message);
    fprintf(G_FP, " [");
    fflush(G_FP);
  }
  if(acurrent%fivePerc == fivePerc-1) {
    fprintf(G_FP, "#");
    fflush(G_FP);
  }
  if(acurrent == asize-1) {
    fprintf(G_FP, "] ");
    if(apch_message != NULL)
      fprintf(G_FP, "%*s\n", 1, "[ ok ]");
  }
  return 1;
}

/*!
  \fn b_outCurvFile_write(char* apch_operator)
  \brief A simple function that determines if an output file should be saved.
  \param void
  \see
  \return Internal "class" global variables are set by this process.
*/
short
b_outCurvFile_write(char* apch_operator)
{
    // PRECONDITIONS
    //	pch_operator		in		Action to perform
    //
    // POSTCONDITIONS
    //	Internally an output curvature data structure is always generated.
    //	Saving this output file to disk is only meaningful for a given
    //	set of actions -- usually the typical mathematical operators.
    //	Other actions don't require an output file to be saved, e.g the
    //	'stats' action.
    //

    short	b_ret	= 0;

    if( (!strcmp(apch_operator, "mul"))		|| 	
    	(!strcmp(apch_operator, "div")) 	||	
    	(!strcmp(apch_operator, "add")) 	||	
    	(!strcmp(apch_operator, "sub")) 	||
    	(!strcmp(apch_operator, "set")) 	||
    	(!strcmp(apch_operator, "abs")) )
	b_ret = 1;

    return b_ret;	
}

/*!
  \fn CURV_process(void)
  \brief The main entry point for processing the curvature operations.
  \param void
  \see
  \return Internal "class" global variables are set by this process.
*/
short
CURV_process(void)
{
  // PRECONDITIONS
  //  o The following internal "class" variables are extant and valid:
  //    - G_pf_arrayCurv1, G_pf_arrayCurv2, G_pf_arrayCurv3
  //    - G_pch_operator
  //
  // POSTCONDITIONS
  //  o Depending on <G_pch_operator>, a simple calculation is performed
  //    to generate G_pf_arrayCurv3.
  //  o G_pf_arrayCurv3 is saved to G_pch_curvFile3
  //

  float f_min           = 0.;
  float f_max           = 0.;
  int   mini            = -1;
  int   maxi            = -1;
  int   b_canWrite      = 0;
  float f_mean  = 0.;
  float f_std   = 0.;
  float f_dev   = 0.;

  b_canWrite	= b_outCurvFile_write(G_pch_operator);

  if(     !strcmp(G_pch_operator, "mul")) {CURV_functionRunABC(fn_mul);}
  else if(!strcmp(G_pch_operator, "div")) {CURV_functionRunABC(fn_div);}
  else if(!strcmp(G_pch_operator, "add")) {CURV_functionRunABC(fn_add);}
  else if(!strcmp(G_pch_operator, "sub")) {CURV_functionRunABC(fn_sub);}
  else if(!strcmp(G_pch_operator, "set")) {CURV_functionRunABC(fn_set);}
  else if(!strcmp(G_pch_operator, "abs")) {CURV_functionRunAC( fn_abs);}

  if(!strcmp(G_pch_operator, "ascii")) {
    if(!Gb_file3) strcat(G_pch_curvFile3, ".ascii");
	printf("Write : %s", G_pch_curvFile3);

    if((ascii_fileWrite(G_pch_curvFile3, G_pf_arrayCurv1))==e_WRITEACCESSERROR)
	printf("Write error: %s", G_pch_curvFile3);
  }

  if(!strcmp(G_pch_operator, "size") || !strcmp(G_pch_operator, "stats")) {
    cprintd("Size", G_sizeCurv1);
  }
  if(!strcmp(G_pch_operator, "min") || !strcmp(G_pch_operator, "stats")) {
    f_min       = CURV_functionRunAC(fn_min);
    mini        = (int) CURV_functionRunAC(fn_mini);    
    cprintf("Min", f_min);
    cprintd("Index", mini);
  }
  if(!strcmp(G_pch_operator, "max") || !strcmp(G_pch_operator, "stats")) {
    f_max       = CURV_functionRunAC(fn_max);
    maxi        = (int) CURV_functionRunAC(fn_maxi);
    cprintf("Max", f_max);
    cprintd("Index", maxi);
  }
  if(!strcmp(G_pch_operator, "mean") || !strcmp(G_pch_operator, "stats")) {
    f_mean      = CURV_functionRunAC(fn_mean);
    cprintf("Mean", f_mean);
  }
  if(!strcmp(G_pch_operator, "std") || !strcmp(G_pch_operator, "stats")) {
    f_dev       = CURV_functionRunAC(fn_dev);
    f_std       = sqrt(f_dev/(G_sizeCurv1-1));
    cprintf("Std", f_std);
  }

  if(b_canWrite) fileWrite(G_pch_curvFile3, G_sizeCurv3, G_pf_arrayCurv3);

  return 1;
}

/*!
  \fn CURV_functionRunABC( (*F)(float f_A, float f_B) )
  \brief Loops over the internal curvature arrays and applies (*F) at each index
  \param (*F) A function of two floats that is applied at each curvature index.
  \see
  \return Internal "class" global variables are set by this process.
*/
short
CURV_functionRunABC( double (*F)(float f_A, float f_B) )
{
  // PRECONDITIONS
  //  o The following internal "class" variables are extant and valid:
  //    - G_pf_arrayCurv1, G_pf_arrayCurv2, G_pf_arrayCurv3
  //    - G_pch_operator
  //
  // POSTCONDITIONS
  //  o Depending on <G_pch_operator>, a simple calculation is performed
  //    to generate G_pf_arrayCurv3.
  //  o G_pf_arrayCurv3 is saved to G_pch_curvFile3
  //
  int   i;
  double f_a = 0.;
  double f_b = 0.;
  double f_c = 0.;

  for(i=0; i<G_sizeCurv1; i++) {
    f_a                 = G_pf_arrayCurv1[i];
    f_b                 = G_pf_arrayCurv2[i];
    f_c                 = (F)(f_a, f_b);
    G_pf_arrayCurv3[i]  = f_c;
  }
  return 1;
}

/*!
  \fn CURV_functionRunAC( (*F)(float f_A) )
  \brief Loops over the internal curvature arrays and applies (*F) at each index
  \param (*F) A function of two floats that is applied at each curvature index.
  \see
  \return Internal "class" global variables are set by this process.
*/
double
CURV_functionRunAC( double (*F)(float f_A) )
{
  // PRECONDITIONS
  //  o The following internal "class" variables are extant and valid:
  //    - G_pf_arrayCurv1, G_pf_arrayCurv3
  //    - G_pch_operator
  //
  // POSTCONDITIONS
  //  o Depending on <G_pch_operator>, a simple calculation is performed
  //    to generate G_pf_arrayCurv3.
  //  o G_pf_arrayCurv3 is saved to G_pch_curvFile3
  //
  int   i;
  double f_a = 0.;
  double f_c = 0.;

  for(i=0; i<G_sizeCurv1; i++) {
    f_a                 = G_pf_arrayCurv1[i];
    f_c                 = (F)(f_a);
    G_pf_arrayCurv3[i]  = f_c;
  }
  return f_c;
}
