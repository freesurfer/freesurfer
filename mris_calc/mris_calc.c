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
 *    $Date: 2008/09/24 20:34:54 $
 *    $Revision: 1.8 $
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

#define  STRBUF   65536
#define  MAX_FILES    1000
#define  CO( x )    fprintf(stdout, ( x ))
#define  CE( x )    fprintf(stderr, ( x ))
#define  START_i    3

static const char vcid[] =
"$Id: mris_calc.c,v 1.8 2008/09/24 20:34:54 rudolph Exp $";

// ----------------------------------------------------------------------------
// DECLARATION
// ----------------------------------------------------------------------------

// ------------------
// typedefs and enums
// ------------------
typedef enum _operation {
  e_mul,
  e_div,
  e_add,
  e_sub,
  e_set
} e_operation;

const char* Gppch_operation[] = {
  "multiply",
  "divide",
  "add",
  "subtract",
  "set",
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

char*           G_pch_progname ;
char*           Progname ;

static int      G_verbosity             = 0;
static FILE*    G_FP                    = NULL;

// Curvature 1
static int      G_sizeCurv1             = 0;
static char*    G_pch_curvFile1         = NULL;
static float*   G_pf_arrayCurv1         = NULL;
static int      G_nfaces                = 0;
static int      G_valsPerVertex         = 0;

// Curvature 2
static int      Gb_curvFile2            = 0;            //  The second curvature file
static char*    G_pch_curvFile2         = NULL;         //+ is optional.
static int      G_sizeCurv2             = 0;
static float*   G_pf_arrayCurv2         = NULL;

// Operation to perform on Curvature1 and Curvature2
static char*    G_pch_operator          = NULL;

// Output curvature file
static int      G_sizeCurv3             = 0;
static char*    G_pch_curvFile3         = "out.crv";
static float*   G_pf_arrayCurv3         = NULL;

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

// I/O functions
short CURV_arrayProgress_print(
  int   asize,
  int   acurrent,
  char* apch_message
  );

short   CURV_fileRead(
  char* apch_curvFileName,
  int*  ap_vectorSize,
  float*  apf_curv[]
  );

short   CURV_fileWrite(
  char* apch_curvFileName,
  int*  ap_vectorSize,
  float*  apf_curv
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

  sprintf(pch_synopsis, " \n\
    NAME \n\
 \n\
          mris_calc \n\
 \n\
    SYNOPSIS \n\
 \n\
          mris_calc [OPTIONS]					\\ \n\
          	<curvFile1> <ACTION> [<curvFile2> | <floatNumber>] \n\
 \n\
    DESCRIPTION \n\
 \n\
	'mris_calc' is a simple calculator that operates on FreeSurfer \n\
	curvature files. \n\
 \n\
	In most cases, the calculator functions with three arguments: \n\
        two inputs and an <ACTION> linking them. Some actions, however, \n\
        operate with only one input <curvFile1>. \n\
 \n\
        In all cases, the first input <curvFile1> is the name of a FreeSurfer \n\
        curvature file. \n\
 \n\
        For two inputs, the calculator first assumes that the second input \n\
        is the name of a second curvature file. If, however, this second input \n\
        is not found on the filesystem, the calculator attempts to parse this \n\
        input as a float number, which is then processed according to <ACTION>. \n\
 \n\
    OPTIONS \n\
 \n\
    	--output <outputCurvFile> \n\
   	 -o <outputCurvFile> \n\
 \n\
    	By default, 'mris_calc' will save the output curvature to a file \n\
    	in the current working directory called 'out.crv'. This name can be \n\
	set to something more meaningful with the '--output' option. \n\
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
	i.e. <curvFile1>[n] <ACTION> <curvFile2>[n] where 'n' is an index \n\
	counter. \n\
 \n\
	ACTION  INPUTS OUTPUTS	                EFFECT \n\
	  mul	   2      1     <outputCurvFile> = <curvFile1> * <curvFile2> \n\
	  div	   2	  1     <outputCurvFile> = <curvFile1> / <curvFile2> \n\
	  add      2	  1	<outputCurvFile> = <curvFile1> + <curvFile2> \n\
	  sub	   2	  1     <outputCurvFile> = <curvFile1> - <curvFile2> \n\
	  set	   2	  1     <outputCurvFile> = <curvFile2> \n\
 \n\
          ascii    1      1     <outputCurvFile> = ascii <curvFile1> \n\
 \n\
          size     1      0     print the size of <curvFile1> \n\
          min      1      0     print the min value (and index) of <curvFile1> \n\
          max      1      0     print the max value (and index) of <curvFile1> \n\
          mean     1      0     print the mean value of <curvFile1> \n\
          std      1      0     print the standard deviation of <curvFile1> \n\
          stats    1      0     process 'size', 'min', 'max', 'mean', 'std' \n\
 \n\
	The 'set' command is somewhat different in that for practical purposes \n\
	the contents of <curvFile1> are ignored. It is still important to \n\
	specifiy a valid <curvFile1> since it is parsed by 'mris_calc' \n\
	to determine the size of output curvature file to create. In most \n\
	instances, <curvFile2> will denote a float value, and not an actual \n\
	curvature file, i.e. 'mris_calc set rh.area 0.005' will create \n\
	an output curvature, 'out.crv' of the same size as rh.area, and with \n\
	each element set to 0.005. \n\
 \n\
        The 'ascii' command converts <curvFile1> to a text format file, \n\
        suitable for reading into MatLAB, for example. \n\
 \n\
        Note also that the standard deviation can suffer from float rounding \n\
        errors and is only accurate to 4 digits of precision. \n\
 \n\
    ARBITRARY FLOATS AS SECOND INPUT ARGUMENT \n\
 \n\
	If a second input argument is specified, 'mris_calc' will attempt to \n\
        open the argument following <ACTION> as if it were a curvature file. \n\
        Should this file not exist, 'mricurc_calc' will attempt to parse the \n\
        argument as if it were a float value. \n\
 \n\
	In such a case, 'mris_calc' will create a dummy internal \n\
	curvature file and set all its elements to this float value. \n\
 \n\
    NOTES \n\
 \n\
	<curvFile1> and <curvFile2> should typically be generated on the \n\
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
          "%s\nYou should specify '<curvFile1> <ACTION> [<curvFile2> | <floatNumber>]'",
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
error_unmatchedSizes(void) {
  char  pch_errorMessage[STRBUF];
  strcpy(pch_errorMessage, "");
  sprintf(pch_errorMessage,
          "%s\nlength(<curvFile1>) != length(<curvFile2>)",
          pch_errorMessage);
  sprintf(pch_errorMessage,
          "%s\nCurvature files should have the same size.",
          pch_errorMessage);
  ErrorExit(10, "%s: %s", G_pch_progname, pch_errorMessage);
}

void
output_init(void) {
  int i;
  G_sizeCurv3   = G_sizeCurv1;
  G_pf_arrayCurv3 = (float*) xmalloc(G_sizeCurv1 * sizeof(float));
  for(i=0; i<G_sizeCurv3; i++)
    G_pf_arrayCurv3[i] = 0.0;
}

void
curv2_fileRead() {
  float f_curv2 = 0.0;
  int   i;
  if(!CURV_fileRead(G_pch_curvFile2, &G_sizeCurv2, &G_pf_arrayCurv2)) {
    f_curv2   = atof(G_pch_curvFile2);
    G_sizeCurv2 = G_sizeCurv1;
    G_pf_arrayCurv2 = (float*) malloc(G_sizeCurv1 * sizeof(float));
    for(i=0; i<G_sizeCurv2; i++)
      G_pf_arrayCurv2[i] = f_curv2;
  }
  if(G_verbosity) cprintd("Size of curvFile1: ", G_sizeCurv1);
  if(G_verbosity) cprintd("Size of curvFile2: ", G_sizeCurv2);
  if(G_sizeCurv1 != G_sizeCurv2) error_unmatchedSizes();
}

int
main(
  int   argc,
  char  *argv[]
  ) {

  int   nargs;
  short ret       = 0;

  init();
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_calc.c,v 1.8 2008/09/24 20:34:54 rudolph Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  G_pch_progname = argv[0] ;

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

  //    Read in relevant curvature files -
  //+   the second curvature file is optional.
  ret     = CURV_fileRead(G_pch_curvFile1, &G_sizeCurv1, &G_pf_arrayCurv1);
  if(Gb_curvFile2) curv2_fileRead();
  output_init();
  ret     = CURV_process();

  return(0);
}

static int
options_print(void) {
  cprints("Curvature file 1:", G_pch_curvFile1);
  cprints("ACTION:", G_pch_operator);
  if(G_pch_curvFile2) cprints("Curvature file 2:", G_pch_curvFile2);
  cprints("Output curvature file:", G_pch_curvFile3);

  return 1;
}

static int
options_parse(int argc, char *argv[]) {
  int    nargs    = 0;
  char*  option;
  char*  pch_text;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-output") || (toupper(*option) == 'O')) {
    G_pch_curvFile3   = (argv[2]);
    nargs     = 1;
  } else if (!stricmp(option, "-version") || (toupper(*option) == 'V')) {
    version_print();
  } else if (!stricmp(option, "-verbosity")) {
    G_verbosity     = atoi(argv[2]);
    nargs     = 1;
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
  \fn CURV_fileRead(char* apch_curvFileName, int* ap_vectorSize, float* apf_curv)
  \brief Read a FreeSurfer curvature file into a float array
  \param apch_curvFileName The name of a FreeSurfer curvature file.
  \param ap_vectorSize Pointer to the size (i.e. number of elements) in file.
  \param apf_curv Array containing the curvature values.
  \return If curvature file is successfully read, return a 1, else return 0.
*/
#define   NEW_VERSION_MAGIC_NUMBER  16777215
short
CURV_fileRead(
  char* apch_curvFileName,
  int*  ap_vectorSize,
  float*  apf_curv[]
  ) {
  FILE* FP_curv;
  int   vnum;
  int   nvertices;
  int   i;
  char  pch_readMessage[STRBUF];
  float*  pf_data=NULL;

  if((FP_curv = fopen(apch_curvFileName, "r")) == NULL) return(0);
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
  } else
    ErrorExit(ERROR_BADPARM,
              "\n%s: curvature file '%s' has wrong magic number.\n",
              G_pch_progname, apch_curvFileName);
  *apf_curv   = pf_data;
  fclose(FP_curv);
  return(1);
}

/*!
  \fn ascii_fileWrite(char* apch_curvFileName, float* apf_curv)
  \brief Write a FreeSurfer curvature array to an ascii text file
  \param apch_curvFileName The name of a FreeSurfer curvature file.
  \param apf_curv Array containing the curvature values.
  \return If curvature file is successfully written, return a 1, else return 0.
*/
short
ascii_fileWrite(
  char* apch_curvFileName,
  float*  apf_curv
  ) {
  FILE* FP_curv;
  int   i;
  char  pch_readMessage[STRBUF];

  if((FP_curv = fopen(apch_curvFileName, "w")) == NULL)
    return(0);
  sprintf(pch_readMessage, "Writing %s", apch_curvFileName);
  for(i=0; i<G_sizeCurv1; i++) {
    CURV_arrayProgress_print(G_sizeCurv1, i, pch_readMessage);
    fprintf(FP_curv, "%f\n", apf_curv[i]);
  }
  fclose(FP_curv);
  return(1);
}


/*!
  \fn CURV_fileWrite(char* apch_curvFileName, int* ap_vectorSize, float* apf_curv)
  \brief Write a FreeSurfer curvature array to a file
  \param apch_curvFileName The name of a FreeSurfer curvature file.
  \param ap_vectorSize Pointer to the size (i.e. number of elements) in file.
  \param apf_curv Array containing the curvature values.
  \return If curvature file is successfully written, return a 1, else return 0.
*/
short
CURV_fileWrite(
  char* apch_curvFileName,
  int*  ap_vectorSize,
  float*  apf_curv
  ) {
  FILE* FP_curv;
  int   i;
  char  pch_readMessage[STRBUF];

  if((FP_curv = fopen(apch_curvFileName, "w")) == NULL)
    return(0);
  fwrite3(NEW_VERSION_MAGIC_NUMBER, FP_curv);
  fwriteInt(G_sizeCurv1, FP_curv);
  fwriteInt(G_nfaces, FP_curv);
  fwriteInt(G_valsPerVertex, FP_curv);
  sprintf(pch_readMessage, "Writing %s", apch_curvFileName);
  for(i=0; i<G_sizeCurv1; i++) {
    CURV_arrayProgress_print(G_sizeCurv1, i, pch_readMessage);
    fwriteFloat(apf_curv[i], FP_curv);
  }
  fclose(FP_curv);
  return(1);
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
  if(acurrent%fivePerc == fivePerc-1) fprintf(G_FP, "#");
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
    	(!strcmp(apch_operator, "ascii")) )
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

  if(!strcmp(G_pch_operator, "ascii"))
    ascii_fileWrite(G_pch_curvFile3, G_pf_arrayCurv1);

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

  if(b_canWrite) CURV_fileWrite(G_pch_curvFile3, &G_sizeCurv3, G_pf_arrayCurv3);

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
