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
 *    $Author: nicks $
 *    $Date: 2009/04/03 17:33:00 $
 *    $Revision: 1.2.2.2 $
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

#define  STRBUF   	65536
#define  MAX_FILES    	1000
#define  CO( x )    	fprintf(stdout, ( x ))
#define  CE( x )    	fprintf(stderr, ( x ))
#define  START_i    	3

static const char vcid[] =
"$Id: mris_calc.c,v 1.2.2.2 2009/04/03 17:33:00 nicks Exp $";

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
  e_sqd,
  e_sqr,
  e_sqrt,
  e_set,
  e_abs,
  e_sign,
  e_lt,
  e_lte,
  e_gt,
  e_gte,
  e_masked,
  e_and,
  e_or,
  e_andbw,
  e_orbw,
  e_upperlimit,
  e_lowerlimit,
  e_ascii,
  e_stats,
  e_min,
  e_max,
  e_size,
  e_mini,
  e_maxi,
  e_mean,
  e_std,
  e_sum,
  e_prod,
  e_norm
} e_operation;

const char* Gppch_operation[] = {
  "multiply",
  "divide",
  "add",
  "subtract",
  "square difference",
  "square",
  "square root"
  "set",
  "abs",
  "sign",
  "less than",
  "less than or equal to",
  "greated than",
  "greater than or equal to",
  "masked by",
  "logical and",
  "logical or",
  "bitwise and",
  "bitwise or",
  "upper limit",
  "lower limit",
  "ascii",
  "stats",
  "min",
  "max",
  "size",
  "mini",
  "maxi",
  "mean",
  "std",
  "sum",
  "product",
  "normalize"
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
static e_operation      Ge_operation            = e_stats;

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
static e_operation operation_lookup(char* apch_operation);

// Simple functions on two float arguments
double fn_mul(float af_A, float af_B)   {return (af_A * af_B);}
double fn_div(float af_A, float af_B)   {return (af_B != 0 ?
                                        (af_A / af_B) : 0.0);}
double fn_add(float af_A, float af_B)   {return (af_A + af_B);}
double fn_sub(float af_A, float af_B)   {return (af_A - af_B);}
double fn_sqd(float af_A, float af_B)   {return (af_A - af_B)*(af_A - af_B);}
double fn_set(float af_A, float af_B)   {return (af_B);}

// Simple relational functions on two float arguments
double  fn_lt(float af_A,  float af_B)  {return (af_A < af_B ? af_A : 0.0);}
double  fn_lte(float af_A, float af_B)  {return (af_A <= af_B? af_A : 0.0);}
double  fn_gt(float af_A,  float af_B)  {return (af_A > af_B ? af_A : 0.0);}
double  fn_gte(float af_A, float af_B)  {return (af_A >= af_B? af_A : 0.0);}
double  fn_upl(float af_A, float af_B)  {return (af_A >= af_B ? af_B : af_A);}
double  fn_lrl(float af_A, float af_B)  {return (af_A <= af_B ? af_B : af_A);}
  
// Simple logical functions on two float arguments
double  fn_and(float af_A,    float af_B)       {return (af_A && af_B);}
double  fn_or(float af_A,     float af_B)       {return (af_A || af_B);}
double  fn_andbw(float af_A,  float af_B)       {return (
                                                  (int) af_A & (int) af_B);}
double  fn_orbw(float af_A,   float af_B)       {return (
                                                  (int) af_A | (int) af_B);}
double  fn_masked(float af_A, float af_B)       {return(af_B ? af_A : af_B);}

                

// Simple functions on one argument
double fn_abs(float af_A)		{return fabs(af_A);}
double fn_sqr(float af_A)		{return (af_A*af_A);}
double fn_sqrt(float af_A)		{return sqrt(af_A);}

double fn_sign(float af_A) {
    float f_ret=0.0;
    if(af_A < 0.0)	f_ret 	= -1.0;
    if(af_A == 0.0)	f_ret	= 0.0;
    if(af_A > 0.0)	f_ret	= 1.0;
    return f_ret;
}

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

// Additional 'sum' and 'prod' functions, but using a global vars
//+ as opposed to local static. The caller must make sure to zero the global
//+ var where appropriate.
double Gf_sum		= 0.;
double Gf_prod		= 1.;
double Gfn_sum(float af_A) {
    Gf_sum += af_A;
    return Gf_sum;
}
double Gfn_prod(float af_A) {
    Gf_prod *= af_A;
    return Gf_prod;
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

short   b_outCurvFile_write(e_operation e_op);
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
 \n\
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
	The action to be perfomed on the two input files. This is a \n\
	text string that defines the mathematical operation to execute. For \n\
	two inputs, this action is applied in an indexed element-by-element \n\
	fashion, i.e. \n\
	 \n\
			<file3>[n] = <file1>[n] <ACTION> <file2>[n] \n\
 \n\
	where 'n' is an index counter into the data space. \n\
 \n\
	ACTION  INPUTS OUTPUTS	                EFFECT \n\
     MATHEMATICAL \n\
	  mul	   2      1     <outputFile> = <file1> * <file2> \n\
	  div	   2	  1     <outputFile> = <file1> / <file2> \n\
	  add      2	  1	<outputFile> = <file1> + <file2> \n\
	  sub	   2	  1     <outputFile> = <file1> - <file2> \n\
	  sqd	   2	  1     <outputFile> = (<file1> - <file2>)^2 \n\
	  set	   2	  1     <file1>	     = <file2> \n\
	  sqr	   1      1	<outputFile> = <file1> * <file1> \n\
	  sqrt	   1      1     <outputFile> = sqrt(<file1>) \n\
          abs      1      1     <outputFile> = abs(<file1>) \n\
	  sign	   1	  1	<outputFile> = sign(<file1>) \n\
          norm     1      1     <outputFile> = norm(<file1>) \n\
 \n\
      RELATIONAL \n\
          lt       2      1     <outputFile> = <file1> <  <file2> \n\
          lte      2      1     <outputFile> = <file1> <= <file2> \n\
          gt       2      1     <outputFile> = <file1> >  <file2> \n\
          gte      2      1     <outputFile> = <file1> >= <file2> \n\
          upl      2      1     <outputFile> = <file1> upperlimit <file2> \n\
          lrl      2      1     <outputFile> = <file1> lowerlimit <file2> \n\
 \n\
        LOGICAL \n\
          and      2      1     <outputFile> = <file1> && <file2> \n\
           or      2      1     <outputFile> = <file1> || <file2> \n\
         andbw     2      1     <outputFile> = (int)<file1> & (int)<file2> \n\
          orbw     2      1     <outputFile> = (int)<file1> | (int)<file2> \n\
         masked    2      1     <outputFile> = <file1> maskedby <file2> \n\
 \n\
     DATA CONVERSION \n\
          ascii    1      1     <outputFile> = ascii <file1> \n\
 \n\
      STATISTICAL \n\
          size     1      0     print the size (number of elements) of <file1> \n\
          min      1      0     print the min value (and index) of <file1> \n\
          max      1      0     print the max value (and index) of <file1> \n\
          mean     1      0     print the mean value of <file1> \n\
          std      1      0     print the standard deviation of <file1> \n\
	  sum	   1      0	print the sum across all values of <file1> \n\
	  prod	   1      0     print the inner product across <file1> \n\
          stats    1      0     process 'size', 'min', 'max', 'mean', 'std' \n\
 \n\
    NOTES ON ACTIONS \n\
 \n\
      MATHEMATICAL \n\
        The 'add', 'sub', 'div', 'mul' operations all function as one would \n\
        expect. The 'norm' creates an output file such that all values are \n\
        constrained (normalized) between 0.0 and 1.0. The 'sqd' stores the \n\
        square difference between two inputs. \n\
 \n\
        The 'sqr' and 'sqrt' return the square and square-root of an input \n\
        file. \n\
	 \n\
	The 'sign' function returns at each index of an input file: \n\
 \n\
		-1 if <file1>[n] <  0 \n\
		 0 if <file1>[n] == 0 \n\
		 1 if <file1>[n] >  0 \n\
 \n\
        NOTE: For volume files, be very careful about data types! If the input \n\
        volume is has data of type UCHAR (i.e. integers between zero and 255), \n\
        the output will be constrained to this range and type as well! That \n\
        means, simply, that if type UCHAR vols are multiplied together, the \n\
        resultant output will itself be a UCHAR volume. This is probably not \n\
        what you want. In order for calculations to evaluate correctly, \n\
        especially 'mul', 'div', and 'norm', convert the input volumes to \n\
        float format, i.e.: \n\
 \n\
                $>mri_convert -odt float input.mgz input.f.mgz \n\
                $>mris_calc -o input_norm.mgz input.f.mgz norm \n\
 \n\
        will give correct results, while \n\
 \n\
                $>mris_calc -o input_norm.mgz input.mgz norm \n\
 \n\
        most likely be *not* what you are looking for. \n\
 \n\
	The 'set' command overwrites its input data. It still requires a valid \n\
	<file1> -- since in most instances the 'set' command is used to set \n\
	input data values to a single float constant, i.e. \n\
 \n\
		$>mris_calc rh.area set 0.005 \n\
 \n\
	will set all values of rh.area to 0.005. It might be more meaningful \n\
	to first make a copy of the input file, and set this \n\
 \n\
		$>cp rh.area rh-0.005 \n\
		$>mris_calc rh-0.005 set 0.005 \n\
 \n\
	Similarly for volumes \n\
 \n\
		$>cp orig.mgz black.mgz \n\
		$>mris_calc black.mgz set 0 \n\
 \n\
	will result in the 'black.mgz' volume having all its intensity values \n\
	set to 0. \n\
 \n\
    RELATIONAL \n\
        The relational operations apply the relevant evaluation at each \n\
        index in the input data sets. \n\
 \n\
              'lt'      -- less than \n\
              'gt'      -- greater than \n\
              'lte'     -- less than or equal to \n\
              'gte'     -- greater than or equal to \n\
 \n\
        If the comparison is valid for a point in <file1> compared to \n\
        corresponding point in <file2>, the <file1> value is retained; \n\
        otherwise the <file1> value is set to zero. Thus, if we run \n\
        'mris_calc input1.mgz lte input2.mgz', the output volume 'out.mgz' \n\
        will have all input1.mgz values that are not less than or equal to \n\
        to input2.mgz set to zero. \n\
 \n\
        The 'upl' and 'lrl' are hardlimiting filters. In the case of 'upl', \n\
        any values in <file1> that are greater or equal to corresponding data \n\
        points in <file2> are limited to the values in <file2>. Similarly, for \n\
        'lpl', any <file1> values that are less than corresponding <file2> \n\
        values are set to these <file2> values. \n\
 \n\
        Essentially, 'upl' guarantees that all values in <file1> are less than \n\
        or at most equal to corresponding values in <file2>; 'lpl' guarantees \n\
        that all values in <file1> are greater than or at least equal to \n\
        corresponding values in <file2>. \n\
 \n\
    LOGICAL \n\
        The logical operations follow C convention, i.e. and is a logical 'and' \n\
        equivalent to the C '&&' operator, similarly for 'or' which is \n\
        evaluated with the C '||'. The 'andbw' and 'orbw' are bit-wise \n\
        operators, evaluted with the C operators '&' and '|' respectively. \n\
 \n\
    DATA CONVERSION \n\
        The 'ascii' command converts <file1> to a text format file, \n\
        suitable for reading into MatLAB, for example. Note that for volumes \n\
	data values are written out as a 1D linear array with looping order \n\
	(slice, height, width). \n\
 \n\
    STATISTICAL \n\
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
	Multiply each value in <rh.area> with the corresponding value \n\
	in <rh.thickness>, creating a new file called 'out.crv' that \n\
	contains the result. \n\
 \n\
    	$>mris_calc --output rh.weightedCortex rh.area mul rh.thickness \n\
	Same as above, but give the ouput file the more meaningful name \n\
	of 'rh.weightedCortex'. \n\
 \n\
        $>mris_calc rh.area max \n\
        Determine the maximum value in 'rh.area' and print to stdout. In \n\
        addition to the max value, the index offset in 'rh.area' containing \n\
        this value is also printed. \n\
 \n\
        $>mris_calc rh.area stats \n\
        Determine the size, min, max, mean, and std of 'rh.area'. \n\
 \n\
	$>mris_calc orig.mgz sub brainmask.mgz \n\
	Subtract the brainmask.mgz volume from the orig.mgz volume. Result is \n\
	saved by default to out.mgz. \n\
 \n\
        $>mris_calc -o ADC_masked.nii ADC.nii masked B0_mask.img \n\
        Mask a volume 'ADC.nii' with 'B0_mask.img', saving the output in \n\
        'ADC_masked.nii'. Note that the input volumes are different formats, \n\
        but the same logical size. \n\
 \n\
    ADVANCED EXAMPLES \n\
 \n\
	Consider the case when calculating the right hemisphere pseudo volume \n\
	formed by the FreeSurfer generated white matter 'rh.area' curvature \n\
	file, and the cortical thickness, 'rh.thickness'. Imagine this is to \n\
	be expressed as a percentage of intercranial volume. \n\
 \n\
	First, calculate the volume and store in a curvature format: \n\
 \n\
		$>mris_calc -o rh.cortexVol rh.area mul rh.thickness \n\
 \n\
	Now, find the intercranial volume (ICV) in the corresponding output \n\
	file generated by FreeSurfer for this subject. Assume ICV = 100000. \n\
 \n\
		$>mris_calc -o rh.cortexVolICV rh.cortexVol div 100000 \n\
	 \n\
	Here the second <ACTION> argument is a number and not a curvature file. \n\
 \n\
	We could have achieved the same effect by first creating an \n\
	intermediate curvature file, 'rh.ICV' with each element set to \n\
	the ICV, and then divided by this curvature: \n\
	 \n\
		$>cp rh.area rh.ICV \n\
		$>mris_calc rh.ICV set 100000 \n\
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
  fprintf(stderr, "\n");
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
error_unknownFileType(char* apch_inputFile) {
  char  pch_action[65536];
  sprintf(pch_action, "accessing file '%s'", apch_inputFile);
  error_exit(   pch_action,
                "I was unable to identify this file type. Does it exist?",
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

  if(type == -1)
    error_unknownFileType(apch_inputFile);

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

short
CURV_set(
  float*        apf_array[], 
  int           a_size, 
  float         af_val) {
  //
  // PRECONDITIONS
  // o apf_array does not need to have been allocated... if non-NULL will
  //   first free.
  //   
  // POSTCONDITIONS
  // o apf_array is set to have all values set to <af_val>.
  // 
  int           i;
  float*        pf_array        = NULL;

//   if(*apf_array) free(*apf_array);
  pf_array = (float*) malloc(a_size * sizeof(float));
  for(i=0; i<a_size; i++)
    pf_array[i] = af_val;
  *apf_array    = pf_array;
  return 1;
}

short
CURV_copy(
  float         apf_arraySrc[],
  float         apf_arrayDst[],
  int           a_size
) {
  //
  // PRECONDITIONS
  // o <a_size> must be correct for both array inputs.
  // 
  // POSTCONDITIONS
  // o contents of <apf_arrayDst> are copied to <apf_arraySrc>, i.e.
  //   <apf_arraySrc> = <apf_arrayDst>
  // 
  int   i;
  for(i=0; i<a_size; i++)
    apf_arraySrc[i]     = apf_arrayDst[i];
  return 1;
}


int
main(
  int   argc,
  char  *argv[]
  ) {

  int   	nargs;
  short 	ret		= 0;
  float		f_curv2;
  e_FILEACCESS	eACCESS;

  init();
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_calc.c,v 1.2.2.2 2009/04/03 17:33:00 nicks Exp $",
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
        CURV_set(&G_pf_arrayCurv2, G_sizeCurv2, f_curv2);
/*    	G_pf_arrayCurv2 = (float*) malloc(G_sizeCurv1 * sizeof(float));
    	for(i=0; i<G_sizeCurv2; i++)
      	    G_pf_arrayCurv2[i] = f_curv2;*/
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
  \fn static e_operation operation_lookup(char* apch_operation)
  \brief Convert a string 'operation' to an enumerated e_operation
  \param apch_operation The mathematical operation to perform..
  \return The internal lookup enumerated e_operation type is returned.
*/static e_operation
operation_lookup(
  char* apch_operation
) {

  e_operation   e_op = e_stats;

  if(     !strcmp(apch_operation, "mul"))       e_op    = e_mul;
  else if(!strcmp(apch_operation, "div"))       e_op    = e_div;
  else if(!strcmp(apch_operation, "add"))       e_op    = e_add;
  else if(!strcmp(apch_operation, "sub"))       e_op    = e_sub;
  else if(!strcmp(apch_operation, "sqd"))       e_op    = e_sqd;
  else if(!strcmp(apch_operation, "sqr"))       e_op    = e_sqr;
  else if(!strcmp(apch_operation, "sqrt"))	e_op    = e_sqrt;
  else if(!strcmp(apch_operation, "set"))       e_op    = e_set;
  else if(!strcmp(apch_operation, "abs"))       e_op    = e_abs;
  else if(!strcmp(apch_operation, "sign"))      e_op    = e_sign;

  else if(!strcmp(apch_operation, "lt"))        e_op    = e_lt;
  else if(!strcmp(apch_operation, "lte"))       e_op    = e_lte;
  else if(!strcmp(apch_operation, "gt"))        e_op    = e_gt;
  else if(!strcmp(apch_operation, "gte"))       e_op    = e_gte;
  else if(!strcmp(apch_operation, "masked"))    e_op    = e_masked;
  else if(!strcmp(apch_operation, "and"))       e_op    = e_and;
  else if(!strcmp(apch_operation, "or"))        e_op    = e_or;
  else if(!strcmp(apch_operation, "andbw"))     e_op    = e_andbw;
  else if(!strcmp(apch_operation, "orbw"))      e_op    = e_orbw;
  else if(!strcmp(apch_operation, "upperlimit"))e_op    = e_upperlimit;
  else if(!strcmp(apch_operation, "lowerlimit"))e_op    = e_lowerlimit;

  else if(!strcmp(apch_operation, "stats"))     e_op    = e_stats;
  else if(!strcmp(apch_operation, "size"))      e_op    = e_size;
  else if(!strcmp(apch_operation, "min"))       e_op    = e_min;
  else if(!strcmp(apch_operation, "max"))       e_op    = e_max;
  else if(!strcmp(apch_operation, "mean"))      e_op    = e_mean;
  else if(!strcmp(apch_operation, "std"))       e_op    = e_std;
  else if(!strcmp(apch_operation, "sum"))       e_op    = e_sum;
  else if(!strcmp(apch_operation, "prod"))      e_op    = e_prod;
  else if(!strcmp(apch_operation, "norm"))      e_op    = e_norm;

  else if(!strcmp(apch_operation, "ascii"))     e_op    = e_ascii;

  return e_op;
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
	pf_data[I++]	= (float) MRIgetVoxVal(pMRI, i, j, k, 0);
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
  // o Gp_MRI saved to <apchq_volFileName>.
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
	MRIsetVoxVal(Gp_MRI, i, j, k, 0, (float) apf_data[I++]);
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
  \fn b_outCurvFile_write(e_operation e_op)
  \brief A simple function that determines if an output file should be saved.
  \param e_op The enumerated operation to perform.
  \see
  \return A simple TRUE or FALSE indicating if for the passed operation an output file is generated.
*/
short
b_outCurvFile_write(e_operation e_op)
{
    // PRECONDITIONS
    //	e_op		in		Action to perform
    //
    // POSTCONDITIONS
    //	Internally an output curvature data structure is always generated.
    //	Saving this output file to disk is only meaningful for a given
    //	set of actions -- usually the typical mathematical operators.
    //	Other actions don't require an output file to be saved, e.g the
    //	'stats' action.
    //

    short	b_ret	= 0;

    if(
        e_op == e_mul           ||
        e_op == e_div           ||
        e_op == e_add           ||
        e_op == e_sub           ||
        e_op == e_sqd           ||
        e_op == e_sqr           ||
        e_op == e_sqrt          ||
        e_op == e_set           ||
        e_op == e_abs           ||
	e_op == e_sign		||
        e_op == e_lt            ||
        e_op == e_lte           ||
        e_op == e_gt            ||
        e_op == e_gte           ||
        e_op == e_masked        ||
        e_op == e_and           ||
        e_op == e_or            ||
        e_op == e_andbw         ||
        e_op == e_orbw          ||
        e_op == e_upperlimit    ||
        e_op == e_lowerlimit    ||
        e_op == e_norm
  ) b_ret = 1;

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
  float f_range         = 0.;
  int   mini            = -1;
  int   maxi            = -1;
  int   b_canWrite      = 0;
  float f_mean   	= 0.;
  float f_std   	= 0.;
  float f_dev   	= 0.;
  char	pch_text[STRBUF];

  Ge_operation  = operation_lookup(G_pch_operator);
  b_canWrite    = b_outCurvFile_write(Ge_operation);
  switch(Ge_operation) {
    case  e_mul:        CURV_functionRunABC(fn_mul);    break;
    case  e_div:        CURV_functionRunABC(fn_div);    break;
    case  e_add:        CURV_functionRunABC(fn_add);    break;
    case  e_sub:        CURV_functionRunABC(fn_sub);    break;
    case  e_set:        CURV_functionRunABC(fn_set);    break;
    case  e_sqd:        CURV_functionRunABC(fn_sqd);    break;
    case  e_abs:        CURV_functionRunAC( fn_abs);    break;
    case  e_sign:	CURV_functionRunAC( fn_sign);	break;
    case  e_sqr:        CURV_functionRunAC( fn_sqr);    break;
    case  e_sqrt:       CURV_functionRunAC( fn_sqrt);   break;
    case  e_lt:         CURV_functionRunABC(fn_lt);     break;
    case  e_lte:        CURV_functionRunABC(fn_lte);    break;
    case  e_gt:         CURV_functionRunABC(fn_gt);     break;
    case  e_gte:        CURV_functionRunABC(fn_gte);    break;
    case  e_masked:     CURV_functionRunABC(fn_masked); break;
    case  e_and:        CURV_functionRunABC(fn_and);    break;
    case  e_or:         CURV_functionRunABC(fn_or);     break;
    case  e_andbw:      CURV_functionRunABC(fn_andbw);  break;
    case  e_orbw:       CURV_functionRunABC(fn_orbw);   break;
    case  e_upperlimit: CURV_functionRunABC(fn_upl);    break;
    case  e_lowerlimit: CURV_functionRunABC(fn_lrl);    break;
    case  e_ascii:                                      break;
    case  e_stats:                                      break;
    case  e_size:                                       break;
    case  e_min:                                        break;
    case  e_mini:                                       break;
    case  e_maxi:                                       break;
    case  e_max:                                        break;
    case  e_mean:                                       break;
    case  e_std:                                        break;
    case  e_sum:                                        break;
    case  e_prod:                                       break;
    case  e_norm:
                f_min   = CURV_functionRunAC(fn_min);
                f_max   = CURV_functionRunAC(fn_max);
                f_range = f_max - f_min;
                // v2   = f_min
                CURV_set(&G_pf_arrayCurv2, G_sizeCurv1, f_min);
                // v3 = v1 - v2
                CURV_functionRunABC(fn_sub);
                // v1 = v3
                CURV_copy(G_pf_arrayCurv1, G_pf_arrayCurv3, G_sizeCurv1);
                // v2 = f_range
                CURV_set(&G_pf_arrayCurv2, G_sizeCurv1, f_range);
                // v3 = v1 / v2
                CURV_functionRunABC(fn_div);
                break;
  }

  if(Ge_operation == e_set) strcpy(G_pch_curvFile3, G_pch_curvFile1);

  if(Ge_operation == e_ascii) {
    if(!Gb_file3) strcat(G_pch_curvFile3, ".ascii");
	printf("Write : %s", G_pch_curvFile3);

    if((ascii_fileWrite(G_pch_curvFile3, G_pf_arrayCurv1))==e_WRITEACCESSERROR)
	printf("Write error: %s", G_pch_curvFile3);
  }

  if(Ge_operation == e_size || Ge_operation == e_stats) {
    cprintd("Size", G_sizeCurv1);
  }

  if(Ge_operation == e_min || Ge_operation == e_stats) {
    f_min       = CURV_functionRunAC(fn_min);
    mini        = (int) CURV_functionRunAC(fn_mini);
    sprintf(pch_text, "%f (%d)", f_min, mini);    
    cprints("Min@(index)", pch_text);
  }
  if(Ge_operation == e_max || Ge_operation == e_stats) {
    f_max       = CURV_functionRunAC(fn_max);
    maxi        = (int) CURV_functionRunAC(fn_maxi);
    sprintf(pch_text, "%f (%d)", f_max, maxi);    
    cprints("Max@(index)", pch_text);
  }
  if(Ge_operation == e_mean || Ge_operation == e_stats) {
    f_mean      = CURV_functionRunAC(fn_mean);
    cprintf("Mean", f_mean);
  }
  if(Ge_operation == e_std || Ge_operation == e_stats) {
    f_dev       = CURV_functionRunAC(fn_dev);
    f_std       = sqrt(f_dev/(G_sizeCurv1-1));
    cprintf("Std", f_std);
  }
  if(Ge_operation == e_sum || Ge_operation == e_stats) {
    Gf_sum	= 0.0;
    Gf_sum	= CURV_functionRunAC(Gfn_sum);
    cprintf("Sum", Gf_sum);
  }
  if(Ge_operation == e_prod || Ge_operation == e_stats) {
    Gf_prod	= 1.0;
    Gf_prod	= CURV_functionRunAC(Gfn_prod);
    cprintf("Prod", Gf_prod);
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
