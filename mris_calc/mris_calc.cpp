/**
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
#define _ISOC99_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/stat.h>
#include <libgen.h>

#include <assert.h>
#include <errno.h>

#include <getopt.h>
#include <stdarg.h>
#include <float.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "fio.h"
#include "mri_identify.h"
#include "label.h"

#include "romp_support.h"


#define  STRBUF         65536
#define  MAX_FILES      1000
#define  CO( x )        fprintf(stdout, ( x ))
#define  CE( x )        fprintf(stderr, ( x ))
#define  START_i        3

double fn_sign(float af_A);

// ----------------------------------------------------------------------------
// DECLARATION
// ----------------------------------------------------------------------------

// ------------------
// typedefs and enums
// ------------------

typedef enum _FILETYPE
{
  e_Unknown, e_VolumeFile, e_CurvatureFile, e_SurfaceFile, e_FloatArg
} e_FILETYPE;

const char*   Gppch_filetype[] =
{
  "Unknown",
  "Volume",
  "Curvature",
  "Surface",
  "FloatArg"
};

const char* Gppch_fileExt[] =
{
  "null",
  "mgz",
  "crv",
  "surf",
  "dummy"
};

const char* Gppch_fileDotExt[] =
{
  ".null",
  ".mgz",
  ".crv",
  ".surf",
  ".dummy"
};

typedef enum _FILEACCESS
{
  e_UNSPECIFIED     = -10,
  e_WRONGMAGICNUMBER    = -1,
  e_OK          =  0,
  e_READACCESSERROR   =  1,
  e_WRITEACCESSERROR    =  2
} e_FILEACCESS;

typedef enum _operation
{
  e_mul,
  e_div,
  e_mod,
  e_add,
  e_pow,
  e_sub,
  e_sub0,
  e_sratio,
  e_pctdiff,
  e_pctdiff0,
  e_sqd,
  e_sqr,
  e_not,
  e_sqrt,
  e_set,
  e_atan2,
  e_bcor,
  e_sig2p,
  e_mag,
  e_log10,
  e_abs,
  e_inv,
  e_sign,
  e_eq,
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
  e_norm,
  e_log,
  e_unknown
} e_operation;

const char* Gppch_operation[] =
{
  "multiply",
  "divide",
  "modulus",
  "add",
  "pow",
  "subtract",
  "subtract0",
  "sratio",
  "pctdiff",
  "pctdiff0",
  "square difference",
  "square",
  "square root"
  "set",
  "atan2",
  "bcor",
  "sig2p",
  "mag",
  "log10",
  "abs",
  "inv",
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
  "normalize",
  "log",
  "unknown operation"
};

// -------------------------------
// Global "class" member variables
//--------------------------------

char*             G_pch_progname ;
const char*             Progname ;

static int        G_verbosity             = 0;
static FILE*      G_FP                    = NULL;
static int    Gb_strictExtensions   = 0;

// Input 1
static unsigned long  G_sizeCurv1             = 0;
static char*      G_pch_curvFile1         = NULL;
static float*     G_pf_arrayCurv1         = NULL;
static float*     G_pf_arrayCurv1Copy     = NULL;
static int        G_sizeCurv1Copy         = 0;
static int        G_nfaces                = 0;
static int        G_valsPerVertex         = 0;
static int        G_FSFILETYPE1           = -1; // FreeSurfer type designator
static e_FILETYPE G_eFILETYPE1            = e_Unknown;

// Input 2
static int        Gb_curvFile2            = 0;  //  The second input
static char*      G_pch_curvFile2         = NULL; //+ file is optional.
static unsigned long        G_sizeCurv2             = 0;
static float*     G_pf_arrayCurv2         = NULL;
static float*     G_pf_arrayCurv2Copy     = NULL;
static unsigned long       G_sizeCurv2Copy         = 0;
static int        G_FSFILETYPE2           = -1; // FreeSurfer type designator
static e_FILETYPE G_eFILETYPE2            = e_Unknown;

// "Helper" pointers
static MRI*       Gp_MRI            = NULL; // Pointer to most
//+ recently read
//+ MRI_VOLUME struct
//+ and used in volume
//+ handling... mostly
//+ for volume size.

// Operation to perform on input1 and input2
static char*            G_pch_operator          = NULL;
static e_operation      Ge_operation            = e_unknown;

// Output file
static unsigned long    G_sizeCurv3             = 0;
static char             G_pch_curvFile3[STRBUF];
static float*           G_pf_arrayCurv3         = NULL;
static float*           G_pf_arrayCurv3Copy     = NULL;
static unsigned long    G_sizeCurv3Copy         = 0;
static short            Gb_file3                = 0;
static short            Gb_canWrite             = 0;
static int              G_FSFILETYPE3           = -1; // FreeSurfer type designator
static e_FILETYPE       G_eFILETYPE3            = e_Unknown;

// Label file and data
static short            Gb_labelMask            = 0;
static char             G_pch_labelFile[STRBUF];
static LABEL*           G_plabel                = NULL;
static float*           G_pf_labelMask          = NULL;
static float*           G_pf_labelMaskIndex     = NULL;
static int              G_sizeLabel             = 0;

//----------------
// "class" methods
//----------------

// Housekeeping functions
static void synopsis_show(void);
static void version_print(void);

// Setup functions
static int  options_parse(
  int   argc,
  char* apchv[]
);
static int  options_print(void);
static e_operation operation_lookup(char* apch_operation);

// Simple functions on two float arguments
double fn_mul(float af_A, float af_B)
{
  return (af_A * af_B);
}
double fn_div(float af_A, float af_B)
{
  return (af_B != 0 ?
          (af_A / af_B) : 0.0);
}
double fn_mod(float af_A, float af_B)
{
  return fmod(af_A, af_B);
}
double fn_add(float af_A, float af_B)
{
  return (af_A + af_B);
}
double fn_pow(float af_A, float af_B)
{
  return ( pow(af_A,af_B) );
}
double fn_sub(float af_A, float af_B)
{
  return (af_A - af_B);
}
double fn_sub0(float af_A, float af_B)
{
  // Return 0 if either A or B are 0
  if(fabs(af_A) < 2.0*FLT_MIN) return(0);
  if(fabs(af_B) < 2.0*FLT_MIN) return(0);
  return (af_A - af_B);
}
double fn_sratio(float af_A, float af_B)
{
  if(fabs(af_A) < 2.0*FLT_MIN) return(0);
  if(fabs(af_B) < 2.0*FLT_MIN) return(0);
  if(af_A > af_B) return(af_A/(af_B + FLT_MIN));
  return(-af_B/(af_A + FLT_MIN));
}
double fn_pctdiff(float af_A, float af_B)
{
  return (100*(af_A - af_B)/( (af_A + af_B + FLT_MIN)/2) );
}
double fn_pctdiff0(float af_A, float af_B)
{
  // Return 0 if either A or B are 0
  if(fabs(af_A) < 2.0*FLT_MIN) return(0);
  if(fabs(af_B) < 2.0*FLT_MIN) return(0);
  return (100*(af_A - af_B)/( (af_A + af_B + FLT_MIN)/2) );
}
double fn_sqd(float af_A, float af_B)
{
  return (af_A - af_B)*(af_A - af_B);
}
double fn_set(float af_A, float af_B)
{
  return (af_B);
}
double fn_atan2(float af_A, float af_B)
{
  return (atan2(af_A,af_B));
}
double fn_bcor(float af_A, float af_B)
{
  // bonferroni correction, assumes input is sig (ie, -log10(p))
  double p;
  p = pow(10,-fabs(af_A));
  p = 1 - pow(1-p,af_B);
  return (-log10(p)*fn_sign(af_A));
}
double fn_sig2p(float af_A)
{
  return (pow(10,-fabs(af_A)));
}
double fn_mag(float af_A, float af_B)
{
  return (sqrt(af_A*af_A + af_B*af_B));
}

// Simple relational functions on two float arguments
double  fn_eq(float af_A,  float af_B)
{
  return (af_A == af_B ? af_A : 0.0);
}
double  fn_lt(float af_A,  float af_B)
{
  return (af_A < af_B ? af_A : 0.0);
}
double  fn_lte(float af_A, float af_B)
{
  return (af_A <= af_B? af_A : 0.0);
}
double  fn_gt(float af_A,  float af_B)
{
  return (af_A > af_B ? af_A : 0.0);
}
double  fn_gte(float af_A, float af_B)
{
  return (af_A >= af_B? af_A : 0.0);
}
double  fn_upl(float af_A, float af_B)
{
  return (af_A >= af_B ? af_B : af_A);
}
double  fn_lrl(float af_A, float af_B)
{
  return (af_A <= af_B ? af_B : af_A);
}

// Simple logical functions on two float arguments
double  fn_and(float af_A,    float af_B)
{
  return (af_A && af_B);
}
double  fn_or(float af_A,     float af_B)
{
  return (af_A || af_B);
}
double  fn_andbw(float af_A,  float af_B)
{
  return (
           (int) af_A & (int) af_B);
}
double  fn_orbw(float af_A,   float af_B)
{
  return (
           (int) af_A | (int) af_B);
}
double  fn_masked(float af_A, float af_B)
{
  return(af_B ? af_A : af_B);
}

// Simple functions on one argument
double fn_log10(float af_A)
{
  if(af_A <= 0) return(0);
  return(log10(af_A));
}
double fn_log(float af_A)
{
  if(af_A <= 0) return(0);
  return(log(af_A));
}
double fn_abs(float af_A)
{
  return fabs(af_A);
}
double fn_inv(float af_A)
{
  return(1.0/af_A);
}
double fn_sqr(float af_A)
{
  return (af_A*af_A);
}
double fn_not(float af_A)
{
  if(af_A > 0.5)
  {
    return(0);
  }
  return(1);
}
double fn_sqrt(float af_A)
{
  return sqrt(af_A);
}
double fn_sign(float af_A)
{
  float f_ret=0.0;
  if(af_A < 0.0)
  {
    f_ret   = -1.0;
  }
  if(af_A == 0.0)
  {
    f_ret = 0.0;
  }
  if(af_A > 0.0)
  {
    f_ret = 1.0;
  }
  return f_ret;
}

double fn_zeroCount(float af_A)
{
  static int          count = 0;
  if(af_A == 0.0)
  {
    count++;
  }
  return count;
}
double fn_negCount(float af_A)
{
  static int          count = 0;
  if(af_A < 0.0)
  {
    count++;
  }
  return count;
}
double fn_posCount(float af_A)
{
  static int          count = 0;
  if(af_A > 0.0)
  {
    count++;
  }
  return count;
}

double fn_min(float af_A)
{
  static float        f_min  = 0;
  static int          count  = 0;
  if(!count)
  {
    f_min = af_A;
  }
  if(af_A <= f_min)
  {
    f_min = af_A;
  }
  count++;
  return f_min;
}
double fn_mini(float af_A)
{
  static float        f_min  = 0;
  static int          mini   = -1;
  static int          count  = 0;
  if(!count)
  {
    f_min = af_A;
  }
  if(af_A <= f_min)
  {
    f_min   = af_A;
    mini    = count;
  }
  count++;
  return (float) mini;
}
double fn_max(float af_A)
{
  static float        f_max  = 0;
  static int          count  = 0;
  if(!count)
  {
    f_max = af_A;
  }
  if(af_A >= f_max)
  {
    f_max = af_A;
  }
  count++;
  return f_max;
}
double fn_maxi(float af_A)
{
  static float        f_max  = 0;
  static int          maxi   = -1;
  static int          count  = 0;
  if(!count)
  {
    f_max = af_A;
  }
  if(af_A >= f_max)
  {
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
double fn_sum(float af_A)
{
  static double f_sum  = 0.;
  f_sum += af_A;
  return f_sum;
}

double fn2_sum(float af_A)
{
  static double f_sum  = 0.;
  f_sum += af_A;
  return f_sum;
}

// Additional 'sum' and 'prod' functions, but using a global vars
//+ as opposed to local static. The caller must make sure to zero the global
//+ var where appropriate.
double Gf_sum   = 0.;
double Gf_prod    = 1.;
double Gfn_sum(float af_A)
{
  Gf_sum += af_A;
  return Gf_sum;
}
double Gfn_prod(float af_A)
{
  Gf_prod *= af_A;
  return Gf_prod;
}

double fn_sum2(float af_A)
{
  static double f_sum2 = 0.;
  f_sum2  += af_A * af_A;
  return f_sum2;
}

double fn_mean(float af_A)
{
  static int  count   = 0;
  float       f_sum   = 0.;
  float       f_mean  = 0.;
  f_sum       = fn_sum(af_A);
  f_mean      = f_sum / ++count;
  return(f_mean);
}

double fn_dev(float af_A)
{
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
  char*   apch_inputFile,
  int*    ap_FSFILETYPE
);

// I/O functions
short CURV_arrayProgress_print(
  int     asize,
  int     acurrent,
  char*   apch_message
);

e_FILEACCESS
VOL_fileRead(
  char*   apch_VolFileName,
  unsigned long*    ap_vectorSize,
  float*    apf_volData[]
);

e_FILEACCESS
VOL_fileWrite(
  char*   apch_VolFileName,
  unsigned long     a_vectorSize,
  float*    apf_volData
);

e_FILEACCESS
CURV_fileRead(
  char*   apch_curvFileName,
  unsigned long*    ap_vectorSize,
  float*    apf_curv[]
);

e_FILEACCESS
CURV_fileWrite(
  char*   apch_curvFileName,
  unsigned long     a_vectorSize,
  float*    apf_curv
);

short   b_outCurvFile_write(e_operation e_op);
short   CURV_process(void);
short   CURV_functionRunABC( double (*F)(float f_A, float f_B) );
double  CURV_functionRunAC( double (*F)(float f_A) );
int PrintStatToFile(char *fname, double stat);

int main(int argc, char *argv[]) ;

// ----------------------------------------------------------------------------
// IMPLEMENTATION
// ----------------------------------------------------------------------------

#include "mris_calc.help.xml.h"
static void
synopsis_show(void)
{
  outputHelpXml(mris_calc_help_xml,
                mris_calc_help_xml_len);
  exit(1) ;
}

void
simpleSynopsis_show(void)
{
  char  pch_errorMessage[STRBUF];

  sprintf(pch_errorMessage,
          "Insufficient number of arguments.\n"
          "You should specify '<input1> <ACTION> [<input2> | <floatNumber>]'\n"
          "Use a '-u' for full usage instructions.");

  ErrorExit(10, "%s: %s", G_pch_progname, pch_errorMessage);
}

void *
xmalloc (size_t size)
{
  void *value = malloc (size);
  if (value == 0)
  {
    ErrorExit(10, "%s: virtual memory exhausted.", G_pch_progname);
  }
  return value;
}

void
verbosity_set(void)
{
  if(G_verbosity)
  {
    G_FP  = stdout;
    options_print();
  }
  else if((G_FP = fopen("/dev/null", "w")) == NULL)
    ErrorExit(  ERROR_NOFILE,
                "%s: Could not open /dev/null for console sink!",
                G_pch_progname);
}

void
init(void)
{
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
}

void
error_exit(
  const char*   apch_action,
  const char*   apch_error,
  int   exitCode)
{

  char  pch_errorMessage[STRBUF];
  strcpy(pch_errorMessage, "");

  sprintf(pch_errorMessage, "\n%s:\n"
  "\tSorry, but I seem to have encountered an error.\n"
  "\tWhile %s,\n"
  "\t%s\n",
  G_pch_progname, apch_action, apch_error);

  fprintf(stderr, "%s", pch_errorMessage);
  fprintf(stderr, "\n");
  exit(exitCode);
}

void
error_unmatchedSizes(void)
{
  error_exit( "checking on input filetype sizes",
              "I found a size mismatch, i.e. len(input1)!=len(input2)",
              10);
}

void
error_noVolumeStruct(void)
{
  error_exit( "checking on output volume",
              "it seems that internal volume structure is invalid.",
              10);
}

void
error_volumeWriteSizeMismatch(void)
{
  error_exit( "checking on output volume",
              "I found a size mismatch, internal data does not fit into vol.",
              10);
}

void
error_incompatibleInputFileTypes(void)
{
  error_exit( "checking on input filetypes",
              "it seems that you specified incompatible types.",
              10);
}

void
error_incompatibleOutputFileType(void)
{
  error_exit( "checking on output filetype",
              "it seems there is a type mismatch to the input type.",
              10);
}

void
error_zeroLengthInput(void)
{
  error_exit( "checking input 1",
              "it seems that the input file has no contents.",
              10);
}

void
error_surfacesNotHandled(void)
{
  error_exit( "checking inputs",
              "FreeSurfer surface files are not handled yet.",
              10);
}

void
error_unknownFileType(char* apch_inputFile)
{
  char  pch_action[65536];
  sprintf(pch_action, "accessing file '%s'", apch_inputFile);
  error_exit(   pch_action,
                "I was unable to identify this file type. Does it exist?",
                10);
}

void
error_primaryCurvsBackup(void)
{
  error_exit( "making backup of internal data arrays",
              "it seems that some of the backups already exist.",
              11);
}

void
error_primaryCurvsCompress(void)
{
  error_exit( "compressing primary internal data arrays about labels",
              "an internal memory error occurred.",
              11);
}

void
error_expandCurv3(void)
{
  error_exit( "expanding output data array",
              "an internal memory error occurred.",
              11);
}

void
output_init(void)
{
  unsigned long i;
  G_sizeCurv3   = G_sizeCurv1;
  G_pf_arrayCurv3 = (float*) xmalloc(G_sizeCurv1 * sizeof(float));
  for(i=0; i<G_sizeCurv3; i++)
  {
    G_pf_arrayCurv3[i] = 0.0;
  }
  if(!Gb_file3)
  {
    strcat(G_pch_curvFile3, Gppch_fileDotExt[G_eFILETYPE1]);
  }
}

e_FILETYPE
fileType_find(
  char*   apch_inputFile,
  int*    ap_FSFILETYPE)
{

  int       len, type, filestat;
  float       f           = -1.0;
  char*     pch_end;
  char      pch_floatArg[STRBUF];
  struct stat   stFileInfo;

  // The acph_inputFile can be either a file on the filesystem
  // or a float argument. Check if <apch_inputFile> is an actual
  // file, if not, check if it converts to a float.
//  filestat = stat(basename(apch_inputFile), &stFileInfo);
  filestat = stat((apch_inputFile), &stFileInfo);
  if(filestat)
  {
    // <apch_inputFile> does not seem to refer to a valid file.
    // Check if it converts to a float.
    strcpy(pch_floatArg, apch_inputFile);
    f = strtof(pch_floatArg, &pch_end);
    len = strlen(pch_end);
    if(!len)
    {
      return(e_FloatArg);
    }
  }

  // Check if input is a volume file...
//  type     = mri_identify(basename(apch_inputFile));
  type     = mri_identify((apch_inputFile));
  *ap_FSFILETYPE = type;
  if(type != MRI_VOLUME_TYPE_UNKNOWN && type != MRI_CURV_FILE)
  {
    return(e_VolumeFile);
  }

  // Check if input is a curvature file...
  if(type == MRI_CURV_FILE)
  {
    return(e_CurvatureFile);
  }

  if(type == -1)
  {
    return(e_Unknown);
  }

  // Assume that the input is a surface file...
  return(e_SurfaceFile);
}

void
fileIO_errorHander(
  char*   pch_filename,
  e_FILEACCESS eERROR)
{
  switch(eERROR)
  {
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
  char*     apch_fileName,
  unsigned long*      ap_vectorSize,
  float*      apf_curv[],
  e_FILETYPE*   aeFILETYPE,
  int*    a_FSFILETYPE)
{

  e_FILEACCESS  eACCESS   = e_UNSPECIFIED;

  *aeFILETYPE = fileType_find(apch_fileName, a_FSFILETYPE);
  switch(*aeFILETYPE)
  {
  case e_Unknown:
    error_unknownFileType(apch_fileName);
    break;
  case e_FloatArg:
    break;
  case e_VolumeFile:
    eACCESS = VOL_fileRead( apch_fileName, ap_vectorSize, apf_curv);
    break;
  case e_CurvatureFile:
    eACCESS = CURV_fileRead(apch_fileName, ap_vectorSize, apf_curv);
    break;
  case e_SurfaceFile:
    error_surfacesNotHandled();
    break;
  }
  return eACCESS;
}

e_FILEACCESS
fileWrite(
  char*     apch_fileName,
  unsigned long  a_vectorSize,
  float*      apf_curv
)
{

  // HISTORY
  // 09 March 2011
  // o Fixed output file name type handling. If no
  //   output extension is specified, append extension
  //   corresponding to input type.
  // o Correctly handle curvature file name extension spec.
  //

  e_FILEACCESS  eACCESS     = e_UNSPECIFIED;
  char*   pch_fileExt;
  int     b_fileID_OK = 0;

  pch_fileExt = fio_extension(apch_fileName);
  if(pch_fileExt)
  {
    if(!strcmp(pch_fileExt, Gppch_fileExt[e_CurvatureFile]))
    {
      G_eFILETYPE3  = e_CurvatureFile;
      G_FSFILETYPE3 = MRI_CURV_FILE;
      b_fileID_OK   = 1;
    }
    else
    {
      G_eFILETYPE3 = fileType_find(apch_fileName, &G_FSFILETYPE3);
      if(G_eFILETYPE3 == e_Unknown)
      {
        if(Gb_strictExtensions)
        {
          strcat(apch_fileName, Gppch_fileDotExt[G_eFILETYPE1]);
        }
      }
      else
      {
        if(strcmp(pch_fileExt, Gppch_fileExt[G_eFILETYPE3]) &&
            Gb_strictExtensions)
        {
          strcat(apch_fileName, Gppch_fileDotExt[G_eFILETYPE3]);
          G_eFILETYPE3  = G_eFILETYPE1;
        }
        b_fileID_OK = 1;
//        G_FSFILETYPE3 = G_FSFILETYPE1;
      }
    }
  }
  else
  {
    strcat(apch_fileName, Gppch_fileDotExt[G_eFILETYPE1]);
  }

  if(!b_fileID_OK)
  {
    G_eFILETYPE3  = G_eFILETYPE1;
    G_FSFILETYPE3 = G_FSFILETYPE1;
  }
  if(G_eFILETYPE3 != G_eFILETYPE1)
  {
    error_incompatibleOutputFileType();
  }
  switch(G_eFILETYPE1)
  {
  case e_Unknown:
    break;
  case e_FloatArg:
    break;
  case e_VolumeFile:
    eACCESS = VOL_fileWrite( apch_fileName, a_vectorSize, apf_curv);
    break;
  case e_CurvatureFile:
    eACCESS = CURV_fileWrite(apch_fileName, a_vectorSize, apf_curv);
    break;
  case e_SurfaceFile:
    error_surfacesNotHandled();
    break;
  }
  return eACCESS;
}

void
debuggingInfo_display()
{
  cprintd("Size of input1: ",   G_sizeCurv1);
  cprints("Type of input1: ",   Gppch_filetype[G_eFILETYPE1]);
  cprints("ACTION:",      G_pch_operator);
  cprintd("Size of input2: ",   G_sizeCurv2);
  cprints("Type of input2: ",   Gppch_filetype[G_eFILETYPE2]);
}

short
CURV_set(
  float*        apf_array[],
  unsigned long a_size,
  float         af_val)
{
  //
  // PRECONDITIONS
  // o apf_array does not need to have been allocated... if non-NULL will
  //   first free.
  //
  // POSTCONDITIONS
  // o apf_array is set to have all values set to <af_val>.
  //
  unsigned long          i;
  float*        pf_array        = NULL;

//   if(*apf_array) free(*apf_array);
  pf_array = (float*) malloc(a_size * sizeof(float));
  for(i=0; i<a_size; i++)
  {
    pf_array[i] = af_val;
  }
  *apf_array    = pf_array;
  return 1;
}

short
CURV_copy(
  float         apf_arraySrc[],
  float         apf_arrayDst[],
  int           a_size
)
{
  //
  // PRECONDITIONS
  // o <a_size> must be correct for both array inputs.
  //
  // POSTCONDITIONS
  // o contents of <apf_arrayDst> are copied to <apf_arraySrc>, i.e.
  //   <apf_arraySrc> = <apf_arrayDst>
  //
  unsigned long   i;
  for(i=0; i<a_size; i++)
  {
    apf_arraySrc[i]     = apf_arrayDst[i];
  }
  return 1;
}

/*!
  \fn array_compressUsingMask(float* apf_input, int a_sizeInput, float* apf_output[], int a_sizeOutput)
  \brief Using the internal mask label, compress the input array to masked values
  \param apf_input The input source array to compress.
  \param a_sizeInput The size of the input source array.
  \param apf_output The output array containing the compression. This array is destroyed and rebuilt
  \param a_sizeOutput The size of the output source array.
  \return TRUE -> OK; FALSE -> Error
*/
short
array_compressUsingLabelMask(
  float*      apf_input,
  unsigned long         a_sizeInput,
  float*      apf_output[],
  unsigned long a_sizeOutput
)
{
  //
  // PRECONDITIONS
  // o Label mask array and size must be non-zero (i.e. defined)
  // o apf_input and apf_output must be allocated and non-NULL!
  //
  // POSTCONDITIONS
  // o apf_output is freed, and then re-allocated to mask size.
  // o apf_output will only contain values that are masked by the
  //   internal label array G_pf_labelMask.
  //

  short       ret             = 0;
  unsigned long         n               = 0;
  unsigned long         m               = 0;
  float*      pf_array        = NULL;

  if(!apf_input || !apf_output)
  {
    ret     = 0;
  }
  else
  {
    free(*apf_output);
    CURV_set(&pf_array, G_sizeLabel, 0.0);
    for(n = 0; n < a_sizeInput; n++)
    {
      if(G_pf_labelMask[n])
      {
        pf_array[m++] = apf_input[n];
      }
    }
    *apf_output     = pf_array;
    ret     = 1;
  }
  return ret;
}

/*!
  \fn short primaryCurvs_backup(void)
  \brief Make backup copies of the primary curvs
  \return TRUE -> OK; FALSE -> Error
*/
short
primaryCurvs_backup(void)
{
  //
  // PRECONDITIONS
  // o The data arrays pf_arrayCurvN (N in 1, 2, 3) exist.
  // o The copy arrays MUST initially be NULL.
  // o Gb_canWrite must have been evaluated!
  //
  // POSTCONDITIONS
  // o The contents of the primary data arrays are copied
  //   to pf_arrayCurvNCopy.
  //

  short ret   = 1;

  if(!G_pf_arrayCurv1Copy)
  {
    CURV_set(&G_pf_arrayCurv1Copy, G_sizeCurv1, 0.0);
    CURV_copy(G_pf_arrayCurv1Copy, G_pf_arrayCurv1, G_sizeCurv1);
    G_sizeCurv1Copy = G_sizeCurv1;
  }
  else
  {
    ret = 0;
  }
  if(Gb_curvFile2)
  {
    if(!G_pf_arrayCurv2Copy)
    {
      CURV_set(&G_pf_arrayCurv2Copy, G_sizeCurv2, 0.0);
      CURV_copy(G_pf_arrayCurv2Copy, G_pf_arrayCurv2, G_sizeCurv2);
      G_sizeCurv2Copy = G_sizeCurv2;
    }
    else
    {
      ret = 0;
    }
  }
  if(Gb_canWrite)
  {
    if(!G_pf_arrayCurv3Copy)
    {
      CURV_set(&G_pf_arrayCurv3Copy, G_sizeCurv3, 0.0);
      CURV_copy(G_pf_arrayCurv3Copy, G_pf_arrayCurv3, G_sizeCurv3);
      G_sizeCurv3Copy = G_sizeCurv3;
    }
    else
    {
      ret = 0;
    }
  }
  return ret;
}

/*!
  \fn short primaryCurvs_compress(void)
  \brief Compress the primary curvs
  \return TRUE -> OK; FALSE -> Error
*/
short
primaryCurvs_compress(void)
{
  //
  // PRECONDITIONS
  // o The data arrays pf_arrayCurvN (N in 1, 2, 3) exist.
  // o The backup arrays must exist.
  //
  // POSTCONDITIONS
  // o The contents of the primary data arrays are compressed about
  //   the labels.
  //

  short ret   = 1;

  if(!array_compressUsingLabelMask(G_pf_arrayCurv1Copy, G_sizeCurv1,
                                   &G_pf_arrayCurv1, G_sizeLabel))
  {
    ret = 0;
  }
  G_sizeCurv1 = G_sizeLabel;
  if(Gb_curvFile2)
  {
    if(!array_compressUsingLabelMask(G_pf_arrayCurv2Copy, G_sizeCurv2,
                                     &G_pf_arrayCurv2, G_sizeLabel))
    {
      ret = 0;
    }
    G_sizeCurv2 = G_sizeLabel;
  }
  if(Gb_canWrite)
  {
    if(!array_compressUsingLabelMask(G_pf_arrayCurv3Copy, G_sizeCurv3,
                                     &G_pf_arrayCurv3, G_sizeLabel))
    {
      ret = 0;
    }
    G_sizeCurv3 = G_sizeLabel;
  }
  return ret;
}

/*!
  \fn label_initialize(void)
  \brief Initializes internal label structure and label mask array.
  \return TRUE -> OK; FALSE -> Error
*/
short
label_initialize(void)
{
  //
  // PRECONDITIONS
  // o Gb_labelMask is true
  //
  // POSTCONDITIONS
  // o Reads the FreeSurfer label specified by G_pch_labelFile
  // o Populates the G_plabel
  // o Generates the G_pf_labelMask which has TRUE (1) values
  //   at each index corresponding to the label.
  // o Generates the G_pf_labelMaskIndex which contains the
  //   offset indexes into the original array of each label element.
  //

  short       ret             = 0;
  unsigned long n = 0;
  unsigned long i  = 0;

  G_sizeLabel = 0;
  if(!Gb_labelMask)
  {
    ret     = 0;
  }
  else
  {
    if(G_verbosity)
    {
      cprints("Reading label...", "");
    }
    G_plabel    = LabelRead(NULL, G_pch_labelFile) ;
    if(G_verbosity)
    {
      cprints("", "ok");
    }
    if (!G_plabel)
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s",
                G_pch_progname, G_pch_labelFile) ;
    CURV_set(&G_pf_labelMask, G_sizeCurv1, 0.0);
    for(n = 0; n < G_plabel->n_points; n++)
    {
      G_pf_labelMask[G_plabel->lv[n].vno] = 1.0;
    }
    for(n=0; n < G_sizeCurv1; n++)
      if(G_pf_labelMask[n])
      {
        G_sizeLabel++;
      }
    CURV_set(&G_pf_labelMaskIndex, G_sizeLabel, 0.0);
    i=0;
    for(n=0; n < G_sizeCurv1; n++)
      if(G_pf_labelMask[n])
      {
        G_pf_labelMaskIndex[i++] = n;
      }
    ret = 1;
  }
  return ret;
}

/*!
  \fn array_expandUsingLabelMask(float* apf_compressedIn, int a_sizeInput, float* apf_output[], int a_sizeOutput)
  \brief Using the outputTemplate, expand the compressedInput array.
  \param apf_compressedIn The input source array to compress.
  \param a_sizeInput The size of the input source array.
  \param apf_outputTemplate The output array containing the compression. This array is destroyed and rebuilt
  \param a_sizeOutput The size of the output source array.
  \return TRUE -> OK; FALSE -> Error
*/
short
array_expandUsingLabelMask(
  float*      apf_compressedIn[],
  unsigned long         a_sizeInput,
  float*      apf_outputTemplate,
  unsigned long         a_sizeOutput
)
{
  //
  // PRECONDITIONS
  // o Label mask array and size must be non-zero (i.e. defined)
  // o apf_compressedIn and apf_outputTemplate must be allocated and non-NULL!
  // o a_sizeInput corresponds to compressed size (i.e. labelSize)
  // o a_sizeOutput corresponds to uncompressed size (i.e. original size)
  //
  // POSTCONDITIONS
  // o apf_compressedIn is expanded to the same size as apf_outputTemplate
  //   with the contents of compressedIn written to labelIndexed positions.
  //

  short       ret             = 0;
  unsigned long  n               = 0;
  unsigned long  m               = 0;
  float*      pf_inputCopy    = NULL;
  float*      pf_compressed   = NULL;

  if(!apf_compressedIn || !apf_outputTemplate)
  {
    ret     = 0;
  }
  else
  {
    CURV_set(&pf_inputCopy, a_sizeInput, 0.0);
    CURV_copy(pf_inputCopy, *apf_compressedIn, a_sizeInput);
    free(*apf_compressedIn);
    CURV_set(&pf_compressed, a_sizeOutput, 0.0);
    CURV_copy(pf_compressed, apf_outputTemplate, a_sizeOutput);
    for(n=0; n<a_sizeInput; n++)
    {
      pf_compressed[(int)G_pf_labelMaskIndex[n]] = pf_inputCopy[m++];
    }
    *apf_compressedIn       = pf_compressed;
    ret     = 1;
  }
  return ret;
}

int
main(
  int   argc,
  char  *argv[]
)
{

  int     nargs;
  short   ret   = 0;
  float   f_curv2;
  e_FILEACCESS  eACCESS;

  init();
  nargs = handleVersionOption(argc, argv, "mris_calc");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname=argv[0];
  G_pch_progname = argv[0] ;
  strcpy(G_pch_curvFile3, "out");

  for (; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = options_parse(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  // Process command line surface options and operator
  if((argc != 4) && (argc != 3))
  {
    simpleSynopsis_show();
  }
  if(argc==4)
  {
    Gb_curvFile2              = 1;
  }
  G_pch_curvFile1                       = argv[1];
  G_pch_operator                        = argv[2];
  if(Gb_curvFile2)
  {
    G_pch_curvFile2      = argv[3];
  }
  verbosity_set();

  // Read in relevant input files -
  eACCESS = fileRead(G_pch_curvFile1,
                     &G_sizeCurv1,
                     &G_pf_arrayCurv1,
                     &G_eFILETYPE1,
                     &G_FSFILETYPE1);
  if(eACCESS != e_OK)
  {
    fileIO_errorHander(G_pch_curvFile1, eACCESS);
  }

  // Second input file is optional, and, if specified could in fact
  //+ denote a float value and not an actual file per se.
  if(Gb_curvFile2)
  {
    eACCESS   = fileRead(G_pch_curvFile2,
                         &G_sizeCurv2,
                         &G_pf_arrayCurv2,
                         &G_eFILETYPE2,
                         &G_FSFILETYPE2);
    if(G_eFILETYPE2 == e_FloatArg)
    {
      f_curv2     = atof(G_pch_curvFile2);
      G_sizeCurv2   = G_sizeCurv1;
      CURV_set(&G_pf_arrayCurv2, G_sizeCurv2, f_curv2);
      /*      G_pf_arrayCurv2 = (float*) malloc(G_sizeCurv1 * sizeof(float));
            for(i=0; i<G_sizeCurv2; i++)
                  G_pf_arrayCurv2[i] = f_curv2;*/
    }
  }
  if(G_verbosity)
  {
    debuggingInfo_display();
  }
  if(!G_sizeCurv1)
  {
    error_zeroLengthInput();
  }
  if(G_sizeCurv1 != G_sizeCurv2 && argc==4)
  {
    error_unmatchedSizes();
  }
  if(G_eFILETYPE1 != G_eFILETYPE2 && G_eFILETYPE2 != e_FloatArg && argc==4)
  {
    error_incompatibleInputFileTypes();
  }

  output_init();
  ret     = CURV_process();
  if(G_verbosity)
  {
    printf("\n");
  }
  if(!ret)
  {
    ErrorExit(ERROR_BADPARM, "Unknown operation %s failed.\n", G_pch_operator);
  }
  printf("#VMPC# mris_calc VmPeak  %d\n",GetVmPeak());

  return NO_ERROR;
}

static int
options_print(void)
{
  cprints("Input file 1:", G_pch_curvFile1);
  cprints("ACTION:", G_pch_operator);
  if(G_pch_curvFile2)
  {
    cprints("Input file 2:", G_pch_curvFile2);
  }
  cprints("Output file:", G_pch_curvFile3);

  return 1;
}

static int
options_parse(int argc, char *argv[])
{
  int    nargs    = 0;
  char*  option;
  const char*  pch_text;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-output") || (toupper(*option) == 'O'))
  {
    strcpy(G_pch_curvFile3, argv[2]);
    Gb_file3      = 1;
    nargs           = 1;
  }
  else if(!stricmp(option, "threads") || !stricmp(option, "nthreads") ){
    int nthreads;
    sscanf(argv[2],"%d",&nthreads);
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    printf("nthreads %d\n",nthreads);
    nargs = 1;
  } 

  else if (!stricmp(option, "-strictExtensions") || (toupper(*option) == 'E'))
  {
    Gb_strictExtensions = 1;
  }
  else if (!stricmp(option, "-version") || (toupper(*option) == 'V'))
  {
    version_print();
  }
  else if (!stricmp(option, "-verbosity"))
  {
    G_verbosity   = atoi(argv[2]);
    nargs         = 1;
  }
  else if (!stricmp(option, "-label") || (toupper(*option) == 'L'))
  {
    strcpy(G_pch_labelFile, argv[2]);
    Gb_labelMask        = 1;
    nargs               = 1;
  }
  else if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    synopsis_show() ;
    exit(1) ;
  }
  else switch (toupper(*option))
    {
    case 'T':
      pch_text  = "void";
      break ;
    case '?':
    case 'U':
    case 'H':
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
version_print(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

/*!
  \fn static e_operation operation_lookup(char* apch_operation)
  \brief Convert a string 'operation' to an enumerated e_operation
  \param apch_operation The mathematical operation to perform..
  \return The internal lookup enumerated e_operation type is returned.
*/
static e_operation
operation_lookup(
  char* apch_operation
)
{

  e_operation   e_op = e_unknown;

  if(     !strcmp(apch_operation, "mul"))
  {
    e_op    = e_mul;
  }
  else if(!strcmp(apch_operation, "div"))
  {
    e_op    = e_div;
  }
  else if(!strcmp(apch_operation, "mod"))
  {
    e_op    = e_mod;
  }
  else if(!strcmp(apch_operation, "add"))
  {
    e_op    = e_add;
  }
  else if(!strcmp(apch_operation, "pow"))
  {
    e_op    = e_pow;
  }
  else if(!strcmp(apch_operation, "sub"))
  {
    e_op    = e_sub;
  }
  else if(!strcmp(apch_operation, "sub0"))
  {
    e_op    = e_sub0;
  }
  else if(!strcmp(apch_operation, "sratio"))
  {
    e_op    = e_sratio;
  }
  else if(!strcmp(apch_operation, "pctdiff"))
  {
    e_op    = e_pctdiff;
  }
  else if(!strcmp(apch_operation, "pctdiff0"))
  {
    e_op    = e_pctdiff0;
  }
  else if(!strcmp(apch_operation, "sqd"))
  {
    e_op    = e_sqd;
  }
  else if(!strcmp(apch_operation, "sqr"))
  {
    e_op    = e_sqr;
  }
  else if(!strcmp(apch_operation, "sqrt"))
  {
    e_op    = e_sqrt;
  }
  else if(!strcmp(apch_operation, "not"))
  {
    e_op    = e_not;
  }
  else if(!strcmp(apch_operation, "set"))
  {
    e_op    = e_set;
  }
  else if(!strcmp(apch_operation, "atan2"))
  {
    e_op    = e_atan2;
  }
  else if(!strcmp(apch_operation, "bcor"))
  {
    e_op    = e_bcor;
  }
  else if(!strcmp(apch_operation, "sig2p"))
  {
    e_op    = e_sig2p;
  }
  else if(!strcmp(apch_operation, "mag"))
  {
    e_op    = e_mag;
  }
  else if(!strcmp(apch_operation, "log10"))
  {
    e_op    = e_log10;
  }
  else if(!strcmp(apch_operation, "abs"))
  {
    e_op    = e_abs;
  }
  else if(!strcmp(apch_operation, "inv"))
  {
    e_op    = e_inv;
  }
  else if(!strcmp(apch_operation, "sign"))
  {
    e_op    = e_sign;
  }

  else if(!strcmp(apch_operation, "eq"))
  {
    e_op    = e_eq;
  }
  else if(!strcmp(apch_operation, "lt"))
  {
    e_op    = e_lt;
  }
  else if(!strcmp(apch_operation, "lte"))
  {
    e_op    = e_lte;
  }
  else if(!strcmp(apch_operation, "gt"))
  {
    e_op    = e_gt;
  }
  else if(!strcmp(apch_operation, "gte"))
  {
    e_op    = e_gte;
  }
  else if(!strcmp(apch_operation, "masked"))
  {
    e_op    = e_masked;
  }
  else if(!strcmp(apch_operation, "and"))
  {
    e_op    = e_and;
  }
  else if(!strcmp(apch_operation, "or"))
  {
    e_op    = e_or;
  }
  else if(!strcmp(apch_operation, "andbw"))
  {
    e_op    = e_andbw;
  }
  else if(!strcmp(apch_operation, "orbw"))
  {
    e_op    = e_orbw;
  }
  else if(!strcmp(apch_operation, "upperlimit"))
  {
    e_op    = e_upperlimit;
  }
  else if(!strcmp(apch_operation, "lowerlimit"))
  {
    e_op    = e_lowerlimit;
  }

  else if(!strcmp(apch_operation, "stats"))
  {
    e_op    = e_stats;
  }
  else if(!strcmp(apch_operation, "size"))
  {
    e_op    = e_size;
  }
  else if(!strcmp(apch_operation, "min"))
  {
    e_op    = e_min;
  }
  else if(!strcmp(apch_operation, "max"))
  {
    e_op    = e_max;
  }
  else if(!strcmp(apch_operation, "mean"))
  {
    e_op    = e_mean;
  }
  else if(!strcmp(apch_operation, "std"))
  {
    e_op    = e_std;
  }
  else if(!strcmp(apch_operation, "sum"))
  {
    e_op    = e_sum;
  }
  else if(!strcmp(apch_operation, "prod"))
  {
    e_op    = e_prod;
  }
  else if(!strcmp(apch_operation, "norm"))
  {
    e_op    = e_norm;
  }
  else if(!strcmp(apch_operation, "log"))
  {
    e_op    = e_log;
  }

  else if(!strcmp(apch_operation, "ascii"))
  {
    e_op    = e_ascii;
  }

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
  char*   apch_volFileName,
  unsigned long*    ap_vectorSize,
  float*    apf_data[]
)
{

  //char      pch_readMessage[STRBUF];
  int     i, j, k, f;
  unsigned long I = 0;
  MRI*      pMRI            = NULL;
  float*      pf_data           = NULL;

  if(G_verbosity)
  {
    cprints("Reading...", "");
  }
  if( (pMRI = MRIread(apch_volFileName)) == NULL)
  {
    return e_READACCESSERROR;
  }
  if(G_verbosity)
  {
    cprints("", "ok");
  }
  *ap_vectorSize  = (unsigned long)pMRI->width*pMRI->height*pMRI->depth*pMRI->nframes;
  pf_data       = (float*) xmalloc(*ap_vectorSize * sizeof(float));
  if(!pf_data){
    printf("Alloc failed\n");
    exit(1);
  }

  printf("Packing\n"); fflush(stdout);
  for(f=0; f<pMRI->nframes; f++)        // number of frames
    for(i=0; i<pMRI->width; i++)  // 'x', i.e. columns in slice
      for(j=0; j<pMRI->height; j++) // 'y', i.e. rows in slice
        for(k=0; k<pMRI->depth; k++)    // 'z', i.e. # of slices
        {
          //CURV_arrayProgress_print(*ap_vectorSize, I, pch_readMessage);
          pf_data[I++]  = (float) MRIgetVoxVal(pMRI, i, j, k, f);
        }
  *apf_data = pf_data;
  printf("Done packing\n"); fflush(stdout);

  // Only need to keep the header around
  //Gp_MRI = pMRI;
  Gp_MRI = MRIcopyHeader(pMRI,NULL);   // Global pointer.
  MRIfree(&pMRI); 
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
  char*   apch_volFileName,
  unsigned long  a_vectorSize,
  float*    apf_data
)
{
  //
  // PRECONDITIONS
  // o The Gp_MRI *must* have been set with a previous call to VOL_fileRead.
  //
  // POSTCONDITIONS
  // o Gp_MRI will have its voxel data set element-by-element to apf_data.
  // o Gp_MRI saved to <apchq_volFileName>.
  //

  unsigned long  volSize;
  int     i, j, k, f;
  unsigned long I=0;
  char    pch_readMessage[STRBUF];
  int     ret;
  e_FILEACCESS retfa;
  MRI     *out;

  if(!Gp_MRI)
  {
    error_noVolumeStruct();
  }
  volSize = (unsigned long)Gp_MRI->width*Gp_MRI->height*Gp_MRI->depth*Gp_MRI->nframes;
  if(volSize != a_vectorSize)
  {
    error_volumeWriteSizeMismatch();
  }

  // Free before MRI alloc to reduce max
  if(G_pf_arrayCurv1) free(G_pf_arrayCurv1);
  if(G_pf_arrayCurv2) free(G_pf_arrayCurv2);
  // this is where the output is stored, so don't free
  //if(G_pf_arrayCurv3) free(G_pf_arrayCurv3);

  sprintf(pch_readMessage, "Repacking %s", apch_volFileName);
  fflush(stdout);  fflush(stderr);
  out = MRIallocSequence(Gp_MRI->width, Gp_MRI->height, Gp_MRI->depth, MRI_FLOAT, Gp_MRI->nframes);
  MRIcopyHeader(Gp_MRI,out);
  MRIcopyPulseParameters(Gp_MRI,out);

  printf("Unpacking\n");  fflush(stdout);
  for(f=0; f<Gp_MRI->nframes; f++)              // number of frames
    for(i=0; i<Gp_MRI->width; i++)        // 'x', i.e. columns in slice
      for(j=0; j<Gp_MRI->height; j++)     // 'y', i.e. rows in slice
        for(k=0; k<Gp_MRI->depth; k++)      // 'z', i.e. # of slices
        {
          //CURV_arrayProgress_print(a_vectorSize, I, pch_readMessage);
          MRIsetVoxVal(out, i, j, k, f, (float) apf_data[I++]);
        }

  sprintf(pch_readMessage, "Saving result to '%s' (type = %s )",
          apch_volFileName, type_to_string (G_FSFILETYPE3));
  cprints(pch_readMessage, "");
  fflush(stdout);  fflush(stderr);
  ret = MRIwrite(out, apch_volFileName);
  if(!ret)
  {
    retfa = e_OK;
    cprints("", "ok");
  }
  else
  {
    retfa = e_WRITEACCESSERROR;
    cprints("", "error");
  }
  MRIfree(&out);
  return(retfa);
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
  char*   apch_curvFileName,
  unsigned long*    ap_vectorSize,
  float*    apf_curv[]
)
{

  FILE*     FP_curv;
  int       vnum;
  int       nvertices;
  unsigned long   i;
  char      pch_readMessage[STRBUF];
  float*      pf_data           = NULL;

  if((FP_curv = fopen(apch_curvFileName, "r")) == NULL)
  {
    return(e_READACCESSERROR);
  }
  fread3(&vnum, FP_curv);
  if(vnum == NEW_VERSION_MAGIC_NUMBER)
  {
    nvertices = freadInt(FP_curv);
    G_nfaces  = freadInt(FP_curv);
    G_valsPerVertex = freadInt(FP_curv);
    *ap_vectorSize  = nvertices;
    pf_data   = (float*) xmalloc(nvertices * sizeof(float));
    sprintf(pch_readMessage, "Reading %s", apch_curvFileName);
    for(i=0; i<nvertices; i++)
    {
      CURV_arrayProgress_print(nvertices, i, pch_readMessage);
      pf_data[i]  = freadFloat(FP_curv);
    }
  }
  else
  {
    return(e_WRONGMAGICNUMBER);
  }
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
  char*   apch_fileName,
  float*    apf_data
)
{
  FILE* FP_curv;
  unsigned long i;
  char  pch_readMessage[STRBUF];

  if((FP_curv = fopen(apch_fileName, "w")) == NULL)
  {
    return(e_WRITEACCESSERROR);
  }
  sprintf(pch_readMessage, "Writing %s", apch_fileName);
  sprintf(pch_readMessage, "Saving result to '%s' (type = ASCII)",
          apch_fileName);
  cprints(pch_readMessage, "");
  for(i=0; i<G_sizeCurv1; i++)
  {
    CURV_arrayProgress_print(G_sizeCurv1, i, pch_readMessage);
    fprintf(FP_curv, "%f\n", apf_data[i]);
  }
  cprints("", "ok");
  fclose(FP_curv);
  return(e_OK);
}


/*!
  \fn CURV_fileWrite(char* apch_curvFileName, int* a_vectorSize, float* apf_curv)
  \brief Write a FreeSurfer curvature array to a file
  \param apch_curvFileName The name of a FreeSurfer curvature file.
  \param a_vectorSize Size (i.e. number of elements) of data array.
  \param apf_curv Array containing the curvature values.
  \return If curvature file is successfully written, return e_OK, else return e_WRITEACCESSERROR.
*/
e_FILEACCESS
CURV_fileWrite(
  char*   apch_curvFileName,
  unsigned long a_vectorSize,
  float*    apf_curv
)
{
  FILE* FP_curv;
  unsigned long i;
  char  pch_readMessage[STRBUF];

  if((FP_curv = fopen(apch_curvFileName, "w")) == NULL)
  {
    return(e_WRITEACCESSERROR);
  }
  fwrite3(NEW_VERSION_MAGIC_NUMBER, FP_curv);
  fwriteInt((int)a_vectorSize, FP_curv);
  fwriteInt(G_nfaces, FP_curv);
  fwriteInt(G_valsPerVertex, FP_curv);
  sprintf(pch_readMessage, "Writing %s", apch_curvFileName);
  sprintf(pch_readMessage, "Saving result to '%s' (type = MRI_CURV_FILE)",
          apch_curvFileName);
  cprints(pch_readMessage, "");
  for(i=0; i<a_vectorSize; i++)
  {
    CURV_arrayProgress_print(a_vectorSize, i, pch_readMessage);
    fwriteFloat(apf_curv[i], FP_curv);
  }
  cprints("", "ok");
  fclose(FP_curv);
  return(e_OK);
}

short
CURV_arrayProgress_print(
  int   asize,
  int   acurrent,
  char* apch_message
)
{
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

  if(!acurrent)
  {
    if(apch_message != NULL)
    {
      fprintf(G_FP, "%*s", G_LC, apch_message);
    }
    fprintf(G_FP, " [");
    fflush(G_FP);
  }
  else if(fivePerc>0)
  {
    if (acurrent%fivePerc == fivePerc-1)
    {
      fprintf(G_FP, "#");
      fflush(G_FP);
    }
  }

  if(acurrent == asize-1)
  {
    fprintf(G_FP, "] ");
    if(apch_message != NULL)
    {
      fprintf(G_FP, "%*s\n", 1, "[ ok ]");
    }
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
  //  e_op    in    Action to perform
  //
  // POSTCONDITIONS
  //  Internally an output curvature data structure is always generated.
  //  Saving this output file to disk is only meaningful for a given
  //  set of actions -- usually the typical mathematical operators.
  //  Other actions don't require an output file to be saved, e.g the
  //  'stats' action.
  //

  short b_ret = 0;

  if(
    e_op == e_mul           ||
    e_op == e_div           ||
    e_op == e_mod           ||
    e_op == e_pow           ||
    e_op == e_add           ||
    e_op == e_sub           ||
    e_op == e_sub0          ||
    e_op == e_sratio        ||
    e_op == e_pctdiff       ||
    e_op == e_pctdiff0      ||
    e_op == e_sqd           ||
    e_op == e_sqr           ||
    e_op == e_sqrt          ||
    e_op == e_not           ||
    e_op == e_set           ||
    e_op == e_atan2         ||
    e_op == e_bcor          ||
    e_op == e_sig2p         ||
    e_op == e_mag           ||
    e_op == e_log10         ||
    e_op == e_abs           ||
    e_op == e_inv           ||
    e_op == e_sign          ||
    e_op == e_eq            ||
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
    e_op == e_norm          ||
    e_op == e_log
  )
  {
    b_ret = 1;
  }

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
  //  o If the operation to perform is unknown, then return zero,
  //    else return 1.
  //

  float f_min           = 0.;
  float f_max           = 0.;
  float f_range         = 0.;
  int   mini            = -1;
  int   maxi            = -1;
  float f_mean      = 0.;
  float f_std       = 0.;
  float f_dev       = 0.;
  char  pch_text[STRBUF];

  Ge_operation  = operation_lookup(G_pch_operator);
  Gb_canWrite   = b_outCurvFile_write(Ge_operation);

  if(Gb_labelMask)
  {
    // Backup the primary curvs
    if(!primaryCurvs_backup())
    {
      error_primaryCurvsBackup();
    }
    // Read in and initialize the labels
    label_initialize();
    // Compress/collapse the input arrays about the label indices
    if(!primaryCurvs_compress())
    {
      error_primaryCurvsCompress();
    }
  }

  switch(Ge_operation)
  {
  case  e_mul:
    CURV_functionRunABC(fn_mul);
    break;
  case  e_div:
    CURV_functionRunABC(fn_div);
    break;
  case  e_mod:
    CURV_functionRunABC(fn_mod);
    break;
  case  e_pow:
    CURV_functionRunABC(fn_pow);
    break;
  case  e_add:
    CURV_functionRunABC(fn_add);
    break;
  case  e_sub:
    CURV_functionRunABC(fn_sub);
    break;
  case  e_sub0:
    CURV_functionRunABC(fn_sub0);
    break;
  case  e_sratio:
    CURV_functionRunABC(fn_sratio);
    break;
  case  e_pctdiff:
    CURV_functionRunABC(fn_pctdiff);
    break;
  case  e_pctdiff0:
    CURV_functionRunABC(fn_pctdiff0);
    break;
  case  e_set:
    CURV_functionRunABC(fn_set);
    break;
  case  e_atan2:
    CURV_functionRunABC(fn_atan2);
    break;
  case  e_bcor:
    CURV_functionRunABC(fn_bcor);
    break;
  case  e_sig2p:
    CURV_functionRunAC(fn_sig2p);
    break;
  case  e_mag:
    CURV_functionRunABC(fn_mag);
    break;
  case  e_sqd:
    CURV_functionRunABC(fn_sqd);
    break;
  case  e_abs:
    CURV_functionRunAC( fn_abs);
    break;
  case  e_log10:
    CURV_functionRunAC( fn_log10);
    break;
  case  e_inv:
    CURV_functionRunAC( fn_inv);
    break;
  case  e_sign:
    CURV_functionRunAC( fn_sign);
    break;
  case  e_sqr:
    CURV_functionRunAC( fn_sqr);
    break;
  case  e_sqrt:
    CURV_functionRunAC( fn_sqrt);
    break;
  case  e_not:
    CURV_functionRunAC( fn_not);
    break;
  case  e_eq:
    CURV_functionRunABC(fn_eq);
    break;
  case  e_lt:
    CURV_functionRunABC(fn_lt);
    break;
  case  e_lte:
    CURV_functionRunABC(fn_lte);
    break;
  case  e_gt:
    CURV_functionRunABC(fn_gt);
    break;
  case  e_gte:
    CURV_functionRunABC(fn_gte);
    break;
  case  e_masked:
    CURV_functionRunABC(fn_masked);
    break;
  case  e_and:
    CURV_functionRunABC(fn_and);
    break;
  case  e_or:
    CURV_functionRunABC(fn_or);
    break;
  case  e_andbw:
    CURV_functionRunABC(fn_andbw);
    break;
  case  e_orbw:
    CURV_functionRunABC(fn_orbw);
    break;
  case  e_upperlimit:
    CURV_functionRunABC(fn_upl);
    break;
  case  e_lowerlimit:
    CURV_functionRunABC(fn_lrl);
    break;
  case  e_ascii:
    break;
  case  e_stats:
    break;
  case  e_size:
    break;
  case  e_min:
    break;
  case  e_mini:
    break;
  case  e_maxi:
    break;
  case  e_max:
    break;
  case  e_mean:
    break;
  case  e_std:
    break;
  case  e_sum:
    break;
  case  e_prod:
    break;
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
  case  e_log:
    CURV_functionRunAC(fn_log);
    break;
  case  e_unknown:
    // The operand is unknown
    // so return to caller with zero (error)
    return 0;
  }

  if(Ge_operation == e_set)
  {
    strcpy(G_pch_curvFile3, G_pch_curvFile1);
  }

  if(Ge_operation == e_size || Ge_operation == e_stats)
  {
    cprintd("Size", G_sizeCurv1);
  }

  if(Ge_operation == e_min || Ge_operation == e_stats)
  {
    f_min       = CURV_functionRunAC(fn_min);
    mini        = (int) CURV_functionRunAC(fn_mini);
    sprintf(pch_text, "%f (%d)", f_min, mini);
    cprints("Min@(index)", pch_text);
    if(Gb_file3)
    {
      PrintStatToFile(G_pch_curvFile3,f_min);
    }
  }
  if(Ge_operation == e_max || Ge_operation == e_stats)
  {
    f_max       = CURV_functionRunAC(fn_max);
    maxi        = (int) CURV_functionRunAC(fn_maxi);
    sprintf(pch_text, "%f (%d)", f_max, maxi);
    cprints("Max@(index)", pch_text);
    if(Gb_file3)
    {
      PrintStatToFile(G_pch_curvFile3,f_max);
    }
  }
  if(Ge_operation == e_mean || Ge_operation == e_stats)
  {
    f_mean      = CURV_functionRunAC(fn_mean);
    cprintf("Mean", f_mean);
    if(Gb_file3)
    {
      PrintStatToFile(G_pch_curvFile3,f_mean);
    }
  }
  if(Ge_operation == e_std || Ge_operation == e_stats)
  {
    f_dev       = CURV_functionRunAC(fn_dev);
    f_std       = sqrt(f_dev/(G_sizeCurv1-1));
    cprintf("Std", f_std);
    if(Gb_file3)
    {
      PrintStatToFile(G_pch_curvFile3,f_std);
    }
  }
  if(Ge_operation == e_sum || Ge_operation == e_stats)
  {
    Gf_sum  = 0.0;
    Gf_sum  = CURV_functionRunAC(Gfn_sum);
    cprintf("Sum", Gf_sum);
    if(Gb_file3)
    {
      PrintStatToFile(G_pch_curvFile3,Gf_sum);
    }
    int negCount        = CURV_functionRunAC(fn_negCount);
    int zeroCount       = CURV_functionRunAC(fn_zeroCount);
    int posCount        = CURV_functionRunAC(fn_posCount);
    cprintf("negCount",         negCount);
    cprintf("zeroCount",        zeroCount);
    cprintf("posCount",         posCount);
  }
  if(Ge_operation == e_prod || Ge_operation == e_stats)
  {
    Gf_prod = 1.0;
    Gf_prod = CURV_functionRunAC(Gfn_prod);
    cprintf("Prod", Gf_prod);
    if(Gb_file3)
    {
      PrintStatToFile(G_pch_curvFile3,Gf_prod);
    }
  }

  if(Ge_operation == e_ascii)
  {
    if(!Gb_file3)
    {
      strcat(G_pch_curvFile3, ".ascii");
    }

    if((ascii_fileWrite(G_pch_curvFile3, G_pf_arrayCurv1))==e_WRITEACCESSERROR)
    {
      printf("Write error: %s", G_pch_curvFile3);
    }
  }

  if(Gb_labelMask)
  {
    // Expand any output data in preparation for saving...
    if(Gb_canWrite)
      if(!array_expandUsingLabelMask( &G_pf_arrayCurv3, G_sizeLabel,
                                      G_pf_arrayCurv3Copy, G_sizeCurv3Copy))
      {
        error_expandCurv3();
      }
    G_sizeCurv3     = G_sizeCurv3Copy;
  }
  if(Gb_canWrite)
  {
    fileWrite(G_pch_curvFile3, G_sizeCurv3, G_pf_arrayCurv3);
  }

  return 1;
}

/*!
  \fn int PrintStatToFile(char *fname, double stat)
  \brief Simply prints a number to a file.
  \param fname - file name
  \param stat - value to print
  \return 0 if successful, 1 if failed
*/
int PrintStatToFile(char *fname, double stat)
{
  FILE *fp;
  fp = fopen(fname,"w");
  if(fp == NULL)
  {
    printf("ERROR: opening %s\n",fname);
    return(1);
  }
  fprintf(fp,"%lf\n",stat);
  fclose(fp);
  return(0);
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
  unsigned long i;

  ROMP_PF_begin
  #ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) 
  #endif
  for(i=0; i<G_sizeCurv1; i++)
  {
    ROMP_PFLB_begin
    double f_a = 0.;
    double f_b = 0.;
    double f_c = 0.;
    f_a                 = G_pf_arrayCurv1[i];
    f_b                 = G_pf_arrayCurv2[i];
    f_c                 = (F)(f_a, f_b);
    G_pf_arrayCurv3[i]  = f_c;
    ROMP_PFLB_end
  }
  ROMP_PF_end
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
  unsigned long i;
  double f_a = 0.;
  double f_c = 0.;

  for(i=0; i<G_sizeCurv1; i++)
  {
    f_a                 = G_pf_arrayCurv1[i];
    f_c                = (F)(f_a);
    G_pf_arrayCurv3[i]  = f_c;
  }
  return f_c;
}
