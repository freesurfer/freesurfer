/**
 * @brief Determine the mean and std of curvature files, and much much more.
 *
 * 'mris_curvature_stats' calculates statistics on curvature across a surface.
 * It also creates simple histograms on curvature distribution.
 * Additionally, it can also be used to analyze a surface and determine the
 * principal curvatures and save to extra curvature files, e.g. K, H, k1, k2.
 *
 */
/*
 * Original Author: Bruce Fischl / heavily hacked by Rudolph Pienaar
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
//#include "xDebug.h"
#include "label.h"

#define  STRBUF   65536
#define  MAX_FILES    1000
#define  CO( x )    fprintf(stdout, ( x ))
#define  CE( x )    fprintf(stderr, ( x ))
#define  START_i    3

// Calculations performed on the curvature surface
typedef enum _secondOrderType {
  e_Raw     = 0,  // "Raw" (native) curvature - no calcuation
  e_Gaussian    = 1,  // Gaussian curvature   = k1*k2
  e_Mean      = 2,  // Mean curvature = 0.5*(k1+k2))
  e_K1      = 3,  // k1 curvature
  e_K2      = 4,  // k2 curvature
  e_S     = 5,  // "sharpness"    = (k1-k2)^2
  e_C     = 6,  // "curvedness"   = sqrt(0.5*(k1^2+k2^2))
  e_SI      = 7,  // "shape index"  = atan((k1+k2)/(k2-k1))
  e_BE      = 8,  // "bending energy" = k1^2 + k2^2
  e_Normal    = 9,  // Normalised curvature
  e_Scaled    = 10,   // "Raw" scaled curvature
  e_ScaledTrans   = 11,   // Scaled and translated curvature
  e_FI      = 12  // Folding Index
} e_secondOrderType;

typedef enum _surfaceIntegrals {
  e_natural   = 0,  // "Natural" (native) integral - no conditions
  e_rectified   = 1,  // curvature is first rectified
  e_pos     = 2,  // only positive curvatures are considered
  e_neg     = 3     // only negative curvatures are considered
} e_surfaceIntegral;

// Set of possible output curvature files
const char* Gppch[] = {
  "Raw",
  "K",
  "H",
  "k1",
  "k2",
  "S",
  "C",
  "SI",
  "BE",
  "Normalised",
  "Scaled",
  "ScaledTrans",
  "FI"
};

// Output file prefixes
typedef enum _OFSP {
  e_None,         // No prefix
  e_Partial,        // A partial prefix (does not include the stem)
  e_Full,         // Full prefix - including stem
}   e_OFSP;

// This structure is used to house data on the min/max curvature of a
//  particular type (raw, K, H, etc.) across a surface. It is also
//  used to convey information between the various functions that
//  use it.
typedef struct _minMax {
  float f_min;      // The minimum curvature value
  float f_max;      // The maximum curvature value
  int vertexMin;    // The vertex where min lives
  int vertexMax;    // The vertex where max lives
  short b_minTest;    // If true, perform explicit min test
  short b_maxTest;    // If true, perform explicit max test
  short b_minViolation;   // If true, min value != explicit min
  short b_maxViolation;     // If true, max value != explicit max
} s_MINMAX;


int   main(int argc, char *argv[]) ;

static int   get_option(int argc, char *argv[]) ;
static void  usage_exit(void) ;
static void  print_usage(void) ;
static void  print_help(void) ;
static void  print_version(void) ;

//----------------------------//
// Function Prototypes: START //
//----------------------------//

// Simple functions of principal curvatures
float   f_sharpnessCurvature(float af_k1, float af_k2)
{
  return ((af_k1 - af_k2)*(af_k1 - af_k2));
}
float   f_bendingEnergyCurvature(float af_k1, float af_k2)
{
  return (af_k1*af_k1 + af_k2*af_k2);
}
float   f_curvednessCurvature(float af_k1, float af_k2)
{
  return (sqrt(0.5*(af_k1*af_k1 + af_k2*af_k2)));
}
float   f_shapeIndexCurvature(float af_k1, float af_k2)
{
  return (af_k1 == af_k2 ? 0 : atan((af_k1+af_k2)/(af_k2 - af_k1)));
}
float   f_foldingIndexCurvature(float af_k1, float af_k2)
{
  return (fabs(af_k1)*(fabs(af_k1) - fabs(af_k2)));
}

// Simple functions on specific vertices
float   f_absCurv(VERTEX*   pv)
{
  return (fabs(pv->curv));
}
float   f_pass(VERTEX*      pv)
{
  return (pv->curv);
}
short   f_allVertices(VERTEX*   pv)
{
  return 1;
}
short   f_greaterEqualZero(VERTEX*  pv)
{
  return (pv->curv >= 0 ? 1 : 0);
}
short   f_lessThanZero(VERTEX*    pv)
{
  return (pv->curv <0 ? 1 : 0);
}

int   MRISusePrincipalCurvatureFunction(
  MRI_SURFACE*    pmris,
  float       (*f)(float k1, float k2)
);

short   MRIS_surfaceIntegral_compute(
  MRI_SURFACE*    pmris,
  int*      p_verticesCounted,
  float*      pf_areaCounted,
  short     (*fcond)(VERTEX*  pv),
  float     (*fv) (VERTEX*  pv),
  float*      pf_surfaceIntegral,
  float*      pf_meanIntegral,
  float*      pf_areaNormalizedIntegral
);

int comp(
  const void *p,
  const void *q
);          // Compare p and q for qsort()

void  histogram_wrapper(
  MRIS*     amris,
  e_secondOrderType aesot
);

void  histogram_create(
  MRIS*       amris_curvature,
  float       af_minCurv,
  double      af_binSize,
  int         abins,
  float*      apf_histogram,
  e_secondOrderType   aesot
);

void  OFSP_create(
  char*       apch_prefix,
  char*       apch_suffix,
  e_OFSP      ae_OFSP
);

void  outputFileNames_create(void);
void  outputFiles_open(void);
void  outputFiles_close(void);

int MRIS_curvatureStats_analyze(
  MRIS*       apmris,
  e_secondOrderType   aesot
);

int   MRISminMaxCurvaturesSearchSOT(
  MRI_SURFACE*      apmris,
  int*        ap_vertexMin,
  int*        ap_vertexMax,
  float*        apf_min,
  float*        apf_max,
  e_secondOrderType     aesot
);

int   MRISminMaxCurvaturesSearch(
  MRI_SURFACE*      apmris,
  int*        ap_vertexMin,
  int*        ap_vertexMax,
  float*        apf_min,
  float*        apf_max
)
{
  MRISminMaxCurvaturesSearchSOT(  apmris,
                                  ap_vertexMin,
                                  ap_vertexMax,
                                  apf_min,
                                  apf_max,
                                  e_Raw);
  return(NO_ERROR);
};

int MRISvertexAreaPostProcess(
  MRI_SURFACE*    pmris
);

int MRIS_surfaceRipFlags_filter(
  MRI_SURFACE*  apmris,
  float*    apf_notRippedArea
);


//----------------------------//
// Function Prototypes: END   //
//----------------------------//


// Global variables

const char*  Progname ;
char*  hemi;

static int    navgs         = 0;
static int    normalize_flag      = 0;
static int    condition_no        = 0;
static int    stat_flag         = 0;
static char*  label_name        = NULL;
static char*  output_fname        = NULL;
static char*  curv_fname      = NULL;
static char   surf_name[STRBUF];
// Additional global variables (prefixed by 'G') added by RP.
// Flags are prefixed with Gb_ and are used to track
// user spec'd command line flags.
static int  Gb_shapeIndex     = 0;
static int  Gb_discreteCurvaturesUse  = 1;
static int    Gb_minMaxShow       = 0;
static int  Gb_histogram        = 0;
static int  Gb_histogramPercent     = 0;
static int  Gb_binSizeOverride    = 0;
static double   Gf_binSize        = 0.;
static int  Gb_histStartOverride    = 0;
static float  Gf_histStart        = 0.;
static int  Gb_histEndOverride    = 0;
static float  Gf_histEnd        = 0.;
static int  Gb_gaussianAndMean    = 0;
static int  Gb_output2File      = 0;
static int  Gb_scale        = 0;
static int  Gb_scaleMin       = 0;
static int  Gb_scaleMax       = 0;
static int  Gb_zeroVertex       = 0;
static int  G_zeroVertex        = 0;
static int    G_nbrs          = 2;
static int  G_bins        = 1;
static int  Gb_maxUlps        = 0;
static int  G_maxUlps       = 0;

static float  Gf_regionalTotalSurfaceArea = 0;
static int  G_regionalTotalVertexCount  = 0;
static int  Gb_regionalPercentages    = 0;
static int  Gb_vertexAreaNormalize    = 0;
static int  Gb_vertexAreaWeigh    = 0;
static int      Gb_vertexAreaNormalizeFrac      = 0;
static int      Gb_vertexAreaWeighFrac          = 0;
static int  Gb_postScale      = 0;
static float  Gf_postScale      = 0.;
static int  Gb_filter     = 0;
static int  Gb_lowPassFilter    = 0;
static float  Gf_lowPassFilter    = 0.;
static int  Gb_lowPassFilterGaussian  = 0;
static float  Gf_lowPassFilterGaussian  = 0.;
static int  Gb_highPassFilter   = 0;
static float  Gf_highPassFilter   = 0.;
static int  Gb_highPassFilterGaussian = 0;
static float  Gf_highPassFilterGaussian = 0.;

static short  Gb_signedPrincipals   = 0;
static char Gpch_filterLabel[STRBUF];
static short  Gb_filterLabel      = 0;

static float  Gf_foldingIndex     = 0.;
static float  Gf_intrinsicCurvaturePos  = 0.;
static float  Gf_intrinsicCurvatureNeg  = 0.;
static float  Gf_intrinsicCurvatureNat  = 0.;

// All possible output file name and suffixes
static  FILE* GpSTDOUT      = NULL;
static  short Gb_writeCurvatureFiles    = 0;
static  char  Gpch_log[STRBUF];
static  char  Gpch_logS[]       = "log";
static  FILE*   GpFILE_log        = NULL;
static  char  Gpch_rawHist[STRBUF];
static  char  Gpch_rawHistS[]     = "raw.hist";
static  FILE*   GpFILE_rawHist      = NULL;
static  char  Gpch_normCurv[STRBUF];
static  char  Gpch_normHist[STRBUF];
static  char  Gpch_normHistS[]    = "norm.hist";
static  FILE*   GpFILE_normHist     = NULL;
static  char  Gpch_normCurvS[]    = "norm.crv";
static  char  Gpch_KHist[STRBUF];
static  char  Gpch_KHistS[]       = "K.hist";
static  FILE*   GpFILE_KHist        = NULL;
static  char  Gpch_KCurv[STRBUF];
static  char  Gpch_KCurvS[]       = "K.crv";
static  char  Gpch_HHist[STRBUF];
static  char  Gpch_HHistS[]       = "H.hist";
static  FILE*   GpFILE_HHist        = NULL;
static  char  Gpch_HCurv[STRBUF];
static  char  Gpch_HCurvS[]       = "H.crv";
static  char  Gpch_scaledHist[STRBUF];
static  char  Gpch_scaledHistS[]    = "scaled.hist";
static  FILE*   GpFILE_scaledHist     = NULL;
static  char  Gpch_scaledCurv[STRBUF];
static  char  Gpch_scaledCurvS[]    = "scaled.crv";

static char   Gpch_K1Hist[STRBUF];
static char   Gpch_K1HistS[]      = "K1.hist";
static FILE*  GpFILE_K1Hist       = NULL;
static char   Gpch_K1Curv[STRBUF];
static char   Gpch_K1CurvS[]      = "K1.crv";

static char   Gpch_K2Hist[STRBUF];
static char   Gpch_K2HistS[]      = "K2.hist";
static FILE*  GpFILE_K2Hist       = NULL;
static char   Gpch_K2Curv[STRBUF];
static char   Gpch_K2CurvS[]      = "K2.crv";

static char   Gpch_SHist[STRBUF];
static char   Gpch_SHistS[]       = "S.hist";
static FILE*  GpFILE_SHist        = NULL;
static char   Gpch_SCurv[STRBUF];
static char   Gpch_SCurvS[]       = "S.crv";

static char   Gpch_CHist[STRBUF];
static char   Gpch_CHistS[]       = "C.hist";
static FILE*  GpFILE_CHist        = NULL;
static char   Gpch_CCurv[STRBUF];
static char   Gpch_CCurvS[]       = "C.crv";

static char   Gpch_BEHist[STRBUF];
static char   Gpch_BEHistS[]      = "BE.hist";
static FILE*  GpFILE_BEHist       = NULL;
static char   Gpch_BECurv[STRBUF];
static char   Gpch_BECurvS[]      = "BE.crv";

static char   Gpch_SIHist[STRBUF];
static char   Gpch_SIHistS[]      = "SI.hist";
static FILE*  GpFILE_SIHist       = NULL;
static char   Gpch_SICurv[STRBUF];
static char   Gpch_SICurvS[]      = "SI.crv";

static char   Gpch_FIHist[STRBUF];
static char   Gpch_FIHistS[]      = "FI.hist";
static FILE*  GpFILE_FIHist       = NULL;
static char   Gpch_FICurv[STRBUF];
static char   Gpch_FICurvS[]      = "FI.crv";

// These are used for tabular output
const int   G_leftCols        = 10;
const int   G_rightCols       = 40;

// Mean / sigma tracking and scaling
static double   Gpf_means[MAX_FILES] ;
static double   Gf_mean         = 0.;
static double   Gf_sigma        = 0.;
static double   Gf_n          = 0.;
static double   Gf_total        = 0.;
static double   Gf_total_sq       = 0.;
static double   Gf_scaleFactor      = 1.;
static double   Gf_scaleMin       = 0.;
static double   Gf_scaleMax       = 0.;

static int which_norm = NORM_MEAN;

char    Gpch_calc[1024];

int
stats_update(
  int ai
)
{
  //
  // PRECONDITIONS
  //  o Typically called if the mean and sigma for a surface
  //    has been determined.
  //  o <ai> denotes the position in the global 'mean' array
  //    to update.
  //
  // POSTCONDITIONS
  //  o Updates several globals with new mean/sigma data.
  //

  Gpf_means[ai-START_i]    = Gf_mean;
  Gf_total          += Gf_mean;
  Gf_total_sq       += Gf_mean*Gf_mean ;
  Gf_n        += 1.0;
  return 1;
}

void
surface_applyFilter(MRI_SURFACE*  apmris)
{
  cprints("Scanning surface for labels and/or filters...", "");
  G_regionalTotalVertexCount = MRIS_surfaceRipFlags_filter(
                                 apmris, &Gf_regionalTotalSurfaceArea);
  cprints("", "ok");
  cprintd("Filtered vertices",  G_regionalTotalVertexCount);
  cprintf("Filtered Area (mm^2)",   Gf_regionalTotalSurfaceArea);
}

int
LABEL_RipSurface(
  MRI_SURFACE*  apmris,
  LABEL*    aplbl
)
{
  //
  // PRECONDITIONS
  //  o Initialized surface <apmris>
  //  o Initialized label <aplbl>
  //
  // POSTCONDITIONS
  //  o Surface vertices that correspond to label positions
  //    will have their "ripflag" set false.
  //  o This function is a simplified version of LabelRipRestOfSurface(...)
  //    and specifically does not have any surface side effects other than
  //    setting the ripflag value.
  //

  int    vno, n ;
  VERTEX *v ;

  LabelToFlat(aplbl, apmris) ;
  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    v = &apmris->vertices[vno] ;
    v->ripflag = 1 ;
  }

  for (n = 0 ; n < aplbl->n_points ; n++) {
    vno = aplbl->lv[n].vno ;
    if (vno < 0 || vno >= apmris->nvertices) {
      continue ;
    }
    if (vno == Gdiag_no) {
      DiagBreak() ;
    }
    v = &apmris->vertices[vno] ;
    v->ripflag = 0 ;
  }
  return(NO_ERROR) ;
}

int
main(int argc, char *argv[])
{
  char          **av, fname[STRBUF], *sdir ;
  char          *subject_name;
  char    pch_surface[16384];
  char    pch_tmp[1024];
  int           ac, nargs;
  int   i       = START_i;
  MRI_SURFACE   *mris;
  int   b_surfaceFiltered   = 0;

  GpSTDOUT  = stdout;
//  InitDebugging( "mris_curvature_stats" );
  nargs = handleVersionOption(argc, argv, "mris_curvature_stats");
  if (nargs && argc - nargs == 1) {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  sdir = getenv("SUBJECTS_DIR") ;
  if (!sdir) {
    ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in environment.\n",Progname);
  }
  ac = argc ;
  av = argv ;
  strcpy(surf_name, "-x");
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  if(Gb_discreteCurvaturesUse) {
    strcpy(Gpch_calc, "discrete");
  } else {
    strcpy(Gpch_calc, "2nd order form");
  }

  if (!strcmp(surf_name, "-x")) {
    strcpy(surf_name, "smoothwm") ;
  }

  if (argc < 3) {
    usage_exit() ;
  }

  //
  // Once all the options have been processed, the command line is
  // mris_curvature_stats <str_subjectName> <str_hemi> [<str_curv0> ... <str_curvN>]
  //    0   1   2     3     ...     N
  //
  // The <str_curvN> are not necessary if the Gaussian curvatures are calculated.

  subject_name  = argv[1] ;
  hemi    = argv[2] ;
  int req = snprintf(fname, STRBUF, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, surf_name) ;
  if( req >= STRBUF ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  req = snprintf(pch_surface, 16384, "%s/%s.%s", subject_name, hemi, surf_name);
  if( req >= STRBUF ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  cprints("Setting surface", pch_surface);
  cprints("Reading surface...", "");
  mris = MRISread(fname) ;
  cprints("", "ok");
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;


  if (label_name) {
    cprints("Reading label...", "");
    LABEL *area ;
    area = LabelRead(subject_name, label_name) ;
    cprints("", "ok");
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s",
                Progname, label_name) ;
    LABEL_RipSurface(mris, area);
    LabelFree(&area) ;
  }

  outputFiles_close();
  outputFileNames_create();
  outputFiles_open();
  MRIScomputeMetricProperties(mris);

  // Process all the command-line spec'd curvature files
  if(argc>=4)
    for (Gf_n= Gf_total_sq = Gf_total = 0.0, i = START_i ; i < argc ; i++) {
      curv_fname = argv[i] ;
      cprints("Setting texture", curv_fname);
      cprints("Reading texture...", "");
      if (MRISreadCurvatureFile(mris, curv_fname) != NO_ERROR) {
        fprintf(stderr,
                "\n\t\t***WARNING!***\n");
        fprintf(stderr,
                "\tSome error has occurred while reading '%s.%s'.\n",
                hemi, curv_fname);
        fprintf(stderr,
                "\tThis might be due a vertex incompatibility\n");
        fprintf(stderr,
                "\tbetween the surface '%s.%s' and curv '%s.%s'.\n",
                hemi, surf_name, hemi, curv_fname);
        fprintf(stderr,
                "\n\tYou might be able to correct this by re-running\n");
        fprintf(stderr,
                "\t'mris_make_surfaces' on this dataset.\n\n");
        fprintf(stderr,
                "\tSkipping this (and any remaining) curvature files.\n");
        fprintf(stderr,
                "\tAny measurements / calcuations that depend on the\n");
        fprintf(stderr,
                "\tcurvature file will be skipped.\n");
        fprintf(stderr, "\t\t*************\n\n");
        break;
      }
      cprints("", "ok");

      if (Gb_zeroVertex) {
        MRISvertexCurvature_set(mris, G_zeroVertex, 0);
      }

      MRISaverageCurvatures(mris, navgs) ;

      if (Gb_scale) {
        MRISscaleCurvature(mris, Gf_scaleFactor);
        if (Gpch_scaledCurv[0] && Gb_writeCurvatureFiles) {
          MRISwriteCurvature(mris, Gpch_scaledCurv);
        }
      }

      printf("Gb_filter = %d\n", Gb_filter);

      if(Gb_filter || label_name) {
        surface_applyFilter(mris);
        b_surfaceFiltered = 1;
      }

      if (Gb_scaleMin && Gb_scaleMax) {
        MRISscaleCurvatures(mris, Gf_scaleMin, Gf_scaleMax);
        if (Gpch_scaledCurv[0] && Gb_writeCurvatureFiles) {
          MRISwriteCurvature(mris, Gpch_scaledCurv);
        }
      }

      if (normalize_flag) {
        MRISnormalizeCurvature(mris,which_norm);
        if (Gpch_normCurv[0] && Gb_writeCurvatureFiles) {
          MRISwriteCurvature(mris, Gpch_normCurv);
        }
      }
      MRIS_curvatureStats_analyze(mris, e_Raw);
      stats_update(i);
    }

  // Should we calculate all the principal curvature based curvatures?
  // This is a surface-based calculation, and does not depend on the
  // curvature processing loop - thus the 'curv' input is irrelevant
  if (Gb_gaussianAndMean) {
    MRISsetNeighborhoodSizeAndDist(mris, G_nbrs);
    if(!Gb_discreteCurvaturesUse) {
      cprints("Calculating Continuous Principal Curvatures...", "");
      MRIScomputeSecondFundamentalForm(mris);
      cprints("", "ok");
    } else {
      cprints("Calculating Discrete Principal Curvatures...", "");
      MRIScomputeSecondFundamentalFormDiscrete(mris, Gb_signedPrincipals);
    }

    // Apply any filters -- the pattern of ripMarks is used to limit
    // analysis to tagged vertices
    if((Gb_filter || label_name) && !b_surfaceFiltered) {
      surface_applyFilter(mris);
    }

    MRIS_curvatureStats_analyze(mris, e_Gaussian);
    stats_update(i++);
    if (Gpch_KCurv[0] && Gb_writeCurvatureFiles) {
      MRISwriteCurvature(mris,  Gpch_KCurv);
    }
    MRIS_curvatureStats_analyze(mris, e_Mean);
    stats_update(i++);
    if (Gpch_HCurv[0] && Gb_writeCurvatureFiles) {
      MRISwriteCurvature(mris,  Gpch_HCurv);
    }
    MRIS_curvatureStats_analyze(mris, e_K1);
    stats_update(i++);
    if (Gpch_K1Curv[0] && Gb_writeCurvatureFiles) {
      MRISwriteCurvature(mris,  Gpch_K1Curv);
    }
    MRIS_curvatureStats_analyze(mris, e_K2);
    stats_update(i++);
    if (Gpch_K2Curv[0] && Gb_writeCurvatureFiles) {
      MRISwriteCurvature(mris,  Gpch_K2Curv);
    }
    MRIS_curvatureStats_analyze(mris, e_S);
    stats_update(i++);
    if (Gpch_SCurv[0] && Gb_writeCurvatureFiles) {
      MRISwriteCurvature(mris,  Gpch_SCurv);
    }
    MRIS_curvatureStats_analyze(mris, e_C);
    stats_update(i++);
    if (Gpch_CCurv[0] && Gb_writeCurvatureFiles) {
      MRISwriteCurvature(mris,  Gpch_CCurv);
    }
    MRIS_curvatureStats_analyze(mris, e_BE);
    stats_update(i++);
    if (Gpch_BECurv[0] && Gb_writeCurvatureFiles) {
      MRISwriteCurvature(mris,  Gpch_BECurv);
    }
    MRIS_curvatureStats_analyze(mris, e_FI);
    stats_update(i++);
    if (Gpch_FICurv[0] && Gb_writeCurvatureFiles) {
      MRISwriteCurvature(mris,  Gpch_FICurv);
    }

    // NOTE:
    // The "Shape Index" can be problematic due to the atan calculations.
    // To analyze the "Shape Index", you must pass an explicit
    // '--shapeIndex' on the command line
    if(Gb_shapeIndex) {
      MRIS_curvatureStats_analyze(mris, e_SI);
      stats_update(i++);
      if (Gpch_SICurv[0] && Gb_writeCurvatureFiles) {
        MRISwriteCurvature(mris,  Gpch_SICurv);
      }
    }

    strcpy(pch_tmp, "");
    fprintf(GpSTDOUT, "\n");
    fprintf(GpSTDOUT, "%-55s%10s\n", "curv -- calculation type:", Gpch_calc);
    sprintf(pch_tmp, "curv -- Folding Index (FI):");
    fprintf(GpSTDOUT, "%-55s%10.5f\n", pch_tmp , Gf_foldingIndex);
    sprintf(pch_tmp, "curv -- Intrinsic Curvature Index - positive (ICIp):");
    fprintf(GpSTDOUT, "%-55s%10.5f\n", pch_tmp, Gf_intrinsicCurvaturePos);
    sprintf(pch_tmp, "curv -- Intrinsic Curvature Index - negative (ICIn):");
    fprintf(GpSTDOUT, "%-55s%10.5f\n", pch_tmp, Gf_intrinsicCurvatureNeg);
    sprintf(pch_tmp, "curv -- Intrinsic Curvature Index - natural  (ICIt):");
    fprintf(GpSTDOUT, "%-55s%10.5f\n",pch_tmp, Gf_intrinsicCurvatureNat);
  }

  if (Gf_n> 1.8) {
    Gf_mean   = Gf_total / Gf_n;
    Gf_sigma  = sqrt(Gf_total_sq/Gf_n- Gf_mean*Gf_mean) ;
    fprintf(GpSTDOUT, "\nMean across %d curvatures: %8.4e +- %8.4e\n",
            (int) Gf_n, Gf_mean, Gf_sigma) ;
  }

  MRISfree(&mris) ;
  fprintf(GpSTDOUT, "\n\n");
  outputFiles_close();
  exit(0) ;
  return(0) ;  /* for ansi */
}

short
surfaceIntegrals_compute(
  MRIS*     apmris,
  e_surfaceIntegral aeSI_type,
  float*      apf_surfaceIntegral,
  float*      apf_meanIntegral,
  int*      ap_vertices,
  float*      apf_areaNormalizedIntegral,
  float*      apf_area
)
{
  //
  // PRECONDITIONS
  // o amris must have been prepared with a call to
  //   MRIScomputeSecondFundamentalForm(...)
  // o ae_type defines the integral type to perform
  //
  // POSTCONDITIONS
  // o integral and vertex count are returned in pointers
  // o OK / not OK returned in function name
  //

  short b_ret = 0;

  switch(aeSI_type) {
  case e_natural:
    b_ret   = MRIS_surfaceIntegral_compute(
                apmris, ap_vertices, apf_area,
                f_allVertices,
                f_pass,
                apf_surfaceIntegral,
                apf_meanIntegral,
                apf_areaNormalizedIntegral);
    break;
  case e_rectified:
    b_ret   = MRIS_surfaceIntegral_compute(
                apmris, ap_vertices, apf_area,
                f_allVertices,
                f_absCurv,
                apf_surfaceIntegral,
                apf_meanIntegral,
                apf_areaNormalizedIntegral);
    break;
  case e_pos:
    b_ret   = MRIS_surfaceIntegral_compute(
                apmris, ap_vertices, apf_area,
                f_greaterEqualZero,
                f_absCurv,
                apf_surfaceIntegral,
                apf_meanIntegral,
                apf_areaNormalizedIntegral);
    break;
  case e_neg:
    b_ret   = MRIS_surfaceIntegral_compute(
                apmris, ap_vertices, apf_area,
                f_lessThanZero,
                f_absCurv,
                apf_surfaceIntegral,
                apf_meanIntegral,
                apf_areaNormalizedIntegral);
    break;
  }
  return b_ret;
}

int
MRIS_minMaxCurvature_analyze(
  MRIS*     apmris,
  e_secondOrderType aesot,
  s_MINMAX*   aps_minMax
)
{
  //
  // PRECONDITIONS
  //  o Typically called only by MRIS_minMaxCurve_report(...)
  //  o Either surfFunc or k1k2Func should be NULL -- this function
  //    will call whichever of the above is non-NULL.
  //  o The <s_MINMAX> structure is used to communicate back to the
  //    calling function. If the boolean b_{min,max}Violation fields are set
  //    to TRUE, the caller should interpret this to mean that an explicit
  //    lookup failed.
  //  o Similarly, if the <s_MINMAX> b_{min,max}Test fields are true as
  //    set by caller, this function will perform an explicit min/max
  //    search across the surface.
  //
  // POSTCONDITIONS
  //  o The <s_MINMAX> structure is populated, and pre- post- behaviour
  //    is set by the <s_MINMAX> boolean control flags.
  //

  int   vmin    = 0;
  int   vmax    = 0;
  float f_minExplicit = 0.;
  float f_maxExplicit = 0.;
//  char* pch_function  = "MRIS_minMaxCurvature_analyze";

//  DebugEnterFunction (( pch_function ));
  aps_minMax->b_minViolation  = 0;  // Assume no "violations" for min
  aps_minMax->b_maxViolation  = 0;  // Assume no "violations" for max
  MRISminMaxCurvaturesSearch( apmris,
                              &vmin,    &vmax,
                              &f_minExplicit, &f_maxExplicit);
  aps_minMax->vertexMin = vmin;
  aps_minMax->vertexMax = vmax;
  if (aps_minMax->b_minTest && aps_minMax->f_min != f_minExplicit) {
    fprintf(stderr, "\nWARN:%5s%-40s%f\n", Gppch[aesot],
            " lookup   min:", aps_minMax->f_min);
    fprintf(stderr, "WARN:%5s%-40s%f\tvertex = %d\n", Gppch[aesot],
            " explicit min:", f_minExplicit, vmin);
    aps_minMax->b_minViolation  = 1;
    aps_minMax->f_min     = f_minExplicit;
  }
  if (aps_minMax->b_maxTest && aps_minMax->f_max != f_maxExplicit) {
    fprintf(stderr, "\nWARN:%5s%-40s%f\n", Gppch[aesot],
            " lookup   max:", aps_minMax->f_max);
    fprintf(stderr, "WARN:%5s%-40s%f\tvertex = %d\n", Gppch[aesot],
            " explicit max:", f_maxExplicit, vmax);
    aps_minMax->b_maxViolation  = 1;
    aps_minMax->f_max     = f_maxExplicit;
  }
  return 1;
}

short
MRIS_curvatures_prepare(
  MRIS*     apmris,
  e_secondOrderType aesot
)
{
  //
  // PRECONDITIONS
  //  o For second order type curves, <aprmis> *must* have been
  //    pre-processed with a call to one of the second order
  //    functions:
  //
  //    MRIScomputeSecondFundamentalForm(...)
  //    MRIScomputeSecondFundamentalFormDiscrete(...)
  //
  //    Note that there is no way for this function to know that this
  //    preprocessing has occurred.
  //  o The <aesot> "second order type" determines which particular curvature
  //    value to process.
  //  o This is the main entry point for preparing a surface for subsequent
  //    analysis.
  //
  // POSTCONDITIONS
  //  o Depending on <aesot>, the <apmris> surface is prepared for subsequent
  //    analysis.
  //  o The return value from the surface prepation function is returned, or
  //    '-1' if no preparation was performed.
  //

//  char* pch_function  = "MRIS_curvatures_prepare";
  short ret   = -1;

//  DebugEnterFunction (( pch_function ));

  switch (aesot) {
  case e_Gaussian:
    ret = MRISuseGaussianCurvature(apmris);
    break;
  case e_Mean:
    ret = MRISuseMeanCurvature(apmris);
    break;
  case e_K1:
    ret = MRISuseK1Curvature(apmris);
    break;
  case e_K2:
    ret = MRISuseK2Curvature(apmris);
    break;
  case e_S:
    ret = MRISusePrincipalCurvatureFunction(apmris,
                                            f_sharpnessCurvature);
    break;
  case e_C:
    ret = MRISusePrincipalCurvatureFunction(apmris,
                                            f_curvednessCurvature);
    break;
  case e_BE:
    ret = MRISusePrincipalCurvatureFunction(apmris,
                                            f_bendingEnergyCurvature);
    break;
  case e_SI:
    ret = MRISusePrincipalCurvatureFunction(apmris,
                                            f_shapeIndexCurvature);
    break;
  case e_FI:
    ret = MRISusePrincipalCurvatureFunction(apmris,
                                            f_foldingIndexCurvature);
    break;
  case e_Raw:
  case e_Normal:
  case e_Scaled:
  case e_ScaledTrans:
    break;
  }

  // Check for possible normalizing / weighing factors
  if(   Gb_vertexAreaNormalize          ||
        Gb_vertexAreaWeigh              ||
        Gb_vertexAreaNormalizeFrac      ||
        Gb_vertexAreaWeighFrac) {
    MRISvertexAreaPostProcess(apmris);
  }
  return ret;
}

short
MRIS_minMaxCurve_report(
  MRIS*       apmris,
  e_secondOrderType   aesot,
  char*     apch_report
)
{
  //
  // PRECONDITIONS
  //  o MRIS_curvatures_prepare(...) *must* have been called for <aesot>.
  //  o The <aesot> "second order type" determines which particular curvature
  //    value to process.
  //  o This is the main entry point for analyzing the min/max curves. It is
  //    typically called by any higher-level function that wants to determine
  //    the min/max for a curvature function type (as specified in <aesot>).
  //
  // POSTCONDITIONS
  //  o Calls MRIS_minMaxCurvature_analyze(...) to find the min/max
  //    curvatures and vertex "occurrences".
  //  o Depending on <aesot>, data is printed to stdout (and
  //    output file, if spec'd).
  //  o Essentially, function determines the minimum and maximum values of
  //    a particular curvature by explicit search. It also checks, for
  //    curvature cases that record min/max in the <apmris> structure, that
  //    the embedded min/max matches the searched min/max.
  //  o The report itself is returned in the apch_report string.
  //

  s_MINMAX  s_minMax;
//  char* pch_function  = "MRIS_minMaxCurve_report";
  char  tmp[1024];

//  DebugEnterFunction (( pch_function ));

  s_minMax.f_min    = apmris->min_curv;
  s_minMax.f_max    = apmris->max_curv;
  s_minMax.vertexMin    = -1;
  s_minMax.vertexMax    = -1;
  s_minMax.b_minTest    = 1;
  s_minMax.b_maxTest    = 1;
  s_minMax.b_minViolation = 0;
  s_minMax.b_maxViolation = 0;

  if(aesot == e_Gaussian) {
    s_minMax.f_min = apmris->Kmin;
    s_minMax.f_max  = apmris->Kmax;
  }
  if(aesot == e_Mean) {
    s_minMax.f_min = apmris->Hmin;
    s_minMax.f_max  = apmris->Hmax;
  }

  MRIS_minMaxCurvature_analyze(apmris, aesot, &s_minMax);

  if(aesot == e_Gaussian) {
    if(s_minMax.b_minViolation) {
      apmris->Kmin = s_minMax.f_min;
    }
    if(s_minMax.b_maxViolation) {
      apmris->Kmax = s_minMax.f_max;
    }
  }
  if(aesot == e_Mean) {
    if(s_minMax.b_minViolation) {
      apmris->Hmin = s_minMax.f_min;
    }
    if(s_minMax.b_maxViolation) {
      apmris->Hmax = s_minMax.f_max;
    }
  }

  if(s_minMax.b_minViolation) {
    apmris->min_curv = s_minMax.f_min;
  }
  if(s_minMax.b_maxViolation) {
    apmris->max_curv = s_minMax.f_max;
  }

  strcpy(apch_report, "");
  sprintf(apch_report, "%*s%-*s",   G_leftCols,   Gppch[aesot],
          G_rightCols,  " Min:");
  sprintf(tmp, "%12.5f at vertex %-8d\n", s_minMax.f_min, s_minMax.vertexMin);
  strcat(apch_report,tmp);
  sprintf(tmp, "%*s%-*s",
          G_leftCols,   Gppch[aesot],
          G_rightCols,  " Max:");
  strcat(apch_report,tmp);
  sprintf(tmp, "%12.5f at vertex %-8d\n", s_minMax.f_max, s_minMax.vertexMax);
  strcat(apch_report,tmp);

  return 1;
}

short
MRIS_surfaceIntegrals_report(
  MRIS*       apmris,
  e_secondOrderType aesot,
  char*     apch_report
)
{
  //
  // PRECONDITIONS
  //  o The 'curv' field of <apmris> must contain the particular curvature
  //    map value to integrate.
  //  o The <apch_curvName> is a string that is prefixed to each output line
  //    denoting the curve being processed.
  //
  // POSTCONDITIONS
  //  o This function is the typical entry point to performing a surface
  //    integral. The 'report' in the function name is meant to convey
  //    that this function performs the integral, and *prints* the results
  //    to the output device.
  //  o The function performs/prints the following surface integral
  //    functions:
  //    - "natural":  no change/filtering on the vertex curv values
  //    - "abs":  the fabs(...) of each vertex curv value
  //    - "pos":  process only the positive vertex curv values
  //    - "neg":  process only the negative vertex curv values
  //  o In addition, the following qualifiers are also evaluated:
  //    - "Mean": the integral normalized to number of vertices
  //    - "AreaNorm": the integral value normalized to unit surface
  //  o The report itself is returned in the apch_report string.
  //

  char  pch_processedArea[65536];
  char  pch_misc[65536];
  // Surface integral variables
  float f_SInatural   = 0.0;
  float f_SIabs   = 0.0;
  float f_SIpos   = 0.0;
  float f_SIneg   = 0.0;
  float f_SInaturalMean = 0.0;
  int SInaturalVertices = 0;
  float f_SIabsMean = 0.0;
  int SIabsVertices = 0;
  float f_SIposMean = 0.0;
  int SIposVertices = 0;
  float f_SInegMean = 0.0;
  int SInegVertices = 0;

  float f_SInaturalAreaNorm = 0.0;
  float f_SInaturalArea = 0.0;
  float f_SIabsAreaNorm = 0.0;
  float f_SIabsArea = 0.0;
  float f_SIposAreaNorm = 0.0;
  float f_SIposArea = 0.0;
  float f_SInegAreaNorm = 0.0;
  float f_SInegArea = 0.0;
  char  pch_curveName[1024];
//  char* pch_function    = "MRIS_surfaceIntegrals_report";
  char  tmp[1024];

  float f_totalSurfaceArea  = apmris->total_area;
  int   totalVertices   = apmris->nvertices;

//  DebugEnterFunction (( pch_function ));
  sprintf(pch_curveName, "%s", Gppch[aesot]);

  surfaceIntegrals_compute(apmris, e_natural,
                           &f_SInatural,
                           &f_SInaturalMean,    &SInaturalVertices,
                           &f_SInaturalAreaNorm,  &f_SInaturalArea);

  surfaceIntegrals_compute(apmris, e_rectified,
                           &f_SIabs,
                           &f_SIabsMean,    &SIabsVertices,
                           &f_SIabsAreaNorm,    &f_SIabsArea);

  surfaceIntegrals_compute(apmris, e_pos,
                           &f_SIpos,
                           &f_SIposMean,    &SIposVertices,
                           &f_SIposAreaNorm,    &f_SIposArea);

  surfaceIntegrals_compute(apmris, e_neg,
                           &f_SIneg,
                           &f_SInegMean,    &SInegVertices,
                           &f_SInegAreaNorm,    &f_SInegArea);

  strcpy(apch_report, "");

  if(aesot == e_Gaussian) {
    Gf_intrinsicCurvaturePos = f_SIpos / 4 / M_PI;
    Gf_intrinsicCurvatureNeg = f_SIneg / 4 / M_PI;
    Gf_intrinsicCurvatureNat = f_SInatural / 4 / M_PI;
  }
  if(aesot == e_FI) {
    Gf_foldingIndex    = f_SInatural / 4 / M_PI;
  }

  if(Gb_filter) {
    // Recompute metric properties since filtering might change
    // vertices that satisfy filter...
    MRIScomputeMetricProperties(apmris);
  }

  if(Gb_filter || label_name) {
    strcpy(pch_processedArea, "ROI Surface");
  } else {
    strcpy(pch_processedArea, "Whole Surface");
  }

  int req = snprintf(pch_misc, STRBUF, " Mean Vertex Separation (%s):", pch_processedArea);
  if( req >= STRBUF ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Curvature Calculation Type:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%12s\n", Gpch_calc);
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s", pch_curveName, pch_misc);
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f +- %2.5f mm\n",
          apmris->avg_vertex_dist, apmris->std_vertex_dist);
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s", pch_curveName, " Total Surface Area:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f mm^2\n", apmris->total_area);
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s", pch_curveName, " Total Number of Vertices:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5d\n", apmris->nvertices);
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Average Vertex Area (Whole Surface):");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f mm^2\n", apmris->avg_vertex_area);
  strcat(apch_report,tmp);

  if(Gb_filter || label_name) {
    int req = snprintf(tmp, 1024, "%10s%-40s", pch_curveName, " ROI Surface Area:");
    if( req >= 1024 ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    strcat(apch_report,tmp);

    sprintf(tmp, "%12.5f mm^2\n", f_SInaturalArea);
    strcat(apch_report,tmp);

    req = snprintf(tmp, 1024, "%10s%-40s", pch_curveName, " ROI Number of Vertices:");
    if( req >= 1024 ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    strcat(apch_report,tmp);

    sprintf(tmp, "%12.5d\n", SInaturalVertices);
    strcat(apch_report,tmp);

    req = snprintf(tmp, 1024, "%10s%-40s",
            pch_curveName, " ROI Surface Area Percentage:");
    if( req >= STRBUF ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    strcat(apch_report,tmp);

    sprintf(tmp, "%12.2f%s\n",
            100 * (float)f_SInaturalArea / apmris->total_area, "%");
    strcat(apch_report,tmp);

    req = snprintf(tmp, 1024, "%10s%-40s",
            pch_curveName, " ROI Surface Vertex Percentage:");
    if( req >= 1024 ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    strcat(apch_report,tmp);

    sprintf(tmp, "%12.2f%s\n",
            100 * (float)SInaturalVertices / apmris->nvertices, "%");
    strcat(apch_report,tmp);


    req = snprintf(tmp, 1024, "%10s%-40s",
            pch_curveName, " Average Vertex Area (ROI Surface):");
    if( req >= 1024 ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    strcat(apch_report,tmp);

    sprintf(tmp, "%12.5f mm^2\n", f_SInaturalArea / SInaturalVertices);
    strcat(apch_report,tmp);

    if(Gb_regionalPercentages) {
      totalVertices     = G_regionalTotalVertexCount;
      f_totalSurfaceArea  = Gf_regionalTotalSurfaceArea;
    }
  }

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Natural Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f\n",
          f_SInatural);
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Rectified Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f\n",
          f_SIabs);
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Positive Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f\n",
          f_SIpos);
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Negative Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp,"%12.5f\n",
          f_SIneg);
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Mean Natural Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f across %d (%05.2f%s) vertices\n",
          f_SInaturalMean, SInaturalVertices,
          100 * (float)SInaturalVertices / totalVertices,
          "%");
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Mean Rectified Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f across %d (%05.2f%s) vertices\n",
          f_SIabsMean, SIabsVertices,
          100 * (float)SIabsVertices / totalVertices, "%");
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Mean Positive Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f across %d (%05.2f%s) vertices\n",
          f_SIposMean, SIposVertices,
          100 * (float)SIposVertices / totalVertices, "%");
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " Mean Negative Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f across %d (%05.2f%s) vertices\n",
          f_SInegMean, SInegVertices,
          100 * (float)SInegVertices / totalVertices, "%");
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " AreaNorm Natural Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f across %f (%05.2f%s) mm^2\n",
          f_SInaturalAreaNorm, f_SInaturalArea,
          100 * f_SInaturalArea / f_totalSurfaceArea, "%");
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " AreaNorm Rectified Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f across %f (%05.2f%s) mm^2\n",
          f_SIabsAreaNorm, f_SIabsArea,
          100 * f_SIabsArea / f_totalSurfaceArea, "%");
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " AreaNorm Positive Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f across %f (%05.2f%s) mm^2\n",
          f_SIposAreaNorm, f_SIposArea,
          100 * f_SIposArea / f_totalSurfaceArea, "%");
  strcat(apch_report,tmp);

  req = snprintf(tmp, 1024, "%10s%-40s",
          pch_curveName, " AreaNorm Negative Surface Integral:");
  if( req >= 1024 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  strcat(apch_report,tmp);

  sprintf(tmp, "%12.5f across %f (%05.2f%s) mm^2\n",
          f_SInegAreaNorm, f_SInegArea,
          100 * f_SInegArea / f_totalSurfaceArea, "%");
  strcat(apch_report,tmp);

  return 1;
}

int
MRIS_curvatureStats_analyze(
  MRIS*     apmris,
  e_secondOrderType aesot
)
{
  //
  // PRECONDITIONS
  //  o <aprmis> is valid
  //  o <aesot> denotes the enumerated surface type to analyze.
  //
  // POSTCONDITIONS
  //  o For the specific <aesot>, the following are analyzed:
  //    - min/max values and vertices (*)
  //    - Surface intergrals are performed and reported
  //    - histograms are prepared (*)
  //    (*) if the appropriate flag has been set by the main
  //    user on the command line.
  //
  // SIDE EFFECTS
  //  o Modifies the current global Gf_mean and Gf_sigma
  //

  char  pch_text[65536];
  char  pch_minMaxReport[65536];
  char  pch_surfaceIntegralReport[65536];
//  char* pch_function  = "MRIS_curvatureStats_analyze";
  const char* pch_unitsmm = "mm^-1";
  const char* pch_unitsmm2  = "mm^-2";
  char  pch_units[256];

//  DebugEnterFunction (( pch_function ));

  strcpy(pch_text, "");
  strcpy(pch_minMaxReport, "");
  strcpy(pch_surfaceIntegralReport, "");

  // Perform the analyses and prepare reports
  MRIS_curvatures_prepare(apmris, aesot);
  if(Gb_minMaxShow) {
    MRIS_minMaxCurve_report(apmris, aesot, pch_minMaxReport);
  }
  MRIS_surfaceIntegrals_report(apmris, aesot, pch_surfaceIntegralReport);

  strcpy(pch_units, "");
  switch(aesot) {
  case e_K1:
  case e_K2:
  case e_Mean:
    strcpy(pch_units, pch_unitsmm);
    break;
  default:
    strcpy(pch_units, pch_unitsmm2);
    break;
  }

  // Now, dump the reports to stdout
  //  First the mean/sigma results
  Gf_mean = MRIScomputeAverageCurvature(apmris, &Gf_sigma);
  if(aesot == e_Raw) {
    int req = snprintf(pch_text, STRBUF,
            "\n%s <mean> +- <std> (using '%s.%s'):",
            Gppch[aesot], hemi, curv_fname);
    if( req >= STRBUF ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    strcpy(pch_units, "");
  } else {
    int req = snprintf(pch_text, STRBUF,
              "\n%s <mean> +- <std> (using '%s.%s'):",
              Gppch[aesot], hemi, surf_name);
    if( req >= STRBUF ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }
  fprintf(GpSTDOUT, "%-50s", pch_text);
  fprintf(GpSTDOUT, " %12.5f +- %2.4f %s\n", Gf_mean, Gf_sigma, pch_units);

  // Now the min/max report
  if(Gb_minMaxShow) {
    fprintf(GpSTDOUT, "%s", pch_minMaxReport);
  }

  // The surface integral report
  fprintf(GpSTDOUT, "%s", pch_surfaceIntegralReport);

  // and any histograms...
  if(Gb_histogram) {
    histogram_wrapper(apmris, aesot);
  }

  return 1;
}

int
MRISvertexCurvature_set(
  MRI_SURFACE*    apmris,
  int       aindex,
  float     af_val
)
{
  //
  // PRECONDITIONS
  //  o <af_val> is typically zero.
  //  o MRIScomputeSecondFundamentalForm() should have been called
  //   if mean/gaussian values are important.
  //
  // POSTCONDITIONS
  //  o The curvature of the vertex at aindex is
  //    set to <af_val>. The Gaussian and Mean are also set
  //    to <af_val>.
  //  o apmris max_curv and min_curv values are recomputed across
  //   the "raw", gaussian, and mean curvatures.
  //

  VERTEX*   pvertex;
  int     vno;
  float   f_maxCurv = 0.;
  float   f_minCurv = 0.;
  float   f_maxCurvK = 0.;
  float   f_minCurvK = 0.;
  float   f_maxCurvH = 0.;
  float   f_minCurvH = 0.;

  if (aindex > apmris->nvertices) {
    ErrorExit(ERROR_SIZE, "%s: target vertex is out of range.", Progname);
  }

  apmris->vertices[aindex].curv = af_val;
  apmris->vertices[aindex].K  = af_val;
  apmris->vertices[aindex].H  = af_val;

  f_maxCurv  = f_minCurv  = apmris->vertices[0].curv;
  f_maxCurvK = f_minCurvK = apmris->vertices[0].K;
  f_maxCurvH = f_minCurvH = apmris->vertices[0].H;
  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    pvertex = &apmris->vertices[vno] ;
    if (pvertex->ripflag) {
      continue;
    }
    if (pvertex->curv > f_maxCurv) {
      f_maxCurv  = pvertex->curv;
    }
    if (pvertex->curv < f_minCurv) {
      f_minCurv  = pvertex->curv;
    }
    if (pvertex->K    > f_maxCurvK) {
      f_maxCurvK = pvertex->K;
    }
    if (pvertex->K    < f_minCurvK) {
      f_minCurvK = pvertex->K;
    }
    if (pvertex->H    > f_maxCurvH) {
      f_maxCurvH = pvertex->H;
    }
    if (pvertex->H    < f_minCurvH) {
      f_minCurvH = pvertex->H;
    }
  }
  apmris->max_curv = f_maxCurv;
  apmris->min_curv = f_minCurv;
  apmris->Kmax = f_maxCurvK;
  apmris->Kmin = f_minCurvK;
  apmris->Hmax = f_maxCurvH;
  apmris->Hmin = f_minCurvH;

  return(NO_ERROR);
}

int
MRISminMaxCurvaturesSearchSOT(
  MRI_SURFACE*    apmris,
  int*      ap_vertexMin,
  int*      ap_vertexMax,
  float*    apf_min,
  float*    apf_max,
  e_secondOrderType   aesot
)
{
  //
  // PRECONDITIONS
  //  o apmris should have its curvature fields applicably set (i.e. Gaussian
  //   and Mean usage has already been called).
  //
  // POSTCONDITIONS
  // o Return the min and max curvatures in apf_min and apf_max respectively,
  //   and similarly the indices in ap_vertex{Min/Max}.
  // o This is an explicit search.
  //

  VERTEX* pvertex;
  int  vno;

  float f_min = 1e6;
  float f_max = -1e6;
  float f_curv = 0.;

  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    pvertex = &apmris->vertices[vno] ;
    switch (aesot) {
    case e_Normal:
    case e_Scaled:
    case e_ScaledTrans:
    case e_C:
    case e_S:
    case e_SI:
    case e_FI:
    case e_BE:
    case e_Raw:
      f_curv  = pvertex->curv;
      break;
    case e_Gaussian:
      f_curv = pvertex->K;
      break;
    case e_Mean:
      f_curv = pvertex->H;
      break;
    case e_K1:
      f_curv = pvertex->k1;
      break;
    case e_K2:
      f_curv = pvertex->k2;
      break;
    }
    if (!vno) {
      f_min = f_max = f_curv;
      *ap_vertexMin = *ap_vertexMax = 0;
    }
    if (pvertex->ripflag) {
      continue;
    }
    if (f_curv > f_max) {
      f_max   = f_curv;
      *ap_vertexMax = vno;
    }
    if (f_curv < f_min) {
      f_min = f_curv;
      *ap_vertexMin = vno;
    }
  }
  *apf_max = f_max;
  *apf_min = f_min;
  return(NO_ERROR);
}

int
MRISminMaxCurvatureIndicesLookup(
  MRI_SURFACE*  apmris,
  int*   ap_vertexMin,
  int*   ap_vertexMax
)
{
  //
  // PRECONDITIONS
  //  o apmris should already have its max_curv and min_curv fields
  //   defined.
  //  o for second order min/max, make sure that apmris has its Kmin/Kmax
  //   and Hmin/Hmax fields defined.
  //
  // POSTCONDITIONS
  // o Return the actual vertex index number corresponding to the
  //   minimum and maximum curvature.
  //

  VERTEX* pvertex;
  int  vno;
  float f_min = apmris->min_curv;
  float f_max = apmris->max_curv;

  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    pvertex = &apmris->vertices[vno] ;
    if (pvertex->ripflag) {
      continue;
    }
    if (pvertex->curv == f_min) {
      *ap_vertexMin = vno;
    }
    if (pvertex->curv == f_max) {
      *ap_vertexMax = vno;
    }
  }
  return(NO_ERROR);
}



int
MRISvertexAreaPostProcess(
  MRI_SURFACE*    pmris
)
{
  //
  // PRECONDITIONS
  //  o The 'curv' member at each vertex contains the curvature
  //    function value of interest.
  //
  // POSTCONDITIONS
  //  o If the Gb_vertexAreaNormalize{Frac} or Gb_vertexAreaWeigh{Frac}
  //    flags are set, process the 'curv' member accordingly.
  //  o Reset the min_curv and max_curv values if necessary.
  //

  int     vno ;
  VERTEX* pvertex ;
  float         f_surfaceArea   = pmris->total_area;

  float   f_min =  pmris->vertices[0].curv;
  float   f_max = f_min;

  for (vno = 0 ; vno < pmris->nvertices ; vno++) {
    pvertex = &pmris->vertices[vno] ;
    if (pvertex->ripflag) {
      continue ;
    }
    if (Gb_vertexAreaNormalize) {
      pvertex->curv /= pvertex->area;
    }
    if (Gb_vertexAreaWeigh) {
      pvertex->curv *= pvertex->area;
    }
    if (Gb_vertexAreaNormalizeFrac) {
      pvertex->curv /= (pvertex->area/f_surfaceArea);
    }
    if (Gb_vertexAreaWeighFrac) {
      pvertex->curv *= (pvertex->area/f_surfaceArea);
    }
    if (pvertex->curv < f_min) {
      f_min = pvertex->curv;
    }
    if (pvertex->curv > f_max) {
      f_max = pvertex->curv;
    }
  }

  pmris->min_curv = f_min ;
  pmris->max_curv = f_max;
  return(NO_ERROR) ;
}

void *
xmalloc (size_t size)
{
  void *value = malloc (size);
  if (value == 0) {
    ErrorExit(10, "%s: virtual memory exhausted.", Progname);
  }
  return value;
}

short
LABEL_save(
  MRI_SURFACE*  apmris,
  int   a_labelSize,
  int*    ap_labelMark,
  char*   apch_filename
)
{
  //
  // DESCRIPTION
  //  This function saves a FreeSurfer format 'LABEL' file
  //  to the file <apch_filename>. The integer array <ap_labelMark>
  //  has length <aprmis->nvertices> and denotes the vertices to
  //  label by a non-zero element value.
  //
  // PRECONDITIONS
  //  o ap_labelMark denotes the vertices to label -- marked
  //    vertices are indicated by a non-zero element value and
  //    index counter corresponds to vertex number.
  //
  // POSTCONDITIONS
  //    o The file <apch_filename> contains the marked vertices, saved
  //      as a FreeSurfer label file.
  //
  // RETURN
  //  o TRUE, i.e. 1.
  //

  LABEL*  pLBL      = NULL;
  int   vno     = 0;

  pLBL = LabelAlloc(a_labelSize, "", apch_filename);
  for(vno=0; vno<apmris->nvertices; vno++) {
    if(ap_labelMark[vno]) {
      pLBL->lv[pLBL->n_points].vno  = vno;
      pLBL->lv[pLBL->n_points].x    = apmris->vertices[vno].x;
      pLBL->lv[pLBL->n_points].y    = apmris->vertices[vno].y;
      pLBL->lv[pLBL->n_points].z    = apmris->vertices[vno].z;
      pLBL->n_points++;
    }
  }
  LabelWrite(pLBL, apch_filename);
  LabelFree(&pLBL);
  return 1;
}

int
MRIS_surfaceRipFlags_filterVertex(
  MRI_SURFACE*  apmris,
  int   avno
)
{
  //
  // DESCRIPTION
  //  Depending on a user-specified filter pattern, this
  //  function examines the passed vertex, and sets its
  //  ripflag accordingly.
  //
  //  This "ripping" is performed according to the pattern
  //  of 'filter' variables that might have been set at the
  //  command line.
  //
  // PRECONDITIONS
  //  o A valid surface.
  //  o The 'curv' value at each vertex MUST be set for current
  //    processing regime.
  //
  // POSTCONDITIONS
  //  o If passed vertex index <vno> does not satisfy
  //    filter contraints, it is "ripped"!
  //
  // RETURN
  //  o If vertex is within a label region (or whole hemisphere) it is
  //    not ripped, and thus return 0; else return 1.
  //

  VERTEX* pv ;
  short b_canProcess  = 1;

  pv = &apmris->vertices[avno] ;
  // Check if vertex might have been ripped by a label load...
  if(pv->ripflag) {
    return pv->ripflag;
  }
  // Process filter patterns...
  if(Gb_lowPassFilter) {
    if(fabs(pv->curv)<fabs(Gf_lowPassFilter)) {
      b_canProcess  &= 1;
    } else {
      b_canProcess  &= 0;
    }
  }
  if(Gb_highPassFilter) {
    if(fabs(pv->curv)>=fabs(Gf_highPassFilter)) {
      b_canProcess  &= 1;
    } else {
      b_canProcess  &= 0;
    }
  }
  if(Gb_lowPassFilterGaussian) {
    if(fabs(pv->K)<fabs(Gf_lowPassFilterGaussian)) {
      b_canProcess  &= 1;
    } else {
      b_canProcess  &= 0;
    }
  }
  if(Gb_highPassFilterGaussian) {
    if(fabs(pv->K)>=fabs(Gf_highPassFilterGaussian)) {
      b_canProcess  &= 1;
    } else {
      b_canProcess  &= 0;
    }
  }
  if(b_canProcess) {
    pv->ripflag = 0;
  } else {
    pv->ripflag = 1;
  }
  return 0;
}

int
MRIS_surfaceRipFlags_filter(
  MRI_SURFACE*  apmris,
  float*    apf_notRippedArea
)
{
  //
  // DESCRIPTION
  //  This function "rips" a surface (i.e. sets the 'ripgflag'
  //  field at each vertex) to indicate whether a vertex is a
  //  candidate for further processing. A "ripped" vertex
  //  should not be processed.
  //
  //  This "ripping" is performed according to the pattern
  //  of 'filter' variables that might have been set at the
  //  command line.
  //
  // PRECONDITIONS
  //  o A valid surface.
  //
  // POSTCONDITIONS
  //  o All vertices that satisfy the command line specified
  //    filter constraints have their ripflags set to 0, those
  //    that do not satisfy filter constraints have their
  //    ripflags set to 1.
  //  o Any existing rip pattern is overwritten!
  //
  // RETURN
  //  o Returns the number of vertices that can be processed,
  //    i.e. number of vertices that were *not* ripped.
  //  o Returns the sum area of valid vertices in <apf_notRippedArea>.
  //

  int         vno ;
  int   ripped    = 0;
  int   notRipped = 0;
  float f_rippedArea  = 0;

  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    ripped = MRIS_surfaceRipFlags_filterVertex(apmris, vno);
    if(!ripped) {
      notRipped++;
      f_rippedArea += apmris->vertices[vno].area;
    }
  }
  *apf_notRippedArea  = f_rippedArea;
  return notRipped;
}

short
MRIS_surfaceIntegral_compute(
  MRI_SURFACE*  pmris,
  int*    p_verticesCounted,
  float*    pf_areaCounted,
  short   (*fcond)(VERTEX*  pv),
  float   (*fv) (VERTEX*  pv),
  float*    pf_surfaceIntegral,
  float*    pf_meanIntegral,
  float*    pf_areaNormalizedIntegral
)
{
  //
  // DESCRIPTION
  //  This function computes a surface integral across the
  //  vertex curvature. The resultant sum is normalized to the
  //  number of the vertices that were counted and also normalized
  //  to the total area of the counted vertices.
  //
  //  Before processing a particular vertex curvature, the vertex
  //  is first checked by the (*fcond) function for eligibility.
  //
  //  The actual curvature that is integrated can also be shaped
  //  by the (*fv) function.
  //
  //  Global low pass filtering (on actual curvature values and
  //  Gaussian test) is always implemented.
  //
  // PRECONDITIONS
  //  o The vertex curvature 'curv' contains the particular
  //    function curvature of interest.
  //  o The (*fcond) is a conditional function that checks the
  //    curvature value for processing. Such tests can for example
  //    only trigger on negative curvatures, or only positive, etc.
  //  o The (*fcurv) is an additional function of the vertex curvature
  //    itself - often this will return the absolute value of the
  //    curvature.
  //
  // POSTCONDITIONS
  //  o The integral of the curvature value over the surface is returned
  //    according to various criteria.
  //  o <pf_surfaceIntegral> is the actual integral for the given curvatures
  //    in the vertex 'curv'.
  //  o The <pf_meanIntegral> and the <pf_areaNormalizedIntegral> are the
  //    surface intergrals normalized to either the number of vertices or
  //    the integral area.
  //  o If the integral could not be computed, -1 is returned.
  //  o If the global <Gb_lowPassFilter> is set, only count vertices
  //    where abs 'curv' is less than <Gb_lowPassFilter>.
  //  o If the global <Gb_lowPassFilterGaussian> is set, only count
  //    vertices where the abs Gaussian at that vertex is less than
  //    <Gb_lowPassFilterGaussian>.
  //

  VERTEX* pv ;
  int         vno ;
  double      f_total, f_n, f_totalArea ;
  short   b_ret     = 0;
  int*    p_labelMark   = NULL;
  static short  b_labelFileCreated  = 0;

  *p_verticesCounted  = 0;
  p_labelMark   = (int*) xmalloc(pmris->nvertices * sizeof(int));
  if(Gb_filterLabel && !b_labelFileCreated)
    for(vno=0; vno<pmris->nvertices; vno++) {
      p_labelMark[vno] = 0;
    }
  for (f_n = f_total =f_totalArea = 0.0, vno = 0 ; vno < pmris->nvertices ; vno++) {
    pv = &pmris->vertices[vno] ;
    // This call to MRIS_surfaceRipFlags_filterVertex is needed
    // since Gaussian curvatures are now available.
//    if(Gb_filter)   MRIS_surfaceRipFlags_filterVertex(pmris, vno);
    if(pv->ripflag) {
      continue ;
    }
    if( (*fcond)(pv)) {
      f_total       += ((*fv)(pv) * pv->area) ;
      f_n       += 1.0 ;
      (*pf_areaCounted)   += pv->area;
      (*p_verticesCounted)++;
      if(Gb_filterLabel) {
        p_labelMark[vno] = 1;
      }
    }
  }
  if(f_n > 1) {
    *pf_meanIntegral    = f_total / f_n;
    *pf_areaNormalizedIntegral  = f_total / (*pf_areaCounted);
    b_ret       = 1;
  }
  if(Gb_postScale) {
    *pf_meanIntegral    *= Gf_postScale;
    *pf_areaNormalizedIntegral  *= Gf_postScale;
  }
  *pf_surfaceIntegral   = f_total;
  if(Gb_filterLabel && !b_labelFileCreated) {
    LABEL_save(pmris, *p_verticesCounted, p_labelMark, Gpch_filterLabel);
    cprintd("Label size", *p_verticesCounted);
    free(p_labelMark);
    b_labelFileCreated  = 1;
  }
  return(b_ret);
}

int
MRISzeroCurvature(
  MRI_SURFACE* apmris
)
{
  //
  // POSTCONDITIONS
  // o Each curvature (as well as Gaussian and Mean )value in
  //   apmris is simply set to zero.
  //

  VERTEX* pvertex;
  int  vno;

  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    pvertex = &apmris->vertices[vno] ;
    if (pvertex->ripflag) {
      continue;
    }
    pvertex->curv  = 0.;
    pvertex->K = 0.;
    pvertex->H = 0.;
  }
  return(NO_ERROR);
}

void
OFSP_create(
  char*  apch_prefix,
  char*  apch_suffix,
  e_OFSP  ae_OFSP
)
{
  //
  // PRECONDITIONS
  // o pch_prefix / pch_suffix must be large enough to contain
  //   required text info.
  //
  // POSTCONDITIONS
  //  o a full pch_prefix is composed of:
  //  <output_fname>.<hemi>.<surface>.<pch_suffix>
  //  o a partial pch_prefix is compsed of:
  //  <hemi>.<surface>.<pch_suffix>
  //

  switch (ae_OFSP) {
  case e_None:
    if(output_fname != NULL) {
      sprintf(apch_prefix, "%s", output_fname);
    }
    break;
  case e_Partial:
    sprintf(apch_prefix, "%s.%s.%s",
            hemi,
            surf_name,
            apch_suffix
           );
    break;
  case e_Full:
    sprintf(apch_prefix, "%s.%s.%s.%s",
            output_fname,
            hemi,
            surf_name,
            apch_suffix
           );
    break;
  }
}

void
histogram_wrapper(
  MRIS*   amris,
  e_secondOrderType aesot
)
{
  //
  // PRECONDITIONS
  // o amris_curvature must be valid and populated with curvature values
  //   of correct type (i.e. normalised, Gaussian/Mean, scaled)
  // o apf_histogram memory must be handled (i.e. created/freed) by the
  //   caller.
  //
  // POSTCONDITIONS
  //  o max/min_curvatures and binSize are checked for validity.
  // o histogram memory is created / freed.
  //  o histogram_create() is called.
  //

  int     i;
  float*  pf_histogram;
  float   f_maxCurv   = amris->max_curv;
  float   f_minCurv   = amris->min_curv;
  double    f_binSize   = 0.;
  int     b_error   = 0;
  short   b_writeToFile = 0;
  FILE*   pFILE   = NULL;

  if (Gb_binSizeOverride) {
    f_binSize = Gf_binSize;
  }

  if (Gb_histStartOverride) {
    f_minCurv = Gf_histStart;
  }

  if (f_minCurv < amris->min_curv) {
    f_minCurv = amris->min_curv;
  }

  if (Gb_histEndOverride) {
    f_maxCurv = Gf_histEnd;
  }

  if (f_maxCurv > amris->max_curv) {
    f_maxCurv = amris->max_curv;
  }

  f_binSize  = ((double)f_maxCurv - (double)f_minCurv) / (double)G_bins;
  if (f_maxCurv <= f_minCurv)
    ErrorExit(ERROR_SIZE, "%s: f_maxCurv (%f) < f_minCurv (%f)",
              Progname, f_maxCurv, f_minCurv);

  b_error =   (long)(((double)f_minCurv+(double)G_bins*(double)f_binSize)*1e10) >
              (long)(((double)f_maxCurv)*1e10);

  if (b_error)
    ErrorExit(ERROR_SIZE, "%s: Invalid <binSize> and <bins> combination",
              Progname);

  pf_histogram = (float *)calloc(G_bins, sizeof(float));
  histogram_create( amris,
                    f_minCurv,
                    f_binSize,
                    G_bins,
                    pf_histogram,
                    aesot);

  printf("%*s%*s%*s\n", G_leftCols, "\tbin start",
         G_leftCols*2, "\tbin end",
         G_leftCols*2, "\tcount");
  for (i=0; i<G_bins; i++) {
    printf("%*.4f%*.4f%*.4f ",
           G_leftCols*2, (i*f_binSize)+f_minCurv,
           G_leftCols*2, ((i+1)*f_binSize)+f_minCurv,
           G_leftCols*2, pf_histogram[i]);
    printf("\n");
    switch (aesot) {
    case e_Raw:
      if (GpFILE_rawHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_rawHist;
      }
      break;
    case e_Gaussian:
      if (GpFILE_KHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_KHist;
      }
      break;
    case e_Mean:
      if (GpFILE_HHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_HHist;
      }
      break;
    case e_K1:
      if (GpFILE_K1Hist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_K1Hist;
      }
      break;
    case e_K2:
      if (GpFILE_K2Hist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_K2Hist;
      }
      break;
    case e_S:
      if (GpFILE_SHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_SHist;
      }
      break;
    case e_C:
      if (GpFILE_CHist !=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_CHist;
      }
      break;
    case e_BE:
      if (GpFILE_BEHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_BEHist;
      }
      break;
    case e_SI:
      if (GpFILE_SIHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_SIHist;
      }
      break;
    case e_FI:
      if (GpFILE_FIHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_FIHist;
      }
      break;
    case e_Normal:
      if (GpFILE_normHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_normHist;
      }
      break;
    case e_ScaledTrans:
    case e_Scaled:
      if (GpFILE_scaledHist!=NULL) {
        b_writeToFile = 1;
        pFILE = GpFILE_scaledHist;
      }
      break;
    }
    if(b_writeToFile)
      fprintf(pFILE, "%f\t%f\t%f\n",
              (i*f_binSize)+f_minCurv,
              ((i+1)*f_binSize)+f_minCurv,
              pf_histogram[i]);
  }
  free(pf_histogram);
}

// Usable AlmostEqual function
// Bruce Dawson, see
//  http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
int AlmostEqual2sComplement(float A, float B, int maxUlps)
{
  int aInt;
  int bInt;
  int intDiff;
  void* ptrA = (void*)&A;
  void* ptrB = (void*)&B;

  // Make sure maxUlps is non-negative and small enough that the
  // default NAN won't compare as equal to anything.
  assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
  aInt = *(int*)ptrA;
  // Make aInt lexicographically ordered as a twos-complement int
  if (aInt < 0) {
    aInt = 0x80000000 - aInt;
  }
  // Make bInt lexicographically ordered as a twos-complement int
  bInt = *(int*)ptrB;
  if (bInt < 0) {
    bInt = 0x80000000 - bInt;
  }
  intDiff = abs(aInt - bInt);
  if (intDiff <= maxUlps) {
    return true;
  }
  return false;
}

void
histogram_create(
  MRIS*   amris_curvature,
  float   af_minCurv,
  double   af_binSize,
  int   abins,
  float*   apf_histogram,
  e_secondOrderType aesot
)
{
  //
  // PRECONDITIONS
  // o amris_curvature must be valid and populated with curvature values
  //   of correct type (i.e. normalised, Gaussian/Mean, scaled)
  // o apf_histogram memory must be handled (i.e. created/freed) by the
  //   caller.
  //  o af_minCurv < af_maxCurv.
  //
  // POSTCONDITIONS
  // o curvature values in amris_curvature will be sorted into abins.
  //
  // HISTORY
  // 12 September 2005
  // o Initial design and coding.
  //

  //char* pch_function = "histogram_create()";

  int     vno     = 0;
  float*  pf_curvature  = NULL;
  int     i     = 0;
  int     j     = 0;
  int     start     = 0;
  int     count     = 0;
  int     totalCount  = 0;
  int     nvertices   = 0;

  int     b_onLeftBound = 0;  // Don't trigger higher order float
  int     b_onRightBound = 0;   // comparison on left/right bounds

  double  l_curvature;    // These three variables
  double  l_leftBound;    // are scaled up and truncated
  double  l_rightBound;     // to minimise rounding errors

  char pch_sot[64];
  strcpy(pch_sot, Gppch[aesot]);

  nvertices  = amris_curvature->nvertices;
  fprintf(GpSTDOUT, "\n%*s%s = %f\n",
          G_leftCols, pch_sot, " bin size", af_binSize);
  fprintf(GpSTDOUT, "%*s%s = %d\n",
          G_leftCols, pch_sot, " surface vertices", nvertices);
  pf_curvature = (float *)calloc(nvertices, sizeof(float));

  for (vno=0; vno < nvertices; vno++) {
    pf_curvature[vno] = amris_curvature->vertices[vno].curv;
  }

  qsort(pf_curvature, nvertices, sizeof(float), comp);

  for (i=0; i<abins; i++) {
    count = 0;
    for (j=start; j<nvertices; j++) {
      l_curvature  = (double)(pf_curvature[j]);
      l_leftBound  = (double)((i*af_binSize+af_minCurv));
      l_rightBound = (double)((((i+1)*af_binSize)+af_minCurv));
      if (Gb_maxUlps) {
        b_onRightBound =
          AlmostEqual2sComplement(l_curvature,l_rightBound,G_maxUlps);
        b_onLeftBound =
          AlmostEqual2sComplement(l_curvature,l_leftBound,G_maxUlps);
      }
      if ( (l_curvature >= l_leftBound && l_curvature <= l_rightBound)
           ||
           b_onLeftBound || b_onRightBound ) {
        count++;
        totalCount++;
      }
      if (l_curvature > l_rightBound && !b_onRightBound) {
        start = j;
        break;
      }
    }
    apf_histogram[i] = count;
    if (Gb_histogramPercent) {
      apf_histogram[i] /= nvertices;
    }
    if (totalCount == nvertices) {
      break;
    }
  }
  fprintf(GpSTDOUT, "%*s%s = %d\n",
          G_leftCols, pch_sot, " sorted vertices", totalCount);
  fprintf(GpSTDOUT, "%*s%s = %f\n",
          G_leftCols, pch_sot, " ratio", (float)totalCount / (float)nvertices);
  free(pf_curvature);
}

static int
get_option(int argc, char *argv[])
{
  int    nargs   = 0;
  char*  option;
  char*  pch_text;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-lowPassFilter")) {
    Gb_filter     = 1;
    Gb_lowPassFilter    = 1;
    Gf_lowPassFilter    = atof(argv[2]);
    nargs     = 1;
    cprintf("Setting rectified low pass filter to", Gf_lowPassFilter);
  } else if (!stricmp(option, "-regionalPercentages")) {
    Gb_regionalPercentages  = 1;
    cprints("Toggling regional percentages", "ok");
  } else if (!stricmp(option, "-shapeIndex")) {
    Gb_shapeIndex   = 1;
    cprints("Toggling shape index map on", "ok");
  } else if (!stricmp(option, "-vertexAreaNormalize")) {
    Gb_vertexAreaNormalize  = 1;
    cprints("Toggling vertex area normalize on", "ok");
  } else if (!stricmp(option, "-vertexAreaWeigh")) {
    Gb_vertexAreaWeigh    = 1;
    cprints("Toggling vertex area weighing on", "ok");
  } else if (!stricmp(option, "-vertexAreaNormalizeFrac")) {
    Gb_vertexAreaNormalizeFrac   = 1;
    cprints("Toggling fractional vertex area normalize on", "ok");
  } else if (!stricmp(option, "-vertexAreaWeighFrac")) {
    Gb_vertexAreaWeighFrac       = 1;
    cprints("Toggling fractional vertex area weighing on", "ok");
  } else if (!stricmp(option, "-lowPassFilterGaussian")) {
    Gb_filter     = 1;
    Gb_lowPassFilterGaussian  = 1;
    Gf_lowPassFilterGaussian  = atof(argv[2]);
    nargs     = 1;
    cprintf("Setting rectified low pass Gaussian filter",
            Gf_lowPassFilterGaussian);
  } else if (!stricmp(option, "-highPassFilterGaussian")) {
    Gb_filter     = 1;
    Gb_highPassFilterGaussian = 1;
    Gf_highPassFilterGaussian = atof(argv[2]);
    nargs     = 1;
    cprintf("Setting rectified high pass Gaussian filter",
            Gf_highPassFilterGaussian);
  } else if (!stricmp(option, "-highPassFilter")) {
    Gb_filter     = 1;
    Gb_highPassFilter   = 1;
    Gf_highPassFilter   = atof(argv[2]);
    nargs     = 1;
    cprintf("Setting rectified high pass filter", Gf_highPassFilter);
  } else if (!stricmp(option, "-postScale")) {
    Gb_postScale    = 1;
    Gf_postScale    = atof(argv[2]);
    nargs     = 1;
    cprintf("Setting post scale factor", Gf_postScale);
  } else if (!stricmp(option, "-filterLabel")) {
    Gb_filterLabel    = 1;
    pch_text      = argv[2];
    strcpy(Gpch_filterLabel, pch_text);
    nargs     = 1;
    cprints("Filter hits saved to ", Gpch_filterLabel);
  } else if (!stricmp(option, "-discrete")) {
    Gb_discreteCurvaturesUse  = 1;
    cprints("Using discrete curvature calculations", "ok");
  } else if (!stricmp(option, "-continuous")) {
    Gb_discreteCurvaturesUse  = 0;
    cprints("Using continuous curvature calculations", "ok");
  } else if (!stricmp(option, "-writeCurvatureFiles")) {
    Gb_writeCurvatureFiles  = 1;
    cprints("Toggling save flag on curvature files", "ok");
  } else if (!stricmp(option, "-signedPrincipals")) {
    Gb_signedPrincipals = 1;
    cprints("Toggling 'signedPrincipals' flag ON", "ok");
  } else switch (toupper(*option)) {
    case 'O':
      output_fname      = argv[2] ;
      Gb_output2File    = 1;
      Gb_writeCurvatureFiles  = 1;
      nargs         = 1 ;
      cprints("Outputting results using filestem", output_fname) ;
      cprints("Toggling save flag on curvature files", "ok");
      break ;
    case 'F':
      pch_text = argv[2];
      strcpy(surf_name, pch_text);
      nargs = 1 ;
      break;
    case 'L':
      label_name = argv[2] ;
      cprints("Using label", label_name) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      cprintd("Averaging curvature times", navgs);
      nargs = 1 ;
      break ;
    case 'C':
      Gb_scale  = 1;
      Gf_scaleFactor  = atof(argv[2]) ;
      cprintf("Setting raw scale factor", Gf_scaleFactor);
      nargs   = 1 ;
      break ;
    case 'D':
      Gb_scaleMin  = 1;
      Gf_scaleMin  = atof(argv[2]);
      cprintf("Setting scale min factor", Gf_scaleMin);
      nargs  = 1;
      break;
    case 'E':
      Gb_scaleMax  = 1;
      Gf_scaleMax  = atof(argv[2]);
      cprintf("Setting scale max factor", Gf_scaleMax);
      nargs  = 1;
      break;
    case 'G':
      Gb_gaussianAndMean = 1;
      break ;
    case 'S':   /* write out stats */
      stat_flag   = 1 ;
      condition_no  = atoi(argv[2]) ;
      nargs   = 1 ;
      cprintd("Setting out summary statistics as condition", condition_no);
      break ;
    case 'H':   /* histogram */
      if (argc == 2) {
        print_usage();
      }
      Gb_histogram = 1;
      Gb_histogramPercent = 0;
      Gb_minMaxShow = 1;
      G_bins   = atoi(argv[2]);
      nargs   = 1 ;
      cprintf("Setting curvature histogram bins", G_bins);
      if (G_bins <=0 ) {
        ErrorExit(ERROR_BADPARM, "%s: Invalid bin number.\n", Progname);
      }
      break ;
    case 'P':   /* percentage histogram */
      Gb_histogram = 1;
      Gb_histogramPercent = 1;
      Gb_minMaxShow = 1;
      G_bins   = atoi(argv[2]);
      nargs   = 1 ;
      cprintd("Creating percentage curvature histogram bins", G_bins);
      break ;
    case 'B':   /* binSize override*/
      Gb_binSizeOverride = 1;
      Gf_binSize   = (double) atof(argv[2]);
      nargs   = 1 ;
      cprintf("Setting curvature histogram binSize", Gf_binSize);
      break;
    case 'I':   /* histogram start */
      Gb_histStartOverride = 1;
      Gf_histStart   = atof(argv[2]);
      nargs    = 1 ;
      cprintf("Setting histogram start point", Gf_histStart);
      break;
    case 'J':   /* histogram end */
      Gb_histEndOverride = 1;
      Gf_histEnd   = atof(argv[2]);
      nargs   = 1 ;
      cprintf("Setting histogram end point", Gf_histEnd);
      break;
    case 'Z':   /* zero a target vertex */
      Gb_zeroVertex = 1;
      G_zeroVertex = atoi(argv[2]);
      nargs   = 1 ;
      cprintd("Setting zero vertex index", G_zeroVertex);
      break;
    case 'Q':   /* set maxUlps */
      Gb_maxUlps  = 1;
      G_maxUlps  = atoi(argv[2]);
      nargs   = 1 ;
      cprintd("Setting maxUlps", G_maxUlps);
      break;
    case 'N':
      normalize_flag  = 1 ;
      cprints("Setting normalisation", "ON");
      break ;
    case 'M':
      Gb_minMaxShow = 1;
      break;
    case 'V':
      print_version() ;
      exit(1);
      break;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case '-':
      break;
    default:
      fprintf(stderr, "Unknown option '%s'. Looking for help? Use '-u' instead.\n", argv[1]) ;
      exit(1) ;
      break ;
    }
  return(nargs) ;
}

void
outputFileNames_create(void)
{
  //
  // POSTCONDITIONS
  // o All necessary (depending on user flags) file names are created.
  //
  OFSP_create(Gpch_log,   Gpch_logS,    e_None);
  OFSP_create(Gpch_rawHist,   Gpch_rawHistS,    e_Full);
  OFSP_create(Gpch_normHist,  Gpch_normHistS,   e_Full);
  OFSP_create(Gpch_normCurv,  Gpch_normCurvS,   e_Partial);
  OFSP_create(Gpch_KHist,   Gpch_KHistS,    e_Full);
  OFSP_create(Gpch_KCurv,   Gpch_KCurvS,    e_Partial);
  OFSP_create(Gpch_HHist,   Gpch_HHistS,    e_Full);
  OFSP_create(Gpch_HCurv,   Gpch_HCurvS,    e_Partial);
  OFSP_create(Gpch_K1Hist,    Gpch_K1HistS,     e_Full);
  OFSP_create(Gpch_K1Curv,    Gpch_K1CurvS,     e_Partial);
  OFSP_create(Gpch_K2Hist,    Gpch_K2HistS,     e_Full);
  OFSP_create(Gpch_K2Curv,    Gpch_K2CurvS,     e_Partial);
  OFSP_create(Gpch_SHist,   Gpch_SHistS,    e_Full);
  OFSP_create(Gpch_SCurv,   Gpch_SCurvS,    e_Partial);
  OFSP_create(Gpch_CHist,   Gpch_CHistS,    e_Full);
  OFSP_create(Gpch_CCurv,   Gpch_CCurvS,    e_Partial);
  OFSP_create(Gpch_BEHist,    Gpch_BEHistS,     e_Full);
  OFSP_create(Gpch_BECurv,    Gpch_BECurvS,     e_Partial);
  OFSP_create(Gpch_SIHist,    Gpch_SIHistS,     e_Full);
  OFSP_create(Gpch_SICurv,    Gpch_SICurvS,     e_Partial);
  OFSP_create(Gpch_FIHist,    Gpch_FIHistS,     e_Full);
  OFSP_create(Gpch_FICurv,    Gpch_FICurvS,     e_Partial);
  OFSP_create(Gpch_scaledHist,  Gpch_scaledHistS, e_Full);
  OFSP_create(Gpch_scaledCurv,  Gpch_scaledCurvS, e_Partial);
}

void
outputFiles_open(void)
{
  //
  // POSTCONDITIONS
  // o All necessary (depending on user flags) files are opened.
  //

  if (Gb_output2File) {
    if ((GpFILE_log=fopen(Gpch_log, "a"))==NULL)
      ErrorExit(ERROR_NOFILE, "%s: Could not open file '%s' for appending.\n",
                Progname, Gpch_log);
    GpSTDOUT  = GpFILE_log;
    if (Gb_histogram && curv_fname) {
      if ((GpFILE_rawHist=fopen(Gpch_rawHist, "w"))==NULL)
        ErrorExit(ERROR_NOFILE, "%s: Could not open file '%s' for writing.\n",
                  Progname, Gpch_rawHist);
    }
    if (normalize_flag) {
      if (Gb_histogram) {
        if ((GpFILE_normHist=fopen(Gpch_normHist, "w"))==NULL)
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_normHist);
      }
    }
    if (Gb_gaussianAndMean) {
      if (Gb_histogram) {
        if ((GpFILE_KHist=fopen(Gpch_KHist, "w"))==NULL) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_KHist);
        }
        if ((GpFILE_HHist=fopen(Gpch_HHist, "w"))==NULL) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_HHist);
        }
        if ((GpFILE_K1Hist=fopen(Gpch_K1Hist, "w"))==NULL) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_K1Hist);
        }
        if ((GpFILE_K2Hist=fopen(Gpch_K2Hist, "w"))==NULL) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_K2Hist);
        }
        if ((GpFILE_SHist=fopen(Gpch_SHist, "w"))==NULL) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_SHist);
        }
        if ((GpFILE_CHist=fopen(Gpch_CHist, "w"))==NULL) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_CHist);
        }
        if ((GpFILE_BEHist=fopen(Gpch_BEHist, "w"))==NULL) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_BEHist);
        }
        if ((GpFILE_SIHist=fopen(Gpch_SIHist, "w"))==NULL && Gb_shapeIndex) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_SIHist);
        }
        if ((GpFILE_FIHist=fopen(Gpch_FIHist, "w"))==NULL) {
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_FIHist);
        }
      }
    }
    if (Gb_scale || (Gb_scaleMin && Gb_scaleMax)) {
      if (Gb_histogram) {
        if ((GpFILE_scaledHist=fopen(Gpch_scaledHist, "w"))==NULL)
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_scaledHist);
      }
    }
  }
}

void
outputFiles_close(void)
{
  //
  // POSTCONDITIONS
  // o Checks for any open files and closes them.
  //

  if (GpFILE_log) {
    fclose(GpFILE_log);
  }
  GpFILE_log     = NULL;
  if (GpFILE_rawHist) {
    fclose(GpFILE_rawHist);
  }
  GpFILE_rawHist   = NULL;
  if (GpFILE_normHist) {
    fclose(GpFILE_normHist);
  }
  GpFILE_normHist  = NULL;
  if (GpFILE_KHist) {
    fclose(GpFILE_KHist);
  }
  GpFILE_KHist   = NULL;
  if (GpFILE_HHist) {
    fclose(GpFILE_HHist);
  }
  GpFILE_HHist   = NULL;
  if (GpFILE_K1Hist) {
    fclose(GpFILE_K1Hist);
  }
  GpFILE_K1Hist  = NULL;
  if (GpFILE_K2Hist) {
    fclose(GpFILE_K2Hist);
  }
  GpFILE_K2Hist  = NULL;
  if (GpFILE_SHist) {
    fclose(GpFILE_SHist);
  }
  GpFILE_SHist   = NULL;
  if (GpFILE_CHist) {
    fclose(GpFILE_CHist);
  }
  GpFILE_CHist   = NULL;
  if (GpFILE_BEHist) {
    fclose(GpFILE_BEHist);
  }
  GpFILE_BEHist  = NULL;
  if (GpFILE_SIHist) {
    fclose(GpFILE_SIHist);
  }
  GpFILE_SIHist  = NULL;
  if (GpFILE_FIHist) {
    fclose(GpFILE_FIHist);
  }
  GpFILE_FIHist  = NULL;
  if (GpFILE_scaledHist) {
    fclose(GpFILE_scaledHist);
  }
  GpFILE_scaledHist  = NULL;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  print_help() ;
}

static void
print_help(void)
{
  char  pch_synopsis[65536];

  sprintf(pch_synopsis, "\n\
 \n\
 \n\
    NAME \n\
 \n\
          mris_curvature_stats \n\
 \n\
    SYNOPSIS \n\
 \n\
          mris_curvature_stats [options] 			\\ \n\
          	<subjectName> <hemi> [<curvFile1> ... <curvFileN] \n\
 \n\
    DESCRIPTION \n\
 \n\
          In its simplest usage, 'mris_curvature_stats' will compute a set \n\
	  of statistics on its input <curvFile>. These statistics are the \n\
	  mean and standard deviation of the particular curvature on the \n\
	  surface, as well as the results from several surface-based \n\
	  integrals. \n\
 \n\
          Additionally, 'mris_curvature_stats' can report the max/min \n\
          curvature values, and compute a simple histogram based on \n\
          all curvature values. \n\
 \n\
          Curvatures can also be normalised and constrained to a given \n\
          range before computation. \n\
 \n\
          Principal curvature (K, H, k1 and k2) calculations on a surface \n\
	  structure can also be performed, as well as several functions \n\
	  derived from k1 and k2. \n\
 \n\
          Finally, all output to the console, as well as any new \n\
          curvatures that result from the above calculations can be \n\
          saved to a series of text and binary-curvature files. \n\
 \n\
	PRINCIPAL CURVATURES AND FUNCTIONS \n\
 \n\
	  Given a surface file, 'mris_curvature_stats' can also compute \n\
	  all the principal curvatures relating to each point on the \n\
	  surface as well as several functions of these principal \n\
	  curvatures. \n\
 \n\
	  To calculate these principal curvatures, use a '-G' flag on the \n\
	  command line. In such an instance, you do not need (nor probably even \n\
	  want) any <curvFile> arguments. The following principal curvatures \n\
	  and derived measures are calculated (and statistics on each are \n\
	  presented) for each vertex on the surface: \n\
 \n\
	  		k1 	maximum curvature \n\
	  		k2 	minimum curvature	 \n\
	  		K 	Gaussian 	= k1*k2 \n\
	  		H	Mean 		= 0.5*(k1+k2) \n\
          		C	Curvedness 	= sqrt(0.5*(k1*k1+k2*k2)) \n\
			S	Sharpness 	= (k1 - k2)^2 \n\
	  		BE	Bending Energy 	= k1*k1 + k2*k2 \n\
			SI	Shape Index	= atan((k1+k2)/(k2-k1)) \n\
			FI	Folding Index	= |k1|*(|k1| - |k2|) \n\
	 \n\
	  Note that the 'BE' is not the same as the Willmore Bending Energy. \n\
	  The Willmore Bending Energy is the same as 'S' (sharpness). The \n\
	  BE in this case only considers the square of the principal \n\
	  curvatures. \n\
 \n\
	  Also note that the SI is not calculated by default due to some issues \n\
	  with atan calculations. Use the '--shapeIndex' flag on the command \n\
	  line to force SI. \n\
 \n\
	  In addition, if a '--writeCurvatureFiles' is specified, each of the \n\
	  above are written to a FreeSurfer format curvature file. \n\
 \n\
	  Note that there are two methods to calculate principal curvatures. \n\
	  The first method is based on calculating the Second Fundamental \n\
	  Form, and is the method used by 'mris_anatomical_stats'. This \n\
	  approach can be selected by the '--continuous' command line \n\
	  argument. This method, however, suffers occassionally from \n\
	  accuracy problems and can have very large curvature outliers. \n\
 \n\
	  A slower, more accurate method using discrete methods is available \n\
       	  and can be selected with the '--discrete' command line argument. \n\
	  This is in fact the default mode for 'mris_curvature_stats'. \n\
 \n\
	NOTES ON CURVATURE / FOLDING INDICES \n\
	 \n\
 	  If a '-G' has been specified, 'mris_curvature_stats' will also \n\
 	  output some folding/curvature index values. The 'Folding Index' \n\
	  is determined by dividing the integral of the FI curvature function \n\
	  by 4*pi, and can also be reported by 'mris_anatomical_stats'. \n\
	 \n\
	  Three 'Intrinsic Curvuture Index' measures are reported. These are \n\
	  simply found by dividing the Gaussian curvature surface integrals \n\
	  by 4*pi. The ICIp is also reported by 'mris_anatomical_stats'. The \n\
	  INIn and ICIt are extensions of the Intrinsic Curvature, derived from \n\
	  the negative and natural Gaussian surface integrals -- again by \n\
	  dividing these integral results by 4*pi. \n\
	 \n\
	  Note that the ICIt (the integral of the Gaussian curvature divided \n\
	  by 4*pi) should be close to 1.000 as per the Gauss-Bonnet theorem \n\
	  which states that the total Gaussian curvature of any closed surface \n\
	  is 2*pi*Euler_number. For topologically correct surfaces, the \n\
	  Euler number is 2; thus the Gaussian integral is 4*pi and the ICIt \n\
	  is (in the ideal case) 1. \n\
 \n\
	FILTERING SURFACE DOMAINS \n\
 \n\
	  There are two main mechanisms for restricting the domain (or regions \n\
	  on the surface) where calculations occur. The simplest technique is \n\
	  by passing a label file with [-l <labelFileName>]. This \n\
	  <labelFileName> defines a surface region of interest, and all \n\
	  calculations will be constrained to this region. \n\
	 \n\
	  A second method is by specifying upper and lower filter bounds on the \n\
	  command line. Two techniques are available. The first filters accord- \n\
	  ing to the current curvature function value at a vertex of interest, \n\
	  and the second filters according to the Gaussian curvature value at \n\
	  a vertex of interest. \n\
 \n\
	  Filtering by the current curvature function is controlled by the \n\
	  [--lowPassFilter <value>] and [--highPassFilter <value>] flags, while \n\
	  filtering by vertex Gaussian curvature value is similarly specified \n\
	  by [--lowPassFilterGaussian <value>] and \n\
	  [--highPassFilterGaussian <value>]. If an additional argument, \n\
	  [--filterLabel <labelFileName>] is also given, the vertices tagged \n\
	  by the filter pattern will be saved to the <labelFileName>. \n\
 \n\
	  Note that the [-l <labelFileName>] and the high/low pass filters \n\
	  can be specified concurrently. If a <labelFileName> is given, then \n\
	  the low/high pass filters will operate only on the surface defined \n\
	  by the <labelFileName>. If a '--filterLabel <filename>' is also \n\
	  given, the filter <filename> will only contain vertices that lie \n\
	  within the original <labelFileName> *and* satisfy the filter pattern. \n\
 \n\
	  Specifying either a '-l <labelFileName>' or adding a filter pattern \n\
	  will result in some additional information being presented with the \n\
	  curvature report, specifically the ROI surface area, the ROI number \n\
	  number of vertices, and the ROI surface area percentage. \n\
		 	 \n\
	SURFACE INTEGRALS \n\
 \n\
	  The surface integrals for a given curvature map are filtered/modified \n\
	  in several ways. \n\
 \n\
    		- 'natural':	no change/filtering on the vertex curv values \n\
    		- 'abs':	the fabs(...) of each vertex curv value \n\
    		- 'pos':	process only the positive vertex curv values \n\
    		- 'neg':	process only the negative vertex curv values \n\
    	 \n\
	  In addition, the following qualifiers are also evaluated: \n\
    		- 'Mean':	the integral normalized to number of vertices \n\
    		- 'AreaNorm':	the integral value normalized to unit surface \n\
 \n\
	  The surface integral is defined as the sum across all the vertices \n\
	  of the particular curvature function at that vertex, weighted by the \n\
	  area of the vertex. \n\
 \n\
    OPTIONS \n\
 \n\
    [-a <numberOfAverages] \n\
 \n\
          Average the curvature <numberOfAverages> times. To save the resultant \n\
          curvature file, also specify '--writeCurvatureFiles -c 1', which \n\
          scales the resultant curvature by 1.0 and saves to *.scaled.crv. \n\
 \n\
    [-G] [--discrete] [--continuous] [--signedPrincipals] \n\
 \n\
          The '-G' flag triggers a series of derived curvature values that \n\
          are calculated from the current <FreeSurferSurface>. These curvatures \n\
          include the Gaussian, K, the mean, H, and the principal curvatures \n\
          K1 and K2. Additionally, several curvature maps based off these \n\
          values are also determined. \n\
 \n\
          Two calculation methods are available, a '--discrete' method and a \n\
          '--continuous' method. By default, '--discrete' is assumed. This \n\
          method is slower than the '--continuous' method, but more accurate. \n\
          The '--continuous' method is also used by 'mris_anatomical_stats', \n\
          however, it does suffer from occasionally large Gaussian outliers. \n\
 \n\
          Note that if both '--discrete' and '--continuous' are specified \n\
          together, the system will assume '--discrete'. \n\
 \n\
          In the '--discrete' mode, 'mris_curvature_stats' solves first for \n\
          K and H, and then solves a quadratic function to determine K1 and K2. \n\
          By default, K1 is assigned the fabs(max) solution, and K2 is \n\
          assigned the fabs(min) solution. By specifying '--signedPrincipals', \n\
          K1 is assigned the (signed) max and K2 the (signed) min. \n\
 \n\
          Thus, if the solution principal curves at a vertex are -10 and 0.5, \n\
          the default assignment would be: \n\
 \n\
                K1 = -10 and K2 = 0.5 \n\
 \n\
           but with '--signedPrincipals', this becomes: \n\
 \n\
                K1 = 0.5 and K2 = -10 \n\
 \n\
    [--vertexAreaWeigh]         [--vertexAreaNormalize] \n\
    [--vertexAreaWeighFrac]     [--vertexAreaNormalizeFrac] \n\
 \n\
	  If specified, will change the value of the curvature value \n\
	  at each point. The '--vertexAreaWeigh' will multiply each \n\
	  curvature value by the area of its vertex, and \n\
	  '--vertexAreaNormalize' will divide each curvature value by \n\
	  the area of its vertex. \n\
 \n\
          The {..}Frac flags use fractional operations, weighing or \n\
          normalizing by the fractional vertex area which is defined by \n\
          fractionalArea = (vertexArea/TotalSurfaceArea). \n\
 \n\
    [--postScale <scaleFactor>] \n\
	 \n\
	  If specified, scale the mean and areaNorm integrals by the amount \n\
	  <scaleFactor>. This can be useful in normalizing integral values. \n\
 \n\
    [--shapeIndex] \n\
 \n\
	  The 'Shape Index' curvature map can be problematic due to \n\
	  issues with some atan calculations. By default the shape index \n\
	  is not calculated. Use this flag to force shape index \n\
	  calculations. \n\
 \n\
    [-o <outputFileStem>] \n\
 \n\
          Save processing results to a series of files. This includes \n\
          condensed text output, histogram files (MatLAB friendly) \n\
          and curvature files. \n\
 \n\
          The actual files that are saved depends on which additional \n\
          calculation flags have been specified (i.e. normalisation, \n\
          Gaussian / Mean, scaling). \n\
 \n\
          In the case when a Gaussian/Mean calculation has been \n\
          performed, 'mris_curvature_stats' will act in a manner \n\
          similar to 'mris_curvature -w'.  Note though that though the \n\
          name of the curvature files that are created are different, \n\
          the contents are identical to 'mris_curvature -w' created files. \n\
 \n\
          All possible files that can be saved are: \n\
	  (where	O	= <outputFileStem> \n\
			H	= <hemisphere> \n\
			S	= <surface> \n\
			C	= <curvature>) \n\
 \n\
          <O> 			Log only a single  mean+-sigma. If several \n\
          			curvature files have been specified, log the \n\
				mean+-sigma across all the curvatures. Note \n\
				also that this file is *appended* for each \n\
				new run. \n\
          <OHSC>.log		Full output, i.e the output of each curvature \n\
				file mean +- sigma, as well as min/max and \n\
				surface integrals \n\
 \n\
          <OHS>.raw.hist  	Raw histogram file. By 'raw' is implied that \n\
				the curvature has not been further processed \n\
          			in any manner. \n\
          <OHS>.norm.hist  	Normalised histogram file \n\
          <OHS>.scaled.hist 	Scaled histogram file \n\
          <OHS>.K.hist  	Gaussian histogram file \n\
          <OHS>.H.hist  	Mean histogram file \n\
	  <OHS>.k1.hist		k1 histogram file \n\
	  <OHS>.k2.hist		k2 histogram file \n\
	  <OHS>.C.hist		C 'curvedness' histogram file \n\
	  <OHS>.S.hist		S 'sharpness' histogram file \n\
	  <OHS>.BE.hist		BE 'bending energy' histogram file \n\
	  <OHS>.SI.hist		SI 'shape index' histogram file \n\
	  <OHS>.FI.hist		FI 'folding index' histogram file \n\
 \n\
          <HS>.norm.crv   	Normalised curv file	('-n') \n\
          <HS>.scaled.crv  	Scaled curv file	('-c' or '-d' '-e') \n\
          <HS>.K.crv   		Gaussian curv file  	('-G') \n\
          <HS>.H.crv   		Mean curv file  	('-G') \n\
          <HS>.k1.crv   	k1 curv file  		('-G') \n\
          <HS>.k2.crv   	k2 curv file  		('-G') \n\
          <HS>.S.crv   		S curv file  		('-G') \n\
          <HS>.C.crv   		C curv file  		('-G') \n\
          <HS>.BE.crv   	BE curv file  		('-G') \n\
          <HS>.SI.crv   	SI curv file  		('-G') \n\
          <HS>.FI.crv   	FI curv file  		('-G') \n\
 \n\
	  (The above *.crv files can also be created with a \n\
	  '--writeCurvatureFiles' flag) \n\
 \n\
          Note that curvature files are saved to \n\
 \n\
			$SUBJECTS_DIR/<subjectname>/surf \n\
 \n\
	  and *not* to the current working directory. \n\
 \n\
    [-h <numberOfBins>] [-p <numberOfBins] \n\
 \n\
          If specified, prepare a histogram over the range of curvature \n\
          values, using <numberOfBins> buckets. These are dumped to \n\
          stdout. \n\
 \n\
          If '-p' is used, then the histogram is expressed as a \n\
          percentage. \n\
 \n\
          Note that this histogram, working off float values and float \n\
          boundaries, can suffer from rounding errors! There might be \n\
          instances when a very few (on average) curvature values might \n\
          not be sorted. \n\
 \n\
          The histogram behaviour can be further tuned with the \n\
          following: \n\
 \n\
    [-b <binSize>] [-i <binStartCurvature] [-j <binEndCurvature] \n\
 \n\
          These arguments are only processed iff a '-h <numberOfBins>' \n\
          has also been specified. By default, <binSize> is defined as \n\
 \n\
          (maxCurvature - minCurvature) / <numberOfBins> \n\
 \n\
          The '-b' option allows the user to specify an arbitrary \n\
          <binSize>. This is most useful when used in conjunction with \n\
          the '-i <binStartCurvature>' option, which starts the histogram \n\
          not at (minCurvature), but at <binStartCurvature>. So, if \n\
          a histogram reveals that most values seem confined to a very \n\
          narrow range, the '-b' and '-i' allow the user to 'zoom in' \n\
          to this range and expand. \n\
 \n\
          If <binStartCurvature> < (minCurvature), then regardless \n\
          of its current value, <binStartCurvature> = (minCurvature). \n\
          Also, if (<binStartCurvature> + <binSize>*<numberOfBins> >) \n\
          (maxCurvature), an error is raised and processing aborts. \n\
 \n\
          The '-j' allows the user to specify an optional end \n\
          value for the histogram. Using '-i' and '-j' together \n\
          are the most convenient ways to zoom into a region of interest \n\
          in a histogram. \n\
 \n\
    [-l <labelFileName>] [--regionalPercentages] \n\
 \n\
          Constrain statistics to the region defined in <labelFileName>. \n\
          If additionally, the [--regionalPercentages] flag is passed, \n\
	  the integral percentages are reported relative to the region, and \n\
	  not the whole brain surface (the default). \n\
 \n\
    [--highPassFilter <HPvalue>] [--lowPassFilter <LPvalue>] \n\
 \n\
	  When computing surface properties, only consider vertices where \n\
	  the current curvature function map satisfies the filter constraints, \n\
	  i.e. if [--highPassFilter <HPvalue>], then only process vertex if \n\
	  f_curv >= HPvalue. Similarly, if [--lowPassFilter <LPvalue] is passed \n\
	  then only process vertex if f_curv <= LPvalue. \n\
 \n\
    [--highPassFilterGaussian <HPvalue>] [--lowPassFilterGaussian <LPvalue>] \n\
 \n\
	  Same as above, but filter vertices according to their Gaussian \n\
	  curvature values. \n\
	 \n\
    [--filterLabel <labelFile>] \n\
 \n\
	  If any command line filters are specified, adding a [--filterLabel] \n\
	  argument will store the surface vertices that were processed in \n\
	  the FreeSurfer label file, <labelFile>. \n\
 \n\
    [-m] \n\
 \n\
          Output min / max information for the processed curvature. \n\
 \n\
    [-n] \n\
 \n\
          Normalise the curvature before computation. Normalisation \n\
          takes precedence over scaling, so if '-n' is specified \n\
          in conjunction with '-c' or '-smin'/'-smax' it will \n\
          override the effects of the scaling. \n\
 \n\
          If specified in conjunction with '-o <outputFileStem>' \n\
          will also create a curvature file with these values. \n\
 \n\
    [-s <summaryCondition>] \n\
 \n\
          Write out stats as <summaryCondition>. \n\
 \n\
    [-d <minCurvature> -e <maxCurvature>] \n\
 \n\
          Scale curvature values between <minCurvature> and \n\
          <maxCurvature>. If the minimum curvature is greater \n\
          than the maximum curvature, or if either is \n\
          unspecified, these flags are ignored. \n\
 \n\
          This scale computation takes precedence over '-c' scaling. \n\
 \n\
          Note also that the final scaling bounds might not correspond \n\
          to <minCurvature>... <maxCurvature> since values are scaled \n\
          across this range so as to preserve the original mean profile. \n\
 \n\
          If specified in conjunction with '-o <outputFileStem>' \n\
          will also create a curvature file with these values. \n\
 \n\
    [-c <factor>] \n\
 \n\
          Scale curvature values with <factor>. The mean of the \n\
          original curvature is conserved (and the sigma increases \n\
          with <factor>). \n\
 \n\
    [-version] \n\
 \n\
          Print out version number. \n\
 \n\
    [-z <vertexIndex>] \n\
 \n\
          Sets the curvature values at that index to zero. The \n\
          'raw' curvature, as well as the Gaussian and Mean curvatures \n\
          are set to zero, and min/max values are recomputed. \n\
 \n\
          This is useful in cases when outliers in the data (particularly \n\
          evident in Gaussian calcuations) skew mean and sigma values. \n\
 \n\
    [-q <maxUlps>] \n\
 \n\
          The <maxUlps> is used to toggle a more rigorous floating point \n\
          comparison operation in the histogram function. Comparing \n\
          float values for sorting into bins can at times fail due to \n\
          number precision issues. If, over the range of comparison \n\
          some curvature values are not sorted, add <maxUlps>. \n\
 \n\
          This adds extra function calls to AlmostEqual2sComplement(..) \n\
          for float comparisons and improves the general accuracy, at \n\
          a very slight performance penalty. \n\
 \n\
          You will most likely never have to use this argument, and is \n\
          for advanced use only. \n\
 \n\
    NOTES \n\
 \n\
          It is important to note that some combinations of the command \n\
          line parameters are somewhat meaningless, such as normalising \n\
          a 'sulc' curvature file (since it's normalised by definition). \n\
 \n\
    EXAMPLES \n\
 \n\
    $>mris_curvature_stats 801_recon rh curv \n\
 \n\
          For subject '801_recon', determine the mean+-sigma and surface \n\
	  integrals for the curvature file 'curv' on the right hemisphere. \n\
 \n\
    $>mris_curvature_stats -m 801_recon rh curv \n\
 \n\
          Same as above, but print the min/max curvature values \n\
          across the surface. \n\
 \n\
    $>mris_curvature_stats -h 20 -m 801_recon rh curv \n\
 \n\
          Same as above, and also print a histogram of curvature \n\
          values over the min/max range, using 20 bins. By replacing \n\
          the '-h' with '-p', print the histogram as a percentage. \n\
 \n\
    $>mris_curvature_stats -h 20 -b 0.01 -i -0.1 -m 801_recon rh curv \n\
    $>mris_curvature_stats -h 20 -i -0.1  -j 0.1 -m 801_recon rh curv \n\
 \n\
          Same as above, but this time constrain the histogram to the 20 \n\
          bins from -0.1 to 0.1, with a bin size of 0.01. \n\
 \n\
          Note that the count / percentage values are taken across the \n\
          total curvature range and not the constrained window defined \n\
          by the '-i' and '-b' arguments. \n\
 \n\
    $>mris_curvature_stats -G 801_recon rh \n\
 \n\
	  Calculate all the second order measures for the right hemisphere \n\
	  of subject '801_recon'. By default, the 'smoothwm' surface is \n\
	  selected. \n\
 \n\
    $>mris_curvature_stats -G -F inflated 801_recon rh \n\
 \n\
	  Same as above, but use the 'inflated' surface. \n\
 \n\
    $>mris_curvature_stats -h 10 -G -m 801_recon rh \n\
 \n\
          Same as above, with the addition of a histogram and min/max \n\
	  curvatures for all the second order measures. \n\
 \n\
    $>mris_curvature_stats -h 10 -G -m -o foo 801_recon rh \n\
 \n\
          Generate several output text files using the stem 'foo' that \n\
	  capture the min/max and histograms for each curvature processed. \n\
	  Also create new second order curvature files. \n\
 \n\
          In this case, the curvature files created are called: \n\
 \n\
          		rh.smoothwm.K.crv \n\
          		rh.smoothwm.H.crv \n\
          		rh.smoothwm.k1.crv \n\
          		rh.smoothwm.k2.crv \n\
          		rh.smoothwm.S.crv \n\
          		rh.smoothwm.C.crv \n\
          		rh.smoothwm.BE.crv \n\
 \n\
          and are saved to the $SUBJECTS_DIR/<subjectname>/surf directory. \n\
          These can be re-read by 'mris_curvature_stats' using \n\
 \n\
    $>mris_curvature_stats -m 801_recon rh 		\\ \n\
			smoothwm.K.crv			\\ \n\
			smoothwm.H.crv			\\ \n\
			smoothwm.k1.crv			\\ \n\
			smoothwm.k2.crv   		\\ \n\
			smoothwm.S.crv			\\ \n\
			smoothwm.C.crv			\\ \n\
			smoothwm.BE.crv			\\ \n\
 \n\
    ADVANCED EXAMPLES \n\
 \n\
          'mris_curvature_stats' can also provide some useful side \n\
          effects. Reading in a curvature, and applying any calculation \n\
          to it (scaling, gaussian, etc.) can result in data that \n\
          can be visualised quite well in a tool such as 'tksurfer'. \n\
 \n\
          Consider the normal curvature display in 'tksurfer', which \n\
          is usually quite dark due to the dynamic range of the display. \n\
          We can considerably improve the brightness by scaling a \n\
          curvature file and rendering the resultant in 'tksurfer'. \n\
 \n\
          First, take an arbitrary curvature, apply a scale factor, \n\
          and an output filestem: \n\
 \n\
    $>mris_curvature_stats --writeCurvatureFiles -c 10 801_recon rh curv \n\
 \n\
          This scales each curvature value by 10. A new curvature file \n\
          is saved in \n\
 \n\
          	$SUBJECTS_DIR/801_recon/surf/rh.smoothwm.scaled.crv \n\
 \n\
          Comparing the two curvatures in 'tksurfer' will clearly show \n\
          the scaled file as much brighter. \n\
 \n\
          Similarly, the Gaussian curvature can be processed, scaled, and \n\
          displayed, yielding very useful visual information. First \n\
          create and save the principal curvature based files: \n\
 \n\
    $>mris_curvature_stats --writeCurvatureFiles -G 801_recon rh curv \n\
 \n\
	  This command will create Gaussian and Mean curvature files in the \n\
          $SUBJECTS_DIR/<subjectName>/surf directory: \n\
 \n\
          		rh.smoothwm.K.crv \n\
          		rh.smoothwm.H.crv \n\
 \n\
          Now, process the created Gaussian with the scaled curvature: \n\
 \n\
    $>mris_curvature_stats --writeCurvatureFiles -c 10 	\\ \n\
			801_recon rh smoothwm.K.crv \n\
 \n\
          The final scaled Gaussian curvature is saved to (again in the \n\
	  $SUBJECTS_DIR/801_recon/surf directory): \n\
 \n\
          		rh.smoothwm.scaled.crv \n\
 \n\
          which is a much better candidate to view in 'tksurfer' than \n\
          the original Gaussian curvature file. \n\
 \n\
\n");

  fprintf(stdout, "%s", pch_synopsis);
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

int  comp(const void *p, const void *q)
{
  float *ptr1 = (float *)(p);
  float *ptr2 = (float *)(q);

  if (*ptr1 < *ptr2) {
    return -1;
  } else if (*ptr1 == *ptr2) {
    return 0;
  } else {
    return 1;
  }
}
