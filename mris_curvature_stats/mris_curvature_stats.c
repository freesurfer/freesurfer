/**
 * @file  mris_curvature_stats.c
 * @brief Determine the mean and std of curvature files, and much much more. 
 *
 * 'mris_curvature_stats' calculates statistics on curvature across a surface. It also
 * creates simple histograms on curvature distribution. Additionally, it can also
 * be used to analyze a surface and determine the principle curvatures and save to
 * extra curvature files, e.g. K, H, k1, k2.
 *
 */
/*
 * Original Author: Bruce Fischl / heavily hacked by Rudolph Pienaar
 * CVS Revision Info:
 *    $Author: rudolph $
 *    $Date: 2007/04/17 16:22:22 $
 *    $Revision: 1.26 $
 *
 * Copyright (C) 2002-2007,
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

#include <getopt.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"

#define  STRBUF  65536
#define  MAX_FILES  1000
#define  CE( x )  fprintf(stdout, ( x ))
#define  START_i   3

// Calculations performed on the curvature surface
typedef enum _secondOrderType {
  e_Raw  		= 0, 	// "Raw" (native) curvature - no calcuation
  e_Gaussian 		= 1, 	// Gaussian curvature 	= k1*k2
  e_Mean  		= 2, 	// Mean curvature	= 0.5*(k1+k2))
  e_K1  		= 3, 	// k1 curvature
  e_K2  		= 4, 	// k2 curvature
  e_S			= 5,	// "sharpness" 		= (k1-k2)^2
  e_C			= 6,	// "curvedness" 	= sqrt(0.5*(k1^2+k2^2))
  e_SI			= 7,	// "shape index"	= atan((k1+k2)/(k2-k1))
  e_BE			= 8,	// "bending energy"	= k1^2 + k2^2
  e_Normal 		= 9, 	// Normalised curvature
  e_Scaled 		= 10, 	// "Raw" scaled curvature
  e_ScaledTrans 	= 11 	// Scaled and translated curvature
} e_secondOrderType;

// Set of possible output curvature files
char* Gppch[] = {
                  "Raw",
                  "Gaussian",
                  "Mean",
                  "k1",
                  "k2",
		  "sharpness",
		  "curvedness",
		  "shapeIndex",
		  "bending energy",
                  "Normal",
                  "Scaled",
                  "ScaledTrans"
                };

// Output file prefixes
typedef enum _OFSP {
  e_None,   // No prefix
  e_Partial,  // A partial prefix (does not include the stem)
  e_Full,   // Full prefix - including stem
} e_OFSP;

static char vcid[] =
  "$Id: mris_curvature_stats.c,v 1.26 2007/04/17 16:22:22 rudolph Exp $";

int   main(int argc, char *argv[]) ;

static int   get_option(int argc, char *argv[]) ;
static void  usage_exit(void) ;
static void  print_usage(void) ;
static void  print_help(void) ;
static void  print_version(void) ;

// Function Prototypes

// Simple functions of principle curvatures
float f_sharpnessCurvature(float af_k1, float af_k2) 
	{return ((af_k1 - af_k2)*(af_k1 - af_k2));}
float f_bendingEnergyCurvature(float af_k1, float af_k2) 
	{return (af_k1*af_k1 + af_k2*af_k2);}
float f_curvednessCurvature(float af_k1, float af_k2) 
	{return (sqrt(0.5*(af_k1*af_k1 + af_k2*af_k2)));}
float f_shapeIndexCurvature(float af_k1, float af_k2) 
	{return (af_k1 == af_k2 ? 0 : atan((af_k1+af_k2)/(af_k2 - af_k1)));}

// Simple functions on specific vertices
float f_absCurv(VERTEX*			pv) {return (fabs(pv->curv));}
float f_pass(VERTEX*			pv) {return (pv->curv);}
short f_allVertices(VERTEX*		pv) {return 1;}
short f_greaterEqualZero(VERTEX*	pv) {return (pv->curv >= 0 ? 1 : 0);}
short f_lessThanZero(VERTEX*		pv) {return (pv->curv <0 ? 1 : 0);}

int MRISusePrincipleCurvatureFunction(
	MRI_SURFACE*		pmris, 
	float 			(*f)(float k1, float k2)
);

double MRIScomputeSurfaceIntegral(
	MRI_SURFACE* 	pmris, 
	int* 		p_verticesCounted,
	short		(*fcond)(VERTEX*	pv),
	float		(*fv)	(VERTEX*	pv)	
);

int    comp(const void *p, const void *q);    // Compare p and q for qsort()

void  histogram_wrapper(
  MRIS*   amris,
  e_secondOrderType aesot
);

void  histogram_create(
  MRIS*   amris_curvature,
  float   af_minCurv,
  double   af_binSize,
  int   abins,
  float*   apf_histogram,
  e_secondOrderType aesot
);

void  OFSP_create(
  char*   apch_prefix,
  char*   apch_suffix,
  char*   apch_curv,
  e_OFSP   ae_OFSP
);

void   secondOrderParams_print(
  MRIS*   apmris,
  e_secondOrderType aesot,
  int   ai
);

void  outputFileNames_create(
  char*   apch_curv
);

void  outputFiles_open(void);

void  outputFiles_close(void);

int  MRISminMaxCurvaturesSearchSOT(
  MRI_SURFACE*  apmris,
  int*   ap_vertexMin,
  int*   ap_vertexMax,
  float*   apf_min,
  float*   apf_max,
  e_secondOrderType aesot
);

int  MRISminMaxCurvaturesSearch(
  MRI_SURFACE*  apmris,
  int*   ap_vertexMin,
  int*   ap_vertexMax,
  float*   apf_min,
  float*   apf_max
) {
  MRISminMaxCurvaturesSearchSOT(  apmris,
                                  ap_vertexMin,
                                  ap_vertexMax,
                                  apf_min,
                                  apf_max,
                                  e_Raw);
  return(NO_ERROR);
};

int  MRISminMaxCurvatureIndicesLookup(
  MRI_SURFACE*  apmris,
  int*   ap_vertexMin,
  int*   ap_vertexMax
);

int  MRISvertexCurvature_set(
  MRI_SURFACE*  apmris,
  int   aindex,
  float   af_val
);

int    MRISzeroCurvature(
  MRI_SURFACE*  apmris
);

int  MRISuseK1Curvature(
  MRI_SURFACE*  mris
);

int  MRISuseK2Curvature(
  MRI_SURFACE*  mris
);

// static int	longOpts_process(
//   int argc, 
//   char** ppch_argv
// );


// Global variables

char*  Progname ;
char*  hemi;


// static int Gb_statsTableShow	= 0;
// static struct option Gs_longOptions[] =
// {
//   /* These options set a flag. */
//   {"statsTable",		no_argument,		&Gb_statsTableShow, 1},
//   /* These options don't set a flag.
//      We distinguish them by their indices. */
//   {"lowPassFilterGaussian",	required_argument, 	0, 'L'},
//   {"lowPassFilter",		required_argument, 	0, 'l'},
//   {"highPassFilterGaussian",	required_argument, 	0, 'H'},
//   {"highPassFilter",		required_argument, 	0, 'h'},
//   {0, 0, 0, 0}
// };

static int  navgs    = 0 ;
static int  normalize_flag   = 0 ;
static int  condition_no   = 0 ;
static int  stat_flag   = 0 ;
static char* label_name   = NULL ;
static char* output_fname   = NULL ;
static char surf_name[STRBUF];

// Additional global variables (prefixed by 'G') added by RP.
// Flags are prefixed with Gb_ and are used to track
// user spec'd command line flags.
static int  Gb_minMaxShow  = 0;
static int Gb_histogram  = 0;
static int Gb_histogramPercent = 0;
static int Gb_binSizeOverride = 0;
static double Gf_binSize  = 0.;
static int Gb_histStartOverride = 0;
static float Gf_histStart  = 0.;
static int Gb_histEndOverride = 0;
static float Gf_histEnd  = 0.;
static int Gb_gaussianAndMean = 0;
static int Gb_output2File  = 0;
static int Gb_scale  = 0;
static int Gb_scaleMin  = 0;
static int Gb_scaleMax  = 0;
static int Gb_zeroVertex  = 0;
static int G_zeroVertex  = 0;
static int  G_nbrs    = 2;
static int G_bins   = 1;
static int Gb_maxUlps  = 0;
static int G_maxUlps  = 0;

static int	Gb_lowPassFilter		= 0;
static float	Gf_lowPassFilter		= 0.;
static int	Gb_lowPassFilterGaussian	= 0;
static float	Gf_lowPassFilterGaussian	= 0.;
static int	Gb_highPassFilter		= 0;
static float	Gf_highPassFilter		= 0.;
static int	Gb_highPassFilterGaussian	= 0;
static float	Gf_highPassFilterGaussian	= 0.;

// All possible output file name and suffixes
static	short	Gb_writeCurvatureFiles		= 0;
static 	char 	Gpch_log[STRBUF];
static 	char 	Gpch_logS[]  			= "log";
static 	FILE* 	GpFILE_log  			= NULL;
static 	char 	Gpch_allLog[STRBUF];
static 	char 	Gpch_allLogS[]  		= "log";
static 	FILE* 	GpFILE_allLog  = NULL;
static 	char 	Gpch_rawHist[STRBUF];
static 	char 	Gpch_rawHistS[]  		= "raw.hist";
static 	FILE* 	GpFILE_rawHist  		= NULL;
static 	char 	Gpch_normCurv[STRBUF];
static 	char 	Gpch_normHist[STRBUF];
static 	char 	Gpch_normHistS[] 		= "norm.hist";
static 	FILE* 	GpFILE_normHist  		= NULL;
static 	char 	Gpch_normCurv[STRBUF];
static 	char 	Gpch_normCurvS[] 		= "norm.crv";
static 	char 	Gpch_KHist[STRBUF];
static 	char 	Gpch_KHistS[]  			= "K.hist";
static 	FILE* 	GpFILE_KHist  			= NULL;
static 	char 	Gpch_KCurv[STRBUF];
static 	char 	Gpch_KCurvS[]  			= "K.crv";
static 	char 	Gpch_HHist[STRBUF];
static 	char 	Gpch_HHistS[]  			= "H.hist";
static 	FILE* 	GpFILE_HHist  			= NULL;
static 	char 	Gpch_HCurv[STRBUF];
static 	char 	Gpch_HCurvS[]  			= "H.crv";
static 	char 	Gpch_scaledHist[STRBUF];
static 	char 	Gpch_scaledHistS[] 		= "scaled.hist";
static 	FILE* 	GpFILE_scaledHist 		= NULL;
static 	char 	Gpch_scaledCurv[STRBUF];
static 	char 	Gpch_scaledCurvS[] 		= "scaled.crv";

static char Gpch_K1Hist[STRBUF];
static char Gpch_K1HistS[]  = "K1.hist";
static FILE* GpFILE_K1Hist  = NULL;
static char Gpch_K1Curv[STRBUF];
static char Gpch_K1CurvS[]  = "K1.crv";

static char Gpch_K2Hist[STRBUF];
static char Gpch_K2HistS[]  = "K2.hist";
static FILE* GpFILE_K2Hist  = NULL;
static char Gpch_K2Curv[STRBUF];
static char Gpch_K2CurvS[]  = "K2.crv";

static char Gpch_SHist[STRBUF];
static char Gpch_SHistS[]  = "S.hist";
static FILE* GpFILE_SHist  = NULL;
static char Gpch_SCurv[STRBUF];
static char Gpch_SCurvS[]  = "S.crv";

static char Gpch_CHist[STRBUF];
static char Gpch_CHistS[]  = "C.hist";
static FILE* GpFILE_CHist  = NULL;
static char Gpch_CCurv[STRBUF];
static char Gpch_CCurvS[]  = "C.crv";

static char Gpch_BEHist[STRBUF];
static char Gpch_BEHistS[]  = "BE.hist";
static FILE* GpFILE_BEHist  = NULL;
static char Gpch_BECurv[STRBUF];
static char Gpch_BECurvS[]  = "BE.crv";

static char Gpch_SIHist[STRBUF];
static char Gpch_SIHistS[]  = "SI.hist";
static FILE* GpFILE_SIHist  = NULL;
static char Gpch_SICurv[STRBUF];
static char Gpch_SICurvS[]  = "SI.crv";

// These are used for tabular output
const int G_leftCols  = 20;
const int G_rightCols  = 20;

// Mean / sigma tracking and scaling
static double  Gpf_means[MAX_FILES] ;
static double Gf_mean   = 0.;
static double Gf_sigma  = 0.;
static double Gf_n   = 0.;
static double Gf_total  = 0.;
static double Gf_total_sq  = 0.;
static double Gf_scaleFactor  = 1.;
static double Gf_scaleMin  = 0.;
static double Gf_scaleMax  = 0.;

static int which_norm = NORM_MEAN;

int
main(int argc, char *argv[]) {
  char          **av, fname[STRBUF], *sdir ;
  char          *subject_name, *curv_fname ;
  int           ac, nargs, i ;
  MRI_SURFACE   *mris ;

  char  pch_text[65536];
  int  vmin    = -1;
  int  vmax    = -1;

//   int 	ret	= 0;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv,
                                 "$Id: mris_curvature_stats.c,v 1.26 2007/04/17 16:22:22 rudolph Exp $", "$Name:  $");
//   ret = longOpts_process(argc, argv);
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  sdir = getenv("SUBJECTS_DIR") ;
  if (!sdir)
    ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in environment.\n",Progname);
  ac = argc ;
  av = argv ;
  strcpy(surf_name, "-x");
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (!strcmp(surf_name, "-x"))
    strcpy(surf_name, "orig") ;

  if (argc < 4)
    usage_exit() ;

  /* parameters are
     1. subject name
     2. hemisphere
     3. specific curvature file to analyze
  */
  subject_name = argv[1] ;
  hemi = argv[2] ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, surf_name) ;
  printf("Reading surface from %s...\n", fname);
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;


  if (label_name) {
    LABEL *area ;
    area = LabelRead(subject_name, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s",
                Progname, label_name) ;
    LabelRipRestOfSurface(area, mris) ;
    LabelFree(&area) ;
  }
  if (label_name)
    fprintf(stdout, "%s: ", label_name) ;

  // Process all the command-line spec'd curvature files
  for (Gf_n= Gf_total_sq = Gf_total = 0.0, i = START_i ; i < argc ; i++) {
    curv_fname = argv[i] ;

    outputFiles_close();
    outputFileNames_create(curv_fname);
    outputFiles_open();

    MRIScomputeMetricProperties(mris);
    if (MRISreadCurvatureFile(mris, curv_fname) != NO_ERROR) {
      fprintf(stderr,
              "\n\t\t***WARNING!***\n");
      fprintf(stderr, "\tSome error has occurred while reading '%s.%s'.\n",
              hemi, curv_fname);
      fprintf(stderr, "\tThis might be due a vertex incompatibility\n");
      fprintf(stderr, "\tbetween the surface '%s.%s' and curv '%s.%s'.\n",
              hemi, surf_name, hemi, curv_fname);
      fprintf(stderr, "\n\tYou might be able to correct this by re-running\n");
      fprintf(stderr, "\t'mris_make_surfaces' on this dataset.\n\n");
      fprintf(stderr, "\tSkipping this (and any remaining) curvature files.\n");
      fprintf(stderr, "\tAny measurements / calcuations that depend on the\n");
      fprintf(stderr, "\tcurvature file will be skipped.\n");
      fprintf(stderr, "\t\t*************\n\n");
      break;
      /*      ErrorExit(ERROR_BADFILE,"%s: could not read curvature file %s.\n",
                      Progname, curv_fname);*/
    }

    if (Gb_zeroVertex)
      MRISvertexCurvature_set(mris, G_zeroVertex, 0);

    if (Gb_scale) {
      MRISscaleCurvature(mris, Gf_scaleFactor);
      if (Gpch_scaledCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, Gpch_scaledCurv);
    }

    if (Gb_scaleMin && Gb_scaleMax) {
      MRISscaleCurvatures(mris, Gf_scaleMin, Gf_scaleMax);
      if (Gpch_scaledCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, Gpch_scaledCurv);
    }

    if (normalize_flag) {
      MRISnormalizeCurvature(mris,which_norm);
      if (Gpch_normCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, Gpch_normCurv);
    }

    MRISaverageCurvatures(mris, navgs) ;
    Gf_mean = MRIScomputeAverageCurvature(mris, &Gf_sigma) ;
    sprintf(pch_text, "Raw Curvature     (using '%s.%s'):", hemi, curv_fname);
    fprintf(stdout, "%-50s", pch_text);
    fprintf(stdout, "%20.4f +- %2.4f\n", Gf_mean, Gf_sigma) ;
    sprintf(pch_text, "Vertex Statistics (using '%s.%s'):", hemi, surf_name);
    fprintf(stdout, "%-50s", pch_text);
    fprintf(stdout, "%20.4f +- %2.4f\n",
            mris->avg_vertex_dist, mris->std_vertex_dist);
    if (GpFILE_log)
      fprintf(GpFILE_allLog, "mean/sigma = %20.4f +- %2.4f\n", Gf_mean, Gf_sigma);

    if (Gb_minMaxShow) {
      MRISminMaxCurvatureIndicesLookup(mris, &vmin, &vmax);
      fprintf(stdout,
              "%*s%20.6f\tvertex = %d\n%*s%20.6f\tvertex = %d\n",
              G_leftCols, "Raw min = ", mris->min_curv, vmin,
              G_leftCols, "Raw max = ", mris->max_curv, vmax);
      if (GpFILE_allLog)
        fprintf(GpFILE_allLog, "min = %f\nmax = %f\n",
                mris->min_curv, mris->max_curv);
      if (Gb_histogram) histogram_wrapper(mris, e_Raw);
    }

    Gpf_means[i-START_i] = Gf_mean ;
    Gf_total += Gf_mean ;
    Gf_total_sq += Gf_mean*Gf_mean ;
    Gf_n+= 1.0 ;
  }

  // Should we calculate all the principle curvature based curvatures? 
  // This is a surface-based calculation, and does not depend on the 
  // curvature processing loop - thus the 'curv' input is irrelevant
  if (Gb_gaussianAndMean) {
    fprintf(stderr,
	    "%-50s",
            "Calculating second fundamental form...");

    /* Principle curvature calculations */
    MRISsetNeighborhoodSize(mris, G_nbrs);
    MRIScomputeSecondFundamentalForm(mris);
    fprintf(stderr, "%20s\n", "[ ok ]");

    secondOrderParams_print(mris, e_Gaussian, i);
    if (Gpch_KCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_KCurv);
    secondOrderParams_print(mris, e_Mean,   i);
    if (Gpch_HCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_HCurv);
    secondOrderParams_print(mris, e_K1, i);
    if (Gpch_K1Curv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_K1Curv);
    secondOrderParams_print(mris, e_K2, i);
    if (Gpch_K2Curv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_K2Curv);
    secondOrderParams_print(mris, e_S, i);
    if (Gpch_SCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_SCurv);
    secondOrderParams_print(mris, e_C, i);
    if (Gpch_CCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_CCurv);
    secondOrderParams_print(mris, e_BE, i);
    if (Gpch_BECurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_BECurv);
    secondOrderParams_print(mris, e_SI, i);
    if (Gpch_SICurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_SICurv);
  }


  if (Gf_n> 1.8) {
    Gf_mean = Gf_total / Gf_n;
    Gf_sigma = sqrt(Gf_total_sq/Gf_n- Gf_mean*Gf_mean) ;
    fprintf(stdout, "\nMean across %d curvatures: %8.4e +- %8.4e\n",
            (int) Gf_n, Gf_mean, Gf_sigma) ;
  }
  //else
  //  Gf_mean = Gf_sigma = 0.0 ;


  MRISfree(&mris) ;
  if (output_fname) {
    // This code dumps only the *final* mean/sigma values to the summary
    // log file. For cases where only 1 curvature file has been processed,
    // the mean/sigma of this file is written. For multiple files, the
    // mean/sigma over *all* of the files is written.
    if (label_name)
      fprintf(GpFILE_log, "%s: ", label_name) ;
    fprintf(GpFILE_log, "%8.4e +- %8.4e\n", Gf_mean, Gf_sigma) ;
  }
  outputFiles_close();
  exit(0) ;
  return(0) ;  /* for ansi */
}

void secondOrderParams_print(
  MRIS*   apmris,
  e_secondOrderType aesot,
  int   ai
) {
  //
  // PRECONDITIONS
  // o amris must have been prepared with a call to
  //   MRIScomputeSecondFundamentalForm(...)
  //
  // POSTCONDITIONS
  //  o Depending on <esot>, data is printed to stdout (and
  //   output file, if spec'd).
  //

  char 	pch_out[65536];
  char 	pch_text[65536];
  char 	pch_type[32];
  char 	pch_min[32];
  char 	pch_max[32];
  float f_max   		= 0.;
  float f_min   		= 0.;
  float f_maxExplicit  		= 0.;
  float f_minExplicit  		= 0.;
  int  	vmax   			= -1;
  int  	vmin   			= -1;
  float	f_SInatural		= 0.0;	int SInaturalVertices	= 0;
  float f_SIabs			= 0.0; 	int SIabsVertices	= 0;
  float f_SIpos			= 0.0;	int SIposVertices	= 0;
  float f_SIneg			= 0.0;  int SInegVertices	= 0;

  sprintf(pch_type, " ");
  if (Gb_zeroVertex)
    MRISvertexCurvature_set(apmris, G_zeroVertex, 0);
  switch (aesot) {
  case e_Gaussian:
    sprintf(pch_type, "K");
    MRISuseGaussianCurvature(apmris);
    sprintf(pch_out, "\nK");
    f_min = apmris->Kmin;
    f_max = apmris->Kmax;
    MRISminMaxCurvaturesSearch(apmris, &vmin, &vmax, &f_minExplicit, &f_maxExplicit);
    if (f_min != f_minExplicit) {
      printf("\tLookup   min: %f\n", f_min);
      printf("\tExplicit min: %f\tvertex = %d\n", f_minExplicit, vmin);
      f_min = f_minExplicit;
      apmris->Kmin   = f_minExplicit;
      apmris->min_curv = f_minExplicit;
    }
    if (f_max != f_maxExplicit) {
      printf("\tLookup   max: %f\n", f_max);
      printf("\tExplicit max: %f\tvertex = %d\n", f_maxExplicit, vmax);
      f_max = f_maxExplicit;
      apmris->Kmax   = f_maxExplicit;
      apmris->max_curv = f_maxExplicit;
    }
    break;
  case e_Mean:
    sprintf(pch_type, "H");
    MRISuseMeanCurvature(apmris);
    sprintf(pch_out, "\nH");
    f_min = apmris->Hmin;
    f_max = apmris->Hmax;
    MRISminMaxCurvaturesSearch(apmris, &vmin, &vmax, &f_minExplicit, &f_maxExplicit);
    if (f_min != f_minExplicit) {
      printf("\tLookup   min: %f\n", f_min);
      printf("\tExplicit min: %f\tvertex = %d\n", f_minExplicit, vmin);
      f_min = f_minExplicit;
      apmris->Hmin   = f_minExplicit;
      apmris->min_curv = f_minExplicit;
    }
    if (f_max != f_maxExplicit) {
      printf("\tLookup   max: %f\n", f_max);
      printf("\tExplicit max: %f\tvertex = %d\n", f_maxExplicit, vmax);
      f_max = f_maxExplicit;
      apmris->Hmax   = f_maxExplicit;
      apmris->max_curv = f_maxExplicit;
    }
    break;
  case e_K1:
    sprintf(pch_type, "k1");
    MRISuseK1Curvature(apmris);
    sprintf(pch_out, "\nk1");
    f_min = apmris->min_curv;
    f_max = apmris->max_curv;
    MRISminMaxCurvaturesSearch(apmris, &vmin, &vmax, &f_minExplicit, &f_maxExplicit);
    if (f_min != f_minExplicit) {
      printf("\tLookup   min: %f\n", f_min);
      printf("\tExplicit min: %f\tvertex = %d\n", f_minExplicit, vmin);
      f_min = f_minExplicit;
      apmris->Hmin   = f_minExplicit;
      apmris->min_curv = f_minExplicit;
    }
    if (f_max != f_maxExplicit) {
      printf("\tLookup   max: %f\n", f_max);
      printf("\tExplicit max: %f\tvertex = %d\n", f_maxExplicit, vmax);
      f_max = f_maxExplicit;
      apmris->Hmax   = f_maxExplicit;
      apmris->max_curv = f_maxExplicit;
    }
    break;
  case e_K2:
    sprintf(pch_type, "k2");
    MRISuseK2Curvature(apmris);
    sprintf(pch_out, "\nk2");
    f_min = apmris->min_curv;
    f_max = apmris->max_curv;
    MRISminMaxCurvaturesSearch(apmris, &vmin, &vmax, &f_minExplicit, &f_maxExplicit);
    if (f_min != f_minExplicit) {
      printf("\tLookup   min: %f\n", f_min);
      printf("\tExplicit min: %f\tvertex = %d\n", f_minExplicit, vmin);
      f_min = f_minExplicit;
      apmris->Hmin   = f_minExplicit;
      apmris->min_curv = f_minExplicit;
    }
    if (f_max != f_maxExplicit) {
      printf("\tLookup   max: %f\n", f_max);
      printf("\tExplicit max: %f\tvertex = %d\n", f_maxExplicit, vmax);
      f_max = f_maxExplicit;
      apmris->Hmax   = f_maxExplicit;
      apmris->max_curv = f_maxExplicit;
    }
    break;
  case e_S:
    sprintf(pch_type, "S");
    MRISusePrincipleCurvatureFunction(apmris, f_sharpnessCurvature);
    sprintf(pch_out, "\nS");
    f_min = apmris->min_curv;
    f_max = apmris->max_curv;
    break;
  case e_C:
    sprintf(pch_type, "C");
    MRISusePrincipleCurvatureFunction(apmris, f_curvednessCurvature);
    sprintf(pch_out, "\nC");
    f_min = apmris->min_curv;
    f_max = apmris->max_curv;
    break;
  case e_BE:
    sprintf(pch_type, "BE");
    MRISusePrincipleCurvatureFunction(apmris, f_bendingEnergyCurvature);
    sprintf(pch_out, "\nBE");
    f_min = apmris->min_curv;
    f_max = apmris->max_curv;
    break;
  case e_SI:
    sprintf(pch_type, "SI");
    MRISusePrincipleCurvatureFunction(apmris, f_shapeIndexCurvature);
    sprintf(pch_out, "\nSI");
    f_min = apmris->min_curv;
    f_max = apmris->max_curv;
    break;
  case e_Raw:
  case e_Normal:
  case e_Scaled:
  case e_ScaledTrans:
    break;
  }

  f_SInatural	= MRIScomputeSurfaceIntegral(
				apmris, &SInaturalVertices,
				f_allVertices,
				f_pass);
  f_SIabs	= MRIScomputeSurfaceIntegral(
				apmris, &SIabsVertices,
				f_allVertices,
				f_absCurv);
  f_SIpos	= MRIScomputeSurfaceIntegral(
				apmris, &SIposVertices,
				f_greaterEqualZero,
				f_absCurv);
  f_SIneg	= MRIScomputeSurfaceIntegral(
				apmris, &SInegVertices,
				f_lessThanZero,
				f_absCurv);

  sprintf(pch_min, "%s min =", pch_type);
  sprintf(pch_max, "%s max =", pch_type);

  Gf_mean = MRIScomputeAverageCurvature(apmris, &Gf_sigma);
  sprintf(pch_text, "%s Curvature (using '%s.%s'):",
          pch_out, hemi, surf_name);
  fprintf(stdout, "%-50s", pch_text);
  fprintf(stdout, " %10.4f +- %2.4f\n", Gf_mean, Gf_sigma);
  fprintf(stdout, "%10s%-40s", pch_type, " Total number of vertices:");
  fprintf(stdout, "%10.4d\n", apmris->nvertices);
  fprintf(stdout, "%10s%-40s", pch_type, " Average Vertex Area:");
  fprintf(stdout, "%10.4f\n", apmris->avg_vertex_area);
  fprintf(stdout, "%10s%-40s", pch_type, " Normalized Natural Surface Integral:");
  fprintf(stdout, "%10.4f across %d (%05.2f%s) vertices\n", 
			f_SInatural, SInaturalVertices, 
			100 * (float)SInaturalVertices / apmris->nvertices, "%");
  fprintf(stdout, "%10s%-40s", pch_type, " Normalized Rectified Surface Integral:");
  fprintf(stdout, "%10.4f across %d (%05.2f%s) vertices\n", 
			f_SIabs, SIabsVertices, 
			100 * (float)SIabsVertices / apmris->nvertices, "%");
  fprintf(stdout, "%10s%-40s", pch_type, " Normalized Positive Surface Integral:");
  fprintf(stdout, "%10.4f across %d (%05.2f%s) vertices\n", 
			f_SIpos, SIposVertices, 
			100 * (float)SIposVertices / apmris->nvertices, "%");
  fprintf(stdout, "%10s%-40s", pch_type, " Normalized Negative Surface Integral:");
  fprintf(stdout, "%10.4f across %d (%05.2f%s) vertices\n", 
			f_SIneg, SInegVertices, 
			100 * (float)SInegVertices / apmris->nvertices, "%");

  if (GpFILE_allLog)
    fprintf(GpFILE_allLog, "mean/sigma = %20.4f +- %2.4f\n", Gf_mean, Gf_sigma);
  if (Gb_minMaxShow) {
    fprintf(stdout,
            "%*s%20.6f\tvertex = %d\n%*s%20.6f\tvertex = %d\n",
            G_leftCols, pch_min, f_min, vmin,
            G_leftCols, pch_max, f_max, vmax);
    if (GpFILE_allLog)
      fprintf(GpFILE_allLog, "min = %f\nmax = %f\n",
              f_min, f_max);
  }

  if (Gb_histogram) histogram_wrapper(apmris, aesot);
  Gpf_means[ai-START_i]  = Gf_mean;
  Gf_total    += Gf_mean;
  Gf_total_sq   += Gf_mean*Gf_mean ;
  Gf_n   += 1.0;

}

int
MRISvertexCurvature_set(
  MRI_SURFACE*  apmris,
  int   aindex,
  float   af_val
) {
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

  VERTEX* pvertex;
  int  vno;
  float f_maxCurv = 0.;
  float  f_minCurv = 0.;
  float  f_maxCurvK = 0.;
  float  f_minCurvK = 0.;
  float  f_maxCurvH = 0.;
  float  f_minCurvH = 0.;

  if (aindex > apmris->nvertices)
    ErrorExit(ERROR_SIZE, "%s: target vertex is out of range.", Progname);

  apmris->vertices[aindex].curv = af_val;
  apmris->vertices[aindex].K  = af_val;
  apmris->vertices[aindex].H  = af_val;

  f_maxCurv  = f_minCurv  = apmris->vertices[0].curv;
  f_maxCurvK = f_minCurvK = apmris->vertices[0].K;
  f_maxCurvH = f_minCurvH = apmris->vertices[0].H;
  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    pvertex = &apmris->vertices[vno] ;
    if (pvertex->ripflag)
      continue;
    if (pvertex->curv > f_maxCurv)  f_maxCurv  = pvertex->curv;
    if (pvertex->curv < f_minCurv)  f_minCurv  = pvertex->curv;
    if (pvertex->K    > f_maxCurvK) f_maxCurvK = pvertex->K;
    if (pvertex->K    < f_minCurvK) f_minCurvK = pvertex->K;
    if (pvertex->H    > f_maxCurvH) f_maxCurvH = pvertex->H;
    if (pvertex->H    < f_minCurvH) f_minCurvH = pvertex->H;
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
  MRI_SURFACE*  apmris,
  int*   ap_vertexMin,
  int*   ap_vertexMax,
  float*   apf_min,
  float*   apf_max,
  e_secondOrderType aesot
) {
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
    if (pvertex->ripflag)
      continue;
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
) {
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
    if (pvertex->ripflag)
      continue;
    if (pvertex->curv == f_min)
      *ap_vertexMin = vno;
    if (pvertex->curv == f_max)
      *ap_vertexMax = vno;
  }
  return(NO_ERROR);
}

int
MRISuseK1Curvature(MRI_SURFACE *mris) {
  int    vno ;
  VERTEX *vertex ;

  float  f_min =  mris->vertices[0].curv;
  float  f_max = f_min;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->curv = vertex->k1 ;
    if (vertex->curv < f_min) f_min = vertex->curv;
    if (vertex->curv > f_max) f_max = vertex->curv;
  }

  mris->min_curv = f_min ;
  mris->max_curv = f_max;
  return(NO_ERROR) ;
}

int
MRISuseK2Curvature(MRI_SURFACE *mris) {
  int    vno ;
  VERTEX *vertex ;

  float  f_min =  mris->vertices[0].curv;
  float  f_max = f_min;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->curv = vertex->k2 ;
    if (vertex->curv < f_min) f_min = vertex->curv;
    if (vertex->curv > f_max) f_max = vertex->curv;
  }

  mris->min_curv = f_min ;
  mris->max_curv = f_max;
  return(NO_ERROR) ;
}

int
MRISusePrincipleCurvatureFunction(
	MRI_SURFACE*		pmris, 
	float 			(*f)(float k1, float k2)) {
  //
  // PRECONDITIONS
  //	o The principle curvatures k1 and k2 for each vertex point on the 
  //	  surface have been defined.
  //
  // POSTCONDITIONS
  // 	o Each vertex 'curv' value is replaced with the result from the
  //	  the (*f)(k1, k2) function.
  //	o Surface min and max values are set appropriately.
  //

  int    	vno ;
  VERTEX*	pvertex ;
  float		f_k1;
  float		f_k2;

  float  f_min =  pmris->vertices[0].curv;
  float  f_max = f_min;

  for (vno = 0 ; vno < pmris->nvertices ; vno++) {
    pvertex = &pmris->vertices[vno] ;
    if (pvertex->ripflag)
      continue ;
    f_k1		= pvertex->k1;
    f_k2		= pvertex->k2;
    pvertex->curv 	= (*f)(f_k1, f_k2);
    if (pvertex->curv < f_min) f_min = pvertex->curv;
    if (pvertex->curv > f_max) f_max = pvertex->curv;
  }

  pmris->min_curv = f_min ;
  pmris->max_curv = f_max;
  return(NO_ERROR) ;
}

double
MRIScomputeSurfaceIntegral(
	MRI_SURFACE* 	pmris, 
	int* 		p_verticesCounted,
	short		(*fcond)(VERTEX*	pv),
	float		(*fv)	(VERTEX*	pv)	
) {
  //
  // DESCRIPTION
  //	This function computes a surface integral across the 
  //	vertex curvature. The resultant sum is normalized to the
  //	area of the vertices that were counted.
  //
  //	Before processing a particular vertex curvature, the vertex
  //	is first checked by the (*fcond) function for eligibility. 
  //	
  //	The actual curvature that is integrated can also be shaped
  //	by the (*fcurv) function.
  //
  //	Global low pass filtering (on actual curvature values and
  //	Gaussian test is always implemented). 
  //
  // PRECONDITIONS
  //	o The vertex curvature 'curv' contains the particular
  //	  function curvature of interest.
  //	o The (*fcond) is a conditional function that checks the
  //	  curvature value for processing. Such tests can for example
  //	  only trigger on negative curvatures, or only positive, etc.
  //	o The (*fcurv) is an additional function of the vertex curvature
  //	  itself - often this will return the absolute value of the 
  //	  curvature.
  //
  // POSTCONDITIONS
  //	o The integral of the curvature value over the surface is returned.
  //	o If the integral could not be computed, -1 is returned.
  //	o If the global <Gb_lowPassFilter> is set, only count vertices
  //	  where abs 'curv' is less than <Gb_lowPassFilter>.
  //	o If the global <Gb_lowPassFilterGaussian> is set, only count
  //	  vertices where the abs Gaussian at that vertex is less than
  //	  <Gb_lowPassFilterGaussian>.
  //

  VERTEX*	pv ;
  int       	vno ;
  double    	f_total, f_n ;
  short		b_canCount	= 1;
  float		f_ret		= -1.0;

  *p_verticesCounted = 0;
  for (f_n = f_total = 0.0, vno = 0 ; vno < pmris->nvertices ; vno++) {
    pv = &pmris->vertices[vno] ;
    if (pv->ripflag)
      continue ;
 
    if(Gb_lowPassFilter) {
	if(fabs(pv->curv)<fabs(Gf_lowPassFilter))
	    b_canCount = 1;
	else
	    b_canCount = 0;
    }
    if(Gb_highPassFilter) {
	if(fabs(pv->curv)>=fabs(Gf_highPassFilter))
	    b_canCount = 1;
	else
	    b_canCount = 0;
    }
    if(Gb_lowPassFilterGaussian) {
	if(fabs(pv->K)<fabs(Gf_lowPassFilterGaussian))
	    b_canCount = 1;
	else
	    b_canCount = 0;
    }
    if(Gb_highPassFilterGaussian) {
	if(fabs(pv->K)>=fabs(Gf_highPassFilterGaussian))
	    b_canCount = 1;
	else
	    b_canCount = 0;
    }
    if(b_canCount && (*fcond)(pv)) {
        f_total += ((*fv)(pv) * pv->area) ;
        f_n 	+= 1.0 ;
	(*p_verticesCounted)++;
    }
  }
  if(f_n > 1) 	f_ret = f_total/f_n;
  return(f_ret);
}

int
MRISzeroCurvature(
  MRI_SURFACE* apmris
) {
  //
  // POSTCONDITIONS
  // o Each curvature (as well as Gaussian and Mean )value in
  //   apmris is simply set to zero.
  //

  VERTEX* pvertex;
  int  vno;

  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    pvertex = &apmris->vertices[vno] ;
    if (pvertex->ripflag)
      continue;
    pvertex->curv  = 0.;
    pvertex->K = 0.;
    pvertex->H = 0.;
  }
  return(NO_ERROR);
}

int
MRISscaleCurvature(
  MRI_SURFACE* apmris,
  float   af_scale) {
  //
  // POSTCONDITIONS
  // o Each curvature value in apmris is simply scaled by <af_scale>
  //

  VERTEX* pvertex;
  int  vno;
  int  vtotal;
  double f_mean;

  for (f_mean = 0.0, vtotal = vno = 0 ; vno < apmris->nvertices ; vno++) {
    pvertex = &apmris->vertices[vno] ;
    if (pvertex->ripflag)
      continue ;
    vtotal++ ;
    f_mean += pvertex->curv ;
  }
  f_mean /= (double)vtotal ;

  for (vno = 0 ; vno < apmris->nvertices ; vno++) {
    pvertex = &apmris->vertices[vno] ;
    if (pvertex->ripflag)
      continue;
    pvertex->curv = (pvertex->curv - f_mean) * af_scale + f_mean ;
  }
  return(NO_ERROR);
}

void
OFSP_create(
  char*  apch_prefix,
  char*  apch_suffix,
  char*  apch_curv,
  e_OFSP  ae_OFSP
) {
  //
  // PRECONDITIONS
  // o pch_prefix / pch_suffix must be large enough to contain
  //   required text info.
  //
  // POSTCONDITIONS
  //  o a full pch_prefix is composed of:
  //  <output_fname>.<hemi>.<surface>.<pch_suffix>
  //  o a partial pch_prefix is compsed of:
  //  <output_fname>.<hemi>.<surface>.<pch_suffix>
  //

  switch (ae_OFSP) {
  case e_None:
    sprintf(apch_prefix, output_fname);
    break;
  case e_Partial:
    sprintf(apch_prefix, "%s.%s.%s.%s",
            hemi,
            surf_name,
            apch_curv,
            apch_suffix
           );
    break;
  case e_Full:
    sprintf(apch_prefix, "%s.%s.%s.%s.%s",
            output_fname,
            hemi,
            surf_name,
            apch_curv,
            apch_suffix
           );
    break;
  }
}

void
histogram_wrapper(
  MRIS*   amris,
  e_secondOrderType aesot
) {
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

  int  		i;
  float* 	pf_histogram;
  float 	f_maxCurv 	= amris->max_curv;
  float 	f_minCurv 	= amris->min_curv;
  double  	f_binSize 	= 0.;
  int  		b_error  	= 0;
  short		b_writeToFile	= 0;
  FILE*		pFILE		= NULL;

  if (Gb_binSizeOverride)
    f_binSize = Gf_binSize;

  if (Gb_histStartOverride)
    f_minCurv = Gf_histStart;

  if (f_minCurv < amris->min_curv)
    f_minCurv = amris->min_curv;

  if (Gb_histEndOverride)
    f_maxCurv = Gf_histEnd;

  if (f_maxCurv > amris->max_curv)
    f_maxCurv = amris->max_curv;

  f_binSize  = ((double)f_maxCurv - (double)f_minCurv) / (double)G_bins;
  if (f_maxCurv <= f_minCurv)
    ErrorExit(ERROR_SIZE, "%s: f_maxCurv < f_minCurv.",
              Progname);

  b_error =   (long)(((double)f_minCurv+(double)G_bins*(double)f_binSize)*1e10) >
              (long)(((double)f_maxCurv)*1e10);

  if (b_error)
    ErrorExit(ERROR_SIZE, "%s: Invalid <binSize> and <bins> combination",
              Progname);

  pf_histogram = calloc(G_bins, sizeof(float));
  histogram_create( amris,
                    f_minCurv,
                    f_binSize,
                    G_bins,
                    pf_histogram,
                    aesot);

  printf("%*s%*s%*s\n", G_leftCols, "bin start",
         G_leftCols, "bin end",
         G_leftCols, "count");
  for (i=0; i<G_bins; i++) {
    printf("%*.4f%*.4f%*.4f ",
           G_leftCols, (i*f_binSize)+f_minCurv,
           G_leftCols, ((i+1)*f_binSize)+f_minCurv,
           G_leftCols, pf_histogram[i]);
    printf("\n");
    switch (aesot) {
    case e_Raw:
      if (GpFILE_rawHist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_rawHist; }
      break;
    case e_Gaussian:
      if (GpFILE_KHist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_KHist; }
      break;
    case e_Mean:
      if (GpFILE_HHist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_HHist; }
      break;
    case e_K1:
      if (GpFILE_K1Hist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_K1Hist; }
      break;
    case e_K2:
      if (GpFILE_K2Hist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_K2Hist; }
      break;
    case e_S:
      if (GpFILE_SHist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_SHist; }
      break;
    case e_C:
      if (GpFILE_CHist !=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_CHist; }
      break;
    case e_BE:
      if (GpFILE_BEHist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_BEHist; }
      break;
    case e_SI:
      if (GpFILE_SIHist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_SIHist; }
      break;
    case e_Normal:
      if (GpFILE_normHist!=NULL)
	{ b_writeToFile = 1; pFILE = GpFILE_normHist; }
      break;
    case e_ScaledTrans:
    case e_Scaled:
      if (GpFILE_scaledHist!=NULL) 
	{ b_writeToFile = 1; pFILE = GpFILE_scaledHist; }
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
#define true  1
#define false 0
int AlmostEqual2sComplement(float A, float B, int maxUlps) {
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
  if (aInt < 0)
    aInt = 0x80000000 - aInt;
  // Make bInt lexicographically ordered as a twos-complement int
  bInt = *(int*)ptrB;
  if (bInt < 0)
    bInt = 0x80000000 - bInt;
  intDiff = abs(aInt - bInt);
  if (intDiff <= maxUlps)
    return true;
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
) {
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

  int  vno  = 0;
  float* pf_curvature = NULL;
  int  i  = 0;
  int  j  = 0;
  int  start  = 0;
  int  count  = 0;
  int  totalCount = 0;
  int  nvertices = 0;

  int  b_onLeftBound = 0; // Don't trigger higher order float comp
  int  b_onRightBound = 0; // on left/right bounds

  double l_curvature; // These three variables
  double l_leftBound; // are scaled up and truncated
  double l_rightBound; // to minimise rounding errors

  char pch_sot[64];
  strcpy(pch_sot, Gppch[aesot]);

  nvertices  = amris_curvature->nvertices;
  fprintf(stdout, "\n%*s%s = %f\n",
          G_leftCols, pch_sot, " bin size", af_binSize);
  fprintf(stdout, "%*s%s = %d\n",
          G_leftCols, pch_sot, " surface vertices", nvertices);
  pf_curvature = calloc(nvertices, sizeof(float));

  for (vno=0; vno < nvertices; vno++)
    pf_curvature[vno] = amris_curvature->vertices[vno].curv;

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
    if (Gb_histogramPercent)
      apf_histogram[i] /= nvertices;
    if (totalCount == nvertices)
      break;
  }
  fprintf(stdout, "%*s%s = %d\n",
          G_leftCols, pch_sot, " sorted vertices", totalCount);
  fprintf(stdout, "%*s%s = %f\n",
          G_leftCols, pch_sot, " ratio", (float)totalCount / (float)nvertices);
  free(pf_curvature);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
// static int
// longOpts_process(int argc, char** ppch_argv) {
//     int 	c;
//     int 	option_index	= 0;
// 
//     opterr = 0;
//     c = getopt_long (argc, ppch_argv, "sl:L:h:H:",
//      			Gs_longOptions, &option_index);
// 
//     switch(c) {
//     case 'l':
// 	Gb_lowPassFilter		= 1;
// 	Gf_lowPassFilter		= atof(optarg);
// 	break;
//     case 'L':
// 	Gb_lowPassFilterGaussian	= 1;
// 	Gf_lowPassFilterGaussian	= atof(optarg);
// 	break;
//     case 'h':
// 	Gb_highPassFilter		= 1;
// 	Gf_highPassFilter		= atof(optarg);
// 	break;
//     case 'H':
// 	Gb_highPassFilterGaussian	= 1;
// 	Gf_highPassFilterGaussian	= atof(optarg);
// 	break;
//     }
// 
//     return 1;
// }

static int
get_option(int argc, char *argv[]) {
  int    nargs   = 0;
  char*  option;
  char*  pch_text;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-lowPassFilter")) {
    Gb_lowPassFilter		= 1;
    Gf_lowPassFilter		= atof(argv[2]);
    nargs			= 1;
    fprintf(stderr, "Setting rectified low pass filter to %f\n", 
			Gf_lowPassFilter);
  } else if (!stricmp(option, "-lowPassFilterGaussian")) {
    Gb_lowPassFilterGaussian	= 1;
    Gf_lowPassFilterGaussian	= atof(argv[2]);
    nargs			= 1;
    fprintf(stderr, "Setting rectified low pass Gaussian filter to %f\n", 
			Gf_lowPassFilterGaussian);
  } else if (!stricmp(option, "-highPassFilterGaussian")) {
    Gb_highPassFilterGaussian	= 1;
    Gf_highPassFilterGaussian	= atof(argv[2]);
    nargs			= 1;
    fprintf(stderr, "Setting rectified high pass Gaussian filter to %f\n", 
			Gf_highPassFilterGaussian);
  } else if (!stricmp(option, "-highPassFilter")) {
    Gb_highPassFilter		= 1;
    Gf_highPassFilter		= atof(argv[2]);
    nargs			= 1;
    fprintf(stderr, "Setting rectified high pass filter to %f\n", 
			Gf_highPassFilter);
  } else if (!stricmp(option, "-writeCurvatureFiles")) {
    Gb_writeCurvatureFiles	= 1;
    fprintf(stderr, "Toggling save flag on curvature files\n");
  } else switch (toupper(*option)) {
  case 'O':
    output_fname  = argv[2] ;
    Gb_output2File = 1;
    nargs   = 1 ;
    fprintf(stderr, "Outputting results using filestem '%s'...\n", output_fname) ;
    break ;
  case 'F':
    pch_text = argv[2];
    strcpy(surf_name, pch_text);
    nargs = 1 ;
    fprintf(stderr, "Setting surface to '%s'...\n", surf_name);
    break;
  case 'L':
    label_name = argv[2] ;
    fprintf(stderr, "Using label '%s'...\n", label_name) ;
    nargs = 1 ;
    break ;
  case 'A':
    navgs = atoi(argv[2]) ;
    fprintf(stderr, "Averaging curvature %d times...\n", navgs);
    nargs = 1 ;
    break ;
  case 'C':
    Gb_scale  = 1;
    Gf_scaleFactor  = atof(argv[2]) ;
    fprintf(stderr, "Setting raw scale factor to %f...\n", Gf_scaleFactor);
    nargs   = 1 ;
    break ;
  case 'D':
    Gb_scaleMin  = 1;
    Gf_scaleMin  = atof(argv[2]);
    fprintf(stderr, "Setting scale min factor to %f...\n", Gf_scaleMin);
    nargs  = 1;
    break;
  case 'E':
    Gb_scaleMax  = 1;
    Gf_scaleMax  = atof(argv[2]);
    fprintf(stderr, "Setting scale max factor to %f...\n", Gf_scaleMax);
    nargs  = 1;
    break;
  case 'G':
    Gb_gaussianAndMean = 1;
    break ;
  case 'S':   /* write out stats */
    stat_flag   = 1 ;
    condition_no  = atoi(argv[2]) ;
    nargs   = 1 ;
    fprintf(stderr, "Setting out summary statistics as condition %d...\n",
            condition_no) ;
    break ;
  case 'H':   /* histogram */
    if (argc == 2)
      print_usage();
    Gb_histogram = 1;
    Gb_histogramPercent = 0;
    Gb_minMaxShow = 1;
    G_bins   = atoi(argv[2]);
    nargs   = 1 ;
    fprintf(stderr, "Setting curvature histogram to %d bins...\n",
            G_bins);
    if (G_bins <=0 )
      ErrorExit(ERROR_BADPARM, "%s: Invalid bin number.\n", Progname);
    break ;
  case 'P':   /* percentage histogram */
    Gb_histogram = 1;
    Gb_histogramPercent = 1;
    Gb_minMaxShow = 1;
    G_bins   = atoi(argv[2]);
    nargs   = 1 ;
    fprintf(stderr, "Creating percentage curvature histogram with %d bins...\n",
            G_bins);
    break ;
  case 'B':   /* binSize override*/
    Gb_binSizeOverride = 1;
    Gf_binSize   = (double) atof(argv[2]);
    nargs   = 1 ;
    fprintf(stderr, "Setting curvature histogram binSize to %f...\n",
            Gf_binSize);
    break;
  case 'I':   /* histogram start */
    Gb_histStartOverride = 1;
    Gf_histStart   = atof(argv[2]);
    nargs    = 1 ;
    fprintf(stderr, "Setting histogram start point to %f...\n",
            Gf_histStart);
    break;
  case 'J':   /* histogram end */
    Gb_histEndOverride = 1;
    Gf_histEnd   = atof(argv[2]);
    nargs   = 1 ;
    fprintf(stderr, "Setting histogram end point to %f...\n",
            Gf_histEnd);
    break;
  case 'Z':   /* zero a target vertex */
    Gb_zeroVertex = 1;
    G_zeroVertex = atoi(argv[2]);
    nargs   = 1 ;
    fprintf(stderr, "Setting zero vertex index to %d...\n",
            G_zeroVertex);
    break;
  case 'Q':   /* set maxUlps */
    Gb_maxUlps  = 1;
    G_maxUlps  = atoi(argv[2]);
    nargs   = 1 ;
    fprintf(stderr, "Setting maxUlps to %d...\n",
            G_maxUlps);
    break;
  case 'N':
    normalize_flag  = 1 ;
    fprintf(stderr, "Setting normalisation ON...\n");
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
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }
  return(nargs) ;
}

void
outputFileNames_create(
  char*  apch_curv
) {
  //
  // POSTCONDITIONS
  // o All necessary (depending on user flags) file names are created.
  //
  OFSP_create(Gpch_log,  	Gpch_logS, 		apch_curv, e_None);
  OFSP_create(Gpch_allLog, 	Gpch_allLogS, 		apch_curv, e_Full);
  OFSP_create(Gpch_rawHist, 	Gpch_rawHistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_normHist, 	Gpch_normHistS, 	apch_curv, e_Full);
  OFSP_create(Gpch_normCurv, 	Gpch_normCurvS, 	apch_curv, e_Partial);
  OFSP_create(Gpch_KHist,  	Gpch_KHistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_KCurv,  	Gpch_KCurvS, 		apch_curv, e_Partial);
  OFSP_create(Gpch_HHist, 	Gpch_HHistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_HCurv, 	Gpch_HCurvS, 		apch_curv, e_Partial);
  OFSP_create(Gpch_K1Hist,  	Gpch_K1HistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_K1Curv,  	Gpch_K1CurvS, 		apch_curv, e_Partial);
  OFSP_create(Gpch_K2Hist,  	Gpch_K2HistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_K2Curv,  	Gpch_K2CurvS, 		apch_curv, e_Partial);
  OFSP_create(Gpch_SHist,  	Gpch_SHistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_SCurv,  	Gpch_SCurvS, 		apch_curv, e_Partial);
  OFSP_create(Gpch_CHist,  	Gpch_CHistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_CCurv,  	Gpch_CCurvS, 		apch_curv, e_Partial);
  OFSP_create(Gpch_BEHist,  	Gpch_BEHistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_BECurv,  	Gpch_BECurvS, 		apch_curv, e_Partial);
  OFSP_create(Gpch_SIHist,  	Gpch_SIHistS, 		apch_curv, e_Full);
  OFSP_create(Gpch_SICurv,  	Gpch_SICurvS, 		apch_curv, e_Partial);
  OFSP_create(Gpch_scaledHist,	Gpch_scaledHistS,	apch_curv, e_Full);
  OFSP_create(Gpch_scaledCurv,	Gpch_scaledCurvS,	apch_curv, e_Partial);
}

void
outputFiles_open(void) {
  //
  // POSTCONDITIONS
  // o All necessary (depending on user flags) files are opened.
  //

  if (Gb_output2File) {
    CE("\n\tFiles processed for this curvature:\n");
    printf("%*s\n", G_leftCols*2, Gpch_log);
    if ((GpFILE_log=fopen(Gpch_log, "a"))==NULL)
      ErrorExit(ERROR_NOFILE, "%s: Could not open file '%s' for apending.\n",
                Progname, Gpch_log);
    printf("%*s\n", G_leftCols*2, Gpch_allLog);
    if ((GpFILE_allLog=fopen(Gpch_allLog, "w"))==NULL)
      ErrorExit(ERROR_NOFILE, "%s: Could not open file '%s' for writing.\n",
                Progname, Gpch_allLog);
    if (Gb_histogram) {
      printf("%*s\n", G_leftCols*2, Gpch_rawHist);
      if ((GpFILE_rawHist=fopen(Gpch_rawHist, "w"))==NULL)
        ErrorExit(ERROR_NOFILE, "%s: Could not open file '%s' for writing.\n",
                  Progname, Gpch_rawHist);
    }
    if (normalize_flag) {
      if (Gb_histogram) {
        printf("%*s\n", G_leftCols*2, Gpch_normHist);
        if ((GpFILE_normHist=fopen(Gpch_normHist, "w"))==NULL)
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_normHist);
      }
      printf("%*s\n", G_leftCols*2, Gpch_normCurv);
    }
    if (Gb_gaussianAndMean) {
      if (Gb_histogram) {
        if ((GpFILE_KHist=fopen(Gpch_KHist, "w"))==NULL) {
          printf("%*s\n", G_leftCols*2, Gpch_KHist);
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_KHist);
        }
        if ((GpFILE_HHist=fopen(Gpch_HHist, "w"))==NULL) {
          printf("%*s\n", G_leftCols*2,Gpch_HHist);
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_HHist);
        }
        if ((GpFILE_K1Hist=fopen(Gpch_K1Hist, "w"))==NULL) {
          printf("%*s\n", G_leftCols*2, Gpch_K1Hist);
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_K1Hist);
        }
        if ((GpFILE_K2Hist=fopen(Gpch_K2Hist, "w"))==NULL) {
          printf("%*s\n", G_leftCols*2, Gpch_K2Hist);
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_K2Hist);
        }
        if ((GpFILE_SHist=fopen(Gpch_SHist, "w"))==NULL) {
          printf("%*s\n", G_leftCols*2, Gpch_SHist);
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_SHist);
        }
        if ((GpFILE_CHist=fopen(Gpch_CHist, "w"))==NULL) {
          printf("%*s\n", G_leftCols*2, Gpch_CHist);
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_CHist);
        }
        if ((GpFILE_BEHist=fopen(Gpch_BEHist, "w"))==NULL) {
          printf("%*s\n", G_leftCols*2, Gpch_BEHist);
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_BEHist);
        }
        if ((GpFILE_SIHist=fopen(Gpch_SIHist, "w"))==NULL) {
          printf("%*s\n", G_leftCols*2, Gpch_SIHist);
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_SIHist);
        }
      }
      printf("%*s\n", G_leftCols*2, Gpch_KCurv);
      printf("%*s\n", G_leftCols*2, Gpch_HCurv);
    }
    if (Gb_scale || (Gb_scaleMin && Gb_scaleMax)) {
      if (Gb_histogram) {
        printf("%*s\n", G_leftCols*2, Gpch_scaledHist);
        if ((GpFILE_scaledHist=fopen(Gpch_scaledHist, "w"))==NULL)
          ErrorExit(ERROR_NOFILE,
                    "%s: Could not open file '%s' for writing.\n",
                    Progname, Gpch_scaledHist);
      }
      printf("%*s\n", G_leftCols*2, Gpch_scaledCurv);
    }
    printf("\n");
  }
}

void
outputFiles_close(void) {
  //
  // POSTCONDITIONS
  // o Checks for any open files and closes them.
  //

  if (GpFILE_log) 	 fclose(GpFILE_log);	   GpFILE_log  		= NULL;
  if (GpFILE_allLog) 	 fclose(GpFILE_allLog);    GpFILE_allLog  	= NULL;
  if (GpFILE_rawHist) 	 fclose(GpFILE_rawHist);   GpFILE_rawHist 	= NULL;
  if (GpFILE_normHist) 	 fclose(GpFILE_normHist);  GpFILE_normHist 	= NULL;
  if (GpFILE_KHist) 	 fclose(GpFILE_KHist);	   GpFILE_KHist 	= NULL;
  if (GpFILE_HHist) 	 fclose(GpFILE_HHist);	   GpFILE_HHist 	= NULL;
  if (GpFILE_K1Hist) 	 fclose(GpFILE_K1Hist);    GpFILE_K1Hist 	= NULL;
  if (GpFILE_K2Hist) 	 fclose(GpFILE_K2Hist);    GpFILE_K2Hist 	= NULL;
  if (GpFILE_SHist) 	 fclose(GpFILE_SHist);	   GpFILE_SHist 	= NULL;
  if (GpFILE_CHist) 	 fclose(GpFILE_CHist);	   GpFILE_CHist 	= NULL;
  if (GpFILE_BEHist) 	 fclose(GpFILE_BEHist);	   GpFILE_BEHist 	= NULL;
  if (GpFILE_SIHist) 	 fclose(GpFILE_SIHist);	   GpFILE_SIHist 	= NULL;
  if (GpFILE_scaledHist) fclose(GpFILE_scaledHist);GpFILE_scaledHist 	= NULL;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  print_help() ;
}

static void
print_help(void) {
  char  pch_synopsis[65536];

  sprintf(pch_synopsis,"\n\
          \n\
          NAME         \n\
          \n\
          mris_curvature_stats      \n\
          \n\
          SYNOPSIS        \n\
          \n\
          mris_curvature_stats [options] \\    \n\
          <subjectName> <hemi> <curvFile1> [... <curvFileN] \n\
          \n\
          DESCRIPTION        \n\
          \n\
          'mris_curvature_stats' will primarily compute the mean and \n\
          variances for a curvature file (or a set of curvature files). \n\
          These files can be possibly constrained by a label.  \n\
          \n\
          Additionally, 'mris_curvature_stats' can report the max/min \n\
          curvature values, and compute a simple histogram based on  \n\
          these vales.       \n\
          \n\
          Curvatures can also be normalised and constrained to a given \n\
          range before computation.     \n\
          \n\
          Gaussian and mean calculations on a surface structure can \n\
          be performed.       \n\
          \n\
          Finally, all output to the console, as well as any new  \n\
          curvatures that result from the above calculations can be \n\
          saved to a series of text and binary-curvature files.  \n\
          \n\
          Please see the following options for more details.  \n\
          \n\
          OPTIONS         \n\
          \n\
          [-a <numberOfAverages]      \n\
          \n\
          Average the curvature <numberOfAverages> times.   \n\
          \n\
          [-o <outputFileStem>]      \n\
          \n\
          Save processing results to a series of files. This includes \n\
          condensed text output, histogram files (MatLAB friendly) \n\
          and curvature files.      \n\
          \n\
          The actual files that are saved depends on which additional \n\
          calculation flags have been specified (i.e. normalisation, \n\
          Gaussian / Mean, scaling).     \n\
          \n\
          In the case when a Gaussian/Mean calculation has been   \n\
          performed, 'mris_curvature_stats' will act in a manner  \n\
          similar to 'mris_curvature -w'.  Note though that though the \n\
          name of the curvature files that are created are different,  \n\
          the contents are identical to 'mris_curvature -w' created files.\n\
          \n\
          All possible files that can be saved are:   \n\
          (where  	OHSC  = <outputFileStem>.<hemi>.<surface>.<curvFile>\n\
          		HSC   = <hemi>.<surface>.<curvFile>   \n\
          \n\
          <outputFileStem> 		(Log only a single  \n\
          				mean+-sigma. If several \n\
          				curvature files have been \n\
          				specified, log the mean+-sigma \n\
          				across all the curvatures. \n\
					Note also that this file is \n\
					*appended* for each new run.) \n\
          <OHSC>.log  			(Full output, i.e the output of \n\
          				each curvature file mean +- \n\
          				sigma, as well as min/max \n\
          				as it is processed.)    \n\
          <OHS>.raw.hist  		(Raw histogram file. By 'raw'   \n\
          				is implied that the curvature  \n\
          				has not been further processed \n\
          				in any manner.)  \n\
          <OHS>.norm.hist  		(Normalised histogram file) \n\
          <HS>.norm.crv   		(Normalised curv file)  \n\
          <OHS>.K.hist  		(Gaussian histogram file) \n\
          <HS>.K.crv   			(Gaussian curv file)  \n\
          <OHS>.H.hist  		(Mean histogram file)  \n\
          <HS>.H.crv   			(Mean curv file)  \n\
          <OHS>.scaled.hist 		(Scaled histogram file)  \n\
          <HS>.scaled.crv  		(Scaled curv file)  \n\
          \n\
          Note that curvature files are saved to $SUBJECTS_DIR/surf \n\
          and *not* to the current working directory.   \n\
          \n\
          Note also that file names can become quite long and somewhat \n\
          unyielding.       \n\
          \n\
          [-h <numberOfBins>] [-p <numberOfBins]    \n\
          \n\
          If specified, prepare a histogram over the range of curvature \n\
          values, using <numberOfBins> buckets. These are dumped to \n\
          stdout.        \n\
          \n\
          If '-p' is used, then the histogram is expressed as a   \n\
          percentage.       \n\
          \n\
          Note that this histogram, working off float values and float \n\
          boundaries, can suffer from rounding errors! There might be \n\
          instances when a very few (on average) curvature values might \n\
          not be sorted.       \n\
          \n\
          The histogram behaviour can be further tuned with the   \n\
          following:       \n\
          \n\
          [-b <binSize>] [-i <binStartCurvature] [-j <binEndCurvature] \n\
          \n\
          These arguments are only processed iff a '-h <numberOfBins>' \n\
          has also been specified. By default, <binSize> is defined as \n\
          \n\
          (maxCurvature - minCurvature) / <numberOfBins>  \n\
          \n\
          The '-b' option allows the user to specify an arbitrary  \n\
          <binSize>. This is most useful when used in conjunction with \n\
          the '-i <binStartCurvature>' option, which starts the histogram \n\
          not at (minCurvature), but at <binStartCurvature>. So, if  \n\
          a histogram reveals that most values seem confined to a very \n\
          narrow range, the '-b' and '-i' allow the user to 'zoom in' \n\
          to this range and expand.     \n\
          \n\
          If <binStartCurvature> < (minCurvature), then regardless \n\
          of its current value, <binStartCurvature> = (minCurvature). \n\
          Also, if (<binStartCurvature> + <binSize>*<numberOfBins> >) \n\
          (maxCurvature), an error is raised and processing aborts. \n\
          \n\
          The '-j' allows the user to specify an optional end  \n\
          value for the histogram. Using '-i' and '-j' together  \n\
          are the most convenient ways to zoom into a region of interest \n\
          in a histogram.       \n\
          \n\
          [-l <labelFileName>]      \n\
          \n\
          Constrain statistics to the region defined in <labelFileName>. \n\
          \n\
          [-m]        \n\
          \n\
          Output min / max information for the processed curvature. \n\
          \n\
          [-n]        \n\
          \n\
          Normalise the curvature before computation. Normalisation \n\
          takes precedence over scaling, so if '-n' is specified  \n\
          in conjunction with '-c' or '-smin'/'-smax' it will  \n\
          override the effects of the scaling.    \n\
          \n\
          If specified in conjunction with '-o <outputFileStem>'  \n\
          will also create a curvature file with these values.  \n\
          \n\
          [-s <summaryCondition>]      \n\
          \n\
          Write out stats as <summaryCondition>.    \n\
          \n\
          [-d <minCurvature> -e <maxCurvature>]    \n\
          \n\
          Scale curvature values between <minCurvature> and   \n\
          <maxCurvature>. If the minimum curvature is greater  \n\
          than the maximum curvature, or if either is    \n\
          unspecified, these flags are ignored.    \n\
          \n\
          This scale computation takes precedence over '-c' scaling. \n\
          \n\
          Note also that the final scaling bounds might not correspond \n\
          to <minCurvature>... <maxCurvature> since values are scaled \n\
          across this range so as to preserve the original mean profile. \n\
          \n\
          If specified in conjunction with '-o <outputFileStem>'  \n\
          will also create a curvature file with these values.  \n\
          \n\
          [-c <factor>]       \n\
          \n\
          Scale curvature values with <factor>. The mean of the   \n\
          original curvature is conserved (and the sigma increases \n\
          with <factor>).       \n\
          \n\
          If specified in conjunction with '-o <outputFileStem>'  \n\
          will also create a curvature file with these values.  \n\
          \n\
          [-version]        \n\
          \n\
          Print out version number.     \n\
          \n\
          [-z <vertexIndex>]       \n\
          \n\
          Sets the curvature values at that index to zero. The   \n\
          'raw' curvature, as well as the Gaussian and Mean curvatures \n\
          are set to zero, and min/max values are recomputed.  \n\
          \n\
          This is useful in cases when outliers in the data (particularly \n\
          evident in Gaussian calcuations) skew mean and sigma values. \n\
          \n\
          [-q <maxUlps>]       \n\
          \n\
          The <maxUlps> is used to toggle a more rigorous floating point \n\
          comparison operation in the histogram function. Comparing  \n\
          float values for sorting into bins can at times fail due to \n\
          number precision issues. If, over the range of comparison \n\
          some curvature values are not sorted, add <maxUlps>.  \n\
          \n\
          This adds extra function calls to AlmostEqual2sComplement(..) \n\
          for float comparisons and improves the general accuracy, at  \n\
          a very slight performance penalty.    \n\
          \n\
          You will most likely never have to use this argument, and is  \n\
          for advanced use only.      \n\
          \n\
          NOTES         \n\
          \n\
          It is important to note that some combinations of the command \n\
          line parameters are somewhat meaningless, such as normalising \n\
          a 'sulc' curvature file (since it's normalised by definition). \n\
          \n\
          EXAMPLES        \n\
          \n\
          mris_curvature_stats 801_recon rh curv    \n\
          \n\
          For subject '801_recon', determine the mean and sigma for \n\
          the curvature file on the right hemisphere.   \n\
          \n\
          mris_curvature_stats -m 801_recon rh curv    \n\
          \n\
          Same as above, but print the min/max curvature values  \n\
          across the surface.      \n\
          \n\
          mris_curvature_stats -h 20 -m 801_recon rh curv   \n\
          \n\
          Same as above, and also print a histogram of curvature   \n\
          values over the min/max range, using 20 bins. By replacing \n\
          the '-h' with '-p', print the histogram as a percentage. \n\
          \n\
          mris_curvature_stats -h 20 -b 0.01 -i 0.1 -m 801_recon rh curv \n\
          \n\
          Same as above, but this time constrain the histogram to the 20  \n\
          bins from -0.1 to 0.1, with a bin size of 0.01.   \n\
          \n\
          Note that the count / percentage values are taken across the    \n\
          total curvature range and not the constrained window defined  \n\
          by the '-i' and '-b' arguments.     \n\
          \n\
          mris_curvature_stats -G -m 801_recon rh curv   \n\
          \n\
          Print the min/max curvatures for 'curv', and also calculate \n\
          the Gaussian and Mean curvatures (also printing the min/max     \n\
          for these).       \n\
          \n\
          mris_curvature_stats -G -F smoothwm -m 801_recon rh curv  \n\
          \n\
          By default, 'mris_curvature_stats' reads the 'orig' surface \n\
          for the passed subject. This is not generally the best surface \n\
          for Gaussian determination. The '-F' uses the 'smoothwm'  \n\
          surface, which is a better choice.    \n\
          \n\
          mris_curvature_stats -h 10 -G -F smoothwm -m 801_recon rh curv \n\
          \n\
          Same as above, with the addition of a histogram for the  \n\
          Gaussian and Mean curvatures as well.    \n\
          \n\
          mris_curvature_stats -h 10 -G -F smoothwm -m -o foo \\   \n\
          801_recon rh curv sulc      \n\
          \n\
          Generate several output text files that capture the min/max \n\
          and histograms for each curvature processed. Also create \n\
          new Gaussian and Mean curvature files.    \n\
          \n\
          In this case, the curvature files created are called:  \n\
          \n\
          rh.smoothwm.curv.K.crv     \n\
          rh.smoothwm.curv.H.crv     \n\
          rh.smoothwm.sulc.K.crv     \n\
          rh.smoothwm.sulc.H.crv     \n\
          \n\
          and are saved to the $SUBJECTS_DIR/<subjid>/surf directory. \n\
          These can be re-read by 'mris_curvature_stats' using   \n\
          \n\
          mris_curvature_stats -m 801_recon rh \\   \n\
          smoothwm.curv.K.crv smoothwm.sulc.K.crv  \n\
          \n\
          ADVANCED EXAMPLES       \n\
          \n\
          'mris_curvature_stats' can also provide some useful side  \n\
          effects. Reading in a curvature, and applying any calculation \n\
          to it (scaling, gaussian, etc.) can result in data that  \n\
          can be visualised quite well in a tool such as 'tksurfer'. \n\
          \n\
          Consider the normal curvature display in 'tksurfer', which \n\
          is usually quite dark due to the dynamic range of the display. \n\
          We can considerably improve the brightness by scaling a  \n\
          curvature file and rendering the resultant in 'tksurfer'. \n\
          \n\
          First, take an arbitrary curvature, apply a scale factor, \n\
          and an output filestem:      \n\
          \n\
          mris_curvature_stats -o foo -c 10 801_recon rh curv   \n\
          \n\
          This scales each curvature value by 10. A new curvature file \n\
          is saved in       \n\
          \n\
          $SUBJECTS_DIR/801_recon/surf/rh.orig.curv.scaled.crv \n\
          \n\
          Comparing the two curvatures in 'tksurfer' will clearly show \n\
          the scaled file as much brighter.    \n\
          \n\
          Similarly, the Gaussian curvature can be processed, scaled, and \n\
          displayed, yielding very useful visual information. First  \n\
          create and save a Gaussian curvature file (remember that the \n\
          smoothwm surface is a better choice than the default orig \n\
          surface):       \n\
          \n\
          mris_curvature_stats -o foo -G -F smoothwm 801_recon rh curv \n\
          \n\
          The 'foo' filestem is ignored when saving curvature files, but \n\
          needs to be specified in order to trigger output saving. This \n\
          command will create Gaussian and Mean curvature files in the \n\
          $SUBJECTS_DIR/surf directory:     \n\
          \n\
          rh.smoothwm.curv.K.crv     \n\
          rh.smoothwm.curv.H.crv     \n\
          \n\
          Now, process the created Gaussian with the scaled curvature: \n\
          \n\
          mris_curvature_stats -o foo -c 10 801_recon rh smoothwm.curv.K.crv \n\
          \n\
          Again, the 'foo' filestem is ignored, but needs to be specified \n\
          to trigger the save. The final scaled Gaussian curvature is  \n\
          saved to (again in the $SUBJECTS_DIR/801_recon/surf directory): \n\
          \n\
          rh.orig.smoothwm.curv.K.crv.scaled.crv   \n\
          \n\
          which is a much better candidate to view in 'tksurfer' than \n\
          the original Gaussian curvature file.    \n\
          \n\
          \n\
          \n");

  CE(pch_synopsis);
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

int  comp(const void *p, const void *q)
{
  float *ptr1 = (float *)(p);
  float *ptr2 = (float *)(q);

  if (*ptr1 < *ptr2) return -1;
  else if (*ptr1 == *ptr2) return 0;
  else return 1;
}
