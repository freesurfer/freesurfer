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
 *    $Date: 2007/09/05 14:58:27 $
 *    $Revision: 1.32 $
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
#include "xDebug.h"

#define  STRBUF  	65536
#define  MAX_FILES  	1000
#define  CO( x )  	fprintf(stdout, ( x ))
#define  CE( x )  	fprintf(stderr, ( x ))
#define  START_i   	3

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

typedef enum _surfaceIntegrals {
  e_natural		= 0,	// "Natural" (native) integral - no conditions
  e_rectified		= 1,	// curvature is first rectified
  e_pos			= 2,	// only positive curvatures are considered
  e_neg			= 3     // only negative curvatures are considered
} e_surfaceIntegral;

// Set of possible output curvature files
char* Gppch[] = {
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
                  "ScaledTrans"
                };

// Output file prefixes
typedef enum _OFSP {
    e_None,   			// No prefix
    e_Partial,  		// A partial prefix (does not include the stem)
    e_Full,   			// Full prefix - including stem
}   e_OFSP;

// This structure is used to house data on the min/max curvature of a
//	particular type (raw, K, H, etc.) across a surface. It is also
//	used to convey information between the various functions that
//	use it.
typedef struct _minMax {
    float	f_min;			// The minimum curvature value
    float	f_max;			// The maximum curvature value
    int		vertexMin;		// The vertex where min lives
    int		vertexMax;		// The vertex where max lives
    short	b_minTest;		// If true, perform explicit min test
    short	b_maxTest;		// If true, perform explicit max test
    short	b_minViolation;		// If true, min value != explicit min
    short	b_maxViolation;	  	// If true, max value != explicit max
} s_MINMAX;

static char vcid[] =
  "$Id: mris_curvature_stats.c,v 1.32 2007/09/05 14:58:27 rudolph Exp $";

int   main(int argc, char *argv[]) ;

static int   get_option(int argc, char *argv[]) ;
static void  usage_exit(void) ;
static void  print_usage(void) ;
static void  print_help(void) ;
static void  print_version(void) ;

	//----------------------------//
	// Function Prototypes: START //
	//----------------------------// 


// To mrisurf.h vvvvvvvvvvvvvvvvvv
short	FACES_aroundVertex_reorder(
	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			pv_geometricOrder
);

short	MRIScomputeSecondFundamentalFormDiscrete(
	MRIS*			apmris
);

short	FACE_vertexIndex_find(
    	FACE*			pFace,
    	int 			avertex 
);

short	VECTOR_elementIndex_findNotEqual(
	VECTOR*			apV,
	float			af_searchTerm
);

short	VECTOR_elementIndex_find(
	VECTOR*			apV,
	float			af_searchTerm
);

short	MRIS_vertexProgress_print(
    	MRIS*			apmris,
    	int			avertex,
    	char*			apch_message
);

int	FACE_vertexIndexAtMask_find(
	FACE*			apFACE_I,
	VECTOR*			apv_verticesCommon
);

short	VERTICES_commonInFaces_find(
	FACE*			apFACE_I,
	FACE*			apFACE_J,
	VECTOR*			apv_verticesCommon
);

short	FACES_Hcurvature_determineSign(
    	MRIS*			apmris,
    	FACE*			apFACE_O,
    	FACE*			apFACE_I
);

int	VERTEX_faceAngles_determine(
    	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			apv_angle
);

int	VERTEX_faceMinMaxAngles_determine(
    	MRIS*			apmris,
    	int			avertex,
    	int*			ap_minIndex,
    	float*			apf_minAngle,
    	int*			ap_maxIndex,
    	float*			apf_maxAngle
);

int	MRIS_facesAtVertices_reorder(
    	MRIS*			apmris
);

int	MRIScomputeGeometricProperties(
    	MRIS*			apmris
);

short	FACES_aroundVertex_reorder(
    	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			pv_geometricOrder
);

float	FACES_angleNormal_find(
    	MRIS*			apmris,
    	FACE*			apFACE_I,
    	FACE*			apFACE_J
);

float	FACES_commonEdgeLength_find(
    	MRIS*			apmris,
    	FACE*			apFACE_I,
    	FACE*			apFACE_J
);

short	MRIS_discreteKH_compute(
	MRIS*			apmris
);

short	MRIS_discretek1k2_compute(
	MRIS*			apmris
);

// The discrete curvature calculations are based on the Gauss-Bonnet Scheme.
// 
// For 'K' and 'H', see
// 
// @article{1280456,
//  author = {Evgeni Magid and Octavian Soldea and Ehud Rivlin},
//  title = {A comparison of Gaussian and mean curvature estimation 
// 	  methods on triangular meshes of range image data},
//  journal = {Comput. Vis. Image Underst.},
//  volume = {107},
//  number = {3},
//  year = {2007},
//  issn = {1077-3142},
//  pages = {139--159},
//  doi = {http://dx.doi.org/10.1016/j.cviu.2006.09.007},
//  publisher = {Elsevier Science Inc.},
//  address = {New York, NY, USA},
//  }
// 
// and for k1 and k2 see:
// 
// @misc{ meyer02discrete,
//   author = "M. MEYER and M. DESBRUN and P. SCHR and A. BARR",
//   title = "Discrete DifferentialGeometry Operators for Triangulated 2-Manifolds",
//   text = "MEYER, M., DESBRUN, M., SCHR ODER, P., AND BARR, 
// 	  A. H. Discrete DifferentialGeometry
//     	  Operators for Triangulated 2-Manifolds, 2002. VisMath.",
//   year = "2002",
//   url = "citeseer.ist.psu.edu/meyer02discrete.html" }
// 

short	MRIScomputeSecondFundamentalFormDiscrete(
	MRIS*			apmris
);

int  	MRISminMaxCurvatureIndicesLookup(
  	MRI_SURFACE*  		apmris,
  	int*   			ap_vertexMin,
  	int*   			ap_vertexMax
);

int  	MRISvertexCurvature_set(
  	MRI_SURFACE*  		apmris,
  	int   			aindex,
  	float   		af_val
);

int	MRISzeroCurvature(
  	MRI_SURFACE*  		apmris
);

int	MRISuseK1Curvature(
  	MRI_SURFACE*  		mris
);

int	MRISuseK2Curvature(
  	MRI_SURFACE*  		mris
);


// To mrisurf.h ^^^^^^^^^^^^^^^^^^


// Simple functions of principle curvatures
float 	f_sharpnessCurvature(float af_k1, float af_k2) 
	{return ((af_k1 - af_k2)*(af_k1 - af_k2));}
float 	f_bendingEnergyCurvature(float af_k1, float af_k2) 
	{return (af_k1*af_k1 + af_k2*af_k2);}
float 	f_curvednessCurvature(float af_k1, float af_k2) 
	{return (sqrt(0.5*(af_k1*af_k1 + af_k2*af_k2)));}
float 	f_shapeIndexCurvature(float af_k1, float af_k2) 
	{return (af_k1 == af_k2 ? 0 : atan((af_k1+af_k2)/(af_k2 - af_k1)));}

// Simple functions on specific vertices
float 	f_absCurv(VERTEX*		pv) {return (fabs(pv->curv));}
float 	f_pass(VERTEX*			pv) {return (pv->curv);}
short 	f_allVertices(VERTEX*		pv) {return 1;}
short 	f_greaterEqualZero(VERTEX*	pv) {return (pv->curv >= 0 ? 1 : 0);}
short 	f_lessThanZero(VERTEX*		pv) {return (pv->curv <0 ? 1 : 0);}

int 	MRISusePrincipleCurvatureFunction(
	MRI_SURFACE*		pmris, 
	float 			(*f)(float k1, float k2)
);

short 	MRIS_surfaceIntegral_compute(
	MRI_SURFACE* 		pmris, 
	int* 			p_verticesCounted,
	float*			pf_areaCounted,
	short			(*fcond)(VERTEX*	pv),
	float			(*fv)	(VERTEX*	pv),
      	float*			pf_surfaceIntegral,
	float*			pf_meanIntegral,
	float*			pf_areaNormalizedIntegral	
);

int	comp(
	const void *p, 
	const void *q
	);    			// Compare p and q for qsort()

void	histogram_wrapper(
  	MRIS*   	amris,
  	e_secondOrderType aesot
);

void	histogram_create(
  	MRIS*   		amris_curvature,
  	float   		af_minCurv,
  	double   		af_binSize,
  	int   			abins,
  	float*   		apf_histogram,
  	e_secondOrderType 	aesot
);

void	OFSP_create(
  	char*   		apch_prefix,
  	char*   		apch_suffix,
  	e_OFSP   		ae_OFSP
);

void	outputFileNames_create(void);
void	outputFiles_open(void);
void	outputFiles_close(void);

int	MRIS_curvatureStats_analyze(
	MRIS*			apmris,
	e_secondOrderType	aesot
);

int  	MRISminMaxCurvaturesSearchSOT(
  	MRI_SURFACE*  		apmris,
  	int*   			ap_vertexMin,
  	int*   			ap_vertexMax,
  	float*   		apf_min,
  	float*   		apf_max,
  	e_secondOrderType 	aesot
);

int  	MRISminMaxCurvaturesSearch(
  	MRI_SURFACE*  		apmris,
  	int*   			ap_vertexMin,
  	int*   			ap_vertexMax,
  	float*   		apf_min,
  	float*   		apf_max
) {
  MRISminMaxCurvaturesSearchSOT(  apmris,
                                  ap_vertexMin,
                                  ap_vertexMax,
                                  apf_min,
                                  apf_max,
                                  e_Raw);
  return(NO_ERROR);
};

	//----------------------------//
	// Function Prototypes: END   //
	//----------------------------// 


// Global variables

char*  Progname ;
char*  hemi;

// Column / width for formatted output
static short	G_LC				= 50;
static short	G_RC				= 30;

static int  	navgs    			= 0;
static int  	normalize_flag   		= 0;
static int  	condition_no   			= 0;
static int  	stat_flag   			= 0;
static char* 	label_name   			= NULL;
static char* 	output_fname   			= NULL;
static char*	curv_fname			= NULL;
static char 	surf_name[STRBUF];
// Additional global variables (prefixed by 'G') added by RP.
// Flags are prefixed with Gb_ and are used to track
// user spec'd command line flags.
static int	Gb_shapeIndex			= 0;
static int	Gb_discreteCurvaturesUse	= 1;
static int  	Gb_minMaxShow  			= 0;
static int 	Gb_histogram  			= 0;
static int 	Gb_histogramPercent 		= 0;
static int 	Gb_binSizeOverride 		= 0;
static double 	Gf_binSize  			= 0.;
static int 	Gb_histStartOverride 		= 0;
static float 	Gf_histStart  			= 0.;
static int 	Gb_histEndOverride 		= 0;
static float 	Gf_histEnd  			= 0.;
static int 	Gb_gaussianAndMean 		= 0;
static int 	Gb_output2File  		= 0;
static int 	Gb_scale  			= 0;
static int 	Gb_scaleMin  			= 0;
static int 	Gb_scaleMax  			= 0;
static int 	Gb_zeroVertex  			= 0;
static int 	G_zeroVertex  			= 0;
static int  	G_nbrs    			= 2;
static int 	G_bins   			= 1;
static int 	Gb_maxUlps  			= 0;
static int 	G_maxUlps  			= 0;

static int	Gb_postScale			= 0;
static float	Gf_postScale			= 0.;
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

static char 	Gpch_K1Hist[STRBUF];
static char 	Gpch_K1HistS[] 			= "K1.hist";
static FILE* 	GpFILE_K1Hist  			= NULL;
static char 	Gpch_K1Curv[STRBUF];
static char 	Gpch_K1CurvS[] 			= "K1.crv";

static char 	Gpch_K2Hist[STRBUF];
static char 	Gpch_K2HistS[]  		= "K2.hist";
static FILE* 	GpFILE_K2Hist  			= NULL;
static char 	Gpch_K2Curv[STRBUF];
static char 	Gpch_K2CurvS[]  		= "K2.crv";

static char 	Gpch_SHist[STRBUF];
static char 	Gpch_SHistS[]  			= "S.hist";
static FILE* 	GpFILE_SHist  			= NULL;
static char 	Gpch_SCurv[STRBUF];
static char 	Gpch_SCurvS[]  			= "S.crv";

static char 	Gpch_CHist[STRBUF];
static char 	Gpch_CHistS[]  			= "C.hist";
static FILE* 	GpFILE_CHist  			= NULL;
static char 	Gpch_CCurv[STRBUF];
static char 	Gpch_CCurvS[]  			= "C.crv";

static char 	Gpch_BEHist[STRBUF];
static char 	Gpch_BEHistS[]  		= "BE.hist";
static FILE* 	GpFILE_BEHist  			= NULL;
static char 	Gpch_BECurv[STRBUF];
static char 	Gpch_BECurvS[]  		= "BE.crv";

static char 	Gpch_SIHist[STRBUF];
static char 	Gpch_SIHistS[]  		= "SI.hist";
static FILE* 	GpFILE_SIHist  			= NULL;
static char 	Gpch_SICurv[STRBUF];
static char 	Gpch_SICurvS[]  		= "SI.crv";

// These are used for tabular output
const int 	G_leftCols  			= 10;
const int 	G_rightCols 			= 40;

// Mean / sigma tracking and scaling
static double 	Gpf_means[MAX_FILES] ;
static double 	Gf_mean   			= 0.;
static double 	Gf_sigma  			= 0.;
static double 	Gf_n   				= 0.;
static double 	Gf_total  			= 0.;
static double 	Gf_total_sq  			= 0.;
static double 	Gf_scaleFactor 			= 1.;
static double 	Gf_scaleMin  			= 0.;
static double 	Gf_scaleMax  			= 0.;

static int which_norm = NORM_MEAN;

int 
stats_update(
	int	ai
) {
    //
    // PRECONDITIONS
    //	o Typically called if the mean and sigma for a surface
    //	  has been determined.
    //	o <ai> denotes the position in the global 'mean' array
    //	  to update.
    //
    // POSTCONDITIONS
    //	o Updates several globals with new mean/sigma data.
    //

    Gpf_means[ai-START_i]  	 = Gf_mean;
    Gf_total    		+= Gf_mean;
    Gf_total_sq   		+= Gf_mean*Gf_mean ;
    Gf_n   			+= 1.0;
    return 1;
}

void
cprints(
	char*		apch_left,
	char*		apch_right
) {
    //
    // PRECONDITIONS
    //	o The length of each text string should be such that
    //	  the string will "fit" into its column.
    //
    // POSTCONDITIONS
    //	o Column prints the left (action) and right (status)
    //	  text strings in a formatted manner.
    //

    static char	pch_right[16384];

    sprintf(pch_right, " [ %s ]\n", apch_right);
    if(strlen(apch_left))
        fprintf(stderr, "%*s", 	G_LC, 	apch_left);
    if(strlen(apch_right))
        fprintf(stderr, "%*s",	G_RC, 	pch_right);
    fflush(stderr);
}

void
cprintd(
	char*		apch_left,
	int		a_right
) {
    //
    // PRECONDITIONS
    //	o The length of each text string should be such that
    //	  the string will "fit" into its column.
    //
    // POSTCONDITIONS
    //	o Column prints the left (action) and right (status)
    //	  text strings in a formatted manner.
    //

    static char	pch_right[16384];

    sprintf(pch_right, " [ %d ]\n", a_right);
    if(strlen(apch_left))
        fprintf(stderr, "%*s", 	G_LC, 	apch_left);
    else
	fprintf(stderr, "%*s", 	G_LC,	" ");
    fprintf(stderr, "%*s",	G_RC, 	pch_right);
    fflush(stderr);
}

void
cprintf(
	char*		apch_left,
	float		af_right
) {
    //
    // PRECONDITIONS
    //	o The length of each text string should be such that
    //	  the string will "fit" into its column.
    //
    // POSTCONDITIONS
    //	o Column prints the left (action) and right (status)
    //	  text strings in a formatted manner.
    //

    static char	pch_right[16384];

    sprintf(pch_right, " [ %f ]\n", af_right);
    if(strlen(apch_left))
        fprintf(stderr, "%*s", 	G_LC, 	apch_left);
    else
	fprintf(stderr, "%*s", 	G_LC,	" ");
    fprintf(stderr, "%*s",	G_RC, 	pch_right);
    fflush(stderr);
}

int
main(int argc, char *argv[]) {
  char          **av, fname[STRBUF], *sdir ;
  char          *subject_name;
  char		pch_surface[16384];
  int           ac, nargs;
  int		i = START_i;
  MRI_SURFACE   *mris ;

  InitDebugging( "mris_curvature_stats" );
  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv,
	"$Id: mris_curvature_stats.c,v 1.32 2007/09/05 14:58:27 rudolph Exp $", "$Name:  $");
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
    strcpy(surf_name, "smoothwm") ;

  if (argc < 3)
    usage_exit() ;

  //
  // Once all the options have been processed, the command line is
  // mris_curvature_stats <str_subjectName> <str_hemi> [<str_curv0> ... <str_curvN>]
  //		0		1		2	    3	    ...	    N
  //
  // The <str_curvN> are not necessary if the Gaussian curvatures are calculated.

  subject_name 	= argv[1] ;
  hemi 		= argv[2] ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, surf_name) ;
  sprintf(pch_surface, "%s/%s.%s", subject_name, hemi, surf_name);
  cprints("Setting surface", pch_surface);
  cprints("Reading surface...", "");
  mris = MRISread(fname) ;
  cprints("", "ok");
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
    fprintf(stdout,"%s: ", label_name) ;

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
	MRIS_curvatureStats_analyze(mris, e_Raw); stats_update(i);
      }

  // Should we calculate all the principle curvature based curvatures? 
  // This is a surface-based calculation, and does not depend on the 
  // curvature processing loop - thus the 'curv' input is irrelevant
  if (Gb_gaussianAndMean) {
    MRISsetNeighborhoodSize(mris, G_nbrs);
    if(!Gb_discreteCurvaturesUse) {
	cprints("Calculating Continuous Principle Curvatures...", "");
        MRIScomputeSecondFundamentalForm(mris);
        cprints("", "ok");
    } else {
	cprints("Calculating Discrete Principle Curvatures...", "");
        MRIScomputeSecondFundamentalFormDiscrete(mris);
    }

    MRIS_curvatureStats_analyze(mris, e_Gaussian); stats_update(i++);
    if (Gpch_KCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_KCurv);
    MRIS_curvatureStats_analyze(mris, e_Mean); stats_update(i++);
    if (Gpch_HCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_HCurv);
    MRIS_curvatureStats_analyze(mris, e_K1); stats_update(i++);
    if (Gpch_K1Curv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_K1Curv);
    MRIS_curvatureStats_analyze(mris, e_K2); stats_update(i++);
    if (Gpch_K2Curv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_K2Curv);
    MRIS_curvatureStats_analyze(mris, e_S); stats_update(i++);
    if (Gpch_SCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_SCurv);
    MRIS_curvatureStats_analyze(mris, e_C); stats_update(i++);
    if (Gpch_CCurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_CCurv);
    MRIS_curvatureStats_analyze(mris, e_BE); stats_update(i++);
    if (Gpch_BECurv && Gb_writeCurvatureFiles) 
	MRISwriteCurvature(mris, 	Gpch_BECurv);
    
    // NOTE:
    // The "Shape Index" can be problematic due to the atan calculations.
    // To analyze the "Shape Index", you must pass an explicit
    // '--shapeIndex' on the command line
    if(Gb_shapeIndex) {
        MRIS_curvatureStats_analyze(mris, e_SI); stats_update(i++);
        if (Gpch_SICurv && Gb_writeCurvatureFiles) 
    	    MRISwriteCurvature(mris, 	Gpch_SICurv);
        }
  }

  if (Gf_n> 1.8) {
    Gf_mean 	= Gf_total / Gf_n;
    Gf_sigma 	= sqrt(Gf_total_sq/Gf_n- Gf_mean*Gf_mean) ;
    fprintf(stdout, "\nMean across %d curvatures: %8.4e +- %8.4e\n",
            (int) Gf_n, Gf_mean, Gf_sigma) ;
  }

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

// To mrisurf.c vvvvvvvvvvvvvvvvvv

short
FACE_vertexIndex_find(
    FACE*	pFace,
    int 	avertex 
) {
    //
    // PRECONDITIONS
    // o <avertex> denotes a vertex number to lookup in the <pFace>.
    //
    // POSTCONDITIONS
    // o The vertex index (0, 1, 2) containing the <avertex> is returned
    //	 or -1 if not found.
    //
    int 	vertex			= 0;
    int		ret			= -1;
    for(vertex=0; vertex<VERTICES_PER_FACE; vertex++)
	if(pFace->v[vertex] == avertex) ret = vertex;
    return ret;
}

short
VECTOR_elementIndex_findNotEqual(
	VECTOR*	apV,
	float	af_searchTerm
) {
    //
    // PRECONDITIONS
    //  o The <apV> is a column vector.
    //
    // POSTCONDITIONS
    // 	o The index of the first element in <apV> that is not equal to
    //	  the <af_searchTerm> is returned in the function name. If no
    //	  hits are found, a '0' is returned.
    //  o NOTE that the first element index in the vector is '1', not
    //	  zero (i.e. MatLAB convention).
    //
  	
    int		i	= 0;
    short	b_ret	= 0;
    for(i=0; i<apV->rows; i++)
	if(VECTOR_ELT(apV, i+1) != af_searchTerm) {
	    b_ret		= i;
	    break;
	}
    return b_ret;
}

short
VECTOR_elementIndex_find(
	VECTOR*	apV,
	float	af_searchTerm
) {
    //
    // PRECONDITIONS
    //  o The <apV> is a column vector.
    //
    // POSTCONDITIONS
    // 	o The index of the first element in <apV> that is equal to
    //	  the <af_searchTerm> is returned in the function name. If no
    //	  hits are found, a '0' is returned.
    //  o NOTE that the first element index in the vector is '1', not
    //	  zero (i.e. MatLAB convention).
    //
  	
    int		i	= 0;
    short	b_ret	= 0;
    for(i=1; i<=apV->rows; i++)
	if(VECTOR_ELT(apV, i) == af_searchTerm) {
	    b_ret		= i;
	    break;
	}
    return b_ret;
}

short
MRIS_vertexProgress_print(
    MRIS*	apmris,
    int		avertex,
    char*	apch_message
) {
    //
    // PRECONDITIONS
    //  o <avertex> is the current vertex being processed in a stream.
    //	o If <apch_message> is non-NULL, then prefix the progress bar
    //	  with <apch_message> (and terminate progress bar with [ ok ]).
    //
    // POSTCONDITIONS
    //	o For every 5% of processed vertices a "#" is written to stderr
    //

    static int 		totalVertices	= 0;
    static int		fivePerc	= 0;	
   
    totalVertices		= apmris->nvertices;
    fivePerc			= 0.05 * totalVertices;

    if(!avertex) {
        if(apch_message != NULL)
    	    fprintf(stderr, "%*s", G_LC, apch_message);
	fprintf(stderr, " [");
    }
    if(avertex%fivePerc == fivePerc-1)	fprintf(stderr, "#");
    if(avertex == apmris->nvertices-1) {
	fprintf(stderr, "] ");
        if(apch_message != NULL)
    	    fprintf(stderr, "%*s\n", 1, "[ ok ]");
    }
    return 1;
}

int
FACE_vertexIndexAtMask_find(
	FACE*	apFACE_I,
	VECTOR*	apv_verticesCommon
) {
    //
    // PRECONDITIONS
    //  o Called after <apv_verticesCommon> has been processed by
    //	  VERTICES_commonInFaces_find().
    //
    // POSTCONDITIONS
    //	o The vertex index in <apFACE_I> corresponding to the '-1'
    //	  in <apv_verticesCommon> is returned. For two bordering faces
    //	  that share an edge defined by <apv_verticesCommon>, this function
    //	  determines the index of the vertex that is *not* on this shared
    //	  edge.
    //
    // HISTORY
    // 03 July 2007
    //	o Initial design and coding.
    //

    int 	face		= 0;
    int 	vertex		= 0;
    short	b_inCommon	= 0;
    int		ret		= -1;

    for(face=0; face<3; face++) {
	vertex 		= apFACE_I->v[face];
	b_inCommon	= VECTOR_elementIndex_find(apv_verticesCommon, 
				(float) vertex);
	if(!b_inCommon) {
	    ret	= vertex;
	    break;
	}
    }

    return ret;
}

short
VERTICES_commonInFaces_find(
	FACE*	apFACE_I,
	FACE*	apFACE_J,
	VECTOR*	apv_verticesCommon
) {
    //
    // PRECONDITIONS
    //  o The <apFACE>s must be triangles with 3 vertices each.
    //	o It is assumed (but not mandatory) that the FACES share 
    //	  at least one common vertex - or more often a common 
    //	  edge, i.e. two common vertices.
    // 	o The <apv_uncommon> VECTOR's memory must be managed by the
    //	  caller, i.e. created and freed, and should be a 3x1 VECTOR.
    //
    // POSTCONDITIONS
    //	o The number of vertices that are in common between the two
    //	  faces are returned in the function name. This can be either 
    //	  (0, 1, 2, 3).
    //	o The indices of the common vertices are returned in 
    //	  apv_verticesCommon. This is a three element vector, with each
    //	  element corresponding to a common vertex. Vertex indices that
    //	  are not common between the two faces have a -1.
    //

    int		i		= 0;
    int		j 		= 0;
    int		k		= 0;
    char*	pch_function	= "VERTICES_commonInFaces_find";
    short	b_hit		= 0;
    float	f_val		= -1.;

    if(apv_verticesCommon->rows != 3 || apv_verticesCommon->cols !=1)
	ErrorExit(-1, "%s: Return VECTOR must be 3x1.\n", pch_function);
    V3_LOAD(apv_verticesCommon, -1, -1, -1);

    for(i=0; i<3; i++) {
	b_hit 	= 0;
	f_val	= -1.;
	for(j=0; j<3; j++) {
	    if(apFACE_J->v[j] == apFACE_I->v[i]) b_hit = 1; }
	if(b_hit) {
	    f_val	= apFACE_I->v[i];
	    VECTOR_ELT(apv_verticesCommon, ++k) = f_val;
	} 
    }
    return k;
}

short
FACES_Hcurvature_determineSign(
    MRIS*	apmris,
    FACE*	apFACE_O,
    FACE*	apFACE_I
) {
    //
    // PRECONDITIONS
    //	o Typically called from MRIS_Hcurvature_determineSign()
    //	o apFACE_I and apFACE_J are geometric neighbouring faces
    //	  about a given vertex.
    //
    // POSTCONDITIONS
    //	o If the faces "diverge", i.e. have face normals that point "away"
    //	  from each other, then the curvature sign return is -1.
    //	o If the faces "converge", i.e. have face normals that point "toward"
    //	  each other, then the curvature sign return is +1
    //
    // HISTORY
    // 	02 July 2007
    //	o Initial design and coding.
    //

    static int	calls			= 0;
    char*	pch_function		= "FACES_Hcurvature_determineSign";
    int		ret			= 0;
    int		vertexO			= -1;
    int		vertexI			= -1;
    VERTEX*	pVERTEX_O		= NULL;
    VERTEX*	pVERTEX_I		= NULL;
    VECTOR*	pv_O			= NULL;	// Coords of 1st common vertex
    VECTOR*	pv_normalO		= NULL; // Normal for face O
    VECTOR*	pv_OnO			= NULL; // pv_O + normal
    VECTOR*	pv_I			= NULL;	// Coords of 2nd common vertex
    VECTOR*	pv_normalI		= NULL; // Normal for face I
    VECTOR*	pv_InI			= NULL; // pv_I + normal
    VECTOR*	pv_connectOI		= NULL; // vector connecting normals
    VECTOR*	pv_commonVertices	= NULL; // Vector housing vertices that
						// are common between two
						// neighbouring faces.
    int		commonVertices		= 0;	// number of vertices in common
						// between two faces
    float	f_distNormal		= 0;
    float	f_distAntiNormal	= 0;
    int		sign			= -1;

    pv_commonVertices	= VectorAlloc(3, MATRIX_REAL);
    pv_I		= VectorAlloc(3, MATRIX_REAL);
    pv_O		= VectorAlloc(3, MATRIX_REAL);
    pv_normalI		= VectorAlloc(3, MATRIX_REAL);
    pv_normalO		= VectorAlloc(3, MATRIX_REAL);
    pv_InI		= VectorAlloc(3, MATRIX_REAL);
    pv_OnO		= VectorAlloc(3, MATRIX_REAL);
    pv_connectOI	= VectorAlloc(3, MATRIX_REAL);

    commonVertices	= VERTICES_commonInFaces_find(apFACE_O, apFACE_I,
					pv_commonVertices
					);
    if(commonVertices!=2) 
	ErrorExit(-4, 
	"%s: During call %d, the passed faces do not share an edge",
	pch_function, calls);

    vertexO		= FACE_vertexIndexAtMask_find(apFACE_O, pv_commonVertices);
    vertexI		= FACE_vertexIndexAtMask_find(apFACE_I, pv_commonVertices);
    pVERTEX_O		= &apmris->vertices[vertexO];
    pVERTEX_I		= &apmris->vertices[vertexI];     
    V3_LOAD(pv_O, pVERTEX_O->x, pVERTEX_O->y, pVERTEX_O->z);
    V3_LOAD(pv_I, pVERTEX_I->x, pVERTEX_I->y, pVERTEX_I->z);
    V3_LOAD(pv_normalO, apFACE_O->nx, apFACE_O->ny, apFACE_O->nz);
    V3_LOAD(pv_normalI, apFACE_I->nx, apFACE_I->ny, apFACE_I->nz);

    for(sign = 1; sign >= -1; sign-=2) {
	V3_SCALAR_MUL(pv_normalO, sign, pv_normalO);
	V3_SCALAR_MUL(pv_normalI, sign, pv_normalI);
        V3_ADD(pv_O, pv_normalO, pv_OnO);
        V3_ADD(pv_I, pv_normalI, pv_InI);
        V3_SUBTRACT(pv_OnO, pv_InI, pv_connectOI);

        if(sign == 1) 	f_distNormal	= V3_LEN(pv_connectOI);
	else		f_distAntiNormal = V3_LEN(pv_connectOI);
    }
	
    ret = (f_distNormal < f_distAntiNormal) ? +1 : -1;

    VectorFree(&pv_O);
    VectorFree(&pv_normalO);
    VectorFree(&pv_OnO);
    VectorFree(&pv_I);
    VectorFree(&pv_normalI);
    VectorFree(&pv_InI);
    VectorFree(&pv_connectOI);
    VectorFree(&pv_commonVertices);
    calls++;
    return ret;
}

int
VERTEX_faceAngles_determine(
    MRIS*	apmris,
    int		avertex,
    VECTOR*	apv_angle
) {
    //
    // PRECONDITIONS
    //  o <apmris> is a valid surface.
    //	o <apVERTEX> is a vertex to analyze.
    //
    // POSTCONDITIONS
    //	o The angle between each face in the <apex> normal
    //	  is determined and returned in <apv_angle>.
    //	o The caller is responsible for clearing the memory allocated to
    //	  <apv_angle>!
    //	o The number of faces processed (i.e. size of the <apv_angle
    //	  vector) is returned in the function name.
    //
    // HISTORY
    // 	30 July 2007
    //	o Initial design and coding.
    //

    char*		pch_function	= "VERTEX_faceAngles_determine";
    int			nfaces		= -1;
    float		f_angle		= 0.;
    float		f_lenApexNormal	= 0.;
    float		f_lenFaceNormal	= 0.;
    float		f_lenNormals	= 0.;
    float		f_acosArg	= 0.;
    float		f_dot		= 0.;

    int			face		= 0;
    static int		calls		= 0;
    static VECTOR*	pv_faceNormal	= NULL;
    static VECTOR*	pv_apexNormal	= NULL;

    VERTEX*		pVERTEX_apex	= NULL;
    FACE*		pFACE_side	= NULL;

    if(!calls) {
	pv_faceNormal	= VectorAlloc(3, MATRIX_REAL);
	pv_apexNormal	= VectorAlloc(3, MATRIX_REAL);
    }
    pVERTEX_apex	= &apmris->vertices[avertex];
    nfaces		= pVERTEX_apex->num;
    VECTOR_ELT(pv_apexNormal, 1)	= pVERTEX_apex->nx;
    VECTOR_ELT(pv_apexNormal, 2)	= pVERTEX_apex->ny;
    VECTOR_ELT(pv_apexNormal, 3)	= pVERTEX_apex->nz;
    f_lenApexNormal			= V3_LEN(pv_apexNormal);

    for(face=0; face<nfaces; face++) {
	pFACE_side		= &apmris->faces[pVERTEX_apex->f[face]];
    	VECTOR_ELT(pv_faceNormal, 1)	= pFACE_side->nx;
    	VECTOR_ELT(pv_faceNormal, 2)	= pFACE_side->ny;
    	VECTOR_ELT(pv_faceNormal, 3)	= pFACE_side->nz;
	f_lenFaceNormal			= V3_LEN(pv_faceNormal);
	f_lenNormals			= f_lenApexNormal * f_lenFaceNormal;
 	f_dot				= V3_DOT(pv_apexNormal, pv_faceNormal);
	errno				= 0;
// 	feclearexcept(FE_ALL_EXCEPT);
	f_acosArg			= f_dot / f_lenNormals;
        // Check on the bounds of the acos argument. Without this bounds check, 
        // it is quite possible to have 'nan' acos results, especially on 64-bit
        // builds.
        if(f_acosArg > 1.)  f_acosArg 	= 1.0;
        if(f_acosArg < -1.) f_acosArg 	= -1.0;
    	f_angle				= acos(f_acosArg);
	if(errno) {
	    f_angle			= 0.;
	    printf("%s: acos error - angle set to zero for vertex = %d, face = %d.\n", 
			pch_function, avertex, face);
	}
	VECTOR_ELT(apv_angle, face+1)	= f_angle;
    }
    calls++;
    return face;
}

int
VERTEX_faceMinMaxAngles_determine(
    MRIS*	apmris,
    int		avertex,
    int*	ap_minIndex,
    float*	apf_minAngle,
    int*	ap_maxIndex,
    float*	apf_maxAngle
) {
    //
    // PRECONDITIONS
    //  o <apmris> is a valid surface.
    //	o <apVERTEX> is a vertex to analyze.
    //
    // POSTCONDITIONS
    //	o For the given <avertex>, the minimum and maximum
    //	  face angles and their indices are returned in the
    //	  argument pointers.
    //	o The number of faces is returned in the function name.
    //
    // HISTORY
    // 	31 July 2007
    //	o Initial design and coding.
    //

    char*	pch_function		= "VERTEX_faceMinMaxAngles_determine";
    int		face			= 0;	// Face index counte
    int		nfaces			= 0;	// Number of faces at <avertex>
    float	f_faceAngle		= 0.;	// Actual face angle
    VECTOR* 	pv_faceAngles		= NULL; // vector containing angles
						// between each face normal
						// and apex normal
    VERTEX*	pVERTEX			= NULL;

    // Determine the angles between each face and the vertex normal;
    //	find the min/max angles and indices
    pVERTEX		= &apmris->vertices[avertex];
    nfaces		= pVERTEX->num;
    pv_faceAngles	= VectorAlloc(nfaces, MATRIX_REAL);
    nfaces		= VERTEX_faceAngles_determine(apmris, avertex, 
							pv_faceAngles);
    if(!nfaces)
	ErrorExit(-4, "%s: error with determining face angles.", pch_function);
    f_faceAngle		= VECTOR_ELT(pv_faceAngles, 1);
    *apf_minAngle	= f_faceAngle;
    *apf_maxAngle	= f_faceAngle;
    for(face=1; face<nfaces; face++) {
	f_faceAngle	= VECTOR_ELT(pv_faceAngles, face+1);	// base 1 index
	if(f_faceAngle < *apf_minAngle) {
	    *apf_minAngle	= f_faceAngle;
	    *ap_minIndex	= face;
	}
	if(f_faceAngle > *apf_maxAngle) {
	    *apf_maxAngle	= f_faceAngle;
	    *ap_maxIndex	= face;
	}
    }
    VectorFree(&pv_faceAngles);
    return face;
}

int
signum_eval(
    float	af
) {
    //
    // PRECONDITIONS
    //  o <af> is an input float.
    //
    // POSTCONDITIONS
    // 	o if <af> < 0, a -1 is returned, else +1 is returned
    //
   
    return af<0 ? -1 : 1;
}

int
MRIS_Hcurvature_determineSign(
    MRIS*	apmris
) {
    //
    // NOTE
    //	This function is obsolete and should not be used! Mean curvature (H)
    //	determination is now done directly when processing faces and normal
    //	angles.
    //
    // PRECONDITIONS
    //  o <apmris> is a valid surface.
    //	o MRIS_facesAtVertices_reorder()
    //
    // POSTCONDITIONS
    //	o The face geometry at each vertex is examined, and the orientation
    //	  of each face relative to its neighbor is determined as either
    //	  converging (-) or diverging (+).
    //	o If the sum of each convergence/divergence is determined to be
    //	  negative, the vertex is marked as diverging; otherwise it is
    //	  marked as converging.
    //	o Convergence is indicated with pVERTEX->undefval=-1; divergence
    //	  is marked with pVERTEX->undefval=1
    //
    // HISTORY
    // 	02 July 2007
    //	o Initial design and coding.
    //

    int 	vertex			= 0;
    int		face			= 0;
    int		nfaces			= 0;
    VERTEX*	pVERTEX			= NULL;
    int		ret			= 1;
    int		signSum			= 0;
    FACE*	pFACE_I;
    FACE*	pFACE_J;

    for(vertex=0; vertex<apmris->nvertices; vertex++) {
	MRIS_vertexProgress_print(apmris, vertex, 
		"Determining H sign for vertex faces...");
	pVERTEX			= &apmris->vertices[vertex];
	nfaces			= pVERTEX->num;
	signSum			= 0;
	for(face=0; face<nfaces; face++) {
	    pFACE_I		= &apmris->faces[pVERTEX->f[face]];
	    pFACE_J		= &apmris->faces[pVERTEX->f[(face+1)%nfaces]];
	    signSum 	+= FACES_Hcurvature_determineSign(apmris, pFACE_I, pFACE_J);
	}
	pVERTEX->undefval	= (signSum >= 0) ? 1 : -1;
    }
    return ret;
}

int
MRIS_facesAtVertices_reorder(
    MRIS*	apmris
) {
    //
    // PRECONDITIONS
    //  o <apmris> is a valid surface.
    //
    // POSTCONDITIONS
    //	o The 'f' FACE array at each vertex has its indices reordered
    //	  so that bordering face indices index (i) and index (i+1)
    //	  correspond to the actual geometric order of the faces about
    //	  each vertex.
    //	o Note that the 'f' FACE array is changed at each vertex by 
    //	  this function.
    //
    // HISTORY
    // 	02 July 2007
    //	o Initial design and coding.
    //

    int 	vertex			= 0;
    int		face			= 0;
    int		nfaces			= 0;
    int		orderedIndex		= -1;
    int		orderedFace		= -1;
    VECTOR*	pv_geometricOrderIndx	= NULL;
    VECTOR*	pv_logicalOrderFace	= NULL;
    VERTEX*	pVERTEX			= NULL;
    int		ret			= 1;
    char*	pch_function		= "MRIS_facesAtVertices_reorder";

    DebugEnterFunction(( pch_function ));
    fprintf(stderr, "\n");
    for(vertex=0; vertex<apmris->nvertices; vertex++) {
	MRIS_vertexProgress_print(apmris, vertex,
				"Determining geometric order for vertex faces...");
	pVERTEX			= &apmris->vertices[vertex];
	nfaces			= pVERTEX->num;
	pv_geometricOrderIndx	= VectorAlloc(nfaces, MATRIX_REAL);
	pv_logicalOrderFace	= VectorAlloc(nfaces, MATRIX_REAL);
	FACES_aroundVertex_reorder(apmris, vertex, pv_geometricOrderIndx);
	for(face=0; face<nfaces; face++) {
	    VECTOR_ELT(pv_logicalOrderFace, face+1) = pVERTEX->f[face];
	}
	for(face=0; face<nfaces; face++) {
	    orderedIndex	= VECTOR_ELT(pv_geometricOrderIndx, face+1);
	    orderedFace		= VECTOR_ELT(pv_logicalOrderFace, orderedIndex+1);
	    pVERTEX->f[face]	= orderedFace;
	}
	VectorFree(&pv_geometricOrderIndx);	
	VectorFree(&pv_logicalOrderFace);	
    }
    xDbg_PopStack();
    return ret;
}

int
MRIScomputeGeometricProperties(
    MRIS*	apmris
) {
    //
    // PRECONDITIONS
    //  o Needs to be called before computing discrete curvatures.
    //
    // POSTCONDITIONS
    // 	o The face array at each vertex is re-ordered in a geometric sense.
    //	o Each pair of bordering faces at each vertex are processed to
    //	  to determine overall convexity/concavity of "node".
    //
    // HISTORY
    // 03 July 2007
    //	o Initial design and coding.
    //

    int	ret	= 0;
    ret  	= MRIS_facesAtVertices_reorder(apmris);
    return ret;
}

short
FACES_aroundVertex_reorder(
    MRIS*	apmris,
    int		avertex,
    VECTOR*	pv_geometricOrder
) {
    //
    // PRECONDITIONS
    //  o <avertex> is a valid vertex on the surface.
    //	o <pll_faces> should not be allocated.
    //
    // POSTCONDITIONS
    // 	o The face indices about vertex <avertex>, starting with face 0
    //	  are returned in geometric order in the <pv_geometricOrder> vector.
    //	  By geometric order is implied that the "next" index denotes 
    //	  the next face that directly borders the current face.
    //	o The number of connected faces is returned in the function name, or
    //	  0 is there is some error.
    //
    // HISTORY
    // 	25 June 2007
    //	o Initial design and coding.
    //

    char*	pch_function	= "FACES_aroundVertex_reorder";
    VERTEX*	pVERTEX;
    int		nfaces		= 0;
    int*	pFaceIndex	= NULL;
    int		packedCount	= 1;
    int		i		= 0;
    int		I		= 0;
    int		j		= 0;
    int		k 		= 0;
    FACE*	pFACE_I;
    FACE*	pFACE_J;
    VECTOR*	pv_commonVertices	= NULL; // Vector housing vertices that
						// are common between two
						// neighbouring faces.
    int		commonVertices		= 0;	// number of vertices in common
						// between two faces
    short	b_borderFound		= 0;		

    pv_commonVertices	= VectorAlloc(3, MATRIX_REAL);
    pVERTEX		= &apmris->vertices[avertex];
    nfaces		= pVERTEX->num;
    pFaceIndex		= pVERTEX->f;

    for(i=1; i<=nfaces; i++) {
	VECTOR_ELT(pv_geometricOrder, i)= -1;
	pFACE_I	= &apmris->faces[pFaceIndex[i-1]];
    } 
    VECTOR_ELT(pv_geometricOrder, 1)	= 0;
    for(i=0; i<nfaces; i++) {
	if(packedCount == nfaces) break;
	I	= VECTOR_ELT(pv_geometricOrder, i+1);
	pFACE_I	= &apmris->faces[pFaceIndex[I]];
	for(j=0; j<nfaces; j++) {
	    k 			= (i+j) % nfaces;
	    pFACE_J		= &apmris->faces[pFaceIndex[k]];
	    commonVertices	= VERTICES_commonInFaces_find(pFACE_I, pFACE_J,
					pv_commonVertices
					);
	    if(commonVertices==2) {
		if(!VECTOR_elementIndex_find(pv_geometricOrder, k)) {
		    VECTOR_ELT(pv_geometricOrder, i+2) = k;
		    b_borderFound		= 1;
		    packedCount++;
		    break;
		}
	    }
	}
    }
    if(packedCount != nfaces)
	ErrorExit(-4, "%s: packed / faces mismatch; vertex = %d, faces = %d, packed = %d",
 			pch_function, avertex, nfaces, packedCount);
    VectorFree(&pv_commonVertices);
    return 1;	
}

float
FACES_angleNormal_find(
    MRIS*	apmris,
    FACE*	apFACE_I,
    FACE*	apFACE_J
) {
    //
    // PRECONDITIONS
    //  o The <apFACE>s should be triangles with 3 vertices each.
    //	o It is assumed (but not mandatory) that the FACES share 
    //	  at least one common vertex - or more often a common 
    //	  edge, i.e. two common vertices.
    //
    // POSTCONDITIONS
    // 	o The angle between the normals on each FACE is returned.
    //

    static int 	    calls		= 0;	// Used for vector allocation
    static VECTOR*  pv_faceNormalI	= NULL;	// Normal vector for face I
    static VECTOR*  pv_faceNormalJ	= NULL; // Normal vector for face J
    static VECTOR*  pv_crossIJ		= NULL; // Cross product of input vectors
    float	    f_faceNormalIlen	= 0.;	// Length of face normal I
    float	    f_faceNormalJlen	= 0.;	// Length of face normal J
    float	    f_faceNormalIJlen	= 0.;	// Face normal length product
    float	    f_angleNormalIJ	= 0.;	// Angle between face normals
    float	    f_acosArg		= 0.;	// Dot product arguments
    float	    f_dot		= 0.;	// Dot product
    short	    sign		= 1; 	// Angle "sign"
    char*	    pch_function	= "FACES_angleNormal_find";

    DebugEnterFunction(( "%s", pch_function));
    if(!calls) {
        pv_faceNormalI	= VectorAlloc(3, MATRIX_REAL);
        pv_faceNormalJ	= VectorAlloc(3, MATRIX_REAL);
        pv_crossIJ	= VectorAlloc(3, MATRIX_REAL);
    }
    V3_LOAD(pv_faceNormalI, apFACE_I->nx, apFACE_I->ny, apFACE_I->nz);
    V3_LOAD(pv_faceNormalJ, apFACE_J->nx, apFACE_J->ny, apFACE_J->nz);
    DebugPrint(( "%20s = %5f %20s = %5f %20s = %5f\n",
		 "Inx", apFACE_I->nx, 
		 "Iny", apFACE_I->ny, 
		 "Inz", apFACE_I->nz));
    DebugPrint(( "%20s = %5f %20s = %5f %20s = %5f\n",
		 "Jnx", apFACE_J->nx, 
		 "Jny", apFACE_J->ny, 
		 "Jnz", apFACE_J->nz));
    f_faceNormalIlen	= V3_LEN(pv_faceNormalI);
    f_faceNormalJlen	= V3_LEN(pv_faceNormalJ);
    if(f_faceNormalIlen > 1.0001 || f_faceNormalJlen > 1.0001 )
	ErrorExit(-4, "%s: face normal not unit length -- Ni: %f, Nj: %f\n", 
			pch_function, f_faceNormalIlen, f_faceNormalJlen);
    f_faceNormalIJlen	= f_faceNormalIlen * f_faceNormalJlen;
    f_dot		= V3_DOT(pv_faceNormalI, pv_faceNormalJ);
    sign		= FACES_Hcurvature_determineSign(apmris, apFACE_I, apFACE_J);
    f_acosArg		= f_dot / f_faceNormalIJlen;
    // Check on the bounds of the acos argument. Without this bounds check, 
    //	it is quite possible to have 'nan' acos results, especially on 64-bit
    //	builds.
    if(f_acosArg > 1.)	f_acosArg = 1.0;
    if(f_acosArg < -1.) f_acosArg = -1.0;
    f_angleNormalIJ	= acosf(f_acosArg) * sign;
    calls++;
    xDbg_PopStack();
    return f_angleNormalIJ;
}

float
FACES_commonEdgeLength_find(
    MRIS*	apmris,
    FACE*	apFACE_I,
    FACE*	apFACE_J
) {
    //
    // PRECONDITIONS
    //  o The <apFACE>s should be triangles with 3 vertices each.
    //	o The FACES share a common edge.
    //	o The FACES are not the same.
    //
    // POSTCONDITIONS
    //  o Length of common edge is returned.
    //  o If common vertices != 2, then function ErrorExits.
    //

    static int	      calls		= 0;
    char*	      pch_function	= "FACES_commonEdgeLength_find";
    VERTEX*	      pVERTEX_O		= NULL;	// Common vertex O
    VERTEX*	      pVERTEX_I		= NULL;	// Common vertex I
    static VECTOR*    pVECTOR_O		= NULL; // Common vertex O cart. coords 
    static VECTOR*    pVECTOR_I		= NULL; // Common vertex I cart. coords 
    static VECTOR*    pVECTOR_edgeVoVi	= NULL; // Edge Vo->Vi 
    static VECTOR*    pv_commonVertices	= NULL; // Vector housing vertices that
						// are common between two
						// neighbouring faces.
    int		      commonVertices	= 0;	// number of vertices in common
						// between two faces
    float	      f_edgeLength	= 0.;	// Length of edge v->vI

    DebugEnterFunction(( "%s", pch_function));
    if(!calls) {
        pv_commonVertices	= VectorAlloc(3, MATRIX_REAL);
        pVECTOR_O		= VectorAlloc(3, MATRIX_REAL);
        pVECTOR_I		= VectorAlloc(3, MATRIX_REAL);
        pVECTOR_edgeVoVi	= VectorAlloc(3, MATRIX_REAL);
    }
    commonVertices	= VERTICES_commonInFaces_find(
					apFACE_I,
					apFACE_J,
					pv_commonVertices
				  );
    DebugPrint(( "%20s = %5d %20s = %5d %20s = %5d\n",
		 "cV1", (int) VECTOR_ELT(pv_commonVertices, 1), 
		 "cV2", (int) VECTOR_ELT(pv_commonVertices, 2), 
		 "cV3", (int) VECTOR_ELT(pv_commonVertices, 3) ));
    if(commonVertices != 2)
	ErrorExit(-4, 
			"%s: No common edge found! <commonVertices> = %d\n",
			pch_function, commonVertices);	

    pVERTEX_O		= &apmris->vertices[(int)VECTOR_ELT(pv_commonVertices, 1)];
    pVERTEX_I		= &apmris->vertices[(int)VECTOR_ELT(pv_commonVertices, 2)];

    DebugPrint(( "(1: %d) vx  = %f, vy  = %f, vz  = %f\n", 
		(int)VECTOR_ELT(pv_commonVertices, 1),
		pVERTEX_O->x, pVERTEX_O->y, pVERTEX_O->z));
    DebugPrint(( "(1: %d) vnx = %f, vny = %f, vnz = %f\n", 
		(int)VECTOR_ELT(pv_commonVertices, 1),
		pVERTEX_O->nx, pVERTEX_O->ny, pVERTEX_O->nz));
    DebugPrint(( "(2: %d) vx  = %f, vy  = %f, vz  = %f\n", 
		(int)VECTOR_ELT(pv_commonVertices, 2),
		pVERTEX_I->x, pVERTEX_I->y, pVERTEX_I->z));
    DebugPrint(( "(2: %d) vnx = %f, vny = %f, vnz = %f\n", 
		(int)VECTOR_ELT(pv_commonVertices, 2),
		pVERTEX_I->nx, pVERTEX_I->ny, pVERTEX_I->nz));

    V3_LOAD(pVECTOR_O, pVERTEX_O->x, pVERTEX_O->y, pVERTEX_O->z);
    V3_LOAD(pVECTOR_I, pVERTEX_I->x, pVERTEX_I->y, pVERTEX_I->z);
    V3_SUBTRACT(pVECTOR_I, pVECTOR_O, pVECTOR_edgeVoVi);
    f_edgeLength 	= V3_LEN(pVECTOR_edgeVoVi);
    calls++;
    xDbg_PopStack();
    return(f_edgeLength);
}

short
MRIS_discreteKH_compute(
	MRIS*			apmris
) {
    //
    // PRECONDITIONS
    //  o A valid SURFACE with computed triangle properties.
    //
    // POSTCONDITIONS
    //	o The discrete K and H curvatures at each vertex are
    //	  computed.
    //

    char*	pch_function		= "MRIS_discreteKH_compute";

    VECTOR*	pv_geometricOrder	= NULL; // Geometrically ordered faces
    VERTEX*	pVertex			= NULL;	// Each vertex in surface
    int       	vertex			= 0;	// Vertex index number
    FACE*	pFACE_I			= NULL;	// Face I with vertex apex
    FACE*	pFACE_J			= NULL; // Face I+1 with vertex apex
    int		face			= 0;	// face counter
    int 	faceI			= 0;	// face I index
    int		faceJ			= 0;	// face J index
    int		nfaces			= 0;	// total number of faces
    int*	pFaceIndex		= NULL;	// face index array at vertex
    int		angleIndex		= -1;	// angle index
    float	f_faceAreaSum		= 0.;	// area about vertex
    float	f_angleDeficitSum	= 0.;	// angle deficit about vertex
    float	f_angleNormalIJ		= 0.;	// angle between normals
    float	f_angleNormalIJSum	= 0.;	// sum angle between normals
    float	f_edgeLength		= 0.;	// Length of edge v->vI
    double	f_K			= 0.;	// Gaussian curvature at vertex
    double	f_H			= 0.;	// Mean curvature at vertex
    float	f_Kmin			= 0.;
    float	f_Kmax			= 0.;
    float	f_Hmin			= 0.;
    float	f_Hmax			= 0.;
    float	f_Ktotal		= 0.;

    DebugEnterFunction(( "%s", pch_function));
    for (vertex = 0 ; vertex < apmris->nvertices ; vertex++) {
	MRIS_vertexProgress_print(apmris, vertex, 
				"Determining KH curvatures...");
	f_faceAreaSum		= 0.;
	f_angleDeficitSum	= 0.;
	f_angleNormalIJSum	= 0.;
	pVertex			= &apmris->vertices[vertex];
	nfaces			= pVertex->num;
	pFaceIndex		= pVertex->f;
	pv_geometricOrder	= VectorAlloc(nfaces, MATRIX_REAL);
	for(face=0; face<nfaces; face++) {
	    faceI		= face;
	    faceJ		= (face+1)%nfaces;

	    pFACE_I		=  &apmris->faces[pFaceIndex[faceI]];
	    pFACE_J		=  &apmris->faces[pFaceIndex[faceJ]];

	    f_angleNormalIJ	= FACES_angleNormal_find(apmris, pFACE_I, pFACE_J); 	
	    f_edgeLength	= FACES_commonEdgeLength_find(
					apmris,
					pFACE_I,
				  	pFACE_J
					);
	    f_faceAreaSum 	+= pFACE_I->area;
	    angleIndex		=  FACE_vertexIndex_find(pFACE_I, vertex);
	    if(angleIndex == -1)
		ErrorExit(-4, 
		"%s:\n\tangleIndex lookup failure for vertex %d, face %d", 
		pch_function, vertex, face);
	    f_angleDeficitSum	+= pFACE_I->angle[angleIndex];
	    f_angleNormalIJSum	+= f_angleNormalIJ*f_edgeLength;
	}
	VectorFree(&pv_geometricOrder);
	pv_geometricOrder	= NULL;
	f_K = 3/f_faceAreaSum   * (2*M_PI - f_angleDeficitSum);
	f_H = 0.75/f_faceAreaSum * f_angleNormalIJSum;
	apmris->vertices[vertex].K 	= f_K;
	apmris->vertices[vertex].H	= f_H;	
	if(!vertex) {
	    f_Kmin = f_Kmax = f_K;
	    f_Hmin = f_Hmax = f_H;
	}
	if(f_K > f_Kmax) f_Kmax = f_K;
	if(f_K < f_Kmin) f_Kmin = f_K;
	if(f_H > f_Hmax) f_Hmax = f_H;
	if(f_H < f_Hmin) f_Hmin = f_H;
	f_Ktotal += f_K * pVertex->area;
    }
    apmris->Kmax	= f_Kmax;	apmris->Kmin	= f_Kmin;
    apmris->Hmax	= f_Hmax;	apmris->Hmin	= f_Hmin;
    apmris->Ktotal	= f_Ktotal;
    xDbg_PopStack();
    return(NO_ERROR);
}

short
MRIS_discretek1k2_compute(
	MRIS*			apmris
) {
    //
    // PRECONDITIONS
    //  o A valid SURFACE with computed triangle properties.
    //	o A valid K and H at each vertex.
    //
    // POSTCONDITIONS
    //	o The discrete K and H curvatures at each vertex are
    //	  computed.
    //

    char*	pch_function	= "MRIS_discretek1k2_compute";
    VERTEX*	pVERTEX		= NULL;
    float	f_k1		= 0.;
    float	f_k2		= 0.;
    float	f_A		= 0.;
    float	f_B		= 0.;
    float	f_delta		= 0.;
    float	f_K		= 0.;
    float	f_H		= 0.;
    int		vertex		= 0;
    int		deltaViolations = 0;
    
    for (vertex = 0 ; vertex < apmris->nvertices ; vertex++) {
	MRIS_vertexProgress_print(apmris, vertex, 
				"Determining k1k2 curvatures...");
	pVERTEX			= &apmris->vertices[vertex];
	f_K	= pVERTEX->K;
	f_H	= pVERTEX->H;
	f_delta	= f_H*f_H - f_K;
	if(f_delta<0) {deltaViolations++; f_delta = 0.;}
	if(f_delta < 0)
	    ErrorExit(-4, "%s: f_delta = %f, vertex = %d, f_K = %f, f_H = %f\n",
			pch_function, f_delta, vertex, f_K, f_H
			);
	f_A	= f_H + sqrt(f_delta);
	f_B	= f_H - sqrt(f_delta);
	f_k1	= fabs(f_A) >= fabs(f_B) ? f_A : f_B;
	f_k2	= fabs(f_A) <= fabs(f_B) ? f_A : f_B;
	pVERTEX->k1	= f_k1;
	pVERTEX->k2	= f_k2;
    }
    cprintd("deltaViolations", deltaViolations);
    return(NO_ERROR);
}

short
MRIScomputeSecondFundamentalFormDiscrete(
	MRIS*			apmris
) {
    int 	retKH, retk1k2;
	
    retKH	= 1;
    retk1k2	= 1;
    MRIScomputeTriangleProperties(apmris);
    MRIScomputeGeometricProperties(apmris);
    retKH	= MRIS_discreteKH_compute(apmris);
    retk1k2	= MRIS_discretek1k2_compute(apmris);
    return(retKH | retk1k2);
}

int
MRISscaleCurvature(
	MRI_SURFACE* 	apmris,
  	float   	af_scale) 
{
    //
    // POSTCONDITIONS
    // o Each curvature value in apmris is scaled by:
    //
    //		(curv-f_mean)*<af_scale> + f_mean
    //
    //   where f_mean is the mean of all the surface curvatures
    //

    VERTEX* 	pvertex;
    int  	vno;
    int  	vtotal;
    double 	f_mean;

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

// To mrisurf.c ^^^^^^^^^^^^^^^^^^

short
surfaceIntegrals_compute(
	MRIS*			apmris,
	e_surfaceIntegral	aeSI_type,
	float*			apf_surfaceIntegral,
	float*			apf_meanIntegral,
	int*			ap_vertices,
	float*			apf_areaNormalizedIntegral,
	float*			apf_area
) {
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

  short b_ret	= 0;

  switch(aeSI_type) {
    case e_natural:
	b_ret		= MRIS_surfaceIntegral_compute(
				apmris, ap_vertices, apf_area,
				f_allVertices,
				f_pass,
				apf_surfaceIntegral,
				apf_meanIntegral,
				apf_areaNormalizedIntegral);
    break;
    case e_rectified:
	b_ret		= MRIS_surfaceIntegral_compute(
				apmris, ap_vertices, apf_area,
				f_allVertices,
				f_absCurv,
				apf_surfaceIntegral,
				apf_meanIntegral,
				apf_areaNormalizedIntegral); 
    break;
    case e_pos:
	b_ret		= MRIS_surfaceIntegral_compute(
				apmris, ap_vertices, apf_area,
				f_greaterEqualZero,
				f_absCurv,
				apf_surfaceIntegral,
				apf_meanIntegral,
				apf_areaNormalizedIntegral);
    break;
    case e_neg:
	b_ret		= MRIS_surfaceIntegral_compute(
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
	MRIS*			apmris,
	e_secondOrderType	aesot,
	s_MINMAX*		aps_minMax
) {
    //
    // PRECONDITIONS
    //	o Typically called only by MRIS_minMaxCurve_report(...)
    //	o Either surfFunc or k1k2Func should be NULL -- this function
    //	  will call whichever of the above is non-NULL.
    //	o The <s_MINMAX> structure is used to communicate back to the
    //	  calling function. If the boolean b_{min,max}Violation fields are set
    //	  to TRUE, the caller should interpret this to mean that an explicit
    //	  lookup failed.
    //	o Similarly, if the <s_MINMAX> b_{min,max}Test fields are true as
    //	  set by caller, this function will perform an explicit min/max
    //	  search across the surface.
    //
    // POSTCONDITIONS
    //	o The <s_MINMAX> structure is populated, and pre- post- behaviour
    //	  is set by the <s_MINMAX> boolean control flags.
    //

    int		vmin		= 0;
    int 	vmax		= 0;
    float	f_minExplicit	= 0.;
    float	f_maxExplicit	= 0.;
    char*	pch_function	= "MRIS_minMaxCurvature_analyze";

    DebugEnterFunction (( pch_function ));
    aps_minMax->b_minViolation	= 0;	// Assume no "violations" for min
    aps_minMax->b_maxViolation	= 0;	// Assume no "violations" for max
    MRISminMaxCurvaturesSearch(	apmris, 
				&vmin, 		&vmax, 
				&f_minExplicit, &f_maxExplicit);
    aps_minMax->vertexMin	= vmin;
    aps_minMax->vertexMax	= vmax;
    if (aps_minMax->b_minTest && aps_minMax->f_min != f_minExplicit) {
     	fprintf(stderr, "\n\nWARNING\n%10s%-40s%f\n", Gppch[aesot], 
			" lookup   min:", aps_minMax->f_min);
      	fprintf(stderr, "%10s%-40s%f\tvertex = %d\n", Gppch[aesot],
			" explicit min:", f_minExplicit, vmin);
	aps_minMax->b_minViolation	= 1;
	aps_minMax->f_min 		= f_minExplicit;
    }
    if (aps_minMax->b_maxTest && aps_minMax->f_max != f_maxExplicit) {
     	fprintf(stderr, "%10s%-40s%f\n", Gppch[aesot], 
			" lookup   max:", aps_minMax->f_max);
      	fprintf(stderr, "%10s%-40s%f\tvertex = %d\n", Gppch[aesot],
			" explicit max:", f_maxExplicit, vmax);
	aps_minMax->b_maxViolation	= 1;
	aps_minMax->f_max 		= f_minExplicit;
    }
    xDbg_PopStack();
    return 1;
}

short
MRIS_curvatures_prepare(
	MRIS*			apmris,
	e_secondOrderType	aesot
) {
    //
    // PRECONDITIONS
    //	o For second order type curves, <aprmis> *must* have been 
    //	  pre-processed with a call to one of the second order 
    //	  functions:
    //
    //		MRIScomputeSecondFundamentalForm(...)
    //		MRIScomputeSecondFundamentalFormDiscrete(...)
    //
    //	  Note that there is no way for this function to know that this
    //	  preprocessing has occurred.
    //  o The <aesot> "second order type" determines which particular curvature
    //	  value to process.
    //	o This is the main entry point for preparing a surface for subsequent
    //	  analysis.
    //
    // POSTCONDITIONS
    //  o Depending on <aesot>, the <apmris> surface is prepared for subsequent
    //	  analysis.
    //	o The return value from the surface prepation function is returned, or
    //	  '-1' if no preparation was performed.
    //

    char*	pch_function	= "MRIS_curvatures_prepare";
    short	ret		= -1;

    DebugEnterFunction (( pch_function ));

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
	    ret = MRISusePrincipleCurvatureFunction(apmris, 
			f_sharpnessCurvature);
	break;
	case e_C:
	    ret = MRISusePrincipleCurvatureFunction(apmris, 
			f_sharpnessCurvature);
        break;
	case e_BE:
	    ret = MRISusePrincipleCurvatureFunction(apmris,
			f_bendingEnergyCurvature);
        break;
	case e_SI:
	    ret = MRISusePrincipleCurvatureFunction(apmris, 
			f_shapeIndexCurvature);
        break;
  	case e_Raw:
  	case e_Normal:
  	case e_Scaled:
  	case e_ScaledTrans:
    	break;
    }
    xDbg_PopStack();
    return ret;
}

short
MRIS_minMaxCurve_report(
	MRIS*   		apmris,
  	e_secondOrderType 	aesot,
	char*			apch_report
) {
    //
    // PRECONDITIONS
    //	o MRIS_curvatures_prepare(...) *must* have been called for <aesot>.
    //  o The <aesot> "second order type" determines which particular curvature
    //	  value to process.
    //	o This is the main entry point for analyzing the min/max curves. It is
    //	  typically called by any higher-level function that wants to determine
    //	  the min/max for a curvature function type (as specified in <aesot>).
    //
    // POSTCONDITIONS
    //	o Calls MRIS_minMaxCurvature_analyze(...) to find the min/max
    //	  curvatures and vertex "occurrences".
    //  o Depending on <aesot>, data is printed to stdout (and
    //    output file, if spec'd).
    //	o Essentially, function determines the minimum and maximum values of
    //	  a particular curvature by explicit search. It also checks, for 
    //	  curvature cases that record min/max in the <apmris> structure, that 
    //	  the embedded min/max matches the searched min/max.
    //	o The report itself is returned in the apch_report string.
    //

    s_MINMAX	s_minMax;
    char*	pch_function	= "MRIS_minMaxCurve_report";

    DebugEnterFunction (( pch_function ));

    s_minMax.f_min		= apmris->min_curv;
    s_minMax.f_max		= apmris->max_curv;
    s_minMax.vertexMin		= -1;
    s_minMax.vertexMax		= -1;
    s_minMax.b_minTest		= 1;
    s_minMax.b_maxTest		= 1;
    s_minMax.b_minViolation	= 0;
    s_minMax.b_maxViolation	= 0;

    if(aesot == e_Gaussian) 
	{s_minMax.f_min	= apmris->Kmin; s_minMax.f_max	= apmris->Kmax;}
    if(aesot == e_Mean) 
	{s_minMax.f_min	= apmris->Hmin; s_minMax.f_max	= apmris->Hmax;}

    MRIS_minMaxCurvature_analyze(apmris, aesot, &s_minMax);

    if(aesot == e_Gaussian) {	    
	if(s_minMax.b_minViolation) apmris->Kmin = s_minMax.f_min;
	if(s_minMax.b_maxViolation) apmris->Kmax = s_minMax.f_max;
    }
    if(aesot == e_Mean) {
	if(s_minMax.b_minViolation) apmris->Hmin = s_minMax.f_min;
	if(s_minMax.b_maxViolation) apmris->Hmax = s_minMax.f_max;
    }

    if(s_minMax.b_minViolation) apmris->min_curv = s_minMax.f_min;
    if(s_minMax.b_maxViolation) apmris->max_curv = s_minMax.f_max;

    strcpy(apch_report, "");
    sprintf(apch_report, "%*s%-*s", 	G_leftCols, 	Gppch[aesot], 
					G_rightCols, 	" Min:");
    sprintf(apch_report, "%s%12.5f at vertex %-8d\n", 	apch_report,
		s_minMax.f_min, s_minMax.vertexMin);
    sprintf(apch_report, "%s%*s%-*s",	apch_report,
					G_leftCols, 	Gppch[aesot], 
					G_rightCols, 	" Max:");
    sprintf(apch_report, "%s%12.5f at vertex %-8d\n", 	apch_report, 
		s_minMax.f_max, s_minMax.vertexMax);
    if (GpFILE_allLog)
      fprintf(GpFILE_allLog, "min = %f\nmax = %f\n",
              s_minMax.f_min, s_minMax.f_max);
    xDbg_PopStack();
    return 1;
}

short
MRIS_surfaceIntegrals_report(
	MRIS*   		apmris,
	e_secondOrderType	aesot,
	char*			apch_report
) {
    //
    // PRECONDITIONS
    //  o The 'curv' field of <apmris> must contain the particular curvature
    //    map value to integrate.
    //  o The <apch_curvName> is a string that is prefixed to each output line
    //	  denoting the curve being processed.
    //
    // POSTCONDITIONS
    //  o This function is the typical entry point to performing a surface
    //	  integral. The 'report' in the function name is meant to convey
    //	  that this function performs the integral, and *prints* the results
    //	  to the output device.
    //	o The function performs/prints the following surface integral
    //	  functions:
    //		- "natural":	no change/filtering on the vertex curv values
    //		- "abs":	the fabs(...) of each vertex curv value
    //		- "pos":	process only the positive vertex curv values 
    //		- "neg":	process only the negative vertex curv values
    //	o In addition, the following qualifiers are also evaluated:
    //		- "Mean":	the integral normalized to number of vertices
    //		- "AreaNorm":	the integral value normalized to unit surface
    //	o The report itself is returned in the apch_report string.
    //

    // Surface integral variables
    float f_SInatural		= 0.0;
    float f_SIabs		= 0.0;
    float f_SIpos		= 0.0;
    float f_SIneg		= 0.0;
    float f_SInaturalMean	= 0.0;	int SInaturalVertices	= 0;
    float f_SIabsMean		= 0.0; 	int SIabsVertices	= 0;
    float f_SIposMean		= 0.0;	int SIposVertices	= 0;
    float f_SInegMean		= 0.0;  int SInegVertices	= 0;
  
    float f_SInaturalAreaNorm	= 0.0;	float f_SInaturalArea	= 0.0;
    float f_SIabsAreaNorm	= 0.0; 	float f_SIabsArea	= 0.0;
    float f_SIposAreaNorm	= 0.0;	float f_SIposArea	= 0.0;
    float f_SInegAreaNorm	= 0.0;  float f_SInegArea	= 0.0;
    char* pch_curveName;
    char* pch_function		= "MRIS_surfaceIntegrals_report";

    DebugEnterFunction (( pch_function ));
    pch_curveName		= Gppch[aesot];

    surfaceIntegrals_compute(apmris, e_natural,
				&f_SInatural, 		 	
				&f_SInaturalMean, 	&SInaturalVertices,
				&f_SInaturalAreaNorm,	&f_SInaturalArea);

    surfaceIntegrals_compute(apmris, e_rectified, 	 	
				&f_SIabs, 		 	
				&f_SIabsMean,		&SIabsVertices,
				&f_SIabsAreaNorm,	&f_SIabsArea);
 
    surfaceIntegrals_compute(apmris, e_pos, 		 	
				&f_SIpos, 		 	
				&f_SIposMean,		&SIposVertices,
				&f_SIposAreaNorm,	&f_SIposArea);
  
    surfaceIntegrals_compute(apmris, e_neg, 		 	
				&f_SIneg, 		 	
				&f_SInegMean,		&SInegVertices,
				&f_SInegAreaNorm,	&f_SInegArea);

    strcpy(apch_report, "");

    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Average Vertex Separation:");
    sprintf(apch_report, "%s%12.5f +- %2.5f mm\n", apch_report,
		apmris->avg_vertex_dist, apmris->std_vertex_dist);

    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Total number of vertices:");
    sprintf(apch_report, "%s%12.5d\n", apch_report,
			apmris->nvertices);
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Total surface area:");
    sprintf(apch_report, "%s%12.5f mm^2\n", 	apch_report, 
			apmris->total_area);
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Average Vertex Area:");
    sprintf(apch_report, "%s%12.5f mm^2\n",	apch_report, 
			apmris->avg_vertex_area);

    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Natural Surface Integral:");
    sprintf(apch_report, "%s%12.5f\n", 	apch_report,
			f_SInatural);
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Rectified Surface Integral:");
    sprintf(apch_report, "%s%12.5f\n",	apch_report, 
			f_SIabs);
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Positive Surface Integral:");
    sprintf(apch_report, "%s%12.5f\n", 	apch_report,
			f_SIpos);
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Negative Surface Integral:");
    sprintf(apch_report,"%s%12.5f\n", 	apch_report,
			f_SIneg);

    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Mean Natural Surface Integral:");
    sprintf(apch_report, "%s%12.5f across %d (%05.2f%s) vertices\n",
			apch_report, 
			f_SInaturalMean, SInaturalVertices, 
			100 * (float)SInaturalVertices / apmris->nvertices, 
			"%");
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Mean Rectified Surface Integral:");
    sprintf(apch_report, "%s%12.5f across %d (%05.2f%s) vertices\n", 
			apch_report,
			f_SIabsMean, SIabsVertices, 
			100 * (float)SIabsVertices / apmris->nvertices, "%");
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Mean Positive Surface Integral:");
    sprintf(apch_report, "%s%12.5f across %d (%05.2f%s) vertices\n", 
			apch_report,
			f_SIposMean, SIposVertices, 
			100 * (float)SIposVertices / apmris->nvertices, "%");
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " Mean Negative Surface Integral:");
    sprintf(apch_report, "%s%12.5f across %d (%05.2f%s) vertices\n", 
			apch_report,
			f_SInegMean, SInegVertices, 
			100 * (float)SInegVertices / apmris->nvertices, "%");

    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " AreaNorm Natural Surface Integral:");
    sprintf(apch_report, "%s%12.5f across %f (%05.2f%s) mm^2\n", 
			apch_report,
			f_SInaturalAreaNorm, f_SInaturalArea, 
			100 * f_SInaturalArea / apmris->total_area, "%");
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " AreaNorm Rectified Surface Integral:");
    sprintf(apch_report, "%s%12.5f across %f (%05.2f%s) mm^2\n", 
			apch_report,
			f_SIabsAreaNorm, f_SIabsArea, 
			100 * f_SIabsArea / apmris->total_area, "%");
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " AreaNorm Positive Surface Integral:");
    sprintf(apch_report, "%s%12.5f across %f (%05.2f%s) mm^2\n", 
			apch_report,
			f_SIposAreaNorm, f_SIposArea, 
			100 * f_SIposArea / apmris->total_area, "%");
    sprintf(apch_report, "%s%10s%-40s", apch_report,
			pch_curveName, " AreaNorm Negative Surface Integral:");
    sprintf(apch_report, "%s%12.5f across %f (%05.2f%s) mm^2\n", 
			apch_report,
			f_SInegAreaNorm, f_SInegArea, 
			100 * f_SInegArea / apmris->total_area, "%");
    xDbg_PopStack();
    return 1;
}

int
MRIS_curvatureStats_analyze(
	MRIS*			apmris,
	e_secondOrderType	aesot
) {
    //
    // PRECONDITIONS
    //	o <aprmis> is valid
    // 	o <aesot> denotes the enumerated surface type to analyze.
    //
    // POSTCONDITIONS
    //  o For the specific <aesot>, the following are analyzed:
    //		- min/max values and vertices (*)
    //		- Surface intergrals are performed and reported
    //		- histograms are prepared (*)
    //	  (*) if the appropriate flag has been set by the main
    //	  user on the command line.
    //
    // SIDE EFFECTS
    //	o Modifies the current global Gf_mean and Gf_sigma
    //

    char 	pch_text[65536];
    char	pch_minMaxReport[65536];
    char	pch_surfaceIntegralReport[65536];
    char*	pch_function	= "MRIS_curvatureStats_analyze";

    DebugEnterFunction (( pch_function ));

    strcpy(pch_text, "");	
    strcpy(pch_minMaxReport, "");	
    strcpy(pch_surfaceIntegralReport, "");	

    // Perform the analyses and prepare reports
    MRIS_curvatures_prepare(apmris, aesot);
    if(Gb_minMaxShow) MRIS_minMaxCurve_report(apmris, aesot, pch_minMaxReport);
    MRIS_surfaceIntegrals_report(apmris, aesot, pch_surfaceIntegralReport);

    // Now, dump the reports to stdout
    //	First the mean/sigma results
    Gf_mean = MRIScomputeAverageCurvature(apmris, &Gf_sigma);
    if(aesot == e_Raw)
        sprintf(pch_text, 
	"\n%s <mean> +- <std> (using '%s.%s'):",
          Gppch[aesot], hemi, curv_fname);
    else
        sprintf(pch_text, 
	"\n%s <mean> +- <std> (using '%s.%s'):",
          Gppch[aesot], hemi, surf_name);
    fprintf(stdout, "%-50s", pch_text);
    fprintf(stdout, " %12.5f +- %2.4f mm\n", Gf_mean, Gf_sigma);
    if (GpFILE_allLog)
        fprintf(GpFILE_allLog, "mean/sigma = %20.4f +- %2.4f\n", 
				Gf_mean, Gf_sigma);

    // Now the min/max report
    if(Gb_minMaxShow) fprintf(stdout, "%s", pch_minMaxReport);

    // The surface integral report
    fprintf(stdout, "%s", pch_surfaceIntegralReport);

    // and any histograms...
    if(Gb_histogram) histogram_wrapper(apmris, aesot);

    xDbg_PopStack();
    return 1;
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
MRISusePrincipleCurvature(
	MRIS*		apmris,
	int		(*f)(MRIS* apmris)
) {
    //
    // 	This is a thin wrapper/dispatch about an underlying
    //	MRISuse[Gaussian|Mean|k1|k2]Curvature() function.
    //
    int	ret	= 0;

    ret		= f(apmris);

    return(ret);
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

short
MRIS_surfaceIntegral_compute(
	MRI_SURFACE* 	pmris, 
	int* 		p_verticesCounted,
	float*		pf_areaCounted,
	short		(*fcond)(VERTEX*	pv),
	float		(*fv)	(VERTEX*	pv),
	float*		pf_surfaceIntegral,
	float*		pf_meanIntegral,
	float*		pf_areaNormalizedIntegral	
) {
  //
  // DESCRIPTION
  //	This function computes a surface integral across the 
  //	vertex curvature. The resultant sum is normalized to the
  //	number of the vertices that were counted and also normalized
  //	to the total area of the counted vertices.
  //
  //	Before processing a particular vertex curvature, the vertex
  //	is first checked by the (*fcond) function for eligibility. 
  //	
  //	The actual curvature that is integrated can also be shaped
  //	by the (*fcurv) function.
  //
  //	Global low pass filtering (on actual curvature values and
  //	Gaussian test) is always implemented. 
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
  //	o The integral of the curvature value over the surface is returned 
  //	  according to various criteria.
  //	o <pf_surfaceIntegral> is the actual integral for the given curvatures
  //	  in the vertex 'curv'.
  //	o The <pf_meanIntegral> and the <pf_areaNormalizedIntegral> are the
  //	  surface intergrals normalized to either the number of vertices or
  //	  the integral area.
  //	o If the integral could not be computed, -1 is returned.
  //	o If the global <Gb_lowPassFilter> is set, only count vertices
  //	  where abs 'curv' is less than <Gb_lowPassFilter>.
  //	o If the global <Gb_lowPassFilterGaussian> is set, only count
  //	  vertices where the abs Gaussian at that vertex is less than
  //	  <Gb_lowPassFilterGaussian>.
  //

  VERTEX*	pv ;
  int       	vno ;
  double    	f_total, f_n, f_totalArea ;
  short		b_canCount		= 1;
  short		b_ret			= 0;

  *p_verticesCounted = 0;
  for (f_n = f_total =f_totalArea = 0.0, vno = 0 ; vno < pmris->nvertices ; vno++) {
    b_canCount = 1;
    pv = &pmris->vertices[vno] ;
    if (pv->ripflag)
      continue ;
 
    if(Gb_lowPassFilter) {
	if(fabs(pv->curv)<fabs(Gf_lowPassFilter))
	    b_canCount &= 1;
	else
	    b_canCount &= 0;
    }
    if(Gb_highPassFilter) {
	if(fabs(pv->curv)>=fabs(Gf_highPassFilter))
	    b_canCount &= 1;
	else
	    b_canCount &= 0;
    }
    if(Gb_lowPassFilterGaussian) {
	if(fabs(pv->K)<fabs(Gf_lowPassFilterGaussian))
	    b_canCount &= 1;
	else
	    b_canCount &= 0;
    }
    if(Gb_highPassFilterGaussian) {
	if(fabs(pv->K)>=fabs(Gf_highPassFilterGaussian))
	    b_canCount &= 1;
	else
	    b_canCount &= 0;
    }
    if(b_canCount && (*fcond)(pv)) {
        f_total 		+= ((*fv)(pv) * pv->area) ;
        f_n 			+= 1.0 ;
	(*pf_areaCounted) 	+= pv->area;
	(*p_verticesCounted)++;
    }
  }
  if(f_n > 1) {
    *pf_meanIntegral		= f_total / f_n;
    *pf_areaNormalizedIntegral	= f_total / (*pf_areaCounted);
    b_ret 			= 1;
  }
  if(Gb_postScale){	
    *pf_meanIntegral 		*= Gf_postScale;
    *pf_areaNormalizedIntegral	*= Gf_postScale;
  }
  *pf_surfaceIntegral		= f_total;
  return(b_ret);
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

void
OFSP_create(
  char*  apch_prefix,
  char*  apch_suffix,
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
  //  <hemi>.<surface>.<pch_suffix>
  //

  switch (ae_OFSP) {
  case e_None:
    if(output_fname != NULL)
        sprintf(apch_prefix, output_fname);
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
    ErrorExit(ERROR_SIZE, "%s: f_maxCurv (%f) < f_minCurv (%f)",
              Progname, f_maxCurv, f_minCurv);

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

  int  		vno  		= 0;
  float* 	pf_curvature 	= NULL;
  int  		i  		= 0;
  int  		j  		= 0;
  int  		start  		= 0;
  int  		count  		= 0;
  int  		totalCount 	= 0;
  int  		nvertices 	= 0;

  int  		b_onLeftBound = 0; 	// Don't trigger higher order float
  int  		b_onRightBound = 0; 	// comparison on left/right bounds

  double 	l_curvature; 		// These three variables
  double 	l_leftBound; 		// are scaled up and truncated
  double 	l_rightBound; 		// to minimise rounding errors

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
    cprintf("Setting rectified low pass filter to", Gf_lowPassFilter);
  } else if (!stricmp(option, "-shapeIndex")) {
    Gb_shapeIndex		= 1;
    cprints("Toggling shape index map on", "ok");
  } else if (!stricmp(option, "-lowPassFilterGaussian")) {
    Gb_lowPassFilterGaussian	= 1;
    Gf_lowPassFilterGaussian	= atof(argv[2]);
    nargs			= 1;
    cprintf("Setting rectified low pass Gaussian filter", 
		Gf_lowPassFilterGaussian);
  } else if (!stricmp(option, "-highPassFilterGaussian")) {
    Gb_highPassFilterGaussian	= 1;
    Gf_highPassFilterGaussian	= atof(argv[2]);
    nargs			= 1;
    cprintf("Setting rectified high pass Gaussian filter", 
		Gf_highPassFilterGaussian);
  } else if (!stricmp(option, "-highPassFilter")) {
    Gb_highPassFilter		= 1;
    Gf_highPassFilter		= atof(argv[2]);
    nargs			= 1;
    cprintf("Setting rectified high pass filter", Gf_highPassFilter);
  } else if (!stricmp(option, "-postScale")) {
    Gb_postScale		= 1;
    Gf_postScale		= atof(argv[2]);
    nargs			= 1;
    cprintf("Setting post scale factor", Gf_postScale);
  } else if (!stricmp(option, "-discrete")) {
    Gb_discreteCurvaturesUse	= 1;
    cprints("Using discrete curvature calculations", "ok");
  } else if (!stricmp(option, "-continuous")) {
    Gb_discreteCurvaturesUse	= 0;
    cprints("Using continuous curvature calculations", "ok");
  } else if (!stricmp(option, "-writeCurvatureFiles")) {
    Gb_writeCurvatureFiles	= 1;
    cprints("Toggling save flag on curvature files", "ok");
  } else switch (toupper(*option)) {
  case 'O':
    output_fname  		= argv[2] ;
    Gb_output2File 		= 1;
    Gb_writeCurvatureFiles	= 1;
    nargs   			= 1 ;
    cprints("Outputting results using filestem", output_fname) ;
    cprints("Toggling save flag on curvature files", "ok");
    break ;
  case 'F':
    pch_text = argv[2];
    strcpy(surf_name, pch_text);
    nargs = 1 ;
    cprints("Setting surface to '%s'", surf_name);
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
    if (argc == 2)
      print_usage();
    Gb_histogram = 1;
    Gb_histogramPercent = 0;
    Gb_minMaxShow = 1;
    G_bins   = atoi(argv[2]);
    nargs   = 1 ;
    cprintf("Setting curvature histogram bins", G_bins);
    if (G_bins <=0 )
      ErrorExit(ERROR_BADPARM, "%s: Invalid bin number.\n", Progname);
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
    fprintf(stderr, "Unknown option '%s'. Looking for help? Try with '-u'.\n", argv[1]) ;
    exit(1) ;
    break ;
  }
  return(nargs) ;
}

void
outputFileNames_create(void) {
  //
  // POSTCONDITIONS
  // o All necessary (depending on user flags) file names are created.
  //
  OFSP_create(Gpch_log,  	Gpch_logS, 		e_None);
  OFSP_create(Gpch_allLog, 	Gpch_allLogS, 		e_Full);
  OFSP_create(Gpch_rawHist, 	Gpch_rawHistS, 		e_Full);
  OFSP_create(Gpch_normHist, 	Gpch_normHistS, 	e_Full);
  OFSP_create(Gpch_normCurv, 	Gpch_normCurvS, 	e_Partial);
  OFSP_create(Gpch_KHist,  	Gpch_KHistS, 		e_Full);
  OFSP_create(Gpch_KCurv,  	Gpch_KCurvS, 		e_Partial);
  OFSP_create(Gpch_HHist, 	Gpch_HHistS, 		e_Full);
  OFSP_create(Gpch_HCurv, 	Gpch_HCurvS, 		e_Partial);
  OFSP_create(Gpch_K1Hist,  	Gpch_K1HistS, 		e_Full);
  OFSP_create(Gpch_K1Curv,  	Gpch_K1CurvS, 		e_Partial);
  OFSP_create(Gpch_K2Hist,  	Gpch_K2HistS, 		e_Full);
  OFSP_create(Gpch_K2Curv,  	Gpch_K2CurvS, 		e_Partial);
  OFSP_create(Gpch_SHist,  	Gpch_SHistS, 		e_Full);
  OFSP_create(Gpch_SCurv,  	Gpch_SCurvS, 		e_Partial);
  OFSP_create(Gpch_CHist,  	Gpch_CHistS, 		e_Full);
  OFSP_create(Gpch_CCurv,  	Gpch_CCurvS, 		e_Partial);
  OFSP_create(Gpch_BEHist,  	Gpch_BEHistS, 		e_Full);
  OFSP_create(Gpch_BECurv,  	Gpch_BECurvS, 		e_Partial);
  OFSP_create(Gpch_SIHist,  	Gpch_SIHistS, 		e_Full);
  OFSP_create(Gpch_SICurv,  	Gpch_SICurvS, 		e_Partial);
  OFSP_create(Gpch_scaledHist,	Gpch_scaledHistS,	e_Full);
  OFSP_create(Gpch_scaledCurv,	Gpch_scaledCurvS,	e_Partial);
}

void
outputFiles_open(void) {
  //
  // POSTCONDITIONS
  // o All necessary (depending on user flags) files are opened.
  //

  if (Gb_output2File) {
    fprintf(stdout,"\n\tFiles processed for this curvature:\n");
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

  sprintf(pch_synopsis, "\n\
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
          Principle curvature (K, H, k1 and k2) calculations on a surface \n\
	  structure can also be performed, as well as several functions \n\
	  derived from k1 and k2. \n\
 \n\
          Finally, all output to the console, as well as any new \n\
          curvatures that result from the above calculations can be \n\
          saved to a series of text and binary-curvature files. \n\
 \n\
 	  PRINCIPLE CURVATURES AND FUNCTIONS \n\
 \n\
	  Given a surface file, 'mris_curvature_stats' can also compute \n\
	  all the principle curvatures relating to each point on the \n\
	  surface as well as several functions of these principle \n\
	  curvatures. \n\
 \n\
	  To calculate these principle curvatures, use a '-G' flag on the \n\
	  command line. In such an instance, you do need (nor probably even \n\
	  want) any <curvFile> arguments. The following principle curvatures \n\
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
	 \n\
	  Note that the SI is not calculated by default due to some issues \n\
	  with atan calculations. Use the '--shapeIndex' flag on the command \n\
	  line to force SI. \n\
 \n\
	  In addition, if a '--writeCurvatureFiles' is specified, each of the \n\
	  above are written to a FreeSurfer format curvature file. \n\
 \n\
	  Note that there are two methods to calculate principle curvatures. \n\
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
    OPTIONS \n\
 \n\
    [-a <numberOfAverages] \n\
 \n\
          Average the curvature <numberOfAverages> times. \n\
 \n\
    [-G] [--discrete] [--continuous] \n\
 \n\
	  Calculate the principle curvatures ('-G'). The calculation method \n\
	  used can be specified explicitly as '--discrete' which is slower \n\
	  but more accurate, or as '--continuous' which will result in \n\
	  measures similar to 'mris_anatomical_stats'. \n\
 \n\
	  Note that the default assumes '--discrete'. If you specify both \n\
	  '--continuous' and '--discrete', the system will revert to its \n\
	  default, i.e. '--discrete', behaviour. \n\
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
    [-l <labelFileName>] \n\
 \n\
          Constrain statistics to the region defined in <labelFileName>. \n\
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
          create and save the principle curvature based files: \n\
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

  fprintf(stdout,pch_synopsis);
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
