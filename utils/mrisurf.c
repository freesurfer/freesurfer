#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

#include "diag.h"
#include "error.h"
#include "utils.h"
#include "macros.h"
#include "fio.h"
#include "mri.h"
#include "mri2.h"
#include "mrisurf.h"
#include "matrix.h"
#include "proto.h"
#include "stats.h"
#include "timer.h"
#include "const.h"
#include "mrishash.h"
#include "icosahedron.h"
#include "tritri.h"
#include "timer.h"
#include "chklc.h"
#include "mri_identify.h"
#include "colortab.h"
#include "tags.h"
#include "selxavgio.h"
#include "machine.h"
#include "tags.h"
#include "talairachex.h"
/*---------------------------- STRUCTURES -------------------------*/

/*---------------------------- CONSTANTS -------------------------*/

#define MAX_NBRS 10000
#define REPULSE_K   1.0
#define REPULSE_E   0.5

#define NO_NEG_DISTANCE_TERM  0
#define ONLY_NEG_AREA_TERM    1
#define ANGLE_AREA_SCALE      0.0

#define ORIG_AREAS          0
#define CURRENT_AREAS       1
#define AVERAGE_AREAS       0

#define CURV_SCALE          1000.0
#define MIN_DT_SCALE        0.01
#if 0
#define MAX_SMALL           10
#define TOTAL_SMALL         (4*MAX_SMALL)
#else
#define MAX_SMALL           50000
#define TOTAL_SMALL         15000
#endif

#define METRIC_SCALE        1

#define MAX_NBHD_SIZE       200
#define MAX_NEG_AREA_PCT    0.005f

/* this definition should come from mrisurf.h 
 * #define MRIS_ASCII_FILE     1
*/

static double NEG_AREA_K=20.0 ; /* was 200 */
/* limit the size of the ratio so that the exp() doesn't explode */
#define MAX_NEG_RATIO       (400 / NEG_AREA_K)
#define MAX_ASYNCH_MM       0.3
#define MAX_ASYNCH_NEW_MM   0.3

#define NOT_USED                       0
#define USED_IN_ORIGINAL_TESSELLATION  1
#define USED_IN_NEW_TESSELLATION       2

typedef struct
{
  int   vno1, vno2 ;
  float len ;
  short used ;
} EDGE ;


/*------------------------ STATIC PROTOTYPES -------------------------*/

static int mrisReadAsciiCurvatureFile(MRI_SURFACE *mris, char *fname) ;
static int mrisAverageSignedGradients(MRI_SURFACE *mris, int num_avgs) ;
#if 0
static int mrisAverageWeightedGradients(MRI_SURFACE *mris, int num_avgs) ;
#endif
int MRISrestoreExtraGradients(MRI_SURFACE *mris) ;
static int mrisComputePositioningGradients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static int mrisFindGrayWhiteBorderMean(MRI_SURFACE *mris, MRI *mri) ;
static int mrisDumpDefectiveEdge(MRI_SURFACE *mris, int vno1, int vno2) ;
static int mrisMarkBadEdgeVertices(MRI_SURFACE *mris, int mark) ;
static int mrisCheckSurface(MRI_SURFACE *mris) ;
#if 0
static int mrisComputeCanonicalBasis(MRI_SURFACE *mris, int fno,
                                     double origin[3],double e0[3],
                                     double e1[3]);
#endif
static int mrisInitializeNeighborhood(MRI_SURFACE *mris, int vno) ;
static int mrisSetVertexFaceIndex(MRI_SURFACE *mris, int vno, int fno) ;
static int isFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2) ;
static int findFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2) ;
static int mrisAddFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2) ;

static int mrisComputeCanonicalEdgeBasis(MRI_SURFACE *mris, EDGE *edge1,
                                         EDGE *edge2, double origin[3],
                                         double e0[3], double e1[3]);
static int mrisFindSecondNeighborhood(MRI_SURFACE *mris, int vno, int *nbrs, int *num_nbrs) ;
#if 0
static int mrisDumpTriangle(MRI_SURFACE *mris, int fno) ;
static int mrisDilateAmbiguousVertices(MRI_SURFACE *mris, int mark,int ndil) ;
static int triangleNeighbors(MRI_SURFACE *mris, int fno1, int fno2) ;
#endif
static int triangleMarked(MRI_SURFACE *mris, int fno) ;
static int mrisScaleMaxDimension(MRI_SURFACE *mris, float maxr) ;
static int mrisCalculateOriginalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                     float *px, float *py, float *pz) ;
static int mrisCalculateFaceCentroid(MRI_SURFACE *mris, int fno, 
                                     float *px, float *py, float *pz) ;
static int mrisCalculateCanonicalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                     float *px, float *py, float *pz) ;
static int mrisDirectionTriangleIntersection(MRI_SURFACE *mris, float x0, 
                                             float y0, float z0, float nx, 
                                             float ny, float nz, MHT *mht, 
                                             double *pdist) ;
static int mrisComputeCurvatureValues(MRI_SURFACE *mris) ;
static int 
mrisAllNormalDirectionCurrentTriangleIntersections(MRI_SURFACE *mris,
                                                   VERTEX *v, MHT *mht,
                                                          double *pdist,
                                                          int *flist);
static int  load_triangle_vertices(MRI_SURFACE *mris, int fno, double U0[3], 
                                   double U1[3], double U2[3]) ;
static int  load_orig_triangle_vertices(MRI_SURFACE *mris, int fno, 
                                        double U0[3], double U1[3], 
                                        double U2[3]) ;
static void    mrisDumpFace(MRI_SURFACE *mris, int fno, FILE *fp) ;
static int    mrisAddEdge(MRI_SURFACE *mris, int vno1, int vno2) ;

#if 0
static int mrisNormalDirectionTriangleIntersection(MRI_SURFACE*mris,VERTEX *v,
                                                   MHT *mht, double *pdist,
                                                   int *flist);
static int mrisAllCurrentTriangleIntersections(MRI_SURFACE *mris, float x, 
                                               float y, float z, float nx, 
                                               float ny, float nz, 
                                               MHT *mht, int *flist) ;
static double mrisFindClosestFilledVoxel(MRI_SURFACE *mris, MRI *mri_filled, 
                                         int vno, double max_dist) ;
static int   mrisCheck(MRI_SURFACE *mris) ;
static int   mrisClipGradient(MRI_SURFACE *mris, float max_len) ;
static int   mrisClipMomentumGradient(MRI_SURFACE *mris, float max_len) ;
#endif
static int   mrisComputeSurfaceDimensions(MRI_SURFACE *mris) ;
static int   mrisFindNeighbors(MRI_SURFACE *mris) ;
static void  mrisNormalize(float v[3]) ;
static float mrisTriangleArea(MRIS *mris, int fac, int n) ;
static int   mrisNormalFace(MRIS *mris, int fac,int n,float norm[]) ;
static int   mrisComputeOrigNormal(MRIS *mris, int vno, float norm[]) ;
static int   mrisOrigNormalFace(MRIS *mris, int fac,int n,float norm[]) ;
static int   mrisReadTransform(MRIS *mris, char *mris_fname) ;
static MRI_SURFACE *mrisReadAsciiFile(char *fname) ;
static MRI_SURFACE *mrisReadGeoFile(char *fname) ;
static int         mrisReadGeoFilePositions(MRI_SURFACE *mris,char *fname) ;
static MRI_SURFACE *mrisReadTriangleFile(char *fname, double pct_over) ;
static int         mrisReadTriangleFilePositions(MRI_SURFACE*mris,
                                                  char *fname) ;

static SMALL_SURFACE *mrisReadTriangleFileVertexPositionsOnly(char *fname) ;

/*static int   mrisReadFieldsign(MRI_SURFACE *mris, char *fname) ;*/
static double mrisComputeNonlinearAreaSSE(MRI_SURFACE *mris) ;
static double mrisComputeNonlinearDistanceSSE(MRI_SURFACE *mris) ;
static double mrisComputeSpringEnergy(MRI_SURFACE *mris) ;
static double mrisComputeThicknessSmoothnessEnergy(MRI_SURFACE *mris,
                                         double l_repulse) ;
static double mrisComputeRepulsiveEnergy(MRI_SURFACE *mris,double l_repulse, MHT *mht_v_current) ;
static int    mrisComputeRepulsiveTerm(MRI_SURFACE *mris,
                                         double l_repulse, MHT *mht_v) ;
static double mrisComputeRepulsiveRatioEnergy(MRI_SURFACE *mris,
                                              double l_repulse) ;
static int    mrisComputeRepulsiveRatioTerm(MRI_SURFACE *mris,
                                         double l_repulse, MHT *mht_v) ;
static int    mrisComputeSurfaceRepulsionTerm(MRI_SURFACE *mris, 
                                              double l_repulse, MHT *mht);
static int    mrisComputeThicknessSmoothnessTerm(MRI_SURFACE *mris,
                                         double l_tsmooth) ;
static double mrisComputeTangentialSpringEnergy(MRI_SURFACE *mris) ;
static double mrisComputeIntensityError(MRI_SURFACE *mris, 
                                        INTEGRATION_PARMS *parms);
#if 0
static int    mrisMarkSulcalVertices(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static int    mrisUpdateSulcalGradients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
#endif
static double mrisComputeIntensityGradientError(MRI_SURFACE *mris, 
                                        INTEGRATION_PARMS *parms);
static double mrisComputeSphereError(MRI_SURFACE *mris, 
                                     double l_sphere,double a);
static double mrisComputeDistanceError(MRI_SURFACE *mris) ;
static double mrisComputeCorrelationError(MRI_SURFACE *mris, 
                                          INTEGRATION_PARMS *parms, 
                                          int use_stds) ;
static int    mrisComputeVertexDistances(MRI_SURFACE *mris) ;
static double mrisComputeError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                                float *parea_rms, float *pangle_rms,
                               float *pcurv_rms, float *pdist_rms,
                               float *pcorr_rms);
static int   mrisAverageGradients(MRI_SURFACE *mris, int num_avgs) ;
static int   mrisIntegrationEpoch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                                  int n_avgs);
static int   mrisRemoveNegativeArea(MRI_SURFACE *mris,INTEGRATION_PARMS *parms,
                                    int n_avgs, float min_area_pct,
                                    int max_passes);
static double mrisLineMinimize(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisLineMinimizeSearch(MRI_SURFACE *mris, 
                                     INTEGRATION_PARMS *parms);
static double  mrisAsynchronousTimeStep(MRI_SURFACE *mris, float momentum, 
                                    float dt, MHT *mht, float max_mag) ;
static double  mrisAsynchronousTimeStepNew(MRI_SURFACE *mris, float momentum, 
                                    float dt, MHT *mht, float max_mag) ;
static double mrisAdaptiveTimeStep(MRI_SURFACE *mris,INTEGRATION_PARMS*parms);
static int   mrisOrientEllipsoid(MRI_SURFACE *mris) ;
static int   mrisOrientPlane(MRI_SURFACE *mris) ;
#if AVERAGE_AREAS
static int   mrisAverageAreas(MRI_SURFACE *mris, int num_avgs, int which) ;
#endif
static int   transform(float *xptr, float *yptr, float *zptr, 
                       float nx, float ny, float nz, float d) ;

static int   mrisComputeTangentPlanes(MRI_SURFACE *mris) ;
static int   mrisRemoveLink(MRI_SURFACE *mris, int vno1, int vno2) ;
static int   mrisRemoveEdge(MRI_SURFACE *mris, int vno1, int vno2) ;
static int   mrisRemoveFace(MRI_SURFACE *mris, int fno) ;
static int   mrisCountTotalNeighbors(MRI_SURFACE *mris) ;
static int   mrisCountValidLinks(MRI_SURFACE *mris, int vno1, int vno2) ;
static int   mrisComputeSpringTerm(MRI_SURFACE *mris, double l_spring);
static int   mrisComputeLinkTerm(MRI_SURFACE *mris, double l_spring, int pial);
static int   mrisComputeNormalizedSpringTerm(MRI_SURFACE *mris, 
                                             double l_spring);
static int   mrisComputeIntensityTerm(MRI_SURFACE*mris,double l_intensity,
                                      MRI *mri_brain, MRI *mri_smooth,
                                      double sigma);
static int   mrisComputeIntensityGradientTerm(MRI_SURFACE*mris,
                                              double l_grad,
                                              MRI *mri_brain, MRI *mri_smooth);
static int   mrisComputeSphereTerm(MRI_SURFACE *mris, double l_sphere, 
                                   float radius) ;
static int   mrisComputeConvexityTerm(MRI_SURFACE *mris, double l_convex) ;
static int   mrisComputeExpansionTerm(MRI_SURFACE *mris, double l_expand) ;
static int   mrisComputeDistanceTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeNonlinearDistanceTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeCorrelationTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, 
                                              double l_curv) ;
static double  mrisComputeQuadraticCurvatureSSE(MRI_SURFACE *mris, 
                                              double l_curv) ;
static int   mrisComputePolarCorrelationTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeAngleAreaTerms(MRI_SURFACE *mris, 
                                       INTEGRATION_PARMS *parms) ;
static int   mrisComputeNonlinearAreaTerm(MRI_SURFACE *mris, 
                                       INTEGRATION_PARMS *parms) ;
static int   mrisClearDistances(MRI_SURFACE *mris) ;
static int   mrisClearExtraGradient(MRI_SURFACE *mris) ;
static int   mrisClearMomentum(MRI_SURFACE *mris) ;
static int   mrisValidVertices(MRI_SURFACE *mris) ;
static int   mrisValidFaces(MRI_SURFACE *mris) ;
static int   mrisLabelVertices(MRI_SURFACE *mris, float cx, float cy, 
                               float cz, int label, float radius) ;
static int mrisComputeShrinkwrapTerm(MRI_SURFACE *mris, MRI *mri_brain, double  l_shrinkwrap) ;
static double mrisComputeShrinkwrapError(MRI_SURFACE *mris, MRI *mri_brain, double l_shrinkwrap) ;


#if 0
static double mrisFindNormalDistance(MRI_SURFACE *mris, MHT *mht, int vno, 
                                      double max_dist);
static int    mrisFindNextOutwardFace(MRI_SURFACE *mris, MHT *mht, int vno, 
                                      double max_dist);
static int    mrisFindNextInwardFace(MRI_SURFACE *mris, MHT *mht, int vno, 
                                      double max_dist);
#endif

static int mrisProjectSurface(MRI_SURFACE *mris) ;
static int mrisOrientSurface(MRI_SURFACE *mris) ;
static int   mrisComputeBoundaryNormals(MRI_SURFACE *mris) ;
static int   mrisSmoothBoundaryNormals(MRI_SURFACE *mris, int niter) ;
static int   mrisFlipPatch(MRI_SURFACE *mris) ;

static int    mrisPlaceVertexInOrigFace(MRI_SURFACE *mris, VERTEX *v,int fno);

#if 0
static int    vertexInFace(MRI_SURFACE *mris, int vno, int fno)  ;
/* not currently used */
static int  mrisNeighborAtVoxel(MRI_SURFACE *mris, MRI *mri, int vno, 
                                int xv,int yv,int zv) ;
static int   mrisAverageDs(MRI_SURFACE *mris, int num_avgs) ;
static int mrisComputeAverageNormalTerm(MRI_SURFACE *mris, int navgs, 
                                        double l_normal) ;
static int  mrisFindNormalDistanceLimits(MRI_SURFACE *mris, MRI *mri_filled, 
                                         int vno, float max_dist,
                                         float *pmax_outward_distance,
                                         float *pmax_inward_distance) ;
#endif
  
static int  mrisComputeTangentialSpringTerm(MRI_SURFACE *mris,double l_spring);
static int   mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring) ;
static int   mrisRemoveNeighborGradientComponent(MRI_SURFACE *mris, int vno) ;
static int   mrisRemoveNormalGradientComponent(MRI_SURFACE *mris, int vno) ;

#if 0
static int   mrisSmoothNormalOutliers(MRI_SURFACE *mris, double nlen) ;
static int   mrisDebugVertex(MRI_SURFACE *mris, int vno) ;
static int    mrisComputeBoundaryTerm(MRI_SURFACE *mris, 
                                      INTEGRATION_PARMS *parms) ;
static int   mrisComputeCurvatureTerm(MRI_SURFACE *mris, 
                                      INTEGRATION_PARMS *parms) ;
static int   mrisComputeNegTerm(MRI_SURFACE *mris,INTEGRATION_PARMS *parms);
static int   mrisCountNegativeVertices(MRI_SURFACE *mris) ;
static double mrisComputeAverageHeight(MRI_SURFACE *mris) ;


static int   mrisComputeSethianCurvatureTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisSmoothCurvatures(MRI_SURFACE *mris, int niterations) ;
static int   mrisSmoothNormals(MRI_SURFACE *mris, int niterations) ;
static int   mrisComputeCurvatureGradientTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisStoreCurrentGradient(MRI_SURFACE *mris) ;
static int   mrisFindPoles(MRIS *mris) ;
static int   mrisComputeEllipsoidProperties(MRI_SURFACE *mris) ;
#endif
static int   mrisLogIntegrationParms(FILE *fp, MRI_SURFACE *mris,
                                     INTEGRATION_PARMS *parms) ;
static int   mrisLogStatus(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                           FILE *fp, float dt) ;
static int   mrisWriteSnapshots(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                                int t) ;
static int   mrisWriteSnapshot(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                               int t) ;
static int   mrisTrackTotalDistance(MRI_SURFACE *mris) ;
static int   mrisTrackTotalDistanceNew(MRI_SURFACE *mris) ;
static int  mrisLimitGradientDistance(MRI_SURFACE *mris, MHT *mht, int vno) ;
static int mrisFillFace(MRI_SURFACE *mris, MRI *mri, int fno) ;
static int mrisHatchFace(MRI_SURFACE *mris, MRI *mri, int fno, int on) ;
#if 0
static int mrisEraseFace(MRI_SURFACE *mris, MRI *mri, int fno) ;
static int  mrisRipVertices(MRI_SURFACE *mris) ;
#endif
static double mrisRmsValError(MRI_SURFACE *mris, MRI *mri) ;
static int mrisRemoveVertexLink(MRI_SURFACE *mris, int vno1, int vno2) ;
static int mrisStoreVtotalInV3num(MRI_SURFACE *mris) ;
static int  mrisFindAllOverlappingFaces(MRI_SURFACE *mris, MHT *mht,int fno, 
                                        int *flist) ;
#if 0
static int   mrisAddVertices(MRI_SURFACE *mris, double max_len) ;
#endif
static int mrisDivideEdge(MRI_SURFACE *mris, int vno1, int vno2) ;
static int mrisDivideFace(MRI_SURFACE *mris, int fno, int vno1, int vno2, 
                          int vnew_no) ;

/*--------------------------------------------------------------------*/

/*--------------------- CONSTANTS AND MACROS -------------------------*/

#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_QUAD_FILE_MAGIC_NUMBER  (-3 & 0x00ffffff)
/* 16777215 = 0xFFFFFF */
#define NEW_VERSION_MAGIC_NUMBER  16777215
#define START_Y                   (-128)
#define SLICE_THICKNESS           1

#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y,\
                                               v1->z-v0->z)
#define VERTEX_ORIG_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->origx-v0->origx,v1->origy-v0->origy,\
                                               v1->origz-v0->origz)

/* 77557 is horizontal with positive area */
/* 126906 has dy = 0 with negative area */
/* 102961 has a = 0 */
/* 77115 is > pi/2, face 1 */
/* v 82875, f 0 has negative area and horizontal */
/* v 115365, f 1 has negative area and is vertical */
/* v 75530, f 4 has negative area and is vertical */
#if 0
#define DEBUG_FACE(vno, fno)   (((fno) == 0) && (Gdiag & DIAG_SURFACE) &&\
#define DEBUG_FACE(vno, fno)   (((fno) == 2) && (Gdiag & DIAG_SURFACE) &&\
                                (vno == 79881))
#endif
#define DEBUG_FACE(vno, fno)   (((fno) == 4) && (Gdiag & DIAG_SURFACE) &&\
                                (DEBUG_VERTEX(vno)))
#define VDEBUG_FACE(fno)   (DEBUG_FACE(fno) && 0)
#define DEBUG_VERTEX(v)   (((v) == 75530) && (Gdiag & DIAG_SURFACE) && 1)
#define VDEBUG_VERTEX(v)   (((v) == 77115) && (Gdiag & DIAG_SURFACE) && 0)

/*--------------------------------------------------------------------*/

/*-------------------------- STATIC DATA -------------------------------*/

/*-------------------------- FUNCTIONS -------------------------------*/
double (*gMRISexternalGradient)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL ;
double (*gMRISexternalSSE)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL ;
double (*gMRISexternalRMS)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL ;
int (*gMRISexternalTimestep)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL ;
int (*gMRISexternalRipVertices)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)=NULL;
int (*gMRISexternalClearSSEStatus)(MRI_SURFACE *mris) = NULL ;
int (*gMRISexternalReduceSSEIncreasedGradients)(MRI_SURFACE *mris, double pct) = NULL ;

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISreadOverAlloc(char *fname, double pct_over)
{
  MRI_SURFACE *mris = NULL ;
  int         nquads, nvertices, magic, version, ix, iy, iz, vno, fno, n, m,
              imnr, imnr0, imnr1, type, vertices[VERTICES_PER_FACE+1], num ;
  float       x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  FILE        *fp = NULL ;
  VERTEX      *vertex ;
  FACE        *face ;
  int         tag;

  chklc() ;    /* check to make sure license.dat is present */
  type = MRISfileNameType(fname) ; /* using extension to get type */
  if (type == MRIS_ASCII_TRIANGLE_FILE)  /* .ASC */
  {
    mris = mrisReadAsciiFile(fname) ;
    if (!mris)
      return(NULL) ;
    version = -3 ;
  }
  else if (type == MRIS_ICO_FILE)        /* .TRI, .ICO */
  {
    mris = ICOreadOverAlloc(fname, pct_over) ;
    if (!mris)
      return(NULL) ;
    return(mris) ;
    version = -2 ;
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE) /* .GEO */ 
  {
    mris = mrisReadGeoFile(fname) ;
    if (!mris)
      return(NULL) ;
    version = -4;
  }
  else // default type MRIS_BINARY_QUADRANGLE_FILE ... use magic number
  {
    fp = fopen(fname, "rb") ;
    if (!fp)
      ErrorReturn(NULL,(ERROR_NOFILE,"MRISread(%s): could not open file",
                        fname));

    fread3(&magic, fp) ;
    if (magic == QUAD_FILE_MAGIC_NUMBER) 
    {
      version = -1;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stdout, "new surface file format\n");
    }
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) 
    {
      version = -2 ;
    }
    else if (magic == TRIANGLE_FILE_MAGIC_NUMBER)
    {
      fclose(fp) ;
      mris = mrisReadTriangleFile(fname, pct_over) ;
      if (!mris)
        ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n")) ;
      version = -3 ;
    }
    else /* no magic number assigned */
    {
      rewind(fp);
      version = 0;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        printf("surfer: old surface file format\n");
    }
  }
  /* some type of quadrangle file processing */
  if (version >= -2) 
  {
    fread3(&nvertices, fp);
    fread3(&nquads, fp);   /* # of qaudrangles - not triangles */
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout,"reading %d vertices and %d faces.\n",
      nvertices,2*nquads);
    
    mris = MRISoverAlloc(pct_over*nvertices,pct_over*2*nquads,nvertices, 
       2*nquads);
    mris->type = MRIS_BINARY_QUADRANGLE_FILE ;

    imnr0 = 1000 ;
    imnr1 = 0 ;
    /* read vertices *************************************************/
    for (vno = 0 ; vno < nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (version == -1)        /* QUAD_FILE_MAGIC_NUMBER */
      {
        fread2(&ix,fp);
        fread2(&iy,fp);
        fread2(&iz,fp);
        vertex->x = ix/100.0;
        vertex->y = iy/100.0;
        vertex->z = iz/100.0;
      }
      else  /* version == -2 */ /* NEW_QUAD_FILE_MAGIC_NUMBER */
      {
        vertex->x = freadFloat(fp) ;
        vertex->y = freadFloat(fp) ;
        vertex->z = freadFloat(fp) ;
      }
#if 0
      vertex->label = NO_LABEL ;
#endif
      /* brain-dead code and never used again either */
      imnr = (int)((vertex->y-START_Y)/SLICE_THICKNESS+0.5);
      if (imnr > imnr1)
        imnr1 = imnr ;
      if (imnr < imnr0)
        imnr0 = imnr ;
      if (version == 0)  /* old surface format */
      {
        fread1(&num,fp);   /* # of faces we are part of */
        vertex->num = num ;
        vertex->f = (int *)calloc(vertex->num,sizeof(int));
        if (!vertex->f)
          ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                    vertex->num) ;
        vertex->n = (uchar *)calloc(vertex->num,sizeof(uchar));
        if (!vertex->n)
          ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d nbrs",
                    vertex->n) ;
        for (n=0;n<vertex->num;n++)
          fread3(&vertex->f[n],fp);
      } 
      else 
        vertex->num = 0;   /* will figure it out */
    }
    /* read face vertices *******************************************/
    for (fno = 0 ; fno < mris->nfaces ; fno += 2)
    {
      int which ;

      if (fno == 86)
        DiagBreak() ;
      for (n = 0 ; n < 4 ; n++)   /* read quandrangular face */
      {
        fread3(&vertices[n],fp);
        if (vertices[n] == 22)
          DiagBreak() ;
      }

/* if we're going to be arbitrary, we might as well be really arbitrary */
#define WHICH_FACE_SPLIT(vno0, vno1) \
            (1*nint(sqrt(1.9*vno0) + sqrt(3.5*vno1)))
      /* 
         NOTE: for this to work properly in the write, the first two
         vertices in the first face (EVEN and ODD) must be 0 and 1.
         */
      which = WHICH_FACE_SPLIT(vertices[0], vertices[1]) ;

      /* 1st triangle */
      if (EVEN(which))
      {
        mris->faces[fno].v[0] = vertices[0] ;
        mris->faces[fno].v[1] = vertices[1] ;
        mris->faces[fno].v[2] = vertices[3] ;
        
        /* 2nd triangle */
        mris->faces[fno+1].v[0] = vertices[2] ;
        mris->faces[fno+1].v[1] = vertices[3] ;
        mris->faces[fno+1].v[2] = vertices[1] ;
      }
      else
      {
        mris->faces[fno].v[0] = vertices[0] ;
        mris->faces[fno].v[1] = vertices[1] ;
        mris->faces[fno].v[2] = vertices[2] ;
        
        /* 2nd triangle */
        mris->faces[fno+1].v[0] = vertices[0] ;
        mris->faces[fno+1].v[1] = vertices[2] ;
        mris->faces[fno+1].v[2] = vertices[3] ;
      }
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      {
        mris->vertices[mris->faces[fno].v[n]].num++;
        mris->vertices[mris->faces[fno+1].v[n]].num++;
      }
    }
    // new addition
    if (freadIntEx(&tag, fp))
    {
      if (tag == TAG_USEREALRAS)
	if (!freadIntEx(&mris->useRealRAS,fp)) // set useRealRAS
	  mris->useRealRAS = 0; // if error, set to default
    }
    else // no tag found.
    {
      // mark vertex coordinates are using the conformed (256^3) and c_(r,a,s) = 0.
      mris->useRealRAS = 0;
    }
    fclose(fp);
  }
  /* end of quadrangle file processing */
  /* file is closed now for all types ***********************************/

  /* find out if this surface is lh or rh from fname */
  strcpy(mris->fname, fname) ;
  {
    char *surf_name ;

    surf_name = strrchr(fname, '/') ;
    if (surf_name == NULL)
      surf_name = fname ;
    else
      surf_name++ ;  /* past the last slash */
    if (toupper(*surf_name) == 'R')
      mris->hemisphere = RIGHT_HEMISPHERE ;
    else
      mris->hemisphere = LEFT_HEMISPHERE ;
  }

  /***********************************************************************/
  /* build members of mris structure                                     */
  /***********************************************************************/
  if ((version<0) || type == MRIS_ASCII_TRIANGLE_FILE)
  {
    for (vno = 0 ; vno< mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      mris->vertices[vno].f = 
        (int *)calloc(mris->vertices[vno].num,sizeof(int));
      if (!mris->vertices[vno].f)
        ErrorExit(ERROR_NOMEMORY, 
                  "MRISread(%s): could not allocate %d faces at %dth vertex",
                  fname, vno, mris->vertices[vno].num) ;
      
      mris->vertices[vno].n = 
        (uchar *)calloc(mris->vertices[vno].num,sizeof(uchar));
      if (!mris->vertices[vno].n)
        ErrorExit(ERROR_NOMEMORY, 
                  "MRISread(%s): could not allocate %d indices at %dth vertex",
                  fname, vno, mris->vertices[vno].num) ;
      mris->vertices[vno].num = 0 ;
    }
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        mris->vertices[face->v[n]].f[mris->vertices[face->v[n]].num++] = fno;
    }
  }


  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    mris->vertices[vno].curv = 0;
    mris->vertices[vno].origarea = -1;
    mris->vertices[vno].border = 0;
#if 0
    mris->vertices[vno].origripflag = 0;
    mris->vertices[vno].ripflag = 0;
    mris->vertices[vno].val = 0;
    mris->vertices[vno].dist = 0;
    mris->vertices[vno].mx = 0;
    mris->vertices[vno].my = 0;
    mris->vertices[vno].mz = 0;
    mris->vertices[vno].fieldsign = 0;
    mris->vertices[vno].fsmask = 1;
    mris->vertices[vno].nc = 0;
    mris->vertices[vno].marked = 0;
#endif
    for (n=0;n<mris->vertices[vno].num;n++)
    {
      for (m=0;m<VERTICES_PER_FACE;m++)
      {
        if (mris->faces[mris->vertices[vno].f[n]].v[m] == vno)
          mris->vertices[vno].n[n] = m;
      }
    }
    x = mris->vertices[vno].x;
    y = mris->vertices[vno].y;
    z = mris->vertices[vno].z;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  mris->xlo = xlo ; mris->ylo = ylo ; mris->zlo = zlo ;
  mris->xhi = xhi ; mris->yhi = yhi ; mris->zhi = zhi ;
  mris->xctr = (xhi+xlo)/2;
  mris->yctr = (yhi+ylo)/2;
  mris->zctr = (zhi+zlo)/2;
  mrisFindNeighbors(mris);
  MRIScomputeNormals(mris);
  mrisComputeVertexDistances(mris) ;
  mrisReadTransform(mris, fname) ;
  if (type == MRIS_ASCII_TRIANGLE_FILE || type == MRIS_GEO_TRIANGLE_FILE)
  {
#if 0
    MRISsetNeighborhoodSize(mris, 2) ;
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
#endif
  }
  else 
  {
#if 0
    if (MRISreadBinaryCurvature(mris, fname) != NO_ERROR)
    {
      fprintf(stdout, "computing surface curvature directly...\n") ;
      MRISsetNeighborhoodSize(mris, 2) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
    }
       
    if (MRISreadBinaryAreas(mris, fname) != NO_ERROR)
     fprintf(stdout, "ignoring area file...\n") ; /*return(NULL) ;*/
#endif
  }

  mris->radius = MRISaverageRadius(mris) ;
#if 0
  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
#endif
  MRIScomputeMetricProperties(mris) ;
  /*  mrisFindPoles(mris) ;*/

  MRISstoreCurrentPositions(mris) ;
  return(mris) ;
}
/*-----------------------------------------------------

  MRISfastRead() just calls MRISRead()

        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISfastRead(char *fname)
{
  /********* why you keep the rest ? ******************/
#if 1
  return(MRISread(fname)) ;
#else
  MRI_SURFACE *mris ;
  int         nquads, nvertices, magic, version, ix, iy, iz, vno, fno, n, m,
              imnr, imnr0, imnr1, type, vertices[4], num ;
  float       x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  FILE        *fp ;
  VERTEX      *vertex ;
  FACE        *face ;

  mris = NULL ; fp = NULL ;
  type = MRISfileNameType(fname) ;
  if (type == MRIS_ASCII_TRIANGLE_FILE)
  {
    mris = mrisReadAsciiFile(fname) ;
    if (!mris)
      return(NULL) ;
    version = -3 ;
  }
  else if (type == MRIS_ICO_FILE)
  {
    mris = ICOread(fname) ;
    if (!mris)
      return(NULL) ;
    return(mris) ;
    version = -2 ;
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE)
  {
    mris = mrisReadGeoFile(fname) ;
    if (!mris)
      return(NULL) ;
    version = -4 ;
  }
  else   /* custom binary file - find out which type using magic # */
  {
    fp = fopen(fname, "rb") ;
    if (!fp)
      ErrorReturn(NULL,(ERROR_NOFILE,"MRISread(%s): could not open file",
                        fname));

    fread3(&magic, fp) ;
    if (magic == TRIANGLE_FILE_MAGIC_NUMBER)
    {
      fclose(fp) ;
      mris = mrisReadTriangleFile(fname, 0.0) ;
      if (!mris)
        ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n")) ;
      version = -3 ;
    }
    else if (magic == QUAD_FILE_MAGIC_NUMBER) 
    {
      version = -1;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stdout, "new surface file format\n");
    }
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) 
    {
      version = -2;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stdout, "new surface file format\n");
    }
    else 
    {
      rewind(fp);
      version = 0;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        printf("surfer: old surface file format\n");
    }
  }
  if (version >= -2)
  {
    fread3(&nvertices, fp);
    fread3(&nquads, fp);
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout,
              "reading %d vertices and %d faces.\n",nvertices,2*nquads);
    
    mris = MRISalloc(nvertices, 2*nquads) ;
    mris->type = MRIS_BINARY_QUADRANGLE_FILE ;
    
    imnr0 = 1000 ;
    imnr1 = 0 ;
    for (vno = 0 ; vno < nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (version == -1)
      {
        fread2(&ix,fp);
        fread2(&iy,fp);
        fread2(&iz,fp);
        vertex->x = ix/100.0;
        vertex->y = iy/100.0;
        vertex->z = iz/100.0;
      }
      else  /* version == -2 */
      {
        vertex->x = freadFloat(fp) ;
        vertex->y = freadFloat(fp) ;
        vertex->z = freadFloat(fp) ;
      }
#if 0
      vertex->label = NO_LABEL ;
#endif
      imnr = (int)((vertex->y-START_Y)/SLICE_THICKNESS+0.5);
      if (imnr > imnr1)
        imnr1 = imnr ;
      if (imnr < imnr0)
        imnr0 = imnr ;
      if (version == 0)  /* old surface format */
      {
        fread1(&num,fp);
        vertex->num = num ;
        vertex->f = (int *)calloc(vertex->num,sizeof(int));
        if (!vertex->f)
          ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                    vertex->num) ;
        vertex->n = (uchar *)calloc(vertex->num,sizeof(uchar));
        if (!vertex->n)
          ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d nbrs",
                    vertex->n) ;
        for (n=0;n<vertex->num;n++)
          fread3(&vertex->f[n],fp);
      } else vertex->num = 0;
    }
    
    for (fno = 0 ; fno < mris->nfaces ; fno += 2)
    {
      for (n = 0 ; n < 4 ; n++)   /* read quandrangular face */
        fread3(&vertices[n],fp);

      /* 1st triangle */
      mris->faces[fno].v[0] = vertices[0] ;
      mris->faces[fno].v[1] = vertices[1] ;
      mris->faces[fno].v[2] = vertices[3] ;
      if (version < 0)
        for (n = 0 ; n < VERTICES_PER_FACE ; n++)
          mris->vertices[mris->faces[fno].v[n]].num++;

      /* 2nd triangle */
      mris->faces[fno+1].v[0] = vertices[2] ;
      mris->faces[fno+1].v[1] = vertices[3] ;
      mris->faces[fno+1].v[2] = vertices[1] ;
      if (version < 0)
        for (n = 0 ; n < VERTICES_PER_FACE ; n++)
          mris->vertices[mris->faces[fno+1].v[n]].num++;
    }
    fclose(fp);
  }
  strcpy(mris->fname, fname) ;
  {
    char *surf_name ;

    surf_name = strrchr(fname, '/') ;
    if (surf_name == NULL)
      surf_name = fname ;
    else
      surf_name++ ;  /* past the last slash */
    if (toupper(*surf_name) == 'R')
      mris->hemisphere = RIGHT_HEMISPHERE ;
    else
      mris->hemisphere = LEFT_HEMISPHERE ;
  }
  if ((version<0) || type == MRIS_ASCII_TRIANGLE_FILE)
  {
    for (vno = 0 ; vno< mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      mris->vertices[vno].f = 
        (int *)calloc(mris->vertices[vno].num,sizeof(int));
      if (!mris->vertices[vno].f)
        ErrorExit(ERROR_NOMEMORY, 
                  "MRISread(%s): could not allocate %d faces at %dth vertex",
                  fname, vno, mris->vertices[vno].num) ;
      
      mris->vertices[vno].n = 
        (uchar *)calloc(mris->vertices[vno].num,sizeof(uchar));
      if (!mris->vertices[vno].n)
        ErrorExit(ERROR_NOMEMORY, 
                  "MRISread(%s): could not allocate %d indices at %dth vertex",
                  fname, vno, mris->vertices[vno].num) ;
      mris->vertices[vno].num = 0 ;
    }
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      {
        vertex = &mris->vertices[face->v[n]] ;
        vertex->f[vertex->num++] = fno;
      }
    }
  }


  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    mris->vertices[vno].curv = 0;
    mris->vertices[vno].origarea = -1;
    mris->vertices[vno].border = 0;
#if 0
    mris->vertices[vno].origripflag = 0;
    mris->vertices[vno].ripflag = 0;
    mris->vertices[vno].val = 0;
    mris->vertices[vno].dist = 0;
    mris->vertices[vno].mx = 0;
    mris->vertices[vno].my = 0;
    mris->vertices[vno].mz = 0;
    mris->vertices[vno].fieldsign = 0;
    mris->vertices[vno].fsmask = 1;
    mris->vertices[vno].nc = 0;
    mris->vertices[vno].marked = 0;
#endif
    for (n=0;n<mris->vertices[vno].num;n++)
    {
      for (m=0;m<VERTICES_PER_FACE;m++)
      {
        if (mris->faces[mris->vertices[vno].f[n]].v[m] == vno)
          mris->vertices[vno].n[n] = m;
      }
    }
    x = mris->vertices[vno].x;
    y = mris->vertices[vno].y;
    z = mris->vertices[vno].z;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  mris->xlo = xlo ; mris->ylo = ylo ; mris->zlo = zlo ;
  mris->xhi = xhi ; mris->yhi = yhi ; mris->zhi = zhi ;
  mris->xctr = (xhi+xlo)/2;
  mris->yctr = (yhi+ylo)/2;
  mris->zctr = (zhi+zlo)/2;
  mrisFindNeighbors(mris);
  MRIScomputeNormals(mris);
#if 0
  mrisComputeVertexDistances(mris) ;
  mrisReadTransform(mris, fname) ;
#endif
  if (type == MRIS_ASCII_TRIANGLE_FILE || type == MRIS_GEO_TRIANGLE_FILE)
  {
    MRISsetNeighborhoodSize(mris, 2) ;
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
  }
  else 
  {
    if (MRISreadBinaryCurvature(mris, fname) != NO_ERROR)
      fprintf(stdout, "ignoring curvature file...\n") ; /*return(NULL) ;*/
#if 0
    if (MRISreadBinaryAreas(mris, fname) != NO_ERROR)
      return(NULL) ;
#endif
  }

#if 0
  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
#endif
  MRISstoreCurrentPositions(mris) ;
  return(mris) ;
#endif
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISread(char *fname)
{
  MRI_SURFACE  *mris ;

  mris = MRISreadOverAlloc(fname, 0.0) ;
  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwrite(MRI_SURFACE *mris, char *name)
{
  int   k, type ;
  float x,y,z;
  FILE  *fp;

  char  fname[STRLEN] ;

  chklc() ;
  MRISbuildFileName(mris, name, fname) ;
  type = MRISfileNameType(fname) ;
  if (type == MRIS_ASCII_TRIANGLE_FILE)
    return(MRISwriteAscii(mris, fname)) ;
  else if (type == MRIS_VTK_FILE)
    return(MRISwriteVTK(mris, fname)) ;
  else if (type == MRIS_GEO_TRIANGLE_FILE)
    return(MRISwriteGeo(mris, fname)) ;
  else if (type == MRIS_ICO_FILE)
    return MRISwriteICO(mris, fname);

  if (mris->type == MRIS_TRIANGULAR_SURFACE)
    return(MRISwriteTriangularSurface(mris, fname)) ;

  fp = fopen(fname,"w");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE,"MRISwrite(%s): can't create file\n",fname));
#if 0
  fwrite3(NEW_QUAD_FILE_MAGIC_NUMBER,fp);
#else
  fwrite3(QUAD_FILE_MAGIC_NUMBER,fp);
#endif
  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces/2,fp);   /* # of quadrangles */
  for (k = 0 ; k < mris->nvertices ; k++)
  {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
#if 1
    fwrite2((int)(x*100),fp);
    fwrite2((int)(y*100),fp);
    fwrite2((int)(z*100),fp);
#else
    fwriteFloat(x, fp) ;
    fwriteFloat(y, fp) ;
    fwriteFloat(z, fp) ;
#endif
  }
  for (k = 0 ; k < mris->nfaces ; k+=2)
  {
    int which ;
    FACE *f ;

    f = &mris->faces[k] ;
    {
      int n ;
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      {
        if ((mris->faces[k].v[n] == 22) || (mris->faces[k+1].v[n] == 22))
          DiagBreak() ;
      }
    }
    which = WHICH_FACE_SPLIT(f->v[0], f->v[1]) ;
    if (EVEN(which))
    {
      fwrite3(mris->faces[k].v[0],fp);
      fwrite3(mris->faces[k].v[1],fp);
      fwrite3(mris->faces[k+1].v[0],fp);
      fwrite3(mris->faces[k].v[2],fp);
    }
    else
    {
      fwrite3(mris->faces[k].v[0],fp);
      fwrite3(mris->faces[k].v[1],fp);
      fwrite3(mris->faces[k].v[2],fp);
      fwrite3(mris->faces[k+1].v[2],fp);
    }
  }
  /* write whether vertex data was using the real RAS rather than conformed RAS */
  fwriteInt(TAG_USEREALRAS, fp);
  fwriteInt(mris->useRealRAS, fp);
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE  *
MRISoverAlloc(int max_vertices, int max_faces, int nvertices, int nfaces)
{
  MRI_SURFACE  *mris ;

  if (max_vertices <= 0)
    max_vertices = nvertices ;
  if (max_faces <= 0)
    max_faces = nfaces ;
  mris = MRISalloc(max_vertices, max_faces) ;
  mris->nvertices = nvertices ;
  mris->nfaces = nfaces ;
  mris->max_vertices = max_vertices ;
  mris->max_faces = max_faces ;
  mris->useRealRAS = 0; /* just initialize */
  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISalloc(int nvertices, int nfaces)
{
  MRI_SURFACE   *mris ;

  mris = (MRI_SURFACE *)calloc(1, sizeof(MRI_SURFACE)) ;
  if (!mris)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISalloc(%d, %d): could not allocate mris structure",
              nvertices, nfaces);

  mris->nsize = 1 ;  /* only 1-connected neighbors initially */
  mris->nvertices = nvertices ;
  mris->nfaces = nfaces ;
  mris->vertices = (VERTEX *)calloc(nvertices, sizeof(VERTEX)) ;
  if (!mris->vertices)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISalloc(%d, %d): could not allocate vertices",
              nvertices, sizeof(VERTEX));
  mris->faces = (FACE *)calloc(nfaces, sizeof(FACE)) ;
  if (!mris->faces)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISalloc(%d, %d): could not allocate faces",
              nfaces, sizeof(FACE));
  mris->useRealRAS = 0; /* just initialize */
  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISfree(MRI_SURFACE **pmris)
{
  MRI_SURFACE  *mris ;
  int          vno ;

  mris = *pmris ;
  *pmris = NULL ;

  if (mris->dx2)
    free(mris->dx2) ;
  if (mris->dy2)
    free(mris->dy2) ;
  if (mris->dz2)
    free(mris->dz2) ;
  if (mris->labels)
    free(mris->labels) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (mris->vertices[vno].f)
      free(mris->vertices[vno].f) ;
    if (mris->vertices[vno].n)
      free(mris->vertices[vno].n) ;
    if (mris->vertices[vno].dist)
      free(mris->vertices[vno].dist) ;
    if (mris->vertices[vno].dist_orig)
      free(mris->vertices[vno].dist_orig) ;
    if (mris->vertices[vno].v)
      free(mris->vertices[vno].v) ;
  }

  if (mris->vertices)
    free(mris->vertices) ;
  if (mris->faces)
    free(mris->faces) ;
  free(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISfreeDists(MRI_SURFACE *mris)
{
  int          vno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (mris->vertices[vno].dist)
      free(mris->vertices[vno].dist) ;
    if (mris->vertices[vno].dist_orig)
      free(mris->vertices[vno].dist_orig) ;
    mris->vertices[vno].dist = mris->vertices[vno].dist_orig = NULL ;
    mris->vertices[vno].vtotal = 0 ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          if nsize <=0 then neighborhood size gets reset back to whatever
          it's max was.
------------------------------------------------------*/
int
MRISresetNeighborhoodSize(MRI_SURFACE *mris, int nsize)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    switch (nsize)
    {
    default:   /* reset back to original */
      switch (v->nsize)
      {
      default:
      case 1:
        v->vtotal = v->vnum ;
        break ;
      case 2:
        v->vtotal = v->v2num ;
        break ;
      case 3:
        v->vtotal = v->v3num ;
        break ;
      }
      break ;
    case 1:
      v->vtotal = v->vnum ;
      break ;
    case 2:
      v->vtotal = v->v2num ;
      break ;
    case 3:
      v->vtotal = v->v3num ;
      break ;
    }
  }
  mris->nsize = nsize ;
  return(NO_ERROR) ;
}


#define MAX_4_NEIGHBORS     100
#define MAX_3_NEIGHBORS     70
#define MAX_2_NEIGHBORS     20
#define MAX_1_NEIGHBORS     8
#define MAX_NEIGHBORS       (400)
static int
mrisFindNeighbors(MRI_SURFACE *mris)
{
  int          n0,n1,i,k,m,n, vno, vtotal, ntotal, vtmp[MAX_NEIGHBORS] ;
  FACE         *f;
  VERTEX       *v ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "finding surface neighbors...") ;

  for (k=0;k<mris->nvertices;k++)
  {
    if (k == Gdiag_no)
      DiagBreak() ;
    v = &mris->vertices[k];
    v->vnum = 0;
    for (m=0;m<v->num;m++)
    {
      n = v->n[m];     /* # of this vertex in the mth face that it is in */
      f = &mris->faces[v->f[m]];  /* ptr to the mth face */
      /* index of vertex we are connected to */
      n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;   
      n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;   
      for (i=0;i<v->vnum && vtmp[i]!=f->v[n0];i++);
      if (i==v->vnum)
        vtmp[(int)v->vnum++] = f->v[n0];
      for (i=0;i<v->vnum && vtmp[i]!=f->v[n1];i++);
      if (i==v->vnum)
        vtmp[(int)v->vnum++] = f->v[n1];
    }
    if (mris->vertices[k].v)
      free(mris->vertices[k].v) ;
    mris->vertices[k].v = (int *)calloc(mris->vertices[k].vnum,sizeof(int));
    if (!mris->vertices[k].v)
      ErrorExit(ERROR_NOMEMORY, 
                "mrisFindNeighbors: could not allocate nbr array") ;

    v->vtotal = v->vnum ;
    v->nsize = 1 ;
    for (i=0;i<v->vnum;i++)
    {
      v->v[i] = vtmp[i];
    }

    if (v->dist)
      free(v->dist) ;
    if (v->dist_orig)
      free(v->dist_orig) ;
    
    v->dist = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist)
      ErrorExit(ERROR_NOMEMORY,
                "mrisFindNeighbors: could not allocate list of %d "
                "dists at v=%d", v->vnum, k) ;
    v->dist_orig = (float *)calloc(v->vnum, sizeof(float)) ;
      if (!v->dist_orig)
        ErrorExit(ERROR_NOMEMORY,
                  "mrisFindNeighbors: could not allocate list of %d "
                  "dists at v=%d", v->vnum, k) ;
/*
    if (v->num != v->vnum)
      printf("%d: num=%d vnum=%d\n",k,v->num,v->vnum);
*/
  }
  for (k=0;k<mris->nfaces;k++)
  {
    f = &mris->faces[k];
    for (m=0;m<VERTICES_PER_FACE;m++)
    {
      v = &mris->vertices[f->v[m]];
      for (i=0;i<v->num && k!=v->f[i];i++);
      if (i==v->num)   /* face has vertex, but vertex doesn't have face */
        ErrorExit(ERROR_BADPARM,
                  "%s: face[%d].v[%d] = %d, but face %d not in vertex %d "
                  "face list\n", mris->fname,k,m,f->v[m], k, f->v[m]);
    }
  }

  for (vno = ntotal = vtotal = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    vtotal += v->vtotal ;
    ntotal++ ;
  }

  mris->avg_nbrs = (float)vtotal / (float)ntotal ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int          
MRISsampleAtEachDistance(MRI_SURFACE *mris,int nbhd_size,int nbrs_per_distance)
{
  int  n, nbrs_array[MAX_NBHD_SIZE] ;

  if (!nbhd_size)
    return(NO_ERROR) ;

  if (Gdiag & (DIAG_HEARTBEAT | DIAG_SHOW))
    fprintf(stdout, "sampling long-range distances") ;
  for (n = 0 ; n <= nbhd_size ; n++)
    nbrs_array[n] = nbrs_per_distance ;
  return(MRISsampleDistances(mris, nbrs_array, nbhd_size)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Expand the list of neighbors of each vertex, reallocating
          the v->v array to hold the expanded list.
------------------------------------------------------*/
#define MAX_VERTICES  20000
#define MAX_V         1000  /* max for any one node, actually way too big */
#define TRIANGLE_DISTANCE_CORRECTION  1.09f /*1.1f*/ /*1.066f*/ /*1.12578*/ /* 1.13105f*/
 /*1.1501f  (1.1364f)*/
#define QUADRANGLE_DISTANCE_CORRECTION  ((1+sqrt(2)) / 2) /* 1.2071  */
int
MRISsampleDistances(MRI_SURFACE *mris, int *nbrs, int max_nbhd)
{
  int          i,n, vno, vnum, old_vnum, total_nbrs, max_possible,max_v,vtotal;
  VERTEX       *v, *vn, *vn2 ;
  int          *vnbrs, *vall, found, n2, vnbrs_num, vall_num, nbhd_size,done ;
  float        xd, yd, zd, min_dist, dist, dist_scale, old_dist[MAX_V], 
               old_v[MAX_V], min_angle, angle ;
  VECTOR       *v1, *v2 ;
  float         c[100] ;
  int           nc[100] ;
  
  memset(c, 0, 100*sizeof(float)) ;
  memset(nc, 0, 100*sizeof(int)) ;
  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;

  /* adjust for Manhattan distance */
  if (IS_QUADRANGULAR(mris))
    dist_scale = (1.0 + sqrt(2.0)) / 2.0f ;
  else
    dist_scale = TRIANGLE_DISTANCE_CORRECTION ;

  vnbrs = (int *)calloc(MAX_VERTICES, sizeof(int)) ;
  vall = (int *)calloc(MAX_VERTICES, sizeof(int)) ;
  vtotal = total_nbrs = 0 ;
  for (vtotal = max_possible = 0, n = 1 ; n <= max_nbhd ; n++)
  {
    max_possible += nbrs[n] ;
    if (n > mris->nsize)
      vtotal += nbrs[n] ;
  }
  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stdout, 
            "\nsampling %d dists/vertex (%2.1f at each dist) = %2.1fMB\n",
            vtotal, 
            (float)vtotal/((float)max_nbhd-(float)mris->nsize),
            (float)vtotal*mrisValidVertices(mris)*sizeof(float)*3.0f / 
            (1024.0f*1024.0f)) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((Gdiag & DIAG_HEARTBEAT) && (!(vno % (mris->nvertices/10))))
      fprintf(stdout, "%%%1.0f done\n", 
              100.0f*(float)vno / (float)mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak()  ;

    if (v->ripflag)
      continue ;

    /* small neighborhood is always fixed, don't overwrite them */
    vtotal = v->vtotal ;
    if (v->v3num > 0)
      v->vtotal = v->v3num ;
    else if (v->v2num > 0)
      v->vtotal = v->v2num ;
    else
      v->vtotal = v->vnum ;

    max_v = v->vtotal+max_possible ;
    if (vtotal < max_v)  /* won't fit in current allocation,reallocate stuff */
    {
      /* save and restore neighbor list */
      memmove(old_v, v->v, v->vtotal*sizeof(v->v[0])) ;
      free(v->v) ;
      v->v = (int *)calloc(max_v, sizeof(int)) ;
      if (!v->v)
        ErrorExit(ERROR_NO_MEMORY, 
                  "MRISsampleDistances: could not allocate list of %d "
                  "nbrs at v=%d", max_v, vno) ;
      memmove(v->v, old_v, v->vtotal*sizeof(v->v[0])) ;
      
      /* save and restore distance vector */
      memmove(old_dist, v->dist, v->vtotal*sizeof(v->dist[0])) ;
      free(v->dist) ;
      v->dist = (float *)calloc(max_v, sizeof(float)) ;
      if (!v->dist)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d", max_v, vno) ;
      memmove(v->dist, old_dist, v->vtotal*sizeof(v->dist[0])) ;

      /* save and restore original distance vector */
      memmove(old_dist, v->dist_orig, v->vtotal*sizeof(v->dist_orig[0])) ;
      free(v->dist_orig) ;
      v->dist_orig = (float *)calloc(max_v, sizeof(float)) ;
      if (!v->dist_orig)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d", max_v, vno) ;
      memmove(v->dist_orig, old_dist, v->vtotal*sizeof(v->dist_orig[0])) ;
    }

    /* 
       find all the neighbors at each extent (i.e. 1-neighbors, then 
       2-neighbors, etc..., marking their corrected edge-length distances
       as you go.
    */
    vall[0] = vno ; vall_num = 1 ; old_vnum = 0 ; 
    v->marked = 1 ;  /* a hack - it is a zero neighbor */
    for (nbhd_size = 1 ; vall_num < MAX_VERTICES && nbhd_size <= max_nbhd ; 
         nbhd_size++)
    {
      /* expand neighborhood outward by a ring of vertices */
      vnbrs_num = 0 ;  /* will count neighbors in this ring */
      vnum = vall_num ;
      for (found = 0, n = old_vnum; vall_num<MAX_VERTICES && n < vall_num; n++)
      {
        vn = &mris->vertices[vall[n]] ;
        if (vn->ripflag)
          continue ;

        /* search through vn's neighbors to find an unmarked vertex */
        for (n2 = 0 ; n2 < vn->vnum ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (vn2->ripflag || vn2->marked)
            continue ;

          /* found one, mark it and put it in the vall list */
          found++ ;
          vn2->marked = nbhd_size ;
          vall[vnum++] = vn->v[n2] ;
          if (nbrs[nbhd_size] > 0)  /* want to store this distance */
          {
            vnbrs[vnbrs_num++] = vn->v[n2] ;
          }
        }
      }  /* done with all neighbors at previous distance */

      /* found all neighbors at this extent - calculate distances */
      old_vnum = vall_num ;  /* old_vnum is index of 1st nbr at this distance*/
      vall_num += found ;    /* vall_num is total # of nbrs */
      for (n = old_vnum ; n < vall_num ; n++)
      {
        vn = &mris->vertices[vall[n]] ;
        if (vn->ripflag)
          continue ;
#define UNFOUND_DIST 1000000.0f
        for (min_dist = UNFOUND_DIST, n2 = 0 ; n2 < vn->vnum ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (vn2->ripflag)
            continue ;
          if (!vn2->marked || vn2->marked == nbhd_size)
            continue ;
          xd = vn2->x - vn->x ; yd = vn2->y - vn->y ; zd = vn2->z - vn->z ;
          dist = sqrt(xd*xd + yd*yd + zd*zd) ;
#if !MULTI_DIST_SCALING
          if (nbhd_size > 1)
            dist /= dist_scale ;
          if (vn2->d+dist < min_dist)
            min_dist = vn2->d+dist ;
#else
          dist = 
            (dist + vn2->d * distance_scale[vn2->marked]) / 
            distance_scale[vn2->marked+1] ;
          if (dist < min_dist)
            min_dist = dist ;
#endif
        }
        vn->d = min_dist  ;
        if (nbhd_size <= 2)
        {
          xd = vn->x - v->x ; yd = vn->y - v->y ; zd = vn->z - v->z ;
          dist = sqrt(xd*xd + yd*yd + zd*zd) ;
          vn->d = dist ;
        }
      }

      /* 
         now check each to see if a neighbor at the same 'distance'
         is actually closer than neighbors which are 'nearer' (i.e. maybe
         the path through a 3-neighbor is shorter than that through any
         of the 2-neighbors.
      */
      for (n = old_vnum ; n < vall_num ; n++)
      {
        vn = &mris->vertices[vall[n]] ;
        if (vn->ripflag)
          continue ;
        min_dist = vn->d ;
        for (n2 = 0 ; n2 < vn->vnum ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (!vn2->marked || vn2->marked != nbhd_size || vn2->ripflag)
            continue ;
          xd = vn2->x - vn->x ; yd = vn2->y - vn->y ; zd = vn2->z - vn->z ;
          dist = sqrt(xd*xd + yd*yd + zd*zd) ;
#if !MULTI_DIST_SCALING
          if (nbhd_size > 1)
            dist /= dist_scale ;
          if (vn2->d+dist < min_dist)
            min_dist = vn2->d + dist ;
#else
          dist = 
            (dist + vn2->d * distance_scale[vn2->marked]) /
            distance_scale[vn2->marked+1] ;
          if (dist < min_dist)
            min_dist = dist ;
#endif
        }
        vn->d = min_dist ;
        {
          xd = vn->x - v->x ; yd = vn->y - v->y ; zd = vn->z - v->z ;
          dist = sqrt(xd*xd + yd*yd + zd*zd) ;
          c[nbhd_size] += vn->d / dist ;
          nc[nbhd_size]++ ;
        }
      }

      /* if this set of neighbors are to be stored, sample from them */
      if (nbrs[nbhd_size] <= 0)
        continue ;

      /* make sure the points are not too close together */
      min_angle = 0.9*2.0*M_PI / (float)nbrs[nbhd_size] ;

      /*
        at this point, the vall list contains all the neighbors currently found
        at ALL distances, while the vnbrs list contains ONLY the 
        nbhd_size-neighbors.
        */
      if (found <= nbrs[nbhd_size])  /* just copy them all in */
      {
        for (n = 0 ; n < found ; n++, v->vtotal++)
        {
          v->v[v->vtotal] = vnbrs[n] ;
          v->dist_orig[v->vtotal] = mris->vertices[vnbrs[n]].d ;
        }
      }
      else                   /* randomly sample from them */
      {
        int vstart = v->vtotal ;
        for (n = 0 ; n < nbrs[nbhd_size] ; n++, v->vtotal++)
        {
          int j, niter = 0 ;
          do
          {
            do
            {
              i = nint(randomNumber(0.0, (double)found-1)) ;
            } while (vnbrs[i] < 0) ;
            /* 
               now check to make sure that the angle between this
               point and the others already selected is not too
               small to make sure the points are not bunched.
               */
            vn = &mris->vertices[vnbrs[i]] ;
            VECTOR_LOAD(v1, vn->x-v->x, vn->y-v->y, vn->z-v->z) ;
            done = 1 ;
            for (j = vstart ; done && j < v->vtotal ; j++)
            {
              vn2 = &mris->vertices[v->v[j]] ;
              VECTOR_LOAD(v2, vn2->x-v->x, vn2->y-v->y, vn2->z-v->z) ;
              angle = Vector3Angle(v1, v2) ;
              if (angle < min_angle)
                done = 0 ;
            }
            if (++niter > found)  /* couldn't find enough at this difference */
            {
              min_angle *= 0.75f ;  /* be more liberal */
              niter = 0 ;
            }
          } while (!done && !FZERO(min_angle)) ;
          vn = &mris->vertices[vnbrs[i]] ;
          v->v[v->vtotal] = vnbrs[i] ;
          v->dist_orig[v->vtotal] = vn->d ;
          if (FZERO(vn->d))
            DiagBreak() ;
          vnbrs[i] = -1 ;
        }
      }
    }


    if ((Gdiag_no == vno) && DIAG_VERBOSE_ON)
    {
      FILE  *fp ;
      char  fname[STRLEN] ;

      sprintf(fname, "v%d", vno) ;
      fp = fopen(fname, "w") ;
      fprintf(fp, "%d\n", vall_num) ;
      for (n = 0 ; n < vall_num ; n++)
        fprintf(fp, "%d\n", vall[n]) ;
      fclose(fp) ;

      sprintf(fname, "vn%d", vno) ;
      fp = fopen(fname, "w") ;
      fprintf(fp, "%d\n", v->vtotal) ;
      for (n = 0 ; n < v->vtotal ; n++)
        fprintf(fp, "%d\n", v->v[n]) ;
      fclose(fp) ;
      for (n = 0 ; n < mris->nvertices ; n++)
      {
        vn = &mris->vertices[n] ;
#if 0
        if (vn->ripflag)
          continue ;
#endif
        vn->curv = vn->d ;
      }
      sprintf(fname, "%s.dist", mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh");
      MRISwriteCurvature(mris, fname) ;
    }

    /* 
       done building arrays - allocate distance vectors and
       sample from the found neighbors list.
       */
    /* now unmark them all */
    for (n = 0 ; n < vall_num ; n++)
    {
      mris->vertices[vall[n]].marked = 0 ;
      mris->vertices[vall[n]].d = 0.0 ;
    }

    total_nbrs += v->vtotal ;
  }

  /* now fill in immediate neighborhood(Euclidean) distances */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak()  ;
    if (v->ripflag)
      continue ;
    if (v->v3num > 0)
      vtotal = v->v3num ;
    else if (v->v2num > 0)
      vtotal = v->v2num ;
    else
      vtotal = v->vnum ;
    for (n = 0 ; n < vtotal ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag)
        continue ;
      xd = v->x - vn->x ; yd = v->y - vn->y ; zd = v->z - vn->z ;
      v->dist_orig[n] = sqrt(xd*xd+yd*yd+zd*zd) ;
    }
  }

  /* check reasonableness of distances */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    for (n = 0 ; n < v->vtotal ; n++)
    {
      if (FZERO(v->dist_orig[n]))
        fprintf(stderr, "zero distance at v %d, n %d (vn = %d)\n",
                vno, n, v->v[n]) ;
    }
  }

  mris->avg_nbrs = (float)total_nbrs / (float)mrisValidVertices(mris) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "avg_nbrs = %2.1f\n", mris->avg_nbrs) ;

#if MULTI_DIST_SCALING
  if (Gdiag & DIAG_SHOW)
  {
    for (n = 0 ; n <= max_nbhd ; n++)
    {
      if (nc[n])
        c[n] /= (float)nc[n] ;
      fprintf(stdout, "c[%d] = %2.5f (%d samples)\n", n, c[n], nc[n]) ;
    }
    fprintf(stdout, "c[] = { ") ;
    for (n = 0 ; n <= max_nbhd ; n++)
    {
      fprintf(stdout, "%2.5f", c[n]) ;
      if (n < max_nbhd)
        fprintf(stdout, ", ") ;
    }
  }
#endif
  free(vnbrs) ;
  free(vall) ;
  VectorFree(&v1) ; VectorFree(&v2) ;
  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stdout, " done.\n") ;
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Expand the list of neighbors of each vertex, reallocating
          the v->v array to hold the expanded list.
------------------------------------------------------*/
#define MAX_VERTICES  20000
#define MAX_V         1000  /* max for any one node, actually way too big */
int
MRISsampleDistances(MRI_SURFACE *mris, int *nbrs, int max_nbhd)
{
  int          i,n, vno, vnum, old_vnum, total_nbrs, max_possible,max_v,vtotal;
  VERTEX       *v, *vn, *vn2 ;
  int          *vnbrs, *vall, found, n2, vnbrs_num, vall_num, nbhd_size,done ;
  float        xd, yd, zd, min_dist, dist, dist_scale, old_dist[MAX_V], 
               old_v[MAX_V], min_angle, angle ;
  VECTOR       *v1, *v2 ;
  
  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;

  dist_scale = (1.0 + sqrt(2.0)) / 2 ;  /* adjust for Manhattan distance */
  vnbrs = (int *)calloc(MAX_VERTICES, sizeof(int)) ;
  vall = (int *)calloc(MAX_VERTICES, sizeof(int)) ;
  total_nbrs = 0 ;
  for (vtotal = max_possible = 0, n = 1 ; n <= max_nbhd ; n++)
  {
    max_possible += nbrs[n] ;
    if (n > mris->nsize)
      vtotal += nbrs[n] ;
  }

  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stdout, 
            "\nsampling %d dists/vertex (%2.1f at each dist) = %2.1fMB\n",
            vtotal, 
            (float)vtotal/((float)max_nbhd-(float)mris->nsize),
            (float)vtotal*mris->nvertices*sizeof(float)*3.0f / 
            (1024.0f*1024.0f)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((Gdiag & DIAG_HEARTBEAT) && (!(vno % (mris->nvertices/25))))
      fprintf(stdout, " %%%1.0f", 100.0f*(float)vno / (float)mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak()  ;

    if (v->ripflag)
      continue ;

    /* small neighborhood is always fixed, don't overwrite them */
    vtotal = v->vtotal ;
    if (v->v3num > 0)
      v->vtotal = v->v3num ;
    else if (v->v2num > 0)
      v->vtotal = v->v2num ;
    else
      v->vtotal = v->vnum ;

    max_v = v->vtotal+max_possible ;
    if (vtotal < max_v) /* won't fit in current allocation,reallocate stuff */
    {
      /* save and restore neighbor list */
      memmove(old_v, v->v, v->vtotal*sizeof(v->v[0])) ;
      free(v->v) ;
      v->v = (int *)calloc(max_v, sizeof(int)) ;
      if (!v->v)
        ErrorExit(ERROR_NO_MEMORY, 
                  "MRISsampleDistances: could not allocate list of %d "
                  "nbrs at v=%d", max_v, vno) ;
      memmove(v->v, old_v, v->vtotal*sizeof(v->v[0])) ;
      
      /* save and restore distance vector */
      memmove(old_dist, v->dist, v->vtotal*sizeof(v->dist[0])) ;
      free(v->dist) ;
      v->dist = (float *)calloc(max_v, sizeof(float)) ;
      if (!v->dist)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d", max_v, vno) ;
      memmove(v->dist, old_dist, v->vtotal*sizeof(v->dist[0])) ;

      /* save and restore original distance vector */
      memmove(old_dist, v->dist_orig, v->vtotal*sizeof(v->dist_orig[0])) ;
      free(v->dist_orig) ;
      v->dist_orig = (float *)calloc(max_v, sizeof(float)) ;
      if (!v->dist_orig)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d", max_v, vno) ;
      memmove(v->dist_orig, old_dist, v->vtotal*sizeof(v->dist_orig[0])) ;
    }

    vall[0] = vno ; vall_num = 1 ; old_vnum = 0 ; 
    v->marked = 1 ;  /* a hack - it is a zero neighbor */
    for (nbhd_size = 1 ; vall_num < MAX_VERTICES && nbhd_size <= max_nbhd ; 
         nbhd_size++)
    {
      /* expand neighborhood outward by a ring of vertices */
      vnbrs_num = 0 ;  /* will count neighbors in this ring */
      vnum = vall_num ;
      for (found = 0, n = old_vnum; vall_num<MAX_VERTICES && n < vall_num; n++)
      {
        vn = &mris->vertices[vall[n]] ;
        if (vn->ripflag)
          continue ;

        /* search through vn's neighbors to find an unmarked vertex */
        for (n2 = 0 ; n2 < vn->vnum ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (vn2->ripflag || vn2->marked)
            continue ;

          /* found one, mark it and put it in the vall list */
          found++ ;
          vn2->marked = nbhd_size ;
          vall[vnum++] = vn->v[n2] ;
          if (nbrs[nbhd_size] > 0)  /* want to store this distance */
          {
            vnbrs[vnbrs_num++] = vn->v[n2] ;
          }
        }
      }  /* done with all neighbors at previous distance */

      /* found all neighbors at this extent - calculate distances */
      old_vnum = vall_num ;
      vall_num += found ;
      for (n = old_vnum ; n < vall_num ; n++)
      {
        vn = &mris->vertices[vall[n]] ;
        for (min_dist = 10000.0, n2 = 0 ; n2 < vn->vnum ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (!vn2->marked)
            continue ;
          xd = vn2->x - vn->x ; yd = vn2->y - vn->y ; zd = vn2->z - vn->z ;
          dist = sqrt(xd*xd + yd*yd + zd*zd) ;
          if (nbhd_size > 1)
            dist /= dist_scale ;
          if (vn2->d + dist < min_dist)
            min_dist = vn2->d + dist ;
        }
        vn->d = min_dist  ;
        if (nbhd_size <= 2)
        {
          xd = vn->x - v->x ; yd = vn->y - v->y ; zd = vn->z - v->z ;
          dist = sqrt(xd*xd + yd*yd + zd*zd) ;
          vn->d = dist ;
        }
      }

      /* if this set of neighbors are to be stored, sample from them */
      if (nbrs[nbhd_size] <= 0)
        continue ;

      /* make sure the points are not too close together */
      min_angle = 0.9*2.0*M_PI / (float)nbrs[nbhd_size] ;

      /*
        at this point, the vall list contains all the neighbors currently found
        at ALL distances, while the vnbrs list contains ONLY the 
        nbhd_size-neighbors.
        */
      if (found <= nbrs[nbhd_size])  /* just copy them all in */
      {
        for (n = 0 ; n < found ; n++, v->vtotal++)
        {
          v->v[v->vtotal] = vnbrs[n] ;
          v->dist_orig[v->vtotal] = mris->vertices[vnbrs[n]].d ;
        }
      }
      else                   /* randomly sample from them */
      {
        int vstart = v->vtotal ;
        for (n = 0 ; n < nbrs[nbhd_size] ; n++, v->vtotal++)
        {
          int j, niter = 0 ;
          do
          {
            do
            {
              i = nint(randomNumber(0.0, (double)found-1)) ;
            } while (vnbrs[i] < 0) ;
            /* 
               now check to make sure that the angle between this
               point and the others already selected is not too
               small to make sure the points are not bunched.
               */
            vn = &mris->vertices[vnbrs[i]] ;
            VECTOR_LOAD(v1, vn->x-v->x, vn->y-v->y, vn->z-v->z) ;
            done = 1 ;
            for (j = vstart ; done && j < v->vtotal ; j++)
            {
              vn2 = &mris->vertices[v->v[j]] ;
              VECTOR_LOAD(v2, vn2->x-v->x, vn2->y-v->y, vn2->z-v->z) ;
              angle = Vector3Angle(v1, v2) ;
              if (angle < min_angle)
                done = 0 ;
            }
            if (++niter > found)  /* couldn't find enough at this difference */
            {
              min_angle *= 0.75f ;  /* be more liberal */
              niter = 0 ;
            }
          } while (!done && !FZERO(min_angle)) ;
          vn = &mris->vertices[vnbrs[i]] ;
          v->v[v->vtotal] = vnbrs[i] ;
          v->dist_orig[v->vtotal] = vn->d ;
          vnbrs[i] = -1 ;
        }
      }
    }


    if ((Gdiag_no == vno) && DIAG_VERBOSE_ON)
    {
      FILE  *fp ;
      char  fname[STRLEN] ;

      sprintf(fname, "v%d", vno) ;
      fp = fopen(fname, "w") ;
      fprintf(fp, "%d\n", vall_num) ;
      for (n = 0 ; n < vall_num ; n++)
        fprintf(fp, "%d\n", vall[n]) ;
      fclose(fp) ;

      sprintf(fname, "vn%d", vno) ;
      fp = fopen(fname, "w") ;
      fprintf(fp, "%d\n", v->vtotal) ;
      for (n = 0 ; n < v->vtotal ; n++)
        fprintf(fp, "%d\n", v->v[n]) ;
      fclose(fp) ;
      for (n = 0 ; n < mris->nvertices ; n++)
      {
        vn = &mris->vertices[n] ;
#if 0
        if (vn->ripflag)
          continue ;
#endif
        vn->curv = vn->d ;
      }
      sprintf(fname, "%s.dist", mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh");
      MRISwriteCurvature(mris, fname) ;
    }

    /* 
       done building arrays - allocate distance vectors and
       sample from the found neighbors list.
       */
    /* now unmark them all */
    for (n = 0 ; n < vall_num ; n++)
    {
      mris->vertices[vall[n]].marked = 0 ;
      mris->vertices[vall[n]].d = 0.0 ;
    }

    total_nbrs += v->vtotal ;
  }

  /* now fill in immediate neighborhood(Euclidean) distances */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak()  ;
    if (v->ripflag)
      continue ;
    if (v->v3num > 0)
      vtotal = v->v3num ;
    else if (v->v2num > 0)
      vtotal = v->v2num ;
    else
      vtotal = v->vnum ;
    for (n = 0 ; n < vtotal ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag)
        continue ;
      xd = v->x - vn->x ; yd = v->y - vn->y ; zd = v->z - vn->z ;
      v->dist_orig[n] = sqrt(xd*xd+yd*yd+zd*zd) ;
    }
  }

  mris->avg_nbrs = (float)total_nbrs / (float)mrisValidVertices(mris) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "avg_nbrs = %2.1f\n", mris->avg_nbrs) ;

  free(vnbrs) ;
  free(vall) ;
  VectorFree(&v1) ; VectorFree(&v2) ;
  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stdout, " done.\n") ;
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Expand the list of neighbors of each vertex, reallocating
          the v->v array to hold the expanded list.
------------------------------------------------------*/
int
MRISsetNeighborhoodSize(MRI_SURFACE *mris, int nsize)
{
  int          i,n, vno, neighbors, j, vnum, niter, nb_vnum ;
  VERTEX       *v, *vnb, *vnb2 ;
  int          vtmp[MAX_NEIGHBORS], vtotal, ntotal ;

/*
   now build a list of 2-connected neighbors. After this is done,
   reallocate the v->n list and arrange the 2-connected neighbors
   sequentially after it.
*/
  for (niter = 0 ; niter < nsize-mris->nsize ; niter++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak()  ;
      if (v->nsize >= nsize)
        continue ;
      vnum = v->vtotal ;
      if (v->ripflag || !vnum)
        continue ;
      memmove(vtmp, v->v, vnum*sizeof(int)) ;

      /* mark 1-neighbors so we don't count them twice */
      v->marked = 1 ;
      for (i = 0 ; i < vnum; i++)
        mris->vertices[v->v[i]].marked = 1 ;

      /* count 2-neighbors */
      for (neighbors = vnum, i = 0; neighbors < MAX_NEIGHBORS && i < vnum; i++)
      {
        n = v->v[i] ;
        vnb = &mris->vertices[n] ;
        vnb->marked = 1 ;
        if (vnb->ripflag)
          continue ;

        switch (v->nsize)
        {
        case 1: nb_vnum = vnb->vnum ;    break ;
        case 2: nb_vnum = vnb->v2num ;   break ;
        default:
        case 3: nb_vnum = vnb->v3num ;  break ;
        }

        for (j = 0 ; j < nb_vnum ; j++)
        {
          vnb2 = &mris->vertices[vnb->v[j]] ;
          if (vnb2->ripflag || vnb2->marked)
            continue ;
          vtmp[neighbors] = vnb->v[j] ;
          vnb2->marked = 1 ;
          if (++neighbors >= MAX_NEIGHBORS)
          {
            fprintf(stderr, "vertex %d has too many neighbors!\n",vno) ;
            break ;
          }
        }
      }
      /*
        now reallocate the v->v structure and place the 2-connected neighbors
        suquentially after the 1-connected neighbors.
        */
      free(v->v) ;
      v->v = (int *)calloc(neighbors, sizeof(int)) ;
      if (!v->v)
        ErrorExit(ERROR_NO_MEMORY, 
                  "MRISsetNeighborhoodSize: could not allocate list of %d "
                  "nbrs at v=%d", neighbors, vno) ;

      v->marked = 0 ;
      for (n = 0 ; n < neighbors ; n++)
      {
        v->v[n] = vtmp[n] ;
        mris->vertices[vtmp[n]].marked = 0 ;
      }
      if (v->dist)
        free(v->dist) ;
      if (v->dist_orig)
        free(v->dist_orig) ;
          
      v->dist = (float *)calloc(neighbors, sizeof(float)) ;
      if (!v->dist)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsetNeighborhoodSize: could not allocate list of %d "
                  "dists at v=%d", neighbors, vno) ;
      v->dist_orig = (float *)calloc(neighbors, sizeof(float)) ;
      if (!v->dist_orig)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsetNeighborhoodSize: could not allocate list of %d "
                  "dists at v=%d", neighbors, vno) ;
      v->nsize++ ;
      switch (v->nsize)
      {
      case 2:
        v->v2num = neighbors ;
        break ;
      case 3:
        v->v3num = neighbors ;
        break ;
      default:   /* store old neighborhood size in v3num */
        v->v3num = v->vtotal ;
        break ;
      }
      v->vtotal = neighbors ;
      for (n = 0 ; n < neighbors ; n++)
        for (i = 0 ; i < neighbors ; i++)
          if (i != n && v->v[i] == v->v[n])
            fprintf(stderr, 
                    "warning: vertex %d has duplicate neighbors %d and %d!\n",
                    vno, i, n) ;
      if ((vno == Gdiag_no) && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
      {
        fprintf(stdout, "v %d: vnum=%d, v2num=%d, vtotal=%d\n",
                vno, v->vnum, v->v2num, v->vtotal) ;
        for (n = 0 ; n < neighbors ; n++)
          fprintf(stdout, "v[%d] = %d\n", n, v->v[n]) ;
      }
    }
  }

  for (vno = ntotal = vtotal = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    vtotal += v->vtotal ;
    ntotal++ ;
  }

  mris->avg_nbrs = (float)vtotal / (float)ntotal ;
  mris->nsize = nsize ;
  if (Gdiag & DIAG_SHOW && mris->nsize > 1 && DIAG_VERBOSE_ON)
    fprintf(stdout, "avg_nbrs = %2.1f\n", mris->avg_nbrs) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Remove ripped vertices and faces from the v->v and the
          v->f arrays respectively.
------------------------------------------------------*/
int
MRISremoveRipped(MRI_SURFACE *mris)
{
  int     vno, n, fno, nripped ;
  VERTEX  *v ;
  FACE    *face ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "removing ripped vertices and faces...\n") ;
  do
  {
    nripped = 0 ;
    // go through all vertices
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      // if rip flag set
      if (v->ripflag)
      {
	// remove it
        if (v->dist)
          free(v->dist) ;
        if (v->dist_orig)
          free(v->dist_orig) ;
        v->dist = v->dist_orig = NULL ;
        continue ;
      }

      for (n = 0 ; n < v->vtotal ; n++)
      {
        /* remove this vertex from neighbor list if it is ripped */
        if (mris->vertices[v->v[n]].ripflag)
        {
          if (n < v->vtotal-1)  /* not the last one in the list */
          {
            memmove(v->v+n, v->v+n+1, (v->vtotal-n-1)*sizeof(int)) ;
            memmove(v->dist+n, v->dist+n+1, (v->vtotal-n-1)*sizeof(float)) ;
            memmove(v->dist_orig+n, v->dist_orig+n+1, 
                    (v->vtotal-n-1)*sizeof(float)) ;
          }
          if (n < v->vnum)      /* it was a 1-neighbor */
            v->vnum-- ;
          if (n < v->v2num)     /* it was a 2-neighbor */
            v->v2num-- ;
          if (n < v->v3num)     /* it was a 3-neighbor */
            v->v3num-- ;
          n-- ; v->vtotal-- ;
        }
      }
      for (fno = 0 ; fno < v->num ; fno++)
      {
        /* remove this face from face list if it is ripped */
        if (mris->faces[v->f[fno]].ripflag)
        {
          if (fno < v->num-1)  /* not the last one in the list */
          {
            memmove(v->f+fno, v->f+fno+1, (v->num-fno-1)*sizeof(int)) ;
            memmove(v->n+fno, v->n+fno+1, (v->num-fno-1)*sizeof(uchar)) ;
          }
          v->num-- ; fno-- ;
        }
      }
      if (v->num <= 0 || v->vnum <= 0)  /* degenerate vertex */
      {
        v->ripflag = 1 ;
        nripped++ ;
      }
    }
  } while (nripped > 0) ;

  /* now recompute total original area for scaling */
  mris->orig_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    mris->orig_area += face->orig_area ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          no triangle area in msurfer, no explodeflag here 
------------------------------------------------------*/
#define RAN   0.001   /* one thousandth of a millimeter */

int
MRIScomputeNormals(MRI_SURFACE *mris) 
{
  int       k,n, num, i ;
  VERTEX    *v ;
  FACE      *f;
  float     norm[3],snorm[3], len ;

  i = 0 ;

#if 0
  if (mris->status == MRIS_PLANE)
  {
    mrisComputeBoundaryNormals(mris);
    mrisSmoothBoundaryNormals(mris,10);
  }
#endif
  for (k=0;k<mris->nfaces;k++) if (mris->faces[k].ripflag)
  {
    f = &mris->faces[k];
    for (n=0;n<VERTICES_PER_FACE;n++)
      mris->vertices[f->v[n]].border = TRUE;
  }

  for (k=0;k<mris->nvertices;k++) if (!mris->vertices[k].ripflag)
  {
    v = &mris->vertices[k];
    snorm[0]=snorm[1]=snorm[2]=0;
    v->area = 0;
    for (num = n=0;n<v->num;n++) if (!mris->faces[v->f[n]].ripflag)
    {
      num++ ;
      mrisNormalFace(mris, v->f[n], (int)v->n[n],norm);
      snorm[0] += norm[0];
      snorm[1] += norm[1];
      snorm[2] += norm[2];

      /* Note: overestimates area by *2 !! */
      v->area += mrisTriangleArea(mris, v->f[n], (int)v->n[n]); 
    }
    if (!num)
      continue ;
    mrisNormalize(snorm);

    v->area /= 2 ;
    if (v->origarea<0)        /* has never been set */
      v->origarea = v->area;

    len = sqrt(snorm[0]*snorm[0] + snorm[1]*snorm[1] + snorm[2]*snorm[2]) ;
    if (!FZERO(len))
    {
      v->nx = snorm[0];
      v->ny = snorm[1];
      v->nz = snorm[2];
      i = 0 ;
    }
    else
    {
      if (i++ > 5)
        continue ;

      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stderr, "vertex %d: degenerate normal\n", k) ;

      v->x += (float)randomNumber(-RAN, RAN) ;
      v->y += (float)randomNumber(-RAN, RAN) ;
      /* normal is always (0,0,+-1) anyway */
      if (mris->status == MRIS_PLANE || mris->status == MRIS_CUT)
      {
        v->nx = v->ny = 0.0f ; v->nz = 1.0f ;
        continue ;
      }

      v->z += (float)randomNumber(-RAN, RAN) ;
      for (n=0;n<v->vnum;n++) /*if (!mris->faces[v->f[n]].ripflag)*/
      {
        mris->vertices[v->v[n]].x += (float)randomNumber(-RAN, RAN) ;
        mris->vertices[v->v[n]].y += (float)randomNumber(-RAN, RAN) ;
        mris->vertices[v->v[n]].z += (float)randomNumber(-RAN, RAN) ;
      }
      k-- ;   /* recalculate the normal for this vertex */
    }
  }
#if 0
  mris->vertices[0].nx = mris->vertices[0].ny = 0 ;
  mris->vertices[0].nz = mris->vertices[0].nz / fabs(mris->vertices[0].nz) ;
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        calculate distances between each vertex and all of its neighbors
------------------------------------------------------*/
static int
mrisComputeVertexDistances(MRI_SURFACE *mris)
{
  int     vno, n, vtotal, *pv ;
  VERTEX  *v, *vn ;
  float   d, xd, yd, zd, circumference = 0.0f, angle ;
  VECTOR  *v1, *v2 ;

  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;

  for (vno=0;vno<mris->nvertices;vno++)
  {
    v = &mris->vertices[vno];
    if (v->ripflag || v->dist == NULL)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vtotal = v->vtotal ;
    switch (mris->status)
    {
    default:   /* don't really know what to do in other cases */
    case MRIS_PLANE:
      for (pv = v->v, n = 0 ; n < vtotal ; n++)
      {
        vn = &mris->vertices[*pv++] ;
        if (vn->ripflag)
          continue ;
        xd = v->x - vn->x ; yd = v->y - vn->y ; zd = v->z - vn->z ;
        d = xd*xd + yd*yd + zd*zd ;
        v->dist[n] = sqrt(d) ;
      }
      break ;
    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE:
      VECTOR_LOAD(v1, v->x, v->y, v->z) ;  /* radius vector */
      if (FZERO(circumference))   /* only calculate once */
        circumference = M_PI * 2.0 * V3_LEN(v1) ;
      for (pv = v->v, n = 0 ; n < vtotal ; n++)
      {
        vn = &mris->vertices[*pv++] ;
        if (vn->ripflag)
          continue ;
        VECTOR_LOAD(v2, vn->x, vn->y, vn->z) ;  /* radius vector */
        angle = fabs(Vector3Angle(v1, v2)) ;
#if 0
        xd = v->x - vn->x ; yd = v->y - vn->y ; zd = v->z - vn->z ;
        d = sqrt(xd*xd + yd*yd + zd*zd) ;
#endif
        d = circumference * angle / (2.0 * M_PI) ;
        v->dist[n] = d ;
      }
      break ;
    }
  }
  
  VectorFree(&v1) ; VectorFree(&v2) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static void
mrisNormalize(float v[3])
{
  float d;

  d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (d>0)
  {
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
  }
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static float
mrisTriangleArea(MRIS *mris, int fac, int n)
{
  int n0,n1;
  face_type *f;
  float v0[3],v1[3],d1,d2,d3;

  n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
  v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
  v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
  v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
  v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
  v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
  d1 = -v1[1]*v0[2] + v0[1]*v1[2];
  d2 = v1[0]*v0[2] - v0[0]*v1[2];
  d3 = -v1[0]*v0[1] + v0[0]*v1[1];
  return sqrt(d1*d1+d2*d2+d3*d3)/2;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeOrigNormal(MRIS *mris, int vno, float norm[])
{
  float snorm[3] ;
  VERTEX *v ;
  int    n, num ;

  v = &mris->vertices[vno] ;

  norm[0]=norm[1]=norm[2]=0.0;
  for (num = n=0;n<v->num;n++) if (!mris->faces[v->f[n]].ripflag)
  {
    num++ ;
    mrisOrigNormalFace(mris, v->f[n], (int)v->n[n],snorm);
    norm[0] += snorm[0];
    norm[1] += snorm[1];
    norm[2] += snorm[2];

  }
  if (!num)
    return(ERROR_BADPARM) ;
  mrisNormalize(norm);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisOrigNormalFace(MRIS *mris, int fac,int n,float norm[])
{
  int     n0,n1, *pv ;
  FACE    *f;
  float   v0[3],v1[3];
  register VERTEX  *v, *vn0, *vn1 ;

  n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  pv = f->v ;
  vn0 = &mris->vertices[pv[n0]] ;
  vn1 = &mris->vertices[pv[n1]] ;
  v =  &mris->vertices[pv[n]] ;
  v0[0] = v->origx - vn0->origx; v0[1] = v->origy - vn0->origy; 
  v0[2] = v->origz - vn0->origz;
  v1[0] = vn1->origx - v->origx; v1[1] = vn1->origy - v->origy; 
  v1[2] = vn1->origz - v->origz;
  mrisNormalize(v0);
  mrisNormalize(v1);
  norm[0] = -v1[1]*v0[2] + v0[1]*v1[2];
  norm[1] = v1[0]*v0[2] - v0[0]*v1[2];
  norm[2] = -v1[0]*v0[1] + v0[0]*v1[1];
/*
  printf("[%5.2f,%5.2f,%5.2f] x [%5.2f,%5.2f,%5.2f] = [%5.2f,%5.2f,%5.2f]\n",
         v0[0],v0[1],v0[2],v1[0],v1[1],v1[2],norm[0],norm[1],norm[2]);
*/
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisNormalFace(MRIS *mris, int fac,int n,float norm[])
{
  int     n0,n1, *pv ;
  FACE    *f;
  float   v0[3],v1[3];
  register VERTEX  *v, *vn0, *vn1 ;

  n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  pv = f->v ;
  vn0 = &mris->vertices[pv[n0]] ;
  vn1 = &mris->vertices[pv[n1]] ;
  v =  &mris->vertices[pv[n]] ;
  v0[0] = v->x - vn0->x; v0[1] = v->y - vn0->y; v0[2] = v->z - vn0->z;
  v1[0] = vn1->x - v->x; v1[1] = vn1->y - v->y; v1[2] = vn1->z - v->z;
  mrisNormalize(v0);
  mrisNormalize(v1);
  norm[0] = -v1[1]*v0[2] + v0[1]*v1[2];
  norm[1] = v1[0]*v0[2] - v0[0]*v1[2];
  norm[2] = -v1[0]*v0[1] + v0[0]*v1[1];
/*
  printf("[%5.2f,%5.2f,%5.2f] x [%5.2f,%5.2f,%5.2f] = [%5.2f,%5.2f,%5.2f]\n",
         v0[0],v0[1],v0[2],v1[0],v1[1],v1[2],norm[0],norm[1],norm[2]);
*/
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisReadTransform(MRIS *mris, char *mris_fname)
{
  char transform_fname[STRLEN], fpref[300] ;

  FileNamePath(mris_fname, fpref) ;
  sprintf(transform_fname, "%s/../mri/transforms/talairach.xfm", fpref) ;
  if (!FileExists(transform_fname))
    return(ERROR_NO_FILE) ;
  if (input_transform_file(transform_fname, &mris->transform) != OK)
  {
    ErrorReturn(ERROR_NO_FILE, 
                (ERROR_NOFILE, 
                 "mrisReadTransform: could not read xform file '%s'", 
                 transform_fname)) ;
  }
  else
  {
    mris->linear_transform = get_linear_transform_ptr(&mris->transform) ;
    mris->inverse_linear_transform = 
      get_inverse_linear_transform_ptr(&mris->transform) ;
    mris->free_transform = 1 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadBinaryCurvature(MRI_SURFACE *mris, char *mris_fname)
{
  char   fname[STRLEN], fpref[STRLEN], hemi[20] ;

  FileNamePath(mris_fname, fpref) ;
  strcpy(hemi, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh") ;
  sprintf(fname, "%s/%s.curv", fpref, hemi) ;
  return(MRISreadCurvatureFile(mris, fname)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadCurvatureFile(MRI_SURFACE *mris, char *sname)
{
  int    k,i,vnum,fnum;
  float  curv, curvmin, curvmax;
  FILE   *fp;
  char   *cp, path[STRLEN], fname[STRLEN], type ;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (cp && ((strncmp(cp-2, "lh", 2) == 0) || (strncmp(cp-2, "rh", 2) == 0)))
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explicitly */


  type = MRISfileNameType(fname) ;
  if (type == MRIS_ASCII_TRIANGLE_FILE)
    return(mrisReadAsciiCurvatureFile(mris, fname)) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
    fprintf(stdout, "reading curvature file...") ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadCurvature: could not open %s", 
                 fname)) ;

  fread3(&vnum,fp);
  if (vnum == NEW_VERSION_MAGIC_NUMBER)
  {
    fclose(fp) ;
    return(MRISreadNewCurvatureFile(mris, fname)) ;
  }
  
  fread3(&fnum,fp);
  if (vnum!= mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadBinaryCurvature: incompatible vertex "
                 "number in file %s", fname)) ;
  }
  curvmin = 10000.0f ; curvmax = -10000.0f ;  /* for compiler warnings */
  for (k=0;k<vnum;k++)
  {
    fread2(&i,fp);
    curv = i/100.0;
    if (k==0) curvmin=curvmax=curv;
    if (curv>curvmax) curvmax=curv;
    if (curv<curvmin) curvmin=curv;
    mris->vertices[k].curv = curv;
  }
  mris->max_curv = curvmax ;
  mris->min_curv = curvmin ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "done. min=%2.3f max=%2.3f\n", curvmin, curvmax) ;
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float *
MRISreadCurvatureVector(MRI_SURFACE *mris, char *sname)
{
  int    k,i,vnum,fnum;
  float  *cvec ;
  FILE   *fp;
  char   *cp, path[STRLEN], fname[STRLEN] ;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (cp)
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
    fprintf(stdout, "reading curvature file...") ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    return(NULL) ;

  fread3(&vnum,fp);
  if (vnum == NEW_VERSION_MAGIC_NUMBER)
  {
    fclose(fp) ;
    return(MRISreadNewCurvatureVector(mris, fname)) ;
  }
  
  fread3(&fnum,fp);
  if (vnum!= mris->nvertices)
  {
    fclose(fp) ;
    return(NULL) ;
  }
  cvec = (float *)calloc(mris->nvertices, sizeof(float)) ;
  if (!cvec)
    ErrorExit(ERROR_NOMEMORY, "MRISreadCurvatureVector(%s): calloc failed",
              fname) ;

  for (k=0;k<vnum;k++)
  {
    fread2(&i,fp);
    cvec[k] = i/100.0 ;
  }
  fclose(fp);
  return(cvec) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadFloatFile(MRI_SURFACE *mris, char *sname)
{
  int    k,vnum,fnum;
  float  f, fmin, fmax;
  FILE   *fp;
  char   *cp, path[STRLEN], fname[STRLEN] ;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (cp)
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
    fprintf(stdout, "reading float file...") ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadFloatFile: could not open %s", 
                 fname)) ;

  vnum = freadInt(fp);
  fnum = freadInt(fp);
  if (vnum!= mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadFloatFile: incompatible # of vertices "
                 "in file %s", fname)) ;
  }
  if (fnum!= mris->nfaces)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadFloatFile: incompatible # of faces "
                 "file %s", fname)) ;
  }
  fmin = 10000.0f ; fmax = -10000.0f ;  /* for compiler warnings */
  for (k=0;k<vnum;k++)
  {
    f = freadFloat(fp);
    mris->vertices[k].val = f;
  }
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadBinaryAreas(MRI_SURFACE *mris, char *mris_fname)
{
  int   k,vnum,fnum;
  float f;
  FILE  *fp;
  char  fname[STRLEN], fpref[STRLEN], hemi[20] ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "reading area file...") ;

  FileNamePath(mris_fname, fpref) ;
  strcpy(hemi, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh") ;
  sprintf(fname, "%s/%s.area", fpref, hemi) ;

  /*  mris->orig_area = 0.0f ;*/
  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADPARM, 
              (ERROR_BADPARM,"MRISreadBinaryAreas: no area file %s\n",fname));
  fread3(&vnum,fp);
  fread3(&fnum,fp);
  if (vnum!=mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadBinaryAreas: incompatible vertex "
                 "number in file %s", fname)) ;
  }

  for (k=0;k<vnum;k++)
  {
    f = freadFloat(fp);
    mris->vertices[k].origarea = f ;
    /*    mris->orig_area += f;*/
  }
  fclose(fp);

  /* hack to correct for overestimation of area in compute_normals */
#if 0
  mris->orig_area /= 2; 
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "total area = %2.0f.\n", mris->orig_area) ;
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteArea(MRI_SURFACE *mris, char *sname)
{
  float  *curv_save ;

  curv_save = (float *)calloc(mris->nvertices, sizeof(float)) ;
  if (!curv_save)
    ErrorExit(ERROR_NOMEMORY, "MRISwriteArea: could not alloc %d vertex curv storage",
              mris->nvertices) ;

  MRISextractCurvatureVector(mris, curv_save) ;
  MRISareaToCurv(mris) ;
  MRISwriteCurvature(mris, sname) ;
  MRISimportCurvatureVector(mris, curv_save) ;
  free(curv_save) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisReadFieldsign(MRI_SURFACE *mris, char *mris_fname)
{
  int k,i,vnum;
  float f;
  FILE *fp;

  printf("surfer: read_fieldsign(%s)\n",fname);
  fp = fopen(fname,"r");
  if (fp==NULL) {printf("surfer: ### File %s not found\n",fname);PR return;}
  fread(&vnum,1,sizeof(int),fp);
  printf("surfer: vertex_index = %d, vnum = %d\n",vertex_index,vnum);
  if (vnum!=vertex_index)
    printf("surfer: Warning: incompatible vertex number in file %s\n",fname);
  for (k=0;k<vnum;k++)
  {
    fread(&f,1,sizeof(float),fp);
    vertex[k].fieldsign = f;
  }
  fclose(fp);
  fieldsignflag = TRUE;
  PR
    return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
normalize_binary_curvature(MRI_SURFACE *mris)
{
  int k;
  float curv,min,max;
  float sum,avg,sq,sum_sq,sd,n;
  FILE *fp;

  if (!curvloaded)   { printf("surfer: ### curv not loaded!\n");PR return; }

  sum = 0;
  for (k=0;k<vertex_index;k++)
    sum += vertex[k].curv;
  avg = sum/vertex_index;

  n = (float)vertex_index;
  sum = sum_sq = 0.0;
  for (k=0;k<vertex_index;k++) {
    vertex[k].curv -= avg;
    curv = vertex[k].curv;
    sum += curv;
    sum_sq += curv*curv;
  }
  sd = sqrt((n*sum_sq - sum*sum)/(n*(n-1.0)));

  for (k=0;k<vertex_index;k++) {
    curv = (vertex[k].curv)/sd;
    if (k==0) min=max=curv;
    if (curv<min) min=curv;
    if (curv>max) max=curv;
    if (curv<CURVIM_NORM_MIN) curv = CURVIM_NORM_MIN;
    if (curv>CURVIM_NORM_MAX) curv = CURVIM_NORM_MAX;
    vertex[k].curv = curv;
  }
  curvmin = CURVIM_NORM_MIN;
  curvmax = CURVIM_NORM_MAX;
  printf("surfer: curvature normalized: avg=%f sd=%f\n",avg,sd);
  printf("surfer: min=%f max=%f trunc to (%f,%f)\n",min,max,curvmin,curvmax);PR
}

#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Perform a projection onto a cylinder moving each
          point on the cortical surface to the closest cylindrical
          coordinate.
------------------------------------------------------*/
int
MRISprojectOntoCylinder(MRI_SURFACE *mris, float radius)
{
  VERTEX  *v;
  int     k;
  float   x,y,z,x2,z2,dx,dz ;
  float   d ;


  MRIScenter(mris, mris) ;

  for (k=0;k<mris->nvertices;k++) 
  {
    v = &mris->vertices[k];
    x = v->x;
    y = v->y;
    z = v->z;

    x2 = x*x;
    z2 = z*z;

    d = (-1.0+(float)radius/sqrt(x2+z2)) ;
    if (!finite(d))
    {
      ErrorPrintf(ERROR_BADPARM, 
                  "point (%2.2f,%2.2f,%2.2f) cannot be projected on cylinder",
                  x, y, z) ;
      
    }
    dx = d*x ;
    dz = d*z ;
    v->x = x+dx ;
    v->z = z+dz;

    if (!finite(v->x) || !finite(v->y) || !finite(v->z))
      DiagBreak() ;

  }
  MRISupdateSurface(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Perform a projection onto an sphere moving each
          point on the cortical surface to the closest spherical
          coordinate.
------------------------------------------------------*/
MRI_SURFACE  *
MRISprojectOntoSphere(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, double r)
{
  VERTEX  *v;
  int     vno ;
  double  x, y, z, d, dx, dy, dz, dist, total_dist, x2, y2, z2 ;

  if (FZERO(r))
    r = DEFAULT_RADIUS ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  if ((mris_dst->status != MRIS_SPHERE) &&
      (mris_dst->status != MRIS_PARAMETERIZED_SPHERE))
    MRIScenter(mris_dst, mris_dst) ;

  mris_dst->radius = r ;

  for (total_dist = vno = 0 ; vno < mris_dst->nvertices ; vno++) 
  {
    v = &mris_dst->vertices[vno];
    if (v->ripflag)  /* shouldn't happen */
      continue ;
    if (vno == 118009)
      { DiagBreak() ; }
    x = (double)v->x;
    y = (double)v->y;
    z = (double)v->z;

    x2 = x*x ; y2 = y*y ; z2 = z*z ;
    dist = sqrt(x2+y2+z2) ;
    if (FZERO(dist))
      d = 0 ;
    else
      d = 1 - r / dist ;
    dx = d*x ;
    dy = d*y;
    dz = d*z;
    v->x = x-dx ;
    v->y = y-dy;
    v->z = z-dz;

    if (!finite(v->x) || !finite(v->y) || !finite(v->z))
      DiagBreak() ;

    /*    if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)*/
    {
      dist = sqrt((double)(dx*dx+dy*dy+dz*dz));
      total_dist += dist;
    }
#if 1
x = (double)v->x;
y = (double)v->y;
z = (double)v->z;
x2 = x*x ; y2 = y*y ; z2 = z*z ;
dist = sqrt(x2+y2+z2) ;
#endif

  }
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
    fprintf(stdout, 
            "sphere_project: total dist = %f\n",total_dist);
  MRISupdateEllipsoidSurface(mris_dst) ;
  mris_dst->status = mris_src->status == MRIS_PARAMETERIZED_SPHERE ? 
    MRIS_PARAMETERIZED_SPHERE : MRIS_SPHERE ;
  return(mris_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Perform a projection onto an ellipsoid moving each
          point on the cortical surface to the closest ellipsoidal
          coordinate.
------------------------------------------------------*/
extern double sqrt(double) ;

MRI_SURFACE *
MRISprojectOntoEllipsoid(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, 
                                       float a, float b, float c)
{
  VERTEX  *v;
  int     k;
  float   x,y,z,x2,y2,z2,dx,dy,dz,a2,b2,c2,a4,b4,c4,a6,b6,c6;
  float   f,g,h,d,dist,avgdist=0.0f ;

  if (FZERO(a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  MRIScenter(mris_dst, mris_dst) ;

  mris_dst->a = a ; mris_dst->b = b ; mris_dst->c = c ;

  /*  printf("ellipsoid_project(%f,%f,%f)\n",a,b,c);*/
  a2 = a*a;
  b2 = b*b;
  c2 = c*c;
  a4 = a2*a2;
  b4 = b2*b2;
  c4 = c2*c2;
  a6 = a2*a4;
  b6 = b2*b4;
  c6 = c2*c4;

#if 0
  /* rescale brain so that it is contained within the ellipsoid */
  xscale = mris_dst->xhi / a ;
  yscale = mris_dst->yhi / b ;
  zscale = mris_dst->zhi / c ;
  if ((xscale > yscale) && (xscale > zscale))
    scale = 1.0f / xscale ;
  else if (yscale > zscale)
    scale = 1.0f / yscale ;
  else
    scale = 1.0f / zscale ;

  MRISscaleBrain(mris_dst, mris_dst, scale) ;
#endif

  for (k=0;k<mris_dst->nvertices;k++) 
  {
    v = &mris_dst->vertices[k];
/*
    printf("%6d: before: %6.2f\n",k,SQR(v->x/a)+SQR(v->y/b)+SQR(v->z/c));
*/
    x = v->x;
    y = v->y;
    z = v->z;
#if 0
    if ((fabs(x) > a) || (fabs(y) > b) || (fabs(z) > c))
      return(MRISradialProjectOntoEllipsoid(mris_src, mris_dst, a, b, c)) ;
#endif

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    f = x2/a6+y2/b6+z2/c6;
    g = 2*(x2/a4+y2/b4+z2/c4);
    h = x2/a2+y2/b2+z2/c2-1;
    d = (-g+(float)sqrt((double)(g*g-4*f*h)))/(2*f);
    if (!finite(d))
    {
      ErrorPrintf(ERROR_BADPARM, 
              "point (%2.2f,%2.2f,%2.2f) cannot be projected on ell "
              "(%2.0f,%2.0f,%2.0f...\n",
              x, y, z, a, b, c) ;
      
      return(MRISradialProjectOntoEllipsoid(mris_src, mris_dst, a, b, c)) ;
    }
    dx = d*x/a2;
    dy = d*y/b2;
    dz = d*z/c2;
    v->x = x+dx ;
    v->y = y+dy;
    v->z = z+dz;

    if (!finite(v->x) || !finite(v->y) || !finite(v->z))
      DiagBreak() ;

    if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
    {
      dist = (float)sqrt((double)(dx*dx+dy*dy+dz*dz));
      avgdist += dist;
    }
/*
    printf("%6d: after: %6.2f\n",k,SQR(v->x/a)+SQR(v->y/b)+SQR(v->z/c));
*/
  }
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
    fprintf(stdout, 
            "ellipsoid_project: avgdist = %f\n",avgdist/mris_dst->nvertices);
  MRISupdateEllipsoidSurface(mris_dst) ;
  if (FZERO(a-b) && FZERO(b-c))
    mris_dst->status = MRIS_SPHERE ;
  else
    mris_dst->status = MRIS_ELLIPSOID ;
  return(mris_dst) ;
}


/*
  this one projects along the line from the origin to the ellipsoidal
  surface - not orthographic unless the ellipsoid is a sphere.
  */
MRI_SURFACE *
MRISradialProjectOntoEllipsoid(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, 
                                       float a, float b, float c)
{
  int    vno ;
  VERTEX *vsrc, *vdst ;
  float  x0, y0, z0, x1, y1, z1, denom,
         asq_bsq, asq_csq, bsq_csq, x1sq, y1sq, z1sq, abc ;

  if (FZERO(a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  x0 = mris_dst->xctr ; y0 = mris_dst->yctr ; z0 = mris_dst->zctr ;
  asq_bsq = a*a*b*b ; bsq_csq= b*b*c*c ; asq_csq = a*a*c*c ; abc = a * b * c ;

  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vsrc = &mris_src->vertices[vno] ;
    vdst = &mris_dst->vertices[vno] ;
    x1 = (vsrc->x-x0) ; y1 = (vsrc->y-y0) ; z1 = (vsrc->z-z0) ;
    x1sq = x1*x1 ; y1sq = y1*y1 ; z1sq = z1*z1 ;

    /* right out of mathematica (almost) */
    denom = sqrt(bsq_csq*x1sq + asq_csq*y1sq + asq_bsq*z1sq) ;

    vdst->x = abc*x1 / denom /* + x0 */ ;
    vdst->y = abc*y1 / denom /* + y0 */ ;
    vdst->z = abc*z1 / denom /* + z0 */ ;
  }
  
  x0 = y0 = z0 = 0 ;   /* set center of ellipsoid at origin */
#if 0
  if (mris_dst->v_temporal_pole)
  {
    mris_dst->v_temporal_pole->x = x0 ;
    mris_dst->v_temporal_pole->y = y0 ;
    mris_dst->v_temporal_pole->z = -c+z0 ;
    mris_dst->v_temporal_pole->tethered = TETHERED_TEMPORAL_POLE ;
  }
  if (mris_dst->v_frontal_pole)
  {
    mris_dst->v_frontal_pole->x = x0 ;
    mris_dst->v_frontal_pole->y = b+y0 ;
    mris_dst->v_frontal_pole->z = z0 ;
    mris_dst->v_frontal_pole->tethered = TETHERED_FRONTAL_POLE ;
  }
  if (mris_dst->v_occipital_pole)
  {
    mris_dst->v_occipital_pole->x = x0 ;
    mris_dst->v_occipital_pole->y = -b+y0 ;
    mris_dst->v_occipital_pole->z = z0 ;
    mris_dst->v_occipital_pole->tethered = TETHERED_OCCIPITAL_POLE ;
  }
#endif


  MRISupdateEllipsoidSurface(mris_dst) ;
  if (FZERO(a-b) && FZERO(b-c))
    mris_dst->status = MRIS_SPHERE ;
  else
    mris_dst->status = MRIS_ELLIPSOID ;
  return(mris_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISclone(MRI_SURFACE *mris_src)
{
  MRI_SURFACE *mris_dst ;
  int         vno, fno, n ;
  VERTEX      *vsrc, *vdst ;
  FACE        *fsrc, *fdst ;

  mris_dst = MRISalloc(mris_src->nvertices, mris_src->nfaces) ;

  mris_dst->hemisphere = mris_src->hemisphere ;
  mris_dst->xctr = mris_src->xctr ;
  mris_dst->yctr = mris_src->yctr ;
  mris_dst->zctr = mris_src->zctr ;
  mris_dst->xlo = mris_src->xlo ;
  mris_dst->ylo = mris_src->ylo ;
  mris_dst->zlo = mris_src->zlo ;
  mris_dst->xhi = mris_src->xhi ;
  mris_dst->yhi = mris_src->yhi ;
  mris_dst->zhi = mris_src->zhi ;
  mris_dst->min_curv = mris_src->min_curv ;
  mris_dst->max_curv = mris_src->max_curv ;
  mris_dst->total_area = mris_src->total_area ;
  mris_dst->orig_area = mris_src->orig_area ;
  mris_dst->linear_transform = mris_src->linear_transform ;
  mris_dst->inverse_linear_transform = mris_src->inverse_linear_transform ;
  mris_dst->free_transform = 0 ;
  if (mris_src->v_frontal_pole)
    mris_dst->v_frontal_pole = 
      &mris_dst->vertices[mris_src->v_frontal_pole - mris_src->vertices] ;
  if (mris_src->v_occipital_pole)
    mris_dst->v_occipital_pole = 
      &mris_dst->vertices[mris_src->v_occipital_pole - mris_src->vertices] ;
  if (mris_src->v_temporal_pole)
    mris_dst->v_temporal_pole = 
      &mris_dst->vertices[mris_src->v_temporal_pole - mris_src->vertices] ;
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vsrc = &mris_src->vertices[vno] ;
    vdst = &mris_dst->vertices[vno] ;
    vdst->x = vsrc->x ;
    vdst->y = vsrc->y ;
    vdst->z = vsrc->z ;
    vdst->nx = vsrc->nx ;
    vdst->ny = vsrc->ny ;
    vdst->nz = vsrc->nz ;
    vdst->cx = vsrc->cx ;
    vdst->cy = vsrc->cy ;
    vdst->cz = vsrc->cz ;
#if 0
    vdst->ox = vsrc->ox ;
    vdst->oy = vsrc->oy ;
    vdst->oz = vsrc->oz ;
#endif
    vdst->curv = vsrc->curv ;
    vdst->num = vsrc->num ;

    if (vdst->num)
    {
      vdst->f = (int *)calloc(vdst->num,sizeof(int));
      if (!vdst->f)
        ErrorExit(ERROR_NO_MEMORY, "MRISclone: could not allocate %d faces",
                  vdst->num) ;
      vdst->n = (uchar *)calloc(vdst->num,sizeof(uchar));
      if (!vdst->n)
        ErrorExit(ERROR_NO_MEMORY, "MRISclone: could not allocate %d num",
                  vdst->n) ;
      for (n = 0; n < vdst->num; n++)
      {
        vdst->n[n] = vsrc->n[n] ;
        vdst->f[n] = vsrc->f[n] ;
      }
    }

    vdst->vnum = vsrc->vnum ;
    vdst->v2num = vsrc->v2num ;
    vdst->v3num = vsrc->v3num ;
    vdst->vtotal = vsrc->vtotal ;
    if (vdst->vnum)
    {
      vdst->v = (int *)calloc(vdst->vtotal, sizeof(int)) ;
      if (!vdst->v)
        ErrorExit(ERROR_NO_MEMORY, "MRISclone: could not allocate %d nbrs",
                  vdst->vtotal) ;
      for (n = 0; n < vdst->vtotal; n++)
        vdst->v[n] = vsrc->v[n] ;
    }

    vdst->ripflag = vsrc->ripflag ;
#if 0
    vdst->oripflag = vsrc->oripflag ;
    vdst->origripflag = vsrc->origripflag ;
    memcpy(vdst->coord, vsrc->coord, sizeof(vsrc->coord)) ;
#endif
    vdst->border = vsrc->border ;
    vdst->area = vsrc->area ;
    vdst->origarea = vsrc->origarea ;
  }

  for (fno = 0 ; fno < mris_src->nfaces ; fno++)
  {
    fsrc = &mris_src->faces[fno] ;
    fdst = &mris_dst->faces[fno] ;
    memmove(fdst, fsrc, sizeof(FACE)) ;
  }
  return(mris_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRIScenter(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int         vno ;
  VERTEX      *vdst ;
  float       x, y, z, x0, y0, z0, xlo, xhi, zlo, zhi, ylo, yhi ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  x = y = z = 0 ;   /* silly compiler warning */
  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->ripflag)
      continue ;
    x = vdst->x;
    y = vdst->y;
    z = vdst->z;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  x0 = (xlo+xhi)/2.0f ; y0 = (ylo+yhi)/2.0f ; z0 = (zlo+zhi)/2.0f ;
  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->ripflag)
      continue ;
    vdst->x -= x0 ; vdst->y -= y0 ; vdst->z -= z0 ;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }

  mris_dst->xctr = mris_dst->yctr = mris_dst->zctr = 0 ;
  mris_dst->xlo = xlo ; mris_dst->ylo = ylo ; mris_dst->zlo = zlo ;
  mris_dst->xhi = xhi ; mris_dst->yhi = yhi ; mris_dst->zhi = zhi ;

  return(mris_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRIStalairachTransform(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int         vno ;
  VERTEX      *v ;
  Real        x, y, z, xt, yt, zt ;
  float       xlo, ylo, zlo, xhi, yhi, zhi ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  if (!mris_src->linear_transform)
    return(mris_dst) ;
#if 0
    ErrorReturn(mris_dst, 
                (ERROR_BADPARM, "MRIStalairachTransform: no xform loaded")) ;
#endif

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    v = &mris_dst->vertices[vno] ;
    x = v->x ; y = v->y ; z = v->z ;
    transform_point(mris_src->linear_transform, -x, z, y, &xt, &yt, &zt) ;
    v->x = -xt ; v->y = zt ; v->z = yt ;
    if (v->x > xhi) xhi = v->x;
    if (v->x < xlo) xlo = v->x;
    if (v->y > yhi) yhi = v->y;
    if (v->y < ylo) ylo = v->y;
    if (v->z > zhi) zhi = v->z;
    if (v->z < zlo) zlo = v->z;
  }

  mris_dst->xlo = xlo ; mris_dst->ylo = ylo ; mris_dst->zlo = zlo ;
  mris_dst->xctr = (xhi + xlo)/2 ;
  mris_dst->yctr = (yhi + ylo)/2 ;
  mris_dst->zctr = (zhi + zlo)/2 ;
  return(mris_dst) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisCountNegativeVertices(MRI_SURFACE *mris)
{
  int     vno, neg ;
  VERTEX  *v ;

  for (neg = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->neg)
      neg++ ;
  }

  return(neg) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISremoveNegativeVertices(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                           int min_neg, float min_neg_pct)
{
  int   t, niterations, write_iterations, neg, total_vertices, base_avgs ;
  float pct_neg, delta_t, scale, pct_neg_area, l_dist ;

  if (min_neg < 0)
    min_neg = 0 ;
  if (min_neg_pct < 0.0f)
    min_neg_pct = 0.0f ;

  if (Gdiag & DIAG_WRITE && parms->fp == NULL)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    parms->fp = fopen(fname, "w") ;
    if (!parms->fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
    Progname, fname) ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  parms->start_t = 0 ;
  mrisProjectSurface(mris) ;
  if (Gdiag & DIAG_WRITE)
    mrisLogStatus(mris, parms, parms->fp, 0.0f) ;
  mrisClearMomentum(mris) ;
  niterations = parms->niterations ;
  write_iterations = parms->write_iterations ;
  if (Gdiag & DIAG_WRITE && write_iterations > 0)
    mrisWriteSnapshot(mris, parms, 0) ;
  total_vertices = mrisValidVertices(mris) ;
  neg = mrisCountNegativeVertices(mris) ;
  pct_neg = (float)neg / (float)total_vertices ;
  l_dist = parms->l_dist ;
  pct_neg_area = 
    (float)mris->neg_area / (float)(mris->total_area+mris->neg_area) ;
  base_avgs = parms->n_averages ;
  for (t = 0 ; 
       (t < niterations) && (neg > min_neg) && (pct_neg_area > min_neg_pct) ; 
       t++)
  {
    if (pct_neg_area < 0.001)  /* hack!!, but it speeds things up */
    {
      parms->l_dist *= 1.1 ;
      if (parms->l_dist > 10*l_dist)
        parms->l_dist = 10*l_dist ;
      if (parms->l_dist > 1.0)
        parms->l_dist = 1.0 ;
    }
    if (pct_neg_area < 0.001)  /* another hack!!, but it speeds things up */
    {
      static int first = 1 ;
      /* don't want big steps or momentum for fine-scale stuff */
      parms->momentum = 0.0f ; 
      parms->dt = 0.1 ;
      if (Gdiag & DIAG_SHOW && first)
        fprintf(stdout, "setting momentum=%2.1f, dt=%2.1f, l_dist=%2.2f\n",
                parms->momentum, parms->dt, parms->l_dist) ;
      first = 0 ;
    }
    
    if (mris->patch)  /* area is constant so spring force doesn't decrease */
    {
      scale = sqrt(mris->orig_area / (mris->total_area+mris->neg_area)) ;
      MRISscaleBrain(mris, mris, scale) ;
      MRIScomputeMetricProperties(mris) ;
    }
    else
      mrisComputeVertexDistances(mris) ;
    MRISclearGradient(mris) ;      
    mrisComputeSpringTerm(mris, parms->l_spring) ;
    mrisComputeDistanceTerm(mris, parms) ;
    mrisComputeAngleAreaTerms(mris, parms) ;
/*    mrisAverageGradient(mris, parms->n_averages) ;*/

    switch (parms->integration_type)
    {
    case INTEGRATE_LM_SEARCH:
      delta_t = mrisLineMinimizeSearch(mris, parms) ;
      break ;
    default:
    case INTEGRATE_LINE_MINIMIZE:
      delta_t = mrisLineMinimize(mris, parms) ;
      break ;
    case INTEGRATE_MOMENTUM:
      delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, 
                                     parms->tol, parms->n_averages) ;
      break ;
    case INTEGRATE_ADAPTIVE:
      delta_t = mrisAdaptiveTimeStep(mris, parms);
      break ;
    }
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    neg = mrisCountNegativeVertices(mris) ;
    pct_neg = (float)neg / (float)total_vertices ;
    pct_neg_area = 
      (float)mris->neg_area / (float)(mris->total_area+mris->neg_area) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "%3.3d: count: %d (%2.2f%%), area: %2.2f (%2.2f%%)   \n",
              t, neg, 100.0f*pct_neg, mris->neg_area, 100.0f*pct_neg_area) ;
    if ((write_iterations > 0) &&!((t+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
      mrisWriteSnapshot(mris, parms, t+1) ;
    if (parms->n_averages == 0)
      parms->n_averages = base_avgs ;
    else
      parms->n_averages /= 2 ;
  }
  
  parms->n_averages = base_avgs ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "\n") ;
    if (Gdiag & DIAG_WRITE)
    {
      fclose(parms->fp) ;
      parms->fp = NULL ;
    }
  }
  mrisProjectSurface(mris) ;
  return(mris) ;
}
#endif
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE  *
MRISflatten(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     base_averages, n_averages, done, steps, total_steps ;
  float   base_tol ;
  double  starting_sse, ending_sse ;

  base_averages = parms->n_averages ;


  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s.out", parms->base_name) ;
    parms->fp = fopen(fname, "w") ;
    if (!parms->fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
    Progname, fname) ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;


  parms->start_t = 0 ;
  base_tol = parms->tol ;
  do
  {
    done = 0 ;
    mrisClearMomentum(mris) ;
    starting_sse = MRIScomputeSSE(mris, parms) ;
    for (total_steps = 0, n_averages = base_averages; !done ;n_averages /= 2)
    {
      steps = MRISintegrate(mris, parms, n_averages) ;
      parms->start_t += steps ;
      total_steps += steps ;
      done = n_averages == 0 ;   /* finished integrating at smallest scale */
    }
    parms->dt = parms->base_dt ;         /* reset time step */
    ending_sse = MRIScomputeSSE(mris, parms) ;
  } while (!FZERO(ending_sse) && ((starting_sse-ending_sse) > parms->tol)) ;


  if (Gdiag & DIAG_WRITE)
    fclose(parms->fp) ;

  mrisProjectSurface(mris) ;
  return(mris) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static float sigmas[] = { 4.0f, 2.0f, 1.0f, 0.5f } ;
#define NSIGMAS  (sizeof(sigmas)  / sizeof(sigmas[0]))

static char *surface_names[] = 
{
  "inflated",
  "smoothwm",
  "smoothwm"
} ;

static char *curvature_names[] = 
{
  NULL,
  "sulc",
  NULL
} ;

int
MRISsetOriginalFile(char *orig_name)
{
  surface_names[1] = surface_names[2] = orig_name ;
  return(NO_ERROR) ;
}


#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])
#define PARAM_IMAGES         (IMAGES_PER_SURFACE*SURFACES)

int
MRISregister(MRI_SURFACE *mris, MRI_SP *mrisp_template, 
             INTEGRATION_PARMS *parms, int max_passes, float min_degrees, float max_degrees, int nangles)
{
  float   sigma ;
  int     i, /*steps,*/ done, sno, ino, msec ;
  MRI_SP  *mrisp ;
  char    fname[STRLEN], base_name[STRLEN], path[STRLEN] ;
  double  base_dt ;
  struct  timeb start ;
  static  int first = 1 ;
  
  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  TimerStart(&start) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  FileNamePath(mris->fname, path) ;
  sprintf(base_name, "%s/%s.%s", path, 
          mris->hemisphere == LEFT_HEMISPHERE ? "lh":"rh", parms->base_name);

  base_dt = parms->dt ;
  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    if (!parms->start_t)      
    {
      parms->fp = fopen(fname, "w") ;
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
      Progname, fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris,parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris,parms) ;

  MRISuseMeanCurvature(mris) ;
  MRISnormalizeCurvature(mris) ;
  MRISstoreMeanCurvature(mris) ;

  if (parms->nbhd_size > 0)  /* compute long-range distances */
  {
    int i, nbrs[MAX_NBHD_SIZE] ;
    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      nbrs[i] = parms->max_nbrs ;
  }

  for (sno = 1 ; sno < SURFACES ; sno++)
  {
    if (!first && ((parms->flags & IP_USE_CURVATURE) == 0))
      break ;

    ino = parms->frame_no = sno*IMAGES_PER_SURFACE ;
    if (curvature_names[sno])  /* read in precomputed curvature file */
    {
      sprintf(fname, "%s.%s", 
              mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh", 
              curvature_names[sno]) ;
      if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read curvature file '%s'\n",
                  "MRISregister", fname) ;
      MRISnormalizeCurvature(mris) ;
    }
    else                       /* compute curvature of surface */
    {
      sprintf(fname, "%s", surface_names[sno]) ;
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  "MRISregister", fname) ;
      
      MRISsetNeighborhoodSize(mris, -1) ;  /* back to max */
      MRIScomputeMetricProperties(mris) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
      MRISnormalizeCurvature(mris) ;
      MRISresetNeighborhoodSize(mris,1);/*only use nearest neighbor distances*/
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    }
    MRISstoreMeanCurvature(mris) ;

    if (Gdiag & DIAG_SHOW)
    {
      if (curvature_names[sno])
        fprintf(stdout, "reading precomputed curvature from %s\n",fname) ;
      else
        fprintf(stdout, "calculating curvature of %s surface\n",fname) ;
    }

    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp,"calculating curvature of %s surface\n",fname);

    if (!first && parms->flags & IP_USE_CURVATURE)
    {
      /* only small adjustments needed after 1st time around */
      parms->tol *= 2.0f ;
      parms->l_corr /= 20.0f ;  /* should be more adaptive */
      if (Gdiag & DIAG_WRITE)
        mrisLogIntegrationParms(parms->fp, mris, parms) ;
      if (Gdiag & DIAG_SHOW)
        mrisLogIntegrationParms(stderr, mris, parms) ;
    }
    else
      if (!first) /* don't do curvature alignment */
        break ;   /* finished */

    for (i = 0 ; i < NSIGMAS ; i++)  /* for each spatial scale (blurring) */
    {
      parms->sigma = sigma = sigmas[i] ;
      parms->dt = base_dt ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "\nblurring surfaces with sigma=%2.2f...\n", sigma) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,"\ncorrelating surfaces with with sigma=%2.2f\n",
                sigma) ;
      if (Gdiag & DIAG_WRITE && !i && !parms->start_t)
      {
        MRISfromParameterization(mrisp_template, mris, ino);
				MRISnormalizeCurvature(mris) ;
        sprintf(fname, "%s/%s.target", path, mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh") ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "writing curvature file %s...\n", fname) ;
        MRISwriteCurvature(mris, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "done.\n") ;
      }
      MRISuseMeanCurvature(mris) ;
      mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
      parms->mrisp = MRISPblur(mrisp, NULL, sigma, 0) ;
      parms->mrisp_template = MRISPblur(mrisp_template, NULL, sigma, ino) ;
      MRISPblur(parms->mrisp_template, NULL, sigma, ino+1) ; /* variances */
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "done.\n") ;
      /* normalize curvature intensities for both source and target */
      MRISfromParameterization(parms->mrisp_template, mris, ino);
      MRISnormalizeCurvature(mris) ;
      MRIStoParameterization(mris, parms->mrisp_template, 1, ino) ;

#if 0
      /* normalize variances for both source and target */
      MRISfromParameterization(parms->mrisp_template, mris, ino+1);
      MRISnormalizeCurvature(mris) ;
      MRIStoParameterization(mris, parms->mrisp_template, 1, ino+1) ;
#endif

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        sprintf(fname, "%s/%s.%4.4dtarget%2.2f", 
                path, mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",
                parms->start_t, sigma) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "writing curvature file %s...", fname) ;
        MRISwriteCurvature(mris, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "done.\n") ;
      }
      
      MRISfromParameterization(parms->mrisp, mris, 0);
      MRISnormalizeCurvature(mris) ;
      MRIStoParameterization(mris, parms->mrisp, 1, 0) ;
      MRISPfree(&mrisp) ;

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        MRISPwrite(parms->mrisp, "mrisp_blur.hipl") ;
        MRISPwrite(parms->mrisp_template, "mrisp_template_blur.hipl") ;
      }
      mris->vp = (void *)parms->mrisp ;  /* hack to get it to projectSurface */

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        sprintf(fname, "%s/%s.%4.4dblur%2.2f", 
                path, mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",
                parms->start_t, sigma) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "writing curvature file %s...", fname) ;
        MRISwriteCurvature(mris, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "done.\n") ;
        sprintf(fname, "target.%s.%4.4d.hipl",parms->base_name,parms->start_t);
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "writing parameterization file %s...", fname) ;
        MRISPwrite(parms->mrisp_template, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "done.\n") ;
      }
      if (first)  /* only do rigid alignment first time through */
      {
        first = 0 ;
        if ((parms->flags & IP_NO_RIGID_ALIGN) == 0)
        {
          if (Gdiag & DIAG_SHOW)
            fprintf(stdout, "finding optimal rigid alignment\n") ;
          if (Gdiag & DIAG_WRITE)
            fprintf(parms->fp, "finding optimal rigid alignment\n") ;
          MRISrigidBodyAlignGlobal(mris, parms, min_degrees, max_degrees, nangles) ;
					/*          MRISrigidBodyAlignGlobal(mris, parms, 0.5f, 32.0f, 8) ;*/
           if (Gdiag & DIAG_WRITE && parms->write_iterations != 0)
             MRISwrite(mris, "rotated") ;
        }
      }

      mrisClearMomentum(mris) ;
      done = 0 ;
      mrisIntegrationEpoch(mris, parms, parms->n_averages) ;
    }
  }

  parms->tol /= 10 ;  /* remove everything possible pretty much */
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "removing remaining folds...\n") ;
  if (Gdiag & DIAG_WRITE)
    fprintf(parms->fp, "removing remaining folds...\n") ;
#if 1
  parms->l_nlarea *= 5 ;
  mrisIntegrationEpoch(mris, parms, parms->n_averages) ;
#else
  parms->l_nlarea = 1 ; parms->l_corr /= 10.0 ;
  parms->l_area = parms->l_parea = 0 ;
  mrisRemoveNegativeArea(mris,parms,parms->n_averages,MAX_NEG_AREA_PCT,3);
#endif
  MRISPfree(&parms->mrisp) ; MRISPfree(&parms->mrisp_template) ;
  msec = TimerStop(&start) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "registration took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
  if (Gdiag & DIAG_WRITE)
    fprintf(parms->fp, "registration took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static float area_coefs[] = { 1.0f,  1.0f, 0.1f } ;
static float dist_coefs[] = { 0.1f,  1.0f, 1.0f } ;

#define NCOEFS  sizeof(area_coefs) / sizeof(area_coefs[0])

#define MAX_NBHD_SIZE  200
#define NBR_COEF       (M_PI*1.0f)

MRI_SURFACE *
MRISunfold(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int max_passes)
{
  int     base_averages, i, nbrs[MAX_NBHD_SIZE], niter, passno, msec ;
  double  starting_sse, ending_sse, l_area, pct_error ;
  struct  timeb start ;
  
  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  TimerStart(&start) ;
  starting_sse = ending_sse = 0.0f ;   /* compiler warning */
  memset(nbrs, 0, MAX_NBHD_SIZE*sizeof(nbrs[0])) ;
#if 0
  if (mris->nsize < 2)
    nbrs[2] = nint(NBR_COEF*2.0) ;
  for (i = 4 ; i <= parms->nbhd_size ; i*= 2)
    nbrs[i] = nint(NBR_COEF*(float)i) ;
#else
  for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
    nbrs[i] = parms->max_nbrs ;
#endif
  
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    
    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
            parms->base_name);
    if (!parms->fp)
    {
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      else
        parms->fp = fopen(fname, "a") ;

      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "MRISunfold: could not open log file %s\n",
                  fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris,parms) ;
    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      if (nbrs[i])
        fprintf(parms->fp, "%d: %d | ", i, nbrs[i]) ;
    fprintf(parms->fp, "\n") ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      if (nbrs[i])
        fprintf(stdout, "%d: %d | ", i, nbrs[i]) ;
    fprintf(stdout, "\n") ;
    mrisLogIntegrationParms(stderr, mris, parms) ;
  }

  /*  parms->start_t = 0 ;*/
/*
   integrate until no improvement can be made at ANY scale, or until
   the error is effectively zero.
*/
  base_averages = parms->n_averages ;
  l_area = parms->l_area ;
  niter = parms->niterations ;
  passno = 0 ;
  do
  {
    if (mris->nsize < parms->nbhd_size)  /* resample distances on surface */
    {
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "resampling long-range distances") ;
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
      MRISsampleDistances(mris, nbrs, parms->nbhd_size) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      mrisClearMomentum(mris) ;
    }

    {
      char   *cp ;
      int    vno, n ;
      VERTEX *v, *vn ;
      float  d ;
      FILE   *fp ;
      
      cp = getenv("MEASURE_DISTANCES") ;
      if (cp)
      {
        fprintf(stdout, "outputting distance errors to distance.log...\n") ;
        fp = fopen("distance.log", "w") ;
        for (vno = 0 ; vno < mris->nvertices ; vno++)
        {
          v = &mris->vertices[vno] ;
          if (v->ripflag)
            continue ;
          for (n = 0 ; n < v->vtotal ; n++)
          {
            vn = &mris->vertices[v->v[n]] ;
            if (vn->ripflag)
              continue ;
            d = sqrt(SQR(vn->origx-v->origx)+SQR(vn->origy-v->origy)+
                     SQR(vn->origz-v->origz)) ;
            fprintf(fp, "%2.4f  %2.4f\n", d, v->dist_orig[n]) ;
          }
        }
        fclose(fp) ;
        exit(1) ;
      }
    }
    {
      char   *cp ;
      double max_pct ;

      cp = getenv("DISTURB_DISTANCES") ;
      if (cp)
      {
        max_pct = atof(cp) ;
        fprintf(stdout, "disturbing distances by %%%2.1f\n", (float)max_pct) ;
        if (Gdiag & DIAG_WRITE)
          fprintf(parms->fp, "disturbing distances by %%%2.1f\n", 
                  (float)max_pct);
        MRISdisturbOriginalDistances(mris, max_pct) ;
      }
    }

    {
      char   *cp ;

      cp = getenv("SPHERE") ;
      if (cp)
        MRISstoreAnalyticDistances(mris, MRIS_SPHERE) ;
      cp = getenv("PLANE") ;
      if (cp)
        MRISstoreAnalyticDistances(mris, MRIS_PLANE) ;
    }

    if (!passno && ((parms->flags & IPFLAG_QUICK) == 0))
    {
      double tol = parms->tol ;
      parms->tol = 0.5 ;
      if (niter > 30)
        parms->niterations = 30 ;
      mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 2);
      parms->niterations = niter ; parms->tol = tol ;
    }

    
    for (i = 0 ; i < NCOEFS ; i++)
    {
      if (mris->status == MRIS_SPHERE && i == NCOEFS-1)
        continue ;

      pct_error = MRISpercentDistanceError(mris) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp, 
                "pass %d: epoch %d of %d starting distance error %%%2.2f\n",
                passno, i+1, (int)(NCOEFS), (float)pct_error);
      fprintf(stdout, 
              "pass %d: epoch %d of %d starting distance error %%%2.2f\n",
              passno+1, i+1, (int)(NCOEFS), (float)pct_error);
    
      parms->l_dist = dist_coefs[i] ;
#if 1
      parms->l_area = area_coefs[i] ;
#else
      parms->l_nlarea = area_coefs[i] ;
#endif
      parms->l_angle = ANGLE_AREA_SCALE * parms->l_area ;
      if (i == NCOEFS-1)  /* see if distance alone can make things better */
        starting_sse = MRIScomputeSSE(mris, parms) ;
      mrisIntegrationEpoch(mris, parms, base_averages) ;
    }

#if 1
    parms->l_area = area_coefs[NCOEFS-1] ;
#else
    parms->l_nlarea = area_coefs[NCOEFS-1] ;
#endif
    parms->l_dist = dist_coefs[NCOEFS-1] ;
    ending_sse = MRIScomputeSSE(mris, parms) ;
    if (Gdiag & DIAG_SHOW)
    {
#if 0
      fprintf(stdout, "pass %d: start=%2.1f, end=%2.1f, ratio=%2.3f\n",
              passno+1, starting_sse, ending_sse, 
              (starting_sse-ending_sse)/starting_sse) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp, "pass %d: start=%2.4f, end=%2.4f, ratio=%2.4f\n",
                passno+1, starting_sse, ending_sse, 
                (starting_sse-ending_sse)/starting_sse) ;
#endif
    }
  } while (
           !FZERO(ending_sse) && 
           (((starting_sse-ending_sse)/starting_sse) > parms->tol) &&
           (++passno < max_passes)
           ) ;


  fprintf(stdout, "unfolding complete - removing small folds...\n") ;
  pct_error = MRISpercentDistanceError(mris) ;
  if (Gdiag & DIAG_WRITE)
    fprintf(parms->fp, 
            "starting distance error %%%2.2f\n", (float)pct_error);
  fprintf(stdout, 
          "starting distance error %%%2.2f\n", (float)pct_error);

  /* finally, remove all the small holes */
  parms->l_nlarea = 1.0f ; parms->l_area = 0.0 ;
  parms->l_dist = 0.1f ;  /* was 0.001 */
  parms->l_angle = ANGLE_AREA_SCALE * parms->l_nlarea ;
  parms->niterations = niter ;
#if 1
  parms->tol = 1e-2 ;  /* try and remove as much negative stuff as possible */
#else
  parms->tol = 1e-1 ;  /* try and remove as much negative stuff as possible */
#endif
  mrisStoreVtotalInV3num(mris) ;  /* hack to speed up neg. area removal */
  MRISresetNeighborhoodSize(mris, 1) ;
  fprintf(stdout, "removing remaining folds...\n") ;
  mrisRemoveNegativeArea(mris, parms, base_averages > 32 ? 32 : base_averages,
                         MAX_NEG_AREA_PCT, 3);
  MRISresetNeighborhoodSize(mris, 3) ;

  if (mris->status == MRIS_PLANE)  /* smooth out remaining folds */
  {
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "smoothing final surface...\n") ;
    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp, "smoothing final surface...\n") ;
    parms->l_spring = 1.0f ;
    parms->l_area = parms->l_nlarea = 0.0f ;
    parms->niterations = 5 ;
    parms->integration_type = INTEGRATE_MOMENTUM ;
    parms->dt = 0.5f ; parms->momentum = 0.0f ;
    parms->n_averages = 0 ;
    MRISintegrate(mris, parms, 0) ;
    /*    mrisRemoveNegativeArea(mris, parms, 0, MAX_NEG_AREA_PCT, 1);*/
  }
 

  pct_error = MRISpercentDistanceError(mris) ;
  fprintf(stdout, "final distance error %%%2.2f\n", (float)pct_error);
  mrisProjectSurface(mris) ;
  msec = TimerStop(&start) ;
  if (Gdiag & DIAG_SHOW)
  {
    mrisLogStatus(mris, parms, stderr, 0) ;
    fprintf(stdout, "optimization complete.\n") ;
    fprintf(stdout, "unfolding took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
  }
  if (Gdiag & DIAG_WRITE)
  {
    fprintf(parms->fp, "unfolding took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
    mrisLogStatus(mris, parms, parms->fp, 0) ;
    fprintf(parms->fp, "final distance error %%%2.2f\n", pct_error);
    fclose(parms->fp) ; parms->fp = NULL ;
  }

  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISquickSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int max_passes)
{
  int     niter, passno, msec, nbrs[MAX_NBHD_SIZE],i,use_dists, base_averages ;
  double  pct_error, orig_k, last_sse, sse, pct_change ;
  struct  timeb start ;

  orig_k = NEG_AREA_K ;

  TimerStart(&start) ;

  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;

  use_dists = (!FZERO(parms->l_dist) || !FZERO(parms->l_nldist)) &&
    (parms->nbhd_size > mris->nsize) ;

  memset(nbrs, 0, MAX_NBHD_SIZE*sizeof(nbrs[0])) ;
  for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
    nbrs[i] = parms->max_nbrs ;

  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    
    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
            parms->base_name);
    if (!parms->fp)
    {
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      else
        parms->fp = fopen(fname, "a") ;
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE,"MRISquickSphere: could not open log file %s\n",
                  fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris,parms) ;
    if (use_dists)
    {
      for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
        if (nbrs[i])
          fprintf(parms->fp, "%d: %d | ", i, nbrs[i]) ;
      fprintf(parms->fp, "\n") ;
    }
  }
  if (Gdiag & DIAG_SHOW)
  {
    if (use_dists)
    {
      for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
        if (nbrs[i])
          fprintf(stdout, "%d: %d | ", i, nbrs[i]) ;
      fprintf(stdout, "\n") ;
    }
    mrisLogIntegrationParms(stderr, mris, parms) ;
  }

  /* resample distances on surface */
  if (use_dists && mris->nsize < parms->nbhd_size)
  {
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "resampling long-range distances") ;
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISsampleDistances(mris, nbrs, parms->nbhd_size) ;
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    mrisClearMomentum(mris) ;
  }

/*
   integrate until no improvement can be made at ANY scale, or until
   the error has asymptoted.
*/
  base_averages = parms->n_averages ;
  parms->flags |= IP_RETRY_INTEGRATION ;
  niter = parms->niterations ;
  passno = 0 ;
#if 1
  if ((parms->flags & IPFLAG_QUICK) == 0)
    parms->tol = parms->tol * 1024 / (sqrt((double)base_averages+1)) ;
#endif	
for (i = 0, NEG_AREA_K = orig_k ; i < 4 ; NEG_AREA_K *= 4, i++)
  {
		passno = 0 ;
		do
		{
			last_sse = MRIScomputeSSE(mris, parms) ;
			printf("epoch %d (K=%2.1f), pass %d, starting sse = %2.2f\n",
						 i+1, NEG_AREA_K, passno+1, last_sse) ;
			niter = mrisIntegrationEpoch(mris, parms, base_averages) ;
			sse = MRIScomputeSSE(mris, parms) ;
			pct_change = (last_sse - sse) / (last_sse*niter) ; /* per timestep */
			passno++ ;
			printf("pass %d complete, delta sse/iter = %2.2f/%d = %2.2f\n",
						 passno, (last_sse-sse)/last_sse, niter, pct_change) ;
		} while (pct_change > parms->tol) ;
#if 0
		if (passno == 1)   /* couldn't make any progress at all */
			break ;
#endif
  } 

	NEG_AREA_K = orig_k ;
  pct_error = MRISpercentDistanceError(mris) ;
  fprintf(stdout, "final distance error %%%2.2f\n", (float)pct_error);
  mrisProjectSurface(mris) ;
  msec = TimerStop(&start) ;
  if (Gdiag & DIAG_SHOW)
  {
    mrisLogStatus(mris, parms, stderr, 0) ;
    fprintf(stdout, "optimization complete.\n") ;
    fprintf(stdout, "unfolding took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
  }
  if (Gdiag & DIAG_WRITE)
  {
    fprintf(parms->fp, "unfolding took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
    mrisLogStatus(mris, parms, parms->fp, 0) ;
    fprintf(parms->fp, "final distance error %%%2.2f\n", pct_error);
#if 0
    fclose(parms->fp) ;
#endif
  }

  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
#if 0
static float area_scoefs[] = 
    { 1.0f,    1.0f,   1.0f,  1.0f, 1.0f, 0.1f, 0.01f, 0.001f };
static float dist_scoefs[] = 
    { 0.0001f, 0.001f, 0.01f, 0.1f, 1.0f, 1.0f, 1.0f,  1.0f } ;
#else
static float area_scoefs[] = { 1.0f,  0.1f, 0.01f } ;
static float dist_scoefs[] = { 1.0f,  1.0f, 1.0f } ;
#endif

#define NSCOEFS  sizeof(area_scoefs) / sizeof(area_scoefs[0])


MRI_SURFACE *
MRISunfoldOnSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int max_passes)
{
  int     base_averages, i, nbrs[MAX_NBHD_SIZE], niter, passno, msec ;
  double  starting_sse, ending_sse ;
  struct  timeb start ;

  starting_sse = ending_sse = 0.0f ;   /* compiler warning */
  memset(nbrs, 0, MAX_NBHD_SIZE*sizeof(nbrs[0])) ;
  for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
    nbrs[i] = parms->max_nbrs ;
  
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    
    sprintf(fname, "%s.out", parms->base_name) ;
    if (!parms->fp)
    {
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      else
        parms->fp = fopen(fname, "a") ;
      if (!parms->fp)
  ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
      Progname, fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris,parms) ;
    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      if (nbrs[i])
        fprintf(parms->fp, "%d: %d | ", i, nbrs[i]) ;
    fprintf(parms->fp, "\n") ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      if (nbrs[i])
        fprintf(stdout, "%d: %d | ", i, nbrs[i]) ;
    fprintf(stdout, "\n") ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

/*
   integrate until no improvement can be made at ANY scale, or until
   the error has asymptoted.
*/
  base_averages = parms->n_averages ;
  niter = parms->niterations ;
  passno = 0 ;
  mrisProjectSurface(mris) ;
  MRIScomputeMetricProperties(mris) ;
  do
  {
    if (passno++ >= max_passes)
      break ;

    /* first time through only - use big ratio to remove folds */
    TimerStart(&start) ;
    for (i = 0 ; i < NSCOEFS ; i++)
    {
      if (mris->nsize < parms->nbhd_size)  /* resample distances on surface */
      {
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "resampling long-range distances") ;
        MRISsaveVertexPositions(mris, TMP_VERTICES) ;
        MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
        MRISsampleDistances(mris, nbrs, parms->nbhd_size) ;
        MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
        mrisClearMomentum(mris) ;
      }
      if (!i)  /* record starting error */
      {
        parms->l_nlarea = area_scoefs[NSCOEFS-1] ;
        parms->l_dist = dist_scoefs[NSCOEFS-1] ;
        starting_sse = MRIScomputeSSE(mris, parms) ;
      }

      /* remove any folds in the surface */
      mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 2) ;

      parms->l_dist = dist_scoefs[i] ;
      parms->l_nlarea = area_scoefs[i] ;
      parms->l_angle = ANGLE_AREA_SCALE * area_scoefs[i] ;
      mrisIntegrationEpoch(mris, parms, base_averages) ;
    }

    parms->l_area = area_scoefs[NSCOEFS-1] ;
    parms->l_dist = dist_scoefs[NSCOEFS-1] ;
    ending_sse = MRIScomputeSSE(mris, parms) ;
    if (Gdiag & DIAG_SHOW)
    {
#if 0
      fprintf(stdout, "pass %d: start=%2.1f, end=%2.1f, ratio=%2.3f\n",
              passno, starting_sse, ending_sse, 
              (starting_sse-ending_sse)/ending_sse) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp, "pass %d: start=%2.4f, end=%2.4f, ratio=%2.4f\n",
                passno, starting_sse, ending_sse, 
                (starting_sse-ending_sse)/starting_sse) ;
#endif
    }
    msec = TimerStop(&start) ;
    if (Gdiag & DIAG_SHOW)
    {
      fprintf(stdout,"epoch took %2.2f minutes\n",(float)msec/(1000.0f*60.0f));
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,
                "epoch took %2.2f minutes\n",(float)msec/(1000.0f*60.0f));
    }
  } while (!FZERO(ending_sse) && 
           (((starting_sse-ending_sse)/starting_sse) > parms->tol)) ;

#if 0
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, 
            "initial unfolding complete - settling to equilibrium...\n") ;
  parms->niterations = 1000 ;  /* let it go all the way to equilibrium */
  mrisIntegrationEpoch(mris, parms, parms->n_averages = 0) ; 
#endif

  /* finally, remove all the small holes */
  parms->l_nlarea = 1.0f ; parms->l_area = 0.0 ;
  parms->l_dist = 0.001f ;
  parms->l_angle = ANGLE_AREA_SCALE * area_scoefs[0] ;
  parms->niterations = niter ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "removing remaining folds...\n") ;
  mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 3);
  if (Gdiag & DIAG_SHOW)
    mrisLogStatus(mris, parms, stderr, 0) ;
  if (Gdiag & DIAG_WRITE)
    mrisLogStatus(mris, parms, parms->fp, 0) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "unfolding complete.\n") ;
  if (Gdiag & DIAG_WRITE)
  {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }

  mrisProjectSurface(mris) ;
  return(mris) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
static float neg_area_ratios[] = { 0.01f, 0.001f, 0.0001f } ;
#define MAX_PASSES (sizeof(neg_area_ratios) / sizeof(neg_area_ratios[0]))
static int
mrisRemoveNegativeArea(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, 
                       int base_averages, float min_area_pct, int max_passes)
{
  int    total_steps, done, steps, n_averages, old_averages, npasses, niter ;
  float  pct_neg, ratio ;
  char   *snum, *sdenom ;
  float  l_area, l_parea, l_corr, l_spring, l_dist, l_nlarea, 
         *pnum, *pdenom, cmod ;
  double tol ;

  if (Gdiag & DIAG_WRITE && parms->fp == NULL)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    if (!parms->start_t)
      parms->fp = fopen(fname, "w") ;
    else
      parms->fp = fopen(fname, "a") ;
    if (!parms->fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
    Progname, fname) ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  pct_neg = 100.0*mris->neg_area/(mris->neg_area+mris->total_area) ;
  if (pct_neg <= min_area_pct)
    return(0) ;   /* no steps */

  tol = parms->tol ;
#if 0
  parms->tol = 1e-2 ;
#endif
  niter = parms->niterations ;
  old_averages = parms->n_averages ;
  /*  parms->integration_type = INTEGRATE_LINE_MINIMIZE ;*/
  /* parms->niterations = 25 ;*/
  /*  base_averages = 1024 ;*/
  l_area = parms->l_area ; l_parea = parms->l_parea ; 
  l_nlarea = parms->l_nlarea ;
  l_spring = parms->l_spring ;
  l_dist = parms->l_dist ; l_corr = parms->l_corr ;
  parms->l_area = parms->l_parea = parms->l_dist = 
    parms->l_corr = parms->l_spring = parms->l_nlarea = 0.0 ;

  /* there is one negative area removing term (area, nlarea, parea, spring), 
     and one term we are seaking to retain (corr, dist).
     */
  cmod = 1.0f ;
  if (!FZERO(l_corr))
  { sdenom = "corr" ; pdenom = &parms->l_corr  ; /*cmod = 10.0f ;*/ }
  else
  { sdenom = "dist" ; pdenom = &parms->l_dist  ; }

  if (!FZERO(l_area))
  { snum = "area" ;   pnum = &parms->l_area ; }
  else if (!FZERO(l_parea))
#if 0
  { snum = "parea" ;  pnum = &parms->l_parea  ; }
#else
  { snum = "area" ;  pnum = &parms->l_area  ; }
#endif
  else if (!FZERO(l_nlarea))
  { snum = "nlarea" ; pnum = &parms->l_nlarea  ; }
  else
  { snum = "spring" ; pnum = &parms->l_spring  ; }

  npasses = 0 ;
  for (done=total_steps=0, n_averages = base_averages; !done ; n_averages /= 4)
  {
    *pnum = 1.0 ;
    *pdenom = cmod * (npasses >= MAX_PASSES ? 
        neg_area_ratios[MAX_PASSES-1] : neg_area_ratios[npasses]) ; 

    ratio = *pnum / *pdenom ;

    if (Gdiag & DIAG_SHOW && (n_averages == base_averages))
      fprintf(stdout, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    if (Gdiag & DIAG_WRITE && (n_averages == base_averages))
      fprintf(parms->fp, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    parms->n_averages = n_averages ;
    steps = MRISintegrate(mris, parms, n_averages) ;
    parms->start_t += steps ;
    total_steps += steps ;
    pct_neg = 100.0*mris->neg_area/(mris->neg_area+mris->total_area) ;
    if (pct_neg < min_area_pct)
      break ;
    done = n_averages == 0 ;    /* finished integrating at smallest scale */
    if (done && pct_neg > min_area_pct)  /* still too much negative area */
    {
      if (++npasses >= max_passes)
        break ;     /* can't get there from here (at least not quickly) */
      done = 0 ;     
      n_averages = base_averages*4 ;  /* try, try again */
    }
  }
#if 0
  MRIScomputeNormals(mris) ;
  mrisComputeVertexDistances(mris) ;
  MRIScomputeTriangleProperties(mris) ;  /* compute areas and normals */
  mrisOrientSurface(mris) ;
#endif
  parms->n_averages = old_averages  ;   /* hack, but no time to clean up now */
  parms->l_area = l_area ; parms->l_parea = l_parea ; 
  parms->l_spring = l_spring ;
  parms->l_dist = l_dist ; parms->l_corr = l_corr ;
  parms->niterations = niter ;
  parms->tol = tol ;
  return(total_steps) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
static int
mrisIntegrationEpoch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,int base_averages)
{
  int   total_steps, done, steps, n_averages, old_averages ;
  char  *snum, *sdenom ;
  float ratio, *pdenom, *pnum ;

  if (!FZERO(parms->l_corr))
  { sdenom = "corr" ; pdenom = &parms->l_corr  ; }
  else
  { sdenom = "dist" ; pdenom = &parms->l_dist  ; }

  if (!FZERO(parms->l_area))
  { snum = "area" ;   pnum = &parms->l_area ; }
  else if (!FZERO(parms->l_parea))
  { snum = "parea" ;  pnum = &parms->l_parea  ; }
  else if (!FZERO(parms->l_nlarea))
  { snum = "nlarea" ;  pnum = &parms->l_nlarea  ; }
  else
  { snum = "spring" ; pnum = &parms->l_spring  ; }

  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;
  if (Gdiag & DIAG_WRITE)
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  if (!FZERO(*pdenom))
  {
    ratio = *pnum / *pdenom ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN] ;
      if (!parms->fp)
      {
        sprintf(fname, "%s.%s.out", 
                mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",
                parms->base_name);
        if (!parms->start_t)
          parms->fp = fopen(fname, "w") ;
        else
          parms->fp = fopen(fname, "a") ;
  if (!parms->fp)
    ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
        Progname, fname) ;
      }
      fprintf(parms->fp, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    }
  }

  old_averages = parms->n_averages ;
  for (done = total_steps = 0, n_averages = base_averages ; !done ; 
       n_averages /= 4)
  {
    parms->n_averages = n_averages ;
    steps = MRISintegrate(mris, parms, n_averages) ;
    if (n_averages > 0 && parms->flags & IP_RETRY_INTEGRATION && 
        ((parms->integration_type == INTEGRATE_LINE_MINIMIZE) ||
        (parms->integration_type == INTEGRATE_LM_SEARCH)))
    {
      int niter = parms->niterations ;
      int integration_type = parms->integration_type ;

      fprintf(stdout, "taking momentum steps...\n") ;
      parms->integration_type = INTEGRATE_MOMENTUM ; parms->niterations = 10 ;
      parms->start_t += steps ;
      total_steps += steps ;
      steps = MRISintegrate(mris, parms, n_averages) ;
      parms->integration_type = integration_type ;
      parms->niterations = niter ;
      parms->start_t += steps ;
      total_steps += steps ;
      steps = MRISintegrate(mris, parms, n_averages) ;
    }
    parms->start_t += steps ;
    total_steps += steps ;
    done = n_averages == parms->min_averages ;
    if (mris->status == MRIS_SPHERE)
    {
      if (Gdiag & DIAG_SHOW)
        MRISprintTessellationStats(mris, stderr) ;
      parms->scale *= parms->dt_decrease ;
      if (parms->scale < 1.0f)
        parms->scale = 1.0f ;
    }
  }
#if 0
  MRIScomputeNormals(mris) ;
  mrisComputeVertexDistances(mris) ;
  MRIScomputeTriangleProperties(mris) ;  /* compute areas and normals */
  mrisOrientSurface(mris) ;
#endif
  parms->n_averages = old_averages  ;  /* hack, but no time to clean up now */
  return(total_steps) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
MRISintegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_averages)
{
  int     t, write_iterations, niterations, nsmall, neg ;
  double  l_dist, l_area, l_spring, sse, old_sse, delta_t, total_small = 0.0, 
          sse_thresh, pct_neg, pct_neg_area, total_vertices, tol
          /*, scale, last_neg_area */ ;
  MHT     *mht_v_current = NULL ;

  if (Gdiag & DIAG_WRITE && parms->fp == NULL)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    if (!parms->start_t)
      parms->fp = fopen(fname, "w") ;
    else
      parms->fp = fopen(fname, "a") ;
    if (!parms->fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
    Progname, fname) ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  l_spring = parms->l_spring ;
  l_dist = parms->l_dist ;
  l_area = parms->l_area ;
  write_iterations = parms->write_iterations ;
  niterations = parms->niterations ;
  if ((parms->flags & IPFLAG_QUICK) == 0)
    tol = parms->tol * sqrt(((double)n_averages + 1.0) / 1024.0);
  else
    tol = parms->tol ;
  sse_thresh = tol ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout,"integrating with navgs=%d and tol=%2.3e\n",n_averages,tol);
    
  mrisProjectSurface(mris) ;
  MRIScomputeMetricProperties(mris) ;

#if AVERAGE_AREAS
  MRISreadTriangleProperties(mris, mris->fname) ;
  mrisAverageAreas(mris, n_averages, ORIG_AREAS) ;
#endif

  parms->starting_sse = sse = old_sse = MRIScomputeSSE(mris, parms) ;

  delta_t = 0.0 ;
  niterations += parms->start_t ;
  parms->t = parms->start_t ;
  if (!parms->start_t)
  {
    mrisLogStatus(mris, parms, stderr, 0.0f) ;
    if (Gdiag & DIAG_WRITE)
    {
      mrisLogStatus(mris, parms, parms->fp, 0.0f) ;
      if (write_iterations > 0)
        mrisWriteSnapshot(mris, parms, 0) ;
    }
  }

  if (Gdiag_no >= 0)
    fprintf(stdout,
            "v %d curvature = %2.5f, position = (%2.3f,%2.3f,%2.3f)\n",
        Gdiag_no, mris->vertices[Gdiag_no].H,
        mris->vertices[Gdiag_no].x,
        mris->vertices[Gdiag_no].y,
        mris->vertices[Gdiag_no].x) ;
  total_small = 0.0 ; nsmall = 0 ;

  total_vertices = (double)mrisValidVertices(mris) ;
  neg = MRIScountNegativeTriangles(mris) ;
  pct_neg = (double)neg / total_vertices ;
  pct_neg_area = 
    (float)mris->neg_area / (float)(mris->total_area+mris->neg_area) ;
  /*  mrisClearMomentum(mris) ;*/
  for (parms->t = t = parms->start_t ; t < niterations ; t++)
  {
    if (!FZERO(parms->l_repulse_ratio))
      mht_v_current = 
        MHTfillVertexTableRes(mris, mht_v_current,CURRENT_VERTICES, 3.0);

    if (!FZERO(parms->l_curv))
      MRIScomputeSecondFundamentalForm(mris) ;

    MRISclearGradient(mris) ;      /* clear old deltas */
    mrisComputeDistanceTerm(mris, parms) ;
    mrisComputeAngleAreaTerms(mris, parms) ;
    mrisComputeCorrelationTerm(mris, parms) ;
    mrisComputePolarCorrelationTerm(mris, parms) ;
    /*    mrisComputeSpringTerm(mris, parms->l_spring) ;*/
#if 0
    mrisComputeCurvatureTerm(mris, parms) ;
    mrisComputeNegTerm(mris, parms) ;
    mrisComputeBoundaryTerm(mris, parms) ;
#endif
    mrisComputeRepulsiveRatioTerm(mris, parms->l_repulse_ratio, mht_v_current);

    mrisAverageGradients(mris, n_averages) ;
    mrisComputeSpringTerm(mris, parms->l_spring) ;
    mrisComputeTangentialSpringTerm(mris, parms->l_tspring) ;
    mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv) ;
    switch (parms->integration_type)
    {
    case INTEGRATE_LM_SEARCH:
      delta_t = mrisLineMinimizeSearch(mris, parms) ;
      break ;
    default:
    case INTEGRATE_LINE_MINIMIZE:
      delta_t = mrisLineMinimize(mris, parms) ;
      break ;
    case INTEGRATE_MOMENTUM:
      delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, 
                                     tol, parms->n_averages) ;
      break ;
    case INTEGRATE_ADAPTIVE:
      delta_t = mrisAdaptiveTimeStep(mris, parms);
      break ;
    }
      
    if (!FZERO(parms->l_pcorr) && (Gdiag & DIAG_SHOW))
    {
      float  alpha, beta, gamma ;

      alpha = DEGREES(delta_t*mris->alpha) ;
      beta = DEGREES(delta_t*mris->beta) ;
      gamma = DEGREES(delta_t*mris->gamma) ;
      fprintf(stdout, "rotating brain by (%2.1f, %2.1f, %2.1f)\n",
              alpha, beta, gamma) ;
    }
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    if (Gdiag_no >= 0)
      fprintf(stdout, "v %d curvature = %2.3f\n",
              Gdiag_no, mris->vertices[Gdiag_no].H) ;
    /* only print stuff out if we actually took a step */
    sse = MRIScomputeSSE(mris, parms) ;
    if (!FZERO(old_sse) && ((old_sse-sse)/(old_sse) < sse_thresh))
    {
      if (++nsmall > MAX_SMALL)
        break ;
      if (++total_small > TOTAL_SMALL)
        break ;
    }
    else
    {
      if (total_small > 0.0)  /* if error increases more than 1/4 time quit */
        total_small -= .25 ;
      nsmall = 0 ;
    }

    parms->t++ ;
    mrisLogStatus(mris, parms, stderr, delta_t) ;
    if (Gdiag & DIAG_WRITE)
      mrisLogStatus(mris, parms, parms->fp, delta_t) ;

    if ((Gdiag & DIAG_SHOW) && !((t+1)%10))
      MRISprintTessellationStats(mris, stderr) ;

    if ((write_iterations > 0) &&!((t+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
      mrisWriteSnapshot(mris, parms, t+1) ;
    if (mris->status == MRIS_PLANE && mris->neg_area > 4*mris->total_area)
    {
      fprintf(stdout, "flipping flattened patch...\n") ;
      mrisClearMomentum(mris) ;
      mrisFlipPatch(mris) ;
      MRIScomputeMetricProperties(mris) ;
    }

    if (FZERO(sse))
      break ;
    if ((parms->integration_type == INTEGRATE_LINE_MINIMIZE) ||
        (parms->integration_type == INTEGRATE_LM_SEARCH))
    {
      if ((100*(old_sse - sse) / sse) < tol)
        break ;
    }
    old_sse = sse ;
    if (FZERO(delta_t))   /* reached the minimum */
      break ;
  }

  if (!FZERO(parms->l_repulse))
    MHTfree(&mht_v_current) ;

  parms->ending_sse = MRIScomputeSSE(mris, parms) ;
  /*  mrisProjectSurface(mris) ;*/

  return(parms->t-parms->start_t) ;  /* return actual # of steps taken */
}
/*-----------------------------------------------------
        Parameters:


        Returns value:

        Description
          This routine is solely for reporting purposes - it is
          not used by any of the numerical integration routines.
------------------------------------------------------*/
static double
mrisComputeError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                 float *parea_rms, float *pangle_rms, 
                 float *pcurv_rms, float *pdist_rms,
                 float *pcorr_rms)
{
  double  rms, sse_area, sse_angle, sse_curv, delta, area_scale, sse_dist,
          sse_corr ;
  int     ano, fno, ntriangles, total_neighbors ;
  FACE    *face ;
  float   nv ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  sse_angle = sse_area = 0.0 ;
  for (ntriangles = fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    ntriangles++ ;
    delta = (double)(area_scale * face->area - face->orig_area) ;
    sse_area += delta*delta ;
    for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
    {
      delta = deltaAngle(face->angle[ano], face->orig_angle[ano]);
      sse_angle += delta*delta ;
    }
    if (!finite(sse_area) || !finite(sse_angle))
      ErrorExit(ERROR_BADPARM, "sse is not finite at face %d!\n",fno);
  }

  sse_corr = mrisComputeCorrelationError(mris, parms, 1) ;
  sse_dist = mrisComputeDistanceError(mris) ;
#if 0
  if (!FZERO(parms->l_spring))
    sse_spring = mrisComputeSpringEnergy(mris) ;
  if (!FZERO(parms->l_tspring))
    sse_tspring = mrisComputeTangentialSpringEnergy(mris) ;
#endif
  sse_curv = mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv) ;

  total_neighbors = mrisCountTotalNeighbors(mris) ;

  nv = (float)mrisValidVertices(mris) ;
  if  (mris->status != MRIS_PLANE)
    *pcurv_rms = (float)sqrt(sse_curv / nv) ;
  else
    *pcurv_rms = 0.0f ;
  *pdist_rms = (float)sqrt(sse_dist / (double)total_neighbors) ;
  *parea_rms = (float)sqrt(sse_area/(double)ntriangles) ;
  *pangle_rms =(float)sqrt(sse_angle/(double)(ntriangles*ANGLES_PER_TRIANGLE));
  *pcorr_rms = (float)sqrt(sse_corr / (double)nv) ;
  rms = MRIScomputeSSE(mris, parms) ;
#if 0
  rms = 
    *pdist_rms * parms->l_dist +
    *parea_rms * parms->l_area +
    *parea_rms * parms->l_parea +
    *pangle_rms * parms->l_angle +
    *pcorr_rms * (parms->l_corr+parms->l_pcorr) +
    sqrt(sse_spring/nv) * parms->l_spring ;
#endif
  return(rms) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          this is the routine actually used for the numerical integration.
          As such, it must represent the exact error function being
          minimized (as opposed to computeError above).
------------------------------------------------------*/
double
MRIScomputeSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double  sse, sse_area, sse_angle, delta, sse_curv, sse_spring, sse_dist,
          area_scale, sse_corr, sse_neg_area, l_corr, sse_val, sse_sphere,
          sse_grad, sse_nl_area, sse_tspring, sse_repulse, sse_tsmooth,
          sse_repulsive_ratio, sse_shrinkwrap;
  int     ano, fno ;
  FACE    *face ;
  MHT     *mht_v_current = NULL ;

#if METRIC_SCALE
  if (mris->patch || mris->noscale)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  sse_repulse = sse_nl_area = sse_corr = sse_angle = sse_neg_area = sse_val = sse_sphere =
    sse_shrinkwrap = sse_area = sse_spring = sse_curv = sse_dist = sse_tspring = sse_grad = 0.0;

  if (!FZERO(parms->l_repulse))
    mht_v_current = MHTfillVertexTable(mris, mht_v_current,CURRENT_VERTICES);

  if (!FZERO(parms->l_angle)||!FZERO(parms->l_area)||(!FZERO(parms->l_parea)))
  {
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      if (face->ripflag)
        continue ;

      delta = (double)(area_scale*face->area - face->orig_area) ;
#if ONLY_NEG_AREA_TERM
      if (face->area < 0.0f)
        sse_neg_area += delta*delta ;
#endif
      sse_area += delta*delta ;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      {
        delta = deltaAngle(face->angle[ano],face->orig_angle[ano]);
#if ONLY_NEG_AREA_TERM
        if (face->angle[ano] >= 0.0f)
          delta = 0.0f ;
#endif
        sse_angle += delta*delta ;
      }
      if (!finite(sse_area) || !finite(sse_angle))
        ErrorExit(ERROR_BADPARM, "sse not finite at face %d!\n",fno);
    }
  }
  if (parms->l_repulse > 0)
    sse_repulse = mrisComputeRepulsiveEnergy(mris, parms->l_repulse, mht_v_current) ;
  sse_repulsive_ratio = 
    mrisComputeRepulsiveRatioEnergy(mris, parms->l_repulse_ratio) ;
  sse_tsmooth = mrisComputeThicknessSmoothnessEnergy(mris, parms->l_tsmooth) ;
  if (!FZERO(parms->l_nlarea))
    sse_nl_area = mrisComputeNonlinearAreaSSE(mris) ;
  if (!FZERO(parms->l_nldist))
    sse_nl_area = mrisComputeNonlinearDistanceSSE(mris) ;
  if (!FZERO(parms->l_dist))
    sse_dist = mrisComputeDistanceError(mris) ;
  if (!FZERO(parms->l_spring))
    sse_spring = mrisComputeSpringEnergy(mris) ;
  if (!FZERO(parms->l_tspring))
    sse_tspring = mrisComputeTangentialSpringEnergy(mris) ;
  if (!FZERO(parms->l_curv))
    sse_curv = mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv) ;
  l_corr = (double)(parms->l_corr + parms->l_pcorr) ;
  if (!FZERO(l_corr))
    sse_corr = mrisComputeCorrelationError(mris, parms, 1) ;
  if (!FZERO(parms->l_intensity))
    sse_val = mrisComputeIntensityError(mris, parms) ;
  if (!FZERO(parms->l_grad))
    sse_grad = mrisComputeIntensityGradientError(mris, parms) ;
  if (!FZERO(parms->l_sphere))
    sse_sphere = mrisComputeSphereError(mris, parms->l_sphere, parms->a) ;
  if (!FZERO(parms->l_shrinkwrap))
    sse_shrinkwrap = mrisComputeShrinkwrapError(mris, parms->mri_brain, parms->l_shrinkwrap) ;

  sse = 0 ;

  sse += 
    (double)parms->l_area      * sse_neg_area + sse_repulse + sse_tsmooth +
    (double)parms->l_sphere    * sse_sphere + sse_repulsive_ratio +
    (double)parms->l_intensity * sse_val + 
    (double)parms->l_shrinkwrap * sse_shrinkwrap + 
    (double)parms->l_grad      * sse_grad + 
    (double)parms->l_parea     * sse_area + 
    (double)parms->l_nlarea    * sse_nl_area + 
    (double)parms->l_angle     * sse_angle + 
    (double)parms->l_dist      * sse_dist + 
    (double)parms->l_spring    * sse_spring + 
    (double)parms->l_tspring   * sse_tspring + 
    (double)l_corr             * sse_corr + 
    (double)parms->l_curv      * CURV_SCALE * sse_curv ;

  if (mht_v_current)
    MHTfree(&mht_v_current) ;

  return(sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          this is the routine actually used for the numerical integration.
          As such, it must represent the exact error function being
          minimized (as opposed to computeError above).
------------------------------------------------------*/
double
MRIScomputeSSEExternal(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                       double *ext_sse)
{
  double  sse ;

  if (gMRISexternalSSE)
    sse = (*gMRISexternalSSE)(mris, parms) ;
  else
    sse = 0 ;
  *ext_sse = sse ;
  sse += MRIScomputeSSE(mris, parms) ;

  return(sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeNonlinearAreaSSE(MRI_SURFACE *mris)
{
  double  sse, area_scale, error, ratio ;
  int     fno ;
  FACE    *face ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  sse = 0.0 ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
#define SCALE_NONLINEAR_AREA 0
#if SCALE_NONLINEAR_AREA
    if (!FZERO(face->orig_area))
      ratio = area_scale*face->area / face->orig_area ;
    else
      ratio = 0.0f ;
#else
    ratio = area_scale*face->area ;
#endif
    if (ratio > MAX_NEG_RATIO)
      ratio = MAX_NEG_RATIO ;
    else if (ratio < -MAX_NEG_RATIO)
      ratio = -MAX_NEG_RATIO ;
#if 0
    error = (1.0 / NEG_AREA_K) * log(1.0+exp(-NEG_AREA_K*ratio)) ;
#else
    error = (log(1.0+exp(NEG_AREA_K*ratio)) / NEG_AREA_K) - ratio ;
#endif

    sse += error ;
    if (!finite(sse) || !finite(error))
      ErrorExit(ERROR_BADPARM, "nlin area sse not finite at face %d!\n",fno);
  }
  return(sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISscaleBrainArea(MRI_SURFACE *mris)
{
  float   scale ;

  scale = sqrt(mris->orig_area / (mris->total_area+mris->neg_area)) ;
  MRISscaleBrain(mris, mris, scale) ;
  MRIScomputeMetricProperties(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Find the 3 poles (temporal, occipital, and frontal) of
          the cortical surface.
------------------------------------------------------*/
#if 0
#define MIN_Z_DISTANCE       30.0f
#define MIN_Y_DISTANCE       30.0f
#define MIN_ANGLE_VARIATION  RADIANS(30.0f)

static int
mrisFindPoles(MRIS *mris)
{
  int     vno, n, neigh, local_max ;
  VERTEX  *vertex, *v_neigh ;
  Real    x, y, z, xt, yt, zt ;
  float   temporal_y_hi = -1000, zfront, yfront, angle ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "(%2.0f, %2.0f, %2.0f) --> (%2.0f, %2.0f, %2.0f), ctr "
            "(%2.0f, %2.0f, %2.0f)\n",
            mris->xlo, mris->ylo, mris->zlo, mris->xhi, mris->yhi, mris->zhi,
            mris->xctr, mris->yctr, mris->zctr);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "finding cortical poles...") ;

  /* first find frontal and occipital poles */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->y >= mris->yhi)
      mris->v_frontal_pole = vertex ;
    else if (vertex->y <= mris->ylo)
      mris->v_occipital_pole = vertex ;
  }

  /* now find temporal pole */
  zfront = mris->v_frontal_pole->z ; yfront = mris->v_frontal_pole->y ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = (Real)vertex->x ; y = (Real)vertex->y ; z = (Real)vertex->z ;
    if (mris->linear_transform)
      transform_point(mris->linear_transform, -x, z, y, &xt, &yt, &zt) ;
    else
    { xt = -x ; yt = z ; zt = y ; }

/*
  some rules for finding the temporal pole:

  1. must be a local max in y (posterior-anterior) direction.
  2. Must be below some absolute y talairach coordinate.
  3. Must be a minimum distance from the frontal pole in y and z directions.
  4. Must have a normal vector within 30 degrees of (0,1,0).
*/
    if ((yt < MAX_TALAIRACH_Y) &&  /* low enough to be temporal pole */
        ((zfront - vertex->z) > MIN_Z_DISTANCE) &&
        ((yfront - vertex->y) > MIN_Y_DISTANCE))
    {
      local_max = 1 ;
      if (vertex->y > temporal_y_hi)  /* check neighbors positions */
      {
        for (n = 0 ; n < vertex->vnum ; n++)
        {
          neigh = vertex->v[n] ;
          v_neigh = &mris->vertices[neigh] ;
          if (v_neigh->y > vertex->y)
          {
            local_max = 0 ;
            break ;
          }
        }

        angle = acos(vertex->ny) ;
        if (local_max && (angle < MIN_ANGLE_VARIATION))
        {
          mris->v_temporal_pole = vertex ;
          temporal_y_hi = vertex->y ;
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    if (mris->v_temporal_pole)
      fprintf(stdout, "F: (%2.0f,%2.0f,%2.0f), T: (%2.0f,%2.0f,%2.0f) "
              "O: (%2.0f,%2.0f,%2.0f).\n",
              mris->v_frontal_pole->x, mris->v_frontal_pole->y, 
              mris->v_frontal_pole->z,
              mris->v_temporal_pole->x, mris->v_temporal_pole->y, 
              mris->v_temporal_pole->z,
              mris->v_occipital_pole->x, mris->v_occipital_pole->y, 
              mris->v_occipital_pole->z) ;
    else
      fprintf(stdout, "F: (%2.0f,%2.0f,%2.0f), T: (NOT FOUND), "
              "O: (%2.0f,%2.0f,%2.0f).\n",
              mris->v_frontal_pole->x, mris->v_frontal_pole->y, 
              mris->v_frontal_pole->z,
              mris->v_occipital_pole->x, mris->v_occipital_pole->y, 
              mris->v_occipital_pole->z) ;
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisOrientEllipsoid(MRI_SURFACE *mris)
{
  int     fno, ano ;
  VERTEX  *v ;
  FACE    *face ;
  float   dot ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
      
    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
       */
    v = &mris->vertices[face->v[0]] ;
    dot = v->x * face->nx + v->y * face->ny + v->z * face->nz;
    if (dot < 0.0f)   /* not in same direction, area < 0 and reverse n */
    {
      face->area *= -1.0f ;
      face->nx *= -1.0f; face->ny *= -1.0f; face->nz *= -1.0f;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        face->angle[ano] *= -1.0f ;
    }
  }

  /* now recompute the total surface area, ignoring negative areas */
#if 0
  if ((mris->status != MRIS_PARAMETERIZED_SPHERE) || (!mris->total_area))
#endif
  {
    mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f ;
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      if (face->ripflag)
        continue ;
      if (face->area >= 0.0f)
        mris->total_area += face->area ;
      else
      {
        mris->neg_area += -face->area ;
        mris->neg_orig_area += face->orig_area ;
      }
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisOrientPlane(MRI_SURFACE *mris)
{
  int     fno, ano, vno ;
  FACE    *face ;
  VERTEX  *v ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
      
    /* now give the area an orientation: if the unit normal is pointing
       downwards in the plane then the area should be negative.
       */
    if (face->nz < 0.0f)   
    {
      /* not in same direction, area < 0 and reverse n */
      face->area *= -1.0f ;
      face->nx *= -1.0f; face->ny *= -1.0f; face->nz *= -1.0f;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        face->angle[ano] *= -1.0f ;
    }
  }

  /* now recompute the total surface area, ignoring negative areas */
  mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    if (face->area >= 0.0f)
      mris->total_area += face->area ;
    else
    {
      mris->neg_area += -face->area ;
      mris->neg_orig_area += face->orig_area ;
    }
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    if (v->nz < 0)
    {
      v->nz *= -1 ;
      v->neg = 1 ;
    }
    else
      v->neg = 0 ;
    v->area = 0 ;
    for (fno = 0 ; fno < v->num ; fno++)
    {
      face = &mris->faces[v->f[fno]] ;
      v->area += face->area ;
    }
    v->area /= 2 ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScountNegativeTriangles(MRI_SURFACE *mris)
{
  int     fno, negative ;
  FACE    *face ;

  negative = 0 ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    if (face->area < 0.0f)
      negative++ ;
  }

  return(negative) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISupdateEllipsoidSurface(MRI_SURFACE *mris)
{
#if 0
  MRIScomputeTriangleProperties(mris) ;  /* recompute areas and normals */
#endif
  mrisOrientEllipsoid(mris) ;      /* orient the normals and angles */
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisAverageGradients(MRI_SURFACE *mris, int num_avgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  dx, dy, dz, num, sigma ;
  VERTEX *v, *vn ;
  MRI_SP *mrisp, *mrisp_blur ;

  if (num_avgs <= 0)
    return(NO_ERROR) ;

  if (Gdiag_no >= 0)
  {
    v = &mris->vertices[Gdiag_no] ;
    fprintf(stdout, "before averaging dot = %2.2f ",
            v->dx*v->nx+v->dy*v->ny+v->dz*v->nz) ;
  }
  if (0 && mris->status == MRIS_PARAMETERIZED_SPHERE)  /* use convolution */
  {
    sigma = sqrt((float)num_avgs) / 4.0 ;
    mrisp = MRISgradientToParameterization(mris, NULL, 1.0) ;
    mrisp_blur = MRISPblur(mrisp, NULL, sigma, -1) ;
    MRISgradientFromParameterization(mrisp_blur, mris) ;
    MRISPfree(&mrisp) ; MRISPfree(&mrisp_blur) ;
  }
  else for (i = 0 ; i < num_avgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      dx = v->dx ; dy = v->dy ; dz = v->dz ;
      pnb = v->v ;
      /*      vnum = v->v2num ? v->v2num : v->vnum ;*/
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        dx += vn->dx ; dy += vn->dy ; dz += vn->dz ;
#if 0
        if (vno == Gdiag_no)
        {
          float dot ;
          dot = vn->dx*v->dx + vn->dy*v->dy + vn->dz*v->dz ;
          if (dot < 0)
            fprintf(stdout, "vn %d: dot = %2.3f, dx = (%2.3f, %2.3f, %2.3f)\n",
                    v->v[vnb], dot, vn->dx, vn->dy, vn->dz) ;
        }
#endif
      }
      num++ ;
      v->tdx = dx / num ;
      v->tdy = dy / num ;
      v->tdz = dz / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->dx = v->tdx ; v->dy = v->tdy ; v->dz = v->tdz ;
    }
  }
  if (Gdiag_no >= 0)
  {
    float dot ;
    v = &mris->vertices[Gdiag_no] ;
    dot = v->nx*v->dx + v->ny*v->dy + v->nz*v->dz ;
    fprintf(stdout, " after dot = %2.2f\n",dot) ;
    if (fabs(dot) > 50)
      DiagBreak() ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadTriangleProperties(MRI_SURFACE *mris, char *mris_fname)
{
  int     ano, vnum,fnum, fno, vno ;
  FACE    *face ;
  VERTEX  *v ;
  float   f;
  FILE    *fp;
  char    fname[STRLEN], fpref[STRLEN], hemi[20], *cp ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "reading triangle files...") ;

  cp = strrchr(mris_fname, '.') ;
  if (!cp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM,"MRISreadTriangleProperties(%s): could not scan"
                 " hemisphere from fname", mris_fname)) ;
  strncpy(hemi, cp-2, 2) ; hemi[2] = 0 ;
  FileNamePath(mris_fname, fpref) ;

  sprintf(fname, "%s/%s.triangle_area", fpref, hemi) ;

  fp = fopen(fname,"r");
  if (fp==NULL)   
  {
    fprintf(stdout, 
            "\nno precomputed triangle areas and angles - computing...\n");
    return(1) ;  /* doesn't exist */
  }

  fread4((float *)&vnum,fp);
  fread4((float *)&fnum,fp);
  if (vnum!=mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE, "MRISreadTriangleProperties: incompatible vertex "
                 "number in file %s", fname)) ;
  }
  
  mris->orig_area = 0.0f ;
  for (fno=0;fno<fnum;fno++)
  {
    face = &mris->faces[fno] ;
    f = freadFloat(fp);
    face->orig_area = f ;
    mris->orig_area += f;
  }

  /* compute original vertex areas from faces */
  for (vno=0;vno<vnum;vno++)
  {
    v = &mris->vertices[vno] ;
    v->origarea = 0.0f ;
    for (fno = 0 ; fno < v->num ; fno++)
    {
      face = &mris->faces[v->f[fno]] ;
      v->origarea += face->orig_area ;
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "total area = %2.0f.\n", mris->orig_area) ;


  /* now open and read the angle file */
  sprintf(fname, "%s/%s.triangle_angle", fpref, hemi) ;
  fp = fopen(fname,"r");
  if (fp==NULL)   
    return(1) ;  /* doesn't exist */

  fread4((float *)&vnum,fp);
  fread4((float *)&fnum,fp);
  if (vnum!=mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE, "MRISreadTriangleProperties: incompatible vertex "
                 "number in file %s", fname)) ;
  }
  
  for (fno=0;fno<fnum;fno++)
  {
    face = &mris->faces[fno] ;
    for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
    {
      f = freadFloat(fp);
      face->orig_angle[ano] = f ;
    }
  }

#if 0
  /* read in the distances to all neighboring vertices */
  sprintf(fname, "%s/%s.dist", fpref, hemi) ;
  fp = fopen(fname,"rb");
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, 
                              "MRISreadTriangleProperties: could not open %s",
                                fname)) ;

  fread4((float *)&vnum,fp);
  fread4((float *)&fnum,fp);
  for (vno = 0; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->vtotal ; n++)
      v->dist_orig[n] = freadFloat(fp) ;
  }

  fclose(fp);
#endif

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteTriangleProperties(MRI_SURFACE *mris, char *mris_fname)
{
  int     fno, ano, vno ;
  FACE    *face ;
  FILE    *fp;
  char    fname[STRLEN], fpref[STRLEN], hemi[20], *cp ;

  MRIScomputeTriangleProperties(mris) ;

  cp = strrchr(mris_fname, '.') ;
  if (!cp)
    ErrorReturn(ERROR_BADPARM, 
              (ERROR_BADPARM, "MRISwriteTriangleProperties(%s): could not scan"
                 "hemisphere from fname", mris_fname)) ;
  strncpy(hemi, cp-2, 2) ; hemi[2] = 0 ;
  FileNamePath(mris_fname, fpref) ;


  sprintf(fname, "%s/%s.triangle_area", fpref, hemi) ;
  fp = fopen(fname,"wb");
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, 
                              "MRISwriteTriangleProperties: could not open %s",
                                fname)) ;

  /* write out the distances to all neighboring vertices */
  fwrite4(mris->nvertices,fp);
  fwrite4(mris->nfaces,fp);
  for (fno=0;fno<mris->nfaces;fno++)
  {
    face = &mris->faces[fno] ;
    fwriteFloat(face->area, fp) ;
  }

  fclose(fp);


  sprintf(fname, "%s/%s.triangle_angle", fpref, hemi) ;
  fp = fopen(fname,"wb");
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, 
                              "MRISwriteTriangleProperties: could not open %s",
                                fname)) ;

  /* write out the area of all the triangles */
  fwrite4(mris->nvertices,fp);
  fwrite4(mris->nfaces,fp);
  for (vno=0;vno<mris->nfaces;vno++)
  {
    face = &mris->faces[fno] ;
    for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      fwriteFloat(face->angle[ano], fp) ;
  }

  fclose(fp);

#if 0
  /* write out the distances to all neighboring vertices */
  sprintf(fname, "%s/%s.dist", fpref, hemi) ;
  fp = fopen(fname,"wb");
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, 
                              "MRISwriteTriangleProperties: could not open %s",
                                fname)) ;

  fwrite4(mris->nvertices,fp);
  fwrite4(mris->nfaces,fp);
  for (vno = 0; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->v2num ; n++)
      fwriteFloat(v->dist[n], fp) ;
  }

  fclose(fp);
#endif

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        each face has 2 triangles defined by it:

       V0       d      V3
        o--------------o
        |              |
        | A0           |
      a |              | c
        |              |
        |           A1 |
        o--------------o
       V1      b        V2        

       a = V1 - V0
       d = V3 - V0
       e = V3 - V1
       A0 = 0.5 (a x d) . n

       b = V1 - V2
       c = V3 - V2
       A1 = 0.5 (c x b) . n


        each face has 1 triangle defined by it:
       V0    b     V2
        o----------o
        |         /      
        | A0    /        
      a |     /          
        |   /            
        | /
        o
       V1

       a = V1 - V0
       b = V2 - V0
       A0 = 0.5 (a x b) . n

       
------------------------------------------------------*/
int
MRIScomputeTriangleProperties(MRI_SURFACE *mris)
{
  VECTOR  *v_a, *v_b, *v_n ;
  VERTEX  *v0, *v1, *v2, *va, *vb, *vo, *v ;
  FACE    *face ;
  int     fno, ano, vno  ;
  float   area, angle, dot, cross, dz ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */

  mris->total_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;
    VERTEX_EDGE(v_a, v0, v1) ;  
    VERTEX_EDGE(v_b, v0, v2) ;

    /* compute metric properties of first triangle */
    V3_CROSS_PRODUCT(v_a, v_b, v_n) ;
    area = V3_LEN(v_n) * 0.5f ;
    dot = V3_DOT(v_a, v_b) ;
    face->area = area ;
    V3_NORMALIZE(v_n, v_n) ;             /* make it a unit vector */
    face->nx = V3_X(v_n); face->ny = V3_Y(v_n); face->nz = V3_Z(v_n);
    mris->total_area += area ;

    /* now compute angles */
    VECTOR_LOAD(v_n, face->nx, face->ny, face->nz) ;
    if ((V3_X(v_n) < V3_Y(v_n)) && (V3_X(v_n) < V3_Z(v_n)))
      dz = fabs(V3_X(v_n)) ;
    else if (V3_Y(v_n) < V3_Z(v_n))
      dz = fabs(V3_Y(v_n)) ;
    else
      dz = fabs(V3_Z(v_n)) ;
    for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
    {
      switch (ano)   /* vertices for triangle 1 */
      {
      default:
      case 0: vo = v0 ; va = v2 ; vb = v1 ; break ;
      case 1: vo = v1 ; va = v0 ; vb = v2 ; break ;
      case 2: vo = v2 ; va = v1 ; vb = v0 ; break ;
      }
      
      VERTEX_EDGE(v_a, vo, va) ;VERTEX_EDGE(v_b, vo, vb) ;
      cross = VectorTripleProduct(v_b, v_a, v_n) ;
      dot = V3_DOT(v_a, v_b) ;
      angle = atan2(cross, dot) ;
      face->angle[ano] = angle ;
      
#if 0
      if (angle < 0.0f || angle >= M_PI)
        fprintf(stdout, "angle [%d][%d] = %2.1f\n",
                fno,ano,(float)DEGREES(angle)) ;
#endif
    }
  }

  /* calculate the "area" of the vertices */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->area = 0.0 ;
    for (fno = 0 ; fno < v->num ; fno++)
      v->area += mris->faces[v->f[fno]].area ;
    v->area /= 2.0 ;
  }

  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_n) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
#define MAX_MM           MAX_DIM/10.0f
#else
#define MAX_MM           MAX_DIM/10.0f
#endif
#define MAX_PLANE_MM     100*10.0f
#define MAX_MOMENTUM_MM  1
#define MIN_MM           0.001

static double
mrisAdaptiveTimeStep(MRI_SURFACE *mris,INTEGRATION_PARMS *parms)
{
  double  delta_t, sse, starting_sse ;

  starting_sse = MRIScomputeSSE(mris, parms) ;

  MRISstoreCurrentPositions(mris) ;
  delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, parms->tol,
                                 parms->n_averages) ;

  sse = MRIScomputeSSE(mris, parms) ;

  if (sse > starting_sse)  /* error increased - turn off momentum */
  {
    mrisClearMomentum(mris) ;
    parms->dt *= parms->dt_decrease ;
    if (parms->dt <= parms->base_dt)
      parms->dt = parms->base_dt ;
      
    if (sse / starting_sse > parms->error_ratio)  /* undo time step */
    {
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "sse increased by %2.0f%%, undoing time step...\n",
                (float)sse/starting_sse * 100.0f) ;
      if (parms->dt > parms->base_dt)    /* reset step size */
        parms->dt = parms->base_dt ;

      /* undo the time step */
      MRISrestoreOldPositions(mris) ;
      mrisProjectSurface(mris) ;
    }
  }
  else   /* error decreased */
    parms->dt *= parms->dt_increase ;

  return(delta_t) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double  
mrisAsynchronousTimeStep(MRI_SURFACE *mris, float momentum, 
                                    float delta_t, MHT *mht, float max_mag)
{
  static int direction = 1 ;
  double  mag ;
  int     vno, i ;
  VERTEX  *v ;

  /* take a step in the gradient direction modulated by momentum */
  if (mris->status == MRIS_RIGID_BODY)
  {
    mris->da = delta_t * mris->alpha + momentum * mris->da ;
    mris->db = delta_t * mris->beta + momentum * mris->db ;
    mris->dg = delta_t * mris->gamma + momentum * mris->dg ;
    MRISrotate(mris, mris, mris->da, mris->db, mris->dg) ;
  }
  else for (i = 0 ; i < mris->nvertices ; i++)
  {
    if (direction < 0)
      vno = mris->nvertices - i - 1 ;
    else
      vno = i ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->odx = delta_t * v->dx + momentum*v->odx ;
    v->ody = delta_t * v->dy + momentum*v->ody ;
    v->odz = delta_t * v->dz + momentum*v->odz ;
    mag = sqrt(v->odx*v->odx + v->ody*v->ody + v->odz*v->odz) ;
    if (mag > max_mag) /* don't let step get too big */
    {
      mag = max_mag / mag ;
      v->odx *= mag ; v->ody *= mag ; v->odz *= mag ;
    }
    if (vno == Gdiag_no)
    {
      float dist, dot, dx, dy, dz ;

      dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
      dist = sqrt(dx*dx+dy*dy+dz*dz) ;
      dot = dx*v->nx + dy*v->ny + dz*v->nz ;
      fprintf(stdout, "moving v %d by (%2.2f, %2.2f, %2.2f) dot=%2.2f-->"
              "(%2.1f, %2.1f, %2.1f)\n", vno, v->odx, v->ody, v->odz, 
              v->odx*v->nx+v->ody*v->ny+v->odz*v->nz,
              v->x, v->y, v->z) ;
      fprintf(stdout, "n = (%2.1f,%2.1f,%2.1f), total dist=%2.3f, total dot = %2.3f\n", 
              v->nx, v->ny, v->nz, dist, dot) ;
    }

    /* erase the faces this vertex is part of */
#if 0
    for (fno = 0 ; fno < v->num ; fno++)
      mrisEraseFace(mris, mri_filled, v->f[fno]) ;
#else
    if (mht)
      MHTremoveAllFaces(mht, mris, v) ;
#endif

    if (mht)
      mrisLimitGradientDistance(mris, mht, vno) ;

    v->x += v->odx ; 
    v->y += v->ody ;
    v->z += v->odz ;

    if ((fabs(v->x) > 128.0f) ||
        (fabs(v->y) > 128.0f) ||
        (fabs(v->z) > 128.0f))
      DiagBreak() ;

    /* should this be done here????? (BRF) what about undoing step??? */
    v->dx = v->odx ;  /* for mrisTrackTotalDistances */
    v->dy = v->ody ;
    v->dz = v->odz ;

#if 0
    /* write the new face positions into the filled volume */
    for (fno = 0 ; fno < v->num ; fno++)
      mrisFillFace(mris, mri_filled, v->f[fno]) ;
#else
    if (mht)
      MHTaddAllFaces(mht, mris, v) ;
#endif

  }

  direction *= -1 ;
  return(delta_t) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double  
mrisAsynchronousTimeStepNew(MRI_SURFACE *mris, float momentum, 
                            float delta_t, MHT *mht, float max_mag)
{
  static int direction = 1 ;
  double  mag ;
  int     vno, i ;
  VERTEX  *v ;

  /* take a step in the gradient direction modulated by momentum */
  if (mris->status == MRIS_RIGID_BODY)
  {
    mris->da = delta_t * mris->alpha + momentum * mris->da ;
    mris->db = delta_t * mris->beta + momentum * mris->db ;
    mris->dg = delta_t * mris->gamma + momentum * mris->dg ;
    MRISrotate(mris, mris, mris->da, mris->db, mris->dg) ;
  }
  else for (i = 0 ; i < mris->nvertices ; i++)
  {
    if (direction < 0)
      vno = mris->nvertices - i - 1 ;
    else
      vno = i ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->odx = delta_t * v->dx + momentum*v->odx ;
    v->ody = delta_t * v->dy + momentum*v->ody ;
    v->odz = delta_t * v->dz + momentum*v->odz ;
    mag = sqrt(v->odx*v->odx + v->ody*v->ody + v->odz*v->odz) ;
    if (mag > max_mag) /* don't let step get too big */
    {
      mag = max_mag / mag ;
      v->odx *= mag ; v->ody *= mag ; v->odz *= mag ;
    }
    if (vno == Gdiag_no)
    {
      float dist, dot, dx, dy, dz ;

      dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
      dist = sqrt(dx*dx+dy*dy+dz*dz) ;
      dot = dx*v->nx + dy*v->ny + dz*v->nz ;
      fprintf(stdout, "moving v %d by (%2.2f, %2.2f, %2.2f) dot=%2.2f-->"
              "(%2.1f, %2.1f, %2.1f)\n", vno, v->odx, v->ody, v->odz, 
              v->odx*v->nx+v->ody*v->ny+v->odz*v->nz,
              v->x, v->y, v->z) ;
#if 0
      fprintf(stdout, "n = (%2.1f,%2.1f,%2.1f), total dist=%2.3f, total dot = %2.3f\n", 
              v->nx, v->ny, v->nz, dist, dot) ;
#endif
    }

    /* erase the faces this vertex is part of */
#if 0
    for (fno = 0 ; fno < v->num ; fno++)
      mrisEraseFace(mris, mri_filled, v->f[fno]) ;
#else
    if (mht)
      MHTremoveAllFaces(mht, mris, v) ;
#endif

    if (mht)
      mrisLimitGradientDistance(mris, mht, vno) ;

    v->x += v->odx ; 
    v->y += v->ody ;
    v->z += v->odz ;

    if ((fabs(v->x) > 128.0f) ||
        (fabs(v->y) > 128.0f) ||
        (fabs(v->z) > 128.0f))
      DiagBreak() ;

#if 0
    /* should this be done here????? (BRF) what about undoing step??? */
    v->dx = v->odx ;  /* for mrisTrackTotalDistances */
    v->dy = v->ody ;
    v->dz = v->odz ;
#endif

#if 0
    /* write the new face positions into the filled volume */
    for (fno = 0 ; fno < v->num ; fno++)
      mrisFillFace(mris, mri_filled, v->f[fno]) ;
#else
    if (mht)
      MHTaddAllFaces(mht, mris, v) ;
#endif

  }

  direction *= -1 ;
  return(delta_t) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRISmomentumTimeStep(MRI_SURFACE *mris, float momentum, float dt, float tol, 
                     float n_averages)
{
  double  delta_t, mag ;
  int     vno ;
  VERTEX  *v ;
#if 0
  double  max_delta ;
  float   dx, dy, dz ;
#endif

  delta_t = dt * sqrt((double)n_averages+1.0) ; ;


#if 0
  /* find the largest delta, and scale the gradient by it */
  max_delta = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dx = v->dx ; dy = v->dy ; dz = v->dz ;
    mag = sqrt(dx*dx+dy*dy+dz*dz) ;
    if (mag > max_delta)
      max_delta = mag ;
  }
  if (FZERO(max_delta))
    max_delta = tol ;

  if (delta_t > MAX_MOMENTUM_MM / max_delta)   /* no bigger than 1mm */
    delta_t = MAX_MOMENTUM_MM / max_delta ;
#endif

  /* take a step in the gradient direction modulated by momentum */
  if (mris->status == MRIS_RIGID_BODY)
  {
    mris->da = delta_t * mris->alpha + momentum * mris->da ;
    mris->db = delta_t * mris->beta + momentum * mris->db ;
    mris->dg = delta_t * mris->gamma + momentum * mris->dg ;
    MRISrotate(mris, mris, mris->da, mris->db, mris->dg) ;
  }
  else for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->odx = delta_t * v->dx + momentum*v->odx ;
    v->ody = delta_t * v->dy + momentum*v->ody ;
    v->odz = delta_t * v->dz + momentum*v->odz ;
    mag = 
      sqrt(v->odx*v->odx + 
           v->ody*v->ody +
           v->odz*v->odz) ;
    if (mag > MAX_MOMENTUM_MM) /* don't let step get too big */
    {
      mag = MAX_MOMENTUM_MM / mag ;
      v->odx *= mag ; v->ody *= mag ; v->odz *= mag ;
    }
    if (vno == Gdiag_no)
    {
      float dist, dot, dx, dy, dz ;

      dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
      dist = sqrt(dx*dx+dy*dy+dz*dz) ;
      dot = dx*v->nx + dy*v->ny + dz*v->nz ;
      fprintf(stdout, "moving v %d by (%2.2f, %2.2f, %2.2f) dot=%2.2f-->"
              "(%2.1f, %2.1f, %2.1f)\n", vno, v->odx, v->ody, v->odz, 
              v->odx*v->nx+v->ody*v->ny+v->odz*v->nz,
              v->x, v->y, v->z) ;
      fprintf(stdout, "n = (%2.1f,%2.1f,%2.1f), total dist=%2.3f, total dot = %2.3f\n", 
              v->nx, v->ny, v->nz, dist, dot) ;
    }
    v->x += v->odx ; 
    v->y += v->ody ;
    v->z += v->odz ;
  }

  mrisProjectSurface(mris) ;
  return(delta_t) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Fit a quadratic form to the error function in the gradient
          direction and use it to predict the location of the minimum.
          Pick the dt which minimizes the error function among the
          sampled points, including the predicted one.
------------------------------------------------------*/
#define MAX_ENTRIES  100
static double
mrisLineMinimize(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  char    fname[STRLEN] ;
  FILE    *fp = NULL ;
  double  starting_sse, sse, min_sse, max_dt, min_delta,
          max_delta, mag, grad, delta_t, min_dt, mean_delta ;
  float   dx, dy, dz ;
  int     vno, n ;
  VERTEX  *vertex ;
  VECTOR  *vY ;
  MATRIX  *mX, *m_xTx, *m_xTx_inv, *m_xTy, *mP, *m_xT ;
  int     i, N, mini ;
  double  a, b, c, sse0, sse2, dt0, dt2, dt_in[MAX_ENTRIES],
          sse_out[MAX_ENTRIES] ;

  
  if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "%s%4.4d.dat", FileName(parms->base_name), parms->t+1);
    fp = fopen(fname, "w") ;
  }

  min_sse = starting_sse = MRIScomputeSSE(mris, parms) ;

  /* compute the magnitude of the gradient, and the max delta */
  max_delta = grad = mean_delta = 0.0f ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    dx = vertex->dx ; dy = vertex->dy ; dz = vertex->dz ;
    mag = sqrt(dx*dx+dy*dy+dz*dz) ;
    grad += dx*dx+dy*dy+dz*dz ;
    mean_delta += mag ;
#if 0
    if (!FZERO(mag))
#endif
      n++ ;
    if (mag > max_delta)
      max_delta = mag ;
  }
  mean_delta /= (float)n ;
  grad = sqrt(grad) ;

  if (FZERO(max_delta))
    return(0.0) ;       /* at a local minimum */

  /* limit the size of the largest time step */
  switch (parms->projection)
  {
  case PROJECT_SPHERE:
  case PROJECT_ELLIPSOID:
    max_dt = MAX_MM / mean_delta ;
    break ;
  case NO_PROJECTION:
    max_dt = 100.0*MAX_MM / max_delta ;
    break ;
  default:
  case PROJECT_PLANE:
    max_dt = MAX_PLANE_MM / max_delta ;
    break ;
  }
  min_dt = MIN_MM / mean_delta ;

  /* write out some data on supposed quadratic form */
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    double  delta ;
    float   predicted_sse ;
    FILE    *fp2 ;

    sprintf(fname, "nn%s%4.4d.dat",FileName(parms->base_name), parms->t+1) ;
    fp2 = fopen(fname, "w") ;

    delta = max_dt / 100.0 ;
    for (delta_t = delta ; delta_t <= max_dt ; delta_t += delta)
    {
      predicted_sse = starting_sse - grad * delta_t ;
      MRISapplyGradient(mris, delta_t) ;
      mrisProjectSurface(mris) ;
      MRIScomputeMetricProperties(mris) ;
        
      sse = MRIScomputeSSE(mris, parms) ;
      fprintf(fp2, "%f  %f  %f\n", delta_t, sse, predicted_sse) ;
      mrisProjectSurface(mris) ;
      sse = MRIScomputeSSE(mris, parms) ;
      fprintf(fp, "%f  %f  %f\n", delta_t, sse, predicted_sse) ;
      fflush(fp) ;
      MRISrestoreOldPositions(mris) ;
    }
  }

  /* pick starting step size */
  min_delta = 0.0f ; /* to get rid of compiler warning */
  for (delta_t = min_dt ; delta_t < max_dt ; delta_t *= 10.0)
  {
    MRISapplyGradient(mris, delta_t) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    sse = MRIScomputeSSE(mris, parms) ;
    if (sse <= min_sse)   /* new minimum found */
    {
      min_sse = sse ;
      min_delta = delta_t ;
    }

    /* undo step */
    MRISrestoreOldPositions(mris) ;
  }

  if (FZERO(min_delta))  /* dt=0 is min starting point, look mag smaller */
  {
    min_delta = min_dt/10.0 ;  /* start at smallest step */
    MRISapplyGradient(mris, min_delta) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    min_sse = MRIScomputeSSE(mris, parms) ;
    MRISrestoreOldPositions(mris) ;
  }

  delta_t = min_delta ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout,"grad=%2.3f, max_del=%2.3f, mean=%2.3f, max_dt=%2.1f, "
            "starting dt=%2.3f, min_dt=%2.3f\n", (float)grad, (float)max_delta, 
            mean_delta,(float)max_dt, (float)delta_t,min_dt) ;

  /* fit a quadratic form to it, and predict location of minimum */
  /* bracket the minimum by sampling on either side */
  N = 3 ;
  dt0 = min_delta - (min_delta/2) ;    dt2 = min_delta + (min_delta/2) ;
  MRISapplyGradient(mris, dt0) ;
  mrisProjectSurface(mris) ;           MRIScomputeMetricProperties(mris) ;
  sse0 = MRIScomputeSSE(mris, parms) ; MRISrestoreOldPositions(mris) ;
  
  MRISapplyGradient(mris, dt2) ;
  mrisProjectSurface(mris) ;           MRIScomputeMetricProperties(mris) ;
  sse2 = MRIScomputeSSE(mris, parms) ; MRISrestoreOldPositions(mris) ;
  
  /* now fit a quadratic form to these values */
  sse_out[0] = sse0 ; sse_out[1] = min_sse ; sse_out[2] = sse2 ;
  dt_in[0] = dt0 ; dt_in[1] = min_delta ; dt_in[2] = dt2 ;
  
  mX = MatrixAlloc(N, 3, MATRIX_REAL) ;
  vY = VectorAlloc(N, MATRIX_REAL) ;   
  
  for (i = 1 ; i <= N ; i++)
  {
    *MATRIX_RELT(mX, i, 1) = dt_in[i-1] * dt_in[i-1] ;
    *MATRIX_RELT(mX, i, 2) = 2*dt_in[i-1] ;
    *MATRIX_RELT(mX, i, 3) = 1.0f ;
    
    VECTOR_ELT(vY, i) = sse_out[i-1] ;
  }
  
  m_xT = MatrixTranspose(mX, NULL) ;
  m_xTx = MatrixMultiply(m_xT, mX, NULL) ;
  m_xTx_inv = MatrixInverse(m_xTx, NULL) ;
  if (m_xTx_inv)
  {
    m_xTy = MatrixMultiply(m_xT, vY, NULL) ;
    mP = MatrixMultiply(m_xTx_inv, m_xTy, NULL) ;
    a = RVECTOR_ELT(mP, 1) ; b = RVECTOR_ELT(mP, 2) ; c = RVECTOR_ELT(mP, 3);
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout,
              "(a,b,c) = (%2.3f, %2.3f, %2.3f), predicted min at %2.3f\n",
              a, b, c, -b/a) ;
    if (!finite(a))
      DiagBreak() ;
    MatrixFree(&mX) ; MatrixFree(&mP) ; VectorFree(&vY) ;
    MatrixFree(&m_xT) ; MatrixFree(&m_xTx) ; MatrixFree(&m_xTx_inv) ;
    MatrixFree(&m_xTy) ;
    
    dt_in[N] = 0 ;
    sse_out[N++] = starting_sse ;
    if (finite(a) && !FZERO(a))
    {
      float new_min_delta ;
      
      new_min_delta = -b/a ;
      if (new_min_delta < 10.0f*min_delta && new_min_delta > min_delta/10.0f)
      {
        MRISapplyGradient(mris, new_min_delta) ;
        mrisProjectSurface(mris) ;           
        MRIScomputeMetricProperties(mris) ;
        sse = MRIScomputeSSE(mris, parms) ;
        MRISrestoreOldPositions(mris) ;         
        dt_in[N] = new_min_delta ;
        sse_out[N++] = sse ;
      }
    }
  }
  else   /* couldn't invert matrix */
  {
    fprintf(stderr, "singular matrix in quadratic form\n") ;
    MatrixFree(&mX) ; VectorFree(&vY) ;
    MatrixFree(&m_xT) ; MatrixFree(&m_xTx) ; 
  }
  
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "sses: %2.2f  ", sse_out[0]) ;
  mini = 0 ; min_sse = sse_out[mini] ; 
  for (i = 1 ; i < N ; i++)
  {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "%2.2f  ", sse_out[i]) ;
    if (sse_out[i] < min_sse)
    {
      min_sse = sse_out[i] ;
      mini = i ;
    }
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "min %d (%2.3f)\n", mini, dt_in[mini]) ;
  MRISapplyGradient(mris, dt_in[mini]) ;
  return(dt_in[mini]) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Use a binary search in the gradient direction to find the
          location of the minimum.
------------------------------------------------------*/
static double
mrisLineMinimizeSearch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  char    fname[STRLEN] ;
  FILE    *fp = NULL ;
  double  starting_sse, sse, min_sse, max_dt, total_delta, min_delta,
          max_delta, mag, grad, delta_t, min_dt, mean_delta ;
  float   dx, dy, dz ;
  int     vno, done = 0, increasing, n ;
  VERTEX  *vertex ;
  
  if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "%s%4.4d.dat", FileName(parms->base_name), parms->t+1);
    fp = fopen(fname, "w") ;
  }

  min_sse = starting_sse = MRIScomputeSSE(mris, parms) ;

  /* compute the magnitude of the gradient, and the max delta */
  max_delta = grad = mean_delta = 0.0f ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    dx = vertex->dx ; dy = vertex->dy ; dz = vertex->dz ;
    mag = sqrt(dx*dx+dy*dy+dz*dz) ;
    grad += dx*dx+dy*dy+dz*dz ;
    mean_delta += mag ;
    if (!FZERO(mag))
      n++ ;
    if (mag > max_delta)
      max_delta = mag ;
  }
  mean_delta /= (float)n ;
  grad = sqrt(grad) ;

  if (FZERO(max_delta))
    return(0.0) ;       /* at a local minimum */

  /* limit the size of the largest time step */
  switch (parms->projection)
  {
  case PROJECT_SPHERE:
  case PROJECT_ELLIPSOID:
    max_dt = MAX_MM / mean_delta ;
    break ;
  case NO_PROJECTION:
    max_dt = 100.0*MAX_MM / max_delta ;
    break ;
  default:
  case PROJECT_PLANE:
    max_dt = MAX_PLANE_MM / max_delta ;
    break ;
  }
  min_dt = MIN_MM / mean_delta ;

  /* pick starting step size */
  min_delta = 0.0f ; /* to get rid of compiler warning */
  for (delta_t = min_dt ; delta_t < max_dt ; delta_t *= 10.0)
  {
    MRISapplyGradient(mris, delta_t) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    sse = MRIScomputeSSE(mris, parms) ;
    if (sse <= min_sse)   /* new minimum found */
    {
      min_sse = sse ;
      min_delta = delta_t ;
    }

    /* undo step */
    MRISrestoreOldPositions(mris) ;
  }

  if (FZERO(min_delta))  /* dt=0 is min starting point, look mag smaller */
  {
    min_delta = min_dt/10.0 ;  /* start at smallest step */
    MRISapplyGradient(mris, min_delta) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    min_sse = MRIScomputeSSE(mris, parms) ;
    MRISrestoreOldPositions(mris) ;
  }

  delta_t = min_delta ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout,"grad=%2.3f, max_del=%2.3f, mean=%2.3f, max_dt=%2.1f, "
            "starting dt=%2.3f, min_dt=%2.3f\n", (float)grad, (float)max_delta, 
            mean_delta,(float)max_dt, (float)delta_t,min_dt) ;

  /* now search for minimum in gradient direction */
  increasing = 1 ;
  total_delta = 0.0 ;
  min_sse = starting_sse ;
  while (!done)
  {
    MRISapplyGradient(mris, delta_t) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    sse = MRIScomputeSSE(mris, parms) ;
#if 0
    if (Gdiag & DIAG_WRITE)
      fprintf(fp, "%2.8f   %2.8f\n", total_delta+delta_t, sse) ;
#endif
    if (sse <= min_sse)   /* new minimum found */
    {
      if ((parms->projection == PROJECT_ELLIPSOID) &&  
          (total_delta+delta_t > max_dt))
        increasing = 0 ;  /* limit size of largest time step */
      min_sse = sse ;
      total_delta += delta_t ;           /* keep track of total time step */
    }
    else                 /* error increased - undo it and decrease time step */
    {
      if (increasing)    /* went too far - reverse delta_t change */
        increasing = 0 ;

      MRISrestoreOldPositions(mris) ;
      mrisProjectSurface(mris) ;
      MRIScomputeMetricProperties(mris) ;
    }
    if (total_delta + delta_t >= 10.0*min_delta)
      increasing = 0 ;
    if (increasing)
      delta_t *= 2.0 ;    /* increase time step and search further out */
    else               /* decreasing - reduce time step */
      delta_t *= 0.5 ;
    done = delta_t < min_dt ;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON && fp)
    fclose(fp) ;

  return(total_delta) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISscaleBrain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, float scale)
{
  VERTEX  *v;
  int     k;
  float   xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0 ;

  /*  mris_dst = MRIScenter(mris_src, mris_dst) ;*/
  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;

  /* find the center */
  for (k=0;k<mris_dst->nvertices;k++) 
  {
    v = &mris_dst->vertices[k];
    if (v->ripflag)
      continue ;
    if (v->x > xhi) xhi = v->x;
    if (v->x < xlo) xlo = v->x;
    if (v->y > yhi) yhi = v->y;
    if (v->y < ylo) ylo = v->y;
    if (v->z > zhi) zhi = v->z;
    if (v->z < zlo) zlo = v->z;
  }

  /* scale around the center */
  x0 = (xlo+xhi)/2.0f ; y0 = (ylo+yhi)/2.0f ; z0 = (zlo+zhi)/2.0f ;
  for (k=0;k<mris_dst->nvertices;k++) 
  {
    v = &mris_dst->vertices[k];
    if (v->ripflag)
      continue ;
    v->x = (v->x - x0) * scale + x0 ; 
    v->y = (v->y - y0) * scale + y0 ; 
    v->z = (v->z - z0) * scale + z0 ;
    if (v->x > xhi) xhi = v->x;
    if (v->x < xlo) xlo = v->x;
    if (v->y > yhi) yhi = v->y;
    if (v->y < ylo) ylo = v->y;
    if (v->z > zhi) zhi = v->z;
    if (v->z < zlo) zlo = v->z;
  }

  /* recompute the dimensions */
  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (k=0;k<mris_dst->nvertices;k++) 
  {
    v = &mris_dst->vertices[k];
    if (v->ripflag)
      continue ;
    if (v->x > xhi) xhi = v->x;
    if (v->x < xlo) xlo = v->x;
    if (v->y > yhi) yhi = v->y;
    if (v->y < ylo) ylo = v->y;
    if (v->z > zhi) zhi = v->z;
    if (v->z < zlo) zlo = v->z;
  }
  mris_dst->xlo = xlo ; mris_dst->ylo = ylo ; mris_dst->zlo = zlo ;
  mris_dst->xctr = (xhi + xlo)/2 ;
  mris_dst->yctr = (yhi + ylo)/2 ;
  mris_dst->zctr = (zhi + zlo)/2 ;
  /*  MRISupdateSurface(mris_dst) ;*/
  return(mris_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteCurvature(MRI_SURFACE *mris, char *sname)
{
  int    k ;
  float  curv;
  char   fname[STRLEN], *cp, path[STRLEN], name[STRLEN] ;
  FILE   *fp;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
    cp = strchr(sname, '.') ;
    if (!cp)
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
    else
      sprintf(fname, "%s/%s", path, sname) ;
  }
  else
  {
    FileNamePath(sname, path) ;
    FileNameOnly(sname, name) ;
    cp = strchr(sname, '.') ;
    if (!cp)
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", name) ;
    else
      sprintf(fname, "%s/%s", path, name) ;
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "writing curvature file %s\n", fname) ;

  fp = fopen(fname,"wb") ;
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISwriteCurvature: could not open %s", 
                 fname)) ;

  fwrite3(-1, fp) ;   /* same old trick - mark it as new format */
  fwriteInt(mris->nvertices,fp);
  fwriteInt(mris->nfaces,fp);
  fwriteInt(1, fp) ;    /* 1 value per vertex */
  for (k=0;k<mris->nvertices;k++)
  {
    curv = mris->vertices[k].curv ;
    fwriteFloat(curv, fp) ;
  }
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteDists(MRI_SURFACE *mris, char *sname)
{
  int    k,i ;
  float  dist ;
  char   fname[STRLEN], *cp, path[STRLEN] ;
  FILE   *fp;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
    cp = strchr(sname, '.') ;
    if (!cp)
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
    else
      sprintf(fname, "%s/%s", path, sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */
  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISwriteDists: could not open %s", 
                 fname)) ;

  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    dist = mris->vertices[k].d ;
    i = nint(dist * 100.0) ;
    fwrite2((int)i,fp);
  }
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRISwriteIntoVolume(MRI_SURFACE *mris, MRI *mri)
{
  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISreadFromVolume(MRI *mri, MRI_SURFACE *mris)
{
  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISrotate(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, 
           float alpha, float beta, float gamma)
{
  int      vno ;
  VERTEX   *vertex ;
  float    x, y, z, ca, cb, cg, sa, sb, sg, xp, yp, zp ;
  float    cacb, cacgsb, sasg, cgsa ;
  float    casbsg, cbsa, cgsasb, casg ;
  float    cacg, sasbsg, cbcg, cbsg ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  sa = sin(alpha) ; sb = sin(beta) ; sg = sin(gamma) ;
  ca = cos(alpha) ; cb = cos(beta) ; cg = cos(gamma) ;
  cacb = ca*cb ; cacgsb = ca*cg*sb ; sasg = sa*sg ; cgsa = cg*sa ;
  casbsg = ca*sb*sg ; cbsa = cb*sa ; cgsasb = cg*sa*sb ; casg = ca*sg ;
  cacg = ca*cg ; sasbsg = sa*sb*sg ; cbcg = cb*cg ; cbsg = cb*sg ;
  
  
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;

    vertex = &mris_src->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    xp = x*cacb + z*(-cacgsb - sasg) + y*(cgsa-casbsg) ;
    yp = -x*cbsa + z*(cgsasb-casg) + y*(cacg+sasbsg) ;
    zp = z*cbcg + x*sb + y*cbsg ;
    vertex->x = xp ; vertex->y = yp ; vertex->z = zp ;
  }
  return(mris_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteAreaError(MRI_SURFACE *mris, char *name)
{
  int    vno, fno, i ;
  float  area, orig_area ;
  FACE   *face ;
  VERTEX *vertex ;
  FILE   *fp;
  char   fname[STRLEN] ;
  
  MRISbuildFileName(mris, name, fname) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "writing area error file %s...", fname) ;

  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISwriteAreaError: could not open %s", 
                 fname)) ;

  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);
  for (vno = 0 ; vno < mris->nvertices; vno++)
  {
    vertex = &mris->vertices[vno] ;
    area = orig_area = 0.0f ;

    /* use average area of all faces this vertex is part of -
       this is not really correct, but should be good enough for
       visualization purposes.
       */
    for (fno = 0 ; fno < vertex->num ; fno++)
    {
      face = &mris->faces[vertex->f[fno]] ;
      area += face->area ;
      orig_area += face->orig_area ;
    }
    i = nint((area-orig_area) * 100.0f / (float)(vertex->num)) ;
    fwrite2((int)i,fp);
  }
  fclose(fp);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "done.\n") ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteAreaErrorToValFile(MRI_SURFACE *mris, char *name)
{
  int    vno, fno ;
  float  area, orig_area ;
  FACE   *face ;
  VERTEX *v ;
  
  for (vno = 0 ; vno < mris->nvertices; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    area = orig_area = 0.0f ;

    /* use average area of all faces this vertex is part of -
       this is not really correct, but should be good enough for
       visualization purposes.
       */
    for (fno = 0 ; fno < v->num ; fno++)
    {
      face = &mris->faces[v->f[fno]] ;
      area += face->area ;
      orig_area += face->orig_area ;
    }
    v->val = (area-orig_area) / (float)(v->num) ;
  }

  MRISwriteValues(mris, name) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteAngleError(MRI_SURFACE *mris, char *fname)
{
  int    vno, fno, ano, i ;
  float  error ;
  FILE   *fp;
  FACE   *face ;
  VERTEX *v ;
  
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "writing angular error file %s...", fname) ;

  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISwriteAngleError: could not open %s", 
                 fname)) ;

  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    error = 0.0f ;
    for (fno = 0 ; fno < v->num ; fno++)
    {
      face = &mris->faces[v->f[fno]] ;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      {
        error += fabs(deltaAngle(face->angle[ano],face->orig_angle[ano]));
      }
      error /= (float)(v->num*ANGLES_PER_TRIANGLE) ;
    }
    i = DEGREES(error) * 100.0 ;
    fwrite2((int)i,fp);
  }
  fclose(fp);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "done.\n") ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if AVERAGE_AREAS
static int
mrisAverageAreas(MRI_SURFACE *mris, int num_avgs, int which)
{
  int    i, vno, vnb, *pnb, fno, num, vf, nfno ;
  float  area ;
  VERTEX *v, *vn ;
  FACE   *f ;

  for (i = 0 ; i < num_avgs ; i++)
  {
    switch (which)
    {
    case ORIG_AREAS:
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;

        for (fno = 0 ; fno < v->num ; fno++) /* each face of this vertex */
        {
          f = &mris->faces[v->f[fno]] ;  /* pointer to the face */
          if (f->ripflag)
            continue ;
          area = 0.0f ;

          /* now go through each vertex associated with this face */
          for (vf = 0 ; vf < VERTICES_PER_FACE ; vf++)
          {
            vn = &mris->vertices[f->v[vf]] ;
            num += vn->num ;
            for (nfno = 0 ; nfno < vn->num ; nfno++)
              area += vn->orig_tri_area[nfno] ;
          }
          area /= (float)num ;
          v->orig_tri_area[fno] = area ;
        }
      }
      break ;
    case CURRENT_AREAS:
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;

        for (fno = 0 ; fno < v->num ; fno++) /* each face of this vertex */
        {
          f = &mris->faces[v->f[fno]] ;  /* pointer to the face */
          if (f->ripflag)
            continue ;
          area = 0.0f ;

          /* now go through each vertex associated with this face */
          for (vf = 0 ; vf < VERTICES_PER_FACE ; vf++)
          {
            vn = &mris->vertices[f->v[vf]] ;
            num += vn->num ;
            for (nfno = 0 ; nfno < vn->num ; nfno++)
              area += vn->tri_area[nfno] ;
          }
          area /= (float)num ;
          v->tri_area[fno] = area ;
        }
      }
      break ;
    }
  }
  return(NO_ERROR) ;
}

#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#include "stats.h"
int
MRISsampleStatVolume(MRI_SURFACE *mris, void *vsv,int time_point,
                     int coords)
{
  VERTEX   *v ;
  int      vno, xv, yv, zv, width, height, depth ;
  Real     x, y, z, xt, yt, zt ;
  STAT_VOLUME *sv = (STAT_VOLUME *)vsv ;

  if (time_point >= sv->mri_pvals[0]->nframes)
    ErrorExit(ERROR_BADPARM, 
              "MRISsampleStatVolume: time point (%d) out of bounds [%d, %d]\n",
              time_point, 0, sv->mri_pvals[0]->nframes-1) ;
  width  = sv->mri_pvals[0]->width ;
  height  = sv->mri_pvals[0]->height ;
  depth  = sv->mri_pvals[0]->depth ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == 47)
      DiagBreak() ;
    v = &mris->vertices[vno] ;
    x = (Real)v->x ; y = (Real)v->y ; z = (Real)v->z ;

    /* now convert them into talairach space */
    switch (coords)
    {
    case TALAIRACH_COORDS:
      MRIworldToTalairachVoxel(sv->mri_pvals[0], x, y, z, &xt, &yt, &zt) ;
      break ;
    case SPHERICAL_COORDS:
      x = (Real)v->cx ; y = (Real)v->cy ; z = (Real)v->cz ;
      MRIworldToVoxel(sv->mri_pvals[0], x, y, z, &xt, &yt, &zt) ;
      break ;
    default:
      MRIworldToVoxel(sv->mri_pvals[0], x, y, z, &xt, &yt, &zt) ;
      break ;
    }
    xv = nint(xt) ; yv = nint(yt) ; zv = nint(zt) ;
    if (xv >= 0 && xv < width && yv >= 0 && yv < height && zv>=0&&zv<depth)
      v->val = MRIFseq_vox(sv->mri_pvals[0], xv, yv, zv, time_point) ;
    if (vno == 1446)
      DiagBreak() ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteCurvatureToWFile(MRI_SURFACE *mris, char *fname)
{
  int k,num;                   /* loop counters */
  float f;
  FILE *fp;
  double sum=0,sum2=0,max= -1000,min=1000;

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "writing out surface values to %s.\n", fname) ;

  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorExit(ERROR_NOFILE, "Can't create file %s\n",fname) ;

  num = mris->nvertices ;
  fwrite2(0,fp);
  fwrite3(num,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    fwrite3(k,fp);
    f = mris->vertices[k].curv;
    if (!finite(f))
      ErrorPrintf(ERROR_BADPARM, 
                  "MRISwriteCurvatureToWFile(%s): val at vertex %d is not"
                  "finite", fname, k) ;
    
    fwriteFloat(f, fp) ;
    sum += f;
    sum2 += f*f;
    if (f>max) max=f;
    if (f<min) min=f;
  }
  fclose(fp);
  sum /= num;
  sum2 = sqrt(sum2/num-sum*sum);
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "avg = %2.3f, stdev = %2.3f, min = %2.3f, max = %2.3f\n",
            sum,sum2,min,max);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteValues(MRI_SURFACE *mris, char *sname)
{
  int k,num;                   /* loop counters */
  float f;
  char  fname[STRLEN], *cp ;
  FILE *fp;
  double sum=0,sum2=0,max= -1000,min=1000;

#if 1
  MRISbuildFileName(mris, sname, fname) ;
#else
  strcpy(fname, sname) ;
#endif
  cp = strrchr(fname, '.') ;
  if (!cp || *(cp+1) != 'w')
    strcat(fname, ".w") ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "writing surface values to %s.\n", fname) ;

  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorExit(ERROR_NOFILE, "Can't create file %s\n",fname) ;

  for (k=0,num=0;k<mris->nvertices;k++) 
    if (mris->vertices[k].val!=0) num++;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("num = %d\n",num);
  fwrite2(0,fp);
  fwrite3(num,fp);

  for(k=0;k<mris->nvertices;k++){
    if(mris->vertices[k].val != 0){
      fwrite3(k,fp);
      f = mris->vertices[k].val;
      if (!finite(f))
        ErrorPrintf(ERROR_BADPARM, 
                    "MRISwriteValues(%s): val at vertex %d is not finite",
                    fname, k) ;

      fwriteFloat(f, fp) ;
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
  }
  fclose(fp);

  if(num > 0){
    sum /= num;
    sum2 = (sum2/num-sum*sum);
    if(sum2 > 0) sum2 = sqrt(sum2) ;
    else          sum2 = 0 ;
    printf("avg = %2.3f, stdev = %2.3f, min = %2.3f, max = %2.3f\n",
            sum,sum2,min,max);
  }
  else
    printf("Warning: all vertex values are zero\n");
    

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadAnnotation(MRI_SURFACE *mris, char *sname)
{
  int   i,j,vno,num, need_hemi ;
  FILE  *fp;
  char  *cp, fname[STRLEN], path[STRLEN], fname_no_path[STRLEN] ;
#if 0
  int   numannothist;
  float f;
  char  histfname[STRLEN], freqfname[STRLEN];
#endif

  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNameOnly(sname, fname_no_path) ;
    cp = strstr(fname_no_path, ".annot") ;
    if (!cp)
      strcat(fname_no_path, ".annot") ;

    need_hemi = stricmp(fname_no_path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh") ;

    FileNamePath(mris->fname, path) ;
    if (!need_hemi)
      sprintf(fname, "%s/../label/%s", path, fname_no_path) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/../label/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh",fname_no_path);
  }
  else
  {
    strcpy(fname, sname) ;  /* full path specified */
    cp = strstr(fname, ".annot") ;
    if (!cp)
      strcat(fname, ".annot") ;
  }

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "could not read annot file %s",
                                fname)) ;
  MRISclearAnnotations(mris) ;
  num = freadInt(fp) ;
  for (j=0;j<num;j++)
  {
    vno = freadInt(fp) ; i = freadInt(fp) ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (vno>=mris->nvertices||vno<0)
      fprintf(stderr, "MRISreadAnnotation: vertex index out of range: %d i=%d\n",vno,i);
    else
      mris->vertices[vno].annotation = i;
  }

	while (!feof(fp))
	{
		int tag ;

		tag = freadInt(fp) ;
		switch (tag)
		{
		case TAG_COLORTABLE:
			fprintf(stderr, "reading colortable from annotation file...\n") ;
			mris->ct = CTABreadFrom(fp) ;
			fprintf(stderr, "colortable with %d entries read (originally %s)\n", mris->ct->nbins, mris->ct->fname) ;
			break ;
		default:
			break ;
		}
	}

  fclose(fp);

#if 0
  for (vno=0;vno<vertex_index;vno++)
    vertex[vno].annotfreq=1;

  sprintf(freqfname,"%s.freq",fname);
  fp = fopen(freqfname,"r");
  if (fp!=NULL)
  {
    printf("file %s read\n",freqfname);
    for (vno=0;vno<vertex_index;vno++)
      vertex[vno].annotfreq=0;
    fread(&num,1,sizeof(int),fp);
    printf("surfer: num=%d\n",num);
    for (j=0;j<num;j++)
    {
      fread(&vno,1,sizeof(int),fp);
      fread(&f,1,sizeof(float),fp);
      if (vno>=vertex_index||vno<0)
        printf("surfer: vertex index out of range: %d f=%f\n",vno,f);
      else
        vertex[vno].annotfreq = f;
    }
    fclose(fp);
  }

  sprintf(histfname,"%s.hist",fname);
  fp = fopen(histfname,"r");
  if (fp!=NULL)
  {
    printf("file %s read\n",histfname);
    for (vno=0;vno<vertex_index;vno++)
      vertex[vno].numannothist=0;
    fread(&num,1,sizeof(int),fp);
    printf("surfer: num=%d\n",num);
    for (j=0;j<num;j++)
    {
      fread(&vno,1,sizeof(int),fp);
      fread(&numannothist,1,sizeof(int),fp);
      if (vno>=vertex_index||vno<0)
        printf("surfer: vertex index out of range: %d f=%f\n",vno,f);
      else
      {
        vertex[vno].numannothist = numannothist;
        vertex[vno].annothistlabel = calloc(numannothist,sizeof(int));
        vertex[vno].annothistcount = calloc(numannothist,sizeof(int));
        for (i=0;i<numannothist;i++)
          fread(&vertex[vno].annothistlabel[i],1,sizeof(int),fp);
        for (i=0;i<numannothist;i++)
          fread(&vertex[vno].annothistcount[i],1,sizeof(int),fp);
      }
    }
    fclose(fp);
  }
#endif

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteAnnotation(MRI_SURFACE *mris, char *sname)
{
  int   i,vno, need_hemi ;
  FILE  *fp;
  char  *cp, fname[STRLEN], path[STRLEN], fname_no_path[STRLEN];

  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNameOnly(sname, fname_no_path) ;
    cp = strstr(fname_no_path, ".annot") ;
    if (!cp)
      strcat(fname_no_path, ".annot") ;

    need_hemi = 
      !stricmp(fname, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh") ;

    FileNamePath(mris->fname, path) ;
    if (!need_hemi)
      sprintf(fname, "%s/../label/%s", path, fname_no_path) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/../label/%s_%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh",fname_no_path);
  }
  else
  {
    strcpy(fname, sname) ;  /* full path specified */
    cp = strstr(fname, ".annot") ;
    if (!cp)
      strcat(fname, ".annot") ;
  }

  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "could not write annot file %s",
                                fname)) ;
  fwriteInt(mris->nvertices, fp) ;
  for (vno=0;vno<mris->nvertices;vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    i = mris->vertices[vno].annotation ;
    fwriteInt(vno,fp) ; i = fwriteInt(i,fp) ;
  }

	if (mris->ct)   /* also write annotation in */
	{
		printf("writing colortable into annotation file...\n") ;
		fwriteInt(TAG_COLORTABLE, fp) ;
		CTABwriteInto(fp, mris->ct);
	}

  fclose(fp);

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadValues(MRI_SURFACE *mris, char *sname)
{
  int   i,k,num,ilat, vno ;
  float f;
  float lat, *cvec;
  FILE  *fp;
  char  *cp, fname[STRLEN] ;

  cvec = MRISreadCurvatureVector(mris, sname) ;
  if (cvec)
  {
    printf("reading values from curvature-format file...\n") ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      if (vno == Gdiag_no)
        DiagBreak() ;
      mris->vertices[vno].val = cvec[vno];
    }
    free(cvec) ;
  }
  else
  {
    strcpy(fname, sname) ;
    cp = strrchr(fname, '.') ;
    if (!cp || *(cp+1) != 'w')
      strcat(fname, ".w") ;
    fp = fopen(fname,"rb");
    if (fp==NULL) 
      ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE,
                                 "MRISreadValues: File %s not found\n",fname));
    fread2(&ilat,fp);
    lat = ilat/10.0;
    
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].val=0;
    if (fread3(&num,fp) < 1)
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE, 
                   "MRISreadValues(%s): could not read # of vertices",
                   fname)) ;
    for (i=0;i<num;i++)
    {
      if (fread3(&k,fp) < 1)
        ErrorReturn(ERROR_BADFILE,
                    (ERROR_BADFILE, 
                     "MRISreadValues(%s): could not read %dth vno",
                     fname, i)) ;
      f = freadFloat(fp) ;
      if (k>=mris->nvertices||k<0)
        printf("MRISreadValues: vertex index out of range: %d f=%f\n",k,f);
      else
      {
        if (k == Gdiag_no)
          DiagBreak() ;
        mris->vertices[k].val = f;
      }
    }
    fclose(fp);
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadValuesBak(MRI_SURFACE *mris, char *fname)
{
  int i,k,num,ilat;
  float f;
  float lat;
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE,
                               "MRISreadValuesBak: File %s not found\n",fname));
  fread2(&ilat,fp);
  lat = ilat/10.0;

  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].valbak=0;
  fread3(&num,fp);
  for (i=0;i<num;i++)
    {
    fread3(&k,fp);
    f = freadFloat(fp) ;
    if (k>=mris->nvertices||k<0)
      printf("MRISreadValuesBak: vertex index out of range: %d f=%f\n",k,f);
    else
      {
      mris->vertices[k].valbak = f;
      }
    }
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadImagValues(MRI_SURFACE *mris, char *fname)
{
  int i,k,num,ilat;
  float f;
  float lat;
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE,"MRISreadImagValues: File %s not found\n",fname));
  fread2(&ilat,fp);
  lat = ilat/10.0;

  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].imag_val=0;
  fread3(&num,fp);
  for (i=0;i<num;i++)
  {
    fread3(&k,fp);
    f = freadFloat(fp);
    if (k>=mris->nvertices||k<0)
      printf("MRISreadImagValues: vertex index out of range: %d f=%f\n",k,f);
/*
    else if (mris->vertices[k].dist!=0)
      printf("MRISreadImagValues: subsample and data file mismatch\n");
*/
    else
    {
      mris->vertices[k].imag_val = f;
      /*      mris->vertices[k].dist=0;*/
    }
  }
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyValuesToImagValues(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->imag_val = v->val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadInflatedCoordinates(MRI_SURFACE *mris, char *sname)
{
  if (!sname)
    sname = "inflated" ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  if (MRISreadVertexPositions(mris, sname) != NO_ERROR)
    return(Gerror) ;
  MRISsaveVertexPositions(mris, INFLATED_VERTICES) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadFlattenedCoordinates(MRI_SURFACE *mris, char *sname)
{
  int    vno, fno ;
  VERTEX *v ;
  FACE   *f ;

  if (!sname)
    sname = "patch" ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  if (MRISreadPatchNoRemove(mris, sname) != NO_ERROR)
    return(Gerror) ;
  MRISsaveVertexPositions(mris, FLATTENED_VERTICES) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      v->z = -10000 ;
      v->ripflag = 0 ;
    }
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
    {
      f->ripflag = 0 ;
    }
  }
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadCanonicalCoordinates(MRI_SURFACE *mris, char *sname)
{
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  if (MRISreadVertexPositions(mris, sname) != NO_ERROR)
    return(Gerror) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeCanonicalCoordinates(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadPatchNoRemove(MRI_SURFACE *mris, char *pname)
{
  int         ix, iy, iz, k, i, j, npts ;
  double        rx, ry, rz;
  FILE        *fp ;
  char        fname[STRLEN] ;
  int         type;
  char        line[256];
  char        *cp;

#if 0
  char        path[STRLEN], *cp ;
  cp = strchr(pname, '/') ;
  if (cp)
    strcpy(fname, pname) ;    /* path already specified */
  else                        /* no path - use same as was used in MRISread */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/%s", path, pname) ;
  }
#else
  MRISbuildFileName(mris, pname, fname) ;
#endif

  // check whether the patch file is ascii or binary
  type = MRISfileNameType(fname) ; /* using extension to get type */
  if (type == MRIS_ASCII_TRIANGLE_FILE)  /* .ASC */
  {
    fp = fopen(fname, "r");
    if (!fp)
      ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,
				"MRISreadPatch(%s): could not open file", fname));
    cp = fgetl(line, 256, fp);  // this would skip # lines
    sscanf(cp, "%d %*s", &npts);   // get points
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "reading patch %s with %d vertices (%2.1f%% of total)\n",
	      pname, npts, 100.0f*(float)npts/(float)mris->nvertices) ;

    // set all vertices ripflag to be true
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].ripflag = TRUE;
    
    // go through points
    for (j=0;j<npts;j++)
    {
      // read int 
      if((cp = fgetl(line, 256, fp)))
	sscanf(cp, "%d", &i);
      else
	ErrorReturn(ERROR_BADPARM,
		    (ERROR_BADPARM, "MRISreadPatch(%s): could not read line for point %d\n", fname, j));

      // if negative, flip it
      if (i<0)
      {
	k = -i-1; // convert it to zero based number
	if (k < 0 || k >= mris->nvertices)
	  ErrorExit(ERROR_BADFILE, 
		    "MRISreadPatch: bad vertex # (%d) found in patch file", k) ;
	// negative -> set the border to be true
	mris->vertices[k].border = TRUE;
      } 
      // if positive
      else
      {
	k = i-1; // vertex number is zero based
	if (k < 0 || k >= mris->nvertices)
	  ErrorExit(ERROR_BADFILE, 
		    "MRISreadPatch: bad vertex # (%d) found in patch file", k) ;
	// positive -> set the border to be false
	mris->vertices[k].border = FALSE;
      }
      // rip flag for this vertex to be false
      mris->vertices[k].ripflag = FALSE;
      // read 3 positions
      if ((cp = fgetl(line, 256, fp)))
	sscanf(cp, "%lf %lf %lf", &rx, &ry, &rz);
      else
	ErrorReturn(ERROR_BADPARM,
		    (ERROR_BADPARM, 
		     "MRISreadPatch(%s): could not read 3 floating values line for point %d\n", 
		     fname, j));

      // convert it to mm, i.e. change the vertex position
      mris->vertices[k].x = rx;
      mris->vertices[k].y = ry;
      mris->vertices[k].z = rz;
      // change the hi, lo values
      if (mris->vertices[k].x > mris->xhi) mris->xhi = mris->vertices[k].x;
      if (mris->vertices[k].x < mris->xlo) mris->xlo = mris->vertices[k].x;
      if (mris->vertices[k].y > mris->yhi) mris->yhi = mris->vertices[k].y;
      if (mris->vertices[k].y < mris->ylo) mris->ylo = mris->vertices[k].y;
      if (mris->vertices[k].z > mris->zhi) mris->zhi = mris->vertices[k].z;
      if (mris->vertices[k].z < mris->zlo) mris->zlo = mris->vertices[k].z;
      if (k == Gdiag_no && Gdiag & DIAG_SHOW)
	fprintf(stdout, "vertex %d read @ (%2.2f, %2.2f, %2.2f)\n",k,
		mris->vertices[k].x,mris->vertices[k].y,mris->vertices[k].z) ;
    }

  }
  /////////////////////////////////////////////////////////////////////////
  // here file was binary
  /////////////////////////////////////////////////////////////////////////
  else
  {
    fp = fopen(fname, "rb") ;
    if (!fp)
      ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,
				"MRISreadPatch(%s): could not open file", fname));
    
    // read number of vertices
    npts = freadInt(fp) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "reading patch %s with %d vertices (%2.1f%% of total)\n",
	      pname, npts, 100.0f*(float)npts/(float)mris->nvertices) ;
    // set all vertices ripflag to be true
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].ripflag = TRUE;
    
    // go through points
    for (j=0;j<npts;j++)
    {
      // read int 
      i = freadInt(fp) ;
      // if negative, flip it
      if (i<0)
      {
	k = -i-1; // convert it to zero based number
	if (k < 0 || k >= mris->nvertices)
	  ErrorExit(ERROR_BADFILE, 
		    "MRISreadPatch: bad vertex # (%d) found in patch file", k) ;
	// negative -> set the border to be true
	mris->vertices[k].border = TRUE;
      } 
      // if positive
      else
      {
	k = i-1; // vertex number is zero based
	if (k < 0 || k >= mris->nvertices)
	  ErrorExit(ERROR_BADFILE, 
		    "MRISreadPatch: bad vertex # (%d) found in patch file", k) ;
	// positive -> set the border to be false
	mris->vertices[k].border = FALSE;
      }
      // rip flag for this vertex to be false
      mris->vertices[k].ripflag = FALSE;
      // read 3 positions
      fread2(&ix,fp);
      fread2(&iy,fp);
      fread2(&iz,fp);
      // convert it to mm, i.e. change the vertex position
      mris->vertices[k].x = ix/100.0;
      mris->vertices[k].y = iy/100.0;
      mris->vertices[k].z = iz/100.0;
      // change the hi, lo values
      if (mris->vertices[k].x > mris->xhi) mris->xhi = mris->vertices[k].x;
      if (mris->vertices[k].x < mris->xlo) mris->xlo = mris->vertices[k].x;
      if (mris->vertices[k].y > mris->yhi) mris->yhi = mris->vertices[k].y;
      if (mris->vertices[k].y < mris->ylo) mris->ylo = mris->vertices[k].y;
      if (mris->vertices[k].z > mris->zhi) mris->zhi = mris->vertices[k].z;
      if (mris->vertices[k].z < mris->zlo) mris->zlo = mris->vertices[k].z;
      if (k == Gdiag_no && Gdiag & DIAG_SHOW)
	fprintf(stdout, "vertex %d read @ (%2.2f, %2.2f, %2.2f)\n",k,
		mris->vertices[k].x,mris->vertices[k].y,mris->vertices[k].z) ;
    }
  }
  fclose(fp);
  // remove ripflag set vertices
  MRISripFaces(mris);
  // set the patch flag
  mris->patch = 1 ;
  mris->status = MRIS_CUT ;
  // recalculate properties
  mrisComputeBoundaryNormals(mris);
  mrisSmoothBoundaryNormals(mris,10);
  MRISupdateSurface(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadPatch(MRI_SURFACE *mris, char *pname)
{
  int ret  ;

  // update the vertices in patch file
  ret = MRISreadPatchNoRemove(mris, pname) ;
  if (ret != NO_ERROR)
    return(ret) ;
  // remove ripflag set vertices
  MRISremoveRipped(mris) ;
  // recalculate properties
  MRISupdateSurface(mris) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISripFaces(MRI_SURFACE *mris)
{
  int n,k;
  face_type *f;

  for (k=0;k<mris->nfaces;k++)
    mris->faces[k].ripflag = FALSE;
  for (k=0;k<mris->nfaces;k++)
  {
    f = &mris->faces[k];
    for (n=0;n<VERTICES_PER_FACE;n++)
      if (mris->vertices[f->v[n]].ripflag)
        f->ripflag = TRUE;
  }
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].border = FALSE;
  for (k=0;k<mris->nfaces;k++)
  if (mris->faces[k].ripflag)
  {
    f = &mris->faces[k];
    for (n=0;n<VERTICES_PER_FACE;n++)
        mris->vertices[f->v[n]].border = TRUE;
  }
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisRipVertices(MRI_SURFACE *mris)
{
  int     fno, n ;
  VERTEX  *v ;
  FACE    *f ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag == 0)
      continue ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;
      v->ripflag = 1 ;
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwritePatch(MRI_SURFACE *mris, char *fname)
{
  int    k,i,npts, type ;
  float  x,y,z;
  FILE   *fp;

  type = MRISfileNameType(fname) ;
  if (type == MRIS_ASCII_TRIANGLE_FILE) // extension is ASC
    return(MRISwritePatchAscii(mris, fname)) ;
  else if (type == MRIS_GEO_TRIANGLE_FILE) // extension is GEO
    return(MRISwriteGeo(mris, fname)) ;

  // binary file write
  // count number of points
  npts = 0;
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag) npts++;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "writing surface patch with npts=%d to %s\n",npts,fname);
  // binary write
  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, 
                 "MRISwritePatch: can't create file %s\n",fname)) ;
  // write num points 
  fwriteInt(npts, fp) ;
  // go through all points
  for (k=0;k<mris->nvertices;k++)
  if (!mris->vertices[k].ripflag)
  {
    i = (mris->vertices[k].border)? (-(k+1)): (k+1);
    fwriteInt(i, fp) ;
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
    fwrite2((int)(x*100),fp);
    fwrite2((int)(y*100),fp);
    fwrite2((int)(z*100),fp);
/*
    printf("k=%d, i=%d\n",k,i);
*/
  }
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsmoothSurfaceNormals(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  nx, ny, nz, num, len ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      nx = v->nx ; ny = v->ny ; nz = v->nz ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        nx += vn->nx ; ny += vn->ny ; nz += vn->nz ;
      }
      num++ ;   /* account for central vertex */
      v->tdx = nx / num ;
      v->tdy = ny / num ;
      v->tdz = nz / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
			len = sqrt(v->tdx*v->tdx+v->tdy*v->tdy+v->tdz*v->tdz) ;
			if (FZERO(len))
				len = 1 ;
      v->nx = v->tdx/len ; v->ny = v->tdy/len ; v->nz = v->tdz/len ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisSmoothBoundaryNormals(MRI_SURFACE *mris, int niter)
{
#if 0
  int iter,k,m,n;
  vertex_type *v;
  float sumx,sumy,r;

  for (iter=0;iter<niter;iter++)
  {
    for (k=0;k<mris->nvertices;k++)
    if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
    {
      mris->vertices[k].obnx = mris->vertices[k].bnx;
      mris->vertices[k].obny = mris->vertices[k].bny;
    }
    for (k=0;k<mris->nvertices;k++)
    if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
    {
      v = &mris->vertices[k];
      n = 1;
      sumx = v->obnx;
      sumy = v->obny;
      for (m=0;m<v->vnum;m++)
      if ((!mris->vertices[v->v[m]].ripflag)&&mris->vertices[v->v[m]].border)
      {
          sumx += mris->vertices[v->v[m]].obnx;
          sumy += mris->vertices[v->v[m]].obny;
          n++;
      }
      v->bnx = (n>0)?sumx/n:0;
      v->bny = (n>0)?sumy/n:0;
    }
  }
  for (k=0;k<mris->nvertices;k++)
  if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
  {
    r = sqrt(SQR(mris->vertices[k].bnx)+SQR(mris->vertices[k].bny));
    if (r>0)
    {
      mris->vertices[k].bnx /= r;
      mris->vertices[k].bny /= r;
    }
  }
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeBoundaryNormals(MRI_SURFACE *mris)
{
#if 0
  int      k,m,n;
  VERTEX   *v;
  float    sumx,sumy,r,nx,ny,f;

  for (k=0;k<mris->nvertices;k++)
  if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
  {
    v = &mris->vertices[k];
    n = 0;
    sumx = 0;
    sumy = 0;
    for (m=0;m<v->vnum;m++)
    if (!mris->vertices[v->v[m]].ripflag)
    {
      sumx += v->x-mris->vertices[v->v[m]].x;
      sumy += v->y-mris->vertices[v->v[m]].y;
      n++;
    }
    v->bnx = (n>0)?sumx/n:0;
    v->bny = (n>0)?sumy/n:0;
  }
  for (k=0;k<mris->nvertices;k++)
  if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
  {
    v = &mris->vertices[k];
    n = 0;
    sumx = 0;
    sumy = 0;
    for (m=0;m<v->vnum;m++)
    if ((!mris->vertices[v->v[m]].ripflag)&&mris->vertices[v->v[m]].border)
    {
      nx = -(v->y-mris->vertices[v->v[m]].y); 
      ny = v->x-mris->vertices[v->v[m]].x; 
      f = nx*v->bnx+ny*v->bny;
/*
      f = (f<0)?-1.0:(f>0)?1.0:0.0;
*/
      sumx += f*nx;
      sumy += f*ny;
      n++;
    }
    v->bnx = (n>0)?sumx/n:0;
    v->bny = (n>0)?sumy/n:0;
  }
  for (k=0;k<mris->nvertices;k++)
  if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
  {
    r = sqrt(SQR(mris->vertices[k].bnx)+SQR(mris->vertices[k].bny));
    if (r>0)
    {
      mris->vertices[k].bnx /= r;
      mris->vertices[k].bny /= r;
    }
  }
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISflattenPatchRandomly(MRI_SURFACE *mris)
{
  float   extent ;
  int     vno ;

  extent = sqrt(mris->total_area) ;
  for (vno=0;vno<mris->nvertices;vno++)
  {
    mris->vertices[vno].x = randomNumber(-extent, extent) ;
    mris->vertices[vno].y = randomNumber(-extent, extent) ;
    mris->vertices[vno].z = 0;
    if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stdout, "vertex %d flattened @ (%2.2f, %2.2f, %2.2f)\n",vno,
              mris->vertices[vno].x,mris->vertices[vno].y,
              mris->vertices[vno].z) ;
  }
  mris->status = MRIS_PLANE ;
  MRIScomputeMetricProperties(mris) ;
  mrisOrientPlane(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISflattenPatch(MRI_SURFACE *mris)
{
  float       x,y,z,d,d1,d2;
  float       nx,ny,nz ;
  VERTEX      *v;
  int         k, an;

  x = y = z = nx = ny = nz = 0; an = 0;

  /* calculate average normal and vertex position */
  MRIScomputeNormals(mris) ;
  for (k=0;k<mris->nvertices;k++)  
  if (!mris->vertices[k].ripflag)
  {
    v = &mris->vertices[k];
    x += v->x;
    y += v->y;
    z += v->z;
#if 0
    if (!FZERO(v->nx))
      fprintf(stdout, "vertex %d, normal = (%2.3f, %2.3f, %2.3f)\n",
              k, v->nx, v->ny, v->nz) ;
#endif
    nx += v->nx;
    ny += v->ny;
    nz += v->nz;
    an++;
  }
  x /= an; y /= an; z /= an;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "flatten: avg p = {%2.1f, %2.1f, %2.1f}\n",x,y,z);
    fprintf(stdout, "flatten: sum n = {%2.2f, %2.2f, %2.2f}\n",nx,ny,nz);
  }
#if 0
  /* or override with direct front,back */
  if (project==POSTERIOR) { nx = nz = 0.0; ny = -1.0; }
  if (project==ANTERIOR)  { nx = nz = 0.0; ny = 1.0;  }
#endif

  /* make the average normal unit length */
  d = sqrt(nx*nx+ny*ny+nz*nz);
  nx /= d; ny /= d; nz /= d;
  d = sqrt(nx*nx+ny*ny);
  if (!FZERO(d))  /* not already in a plane */
  {
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "flatten: norm n = {%2.2f, %2.2f, %2.2f}\n",nx,ny,nz);
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        v = &mris->vertices[k];
        v->x -= x;
        v->y -= y;
        v->z -= z;
        d1 = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
        transform(&v->x,&v->y,&v->z,nx,ny,nz,d);
        d2 = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
        if (fabs(d1-d2)>0.0001) 
          printf("flatten: d1=%f, d2=%f\n",d1,d2);
        transform(&v->nx,&v->ny,&v->nz,nx,ny,nz,d);
      }

    /* print transform matrix in tmp dir */
#if 0
    sprintf(fname,"%s/surfer.mat",dir);
    fp = fopen(fname,"w");
    if (fp==NULL)
      ErrorPrintf(ERROR_NOFILE, "flatten: can't create file %s\n",fname);
    else 
    {
      fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",  nx*nz/d,  -nx,  ny/d,  0.0);
      fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",    d,       nz,   0.0, 0.0);
      fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",  -ny*nz/d,  ny,   nx/d, 0.0);
      fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",   0.0,      0.0,  0.0, 1.0);
      fclose(fp);
      if (Gdiag & DIAG_SHOW)
        printf("flatten: file %s written\n",fname);
    }
#endif

    transform(&nx,&ny,&nz,nx,ny,nz,d);
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "flatten: transformed n = {%2.1f, %2.1f, %2.1f}\n",
              nx,ny,nz);
  }
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].z = 0;
    if (k == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stdout, "vertex %d flattened @ (%2.2f, %2.2f, %2.2f)\n",k,
              mris->vertices[k].x,mris->vertices[k].y,mris->vertices[k].z) ;
  }

  mris->status = MRIS_PLANE ;
  MRIScomputeMetricProperties(mris) ;
  if (Gdiag & DIAG_SHOW && Gdiag_no >= 0)
  {
    int    n ;
    VERTEX *v, *vn ;
    
    v = &mris->vertices[Gdiag_no] ;
    fprintf(stdout, 
            "%d @ (%2.1f, %2.1f, %2.1f): area %2.3f, oa %2.3f, nz=%2.3f, vnum=%d\n",
            Gdiag_no, v->x, v->y, v->z, v->area, v->origarea, v->nz, v->vnum) ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      fprintf(stdout, 
              "%d @ (%2.1f, %2.1f, %2.1f): area %2.3f, oa %2.3f, nz=%2.3f, "
              "vnum=%d, d=%2.2f, od=%2.2f\n", 
              v->v[n], vn->x, vn->y, vn->z, vn->area, v->origarea, 
              vn->nz, vn->vnum, v->dist[n], v->dist_orig[n]) ;
    }
    fprintf(stdout, "\n") ;
  }


  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          2 vects ortho to summed normal 
------------------------------------------------------*/
static int
transform(float *xptr, float *yptr, float *zptr, 
          float nx, float ny, float nz, float d)  
{
  float x = *xptr, y = *yptr, z = *zptr;

  *zptr = nx*x + ny*y + nz*z;
  *yptr = -ny/d*x + nx/d*y;
  *xptr = nx*nz/d*x + ny*nz/d*y - d*z;
/*
  printf("transform {%f,%f,%f} -> {%f,%f,%f}\n",
         x,y,z,*xptr,*yptr,*zptr);
*/
  return(NO_ERROR) ;
}
static int
mrisComputeLinkTerm(MRI_SURFACE *mris, double l_link, int pial)
{
  int     vno ;
  VERTEX  *v ;
  float   dx, dy, dz, lx, ly, lz, len ;

  if (FZERO(l_link))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

		lx = v->pialx-v->origx ; ly = v->pialy-v->origy ; lz = v->pialz-v->origz ; 
		len = sqrt(lx*lx + ly*ly + lz*lz) ;
		if (len < .25)  /* can't accurately estimate vector connecting white and pial */
			continue ;
		lx /= len ; ly /= len ; lz /= len ;

		dx = l_link*(v->nx - lx) ;
		dy = l_link*(v->ny - ly) ;
		dz = l_link*(v->nz - lz) ;

		if (pial == 0)
		{
			dx *= -1 ; dy *= -1 ; dz *= -1 ;
		}
		
    v->dx += dx ;
    v->dy += dy ;
    v->dz += dz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d link %s term: (%2.3f, %2.3f, %2.3f), Nl=(%2.1f, %2.1f, %2.1f), Ns=(%2.1f, %2.1f, %2.1f), dot=%2.3f\n",
              vno, pial ? "pial" : "white", dx, dy, dz, lx, ly, lz, v->nx, v->ny, v->nz, lx*v->nx+ly*v->ny+lz*v->nz) ;
  }
  

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeTangentPlanes(MRI_SURFACE *mris)
{
  VECTOR  *v_n, *v_e1, *v_e2, *v ;
  int     vno ;
  VERTEX  *vertex ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    /* now find some other non-parallel vector */
#if 0
    if (!FZERO(vertex->nx) || !FZERO(vertex->ny))
    {VECTOR_LOAD(v, 0.0, 0.0, 1.0) ; }
    else
    {VECTOR_LOAD(v, 0.0, 1.0, 0.0) ; }
#else
    VECTOR_LOAD(v, vertex->ny, vertex->nz, vertex->nx) ;
#endif
    V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    if ((V3_LEN_IS_ZERO(v_e1)))  /* happened to pick a parallel vector */
    {
      VECTOR_LOAD(v, vertex->ny, -vertex->nz, vertex->nx) ;
      V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    }

    if ((V3_LEN_IS_ZERO(v_e1)) && DIAG_VERBOSE_ON)  /* happened to pick a parallel vector */
      fprintf(stderr, "vertex %d: degenerate tangent plane\n", vno) ;
    V3_CROSS_PRODUCT(v_n, v_e1, v_e2) ;
    V3_NORMALIZE(v_e1, v_e1) ;
    V3_NORMALIZE(v_e2, v_e2) ;
    vertex->e1x = V3_X(v_e1) ; vertex->e2x = V3_X(v_e2) ;
    vertex->e1y = V3_Y(v_e1) ; vertex->e2y = V3_Y(v_e2) ;
    vertex->e1z = V3_Z(v_e1) ; vertex->e2z = V3_Z(v_e2) ;
  }

  VectorFree(&v) ;
  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScomputeMeanCurvature(MRI_SURFACE *mris)
{
  VECTOR  *v_n, *v_e1, *v_e2, *v_i ;
  int     vno, i, N ;
  VERTEX  *vertex, *vnb ;
  float   rsq, z, H, u, v, Hmin, Hmax ;

  mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_i = VectorAlloc(3, MATRIX_REAL) ;

  Hmin = 10000.0f ; Hmax = -Hmin ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z) ;
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z) ;

    H = 0.0 ; N = 0 ;
    for (i = 0 ; i < vertex->vnum ; i++)  /* for each neighbor */
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag)
        continue ;
      VECTOR_LOAD(v_i, vnb->x-vertex->x, vnb->y-vertex->y, vnb->z-vertex->z) ;
      
      /* calculate projection onto tangent plane */
      u = V3_DOT(v_i, v_e1) ; v = V3_DOT(v_i, v_e2) ;
      rsq = u*u + v*v ;
      z = V3_DOT(v_i, v_n) ;   /* height above tangent plane */
      if (!FZERO(rsq))
      {
        H += z / rsq ;
        N++ ;
        if ((fabs(z/rsq) > 5.0) && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
          fprintf(stdout, "%d --> %d: curvature = %2.1f\n", vno, vertex->v[i],
                  z/rsq) ;
      }
    }
    if (N > 0)
      vertex->H = 0.5f * H / (float)N ;
    else
      vertex->H = 0 ;
    if (vertex->H > Hmax)
      Hmax = vertex->H ;
    if (vertex->H < Hmin)
      Hmin = vertex->H ;
  }
  
  mris->min_curv = Hmin ; mris->max_curv = Hmax ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "mean curvature range [%2.3f --> %2.3f]\n",Hmin, Hmax) ;
  VectorFree(&v_i) ;
  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 1
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Fit a 1-d quadratic to the surface locally and move the
          vertex in the normal direction to improve the fit.
------------------------------------------------------*/
static int
mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv)
{
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n ;
  VERTEX   *v, *vn ;
  float    ui, vi, rsq, a, b ;

  if (FZERO(l_curv))
    return(NO_ERROR) ;

  mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(2, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v_Y = VectorAlloc(v->vtotal, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(v->vtotal, 2, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;
    for (n = 0 ; n < v->vtotal ; n++)  /* build data matrices */
    {
      vn = &mris->vertices[v->v[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ; vi = V3_DOT(v_e2, v_nbr) ; 
      rsq = ui*ui + vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsq ;
      *MATRIX_RELT(m_R, n+1, 2) = 1 ;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL) ;
    if (!m_R_inv)
    {
      MatrixFree(&m_R) ; VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ; b *= l_curv ;
    v->dx += b * v->nx ; v->dy += b * v->ny ; v->dz += b * v->nz ; 

    if (vno == Gdiag_no)
      fprintf(stdout, "v %d curvature term:      (%2.3f, %2.3f, %2.3f), "
              "a=%2.2f, b=%2.1f\n",
              vno, b*v->nx, b*v->ny, b*v->nz, a, b) ;
    MatrixFree(&m_R) ; VectorFree(&v_Y) ; MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ; VectorFree(&v_e1) ; VectorFree(&v_e2) ; 
  VectorFree(&v_nbr) ; VectorFree(&v_A) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Fit a 1-d quadratic to the surface locally and move the
          vertex in the normal direction to improve the fit.
------------------------------------------------------*/
static double
mrisComputeQuadraticCurvatureSSE(MRI_SURFACE *mris, double l_curv)
{
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n ;
  VERTEX   *v, *vn ;
  float    ui, vi, rsq, a, b ;
  double   sse = 0.0 ;

  if (FZERO(l_curv))
    return(0.0) ;

  mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(2, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v_Y = VectorAlloc(v->vtotal, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(v->vtotal, 2, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;
    for (n = 0 ; n < v->vtotal ; n++)  /* build data matrices */
    {
      vn = &mris->vertices[v->v[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ; vi = V3_DOT(v_e2, v_nbr) ; 
      rsq = ui*ui + vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsq ;
      *MATRIX_RELT(m_R, n+1, 2) = 1 ;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL) ;
    if (!m_R_inv)
    {
      MatrixFree(&m_R) ; VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ;
    sse += b*b ;
    if (vno == Gdiag_no)
      printf("v %d: curvature sse %2.2f\n", vno, b*b) ;
    MatrixFree(&m_R) ; VectorFree(&v_Y) ; MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ; VectorFree(&v_e1) ; VectorFree(&v_e2) ; 
  VectorFree(&v_nbr) ; VectorFree(&v_A) ;
  return(sse) ;
}
#else
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Fit a quadric to the surface locally and move the
          vertex in the normal direction to improve the fit.
------------------------------------------------------*/
static int
mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv)
{
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n ;
  VERTEX   *v, *vn ;
  float    ui, vi, a, b, c ;

  /* will set tangent plane basis to be principal directions */
  MRIScomputeSecondFundamentalForm(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v_Y = VectorAlloc(v->vtotal, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(v->vtotal, 3, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;
    for (n = 0 ; n < v->vtotal ; n++)  /* build data matrices */
    {
      vn = &mris->vertices[v->v[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ; vi = V3_DOT(v_e2, v_nbr) ; 
      *MATRIX_RELT(m_R, n+1, 1) = ui*ui ;
      *MATRIX_RELT(m_R, n+1, 2) = vi*vi ;
      *MATRIX_RELT(m_R, n+1, 3) = 1 ;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL) ;
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ; 
    c = VECTOR_ELT(v_A, 3) ; c *= l_curv ;
    v->dx += c * v->nx ; v->dy += c * v->ny ; v->dz += c * v->nz ; 

    MatrixFree(&m_R) ; VectorFree(&v_Y) ; MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ; VectorFree(&v_e1) ; VectorFree(&v_e2) ; 
  VectorFree(&v_nbr) ; VectorFree(&v_A) ;
  return(NO_ERROR) ;
}
#endif
#if 0
static int
mrisAverageDs(MRI_SURFACE *mris, int num_avgs)
{
  VERTEX  *v, *vn ;
  double  d, vnum ;
  int     n, vno, i, marked ;

  /* now average them */
  for (i = 0 ; i < num_avgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      marked = v->marked ;
      d = v->d ;
      for (vnum = 1.0, n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked != marked)
          continue ;
        vnum++ ;
        d += vn->d ;
      }
      if (vnum > 0.0)
        d /= vnum ;
      else
        d = 0.0 ;
      v->tdx = d ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->d = v->tdx ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define USE_TANGENTIAL_TERM  0
#define SCALE_BY_N           1
static int
mrisComputeCurvatureTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_n, *v_delta ;
  int     vno ;
  VERTEX  *vertex ;
  double  l_curv, deltaH, Hdesired ;

  l_curv = parms->l_curv ;
  if (FZERO(l_curv))
    return(NO_ERROR) ;

  Hdesired = parms->Hdesired ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;

    deltaH = (Hdesired - vertex->H) ;  /* sign will be reversed in MUL */
#if 1
#define TAN_SCALE 0.5
    deltaH = tanh(deltaH * TAN_SCALE) ;
#endif
    V3_SCALAR_MUL(v_n, -deltaH, v_delta) ;

    vertex->dx += l_curv * V3_X(v_delta) ; 
    vertex->dy += l_curv * V3_Y(v_delta) ; 
    vertex->dz += l_curv * V3_Z(v_delta) ;
    if (vno == Gdiag_no && DIAG_VERBOSE_ON)
      fprintf(stdout, "Hdes=%2.3f, dH = %2.3f, tanh= %2.3f, dx=%2.3f\n",
              Hdesired,Hdesired - vertex->H, deltaH, vertex->dx) ;
  }
  
  VectorFree(&v_delta) ;
  VectorFree(&v_n) ;

  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisComputeSethianCurvatureTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_n, *v_delta ;
  int     vno ;
  VERTEX  *vertex ;
  double  l_curv, deltaH, Hdesired ;

  Hdesired = parms->Hdesired ;
  l_curv = parms->l_curv ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;


    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    /* Osher and Sethian */
#if 1
    deltaH = (1 - parms->epsilon*tanh(Hdesired-vertex->H)*.5) ;
#else
    deltaH = (tanh(vertex->H-Hdesired)*.5) ;
#endif
    V3_SCALAR_MUL(v_n, deltaH, v_delta) ;

    vertex->dx += l_curv * V3_X(v_delta) ; 
    vertex->dy += l_curv * V3_Y(v_delta) ; 
    vertex->dz += l_curv * V3_Z(v_delta) ;
  }
  
  VectorFree(&v_delta) ;
  VectorFree(&v_n) ;

  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisComputeCurvatureGradientTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_n, *v_g1, *v_g2, *v_y, *v_delta, *v_tmp, *v_n_r2, *v_u_g1,*v_v_g2;
  int     vno, i, N ;
  VERTEX  *vertex, *vnb ;
  double  r2, r3, z, u, v, l_curv, deltaH, Hdesired ;
  FILE    *fp = NULL ;

  Hdesired = parms->Hdesired ;
  l_curv = parms->l_curv ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_g1 = VectorAlloc(3, MATRIX_REAL) ;
  v_g2 = VectorAlloc(3, MATRIX_REAL) ;
  v_y = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;
  v_tmp = VectorAlloc(3, MATRIX_REAL) ;
  v_n_r2 = VectorAlloc(3, MATRIX_REAL) ;
  v_u_g1 = VectorAlloc(3, MATRIX_REAL) ;
  v_v_g2 = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

/* 
   first compute term which comes from moving the this vertex in its
   own tangent plane.
 */
#if 1
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    VECTOR_LOAD(v_g1, vertex->e1x, vertex->e1y, vertex->e1z) ;
    VECTOR_LOAD(v_g2, vertex->e2x, vertex->e2y, vertex->e2z) ;
    V3_CLEAR(v_delta) ;

    if (Gdiag & DIAG_SHOW && vno == Gdiag_no && DIAG_VERBOSE_ON)
    {
      char fname[STRLEN] ;

      sprintf(fname, "v%d_%d.m", vno, parms->t) ;
      fp = fopen(fname, "w") ;
      fprintf(fp, "U%d_%d = [", vno, parms->t) ;
    }


    deltaH = (Hdesired - vertex->H) ;
    N = 0 ;
    for (i = 0 ; i < vertex->v2num ; i++)  /* for each neighbor */
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag)
        continue ;
      VECTOR_LOAD(v_y, vnb->x-vertex->x, vnb->y-vertex->y, vnb->z-vertex->z) ;
      
      /* calculate projection onto tangent plane */
      u = V3_DOT(v_y, v_g1) ; v = V3_DOT(v_y, v_g2) ;
      r2 = u*u + v*v ; r3 = r2 * sqrt(r2) ;
      if (FZERO(r3))
        continue ;
      z = V3_DOT(v_y, v_n) ;                /* height above tangent plane */
      if (Gdiag & DIAG_SHOW && vno == Gdiag_no && DIAG_VERBOSE_ON)
      {
        fprintf(fp, "%2.3f  %2.3f  %2.3f; ", u, v, z) ;
      }
      V3_SCALAR_MUL(v_n, -1.0/r2, v_n_r2) ;  /* -n/r^2 */
      V3_SCALAR_MUL(v_g1, u, v_u_g1) ;       /* ui g1 */
      V3_SCALAR_MUL(v_g2, v, v_v_g2) ;       /* vi g2 */
      V3_ADD(v_u_g1, v_v_g2, v_tmp) ;        /* ui g1 + vi g2 */
      V3_SCALAR_MUL(v_tmp, 4*z/r3, v_tmp) ;  /*  4 z / n^3 (ui g1 + vi g2) */
#if USE_TANGENTIAL_TERM
      V3_ADD(v_tmp, v_delta, v_delta) ;      /* add it into total delta */
#endif
      V3_ADD(v_n_r2, v_delta, v_delta) ;     /* add it into total delta */
      N++ ;
    }
    if (N > 0)
#if SCALE_BY_N
      V3_SCALAR_MUL(v_delta, deltaH * 2.0/N, v_delta) ;
#else
      V3_SCALAR_MUL(v_delta, deltaH * 2.0, v_delta) ;
#endif


#endif
      
      if (Gdiag & DIAG_SHOW && vno == Gdiag_no)
      {
        fprintf(fp, "] ;\n") ;
        fclose(fp) ;
      }

/* 
   now add terms which come from this vertex's appearance in
   neighboring tangent planes.
   */
    for (i = 0 ; i < vertex->v2num ; i++)  /* for each neighbor */
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag || !vnb->v2num)
        continue ;

      /* load coordinate system for neighbor's tangent plane */
      VECTOR_LOAD(v_n, vnb->nx, vnb->ny, vnb->nz) ;
      VECTOR_LOAD(v_g1, vnb->e1x, vnb->e1y, vnb->e1z) ;
      VECTOR_LOAD(v_g2, vnb->e2x, vnb->e2y, vnb->e2z) ;

      deltaH = (Hdesired - vnb->H) ;
      VECTOR_LOAD(v_y, vertex->x-vnb->x, vertex->y-vnb->y, vertex->z-vnb->z) ;

      /* calculate projection onto tangent plane */
      u = V3_DOT(v_y, v_g1) ; v = V3_DOT(v_y, v_g2) ;
      r2 = u*u + v*v ; r3 = r2 * sqrt(r2) ;
      if (FZERO(r3))
        continue ;
      z = V3_DOT(v_y, v_n) ;                 /* height above tangent plane */
      V3_SCALAR_MUL(v_n, 1.0/r2, v_n_r2) ;   /* n/r^2 */
      V3_SCALAR_MUL(v_g1, u, v_u_g1) ;       /* ui g1 */
      V3_SCALAR_MUL(v_g2, v, v_v_g2) ;       /* vi g2 */
      V3_ADD(v_u_g1, v_v_g2, v_tmp) ;        /* ui g1 + vi g2 */

#if USE_TANGENTIAL_TERM
      V3_SCALAR_MUL(v_tmp, -4*z/r3, v_tmp) ;  /*  -4z / n^3 (ui g1 + vi g2) */
      V3_ADD(v_n_r2, v_tmp, v_tmp) ;
#else
      V3_SCALAR_MUL(v_n_r2, 1.0, v_tmp) ;
#endif
#if SCALE_BY_N
      V3_SCALAR_MUL(v_tmp, deltaH*2/(double)vnb->v2num, v_tmp) ;
#else
      V3_SCALAR_MUL(v_tmp, deltaH*2, v_tmp) ;
#endif
      V3_ADD(v_tmp, v_delta, v_delta) ;
    }

    vertex->dx += l_curv * V3_X(v_delta) ; 
    vertex->dy += l_curv * V3_Y(v_delta) ; 
    vertex->dz += l_curv * V3_Z(v_delta) ;
    if (vno == Gdiag_no)
      fprintf(stdout, "moving v %d by (%2.3f, %2.3f, %2.3f)\n",
              vno, vertex->dx, vertex->dy, vertex->dz) ;
  }
  
  VectorFree(&v_tmp) ;
  VectorFree(&v_delta) ;
  VectorFree(&v_y) ;
  VectorFree(&v_n) ;
  VectorFree(&v_g1) ;
  VectorFree(&v_g2) ;
  VectorFree(&v_n_r2) ;
  VectorFree(&v_u_g1) ;
  VectorFree(&v_v_g2) ;

  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISupdateSurface(MRI_SURFACE *mris)
{
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScomputeEulerNumber(MRI_SURFACE *mris, int *pnvertices, 
                                    int *pnfaces, int *pnedges)
{
  int     eno, nfaces, nedges, nvertices, vno, fno, vnb, i ;
  VERTEX  *v1 ;

  /*  MRISripFaces(mris) ;*/
  for (nfaces = fno = 0 ; fno < mris->nfaces ; fno++)
    if (!mris->faces[fno].ripflag)
      nfaces++ ;

  for (nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
    if (!mris->vertices[vno].ripflag)
      nvertices++ ;

  for (nedges = vno = 0 ; vno < mris->nvertices ; vno++)
    if (!mris->vertices[vno].ripflag)
    {
      v1 = &mris->vertices[vno] ;
      for (i = 0 ; i < v1->vnum ; i++)
      {
        vnb = v1->v[i] ;
        /* already counted */
        if ((vnb > vno) && !mris->vertices[vnb].ripflag)
          nedges++ ;
      }
    }

  *pnfaces = nfaces ;
  *pnvertices = nvertices ;
  *pnedges = nedges ; 
  eno = nvertices - nedges + nfaces ;
  return(eno) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIStopologicalDefectIndex(MRI_SURFACE *mris)
{
  int     eno, nfaces, nedges, nvertices, vno, fno, vnb, i, dno ;
  VERTEX  *v1 ;

  for (nfaces = fno = 0 ; fno < mris->nfaces ; fno++)
    if (!mris->faces[fno].ripflag)
      nfaces++ ;

  for (nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
    if (!mris->vertices[vno].ripflag)
      nvertices++ ;

  for (nedges = vno = 0 ; vno < mris->nvertices ; vno++)
    if (!mris->vertices[vno].ripflag)
    {
      v1 = &mris->vertices[vno] ;
      for (i = 0 ; i < v1->vnum ; i++)
      {
        vnb = v1->v[i] ;
        /* already counted */
        if ((vnb > vno) && !mris->vertices[vnb].ripflag)
          nedges++ ;
      }
    }
#if 0
  nedges += nfaces ;        /* one additional edge added for each triangle */
  nfaces *= 2 ;             /* two triangular faces per face */
#endif
  eno = nvertices - nedges + nfaces ;

  dno = abs(2-eno) + abs(2*nedges-3*nfaces) ;
  return(dno) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISremoveTopologicalDefects(MRI_SURFACE *mris,float curv_thresh)
{
  VECTOR  *v_n, *v_e1, *v_e2, *v_i ;
  int     vno, i ;
  VERTEX  *vertex, *vnb ;
  float   rsq, z, u, v ;

  mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_i = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z) ;
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z) ;

    for (i = 0 ; i < vertex->vnum ; i++)  /* for each neighbor */
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag)
        continue ;
      VECTOR_LOAD(v_i, vnb->x-vertex->x, vnb->y-vertex->y, vnb->z-vertex->z) ;
      
      /* calculate projection onto tangent plane */
      u = V3_DOT(v_i, v_e1) ; v = V3_DOT(v_i, v_e2) ;
      rsq = u*u + v*v ;
      z = V3_DOT(v_i, v_n) ;   /* height above tangent plane */
      if (!FZERO(rsq))
      {
        if (fabs(z/rsq) > curv_thresh)
          mrisRemoveLink(mris, vno, vertex->v[i--]) ;
      }
    }
  }
  

  /*  MRISripFaces(mris) ;*/
  VectorFree(&v_i) ;
  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisRemoveLink(MRI_SURFACE *mris, int vno1, int vno2)
{
  FACE    *face ;
  int     vno, fno, nvalid ;
  VERTEX  *v1, *v2 ;

  v1 = &mris->vertices[vno1] ; v2 = &mris->vertices[vno2] ;
  mrisRemoveEdge(mris, vno1, vno2) ;
  mrisRemoveEdge(mris, vno2, vno1) ;

  /* now remove all the faces which contain both edges */
  for (fno = 0 ; fno < v1->num ; fno++)
  {
    face = &mris->faces[v1->f[fno]] ;
    for (vno = 0 ; vno < VERTICES_PER_FACE ; vno++)
      if (face->v[vno] == vno2)   /* found a face with both vertices */
        face->ripflag = 1 ;
    if (face->ripflag)
      mrisRemoveFace(mris, v1->f[fno]) ;
  }

  /* make sure vno1 and vno2 are still part of at least 1 valid face */
  for (nvalid = fno = 0 ; fno < v1->num ; fno++)
    if (!mris->faces[v1->f[fno]].ripflag)
      nvalid++ ;
  if (nvalid <= 0)
    v1->ripflag = 1 ;

  for (nvalid = fno = 0 ; fno < v2->num ; fno++)
    if (!mris->faces[v2->f[fno]].ripflag)
      nvalid++ ;
  if (nvalid <= 0)
    v2->ripflag = 1 ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisRemoveEdge(MRI_SURFACE *mris, int vno1, int vno2)
{
  int     i ;
  VERTEX  *v1, *v2 ;

  v1 = &mris->vertices[vno1] ; v2 = &mris->vertices[vno2] ;
  for (i = 0 ; i < v1->vnum ; i++)
    if (v1->v[i] == vno2)
    {
      v1->vnum-- ;
      if (i < v1->vnum)  /* not the (previous) last index */
        memmove(&v1->v[i], &v1->v[i+1], (v1->vnum-i)*sizeof(v1->v[0])) ;
      return(NO_ERROR) ;
    }
  return(ERROR_BADPARM) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          remove this face, as well as the link two vertices if
          they only exist through this face.
------------------------------------------------------*/
static int
mrisRemoveFace(MRI_SURFACE *mris, int fno)
{
  int   vno1, vno2, vno ;
  FACE  *face ;
  
  face = &mris->faces[fno] ;
  face->ripflag = 1 ;
  for (vno = 0 ; vno < VERTICES_PER_FACE ; vno++)
  {
    vno1 = face->v[vno] ; 
    vno2 = face->v[vno < VERTICES_PER_FACE-1 ? vno+1 : 0] ;
    if (!mrisCountValidLinks(mris, vno1, vno2))
    {
      mrisRemoveEdge(mris, vno1, vno2) ;
      mrisRemoveEdge(mris, vno2, vno1) ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Go through all the (valid) faces in vno1 and see how man
           valid links to vno2 exists.
------------------------------------------------------*/
static int
mrisCountValidLinks(MRI_SURFACE *mris, int vno1, int vno2)
{
  int    nvalid, fno, vno ;
  FACE   *face ;
  VERTEX *v ;

  v = &mris->vertices[vno1] ;
  for (nvalid = fno = 0 ; fno < v->num ; fno++)
  {
    face = &mris->faces[v->f[fno]] ;
    if (face->ripflag)
      continue ;
    for (vno = 0 ; vno < VERTICES_PER_FACE ; vno++)
      if (face->v[vno] == vno1)
      {
        if ((vno == VERTICES_PER_FACE-1) && (face->v[0] == vno2))
          nvalid++ ;
        else
          if (face->v[vno+1] == vno2)
            nvalid++ ;
      }
    
  }
  return(nvalid) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define ILL_CONDITIONED   500000.0
int
MRIScomputeSecondFundamentalForm(MRI_SURFACE *mris)
{
  int    vno, i, n, vmax, nbad = 0 ;
  VERTEX *vertex, *vnb ;
  MATRIX *m_U, *m_Ut, *m_tmp1, *m_tmp2, *m_inverse, *m_eigen, *m_Q ;
  VECTOR *v_c, *v_z, *v_n, *v_e1, *v_e2, *v_yi ;
  float  k1, k2, evalues[3], a11, a12, a21, a22, cond_no, kmax, kmin, rsq, k ;
  double ui, vi, total_area = 0.0, max_error ;

  if (mris->status == MRIS_PLANE)
    return(NO_ERROR) ;

  mrisComputeTangentPlanes(mris) ;

  v_c = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_yi = VectorAlloc(3, MATRIX_REAL) ;
  m_Q = MatrixAlloc(2, 2, MATRIX_REAL) ;   /* the quadratic form */
  m_eigen = MatrixAlloc(2, 2, MATRIX_REAL) ;

  mris->Kmin = mris->Hmin = 10000.0f ;
  mris->Kmax = mris->Hmax = -10000.0f ;
  mris->Ktotal = 0.0f ;
  vmax = -1 ; max_error = -1.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;

    if (vno == 142915)
      DiagBreak() ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z) ;
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z) ;

    if (vertex->vtotal <= 0)
      continue ;
    
    m_U = MatrixAlloc(vertex->vtotal, 3, MATRIX_REAL) ;
    v_z = VectorAlloc(vertex->vtotal, MATRIX_REAL) ;

    if (vno == Gdiag_no)
      DiagBreak() ;

    /* fit a quadratic form to the surface at this vertex */
    kmin = 10000.0f ; kmax = -kmin ;
    for (n = i = 0 ; i < vertex->vtotal ; i++)
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag)
        continue ;
/* 
   calculate the projection of this vertex onto the local tangent plane 
*/
      VECTOR_LOAD(v_yi, vnb->x-vertex->x, vnb->y-vertex->y,vnb->z-vertex->z);
      ui = V3_DOT(v_yi, v_e1) ; vi = V3_DOT(v_yi, v_e2) ;
      *MATRIX_RELT(m_U, n+1, 1) = ui*ui ;
      *MATRIX_RELT(m_U, n+1, 2) = 2*ui*vi ;
      *MATRIX_RELT(m_U, n+1, 3) = vi*vi ;
      VECTOR_ELT(v_z, n+1) = V3_DOT(v_n, v_yi) ;  /* height above TpS */
      rsq = ui*ui + vi*vi ;
      if (!FZERO(rsq))
      {
        k = VECTOR_ELT(v_z, n+1) / rsq ;
        if (k > kmax)
          kmax = k ;
        if (k < kmin)
          kmin = k ;
      }
      n++ ;
    }

    m_Ut = MatrixTranspose(m_U, NULL) ;          /* Ut */
    m_tmp2 = MatrixMultiply(m_Ut, m_U, NULL) ;   /* Ut U */
    cond_no = MatrixConditionNumber(m_tmp2) ;
#if 0
    m_inverse = MatrixInverse(m_tmp2, NULL) ;    /* (Ut U)^-1 */
#else
    m_inverse = MatrixSVDInverse(m_tmp2, NULL) ;    /* (Ut U)^-1 */
#endif
    if (!m_inverse)   /* singular matrix - must be planar?? */
    {
      nbad++ ;
      evalues[0] = evalues[1] = 0.0 ;
    }
    else
    {
      m_tmp1 = MatrixMultiply(m_Ut, v_z, NULL) ;   /* Ut z */
      MatrixMultiply(m_inverse, m_tmp1, v_c) ;     /* (Ut U)^-1 Ut z */

      /* now build Hessian matrix */
      *MATRIX_RELT(m_Q,1,1) = 2*VECTOR_ELT(v_c, 1) ;
      *MATRIX_RELT(m_Q,1,2) = *MATRIX_RELT(m_Q,2,1) = 2*VECTOR_ELT(v_c, 2) ;
      *MATRIX_RELT(m_Q,2,2) = 2*VECTOR_ELT(v_c, 3) ;

      if (cond_no >= ILL_CONDITIONED)
      {
#if 0
        MatrixSVDEigenValues(m_Q, evalues) ;
        vertex->k1 = k1 = evalues[0] ;
        vertex->k2 = k2 = evalues[1] ;
#else
        vertex->k1 = k1 = kmax ;
        vertex->k2 = k2 = kmin ;
#endif
        vertex->K = k1*k2 ; vertex->H = (k1+k2)/2 ;
        MatrixFree(&m_Ut) ;
        MatrixFree(&m_tmp2) ;
        MatrixFree(&m_U) ;
        VectorFree(&v_z) ;
        MatrixFree(&m_tmp1) ;
        MatrixFree(&m_inverse) ;
        continue ;
      }

      /* the columns of m_eigen will be the eigenvectors of m_Q */
      if (MatrixEigenSystem(m_Q, evalues, m_eigen) == NULL)
      {
        nbad++ ;
        MatrixSVDEigenValues(m_Q, evalues) ;
        vertex->k1 = k1 = evalues[0] ;
        vertex->k2 = k2 = evalues[1] ;
        vertex->K = k1*k2 ; vertex->H = (k1+k2)/2 ;
        MatrixFree(&m_Ut) ;
        MatrixFree(&m_tmp2) ;
        MatrixFree(&m_U) ;
        VectorFree(&v_z) ;
        MatrixFree(&m_tmp1) ;
        MatrixFree(&m_inverse) ;
        continue ;
      }

      MatrixFree(&m_tmp1) ;
      MatrixFree(&m_inverse) ;
    }
    k1 = evalues[0] ; k2 = evalues[1] ;
    vertex->k1 = k1 ; vertex->k2 = k2 ;
    vertex->K = k1 * k2 ;
    vertex->H = (k1 + k2) / 2 ;
    if (vno == Gdiag_no && (Gdiag & DIAG_SHOW))
      fprintf(stdout, "v %d: k1=%2.3f, k2=%2.3f, K=%2.3f, H=%2.3f\n",
              vno, vertex->k1, vertex->k2, vertex->K, vertex->H) ;
    if (vertex->K < mris->Kmin)
      mris->Kmin = vertex->K ;
    if (vertex->H < mris->Hmin)
      mris->Hmin = vertex->H ;
    if (vertex->K > mris->Kmax)
      mris->Kmax = vertex->K ;
    if (vertex->H > mris->Hmax)
      mris->Hmax = vertex->H ;
    mris->Ktotal += (double)k1 * (double)k2 * (double)vertex->area ;
    total_area += (double)vertex->area ;

    /* now update the basis vectors to be the principal directions */
    a11 = *MATRIX_RELT(m_eigen,1,1) ; a12 = *MATRIX_RELT(m_eigen,1,2) ;
    a21 = *MATRIX_RELT(m_eigen,2,1) ; a22 = *MATRIX_RELT(m_eigen,2,2) ;
    vertex->e1x = V3_X(v_e1) * a11 + V3_X(v_e2) * a21 ;
    vertex->e1y = V3_Y(v_e1) * a11 + V3_Y(v_e2) * a21 ;
    vertex->e1z = V3_Z(v_e1) * a11 + V3_Z(v_e2) * a21 ;
    vertex->e2x = V3_X(v_e1) * a12 + V3_X(v_e2) * a22 ;
    vertex->e2y = V3_Y(v_e1) * a12 + V3_Y(v_e2) * a22 ;
    vertex->e2z = V3_Z(v_e1) * a12 + V3_Z(v_e2) * a22 ;

    MatrixFree(&m_Ut) ;
    MatrixFree(&m_tmp2) ;
    MatrixFree(&m_U) ;
    VectorFree(&v_z) ;

  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "max H error=%2.3f at %d\n", max_error, vmax) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "total area = %2.3f\n", total_area);

  if (Gdiag & DIAG_SHOW && (nbad > 0))
    fprintf(stdout, "%d ill-conditioned points\n", nbad) ;
  MatrixFree(&m_eigen) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_c) ;
  VectorFree(&v_n) ;
  VectorFree(&v_yi) ;
  MatrixFree(&m_Q) ;
  return(NO_ERROR) ;
}
int
MRIScomputeSecondFundamentalFormAtVertex(MRI_SURFACE *mris, int vno, 
                                         int *vertices, int vnum)
{
  int    i, n, nbad = 0 ;
  VERTEX *vertex, *vnb ;
  MATRIX *m_U, *m_Ut, *m_tmp1, *m_tmp2, *m_inverse ;
  VECTOR *v_z ;
  static MATRIX *m_Q, *m_eigen ;
  static VECTOR *v_c = NULL, *v_n, *v_e1, *v_e2, *v_yi ;
  float  k1, k2, evalues[3], a11, a12, a21, a22, cond_no, rsq, k, kmin, kmax ;
  double ui, vi ;

  if (mris->status == MRIS_PLANE)
    return(NO_ERROR) ;

  if (v_c == NULL)
  {
    v_c = VectorAlloc(3, MATRIX_REAL) ;
    v_n = VectorAlloc(3, MATRIX_REAL) ;
    v_e1 = VectorAlloc(3, MATRIX_REAL) ;
    v_e2 = VectorAlloc(3, MATRIX_REAL) ;
    v_yi = VectorAlloc(3, MATRIX_REAL) ;
    m_Q = MatrixAlloc(2, 2, MATRIX_REAL) ;   /* the quadratic form */
    m_eigen = MatrixAlloc(2, 2, MATRIX_REAL) ;
  }

  vertex = &mris->vertices[vno] ;
  if (vertex->ripflag)
    return(ERROR_BADPARM) ;
  
  if (vno == 142915)
    DiagBreak() ;
  VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
  VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z) ;
  VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z) ;
  
  if (vnum <= 0)
    return(ERROR_BADPARM) ;
  
  m_U = MatrixAlloc(vnum, 3, MATRIX_REAL) ;
  v_z = VectorAlloc(vnum, MATRIX_REAL) ;
  
  if (vno == Gdiag_no)
    DiagBreak() ;
  
  /* fit a quadratic form to the surface at this vertex */
  kmin = 10000.0f ; kmax = -kmin ;
  for (n = i = 0 ; i < vnum ; i++)
  {
    vnb = &mris->vertices[vertices[i]] ;
    if (vnb->ripflag)
      continue ;
    /* 
       calculate the projection of this vertex onto the local tangent plane 
    */
    VECTOR_LOAD(v_yi, vnb->x-vertex->x, vnb->y-vertex->y,vnb->z-vertex->z);
    ui = V3_DOT(v_yi, v_e1) ; vi = V3_DOT(v_yi, v_e2) ;
    *MATRIX_RELT(m_U, n+1, 1) = ui*ui ;
    *MATRIX_RELT(m_U, n+1, 2) = 2*ui*vi ;
    *MATRIX_RELT(m_U, n+1, 3) = vi*vi ;
    VECTOR_ELT(v_z, n+1) = V3_DOT(v_n, v_yi) ;  /* height above TpS */
    rsq = ui*ui + vi*vi ;
    if (!FZERO(rsq))
    {
      k = VECTOR_ELT(v_z, n+1) / rsq ;
      if (k > kmax)
        kmax = k ;
      if (k < kmin)
        kmin = k ;
    }
    n++ ;
  }
  
  m_Ut = MatrixTranspose(m_U, NULL) ;          /* Ut */
  m_tmp2 = MatrixMultiply(m_Ut, m_U, NULL) ;   /* Ut U */
  cond_no = MatrixConditionNumber(m_tmp2) ;
#if 0
  m_inverse = MatrixInverse(m_tmp2, NULL) ;    /* (Ut U)^-1 */
#else
  m_inverse = MatrixSVDInverse(m_tmp2, NULL) ;    /* (Ut U)^-1 */
#endif
  if (!m_inverse)   /* singular matrix - must be planar?? */
  {
    nbad++ ;
    evalues[0] = evalues[1] = 0.0 ;
  }
  else
  {
    m_tmp1 = MatrixMultiply(m_Ut, v_z, NULL) ;   /* Ut z */
    MatrixMultiply(m_inverse, m_tmp1, v_c) ;     /* (Ut U)^-1 Ut z */
    
    /* now build Hessian matrix */
    *MATRIX_RELT(m_Q,1,1) = 2*VECTOR_ELT(v_c, 1) ;
    *MATRIX_RELT(m_Q,1,2) = *MATRIX_RELT(m_Q,2,1) = 2*VECTOR_ELT(v_c, 2) ;
    *MATRIX_RELT(m_Q,2,2) = 2*VECTOR_ELT(v_c, 3) ;
    
    if (cond_no >= ILL_CONDITIONED)
    {
#if 0
      MatrixSVDEigenValues(m_Q, evalues) ;
      vertex->k1 = k1 = evalues[0] ;
      vertex->k2 = k2 = evalues[1] ;
#else
      vertex->k1 = k1 = kmax ;
      vertex->k2 = k2 = kmin ;
#endif
      vertex->K = k1*k2 ; vertex->H = (k1+k2)/2 ;
      MatrixFree(&m_Ut) ;
      MatrixFree(&m_tmp2) ;
      MatrixFree(&m_U) ;
      VectorFree(&v_z) ;
      MatrixFree(&m_tmp1) ;
      MatrixFree(&m_inverse) ;
      return(ERROR_BADPARM) ;
    }
    
    /* the columns of m_eigen will be the eigenvectors of m_Q */
    if (MatrixEigenSystem(m_Q, evalues, m_eigen) == NULL)
    {
      nbad++ ;
      MatrixSVDEigenValues(m_Q, evalues) ;
      vertex->k1 = k1 = evalues[0] ;
      vertex->k2 = k2 = evalues[1] ;
      vertex->K = k1*k2 ; vertex->H = (k1+k2)/2 ;
      MatrixFree(&m_Ut) ;
      MatrixFree(&m_tmp2) ;
      MatrixFree(&m_U) ;
      VectorFree(&v_z) ;
      MatrixFree(&m_tmp1) ;
      MatrixFree(&m_inverse) ;
      return(ERROR_BADPARM) ;
    }
    
    MatrixFree(&m_tmp1) ;
    MatrixFree(&m_inverse) ;
  }
  k1 = evalues[0] ; k2 = evalues[1] ;
  vertex->k1 = k1 ; vertex->k2 = k2 ;
  vertex->K = k1 * k2 ;
  vertex->H = (k1 + k2) / 2 ;
  if (vno == Gdiag_no && (Gdiag & DIAG_SHOW))
    fprintf(stdout, "v %d: k1=%2.3f, k2=%2.3f, K=%2.3f, H=%2.3f\n",
            vno, vertex->k1, vertex->k2, vertex->K, vertex->H) ;
  if (vertex->K < mris->Kmin)
    mris->Kmin = vertex->K ;
  if (vertex->H < mris->Hmin)
    mris->Hmin = vertex->H ;
  if (vertex->K > mris->Kmax)
    mris->Kmax = vertex->K ;
  if (vertex->H > mris->Hmax)
    mris->Hmax = vertex->H ;
  mris->Ktotal += (double)k1 * (double)k2 * (double)vertex->area ;
  
  /* now update the basis vectors to be the principal directions */
  a11 = *MATRIX_RELT(m_eigen,1,1) ; a12 = *MATRIX_RELT(m_eigen,1,2) ;
  a21 = *MATRIX_RELT(m_eigen,2,1) ; a22 = *MATRIX_RELT(m_eigen,2,2) ;
  vertex->e1x = V3_X(v_e1) * a11 + V3_X(v_e2) * a21 ;
  vertex->e1y = V3_Y(v_e1) * a11 + V3_Y(v_e2) * a21 ;
  vertex->e1z = V3_Z(v_e1) * a11 + V3_Z(v_e2) * a21 ;
  vertex->e2x = V3_X(v_e1) * a12 + V3_X(v_e2) * a22 ;
  vertex->e2y = V3_Y(v_e1) * a12 + V3_Y(v_e2) * a22 ;
  vertex->e2z = V3_Z(v_e1) * a12 + V3_Z(v_e2) * a22 ;
  
  MatrixFree(&m_Ut) ;
  MatrixFree(&m_tmp2) ;
  MatrixFree(&m_U) ;
  VectorFree(&v_z) ;

  if (Gdiag & DIAG_SHOW && (nbad > 0))
    fprintf(stdout, "%d ill-conditioned points\n", nbad) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISuseAreaErrors(MRI_SURFACE *mris)
{
  int    vno, fno, n ;
  float  area, orig_area ;
  FACE   *face ;
  VERTEX *vertex ;

  for (vno = 0 ; vno < mris->nvertices; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    area = orig_area = 0.0f ;

    /* use average area of all faces this vertex is part of -
       this is not really correct, but should be good enough for
       visualization purposes.
       */
    for (n = fno = 0 ; fno < vertex->num ; fno++)
    {
      face = &mris->faces[vertex->f[fno]] ;
      area += face->area ;
      orig_area += face->orig_area ;
    }
    vertex->curv = (area-orig_area) / (float)n ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISuseGaussianCurvature(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *vertex ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->curv = vertex->K ;
  }

  mris->min_curv = mris->Kmin ; mris->max_curv = mris->Kmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISuseMeanCurvature(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *vertex ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->curv = vertex->H ;
  }

  mris->min_curv = mris->Hmin ; mris->max_curv = mris->Hmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISusePrincipalCurvature(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *vertex ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->curv = vertex->k1 ;
  }

  mris->min_curv = mris->Hmin ; mris->max_curv = mris->Hmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the folding (FI) and intrinsic curvature (ICI)
          indices of the surface.

          Drury and Van Essen:  ICI = 56  (lh) and 54  (rh)
                                 FI = 500 (lh) and 520 (rh)
------------------------------------------------------*/
int
MRIScomputeCurvatureIndices(MRI_SURFACE *mris, double *pici, double *pfi)
{
  int    vno ;
  VERTEX *vertex ;
  double k1, k2, ici, fi, area, Kmax, Kmin ;

  ici = fi = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    k1 = (double)vertex->k1 ; k2 = (double)vertex->k2 ;
    area = (double)vertex->area ;
    if (vertex->K > 0)
      ici += area * (double)vertex->K ;
    Kmax = (double)fabs(k1) ; Kmin = (double)fabs(k2) ;
    fi += area * Kmax * (Kmax - Kmin) ;
  }

  *pfi = fi / (4.0*M_PI) ; *pici = ici / (4.0 * M_PI) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRISmaxRadius(MRI_SURFACE *mris)
{
  double  radius ;
  int    vno, n ;
  VERTEX *vertex ;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0, r ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
#if 0
    if (vertex->ripflag)
      continue ;
#endif
    x = (double)vertex->x ; y = (double)vertex->y ; z = (double)vertex->z ;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  mris->xhi = xhi ; mris->xlo = xlo ; 
  mris->yhi = yhi ; mris->ylo = ylo ; 
  mris->zhi = zhi ; mris->zlo = zlo ; 
  x0 = (xlo+xhi)/2.0f ; y0 = (ylo+yhi)/2.0f ; z0 = (zlo+zhi)/2.0f ;
  for (radius = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
#if 0
    if (vertex->ripflag)
      continue ;
#endif
    n++ ;
    x = (double)vertex->x-x0 ; 
    y = (double)vertex->y-y0 ; 
    z = (double)vertex->z-z0 ;
    r = sqrt(x*x + y*y + z*z) ;
    if (r > radius)
      radius = r ;
  }

  return(radius) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeSurfaceDimensions(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *vertex ;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
#if 0
    if (vertex->ripflag)
      continue ;
#endif
    x = (double)vertex->x ; y = (double)vertex->y ; z = (double)vertex->z ;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  mris->xlo = xlo ; mris->xhi = xhi ;
  mris->ylo = ylo ; mris->yhi = yhi ;
  mris->zlo = zlo ; mris->zhi = zhi ;
  mris->xctr = (xlo+xhi)/2.0f ; 
  mris->yctr = (ylo+yhi)/2.0f ; 
  mris->zctr = (zlo+zhi)/2.0f ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRISaverageRadius(MRI_SURFACE *mris)
{
  double  radius ;
  int    vno, n ;
  VERTEX *vertex ;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0 ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
#if 0
    if (vertex->ripflag)
      continue ;
#endif
    x = (double)vertex->x ; y = (double)vertex->y ; z = (double)vertex->z ;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) 
      zhi=z;
    if (z<zlo) 
      zlo=z;
    if (zlo < -1000)
      DiagBreak() ;
  }
  x0 = (xlo+xhi)/2.0f ; y0 = (ylo+yhi)/2.0f ; z0 = (zlo+zhi)/2.0f ;
  for (radius = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
#if 0
    if (vertex->ripflag)
      continue ;
#endif
    n++ ;
    x = (double)vertex->x-x0 ; 
    y = (double)vertex->y-y0 ; 
    z = (double)vertex->z-z0 ;
    radius += sqrt(x*x + y*y + z*z) ;
  }

  radius /= (double)n ;
  return(radius) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Descripton
------------------------------------------------------*/
double
MRISaverageCanonicalRadius(MRI_SURFACE *mris)
{
  double  radius ;
  int    vno, n ;
  VERTEX *vertex ;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0 ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
#if 0
    if (vertex->ripflag)
      continue ;
#endif
    x = (double)vertex->cx ; y = (double)vertex->cy ; z = (double)vertex->cz ;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  x0 = (xlo+xhi)/2.0f ; y0 = (ylo+yhi)/2.0f ; z0 = (zlo+zhi)/2.0f ;
  for (radius = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
#if 0
    if (vertex->ripflag)
      continue ;
#endif
    n++ ;
    x = (double)vertex->cx-x0 ; 
    y = (double)vertex->cy-y0 ; 
    z = (double)vertex->cz-z0 ;
    radius += sqrt(x*x + y*y + z*z) ;
  }

  radius /= (double)n ;
  return(radius) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISinflateBrain(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     n_averages, n, write_iterations, niterations ;
  double  delta_t = 0.0, rms_height, desired_rms_height, sse, l_dist ;

  write_iterations = parms->write_iterations ;

  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s.out", parms->base_name) ;
    if (!parms->start_t)
      parms->fp = fopen(fname, "w") ;
    else
      parms->fp = fopen(fname, "a") ;
    if (!parms->fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
    Progname, fname) ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  MRIScomputeMetricProperties(mris) ;
  niterations = parms->niterations ;
  rms_height = MRISrmsTPHeight(mris) ;
  desired_rms_height = parms->desired_rms_height ;
  fprintf(stdout, "inflating to desired rms height = %2.3f\n", 
          desired_rms_height);

  /* write out initial surface */
  if (!parms->start_t && (parms->write_iterations > 0) && (Gdiag&DIAG_WRITE))
    mrisWriteSnapshot(mris, parms, 0) ;
  
  sse = MRIScomputeSSE(mris, parms) ;
  if (!parms->start_t)
  {
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout,"%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 
              0, 0.0f, (float)rms_height, parms->n_averages) ;
    else
      fprintf(stdout, "\rstep %3.3d: RMS=%2.3f (target=%2.3f)   ", 0, 
              rms_height, desired_rms_height);
    if (Gdiag & DIAG_WRITE)
    {
      fprintf(parms->fp, 
              "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 
              0, 0.0f, (float)rms_height, parms->n_averages) ;
      fflush(parms->fp) ;
    }
    
    MRISclearCurvature(mris) ;  /* curvature will be used to calculate sulc */
  }

  l_dist = parms->l_dist ;
  for (n_averages = parms->n_averages ; n_averages >= 0 ; n_averages /= 2)
  {
    parms->l_dist = l_dist * sqrt(n_averages) ;
    for (n = parms->start_t ; n < parms->start_t+niterations ; n++)
    {
      MRISclearGradient(mris) ;
      mrisComputeDistanceTerm(mris, parms) ;
      mrisComputeSphereTerm(mris, parms->l_sphere, parms->a) ;
      mrisComputeExpansionTerm(mris, parms->l_expand) ;
      
      mrisAverageGradients(mris, n_averages) ;
      mrisComputeNormalSpringTerm(mris, parms->l_nspring) ;
      mrisComputeTangentialSpringTerm(mris, parms->l_tspring) ;
      mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv) ;
      mrisComputeSpringTerm(mris, parms->l_spring) ;
      mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm) ;
      switch (parms->integration_type)
      {
      case INTEGRATE_LM_SEARCH:
        delta_t = mrisLineMinimizeSearch(mris, parms) ;
        break ;
      default:
      case INTEGRATE_LINE_MINIMIZE:
        delta_t = mrisLineMinimize(mris, parms) ;
        break ;
      case INTEGRATE_MOMENTUM:
        delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, 
                                       parms->tol, 0/*parms->n_averages*/) ;
        break ;
      case INTEGRATE_ADAPTIVE:
        mrisAdaptiveTimeStep(mris, parms);
        break ;
      }
      mrisTrackTotalDistance(mris) ;  /* update sulc */
      MRIScomputeMetricProperties(mris) ; 
      sse = MRIScomputeSSE(mris, parms) ;
      rms_height = MRISrmsTPHeight(mris) ;
      if (!((n+1) % 5))     /* print some diagnostics */
      {
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, 
                "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d, l_dist=%2.2f\n",
                  n+1,(float)delta_t, (float)rms_height, n_averages,
                  parms->l_dist);
        else
          fprintf(stdout, "\rstep %3.3d: RMS=%2.3f (target=%2.3f)   ", 
                  n+1, rms_height, desired_rms_height) ;
        if (Gdiag & DIAG_WRITE)
        {
          fprintf(parms->fp, 
                "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d, l_dist=%2.2f\n",
                  n+1,(float)delta_t, (float)rms_height, n_averages, 
                  parms->l_dist);
          fflush(parms->fp) ;
        }
      }

      if (parms->scale > 0)
        MRISscaleBrainArea(mris) ;
      if ((parms->write_iterations > 0) &&
          !((n+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
        mrisWriteSnapshot(mris, parms, n+1) ;
      if (rms_height < desired_rms_height)
        break ;
    }

    parms->start_t = n ;
    if (!n_averages || rms_height < desired_rms_height)
      break ;
  }

  fprintf(stdout, "\ninflation complete.\n") ;
  if (Gdiag & DIAG_WRITE)
  {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISinflateToSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     n_averages, n, write_iterations, niterations, base_averages ;
  double  delta_t = 0.0, rms_radial_error, sse, base_dt ;
  MHT     *mht_v_current = NULL ;

  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  write_iterations = parms->write_iterations ;
  n_averages = parms->n_averages ;

#if 1
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    
    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
            parms->base_name);
    if (!parms->fp)
    {
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      else
        parms->fp = fopen(fname, "a") ;

      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "MRISunfold: could not open log file %s\n",
                  fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris,parms) ;
  }
#else
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s.out", parms->base_name) ;
    if (!parms->start_t)
      parms->fp = fopen(fname, "w") ;
    else
      parms->fp = fopen(fname, "a") ;
    if (!parms->fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
    Progname, fname) ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
#endif
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  MRIScomputeMetricProperties(mris) ;
#if 0
  MRIScomputeSecondFundamentalForm(mris) ;
#endif
  /*  parms->start_t = 0 ;*/
  niterations = parms->niterations ;
  MRISstoreMetricProperties(mris) ;
  rms_radial_error = 
    sqrt(mrisComputeSphereError(mris, 1.0, parms->a)/mris->nvertices);
  fprintf(stdout, "inflating to sphere (rms error < %2.2f)\n",
          parms->tol*parms->a / 100.0f) ;

  /* write out initial surface */
  if (!parms->start_t && (parms->write_iterations > 0) && (Gdiag&DIAG_WRITE))
    mrisWriteSnapshot(mris, parms, 0) ;
  
  sse = MRIScomputeSSE(mris, parms) ;
  if (!parms->start_t)
  {
    fprintf(stdout,"%3.3d: dt: %2.4f, rms radial error=%2.3f, avgs=%d\n", 
            0, 0.0f, (float)rms_radial_error, n_averages) ;
    if (Gdiag & DIAG_WRITE)
    {
      fprintf(parms->fp,"%3.3d: dt: %2.4f, rms radial error=%2.3f, avgs=%d\n", 
            0, 0.0f, rms_radial_error, n_averages) ;
      fflush(parms->fp) ;
    }
    
    MRISclearCurvature(mris) ;   /* curvature will be used to calculate sulc */
  }

  base_averages = parms->n_averages ; base_dt = parms->dt ;
  for (n_averages = base_averages ; n_averages >= 0 ; n_averages /= 4)
  {
    parms->n_averages = n_averages ;
    /*    parms->dt = (sqrt((float)n_averages)+1)*base_dt ;*/
    for (n = parms->start_t ; n < parms->start_t+niterations ; n++)
    {
      if (!FZERO(parms->l_repulse_ratio))
        mht_v_current =
          MHTfillVertexTableRes(mris,mht_v_current,CURRENT_VERTICES, 3.0f);
      MRISclearGradient(mris) ;
      mrisComputeDistanceTerm(mris, parms) ;
      mrisComputeSphereTerm(mris, parms->l_sphere, parms->a) ;
      mrisComputeExpansionTerm(mris, parms->l_expand) ;
      mrisComputeRepulsiveRatioTerm(mris,parms->l_repulse_ratio,mht_v_current);
      mrisComputeConvexityTerm(mris, parms->l_convex) ;
      
      mrisAverageGradients(mris, n_averages) ;
      mrisComputeSpringTerm(mris, parms->l_spring) ;
      mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm) ;
      switch (parms->integration_type)
      {
      case INTEGRATE_LM_SEARCH:
        delta_t = mrisLineMinimizeSearch(mris, parms) ;
        break ;
      default:
      case INTEGRATE_LINE_MINIMIZE:
        delta_t = mrisLineMinimize(mris, parms) ;
        break ;
      case INTEGRATE_MOMENTUM:
        delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, 
                                       parms->tol, 0/*parms->n_averages*/) ;
        break ;
      case INTEGRATE_ADAPTIVE:
        mrisAdaptiveTimeStep(mris, parms);
        break ;
      }
      mrisTrackTotalDistance(mris) ;  /* update sulc */
      MRIScomputeMetricProperties(mris) ; 
      sse = MRIScomputeSSE(mris, parms) ;
      rms_radial_error = 
        sqrt(mrisComputeSphereError(mris, 1.0, parms->a)/mris->nvertices);
      if (!((n+1) % 5))     /* print some diagnostics */
      {
        fprintf(stdout, 
                "%3.3d: dt: %2.4f, rms radial error=%2.3f, avgs=%d\n", 
                n+1,(float)delta_t, (float)rms_radial_error, n_averages);
        if (Gdiag & DIAG_WRITE)
        {
          fprintf(parms->fp, 
                  "%3.3d: dt: %2.4f, rms radial error=%2.3f, avgs=%d\n", 
                  n+1,(float)delta_t, (float)rms_radial_error, n_averages);
          fflush(parms->fp) ;
        }
      }
      
      if ((parms->write_iterations > 0) &&
          !((n+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
        mrisWriteSnapshot(mris, parms, n+1) ;
      if (100.0*rms_radial_error/parms->a < parms->tol)
        break ;
    }
    parms->start_t = n ;
    if (!n_averages || (100.0*rms_radial_error/parms->a < parms->tol))
      break ;
  }

  fprintf(stdout, "\nspherical inflation complete.\n") ;
  if (Gdiag & DIAG_WRITE)
  {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }
  if (!FZERO(parms->l_repulse))
    MHTfree(&mht_v_current) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISapplyGradient(MRI_SURFACE *mris, double dt)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  MRISstoreCurrentPositions(mris) ;
  if (mris->status == MRIS_RIGID_BODY)
    MRISrotate(mris, mris, dt*mris->alpha, dt*mris->beta, dt*mris->gamma) ;
  else for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (!finite(v->x) || !finite(v->y) || !finite(v->z))
      ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n",vno) ;
    if (!finite(v->dx) || !finite(v->dy) || !finite(v->dz))
      ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n",vno) ;
    v->x += dt*v->dx ; v->y += dt*v->dy; v->z += dt*v->dz;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisClearDistances(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->d = 0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISclearGradient(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->dx = 0 ; v->dy = 0 ; v->dz = 0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisClearExtraGradient(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (mris->dx2)
      mris->dx2[vno] = mris->dy2[vno] = mris->dz2[vno] = 0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisClearMomentum(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->odx = 0 ; v->ody = 0 ; v->odz = 0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRIStotalVariationDifference(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;
  double  var, delta ;

  nvertices = mris->nvertices ;
  for (var = 0.0, vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    delta = fabs(v->k1 - v->k2) ; delta *= delta ;
    var += v->area * delta ;
    if (!finite(var))
      ErrorPrintf(ERROR_BADPARM, "curvature at vertex %d is not finite!\n",
                  vno) ;
  }
  return(var) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRIStotalVariation(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;
  double  var, delta ;

  nvertices = mris->nvertices ;
  for (var = 0.0, vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    delta = v->k1 * v->k1 + v->k2 * v->k2 ;
    var += v->area * delta ;
    if (!finite(var))
      ErrorPrintf(ERROR_BADPARM, "curvature at vertex %d is not finite!\n",
                  vno) ;
  }
  return(var) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRIScurvatureError(MRI_SURFACE *mris, double Kd)
{
  int     vno, nvertices, n ;
  VERTEX  *v ;
  double  Kerror, deltaK ;

  nvertices = mris->nvertices ;
  for (Kerror = 0.0, n = vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    deltaK = (v->k1-Kd) ; deltaK *= deltaK ; Kerror += v->area * deltaK ; n++;
    deltaK = (v->k2-Kd) ; deltaK *= deltaK ; Kerror += v->area * deltaK ; n++ ;
  }
/*  return(sqrt(Kerror/(double)n)) ;*/
  return(Kerror) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISrestoreOldPositions(MRI_SURFACE *mris)
{
#if 1
  return(MRISrestoreVertexPositions(mris, TMP_VERTICES)) ;
#else
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->x = v->ox ; v->y = v->oy ; v->z = v->oz ;
  }
  return(NO_ERROR) ;
#endif
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISstoreCurrentPositions(MRI_SURFACE *mris)
{
#if 1
  return(MRISsaveVertexPositions(mris, TMP_VERTICES)) ;
#else
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->ox = v->x ; v->oy = v->y ; v->oz = v->z ;
  }
  return(NO_ERROR) ;
#endif
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisStoreCurrentGradient(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->odx = v->dx ; v->ody = v->dy ; v->odz = v->dz ;
  }
  return(NO_ERROR) ;
}

#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISstoreMeanCurvature(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->H = v->curv ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISstoreMetricProperties(MRI_SURFACE *mris)
{
  int     vno, nvertices, fno, ano, n ;
  VERTEX  *v ;
  FACE    *f ;

#if 0
  MRIScomputeNormals(mris);              /* update vertex areas */
  MRIScomputeTriangleProperties(mris) ;  /* update triangle properties */
#endif
  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->origarea = v->area ;
#if 1
    if (v->dist && v->dist_orig)
      for (n = 0 ; n < v->vtotal ; n++)
        v->dist_orig[n] = v->dist[n] ;
#endif
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    f->orig_area = f->area ;
    for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      f->orig_angle[ano] = f->angle[ano] ;
  }
  mris->orig_area = mris->total_area ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISrestoreMetricProperties(MRI_SURFACE *mris)
{
  int     vno, nvertices, fno, ano, n ;
  VERTEX  *v ;
  FACE    *f ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->area = v->origarea ;
    for (n = 0 ; n < v->vtotal ; n++)
      v->dist[n] = v->dist_orig[n] ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    f->area = f->orig_area ;
    for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      f->angle[ano] = f->orig_angle[ano] ;
  }
  mris->total_area = mris->orig_area ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisValidVertices(MRI_SURFACE *mris)
{
  int vno, nvertices, nvalid ;

  nvertices = mris->nvertices ;
  for (vno = nvalid = 0 ; vno < nvertices ; vno++)
    if (!mris->vertices[vno].ripflag)
      nvalid++ ;

  return(nvalid) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisValidFaces(MRI_SURFACE *mris)
{
  int fno, nfaces, nvalid ;

  nfaces = mris->nfaces ;
  for (fno = nvalid = 0 ; fno < nfaces ; fno++)
    if (!mris->faces[fno].ripflag)
      nvalid++ ;

  return(nvalid) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
int
MRISreadTetherFile(MRI_SURFACE *mris, char *fname, float radius)
{
  int    l ;
  float  cx, cy, cz ;
  FILE   *fp ;
  char   line[STRLEN], *cp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, 
                               "MRISreadTetherFile(%s): could not open file",
                               fname)) ;

  /* first count # of labels in file */
  mris->nlabels = 0 ;
  while ((cp = fgetl(line, 199, fp)) != NULL)
    mris->nlabels++ ;
  mris->labels = 
    (MRIS_AREA_LABEL *)calloc(mris->nlabels, sizeof(MRIS_AREA_LABEL)) ;
  if (!mris->labels)
    ErrorExit(ERROR_NOMEMORY, 
              "MRISreadTetherFile: could not allocate %d labels",
              mris->nlabels) ;

  for (l = 0 ; l < mris->nlabels ; l++)
  {
    cp = fgetl(line, 199, fp) ;
    if (!sscanf(cp, "%s %f %f %f", mris->labels[l].name,&cx, &cy, &cz) != 4)
    {
      fclose(fp) ;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE,
                  "MRISreadTetherFile(%s): could not scan parms from %s",
                                  fname, line)) ;
    }
    mris->labels[l].cx = cx; mris->labels[l].cy = cy; mris->labels[l].cz = cz;
    mrisLabelVertices(mris, cx, cy, cz, l, radius) ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
static int
mrisLabelVertices(MRI_SURFACE *mris, float cx, float cy, float cz, int label, 
                  float radius)
{
  VERTEX *v ;
  float  xd, yd, zd, d ;
  int    vno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    xd = cx - v->x ; yd = cy - v->y ; zd = cz - v->z ;
    d = sqrt(xd*xd + yd*yd + zd*zd) ;
#if 0
    if (d <= radius)
      v->label = label ;
#endif
  }
  
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
#define K_A         0.4f
static float kernel[] = { K_A, 0.25f, 0.25f-K_A/2.0f } ;
static int
mrisSmoothCurvatures(MRI_SURFACE *mris, int niterations)
{
  int     vno, i, vn ;
  double  g, H, norm ;
  VERTEX  *vertex ;

  for (i = 0 ; i < niterations ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
        continue ;
      H = kernel[0]*vertex->H ;
      g = kernel[1] ;
      for (vn = 0 ; vn < vertex->vnum ; vn++)
        H += g * mris->vertices[vertex->v[vn]].H ;
      g = kernel[2] ;
      for ( ; vn < vertex->v2num ; vn++)
        H += g * mris->vertices[vertex->v[vn]].H ;
      norm = 
        kernel[0] + 
        vertex->vnum*kernel[1] + 
        (vertex->v2num-vertex->vnum) * kernel[2] ;
      vertex->d = H/norm ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
        continue ;
      vertex->H = vertex->d ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
static int
mrisSmoothNormals(MRI_SURFACE *mris, int niterations)
{
  int     vno, i, vn ;
  double  g ;
  VERTEX  *vertex, *vnb ;
  VECTOR  *v_n, *v_n2 ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_n2 = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */

  for (i = 0 ; i < niterations ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
        continue ;
      VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
      V3_SCALAR_MUL(v_n, kernel[0], v_n) ;
      g = kernel[1] ;
      for (vn = 0 ; vn < vertex->vnum ; vn++)
      {
        vnb = &mris->vertices[vertex->v[vn]] ;
        VECTOR_LOAD(v_n2, vnb->nx, vnb->ny, vnb->nz) ;
        V3_SCALAR_MUL(v_n2, g, v_n2) ;
        V3_ADD(v_n, v_n2, v_n) ;
      }
      g = kernel[2] ;
      for ( ; vn < vertex->v2num ; vn++)
      {
        vnb = &mris->vertices[vertex->v[vn]] ;
        VECTOR_LOAD(v_n2, vnb->nx, vnb->ny, vnb->nz) ;
        V3_SCALAR_MUL(v_n2, g, v_n2) ;
        V3_ADD(v_n, v_n2, v_n) ;
      }
      V3_NORMALIZE(v_n, v_n) ;
      vertex->tdx = V3_X(v_n) ; 
      vertex->tdy = V3_Y(v_n) ;
      vertex->tdz = V3_Z(v_n) ; 
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
        continue ;
      vertex->nx = vertex->tdx ; 
      vertex->ny = vertex->tdy ;
      vertex->nz = vertex->tdz ;
    }
  }

  VectorFree(&v_n) ;
  VectorFree(&v_n2) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisSmoothNormalOutliers(MRI_SURFACE *mris, double ndist)
{
  int     vno, n, m, smooth, nsmoothed ;
  VERTEX  *v, *vn ;
  float   sx, sy, sz, nx, ny, nz, nc, x, y, z, dist, dx, dy, dz ;


  for (nsmoothed = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    x = v->x ;    y = v->y ;   z = v->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    /* see if we are futher than ndist from all our neighbors in
       the normal direction. If so, smooth.
    */
    for (smooth = 1, m = 0 ; smooth && m < v->vnum ; m++)
    {
      vn = &mris->vertices[v->v[m]] ;
      if (!vn->ripflag)
      {
        dx = vn->x - x;
        dy = vn->y - y;
        dz = vn->z - z;
        nc = dx*nx+dy*ny+dz*nz;   /* projection onto normal */
        dist = sqrt(nc) ;         /* distance in normal direction */
        if (dist < ndist)
          smooth = 0 ;
        sx += dx ;
        sy += dy ;
        sz += dz ;
        n++;
      }
    }
    if (!smooth)
      continue ;
    nsmoothed++ ;
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
#if 0
    nc = sx*nx+sy*ny+sz*nz;   /* projection onto normal */
    sx = nc*nx ;              /* move in normal direction */
    sy = nc*ny ;
    sz = nc*nz;
#endif
    
    v->dx += sx ;
    v->dy += sy ;
    v->dz += sz ;
  }
  
  fprintf(stdout, "%d smoothed (%2.2f%%)\n",
          nsmoothed, 100.0f*(float)nsmoothed / (float)mris->nvertices) ;
  return(NO_ERROR) ;
}
static int
mrisComputeAverageNormalTerm(MRI_SURFACE *mris, int navgs, double l_normal)
{
  VERTEX  *v, *vn ;
  double  nc_avg, nc, vnum, delta ;
  int     n, vno, marked ;
  float   x, y, z, dx, dy, dz, nx, ny, nz, s ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    marked = v->marked ;

    /* compute projection onto normal */
    x = v->x ; y = v->y ; z = v->z ; nx = v->nx ; ny = v->ny ; nz = v->nz ;
    for (vnum = nc_avg = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked != marked)
        continue ;
      vnum++ ;
      dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ;
      nc_avg += dx*nx+dy*ny+dz*nz;   /* projection onto normal */
    }
    if (vnum > 0.0)
      nc_avg /= vnum ;
    else
      nc_avg = 0.0 ;
    v->d = nc_avg ;
  }

  mrisAverageDs(mris, navgs) ;

  /* now move each vertex in the direction of the local average */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    marked = v->marked ;

    /* compute projection onto normal */
    nc_avg = v->d ;
    x = v->x ; y = v->y ; z = v->z ; nx = v->nx ; ny = v->ny ; nz = v->nz ;
    for (nc = vnum = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked != marked)
        continue ;
      vnum++ ;
      dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ;
      nc += dx*nx+dy*ny+dz*nz;   /* projection onto normal */
    }
    if (vnum > 0.0)
      nc /= (float)vnum ;
    else
      nc = 0.0 ;

    s = nc_avg < 0.0 ? -1.0 : 1.0 ;
    nc_avg = sqrt(fabs(nc_avg)) * s ;
    s = nc < 0.0 ? -1.0 : 1.0 ;
    nc = sqrt(fabs(nc)) * s ;
    delta = -l_normal * (nc_avg - nc) ;
    dx = nx * delta ; dy = ny * delta ; dz = nz * delta ;
    v->dx += dx ; v->dy += dy ; v->dz += dz ;
  }
  return(NO_ERROR) ;
}
static int
mrisComputeCurvatureTerm(MRI_SURFACE *mris, double l_curv)
{
  int     vno, n ;
  VERTEX  *v, *vn ;
  float   nc, nx, ny, nz, x, y, z, ndx, ndy, ndz, tdx, tdy, tdz, dx, dy, dz, 
          tdist, ndist, nc_avg ;
  double  curv ;

  if (FZERO(l_curv))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    x = v->x ;    y = v->y ;   z = v->z ;

    /* first mean tangential and normal spacing */
    nc_avg = ndx = ndy = ndz = tdist = ndist = 0.0f ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ;
      nc = (dx*nx + dy*ny + dz*nz) ;
      nc_avg += nc ;
#if 0
      ndx += nc*nx ; ndy += nc*ny ; ndz += nc*nz ; 
#endif
      tdx = dx-nc*nx ; tdy = dy-nc*ny ; tdz = dz-nc*nz ; 
      tdist += sqrt(tdx*tdx+tdy*tdy+tdz*tdz) ;
    }
#if 0
    ndx /= (float)v->vnum ; ndy /= (float)v->vnum ; ndz /= (float)v->vnum ;
    ndist = sqrt(ndx*ndx+ndy*ndy+ndz*ndz) ;
    if (nc_avg < 0.0f)
      ndist *= -1 ;       /* neighbors are predominantly below tangent plane */
#else
    ndist = nc_avg ;
#endif
    tdist /= (float)v->vnum ; ndist /= (float)v->vnum ;
    if (FZERO(tdist))
      continue ;
    curv = ndist / tdist ;

    if (fabs(curv) < 0.25f)
      continue ;
    curv *= l_curv ;
    v->dx += curv * nx ;
    v->dy += curv * ny ;
    v->dz += curv * nz ;
    if (vno == Gdiag_no)
      fprintf(stdout,"vertex %d normal curvature term: (%2.3f, %2.3f, %2.3f)\n"
              , vno, curv*nx, curv*ny, curv*nz) ;
  }
  
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
static int
mrisComputeConvexityTerm(MRI_SURFACE *mris, double l_convex)
{
  int     vno, n, m ;
  VERTEX  *vertex, *vn ;
  float   sx, sy, sz, nx, ny, nz, nc, x, y, z ;

  if (FZERO(l_convex))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    nx = vertex->nx ; ny = vertex->ny ; nz = vertex->nz ;
    x = vertex->x ;    y = vertex->y ;   z = vertex->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < vertex->vnum ; m++)
    {
      vn = &mris->vertices[vertex->v[m]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
    nc = sx*nx+sy*ny+sz*nz;   /* projection onto normal */
    if (nc < 0)
      nc = 0 ;
    sx = nc*nx ;              /* move in normal direction */
    sy = nc*ny ;
    sz = nc*nz;
    
    vertex->dx += l_convex * sx ;
    vertex->dy += l_convex * sy ;
    vertex->dz += l_convex * sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d convexity term: (%2.3f, %2.3f, %2.3f)\n",
              vno, vertex->dx, vertex->dy, vertex->dz) ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
static int
mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int     vno, n, m ;
  VERTEX  *vertex, *vn ;
  float   sx, sy, sz, nx, ny, nz, nc, x, y, z ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    nx = vertex->nx ; ny = vertex->ny ; nz = vertex->nz ;
    x = vertex->x ;    y = vertex->y ;   z = vertex->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < vertex->vnum ; m++)
    {
      vn = &mris->vertices[vertex->v[m]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
    nc = sx*nx+sy*ny+sz*nz;   /* projection onto normal */
    sx = l_spring*nc*nx ;              /* move in normal direction */
    sy = l_spring*nc*ny ;
    sz = l_spring*nc*nz ;
    
    vertex->dx += sx ;
    vertex->dy += sy ;
    vertex->dz += sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring normal term:  (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
  }
  

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
static int
mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int     vno, n, m ;
  VERTEX  *v, *vn ;
  float   sx, sy, sz, x, y, z, nc ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    if (v->border && !v->neg)
      continue ;

    x = v->x ;    y = v->y ;   z = v->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < v->vnum ; m++)
    {
      vn = &mris->vertices[v->v[m]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
#if 0
    n = 4 ;  /* avg # of nearest neighbors */
#endif
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
    
    nc = sx*v->nx+sy*v->ny+sz*v->nz;   /* projection onto normal */
    sx -= l_spring*nc*v->nx ;                   /* remove  normal component */
    sy -= l_spring*nc*v->ny ;
    sz -= l_spring*nc*v->nz;

    v->dx += sx ;
    v->dy += sy ;
    v->dz += sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring tangent term: (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
  }
  

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
#define NORMAL_MOVEMENT   0.1
#define NSAMPLES          15
#define SAMPLE_DISTANCE   0.1

static int
mrisComputeIntensityTerm(MRI_SURFACE *mris, double l_intensity, MRI *mri_brain,
                         MRI *mri_smooth, double sigma_global)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z, nx, ny, nz, dx, dy, dz ;
  Real    val0, xw, yw, zw, del, val_outside, val_inside, delI, delV, k,
          ktotal ;
	double  sigma ;

  if (FZERO(l_intensity))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ; y = v->y ; z = v->z ;

    // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_brain, xw, yw, zw, &val0) ;
    sigma = v->val2 ;
		if (FZERO(sigma))
			sigma = sigma_global ;

    nx = v->nx ; ny = v->ny ; nz = v->nz ;

    /* compute intensity gradient using smoothed volume */

#if 1
    {
      Real dist, val, step_size ;
      int  n ;

      step_size = MIN(sigma/2, 0.5) ;
      ktotal = 0.0 ;
      for (n = 0, val_outside = val_inside = 0.0, dist = step_size ; 
           dist <= 2*sigma; 
           dist += step_size, n++)
      {
        k = exp(-dist*dist/(2*sigma*sigma)) ; ktotal += k ;
        xw = x + dist*nx ; yw = y + dist*ny ; zw = z + dist*nz ;
        // MRIworldToVoxel(mri_brain, xw, yw, zw, &xw, &yw, &zw) ;
				MRIsurfaceRASToVoxel(mri_brain, xw, yw, zw, &xw, &yw, &zw);
        MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;
        val_outside += k*val ;
				
        xw = x - dist*nx ; yw = y - dist*ny ; zw = z - dist*nz ;
        // MRIworldToVoxel(mri_brain, xw, yw, zw, &xw, &yw, &zw) ;
				MRIsurfaceRASToVoxel(mri_brain, xw, yw, zw, &xw, &yw, &zw);
        MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;
        val_inside += k*val ;
      }
      val_inside /= (double)ktotal ; val_outside /= (double)ktotal ;
    }
#else
    /* sample outward from surface */
    xw = x + nx ; yw = y + ny ; zw = z + nz ;
    // MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw);
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_outside) ;

    /* sample inward from surface */
    xw = x - nx ; yw = y - ny ; zw = z - nz ;
    // MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw);
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_inside) ;
#endif

#if 1
    /* this stuff was used for the shrink-wrapping */

    /* if everthing is 0 force vertex to move inwards */
    if (val0 < 25 && v->marked)  /* vertex has never been higher than 25 */
    {
      /* just noise - set everything to 0 */
      val0 = val_inside = val_outside = 0 ;
    }
    if (val0 > 25)
    {
      if (vno == 0)
        DiagBreak() ;
      v->marked = 0 ;
    }
    if (val_inside == 0 && val_outside == 0 && val0 == 0)
    {
      val_outside = val0-5 ; val_inside = val0+5 ;
    }
#endif
    delV = v->val - val0 ;
    delI = (val_outside - val_inside) / 2.0 ;
#if 1
    if (!FZERO(delI))
      delI /= fabs(delI) ;
    else
      delI = -1 ;   /* intensities tend to increase inwards */
#if 0
		if (delI > 0)
			delI = 0 ;
#endif
#else
    delI = -1 ;  /* ignore gradient and assume that too bright means move out */
#endif

    if (delV > 5)
      delV = 5 ;
    else if (delV < -5)
      delV = -5 ;
    
    del = l_intensity * delV * delI ;
#if 0
    if (delV > 0)  /* in a sulcus? */
    {
      del *= 10 ; 
      if (Gdiag_no == vno)
	printf("v %d: augmenting intensity term to prevent sulcal crossing\n",vno) ;
    }
#endif

    dx = nx * del ; dy = ny * del ; dz = nz * del ;

    v->dx += dx ;   
    v->dy += dy ;
    v->dz += dz ;

#if 0
    {
      int n ;
      for (n = 0 ; n < v->vnum ; n++)
        if (v->v[n] == Gdiag_no)
        {
          Real xwi, ywi, zwi, xwo, ywo, zwo ;
          
          x = v->x ; y = v->y ; z = v->z ;
          
          /* sample outward from surface */
          xw = x + nx ; yw = y + ny ; zw = z + nz ;
          // MRIworldToVoxel(mri_smooth, xw, yw, zw, &xwo, &ywo, &zwo) ;
	  MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xwo, &ywo, &zwo);
          /* sample inward from surface */
          xw = x - nx ; yw = y - ny ; zw = z - nz ;
          // MRIworldToVoxel(mri_smooth, xw, yw, zw, &xwi, &ywi, &zwi) ;
	  MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xwi, &ywi, &zwi);
          
          // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
	  MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw);
          fprintf(stdout, 
                  "I(%2.1f,%2.1f,%2.1f)=%2.1f, Io(%2.1f,%2.1f,%2.1f)=%2.1f, "
                  "Ii(%2.1f,%2.1f,%2.1f)=%2.1f\n", 
                  xw,yw,zw,val0, xwo, ywo,zwo,val_outside,xwi,ywi,zwi,val_inside);
          fprintf(stdout, "v %d intensity term:      (%2.3f, %2.3f, %2.3f), "
                  "delV=%2.1f, delI=%2.0f\n", vno, dx, dy, dz, delV, delI) ;
        }
    }
#endif

    if (vno == Gdiag_no)
    {
      Real xwi, ywi, zwi, xwo, ywo, zwo ;

      x = v->x ; y = v->y ; z = v->z ;

      /* sample outward from surface */
      xw = x + nx ; yw = y + ny ; zw = z + nz ;
      // MRIworldToVoxel(mri_smooth, xw, yw, zw, &xwo, &ywo, &zwo) ;
      MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xwo, &ywo, &zwo);
      /* sample inward from surface */
      xw = x - nx ; yw = y - ny ; zw = z - nz ;
      // MRIworldToVoxel(mri_smooth, xw, yw, zw, &xwi, &ywi, &zwi) ;
      MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xwi, &ywi, &zwi);

      // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw);
      fprintf(stdout, 
              "I(%2.1f,%2.1f,%2.1f)=%2.1f, Io(%2.1f,%2.1f,%2.1f)=%2.1f, "
              "Ii(%2.1f,%2.1f,%2.1f)=%2.1f\n", 
              xw,yw,zw,val0, xwo, ywo,zwo,val_outside,xwi,ywi,zwi,val_inside);
      fprintf(stdout, "v %d intensity term:      (%2.3f, %2.3f, %2.3f), "
              "delV=%2.1f, delI=%2.0f\n", vno, dx, dy, dz, delV, delI) ;
    }
  }
  

  return(NO_ERROR) ;
}
#if 1
static int
mrisComputeIntensityGradientTerm(MRI_SURFACE*mris, double l_grad,
                                 MRI *mri_brain, MRI *mri_smooth)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z, nx, ny, nz ;
  Real    val0, mag0, xw, yw, zw, del, mag_outside, mag_inside, delI, delV,
          dx, dy, dz, val_outside, val_inside, val_dist, dn, xw1, yw1, zw1 ;

  if (FZERO(l_grad))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    /*
      sample intensity value and derivative in normal direction
      at current point.
    */
    x = v->x+v->nx ; y = v->y+v->ny ; z = v->z+v->nz ;
    //MRIworldToVoxel(mri_brain, x, y, z, &xw1, &yw1, &zw1) ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw1, &yw1, &zw1) ;
    x = v->x ; y = v->y ; z = v->z ;
    //MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    nx = xw1-xw ; ny = yw1-yw ; nz = zw1-zw ; 
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    MRIsampleVolume(mri_brain, xw, yw, zw, &val0) ;
    mag0 = sqrt(dx*dx+dy*dy+dz*dz) ;
    MRIsampleVolumeDerivative(mri_smooth, xw, yw, zw, nx, ny, nz, &dn) ;

    /* compute intensity gradient using smoothed volume */

    /* sample outward from surface */
    xw = x + nx ; yw = y + ny ; zw = z + nz ;
    //MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    mag_outside = sqrt(dx*dx+dy*dy+dz*dz) ;
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_outside) ;

    /* sample inward from surface */
    xw = x - nx ; yw = y - ny ; zw = z - nz ;
    //MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    mag_inside = sqrt(dx*dx+dy*dy+dz*dz) ;
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_inside) ;

    if (mag_outside > mag_inside)   /* gradients suggest edge is outwards */
      val_dist = val_outside - v->val ;
    else   /* gradients suggest edge is inwards */
      val_dist = val_inside - v->val ;

#if 0
    /* diminish the effect of gradients that are in locations whose
       intensity values are far from the target.
       */
    val_dist = 1 / (1 + val_dist*val_dist) ;
#else
    /* check to make sure that gradient is changing in right direction */
    val_dist = 1 ;
    /* is this right??? Used to be val0 > v->val, what about dn < 0?? */
    if (val0 > v->val)   /* current point is brighter than target */
    {
      /* dn > 0 implies should move in, but gradient mag points out */
      if (((mag_inside < mag_outside) && dn > 0) ||
          ((mag_inside > mag_outside) && dn < 0))
        val_dist = 0 ;
    }
    else                /* current point is dimmer than target */
    {
      /* dn > 0 implies should move out, but gradient mag points in */
      if (((mag_inside > mag_outside) && dn > 0) ||
          ((mag_inside < mag_outside) && dn < 0))
        val_dist = 0 ;
    }
#endif
    
    delV = 1.0f /*v->mean - mag0*/ ;
    delI = (mag_outside - mag_inside) / 2.0 ;
#if 1
    if (!FZERO(delI))
      delI /= fabs(delI) ;
#endif
    del = val_dist * l_grad * delV * delI ;
    dx = v->nx * del ; dy = v->ny * del ; dz = v->nz * del ;

    v->dx += dx ;   
    v->dy += dy ;
    v->dz += dz ;
    if (vno == Gdiag_no)
      fprintf(stdout, 
        "v %d intensity gradient term: (%2.3f, %2.3f, %2.3f) "
              "(mag = %2.1f, [%2.1f,%2.1f])\n",
              vno, v->dx, v->dy, v->dz, mag0, mag_inside, mag_outside) ;
  }
  

  return(NO_ERROR) ;
}
#else
static int
mrisComputeIntensityGradientTerm(MRI_SURFACE*mris, double l_grad,
                                 MRI *mri_brain, MRI *mri_smooth)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z, nx, ny, nz ;
  Real    xw, yw, zw, dx, dy, dz, val, mag, mag_next, scale ;

  if (FZERO(l_grad))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ; y = v->y ; z = v->z ;
    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    //MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    mag = sqrt(dx*dx + dy*dy + dz*dz) ;
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val) ;

    /* sample outward from surface */
    xw = x + v->dx ; yw = y + v->dy ; zw = z + v->dz ;
    //MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    mag_next = sqrt(dx*dx + dy*dy + dz*dz) ;

    /* gradient decreasing in outwards direction */
    scale = 1.0+l_grad ;
    if (fabs(mag) > fabs(mag_next))
      scale /= 1.0f ;

    v->dx *= scale ;   
    v->dy *= scale ;
    v->dz *= scale ;
    if (vno == Gdiag_no)
      fprintf(stdout, 
        "v %d intensity gradient term: (%2.3f, %2.3f, %2.3f) "
              "(mag = %2.1f, mag_next=%2.1f)\n",
              vno, v->dx, v->dy, v->dz, mag, mag_next) ;
  }
  

  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Apply a uniform outward expansion force.
------------------------------------------------------*/
static int
mrisComputeSphereTerm(MRI_SURFACE *mris, double l_sphere, float radius)
{
  int     vno ;
  VERTEX  *v ;
  float   r, x, y, z, x0, y0, z0 ;

  if (FZERO(l_sphere))
    return(0.0f) ;

#if 0
#if 0
  radius = MRISaverageRadius(mris) ;
#else
  radius = MRISmaxRadius(mris) ;
#endif
  fprintf(stdout, "max radius = %2.1f\n", radius) ;
#endif
  x0 = (mris->xlo+mris->xhi)/2.0f ; 
  y0 = (mris->ylo+mris->yhi)/2.0f ; 
  z0 = (mris->zlo+mris->zhi)/2.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x-x0 ; y = v->y-y0 ; z = v->z-z0 ;
    r = sqrt(x*x+y*y+z*z) ;
    x /= r ; y /= r ; z /= r ;  /* normal direction */
    r = (radius - r) / radius ;
    v->dx += r*l_sphere * x ;
    v->dy += r*l_sphere * y ;
    v->dz += r*l_sphere * z ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d sphere    term: (%2.3f, %2.3f, %2.3f), r=%2.2f\n",
              vno, v->dx, v->dy, v->dz, r) ;
  }
  

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Apply a uniform outward expansion force.
------------------------------------------------------*/
static int
mrisComputeExpansionTerm(MRI_SURFACE *mris, double l_expand)
{
  int     vno ;
  VERTEX  *v ;

  if (FZERO(l_expand))
    return(0.0f) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    v->dx += l_expand * v->nx ;
    v->dy += l_expand * v->ny ;
    v->dz += l_expand * v->nz ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        File Format is:

        name x y z
------------------------------------------------------*/
static int
mrisComputeSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int     vno, n, m ;
  VERTEX  *v, *vn ;
  float   sx, sy, sz, x, y, z, dist_scale ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

#if METRIC_SCALE
  if (mris->patch)
    dist_scale = 1.0 ;
  else
    dist_scale = sqrt(mris->orig_area / mris->total_area) ;
#else
  dist_scale = 1.0 ;
#endif
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    if (v->border && !v->neg)
      continue ;

    x = v->x ;    y = v->y ;   z = v->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < v->vnum ; m++)
    {
      vn = &mris->vertices[v->v[m]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
#if 0
    n = 4 ;  /* avg # of nearest neighbors */
#endif
    if (n>0)
    {
      sx = dist_scale*sx/n;
      sy = dist_scale*sy/n;
      sz = dist_scale*sz/n;
    }
    
    sx *= l_spring ; sy *= l_spring ; sz *= l_spring ; 
    v->dx += sx ;
    v->dy += sy ;
    v->dz += sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring term:         (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
  }
  

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute a spring term, and normalize it by removing the
          average normal component which typically forces the surface
          to shrink.
------------------------------------------------------*/
#if 0
static int
mrisComputeNormalizedSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int     vno, n, m ;
  VERTEX  *v, *vn ;
  float   sx, sy, sz, x, y, z, dist_scale, nx, ny, nz, dx, dy, dz ;
  double  dot_total, avg_len, std_len, len, num ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

#if METRIC_SCALE
  if (mris->patch)
    dist_scale = 1.0 ;
  else
    dist_scale = sqrt(mris->orig_area / mris->total_area) ;
#else
  dist_scale = 1.0 ;
#endif

  std_len = avg_len = dot_total = 0.0 ;
  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ;      y = v->y ;     z = v->z ;
    nx = v->nx ;    ny = v->ny ;   nz = v->nz ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < v->vnum ; m++)
    {
      vn = &mris->vertices[v->v[m]] ;
      if (!vn->ripflag)
      {
        dx = vn->x - x;
        dy = vn->y - y;
        dz = vn->z - z;
        sx += dx ;
        sy += dy ;
        sz += dz ;
        n++;
        len = sqrt(dx*dx + dy*dy + dz*dz) ;
        avg_len += len ;
        std_len += len*len ;
      }
    }
    num += n ;
    if (n>0)
    {
      sx = dist_scale*sx/n;
      sy = dist_scale*sy/n;
      sz = dist_scale*sz/n;
    }
  }
  avg_len /= num ;
  std_len = sqrt(std_len/num-avg_len*avg_len) ;

  num = (double)mrisValidVertices(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ;      y = v->y ;     z = v->z ;
    nx = v->nx ;    ny = v->ny ;   nz = v->nz ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < v->vnum ; m++)
    {
      vn = &mris->vertices[v->v[m]] ;
      if (!vn->ripflag)
      {
        dx = vn->x - x;
        dy = vn->y - y;
        dz = vn->z - z;
        len = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (vno == 91951 && v->v[m] == 90994)
        {
          fprintf(stdout, "v %d-->%d: len %2.3f, avg %2.3f (ratio=%2.3f)\n",
                  vno, v->v[m], len, avg_len, len/avg_len) ;
        }
#if 1
        dx /= len ; dy /= len ; dz /= len ;
        if (len/avg_len > 1.5)
          len = avg_len * 1.0 * dist_scale ;
        dx *= len ; dy *= len ; dz *= len ;
#endif
        sx += dx ;
        sy += dy ;
        sz += dz ;
        n++;
      }
    }
    if (n>0)
    {
      sx = dist_scale*sx/n;
      sy = dist_scale*sy/n;
      sz = dist_scale*sz/n;
    }

    dot_total += l_spring*(nx*sx + ny*sy + nz*sz) ;
    v->dx += l_spring * sx ;
    v->dy += l_spring * sy ;
    v->dz += l_spring * sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring nm term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
  }
  dot_total /= num ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    nx = v->nx ;    ny = v->ny ;   nz = v->nz ;

    v->dx -= dot_total * nx ;
    v->dy -= dot_total * ny ;
    v->dz -= dot_total * nz ;
  }
  

  return(NO_ERROR) ;
}
#else
static int
mrisComputeNormalizedSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int     vno, n, m ;
  VERTEX  *v, *vn ;
  float   sx, sy, sz, x, y, z, dist_scale, nx, ny, nz, dx, dy, dz ;
  double  dot_total, num ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

#if METRIC_SCALE
  if (mris->patch)
    dist_scale = 1.0 ;
  else
    dist_scale = sqrt(mris->orig_area / mris->total_area) ;
#else
  dist_scale = 1.0 ;
#endif

  dot_total = 0.0 ;
  num = (double)mrisValidVertices(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ;      y = v->y ;     z = v->z ;
    nx = v->nx ;    ny = v->ny ;   nz = v->nz ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < v->vnum ; m++)
    {
      vn = &mris->vertices[v->v[m]] ;
      if (!vn->ripflag)
      {
        dx = vn->x - x;
        dy = vn->y - y;
        dz = vn->z - z;
        sx += dx ;
        sy += dy ;
        sz += dz ;
        n++;
      }
    }
    if (n>0)
    {
      sx = dist_scale*sx/n;
      sy = dist_scale*sy/n;
      sz = dist_scale*sz/n;
    }

    dot_total += l_spring*(nx*sx + ny*sy + nz*sz) ;
    v->dx += l_spring * sx ;
    v->dy += l_spring * sy ;
    v->dz += l_spring * sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring norm term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
  }
  dot_total /= num ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    nx = v->nx ;    ny = v->ny ;   nz = v->nz ;

    v->dx -= dot_total * nx ;
    v->dy -= dot_total * ny ;
    v->dz -= dot_total * nz ;
  }
  
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the folding of the surface.
------------------------------------------------------*/
double
MRIScomputeFolding(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *vertex ;
  double k1, k2, folding, area ;

  folding = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    k1 = (double)vertex->k1 ; k2 = (double)vertex->k2 ;
    area = (double)vertex->area ;
    folding += area * (k1 - k2) * (k1 - k2) ;
  }

  return(folding) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the effects of the gradient of the distance term
------------------------------------------------------*/
static int
mrisComputeDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_y, *v_delta, *v_n ;
  float   l_dist, d0, dt, delta, nc, scale, norm, max_del ;
  VERTEX  *v, *vn ;
  int     vno, n, vnum, max_v, max_n ;

  if (!FZERO(parms->l_nldist))
    mrisComputeNonlinearDistanceTerm(mris, parms) ;

  l_dist = parms->l_dist ;
  if (FZERO(l_dist))
    return(NO_ERROR) ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_y = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;
  norm = 1.0f / mris->avg_nbrs ;

#if METRIC_SCALE
  if (mris->patch)
    scale = 1.0f ;
  else
    if (mris->status == MRIS_PARAMETERIZED_SPHERE)
      scale = sqrt(mris->orig_area / mris->total_area) ;
  else
    scale = mris->neg_area < mris->total_area ? 
      sqrt(mris->orig_area / (mris->total_area-mris->neg_area)) :
      sqrt(mris->orig_area / mris->total_area) ;
#else
  scale = 1.0f ;
#endif
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "distance scale = %2.3f\n", scale) ;
max_del = 10000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vnum = v->vtotal ;
    if (v->ripflag || vnum <= 0)
      continue ;

    if (v->border)
      DiagBreak() ;
    V3_CLEAR(v_delta) ;
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;

    if (vno == Gdiag_no && Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, 
            "computing distance term for v %d @ (%2.2f, %2.2f, %2.2f)\n",
              vno, v->x, v->y, v->z) ;

    for (n = 0 ; n < vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag)
        continue ;
      d0 = v->dist_orig[n] ; dt = scale * v->dist[n] ; delta = dt - d0 ;
if (fabs(delta) > max_del)
{
  max_del = delta ;
  max_v = vno ;
  max_n = n ;
}
      VECTOR_LOAD(v_y, vn->x - v->x, vn->y - v->y, vn->z - v->z) ;
      if ((V3_LEN_IS_ZERO(v_y)))
        continue ;
      V3_NORMALIZE(v_y, v_y) ;   /* make it a unit vector */
      V3_SCALAR_MUL(v_y, delta, v_y) ;
      V3_ADD(v_y, v_delta, v_delta) ;
      if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
        fprintf(stdout, 
                "nbr %d (%6.6d) @ (%2.2f, %2.2f, %2.2f), "
                "d0 %2.2f, dt %2.2f, delta %2.3f\n\ty=%2.3f, %2.3f, %2.3f)\n",
                n, v->v[n], vn->x, vn->y, vn->z, d0, dt,
                delta, V3_X(v_y), V3_Y(v_y), V3_Z(v_y)) ;
    }

    V3_SCALAR_MUL(v_delta, norm, v_delta) ;

    if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stdout, 
                "total delta=(%2.3f, %2.3f, %2.3f)\n",
                V3_X(v_delta), V3_Y(v_delta), V3_Z(v_delta)) ;
    /* take out normal component */
    nc = V3_DOT(v_n, v_delta) ; V3_SCALAR_MUL(v_n, -nc, v_n) ;
    V3_ADD(v_delta, v_n, v_delta) ;

    v->dx += l_dist * V3_X(v_delta) ;
    v->dy += l_dist * V3_Y(v_delta) ;
    v->dz += l_dist * V3_Z(v_delta) ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d, distance term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
  }

  VectorFree(&v_n) ;
  VectorFree(&v_y) ;
  VectorFree(&v_delta) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the effects of the gradient of the distance term
------------------------------------------------------*/
static int
mrisComputeNonlinearDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_y, *v_delta, *v_n ;
  float   l_dist, d0, dt, delta, nc, scale, norm, ratio ;
  VERTEX  *v, *vn ;
  int     vno, n, vnum ;

  l_dist = parms->l_nldist ;
  if (FZERO(l_dist))
    return(NO_ERROR) ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_y = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;
  norm = 1.0f / mris->avg_nbrs ;

#if METRIC_SCALE
  if (mris->patch)
    scale = 1.0f ;
  else
    if (mris->status == MRIS_PARAMETERIZED_SPHERE)
      scale = sqrt(mris->orig_area / mris->total_area) ;
  else
    scale = mris->neg_area < mris->total_area ? 
      sqrt(mris->orig_area / (mris->total_area-mris->neg_area)) :
      sqrt(mris->orig_area / mris->total_area) ;
#else
  scale = 1.0f ;
#endif
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "distance scale = %2.3f\n", scale) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vnum = v->vtotal ;
    if (v->ripflag || vnum <= 0)
      continue ;

    if (v->border)
      DiagBreak() ;
    V3_CLEAR(v_delta) ;
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;

    if (vno == Gdiag_no && Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, 
            "computing distance term for v %d @ (%2.2f, %2.2f, %2.2f)\n",
              vno, v->x, v->y, v->z) ;

    for (n = 0 ; n < vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag)
        continue ;
      d0 = v->dist_orig[n] ; dt = scale * v->dist[n] ; delta = dt - d0 ;
      VECTOR_LOAD(v_y, vn->x - v->x, vn->y - v->y, vn->z - v->z) ;
      if ((V3_LEN_IS_ZERO(v_y)))
        continue ;
      V3_NORMALIZE(v_y, v_y) ;   /* make it a unit vector */
      if (!FZERO(d0))
      {
        ratio = dt / d0 ;
        delta *= 1 / (1 + exp(-1*ratio)) ;
      }
      V3_SCALAR_MUL(v_y, delta, v_y) ;
      V3_ADD(v_y, v_delta, v_delta) ;
      if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
        fprintf(stdout, 
                "nbr %d (%6.6d) @ (%2.2f, %2.2f, %2.2f), "
                "d0 %2.2f, dt %2.2f, delta %2.3f\n\ty=%2.3f, %2.3f, %2.3f)\n",
                n, v->v[n], vn->x, vn->y, vn->z, d0, dt,
                delta, V3_X(v_y), V3_Y(v_y), V3_Z(v_y)) ;
    }

    V3_SCALAR_MUL(v_delta, norm, v_delta) ;

    if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stdout, 
                "total delta=(%2.3f, %2.3f, %2.3f)\n",
                V3_X(v_delta), V3_Y(v_delta), V3_Z(v_delta)) ;
    /* take out normal component */
    nc = V3_DOT(v_n, v_delta) ; V3_SCALAR_MUL(v_n, -nc, v_n) ;
    V3_ADD(v_delta, v_n, v_delta) ;

    v->dx += l_dist * V3_X(v_delta) ;
    v->dy += l_dist * V3_Y(v_delta) ;
    v->dz += l_dist * V3_Z(v_delta) ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d, distance term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
  }

  VectorFree(&v_n) ;
  VectorFree(&v_y) ;
  VectorFree(&v_delta) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisLogIntegrationParms(FILE *fp, MRI_SURFACE *mris,INTEGRATION_PARMS *parms)
{
  char  *cp, host_name[STRLEN] ;

  if (!fp)
    return(NO_ERROR) ;

  cp = getenv("HOST") ;
  if (cp)
    strcpy(host_name, cp) ;
  else
    strcpy(host_name, "unknown") ;

  fprintf(fp, "tol=%2.1e, host=%5.5s, nav=%d, nbrs=%d",
          (float)parms->tol, host_name, parms->n_averages, mris->nsize) ;
  if (!FZERO(parms->l_area))
    fprintf(fp, ", l_area=%2.3f", parms->l_area) ;
  if (!FZERO(parms->l_shrinkwrap))
    fprintf(fp, ", l_shrinkwrap=%2.3f", parms->l_shrinkwrap) ;
  if (!FZERO(parms->l_external))
    fprintf(fp, ", l_extern=%2.3f", parms->l_external) ;
  if (!FZERO(parms->l_parea))
    fprintf(fp, ", l_parea=%2.3f", parms->l_parea) ;
  if (!FZERO(parms->l_nlarea))
    fprintf(fp, ", l_nlarea=%2.3f", parms->l_nlarea) ;
  if (!FZERO(parms->l_nldist))
    fprintf(fp, ", l_nldist=%2.3f", parms->l_nldist) ;
  if (!FZERO(parms->l_angle))
    fprintf(fp, ", l_angle=%2.3f", parms->l_angle) ;
  if (!FZERO(parms->l_repulse))
    fprintf(fp, ", l_repulse=%2.3f", parms->l_repulse) ;
  if (!FZERO(parms->l_repulse_ratio))
    fprintf(fp, ", l_repulse_ratio=%2.3f", parms->l_repulse_ratio) ;
  if (!FZERO(parms->l_surf_repulse))
    fprintf(fp, ", l_surf_repulse=%2.3f", parms->l_surf_repulse) ;
  if (!FZERO(parms->l_corr))
    fprintf(fp, ", l_corr=%2.3f", parms->l_corr) ;
  if (!FZERO(parms->l_spring))
    fprintf(fp, ", l_spring=%2.3f", parms->l_spring) ;
  if (!FZERO(parms->l_spring_norm))
    fprintf(fp, ", l_spring_norm=%2.3f", parms->l_spring_norm) ;
  if (!FZERO(parms->l_tspring))
    fprintf(fp, ", l_tspring=%2.3f", parms->l_tspring) ;
  if (!FZERO(parms->l_nspring))
    fprintf(fp, ", l_nspring=%2.3f", parms->l_nspring) ;
  if (!FZERO(parms->l_dist))
    fprintf(fp, ", l_dist=%2.3f", parms->l_dist) ;
  if (!FZERO(parms->l_intensity))
    fprintf(fp, ", l_intensity=%2.3f", parms->l_intensity) ;
  if (!FZERO(parms->l_grad))
    fprintf(fp, ", l_grad=%2.3f", parms->l_grad) ;
  if (!FZERO(parms->l_sphere))
    fprintf(fp, ", l_sphere=%2.3f", parms->l_sphere) ;
  if (!FZERO(parms->l_expand))
    fprintf(fp, ", l_expand=%2.3f", parms->l_expand) ;
  if (!FZERO(parms->l_curv))
    fprintf(fp, ", l_curv=%2.3f", parms->l_curv) ;
  if (!FZERO(parms->l_convex))
    fprintf(fp, ", l_convex=%2.3f", parms->l_convex) ;
  if (!FZERO(parms->l_boundary))
    fprintf(fp, ", l_boundary=%2.3f", parms->l_boundary) ;
  if (!FZERO(parms->l_neg))
    fprintf(fp, ", l_neg=%2.3f", parms->l_neg) ;
  if (!FZERO(parms->l_tsmooth))
    fprintf(fp, ", l_tsmooth=%2.3f", parms->l_tsmooth) ;
  fprintf(fp, "\n") ;
  switch (parms->integration_type)
  {
  case INTEGRATE_LM_SEARCH:
    fprintf(fp, "using binary search line minimization\n") ;
    break ;
  case INTEGRATE_LINE_MINIMIZE:
    fprintf(fp, "using quadratic fit line minimization\n") ;
    break ;
  case INTEGRATE_ADAPTIVE:
    fprintf(fp, 
            "mom=%2.2f, dt=%2.2f, base_dt=%2.3f, dt_inc=%2.2f, "
            "dt_dec=%2.2f, err_rat=%2.2f\n", 
            (float)parms->momentum, (float)parms->dt,
            (float)parms->base_dt, (float)parms->dt_increase, 
            (float)parms->dt_decrease, (float)parms->error_ratio) ;
    break ;
  default:
  case INTEGRATE_MOMENTUM:
    fprintf(fp, 
            "mom=%2.2f, dt=%2.2f\n",(float)parms->momentum, (float)parms->dt);
    break ;
  }
#if 0
  fprintf(fp, "nbhd_size=%d, max_nbrs=%d ", parms->nbhd_size,parms->max_nbrs);
#endif
  if (parms->desired_rms_height > 0.0)
    fprintf(fp, "desired rms height=%2.3f", parms->desired_rms_height) ;
  fprintf(fp, "\n") ;
  fflush(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisWriteSnapshots(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int t)
{
  char base_name[STRLEN] ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;

  strcpy(base_name, parms->base_name) ;

  MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
  sprintf(parms->base_name, "%s_pial", base_name) ;
  mrisWriteSnapshot(mris, parms, t) ;

  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  sprintf(parms->base_name, "%s_white", base_name) ;
  mrisWriteSnapshot(mris, parms, t) ;

  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;

  strcpy(parms->base_name, base_name) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisWriteSnapshot(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int t)
{
  char fname[STRLEN], path[STRLEN], base_name[STRLEN], *cp ;

  FileNamePath(mris->fname, path) ;
  sprintf(base_name, "%s/%s.%s", path, 
          mris->hemisphere == LEFT_HEMISPHERE ? "lh":"rh", parms->base_name);
  if ((cp = strstr(base_name, ".geo")) != NULL)
  {
    *cp = 0;
    sprintf(fname, "%s%4.4d.geo", base_name, t) ;
    *cp = '.' ;
  }
  else
    sprintf(fname, "%s%4.4d", base_name, t) ;
#if 1
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "writing %s...", fname) ;
  if (mris->patch)
    MRISwritePatch(mris, fname) ;
  else
    MRISwrite(mris, fname) ;
  
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "done.\n") ;
    fflush(stderr) ;
  }
#endif

  if (mris->status == MRIS_PARAMETERIZED_SPHERE && DIAG_VERBOSE_ON)
  {
    MRI_SP *mrisp = (MRI_SP *)mris->vp ;

    sprintf(fname, "%s%4.4d.hipl", parms->base_name, t) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "writing %s\n", fname) ;
    MRIStoParameterization(mris, mrisp, 1, 0) ;
    MRISPwrite(mrisp, fname) ;
  }
  if (!FZERO(parms->l_area) && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "%s%4.4d.area_error", base_name, t) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, " %s...", fname) ;
    MRISwriteAreaError(mris, fname) ;
  }

  if (!FZERO(parms->l_corr) && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "%s%4.4d.angle_error", base_name, t) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, " %s...", fname) ;
    MRISwriteAngleError(mris, fname) ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeDistanceError(MRI_SURFACE *mris)
{
  VERTEX  *v ;
  int     vno, n, nvertices, max_v, max_n ;
  double  dist_scale, sse_dist, delta, v_sse, max_del ;

#if METRIC_SCALE
  if (mris->patch)
    dist_scale = 1.0 ;
  else
    if (mris->status == MRIS_PARAMETERIZED_SPHERE)
      dist_scale = sqrt(mris->orig_area / mris->total_area) ;
  else
    dist_scale = mris->neg_area < mris->total_area ? 
      sqrt(mris->orig_area / (mris->total_area-mris->neg_area)) :
      sqrt(mris->orig_area / mris->total_area) ;
#else
  dist_scale = 1.0 ;
#endif
max_del = -1.0 ; max_v = max_n = -1 ;
  for (sse_dist = 0.0, nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
#if NO_NEG_DISTANCE_TERM
    if (v->neg)
      continue ;
#endif
    nvertices++ ;
    for (v_sse = 0.0, n = 0 ; n < v->vtotal ; n++)
    {
#if NO_NEG_DISTANCE_TERM
    if (mris->vertices[v->v[n]].neg)
      continue ;
#endif
		if (v->dist_orig[n] >= UNFOUND_DIST)
			continue ;
		delta = dist_scale*v->dist[n] - v->dist_orig[n] ;
			
		if (fabs(delta) > fabs(max_del))
		{
			max_del = delta ;
			max_v = vno ;
			max_n = n ;
		}
		v_sse += delta*delta ;
		if (!finite(delta) || !finite(v_sse))
			DiagBreak() ;
    }
		if (v_sse > 10000)
			DiagBreak() ;
    sse_dist += v_sse ;
    if (!finite(sse_dist) || !finite(v_sse))
      DiagBreak() ;
  }
	
  /*fprintf(stdout, "max_del = %f at v %d, n %d\n", max_del, max_v, max_n) ;*/
  return(sse_dist) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeNonlinearDistanceSSE(MRI_SURFACE *mris)
{
  VERTEX  *v ;
  int     vno, n, nvertices, max_v, max_n ;
  double  dist_scale, sse_dist, delta, v_sse, max_del, ratio ;

#if METRIC_SCALE
  if (mris->patch)
    dist_scale = 1.0 ;
  else
    if (mris->status == MRIS_PARAMETERIZED_SPHERE)
      dist_scale = sqrt(mris->orig_area / mris->total_area) ;
  else
    dist_scale = mris->neg_area < mris->total_area ? 
      sqrt(mris->orig_area / (mris->total_area-mris->neg_area)) :
      sqrt(mris->orig_area / mris->total_area) ;
#else
  dist_scale = 1.0 ;
#endif
max_del = -1.0 ; max_v = max_n = -1 ;
  for (sse_dist = 0.0, nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nvertices++ ;
    for (v_sse = 0.0, n = 0 ; n < v->vtotal ; n++)
    {
#if 0
      delta = dist_scale*v->dist[n] - v->dist_orig[n] ;
      delta *= delta ;
#endif
      if (FZERO(v->dist_orig[n]))
        continue ;
      ratio = dist_scale*v->dist[n] / v->dist_orig[n] ;
      delta = log(1+exp(ratio)) ;
      v_sse += delta ;
      if (!finite(delta) || !finite(v_sse))
        DiagBreak() ;
    }
    sse_dist += v_sse ;
    if (!finite(sse_dist) || !finite(v_sse))
      DiagBreak() ;
  }

  return(sse_dist) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeThicknessSmoothnessEnergy(MRI_SURFACE *mris, double l_tsmooth)
{
  int     vno, n ;
  double  sse_tsmooth, v_sse, dn, dx, dy, dz, x, y, z, d0 ;
  VERTEX  *v, *vn ;

  if (FZERO(l_tsmooth))
    return(0.0) ;

  for (sse_tsmooth = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    x = v->x ; y = v->y ; z = v->z ; 
    d0 = SQR(x-v->origx)+SQR(y-v->origy)+SQR(z-v->origz) ;
    for (v_sse = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (!vn->ripflag)
      {
        dx = vn->x-vn->origx ; dy = vn->y-vn->origy ; dz = vn->z-vn->origz ; 
        dn = (dx*dx+dy*dy+dz*dz) ;
        v_sse += (dn - d0) * (dn - d0) ;
      }
    }
    sse_tsmooth += v_sse ;
  }
  return(l_tsmooth * sse_tsmooth) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static double
mrisComputeRepulsiveEnergy(MRI_SURFACE *mris, double l_repulse)
{
  int     vno, n ;
  double  sse_repulse, v_sse, dist, dx, dy, dz, x, y, z ;
  VERTEX  *v, *vn ;

  if (FZERO(l_repulse))
    return(0.0) ;

  for (sse_repulse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    x = v->x ; y = v->y ; z = v->z ; 
    for (v_sse = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (!vn->ripflag)
      {
        dx = x - vn->x ; dy = y - vn->y ; dz = z - vn->z ; 
        dist = sqrt(dx*dx+dy*dy+dz*dz) + REPULSE_E ;
        v_sse += REPULSE_K / (dist*dist*dist*dist) ;
      }
    }
    sse_repulse += v_sse ;
  }
  return(l_repulse * sse_repulse) ;
}
#else
static double
mrisComputeRepulsiveEnergy(MRI_SURFACE *mris, double l_repulse, MHT *mht)
{
  int     vno, num, min_vno, i, n ;
  float   dist, dx, dy, dz, x, y, z, min_d ;
  double  sse_repulse, v_sse ;
  VERTEX  *v, *vn ;
  MHBT    *bucket ;
  MHB     *bin ;

  if (FZERO(l_repulse))
    return(NO_ERROR) ;

  min_d = 1000.0 ; min_vno = 0 ;
  for (sse_repulse = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    x = v->x ; y = v->y ; z = v->z ; 
    bucket = MHTgetBucket(mht, x, y, z) ;
    if (!bucket)
      continue ;
    for (v_sse = 0.0, bin = bucket->bins, num = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      if (bin->fno == vno)
        continue ;  /* don't be repelled by myself */
      for (n = 0 ; n < v->vtotal ; n++)
        if (v->v[n] == bin->fno)
          break ;
      if (n < v->vtotal)   /* don't be repelled by a neighbor */
        continue ;
      vn = &mris->vertices[bin->fno] ;
      if (!vn->ripflag)
      {
        dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ; 
        dist = sqrt(dx*dx+dy*dy+dz*dz) + REPULSE_E ;
        if (vno == Gdiag_no)
        {
          if (dist-REPULSE_E < min_d)
          {
            min_vno = bin->fno ;
            min_d = dist-REPULSE_E ;
          }
        }
        dist = dist*dist*dist ; dist *= dist ; /* dist^6 */
        v_sse += REPULSE_K / dist ;
      }
    }
    sse_repulse += v_sse ;

    if (vno == Gdiag_no && !FZERO(v_sse))
    {
      printf("v %d: repulse sse:    min_dist=%2.4f, v_sse %2.4f\n", vno, 
              min_d, v_sse) ;
    }
  }
  return(sse_repulse) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/

static double
mrisComputeRepulsiveRatioEnergy(MRI_SURFACE *mris, double l_repulse)
{
  int     vno, n ;
  double  sse_repulse, v_sse, dist, dx, dy, dz, x, y, z, canon_dist,
          cdx, cdy, cdz ;
  VERTEX  *v, *vn ;

  if (FZERO(l_repulse))
    return(0.0) ;

  for (sse_repulse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    x = v->x ; y = v->y ; z = v->z ; 
    for (v_sse = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (!vn->ripflag)
      {
        dx = x - vn->x ; dy = y - vn->y ; dz = z - vn->z ; 
        dist = sqrt(dx*dx+dy*dy+dz*dz) ;
        cdx = vn->cx - v->cx ; cdy = vn->cy - v->cy ; cdz = vn->cz - v->cz ; 
        canon_dist = sqrt(cdx*cdx+cdy*cdy+cdz*cdz) + REPULSE_E ;
        dist /= canon_dist ;
        dist += REPULSE_E ;
#if 0
        v_sse += REPULSE_K / (dist*dist*dist*dist) ;
#else
        v_sse += REPULSE_K / (dist*dist) ;
#endif
      }
    }
    sse_repulse += v_sse ;
  }
  return(l_repulse * sse_repulse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeThicknessSmoothnessTerm(MRI_SURFACE *mris, double l_tsmooth)
{
  int     vno, n, num ;
  float   dx, dy, dz, x, y, z, dn, d0, vx, vy, vz, delta ;
  VERTEX  *v, *vn ;

  if (FZERO(l_tsmooth))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    x = v->x ; y = v->y ; z = v->z ; 
    vx = v->x - v->origx ; vy = v->y - v->origy ; vz = v->z - v->origz ;
    d0 = vx*vx + vy*vy + vz*vz ;
    dx = dy = dz = 0.0 ;
    for (num = n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (!vn->ripflag)
      {
        dn = SQR(vn->x-vn->origx)+SQR(vn->origy-vn->y)+SQR(vn->origz-vn->z) ;
        delta = d0 - dn ;
        dx -= delta * vx ; dy -= delta * vy ; dz -= delta * vz ; 
        num++ ;
      }
    }
    if (num)
    {
      dx /= num ; dy /= num ; dz /= num ;
    }
    v->dx += dx ; v->dy += dy ; v->dz += dz ;
    if (vno == Gdiag_no)
    {
      fprintf(stdout, "v %d tsmooth term:        (%2.3f, %2.3f, %2.3f)\n",
              vno, dx, dy, dz) ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeRepulsiveTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht)
{
  int     vno, num, min_vno, i, n ;
  float   dist, dx, dy, dz, x, y, z, sx, sy, sz, min_d, min_scale, norm ;
  double  scale ;
  VERTEX  *v, *vn ;
  MHBT    *bucket ;
  MHB     *bin ;

  if (FZERO(l_repulse))
    return(NO_ERROR) ;

  min_d = 100000.0 ; min_scale = 1.0 ; min_vno = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
			DiagBreak() ;
    x = v->x ; y = v->y ; z = v->z ; 
    bucket = MHTgetBucket(mht, x, y, z) ;
    if (!bucket)
      continue ;
    sx = sy = sz = 0.0 ;
    for (bin = bucket->bins, num = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      if (bin->fno == vno)
        continue ;  /* don't be repelled by myself */
      for (n = 0 ; n < v->vtotal ; n++)
        if (v->v[n] == bin->fno)
          break ;
      if (n < v->vtotal)   /* don't be repelled by a neighbor */
        continue ;
      vn = &mris->vertices[bin->fno] ;
      if (!vn->ripflag)
      {
        dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ; 
        dist = sqrt(dx*dx+dy*dy+dz*dz) + REPULSE_E ;
        scale = -4*REPULSE_K / (dist*dist*dist*dist*dist*dist*dist) ;  /* ^-7 */
        if (vno == Gdiag_no)
        {
          if (dist-REPULSE_E < min_d)
          {
            min_vno = bin->fno ;
            min_d = dist-REPULSE_E ;
            min_scale = scale ;
          }
        }
        norm = sqrt(dx*dx+dy*dy+dz*dz) ; 
				if (FZERO(norm))
					norm = 1.0 ;
        dx /= norm ; dy /= norm ; dz /= norm ;
				if (!finite(dx) || !finite(dy) || !finite(dz))
					DiagBreak() ;
        sx += scale * dx ; sy += scale * dy ; sz += scale*dz ;
        num++ ;
      }
    }
    if (num)
    {
      scale = l_repulse / (double)num ;
      sx *= scale ; sy *= scale ; sz *= scale ;
    }
    v->dx += sx ; v->dy += sy ; v->dz += sz ;
    if ((vno == Gdiag_no) && min_d < 1000)
    {
      vn = &mris->vertices[min_vno] ;
      dx = x - vn->x ; dy = y - vn->y ; dz = z - vn->z ; 
      
      fprintf(stdout, "v %d self repulse term:   (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
      fprintf(stdout, "min_dist @ %d = %2.2f, scale = %2.1f\n",
              min_vno, min_d, min_scale) ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeRepulsiveRatioTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht)
{
  int     vno, num, min_vno, i ;
  float   dist, dx, dy, dz, x, y, z, sx, sy, sz, min_d, min_scale, canon_dist,
          cdx, cdy, cdz ;
  double  scale ;
  VERTEX  *v, *vn ;
  MHBT    *bucket ;
  MHB     *bin ;

  if (FZERO(l_repulse))
    return(NO_ERROR) ;

  min_d = 1000.0 ; min_scale = 1.0 ; min_vno = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    x = v->x ; y = v->y ; z = v->z ; 
    bucket = MHTgetBucket(mht, x, y, z) ;
    if (!bucket)
      continue ;
    sx = sy = sz = 0.0 ;
    for (bin = bucket->bins, num = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      if (bin->fno == vno)
        continue ;  /* don't be repelled by myself */
      vn = &mris->vertices[bin->fno] ;
      if (!vn->ripflag)
      {
        dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ; 
        dist = sqrt(dx*dx+dy*dy+dz*dz) ;
        cdx = vn->cx - v->cx ; cdy = vn->cy - v->cy ; cdz = vn->cz - v->cz ; 
        canon_dist = sqrt(cdx*cdx+cdy*cdy+cdz*cdz) + REPULSE_E ;
        dist /= canon_dist ;
        dist += REPULSE_E ;
#if 0
        scale = -4*REPULSE_K / (dist*dist*dist*dist*dist) ;
#else
        scale = -4*REPULSE_K / (dist*dist*dist) ;
#endif
        if (vno == Gdiag_no)
        {
          if (dist-REPULSE_E < min_d)
          {
            min_vno = bin->fno ;
            min_d = dist-REPULSE_E ;
            min_scale = scale ;
          }
        }
        sx += scale * dx ; sy += scale * dy ; sz += scale*dz ;
        num++ ;
      }
    }
    if (num)
    {
      scale = l_repulse / (double)num ;
      sx *= scale ; sy *= scale ; sz *= scale ;
    }
    v->dx += sx ; v->dy += sy ; v->dz += sz ;
    if (vno == Gdiag_no)
    {
      vn = &mris->vertices[min_vno] ;
      dx = x - vn->x ; dy = y - vn->y ; dz = z - vn->z ; 
      
      fprintf(stdout, "v %d repulse term:        (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
      fprintf(stdout, "min_dist @ %d = %2.2f, scale = %2.1f\n",
              min_vno, min_d, min_scale) ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeSpringEnergy(MRI_SURFACE *mris)
{
  int     vno, n ;
  double  area_scale, sse_spring, v_sse ;
  VERTEX  *v ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  for (sse_spring = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    for (v_sse = 0.0, n = 0 ; n < v->vnum ; n++)
      v_sse += (v->dist[n]*v->dist[n]) ;
    sse_spring += area_scale * v_sse ;
  }
  return(sse_spring) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeTangentialSpringEnergy(MRI_SURFACE *mris)
{
  int     vno, n ;
  double  area_scale, sse_spring, v_sse ;
  VERTEX  *v, *vn ;
  float   dx, dy, dz, x, y, z, nc, dist_sq ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  for (sse_spring = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    x = v->x ;    y = v->y ;   z = v->z ;

    for (v_sse = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ; 
      nc = dx * v->nx + dy*v->ny + dz*v->nz ;
      dx -= nc*v->nx ; dy -= nc*v->ny ; dz -= nc*v->nz ; 
      dist_sq = dx*dx+dy*dy+dz*dz ;
      v_sse += dist_sq ;
    }
    sse_spring += area_scale * v_sse ;
  }
  return(sse_spring) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        each face has 2 triangles defined by it:

       V0       d      V3
        o--------------o
        |              |
        | A0           |
      a |              | c
        |              |
        |           A1 |
        o--------------o
       V1      b        V2        

       a = V1 - V0
       d = V3 - V0
       A0 = 0.5 (a x d) . n

       b = V1 - V2
       c = V3 - V2
       A1 = 0.5 (c x b) . n

        each face has 1 triangle defined by it:

       V0    b     V2
        o----------o
        |         /      
        | A0    /        
      a |     /          
        |   /            
        | /
        o
       V1

       a = V1 - V0
       b = V2 - V0
       A0 = 0.5 (a x b) . n

------------------------------------------------------*/
static int
mrisComputeAngleAreaTerms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     fno, ano ;
  VERTEX  *v0, *v1, *v2, *va, *vb, *vo ;
  VECTOR  *v_a, *v_b, *v_a_x_n, *v_b_x_n, *v_n, *v_tmp, *v_sum ;
  FACE    *face ;
  float   orig_area, area, l_parea, l_area, l_angle, delta, len, area_scale ;

#if METRIC_SCALE
  if (mris->patch || 
      (mris->status != MRIS_SPHERE && 
       mris->status != MRIS_PARAMETERIZED_SPHERE))
    area_scale = 1.0f ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0f ;
#endif

  l_angle = parms->l_angle ;
  l_area = parms->l_area ;
  l_parea = parms->l_parea ;
  if (!FZERO(parms->l_nlarea))
    mrisComputeNonlinearAreaTerm(mris, parms) ;

  if (FZERO(l_area) && FZERO(l_angle) && FZERO(l_parea))
    return(NO_ERROR) ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;

  v_tmp = VectorAlloc(3, MATRIX_REAL) ;   
  v_sum = VectorAlloc(3, MATRIX_REAL) ;   
  v_a_x_n = VectorAlloc(3, MATRIX_REAL) ;      
  v_b_x_n = VectorAlloc(3, MATRIX_REAL) ;      

  /* calculcate movement of each vertex caused by each triangle */
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    VECTOR_LOAD(v_n, face->nx, face->ny, face->nz) ;
    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;
    VERTEX_EDGE(v_a, v0, v1) ;  
    VERTEX_EDGE(v_b, v0, v2) ;
    orig_area = face->orig_area ; area = area_scale * face->area ;
    delta = 0.0 ;
    if (!FZERO(l_parea))
      delta += l_parea * (area - orig_area) ; 
    if (!FZERO(l_area))
    {
      if (area <= 0.0f)
        delta += l_area * (area - orig_area) ; 
    }

    V3_CROSS_PRODUCT(v_a, v_n, v_a_x_n) ;
    V3_CROSS_PRODUCT(v_b, v_n, v_b_x_n) ;

    /* calculate movement of vertices in order, 0-3 */
    
    /* v0 */
    V3_SCALAR_MUL(v_a_x_n, -1.0f, v_sum) ;
    V3_ADD(v_sum, v_b_x_n, v_sum) ;
    V3_SCALAR_MUL(v_sum, delta, v_sum) ;
    v0->dx += V3_X(v_sum) ; v0->dy += V3_Y(v_sum) ; v0->dz += V3_Z(v_sum) ;
    
    /* v1 */
    V3_SCALAR_MUL(v_b_x_n, -delta, v_sum) ;
    v1->dx += V3_X(v_sum) ; v1->dy += V3_Y(v_sum) ; v1->dz += V3_Z(v_sum) ;
    
    /* v2 */
    V3_SCALAR_MUL(v_a_x_n, delta, v_sum) ;
    v2->dx += V3_X(v_sum) ; v2->dy += V3_Y(v_sum) ; v2->dz += V3_Z(v_sum) ;
    
    /* now calculate the angle contributions */
    if (!FZERO(l_angle)) 
    {
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      {
        switch (ano)
        {
        default:
        case 0: vo = v0 ; va = v2 ; vb = v1 ; break ;
        case 1: vo = v1 ; va = v0 ; vb = v2 ; break ;
        case 2: vo = v2 ; va = v1 ; vb = v0 ; break ;
        }
        delta = deltaAngle(face->angle[ano],face->orig_angle[ano]);
#if ONLY_NEG_AREA_TERM
        if (face->angle[ano] >= 0.0f)
          delta = 0.0f ;
#endif
        delta *= parms->l_angle ;
        VERTEX_EDGE(v_a, vo, va) ; VERTEX_EDGE(v_b, vo, vb) ;
        
        /* this angle's contribution to va */
        V3_CROSS_PRODUCT(v_a, v_n, v_tmp) ;
        len = V3_DOT(v_a,v_a) ;
        if (!FZERO(len))
          V3_SCALAR_MUL(v_tmp, delta/len, v_tmp) ;
        else
          V3_SCALAR_MUL(v_tmp, 0.0f, v_tmp) ;
        va->dx += V3_X(v_tmp) ; va->dy += V3_Y(v_tmp); va->dz += V3_Z(v_tmp);
        
        /* this angle's contribution to vb */
        V3_CROSS_PRODUCT(v_n, v_b, v_sum) ;
        len = V3_DOT(v_b,v_b) ;
        if (!FZERO(len))
          V3_SCALAR_MUL(v_sum, delta/len, v_sum) ;
        else
          V3_SCALAR_MUL(v_sum, 0.0f, v_sum) ;
        vb->dx += V3_X(v_sum) ; vb->dy += V3_Y(v_sum); vb->dz += V3_Z(v_sum);
        
        /* this angle's contribution to vo */
        V3_ADD(v_tmp, v_sum, v_sum) ;
        vo->dx -= V3_X(v_sum); vo->dy -= V3_Y(v_sum); vo->dz -= V3_Z(v_sum) ;
      }
    }
  }    /* done with all faces */
  
  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_tmp) ;
  VectorFree(&v_sum) ;
  VectorFree(&v_n) ;

  VectorFree(&v_a_x_n) ;
  VectorFree(&v_b_x_n) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeNonlinearAreaTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     fno ;
  VERTEX  *v0, *v1, *v2 ;
  VECTOR  *v_a, *v_b, *v_a_x_n, *v_b_x_n, *v_n, *v_tmp, *v_sum ;
  FACE    *face ;
  double  orig_area, area, delta, area_scale, scale, l_nlarea, ratio ;

#if METRIC_SCALE
  if (mris->patch || 
      (mris->status != MRIS_SPHERE && 
       mris->status != MRIS_PARAMETERIZED_SPHERE))
    area_scale = 1.0f ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0f ;
#endif

  l_nlarea = parms->l_nlarea ;

  if (FZERO(l_nlarea))
    return(NO_ERROR) ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;

  v_tmp = VectorAlloc(3, MATRIX_REAL) ;   
  v_sum = VectorAlloc(3, MATRIX_REAL) ;   
  v_a_x_n = VectorAlloc(3, MATRIX_REAL) ;      
  v_b_x_n = VectorAlloc(3, MATRIX_REAL) ;      

  /* calculcate movement of each vertex caused by each triangle */
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    VECTOR_LOAD(v_n, face->nx, face->ny, face->nz) ;
    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;
    VERTEX_EDGE(v_a, v0, v1) ;  
    VERTEX_EDGE(v_b, v0, v2) ;
    orig_area = face->orig_area ; area = area_scale * face->area ;
#if SCALE_NONLINEAR_AREA
    if (!FZERO(orig_area))
      ratio = area / orig_area ;
    else
      ratio = 0.0f ;
#else
    ratio = area ;
#endif

    if (ratio > MAX_NEG_RATIO)
      ratio = MAX_NEG_RATIO ;
    else if (ratio < -MAX_NEG_RATIO)
      ratio = -MAX_NEG_RATIO ;
#if 0
    scale = l_nlarea * (1 - (1/(1.0+exp(-NEG_AREA_K*ratio)))) ;
#else
    scale = l_nlarea / (1.0+exp(NEG_AREA_K*ratio)) ;
#endif
    delta = scale * (area - orig_area) ; 

    V3_CROSS_PRODUCT(v_a, v_n, v_a_x_n) ;
    V3_CROSS_PRODUCT(v_b, v_n, v_b_x_n) ;

    /* calculate movement of vertices in order, 0-3 */
    
    /* v0 */
    V3_SCALAR_MUL(v_a_x_n, -1.0f, v_sum) ;
    V3_ADD(v_sum, v_b_x_n, v_sum) ;
    V3_SCALAR_MUL(v_sum, delta, v_sum) ;
    v0->dx += V3_X(v_sum) ; v0->dy += V3_Y(v_sum) ; v0->dz += V3_Z(v_sum) ;
    
    /* v1 */
    V3_SCALAR_MUL(v_b_x_n, -delta, v_sum) ;
    v1->dx += V3_X(v_sum) ; v1->dy += V3_Y(v_sum) ; v1->dz += V3_Z(v_sum) ;
    
    /* v2 */
    V3_SCALAR_MUL(v_a_x_n, delta, v_sum) ;
    v2->dx += V3_X(v_sum) ; v2->dy += V3_Y(v_sum) ; v2->dz += V3_Z(v_sum) ;
  }    /* done with all faces */
  
  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_tmp) ;
  VectorFree(&v_sum) ;
  VectorFree(&v_n) ;

  VectorFree(&v_a_x_n) ;
  VectorFree(&v_b_x_n) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeSurfaceRepulsionTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht)
{
  int     vno, max_vno, i ;
  float   dx, dy, dz, x, y, z, sx, sy, sz,norm[3],dot;
  float   max_scale, max_dot ;
  double  scale ;
  VERTEX  *v, *vn ;
  MHBT    *bucket ;
  MHB     *bin ;

  if (FZERO(l_repulse))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    x = v->x ; y = v->y ; z = v->z ; 
    bucket = MHTgetBucket(mht, x, y, z) ;
    if (!bucket)
      continue ;
    bin = bucket->bins ;
    sx = sy = sz = 0.0 ;
    max_dot = max_scale = 0.0 ; max_vno = 0 ;
    for (i = 0 ; i < bucket->nused ; i++, bin++)
    {
      vn = &mris->vertices[bin->fno] ;
      if (bin->fno == Gdiag_no)
        DiagBreak() ;
      if (vn->ripflag)
        continue ;
      dx = x - vn->origx ; dy = y - vn->origy ; dz = z - vn->origz ; 
      mrisComputeOrigNormal(mris, bin->fno, norm) ;
      dot = dx*norm[0] + dy*norm[1] + dz*norm[2] ;
			if (dot > 1)
				continue ;
      if (dot < 0 && vno == Gdiag_no)
        DiagBreak() ;
      if (dot > MAX_NEG_RATIO)
        dot = MAX_NEG_RATIO ;
      else if (dot < -MAX_NEG_RATIO)
        dot = -MAX_NEG_RATIO ;
#if 0
      scale = l_repulse / (1.0+exp(NEG_AREA_K*dot)) ;
#else
			scale = l_repulse*pow(1.0-(double)dot,4.0) ;
#endif
      if (scale > max_scale)
      {
        max_scale = scale ;
        max_vno = bin->fno ;
        max_dot = dot ;
      }
      sx += (scale*v->nx) ; sy += (scale*v->ny) ; sz += (scale*v->nz) ;
    }
    
    v->dx += sx ; v->dy += sy ; v->dz += sz ;
    if (vno == Gdiag_no)
    {
      vn = &mris->vertices[max_vno] ;
      dx = x - vn->x ; dy = y - vn->y ; dz = z - vn->z ; 

      fprintf(stdout, "v %d inside repulse term:  (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
      fprintf(stdout, "max_scale @ %d = %2.2f, max dot = %2.2f\n", 
              max_vno, max_scale, max_dot) ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static double
mrisComputeAverageHeight(MRI_SURFACE *mris)
{
  int    vno, n, nv ;
  VERTEX *vertex, *vn ;
  double height ;
  float  dx, dy, dz, nx, ny, nz, x, y, z ;

  for (height = 0.0, nv = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    nx = vertex->nx ; ny = vertex->ny ; nz = vertex->nz ;
    if (vertex->ripflag)
      continue ;
    nv += vertex->vtotal ;
    for (n = 0 ; n < vertex->vtotal ; n++)
    {
      vn = &mris->vertices[vertex->v[n]] ;
      dx = x - vn->x ; dy = y - vn->y ; dz = z - vn->z ;
      height += fabs(dx * nx + dy * ny + dz * nz) ;
    }
  }

  return(sqrt(height/(double)nv)) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisProjectSurface(MRI_SURFACE *mris)
{

  /*  MRISupdateSurface(mris) ;*/
  switch (mris->status)
  {
  case MRIS_PARAMETERIZED_SPHERE:
    MRISprojectOntoSphere(mris, mris, mris->radius) ;
    break ;
  case MRIS_SPHERE:
    MRISprojectOntoSphere(mris, mris, mris->radius) ;
    break ;
  case MRIS_ELLIPSOID:
    MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
    break ;
  case PROJECT_PLANE:
    /*    mrisOrientPlane(mris) ;*/
    break ;
  case MRIS_RIGID_BODY:
    /*    MRISprojectOntoSphere(mris, mris, mris->radius) ;*/
    mris->status = MRIS_RIGID_BODY ;
    break ;
  default:
    break ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisOrientSurface(MRI_SURFACE *mris)
{
  switch (mris->status)
  {
  case MRIS_RIGID_BODY:
  case MRIS_PARAMETERIZED_SPHERE:
  case MRIS_SPHERE:
  case MRIS_ELLIPSOID:
    MRISupdateEllipsoidSurface(mris) ;
    break ;
  case MRIS_PLANE:
    mrisOrientPlane(mris) ;
    break ;
  default:
/*    MRISupdateSurface(mris) ;*/
    break ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int   
mrisLogStatus(MRI_SURFACE *mris,INTEGRATION_PARMS *parms,FILE *fp, float dt)
{
  float  area_rms, angle_rms, curv_rms, sse, dist_rms, corr_rms ;
  int    negative ;

  if (!(Gdiag & DIAG_SHOW))
    return(NO_ERROR) ;

  negative = MRIScountNegativeTriangles(mris) ;
  sse  = mrisComputeError(mris, parms,&area_rms,&angle_rms,&curv_rms,&dist_rms,
                   &corr_rms);
#if  0
  sse = MRIScomputeSSE(mris, parms) ;
#endif
#if 0
  sse /= (float)mrisValidVertices(mris) ;
  sse = sqrt(sse) ;
#endif
  if (FZERO(parms->l_corr) && FZERO(parms->l_pcorr))
    fprintf(fp, "%3.3d: dt: %2.2f, sse: %2.1f (%2.3f, %2.1f, %2.3f), "
            "neg: %d (%%%2.3f:%%%2.2f), avgs: %d\n", 
            parms->t, dt, sse, area_rms, (float)DEGREES(angle_rms), dist_rms,
            negative, 100.0*mris->neg_area/(mris->neg_area+mris->total_area),
            100.0*mris->neg_orig_area/(mris->orig_area),
            parms->n_averages);
  else
    fprintf(fp, "%3.3d: dt: %2.3f, sse: %2.1f (%2.3f, %2.1f, %2.3f, %2.3f), "
            "neg: %d (%%%2.2f:%%%2.2f), avgs: %d\n", 
            parms->t, dt, sse, area_rms, (float)DEGREES(angle_rms), dist_rms,
            corr_rms, negative, 
            100.0*mris->neg_area/(mris->neg_area+mris->total_area),
            100.0*mris->neg_orig_area/(mris->orig_area),
            parms->n_averages);

  fflush(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadVertexPositions(MRI_SURFACE *mris, char *name)
{
  char    fname[STRLEN] ;
  int     vno, nvertices, nfaces, magic, version, tmp, ix, iy, iz, n, type ;
  VERTEX  *vertex ;
  FILE    *fp ;

  type = MRISfileNameType(name) ;
  MRISbuildFileName(mris, name, fname) ;
  if (type == MRIS_GEO_TRIANGLE_FILE)
    return(mrisReadGeoFilePositions(mris, fname)) ;
  else if (type == MRIS_ICO_FILE)
    return(ICOreadVertexPositions(mris, fname, CURRENT_VERTICES)) ;
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadVertexPosition(%s): could not open file %s", 
                 name, fname));

  fread3(&magic, fp) ;
  if (magic == QUAD_FILE_MAGIC_NUMBER) 
  {
    version = -1;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "new surface file format\n");
  }
  else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER)
  {
    version = -2 ;
  }
  else if (magic == TRIANGLE_FILE_MAGIC_NUMBER)
  {
    fclose(fp) ;
    if (mrisReadTriangleFilePositions(mris,fname)  != NO_ERROR)
      ErrorReturn(Gerror, (Gerror, "mrisReadTriangleFile failed.\n")) ;
    version = -3 ;
  }
  else 
  {
    rewind(fp);
    version = 0;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      printf("surfer: old surface file format\n");
  }

  if (version >= -2)
  {
    fread3(&nvertices, fp);
    fread3(&nfaces, fp); nfaces *= 2 ;
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "reading %d vertices and %d faces.\n",nvertices,nfaces);
    
    if (nvertices != mris->nvertices || nfaces != mris->nfaces)
    {
      fclose(fp) ;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, 
                "MRISreadVertexPositions(%s): surfaces differ.\n",
                fname)) ;
    }

    for (vno = 0 ; vno < nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (version == -1)
      {
        fread2(&ix,fp);
        fread2(&iy,fp);
        fread2(&iz,fp);
        vertex->x = ix/100.0;
        vertex->y = iy/100.0;
        vertex->z = iz/100.0;
      }
      else  /* version == -2 */
      {
        vertex->x = freadFloat(fp) ;
        vertex->y = freadFloat(fp) ;
        vertex->z = freadFloat(fp) ;
      }
      if (version == 0)  /* old surface format */
      {
        fread1(&tmp, fp);
        for (n=0;n<vertex->num;n++)
          fread3(&tmp,fp);
      } 
    }
    
    fclose(fp);
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScomputeMetricProperties(MRI_SURFACE *mris)
{
  MRIScomputeNormals(mris) ;
  mrisComputeVertexDistances(mris) ;
  mrisComputeSurfaceDimensions(mris) ;
  MRIScomputeTriangleProperties(mris) ;  /* compute areas and normals */
  mrisOrientSurface(mris) ;
  if (mris->status == MRIS_PARAMETERIZED_SPHERE || 
      mris->status == MRIS_RIGID_BODY)
  {
    double old_area ;
    old_area = mris->total_area ;
    mris->total_area = M_PI * mris->radius * mris->radius * 4.0 ;
#if 0
    fprintf(stdout, "computed %.1f, analytic %.1f\n",
            old_area, mris->total_area) ;
    mris->total_area = old_area ;
#endif
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisCountTotalNeighbors(MRI_SURFACE *mris)
{
  int     vno, total ;
  VERTEX  *v ;

  for (total = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    
    total += v->vtotal+1 ;  /* include this vertex in count */
  }
  return(total) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisComputeBoundaryTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
#if 0
  int    vno ;
  VERTEX *v ;
  double l_boundary ;

  l_boundary = parms->l_boundary ;

  if ((mris->status != MRIS_PLANE) || FZERO(l_boundary))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || !v->border)
      continue ;
    if (v->neg)
    {
#if 0
      v->dx -= l_boundary * v->bnx ;
      v->dy -= l_boundary * v->bny ;
#endif
    }
    else
    {
#if 1
      v->dx += l_boundary * v->bnx ;
      v->dy += l_boundary * v->bny ;
#endif
    }
  }

#endif
  return(NO_ERROR) ;
}
#endif
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeNegTerm(MRI_SURFACE *mris,INTEGRATION_PARMS *parms)
{
  int    vno, n, neg ;
  VERTEX *v, *vn ;
  double l_neg, dx, dy, dz, len ;

  l_neg = parms->l_neg ;

  if (FZERO(l_neg))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->neg)
      neg = 1 ;
    else
      neg = 0 ;
    dx = dy = dz = 0.0 ;
    for (n = 0 ; n < v->vtotal ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->neg)
      { neg++ ; continue ; }
      dx += vn->x - v->x ; dy += vn->y - v->y ; dz += vn->z - v->z ;
    }
    len = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (!FZERO(len) && neg > 0)
    {
      dx /= len ; dy /= len ; dz /= len ;
      v->dx += l_neg*dx ; v->dy += l_neg*dy ; v->dz += l_neg*dz ;
    }
  }

  return(NO_ERROR) ;
}
#endif
static int
mrisFlipPatch(MRI_SURFACE *mris)
{
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->x = -v->x ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           read in vertex positions from the original file,
           and compute the metric properties of that surface.
           Note we must change the status of the mris as the
           original surface has no externally defined orientation.
------------------------------------------------------*/
int
MRISreadOriginalProperties(MRI_SURFACE *mris, char *sname)
{
  int  old_status ;
   
  if (!sname)
    sname = "smoothwm" ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;

  old_status = mris->status ;
  mris->status = MRIS_PATCH ;  /* so no orientating will be done */
  if (MRISreadVertexPositions(mris, sname) != NO_ERROR)
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE, 
                 "MRISreadOriginalProperties: could not read surface file %s",
                 sname)) ;
  
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeTriangleProperties(mris) ;
  MRISstoreMetricProperties(mris) ;
  mris->status = old_status ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeTriangleProperties(mris) ;
  mrisOrientSurface(mris) ;
  mris->orig_area = mris->total_area ;
        
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsaveVertexPositions(MRI_SURFACE *mris, int which)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
#if 0
    if (v->ripflag)
      continue ;
#endif
    switch (which)
    {
		case WHITE_VERTICES:
      v->whitex = v->x ; v->whitey = v->y ; v->whitez = v->z ;
			break ;
    case PIAL_VERTICES:
      v->pialx = v->x ; v->pialy = v->y ; v->pialz = v->z ;
      break ;
    case INFLATED_VERTICES:
      v->infx = v->x ; v->infy = v->y ; v->infz = v->z ;
      break ;
    case FLATTENED_VERTICES:
      v->fx = v->x ; v->fy = v->y ; v->fz = v->z ;
      break ;
    case CANONICAL_VERTICES:
      v->cx = v->x ; v->cy = v->y ; v->cz = v->z ;
      break ;
    case ORIGINAL_VERTICES:
      v->origx = v->x ; v->origy = v->y ; v->origz = v->z ;
      break ;
    default:
    case TMP_VERTICES:
      v->tx = v->x ; v->ty = v->y ; v->tz = v->z ;
      break ;
    }
  }
  if (which == CANONICAL_VERTICES)
    MRIScomputeCanonicalCoordinates(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScomputeCanonicalCoordinates(MRI_SURFACE *mris)
{
  float   theta, phi, r, d, x, y, z ;
  VERTEX  *v ;
  int     vno ;

  r = mris->radius = MRISaverageCanonicalRadius(mris) ;
  r = mris->radius = (float)nint(mris->radius) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    x = v->cx ; y = v->cy ; z = v->cz ;
    theta = atan2(y/r, x/r) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = r*r-z*z ; if (d < 0.0) d = 0.0 ;
    phi = atan2(sqrt(d), z) ;
    v->theta = theta ; v->phi = phi ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISrestoreVertexPositions(MRI_SURFACE *mris, int which)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
#if 0
    if (v->ripflag)
      continue ;
#endif
    switch (which)
    {
    case WHITE_VERTICES:
      v->x = v->whitex ; v->y = v->whitey ; v->z = v->whitez ;
      break ;
    case PIAL_VERTICES:
      v->x = v->pialx ; v->y = v->pialy ; v->z = v->pialz ;
      break ;
    case INFLATED_VERTICES:
      v->x = v->infx ; v->y = v->infy ; v->z = v->infz ;
      break ;
    case FLATTENED_VERTICES:
      v->x = v->fx ; v->y = v->fy ; v->z = v->fz ;
      break ;
    case CANONICAL_VERTICES:
      v->x = v->cx ; v->y = v->cy ; v->z = v->cz ;
      break ;
    case ORIGINAL_VERTICES:
      v->x = v->origx ; v->y = v->origy ; v->z = v->origz ;
      break ;
    default:
    case TMP_VERTICES:
      v->x = v->tx ; v->y = v->ty ; v->z = v->tz ;
      break ;
    }
  }
  mrisComputeSurfaceDimensions(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRISpercentDistanceError(MRI_SURFACE *mris)
{
  VERTEX  *v ;
  int     vno, n, nvertices ;
  double  dist_scale, pct, dist, odist, mean, mean_error ;

  if (mris->patch)
    dist_scale = 1.0 ;
  else
    dist_scale = sqrt(mris->orig_area / mris->total_area) ;

  mean = 0.0 ;
  for (pct = 0.0, nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    for (n = 0 ; n < v->vtotal ; n++)
    {
      nvertices++ ;
      dist = dist_scale*v->dist[n] ;
      odist = v->dist_orig[n] ;
      if (!FZERO(odist))
        pct += fabs(dist - odist) / odist ;
      mean += odist ;
    }
  }

  for (mean_error = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    for (n = 0 ; n < v->vtotal ; n++)

    {
      dist = dist_scale*v->dist[n] ;
      odist = v->dist_orig[n] ;
      mean_error += fabs(dist-odist) ;
    }
  }
  mean /= (double)nvertices ;
  mean_error /= (double)nvertices ;
#if 0
  pct /= (double)nvertices ;
#else
  if (!FZERO(mean))
    pct = mean_error / mean ;
  else
    pct = 1000.0 ;  /* should never happen */
#endif
  return(100.0*pct) ;
}
/*-----------------------------------------------------
  int MRISfileNameType()

        Parameters: char *fname

        Returns value: int 

        Description: return filetype using the extension
	             default is MRIS_BINARY_QUADRANGLE_FILE
------------------------------------------------------*/
int   
MRISfileNameType(char *fname)
{
  int   type ;
  char  *dot, ext[STRLEN], str[STRLEN] ;

	FileNameOnly(fname, str) ; fname = str ; /* remove path */
  ext[0] = 0 ;
  dot = strrchr(fname, '@') ;   /* forces the type of the file */
  if (dot)
  {
    *dot = 0 ;   /* remove it from file name so that fopen will work */
    strcpy(ext, dot+1) ;
  }
  else     /* no explicit type - check extension */
  {
    dot = strrchr(fname, '.') ;
    if (dot)
      strcpy(ext, dot+1) ;
  }
  StrUpper(ext) ;
  if (!strcmp(ext, "ASC"))
    type = MRIS_ASCII_TRIANGLE_FILE ;
  else if (!strcmp(ext, "GEO"))
    type = MRIS_GEO_TRIANGLE_FILE ;
  else if (!strcmp(ext, "TRI") || !strcmp(ext, "ICO"))
    type = MRIS_ICO_FILE ;
  else if (!strcmp(ext, "VTK"))
    type = MRIS_VTK_FILE ;
  else
    type = MRIS_BINARY_QUADRANGLE_FILE ;

  return(type) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteAscii(MRI_SURFACE *mris, char *fname)
{
  int     vno, fno, n ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE, "MRISwriteAscii: could not open file %s",fname));
                 
  fprintf(fp, "#!ascii version of %s\n", mris->fname) ;
  fprintf(fp, "%d %d\n", mris->nvertices, mris->nfaces) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    fprintf(fp, "%f  %f  %f  %d\n", v->x, v->y, v->z, v->ripflag) ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      fprintf(fp, "%d ", face->v[n]) ;
    fprintf(fp, "%d\n", face->ripflag) ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteVTK(MRI_SURFACE *mris, char *fname)
{
  int     vno, fno, n ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE, "MRISwriteAscii: could not open file %s",fname));
                 
  fprintf(fp, "# vtk DataFile Version 1.0\n") ;
  fprintf(fp, "vtk output\nASCII\nDATASET POLYDATA\nPOINTS %d float\n",
          mris->nvertices) ;
  /*  fprintf(fp, "%d %d\n", mris->nvertices, mris->nfaces) ;*/

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    fprintf(fp, "%f  %f  %f\n", v->x, v->y, v->z) ;
  }
  fprintf(fp, "POLYGONS %d %d\n", mris->nfaces, mris->nfaces*4) ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    fprintf(fp,"%d ", VERTICES_PER_FACE) ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      fprintf(fp, "%d ", face->v[n]) ;
    fprintf(fp, "\n") ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteGeo(MRI_SURFACE *mris, char *fname)
{
  int     vno, fno, n, actual_vno, toggle, nfaces, nvertices, vnos[300000] ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  nfaces = mrisValidFaces(mris) ;
  nvertices = mrisValidVertices(mris) ;
  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE, "MRISwriteGeo: could not open file %s",fname));

  fprintf(fp, "    1 %d %d %d\n", nvertices, 
          nfaces, VERTICES_PER_FACE*nfaces) ;
  fprintf(fp, "    1 %d\n", nfaces) ;
  for (actual_vno = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    vnos[vno] = actual_vno++ ;
  }

  for (toggle = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    fprintf(fp, "%12.5e %12.5e %12.5e", v->x, v->y, v->z) ;
    if (toggle)
    {
      toggle = 0 ;
      fprintf(fp, "\n") ;
    }
    else
    {
      toggle = 1 ;
      fprintf(fp, " ") ;
    }
  }
  if (toggle)
    fprintf(fp, "\n") ;
  
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    for (n = 0 ; n < VERTICES_PER_FACE-2 ; n++)
      fprintf(fp, "%d ", vnos[face->v[n]]+1); /* 1-based */

    /* swap order on output to conform to movie.byu convention */
    fprintf(fp, "%d ",vnos[face->v[VERTICES_PER_FACE-1]]+1) ;
    fprintf(fp, "-%d\n",vnos[face->v[VERTICES_PER_FACE-2]]+1);
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters: MRI_SURFACE *mris, char *fname

        Returns value: int

        Description: write ascii icosahedron data (vertices and face vertices info)
------------------------------------------------------*/
/* note that .tri or .ico file.  numbering is 1-based output.*/
int
MRISwriteICO(MRI_SURFACE *mris, char *fname)
{
  int     vno, fno, nfaces, nvertices;
  int     actual_fno, actual_vno;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;
  // get the valid faces and vertices numbers
  nfaces = mrisValidFaces(mris) ;
  nvertices = mrisValidVertices(mris) ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE, "MRISwriteICO: could not open file %s",fname));
  
  // write number of vertices
  fprintf(fp, "%8d\n", nvertices);
  // go over all vertices and ignoring bad ones
  actual_vno = 1; // count from 1 (1-based output)
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    fprintf(fp, "%8d %8.4f %8.4f %8.4f\n", actual_vno, v->x, v->y, v->z) ;
    actual_vno++;
  }
  // write number of faces
  fprintf(fp, "%8d\n", nfaces);
  // go over all faces and ignoring bad ones
  actual_fno = 1; // count from 1 (1-based)
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    // make the vertex number 1 based
    // the vertex ordering flipped to clockwise (see icosahedron.c)
    fprintf(fp, "%8d %8d %8d %8d\n", actual_fno,face->v[0]+1,face->v[2]+1,face->v[1]+1);
    actual_fno++;
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwritePatchAscii(MRI_SURFACE *mris, char *fname)
{
  FILE    *fp ;
  int     vno, fno, n, nvertices, nfaces, type ;
  VERTEX  *v ;
  FACE    *face ;

  type = MRISfileNameType(fname) ;
#if 0
  if (type == MRIS_ASCII_TRIANGLE_FILE)
    return(MRISwriteAscii(mris, fname)) ;
#endif

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISwritePatchAscii: could not open file %s",
                 fname));
                 
  for (nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nvertices++ ;
  }
  for (nfaces = fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    nfaces++ ;
  }
  fprintf(fp, "#!ascii version of patch %s\n", mris->fname) ;
  fprintf(fp, "%d %d\n", nvertices, nfaces) ;
  fprintf(stdout, "nvertices=%d (valid=%d) nfaces=%d\n", nvertices, 
          mrisValidVertices(mris), nfaces) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    fprintf(fp, "%d\n", vno) ;
    fprintf(fp, "%f  %f  %f\n", v->x, v->y, v->z) ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    fprintf(fp, "%d\n", fno) ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      fprintf(fp, "%d ", face->v[n]) ;
    fprintf(fp, "\n") ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static MRI_SURFACE *
mrisReadAsciiFile(char *fname)
{
  MRI_SURFACE   *mris ;
  char    line[STRLEN], *cp ;
  int     vno, fno, n, nvertices, nfaces, patch, rip ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
              (ERROR_NOFILE, 
               "MRISreadAsciiFile: could not open file %s",fname));

  patch = 0 ;
  cp = fgetl(line, STRLEN, fp) ;
  sscanf(cp, "%d %d\n", &nvertices, &nfaces) ;
  mris = MRISalloc(nvertices, nfaces) ;
#if 0
  mris->type = MRIS_ASCII_TRIANGLE_FILE ;
#else
  mris->type = MRIS_TRIANGULAR_SURFACE ;
#endif
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    fscanf(fp, "%f  %f  %f  %d\n", &v->x, &v->y, &v->z, &rip) ;
    v->ripflag = rip ;
    if (v->ripflag)
      patch = 1 ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      fscanf(fp, "%d ", &face->v[n]) ;
      mris->vertices[face->v[n]].num++;
    }
    fscanf(fp, "%d\n", &rip) ;
    face->ripflag = rip ;
  }

  mris->patch = patch ;
  if (patch)
    mris->status = MRIS_PLANE ;
  fclose(fp) ;
  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static MRI_SURFACE *
mrisReadGeoFile(char *fname)
{
  MRI_SURFACE   *mris ;
  char    line[202], *cp ;
  int     vno, fno, n, nvertices, nfaces, patch, vertices_per_face, nedges ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
              (ERROR_NOFILE, 
               "mrisReadGeoFile: could not open file %s",fname));

  patch = 0 ;
  cp = fgetl(line, 100, fp) ;
  if (sscanf(cp, "%*d %d %d %d\n", &nvertices, &nfaces, &nedges) != 3)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "mrisReadGeoFile(%s): could not scan "
                       "dimensions from '%s'", fname, cp)) ;
  }
  vertices_per_face = nedges / nfaces ;
  if (vertices_per_face != VERTICES_PER_FACE)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, 
                       "mrisReadGeoFile(%s): unsupported vertices/face %d.",
                       fname, vertices_per_face)) ;
  }

  cp = fgetl(line, 200, fp) ;   /* nfaces again */
  mris = MRISalloc(nvertices, nfaces) ;
  mris->type = MRIS_GEO_TRIANGLE_FILE ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    fscanf(fp, "%e %e %e", &v->x, &v->y, &v->z) ;
    if (ISODD(vno))
      fscanf(fp, "\n") ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    int tmp ;
    face = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE-1 ; n++)
    {
      fscanf(fp, "%d ", &face->v[n]) ; face->v[n]-- ;   /* make it 0-based */
      if (face->v[n] < 0)
        DiagBreak() ;
      mris->vertices[face->v[n]].num++;
    }
    n = VERTICES_PER_FACE-1 ;   /* already true - but make it explicit */
    fscanf(fp, "-%d\n", &face->v[n]) ; face->v[n]-- ;   /* make it 0-based */
    mris->vertices[face->v[n]].num++; 

    /* swap positions so normal (via cross-product) will point outwards */
    tmp = face->v[1] ;
    face->v[1] = face->v[2] ; face->v[2] = tmp ;
  }

  fclose(fp) ;
  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisReadGeoFilePositions(MRI_SURFACE *mris, char *fname)
{
  char    line[202], *cp ;
  int     vno, nvertices, nfaces, patch, vertices_per_face, nedges ;
  VERTEX  *v ;
  FILE    *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE, 
               "mrisReadGeoFile: could not open file %s",fname));

  patch = 0 ;
  cp = fgetl(line, 100, fp) ;
  if (sscanf(cp, "%*d %d %d %d\n", &nvertices, &nfaces, &nedges) != 3)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_BADFILE, 
                (ERROR_BADFILE, "mrisReadGeoFile(%s): could not scan "
                 "dimensions from '%s'", fname, cp)) ;
  }
  vertices_per_face = nedges / nfaces ;
  if (vertices_per_face != VERTICES_PER_FACE)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, 
                       "mrisReadGeoFile(%s): unsupported vertices/face %d.",
                       fname, vertices_per_face)) ;
  }

  cp = fgetl(line, 200, fp) ;   /* nfaces again */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    fscanf(fp, "%e %e %e", &v->x, &v->y, &v->z) ;
    if (ISODD(vno))
      fscanf(fp, "\n") ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static MRI_SURFACE *
mrisReadAsciiPatchFile(char *fname)
{
  MRI_SURFACE   *mris ;
  char    line[STRLEN], *cp ;
  int     vno, fno, n, nvertices, nfaces ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
              (ERROR_NOFILE, 
               "MRISreadAsciiFile: could not open file %s",fname));

  cp = fgetl(line, 100, fp) ;
  sscanf(cp, "%d %d\n", &nvertices, &nfaces) ;
  mris = MRISalloc(nvertices, nfaces) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    fscanf(fp, "%f  %f  %f\n", &v->x, &v->y, &v->z) ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      fscanf(fp, "%d ", &face->v[n]) ;
      mris->vertices[face->v[n]].num++;
    }
    fscanf(fp, "\n") ;
  }

  fclose(fp) ;
  return(mris) ;
}
#endif
double
MRIScomputeCorrelationError(MRI_SURFACE *mris,MRI_SP *mrisp_template,int fno)
{
  INTEGRATION_PARMS  parms ;
  float              error ;

  if (!mrisp_template)
    return(0.0) ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.mrisp_template = mrisp_template ;
  parms.l_corr = 1.0f ;
  parms.frame_no = fno ;
  error = mrisComputeCorrelationError(mris, &parms, 1) ;
  return(sqrt(error / (double)mrisValidVertices(mris))) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeCorrelationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                            int use_stds)
{
  double   src, target, sse, delta, std ;
  VERTEX   *v ;
  int      vno ;
  float    x, y, z, l_corr ;

  l_corr = parms->l_corr + parms->l_pcorr ;  /* only one will be nonzero */
  if (FZERO(l_corr))
    return(0.0) ;

  for (sse = 0.0f, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    x = v->x ; y = v->y ; z = v->z ;
#if 0
    src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
    src = v->curv ;
#endif
    target = MRISPfunctionVal(parms->mrisp_template, mris, x, y, z, 
                              parms->frame_no) ;
#define DEFAULT_STD  4.0f
#define DISABLE_STDS  0
#if DISABLE_STDS
std = 1.0f ;
#else
    std = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,parms->frame_no+1);
    std = sqrt(std) ;
    if (FZERO(std))
      std = DEFAULT_STD /*FSMALL*/ ;
    if (!use_stds)
      std = 1.0f ;
#endif
    delta = (src - target) / std ;
    if (!finite(target) || !finite(delta))
      DiagBreak() ;
    sse += delta * delta ;
  }
  return(sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define D_DIST  0.25

static int
mrisComputeCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double   du, dv, up1, um1, vp1, vm1, delta, src, target, mag, max_mag,l_corr;
  VERTEX   *v ;
  int      vno, fno ;
  float    x, y, z, e1x, e1y, e1z, e2x, e2y, e2z, ux, uy, uz, vx, vy, vz, 
           std, coef;

  l_corr = parms->l_corr ;
  if (FZERO(l_corr))
    return(NO_ERROR) ;
  fno = parms->frame_no ;
  mrisComputeTangentPlanes(mris) ;
  max_mag = 0.0f ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRISPwrite(parms->mrisp_template, "temp.hipl") ;
    MRISPwrite(parms->mrisp, "srf.hipl") ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    e1x = v->e1x ;  e1y = v->e1y ; e1z = v->e1z ;
    e2x = v->e2x ;  e2y = v->e2y ; e2z = v->e2z ;
    x = v->x ; y = v->y ; z = v->z ;
#if 0
    src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
    src = v->curv ;
#endif
    target = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,
                              parms->frame_no);
    std = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,parms->frame_no+1);
    std = sqrt(std) ;
    if (FZERO(std))
      std = DEFAULT_STD /*FSMALL*/ ;
#if DISABLE_STDS
std = 1.0f ;
#endif
    delta = (target-src) ; coef = delta * l_corr / std ;

    /* now compute gradient of template w.r.t. a change in vertex position */

    /*
      sample the curvature functions along the tangent plane axes and
      compute the derivates using them.
      */
    ux = e1x*D_DIST ; uy = e1y*D_DIST ; uz = e1z*D_DIST ;
    vx = e2x*D_DIST ; vy = e2y*D_DIST ; vz = e2z*D_DIST ;

#if 0
    /* compute src term */
    up1 = MRISPfunctionVal(parms->mrisp, mris, x+ux, y+uy, z+uz, fno) ;
    um1 = MRISPfunctionVal(parms->mrisp, mris, x-ux, y-uy, z-uz, fno) ;
    vp1 = MRISPfunctionVal(parms->mrisp, mris, x+vx, y+vy, z+vz, fno) ;
    vm1 = MRISPfunctionVal(parms->mrisp, mris, x-vx, y-vy, z-vz, fno) ;
    du = (up1 - um1) / (2 * D_DIST) ;
    dv = (vp1 - vm1) / (2 * D_DIST) ;
    v->dx += coef * (du*e1x + dv*e2x) ;  /* in negative of grad. direction */
    v->dy += coef * (du*e1y + dv*e2y) ;
    v->dz += coef * (du*e1z + dv*e2z) ;
#endif

    /* compute target term */
    up1 = MRISPfunctionVal(parms->mrisp_template, mris, x+ux, y+uy, z+uz, fno);
    um1 = MRISPfunctionVal(parms->mrisp_template, mris, x-ux, y-uy, z-uz, fno);
    vp1 = MRISPfunctionVal(parms->mrisp_template, mris, x+vx, y+vy, z+vz, fno);
    vm1 = MRISPfunctionVal(parms->mrisp_template, mris, x-vx, y-vy, z-vz, fno);
    du = (up1 - um1) / (2 * D_DIST) ;
    dv = (vp1 - vm1) / (2 * D_DIST) ;
    v->dx -= coef * (du*e1x + dv*e2x) ;  
    v->dy -= coef * (du*e1y + dv*e2y) ;
    v->dz -= coef * (du*e1z + dv*e2z) ;

    mag = sqrt(v->dx*v->dx + v->dy*v->dy + v->dz*v->dz) ;
    if (mag > max_mag)
      max_mag = mag ;
    if (!finite(v->dx) || !finite(v->dy) || !finite(v->dz))
    {
      DiagBreak() ;
      ErrorExit(ERROR_BADPARM, 
                "mrisComputeCorrelationTerm: delta is not finite at vno %d", 
                vno) ;
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "max gradient magnitude = %2.5f\n", max_mag) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define D_ANGLE  RADIANS(0.25)

static int
mrisComputePolarCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double   ap1, am1, da, bp1, bm1, db, gp1, gm1, dg, delta, src, target, mag, 
           max_mag ;
  VERTEX   *v ;
  int      vno, fno ;
  float    x, y, z, std, coef, dx, dy, dz, nv, r ;

  if (FZERO(parms->l_pcorr))
    return(NO_ERROR) ;
  fno = parms->frame_no ;
  max_mag = 0.0f ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRISPwrite(parms->mrisp_template, "temp.hipl") ;
    /*    MRISPwrite(parms->mrisp, "srf.hipl") ;*/
  }
  mris->gamma = mris->beta = mris->alpha = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    x = v->x ; y = v->y ; z = v->z ;
#if 0
    src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
    src = v->curv ;
#endif
    target = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,
                              parms->frame_no);
    std = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,parms->frame_no+1);
    std = sqrt(std) ;
    if (FZERO(std))
      std = DEFAULT_STD /*FSMALL*/ ;
#if DISABLE_STDS
std = 1.0f ;
#endif
    delta = (target-src) ; coef = delta / std ;

    /* now compute gradient of template w.r.t. a change in vertex position */

    /* 
       compute the gradient by using differential rotations in the 3 
       rotational directions using the associated skew-symmetric
       differential rotation matrices
       */
    /* compute alpha term - differential rotation around z axis */
    dx = y*D_ANGLE ; dy = -x*D_ANGLE ; dz = 0 ;
    am1 = MRISPfunctionVal(parms->mrisp_template, mris, x-dx, y-dy, z-dz, 0) ;
    ap1 = MRISPfunctionVal(parms->mrisp_template, mris, x+dx, y+dy, z+dz, 0) ;
    da = (ap1 - am1) / (2*D_ANGLE);

    /* compute beta term - differential rotation around y axis */
    dx = -z*D_ANGLE ; dy = 0 ; dz = x*D_ANGLE ;
    bm1 = MRISPfunctionVal(parms->mrisp_template, mris, x-dx, y-dy, z-dz, 0) ;
    bp1 = MRISPfunctionVal(parms->mrisp_template, mris, x+dx, y+dy, z+dz, 0) ;
    db = (bp1 - bm1) / (2*D_ANGLE);

    /* compute gamma term - differential rotation around x axis */
    dx = 0 ; dy = -z*D_ANGLE ; dz = y*D_ANGLE ;
    gm1 = MRISPfunctionVal(parms->mrisp_template, mris, x-dx, y-dy, z-dz, 0) ;
    gp1 = MRISPfunctionVal(parms->mrisp_template, mris, x+dx, y+dy, z+dz, 0) ;
    dg = (gp1 - gm1) / (2*D_ANGLE);

    mris->gamma -= coef * dg ;   /* around x-axis */
    mris->beta  -= coef * db ;   /* around y-axis */
    mris->alpha -= coef * da ;   /* around z-axis */

    mag = sqrt(v->dx*v->dx + v->dy*v->dy + v->dz*v->dz) ;
    if (mag > max_mag)
      max_mag = mag ;
    if (!finite(v->dx) || !finite(v->dy) || !finite(v->dz))
    {
      DiagBreak() ;
      ErrorExit(ERROR_BADPARM, 
             "mrisComputePolarCorrelationTerm: delta is not finite at vno %d", 
                vno) ;
    }
  }

  nv = mrisValidVertices(mris) ;
  r = mris->radius ;
  mris->alpha /= nv ;
  mris->beta /= nv ;
  mris->gamma /= nv ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->dx = r*mris->alpha ; v->dy = r*mris->beta ; v->dz = r*mris->gamma ;
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "max gradient magnitude = %2.5f\n", max_mag) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISzeroMeanCurvature(MRI_SURFACE *mris)
{
  double    mean ;
  int       vno, vtotal ;
  VERTEX    *v ;

  for (mean = 0.0f, vtotal = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    mean += v->curv ;
    vtotal++ ;
  }

  mean /= (double)vtotal ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "curvature mean = %2.3f\n", mean) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv -= mean ;
  }

  mrisComputeCurvatureValues(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISnormalizeCurvature(MRI_SURFACE *mris)
{
  double    mean, var, std ;
  int       vno, vtotal ;
  VERTEX    *v ;

  for (mean = 0.0f, vtotal = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    mean += v->curv ;
    vtotal++ ;
  }

  mean /= (double)vtotal ;

  for (var = 0.0f, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    std = v->curv - mean ;
    var += std * std ;
  }

  var /= (double)vtotal ;
  std = sqrt(var) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "curvature mean = %2.3f, std = %2.3f\n", mean, std) ;

  /* now normalize the curvatures so they have unit standard deviation, but
     leave the mean alone */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = (v->curv - mean) / std /* + mean*/ ;
  }

  mrisComputeCurvatureValues(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISnormalizeCurvatureVariance(MRI_SURFACE *mris)
{
  double    mean, var, std ;
  int       vno, vtotal ;
  VERTEX    *v ;

  for (mean = 0.0f, vtotal = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    mean += v->curv ;
    vtotal++ ;
  }

  mean /= (double)vtotal ;

  for (var = 0.0f, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    std = v->curv - mean ;
    var += std * std ;
  }

  var /= (double)vtotal ;
  std = sqrt(var) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "curvature mean = %2.3f, std = %2.3f\n", mean, std) ;

  /* now normalize the curvatures so they have unit standard deviation, but
     leave the mean alone */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = (v->curv - mean) / std  + mean ;
  }

  mrisComputeCurvatureValues(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISscaleDistances(MRI_SURFACE *mris, float scale)
{
  int       vno, n ;
  VERTEX    *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    for (n = 0 ; n < v->vtotal ; n++)
      v->dist[n] *= scale ;
  }
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisCheck(MRI_SURFACE *mris)
{
  int       vno ; 
  VERTEX    *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (fabs(v->curv) > 10000.0)
    {
      DiagBreak() ;
      return(ERROR_BADPARM) ;
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISminFilterCurvatures(MRI_SURFACE *mris, int niter)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  curv ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < niter ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      curv = v->curv ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        if (vn->curv < curv)
          curv = vn->curv ;
      }
      v->tdx = curv ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->curv = v->tdx ;
    }
  }
  mrisComputeCurvatureValues(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmaxFilterCurvatures(MRI_SURFACE *mris, int niter)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  curv ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < niter ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      curv = v->curv ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        if (vn->curv > curv)
          curv = vn->curv ;
      }
      v->tdx = curv ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->curv = v->tdx ;
    }
  }
  mrisComputeCurvatureValues(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaverageCurvatures(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  curv, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      curv = v->curv ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        curv += vn->curv ;
      }
      num++ ;  /*  account for central vertex */
      v->tdx = curv / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->curv = v->tdx ;
    }
  }
  mrisComputeCurvatureValues(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaverageMarkedCurvatures(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  curv, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || !v->marked)
        continue ;
      curv = v->curv ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        curv += vn->curv ;
      }
      num++ ;  /*  account for central vertex */
      v->tdx = curv / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 0)
        continue ;
      v->curv = v->tdx ;
    }
  }
  mrisComputeCurvatureValues(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaverageMarkedVals(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  val, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 0)
        continue ;
      val = v->val ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag || vn->marked == 0)
          continue ;
        num++ ;
        val += vn->val ;
      }
      num++ ;  /*  account for central vertex */
      v->tdx = val / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 0)
        continue ;
      v->val = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaverageVals(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  val, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      val = v->val ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        val += vn->val ;
      }
      num++ ;  /*  account for central vertex */
      v->tdx = val / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->val = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
static int
compare_sort_vals(const void *pc1, const void *pc2)
{
  register float c1, c2 ;

  c1 = *(float *)pc1 ;
  c2 = *(float *)pc2 ;

/*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (c1 > c2)
    return(1) ;
  else if (c1 < c2)
    return(-1) ;

  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmedianFilterVals(MRI_SURFACE *mris, int nmedians)
{
  int    i, vno, vnb, *pnb, vnum, num ;
  float  val_list[MAX_NEIGHBORS] ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < nmedians ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      pnb = v->v ;
      vnum = v->vnum ;
			val_list[0] = v->val ;
      for (num = 1, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;

				val_list[num++] = vn->val ;
      }
			qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals) ;
			if (ISODD(num))
				v->tdx = val_list[(num-1)/2] ;
			else
				v->tdx = (val_list[num/2] + val_list[num/2-1])/2 ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->val = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmedianFilterCurvature(MRI_SURFACE *mris, int nmedians)
{
  int    i, vno, vnb, *pnb, vnum, num ;
  float  val_list[MAX_NEIGHBORS] ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < nmedians ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      pnb = v->v ;
      vnum = v->vnum ;
			val_list[0] = v->curv ;
      for (num = 1, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;

				val_list[num++] = vn->curv ;
      }
			qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals) ;
			if (ISODD(num))
				v->tdx = val_list[(num-1)/2] ;
			else
				v->tdx = (val_list[num/2] + val_list[num/2-1])/2 ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->curv = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmedianFilterVal2s(MRI_SURFACE *mris, int nmedians)
{
  int    i, vno, vnb, *pnb, vnum, num ;
  float  val_list[MAX_NEIGHBORS] ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < nmedians ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      pnb = v->v ;
      vnum = v->vnum ;
			val_list[0] = v->val2 ;
      for (num = 1, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;

				val_list[num++] = vn->val2 ;
      }
			qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals) ;
			if (ISODD(num))
				v->tdx = val_list[(num-1)/2] ;
			else
				v->tdx = (val_list[num/2] + val_list[num/2-1])/2 ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->val2 = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmedianFilterVal2baks(MRI_SURFACE *mris, int nmedians)
{
  int    i, vno, vnb, *pnb, vnum, num ;
  float  val_list[MAX_NEIGHBORS] ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < nmedians ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      pnb = v->v ;
      vnum = v->vnum ;
			val_list[0] = v->val2bak ;
      for (num = 1, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;

				val_list[num++] = vn->val2bak ;
      }
			qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals) ;
			if (ISODD(num))
				v->tdx = val_list[(num-1)/2] ;
			else
				v->tdx = (val_list[num/2] + val_list[num/2-1])/2 ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->val2bak = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaverageVal2s(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  val, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      val = v->val2 ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        val += vn->val2 ;
      }
      num++ ;  /*  account for central vertex */
      v->tdx = val / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->val2 = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaverageVal2baks(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  val, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      val = v->val2bak ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        val += vn->val2bak ;
      }
      num++ ;  /*  account for central vertex */
      v->tdx = val / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->val2bak = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISscaleCurvatures(MRI_SURFACE *mris, float min_curv, float max_curv)
{
  double    old_min_curv, old_max_curv, mean, scale ;
  int       vno, vtotal ;
  VERTEX    *v ;

  old_min_curv = 100000.0 ; old_max_curv = -100000.0 ;
  for (mean = 0.0, vtotal = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    vtotal++ ;
    mean += v->curv ;
    if (v->curv > old_max_curv)
      old_max_curv = v->curv ;
    if (v->curv < old_min_curv)
      old_min_curv = v->curv ;
  }

  mean /= (double)vtotal ;
  scale = (max_curv - min_curv) / (old_max_curv - old_min_curv) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = (v->curv - mean) * scale + mean ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISnonmaxSuppress(MRI_SURFACE *mris)
{
  double   du, dv, up1, um1, vp1, vm1, src, dx, dy, dz, fp1, fm1, mag ;
  VERTEX   *v ;
  int      vno ;
  float    x, y, z, e1x, e1y, e1z, e2x, e2y, e2z, ux, uy, uz, vx, vy, vz ;
  MRI_SP   *mrisp, *mrisp_blur ;

  mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
  mrisp_blur = MRISPblur(mrisp, NULL, 20.0, 0) ;
  mrisComputeTangentPlanes(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    e1x = v->e1x ;  e1y = v->e1y ; e1z = v->e1z ;
    e2x = v->e2x ;  e2y = v->e2y ; e2z = v->e2z ;
    x = v->x ; y = v->y ; z = v->z ;
#if 0
    src = v->curv ;
#else
    src = MRISPfunctionVal(mrisp, mris, x, y, z, 0) ;
#endif


    /* now compute gradient of template w.r.t. a change in vertex position */

    /*
      sample the curvature functions along the tangent plane axes and
      compute the derivates using them.
      */
    ux = 10*e1x ; uy = 10*e1y ; uz = 10*e1z ;
    vx = 10*e2x ; vy = 10*e2y ; vz = 10*e2z ;
    ux = e1x*D_DIST ; uy = e1y*D_DIST ; uz = e1z*D_DIST ;
    vx = e2x*D_DIST ; vy = e2y*D_DIST ; vz = e2z*D_DIST ;

    /* compute gradient usnig blurred image */
    up1 = MRISPfunctionVal(mrisp_blur, mris, x+ux, y+uy, z+uz, 0) ;
    um1 = MRISPfunctionVal(mrisp_blur, mris, x-ux, y-uy, z-uz, 0) ;
    vp1 = MRISPfunctionVal(mrisp_blur, mris, x+vx, y+vy, z+vz, 0) ;
    vm1 = MRISPfunctionVal(mrisp_blur, mris, x-vx, y-vy, z-vz, 0) ;
    du = (up1 - um1) / (2 * D_DIST) ;
    dv = (vp1 - vm1) / (2 * D_DIST) ;

    /* calculate curvature gradient */
    dx = (du*e1x + dv*e2x) ;  
    dy = (du*e1y + dv*e2y) ;
    dz = (du*e1z + dv*e2z) ;
    mag = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (FZERO(mag))   /* zero gradient */
      v->curv = 0 ;
    else
    { mag *= .1 ; dx = dx / mag ; dy = dy / mag ; dz = dz / mag ; }
    fp1 = MRISPfunctionVal(mrisp, mris, x+dx, y+dy, z+dz, 0) ;
    fm1 = MRISPfunctionVal(mrisp, mris, x-dx, y-dy, z-dz, 0) ;

    if ((src >= fp1) && (src >= fm1))       /* local max */
      v->curv = 1 ;
    else if ((src <= fp1) && (src <= fm1))  /* local min */
      v->curv = -1 ;
    else
      v->curv = 0 ;
    v->curv = src ;
  }
  MRISPfree(&mrisp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
				same as  below, but tracks odx,ody,odz fields
------------------------------------------------------*/
static int
mrisTrackTotalDistanceNew(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;
  float  nc ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nc = v->odx*v->nx + v->ody*v->ny + v->odz*v->nz ;
    v->curv += nc ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisTrackTotalDistance(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;
  float  nc ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nc = v->dx*v->nx + v->dy*v->ny + v->dz*v->nz ;
    v->curv += nc ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISclearCurvature(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = 0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISrigidBodyAlignLocal(MRI_SURFACE *mris, INTEGRATION_PARMS *old_parms)
{
  int                old_status, steps ;
  INTEGRATION_PARMS  parms ;

  /* dx,dy,dz interpreted as rotations in applyGradient when status is rigid */
  old_status = mris->status ;   /* okay, okay, this is a hack too... */
  mris->status = MRIS_RIGID_BODY ; 
  memset(&parms, 0, sizeof(parms)) ;
  parms.integration_type = INTEGRATE_LM_SEARCH ;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;

  parms.mrisp_template = old_parms->mrisp_template ;
  parms.fp = old_parms->fp ;
  parms.niterations =  25 ;
  parms.frame_no = old_parms->frame_no ;
  parms.mrisp = old_parms->mrisp ;
  parms.tol = old_parms->tol ;
  parms.l_pcorr = 1.0f ;
  parms.dt = old_parms->dt ;
  /*  parms.integration_type = old_parms->integration_type ;*/
  parms.momentum = old_parms->momentum ;
  parms.write_iterations = old_parms->write_iterations ;
  parms.start_t = old_parms->start_t ;
  strcpy(parms.base_name, old_parms->base_name) ;

  steps = MRISintegrate(mris, &parms, 0) ;
  old_parms->start_t += steps ;
  mris->status = old_status ;
  if (Gdiag & DIAG_WRITE)
    fprintf(old_parms->fp, "rigid alignment complete, sse = %2.3f\n",
            MRIScomputeSSE(mris, old_parms)) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    float area_rms, angle_rms, curv_rms, dist_rms, corr_rms, rms ;

    rms = 
      mrisComputeError(mris, &parms,&area_rms,&angle_rms,&curv_rms,&dist_rms,
                       &corr_rms);
    fprintf(stdout, "rms = %2.3f, corr_rms = %2.3f ", rms, corr_rms) ;
    rms = 
      mrisComputeError(mris, old_parms,&area_rms,&angle_rms,&curv_rms,
                       &dist_rms, &corr_rms);
    fprintf(stdout, "(%2.3f, %2.3f)\n", rms, corr_rms) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define STARTING_ANGLE   RADIANS(16.0f)
#define ENDING_ANGLE     RADIANS(4.0f)
#define NANGLES          8

int
MRISrigidBodyAlignGlobal(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                         float min_degrees, float max_degrees, int nangles)
{
  double   alpha, beta, gamma, degrees, delta, mina, minb, ming,
           sse, min_sse ;
  int      old_status = mris->status ;

  min_degrees = RADIANS(min_degrees) ; max_degrees = RADIANS(max_degrees) ;
  mrisOrientSurface(mris) ;
  mris->status = MRIS_RIGID_BODY ; 
  if (!parms->start_t)
  {
    mrisLogStatus(mris, parms, stderr, 0.0f) ;
    if (Gdiag & DIAG_WRITE)
    {
      mrisLogStatus(mris, parms, parms->fp, 0.0f) ;
      if (parms->write_iterations > 0)
        mrisWriteSnapshot(mris, parms, 0) ;
    }
  }
  for (degrees = max_degrees ; degrees >= min_degrees ; degrees /= 2.0f)
  {
    mina = minb = ming = 0.0 ;
    min_sse = mrisComputeCorrelationError(mris, parms, 1) ;  /* was 0 !!!! */
    delta = 2*degrees / (float)nangles ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "scanning %2.2f degree nbhd, min sse = %2.2f\n", 
              (float)DEGREES(degrees), (float)min_sse) ;
    for (alpha = -degrees ; alpha <= degrees ; alpha += delta)
    {
      for (beta = -degrees ; beta <= degrees ; beta += delta)
      {
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "\r(%+2.2f, %+2.2f, %+2.2f), "
                  "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                  (float)DEGREES(alpha), (float)DEGREES(beta), (float)
                  DEGREES(-degrees), (float)DEGREES(mina), 
                  (float)DEGREES(minb), (float)DEGREES(ming),(float)min_sse);

        for (gamma = -degrees ; gamma <= degrees ; gamma += delta)
        {
          MRISsaveVertexPositions(mris, TMP_VERTICES) ;
          MRISrotate(mris, mris, alpha, beta, gamma) ;
          sse = mrisComputeCorrelationError(mris, parms, 1) ;  /* was 0 !!!! */
          MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
          if (sse < min_sse)
          {
            mina = alpha ; minb = beta ; ming = gamma ;
            min_sse = sse ;
          }
#if 0
          if (Gdiag & DIAG_SHOW)
            fprintf(stdout, "\r(%+2.2f, %+2.2f, %+2.2f), "
                    "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                    (float)DEGREES(alpha), (float)DEGREES(beta), (float)
                    DEGREES(gamma), (float)DEGREES(mina), 
                    (float)DEGREES(minb), (float)DEGREES(ming),(float)min_sse);
#endif
        }
      }
    }
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "\n") ;
    if (!FZERO(mina) || !FZERO(minb) || !FZERO(ming))
    {
      MRISrotate(mris, mris, mina, minb, ming) ;
      sse = mrisComputeCorrelationError(mris, parms, 1) ;  /* was 0 !!!! */
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "min sse = %2.2f at (%2.2f, %2.2f, %2.2f)\n",
                sse, (float)DEGREES(mina), (float)DEGREES(minb), 
                (float)DEGREES(ming)) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp, 
                "rotating brain by (%2.2f, %2.2f, %2.2f), sse: %2.2f\n",
                (float)DEGREES(mina), (float)DEGREES(minb), 
                (float)DEGREES(ming), (float)sse) ;
      parms->start_t += 1.0f ;
      parms->t += 1.0f ;
      if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
        mrisWriteSnapshot(mris, parms, parms->start_t) ;
      if (Gdiag & DIAG_WRITE)
      mrisLogStatus(mris, parms, parms->fp, 0.0f) ;
      if (Gdiag & DIAG_SHOW)
        mrisLogStatus(mris, parms, stderr, 0.0f) ;
    }
  }

  mris->status = old_status ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRISrmsTPHeight(MRI_SURFACE *mris)
{
  int    vno, i, total_nbrs ;
  VERTEX *vertex, *vnb ;
  double avg_height, dot, nx, ny, nz, x, y, z, d ;

  if (mris->status == MRIS_PLANE)
    return(NO_ERROR) ;

  avg_height = 0.0 ; total_nbrs = 0 ;
  mrisComputeTangentPlanes(mris) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;

    if (vertex->vtotal <= 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    
    nx = vertex->nx ; ny = vertex->ny ; nz = vertex->nz ;

    for (i = 0 ; i < vertex->vtotal ; i++)
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag)
        continue ;
      x = vnb->x-vertex->x ; y = vnb->y-vertex->y ; z = vnb->z-vertex->z ;
      d = (x*x+y*y+z*z) ;
      if (FZERO(d))
        continue ;
/* 
   calculate the projection of this vertex onto the surface normal
*/
      dot = nx*x + ny*y + nz*z ;  /* height above TpS */
      avg_height += dot*dot / d ;
      total_nbrs++ ;
    }
  }

  return(sqrt(avg_height / (double)total_nbrs)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISzeroNegativeAreas(MRI_SURFACE *mris)
{
  int     vno, fno ;
  VERTEX  *v ;
  FACE    *face ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->area < 0)
      v->area = 0 ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    if (face->area < 0.0f)
      face->area = 0.0f ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISfindClosestVertex(MRI_SURFACE *mris, float x, float y, float z)
{
  int    vno, min_v = -1 ;
  VERTEX *v ;
  float  d, min_d, dx, dy, dz ;

  min_d = 10000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dx = v->x - x ; dy = v->y - y ; dz = v->z - z ;
    d = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (d < min_d)
    {
      min_d = d ;
      min_v = vno ;
    }
  }

  return(min_v) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISfindClosestOriginalVertex(MRI_SURFACE *mris, float x, float y, float z)
{
  int    vno, min_v = -1 ;
  VERTEX *v ;
  float  d, min_d, dx, dy, dz ;

  min_d = 10000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == 91007 || vno == 91814)
      DiagBreak() ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dx = v->origx - x ; dy = v->origy - y ; dz = v->origz - z ;
    d = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (d < min_d)
    {
      min_d = d ;
      min_v = vno ;
    }
  }

  return(min_v) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISfindClosestCanonicalVertex(MRI_SURFACE *mris, float x, float y, float z)
{
  int    vno, min_v = -1 ;
  VERTEX *v ;
  float  d, min_d, dx, dy, dz ;

  min_d = 10000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dx = v->cx - x ; dy = v->cy - y ; dz = v->cz - z ;
    d = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (d < min_d)
    {
      min_d = d ;
      min_v = vno ;
    }
  }

  return(min_v) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISuseCurvatureDifference(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *vertex ;
  float   kmin, kmax ;

  kmin = 100000.0f ; kmax = -100000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->curv = fabs(vertex->k1 - vertex->k2) ;
    if (vertex->curv > kmax)
      kmax = vertex->curv ;
    if (vertex->curv < kmin)
      kmin = vertex->curv ;
  }

  mris->min_curv = kmin ; mris->max_curv = kmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISuseCurvatureMax(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;
  float   kmin, kmax ;

  kmin = 100000.0f ; kmax = -100000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (fabs(v->k1) > fabs(v->k2))
      v->curv = v->k1 ;
    else
      v->curv = v->k2 ;
    /*    v->curv = MAX(fabs(v->k1), fabs(v->k2)) ;*/
    if (v->curv > kmax)
      kmax = v->curv ;
    if (v->curv < kmin)
      kmin = v->curv ;
    if (v->curv < 0)
      DiagBreak() ;
  }

  fprintf(stdout, "kmin = %2.2f, kmax = %2.2f\n", kmin, kmax) ;
  mris->min_curv = kmin ; mris->max_curv = kmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISuseCurvatureMin(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;
  float   kmin, kmax ;

  kmin = 100000.0f ; kmax = -100000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (fabs(v->k1) > fabs(v->k2))
      v->curv = v->k2 ;
    else
      v->curv = v->k1 ;
    /*    v->curv = MIN(v->k1, v->k2) ;*/
    if (v->curv > kmax)
      kmax = v->curv ;
    if (v->curv < kmin)
      kmin = v->curv ;
    if (v->curv < 0)
      DiagBreak() ;
  }

  mris->min_curv = kmin ; mris->max_curv = kmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISuseCurvatureStretch(MRI_SURFACE *mris)
{
  int    vno, n ;
  VERTEX *v ;
  float   kmin, kmax, dist, dist_orig, curv, dist_scale, max_stretch, stretch ;

  dist_scale = sqrt(mris->orig_area / mris->total_area) ;
  kmin = 100000.0f ; kmax = -100000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    max_stretch = 0.0f ;
    for (curv = 0.0f, n = 0 ; n < v->vtotal ; n++)
    {
      dist = dist_scale * v->dist[n] ; dist_orig = v->dist_orig[n] ;
      stretch = dist - dist_orig ;
      if (stretch > max_stretch)
        max_stretch = stretch ;
    }
    v->curv = max_stretch ;
    if (v->curv > kmax)
      kmax = v->curv ;
    if (v->curv < kmin)
      kmin = v->curv ;
    if (v->curv < 0)
      DiagBreak() ;
  }

  mris->min_curv = kmin ; mris->max_curv = kmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISuseNegCurvature(MRI_SURFACE *mris)
{
  int    vno, fno ;
  VERTEX *v ;
  FACE   *f ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = 0 ;
    for (fno = 0 ; fno < v->num ; fno++)
    {
      f = &mris->faces[v->f[fno]] ;
      if (f->area < 0.0f)
        v->curv = 1.0f ;
    }
  }

  mris->min_curv = 0 ; mris->max_curv = 1 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIStalairachToVertex(MRI_SURFACE *mris, Real xt, Real yt, Real zt)
{
  int     vno ;
  Real    xw, yw, zw ;

  transform_point(mris->inverse_linear_transform, xt, yt, zt, &xw, &yw, &zw) ;
  vno = MRISfindClosestOriginalVertex(mris, xw, yw, zw) ;
  return(vno) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScanonicalToVertex(MRI_SURFACE *mris, Real phi, Real theta)
{
  int     vno ;
  Real    xw, yw, zw ;

  MRIScanonicalToWorld(mris, phi, theta, &xw, &yw, &zw);
  vno = MRISfindClosestCanonicalVertex(mris, xw, yw, zw) ;
  return(vno) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISworldToTalairachVoxel(MRI_SURFACE *mris, MRI *mri, Real xw, Real yw, 
                          Real zw, Real *pxv, Real *pyv, Real *pzv)
{
  Real  xt, yt, zt ;

  transform_point(mris->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  MRIworldToVoxel(mri, xt, yt, zt, pxv, pyv, pzv) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISvertexToVoxel(VERTEX *v, MRI *mri,Real *pxv, Real *pyv, Real *pzv)
{
  Real  xw, yw, zw ;

  xw = v->x ; yw = v->y ; zw = v->z ;
  // MRIworldToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;
  MRIsurfaceRASToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISorigVertexToVoxel(VERTEX *v, MRI *mri,Real *pxv, Real *pyv, Real *pzv)
{
  Real  xw, yw, zw ;

  xw = v->origx ; yw = v->origy ; zw = v->origz ;
  // MRIworldToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;
  MRIsurfaceRASToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaverageVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  x, y, z, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      x = v->x ; y = v->y ; z = v->z ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        x += vn->x ; y += vn->y ; z += vn->z ;
      }
      num++ ;   /* account for central vertex */
      v->tdx = x / num ;
      v->tdy = y / num ;
      v->tdz = z / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->x = v->tdx ; v->y = v->tdy ; v->z = v->tdz ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaverageEveryOtherVertexPositions(MRI_SURFACE *mris, int navgs, int which)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  x, y, z, num ;
  VERTEX *v, *vn ;

  which = ISODD(which) ;
  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = which ; vno < mris->nvertices ; vno += 2)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      x = v->x ; y = v->y ; z = v->z ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        x += vn->x ; y += vn->y ; z += vn->z ;
      }
      num++ ;   /* account for central vertex */
      v->tdx = x / num ;
      v->tdy = y / num ;
      v->tdz = z / num ;
    }
    for (vno = which ; vno < mris->nvertices ; vno += 2)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->x = v->tdx ; v->y = v->tdy ; v->z = v->tdz ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/

#define MAX_TOTAL_MOVEMENT 5.0
#define MAX_MOVEMENT       .2 
#define DELTA_M            (MAX_MOVEMENT/2.0)


#define DEBUG_V            33100

#define MAX_REDUCTIONS     2
#define REDUCTION_PCT      0.5

int   
MRISpositionSurfaces(MRI_SURFACE *mris, MRI **mri_flash, int nvolumes, 
										 INTEGRATION_PARMS *parms)
{
  /*  char   *cp ;*/
  int    niterations, n, write_iterations, nreductions = 0, ripped = 0, increased = 0 ;
  double pial_sse, sse, wm_sse, delta_t = 0.0, dt, l_intensity, base_dt, last_sse, rms, mle_sse, last_mle_sse, 
         pct_sse_decrease, l_repulse, l_surf_repulse, last_wm_sse, last_pial_sse ;
  MHT    *mht = NULL ;
  struct timeb  then ; int msec ;
  MRI    *mri_brain = mri_flash[0] ;

  base_dt = parms->dt ;
  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  TimerStart(&then) ;
  parms->mri_smooth = parms->mri_brain = mri_brain ;
  niterations = parms->niterations ; write_iterations = parms->write_iterations ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;

    if (!parms->fp)
    {
      sprintf(fname, "%s.%s.out", 
              mris->hemisphere==RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      else
        parms->fp = fopen(fname, "a") ;
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
      Progname, fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  mrisClearMomentum(mris) ;
  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;

  MRIScomputeNormals(mris) ;
  mrisClearDistances(mris) ;

  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  wm_sse = MRIScomputeSSEExternal(mris, parms, &mle_sse) ;
  MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  pial_sse = MRIScomputeSSE(mris, parms) ;
  sse = last_sse = wm_sse + pial_sse ; last_mle_sse = mle_sse ;
#if 0
  rms = (*gMRISexternalRMS)(mris, parms) ;
#else
  rms = sqrt(mle_sse/(float)mris->nvertices) ;
#endif

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "%3.3d: dt: %2.4f, sse=%2.2f, pial sse=%2.2f, wm sse=%2.2f, rms=%2.2f\n", 
            0, 0.0f, (float)sse/(float)mris->nvertices, 
            (float)pial_sse/(float)mris->nvertices, 
            (float)wm_sse/(float)mris->nvertices, (float)rms);
																/*  */
  if (Gdiag & DIAG_WRITE)
  {
    fprintf(parms->fp, "%3.3d: dt: %2.4f, sse=%2.2f, pial sse=%2.2f, wm sse=%2.2f, rms=%2.2f\n", 
            0, 0.0f, 
            (float)sse/(float)mris->nvertices, 
            (float)pial_sse/(float)mris->nvertices, 
            (float)wm_sse/(float)mris->nvertices, rms);
    fflush(parms->fp) ;
  }

  /* write out initial surface */
  if ((parms->write_iterations > 0) && (Gdiag&DIAG_WRITE) && !parms->start_t)
    mrisWriteSnapshots(mris, parms, 0) ;

  dt = parms->dt ; l_intensity = parms->l_intensity ;
  mris->noscale = TRUE ; l_repulse = parms->l_repulse ;
  l_surf_repulse = parms->l_surf_repulse ;
  for (n = parms->start_t ; n < parms->start_t+niterations ; n++)
  {
    /* compute and apply wm derivatative */
    MRISclearGradient(mris) ; mrisClearExtraGradient(mris) ;
		if (!increased)
		{
			if (gMRISexternalClearSSEStatus)
				(*gMRISexternalClearSSEStatus)(mris) ;
		}
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;  /* make wm positions current */
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;          
    MRIScomputeMetricProperties(mris) ;
    if (gMRISexternalGradient)
      mle_sse = (*gMRISexternalGradient)(mris, parms) ;   /* this computes the external sse for both wm and pial */
		if (increased && gMRISexternalReduceSSEIncreasedGradients)
		{
			printf("decreasing gradient at vertices with delta SSE>0 by %2.2f\n", 0.5/increased) ;
			(*gMRISexternalReduceSSEIncreasedGradients)(mris, 0.5/increased) ;
			if (gMRISexternalClearSSEStatus)
				(*gMRISexternalClearSSEStatus)(mris) ;
		}

    parms->l_repulse = l_repulse ;  /* use self-repulsion for wm surface */
    parms->l_surf_repulse = 0    ;  /* don't repel wm surface outwards from itself */
    mrisComputePositioningGradients(mris, parms) ;
		if (!FZERO(parms->l_link))
			mrisComputeLinkTerm(mris, parms->l_link, 0) ;
    if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
      mht = MHTfillTable(mris, mht) ;
    last_wm_sse = MRIScomputeSSEExternal(mris,parms, &last_mle_sse) ;

    delta_t = mrisAsynchronousTimeStepNew(mris, 0, dt, mht, MAX_ASYNCH_NEW_MM) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    if (gMRISexternalTimestep)
      (*gMRISexternalTimestep)(mris, parms) ;
    wm_sse = MRIScomputeSSE(mris, parms) ;  /* needs update orig to compute sse - will compute external sse later */

    /* store current wm positions in WHITE vertices, and pial in INFLATED vertices for undo */
    MRISclearGradient(mris) ;
    MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;  /* make pial positions current */
    MRISsaveVertexPositions(mris, INFLATED_VERTICES) ; /* pial->inflated */
    MRIScomputeMetricProperties(mris) ; 
    MRISrestoreExtraGradients(mris) ;    /* put pial deltas into v->d[xyz] */
    parms->l_repulse = 0 ;  /* don't use self-repulsion for pial surface */
    parms->l_surf_repulse = l_surf_repulse ;  /* repel pial surface out from wm */
    mrisComputePositioningGradients(mris, parms) ;
		if (!FZERO(parms->l_link))
			mrisComputeLinkTerm(mris, parms->l_link, 1) ;
    if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
      mht = MHTfillTable(mris, mht) ;
    last_pial_sse = MRIScomputeSSE(mris,parms) ;
    delta_t += mrisAsynchronousTimeStepNew(mris, 0, dt, mht, MAX_ASYNCH_NEW_MM) ;
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    if (gMRISexternalTimestep)
      (*gMRISexternalTimestep)(mris, parms) ;
    delta_t /= 2 ;
    pial_sse = MRIScomputeSSEExternal(mris, parms, &mle_sse) ;  /* needs update pial to compute sse. mle_sse includes wm and pial */
		printf("MLE sse %2.3f --> %2.3f, delta = %2.3f\n", last_mle_sse, mle_sse, mle_sse-last_mle_sse) ;
    last_sse = last_wm_sse + last_pial_sse + last_mle_sse ;
    sse = wm_sse + pial_sse ; /* pial sse includes current mle_sse */

    pct_sse_decrease = 1 - sse/last_sse ;
		pct_sse_decrease = 1 - mle_sse / last_mle_sse ;  /* only terminate if surfaces have asymptoted to desired positions */
    if (pct_sse_decrease < parms->tol)  /* error didn't decrease much */
    {
      nreductions++ ;
      dt *= .5 ;

      if (pct_sse_decrease < 0)  /* error increased - reject time step */
      {
				increased++ ;
        printf("error increased by %2.3f%% - time step reduction #%d: dt=%2.3f, undoing step...\n",
               -100.0f*pct_sse_decrease, nreductions, dt) ;
        MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
        MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
        MRISrestoreVertexPositions(mris, INFLATED_VERTICES) ;
        MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
        if (gMRISexternalTimestep)
          (*gMRISexternalTimestep)(mris, parms) ;
        sse = last_sse ; mle_sse = last_mle_sse ;
				/*        nreductions = MAX_REDUCTIONS+1 ;*/
        if (ripped)
          break ;
        n-- ; /* don't count this as a time step */
      }
      else
      {
        printf("error decreased by %2.4f%% - %dth time step reduction: dt=%2.3f\n",
               100.0f*pct_sse_decrease, nreductions, dt) ;
				increased = 0 ;
      }
      if ((nreductions > MAX_REDUCTIONS))
      {
#if 0
        if (ripped == 0)
        {
          nreductions = 0 ;
          dt = parms->dt ;
          ripped = 1 ;
          nreductions = 0 ;
          printf("****** ripping vertices that have asymptoted *****\n") ;
          if (gMRISexternalRipVertices)
            (*gMRISexternalRipVertices)(mris, parms) ;
          continue ;
        }
#endif
        n++ ; break ;
      }
    }
		else
		{
			last_mle_sse = mle_sse ;
			increased = 0 ;
		}


    if (parms->flags & IPFLAG_ADD_VERTICES)
    {
      float max_len ;
      
      MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
      for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
        while (MRISdivideLongEdges(mris, max_len) > 0)
        {}

      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
      for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
        while (MRISdivideLongEdges(mris, max_len) > 0)
        {}

      if (gMRISexternalTimestep)
        (*gMRISexternalTimestep)(mris, parms) ;
    }

    /* recompute sse after external timestep, since outward and inward
       distances will have changed */
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
#if 0
    wm_sse = MRIScomputeSSE(mris, parms) ;
    MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    pial_sse = MRIScomputeSSE(mris, parms) ;
    sse = last_sse = wm_sse + pial_sse ;
#endif
#if 0
    rms = (*gMRISexternalRMS)(mris, parms) ;
#else
    rms = sqrt(mle_sse/(float)mris->nvertices) ;
#endif
    if (Gdiag & DIAG_SHOW)
      printf("%3.3d: dt: %2.4f, sse=%2.2f, pial sse=%2.2f, wm sse=%2.2f, rms=%2.2f, (%2.2f%%)\n",
              n+1,(float)delta_t, (float)sse/(float)mris->nvertices, 
              (float)pial_sse/(float)mris->nvertices, 
              (float)wm_sse/(float)mris->nvertices, (float)rms,
             100.0f*pct_sse_decrease);

    if (Gdiag & DIAG_WRITE)
    {
      fprintf(parms->fp, "%3.3d: dt: %2.4f, sse=%2.2f, pial sse=%2.2f, "
              "wm sse=%2.2f, rms=%2.2f (%2.2f%%)\n",
              n+1,(float)delta_t, (float)sse/(float)mris->nvertices, 
              (float)pial_sse/(float)mris->nvertices, 
              (float)wm_sse/(float)mris->nvertices, (float)rms,
              100.0f*pct_sse_decrease) ;

      fflush(parms->fp) ;
    }
    if ((parms->write_iterations > 0) &&
        !((n+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
      mrisWriteSnapshots(mris, parms, n+1) ;

    if ((Gdiag & DIAG_SHOW) && !((n+1)%5))
      MRISprintTessellationStats(mris, stderr) ;
    if ((nreductions > MAX_REDUCTIONS) && ripped)
    {
      n++ ;  /* count this step */
      break ;
    }
  }

  MRISunrip(mris) ;

  parms->start_t = n ; parms->dt = base_dt ;
  if (Gdiag & DIAG_SHOW)
  {
    msec = TimerStop(&then) ;
    fprintf(stdout,"positioning took %2.1f minutes\n", 
            (float)msec/(60*1000.0f));
  }
  if (Gdiag & DIAG_WRITE)
  {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }

  /*  MHTcheckSurface(mris, mht) ;*/
  if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
    MHTfree(&mht) ;
  return(NO_ERROR) ;
}


int
MRISpositionSurface(MRI_SURFACE *mris, MRI *mri_brain, MRI *mri_smooth,
                    INTEGRATION_PARMS *parms)
{
  /*  char   *cp ;*/
  int    avgs, niterations, n, write_iterations, nreductions = 0, done ;
  double sse, delta_t = 0.0, rms, dt, l_intensity, base_dt, last_sse,last_rms;
  MHT    *mht = NULL, *mht_v_orig = NULL, *mht_v_current = NULL ;
  struct timeb  then ; int msec ;

  if (!FZERO(parms->l_surf_repulse))
    mht_v_orig = MHTfillVertexTable(mris, NULL, ORIGINAL_VERTICES) ;
    
  base_dt = parms->dt ;
  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  TimerStart(&then) ;
  parms->mri_brain = mri_brain ;
  parms->mri_smooth = mri_smooth ;
  niterations = parms->niterations ;
  write_iterations = parms->write_iterations ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;

    if (!parms->fp)
    {
      sprintf(fname, "%s.%s.out", 
              mris->hemisphere==RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      else
        parms->fp = fopen(fname, "a") ;
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
      Progname, fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  mrisClearMomentum(mris) ;
  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;

  MRIScomputeNormals(mris) ;
  mrisClearDistances(mris) ;

	MRISclearCurvature(mris) ;  /* curvature will be used to calculate sulc */

  /* write out initial surface */
  if ((parms->write_iterations > 0) && (Gdiag&DIAG_WRITE) && !parms->start_t)
    mrisWriteSnapshot(mris, parms, 0) ;

  avgs = parms->n_averages ;
  last_rms = rms = mrisRmsValError(mris, mri_brain) ;
  last_sse = sse = MRIScomputeSSE(mris, parms) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "%3.3d: dt: %2.4f, sse=%2.1f, rms=%2.2f\n", 
            0, 0.0f, (float)sse, (float)rms);

  if (Gdiag & DIAG_WRITE)
  {
    fprintf(parms->fp, "%3.3d: dt: %2.4f, sse=%2.1f, rms=%2.2f\n", 
            0, 0.0f, (float)sse, (float)rms);
    fflush(parms->fp) ;
  }

  dt = parms->dt ; l_intensity = parms->l_intensity ;
  for (n = parms->start_t ; n < parms->start_t+niterations ; n++)
  {
#if 0
    if (n == parms->start_t+niterations-5)
    {
      dt = parms->dt / 4.0f ;  /* take some small steps at the end */
      /*      l_intensity = 0.0f ;*/
    }
    if (n == parms->start_t+niterations-2)
    {
      dt = dt / 4.0f ;        /* take some really small steps at the end */
      /*      l_intensity = 0.0f ;*/
    }
#endif
    if (!FZERO(parms->l_repulse))
      mht_v_current = MHTfillVertexTable(mris, mht_v_current,CURRENT_VERTICES);
    if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
      mht = MHTfillTable(mris, mht) ;
    MRISclearGradient(mris) ;
    mrisComputeIntensityTerm(mris, l_intensity, mri_brain, mri_smooth,
                             parms->sigma);
		mrisComputeShrinkwrapTerm(mris, mri_brain, parms->l_shrinkwrap) ;
    mrisComputeIntensityGradientTerm(mris, parms->l_grad,mri_brain,mri_smooth);
    mrisComputeSurfaceRepulsionTerm(mris, parms->l_surf_repulse, mht_v_orig);
    if (gMRISexternalGradient)
      (*gMRISexternalGradient)(mris, parms) ;

		/*		mrisMarkSulcalVertices(mris, parms) ;*/
#if 1
    mrisAverageSignedGradients(mris, avgs) ;
#else
    mrisAverageWeightedGradients(mris, avgs) ;
#endif
		/*		mrisUpdateSulcalGradients(mris, parms) ;*/

    /* smoothness terms */
    mrisComputeSpringTerm(mris, parms->l_spring) ;
    mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm) ;
    mrisComputeRepulsiveTerm(mris, parms->l_repulse, mht_v_current) ;
    mrisComputeThicknessSmoothnessTerm(mris, parms->l_tsmooth) ;
    mrisComputeNormalSpringTerm(mris, parms->l_nspring) ;
    mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv) ;
    /*    mrisComputeAverageNormalTerm(mris, avgs, parms->l_nspring) ;*/
    /*    mrisComputeCurvatureTerm(mris, parms->l_curv) ;*/
    mrisComputeTangentialSpringTerm(mris, parms->l_tspring) ;
    
#if 0
    switch (parms->integration_type)
    {
    case INTEGRATE_LM_SEARCH:
      delta_t = mrisLineMinimizeSearch(mris, parms) ;
      break ;
    default:
    case INTEGRATE_LINE_MINIMIZE:
      delta_t = mrisLineMinimize(mris, parms) ;
      break ;
    case INTEGRATE_MOMENTUM:
      delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, 
                                     parms->tol, avgs) ;
      break ;
    case INTEGRATE_ADAPTIVE:
      mrisAdaptiveTimeStep(mris, parms);
      break ;
    }
#else
    do
    {
      MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
      delta_t = mrisAsynchronousTimeStep(mris, parms->momentum, dt,mht,
                                       MAX_ASYNCH_MM) ;
      if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
        MHTcheckFaces(mris, mht) ;
      MRIScomputeMetricProperties(mris) ; 
      rms = mrisRmsValError(mris, mri_brain) ;
      sse = MRIScomputeSSE(mris, parms) ;
      done = 1 ;
      if (rms > last_rms-0.05)  /* error increased - reduce step size */
      {
        nreductions++ ;
        parms->dt *= REDUCTION_PCT ; dt = parms->dt ;
        fprintf(stdout, 
                "rms = %2.2f, time step reduction %d of %d to %2.3f...\n",
                rms, nreductions, MAX_REDUCTIONS+1, dt) ;
        mrisClearMomentum(mris) ;
        if (rms > last_rms)  /* error increased - reject step */
        {
          MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
          MRIScomputeMetricProperties(mris) ; 

          /* if error increased and we've only reduced the time
             step a few times, try taking a smaller step (done=0).
          */
          done = (nreductions > MAX_REDUCTIONS) ;
        }
      }
    } while (!done) ;
    last_sse = sse ; last_rms = rms ;
#endif
		mrisTrackTotalDistanceNew(mris) ;  /* computes signed  deformation amount */
    rms = mrisRmsValError(mris, mri_brain) ;
    sse = MRIScomputeSSE(mris, parms) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "%3.3d: dt: %2.4f, sse=%2.1f, rms=%2.2f\n", 
              n+1,(float)delta_t, (float)sse, (float)rms);

    if (Gdiag & DIAG_WRITE)
    {
      fprintf(parms->fp, "%3.3d: dt: %2.4f, sse=%2.1f, rms=%2.2f\n", 
              n+1,(float)delta_t, (float)sse, (float)rms);
      fflush(parms->fp) ;
    }
 
    if ((parms->write_iterations > 0) &&
        !((n+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
      mrisWriteSnapshot(mris, parms, n+1) ;

    if ((Gdiag & DIAG_SHOW) && !((n+1)%5))
      MRISprintTessellationStats(mris, stderr) ;
    if (nreductions > MAX_REDUCTIONS)
    {
      n++ ;  /* count this step */
      break ;
    }
  }

  parms->start_t = n ; parms->dt = base_dt ;
  if (Gdiag & DIAG_SHOW)
  {
    msec = TimerStop(&then) ;
    fprintf(stdout,"positioning took %2.1f minutes\n", 
            (float)msec/(60*1000.0f));
  }
  if (Gdiag & DIAG_WRITE)
  {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }

  /*  MHTcheckSurface(mris, mht) ;*/
  if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
    MHTfree(&mht) ;
  if (mht_v_current)
    MHTfree(&mht_v_current) ;
  if (mht_v_orig)
    MHTfree(&mht_v_orig) ;
  return(NO_ERROR) ;
}
int
MRISmoveSurface(MRI_SURFACE *mris, MRI *mri_brain, MRI  *mri_smooth,
                INTEGRATION_PARMS *parms)
{
  /*  char   *cp ;*/
  double sse_before, sse_after, rms_before, rms_after ;
  MHT    *mht = NULL ;
  int     vno ;
  VERTEX  *v ;

  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  parms->mri_brain = mri_brain ;
  parms->mri_smooth = mri_smooth ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;

    if (!parms->fp)
    {
      sprintf(fname, "%s.%s.out", 
              mris->hemisphere==RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      else
        parms->fp = fopen(fname, "a") ;
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
                  Progname, fname) ;
    }
  }

  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;

  MRIScomputeNormals(mris) ;

  rms_before = mrisRmsValError(mris, mri_brain) ;
  sse_before = MRIScomputeSSE(mris, parms) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "before expansion, sse = %2.3f, rms = %2.3f\n",
            (float)sse_before, (float)rms_before) ;

  if (Gdiag & DIAG_WRITE)
  {
    /* write out initial surface */
    if (parms->write_iterations > 0)
    {
      fprintf(stdout, "writing out pre expansion surface.\n") ;
      MRISwrite(mris, "pre") ;
    }
    fprintf(parms->fp, "before expansion, sse = %2.1f, rms = %2.1f\n",
            (float)sse_before, (float)rms_before) ;
    fflush(parms->fp) ;
  }

  if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
    mht = MHTfillTable(mris, mht) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->dx = v->nx * v->d ;
    v->dy = v->ny * v->d ;
    v->dz = v->nz * v->d ;
  }
  mrisAsynchronousTimeStep(mris, 0.0, 1.0,mht, 3.0f) ;
  MRIScomputeMetricProperties(mris) ; 
  rms_after = mrisRmsValError(mris, mri_brain) ;
  sse_after = MRIScomputeSSE(mris, parms) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "after expansion, sse = %2.1f, rms = %2.1f\n",
            (float)sse_after, (float)rms_after) ;

  if (Gdiag & DIAG_WRITE)
  {
    if (parms->write_iterations > 0)
    {
      fprintf(stdout, "writing post expansion surface...\n") ;
      MRISwrite(mris, "post") ;
    }
    fprintf(parms->fp, "after expansion, sse = %2.3f, rms = %2.3f\n",
            (float)sse_after, (float)rms_after) ;
    fflush(parms->fp) ;
  }

  if (Gdiag & DIAG_WRITE)
  {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }

  /*  MHTcheckSurface(mris, mht) ;*/
  if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
    MHTfree(&mht) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRISwriteSurfaceIntoVolume(MRI_SURFACE *mris, MRI *mri_template, MRI *mri)
{
  int fno ;

  if (!mri)
  {
    mri = MRIalloc(256*2, 256*2, 256*2, MRI_BITMAP) ;
    MRIcopyHeader(mri_template, mri) ;
    MRIsetResolution(mri, 0.5, 0.5, 0.5) ;
    mri->xstart = mri_template->xstart ; mri->xend = mri_template->xend ;
    mri->ystart = mri_template->ystart ; mri->yend = mri_template->yend ;
    mri->zstart = mri_template->zstart ; mri->zend = mri_template->zend ;
  }
  else
    MRIclear(mri) ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
    mrisFillFace(mris, mri, fno) ;
    
  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Fill in a triangle to prevent self-intersection at SAMPLE_DIST
          intervals (should be 1/2 resolution of mri volume).

       V0    b     V2
        o----------o
        |        /       
        |      /         
      a |    /            
        |  /             
        |/      
        o
       V1      b        V2        
------------------------------------------------------*/
static int
mrisFillFace(MRI_SURFACE *mris, MRI *mri, int fno)
{
  return(mrisHatchFace(mris, mri, fno, 1)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        Description
        each face has 2 triangles defined by it:

       V0    b     V2
        o----------o
        |        /       
        |      /         
      a |    /            
        |  /             
        |/      
        o
       V1      b        V2        
------------------------------------------------------*/
#define SAMPLE_DIST   0.25

static int
mrisHatchFace(MRI_SURFACE *mris, MRI *mri, int fno, int on)
{
  Real   x, y, z, xa, ya, za, xc, yc, zc, t0, t1, adx, ady, adz, dx, dy, dz, 
         cdx, cdy, cdz, alen, clen, delta_t0, delta_t1, len ;
  int    xv, yv, zv, i ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *face ;

  face = &mris->faces[fno] ;
  if (face->ripflag)
    return(NO_ERROR) ;
  
  for (i = 0 ; i < 1 ; i++)
  {
    switch (i)
    {
    default:
    case 0:
      v0 = &mris->vertices[face->v[0]] ;
      v1 = &mris->vertices[face->v[1]] ;
      v2 = &mris->vertices[face->v[2]] ;
      break ;
    case 1:
      v0 = &mris->vertices[face->v[1]] ;
      v1 = &mris->vertices[face->v[2]] ;
      v2 = &mris->vertices[face->v[0]] ;
      break ;
    case 2:
      v0 = &mris->vertices[face->v[2]] ;
      v1 = &mris->vertices[face->v[0]] ;
      v2 = &mris->vertices[face->v[1]] ;
      break ;
    }

    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;
    adx = v1->x - v0->x ; ady = v1->y - v0->y ; adz = v1->z - v0->z ;
    alen = sqrt(SQR(adx)+SQR(ady)+SQR(adz)) ;
    cdx = v2->x - v0->x ; cdy = v2->y - v0->y ; cdz = v2->z - v0->z ;
    clen = sqrt(SQR(cdx)+SQR(cdy)+SQR(cdz)) ;
    
    /*
      sample along legs of the triangle making sure the maximum spacing
      between samples (along the longer leg) is SAMPLE_DIST.
    */
    
    /*
      move along v0->v1 and v3->v2 lines and draw in crossing line to fill face
      t0 parameterizes lines from v0->v1 and v0->v2 
    */
    if (FZERO(alen) && FZERO(clen))
      delta_t0 = 0.99 ;
    else
      delta_t0 = (alen > clen) ? (SAMPLE_DIST / alen) : (SAMPLE_DIST / clen ) ;
    if (FZERO(delta_t0))
      ErrorReturn(ERROR_BADPARM, 
                  (ERROR_BADPARM,
                   "mrisFillFace: face %d has infinite leg (%d, %d)\n",
                   fno, alen, clen)) ;
    
    if (delta_t0 >= 1.0)
      delta_t0 = 0.99 ;
    
    /* delta_t0 is % of alen or clen (whichever is bigger) of SAMPLE_DIST */
    for (t0 = 0 ; t0 <= 1.0f ; t0 += delta_t0)
    {
      /* compute points (xa,ya,za) and (xc,yc,zc) on the a and c lines resp. */
      xa = v0->x + t0*adx ; ya = v0->y + t0*ady ; za = v0->z + t0*adz ;
      xc = v0->x + t0*cdx ; yc = v0->y + t0*cdy ; zc = v0->z + t0*cdz ;
      dx = xc-xa ; dy = yc-ya ; dz = zc-za ;
      len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
      if (FZERO(len))
        delta_t1 = 0.99 ;
      else
      {
        delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
        if (delta_t1 >= 1.0f)
          delta_t1 = 0.99 ;
      }
      
      /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
      for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1)
      {
        /* compute a point on the line connecting a and c */
        x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
        // MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
        MRIsurfaceRASToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
        xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;  /* voxel coordinate */
        if (on)
          MRIset_bit(mri, xv, yv, zv) ;                 /* mark it filled */
        else
          MRIclear_bit(mri, xv, yv, zv) ;               /* mark it empty */
      }
      /* compute last point on line */
      t1 = 1.0f ;
      x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
      // MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
      MRIsurfaceRASToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;  /* voxel coordinate */
      if (on)
        MRIset_bit(mri, xv, yv, zv) ;                 /* mark it filled */
      else
        MRIclear_bit(mri, xv, yv, zv) ;               /* mark it empty */
    }
  
    /* compute last line on the a and c lines resp. */
    t0 = 1.0f ;
    xa = v0->x + t0*adx ; ya = v0->y + t0*ady ; za = v0->z + t0*adz ;
    xc = v0->x + t0*cdx ; yc = v0->y + t0*cdy ; zc = v0->z + t0*cdz ;
    dx = xc-xa ; dy = yc-ya ; dz = zc-za ;
    len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
    if (FZERO(len))
      delta_t1 = 0.99 ;
    else
    {
      delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
      if (delta_t1 >= 1.0f)
        delta_t1 = 0.99 ;
    }

    /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
    for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1)
    {
      /* compute a point on the line connecting a and c */
      x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
      // MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
      MRIsurfaceRASToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;  /* voxel coordinate */
      if (on)
        MRIset_bit(mri, xv, yv, zv) ;                 /* mark it filled */
      else
        MRIclear_bit(mri, xv, yv, zv) ;               /* mark it empty */
    }
    /* compute last point on line */
    t1 = 1.0f ;
    x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
    // MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
    MRIsurfaceRASToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;  /* voxel coordinate */
    if (on)
      MRIset_bit(mri, xv, yv, zv) ;                 /* mark it filled */
    else
      MRIclear_bit(mri, xv, yv, zv) ;               /* mark it empty */
  }

  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisEraseFace(MRI_SURFACE *mris, MRI *mri, int fno)
{
  return(mrisHatchFace(mris, mri, fno, 0)) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the cortical thickness at each vertex by measuring the
          distance from it to the pial surface. 

          This routine assumes that the white matter surface is stored in
          ORIGINAL_VERTICES, and that the current vertex positions reflect
          the pial surface.
------------------------------------------------------*/

#if 1
int
MRISfindClosestOrigVertices(MRI_SURFACE *mris, int nbhd_size)
{
  int     vno, n, vlist[100000], vtotal, ns, i, vnum, nbr_count[100], min_n, min_vno ;
  VERTEX  *v, *vn, *vn2 ;
  float   dx, dy, dz, dist, min_dist, nx, ny, nz, dot ;

  memset(nbr_count, 0, 100*sizeof(int)) ;

  /* current vertex positions are gray matter, orig are white matter */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
    min_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    v->marked = 1 ; vtotal = 1 ; vlist[0] = vno ;
    min_n = 0 ; min_vno = vno ;
    for (ns = 1 ; ns <= nbhd_size ; ns++)
    {
      vnum = 0 ;  /* will be # of new neighbors added to list */
      for (i = 0 ; i < vtotal ; i++)
      {
        vn = &mris->vertices[vlist[i]] ;
        if (vn->ripflag)
          continue ;
        if (vn->marked && vn->marked < ns-1)
          continue ;
        for (n = 0 ; n < vn->vnum ; n++)
        {
          vn2 = &mris->vertices[vn->v[n]] ;
          if (vn2->ripflag || vn2->marked)  /* already processed */
            continue ;
          vlist[vtotal+vnum++] = vn->v[n] ;
          vn2->marked = ns ;
          dx = vn2->x-v->origx ; dy = vn2->y-v->origy ; dz = vn2->z-v->origz ;
          dot = dx*nx + dy*ny + dz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dot = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist)
          {
            min_n = ns ;
            min_dist = dist ;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON)
              fprintf(stdout, "%d --> %d = %2.3f\n",
                      vno,vn->v[n], dist) ;
						min_vno = vn->v[n] ;
          }
        }
      }
      vtotal += vnum ;
    }

    nbr_count[min_n]++ ;
    for (n = 0 ; n < vtotal ; n++)
    {
      vn = &mris->vertices[vlist[n]] ;
      if (vn->ripflag)
        continue ;
      vn->marked = 0 ;
    }
		v->curv = min_vno ;
  }


  for (n = 0 ; n <= nbhd_size ; n++)
    fprintf(stdout, "%d vertices at %d distance\n", nbr_count[n], n) ;
  return(NO_ERROR) ;
}
int
MRISmeasureCorticalThickness(MRI_SURFACE *mris, int nbhd_size, float max_thick)
{
  int     vno, n, vlist[100000], vtotal, ns, i, vnum, nbr_count[100], min_n,
          nwg_bad, ngw_bad ;
  VERTEX  *v, *vn, *vn2 ;
  float   dx, dy, dz, dist, min_dist, nx, ny, nz, dot ;

  memset(nbr_count, 0, 100*sizeof(int)) ;
  nwg_bad = ngw_bad = 0 ;

  /* current vertex positions are gray matter, orig are white matter */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (!(vno % 25000))
      fprintf(stdout, "%d of %d vertices processed\n", vno,mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
    min_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    v->marked = 1 ; vtotal = 1 ; vlist[0] = vno ;
    min_n = 0 ;
    for (ns = 1 ; ns <= nbhd_size ; ns++)
    {
      vnum = 0 ;  /* will be # of new neighbors added to list */
      for (i = 0 ; i < vtotal ; i++)
      {
        vn = &mris->vertices[vlist[i]] ;
        if (vn->ripflag)
          continue ;
        if (vn->marked && vn->marked < ns-1)
          continue ;
        for (n = 0 ; n < vn->vnum ; n++)
        {
          vn2 = &mris->vertices[vn->v[n]] ;
          if (vn2->ripflag || vn2->marked)  /* already processed */
            continue ;
          vlist[vtotal+vnum++] = vn->v[n] ;
          vn2->marked = ns ;
          dx = vn2->x-v->origx ; dy = vn2->y-v->origy ; dz = vn2->z-v->origz ;
          dot = dx*nx + dy*ny + dz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dot = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist)
          {
            min_n = ns ;
            min_dist = dist ;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON)
              fprintf(stdout, "%d --> %d = %2.3f\n",
                      vno,vn->v[n], dist) ;
          }
        }
      }
      vtotal += vnum ;
    }

    nbr_count[min_n]++ ;
    for (n = 0 ; n < vtotal ; n++)
    {
      vn = &mris->vertices[vlist[n]] ;
      if (vn->ripflag)
        continue ;
      vn->marked = 0 ;
    }
    if (min_dist > max_thick)
    {
      nwg_bad++ ;
      min_dist = max_thick ;
    }
    v->curv = min_dist ;
  }


  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (!(vno % 25000))
      fprintf(stdout, "%d of %d vertices processed\n", vno,mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
    min_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    v->marked = 1 ; vtotal = 1 ; vlist[0] = vno ;
    min_n = 0 ;
    for (ns = 1 ; ns <= nbhd_size ; ns++)
    {
      vnum = 0 ;  /* will be # of new neighbors added to list */
      for (i = 0 ; i < vtotal ; i++)
      {
        vn = &mris->vertices[vlist[i]] ;
        if (vn->ripflag)
          continue ;
        if (vn->marked && vn->marked < ns-1)
          continue ;
        for (n = 0 ; n < vn->vnum ; n++)
        {
          vn2 = &mris->vertices[vn->v[n]] ;
          if (vn2->ripflag || vn2->marked)  /* already processed */
            continue ;
          vlist[vtotal+vnum++] = vn->v[n] ;
          vn2->marked = ns ;
          dx = v->x-vn2->origx ; dy = v->y-vn2->origy ; dz = v->z-vn2->origz ;
          dot = dx*nx + dy*ny + dz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dot = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist)
          {
            min_n = ns ;
            min_dist = dist ;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON)
              fprintf(stdout, "%d --> %d = %2.3f\n",
                      vno,vn->v[n], dist) ;
          }
        }
      }
      vtotal += vnum ;
    }

    nbr_count[min_n]++ ;
    for (n = 0 ; n < vtotal ; n++)
    {
      vn = &mris->vertices[vlist[n]] ;
      if (vn->ripflag)
        continue ;
      vn->marked = 0 ;
    }
    if (DIAG_VERBOSE_ON && fabs(v->curv - min_dist) > 4.0)
      fprintf(stdout, "v %d, white->gray=%2.2f, gray->white=%2.2f\n",
              vno, v->curv, min_dist) ;
    if (min_dist > max_thick)
    {
      min_dist = max_thick ;
      ngw_bad++ ;
    }
    v->curv = (v->curv+min_dist)/2 ;
  }

  fprintf(stdout, "thickness calculation complete, %d:%d truncations.\n",
          nwg_bad, ngw_bad) ;
  for (n = 0 ; n <= nbhd_size ; n++)
    fprintf(stdout, "%d vertices at %d distance\n", nbr_count[n], n) ;
  return(NO_ERROR) ;
}
#else
#define MAX_THICKNESS 6.0f
int
MRISmeasureCorticalThickness(MRI_SURFACE *mris)
{
  int     vno ;
  VERTEX  *v ;
  float   max_out_dist ;
  MHT     *mht ;

  mht = MHTfillTable(mris, NULL) ;  /* fill table with pial surface position */
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;

  /* compute white matter vertex normals */
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeNormals(mris);
  MRISsmoothSurfaceNormals(mris, 10) ;

  /* must restore gray matter surface so self-intersection will work */
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ; 

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (!(vno%50))
      DiagHeartbeat((float)vno / (float)(mris->nvertices-1)) ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    max_out_dist = 
      mrisFindNormalDistance(mris, mht, vno, 1.5*MAX_THICKNESS);
    if (max_out_dist > MAX_THICKNESS)   /* can't compute it properly */
    {
      float  dx, dy, dz ;
      dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
      v->curv = sqrt(dx*dx+dy*dy+dz*dz) ;
      v->marked = 1 ;
    }
    else
      v->curv = max_out_dist ;
  }

  MHTfree(&mht) ;
  return(NO_ERROR) ;
}
#endif
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisFindNextOutwardFace(MRI_SURFACE *mris, MHT *mht, int vno, double max_dist)
{
  int     nfound, flist[1000], *fptr, total_found, i, j ;
  double  dist, d ;
  VERTEX  *v ;

  v = &mris->vertices[vno] ;

  for (total_found = nfound = 0, dist = 0.0 ; dist <= max_dist ; dist += .25)
  {
    d = dist ; fptr = &flist[total_found] ;
    nfound = 
      mrisAllNormalDirectionCurrentTriangleIntersections(mris, v,mht,&d,fptr);
    if (nfound > 0)
    {
      for (i = 0 ; i < total_found ; i++)
      {
        for (j = 0 ; j < nfound ; j++)
        {
          if (flist[i] == fptr[j])
            fptr[j] = -1 ;   /* was already found */
        }
      }
      for (j = 0 ; j < nfound ; j++)
      {
        if (fptr[j] >= 0)
          flist[total_found++] = fptr[j] ;
      }
      for (i = 0 ; i < total_found ; i++)
        if (vertexInFace(mris, vno, flist[i]))
          flist[i] = -1 ;
      nfound = total_found ;
      for (fptr = flist,total_found = 0, j = 0 ; j < nfound ; j++, fptr++)
      {
        if (*fptr >= 0)
          flist[total_found++] = *fptr ;
      }
    }
    if (total_found > 0)
    {
      if (vno == Gdiag_no)
      {
        fprintf(stdout, "v %d @ (%2.2f, %2.2f, %2.2f), f ", 
                vno, v->x,v->y,v->z);
        for (i = 0 ; i < v->num ; i++)
          fprintf(stdout, "%d ", v->f[i]) ;
        fprintf(stdout, "\n") ;
        for (i = 0 ; i < total_found ; i++)
        {
          FACE *f = &mris->faces[flist[i]] ;
          fprintf(stdout, "\tface %d with vertices (%d, %d, %d)\n",
                  flist[i], f->v[0], f->v[1], f->v[2]) ;
        }
      }
      return(flist[0]) ;
    }
  }
  return(-1) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisFindNextInwardFace(MRI_SURFACE *mris, MHT *mht, int vno, double max_dist)
{
  int     nfound, flist[1000], *fptr, total_found, i, j ;
  double  dist, d ;
  VERTEX  *v ;

  v = &mris->vertices[vno] ;
  v->nx *= -1 ; v->ny *= -1 ; v->nz *= -1 ; 
  for (total_found = nfound = 0, dist = 0.0 ; dist <= max_dist ; dist += .25)
  {
    d = dist ; fptr = &flist[total_found] ;
    nfound = 
      mrisAllNormalDirectionCurrentTriangleIntersections(mris, v,mht,&d,fptr);
    if (nfound > 0)
    {
      for (i = 0 ; i < total_found ; i++)
      {
        for (j = 0 ; j < nfound ; j++)
        {
          if (flist[i] == fptr[j])
            fptr[j] = -1 ;   /* was already found */
        }
      }
      for (j = 0 ; j < nfound ; j++)
      {
        if (fptr[j] >= 0)
          flist[total_found++] = fptr[j] ;
      }
      for (i = 0 ; i < total_found ; i++)
        if (vertexInFace(mris, vno, flist[i]))
          flist[i] = -1 ;
      nfound = total_found ;
      for (fptr = flist,total_found = 0, j = 0 ; j < nfound ; j++, fptr++)
      {
        if (*fptr >= 0)
          flist[total_found++] = *fptr ;
      }
    }
    if (total_found > 0)
    {
      v->nx *= -1 ; v->ny *= -1 ; v->nz *= -1 ; 
      return(flist[0]) ;
    }
  }
  v->nx *= -1 ; v->ny *= -1 ; v->nz *= -1 ; 
  return(-1) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisFindNormalDistance(MRI_SURFACE *mris, MHT *mht, int vno, double max_dist)
{
  double   dist ;
  VERTEX   *v ;

  v = &mris->vertices[vno] ;
  if (v->ripflag)
    return(0.0) ;

  for (dist = 0.0f ; dist < max_dist ; dist += .25)
  {
    if (mrisNormalDirectionTriangleIntersection(mris, v, mht, &dist,NULL) > 0)
      return(dist) ;
  }
  
  return(dist) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           See if the line in the normal direction passing
           through vertex v intersects any of the triangles at the
           given location.
------------------------------------------------------*/
static int
mrisDirectionTriangleIntersection(MRI_SURFACE *mris, float x0, float y0, 
                                  float z0, float nx, float ny,
                                  float nz, MHT *mht, double *pdist)
{
  double  dist, min_dist, U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3], dot ;
  float   x, y, z, dx, dy, dz ;
  MHBT    *bucket ;
  MHB     *bin ;
  int     i, found, fno, ret ;
  static MHBT *last_bucket = NULL ;
  static float lastx, lasty, lastz = -1 ;

  dist = *pdist ;
  dir[0] = nx ; dir[1] = ny ; dir[2] = nz ;
  pt[0] = x0  ; pt[1] = y0  ; pt[2] = z0  ;
  x = x0 + nx * dist ;
  y = y0 + ny * dist ; 
  z = z0 + nz * dist ;

  min_dist = 10000.0f ;
#if 1
  bucket = MHTgetBucket(mht, x, y, z) ;
  if (bucket == NULL)
    return(0) ;

#if 0
  if (lastx == x0 && lasty == y0 && lastz == z0 && bucket == last_bucket)
    return(0) ;
#endif

  lastx = x0 ; lasty = y0 ; lastz = z0 ; last_bucket = bucket ;

  for (bin = bucket->bins, found = i = 0 ; i < bucket->nused ; i++, bin++)
  {
    fno = bin->fno ;
    if (fno == 1287 || fno == 5038)
      DiagBreak() ;
#else
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
#endif
    load_triangle_vertices(mris, fno, U0, U1, U2) ;
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt) ;
    if (ret)
    {
      dx = int_pt[0] - x0 ; 
      dy = int_pt[1] - y0 ; 
      dz = int_pt[2] - z0 ; 
      dist = sqrt(dx*dx + dy*dy + dz*dz) ;
      dot = dx*nx + dy*ny + dz*nz ;
      if (dot >= 0 && dist < min_dist)
      {
        found = 1 ;
        *pdist = min_dist = dist ;
      }
    }
  }
  return(found) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           See if the line in the normal direction passing
           through vertex v intersects any of the triangles at the
           given location.
------------------------------------------------------*/
#if 0
static int
mrisNormalDirectionTriangleIntersection(MRI_SURFACE *mris, VERTEX *v, 
                                        MHT *mht, double *pdist, int *flist)
{
  double  dist, min_dist, U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3], dot ;
  float   nx, ny, nz, x, y, z, dx, dy, dz ;
  MHBT    *bucket ;
  MHB     *bin ;
  int     i, found, fno, ret ;
  static MHBT *last_bucket = NULL ;
  static VERTEX *last_v = NULL ;


  dist = *pdist ;
  nx = v->nx ; ny = v->ny ; nz = v->nz ; 
  dir[0] = v->nx ; dir[1] = v->ny ; dir[2] = v->nz ;
  pt[0] = v->origx  ; pt[1] = v->origy  ; pt[2] = v->origz  ;
  x = v->origx + nx * dist ;
  y = v->origy + ny * dist ; 
  z = v->origz + nz * dist ;

  bucket = MHTgetBucket(mht, x, y, z) ;
  if (bucket == NULL)
    return(-1) ;

#if 0
  v->origx = v->origy = .5 ; v->origz = 1 ;
  v->nx = v->ny = 0 ; v->nz = 1 ;
  nx = v->nx ; ny = v->ny ; nz = v->nz ; 
  dir[0] = v->nx ; dir[1] = v->ny ; dir[2] = v->nz ;
  pt[0] = v->origx  ; pt[1] = v->origy  ; pt[2] = v->origz  ;
  x = v->origx + nx * dist ;
  y = v->origy + ny * dist ; 
  z = v->origz + nz * dist ;
#endif

  if (last_v == v && bucket == last_bucket)
    return(-2) ;
  last_v = v ; last_bucket = bucket ;

  min_dist = 10000.0f ;
  for (bin = bucket->bins, found = i = 0 ; i < bucket->nused ; i++, bin++)
  {
    fno = bin->fno ;

#if 0
    v->origx = v->origy = .5 ; v->origz = 1 ;
    {
      FACE   *f = &mris->faces[fno] ;
      VERTEX *v ;
      v = &mris->vertices[f->v[0]] ;
      v->x = v->y = v->z = 0 ;
      v = &mris->vertices[f->v[1]] ;
      v->x = 1 ; v->y = v->z = 0 ;
      v = &mris->vertices[f->v[2]] ;
      v->x = 0 ; v->y = 1 ; v->z = 0 ;
    }
#endif

    load_triangle_vertices(mris, fno, U0, U1, U2) ;
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt) ;
    if (ret)
    {
      dx = int_pt[0] - v->origx ; 
      dy = int_pt[1] - v->origy ; 
      dz = int_pt[2] - v->origz ; 
      dist = sqrt(dx*dx + dy*dy + dz*dz) ;
      dot = dx*nx + dy*ny + dz*nz ;
      if (dot >= 0 && dist < min_dist)
      {
        if (flist)
          flist[found] = fno ;
        found++ ;
        *pdist = min_dist = dist ;
      }
    }
  }
  return(found) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           See if the line in the normal direction passing
           through vertex v intersects any of the triangles at the
           given location.
------------------------------------------------------*/
static int
mrisAllNormalDirectionCurrentTriangleIntersections(MRI_SURFACE *mris, 
                                                   VERTEX *v, MHT *mht, 
                                                   double *pdist, int *flist)
{
  double  dist, min_dist, U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3] ;
  float   nx, ny, nz, x, y, z, dx, dy, dz, dot ;
  MHBT    *bucket ;
  MHB     *bin ;
  int     i, found, fno, ret ;
  static MHBT *last_bucket = NULL ;
  static VERTEX *last_v = NULL ;


  dist = *pdist ;
  nx = v->nx ; ny = v->ny ; nz = v->nz ; 
  dir[0] = v->nx ; dir[1] = v->ny ; dir[2] = v->nz ;
  pt[0] = v->x  ; pt[1] = v->y  ; pt[2] = v->z  ;
  x = v->x + nx * dist ;
  y = v->y + ny * dist ; 
  z = v->z + nz * dist ;

  bucket = MHTgetBucket(mht, x, y, z) ;
  if (bucket == NULL)
    return(-1) ;

  if (last_v == v && bucket == last_bucket)
    return(-2) ;
  last_v = v ; last_bucket = bucket ;

  min_dist = 10000.0f ;
  for (bin = bucket->bins, found = i = 0 ; i < bucket->nused ; i++, bin++)
  {
    fno = bin->fno ;

    load_triangle_vertices(mris, fno, U0, U1, U2) ;
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt) ;
    if (ret)
    {
      dx = int_pt[0] - v->x ; 
      dy = int_pt[1] - v->y ; 
      dz = int_pt[2] - v->z ; 
      dot = dx*nx + dy*ny + dz*nz ;
      if (dot < 0)   /* in direciton antiparallel to normal direction */
        continue ;
      dist = sqrt(dot) ;
      if (flist)
        flist[found] = fno ;
      found++ ;
      if (dist < min_dist)
        *pdist = min_dist = dist ;
    }
  }
  return(found) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           See if the line in the normal direction passing
           through vertex v intersects any of the triangles at the
           given location.
------------------------------------------------------*/
static int
mrisAllCurrentTriangleIntersections(MRI_SURFACE *mris, float x, float y,
                                    float z, float nx, float ny, float nz, 
                                    MHT *mht, int *flist)
{
  double  U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3] ;
  MHBT    *bucket ;
  MHB     *bin ;
  int     i, found, fno, ret ;


  dir[0] = nx ; dir[1] = ny ; dir[2] = nz ;
  pt[0] = x  ; pt[1] = y  ; pt[2] = z  ;

  bucket = MHTgetBucket(mht, x, y, z) ;
  if (bucket == NULL)
    return(-1) ;

  for (bin = bucket->bins, found = i = 0 ; i < bucket->nused ; i++, bin++)
  {
    fno = bin->fno ;

    load_triangle_vertices(mris, fno, U0, U1, U2) ;
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt) ;
    if (ret)
    {
      flist[found] = fno ;
      found++ ;
    }
  }
  return(found) ;
}
#endif

static int
load_orig_triangle_vertices(MRI_SURFACE *mris, int fno, double U0[3], 
                            double U1[3], double U2[3])
{
  VERTEX *v ;
  FACE   *face ;
  
  face = &mris->faces[fno] ;
  v = &mris->vertices[face->v[0]] ;
  U0[0] = v->origx ; U0[1] = v->origy ; U0[2] = v->origz ; 
  v = &mris->vertices[face->v[1]] ;
  U1[0] = v->origx ; U1[1] = v->origy ; U1[2] = v->origz ; 
  v = &mris->vertices[face->v[2]] ;
  U2[0] = v->origx ; U2[1] = v->origy ; U2[2] = v->origz ; 
  return(NO_ERROR) ;
}
static int
load_triangle_vertices(MRI_SURFACE *mris, int fno, double U0[3], double U1[3], 
                       double U2[3])
{
  VERTEX *v ;
  FACE   *face ;
  
  face = &mris->faces[fno] ;
  v = &mris->vertices[face->v[0]] ;
  U0[0] = v->x ; U0[1] = v->y ; U0[2] = v->z ; 
  v = &mris->vertices[face->v[1]] ;
  U1[0] = v->x ; U1[1] = v->y ; U1[2] = v->z ; 
  v = &mris->vertices[face->v[2]] ;
  U2[0] = v->x ; U2[1] = v->y ; U2[2] = v->z ; 
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisFindClosestFilledVoxel(MRI_SURFACE *mris, MRI *mri_filled, int vno, 
                           double max_dist)
{
  VERTEX  *v ;
  double  min_dist, dist, xd, yd, zd ;
  int     xoff, yoff, zoff, whalf, x, y, z, xi, yi, zi, window_size ;
  Real    xw, yw, zw ;

  v = &mris->vertices[vno] ;
  if (v->ripflag)
    return(0.0) ;
  MRISvertexToVoxel(v, mri_filled, &xw, &yw, &zw) ;
  x = nint(xw) ; y = nint(yw) ; ; z = nint(zw) ;

  window_size = nint(max_dist/mri_filled->xsize + 0.5) ;
  whalf = (window_size-1)/2 ;

  min_dist = 10.0*max_dist*max_dist ;
  for (zoff = -whalf ; zoff <= whalf ; zoff++)
  {
    zi = mri_filled->zi[zoff+z] ;
    zd = zi - z ;
    for (yoff = -whalf ; yoff <= whalf ; yoff++)
    {
      yi = mri_filled->yi[yoff+y] ;
      yd = yi-y ;
      for (xoff = -whalf ; xoff <= whalf ; xoff++)
      {
        xi = mri_filled->xi[xoff+x] ;
        if (MRItest_bit(mri_filled, xi, yi, zi))
        {
          xd = xi-x ; 
          dist = xd*xd+yd*yd+zd*zd ;
          if (dist < min_dist)
            min_dist = dist ;
        }
      }
    }
  }

  min_dist = sqrt(min_dist) * mri_filled->xsize ;

  return(min_dist) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  x, y, z, num ;
  VERTEX *v, *vn ;
  int    nmarked ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (nmarked = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 1)
        continue ;
      x = y = z = 0;
      num = 0;
      /*      if (v->marked == 2)*/
      {
        x = v->x ; y = v->y ; z = v->z ;
        num++ ;   /* account for central vertex */
      }
      pnb = v->v ;
      vnum = v->vnum ;
      for (vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
#if 0
        if (vn->ripflag || !vn->marked || vn->marked > 2) /* no valid data */
#else
        if (vn->ripflag) /* no valid data */
#endif
          continue ;
        num++ ;
        x += vn->x ; y += vn->y ; z += vn->z ;
      }
      if (num>0)
      {
        v->tdx = x / num ;
        v->tdy = y / num ;
        v->tdz = z / num ;
        if (!v->marked)
          nmarked++ ;
        v->marked = 3 ;  /* need modification */
      }
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 1)
        continue ;
      if (v->marked)
      {
        v->x = v->tdx ; v->y = v->tdy ; v->z = v->tdz ;
      }
      if (v->marked == 3)  /* needs modification */
        v->marked = 2 ;    /* modified, but not fixed */
    }
    if (nmarked && (Gdiag & DIAG_SHOW))
      printf("%d: %d vertices marked\n", i,nmarked);
#if 0
    if (!nmarked)
      break ;
#endif
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsoapBubbleOrigVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  x, y, z, num ;
  VERTEX *v, *vn ;
  int    nmarked, num_none_marked = 0 ;

  /*
    v->marked:
    
    0 -  never processed
    1 -  fixed value
    2 -  has had value computed via previous soap bubble.
    3 -  has had value computed on this soap bubble, but not yet usable.
  */
  for (i = 0 ; i < navgs ; i++)
  {
    for (nmarked = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 1)
        continue ;
      x = y = z = 0;
      num = 0;
      if (v->marked == 2)  /* computed on previous iteration, use value */
      {
        x = v->origx ; y = v->origy ; z = v->origz ;
        num++ ;   /* account for central vertex */
      }
      pnb = v->v ;
      vnum = v->vnum ;
      for (vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag || !vn->marked || vn->marked > 2) /* no valid data */
          continue ;
        num++ ;
        x += vn->origx ; y += vn->origy ; z += vn->origz ;
      }
      if (num>0)
      {
        v->tdx = x / num ;
        v->tdy = y / num ;
        v->tdz = z / num ;
        if (!v->marked)
          nmarked++ ;
        v->marked = 3 ;  /* need modification */
      }
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 1)
        continue ;
      if (v->marked)        /* update value */
      {
        v->origx = v->tdx ; v->origy = v->tdy ; v->origz = v->tdz ;
      }
      if (v->marked == 3)  /* needs modification */
        v->marked = 2 ;    /* modified, but not fixed */
    }
    if (Gdiag & DIAG_SHOW)
      printf("%d: %d vertices marked\n", i,nmarked);
    if (!nmarked && ++num_none_marked > 5)
      break ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->marked == 1)
      continue ;
    v->marked = 0 ;
  }
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs, float pct_fixed)
{
  int    i, vno, vnb, *pnb, vnum, j ;
  float  x, y, z, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    MRISclearMarks(mris) ;
    MRISmarkRandomVertices(mris, pct_fixed) ;
    for (j = 0 ; j < 10 ; j++)
    {
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag || v->marked)
          continue ;
        x = v->x ; y = v->y ; z = v->z ;
        pnb = v->v ;
        vnum = v->vnum ;
        for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
        {
          vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
          if (vn->ripflag)
            continue ;
          num++ ;
          x += vn->x ; y += vn->y ; z += vn->z ;
        }
        num++ ;   /* account for central vertex */
        v->tdx = x / num ;
        v->tdy = y / num ;
        v->tdz = z / num ;
      }
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag || v->marked)
          continue ;
        v->x = v->tdx ; v->y = v->tdy ; v->z = v->tdz ;
      }
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmarkRandomVertices(MRI_SURFACE *mris, float prob_marked)
{
  int    vno ;
  VERTEX *v ;
  float  r ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    r = randomNumber(0.0, 1.0) ;
    if (r < prob_marked)
      v->marked = 1 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISclearMarks(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->marked = 0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISclearFixedValFlags(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->fixedval = FALSE ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyFixedValFlagsToMarks(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->marked = v->fixedval ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISclearAnnotations(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->annotation = 0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsetMarks(MRI_SURFACE *mris, int mark)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->marked = mark ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsequentialAverageVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  x, y, z, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked)
        continue ;
      x = v->x ; y = v->y ; z = v->z ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        x += vn->x ; y += vn->y ; z += vn->z ;
      }
      num++ ;   /* account for central vertex */
      v->x = x / num ;
      v->y = y / num ;
      v->z = z / num ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
#define WHALF                 (5-1)/2
static int
mrisComputeWhiteSurfaceValues(MRI_SURFACE *mris, MRI *mri_brain, 
                              MRI *mri_wm, float nsigma)
{
  Real    val, x, y, z ;
  int     total_vertices, vno, xv, yv, zv, xo, yo, zo, xi, yi, zi, nvox ;
  float   total, total_sq, sigma, mean_wm, mean_gray, mean ;
  VERTEX  *v ;

  /* first compute intensity of local gray/white boundary */
  mean_wm = mean_gray = 0.0f ;
  
  for (total_vertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->d = 0.0f ;
    MRISvertexToVoxel(v, mri_wm, &x, &y, &z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    
    /* compute mean and variance of wm in a neighborhood of voxel */
    total = total_sq = 0.0f ; nvox = 0 ;
    for (zo = zv-WHALF ; zo <= zv + WHALF ; zo++)
    {
      zi = mri_wm->zi[zo] ;
      for (yo = yv-WHALF ; yo <= yv + WHALF ; yo++)
      {
        yi = mri_wm->yi[yo] ;
        for (xo = xv-WHALF ; xo <= xv + WHALF ; xo++)
        {
          xi = mri_wm->xi[xo] ;
          val = (Real)MRIvox(mri_wm, xi, yi, zi) ;
          if (val > WM_MIN_VAL)
          {
#if 0
            if (MRIneighborsOff(mri_wm, xi, yi, zi, WM_MIN_VAL) == 0)
              continue ;   /* not a border voxel */
#endif
            val = (Real)MRIvox(mri_brain, xi, yi, zi) ;  /* use smoothed val */
            total += val ;
            total_sq += val * val ;
            nvox++ ;
          }
        }
      }
    }
    if (!nvox)
      v->val = 0.0f ;
    else
    {
      mean = total / (float)nvox ;
      sigma = sqrt(total_sq / (float)nvox - mean*mean) ;
      MRISvertexToVoxel(v, mri_wm, &x, &y, &z) ;
      MRIsampleVolume(mri_brain, x, y, z, &val) ;
      v->val = mean - nsigma * sigma ;
      mean_gray += v->val ;
      mean_wm += mean ;
      total_vertices++ ;
    }
    
  }
  mean_wm /= (float)total_vertices ; mean_gray /= (float)total_vertices ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (FZERO(v->val))   /* no border voxels nearby */
      v->val = mean_gray ;
  }

#if 0
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = mean_gray ;
  }
#endif

  fprintf(stdout, "mean wm=%2.1f, gray=%2.1f, averaging targets and "
            "smoothing surface...", mean_wm, mean_gray) ;
  return(NO_ERROR) ;
}
#else
#if 1
#define MAX_CSF   55.0f
#define STEP_SIZE 0.1
int
MRIScomputeWhiteSurfaceValues(MRI_SURFACE *mris,MRI *mri_brain,MRI *mri_smooth)
{
  Real    val, x, y, z, min_val, xw, yw, zw,mag,max_mag, xw1, yw1, zw1,
          previous_val, next_val ;
  int     total_vertices, vno, nmissing = 0 ;
  float   mean_white, dist, nx, ny, nz ;
  VERTEX  *v ;

  /* first compute intensity of local gray/white boundary */
  mean_white = 0.0f ;

  MRISclearMarks(mris) ;  /* for soap bubble smoothing later */
  for (total_vertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

  /*
    search outwards and inwards and find the local gradient maximum
    at a location with a reasonable MR intensity value. This will
    be the location of the edge.
    */

    /* search in the normal direction to find the min value */
    min_val = -10.0f ; mag = 5.0f ; max_mag = 0.0f ;
    nx = v->nx ; ny = v->ny ; nz = v->nz ; 
    x = v->x ; y = v->y ; z = v->z ;
    // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    x = v->x+nx ; y = v->y + ny ; z = v->z + nz ;
    // MRIworldToVoxel(mri_brain, x, y, z, &xw1, &yw1, &zw1) ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw1, &yw1, &zw1) ;
    nx = xw1 - xw ; ny = yw1 - yw ; nz = zw1 - zw ; 
    for (dist = -3.0f ; dist < 10.0f ; dist += STEP_SIZE)
    {
      x = v->x+v->nx*(dist-1) ; y = v->y + v->ny*(dist-1) ; 
      z = v->z + v->nz*(dist-1) ;
      // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
      MRIsampleVolume(mri_brain, xw, yw, zw, &previous_val) ;
      if (previous_val < 120 && previous_val > 95)  /* in right range */
      {
        x = v->x + v->nx*dist ; y = v->y + v->ny*dist ; z = v->z + v->nz*dist ;
        // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;

        /* see if we are at a local maximum in the gradient magnitude */
        MRIsampleVolumeDerivative(mri_smooth, xw, yw, zw, nx, ny, nz, &mag) ;
        MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;

        /* if gradient is big and pointing towards wm */
        if ((previous_val > val) && (fabs(mag) > max_mag))
        {
          x = v->x + v->nx*(dist+1) ;
          y = v->y + v->ny*(dist+1) ;
          z = v->z + v->nz*(dist+1) ;
          // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
          MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
          MRIsampleVolume(mri_brain, xw, yw, zw, &next_val) ;
          if (next_val > 60 && next_val < 95)
          {          
            max_mag = fabs(mag) ;
            min_val = val ;
          }
        }
      }
    }

    if (min_val > 0)
    {
      v->val = min_val ;
      v->mean = max_mag ;
      mean_white += min_val ; total_vertices++ ;
      v->marked = 1 ;
    }
    else
      nmissing++ ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d, target value = %2.1f, mag = %2.1f\n",
              Gdiag_no, v->val, v->mean) ;
  }
  mean_white /= (float)total_vertices ; 
  MRISsoapBubbleVals(mris, 100) ; MRISclearMarks(mris) ;

  /*  MRISaverageVals(mris, 3) ;*/
  fprintf(stdout, "mean white matter surface=%2.1f, %d missing vertices\n", 
          mean_white, nmissing) ;
  return(NO_ERROR) ;
}
#else
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScomputeWhiteSurfaceValues(MRI_SURFACE *mris, MRI *mri_brain, 
                              MRI *mri_wm, float nsigma)
{
  Real    val, x, y, z ;
  int     total_vertices, vno, xv, yv, zv, xo, yo, 
          zo, xi, yi, zi, nwhite_vox, ngray_vox ;
  float   total_white, total_gray, total_sq_white, total_sq_gray, std_white,
          mean_white, mean_gray, total_mean_gray, total_mean_white, std_gray ;
  VERTEX  *v ;

  /* first compute intensity of local gray/white boundary */
  total_vertices = 0 ;
  total_mean_white = total_mean_gray = 0.0f ;
  total_sq_white = total_sq_gray = 0.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->d = 0.0f ;
    MRISvertexToVoxel(v, mri_wm, &x, &y, &z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    
    /* compute mean and variance of white in a neighborhood of voxel */
    total_white = 0.0f ; nwhite_vox = 0 ;
    total_gray = 0.0f ; ngray_vox = 0 ;
    for (zo = zv-WHALF ; zo <= zv + WHALF ; zo++)
    {
      zi = mri_wm->zi[zo] ;
      for (yo = yv-WHALF ; yo <= yv + WHALF ; yo++)
      {
        yi = mri_wm->yi[yo] ;
        for (xo = xv-WHALF ; xo <= xv + WHALF ; xo++)
        {
          xi = mri_wm->xi[xo] ;
          val = (Real)MRIvox(mri_wm, xi, yi, zi) ;
          if (val > WM_MIN_VAL)  /* white matter */
          {
            if (MRIneighborsOff(mri_wm, xi, yi, zi, WM_MIN_VAL) == 0)
              continue ;   /* not a border voxel */
            val = (Real)MRIvox(mri_brain, xi, yi, zi) ;  /* use smoothed val */
            total_white += val ;
            total_sq_white += val*val ;
            nwhite_vox++ ;
          }
          else    /* gray matter */
          {
            if (MRIneighborsOn(mri_wm, xi, yi, zi, WM_MIN_VAL+1) == 0)
              continue ;   /* not a border voxel */
            val = (Real)MRIvox(mri_brain, xi, yi, zi) ;  /* use smoothed val */
            total_gray += val ;
            total_sq_gray += val*val ;
            ngray_vox++ ;
          }
        }
      }
    }
    if (!nwhite_vox || !ngray_vox)
      v->val = 0.0f ;
    else
    {
      mean_white = total_white / (float)nwhite_vox ;
      mean_gray = total_gray / (float)ngray_vox ;
      std_white = sqrt(total_sq_white/(float)nwhite_vox-mean_white*mean_white);
      std_gray = sqrt(total_sq_gray / (float)ngray_vox-mean_gray*mean_gray);
      std_white = 0.0f ;   /* only use white matter */
      if (mean_gray > mean_white)   /* shouldn't happen */
      {
        if (DIAG_VERBOSE_ON)
          fprintf(stdout, "mean gray (%2.1f) > mean white (%2.1f) at v %d!\n",
                  mean_gray, mean_white, vno) ;
        v->val = mean_white ;
      }
      else
        v->val = (std_gray*mean_white+std_white*mean_gray)/(std_gray+std_white);
      total_vertices++ ;
      total_mean_white += mean_white ;
      total_mean_gray += mean_gray ;
    }
    v->mean = 20.0f ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d, target value = %2.1f, mag = %2.1f\n",
              Gdiag_no, v->val, v->mean) ;
  }
  total_mean_white /= (float)total_vertices ; 
  total_mean_gray /= (float)total_vertices ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (FZERO(v->val))   /* no border voxels nearby */
      v->val = (total_mean_gray+total_mean_white)/2 ;
  }

#if 0
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = mean_gray ;
  }
#endif

  fprintf(stdout, "mean white=%2.1f, gray=%2.1f\n", total_mean_white, 
          total_mean_gray) ;
  return(NO_ERROR) ;
}
#endif
#endif

int
MRIScomputeBorderValues(MRI_SURFACE *mris,MRI *mri_brain,
                        MRI *mri_smooth, Real inside_hi, Real border_hi,
                        Real border_low, Real outside_low, Real outside_hi,
                        double sigma,
                        float max_thickness, FILE *log_fp, int which)
{
  Real    val, x, y, z, max_mag_val, xw, yw, zw,mag,max_mag, max_mag_dist=0.0f,
          previous_val, next_val, min_val,inward_dist,outward_dist,xw1,yw1,zw1,
          min_val_dist, orig_dist, dx, dy, dz, previous_mag, next_mag ;
  int     total_vertices, vno, nmissing = 0, nout = 0, nin = 0, nfound = 0,
          nalways_missing = 0, local_max_found, ngrad_max, ngrad, nmin, num_changed=0 ;
  float   mean_border, mean_in, mean_out, dist, nx, ny, nz, mean_dist ;
  double  current_sigma ;
  VERTEX  *v ;
  FILE    *fp = NULL ;

	MRI     *mri_tmp ;

  mri_tmp = MRIreplaceValues(mri_brain, NULL, 255, 0) ;

  /* first compute intensity of local gray/white boundary */
  mean_dist = mean_in = mean_out = mean_border = 0.0f ;
  ngrad_max = ngrad = nmin = 0 ;
  MRISclearMarks(mris) ;  /* for soap bubble smoothing later */
  for (total_vertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ; y = v->y ; z = v->z ;
    // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    x = v->x + v->nx ; y = v->y + v->ny ; z = v->z + v->nz ;
    // MRIworldToVoxel(mri_brain, x, y, z, &xw1, &yw1, &zw1) ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw1, &yw1, &zw1) ;
    nx = xw1 - xw ; ny = yw1 - yw ; nz = zw1 - zw ; 


    /* 
       find the distance in the directions parallel and anti-parallel to
       the surface normal in which the gradient is pointing 'inwards'.
       The border will then be constrained to be within that region.
    */
    inward_dist = 1.0 ; outward_dist = -1.0 ;
    for (current_sigma = sigma; current_sigma <= 10*sigma; current_sigma *= 2)
    {
      for (dist = 0 ; dist > -max_thickness ; dist -= 0.5)
      {
        dx = v->x-v->origx ; dy = v->y-v->origy ; dz = v->z-v->origz ; 
        orig_dist = fabs(dx*v->nx + dy*v->ny + dz*v->nz) ;
        if (fabs(dist)+orig_dist > max_thickness)
          break ;
        x = v->x + v->nx*dist ; y = v->y + v->ny*dist ; z = v->z + v->nz*dist ;
        // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny,nz,&mag,
                                       current_sigma);
        if (mag >= 0.0)
          break ;
        MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;
        if (val > border_hi)
          break ;
      }
      inward_dist = dist+.25 ; 
      for (dist = 0 ; dist < max_thickness ; dist += 0.5)
      {
        dx = v->x-v->origx ; dy = v->y-v->origy ; dz = v->z-v->origz ; 
        orig_dist = fabs(dx*v->nx + dy*v->ny + dz*v->nz) ;
        if (fabs(dist)+orig_dist > max_thickness)
          break ;
        x = v->x + v->nx*dist ; y = v->y + v->ny*dist ; z = v->z + v->nz*dist ;
        // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny,nz, &mag,
                                       current_sigma);
        if (mag >= 0.0)
          break ;
        MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;
        if (val < border_low)
          break ;
      }
      outward_dist = dist-.25 ; 
      if (!finite(outward_dist))
        DiagBreak() ;
      if (inward_dist <= 0 || outward_dist >= 0)
        break ;
    }
    
    if (inward_dist > 0 && outward_dist < 0)
      current_sigma = sigma ;  /* couldn't find anything */

    if (vno == Gdiag_no)
    {
      char fname[STRLEN] ;
      sprintf(fname, "v%d.%2.0f.log", Gdiag_no, sigma*100) ;
      fp = fopen(fname, "w") ;
      fprintf(stdout, 
              "v %d: inward dist %2.2f, outward dist %2.2f, sigma %2.1f\n",
              vno, inward_dist, outward_dist, current_sigma) ;
    }
              
    v->val2 = current_sigma ;
    /*
      search outwards and inwards and find the local gradient maximum
      at a location with a reasonable MR intensity value. This will
      be the location of the edge.
    */

    /* search in the normal direction to find the min value */
    max_mag_val = -10.0f ; mag = 5.0f ; max_mag = 0.0f ; min_val = 10000.0 ;
    min_val_dist = 0.0f ; local_max_found = 0 ;
    for (dist = inward_dist ; dist <= outward_dist ; dist += STEP_SIZE)
    {
#if 1
      x = v->x + v->nx*(dist-1) ; 
      y = v->y + v->ny*(dist-1) ; 
      z = v->z + v->nz*(dist-1) ;
      // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
      MRIsampleVolume(mri_brain, xw, yw, zw, &previous_val) ;
#else
      /* find max val within 1 mm in inwards direction */
      {
        float  d ;
        Real   tmp_val ;

        previous_val = 0 ;
        for (d = 0.25 ; d <= 1.5 ; d += 0.25)
        {
          x = v->x + v->nx*(d-1) ; 
          y = v->y + v->ny*(d-1) ; 
          z = v->z + v->nz*(d-1) ;
          // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
          MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
          MRIsampleVolume(mri_brain, xw, yw, zw, &tmp_val) ;
          if (tmp_val > previous_val)
            previous_val = tmp_val ;
        }
      }
#endif

      /* the previous point was inside the surface */
      if (previous_val < inside_hi && previous_val >= border_low)
      {
        /* see if we are at a local maximum in the gradient magnitude */
        x = v->x + v->nx*dist ; y = v->y + v->ny*dist ; z = v->z + v->nz*dist ;
        // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;

        x = v->x + v->nx*(dist+STEP_SIZE) ;
        y = v->y + v->ny*(dist+STEP_SIZE) ;
        z = v->z + v->nz*(dist+STEP_SIZE) ;
        // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny, nz,
                                       &next_mag, sigma);

        x = v->x + v->nx*(dist-STEP_SIZE) ;
        y = v->y + v->ny*(dist-STEP_SIZE) ;
        z = v->z + v->nz*(dist-STEP_SIZE) ;
        // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny, nz,
                                       &previous_mag, sigma);

        if (val < min_val)
        {
          min_val = val ;  /* used if no gradient max is found */
          min_val_dist = dist ;
        }

        /* if gradient is big and val is in right range */
        x = v->x + v->nx*dist;
        y = v->y + v->ny*dist;
        z = v->z + v->nz*dist;
        // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny, nz,&mag,
                                       sigma);
        if (which == GRAY_CSF)
        {
          /* 
             sample the next val we would process. If it is too low, then we 
             have definitely reached the border, and the current gradient 
             should be considered a local max.

             Don't want to do this for gray/white, as the gray/white gradient 
             often continues seemlessly into the gray/csf.
          */
          x = v->x + v->nx*(dist+STEP_SIZE) ;
          y = v->y + v->ny*(dist+STEP_SIZE) ;
          z = v->z + v->nz*(dist+STEP_SIZE) ;
          // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
          MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
          MRIsampleVolume(mri_brain, xw, yw, zw, &next_val) ;
          if (next_val < border_low)
            next_mag = 0 ;
        }

        if (vno == Gdiag_no)
          fprintf(fp, "%2.3f  %2.3f  %2.3f  %2.3f  %2.3f\n",
                 dist, val, mag, previous_mag, next_mag) ;

        /*
          if no local max has been found, or this one has a greater magnitude,
          and it is in the right intensity range....
          */
        if (
            /*            (!local_max_found || (fabs(mag) > max_mag)) && */
            (fabs(mag) > fabs(previous_mag)) &&
            (fabs(mag) > fabs(next_mag)) &&
            (val <= border_hi) && (val >= border_low)
            )
        {
          x = v->x + v->nx*(dist+1) ;
          y = v->y + v->ny*(dist+1) ;
          z = v->z + v->nz*(dist+1) ;
          // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
          MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
          MRIsampleVolume(mri_brain, xw, yw, zw, &next_val) ;
          /*
            if next val is in the right range, and the intensity at
            this local max is less than the one at the previous local
            max, assume it is the correct one.
          */
          if ((next_val >= outside_low) && 
              (next_val <= border_hi) &&
              (next_val <= outside_hi) &&
#if 0
              (!local_max_found || (val < max_mag_val)))
#else
              (!local_max_found || (max_mag < fabs(mag))))
#endif
          {          
            local_max_found = 1 ;
            max_mag_dist = dist ;
            max_mag = fabs(mag) ;
            max_mag_val = val ;
          }
        }
        else
        {
          /*
            if no local max found yet, just used largest gradient
            if the intensity is in the right range.
            */
          if ((local_max_found == 0) &&
              (fabs(mag) > max_mag) && 
              (val <= border_hi) && 
              (val >= border_low)
              )
          {
            x = v->x + v->nx*(dist+1) ;
            y = v->y + v->ny*(dist+1) ;
            z = v->z + v->nz*(dist+1) ;
            // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
            MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
            MRIsampleVolume(mri_brain, xw, yw, zw, &next_val) ;
            if (next_val >= outside_low && next_val <= border_hi &&
                next_val < outside_hi)
            {          
              max_mag_dist = dist ;
              max_mag = fabs(mag) ;
              max_mag_val = val ;
            }
          }
        }

      }
    }

    if (vno == Gdiag_no)
      fclose(fp) ;

    if (which == GRAY_CSF && local_max_found == 0 && max_mag_dist > 0)
    {
      float outlen ;
      int   allgray = 1 ;
      
      /* check to make sure it's not ringing near the gray white boundary,
	 by seeing if there is uniform stuff outside that could be gray matter.
      */
      for (outlen = max_mag_dist ; outlen < max_mag_dist+2 ; outlen += STEP_SIZE)
      {
	x = v->x + v->nx*outlen ;
	y = v->y + v->ny*outlen ;
	z = v->z + v->nz*outlen ;
	// MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
	MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
	MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;
	if ((val < outside_hi /*border_low*/) || (val > border_hi))
	{
	  allgray = 0 ;
	  break ;
	}
			}
      if (allgray)
      {
	if (Gdiag_no == vno)
	  printf("v %d: exterior gray matter detected, "
		 "ignoring large gradient at %2.3f (I=%2.1f)\n",
		 vno, max_mag_dist, max_mag_val) ;
	max_mag_val = -10 ;   /* don't worry about largest gradient */
				max_mag_dist = 0 ;
				num_changed++ ;
      }
    }

    if (max_mag_val > 0)   /* found the border value */
    {
      if (local_max_found)
        ngrad_max++ ;
      else
        ngrad++ ;
      if (max_mag_dist > 0)
      {
        nout++ ; nfound++ ;
        mean_out += max_mag_dist ;
      }
      else
      {
        nin++ ; nfound++ ;
        mean_in -= max_mag_dist ;
      }
      
      if (max_mag_val < border_low)
        max_mag_val = border_low ;
      mean_dist += max_mag_dist ;
      v->val = max_mag_val ;
      v->mean = max_mag ;
      mean_border += max_mag_val ; total_vertices++ ;
      v->d = max_mag_dist ;
      v->marked = 1 ;
    }
    else         /* couldn't find the border value */
    {
      if (min_val < 1000)
      {
        nmin++ ;
        v->d = min_val_dist ;
#if 0
        if (min_val > border_hi)  /* found a low value, but not low enough */
          min_val = border_hi ;
        else if (min_val < border_low)
          min_val = border_low ;
#else
        if (min_val < border_low)
          min_val = border_low ;
#endif
        v->val = min_val ;
        mean_border += min_val ; total_vertices++ ;
        v->marked = 1 ;
      }
      else
      {
        /* don't overwrite old target intensity if it was there */
        /*        v->val = -1.0f ;*/
        v->d = 0 ;
        if (v->val < 0)
        {
          nalways_missing++ ;
#if 0
          v->val = (border_low+border_hi)/2 ;
#endif
          v->marked = 0 ;
        }
        else
          v->marked = 1 ;
        nmissing++ ;
      }
    }
    if (vno == Gdiag_no)
      fprintf(stdout, 
              "v %d, target value = %2.1f, mag = %2.1f, dist = %2.2f, %s\n",
              Gdiag_no, v->val, v->mean, v->d,
              local_max_found ? "local max" : max_mag_val > 0 ? "grad":"min");
#if 0
    if (vno == 44289 || vno == 91080 || vno == 92286 || vno == 46922)
      fprintf(stdout, "v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
              Gdiag_no, v->val, v->mean, v->d) ;
#endif
  }
  mean_dist /= (float)(total_vertices-nmissing) ;
  mean_border /= (float)total_vertices ; 
#if 0
  mean_in /= (float)total_vertices ; 
  mean_out /= (float)total_vertices ; 
#else
  if (nin > 0)
    mean_in /= (float)nin ;
  if (nout > 0)
    mean_out /= (float)nout ;
#endif

#if 0
  MRISsoapBubbleVals(mris, 100) ; 
#endif

  /*  MRISaverageVals(mris, 3) ;*/
  fprintf(stdout, 
          "mean border=%2.1f, %d (%d) missing vertices, mean dist %2.1f "
          "[%2.1f (%%%2.1f)->%2.1f (%%%2.1f))]\n",
          mean_border, nmissing, nalways_missing, mean_dist, 
          mean_in, 100.0f*(float)nin/(float)nfound, 
          mean_out, 100.0f*(float)nout/(float)nfound) ;
  fprintf(stdout, "%%%2.0f local maxima, %%%2.0f large gradients "
          "and %%%2.0f min vals, %d gradients ignored\n",
          100.0f*(float)ngrad_max/(float)mris->nvertices,
          100.0f*(float)ngrad/(float)mris->nvertices,
          100.0f*(float)nmin/(float)mris->nvertices, num_changed) ;
  if (log_fp)
  {
    fprintf(log_fp, 
            "mean border=%2.1f, %d (%d) missing vertices, mean dist %2.1f "
            "[%2.1f (%%%2.1f)->%2.1f (%%%2.1f))]\n",
            mean_border, nmissing, nalways_missing, mean_dist, 
            mean_in, 100.0f*(float)nin/(float)nfound, 
            mean_out, 100.0f*(float)nout/(float)nfound) ;
    fprintf(log_fp, "%%%2.0f local maxima, %%%2.0f large gradients "
            "and %%%2.0f min vals, %d gradients ignored\n",
            100.0f*(float)ngrad_max/(float)mris->nvertices,
            100.0f*(float)ngrad/(float)mris->nvertices,
            100.0f*(float)nmin/(float)mris->nvertices, num_changed) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/

int
MRIScomputeGraySurfaceValues(MRI_SURFACE *mris,MRI *mri_brain,MRI *mri_smooth,
                             float gray_surface)
{
  Real    val, x, y, z, min_val, xw, yw, zw, dx, dy, dz, mag, max_mag ;
  int     total_vertices, vno, nmissing ;
  float   mean_gray, dist ;
  VERTEX  *v ;

  if (gray_surface <= 0.0f)
    gray_surface = MAX_CSF ;

  /* first compute intensity of local gray/white boundary */
  mean_gray = 0.0f ;

  MRISclearMarks(mris) ;   /* for use in soap bubble later */
  for (nmissing = total_vertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

  /*
    search outwards and inwards and find the local gradient maximum
    at a location with a reasonable MR intensity value. This will
    be the location of the edge.
    */

    /* search in the normal direction to find the min value */
    min_val = -10.0f ; mag = 5.0f ; max_mag = 0.0f ;
    for (dist = 0.0f ; dist < 6.0f ; dist += STEP_SIZE)
    {
      x = v->x + v->nx*dist ; y = v->y + v->ny*dist ; z = v->z + v->nz*dist ; 
      // MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
      MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;
      if (val < 70 && val > gray_surface)  /* in right range */
      {
        /* see if we are at a local maximum in the gradient magnitude */
        MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
        mag = sqrt(dx*dx + dy*dy + dz*dz) ;

        if (mag > max_mag)
        {
          max_mag = mag ;
          min_val = val ;
        }
      }
    }
    if (!finite(min_val) || !finite(max_mag) || !finite(mag))
      DiagBreak() ;
    if (min_val > 0)
    {
      v->marked = 1 ;
      v->val = min_val ;
      v->mean = max_mag ;
      mean_gray += min_val ; total_vertices++ ;
    }
    else
    {
      nmissing++ ;
      v->val = 0.0f ;
    }

    if (vno == Gdiag_no)
      fprintf(stdout, "v %d, target value = %2.1f, mag = %2.1f\n",
              Gdiag_no, v->val, v->mean) ;
  }
  mean_gray /= (float)total_vertices ; 
  MRISsoapBubbleVals(mris, 100) ; MRISclearMarks(mris) ;
  /*  MRISaverageVals(mris, 3) ;*/
  fprintf(stdout, "mean pial surface=%2.1f, %d missing\n", mean_gray,nmissing);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreverse(MRI_SURFACE *mris, int which)
{
  int    vno ;
  float  x, y, z ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    switch (which)
    {
    default:
    case REVERSE_X: x = -x ; break ;
    case REVERSE_Y: y = -y ; break ;
    case REVERSE_Z: z = -z ; break ;
    }
    v->x = x ;
    v->y = y ;
    v->z = z ;
  }
  if (which == REVERSE_X)   /* swap order of faces */
  {
    int  fno, vno0, vno1, vno2 ;
    FACE *f ;

    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      f = &mris->faces[fno] ;
      vno0 = f->v[0] ; vno1 = f->v[1] ; vno2 = f->v[2] ;
      f->v[0] = vno2 ; f->v[1] = vno1 ; f->v[2] = vno0 ;
      mrisSetVertexFaceIndex(mris, vno0, fno) ;
      mrisSetVertexFaceIndex(mris, vno1, fno) ;
      mrisSetVertexFaceIndex(mris, vno2, fno) ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
#define X_DOM   0
#define Y_DOM   1
#define Z_DOM   2


static int
mrisFindNormalDistanceLimits(MRI_SURFACE *mris, MRI *mri_filled, int vno, 
                             float max_dist,
                             float *pmax_outward_distance,
                             float *pmax_inward_distance)
{
  VERTEX   *v, *vn ;
  float    dx, dy, dz, len, max_outward, max_inward, dot, dist ;
  Real     nx, ny, nz, xw, yw, zw, x, y, z, x0, y0, z0,
           ax, ay, az ;
  int      n, xv, yv, zv, xv0, yv0, zv0, dom, xsign, ysign, zsign ;

  max_inward = max_outward = max_dist ;
  v = &mris->vertices[vno] ;
  
  /* compute normal vector in mri_filled coordinate system */
  MRIvoxelToWorld(mri_filled, 0, 0, 0, &xw, &yw, &zw) ;  /* origin */
  nx = xw + (Real)v->nx ; ny = yw + (Real)v->ny ; nz = zw + (Real)v->nz ;
  MRIworldToVoxel(mri_filled, nx, ny, nz, &nx, &ny, &nz) ;

  MRISvertexToVoxel(v, mri_filled, &x0, &y0, &z0) ;
  xv0 = nint(x0) ; yv0 = nint(y0) ; zv0 = nint(z0) ;
  xsign = nx > 0 ? 1 : -1 ; ysign = ny > 0 ? 1 : -1 ; zsign = nz > 0 ? 1 : -1 ;

  /* compute max outward movement without self-intersection */
  
  /* 
     first check to make sure no neihgbors occupy the same voxel in
     the filled volume and are in the normal direction. If they are,
     don't allow any movement in that direction.
  */
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    MRISvertexToVoxel(vn, mri_filled, &x, &y, &z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    dx = v->x - vn->x ; dy = v->y - vn->y ; dz = v->z - vn->z ;
    len = sqrt(SQR(dx*dx)+SQR(dy*dy)+SQR(dz*dz)) ;
    /*
      don't let vertices within the same filled voxel move towards each other
      and cross.
      */
    if ((len <= 2.5*MAX_MOVEMENT) && (xv == xv0 && yv == yv0 && zv == zv0))
    {
      dot = v->nx * dx + v->ny * dy + v->nz * dz ;
      if (dot < 0)   
        max_inward = 0 ;    /* don't allow inward movement */
      else
        max_outward = 0 ;   /* don't allow outward movement */
    }
  }

  /*
    now compute differentials, and use Bresenham type algorithm to
    move in 6-connected fashion in order to check for self-intersection.
  */
  ax = fabs(nx) ; ay = fabs(ny) ; az = fabs(nz) ;
  if (ax > ay && ax > az)    /* nx biggest - set it to unit movement */
  {
    len = nx ; nx = (Real)xsign ; ny /= len ; nz /= len ; dom = X_DOM ;
  }
  else if (ay > az)          /* ny biggest - set it to unit movement */
  {
    len = ny ; ny = (Real)ysign ; nx /= len ; nz /= len ; dom = Y_DOM ;
  }
  else                       /* nz biggest - set it to unit movement */
  {
    len = nz ; nz = (Real)zsign ; ny /= len ; nx /= len ; dom = Z_DOM ;
  }
  ax = fabs(nx) ; ay = fabs(ny) ; az = fabs(nz) ;

  dx = ax ; dy = ay ; dz = az ;
  xv = xv0 ; yv = yv0 ; zv = zv0 ; x = x0 ; y = y0 ; z = z0 ;
  dist = 0.0f ;
  do
  {
    if (dom == X_DOM)   /* take unit step in x direction, and possibly other */
    {
      if ((dy > dz) && (dy >= 1.0f))   /* take step in y direction */
      {
        yv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y += ny ;
        dy -= 1.0 ;
      }
      if (dz >= 1.0f)  /* take unit step in z direction */
      {
        zv += zsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        z += nz ;
        dz -= 1.0 ;
      }
      if (dy >= 1.0f)  /* take unit step in z direction */
      {
        yv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y += ny ;
        dy -= 1.0 ;
      }
      xv += xsign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dx -= 1.0 ;
      x += nx ;
    }
    else if (dom == Y_DOM)
    {
      if ((dx > dz) && (dx >= 1.0f))   /* take step in x direction */
      {
        xv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x += nx ;
        dx -= 1.0 ;
      }
      if (dz >= 1.0f)  /* take unit step in z direction */
      {
        zv += zsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        z += nz ;
        dz -= 1.0 ;
      }
      if (dx >= 1.0f)  /* take unit step in x direction */
      {
        xv += xsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x += ny ;
        dx -= 1.0 ;
      }
      yv += ysign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dy -= 1.0 ;
      y += ny ;
    }
    else    /* z dominant - take unit step in z direction (at least) */
    {
      if ((dy > dx) && (dy >= 1.0f))   /* take step in y direction */
      {
        yv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y += ny ;
        dy -= 1.0 ;
      }
      if (dx >= 1.0f)  /* take unit step in x direction */
      {
        xv += xsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x += nx ;
        dx -= 1.0 ;
      }
      if (dy >= 1.0f)  /* take unit step in y direction */
      {
        yv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y += ny ;
        dy -= 1.0 ;
      }
      zv += zsign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dz -= 1.0 ;
      z += nz ;
    }
    dist = sqrt(SQR(mri_filled->xsize*(x-x0))+
                SQR(mri_filled->ysize*(y-y0))+
                SQR(mri_filled->zsize*(z-z0))) ;
    dx += ax ; dy += ay ; dz += az ;
  } while (dist < max_outward) ;

  if (dist > max_outward)
    dist = max_outward ;
  *pmax_outward_distance = max_outward = dist ;

  /*
    now do the same thing in the inward direction
  */
  dx = ax ; dy = ay ; dz = az ;
  xv = xv0 ; yv = yv0 ; zv = zv0 ; x = x0 ; y = y0 ; z = z0 ;
  dist = 0.0f ;
  do
  {
    if (dom == X_DOM)   /* take unit step in x direction, and possibly other */
    {
      if ((dy > dz) && (dy >= 1.0f))   /* take step in y direction */
      {
        yv -= ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y -= ny ;
        dy -= 1.0 ;
      }
      if (dz >= 1.0f)  /* take unit step in z direction */
      {
        zv -= zsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        z -= nz ;
        dz -= 1.0 ;
      }
      if (dy >= 1.0f)  /* take unit step in z direction */
      {
        yv -= ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y -= ny ;
        dy -= 1.0 ;
      }
      xv -= xsign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dx -= 1.0 ;
      x -= nx ;
    }
    else if (dom == Y_DOM)
    {
      if ((dx > dz) && (dx >= 1.0f))   /* take step in x direction */
      {
        xv -= ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x -= nx ;
        dx -= 1.0 ;
      }
      if (dz >= 1.0f)  /* take unit step in z direction */
      {
        zv -= zsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        z -= nz ;
        dz -= 1.0 ;
      }
      if (dx >= 1.0f)  /* take unit step in x direction */
      {
        xv -= xsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x -= ny ;
        dx -= 1.0 ;
      }
      yv -= ysign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dy -= 1.0 ;
      y -= ny ;
    }
    else    /* z dominant - take unit step in z direction (at least) */
    {
      if ((dy > dx) && (dy >= 1.0f))   /* take step in y direction */
      {
        yv -= ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y -= ny ;
        dy -= 1.0 ;
      }
      if (dx >= 1.0f)  /* take unit step in x direction */
      {
        xv -= xsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x -= nx ;
        dx -= 1.0 ;
      }
      if (dy >= 1.0f)  /* take unit step in y direction */
      {
        yv -= ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y -= ny ;
        dy -= 1.0 ;
      }
      zv -= zsign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dz -= 1.0 ;
      z -= nz ;
    }
    dist = sqrt(SQR(mri_filled->xsize*(x-x0))+
                SQR(mri_filled->ysize*(y-y0))+
                SQR(mri_filled->zsize*(z-z0))) ;
    dx += ax ; dy += ay ; dz += az ;
  } while (dist < max_inward) ;

  if (dist > max_inward)
    dist = max_inward ;
  if (!FZERO(dist))
    DiagBreak() ;
  *pmax_inward_distance = max_inward = dist ;


  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisLimitGradientDistance(MRI_SURFACE *mris, MHT *mht, int vno)
{
  VERTEX   *v ;

  v = &mris->vertices[vno] ;

  mrisRemoveNeighborGradientComponent(mris, vno) ;
  if (MHTisVectorFilled(mht, mris, vno, v->odx, v->ody, v->odz))
  {
    mrisRemoveNormalGradientComponent(mris, vno) ;
    if (MHTisVectorFilled(mht, mris, vno, v->odx, v->ody, v->odz))
    {
      v->odx = v->ody = v->odz = 0.0 ;
      return(NO_ERROR) ;
    }
  }

  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisLimitGradientDistance(MRI_SURFACE *mris, MRI *mri_filled, int vno)
{
  VERTEX   *v, *vn ;
  float    dx, dy, dz, len, max_outward, dot, dist ;
  Real     nx, ny, nz, xw, yw, zw, x, y, z, x0, y0, z0,
           ax, ay, az ;
  int      n, xv, yv, zv, xv0, yv0, zv0, dom, xsign, ysign, zsign ;

  v = &mris->vertices[vno] ;
  max_outward = sqrt(SQR(v->odx)+SQR(v->ody)+SQR(v->odz)) ;

  /* compute gradient in mri_filled coordinate system */
  MRIvoxelToWorld(mri_filled, 0, 0, 0, &xw, &yw, &zw) ;  /* origin */
  dx = v->odx ; dy = v->ody ; dz = v->odz ;
  len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
  if (FZERO(len))
    return(NO_ERROR) ;
  dx /= len ; dy /= len ; dz /= len ;
  nx = xw + (Real)dx ; ny = yw + (Real)dy ; nz = zw + (Real)dz ;
  MRIworldToVoxel(mri_filled, nx, ny, nz, &nx, &ny, &nz) ;

  MRISvertexToVoxel(v, mri_filled, &x0, &y0, &z0) ;
  xv0 = nint(x0) ; yv0 = nint(y0) ; zv0 = nint(z0) ;
  xsign = nx > 0 ? 1 : -1 ; ysign = ny > 0 ? 1 : -1 ; zsign = nz > 0 ? 1 : -1 ;

  /* compute max outward movement without self-intersection */
  
  /* 
     first check to make sure no neighbors occupy the same voxel in
     the filled volume and are in the gradient direction. If they are,
     don't allow any movement in that direction.
  */
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    MRISvertexToVoxel(vn, mri_filled, &x, &y, &z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    dx = v->x - vn->x ; dy = v->y - vn->y ; dz = v->z - vn->z ;
    len = sqrt(SQR(dx*dx)+SQR(dy*dy)+SQR(dz*dz)) ;
    /*
      don't let vertices within the same filled voxel move towards each other
      and cross.
      */
    if ((len <= 2.5*MAX_MOVEMENT) && (xv == xv0 && yv == yv0 && zv == zv0))
    {
      dot = v->odx * dx + v->ody * dy + v->odz * dz ;
      if (dot >= 0)   
      {
        /* take out component toward neighbor */
#if 1
        v->odx -= dot * dx ; v->ody -= dot * dy ; v->odz -= dot * dz ;
#else
        max_outward = 0 ;   /* don't allow outward movement */
#endif
      }
    }
  }

  /*
    now compute differentials, and use Bresenham type algorithm to
    move in 6-connected fashion in order to check for self-intersection.
  */
  ax = fabs(nx) ; ay = fabs(ny) ; az = fabs(nz) ;
  if (ax > ay && ax > az)    /* nx biggest - set it to unit movement */
  {
    len = nx ; nx = (Real)xsign ; ny /= len ; nz /= len ; dom = X_DOM ;
  }
  else if (ay > az)          /* ny biggest - set it to unit movement */
  {
    len = ny ; ny = (Real)ysign ; nx /= len ; nz /= len ; dom = Y_DOM ;
  }
  else                       /* nz biggest - set it to unit movement */
  {
    len = nz ; nz = (Real)zsign ; ny /= len ; nx /= len ; dom = Z_DOM ;
  }
  ax = fabs(nx) ; ay = fabs(ny) ; az = fabs(nz) ;

  dx = ax ; dy = ay ; dz = az ;
  xv = xv0 ; yv = yv0 ; zv = zv0 ; x = x0 ; y = y0 ; z = z0 ;
  dist = 0.0f ;
  do
  {
    if (dom == X_DOM)   /* take unit step in x direction, and possibly other */
    {
      if ((dy > dz) && (dy >= 1.0f))   /* take step in y direction */
      {
        yv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y += ny ;
        dy -= 1.0 ;
      }
      if (dz >= 1.0f)  /* take unit step in z direction */
      {
        zv += zsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        z += nz ;
        dz -= 1.0 ;
      }
      if (dy >= 1.0f)  /* take unit step in z direction */
      {
        yv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y += ny ;
        dy -= 1.0 ;
      }
      xv += xsign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dx -= 1.0 ;
      x += nx ;
    }
    else if (dom == Y_DOM)
    {
      if ((dx > dz) && (dx >= 1.0f))   /* take step in x direction */
      {
        xv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x += nx ;
        dx -= 1.0 ;
      }
      if (dz >= 1.0f)  /* take unit step in z direction */
      {
        zv += zsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        z += nz ;
        dz -= 1.0 ;
      }
      if (dx >= 1.0f)  /* take unit step in x direction */
      {
        xv += xsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x += ny ;
        dx -= 1.0 ;
      }
      yv += ysign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dy -= 1.0 ;
      y += ny ;
    }
    else    /* z dominant - take unit step in z direction (at least) */
    {
      if ((dy > dx) && (dy >= 1.0f))   /* take step in y direction */
      {
        yv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y += ny ;
        dy -= 1.0 ;
      }
      if (dx >= 1.0f)  /* take unit step in x direction */
      {
        xv += xsign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        x += nx ;
        dx -= 1.0 ;
      }
      if (dy >= 1.0f)  /* take unit step in y direction */
      {
        yv += ysign ;
        if (MRItest_bit(mri_filled, xv, yv, zv))
        {
          if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
            break ;
        }
        y += ny ;
        dy -= 1.0 ;
      }
      zv += zsign ;
      if (MRItest_bit(mri_filled, xv, yv, zv))
      {
        if (!mrisNeighborAtVoxel(mris, mri_filled, vno, xv, yv, zv))
          break ;
      }
      dz -= 1.0 ;
      z += nz ;
    }
    dist = sqrt(SQR(mri_filled->xsize*(x-x0))+
                SQR(mri_filled->ysize*(y-y0))+
                SQR(mri_filled->zsize*(z-z0))) ;
    dx += ax ; dy += ay ; dz += az ;
  } while (dist < max_outward) ;

  if (dist > max_outward)
    dist = max_outward ;

  if (FZERO(dist))
    v->odx = v->ody = v->odz = 0.0f ;
  else   /* limit gradient to max_outward length */
  { 
    dx = v->odx ; dy = v->ody ; dz = v->odz ;
    len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
    v->odx = max_outward * v->odx / len ;
    v->ody = max_outward * v->ody / len ;
    v->odz = max_outward * v->odz / len ;
  }
  return(NO_ERROR) ;
}
#endif
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisDebugVertex(MRI_SURFACE *mris, int vno)
{
  int     n ;
  VERTEX  *v, *vn ;
  float   d, dx, dy, dz ;

  v = &mris->vertices[vno] ;
  fprintf(stdout, 
          "vertex #%d @ (%2.2f, %2.2f, %2.2f), n = (%2.2f, %2.2f, %2.2f) "
          "(%2.2f, %2.2f, %2.2f), val=%2.1f\n",
          vno, v->x, v->y, v->z, v->nx, v->ny, v->nz, v->dx, v->dy, v->dz,
          v->val) ;

  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    dx = vn->x - v->x ; dy = vn->y - v->y ; dz = vn->z - v->z ;
    d = sqrt(dx*dx + dy*dy + dz*dz) ;
    fprintf(stdout, 
            "\tn %d: %6.6d, delta = (%2.3f, %2.3f, %2.3f), dist = %2.3f "
            "(%2.2f, %2.2f, %2.2f), val=%2.1f\n",
            n, v->v[n], dx, dy, dz, d, vn->dx, vn->dy, vn->dz, vn->val) ;
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisRmsValError(MRI_SURFACE *mris, MRI *mri)
{
  int     vno, n, xv, yv, zv ;
  Real    val, total, delta, x, y, z ;
  VERTEX  *v ;

  for (total = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    n++ ;
    MRISvertexToVoxel(v, mri, &x, &y, &z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    MRIsampleVolume(mri, x, y, z, &val) ;
    delta = (val - v->val) ;
    total += delta*delta ;
  }
  return(sqrt(total / (double)n)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisNeighborAtVoxel(MRI_SURFACE *mris, MRI *mri, int vno, int xv,int yv,int zv)
{
  int      n, xnv, ynv, znv ;
  VERTEX   *v, *vn ;
  Real     xn, yn, zn ;

  v = &mris->vertices[vno] ;
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    MRISvertexToVoxel(v, mri, &xn, &yn, &zn) ;
    xnv = nint(xn) ; ynv = nint(yn) ; znv = nint(zn) ;
    if (xnv == xv && ynv == yv && znv == zv)
      return(1) ;
  }
  return(0) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRIScomputeAnalyticDistanceError(MRI_SURFACE *mris, int which, FILE *fp)
{
  int     vno, n, vtotal, *pv, ndists ;
  VERTEX  *v, *vn ;
  float   d, xd, yd, zd, circumference = 0.0f, angle, odist ;
  double  pct_orig, pct, mean, mean_orig_error, mean_error,
          smean_error, smean_orig_error ;
  VECTOR  *v1, *v2 ;

  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;

  mean_orig_error = mean_error = pct_orig = pct=  mean = 0.0 ;
  smean_orig_error = smean_error = 0.0 ;
  ndists = 0 ;
  for (vno=0;vno<mris->nvertices;vno++)
  {
    v = &mris->vertices[vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vtotal = v->vtotal ;
    switch (which)
    {
    default:   /* don't really know what to do in other cases */
    case MRIS_PLANE:
      for (pv = v->v, n = 0 ; n < vtotal ; n++)
      {
        vn = &mris->vertices[*pv++] ;
        if (vn->ripflag)
          continue ;
        xd = v->origx - vn->origx ; yd = v->origy - vn->origy ; 
        zd = v->origz - vn->origz ;
        d = xd*xd + yd*yd + zd*zd ;
        odist = sqrt(d) ;
        mean_orig_error += fabs(v->dist_orig[n] - odist) ;
        mean_error += fabs(v->dist[n] - odist) ;
        smean_orig_error += (v->dist_orig[n] - odist) ;
        smean_error += (v->dist[n] - odist) ;
        mean += odist ;
        ndists++ ;
      }
      break ;
    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE:
      VECTOR_LOAD(v1, v->origx, v->origy, v->origz) ;  /* radius vector */
      if (FZERO(circumference))   /* only calculate once */
        circumference = M_PI * 2.0 * V3_LEN(v1) ;
      for (pv = v->v, n = 0 ; n < vtotal ; n++)
      {
        vn = &mris->vertices[*pv++] ;
        if (vn->ripflag)
          continue ;
        VECTOR_LOAD(v2, vn->origx, vn->origy, vn->origz) ;  /* radius vector */
        angle = fabs(Vector3Angle(v1, v2)) ;
        d = circumference * angle / (2.0 * M_PI) ;
        odist = d ;
        mean_orig_error += fabs(v->dist_orig[n] - odist) ;
        mean_error += fabs(v->dist[n] - odist) ;
        smean_orig_error += (v->dist_orig[n] - odist) ;
        smean_error += (v->dist[n] - odist) ;
        mean += fabs(odist) ;
        ndists++ ;
      }
      break ;
    }
  }

  mean /= (double)ndists ; mean_error /= (double)ndists ;
  mean_orig_error /= (double)ndists ; 
  smean_orig_error /= (double)ndists ; smean_error /= (double)ndists ; 
  pct = mean_error / mean ; pct_orig = mean_orig_error / mean ;
  fprintf(stdout, 
      "mean orig = %2.3f mm (%%%2.2f), final = %2.3f mm (%%%2.2f)\n",
          mean_orig_error, 100.0*pct_orig, mean_error, 100.0*pct) ;
  fprintf(stdout, "signed mean orig error = %2.3f, final mean error = %2.3f\n",
          smean_orig_error, smean_error) ;
  if (fp)
  {
    char  *cp ;
    float measured_error, disturb_pct ;

    cp = getenv("DISTURB_DISTANCES") ;
    if (cp)
      disturb_pct = atof(cp) ;
    else
      disturb_pct = 0.0 ;
    measured_error = MRISpercentDistanceError(mris) ;
    fprintf(fp, "%2.3f  %2.3f  %2.3f  %2.3f  %2.3f  %2.3f  %2.3f\n",
            100.0f*(float)mrisValidVertices(mris) / (float)mris->nvertices,
            disturb_pct,
            100.0*pct_orig, mean_orig_error, 100.0*pct, mean_error, 
            measured_error) ;
#if 0
    fprintf(fp, 
            "mean orig = %2.3f mm (%%%2.2f), final = %2.3f mm (%%%2.2f)\n",
            mean_orig_error, 100.0*pct_orig, mean_error, 100.0*pct) ;
    fprintf(fp, "signed mean orig error = %2.3f, final mean error = %2.3f\n",
            smean_orig_error, smean_error) ;
#endif
  }
  VectorFree(&v1) ; VectorFree(&v2) ;
  return(pct) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRISstoreAnalyticDistances(MRI_SURFACE *mris, int which)
{
  int     vno, n, vtotal, *pv ;
  VERTEX  *v, *vn ;
  float   d, xd, yd, zd, circumference = 0.0f, angle, odist ;
  double  pct_orig, pct, mean, mean_orig_error, mean_error,
          smean_error, smean_orig_error ;
  VECTOR  *v1, *v2 ;

  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;

  mean_orig_error = mean_error = pct_orig = pct=  mean = 0.0 ;
  smean_orig_error = smean_error = 0.0 ;
  for (vno=0;vno<mris->nvertices;vno++)
  {
    v = &mris->vertices[vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vtotal = v->vtotal ;
    switch (which)
    {
    default:   /* don't really know what to do in other cases */
    case MRIS_PLANE:
      for (pv = v->v, n = 0 ; n < vtotal ; n++)
      {
        vn = &mris->vertices[*pv++] ;
        if (vn->ripflag)
          continue ;
        xd = v->origx - vn->origx ; yd = v->origy - vn->origy ; 
        zd = v->origz - vn->origz ;
        d = xd*xd + yd*yd + zd*zd ;
        odist = sqrt(d) ;
        v->dist_orig[n] = odist ;
      }
      break ;
    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE:
      VECTOR_LOAD(v1, v->origx, v->origy, v->origz) ;  /* radius vector */
      if (FZERO(circumference))   /* only calculate once */
        circumference = M_PI * 2.0 * V3_LEN(v1) ;
      for (pv = v->v, n = 0 ; n < vtotal ; n++)
      {
        vn = &mris->vertices[*pv++] ;
        if (vn->ripflag)
          continue ;
        VECTOR_LOAD(v2, vn->origx, vn->origy, vn->origz) ;  /* radius vector */
        angle = fabs(Vector3Angle(v1, v2)) ;
        d = circumference * angle / (2.0 * M_PI) ;
        v->dist_orig[n] = d ;
      }
      break ;
    }
  }

  VectorFree(&v1) ; VectorFree(&v2) ;
  return(pct) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISdisturbOriginalDistances(MRI_SURFACE *mris, double max_pct)
{
  int    vno, n ;
  VERTEX *v ;
  

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    for (n = 0 ; n < v->vtotal ; n++)
    {
      v->dist_orig[n] *= (1.0 + randomNumber(-max_pct/100.0f, max_pct/100.0f));
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISnegateValues(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val *= -1.0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyMeansToValues(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = v->mean ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyImaginaryMeansToValues(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = v->mean_imag ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyStandardErrorsToValues(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = v->std_error ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaccumulateMeansInVolume(MRI_SURFACE *mris, MRI *mri, int mris_dof, 
                            int mri_dof, int coordinate_system, int sno)
{
  VERTEX    *vertex ;
  Real      ndof, x, y, z, mean ;
  int       vno, xv, yv, zv ;
  
  ndof = (Real)(mris_dof + mri_dof) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->val != 0.0f)
    {
      x = vertex->x ; y = vertex->y ; z = vertex->z ;
      switch (coordinate_system)
      {
      case TALAIRACH_COORDS:
        MRISworldToTalairachVoxel(mris, mri, x, y, z, &x, &y, &z) ;
        break ;
      default:  /* surface-based */
        // MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;
        MRIsurfaceRASToVoxel(mri, x, y, z, &x, &y, &z) ;
        break ;
      }
      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
      if ((xv >= 0 && xv < mri->width) &&
          (yv >= 0 && yv < mri->height) &&
          (zv >= 0 && zv < mri->depth))
      {
        mean = MRIFseq_vox(mri, xv, yv, zv, sno) ;
        mean = (mris_dof * vertex->val + mri_dof * mean) / ndof ;
        MRIFseq_vox(mri, xv, yv, zv, sno) = mean ;
      }
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        actually these are squared standard errors 
------------------------------------------------------*/
int
MRISaccumulateStandardErrorsInVolume(MRI_SURFACE *mris, MRI *mri, 
                                         int mris_dof, int mri_dof, 
                                         int coordinate_system, int sno)
{
  VERTEX    *vertex ;
  Real      ndof, x, y, z, mris_sigma, mri_sigma ;
  int       vno, xv, yv, zv ;
  
  ndof = (Real)(mris_dof + mri_dof) ;

/* 
   now that we have the values read in, go through the surface, and
   map each vertex into the structural volume via its ellipsoidal
   coordinate.
   */

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->val != 0.0f)
    {
      x = vertex->x ; y = vertex->y ; z = vertex->z ;
      switch (coordinate_system)
      {
      case TALAIRACH_COORDS:
        MRISworldToTalairachVoxel(mris, mri, x, y, z, &x, &y, &z) ;
        break ;
      default:  /* surface-based */
        // MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;
        MRIsurfaceRASToVoxel(mri, x, y, z, &x, &y, &z) ;
        break ;
      }
      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
      
      if ((xv < 0 || xv >= mri->width) ||
          ((yv < 0) || (yv >= mri->height)) ||
          ((zv < 0) || (zv >= mri->depth)))
        continue ;
      
      mri_sigma = MRIFseq_vox(mri, xv, yv, zv, sno) ; /* variance */
      mris_sigma = vertex->val  ;
      mri_sigma = mris_sigma*SQR(mris_dof) + mri_sigma*SQR(mri_dof) ;
      mri_sigma /= SQR(ndof) ;
      if (!finite(mri_sigma))
      {
        fprintf(stderr, "variance not finite at vno %d!\n", vno) ;
        DiagBreak() ;
      }
      MRIFseq_vox(mri, xv, yv, zv, sno) = mri_sigma ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaccumulateMeansOnSurface(MRI_SURFACE *mris, int total_dof,int new_dof)
{
  int    vno, ndof ;
  VERTEX *v ;

  ndof = total_dof + new_dof ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    v->mean = (v->mean*total_dof + v->val*new_dof) / ndof ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISaccumulateImaginaryMeansOnSurface(MRI_SURFACE *mris, int total_dof,
                                      int new_dof)
{
  int    vno, ndof ;
  VERTEX *v ;

  ndof = total_dof + new_dof ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    v->mean_imag = (v->mean_imag*total_dof + v->val*new_dof) / ndof ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        actually these are squared standard errors 
------------------------------------------------------*/
int
MRISaccumulateStandardErrorsOnSurface(MRI_SURFACE *mris,
                                      int total_dof,int new_dof)
{
  int    vno, ndof ;
  VERTEX *v ;
  double var ;

  ndof = total_dof + new_dof ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    var = (v->std_error * SQR(total_dof) + v->val * SQR(new_dof)) / SQR(ndof);
    v->std_error = var ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        actually these are squared standard errors 
------------------------------------------------------*/
int
MRIScomputeAverageCircularPhaseGradient(MRI_SURFACE *mris, LABEL *area,
                                        float *pdx, float *pdy, float *pdz)
{
  int     N, vno, n, i ;
  VERTEX  *v, *vn ;
  VECTOR  *vdf, *vfz ;
  MATRIX  *mz, *mzt, *mztz, *mztz_inv, *mztz_inv_zt ;
  double  x0, y0, z0, f0, x1, y1, z1, f1, dx, dy, dz, df ;

  dx = dy = dz = 0.0 ;
  for (i = 0 ; i < area->n_points ; i++)
  {
    vno = area->lv[i].vno ;
    v = &mris->vertices[vno] ;
    x0 = v->x ; y0 = v->y ; z0 = v->z ; f0 = atan2(v->imag_val, v->val) ;

    /* first count # of valid neighbors */
    for (N = n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag)
        continue ;
      N++ ;
    }
    /* now allocate vectors and matrix */
    vdf = VectorAlloc(N, MATRIX_REAL) ;    /* function deltas */
    if (mris->patch)
      mz = MatrixAlloc(N, 2, MATRIX_REAL) ;  /* vertex spacing deltas */
    else
      mz = MatrixAlloc(N, 3, MATRIX_REAL) ;  /* vertex spacing deltas */

    /* now fill in matrix and vector entries */
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag)
        continue ;
      x1 = vn->x ; y1 = vn->y ; z1 = vn->z ; f1 = atan2(vn->imag_val,vn->val);
      df = f1 - f0 ;
      while (df > M_PI)
        df -= 2*M_PI ;
      while (df < -M_PI)
        df += 2*M_PI ;
      VECTOR_ELT(vdf, n+1) = df ;
      *MATRIX_RELT(mz, n+1, 1) = x1 - x0 ;
      *MATRIX_RELT(mz, n+1, 2) = y1 - y0 ;
      if (!mris->patch)
        *MATRIX_RELT(mz, n+1, 3) = z1 - z0 ;
    }
    mzt = MatrixTranspose(mz, NULL) ;
    mztz = MatrixMultiply(mzt, mz, NULL) ;
    mztz_inv = MatrixSVDInverse(mztz, NULL) ;
    if (mztz_inv)
    {
      mztz_inv_zt = MatrixMultiply(mztz_inv, mzt, NULL) ;
      vfz = MatrixMultiply(mztz_inv_zt, vdf, NULL) ;
      v->dx = VECTOR_ELT(vfz, 1) ;
      v->dy = VECTOR_ELT(vfz, 2) ;
      if (!mris->patch)
        v->dz = VECTOR_ELT(vfz, 3) ;
      else
        v->dz = 0.0f ;
      dx += v->dx ; dy += v->dy ; dz += v->dz ;
    }

    VectorFree(&vdf) ; 
    MatrixFree(&mz) ;
    MatrixFree(&mzt) ;
    MatrixFree(&mztz) ;
    if (mztz_inv)
    {
      VectorFree(&vfz) ;
      MatrixFree(&mztz_inv) ;
      MatrixFree(&mztz_inv_zt) ;
    }
  }
  
  *pdx = dx /= (float)area->n_points ; *pdy = dy /= (float)area->n_points ;
  *pdz = dz /= (float)area->n_points ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeIntensityError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z ;
  Real    val0, xw,yw,zw ;
  double  sse, del0 ;

  if (FZERO(parms->l_intensity))
    return(0.0f) ;

  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ; y = v->y ; z = v->z ;
    
    // MRIworldToVoxel(parms->mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(parms->mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolume(parms->mri_brain, xw, yw, zw, &val0) ;

    del0 = v->val - val0 ;
    sse += (del0 * del0) ;
  }
  
  return(parms->l_intensity * sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeIntensityGradientError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z ;
  Real    mag0, xw, yw, zw, dx, dy, dz ;
  double  sse, del0 ;

  if (FZERO(parms->l_grad))
    return(0.0f) ;

  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ; y = v->y ; z = v->z ;
    
    // MRIworldToVoxel(parms->mri_smooth, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(parms->mri_smooth, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolumeGradient(parms->mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    mag0 = sqrt(dx*dx + dy*dy + dz*dz) ;

    del0 = v->mean - mag0 ;
    sse += (del0 * del0) ;
  }
  
  return(parms->l_grad * sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Compute the ratio of the two principal curvatures and
           store it in the vertex->curv variable.
------------------------------------------------------*/
int
MRISuseCurvatureRatio(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;
  float  min_curv, max_curv, k1, k2, curv ;

  /*  MRIScomputeSecondFundamentalForm(mris) ;*/

  min_curv = 10000.0f ; max_curv = -min_curv ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (fabs(v->k1) > fabs(v->k2))
    { k1 = v->k1 ; k2 = v->k2 ; }
    else
    { k1 = v->k2 ; k2 = v->k1 ; }

    if (!FZERO(k2))
    {
      curv = fabs(k1 / k2) ;
      if (curv < min_curv)
        min_curv = curv ;
      if (curv > max_curv)
        max_curv = curv ;

      v->curv = curv ;
    }
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (fabs(v->k1) > fabs(v->k2))
    { k1 = v->k1 ; k2 = v->k2 ; }
    else
    { k1 = v->k2 ; k2 = v->k1 ; }

    if (FZERO(k2))
    {
      if (FZERO(k1))
        curv = 0.0 ;
      else
        curv = k1 < 0 ? min_curv : max_curv ;
      v->curv = curv ;
    }
  }

  mris->min_curv = min_curv ; mris->max_curv = max_curv ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Compute the contrast of the two principal curvatures and
           store it in the vertex->curv variable.
------------------------------------------------------*/
int
MRISuseCurvatureContrast(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;
  float  min_curv, max_curv, k1, k2, curv, min_k ;

  /*  MRIScomputeSecondFundamentalForm(mris) ;*/

  min_k = 10000.0f ; 
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->k1 < min_k)
      min_k = v->k1 ;
    if (v->k2 < min_k)
      min_k = v->k2 ;
  }

  min_curv = 10000.0f ; max_curv = -min_curv ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->k1 > v->k2)
    { k1 = v->k1 ; k2 = v->k2 ; }
    else
    { k1 = v->k2 ; k2 = v->k1 ; }
    k1 -= min_k ; k2 -= min_k ;
    curv = (k1 - k2) / (k1+k2) ;

    if (curv < min_curv)
      min_curv = curv ;
    if (curv > max_curv)
      max_curv = curv ;

    v->curv = curv ;
  }

  mris->min_curv = min_curv ; mris->max_curv = max_curv ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeSphereError(MRI_SURFACE *mris, double l_sphere, double r0)
{
  int     vno ;
  VERTEX  *v ;
  double  del, sse, x, y, z, x0, y0, z0, r ;

  if (FZERO(l_sphere))
    return(0.0f) ;

  x0 = (mris->xlo+mris->xhi)/2.0f ; 
  y0 = (mris->ylo+mris->yhi)/2.0f ; 
  z0 = (mris->zlo+mris->zhi)/2.0f ;
  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = (double)v->x-x0 ; 
    y = (double)v->y-y0 ; 
    z = (double)v->z-z0 ;
    r = sqrt(x*x + y*y + z*z) ;
    
    del = r0 - r ;
    sse += (del * del) ;
#if 0
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d sphere term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
#endif
  }
  
  return(l_sphere * sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScanonicalToWorld(MRI_SURFACE *mris, Real phi, Real theta,
                                  Real *pxw, Real *pyw, Real *pzw)
{
  Real x, y, z, radius ;
  
  radius = mris->radius ;
  *pxw = x = radius * sin(phi) * cos(theta) ;
  *pyw = y = radius * sin(phi) * sin(theta) ;
  *pzw = z = radius * cos(phi) ;
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisClipGradient(MRI_SURFACE *mris, float max_len)
{
  int     vno ;
  VERTEX  *v ;
  float   dx, dy, dz, len ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dx = v->dx ; dy = v->dy ; dz = v->dz ; 
    len = sqrt(dx*dx+dy*dy+dz*dz) ;
    if (len > max_len)
    {
      len = max_len / len ;
      v->dx *= len ; v->dy *= len ; v->dz *= len ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisClipMomentumGradient(MRI_SURFACE *mris, float max_len)
{
  int     vno ;
  VERTEX  *v ;
  float   dx, dy, dz, len ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dx = v->odx ; dy = v->ody ; dz = v->odz ; 
    len = sqrt(dx*dx+dy*dy+dz*dz) ;
    if (len > max_len)
    {
      len = max_len / len ;
      v->odx *= len ; v->ody *= len ; v->odz *= len ;
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyCurvatureToValues(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = v->curv ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyCurvatureToImagValues(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->imag_val = v->curv ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyValuesToCurvature(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = v->val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyImagValuesToCurvature(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = v->imag_val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyCurvatureFromValues(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = v->val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyCurvatureFromImagValues(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = v->imag_val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScomputeDistanceErrors(MRI_SURFACE *mris, int nbhd_size, int max_nbrs)
{
  VERTEX  *v ;
  int     vno, n, nvertices ;
  double  dist_scale, pct, dist, odist, mean, mean_error, smean,
          total_mean_error, total_mean ;

  MRIScomputeMetricProperties(mris) ;
  if (mris->patch)
    dist_scale = 1.0 ;
  else
    dist_scale = sqrt(mris->orig_area / mris->total_area) ;

  total_mean = total_mean_error = mean = 0.0 ;
  for (pct = 0.0, nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      v->val = 0.0f ;
      continue ;
    }
    for (smean = mean_error = mean = 0.0, n = 0 ; n < v->vtotal ; n++)
    {
      nvertices++ ;
      dist = dist_scale*v->dist[n] ;
      odist = v->dist_orig[n] ;
#if 0
      mean += (dist - odist) * (dist - odist) ;
#else
      mean += odist ;
#endif
      total_mean += odist ;
      smean += dist - odist ;
#define USE_FABS  1
#if USE_FABS
      mean_error += fabs(dist-odist) ;
      total_mean_error += fabs(dist-odist) ;
#else
      mean_error += (dist-odist)*(dist-odist) ;
      total_mean_error += (dist-odist)*(dist-odist) ;
#endif
      if (!FZERO(odist))
        pct += fabs(dist - odist) / odist ;
    }
    mean /= (double)v->vtotal ;
#if USE_FABS
    mean_error /= (double)v->vtotal ;
#else
    mean_error = sqrt(mean_error / (double)v->vtotal) ;
#endif
#if 0
    if (smean < 0.0f)
      mean_error *= -1.0f ;
#endif
    v->val = mean_error / mean ;
  }

#if USE_FABS
  total_mean_error /= (double)nvertices ;
#else
  total_mean_error = sqrt(total_mean_error / (double)nvertices) ;
#endif
  total_mean /= (double)nvertices ;
  total_mean_error /= total_mean ;
  fprintf(stdout, "mean dist = %2.3f, rms error = %2.2f%%\n",
          total_mean, 100.0*total_mean_error) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteTriangularSurface(MRI_SURFACE *mris, char *fname)
{
  int      k, n ;
  FILE     *fp;
  time_t   tt ;
  char     *user, *time_str ;
  
  user = getenv("USER") ;
  if (!user)
    user = getenv("LOGNAME") ;
  if (!user)
    user = "UNKNOWN" ;
  tt = time(&tt) ;
  time_str = ctime(&tt) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "writing surface file %s, created by %s on %s.\n",
            fname, user, time_str) ;
  fp = fopen(fname,"w");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE,"MRISwrite(%s): can't create file\n",fname));
  fwrite3(TRIANGLE_FILE_MAGIC_NUMBER,fp);
  fprintf(fp, "created by %s on %s\n", user, time_str) ;
  fwriteInt(mris->nvertices,fp);
  fwriteInt(mris->nfaces,fp);   /* # of triangles */
  for (k = 0 ; k < mris->nvertices ; k++)
  {
    fwriteFloat(mris->vertices[k].x, fp) ;
    fwriteFloat(mris->vertices[k].y, fp) ;
    fwriteFloat(mris->vertices[k].z, fp) ;
  }
  for (k = 0 ; k < mris->nfaces ; k++)
  {
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      fwriteInt(mris->faces[k].v[n],fp);
  }
  /* write whether vertex data was using the real RAS rather than conformed RAS */
  fwriteInt(TAG_USEREALRAS, fp);
  fwriteInt(mris->useRealRAS, fp);

  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static SMALL_SURFACE *
mrisReadTriangleFileVertexPositionsOnly(char *fname)
{
  SMALL_VERTEX   *v ;
  int            nvertices, nfaces, magic, vno ;
  char           line[STRLEN] ;
  FILE           *fp ;
  SMALL_SURFACE  *mriss ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL,(ERROR_NOFILE,
                      "mrisReadTriangleFile(%s): could not open file",fname));

  fread3(&magic, fp) ;
  fgets(line, 200, fp) ;
  fscanf(fp, "\n") ;
  /*  fscanf(fp, "\ncreated by %s on %s\n", user, time_str) ;*/
  nvertices = freadInt(fp);
  nfaces = freadInt(fp);  
  
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout,"surface %s: %d vertices and %d faces.\n",
            fname, nvertices,nfaces);
  
  mriss = (SMALL_SURFACE *)calloc(1, sizeof(SMALL_SURFACE)) ;
  if (!mriss)
    ErrorReturn(NULL, 
                (ERROR_NOMEMORY, 
                 "MRISreadVerticesOnly: could not allocate surface")) ;
  
  mriss->nvertices = nvertices ;
  mriss->vertices = (SMALL_VERTEX *)calloc(nvertices, sizeof(SMALL_VERTEX)) ;
  if (!mriss->nvertices)
  {
    free(mriss) ;
    ErrorReturn(NULL, 
                (ERROR_NOMEMORY, 
                 "MRISreadVerticesOnly: could not allocate surface")) ;
  }
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mriss->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->x = freadFloat(fp);
    v->y = freadFloat(fp);
    v->z = freadFloat(fp);
#if 0
    v->label = NO_LABEL ;
#endif
    if (fabs(v->x) > 10000 || !finite(v->x))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d x coordinate %f!",
                Progname, vno, v->x) ;
    if (fabs(v->y) > 10000 || !finite(v->y))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d y coordinate %f!",
                Progname, vno, v->y) ;
    if (fabs(v->z) > 10000 || !finite(v->z))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d z coordinate %f!",
                Progname, vno, v->z) ;
  }
  fclose(fp);
  return(mriss) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisReadTriangleFilePositions(MRI_SURFACE *mris, char *fname)
{
  VERTEX      *v ;
  int         nvertices, nfaces, magic, vno ;
  char        line[STRLEN] ;
  FILE        *fp ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,
                      "mrisReadTriangleFile(%s): could not open file",fname));

  fread3(&magic, fp) ;
  fgets(line, 200, fp) ;
  fscanf(fp, "\n") ;
  /*  fscanf(fp, "\ncreated by %s on %s\n", user, time_str) ;*/
  nvertices = freadInt(fp);
  nfaces = freadInt(fp);  

  if (nvertices != mris->nvertices || nfaces != mris->nfaces)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "mrisReadTriangleFile(%s): surface doesn't match %s\n",
                 fname, mris->fname)) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout,"surface %s: %d vertices and %d faces.\n",
            fname, nvertices,nfaces);
  
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->x = freadFloat(fp);
    v->y = freadFloat(fp);
    v->z = freadFloat(fp);
  }
  
  fclose(fp);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static MRI_SURFACE *
mrisReadTriangleFile(char *fname, double pct_over)
{
  VERTEX      *v ;
  FACE        *f ;
  int         nvertices, nfaces, magic, vno, fno, n ;
  char        line[STRLEN] ;
  FILE        *fp ;
  MRI_SURFACE *mris ;
  int         tag;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL,(ERROR_NOFILE,
                      "mrisReadTriangleFile(%s): could not open file",fname));

  fread3(&magic, fp) ;
  fgets(line, 200, fp) ;
  fscanf(fp, "\n") ;
  /*  fscanf(fp, "\ncreated by %s on %s\n", user, time_str) ;*/
  nvertices = freadInt(fp);
  nfaces = freadInt(fp);  
  
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout,"surface %s: %d vertices and %d faces.\n",
            fname, nvertices,nfaces);
  
  mris = MRISoverAlloc(pct_over*nvertices, pct_over*nfaces,nvertices,nfaces) ;
  mris->type = MRIS_TRIANGULAR_SURFACE ;
  
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->x = freadFloat(fp);
    v->y = freadFloat(fp);
    v->z = freadFloat(fp);
#if 0
    v->label = NO_LABEL ;
#endif
    v->num = 0;   /* will figure it out */
    if (fabs(v->x) > 10000 || !finite(v->x))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d x coordinate %f!",
                Progname, vno, v->x) ;
    if (fabs(v->y) > 10000 || !finite(v->y))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d y coordinate %f!",
                Progname, vno, v->y) ;
    if (fabs(v->z) > 10000 || !finite(v->z))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d z coordinate %f!",
                Progname, vno, v->z) ;
  }
  
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)  
    {
      f->v[n] = freadInt(fp);
      if (f->v[n] >= mris->nvertices || f->v[n] < 0)
        ErrorExit(ERROR_BADFILE, "f[%d]->v[%d] = %d - out of range!\n",
                  fno, n, f->v[n]) ;
    }
    
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      mris->vertices[mris->faces[fno].v[n]].num++;
  }
  // new addition
  if (freadIntEx(&tag, fp))
  {
    if (tag == TAG_USEREALRAS)
      if (!freadIntEx(&mris->useRealRAS,fp)) // set useRealRAS
	mris->useRealRAS = 0; // if error, set to default
  }
  else // no tag found.
  {
    // mark vertex coordinates are using the conformed (256^3) and c_(r,a,s) = 0.
    mris->useRealRAS = 0;
  }
  fclose(fp);
  return(mris) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisRemoveNormalGradientComponent(MRI_SURFACE *mris, int vno)
{
  VERTEX   *v ;
  float    dot ;

  v = &mris->vertices[vno] ;
  if (v->ripflag)
    return(NO_ERROR) ;
  
  dot = v->nx*v->odx + v->ny*v->ody + v->nz*v->odz ;
  v->odx -= dot*v->nx ;
  v->ody -= dot*v->ny ;
  v->odz -= dot*v->nz ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MIN_NBR_DIST  (0.25)

static int
mrisRemoveNeighborGradientComponent(MRI_SURFACE *mris, int vno)
{
  VERTEX   *v, *vn ;
  int      n ;
  float    dx, dy, dz, dot, x, y, z, dist ;

  v = &mris->vertices[vno] ;
  if (v->ripflag)
    return(NO_ERROR) ;
  
  x = v->x ; y = v->y ; z = v->z ; 
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;

    /* too close - take out gradient component in this dir. */
    if (dist <= MIN_NBR_DIST)  
    {
      dx /= dist ; dy /= dist ; dz /= dist ;
      dot = dx*v->odx + dy*v->ody + dz*v->odz ;
      if (dot > 0.0)
      {
        v->odx -= dot*dx ;
        v->ody -= dot*dy ;
        v->odz -= dot*dz ;
      }
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISbuildFileName(MRI_SURFACE *mris, char *sname, char *fname)
{
  char   path[STRLEN], *slash, *dot ;
  
  slash = strchr(sname, '/') ;
  if (!slash)              /* no path - use same one as mris was read from */
  {
    dot = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (dot && (*(dot-1) == 'h') && (*(dot-2) == 'l' || *(dot-2) == 'r'))
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explicitly */
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmodeFilterVals(MRI_SURFACE *mris, int niter)
{
  int    histo[256], i, n, vno, ino, index, max_histo, max_index, nchanged, nzero ;
  VERTEX *v, *vn ;

  MRISclearMarks(mris) ;  /* v->marked = 0 means it hasn't converged yet */
  for (ino  = 0 ; ino < niter ; ino++)
  {
    nzero = nchanged = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->marked)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;

      if (nint(v->val) == 0)
        nzero++ ;
      memset(histo, 0, sizeof(histo)) ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        index = (int)nint(vn->val) ;
        if (index < 0 || index > 255)
          continue ;
        histo[index]++ ;
      }
      max_histo = histo[0] ; max_index = 0 ;
      for (i = 1 ; i < 256 ; i++)
      {
        if (histo[i] > max_histo)
        {
          max_histo = histo[i] ;
          max_index = i ;
        }
      }
      v->valbak = max_index ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->val != v->valbak)
      {
        v->marked = 0 ;  /* process it again */
        nchanged++ ;
      }
      else
        v->marked = 1 ;  /* didn't change */
      v->val = v->valbak ;
    }

    /* unmark all nbrs of unmarked vertices */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 1)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        vn->marked = 0 ;   /* process it again */
      }
    }
    printf("iter %d: %d changed, %d zero\n", ino, nchanged, nzero) ;
    if (!nchanged)
      break ;
  }
  MRISclearMarks(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmodeFilterZeroVals(MRI_SURFACE *mris)
{
  int    histo[256], i, n, vno, ino, index, max_histo, max_index, nchanged, nzero ;
  VERTEX *v, *vn ;

  MRISclearMarks(mris) ;  /* v->marked = 0 means it hasn't converged yet */
  ino = 0 ;
  do
  {
    nzero = nchanged = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->marked)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;

      if (nint(v->val) == 0)
        nzero++ ;
      memset(histo, 0, sizeof(histo)) ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        index = (int)nint(vn->val) ;
        if (index < 0 || index > 255)
          continue ;
        histo[index]++ ;
      }
      max_histo = 0 ; max_index = 0 ;
      for (i = 1 ; i < 256 ; i++)
      {
        if (histo[i] > max_histo)
        {
          max_histo = histo[i] ;
          max_index = i ;
        }
      }
      v->valbak = max_index ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->val != v->valbak)
      {
        v->marked = 0 ;  /* process it again */
        nchanged++ ;
      }
      else
        v->marked = 1 ;  /* didn't change */
      v->val = v->valbak ;
    }

    /* unmark all nbrs of unmarked vertices */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 1)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        vn->marked = 0 ;   /* process it again */
      }
    }
    printf("iter %d: %d changed, %d zero\n", ino++, nchanged, nzero) ;
    if (!nchanged)
      break ;
  } while (nchanged > 0 && nzero > 0) ;
  MRISclearMarks(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_ANNOTATION 1000
extern int annotation_to_index(int annotation) ;
int
MRISmodeFilterAnnotations(MRI_SURFACE *mris, int niter)
{
  int    histo[MAX_ANNOTATION], i, n, vno, ino, index, max_histo, 
         max_annotation, annotations[MAX_ANNOTATION], nchanged = 0 ;
  VERTEX *v, *vn ;

  for (ino  = 0 ; ino < niter ; ino++)
  {
    for (nchanged = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;

      memset(histo, 0, sizeof(histo)) ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        index = annotation_to_index(vn->annotation) ;
        if (index < 0 || index >= MAX_ANNOTATION)
          continue ;
        histo[index]++ ;
        annotations[index] = vn->annotation ;
      }
      index = annotation_to_index(v->annotation) ;
      max_histo = histo[index] ; max_annotation = v->annotation ;
      for (i = 1 ; i < MAX_ANNOTATION ; i++)
      {
        if (histo[i] > max_histo)
        {
          max_histo = histo[i] ;
          max_annotation = annotations[i] ;
        }
      }
      v->undefval = max_annotation ;

    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->annotation != v->undefval)
        nchanged++ ;
      v->annotation = v->undefval ;
    }
    if (nchanged == 0)
      break ;
  }
  printf("%d filter iterations complete (%d requested, %d changed)\n", ino, niter, nchanged) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsoapBubbleVals(MRI_SURFACE *mris, int navgs)
{
  int     vno, n, i, cmpt, nmarked ;
  VERTEX  *v, *vn ;
  double  mean ;

  fprintf(stdout, 
          "performing soap bubble smoothing of vals for %d iterations.\n",
          navgs) ;
  for (i = 0 ; i < navgs ; i++)
  {
    for (nmarked = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || (v->marked==1))
        continue ;

      /* compute average of self and neighbors */
      mean = 0.0 ; cmpt = 0 ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked) 
        {
          mean += vn->val ;
          cmpt++ ;
        }
      }
      if (cmpt>0)   /* had some neighbors with real values */
      {
        v->val = mean / (double)(cmpt) ;
        if (!v->marked)  /* has never been computed before */
          nmarked++ ;
        v->marked = 2 ;  /* has a real value, but is not fixed */
      }
    }
#if 0
    if (!nmarked)
      break ;
#endif
  }
  
  /*  fprintf(stdout, "\n") ;*/
  return(NO_ERROR) ;
}
int
MRISremoveTriangleLinks(MRI_SURFACE *mris)
{
  int    fno, which ;
  FACE   *f ;

  if (!IS_QUADRANGULAR(mris))
    return(NO_ERROR) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "removing non-quadrangular links.\n") ;

  for (fno = 0 ; fno < mris->nfaces ; fno += 2)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    which = WHICH_FACE_SPLIT(f->v[0], f->v[1]) ;
    if (EVEN(which))
    {
      mrisRemoveVertexLink(mris, f->v[1], f->v[2]) ;
      mrisRemoveVertexLink(mris, f->v[2], f->v[1]) ;
    }
    else
    {
      mrisRemoveVertexLink(mris, f->v[0], f->v[2]) ;
      mrisRemoveVertexLink(mris, f->v[2], f->v[0]) ;
    }
  }
  return(NO_ERROR) ;
}
static int
mrisRemoveVertexLink(MRI_SURFACE *mris, int vno1, int vno2)
{
  int    n ;
  VERTEX *v ;

  v = &mris->vertices[vno1] ;
  for (n = 0 ; n < v->vnum ; n++)
    if (v->v[n] == vno2)
      break ;

  if (n < v->vnum)
  {
    memmove(v->v+n, v->v+n+1, (v->vtotal-(n+1))*sizeof(int)) ;
    v->vnum-- ; v->vtotal-- ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISdivideLongEdges(MRI_SURFACE *mris, double thresh)
{
  double   dist ;
  int      vno, nadded, n /*,nvertices, nfaces, nedges, eno*/ ;
  VERTEX   *v, *vn ;
  float    x, y, z ;

  /* make it squared so we don't need sqrts later */
  for (nadded = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    x = v->x ; y = v->y ; z = v->z ;

    /* 
       only add vertices if average neighbor vector is in
       normal direction, that is, if the region is concave or sulcal.
    */
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dist = sqrt(SQR(vn->x-x) + SQR(vn->y - y) + SQR(vn->z - z));
      if (dist > thresh)
      {
        if (mrisDivideEdge(mris, vno, v->v[n]) == NO_ERROR)
          nadded++ ;
      }
    }
  }

  if (Gdiag & DIAG_SHOW && nadded > 0)
  {
    fprintf(stdout, 
            "%2.2f mm: %d vertices added: # of vertices=%d, # of faces=%d.\n", 
            thresh, nadded, mris->nvertices, mris->nfaces) ;
#if 0
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
#endif
  }
  return(nadded) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisAddVertices(MRI_SURFACE *mris, double thresh)
{
  double   dist ;
  int      vno, nadded, n,nvertices, nfaces, nedges, eno ;
  VERTEX   *v, *vn ;
  float    x, y, z ;

  /* make it squared so we don't need sqrts later */
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "dividing edges more than %2.2f mm long.\n", thresh) ;
  for (nadded = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    x = v->origx ; y = v->origy ; z = v->origz ;

    /* 
       only add vertices if average neighbor vector is in
       normal direction, that is, if the region is concave or sulcal.
    */
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dist = sqrt(SQR(vn->origx-x) + SQR(vn->origy - y) + SQR(vn->origz - z));
      if (dist > thresh)
      {
        if (mrisDivideEdge(mris, vno, v->v[n]) == NO_ERROR)
          nadded++ ;
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "%d vertices added: # of vertices=%d, # of faces=%d.\n", 
            nadded, mris->nvertices, mris->nfaces) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
  }
  return(nadded) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_VERTEX_NEIGHBORS 50
#define MAX_FACES            50
static int
mrisDivideEdge(MRI_SURFACE *mris, int vno1, int vno2)
{
  VERTEX   *v1, *v2, *vnew ;
  int      vnew_no, n, m, fno, n1, n2, flist[100] ;
  FACE     *face ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "dividing edge %d --> %d\n", vno1, vno2) ;

  if (vno1 == Gdiag_no || vno2 == Gdiag_no || mris->nvertices == Gdiag_no)
  {
    printf("dividing edge %d --> %d, adding vertex number %d\n",
           vno1, vno2, mris->nvertices) ;
    DiagBreak() ;
  }
  v1 = &mris->vertices[vno1] ;
  v2 = &mris->vertices[vno2] ;
  if (v1->ripflag || v2->ripflag || mris->nvertices >= mris->max_vertices ||
      mris->nfaces >= (mris->max_faces-1))
    return(ERROR_NO_MEMORY) ;

  /* check to make sure these vertices or the faces they are part of
     have enough room to expand.
   */
  if (v1->vnum >= MAX_VERTEX_NEIGHBORS || v2->vnum >= MAX_VERTEX_NEIGHBORS ||
      v1->num  >= MAX_FACES ||     v2->num >= MAX_FACES)
    return(ERROR_NO_MEMORY) ;


  /* add 1 new vertex, 2 new faces, and 2 new edges */
  vnew_no = mris->nvertices ;
  vnew = &mris->vertices[vnew_no] ;
  vnew->x = (v1->x + v2->x) / 2 ;
  vnew->y = (v1->y + v2->y) / 2 ;
  vnew->z = (v1->z + v2->z) / 2 ;
  vnew->tx = (v1->tx + v2->tx) / 2 ;
  vnew->ty = (v1->ty + v2->ty) / 2 ;
  vnew->tz = (v1->tz + v2->tz) / 2 ;

  vnew->infx = (v1->infx + v2->infx) / 2 ;
  vnew->infy = (v1->infy + v2->infy) / 2 ;
  vnew->infz = (v1->infz + v2->infz) / 2 ;

  vnew->pialx = (v1->pialx + v2->pialx) / 2 ;
  vnew->pialy = (v1->pialy + v2->pialy) / 2 ;
  vnew->pialz = (v1->pialz + v2->pialz) / 2 ;

  vnew->cx = (v1->cx + v2->cx) / 2 ;
  vnew->cy = (v1->cy + v2->cy) / 2 ;
  vnew->cz = (v1->cz + v2->cz) / 2 ;
  vnew->x = (v1->x + v2->x) / 2 ;
  vnew->y = (v1->y + v2->y) / 2 ;
  vnew->z = (v1->z + v2->z) / 2 ;
  vnew->odx = (v1->odx + v2->odx) / 2 ;
  vnew->ody = (v1->ody + v2->ody) / 2 ;
  vnew->odz = (v1->odz + v2->odz) / 2 ;
  vnew->val = (v1->val + v2->val) / 2 ;
  vnew->origx = (v1->origx + v2->origx) / 2 ;
  vnew->origy = (v1->origy + v2->origy) / 2 ;
  vnew->origz = (v1->origz + v2->origz) / 2 ;
  vnew->vnum = 2 ;    /* at least connected to two bisected vertices */

  /* count the # of faces that both vertices are part of */
  for (n = 0 ; n < v1->num ; n++)
  {
    fno = v1->f[n] ; face = &mris->faces[fno] ;
    for (m = 0 ; m < VERTICES_PER_FACE ; m++)
      if (face->v[m] == vno2)  
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, " face %d shared.\n", fno) ;
        flist[vnew->num++] = fno ;
        vnew->vnum++ ; vnew->vtotal = vnew->vnum ;
      }
  }

  mris->nvertices++ ;

  /* will be part of two new faces also */
  for (n = 0 ; n < vnew->num ; n++)
    flist[vnew->num+n] = mris->nfaces+n ;
  vnew->num *= 2 ;
#if 0
  flist[vnew->num++] = mris->nfaces ; flist[vnew->num++] = mris->nfaces+1 ;
  vnew->num = 4 ; vnew->vnum = 4 ;
#endif
  vnew->f = (int *)calloc(vnew->num, sizeof(int)) ;
  if (!vnew->f)
    ErrorExit(ERROR_NOMEMORY, "could not allocate %dth face list.\n", vnew_no);
  vnew->n = (uchar *)calloc(vnew->num, sizeof(uchar)) ;
  if (!vnew->n)
    ErrorExit(ERROR_NOMEMORY, "could not allocate %dth face list.\n", vnew_no);
  vnew->v = (int *)calloc(vnew->vnum, sizeof(int)) ;
  if (!vnew->v)
    ErrorExit(ERROR_NOMEMORY, "could not allocate %dth vertex list.\n", 
              vnew_no);

#if 0  
  vnew->v[0] = vno1 ; vnew->v[0] = vno2 ; vnew->vnum = 2 ; vnew->num = 0;
#else
  vnew->num = vnew->vnum = 0 ;
#endif

  /* divide every face that both vertices are part of in two */
  for (n = 0 ; n < v1->num ; n++)
  {
    fno = v1->f[n] ; face = &mris->faces[fno] ;
    for (m = 0 ; m < VERTICES_PER_FACE ; m++)
      if (face->v[m] == vno2)  
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "dividing face %d along edge %d-->%d.\n", 
                  fno, vno1, vno2) ;
        if (face->v[m] == Gdiag_no || vno2 == Gdiag_no)
          DiagBreak() ;
        mrisDivideFace(mris, fno, vno1, vno2, vnew_no) ;
      }
  }

  /* build vnew->f and vnew->n lists by going through all faces v1 and
     v2 are part of */
  for (fno = 0 ; fno < vnew->num ; fno++)
  {
    vnew->f[fno] = flist[fno] ;
    face = &mris->faces[flist[fno]] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      if (face->v[n] == vnew_no)
        vnew->n[fno] = (uchar)n ;
  }

  /* remove vno1 from vno2 list and visa-versa */
  for (n = 0 ; n < v1->vnum ; n++)
    if (v1->v[n] == vno2)
    {
      v1->v[n] = vnew_no ;
      break ;
    }
  for (n = 0 ; n < v2->vnum ; n++)
    if (v2->v[n] == vno1)
    {
      v2->v[n] = vnew_no ;
      break ;
    }
  /* build vnew->v list by going through faces it is part of and
     rejecting duplicates
  */
  for (fno = 0 ; fno < vnew->num ; fno++)
  {
    face = &mris->faces[vnew->f[fno]] ;
    n1 = vnew->n[fno] == 0 ? VERTICES_PER_FACE-1 : vnew->n[fno]-1 ;
    n2 = vnew->n[fno] == VERTICES_PER_FACE-1 ? 0 : vnew->n[fno]+1 ;
    vno1 = face->v[n1] ; vno2 = face->v[n2] ;

    /* go through this faces vertices and see if they should be added to v[] */
    for (n = 0 ; n < vnew->vnum ; n++)
    {
      if (vnew->v[n] == vno1)
        vno1 = -1 ;
      if (vnew->v[n] == vno2)
        vno2 = -1 ;
    }
    if (vno1 >= 0)
      vnew->v[vnew->vnum++] = vno1 ; 
    if (vno2 >= 0)
      vnew->v[vnew->vnum++] = vno2 ;
    vnew->vtotal = vnew->vnum ;
  }
  if (0 && Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "%d edges and %d faces.\n", vnew->vnum, vnew->num) ;

  if (!vnew->vnum || !v1->vnum || !v2->vnum)
  {
    fprintf(stderr, "empty vertex (%d <-- %d --> %d!\n",
            vno1, vnew_no, vno2) ;
    DiagBreak() ;
  }
  if (vnew->vnum != 4 || vnew->num != 4)
    DiagBreak() ;
  mrisInitializeNeighborhood(mris, vnew_no) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisDivideFace(MRI_SURFACE *mris, int fno, int vno1, int vno2, int vnew_no)
{
  FACE   *f1, *f2 ;
  VERTEX *v1, *v2, *v3, *vnew ;
  int    fnew_no, n, vno3, flist[5000], vlist[5000], nlist[5000] ;

  if (vno1 == Gdiag_no || vno2 == Gdiag_no || vnew_no == Gdiag_no)
    DiagBreak() ;

  if (fno ==  10919 || mris->nfaces == 10919)
    DiagBreak() ;

  /* divide this face in two, reusing one of the face indices, and allocating
     one new one
  */
  if (mris->nfaces >= mris->max_faces)
    return(ERROR_NO_MEMORY) ;

  fnew_no = mris->nfaces++ ;
  if (fnew_no == 12452)
    DiagBreak() ;
  f1 = &mris->faces[fno] ;
  f2 = &mris->faces[fnew_no] ;
  v1 = &mris->vertices[vno1] ;
  v2 = &mris->vertices[vno2] ;
  vnew = &mris->vertices[vnew_no] ;
  memmove(f2->v, f1->v, VERTICES_PER_FACE * sizeof(int));
  
  /* set v3 to be other vertex in face being divided */

  /* 1st construct f1 by replacing vno2 with vnew_no */
  for (vno3 = -1, n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    if (f1->v[n] == vno2)   /* replace it with vnew */
    {
      f1->v[n] = vnew_no ;
      vnew->f[vnew->num] = fno ;
      vnew->n[vnew->num++] = (uchar)n ;
    }
    else if (f1->v[n] != vno1)
      vno3 = f1->v[n] ;
  }
  v3 = &mris->vertices[vno3] ;

  if (vno1 == Gdiag_no || vno2 == Gdiag_no || vno3 == Gdiag_no)
    DiagBreak() ;

  /* now construct f2 */

  /*  replace vno1 with vnew_no in f2 */
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    if (f2->v[n] == vno1)   /* replace it with vnew */
    {
      f2->v[n] = vnew_no ;
      vnew->f[vnew->num] = fnew_no ;
      vnew->n[vnew->num++] = (uchar)n ;
    }
  }

  /* now replace f1 in vno2 with f2 */
  for (n = 0 ; n < v2->num ; n++)
    if (v2->f[n] == fno)
      v2->f[n] = fnew_no ;


  /* add new face and edge connected to new vertex to v3 */
  memmove(flist, v3->f, v3->num*sizeof(v3->f[0])) ;
  memmove(vlist, v3->v, v3->vnum*sizeof(v3->v[0])) ;
  memmove(nlist, v3->n, v3->num*sizeof(v3->n[0])) ;
  free(v3->f) ; 
  free(v3->v) ; 
  free(v3->n) ;
  v3->v = (int *)calloc(v3->vnum+1,sizeof(int));
  if (!v3->v)
    ErrorExit(ERROR_NO_MEMORY,"mrisDivideFace: could not allocate %d vertices",
              v3->vnum) ;
  v3->f = (int *)calloc(v3->num+1,sizeof(int));
  if (!v3->f)
    ErrorExit(ERROR_NO_MEMORY, "mrisDivideFace: could not allocate %d faces",
              v3->num) ;
  v3->n = (uchar *)calloc(v3->num+1,sizeof(uchar));
  if (!v3->n)
    ErrorExit(ERROR_NO_MEMORY, "mrisDivideFace: could not allocate %d nbrs",
              v3->n) ;
  memmove(v3->f, flist, v3->num*sizeof(v3->f[0])) ;
  memmove(v3->n, nlist, v3->num*sizeof(v3->n[0])) ;
  memmove(v3->v, vlist, v3->vnum*sizeof(v3->v[0])) ;
  v3->v[v3->vnum++] = vnew_no ; v3->vtotal = v3->vnum ;
  v3->f[v3->num] = fnew_no ;

  /*  find position of v3 in new face f2 */
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    if (f2->v[n] == vno3)   
    {
      v3->n[v3->num] = n ;
      break ;
    }
  }
  if (n >= VERTICES_PER_FACE)
    fprintf(stderr, "could not find v3 (%d) in new face %d!\n",
            vno3, fnew_no) ;
  v3->num++ ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    fprintf(stdout, "face %d: (%d, %d, %d)\n",
            fno, f1->v[0], f1->v[1], f1->v[2]);
    fprintf(stdout, "face %d: (%d, %d, %d)\n",
            fnew_no, f2->v[0], f2->v[1], f2->v[2]);
  }
#if 1
  mrisInitializeNeighborhood(mris, vno3) ;
#else
  /* 
     NOTE!!!!!! This won't work if the neighborhood size is 1
  */
  if (v3->dist)
  {
    memmove(dlist, v3->dist, v3->vtotal*sizeof(v3->dist[0])) ;
    free(v3->dist) ;
    v3->dist = (float *)calloc(v3->vtotal+1, sizeof(float)) ;
    if (!v3->dist)
      ErrorExit(ERROR_NOMEMORY, "mrisDivideFace: could not allocate %d dists",
                v3->vtotal+1) ;
    memmove(v3->dist, dlist, v3->vtotal*sizeof(v3->dist[0])) ;
  }
  if (v3->dist_orig)
  {
    memmove(dlist, v3->dist_orig, v3->vtotal*sizeof(v3->dist_orig[0])) ;
    free(v3->dist_orig) ;
    v3->dist_orig = (float *)calloc(v3->vtotal+1, sizeof(float)) ;
    if (!v3->dist_orig)
      ErrorExit(ERROR_NOMEMORY, 
                "mrisDivideFace: could not allocate %d dist_origs",
                v3->vtotal+1) ;
    memmove(v3->dist_orig, dlist, v3->vtotal*sizeof(v3->dist_orig[0])) ;
  }
#endif
  return(NO_ERROR) ;
}
#if 0

static double mrisComputeVertexNormalSpacingStats(MRI_SURFACE *mris, 
static int mrisTessellateFace(MRI_SURFACE *mris, int fno) ;
static int VertexReplaceNeighbor(VERTEX *v, int vno_old, int vno_new) ;
                                                  double *psigma);

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisAddVertices(MRI_SURFACE *mris, double nsigma)
{
  double   mean, sigma, dist, thresh, dot ;
  int      vno, nadded, n,nvertices, nfaces, nedges, eno, added ;
  VERTEX   *v, *vn ;
  float    x, y, z, sx, sy, sz, nx, ny, nz ;

#if 1
  mean = MRIScomputeVertexSpacingStats(mris, &sigma, NULL, NULL, NULL, NULL) ;
#else
  mean = mrisComputeVertexNormalSpacingStats(mris, &sigma) ;
  thresh *= thresh ;   /* make it squared so we don't need sqrts later */
#endif
  thresh = mean + sigma * nsigma ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "dividing edges more than %2.2f mm long.\n", thresh) ;
  for (nadded = vno = 0 ; vno < mris->nvertices ; vno++)
   {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    /* calculate average normal displacement to neighbors */
    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    x = v->x ;    y = v->y ;   z = v->z ;
    sx = sy = sz = 0.0 ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
      }
    }
    dot = sx*nx+sy*ny+sz*nz;   /* projection onto normal */

    /* 
       only add vertices if average neighbor vector is in
       normal direction, that is, if the region is concave or sulcal.
    */
    for (added = n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dist = sqrt(SQR(vn->x - x) + SQR(vn->y - y) + SQR(vn->z - z)) ;
      if (dist > thresh)
      {
        if (mrisDivideEdge(mris, vno, v->v[n]) == NO_ERROR)
          nadded++ ;
      }
    }

    /* check for sulcal vertices that have asymptoted */
    if (!added && v->marked && dot >= 0.0f)
    {
      dot = v->odx * nx + v->ody * ny + v->odz * nz ;
      if (dot > 0.0f)
      {
        for (n = 0 ; n < v->vnum ; n++)
        {
          vn = &mris->vertices[v->v[n]] ;
          dist = sqrt(SQR(vn->x - x) + SQR(vn->y - y) + SQR(vn->z - z)) ;
          if (dist > mean)
          {
            if (mrisDivideEdge(mris, vno, v->v[n]) == NO_ERROR)
              added++ ;
          }
        }
      }
    }
    nadded += added ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "%d vertices added: # of vertices=%d, # of faces=%d.\n", 
            nadded, mris->nvertices, mris->nfaces) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
  }
  return(mean) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        NOT FINISHED!!!!!!
------------------------------------------------------*/
static  int
mrisTessellateFace(MRI_SURFACE *mris, int fno)
{
  VERTEX  *vc, *v1, *v2, *v, *vnew[3] ;
  int     n, vno1, vno2, vnew_no, nv, fnew_no, vc_no, vlist[15] ;
  float   x, y, z, dx, dy, dz, ox, oy, oz, val ;
  FACE    *f, *fnew ;

  if (mris->nvertices + 4 >= mris->max_vertices)
    return(ERROR_NO_MEMORY) ;
  if (mris->nfaces + 5 >= mris->max_faces)
    return(ERROR_NO_MEMORY) ;

  f = &mris->faces[fno] ;

  /* find centroid of current face and put a vertex there */
  ox = oy = oz = x = y = z = dx = dy = dz = val = 0.0 ;
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    v = &mris->vertices[f->v[n]] ;
    x += v->x ; y += v->y ; z += v->z ;
    dx += v->odx ; dy += v->ody ; dz += v->odz ;
    ox += v->origx ; oy += v->origy ; oz += v->origz ;
    val += v->val ;
  }

  vc_no = mris->nvertices++ ;
  vc = &mris->vertices[vc_no] ;
  vc->val = val / (float)VERTICES_PER_FACE ;
  vc->x = x / (float)VERTICES_PER_FACE ;
  vc->y = y / (float)VERTICES_PER_FACE ;
  vc->z = z / (float)VERTICES_PER_FACE ;
  vc->odx = dx / (float)VERTICES_PER_FACE ;
  vc->ody = dy / (float)VERTICES_PER_FACE ;
  vc->odz = dz / (float)VERTICES_PER_FACE ;
  vc->origx = ox / (float)VERTICES_PER_FACE ;
  vc->origy = oy / (float)VERTICES_PER_FACE ;
  vc->origz = oz / (float)VERTICES_PER_FACE ;
  vc->vnum = 0 ;
  vc->num = 6 ;

  /* now allocate 3 vertices which bisect each edge of the face */
  vlist[0] = f->v[0] ; nv = 1 ;
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    vno1 = f->v[n] ; vno2 = n < VERTICES_PER_FACE-1 ? f->v[n+1] : f->v[0] ;
    v1 = &mris->vertices[vno1] ; v2 = &mris->vertices[vno2] ;
    vnew_no = mris->nvertices++ ;
    v = vnew[n] = &mris->vertices[vnew_no] ;
    v->val = (v1->val+v2->val) / 2.0f ;

    v->x = (v1->x+v2->x) / 2.0f ;
    v->y = (v1->y+v2->y) / 2.0f ;
    v->z = (v1->z+v2->z) / 2.0f ;

    v->odx = (v1->odx+v2->odx) / 2.0f ;
    v->ody = (v1->ody+v2->ody) / 2.0f ;
    v->odz = (v1->odz+v2->odz) / 2.0f ;

    v->origx = (v1->origx+v2->origx) / 2.0f ;
    v->origy = (v1->origy+v2->origy) / 2.0f ;
    v->origz = (v1->origz+v2->origz) / 2.0f ;
    v->num = 0 ;
    VertexReplaceNeighbor(v1, vno2, vnew_no) ;
    VertexReplaceNeighbor(v2, vno1, vnew_no) ;
    vlist[nv++] = vnew_no ; vlist[nv++] = vno2 ;

    /* now build the new vertex's neighbor list */
    v->vnum = 3 ; v->v[0] = vno1 ; v->v[1] = vno2 ; v->v[2] = vc_no ;
    vc->v[vc->vnum++] = vno1 ;
    vc->v[vc->vnum++] = vnew_no ;
    vc->vtotal = vc->vnum ;
  }

  /* 
     at this point all the vertices and edges are in place. Now
     put in new faces, reusing the one we are supertessellating.
  */
  for (n = 0 ; n < nv-1 ; n++)
  {
    if (!n)
      fnew_no = fno ;
    else
      fnew_no = mris->nfaces++ ;
    fnew = &mris->faces[fnew_no] ;
    fnew->v[0] = vlist[n] ;
    fnew->v[1] = vlist[n+1] ;
    fnew->v[2] = vc_no ;
  }

  return(NO_ERROR) ;
}
static int
VertexReplaceNeighbor(VERTEX *v, int vno_old, int vno_new)
{
  int n ;

  for (n = 0 ; n < v->vnum ; n++)
  {
    if (v->v[n] == vno_old)
    {
      v->v[n] = vno_new ;
      break ;
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRIScomputeFaceAreaStats(MRI_SURFACE *mris, double *psigma,
                         double *pmin, double *pmax)
{
  double   total_area, mean, var, nf, sigma, min_area, max_area, area,
           area_scale ;
  int      fno ;
  FACE     *f ;

  MRIScomputeMetricProperties(mris) ;

  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;

  min_area = 1000 ; max_area = -1 ;
  for (var = nf = total_area = 0.0, fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    area = area_scale * f->area ;
    total_area += area ;
    nf++ ;
    var += area*area ;
    if (area > max_area)
      max_area = area ;
    if (area < min_area)
      min_area = area ;
  }
  mean = total_area / nf ;
  *psigma = sigma = sqrt(var / nf - mean*mean) ;
  if (pmin)
    *pmin = min_area ; 
  if (pmax)
    *pmax = max_area ;
  return(mean) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRIScomputeVertexSpacingStats(MRI_SURFACE *mris, double *psigma,
                              double *pmin, double *pmax, int *pvno,int *pvno2)
{
  double   total_dist, mean, var, nv, dist, sigma, min_dist, max_dist,
           dist_scale ;
  int      vno, n ;
  VERTEX   *v, *vn ;

  MRIScomputeMetricProperties(mris) ;
  if (mris->patch)
    dist_scale = 1.0 ;
  else
    dist_scale = sqrt(mris->orig_area / mris->total_area) ;
  dist_scale = 1.0f ;
  min_dist = 1000 ; max_dist = -1 ;
  for (var = nv = total_dist = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      nv++ ;
      dist = sqrt(SQR(vn->x - v->x) + SQR(vn->y - v->y) + SQR(vn->z - v->z));
      dist *= dist_scale ;
      if (dist > max_dist)
      {
        if (pvno)
          *pvno = vno ;
        if (pvno2)
          *pvno2 = v->v[n] ;
        max_dist = dist ;
      }
      if (dist < min_dist)
        min_dist = dist ;
      total_dist += dist ;
      var += dist*dist ;
    }
  }
  mean = total_dist / nv ;
  if (psigma)
    *psigma = sigma = sqrt(var / nv - mean*mean) ;
  if (pmin)
    *pmin = min_dist ; 
  if (pmax)
    *pmax = max_dist ;
  return(mean) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeCurvatureValues(MRI_SURFACE *mris)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->curv > mris->max_curv)
      mris->max_curv = v->curv ;
    if (v->curv < mris->min_curv)
      mris->min_curv = v->curv ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsetVals(MRI_SURFACE *mris, float val)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISscaleVals(MRI_SURFACE *mris, float scale)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val *= scale ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float
MRISdistanceToSurface(MRI_SURFACE *mris, MHT *mht,float x0,float y0,float z0,
                      float nx, float ny, float nz)
{
  double   dist, len ;

#if 0
  {
    FACE   *f = &mris->faces[0] ;
    VERTEX *v0, *v1, *v2 ;

    v0 = &mris->vertices[f->v[0]] ;
    v1 = &mris->vertices[f->v[1]] ;
    v2 = &mris->vertices[f->v[2]] ;
    
    x0 = (v1->x + v0->x) / 2 ;
    y0 = (v1->y + v0->y) / 2 ;
    z0 = (v1->z + v0->z) / 2 ;
    x0 += (v2->x - x0) / 2 ;
    y0 += (v2->y - y0) / 2 ;
    z0 += (v2->z - z0) / 2 ;
    x0 -= f->nx ; y0 -= f->ny ; z0 -= f->nz ;
    nx = f->nx ; ny = f->ny ; nz = f->nz ;
  }
#endif

  len = sqrt(nx*nx+ny*ny+nz*nz) ; nx /= len ; ny /= len ; nz /= len ;
  for (dist = 0.0f ; dist < 128 ; dist += .25)
  {
    if (mrisDirectionTriangleIntersection(mris, x0, y0, z0, 
                                          nx, ny, nz, mht, &dist))
      return(dist) ;
  }
  
  return(0.0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISexpandSurface(MRI_SURFACE *mris, float distance, INTEGRATION_PARMS *parms)
{
  int    vno, n, niter, avgs ;
  VERTEX *v ;
  MHT    *mht = NULL ;

  if (parms == NULL)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->x += distance*v->nx ; 
      v->y += distance*v->ny ; 
      v->z += distance*v->nz ; 
    }
  }
  else
  {
#define EXPANSION_STEP_SIZE 0.25
    if ((parms->write_iterations > 0) && (Gdiag&DIAG_WRITE) && !parms->start_t)
      mrisWriteSnapshot(mris, parms, 0) ;
    mrisClearMomentum(mris) ;
    niter = nint(distance / EXPANSION_STEP_SIZE) ;
    avgs = parms->n_averages ;
    for (n = parms->start_t ; n < parms->start_t+niter ; n++)
    {
      MRIScomputeMetricProperties(mris) ;
      if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
        mht = MHTfillTable(mris, mht) ;
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;
        v->dx = v->nx*EXPANSION_STEP_SIZE ; 
        v->dy = v->ny*EXPANSION_STEP_SIZE ; 
        v->dz = v->nz*EXPANSION_STEP_SIZE ; 
      }
      /*      mrisAverageGradients(mris, avgs) ;*/
      mrisComputeSpringTerm(mris, parms->l_spring) ;
      mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm) ;
      mrisComputeThicknessSmoothnessTerm(mris, parms->l_tsmooth) ;
      mrisComputeNormalSpringTerm(mris, parms->l_nspring) ;
      mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv) ;

      mrisComputeTangentialSpringTerm(mris, parms->l_tspring) ;
      mrisAsynchronousTimeStep(mris, parms->momentum, parms->dt,mht,
                               MAX_ASYNCH_MM) ;
      
      if ((parms->write_iterations > 0) &&
          !((n+1)%parms->write_iterations)&&(Gdiag&DIAG_WRITE))
        mrisWriteSnapshot(mris, parms, n+1) ;
    }
    parms->start_t += n ;
  }
  if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
    MHTfree(&mht) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Translate a surface by (dx, dy, dz) 
------------------------------------------------------*/
int
MRIStranslate(MRI_SURFACE *mris, float dx, float dy, float dz)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->x += dx ;
    v->y += dy ;
    v->z += dz ;
  }
  mrisComputeSurfaceDimensions(mris) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Apply a linear transform (possibly octree) to a surface.
           Note that the LT is given in MRI coordinates, so the
           surface must be transformed into that coordinate system
           before applying the linear transform
------------------------------------------------------*/
int
MRIStransform(MRI_SURFACE *mris, MRI *mri, LTA *lta, MRI *mri_dst)
{
  int    vno ;
  VERTEX *v ;
  Real   xw, yw, zw ;
  MATRIX *m;
  // for ras-to-ras transform
  MATRIX *RASFromSurfaceRAS = 0;
  MATRIX *surfaceRASFromRAS = 0;
  // for voxel-to-voxel transform
  MATRIX *voxelFromSurfaceRAS = 0;
  MATRIX *surfaceRASFromVoxel = 0;
  //
  MATRIX *surfaceRASFromSurfaceRAS = 0;
  LT *lt = 0;
  int srcPresent = 1;
  int dstPresent = 1;
  int error = 0;
  char errMsg[256];
  // depend on the type of transform you have to handle differently
  //
  //       orig              ------>      RAS   c_(ras) != 0
  //        |                              |
  //        |                              | identity
  //        V                              V 
  //    conformed vol        ------>      RAS   c_(ras) != 0
  //        |                              |
  //        | identity                     | surfaceRASFromRAS
  //        V                              V
  //    conformed vol        ------>   surfaceRAS  c_(ras) = 0
  //
  // given a volume transform you have to create a surfaceRAS transform
  //
  // Note that vertices are given by surfaceRAS coordinates
  //
  // RAS-to-RAS transform
  //                  orig                 dst
  //    surfaceRAS--->RAS --(ras-to-ras)-->RAS -->surfaceRAS
  //
  // VOX-to-Vox transform
  //
  //    surfaceRAS--->Vox---(vox-to-vox)-->Vox -->surfaceRAS
  //
  //
  if (lta->num_xforms > 1)
    ErrorExit(ERROR_BADPARM, "we cannot handle multiple transforms\n");
  if (lta->num_xforms == 0)
    ErrorExit(ERROR_BADPARM, "transform does not have transform ;-) \n");

  // if volumes are not given, then try to get them from transform
  lt = &lta->xforms[0];

  // if mri is given, then override the one stored in the transform
  if (!mri && lt->src.valid == 1)
  {
    srcPresent = 0;
    fprintf(stderr, "INFO:try to get src info from transform.\n");
    mri = MRIallocHeader(lt->src.width, lt->src.height, lt->src.depth, MRI_UCHAR);
    mri->x_r = lt->src.x_r; mri->y_r = lt->src.y_r; mri->z_r = lt->src.z_r; mri->c_r = lt->src.c_r;
    mri->x_a = lt->src.x_a; mri->y_a = lt->src.y_a; mri->z_a = lt->src.z_a; mri->c_a = lt->src.c_a;
    mri->x_s = lt->src.x_s; mri->y_s = lt->src.y_s; mri->z_s = lt->src.z_s; mri->c_s = lt->src.c_s;
    mri->xsize = lt->src.xsize; mri->ysize = lt->src.ysize; mri->zsize = lt->src.zsize;
    mri->ras_good_flag = 1;
  }
  else if (!mri)
  {
    error = 1;
    strcpy(errMsg, "When mri == NULL, the transform must have the valid src info.\n");
    goto mristransform_cleanup;
  }
  // mri_dst side
  // Note: if mri_dst is given, override the one stored in the transform
  if (!mri_dst && lt->dst.valid == 1)
  {
    dstPresent = 0;
    fprintf(stderr, "INFO:try to get dst info from transform.\n");
    lt = &lta->xforms[0];
    mri_dst = MRIallocHeader(lt->dst.width, lt->dst.height, lt->dst.depth, MRI_UCHAR);
    mri_dst->x_r = lt->dst.x_r; mri_dst->y_r = lt->dst.y_r; mri_dst->z_r = lt->dst.z_r; 
    mri_dst->c_r = lt->dst.c_r;
    mri_dst->x_a = lt->dst.x_a; mri_dst->y_a = lt->dst.y_a; mri_dst->z_a = lt->dst.z_a; 
    mri_dst->c_a = lt->dst.c_a;
    mri_dst->x_s = lt->dst.x_s; mri_dst->y_s = lt->dst.y_s; mri_dst->z_s = lt->dst.z_s; 
    mri_dst->c_s = lt->dst.c_s;
    mri_dst->xsize = lt->dst.xsize; mri_dst->ysize = lt->dst.ysize; mri_dst->zsize = lt->dst.zsize;
    mri_dst->ras_good_flag = 1;
  }
  else if (!mri_dst)
  {
    error = 1;
    strcpy(errMsg, "INFO:When mri_dst == NULL, the transform must have the valid dst info.\n");
    strcpy(errMsg, "INFO:If your target is average_305 and the transform is RAS-to-RAS,\n");
    strcpy(errMsg, "INFO:then you can set environmental variable USE_AVERAGE305 to be true\n");
    strcpy(errMsg, "INFO:and try again.\n");
    goto mristransform_cleanup;
  }
  /////////////////////////////////////////////////////////////////////////////
  // Now we can calculate
  if (lta->type == LINEAR_RAS_TO_RAS)
  {
    RASFromSurfaceRAS = RASFromSurfaceRAS_(mri); // needs only c_(ras) info
    surfaceRASFromRAS = surfaceRASFromRAS_(mri_dst);  // need only c_(ras) info
    m = MatrixMultiply(lta->xforms[0].m_L, RASFromSurfaceRAS, NULL);
    surfaceRASFromSurfaceRAS = MatrixMultiply(surfaceRASFromRAS, m, NULL);
  }
  else if (lta->type == LINEAR_VOX_TO_VOX)
  {
    if (mri->width != mri_dst->width 
	|| mri->height != mri_dst->height
	|| mri->depth != mri_dst->depth)
    {
      fprintf(stderr, "WARNING:********************************************************\n");
      fprintf(stderr, "WARNING:voxel-to-voxel transform must have the same volume sizes.\n");
      fprintf(stderr, "WARNING:You gave src (%d, %dm, %d) vs. dst (%d, %d, %d).\n",
	      mri->width, mri->height, mri->depth, mri_dst->width, mri_dst->height, mri_dst->depth);
      fprintf(stderr, "WARNING:The result of this transform is most likely wrong.\n");
      fprintf(stderr, "WARNING:********************************************************\n");
    }
    voxelFromSurfaceRAS = voxelFromSurfaceRAS_(mri);
    surfaceRASFromVoxel = surfaceRASFromVoxel_(mri_dst);
    m = MatrixMultiply(lta->xforms[0].m_L, voxelFromSurfaceRAS, NULL);
    surfaceRASFromSurfaceRAS = MatrixMultiply(surfaceRASFromVoxel, m, NULL);
  }
  // now apply the transform
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    // transform vertex point to actual voxel point
    TransformWithMatrix(surfaceRASFromSurfaceRAS, v->x, v->y, v->z, &xw, &yw, &zw);
    v->x = xw ; v->y = yw ; v->z = zw ;
  }
  mrisComputeSurfaceDimensions(mris) ;

 mristransform_cleanup:
  // free memory //////////////////////////////
  if (RASFromSurfaceRAS)
    MatrixFree(&RASFromSurfaceRAS);
  if (surfaceRASFromRAS)
    MatrixFree(&surfaceRASFromRAS);
  if (m)
    MatrixFree(&m);
  if (voxelFromSurfaceRAS)
    MatrixFree(&voxelFromSurfaceRAS);
  if (surfaceRASFromVoxel)
    MatrixFree(&surfaceRASFromVoxel);
  if (surfaceRASFromSurfaceRAS)
    MatrixFree(&surfaceRASFromSurfaceRAS);
  if (!srcPresent && mri)
    MRIfree(&mri);
  if (!dstPresent && mri_dst)
    MRIfree(&mri_dst);

  if (error)
  {
    ErrorExit(ERROR_BADPARM, errMsg);
    return -1; // just to satisfy compiler
  }
  else
    return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale a surface anisotropically.
------------------------------------------------------*/
int
MRISanisotropicScale(MRI_SURFACE *mris, float sx, float sy, float sz)
{
  VERTEX  *v;
  int     k;
  float   x0, y0, z0 ;

  mrisComputeSurfaceDimensions(mris) ;
  /* scale around the center */
  x0 = mris->xctr ; y0 = mris->yctr ; z0 = mris->zctr ;

  for (k=0;k<mris->nvertices;k++) 
  {
    v = &mris->vertices[k];
    if (v->ripflag)
      continue ;
    v->x = (v->x - x0) * sx + x0 ; 
    v->y = (v->y - y0) * sy + y0 ; 
    v->z = (v->z - z0) * sz + z0 ;
  }

  mrisComputeSurfaceDimensions(mris) ;
  return(NO_ERROR) ;
}
static int mrisChooseFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)  ;
static int mrisFindUnambiguousFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v, 
                                   int *pnum) ;

#if 1
int
MRISinverseSphericalMap(MRI_SURFACE *mris, MRI_SURFACE *mris_ico)
{
  double r ;
  int    fno, vno, num_ambiguous = 0, nfound, i, max_count ;
  VERTEX *v ;
  MHT    *mht ;
  short  *vcount ;

  vcount = (short *)calloc(mris->nfaces, sizeof(short)) ;

  MRISfreeDists(mris) ;

  /* make sure they are they same size */
  r = MRISaverageRadius(mris_ico) ; 
  MRISscaleBrain(mris_ico,mris_ico, 100.0/r);
  r = MRISaverageRadius(mris) ;
  MRISscaleBrain(mris, mris, 100.0/r);
  MRISstoreMetricProperties(mris) ;

  /*
    orig       positions are on cortical surface
    current    positions are on sphere.
  */
  mht = MHTfillTable(mris, NULL) ;

  /*
    for each vertex on the icosahedral surface, find what face it lies
    in on the spherical representation of the cortex. If it only lies
    within a single face, position it at the centroid of the original
    position of the face, and mark it as positioned.
  */
  MRISclearMarks(mris_ico) ; MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    v = &mris_ico->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    fno = mrisFindUnambiguousFace(mris, mht, v, &nfound) ;
    if (fno >= 0)
    {
      vcount[fno]++ ;
#if 0
      mrisCalculateOriginalFaceCentroid(mris,fno,
                                        &v->origx,&v->origy,&v->origz);
#else
      mrisPlaceVertexInOrigFace(mris, v, fno) ;
#endif
      if (vno == Gdiag_no)
      {
        fprintf(stdout, "vertex %d maps to face %d at (%2.1f, %2.1f, %2.1f)\n",
                vno, fno, v->origx, v->origy, v->origz) ;
        mrisDumpFace(mris, fno, stderr) ;
      }
      v->marked = 1 ;
    }
    else
    {
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "v %d maps to %d faces\n", vno, nfound) ;
      num_ambiguous++ ;
    }
  }
  fprintf(stdout, "%d non-invertible locations found - resolving ambiguity\n",
          num_ambiguous) ;
  
  MRISsoapBubbleOrigVertexPositions(mris_ico, 100) ;
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    v = &mris_ico->vertices[vno] ;
    if (v->marked || v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    fno = mrisChooseFace(mris, mht, v) ;
    if (fno < 0)
      ErrorPrintf(ERROR_BADPARM, "unable to find face for ico vertex %d!!!\n",
                  vno) ;
    else
    {
#if 0
      mrisCalculateOriginalFaceCentroid(mris, fno, &v->x, &v->y, &v->z) ;
#else
      mrisPlaceVertexInOrigFace(mris, v, fno) ;
#endif
      vcount[fno]++ ;
    }
  }
  for (max_count = i = 0 ; i < 50 ; i++)
  {
    for (nfound = fno = 0 ; fno < mris->nfaces ; fno++)
    {
      if (vcount[fno] == i)
        nfound++ ;
    }
    if (nfound)
    {
      if (i > max_count)
        max_count = i ;
      fprintf(stdout, "%d mappings to a single face %d times.\n", i, nfound) ;
    }
  }

  fprintf(stdout, "faces mapped to %d times: \n", max_count) ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    if (vcount[fno] == max_count)
      fprintf(stdout, "\t%d (%d, %d, %d)\n",
              fno, mris->faces[fno].v[0],mris->faces[fno].v[1],
              mris->faces[fno].v[2]) ;
  }
  MHTfree(&mht) ;
  free(vcount) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Given an intersection point int_pt in which a vertex intersects
          a face, compute the analagous location that the vertex should be
          placed in the faces original coordnates.

          note that fno refers to a face in mris and
          v is a vertex in mris_ico (not given), NOT mris.
------------------------------------------------------*/
static int
mrisPlaceVertexInOrigFace(MRI_SURFACE *mris, VERTEX *v, int fno)
{
  double  U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3], l1[3], l2[3], l2_len,
          l1_len, l_len, P[3], theta1, theta2, dot, theta_ratio, len_scale,
          e1[3], e2[3], etmp[3], x, y ;
  int     ret ;

  /* first compute point where normal to vertex intersects face */
  dir[0] = v->nx ; dir[1] = v->ny ; dir[2] = v->nz ;
  pt[0] = v->x  ; pt[1] = v->y  ; pt[2] = v->z  ;
  load_triangle_vertices(mris, fno, U0, U1, U2) ;
  ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt) ;
  if (ret == 0)   /* try in negative of normal direction */
  {
    dir[0] = -v->nx ; dir[1] = -v->ny ; dir[2] = -v->nz ;
    pt[0] = v->x  ; pt[1] = v->y  ; pt[2] = v->z  ;
    load_triangle_vertices(mris, fno, U0, U1, U2) ;
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt) ;
    if (ret == 0)
      ErrorReturn(ERROR_BADPARM, 
                  (ERROR_BADPARM, 
                   "mrisPlaceVertexInOrigFace: v does not intersect face!")) ;
  }

  /* now normalize the edges (l1 and l2) of the current triangle */
  SUB(l1, U1, U0) ; SUB(l2, U2, U0) ; SUB(P, int_pt, U0) ;
  l1_len = VLEN(l1) ; SCALAR_MUL(l1, 1.0/l1_len, l1) ;
  l2_len = VLEN(l2) ; SCALAR_MUL(l2, 1.0/l2_len, l2) ;
  l_len =  VLEN(P) ;  SCALAR_MUL(P,  1.0/l_len, P) ;

  /*
    compute the angle between the two legs, and between P and l1.
    The ratio of these two angles will be used to place the vertex
    in the original triangle.
  */
  dot = DOT(l1, P) ;  theta1 = acos(dot) ;
  if (theta1 < 0)
    theta1 += 2*PI ;
  dot = DOT(l1, l2) ; theta2 = acos(dot) ;
  if (theta2 < 0)
    theta2 += 2*PI ;
  if (!DZERO(theta2))  
    theta_ratio = theta1 / theta2 ;
  else   /* degenerate triangle */
    theta_ratio = 0 ;

  /*
    express the ratio of the length of the line segment P-U0 as
    a scaled linear combination of the l1 and l2 (the legs of the
    triangle), where the relative weighting is based on the theta
    ratio. This will allow us to use the ratio and the original
    lengths of the legs to place the point in the corresponding location
    in the original triangle.
  */
  len_scale = l_len / (((1-theta_ratio)*l1_len)+theta_ratio*l2_len) ;
  load_orig_triangle_vertices(mris, fno, U0, U1, U2) ;
  SUB(l1, U1, U0) ; SUB(l2, U2, U0) ;
  l1_len = VLEN(l1) ; SCALAR_MUL(l1, 1.0/l1_len, l1) ;
  l2_len = VLEN(l2) ; SCALAR_MUL(l2, 1.0/l2_len, l2) ;
  l_len =  VLEN(P) ;  SCALAR_MUL(P,  1.0/l_len,  P) ;

  /* compute angle between original legs */
  dot = DOT(l1, l2) ; theta2 = acos(dot) ;
  if (theta2 < 0)
    theta2 += 2*PI ;

  theta1 = theta_ratio * theta2 ;   /* analogous angle in orig triangle */
  
  /* construct basis vector for plane defined by original triangle */
  SCALAR_MUL(l1, 1.0, e1) ;  /* 1st basis vector is just l1 */
  CROSS(l1, l2, etmp) ;      /* vector orthogonal to l1 and l2 */
  CROSS(etmp, l1, e2) ;
  SCALAR_MUL(e2, (1.0/VLEN(e2)), e2) ;

  /* 
     express length of line segment in original triangle as a linear
     combination of the original leg lengths, using the same weighting.
  */
  l_len = len_scale * (((1-theta_ratio)*l1_len)+theta_ratio*l2_len) ;

  /* rotate e1 by theta1 and scale it by l_len */
  x = l_len * cos(theta1) ;
  y = l_len * sin(theta1) ;

  /* express it in the global coordinate system */
  SCALAR_MUL(e1, x, e1) ; SCALAR_MUL(e2, y, e2) ; ADD(e1, e2, P) ;
  ADD(P, U0, P) ;
  
  v->origx = P[0] ; v->origy = P[1] ; v->origz = P[2] ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Find the face in which v lies. If it lies in more than one
          face, return not found.
------------------------------------------------------*/
static int
mrisFindUnambiguousFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v, int *pnfound)
{
  int     nfound, flist[1000], *fptr, total_found, i, j ;
  double  dist, d ;

  for (total_found = nfound = 0, dist = -.25 ; dist <=.25 ; dist += .25)
  {
    d = dist ; fptr = &flist[total_found] ;
    nfound = 
      mrisAllNormalDirectionCurrentTriangleIntersections(mris, v,mht,&d,fptr);
    if (nfound > 0)
    {
      for (i = 0 ; i < total_found ; i++)
      {
        for (j = 0 ; j < nfound ; j++)
        {
          if (flist[i] == fptr[j])
            fptr[j] = -1 ;   /* was already found */
        }
      }
      for (j = 0 ; j < nfound ; j++)
      {
        if (fptr[j] >= 0)
          flist[total_found++] = fptr[j] ;
      }
    }
  }

  *pnfound = total_found ;
  return(total_found == 1 ? flist[0] : -1) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Choose the face whose centroid is clostest to the orig position 
          of v.
------------------------------------------------------*/
#if 1
static int 
mrisChooseFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)
{
  int     fno, nfound, flist[1000], min_fno, i, j, total_found, *fptr ;
  double  dist, d ;
  float   dx, dy, dz, cx[1000], cy[1000], cz[1000], total_dist,max_dist;

  for (total_found = nfound = 0, dist = -.25 ; dist <=.25 ; dist += .25)
  {
    d = dist ; fptr = &flist[total_found] ;
    nfound = 
      mrisAllNormalDirectionCurrentTriangleIntersections(mris, v,mht,&d,fptr);
    if (nfound > 0)
    {
      for (i = 0 ; i < total_found ; i++)
      {
        for (j = 0 ; j < nfound ; j++)
        {
          if (flist[i] == fptr[j])
            fptr[j] = -1 ;   /* was already found */
        }
      }
      for (j = 0 ; j < nfound ; j++)
      {
        if (fptr[j] >= 0)
          flist[total_found++] = fptr[j] ;
      }
    }
  }

  if (total_found <= 0)
    return(-1) ;

  /*
    use face which is furthest distance from negative faces.
  */
  max_dist = -10.0f ; min_fno = -1 ;
  
  for (i = 0 ; i < total_found ; i++)
  {
    fno = flist[i] ;
    mrisCalculateOriginalFaceCentroid(mris, fno, cx+i, cy+i, cz+i) ;
  }
  for (i = 0 ; i < total_found ; i++)
  {
    fno = flist[i] ;
    if (mris->faces[fno].area < 0)
      continue ;   /* don't map it to one with negative area */
    
    for (total_dist = 0.0, j = 0 ; j < total_found ; j++)
    {
      if (mris->faces[flist[j]].area > 0)
        continue ;
      dx = cx[j] - cx[i] ; dy = cy[j] - cy[i] ; dz = cz[j] - cz[i] ;
      dist = sqrt(dx*dx + dy*dy + dz*dz) ;
      total_dist += dist ;
    }
    if (total_dist > max_dist)
    {
      max_dist = dist ;
      min_fno = fno ;
    }
  }
  return(min_fno) ;
}
#else
static int 
mrisChooseFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)
{
  int     fno, nfound, flist[1000], min_fno, i, j, total_found, *fptr ;
  float   x, y, z, dx, dy, dz ;
  double  dist, min_dist, d ;

  for (total_found = nfound = 0, dist = -.25 ; dist <=.25 ; dist += .25)
  {
    d = dist ; fptr = &flist[total_found] ;
    nfound = 
      mrisAllNormalDirectionCurrentTriangleIntersections(mris, v,mht,&d,fptr);
    if (nfound > 0)
    {
      for (i = 0 ; i < total_found ; i++)
      {
        for (j = 0 ; j < nfound ; j++)
        {
          if (flist[i] == fptr[j])
            fptr[j] = -1 ;   /* was already found */
        }
      }
      for (j = 0 ; j < nfound ; j++)
      {
        if (fptr[j] >= 0)
          flist[total_found++] = fptr[j] ;
      }
    }
  }

  if (total_found <= 0)
    return(-1) ;
  min_dist = 10000.0f ; min_fno = -1 ;

  for (i = 0 ; i < total_found ; i++)
  {
    fno = flist[i] ;
    mrisCalculateOriginalFaceCentroid(mris, fno, &x, &y, &z) ;
    dx = x - v->origx ; dy = y - v->origy ; dz = z - v->origz ; 
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (dist < min_dist)
    {
      min_dist = dist ;
      min_fno = fno ;
    }
  }
  return(min_fno) ;
}
#endif
#else
int
MRISinverseSphericalMap(MRI_SURFACE *mris, MRI_SURFACE *mris_ico)
{
  double r ;
  int    vno /*, fno, n, num*/ ;
#if 0
  VERTEX *v, *vn ;
  MHT    *mht ;
  float  x, y, z ;
#endif
  float orient ;
  static char ***flagvol, ***numvol ; /* volume of flags for neg area test */
        static float ***xvol, ***yvol, ***zvol ;
        static int allocated = 0;
  float flagvolres = 2.0, flagvolfov = 300;
  int flagvoldim,i,j;
        int ix, iy, iz ;

  MRISclearCurvature(mris) ;   /* curvature will be used to calculate sulc */
  r = MRISaverageRadius(mris) ; MRISscaleBrain(mris, mris, 100.0f/r) ;
/*
  r = MRISaverageRadius(mris_ico) ; MRISscaleBrain(mris_ico,mris_ico,100.0f/r);
*/

  /*
    orig       positions are on orig
    cx,cy,cz   positions are on inflated surface
    current    positions are on sphere.
  */
/*
  printf("begin MHTfillTableWithFaces\n");

  mht = MHTfillTableWithFaces(mris, NULL) ;

  printf("end MHTfillTableWithFaces\n");
*/

        flagvoldim = ceil(flagvolfov/flagvolres);
        if (!allocated)
        {
          fprintf(stdout, "allocating flagvol...\n") ;
          flagvol = (char ***)calloc(flagvoldim, sizeof(char **));
          numvol = (char ***)calloc(flagvoldim, sizeof(char **));
          xvol = (float ***)calloc(flagvoldim, sizeof(float **));
          yvol = (float ***)calloc(flagvoldim, sizeof(float **));
          zvol = (float ***)calloc(flagvoldim, sizeof(float **));
          for (i=0;i<flagvoldim;i++)
          {
            flagvol[i] = (char **)calloc(flagvoldim, sizeof(char *));
            numvol[i] = (char **)calloc(flagvoldim, sizeof(char *));
            xvol[i] = (float **)calloc(flagvoldim, sizeof(float *));
            yvol[i] = (float **)calloc(flagvoldim, sizeof(float *));
            zvol[i] = (float **)calloc(flagvoldim, sizeof(float *));
            for (j=0;j<flagvoldim;j++)
            {
              flagvol[i][j] = (char *)calloc(flagvoldim, sizeof(char));
              numvol[i][j] = (char *)calloc(flagvoldim, sizeof(char));
              xvol[i][j] = (float *)calloc(flagvoldim, sizeof(float));
              yvol[i][j] = (float *)calloc(flagvoldim, sizeof(float));
              zvol[i][j] = (float *)calloc(flagvoldim, sizeof(float));
            }
          }
          allocated = 1;
        }
  
        for (ix=0;ix<flagvoldim;ix++)
        for (iy=0;iy<flagvoldim;iy++)
        for (iz=0;iz<flagvoldim;iz++)
          numvol[ix][iy][iz] = xvol[ix][iy][iz] = yvol[ix][iy][iz] = zvol[ix][iy][iz] = flagvol[ix][iy][iz] = 0;

        for (vno = 0 ; vno < mris->nvertices ; vno++)
        {
          orient = mris->vertices[vno].x*mris->vertices[vno].nx+
                   mris->vertices[vno].y*mris->vertices[vno].ny+
                   mris->vertices[vno].z*mris->vertices[vno].nz ;
          if (orient < 0)
          {
            printf("vertex %d inside out (orient = %f)\n",vno,orient);
/*
            printf("x = (%3.1f, %3.1f,%3.1f) n = (%3.1f,%3.1f,%3.1f)\n",
                    mris->vertices[vno].x,mris->vertices[vno].y,mris->vertices[vno].z,
                    mris->vertices[vno].nx,mris->vertices[vno].ny,mris->vertices[vno].nz);
*/
            mris->vertices[vno].curv = 1;
          }
          else
            mris->vertices[vno].curv = 0;
        }

        for (vno = 0 ; vno < mris->nvertices ; vno++)
        {
          ix = floor(0.5+(mris->vertices[vno].x+flagvolfov/2)/flagvolres);
          iy = floor(0.5+(mris->vertices[vno].y+flagvolfov/2)/flagvolres);
          iz = floor(0.5+(mris->vertices[vno].z+flagvolfov/2)/flagvolres);
          numvol[ix][iy][iz]++;
          xvol[ix][iy][iz] += mris->vertices[vno].cx;
          yvol[ix][iy][iz] += mris->vertices[vno].cy;
          zvol[ix][iy][iz] += mris->vertices[vno].cz;
          if (mris->vertices[vno].curv != 0) /* inverted vertex? */
            flagvol[ix][iy][iz] = mris->vertices[vno].curv;
        }

        for (ix=0;ix<flagvoldim;ix++)
        for (iy=0;iy<flagvoldim;iy++)
        for (iz=0;iz<flagvoldim;iz++)
        if (numvol[ix][iy][iz]>0)
        {
          xvol[ix][iy][iz] /= numvol[ix][iy][iz];
          yvol[ix][iy][iz] /= numvol[ix][iy][iz];
          zvol[ix][iy][iz] /= numvol[ix][iy][iz];
        }

        /* restore orig. vertex coords */
        for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
        {
          mris_ico->vertices[vno].x = mris_ico->vertices[vno].origx;
          mris_ico->vertices[vno].y = mris_ico->vertices[vno].origy;
          mris_ico->vertices[vno].z = mris_ico->vertices[vno].origz;
        }

        for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
        {
          ix = floor(0.5+(mris_ico->vertices[vno].x+flagvolfov/2)/flagvolres);
          iy = floor(0.5+(mris_ico->vertices[vno].y+flagvolfov/2)/flagvolres);
          iz = floor(0.5+(mris_ico->vertices[vno].z+flagvolfov/2)/flagvolres);
            mris_ico->vertices[vno].curv = flagvol[ix][iy][iz];
          if (mris_ico->vertices[vno].curv != 0)
          {
            mris_ico->vertices[vno].marked = 0;
            printf("ambiguous ico vertex %d\n",vno);
          }
          else
            mris_ico->vertices[vno].marked = 1;
          if (numvol[ix][iy][iz]>0)
          {
            mris_ico->vertices[vno].x = xvol[ix][iy][iz];
            mris_ico->vertices[vno].y = yvol[ix][iy][iz];
            mris_ico->vertices[vno].z = zvol[ix][iy][iz];
          }
          else
          {
            printf("### ico vertex %d missed volume\n",vno);
            mris_ico->vertices[vno].marked = 0;
          }
        }


  MRISsoapBubbleVertexPositions(mris_ico, 100) ;

/*
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    v = &mris_ico->vertices[vno] ;
    if (v->marked || v->ripflag)
      continue ;
    fno = mrisChooseFace(mris, mht, v) ;
    mrisCalculateCanonicalFaceCentroid(mris, fno, &v->x, &v->y, &v->z) ;
  }
  MHTfree(&mht) ;
*/
  return(NO_ERROR) ;
}
#endif
#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define VERTEX_DIF(leg, v0, v1)   leg[0] = v1->x-v0->x, \
                                  leg[1] = v1->y-v0->y,\
                                  leg[2] = v1->z-v0->z ;

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#if 0
static float mrisComputeFaceStretch(MRI_SURFACE *mris, int fno) ;
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        mris is the original (hi-res) surface (from smoothwm).
        v is a vertex on the icosahedron.
------------------------------------------------------*/
static int
mrisChooseFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)
{
  int    i, /*fno,*/ flist[10000], nfaces, min_fno = 0, l ;
  MHBT   *bucket, bucket_bak  ;
  MHB    *bin ;
  FACE   *f ;
  VERTEX *v1, *v2, *v3 ;
  float  /*ut, vt, stretch,*/ n[3], leg[3], p[3], d[3], dot, ldot, leg2[3] ;

  bucket = MHTgetBucket(mht, v->x, v->y, v->z) ;
  bucket_bak = *bucket ;
  bin = bucket->bins ;
  for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)  /* check each face */
  {
    n[0] = n[1] = n[2] = 0.0f ;
    f = &mris->faces[bin->fno] ;
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      n[0] += v1->nx ; n[1] += v1->ny ; n[2] += v1->nz ; 
    }
    n[0] /= VERTICES_PER_FACE ; n[1] /= VERTICES_PER_FACE ; 
    n[2] /= VERTICES_PER_FACE ; 
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      switch (l)
      {
      default:
      case 0:
        v2 = &mris->vertices[f->v[1]] ;
        v3 = &mris->vertices[f->v[2]] ;
        break ;
      case 1:
        v2 = &mris->vertices[f->v[2]] ;
        v3 = &mris->vertices[f->v[0]] ;
        break ;
      case 2:
        v2 = &mris->vertices[f->v[0]] ;
        v3 = &mris->vertices[f->v[1]] ;
        break ;
      }

      VERTEX_DIF(leg, v1, v2) ;   /* leg of triangle */
      VERTEX_DIF(leg2, v1, v3) ;  /* leg of triangle */

      /* express p as point in triangle plane */
      VERTEX_DIF(p, v1, v) ;     /* vector from vertex to point in question */
      dot = DOT(p,n) ;
      p[0] -= dot*n[0] ; p[1] -= dot*n[1] ; p[2] -= dot*n[2] ; 
#if 0
      p[0] = ut*leg[0] + vt*leg2[0] ;
      p[1] = ut*leg[1] + vt*leg2[1] ;
      p[2] = ut*leg[2] + vt*leg2[2] ;
      ut = DOT(p, leg) ; vt = DOT(p, leg2) ;
#endif

      CROSS(d, leg, n) ;
      dot = DOT(d, p) ; ldot = DOT(d, leg2) ;

      /* on different side of leg from 3rd vertex */
      if (!FZERO(ldot) && !FZERO(dot) && dot*ldot < 0)
        break ;
    }
    if (l >= VERTICES_PER_FACE)
      flist[nfaces++] = bin->fno ;
  }
  if (!nfaces)  /* something went wrong, but Anders will fix it */
  {
    float dist, min_dist ;

    fprintf(stderr, "no faces found on sphere!\n") ;
    bin = bucket->bins ;
    min_dist = 1000000.0f ; min_fno = 0 ;
    for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      f = &mris->faces[bin->fno] ;
      v1 = &mris->vertices[f->v[0]] ;
      v2 = &mris->vertices[f->v[1]] ;
      v3 = &mris->vertices[f->v[2]] ;
#define VDIST(v1,v2) (sqrt(SQR(v1->x-v2->x)+SQR(v1->y-v2->y)+SQR(v1->z-v2->z)))
      dist = VDIST(v1,v) + VDIST(v2, v) + VDIST(v3,v) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = bin->fno ;
      }
    }
  }
  else   /* pix the face that is closest to the soap bubble location */
  {
    float min_dist, x, y, z, dist, curv ;

    if (nfaces > 1)
    {
      DiagBreak() ;
      curv = 1.0f ;
    }
    else
      curv = 0.0f ;
    for ( i = 0 ; i < nfaces ; i++)/* check each face */
    {
      f = &mris->faces[flist[i]] ;
      for (l = 0 ; l < VERTICES_PER_FACE ; l++)
        mris->vertices[f->v[l]].curv = curv ;
    }

    min_dist = 100000.0f ;
    for (i = 0 ; i < nfaces ; i++)
    {
      mrisCalculateCanonicalFaceCentroid(mris, flist[i], &x, &y, &z) ;
      dist = sqrt(SQR(v->x-x)+SQR(v->y-y)+SQR(v->z-z)) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = flist[i] ;
      }
    }
  }
  return(min_fno) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale a surface anisotropically.
------------------------------------------------------*/
static int
mrisFindUnambiguousFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)
{
  int    i, /*fno,*/ flist[10000], nfaces, min_fno = 0, l ;
  MHBT   *bucket ;
  MHB    *bin ;
  FACE   *f ;
  VERTEX *v1, *v2, *v3 ;
  float  /*ut, vt, */stretch, n[3], leg[3], p[3], d[3], dot, ldot, leg2[3] ;

  bucket = MHTgetBucket(mht, v->x, v->y, v->z) ;
  bin = bucket->bins ;
  for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)  /* check each face */
  {
    n[0] = n[1] = n[2] = 0.0f ;
    f = &mris->faces[bin->fno] ;
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      n[0] += v1->nx ; n[1] += v1->ny ; n[2] += v1->nz ; 
    }
    n[0] /= VERTICES_PER_FACE ; n[1] /= VERTICES_PER_FACE ; 
    n[2] /= VERTICES_PER_FACE ; 
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      switch (l)
      {
      default:
      case 0:
        v2 = &mris->vertices[f->v[1]] ;
        v3 = &mris->vertices[f->v[2]] ;
        break ;
      case 1:
        v2 = &mris->vertices[f->v[2]] ;
        v3 = &mris->vertices[f->v[0]] ;
        break ;
      case 2:
        v2 = &mris->vertices[f->v[0]] ;
        v3 = &mris->vertices[f->v[1]] ;
        break ;
      }

      VERTEX_DIF(leg, v1, v2) ;   /* leg of triangle */
      VERTEX_DIF(leg2, v1, v3) ;  /* leg of triangle */

      /* express p as point in triangle plane */
      VERTEX_DIF(p, v1, v) ;     /* vector from vertex to point in question */
      dot = DOT(p,n) ;
      p[0] -= dot*n[0] ; p[1] -= dot*n[1] ; p[2] -= dot*n[2] ; 
#if 0
      p[0] = ut*leg[0] + vt*leg2[0] ;
      p[1] = ut*leg[1] + vt*leg2[1] ;
      p[2] = ut*leg[2] + vt*leg2[2] ;
      ut = DOT(p, leg) ; vt = DOT(p, leg2) ;
#endif

      CROSS(d, leg, n) ;
      dot = DOT(d, p) ; ldot = DOT(d, leg2) ;

      /* on different side of leg from 3rd vertex */
      if (!FZERO(ldot) && !FZERO(dot) && dot*ldot < 0)
        break ;
    }
    if (l >= VERTICES_PER_FACE)
      flist[nfaces++] = bin->fno ;
  }
  if (!nfaces)
  {
    float dist, min_dist ;

    fprintf(stderr, "no faces found on sphere!\n") ;
    bin = bucket->bins ;
    min_dist = 1000000.0f ; min_fno = 0 ;
    for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      f = &mris->faces[bin->fno] ;
      v1 = &mris->vertices[f->v[0]] ;
      v2 = &mris->vertices[f->v[1]] ;
      v3 = &mris->vertices[f->v[2]] ;
#define VDIST(v1,v2) (sqrt(SQR(v1->x-v2->x)+SQR(v1->y-v2->y)+SQR(v1->z-v2->z)))
      dist = VDIST(v1,v) + VDIST(v2, v) + VDIST(v3,v) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = bin->fno ;
      }
    }
                printf("min_dist = %f (min_fno=%d)\n",min_dist,min_fno);
                
  }
  else
  {
    float min_stretch, curv ;

    if (nfaces > 1)
    {
      DiagBreak() ;
      curv = 1.0f ;
    }
    else
      curv = 0.0f ;
    for ( i = 0 ; i < nfaces ; i++)/* check each face */
    {
      f = &mris->faces[flist[i]] ;
      for (l = 0 ; l < VERTICES_PER_FACE ; l++)
        mris->vertices[f->v[l]].curv = curv ;
    }

    min_stretch = 100000.0f ;
    for (i = 0 ; i < nfaces ; i++)
    {
      stretch = mrisComputeFaceStretch(mris, flist[i]) ;
      if (stretch < min_stretch)
      {
        min_stretch = stretch ;
        min_fno = flist[i] ;
      }
    }
  }
  if (nfaces <= 1)
    return(min_fno) ;
  else
    return(-1) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale a surface anisotropically.
------------------------------------------------------*/
static float
mrisComputeFaceStretch(MRI_SURFACE *mris, int fno)
{
  int    fv, n0, n1 ;
  float  stretch, max_stretch, dist, inflated_dist ;
  FACE   *f ;
  VERTEX *v0, *v1 ;

  f = &mris->faces[fno] ;
  max_stretch = -1.0f ;
  for (fv = 0 ; fv < VERTICES_PER_FACE ; fv++)
  {
    n0 = f->v[fv] ;
    n1 = fv < VERTICES_PER_FACE - 1 ? f->v[fv+1] : f->v[0] ;
    v0 = &mris->vertices[n0] ; v1 = &mris->vertices[n1] ;
    inflated_dist = 
      SQR(v0->cx-v1->cx) + SQR(v0->cy-v1->cy) + SQR(v0->cz-v1->cz);
    dist = 
      SQR(v0->origx-v1->origx) + 
      SQR(v0->origy-v1->origy) + SQR(v0->origz-v1->origz);
    if (!FZERO(dist))
    {
      stretch = inflated_dist / dist ;
      if (stretch > max_stretch)
        max_stretch = stretch ;
    }
  }
  return(max_stretch) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale a surface anisotropically.
------------------------------------------------------*/
static int
mrisCalculateCanonicalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                  float *px, float *py, float *pz)
{
  float  x, y, z ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *f ;

  f = &mris->faces[fno] ;
  v0 = &mris->vertices[f->v[0]] ; v1 = &mris->vertices[f->v[1]] ; 
  v2 = &mris->vertices[f->v[2]] ;

  /* first bisect v1->v2 line */

  x = (v1->cx + v2->cx) / 2.0f ;
  y = (v1->cy + v2->cy) / 2.0f ;
  z = (v1->cz + v2->cz) / 2.0f ;
  
  /* now bisect v0->bisector line */
  *px = (v0->cx + x) / 2.0f ;
  *py = (v0->cy + y) / 2.0f ;
  *pz = (v0->cz + z) / 2.0f ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale a surface anisotropically.
------------------------------------------------------*/
static int
mrisCalculateOriginalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                  float *px, float *py, float *pz)
{
  float  x, y, z ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *f ;

  f = &mris->faces[fno] ;
  v0 = &mris->vertices[f->v[0]] ; v1 = &mris->vertices[f->v[1]] ; 
  v2 = &mris->vertices[f->v[2]] ;

  /* first bisect v1->v2 line */

  x = (v1->origx + v2->origx) / 2.0f ;
  y = (v1->origy + v2->origy) / 2.0f ;
  z = (v1->origz + v2->origz) / 2.0f ;
  
  /* now bisect v0->bisector line */
  *px = (v0->origx + x) / 2.0f ;
  *py = (v0->origy + y) / 2.0f ;
  *pz = (v0->origz + z) / 2.0f ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale a surface anisotropically.
------------------------------------------------------*/
static int
mrisCalculateFaceCentroid(MRI_SURFACE *mris, int fno, float *px, float *py, 
                          float *pz)
{
  float  x, y, z ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *f ;

  f = &mris->faces[fno] ;
  v0 = &mris->vertices[f->v[0]] ; v1 = &mris->vertices[f->v[1]] ; 
  v2 = &mris->vertices[f->v[2]] ;

  /* first bisect v1->v2 line */

  x = (v1->x + v2->x) / 2.0f ;
  y = (v1->y + v2->y) / 2.0f ;
  z = (v1->z + v2->z) / 2.0f ;
  
  /* now bisect v0->bisector line */
  *px = (v0->x + x) / 2.0f ;
  *py = (v0->y + y) / 2.0f ;
  *pz = (v0->z + z) / 2.0f ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale a surface anisotropically.
------------------------------------------------------*/
int
MRISprintTessellationStats(MRI_SURFACE *mris, FILE *fp)
{
  double  mean, dsigma, dmin, dmax ;
  int     vno, vno2 ;

  mean = MRIScomputeVertexSpacingStats(mris, &dsigma, &dmin,&dmax,&vno,&vno2) ;
  fprintf(fp, "vertex spacing %2.2f +- %2.2f (%2.2f-->%2.2f) "
          "(max @ vno %d --> %d)\n",
          mean, dsigma, dmin, dmax, vno, vno2) ;
  mean = MRIScomputeFaceAreaStats(mris, &dsigma, &dmin, &dmax) ;
  fprintf(fp, "face area %2.2f +- %2.2f (%2.2f-->%2.2f)\n",
          mean, dsigma, dmin, dmax) ;

  if (dmax > 20)
  {
    VERTEX  *v, *vn ;
    int     n ;
    float   dist ;

    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dist = sqrt(SQR(vn->x - v->x) + SQR(vn->y - v->y) + SQR(vn->z - v->z));
      if (dist > 20)
        fprintf(stdout, "\t%d --> %d = %2.1f mm\n", vno, v->v[n], dist) ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static void
mrisDumpFace(MRI_SURFACE *mris, int fno, FILE *fp) 
{
  FACE   *f ;
  VERTEX *v ;
  int    n ;

  f = &mris->faces[fno] ;
  fprintf(fp, "face %d, area %2.1f, orig area %2.1f\n",
          fno, f->area, f->orig_area) ;
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    v = &mris->vertices[f->v[n]] ;
    fprintf(fp,"\tv %d (%d) @ (%2.1f, %2.1f, %2.1f) o (%2.1f, %2.1f, %2.1f)\n",
            n, f->v[n], v->x, v->y, v->z, v->origx, v->origy, v->origz) ;
  }
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisStoreVtotalInV3num(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->v3num = v->vtotal ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISripDefectiveFaces(MRI_SURFACE *mris)
{
  FACE   *f ;
  int    fno, flist[1000], retained_i, i, j, nfaces, nripped, n ;
  double dist, min_face_dist, max_dist, r ;
  float  dx, dy, dz, cx[1000], cy[1000], cz[1000] ;
  MHT    *mht ;

  MRISclearCurvature(mris) ;
  r = MRISaverageRadius(mris) ;
  MRISscaleBrain(mris, mris, 100.0/r) ;
  mht = MHTfillTable(mris, NULL) ;

  /*
    first remove all faces with negative areas as they will be
    the 'underside' of a defect.
    */
#if 0
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    if (mris->faces[177616].v[2] > mris->nfaces)
      DiagBreak() ;
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    if (f->area < 0)
    {
      f->ripflag = 1 ;  /* part of a defect! */
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "ripping face %d\n", fno) ;
      nripped++ ;
      for (i = 0 ; i < VERTICES_PER_FACE ; i++)
        mris->vertices[f->v[i]].curv = -2.0f ;
    }
  }
#endif
  for (nripped = fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    nfaces = mrisFindAllOverlappingFaces(mris, mht, fno, flist) ;
    if (nfaces > 1)   /* overlapping faces - rip all but one of them */
    {
      /* retain the face that is furthest from any negative face */
      max_dist = 0.0f ;
      for (i = 0 ; i < nfaces ; i++)  
        mrisCalculateFaceCentroid(mris, flist[i], cx+i, cy+i, cz+i) ;
      for (retained_i = -1, i = 0 ; i < nfaces ; i++)
      {
        if (mris->faces[flist[i]].area < 0)
          continue ;   /* don't ever retain a negative face */

        /* find the distance to the closest negative face */
        for (min_face_dist = 1000000.0, j = 0 ; j < nfaces ; j++)
        {
          if (mris->faces[flist[j]].area > 0)
            continue ;   /* only consider distances to negative faces */
          dx = cx[j] - cx[i] ; dy = cy[j] - cy[i] ; dz = cz[j] - cz[i] ; 
          dist = (dx*dx+dy*dy+dz*dz) ;
          if (dist < min_face_dist)
            min_face_dist = dist ;
        }

        /* 
           if this face is more distant than any other face to a negative
           face (so far), tentatively mark it as the one to keep.
           */
        if (min_face_dist > max_dist)
        {
          max_dist = min_face_dist ;
          retained_i = i ;
        }
      }

      if (retained_i >= 0) for (i = 0 ; i < nfaces ; i++)
      {
        VERTEX *v ;
        FACE   *f ;

        f = &mris->faces[flist[i]] ;
        if (i == retained_i)
        {
          for (n = 0 ; n < VERTICES_PER_FACE ; n++)
          {
            v = &mris->vertices[f->v[n]] ;
            if (FZERO(v->curv))
              v->curv = 1.0 ; /* good part of a defect! */
          }
        }

        if (i == retained_i || f->ripflag)
          continue ;
        f->ripflag = 1 ;  /* part of a defect! */
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "ripping face %d\n", flist[i]) ;
        nripped++ ;
        for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        {
          v = &mris->vertices[f->v[n]] ;
          if (v->curv >= 0.0)
            v->curv = f->area > 0 ? -1.0 : -2.0 ; 
        }
      }
    }
  }
  
  MRISscaleBrain(mris, mris, r/100.0) ;
  fprintf(stdout, "removing %d ripped faces from tessellation\n", nripped) ;
#if 0
  mrisRipVertices(mris) ;
#else
  {
#if 0
    int vno ;
    
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->curv = 0 ;
      for (fno = 0 ; fno < v->num ; fno++)
      {
        f = &mris->faces[v->f[fno]] ;
        if (f->ripflag == 1)
        {
          v->curv = -1 ;
          break ;
        }
        else if (f->ripflag == 2)
          v->curv = 1.0 ;
      }
    }
#endif
    MRISwriteCurvature(mris, "defects") ;
  }
#endif
  /*  MRISremoveRipped(mris) ;*/
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int edgesIntersect(MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2) ;
static int
edgesIntersect(MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2)
{
  VERTEX *v1, *v2 ;
  double  b1, b2, m1, m2, x1_start, x1_end, y1_start, y1_end,
         x2_start, x2_end, y2_start, y2_end, x, y, x1min, x1max,
         y1min, y1max, x2min, x2max, y2min, y2max, cx, cy, cz ;
  double origin[3], e0[3], e1[3], e0_0, e0_1, e0_2, e1_0, e1_1, e1_2 ;

	if (edge1->vno1 == edge2->vno1 || edge1->vno1 == edge2->vno2 ||
			edge1->vno2 == edge2->vno1 || edge1->vno2 == edge2->vno2)
		return(0) ;

  mrisComputeCanonicalEdgeBasis(mris, edge1, edge2,origin,e0,e1);
	e0_0 = e0[0] ; e0_1 = e0[1] ; e0_2 = e0[2] ; 
	e1_0 = e1[0] ; e1_1 = e1[1] ; e1_2 = e1[2] ; 
  if (edge1->vno1 == Gdiag_no || edge2->vno1 == Gdiag_no ||
      edge1->vno2 == Gdiag_no || edge2->vno2 == Gdiag_no)
    DiagBreak() ;
  v1 = &mris->vertices[edge1->vno1] ; v2 = &mris->vertices[edge1->vno2] ;
  cx = v1->cx - origin[0] ; cy = v1->cy - origin[1] ; cz = v1->cz - origin[2] ;
  x1_start = cx*e0_0 + cy*e0_1 + cz*e0_2 ;
  y1_start = cx*e1_0 + cy*e1_1 + cz*e1_2 ;
  cx = v2->cx - origin[0] ; cy = v2->cy - origin[1] ; cz = v2->cz - origin[2] ;
  x1_end = cx*e0_0 + cy*e0_1 + cz*e0_2 ;
  y1_end = cx*e1_0 + cy*e1_1 + cz*e1_2 ;
  x1min = MIN(x1_start, x1_end) ; x1max = MAX(x1_start, x1_end) ;
  y1min = MIN(y1_start, y1_end) ; y1max = MAX(y1_start, y1_end) ;
  if (!FEQUAL(x1_start, x1_end))
  {
    m1 = (y1_end - y1_start) / (x1_end - x1_start) ;
    b1 = y1_end - m1 * x1_end ;
  }
  else
    m1 = b1 = 0 ;   /* will be handled differently */

  v1 = &mris->vertices[edge2->vno1] ; v2 = &mris->vertices[edge2->vno2] ;
  cx = v1->cx-origin[0] ; cy = v1->cy-origin[1] ; cz = v1->cz-origin[2];
  x2_start = cx*e0_0 + cy*e0_1 + cz*e0_2 ;
  y2_start = cx*e1_0 + cy*e1_1 + cz*e1_2 ;
  cx = v2->cx-origin[0] ; cy = v2->cy-origin[1] ; cz = v2->cz-origin[2];
  x2_end = cx*e0_0 + cy*e0_1 + cz*e0_2 ;
  y2_end = cx*e1_0 + cy*e1_1 + cz*e1_2 ;
  x2min = MIN(x2_start, x2_end) ; x2max = MAX(x2_start, x2_end) ;
  y2min = MIN(y2_start, y2_end) ; y2max = MAX(y2_start, y2_end) ;
#if 0
#define INT_EPSILON 0.000
  x1max += INT_EPSILON ; y1max += INT_EPSILON ; 
  x2max += INT_EPSILON ; y2max += INT_EPSILON ; 
  x1min -= INT_EPSILON ; y1min -= INT_EPSILON ; 
  x2min -= INT_EPSILON ; y2min -= INT_EPSILON ; 
#endif
  
#if 1
  /* non-overlapping intervals can't intersect */
  if (x1min > x2max || x1max < x2min || y1min > y2max || y1max < y2min)
    return(0) ;
#endif
  
  /* handle special cases first */
  if (FEQUAL(x1_start, x1_end))        /* l1 is vertical */
  {
    if (FEQUAL(x2_start, x2_end))      /* both vertical */
      return (FEQUAL(x1_start, x2_start)) ;  /* same line */
    m2 = (y2_end - y2_start) / (x2_end - x2_start) ;
    b2 = y2_end - m2 * x2_end ;
    y = m2 * x1_start + b2 ;           /* y coordinate of intersection */
    if (y >= y1min && y <= y1max)
      return(1) ;
  }
  else if (FEQUAL(x2_start, x2_end))   /* l2 is vertical */
  {
    if (FEQUAL(x1_start, x1_end))      /* both vertical */
      return (FEQUAL(x1_start, x2_start)) ;  /* same line */
    y = m1 * x2_start + b1 ;  /* y coord of intersection */
    if (y >= y2min && y <= y2max)
      return(1) ;
  }
  else    /* neither line is vertical, compute intersection point */
  {

    m2 = (y2_end - y2_start) / (x2_end - x2_start) ;
    b2 = y2_end - m2 * x2_end ;
    if (FEQUAL(m1, m2))                /* parallel lines */
    {
      if (FEQUAL(b1, b2))  /* same line, see if segments overlap */
      {
        return(x2max >= x1min && x2min <= x1max) ;
      }
      return(0) ;
    }
    x = (b2 - b1) / (m1-m2) ;            /* intersection point */
    if ((x <= x1max && x >= x1min) &&    /* is it in l1 interval? */
        (x <= x2max && x >= x2min))      /* is it in l2 interval? */
      return(1) ;                        /* both true - intersection */
  }
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int 
mrisFindAllOverlappingFaces(MRI_SURFACE *mris, MHT *mht,int fno, int *flist)
{
  double   x0, x1, y0, y1, z0, z1, x, y, z ;
  int     i, n, m, total_found, all_faces[100000], nfaces ;
  MHBT    *bucket, *last_bucket ;
  MHB     *bin ;
  EDGE    edge1, edge2 ;
  FACE    *f1, *f2 ;
  VERTEX  *v ;


  f1 = &mris->faces[fno] ;
  if (fno == Gdiag_no)
    DiagBreak() ;
  x0 = y0 = z0 = 100000.0 ; x1 = y1 = z1 = -x0 ;
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    v = &mris->vertices[f1->v[n]] ;
    x = v->x ; y = v->y ; z = v->z ;
    x0 = MIN(x, x0) ; y0 = MIN(y, y0) ; z0 = MIN(z, z0) ;
    x1 = MAX(x, x1) ; y1 = MAX(y, y1) ; z1 = MAX(z, z1) ;
  }

  nfaces = total_found = 0 ; flist[total_found++] = fno ;
  last_bucket = NULL ;
  for (x = x0 ; x <= x1 ; x += 0.5)
  for (y = y0 ; y <= y1 ; y += 0.5)
  for (z = z0 ; z <= z1 ; z += 0.5)
  {
    bucket = MHTgetBucket(mht, x, y, z) ;
    if (!bucket || bucket == last_bucket)
      continue ;
    last_bucket = bucket ;
    for (bin = bucket->bins, i = 0 ; i < bucket->nused ; i++, bin++)
    {
      f2 = &mris->faces[bin->fno] ;
      if ((bin->fno == 114877 && fno == 114875) ||
          (fno == 114877 && bin->fno == 114875))
        DiagBreak() ;
      if (!f2->ripflag)  /* only add it once */
      {
        all_faces[nfaces++] = bin->fno ;
        f2->ripflag = 1 ;
      }
    }
  }
  for (i = 0 ; i < nfaces ; i++)     /* reset ripflag */
    mris->faces[all_faces[i]].ripflag = 0 ;

  for (i = 0 ; i < nfaces ; i++)
  {
    if (i == Gdiag_no)
      DiagBreak() ;
    if (all_faces[i] < 0)  /* duplicate */
      continue ;
    f2 = &mris->faces[all_faces[i]] ;
    if (all_faces[i] == Gdiag_no)
      DiagBreak() ;
    if ((all_faces[i] == 117486 && fno == 114877) ||
        (fno == 117486 && all_faces[i] == 114877))
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      edge1.vno1 = f1->v[n] ; 
      edge1.vno2 = f1->v[n == VERTICES_PER_FACE-1 ? 0 : n+1] ;
      for (m = 0 ; m < VERTICES_PER_FACE ; m++)
      {
        if (f2->ripflag)
          continue ;
        edge2.vno1 = f2->v[m] ;
        edge2.vno2 = f2->v[m == VERTICES_PER_FACE-1 ? 0 : m+1] ;
        if (edge2.vno1 == edge1.vno1 || edge2.vno1 == edge1.vno2 ||
            edge2.vno2 == edge1.vno1 || edge2.vno2 == edge1.vno2)
          continue ;
        if (edgesIntersect(mris, &edge1, &edge2))
        {
          f2->ripflag = 1 ;  /* use ripflag as a mark */
          flist[total_found++] = all_faces[i] ;
        }
      }
    }
  }

  for (i = 0 ; i < total_found ; i++)
  {
    f1 = &mris->faces[flist[i]] ;
    f1->ripflag = 0 ;
  }
  return(total_found) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          restore ripflag to 0 (NOTE: Won't undo MRISremoveRipped!!)
------------------------------------------------------*/
int
MRISunrip(MRI_SURFACE *mris)
{
  int    vno, fno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].ripflag = 0 ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
    mris->faces[fno].ripflag = 0 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MARK_AMBIGUOUS_RETAIN   1
#define MARK_RETAIN             MARK_AMBIGUOUS_RETAIN
#define MARK_AMBIGUOUS          MARK_AMBIGUOUS_RETAIN   
#define MARK_AMBIGUOUS_DISCARD  2
#define MARK_DISCARD            MARK_AMBIGUOUS_DISCARD
#define MARK_SEGMENTED          3
#define MAX_DEFECTS             25000


typedef struct
{
  int **faces ;   /* array of pointer to list of ambiguous faces */
  int *nfaces ;   /* array of ints specifying how many faces are ambigous */
} FACE_DEFECT_LIST, FDL ;

#define KEEP_VERTEX     0
#define DISCARD_VERTEX  1

typedef struct
{
  float  nx, ny, nz ;    /* average normal in the defect */
  float  cx, cy, cz ;    /* centroid of the defect */
  float  area ;          /* total (original) area of the defect */
  int    *vertices ;     /* vertices in the defect */    
  char   *status ;       /* keep or discard */
  int    nvertices ;     /* # of vertices in the defect */
  int    *border ;
  int    nborder ;
  int    *chull ;        /* vertices in the convex hull */
  int    nchull ;        /* # of vertices in the convex hull */
} DEFECT ;

typedef struct
{
  int    ndefects ;
  DEFECT defects[MAX_DEFECTS] ;
} DEFECT_LIST, DL ;

#define ET_OVERLAP_LIST_INCOMPLETE  0x0001

typedef struct
{
	int  nedges ;
	EDGE *edges ;
	int  **overlapping_edges ;  /* indices of all edges overlapping this one (i.e. [i][j]) */
	int  *noverlap ;
	unsigned char *flags ;
} EDGE_TABLE ;

typedef struct
{
	double        fitness ;
	int           nedges ;
	int           *ordering ;  /* order of edges in EDGE_TABLE for retessellation */
	DEFECT        *defect ;
	EDGE_TABLE    *etable ;
	int           rank ;
} DEFECT_PATCH, DP ;

typedef struct
{
	int  vno ;     /* vertex # in surface */
	int  vnum ;    /* original # of 1 neighbors */
	int  vtotal ;
	int  *v ;      /* list of original neighbors */
} VERTEX_STATE, VS ;

typedef struct
{
	DEFECT       *defect ;
	VERTEX_STATE *vs ;
	int          nvertices ;
	int          *vertex_trans ;  /* not allocated - pointer to preexisting table */
} DEFECT_VERTEX_STATE, DVS ;

static DEFECT_VERTEX_STATE *mrisRecordVertexState(MRI_SURFACE *mris, DEFECT *defect,
																									int *vertex_trans) ;
static int mrisRestoreVertexState(MRI_SURFACE *mris, DEFECT_VERTEX_STATE *dvs) ;

static int mrisAddAllDefectFaces(MRI_SURFACE *mris, DEFECT_LIST *dl, 
                                 int *vertex_trans) ;
static int mrisFindDefectConvexHull(MRI_SURFACE *mris, DEFECT *defect) ;
static int mrisOrientRetessellatedSurface(MRI_SURFACE *mris, DEFECT_LIST *dl,
                                          int *vtrans) ;
static int containsAnotherVertex(MRI_SURFACE *mris, int vno0, int vno1, 
                                 int vno2, double e0[3], double e1[3], 
                                 double origin[3]) ;
static int       mrisMarkDefect(MRI_SURFACE *mris, DEFECT *defect, int mark) ;
static int       mrisMarkDefectBorder(MRI_SURFACE *mris, DEFECT *defect,
                                      int mark);
static int       mrisMarkDefectConvexHull(MRI_SURFACE *mris, DEFECT *defect,
                                      int mark);
#if 0
static int mrisAddDefectFaces(MRI_SURFACE *mris, double e0[3], double e1[3], 
                              double origin[3], EDGE *et, int nedges);
static int       mrisOrientFaces(MRI_SURFACE *mris,DEFECT *defect,int *vtrans);
static int       mrisRestoreDefectPositions(MRI_SURFACE *mris, DEFECT *defect,
                                            int which) ;
#endif
FACE_DEFECT_LIST *MRISmarkAmbiguousVertices(MRI_SURFACE *mris, int mark) ;
DEFECT_LIST      *MRISsegmentDefects(MRI_SURFACE *mris, int mark_ambiguous,
                                     int mark_segmented) ;
static int       mrisSegmentDefect(MRI_SURFACE *mris,int vno,DEFECT *defect,
                                   int mark_ambiguous, int mark_segmented);
static int       mrisMarkRetainedPartOfDefect(MRI_SURFACE *mris, 
                                              DEFECT *defect, 
                                              FACE_DEFECT_LIST *fdl,
                                              float area_threshold,
                                              int mark_retain, 
                                              int mark_discard,
                                              MHT *mht) ;

static int       mrisTessellateDefect(MRI_SURFACE *mris, 
                                      MRI_SURFACE *mris_corrected, 
                                      DEFECT *defect, int *vertex_trans,
                                      MRI *mri,
																			HISTOGRAM *h_k1,
																			HISTOGRAM *h_k2,
																			HISTOGRAM *h_white,
																			HISTOGRAM *h_gray,
																			HISTOGRAM *h_border,
																			HISTOGRAM *h_grad,
																			MRI *mri_gray_white,
																			HISTOGRAM *h_dot, 
																			TOPOLOGY_PARMS *parms
																			) ;
static int      mrisDefectRemoveDegenerateVertices(MRI_SURFACE *mris, 
                                                   float min_sphere_dist,
                                                   DEFECT *defect) ;
static int       mrisDefectRemoveProximalVertices(MRI_SURFACE *mris, 
                                                  float min_orig_dist,
                                                  DEFECT *defect) ;
static int       mrisDefectRemoveNegativeVertices(MRI_SURFACE *mris, 
                                                  DEFECT *defect) ;
static int       mrisComputeDefectTangentPlanes(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans) ;
static int       intersectDefectEdges(MRI_SURFACE *mris, DEFECT *defect, 
                               EDGE *e, int *vertex_trans);
static int       intersectDefectConvexHullEdges(MRI_SURFACE *mris, DEFECT *defect, 
																								EDGE *e, int *vertex_trans);
#if 0
static int       colinearDefectEdges(MRI_SURFACE *mris, DEFECT *defect, 
																		 EDGE *e, int *vertex_trans);
#endif
#if 0
static int       mrisCheckDefectEdges(MRI_SURFACE *mris, DEFECT *defect,
                                      int vno, int *vertex_trans) ;
#endif
static int       vertexNeighbor(MRI_SURFACE *mris, int vno1, int vno2) ;
static int       mrisRetessellateDefect(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, 
																				DEFECT *defect, int *vertex_trans, EDGE *et, 
																				int nedges, int *ordering, EDGE_TABLE *etable) ;
static double    mrisDefectPatchFitness(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, 
																				MRI *mri, DEFECT_PATCH *dp, int *vertex_trans,
																				DEFECT_VERTEX_STATE *dvs,
																				HISTOGRAM *h_k1,
																				HISTOGRAM *h_k2,
																				HISTOGRAM *h_white,
																				HISTOGRAM *h_gray,
																				HISTOGRAM *h_border,
																				HISTOGRAM *h_grad,
																				MRI *mri_gray_white,
																				HISTOGRAM *h_dot, TOPOLOGY_PARMS *parms
																				) ;
static double mrisComputeDefectLogLikelihood(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, 
																						 int *vertex_trans, DEFECT_PATCH *dp,
																						 HISTOGRAM *h_k1,
																						 HISTOGRAM *h_k2,
																						 HISTOGRAM *h_white,
																						 HISTOGRAM *h_gray,
																						 HISTOGRAM *h_border,
																						 HISTOGRAM *h_grad,
																						 MRI *mri_gray_white,
																						 HISTOGRAM *h_dot, TOPOLOGY_PARMS *parms) ;
#if 0
static double mrisComputeDefectEnergy(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, 
																			int *vertex_trans) ;
static double mrisComputeDefectQuadraticCurvatureEnergy(MRI_SURFACE *mris, DEFECT *defect, 
																												int *vertex_trans) ;
static double mrisComputeDefectCurvatureEnergy(MRI_SURFACE *mris, DEFECT *defect, 
																							 int *vertex_trans) ;
static double mrisComputeDefectMRIEnergy(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, 
																				 int *vertex_trans) ;
#endif
static double mrisComputeDefectMRILogLikelihood(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, 
																								int *vertex_trans, HISTOGRAM *h_white,
																								HISTOGRAM *h_gray, HISTOGRAM *h_grad,
																								MRI *mri_gray_white) ;
static double mrisComputeDefectMRILogUnlikelihood(MRI_SURFACE *mris, MRI *mri, DEFECT_PATCH *dp, 
																									int *vertex_trans, HISTOGRAM *h_border) ;
static double mrisComputeDefectCurvatureLogLikelihood(MRI_SURFACE *mris, DEFECT *defect, 
																											int *vertex_trans, HISTOGRAM *h_k1, HISTOGRAM *h_k2) ;
static int    mrisComputeOptimalRetessellation(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected,
																							 MRI *mri, DEFECT *defect, int *vertex_trans,
																							 EDGE *et, int nedges, HISTOGRAM *h_k1, HISTOGRAM *h_k2,
																							 HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_border,
																							 HISTOGRAM *h_grad, 
																							 MRI *mri_gray_white,
																							 HISTOGRAM *h_dot, TOPOLOGY_PARMS *parms) ;



#define AREA_THRESHOLD     35.0f


static int mrisFreeDefectVertexState(DEFECT_VERTEX_STATE *dvs) ;

static int
mrisFreeDefectVertexState(DEFECT_VERTEX_STATE *dvs)
{
	int  i ;

	for (i = 0 ; i < dvs->nvertices ; i++)
		free(dvs->vs[i].v) ;

	free(dvs->vs) ;
	free(dvs) ;
	return(NO_ERROR) ;
}
static DEFECT_VERTEX_STATE *
mrisRecordVertexState(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans)
{
	DEFECT_VERTEX_STATE *dvs ;
	int                 i, n, vno ;
	VERTEX              *v ;
	VERTEX_STATE        *vs ;

	dvs = calloc(1, sizeof(DVS)) ;
	if (!dvs)
		ErrorExit(ERROR_NOMEMORY, "mrisRecordVertexState: could not allocate dvs") ;

	dvs->defect = defect ;
	dvs->vertex_trans = vertex_trans ;
	dvs->nvertices = defect->nvertices+defect->nborder ;
	dvs->vs = (VS *)calloc(dvs->nvertices, sizeof(VS)) ;
	if (!dvs->vs)
		ErrorExit(ERROR_NOMEMORY, "mrisRecordVertexState: could not allocate %d dvs->vs",
							dvs->nvertices) ;
	for (n = 0 ; n < defect->nvertices ; n++)
		dvs->vs[n].vno = vertex_trans[defect->vertices[n]] ;
	for (n = 0 ; n < defect->nborder ; n++)
		dvs->vs[defect->nvertices+n].vno = vertex_trans[defect->border[n]] ;

	for (i = 0 ; i < dvs->nvertices ; i++)
	{
		vs = &dvs->vs[i] ;
		vno = vs->vno ;
		if (vno < 0)
			continue ;
		v = &mris->vertices[vno] ;
		vs->vtotal = v->vtotal ; vs->vnum = v->vnum ;
		if (!v->vtotal)
			continue ;
		vs->v = (int *)calloc(vs->vtotal, sizeof(int)) ;
		if (!vs->v)
			ErrorExit(ERROR_NOMEMORY, "mrisRecordVertexState: could not allocate %dth array of %d elts", i, vs->vtotal) ;
		for (n = 0 ; n < v->vtotal ; n++)
			vs->v[n] = v->v[n] ;
	}

	return(dvs) ;
}

static int
mrisRestoreVertexState(MRI_SURFACE *mris, DEFECT_VERTEX_STATE *dvs)
{
	int                 i, n, vno ;
	VERTEX              *v ;
	VERTEX_STATE        *vs ;

	for (i = 0 ; i < dvs->nvertices ; i++)
	{
		vs = &dvs->vs[i] ;
		vno = vs->vno ;
		if (vno < 0)
			continue ;
		v = &mris->vertices[vno] ;
		free(v->v) ;
		v->v = NULL ;
		v->vtotal = vs->vtotal ; v->vnum = vs->vnum ;
		if (!v->vtotal)
			continue ;
		v->v = (int *)calloc(vs->vtotal, sizeof(int)) ;
		if (!v->v)
			ErrorExit(ERROR_NOMEMORY, 
								"mrisRestoreVertexState: could not allocate %dth array of %d elts", 
								i, vs->vtotal) ;
		for (n = 0 ; n < v->vtotal ; n++)
			v->v[n] = vs->v[n] ;
	}

	return(NO_ERROR) ;
}

static int mrisComputeJointGrayWhiteBorderDistributions(MRI_SURFACE *mris, MRI *mri, 
																												MRI *mri_gray_white, MRI *mri_wm) ;
static int mrisComputeGrayWhiteBorderDistributions(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, 
																									 HISTOGRAM *h_white, HISTOGRAM *h_gray, 
																									 HISTOGRAM *h_border, HISTOGRAM *h_grad);
static int mrisComputePrincipalCurvatureDistributions(MRI_SURFACE *mris, HISTOGRAM *h_k1, HISTOGRAM *h_k2) ;
static int mrisComputeNormalDotDistribution(MRI_SURFACE *mris, HISTOGRAM *h_dot) ;
static int
mrisComputeNormalDotDistribution(MRI_SURFACE *mris, HISTOGRAM *h_dot)
{
	int    vno, bin, n, num ;
	VERTEX *v, *vn ;
	float  bin_size, min_dot, max_dot, bin_val, dot, dx, dy, dz, nx, ny, nz, x, y, z ;
	HISTOGRAM *h_raw;

	MRISsaveVertexPositions(mris, TMP_VERTICES) ; MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
	MRIScomputeMetricProperties(mris) ;
	min_dot = 100000 ; 
	max_dot = -100000 ; 

	/* first compute min and max */
	for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		if (v->ripflag)
			continue ;

    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    x = v->x ;    y = v->y ;   z = v->z ;
		for (n = 0 ; n < v->vnum ; n++)
		{
			vn = &mris->vertices[v->v[n]] ;
			dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ; 
			dot = dx*nx + dy*ny * dz*nz ;
			if (dot < min_dot)
				min_dot = dot ;
			if (dot > max_dot)
				max_dot = dot ;
		}
	}

	/* add one bin at either end for almost zero probability events */
	bin_size = (max_dot - min_dot) / (h_dot->nbins-2) ;
	h_dot->bin_size = bin_size ;
	for (bin_val = min_dot-bin_size, bin = 0 ; bin < h_dot->nbins ; bin++, bin_val += bin_size)
		h_dot->bins[bin] = bin_val ;

	min_dot = h_dot->bins[0] ;

	/* now fill in distribution */
	for (num = vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		if (v->ripflag)
			continue ;

    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    x = v->x ;    y = v->y ;   z = v->z ;
		for (n = 0 ; n < v->vnum ; n++)
		{
			num++ ;
			vn = &mris->vertices[v->v[n]] ;
			dx = vn->x - x ; dy = vn->y - y ; dz = vn->z - z ; 
			dot = dx*nx + dy*ny * dz*nz ;
			bin = (int)((dot - min_dot) / bin_size) ;
			if (bin == 0)
				DiagBreak() ;
			h_dot->counts[bin]++ ;
		}
	}

	for (bin = 0 ; bin < h_dot->nbins ; bin++)
	{
		if (h_dot->counts[bin] == 0)
			h_dot->counts[bin] = 0.01 ;
		h_dot->counts[bin] /= (float)num ;
	}
	h_raw = HISTOcopy(h_dot, NULL) ;
	HISTOsmooth(h_raw, h_dot, 2.0) ;
	
	/*	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)*/
	{
		HISTOplot(h_dot, "ndot.plt") ;
		HISTOplot(h_raw, "rdot.plt") ;
	}
	HISTOfree(&h_raw) ;
	MRISrestoreVertexPositions(mris, TMP_VERTICES) ; MRIScomputeMetricProperties(mris) ;
	return(NO_ERROR) ;
}
static int
mrisComputePrincipalCurvatureDistributions(MRI_SURFACE *mris, HISTOGRAM *h_k1, HISTOGRAM *h_k2)
{
	int    vno, bin ;
	VERTEX *v ;
	float  k1_bin_size, k2_bin_size, min_k1, max_k1, min_k2, max_k2, bin_val ;
#if 0
	HISTOGRAM *h_k1_raw, *h_k2_raw ;
#endif

	MRIScomputeSecondFundamentalForm(mris) ;
	
	min_k1 = min_k2 = 100000 ; 
	max_k1 = max_k2 = -100000 ; 
	for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		if (v->ripflag)
			continue ;
		if (v->k1 < min_k1)
			min_k1 = v->k1 ;
		if (v->k2 < min_k2)
			min_k2 = v->k2 ;
		if (v->k1 > max_k1)
			max_k1 = v->k1 ;
		if (v->k2 > max_k2)
			max_k2 = v->k2 ;
	}

	
	k1_bin_size = (max_k1 - min_k1) / h_k1->nbins ;
	k2_bin_size = (max_k2 - min_k2) / h_k2->nbins ;

	h_k1->bin_size = k1_bin_size ; h_k2->bin_size = k2_bin_size ; 
	for (bin_val = min_k1, bin = 0 ; bin < h_k1->nbins ; bin++, bin_val += k1_bin_size)
		h_k1->bins[bin] = bin_val ;
	for (bin_val = min_k2, bin = 0 ; bin < h_k2->nbins ; bin++, bin_val += k2_bin_size)
		h_k2->bins[bin] = bin_val ;

	for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		if (v->ripflag)
			continue ;
		bin = (int)((v->k1 - min_k1) / k1_bin_size) ;
		h_k1->counts[bin]++ ;
		bin = (int)((v->k2 - min_k2) / k2_bin_size) ;
		h_k2->counts[bin]++ ;
	}
	for (bin = 0 ; bin < h_k1->nbins ; bin++)
	{
		if (h_k1->counts[bin] == 0)
			h_k1->counts[bin] = 0.01 ;
		h_k1->counts[bin] /= (float)mris->nvertices ; ;
	}
	for (bin = 0 ; bin < h_k2->nbins ; bin++)
	{
		if (h_k2->counts[bin] == 0)
			h_k2->counts[bin] = 0.01 ;
		h_k2->counts[bin] /= (float)mris->nvertices ; ;
	}
#if 0
	h_k1_raw = HISTOcopy(h_k1, NULL) ;
	h_k2_raw = HISTOcopy(h_k2, NULL) ;
	HISTOsmooth(h_k1_raw, h_k1, 2.0) ;
	HISTOsmooth(h_k2_raw, h_k2, 2.0) ;
	HISTOplot(h_k1_raw, "k1r.plt") ;
	HISTOplot(h_k2_raw, "k2r.plt") ;
#endif
	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
	{
		HISTOplot(h_k1, "k1.plt") ;
		HISTOplot(h_k2, "k2.plt") ;
	}
	return(NO_ERROR) ;
}


static long ncross = 0 ;
static long nmut = 0 ;

static int mrisMarkAllDefects(MRI_SURFACE *mris, DEFECT_LIST *dl, int flag) ;
static int mrisRipAllDefects(MRI_SURFACE *mris, DEFECT_LIST *dl, int ripflag) ;
static int mrisRipDefect(MRI_SURFACE *mris, DEFECT *defect, int ripflag) ;

MRI_SURFACE *
MRIScorrectTopology(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, MRI *mri, MRI *mri_wm,
                    int nsmooth, TOPOLOGY_PARMS *parms)
{
  FACE_DEFECT_LIST   *fdl ;
  DEFECT_LIST        *dl ;
  DEFECT             *defect ;
  int                fno, i, n, vno, kept_vertices, *face_trans,*vertex_trans;
  MHT                *mht ;
  VERTEX             *v, *vdst ;
  FACE               *f, *fdst ;
	HISTOGRAM          *h_k1, *h_k2, *h_gray, *h_white, *h_dot, *h_border, *h_grad ;
	MRI                *mri_gray_white ;
#if 0
  float              max_len ;
#endif

  /*  mrisCheckSurface(mris) ;*/
  fdl = MRISmarkAmbiguousVertices(mris, MARK_AMBIGUOUS) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "segmenting defects...\n") ;
  dl = MRISsegmentDefects(mris, MARK_AMBIGUOUS, MARK_SEGMENTED) ;
  MRISsetVals(mris, 0.0f) ;
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    defect = &dl->defects[i] ;
    for (n = 0 ; n < defect->nvertices ; n++)
    {
      mris->vertices[defect->vertices[n]].val = defect->area ;
    }
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRISwriteValues(mris, "defect_area") ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "%d defects found, arbitrating ambiguous regions...\n",
            dl->ndefects) ;

#if 0
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISreadVertexPositions(mris, "inflated") ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
#endif
#if 0
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
#endif
  MRIScomputeMetricProperties(mris) ;
  mrisScaleMaxDimension(mris, FIELD_OF_VIEW*.9f) ;
  MRISclearMarks(mris) ; MRISclearCurvature(mris) ;
  mht = MHTfillTable(mris, NULL) ;
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    defect = &dl->defects[i] ;
    if (i == 47)
      DiagBreak() ;
#if 0
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "\rdefect %d (area %2.3f)      ", 
              i, defect->area) ;
#endif
    mrisMarkRetainedPartOfDefect(mris, defect, fdl, AREA_THRESHOLD,
                                 MARK_RETAIN, MARK_DISCARD, mht) ;
    mrisFindDefectConvexHull(mris, defect) ;
  }
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    int    vno2, n2 ;
    VERTEX *vn ;

    defect = &dl->defects[i] ;
    for (n = 0 ; n < defect->nvertices+defect->nborder ; n++)
    {
      if (n < defect->nvertices)
      {
        if (defect->status[n] == DISCARD_VERTEX)
          continue ;
        vno = defect->vertices[n] ;
      }
      else
        vno = defect->border[n-defect->nvertices] ;
      v = &mris->vertices[vno] ;
      for (n2 = n+1 ; n2 < defect->nvertices+defect->nborder ; n2++)
      {
        if (n2 < defect->nvertices)
        {
          if (defect->status[n2] == DISCARD_VERTEX)
            continue ;
          vno2 = defect->vertices[n2] ;
        }
        else
          vno2 = defect->border[n2-defect->nvertices] ;
        if (vno == vno2)
          continue ;
        vn = &mris->vertices[vno2] ;
        if (FEQUAL(vn->x,v->x) && FEQUAL(vn->y,v->y) && FEQUAL(vn->z,v->z))
        {
          fprintf(stdout, "defect %d, vertices %d and %d coincident!\n",
                  i, vno, vno2) ;
        }
      }
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "\n") ;
  MHTfree(&mht) ;

  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
	/* v->val = border, v->val2 = white, v->val2bak = gray */
	mrisRipAllDefects(mris, dl, 1) ;
  mrisFindGrayWhiteBorderMean(mris, mri) ;
	mrisRipAllDefects(mris, dl, 0) ;

#if 1
	MRISsaveVertexPositions(mris, TMP_VERTICES) ;
	MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
	MRISaverageVertexPositions(mris, 2) ;
	MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
	MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
#endif
	MRISsetNeighborhoodSize(mris, 2) ;
	h_k1 = HISTOalloc(100) ; h_k2 = HISTOalloc(100) ;
	h_dot = HISTOalloc(100) ;
	mrisComputePrincipalCurvatureDistributions(mris, h_k1, h_k2) ;
	mrisComputeNormalDotDistribution(mris, h_dot) ;


  /* now knit each defect together by retessellating the surface,
     using the spherical space for topology (i.e. edge intersection),
     and the original space for geometry (i.e. edge length).
  */
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;  /* bring inflated back */


#if 0
  MRIScopyValuesToCurvature(mris) ; MRISwriteCurvature(mris, "defect_types");
  MRIScopyImagValuesToCurvature(mris) ; 
#endif
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    defect = &dl->defects[i] ;
    for (n = 0 ; n < defect->nvertices ; n++)
    {
      mris->vertices[defect->vertices[n]].curv = 
        defect->status[n] == DISCARD_VERTEX ? -1 : 1 ;
    }
  }
  if (((Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)) || (Gdiag & DIAG_SAVE_DIAGS))
    MRISwriteCurvature(mris, "defect_status");

  /*  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;*/
  MRIScomputeMetricProperties(mris) ;
  MRISclearCurvature(mris) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    for (i = 0 ; i < dl->ndefects ; i++)
    {
      int   total_defective_vertices ;
      float total_defective_area ;
      FILE *fp ;
      char fname[STRLEN] ;
      
      sprintf(fname, "%s.%s.defect%d.log", 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", 
              mris->subject_name,i) ;
      fp = fopen(fname, "wb") ;
      fprintf(fp, "%d %2.3f\n", dl->defects[i].nvertices,dl->defects[i].area);
      for (n = 0 ; n < dl->defects[i].nvertices ; n++)
      {
        fprintf(fp, "%d\n", dl->defects[i].vertices[n]) ;
        mris->vertices[dl->defects[i].vertices[n]].curv = (float)i+1 ;
        if (dl->defects[i].vertices[n] == Gdiag_no)
          DiagBreak() ;
      }
      fprintf(fp, "\nborder (%d)\n", dl->defects[i].nborder) ;
      for (n = 0 ; n < dl->defects[i].nborder ; n++)
      {
        fprintf(fp, "%d\n", dl->defects[i].border[n]) ;
        if (dl->defects[i].border[n] == Gdiag_no)
          DiagBreak() ;
      }
#if 0
      if (i != Gdiag_no)
        continue ;
#endif
      fprintf(fp, "\nconvex hull (%d)\n", dl->defects[i].nchull) ;
      for (n = 0 ; n < dl->defects[i].nchull ; n++)
      {
        fprintf(fp, "%d\n", dl->defects[i].chull[n]) ;
        if (dl->defects[i].chull[n] == Gdiag_no)
          DiagBreak() ;
      }
      fclose(fp) ;
      
      sprintf(fname, "%s.%s.defects.log", 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh",
              mris->subject_name) ;
      fp = fopen(fname, "wb") ;
      for (total_defective_area = 0.0f, total_defective_vertices = i = 0 ; 
           i < dl->ndefects ; i++)
      {
        total_defective_vertices += dl->defects[i].nvertices ;
        total_defective_area += dl->defects[i].area ;
      }
      fprintf(fp, "%d %2.3f\n", total_defective_vertices,total_defective_area);
      for (i = 0 ; i < dl->ndefects ; i++)
      {
        for (n = 0 ; n < dl->defects[i].nvertices ; n++)
          fprintf(fp, "%d\n", dl->defects[i].vertices[n]) ;
      }
      fclose(fp) ;
    }
  }
  if (Gdiag & DIAG_WRITE)
  {
    MRISclearCurvature(mris) ;
    for (i = 0 ; i < dl->ndefects ; i++)
    {
      defect = &dl->defects[i] ;
      for (n = 0 ; n < defect->nvertices ; n++)
      {
        mris->vertices[defect->vertices[n]].curv = i+1 ;
      }
    }
    if (DIAG_VERBOSE_ON || (Gdiag & DIAG_SAVE_DIAGS))
      MRISwriteCurvature(mris, "defect_labels") ;
  }

  /* now start building the target surface */
  MRISclearMarks(mris) ;
  kept_vertices = mris->nvertices ;
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    mrisMarkDefect(mris, &dl->defects[i], 1) ;
    for (n = 0 ; n < dl->defects[i].nvertices ; n++)
      if (dl->defects[i].status[n] == DISCARD_VERTEX)
        kept_vertices-- ;
  }

  face_trans  =(int *)calloc(mris->nfaces, sizeof(int)) ;
  vertex_trans = (int *)calloc(mris->nvertices, sizeof(int)) ;
  memset(vertex_trans, -1, mris->nvertices*sizeof(int)) ;
  memset(face_trans, -1, mris->nfaces*sizeof(int)) ;
  mris_corrected = MRISoverAlloc(mris->nvertices+10, 2*mris->nfaces,
                                 kept_vertices, 2*mris->nfaces) ;
                                 
  mris_corrected->type = MRIS_TRIANGULAR_SURFACE ;
#if 0
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
#endif
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;  /* inflated */
  MRIScomputeMetricProperties(mris) ;  
  for (mris_corrected->nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->marked)
      continue ;
    vdst = &mris_corrected->vertices[mris_corrected->nvertices] ;
    if (mris_corrected->nvertices == Gdiag_no)
      DiagBreak() ;
    
    vdst->nsize = v->nsize ;
    vdst->x = v->x ; vdst->y = v->y ; vdst->z = v->z ; 
    vdst->origx = v->origx ; vdst->origy = v->origy ; vdst->origz = v->origz ; 
    vdst->tx = v->tx ; vdst->ty = v->ty ; vdst->tz = v->tz ; 
    vdst->nx = v->nx ; vdst->ny = v->ny ; vdst->nz = v->nz ; 
    vdst->cx = v->cx ; vdst->cy = v->cy ; vdst->cz = v->cz ; 
    vdst->num = v->num ; 
		vdst->val = v->val ;  vdst->val2 = v->val2 ; 
		vdst->valbak = v->valbak ; vdst->val2bak = v->val2bak ; 
		vdst->imag_val = v->imag_val ; 
		vdst->curv = v->curv ;  vdst->curvbak = v->curvbak ; vdst->stat = v->stat ;
		vdst->mean = v->mean ; vdst->mean_imag = v->mean_imag ; vdst->std_error = v->std_error;
		vdst->H = v->H ; vdst->K = v->K ; vdst->k1 = v->k1 ; vdst->k2 = v->k2 ; 
    vertex_trans[vno] = mris_corrected->nvertices++ ;
  }
  /* now add all the retained vertices in the defects */
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    defect = &dl->defects[i] ;
    for (n = 0 ; n < defect->nvertices ; n++)
    {
      if (defect->vertices[n] == Gdiag_no)
        DiagBreak() ;
      if (defect->status[n] == KEEP_VERTEX)
      {
        vno = defect->vertices[n] ; v = &mris->vertices[vno] ;
        if (vno == Gdiag_no)
          DiagBreak() ;
        vdst = &mris_corrected->vertices[mris_corrected->nvertices] ;
        if (mris_corrected->nvertices == Gdiag_no)
          DiagBreak() ;

        vdst->nsize = v->nsize ;
        vdst->x = v->x ; vdst->y = v->y ; vdst->z = v->z ; 
        vdst->origx = v->origx; vdst->origy = v->origy; vdst->origz = v->origz;
        vdst->tx = v->tx ; vdst->ty = v->ty ; vdst->tz = v->tz ; 
        vdst->cx = v->cx ; vdst->cy = v->cy ; vdst->cz = v->cz ; 
        vdst->nx = v->nx ; vdst->ny = v->ny ; vdst->nz = v->nz ; 
				vdst->val = v->val ;  vdst->val2 = v->val2 ; 
				vdst->valbak = v->valbak ; vdst->val2bak = v->val2bak ; 
				vdst->imag_val = v->imag_val ; 
				vdst->curv = v->curv ;  vdst->curvbak = v->curvbak ; vdst->stat = v->stat ;
				vdst->mean = v->mean ; vdst->mean_imag = v->mean_imag ; vdst->std_error = v->std_error;
				vdst->H = v->H ; vdst->K = v->K ; vdst->k1 = v->k1 ; vdst->k2 = v->k2 ; 
        vdst->num = vdst->vnum = 0 ; 
        vertex_trans[vno] = mris_corrected->nvertices++ ;
      }
    }
  }

  for (mris_corrected->nfaces = fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (triangleMarked(mris, fno))
      continue ;
    fdst = &mris_corrected->faces[mris_corrected->nfaces] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      fdst->v[n] = vertex_trans[f->v[n]] ;
    face_trans[fno] = mris_corrected->nfaces++ ;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    FILE  *fp ;
    char  fname[STRLEN] ;
    sprintf(fname, "%s.%s.vtrans.log", 
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh",
            mris->subject_name) ;
    fp = fopen(fname, "wb") ;
    if (!fp)
      DiagBreak() ;

    for (vno = 0 ; vno < mris->nvertices ; vno++)
      fprintf(fp, "%6d --> %6d\n", vno, vertex_trans[vno]) ;
    fclose(fp);
    sprintf(fname, "%s.%s.ftrans.log", 
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh",
            mris->subject_name) ;
    fp = fopen(fname, "wb") ;

    for (vno = 0 ; vno < mris->nfaces ; vno++)
      fprintf(fp, "%6d --> %6d\n", vno, face_trans[vno]) ;
    fclose(fp);
  }

  /* now allocate face and neighbor stuff in mris_corrected */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked)
      continue ;
    if (vertex_trans[vno]<0 || vertex_trans[vno] >= mris_corrected->nvertices)
      continue ;
    vdst = &mris_corrected->vertices[vertex_trans[vno]] ;
    
    /* count # of good trangles attached to this vertex */
    for (vdst->num = n = 0 ; n < v->num ; n++)
      if (triangleMarked(mris, v->f[n]) == 0)
        vdst->num++ ;
    vdst->f = (int *)calloc(vdst->num, sizeof(int)) ;
    vdst->n = (uchar *)calloc(vdst->num, sizeof(uchar)) ;
    for (i = n = 0 ; n < v->num ; n++)
    {
      if (triangleMarked(mris, v->f[n]))
        continue ;
      vdst->n[i] = v->n[n] ; vdst->f[i] = face_trans[v->f[n]] ; i++ ;
    }
    /* count # of valid neighbors */
    for (n = vdst->vnum = 0 ; n < v->vnum ; n++)
      if (mris->vertices[v->v[n]].marked == 0)
        vdst->vnum++ ;
    vdst->vtotal = vdst->vnum ;
    vdst->v = (int *)calloc(vdst->vnum, sizeof(int)) ;
    for (i = n = 0 ; n < v->vnum ; n++)
      if (mris->vertices[v->v[n]].marked == 0)
        vdst->v[i++] = vertex_trans[v->v[n]] ;
  }

  MRISclearMarks(mris) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MRISsmoothSurfaceNormals(mris, 10) ;
#if 0
	{
		int  min_vertices, min_i ;

		min_vertices = dl->defects[0].nvertices ; min_i = 0 ;
		for (i = 0 ; i < dl->ndefects ; i++)
		{
			defect = &dl->defects[i] ;
			if (defect->nvertices < min_vertices)
			{
				min_vertices = defect->nvertices ;
				min_i = i ;
			}
#if 0
			fprintf(stdout, 
							"\rretessellating defect %d with %d vertices (chull=%d).    ", 
							i, defect->nvertices+defect->nborder, defect->nchull) ;
#endif
		}
		if (Gdiag & 0x1000000)
			mrisTessellateDefect(mris, mris_corrected, &dl->defects[min_i], vertex_trans, mri,) ;
	}
#endif
	h_gray = HISTOalloc(256) ; h_white = HISTOalloc(256) ; h_border = HISTOalloc(256) ;
	h_grad = HISTOalloc(256) ; 
	mri_gray_white = MRIalloc(256, 256, 1, MRI_FLOAT) ;
	mrisMarkAllDefects(mris, dl, 1) ;
	mrisComputeJointGrayWhiteBorderDistributions(mris, mri, mri_gray_white, mri_wm) ;
	mrisMarkAllDefects(mris, dl, 0) ;
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    defect = &dl->defects[i] ;
		if (i == Gdiag_no)
			DiagBreak() ;
#if 0
    fprintf(stdout, 
            "\rretessellating defect %d with %d vertices (chull=%d).    ", 
            i, defect->nvertices+defect->nborder, defect->nchull) ;
#endif
		mrisMarkAllDefects(mris, dl, 1) ;
		mrisComputeGrayWhiteBorderDistributions(mris, mri, defect, h_white,h_gray, h_border, h_grad) ;
		mrisMarkAllDefects(mris, dl, 0) ;
    mrisTessellateDefect(mris, mris_corrected, defect, vertex_trans, mri,
												 h_k1,h_k2,h_white,h_gray, h_border, h_grad, mri_gray_white, h_dot, parms) ;
  }
	HISTOfree(&h_white) ; HISTOfree(&h_gray) ; HISTOfree(&h_border) ; HISTOfree(&h_dot) ;
	HISTOfree(&h_k1) ; HISTOfree(&h_k2) ; HISTOfree(&h_grad) ; MRIfree(&mri_gray_white) ;
  mrisAddAllDefectFaces(mris_corrected, dl, vertex_trans) ;
  mrisCheckSurface(mris_corrected) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    MHT *mht ;
    
    fprintf(stdout, "checking corrected surface for self-intersection...\n") ;
    MRISsaveVertexPositions(mris_corrected, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_corrected, ORIG_VERTICES) ;
    mht = MHTfillTable(mris_corrected, NULL) ;
    MHTcheckSurface(mris_corrected, mht) ;
    MHTfree(&mht) ;
    MRISrestoreVertexPositions(mris_corrected, TMP_VERTICES) ;
  }

#if 0
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    defect = &dl->defects[i] ;
    for (n = 0 ; n < defect->nvertices ; n++)
    {
      vno = vertex_trans[defect->vertices[n]] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (vno < 0)
        continue ;
      v = &mris_corrected->vertices[vno] ;
      if (v->vnum < 2)
      {
        fprintf(stderr, "Warning: vertex %d has only %d neighbors!\n",
                vno, v->vnum) ;
        DiagBreak() ;
      }
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "\n") ;
#endif
  for (vno = 0 ; vno < mris_corrected->nvertices ; vno++)
    mris_corrected->vertices[vno].vtotal = mris_corrected->vertices[vno].vnum ;


#if 0  
  fprintf(stdout, "tessellation finished, orienting corrected surface...\n") ;
  /*  MRISprojectOntoSphere(mris, mris, MRISaverageRadius(mris)) ;*/
  mrisOrientRetessellatedSurface(mris_corrected, dl, vertex_trans) ;
#endif

  fprintf(stdout, "computing original vertex metric properties...\n") ;
  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris_corrected) ;
  fprintf(stdout, "storing new metric properties...\n") ;
  /*  MRISstoreMetricProperties(mris_corrected) ;*/
  fprintf(stdout, "computing tessellation statistics...\n") ;
  MRISprintTessellationStats(mris_corrected, stderr) ;
  MRISrestoreVertexPositions(mris_corrected, CANONICAL_VERTICES) ;

	MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
	MRISsetMarks(mris_corrected, 1) ;
	for (i = 0 ; i < dl->ndefects ; i++)
	{
		defect = &dl->defects[i] ;
		for (n = 0 ; n < defect->nvertices ; n++)
		{
			vno = vertex_trans[defect->vertices[n]] ;
			if (vno < 0 || vno >= mris_corrected->nvertices)
				continue ;
			v = &mris_corrected->vertices[vno] ;
			v->marked = 0 ;
		}
		for (n = 0 ; n < defect->nborder ; n++)
		{
			vno = vertex_trans[defect->border[n]] ;
			if (vno < 0 || vno >= mris_corrected->nvertices)
				continue ;
			v = &mris_corrected->vertices[vno] ;
			v->marked = 0 ;
		}
	}
	fprintf(stdout, 
					"performing soap bubble on retessellated vertices for %d "
					"iterations...\n", nsmooth) ;
	MRISsoapBubbleVertexPositions(mris_corrected, nsmooth) ;
	MRISsaveVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
	MRISclearMarks(mris_corrected) ;

  MRISprintTessellationStats(mris_corrected, stderr) ;
  MRISrestoreVertexPositions(mris_corrected, TMP_VERTICES) ;
  fprintf(stdout, "tessellation finished, orienting corrected surface...\n") ;

  /*  MRISprojectOntoSphere(mris, mris, MRISaverageRadius(mris)) ;*/
  mrisOrientRetessellatedSurface(mris_corrected, dl, vertex_trans) ;
  fprintf(stdout, "final surface representation: %d vertices, %d triangles\n",
          mris->nvertices, mris->nfaces) ;

  free(face_trans) ; free(vertex_trans) ;

  /* free structures */
  for (fno = 0 ; fno < mris->nfaces ; fno++)
    if (fdl->nfaces[fno] > 0)
      free(fdl->faces[fno]) ;
  free(fdl->faces) ; free(fdl->nfaces) ; free(fdl) ;

  for (i = 0 ; i < dl->ndefects ; i++)
  {
    free(dl->defects[i].vertices) ;
    free(dl->defects[i].status) ;
    free(dl->defects[i].border) ;
  }
  free(dl) ;
	if (nmut + ncross > 0)
		printf("%ld mutations (%2.1f%%), %ld crossovers (%2.1f%%)\n",
					 nmut, (float)nmut*100/(nmut+ncross), ncross, (float)ncross*100/(nmut+ncross)) ;
  return(mris_corrected) ;
}
static int
mrisMarkAllDefects(MRI_SURFACE *mris, DEFECT_LIST *dl, int flag)
{
	int j ;

	for (j = 0 ; j < dl->ndefects ; j++)
		mrisMarkDefect(mris, &dl->defects[j], flag) ;
	return(NO_ERROR) ;
}
static int
mrisRipAllDefects(MRI_SURFACE *mris, DEFECT_LIST *dl, int flag)
{
	int j ;

	for (j = 0 ; j < dl->ndefects ; j++)
		mrisRipDefect(mris, &dl->defects[j], flag) ;
	return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Mark all the vertices in the tessellation that are part
          of faces whose centroid intersects other faces. These are
          regions in which the topology is broken as the spherical
          homeomorphism is non-invertible.
------------------------------------------------------*/
FACE_DEFECT_LIST *
MRISmarkAmbiguousVertices(MRI_SURFACE *mris, int mark)
{
  FACE   *f ;
  VERTEX *v ;
  int    fno, flist[1000], i, nfaces, nmarked, n /*, vno, neg*/ ;
  double r, area_scale ;
  MHT    *mht ;
  FILE   *fp = NULL ;
  FDL    *fdl ;

  fdl = (FACE_DEFECT_LIST *)calloc(1, sizeof(FDL)) ;
  if (!fdl)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISmarkAmbiguousFaces: could allocate face defect list") ;
  fdl->nfaces = (int *)calloc(mris->nfaces, sizeof(int)) ;
  if (!fdl->nfaces)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISmarkAmbiguousFaces: could allocate face defect list") ;
  fdl->faces = (int **)calloc(mris->nfaces, sizeof(int *)) ;
  if (!fdl->faces)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISmarkAmbiguousFaces: could allocate face defect list") ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.%s.topology.log", 
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh",
            mris->subject_name) ;
    fp = fopen(fname, "w") ;
  }

  /* the curvature is for diagnostic purposes so I can write it out */
  MRISclearMarks(mris) ; MRISclearCurvature(mris) ;
  mrisMarkBadEdgeVertices(mris, mark) ;
  r = MRISaverageRadius(mris) ;
  MRISscaleBrain(mris, mris, 100.0/r) ;
  mht = MHTfillTable(mris, NULL) ;

  area_scale = mris->orig_area / mris->total_area ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "marking ambiguous vertices...\n") ;
  for (nmarked = fno = 0 ; fno < mris->nfaces ; fno++)
  {
    if (Gdiag & DIAG_SHOW && !(fno % 25000) && fno)
      fprintf(stdout, "%d of %d faces processed, %d ambiguous\n", 
              fno, mris->nfaces-1, nmarked) ;
    f = &mris->faces[fno] ;
    if (fno == Gdiag_no)
      DiagBreak() ;
    if (f->ripflag)
      continue ;
    nfaces = mrisFindAllOverlappingFaces(mris, mht, fno, flist) ;

    /* make sure fno is in list, and add it if it isn't (it should be) */
    for (i = 0 ; i < nfaces ; i++)
      if (flist[i] == fno)
        break ;
    if (i >= nfaces)
      flist[nfaces++] = fno ;

    if ((nfaces > 1 || area_scale*f->area < 0.001) || ((fno <= 5) && Gdiag & DIAG_SAVE_DIAGS))  /* part of a defect */
    {
      nmarked++ ;
#if 0
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "\r%d of %d faces processed, %d ambiguous", 
                fno, mris->nfaces-1, nmarked) ;
#endif


      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      {
        fprintf(stdout, "%d faces @ fno %d\n", nfaces, fno) ;
        for (i = 0 ; i < nfaces ; i++)
        {
          f = &mris->faces[flist[i]] ;
          fprintf(stdout, "\tface %d area %2.4f (%d, %d, %d)\n",
                  flist[i], f->area, f->v[0], f->v[1], f->v[2]) ;
        }
        fprintf(stdout, "\n") ;
      }
      if (Gdiag & DIAG_WRITE && fp != NULL && DIAG_VERBOSE_ON)
      {
        fprintf(fp, "%d faces @ fno %d\n", nfaces, fno) ;
        for (i = 0 ; i < nfaces ; i++)
        {
          f = &mris->faces[flist[i]] ;
          fprintf(fp, "\tface %d area %2.4f (%d, %d, %d)\n",
                  flist[i], f->area, f->v[0], f->v[1], f->v[2]) ;
        }
        fprintf(fp, "\n") ;
        fflush(fp) ;
      }

      fdl->nfaces[fno] = nfaces ;
      fdl->faces[fno] = (int *)calloc(nfaces, sizeof(int)) ;
      if (!fdl->faces[fno])
        ErrorExit(ERROR_NO_MEMORY, 
                  "MRISmarkAmbiguousFaces: could allocate %d defect list",
                  fno) ;
      for (i = 0 ; i < nfaces ; i++)
      {
        fdl->faces[fno][i] = flist[i] ;
        f = &mris->faces[flist[i]] ;
        for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        {
          v = &mris->vertices[f->v[n]] ;
          if (f->v[n] == Gdiag_no)
            DiagBreak() ;
          v->marked = mark ;
        }
      }
    }
  }

#if 0
  /* expand defective vertices outwards by one */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked == mark)
    {
      for (n = 0 ; n < v->vnum ; n++)
        mris->vertices[v->v[n]].marked = mark+1 ;
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked == mark+1)
      v->marked = mark ;
  }
#endif

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "\n") ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON && fp)
    fclose(fp) ;
  MRISscaleBrain(mris, mris, r/100.0) ;
  fprintf(stdout, " %d ambiguous faces found in tessellation\n", nmarked) ;
  MHTfree(&mht) ;
  return(fdl) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Find all connected components of a defect.
------------------------------------------------------*/
DEFECT_LIST *
MRISsegmentDefects(MRI_SURFACE *mris, int mark_ambiguous, int mark_segmented)
{
  DEFECT_LIST  *dl ;
  int          vno ;
  VERTEX       *v ;
  FILE         *fp = NULL ;

  dl = (DEFECT_LIST *)calloc(1, sizeof(DEFECT_LIST)) ;
  if (!dl)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISsegmentDefects: could allocate defect list") ;

  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.%s.topology.log", 
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh",
            mris->subject_name) ;
    fp = fopen(fname, "a") ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked != mark_ambiguous)
      continue ;
    if (dl->ndefects == Gdiag_no)
      DiagBreak() ;
    mrisSegmentDefect(mris, vno, &dl->defects[dl->ndefects++], 
                      mark_ambiguous, mark_segmented) ;
#if 0
    if (Gdiag & DIAG_SHOW)
    {
      DEFECT *defect = &dl->defects[dl->ndefects-1] ;
#if 0
      fprintf(stdout, "\rdefect %3d: %4d vertices, area %2.2f c (%2.1f,%2.1f,%2.1f), n (%2.1f,%2.1f,%2.1f)  ",
              dl->ndefects, defect->nvertices, defect->area, 
              defect->cx, defect->cy, defect->cz,
              defect->nx, defect->ny, defect->nz) ;
#else
      fprintf(stdout, "\rdefect %3d: %4d vertices, area %2.2f, "
              "c (%2.1f,%2.1f,%2.1f)  ", dl->ndefects, defect->nvertices, 
              defect->area, defect->cx, defect->cy, defect->cz) ;
#endif
    }
#endif
    if (Gdiag & DIAG_WRITE && fp)
    {
      int n ;
      DEFECT *defect = &dl->defects[dl->ndefects-1] ;
      fprintf(fp, "defect %d found with %d vertices, area %2.2f\n"
              "\tcentroid (%2.1f,%2.1f,%2.1f)\n",
              dl->ndefects, defect->nvertices, defect->area, 
              defect->cx, defect->cy, defect->cz) ;
      for (n = 0 ; n < defect->nvertices ; n++)
        fprintf(fp, "\t%d\n", defect->vertices[n]) ;
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "\n") ;
  if (Gdiag & DIAG_WRITE && fp)
    fclose(fp) ;
  return(dl) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Segment the connected region of a defect, starting with vno
           and spreading outwards.
------------------------------------------------------*/
static int
mrisSegmentDefect(MRI_SURFACE *mris, int vno, DEFECT *defect,
                  int mark_ambiguous, int mark_segmented)
{
  int    vlist[200000], i, n, nfilled ;
  VERTEX *v, *vn ;
  float  len, nx, ny, nz ;

  vlist[defect->nvertices++] = vno ;  /* start the list */

  v = &mris->vertices[vno] ;
  v->marked = mark_segmented ;
  defect->cx = v->x ;  defect->cy = v->y ;  defect->cz = v->z ; 
  defect->area = v->origarea ;
  do
  {
    nfilled = 0 ;

    for (i = 0 ; i < defect->nvertices ; i++)
    {
      v = &mris->vertices[vlist[i]] ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (v->v[n] == Gdiag_no)
          DiagBreak() ;
        if (vn->marked == mark_ambiguous)
        {
          vlist[defect->nvertices++] = v->v[n] ;  /* add it to list */

          vn->marked = mark_segmented ;
          defect->cx += vn->x ;  defect->cy += vn->y ;  defect->cz += vn->z ; 
          defect->area += vn->origarea ;
          nfilled++ ;
        }
      }
    }
  } while (nfilled > 0) ;

  defect->cx /= (float)defect->nvertices ;
  defect->cy /= (float)defect->nvertices ;
  defect->cz /= (float)defect->nvertices ;
  defect->vertices = (int *)calloc(defect->nvertices, sizeof(int)) ;
  if (!defect->vertices)
    ErrorExit(ERROR_NO_MEMORY, 
              "mrisSegmentDefect: could allocate defect vertex list") ;
  defect->status = (char *)calloc(defect->nvertices, sizeof(char)) ;
  if (!defect->status)
    ErrorExit(ERROR_NO_MEMORY, 
              "mrisSegmentDefect: could allocate defect status list") ;
  memcpy(defect->vertices, vlist, defect->nvertices*sizeof(int)) ;
  for (nfilled = i = 0 ; i < defect->nvertices ; i++)
  {
    v = &mris->vertices[defect->vertices[i]] ;
    if (defect->vertices[i] == Gdiag_no)
      DiagBreak() ;
    v->val = defect->area ;
    defect->status[i] = KEEP_VERTEX ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      if (v->v[n] == Gdiag_no)
        DiagBreak() ;
      if (mris->vertices[v->v[n]].marked == 0) /* border vertex */
      {
        mris->vertices[v->v[n]].marked = 2 ;
        vlist[nfilled++] = v->v[n] ;
      }
    }
  }

  defect->border = (int *)calloc(nfilled, sizeof(int)) ;
  defect->nborder = nfilled ;
  memcpy(defect->border, vlist, defect->nborder*sizeof(int)) ;
  mrisMarkDefectBorder(mris, defect, 0) ;

  nx = ny = nz = 0.0f ;
  for (n = 0 ; n < defect->nborder ; n++)
  {
    v = &mris->vertices[defect->border[n]] ;
    nx += v->nx ; ny += v->ny ; nz += v->nz ;
  }
  len = sqrt(nx*nx + ny*ny + nz*nz) ;
  if (FZERO(len))
    len = 1.0f ;
  defect->nx = nx/len ; defect->ny = ny/len ; defect->nz = nz/len ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Mark the vertices in a defect as either retained or
          discarded. The decision will be made based on the surface 
          area in the defect. If it is above some threshold, then
          fill it by marking the "outside" faces as kept and all others
          as discarded. Otherwise cut it by marking the "inside" face
          as kept and all others as discarded. The inside/outside
          decision will be made by using the dot product of the average
          inflated surface normal with the inflated face centroid.
------------------------------------------------------*/
#define MIN_SPHERE_DIST   .01
#define MIN_ORIG_DIST     .75

#if 1
int
mrisMarkRetainedPartOfDefect(MRI_SURFACE *mris, DEFECT *defect, 
                             FACE_DEFECT_LIST *fdl, float area_threshold, 
                             int mark_retain, int mark_discard, MHT *mht)
{
#if 0
  int      n, i, j, nfaces, fno, flist[100000], n2, vno, retain ;
  FACE     *f ;
  VERTEX   *v, *vn ;
  float    dot, x0, y0, z0, x, y, z, dx, dy, dz, dist, fn, dot0, len ;
#endif

  mrisMarkDefect(mris, defect, 0) ;
  mrisDefectRemoveDegenerateVertices(mris, MIN_SPHERE_DIST, defect) ;
  mrisDefectRemoveProximalVertices(mris, MIN_ORIG_DIST, defect) ;
  mrisDefectRemoveNegativeVertices(mris, defect) ;

#if 0
  /* throw out anything in a negative face */
  for (i = 0 ; i < defect->nvertices ; i++)
  {
    if (defect->status[i] == DISCARD_VERTEX)
      continue ;
    v = &mris->vertices[defect->vertices[i]] ;
    for (n = 0 ; n < v->num ; n++)
      if (mris->faces[v->f[n]].area < 0.05)
        defect->status[i] = DISCARD_VERTEX ;
  }

  /* compute centroid and average normal of defect using border vertices */
  defect->cx = defect->cy = defect->cz = 0.0f ;
  defect->nx = defect->ny = defect->nz = 0.0f ;
  for (fn = 0.0f, i = 0 ; i < defect->nborder ; i++, fn += 1.0f)
  {
    v = &mris->vertices[defect->border[i]] ;
    defect->nx += v->nx ; defect->ny += v->ny ; defect->nz += v->nz ;
    defect->cx += v->x  ; defect->cy += v->y  ; defect->cz += v->z ;
  }
  len = sqrt(SQR(defect->nx) + SQR(defect->ny) + SQR(defect->nz)) ;
  defect->nx /= len ; defect->ny /= len ; defect->nz /= len ;
  defect->cx /= fn  ; defect->cy /= fn  ; defect->cz /= fn  ;

  /* discard vertices that are too close to another vertex */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    float  dx, dy, dz ;

    if (i < defect->nvertices)
    {
      if (defect->status[i] == DISCARD_VERTEX)
        continue ;
      v = &mris->vertices[defect->vertices[i]] ;
    }
    else
      v = &mris->vertices[defect->border[i-defect->nvertices]] ;
    for (j = i+1 ; j < defect->nvertices ; j++)
    {
      if (defect->status[j] == DISCARD_VERTEX)
        continue ;
      vn = &mris->vertices[defect->vertices[j]] ;
      dx = vn->origx-v->origx ;
      dy = vn->origy-v->origy ;
      dz = vn->origz-v->origz ;
      dist = sqrt(dx*dx+dy*dy+dz*dz) ;
      if (dist <= 0.75)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "discarding proximal vertex %d\n",
                  defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX ;
        vn->imag_val = -1.0f ;
      }
    }
  }

  /* build a list of all faces in this defect */
  for (nfaces = n = 0 ; n < defect->nvertices ; n++)
  {
    if (defect->status[n] == DISCARD_VERTEX)
      continue ;

    vno = defect->vertices[n] ; v = &mris->vertices[vno] ;

    /* build a list of faces */
    for (n2 = 0 ; n2 < v->num ; n2++)
    {
      if (mris->faces[v->f[n2]].ripflag == 0)
      {
        if (nfaces == 7127)
          DiagBreak() ;
        flist[nfaces++] = v->f[n2] ;
        mris->faces[v->f[n2]].ripflag = 1 ;  /* temporary */
      }
    }
  }

  /* for really big defects throw out 'inside' vertices */
  for (n = 0 ; n < nfaces ; n++)
  {
    fno = flist[n] ; f = &mris->faces[fno] ;
    mrisCalculateFaceCentroid(mris, fno, &x0, &y0, &z0) ;
    dx = x0 - defect->cx ; dy = y0 - defect->cy ; dz = z0 - defect->cz ; 
    dot0 = dx*defect->nx + dy*defect->ny + dz*defect->nz ;

    /* see if there are any faces inside (outside) of this one */
    retain = 1 ;
    for (n2 = 0 ; n2 < fdl->nfaces[fno] ; n2++)
    {
      if (triangleNeighbors(mris, fno, fdl->faces[fno][n2]) >= 1)
        continue ;
      mrisCalculateFaceCentroid(mris, fdl->faces[fno][n2], &x, &y, &z);
      dx = x - defect->cx ; dy = y - defect->cy ; dz = z - defect->cz ; 
      dot = dx*defect->nx + dy*defect->ny + dz*defect->nz ;
#define HUGE_DEFECT 10000
#define BIG_DEFECT   5000
      if ((defect->nvertices > HUGE_DEFECT) && (dot > dot0))
      {
        retain = 0 ;
        break ;   /* found a face outside of this one - discard it */
      }
      if (defect->nvertices > BIG_DEFECT && defect->nvertices < HUGE_DEFECT)
      {    
        if (dot < dot0)  /* found a face inside this one - keep it */
        {
          retain = 1 ;
          break ;
        }
        else
          retain =  0 ;
      }
    }
    if (!retain)  /* no faces outside of this one */
    {
      for (n2 = 0 ; n2 < VERTICES_PER_FACE ; n2++)
      {
        if (f->v[n2] == 1245 || f->v[n2] == Gdiag_no)
          DiagBreak() ;
        mris->vertices[f->v[n2]].marked = 1 ;
        mris->vertices[f->v[n2]].imag_val = -1 ;
      }
    }
  }

  /* discard all marked vertices */
  for (i = 0 ; i < defect->nvertices ; i++)
  {
    if (defect->status[i] == DISCARD_VERTEX)
      continue ;
    v = &mris->vertices[defect->vertices[i]] ;
    if (v->marked)
    {
      v->marked = 0 ;
      defect->status[i] = DISCARD_VERTEX ;
    }
  }

  /* unmark the faces */
  for (n = 0 ; n < nfaces ; n++)
  {
    fno = flist[n] ; f = &mris->faces[fno] ;
    f->ripflag = 0 ;
  }

  /* discard vertices that are too close to another vertex */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    float  dx, dy, dz ;

    if (i < defect->nvertices)
    {
      if (defect->status[i] == DISCARD_VERTEX)
        continue ;
      v = &mris->vertices[defect->vertices[i]] ;
    }
    else
      v = &mris->vertices[defect->border[i-defect->nvertices]] ;
    for (j = i+1 ; j < defect->nvertices ; j++)
    {
      if (defect->status[j] == DISCARD_VERTEX)
        continue ;
      vn = &mris->vertices[defect->vertices[j]] ;
      dx = vn->cx-v->cx ;
      dy = vn->cy-v->cy ;
      dz = vn->cz-v->cz ;
      dist = (dx*dx+dy*dy+dz*dz) ;  /* no sqrt */
      if (dist < MIN_SPHERE_DIST)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "discarding proximal vertex %d\n",
                  defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX ;
        vn->imag_val = -1.0f ;
      }
    }
  }
#endif
  return(NO_ERROR) ;
}
#else
static int
mrisMarkRetainedPartOfDefect(MRI_SURFACE *mris, DEFECT *defect, 
                             FACE_DEFECT_LIST *fdl, float area_threshold, 
                             int mark_retain, int mark_discard, MHT *mht)
{
#define USING_CUBE 0
#if USING_CUBE || 1
  /* based on faces */
  int     n, n2, inside, vno, flist[100000], nfaces, fno, i, j ;
  VERTEX  *v, *vn ;
  FACE    *f ;
  float   dot, x0, y0, z0, x, y, z, dx, dy, dz, fn, len, dist ;

  defect->cx = defect->cy = defect->cz = 0.0f ;
  inside = defect->area < area_threshold ;
  mrisMarkDefect(mris, defect, 1) ;
  for (fn = 0.0, nfaces = n = 0 ; n < defect->nvertices ; n++)
  {
    vno = defect->vertices[n] ; v = &mris->vertices[vno] ;

    /* build a list of faces */
    for (n2 = 0 ; n2 < v->num ; n2++)
    {
      if (mris->faces[v->f[n2]].ripflag == 0)
      {
        if (nfaces == 7127)
          DiagBreak() ;
        flist[nfaces++] = v->f[n2] ;
        mris->faces[v->f[n2]].ripflag = 1 ;  /* temporary */
      }
    }
    v->val = inside ? -1.0f : 1.0 ;
    v->imag_val = 1.0f ;  /* assume kept until found otherwise */
    defect->cx += v->x ;  defect->cy += v->y ;  defect->cz += v->z ;

  }
  defect->nx = defect->ny = defect->nz = 0.0f ;
  for (n = 0 ; n < defect->nborder ; n++)
  {
    vno = defect->border[n] ; v = &mris->vertices[vno] ;
    defect->nx += v->nx ; defect->ny += v->ny ; defect->nz += v->nz ;
  }
  mrisMarkDefect(mris, defect, 0) ;
  len =sqrt(defect->nx*defect->nx+defect->ny*defect->ny+defect->nz*defect->nz);
  if (FZERO(len))
    len = 1.0f ;
  defect->nx /= len ; defect->ny /= len ; defect->nz /= len ;
  fn = (float)defect->nvertices ;
  defect->cx /= fn ; defect->cy /= fn ; defect->cz /= fn ;

  /* unrip the faces (used really as a mark, but don't want to add to
     face struct).
  */
  for (n = 0 ; n < defect->nvertices ; n++)
  {
    v = &mris->vertices[defect->vertices[n]] ;

    for (n2 = 0 ; n2 < v->num ; n2++)
      mris->faces[v->f[n2]].ripflag = 0 ;
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "%d faces found in defect\n", nfaces) ;
#if USING_CUBE
  for (n = 0 ; n < nfaces ; n++)
  {
    fno = flist[n] ; f = &mris->faces[fno] ;
    mrisCalculateFaceCentroid(mris, fno, &x0, &y0, &z0) ;

    /* see if there are any faces inside (outside) of this one */
    for (dot = 0.0f, n2 = 0 ; n2 < fdl->nfaces[fno] ; n2++)
    {
      if (triangleNeighbors(mris, fno, fdl->faces[fno][n2]) >= 1)
        continue ;
      mrisCalculateFaceCentroid(mris, fdl->faces[fno][n2], &x, &y, &z);
      dx = -15.5 - x ; dot = dx*dx ;
      if (dot < ((-15.5-x0)*(-15.5-x0)))
        break ;
    }
    if (n2 < fdl->nfaces[fno])  /* found a face inside (outside) of this one */
    {
      for (n2 = 0 ; n2 < VERTICES_PER_FACE ; n2++)
      {
        if (f->v[n2] == 1245 || f->v[n2] == Gdiag_no)
          DiagBreak() ;
        mris->vertices[f->v[n2]].marked = 1 ;
        mris->vertices[f->v[n2]].imag_val = -1 ;
      }
    }
  }
#else
  for (n = 0 ; n < nfaces ; n++)
  {
    fno = flist[n] ; f = &mris->faces[fno] ;
    mrisCalculateFaceCentroid(mris, fno, &x0, &y0, &z0) ;

    /* see if there are any faces inside (outside) of this one */
    for (dot = 0.0f, n2 = 0 ; n2 < fdl->nfaces[fno] ; n2++)
    {
      if (triangleNeighbors(mris, fno, fdl->faces[fno][n2]) >= 1)
        continue ;
      mrisCalculateFaceCentroid(mris, fdl->faces[fno][n2], &x, &y, &z);
      dx = x - x0 ; dy = y - y0 ; dz = z - z0 ; 
      dot = dx*defect->nx + dy*defect->ny + dz*defect->nz ;
      if ((inside && dot < 0) || (!inside && dot > 0))
        break ;
    }
    if (n2 < fdl->nfaces[fno])  /* found a face inside (outside) of this one */
    {
      for (n2 = 0 ; n2 < VERTICES_PER_FACE ; n2++)
      {
        if (f->v[n2] == 1245 || f->v[n2] == Gdiag_no)
          DiagBreak() ;
        mris->vertices[f->v[n2]].marked = 1 ;
        mris->vertices[f->v[n2]].imag_val = -1 ;
      }
    }
  }
#endif
  for (n = 0 ; n < defect->nvertices ; n++)
  {
    vno = defect->vertices[n] ; v = &mris->vertices[vno] ;
    defect->status[n] = v->marked ? DISCARD_VERTEX : KEEP_VERTEX ;
  }      
  mrisMarkDefect(mris, defect, 0) ;
#else
  int     n, n2, inside, vno ;
  VERTEX  *v, *v2 ;
  float   fn, len ;

  defect->cx=defect->cy=defect->cz=defect->nx=defect->ny=defect->nz = 0.0 ;
  inside = defect->area < area_threshold ;
  mrisMarkDefect(mris, defect, 1) ;
  for (fn = 0.0, n = 0 ; n < defect->nvertices ; n++)
  {
    vno = defect->vertices[n] ; v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (inside)
      v->val = -1.0f ;
    else
      v->val = 1.0f ;
    defect->cx += v->x ;  defect->cy += v->y ;  defect->cz += v->z ;

    /* use surrounding unmarked vertices as estimate of local normal */
    for (n2 = 0 ; n2 < v->vnum ; n2++)
    {
      v2 = &mris->vertices[v->v[n2]] ;
      if (!v2->marked)
      {
        defect->nx += v2->nx ; defect->ny += v2->ny ; defect->nz += v2->nz ;
        fn += 1.0 ;
      }
    }
  }
  mrisMarkDefect(mris, defect, 0) ;
  defect->nx /= fn ; defect->ny /= fn ; defect->nz /= fn ;
  len =sqrt(defect->nx*defect->nx+defect->ny*defect->ny+defect->nz*defect->nz);
  if (FZERO(len))
    len = 1.0f ;
  defect->nx /= len ; defect->ny /= len ; defect->nz /= len ;
  fn = (float)defect->nvertices ;
  defect->cx /= fn ; defect->cy /= fn ; defect->cz /= fn ;

  for (n = 0 ; n < defect->nvertices ; n++)  /* for each vertex in defect */
  {
    vno = defect->vertices[n] ; v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    defect->status[n] = KEEP_VERTEX ;
    v->imag_val = 1.0 ; 
    v->nx = defect->nx ; v->ny = defect->ny ; v->nz = defect->nz ; 
    if (inside)
    { 
      if (mrisFindNextInwardFace(mris, mht, vno, 20.0f) >= 0)
      {
        defect->status[n] = DISCARD_VERTEX ;
        defect->status[n] = KEEP_VERTEX ;
        v->imag_val = -1.0 ;
      }
    }
    else 
    { 
      if (mrisFindNextOutwardFace(mris, mht, vno, 20.0f) >= 0)
      {
        defect->status[n] = DISCARD_VERTEX ;
        v->imag_val = -1.0 ;
      }
    }
  }
#endif

  /* throw out anything in a negative face */
  for (i = 0 ; i < defect->nvertices ; i++)
  {
    if (defect->status[i] == DISCARD_VERTEX)
      continue ;
    v = &mris->vertices[defect->vertices[i]] ;
    for (n = 0 ; n < v->num ; n++)
      if (mris->faces[v->f[n]].area < 0)
        defect->status[i] = DISCARD_VERTEX ;
  }

  /* discard vertices that are too close to another vertex */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    float  dx, dy, dz ;

    if (i < defect->nvertices)
    {
      if (defect->status[i] == DISCARD_VERTEX)
        continue ;
      v = &mris->vertices[defect->vertices[i]] ;
    }
    else
      v = &mris->vertices[defect->border[i-defect->nvertices]] ;
    for (j = i+1 ; j < defect->nvertices ; j++)
    {
      if (defect->status[j] == DISCARD_VERTEX)
        continue ;
      vn = &mris->vertices[defect->vertices[j]] ;
      dx = vn->origx-v->origx ;
      dy = vn->origy-v->origy ;
      dz = vn->origz-v->origz ;
      dist = sqrt(dx*dx+dy*dy+dz*dz) ;
      if (dist <= 0.5)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "discarding proximal vertex %d\n",
                  defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX ;
        vn->imag_val = -1.0f ;
      }
    }
  }
  /* discard vertices that are too close to another vertex on sphere */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    float  dx, dy, dz ;

    if (i < defect->nvertices)
    {
      if (defect->status[i] == DISCARD_VERTEX)
        continue ;
      v = &mris->vertices[defect->vertices[i]] ;
    }
    else
      v = &mris->vertices[defect->border[i-defect->nvertices]] ;
    for (j = i+1 ; j < defect->nvertices ; j++)
    {
      if (defect->status[j] == DISCARD_VERTEX)
        continue ;
      vn = &mris->vertices[defect->vertices[j]] ;
      dx = vn->cx-v->cx ;
      dy = vn->cy-v->cy ;
      dz = vn->cz-v->cz ;
      dist = (dx*dx+dy*dy+dz*dz) ;  /* no sqrt */
      if (dist < MIN_SPHERE_DIST)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "discarding proximal vertex %d\n",
                  defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX ;
        vn->imag_val = -1.0f ;
      }
    }
  }
  return(NO_ERROR) ;
}
#endif
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Is vno a vertex of face fno?
------------------------------------------------------*/
static int
vertexInFace(MRI_SURFACE *mris, int vno, int fno)
{
  VERTEX   *v ;
  int      n ;

  v = &mris->vertices[vno] ;
  for (n = 0 ; n < v->num ; n++)
    if (v->f[n] == fno)
      return(1) ;
  return(0) ;
}
#endif
static int
mrisRipDefect(MRI_SURFACE *mris, DEFECT *defect, int ripflag)
{
  int   n ;

  for (n = 0 ; n < defect->nvertices ; n++)
    mris->vertices[defect->vertices[n]].ripflag = ripflag ;

  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Mark all the vertices in the given defect.
------------------------------------------------------*/
static int
mrisMarkDefect(MRI_SURFACE *mris, DEFECT *defect, int mark)
{
  int   n ;

  for (n = 0 ; n < defect->nvertices ; n++)
    mris->vertices[defect->vertices[n]].marked = mark ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Mark all the border vertices in the given defect.
------------------------------------------------------*/
static int
mrisMarkDefectConvexHull(MRI_SURFACE *mris, DEFECT *defect, int mark)
{
  int   n ;

  for (n = 0 ; n < defect->nchull ; n++)
    mris->vertices[defect->chull[n]].marked = mark ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Mark all the border vertices in the given defect.
------------------------------------------------------*/
static int       
mrisMarkDefectBorder(MRI_SURFACE *mris, DEFECT *defect,int mark)
{
  int   n ;

  for (n = 0 ; n < defect->nborder ; n++)
    mris->vertices[defect->border[n]].marked = mark ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale the surface so that it's max dimension has length maxr
------------------------------------------------------*/
static int
mrisScaleMaxDimension(MRI_SURFACE *mris, float maxr)
{
  maxr /= 2.0 ;
  mrisComputeSurfaceDimensions(mris) ;
  if (mris->xhi >= maxr)
    MRISscaleBrain(mris, mris, maxr / mris->xhi) ;
  if (mris->yhi >= maxr)
    MRISscaleBrain(mris, mris, maxr / mris->yhi) ;
  if (mris->yhi >= maxr)
    MRISscaleBrain(mris, mris, maxr / mris->yhi) ;
  if (fabs(mris->xlo) >= maxr)
    MRISscaleBrain(mris, mris, maxr / fabs(mris->xlo)) ;
  if (fabs(mris->ylo) >= maxr)
    MRISscaleBrain(mris, mris, maxr / fabs(mris->ylo)) ;
  if (fabs(mris->zlo) >= maxr)
    MRISscaleBrain(mris, mris, maxr / fabs(mris->zlo)) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          See if any of the vertices in a triangle are marked
------------------------------------------------------*/
static int
triangleMarked(MRI_SURFACE *mris, int fno)
{
  int  n ;
  FACE *f ;

  f = &mris->faces[fno] ;
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    if (mris->vertices[f->v[n]].marked != 0)
      return(1) ;
  }
  return(0) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Return count of # of vertices shared by 2 triangles
------------------------------------------------------*/
static int
triangleNeighbors(MRI_SURFACE *mris, int fno1, int fno2)
{
  int  n1, n2, num ;
  FACE *f1, *f2 ;

  f1 = &mris->faces[fno1] ; f2 = &mris->faces[fno2] ;
  for (num = n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
  {
    for (n2 = 0 ; n2 < VERTICES_PER_FACE ; n2++)
      if (f1->v[n1] == f2->v[n2])
        num++ ;
  }
  return(num) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Tessellate a defect using the spherical space to compute
          edge intersection, and the original space to compute edge length.
------------------------------------------------------*/
static int defect_no = 0 ;

static int compare_edge_length(const void *pe1, const void *pe2) ;
static int edgeExists(MRI_SURFACE *mris, int vno1, int vno2) ;
static int mrisComputeDefectVertexNormal(MRI_SURFACE *mris, int vno, 
																				 double *pnx, double *pny, double *pnz) ;

static int
edgeExists(MRI_SURFACE *mris, int vno1, int vno2)
{
  int    n ;
  VERTEX *v ;
  
  v = &mris->vertices[vno1] ;
  for (n = 0 ; n < v->vnum ; n++)
    if (v->v[n] == vno2)
      return(1) ;
  return(0) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisComputeDefectTangentPlanes(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans)
{
  VECTOR  *v_n, *v_e1, *v_e2, *v ;
  int     vno, i ;
  VERTEX  *vertex ;
	double  nx, ny, nz ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v = VectorAlloc(3, MATRIX_REAL) ;

  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
		if (i < defect->nvertices)
			vno = vertex_trans[defect->vertices[i]] ;
		else
			vno = vertex_trans[defect->border[i-defect->nvertices]] ;
		if (vno < 0)
			continue ;
    vertex = &mris->vertices[vno] ;
		mrisComputeDefectVertexNormal(mris, vno, &nx, &ny, &nz) ;
		vertex->nx = nx ; vertex->ny = ny ; vertex->nz = nz ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    /* now find some other non-parallel vector */
#if 0
    if (!FZERO(vertex->nx) || !FZERO(vertex->ny))
    {VECTOR_LOAD(v, 0.0, 0.0, 1.0) ; }
    else
    {VECTOR_LOAD(v, 0.0, 1.0, 0.0) ; }
#else
    VECTOR_LOAD(v, vertex->ny, vertex->nz, vertex->nx) ;
#endif
    V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    if ((V3_LEN_IS_ZERO(v_e1)))  /* happened to pick a parallel vector */
    {
      VECTOR_LOAD(v, vertex->ny, -vertex->nz, vertex->nx) ;
      V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    }

    if ((V3_LEN_IS_ZERO(v_e1)) && DIAG_VERBOSE_ON)  /* happened to pick a parallel vector */
      fprintf(stderr, "vertex %d: degenerate tangent plane\n", vno) ;
    V3_CROSS_PRODUCT(v_n, v_e1, v_e2) ;
    V3_NORMALIZE(v_e1, v_e1) ;
    V3_NORMALIZE(v_e2, v_e2) ;
    vertex->e1x = V3_X(v_e1) ; vertex->e2x = V3_X(v_e2) ;
    vertex->e1y = V3_Y(v_e1) ; vertex->e2y = V3_Y(v_e2) ;
    vertex->e1z = V3_Z(v_e1) ; vertex->e2z = V3_Z(v_e2) ;
  }

  VectorFree(&v) ;
  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  return(NO_ERROR) ;
}
static double 
mrisComputeDefectCurvatureLogLikelihood(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans, 
																				HISTOGRAM *h_k1, HISTOGRAM *h_k2)
{
	double   total_ll = 0.0 ;
  int      vno, i ;
  VERTEX   *v ;
	int      nbrs[MAX_NBRS], num_nbrs, bin ;

	MRISsaveVertexPositions(mris, TMP_VERTICES) ;
	MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
	/*	MRIScomputeMetricProperties(mris) ;*/

  mrisComputeDefectTangentPlanes(mris, defect, vertex_trans) ;
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
		if (i < defect->nvertices)
			vno = vertex_trans[defect->vertices[i]] ;
		else
			vno = vertex_trans[defect->border[i-defect->nvertices]] ;
		if (vno < 0)
			continue ;
		if (vno == Gdiag_no)
			DiagBreak() ;
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->vnum <= 0)
      continue ;
		mrisFindSecondNeighborhood(mris, vno, nbrs, &num_nbrs) ;
		if (num_nbrs < 3)
			continue ;
		MRIScomputeSecondFundamentalFormAtVertex(mris, vno, nbrs, num_nbrs) ;
		bin = nint((v->k1-h_k1->bins[0]) / h_k1->bin_size) ;
		if (bin < 0)
			bin = 0 ;
		else if (bin >= h_k1->nbins)
			bin = h_k1->nbins-1 ;
		total_ll += log(h_k1->counts[bin]) ;
		bin = nint((v->k2-h_k2->bins[0]) / h_k2->bin_size) ;
		if (bin < 0)
			bin = 0 ;
		else if (bin >= h_k2->nbins)
			bin = h_k2->nbins-1 ;
		total_ll += log(h_k2->counts[bin]) ;
		if (!finite(total_ll))
			DiagBreak() ;
  }

	MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
	return(total_ll) ;
}

static double 
mrisComputeDefectNormalDotLogLikelihood(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans, 
																				HISTOGRAM *h_dot)
{
	double   total_ll = 0.0, nx, ny, nz, x, y, z, dx, dy, dz, ll, dot ;
  int      vno, n, i, bin, nvertices ;
  VERTEX   *v, *vn ;

  for (nvertices = i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
		if (i < defect->nvertices)
			vno = vertex_trans[defect->vertices[i]] ;
		else
			vno = vertex_trans[defect->border[i-defect->nvertices]] ;
		if (vno < 0)
			continue ;
		if (vno == Gdiag_no)
			DiagBreak() ;
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->vnum <= 0)
      continue ;

		nvertices++ ;
		mrisComputeDefectVertexNormal(mris, vno, &nx, &ny, &nz) ;
		x = v->origx ; y = v->origy ; z = v->origz ;
		for (ll = 0.0, n = 0 ; n < v->vnum ; n++)
		{
			vn = &mris->vertices[v->v[n]] ;
			dx = vn->origx - x ; dy = vn->origy - y ; dz = vn->origz - z ; 
			dot = dx*nx + dy*nz * dz*nz ;
			bin = nint((dot-h_dot->bins[0]) / h_dot->bin_size) ;
			if (bin < 0)
				bin = 0 ;
			else if (bin >= h_dot->nbins)
				bin = h_dot->nbins-1 ;
			ll += log(h_dot->counts[bin]) ;
		}
		total_ll += ll / v->vnum ;
		if (!finite(total_ll))
			DiagBreak() ;
  }

	return(total_ll/*/(float)nvertices*/) ;
}

static double 
mrisComputeDefectMRILogUnlikelihood(MRI_SURFACE *mris, MRI *mri, DEFECT_PATCH *dp, 
																		int *vertex_trans, HISTOGRAM *h_border)
{
  double dx, dy, dz, d, len, ll, total_ll ;
  int    i, nedges = 0, nsamples ;
  VERTEX *v, *vn ;
  Real   val, xv, yv, zv, x, y, z ;
  EDGE    *edge ;
  EDGE_TABLE *etable = dp->etable ;

  for (total_ll = 0.0, i = 0 ; i < dp->nedges ; i++)
  {
    edge = &etable->edges[i] ;
    if (edge->used)   /* only count edges not in current tessellation */
      continue ;
    v = &mris->vertices[edge->vno1] ;
    vn = &mris->vertices[edge->vno2] ;
		
    if (v->ripflag || vn->ripflag)
      continue ;

    if (edge->vno1 == Gdiag_no || edge->vno2 == Gdiag_no)
      DiagBreak() ;
    nedges++ ;

    /* sample MR values along line and build estimate of log likelihood
       as distance.
    */
    
    /* traverse the edge connecting the two vertices in equal increments, sampling
       MRI volume outside and inside */
    dx = vn->origx - v->origx ; dy = vn->origy - v->origy ; dz = vn->origz - v->origz ;
    len = sqrt(dx*dx + dy*dy + dz*dz) ;
    len = .5 / len ;  /* sample every 1/2 mm */
    if (len > 1)
      len = 1 ;
    
    ll = 0.0 ;
    for (nsamples = 0, d = 0 ; d <= 1 ; d += len, nsamples++)
    {
      x = v->origx+d*dx ; y = v->origy+d*dy ; z = v->origz+d*dz ;
      //MRIworldToVoxel(mri, x, y, z, &xv, &yv, &zv) ;
      MRIsurfaceRASToVoxel(mri, x, y, z, &xv, &yv, &zv) ;
      MRIsampleVolume(mri, xv, yv, zv, &val) ;
      ll += log(1-h_border->counts[nint(val)]) ;
    }
    
    if (!finite(ll))
      DiagBreak() ;
    total_ll += ll / (float)nsamples ;
  }
  return(total_ll/*/(double)nedges*/) ;
}
		 
static float mrisDefectFaceMRILogLikelihood(MRI_SURFACE *mris, int vno1, int vno2, 
					    int vno3, MRI *mri, DEFECT *defect, 
					    HISTOGRAM *h_white, HISTOGRAM *h_gray, 
					    HISTOGRAM *h_grad, MRI *mri_gray_white) ;
static float
mrisDefectFaceMRILogLikelihood(MRI_SURFACE *mris, int vno0, int vno1, int vno2, 
			       MRI *mri, DEFECT *defect, HISTOGRAM *h_white, 
			       HISTOGRAM *h_gray, HISTOGRAM *h_grad, MRI *mri_gray_white)
{
  Real   x, y, z, xa, ya, za, xc, yc, zc, t0, t1, adx, ady, adz, dx, dy, dz, grad,
         cdx, cdy, cdz, alen, clen, delta_t0, delta_t1, len, nx, ny, nz, xv, yv, zv,
         white_val, gray_val, cnx, cny, cnz, dot, val ;
	double ll = 0.0, jll = 0.0 ; 
  int    bin, nsamples ;
  VERTEX *v0, *v1, *v2 ;
  float   l0[3],l1[3];

	
  v0 = &mris->vertices[vno0] ; 
  v1 = &mris->vertices[vno1] ; 
  v2 = &mris->vertices[vno2] ; 
  l0[0] = v0->origx - v1->origx; l0[1] = v0->origy - v1->origy; 
  l0[2] = v0->origz - v1->origz;
  l1[0] = v2->origx - v0->origx; l1[1] = v2->origy - v0->origy; 
  l1[2] = v2->origz - v0->origz;
  mrisNormalize(l0); mrisNormalize(l1);
  nx = -l1[1]*l0[2] + l0[1]*l1[2];
  ny = l1[0]*l0[2] - l0[0]*l1[2];
  nz = -l1[0]*l0[1] + l0[0]*l1[1];
  alen = nx*nx + ny*nz + nz*nz ; if (FZERO(alen)) alen = 1 ;
  nx /= alen ; ny /= alen ; nz /= alen ; 

	/* now use sphere to orient normal vector */
  l0[0] = v0->cx - v1->cx; l0[1] = v0->cy - v1->cy; 
  l0[2] = v0->cz - v1->cz;
  l1[0] = v2->cx - v0->cx; l1[1] = v2->cy - v0->cy; 
  l1[2] = v2->cz - v0->cz;
  mrisNormalize(l0); mrisNormalize(l1);
  cnx = -l1[1]*l0[2] + l0[1]*l1[2];
  cny = l1[0]*l0[2] - l0[0]*l1[2];
  cnz = -l1[0]*l0[1] + l0[0]*l1[1];
  dot = v0->cx*cnx + v0->cy*cny + v0->cz*cnz ;
  if (dot < 0)
  { nx *= -1 ; ny *= -1 ; nz *= -1 ; }
  
  adx = v1->x - v0->x ; ady = v1->y - v0->y ; adz = v1->z - v0->z ;
  alen = sqrt(SQR(adx)+SQR(ady)+SQR(adz)) ;
  cdx = v2->x - v0->x ; cdy = v2->y - v0->y ; cdz = v2->z - v0->z ;
  clen = sqrt(SQR(cdx)+SQR(cdy)+SQR(cdz)) ;
  
  /*
    sample along legs of the triangle making sure the maximum spacing
    between samples (along the longer leg) is SAMPLE_DIST.
  */
  
  /*
    move along v0->v1 and v3->v2 lines and draw in crossing line to fill face
    t0 parameterizes lines from v0->v1 and v0->v2 
  */
  if (FZERO(alen) && FZERO(clen))
    delta_t0 = 0.99 ;
  else
    delta_t0 = (alen > clen) ? (SAMPLE_DIST / alen) : (SAMPLE_DIST / clen ) ;
  if (FZERO(delta_t0))
    ErrorReturn(ERROR_BADPARM, 
		(ERROR_BADPARM,
		 "mrisDefectFaceMRILogLikelihood: face has infinite leg (%d, %d)\n",
		 alen, clen)) ;
	
  if (delta_t0 >= 1.0)
    delta_t0 = 0.99 ;
  
  /* delta_t0 is % of alen or clen (whichever is bigger) of SAMPLE_DIST */
  for (nsamples = 0, ll = 0.0f, t0 = 0 ; t0 <= 1.0f ; t0 += delta_t0)
  {
    /* compute points (xa,ya,za) and (xc,yc,zc) on the a and c lines resp. */
    xa = v0->x + t0*adx ; ya = v0->y + t0*ady ; za = v0->z + t0*adz ;
    xc = v0->x + t0*cdx ; yc = v0->y + t0*cdy ; zc = v0->z + t0*cdz ;
    dx = xc-xa ; dy = yc-ya ; dz = zc-za ;
    len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
    if (FZERO(len))
      delta_t1 = 0.99 ;
    else
    {
      delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
      if (delta_t1 >= 1.0f)
	delta_t1 = 0.99 ;
    }
    
    /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
    for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1, nsamples++)
    {
      /* compute a point on the line connecting a and c */
      x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
      // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
      MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
      MRIsampleVolume(mri, xv, yv, zv, &gray_val) ;
      ll += log(h_gray->counts[nint(gray_val)]) ;
      // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
      MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
      MRIsampleVolume(mri, xv, yv, zv, &white_val) ;
      ll += log(h_white->counts[nint(white_val)]) ;
      grad = white_val - gray_val ;
      bin = nint((grad-h_grad->bins[0]) / h_grad->bin_size) ;
      if (bin < 0)
	bin = 0 ;
      else if (bin >= h_grad->nbins)
	bin = h_grad->nbins-1 ;
      ll += log(h_grad->counts[bin]) ;
      MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val) ;
      jll += log(val) ;
    }
    /* compute last point on line */
    t1 = 1.0f ;
    x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
    // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &gray_val) ;
    ll += log(h_gray->counts[nint(gray_val)]) ;
    // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &white_val) ;
    ll += log(h_white->counts[nint(white_val)]) ;
    grad = white_val - gray_val ;
    bin = nint((grad-h_grad->bins[0]) / h_grad->bin_size) ;
    if (bin < 0)
      bin = 0 ;
    else if (bin >= h_grad->nbins)
      bin = h_grad->nbins-1 ;
    ll += log(h_grad->counts[bin]) ;
    nsamples++ ;
    MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val) ;
    jll += log(val) ;
  }
  
  /* compute last line on the a and c lines resp. */
  t0 = 1.0f ;
  xa = v0->x + t0*adx ; ya = v0->y + t0*ady ; za = v0->z + t0*adz ;
  xc = v0->x + t0*cdx ; yc = v0->y + t0*cdy ; zc = v0->z + t0*cdz ;
  dx = xc-xa ; dy = yc-ya ; dz = zc-za ;
  len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
  if (FZERO(len))
    delta_t1 = 0.99 ;
  else
  {
    delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
    if (delta_t1 >= 1.0f)
      delta_t1 = 0.99 ;
  }
  
  /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
  for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1, nsamples++)
  {
    /* compute a point on the line connecting a and c */
    x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
    // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &gray_val) ;
    ll += log(h_gray->counts[nint(gray_val)]) ;
    // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &white_val) ;
    ll += log(h_white->counts[nint(white_val)]) ;
    grad = white_val - gray_val ;
    bin = nint((grad-h_grad->bins[0]) / h_grad->bin_size) ;
    if (bin < 0)
      bin = 0 ;
    else if (bin >= h_grad->nbins)
      bin = h_grad->nbins-1 ;
    ll += log(h_grad->counts[bin]) ;
    MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val) ;
    jll += log(val) ;
  }
  /* compute last point on line */
  t1 = 1.0f ;
  x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
  // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
  MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
  MRIsampleVolume(mri, xv, yv, zv, &gray_val) ;
  ll += log(h_gray->counts[nint(gray_val)]) ;
  // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
  MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
  MRIsampleVolume(mri, xv, yv, zv, &white_val) ;
  ll += log(h_white->counts[nint(white_val)]) ;
  grad = white_val - gray_val ;
  bin = nint((grad-h_grad->bins[0]) / h_grad->bin_size) ;
  if (bin < 0)
    bin = 0 ;
  else if (bin >= h_grad->nbins)
    bin = h_grad->nbins-1 ;
  ll += log(h_grad->counts[bin]) ;
  MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val) ;
  jll += log(val) ;
  nsamples++ ;

#if 0
  return(ll/nsamples) ;
#else
  return(jll/nsamples) ;
#endif
}

static double 
mrisComputeDefectMRILogLikelihood(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, 
				  int *vertex_trans, HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_grad, MRI *mri_gray_white)
{
  double sse = 0.0, nx, ny, nz, n1x, n1y, n1z, n2x, n2y, n2z, dx, dy, dz, d, len, total, ll, total_ll, jll, total_jll ;
  int    i, vno, n, nvertices = 0, nsamples, m, vno2, bin ;
  VERTEX *v, *vn ;
  Real   white_val, gray_val, val0, wval, gval, xv, yv, zv, x, y, z, vtotal, val ;
  double  origin[3], e0[3], e1[3], grad ;
  EDGE    edge ;

  for (total_jll = total_ll = 0.0, i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
      vno = vertex_trans[defect->vertices[i]] ;
    else
      vno = vertex_trans[defect->border[i-defect->nvertices]] ;
    if (vno < 0)
      continue ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    nvertices++ ;
    mrisComputeDefectVertexNormal(mris, vno, &n1x, &n1y, &n1z) ;
    vtotal = 0.0f ;
    edge.vno1 = vno ;
    for (nsamples = n = 0 ; n < v->vnum ; n++)
    {
      if (((vno == Gdiag_no) && (v->v[n] == Gx)) ||
	  ((vno == Gx) && (v->v[n] == Gdiag_no)))
	DiagBreak() ;
      vn = &mris->vertices[v->v[n]] ;
      mrisComputeDefectVertexNormal(mris, v->v[n], &n2x, &n2y, &n2z) ;
      nx = (n1x + n2x) / 2 ; ny = (n1y + n2y) / 2 ; nz = (n1z + n2z) / 2 ;
      len = sqrt(nx*nx + ny*ny + nz*nz) ;
      if (FZERO(len))
	len = 1 ;
      nx /= len ; ny /= len ; nz /= len ;
      
      /* sample MR values along line and build estimate of log likelihood
         as distance.
      */
      val0 = (v->val + vn->val) / 2 ;  /* gray/white border value */
      wval = (v->val2 + vn->val2) / 2 ;  /* white matter mean */
      gval = (v->val2bak + vn->val2bak) / 2 ;  /* gray matter mean */
      
      /* traverse the edge connecting the two vertices in equal increments, sampling
	 MRI volume outside and inside */
      dx = vn->origx - v->origx ; dy = vn->origy - v->origy ; dz = vn->origz - v->origz ;
      len = sqrt(dx*dx + dy*dy + dz*dz) ;
#define SAMPLE_DIST 0.25
      len = SAMPLE_DIST / len ;  /* sample every 1/2 mm */
      if (len > 1)
	len = 1 ;
      
      jll = ll = total = 0.0 ;
      for (d = 0 ; d <= 1 ; d += len, nsamples++)
      {
        x = v->origx+d*dx ; y = v->origy+d*dy ; z = v->origz+d*dz ;
	// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
	MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &gray_val) ;
	/*				ll += log(h_gray->counts[nint(gray_val)]) ;*/
        total += fabs(gray_val-gval) ;
	// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
	MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &white_val) ;
        total += fabs(white_val-wval) ;
	/*				ll += log(h_white->counts[nint(white_val)]) ;*/
	grad = white_val - gray_val ;
	bin = nint((grad-h_grad->bins[0]) / h_grad->bin_size) ;
	if (bin < 0)
	  bin = 0 ;
	else if (bin >= h_grad->nbins)
					bin = h_grad->nbins-1 ;
	/*				ll += log(h_grad->counts[bin]) ;*/
	if (!finite(ll))
	  DiagBreak() ;
	MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val) ;
	jll += log(val) ;
      }

      /* add one sample from the centroid of any face that will get added */
      edge.vno2 = vno2 = v->v[n] ;
      mrisComputeCanonicalEdgeBasis(mris, &edge, &edge, origin, e0, e1) ;
      for (m = 0 ; m < v->vnum ; m++)
      {
	if (v->v[m] == vno2)
	  continue ;
	if (vertexNeighbor(mris, vno2, v->v[m]) && 
	    !isFace(mris,vno, vno2, v->v[m]) &&
	    !containsAnotherVertex(mris,vno,vno2,v->v[m],e0,e1,origin))
	{
#if 1
	  jll += mrisDefectFaceMRILogLikelihood(mris, vno, vno2, v->v[m],
						mri, defect, 
						h_white, h_gray, h_grad, mri_gray_white) ;
#else
	  VERTEX *v0, *v1, *v2 ;
	  float  x, y, z ;
	  
	  v0 = v ; v1 = &mris->vertices[vno2] ; v2 = &mris->vertices[v->v[m]] ;
	  x = (v1->origx + v2->origx) / 2.0f ;
	  y = (v1->origy + v2->origy) / 2.0f ;
	  z = (v1->origz + v2->origz) / 2.0f ;
	  
	  /* now bisect v0->bisector line */
	  x = (v0->origx + x) / 2.0f ;
	  y = (v0->origy + y) / 2.0f ;
	  z = (v0->origz + z) / 2.0f ;
	  nsamples++ ;
	  
	  // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
	  MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
	  MRIsampleVolume(mri, xv, yv, zv, &gray_val) ;
	  ll += log(h_gray->counts[nint(gray_val)]) ;
	  total += fabs(gray_val-gval) ;
	  // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
	  MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
	  MRIsampleVolume(mri, xv, yv, zv, &white_val) ;
	  ll += log(h_white->counts[nint(white_val)]) ;
	  total += fabs(white_val-wval) ;
	  grad = white_val - gray_val ;
	  bin = nint((grad-h_grad->bins[0]) / h_grad->bin_size) ;
	  if (bin < 0)
	    bin = 0 ;
	  else if (bin >= h_grad->nbins)
	    bin = h_grad->nbins-1 ;
	  ll += log(h_grad->counts[bin]) ;
#endif
	}
      }

      if (!finite(total))
	DiagBreak() ;
      vtotal += total / (float)nsamples ;  /* 11 samples */
      total_ll += ll / (float)nsamples ;
      total_jll += jll / (float)nsamples ;
    }
    if (v->vnum == 0)
      continue ;
    sse += (vtotal / (float)v->vnum) ;
  }
#if 0
  return(total_ll/*/(double)nvertices*/) ;
#else
  return(total_jll/*/(double)nvertices*/) ;
#endif
}
#if 0
static double
mrisComputeDefectMRIEnergy(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, int *vertex_trans)
{
  double sse = 0.0, nx, ny, nz, n1x, n1y, n1z, n2x, n2y, n2z, dx, dy, dz, d, len, total ;
  int    i, vno, n, nvertices = 0, nsamples, m, vno2 ;
  VERTEX *v, *vn ;
  Real   val, val0, wval, gval, xv, yv, zv, x, y, z, vtotal ;
  double  origin[3], e0[3], e1[3] ;
  EDGE    edge ;

  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
      vno = vertex_trans[defect->vertices[i]] ;
    else
      vno = vertex_trans[defect->border[i-defect->nvertices]] ;
    if (vno < 0)
      continue ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    nvertices++ ;
    mrisComputeDefectVertexNormal(mris, vno, &n1x, &n1y, &n1z) ;
    vtotal = 0.0f ;
    edge.vno1 = vno ;
    for (nsamples = n = 0 ; n < v->vnum ; n++)
    {
      if (((vno == Gdiag_no) && (v->v[n] == Gx)) ||
	  ((vno == Gx) && (v->v[n] == Gdiag_no)))
	DiagBreak() ;
      vn = &mris->vertices[v->v[n]] ;
      mrisComputeDefectVertexNormal(mris, v->v[n], &n2x, &n2y, &n2z) ;
      nx = (n1x + n2x) / 2 ; ny = (n1y + n2y) / 2 ; nz = (n1z + n2z) / 2 ;
      len = sqrt(nx*nx + ny*ny + nz*nz) ;
      if (FZERO(len))
	len = 1 ;
      nx /= len ; ny /= len ; nz /= len ;
      
      /* sample MR values along line and build estimate of log likelihood
         as distance.
      */
      val0 = (v->val + vn->val) / 2 ;  /* gray/white border value */
      wval = (v->val2 + vn->val2) / 2 ;  /* white matter mean */
      gval = (v->val2bak + vn->val2bak) / 2 ;  /* gray matter mean */

      /* traverse the edge connecting the two vertices in equal increments, sampling
	 MRI volume outside and inside */
      dx = vn->origx - v->origx ; dy = vn->origy - v->origy ; dz = vn->origz - v->origz ;
      len = sqrt(dx*dx + dy*dy + dz*dz) ;
      len = .5 / len ;  /* sample every 1/2 mm */

      total = 0.0 ;
      for (d = 0 ; d <= 1 ; d += len, nsamples++)
      {
        x = v->origx+d*dx ; y = v->origy+d*dy ; z = v->origz+d*dz ;
	// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
	MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        total += fabs(val-gval) ;
	// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
	MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        total += fabs(val-wval) ;
      }

      /* add one sample from the centroid of any face that will get added */
      edge.vno2 = vno2 = v->v[n] ;
      mrisComputeCanonicalEdgeBasis(mris, &edge, &edge, origin, e0, e1) ;
      for (m = 0 ; m < v->vnum ; m++)
      {
	if (v->v[m] == vno2)
	  continue ;
	if (vertexNeighbor(mris, vno2, v->v[m]) && 
	    !isFace(mris,vno, vno2, v->v[m]) &&
	    !containsAnotherVertex(mris,vno,vno2,v->v[m],e0,e1,origin))
	{
	  VERTEX *v0, *v1, *v2 ;
	  float  x, y, z ;
						
	  v0 = v ; v1 = &mris->vertices[vno2] ; v2 = &mris->vertices[v->v[m]] ;
	  x = (v1->origx + v2->origx) / 2.0f ;
	  y = (v1->origy + v2->origy) / 2.0f ;
	  z = (v1->origz + v2->origz) / 2.0f ;
	  
	  /* now bisect v0->bisector line */
	  x = (v0->origx + x) / 2.0f ;
	  y = (v0->origy + y) / 2.0f ;
	  z = (v0->origz + z) / 2.0f ;
	  nsamples++ ;
	  
	  // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
	  MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
	  MRIsampleVolume(mri, xv, yv, zv, &val) ;
	  total += fabs(val-gval) ;
	  // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
	  MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
	  MRIsampleVolume(mri, xv, yv, zv, &val) ;
	  total += fabs(val-wval) ;
	}
      }
      
      if (!finite(total))
	DiagBreak() ;
      vtotal += total / (float)nsamples ;  /* 11 samples */
    }
    if (v->vnum == 0)
      continue ;
    sse += (vtotal / (float)v->vnum) ;
  }
  return(sse) ;
}
#endif

#include "tritri.h"
static int
mrisComputeDefectVertexNormal(MRI_SURFACE *mris, int vno, double *pnx, double *pny, double *pnz)
{
  VERTEX  *v, *vn ;
  int     n ;
  double  a0[3], an[3], ac0[3], acn[3], a_cross[3], ac_cross[3], dot, xc[3], nx, ny, nz, 
    ax, ay, az, len ;
  
  v = &mris->vertices[vno] ;
  if (v->vnum < 2)
  {
    *pnx = *pny = *pnz = 0 ;
    return(ERROR_BADPARM) ;
  }

  /* in original space */
  vn = &mris->vertices[v->v[0]] ;
  ax = vn->origx - v->origx ; ay = vn->origy - v->origy ; az = vn->origz - v->origz ;
  a0[0] = ax ; a0[1] = ay ; a0[2] = az ;

  /* on sphere */
  ax = vn->cx - v->cx ; ay = vn->cy - v->cy ; az = vn->cz - v->cz ;
  ac0[0] = ax ; ac0[1] = ay ; ac0[2] = az ;
  
  xc[0] = v->cx ; xc[1] = v->cy ; xc[2] = v->cz ;   /* vector for orienting cross products */
  
  nx = ny = nz = 0.0 ;
  for (n = 1 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    
    /* do original surface */
    ax = vn->origx - v->origx ; ay = vn->origy - v->origy ; az = vn->origz - v->origz ;
    an[0] = ax ; an[1] = ay ; an[2] = az ;
    CROSS(a_cross, a0, an) ;

    /* do it on sphere, to orient things */
    ax = vn->cx - v->cx ; ay = vn->cy - v->cy ; az = vn->cz - v->cz ;
    acn[0] = ax ; acn[1] = ay ; acn[2] = az ;
    CROSS(ac_cross, ac0, acn) ;
    dot = DOT(ac_cross, xc) ;
    if (dot < 0)
    { 
      a_cross[0] *= -1 ; a_cross[1] *= -1 ; a_cross[2] *= -1 ; 
    }
    nx += a_cross[0] ; ny += a_cross[1] ; nz += a_cross[2] ; 
  }
  
  len = sqrt(nx*nx + ny*ny + nz*nz) ;
  if (FZERO(len))
    len = 1 ;
  nx /= len ; ny /= len ; nz /= len ;
  *pnx = nx ; *pny = ny ; *pnz = nz ;
  return(NO_ERROR) ;
}

static double  l_mri   = 1.0 ;
static double  l_unmri = 0.0 ;
static double  l_curv  = 1.0 ;
static double  l_qcurv = 0.0 ;

static double
mrisComputeDefectLogLikelihood(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, int *vertex_trans, DEFECT_PATCH *dp,
			       HISTOGRAM *h_k1, HISTOGRAM *h_k2, HISTOGRAM *h_white, HISTOGRAM *h_gray,
			       HISTOGRAM *h_border, HISTOGRAM *h_grad,MRI *mri_gray_white,
			       HISTOGRAM *h_dot, TOPOLOGY_PARMS *parms)
{
  static int first_time = 1 ;
  double ll = 0.0 ;
  
  if (first_time)
  {
    char *cp ;
    l_mri = parms->l_mri ; 
    l_unmri = parms->l_unmri ; 
    l_curv = parms->l_curv ; 
    if ((cp = getenv("QCURV")) != NULL)
    {
      l_qcurv = atof(cp) ;
      printf("setting qcurv = %2.3f\n", l_qcurv) ;
    }
    if ((cp = getenv("CURV")) != NULL)
    {
      l_curv = atof(cp) ;
      printf("setting curv = %2.3f\n", l_curv) ;
    }
    if ((cp = getenv("MRI")) != NULL)
    {
      l_mri = atof(cp) ;
      printf("setting mri = %2.3f\n", l_mri) ;
    }
    if ((cp = getenv("UNMRI")) != NULL)
    {
      l_unmri = atof(cp) ;
      printf("setting unmri = %2.3f\n", l_unmri) ;
    }
    first_time = 0 ;
    if (!FZERO(l_mri))
      printf("l_mri = %2.2f ", l_mri) ;
    if (!FZERO(l_unmri))
      printf("l_unmri = %2.2f ", l_unmri) ;
    if (!FZERO(l_curv))
      printf("l_curv = %2.2f ", l_curv) ;
    printf("\n") ;
  }
  if (!FZERO(l_mri))
    ll += l_mri * mrisComputeDefectMRILogLikelihood(mris, mri, defect, vertex_trans, h_white, h_gray,h_grad, mri_gray_white) ;
  if (!FZERO(l_unmri))
    ll += l_unmri * mrisComputeDefectMRILogUnlikelihood(mris, mri, dp, vertex_trans, h_border) ;
  if (!FZERO(l_qcurv))
    ll += l_qcurv * mrisComputeDefectCurvatureLogLikelihood(mris, defect, vertex_trans, h_k1, h_k2) ;
  if (!FZERO(l_curv))
    ll += l_curv * mrisComputeDefectNormalDotLogLikelihood(mris, defect, vertex_trans, h_dot) ;
  return(ll) ;
}
#if 0
static double
mrisComputeDefectCurvatureEnergy(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans)
{
  double sse = 0.0, nx, ny, nz, dx, dy, dz, vtotal, dot ;
  int    i, vno, n, nvertices = 0 ;
  VERTEX *v, *vn ;
  
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
      vno = vertex_trans[defect->vertices[i]] ;
    else
      vno = vertex_trans[defect->border[i-defect->nvertices]] ;
    if (vno < 0)
      continue ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nvertices++ ;
    mrisComputeDefectVertexNormal(mris, vno, &nx, &ny, &nz) ;
    vtotal = 0.0f ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dx = vn->origx - v->origx ; dy = vn->origy - v->origy ; dz = vn->origz - v->origz ; 
      dot = dx*nx + dy*ny + dz*nz ;
      vtotal += (dot*dot) ;
      if (!finite(vtotal))
	DiagBreak() ;
    }
    if (v->vnum == 0)
      continue ;
    sse += (vtotal / (float)v->vnum) ;
  }
  return(sse) ;
}
#endif
static int
mrisFindSecondNeighborhood(MRI_SURFACE *mris, int vno, int *nbrs, int *num_nbrs)
{
  int     n, n2 ;
  VERTEX  *v, *vn, *vn2 ;
  
  *num_nbrs = 0 ;
  v = &mris->vertices[vno] ;  v->marked = 1 ;
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ; vn->marked = 1 ;
    nbrs[*num_nbrs] = v->v[n] ; *num_nbrs += 1 ;
  }

  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    for (n2 = 0 ; n2 < vn->vnum ; n2++)
    {
      vn2 = &mris->vertices[vn->v[n2]] ;
      if (vn2->marked)
	continue ;
      vn2->marked = 1 ;
      nbrs[*num_nbrs] = vn->v[n2] ; *num_nbrs += 1 ;
    }
  }
  
  v->marked = 0 ;
  for (n = 0 ; n < *num_nbrs ; n++)
    mris->vertices[nbrs[n]].marked = 0 ;
  return(NO_ERROR) ;
}
#if 0
static double
mrisComputeDefectQuadraticCurvatureEnergy(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans)
{
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n, i ;
  VERTEX   *v, *vn ;
  float    ui, vi, rsq, a, b ;
  double   sse = 0.0 ;
  int      nbrs[MAX_NBRS], num_nbrs ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  
  mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(2, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
      vno = vertex_trans[defect->vertices[i]] ;
    else
      vno = vertex_trans[defect->border[i-defect->nvertices]] ;
    if (vno < 0)
      continue ;
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->vnum <= 0)
      continue ;
    mrisFindSecondNeighborhood(mris, vno, nbrs, &num_nbrs) ;
    if (num_nbrs < 3)
      continue ;
    MRIScomputeSecondFundamentalFormAtVertex(mris, vno, nbrs, num_nbrs) ;
    v_Y = VectorAlloc(num_nbrs, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(num_nbrs, 2, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;
    for (n = 0 ; n < num_nbrs ; n++)  /* build data matrices */
    {
      vn = &mris->vertices[nbrs[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ; vi = V3_DOT(v_e2, v_nbr) ; 
      rsq = ui*ui + vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsq ;
      *MATRIX_RELT(m_R, n+1, 2) = 1 ;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL) ;
    if (!m_R_inv)
    {
      MatrixFree(&m_R) ; VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ;
		if (!finite(b))
			DiagBreak() ;
    sse += b*b ;
    if (vno == Gdiag_no)
      printf("v %d: curvature sse %2.2f\n", vno, b*b) ;
    MatrixFree(&m_R) ; VectorFree(&v_Y) ; MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ; VectorFree(&v_e1) ; VectorFree(&v_e2) ; 
  VectorFree(&v_nbr) ; VectorFree(&v_A) ;

  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  return(sse) ;
}
#endif


static int
mrisTessellateDefect(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, DEFECT *defect, int *vertex_trans, 
		     MRI *mri, HISTOGRAM *h_k1,HISTOGRAM *h_k2,HISTOGRAM *h_white, HISTOGRAM *h_gray,
		     HISTOGRAM *h_border, HISTOGRAM *h_grad, MRI *mri_gray_white,
		     HISTOGRAM *h_dot, TOPOLOGY_PARMS *parms)
{
#define MAX_DEFECT_VERTICES 200000
  int    i, j, vlist[MAX_DEFECT_VERTICES], n, nvertices, nedges, ndiscarded ;
  VERTEX *v, *v2 ;
  EDGE   *et ;
	/*  double  cx, cy, cz, max_len ;*/
  static int dno = 0 ;
  Real    x, y, z, xv, yv, zv, val0, val, total, dx, dy, dz, d, wval, gval, Ix, Iy, Iz ;
	float   norm1[3], norm2[3], nx, ny, nz ;


  /* first build table of all possible edges among vertices in the defect
     and on its border.
  */
  for (nvertices = i = 0 ; i < defect->nvertices ; i++)
  {
    if (nvertices >= MAX_DEFECT_VERTICES)
      ErrorExit(ERROR_NOMEMORY, "mrisTessellateDefect: too many vertices in defect (%d)",
								MAX_DEFECT_VERTICES) ;
    if (defect->vertices[i] == Gdiag_no)
      DiagBreak() ;
    if (vertex_trans[defect->vertices[i]] == Gdiag_no)
      DiagBreak() ;
    if (defect->status[i] == KEEP_VERTEX)
      vlist[nvertices++] = defect->vertices[i] ;
  }
  for (i = 0 ; i < defect->nborder ; i++)
  {
    if (defect->border[i] == Gdiag_no)
      DiagBreak() ;
    if (vertex_trans[defect->border[i]] == Gdiag_no)
      DiagBreak() ;
    vlist[nvertices++] = defect->border[i] ;
  }
  if (nvertices > 250)
    fprintf(stdout, 
            "retessellating defect %d with %d vertices (convex hull=%d).\n", 
            dno, nvertices, defect->nchull) ;
  dno++ ;
  if (nvertices == 0)  /* should never happen */
    return(NO_ERROR) ;

  nedges = (nvertices * (nvertices-1)) / 2 ;  /* won't be more than this */

  et = (EDGE *)calloc(nedges, sizeof(EDGE)) ;
  if (!et)
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate %d edges for retessellation", nvertices) ;

  for (n = i = 0 ; i < nvertices ; i++)
  {
    v = &mris->vertices[vlist[i]] ;
    if (vlist[i] == Gdiag_no)
      DiagBreak() ;
    if (vertex_trans[vlist[i]] == Gdiag_no)
      DiagBreak() ;
    for (j = i+1 ; j < nvertices ; j++, n++)
    {
      if (vlist[j] == Gdiag_no)
        DiagBreak() ;
      if (vlist[j] == Gdiag_no || vlist[i] == Gdiag_no)
        DiagBreak() ;
      if (vertex_trans[vlist[j]] == Gdiag_no || vertex_trans[vlist[i]] == Gdiag_no)
        DiagBreak() ;
      mrisComputeOrigNormal(mris, vlist[i], norm1) ; 
      mrisComputeOrigNormal(mris, vlist[j], norm2) ;
      nx = (norm1[0] + norm2[0]) / 2 ; ny = (norm1[1] + norm2[1]) / 2 ;
      nz = (norm1[2] + norm2[2]) / 2 ;
      total = sqrt(nx*nx + ny*ny + nz*nz) ;
      if (FZERO(total))
	total = 1 ;
      nx /= total ; ny /= total ; nz /= total ;
      v2 = &mris->vertices[vlist[j]] ;
      x = (v->origx+v2->origx)/2 ; y = (v->origy+v2->origy)/2 ; 
			z = (v->origz+v2->origz)/2 ;
      // MRIworldToVoxel(mri, x, y, z, &xv, &yv, &zv) ;
      MRIsurfaceRASToVoxel(mri, x, y, z, &xv, &yv, &zv) ;
      MRIsampleVolumeGradient(mri, xv, yv, zv, &Ix, &Iy, &Iz) ;
      total = sqrt(Ix*Ix + Iy*Iy + Iz*Iz) ;
      if (FZERO(total))
	total = 1 ;
      Ix /= total ; Iy /= total ; Iz /= total ;
      et[n].vno1 = vertex_trans[vlist[i]]; et[n].vno2 = vertex_trans[vlist[j]];
      if ((et[n].vno1 == 141823 && et[n].vno2 == 141908) ||
          (et[n].vno2 == 141823 && et[n].vno1 == 141908))
        DiagBreak() ;
      if ((vlist[i] == Gdiag_no && vlist[j] == Gx) ||
          (vlist[j] == Gdiag_no && vlist[i] == Gx))
        DiagBreak() ;
      if ((vertex_trans[vlist[i]] == Gdiag_no && vertex_trans[vlist[j]] == Gx) ||
          (vertex_trans[vlist[j]] == Gdiag_no && vertex_trans[vlist[i]] == Gx))
        DiagBreak() ;
      
      /* sample MR values along line and build estimate of log likelihood
         as distance.
      */
      val0 = (v->val + v2->val) / 2 ;  /* gray/white border value */
      wval = (v->val2 + v2->val2) / 2 ;  /* white matter mean */
      gval = (v->val2bak + v2->val2bak) / 2 ;  /* gray matter mean */

      /* sample one end point */
      x = v->origx ; y = v->origy ; z = v->origz ;
      // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
      MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
      MRIsampleVolume(mri, xv, yv, zv, &val) ;
      total = fabs(val-gval) ;
      // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
      MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
      MRIsampleVolume(mri, xv, yv, zv, &val) ;
      total += fabs(val-wval) ;

      /* sample the other end point */
      x = v2->origx ; y = v2->origy ; z = v2->origz ;
      // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
      MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
      MRIsampleVolume(mri, xv, yv, zv, &val) ;
      total += fabs(val-gval) ;
      // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
      MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
      MRIsampleVolume(mri, xv, yv, zv, &val) ;
      total += fabs(val-wval) ;

      dx = v2->origx - v->origx ;
      dy = v2->origy - v->origy ;
      dz = v2->origz - v->origz ;
      for (d = .1 ; d <= .9 ; d += .1)
      {
        /* sample the midpoint end point */
        x = v->origx+d*dx ;
        y = v->origy+d*dy ;
        z = v->origz+d*dz ;
	// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
	MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        total += fabs(val-gval) ;
	// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
	MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        total += fabs(val-wval) ;
      }

      et[n].len = total / (4.0+18.0) ;
      if (et[n].vno1 == 120811 && et[n].vno2 == 120951)
      {
        VERTEX *v1, *v2 ;
        v1 = &mris_corrected->vertices[et[n].vno1] ;
        v2 = &mris_corrected->vertices[et[n].vno2] ;
        fprintf(stdout, "v %d (%d) --> %d (%d), len = %2.3f\n",
                et[n].vno1, vlist[i], et[n].vno2, vlist[j], et[n].len) ;
        fprintf(stdout, "INFLATED:  (%2.1f, %2.1f, %2.1f) --> "
                "(%2.1f, %2.1f, %2.1f), len = %2.2f\n",
                v1->tx, v1->ty, v1->tz,
                v2->tx, v2->ty, v2->tz,
                sqrt(SQR(v1->tx-v2->tx)+SQR(v1->ty-v2->ty)+SQR(v1->tz-v2->tz)));
        fprintf(stdout, "CANON:  (%2.1f, %2.1f, %2.1f) --> "
                "(%2.1f, %2.1f, %2.1f), len = %2.2f\n",
                v1->cx, v1->cy, v1->cz,
                v2->cx, v2->cy, v2->cz,
                sqrt(SQR(v1->cx-v2->cx)+SQR(v1->cy-v2->cy)+SQR(v1->cz-v2->cz)));
        fprintf(stdout, "CURRENT:  (%2.1f, %2.1f, %2.1f) --> "
                "(%2.1f, %2.1f, %2.1f), len = %2.2f\n",
                v1->x, v1->y, v1->z,
                v2->x, v2->y, v2->z,
                sqrt(SQR(v1->x-v2->x)+SQR(v1->y-v2->y)+SQR(v1->z-v2->z)));
        fprintf(stdout, "ORIG:  (%2.1f, %2.1f, %2.1f) --> "
                "(%2.1f, %2.1f, %2.1f), len = %2.2f\n",
                v1->origx, v1->origy, v1->origz,
                v2->origx, v2->origy, v2->origz,
                sqrt(SQR(v1->origx-v2->origx)+SQR(v1->origy-v2->origy)+SQR(v1->origz-v2->origz)));
        DiagBreak() ;
      }
      if (edgeExists(mris_corrected, et[n].vno1, et[n].vno2))
      {
        et[n].used = USED_IN_ORIGINAL_TESSELLATION ;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "excluding existing edge %d <--> %d\n",
                  vlist[i],vlist[j]) ;
      }
#if 1
      if (!edgeExists(mris, vlist[i], vlist[j]))  /* prioritize edges in original tessellation */
	et[n].len += 100 ;
#endif
    }
  }

  /* find and discard all edges that intersect one that is already in the
     tessellation.
  */
  for (ndiscarded = i = 0 ; i < nedges ; i++)
  {
    if (et[i].used == 0)
      continue ;
    for (j = i+1 ; j < nedges ; j++)
    {
      if (et[j].used)
	continue ;
      if (edgesIntersect(mris_corrected, &et[i], &et[j]))
      {
	ndiscarded++ ;
	if (j < nedges-1)
	  memmove(&et[j], &et[j+1], (nedges-j-1)*sizeof(EDGE)) ;
	nedges-- ; j-- ;
      }
    }
  }

  printf("%d of %d overlapping edges discarded\n", ndiscarded, nedges) ;

  /* sort the edge list by edge length */
  qsort(et, nedges, sizeof(EDGE), compare_edge_length) ;
#if 0
  {
    char fname[STRLEN] ;
    FILE *fp ;

    sprintf(fname, "alledges%d.log", defect_no) ;
    fp = fopen(fname, "w") ;
    if (fp)
    {
      for (n = 0 ; n < nedges ; n++)
      {
        fprintf(fp, "et[%d]: v %d <--> %d, len %2.2f\n",
                n, et[n].vno1, et[n].vno2, et[n].len) ;
      }
      fclose(fp) ;
    }
  }
#endif

#if 0
  /* compute centroid of spherical points and use it to build planar basis */
  cx = cy = cz = 0.0f ;
  for (n = i = 0 ; i < defect->nvertices ; i++)
  {
    if (defect->status[i] == KEEP_VERTEX)
    {
      v = &mris->vertices[defect->vertices[i]] ;
      cx += v->cx ; cy += v->cy ; cz += v->cz ;
      n++ ;
    }
  }
  for (i = 0 ; i < defect->nborder ; i++)
  {
    v = &mris->vertices[defect->border[i]] ;
    cx += v->cx ; cy += v->cy ; cz += v->cz ;
    n++ ;
  }
#endif

  if (!n)   /* should never happen */
    return(NO_ERROR) ;

#if 0
  {
    VERTEX  *v ;
    FILE    *fp ;
    char    fname[STRLEN] ;
    float   cx, cy, cz, x, y ;
    int     vno ;

    sprintf(fname, "points%d.log", defect_no) ;
    fp  = fopen(fname, "w") ;
    for (i = 0 ; i < defect->nvertices ; i++)
    {
      if (defect->status[i] != KEEP_VERTEX)
        continue ;
      vno = vertex_trans[defect->vertices[i]] ; 
      if (vno < 0)
        DiagBreak() ;
      v = &mris_corrected->vertices[vno] ;
      cx = v->cx - origin[0] ; cy = v->cy - origin[1] ; cz = v->cz - origin[2];
      x = cx*e0[0] + cy*e0[1] + cz*e0[2] ; y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
      fprintf(fp, "%d %2.3f  %2.3f\n", vno, x, y);
    }
    for (i = 0 ; i < defect->nborder ; i++)
    {
      vno = vertex_trans[defect->border[i]] ; 
      if (vno < 0)
        DiagBreak() ;
      v = &mris_corrected->vertices[vno] ;
      cx = v->cx - origin[0] ; cy = v->cy - origin[1] ; cz = v->cz - origin[2];
      x = cx*e0[0] + cy*e0[1] + cz*e0[2] ; y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
      fprintf(fp, "%d %2.3f  %2.3f\n", vno, x, y);
    }
    fclose(fp) ;
  }
#endif

#if 0
  /* remove border edges which intersect in the existing tessellation */
  fprintf(stdout, "removing intersecting border edges...\n") ;
  for (n = 0 ; n < nedges ; n++)
  {
    if (et[n].used == NOT_USED)
      continue ;
    if (intersectDefectEdges(mris_corrected, defect, &et[n], vertex_trans))
    {
      et[n].used = NOT_USED ;
      fprintf(stdout,"removing edge %d (%d-->%d)\n",n,et[n].vno1,et[n].vno1);
    }
  }
  fprintf(stdout, "retessellating planar representation...\n") ;
#endif

  if (getenv("USE_GA_TOPOLOGY_CORRECTION") != NULL)
    mrisComputeOptimalRetessellation(mris, mris_corrected,mri, defect, vertex_trans,et, nedges,
				     h_k1,h_k2,h_white, h_gray, h_border, h_grad, mri_gray_white,
				     h_dot, parms) ;
  else
    mrisRetessellateDefect(mris, mris_corrected, defect, vertex_trans, et, nedges, NULL, NULL);
  
  free(et) ;
  defect_no++ ;     /* for diagnostics */
  return(NO_ERROR) ;
}

static int    mrisMutateDefectPatch(DEFECT_PATCH *dp, EDGE_TABLE *etable, double pmutation) ;
static int    mrisCrossoverDefectPatches(DEFECT_PATCH *dp1, DEFECT_PATCH *dp2,
																				 DEFECT_PATCH *dp_dst, EDGE_TABLE *etable) ;


#define NUM_TO_ADD_FROM_ONE_PARENT 1

static int
mrisCrossoverDefectPatches(DEFECT_PATCH *dp1, DEFECT_PATCH *dp2,
													 DEFECT_PATCH *dp_dst, EDGE_TABLE *etable)
{
  int          i1, i2, *added, i, isrc, j, nadded ;
  double       p ;
  DEFECT_PATCH *dp_src ;
  
  added = (int *)calloc(dp1->nedges, sizeof(int)) ;
  p = randomNumber(0.0, 1.0) ;
  if (p < 0.5)  /* add from first defect */
  {
    dp_src = dp1 ;
  }
  else
  {
    dp_src = dp2 ;
  }
  
  for (nadded = isrc = i1 = i2 = i = 0 ; i < dp_dst->nedges ; i++)
  {
    if (nadded >= NUM_TO_ADD_FROM_ONE_PARENT)
    {
      nadded = 0 ;
      
      if ((dp_src == dp1 && i2 < dp2->nedges) || i1 >= dp1->nedges)  /* use dp2 */
      {
	dp_src = dp2 ; isrc = i2 ;
      }
      else if (i1 < dp1->nedges)
      {
	dp_src = dp1 ; isrc = i1 ;
      }
    }
    else  /* keep adding from the same parent */
    {
      nadded++ ;
      if (dp_src == dp1)
	isrc = i1 ;
      else
	isrc = i2 ;
    }
    if (isrc >= dp_src->nedges) /* shouldn't happen */
    {
      i-- ;
      continue ;
    }
    
    /* find the next index in the src dp that hasn't been added yet */
    while (added[dp_src->ordering[isrc]])
      if (++isrc >= dp_src->nedges)
	break ;
    if (isrc >= dp_src->nedges)
    {
      i-- ;      /* process this one again */
      continue ;
    }
    if (dp_src->ordering[isrc] == Gdiag_no)
      DiagBreak() ;
    if (isrc == Gdiag_no)
      DiagBreak() ;
    
    if (dp_src == dp1)    /* update source index for next iteration */
      i1 = isrc+1 ;
    else
      i2 = isrc+1 ;
    
    dp_dst->ordering[i] = dp_src->ordering[isrc] ;
    if (dp_dst->ordering[i] >= dp_dst->nedges)
      DiagBreak() ;
    
    added[dp_src->ordering[isrc]] = 1 ; /* make sure every edge is represented */
  }
  
  for (i = 0 ; i < dp_dst->nedges ; i++)
  {
    if ((dp_dst->ordering[i] >= dp_dst->nedges) ||
	(dp_dst->ordering[i] < 0))
    {
      DiagBreak() ;
      dp_dst->ordering[i] = 0 ;
    }
    if (added[i] == 0)
    {
      for (j = 0 ; j < dp_dst->nedges ; j++)
	if (dp_dst->ordering[j] < 0)  /* nothing in this slot */
	{
	  dp_dst->ordering[j] = i ;
	  added[i] = 1 ; /* make sure every edge is represented */
	  break ;
	}
    }
  }

  free(added) ;
  
  return(NO_ERROR) ;
}
static int
mrisMutateDefectPatch(DEFECT_PATCH *dp, EDGE_TABLE *etable, double pmutation)
{
  int   i, j, eti, etj, tmp, *dp_indices ;
  double p ;

  dp_indices = (int *)calloc(dp->nedges, sizeof(int)) ;
  for (i = 0 ; i < dp->nedges ; i++)
    dp_indices[dp->ordering[i]] = i ;
  
  for (i = 0 ; i < dp->nedges ; i++)
  {
    p = randomNumber(0.0, 1.0) ;
    eti = dp->ordering[i] ;
    if (p < pmutation && etable->noverlap[eti] > 0)
    {
      if (etable->flags[eti] & ET_OVERLAP_LIST_INCOMPLETE)  /* swap any two */
      {
	j = (int)randomNumber(0.0, dp->nedges-.1) ;
	tmp = dp->ordering[i] ; dp->ordering[i] = dp->ordering[j] ; dp->ordering[j] = tmp ;
      }
      else   /* swap two edges that intersect */
      {
	j = (int)randomNumber(0.0, etable->noverlap[eti]-0.0001) ;
	etj = etable->overlapping_edges[eti][j] ;  /* index of jth overlapping edge */
	j = dp_indices[etj] ;  /* find where it is in this defect patch ordering */
	
	tmp = dp->ordering[i] ; dp->ordering[i] = dp->ordering[j] ; dp->ordering[j] = tmp ;
	dp_indices[dp->ordering[i]] = i ;
	dp_indices[dp->ordering[j]] = j ;
      }
    }
  }
  
  free(dp_indices) ;
  return(NO_ERROR) ;
}
static double
mrisDefectPatchFitness(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, MRI *mri,
		       DEFECT_PATCH *dp, int *vertex_trans, DEFECT_VERTEX_STATE *dvs,
		       HISTOGRAM *h_k1, HISTOGRAM *h_k2, HISTOGRAM *h_white, HISTOGRAM *h_gray,
		       HISTOGRAM *h_border, HISTOGRAM *h_grad, MRI *mri_gray_white,
		       HISTOGRAM *h_dot, TOPOLOGY_PARMS *parms)
{
  int i ;
  
  mrisRetessellateDefect(mris, mris_corrected, dp->defect, 
			 vertex_trans, dp->etable->edges, dp->nedges, dp->ordering,
			 dp->etable) ;
#if 0
  dp->fitness = -mrisComputeDefectEnergy(mris_corrected, mri, dp->defect, vertex_trans) ;
#else
  dp->fitness = mrisComputeDefectLogLikelihood(mris_corrected, mri, dp->defect, vertex_trans, dp,
					       h_k1, h_k2, h_white, h_gray, h_border, h_grad, mri_gray_white, h_dot, parms) ;
#endif
  mrisRestoreVertexState(mris_corrected, dvs) ; 
  
  /* reset the edges to the unused state (unless they were in the original tessellation */
  for (i = 0 ; i < dp->nedges ; i++)
    if (dp->etable->edges[i].used == USED_IN_NEW_TESSELLATION)
      dp->etable->edges[i].used = NOT_USED ;
  
  return(dp->fitness) ;
}

#define MAX_PATCHES         1000
#define SELECTION_PCT       0.50
#define ELITISM_PCT         0.10
#define MUTATION_PCT        0.1
static double MUTATION_PCT_INIT =  (MUTATION_PCT*1.0) ;
#define REPLACEMENT_PCT     0.1   /* replace this many with mutated versions of best */
#define MAX_UNCHANGED       3
static int max_unchanged = MAX_UNCHANGED ;


static int defectPatchRank(DEFECT_PATCH *dps, int index, int npatches) ;
static int mrisCopyDefectPatch(DEFECT_PATCH *dp_src, DEFECT_PATCH *dp_dst) ;

static int
mrisCopyDefectPatch(DEFECT_PATCH *dp_src, DEFECT_PATCH *dp_dst)
{
  int i ;
  
  dp_dst->etable = dp_src->etable ; dp_dst->defect = dp_src->defect ;
  dp_dst->nedges = dp_src->nedges ;
  dp_dst->fitness = dp_src->fitness ; dp_dst->rank = dp_src->rank ;
  
  for (i = 0 ; i < dp_src->nedges ; i++)
    dp_dst->ordering[i] = dp_src->ordering[i] ;
  
  return(NO_ERROR) ;
}


static int
defectPatchRank(DEFECT_PATCH *dps, int index, int npatches)
{
  int   i, rank = 0 ;

  for (i = 0 ; i < npatches ; i++)
  {
    if (i == index)
      continue ;
    
    if (FEQUAL(dps[i].fitness, dps[index].fitness))
    {
      if (index > i)
	rank++ ;
    }
    else if (dps[i].fitness > dps[index].fitness)
      rank++ ;
  }
  return(rank) ;
}

#define MAX_EDGES 1000

static int
mrisComputeOptimalRetessellation(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected,
				 MRI *mri, DEFECT *defect, int *vertex_trans,
				 EDGE *et, int nedges, HISTOGRAM *h_k1, HISTOGRAM *h_k2,
				 HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_border, 
				 HISTOGRAM *h_grad, MRI *mri_gray_white, HISTOGRAM *h_dot,
				 TOPOLOGY_PARMS *parms)
{
  DEFECT_VERTEX_STATE *dvs ;
  DEFECT_PATCH        dps1[MAX_PATCHES], dps2[MAX_PATCHES], *dps, *dp, *dps_next_generation ;
  int                 i, best_i, j, g, nselected, nreplacements,rank, nunchanged = 0,
    nelite, ncrossovers, k, l, noverlap ;
  int                 *overlap ;
  double              fitness, best_fitness, last_best, fitness_mean, fitness_sigma, fitness_norm, pfitness,
    two_sigma_sq ;
  static int dno = 0 ;   /* for debugging */
  EDGE_TABLE          etable ;
  int                 max_patches = MAX_PATCHES, ranks[MAX_PATCHES], next_gen_index, 
    selected[MAX_PATCHES], nzero, sno=0, max_edges, debug_patch_n=-1, nbest = 0 ;

  max_patches = parms->max_patches ; max_unchanged = parms->max_unchanged ;
  max_edges = MAX_EDGES ;
  
  if (dno == Gdiag_no)
    DiagBreak() ;
  
  if (getenv("DEBUG_PATCH") != NULL)
  {
    int debug_patch = atoi(getenv("DEBUG_PATCH"));
    if (debug_patch != dno)
      max_patches = 0 ;
    else
    {
      if (getenv("DEBUG_PATCH_N") != NULL)
      {
	debug_patch_n = atoi(getenv("DEBUG_PATCH_N")) ;
	printf("terminating after %dth best tessellation\n", debug_patch_n) ;
      }
    }
  }
  
  if (!max_patches)
  {
    dno++ ;  /* for debugging */
    mrisRetessellateDefect(mris, mris_corrected, defect, 
			   vertex_trans, et, nedges, NULL, NULL) ;
    return(NO_ERROR) ;
  }
  
#if 0
  if (dno != 5)
  {
    max_patches = 5 ;
    max_unchanged = 2 ;
    dno++ ;  /* for debugging */
    mrisRetessellateDefect(mris, mris_corrected, defect, 
			   vertex_trans, et, nedges, NULL, NULL) ;
    return(NO_ERROR) ;
  }
  else
  {
#if 0
    int bin ;
    
    for (bin = 60 ; bin <= 80 ; bin++)
      h_gray->counts[bin] = 0.000001 ;
#endif
    DiagBreak() ;
  }
#endif
  dno++ ;  /* for debugging */
  
  if (nedges > 200000)
  {
    mrisRetessellateDefect(mris, mris_corrected, defect, 
			   vertex_trans, et, nedges, NULL, NULL) ;
    return(NO_ERROR) ;
  }
  else if (nedges > 100000)
  {
    max_unchanged = MIN(max_unchanged, 1) ;
    max_patches = MIN(max_patches, 10) ;
    max_edges = MIN(max_edges, 100) ;
  }
  else if (nedges > 50000)
  {
    max_patches = max_patches/2 ;
    max_edges = max_edges/5 ;
    max_unchanged = max_unchanged / 5 ;
  }
  
  etable.nedges = nedges ;
  etable.edges = (EDGE *)calloc(nedges, sizeof(EDGE)) ;
  etable.overlapping_edges = (int **)calloc(nedges, sizeof(int *)) ;
  etable.noverlap = (int *)calloc(nedges, sizeof(int)) ;
  etable.flags = (unsigned char *)calloc(nedges, sizeof(unsigned char)) ;
  overlap = (int *)calloc(nedges, sizeof(int)) ;
  if (!etable.edges || !etable.overlapping_edges || !etable.noverlap || !overlap)
    ErrorExit(ERROR_NOMEMORY, "mrisComputeOptimalRetessellation: could not allocate %d edge table",
	      nedges) ;
  
  memmove(etable.edges, et, nedges*sizeof(EDGE)) ;
  
  for (nzero = i = 0 ; i < nedges ; i++)  /* compute overlapping for each edge */
  {
    if (nedges > 50000 && !(i % 25000))
      printf("%d of %d edges processed\n", i, nedges) ;
    etable.noverlap[i] = 0 ;
    for (noverlap = j = 0 ; j < nedges ; j++)
    {
      if (j == i)
	continue ;
      if (edgesIntersect(mris_corrected, &et[i], &et[j]))
      {
	overlap[noverlap] = j ;
	noverlap++ ;
      }
      if (noverlap > MAX_EDGES)
	break ;
    }
    if (noverlap > 0)
    {
      if (noverlap > MAX_EDGES)
      {
	etable.noverlap[i] = MAX_EDGES ;
	etable.flags[i] |= ET_OVERLAP_LIST_INCOMPLETE ;
      }
      else
	etable.noverlap[i] = noverlap ;
      
      etable.overlapping_edges[i] = (int *)calloc(etable.noverlap[i], sizeof(int)) ;
      if (!etable.overlapping_edges[i])
	ErrorExit(ERROR_NOMEMORY, 
		  "mrisComputeOptimalRetessellation: could not allocate overlap list %d "
		  "with %d elts",i, etable.noverlap[i]) ;
      memmove(etable.overlapping_edges[i], overlap, etable.noverlap[i]*sizeof(int)) ;
    }
    else
      nzero++ ;
  }
  
  free(overlap) ;
  
  dvs = mrisRecordVertexState(mris_corrected, defect, vertex_trans) ;
  dps = dps1 ;
  
  /* generate initial population of patches */
  best_fitness = -1000000 ; best_i = 0 ;
  for (i = 0 ; i < max_patches ; i++)
  {
    dp = &dps2[i] ;
    dp->nedges = nedges ; dp->defect = defect ; dp->etable = &etable ;
    dp->ordering = (int *)calloc(nedges, sizeof(int)) ;
    if (!dp->ordering)
      ErrorExit(ERROR_NOMEMORY, "could not allocate %dth defect patch with %d indices",
		i, nedges) ;
    for (j = 0 ; j < nedges ; j++)
      dp->ordering[j] = j ;   /* initial in same order - will change later */
    
    dp = &dps1[i] ;
    dp->nedges = nedges ; dp->defect = defect ; dp->etable = &etable ;
    dp->ordering = (int *)calloc(nedges, sizeof(int)) ;
    if (!dp->ordering)
      ErrorExit(ERROR_NOMEMORY, "could not allocate %dth defect patch with %d indices",
		i, nedges) ;
    for (j = 0 ; j < nedges ; j++)
      dp->ordering[j] = j ;   /* initial in same order - will change later */
    
    if (i)  /* first one is in same order as original edge table */
      mrisMutateDefectPatch(dp, &etable, MUTATION_PCT_INIT) ;
    fitness = mrisDefectPatchFitness(mris, mris_corrected, mri, dp, vertex_trans, dvs, 
				     h_k1,h_k2,h_white,h_gray,h_border,h_grad,mri_gray_white,
				     h_dot, parms) ;
    if (i == 0 && Gdiag & 0x1000000)
    {
      int i ;
      char fname[STRLEN] ;
      sprintf(fname, "%s_defect%d_%03d", mris->fname, dno-1, sno++) ;
      dp = &dps[best_i] ;
      mrisRetessellateDefect(mris, mris_corrected, dp->defect, 
			     vertex_trans, dp->etable->edges, dp->nedges, dp->ordering,
			     dp->etable) ;
      MRISsaveVertexPositions(mris_corrected, TMP_VERTICES) ; 
      MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
      printf("writing surface snapshow to %s...\n",fname) ;
      MRISwrite(mris_corrected, fname) ;
      MRISrestoreVertexPositions(mris_corrected, TMP_VERTICES) ;
      mrisRestoreVertexState(mris_corrected, dvs) ; 
      /* reset the edges to the unused state (unless they were in the original tessellation */
      for (i = 0 ; i < dp->nedges ; i++)
	if (dp->etable->edges[i].used == USED_IN_NEW_TESSELLATION)
	  dp->etable->edges[i].used = NOT_USED ;
    }
    
    if (!i)
    {
      printf("defect %d: initial fitness = %2.4e, nvertices=%d, nedges=%d, max patches=%d\n", dno-1, fitness,
	     defect->nvertices, nedges, max_patches) ;
      best_fitness = fitness ; best_i = 0 ;
      if (++nbest == debug_patch_n)
	goto debug_use_this_patch ;
    }
    
    if (fitness > best_fitness)
    {
      best_fitness = fitness ; best_i = i ;
      printf("new optimal fitness found at %d: %2.4e\n", i, fitness) ;
      if (++nbest == debug_patch_n)
	goto debug_use_this_patch ;
    }
  }
  
  nelite = nint(ELITISM_PCT*max_patches) ;  /* # to keep in next generation */
  if (nelite < 1)
    nelite = 1 ;
  nselected = nint(SELECTION_PCT*max_patches) ;  /* # to allow to crossover */
  if (nselected < 2)
    nselected = 2 ;
  nreplacements = nint(REPLACEMENT_PCT*max_patches) ; /* # to replace with mutated versions of best */
  if (nreplacements + nelite > max_patches)
    nreplacements = max_patches-nelite ;
  g = 0 ; dps = dps1 ;
  while (nunchanged < max_unchanged)
  {
    if (dps == dps1)
      dps_next_generation = dps2 ;
    else
      dps_next_generation = dps1 ;
    
    last_best = best_fitness ;
    for (i = 0 ; i < max_patches ; i++)
    {
      dp = &dps[i] ;
      dp->rank = rank = defectPatchRank(dps, i, max_patches) ;
      ranks[rank] = i ;
    }
    
    /* first add the 'elite' group that are retained unchanged */
    next_gen_index = 0 ;
    for (i = 0 ; i < nelite ; i++)
      mrisCopyDefectPatch(&dps[ranks[i]], &dps_next_generation[next_gen_index++]) ;
    
    /* now replace the worst ones with mutated copies of the best */
    for (i = 0 ; i < nreplacements ; i++)
    {
      dp = &dps_next_generation[next_gen_index++] ;
      mrisCopyDefectPatch(&dps[ranks[i]], dp) ;
      mrisMutateDefectPatch(dp, &etable, MUTATION_PCT) ;
      fitness = mrisDefectPatchFitness(mris, mris_corrected, mri, dp, vertex_trans, dvs,
				       h_k1,h_k2,h_white,h_gray, h_border,h_grad, mri_gray_white, h_dot, parms) ;
      if (fitness > best_fitness)
      {
	nunchanged = 0 ;
	best_fitness = fitness ; best_i = next_gen_index-1 ;
	if  (Gdiag & DIAG_SHOW)
	  printf("replacement %d MUTATION: new optimal fitness found at %d: %2.4e\n", 
		 i, best_i, fitness) ;
	nmut++ ;
	if (++nbest == debug_patch_n)
	{
	  dps = dps_next_generation ;
	  goto debug_use_this_patch ;
	}
      }
    }
    
    for (fitness_mean = fitness_sigma = 0.0, i = 0 ; i < max_patches ; i++)
    {
      dp = &dps[i] ;
      fitness_mean += dp->fitness ;
      fitness_sigma += dp->fitness*dp->fitness ;
    }
    
    fitness_mean /= (float)max_patches ;
    fitness_sigma = (fitness_sigma/max_patches - (fitness_mean*fitness_mean)) ;
    if (fitness_sigma < 0)
      fitness_sigma = 0 ;
    else
      fitness_sigma = sqrt(fitness_sigma) ;
    if (!finite(fitness_sigma))
      DiagBreak() ;
    
    two_sigma_sq = (dps[ranks[0]].fitness-dps[ranks[nselected-1]].fitness) ;
    if (FZERO(two_sigma_sq))
      two_sigma_sq = 1 ;
    for (fitness_norm = 0.0, j = 0 ; j < nselected ; j++)
    {
      i = ranks[j] ; dp = &dps[i] ;
      fitness_norm += exp((dp->fitness-dps[ranks[0]].fitness)/two_sigma_sq) ;  /* make them positive and increasing */
    }
    if (FZERO(fitness_norm))  /* something wrong */
    {
      for (i = 0 ; i < max_patches ; i++)
      {
	dp = &dps[i] ;
	dp->rank = rank = defectPatchRank(dps, i, max_patches) ;
	if (dp->fitness >= best_fitness)
	{
	  best_fitness = dp->fitness ;
	  best_i = i ;
	}
	ranks[rank] = i ;
      }
      break ;
    }
    
    ncrossovers = max_patches-(nelite+nreplacements) ;
    for (l = k = j = 0 ; j < nselected ; j++)
    {
      int  nadd ;
      
      i = ranks[j] ; dp = &dps[i] ;
      pfitness = exp((dp->fitness-dps[ranks[0]].fitness)/two_sigma_sq)/fitness_norm ;
      nadd = nint(pfitness * ncrossovers) ;
      if (nadd >= ncrossovers)
	nadd = ncrossovers-1 ;
      else if (nadd == 0)
	nadd = 1 ;
      
      for (k = 0 ; l < ncrossovers && k < nadd ; l++, k++)
	selected[l] = i ;
    }
    for ( ; l < ncrossovers ; l++)   /* fill out rest of list */
    {
      double p ;
      p = randomNumber(0.0, 1.0)  ;
      for (fitness = 0.0, j = 0 ; j < nselected ; j++)
      {
	i = ranks[j] ; dp = &dps[i] ;
	pfitness = exp(dp->fitness/two_sigma_sq)/fitness_norm ;
	fitness += pfitness ;
	if (fitness > p)
	  break ;
      }
      selected[l] = i ;
    }
    
    for (i = 0 ; i < ncrossovers ; i++)
    {
      int   p1, p2 ;
      
      p1 = selected[i] ;
      do   /* select second parent at random */
      {
	p2 = selected[(int)randomNumber(0, ncrossovers-.001)] ;
      } while (p2 == p1) ;
      
      dp = &dps_next_generation[next_gen_index++] ;
      mrisCrossoverDefectPatches(&dps[p1], &dps[p2], dp, &etable) ;
      fitness = mrisDefectPatchFitness(mris, mris_corrected, mri, dp, vertex_trans, dvs,
				       h_k1,h_k2,h_white,h_gray, h_border,h_grad, mri_gray_white, h_dot, parms) ;
      if (fitness > best_fitness)
      {
	nunchanged = 0 ;
	best_fitness = fitness ; best_i = next_gen_index-1 ;
	if  (Gdiag & DIAG_SHOW)
	  printf("CROSSOVER (%d x %d): new optimal fitness found at %d: %2.4e\n", 
		 dps[p1].rank, dps[p2].rank, best_i, fitness) ;
	ncross++ ;
	if (++nbest == debug_patch_n)
	{
	  dps = dps_next_generation ;
	  goto debug_use_this_patch ;
	}
      }
      else   /* mutate it also */
      {
	mrisMutateDefectPatch(dp, &etable, MUTATION_PCT) ;
	fitness = mrisDefectPatchFitness(mris, mris_corrected, mri, dp, vertex_trans, dvs,
					 h_k1,h_k2,h_white,h_gray, h_border, h_grad, mri_gray_white, h_dot, parms) ;
	if (fitness > best_fitness)
	{
	  nunchanged = 0 ;
	  best_fitness = fitness ; best_i = next_gen_index-1 ;
	  if (Gdiag & DIAG_SHOW)
	    printf("MUTATION: new optimal fitness found at %d: %2.4e\n", best_i, fitness) ;
	  if (++nbest == debug_patch_n)
	  {
	    dps = dps_next_generation ;
	    goto debug_use_this_patch ;
	  }
	  nmut++ ;
	}
      }
    }
    
    /* make next generation current */
    if (dps == dps1)
    {
      dps = dps2 ; dps_next_generation = dps1 ;
    }
    else
    {
      dps = dps1 ; dps_next_generation = dps2 ;
    }
    
    /* regenerate rankings */
#if 0
    for (fitness_mean = fitness_sigma = 0.0, i = 0 ; i < max_patches ; i++)
    {
      dp = &dps[i] ;
      dp->fitness = fitness = mrisDefectPatchFitness(mris, mris_corrected, mri, dp, vertex_trans, dvs,
						     h_k1,h_k2,h_white,h_gray, h_border,h_grad,mri_gray_white, h_dot, parms) ;
      if (!i)
      {
	best_fitness = dps[0].fitness ; best_i = 0 ;
      }
    }
#endif
    
    best_fitness = dps[0].fitness ; best_i = 0 ;
    for (fitness_mean = fitness_sigma = 0.0, i = 0 ; i < max_patches ; i++)
    {
      dp = &dps[i] ;
      dp->rank = rank = defectPatchRank(dps, i, max_patches) ;
      if (dp->fitness >= best_fitness)
      {
	best_fitness = dp->fitness ;
	best_i = i ;
      }
      ranks[rank] = i ;
      fitness_mean += dp->fitness ;
      fitness_sigma += dp->fitness*dp->fitness ;
    }
    
    fitness_mean /= (float)max_patches ;
    fitness_sigma = sqrt(fitness_sigma/max_patches - (fitness_mean*fitness_mean)) ;
    printf("generation %d complete, optimal fitness = %2.4e (%2.4e +- %2.4e)\n", ++g,best_fitness,
	   fitness_mean, fitness_sigma);
    if (FEQUAL(last_best, best_fitness))
      nunchanged++ ;
    else
      nunchanged = 0 ;
  }

#if 0
  dp = &dps[best_i] ;
  fitness = mrisDefectPatchFitness(mris, mris_corrected, mri, dp, vertex_trans,dvs,h_k1,h_k2,h_white,h_gray,
				   h_border, h_grad, mri_gray_white, h_dot, parms);
  best_fitness = dps[0].fitness ; best_i = 0 ;
  for (i = 0 ; i < max_patches ; i++)
  {
    dp = &dps[i] ;
    dp->fitness = fitness = mrisDefectPatchFitness(mris, mris_corrected, mri, dp, vertex_trans, dvs,
						   h_k1,h_k2,h_white,h_gray, h_border, h_grad, mri_gray_white, h_dot, parms) ;
    if (fitness > best_fitness)
    {
      best_fitness = fitness ; best_i = i ;
    }
  }
#endif
 debug_use_this_patch:
  dp = &dps[best_i] ;
  fitness = mrisDefectPatchFitness(mris, mris_corrected, mri, dp, vertex_trans,dvs,h_k1,h_k2,h_white,h_gray,
				   h_border, h_grad, mri_gray_white, h_dot, parms);
  if (fitness != best_fitness)
    fprintf(stderr, "Warning - incorrect dp selected!!!!\n") ;
  mrisRetessellateDefect(mris, mris_corrected, dp->defect, vertex_trans, dp->etable->edges, nedges, 
			 dp->ordering, dp->etable) ;
  
  /* free everything */
  mrisFreeDefectVertexState(dvs) ;
  for (i = 0 ; i < max_patches ; i++)
  {
    free(dps1[i].ordering) ;
    free(dps2[i].ordering) ;
  }
  
  for (i = 0 ; i < nedges ; i++)
  {
    if (etable.overlapping_edges[i])
      free(etable.overlapping_edges[i]) ;
  }
  free(etable.overlapping_edges) ;
  free(etable.edges) ;
  free(etable.noverlap) ;
  free(etable.flags) ;
  return(NO_ERROR) ;
}

static int
mrisRetessellateDefect(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, DEFECT *defect, 
		       int *vertex_trans, EDGE *et, int nedges, int *ordering, EDGE_TABLE *etable)
{
  double  max_len ;
  int     i, j, max_i, max_added, nadded, index ;
  int     (*intersection_function)(MRI_SURFACE *mris, DEFECT *defect, EDGE *e, int *vertex_trans) ;
  
  max_len = 0 ; max_i = 0 ; max_added = nadded = 0 ;
  intersection_function = intersectDefectEdges ;
  for (index = 0 ; index < nedges ; index++)
  {
    if (ordering)
      i = ordering[index] ;
    else
      i = index ;
    
    if (i == 1662 || i == 1685)
      DiagBreak() ;
    if ((et[i].vno1 == 52022 || et[i].vno2 == 52022) &&
        (et[i].vno1 == 53156 || et[i].vno2 == 53156))
      DiagBreak() ;
    if (et[i].vno1 == Gdiag_no || et[i].vno2 == Gdiag_no)
      DiagBreak() ;
    if (et[i].used)  /* already exists in tessellation - don't add it again */
      continue ;
    
    if (etable)   /* use pre-computed intersection table */
    {
      int intersects = 0 ;
      
      for (j = 0 ; j < etable->noverlap[i] ; j++)
	if (et[etable->overlapping_edges[i][j]].used)
	{
	  intersects = 1 ;
	  break ;
	}
      if (intersects)
	continue ;
      if (etable->flags[i] & ET_OVERLAP_LIST_INCOMPLETE)
	intersection_function = intersectDefectEdges ;
      else
	intersection_function = intersectDefectConvexHullEdges ;
    }
    
    if ((*intersection_function)(mris_corrected, defect, &et[i], vertex_trans) == 0)
    {
      mrisAddEdge(mris_corrected, et[i].vno1, et[i].vno2) ;
      et[i].used = USED_IN_NEW_TESSELLATION ; nadded++ ;
      if (et[i].len > max_len)
      {
	max_len = et[i].len ; max_added = nadded-1 ;
	max_i = i ;
      }
    }
  }

  return(NO_ERROR) ;
}

static int
compare_edge_length(const void *pe0, const void *pe1)
{
  register EDGE *e0, *e1 ;

  e0 = (EDGE *)pe0 ;
  e1 = (EDGE *)pe1 ;

/*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (e0->len > e1->len)
    return(1) ;
  else if (e0->len < e1->len)
    return(-1) ;

  return(0) ;
}
#if 0
static int
mrisCheckDefectEdges(MRI_SURFACE *mris, DEFECT *defect, int vno, 
                     int *vertex_trans)
{
  int    i, n1, n2, n3, vno2, fno1, fno2 ;
  EDGE   edge1, edge2 ;
  VERTEX *v, *vn ;
  FACE   *f ;

  v = &mris->vertices[vno] ;
  edge1.vno1 = vno ;
  for (n1 = 0 ; n1 < v->num ; n1++)
  {
    n2 = v->n[n1] == VERTICES_PER_FACE-1 ? 0 : v->n[n1]+1 ;
    fno1 = v->f[n1] ;
    edge1.vno2 = mris->faces[fno1].v[n2] ;
    for (i = 0 ; i < defect->nborder ; i++)
    {
      vno2 = vertex_trans[defect->border[i]] ;
      if (vno2 < 0)
        continue ;
      vn = &mris->vertices[vno2] ;
      for (n2 = 0 ; n2 < vn->num ; n2++)
      {
        fno2 = vn->f[n2] ;
        f = &mris->faces[fno2] ;
        if (fno2 == Gdiag_no)
          DiagBreak() ;
        if ((fno2 == 103815 && fno1 == 101608) ||
            (fno1 == 103815 && fno2 == 101608))
          DiagBreak() ;
        for (n3 = 0 ; n3 < VERTICES_PER_FACE ; n3++)
        {
          edge2.vno1 = f->v[n3] ; 
          edge2.vno2 = f->v[n3 == VERTICES_PER_FACE-1 ? 0 : n3+1] ;
          if (edge2.vno1 == edge1.vno1 || edge2.vno1 == edge1.vno2 ||
              edge2.vno2 == edge1.vno1 || edge2.vno2 == edge1.vno2)
            continue ;
          if (edgesIntersect(mris, &edge1, &edge2))
          {
            fprintf(stdout, "triangles %d (a=%f) and %d (a=%f) intersect!\n",
                    fno1, mris->faces[fno1].area, fno2,
                    mris->faces[fno2].area) ;
            mrisDumpTriangle(mris, fno1) ;
            mrisDumpTriangle(mris, fno2) ;
          }
        }
      }
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          See if the edge e intersects any edges in the defect or
          it's border. Sorry this code is such a hatchet job. I'm sure
          there are far more elegant ways of doing intersection
          (e.g. sorting).
------------------------------------------------------*/
static int
intersectDefectEdges(MRI_SURFACE *mris, DEFECT *defect, EDGE *e, 
                     int *vertex_trans)
{
  VERTEX *v ;
  int    i, n, vno ;
  FILE   *fp = NULL ;
  EDGE   edge2 ;

  if ((e->vno1 == 109259 && e->vno2 == 108332) ||
      (e->vno2 == 109259 && e->vno1 == 108332))
    DiagBreak() ;

  for (i = 0 ; i < defect->nvertices+defect->nchull ; i++)
  {
    if (i < defect->nvertices)
      vno = vertex_trans[defect->vertices[i]] ;
    else
      vno = vertex_trans[defect->chull[i-defect->nvertices]] ;

    if (vno == e->vno1 || vno == e->vno2 || vno < 0)
      continue ;

    edge2.vno1 = vno ;
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      if (vno == Gdiag_no)
        DiagBreak() ;

      if ((vno == 53153 && v->v[n] == 53158) ||
          (v->v[n] == 53153 && vno == 53158))
        DiagBreak() ;
      if ((vno == 52024 && v->v[n] == 53158) ||
          (v->v[n] == 52024 && vno == 53158))
        DiagBreak() ;
      if (v->v[n] == e->vno1 || v->v[n] == e->vno2)
        continue ;
      edge2.vno2 = v->v[n] ;
			if (edgesIntersect(mris, e, &edge2))
				return(1) ;
    }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    double     x1_start, x1_end, y1_start, y1_end, x2_start, x2_end, 
               y2_start,y2_end, x2min, x2max, y2min, y2max, cx, cy, cz, 
               origin[3], e0[3], e1[3] ;
    char       fname[STRLEN] ;
    VERTEX     *v2 ;
    static int dno = -1 ;

    mrisComputeCanonicalEdgeBasis(mris, e, e, origin, e0, e1) ;
    sprintf(fname, "lines%d.log", defect_no) ;
    v = &mris->vertices[e->vno1] ; v2 = &mris->vertices[e->vno2] ;
    cx = v->cx-origin[0] ; cy = v->cy - origin[1] ; cz = v->cz-origin[2];
    x1_start = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
    y1_start = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
    cx = v2->cx-origin[0]; cy = v2->cy - origin[1] ; cz = v2->cz-origin[2];
    x1_end = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
    y1_end = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
    if (dno != defect_no)  /* first time for this defect */
    {
      fp = fopen(fname, "w") ;

      for (i = 0 ; i < defect->nvertices+defect->nchull ; i++)
      {
        if (i < defect->nvertices)
          vno = vertex_trans[defect->vertices[i]] ;
        else
          vno = vertex_trans[defect->chull[i-defect->nvertices]] ;
        
        if (vno == e->vno1 || vno == e->vno2 || vno < 0)
          continue ;
        v = &mris->vertices[vno] ;
        cx = v->cx-origin[0] ; cy = v->cy-origin[1] ; cz = v->cz-origin[2];
        x2_start = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
        y2_start = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
        for (n = 0 ; n < v->vnum ; n++)
        {
          if (vno == 6167)
            DiagBreak() ;
          if (v->v[n] == e->vno1 || v->v[n] == e->vno2)
            continue ;
          v2 = &mris->vertices[v->v[n]] ;
          cx = v2->cx-origin[0] ; cy = v2->cy-origin[1]; cz = v2->cz-origin[2];
          x2_end = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
          y2_end = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
          x2min = MIN(x2_start, x2_end) ; x2max = MAX(x2_start, x2_end) ;
          y2min = MIN(y2_start, y2_end) ; y2max = MAX(y2_start, y2_end) ;
          fprintf(fp, "%2.3f  %2.3f\n%2.3f  %2.3f\n\n",
                  x2_start, y2_start, x2_end, y2_end) ;
        }
      }

    }
    else
      fp = fopen(fname, "a") ;
    fprintf(fp, "%2.3f  %2.3f\n%2.3f  %2.3f\n\n",
            x1_start, y1_start, x1_end, y1_end) ;
    if (fp)
      fclose(fp) ;


    sprintf(fname, "edges%d.log", defect_no) ;
    if (dno != defect_no)  /* first time for this defect */
    {
      dno = defect_no ;
      fp = fopen(fname, "w") ;

      for (i = 0 ; i < defect->nvertices+defect->nchull ; i++)
      {
        if (i < defect->nvertices)
          vno = vertex_trans[defect->vertices[i]] ;
        else
          vno = vertex_trans[defect->chull[i-defect->nvertices]] ;
        
        if (vno == e->vno1 || vno == e->vno2 || vno < 0)
          continue ;
        v = &mris->vertices[vno] ;
        cx = v->cx-origin[0] ; cy = v->cy-origin[1] ; cz = v->cz-origin[2];
        x2_start = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
        y2_start = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
        for (n = 0 ; n < v->vnum ; n++)
        {
          if (vno == 6167)
            DiagBreak() ;
          if (v->v[n] == e->vno1 || v->v[n] == e->vno2)
            continue ;
          v2 = &mris->vertices[v->v[n]] ;
          cx = v2->cx-origin[0] ; cy = v2->cy-origin[1]; cz = v2->cz-origin[2];
          x2_end = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
          y2_end = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
          x2min = MIN(x2_start, x2_end) ; x2max = MAX(x2_start, x2_end) ;
          y2min = MIN(y2_start, y2_end) ; y2max = MAX(y2_start, y2_end) ;
          fprintf(fp, "%d  %d  %2.3f  %2.3f  %2.3f  %2.3f\n", vno, v->v[n],
                  x2_start, y2_start, x2_end, y2_end) ;
        }
      }

    }
    else
      fp = fopen(fname, "a") ;
    fprintf(fp, "%d  %d  %2.3f  %2.3f  %2.3f  %2.3f\n", e->vno1, e->vno2,
            x1_start, y1_start, x1_end, y1_end) ;
    if (fp)
      fclose(fp) ;
  }
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          See if the edge e intersects any edges in the defect or
          it's border. Sorry this code is such a hatchet job. I'm sure
          there are far more elegant ways of doing intersection
          (e.g. sorting).
------------------------------------------------------*/
static int
intersectDefectConvexHullEdges(MRI_SURFACE *mris, DEFECT *defect, EDGE *e, 
															 int *vertex_trans)
{
  VERTEX *v ;
  int    i, n, vno ;
  FILE   *fp = NULL ;
  EDGE   edge2 ;

  if ((e->vno1 == 109259 && e->vno2 == 108332) ||
      (e->vno2 == 109259 && e->vno1 == 108332))
    DiagBreak() ;

  for (i = defect->nborder ; i < defect->nchull ; i++)
  {
		vno = vertex_trans[defect->chull[i]] ;

    if (vno == e->vno1 || vno == e->vno2 || vno < 0)
      continue ;

    edge2.vno1 = vno ;
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      if (vno == Gdiag_no)
        DiagBreak() ;

      if ((vno == 53153 && v->v[n] == 53158) ||
          (v->v[n] == 53153 && vno == 53158))
        DiagBreak() ;
      if ((vno == 52024 && v->v[n] == 53158) ||
          (v->v[n] == 52024 && vno == 53158))
        DiagBreak() ;
      if (v->v[n] == e->vno1 || v->v[n] == e->vno2)
        continue ;
      edge2.vno2 = v->v[n] ;
			if (edgesIntersect(mris, e, &edge2))
				return(1) ;
    }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    double     x1_start, x1_end, y1_start, y1_end, x2_start, x2_end, 
               y2_start,y2_end, x2min, x2max, y2min, y2max, cx, cy, cz, 
               origin[3], e0[3], e1[3] ;
    char       fname[STRLEN] ;
    VERTEX     *v2 ;
    static int dno = -1 ;

    mrisComputeCanonicalEdgeBasis(mris, e, e, origin, e0, e1) ;
    sprintf(fname, "lines%d.log", defect_no) ;
    v = &mris->vertices[e->vno1] ; v2 = &mris->vertices[e->vno2] ;
    cx = v->cx-origin[0] ; cy = v->cy - origin[1] ; cz = v->cz-origin[2];
    x1_start = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
    y1_start = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
    cx = v2->cx-origin[0]; cy = v2->cy - origin[1] ; cz = v2->cz-origin[2];
    x1_end = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
    y1_end = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
    if (dno != defect_no)  /* first time for this defect */
    {
      fp = fopen(fname, "w") ;

      for (i = 0 ; i < defect->nvertices+defect->nchull ; i++)
      {
        if (i < defect->nvertices)
          vno = vertex_trans[defect->vertices[i]] ;
        else
          vno = vertex_trans[defect->chull[i-defect->nvertices]] ;
        
        if (vno == e->vno1 || vno == e->vno2 || vno < 0)
          continue ;
        v = &mris->vertices[vno] ;
        cx = v->cx-origin[0] ; cy = v->cy-origin[1] ; cz = v->cz-origin[2];
        x2_start = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
        y2_start = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
        for (n = 0 ; n < v->vnum ; n++)
        {
          if (vno == 6167)
            DiagBreak() ;
          if (v->v[n] == e->vno1 || v->v[n] == e->vno2)
            continue ;
          v2 = &mris->vertices[v->v[n]] ;
          cx = v2->cx-origin[0] ; cy = v2->cy-origin[1]; cz = v2->cz-origin[2];
          x2_end = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
          y2_end = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
          x2min = MIN(x2_start, x2_end) ; x2max = MAX(x2_start, x2_end) ;
          y2min = MIN(y2_start, y2_end) ; y2max = MAX(y2_start, y2_end) ;
          fprintf(fp, "%2.3f  %2.3f\n%2.3f  %2.3f\n\n",
                  x2_start, y2_start, x2_end, y2_end) ;
        }
      }

    }
    else
      fp = fopen(fname, "a") ;
    fprintf(fp, "%2.3f  %2.3f\n%2.3f  %2.3f\n\n",
            x1_start, y1_start, x1_end, y1_end) ;
    if (fp)
      fclose(fp) ;


    sprintf(fname, "edges%d.log", defect_no) ;
    if (dno != defect_no)  /* first time for this defect */
    {
      dno = defect_no ;
      fp = fopen(fname, "w") ;

      for (i = 0 ; i < defect->nvertices+defect->nchull ; i++)
      {
        if (i < defect->nvertices)
          vno = vertex_trans[defect->vertices[i]] ;
        else
          vno = vertex_trans[defect->chull[i-defect->nvertices]] ;
        
        if (vno == e->vno1 || vno == e->vno2 || vno < 0)
          continue ;
        v = &mris->vertices[vno] ;
        cx = v->cx-origin[0] ; cy = v->cy-origin[1] ; cz = v->cz-origin[2];
        x2_start = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
        y2_start = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
        for (n = 0 ; n < v->vnum ; n++)
        {
          if (vno == 6167)
            DiagBreak() ;
          if (v->v[n] == e->vno1 || v->v[n] == e->vno2)
            continue ;
          v2 = &mris->vertices[v->v[n]] ;
          cx = v2->cx-origin[0] ; cy = v2->cy-origin[1]; cz = v2->cz-origin[2];
          x2_end = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
          y2_end = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
          x2min = MIN(x2_start, x2_end) ; x2max = MAX(x2_start, x2_end) ;
          y2min = MIN(y2_start, y2_end) ; y2max = MAX(y2_start, y2_end) ;
          fprintf(fp, "%d  %d  %2.3f  %2.3f  %2.3f  %2.3f\n", vno, v->v[n],
                  x2_start, y2_start, x2_end, y2_end) ;
        }
      }

    }
    else
      fp = fopen(fname, "a") ;
    fprintf(fp, "%d  %d  %2.3f  %2.3f  %2.3f  %2.3f\n", e->vno1, e->vno2,
            x1_start, y1_start, x1_end, y1_end) ;
    if (fp)
      fclose(fp) ;
  }
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          See if the edge e intersects any edges in the defect or
          it's border. Sorry this code is such a hatchet job. I'm sure
          there are far more elegant ways of doing intersection
          (e.g. sorting).
------------------------------------------------------*/
#if 0
static int
colinearDefectEdges(MRI_SURFACE *mris, DEFECT *defect, EDGE *e, 
                     int *vertex_trans)
{
  VERTEX *v ;
  int    i, n, vno ;
  EDGE   edge2 ;

  if ((e->vno1 == 109259 && e->vno2 == 108332) ||
      (e->vno2 == 109259 && e->vno1 == 108332))
    DiagBreak() ;

  for (i = 0 ; i < defect->nvertices+defect->nchull ; i++)
  {
    if (i < defect->nvertices)
      vno = vertex_trans[defect->vertices[i]] ;
    else
      vno = vertex_trans[defect->chull[i-defect->nvertices]] ;

		if (vno < 0)
			continue ;

    edge2.vno1 = vno ;
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      if (vno == Gdiag_no)
        DiagBreak() ;

      if ((vno == 53153 && v->v[n] == 53158) ||
          (v->v[n] == 53153 && vno == 53158))
        DiagBreak() ;
      if ((vno == 52024 && v->v[n] == 53158) ||
          (v->v[n] == 52024 && vno == 53158))
        DiagBreak() ;
			
      edge2.vno2 = v->v[n] ;
					
      if ((v->v[n] == e->vno1 || v->v[n] == e->vno2) ||
					(vno == e->vno1 || vno == e->vno2))  /* check for colinearity */
			{
				int vno0, vno1, vno2, ncoords ;
				VERTEX *v0, *v1, *v2 ;
				float  dx01, dy01, dz01, dx02, dy02, dz02, len01, len02 ;

				/*
					set vno0 and vno1 so that they are the existing edge, with
					vno0 being the shared vertex, and vno2 is the vertex for the
					putative edge
				*/
				if (vno == e->vno1 || vno == e->vno2)  /* vno is shared vertex */
				{
					vno0 = vno ; vno1 = v->v[n] ;
				}
				else   /* v->v[n] is shared vertex */
				{
					vno0 = v->v[n] ; vno1 = vno ;
				}
				if (e->vno1 == vno0)
					vno2 = e->vno2 ;
				else
					vno2 = e->vno1 ;
				v0 = &mris->vertices[vno0] ; v1 = &mris->vertices[vno1] ;
				v2 = &mris->vertices[vno2] ; 
				dx01 = v1->origx - v0->origx ; dy01 = v1->origy - v0->origy ;
				dz01 = v1->origz - v0->origz ; len01 = sqrt(dx01*dx01+dy01*dy01+dz01*dz01) ;
				if (FZERO(len01))
					len01 = 1 ;
				dx01 /= len01 ; dy01 /= len01 ; dz01 /= len01 ;
				dx02 = v2->origx - v0->origx ; dy02 = v2->origy - v0->origy ;
				dz02 = v2->origz - v0->origz ; len02 = sqrt(dx02*dx02+dy02*dy02+dz02*dz02) ;
				if (FZERO(len02))
					len02 = 1 ;
				dx02 /= len02 ; dy02 /= len02 ; dz02 /= len02 ;
				

				/* see if vno1 and vno2 are colinear. If not, no problemo */
				ncoords = FZERO(dx01-dx02)+FZERO(dy01-dy02)+FZERO(dz01-dz02);
				if (ncoords == 3)
				{
					if (e->vno1 == 16968 || e->vno2 == 16968)
						DiagBreak() ;
					return(1) ;   /* colinear */
				}
			}
    }
  }

  return(0) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Is vno1 is a neighbor of vno2
------------------------------------------------------*/
static int
vertexNeighbor(MRI_SURFACE *mris, int vno1, int vno2)
{
  VERTEX *v ;
  int    n ;

  v = &mris->vertices[vno1] ;
  for (n = 0 ; n < v->vnum ; n++)
    if (v->v[n] == vno2)
      return(1) ;
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Add the edge vnot <--> vno2 to the tessellation
------------------------------------------------------*/
#define MAX_VLIST 10000
static int
mrisAddEdge(MRI_SURFACE *mris, int vno1, int vno2)
{
  int    vlist[MAX_VLIST] ;
  VERTEX *v ;

  if (vno1 < 0 || vno2 < 0)
    DiagBreak() ;
  if (vno1 == Gdiag_no || vno2 == Gdiag_no)
    DiagBreak() ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "adding edge %d <--> %d\n", vno1, vno2) ;

  if (v->vnum >= MAX_VLIST-1)
    ErrorExit(ERROR_NOMEMORY, "mrisAddEdge: too many edges (%d)",v->vnum) ;

  /* add v2 link to v1 struct */
  v = &mris->vertices[vno1] ;
  memcpy(vlist, v->v, v->vnum*sizeof(int)) ;
  vlist[(unsigned int)v->vnum++] = vno2 ; v->vtotal = v->vnum ;
  if (v->v)
    free(v->v) ;
  v->v = (int *)calloc(v->vnum, sizeof(int)) ;
  if (!v->v)
    ErrorExit(ERROR_NO_MEMORY, 
              "mrisAddEdge(%d, %d): could not allocate %d len vlist", v->vnum);

  memcpy(v->v, vlist, v->vnum*sizeof(int)) ;

  /* add v1 link to v2 struct */
  v = &mris->vertices[vno2] ;
  memcpy(vlist, v->v, v->vnum*sizeof(int)) ;
  vlist[(unsigned int)v->vnum++] = vno1 ; v->vtotal = v->vnum ;
  if (v->v)
    free(v->v) ;
  v->v = (int *)calloc(v->vnum, sizeof(int)) ;
  if (!v->v)
    ErrorExit(ERROR_NO_MEMORY, 
              "mrisAddEdge(%d, %d): could not allocate %d len vlist", v->vnum);

  memcpy(v->v, vlist, v->vnum*sizeof(int)) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Add the triangle vno0 --> vno1 --> vno2 to the tessellation
------------------------------------------------------*/
static int
mrisAddFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2)
{
  int    n, fno, ilist[1000], n0, n1 ;
  FACE   *f ;
  VERTEX *v ;
  uchar  ulist[1000] ;
  float  norm[3], dot, cx, cy, cz ;

  if (vno0 < 0 || vno1 < 0 || vno2 < 0)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mrisAddFace(%d,%d,%d)!\n",
                                vno0, vno1, vno2)) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "adding face (%d, %d, %d)\n", vno0, vno1, vno2) ;
  if (vno0 == 6160 && vno1 == 6189 && vno2 == 6176)
    DiagBreak() ;
  if (mris->nfaces >= mris->max_faces)
    ErrorExit(ERROR_NOMEMORY, "mrisAddFace: max faces reached") ;

  fno = mris->nfaces++ ;
  if (fno == Gdiag_no)
    DiagBreak() ;

  f = &mris->faces[fno] ;
  f->v[0] = vno0 ; f->v[1] = vno1 ; f->v[2] = vno2 ; 

  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    v = &mris->vertices[f->v[n]] ;
    if (v->num >= 255)
      continue ;
    memcpy(ilist, v->f, v->num*sizeof(int)) ;
    ilist[v->num++] = fno ;
    if (v->f)
      free(v->f) ;
    v->f = (int *)calloc(v->num, sizeof(int)) ;
    if (!v->f)
      ErrorExit(ERROR_NOMEMORY, "mrisAddFace: could not allocate face list") ;
    memcpy(v->f, ilist, v->num*sizeof(int)) ;

    memcpy(ulist, v->n, (v->num-1)*sizeof(uchar)) ;
    ulist[v->num-1] = n ;
    if (v->n)
      free(v->n) ;
    v->n = (uchar *)calloc(v->num, sizeof(uchar)) ;
    if (!v->n)
      ErrorExit(ERROR_NOMEMORY, "mrisAddFace: could not allocate n list") ;
    memcpy(v->n, ulist, v->num*sizeof(uchar)) ;
  }

  mrisCalculateCanonicalFaceCentroid(mris, fno, &cx, &cy, &cz);
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    mrisNormalFace(mris, fno, n, norm) ;  /* compute face normal */
    dot = norm[0]*cx + norm[1]*cy + norm[2]*cz ;
    
    if (dot < 0) /* they disagree - change order of vertices in face */
    {
      n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
      n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
      vno0 = f->v[n0] ; vno1 = f->v[n1] ;
      f->v[n0] = vno1 ; f->v[n1] = vno0 ;
      mrisSetVertexFaceIndex(mris, vno0, fno) ;
      mrisSetVertexFaceIndex(mris, vno1, fno) ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Is the triangle vno0-->vno1-->vno2 already in the tessellation
------------------------------------------------------*/
static int
isFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2)
{
  VERTEX  *v ;
  FACE    *f ;
  int     n, n1, vno ;

  v = &mris->vertices[vno0] ;
  for (n = 0 ; n < v->num ; n++)
  {
    f = &mris->faces[v->f[n]] ;
    for (n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
    {
      vno = f->v[n1] ;
      if (vno != vno0 && vno != vno1 && vno != vno2)
        break ;
    }
    if (n1 >= VERTICES_PER_FACE)
      return(1) ;
  }
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Is the triangle vno0-->vno1-->vno2 already in the tessellation
------------------------------------------------------*/
static int
findFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2)
{
  VERTEX  *v ;
  FACE    *f ;
  int     n, n1, vno, fno ;

  v = &mris->vertices[vno0] ;
  for (n = 0 ; n < v->num ; n++)
  {
    f = &mris->faces[v->f[n]] ;
    for (n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
    {
      vno = f->v[n1] ;
      fno  = v->f[n] ;
      if (vno != vno0 && vno != vno1 && vno != vno2)
        break ;
    }
    if (n1 >= VERTICES_PER_FACE)
      return(fno) ;
  }
  return(-1) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Orient the faces of the tessellation so that the
          point outward (i.e. in the same direction as the original
          surface).
------------------------------------------------------*/
static int
mrisOrientFaces(MRI_SURFACE *mris, DEFECT *defect, int *vtrans)
{
  int    vno, vno0, vno1, i, n, fno, oriented = 0 ;
  FACE   *f ;
  VERTEX *v, *vn ;
  float  norm[3], dot ;

  /* first orient vertices */
  for (i = 0 ; i < defect->nvertices ; i++)
  {
    vno = vtrans[defect->vertices[i]] ;
    if (vno < 0)   /* not in new tessellation */
      continue ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (dot = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dot += v->nx*vn->nx + v->ny*vn->ny + v->nz*vn->nz ;
    }
    if (dot < 0)
    { 
      oriented++ ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "reversing normal for vertex %d\n", vno) ;
      v->nx *= -1.0f ; v->ny *= -1.0f ; v->nz *= -1.0f ; 
    }
  }
  /*  mrisRestoreDefectPositions(mris, defect, ORIGINAL_VERTICES) ;*/
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
      vno = vtrans[defect->vertices[i]] ;
    else
      vno = vtrans[defect->border[i-defect->nvertices]] ;
    if (vno < 0)   /* not in new tessellation */
      continue ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < v->num ; n++)
    {
      fno = v->f[n] ;
      if (mris->faces[fno].ripflag)  /* only consider it once */
        continue ;
      if (fno == Gdiag_no)
        DiagBreak() ;
      mris->faces[fno].ripflag = 1 ;
      mrisNormalFace(mris, fno, v->n[n], norm) ;
      dot = norm[0]*v->nx + norm[1]*v->ny + norm[2]*v->nz ;
      if (dot < 0)   /* change order of vertices in face */
      {
        oriented++ ;
        f = &mris->faces[fno] ;
        vno0 = f->v[0] ; vno1 = f->v[1] ;
        f->v[0] = vno1 ; f->v[1] = vno0 ;
        mrisSetVertexFaceIndex(mris, vno0, fno) ;
        mrisSetVertexFaceIndex(mris, vno1, fno) ;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "reversing face %d orientation\n", fno) ;
      }
    }
  }
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
      vno = vtrans[defect->vertices[i]] ;
    else
      vno = vtrans[defect->border[i-defect->nvertices]] ;
    if (vno < 0)   /* not in new tessellation */
      continue ;
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->num ; n++)
    {
      fno = v->f[n] ;
      mris->faces[fno].ripflag = 0 ;
    }
  }
  /*  mrisRestoreDefectPositions(mris, defect, TMP_VERTICES) ;*/
  return(oriented) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Set the positions of all vertices in a defect to the
          original, canonical, or tmp vertices.
------------------------------------------------------*/
static int
mrisRestoreDefectPositions(MRI_SURFACE *mris, DEFECT *defect, int which)
{
  int     vno, i ;
  VERTEX  *v ;

  for (i = 0 ; i < defect->nvertices ; i++)
  {
    vno = defect->vertices[i] ; v = &mris->vertices[vno] ;
#if 0
    if (v->ripflag)
      continue ;
#endif
    switch (which)
    {
    case CANONICAL_VERTICES:
      v->x = v->cx ; v->y = v->cy ; v->z = v->cz ;
      break ;
    case ORIGINAL_VERTICES:
      v->x = v->origx ; v->y = v->origy ; v->z = v->origz ;
      break ;
    default:
    case TMP_VERTICES:
      v->x = v->tx ; v->y = v->ty ; v->z = v->tz ;
      break ;
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Search the face for vno and set the v->n[] field
          appropriately.
------------------------------------------------------*/
static int
mrisSetVertexFaceIndex(MRI_SURFACE *mris, int vno, int fno)
{
  VERTEX  *v ;
  FACE    *f ;
  int     n, i ;

  v = &mris->vertices[vno] ;
  f = &mris->faces[fno] ;

  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    if (f->v[n] == vno)
      break ;
  }
  if (n >= VERTICES_PER_FACE)
    return(ERROR_BADPARM) ;

  for (i = 0 ; i < v->num ; i++)
    if (v->f[i] == fno)
      v->n[i] = n ;

  return(n) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Add the appropriate face to the tessellation.
------------------------------------------------------*/
static int
mrisAddDefectFaces(MRI_SURFACE *mris, double e0[3], double e1[3], 
                   double origin[3], EDGE *et, int nedges)
{
  int     i, vno1, vno2, nfaces, n ;
  VERTEX  *v ;

  for (i = 0 ; i < nedges ; i++)
  {
    if (et[i].used == 0)
      continue ;
    vno1 = et[i].vno1 ;
    vno2 = et[i].vno2 ;
    mrisComputeCanonicalEdgeBasis(mris, et+i, et+i, origin, e0, e1) ;
    if (vno1 == 108332 && vno2 == 109240)
      DiagBreak() ;

    /* for every vertex which is a neighbor of both of these, add 1 triangle */
    v = &mris->vertices[vno1] ;
    for (nfaces = n = 0 ; n < v->vnum ; n++)
    {
      if (v->v[n] == vno2)
        continue ;
      if (vertexNeighbor(mris, vno2, v->v[n]) && 
          !isFace(mris,vno1, vno2, v->v[n]) &&
          !containsAnotherVertex(mris,vno1,vno2,v->v[n],e0,e1,origin))
      {
        if (nfaces++ > 1)
          DiagBreak() ;
        if (mris->nfaces == Gdiag_no)
          DiagBreak() ;
#if 0
        if (((vno1 != 110718) && (vno1 != 110732) && (vno1 != 110748)) ||
            ((vno2 != 110718) && (vno2 != 110732) && (vno2 != 110748)) ||
            ((v->v[n] != 110718) && (v->v[n] != 110732) && (v->v[n]!=110748)))
#endif
          mrisAddFace(mris, vno1, vno2, v->v[n]) ;
      }
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Add the appropriate face to the tessellation.
------------------------------------------------------*/
static int
mrisAddAllDefectFaces(MRI_SURFACE *mris, DEFECT_LIST *dl, int *vertex_trans)
{
  DEFECT  *defect ;
  int     i, vno1, vno2, nfaces, n, m, dno ;
  VERTEX  *v ;
  double  origin[3], e0[3], e1[3] ;
  EDGE    edge ;
  
  for (dno = 0 ; dno < dl->ndefects ; dno++)
  {
    defect = &dl->defects[dno] ;
    for (i = 0 ; i < defect->nborder+defect->nvertices ; i++)
    {
      if (i < defect->nvertices)
        vno1 = vertex_trans[defect->vertices[i]] ;
      else
        vno1 = vertex_trans[defect->border[i-defect->nvertices]] ;
      if (vno1 < 0)
        continue ;

      v = &mris->vertices[vno1] ;
      edge.vno1 = vno1 ;
      for (m = 0 ; m < v->vnum ; m++)
      {
        edge.vno2 = vno2 = v->v[m] ;
        mrisComputeCanonicalEdgeBasis(mris, &edge, &edge, origin, e0, e1) ;
        if (vno1 == 108332 && vno2 == 109240)
          DiagBreak() ;

        /* 
           for every vertex which is a neighbor of both of these, 
           add 1 triangle */
        for (nfaces = n = 0 ; n < v->vnum ; n++)
        {
          if (v->v[n] == vno2)
            continue ;
          if (vertexNeighbor(mris, vno2, v->v[n]) && 
              !isFace(mris,vno1, vno2, v->v[n]) &&
              !containsAnotherVertex(mris,vno1,vno2,v->v[n],e0,e1,origin))
          {
            if (nfaces++ > 1)
              DiagBreak() ;
            if (mris->nfaces == Gdiag_no)
              DiagBreak() ;
            mrisAddFace(mris, vno1, vno2, v->v[n]) ;
          }
        }
      }
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          See if this triangle contains another vertex within it.
          If so, it should not be added to the tessellation.
------------------------------------------------------*/
static int
containsAnotherVertex(MRI_SURFACE *mris, int vno0, int vno1, int vno2,
                    double e0[3], double e1[3], double origin[3])
#if 1
{
  int     n, vno, i, n1 ;
  VERTEX  *vn, *v0, *v1, *v2, *v ;
  double  cx, cy, cz, x0, y0, x1, y1, x2, y2, xn, yn, m, b ;
  FILE    *fp ;


  v = &mris->vertices[vno0] ;
  
  if (((vno1 == 110718) || (vno1 == 110732) || (vno1 == 110748)) &&
      ((vno2 == 110718) || (vno2 == 110732) || (vno2 == 110748)) &&
      ((vno0 == 110718) || (vno0 == 110732) || (vno0 == 110748)))
  {
    fp = fopen("triangle.log", "w") ;
    for (n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
    {
      switch (n1)
      {
      default:
      case 0: v = &mris->vertices[vno0] ; break ;
      case 1: v = &mris->vertices[vno1] ; break ;
      case 2: v = &mris->vertices[vno2] ; break ;
      }
      for (n = 0 ; n < v->vnum ; n++)
      {
        vno = v->v[n] ;
        if (vno == vno0 || vno == vno1 || vno == vno2)
          continue ;
        vn = &mris->vertices[vno] ;
        cx = vn->cx-origin[0]; cy = vn->cy - origin[1]; cz = vn->cz -origin[2];
        xn = cx*e0[0] + cy*e0[1] + cz*e0[2]; yn = cx*e1[0] + cy*e1[1]+cz*e1[2];
        fprintf(fp, "# vertex %d\n", vno) ;
        fprintf(fp, "%f %f\n\n", xn, yn) ;
      }
    }
    fprintf(fp, "# vertices %d %d %d\n", vno0, vno1, vno2) ;
    v0 = &mris->vertices[vno0] ;
    v1 = &mris->vertices[vno1] ;
    v2 = &mris->vertices[vno2] ;
    
    cx = v0->cx-origin[0]; cy = v0->cy - origin[1]; cz = v0->cz - origin[2];
    x0 = cx*e0[0] + cy*e0[1] + cz*e0[2]; y0 = cx*e1[0] + cy*e1[1] + cz*e1[2];
    
    cx = v1->cx-origin[0]; cy = v1->cy - origin[1]; cz = v1->cz - origin[2];
    x1 = cx*e0[0] + cy*e0[1] + cz*e0[2]; y1 = cx*e1[0] + cy*e1[1] + cz*e1[2];
    
    cx = v2->cx-origin[0]; cy = v2->cy - origin[1]; cz = v2->cz - origin[2];
    x2 = cx*e0[0] + cy*e0[1] + cz*e0[2]; y2 = cx*e1[0] + cy*e1[1] + cz*e1[2];
      
    fprintf(fp, "%f %f\n%f %f\n%f %f\n%f %f\n\n", 
            x0, y0, x1, y1, x2, y2, x0, y0) ;
    fclose(fp) ;
  }

  for (n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
  {
    switch (n1)
    {
    default:
    case 0: v = &mris->vertices[vno0] ; break ;
    case 1: v = &mris->vertices[vno1] ; break ;
    case 2: v = &mris->vertices[vno2] ; break ;
    }
    for (n = 0 ; n < v->vnum ; n++)
    {
      vno = v->v[n] ;
      if (vno == vno0 || vno == vno1 || vno == vno2)
        continue ;
      
      vn = &mris->vertices[v->v[n]] ;
      cx = vn->cx-origin[0]; cy = vn->cy - origin[1]; cz = vn->cz - origin[2];
      xn = cx*e0[0] + cy*e0[1] + cz*e0[2]; yn = cx*e1[0] + cy*e1[1] + cz*e1[2];
      for (i = 0 ; i < VERTICES_PER_FACE ; i++)
      {
        switch (i)
        {
        default:
        case 0:
          v0 = &mris->vertices[vno0] ;
          v1 = &mris->vertices[vno1] ;
          v2 = &mris->vertices[vno2] ;
          break ;
        case 1:
          v0 = &mris->vertices[vno1] ;
          v1 = &mris->vertices[vno2] ;
          v2 = &mris->vertices[vno0] ;
          break ;
        case 2:
          v0 = &mris->vertices[vno2] ;
          v1 = &mris->vertices[vno0] ;
          v2 = &mris->vertices[vno1] ;
          break ;
        }
        cx = v0->cx-origin[0]; cy = v0->cy - origin[1]; cz = v0->cz-origin[2];
        x0 = cx*e0[0] + cy*e0[1] + cz*e0[2]; y0 = cx*e1[0]+cy*e1[1] + cz*e1[2];
        
        cx = v1->cx-origin[0]; cy = v1->cy - origin[1]; cz = v1->cz-origin[2];
        x1 = cx*e0[0] + cy*e0[1] + cz*e0[2]; y1 = cx*e1[0]+cy*e1[1] + cz*e1[2];
        
        cx = v2->cx-origin[0]; cy = v2->cy - origin[1]; cz = v2->cz-origin[2];
        x2 = cx*e0[0] + cy*e0[1] + cz*e0[2]; y2 = cx*e1[0]+cy*e1[1] + cz*e1[2];
        
        if (FEQUAL(x1, x0))   /* vertical leg */
        {
          if (((x2 - x1) * (xn - x1)) < 0)  /* on opposite sides of leg */
            break ;
        }
        else
        {
          m = (y1 - y0) / (x1 - x0) ;
          b = y1 - m * x1 ;
          if ((y2 - (m * x2 + b)) * (yn - (m * xn + b)) < 0)
            break ;
        }
      }
      if (i >= VERTICES_PER_FACE)  /* all inside */
        return(1) ;
    }
  }
  return(0) ;
}
#else
{
  int     n, vno, i ;
  VERTEX  *v0, *v1, *v2, *v ;
  double   U0[3], U1[3], U2[3], dot, L0[3], L1[3], *V0, *V1, *V2,
          desc[3], cx, cy, cz, Point[3], len, norm_proj[3], U[3] ;

  /* compute planar coordinate representation of triangle vertices */
  v0 = &mris->vertices[vno0] ;
  v1 = &mris->vertices[vno1] ;
  v2 = &mris->vertices[vno2] ;
  cx = v0->cx-origin[0] ; cy = v0->cy-origin[1] ; cz = v0->cz-origin[2];
  U0[0] = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
  U0[1] = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  U0[2] = 0 ;
  cx = v1->cx-origin[0] ; cy = v1->cy-origin[1] ; cz = v1->cz-origin[2];
  U1[0] = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
  U1[1] = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  U1[2] = 0 ;
  cx = v2->cx-origin[0] ; cy = v2->cy-origin[1] ; cz = v2->cz-origin[2];
  U2[0] = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
  U2[1] = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  U2[2] = 0 ;

  for (n = 0 ; n < v0->vnum ; n++)
  {
    vno = v0->v[n] ;
    if (vno == vno1 || vno == vno2)
      continue ;
    v = &mris->vertices[vno] ;
    cx = v->cx-origin[0] ; cy = v->cy-origin[1] ; cz = v->cz-origin[2];
    U[0] = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
    U[1] = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
    U[2] = 0 ;

    for (i = 0 ; i < 3 ; i++)
    {
      /* 
         build a coordinate system with V0 as the origin, then construct
         the vector connecting V2 with it's normal projection onto V0->V1.
         This will be a descriminant vector for dividing the plane by the
         V0->V1 line. A positive dot product with the desc. vector indicates
         that the point is on the positive side of the plane and therefore
         may be contained within the triangle. Doing this for each of the
         legs in sequence gives a test for being inside the triangle.
      */
      
      switch (i)
      {
      default:
      case 0:   V0 = U0 ; V1 = U1 ; V2 = U2 ; break ;
      case 1:   V0 = U1 ; V1 = U2 ; V2 = U0 ; break ;
      case 2:   V0 = U2 ; V1 = U0 ; V2 = U1 ; break ;
      }
      SUB(L0, V1, V0) ; SUB(L1, V2, V0) ;
      
      /* compute normal projection onto base of triangle */
      len = VLEN(L0) ; L0[0] /= len ; L0[1] /= len ; L0[2] /= len ;
      dot = DOT(L0,L1) ; 
      SCALAR_MUL(norm_proj, dot, L0) ;
      
      /* build descriminant vector */
      SUB(desc, L1, norm_proj) ;
      
      /* 
         transform point in question into local coordinate system and build
         the vector from the point in question to the normal projection point.
         The dot product of this vector with the descrimant vector will then
         indicate which side of the V0->V1 line the point is on.
      */
      SUB(Point, U, V0) ; SUB(Point, Point, norm_proj) ;
      dot = DOT(desc, Point) ;
      if (dot < 0 && !DZERO(dot))  /* not in triangle */
        break ;
    }
    if (i >= 3)  /* contained in triangle */
      return(1) ;
  }

  return(0) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Orient the faces of the tessellation so that the
          point outward (i.e. in the same direction as the original
          surface).
------------------------------------------------------*/
#if 1
static int
mrisOrientRetessellatedSurface(MRI_SURFACE *mris, DEFECT_LIST *dl,int *vtrans)
{
  int     vno, n, fno, m, n0, n1, vno0, vno1 ;
  VERTEX  *v ;
  FACE    *f ;
  float   dot, norm[3], cx, cy, cz ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISprojectOntoSphere(mris, mris, MRISaverageRadius(mris)) ;

  MRIScomputeMetricProperties(mris) ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    if (fno == Gdiag_no)
      DiagBreak() ;
    f = &mris->faces[fno] ;
    mrisCalculateFaceCentroid(mris, fno, &cx, &cy, &cz) ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      mrisNormalFace(mris, fno, n, norm) ;  /* compute face normal */
      dot = norm[0]*cx + norm[1]*cy + norm[2]*cz ;
      
      if (dot < 0) /* they disagree - change order of vertices in face */
      {
        n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
        n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
        vno0 = f->v[n0] ; vno1 = f->v[n1] ;
        f->v[n0] = vno1 ; f->v[n1] = vno0 ;
        mrisSetVertexFaceIndex(mris, vno0, fno) ;
        mrisSetVertexFaceIndex(mris, vno1, fno) ;
      }
    }
  }

  MRIScomputeMetricProperties(mris) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dot = v->nx * v->x + v->ny*v->y + v->nz*v->z ;
    if (dot < 0)
    {
      fprintf(stdout, "vertex %d has inverted normal!\n", vno) ;
      DiagBreak() ;
    }

    for (m = 0 ; m < v->num ; m++)
    {
      fno = v->f[m] ; f = &mris->faces[fno] ;
      if (fno == Gdiag_no)
        DiagBreak() ;
      if (f->area < 0)
      {
        fprintf(stderr, "face %d has negative area!\n", fno) ;
        DiagBreak() ;
      }
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      {
        mrisNormalFace(mris, fno, n, norm) ;  /* compute face normal */
        dot = norm[0]*v->x + norm[1]*v->y + norm[2]*v->z ;
        
        if (dot < 0) /* they disagree - change order of vertices in face */
        {
          fprintf(stderr, "face %d has inverted normal!\n", fno) ;
          DiagBreak() ;
        }
      }
    }
  }

  MRISclearMarks(mris) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  return(NO_ERROR) ;
}
#else
static int
mrisOrientRetessellatedSurface(MRI_SURFACE *mris, DEFECT_LIST *dl,int *vtrans)
{
#if 1
  int     dno, fno, vno, i, n, m, vno0, vno1, n0, n1, oriented, nreversed ;
  VERTEX  *v, *vn ;
  FACE    *f ;
  DEFECT  *defect ;
  float   dot, norm[3], len, nx, ny, nz ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIG_VERTICES) ;
  MRISaverageVertexPositions(mris, 200) ;
  MRIScomputeMetricProperties(mris) ;

  /* compute average boundary normal in the smoothed space */
  for (dno = 0 ; dno < dl->ndefects ; dno++)
  {
    defect = &dl->defects[dno] ;
    nx = ny = nz = 0.0f ;
    for (i = 0 ; i < defect->nborder ; i++)
    {
      vno = vtrans[defect->border[i]] ;
      if (vno < 0)        /* not in new tessellation */
        continue ;
      v = &mris->vertices[vno] ;
      nx += v->nx ; ny += v->ny ; nz += v->nz ; 
    }
    len = sqrt(nx*nx + ny*ny + nz*nz) ;
    if (FZERO(len))
      len = 1.0f ;
    defect->nx = nx/len ; defect->ny = ny/len ; defect->nz = nz/len ;
  }

  /* orient defect faces so that the outward direction agrees with boundary */
  for (oriented = dno = 0 ; dno < dl->ndefects ; dno++)
  {
    defect = &dl->defects[dno] ;
    for (i = 0 ; i < defect->nvertices ; i++)
    {
      vno = vtrans[defect->vertices[i]] ;
      if (vno < 0)     /* not in new tessellation */
        continue ;
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;

      /* go through each face and orient it to agree with defect normal */
      for (m = 0 ; m < v->num ; m++)
      {
        fno = v->f[m] ; f = &mris->faces[fno] ;
        if (fno == Gdiag_no)
          DiagBreak() ;
        for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        {
          mrisNormalFace(mris, fno, n, norm) ;  /* compute face normal */
          dot = norm[0]*defect->nx + norm[1]*defect->ny + norm[2]*defect->nz ;
          if (dot < 0)   /* they disagree - change order of vertices in face */
          {
            oriented++ ;
            n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
            n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
            vno0 = f->v[n0] ; vno1 = f->v[n1] ;
            f->v[n0] = vno1 ; f->v[n1] = vno0 ;
            mrisSetVertexFaceIndex(mris, vno0, fno) ;
            mrisSetVertexFaceIndex(mris, vno1, fno) ;
            if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
              fprintf(stdout, "reversing face %d orientation\n", fno) ;
          }
        }
      }
    }
  }

  MRISsetNeighborhoodSize(mris, 2) ;

  i = 0 ;
  do
  {
    MRIScomputeMetricProperties(mris) ;
    nreversed = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      for (nx = ny = nz = 0.0, n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        dot = vn->nx*v->nx + vn->ny*v->ny + vn->nz*v->nz ;
        if (dot < 0)
          DiagBreak() ;
        nx += vn->nx ; ny += vn->ny ; nz += vn->nz ; 
      }
      dot = nx*v->nx + ny*v->ny + nz*v->nz ;
      if (dot < 0)
      {
        v->nx *= -1.0 ; v->ny *= -1.0 ; v->nz *= -1.0 ; 
        DiagBreak() ;
      }
      
      for (m = 0 ; m < v->num ; m++)
      {
        fno = v->f[m] ; f = &mris->faces[fno] ;
        if (fno == Gdiag_no)
          DiagBreak() ;
        dot = nx*f->nx + ny*f->ny + nz*f->nz ;
        if (dot < 0)   /* they disagree - change order of vertices in face */
          DiagBreak() ;
        for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        {
          mrisNormalFace(mris, fno, n, norm) ;  /* compute face normal */
          dot = norm[0]*nx + norm[1]*ny + norm[2]*nz ;
          if (dot < 0)   /* they disagree - change order of vertices in face */
          {
            nreversed++ ;
            n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
            n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
            vno0 = f->v[n0] ; vno1 = f->v[n1] ;
            f->v[n0] = vno1 ; f->v[n1] = vno0 ;
            mrisSetVertexFaceIndex(mris, vno0, fno) ;
            mrisSetVertexFaceIndex(mris, vno1, fno) ;
            if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
              fprintf(stdout, "reversing face %d orientation\n", fno) ;
          }
        }
      }
    }
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "\rpass %d: %d oriented   ", i+1, nreversed) ;

    oriented += nreversed ;
    if (++i > 20)  /* shouldn't happen, but... */
      break ;
  } while (nreversed > 0) ;


  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "orientation complete in %d passes\n", i) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (m = 0 ; m < v->num ; m++)
    {
      fno = v->f[m] ; f = &mris->faces[fno] ;
      if (fno == Gdiag_no)
        DiagBreak() ;
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      {
        mrisNormalFace(mris, fno, n, norm) ;  /* compute face normal */
        dot = norm[0]*v->nx + norm[1]*v->ny + norm[2]*v->nz ;
        if (dot < 0)   /* they disagree - change order of vertices in face */
          DiagBreak() ;
        dot = norm[0]*f->nx + norm[1]*f->ny + norm[2]*f->nz ;
        if (dot < 0)   /* they disagree - change order of vertices in face */
          DiagBreak() ;
      }
    }
  }

#if 0
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;  /* back to inflated */
#endif
  MRIScomputeMetricProperties(mris) ;
  return(oriented) ;
#else
  int    vno, vno0, vno1, i, n, fno, oriented, num, dno, m, blist[200000], nb,
         tmp[200000], nbnew, n0, n1 ;
  FACE   *f ;
  VERTEX *v, *vn ;
  float  norm[3], dot, len ;
  DEFECT *defect ;

  oriented = 0 ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MRISsetNeighborhoodSize(mris, 2) ;

  /* first orient vertices */

  /* mark all defective vertices */
  for (nb = dno = 0 ; dno < dl->ndefects ; dno++)
  {
    defect = &dl->defects[dno] ;
    for (n = 0 ; n < defect->nvertices ; n++)
    {
      vno = vtrans[defect->vertices[n]] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (vno < 0)
        continue ;
      mris->vertices[vno].marked = 1 ;
    }
    for (n = 0 ; n < defect->nborder ; n++)
      mris->vertices[vtrans[defect->border[n]]].marked = 1 ;
    memcpy(blist+nb, defect->border, defect->nborder*sizeof(int)) ;
    nb += defect->nborder ;
  }

  /* starts out as a list of border vertices and will grow inwards */
  for (i = 0 ; i < nb ; i++)
    blist[i] = vtrans[blist[i]] ;

  do   /* grow border inwards one edge length at each iteration */
  {
    for (nbnew = i = 0 ; i < nb ; i++)
    {
      vno = blist[i] ;
      if (vno < 0)   /* not in new tessellation - shouldn't happen */
        continue ;
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no || vno == 135681)
        DiagBreak() ;
      v->nx = v->ny = v->nz = 0 ;
      for (dot = 0.0, num = n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
          
        /* only use already oriented (or non-defective) vertices */
        if (vn->marked)
          continue ; 
        v->nx += vn->nx ; v->ny += vn->ny ; v->nz += vn->nz ; 
        num++ ;
      }
      if (!num)   /* surrounded by unoriented defective vertices */
      {
        tmp[nbnew++] = vno ;   /* so it will be processed next time again */
        continue ;
      }
      len = sqrt(v->nx*v->nx + v->ny*v->ny + v->nz*v->nz) ;
      if (FZERO(len))
        len = 1.0f ;
      v->nx /= len ; v->ny /= len ; v->nz /= len ;
      v->marked = 3 ;   /* it's proper orientation has been established */
    }
    for (i = 0 ; i < nb ; i++)
    {
      vno = blist[i] ;
      if (vno < 0)   /* not in new tessellation - shouldn't happen */
        continue ;
      v = &mris->vertices[vno] ;
      if (v->marked == 3)  /* was oriented properly */
        v->marked = 0 ;
    }

    /* now build new list of vertices, moving inward by one vertex */
    for (i = 0 ; i < nb ; i++)
    {
      vno = blist[i] ;
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
          
        /* only use already oriented (or non-defective) vertices */
        if (vn->marked == 1)
        {
          vn->marked = 2 ;  /* don't add it more than once */
          tmp[nbnew++] = v->v[n] ;
        }
      }
    }
    nb = nbnew ;
    if (nb > 0)
      memcpy(blist, tmp, nb*sizeof(int)) ;
  } while (nb > 0) ;


  MRISclearMarks(mris) ;

  /* now orient faces to agree with their vertices */
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      mrisNormalFace(mris, fno, n, norm) ;  /* how about vertex 2 ???? */
      vno = f->v[n] ; v = &mris->vertices[vno] ;
      dot = norm[0]*v->nx + norm[1]*v->ny + norm[2]*v->nz ;
      if (dot < 0)   /* change order of vertices in face */
      {
        oriented++ ;
        n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
        n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
        vno0 = f->v[n0] ; vno1 = f->v[n1] ;
        f->v[n0] = vno1 ; f->v[n1] = vno0 ;
        mrisSetVertexFaceIndex(mris, vno0, fno) ;
        mrisSetVertexFaceIndex(mris, vno1, fno) ;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "reversing face %d orientation\n", fno) ;
      }
    }
  }

  MRIScomputeMetricProperties(mris) ;
  MRISclearMarks(mris) ;

  /* mark all vertices that have a normal which disagrees with any neighbor */
  for (dno = 0 ; dno < dl->ndefects ; dno++)
  {
    defect = &dl->defects[dno] ;
    for (i = 0 ; i < defect->nvertices ; i++)
    {
      vno = vtrans[defect->vertices[i]] ;
      if (vno < 0)
        continue ;
      v = &mris->vertices[vno] ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        dot = vn->nx * v->nx + vn->ny * v->ny + vn->nz * v->nz ;
        if (dot < 0)
        { 
          v->marked = 1 ;
          if (!vn->marked)
            vn->marked = 1 ; 
          break ;
        }
      }
    }
  }

  /* go back and orient the ambiguous vertices based on the normal of
     the defect.
  */
  for (dno = 0 ; dno < dl->ndefects ; dno++)
  {
    defect = &dl->defects[dno] ;
    for (i = 0 ; i < defect->nvertices ; i++)
    {
      vno = vtrans[defect->vertices[i]] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (vno < 0)
        continue ;
      v = &mris->vertices[vno] ;
      if (!v->marked)
        continue ;
      dot = defect->nx * v->nx + defect->ny * v->ny + defect->nz * v->nz ;
      if (dot < 0)
      {
        v->nx *= -1 ; v->ny *= -1 ; v->nz *= -1 ;
      }
    }
  }
  MRISclearMarks(mris) ;

  /* now orient faces to agree with their vertices */
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      mrisNormalFace(mris, fno, n, norm) ;  /* how about vertex 2 ???? */
      vno = f->v[n] ; v = &mris->vertices[vno] ;
      dot = norm[0]*v->nx + norm[1]*v->ny + norm[2]*v->nz ;
      if (dot < 0)   /* change order of vertices in face */
      {
        oriented++ ;
        m = (n+1) >= VERTICES_PER_FACE ? 0 : n+1 ;
        vno0 = f->v[n] ; vno1 = f->v[m] ;
        f->v[n] = vno1 ; f->v[m] = vno0 ;
        mrisSetVertexFaceIndex(mris, vno0, fno) ;
        mrisSetVertexFaceIndex(mris, vno1, fno) ;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "reversing face %d orientation\n", fno) ;
        mrisNormalFace(mris, fno, n, norm) ;  /* how about vertex 2 ???? */
        dot = norm[0]*v->nx + norm[1]*v->ny + norm[2]*v->nz ;
      }
    }
  }

  MRIScomputeMetricProperties(mris) ;

  if (Gdiag_no >= 0)
  {
    v = &mris->vertices[Gdiag_no] ;
    for (dot = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
          
      /* only use already oriented (or non-defective) vertices */
      dot += v->nx*vn->nx + v->ny*vn->ny + v->nz*vn->nz ;
    }
    if (dot < 0)
      DiagBreak() ;
    for (dot = 0.0, n = 0 ; n < v->num ; n++)
    {
      fno = v->f[n] ; f = &mris->faces[fno] ;
          
      /* only use already oriented (or non-defective) vertices */
      dot += v->nx*f->nx + v->ny*f->ny + v->nz*f->nz ;
    }
    if (dot < 0)
      DiagBreak() ;
    MRISrestoreVertexPositions(mris, ORIG_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    v = &mris->vertices[Gdiag_no] ;
    for (dot = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
          
      /* only use already oriented (or non-defective) vertices */
      dot += v->nx*vn->nx + v->ny*vn->ny + v->nz*vn->nz ;
    }
    if (dot < 0)
      DiagBreak() ;
    
  }

#if 0
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
#endif
  return(oriented) ;
#endif
}
#endif
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int 
mrisComputeCanonicalBasis(MRI_SURFACE *mris, int fno, double origin[3], 
                          double e0[3], double e1[3])
{
  FACE    *f ;
  double   len, normal[3] ;
  float    fx, fy, fz ;

  f = &mris->faces[fno] ;
  mrisCalculateCanonicalFaceCentroid(mris, fno, &fx, &fy, &fz) ;
  origin[0] = (double)fx ; origin[1] = (double)fy ; origin[2] = (double)fz ;
  normal[0] = origin[0] ; normal[1] = origin[1] ; normal[2] = origin[2] ;
  len = 1.0f/VLEN(normal) ; SCALAR_MUL(normal, len, normal) ;

  /* pick any non-parallel vector and cross it with normal */
  e1[0] = normal[1] ; e1[1] = normal[2] ; e1[2] = normal[0] ; 
  CROSS(e0, normal, e1) ;
  if ((VZERO(e0)))  /* happened to pick parallel vector */
  {
    e1[0] = normal[1] ; e1[1] = -normal[2] ; e1[2] = normal[0] ; 
    CROSS(e0, normal, e1) ;
  }
  CROSS(e1, e0, normal) ;
  len = 1.0f/VLEN(e0) ; SCALAR_MUL(e0, len, e0) ;
  len = 1.0f/VLEN(e1) ; SCALAR_MUL(e1, len, e1) ;
  len = DOT(e0, e1) ;
  if ((VZERO(e0)) || (VZERO(e1)))
  {
    fprintf(stderr, "face %d, canonical basis degenerate!\n", fno) ;
  }
  if (fabs(len) > 0.001)
  {
    fprintf(stderr, "face %d, canonical basis not orthogonal!\n", fno) ;
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int 
mrisComputeCanonicalEdgeBasis(MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2,
                              double origin[3], double e0[3], double e1[3])
{
  VERTEX   *v0, *v1, *v2, *v3 ;
  double   len, normal[3] ;
  float    fx, fy, fz ;

  v0 = &mris->vertices[edge1->vno1] ;
  v1 = &mris->vertices[edge1->vno2] ;
  v2 = &mris->vertices[edge2->vno1] ;
  v3 = &mris->vertices[edge2->vno2] ;
  fx = (v0->cx + v1->cx + v2->cx + v3->cx) / 4 ; 
  fy = (v0->cy + v1->cy + v2->cy + v3->cy) / 4 ; 
  fz = (v0->cz + v1->cz + v2->cz + v3->cz) / 4 ;
  normal[0] = origin[0] = (double)fx ; 
  normal[1] = origin[1] = (double)fy ; 
  normal[2] = origin[2] = (double)fz ;
  len = 1.0f/VLEN(normal) ; SCALAR_MUL(normal, len, normal) ;

  /* pick any non-parallel vector and cross it with normal */
  e1[0] = normal[1] ; e1[1] = normal[2] ; e1[2] = normal[0] ; 
  CROSS(e0, normal, e1) ;
  if ((VZERO(e0)))  /* happened to pick parallel vector */
  {
    e1[0] = normal[1] ; e1[1] = -normal[2] ; e1[2] = normal[0] ; 
    CROSS(e0, normal, e1) ;
  }
  CROSS(e1, e0, normal) ;
  len = 1.0f/VLEN(e0) ; SCALAR_MUL(e0, len, e0) ;
  len = 1.0f/VLEN(e1) ; SCALAR_MUL(e1, len, e1) ;
  len = DOT(e0, e1) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        find the convex hull of the defect on the sphere (actually,
        just a disk, but at least it's convex....) 
------------------------------------------------------*/
static int
mrisFindDefectConvexHull(MRI_SURFACE *mris, DEFECT *defect)
{
  float  xmin, xmax, ymin, ymax, zmin, zmax ;
  VERTEX *v, *vn ;
  int    chull[200000], nfound, n, i, vno ;


  xmin = ymin = zmin = 100000 ;
  xmax = ymax = zmax = 0.0f ;


  /* now compute max radius on surface of sphere */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
    {
      vno = defect->vertices[i] ; 
      if (defect->status[i] == DISCARD_VERTEX)
        continue ;
    }
    else
      vno = defect->border[i-defect->nvertices] ; 
    v = &mris->vertices[vno] ;
    if (v->cx < xmin)
      xmin = v->cx ;
    if (v->cy < ymin)
      ymin = v->cy ;
    if (v->cz < zmin)
      zmin = v->cz ;
    if (v->cx > xmax)
      xmax = v->cx ;
    if (v->cy > ymax)
      ymax = v->cy ;
    if (v->cz > zmax)
      zmax = v->cz ;
  }
  defect->chull = chull ;
  defect->nchull = defect->nborder ;
  memcpy(chull, defect->border, defect->nborder*sizeof(int)) ;

  MRISclearMarks(mris) ;
  mrisMarkDefectConvexHull(mris, defect, 1) ;
  mrisMarkDefect(mris, defect, 1) ;

  do
  {
    nfound = 0 ;
    for (i = 0 ; i < defect->nchull ; i++)
    {
      v = &mris->vertices[defect->chull[i]] ;
      if (defect->chull[i] == Gdiag_no)
        DiagBreak() ;
      if ((v->cx >= xmin && v->cx <= xmax) &&
          (v->cy >= ymin && v->cy <= ymax) &&
          (v->cz >= zmin && v->cz <= zmax))
      {   /* vertex inside convex hull - add all its nbrs */
        for (n = 0 ; n < v->vnum ; n++)
        {
          vn = &mris->vertices[v->v[n]] ;
          if (vn->marked)   /* already in defect or convex hull */
            continue ;
          chull[defect->nchull+nfound++] = v->v[n] ;
          vn->marked = 1 ;
        }
      }
    }
    defect->nchull += nfound ;
  } while (nfound > 0) ;

  MRISclearMarks(mris) ;

  defect->chull = (int *)calloc(defect->nchull, sizeof(int)) ;
  if (!defect->chull)
    ErrorExit(ERROR_NO_MEMORY, 
              "mrisFindConvexHull: could not allocate %d vlist\n",
              defect->nchull) ;
  memcpy(defect->chull, chull, defect->nchull*sizeof(int)) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int 
mrisCheckSurface(MRI_SURFACE *mris)
{
  int     vno, n, nfaces, m, vno2, nbad, flist[100] ;
  VERTEX  *v ;

  /*  fprintf(stdout, "\n") ;*/
  for (nbad = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vno2 = v->v[n] ;
      if (vno2 < vno)
        continue ;
      for (nfaces = m = 0 ; m < v->vnum ; m++)
      {
        if (v->v[m] == vno2)
          continue ;
        if (vertexNeighbor(mris, vno2, v->v[m]) && 
            isFace(mris,vno, vno2, v->v[m]))
        {
          flist[nfaces] = findFace(mris, vno, vno2, v->v[m]) ;
          nfaces++ ;
        }
      }
      if (nfaces != 2)
      {
        int i ;

        nbad++ ;
        fprintf(stdout, "edge %d <--> %d has %d face(s)! ",
                vno, vno2, nfaces) ;
        fprintf(stdout, "(") ;
        for (i = 0 ; i < nfaces ; i++)
          fprintf(stdout, "%d%s", flist[i], i < nfaces-1 ? "," : "") ;
        fprintf(stdout, ")\n") ;
        mrisDumpDefectiveEdge(mris, vno, vno2) ;
        DiagBreak() ;
      }
    }
  }
  fprintf(stdout, "%d defective edges\n", nbad) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int 
mrisMarkBadEdgeVertices(MRI_SURFACE *mris, int mark)
{
  int     vno, n, nfaces, m, vno2, nmarked ;
  VERTEX  *v ;

  for (nmarked = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vno2 = v->v[n] ;
      if (vno2 < vno)
        continue ;
      for (nfaces = m = 0 ; m < v->vnum ; m++)
      {
        if (v->v[m] == vno2)
          continue ;
        if (vertexNeighbor(mris, vno2, v->v[m]) && 
            isFace(mris,vno, vno2, v->v[m]))
          nfaces++ ;
      }
      if (nfaces != 2)
      {
        v->marked = mark ;
        mris->vertices[vno2].marked = mark ;
        nmarked += 2 ;
        break ;
      }
    }
  }
  return(nmarked) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisDilateAmbiguousVertices(MRI_SURFACE *mris, int mark, int ndil)
{
  int    vno, i, n ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < ndil ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->marked != mark)
        continue ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked == mark) 
          continue ;
        vn->marked = mark+1 ;
      }
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->marked == mark+1)
        v->marked = mark ;
    }
  }
  return(NO_ERROR) ;
}

#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisDumpDefectiveEdge(MRI_SURFACE *mris, int vno1, int vno2)
{
#if 0
  FILE   *fp ;
  char   fname[STRLEN] ;
  int    n, m, fno, first = 1 ;
  VERTEX *v1, *v2, *vn ;
  double origin[3], e0[3], e1[3], cx, cy, cz, x, y ;
  FACE   *f ;

  sprintf(fname, "edge%d_%d.log", vno1, vno2) ;
  fp = fopen(fname, "w") ;
  

  v1 = &mris->vertices[vno1] ; v2 = &mris->vertices[vno2] ;
  for (n = 0 ; n < v1->vnum ; n++)
  {
    if (v1->v[n] == vno2)
      continue ;
    fno = findFace(mris, vno1, vno2, v1->v[n]) ;
    if ((fno >= 0) && vertexNeighbor(mris, vno2, v1->v[n]))
    {
      f = &mris->faces[fno] ;
      if (first)
      {
        mrisComputeCanonicalBasis(mris, fno, origin, e0, e1) ;
        first = 0 ;
      }
      fprintf(fp, "# triangle %d\n", fno) ;
      for (m = 0 ; m < VERTICES_PER_FACE ; m++)
      {
        vn = &mris->vertices[f->v[m]] ;
        cx = vn->cx-origin[0]; cy = vn->cy-origin[1]; cz = vn->cz-origin[2];
        x = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
        y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
        fprintf(fp, "# vertex %d\n", f->v[m]) ;
        fprintf(fp, "%f %f\n", x, y) ;
      }
      vn = &mris->vertices[f->v[0]] ;
      cx = vn->cx-origin[0]; cy = vn->cy-origin[1]; cz = vn->cz-origin[2];
      x = cx*e0[0] + cy*e0[1] + cz*e0[2] ;
      y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
        fprintf(fp, "#%d\n", f->v[0]) ;
      fprintf(fp, "%f %f\n", x, y) ;
      fprintf(fp, "\n") ;
    }
  }
  cx = v1->cx-origin[0]; cy = v1->cy-origin[1]; cz = v1->cz-origin[2];
  x = cx*e0[0] + cy*e0[1] + cz*e0[2] ; y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  fprintf(fp, "%f %f\n", x, y) ;
  cx = v2->cx-origin[0]; cy = v2->cy-origin[1]; cz = v2->cz-origin[2];
  x = cx*e0[0] + cy*e0[1] + cz*e0[2] ; y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  fprintf(fp, "%f %f\n", x, y) ;
  fclose(fp) ;
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisDumpTriangle(MRI_SURFACE *mris, int fno)
{
  char   fname[STRLEN] ; 
  VERTEX *v0, *v1, *v2 ;
  FACE   *f ;
  FILE   *fp ;
  double cx, cy, cz, x, y, origin[3], e0[3], e1[3] ;

  mrisComputeCanonicalBasis(mris, fno, origin, e0, e1) ;
  f = &mris->faces[fno] ;
  sprintf(fname, "triangle%d.log", fno) ;
  fp = fopen(fname, "w") ;

  v0 = &mris->vertices[f->v[0]] ;
  v1 = &mris->vertices[f->v[1]] ;
  v2 = &mris->vertices[f->v[2]] ;
  fprintf(fp, "# triangle %d, vertices %d, %d, %d\n",
          fno, f->v[0], f->v[1], f->v[2]) ;

  cx = v0->cx-origin[0]; cy = v0->cy-origin[1]; cz = v0->cz-origin[2];
  x = cx*e0[0] + cy*e0[1] + cz*e0[2] ; y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  fprintf(fp, "%f  %f\n", x, y) ;
  cx = v1->cx-origin[0]; cy = v1->cy-origin[1]; cz = v1->cz-origin[2];
  x = cx*e0[0] + cy*e0[1] + cz*e0[2] ; y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  fprintf(fp, "%f  %f\n", x, y) ;
  cx = v2->cx-origin[0]; cy = v2->cy-origin[1]; cz = v2->cz-origin[2];
  x = cx*e0[0] + cy*e0[1] + cz*e0[2] ; y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  fprintf(fp, "%f  %f\n", x, y) ;
  cx = v0->cx-origin[0]; cy = v0->cy-origin[1]; cz = v0->cz-origin[2];
  x = cx*e0[0] + cy*e0[1] + cz*e0[2] ; y = cx*e1[0] + cy*e1[1] + cz*e1[2] ;
  fprintf(fp, "%f  %f\n", x, y) ;

  fclose(fp) ;
  return(NO_ERROR) ;
}
#endif
static int
mrisDefectRemoveNegativeVertices(MRI_SURFACE *mris, DEFECT *defect)
{
  int    i, n ;
  VERTEX *v ;

  for (i = 0 ; i < defect->nvertices ; i++)
  {
    if (defect->status[i] == DISCARD_VERTEX)
      continue ;
    v = &mris->vertices[defect->vertices[i]] ;
    for (n = 0 ; n < v->num ; n++)
      if (mris->faces[v->f[n]].area < 0.0)
        defect->status[i] = DISCARD_VERTEX ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisDefectRemoveDegenerateVertices(MRI_SURFACE *mris, float min_sphere_dist,
                                   DEFECT *defect)
{
  float  dx, dy, dz, dist ;
  int    i, j ;
  VERTEX *v, *vn ;

  /* discard vertices that are too close to another vertex on sphere */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
    {
      if (defect->status[i] == DISCARD_VERTEX)
        continue ;
      v = &mris->vertices[defect->vertices[i]] ;
    }
    else
      v = &mris->vertices[defect->border[i-defect->nvertices]] ;
    for (j = i+1 ; j < defect->nvertices ; j++)
    {
      if (defect->status[j] == DISCARD_VERTEX)
        continue ;
      vn = &mris->vertices[defect->vertices[j]] ;
      dx = vn->cx-v->cx ;
      dy = vn->cy-v->cy ;
      dz = vn->cz-v->cz ;
      dist = (dx*dx+dy*dy+dz*dz) ;  /* no sqrt */
      if (dist < min_sphere_dist)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "discarding proximal vertex %d\n",
                  defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX ;
      }
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisDefectRemoveProximalVertices(MRI_SURFACE *mris, float min_orig_dist,
                                 DEFECT *defect)
{
  float  dx, dy, dz, dist ;
  int    i, j ;
  VERTEX *v, *vn ;

  /* discard vertices that are too close to another vertex on sphere */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
    {
      if (defect->status[i] == DISCARD_VERTEX)
        continue ;
      v = &mris->vertices[defect->vertices[i]] ;
    }
    else
      v = &mris->vertices[defect->border[i-defect->nvertices]] ;
    for (j = i+1 ; j < defect->nvertices ; j++)
    {
      if (defect->status[j] == DISCARD_VERTEX)
        continue ;
      vn = &mris->vertices[defect->vertices[j]] ;
      dx = vn->origx-v->origx ;
      dy = vn->origy-v->origy ;
      dz = vn->origz-v->origz ;
      dist = (dx*dx+dy*dy+dz*dz) ;  /* no sqrt */
      if (dist < min_orig_dist)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "discarding proximal vertex %d\n",
                  defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX ;
      }
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisInitializeNeighborhood(MRI_SURFACE *mris, int vno)
{
  VERTEX  *v, *vnb, *vnb2 ;
  int     vtmp[MAX_NEIGHBORS], vnum, i, j, n, neighbors, nsize ;

  v = &mris->vertices[vno] ;
  if (vno == Gdiag_no)
    DiagBreak()  ;

  v->nsize = mris->nsize ;
  if (v->ripflag || !v->vnum)
    return(ERROR_BADPARM) ;
  memmove(vtmp, v->v, v->vnum*sizeof(int)) ;

  /* mark 1-neighbors so we don't count them twice */
  v->marked = 1 ;

  vnum = neighbors = v->vnum ;
  for (nsize = 2 ; nsize <= v->nsize ; nsize++)
  {
    /* mark all current neighbors */
    vnum = neighbors ;  /* neighbors will be incremented during loop */
    for (i = 0 ; i < neighbors ; i++)
      mris->vertices[vtmp[i]].marked = 1 ;
    for (i = 0; neighbors < MAX_NEIGHBORS && i < vnum; i++)
    {
      n = vtmp[i] ;
      vnb = &mris->vertices[n] ;
      if (vnb->ripflag)
        continue ;
    
      for (j = 0 ; j < vnb->vnum ; j++)
      {
        vnb2 = &mris->vertices[vnb->v[j]] ;
        if (vnb2->ripflag || vnb2->marked)
          continue ;
        vtmp[neighbors] = vnb->v[j] ;
        vnb2->marked = 1 ;
        if (++neighbors >= MAX_NEIGHBORS)
        {
          fprintf(stderr, "vertex %d has too many neighbors!\n",vno) ;
          break ;
        }
      }
    }
  }
  /*
    now reallocate the v->v structure and place the 2-connected neighbors
    suquentially after the 1-connected neighbors.
  */
  free(v->v) ;
  v->v = (int *)calloc(neighbors, sizeof(int)) ;
  if (!v->v)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISsetNeighborhoodSize: could not allocate list of %d "
              "nbrs at v=%d", neighbors, vno) ;
  
  v->marked = 0 ;
  for (n = 0 ; n < neighbors ; n++)
  {
    v->v[n] = vtmp[n] ;
    mris->vertices[vtmp[n]].marked = 0 ;
  }
  if (v->dist)
    free(v->dist) ;
  if (v->dist_orig)
    free(v->dist_orig) ;
  
  v->dist = (float *)calloc(neighbors, sizeof(float)) ;
  if (!v->dist)
    ErrorExit(ERROR_NOMEMORY,
              "MRISsetNeighborhoodSize: could not allocate list of %d "
              "dists at v=%d", neighbors, vno) ;
  v->dist_orig = (float *)calloc(neighbors, sizeof(float)) ;
  if (!v->dist_orig)
    ErrorExit(ERROR_NOMEMORY,
              "MRISsetNeighborhoodSize: could not allocate list of %d "
              "dists at v=%d", neighbors, vno) ;
  switch (v->nsize)
  {
  case 2:
    v->v2num = neighbors ;
    break ;
  case 3:
    v->v3num = neighbors ;
    break ;
  default:   /* store old neighborhood size in v3num */
    v->v3num = v->vtotal ;
    break ;
  }
  v->vtotal = neighbors ;
  for (n = 0 ; n < neighbors ; n++)
    for (i = 0 ; i < neighbors ; i++)
      if (i != n && v->v[i] == v->v[n])
        fprintf(stderr, 
                "warning: vertex %d has duplicate neighbors %d and %d!\n",
                vno, i, n) ;
  if ((vno == Gdiag_no) && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
  {
    fprintf(stdout, "v %d: vnum=%d, v2num=%d, vtotal=%d\n",
            vno, v->vnum, v->v2num, v->vtotal) ;
    for (n = 0 ; n < neighbors ; n++)
      fprintf(stdout, "v[%d] = %d\n", n, v->v[n]) ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmarkNegativeVertices(MRI_SURFACE *mris, int mark) 
{
  int    fno, n ;
  FACE   *f ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->area < 0)
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        mris->vertices[f->v[n]].marked = mark ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISripNegativeVertices(MRI_SURFACE *mris) 
{
  int    fno, n ;
  FACE   *f ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->area < 0)
    {
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        mris->vertices[f->v[n]].ripflag = 1 ;
      f->ripflag = 1 ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRIScomputeAverageCurvature(MRI_SURFACE *mris, double *psigma)
{
  double mean, var, total, total_sq, nv, d;
  int    vno ;
  VERTEX *v ;

  for (vno = 0, nv = total_sq = total = 0.0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    d = (double)v->curv ;
    total_sq += d*d ; total += d ; nv += 1.0 ;
  }
  if (nv)
  {
    mean = total / nv ;
    var = total_sq / nv - (mean*mean) ;
  }
  else
    var = mean = 0.0 ;
  if (psigma)
    *psigma = sqrt(var) ;
  return(mean) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int          
MRIScopyValToVal2(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val2 = v->val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyValuesToVal2Bak(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val2bak = v->val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyValToValBak(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->valbak = v->val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIScopyValToVal2Bak(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val2bak = v->val ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------

        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISsqrtVal(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->val > 0)
      v->val = sqrt(v->val) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------

        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISmulVal(MRI_SURFACE *mris, float mul)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val *= mul ;
  }
  return(NO_ERROR) ;
}


SMALL_SURFACE *
MRISreadVerticesOnly(char *fname)
{
  SMALL_SURFACE *mriss = NULL ;
  int           type, magic, version, ix, iy, iz, nquads, nvertices, vno ;
  SMALL_VERTEX  *vertex ;
  FILE          *fp ;

  type = MRISfileNameType(fname) ;
  switch (type)
  {
  case MRIS_ASCII_TRIANGLE_FILE:
  case MRIS_ICO_FILE:
  case MRIS_GEO_TRIANGLE_FILE:
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, 
                 "MRISreadVerticesOnly: file type %d not supported",type)) ;
    break ;  /* not used */
  default:
    break ;
  }
    
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL,(ERROR_NOFILE,"MRISread(%s): could not open file",
                      fname));
  
  fread3(&magic, fp) ;
  if (magic == QUAD_FILE_MAGIC_NUMBER) 
  {
    version = -1;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "new surface file format\n");
  }
  else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) 
  {
    version = -2 ;
  }
  else if (magic == TRIANGLE_FILE_MAGIC_NUMBER)
  {
    fclose(fp) ;
    mriss = mrisReadTriangleFileVertexPositionsOnly(fname) ;
    if (!mriss)
      ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n")) ;
    version = -3 ;
  }
  else
  {
    rewind(fp);
    version = 0;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      printf("surfer: old surface file format\n");
  }
  if (version >= -2)  /* some type of quadrangle file */
  {
    fread3(&nvertices, fp);
    fread3(&nquads, fp);   /* # of qaudrangles - not triangles */
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout,"reading %d vertices and %d faces.\n",
      nvertices,2*nquads);
    
    mriss = (SMALL_SURFACE *)calloc(1, sizeof(SMALL_SURFACE)) ;
    if (!mriss)
      ErrorReturn(NULL, 
                  (ERROR_NOMEMORY, 
                   "MRISreadVerticesOnly: could not allocate surface")) ;
  
    mriss->nvertices = nvertices ;
    mriss->vertices = (SMALL_VERTEX *)calloc(nvertices, sizeof(SMALL_VERTEX)) ;
    if (!mriss->nvertices)
    {
      free(mriss) ;
      ErrorReturn(NULL, 
                  (ERROR_NOMEMORY, 
                   "MRISreadVerticesOnly: could not allocate surface")) ;
    }
    for (vno = 0 ; vno < nvertices ; vno++)
    {
      vertex = &mriss->vertices[vno] ;
      if (version == -1)
      {
        fread2(&ix,fp);
        fread2(&iy,fp);
        fread2(&iz,fp);
        vertex->x = ix/100.0;
        vertex->y = iy/100.0;
        vertex->z = iz/100.0;
      }
      else  /* version == -2 */
      {
        vertex->x = freadFloat(fp) ;
        vertex->y = freadFloat(fp) ;
        vertex->z = freadFloat(fp) ;
      }
    }
  }

  return(mriss) ;
}
int
MRISSfree(SMALL_SURFACE **pmriss)
{
  SMALL_SURFACE *mriss ;

  mriss = *pmriss ;
  *pmriss = NULL ;
  free(mriss->vertices) ;
  free(mriss) ;
  return(NO_ERROR) ;
}
int
MRISextractCurvatureVector(MRI_SURFACE *mris, float *curvs)
{
  int     vno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    curvs[vno] = mris->vertices[vno].curv ;

  return(NO_ERROR) ;
}
int
MRISextractCurvatureDoubleVector(MRI_SURFACE *mris, double *curvs)
{
  int     vno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    curvs[vno] = (double)mris->vertices[vno].curv ;

  return(NO_ERROR) ;
}
int
MRISimportCurvatureVector(MRI_SURFACE *mris, float *curvs)
{
  int     vno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].curv = curvs[vno] ;

  return(NO_ERROR) ;
}

int
MRISimportValVector(MRI_SURFACE *mris, float *vals)
{
  int     vno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].val = vals[vno] ;

  return(NO_ERROR) ;
}
int
MRISexportValVector(MRI_SURFACE *mris, float *vals)
{
  int     vno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    vals[vno] = mris->vertices[vno].val ;

  return(NO_ERROR) ;
}
int
MRISmaskLabel(MRI_SURFACE *mris, LABEL *area)
{
  int     i ;
  VERTEX  *v ;

  for (i = 0 ; i < area->n_points ; i++)
  {
    v = &mris->vertices[area->lv[i].vno] ;
    v->curv = v->stat = v->val = v->imag_val=v->val2=v->valbak=v->val2bak = 0.0;
  }
  return(NO_ERROR) ;
}
int
MRISmaskNotLabel(MRI_SURFACE *mris, LABEL *area)
{
  int     i, vno ;
  VERTEX  *v ;

  for (i = 0 ; i < area->n_points ; i++)
  {
    v = &mris->vertices[area->lv[i].vno] ;
    v->marked = 1 ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked)
      continue ;
    v->curv = v->stat = v->val = v->imag_val=v->val2=v->valbak=v->val2bak = 0.0;
  }
  for (i = 0 ; i < area->n_points ; i++)
  {
    v = &mris->vertices[area->lv[i].vno] ;
    v->marked = 0 ;
  }
  return(NO_ERROR) ;
}
int
MRISripLabel(MRI_SURFACE *mris, LABEL *area)
{
  int     i ;
  VERTEX  *v ;

  for (i = 0 ; i < area->n_points ; i++)
  {
    v = &mris->vertices[area->lv[i].vno] ;
    v->ripflag = 1 ;
  }
  return(NO_ERROR) ;
}
int
MRISripNotLabel(MRI_SURFACE *mris, LABEL *area)
{
  int     i, vno ;
  VERTEX  *v ;

  for (i = 0 ; i < area->n_points ; i++)
  {
    v = &mris->vertices[area->lv[i].vno] ;
    v->marked = 1 ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked)
      continue ;
    v->ripflag = 1 ;
  }
  for (i = 0 ; i < area->n_points ; i++)
  {
    v = &mris->vertices[area->lv[i].vno] ;
    v->marked = 0 ;
  }
  return(NO_ERROR) ;
}
int
MRISsegmentMarked(MRI_SURFACE *mris, LABEL ***plabel_array, int *pnlabels,
                  float min_label_area)
{
  int     vno, nfound, n, nlabels, *marks ;
  VERTEX  *v ;
  LABEL   *area = NULL, **tmp, **label_array ;

  marks = (int *)calloc(mris->nvertices, sizeof(int)) ;
  label_array = (LABEL **)calloc(mris->nvertices, sizeof(LABEL *)) ;
  if (!label_array || !marks)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: MRISsegmentMarked could not allocate tmp storage",
              Progname) ;

  /* save current marks */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    marks[vno] = v->marked ;
    if (v->marked != 0)
      v->marked = 1 ;
  }

  nlabels = 0 ;
  do
  {
    nfound = 0 ;

    v = &mris->vertices[0] ;

    /* find a marked vertex */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked != 1)
        continue ;
      break ;
    }
    if (vno < mris->nvertices)
    {
      area = LabelAlloc(mris->nvertices, NULL, NULL) ;
      area->n_points = 1 ;
      area->lv[0].x = v->x ; area->lv[0].y = v->y ;area->lv[0].z = v->z ;
      area->lv[0].vno = vno ;
      LabelFillMarked(area, mris) ;
      if (LabelArea(area, mris) >= min_label_area)
        label_array[nlabels++] = LabelCopy(area, NULL) ;
      LabelFree(&area) ;
      nfound = 1 ;
    }
    else
      nfound = 0 ;

  } while (nfound > 0) ;

  /* restore original marks */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->marked = marks[vno] ;
  }

  free(marks) ;

  /* crunch label array down to a reasonable size */
  tmp = label_array ;
  label_array = (LABEL **)calloc(mris->nvertices, sizeof(LABEL *)) ;
  if (!label_array)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: MRISsegmentMarked could not allocate tmp storage",
              Progname) ;
  for (n = 0 ; n < nlabels ; n++)
    label_array[n] = tmp[n] ;
  free(tmp) ;
  *plabel_array = label_array ;
  *pnlabels = nlabels ;
  return(NO_ERROR) ;
}

int
MRISsegmentAnnotated(MRI_SURFACE *mris, LABEL ***plabel_array, int *pnlabels,
                  float min_label_area)
{
  int     vno, nfound, n, nlabels ;
  VERTEX  *v ;
  LABEL   *area = NULL, **tmp, **label_array ;

  label_array = (LABEL **)calloc(mris->nvertices, sizeof(LABEL *)) ;
  if (!label_array)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: MRISsegmentMarked could not allocate tmp storage",
              Progname) ;

  MRISclearMarks(mris) ;

  nlabels = 0 ;
  do
  {
    nfound = 0 ;

    v = &mris->vertices[0] ;

    /* find an un-marked vertex */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->annotation == 0 || v->marked)
        continue ;
      break ;
    }
    if (vno < mris->nvertices)
    {
      area = LabelAlloc(mris->nvertices, NULL, NULL) ;
      area->n_points = 1 ;
      area->lv[0].x = v->x ; area->lv[0].y = v->y ;area->lv[0].z = v->z ;
      area->lv[0].vno = vno ;
      LabelFillAnnotated(area, mris) ;
      if (LabelArea(area, mris) >= min_label_area)
      {
        label_array[nlabels++] = LabelCopy(area, NULL) ;
      }
      LabelMarkSurface(area, mris) ;
      LabelFree(&area) ;
      nfound = 1 ;
    }
    else
      nfound = 0 ;

  } while (nfound > 0) ;


  /* crunch label array down to a reasonable size */
  tmp = label_array ;
  label_array = (LABEL **)calloc(mris->nvertices, sizeof(LABEL *)) ;
  if (!label_array)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: MRISsegmentAnnotated could not allocate tmp storage",
              Progname) ;
  for (n = 0 ; n < nlabels ; n++)
    label_array[n] = tmp[n] ;
  free(tmp) ;
  *plabel_array = label_array ;
  *pnlabels = nlabels ;
  return(NO_ERROR) ;
}

int
MRISsubsampleDist(MRI_SURFACE *mris, float spacing)
{
  int k,m,n, sub_num;
  VERTEX *v;


  sub_num = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->d = 10000;
    v->val = 0;
  }
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    for (m=0;m<v->vnum;m++)
    {
      if (mris->vertices[v->v[m]].d+1 < v->d)
        v->d = mris->vertices[v->v[m]].d+1;
    }
    if (v->d>=spacing)
    {
      v->d = 0;
      v->val = 1;
      sub_num++;
    }
    for (m=0;m<v->vnum;m++)
    {
      if (mris->vertices[v->v[m]].d > v->d+1)
        mris->vertices[v->v[m]].d = v->d+1;
    }
  }
  for (k=mris->nvertices-1;k>=0;k--)
  {
    v = &mris->vertices[k];
    for (m=0;m<v->vnum;m++)
    {
      if (mris->vertices[v->v[m]].d+1 < v->d)
        v->d = mris->vertices[v->v[m]].d+1;
      if (mris->vertices[v->v[m]].d > v->d+1)
        mris->vertices[v->v[m]].d = v->d+1;
    }
  }

  if (spacing==2)
  for (k=0;k<mris->nvertices;k++)
  if (mris->vertices[k].d > 0)
  {
    v = &mris->vertices[k];
    n = 0;
    for (m=0;m<v->vnum;m++)
    {
      if (mris->vertices[v->v[m]].d == 0)
        n++;
    }
    if (n <= 2)
    {
      v->d = 0;
      v->val = 1;
      v->fixedval = TRUE;
      sub_num++;
    }
    for (m=0;m<v->vnum;m++)
    {
      if (mris->vertices[v->v[m]].d > v->d+1)
        mris->vertices[v->v[m]].d = v->d+1;
    }
  }

  return(sub_num) ;
}
int
MRISwriteDecimation(MRI_SURFACE *mris, char *fname)
{
  int k;
  FILE *fptr;
  
  fptr = fopen(fname,"w");
  if (fptr==NULL)
    ErrorReturn(ERROR_BADFILE, 
                (ERROR_BADFILE, "MRISwriteDecimation: could not create %s",
                 fname)) ;
  fputc('\0',fptr);
  fwriteInt(mris->nvertices,fptr);
  for (k=0;k<mris->nvertices;k++)
  {
    if (mris->vertices[k].d==0)
      fputc('\1',fptr);
    else
      fputc('\0',fptr);
  }
  fclose(fptr);
  return(NO_ERROR) ;
}
int
MRISreadDecimation(MRI_SURFACE *mris, char *fname)
{
  int k,d, ndec;
  char c;
  FILE *fptr;


  ndec = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].undefval = TRUE;
    mris->vertices[k].fixedval = FALSE;
  }
  fptr = fopen(fname,"r");
  if (fptr==NULL) 
    ErrorReturn(ERROR_BADFILE, 
                (ERROR_BADFILE, "MRISreadDecimation: could not create %s",
                 fname)) ;
  c = fgetc(fptr);
  if (c=='#')
  {
    fscanf(fptr,"%*s");
    fscanf(fptr,"%d",&d);
    if (d!=mris->nvertices) 
      ErrorReturn(0,
                  (ERROR_BADFILE,
                   "%s: decimation file %s has wrong # of vertices\n",
                   Progname, fname, d)) ;
    for (k=0;k<mris->nvertices;k++)
    {
      fscanf(fptr,"%d",&d);
      if (d!=0)
      {
        mris->vertices[k].d=0;
        mris->vertices[k].fixedval=TRUE;
        mris->vertices[k].undefval=FALSE;
        ndec++;
      }
    }
  }
  else
  {
    d = freadInt(fptr);
    if (d!=mris->nvertices) 
      ErrorReturn(0,
                  (ERROR_BADFILE,
                   "%s: decimation file %s has wrong # of vertices\n",
                   Progname, fname, d)) ;
    for (k=0;k<mris->nvertices;k++)
    {
      c = fgetc(fptr);
      if (c!='\0')
      {
        mris->vertices[k].d=0;
        mris->vertices[k].fixedval=TRUE;
        mris->vertices[k].undefval=FALSE;
        ndec++;
      }
    }
  }
  fclose(fptr);
  return(ndec) ;
}

int
MRIScombine(MRI_SURFACE *mris_src, MRI_SURFACE *mris_total, 
            MRIS_HASH_TABLE *mht, int which)
{
  int    vno ;
  VERTEX *v, *vdst ;
  MHT    *mht_src = NULL ;
  double max_len, mean ;

  MRISclearMarks(mris_total) ;
  for (vno = 0 ; vno < mris_total->nvertices ; vno++)
  {
    vdst = &mris_total->vertices[vno] ;
    vdst->d = 0 ;
  }

  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    v = &mris_src->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    vdst = MHTfindClosestVertex(mht, mris_total, v) ;
    if (!vdst)
    {
      ErrorPrintf(ERROR_BADPARM, "MRIScombine: cannot map vno %d", vno) ;
      continue ;
    }
    if (vdst-mris_total->vertices == Gdiag_no)
      DiagBreak() ;
    vdst->marked++ ;
    switch (which)
    {
		case VERTEX_COORDS:
			vdst->origx += v->origx ; vdst->origy += v->origy ; vdst->origz += v->origz ; 
			break ;
    case VERTEX_AREA:
      vdst->d += v->origarea ;
      break ;
    case VERTEX_CURV:
      vdst->d += v->curv ;
      break ;
    case VERTEX_VALS:
      vdst->d += v->val ;
      break ;
		case VERTEX_ANNOTATION:
			vdst->annotation = v->annotation ;
			break ;
    }
  }

  /* normalize by # of times a vertex is mapped to */
  for (vno = 0 ; vno < mris_total->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    vdst = &mris_total->vertices[vno] ;
    if (vdst->ripflag || !vdst->marked)
      continue ;
    mean = vdst->d / (float)vdst->marked ;
    switch (which)
    {
		case VERTEX_COORDS:
			vdst->origx /= (float)vdst->marked ; 
			vdst->origy /= (float)vdst->marked ; 
			vdst->origz /= (float)vdst->marked ;
			break ;
    case VERTEX_AREA:  /* don't normalize by # of vertices mapped!! */
      vdst->origarea += vdst->d ;
      vdst->val2 += vdst->d * vdst->d ;
      break ;
    case VERTEX_CURV:
      vdst->curv += mean ;
      vdst->val2 += mean*mean ;
      break ;
    case VERTEX_VALS:
      vdst->val += mean ;
      vdst->val2 += mean*mean ;
      break ;
    }
  }

  /* sample from dst to source to fill holes */
  for (vno = 0 ; vno < mris_total->nvertices ; vno++)
  {
    vdst = &mris_total->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (vdst->marked || vdst->ripflag)
      continue ;
    if (!mht_src)
    {
      MRIScomputeVertexSpacingStats(mris_src, NULL, NULL, &max_len, NULL,NULL);
      mht_src = 
        MHTfillVertexTableRes(mris_src, NULL, CURRENT_VERTICES, 2*max_len);
    }
    v = MHTfindClosestVertex(mht_src, mris_src, vdst) ;
    if (!v)
    {
      ErrorPrintf(ERROR_BADPARM, "MRIScombine: cannot map dst vno %d", vno) ;
      continue ;
    }
    if (v - mris_src->vertices == Gdiag_no)
      DiagBreak() ;
    vdst->marked++ ;
    switch (which)
    {
		case VERTEX_COORDS:
			vdst->origx = v->origx ; vdst->origy = v->origy ; vdst->origz = v->origz ; 
			break ;
		case VERTEX_ANNOTATION:
			vdst->annotation = v->annotation ;
			break ;
    case VERTEX_AREA:
      vdst->origarea += v->origarea ;
      vdst->val2 += (v->origarea*v->origarea) ;
      break ;
    case VERTEX_CURV:
      vdst->curv += v->curv ;
      vdst->val2 += (v->curv*v->curv) ;
      break ;
    case VERTEX_VALS:
      vdst->val += v->val ;
      vdst->val2 += (v->val*v->val) ;
      break ;
    }
  }
  if (mht_src) 
    MHTfree(&mht_src) ;

  return(NO_ERROR) ;
}

#if 1
int
MRISsphericalCopy(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, 
            MRIS_HASH_TABLE *mht, int which)
{
  int    vno ;
  VERTEX *v, *vdst ;
  MHT    *mht_src = NULL ;
  double max_len ;

  MRISclear(mris_dst, which) ;
  for (vno = 0 ; vno < mris_dst->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    vdst->d = 0 ; vdst->val2 = 0 ;
  }


  MRIScomputeVertexSpacingStats(mris_src, NULL, NULL, &max_len, NULL,NULL);
  mht_src = MHTfillVertexTableRes(mris_src, NULL, CURRENT_VERTICES, 2*max_len);

  /* sample from dst to source to fill holes */
  for (vno = 0 ; vno < mris_dst->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->ripflag)
      continue ;
    v = MHTfindClosestVertex(mht_src, mris_src, vdst) ;
    if (!v)
    {
      ErrorPrintf(ERROR_BADPARM, 
                  "MRISsphericalCopy: cannot map dst vno %d", vno) ;
      continue ;
    }
    if (v-mris_src->vertices == Gdiag_no)
      DiagBreak() ;
    vdst->marked++ ;
    vdst->val2 = v->val2 ;
    switch (which)
    {
		case VERTEX_COORDS:
			vdst->origx = v->origx ; vdst->origy = v->origy ; vdst->origz = v->origz ; 
			break ;
		case VERTEX_ANNOTATION:
			vdst->annotation = v->annotation ;
			break ;
    case VERTEX_AREA:
      vdst->origarea = v->origarea ;
      break ;
    case VERTEX_CURV:
      vdst->curv = v->curv ;
      break ;
    case VERTEX_VALS:
      vdst->val = v->val ;
      break ;
    }
  }
  if (mht_src) 
    MHTfree(&mht_src) ;

  return(NO_ERROR) ;
}
#else
int
MRISsphericalCopy(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, 
            MRIS_HASH_TABLE *mht, int which)
{
  int    vno ;
  VERTEX *v, *vdst ;
  MHT    *mht_src = NULL ;
  double max_len, mean ;

  MRISclearMarks(mris_dst) ; MRISclear(mris_dst, which) ;
  MRISclearMarks(mris_src) ;
  for (vno = 0 ; vno < mris_dst->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    vdst->d = 0 ; vdst->val2 = 0 ;
  }

  /*
    First determine how much of a fan in there is in the mapping.
  */

  /*
    go through each vertex in the source and see what destination
    vertex it maps to.
  */
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    v = &mris_src->vertices[vno] ;
    if (v->ripflag)
      continue ;
    vdst = MHTfindClosestVertex(mht, mris_dst, v) ;
    if (!vdst)
    {
      ErrorPrintf(ERROR_BADPARM, "MRISsphericalCopy: cannot map vno %d", vno) ;
      continue ;
    }
    if (vdst-mris_dst->vertices == Gdiag_no)
      DiagBreak() ;
    vdst->marked++ ; v->marked++ ;
  }

  /* sample from dst to source  */
  for (vno = 0 ; vno < mris_dst->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->marked || vdst->ripflag)
      continue ;
    if (!mht_src)
    {
      MRIScomputeVertexSpacingStats(mris_src, NULL, NULL, &max_len, NULL,NULL);
      mht_src = 
        MHTfillVertexTableRes(mris_src, NULL, CURRENT_VERTICES, 2*max_len);
    }
    v = MHTfindClosestVertex(mht_src, mris_src, vdst) ;
    if (!v)
    {
      ErrorPrintf(ERROR_BADPARM, 
                  "MRISsphericalCopy: cannot map dst vno %d", vno) ;
      continue ;
    }
    if (v-mris_src->vertices == Gdiag_no)
      DiagBreak() ;
    vdst->marked++ ; v->marked++ ;
  }

  MRISclearMarks(mris_dst) ;
  /*
    go through each vertex in the source and sample it onto
    the destination surface.
  */
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    v = &mris_src->vertices[vno] ;
    if (v->ripflag)
      continue ;
    vdst = MHTfindClosestVertex(mht, mris_dst, v) ;
    if (!vdst)
    {
      ErrorPrintf(ERROR_BADPARM, "MRISsphericalCopy: cannot map vno %d", vno) ;
      continue ;
    }
    if (vdst-mris_dst->vertices == Gdiag_no)
      DiagBreak() ;
    vdst->marked++ ;
    vdst->val2 += v->val2/(float)v->marked ;  /* variances */
    switch (which)
    {
    case VERTEX_AREA:
      vdst->d += v->origarea/(float)v->marked ;
      break ;
    case VERTEX_CURV:
      vdst->d += v->curv ;
      break ;
    case VERTEX_VALS:
      vdst->d += v->val ;
      break ;
    }
  }

  /* normalize by # of times a vertex is mapped to */
  for (vno = 0 ; vno < mris_dst->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->ripflag || !vdst->marked)
      continue ;
    vdst->val2 /= (float)vdst->marked ;
    mean = vdst->d / (float)vdst->marked ;
    switch (which)
    {
    case VERTEX_AREA:
      vdst->origarea = mean ;
      break ;
    case VERTEX_CURV:
      vdst->curv = mean ;
      break ;
    case VERTEX_VALS:
      vdst->val = mean ;
      break ;
    }
  }

  /* sample from dst to source to fill holes */
  for (vno = 0 ; vno < mris_dst->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->marked || vdst->ripflag)
      continue ;
    v = MHTfindClosestVertex(mht_src, mris_src, vdst) ;
    if (!v)
    {
      ErrorPrintf(ERROR_BADPARM, 
                  "MRISsphericalCopy: cannot map dst vno %d", vno) ;
      continue ;
    }
    if (v-mris_src->vertices == Gdiag_no)
      DiagBreak() ;
    vdst->marked++ ;
    vdst->val2 = v->val2/(float)v->marked ;
    switch (which)
    {
    case VERTEX_AREA:
      vdst->origarea = v->origarea/(float)v->marked ;
      break ;
    case VERTEX_CURV:
      vdst->curv = v->curv ;
      break ;
    case VERTEX_VALS:
      vdst->val = v->val ;
      break ;
    }
  }
  if (mht_src) 
    MHTfree(&mht_src) ;

  return(NO_ERROR) ;
}
#endif

int
MRISorigAreaToCurv(MRI_SURFACE *mris)
{
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = v->origarea ;
  }
  return(NO_ERROR) ;
}

int
MRISareaToCurv(MRI_SURFACE *mris)
{
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = v->area ;
  }
  return(NO_ERROR) ;
}

int
MRISclearOrigArea(MRI_SURFACE *mris)
{
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->origarea = 0 ;
  }
  return(NO_ERROR) ;
}


int
MRISclear(MRI_SURFACE *mris, int which)
{
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    switch (which)
    {
    case VERTEX_AREA:
      v->origarea = 0 ;
      break ;
    case VERTEX_CURV:
      v->curv = 0 ;
      break ;
    case VERTEX_VAL:
      v->val = 0 ;
      break ;
    }
    v->val2 = 0 ;
  }
  return(NO_ERROR) ;
}
int
MRISnormalize(MRI_SURFACE *mris, int dof, int which)
{
  int     vno ;
  VERTEX  *v ;
  float   fdof = (float)dof, mean ;

  if (dof <= 0)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "MRISnormalize: invalid dof %d", dof)) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    switch (which)
    {
    case VERTEX_AREA:
      v->origarea /= fdof ; mean = v->origarea ;
      break ;
    case VERTEX_CURV:
      v->curv /= fdof ;     mean = v->curv ;
      break ;
    default:
    case VERTEX_VAL:
      v->val /= fdof ;      mean = v->val ;
      break ;
    }
    v->val2 = v->val2/fdof - mean*mean ;
  }
  return(NO_ERROR) ;
}

static int
mrisComputeGrayWhiteBorderDistributions(MRI_SURFACE *mris, MRI *mri, DEFECT *defect, 
																				HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_border,
																				HISTOGRAM *h_grad)
{
	int    vno, n2, n, i, nvertices, bin ;
	VERTEX *v, *vn, *vn2 ;
	HISTOGRAM *h_white_raw, *h_gray_raw, *h_border_raw, *h_grad_raw ;
	double    grad, min_grad, max_grad, bin_val, bin_size ;

	HISTOclear(h_gray, h_gray) ; HISTOclear(h_white, h_white) ; 
	HISTOclear(h_border, h_border) ;
	HISTOclear(h_grad, h_grad) ;

	mrisMarkDefect(mris, defect, 1);  /* avoid vertices in the defect */
	for (bin = 0 ; bin < h_gray->nbins ; bin++)
		h_gray->bins[bin] = h_white->bins[bin] = h_border->bins[bin] = bin ;

	min_grad = 100000 ; max_grad = -min_grad ;

	for (nvertices = i = 0 ; i < defect->nchull  ; i++)
	{
		vno = defect->chull[i] ;
		v = &mris->vertices[vno] ;
		for (n = 0 ; n < v->vtotal ; n++)
		{
			vn = &mris->vertices[v->v[n]] ;
			for (n2 = 0 ; n2 < vn->vtotal ; n2++)
			{
				vn2 = &mris->vertices[vn->v[n2]] ;
				if (vn2->marked)  /* already processed */
					continue ;
				grad = vn2->val2 - vn2->val2bak ;
				if (grad < min_grad)
					min_grad = grad ;
				if (grad > max_grad)
					max_grad = grad ;
			}
		}
	}
	/* add one bin at either end for almost zero probability events */
	bin_size = (max_grad - min_grad) / (h_grad->nbins-2) ;
	h_grad->bin_size = bin_size ;
	for (bin_val = min_grad-bin_size, bin = 0 ; bin < h_grad->nbins ; bin++, bin_val += bin_size)
		h_grad->bins[bin] = bin_val ;

	min_grad = h_grad->bins[0] ;

	for (nvertices = i = 0 ; i < defect->nchull  ; i++)
	{
		vno = defect->chull[i] ;
		v = &mris->vertices[vno] ;
		for (n = 0 ; n < v->vtotal ; n++)
		{
			vn = &mris->vertices[v->v[n]] ;
			for (n2 = 0 ; n2 < vn->vtotal ; n2++)
			{
				vn2 = &mris->vertices[vn->v[n2]] ;
				if (vn2->marked)  /* already processed */
					continue ;
				nvertices++ ;

				if (vn2->val2bak < 70)
					DiagBreak() ;
				bin = nint(vn2->val2) ;     h_white->counts[bin]++ ;     /* wm value */
				bin = nint(vn2->val2bak) ;  h_gray->counts[bin]++ ;      /* gray value */
				bin = nint(vn2->val) ;      h_border->counts[bin]++ ;      /* border value */
				grad = vn2->val2 - vn2->val2bak ;
				bin = (int)((grad - min_grad) / bin_size) ; 
				h_grad->counts[bin]++ ;

				vn2->marked = 1 ;   /* don't process it twice */
			}
		}
	}

	/* unmark them all */
	for (i = 0 ; i < defect->nchull  ; i++)
	{
		vno = defect->chull[i] ;
		v = &mris->vertices[vno] ;
		for (n = 0 ; n < v->vtotal ; n++)
		{
			vn = &mris->vertices[v->v[n]] ;
			for (n2 = 0 ; n2 < vn->vtotal ; n2++)
			{
				vn2 = &mris->vertices[vn->v[n2]] ;
				vn2->marked = 0 ;
			}
		}
	}
	for (bin = 0 ; bin < h_gray->nbins ; bin++)
	{
		if (h_gray->counts[bin] == 0)
			h_gray->counts[bin] = 0.1 ;
		if (h_white->counts[bin] == 0)
			h_white->counts[bin] = 0.1 ;
		if (h_border->counts[bin] == 0)
			h_border->counts[bin] = 0.1 ;
		if (h_grad->counts[bin] == 0)
			h_grad->counts[bin] = 0.1 ;
		h_grad->counts[bin] /= (float)nvertices ;
		h_gray->counts[bin] /= (float)nvertices ;
		h_white->counts[bin] /= (float)nvertices ;
		h_border->counts[bin] /= (float)nvertices ;
	}
	h_grad_raw = HISTOcopy(h_grad, NULL) ;
	h_gray_raw = HISTOcopy(h_gray, NULL) ;
	h_white_raw = HISTOcopy(h_white, NULL) ;
	h_border_raw = HISTOcopy(h_border, NULL) ;
	HISTOsmooth(h_grad_raw, h_grad, 2) ;
	HISTOsmooth(h_gray_raw, h_gray, 2) ;
	HISTOsmooth(h_white_raw, h_white, 2) ;
	HISTOsmooth(h_border_raw, h_border, 2) ;
	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
	{
		HISTOplot(h_gray, "g.plt") ; HISTOplot(h_white, "w.plt") ; HISTOplot(h_border, "b.plt") ;
		HISTOplot(h_gray_raw, "gr.plt") ; HISTOplot(h_white_raw, "wr.plt") ; HISTOplot(h_border_raw, "br.plt") ;
		HISTOplot(h_grad, "d.plt") ; HISTOplot(h_grad_raw, "dr.plt") ;
	}
	mrisMarkDefect(mris, defect, 0);
	HISTOfree(&h_gray_raw) ; HISTOfree(&h_white_raw) ; HISTOfree(&h_border_raw) ; HISTOfree(&h_grad_raw) ;
	return(NO_ERROR) ;
}
static int
mrisComputeJointGrayWhiteBorderDistributions(MRI_SURFACE *mris, MRI *mri, 
					     MRI *mri_gray_white, MRI *mri_wm)
{
  int    vno, x, y ;
  VERTEX *v ;
  float  norm ;
  Real   nx, ny, nz, xv, yv, zv, xw, yw, zw, white_val, gray_val ;
  
  MRIScomputeMetricProperties(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked || v->ripflag)  /* part of a defect - ignore it */
      continue ;
    
    nx = v->nx ; ny = v->ny ; nz = v->nz ; 
    xw = v->x ; yw = v->y ; zw = v->z ; 
    
    // MRIworldToVoxel(mri_wm, xw+.5*nx, yw+.5*ny, zw+.5*nz, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri_wm, xw+.5*nx, yw+.5*ny, zw+.5*nz, &xv, &yv, &zv) ;
    MRIsampleVolumeType(mri_wm, xv, yv, zv, &gray_val, SAMPLE_NEAREST) ;
    // MRIworldToVoxel(mri_wm, xw-.5*nx, yw-.5*ny, zw-.5*nz, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri_wm, xw-.5*nx, yw-.5*ny, zw-.5*nz, &xv, &yv, &zv) ;
    MRIsampleVolumeType(mri_wm, xv, yv, zv, &white_val, SAMPLE_NEAREST) ;
    if (gray_val >= MIN_WM_VAL && white_val >= MIN_WM_VAL)  /* white on both sides */
      continue ;
    
    MRIFvox(mri_gray_white, nint(v->val2), nint(v->val2bak), 0) += 1.0f ;
    
    if ((nint(nint(v->val2)) == 110) && (nint(v->val2bak) == 110))
      DiagBreak() ;
  }
  
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
    {
      if (FZERO(MRIFvox(mri_gray_white, x, y, 0)))
      {
	MRIFvox(mri_gray_white, x, y, 0) = 0.1 ;
      }
    }
  for (norm = 0.0, x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
    {
      norm += MRIFvox(mri_gray_white, x, y, 0) ;
    }
  
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
    {
      MRIFvox(mri_gray_white, x, y, 0) = MRIFvox(mri_gray_white, x, y, 0) / norm ;
    }
  
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_gray_white, "gw.mgh") ;
  }
  return(NO_ERROR) ;
}

static int
mrisFindGrayWhiteBorderMean(MRI_SURFACE *mris, MRI *mri)
{
  Real    x, y, z, xv, yv, zv, gray_val, white_val, nx, ny, nz ;
  int     vno ;
  VERTEX  *v ;

  MRIScomputeNormals(mris) ;
	/*  MRISsmoothSurfaceNormals(mris, 10) ;*/
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->marked)
      continue ;
		if (vno == Gdiag_no)
			DiagBreak() ;

    nx = v->nx ; ny = v->ny ; nz = v->nz ; 
    x = v->x ; y = v->y ; z = v->z ; 

    // MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &gray_val) ;
    // MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &white_val) ;
		v->val2 = white_val ; v->val2bak = gray_val ;
    v->val = (white_val + gray_val) / 2 ;
  }
#if 0
  MRISaverageVals(mris, 10) ;
  MRISaverageVal2s(mris, 10) ;
  MRISaverageVal2baks(mris, 10) ;
#else
  MRISmedianFilterVals(mris, 2) ;
  MRISmedianFilterVal2s(mris, 2) ;
  MRISmedianFilterVal2baks(mris, 2) ;
#endif
  return(NO_ERROR) ;
}

int
MRISreadNewCurvatureFile(MRI_SURFACE *mris, char *sname)
{
  int    k,vnum,fnum, vals_per_vertex ;
  float  curv, curvmin, curvmax;
  FILE   *fp;
  char   *cp, path[STRLEN], fname[STRLEN] ;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (cp)
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
    fprintf(stdout, "reading curvature file...") ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadBinaryCurvature: could not open %s", 
                 fname)) ;

  fread3(&vnum,fp);
  if (vnum != NEW_VERSION_MAGIC_NUMBER)
  {
    fclose(fp) ;
    return(MRISreadCurvatureFile(mris, fname)) ;
  }
  
  vnum = freadInt(fp);
  fnum = freadInt(fp);
  if (vnum!= mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadNewCurvature: incompatible vertex "
                 "number in file %s", fname)) ;
  }
  vals_per_vertex = freadInt(fp) ;
  if (vals_per_vertex != 1)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadNewCurvature(%s): vals/vertex %d unsupported (must be 1) ",
                 fname, vals_per_vertex)) ;
  }

  curvmin = 10000.0f ; curvmax = -10000.0f ;  /* for compiler warnings */
  for (k=0;k<vnum;k++)
  {
    curv = freadFloat(fp) ;
    if (k==0) curvmin=curvmax=curv;
    if (curv>curvmax) curvmax=curv;
    if (curv<curvmin) curvmin=curv;
    mris->vertices[k].curv = curv;
  }
  mris->max_curv = curvmax ;
  mris->min_curv = curvmin ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "done. min=%2.3f max=%2.3f\n", curvmin, curvmax) ;
  fclose(fp);
  return(NO_ERROR) ;
}
float *
MRISreadNewCurvatureVector(MRI_SURFACE *mris, char *sname)
{
  int    k,vnum,fnum, vals_per_vertex ;
  float  *cvec ;
  FILE   *fp;
  char   *cp, path[STRLEN], fname[STRLEN] ;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (cp)
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
    fprintf(stdout, "reading curvature file...") ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    return(NULL) ;

  fread3(&vnum,fp);
  if (vnum != NEW_VERSION_MAGIC_NUMBER)
  {
    fclose(fp) ;
    return(MRISreadCurvatureVector(mris, fname)) ;
  }
  
  vnum = freadInt(fp);
  fnum = freadInt(fp);
  if (vnum!= mris->nvertices)
  {
    fclose(fp) ;
    return(NULL) ;
  }
  vals_per_vertex = freadInt(fp) ;
  if (vals_per_vertex != 1)
  {
    fclose(fp) ;
    return(NULL) ;
  }

  cvec = (float *)calloc(mris->nvertices, sizeof(float)) ;
  if (!cvec)
    ErrorExit(ERROR_NOMEMORY, "MRISreadNewCurvatureVector(%s): calloc failed",
              fname) ;
  for (k=0;k<vnum;k++)
  {
    cvec[k] = freadFloat(fp) ;
  }
  fclose(fp);
  return(cvec) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISpaintVolume(MRI_SURFACE *mris, LTA *lta, MRI *mri)
{
  VERTEX   *v ;
  int      vno, width, height, depth ;
  Real     x, y, z, val ;
  MATRIX   *m_L, *m_ras_to_voxel ;
  VECTOR   *v_surf, *v_vol ;

  if (lta->type != LINEAR_RAS_TO_RAS)
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, 
                                    "MRISsampleVolume: unsupported transform type %d",
                                    lta->type)) ;

  v_surf = VectorAlloc(4, MATRIX_REAL) ;
  v_vol = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v_surf, 4, 1) = 1.0 ;
  *MATRIX_RELT(v_vol, 4, 1) = 1.0 ;
  m_ras_to_voxel = MRIgetRasToVoxelXform(mri) ;

  m_L = MatrixMultiply(m_ras_to_voxel, lta->xforms[0].m_L, NULL) ;
  width  = mri->width ; height  = mri->height ; depth  = mri->depth ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    v = &mris->vertices[vno] ;
    V3_X(v_surf) = v->x+.5*v->curv*v->nx ; 
    V3_Y(v_surf) = v->y+.5*v->curv*v->ny ; 
    V3_Z(v_surf) = v->z+.5*v->curv*v->nz ; 

    MatrixMultiply(m_L, v_surf, v_vol) ;
    x = V3_X(v_vol) ; y = V3_Y(v_vol) ; z = V3_Z(v_vol) ;

    MRIsampleVolume(mri, x, y, z, &val) ;
    v->val = val ;
    if (Gdiag_no == vno)
      printf("vertex %d at (%2.1f, %2.1f, %2.1f) --> voxel (%2.1f, %2.1f, %2.1f) = %2.2f\n",
             vno, v->x+.5*v->curv*v->nx, 
             v->y+.5*v->curv*v->ny, v->z+.5*v->curv*v->nz, x, y, z, val) ;
  }

  MatrixFree(&v_surf) ; MatrixFree(&v_vol) ; MatrixFree(&m_L) ; 
  MatrixFree(&m_ras_to_voxel) ;
  return(NO_ERROR) ;
}
static int
mrisComputePositioningGradients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  MHT  *mht_v_orig = NULL, *mht_v_current = NULL ;
  MRI  *mri_brain = parms->mri_brain ;
  int  avgs ;

  avgs = parms->n_averages ;
  if (!FZERO(parms->l_surf_repulse))
    mht_v_orig = MHTfillVertexTable(mris, NULL, ORIGINAL_VERTICES) ;
  if (!FZERO(parms->l_repulse))
    mht_v_current = MHTfillVertexTable(mris, mht_v_current,CURRENT_VERTICES);

  mrisComputeIntensityTerm(mris, parms->l_intensity, mri_brain, mri_brain,
                           parms->sigma);
  mrisComputeIntensityGradientTerm(mris, parms->l_grad,mri_brain,mri_brain);
  mrisComputeSurfaceRepulsionTerm(mris, parms->l_surf_repulse, mht_v_orig);
  
  mrisAverageGradients(mris, avgs) ;

  /* smoothness terms */
  mrisComputeSpringTerm(mris, parms->l_spring) ;
  mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm) ;
  mrisComputeRepulsiveTerm(mris, parms->l_repulse, mht_v_current) ;
  mrisComputeThicknessSmoothnessTerm(mris, parms->l_tsmooth) ;
  mrisComputeNormalSpringTerm(mris, parms->l_nspring) ;
  mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv) ;
  /*    mrisComputeAverageNormalTerm(mris, avgs, parms->l_nspring) ;*/
  /*    mrisComputeCurvatureTerm(mris, parms->l_curv) ;*/
  mrisComputeTangentialSpringTerm(mris, parms->l_tspring) ;
  
  
  if (mht_v_orig)
    MHTfree(&mht_v_orig) ;
  if (mht_v_current)
    MHTfree(&mht_v_current) ;
  return(NO_ERROR) ;
}
int
MRISallocExtraGradients(MRI_SURFACE *mris)
{
  if (mris->dx2)
    return(NO_ERROR) ;  /* already allocated */

  mris->dx2 = (float *)calloc(mris->nvertices, sizeof(float)) ;
  mris->dy2 = (float *)calloc(mris->nvertices, sizeof(float)) ;
  mris->dz2 = (float *)calloc(mris->nvertices, sizeof(float)) ;
  if (!mris->dx2 || !mris->dy2 || !mris->dz2)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISallocExtraGradients: could allocate gradient vectors") ;
  return(NO_ERROR) ;
}

int
MRISrestoreExtraGradients(MRI_SURFACE *mris)
{
  int   vno ;
  VERTEX *v ;

  if (!mris->dx2)
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->dx = mris->dx2[vno] ;
    v->dy = mris->dy2[vno] ;
    v->dz = mris->dz2[vno] ;
  }
  
  return(NO_ERROR) ;
}

int
MRISclearDistances(MRI_SURFACE *mris)
{
  int   vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->d = 0.0 ;
  }
  
  return(NO_ERROR) ;
}
/*-----------------------------------------------------------------
  MRISloadSurfVals() - loads surfaces values directly into an MRI
  structure. The surface value file can be any format read by
  MRIread. In addition, it can be a curv or paint file. If
  Surf is non-null, then it is used as a template; otherwise,
  the caller must spec the subject and hemisphere, and then the
  ?h.white surface is used as the template (and then freed).

  If the source file is neither curv nor paint, then MRIreadType is
  used to read the file in as a "volume", and the "volume" is reshaped
  to be nvertices X 1 X 1 X nframes (which is the final size
  regardless).

  If the source is curv format, then the given file is read from the
  subject's surf directory. If Surf is non-null, then the sujbect
  name and hemisphere in the MRI_SURFACE struct are used; otherwise
  they must be passed.

  If the subjectsdir string is NULL, then it reads SUBJECTS_DIR
  from the environment.
  -----------------------------------------------------------------*/
MRI *MRISloadSurfVals(char *srcvalfile, char *typestring, MRI_SURFACE *Surf, 
         char *subject, char *hemi, char *subjectsdir)
{
  MRI *SrcVals, *mritmp;
  char fname[2000];
  int srctype,reshapefactor=0,f;
  float *framepower = NULL;
  SXADAT *sxa;
  int freesurface = 0, err=0;

  if(Surf == NULL){
    /*-------- set SUBJECTS DIR -------------*/
    if(subjectsdir == NULL){
      subjectsdir = getenv("SUBJECTS_DIR");
      if(subjectsdir==NULL){
  fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
  return(NULL);
      }
    }
    /*------- load the surface -------------*/
    sprintf(fname,"%s/%s/surf/%s.white",subjectsdir,subject,hemi);
    printf("INFO: loading surface %s\n",fname);
    Surf = MRISread(fname);
    if(Surf == NULL){
      fprintf(stderr,"ERROR: could not read %s\n",fname);
      return(NULL);
    }
    freesurface = 1;
  }
  else{
    subject = Surf->subject_name;
    if(Surf->hemisphere == LEFT_HEMISPHERE)  hemi = "lh";
    if(Surf->hemisphere == RIGHT_HEMISPHERE) hemi = "rh";
  }
  
  /* ------------------ load the source data ----------------------------*/
  printf("Loading surface source data %s as %s\n",srcvalfile,typestring);
  if(!strcmp(typestring,"curv")){ /* curvature file */
    sprintf(fname,"%s/%s/surf/%s.%s",subjectsdir,subject,hemi,srcvalfile);
    printf("Reading curvature file %s\n",fname);
    err = MRISreadCurvatureFile(Surf, fname);
    if(err){
      printf("ERROR: reading curvature\n");
      return(NULL);
    }
    SrcVals = MRIcopyMRIS(NULL, Surf, 0, "curv");
    if(SrcVals == NULL){
      printf("ERROR: converting surface curv to MRI\n");
      return(NULL);
    }

    //SrcVals = MRIallocSequence(Surf->nvertices, 1, 1,MRI_FLOAT,1);
    //for(vtx = 0; vtx < Surf->nvertices; vtx++){
    //  MRIFseq_vox(SrcVals,vtx,0,0,0) = Surf->vertices[vtx].curv;
    //}
  }
  else if(!strcmp(typestring,"paint") || !strcmp(typestring,"w")){
    MRISreadValues(Surf,srcvalfile);
    SrcVals = MRIcopyMRIS(NULL, Surf, 0, "val");
    //SrcVals = MRIallocSequence(Surf->nvertices, 1, 1,MRI_FLOAT,1);
    //for(vtx = 0; vtx < Surf->nvertices; vtx++)
    //  MRIFseq_vox(SrcVals,vtx,0,0,0) = Surf->vertices[vtx].val;
  }
  else { /* Use MRIreadType */
    srctype = string_to_type(typestring);
    if(srctype == MRI_VOLUME_TYPE_UNKNOWN) {
      printf("ERROR: typestring %s unrecognized\n",typestring);
      return(NULL);
    }
    SrcVals =  MRIreadType(srcvalfile,srctype);
    if(SrcVals == NULL){
      printf("ERROR: could not read %s as type %d\n",srcvalfile,srctype);
      return(NULL);
    }
    if(SrcVals->height != 1 || SrcVals->depth != 1){
      reshapefactor = SrcVals->height * SrcVals->depth;
      printf("Reshaping %d\n",reshapefactor);
      mritmp = mri_reshape(SrcVals, reshapefactor*SrcVals->width, 
         1, 1, SrcVals->nframes);
      MRIfree(&SrcVals);
      SrcVals = mritmp;
    }

    if(SrcVals->width != Surf->nvertices){
      fprintf(stdout,"ERROR: dimesion inconsitency in source data\n");
      fprintf(stdout,"       Number of surface vertices = %d\n",
        Surf->nvertices);
      fprintf(stdout,"       Number of value vertices = %d\n",SrcVals->width);
      return(NULL);
    }
    if(is_sxa_volume(srcvalfile)){
      printf("INFO: Source volume detected as selxavg format\n");
      sxa = ld_sxadat_from_stem(srcvalfile);
      if(sxa == NULL) return(NULL);
      framepower = sxa_framepower(sxa,&f);
      if(f != SrcVals->nframes){
  fprintf(stderr," number of frames is incorrect (%d,%d)\n",
    f,SrcVals->nframes);
  return(NULL);
      }
      printf("INFO: Adjusting Frame Power\n");  fflush(stdout);
      mri_framepower(SrcVals,framepower);
    }
  }
  if(SrcVals == NULL){
    fprintf(stderr,"ERROR loading source values from %s\n",srcvalfile);
    return(NULL);
  }
  printf("Done Loading %s\n",srcvalfile);
  
  if(freesurface) MRISfree(&Surf);

  return(SrcVals);
}
/*-----------------------------------------------------------------
  MRIScopyMRI() - copies the data from an MRI_VOLUME struct into a
  field in the MRI_SURFACE vertex struct. The MRI_VOLUME struct is
  assumed to have the dimension: nvertices X 1 X 1 X nframes.  Frame
  is the zero-based frame number to copy. Field is a string that
  indicates which field of the vertex structure the data should be
  copied to. For example, "val" indicates the val field.  Other
  supported fields are: val, stat, valbak, val2, val2bak, imag_val,
  curv, curvbak, fsmask, nc. Others can be easily added. If there
  is an error, a 1 is returned; otherwise 0.
  -----------------------------------------------------------------*/
int MRIScopyMRI(MRIS *Surf, MRI *Src, int Frame, char *Field)
{
  int vtx, useval=0, usecurv=0;
  float val;

  if(Surf->nvertices != Src->width){
    printf("ERROR: MRIScopyMRI: Surf/Src dimension mismatch.\n");
    return(1);
  }

  if(Frame >= Src->nframes){
    printf("ERROR: MRIScopyMRI: requested frame number is too large.\n");
    printf("ERROR:   requested = %d, max = %d\n",Frame,Src->nframes);
    return(1);
  }

  /* A separate variable is used for val and curv for speed purposes */
  if(!strcmp(Field,"val")) useval = 1;
  else                     useval = 0;
  if(!strcmp(Field,"curv")) usecurv = 1;
  else                      usecurv = 0;

  /*------------------------------------------------*/
  for(vtx = 0; vtx < Surf->nvertices; vtx++){
    val = MRIgetVoxVal(Src, vtx, 0, 0,Frame);

    if(useval)                         Surf->vertices[vtx].val = val;
    else if(usecurv)                   Surf->vertices[vtx].curv = val;
    else if(!strcmp(Field,"stat"))     Surf->vertices[vtx].stat = val;
    else if(!strcmp(Field,"valbak"))   Surf->vertices[vtx].valbak = val;
    else if(!strcmp(Field,"val2"))     Surf->vertices[vtx].val2 = val;
    else if(!strcmp(Field,"val2bak"))  Surf->vertices[vtx].val2bak = val;
    else if(!strcmp(Field,"imag_val")) Surf->vertices[vtx].imag_val = val;
    else if(!strcmp(Field,"curvbak"))  Surf->vertices[vtx].curvbak = val;
    else if(!strcmp(Field,"fsmask"))   Surf->vertices[vtx].fsmask = val;
    else if(!strcmp(Field,"nc"))       Surf->vertices[vtx].nc = val;
    else printf("ERROR: MRIScopyMRI(): Field %s not supported\n",Field);

  }

  return(0);
}
/*-----------------------------------------------------------------
  MRIcopyMRIS() - copies the data from the given field of an
  MRI_SURFACE struct into a given frame of an MRI_VOLUME struct. The
  MRI_VOLUME should have the dimension: nvertices X 1 X 1 X nframes.
  Frame is the zero-based frame number to copy to. Field is a string
  that indicates which field of the vertex structure the data should
  be copied from. For example, "val" indicates the val field.  Other
  supported fields are: val, stat, valbak, val2, val2bak, imag_val,
  curv, curvbak, fsmask, nc. If mri is NULL, it will be allocated
  with nframes=Frame+1 (ie, just enough frames) and type will be
  MRI_FLOAT. A pointer to mri is returned. If an error occurs, 
  NULL is returned.
  -----------------------------------------------------------------*/
MRI *MRIcopyMRIS(MRI *mri, MRIS *surf, int Frame, char *Field)
{
  int vtx, useval=0, usecurv=0;
  float val;

  if(mri == NULL){
    mri = MRIallocSequence(surf->nvertices, 1, 1, MRI_FLOAT, Frame+1);
    if(mri==NULL){
      printf("ERROR: MRIcopyMRIS: could not alloc\n");
      return(NULL);
    }
  }
  else{
    if(surf->nvertices != mri->width){
      printf("ERROR: MRIcopyMRIS: surf/mri dimension mismatch.\n");
      return(NULL);
    }
    if(Frame >= mri->nframes){
      printf("ERROR: MRIScopyMRI: requested frame number is too large.\n");
      printf("ERROR:   requested = %d, max = %d\n",Frame,mri->nframes);
      return(NULL);
    }
  }

  /* A separate variable is used for val and curv for speed purposes */
  if(!strcmp(Field,"val")) useval = 1;
  else                     useval = 0;
  if(!strcmp(Field,"curv")) usecurv = 1;
  else                      usecurv = 0;

  /*------------------------------------------------*/
  for(vtx = 0; vtx < surf->nvertices; vtx++){

    if(useval)                         val = surf->vertices[vtx].val;
    else if(usecurv)                   val = surf->vertices[vtx].curv;
    else if(!strcmp(Field,"stat"))     val = surf->vertices[vtx].stat;
    else if(!strcmp(Field,"valbak"))   val = surf->vertices[vtx].valbak;
    else if(!strcmp(Field,"val2"))     val = surf->vertices[vtx].val2;
    else if(!strcmp(Field,"val2bak"))  val = surf->vertices[vtx].val2bak;
    else if(!strcmp(Field,"imag_val")) val = surf->vertices[vtx].imag_val;
    else if(!strcmp(Field,"curvbak"))  val = surf->vertices[vtx].curvbak;
    else if(!strcmp(Field,"fsmask"))   val = surf->vertices[vtx].fsmask;
    else if(!strcmp(Field,"nc"))       val = surf->vertices[vtx].nc;
    else {
      printf("ERROR: MRIScopyMRI(): Field %s not supported\n",Field);
      return(NULL);
    }
    MRIsetVoxVal(mri, vtx, 0, 0, Frame, val);
  }

  return(mri);
}
/*-------------------------------------------------------------------
  MRISsmoothMRI() - smooths values on the surface when the surface
  values are stored in an MRI_VOLUME structure with the number of
  columns (ie width) equal to the number of nvertices on the
  surface. The number of rows (height) and slices (depth) in the MRI
  struct should be 1. Can handle multiple frames. Can be performed
  in-place. If Targ is NULL, it will automatically allocate a new MRI
  strucutre.
  -------------------------------------------------------------------*/
MRI *MRISsmoothMRI(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *Targ)
{
  int nnbrs, nthstep, frame, vtx, nbrvtx, nthnbr;
  float val;
  MRI *SrcTmp;

  if(Surf->nvertices != Src->width){
    printf("ERROR: MRISsmooth: Surf/Src dimension mismatch\n");
    return(NULL);
  }

  if(Targ == NULL){
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth, 
            MRI_FLOAT, Src->nframes);
    if(Targ==NULL){
      printf("ERROR: MRISsmooth: could not alloc\n");
      return(NULL);
    }
  }
  else{
    if(Src->width   != Targ->width  || 
       Src->height  != Targ->height || 
       Src->depth   != Targ->depth  ||
       Src->nframes != Targ->nframes){
      printf("ERROR: MRISsmooth: output dimension mismatch\n");
      return(NULL);
    }
    if(Targ->type != MRI_FLOAT){
      printf("ERROR: MRISsmooth: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  SrcTmp = MRIcopy(Src,NULL);
  for(nthstep = 0; nthstep < nSmoothSteps; nthstep ++){
    printf("Step = %d\n",nthstep); fflush(stdout);

    for(vtx = 0; vtx < Surf->nvertices; vtx++){
      nnbrs = Surf->vertices[vtx].vnum;

      for(frame = 0; frame < Targ->nframes; frame ++){
  val = MRIFseq_vox(SrcTmp,vtx,0,0,frame);

  for(nthnbr = 0; nthnbr < nnbrs; nthnbr++){
    nbrvtx = Surf->vertices[vtx].v[nthnbr];
    val += MRIFseq_vox(SrcTmp,nbrvtx,0,0,frame) ;
  }/* end loop over neighbor */

  MRIFseq_vox(Targ,vtx,0,0,frame) = (val/(nnbrs+1));
      }/* end loop over frame */

    } /* end loop over vertex */

    MRIcopy(Targ,SrcTmp);
  }/* end loop over smooth step */

  MRIfree(&SrcTmp);

  return(Targ);
}

int
MRISrectifyCurvature(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->curv = fabs(v->curv) ;
  }
  mrisComputeCurvatureValues(mris) ;
  return(NO_ERROR) ;
}


#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double  
mrisAsynchronousTimeStep(MRI_SURFACE *mris, float momentum, 
                                    float delta_t, MHT *mht, float max_mag)
{
  static int direction = 1 ;
  double  mag, biggest_norm, scale ;
  int     vno, i, j ;
  VERTEX  *v ;
	static  int *vlist = NULL, *cropped ;
	static int nvertices = 0, ncropped = 0 ;

	if (mris->nvertices > nvertices)
	{
		if (vlist != NULL)
			free(vlist) ;
		vlist = (int *)calloc(mris->nvertices, sizeof(int)) ;
		cropped = (int *)calloc(mris->nvertices, sizeof(int)) ;
		nvertices = mris->nvertices ;
		if (!vlist || !cropped)
			ErrorExit(ERROR_NOMEMORY, "mrisAsynchronousTimeStep: could not allocate %d len buffers",
								nvertices) ;
	}

  /* take a step in the gradient direction modulated by momentum */
  if (mris->status == MRIS_RIGID_BODY)
  {
    mris->da = delta_t * mris->alpha + momentum * mris->da ;
    mris->db = delta_t * mris->beta + momentum * mris->db ;
    mris->dg = delta_t * mris->gamma + momentum * mris->dg ;
    MRISrotate(mris, mris, mris->da, mris->db, mris->dg) ;
  }
  else 
  {
		/* build a permutation of the vertex numbers */
		MRISclearMarks(mris) ;
		for (i = 0 ; i < ncropped ; i++)
		{
			v = &mris->vertices[cropped[i]] ; 
			v->marked = 1 ;
			vlist[i] = cropped[i] ;
		}

		markAllDistantConvex(mris, 4*max_mag) ; /* if distant, v->marked += 2 */
		for (vno = 0, i = ncropped ; vno < mris->nvertices ; vno++)
		{
			v = &mris->vertices[vno] ;
			if (v->marked != 2)
				continue ;   /* was cropped - already added */
			vlist[ncropped] = vno ; ncropped++ ;
		}

		for (vno = 0, i = ncropped ; vno < mris->nvertices ; vno++)
		{
			v = &mris->vertices[vno] ;
			if (v->marked)
				continue ;   /* was cropped - already added */
			vlist[i] = vno ; i++ ;
		}

		for (i = 0 ; i < mris->nvertices ; i++)
		{
			v = &mris->vertices[vlist[i]] ;
			if (v->marked)
				continue ;  /* don't let this one be swapped */
			j = randomNumber(0, mris->nvertices-1) ;
			if (mris->vertices[vlist[j]].marked)
				continue ;
			vno = vlist[i] ; vlist[i] = vlist[j] ; vlist[j] = vno ;
		}

		ncropped = 0 ;
		biggest_norm = 0.0 ;
		for (i = 0 ; i < mris->nvertices ; i++)
		{
			vno = vlist[i] ;
			v = &mris->vertices[vno] ;
			if (v->ripflag)
				continue ;
			if (vno == Gdiag_no)
				DiagBreak() ;
			v->odx = delta_t * v->dx + momentum*v->odx ;
			v->ody = delta_t * v->dy + momentum*v->ody ;
			v->odz = delta_t * v->dz + momentum*v->odz ;
			mag = sqrt(v->odx*v->odx + v->ody*v->ody + v->odz*v->odz) ;
			if (mag > biggest_norm)
				biggest_norm = mag ;
		}
		scale = max_mag / biggest_norm ;
#define MIN_SCALE (1.0/15.0)
		if (scale < MIN_SCALE)
			scale = MIN_SCALE ;
		for (i = 0 ; i < mris->nvertices ; i++)
		{
			vno = vlist[i] ;
			if (vno == Gdiag_no)
				DiagBreak() ;
			v = &mris->vertices[vno] ;
			mag = sqrt(v->odx*v->odx + v->ody*v->ody + v->odz*v->odz) ;
#if 1
			if (mag > max_mag) /* don't let step get too big */
			{
				mag = max_mag / mag ;
				if (v->marked)
					mag *= 2 ;   /* move distant convex vertices twice as far to avoid pinches */
			}
			else
				mag = 1 ;
#else
			if (mag * scale > max_mag)
				mag = max_mag/mag ;
			else
				mag = scale ;
#endif
      v->odx *= mag ; v->ody *= mag ; v->odz *= mag ;

			/* erase the faces this vertex is part of */
#if 0
			for (fno = 0 ; fno < v->num ; fno++)
				mrisEraseFace(mris, mri_filled, v->f[fno]) ;
#else
			if (mht)
				MHTremoveAllFaces(mht, mris, v) ;
#endif
			
			if (mht)
			{
				if (mrisLimitGradientDistance(mris, mht, vno) > 0)
				{
					if (vno == Gdiag_no)
						DiagBreak() ;
					cropped[ncropped++] = vno ;
				}
			}
			
			mag = sqrt(v->odx*v->odx + v->ody*v->ody + v->odz*v->odz) ;
			v->x += v->odx ; v->y += v->ody ; v->z += v->odz ;

			/* update distance-to-move (just approximate) */
			v->d -= mag ;
			if (v->d < 0)
				v->d = 0 ;
			
			if ((fabs(v->x) > 128.0f) ||
					(fabs(v->y) > 128.0f) ||
					(fabs(v->z) > 128.0f))
				DiagBreak() ;
			
			if (vno == Gdiag_no)
			{
				float dist, dot, dx, dy, dz ;
				
				dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
				dist = sqrt(dx*dx+dy*dy+dz*dz) ;
				dot = dx*v->nx + dy*v->ny + dz*v->nz ;
				fprintf(stdout, "%d: moving v %d by (%2.2f, %2.2f, %2.2f) dot=%2.2f-->"
								"(%2.1f, %2.1f, %2.1f)%s\n", i, vno, v->odx, v->ody, v->odz, 
								v->odx*v->nx+v->ody*v->ny+v->odz*v->nz,
								v->x, v->y, v->z, v->marked ? " CROPPED" : "") ;
				fprintf(stdout, "n = (%2.1f,%2.1f,%2.1f), total dist=%2.3f, total dot = %2.3f\n", 
								v->nx, v->ny, v->nz, dist, dot) ;
			}
			
			/* should this be done here????? (BRF) what about undoing step??? */
			v->dx = v->odx ;  /* for mrisTrackTotalDistances */
			v->dy = v->ody ;
			v->dz = v->odz ;

#if 0
			/* write the new face positions into the filled volume */
			for (fno = 0 ; fno < v->num ; fno++)
				mrisFillFace(mris, mri_filled, v->f[fno]) ;
#else
			if (mht)
				MHTaddAllFaces(mht, mris, v) ;
#endif
			
		}
	}

  direction *= -1 ;
  return(delta_t) ;
}

static int
markAllDistantConvex(MRI_SURFACE *mris, float min_dist)
{
	int   vno, n, num ;
	VERTEX *v, *vn ;
#if 0
	double nx, ny, nz, x, y, z, sx, sy, sz, nc ;
#endif

	for (num = vno = 0; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		if (v->ripflag || v->d < min_dist)
			continue ;
		if (vno == Gdiag_no)
			DiagBreak() ;

#if 0
    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    x = v->x ;    y = v->y ;   z = v->z ;
    sx = sy = sz = 0.0 ;

		for (n = 0; n < v->vnum ; n++)
		{
			vn = &mris->vertices[v->v[n]] ;
			sx += vn->x - x;
			sy += vn->y - y;
			sz += vn->z - z;
		}
    nc = sx*nx+sy*ny+sz*nz;   /* projection onto normal */
		if (nc > 0)  /* convex */
		{
			v->marked = 2 ;
			num++ ;
		}
		else
			v->marked = 0 ;
#else
		v->marked += 2 ;
		for (n = 0 ; n < v->vnum ; n++)
		{
			vn = &mris->vertices[v->v[n]] ;
			if (vn->d > v->d)
			{
				v->marked -= 2 ;
				break ;
			}
		}

		if (v->marked >= 2)
			num++ ;
#endif

		if (vno == Gdiag_no)
			printf("vertex %d %sdistant and convex\n",
						 Gdiag_no, v->marked ? "" : "NOT ") ;
	}

	printf("%d distant convex vertices found...\n", num) ;
	return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisAverageSignedGradients(MRI_SURFACE *mris, int num_avgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  dx, dy, dz, num, sigma, dot ;
  VERTEX *v, *vn ;
  MRI_SP *mrisp, *mrisp_blur ;

  if (num_avgs <= 0)
    return(NO_ERROR) ;

  if (Gdiag_no >= 0)
  {
    v = &mris->vertices[Gdiag_no] ;
    fprintf(stdout, "before averaging dot = %2.2f ",
            v->dx*v->nx+v->dy*v->ny+v->dz*v->nz) ;
  }
  if (0 && mris->status == MRIS_PARAMETERIZED_SPHERE)  /* use convolution */
  {
    sigma = sqrt((float)num_avgs) / 4.0 ;
    mrisp = MRISgradientToParameterization(mris, NULL, 1.0) ;
    mrisp_blur = MRISPblur(mrisp, NULL, sigma, -1) ;
    MRISgradientFromParameterization(mrisp_blur, mris) ;
    MRISPfree(&mrisp) ; MRISPfree(&mrisp_blur) ;
  }
  else for (i = 0 ; i < num_avgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      dx = v->dx ; dy = v->dy ; dz = v->dz ;
      pnb = v->v ;
      /*      vnum = v->v2num ? v->v2num : v->vnum ;*/
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
				dot = vn->dx * v->dx + vn->dy * v->dy + vn->dz*v->dz ;
				if (dot < 0)
					continue ;  /* pointing in opposite directions */

        num++ ;
        dx += vn->dx ; dy += vn->dy ; dz += vn->dz ;
#if 0
        if (vno == Gdiag_no)
        {
          float dot ;
          dot = vn->dx*v->dx + vn->dy*v->dy + vn->dz*v->dz ;
          if (dot < 0)
            fprintf(stdout, "vn %d: dot = %2.3f, dx = (%2.3f, %2.3f, %2.3f)\n",
                    v->v[vnb], dot, vn->dx, vn->dy, vn->dz) ;
        }
#endif
      }
      num++ ;
      v->tdx = dx / num ;
      v->tdy = dy / num ;
      v->tdz = dz / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->dx = v->tdx ; v->dy = v->tdy ; v->dz = v->tdz ;
    }
  }
  if (Gdiag_no >= 0)
  {
    float dot ;
    v = &mris->vertices[Gdiag_no] ;
    dot = v->nx*v->dx + v->ny*v->dy + v->nz*v->dz ;
    fprintf(stdout, " after dot = %2.2f\n",dot) ;
    if (fabs(dot) > 50)
      DiagBreak() ;
  }
  return(NO_ERROR) ;
}
#if 0
static int
mrisAverageWeightedGradients(MRI_SURFACE *mris, int num_avgs)
{
  int    vno, vlist[MAX_NBRS], nbrs, n, n2 ;
  float  nx, ny, nz, dx, dy, dz, sigma, wts[MAX_NBRS], total_wt, wt ;
  VERTEX *v, *vn, *vn2 ;
  MRI_SP *mrisp, *mrisp_blur ;

  if (num_avgs <= 0)
    return(NO_ERROR) ;

	sigma = sqrt((double)num_avgs) * M_PI / 2.0 ;
  if (Gdiag_no >= 0)
  {
    v = &mris->vertices[Gdiag_no] ;
    fprintf(stdout, "before averaging dot = %2.2f ",
            v->dx*v->nx+v->dy*v->ny+v->dz*v->nz) ;
  }
  if (0 && mris->status == MRIS_PARAMETERIZED_SPHERE)  /* use convolution */
  {
    sigma = sqrt((float)num_avgs) / 4.0 ;
    mrisp = MRISgradientToParameterization(mris, NULL, 1.0) ;
    mrisp_blur = MRISPblur(mrisp, NULL, sigma, -1) ;
    MRISgradientFromParameterization(mrisp_blur, mris) ;
    MRISPfree(&mrisp) ; MRISPfree(&mrisp_blur) ;
  }
  else
  {
		MRISclearMarks(mris) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
			if (vno == Gdiag_no)
				DiagBreak() ;

      if (v->ripflag)
        continue ;

			nx = v->nx ; ny = v->ny ; nz = v->nz ;
      dx = v->dx ; dy = v->dy ; dz = v->dz ;

			/* find all 1-neighbors */
			nbrs = 1 ; vlist[0] = vno ; wts[0] = 1.0 ;
			for (n = 0  ; n < v->vnum ; n++)
			{
				vn = &mris->vertices[v->v[n]] ;
				if (vn->marked)
					continue ;
				vn->marked = 1 ;
				wt = vn->nx*nx + vn->ny*ny + vn->nz*nz ; if (wt < 0) wt = 0 ;
				vlist[nbrs] = v->v[n] ; wts[nbrs] = exp(-1.0/(2.0f*sigma*sigma))*wt ; nbrs++ ;
			}
			/* find all 2-neighbors */
			for (n = v->vnum ; n < v->vtotal ; n++)
			{
				vn = &mris->vertices[v->v[n]] ;
				if (vn->marked || v->v[n] == vno)
					continue ;
				vn->marked = 2;
				wt = vn->nx*nx + vn->ny*ny + vn->nz*nz ; if (wt < 0) wt = 0 ;
				vlist[nbrs] = v->v[n] ; wts[nbrs] = exp(-4.0/(2.0f*sigma*sigma))*wt ; nbrs++ ;
			}
			/* find all 3-neighbors */
			for (n = v->vnum ; n < v->vtotal ; n++)
			{
				vn = &mris->vertices[v->v[n]] ;
				if (vn->marked != 2)  /* a two-neighbor */
					continue ;
				for (n2 = 0 ; n2 < vn->vnum ; n2++)
				{
					vn2 = &mris->vertices[vn->v[n2]] ;
					if (vn2->marked || vn->v[n2] == vno)
						continue ;
					vn2->marked = 3 ;
					wt = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ; if (wt < 0) wt = 0 ;
					vlist[nbrs] = vn->v[n2] ; wts[nbrs] = exp(-9.0/(2.0f*sigma*sigma))*wt ; nbrs++ ;
				}
				for (n2 = vn->vnum ; n2 < vn->vtotal ; n2++)
				{
					vn2 = &mris->vertices[vn->v[n2]] ;
					if (vn2->marked || vn->v[n2] == vno)
						continue ;
					vn2->marked = 4 ;
					wt = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ; if (wt < 0) wt = 0 ;
					vlist[nbrs] = vn->v[n2] ; wts[nbrs] = exp(-16.0/(2.0f*sigma*sigma))*wt ; nbrs++ ;
				}
			}

			v->tdx = v->tdy = v->tdz = 0.0 ;
			for (total_wt = 0.0, n = 0 ; n < nbrs ; n++)
			{
				wt = wts[n] ; total_wt += wt ;
				vn = &mris->vertices[vlist[n]] ;
				if (vlist[n] == Gdiag_no || vlist[n] == Gx)
					DiagBreak() ;
				vn->marked = 0 ;
				v->tdx += wt * vn->dx ; v->tdy += wt * vn->dy ; v->tdz += wt * vn->dz ;
			}
			v->tdx /= total_wt ; v->tdy /= total_wt ; v->tdz /= total_wt ; 
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->dx = v->tdx ; v->dy = v->tdy ; v->dz = v->tdz ;
    }
  }
  if (Gdiag_no >= 0)
  {
    float dot ;
    v = &mris->vertices[Gdiag_no] ;
    dot = v->nx*v->dx + v->ny*v->dy + v->nz*v->dz ;
    fprintf(stdout, " after dot = %2.2f\n",dot) ;
    if (fabs(dot) > 50)
      DiagBreak() ;
  }
  return(NO_ERROR) ;
}
#endif

#if 0
static int
mrisMarkSulcalVertices(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z ;
  Real    val0, xw,yw,zw ;
  double  del0, dot ;

	MRISclearFlags(mris, VERTEX_SULCAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ; y = v->y ; z = v->z ;
    
    // MRIworldToVoxel(parms->mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(parms->mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolume(parms->mri_brain, xw, yw, zw, &val0) ;
		dot = v->dx * v->nx + v->dy * v->ny + v->dz * v->nz ;
		
    del0 = v->val - val0 ;
		if (dot < 0 && del0 > 0)  /* too bright and moving inward */
		{
			v->flags |= VERTEX_SULCAL ;
			if (vno == Gdiag_no)
				printf("v %d: intensity %2.1f darker than target %2.1f - marked as sulcal\n",
							 vno, val0, v->val) ;
		}
  }
  
	return(NO_ERROR) ;
}
static int
mrisUpdateSulcalGradients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     vno, num ;
  VERTEX  *v ;
  double  dot ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    if (v->flags & VERTEX_SULCAL)
    {
      dot = v->dx * v->nx + v->dy * v->ny + v->dz * v->nz ;
      if (dot > 0)  /* now moving outwards - take out normal component */
      {
	num++ ;
	v->dx -= dot*v->nx ;
	v->dy -= dot*v->ny ;
	v->dz -= dot*v->nz ;
	if (vno == Gdiag_no)
	  printf("v %d: removing normal component %2.3f to prevent sulcal crossing\n",
		 vno, dot) ;
      }
    }
  }

  printf("%d vertices detected in sulcal-crossing\n", num) ;
	return(NO_ERROR) ;
}
#endif

int
MRISsetFlags(MRI_SURFACE *mris, int flags)
{
  int    vno ;
  
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].flags |= flags ;
  return(NO_ERROR) ;
}

int
MRISclearFlags(MRI_SURFACE *mris, int flags)
{
  int    vno ;
  
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].flags &= (~flags) ;
  return(NO_ERROR) ;
}

static int
mrisReadAsciiCurvatureFile(MRI_SURFACE *mris, char *fname)
{
  FILE   *fp ;
  int    vno ;
  char   line[STRLEN], *cp ;
  VERTEX *v ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "%s could not open output file %s.\n",
				mrisReadAsciiCurvatureFile, fname)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    cp = fgetl(line, 100, fp) ;
    if (!cp)
      break ;
    if (sscanf(line, "%*d %*f %*f %*f %f\n", &v->curv) != 1)
      ErrorReturn(ERROR_BADFILE, 
		  (ERROR_BADFILE, 
		   "mrisReadAsciiCurvatureFile(%s): could not scan curvature from line '%s'",fname,line));
  }
  
  fclose(fp) ;
  return(NO_ERROR) ;
}

// expand the surface by "h" and create a volume which has "val" outside of this surface
unsigned long
MRISeraseOutsideOfSurface(float h,MRI* mri_dst,MRIS *mris,unsigned char val)
{
  int i,j,k,imnr; 
  float x0,y0,z0,x1,y1,z1,x2,y2,z2,d0,d1,d2,dmax,u,v;
  float px,py,pz,px0,py0,pz0,px1,py1,pz1;
  int numu,numv,totalfilled,newfilled;
  double tx,ty,tz;
  unsigned long brainsize;

  int width, height,depth;
  MRI *mri_buff;

  width=mri_dst->width;
  height=mri_dst->height;
  depth=mri_dst->depth;

  mri_buff= MRIalloc(width, height, depth, MRI_UCHAR) ;

  for (k=0;k<mris->nvertices;k++)
  {
    // cache the values
   mris->vertices[k].tx=mris->vertices[k].x;
   mris->vertices[k].ty=mris->vertices[k].y;
   mris->vertices[k].tz=mris->vertices[k].z;

   // expand by h using normal
   mris->vertices[k].x +=h*mris->vertices[k].nx;
   mris->vertices[k].y +=h*mris->vertices[k].ny;
   mris->vertices[k].z +=h*mris->vertices[k].nz;
  }


  for (k=0;k<mris->nfaces;k++)
  {
    // calculate three vertices
    x0 =mris->vertices[mris->faces[k].v[0]].x;    
    y0 =mris->vertices[mris->faces[k].v[0]].y;    
    z0 =mris->vertices[mris->faces[k].v[0]].z;    
    x1 =mris->vertices[mris->faces[k].v[1]].x;    
    y1 =mris->vertices[mris->faces[k].v[1]].y;    
    z1 =mris->vertices[mris->faces[k].v[1]].z;    
    x2 =mris->vertices[mris->faces[k].v[2]].x;    
    y2 =mris->vertices[mris->faces[k].v[2]].y;    
    z2 =mris->vertices[mris->faces[k].v[2]].z;
    // calculate the sides
    d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;
    numu = (int)(ceil(2*d0));
    numv = (int)(ceil(2*dmax));

      
    for (v=0;v<=numv;v++)
    {
      px0 = x0 + (x2-x0)*v/numv;
      py0 = y0 + (y2-y0)*v/numv;
      pz0 = z0 + (z2-z0)*v/numv;
      px1 = x1 + (x2-x1)*v/numv;
      py1 = y1 + (y2-y1)*v/numv;
      pz1 = z1 + (z2-z1)*v/numv;
      for (u=0;u<=numu;u++)
      {
        px = px0 + (px1-px0)*u/numu;
        py = py0 + (py1-py0)*u/numu;
        pz = pz0 + (pz1-pz0)*u/numu;

	// MRIworldToVoxel(mri_dst,px,py,pz,&tx,&ty,&tz);
	MRIsurfaceRASToVoxel(mri_dst,px,py,pz,&tx,&ty,&tz);
	
	imnr=(int)(tz+0.5);
	j=(int)(ty+0.5);
	i=(int)(tx+0.5);
	if (i>=0 && i<width && j>=0 && j<height && imnr>=0 && imnr<depth)
	  MRIvox(mri_buff,i,j,imnr) = 255;
                                
      }  
    }
  }

  MRIvox(mri_buff,1,1,1)= 64;
  totalfilled = newfilled = 1;
  while (newfilled>0)
  {
    newfilled = 0;
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++)
          if (MRIvox(mri_buff,i,j,k)==0)
            if (MRIvox(mri_buff,i,j,mri_buff->zi[k-1])==64||MRIvox(mri_buff,i,mri_buff->yi[j-1],k)==64||
                MRIvox(mri_buff,mri_buff->xi[i-1],j,k)==64)
            {
              MRIvox(mri_buff,i,j,k)= 64;
              newfilled++;
            }
    for (k=depth-1;k>=0;k--)
      for (j=height-1;j>=0;j--)
        for (i=width-1;i>=0;i--)
          if (MRIvox(mri_buff,i,j,k)==0)
            if (MRIvox(mri_buff,i,j,mri_buff->zi[k+1])==64||MRIvox(mri_buff,i,mri_buff->yi[j+1],k)==64||
                MRIvox(mri_buff,mri_buff->xi[i+1],j,k)==64)
	    {
	      MRIvox(mri_buff,i,j,k) = 64;
	      newfilled++;
	    }
    totalfilled += newfilled;
  }

  // modify mri_dst so that outside = 0
  brainsize=0;
  if(val==0)
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++)
        {
          if (MRIvox(mri_buff,i,j,k)==64)
            MRIvox(mri_dst,i,j,k) = 0;
	  else
	    brainsize++;
        }
  else{
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++)
        {
          if (MRIvox(mri_buff,i,j,k)!=64)
            MRIvox(mri_dst,i,j,k) = val;
	  else
	    brainsize++;
        }
  }
  // restore the surface
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].x=mris->vertices[k].tx;
    mris->vertices[k].y=mris->vertices[k].ty;
    mris->vertices[k].z=mris->vertices[k].tz;
  }
  // calculate the normals
  MRIScomputeNormals(mris);
        
  MRIfree(&mri_buff);
  return brainsize;
}

int
MRISspringTermWithGaussianCurvature(MRI_SURFACE *mris, double gaussian_norm, double l_spring)
{
  int     vno, n, m ;
  VERTEX  *vertex, *vn ;
  float   sx, sy, sz, x, y, z, scale ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = vertex->x ;    y = vertex->y ;   z = vertex->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < vertex->vnum ; m++)
    {
      vn = &mris->vertices[vertex->v[m]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
    scale = pow(vertex->K, gaussian_norm) ;
    if (scale > 1)
      scale = 1 ;
    scale *= l_spring ;
    sx *= scale ;              /* move in normal direction */
    sy *= scale ;
    sz *= scale ;
    
    vertex->dx += sx ;
    vertex->dy += sy ;
    vertex->dz += sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring normal term:  (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
  }

  return(NO_ERROR) ;
}

static int
mrisComputeShrinkwrapTerm(MRI_SURFACE *mris, MRI *mri_brain, double  l_shrinkwrap)
{
	int    vno ;
	Real   xw, yw, zw, x, y, z, val, dx, dy, dz ;
	VERTEX *v ;
	float  min_val, max_val, target_val, delta ;

	MRIvalRange(mri_brain, &min_val, &max_val) ;
	target_val = (min_val + max_val) / 2 ;
	for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		target_val = v->val ;
		x = v->x ; y = v->y ; z = v->z ;
    MRIsurfaceRASToVoxel(mri_brain, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;
		delta = (val - target_val) ;
		dx = delta * v->nx * l_shrinkwrap ;
		dy = delta * v->ny * l_shrinkwrap ;
		dz = delta * v->nz * l_shrinkwrap ;

		v->dx += dx ; v->dy += dy ; v->dz += dz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d shrinkwrap term: (%2.3f, %2.3f, %2.3f), target %2.1f, MRI %2.1f, del=%2.1f, N=(%2.1f, %2.1f, %2.1f)\n",
              vno, dx, dy, dz, target_val, val, delta, v->nx, v->ny, v->nz) ;
	}
	return(NO_ERROR) ;
}

static double
mrisComputeShrinkwrapError(MRI_SURFACE *mris, MRI *mri_brain, double l_shrinkwrap)
{
#if 0
	static int iter = 100 ;
	int    vno ;
	Real   xw, yw, zw, x, y, z, val ;
	VERTEX *v ;
	float  min_val, max_val, target_val, error ;
	double sse ;

	MRIvalRange(mri_brain, &min_val, &max_val) ;
	target_val = (min_val + max_val) / 2 ;
	sse = 0 ;
	for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		sse += iter ;
	}
	iter-- ;
	return(sse) ;
#else
	return(0.0) ;
#endif
}
/*-------------------------------------------------------------
  MRISavgInterVetexDist() - computes the average and stddev of
  the distance between neighboring vertices. If StdDev is NULL,
  it is ignored.
  -------------------------------------------------------------*/
double MRISavgInterVetexDist(MRIS *Surf, double *StdDev)
{
  double Avg, Sum, Sum2, d;
  VERTEX *vtx1,*vtx2;
  int nNNbrs, nthNNbr, NbrVtxNo, VtxNo;
  long N;

  Sum = 0;
  Sum2 = 0;
  N = 0;
  for(VtxNo = 0; VtxNo < Surf->nvertices; VtxNo++){
    vtx1 = &Surf->vertices[VtxNo] ;
    nNNbrs = Surf->vertices[VtxNo].vnum;
    for(nthNNbr = 0; nthNNbr < nNNbrs; nthNNbr++){
      NbrVtxNo = Surf->vertices[VtxNo].v[nthNNbr];
      vtx2 = &Surf->vertices[NbrVtxNo] ;
      d = sqrt( (vtx1->x-vtx2->x)*(vtx1->x-vtx2->x) +
		(vtx1->y-vtx2->y)*(vtx1->y-vtx2->y) +
		(vtx1->z-vtx2->z)*(vtx1->z-vtx2->z) );
      Sum  += d;
      Sum2 += (d*d);
      N++;
    }
  }
  Avg = Sum/N;
  if(StdDev != NULL) *StdDev = sqrt( N*(Sum2/N - Avg*Avg)/(N-1) );

  //printf("\n\nN = %ld, Sum = %g, Sum2 = %g, Avg=%g, Std = %g\n\n",
  // N,Sum,Sum2,Avg,*StdDev);

  return(Avg);
}

/*-------------------------------------------------------------
  MRISavgVetexRadius() - computes the average and stddev of
  the distance from the origin to each vertex. If StdDev is NULL,
  it is ignored.
  -------------------------------------------------------------*/
double MRISavgVetexRadius(MRIS *Surf, double *StdDev)
{
  double Avg, Sum, Sum2, d;
  VERTEX *vtx;
  int VtxNo, N;

  Sum = 0;
  Sum2 = 0;
  for(VtxNo = 0; VtxNo < Surf->nvertices; VtxNo++){
    vtx = &Surf->vertices[VtxNo] ;
    d = sqrt( vtx->x*vtx->x + vtx->y*vtx->y + vtx->z*vtx->z );
    Sum  += d;
    Sum2 += (d*d);
  }

  N = Surf->nvertices;
  Avg = Sum/N;
  if(StdDev != NULL) *StdDev = sqrt( N*(Sum2/N - Avg*Avg)/(N-1) );

  //printf("\n\nN = %ld, Sum = %g, Sum2 = %g, Avg=%g, Std = %g\n\n",
  // N,Sum,Sum2,Avg,*StdDev);

  return(Avg);
}

