#include <stdio.h>
#include <math.h>

#include "diag.h"
#include "error.h"
#include "utils.h"
#include "macros.h"
#include "fio.h"
#include "mri.h"
#include "mrisurf.h"
#include "matrix.h"
#include "proto.h"
#include "stats.h"
#include "timer.h"
#include "const.h"
#include "mrishash.h"

/*---------------------------- STRUCTURES -------------------------*/

/*---------------------------- CONSTANTS -------------------------*/

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
#define MAX_NEG_AREA_PCT    0.01f

#define MRIS_BINARY_FILE    0
#define MRIS_ASCII_FILE     1
#define MRIS_GEO_FILE       2    /* movie.byu format */

#define NEG_AREA_K          200.0
/* limit the size of the ratio so that the exp() doesn't explode */
#define MAX_NEG_RATIO       (400 / NEG_AREA_K)

/*------------------------ STATIC PROTOTYPES -------------------------*/

static int mrisComputeCurvatureValues(MRI_SURFACE *mris) ;
static int mrisNormalDirectionTriangleIntersection(MRI_SURFACE *mris,VERTEX *v,
                                                   MHT *mht, double *pdist);
static int  load_triangle_vertices(MRI_SURFACE *mris, int fno, double U0[3], 
                                   double U1[3], double U2[3]) ;
#if 0
static double mrisFindClosestFilledVoxel(MRI_SURFACE *mris, MRI *mri_filled, 
                                         int vno, double max_dist) ;
#endif
static int   mrisCheck(MRI_SURFACE *mris) ;
static int   mrisClipGradient(MRI_SURFACE *mris, float max_len) ;
static int   mrisClipMomentumGradient(MRI_SURFACE *mris, float max_len) ;
static int   mrisFileNameType(char *fname) ;
static int   mrisComputeNormals(MRI_SURFACE *mris) ;
static int   mrisComputeSurfaceDimensions(MRI_SURFACE *mris) ;
static int   mrisFindNeighbors(MRI_SURFACE *mris) ;
static void  mrisNormalize(float v[3]) ;
static float mrisTriangleArea(MRIS *mris, int fac, int n) ;
static int   mrisNormalFace(MRIS *mris, int fac,int n,float norm[]) ;
static int   mrisReadTransform(MRIS *mris, char *mris_fname) ;
static MRI_SURFACE *mrisReadAsciiFile(char *fname) ;
static MRI_SURFACE *mrisReadGeoFile(char *fname) ;
static int         mrisReadGeoFilePositions(MRI_SURFACE *mris,char *fname) ;
static MRI_SURFACE *mrisReadTriangleFile(char *fname) ;
static int         mrisReadTriangleFilePositions(MRI_SURFACE*mris,
                                                  char *fname) ;

/*static int   mrisReadFieldsign(MRI_SURFACE *mris, char *fname) ;*/
static double mrisComputeSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double mrisComputeNonlinearAreaSSE(MRI_SURFACE *mris) ;
static double mrisComputeSpringEnergy(MRI_SURFACE *mris) ;
static double mrisComputeIntensityError(MRI_SURFACE *mris, 
                                        INTEGRATION_PARMS *parms);
static double mrisComputeIntensityGradientError(MRI_SURFACE *mris, 
                                        INTEGRATION_PARMS *parms);
static double mrisComputeSphereError(MRI_SURFACE *mris, 
                                     INTEGRATION_PARMS *parms);
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
static int   mrisIntegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, 
                           int n_avgs);
static int   mrisIntegrationEpoch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                                  int n_avgs);
static int   mrisRemoveNegativeArea(MRI_SURFACE *mris,INTEGRATION_PARMS *parms,
                                    int n_avgs, float min_area_pct,
                                    int max_passes);
static double mrisLineMinimize(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisLineMinimizeSearch(MRI_SURFACE *mris, 
                                     INTEGRATION_PARMS *parms);
static double  mrisMomentumTimeStep(MRI_SURFACE *mris, float momentum, 
                                    float dt, float tol, float n_averages) ;
static double  mrisAsynchronousTimeStep(MRI_SURFACE *mris, float momentum, 
                                    float dt, MHT *mht) ;
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
static int   mrisComputeIntensityTerm(MRI_SURFACE*mris,double l_intensity,
                                      MRI *mri_brain, MRI *mri_smooth);
static int   mrisComputeIntensityGradientTerm(MRI_SURFACE*mris,
                                              double l_grad,
                                              MRI *mri_brain, MRI *mri_smooth);
static int   mrisComputeSphereTerm(MRI_SURFACE *mris, double l_sphere, 
                                   float radius) ;
static int   mrisComputeExpansionTerm(MRI_SURFACE *mris, double l_expand) ;
static int   mrisComputeDistanceTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeCorrelationTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, 
                                              double l_curv) ;
static int   mrisComputePolarCorrelationTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeAngleAreaTerms(MRI_SURFACE *mris, 
                                       INTEGRATION_PARMS *parms) ;
static int   mrisComputeNonlinearAreaTerm(MRI_SURFACE *mris, 
                                       INTEGRATION_PARMS *parms) ;
static int   mrisClearDistances(MRI_SURFACE *mris) ;
static int   mrisClearGradient(MRI_SURFACE *mris) ;
static int   mrisClearMomentum(MRI_SURFACE *mris) ;
static int   mrisApplyGradient(MRI_SURFACE *mris, double dt) ;
static int   mrisValidVertices(MRI_SURFACE *mris) ;
static int   mrisValidFaces(MRI_SURFACE *mris) ;
static int   mrisLabelVertices(MRI_SURFACE *mris, float cx, float cy, 
                               float cz, int label, float radius) ;


static double mrisFindNormalDistance(MRI_SURFACE *mris, MHT *mht, int vno, 
                                      double max_dist);

static int mrisProjectSurface(MRI_SURFACE *mris) ;
static int mrisOrientSurface(MRI_SURFACE *mris) ;
static int   mrisComputeBoundaryNormals(MRI_SURFACE *mris) ;
static int   mrisSmoothBoundaryNormals(MRI_SURFACE *mris, int niter) ;
static int   mrisFlipPatch(MRI_SURFACE *mris) ;

/* not currently used */
static int  mrisNeighborAtVoxel(MRI_SURFACE *mris, MRI *mri, int vno, 
                                int xv,int yv,int zv) ;


#if 0
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
                                             double l_spring);
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
static int   mrisWriteSnapshot(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                               int t) ;
static int   mrisTrackTotalDistance(MRI_SURFACE *mris) ;
static int  mrisLimitGradientDistance(MRI_SURFACE *mris, MHT *mht, int vno) ;
static int mrisFillFace(MRI_SURFACE *mris, MRI *mri, int fno) ;
static int mrisHatchFace(MRI_SURFACE *mris, MRI *mri, int fno, int on) ;
#if 0
static int mrisEraseFace(MRI_SURFACE *mris, MRI *mri, int fno) ;
#endif
static double mrisRmsValError(MRI_SURFACE *mris, MRI *mri) ;
static int mrisRemoveTriangleLinks(MRI_SURFACE *mris) ;
static int mrisRemoveVertexLink(MRI_SURFACE *mris, int vno1, int vno2) ;

/*--------------------------------------------------------------------*/

/*--------------------- CONSTANTS AND MACROS -------------------------*/

#if 0
#define QUAD_FILE_MAGIC_NUMBER      16777215
#else
#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#endif
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_VERSION_MAGIC_NUMBER  16777215
#define START_Y                   (-128)
#define SLICE_THICKNESS           1

#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y,\
                                               v1->z-v0->z)

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


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISread(char *fname)
{
  MRI_SURFACE *mris = NULL ;
  int         nquads, nvertices, magic, version, ix, iy, iz, vno, fno, n, m,
              imnr, imnr0, imnr1, type, vertices[VERTICES_PER_FACE+1] ;
  float       x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  FILE        *fp = NULL ;
  VERTEX      *vertex ;
  FACE        *face ;

  type = mrisFileNameType(fname) ;
  if (type == MRIS_ASCII_QUADRANGLE_FILE)
  {
    mris = mrisReadAsciiFile(fname) ;
    if (!mris)
      return(NULL) ;
    version = -1 ;
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE)
  {
    mris = mrisReadGeoFile(fname) ;
    if (!mris)
      return(NULL) ;
    version = -1 ;
  }
  else
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
        fprintf(stderr, "new surface file format\n");
    }
    else if (magic == TRIANGLE_FILE_MAGIC_NUMBER)
    {
      fclose(fp) ;
      mris = mrisReadTriangleFile(fname) ;
      if (!mris)
        ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n")) ;
      version = -2 ;
    }
    else
    {
      rewind(fp);
      version = 0;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        printf("surfer: old surface file format\n");
    }
  }
  if (version >= -1)  /* some type of quadrangle file */
  {
    fread3(&nvertices, fp);
    fread3(&nquads, fp);   /* # of qaudrangles - not triangles */
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr,"reading %d vertices and %d faces.\n",
      nvertices,2*nquads);
    
    mris = MRISalloc(nvertices, 2*nquads) ;
    mris->type = MRIS_BINARY_QUADRANGLE_FILE ;

    imnr0 = 1000 ;
    imnr1 = 0 ;
    for (vno = 0 ; vno < nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      fread2(&ix,fp);
      fread2(&iy,fp);
      fread2(&iz,fp);
      vertex->x = ix/100.0;
      vertex->y = iy/100.0;
      vertex->z = iz/100.0;
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
        fread1(&vertex->num,fp);   /* # of faces we are part of */
        vertex->f = (int *)calloc(vertex->num,sizeof(int));
        if (!vertex->f)
          ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                    vertex->num) ;
        vertex->n = (int *)calloc(vertex->num,sizeof(int));
        if (!vertex->n)
          ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d nbrs",
                    vertex->n) ;
        for (n=0;n<vertex->num;n++)
          fread3(&vertex->f[n],fp);
      } 
      else 
        vertex->num = 0;   /* will figure it out */
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
  if (strstr(fname, "rh"))
    mris->hemisphere = RIGHT_HEMISPHERE ;
  else
    mris->hemisphere = LEFT_HEMISPHERE ;
  if ((version<0) || type == MRIS_ASCII_QUADRANGLE_FILE)
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
        (int *)calloc(mris->vertices[vno].num,sizeof(int));
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
  mrisComputeNormals(mris);
  mrisComputeVertexDistances(mris) ;
  mrisReadTransform(mris, fname) ;
  if (type == MRIS_ASCII_QUADRANGLE_FILE || type == MRIS_GEO_TRIANGLE_FILE)
  {
    MRISsetNeighborhoodSize(mris, 2) ;
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
  }
  else 
  {
    if (MRISreadBinaryCurvature(mris, fname) != NO_ERROR)
    {
      fprintf(stderr, "computing surface curvature directly...\n") ;
      MRISsetNeighborhoodSize(mris, 2) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
    }
       
    if (MRISreadBinaryAreas(mris, fname) != NO_ERROR)
     fprintf(stderr, "ignoring area file...\n") ; /*return(NULL) ;*/
  }

  mris->radius = MRISaverageRadius(mris) ;
  if (IS_QUADRANGULAR(mris))
    mrisRemoveTriangleLinks(mris) ;
  MRIScomputeMetricProperties(mris) ;
  /*  mrisFindPoles(mris) ;*/

  MRISstoreCurrentPositions(mris) ;
  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISfastRead(char *fname)
{
  MRI_SURFACE *mris ;
  int         nquads, nvertices, magic, version, ix, iy, iz, vno, fno, n, m,
              imnr, imnr0, imnr1, type, vertices[4] ;
  float       x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  FILE        *fp ;
  VERTEX      *vertex ;
  FACE        *face ;

  mris = NULL ; fp = NULL ;
  type = mrisFileNameType(fname) ;
  if (type == MRIS_ASCII_QUADRANGLE_FILE)
  {
    mris = mrisReadAsciiFile(fname) ;
    if (!mris)
      return(NULL) ;
    version = -3 ;
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
      mris = mrisReadTriangleFile(fname) ;
      if (!mris)
        ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n")) ;
      version = -2 ;
    }
    else if (magic == QUAD_FILE_MAGIC_NUMBER) 
    {
      version = -1;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stderr, "new surface file format\n");
    }
    else 
    {
      rewind(fp);
      version = 0;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        printf("surfer: old surface file format\n");
    }
  }
  if (version >= -1)
  {
    fread3(&nvertices, fp);
    fread3(&nquads, fp);
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr,
              "reading %d vertices and %d faces.\n",nvertices,2*nquads);
    
    mris = MRISalloc(nvertices, 2*nquads) ;
    mris->type = MRIS_BINARY_QUADRANGLE_FILE ;
    
    imnr0 = 1000 ;
    imnr1 = 0 ;
    for (vno = 0 ; vno < nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      fread2(&ix,fp);
      fread2(&iy,fp);
      fread2(&iz,fp);
      vertex->x = ix/100.0;
      vertex->y = iy/100.0;
      vertex->z = iz/100.0;
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
        fread1(&vertex->num,fp);
        vertex->f = (int *)calloc(vertex->num,sizeof(int));
        if (!vertex->f)
          ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                    vertex->num) ;
        vertex->n = (int *)calloc(vertex->num,sizeof(int));
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
  if (strstr(fname, "rh"))
    mris->hemisphere = RIGHT_HEMISPHERE ;
  else
    mris->hemisphere = LEFT_HEMISPHERE ;
  if ((version<0) || type == MRIS_ASCII_QUADRANGLE_FILE)
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
        (int *)calloc(mris->vertices[vno].num,sizeof(int));
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
  mrisComputeNormals(mris);
#if 0
  mrisComputeVertexDistances(mris) ;
  mrisReadTransform(mris, fname) ;
#endif
  if (type == MRIS_ASCII_QUADRANGLE_FILE || type == MRIS_GEO_TRIANGLE_FILE)
  {
    MRISsetNeighborhoodSize(mris, 2) ;
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
  }
  else 
  {
    if (MRISreadBinaryCurvature(mris, fname) != NO_ERROR)
      fprintf(stderr, "ignoring curvature file...\n") ; /*return(NULL) ;*/
#if 0
    if (MRISreadBinaryAreas(mris, fname) != NO_ERROR)
      return(NULL) ;
#endif
  }

  if (IS_QUADRANGULAR(mris))
    mrisRemoveTriangleLinks(mris) ;
  MRISstoreCurrentPositions(mris) ;
  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwrite(MRI_SURFACE *mris, char *fname)
{
  int   k, type ;
  float x,y,z;
  FILE  *fp;

  type = mrisFileNameType(fname) ;
  if (type == MRIS_ASCII_QUADRANGLE_FILE)
    return(MRISwriteAscii(mris, fname)) ;
  else if (type == MRIS_GEO_TRIANGLE_FILE)
    return(MRISwriteGeo(mris, fname)) ;
  if (mris->type == MRIS_TRIANGULAR_SURFACE)
    return(MRISwriteTriangularSurface(mris, fname)) ;
  fp = fopen(fname,"w");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE,"MRISwrite(%s): can't create file\n",fname));
  fwrite3(QUAD_FILE_MAGIC_NUMBER,fp);
  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces/2,fp);   /* # of quadrangles */
  for (k = 0 ; k < mris->nvertices ; k++)
  {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
    fwrite2((int)(x*100),fp);
    fwrite2((int)(y*100),fp);
    fwrite2((int)(z*100),fp);
  }
  for (k = 0 ; k < mris->nfaces ; k+=2)
  {
    fwrite3(mris->faces[k].v[0],fp);
    fwrite3(mris->faces[k].v[1],fp);
    fwrite3(mris->faces[k+1].v[0],fp);
    fwrite3(mris->faces[k].v[2],fp);
  }
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
              "MRISalloc(%d, %d): could not allocate mris structure");

  mris->nsize = 1 ;  /* only 1-connected neighbors initially */
  mris->nvertices = nvertices ;
  mris->nfaces = nfaces ;
  mris->vertices = (VERTEX *)calloc(nvertices, sizeof(VERTEX)) ;
  if (!mris->vertices)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISalloc(%d, %d): could not allocate vertices");
  mris->faces = (FACE *)calloc(nfaces, sizeof(FACE)) ;
  if (!mris->faces)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRISalloc(%d, %d): could not allocate faces");
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
    fprintf(stderr, "finding surface neighbors...") ;

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
        vtmp[v->vnum++] = f->v[n0];
      for (i=0;i<v->vnum && vtmp[i]!=f->v[n1];i++);
      if (i==v->vnum)
        vtmp[v->vnum++] = f->v[n1];
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
        printf("face[%d].v[%d] = %d, but face %d not in vertex %d face list\n",
               k,m,f->v[m], k, f->v[m]);
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
    fprintf(stderr, "sampling long-range distances") ;
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
    fprintf(stderr, 
            "\nsampling %d dists/vertex (%2.1f at each dist) = %2.1fMB\n",
            vtotal, 
            (float)vtotal/((float)max_nbhd-(float)mris->nsize),
            (float)vtotal*mrisValidVertices(mris)*sizeof(float)*3.0f / 
            (1024.0f*1024.0f)) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((Gdiag & DIAG_HEARTBEAT) && (!(vno % (mris->nvertices/25))))
      fprintf(stderr, " %%%1.0f", 100.0f*(float)vno / (float)mris->nvertices) ;
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
        for (min_dist = 1000000.0f, n2 = 0 ; n2 < vn->vnum ; n2++)
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
          vnbrs[i] = -1 ;
        }
      }
    }


    if ((Gdiag_no == vno) && DIAG_VERBOSE_ON)
    {
      FILE  *fp ;
      char  fname[200] ;

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
    for (n = 0 ; n < vtotal ; n++)
    {
      if (FZERO(v->dist_orig[n]))
        fprintf(stderr, "zero distance at v %d, n %d (vn = %d)\n",
                vno, n, v->v[n]) ;
    }
  }

  mris->avg_nbrs = (float)total_nbrs / (float)mrisValidVertices(mris) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "avg_nbrs = %2.1f\n", mris->avg_nbrs) ;

#if MULTI_DIST_SCALING
  if (Gdiag & DIAG_SHOW)
  {
    for (n = 0 ; n <= max_nbhd ; n++)
    {
      if (nc[n])
        c[n] /= (float)nc[n] ;
      fprintf(stderr, "c[%d] = %2.5f (%d samples)\n", n, c[n], nc[n]) ;
    }
    fprintf(stderr, "c[] = { ") ;
    for (n = 0 ; n <= max_nbhd ; n++)
    {
      fprintf(stderr, "%2.5f", c[n]) ;
      if (n < max_nbhd)
        fprintf(stderr, ", ") ;
    }
  }
#endif
  free(vnbrs) ;
  free(vall) ;
  VectorFree(&v1) ; VectorFree(&v2) ;
  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stderr, " done.\n") ;
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
    fprintf(stderr, 
            "\nsampling %d dists/vertex (%2.1f at each dist) = %2.1fMB\n",
            vtotal, 
            (float)vtotal/((float)max_nbhd-(float)mris->nsize),
            (float)vtotal*mris->nvertices*sizeof(float)*3.0f / 
            (1024.0f*1024.0f)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((Gdiag & DIAG_HEARTBEAT) && (!(vno % (mris->nvertices/25))))
      fprintf(stderr, " %%%1.0f", 100.0f*(float)vno / (float)mris->nvertices) ;
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
      char  fname[200] ;

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
    fprintf(stderr, "avg_nbrs = %2.1f\n", mris->avg_nbrs) ;

  free(vnbrs) ;
  free(vall) ;
  VectorFree(&v1) ; VectorFree(&v2) ;
  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stderr, " done.\n") ;
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
        fprintf(stderr, "v %d: vnum=%d, v2num=%d, vtotal=%d\n",
                vno, v->vnum, v->v2num, v->vtotal) ;
        for (n = 0 ; n < neighbors ; n++)
          fprintf(stderr, "v[%d] = %d\n", n, v->v[n]) ;
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
    fprintf(stderr, "avg_nbrs = %2.1f\n", mris->avg_nbrs) ;
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
    fprintf(stderr, "removing ripped vertices and faces...") ;
  do
  {
    nripped = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
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
            memmove(v->n+fno, v->n+fno+1, (v->num-fno-1)*sizeof(int)) ;
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

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;

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

static int
mrisComputeNormals(MRI_SURFACE *mris) 
{
  int       k,n, num ;
  VERTEX    *v ;
  FACE      *f;
  float     norm[3],snorm[3], len ;

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
      mrisNormalFace(mris, v->f[n],v->n[n],norm);
      snorm[0] += norm[0];
      snorm[1] += norm[1];
      snorm[2] += norm[2];

      /* Note: overestimates area by *2 !! */
      v->area += mrisTriangleArea(mris, v->f[n],v->n[n]); 
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
    }
    else
    {
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
    if (v->ripflag)
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
  char transform_fname[200], fpref[300] ;

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
  char   fname[200], fpref[200], hemi[20] ;

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
  char   *cp, path[200], fname[200] ;
  
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
    fprintf(stderr, "reading curvature file...") ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISreadBinaryCurvature: could not open %s", 
                 fname)) ;

  fread3(&vnum,fp);
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
    fprintf(stderr, "done. min=%2.3f max=%2.3f\n", curvmin, curvmax) ;
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
  char  fname[200], fpref[200], hemi[20] ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "reading area file...") ;

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
    fprintf(stderr, "total area = %2.0f.\n", mris->orig_area) ;
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
  int   k;
  FILE  *fp;
  char  *cp, fname[200], path[200] ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "writing area file...") ;

  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/%s", path, sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */

  /*  mris->orig_area = 0.0f ;*/
  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADPARM, 
              (ERROR_BADPARM,"mrisWriteBinaryAreas: no area file %s\n",fname));
  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);

  for (k=0;k<mris->nvertices;k++)
    fwriteFloat(mris->vertices[k].origarea, fp);
  fclose(fp);

  /* hack to correct for overestimation of area in compute_normals */
#if 0
  mris->orig_area /= 2; 
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "total area = %2.0f.\n", mris->orig_area) ;
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
    fprintf(stderr, 
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
    fprintf(stderr, 
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
      vdst->n = (int *)calloc(vdst->num,sizeof(int));
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
    char fname[200] ;

    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    parms->fp = fopen(fname, "w") ;
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
        fprintf(stderr, "setting momentum=%2.1f, dt=%2.1f, l_dist=%2.2f\n",
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
    mrisClearGradient(mris) ;      
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
      delta_t = mrisMomentumTimeStep(mris, parms->momentum, parms->dt, 
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
      fprintf(stderr, "%3.3d: count: %d (%2.2f%%), area: %2.2f (%2.2f%%)   \n",
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
    fprintf(stderr, "\n") ;
    if (Gdiag & DIAG_WRITE)
      fclose(parms->fp) ;
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
    char fname[200] ;

    sprintf(fname, "%s.out", parms->base_name) ;
    parms->fp = fopen(fname, "w") ;
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
    starting_sse = mrisComputeSSE(mris, parms) ;
    for (total_steps = 0, n_averages = base_averages; !done ;n_averages /= 2)
    {
      steps = mrisIntegrate(mris, parms, n_averages) ;
      parms->start_t += steps ;
      total_steps += steps ;
      done = n_averages == 0 ;   /* finished integrating at smallest scale */
    }
    parms->dt = parms->base_dt ;         /* reset time step */
    ending_sse = mrisComputeSSE(mris, parms) ;
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
static float sigmas[] = { 2.0f, 1.0f, 0.5f, 0.25 } ;
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


#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])
#define PARAM_IMAGES         (IMAGES_PER_SURFACE*SURFACES)

int
MRISregister(MRI_SURFACE *mris, MRI_SP *mrisp_template, 
             INTEGRATION_PARMS *parms, int max_passes)
{
  float   sigma ;
  int     i, /*steps,*/ done, sno, ino, msec ;
  MRI_SP  *mrisp ;
  char    fname[200], base_name[200], path[200] ;
  double  base_dt ;
  struct  timeb start ;
  static  int first = 1 ;
  
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
      parms->fp = fopen(fname, "w") ;
    mrisLogIntegrationParms(parms->fp, mris,parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris,parms) ;

  MRISuseMeanCurvature(mris) ;
  MRISnormalizeCurvature(mris) ;
  MRISstoreMeanCurvature(mris) ;

  for (sno = 1 ; sno < SURFACES ; sno++)
  {
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
      fprintf(stderr, "calculating curvature of %s surface\n",fname) ;
    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp,"calculating curvature of %s surface\n",fname);

    if (first)  /* only do rigid alignment first time through */
    {
      first = 0 ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "finding optimal rigid alignment\n") ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp, "finding optimal rigid alignment\n") ;
      parms->mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
      parms->mrisp_template = mrisp_template ;
      MRISrigidBodyAlignGlobal(mris, parms, 4.0f, 16.0f, 8) ;
      MRISPfree(&parms->mrisp) ;
    }
    else if (parms->flags & IP_USE_CURVATURE)
    {
      /* only small adjustments needed after 1st time around */
      parms->tol *= 2.0f ;
      parms->l_corr /= 2.0f ;
      mrisLogIntegrationParms(parms->fp, mris, parms) ;
      mrisLogIntegrationParms(stderr, mris, parms) ;
    }
    else
      break ;   /* finished */

    for (i = 0 ; i < NSIGMAS ; i++)
    {
      parms->sigma = sigma = sigmas[i] ;
      parms->dt = base_dt ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "blurring surfaces with sigma=%2.2f...", sigma) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,"correlating surfaces with with sigma=%2.2f\n",
                sigma) ;
      MRISuseMeanCurvature(mris) ;
      mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
      parms->mrisp = MRISPblur(mrisp, NULL, sigma, 0) ;
      parms->mrisp_template = MRISPblur(mrisp_template, NULL, sigma, ino) ;
      MRISPblur(parms->mrisp_template, NULL, sigma, ino+1) ; /* variances */
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "done.\n") ;
      /* normalize curvature intensities for both source and target */
      MRISfromParameterization(parms->mrisp_template, mris, ino);
      MRISnormalizeCurvature(mris) ;
      MRIStoParameterization(mris, parms->mrisp_template, 1, ino) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        sprintf(fname, "%s/%s.%4.4dtarget%2.2f", 
                path, mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",
                parms->start_t, sigma) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, "writing curvature file %s...", fname) ;
        MRISwriteCurvature(mris, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, "done.\n") ;
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
          fprintf(stderr, "writing curvature file %s...", fname) ;
        MRISwriteCurvature(mris, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, "done.\n") ;
        sprintf(fname, "target.%s.%4.4d.hipl",parms->base_name,parms->start_t);
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, "writing parameterization file %s...", fname) ;
        MRISPwrite(parms->mrisp_template, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, "done.\n") ;
      }
      mrisClearMomentum(mris) ;
      done = 0 ;
      mrisIntegrationEpoch(mris, parms, parms->n_averages) ;
    }
  }

  parms->tol = 1e-2 ;  /* remove everything possible pretty much */
  mrisRemoveNegativeArea(mris,parms,parms->n_averages,MAX_NEG_AREA_PCT,3);
  MRISPfree(&parms->mrisp) ; MRISPfree(&parms->mrisp_template) ;
  msec = TimerStop(&start) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "registration took %2.2f hours\n",
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
    char fname[200] ;
    
    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
            parms->base_name);
    if (!parms->start_t)
    {
      parms->fp = fopen(fname, "w") ;
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
        fprintf(stderr, "%d: %d | ", i, nbrs[i]) ;
    fprintf(stderr, "\n") ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

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
        fprintf(stderr, "resampling long-range distances") ;
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
        fprintf(stderr, "outputting distance errors to distance.log...\n") ;
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
        fprintf(stderr, "disturbing distances by %%%2.1f\n", (float)max_pct) ;
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

    if (!passno)
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
      fprintf(stderr, 
              "pass %d: epoch %d of %d starting distance error %%%2.2f\n",
              passno+1, i+1, (int)(NCOEFS), (float)pct_error);
    
      parms->l_dist = dist_coefs[i] ;
      parms->l_area = area_coefs[i] ;
      parms->l_angle = ANGLE_AREA_SCALE * parms->l_area ;
      if (i == NCOEFS-1)  /* see if distance alone can make things better */
        starting_sse = mrisComputeSSE(mris, parms) ;
      mrisIntegrationEpoch(mris, parms, base_averages) ;
    }

    parms->l_area = area_coefs[NCOEFS-1] ;
    parms->l_dist = dist_coefs[NCOEFS-1] ;
    ending_sse = mrisComputeSSE(mris, parms) ;
    if (Gdiag & DIAG_SHOW)
    {
#if 0
      fprintf(stderr, "pass %d: start=%2.1f, end=%2.1f, ratio=%2.3f\n",
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


  /* finally, remove all the small holes */
  parms->l_area = 1.0f ;
  parms->l_dist = 0.1f ;  /* was 0.001 */
  parms->l_angle = ANGLE_AREA_SCALE * parms->l_area ;
  parms->niterations = niter ;
  parms->tol = 1e-2 ;   /* try and remove as much negative stuff as possible */
  fprintf(stderr, "removing remaining folds...\n") ;
  mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 3);

  if (mris->status == MRIS_PLANE)  /* smooth out remaining folds */
  {
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "smoothing final surface...\n") ;
    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp, "smoothing final surface...\n") ;
    parms->l_spring = 1.0f ;
    parms->l_area = 0.0f ;
    parms->niterations = 5 ;
    parms->integration_type = INTEGRATE_MOMENTUM ;
    parms->dt = 0.5f ; parms->momentum = 0.0f ;
    parms->n_averages = 0 ;
    mrisIntegrate(mris, parms, 0) ;
    /*    mrisRemoveNegativeArea(mris, parms, 0, MAX_NEG_AREA_PCT, 1);*/
  }
 

  pct_error = MRISpercentDistanceError(mris) ;
  if (Gdiag & DIAG_SHOW)
    mrisLogStatus(mris, parms, stderr, 0) ;
  fprintf(stderr, "final distance error %%%2.2f\n", (float)pct_error);
  mrisProjectSurface(mris) ;
  msec = TimerStop(&start) ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "optimization complete.\n") ;
    fprintf(stderr, "unfolding took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp, "unfolding took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
  }
  if (Gdiag & DIAG_WRITE)
  {
    mrisLogStatus(mris, parms, parms->fp, 0) ;
    fprintf(parms->fp, "final distance error %%%2.2f\n", pct_error);
    fclose(parms->fp) ;
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
  int     base_averages, niter, passno, msec ;
  double  pct_error ;
  struct  timeb start ;
  
  TimerStart(&start) ;
  
  if (Gdiag & DIAG_WRITE)
  {
    char fname[200] ;
    
    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
            parms->base_name);
    if (!parms->start_t)
    {
      parms->fp = fopen(fname, "w") ;
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "MRISunfold: could not open log file %s\n",
                  fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris,parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  /*  parms->start_t = 0 ;*/
/*
   integrate until no improvement can be made at ANY scale, or until
   the error is effectively zero.
*/
  base_averages = parms->n_averages ;
  niter = parms->niterations ;
  passno = 0 ;
  do
  {
    {
      char   *cp ;

      cp = getenv("SPHERE") ;
      if (cp)
        MRISstoreAnalyticDistances(mris, MRIS_SPHERE) ;
      cp = getenv("PLANE") ;
      if (cp)
        MRISstoreAnalyticDistances(mris, MRIS_PLANE) ;
    }

    mrisIntegrationEpoch(mris, parms, base_averages) ;
  } while (++passno < max_passes) ;



  pct_error = MRISpercentDistanceError(mris) ;
  if (Gdiag & DIAG_SHOW)
    mrisLogStatus(mris, parms, stderr, 0) ;
  fprintf(stderr, "final distance error %%%2.2f\n", (float)pct_error);
  mrisProjectSurface(mris) ;
  msec = TimerStop(&start) ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "optimization complete.\n") ;
    fprintf(stderr, "unfolding took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp, "unfolding took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
  }
  if (Gdiag & DIAG_WRITE)
  {
    mrisLogStatus(mris, parms, parms->fp, 0) ;
    fprintf(parms->fp, "final distance error %%%2.2f\n", pct_error);
    fclose(parms->fp) ;
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
    char fname[200] ;
    
    sprintf(fname, "%s.out", parms->base_name) ;
    if (!parms->start_t)
      parms->fp = fopen(fname, "w") ;
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
        fprintf(stderr, "%d: %d | ", i, nbrs[i]) ;
    fprintf(stderr, "\n") ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  /*  parms->start_t = 0 ;*/
/*
   integrate until no improvement can be made at ANY scale, or until
   the error is effectively zero.
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
          fprintf(stderr, "resampling long-range distances") ;
        MRISsaveVertexPositions(mris, TMP_VERTICES) ;
        MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
        MRISsampleDistances(mris, nbrs, parms->nbhd_size) ;
        MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
        mrisClearMomentum(mris) ;
      }
      if (!i)  /* record starting error */
      {
        parms->l_area = area_scoefs[NSCOEFS-1] ;
        parms->l_dist = dist_scoefs[NSCOEFS-1] ;
        starting_sse = mrisComputeSSE(mris, parms) ;
      }


      /* remove any folds in the surface */
      mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 2) ;

      parms->l_dist = dist_scoefs[i] ;
      parms->l_area = area_scoefs[i] ;
      parms->l_angle = ANGLE_AREA_SCALE * area_scoefs[i] ;
      mrisIntegrationEpoch(mris, parms, base_averages) ;
    }

    parms->l_area = area_scoefs[NSCOEFS-1] ;
    parms->l_dist = dist_scoefs[NSCOEFS-1] ;
    ending_sse = mrisComputeSSE(mris, parms) ;
    if (Gdiag & DIAG_SHOW)
    {
#if 0
      fprintf(stderr, "pass %d: start=%2.1f, end=%2.1f, ratio=%2.3f\n",
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
      fprintf(stderr,"epoch took %2.2f minutes\n",(float)msec/(1000.0f*60.0f));
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,
                "epoch took %2.2f minutes\n",(float)msec/(1000.0f*60.0f));
    }
  } while (!FZERO(ending_sse) && 
           (((starting_sse-ending_sse)/starting_sse) > parms->tol)) ;

#if 0
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, 
            "initial unfolding complete - settling to equilibrium...\n") ;
  parms->niterations = 1000 ;  /* let it go all the way to equilibrium */
  mrisIntegrationEpoch(mris, parms, parms->n_averages = 0) ; 
#endif

  /* finally, remove all the small holes */
  parms->l_area = 1.0f ;
  parms->l_dist = 0.001f ;
  parms->l_angle = ANGLE_AREA_SCALE * area_scoefs[0] ;
  parms->niterations = niter ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "removing remaining folds...\n") ;
  mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 3);
  if (Gdiag & DIAG_SHOW)
    mrisLogStatus(mris, parms, stderr, 0) ;
  if (Gdiag & DIAG_WRITE)
    mrisLogStatus(mris, parms, parms->fp, 0) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "unfolding complete.\n") ;
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
static float neg_area_ratios[] = { 0.01f, 0.001f, 0.0001f } ;
#define MAX_PASSES (sizeof(neg_area_ratios) / sizeof(neg_area_ratios[0]))
static int
mrisRemoveNegativeArea(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, 
                       int base_averages, float min_area_pct, int max_passes)
{
  int    total_steps, done, steps, n_averages, old_averages, npasses, niter ;
  float  pct_neg, ratio ;
  char   *snum, *sdenom ;
  float  l_area, l_parea, l_corr, l_spring, l_dist, *pnum, *pdenom, cmod ;
  double tol ;
  
  if (Gdiag & DIAG_WRITE && parms->fp == NULL)
  {
    char fname[200] ;

    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    parms->fp = fopen(fname, "w") ;
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
  l_spring = parms->l_spring ;
  l_dist = parms->l_dist ; l_corr = parms->l_corr ;
  parms->l_area = parms->l_parea = parms->l_dist = 
    parms->l_corr = parms->l_spring = 0.0 ;

  /* there is one negative area removing term (area, parea, spring), and one
     term we are seaking to retain (corr, dist).
     */
  cmod = 1.0f ;
  if (!FZERO(l_corr))
  { sdenom = "corr" ; pdenom = &parms->l_corr  ; cmod = 10.0f ; }
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
      fprintf(stderr, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    if (Gdiag & DIAG_WRITE && (n_averages == base_averages))
      fprintf(parms->fp, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    parms->n_averages = n_averages ;
    steps = mrisIntegrate(mris, parms, n_averages) ;
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
  mrisComputeNormals(mris) ;
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
  else
  { snum = "spring" ; pnum = &parms->l_spring  ; }

  if (!FZERO(*pdenom))
  {
    ratio = *pnum / *pdenom ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
  }

  old_averages = parms->n_averages ;
  for (done = total_steps = 0, n_averages = base_averages ; !done ; 
       n_averages /= 4)
  {
    parms->n_averages = n_averages ;
    steps = mrisIntegrate(mris, parms, n_averages) ;
    parms->start_t += steps ;
    total_steps += steps ;
    done = n_averages == 0 ;    /* finished integrating at smallest scale */
    if (mris->status == MRIS_SPHERE)
    {
      parms->scale *= parms->dt_decrease ;
      if (parms->scale < 1.0f)
        parms->scale = 1.0f ;
    }
  }
#if 0
  mrisComputeNormals(mris) ;
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
static int
mrisIntegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_averages)
{
  int     t, write_iterations, niterations, nsmall, neg ;
  double  l_dist, l_area, l_spring, sse, old_sse, delta_t, total_small = 0.0, 
          sse_thresh, pct_neg, pct_neg_area, total_vertices, tol
          /*, scale, last_neg_area */ ;

  l_spring = parms->l_spring ;
  l_dist = parms->l_dist ;
  l_area = parms->l_area ;
  write_iterations = parms->write_iterations ;
  niterations = parms->niterations ;
  tol = parms->tol * sqrt(((double)n_averages + 1.0) / 1024.0);
  sse_thresh = tol ;
    
  mrisProjectSurface(mris) ;
  MRIScomputeMetricProperties(mris) ;

#if AVERAGE_AREAS
  MRISreadTriangleProperties(mris, mris->fname) ;
  mrisAverageAreas(mris, n_averages, ORIG_AREAS) ;
#endif

  parms->starting_sse = sse = old_sse = mrisComputeSSE(mris, parms) ;

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
    fprintf(stderr,
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
    if (!FZERO(parms->l_curv))
      MRIScomputeSecondFundamentalForm(mris) ;

    mrisClearGradient(mris) ;      /* clear old deltas */
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

    mrisAverageGradients(mris, n_averages) ;
    mrisComputeSpringTerm(mris, parms->l_spring) ;
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
      delta_t = mrisMomentumTimeStep(mris, parms->momentum, parms->dt, 
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
      fprintf(stderr, "rotating brain by (%2.1f, %2.1f, %2.1f)\n",
              alpha, beta, gamma) ;
    }
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    if (Gdiag_no >= 0)
      fprintf(stderr, "v %d curvature = %2.3f\n",
              Gdiag_no, mris->vertices[Gdiag_no].H) ;
    /* only print stuff out if we actually took a step */
    sse = mrisComputeSSE(mris, parms) ;
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

    if ((write_iterations > 0) &&!((t+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
      mrisWriteSnapshot(mris, parms, t+1) ;
    if (mris->status == MRIS_PLANE && mris->neg_area > 4*mris->total_area)
    {
      fprintf(stderr, "flipping flattened patch...\n") ;
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

  parms->ending_sse = mrisComputeSSE(mris, parms) ;
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
          sse_spring, sse_corr ;
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
  sse_spring = mrisComputeSpringEnergy(mris) ;
  sse_curv = MRIStotalVariation(mris) ;

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
  rms = mrisComputeSSE(mris, parms) ;
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
mrisComputeSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double  sse, sse_area, sse_angle, delta, sse_curv, sse_spring, sse_dist,
          area_scale, sse_corr, sse_neg_area, l_corr, sse_val, sse_sphere,
          sse_grad, sse_nl_area ;
  int     ano, fno ;
  FACE    *face ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  sse_nl_area = sse_corr = sse_angle = sse_neg_area = sse_val = sse_sphere =
    sse_area = sse_spring = sse_curv = sse_dist = sse_grad = 0.0 ;

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
  if (!FZERO(parms->l_narea))
    sse_nl_area = mrisComputeNonlinearAreaSSE(mris) ;
  if (!FZERO(parms->l_dist))
    sse_dist = mrisComputeDistanceError(mris) ;
  if (!FZERO(parms->l_spring))
    sse_spring = mrisComputeSpringEnergy(mris) ;
  if (!FZERO(parms->l_curv))
    sse_curv = MRIScomputeFolding(mris) ;
  l_corr = (double)(parms->l_corr + parms->l_pcorr) ;
  if (!FZERO(l_corr))
    sse_corr = mrisComputeCorrelationError(mris, parms, 1) ;
  if (!FZERO(parms->l_intensity))
    sse_val = mrisComputeIntensityError(mris, parms) ;
  if (!FZERO(parms->l_grad))
    sse_grad = mrisComputeIntensityGradientError(mris, parms) ;
  if (!FZERO(parms->l_sphere))
    sse_sphere = mrisComputeSphereError(mris, parms) ;

  sse = 
    (double)parms->l_area   * sse_neg_area + 
    (double)parms->l_sphere * sse_sphere + 
    (double)parms->l_intensity    * sse_val + 
    (double)parms->l_grad    * sse_grad + 
    (double)parms->l_parea  * sse_area + 
    (double)parms->l_narea  * sse_nl_area + 
    (double)parms->l_angle  * sse_angle + 
    (double)parms->l_dist   * sse_dist + 
    (double)parms->l_spring * sse_spring + 
    (double)l_corr          * sse_corr + 
    (double)parms->l_curv   * CURV_SCALE * sse_curv ;
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
    if (!FZERO(face->orig_area))
      ratio = area_scale*face->area / face->orig_area ;
    else
      ratio = 0.0f ;
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
    fprintf(stderr, "(%2.0f, %2.0f, %2.0f) --> (%2.0f, %2.0f, %2.0f), ctr "
            "(%2.0f, %2.0f, %2.0f)\n",
            mris->xlo, mris->ylo, mris->zlo, mris->xhi, mris->yhi, mris->zhi,
            mris->xctr, mris->yctr, mris->zctr);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "finding cortical poles...") ;

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
      fprintf(stderr, "F: (%2.0f,%2.0f,%2.0f), T: (%2.0f,%2.0f,%2.0f) "
              "O: (%2.0f,%2.0f,%2.0f).\n",
              mris->v_frontal_pole->x, mris->v_frontal_pole->y, 
              mris->v_frontal_pole->z,
              mris->v_temporal_pole->x, mris->v_temporal_pole->y, 
              mris->v_temporal_pole->z,
              mris->v_occipital_pole->x, mris->v_occipital_pole->y, 
              mris->v_occipital_pole->z) ;
    else
      fprintf(stderr, "F: (%2.0f,%2.0f,%2.0f), T: (NOT FOUND), "
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
    mris->total_area = mris->neg_area = 0.0f ;
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      if (face->ripflag)
        continue ;
      if (face->area >= 0.0f)
        mris->total_area += face->area ;
      else
        mris->neg_area += -face->area ;
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
  mris->total_area = mris->neg_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    if (face->area >= 0.0f)
      mris->total_area += face->area ;
    else
      mris->neg_area += -face->area ;
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
  char    fname[200], fpref[200], hemi[20], *cp ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "reading triangle files...") ;

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
    fprintf(stderr, 
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
    fprintf(stderr, "total area = %2.0f.\n", mris->orig_area) ;


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
  char    fname[200], fpref[200], hemi[20], *cp ;

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
        fprintf(stderr, "angle [%d][%d] = %2.1f\n",
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
#define MAX_ASYNCH_MM    0.2
#define MIN_MM           0.001

static double
mrisAdaptiveTimeStep(MRI_SURFACE *mris,INTEGRATION_PARMS *parms)
{
  double  delta_t, sse, starting_sse ;

  starting_sse = mrisComputeSSE(mris, parms) ;

  MRISstoreCurrentPositions(mris) ;
  delta_t = mrisMomentumTimeStep(mris, parms->momentum, parms->dt, parms->tol,
                                 parms->n_averages) ;

  sse = mrisComputeSSE(mris, parms) ;

  if (sse > starting_sse)  /* error increased - turn off momentum */
  {
    mrisClearMomentum(mris) ;
    parms->dt *= parms->dt_decrease ;
    if (parms->dt <= parms->base_dt)
      parms->dt = parms->base_dt ;
      
    if (sse / starting_sse > parms->error_ratio)  /* undo time step */
    {
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "sse increased by %2.0f%%, undoing time step...\n",
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
                                    float delta_t, MHT *mht)
{
  static int direction = 1 ;
  double  mag ;
  int     vno, i ;
  VERTEX  *vertex ;

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
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vertex->odx = delta_t * vertex->dx + momentum*vertex->odx ;
    vertex->ody = delta_t * vertex->dy + momentum*vertex->ody ;
    vertex->odz = delta_t * vertex->dz + momentum*vertex->odz ;
    mag = 
      sqrt(vertex->odx*vertex->odx + 
           vertex->ody*vertex->ody +
           vertex->odz*vertex->odz) ;
    if (mag > MAX_ASYNCH_MM) /* don't let step get too big */
    {
      mag = MAX_ASYNCH_MM / mag ;
      vertex->odx *= mag ; vertex->ody *= mag ; vertex->odz *= mag ;
    }
    if (vno == Gdiag_no)
      fprintf(stderr, "moving v %d by (%2.3f, %2.3f, %2.3f) --> "
              "(%2.1f, %2.1f, %2.1f), nz=%2.2f, a=%2.3f\n",
              vno, vertex->odx, vertex->ody, vertex->odz,
              vertex->x, vertex->y, vertex->z, vertex->nz,vertex->area) ;

    /* erase the faces this vertex is part of */
#if 0
    for (fno = 0 ; fno < vertex->num ; fno++)
      mrisEraseFace(mris, mri_filled, vertex->f[fno]) ;
#else
    if (mht)
      MHTremoveAllFaces(mht, mris, vertex) ;
#endif

    if (mht)
      mrisLimitGradientDistance(mris, mht, vno) ;

    vertex->x += vertex->odx ; 
    vertex->y += vertex->ody ;
    vertex->z += vertex->odz ;

    if ((fabs(vertex->x) > 128.0f) ||
        (fabs(vertex->y) > 128.0f) ||
        (fabs(vertex->z) > 128.0f))
      DiagBreak() ;
    vertex->dx = vertex->odx ;  /* for mrisTrackTotalDistances */
    vertex->dy = vertex->ody ;
    vertex->dz = vertex->odz ;

#if 0
    /* write the new face positions into the filled volume */
    for (fno = 0 ; fno < vertex->num ; fno++)
      mrisFillFace(mris, mri_filled, vertex->f[fno]) ;
#else
    if (mht)
      MHTaddAllFaces(mht, mris, vertex) ;
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
mrisMomentumTimeStep(MRI_SURFACE *mris, float momentum, float dt, float tol, 
                     float n_averages)
{
  double  delta_t, mag ;
  int     vno ;
  VERTEX  *vertex ;
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
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    dx = vertex->dx ; dy = vertex->dy ; dz = vertex->dz ;
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
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vertex->odx = delta_t * vertex->dx + momentum*vertex->odx ;
    vertex->ody = delta_t * vertex->dy + momentum*vertex->ody ;
    vertex->odz = delta_t * vertex->dz + momentum*vertex->odz ;
    mag = 
      sqrt(vertex->odx*vertex->odx + 
           vertex->ody*vertex->ody +
           vertex->odz*vertex->odz) ;
    if (mag > MAX_MOMENTUM_MM) /* don't let step get too big */
    {
      mag = MAX_MOMENTUM_MM / mag ;
      vertex->odx *= mag ; vertex->ody *= mag ; vertex->odz *= mag ;
    }
    if (vno == Gdiag_no)
      fprintf(stderr, "moving v %d by (%2.3f, %2.3f, %2.3f) --> "
              "(%2.1f, %2.1f, %2.1f), nz=%2.2f, a=%2.3f\n",
              vno, vertex->odx, vertex->ody, vertex->odz,
              vertex->x, vertex->y, vertex->z, vertex->nz,vertex->area) ;
    vertex->x += vertex->odx ; 
    vertex->y += vertex->ody ;
    vertex->z += vertex->odz ;
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
  char    fname[200] ;
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

  min_sse = starting_sse = mrisComputeSSE(mris, parms) ;

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
      mrisApplyGradient(mris, delta_t) ;
      mrisProjectSurface(mris) ;
      MRIScomputeMetricProperties(mris) ;
        
      sse = mrisComputeSSE(mris, parms) ;
      fprintf(fp2, "%f  %f  %f\n", delta_t, sse, predicted_sse) ;
      mrisProjectSurface(mris) ;
      sse = mrisComputeSSE(mris, parms) ;
      fprintf(fp, "%f  %f  %f\n", delta_t, sse, predicted_sse) ;
      fflush(fp) ;
      MRISrestoreOldPositions(mris) ;
    }
  }

  /* pick starting step size */
  min_delta = 0.0f ; /* to get rid of compiler warning */
  for (delta_t = min_dt ; delta_t < max_dt ; delta_t *= 10.0)
  {
    mrisApplyGradient(mris, delta_t) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    sse = mrisComputeSSE(mris, parms) ;
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
    mrisApplyGradient(mris, delta_t) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    min_sse = mrisComputeSSE(mris, parms) ;
    MRISrestoreOldPositions(mris) ;
  }

  delta_t = min_delta ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr,"grad=%2.3f, max_del=%2.3f, mean=%2.3f, max_dt=%2.1f, "
            "starting dt=%2.3f, min_dt=%2.3f\n", (float)grad, (float)max_delta, 
            mean_delta,(float)max_dt, (float)delta_t,min_dt) ;

  /* fit a quadratic form to it, and predict location of minimum */
  /* bracket the minimum by sampling on either side */
  N = 3 ;
  dt0 = min_delta - (min_delta/2) ;    dt2 = min_delta + (min_delta/2) ;
  mrisApplyGradient(mris, dt0) ;
  mrisProjectSurface(mris) ;           MRIScomputeMetricProperties(mris) ;
  sse0 = mrisComputeSSE(mris, parms) ; MRISrestoreOldPositions(mris) ;
  
  mrisApplyGradient(mris, dt2) ;
  mrisProjectSurface(mris) ;           MRIScomputeMetricProperties(mris) ;
  sse2 = mrisComputeSSE(mris, parms) ; MRISrestoreOldPositions(mris) ;
  
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
      fprintf(stderr,
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
        mrisApplyGradient(mris, new_min_delta) ;
        mrisProjectSurface(mris) ;           
        MRIScomputeMetricProperties(mris) ;
        sse = mrisComputeSSE(mris, parms) ;
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
    fprintf(stderr, "sses: %2.2f  ", sse_out[0]) ;
  mini = 0 ; min_sse = sse_out[mini] ; 
  for (i = 1 ; i < N ; i++)
  {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%2.2f  ", sse_out[i]) ;
    if (sse_out[i] < min_sse)
    {
      min_sse = sse_out[i] ;
      mini = i ;
    }
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "min %d (%2.3f)\n", mini, dt_in[mini]) ;
  mrisApplyGradient(mris, dt_in[mini]) ;
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
  char    fname[200] ;
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

  min_sse = starting_sse = mrisComputeSSE(mris, parms) ;

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
    mrisApplyGradient(mris, delta_t) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    sse = mrisComputeSSE(mris, parms) ;
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
    mrisApplyGradient(mris, delta_t) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    min_sse = mrisComputeSSE(mris, parms) ;
    MRISrestoreOldPositions(mris) ;
  }

  delta_t = min_delta ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr,"grad=%2.3f, max_del=%2.3f, mean=%2.3f, max_dt=%2.1f, "
            "starting dt=%2.3f, min_dt=%2.3f\n", (float)grad, (float)max_delta, 
            mean_delta,(float)max_dt, (float)delta_t,min_dt) ;

  /* now search for minimum in gradient direction */
  increasing = 1 ;
  total_delta = 0.0 ;
  min_sse = starting_sse ;
  while (!done)
  {
    mrisApplyGradient(mris, delta_t) ;
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    sse = mrisComputeSSE(mris, parms) ;
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

  if (Gdiag & DIAG_WRITE)
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
  int    k,i ;
  float  curv;
  char   fname[200], *cp, path[200], name[100] ;
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
  fp = fopen(fname,"wb");  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISwriteCurvature: could not open %s", 
                 fname)) ;

  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    curv = mris->vertices[k].curv ;
    i = nint(curv * 100.0) ;
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
int
MRISwriteDists(MRI_SURFACE *mris, char *sname)
{
  int    k,i ;
  float  dist ;
  char   fname[200], *cp, path[200] ;
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
  char   fname[200] ;
  
  MRISbuildFileName(mris, name, fname) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "writing area error file %s...", fname) ;

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
    fprintf(stderr, "done.\n") ;

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
    fprintf(stderr, "writing angular error file %s...", fname) ;

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
    fprintf(stderr, "done.\n") ;
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
int
MRISsampleStatVolume(MRI_SURFACE *mris, STAT_VOLUME *sv,int time_point,
                     int coords)
{
  VERTEX   *v ;
  int      vno, xv, yv, zv, width, height, depth ;
  Real     x, y, z, xt, yt, zt ;

  if (time_point >= sv->mri_pvals[0]->nframes)
    ErrorExit(ERROR_BADPARM, 
              "MRISsampleStatVolume: time point (%d) out of bounds [%d, %d]\n",
              time_point, 0, sv->mri_pvals[0]->nframes-1) ;
  width  = sv->mri_pvals[0]->width ;
  height  = sv->mri_pvals[0]->height ;
  depth  = sv->mri_pvals[0]->depth ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == 161)
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
    if (xv >= 0 && xv < width && yv >= 0 && yv <= height && zv>=0&&zv<=depth)
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
    fprintf(stderr, "writing out surface values to %s.\n", fname) ;

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
    fprintf(stderr, "avg = %2.3f, stdev = %2.3f, min = %2.3f, max = %2.3f\n",
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
  char  fname[200], *cp ;
  FILE *fp;
  double sum=0,sum2=0,max= -1000,min=1000;

#if 1
  MRISbuildFileName(mris, sname, fname) ;
#else
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
#endif
  cp = strrchr(fname, '.') ;
  if (!cp || *(cp+1) != 'w')
    strcat(fname, ".w") ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing surface values to %s.\n", fname) ;

  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorExit(ERROR_NOFILE, "Can't create file %s\n",fname) ;

  for (k=0,num=0;k<mris->nvertices;k++) 
    if (mris->vertices[k].val!=0) num++;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("num = %d\n",num);
  fwrite2(0,fp);
  fwrite3(num,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    if (mris->vertices[k].val!=0)
    {
      fwrite3(k,fp);
      f = mris->vertices[k].val;
      if (!finite(f))
        ErrorPrintf(ERROR_BADPARM, 
                    "MRISwriteValues(%s): val at vertex %d is not finite",
                    fname, k) ;

#if 0
      fwrite(&f,1,sizeof(float),fp);
#else
      fwriteFloat(f, fp) ;
#endif
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
  }
  fclose(fp);
  sum /= num;
  sum2 = sqrt(sum2/num-sum*sum);
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "avg = %2.3f, stdev = %2.3f, min = %2.3f, max = %2.3f\n",
            sum,sum2,min,max);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadValues(MRI_SURFACE *mris, char *fname)
{
  int i,k,num,ilat;
  float f;
  float lat;
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE,
                               "MRISreadValues: File %s not found\n",fname));
  fread2(&ilat,fp);
  lat = ilat/10.0;

  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].val=0;
  fread3(&num,fp);
  for (i=0;i<num;i++)
  {
    fread3(&k,fp);
    f = freadFloat(fp) ;
    /*    fread(&f,1,sizeof(float),fp);*/
    if (k>=mris->nvertices||k<0)
      printf("MRISreadValues: vertex index out of range: %d f=%f\n",k,f);
/*
    else if (mris->vertices[k].dist!=0)
      printf("MRISreadValues: subsample and data file mismatch\n");
*/
    else
    {
      mris->vertices[k].val = f;
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
MRISreadValuesScale(MRI_SURFACE *mris, char *fname)
{
  int i,k,num,ilat;
  float f;
  float lat;
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE,
                               "MRISreadValuesScale: File %s not found\n",fname));
  fread2(&ilat,fp);
  lat = ilat/10.0;

  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].val=0;
  fread3(&num,fp);
  for (i=0;i<num;i++)
    {
    fread3(&k,fp);
    f = freadFloat(fp) ;
    if (k>=mris->nvertices||k<0)
      printf("MRISreadValuesScale: vertex index out of range: %d f=%f\n",k,f);
    else
      {
      mris->vertices[k].val *= f;
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
MRISreadCanonicalCoordinates(MRI_SURFACE *mris, char *sname)
{
  int         nvertices, magic, version, ix, iy, iz, vno, n, nfaces, crap ;
  FILE        *fp ;
  VERTEX      *vertex ;
  char        fname[200], path[200], *cp ;
  float       d, x, y, z, r, theta, phi ;

  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
    cp = strstr(sname, "rh.") ;
    if (!cp)
      cp = strstr(sname, "lh.") ;
    if (!cp)   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", sname) ;
    else
      sprintf(fname, "%s/%s", path, sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadCanonicalCoordinates(%s): could not open file",
                 fname));

  fread3(&magic, fp) ;
  if (magic == NEW_VERSION_MAGIC_NUMBER) 
  {
    version = -1;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "new surface file format\n");
  }
  else 
  {
    rewind(fp);
    version = 0;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      printf("surfer: old surface file format\n");
  }
  fread3(&nvertices, fp);
  fread3(&nfaces, fp); nfaces *= 2 ;

  if ((nvertices != mris->nvertices) || (nfaces != mris->nfaces))
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM,
                "MRISreadCanonicalSurface(%s): nfaces %d (%d) or nvertices "
                "%d (%d) mismatch", fname,
                nfaces, mris->nfaces, nvertices, mris->nvertices)) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "reading %d vertices and %d faces.\n",nvertices,nfaces);

  for (vno = 0 ; vno < nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    fread2(&ix,fp);
    fread2(&iy,fp);
    fread2(&iz,fp);
    vertex->cx = ix/100.0;
    vertex->cy = iy/100.0;
    vertex->cz = iz/100.0;
    if (version == 0)  /* old surface format */
    {
      fread1(&crap,fp);
      for (n=0;n<vertex->num;n++)   /* skip over face data */
        fread3(&crap,fp);
    } 
#if 0
    else 
      vertex->num = 0;
#endif
  }
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  r = mris->radius = MRISaverageRadius(mris) ;
  r = mris->radius = (float)nint(mris->radius) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->cx ;
    y = vertex->cy ;
    z = vertex->cz ;
    theta = atan2(y/r, x/r) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = r*r-z*z ; if (d < 0.0) d = 0.0 ;
    phi = atan2(sqrt(d), z) ;
    vertex->theta = theta ; vertex->phi = phi ;
  }
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
  int         ix, iy, iz, k, i, j, npts ;
  FILE        *fp ;
  char        fname[200], path[200], *cp ;

  cp = strchr(pname, '/') ;
  if (cp)
    strcpy(fname, pname) ;    /* path already specified */
  else                        /* no path - use same as was used in MRISread */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/%s", path, pname) ;
  }
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,
                      "MRISreadPatch(%s): could not open file", fname));


  npts = freadInt(fp) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading patch %s with %d vertices (%2.1f%% of total)\n",
            pname, npts, 100.0f*(float)npts/(float)mris->nvertices) ;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].ripflag = TRUE;
  for (j=0;j<npts;j++)
  {
    i = freadInt(fp) ;
    if (i<0)
    {
      k = -i-1;
      mris->vertices[k].border = TRUE;
    } 
    else
    {
      k = i-1;
      mris->vertices[k].border = FALSE;
    }
    if (k >= mris->nvertices)
      ErrorExit(ERROR_BADFILE, 
                "MRISreadPatch: bad vertex # (%d) found in patch file", k) ;
    mris->vertices[k].ripflag = FALSE;
    fread2(&ix,fp);
    fread2(&iy,fp);
    fread2(&iz,fp);
    mris->vertices[k].x = ix/100.0;
    mris->vertices[k].y = iy/100.0;
    mris->vertices[k].z = iz/100.0;
    if (mris->vertices[k].x > mris->xhi) mris->xhi = mris->vertices[k].x;
    if (mris->vertices[k].x < mris->xlo) mris->xlo = mris->vertices[k].x;
    if (mris->vertices[k].y > mris->yhi) mris->yhi = mris->vertices[k].y;
    if (mris->vertices[k].y < mris->ylo) mris->ylo = mris->vertices[k].y;
    if (mris->vertices[k].z > mris->zhi) mris->zhi = mris->vertices[k].z;
    if (mris->vertices[k].z < mris->zlo) mris->zlo = mris->vertices[k].z;
    if (k == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stderr, "vertex %d read @ (%2.2f, %2.2f, %2.2f)\n",k,
              mris->vertices[k].x,mris->vertices[k].y,mris->vertices[k].z) ;
  }
  fclose(fp);
  MRISripFaces(mris);
  mris->patch = 1 ;
  mris->status = MRIS_CUT ;

  MRISremoveRipped(mris) ;
  MRISupdateSurface(mris) ;

  mrisComputeBoundaryNormals(mris);
  mrisSmoothBoundaryNormals(mris,10);
#if 0
  computed_shear_flag = FALSE;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].origripflag = mris->vertices[k].ripflag;
  flag2d = TRUE;
#endif

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

  type = mrisFileNameType(fname) ;
  if (type == MRIS_ASCII_QUADRANGLE_FILE)
    return(MRISwritePatchAscii(mris, fname)) ;
  else if (type == MRIS_GEO_TRIANGLE_FILE)
    return(MRISwriteGeo(mris, fname)) ;

  npts = 0;
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag) npts++;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "writing surface patch with npts=%d to %s\n",npts,fname);
  fp = fopen(fname,"w");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, 
                 "MRISwritePatch: can't create file %s\n",fname)) ;

  fwriteInt(npts, fp) ;
  for (k=0;k<mris->nvertices;k++)
  if (!mris->vertices[k].ripflag)
  {
    i = (mris->vertices[k].border)?-(k+1):k+1;
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
  float  nx, ny, nz, num ;
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
      v->nx = v->tdx ; v->ny = v->tdy ; v->nz = v->tdz ;
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
      fprintf(stderr, "vertex %d flattened @ (%2.2f, %2.2f, %2.2f)\n",vno,
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
  mrisComputeNormals(mris) ;
  for (k=0;k<mris->nvertices;k++)  
  if (!mris->vertices[k].ripflag)
  {
    v = &mris->vertices[k];
    x += v->x;
    y += v->y;
    z += v->z;
#if 0
    if (!FZERO(v->nx))
      fprintf(stderr, "vertex %d, normal = (%2.3f, %2.3f, %2.3f)\n",
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
    fprintf(stderr, "flatten: avg p = {%2.1f, %2.1f, %2.1f}\n",x,y,z);
    fprintf(stderr, "flatten: sum n = {%2.2f, %2.2f, %2.2f}\n",nx,ny,nz);
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
      fprintf(stderr, "flatten: norm n = {%2.2f, %2.2f, %2.2f}\n",nx,ny,nz);
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
      fprintf(stderr, "flatten: transformed n = {%2.1f, %2.1f, %2.1f}\n",
              nx,ny,nz);
  }
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].z = 0;
    if (k == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stderr, "vertex %d flattened @ (%2.2f, %2.2f, %2.2f)\n",k,
              mris->vertices[k].x,mris->vertices[k].y,mris->vertices[k].z) ;
  }

  mris->status = MRIS_PLANE ;
  MRIScomputeMetricProperties(mris) ;
  if (Gdiag & DIAG_SHOW && Gdiag_no >= 0)
  {
    int    n ;
    VERTEX *v, *vn ;
    
    v = &mris->vertices[Gdiag_no] ;
    fprintf(stderr, 
            "%d @ (%2.1f, %2.1f, %2.1f): area %2.3f, oa %2.3f, nz=%2.3f, vnum=%d\n",
            Gdiag_no, v->x, v->y, v->z, v->area, v->origarea, v->nz, v->vnum) ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      fprintf(stderr, 
              "%d @ (%2.1f, %2.1f, %2.1f): area %2.3f, oa %2.3f, nz=%2.3f, "
              "vnum=%d, d=%2.2f, od=%2.2f\n", 
              v->v[n], vn->x, vn->y, vn->z, vn->area, v->origarea, 
              vn->nz, vn->vnum, v->dist[n], v->dist_orig[n]) ;
    }
    fprintf(stderr, "\n") ;
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
    if (FZERO(V3_LEN(v_e1)))  /* happened to pick a parallel vector */
    {
      VECTOR_LOAD(v, vertex->ny, -vertex->nz, vertex->nx) ;
      V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    }

    if (FZERO(V3_LEN(v_e1)))  /* happened to pick a parallel vector */
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
          fprintf(stderr, "%d --> %d: curvature = %2.1f\n", vno, vertex->v[i],
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
    fprintf(stderr, "mean curvature range [%2.3f --> %2.3f]\n",Hmin, Hmax) ;
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
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ; b *= l_curv ;
    v->dx += b * v->nx ; v->dy += b * v->ny ; v->dz += b * v->nz ; 
    if (DIAG_VERBOSE_ON)
    {
      v->x += b * v->nx ; v->y += b * v->ny ; v->z += b * v->nz ; 
      vno-- ;
    }

    MatrixFree(&m_R) ; VectorFree(&v_Y) ; MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ; VectorFree(&v_e1) ; VectorFree(&v_e2) ; 
  VectorFree(&v_nbr) ; VectorFree(&v_A) ;
  return(NO_ERROR) ;
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
      fprintf(stderr, "Hdes=%2.3f, dH = %2.3f, tanh= %2.3f, dx=%2.3f\n",
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
      char fname[200] ;

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
      fprintf(stderr, "moving v %d by (%2.3f, %2.3f, %2.3f)\n",
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

  *pnfaces = 2*nfaces ;       /* two triangular faces per face */
  *pnvertices = nvertices ;
  *pnedges = nedges+nfaces ;  /* one additional edge added for each triangle */
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
      fprintf(stderr, "v %d: k1=%2.3f, k2=%2.3f, K=%2.3f, H=%2.3f\n",
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
    fprintf(stderr, "max H error=%2.3f at %d\n", max_error, vmax) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "total area = %2.3f\n", total_area);

  if (Gdiag & DIAG_SHOW && (nbad > 0))
    fprintf(stderr, "%d ill-conditioned points\n", nbad) ;
  MatrixFree(&m_eigen) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_c) ;
  VectorFree(&v_n) ;
  VectorFree(&v_yi) ;
  MatrixFree(&m_Q) ;
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

        Description
------------------------------------------------------*/
int
MRISinflateBrain(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     n_averages, n, write_iterations, niterations ;
  double  delta_t = 0.0, rms_height, desired_rms_height, sse ;

  write_iterations = parms->write_iterations ;
  n_averages = parms->n_averages ;

  if (Gdiag & DIAG_WRITE)
  {
    char fname[200] ;

    sprintf(fname, "%s.out", parms->base_name) ;
    parms->fp = fopen(fname, "w") ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  MRIScomputeMetricProperties(mris) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  parms->start_t = 0 ;
  niterations = parms->niterations ;
  MRISstoreMetricProperties(mris) ;
  rms_height = MRISrmsTPHeight(mris) ;
#if 0
  if (parms->scale > 0)
    desired_rms_height = parms->scale * rms_height ;
  else
#endif
    desired_rms_height = parms->desired_rms_height ;
  fprintf(stderr, "inflating to desired rms height = %2.3f\n", 
          desired_rms_height);

  /* write out initial surface */
  if ((parms->write_iterations > 0) && (Gdiag&DIAG_WRITE))
    mrisWriteSnapshot(mris, parms, 0) ;
  
  sse = mrisComputeSSE(mris, parms) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,"%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 
            0, 0.0f, (float)rms_height, n_averages) ;
  else
    fprintf(stderr, "\rstep %3.3d: RMS=%2.3f (target=%2.3f)   ", 0, 
            rms_height, desired_rms_height);
  if (Gdiag & DIAG_WRITE)
  {
    fprintf(parms->fp, 
            "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 
            0, 0.0f, (float)rms_height, n_averages) ;
    fflush(parms->fp) ;
  }

  MRISclearCurvature(mris) ;   /* curvature will be used to calculate sulc */
  for (n = 0 ; n < niterations ; n++)
  {
    mrisClearGradient(mris) ;
    mrisComputeDistanceTerm(mris, parms) ;
    mrisComputeSphereTerm(mris, parms->l_sphere, parms->a) ;
    mrisComputeExpansionTerm(mris, parms->l_expand) ;
    mrisComputeSpringTerm(mris, parms->l_spring) ;
    
    mrisAverageGradients(mris, n_averages) ;
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
      delta_t = mrisMomentumTimeStep(mris, parms->momentum, parms->dt, 
                                     parms->tol, parms->n_averages) ;
      break ;
    case INTEGRATE_ADAPTIVE:
      mrisAdaptiveTimeStep(mris, parms);
      break ;
    }
    mrisTrackTotalDistance(mris) ;  /* update sulc */
    MRIScomputeMetricProperties(mris) ; 
    if (!((n+1) % 5))   
    {
      sse = mrisComputeSSE(mris, parms) ;
      rms_height = MRISrmsTPHeight(mris) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, 
                "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 
                n+1,(float)delta_t, (float)rms_height, n_averages);
      else
        fprintf(stderr, "\rstep %3.3d: RMS=%2.3f (target=%2.3f)   ", 
                n+1, rms_height, desired_rms_height) ;
      if (Gdiag & DIAG_WRITE)
      {
        fprintf(parms->fp, 
                "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 
                n+1,(float)delta_t, (float)rms_height, n_averages);
        fflush(parms->fp) ;
      }
      
      if (rms_height < desired_rms_height)
        break ;
    }

    if (parms->scale > 0)
      MRISscaleBrainArea(mris) ;
    if ((parms->write_iterations > 0) &&
        !((n+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
      mrisWriteSnapshot(mris, parms, n+1) ;
  }

  fprintf(stderr, "\ninflation complete.\n") ;
  if (Gdiag & DIAG_WRITE)
    fclose(parms->fp) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisApplyGradient(MRI_SURFACE *mris, double dt)
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
static int
mrisClearGradient(MRI_SURFACE *mris)
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
  mrisComputeNormals(mris);              /* update vertex areas */
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
  char   line[200], *cp ;

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
  
  fprintf(stderr, "%d smoothed (%2.2f%%)\n",
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
      fprintf(stderr,"vertex %d normal curvature term: (%2.3f, %2.3f, %2.3f)\n"
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
    sx = nc*nx ;              /* move in normal direction */
    sy = nc*ny ;
    sz = nc*nz;
    
    vertex->dx += l_spring * sx ;
    vertex->dy += l_spring * sy ;
    vertex->dz += l_spring * sz ;
    if (vno == Gdiag_no)
      fprintf(stderr, "v %d spring normal term: (%2.3f, %2.3f, %2.3f)\n",
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
    sx -= nc*v->nx ;                   /* remove  normal component */
    sy -= nc*v->ny ;
    sz -= nc*v->nz;

    v->dx += l_spring * sx ;
    v->dy += l_spring * sy ;
    v->dz += l_spring * sz ;
    if (vno == Gdiag_no)
      fprintf(stderr, "v %d spring tangent term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
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
                         MRI *mri_smooth)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z, nx, ny, nz, dx, dy, dz ;
  Real    val0, xw, yw, zw, del, val_outside, val_inside, delI, delV ;

  if (FZERO(l_intensity))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ; y = v->y ; z = v->z ;
    
    MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolume(mri_brain, xw, yw, zw, &val0) ;

    nx = v->nx ; ny = v->ny ; nz = v->nz ;

    /* compute intensity gradient using smoothed volume */

    /* sample outward from surface */
    xw = x + nx ; yw = y + ny ; zw = z + nz ;
    MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_outside) ;

    /* sample inward from surface */
    xw = x - nx ; yw = y - ny ; zw = z - nz ;
    MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_inside) ;

    delV = v->val - val0 ;
    delI = (val_outside - val_inside) / 2.0 ;
#if 1
    if (!FZERO(delI))
      delI /= fabs(delI) ;
#endif
    del = l_intensity * delV * delI ;
    dx = nx * del ; dy = ny * del ; dz = nz * del ;

    v->dx += dx ;   
    v->dy += dy ;
    v->dz += dz ;
    if (vno == Gdiag_no)
      fprintf(stderr, "v %d intensity term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
  }
  

  return(NO_ERROR) ;
}
static int
mrisComputeIntensityGradientTerm(MRI_SURFACE*mris, double l_grad,
                                 MRI *mri_brain, MRI *mri_smooth)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z, nx, ny, nz ;
  Real    mag0, xw, yw, zw, del, mag_outside, mag_inside, delI, delV,
          dx, dy, dz, val_outside, val_inside, val_dist ;

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
    MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    mag0 = sqrt(dx*dx+dy*dy+dz*dz) ;
    nx = v->nx ; ny = v->ny ; nz = v->nz ;

    /* compute intensity gradient using smoothed volume */

    /* sample outward from surface */
    xw = x + nx ; yw = y + ny ; zw = z + nz ;
    MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    mag_outside = sqrt(dx*dx+dy*dy+dz*dz) ;
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_outside) ;

    /* sample inward from surface */
    xw = x - nx ; yw = y - ny ; zw = z - nz ;
    MRIworldToVoxel(mri_smooth, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz) ;
    mag_inside = sqrt(dx*dx+dy*dy+dz*dz) ;
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_inside) ;

    if (mag_outside > mag_inside)   /* gradients suggest edge is outwards */
      val_dist = val_outside - v->val ;
    else   /* gradients suggest edge is outwards */
      val_dist = val_inside - v->val ;

    /* diminish the effect of gradients that are in locations whos
       intensity values are far from the target.
       */
    val_dist = 2 / (2 + val_dist*val_dist) ;
    delV = 1.0f /*v->mean - mag0*/ ;
    delI = (mag_outside - mag_inside) / 2.0 ;
#if 1
    if (!FZERO(delI))
      delI /= fabs(delI) ;
#endif
    del = val_dist * l_grad * delV * delI ;
    dx = nx * del ; dy = ny * del ; dz = nz * del ;

    v->dx += dx ;   
    v->dy += dy ;
    v->dz += dz ;
    if (vno == Gdiag_no)
      fprintf(stderr, 
        "v %d intensity gradient term: (%2.3f, %2.3f, %2.3f) "
              "(mag = %2.1f, [%2.1f,%2.1f])\n",
              vno, v->dx, v->dy, v->dz, mag0, mag_inside, mag_outside) ;
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
  fprintf(stderr, "max radius = %2.1f\n", radius) ;
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
    x /= r ; y /= r ; z /= r ;
    r = (radius - r) / radius ;
    v->dx += r*l_sphere * x ;
    v->dy += r*l_sphere * y ;
    v->dz += r*l_sphere * z ;
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
    
    v->dx += l_spring * sx ;
    v->dy += l_spring * sy ;
    v->dz += l_spring * sz ;
    if (vno == Gdiag_no)
      fprintf(stderr, "v %d spring term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
  }
  

  return(NO_ERROR) ;
}
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
          Compute the folding of the surface.
------------------------------------------------------*/
static int
mrisComputeDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_y, *v_delta, *v_n ;
  float   l_dist, d0, dt, delta, nc, scale, norm, max_del ;
  VERTEX  *v, *vn ;
  int     vno, n, vnum, max_v, max_n ;

  l_dist = parms->l_dist ;

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
    fprintf(stderr, "distance scale = %2.3f\n", scale) ;
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
      fprintf(stderr, 
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
      if (FZERO(V3_LEN(v_y)))
        continue ;
      V3_NORMALIZE(v_y, v_y) ;   /* make it a unit vector */
      V3_SCALAR_MUL(v_y, delta, v_y) ;
      V3_ADD(v_y, v_delta, v_delta) ;
      if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
        fprintf(stderr, 
                "nbr %d (%6.6d) @ (%2.2f, %2.2f, %2.2f), "
                "d0 %2.2f, dt %2.2f, delta %2.3f\n\ty=%2.3f, %2.3f, %2.3f)\n",
                n, v->v[n], vn->x, vn->y, vn->z, d0, dt,
                delta, V3_X(v_y), V3_Y(v_y), V3_Z(v_y)) ;
    }

    V3_SCALAR_MUL(v_delta, norm, v_delta) ;

    if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stderr, 
                "total delta=(%2.3f, %2.3f, %2.3f)\n",
                V3_X(v_delta), V3_Y(v_delta), V3_Z(v_delta)) ;
    /* take out normal component */
    nc = V3_DOT(v_n, v_delta) ; V3_SCALAR_MUL(v_n, -nc, v_n) ;
    V3_ADD(v_delta, v_n, v_delta) ;

    v->dx += l_dist * V3_X(v_delta) ;
    v->dy += l_dist * V3_Y(v_delta) ;
    v->dz += l_dist * V3_Z(v_delta) ;
    if (vno == Gdiag_no)
      fprintf(stderr, "v %d, distance term: (%2.3f, %2.3f, %2.3f)\n",
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
  char  *cp, host_name[200] ;

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
    fprintf(fp, ", l_area=%2.4f", parms->l_area) ;
  if (!FZERO(parms->l_parea))
    fprintf(fp, ", l_parea=%2.4f", parms->l_parea) ;
  if (!FZERO(parms->l_narea))
    fprintf(fp, ", l_narea=%2.4f", parms->l_narea) ;
  if (!FZERO(parms->l_angle))
    fprintf(fp, ", l_angle=%2.4f", parms->l_angle) ;
  if (!FZERO(parms->l_corr))
    fprintf(fp, ", l_corr=%2.4f", parms->l_corr) ;
  if (!FZERO(parms->l_spring))
    fprintf(fp, ", l_spring=%2.4f", parms->l_spring) ;
  if (!FZERO(parms->l_tspring))
    fprintf(fp, ", l_tspring=%2.4f", parms->l_tspring) ;
  if (!FZERO(parms->l_nspring))
    fprintf(fp, ", l_nspring=%2.4f", parms->l_nspring) ;
  if (!FZERO(parms->l_dist))
    fprintf(fp, ", l_dist=%2.4f", parms->l_dist) ;
  if (!FZERO(parms->l_intensity))
    fprintf(fp, ", l_intensity=%2.4f", parms->l_intensity) ;
  if (!FZERO(parms->l_grad))
    fprintf(fp, ", l_grad=%2.4f", parms->l_grad) ;
  if (!FZERO(parms->l_sphere))
    fprintf(fp, ", l_sphere=%2.4f", parms->l_sphere) ;
  if (!FZERO(parms->l_expand))
    fprintf(fp, ", l_expand=%2.4f", parms->l_expand) ;
  if (!FZERO(parms->l_curv))
    fprintf(fp, ", l_curv=%2.4f", parms->l_curv) ;
  if (!FZERO(parms->l_boundary))
    fprintf(fp, ", l_boundary=%2.4f", parms->l_boundary) ;
  if (!FZERO(parms->l_neg))
    fprintf(fp, ", l_neg=%2.4f", parms->l_neg) ;
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
mrisWriteSnapshot(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int t)
{
  char fname[200], path[200], base_name[200] ;

  FileNamePath(mris->fname, path) ;
  sprintf(base_name, "%s/%s.%s", path, 
          mris->hemisphere == LEFT_HEMISPHERE ? "lh":"rh", parms->base_name);
  sprintf(fname, "%s%4.4d", base_name, t) ;
#if 1
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing %s...", fname) ;
  if (mris->patch)
    MRISwritePatch(mris, fname) ;
  else
    MRISwrite(mris, fname) ;
  
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "done.\n") ;
    fflush(stderr) ;
  }
#endif

  if (mris->status == MRIS_PARAMETERIZED_SPHERE && DIAG_VERBOSE_ON)
  {
    MRI_SP *mrisp = (MRI_SP *)mris->vp ;

    sprintf(fname, "%s%4.4d.hipl", parms->base_name, t) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "writing %s\n", fname) ;
    MRIStoParameterization(mris, mrisp, 1, 0) ;
    MRISPwrite(mrisp, fname) ;
  }
  if (!FZERO(parms->l_area) && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "%s%4.4d.area_error", base_name, t) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, " %s...", fname) ;
    MRISwriteAreaError(mris, fname) ;
  }

  if (!FZERO(parms->l_corr) && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "%s%4.4d.angle_error", base_name, t) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, " %s...", fname) ;
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
      delta = dist_scale*v->dist[n] - v->dist_orig[n] ;

if (fabs(delta) > max_del)
{
  max_del = delta ;
  max_v = vno ;
  max_n = n ;
}
      v_sse += delta*delta ;
      if (!finite(delta) || !finite(v_sse))
        DiagBreak() ;
    }
    sse_dist += v_sse ;
    if (!finite(sse_dist) || !finite(v_sse))
      DiagBreak() ;
  }

  /*fprintf(stderr, "max_del = %f at v %d, n %d\n", max_del, max_v, max_n) ;*/
  return(sse_dist) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeSpringEnergy(MRI_SURFACE *mris)
{
  int     vno, n, nvertices ;
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

  for (sse_spring = 0.0, nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nvertices++ ;
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
    area_scale = mris->total_area / mris->orig_area ;
#else
  area_scale = 1.0f ;
#endif

  l_angle = parms->l_angle ;
  l_area = parms->l_area ;
  l_parea = parms->l_parea ;
  if (!FZERO(parms->l_narea))
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
    orig_area = area_scale * face->orig_area ; area = face->area ;
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
  double  orig_area, area, delta, area_scale, scale, l_narea, ratio ;

#if METRIC_SCALE
  if (mris->patch || 
      (mris->status != MRIS_SPHERE && 
       mris->status != MRIS_PARAMETERIZED_SPHERE))
    area_scale = 1.0f ;
  else
    area_scale = mris->total_area / mris->orig_area ;
#else
  area_scale = 1.0f ;
#endif

  l_narea = parms->l_narea ;

  if (FZERO(l_narea))
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
    orig_area = area_scale * face->orig_area ; area = face->area ;
    if (!FZERO(orig_area))
      ratio = area / orig_area ;
    else
      ratio = 0.0f ;

    if (ratio > MAX_NEG_RATIO)
      ratio = MAX_NEG_RATIO ;
    else if (ratio < -MAX_NEG_RATIO)
      ratio = -MAX_NEG_RATIO ;
#if 0
    scale = l_narea * (1 - (1/(1.0+exp(-NEG_AREA_K*ratio)))) ;
#else
    scale = l_narea / (1.0+exp(NEG_AREA_K*ratio)) ;
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
  sse = mrisComputeSSE(mris, parms) ;
#endif
#if 0
  sse /= (float)mrisValidVertices(mris) ;
  sse = sqrt(sse) ;
#endif
  if (FZERO(parms->l_corr) && FZERO(parms->l_pcorr))
    fprintf(fp, "%3.3d: dt: %2.1f, sse: %2.3f (%2.3f, %2.1f, %2.3f), "
            "neg: %d (%%%2.3f), avgs: %d\n", 
            parms->t, dt, sse, area_rms, (float)DEGREES(angle_rms), dist_rms,
            negative, 100.0*mris->neg_area/(mris->neg_area+mris->total_area),
            parms->n_averages);
  else
    fprintf(fp, "%3.3d: dt: %2.3f, sse: %2.1f (%2.3f, %2.1f, %2.3f, %2.3f), "
            "neg: %d (%%%2.2f), avgs: %d\n", 
            parms->t, dt, sse, area_rms, (float)DEGREES(angle_rms), dist_rms,
            corr_rms, negative, 
            100.0*mris->neg_area/(mris->neg_area+mris->total_area),
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
  char    fname[200] ;
  int     vno, nvertices, nfaces, magic, version, tmp, ix, iy, iz, n, type ;
  VERTEX  *vertex ;
  FILE    *fp ;

  type = mrisFileNameType(fname) ;
  MRISbuildFileName(mris, name, fname) ;
  if (type == MRIS_GEO_FILE)
    return(mrisReadGeoFilePositions(mris, fname)) ;
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
      fprintf(stderr, "new surface file format\n");
  }
  else if (magic == TRIANGLE_FILE_MAGIC_NUMBER)
  {
    fclose(fp) ;
    if (mrisReadTriangleFilePositions(mris,fname)  != NO_ERROR)
      ErrorReturn(Gerror, (Gerror, "mrisReadTriangleFile failed.\n")) ;
    version = -2 ;
  }
  else 
  {
    rewind(fp);
    version = 0;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      printf("surfer: old surface file format\n");
  }

  if (version >= -1)
  {
    fread3(&nvertices, fp);
    fread3(&nfaces, fp); nfaces *= 2 ;
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "reading %d vertices and %d faces.\n",nvertices,nfaces);
    
    if (nvertices != mris->nvertices || nfaces != mris->nfaces)
      ErrorExit(ERROR_BADFILE, 
                "MRISreadVertexPositions(%s): surfaces differ.\n",
                fname) ;

    for (vno = 0 ; vno < nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      fread2(&ix,fp);
      fread2(&iy,fp);
      fread2(&iz,fp);
      if (!vertex->ripflag)
      {
        vertex->x = ix/100.0;
        vertex->y = iy/100.0;
        vertex->z = iz/100.0;
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
  mrisComputeNormals(mris) ;
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
    fprintf(stderr, "computed %.1f, analytic %.1f\n",
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
  MRISreadVertexPositions(mris, sname) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeTriangleProperties(mris) ;
  MRISstoreMetricProperties(mris) ;
  mris->status = old_status ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeTriangleProperties(mris) ;
  mrisOrientSurface(mris) ;
        
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
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int   
mrisFileNameType(char *fname)
{
  int   type ;
  char  *dot, ext[200] ;

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
    type = MRIS_ASCII_FILE ;
  else if (!strcmp(ext, "GEO"))
    type = MRIS_GEO_FILE ;
  else
    type = MRIS_BINARY_FILE ;

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
MRISwriteGeo(MRI_SURFACE *mris, char *fname)
{
  int     vno, fno, n, actual_vno, toggle, nfaces, nvertices ;
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
    v->marked = actual_vno++ ;
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
      fprintf(fp, "%d ", mris->vertices[face->v[n]].marked+1) ; /* 1-based */

    /* swap order on output to conform to movie.byu convention */
    fprintf(fp, "%d ",mris->vertices[face->v[VERTICES_PER_FACE-1]].marked+1);
    fprintf(fp, "-%d\n",mris->vertices[face->v[VERTICES_PER_FACE-2]].marked+1);
  }

  MRISclearMarks(mris) ;
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

  type = mrisFileNameType(fname) ;
#if 0
  if (type == MRIS_ASCII_QUADRANGLE_FILE)
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
  fprintf(stderr, "nvertices=%d (valid=%d) nfaces=%d\n", nvertices, 
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
  char    line[200], *cp ;
  int     vno, fno, n, nvertices, nfaces, patch ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
              (ERROR_NOFILE, 
               "MRISreadAsciiFile: could not open file %s",fname));

  patch = 0 ;
  cp = fgetl(line, 100, fp) ;
  sscanf(cp, "%d %d\n", &nvertices, &nfaces) ;
  mris = MRISalloc(nvertices, nfaces) ;
  mris->type = MRIS_ASCII_QUADRANGLE_FILE ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    fscanf(fp, "%f  %f  %f  %d\n", &v->x, &v->y, &v->z, &v->ripflag) ;
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
    fscanf(fp, "%d\n", &face->ripflag) ;
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
  char    line[200], *cp ;
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
    fprintf(stderr, "max gradient magnitude = %2.5f\n", max_mag) ;

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
    fprintf(stderr, "max gradient magnitude = %2.5f\n", max_mag) ;

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
    fprintf(stderr, "curvature mean = %2.3f\n", mean) ;

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
    fprintf(stderr, "curvature mean = %2.3f, std = %2.3f\n", mean, std) ;

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
------------------------------------------------------*/
static int
mrisTrackTotalDistance(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;
  float  nc ;

  for (vno = 1 ; vno < mris->nvertices ; vno++)
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

  for (vno = 1 ; vno < mris->nvertices ; vno++)
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

  steps = mrisIntegrate(mris, &parms, 0) ;
  old_parms->start_t += steps ;
  mris->status = old_status ;
  if (Gdiag & DIAG_WRITE)
    fprintf(old_parms->fp, "rigid alignment complete, sse = %2.3f\n",
            mrisComputeSSE(mris, old_parms)) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    float area_rms, angle_rms, curv_rms, dist_rms, corr_rms, rms ;

    rms = 
      mrisComputeError(mris, &parms,&area_rms,&angle_rms,&curv_rms,&dist_rms,
                       &corr_rms);
    fprintf(stderr, "rms = %2.3f, corr_rms = %2.3f ", rms, corr_rms) ;
    rms = 
      mrisComputeError(mris, old_parms,&area_rms,&angle_rms,&curv_rms,
                       &dist_rms, &corr_rms);
    fprintf(stderr, "(%2.3f, %2.3f)\n", rms, corr_rms) ;
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
    min_sse = mrisComputeCorrelationError(mris, parms, 0) ;
    delta = 2*degrees / (float)nangles ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "scanning %2.2f degree nbhd, min sse = %2.2f\n", 
              (float)DEGREES(degrees), (float)min_sse) ;
    for (alpha = -degrees ; alpha <= degrees ; alpha += delta)
    {
      for (beta = -degrees ; beta <= degrees ; beta += delta)
      {
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, "\r(%+2.2f, %+2.2f, %+2.2f), "
                  "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                  (float)DEGREES(alpha), (float)DEGREES(beta), (float)
                  DEGREES(-degrees), (float)DEGREES(mina), 
                  (float)DEGREES(minb), (float)DEGREES(ming),(float)min_sse);

        for (gamma = -degrees ; gamma <= degrees ; gamma += delta)
        {
          MRISsaveVertexPositions(mris, TMP_VERTICES) ;
          MRISrotate(mris, mris, alpha, beta, gamma) ;
          sse = mrisComputeCorrelationError(mris, parms, 0) ;
          MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
          if (sse < min_sse)
          {
            mina = alpha ; minb = beta ; ming = gamma ;
            min_sse = sse ;
          }
#if 0
          if (Gdiag & DIAG_SHOW)
            fprintf(stderr, "\r(%+2.2f, %+2.2f, %+2.2f), "
                    "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                    (float)DEGREES(alpha), (float)DEGREES(beta), (float)
                    DEGREES(gamma), (float)DEGREES(mina), 
                    (float)DEGREES(minb), (float)DEGREES(ming),(float)min_sse);
#endif
        }
      }
    }
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "\n") ;
    if (!FZERO(mina) || !FZERO(minb) || !FZERO(ming))
    {
      MRISrotate(mris, mris, mina, minb, ming) ;
      sse = mrisComputeCorrelationError(mris, parms, 0) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "min sse = %2.2f at (%2.2f, %2.2f, %2.2f)\n",
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
    v->curv = MAX(fabs(v->k1), fabs(v->k2)) ;
    if (v->curv > kmax)
      kmax = v->curv ;
    if (v->curv < kmin)
      kmin = v->curv ;
    if (v->curv < 0)
      DiagBreak() ;
  }

  fprintf(stderr, "kmin = %2.2f, kmax = %2.2f\n", kmin, kmax) ;
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
    v->curv = MIN(v->k1, v->k2) ;
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
  MRIworldToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;
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
  MRIworldToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;
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

#include "timer.h"
int
MRISpositionSurface(MRI_SURFACE *mris, MRI *mri_brain, MRI *mri_smooth,
                    INTEGRATION_PARMS *parms)
{
  char   *cp /*, *nstr*/ ;
  int    avgs, niterations, n, write_iterations ;
  double sse, delta_t = 0.0, rms, dt, l_intensity ;
  MHT    *mht = NULL ;
  struct timeb  then ; int msec ;

  TimerStart(&then) ;
  parms->mri_brain = mri_brain ;
  parms->mri_smooth = mri_smooth ;
  niterations = parms->niterations ;
  write_iterations = parms->write_iterations ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[200] ;

    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    parms->fp = fopen(fname, "w") ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;
  cp = getenv("AVERAGE_VALS") ;
  if (!cp)
    cp = "1" ;
  avgs = atoi(cp) ;
  fflush(stdout) ;
#if 0
  fprintf(stderr, "averaging target values %d times\n", avgs) ;
  if (Gdiag & DIAG_WRITE)
    fprintf(fp, "averaging target values %d times\n", avgs) ;
#endif
  MRISaverageVals(mris, avgs) ;

  mrisComputeNormals(mris) ;

  mrisClearDistances(mris) ;

  /* write out initial surface */
  if ((parms->write_iterations > 0) && (Gdiag&DIAG_WRITE))
    mrisWriteSnapshot(mris, parms, 0) ;

  avgs = parms->n_averages ;
  rms = mrisRmsValError(mris, mri_brain) ;
  sse = mrisComputeSSE(mris, parms) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%3.3d: dt: %2.4f, sse=%2.1f, rms=%2.2f\n", 
            0, 0.0f, (float)sse, (float)rms);

  if (Gdiag & DIAG_WRITE)
  {
    fprintf(parms->fp, "%3.3d: dt: %2.4f, sse=%2.1f, rms=%2.2f\n", 
            0, 0.0f, (float)sse, (float)rms);
    fflush(parms->fp) ;
  }

  MRISclearCurvature(mris) ;   /* curvature will be used to calculate sulc */
  dt = parms->dt ; l_intensity = parms->l_intensity ;
  for (n = 0 ; n < niterations ; n++)
  {
    if (n == niterations-5)
    {
      dt = parms->dt / 4.0f ;  /* take some small steps at the end */
      /*      l_intensity = 0.0f ;*/
    }
    if (n == niterations)
    {
      parms->l_spring *= 2.0f ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "setting l_spring to %2.2f\n", parms->l_spring) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp, "setting l_spring to %2.2f\n", parms->l_spring) ;
    }
    if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
      mht = MHTfillTable(mris, mht) ;
    mrisClearGradient(mris) ;
    mrisComputeIntensityTerm(mris, l_intensity, mri_brain, mri_smooth);
    mrisComputeIntensityGradientTerm(mris, parms->l_grad,mri_brain,mri_smooth);
    mrisAverageGradients(mris, avgs) ;
    mrisComputeSpringTerm(mris, parms->l_spring) ;

    /* smoothness terms */
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
      delta_t = mrisMomentumTimeStep(mris, parms->momentum, parms->dt, 
                                     parms->tol, avgs) ;
      break ;
    case INTEGRATE_ADAPTIVE:
      mrisAdaptiveTimeStep(mris, parms);
      break ;
    }
#else
    /*    mrisClipMomentumGradient(mris, 0.2f) ;*/
    delta_t = mrisAsynchronousTimeStep(mris, parms->momentum, dt,mht) ;
    if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
      MHTcheckFaces(mris, mht) ;
#endif
    mrisTrackTotalDistance(mris) ;  /* update thickness measure */
    MRIScomputeMetricProperties(mris) ; 
    rms = mrisRmsValError(mris, mri_brain) ;
    sse = mrisComputeSSE(mris, parms) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "%3.3d: dt: %2.4f, sse=%2.1f, rms=%2.2f\n", 
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
  }

  fprintf(stderr, "\nsurface positioning complete.\n") ;
  if (Gdiag & DIAG_SHOW)
  {
    msec = TimerStop(&then) ;
    fprintf(stderr,"positioning took %2.1f minutes\n", 
            (float)msec/(60*1000.0f));
  }
  if (Gdiag & DIAG_WRITE)
    fclose(parms->fp) ;

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
        MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
        xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;  /* voxel coordinate */
        if (on)
          MRIset_bit(mri, xv, yv, zv) ;                 /* mark it filled */
        else
          MRIclear_bit(mri, xv, yv, zv) ;               /* mark it empty */
      }
      /* compute last point on line */
      t1 = 1.0f ;
      x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
      MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
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
      MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;  /* voxel coordinate */
      if (on)
        MRIset_bit(mri, xv, yv, zv) ;                 /* mark it filled */
      else
        MRIclear_bit(mri, xv, yv, zv) ;               /* mark it empty */
    }
    /* compute last point on line */
    t1 = 1.0f ;
    x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
    MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;   /* volume coordinate */
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
#define MIN_GRAY              55
#define MAX_GRAY              100
#define WHALF                 (5-1)/2
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
  mrisComputeNormals(mris);

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
    if (max_out_dist > MAX_THICKNESS)
      v->curv = MAX_THICKNESS ;   /* can't compute it properly */
    else
      v->curv = max_out_dist ;
  }
  return(NO_ERROR) ;
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
    if (mrisNormalDirectionTriangleIntersection(mris, v, mht, &dist))
      return(dist) ;
  }
  
  return(0.0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           See if the line in the normal direction passing
           through vertex v intersects any of the triangles at the
           given location.
------------------------------------------------------*/
#include "tritri.h"
static int
mrisNormalDirectionTriangleIntersection(MRI_SURFACE *mris, VERTEX *v, 
                                        MHT *mht, double *pdist)
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
    return(0) ;

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
    return(0) ;
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
        found = 1 ;
        *pdist = min_dist = dist ;
      }
    }
  }
  return(found) ;
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
      if (v->marked == 2)
      {
        x = v->x ; y = v->y ; z = v->z ;
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
    if (Gdiag & DIAG_SHOW)
      printf("%d: %d vertices marked\n", i,nmarked);
    if (!nmarked)
      break ;
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

  for (vno = 1 ; vno < mris->nvertices ; vno++)
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

  for (vno = 1 ; vno < mris->nvertices ; vno++)
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
MRIScomputeWhiteSurfaceValues(MRI_SURFACE *mris,MRI *mri_brain,
                              MRI *mri_smooth, MRI *mri_wm)
{
  Real    val, x, y, z, min_val, xw, yw, zw, dx, dy, dz, mag, max_mag ;
  int     total_vertices, vno, nmissing = 0 ;
  float   mean_white, dist ;
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
    for (dist = -3.0f ; dist < 10.0f ; dist += STEP_SIZE)
    {
      x = v->x + v->nx*dist ; y = v->y + v->ny*dist ; z = v->z + v->nz*dist ; 
      MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
      MRIsampleVolume(mri_brain, xw, yw, zw, &val) ;
      if (val < 100 && val > 90)  /* in right range */
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
      fprintf(stderr, "v %d, target value = %2.1f, mag = %2.1f\n",
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
          fprintf(stderr, "mean gray (%2.1f) > mean white (%2.1f) at v %d!\n",
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
      fprintf(stderr, "v %d, target value = %2.1f, mag = %2.1f\n",
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
      MRIworldToVoxel(mri_brain, x, y, z, &xw, &yw, &zw) ;
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
      fprintf(stderr, "v %d, target value = %2.1f, mag = %2.1f\n",
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
  fprintf(stderr, 
          "vertex #%d @ (%2.2f, %2.2f, %2.2f), n = (%2.2f, %2.2f, %2.2f) "
          "(%2.2f, %2.2f, %2.2f), val=%2.1f\n",
          vno, v->x, v->y, v->z, v->nx, v->ny, v->nz, v->dx, v->dy, v->dz,
          v->val) ;

  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    dx = vn->x - v->x ; dy = vn->y - v->y ; dz = vn->z - v->z ;
    d = sqrt(dx*dx + dy*dy + dz*dz) ;
    fprintf(stderr, 
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
    if (v->ripflag)
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
  fprintf(stderr, 
      "mean orig = %2.3f mm (%%%2.2f), final = %2.3f mm (%%%2.2f)\n",
          mean_orig_error, 100.0*pct_orig, mean_error, 100.0*pct) ;
  fprintf(stderr, "signed mean orig error = %2.3f, final mean error = %2.3f\n",
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
        MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;
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
        MRIworldToVoxel(mri, x, y, z, &x, &y, &z) ;
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

    v->mean_imag = (v->mean*total_dof + v->val*new_dof) / ndof ;
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
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ; y = v->y ; z = v->z ;
    
    MRIworldToVoxel(parms->mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolume(parms->mri_brain, xw, yw, zw, &val0) ;

    del0 = v->val - val0 ;
    sse += (del0 * del0) ;
    if (vno == Gdiag_no)
      fprintf(stderr, "v %d spring term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
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
    
    MRIworldToVoxel(parms->mri_smooth, x, y, z, &xw, &yw, &zw) ;
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
mrisComputeSphereError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     vno ;
  VERTEX  *v ;
  double  del, sse, x, y, z, x0, y0, z0, r, r0 ;

  if (FZERO(parms->l_sphere))
    return(0.0f) ;

#if 0
  r0 = MRISmaxRadius(mris) ;
#else
  r0 = parms->a ;
#endif
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
    if (vno == Gdiag_no)
      fprintf(stderr, "v %d spring term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
  }
  
  return(parms->l_sphere * sse) ;
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
  fprintf(stderr, "mean dist = %2.3f, rms error = %2.2f%%\n",
          total_mean, 100.0*total_mean_error) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
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
    fprintf(stderr, "writing surface file %s, created by %s on %s.\n",
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
  fclose(fp);
  return(NO_ERROR) ;
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
  char        line[200] ;
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
  
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr,"surface %s: %d vertices and %d faces.\n",
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
mrisReadTriangleFile(char *fname)
{
  VERTEX      *v ;
  int         nvertices, nfaces, magic, vno, fno, n ;
  char        line[200] ;
  FILE        *fp ;
  MRI_SURFACE *mris ;

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
    fprintf(stderr,"surface %s: %d vertices and %d faces.\n",
            fname, nvertices,nfaces);
  
  mris = MRISalloc(nvertices, nfaces) ;
  mris->type = MRIS_TRIANGULAR_SURFACE ;
  
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->x = freadFloat(fp);
    v->y = freadFloat(fp);
    v->z = freadFloat(fp);
#if 0
    v->label = NO_LABEL ;
#endif
    v->num = 0;   /* will figure it out */
  }
  
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)  
      mris->faces[fno].v[n] = freadInt(fp);
    
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      mris->vertices[mris->faces[fno].v[n]].num++;
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
  char   path[200], *slash, *dot ;
  
  slash = strchr(sname, '/') ;
  if (!slash)              /* no path - use same one as mris was read from */
  {
    dot = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (dot && ((dot-slash) == 3) && (*(dot-1) == 'h') &&
        (*(dot-2) == 'l' || *(dot-2) == 'r'))
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
MRISsoapBubbleVals(MRI_SURFACE *mris, int navgs)
{
  int     vno, n, i, cmpt, nmarked ;
  VERTEX  *v, *vn ;
  double  mean ;

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
    if (!nmarked)
      break ;
  }
  
  fprintf(stderr, "\n") ;
  return(NO_ERROR) ;
}
static int
mrisRemoveTriangleLinks(MRI_SURFACE *mris)
{
  int    fno ;
  FACE   *f ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "removing non-quadrangular links.\n") ;

  for (fno = 0 ; fno < mris->nfaces ; fno += 2)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    mrisRemoveVertexLink(mris, f->v[1], f->v[2]) ;
    mrisRemoveVertexLink(mris, f->v[2], f->v[1]) ;
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
#if 0
static int mrisDivideEdge(MRI_SURFACE *mris, int vno1, int vno2) ;
static int mrisDivideFace(MRI_SURFACE *mris, int fno, int vno1, int vno2, 
                          int vnew_no) ;

static int   mrisAddVertices(MRI_SURFACE *mris, double nsigma) ;
static double mrisComputeVertexSpacingStats(MRI_SURFACE *mris, double *psigma);
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
  mean = mrisComputeVertexSpacingStats(mris, &sigma) ;
#else
  mean = mrisComputeVertexNormalSpacingStats(mris, &sigma) ;
  thresh *= thresh ;   /* make it squared so we don't need sqrts later */
#endif
  thresh = mean + sigma * nsigma ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "dividing edges more than %2.2f mm long.\n", thresh) ;
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
    fprintf(stderr, "%d vertices added: # of vertices=%d, # of faces=%d.\n", 
            nadded, mris->nvertices, mris->nfaces) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stderr, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
  }
  return(mean) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeVertexSpacingStats(MRI_SURFACE *mris, double *psigma)
{
  double   total_dist, mean, var, nv, dist, sigma ;
  int      vno, n ;
  VERTEX   *v, *vn ;
  
  for (var = nv = total_dist = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
   {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      nv++ ;
      dist = sqrt(SQR(vn->x - v->x) + SQR(vn->y - v->y) + SQR(vn->z - v->z));
      total_dist += dist ;
      var += dist*dist ;
    }
  }
  mean = total_dist / nv ;
  *psigma = sigma = sqrt(var / nv - mean*mean) ;
  return(mean) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_NEIGHBORS 10
#define MAX_FACES     10
static int
mrisDivideEdge(MRI_SURFACE *mris, int vno1, int vno2)
{
  VERTEX   *v1, *v2, *vnew ;
  int      vnew_no, n, m, fno, n1, n2, flist[100] ;
  FACE     *face ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "dividing edge %d --> %d\n", vno1, vno2) ;

  if (vno1 == 2241 && vno2 == 2237)
    DiagBreak() ;
  if (vno1 == 2241 || vno2 == 2237 ||
      vno2 == 2241 || vno1 == 2237)
    DiagBreak() ;
  if (vno1 == Gdiag_no || vno2 == Gdiag_no || mris->nvertices == Gdiag_no)
    DiagBreak() ;
  v1 = &mris->vertices[vno1] ;
  v2 = &mris->vertices[vno2] ;
  if (v1->ripflag || v2->ripflag || mris->nvertices >= mris->max_vertices ||
      mris->nfaces >= (mris->max_faces-1))
    return(ERROR_NO_MEMORY) ;

  /* check to make sure these vertices or the faces they */
  if (v1->vnum >= MAX_NEIGHBORS || v2->vnum >= MAX_NEIGHBORS ||
      v1->num  >= MAX_FACES ||     v2->num >= MAX_FACES)
    return(ERROR_NO_MEMORY) ;


  /* add 1 new vertex, 2 new faces, and 2 new edges */
  vnew_no = mris->nvertices++ ;
  vnew = &mris->vertices[vnew_no] ;
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
          fprintf(stderr, " face %d shared.\n", fno) ;
        flist[vnew->num++] = fno ;
        vnew->vnum++ ;
      }
  }

  /* will be part of two new faces also */
  flist[vnew->num++] = mris->nfaces ; flist[vnew->num++] = mris->nfaces+1 ;
  vnew->num = 4 ; vnew->vnum = 4 ;
  vnew->f = (int *)calloc(vnew->num, sizeof(int)) ;
  if (!vnew->f)
    ErrorExit(ERROR_NOMEMORY, "could not allocate %dth face list.\n", vnew_no);
  vnew->n = (int *)calloc(vnew->num, sizeof(int)) ;
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
          fprintf(stderr, "dividing face %d along edge %d-->%d.\n", 
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
        vnew->n[fno] = n ;
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
  }
  if (0 && Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "%d edges and %d faces.\n", vnew->vnum, vnew->num) ;

  if (!vnew->vnum || !v1->vnum || !v2->vnum)
  {
    fprintf(stderr, "empty vertex!\n") ;
    DiagBreak() ;
  }
  if (vnew->vnum != 4 || vnew->num != 4)
    DiagBreak() ;
  if (mris->nfaces > 12452)
    MHTcheckFaces(mris, mht) ;
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
  int    fnew_no, n, vno3, flist[500], vlist[500], nlist[500] ;

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
      vnew->n[vnew->num++] = n ;
    }
    else if (f1->v[n] != vno1)
      vno3 = f1->v[n] ;
  }
  v3 = &mris->vertices[vno3] ;

  /* now construct f2 */

  /*  replace vno1 with vnew_no in f2 */
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    if (f2->v[n] == vno1)   /* replace it with vnew */
    {
      f2->v[n] = vnew_no ;
      vnew->f[vnew->num] = fnew_no ;
      vnew->n[vnew->num++] = n ;
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
  v3->n = (int *)calloc(v3->num+1,sizeof(int));
  if (!v3->n)
    ErrorExit(ERROR_NO_MEMORY, "mrisDivideFace: could not allocate %d nbrs",
              v3->n) ;
  memmove(v3->f, flist, v3->num*sizeof(v3->f[0])) ;
  memmove(v3->n, nlist, v3->num*sizeof(v3->n[0])) ;
  memmove(v3->v, vlist, v3->vnum*sizeof(v3->v[0])) ;
  v3->v[v3->vnum++] = vnew_no ;
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
    fprintf(stderr, "face %d: (%d, %d, %d)\n",
            fno, f1->v[0], f1->v[1], f1->v[2]);
    fprintf(stderr, "face %d: (%d, %d, %d)\n",
            fnew_no, f2->v[0], f2->v[1], f2->v[2]);
  }
  return(NO_ERROR) ;
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
