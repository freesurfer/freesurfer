#include <stdio.h>
#include <math.h>

#include "diag.h"
#include "error.h"
#include "utils.h"
#include "macros.h"
#include "fio.h"
#include "mrisurf.h"
#include "matrix.h"
#include "proto.h"
#include "stats.h"
#include "timer.h"
#include "const.h"

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
#define MAX_NEG_AREA_PCT    0.02f

#define MRIS_BINARY_FILE    0
#define MRIS_ASCII_FILE     1

/*------------------------ STATIC PROTOTYPES -------------------------*/

static int   mrisCheck(MRI_SURFACE *mris) ;
static int   mrisFileNameType(char *fname) ;
static int   mrisComputeNormals(MRI_SURFACE *mris) ;
static int   mrisFindNeighbors(MRI_SURFACE *mris) ;
static int   mrisRemoveRipped(MRI_SURFACE *mris) ;
static void  mrisNormalize(float v[3]) ;
static float mrisTriangleArea(MRIS *mris, int fac, int n) ;
static int   mrisNormalFace(MRIS *mris, int fac,int n,float norm[]) ;
static int   mrisReadTransform(MRIS *mris, char *mris_fname) ;
static MRI_SURFACE *mrisReadAsciiFile(char *fname) ;

static int   mrisReadBinaryAreas(MRI_SURFACE *mris, char *mris_fname) ;
/*static int   mrisReadFieldsign(MRI_SURFACE *mris, char *fname) ;*/
static double mrisComputeSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double mrisComputeSpringEnergy(MRI_SURFACE *mris) ;
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
static double mrisAdaptiveTimeStep(MRI_SURFACE *mris,INTEGRATION_PARMS *parms);
static int   mrisOrientEllipsoid(MRI_SURFACE *mris) ;
static int   mrisOrientPlane(MRI_SURFACE *mris) ;
#if AVERAGE_AREAS
static int   mrisAverageAreas(MRI_SURFACE *mris, int num_avgs, int which) ;
#endif
static int   mrisRipFaces(MRI_SURFACE *mris) ;
static int   transform(float *xptr, float *yptr, float *zptr, 
                       float nx, float ny, float nz, float d) ;

static int   mrisComputeTangentPlanes(MRI_SURFACE *mris) ;
static int   mrisRemoveLink(MRI_SURFACE *mris, int vno1, int vno2) ;
static int   mrisRemoveEdge(MRI_SURFACE *mris, int vno1, int vno2) ;
static int   mrisRemoveFace(MRI_SURFACE *mris, int fno) ;
static int   mrisCountTotalNeighbors(MRI_SURFACE *mris) ;
static int   mrisCountValidLinks(MRI_SURFACE *mris, int vno1, int vno2) ;
static int   mrisComputeSpringTerm(MRI_SURFACE *mris,INTEGRATION_PARMS *parms);
static int   mrisComputeDistanceTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeCorrelationTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputePolarCorrelationTerm(MRI_SURFACE *mris, 
                                              INTEGRATION_PARMS *parms) ;
static int   mrisComputeAngleAreaTerms(MRI_SURFACE *mris, 
                                       INTEGRATION_PARMS *parms) ;
static int   mrisClearGradient(MRI_SURFACE *mris) ;
static int   mrisClearMomentum(MRI_SURFACE *mris) ;
static int   mrisApplyGradient(MRI_SURFACE *mris, double dt) ;
static int   mrisValidVertices(MRI_SURFACE *mris) ;
static int   mrisLabelVertices(MRI_SURFACE *mris, float cx, float cy, 
                               float cz, int label, float radius) ;



static int mrisProjectSurface(MRI_SURFACE *mris) ;
static int mrisOrientSurface(MRI_SURFACE *mris) ;
static int   mrisComputeBoundaryNormals(MRI_SURFACE *mris) ;
static int   mrisSmoothBoundaryNormals(MRI_SURFACE *mris, int niter) ;
static int   mrisFlipPatch(MRI_SURFACE *mris) ;

/* not currently used */
#if 0
static int    mrisComputeBoundaryTerm(MRI_SURFACE *mris, 
                                      INTEGRATION_PARMS *parms) ;
static int   mrisComputeCurvatureTerm(MRI_SURFACE *mris, 
                                      INTEGRATION_PARMS *parms) ;
static int   mrisComputeNegTerm(MRI_SURFACE *mris,INTEGRATION_PARMS *parms);
static int   mrisCountNegativeVertices(MRI_SURFACE *mris) ;
static int   mrisComputeSpringNormalTerm(MRI_SURFACE *mris,
                                         INTEGRATION_PARMS *parms);
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

/*--------------------------------------------------------------------*/

/*--------------------- CONSTANTS AND MACROS -------------------------*/

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
  MRI_SURFACE *mris ;
  int         nfaces, nvertices, magic, version, ix, iy, iz, vno, fno, n, m,
              imnr, imnr0, imnr1, type ;
  float       x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  FILE        *fp ;
  VERTEX      *vertex ;
  FACE        *face ;

  type = mrisFileNameType(fname) ;
  if (type == MRIS_ASCII_FILE)
  {
    mris = mrisReadAsciiFile(fname) ;
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
    fread3(&nfaces, fp);
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "reading %d vertices and %d faces.\n",nvertices,nfaces);
    
    mris = MRISalloc(nvertices, nfaces) ;
    
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
      vertex->label = NO_LABEL ;
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
    
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      for (n = 0 ; n < 4 ; n++)
      {
        fread3(&mris->faces[fno].v[n],fp);
        if (version < 0)   /* new surface format */
          mris->vertices[mris->faces[fno].v[n]].num++;
      }
    }
    fclose(fp);
  }
  strcpy(mris->fname, fname) ;
  if (strstr(fname, "rh"))
    mris->hemisphere = RIGHT_HEMISPHERE ;
  else
    mris->hemisphere = LEFT_HEMISPHERE ;
  if ((version<0) || type == MRIS_ASCII_FILE)
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
  if (type == MRIS_ASCII_FILE)
  {
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
  }
  else 
  {
    if (MRISreadBinaryCurvature(mris, fname) != NO_ERROR)
      return(NULL) ;
    if (mrisReadBinaryAreas(mris, fname) != NO_ERROR)
      return(NULL) ;
  }

  mris->radius = MRISaverageRadius(mris) ;
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
  int         nfaces, nvertices, magic, version, ix, iy, iz, vno, fno, n, m,
              imnr, imnr0, imnr1, type ;
  float       x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  FILE        *fp ;
  VERTEX      *vertex ;
  FACE        *face ;

  type = mrisFileNameType(fname) ;
  if (type == MRIS_ASCII_FILE)
  {
    mris = mrisReadAsciiFile(fname) ;
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
    fread3(&nfaces, fp);
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "reading %d vertices and %d faces.\n",nvertices,nfaces);
    
    mris = MRISalloc(nvertices, nfaces) ;
    
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
      vertex->label = NO_LABEL ;
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
    
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      for (n = 0 ; n < 4 ; n++)
      {
        fread3(&mris->faces[fno].v[n],fp);
        if (version < 0)   /* new surface format */
          mris->vertices[mris->faces[fno].v[n]].num++;
      }
    }
    fclose(fp);
  }
  strcpy(mris->fname, fname) ;
  if (strstr(fname, "rh"))
    mris->hemisphere = RIGHT_HEMISPHERE ;
  else
    mris->hemisphere = LEFT_HEMISPHERE ;
  if ((version<0) || type == MRIS_ASCII_FILE)
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
#if 0
  mrisComputeVertexDistances(mris) ;
  mrisReadTransform(mris, fname) ;
#endif
  if (type == MRIS_ASCII_FILE)
  {
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
  }
  else 
  {
    if (MRISreadBinaryCurvature(mris, fname) != NO_ERROR)
      fprintf(stderr, "ignoring curvature file...\n") ; /*return(NULL) ;*/
#if 0
    if (mrisReadBinaryAreas(mris, fname) != NO_ERROR)
      return(NULL) ;
#endif
  }


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
  int   k,n, type ;
  float x,y,z;
  FILE  *fp;

  type = mrisFileNameType(fname) ;
  if (type == MRIS_ASCII_FILE)
    return(MRISwriteAscii(mris, fname)) ;

  fp = fopen(fname,"w");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE,"MRISwrite(%s): can't create file\n",fname));
  fwrite3(-1,fp);
  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);
  for (k = 0 ; k < mris->nvertices ; k++)
  {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
    fwrite2((int)(x*100),fp);
    fwrite2((int)(y*100),fp);
    fwrite2((int)(z*100),fp);
  }
  for (k = 0 ; k < mris->nfaces ; k++)
  {
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      fwrite3(mris->faces[k].v[n],fp);
  }
  fclose(fp);
  return(NO_ERROR) ;
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
      n = v->n[m];
      f = &mris->faces[v->f[m]];
      n0 = (n==0)?3:n-1;
      n1 = (n==3)?0:n+1;
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
      if (i==v->num)
        printf("face[%d].v[%d] = %d\n",k,m,f->v[m]);
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
  for (max_possible = 0, n = 1 ; n <= max_nbhd ; n++)
    max_possible += nbrs[n] ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
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
          } while (!done) ;
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
      char  fname[100] ;

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
  return(NO_ERROR) ;
}
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
static int
mrisRemoveRipped(MRI_SURFACE *mris)
{
  int     vno, n, fno, nripped, tno ;
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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
      mris->orig_area += face->orig_area[tno] ;
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
  int       k,n;
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
    for (n=0;n<v->num;n++) if (!mris->faces[v->f[n]].ripflag)
    {
      mrisNormalFace(mris, v->f[n],v->n[n],norm);
      snorm[0] += norm[0];
      snorm[1] += norm[1];
      snorm[2] += norm[2];

      /* Note: overestimates area by *2 !! */
      v->area += mrisTriangleArea(mris, v->f[n],v->n[n]); 
    }
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

  n0 = (n==0)?3:n-1;
  n1 = (n==3)?0:n+1;
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

  n0 = (n==0)?3:n-1;
  n1 = (n==3)?0:n+1;
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
  char transform_fname[100], fpref[300] ;

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
  char   fname[100], fpref[100], hemi[20] ;

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
  char   *cp, path[100], fname[100] ;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/%s", path, sname) ;
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
static int
mrisReadBinaryAreas(MRI_SURFACE *mris, char *mris_fname)
{
  int   k,vnum,fnum;
  float f;
  FILE  *fp;
  char  fname[100], fpref[100], hemi[20] ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "reading area file...") ;

  FileNamePath(mris_fname, fpref) ;
  strcpy(hemi, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh") ;
  sprintf(fname, "%s/%s.area", fpref, hemi) ;

  /*  mris->orig_area = 0.0f ;*/
  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADPARM, 
              (ERROR_BADPARM,"mrisReadBinaryAreas: no area file %s\n",fname));
  fread3(&vnum,fp);
  fread3(&fnum,fp);
  if (vnum!=mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "mrisReadBinaryAreas: incompatible vertex "
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
  char  *cp, fname[100], path[100] ;

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

#if 0
  {
    VERTEX *va, *vb, *vc ;
    FACE   *face ;
    
    vdst = &mris_dst->vertices[0] ;
    face = &mris_dst->faces[vdst->f[0]] ;
    va = &mris_dst->vertices[face->v[3]] ;
    vb = &mris_dst->vertices[face->v[1]] ;
    vc = &mris_dst->vertices[face->v[2]] ;
    vdst->x = 0 ;
    vdst->y = 0 ;
    vdst->z = 0 ;
    va->x = 1 ;
    va->y = 0 ;
    va->z = 0 ;
    vb->x = 0.5 ;
    vb->y = 1 ;
    vb->z = 0 ;
    vc->x = 1 ;
    vc->y = 1 ;
    vc->z = 0 ;
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
    vdst->ox = vsrc->ox ;
    vdst->oy = vsrc->oy ;
    vdst->oz = vsrc->oz ;
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

  if (Gdiag & DIAG_SHOW)
  {
    char fname[100] ;

    sprintf(fname, "%s.%s.out", 
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    if (Gdiag & DIAG_WRITE)
    {
      parms->fp = fopen(fname, "w") ;
      mrisLogIntegrationParms(parms->fp, mris, parms) ;
    }
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
    mrisComputeSpringTerm(mris, parms) ;
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
    char fname[100] ;

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
  char    fname[100], base_name[100], path[100] ;
  double  base_dt ;
  struct  timeb start ;
  
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

  for (sno = 1 ; sno < SURFACES-1 ; sno++)
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
      if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  "MRISregister", fname) ;
      
      MRISsetNeighborhoodSize(mris, -1) ;  /* back to max */
      MRIScomputeMetricProperties(mris) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
      MRISnormalizeCurvature(mris) ;
      MRISresetNeighborhoodSize(mris,1);/*only use nearest neighbor distances*/
    }
    MRISstoreMeanCurvature(mris) ;

    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "calculating curvature of %s surface\n",fname) ;
    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp,"calculating curvature of %s surface\n",fname);
#if 1
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "finding optimal rigid alignment\n") ;
    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp, "finding optimal rigid alignment\n") ;
    parms->mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
    parms->mrisp_template = mrisp_template ;
    MRISrigidBodyAlignGlobal(mris, parms, 4.0f, 16.0f, 8) ;
    MRISPfree(&parms->mrisp) ;
#endif

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
#if 0
      if (i == 10)
      {
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, "finding optimal rigid alignment\n") ;
        if (Gdiag & DIAG_WRITE)
          fprintf(parms->fp, "finding optimal rigid alignment\n") ;
        MRISrigidBodyAlignGlobal(mris, parms, 4.0f, 16.0f, 8) ;
      }
      MRISrigidBodyAlignLocal(mris, parms) ;
#endif
      mrisClearMomentum(mris) ;
      done = 0 ;
#if 0
      steps = mrisIntegrate(mris, parms, parms->n_averages) ;
      parms->start_t += steps ;
#else
       mrisIntegrationEpoch(mris, parms, parms->n_averages) ;
#endif
    }
  }

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
  
  if (Gdiag & DIAG_SHOW)
  {
    char fname[100] ;

    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      if (nbrs[i])
        fprintf(stderr, "%d: %d | ", i, nbrs[i]) ;
    fprintf(stderr, "\n") ;
    if (Gdiag & DIAG_WRITE)
    {
      sprintf(fname, "%s.%s.out", 
              mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
              parms->base_name);
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      mrisLogIntegrationParms(parms->fp, mris,parms) ;
      for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
        if (nbrs[i])
          fprintf(parms->fp, "%d: %d | ", i, nbrs[i]) ;
      fprintf(parms->fp, "\n") ;
    }
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
        fprintf(stderr, "resampling long-range distances...") ;
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
      MRISsampleDistances(mris, nbrs, parms->nbhd_size) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      mrisClearMomentum(mris) ;
    }

    if (!passno)
      mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 2);

    
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
#if 0
      fprintf(stderr, "epoch %d of %d starting distance error %%%2.2f\n",
              (int)(passno*NCOEFS+i+1), (int)(max_passes*NCOEFS), 
              (float)pct_error);
#endif
    
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
      fprintf(stderr, "pass %d: start=%2.1f, end=%2.1f, ratio=%2.3f\n",
              passno+1, starting_sse, ending_sse, 
              (starting_sse-ending_sse)/starting_sse) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp, "pass %d: start=%2.4f, end=%2.4f, ratio=%2.4f\n",
                passno+1, starting_sse, ending_sse, 
                (starting_sse-ending_sse)/starting_sse) ;
    }
  } while (
           !FZERO(ending_sse) && 
           (((starting_sse-ending_sse)/starting_sse) > parms->tol) &&
           (++passno < max_passes)
           ) ;

#if 0
  fprintf(stderr, "initial optimization complete - settling to equilibrium...\n") ;
  parms->niterations = 1000 ;
  mrisIntegrationEpoch(mris, parms, parms->n_averages = 0) ; 
#endif

  /* finally, remove all the small holes */
  parms->l_area = 1.0f ;
  parms->l_dist = 0.1f ;  /* was 0.001 */
  parms->l_angle = ANGLE_AREA_SCALE * parms->l_area ;
  parms->niterations = niter ;
  fprintf(stderr, "removing remaining folds...\n") ;
  mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 3);

  if (mris->patch)  /* smooth out remaining folds */
  {
    parms->l_spring = 1.0f ;
    parms->niterations = 5 ;
    parms->integration_type = INTEGRATE_MOMENTUM ;
    parms->dt = 0.5f ; parms->momentum = 0.0f ;
    parms->n_averages = 0 ;
    mrisRemoveNegativeArea(mris, parms, 0, MAX_NEG_AREA_PCT, 1);
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
  
  if (Gdiag & DIAG_SHOW)
  {
    char fname[100] ;

    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      if (nbrs[i])
        fprintf(stderr, "%d: %d | ", i, nbrs[i]) ;
    fprintf(stderr, "\n") ;
    if (Gdiag & DIAG_WRITE)
    {
      sprintf(fname, "%s.out", parms->base_name) ;
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      mrisLogIntegrationParms(parms->fp, mris,parms) ;
      for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
        if (nbrs[i])
          fprintf(parms->fp, "%d: %d | ", i, nbrs[i]) ;
      fprintf(parms->fp, "\n") ;
    }
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
          fprintf(stderr, "resampling long-range distances...") ;
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
      fprintf(stderr, "pass %d: start=%2.1f, end=%2.1f, ratio=%2.3f\n",
              passno, starting_sse, ending_sse, 
              (starting_sse-ending_sse)/ending_sse) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp, "pass %d: start=%2.4f, end=%2.4f, ratio=%2.4f\n",
                passno, starting_sse, ending_sse, 
                (starting_sse-ending_sse)/starting_sse) ;
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
  
  pct_neg = 100.0*mris->neg_area/(mris->neg_area+mris->total_area) ;
  if (pct_neg <= min_area_pct)
    return(0) ;   /* no steps */

  tol = parms->tol ;
  parms->tol = 1e-2 ;
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
  MRIScomputeTriangleProperties(mris,0) ;  /* compute areas and normals */
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

  ratio = *pnum / *pdenom ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
  if (Gdiag & DIAG_WRITE)
    fprintf(parms->fp, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;

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
  MRIScomputeTriangleProperties(mris,0) ;  /* compute areas and normals */
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
          sse_thresh, pct_neg, pct_neg_area, total_vertices 
          /*, scale, last_neg_area */ ;

  l_spring = parms->l_spring ;
  l_dist = parms->l_dist ;
  l_area = parms->l_area ;
  write_iterations = parms->write_iterations ;
  niterations = parms->niterations ;
  sse_thresh = parms->tol ;
    
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
            "vertex %d curvature = %2.5f, position = (%2.3f,%2.3f,%2.3f)\n",
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
    /*    mrisComputeSpringTerm(mris, parms) ;*/
#if 0
    mrisComputeCurvatureTerm(mris, parms) ;
    mrisComputeNegTerm(mris, parms) ;
    mrisComputeBoundaryTerm(mris, parms) ;
#endif

    mrisAverageGradients(mris, n_averages) ;
    mrisComputeSpringTerm(mris, parms) ;
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
      fprintf(stderr, "vertex %d curvature = %2.3f\n",
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
      if ((100*(old_sse - sse) / sse) < parms->tol)
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
  int     ano, fno, tno, ntriangles, total_neighbors ;
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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
#if 0
if (face->area[tno] >= 0.0f)
  continue ;
#endif
      ntriangles++ ;
      delta = (double)(area_scale * face->area[tno] - face->orig_area[tno]) ;
      sse_area += delta*delta ;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      {
        delta = deltaAngle(face->angle[tno][ano], face->orig_angle[tno][ano]);
        sse_angle += delta*delta ;
      }
      if (!finite(sse_area) || !finite(sse_angle))
        ErrorExit(ERROR_BADPARM, "sse is not finite at face %d:%d!\n",fno,tno);
    }
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
          area_scale, sse_corr, sse_neg_area, l_corr ;
  int     ano, tno, fno ;
  FACE    *face ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  sse_corr = sse_angle = sse_neg_area = 
    sse_area = sse_spring = sse_curv = sse_dist = 0.0 ;
  if (!FZERO(parms->l_angle)||!FZERO(parms->l_area)||(!FZERO(parms->l_parea)))
  {
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      if (face->ripflag)
        continue ;
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
      {
        delta = (double)(area_scale*face->area[tno] - face->orig_area[tno]) ;
#if ONLY_NEG_AREA_TERM
        if (face->area[tno] < 0.0f)
          sse_neg_area += delta*delta ;
#endif
        sse_area += delta*delta ;
        for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        {
          delta = deltaAngle(face->angle[tno][ano],face->orig_angle[tno][ano]);
#if ONLY_NEG_AREA_TERM
        if (face->angle[tno][ano] >= 0.0f)
          delta = 0.0f ;
#endif
          sse_angle += delta*delta ;
        }
        if (!finite(sse_area) || !finite(sse_angle))
          ErrorExit(ERROR_BADPARM, "sse not finite at face %d:%d!\n",fno,tno);
      }
    }
  }
  if (!FZERO(parms->l_dist))
    sse_dist = mrisComputeDistanceError(mris) ;
  if (!FZERO(parms->l_spring))
    sse_spring = mrisComputeSpringEnergy(mris) ;
  if (!FZERO(parms->l_curv))
    sse_curv = MRIScomputeFolding(mris) ;
  l_corr = (double)(parms->l_corr + parms->l_pcorr) ;
  if (!FZERO(l_corr))
    sse_corr = mrisComputeCorrelationError(mris, parms, 1) ;

  sse = 
    (double)parms->l_area   * sse_neg_area + 
    (double)parms->l_parea  * sse_area + 
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
  int     fno, tno, ano ;
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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      v = tno == 0 ? &mris->vertices[face->v[0]]:&mris->vertices[face->v[2]] ;
      dot = v->x * face->nx[tno] + v->y * face->ny[tno] + v->z * face->nz[tno];
      if (dot < 0.0f)   /* not in same direction, area < 0 and reverse n */
      {
        face->area[tno] *= -1.0f ;
        face->nx[tno] *= -1.0f; face->ny[tno] *= -1.0f; face->nz[tno] *= -1.0f;
        for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
          face->angle[tno][ano] *= -1.0f ;
      }
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
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
      {
        if (face->area[tno] >= 0.0f)
          mris->total_area += face->area[tno] ;
        else
          mris->neg_area += -face->area[tno] ;
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
  int     fno, tno, ano, vno ;
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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      if (face->nz[tno] < 0.0f)   
      {
        /* not in same direction, area < 0 and reverse n */
        face->area[tno] *= -1.0f ;
        face->nx[tno] *= -1.0f; face->ny[tno] *= -1.0f; face->nz[tno] *= -1.0f;
        for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
          face->angle[tno][ano] *= -1.0f ;
      }
    }
  }

  /* now recompute the total surface area, ignoring negative areas */
  mris->total_area = mris->neg_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      if (face->area[tno] >= 0.0f)
        mris->total_area += face->area[tno] ;
      else
        mris->neg_area += -face->area[tno] ;
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
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
        v->area += face->area[tno] ;
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
  int     tno, fno, negative ;
  FACE    *face ;

  negative = 0 ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
      if (face->area[tno] < 0.0f)
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
  MRIScomputeTriangleProperties(mris,0) ;  /* recompute areas and normals */
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
  int     ano, tno, vnum,fnum, fno, vno ;
  FACE    *face ;
  VERTEX  *v ;
  float   f;
  FILE    *fp;
  char    fname[100], fpref[100], hemi[20], *cp ;

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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      f = freadFloat(fp);
      face->orig_area[tno] = f ;
      mris->orig_area += f;
    }
  }

  /* compute original vertex areas from faces */
  for (vno=0;vno<vnum;vno++)
  {
    v = &mris->vertices[vno] ;
    v->origarea = 0.0f ;
    for (fno = 0 ; fno < v->num ; fno++)
    {
      face = &mris->faces[v->f[fno]] ;
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
        v->origarea += face->orig_area[tno] ;
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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      {
        f = freadFloat(fp);
        face->orig_angle[tno][ano] = f ;
      }
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
  int     fno, tno, ano, vno ;
  FACE    *face ;
  FILE    *fp;
  char    fname[100], fpref[100], hemi[20], *cp ;

  MRIScomputeTriangleProperties(mris,0) ;

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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
      putf(face->area[tno], fp) ;
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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        putf(face->angle[tno][ano], fp) ;
    }
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
      putf(v->dist[n], fp) ;
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

------------------------------------------------------*/
int
MRIScomputeTriangleProperties(MRI_SURFACE *mris, int no_angles)
{
  VECTOR  *v_a, *v_b, *v_c, *v_d, *v_n ;
  VERTEX  *v0, *v1, *v2, *v3, *va, *vb, *vo, *v ;
  FACE    *face ;
  int     tno, fno, ano, vno  ;
  float   area, angle, dot, cross, dz ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_c = VectorAlloc(3, MATRIX_REAL) ;
  v_d = VectorAlloc(3, MATRIX_REAL) ;
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
    v3 = &mris->vertices[face->v[3]] ;
    VERTEX_EDGE(v_a, v0, v1) ;  
    VERTEX_EDGE(v_d, v0, v3) ;
    VERTEX_EDGE(v_b, v2, v1) ;
    VERTEX_EDGE(v_c, v2, v3) ;

    /* compute metric properties of first triangle */
    V3_CROSS_PRODUCT(v_a, v_d, v_n) ;
    area = V3_LEN(v_n) * 0.5f ;
    dot = V3_DOT(v_a, v_d) ;
    face->area[0] = area ;
    V3_NORMALIZE(v_n, v_n) ;             /* make it a unit vector */
    face->nx[0] = V3_X(v_n); face->ny[0] = V3_Y(v_n); face->nz[0] = V3_Z(v_n);
    mris->total_area += area ;

    /* compute metric properties of second triangle */
    V3_CROSS_PRODUCT(v_c, v_b, v_n) ;
    face->area[1] = area ;
    V3_NORMALIZE(v_n, v_n) ;             /* make it a unit vector */
    face->nx[1] = V3_X(v_n); face->ny[1] = V3_Y(v_n); face->nz[1] = V3_Z(v_n);
    mris->total_area += area ;

    /* now compute angles */
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      if (no_angles)
        break ;
      VECTOR_LOAD(v_n, face->nx[tno], face->ny[tno], face->nz[tno]) ;
      if ((V3_X(v_n) < V3_Y(v_n)) && (V3_X(v_n) < V3_Z(v_n)))
        dz = fabs(V3_X(v_n)) ;
      else if (V3_Y(v_n) < V3_Z(v_n))
        dz = fabs(V3_Y(v_n)) ;
      else
        dz = fabs(V3_Z(v_n)) ;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      {
        if (tno == 0) switch (ano)   /* vertices for triangle 1 */
        {
        default:
        case 0: vo = v0 ; va = v3 ; vb = v1 ; break ;
        case 1: vo = v1 ; va = v0 ; vb = v3 ; break ;
        case 2: vo = v3 ; va = v1 ; vb = v0 ; break ;
        }
        else switch (ano)             /* vertices for triangle 2 */
        {
        default:
        case 0: vo = v1 ; va = v3 ; vb = v2 ; break ;
        case 1: vo = v2 ; va = v1 ; vb = v3 ; break ;
        case 2: vo = v3 ; va = v2 ; vb = v1 ; break ;
        }
        VERTEX_EDGE(v_a, vo, va) ;VERTEX_EDGE(v_b, vo, vb) ;
        
        cross = VectorTripleProduct(v_b, v_a, v_n) ;
        dot = V3_DOT(v_a, v_b) ;
        angle = atan2(cross, dot) ;
        face->angle[tno][ano] = angle ;

#if 0
        if (angle < 0.0f || angle >= M_PI)
          fprintf(stderr, "angle [%d][%d][%d] = %2.1f\n",
                  fno,tno,ano,(float)DEGREES(angle)) ;
#endif
      }
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
    {
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
        v->area += mris->faces[v->f[fno]].area[tno] ;
    }
    v->area /= 2.0 ;
  }

  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_c) ;
  VectorFree(&v_d) ;
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
      fprintf(stderr, "moving vertex %d by (%2.3f, %2.3f, %2.3f) --> "
              "(%2.1f, %2.1f, %2.1f), nz=%2.2f, a=%2.3f\n",
              vno, vertex->odx, vertex->ody, vertex->odz,
              vertex->x, vertex->y, vertex->z, vertex->nz,vertex->area) ;
    vertex->x += vertex->odx ; 
    vertex->y += vertex->ody ;
    vertex->z += vertex->odz ;
  }


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
  char    fname[100] ;
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
  char    fname[100] ;
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
  char   fname[100], *cp, path[100] ;
  FILE   *fp;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/%s", path, sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */
  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "MRISwriteCurvature: could not open %s", 
                 fname)) ;

  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    curv = mris->vertices[k].curv ;
    i = nint(curv * 1000.0) ;
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
MRISwriteAreaError(MRI_SURFACE *mris, char *fname)
{
  int    vno, fno, tno, i, n ;
  float  area, orig_area ;
  FACE   *face ;
  VERTEX *vertex ;
  FILE   *fp;
  
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
    for (n = fno = 0 ; fno < vertex->num ; fno++)
    {
      face = &mris->faces[vertex->f[fno]] ;
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++, n++)
      {
        area += face->area[tno] ;
        orig_area += face->orig_area[tno] ;
      }
    }
    i = nint((area-orig_area) * 100.0f / (float)n) ;
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
MRISwriteAngleError(MRI_SURFACE *mris, char *fname)
{
  int    vno, fno, tno, ano, i ;
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
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
      {
        for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        {
          error += 
            fabs(deltaAngle(face->angle[tno][ano],face->orig_angle[tno][ano]));
        }
      }
      error /= (float)(v->num*TRIANGLES_PER_FACE*ANGLES_PER_TRIANGLE) ;
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
MRISwriteValues(MRI_SURFACE *mris, char *fname)
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
    fprintf(stderr, "avg = %f, stdev = %f, min = %f, max = %f\n",
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
    fread(&f,1,sizeof(float),fp);
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
MRISreadCanonicalCoordinates(MRI_SURFACE *mris, char *sname)
{
  int         nvertices, magic, version, ix, iy, iz, vno, n, nfaces, crap ;
  FILE        *fp ;
  VERTEX      *vertex ;
  char        fname[100], path[100], *cp ;
  float       d, x, y, z, r, theta, phi ;

  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
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
  fread3(&nfaces, fp);

  if ((nvertices != mris->nvertices) || (nfaces != mris->nfaces))
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM,
                "MRISreadCanonicalSurface(%s): nfaces %d (%d) or nvertices "
                "%d (%d) mismatch", nvertices, mris->nvertices,
                nfaces, mris->nfaces)) ;

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
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  r = mris->radius = MRISaverageRadius(mris) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
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
  char        fname[100], path[100], *cp ;

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


#if 0
  fread(&npts,1,sizeof(int),fp);
#else
  npts = freadInt(fp) ;
#endif
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading patch %s with %d vertices (%2.1f%% of total)\n",
            pname, npts, 100.0f*(float)npts/(float)mris->nvertices) ;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].ripflag = TRUE;
  for (j=0;j<npts;j++)
  {
#if 0
    fread(&i,1,sizeof(int),fp);
#else
    i = freadInt(fp) ;
#endif
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
  mrisRipFaces(mris);
  mris->patch = 1 ;
  mris->status = MRIS_CUT ;

  mrisRemoveRipped(mris) ;
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
static int
mrisRipFaces(MRI_SURFACE *mris)
{
  int n,k;
  face_type *f;

  for (k=0;k<mris->nfaces;k++)
    mris->faces[k].ripflag = FALSE;
  for (k=0;k<mris->nfaces;k++)
  {
    f = &mris->faces[k];
    for (n=0;n<4;n++)
      if (mris->vertices[f->v[n]].ripflag)
        f->ripflag = TRUE;
  }
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].border = FALSE;
  for (k=0;k<mris->nfaces;k++)
  if (mris->faces[k].ripflag)
  {
    f = &mris->faces[k];
    for (n=0;n<4;n++)
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
  if (type == MRIS_ASCII_FILE)
    return(MRISwritePatchAscii(mris, fname)) ;

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

#if 0
  fwrite(&npts,1,sizeof(int),fp);
#else
  fwriteInt(npts, fp) ;
#endif
  for (k=0;k<mris->nvertices;k++)
  if (!mris->vertices[k].ripflag)
  {
    i = (mris->vertices[k].border)?-(k+1):k+1;
#if 0
    fwrite(&i,1,sizeof(int),fp);
#else
    fwriteInt(i, fp) ;
#endif
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

    if (vno == 4013)
      DiagBreak() ;
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

#if 0
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
      char fname[100] ;

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
      fprintf(stderr, "moving vertex %d by (%2.3f, %2.3f, %2.3f)\n",
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

  /*  mrisRipFaces(mris) ;*/
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

#if 0
  mris->vertices[414].ripflag = 1 ;
  mris->vertices[3293].ripflag = 1 ;
  mris->vertices[3298].ripflag = 1 ;
  mris->vertices[5025].ripflag = 1 ;
  mris->vertices[8463].ripflag = 1 ;
  mris->vertices[17830].ripflag = 1 ;
  mris->vertices[21865].ripflag = 1 ;
  mris->vertices[25909].ripflag = 1 ;
  mris->vertices[28263].ripflag = 1 ;
  mris->vertices[34104].ripflag = 1 ;
  mris->vertices[102960].ripflag = 1 ;
  mris->vertices[4013].ripflag = 1 ;
  mrisRipFaces(mris) ;
#endif
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

  nedges += nfaces ;        /* one additional edge added for each triangle */
  nfaces *= 2 ;             /* two triangular faces per face */
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

    if (vno == 4013)
      DiagBreak() ;
    if (vno == 3427)
      DiagBreak() ;
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
  

  /*  mrisRipFaces(mris) ;*/
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
int
MRIScomputeSecondFundamentalForm(MRI_SURFACE *mris)
{
  int    vno, i, n, vmax ;
  VERTEX *vertex, *vnb ;
  MATRIX *m_U, *m_Ut, *m_tmp1, *m_tmp2, *m_inverse, *m_eigen, *m_Q ;
  VECTOR *v_c, *v_z, *v_n, *v_e1, *v_e2, *v_yi ;
  float  k1, k2, evalues[3], a11, a12, a21, a22 ;
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
      n++ ;
    }

    m_Ut = MatrixTranspose(m_U, NULL) ;          /* Ut */
    m_tmp2 = MatrixMultiply(m_Ut, m_U, NULL) ;   /* Ut U */
    m_inverse = MatrixInverse(m_tmp2, NULL) ;    /* (Ut U)^-1 */
    if (!m_inverse)   /* singular matrix - must be planar?? */
    {
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

      /* the columns of m_eigen will be the eigenvectors of m_Q */
      if (MatrixEigenSystem(m_Q, evalues, m_eigen) == NULL)
      {
        MatrixFree(&m_Ut) ;
        MatrixFree(&m_tmp2) ;
        MatrixFree(&m_U) ;
        VectorFree(&v_z) ;
        MatrixFree(&m_tmp1) ;
        MatrixFree(&m_inverse) ;
        vertex->K = vertex->H = 0.0 ;
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
#if 0
    vertex->e1x = V3_X(v_e1) * a11 + V3_X(v_e2) * a12 ;
    vertex->e1y = V3_Y(v_e1) * a11 + V3_Y(v_e2) * a12 ;
    vertex->e1z = V3_Z(v_e1) * a11 + V3_Z(v_e2) * a12 ;
    vertex->e2x = V3_X(v_e1) * a21 + V3_X(v_e2) * a22 ;
    vertex->e2y = V3_Y(v_e1) * a21 + V3_Y(v_e2) * a22 ;
    vertex->e2z = V3_Z(v_e1) * a21 + V3_Z(v_e2) * a22 ;
#endif

    MatrixFree(&m_Ut) ;
    MatrixFree(&m_tmp2) ;
    MatrixFree(&m_U) ;
    VectorFree(&v_z) ;

  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "max H error=%2.3f at %d\n", max_error, vmax) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "total area = %2.3f\n", total_area);

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
  int    vno, fno, tno, n ;
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
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++, n++)
      {
        area += face->area[tno] ;
        orig_area += face->orig_area[tno] ;
      }
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
    if (vertex->ripflag)
      continue ;
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
    if (vertex->ripflag)
      continue ;
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
  double  delta_t = 0.0, rms_height, desired_rms_height ;

  write_iterations = parms->write_iterations ;
  n_averages = parms->n_averages ;

  if (Gdiag & DIAG_WRITE)
  {
    char fname[100] ;

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
  
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 
            0, 0.0f, (float)rms_height, n_averages) ;
  else
    fprintf(stderr, "\rstep %3.3d: RMS=%2.3f (target=%2.3f)   ", 0, 
            rms_height, desired_rms_height);
  if (Gdiag & DIAG_WRITE)
  {
    fprintf(parms->fp, "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 
            0, 0.0f, (float)rms_height, n_averages) ;
    fflush(parms->fp) ;
  }

  MRISclearCurvature(mris) ;   /* curvature will be used to calculate sulc */
  for (n = 0 ; n < niterations ; n++)
  {
    mrisClearGradient(mris) ;
    mrisComputeDistanceTerm(mris, parms) ;
    mrisComputeSpringTerm(mris, parms) ;
    
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
    /*
      only compute the second fundamental form (called from 
      MRISupdateSurface) every 5th iteration as it is the most
      expensive part of the inflation.
      */
    MRIScomputeMetricProperties(mris) ; 
    if (!((n+1) % 5))   /* compute curvature also */
    {
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
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISstoreCurrentPositions(MRI_SURFACE *mris)
{
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
  int     vno, nvertices, fno, ano, tno, n ;
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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      f->orig_area[tno] = f->area[tno] ;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        f->orig_angle[tno][ano] = f->angle[tno][ano] ;
    }
  }
  mris->orig_area = mris->total_area ;
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
    if (d <= radius)
      v->label = label ;
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
        File Format is:

        name x y z
------------------------------------------------------*/
static int
mrisComputeSpringNormalTerm(MRI_SURFACE *mris,INTEGRATION_PARMS *parms)
{
  int     vno, n, m ;
  VERTEX  *vertex, *vn ;
  double  l_spring ;
  float   sx, sy, sz, nx, ny, nz, nc, x, y, z ;

  l_spring = parms->l_spring ;

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
      fprintf(stderr, "vertex %d spring term: (%2.3f, %2.3f, %2.3f)\n",
              vno, vertex->dx, vertex->dy, vertex->dz) ;
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
mrisComputeSpringTerm(MRI_SURFACE *mris,INTEGRATION_PARMS *parms)
{
  int     vno, n, m ;
  VERTEX  *vertex, *vn ;
  double  l_spring ;
  float   sx, sy, sz, x, y, z, dist_scale ;

  l_spring = parms->l_spring ;

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
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    if (vertex->border && !vertex->neg)
      continue ;

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
#if 1
    n = 4 ;  /* avg # of nearest neighbors */
#endif
    if (n>0)
    {
      sx = dist_scale*sx/n;
      sy = dist_scale*sy/n;
      sz = dist_scale*sz/n;
    }
    
    vertex->dx += l_spring * sx ;
    vertex->dy += l_spring * sy ;
    vertex->dz += l_spring * sz ;
    if (vno == Gdiag_no)
      fprintf(stderr, "vertex %d spring term: (%2.3f, %2.3f, %2.3f)\n",
              vno, vertex->dx, vertex->dy, vertex->dz) ;
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
            "computing distance term for vertex %d @ (%2.2f, %2.2f, %2.2f)\n",
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
      fprintf(stderr, "vertex %d, distance term: (%2.3f, %2.3f, %2.3f)\n",
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
  char  *cp, host_name[100] ;

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
  if (!FZERO(parms->l_angle))
    fprintf(fp, ", l_angle=%2.4f", parms->l_angle) ;
  if (!FZERO(parms->l_corr))
    fprintf(fp, ", l_corr=%2.4f", parms->l_corr) ;
  if (!FZERO(parms->l_spring))
    fprintf(fp, ", l_spring=%2.4f", parms->l_spring) ;
  if (!FZERO(parms->l_dist))
    fprintf(fp, ", l_dist=%2.4f", parms->l_dist) ;
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
  char fname[100], path[100], base_name[100] ;

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
------------------------------------------------------*/
static int
mrisComputeAngleAreaTerms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     fno, ano, tno ;
  VERTEX  *v0, *v1, *v2, *v3, *va, *vb, *vo ;
  VECTOR  *v_a, *v_b, *v_c, *v_d, *v_a_x_n, *v_b_x_n, *v_c_x_n, *v_d_x_n,
          *v_n0, *v_n, *v_n1, *v_tmp, *v_sum ;
  FACE    *face ;
  float   orig_area0, area0, orig_area1, area1, l_parea,
          delta0, delta1, l_area, l_angle, delta, len, area_scale ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0f ;
  else
    area_scale = mris->total_area / mris->orig_area ;
#else
  area_scale = 1.0f ;
#endif

  l_angle = parms->l_angle ;
  l_area = parms->l_area ;
  l_parea = parms->l_parea ;

  if (FZERO(l_area) && FZERO(l_angle) && FZERO(l_parea))
    return(NO_ERROR) ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_c = VectorAlloc(3, MATRIX_REAL) ;
  v_d = VectorAlloc(3, MATRIX_REAL) ;
  v_n0 = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_n1 = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */

  v_tmp = VectorAlloc(3, MATRIX_REAL) ;   
  v_sum = VectorAlloc(3, MATRIX_REAL) ;   
  v_a_x_n = VectorAlloc(3, MATRIX_REAL) ;      
  v_b_x_n = VectorAlloc(3, MATRIX_REAL) ;      
  v_c_x_n = VectorAlloc(3, MATRIX_REAL) ;      
  v_d_x_n = VectorAlloc(3, MATRIX_REAL) ;      

  /* calculcate movement of each vertex caused by each triangle */
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    VECTOR_LOAD(v_n0, face->nx[0], face->ny[0], face->nz[0]) ;
    VECTOR_LOAD(v_n1, face->nx[1], face->ny[1], face->nz[1]) ;
    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;
    v3 = &mris->vertices[face->v[3]] ;
    VERTEX_EDGE(v_a, v0, v1) ;  
    VERTEX_EDGE(v_d, v0, v3) ;
    VERTEX_EDGE(v_b, v2, v1) ;
    VERTEX_EDGE(v_c, v2, v3) ;
    orig_area0 = area_scale * face->orig_area[0] ; area0 = face->area[0] ;
    orig_area1 = area_scale * face->orig_area[1] ; area1 = face->area[1] ;
    delta0 = delta1 = 0.0 ;
    if (!FZERO(l_parea))
    {
      delta0 += l_parea * (area0 - orig_area0) ; 
      delta1 += l_parea * (area1 - orig_area1) ;
    }
    if (!FZERO(l_area))
    {
      if (area0 <= 0.0f)
        delta0 += l_area * (area0 - orig_area0) ; 
      if (area1 <= 0.0f)
        delta1 += l_area * (area1 - orig_area1) ;
    }

    V3_CROSS_PRODUCT(v_a, v_n0, v_a_x_n) ;
    V3_CROSS_PRODUCT(v_b, v_n1, v_b_x_n) ;
    V3_CROSS_PRODUCT(v_c, v_n1, v_c_x_n) ;
    V3_CROSS_PRODUCT(v_d, v_n0, v_d_x_n) ;

    if (fno == 1235 || fno == 1237)
      DiagBreak() ;
  
    /* calculate movement of vertices in order, 0-3 */
    
    /* v0 */
    V3_SCALAR_MUL(v_a_x_n, -1.0f, v_sum) ;
    V3_ADD(v_sum, v_d_x_n, v_sum) ;
    V3_SCALAR_MUL(v_sum, delta0, v_sum) ;
    v0->dx += V3_X(v_sum) ; v0->dy += V3_Y(v_sum) ; v0->dz += V3_Z(v_sum) ;
    
    /* v1 */
    V3_SCALAR_MUL(v_d_x_n, -delta0, v_sum) ;
    V3_SCALAR_MUL(v_c_x_n, delta1, v_tmp) ;
    V3_ADD(v_tmp, v_sum, v_sum) ;
    v1->dx += V3_X(v_sum) ; v1->dy += V3_Y(v_sum) ; v1->dz += V3_Z(v_sum) ;
    
    /* v2 */
    V3_SCALAR_MUL(v_c_x_n, -1.0f, v_sum) ;
    V3_ADD(v_sum, v_b_x_n, v_sum) ;
    V3_SCALAR_MUL(v_sum, delta1, v_sum) ;
    v2->dx += V3_X(v_sum) ; v2->dy += V3_Y(v_sum) ; v2->dz += V3_Z(v_sum) ;
    
    /* v3 */
    V3_SCALAR_MUL(v_b_x_n, -delta1, v_sum) ;
    V3_SCALAR_MUL(v_a_x_n, delta0, v_tmp) ;
    V3_ADD(v_tmp, v_sum, v_sum) ;
    v3->dx += V3_X(v_sum) ; v3->dy += V3_Y(v_sum) ; v3->dz += V3_Z(v_sum) ;
    
    /* now calculate the angle contributions */
    if (!FZERO(l_angle)) for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      v_n = tno == 0 ? v_n0 : v_n1 ;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      {
        if (tno == 0) switch (ano)   /* vertices for triangle 1 */
        {
        default:
        case 0: vo = v0 ; va = v3 ; vb = v1 ; break ;
        case 1: vo = v1 ; va = v0 ; vb = v3 ; break ;
        case 2: vo = v3 ; va = v1 ; vb = v0 ; break ;
        }
        else switch (ano)             /* vertices for triangle 2 */
        {
        default:
        case 0: vo = v1 ; va = v3 ; vb = v2 ; break ;
        case 1: vo = v2 ; va = v1 ; vb = v3 ; break ;
        case 2: vo = v3 ; va = v2 ; vb = v1 ; break ;
        }
        delta = deltaAngle(face->angle[tno][ano],face->orig_angle[tno][ano]);
#if ONLY_NEG_AREA_TERM
        if (face->angle[tno][ano] >= 0.0f)
          delta = 0.0f ;
#endif
        delta *= parms->l_angle ;
        VERTEX_EDGE(v_a, vo, va) ; VERTEX_EDGE(v_b, vo, vb) ;
        
        /* this angle's contribution to va */
        V3_CROSS_PRODUCT(v_a, v_n0, v_tmp) ;
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
  VectorFree(&v_c) ;
  VectorFree(&v_d) ;
  VectorFree(&v_tmp) ;
  VectorFree(&v_sum) ;
  VectorFree(&v_n0) ;
  VectorFree(&v_n1) ;

  VectorFree(&v_a_x_n) ;
  VectorFree(&v_b_x_n) ;
  VectorFree(&v_c_x_n) ;
  VectorFree(&v_d_x_n) ;

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
            "neg: %d (%%%2.2f), avgs: %d\n", 
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
  char    fname[100], path[100], *cp ;
  int     vno, nvertices, nfaces, magic, version, tmp, ix, iy, iz, n ;
  VERTEX  *vertex ;
  FILE    *fp ;

  cp = strchr(name, '/') ;
  if (cp)
    strcpy(fname, name) ;     /* path already specified */
  else                        /* no path - use same as was used in MRISread */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/%s.%s", path, 
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", name) ;
  }
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadVertexPosition(%s): could not open file", name));

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
  fread3(&nfaces, fp);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "reading %d vertices and %d faces.\n",nvertices,nfaces);

  if (nvertices != mris->nvertices || nfaces != mris->nfaces)
    ErrorExit(ERROR_BADFILE, "MRISreadVertexPositions(%s): surfaces differ.\n",
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
  MRIScomputeTriangleProperties(mris,0) ;  /* compute areas and normals */
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
  MRIScomputeTriangleProperties(mris, 0) ;/* make sure angles are calculated */
  MRISstoreMetricProperties(mris) ;
  mris->status = old_status ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeTriangleProperties(mris, 0) ;/* make sure angles are calculated */
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
    if (v->ripflag)
      continue ;
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
    if (v->ripflag)
      continue ;
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
  char  *dot, ext[100] ;

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
MRISwritePatchAscii(MRI_SURFACE *mris, char *fname)
{
  FILE    *fp ;
  int     vno, fno, n, nvertices, nfaces, type ;
  VERTEX  *v ;
  FACE    *face ;

  type = mrisFileNameType(fname) ;
#if 0
  if (type == MRIS_ASCII_FILE)
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
#define DISABLE_STDS  0
#if DISABLE_STDS
std = 1.0f ;
#else
    std = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,parms->frame_no+1);
    std = sqrt(std) ;
    if (FZERO(std))
      std = FSMALL ;
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
      std = FSMALL ;
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
      std = FSMALL ;
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
        for (gamma = -degrees ; gamma <= degrees ; gamma += delta)
        {
          MRISsaveVertexPositions(mris, TMP_VERTICES) ;
          MRISrotate(mris, mris, alpha, beta, gamma) ;
          sse = mrisComputeCorrelationError(mris, parms, 0) ;
          if (sse < min_sse)
          {
#if 0
            if (Gdiag & DIAG_SHOW)
              fprintf(stderr, 
                      "\nrotating brain by (%+2.2f, %+2.2f, %+2.2f), "
                      "sse: %2.2f\n",
                      (float)DEGREES(alpha), (float)DEGREES(beta), 
                      (float)DEGREES(gamma), (float)sse) ;
#endif
            mina = alpha ; minb = beta ; ming = gamma ;
            min_sse = sse ;
#if 0
            if (Gdiag & DIAG_SHOW)
              fprintf(stderr, " *** ") ;
#endif
          }
          if (Gdiag & DIAG_SHOW)
            fprintf(stderr, "\r(%+2.2f, %+2.2f, %+2.2f), "
                    "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                    (float)DEGREES(alpha), (float)DEGREES(beta), (float)
                    DEGREES(gamma), (float)DEGREES(mina), 
                    (float)DEGREES(minb), (float)DEGREES(ming),(float)min_sse);

          MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
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

    nx = vertex->nx ; ny = vertex->ny ; nz = vertex->nz ;

    if (vertex->vtotal <= 0)
      continue ;
    

    if (vno == Gdiag_no)
      DiagBreak() ;

    /* fit a quadratic form to the surface at this vertex */
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
   calculate the projection of this vertex onto the local tangent plane 
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
  int     vno, fno, tno ;
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
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
      if (face->area[tno] < 0.0f)
        face->area[tno] = 0.0f ;
  }
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISfindClosestCannonicalVertex(MRI_SURFACE *mris, float x, float y, float z)
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
  VERTEX *vertex ;
  float   kmin, kmax ;

  kmin = 100000.0f ; kmax = -100000.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->curv = MAX(fabs(vertex->k1), fabs(vertex->k2)) ;
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
MRISuseNegCurvature(MRI_SURFACE *mris)
{
  int    vno, fno, tno ;
  VERTEX *vertex ;
  FACE   *f ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->curv = 0 ;
    for (fno = 0 ; fno < vertex->num ; fno++)
    {
      f = &mris->faces[vertex->f[fno]] ;
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
        if (f->area[tno] < 0.0f)
          vertex->curv = 1.0f ;
        
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

#define MAX_MOVEMENT    1.0
#define DELTA_M         0.1
#define NSAMPLES        (2*(MAX_MOVEMENT/DELTA_M)+1)

#define WHALF           (5-1)/2
#define MAX_ITER        25
#define DI_DIST         MAX_MOVEMENT
#define WRITE_ITER      10
#define AVERAGE_ITER    50

int
MRISpositionSurface(MRI_SURFACE *mris, MRI *mri_brain, MRI *mri_wm,
                    float nsigma,int where, float dt)
{
  int    vno, xv, yv, zv, xi, yi, zi, xo, yo, zo, nvox, i, max_v, which = 0,
         new_xv, new_yv, new_zv, xf, yf, zf ;
  VERTEX *v ;
  float  mean, sigma, total_sq, total, delta, max_del, mean_wm, mean_gray,
         dist, min_error, error, min_dist ;
  Real   x, y, z, xw, yw, zw, central_val, val, nx, ny, nz, new_x, new_y,new_z;
  FILE   *fp ;
  MRI    *mri_filled = NULL ;

  if (Gdiag & DIAG_WRITE)
    fp = fopen("position.log", "w") ;
  else
    fp = NULL ;

  MRIvoxelToWorld(mri_wm, 0, 0, 0, &xw, &yw, &zw) ;
  
  /* first compute intensity of local gray/white boundary */
  mean_wm = mean_gray = 0.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
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
          if (val > 0)
          {
            val = (Real)MRIvox(mri_brain, xi, yi, zi) ;  /* use smoothed val */
            total += val ;
            total_sq += val * val ;
            nvox++ ;
          }
        }
      }
    }
    mean = total / (float)nvox ;
    sigma = sqrt(total_sq / (float)nvox - mean*mean) ;
    
    MRISvertexToVoxel(v, mri_wm, &x, &y, &z) ;
    MRIsampleVolume(mri_brain, x, y, z, &central_val) ;
    v->val = mean - nsigma * sigma ;
    mean_gray += v->val ;
    mean_wm += mean ;
  }
  mean_wm /= (float)mris->nvertices ; mean_gray /= (float)mris->nvertices ;

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
  if (Gdiag & DIAG_WRITE)
    fprintf(fp, "mean wm=%2.1f, gray=%2.1f\n", mean_wm, mean_gray) ;
  MRISaverageVals(mris, 64) ;
  MRISaverageVertexPositions(mris, 4) ;
  mrisComputeNormals(mris) ;
  fprintf(stdout, "done.\n") ;
  for (i = 0 ; i < MAX_ITER ; i++)
  {
    /*    DiagHeartbeat((float)i / (float)(MAX_ITER-1)) ;*/
    if (Gdiag & DIAG_WRITE && !(i % WRITE_ITER))  /* write snapshot */
    {
      char fname[100], path[100] ;
      FileNamePath(mris->fname, path) ;
      sprintf(fname, "%s/%s.white%3.3d", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", i) ;
      fprintf(stdout, "writing snapshot %s...", fname) ;
      MRISwrite(mris, fname) ;
      fprintf(stdout, "done.\n") ;
    }
    max_v = -1 ;
    mrisClearGradient(mris) ;

    mri_filled = MRISwriteSurfaceIntoVolume(mris, mri_brain, mri_filled) ;

    /* calculate movement for each vertex in normal direction */
    for (max_del = total = 0.0, vno = which ;vno < mris->nvertices ;vno += 2)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      MRISvertexToVoxel(v, mri_filled, &x, &y, &z) ;
      xf = nint(x) ; yf = nint(y) ; zf = nint(z) ;
      MRISvertexToVoxel(v, mri_wm, &x, &y, &z) ;
      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
      MRIsampleVolume(mri_brain, x, y, z, &central_val) ;
      delta = (central_val - v->val) ;
      if (fabs(delta) > max_del)
      {
        max_del = delta ;
        max_v = vno ;
      }
      total += delta*delta ;

      /* calculate normal vector in volume coordinates */
      nx = xw + (Real)v->nx ; ny = yw + (Real)v->ny ; nz = zw + (Real)v->nz ;
      MRIworldToVoxel(mri_brain, nx, ny, nz, &nx, &ny, &nz) ;

      min_error = 10000.0f ; min_dist = 0.0f ;

      /* find optimal position for this vertex in normal direction */

      /* search 'inwards' */
      for (dist = 0 ; dist >= -MAX_MOVEMENT ; dist -= DELTA_M)
      {
        new_x = x+dist*nx ; new_y = y+dist*ny ; new_z = z+dist*nz ;
        MRIsampleVolume(mri_brain, new_x, new_y, new_z, &val) ;

        /* now check for self-intersection */
        new_x = v->x+dist*v->nx ; new_y = v->y+dist*v->ny ; 
        new_z = v->z+dist*v->nz ;
        MRIworldToVoxel(mri_filled,new_x, new_y, new_z, &new_x,&new_y,&new_z);
        new_xv = nint(new_x) ; new_yv = nint(new_y) ; new_zv = nint(new_z) ;

        /* find voxel coords in .5 mm resolution bitmap volume */
        if ((new_xv != xf) || (new_yv != yf) || (new_zv != zf))
        {
          /* use sub-voxel resolution to test surface self-intersection */
          if (MRItest_bit(mri_filled, new_xv, new_yv, new_zv))
            break ;   /* self-intersection - stop searching */
        }

        error = fabs(v->val - val) ;
        if (error < min_error)
        {
          min_error = error ;
          min_dist = dist ;
        }
      }

      /* now search 'outwards' */
      for (dist = DELTA_M ; dist <= MAX_MOVEMENT ; dist += DELTA_M)
      {
        new_x = x+dist*nx ; new_y = y+dist*ny ; new_z = z+dist*nz ;
        MRIsampleVolume(mri_brain, new_x, new_y, new_z, &val) ;

        new_x = v->x+dist*v->nx ; new_y = v->y+dist*v->ny ; 
        new_z = v->z+dist*v->nz ;
        MRIworldToVoxel(mri_filled,new_x, new_y, new_z, &new_x,&new_y,&new_z);
        new_xv = nint(new_x) ; new_yv = nint(new_y) ; new_zv = nint(new_z) ;

        /* find voxel coords in .5 mm resolution bitmap volume */
        if ((new_xv != xf) || (new_yv != yf) || (new_zv != zf))
        {
          /* use sub-voxel resolution to test surface self-intersection */
          if (MRItest_bit(mri_filled, new_xv, new_yv, new_zv))
            break ;   /* self-intersection - stop searching */
        }

        error = fabs(v->val - val) ;
        if (error < min_error)
        {
          min_error = error ;
          min_dist = dist ;
        }
      }

      delta = min_dist ;
      
      /* turn off old bit and turn on new one */
      MRISvertexToVoxel(v, mri_filled, &x, &y, &z) ;
      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
      MRIclear_bit(mri_filled, xv, yv, zv) ;
      new_x = v->x + dist*v->nx ;
      new_y = v->y + dist*v->ny ;
      new_z = v->z + dist*v->nz ;
      MRIworldToVoxel(mri_filled, new_x, new_y,new_z,&new_x,&new_y,&new_z);
      new_xv = nint(new_x) ; new_yv = nint(new_y) ; new_zv = nint(new_z) ;
      MRIset_bit(mri_filled, new_xv, new_yv, new_zv) ;
      v->dx = delta * v->nx ;
      v->dy = delta * v->ny ;
      v->dz = delta * v->nz ;
    }

    mean = sqrt(total / (float)mris->nvertices) ;

#if 1
    mrisApplyGradient(mris, dt) ;
#else
    mrisMomentumTimeStep(mris, 0.5, dt, 1e-4, 0) ;
#endif
    if (!(i % AVERAGE_ITER))
    {
      which = !which ;  /* average ones which weren't just moved */
      MRISaverageEveryOtherVertexPositions(mris, 1, which) ;
    }
    mrisComputeNormals(mris) ;
    MRISvertexToVoxel(&mris->vertices[max_v], mri_wm, &x, &y, &z) ;
    MRIsampleVolume(mri_brain, x, y, z, &central_val) ;
    fprintf(stdout, "%3.3d: max delta = %2.1f, mean = %2.3f, max_v = %3.3d, "
            "val = %2.1f, central val = %2.1f\n", i,
            max_del, mean, max_v, mris->vertices[max_v].val,central_val) ;
    if (Gdiag & DIAG_WRITE)
      fprintf(fp, "%3.3d: max delta = %2.1f, mean = %2.3f, max_v = %3.3d, "
            "val = %2.1f, central val = %2.1f\n", i,
            max_del, mean, max_v, mris->vertices[max_v].val,central_val) ;
  }
  /*  MRISaverageVertexPositions(mris, 1) ;*/
  if (Gdiag & DIAG_WRITE)
    fclose(fp) ;

  MRIfree(&mri_filled) ;

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
  Real   x, y, z ;
  int    xv, yv, zv, vno ;
  VERTEX *v ;

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

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    MRISvertexToVoxel(v, mri, &x, &y, &z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    MRIset_bit(mri, xv, yv, zv) ;
  }
  return(mri) ;
}

