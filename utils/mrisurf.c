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

/*---------------------------- STRUCTURES -------------------------*/

/*---------------------------- CONSTANTS -------------------------*/

#define ORIG_AREAS     0
#define CURRENT_AREAS  1
#define AVERAGE_AREAS  0


/*------------------------ STATIC PROTOTYPES -------------------------*/

static int   mrisComputeNormals(MRI_SURFACE *mris) ;
static int   mrisFindNeighbors(MRI_SURFACE *mris) ;
static void  mrisNormalize(float v[3]) ;
static float mrisTriangleArea(MRIS *mris, int fac, int n) ;
static int   mrisNormalFace(MRIS *mris, int fac,int n,float norm[]) ;
static int   mrisReadTransform(MRIS *mris, char *mris_fname) ;
static int   mrisReadBinaryAreas(MRI_SURFACE *mris, char *mris_fname) ;
/*static int   mrisReadFieldsign(MRI_SURFACE *mris, char *fname) ;*/
static double mrisComputeSSE(MRI_SURFACE *mris, float l_area,float l_angle,
                             float l_curv);
static double mrisComputeError(MRI_SURFACE *mris, 
                               float l_area,float l_angle, float l_curv,
                                double *parea_error, double *pangle_error,
                               double *pcurv_error,
                               int *pn);
static int   mrisScaleEllipsoidArea(MRI_SURFACE *mris) ;
static int   mrisCountNegativeTriangles(MRI_SURFACE *mris) ;
static int   mrisAverageGradients(MRI_SURFACE *mris, int num_avgs) ;
static int   mrisIntegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static float mrisLineMinimize(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static float deltaAngle(float angle1, float angle2) ;
static int   mrisOrientEllipsoid(MRI_SURFACE *mris) ;
static int   mrisOrientPlane(MRI_SURFACE *mris) ;
#if 0
static int   mrisFindPoles(MRIS *mris) ;
static int   mrisComputeEllipsoidProperties(MRI_SURFACE *mris) ;
#endif
#if AVERAGE_AREAS
static int   mrisAverageAreas(MRI_SURFACE *mris, int num_avgs, int which) ;
#endif
#if 0
static int   mrisComputeBoundaryNormals(MRI_SURFACE *mris) ;
static int   mrisSmoothBoundaryNormals(MRI_SURFACE *mris, int niter) ;
#endif
static int   mrisRipFaces(MRI_SURFACE *mris) ;
static int   transform(float *xptr, float *yptr, float *zptr, 
                       float nx, float ny, float nz, float d) ;

static int   mrisComputeTangentPlanes(MRI_SURFACE *mris) ;
static int   mrisComputeCurvatureTerm(MRI_SURFACE *mris, 
                                      INTEGRATION_PARMS *parms) ;
static int   mrisUpdateSurface(MRI_SURFACE *mris) ;
static int   mrisRemoveLink(MRI_SURFACE *mris, int vno1, int vno2) ;
static int   mrisRemoveEdge(MRI_SURFACE *mris, int vno1, int vno2) ;
static int   mrisRemoveFace(MRI_SURFACE *mris, int fno) ;
static int   mrisCountValidLinks(MRI_SURFACE *mris, int vno1, int vno2) ;

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
              imnr, imnr0, imnr1 ;
  float       x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  FILE        *fp ;
  VERTEX      *vertex ;
  FACE        *face ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL,(ERROR_NOFILE,"MRISread(%s): could not open file",fname));

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
  strcpy(mris->fname, fname) ;
  if (strstr(fname, "rh"))
    mris->hemisphere = RIGHT_HEMISPHERE ;
  else
    mris->hemisphere = LEFT_HEMISPHERE ;

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
    vertex->ox = vertex->x ; vertex->oy = vertex->y ; vertex->oz = vertex->z ;
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

  /*  fprintf(stderr, "surfer: imnr0=%d, imnr1=%d\n",imnr0,imnr1);*/
  for (fno=0;fno<nfaces;fno++)
  {
    for (n=0;n<4;n++)
    {
      fread3(&mris->faces[fno].v[n],fp);
      if (version < 0)   /* new surface format */
        mris->vertices[mris->faces[fno].v[n]].num++;
    }
  }
  fclose(fp);
  if (version<0)
  {
    for (vno = 0 ; vno< nvertices ; vno++)
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
    for (fno = 0 ; fno < nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        mris->vertices[face->v[n]].f[mris->vertices[face->v[n]].num++] = fno;
    }
  }


  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < nvertices ; vno++)
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
#if 0
  sub_num=vertex_index;
#endif
  mrisFindNeighbors(mris);
  mrisComputeNormals(mris);
  mrisReadTransform(mris, fname) ;
/*  mrisComputeTangentPlanes(mris) ;*/
  if (MRISreadBinaryCurvature(mris, fname) != NO_ERROR)
    return(NULL) ;

  if (mrisReadBinaryAreas(mris, fname) != NO_ERROR)
    return(NULL) ;

  /*  mrisFindPoles(mris) ;*/
#if 0
  if (mrisReadTriangleProperties(mris, fname))
    mrisComputeTriangleProperties(mris) ;
  MRIScomputeFaceAreas(mris) ;
#endif

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
  int k,n;
  float x,y,z;
  FILE *fp;

  fp = fopen(fname,"w");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE,"MRISwrite(%s): can't create file %s\n",fname));
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

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (mris->vertices[vno].f)
      free(mris->vertices[vno].f) ;
    if (mris->vertices[vno].n)
      free(mris->vertices[vno].n) ;
  }

  if (mris->vertices)
    free(mris->vertices) ;
  if (mris->faces)
    free(mris->faces) ;
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisFindNeighbors(MRI_SURFACE *mris)
{
  int n0,n1,i,k,m,n;
  face_type *f;
  vertex_type *v;
  int vtmp[9];

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "finding surface neighbors...") ;

  for (k=0;k<mris->nvertices;k++)
  {
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
    mris->vertices[k].v = (int *)calloc(mris->vertices[k].vnum,sizeof(int));
    if (!mris->vertices[k].v)
      ErrorExit(ERROR_NOMEMORY, 
                "mrisFindNeighbors: could not allocate nbr array") ;

    for (i=0;i<v->vnum;i++)
    {
      v->v[i] = vtmp[i];
    }
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
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "done.\n") ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          no triangle area in msurfer, no explodeflag here 
------------------------------------------------------*/
static int
mrisComputeNormals(MRI_SURFACE *mris) 
{
  int k,n;
  vertex_type *v;
  face_type *f;
  float norm[3],snorm[3], len;

#if 0
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "computing surface normals...") ;
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

    if (v->origarea<0)        /* has never been set */
      v->origarea = v->area;

    len = sqrt(snorm[0]*snorm[0] + snorm[1]*snorm[1] + snorm[2]*snorm[2]) ;
    if (!FZERO(len))
    {
      v->nx = snorm[0];
      v->ny = snorm[1];
      v->nz = snorm[2];
    }
  }
#if 0
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;
#endif
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
  int n0,n1;
  face_type *f;
  float v0[3],v1[3];

  n0 = (n==0)?3:n-1;
  n1 = (n==3)?0:n+1;
  f = &mris->faces[fac];
  v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
  v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
  v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
  v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
  v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
  v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
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
  char transform_fname[100], fpref[100] ;

  FileNamePath(mris_fname, fpref) ;
  sprintf(transform_fname, "%s/../mri/transforms/talairach.xfm", fpref) ;
  if (input_transform_file(transform_fname, &mris->transform) != OK)
        ErrorPrintf(ERROR_NO_MEMORY, 
                    "mrisReadTransform: could not read xform file '%s'\n", 
                    transform_fname) ;
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
  char   fname[100], fpref[100], hemi[20], *cp ;
  
  cp = strrchr(mris_fname, '.') ;
  if (!cp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "MRISreadBinaryCurvature(%s): could not scan "
                 "hemisphere from fname", mris_fname)) ;
  strncpy(hemi, cp-2, 2) ; hemi[2] = 0 ;
  FileNamePath(mris_fname, fpref) ;
  sprintf(fname, "%s/%s.curv", fpref, hemi) ;
  return(MRISreadCurvatureFile(mris, fname)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISreadCurvatureFile(MRI_SURFACE *mris, char *fname)
{
  int    k,i,vnum,fnum;
  float  curv, curvmin, curvmax;
  FILE   *fp;
  
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
  char  fname[100], fpref[100], hemi[20], *cp ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "reading area file...") ;

  cp = strrchr(mris_fname, '.') ;
  if (!cp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "mrisReadBinaryAreas(%s): could not scan "
                 "hemisphere from fname", mris_fname)) ;
  strncpy(hemi, cp-2, 2) ; hemi[2] = 0 ;
  FileNamePath(mris_fname, fpref) ;
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
#endif
#if 0
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
          Perform a projection onto an ellipsoid moving each
          point on the cortical surface to the closest ellipsoidal
          coordinate.
------------------------------------------------------*/
#if 0
#define DEFAULT_A  44.0f
#define DEFAULT_B  122.0f
#define DEFAULT_C  70.0f
#else
#define DEFAULT_A  122.0f
#define DEFAULT_B  122.0f
#define DEFAULT_C  122.0f
#endif

#define MAX_DIM    DEFAULT_B

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

#if 0
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "projecting onto ellipsoid...") ;
#endif

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
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;
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
    if (vdst->vnum)
    {
      vdst->v = (int *)calloc(vdst->vnum, sizeof(int)) ;
      if (!vdst->v)
        ErrorExit(ERROR_NO_MEMORY, "MRISclone: could not allocate %d nbrs",
                  vdst->vnum) ;
      for (n = 0; n < vdst->vnum; n++)
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
    ErrorReturn(mris_dst, 
                (ERROR_BADPARM, "MRIStalairachTransform: no xform loaded")) ;

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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/

MRI_SURFACE *
MRISunfold(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  char    host_name[100], *cp ;
  int     base_averages, n_averages, done, steps, total_steps ;
  float   base_tol ;
  INTEGRATION_PARMS _parms ;

  if (!parms)   /* use all defaults */
  {
    parms = &_parms ;
    parms->tol = TOL ;
    parms->n_averages = N_AVERAGES ;
    strcpy(parms->base_name, "unfold") ;
    parms->l_angle = L_ANGLE ;
    parms->l_area = L_AREA ;
    parms->niterations = NITERATIONS ;
    parms->write_iterations = WRITE_ITERATIONS ;
  }
  else  /* put in default parameters if not otherwise specified */
  {
    if (parms->tol <= 0.0f)
      parms->tol = TOL ;
    if (parms->n_averages < 0)
      parms->n_averages = N_AVERAGES ;
    if (!strlen(parms->base_name))
      strcpy(parms->base_name, "unfold") ;
    if (parms->l_angle < 0.0f)
      parms->l_angle = L_ANGLE ;
    if (parms->l_area < 0.0f)
      parms->l_area = L_AREA ;
    if (parms->niterations <= 0)
      parms->niterations = NITERATIONS ;
    if (parms->write_iterations < 0)
      parms->write_iterations = WRITE_ITERATIONS ;
  }

  base_averages = parms->n_averages ;
  cp = getenv("HOST") ;
  if (cp)
    strcpy(host_name, cp) ;
  else
    strcpy(host_name, "unknown") ;


  if (Gdiag & DIAG_SHOW)
  {
    char fname[100] ;

    sprintf(fname, "%s.out", parms->base_name) ;
    parms->fp = fopen(fname, "w") ;
  }
  printf("tol=%2.1e, host=%5.5s, nav=%d, l_angle=%2.2f, l_area=%2.2f\n", 
         parms->tol, host_name, parms->n_averages, parms->l_angle, 
         parms->l_area) ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(parms->fp, 
            "tol=%2.1e, host=%5.5s, nav=%d, l_angle=%2.2f, l_area=%2.2f\n", 
            parms->tol, host_name, parms->n_averages, parms->l_angle, 
            parms->l_area) ;
    fflush(parms->fp) ;
  }

  parms->start_t = 0 ;
  parms->niterations = 1 ;
  base_tol = parms->tol ;
  do
  {
    done = 0 ;
    for (total_steps = 0, n_averages = base_averages ; !done ; n_averages /= 2)
    {
      parms->tol = base_tol * sqrt((float)n_averages+1.0f) ;
#if 0
      fprintf(stderr, "integrating at scale = %d, tol = %2.1e\n", n_averages,
              parms->tol) ;
#endif
      parms->n_averages = n_averages ;     /* # of averages == scale */
      steps = mrisIntegrate(mris, parms) ;
      parms->start_t += steps ;
      total_steps += steps ;
      done = n_averages == 0 ;    /* finished integrating at smallest scale */
    }
  } while (total_steps > 0) ;

  if (Gdiag & DIAG_SHOW)
    fclose(parms->fp) ;

  switch (parms->projection)
  {
  case ELLIPSOID_PROJECTION:
    MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
    break ;
  case PROJECT_PLANE:
    mrisOrientPlane(mris) ;
    break ;
  default:
    mrisUpdateSurface(mris) ;
    break ;
  }
  return(mris) ;
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
mrisIntegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  char    base_name[100] ;
  int     vno, t, fno, imin, write_iterations, num, negative, n_averages, 
          niterations, ano, tno ;
  VERTEX  *v0, *v1, *v2, *v3, *va, *vb, *vo ;
  VECTOR  *v_a, *v_b, *v_c, *v_d, *v_a_x_n, *v_b_x_n, *v_c_x_n, *v_d_x_n,
          *v_n0, *v_n, *v_n1, *v_tmp, *v_sum ;
  FACE    *face ;
  double  sse, old_sse, angle_sse, area_sse, curv_sse ;
  float   delta_t, orig_area0, area0, orig_area1, area1, 
          delta0, delta1, l_area, l_corr, l_angle, delta, len ;

  sprintf(base_name, "surf/%s.%s", mris->hemisphere == LEFT_HEMISPHERE ?
              "lh" : "rh", parms->base_name) ;
  n_averages = parms->n_averages ;
  l_angle = parms->l_angle ;
  l_area = parms->l_area ;
  l_corr = parms->l_corr ;
  write_iterations = parms->write_iterations ;
  niterations = parms->niterations ;

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

  switch (parms->projection)
  {
  case ELLIPSOID_PROJECTION:
    MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
    break ;
  case PROJECT_PLANE:
    mrisOrientPlane(mris) ;
    break ;
  default:
    mrisUpdateSurface(mris) ;
    break ;
  }

#if AVERAGE_AREAS
  MRISreadTriangleProperties(mris, mris->fname) ;
  mrisAverageAreas(mris, parms->n_averages, ORIG_AREAS) ;
#endif

  sse = old_sse = mrisComputeError(mris, parms->l_area, parms->l_angle,
                                   parms->l_curv,
                                   &area_sse,&angle_sse, &curv_sse, &num);
  imin = -1 ;
  delta_t = 0.0f ;
  niterations += parms->start_t ;
  if (!parms->start_t)
  {
    char fname[100] ;

    negative = mrisCountNegativeTriangles(mris) ;
    printf("%d: dt: %2.3f, rms: %2.4f (%2.3f, %2.1f, %2.1f), neg: %d (%2.1f)"
           ", avgs: %d\n", parms->start_t, delta_t, 
           sqrt(sse/(float)num), sqrt(area_sse/(float)num), 
           DEGREES(sqrt(angle_sse/(float)(num*3))), curv_sse,
           negative, mris->neg_area, n_averages);
    if (Gdiag & DIAG_SHOW)
    {
      if (Gdiag & DIAG_SHOW)
        fprintf(parms->fp, "%d: dt: %2.3f, rms: %2.4f (%2.3f, %2.1f, %2.1f), "
                "neg: %d (%2.1f) avgs: %d\n", parms->start_t, 
                delta_t, sqrt(sse/(float)num), sqrt(area_sse/(float)num), 
                DEGREES(sqrt(angle_sse/(float)(num*3))), curv_sse,
                negative, mris->neg_area, 
                n_averages);
      fflush(parms->fp) ;
    }
    fflush(stdout) ;
    if (Gdiag & DIAG_WRITE)
    {
      sprintf(fname, "%s%4.4d", base_name, 0) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "writing %s...", fname) ;
      if (parms->projection == PROJECT_PLANE)
        MRISwritePatch(mris, fname) ;
      else
        MRISwrite(mris, fname) ;
      if (!FZERO(parms->l_area) && DIAG_VERBOSE_ON)
      {
        sprintf(fname, "%s%4.4d.area_error", base_name, 0) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, " %s...", fname) ;
        MRISwriteAreaError(mris, fname) ;
      }

      if (!FZERO(parms->l_angle) && DIAG_VERBOSE_ON)
      {
        sprintf(fname, "%s%4.4d.angle_error", base_name, 0) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, " %s...", fname) ;
        MRISwriteAngleError(mris, fname) ;
      }
      if (Gdiag & DIAG_SHOW)
      {
        fprintf(stderr, "done.\n") ;
        fflush(stderr) ;
      }
    }
  }

  for (t = parms->start_t ; t < niterations ; t++)
  {
#if AVERAGE_AREAS
    mrisAverageAreas(mris, parms->n_averages, CURRENT_AREAS) ;
#endif

    if (!FZERO(parms->l_curv))
      MRIScomputeMeanCurvature(mris) ;

    /* clear old deltas */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v0 = &mris->vertices[vno] ;
      v0->dx = v0->dy = v0->dz = 0.0f ;
    }

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
      orig_area0 = face->orig_area[0] ; area0 = face->area[0] ;
      orig_area1 = face->orig_area[1] ; area1 = face->area[1] ;
      delta0 = parms->l_area * (area0 - orig_area0) ; 
      delta1 = parms->l_area * (area1 - orig_area1) ;
      V3_CROSS_PRODUCT(v_a, v_n0, v_a_x_n) ;
      V3_CROSS_PRODUCT(v_b, v_n1, v_b_x_n) ;
      V3_CROSS_PRODUCT(v_c, v_n1, v_c_x_n) ;
      V3_CROSS_PRODUCT(v_d, v_n0, v_d_x_n) ;

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
      for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
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

    if (!FZERO(parms->l_curv))
      mrisComputeCurvatureTerm(mris, parms) ;
#if !AVERAGE_AREAS
    mrisAverageGradients(mris, n_averages) ;
#endif
    parms->t = t ;
    delta_t = mrisLineMinimize(mris, parms) ;

for (vno = 0 ; vno < mris->nvertices ; vno++)
{
  v0 = &mris->vertices[vno] ;
  if (v0->y < -100 || v0->y > 100 || v0->z < -100 || v0->z > 100)
    DiagBreak() ;
}
    if (delta_t < parms->tol)   /* reached the minimum */
      break ;

    /* only print stuff out if we actually took a step */
    sse = mrisComputeError(mris,parms->l_area,parms->l_angle, parms->l_curv,
                           &area_sse, &angle_sse, &curv_sse, &num);
    if (FZERO(sse) || (fabs((sse-old_sse)/(old_sse)) < parms->tol))
      break ;

    old_sse = sse ;
    negative = mrisCountNegativeTriangles(mris) ;
    printf("%d: dt: %2.3f, rms: %2.4f (%2.3f, %2.1f, %2.1f), neg: %d (%2.1f), "
           "avgs: %d\n", t+1,   delta_t, 
           (float)sqrt(sse/(float)num), 
           (float)sqrt(area_sse/(float)num), 
           (float)DEGREES(sqrt(angle_sse/(float)(num*3))), curv_sse, 
           negative, mris->neg_area, n_averages) ;
    fflush(stdout) ;
    if (Gdiag & DIAG_SHOW)
    {
      fprintf(parms->fp, 
             "%d: dt: %2.3f, rms: %2.4f (%2.3f, %2.1f, %2.1f), neg: %d (%2.1f)"
              ", avgs: %d\n", t+1, delta_t, 
              (float)sqrt(sse/(float)num), 
              (float)sqrt(area_sse/(float)num), 
              (float)DEGREES(sqrt(angle_sse/(float)(num*3))),curv_sse, 
              negative,mris->neg_area, n_averages) ;
      fflush(parms->fp) ;
    }

    if ((write_iterations > 0) &&!((t+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
    {
      char fname[100] ;
      sprintf(fname, "%s%4.4d", base_name, t+1) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "writing %s...", fname) ;
      switch (parms->projection)
      {
      case PROJECT_PLANE:
        MRISwritePatch(mris, fname) ;
        break ;
      default:
        MRISwrite(mris, fname) ;
      }
      if (DIAG_VERBOSE_ON)
      {
        sprintf(fname, "%s%4.4d.area_error", base_name, t+1) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr, " %s...", fname) ;
        MRISwriteAreaError(mris, fname) ;
      }
      if (Gdiag & DIAG_SHOW)
      {
        fprintf(stderr, "done.\n") ;
        fflush(stderr) ;
      }
    }
  }

  switch (parms->projection)
  {
  case ELLIPSOID_PROJECTION:
    MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f) ;
    break ;
  case PROJECT_PLANE:
    mrisOrientPlane(mris) ;
    break ;
  default:
    mrisUpdateSurface(mris) ;
    break ;
  }

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

  return(t-parms->start_t) ;  /* return actual # of steps taken */
}
/*-----------------------------------------------------
        Parameters:


        Returns value:

        Description
------------------------------------------------------*/
static double
mrisComputeError(MRI_SURFACE *mris, float l_area,float l_angle,float l_curv,
                  double *parea_error, double *pangle_error, 
                 double *pcurv_error, int *pn)
{
  double  sse, sse_area, sse_angle, sse_curv, delta ;
  int     ano, fno, tno, n, vno, n1 ;
  FACE    *face ;
  VERTEX  *v ;

  sse_angle = sse_area = 0.0 ;
  for (n = fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      n++ ;
      delta = (double)(face->area[tno] - face->orig_area[tno]) ;
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

  sse_curv = 0.0 ;
  for (n1 = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    n1++ ;
    sse_curv += v->curv * v->curv ;
  }

  *pn = n ;
  *pcurv_error = sse_curv/(float)n1 ;
  *parea_error = sse_area ;
  *pangle_error = sse_angle ;
  sse = l_area * sse_area + l_angle * sse_angle + l_curv * sse_curv ; ;
  return(sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
mrisComputeSSE(MRI_SURFACE *mris, float l_area, float l_angle, float l_curv)
{
  double  sse, sse_area, sse_angle, delta, sse_curv ;
  int     ano, tno, fno, vno ;
  FACE    *face ;
  VERTEX  *v ;

  sse_angle = sse_area = 0.0 ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      delta = (double)(face->area[tno] - face->orig_area[tno]) ;
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

  sse_curv = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    sse_curv += v->curv * v->curv ;
  }

  sse = l_area * sse_area + l_angle * sse_angle + l_curv * sse_curv ;
  return(sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisScaleEllipsoidArea(MRI_SURFACE *mris)
{
  float   area_scale ;
  int     tno, vno, fno ;
  VERTEX  *v ;
  FACE    *face ;

  area_scale = mris->orig_area / mris->total_area ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    
    /* scale the area by the ratio of the ellipsoid area to that of the
       original surface.
       */
    v->area *= area_scale ;
  }

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
      face->area[tno] *= area_scale ;
  }

  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Find the 3 poles (temporal, occipital, and frontal) of
          the cortical surface.
------------------------------------------------------*/
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
  int     fno, tno, ano ;
  FACE    *face ;

  MRIScomputeTriangleProperties(mris) ;  /* recompute areas and normals */
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

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisCountNegativeTriangles(MRI_SURFACE *mris)
{
  int     tno, fno, negative ;
  FACE    *face ;

  negative = 0 ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
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
  MRIScomputeTriangleProperties(mris) ;  /* recompute areas and normals */
  mrisOrientEllipsoid(mris) ;      /* orient the normals and angles */
  mrisScaleEllipsoidArea(mris) ;   /* scale it to have same area as orig. */
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
  int    i, vno, vnb, *pnb ;
  float  dx, dy, dz, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < num_avgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      dx = v->dx ; dy = v->dy ; dz = v->dz ;
      pnb = v->v ;
      for (num = 0.0f, vnb = 0 ; vnb < v->vnum ; vnb++)
      {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        num++ ;
        dx += vn->dx ; dy += vn->dy ; dz += vn->dz ;
      }
      num++ ;
      v->odx = dx / num ;
      v->ody = dy / num ;
      v->odz = dz / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->dx = v->odx ; v->dy = v->ody ; v->dz = v->odz ;
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
  int     fno, tno, ano ;
  FACE    *face ;
  FILE    *fp;
  char    fname[100], fpref[100], hemi[20], *cp ;

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

  /* compute the area of all the triangles */
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
  for (fno=0;fno<mris->nfaces;fno++)
  {
    face = &mris->faces[fno] ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        putf(face->angle[tno][ano], fp) ;
    }
  }

  fclose(fp);

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
MRIScomputeTriangleProperties(MRI_SURFACE *mris)
{
  VECTOR  *v_a, *v_b, *v_c, *v_d, *v_n ;
  VERTEX  *v0, *v1, *v2, *v3, *va, *vb, *vo ;
  FACE    *face ;
  int     tno, fno, ano, min_fno, min_tno  ;
  float   area, angle, dot, cross, min_p, dz ;
  static  int first = 1 ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_c = VectorAlloc(3, MATRIX_REAL) ;
  v_d = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */

  min_p = 1000.0f ; min_fno = min_tno = -1 ;
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

    /* compute metric properties of second triangle */
    V3_CROSS_PRODUCT(v_c, v_b, v_n) ;
    face->area[1] = area ;
    V3_NORMALIZE(v_n, v_n) ;             /* make it a unit vector */
    face->nx[1] = V3_X(v_n); face->ny[1] = V3_Y(v_n); face->nz[1] = V3_Z(v_n);
    mris->total_area += area ;

    /* now compute angles */
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      VECTOR_LOAD(v_n, face->nx[tno], face->ny[tno], face->nz[tno]) ;
      if ((V3_X(v_n) < V3_Y(v_n)) && (V3_X(v_n) < V3_Z(v_n)))
        dz = fabs(V3_X(v_n)) ;
      else if (V3_Y(v_n) < V3_Z(v_n))
        dz = fabs(V3_Y(v_n)) ;
      else
        dz = fabs(V3_Z(v_n)) ;
      if (dz < min_p)
      {
        min_p = dz ;
        min_fno = fno ;
        min_tno = tno ;
      }
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
                  fno,tno,ano,DEGREES(angle)) ;
#endif
      }
    }
  }

  if (first)
  {
    first = 0 ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "max planar triangle at face %d, tri %d (%2.3f)\n", 
              min_fno, min_tno, min_p) ;
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
#define MAX_MM    MAX_DIM/10.0f
#define MAX_DT    1000000.0f

#if 0
#define MOMENTUM  0.85

static float
mrisLineMinimize(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  float   delta_t ;
  int     vno ;
  VERTEX  *vertex ;

  delta_t = 1024.0f*parms->tol ;
  delta_t = 16.0f*parms->tol ;

  switch (parms->projection)
  {
  case ELLIPSOID_PROJECTION:
    MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
    break ;
  case PROJECT_PLANE:
    mrisOrientPlane(mris) ;
    break ;
  default:
    mrisUpdateSurface(mris) ;
    break ;
  }

  /* pick starting step size */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    vertex->ox = vertex->x ; vertex->oy = vertex->y ; vertex->oz = vertex->z;
    vertex->odx = (1.0-MOMENTUM)*delta_t * vertex->dx + MOMENTUM*vertex->odx ;
    vertex->ody = (1.0-MOMENTUM)*delta_t * vertex->dy + MOMENTUM*vertex->ody ;
    vertex->odz = (1.0-MOMENTUM)*delta_t * vertex->dz + MOMENTUM*vertex->odz ;
    vertex->x += vertex->odx ;
    vertex->y += vertex->ody ;
    vertex->z += vertex->odz ;
    if (!finite(vertex->x) || !finite(vertex->y) || !finite(vertex->z))
      DiagBreak() ;
  }
  switch (parms->projection)
  {
  case ELLIPSOID_PROJECTION:
    MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
    break ;
  case PROJECT_PLANE:
    mrisOrientPlane(mris) ;
    break ;
  default:
    mrisUpdateSurface(mris);
    break ;
  }

  return(delta_t) ;
}
#else
static float
mrisLineMinimize(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  char    fname[100] ;
  FILE    *fp = NULL ;
  double  starting_sse, sse, min_sse, max_dt ;
  float   delta_t, total_delta, min_delta, max_delta, dx, dy, dz, mag,
          grad ;
  int     vno, done = 0, increasing ;
  VERTEX  *vertex ;
  static float last_max_delta_t = 0.0f ;

  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "%s%4.4d.dat", FileName(parms->base_name), parms->t+1);
    fp = fopen(fname, "w") ;
  }

  /*  delta_t = 128.0f*parms->tol ;*/

  switch (parms->projection)
  {
  case ELLIPSOID_PROJECTION:
    MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
    break ;
  case PROJECT_PLANE:
    mrisOrientPlane(mris) ;
    break ;
  default:
    mrisUpdateSurface(mris) ;
    break ;
  }
  min_sse = starting_sse = 
    mrisComputeSSE(mris, parms->l_area, parms->l_angle, parms->l_curv);

  max_delta = 0.0f ;
  grad = 0.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    dx = vertex->dx ; dy = vertex->dy ; dz = vertex->dz ;
    mag = sqrt(dx*dx+dy*dy+dz*dz) ;
    grad += dx*dx+dy*dy+dz*dz ;
    if (mag > max_delta)
      max_delta = mag ;
  }
  grad = sqrt(grad) ;

  switch (parms->projection)
  {
  case PROJECT_ELLIPSOID:
    max_dt = MAX_MM / max_delta ;
    break ;
  case NO_PROJECTION:
    max_dt = 100.0*MAX_MM / max_delta ;  /* allow at most only 2mm steps */
    break ;
  default:
  case PROJECT_PLANE:
    max_dt = 100*MAX_MM / max_delta ;
    break ;
  }

  /* write out some data on supposed quadratic form */
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    float delta, predicted_sse ;
    FILE  *fp2 ;

    sprintf(fname, "nn%s%4.4d.dat",FileName(parms->base_name), parms->t+1) ;
    fp2 = fopen(fname, "w") ;

    delta = max_dt / 100.0f ;
    for (delta_t = delta ; delta_t <= max_dt ; delta_t += delta)
    {
      predicted_sse = starting_sse - grad * delta_t ;
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        vertex = &mris->vertices[vno] ;
        if (vertex->ripflag)
          continue ;
        vertex->ox = vertex->x ; vertex->oy = vertex->y ; 
        vertex->oz = vertex->z;
        vertex->x += delta_t * vertex->dx ;
        vertex->y += delta_t * vertex->dy ;
        vertex->z += delta_t * vertex->dz ;
      }
      switch (parms->projection)
      {
      case ELLIPSOID_PROJECTION:
        MRISupdateEllipsoidSurface(mris) ;
        break ;
      case PROJECT_PLANE:
        mrisOrientPlane(mris) ;
        break ;
      default:
        mrisUpdateSurface(mris) ;
        break ;
      }
        
      sse = mrisComputeSSE(mris, parms->l_area, parms->l_angle,parms->l_curv) ;
      fprintf(fp2, "%f  %f  %f\n", delta_t, sse, predicted_sse) ;
      switch (parms->projection)
      {
      case ELLIPSOID_PROJECTION:
        fflush(fp2) ;
        MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
        break ;
      case PROJECT_PLANE:
        mrisOrientPlane(mris) ;
        break ;
      default:
        mrisUpdateSurface(mris);
        break ;
      }
      sse = mrisComputeSSE(mris, parms->l_area, parms->l_angle,parms->l_curv) ;
      fprintf(fp, "%f  %f  %f\n", delta_t, sse, predicted_sse) ;
      fflush(fp) ;
      for (vno = 0 ; vno < mris->nvertices ; vno++)  
      {
        vertex = &mris->vertices[vno] ;
        if (vertex->ripflag)
          continue ;
        vertex->x = vertex->ox ; vertex->y = vertex->oy ; vertex->z=vertex->oz;
      }
      switch (parms->projection)
      {
      case ELLIPSOID_PROJECTION:
        MRISupdateEllipsoidSurface(mris) ;
        break ;
      case PROJECT_PLANE:
        mrisOrientPlane(mris) ;
        break ;
      default:
        mrisUpdateSurface(mris);
        break ;
      }
    }
  }



#if 0
  if (max_dt > MAX_DT)
    max_dt = MAX_DT ;
#endif

  /* pick starting step size */
  min_delta = 0.0f ; /* to get rid of compiler warning */
  for (delta_t = parms->tol ; delta_t < max_dt ; delta_t *= 10.0f)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
        continue ;
      vertex->ox = vertex->x ; vertex->oy = vertex->y ; vertex->oz = vertex->z;
      vertex->x += delta_t * vertex->dx ;
      vertex->y += delta_t * vertex->dy ;
      vertex->z += delta_t * vertex->dz ;
      if (!finite(vertex->x) || !finite(vertex->y) || !finite(vertex->z))
        DiagBreak() ;
    }
    switch (parms->projection)
    {
    case ELLIPSOID_PROJECTION:
      MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
      break ;
    case PROJECT_PLANE:
      mrisOrientPlane(mris) ;
      break ;
    default:
      mrisUpdateSurface(mris);
      break ;
    }
    sse = mrisComputeSSE(mris, parms->l_area, parms->l_angle,parms->l_curv) ;
    if (sse <= min_sse)   /* new minimum found */
    {
      min_sse = sse ;
      min_delta = delta_t ;
    }

    /* undo step */
    for (vno = 0 ; vno < mris->nvertices ; vno++)  
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
        continue ;
      vertex->x = vertex->ox ; vertex->y = vertex->oy ; vertex->z=vertex->oz;
      if (!finite(vertex->x) || !finite(vertex->y) || !finite(vertex->z))
        DiagBreak() ;
    }
    switch (parms->projection)
    {
    case ELLIPSOID_PROJECTION:
      MRISupdateEllipsoidSurface(mris) ;
      break ;
    case PROJECT_PLANE:
      mrisOrientPlane(mris) ;
      break ;
    default:
      mrisUpdateSurface(mris);
      break ;
    }
  }

  delta_t = min_delta ;
  fprintf(stderr,"grad=%2.4f, max_delta = %2.3f, max_dt = %2.1f, "
          "starting dt = %2.3f\n", grad, max_delta, max_dt, delta_t) ;

  /* now search for minimum in gradient direction */
  increasing = 1 ;
  total_delta = 0.0f ;
  min_sse = starting_sse ;
  while (!done)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
        continue ;
      vertex->ox = vertex->x ; vertex->oy = vertex->y ; vertex->oz = vertex->z;
      vertex->x += delta_t * vertex->dx ;
      vertex->y += delta_t * vertex->dy ;
      vertex->z += delta_t * vertex->dz ;
      if (!finite(vertex->x) || !finite(vertex->y) || !finite(vertex->z))
        DiagBreak() ;
    }
    switch (parms->projection)
    {
    case ELLIPSOID_PROJECTION:
      MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
      break ;
    case PROJECT_PLANE:
      mrisOrientPlane(mris) ;
      break ;
    default:
      mrisUpdateSurface(mris);
      break ;
    }
    sse = mrisComputeSSE(mris, parms->l_area, parms->l_angle, parms->l_curv) ;
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
      last_max_delta_t = total_delta ;   /* best total time step */
    }
    else                 /* error increased - undo it and decrease time step */
    {
      if (increasing)    /* went too far - reverse delta_t change */
        increasing = 0 ;

      for (vno = 0 ; vno < mris->nvertices ; vno++)  /* undo step */
      {
        vertex = &mris->vertices[vno] ;
        if (vertex->ripflag)
          continue ;
        vertex->x = vertex->ox ; vertex->y = vertex->oy ; vertex->z=vertex->oz;
        if (!finite(vertex->x) || !finite(vertex->y) || !finite(vertex->z))
          DiagBreak() ;
      }
      switch (parms->projection)
      {
      case ELLIPSOID_PROJECTION:
        MRISupdateEllipsoidSurface(mris) ;
        break ;
      case PROJECT_PLANE:
        mrisOrientPlane(mris) ;
        break ;
      default:
        mrisUpdateSurface(mris);
        break ;
      }
    }
    if (increasing)    /* increase time step and search further out */
      delta_t *= 2.0f ;
    else               /* decreasing - reduce time step */
      delta_t *= 0.5f ;
    done = delta_t < parms->tol ;
  }

  if (total_delta < parms->tol)    /* couldn't go anywhere - reset last_max */
    last_max_delta_t = 0.0f ;

  if (Gdiag & DIAG_WRITE)
    fclose(fp) ;

  return(total_delta) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static float
deltaAngle(float angle1, float angle2)
{
  float   delta ;

  delta = angle1 - angle2 ;
  if (delta > M_PI)
    delta = 2.0 * M_PI - delta ;
  else if (delta < -M_PI)
    delta = -2.0f * M_PI - delta ;

  return(delta) ;
}
#if 0
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

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  for (k=0;k<mris_dst->nvertices;k++) 
  {
    v = &mris_dst->vertices[k];
    if (v->ripflag)
      continue ;
    v->x *= scale ; v->y *= scale ; v->z *= scale ;
  }
  return(mris_dst) ;
}
#endif

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISwriteCurvature(MRI_SURFACE *mris, char *fname)
{
  int    k,i ;
  float  curv;
  /*  char   fname[100], fpref[100], hemi[20], *cp ;*/
  FILE   *fp;
  
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing curvature file %s...", fname) ;

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
    i = curv * 100.0 ;
    fwrite2((int)i,fp);
  }
  fclose(fp);
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;
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
#define DEBUG_VNO 79221
#define DEBUG_U  -10  /* -10 25*/
#define DEBUG_V  0  /* 2 */

MRI_SURFACE *
MRISrotate(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, float dphi, 
                float dtheta)
{
  int      vno ;
  VERTEX   *vertex ;
  float    phi, theta, x, y, z, a, b, c ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  if (FZERO(mris_src->a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }
  else
  {
    a = mris_src->a ;
    b = mris_src->b ;
    c = mris_src->c ;
  }
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    if (vno == DEBUG_VNO)
      DiagBreak() ;

    vertex = &mris_src->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2.0f * M_PI + theta ;  /* make it 0 --> 2*PI */
    phi = atan2(sqrt(c*c-z*z), z) ;
    phi += dphi ;
    theta += dtheta ;
    x = a*sin(phi)*cos(theta) ;
    y = b*sin(phi)*sin(theta) ;
    z = c*cos(phi) ;
    vertex->x = x ; vertex->y = y ; vertex->z = z ;
  }
  return(mris_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/

#define BIG            100000.0

/* something of a hack... */
#define UNFILLED_ELT   0
#define FILLED_ELT     1

MRI_SP *
MRIStoParameterization(MRI_SURFACE *mris, MRI_SP *mrisp, float scale)
{
  float     a, b, c, phi, theta, x, y, z, total, uf, vf, d, du, dv, total_d,
            **distances, sigma, two_sigma_sq ;
  int       vno, u, v, uk, vk, n, unfilled, u0, v0, u1, v1, **filled ;
  VERTEX    *vertex ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "computing parameterization...") ;

  if (!mrisp)
    mrisp = MRISPalloc(mris, scale) ;

  if (FZERO(mris->a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }
  else
  {
    a = mris->a ;
    b = mris->b ;
    c = mris->c ;
  }

  filled = (int **)calloc(U_DIM(mrisp), sizeof(int *)) ;
  distances = (float **)calloc(U_DIM(mrisp), sizeof(float *)) ;
  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int)) ;
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float)) ;
    
    for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
      filled[u][v] = UNFILLED_ELT ;
  }

  sigma = scale / 4.0f ;
  two_sigma_sq = 2*sigma*sigma ;

  /* first calculate total distances to a point in parameter space */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    phi = atan2(sqrt(c*c-z*z), z) ;
    if (phi < RADIANS(1))
      DiagBreak() ;
    if (vno == DEBUG_VNO)
      DiagBreak() ;
    vertex->phi = phi ; vertex->theta = theta ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u0 = floor(uf) ; u1 = ceil(uf+0.00001f) ;
    v0 = floor(vf) ; v1 = ceil(vf+0.00001f) ;
    du = uf - (float)u0 ; dv = vf - (float)v0 ;
    if (u0 < 0)  /* enforce spherical topology  */
      u0 += U_DIM(mrisp) ;
    if (u0 >= U_DIM(mrisp))
      u0 -= U_DIM(mrisp) ;
    if (u1 < 0)  /* enforce spherical topology  */
      u1 += U_DIM(mrisp) ;
    if (u1 >= U_DIM(mrisp))
      u1 -= U_DIM(mrisp) ;
    if (v0 < 0)  /* enforce spherical topology  */
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0)  /* enforce spherical topology  */
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;

    if ((u0 == DEBUG_U || u1 == DEBUG_U) && (v0 == DEBUG_V || v1 == DEBUG_V))
      DiagBreak() ;

#if 0
    d = sqrt(du*du+dv*dv) ;  /* 1-distance to vertex 0,0 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+dv*dv ;  /* 1-distance to vertex 0,0 */
    d = exp(-d/two_sigma_sq) ;
#endif
    if (u0 >= U_DIM(mrisp) || v0 >= V_DIM(mrisp))
      DiagBreak() ;
    filled[u0][v0] = vno ;
    distances[u0][v0] += d ;
    if ((u0 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    d = 1.0 - dv ;         /* distance to v1 */
#if 0
    d = sqrt(du*du+d*d) ;  /* distance to vertex 0,1 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+d*d ;  /* distance to vertex 0,1 */
    d = exp(-d/two_sigma_sq) ;
#endif
    distances[u0][v1] += d ;         /* keep track of total distance */
    filled[u0][v1] = vno ;
    if ((u0 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    d = 1.0 - du ;         /* distance to u1 */
#if 0
    d = sqrt(d*d+dv*dv) ;  /* distance to vertex 1,0 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = d*d+dv*dv ;  /* distance to vertex 1,0 */
    d = exp(-d/two_sigma_sq) ;
#endif

    distances[u1][v0] += d ;         /* keep track of total distance */
    filled[u1][v0] = vno ;
    if ((u1 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    du = 1.0 - du ; dv = 1.0 - dv ;
#if 0
    d = sqrt(du*du+dv*dv) ;  /* 1-distance to vertex 1,1 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+dv*dv ;  /* 1-distance to vertex 1,1 */
    d = exp(-d/two_sigma_sq) ;
#endif
    distances[u1][v1] += d ;         /* keep track of total distance */
    filled[u1][v1] = vno ;
    if ((u1 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;
  }

  if (DEBUG_U >= 0)
    fprintf(stderr,"\ndistance[%d][%d] = %2.3f\n\n",
            DEBUG_U, DEBUG_V, distances[DEBUG_U][DEBUG_V]);

  /* now add in curvatures proportional to their distance from the point */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    phi = atan2(sqrt(c*c-z*z), z) ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u0 = floor(uf) ; u1 = ceil(uf+0.00001f) ;
    v0 = floor(vf) ; v1 = ceil(vf+0.00001f) ;
    du = uf - (float)u0 ; dv = vf - (float)v0 ;
    if ((u0 == DEBUG_U || u1 == DEBUG_U) && (v0 == DEBUG_V || v1 == DEBUG_V))
      DiagBreak() ;
    if (u0 < 0)  /* enforce spherical topology  */
      u0 += U_DIM(mrisp) ;
    if (u0 >= U_DIM(mrisp))
      u0 -= U_DIM(mrisp) ;
    if (u1 < 0)  /* enforce spherical topology  */
      u1 += U_DIM(mrisp) ;
    if (u1 >= U_DIM(mrisp))
      u1 -= U_DIM(mrisp) ;
    if (v0 < 0)  /* enforce spherical topology  */
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0)  /* enforce spherical topology  */
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;

    /* 0,0 */
#if 0
    d = sqrt(du*du+dv*dv) ;   /* distance to 0,0 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+dv*dv ;   /* distance to 0,0 */
    d = exp(-d/two_sigma_sq) ;
#endif
    total_d = distances[u0][v0] ;
    d /= total_d ;
    *IMAGEFpix(mrisp->Ip, u0,v0) += d*vertex->curv ;  /* weight by distance */
    if ((u0 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

    /* 1,0 */
    d = 1.0 - du ;          /* distance to u1 */
#if 0
    d = sqrt(d*d+dv*dv) ;   /* distance to u1,v0 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = d*d+dv*dv ;   /* distance to u1,v0 */
    d = exp(-d/two_sigma_sq) ;
#endif
    total_d = distances[u1][v0] ;
    d = d / total_d ;
    *IMAGEFpix(mrisp->Ip, u1, v0) += d*vertex->curv ;  /* weight by distance */
    if ((u1 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

    /* 0,1 */
    d = 1.0 - dv ; 
#if 0
    d = sqrt(du*du+d*d) ;   /* distance to u0,v1 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+d*d ;   /* distance to u0,v1 */
    d = exp(-d/two_sigma_sq) ;
#endif
    total_d = distances[u0][v1] ;
    d = d / total_d ;
    *IMAGEFpix(mrisp->Ip, u0, v1) += d*vertex->curv ;  /* weight by distance */
    if ((u0 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

    /* 1, 1 */
    du = 1.0 - du ; dv = 1.0 - dv ;
#if 0
    d = sqrt(du*du+dv*dv) ;   /* distance to 1,1 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = (du*du+dv*dv) ;   /* distance to 1,1 */
    d = exp(-d/two_sigma_sq) ;
#endif
    total_d = distances[u1][v1] ;
    d = d / total_d ;
    *IMAGEFpix(mrisp->Ip, u1, v1) += d*vertex->curv ;  /* weight by distance */
    if ((u1 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

  }

  if (DEBUG_U >= 0)
    fprintf(stderr,"curv[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, 
            *IMAGEFpix(mrisp->Ip, DEBUG_U, DEBUG_V));

#if 0
  /* normalize by total distances */
  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    for (v = 0 ; v <= V_MAX_INDEX ; v++)
    {
      if (u == DEBUG_U && v == DEBUG_V)
        DiagBreak() ;
      if (filled[u][v] >= 0)
      {
        *IMAGEFpix(mrisp->Ip, u,v) /= distances[u][v] ;
      }
    }
  }
#endif

  /* fill in values which were unmapped */
  do
  {
    unfilled = 0 ;
    for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
    {
      for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
      {
        if (filled[u][v] == UNFILLED_ELT)
        {
          for (total = 0.0f, n = 0, uk = -1 ; uk <= 1 ; uk++)
          {
            u1 = u + uk ;
            if (u1 < 0)  /* enforce spherical topology  */
              u1 += U_DIM(mrisp) ;
            else if (u1 >= U_DIM(mrisp))
              u1 -= U_DIM(mrisp) ;
            for (vk = -1 ; vk <= 1 ; vk++)
            {
              v1 = v + vk ;
              if (v1 < 0)  /* enforce spherical topology  */
                v1 += V_DIM(mrisp) ;
              else if (v1 >= V_DIM(mrisp))
                v1 -= V_DIM(mrisp) ;
              if (filled[u1][v1] != UNFILLED_ELT)
              {
                total += *IMAGEFpix(mrisp->Ip, u1, v1) ;
                n++ ;
              }
            }
          }
          if (n)
          {
            filled[u][v] = FILLED_ELT ;
            *IMAGEFpix(mrisp->Ip, u, v) = total / (float)n ;
          }
          else
            unfilled++ ;
        }
      }
    }
  } while (unfilled > 0) ;

  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    free(filled[u]) ;
    free(distances[u]) ;
  }
  free(filled) ; free(distances) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    ImageWrite(mrisp->Ip, "sphere.hipl") ;

  return(mrisp) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISfromParameterization(MRI_SP *mrisp, MRI_SURFACE *mris)
{
  float     a, b, c, phi, theta, x, y, z, uf, vf, du, dv, curv ;
  int       vno, u0, v0, u1, v1 ;
  VERTEX    *vertex ;

  if (!mris)
    mris = MRISclone(mrisp->mris) ;

  if (FZERO(mris->a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }
  else
  {
    a = mris->a ;
    b = mris->b ;
    c = mris->c ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(vertex->y/b, vertex->x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    phi = atan2(sqrt(c*c-z*z), z) ;
    if (phi < RADIANS(1))
      DiagBreak() ;
    if (vno == 60935)
      DiagBreak() ;

    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u0 = floor(uf) ; u1 = ceil(uf) ;
    v0 = floor(vf) ; v1 = ceil(vf) ;
    du = uf - (float)u0 ; dv = vf - (float)v0 ;

#if 0
    if (vno == 48092)
    {
      DiagBreak() ; 
      fprintf(stderr, "vno %d (%2.1f, %2.1f, %2.1f), u,v = (%2.1f, %2.1f)\n",
             vno, x, y, z, uf, vf) ;
    }
#endif

    /* enforce spherical topology  */
    if (u0 < 0)
      u0 += U_DIM(mrisp) ;
    if (u0 >= U_DIM(mrisp))
      u0 -= U_DIM(mrisp) ;
    if (u1 >= U_DIM(mrisp))
      u1 -= U_DIM(mrisp) ;
    if (v0 < 0) 
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0) 
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;
    
    
    /* do bilinear interpolation */
    curv = 
      du        * dv        * *IMAGEFpix(mrisp->Ip, u1, v1) +
      (1.0f-du) * dv        * *IMAGEFpix(mrisp->Ip, u0, v1) +
      (1.0f-du) * (1.0f-dv) * *IMAGEFpix(mrisp->Ip, u0, v0) +
      du        * (1.0f-dv) * *IMAGEFpix(mrisp->Ip, u1, v0) ;

    vertex->curv = curv ;
  }

  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPclone(MRI_SP *mrisp_src)
{
  MRI_SP   *mrisp_dst ;

  mrisp_dst = (MRI_SP *)calloc(1, sizeof(MRI_SP)) ;
  mrisp_dst->mris = mrisp_src->mris ;
  mrisp_dst->Ip = ImageClone(mrisp_src->Ip) ;

  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_LEN  4

MRI_SP *
MRISPblur(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma)
{
  int         u, v, cart_klen, klen, khalf, uk, vk, u1, v1, no_sphere ;
  double      k, a, b, c, total, ktotal, sigma_sq_inv, udiff, vdiff, sin_sq_u, 
              phi ;
  MRI_SURFACE *mris ;

  no_sphere = getenv("NO_SPHERE") != NULL ;
  if (no_sphere)
    fprintf(stderr, "disabling spherical geometry\n") ;

#if 0
  IMAGE       *Iblur ;

  Iblur = ImageGaussian1d(sigma, 0) ;
  if (!mrisp_dst)
    mrisp_dst = MRISPclone(mrisp_src) ;
  ImageCircularConvolveGaussian(mrisp_src->Ip, Iblur, mrisp_dst->Ip, 0) ;

#else
  if (!mrisp_dst)
    mrisp_dst = MRISPclone(mrisp_src) ;

  /* determine the size of the kernel */
  cart_klen = (int)nint(6.0f * sigma)+1 ;
  if (ISEVEN(cart_klen))   /* ensure it's odd */
    cart_klen++ ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "blurring surface, sigma = %2.3f, cartesian klen = %d\n",
            sigma, cart_klen) ;

  mris = mrisp_src->mris ;
  if (FZERO(mris->a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }
  else
  {
    a = mris->a ;
    b = mris->b ;
    c = mris->c ;
  }

  if (FZERO(sigma))
    sigma_sq_inv = BIG ;
  else
    sigma_sq_inv = 1.0f / (sigma*sigma) ;
  
  for (u = 0 ; u < U_DIM(mrisp_src) ; u++)
  {
    fprintf(stderr, "\r%3.3d of %d     ", u, U_DIM(mrisp_src)) ;
    phi = (double)u*PHI_MAX / PHI_DIM(mrisp_src) ;
    sin_sq_u = sin(phi) ; sin_sq_u *= sin_sq_u ;
    if (!FZERO(sin_sq_u))
    {
      k = cart_klen * cart_klen ;
      klen = sqrt(k + k/sin_sq_u) ;
      if (klen > MAX_LEN*cart_klen)
        klen = MAX_LEN*cart_klen ;
    }
    else
      klen = MAX_LEN*cart_klen ;  /* arbitrary max length */
    if (no_sphere)
      sin_sq_u = 1.0f, klen = cart_klen ;
    khalf = klen/2 ;
    for (v = 0 ; v < V_DIM(mrisp_src) ; v++)
    {
      /*      theta = (double)v*THETA_MAX / THETA_DIM(mrisp_src) ;*/
      if (u == DEBUG_U && v == DEBUG_V)
        DiagBreak() ;

      total = ktotal = 0.0 ;
      for (uk = -khalf ; uk <= khalf ; uk++)
      {
        udiff = (double)(uk*uk) ;  /* distance squared in u */

        u1 = u + uk ;
        if (u1 < 0)  /* enforce spherical topology  */
          u1 += U_DIM(mrisp_src) ;
        else if (u1 >= U_DIM(mrisp_src))
          u1 -= U_DIM(mrisp_src) ;

#if 0       
        phi = (double)u1*PHI_MAX / PHI_DIM(mrisp_src) ;
        sin_sq_u = sin(phi) ; sin_sq_u *= sin_sq_u ;
        if (no_sphere)
          sin_sq_u = 1.0f ;
#endif

        for (vk = -khalf ; vk <= khalf ; vk++)
        {
          vdiff = (double)(vk*vk) ;
          k = exp(-(udiff+sin_sq_u*vdiff)*sigma_sq_inv) ;
          v1 = v + vk ;
          if (v1 < 0)  /* enforce spherical topology */
            v1 += V_DIM(mrisp_src) ;
          else if (v1 >= V_DIM(mrisp_src))
            v1 -= V_DIM(mrisp_src) ;
          ktotal += k ;
          total += k**IMAGEFpix(mrisp_src->Ip, u1, v1) ;
        }
      }
      if (u == DEBUG_U && v == DEBUG_V)
        DiagBreak() ;
      total /= ktotal ;   /* normalize weights to 1 */
      *IMAGEFpix(mrisp_dst->Ip, u, v) = total ;
    }
  }
#endif

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;

  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPalloc(MRI_SURFACE *mris, float scale)
{
  MRI_SP   *mrisp ;
  int      u_dim, v_dim ;

  u_dim = nint(sqrt((float)mris->nvertices/(2*scale))) / 3 ;
  v_dim = 2*u_dim ;
  u_dim = 64 ; v_dim = 128 ; /* for fft */
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "allocating %d by %d parameterization\n",u_dim,v_dim) ;
  mrisp = (MRI_SP *)calloc(1, sizeof(MRI_SP)) ;
  if (!mrisp)
    ErrorExit(ERROR_NOMEMORY, "MRISPalloc(%d, %d): allocation failed",
              u_dim, v_dim) ;
  mrisp->mris = mris ;
  mrisp->Ip = ImageAlloc(v_dim, u_dim, PFFLOAT, 1) ;
  if (!mrisp->Ip)
    ErrorExit(ERROR_NOMEMORY, "MRISPalloc(%d, %d): allocation failed",
              u_dim, v_dim) ;

  return(mrisp) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISPfree(MRI_SP **pmrisp)
{
  MRI_SP  *mrisp ;

  mrisp = *pmrisp ;
  *pmrisp = NULL ;
  ImageFree(&mrisp->Ip) ;
  free(mrisp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPalign(MRI_SP *mrisp_orig, MRI_SP *mrisp_src, MRI_SP *mrisp_tmp, 
           MRI_SP *mrisp_dst)
{
  IMAGE  *Icorr = NULL ;
  int    u_peak, v_peak ;
  float  peak_val, delta_theta, delta_phi ;

  Icorr = ImageCorrelate(mrisp_src->Ip, mrisp_tmp->Ip, 0, NULL) ;
  ImageFindPeak(Icorr, &v_peak, &u_peak, &peak_val) ;
  u_peak -= U_DIM(mrisp_src)/2 ; v_peak -= V_DIM(mrisp_src)/2 ;
  delta_phi = u_peak * PHI_MAX / PHI_DIM(mrisp_src) ;
  delta_theta = v_peak * THETA_MAX / THETA_DIM(mrisp_src) ;
  fprintf(stderr, "peak found at (%d, %d) (%2.2f, %2.2f)\n", 
          u_peak, v_peak, DEGREES(delta_phi), DEGREES(delta_theta)) ;

  mrisp_dst = MRISPtranslate(mrisp_orig, NULL, u_peak, v_peak) ;
  MRISrotate(mrisp_orig->mris, mrisp_orig->mris, delta_phi, delta_theta) ;

#if 0 
  ImageWrite(mrisp_src->Ip, "Isrc.hipl") ;
  ImageWrite(mrisp_tmp->Ip, "Itemplate.hipl") ;
  ImageWrite(mrisp_dst->Ip, "Iorig.hipl") ;
#endif
  ImageWrite(Icorr, "Icorr.hipl") ;
  ImageWrite(mrisp_dst->Ip, "Ialign.hipl") ;

  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPtranslate(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, int du, int dv)
{
  int    u, v, udst, vdst ;

  if (!mrisp_dst)
    mrisp_dst = MRISPclone(mrisp_src) ;

  /* now form the destination as a translated copy of the source */
  for (u = 0 ; u <= U_MAX_INDEX(mrisp_src) ; u++)
  {
    udst = u + du ;
    if (udst < 0)  /* enforce spherical topology  */
      udst += U_DIM(mrisp_src) ;
    if (udst >= U_DIM(mrisp_src))
      udst -= U_DIM(mrisp_src) ;
    for (v = 0 ; v <= V_MAX_INDEX(mrisp_src) ; v++)
    {
      vdst = v + dv ;
      if (vdst < 0)  /* enforce spherical topology  */
        vdst += V_DIM(mrisp_src) ;
      if (vdst >= V_DIM(mrisp_src))
        vdst -= V_DIM(mrisp_src) ;
      *IMAGEFpix(mrisp_dst->Ip, udst, vdst) =
        *IMAGEFpix(mrisp_src->Ip, u, v) ;
    }
  }
  return(mrisp_dst) ;
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
MRISsampleStatVolume(MRI_SURFACE *mris, STAT_VOLUME *sv,int time_point,int tal)
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
    if (tal)
      MRIworldToTalairachVoxel(sv->mri_pvals[0], x, y, z, &xt, &yt, &zt) ;
    else
      MRIworldToVoxel(sv->mri_pvals[0], x, y, z, &xt, &yt, &zt) ;
    xv = nint(xt) ; yv = nint(yt) ; zv = nint(zt) ;
    if (xv >= 0 && xv < width && yv >= 0 && yv <= height && zv >= 0&&zv<=depth)
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

      fwrite(&f,1,sizeof(float),fp);
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
  int         nvertices, magic, version, ix, iy, iz, vno, n, nfaces ;
  FILE        *fp ;
  VERTEX      *vertex ;
  char        fname[100], path[100] ;

  FileNamePath(mris->fname, path) ;
  sprintf(fname, "%s/%s", path, sname) ;
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,
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
      fread1(&vertex->num,fp);
      for (n=0;n<vertex->num;n++)   /* skip over face data */
        fread3(&vertex->f[n],fp);
    } 
    else 
      vertex->num = 0;
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
  int         ix, iy, iz, k, i, j, npts;
  FILE        *fp ;
  char        fname[100], path[100] ;

  FileNamePath(mris->fname, path) ;
  sprintf(fname, "%s/%s", path, pname) ;
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,
                      "MRISreadPatch(%s): could not open file", fname));


  fread(&npts,1,sizeof(int),fp);
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading patch %s with %d vertices (%2.1f%% of total)\n",
            pname, npts, 100.0f*(float)npts/(float)mris->nvertices) ;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].ripflag = TRUE;
  for (j=0;j<npts;j++)
  {
    fread(&i,1,sizeof(int),fp);
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
  }
  fclose(fp);
  mrisRipFaces(mris);
  mrisComputeNormals(mris);
#if 0
  mrisComputeBoundaryNormals(mris);
  mrisSmoothBoundaryNormals(mris,10);
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
  int    k,i,npts;
  float  x,y,z;
  FILE   *fp;

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

  fwrite(&npts,1,sizeof(int),fp);
  for (k=0;k<mris->nvertices;k++)
  if (!mris->vertices[k].ripflag)
  {
    i = (mris->vertices[k].border)?-(k+1):k+1;
    fwrite(&i,1,sizeof(int),fp);
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
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisSmoothBoundaryNormals(MRI_SURFACE *mris, int niter)
{
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
  int k,m,n;
  vertex_type *v;
  float sumx,sumy,r,nx,ny,f;

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
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISflattenPatch(MRI_SURFACE *mris, char *dir)
{
  float       x,y,z,d,d1,d2;
  float       nx,ny,nz ;
  VERTEX      *v;
  int         k, an;
  FILE        *fp;
  char        fname[200];

  x = y = z = nx = ny = nz = 0; an = 0;

  /* calculate average normal and vertex position */
  for (k=0;k<mris->nvertices;k++)  
  if (!mris->vertices[k].ripflag)
  {
    v = &mris->vertices[k];
    x += v->x;
    y += v->y;
    z += v->z;
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
  sprintf(fname,"%s/surfer.mat",dir);
  fp = fopen(fname,"w");
  if (fp==NULL)
    ErrorPrintf(ERROR_NOFILE, "flatten: can't create file %s\n",fname);
  else 
  {
    fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",  nx*nz/d,  -nx,  ny/d,   0.0);
    fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",    d,       nz,   0.0,   0.0);
    fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",  -ny*nz/d,  ny,   nx/d,  0.0);
    fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",   0.0,      0.0,  0.0,   1.0);
    fclose(fp);
    if (Gdiag & DIAG_SHOW)
      printf("flatten: file %s written\n",fname);
  }

  transform(&nx,&ny,&nz,nx,ny,nz,d);
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "flatten: transformed n = {%2.1f, %2.1f, %2.1f}\n",
            nx,ny,nz);
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].z = 0;
  }
  mrisComputeNormals(mris);
  mrisOrientPlane(mris) ;
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
    if (vno == 4013)
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
        if (fabs(z/rsq) > 1.0)
          fprintf(stderr, "%d --> %d: curvaure = %2.1f\n", vno, vertex->v[i],
                  z/rsq) ;
      }
    }
    if (N > 0)
      vertex->curv = 0.5f * H / (float)N ;
    else
      vertex->curv = 0 ;
    if (vertex->curv > Hmax)
      Hmax = vertex->curv ;
    if (vertex->curv < Hmin)
      Hmin = vertex->curv ;
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
static int
mrisComputeCurvatureTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_n, *v_e1, *v_e2, *v_i, *v_delta, *v_tmp ;
  int     vno, i, N ;
  VERTEX  *vertex, *vnb ;
  float   rsq, z, u, v, l_curv ;

  l_curv = parms->l_curv ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_i = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;
  v_tmp = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z) ;
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z) ;
    V3_CLEAR(v_delta) ;
    N = 0 ;
    if (vno == 4013)
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
        N++ ;
        V3_SCALAR_MUL(v_n, -1.0f/rsq, v_tmp) ;
        V3_ADD(v_tmp, v_delta, v_delta) ;    /* -n/ri */
#if 1
        V3_SCALAR_MUL(v_e1, -2*u*z/(rsq*rsq), v_tmp) ;
        V3_ADD(v_tmp, v_delta, v_delta) ;    /* zi/ri^2 * (-2(vi.e1)e1) */
        V3_SCALAR_MUL(v_e2, -2*v*z/(rsq*rsq), v_tmp) ;
        V3_ADD(v_tmp, v_delta, v_delta) ;    /* zi/ri^2 * (-2(vi.e2)e2) */
#endif
      }
    }
    if (N > 0)
    {
      V3_SCALAR_MUL(v_delta, 1.0/N, v_delta) ;
/*
  note that the curvature is positive for convex regions (for compatibility
  with old A&M code, so we must invert the sign here.
*/
      vertex->dx += l_curv * -vertex->curv * V3_X(v_delta) ; 
      vertex->dy += l_curv * -vertex->curv * V3_Y(v_delta) ; 
      vertex->dz += l_curv * -vertex->curv * V3_Z(v_delta) ;
    }
  }
  
  VectorFree(&v_tmp) ;
  VectorFree(&v_delta) ;
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
mrisUpdateSurface(MRI_SURFACE *mris)
{
  mrisComputeNormals(mris);
  MRIScomputeMeanCurvature(mris) ;
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

