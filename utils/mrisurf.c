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


/*------------------------ STATIC PROTOTYPES -------------------------*/

static int   mrisComputeNormals(MRI_SURFACE *mris) ;
static int   mrisFindNeighbors(MRI_SURFACE *mris) ;
static void  mrisNormalize(float v[3]) ;
static float mrisTriangleArea(MRIS *mris, int fac, int n) ;
static int   mrisNormalFace(MRIS *mris, int fac,int n,float norm[]) ;
static int   mrisFindPoles(MRIS *mris) ;
static int   mrisReadTransform(MRIS *mris, char *mris_fname) ;
static int   mrisReadBinaryCurvature(MRI_SURFACE *mris, char *mris_fname) ;
static int   mrisReadBinaryAreas(MRI_SURFACE *mris, char *mris_fname) ;
static int   mrisReadTriangleProperties(MRI_SURFACE *mris, char *mris_fname) ;
/*static int   mrisReadFieldsign(MRI_SURFACE *mris, char *fname) ;*/
float        mrisComputeAreaRMS(MRI_SURFACE *mris) ;
static int   mrisScaleEllipsoidArea(MRI_SURFACE *mris) ;
static int   mrisOrientEllipsoidArea(MRI_SURFACE *mris) ;
static int   mrisCountNegativeVertices(MRI_SURFACE *mris) ;
static int   mrisAverageGradients(MRI_SURFACE *mris, int num_avgs) ;
static int   mrisCalculateTriangleProperties(MRI_SURFACE *mris) ;

/*--------------------------------------------------------------------*/

/*--------------------- CONSTANTS AND MACROS -------------------------*/

#define NEW_VERSION_MAGIC_NUMBER  16777215
#define START_Y                   (-128)
#define SLICE_THICKNESS           1

#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y,\
                                               v1->z-v0->z)

#if 0

/* 77557 is horizontal with positive area */
/* 126906 has dy = 0 with negative area */
#define DEBUG_FACE(fno)   (((fno) == 77540))
#define DEBUG_FACE(fno)   (((fno) == 79379) && 0)
#define DEBUG_FACE(fno)   (((fno) == 77557) && (Gdiag & DIAG_SURFACE))
#else
#define DEBUG_FACE(fno)   (((fno) == 75977) && (Gdiag & DIAG_SURFACE) && 0)
#define VDEBUG_FACE(fno)   (DEBUG_FACE(fno) && 0)
#endif
#define DEBUG_VERTEX(v)   (((v) == 40501) && (Gdiag & DIAG_SURFACE) && 1)

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
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "new surface file format\n");
  }
  else 
  {
    rewind(fp);
    version = 0;
    if (Gdiag & DIAG_SHOW)
      printf("surfer: old surface file format\n");
  }
  fread3(&nvertices, fp);
  fread3(&nfaces, fp);

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading %d vertices and %d faces.\n",nvertices,nfaces);

  mris = MRISalloc(nvertices, nfaces) ;

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
      vertex->tri_angle = (float *)calloc(vertex->num,sizeof(float));
      if (!vertex->tri_angle)
        ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                  vertex->num) ;
      vertex->orig_tri_angle = (float *)calloc(vertex->num,sizeof(float));
      if (!vertex->orig_tri_angle)
        ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                  vertex->num) ;
      vertex->tri_area = (float *)calloc(vertex->num,sizeof(float));
      if (!vertex->tri_area)
        ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                  vertex->num) ;
      vertex->orig_tri_area = (float *)calloc(vertex->num,sizeof(float));
      if (!vertex->orig_tri_area)
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
      vertex->tri_angle = (float *)calloc(vertex->num,sizeof(float));
      if (!vertex->tri_angle)
        ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                  vertex->num) ;
      vertex->orig_tri_angle = (float *)calloc(vertex->num,sizeof(float));
      if (!vertex->orig_tri_angle)
        ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                  vertex->num) ;
      vertex->tri_area = (float *)calloc(vertex->num,sizeof(float));
      if (!vertex->tri_area)
        ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                  vertex->num) ;
      vertex->orig_tri_area = (float *)calloc(vertex->num,sizeof(float));
      if (!vertex->orig_tri_area)
        ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces",
                  vertex->num) ;
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
  mrisReadBinaryCurvature(mris, fname) ;
  mrisReadBinaryAreas(mris, fname) ;
  if (mrisReadTriangleProperties(mris, fname))
    mrisCalculateTriangleProperties(mris) ;
  mrisFindPoles(mris) ;
#if 0
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
    if (mris->vertices[vno].tri_area)
      free(mris->vertices[vno].tri_area) ;
    if (mris->vertices[vno].orig_tri_area)
      free(mris->vertices[vno].orig_tri_area) ;
    if (mris->vertices[vno].tri_angle)
      free(mris->vertices[vno].tri_angle) ;
    if (mris->vertices[vno].orig_tri_angle)
      free(mris->vertices[vno].orig_tri_angle) ;
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

  if (Gdiag & DIAG_SHOW)
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
  if (Gdiag & DIAG_SHOW)
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

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "(%2.0f, %2.0f, %2.0f) --> (%2.0f, %2.0f, %2.0f), ctr "
            "(%2.0f, %2.0f, %2.0f)\n",
            mris->xlo, mris->ylo, mris->zlo, mris->xhi, mris->yhi, mris->zhi,
            mris->xctr, mris->yctr, mris->zctr);
  if (Gdiag & DIAG_SHOW)
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

  if (Gdiag & DIAG_SHOW)
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
static int
mrisReadBinaryCurvature(MRI_SURFACE *mris, char *mris_fname)
{
  int    k,i,vnum,fnum;
  float  curv, curvmin, curvmax;
  char   fname[100], fpref[100], hemi[20], *cp ;
  FILE   *fp;
  
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading curvature file...") ;

  cp = strrchr(mris_fname, '.') ;
  if (!cp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "mrisReadBinaryCurvature(%s): could not scan "
                 "hemisphere from fname", mris_fname)) ;
  strncpy(hemi, cp-2, 2) ; hemi[2] = 0 ;
  FileNamePath(mris_fname, fpref) ;
  sprintf(fname, "%s/%s.curv", fpref, hemi) ;
  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "mrisReadBinaryCurvature: could not open %s", 
                 fname)) ;

  fread3(&vnum,fp);
  fread3(&fnum,fp);
  if (vnum!= mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "mrisReadBinaryCurvature: incompatible vertex "
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
  if (Gdiag & DIAG_SHOW)
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

  if (Gdiag & DIAG_SHOW)
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
#define DEFAULT_A  44.0f
#define DEFAULT_B  122.0f
#define DEFAULT_C  70.0f

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

  /*  printf("ellipsoid_project(%f,%f,%f)\n",a,b,c);*/
  for (k=0;k<mris_dst->nvertices;k++) 
  {
    v = &mris_dst->vertices[k];
/*
    printf("%6d: before: %6.2f\n",k,SQR(v->x/a)+SQR(v->y/b)+SQR(v->z/c));
*/
    x = v->x;
    y = v->y;
    z = v->z;
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    a2 = a*a;
    b2 = b*b;
    c2 = c*c;
    a4 = a2*a2;
    b4 = b2*b2;
    c4 = c2*c2;
    a6 = a2*a4;
    b6 = b2*b4;
    c6 = c2*c4;
    f = x2/a6+y2/b6+z2/c6;
    g = 2*(x2/a4+y2/b4+z2/c4);
    h = x2/a2+y2/b2+z2/c2-1;
    d = (-g+(float)sqrt((double)(g*g-4*f*h)))/(2*f);
    dx = d*x/a2;
    dy = d*y/b2;
    dz = d*z/c2;
    v->x = x+dx ;
    v->y = y+dy;
    v->z = z+dz;
    dist = (float)sqrt((double)(dx*dx+dy*dy+dz*dz));
    avgdist += dist;
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


#if 0
/*
  this one projects along the line from the origin to the ellipsoidal
  surface - not orthographic unless the ellipsoid is a sphere.
  */
MRI_SURFACE *
MRISprojectOntoEllipsoid(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, 
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
#endif
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
        vdst->tri_area[n] = vsrc->tri_area[n] ;
        vdst->orig_tri_area[n] = vsrc->orig_tri_area[n] ;
        vdst->tri_angle[n] = vsrc->tri_angle[n] ;
        vdst->orig_tri_angle[n] = vsrc->orig_tri_angle[n] ;
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
  float       x0, y0, z0 ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  x0 = mris_src->xctr ; y0 = mris_src->yctr ; z0 = mris_src->zctr ;
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    vdst->x -= x0 ; vdst->y -= y0 ; vdst->z -= z0 ;
  }

  mris_src->xctr = mris_src->yctr = mris_src->zctr = 0 ;
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
#define MAX_NEGATIVE        50000
#define DELTA_T             0.25f /* 0.05f */  /* 0.3 */
#define MOMENTUM            0.9f
#define ONE_MINUS_MOMENTUM  (1.0f - MOMENTUM)
#define DT_DECREASE         1.0f /* 0.7f */ /* 0.95f*/
#define DT_INCREASE         1.0f /* 1.04f */ /* 1.1f*/
#define ERROR_RATIO         1.00f
#define MIN_DT              0.000001f
#define MAX_MIN_DTS         30
#define L_SPRING            0.1f   /* coefficient of spring term */

#define WRITE_ITERATIONS    500

#define MAX_VERTICES  500000
#define N_AVERAGES    400

MRI_SURFACE *
MRISunfold(MRI_SURFACE *mris, int niterations, float base_momentum, 
           float l_area, float l_angle,  float l_corr)
{
  char    base_name[100] ;
  int     vno, t, fno, imin, min_dts, write_iterations = WRITE_ITERATIONS,
          negative, index, index_a, index_b, no_projection,
          n_averages, base_averages ;
  VERTEX  *v, *va, *vb ;
  VECTOR  *v_a, *v_b, *v_n, *v_delta, *v_b_x_n, *v_n_x_a, *v_sum, *v_spring ;
  FACE    *face ;
  float   rms, old_rms, delta_t,momentum, one_minus_momentum, min_dz, coef,dz,
          triangle_area, area, delta, dt_increase, dt_decrease, dt_base, 
          l_spring, avg_slope, angle, triangle_angle, a_dot_b, len ;

  dt_increase = DT_INCREASE ; dt_decrease = DT_DECREASE ;

  sprintf(base_name, "surf/%s.tmp", mris->hemisphere == LEFT_HEMISPHERE ?
              "lh" : "rh") ;
  dt_base = delta_t = DELTA_T ;
  base_averages = N_AVERAGES ;
  l_spring = L_SPRING ;
  {
    char *cp ;

    cp = getenv("N_AVERAGES") ;
    if (cp)
      sscanf(cp, "%d", &base_averages) ;

    cp = getenv("BASE_NAME") ;
    if (cp)
      sprintf(base_name, "surf/%s.%s", mris->hemisphere == LEFT_HEMISPHERE ?
              "lh" : "rh", cp) ;
    cp = getenv("L_SPRING") ;
    if (cp)
    {
      sscanf(cp, "%f", &l_spring) ;
      fprintf(stderr, "using l_spring = %2.2f\n", l_spring) ;
    }
    cp = getenv("L_AREA") ;
    if (cp)
    {
      sscanf(cp, "%f", &l_area) ;
      fprintf(stderr, "using l_area = %2.3f\n", l_area) ;
    }
    cp = getenv("DT") ;
    if (cp)
    {
      sscanf(cp, "%f", &dt_base) ;
      delta_t = dt_base ;
    }

    cp = getenv("DT_INCREASE") ;
    if (cp)
      sscanf(cp, "%f", &dt_increase) ;
    cp = getenv("DT_DECREASE") ;
    if (cp)
      sscanf(cp, "%f", &dt_decrease) ;
    cp = getenv("WRITE_ITERATIONS") ;
    if (cp)
    {
      sscanf(cp, "%d", &write_iterations) ;
      fprintf(stderr, 
          "writing out intermediate files every %d iterations to %s%4.4d\n",
          write_iterations, base_name, 0) ;
    }

    cp = getenv("NO_PROJ") ;
    no_projection = (cp != NULL) ;
    if (no_projection)
      fprintf(stderr, "no ellipsoid projection in unfolding\n") ;

  }

  if (base_momentum < 0.0f)
    base_momentum = MOMENTUM ; 
  momentum = base_momentum ;
  one_minus_momentum = 1.0f - momentum ;

  printf(
         "dt_inc=%2.3f, dt_dec=%2.3f, mom=%2.2f, nav=%d, l_sp=%2.2f"
         "l_a=%2.2f\n"
         , dt_increase, dt_decrease, base_momentum, base_averages, l_spring,
         l_area) ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_delta = VectorAlloc(3, MATRIX_REAL) ;   
  v_sum = VectorAlloc(3, MATRIX_REAL) ;   
  v_n_x_a = VectorAlloc(3, MATRIX_REAL) ;      
  v_b_x_n = VectorAlloc(3, MATRIX_REAL) ;      
  v_spring = VectorAlloc(3, MATRIX_REAL) ;

  MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
  rms = old_rms = mrisComputeAreaRMS(mris) ;
  min_dts = 0 ;
  imin = -1 ;
  avg_slope = (float)base_averages / (float)niterations ;
  for (t = 0 ; t < niterations ; t++)
  {
    n_averages = base_averages - nint((float)t*avg_slope) ;
    negative = mrisCountNegativeVertices(mris) ;
    if (!(t%write_iterations))
    {
      char fname[100] ;
      sprintf(fname, "%s%4.4d", base_name, t) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "writing %s...", fname) ;
      MRISwrite(mris, fname) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "done.\n") ;
    }

    delta_t = 4.0f * dt_base * sqrt((float)n_averages+1.0f) ;
    fprintf(stdout,"%d of %d: dt = %2.3f, rms: %2.6f, neg: %d (%2.3f), "
            "avgs = %d\n",
            t, niterations, delta_t, 100.0f*rms, negative,mris->neg_area,
            n_averages) ;
    fflush(stdout) ;
    min_dz = 10000.0f ;


    /* calculcate movement of each vertex */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;

if (DEBUG_VERTEX(vno))
{
  fprintf(stderr, 
          "\n\nt = %d, vertex %d at (%2.3f, %2.3f, %2.3f), n = (%2.3f, %2.3f"
          ", %2.3f)\n area = %2.3f, oarea = %2.3f\n", 
          t, vno, v->x,v->y,v->z,v->nx,v->ny,v->nz,v->area, v->origarea) ;
}

      VectorClear(v_delta) ;
      VectorClear(v_spring) ;
      dz = 0.0 ;
      VECTOR3_LOAD(v_n, v->nx, v->ny, v->nz) ;
      for (fno = 0 ; fno < v->num ; fno++)  /* for each neighboring face */
      {
        triangle_area = v->orig_tri_area[fno] ;    /* original area */
        triangle_angle = v->orig_tri_angle[fno] ;  /* original angle */

        face = &mris->faces[v->f[fno]] ;
        index = v->n[fno] ;
        index_a = index == 0 ? VERTICES_PER_FACE-1 : index-1 ;
        index_b = index == VERTICES_PER_FACE-1 ? 0 : index+1 ;
        va = &mris->vertices[face->v[index_a]] ;
        vb = &mris->vertices[face->v[index_b]] ;
        VERTEX_EDGE(v_a, va, v) ;  /* for some reason A&M invert this */
        VERTEX_EDGE(v_b, v, vb) ;
        dz += (VECTOR_ELT(v_a,3)*VECTOR_ELT(v_a,3) +
               VECTOR_ELT(v_b,3)*VECTOR_ELT(v_b,3)) ;

        area = VectorTripleProduct(v_a, v_b, v_n)*.5f ;
        a_dot_b = VectorDot(v_a, v_b) ;
        angle = atan2(area*2.0f, a_dot_b) ;

        /* calculate movement based on area term */
        delta = area - triangle_area ;
        VectorCrossProduct(v_n, v_a, v_n_x_a) ;   /* n x a */
        VectorCrossProduct(v_n, v_b, v_b_x_n) ;   /* b x n * -1 */
        VectorAdd(v_n_x_a, v_b_x_n, v_sum) ;
        VectorScalarMul(v_sum, delta, v_sum) ;
        VectorAdd(v_sum, v_delta, v_delta) ;

        /* calculate movement based on angle term */
        delta = angle - triangle_angle ;
        len = VectorLen(v_a) ; VectorScalarMul(v_n_x_a, len*len, v_n_x_a) ;
        len = VectorLen(v_b) ; VectorScalarMul(v_b_x_n, len*len, v_b_x_n) ;
        VectorAdd(v_n_x_a, v_b_x_n, v_sum) ;
        VectorScalarMul(v_sum, delta, v_sum) ;
        VectorAdd(v_sum, v_spring, v_spring) ;
      }  /* done with all faces of this vertex */

      dz /= (float)v->num ;
      if (dz < min_dz)
      {
        min_dz = dz ;
        imin = vno ;
      }

#if 0
      /*  calculate spring force - go through all neighbors */
      VectorClear(v_spring) ;
      for (vnb = 0 ; vnb < v->vnum ; vnb++)
      {
        va = &mris->vertices[v->v[vnb]] ;    /* neighboring vertex pointer */
        VERTEX_EDGE(v_a, v, va) ;            /* edge connecting them */
        VectorAdd(v_a, v_spring, v_spring) ;
      }
#endif

      /* now calculate movement of vertex */
      coef = one_minus_momentum * l_spring * delta_t ;
      VectorScalarMul(v_spring, coef, v_spring) ;  /* spring step */


      coef = one_minus_momentum * l_area * delta_t ;
      VectorScalarMul(v_delta, coef, v_delta) ;    /* area step */
      VectorAdd(v_spring, v_delta, v_delta) ;      /* combination */

      if (!finite(VECTOR_ELT(v_delta,1)) || !finite(VECTOR_ELT(v_delta,2)) || 
          !finite(VECTOR_ELT(v_delta,3)))
      { 
        static int first = 1;
      
        if (first)
          fprintf(stderr, "sse is not finite at vertex %d!\n", vno) ;
        first = 0 ;
        DiagBreak() ;
        exit(1) ;
      }

      v->dx = VECTOR_ELT(v_delta,1) + momentum * v->dx ;
      v->dy = VECTOR_ELT(v_delta,2) + momentum * v->dy ;
      v->dz = VECTOR_ELT(v_delta,3) + momentum * v->dz ;
#define SYNCHRONOUS_UPDATE  1

#if !SYNCHRONOUS_UPDATE
      v->x += v->dx ; v->y += v->dy ; v->z += v->dz ;
#endif

      if (!finite(v->x) || !finite(v->y) || !finite(v->z))
      { 
        static int first = 1;
      
        if (first)
          fprintf(stderr, "sse is not finite at vertex %d!\n", vno) ;
        first = 0 ;
        DiagBreak() ;
        exit(1) ;
      }

if (DEBUG_VERTEX(vno))
{
  fprintf(stderr, "moving vertex by (%2.4f, %2.4f, %2.4f)\n",
          v->dx, v->dy, v->dz) ;
}
    }    /* done with all vertices */

    if (!t && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
      fprintf(stderr, "mindz = %2.3f at %d\n",min_dz, imin) ;

#if SYNCHRONOUS_UPDATE
    mrisAverageGradients(mris, n_averages) ;
    /* now apply the deltas */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->x += v->dx ; v->y += v->dy ; v->z += v->dz ;
    }
#endif


    MRISupdateEllipsoidSurface(mris) ;

    rms = mrisComputeAreaRMS(mris) ;
    if (rms > old_rms)   /* error has increased, reduce the time step */
      delta_t *= dt_decrease ;
    else
      delta_t *= dt_increase ;  /* error decreased - increase time step */

#if 0
    if (rms > old_rms * ERROR_RATIO)  /* big increase - kill momentum */
    {
      momentum = 0.0f ; one_minus_momentum = 1.0f ;
      /* restore old vertex positions */
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        v->x -= v->dx ; v->y -= v->dy ; v->z -= v->dz ;
        v->dx = v->dy = v->dz = 0.0f ;
      }
      MRISupdateEllipsoidSurface(mris) ;
      rms = mrisComputeAreaRMS(mris) ;
    }
    else   /* accept the change */
#endif
    {
      if (!no_projection)
        MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
      momentum = base_momentum ; one_minus_momentum = 1.0f - momentum ;
    }
    old_rms = rms ;
    if (delta_t <= MIN_DT)
    {
      delta_t = 2*dt_base ;
      if (min_dts++ >= MAX_MIN_DTS)
        break ;
    }
  }

  MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_delta) ;
  VectorFree(&v_sum) ;
  VectorFree(&v_n) ;
  VectorFree(&v_b_x_n) ;
  VectorFree(&v_n_x_a) ;
  VectorFree(&v_spring) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;

  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float
mrisComputeAreaRMS(MRI_SURFACE *mris)
{
  float   sse, delta, n ;
  int     vno, fno ;
  VERTEX  *v ;

  sse = 0.0f ;
  for (n = 0.0f, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    for (fno = 0 ; fno < v->num ; fno++)
    {
      delta = v->tri_area[fno] - v->orig_tri_area[fno] ;
      sse += delta*delta ;
      if (!finite(sse))
      {
        static int first = 1;
        
        if (first)
          fprintf(stderr, "sse is not finite at vertex %d!\n", vno) ;
        first = 0 ;
        DiagBreak() ;
        exit(1) ;
      }
    }
  }


  return(sqrt(sse) / (float)mris->nvertices) ;
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
  int     vno, fno ;
  VERTEX  *v ;

  area_scale = mris->total_area / mris->orig_area ;
  if (mris->initialized < 2)  /* only do it once */
  {
    mris->initialized = 2 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      
      /* scale the area by the ratio of the ellipsoid area to that of the
         original surface.
         */
      v->origarea *= area_scale ;
      for (fno = 0 ; fno < v->num ; fno++)
        v->orig_tri_area[fno] *= area_scale ;
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
mrisOrientEllipsoidArea(MRI_SURFACE *mris)
{
  int     vno, negative ;
  VERTEX  *v ;
  float   dot ;

  negative = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
      
    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
       */
    dot = v->x * v->nx + v->y * v->ny + v->z * v->nz;
    if (dot < 0.0f)   /* not in same direction, area < 0 and reverse n */
    {
      negative++ ;
      v->area *= -1.0f ;
      v->nx *= -1.0f ; v->ny *= -1.0f ; v->nz *= -1.0f ;
    }
    else
      if (v->area < 0.0f)
        negative++ ;

  }
  return(negative) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the oriented area of each vertex
------------------------------------------------------*/
static int
mrisComputeEllipsoidArea(MRI_SURFACE *mris)
{
  int     vno, fno, index, index_a, index_b ;
  VERTEX  *v, *va, *vb ;
  float   dot, area ;
  VECTOR  *v_a, *v_b, *v_n, *v_n_avg ;
  FACE    *face ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_n_avg = VectorAlloc(3, MATRIX_REAL) ;      

  mris->neg_area = mris->total_area = 0.0f ;
  mris->zeros = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->area = 0.0f ;
    VectorClear(v_n_avg) ;/* will be normal of vertex - avg of all triangles */
if (DEBUG_VERTEX(vno))
{
  fprintf(stderr, "\ncomputing A and N for v %d at "
          "(%2.3f, %2.3f, %2.3f), oarea = %2.2f\n", 
          vno, v->x, v->y, v->z, v->origarea) ;
}
    for (fno = 0 ; fno < v->num ; fno++)  /* all faces of this vertex */
    {
      face = &mris->faces[v->f[fno]] ;
      index = v->n[fno] ;
      index_a = index == 0 ? VERTICES_PER_FACE-1 : index-1 ;
      index_b = index == VERTICES_PER_FACE-1 ? 0 : index+1 ;
      va = &mris->vertices[face->v[index_a]] ;
      vb = &mris->vertices[face->v[index_b]] ;
      VERTEX_EDGE(v_a, va, v) ;
      VERTEX_EDGE(v_b, v, vb) ;
      VectorCrossProduct(v_a, v_b, v_n) ;     /* compute the normal */
      VectorNormalize(v_n, v_n) ;             /* make it a unit vector */

      /* now give the area an orientation: if the unit normal is pointing
         inwards on the ellipsoid then the area should be negative.
         */
      dot = 
        v->x * VECTOR_ELT(v_n,1) + 
        v->y * VECTOR_ELT(v_n,2) + 
        v->z * VECTOR_ELT(v_n,3) ;
      if (dot < 0.0f)     /* pointing inwards - area is negative */
      {
        VectorScalarMul(v_n, -1.0f, v_n) ;
        v->tri_area[fno] *= -1.0f ;   /* area is negative */
        v->tri_angle[fno] *= -1.0f ;  /* angle is negative */
      }
      area = VectorTripleProduct(v_a, v_b, v_n)*.5f ;
      if (FZERO(area))
        mris->zeros++ ;
      v->area += area ;
      VectorAdd(v_n, v_n_avg, v_n_avg) ;  /* keep track of avg normal */

if (DEBUG_VERTEX(vno))
{
  fprintf(stderr, "face %d, area = %2.3f\n", fno, area) ;

  fprintf(stderr, "P1 %d at (%2.3f, %2.3f, %2.3f)\n", 
          va-mris->vertices, va->x, va->y, va->z) ;
  fprintf(stderr, "P2 at (%2.3f, %2.3f, %2.3f)\n", vb->x, vb->y, vb->z) ;
  fprintf(stderr, "a: ") ;
  MatrixPrintTranspose(stderr, v_a) ;
  fprintf(stderr, "b: ") ;
  MatrixPrintTranspose(stderr, v_b) ;
  fprintf(stderr, "n: ") ;
  MatrixPrintTranspose(stderr, v_n) ;
  DiagBreak() ;
}

    }
    if (v->area > 0.0f)
     mris->total_area += v->area ;
    else
     mris->neg_area += -v->area ;
    VectorNormalize(v_n_avg, v_n_avg) ;
    if (!FZERO(VectorLen(v_n_avg)))
    {
      v->nx = VECTOR_ELT(v_n_avg,1) ;
      v->ny = VECTOR_ELT(v_n_avg,2) ;
      v->nz = VECTOR_ELT(v_n_avg,3) ;
    }
    if (!finite(v->nx) || !finite(v->ny) || !finite(v->nz) || !finite(v->area))
    {
      static int first = 1;
      
      if (first)
        fprintf(stderr, "vertex %d parms not finite!\n", vno) ;
      first = 0 ;
      DiagBreak() ;
      exit(1) ;
      
    }
if (DEBUG_VERTEX(vno))
{
  fprintf(stderr, "area = %2.3f, normal = (%2.3f, %2.3f, %2.3f)\n", 
          v->area, v->nx, v->ny, v->nz) ;

}
  }

  mris->total_area /= 2.0f ;  /* every triangle counted twice */
#if 0
  if (Gdiag & DIAG_SHOW)
     fprintf(stderr, "orig = %2.0f, ELL = %2.0f, total area = %2.0f\n",
             mris->orig_area, TOTAL_AREA, mris->total_area) ;
#endif
  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_n) ;
  VectorFree(&v_n_avg) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mrisCountNegativeVertices(MRI_SURFACE *mris)
{
  int     vno, negative ;
  VERTEX  *v ;

  negative = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->area < 0.0f)
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
  mrisComputeNormals(mris);
  mrisOrientEllipsoidArea(mris) ;
#else
  mrisCalculateTriangleProperties(mris) ;
  mrisComputeEllipsoidArea(mris) ;
#endif
  mrisScaleEllipsoidArea(mris) ;
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
  int    i, vno, vnb ;
  float  dx, dy, dz, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < num_avgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      dx = v->dx ; dy = v->dy ; dz = v->dz ;
      for (vnb = 0 ; vnb < v->vnum ; vnb++)
      {
        vn = &mris->vertices[v->v[vnb]] ;    /* neighboring vertex pointer */
        dx += vn->dx ; dy += vn->dy ; dz += vn->dz ;
      }
      num = (float)v->vnum + 1.0f ;
      v->dx = dx / num ;
      v->dy = dy / num ;
      v->dz = dz / num ;
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
mrisReadTriangleProperties(MRI_SURFACE *mris, char *mris_fname)
{
  int     k,vnum,fnum, fno ;
  VERTEX  *v ;
  float   f;
  FILE    *fp;
  char    fname[100], fpref[100], hemi[20], *cp ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading triangle files...") ;

  cp = strrchr(mris_fname, '.') ;
  if (!cp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "mrisReadBinaryAreas(%s): could not scan "
                 "hemisphere from fname", mris_fname)) ;
  strncpy(hemi, cp-2, 2) ; hemi[2] = 0 ;
  FileNamePath(mris_fname, fpref) ;

  sprintf(fname, "%s/%s.triangle_area", fpref, hemi) ;

  fp = fopen(fname,"r");
  if (fp==NULL)   
    return(1) ;  /* doesn't exist */

  fread4((float *)&vnum,fp);
  fread4((float *)&fnum,fp);
  if (vnum!=mris->nvertices)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE, "mrisReadTriangleProperties: incompatible vertex "
                 "number in file %s", fname)) ;
  }
  
  mris->orig_area = 0.0f ;
  for (k=0;k<vnum;k++)
  {
    v = &mris->vertices[k] ;
    for (fno = 0 ; fno < v->num ; fno++)
    {
      f = freadFloat(fp);
      v->orig_tri_area[fno] = f ;
      mris->orig_area += f;
    }
  }

  /* hack to correct for overestimation of area in compute_normals */
  mris->orig_area /= 2; 
  if (Gdiag & DIAG_SHOW)
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
              (ERROR_NOFILE, "mrisReadTriangleProperties: incompatible vertex "
                 "number in file %s", fname)) ;
  }
  
  for (k=0;k<vnum;k++)
  {
    v = &mris->vertices[k] ;
    for (fno = 0 ; fno < v->num ; fno++)
    {
      f = freadFloat(fp);
      v->orig_tri_angle[fno] = f ;
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
  int     k, fno ;
  VERTEX  *v ;
  FILE    *fp;
  char    fname[100], fpref[100], hemi[20], *cp ;

  mrisCalculateTriangleProperties(mris) ;

  cp = strrchr(mris_fname, '.') ;
  if (!cp)
    ErrorReturn(ERROR_BADPARM, 
              (ERROR_BADPARM, "mrisReadTriangleProperties(%s): could not scan "
                 "hemisphere from fname", mris_fname)) ;
  strncpy(hemi, cp-2, 2) ; hemi[2] = 0 ;
  FileNamePath(mris_fname, fpref) ;


  sprintf(fname, "%s/%s.triangle_area", fpref, hemi) ;
  fp = fopen(fname,"wb");
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, 
                              "MRISwriteTriangleProperties: could not open %s",
                                fname)) ;

  /* calculate the area of all the triangles */
  fwrite4(mris->nvertices,fp);
  fwrite4(mris->nfaces,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k] ;
    for (fno = 0 ; fno < v->num ; fno++)
      putf(v->tri_area[fno], fp) ;
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
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k] ;
    for (fno = 0 ; fno < v->num ; fno++)
      putf(v->tri_angle[fno], fp) ;
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
mrisCalculateTriangleProperties(MRI_SURFACE *mris)
{
  VECTOR  *v_a, *v_b, *v_n, *v_b_x_n, *v_n_x_a ;
  VERTEX  *v, *va, *vb ;
  FACE    *face ;
  int     vno, fno, index, index_a, index_b ;
  float   area, angle, a_dot_b ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_n_x_a = VectorAlloc(3, MATRIX_REAL) ;      
  v_b_x_n = VectorAlloc(3, MATRIX_REAL) ;      

  mris->total_area = 0.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    for (fno = 0 ; fno < v->num ; fno++)  /* for each neighboring face */
    {
      face = &mris->faces[v->f[fno]] ;
      index = v->n[fno] ;
      index_a = index == 0 ? VERTICES_PER_FACE-1 : index-1 ;
      index_b = index == VERTICES_PER_FACE-1 ? 0 : index+1 ;
      va = &mris->vertices[face->v[index_a]] ;
      vb = &mris->vertices[face->v[index_b]] ;
      VERTEX_EDGE(v_a, va, v) ;  /* for some reason A&M invert this */
      VERTEX_EDGE(v_b, v, vb) ;
      VectorCrossProduct(v_a, v_b, v_n) ;
      area = VectorLen(v_n) * 0.5f ;
      a_dot_b = VectorDot(v_a, v_b) ;
      angle = atan2(area*2.0f, a_dot_b) ;
      v->tri_angle[fno] = angle ;
      v->tri_area[fno] = area ;
      mris->total_area += area ;
    }
  }

  mris->total_area /= 2.0f ;  /* because each triangle is counted twice */
  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_n) ;
  VectorFree(&v_n_x_a) ;
  VectorFree(&v_b_x_n) ;
  return(NO_ERROR) ;
}



