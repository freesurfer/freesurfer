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
/*static int   mrisReadFieldsign(MRI_SURFACE *mris, char *fname) ;*/
float        mrisComputeAreaRMS(MRI_SURFACE *mris) ;
static int   mrisScaleEllipsoidArea(MRI_SURFACE *mris) ;
static int   mrisOrientEllipsoidArea(MRI_SURFACE *mris) ;

/*--------------------------------------------------------------------*/

/*--------------------- CONSTANTS AND MACROS -------------------------*/

#define NEW_VERSION_MAGIC_NUMBER  16777215
#define START_Y                   (-128)
#define SLICE_THICKNESS           1


#if 0
#define DEBUG_FACE(fno)   (((fno) == 77540))
#else
#define DEBUG_FACE(fno)   (((fno) == 79379) && 0)
#endif
#define DEBUG_VERTEX(v,f)   (((f) == 77557) && 0)

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
      for (n = 0 ; n < 4 ; n++)
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
      for (m=0;m<4;m++)
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
  mrisFindPoles(mris) ;
  MRIScomputeFaceAreas(mris) ;

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
    for (n = 0 ; n < 4 ; n++)
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
    for (m=0;m<4;m++)
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
  float norm[3],snorm[3];

#if 0
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "computing surface normals...") ;
#endif

  for (k=0;k<mris->nfaces;k++) if (mris->faces[k].ripflag)
  {
    f = &mris->faces[k];
    for (n=0;n<4;n++)
      mris->vertices[f->v[n]].border = TRUE;
  }

  for (k=0;k<mris->nvertices;k++) if (!mris->vertices[k].ripflag)
  {
    v = &mris->vertices[k];
    snorm[0]=snorm[1]=snorm[2]=0;
    v->area = 0;
    for (n=0;n<v->num;n++)
    if (!mris->faces[v->f[n]].ripflag)
    {
      mrisNormalFace(mris, v->f[n],v->n[n],norm);
      snorm[0] += norm[0];
      snorm[1] += norm[1];
      snorm[2] += norm[2];

      /* Note: overestimates area by *2 !! */
      v->area += mrisTriangleArea(mris, v->f[n],v->n[n]); 
    }
    mrisNormalize(snorm);

    if (v->origarea<0)
      v->origarea = v->area;

    v->nx = snorm[0];
    v->ny = snorm[1];
    v->nz = snorm[2];
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
    fprintf(stderr, "reading binary curvature file...") ;

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
    fprintf(stderr, "reading binary area file...") ;

  cp = strrchr(mris_fname, '.') ;
  if (!cp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "mrisReadBinaryAreas(%s): could not scan "
                 "hemisphere from fname", mris_fname)) ;
  strncpy(hemi, cp-2, 2) ; hemi[2] = 0 ;
  FileNamePath(mris_fname, fpref) ;
  sprintf(fname, "%s/%s.area", fpref, hemi) ;

  mris->total_area = 0.0f ;
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
    mris->total_area += f;
  }
  fclose(fp);


  /* hack to correct for overestimation of area in compute_normals */
  mris->total_area /= 2; 
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "total area = %2.0f.\n", mris->total_area) ;
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
------------------------------------------------------*/
#define DEFAULT_A  45.0f
#define DEFAULT_B  130.0f
#define DEFAULT_C  75.0f
#define TOTAL_AREA 91860.3  /* area of an ellipsoid with the above A,B and C */

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
          In counter-clockwise fashion, number the vertices 0-3 and the
          edges a-d.
------------------------------------------------------*/
#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y,\
                                               v1->z-v0->z)
int
MRIScomputeFaceAreas(MRI_SURFACE *mris)
{
  int         fno, vno ;
  VERTEX      *v0, *v1, *v2, *v3 ;
  FACE        *face ;
  VECTOR      *v_a, *v_b, *v_c, *v_d, *v_n, *v_n1, *v_n2 ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_c = VectorAlloc(3, MATRIX_REAL) ;
  v_d = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_n1 = VectorAlloc(3, MATRIX_REAL) ;
  v_n2 = VectorAlloc(3, MATRIX_REAL) ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;

    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;
    v3 = &mris->vertices[face->v[3]] ;

    VERTEX_EDGE(v_a, v0, v1) ;    /* these are oriented so that the right- */
    VERTEX_EDGE(v_b, v2, v1) ;    /* hand rule will give an outward pointing */
    VERTEX_EDGE(v_c, v2, v3) ;    /* normal if the triangles have positive */
    VERTEX_EDGE(v_d, v0, v3) ;    /* area */

    VectorCrossProduct(v_a, v_d, v_n1) ;
    VectorCrossProduct(v_c, v_b, v_n2) ;
    VectorAdd(v_n1, v_n2, v_n) ;
    VectorNormalize(v_n, v_n) ;
    face->nx = VECTOR_ELT(v_n,1) ; face->ny = VECTOR_ELT(v_n,2) ;
    face->nz = VECTOR_ELT(v_n,3) ; 
    face->area = VectorTripleProduct(v_a, v_d, v_n)*.5f ;
    face->area += VectorTripleProduct(v_c, v_b, v_n)*.5f ;

    if (face->area < 0.0f)
      DiagBreak() ;
  }

  if (!mris->initialized) 
  {
    mris->total_area = 0.0f ;
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
#if 0
      face->orig_area = face->area ;
#else
      face->orig_area = 0.0 ;
      for (vno = 0 ; vno < VERTICES_PER_FACE ; vno++)
        face->orig_area += mris->vertices[face->v[vno]].origarea ;
      face->orig_area /= (float)(VERTICES_PER_FACE*2) ;
      mris->total_area += face->orig_area ;
#endif
    }
  }
  mris->initialized = 1 ;
  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_c) ;
  VectorFree(&v_d) ;
  VectorFree(&v_n) ;
  VectorFree(&v_n1) ;
  VectorFree(&v_n2) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_NEGATIVE        500
#define DELTA_T             0.05f
#define MOMENTUM            0.2f
#define ONE_MINUS_MOMENTUM  (1.0f - MOMENTUM)
#define RECOMPUTE_NORMALS   10
#define DT_DECREASE         0.7f
#define DT_INCREASE         1.03f
#define ERROR_RATIO         1.0f
#define MIN_DT              0.000001f
#define MAX_MIN_DTS         10

#define WRITE_ITERATIONS    500

MRI_SURFACE *
MRISunfold(MRI_SURFACE *mris, int niterations, float base_momentum, 
           float l_area, float l_angle,  float l_corr)
{
  int         vno, t, fno, imin, min_dts, write_iterations = WRITE_ITERATIONS,
              negative ;
  VERTEX      *v, *v0, *v1, *v2, *v3 ;
  FACE        *face ;
  VECTOR      *v_a, *v_b, *v_c, *v_d, *v_n, *v_delta, *v_sum,
              *v_d_x_n, *v_n_x_c, *v_n_x_a, *v_b_x_n ;
  float       dz, delta, rms, old_rms, delta_t,
              momentum, one_minus_momentum, min_dz, coef ;

  {
    char *cp ;

    cp = getenv("WRITE_ITERATIONS") ;
    if (cp)
    {
      sscanf(cp, "%d", &write_iterations) ;
      fprintf(stderr, "writing out intermediate files every %d iterations\n",
              write_iterations) ;
    }
  }

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_c = VectorAlloc(3, MATRIX_REAL) ;
  v_d = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_sum = VectorAlloc(3, MATRIX_REAL) ;       
  v_delta = VectorAlloc(3, MATRIX_REAL) ;   
  v_n_x_c = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_n_x_a = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_d_x_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_b_x_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */

  MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
  rms = old_rms = mrisComputeAreaRMS(mris) ;
  delta_t = DELTA_T ;
  if (base_momentum < 0.0f)
    base_momentum = MOMENTUM ; 
  momentum = base_momentum ;
  one_minus_momentum = 1.0f - momentum ;
  min_dts = 0 ;
  imin = -1 ;
  for (t = 0 ; t < niterations ; t++)
  {
    negative = mrisOrientEllipsoidArea(mris) ;
    if (!((t+1)%write_iterations))
    {
      char fname[100] ;
      sprintf(fname, "surf/%s.tmp%d", mris->hemisphere == LEFT_HEMISPHERE ?
              "lh" : "rh", t+1) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "writing %s...", fname) ;
      MRISwrite(mris, fname) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "done.\n") ;
    }
    fprintf(stdout, "%d of %d: dt = %2.6f, rms: %2.7f, %d < 0\n", 
            t, niterations, delta_t, 100.0f*rms, negative) ;
    fflush(stdout) ;
    min_dz = 1000.0f ;

    /* keep track of last deltas for momentum, and clear new deltas */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->odx = v->dx ; v->ody = v->dy ; v->odz = v->dz ;
      v->dx = v->dy = v->dz = 0.0f ;
    }

    /* calculcate movement of vertices for each face */
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      if ((negative > MAX_NEGATIVE) && face->area > 0.0f)
        continue ;

      dz = face->nx*face->nx + face->ny*face->ny ;

      v0 = &mris->vertices[face->v[0]] ;
      v1 = &mris->vertices[face->v[1]] ;
      v2 = &mris->vertices[face->v[2]] ;
      v3 = &mris->vertices[face->v[3]] ;
      
      VECTOR_LOAD(v_n, face->nx, face->ny, face->nz) ;

      VERTEX_EDGE(v_a, v0, v1) ;  /* these are oriented so that the right- */
      VERTEX_EDGE(v_b, v2, v1) ;  /* hand rule will give an outward pointing */
      VERTEX_EDGE(v_c, v2, v3) ;  /* normal if the triangles have positive */
      VERTEX_EDGE(v_d, v0, v3) ;  /* area */

      delta = face->area - face->orig_area ;

      VectorCrossProduct(v_b, v_n, v_b_x_n) ;
      VectorCrossProduct(v_d, v_n, v_d_x_n) ;
      VectorCrossProduct(v_n, v_a, v_n_x_a) ;
      VectorCrossProduct(v_n, v_c, v_n_x_c) ;

if (DEBUG_FACE(fno))
{
  DiagBreak() ;
  fprintf(stderr, 
          "\nface %d, normal (%2.1f, %2.1f, %2.1f), area %2.2f, "
          "original area %2.2f\n",
          fno, face->nx, face->ny, face->nz, face->area, face->orig_area) ;
  fprintf(stderr, "P0 %d at (%2.2f, %2.2f, %2.2f)\n",
          face->v[0], v0->x, v0->y, v0->z) ;
  fprintf(stderr, "P1 %d at (%2.2f, %2.2f, %2.2f)\n",
          face->v[1], v1->x, v1->y, v1->z) ;
  fprintf(stderr, "P2 %d at (%2.2f, %2.2f, %2.2f)\n",
          face->v[2], v2->x, v2->y, v2->z) ;
  fprintf(stderr, "P3 %d at (%2.2f, %2.2f, %2.2f)\n",
          face->v[3], v3->x, v3->y, v3->z) ;
  fprintf(stderr, "a = ") ;
  MatrixPrintTranspose(stderr, v_a) ;
  fprintf(stderr, "b = ") ;
  MatrixPrintTranspose(stderr, v_b) ;
  fprintf(stderr, "c = ") ;
  MatrixPrintTranspose(stderr, v_c) ;
  fprintf(stderr, "d = ") ;
  MatrixPrintTranspose(stderr, v_d) ;
  fprintf(stderr, "b x n = ") ;
  MatrixPrintTranspose(stderr, v_b_x_n) ;
  fprintf(stderr, "d x n = ") ;
  MatrixPrintTranspose(stderr, v_d_x_n) ;
  fprintf(stderr, "n x a = ") ;
  MatrixPrintTranspose(stderr, v_n_x_a) ;
  fprintf(stderr, "n x c = ") ;
  MatrixPrintTranspose(stderr, v_n_x_c) ;
}
      /* calculate movements in order 0-3 */
      VectorAdd(v_d_x_n, v_n_x_a, v_sum) ;
      VectorScalarMul(v_sum, delta, v_delta) ;
if (DEBUG_FACE(fno))
{ fprintf(stderr, "delta v0 = ") ; MatrixPrintTranspose(stderr, v_delta);}
      v0->dx += VECTOR_ELT(v_delta, 1) ;
      v0->dy += VECTOR_ELT(v_delta, 2) ;
      v0->dz += VECTOR_ELT(v_delta, 3) ;

      /* P1 */
      VectorAdd(v_d_x_n, v_n_x_c, v_sum) ;
      VectorScalarMul(v_sum, -delta, v_delta) ;
if (DEBUG_FACE(fno))
{ fprintf(stderr, "delta v1 = ") ; MatrixPrintTranspose(stderr, v_delta);}
      v1->dx += VECTOR_ELT(v_delta, 1) ;
      v1->dy += VECTOR_ELT(v_delta, 2) ;
      v1->dz += VECTOR_ELT(v_delta, 3) ;

      /* P2 */
      VectorAdd(v_b_x_n, v_n_x_c, v_sum) ;
      VectorScalarMul(v_sum, delta, v_delta) ;
if (DEBUG_FACE(fno))
{ fprintf(stderr, "delta v2 = ") ; MatrixPrintTranspose(stderr, v_delta);}
      v2->dx += VECTOR_ELT(v_delta, 1) ;
      v2->dy += VECTOR_ELT(v_delta, 2) ;
      v2->dz += VECTOR_ELT(v_delta, 3) ;

      /* P3 */
      VectorAdd(v_b_x_n, v_n_x_a, v_sum) ;
      VectorScalarMul(v_sum, -delta, v_delta) ;
if (DEBUG_FACE(fno))
{ fprintf(stderr, "delta v3 = ") ; MatrixPrintTranspose(stderr, v_delta);}
      v3->dx += VECTOR_ELT(v_delta, 1) ;
      v3->dy += VECTOR_ELT(v_delta, 2) ;
      v3->dz += VECTOR_ELT(v_delta, 3) ;


      if (dz < min_dz && face->area < 0)
      {
        min_dz = dz ;
        imin = fno ;
      }
    }    /* done with all faces */

    /* now apply the deltas */
    coef = delta_t * l_area * one_minus_momentum ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->odx = v->dx = momentum * v->odx + coef * v->dx ; 
      v->ody = v->dy = momentum * v->ody + coef * v->dy ; 
      v->odz = v->dz = momentum * v->odz + coef * v->dz ; 
      v->x += v->dx ; v->y += v->dy ; v->z += v->dz ;
    }
if (DEBUG_FACE(79379))
{
  face = &mris->faces[79379] ;
  v0 = &mris->vertices[face->v[0]] ;
  v1 = &mris->vertices[face->v[1]] ;
  v2 = &mris->vertices[face->v[2]] ;
  v3 = &mris->vertices[face->v[3]] ;
  fprintf(stderr, "moving vertices by:\n") ;
  fprintf(stderr, "v0: (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f)\n",
          v0->dx, v0->dy, v0->dz, v0->x, v0->y, v0->z) ;
  fprintf(stderr, "v1: (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f)\n",
          v1->dx, v1->dy, v1->dz, v1->x, v1->y, v1->z) ;
  fprintf(stderr, "v2: (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f)\n",
          v2->dx, v2->dy, v2->dz, v2->x, v2->y, v2->z) ;
  fprintf(stderr, "v3: (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f)\n",
          v3->dx, v3->dy, v3->dz, v3->x, v3->y, v3->z) ;
}

    MRISupdateEllipsoidSurface(mris) ;

if (DEBUG_FACE(79379))
{
  face = &mris->faces[79379] ;
  fprintf(stderr, "new area: %2.3f\n", face->area) ;
}
    rms = mrisComputeAreaRMS(mris) ;
    if (rms > old_rms)   /* error has increased, reduce the time step */
      delta_t *= DT_DECREASE ;
    else
      delta_t *= DT_INCREASE ;  /* error decreased - increase time step */

    if (rms > old_rms * ERROR_RATIO)  /* big increase - kill momentum */
    {
      momentum = 0.0f ; one_minus_momentum = 1.0f ;
      /* restore old vertex positions */
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        v->x -= v->dx ; v->y -= v->dy ; v->z -= v->dz ;
        v->dx = v->dy = v->dz = v->odx = v->ody = v->odz = 0.0f ;
      }
      MRISupdateEllipsoidSurface(mris) ;
      rms = mrisComputeAreaRMS(mris) ;
    }
    else   /* accept the change */
    {
      MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
      momentum = base_momentum ; one_minus_momentum = 1.0f - momentum ;
    }
    old_rms = rms ;
    if (delta_t <= MIN_DT)
    {
      delta_t = MIN_DT ;
      if (min_dts++ >= MAX_MIN_DTS)
        break ;
    }
#if 0
    else
      min_dts = 0 ;
#endif

#if 1
    if (!t)
      fprintf(stderr, "mindz = %2.3f at %d\n", min_dz, imin) ;
#endif
  }

  MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_c) ;
  VectorFree(&v_d) ;
  VectorFree(&v_delta) ;
  VectorFree(&v_sum) ;
  VectorFree(&v_n) ;

  VectorFree(&v_d_x_n) ;
  VectorFree(&v_b_x_n) ;
  VectorFree(&v_n_x_a) ;
  VectorFree(&v_n_x_c) ;
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
  float sse, delta ;
  int   fno ;
  FACE  *face ;

  sse = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    delta = face->area - face->orig_area ; 
    sse += delta*delta ;
    if (!finite(sse))
    {
      static int first = 1;
      
      if (first)
        fprintf(stderr, "sse is not finite at face %d!\n", fno) ;
      first = 0 ;
        DiagBreak() ;
      exit(1) ;
    }
  }
  return(sqrt(sse) / (float)mris->nfaces) ;
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
  int     fno ;
  FACE    *face ;

  area_scale = mris->total_area / TOTAL_AREA ;
  if (mris->initialized < 2)
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      face = &mris->faces[fno] ;
      
      /* scale the area by the ratio of the ellipsoid area to that of the
         original surface.
         */
      face->area *= area_scale ;
    }
  mris->initialized = 2 ;
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
  int     fno, negative ;
  FACE    *face ;
  VERTEX  *v ;
  float   dot ;

  negative = 0 ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    v = &mris->vertices[face->v[0]] ;  /* arbitrary */
      
    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
       */
    dot = v->x * face->nx + v->y * face->ny + v->z * face->nz;
    if (DEBUG_FACE(fno))
      fprintf(stderr, "face %d v at (%2.3f, %2.3f, %2.3f)\nface normal "
              "(%2.3f, %2.3f, %2.3f), dot = %2.4f\n",
                fno, v->x, v->y, v->z, face->nx, face->ny, face->nz, dot) ;
    if (dot < 0.0f)   /* not in same direction, area < 0 and reverse n */
    {
      negative++ ;
      face->area *= -1.0f ;
      face->nx *= -1.0f ; face->ny *= -1.0f ; face->nz *= -1.0f ;
    }
    else
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
  mrisComputeNormals(mris);
  MRIScomputeFaceAreas(mris) ;
  mrisScaleEllipsoidArea(mris) ;
  return(mrisOrientEllipsoidArea(mris)) ;
}
