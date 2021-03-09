/**
 * @brief utils
 *
 */
/*
 * Original Author: Peng Yu
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

#include "ANN.h"

extern "C"
{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "mri.h"
#include "mrisurf.h"
#include "icosahedron.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "timer.h"
#include "mrinorm.h"
#include "cma.h"
#include "version.h"
#include "error.h"
#include "matrix.h"
}

#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y, v1->z-v0->z)
typedef struct _double_3d
{
  double x;
  double y;
  double z;
}
double3d ;

#define TINY 1.0e-20;
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
 a[k][l]=h+s*(g-h*tau);

int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
const char      *Progname ;
static MRI_SURFACE *center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst);
static MRI_SURFACE *sample_origposition(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
static MRI_SURFACE *sample_origcurvature(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
static double   v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug) ;
static double   brain_volume(MRI_SURFACE *mris);
//static int      sort(double **array, int order, int number, int total_number);

static int      CURV=0;

int
main(int argc, char *argv[])
{
  int           nargs, msec, order, i;
  Timer then ;
  MRIS          *mris_in, *mris_out;
  double        volume; // T[5][5] the transformation matrx (using index 1-4)

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <input sphere> <input surface> <int which> <output wavelets> ", Progname);

  then.reset() ;

  //order = atoi (argv[3]);
  order = 7;
  //fprintf(stdout, "Set %s as the finest scale level\n", argv[3]);
  if (order > 8)
    ErrorExit(ERROR_BADPARM, "the highest order is 7\n");

  if (atoi(argv[3]) == 1) CURV=1;

  /*Spherical Wavelet Analysis*/

  if ( !CURV)
  {
    /* using talairach transformation */
    fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);
    mris_in = MRISread(argv[1]) ;
    if (!mris_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, argv[1]) ;
    MRISreadOriginalProperties(mris_in, argv[2]) ;
    fprintf(stdout, "Reading original surface from %s\n", argv[2]);
    MRISsaveVertexPositions(mris_in, TMP_VERTICES);
    MRISrestoreVertexPositions(mris_in, ORIGINAL_VERTICES);

    /* save the transformed original coordinates */
    center_brain(mris_in, mris_in);
    MRISsaveVertexPositions(mris_in, ORIGINAL_VERTICES);
    MRISrestoreVertexPositions(mris_in, TMP_VERTICES);
    /* sample the surface onto ic7 */
    mris_out = ReadIcoByOrder(order, 100);
    sample_origposition(mris_in, mris_out) ;

    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    center_brain(mris_out, mris_out);
    volume=brain_volume(mris_out);
    //MRISscaleBrain(mris_out, mris_out, cbrt(300000/volume)) ;

    MRISwrite(mris_out,argv[4]) ;
    fprintf(stdout,"Write sampled surface to %s\n",argv[4]);

    MRISfree(&mris_in) ;
    MRISfree(&mris_out) ;
    /*End of Analysis*/
  }
  else if ( CURV)
  {
    mris_in = MRISread(argv[1]) ;
    if (!mris_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, argv[1]) ;
    fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);

    MRISreadCurvatureFile(mris_in, argv[2]) ;
    fprintf(stdout, "Reading input from %s\n", argv[2]);

    mris_out = ReadIcoByOrder(order, 100);
    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    sample_origcurvature(mris_in, mris_out) ;
    MRISwriteCurvature(mris_out,argv[4]);
    fprintf(stdout, "Writing sampled curvature of original surface to %s\n", argv[4]);

    MRISfree(&mris_in) ;
    MRISfree(&mris_out) ;
  }

  msec = then.milliseconds() ;
  fprintf(stdout, "spherical sampling took %2.1f minutes\n", (float)msec/(1000.0f*60.0f));
  exit(0) ;
  return(0) ;
}


double
brain_volume(MRI_SURFACE *mris)
{
  int fno;
  FACE *face;
  double total_volume, face_area;
  VECTOR *v_a, *v_b, *v_n, *v_cen;
  VERTEX  *v0, *v1, *v2;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_cen = VectorAlloc(3, MATRIX_REAL) ;     /* centroid vector */

  total_volume = 0;
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

    /* face normal vector */
    V3_CROSS_PRODUCT(v_a, v_b, v_n) ;
    face_area = V3_LEN(v_n) * 0.5f ;

    V3_NORMALIZE(v_n, v_n) ;             /* make it a unit vector */

    /* compute face centroid */
    V3_X(v_cen) = (v0->x + v1->x + v2->x)/3.0;
    V3_Y(v_cen) = (v0->y + v1->y + v2->y)/3.0;
    V3_Z(v_cen) = (v0->z + v1->z + v2->z)/3.0;

    total_volume += V3_DOT(v_cen, v_n)*face_area;
  }

  total_volume /= 3.0;

  return(total_volume);
}


static MRI_SURFACE *
center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int         fno, vno ;
  FACE        *face ;
  VERTEX      *vdst;
  double       x, y, z, x0, y0, z0 ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  x0 = y0 = z0 = 0 ;   /* silly compiler warning */
  //MRISupdateSurface(mris_dst);
#if 1 //area weighting
  for (fno = 0 ; fno < mris_src->nfaces ; fno++)
  {
    face = &mris_dst->faces[fno] ;
    if (face->ripflag)
      continue ;
    x = mris_dst->vertices[face->v[0]].x;
    y = mris_dst->vertices[face->v[0]].y;
    z = mris_dst->vertices[face->v[0]].z;
    x += mris_dst->vertices[face->v[1]].x;
    y += mris_dst->vertices[face->v[1]].y;
    z += mris_dst->vertices[face->v[1]].z;
    x += mris_dst->vertices[face->v[2]].x;
    y += mris_dst->vertices[face->v[2]].y;
    z += mris_dst->vertices[face->v[2]].z;
    x = x*face->area/3;
    y = y*face->area/3;
    z = z*face->area/3;
    x0 += x/mris_dst->total_area;
    y0 += y/mris_dst->total_area;
    z0 += z/mris_dst->total_area;
    if (face->area == 0) fprintf(stdout, "%d %f\n", fno, face->area);
  }
#else
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno];
    x = vdst->x;
    y = vdst->y;
    z = vdst->z;
    x0 += x/mris_dst->nvertices;
    y0 += y/mris_dst->nvertices;
    z0 += z/mris_dst->nvertices;
  }
#endif
  fprintf(stdout, "brain center is found at %f %f %f\n", x0, y0, z0);
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->ripflag)
      continue ;
    vdst->x -= x0 ;
    vdst->y -= y0 ;
    vdst->z -= z0 ;
  }

  mris_dst->xctr = mris_dst->yctr = mris_dst->zctr = 0 ;
  return(mris_dst) ;
}


static MRI_SURFACE *
sample_origposition(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int index, fno, fnum=0, i;
  VERTEX *vertex;
  double nearest, dist, r, s, t;
  double a, b, c, p;
  ANNpointArray pa = annAllocPts(mris_src->nvertices, 3);

  for (index = 0; index < mris_src->nvertices; index++)
  {
    pa[index][0] = mris_src->vertices[index].x;
    pa[index][1] = mris_src->vertices[index].y;
    pa[index][2] = mris_src->vertices[index].z;
  }

  ANNkd_tree *annkdTree = new ANNkd_tree(pa, mris_src->nvertices, 3);
  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];
  //  ANNpoint query_pt = annAllocPt(3);
  ANNpointArray QueryPt;

  //if(mris_dst == NULL)
  // mris_dst = MRIScopy(mris_src, NULL);

  QueryPt = annAllocPts(1,3);

  for (index = 0; index < mris_dst->nvertices; index++)
  {
    if (mris_dst->vertices[index].border == 1) continue;

    QueryPt[0][0] = mris_dst->vertices[index].x;
    QueryPt[0][1] = mris_dst->vertices[index].y;
    QueryPt[0][2] = mris_dst->vertices[index].z;

    annkdTree->annkSearch( // search
      QueryPt[0],       // query point
      1,   // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,  // distance (returned)
      0);   // error bound

#if 1
    vertex = &mris_src->vertices[annIndex[0]];
    nearest = 100000;
    for (i=0; i<vertex->num; i++)
    {
      fno = vertex->f[i];
      dist = v_to_f_distance(&mris_dst->vertices[index], mris_src, fno, 0);
      if (dist<nearest)
      {
        nearest = dist;
        fnum = fno;
      }
    }

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_src->vertices[mris_src->faces[fnum].v[1]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_src->vertices[mris_src->faces[fnum].v[1]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_src->vertices[mris_src->faces[fnum].v[1]].z),2));
    p = (a+b+c)/2;
    r = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_src->vertices[mris_src->faces[fnum].v[0]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_src->vertices[mris_src->faces[fnum].v[0]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_src->vertices[mris_src->faces[fnum].v[0]].z),2));
    p = (a+b+c)/2;
    s = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_src->vertices[mris_src->faces[fnum].v[1]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_src->vertices[mris_src->faces[fnum].v[1]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_src->vertices[mris_src->faces[fnum].v[1]].z),2));
    p = (a+b+c)/2;
    t = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));
    p= (r+s+t);
    r = r/p;
    s = s/p;
    t = t/p;

    mris_dst->vertices[index].origx = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origx + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origx + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origx;
    mris_dst->vertices[index].origy = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origy + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origy + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origy;
    mris_dst->vertices[index].origz = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origz + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origz + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origz;

#else
    mris_dst->vertices[index].origx = mris_src->vertices[annIndex[0]].origx;
    mris_dst->vertices[index].origy = mris_src->vertices[annIndex[0]].origy;
    mris_dst->vertices[index].origz = mris_src->vertices[annIndex[0]].origz;
#endif
    if (annIndex[0] == 81426)
      printf("src index %d dst index %d face %d r %f s %f t %f\n", index, annIndex[0], fnum, r, s, t);
  }
  return(mris_dst) ;
}

static MRI_SURFACE *
sample_origcurvature(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int index, fno, fnum=0, i;
  VERTEX *vertex;
  double nearest, dist, r, s, t;
  double a, b, c, p;
  ANNpointArray pa = annAllocPts(mris_src->nvertices, 3);

  for (index = 0; index < mris_src->nvertices; index++)
  {
    pa[index][0] = mris_src->vertices[index].x;
    pa[index][1] = mris_src->vertices[index].y;
    pa[index][2] = mris_src->vertices[index].z;
  }

  ANNkd_tree *annkdTree = new ANNkd_tree(pa, mris_src->nvertices, 3);
  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];
  //  ANNpoint query_pt = annAllocPt(3);
  ANNpointArray QueryPt;

  //if(mris_dst == NULL)
  // mris_dst = MRIScopy(mris_src, NULL);

  QueryPt = annAllocPts(1,3);

  for (index = 0; index < mris_dst->nvertices; index++)
  {
    if (mris_dst->vertices[index].border == 1) continue;

    QueryPt[0][0] = mris_dst->vertices[index].x;
    QueryPt[0][1] = mris_dst->vertices[index].y;
    QueryPt[0][2] = mris_dst->vertices[index].z;

    annkdTree->annkSearch( // search
      QueryPt[0],       // query point
      1,   // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,  // distance (returned)
      0);   // error bound

#if 1
    vertex = &mris_src->vertices[annIndex[0]];
    nearest = 100000;
    for (i=0; i<vertex->num; i++)
    {
      fno = vertex->f[i];
      dist = v_to_f_distance(&mris_dst->vertices[index], mris_src, fno, 0);
      if (dist<nearest)
      {
        nearest = dist;
        fnum = fno;
      }
    }

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_src->vertices[mris_src->faces[fnum].v[1]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_src->vertices[mris_src->faces[fnum].v[1]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_src->vertices[mris_src->faces[fnum].v[1]].z),2));
    p = (a+b+c)/2;
    r = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_src->vertices[mris_src->faces[fnum].v[0]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_src->vertices[mris_src->faces[fnum].v[0]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_src->vertices[mris_src->faces[fnum].v[0]].z),2));
    p = (a+b+c)/2;
    s = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_src->vertices[mris_src->faces[fnum].v[1]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_src->vertices[mris_src->faces[fnum].v[1]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_src->vertices[mris_src->faces[fnum].v[1]].z),2));
    p = (a+b+c)/2;
    t = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    p= (r+s+t);
    r = r/p;
    s = s/p;
    t = t/p;

    mris_dst->vertices[index].curv = r*mris_src->vertices[mris_src->faces[fnum].v[0]].curv + s*mris_src->vertices[mris_src->faces[fnum].v[1]].curv + t*mris_src->vertices[mris_src->faces[fnum].v[2]].curv;
#else

    mris_dst->vertices[index].curv = mris_src->vertices[annIndex[0]].curv;
    if (index == 67770)
      printf("src index %d dst index %d face %d r %f s %f t %f\n", index, annIndex[0], fnum, r, s, t);
#endif

  }
  return(mris_dst) ;
}


static double
v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug)
{
  double a, b, c, d, e, f, det, s, t, invDet;
  double numer, denom, tmp0, tmp1;
  VERTEX *V1, *V2, *V3;
  FACE *face;

  VERTEX E0, E1, D;

  face = &mri_surf->faces[face_number];
  V1 = &mri_surf->vertices[face->v[0]];
  V2 = &mri_surf->vertices[face->v[1]];
  V3 = &mri_surf->vertices[face->v[2]];

  E0.x = V2->x - V1->x;
  E0.y = V2->y - V1->y;
  E0.z = V2->z - V1->z;
  E1.x = V3->x - V1->x;
  E1.y = V3->y - V1->y;
  E1.z = V3->z - V1->z;
  D.x = V1->x - P0->x;
  D.y = V1->y - P0->y;
  D.z = V1->z - P0->z;

  a = E0.x *E0.x + E0.y * E0.y + E0.z *E0.z;
  b = E0.x *E1.x + E0.y * E1.y + E0.z *E1.z;
  c = E1.x *E1.x + E1.y * E1.y + E1.z *E1.z;
  d = E0.x *D.x + E0.y * D.y + E0.z *D.z;
  e = E1.x *D.x + E1.y * D.y + E1.z *D.z;
  f = D.x *D.x + D.y * D.y + D.z *D.z;

  det = a*c - b*b;
  s = b*e - c*d;
  t = b*d - a*e;

  if (debug) printf("det = %g\n", det);
  if (s + t <= det)
  {
    if (s < 0)
    {
      if (t<0)
      {
        /* Region 4 */
        tmp0 = b + d;
        tmp1 = c + e;
        if (tmp1 > tmp0)
        {
          numer = tmp1 - tmp0;
          denom = a - b - b +c;
          s = (numer >= denom ? 1 : numer/denom);
          t = 1-s;

        }
        else
        {
          s = 0;
          /* t = (e >= 0 ? 0 : (-e >= c ? 0 > = c + e = tmp1) */
          t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
        }
        if (debug) printf("region 4, s =%g, t =%g\n", s, t);
      }
      else
      {
        /* Region 3 */
        s = 0;
        t = ( e >= 0 ? 0 : (-e >= c ? 1 : (-e/c)));
        if (debug) printf("region 3, s =%g, t =%g\n", s, t);
      }
    }
    else if (t < 0)
    {
      /* Region 5 */
      t = 0;
      s = (d >= 0 ? 0 :(-d >= a ? 1 : (-d/a)));
      if (debug) printf("region 5, s =%g, t =%g\n", s, t);
    }
    else
    {
      /* Region 0 */
      invDet = 1/det;
      s *= invDet;
      t *= invDet;
      if (debug) printf("region 0, s =%g, t =%g\n", s, t);
    }

  }
  else
  {
    if ( s < 0 )
    {
      /* Region 2 */
      if ( d < 0)
      { /* Minimum on edge t = 0 */
        s = (-d >= a ? 1 : -d/a);
        t = 0;
      }
      else if (e < 0)
      { /* Minimum on edge s = 0 */
        t = (-e >= c ? 1 : -e/c);
        s = 0;
      }
      else
      {
        s = 0;
        t = 0;
      }
      if (debug) printf("region 2, s =%g, t =%g\n", s, t);
    }
    else if ( t < 0)
    {
      /* Region 6 */
      tmp0 = b + e;
      tmp1 = a + d;
      if (tmp1 > tmp0)
      { /* Minimum at line s + t = 1 */
        numer = tmp1 - tmp0; /* Positive */
        denom = a + c - b -b;
        t = (numer >= denom ? 1 : (numer/denom));
        s = 1 - t;
      }
      else
      { /* Minimum at line t = 0 */
        s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : -d/a));
        t = 0;
      }
      if (debug) printf("region 6, s =%g, t =%g\n", s, t);
    }
    else
    {
      /* Region 1 */
      numer = c + e - b - d;
      if (numer <= 0)
      {
        s = 0;
      }
      else
      {
        denom = a + c - b - b; /* denom is positive */
        s = (numer >= denom ? 1 : (numer/denom));
      }
      t = 1-s;
      if (debug) printf("region 1, s =%g, t =%g\n", s, t);
    }
  }

  if ( s < 0  || s > 1 || t < 0 || t > 1)
  {
    printf("Error in computing s and t \n");
  }


  /* return (sqrt(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f)); */
  /* The square-root will be taken later to save time */
  return (a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f);
}


#if 0
static int
sort(double **array, int order, int number, int total_number)
{
  int i, j;
  double temp, index;

  if (order)
  {
    for (j=0; j<number; j++)
      for (i=total_number-1; i>j ; i--)
        if (array[0][i]<array[0][i-1])
        {
          temp=array[0][i];
          array[0][i]=array[0][i-1];
          array[0][i-1]=temp;
          index=array[1][i];
          array[1][i]=array[1][i-1];
          array[1][i-1]=index;
        }
  }
  else
  {
    for (j=0; j<number; j++)
      for (i=total_number-1; i>j ; i--)
        if (array[0][i]>array[0][i-1])
        {
          temp=array[0][i];
          array[0][i]=array[0][i-1];
          array[0][i-1]=temp;
          index=array[1][i];
          array[1][i]=array[1][i-1];
          array[1][i-1]=index;
        }
  }
  return(NO_ERROR);
}
#endif
/*----------------------------------------------------------------------

 Parameters:

 Description:
 ----------------------------------------------------------------------*/

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  if (!stricmp(option, "CURV"))
  {
    CURV = 1;
    fprintf(stdout,"Read in the curvature file and decompose it\n");
  }
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      fprintf(stdout, "usage: %s <input spherical  surface> <orig surface> <mode> <output surface> \n", Progname);
      exit(1) ;
      break ;
    default:
      fprintf(stdout, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}






