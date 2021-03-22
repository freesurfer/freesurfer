/*
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


///////////////////////////////////////////
// mris_spharm.c
//
// written by Peng Yu
// date: 12/09/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
////////////////////////////////////////////
#include "ANN.h"

 
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


int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
const char            *Progname ;
static int      SphericalHarmonics(int L, int M, float theta, float phi) ;
static double   legendre(int l, int m, float x) ;
static double   factorial(int L, int M) ;
static MATRIX   *ComplexMatrixTranspose(MATRIX *mIn, MATRIX *mOut);
static MATRIX   *ComplexMatrixPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv);
static MRI_SURFACE *center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
static MRI_SURFACE *sample_origposition(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
static double   v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug) ;
static int      LEVEL = 7;
static double   REAL, IMAG;
static float    threshold = 1;
static int      COMPARE = 0;
static char     *cfname;

int
main(int argc, char *argv[]) {
  char          fname[STRLEN];
  int           nargs, msec, order, i, j, k, nvertices, dimension, count, fno;
  float         phi, d, theta, area;
  Timer then ;
  MRIS          *mris_in, *mris_out;
  //MRI_SP        *mrisp ;
  MATRIX        *m_Z, *m_V, *m_c, *m_Z_inv, *m_V_new;
  VERTEX        *v;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <input surface> <orig surface> <finest order> <output surface>", Progname);

  then.reset() ;

  mris_in = MRISread(argv[1]) ;
  if (!mris_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, argv[1]) ;
  fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);
  MRISreadOriginalProperties(mris_in, argv[2]) ;
  fprintf(stdout, "Reading original surface from %s\n", argv[2]);

  order = atoi (argv[3]);
  fprintf(stdout, "SPHARM decomposition up tp level %s\n", argv[3]);

  /*Sample the original surface*/

  mris_out = ReadIcoByOrder(LEVEL, 100);
  //MRISwrite(mris_out, "/space/xrt/1/users/btquinn/buckner_paper/010223_61223/surf/lh.ic7") ;
  for (i = 0; i<mris_out->nvertices; i++)
    mris_out->vertices[i].nsize=1;
  //mrisp = MRISPalloc(1, 3);
#if 1
  //MRIScoordsToParameterization(mris_in, mrisp, 1, ORIGINAL_VERTICES) ;
  //MRISPblur(mrisp, mrisp, 0.5, 0);
  //MRISPblur(mrisp, mrisp, 0.5, 1);
  //MRISPblur(mrisp, mrisp, 0.5, 2);
  //MRIScoordsFromParameterization(mrisp, mris_out) ;
  sample_origposition(mris_in, mris_out) ;
#else
  MRISreadOriginalProperties(mris_out, argv[2]) ;
#endif
#if 1 /*just to test if the parameterization is correct */
  MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
  MRISupdateSurface(mris_out);
  MRISwrite(mris_out, "/space/xrt/1/users/btquinn/buckner_paper/010223_61223/surf/lh.sampled") ;
  fprintf(stderr, "original area becomes %f\n", mris_out->total_area);
  center_brain(mris_out, mris_out);
  MRISscaleBrain(mris_out, mris_out, sqrt(100000.0f/mris_out->total_area)) ;
  MRISupdateSurface(mris_out);
  for (fno=0; fno<mris_out->nfaces; fno++)
    area += mris_out->faces[fno].area;
  fprintf(stderr, "original area becomes %f\n", area);
  MRISsaveVertexPositions(mris_out, ORIGINAL_VERTICES) ;
  MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
#endif

  /*Initialize Matrix*/
  dimension = (order+1)*(order+1);
  nvertices = mris_out->nvertices;
  m_c = MatrixAlloc(dimension, 3, MATRIX_COMPLEX) ;
  m_Z = MatrixAlloc(nvertices, dimension, MATRIX_COMPLEX) ;
  m_V = MatrixAlloc(nvertices, 3, MATRIX_REAL) ;
  m_Z_inv = MatrixAlloc(dimension, nvertices, MATRIX_COMPLEX) ;
  m_V_new = MatrixAlloc(nvertices, 3, MATRIX_COMPLEX) ;

  for (i=0; i<nvertices; i++) {
    v = &mris_out->vertices[i] ;
    phi = atan2(v->y, v->x) ;
    if (phi < 0.0)  phi = 2 * M_PI + phi ;
    d = 100*100 - v->z*v->z ;
    if (d < 0.0)   d = 0 ;
    theta = atan2(sqrt(d), v->z) ;
    count = 0;
    *MATRIX_RELT(m_V,i+1,1) = v->origx;
    *MATRIX_RELT(m_V,i+1,2) = v->origy;
    *MATRIX_RELT(m_V,i+1,3) = v->origz;

    for (j=0; j<=order; j++)
      for ( k= -1*j; k<=j; k++) {
        count++ ;
        SphericalHarmonics(j, k, theta, phi);
        MATRIX_CELT_REAL(m_Z,i+1,count) = REAL;
        MATRIX_CELT_IMAG(m_Z,i+1,count) = IMAG;
      }
  }

  m_Z_inv=ComplexMatrixPseudoInverse(m_Z, NULL) ;
  MatrixMultiply(m_Z_inv, m_V, m_c) ;
  MatrixPrint(stdout,m_c);

  if (COMPARE) {
    MATRIX *m_cc;
    double r1, r2, m1, m2, diff;
    fprintf(stdout, "%s\n", cfname);
    m_cc = MatrixAlloc(dimension, 3, MATRIX_COMPLEX) ;
    mris_out = ReadIcoByOrder(LEVEL, 100);
    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    //mrisp = MRISPalloc(1, 3);
#if 1
    //MRIScoordsToParameterization(mris_in, mrisp, 1, ORIGINAL_VERTICES) ;
    //MRISPblur(mrisp, mrisp, 0.5, 0);
    //MRISPblur(mrisp, mrisp, 0.5, 1);
    //MRISPblur(mrisp, mrisp, 0.5, 2);
    //MRIScoordsFromParameterization(mrisp, mris_out) ;
    sample_origposition(mris_in, mris_out) ;
#else
    MRISreadOriginalProperties(mris_out, cfname) ;
#endif
#if 1 /*just to test if the parameterization is correct */
    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    MRISupdateSurface(mris_out);
    fprintf(stderr, "original area becomes %f\n", mris_out->total_area);
    center_brain(mris_out, mris_out);
    MRISscaleBrain(mris_out, mris_out, sqrt(100000.0f/mris_out->total_area)) ;
    MRISupdateSurface(mris_out);
    for (fno=0; fno<mris_out->nfaces; fno++)
      area += mris_out->faces[fno].area;
    fprintf(stderr, "original area becomes %f\n", area);
    //MRISwrite(mris_out, "/space/xrt/1/users/btquinn/buckner_paper/011121_vc8048/surf/lh.sampled") ;
    MRISsaveVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
#endif

    for (i=0; i<nvertices; i++) {
      v = &mris_out->vertices[i] ;
      phi = atan2(v->y, v->x) ;
      if (phi < 0.0)  phi = 2 * M_PI + phi ;
      d = 100*100 - v->z*v->z ;
      if (d < 0.0)   d = 0 ;
      theta = atan2(sqrt(d), v->z) ;
      count = 0;
      *MATRIX_RELT(m_V,i+1,1) = v->origx;
      *MATRIX_RELT(m_V,i+1,2) = v->origy;
      *MATRIX_RELT(m_V,i+1,3) = v->origz;

      for (j=0; j<=order; j++)
        for ( k= -1*j; k<=j; k++) {
          count++ ;
          SphericalHarmonics(j, k, theta, phi);
          MATRIX_CELT_REAL(m_Z,i+1,count) = REAL;
          MATRIX_CELT_IMAG(m_Z,i+1,count) = IMAG;
        }
    }

    m_Z_inv=ComplexMatrixPseudoInverse(m_Z, NULL) ;
    MatrixMultiply(m_Z_inv, m_V, m_cc) ;
    MatrixPrint(stdout,m_cc);

    /*check for difference*/
    for (count=0, j=0; j<=order; j++)
      for ( k= -1*j; k<=j; k++) {
        count++;
        for (i=1;i<=3;i++) {
          r1 = MATRIX_CELT_REAL(m_c,count,i);
          r2 = MATRIX_CELT_REAL(m_cc,count,i);
          m1 = MATRIX_CELT_IMAG(m_c,count,i);
          m2 = MATRIX_CELT_IMAG(m_cc,count,i);
          if (sqrt(r1*r1+m1*m1)==0)
            diff = fabs(sqrt(r1*r1+m1*m1)-sqrt(r2*r2+m2*m2));
          else
            diff = fabs(sqrt(r1*r1+m1*m1)-sqrt(r2*r2+m2*m2))/sqrt(r1*r1+m1*m1);
          if (diff>threshold) {
            MATRIX_CELT_REAL(m_c,count,i) = r2;
            MATRIX_CELT_IMAG(m_c,count,i) = m2;
            fprintf(stdout, "%d %d %f\n",count, i, diff);
          }
        }
      }
    MatrixFree(&m_cc) ;
  }

  /*Reconstruct Surface using only part of the coefficients*/
#if 0
  for (i=5; i<=m_c->rows; i++)
    for (j=1; j<=m_c->cols;j++) {
      MATRIX_CELT_REAL(m_c,i,j) = 0;
      MATRIX_CELT_IMAG(m_c,i,j) = 0;
    }
#endif
  MatrixMultiply(m_Z, m_c, m_V_new) ;

  for (i=0; i<nvertices; i++) {
    v = &mris_out->vertices[i] ;
    v->origx = MATRIX_CELT_REAL(m_V_new,i+1,1) ;
    v->origy = MATRIX_CELT_REAL(m_V_new,i+1,2) ;
    v->origz = MATRIX_CELT_REAL(m_V_new,i+1,3) ;
  }

  MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
  fprintf(stdout, "Writing reconstructed surface to %s\n", argv[4]);
  MRISwrite(mris_out,argv[4] ) ;
  MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;


  MatrixFree(&m_V_new) ;
  MatrixFree(&m_Z_inv) ;
  MatrixFree(&m_Z) ;
  MatrixFree(&m_V) ;
  MatrixFree(&m_c) ;
  //MRISPfree(&mrisp) ;
  MRISfree(&mris_in) ;
  MRISfree(&mris_out);
  msec = then.milliseconds() ;
  fprintf(stdout, "spherical wavelet took %2.1f minutes\n", (float)msec/(1000.0f*60.0f));
  exit(0) ;
  return(0) ;
}


static int
SphericalHarmonics(int L, int M, float theta, float phi) {
  double   P, factor;

  P = legendre(L, abs(M), cos(theta));
  factor = sqrt((L+0.5)*factorial(L,abs(M)));
  REAL = factor * P * cos(abs(M)*phi) /sqrt(2 * M_PI);
  IMAG = factor * P * sin(abs(M)*phi) /sqrt(2 * M_PI);
  if (M < 0) {
    REAL = pow(-1,abs(M)) * REAL;
    IMAG = pow(-1,abs(M)) *(-1) * IMAG;
  }
  return(NO_ERROR);
}

static double
legendre(int l, int m, float x) {
  double  fact, pll, pmm, pmmp1, somx2;
  int    i, ll;

  if ( m<0 || m>l || fabs(x)>1.0 )
    fprintf(stdout, "Bad arguments in routine legendre");
  pmm = 1.0;
  if ( m > 0) {
    somx2 = sqrt((1.0-x)*(1.0+x));
    fact = 1.0;
    for ( i=1 ; i<=m ; i++) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }

  if ( l==m )
    return pmm;
  else {
    pmmp1 = x*(2*m+1)*pmm;
    if ( l==(m+1) )
      return pmmp1;
    else {
      for ( ll=m+2 ; ll<=l ; ll++ ) {
        pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
        pmm = pmmp1;
        pmmp1 = pll;
      }
      return pll;
    }
  }
}

static double
factorial(int L, int M) {
  double ratio=1;
  int i;

  for (i=L+M; i>L-M; i--)
    ratio *= (double) 1/i;
  return(ratio);
}

MATRIX *
ComplexMatrixPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv) {
  MATRIX  *mT, *mTm, *mTm_inv ;

  /* build (mT m)-1 mT */
  mT = ComplexMatrixTranspose(m, NULL) ;
  mTm = MatrixMultiply(mT, m, NULL) ;
  mTm_inv = MatrixInverse(mTm, NULL) ;
  if (!mTm_inv) {
    MatrixFree(&mT) ;
    MatrixFree(&mTm) ;
    return(NULL) ;
  }
  m_pseudo_inv = MatrixMultiply(mTm_inv, mT, m_pseudo_inv) ;

  MatrixFree(&mT) ;
  MatrixFree(&mTm) ;
  MatrixFree(&mTm_inv) ;
  return(m_pseudo_inv) ;
}

MATRIX  *
ComplexMatrixTranspose(MATRIX *mIn, MATRIX *mOut) {
  int  row, col, rows, cols ;

  if (!mOut) {
    mOut = MatrixAlloc(mIn->cols, mIn->rows, mIn->type) ;
    if (!mOut)
      return(NULL) ;
  }

  rows = mIn->rows ;
  cols = mIn->cols ;

  for (row = 1 ; row <= rows ; row++) {
    for (col = 1 ; col <= cols ; col++) {
      MATRIX_CELT_REAL(mOut,col,row) = MATRIX_CELT_REAL(mIn,row,col)  ;
      MATRIX_CELT_IMAG(mOut,col,row) = -1*MATRIX_CELT_IMAG(mIn,row,col)  ;
    }
  }

  return(mOut) ;
}


static MRI_SURFACE *
center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) {
  int         fno, vno ;
  FACE        *face ;
  VERTEX      *vdst;
  float       x, y, z, x0, y0, z0 ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  x0 = y0 = z0 = 0 ;   /* silly compiler warning */
  MRISupdateSurface(mris_dst);
  for (fno = 0 ; fno < mris_src->nfaces ; fno++) {
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
    x /= face->area;
    y/= face->area;
    z/= face->area;
    x0 += x;
    y0 += y;
    z0 += z;
    if (face->area == 0) fprintf(stdout, "%d %f\n", fno, face->area);
  }
  x0 /= mris_dst->total_area ;
  y0 /= mris_dst->total_area  ;
  z0 /= mris_dst->total_area ;

  for (vno = 0 ; vno < mris_src->nvertices ; vno++) {
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
sample_origposition(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) {
  int index, vno, fno, fnum;
  VERTEX *vertex;
  double nearest, dist, x1, y1, z1, x2, y2, z2, r, s, t;
  ANNpointArray pa = annAllocPts(mris_src->nvertices, 3);

  for (index = 0; index < mris_src->nvertices; index++) {
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

  for (index = 0; index < mris_dst->nvertices; index++) {
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

    vertex = &mris_src->vertices[annIndex];
    nearest = 100000;
    for (i=0; i<vertex->num; i++) {
      fno = vertex->f[i];
      dist = v_to_f_distance(mris_dst->vertices[index], mris_src, fno, 0);
      if (dist<nearest) {
        nearest = dist;
        fnum = fno;
      }
    }

    x1 = mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x;
    y1 = mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y;
    z1 = mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z;
    x2 = mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x;
    y2 = mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y;
    z2 = mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z;
    r = 0.5*fabs(x1*x2+y1*y2+z1*z2);

    x1 = mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x;
    y1 = mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y;
    z1 = mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z;
    x2 = mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x;
    y2 = mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y;
    z2 = mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z;
    s = 0.5*fabs(x1*x2+y1*y2+z1*z2);

    x1 = mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x;
    y1 = mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y;
    z1 = mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z;
    x2 = mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x;
    y2 = mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y;
    z2 = mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z;
    t = 0.5*fabs(x1*x2+y1*y2+z1*z2);
    r = r/(r+s+t);
    s = s/(r+s+t);
    t = t/(r+s+t);

    mris_dst->vertices[index].origx = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origx + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origx + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origx;
    mris_dst->vertices[index].origy = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origy + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origy + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origy;
    mris_dst->vertices[index].origz = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origz + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origz + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origz;
    if (index == 69894)
      printf("src index %d dst index %d face %d r %f s %f t %f\n", index, annIndex[0], fnum, r, s, t);
  }
  return(mris_dst) ;
}


static double
v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug) {
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
  if (s + t <= det) {
    if (s < 0) {
      if (t<0) {
        /* Region 4 */
        tmp0 = b + d;
        tmp1 = c + e;
        if (tmp1 > tmp0) {
          numer = tmp1 - tmp0;
          denom = a - b - b +c;
          s = (numer >= denom ? 1 : numer/denom);
          t = 1-s;

        } else {
          s = 0;
          /* t = (e >= 0 ? 0 : (-e >= c ? 0 > = c + e = tmp1) */
          t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
        }
        if (debug) printf("region 4, s =%g, t =%g\n", s, t);
      } else {
        /* Region 3 */
        s = 0;
        t = ( e >= 0 ? 0 : (-e >= c ? 1 : (-e/c)));
        if (debug) printf("region 3, s =%g, t =%g\n", s, t);
      }
    } else if (t < 0) {
      /* Region 5 */
      t = 0;
      s = (d >= 0 ? 0 :(-d >= a ? 1 : (-d/a)));
      if (debug) printf("region 5, s =%g, t =%g\n", s, t);
    } else {
      /* Region 0 */
      invDet = 1/det;
      s *= invDet;
      t *= invDet;
      if (debug) printf("region 0, s =%g, t =%g\n", s, t);
    }

  } else {
    if ( s < 0 ) {
      /* Region 2 */
      if ( d < 0) { /* Minimum on edge t = 0 */
        s = (-d >= a ? 1 : -d/a);
        t = 0;
      } else if (e < 0) { /* Minimum on edge s = 0 */
        t = (-e >= c ? 1 : -e/c);
        s = 0;
      } else {
        s = 0;
        t = 0;
      }
      if (debug) printf("region 2, s =%g, t =%g\n", s, t);
    } else if ( t < 0) {
      /* Region 6 */
      tmp0 = b + e;
      tmp1 = a + d;
      if (tmp1 > tmp0) { /* Minimum at line s + t = 1 */
        numer = tmp1 - tmp0; /* Positive */
        denom = a + c - b -b;
        t = (numer >= denom ? 1 : (numer/denom));
        s = 1 - t;
      } else { /* Minimum at line t = 0 */
        s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : -d/a));
        t = 0;
      }
      if (debug) printf("region 6, s =%g, t =%g\n", s, t);
    } else {
      /* Region 1 */
      numer = c + e - b - d;
      if (numer <= 0) {
        s = 0;
      } else {
        denom = a + c - b - b; /* denom is positive */
        s = (numer >= denom ? 1 : (numer/denom));
      }
      t = 1-s;
      if (debug) printf("region 1, s =%g, t =%g\n", s, t);
    }
  }

  if ( s < 0  || s > 1 || t < 0 || t > 1) {
    printf("Error in computing s and t \n");
  }


  /* return (sqrt(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f)); */
  /* The square-root will be taken later to save time */
  return (a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f);
}

/*----------------------------------------------------------------------

 Parameters:

 Description:
 ----------------------------------------------------------------------*/


static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  if (!stricmp(option, "L")) {
    LEVEL = atoi(argv[2]) ;
    fprintf(stdout,"Use %d order icosahedron \n", LEVEL);
    nargs=1;
  } else if (!stricmp(option, "C")) {
    cfname = argv[2] ;
    COMPARE = 1;
    fprintf(stdout,"Compare with %s \n", cfname);
    nargs=1;
  } else if (!stricmp(option, "T")) {
    threshold = atof(argv[2]) ;
    fprintf(stdout,"Set threshold as %f \n", threshold);
    nargs=1;
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      fprintf(stdout,
              "usage: %s <input volumes> <output volume>\n",
              Progname) ;
      exit(1) ;
      break ;
    default:
      fprintf(stdout, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}






