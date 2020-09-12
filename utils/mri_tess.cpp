#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/**
 * @brief tesselation routines
 *
 */
/*
 * Original Author: F. Segonne
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MRIio_old.h" /*.H*/
#include "cma.h"       /* using cma_label_to_name() */
#include "const.h"
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "gca.h"
#include "macros.h"
#include "mri.h"
#include "mrisurf.h"
#include "proto.h"

extern const char *Progname;

#ifndef SQR  // BEVIN - the header files might include the one that defines SQR
#define SQR(x) ((x) * (x))
#endif

#define MAXFACES 3000000
#define MAXVERTICES 1500000
#define FACEINCREASE 1.2
#define VERTEXINCREASE 1.2

typedef struct quad_face_type_
{
  int imnr, i, j, f;
  int num;
  int v[4];
} quad_face_type;

typedef struct quad_vertex_type_
{
  int imnr, i, j;
  int num;
  int f[9];
} quad_vertex_type;

typedef struct tesselation_parms_
{
  /*labeled volume*/
  MRI *mri;

  int imgsize, width, height;

  unsigned char ***im;

  /*label information*/
  int number_of_labels;
  int *label_values;
  int *imin, *imax, *jmin, *jmax, *kmin, *kmax;
  int current_label;
  int ind;
  int xmin, xmax, ymin, ymax, zmin, zmax;

  /*final surfaces*/
  /*The surfaces are saved into a table of surfaces*/
  MRIS **mris_table;

  /*all labels are used for the tesselation if(all_flag) */
  int all_flag;

  /*to compute the tesselation consistent with a 6+ topology*/
  int connectivity;

  /*face information*/
  int face_index;
  quad_face_type *face;
  int maxfaces;
  int *face_index_table0;
  int *face_index_table1;

  /*vertex information*/
  int vertex_index;
  quad_vertex_type *vertex;
  int maxvertices;
  //  int *vertex_index_table;

  /*used for the debugging*/
  char fname[100];

} tesselation_parms;

/*definition of the routines*/
static void initTesselationParms(tesselation_parms *parms);
static int freeTesselationParms(tesselation_parms **parms);
static void allocateTesselation(tesselation_parms *parms);
static void freeTesselation(tesselation_parms *parms);

static int saveTesselation(tesselation_parms *parms);
static void add_face(int imnr, int i, int j, int f, int prev_flag, tesselation_parms *parms);
static int add_vertex(int imnr, int i, int j, tesselation_parms *parms);
static int facep(int im0, int i0, int j0, int im1, int i1, int j1, tesselation_parms *parms);
static void check_face(
    int im0, int i0, int j0, int im1, int i1, int j1, int f, int n, int v_ind, int prev_flag, tesselation_parms *parms);
static void make_surface(tesselation_parms *parms);
static void check_face2(
    int im0, int i0, int j0, int im1, int i1, int j1, int f, int n, int v_ind, int prev_flag, tesselation_parms *parms);

#if 0
static void initImage(tesselation_parms *parms)
{
  int i,j,k,a,b,c;
  int nbh[3][3][3];
  double x,y,z,step,val;
  unsigned char ***im=parms->im;
  MRI *mri=parms->mri;
  int label=parms->current_label;
  double threshold;
  int xmin,ymin,zmin,xmax,ymax,zmax;

  if (parms->connectivity==2)
    threshold=0.444;
  else
    threshold=0.22;

  xmin=parms->xmin/3.;
  ymin=parms->ymin/3.;
  zmin=parms->zmin/3.;
  xmax=parms->xmax/3.;
  ymax=parms->ymax/3.;
  zmax=parms->zmax/3.;

  step=1./3.;
  if (parms->all_flag)
  {
    for (k=xmin;k<=xmax;k++)
      for (j=ymin;j<=ymax;j++)
        for (i=xmin;i<=xmax;i++)
        {
          //load nbh;
          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++)
              for (c=-1;c<2;c++)
                nbh[a+1][b+1][c+1]=MRIvox(mri,a+i,j+b,k+c);

          //modify the neighborhood
          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++)
              for (c=-1;c<2;c++)
                if (MRIvox(mri,a+i,j+b,k+c))
                  MRIvox(mri,a+i,j+b,k+c)=1;
                else
                  MRIvox(mri,a+i,j+b,k+c)=0;

          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++)
              for (c=-1;c<2;c++)
              {
                x=i+step*a;
                y=j+step*b;
                z=k+step*c;
                MRIsampleVolume(mri,x,y,z,&val);
                if (val>=threshold)
                  im[3*k+c+1][3*j+b+1][3*i+a+1]=label;
                else
                  im[3*k+c+1][3*j+b+1][3*i+a+1]=0;
              }

          //restitute the correct neighboohod
          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++)
              for (c=-1;c<2;c++)
                MRIvox(mri,a+i,j+b,k+c)=nbh[a+1][b+1][c+1];
        }

  }
  else
  {
    for (k=xmin;k<=xmax;k++)
      for (j=ymin;j<=ymax;j++)
        for (i=xmin;i<=xmax;i++)
        {
          //load nbh;
          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++)
              for (c=-1;c<2;c++)
                nbh[a+1][b+1][c+1]=MRIvox(mri,a+i,j+b,k+c);

          //modify the neighborhood
          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++)
              for (c=-1;c<2;c++)
                if (MRIvox(mri,a+i,j+b,k+c)==label)
                  MRIvox(mri,a+i,j+b,k+c)=1;
                else
                  MRIvox(mri,a+i,j+b,k+c)=0;

          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++)
              for (c=-1;c<2;c++)
              {
                x=i+step*a;
                y=j+step*b;
                z=k+step*c;
                MRIsampleVolume(mri,x,y,z,&val);
                if (val>=threshold)
                  im[3*k+c+1][3*j+b+1][3*i+a+1]=label;
                else
                  im[3*k+c+1][3*j+b+1][3*i+a+1]=0;
              }

          //restitute the correct neighboohod
          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++)
              for (c=-1;c<2;c++)
                MRIvox(mri,a+i,j+b,k+c)=nbh[a+1][b+1][c+1];
        }
  }


}
#endif

/*initialization of tesselation_parms*/
/*note that not all the fields are allocated*/
static void initTesselationParms(tesselation_parms *parms)
{
  int i, j, k, n;
  int depth, width, height;
  int *xmin, *ymin, *zmin, *xmax, *ymax, *zmax;
  int val, label_nbr, *labels;
  unsigned char ***im;
  MRI *mri = parms->mri;

  depth = mri->depth;
  width = mri->width;
  height = mri->height;

  parms->width = 2 * MAX(width, height) + 1;
  parms->height = parms->width;
  parms->imgsize = parms->width * parms->height;

  /* allocate and initialize the image*/
  im = (unsigned char ***)malloc(sizeof(unsigned char **) * depth);
  if (!im) ErrorExit(ERROR_NO_MEMORY, "MRIStesselate: could not allocate the image");
  for (k = 0; k < depth; k++) {
    im[k] = (unsigned char **)malloc(height * sizeof(char *));
    if (!im[k]) ErrorExit(ERROR_NO_MEMORY, "MRIStesselate: could not allocate the image");
    for (j = 0; j < height; j++) {
      im[k][j] = (unsigned char *)calloc(width, sizeof(char));
      if (!im[k][j]) ErrorExit(ERROR_NO_MEMORY, "MRIStesselate: could not allocate the image");
    }
  }
  parms->im = im;

  label_nbr = parms->number_of_labels;
  labels = parms->label_values;

  xmin = (int *)malloc(label_nbr * sizeof(int));
  ymin = (int *)malloc(label_nbr * sizeof(int));
  zmin = (int *)malloc(label_nbr * sizeof(int));
  xmax = (int *)malloc(label_nbr * sizeof(int));
  ymax = (int *)malloc(label_nbr * sizeof(int));
  zmax = (int *)malloc(label_nbr * sizeof(int));

  if ((!xmin) || (!ymin) || (!zmin) || (!xmax) || (!ymax) || (!zmax))
    ErrorExit(ERROR_NO_MEMORY, "MRIStesselate: table of labels");

  for (n = 0; n < label_nbr; n++) {
    xmin[n] = ymin[n] = zmin[n] = 100000;
    xmax[n] = ymax[n] = zmax[n] = 0;
  }

  for (k = 1; k < depth - 1; k++)
    for (j = 1; j < height - 1; j++)
      for (i = 1; i < width - 1; i++) im[k][j][i] = MRIvox(mri, i, j, k);

  for (k = 1; k < mri->depth - 1; k++)
    for (j = 1; j < mri->height - 1; j++)
      for (i = 1; i < mri->width - 1; i++) {
        val = MRIvox(mri, i, j, k);
        for (n = 0; n < label_nbr; n++)
          if (val == labels[n] || (parms->all_flag && val)) {
            if (i < xmin[n]) xmin[n] = i;
            if (j < ymin[n]) ymin[n] = j;
            if (k < zmin[n]) zmin[n] = k;

            if (i > xmax[n]) xmax[n] = i;
            if (j > ymax[n]) ymax[n] = j;
            if (k > zmax[n]) zmax[n] = k;
          }
      }

  /*be careful: the region is expanded by one voxel in the MRI referential*/
  for (n = 0; n < label_nbr; n++) {
    xmin[n]--;
    ymin[n]--;
    zmin[n]--;
    xmax[n]++;
    ymax[n]++;
    zmax[n]++;
  }

  parms->imin = xmin;
  parms->jmin = ymin;
  parms->kmin = zmin;
  parms->imax = xmax;
  parms->jmax = ymax;
  parms->kmax = zmax;

  parms->all_flag = 0;
}

static int freeTesselationParms(tesselation_parms **parms)
{
  int j, k;
  tesselation_parms *tmp = *parms;
  int depth, width, height;
  MRI *mri = tmp->mri;

  depth = mri->depth;
  height = mri->height;
  width = mri->width;

  for (k = 0; k < depth; k++) {
    for (j = 0; j < height; j++) free(tmp->im[k][j]);
    free(tmp->im[k]);
  }
  free(tmp->im);

  free(tmp->imin);
  free(tmp->jmin);
  free(tmp->kmin);
  free(tmp->imax);
  free(tmp->jmax);
  free(tmp->kmax);

  free(*parms);
  return (NO_ERROR);
}

static void allocateTesselation(tesselation_parms *parms)
{
  int imgsize = parms->imgsize;
  parms->face = (quad_face_type *)lcalloc(MAXFACES, sizeof(quad_face_type));
  parms->maxfaces = MAXFACES;
  parms->face_index_table0 = (int *)lcalloc(6 * imgsize, sizeof(int));
  parms->face_index_table1 = (int *)lcalloc(6 * imgsize, sizeof(int));

  parms->vertex = (quad_vertex_type *)lcalloc(MAXVERTICES, sizeof(quad_vertex_type));
  parms->maxvertices = MAXVERTICES;
  //  parms->vertex_index_table = (int *)lcalloc(8*imgsize,sizeof(int));

  parms->face_index = 0;
  parms->vertex_index = 0;

  if ((!parms->face) || (!parms->face_index_table0) || (!parms->face_index_table1) || (!parms->vertex))
    //     || (!parms->vertex_index_table))
    ErrorExit(ERROR_NO_MEMORY, "MRIStesselate: local tesselation tables");
}

static void freeTesselation(tesselation_parms *parms)
{
  free(parms->face);
  free(parms->face_index_table0);
  free(parms->face_index_table1);

  free(parms->vertex);
  //  free(parms->vertex_index_table);
}

static int saveTesselation(tesselation_parms *parms)
{
  int vno, m, n, fno;
  quad_face_type *face2;
  FACE *face;
  float x, y, z, xhi, xlo, yhi, ylo, zhi, zlo;
  float st, ps, xx1, yy0, zz1;
  float j, i, imnr;
  double xw, yw, zw;

  /*necessary for the coord transformation*/
  ps = parms->mri->ps;
  st = parms->mri->thick;

  yy0 = parms->mri->ystart;
  xx1 = parms->mri->xend;
  zz1 = parms->mri->zend;

  MRIS *mris = MRISoverAlloc(
                    parms->vertex_index, 2 * parms->face_index, 
                    parms->vertex_index, 2 * parms->face_index);
  mris->type = MRIS_BINARY_QUADRANGLE_FILE;

  /*first init vertices*/
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY * const vertex_topology = &mris->vertices_topology[vno];

    quad_vertex_type const * const vertex2 = &parms->vertex[vno];
    
    i = vertex2->i;
    j = vertex2->j;
    imnr = vertex2->imnr;

    if (parms->connectivity >= 1) {
      i = (i + 0.5) / 2.;
      j = (j + 0.5) / 2.;
      imnr = (imnr + 0.5) / 2.;
    }

    //      x = xx1-(j-0.5)*ps;
    // y = yy0+(imnr-0.5)*st;
    // z = zz1-(i-0.5)*ps;

    MRIvoxelToWorld(parms->mri, j - 0.5, i - 0.5, imnr - 0.5, &xw, &yw, &zw);
    x = xw;
    y = yw;
    z = zw;

    MRISsetXYZ(mris,vno, x, y, z);
    vertex_topology->num = 0;
  }
  
  /*then init faces by dividing each quadrant into two triangles*/
  for (m = 0; m < parms->face_index; m++) {
    int which;

    fno = 2 * m;
    face2 = &parms->face[m];

    /* if we're going to be arbitrary, we might as well be really arbitrary */
    /*
       NOTE: for this to work properly in the write, the first two
       vertices in the first face (EVEN and ODD) must be 0 and 1.
    */
    which = WHICH_FACE_SPLIT(face2->v[0], face2->v[1]);
    if (EVEN(which)) {
      /* 1st triangle */
      mris->faces[fno].v[0] = face2->v[0];
      mris->faces[fno].v[1] = face2->v[1];
      mris->faces[fno].v[2] = face2->v[3];

      /* 2nd triangle */
      mris->faces[fno + 1].v[0] = face2->v[2];
      mris->faces[fno + 1].v[1] = face2->v[3];
      mris->faces[fno + 1].v[2] = face2->v[1];
    }
    else {
      /* 1st triangle */
      mris->faces[fno].v[0] = face2->v[0];
      mris->faces[fno].v[1] = face2->v[1];
      mris->faces[fno].v[2] = face2->v[2];

      /* 2nd triangle */
      mris->faces[fno + 1].v[0] = face2->v[0];
      mris->faces[fno + 1].v[1] = face2->v[2];
      mris->faces[fno + 1].v[2] = face2->v[3];
    }
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      mris->vertices_topology[mris->faces[fno    ].v[n]].num++;
      mris->vertices_topology[mris->faces[fno + 1].v[n]].num++;
    }
  }
  /*allocate indices & faces for each vertex*/
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices_topology[vno].f = (int *)calloc(mris->vertices_topology[vno].num, sizeof(int));
    if (!mris->vertices_topology[vno].f)
      ErrorExit(
          ERROR_NOMEMORY, "%s: could not allocate %d faces at %dth vertex", Progname, mris->vertices_topology[vno].num, vno);
    mris->vertices_topology[vno].n = (uchar *)calloc(mris->vertices_topology[vno].num, sizeof(uchar));
    if (!mris->vertices_topology[vno].n)
      ErrorExit(
          ERROR_NOMEMORY, "%s: could not allocate %d indices at %dth vertex", Progname, mris->vertices_topology[vno].num, vno);
    mris->vertices_topology[vno].num = 0;
  }
  /*init faces for each vertex*/
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    for (n = 0; n < VERTICES_PER_FACE; n++) 
    	mris->vertices_topology[face->v[n]].f[mris->vertices_topology[face->v[n]].num++] = fno;
  }

  mrisCheckVertexFaceTopology(mris);

  /*necessary initialization*/
  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].curv = 0;
    mris->vertices[vno].origarea = -1;
    mris->vertices[vno].border = 0;

    for (n = 0; n < mris->vertices_topology[vno].num; n++) {
      for (m = 0; m < VERTICES_PER_FACE; m++) {
        if (mris->faces[mris->vertices_topology[vno].f[n]].v[m] == vno) mris->vertices_topology[vno].n[n] = m;
      }
    }
    x = mris->vertices[vno].x;
    y = mris->vertices[vno].y;
    z = mris->vertices[vno].z;
    if (x > xhi) xhi = x;
    if (x < xlo) xlo = x;
    if (y > yhi) yhi = y;
    if (y < ylo) ylo = y;
    if (z > zhi) zhi = z;
    if (z < zlo) zlo = z;
  }
  mris->xlo = xlo;
  mris->ylo = ylo;
  mris->zlo = zlo;
  mris->xhi = xhi;
  mris->yhi = yhi;
  mris->zhi = zhi;
  mris->xctr = (xhi + xlo) / 2;
  mris->yctr = (yhi + ylo) / 2;
  mris->zctr = (zhi + zlo) / 2;

  mrisCompleteTopology(mris);
  MRIScomputeNormals(mris);

  mris->type = MRIS_TRIANGULAR_SURFACE; /*not so sure about that*/
  MRISsetNeighborhoodSizeAndDist(mris, 2);
  MRIScomputeSecondFundamentalForm(mris);
  MRISuseMeanCurvature(mris);
  mris->radius = MRISaverageRadius(mris);
  MRIScomputeMetricProperties(mris);
  MRISstoreCurrentPositions(mris);
  parms->mris_table[parms->ind] = mris;
  return (NO_ERROR);
}

MRIS *MRISconcatenateQuadSurfaces(int number_of_labels, MRIS **mris_tab)
{
  MRIS *mris;
  int nvertices, nfaces, n, count, countvnbr, vno, fno;

  for (n = 0, nvertices = 0, nfaces = 0; n < number_of_labels; n++) {
    nvertices += mris_tab[n]->nvertices;
    nfaces += mris_tab[n]->nfaces;
  }

  mris = MRISalloc(nvertices, nfaces);
  mris->type = MRIS_BINARY_QUADRANGLE_FILE;

  for (n = 0, count = 0; n < number_of_labels; n++)
    for (vno = 0; vno < mris_tab[n]->nvertices; vno++, count++) {
      MRISsetXYZ(mris,count, 
        mris_tab[n]->vertices[vno].x,
        mris_tab[n]->vertices[vno].y,
        mris_tab[n]->vertices[vno].z);
      mris->vertices[count].curv = mris_tab[n]->vertices[vno].curv;
      mris->vertices[count].val = (float)(n - number_of_labels / 2.) / number_of_labels;
    }
  for (n = 0, count = 0, countvnbr = 0; n < number_of_labels; n++) {
    for (fno = 0; fno < mris_tab[n]->nfaces; fno++, count++) {
      mris->faces[count].v[0] = countvnbr + mris_tab[n]->faces[fno].v[0];
      mris->faces[count].v[1] = countvnbr + mris_tab[n]->faces[fno].v[1];
      mris->faces[count].v[2] = countvnbr + mris_tab[n]->faces[fno].v[2];
    }
    countvnbr += mris_tab[n]->nvertices;
  }

  mrisCheckVertexFaceTopology(mris);

  return mris;
}

static void reallocateVertices(tesselation_parms *parms)
{
  quad_vertex_type *tmp;
  int k, n;
  int vertex_index = parms->vertex_index;

  parms->maxvertices = (int)(parms->maxvertices * VERTEXINCREASE);
  tmp = (quad_vertex_type *)lcalloc(parms->maxvertices, sizeof(quad_vertex_type));
  if (!tmp) ErrorExit(ERROR_NOMEMORY, "%s: max vertices %d exceeded", Progname, parms->maxvertices);
  for (k = 0; k < vertex_index; k++) {
    tmp[k].imnr = parms->vertex[k].imnr;
    tmp[k].i = parms->vertex[k].i;
    tmp[k].j = parms->vertex[k].j;
    tmp[k].num = parms->vertex[k].num;
    for (n = 0; n < 9; n++) tmp[k].f[n] = parms->vertex[k].f[n];
  }
  free(parms->vertex);
  parms->vertex = tmp;
}
static void reallocateFaces(tesselation_parms *parms)
{
  quad_face_type *tmp;
  int k, n;
  int face_index = parms->face_index;

  parms->maxfaces = (int)(parms->maxfaces * FACEINCREASE);
  tmp = (quad_face_type *)lcalloc(parms->maxfaces, sizeof(quad_face_type));
  if (!tmp) ErrorExit(ERROR_NOMEMORY, "%s: max faces %d exceeded", Progname, parms->maxfaces);
  for (k = 0; k < face_index; k++) {
    tmp[k].imnr = parms->face[k].imnr;
    tmp[k].i = parms->face[k].i;
    tmp[k].j = parms->face[k].j;
    tmp[k].f = parms->face[k].f;
    tmp[k].num = parms->face[k].num;
    for (n = 0; n < 4; n++) tmp[k].v[n] = parms->face[k].v[n];
  }

  free(parms->face);
  parms->face = tmp;
}

static void add_face(int imnr, int i, int j, int f, int prev_flag, tesselation_parms *parms)
{
  int pack = f * parms->imgsize + i * parms->width + j;
  int face_index = parms->face_index;

  if (face_index >= parms->maxfaces - 1) reallocateFaces(parms);
  if (prev_flag)
    parms->face_index_table0[pack] = face_index;
  else
    parms->face_index_table1[pack] = face_index;
  parms->face[face_index].imnr = imnr;
  parms->face[face_index].i = i;
  parms->face[face_index].j = j;
  parms->face[face_index].f = f;
  parms->face[face_index].num = 0;
  parms->face_index++;
}

static int add_vertex(int imnr, int i, int j, tesselation_parms *parms)
{
  //  int pack = i*(parms->width+1)+j;
  int vertex_index = parms->vertex_index;

  if (vertex_index >= parms->maxvertices - 1) reallocateVertices(parms);
  //  parms->vertex_index_table[pack] = vertex_index;
  parms->vertex[vertex_index].imnr = imnr;
  parms->vertex[vertex_index].i = i;
  parms->vertex[vertex_index].j = j;
  parms->vertex[vertex_index].num = 0;
  return parms->vertex_index++;
}

static int facep(int im0, int i0, int j0, int im1, int i1, int j1, tesselation_parms *parms)
{
  unsigned char ***im = parms->im;
  int value = parms->current_label;
  return (im0 >= parms->zmin && im0 <= parms->zmax && i0 >= parms->ymin && i0 <= parms->ymax && j0 >= parms->xmin &&
          j0 <= parms->xmax && im1 >= parms->zmin && im1 <= parms->zmax && i1 >= parms->ymin && i1 <= parms->ymax &&
          j1 >= parms->xmin && j1 <= parms->xmax && im[im0][i0][j0] != im[im1][i1][j1] &&
          ((im[im0][i0][j0] == value || im[im1][i1][j1] == value) || parms->all_flag));
}

static void check_face(
    int im0, int i0, int j0, int im1, int i1, int j1, int f, int n, int v_ind, int prev_flag, tesselation_parms *parms)
{
  int f_pack = f * parms->imgsize + i0 * parms->width + j0;
  int f_ind;
  unsigned char ***im = parms->im;
  int value = parms->current_label;

  if ((im0 >= parms->zmin && im0 <= parms->zmax && i0 >= parms->ymin && i0 <= parms->ymax && j0 >= parms->xmin &&
       j0 <= parms->xmax && im1 >= parms->zmin && im1 <= parms->zmax && i1 >= parms->ymin && i1 <= parms->ymax &&
       j1 >= parms->xmin && j1 <= parms->xmax)) {
    if ((parms->all_flag && ((im[im0][i0][j0] != im[im1][i1][j1]) && (im[im1][i1][j1] == 0))) ||
        (((im[im0][i0][j0] == value) && (im[im1][i1][j1] != value)))) {
      if (n == 0) {
        add_face(im0, i0, j0, f, prev_flag, parms);
      }
      if (prev_flag)
        f_ind = parms->face_index_table0[f_pack];
      else
        f_ind = parms->face_index_table1[f_pack];
      parms->face[f_ind].v[n] = v_ind;
      if (parms->vertex[v_ind].num < 9) parms->vertex[v_ind].f[parms->vertex[v_ind].num++] = f_ind;
    }
  }
}

static void check_face2(
    int im0, int i0, int j0, int im1, int i1, int j1, int f, int n, int v_ind, int prev_flag, tesselation_parms *parms)
{
  int f_pack = f * parms->imgsize + i0 * parms->width + j0;
  int f_ind;

  if (n == 0) {
    add_face(im0, i0, j0, f, prev_flag, parms);
  }
  if (prev_flag)
    f_ind = parms->face_index_table0[f_pack];
  else
    f_ind = parms->face_index_table1[f_pack];
  parms->face[f_ind].v[n] = v_ind;
  if (parms->vertex[v_ind].num < 9) parms->vertex[v_ind].f[parms->vertex[v_ind].num++] = f_ind;
}

static void make_surface(tesselation_parms *parms)
{
  int imnr, i, j, f_pack, v_ind, f;
  int n;

  for (n = 0; n < parms->number_of_labels; n++) {
    parms->ind = n;
    parms->current_label = parms->label_values[parms->ind];
    parms->xmin = parms->imin[parms->ind];
    parms->ymin = parms->jmin[parms->ind];
    parms->zmin = parms->kmin[parms->ind];
    parms->xmax = parms->imax[parms->ind];
    parms->ymax = parms->jmax[parms->ind];
    parms->zmax = parms->kmax[parms->ind];

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "working on label %d (%s)\n", parms->current_label, cma_label_to_name(parms->current_label));

    allocateTesselation(parms);

    for (imnr = parms->zmin; imnr <= parms->zmax; imnr++) {
      if ((Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) && (parms->vertex_index || parms->face_index) && !(imnr % 10))
        fprintf(stdout, "slice %d: %d vertices, %d faces\n", imnr, parms->vertex_index, parms->face_index);
      for (i = parms->ymin; i <= parms->ymax; i++)
        for (j = parms->xmin; j <= parms->xmax; j++) {
          if (facep(imnr, i - 1, j - 1, imnr - 1, i - 1, j - 1, parms) ||
              facep(imnr, i - 1, j, imnr - 1, i - 1, j, parms) || facep(imnr, i, j, imnr - 1, i, j, parms) ||
              facep(imnr, i, j - 1, imnr - 1, i, j - 1, parms) ||
              facep(imnr - 1, i, j - 1, imnr - 1, i - 1, j - 1, parms) ||
              facep(imnr - 1, i, j, imnr - 1, i - 1, j, parms) || facep(imnr, i, j, imnr, i - 1, j, parms) ||
              facep(imnr, i, j - 1, imnr, i - 1, j - 1, parms) ||
              facep(imnr - 1, i - 1, j, imnr - 1, i - 1, j - 1, parms) ||
              facep(imnr - 1, i, j, imnr - 1, i, j - 1, parms) || facep(imnr, i, j, imnr, i, j - 1, parms) ||
              facep(imnr, i - 1, j, imnr, i - 1, j - 1, parms)) {
            v_ind = add_vertex(imnr, i, j, parms);
            check_face(imnr, i - 1, j - 1, imnr - 1, i - 1, j - 1, 0, 2, v_ind, 0, parms);
            check_face(imnr, i - 1, j, imnr - 1, i - 1, j, 0, 3, v_ind, 0, parms);
            check_face(imnr, i, j, imnr - 1, i, j, 0, 0, v_ind, 0, parms);
            check_face(imnr, i, j - 1, imnr - 1, i, j - 1, 0, 1, v_ind, 0, parms);
            check_face(imnr - 1, i, j - 1, imnr - 1, i - 1, j - 1, 2, 2, v_ind, 1, parms);
            check_face(imnr - 1, i, j, imnr - 1, i - 1, j, 2, 1, v_ind, 1, parms);
            check_face(imnr, i, j, imnr, i - 1, j, 2, 0, v_ind, 0, parms);
            check_face(imnr, i, j - 1, imnr, i - 1, j - 1, 2, 3, v_ind, 0, parms);
            check_face(imnr - 1, i - 1, j, imnr - 1, i - 1, j - 1, 4, 2, v_ind, 1, parms);
            check_face(imnr - 1, i, j, imnr - 1, i, j - 1, 4, 3, v_ind, 1, parms);
            check_face(imnr, i, j, imnr, i, j - 1, 4, 0, v_ind, 0, parms);
            check_face(imnr, i - 1, j, imnr, i - 1, j - 1, 4, 1, v_ind, 0, parms);

            check_face(imnr - 1, i - 1, j - 1, imnr, i - 1, j - 1, 1, 2, v_ind, 1, parms);
            check_face(imnr - 1, i - 1, j, imnr, i - 1, j, 1, 1, v_ind, 1, parms);
            check_face(imnr - 1, i, j, imnr, i, j, 1, 0, v_ind, 1, parms);
            check_face(imnr - 1, i, j - 1, imnr, i, j - 1, 1, 3, v_ind, 1, parms);
            check_face(imnr - 1, i - 1, j - 1, imnr - 1, i, j - 1, 3, 2, v_ind, 1, parms);
            check_face(imnr - 1, i - 1, j, imnr - 1, i, j, 3, 3, v_ind, 1, parms);
            check_face(imnr, i - 1, j, imnr, i, j, 3, 0, v_ind, 0, parms);
            check_face(imnr, i - 1, j - 1, imnr, i, j - 1, 3, 1, v_ind, 0, parms);
            check_face(imnr - 1, i - 1, j - 1, imnr - 1, i - 1, j, 5, 2, v_ind, 1, parms);
            check_face(imnr - 1, i, j - 1, imnr - 1, i, j, 5, 1, v_ind, 1, parms);
            check_face(imnr, i, j - 1, imnr, i, j, 5, 0, v_ind, 0, parms);
            check_face(imnr, i - 1, j - 1, imnr, i - 1, j, 5, 3, v_ind, 0, parms);
          }
        }

      for (i = 0; i < parms->width; i++)
        for (j = 0; j < parms->height; j++)
          for (f = 0; f < 6; f++) {
            f_pack = f * parms->imgsize + i * parms->width + j;
            parms->face_index_table0[f_pack] = parms->face_index_table1[f_pack];
          }
    }
    saveTesselation(parms);
    freeTesselation(parms);
  }
}

typedef unsigned char CNBH[2][2][2];

static int computeconnectedcomponents(CNBH *tab)
{
  int i, j, k, a, b, c, ik, jk, kk, ct;
  int nvox, label, sum;
  int comp_table[12];
  int min_val;
  int x, y, z;

  memset(comp_table, 0, sizeof(comp_table));

  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        if ((*tab)[i][j][k]) {
          for (nvox = 0, ik = -1; ik <= 1; ik++) {
            a = i + ik;
            if (a < 0 || a >= 2) continue;
            for (jk = -1; jk <= 1; jk++) {
              b = j + jk;
              if (b < 0 || b >= 2) continue;
              for (kk = -1; kk <= 1; kk++) {
                sum = abs(ik) + abs(jk) + abs(kk);
                if (sum > 1 || (!sum)) continue;
                c = k + kk;
                if (c < 0 || c >= 2) continue;
                label = (*tab)[a][b][c];
                if (label > 1) {
                  comp_table[label - 1] = 2;
                  nvox++;
                }
              }
            }
          }
          if (!nvox)  // find new basin!
          {
            for (ct = 1; comp_table[ct] && ct < 10; ct++)
              ;
            (*tab)[i][j][k] = ct + 1;  // label the new basin
            comp_table[ct] = 1;        // note that this number is taken
          }
          else {
            min_val = 11;

            // merging into the smallest value
            for (ct = 1; ct < 10; ct++)
              if (comp_table[ct] == 2) {
                min_val = ct;
                break;
              }

            (*tab)[i][j][k] = min_val + 1;
            comp_table[min_val] = 1;

            // merging of the other neighboring values into the smallest one
            for (ct = min_val + 1; ct < 10; ct++)
              if (comp_table[ct] == 2) {
                for (x = 0; x < 2; x++)
                  for (y = 0; y < 2; y++)
                    for (z = 0; z < 2; z++)
                      if ((*tab)[x][y][z] == ct + 1) (*tab)[x][y][z] = min_val + 1;
                // specify that this basin nbr
                // doesn't exist anymore
                comp_table[ct] = 0;
              }
          }
        }

  for (nvox = 0, ct = 1; ct < 10; ct++)
    if (comp_table[ct]) nvox++;

  return nvox;
}

static void make_surface_with_connectivity(tesselation_parms *parms)
{
  int imnr, i, j, k, f_pack, v_ind, f, a, b, c, label, nlabels, comp, m;
  int n, slcx, slcy, slcz, indx, indy, indz;
  CNBH tab;
  double threshold, val;
  double x, y, z, step;
  MRI *mri, *mri_src = parms->mri;

  switch (parms->connectivity) {
    case 1:
      threshold = 0.510000f;
      break;
    case 2:
      threshold = 0.26000f;
      break;
    case 3:
      threshold = 0.751f;
      break;
    case 4:
      threshold = 0.25f;
      break;
    default:
      threshold = 0.5f;
      break;
  };
  step = 0.5f;

  mri = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_UCHAR);
  for (n = 0; n < parms->number_of_labels; n++) {
    parms->ind = n;
    parms->current_label = parms->label_values[parms->ind];
    parms->xmin = parms->imin[parms->ind];
    parms->ymin = parms->jmin[parms->ind];
    parms->zmin = parms->kmin[parms->ind];
    parms->xmax = parms->imax[parms->ind];
    parms->ymax = parms->jmax[parms->ind];
    parms->zmax = parms->kmax[parms->ind];

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "working on label %d (%s)\n", parms->current_label, cma_label_to_name(parms->current_label));

    // init temporary MRI
    if (parms->all_flag)
      for (k = parms->zmin; k <= parms->zmax; k++)
        for (j = parms->ymin; j <= parms->ymax; j++)
          for (i = parms->xmin; i <= parms->xmax; i++)
            if (MRIvox(mri_src, i, j, k))
              MRIvox(mri, i, j, k) = 1;
            else
              MRIvox(mri, i, j, k) = 0;
    else
      for (k = parms->zmin; k <= parms->zmax; k++)
        for (j = parms->ymin; j <= parms->ymax; j++)
          for (i = parms->xmin; i <= parms->xmax; i++)
            if (MRIvox(mri_src, i, j, k) == parms->current_label)
              MRIvox(mri, i, j, k) = 1;
            else
              MRIvox(mri, i, j, k) = 0;

    allocateTesselation(parms);
    label = parms->current_label;

    for (imnr = parms->zmin; imnr <= parms->zmax + 1; imnr++)
      for (slcz = -1; slcz < 1; slcz++) {
        if ((Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) && (parms->vertex_index || parms->face_index) && !(imnr % 10))
          fprintf(stdout, "slice %d: %d vertices, %d faces\n", imnr, parms->vertex_index, parms->face_index);
        for (i = parms->ymin; i <= parms->ymax + 1; i++)
          for (j = parms->xmin; j <= parms->xmax + 1; j++)
            for (slcy = -1; slcy < 1; slcy++)
              for (slcx = -1; slcx < 1; slcx++) {
                indx = 2 * j + slcx + 1;
                indy = 2 * i + slcy + 1;
                indz = 2 * imnr + slcz + 1;

                // load the 2by2by2 cube and find the connected component
                for (a = 0; a < 2; a++)
                  for (b = 0; b < 2; b++)
                    for (c = 0; c < 2; c++) {
                      x = j + slcx * step + step * a;
                      y = i + slcy * step + step * b;
                      z = imnr + slcz * step + step * c;
                      MRIsampleVolume(mri, x, y, z, &val);
                      if (val >= threshold)
                        tab[c][b][a] = 1;
                      else
                        tab[c][b][a] = 0;
                    }

                // find connected components
                nlabels = computeconnectedcomponents(&tab);

                for (m = 0; m < nlabels; m++) {
                  if ((  tab[0][0][0] * tab[0][0][1] * tab[0][1][0]
		       * tab[0][1][1] * tab[1][0][0] * tab[1][0][1]
		       * tab[1][1][0] * tab[1][1][1])
		      != 0) {
                    continue;
		  }
                  comp = m + 2;

                  /*if(!(((tab[1][0][0]==comp)&&(!tab[0][0][0]))  ||
                     ((tab[1][0][1]==comp)&&(!tab[0][0][1]))  ||
                     ((tab[1][1][1]==comp)&&(!tab[0][1][1]))  ||
                     ((tab[1][1][0]==comp)&&(!tab[0][1][0]))  ||
                     ((tab[0][1][0]==comp)&&(!tab[0][0][0]))  ||
                     ((tab[0][1][1]==comp)&&(!tab[0][0][1]))  ||
                     ((tab[1][1][1]==comp)&&(!tab[1][0][1]))  ||
                     ((tab[1][1][0]==comp)&&(!tab[1][0][0]))  ||
                     ((tab[0][0][1]==comp)&&(!tab[0][0][0]))  ||
                     ((tab[0][1][1]==comp)&&(!tab[0][1][0]))  ||
                     ((tab[1][1][1]==comp)&&(!tab[1][1][0]))  ||
                  ((tab[1][0][1]==comp)&&(!tab[1][0][0]))))
                  continue;*/

                  // for each component, we have to add a vertex
                  v_ind = add_vertex(indz, indy, indx, parms);
                  // fprintf(stderr,"\n lab %d  Nbr = %d ",nlabels,parms->vertex_index);

                  // then check which faces are concerned
                  if ((tab[1][0][0] == comp) && (!tab[0][0][0]))
                    check_face2(indz, indy - 1, indx - 1, indz - 1, indy - 1, indx - 1, 0, 2, v_ind, 0, parms);
                  if ((tab[1][0][1] == comp) && (!tab[0][0][1]))
                    check_face2(indz, indy - 1, indx, indz - 1, indy - 1, indx, 0, 3, v_ind, 0, parms);
                  if ((tab[1][1][1] == comp) && (!tab[0][1][1]))
                    check_face2(indz, indy, indx, indz - 1, indy, indx, 0, 0, v_ind, 0, parms);
                  if ((tab[1][1][0] == comp) && (!tab[0][1][0]))
                    check_face2(indz, indy, indx - 1, indz - 1, indy, indx - 1, 0, 1, v_ind, 0, parms);
                  if ((tab[0][1][0] == comp) && (!tab[0][0][0]))
                    check_face2(indz - 1, indy, indx - 1, indz - 1, indy - 1, indx - 1, 2, 2, v_ind, 1, parms);
                  if ((tab[0][1][1] == comp) && (!tab[0][0][1]))
                    check_face2(indz - 1, indy, indx, indz - 1, indy - 1, indx, 2, 1, v_ind, 1, parms);
                  if ((tab[1][1][1] == comp) && (!tab[1][0][1]))
                    check_face2(indz, indy, indx, indz, indy - 1, indx, 2, 0, v_ind, 0, parms);
                  if ((tab[1][1][0] == comp) && (!tab[1][0][0]))
                    check_face2(indz, indy, indx - 1, indz, indy - 1, indx - 1, 2, 3, v_ind, 0, parms);
                  if ((tab[0][0][1] == comp) && (!tab[0][0][0]))
                    check_face2(indz - 1, indy - 1, indx, indz - 1, indy - 1, indx - 1, 4, 2, v_ind, 1, parms);
                  if ((tab[0][1][1] == comp) && (!tab[0][1][0]))
                    check_face2(indz - 1, indy, indx, indz - 1, indy, indx - 1, 4, 3, v_ind, 1, parms);
                  if ((tab[1][1][1] == comp) && (!tab[1][1][0]))
                    check_face2(indz, indy, indx, indz, indy, indx - 1, 4, 0, v_ind, 0, parms);
                  if ((tab[1][0][1] == comp) && (!tab[1][0][0]))
                    check_face2(indz, indy - 1, indx, indz, indy - 1, indx - 1, 4, 1, v_ind, 0, parms);

                  if ((tab[0][0][0] == comp) && (!tab[1][0][0]))
                    check_face2(indz - 1, indy - 1, indx - 1, indz, indy - 1, indx - 1, 1, 2, v_ind, 1, parms);
                  if ((tab[0][0][1] == comp) && (!tab[1][0][1]))
                    check_face2(indz - 1, indy - 1, indx, indz, indy - 1, indx, 1, 1, v_ind, 1, parms);
                  if ((tab[0][1][1] == comp) && (!tab[1][1][1]))
                    check_face2(indz - 1, indy, indx, indz, indy, indx, 1, 0, v_ind, 1, parms);
                  if ((tab[0][1][0] == comp) && (!tab[1][1][0]))
                    check_face2(indz - 1, indy, indx - 1, indz, indy, indx - 1, 1, 3, v_ind, 1, parms);
                  if ((tab[0][0][0] == comp) && (!tab[0][1][0]))
                    check_face2(indz - 1, indy - 1, indx - 1, indz - 1, indy, indx - 1, 3, 2, v_ind, 1, parms);
                  if ((tab[0][0][1] == comp) && (!tab[0][1][1]))
                    check_face2(indz - 1, indy - 1, indx, indz - 1, indy, indx, 3, 3, v_ind, 1, parms);
                  if ((tab[1][0][1] == comp) && (!tab[1][1][1]))
                    check_face2(indz, indy - 1, indx, indz, indy, indx, 3, 0, v_ind, 0, parms);
                  if ((tab[1][0][0] == comp) && (!tab[1][1][0]))
                    check_face2(indz, indy - 1, indx - 1, indz, indy, indx - 1, 3, 1, v_ind, 0, parms);
                  if ((tab[0][0][0] == comp) && (!tab[0][0][1]))
                    check_face2(indz - 1, indy - 1, indx - 1, indz - 1, indy - 1, indx, 5, 2, v_ind, 1, parms);
                  if ((tab[0][1][0] == comp) && (!tab[0][1][1]))
                    check_face2(indz - 1, indy, indx - 1, indz - 1, indy, indx, 5, 1, v_ind, 1, parms);
                  if ((tab[1][1][0] == comp) && (!tab[1][1][1]))
                    check_face2(indz, indy, indx - 1, indz, indy, indx, 5, 0, v_ind, 0, parms);
                  if ((tab[1][0][0] == comp) && (!tab[1][0][1]))
                    check_face2(indz, indy - 1, indx - 1, indz, indy - 1, indx, 5, 3, v_ind, 0, parms);
                }
              }
        //     for (i=parms->xmin;i<=parms->xmax;i++)
        // for (j=parms->ymin;j<=parms->ymax;j++)
        // for(slcy=-1;slcy<1;slcy++)
        //   for(slcx=-1;slcx<1;slcx++)
        for (i = 0; i < parms->width; i++)
          for (j = 0; j < parms->height; j++)
            for (f = 0; f < 6; f++) {
              // indx=2*j+slcx+1; indy=2*i+slcy+1;
              f_pack = f * parms->imgsize + i * parms->width + j;
              parms->face_index_table0[f_pack] = parms->face_index_table1[f_pack];
            }
      }

    saveTesselation(parms);
    freeTesselation(parms);
    for (k = parms->zmin; k <= parms->zmax; k++)
      for (j = parms->ymin; j <= parms->ymax; j++)
        for (i = parms->xmin; i <= parms->xmax; i++) MRIvox(mri, i, j, k) = 0;
  }
  MRIfree(&mri);
}

MRIS *MRIScreateSurfaceFromVolume(MRI *mri, int label, int connectivity)
{
  tesselation_parms *parms;
  MRIS **mris_table, *mris;

  parms = (tesselation_parms *)calloc(1, sizeof(tesselation_parms));
  if (!parms) ErrorExit(ERROR_NOMEMORY, "tesselation parms\n");
  parms->mri = mri;

  /*init tesselation_parms structure*/
  parms->number_of_labels = 1;  // only one single label
  parms->label_values = (int *)malloc(sizeof(int));
  parms->label_values[0] = label;
  parms->ind = 0;
  mris_table = (MRIS **)malloc(sizeof(MRIS *));  // final surface information
  parms->mris_table = mris_table;
  if ((!parms->label_values) || (!mris_table)) ErrorExit(ERROR_NOMEMORY, "labels/surfaces tables\n");

  parms->connectivity = connectivity;
  if (!label)  // if label==0 then all_flag mode
    parms->all_flag = 1;
  initTesselationParms(parms);  // final parms init

  if (parms->connectivity > 0)
    make_surface_with_connectivity(parms);
  else
    make_surface(parms);

  free(parms->label_values);
  mris = parms->mris_table[0];
  free(parms->mris_table);
  freeTesselationParms(&parms);

  return mris;
}

MRIS **MRIScreateSurfacesFromVolume(MRI *mri, int number_of_labels, int *labelvalues, int connectivity)
{
  tesselation_parms *parms;
  MRIS **mris_table;

  if (!labelvalues) ErrorExit(ERROR_BADPARM, "int *labelvalues\n");

  parms = (tesselation_parms *)calloc(1, sizeof(tesselation_parms));
  parms->mri = mri;
  if (!parms) ErrorExit(ERROR_NOMEMORY, "tesselation parms\n");

  /*init tesselation_parms structure*/
  parms->number_of_labels = number_of_labels;  // label information
  parms->label_values = labelvalues;
  mris_table = (MRIS **)malloc(number_of_labels * sizeof(MRIS *));  // final surfaces information
  if ((!parms->label_values) || (!mris_table)) ErrorExit(ERROR_NOMEMORY, "surfaces tables\n");

  parms->mris_table = mris_table;
  initTesselationParms(parms);  // final parms init

  parms->connectivity = connectivity;
  if (parms->connectivity > 0)
    make_surface_with_connectivity(parms);
  else
    make_surface(parms);

  freeTesselationParms(&parms);

  return mris_table;
}
