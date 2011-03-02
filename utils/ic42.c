/**
 * @file  ic42.c
 * @brief 1st order icosahedral subdivision with 42 vertices.
 *
 * in memory version of ic1.tri
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.2 $
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


#include <stdio.h>
#include "mrisurf.h"
#include "error.h"
#include "diag.h"
#include "icosahedron.h"


#define ICO_NVERTICES 42


IC_VERTEX ic42_vertices[42] =
  {
    {   .0000,   .0000,  1.0000},
    {   .8944,   .0000,   .4472},
    {   .2764,   .8507,   .4472},
    {  -.7236,   .5257,   .4472},
    {  -.7236,  -.5257,   .4472},
    {   .2764,  -.8507,   .4472},
    {   .7236,  -.5257,  -.4472},
    {   .7236,   .5257,  -.4472},
    {  -.2764,   .8507,  -.4472},
    {  -.8944,   .0000,  -.4472},
    {  -.2764,  -.8507,  -.4472},
    {   .0000,   .0000, -1.0000},
    {  -.4253,  -.3090,   .8507},
    {  -.8507,   .0000,   .5257},
    {  -.4253,   .3090,   .8507},
    {   .1625,  -.5000,   .8507},
    {  -.2629,  -.8090,   .5257},
    {   .5257,   .0000,   .8507},
    {   .6882,  -.5000,   .5257},
    {   .1625,   .5000,   .8507},
    {   .6882,   .5000,   .5257},
    {  -.2629,   .8090,   .5257},
    {  -.5878,   .8090,   .0000},
    {   .0000,  1.0000,   .0000},
    {  -.9511,   .3090,   .0000},
    {  -.6882,   .5000,  -.5257},
    {  -.9511,  -.3090,   .0000},
    {  -.5878,  -.8090,   .0000},
    {  -.6882,  -.5000,  -.5257},
    {   .0000, -1.0000,   .0000},
    {   .5878,  -.8090,   .0000},
    {   .2629,  -.8090,  -.5257},
    {   .9511,  -.3090,   .0000},
    {   .9511,   .3090,   .0000},
    {   .8507,   .0000,  -.5257},
    {   .5878,   .8090,   .0000},
    {   .2629,   .8090,  -.5257},
    {  -.5257,   .0000,  -.8507},
    {  -.1625,   .5000,  -.8507},
    {  -.1625,  -.5000,  -.8507},
    {   .4253,  -.3090,  -.8507},
    {   .4253,   .3090,  -.8507}
  } ;

#define ICO_NFACES  80

IC_FACE  ic42_faces[ICO_NFACES] =
  {
    {{ 1,  13,  15 }},
    {{ 13,   5,  14 }},
    {{ 15,  13,  14 }},
    {{ 15,  14,   4 }},
    {{ 1,  16,  13 }},
    {{ 16,   6,  17 }},
    {{ 13,  16,  17 }},
    {{ 13,  17,   5 }},
    {{ 1,  18,  16 }},
    {{ 18,   2,  19 }},
    {{ 16,  18,  19 }},
    {{ 16,  19,   6 }},
    {{ 1,  20,  18 }},
    {{ 20,   3,  21 }},
    {{ 18,  20,  21 }},
    {{ 18,  21,   2 }},
    {{ 1,  15,  20 }},
    {{ 15,   4,  22 }},
    {{ 20,  15,  22 }},
    {{ 20,  22,   3 }},
    {{ 4,  23,  22 }},
    {{ 23,   9,  24 }},
    {{ 22,  23,  24 }},
    {{ 22,  24,   3 }},
    {{ 4,  25,  23 }},
    {{ 25,  10,  26 }},
    {{ 23,  25,  26 }},
    {{ 23,  26,   9 }},
    {{ 4,  14,  25 }},
    {{ 14,   5,  27 }},
    {{ 25,  14,  27 }},
    {{ 25,  27,  10 }},
    {{ 5,  28,  27 }},
    {{ 28,  11,  29 }},
    {{ 27,  28,  29 }},
    {{ 27,  29,  10 }},
    {{ 5,  17,  28 }},
    {{ 17,   6,  30 }},
    {{ 28,  17,  30 }},
    {{ 28,  30,  11 }},
    {{ 6,  31,  30 }},
    {{ 31,   7,  32 }},
    {{ 30,  31,  32 }},
    {{ 30,  32,  11 }},
    {{ 6,  19,  31 }},
    {{ 19,   2,  33 }},
    {{ 31,  19,  33 }},
    {{ 31,  33,   7 }},
    {{ 2,  34,  33 }},
    {{ 34,   8,  35 }},
    {{ 33,  34,  35 }},
    {{ 33,  35,   7 }},
    {{ 2,  21,  34 }},
    {{ 21,   3,  36 }},
    {{ 34,  21,  36 }},
    {{ 34,  36,   8 }},
    {{ 3,  24,  36 }},
    {{ 24,   9,  37 }},
    {{ 36,  24,  37 }},
    {{ 36,  37,   8 }},
    {{ 9,  26,  39 }},
    {{ 26,  10,  38 }},
    {{ 39,  26,  38 }},
    {{ 39,  38,  12 }},
    {{ 10,  29,  38 }},
    {{ 29,  11,  40 }},
    {{ 38,  29,  40 }},
    {{ 38,  40,  12 }},
    {{ 11,  32,  40 }},
    {{ 32,   7,  41 }},
    {{ 40,  32,  41 }},
    {{ 40,  41,  12 }},
    {{ 7,  35,  41 }},
    {{ 35,   8,  42 }},
    {{ 41,  35,  42 }},
    {{ 41,  42,  12 }},
    {{ 8,  37,  42 }},
    {{ 37,   9,  39 }},
    {{ 42,  37,  39 }},
    {{ 42,  39,  12 }}
  } ;
MRI_SURFACE *
ic42_make_surface(int max_vertices, int max_faces)
{
  MRI_SURFACE *mris ;
  int         vno, fno, n, vn, n1, n2 ;
  VERTEX      *v ;
  FACE        *f ;
  static int first_time = 1 ;

  if (first_time)
  {
    first_time = 0 ;
    for (fno = 0 ; fno < ICO_NFACES ; fno++)
    {
      vno = ic42_faces[fno].vno[1] ;
      ic42_faces[fno].vno[1] = ic42_faces[fno].vno[2] ;
      ic42_faces[fno].vno[2] = vno ;
    }
  }

  mris = MRISoverAlloc(max_vertices, max_faces, ICO_NVERTICES, ICO_NFACES) ;

  /* position vertices */
  for (vno = 0 ; vno < ICO_NVERTICES ; vno++)
  {
    v = &mris->vertices[vno] ;

    v->x = 100.0*ic42_vertices[vno].x ;
    v->y = 100.0*ic42_vertices[vno].y ;
    v->z = 100.0*ic42_vertices[vno].z ;
  }

  /* fill in faces, and count # of faces each vertex is part of */
  for (fno = 0 ; fno < ICO_NFACES ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == 15)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      f->v[n] = ic42_faces[fno].vno[n]-1 ;  /* make it zero-based */
      v = &mris->vertices[f->v[n]] ;
      v->num++ ;
      v->vnum += 2 ;   /* will remove duplicates later */
    }
  }

  for (vno = 0 ; vno < ICO_NVERTICES ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->v = (int *)calloc(v->vnum/2, sizeof(int)) ;
    if (!v->v)
      ErrorExit(ERROR_NOMEMORY, "ic42: could not allocate %dth vertex list.",
                vno) ;
    v->vnum = 0 ;
  }

  /* now build list of neighbors */
  for (fno = 0 ; fno < ICO_NFACES ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == 3)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;

      /* now add an edge to other 2 vertices if not already in list */
      for (n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
      {
        if (n1 == n)   /* don't connect vertex to itself */
          continue ;
        vn = ic42_faces[fno].vno[n1]-1 ;  /* make it zero-based */

        /* now check to make sure it's not a duplicate */
        for (n2 = 0 ; n2 < v->vnum ; n2++)
        {
          if (v->v[n2] == vn)
          {
            vn = -1 ; /* mark it as a duplicate */
            break ;
          }
        }
        if (vn >= 0)
          v->v[v->vnum++] = vn ;
      }
    }
  }

  /* now allocate face arrays in vertices */
  for (vno = 0 ; vno < ICO_NVERTICES ; vno++)
  {
    v->vtotal = v->vnum ;
    v = &mris->vertices[vno] ;
    v->f = (int *)calloc(v->num, sizeof(int)) ;
    if (!v->f)
      ErrorExit(ERROR_NO_MEMORY,"ic42: could not allocate %d faces",v->num);
    v->n = (uchar *)calloc(v->num,sizeof(uchar));
    if (!v->n)
      ErrorExit(ERROR_NO_MEMORY, "ic42: could not allocate %d nbrs", v->n);
    v->num = 0 ;   /* for use as counter in next section */
    v->dist = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist)
      ErrorExit(ERROR_NOMEMORY,
                "mrisFindNeighbors: could not allocate list of %d "
                "dists at v=%d", v->vnum, vno) ;
    v->dist_orig = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist_orig)
      ErrorExit(ERROR_NOMEMORY,
                "mrisFindNeighbors: could not allocate list of %d "
                "dists at v=%d", v->vnum, vno) ;
    v->nsize = 1 ;
    v->vtotal = v->vnum ;
  }

  /* fill in face indices in vertex structures */
  for (fno = 0 ; fno < ICO_NFACES ; fno++)
  {
    f = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;
      v->n[v->num] = n ;
      v->f[v->num++] = fno ;
    }
  }

  MRIScomputeMetricProperties(mris) ;
#if 0
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    float dot ;
    int   ano ;

    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;

    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
       */
    v = &mris->vertices[f->v[0]] ;
    dot = v->x * f->nx + v->y * f->ny + v->z * f->nz;
    if (dot < 0.0f)   /* not in same direction, area < 0 and reverse n */
    {
      f->area *= -1.0f ;
      f->nx *= -1.0f;
      f->ny *= -1.0f;
      f->nz *= -1.0f;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        f->angle[ano] *= -1.0f ;
    }
  }
#endif
  mris->type = MRIS_ICO_SURFACE ;
  MRISsetNeighborhoodSize(mris, 1) ;
  return(mris) ;
}
