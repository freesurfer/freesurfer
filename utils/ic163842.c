/**
 * @file  ic163842.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.5 $
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
#include "utils.h" // fgetl

#define ICO_NVERTICES 163842
#define ICO_NFACES  327680


static int static_read_icosahedron(void) ;
IC_VERTEX *ic163842_vertices = NULL ;
IC_FACE  *ic163842_faces = NULL ;

MRI_SURFACE *
ic163842_make_surface(int max_vertices, int max_faces)
{
  MRI_SURFACE *mris ;
  int         vno, fno, n, vn, n1, n2 ;
  VERTEX      *v ;
  FACE        *f ;
  static int first_time = 1 ;

  if (first_time)
  {
    static_read_icosahedron() ;
    first_time = 0 ;
    for (fno = 0 ; fno < ICO_NFACES ; fno++)
    {
      vno = ic163842_faces[fno].vno[1] ;
      ic163842_faces[fno].vno[1] = ic163842_faces[fno].vno[2] ;
      ic163842_faces[fno].vno[2] = vno ;
    }
  }

  mris = MRISoverAlloc(max_vertices, max_faces, ICO_NVERTICES, ICO_NFACES) ;

  /* position vertices */
  for (vno = 0 ; vno < ICO_NVERTICES ; vno++)
  {
    v = &mris->vertices[vno] ;

    v->x = 100.0*ic163842_vertices[vno].x ;
    v->y = 100.0*ic163842_vertices[vno].y ;
    v->z = 100.0*ic163842_vertices[vno].z ;
  }

  /* fill in faces, and count # of faces each vertex is part of */
  for (fno = 0 ; fno < ICO_NFACES ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == 15)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      f->v[n] = ic163842_faces[fno].vno[n]-1 ;  /* make it zero-based */
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
      ErrorExit(ERROR_NOMEMORY, "ic163842: could not allocate %dth vertex list.",
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
        vn = ic163842_faces[fno].vno[n1]-1 ;  /* make it zero-based */

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
    v = &mris->vertices[vno] ;
    v->vtotal = v->vnum ;
    v->f = (int *)calloc(v->num, sizeof(int)) ;
    if (!v->f)
      ErrorExit(ERROR_NO_MEMORY,
                "ic163842: could not allocate %d faces",v->num);
    v->n = (uchar *)calloc(v->num,sizeof(uchar));
    if (!v->n)
      ErrorExit(ERROR_NO_MEMORY, "ic163842: could not allocate %d nbrs", v->n);
    v->num = 0 ;   /* for use as counter in next section */
    v->dist = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist)
      ErrorExit(ERROR_NOMEMORY,
                "ic163842: could not allocate list of %d "
                "dists at v=%d", v->vnum, vno) ;
    v->dist_orig = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist_orig)
      ErrorExit(ERROR_NOMEMORY,
                "ic163842: could not allocate list of %d "
                "dists at v=%d", v->vnum, vno) ;
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
  mris->type = MRIS_ICO_SURFACE ;
  free(ic163842_vertices) ;
  free(ic163842_faces) ;
  first_time = 1 ;
  return(mris) ;
}
static int
static_read_icosahedron(void)
{
  FILE      *fp ;
  char      line[200], *cp ;
  int       vno, fno, vno1, vno2, vno3, n ;
  float     x, y, z ;
  IC_VERTEX *ic_vertices ;
  IC_FACE   *ic_faces ;

  fp = fopen("ic163842.tri", "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "static_read_icosahedron: could not open %s",
              "ic163842.tri") ;

  ic_vertices = ic163842_vertices =
                  (IC_VERTEX *)calloc(ICO_NVERTICES, sizeof(IC_VERTEX)) ;
  if (!ic163842_vertices)
    ErrorExit(ERROR_NOMEMORY, "static_read_ico: could not allocate vertex list") ;
  ic_faces = ic163842_faces =
               (IC_FACE *)calloc(ICO_NFACES, sizeof(IC_FACE)) ;
  if (!ic163842_faces)
    ErrorExit(ERROR_NOMEMORY, "static_read_ico: could not allocate vertex list") ;


  fgetl(line, 150, fp) ;   /* discard # of vertices */

  /* first read vertices */
  n = 0 ;
  while ((cp = fgetl(line, 150, fp)) != NULL)
  {
    if (sscanf(cp, "%d %f %f %f\n", &vno, &x, &y, &z) < 4)
      break ;
    ic_vertices[vno-1].x = x ;
    ic_vertices[vno-1].y = y ;
    ic_vertices[vno-1].z = z ;
    if (++n >= ICO_NVERTICES)
      break ;
  }
  n = 0 ;
  fgetl(line, 150, fp) ;   /* discard # of faces */
  while ((cp = fgetl(line, 150, fp)) != NULL)
  {
    if (sscanf(cp, "%d %d %d %d\n", &fno, &vno1, &vno2, &vno3) < 4)
      break ;
    ic_faces[fno-1].vno[0] = vno1 ;
    ic_faces[fno-1].vno[1] = vno2 ;
    ic_faces[fno-1].vno[2] = vno3 ;
    if (++n >= ICO_NFACES)
      break ;
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

