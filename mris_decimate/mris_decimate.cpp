/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/09 21:10:45 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2010,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

///
/// \file mris_decimate.cpp
///
/// \brief Brief description
/// Reduce the number of vertices and faces in a surface.
///
/// \b NAME
///
/// mris_decimate
///
/// \b SYNPOSIS
///
/// mris_decimate [options] <input surface> <output surface>
///
/// \b DESCRIPTION
///
///  This tool reduces the number of triangles in a surface using the
///  the GNU Triangulated Surface Library documented at:
///
///           http://gts.sourceforge.net/reference/book1.html
///
///  Please see the GTS documentation for details on the decimation algorithm.
///  mris_decimate will read in an existing surface and write out a new one
///  with less triangles.  The decimation level and other options can be provided
///  on the command-line.
///


// $Id: mris_decimate.cpp,v 1.5 2011/02/09 21:10:45 nicks Exp $

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sstream>
#include <gts.h>
#include <iostream>
#include "mris_decimate.h"

extern "C"
{
#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
}


///////////////////////////////////////////////////////////////////////////
//
//  Private Functions
//
//


///
//  Callback function for loading vertex data from GtsSurface in
//  decimateSurface()
//
static void vertexLoad(GtsPoint * p, gpointer * data)
{
  MRI_SURFACE *mris = (MRI_SURFACE *)data;
  mris->vertices[mris->nvertices].x = p->x;
  mris->vertices[mris->nvertices].y = p->y;
  mris->vertices[mris->nvertices].z = p->z;
  unsigned int vertexID = mris->nvertices;
  GTS_OBJECT(p)->reserved = GUINT_TO_POINTER(vertexID);
  mris->nvertices++;
}

///
//  Callbackfunction for loading traingle data from GtsSurface in
//  decimateSurface()
//
static void faceLoad(GtsTriangle * t, gpointer * data)
{
  MRI_SURFACE *mris = (MRI_SURFACE *)data;
  GtsVertex *v1, *v2, *v3;
  gts_triangle_vertices(t, &v1, &v2, &v3);

  mris->faces[mris->nfaces].v[0] = (guint) (glong) GTS_OBJECT (v1)->reserved;
  mris->faces[mris->nfaces].v[1] = (guint) (glong) GTS_OBJECT (v2)->reserved;
  mris->faces[mris->nfaces].v[2] = (guint) (glong) GTS_OBJECT (v3)->reserved;
  mris->nfaces++;
}

///
//  From cleanup.c example in gts, use this callback function in
//  edgeCleanup to remove duplicate edges
//
static void buildList(gpointer data, GSList ** list)
{
  *list = g_slist_prepend (*list, data);
}

///
//  Cleanup any duplicate edges from GtsSurface.  This is needed because
//  the GTS library uses a data structure where each face edge needs
//  to be stored only once (vs. the way Freesurfer stores the edge in
//  each face).  This function goes through and cleans out any duplicate
//  edges from the mesh.  This is applied before the surface is
//  decimated.
//
//  This code was adapted from gts-0.7.6 examples/cleanup.c
//
static void edgeCleanup(GtsSurface * surface)
{
  GSList * edges = NULL;
  GSList * i;

  g_return_if_fail (surface != NULL);

  /* build list of edges */
  gts_surface_foreach_edge (surface, (GtsFunc) buildList, &edges);

  /* remove degenerate and duplicate edges.
   Note: we could use gts_edges_merge() to remove the duplicates and then
   remove the degenerate edges but it is more efficient to do everything
   at once (and it's more pedagogical too ...) */

  /* We want to control manually the destruction of edges */
  gts_allow_floating_edges = TRUE;

  i = edges;
  while (i)
  {
    GtsEdge * e = (GtsEdge*) i->data;
    GtsEdge * duplicate;
    if (GTS_SEGMENT (e)->v1 == GTS_SEGMENT (e)->v2) /* edge is degenerate */
      /* destroy e */
    {
      gts_object_destroy (GTS_OBJECT (e));
    }
    else if ((duplicate = gts_edge_is_duplicate (e)))
    {
      /* replace e with its duplicate */
      gts_edge_replace (e, duplicate);
      /* destroy e */
      gts_object_destroy (GTS_OBJECT (e));
    }
    i = i->next;
  }

  /* don't forget to reset to default */
  gts_allow_floating_edges = FALSE;

  /* free list of edges */
  g_slist_free (edges);
}

///////////////////////////////////////////////////////////////////////////////////
//
//  Public Functions
//
//


///
/// \fn int decimateSurface(MRI_SURFACE *mris, const DECIMATION_OPTIONS &decimationOptions)
/// \brief This function performs decimation on the input surface and outputs the new surface to a file
////     using GTS (GNUT Triangulated Surface Library)
/// \param mris Input loaded MRI_SURFACE to decimate
/// \param decimationOptions Options controlling the decimation arguments (see DECIMATION_OPTIONS)
/// \param decimateProgressFn If non-NULL, provides updates on decimation percentage complete and
///               a status message that can, for example, be displayed in a GUI.
/// \param userData If decimateProgressFn is non-NULL, argument passed into decimateProgressFn
/// \return 0 on success, 1 on failure
///
int decimateSurface(MRI_SURFACE **pmris,
                    const DECIMATION_OPTIONS &decimationOptions,
                    DECIMATE_PROGRESS_FUNC decimateProgressFn,
                    void *userData)
{
  MRI_SURFACE *mris = (*pmris);

  // Special case: if decimation level is 1.0, just make a copy of the
  // surface and return
  if (decimationOptions.decimationLevel == 1.0)
  {
    if (decimateProgressFn != NULL)
    {
      decimateProgressFn(1.0, "No decimation requested, finished.", userData);
    }

    return 0;
  }
  GtsSurface *gtsSurface = NULL;
  GtsVertex **gtsVertices = NULL;
  GtsEdge **gtsEdges = NULL;

  gtsSurface = gts_surface_new( gts_surface_class(),
                                gts_face_class(),
                                gts_edge_class(),
                                gts_vertex_class() );
  gtsVertices = (GtsVertex**) g_malloc( sizeof(GtsVertex*) * 
                                        (mris->nvertices) );

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.10, "Creating vertices...", userData);
  }

  for (int vno = 0; vno < mris->nvertices; vno++)
  {
    gtsVertices[vno] = gts_vertex_new ( gtsSurface->vertex_class,
                                        mris->vertices[vno].x,
                                        mris->vertices[vno].y,
                                        mris->vertices[vno].z );
  }

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.20, "Creating edges...", userData);
  }

  gtsEdges = (GtsEdge**) g_malloc(sizeof(GtsEdge*) * mris->nfaces * 3);
  int edge = 0;
  for (int fno = 0; fno < mris->nfaces; fno++)
  {
    GtsVertex *v1 = gtsVertices[mris->faces[fno].v[0]];
    GtsVertex *v2 = gtsVertices[mris->faces[fno].v[1]];
    GtsVertex *v3 = gtsVertices[mris->faces[fno].v[2]];

    gtsEdges[edge++] = gts_edge_new( gtsSurface->edge_class, v1, v2 );
    gtsEdges[edge++] = gts_edge_new( gtsSurface->edge_class, v2, v3 );
    gtsEdges[edge++] = gts_edge_new( gtsSurface->edge_class, v3, v1 );

  }

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.30, "Adding faces...", userData);
  }

  edge = 0;
  for (int fno = 0; fno < mris->nfaces; fno++)
  {
    gts_surface_add_face( gtsSurface,
                          gts_face_new(gtsSurface->face_class,
                                       gtsEdges[edge],
                                       gtsEdges[edge + 1],
                                       gtsEdges[edge + 2]) );
    edge += 3;
  }

  // The GTS data structure requires the edges to be merged.  Since
  // we did not cull duplicate edges from the Freesurfer data structure,
  // use this GTS function to merge the edges of the surface before
  // decimation.
  edgeCleanup(gtsSurface);
  int numDistinctEdges = gts_surface_edge_number(gtsSurface);

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.40, "Decimating surface...", userData);
  }

  float minimumAngle = 1.0;
  if (decimationOptions.setMinimumAngle)
  {
    minimumAngle = decimationOptions.minimumAngle;
  }

  GtsVolumeOptimizedParams params = { 0.5, 0.5, 0. };
  gdouble fold = minimumAngle * PI/180.;
  int stop = (int)((float)numDistinctEdges * 
                   decimationOptions.decimationLevel);
  gts_surface_coarsen (gtsSurface,
                       (GtsKeyFunc) gts_volume_optimized_cost,
                       &params,
                       (GtsCoarsenFunc) gts_volume_optimized_vertex,
                       &params,
                       (GtsStopFunc) gts_coarsen_stop_number,
                       &stop,
                       fold);

  // Free the vertex and face buffers
  for (int vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (mris->vertices[vno].f)
    {
      free(mris->vertices[vno].f) ;
      mris->vertices[vno].f = NULL ;
    }
    if (mris->vertices[vno].n)
    {
      free(mris->vertices[vno].n) ;
      mris->vertices[vno].n = NULL ;
    }
    if (mris->vertices[vno].dist)
    {
      free(mris->vertices[vno].dist) ;
      mris->vertices[vno].dist = NULL ;
    }
    if (mris->vertices[vno].dist_orig)
    {
      free(mris->vertices[vno].dist_orig) ;
      mris->vertices[vno].dist_orig = NULL ;
    }
    if (mris->vertices[vno].v)
    {
      free(mris->vertices[vno].v) ;
      mris->vertices[vno].v = NULL ;
    }
  }
  free(mris->vertices);
  free(mris->faces);

  GtsSurfaceStats stats;
  gts_surface_stats( gtsSurface, &stats);

  printf("Decimated Surface Number of vertices: %d\n", 
         stats.edges_per_vertex.n);
  printf("Decimated Surface Number of faces: %d\n",
         stats.n_faces);

  mris->nvertices = stats.edges_per_vertex.n;
  mris->nfaces = stats.n_faces;

  // Allocate at the new size
  mris->vertices = (VERTEX *)calloc(stats.edges_per_vertex.n, sizeof(VERTEX)) ;
  if (!mris->vertices)
  {
    ErrorExit(ERROR_NO_MEMORY,
              "decimateSurface(%d, %d): could not allocate vertices",
              mris->nvertices, sizeof(VERTEX));
  }
  mris->faces = (FACE *)calloc(stats.n_faces, sizeof(FACE)) ;
  if (!mris->faces)
  {
    ErrorExit(ERROR_NO_MEMORY,
              "decimateSurface(%d, %d): could not allocate faces",
              mris->nfaces, sizeof(FACE));
  }

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.70, "Copying vertices to mris...", userData);
  }

  gpointer data;
  mris->nvertices = 0;
  data = (gpointer)mris;
  gts_surface_foreach_vertex( gtsSurface, (GtsFunc) vertexLoad, data);


  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.80, "Copying faces to mris...", userData);
  }

  mris->nfaces = 0;
  data = (gpointer)mris;
  gts_surface_foreach_face( gtsSurface, (GtsFunc) faceLoad, data);


  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.90, "Creating new mris...", userData);
  }

  // Create a temporary file to write the output and then read it back in
  // to get the decimated surface. The reason I do this is because there
  // was no clear API call in mrisurf.c that would allow me to muck with
  // the vertex/face array and properly calculate everything else in the
  // structure.  I felt this was the safest way to make sure everything
  // in the surface got recalculated properly.
  char *tmpName = strdup("/tmp/mris_decimateXXXXXX");
  int fd = mkstemp(tmpName);
  char tmp_fpath[STRLEN];

  if (fd == -1)
  {
    std::cerr << "Error creating temporary file: " 
              << std::string(tmpName) << std::endl;
    return -1;
  }

  FileNameAbsolute(tmpName, tmp_fpath);

  MRISwrite(mris, tmp_fpath);
  MRISfree(pmris);
  *pmris = MRISread(tmp_fpath);
  remove(tmp_fpath);

  g_free(gtsVertices);
  g_free(gtsEdges);
  g_free(gtsSurface);

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(1.0, "Done.", userData);
  }

  return 0;
}




