#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
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



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sstream>
#include <gts.h>
#include <iostream>
#include "mris_decimate.h"



#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"



///////////////////////////////////////////////////////////////////////////
//
//  Private Functions
//
//


///
//  Callback function for loading vertex data from GtsSurface in
//  decimateSurface()
//
struct VertexLoadContext {
    MRI_SURFACE *mris;
    int          nextVertex;
};

static void vertexLoad(GtsPoint * p, gpointer * data)
{
  VertexLoadContext* vertexLoadContext = (VertexLoadContext*)data;

  MRIS* mris = vertexLoadContext->mris;

  MRISsetXYZ(mris, vertexLoadContext->nextVertex,
    p->x,
    p->y,
    p->z);
  
  unsigned int vertexID = vertexLoadContext->nextVertex;
  GTS_OBJECT(p)->reserved = GUINT_TO_POINTER(vertexID);
  vertexLoadContext->nextVertex++;
}

///
//  Callbackfunction for loading traingle data from GtsSurface in
//  decimateSurface()
//
struct FaceLoadContext {
    MRI_SURFACE *mris;
    int          nextFace;
};

static void faceLoad(GtsTriangle * t, gpointer * data)
{
  FaceLoadContext* faceLoadContext = (FaceLoadContext*)data;
  MRI_SURFACE *mris = faceLoadContext->mris;
  
  GtsVertex *v1, *v2, *v3;
  gts_triangle_vertices(t, &v1, &v2, &v3);

  mris->faces[faceLoadContext->nextFace].v[0] = (guint) (glong) GTS_OBJECT (v1)->reserved;
  mris->faces[faceLoadContext->nextFace].v[1] = (guint) (glong) GTS_OBJECT (v2)->reserved;
  mris->faces[faceLoadContext->nextFace].v[2] = (guint) (glong) GTS_OBJECT (v3)->reserved;
  faceLoadContext->nextFace++;
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
  MRI_SURFACE *mris0 = MRISclone(mris); // make a copy

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
  // BRB 08/23/18 this seems wrong to me because it doesn't change the counts
  for (int vno = 0 ; vno < mris->nvertices ; vno++)
  {
    freeAndNULL(mris->vertices_topology[vno].f);
    freeAndNULL(mris->vertices_topology[vno].n);
    freeAndNULL(mris->vertices_topology[vno].v);
  }
  MRISfreeDistsButNotOrig(mris);
  MRISfreeDistOrigs      (mris);

  // DNG 7/16/18: these two frees were in the original. I don't know how it ever worked
  // because the realloc below needs these to be valid pointers
  //free(mris->vertices);
  //free(mris->faces);

  GtsSurfaceStats stats;
  gts_surface_stats( gtsSurface, &stats);

  printf("Decimated Surface Number of vertices: %d\n", 
         stats.edges_per_vertex.n);
  printf("Decimated Surface Number of faces: %d\n",
         stats.n_faces);

  
  MRISreallocVerticesAndFaces(mris, stats.edges_per_vertex.n, stats.n_faces);

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.70, "Copying vertices to mris...", userData);
  }

  VertexLoadContext vertexLoadContext;
  	vertexLoadContext.mris     = mris;
	vertexLoadContext.nextVertex = 0;
  gts_surface_foreach_vertex( gtsSurface, (GtsFunc) vertexLoad, (gpointer)&vertexLoadContext);

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.80, "Copying faces to mris...", userData);
  }

  FaceLoadContext faceLoadContext;
  	faceLoadContext.mris     = mris;
	faceLoadContext.nextFace = 0;
  gts_surface_foreach_face( gtsSurface, (GtsFunc) faceLoad, (gpointer)&faceLoadContext);

  if (!mrisCheckVertexFaceTopology(mris)) {
    std::cerr << "Error: surface has invalid topology" << std::endl;
    return 1;
  }

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(0.90, "Creating new mris...", userData);
  }

  MRISreallocVerticesAndFaces(mris, vertexLoadContext.nextVertex, faceLoadContext.nextFace);
  
  // Create a temporary file to write the output and then read it back in
  // to get the decimated surface. The reason I do this is because there
  // was no clear API call in mrisurf.c that would allow me to muck with
  // the vertex/face array and properly calculate everything else in the
  // structure.  I felt this was the safest way to make sure everything
  // in the surface got recalculated properly.
  char *tmpName = strdup("/tmp/mris_decimateXXXXXX");
  int fd = mkstemp(tmpName);
  char tmp_fpath[STRLEN];
  if (fd == -1)  {
    std::cerr << "Error creating temporary file: "<< std::string(tmpName) << std::endl;
    return -1;
  }
  FileNameAbsolute(tmpName, tmp_fpath);
  MRISwrite(mris, tmp_fpath);
  MRISfree(pmris);
  *pmris = MRISread(tmp_fpath);
  remove(tmp_fpath);
  mris = *pmris;

  if(decimationOptions.SortVertices){
    // The gts_surface_coarsen() function will always produce the same
    // surface for the same input in that the xyz of the vertices will
    // be the same. However, when GTS is compiled with hashes, the
    // order of the vertices and faces will be different from run to
    // run. The code below makes it so the output will always be the
    // same. When GTS is compiled with BTREES (which we just changed
    // to in July 2018), it will produce the same output. So, by
    // default, the sorting code is not run. It can be turned on with
    // mris_decimate -q. There may still be some non-deterministic
    // behavior. For example, when the orig.nofix is input. Not sure
    // why, but probably because the lengths of the edges are all
    // either 1 or sqrt(2) thus creating some abiguity which is
    // handled differently on different runs.
    printf("Sorting surface vertices\n");
    MRI_SURFACE *mris2;
    mris2 = MRISsortVertices(mris);
    MRISfree(&mris);
    mris = mris2;
    *pmris = mris;
  }

  g_free(gtsVertices);
  g_free(gtsEdges);
  g_free(gtsSurface);

  if (decimateProgressFn != NULL)
  {
    decimateProgressFn(1.0, "Done.", userData);
  }

  MRISfree(&mris0);
  return 0;
}




