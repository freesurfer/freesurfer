#pragma once
/*
 * @file utilities dealing with the topology 
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2018 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mrisurf_base.h"

// MRIS code dealing with the existence and connectedness of the vertices, edges, and faces
//                   and with their partitioning into sets (ripped, marked, ...)
//                   but not with their placement in the xyz coordinate space
int mrisSetVertexFaceIndex(MRI_SURFACE *mris, int vno, int fno);


int edgeExists(MRI_SURFACE *mris, int vno1, int vno2);
int triangleMarked(MRI_SURFACE *mris, int fno);

// Vertices and Faces interact via edges
//
int mrisCountValidLinks(MRIS *mris, int vno1, int vno2);
int isFace(MRIS *mris, int vno0, int vno1, int vno2);
int findFace(MRIS *mris, int vno0, int vno1, int vno2);
int mrisValidFaces(MRIS *mris);
int vertexInFace(MRIS *mris, int vno, int fno);
int findOtherEdgeFace(MRIS const *mris, int fno, int vno, int vn1);
int mrisRemoveLink(MRIS *mris, int vno1, int vno2);
int findNonMarkedFace(MRIS *mris, int vno, int vn1);

int computeOrientation(MRIS *mris, int f, int v0, int v1);

static int vertexNeighbor(MRI_SURFACE *mris, int vno1, int vno2)
{
  int n;

  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno1];
  for (n = 0; n < v->vnum; n++)
    if (v->v[n] == vno2) {
      return (1);
    }
  return (0);
}

int mrisMarkBadEdgeVertices(MRIS *mris, int mark);

#define MAX_VLIST 255
int mrisAddEdge(MRIS *mris, int vno1, int vno2);
int mrisAddFace(MRIS *mris, int vno0, int vno1, int vno2);

