#pragma once
/*
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mrisurf_aaa.h"
#include "mrisurf.h"


// MRIS code dealing with the existence and connectedness of the vertices, edges, and faces
//                   and with their partitioning into sets (ripped, marked, ...)
//                   but not with their placement in the xyz coordinate space



int edgeExists(MRIS *mris, int vno1, int vno2);
int mrisRemoveLink(MRIS *mris, int vno1, int vno2);


// Each VERTEX has a list of neighbours at varying distances
// Sadly there are several algorithms that maintain this list
// This supports comparing those algorithms
//
#define MRIS_VertexNeighbourInfo_MAX_HOPS 5

typedef struct MRIS_VertexNeighbourInfo {
  size_t hops;
  int    vnum[MRIS_VertexNeighbourInfo_MAX_HOPS + 1];
  int    v   [MAX_NEIGHBORS];
} MRIS_VertexNeighbourInfo;

void MRIS_VertexNeighbourInfo_check            (MRIS_VertexNeighbourInfo* lhs,  MRIS_VertexNeighbourInfo* rhs);
void MRIS_VertexNeighbourInfo_load_from_VERTEX (MRIS_VertexNeighbourInfo* info, MRIS* mris, int vno);

void MRIS_check_vertexNeighbours(MRIS* mris);

void MRIS_setNsizeCur(MRIS *mris, int vno, int nsize);


// Vertices and Faces interact via edges
//
int mrisCountValidLinks(MRIS *mris, int vno1, int vno2);

int findFace(MRIS *mris, int vno0, int vno1, int vno2);
    // Return the index of the triangular face vno0-->vno1-->vno2, or -1
bool isFace(MRIS *mris, int vno0, int vno1, int vno2);
    // Like findFace but returns findFace(...) >= 0
    
int mrisValidFaces(MRIS *mris);
int vertexInFace(MRIS *mris, int vno, int fno);
int findOtherEdgeFace(MRIS const *mris, int fno, int vno, int vn1);

int mrisCountTotalNeighbors(MRIS *mris);

bool triangleMarked(MRIS *mris, int fno);               // are any of the face's vertices marked?
int findNonMarkedFace(MRIS *mris, int vno, int vn1);

int computeOrientation(MRIS *mris, int f, int v0, int v1);

static int vertexNeighbor(MRIS *mris, int vno1, int vno2)
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

int mrisAddFace(MRIS *mris, int vno0, int vno1, int vno2);

int mrisDivideEdgeTopologically(MRIS *mris, int vno1, int vno2);
