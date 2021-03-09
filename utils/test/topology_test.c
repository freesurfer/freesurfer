/**
 * @brief mrisurf_topology tests
 *
 */
/*
 * Original Author: B. R. Brett
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

#include "mrisurf.h"
#include "utils.h"
#include <stdio.h>

#define CHECK(COND) if (!(COND)) { printf("%s failed at line %d\n", #COND, __LINE__); fails++; }

const char *Progname = "topology_test";

int test1()
{
  int fails = 0;

  // Make the src
  //
  int const max_vertices = 4, max_faces = 2, 
               nvertices = 4,    nfaces = 2;

  MRIS mrisObject;
  MRIS* src = &mrisObject;
  
  MRISctr(src, max_vertices, max_faces, nvertices, nfaces);
  CHECK(src->nvertices == 4);
  CHECK(src->nfaces == 2);
   
  mrisAddEdge(src, 0, 1);                   // will get removed
  mrisAddEdge(src, 0, 2);
  mrisAddEdge(src, 0, 3);
  mrisAddEdge(src, 1, 2);
  mrisAddEdge(src, 2, 3);
  mrisAddEdge(src, 3, 1);

  mrisAttachFaceToVertices(src, 0, 0,1,2);  // will get removed
  mrisAttachFaceToVertices(src, 1, 1,2,3);  // will survive

  CHECK(src->nvertices == 4);
  CHECK(src->nfaces == 2);
   
  // Copy the src
  //
  size_t     dst_nvertices;
  int const* dst_mapToDstVno;
  size_t     dst_nfaces;
  int const* dst_mapToDstFno;

  MRIS* dst;

  MRIS_HASH srcHash, dstHash;
  
  MRIScreateSimilarTopologyMapsForNonripped(
    src,
    &dst_nvertices,
    &dst_mapToDstVno,
    &dst_nfaces,
    &dst_mapToDstFno);
  CHECK(dst_nvertices == 4);
  CHECK(dst_nfaces == 2);
  
  dst =
    MRIScreateWithSimilarTopologyAsSubset(
      src,
      dst_nvertices,
      dst_mapToDstVno,
      dst_nfaces,
      dst_mapToDstFno);
  freeAndNULL(dst_mapToDstVno);
  freeAndNULL(dst_mapToDstFno);
  
  CHECK(dst->nvertices == 4);
  CHECK(dst->nfaces == 2);
  
  mris_print_diff(stdout, src, dst);
  mris_hash_init(&srcHash, src);
  mris_hash_init(&dstHash, dst);
  CHECK(srcHash.hash == dstHash.hash);

  MRISfree(&dst);

  // Subset the src
  //
  src->vertices[0].ripflag = 1;
  src->faces   [0].ripflag = 1;
  
  MRIScreateSimilarTopologyMapsForNonripped(
    src,
    &dst_nvertices,
    &dst_mapToDstVno,
    &dst_nfaces,
    &dst_mapToDstFno);
  CHECK(dst_nvertices == 3);
  CHECK(dst_nfaces == 1);
      
  dst =
    MRIScreateWithSimilarTopologyAsSubset(
      src,
      dst_nvertices,
      dst_mapToDstVno,
      dst_nfaces,
      dst_mapToDstFno);
  freeAndNULL(dst_mapToDstVno);
  freeAndNULL(dst_mapToDstFno);

  CHECK(dst->nvertices == 3);
  CHECK(dst->nfaces == 1);

  // TODO verify the shape
  //

  // Done
  //  
  MRISfree(&dst);
  MRISdtr(src);

  return fails;
}

int bitCount(int x) {
  int count = 0;
  int i = 1;
  while (x) {
    if (x&i) { x ^= i; count++; }
    i *= 2;
  }
  return count;
}

int test2Wkr(int verticesLog2)
{
  if (bitCount(5) != 2) {
    printf("bitCount(5) != 2\n");
    return 1;
  }
  
  int fails = 0;

  // Make the src
  //
  int const max_vertices = 1<<verticesLog2, max_faces = 0, 
               nvertices = max_vertices,    nfaces    = 0;

  MRIS mrisObject;
  MRIS* src = &mrisObject;
  
  MRISctr(src, max_vertices, max_faces, nvertices, nfaces);

  { int i,j;
    for (i = 0; i < nvertices; i++) 
    for (j = 0; j < i; j++)
      if (bitCount(i^j) == 1) 
        mrisAddEdge(src, i, j);
  }
  
  size_t const neighborsCapacity = nvertices;
  int* neighbors = (int*)malloc(neighborsCapacity*sizeof(int));
  int* hops      = (int*)malloc(neighborsCapacity*sizeof(int));
  
  int const vno1 = 0; 
  int const nlinks = 4; 
  int const neighborsSize = MRISfindNeighborsAtVertex(src, vno1, nlinks, neighborsCapacity, neighbors, hops);
  
  { int i;
    for (i = 0; i < neighborsSize; i++) {
      src->vertices[neighbors[i]].marked = hops[i];
    }
  }

  int vno2;
  for (vno2 = 0; vno2 < nvertices; vno2++) {
    int expected = bitCount(vno1^vno2);
    if (expected > nlinks) expected = 0;
    if (src->vertices[vno2].marked != expected) {
      fails++;
      if (fails == 1) printf("FAIL when verticesLog2:%d\n", verticesLog2);
      if (fails > 10) goto Done;
      printf("src->vertices[vno2:%d]->marked:%d != expected:%d\n", vno2, src->vertices[vno2].marked, expected);
    }
  }
  
  // Randomly rip a few vertices
  //
  for (vno2 = 0; vno2 < nvertices; vno2++) {
    if (vno2%31 == 0) src->vertices[vno2].ripflag = 1;
  }
  MRISremoveRipped(src);
  
  mrisCheckVertexFaceTopologyWkr(__FILE__,__LINE__,src,true);
  
  // Done
  //  
Done:
  MRISdtr(src);

  return fails;
}


int test2() {
         test2Wkr(2);
         test2Wkr(3);
  return test2Wkr(8);
}


int main(int argc, char *argv[])
{
    return 
      //test1() ? 1 : 
        test2() ? 1 :
        0;
}
