//
// for calculating geodesics on a polyhedral surface
//

#ifndef geodesics
#define geodesics

#include "mrisurf.h"

#define MAX_GEODESICS 8000

typedef struct {
  int vnum;     // number of surrounding vertices within limit
  int v[MAX_GEODESICS];       // surrounding vertices
  float dist[MAX_GEODESICS];  // distances to vertices
} Geodesics;

// computes and returns the nearest geodesics for every vertex in the surface:
Geodesics* computeGeodesics(MRIS* surf, float maxdist);

// save/load geodesics:
void geodesicsWrite(Geodesics* geo, int nvertices, char* fname);
Geodesics* geodesicsRead(char* fname, int *nvertices);
int geodesicsUniquify(Geodesics *geod);

typedef struct {
  int vtxno; // surface vertex no
  int volindex; // voxel index from volume
  float dist;
}  VTXVOLINDEX;

MRI *GeoSmooth(MRI *src, double fwhm, MRIS *surf, Geodesics *geod, MRI *volindex, MRI *out);
int GeoCount(Geodesics *geod, int nvertices);
int GeoDumpVertex(char *fname, Geodesics *geod, int vtxno);
int geodesicsWriteV2(Geodesics* geo, int nvertices, char* fname) ;
Geodesics* geodesicsReadV2(char* fname, int *pnvertices) ;
double geodesicsCheckSphereDist(MRIS *sphere, Geodesics *geod);
double MRISsphereDist(MRIS *sphere, VERTEX *vtx1, VERTEX *vtx2);

int VtxVolIndexCompare(const void *a, const void *b);
int VtxVolIndexSort(VTXVOLINDEX *vvi, int nlist);
int VtxVolIndexPrint(VTXVOLINDEX *vvi, int nlist);
int VtxVolIndexSortTest(int nlist);
VTXVOLINDEX *VtxVolIndexUnique(VTXVOLINDEX *vvi, int nlist, int *nunique);
VTXVOLINDEX *VtxVolIndexPack(Geodesics *geod, int vtxno, MRI *volindex);
Geodesics *VtxVolPruneGeod(Geodesics *geod, int vtxno, MRI *volindex);

#endif
