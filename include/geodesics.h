//
// for calculating geodesics on a polyhedral surface
//

#ifndef geodesics
#define geodesics

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif
