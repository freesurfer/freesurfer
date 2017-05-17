//
// for calculating geodesics on a polyhedral surface
//

#ifndef geodesics
#define geodesics

extern "C"
{

#include "mrisurf.h"

#define MAX_GEODESICS 4000


struct Geodesics {
  int vnum;     // number of surrounding vertices within limit
  int v[MAX_GEODESICS];       // surrounding vertices
  float dist[MAX_GEODESICS];  // distances to vertices
};

// computes and returns the nearest geodesics for every vertex in the surface:
Geodesics* computeGeodesics(MRIS* surf, int maxdist);

// save/load geodesics:
void GeodesicsWrite(Geodesics* geo, int nvertices, char* fname);
Geodesics* GeodesicsRead(char* fname, int nvertices);


}

#endif
