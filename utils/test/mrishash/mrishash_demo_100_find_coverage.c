/*--------------------------------------------
  mrishash_demo_100_find_coverage.c

  Version History
  ----------------
  2007-04-05  GW Tidies
  2007-03-28  GW original

  Notes:
  ------
  Demos the spatial coverage of the find function.  Program creates a small
  surface (1x1x1 cube) then scans a 3-D range of points around this
  surface to probe for "nearest vertex". For each point, the result is
  painted to an MRI.

  If using V1.27 version of mrishash, the "result" is simply a value
  (50) showing that a vertex was found.

  If using GWs0070404 version of mrishash, the "result" is a number indicating
  how many MHT grid boxes were inspected to find the nearest vertex. Since that
  number could be upto 7 x 7 x 7, and this exceeds 255 (the MRI data type max),
  instead, this number is plotted: 50 + boxes/2 (see code).

  Also painted in the MRI is the small surface itself.

  View results in tkmedit, or other MRI viewer.

  Things to adjust:
  -- hashres: Grid resolution of the hash table
  -- position of surface
  -- different find functions

  ----------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "version.h"

#include "gw_utils.h"


#define VERTEXCOUNT 8

GWUTILS_VERTEX vertices[VERTEXCOUNT] = {
  {0, 0, 0},
  {0, 0, 1},
  {0, 1, 0},
  {0, 1, 1},
  {1, 0, 0},
  {1, 0, 1},
  {1, 1, 0},
  {1, 1, 1}
};

#define FACECOUNT 12

GWUTILS_FACE faces[FACECOUNT] = {
  {{3, 2, 6}},
  {{3, 6, 7}},
  {{1, 3, 7}},
  {{1, 7, 5}},
  {{0, 1, 5}},
  {{0, 5, 4}},
  {{2, 0, 4}},
  {{2, 4, 6}},
  {{6, 4, 5}},
  {{6, 5, 7}},
  {{3, 1, 0}},
  {{3, 0, 2}}
};

//---------------- main at bottom ---------------
char * Progname;

#define int_bkg       10
#define int_default   20
#define int_base      50
#define int_cube     255
//-----------------------------------
void CoverageTest() {
//-----------------------------------
  const int xyzmid = 127;
  const float hashres = 16.0;

  MRI * amri;
  MRIS * amris;
  MHT  * amht;
  int rslt;
  int  xsi, ysi, zsi, xvi, yvi, zvi;
  int intens, halfscansize, closest_vertnum, vno;
  Real xs,  ys,  zs,  xv,  yv,  zv;  // WTF is Real?
  VERTEX vtx, *v;
  double ddist;

  int BucketsChecked, BucketsPresent;

  amri = MRIalloc(256, 256, 256, MRI_UCHAR);

  //-------------------------------
  // Set some voxels as background test
  //-------------------------------
  for (zvi = xyzmid - 30; zvi <= xyzmid+30; zvi++) {
    for (yvi = xyzmid - 30; yvi <= xyzmid+30; yvi++) {
      for (xvi = xyzmid - 30; xvi <= xyzmid+30; xvi++) {

        MRIvox(amri, xvi, yvi, zvi) = int_bkg;
      }
    }
  }

  //-------------------------------
  // Make a small cube at 0,0,0
  //-------------------------------
  amris = GWU_make_surface_from_lists(vertices,
                                      VERTEXCOUNT,
                                      faces,
                                      FACECOUNT);

  //-------------------------------
  // Move surface to some test position
  //-------------------------------
  for (vno = 0; vno < amris->nvertices; vno++) {
    v = &amris->vertices[vno];
    v->x += 6;
    v->y += 6;
    v->z += 6;
  }

  //-------------------------------
  // Create MHT using surface
  //-------------------------------
  printf("Creating MHT\n");
  amht = MHTfillVertexTableRes(amris, NULL,CURRENT_VERTICES, hashres);
  if (!amht) {
    printf("MHTfillVertexTableRes failed\n");
    goto done;
  }

  //-------------------------------
  // Test a range of points against that surface, and mark them
  // in volume as to whether a
  // vertex is findable from there. If so give it a value
  // that's some metric, like how
  // many hash boxes had to be examined to find vertex from that point.
  //-------------------------------
  printf("Starting closest-vertex scan\n");
  halfscansize = 64;
  for (zsi =         -halfscansize; zsi <= halfscansize; zsi++) {
    printf("zsi: %d\n", zsi);
    for (ysi =     -halfscansize; ysi <= halfscansize; ysi++) {
      for (xsi = -halfscansize; xsi <= halfscansize; xsi++) {
        // Surface-space coords
        xs = xsi;
        ys = ysi;
        zs = zsi;

        //------- Calc the voxel coords --------
        MRIsurfaceRASToVoxel(amri, xs, ys, zs, &xv, &yv, &zv);
        xvi=xv;
        yvi=yv;
        zvi=zv;

        // Copy surface-space coords to vertex
        vtx.x = xs;
        vtx.y = ys;
        vtx.z = zs;

        intens = int_default;

        //----------------------------------------
        // Here you could exercise any of the find functions to discover
        // their behavior
        //----------------------------------------
        MHTfindClosestVertexGeneric(
          amht, amris,
          xs, ys, zs,
          1000, 3,
          NULL, &closest_vertnum, &ddist
          );
        if (closest_vertnum >= 0) {
          MHTfindReportCounts(&BucketsChecked, &BucketsPresent);
          intens = int_base + (BucketsChecked >> 1);
        }
        MRIvox(amri, xvi, yvi, zvi) = intens;
      }
    }
  }
  printf("Done closest-vertex scan\n");

  //--------------------------------
  // Mark center of cube
  //--------------------------------
  // Note, really should mark all voxels hit by
  // vertices in surface

  for (vno = 0; vno < amris->nvertices; vno++) {
    v = &amris->vertices[vno];

    MRIsurfaceRASToVoxel(amri, v->x, v->y, v->z, &xv, &yv, &zv);

    xvi = xv;
    yvi = yv;
    zvi = zv;

    MRIvox(amri, xvi, yvi, zvi) = int_cube;
  }

  rslt = MRIwrite(amri, "./mrishash_demo_100_find_coverage_mri.mgz");
  printf("MRIWrite: %d\n", rslt);
 done:
  return;
}

//-----------------------------------
int main(int argc, char *argv[]) {
//-----------------------------------
  //   char msg[500];
  char * cp;
  Progname = argv[0];

  printf("Program: %s  V1.07: \n", Progname);

  cp = getenv("FREESURFER_HOME");
  printf("FREESURFER_HOME: %s\n", cp);

  CoverageTest();

  printf("Done.\n");

  return 0;
}

