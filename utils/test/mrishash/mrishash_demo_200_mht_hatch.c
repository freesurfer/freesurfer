/*--------------------------------------------
  mrishash_demo_200_mht_hatch.c

  Version History
  ----------------
  2007-04-05  GW Tidies
  2007-03-28  GW original

  Notes:
  ------
  This demo creates a surface, creates a MRIS_HATCH_TABLE from that and
  then copies the MHT data to an MRI volume for visual inspection.

  This shows the completeness and/or waste for a hatch algorithm.

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

char * Progname;

#define VERTEXCOUNT 6

GWUTILS_VERTEX vertices[VERTEXCOUNT] = {
  { 15,   0,   0 },
  {-20,   0,   0 },
  {  0,  25,   0 },
  {  0, -30,   0 },
  {  0,   0,  35 },
  {  0,   0, -40 }
};

#define FACECOUNT 8

GWUTILS_FACE faces[FACECOUNT] = {
  {{ 0, 5, 2 }},
  {{ 0, 2, 4 }},
  {{ 4, 2, 1 }},
  {{ 1, 2, 5 }},
  {{ 0, 3, 5 }},
  {{ 5, 3, 1 }},
  {{ 1, 3, 4 }},
  {{ 4, 0, 3 }}
};

//-----------------------------------
int main(int argc, char *argv[]) {
//-----------------------------------
  char msg[1000];
  char * cp;
  MRI_SURFACE * mris;
  MHT * mht;
  MRI * mri;
  int rslt;

  Progname = argv[0];
  printf("------------------------------\n");
  sprintf(msg, "Program: %s  V1.08", Progname);
  printf("%s\n", msg);

  cp = getenv("FREESURFER_HOME");
  printf("FREESURFER_HOME: %s\n", cp);

  printf("Make surface\n");

  //------------------------------------------
  // GW manual surface. Possibly change shape by changing
  // tables above.
  //------------------------------------------
  mris = GWU_make_surface_from_lists(vertices, VERTEXCOUNT, faces, FACECOUNT);

  printf("MRISwrite\n");
  MRISwrite(mris, "lh.mrishash_demo_200_mht_hatch.tri");

  //------------------------------------------
  // Hash the surface. Possibly try different
  // resolution arguments here
  //------------------------------------------
  printf("MHTfillTableAtResolution\n");
  mht = NULL;
  mht = MHTfillTableAtResolution(mris, mht, CURRENT_VERTICES, 1.0);

  printf("MRIFromMHTandMRIS\n");
  mri = MRIFromMHTandMRIS(mht, mris, MFMM_Count);

  printf("MRIwrite\n");
  rslt = MRIwrite(mri, "./mrishash_demo_200_mht_hatch_mri.mgz");
  printf("MRIWrite done: %d\n", rslt);

  printf("Done.\n");

  return 0;
}

