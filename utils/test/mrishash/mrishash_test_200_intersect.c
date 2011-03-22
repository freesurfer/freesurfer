/*--------------------------------------------
  Mrishash_test_200_intersect.c

  Version History
  ----------------
  2007-04-05  GW Revs to make compatible with MGH test procedure
  2007-03-30  GW original
  ----------------------------------------------*/

#define TestRepetitions  100

// Probably don't need all these includes...
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "version.h"
#include "icosahedron.h"
#include "gw_utils.h"

char * Progname;
char progver[] = "V1.11";
char logfilepath[1000];
char surffilepath[1000];
char mrifilepath[1000];

//------------------------------
void init_various(char * AProgname) {
//------------------------------
  int rslt;
  sprintf( logfilepath, "%s_log.txt", Progname);
  sprintf(surffilepath, "lh.%s_.tri", Progname);
  sprintf( mrifilepath, "%s_mri.mgz", Progname);

  rslt = gw_log_init(Progname, progver, logfilepath, 1); // empty file
  if (rslt) {
    printf("Couldn't open log file %s", logfilepath);
    exit(-1);
  }

}

//---------------------------------------------
int TestIntersectionDistance(int repnum) {
//---------------------------------------------
  int rslt = 0; // default OK
  char msg[1000];
  double radius1, radius2, rtot=0.0, maxradius, mhtres, pctspheredev;
  double offset;
  double vecx, vecy, vecz, veclen, osx=0.0, osy=0.0, osz=0.0;
  double jump = 70;
  int intersect;
  MRI_SURFACE * mris = NULL;
  msg[0]=0;

#define pctspheredevMAX 0.125

  //---------------------------------------------
  // Random radii
  //---------------------------------------------
  radius1 = (0.95 * ((double) rand() / RAND_MAX) + 0.5) * 40;//3*radius < 127
  radius2 = (0.95 * ((double) rand() / RAND_MAX) + 0.5) * 40;//3*radius < 127
  maxradius = radius1;
  if (maxradius < radius2) maxradius = radius2;

  offset = radius1 + radius2 + 5;
  jump   = offset - maxradius; // don't want one to end up inside the other

  //---------------------------------------------
  // Random mht resolution
  //---------------------------------------------
  mhtres = 1.0 + 3.9 * ((double) rand() / RAND_MAX);
  mhtres = floor(mhtres);// only really interested in integer values from 1..4

  //----------------------------------------
  // Random direction of offset
  //----------------------------------------
  // unit random vector
  vecx = ((double) rand() / RAND_MAX);
  vecy = ((double) rand() / RAND_MAX);
  vecz = ((double) rand() / RAND_MAX);
  veclen = sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
  vecx /= veclen;
  vecy /= veclen;
  vecz /= veclen;

  sprintf(msg, "[%d] Intersect test... Rs: %8.4f %8.4f = %8.4f offset: "
          "%8.4f  (%8.4f  %8.4f  %8.4f) res: %8.4f",
          repnum, radius1, radius2, rtot, offset, vecx, vecy, vecz, mhtres);
  printf("%s\n", msg);

  while (jump > 0.01) {
    osx = vecx * offset;
    osy = vecy * offset;
    osz = vecz * offset;

    if (mris) MRISfree(&mris);

    mris = ic2562_make_two_icos(0,0,0,radius1, osx,osy,osz,radius2 );

    intersect = MHTtestIsMRISselfIntersecting(mris, mhtres);

    jump *= 0.5;

    if (intersect) {
      offset += jump;
    } else {
      offset -= jump;
    }
  }
  rtot = radius1 + radius2;
  //------------------------------------------------------------
  // Calculate percent deviation from touching spheres.
  // Note that collision will happen at the furthest at R1 + R2. But
  // collision may be closer than that because of the flat faces
  // that are inside boundary of sphere.
  //------------------------------------------------------------
  pctspheredev = 100 * (rtot - offset)/rtot;
  if ( (-0.02 > pctspheredev) || (pctspheredev > pctspheredevMAX) ) {
    rslt = 1;
  }

  sprintf(msg, "%d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %d",
          repnum, radius1, radius2, rtot, offset,
          osx, osy, osz, mhtres, pctspheredev, rslt);
  printf("%s\n", msg);
  gw_log_message(msg);

  return rslt;
}

//-----------------------------------
int main(int argc, char *argv[]) {
//-----------------------------------
  char msg[1000];
  int n;
  char * cp;
  int rslt = 0; // default to OK

  if (getenv("SKIP_MRISHASH_TEST")) exit(77); // bypass

  Progname = argv[0];
  init_various(Progname);  // and gw_log_init

  gw_log_begin();

  printf("------------------------------\n");
  printf("Program: %s\n", Progname);

  gw_log_message(msg);

  cp = getenv("FREESURFER_HOME");
  printf("FREESURFER_HOME: %s\n", cp);

  //------------------------------------------
  // Log column heads
  //------------------------------------------
  gw_log_message("repnum radius1 radius2 rtot offset "
                 "osx osy osz mhtres pctspheredev rslt");

  //------------------------------------------
  // Randomize -- may want this or not
  //------------------------------------------
  srand((unsigned int) time((time_t *) NULL) );

  //------------------------------------------
  // Call collision exercise routine.
  //------------------------------------------
  for (n = 0; n < TestRepetitions; n++) {
    rslt = TestIntersectionDistance(n);
    if (rslt) goto done;
  }

 done:

  gw_log_end();

  return rslt;
}

