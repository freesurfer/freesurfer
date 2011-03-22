/*--------------------------------------------
  mrishash_test_100_find_tests.c

  The test code has to test:

  1. That a vertex is returned when it should be.
  2. That the vertex that is returned is the closest.
  3. That no vertex is returned if the MHT code can't tell if it's 
  the closests.

  Version History
  ----------------
  2007-04-05  GW Revs to make compatible with MGH test procedure
  2007-03-30  GW original

  ----------------------------------------------*/

#define TestRepetitions    10
#define ProbeRepetitions  100
#define PURPOSE purpose_Test

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "version.h"

#include "gw_utils.h"
#include "icosahedron.h"

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

typedef enum {
  purpose_Test = 0,
  purpose_ListResults
} purpose_t;

//---------------------------------------
int TestNearestVtxInConcentricIcos(int surfacenum,
                                   purpose_t purpose) {
//---------------------------------------
  int errnum, rslt = 0; // default OK
  int pick, probeix;
  char msg[1000];
  double radius1, radius2, radiusdif, mhtres=0.0, searchrange=0.0;
  double probedistance, probex, probey, probez;
  double vecx, vecy, vecz, veclen;
  MRI_SURFACE * mris = NULL;
  MRIS_HASH_TABLE * mht = NULL;
  VERTEX probe_vtx, *brute_vtx, *fcv_vtx, *fcvit_vtx, *gnrc_vtx;
  int brute_vno, fcvn_vno;
  int LegacyAnswerExpected, GenericAnswerExpected;
  int gnrc_vno;
  double gnrc_dist;
  float brute_dist, fcvn_dist, fcv_dist, fcvit_dist;

  //------------------------------------------
  // pick some initial random dimensions
  //------------------------------------------
  pick = floor(6.9 * ((double) rand() / RAND_MAX));
  switch(pick) {
  case 0:
  case 1:
  case 2: mhtres =  1.0;   break;
  case 3: mhtres =  2.0;   break;
  case 4: mhtres =  4.0;   break;
  case 5: mhtres =  8.0;   break;
  case 6: mhtres = 16.0;   break;
  }

  pick = floor(4.9 * ((double) rand() / RAND_MAX));
  switch(pick) {
  case 0:
  case 1:
  case 2: searchrange = 0.5 * mhtres;   break;
  case 3: searchrange = 1.0 * mhtres;   break;
  case 4: searchrange = 2.0 * mhtres;   break;
  }

  radius1   = 10 + 64 * ((double) rand() / RAND_MAX);
  radiusdif = 2.3 * searchrange;// so sometimes nearest vertex is out-of-range
  radius2   = radius1 + radiusdif;

  if (mris) MRISfree(&mris);
  mris = ic2562_make_two_icos(0,0,0,radius1, 0,0,0, radius2);
  mht  = MHTfillVertexTableRes(mris, mht, CURRENT_VERTICES, mhtres);

  for (probeix = 0; probeix < ProbeRepetitions; probeix++) {
    probedistance = radius1 + radiusdif * ((double) rand() / RAND_MAX);

    // unit random vector
    vecx = ((double) rand() / RAND_MAX);
    vecy = ((double) rand() / RAND_MAX);
    vecz = ((double) rand() / RAND_MAX);
    veclen = sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
    vecx /= veclen;
    vecy /= veclen;
    vecz /= veclen;
    probex = vecx * probedistance;
    probey = vecy * probedistance;
    probez = vecz * probedistance;

    probe_vtx.x = probex;
    probe_vtx.y = probey;
    probe_vtx.z = probez;

    //------------------------------------
    // Init some variables
    //------------------------------------
    brute_vtx  = NULL;
    fcv_vtx    = NULL;
    fcvit_vtx  = NULL;
    gnrc_vtx   = NULL;
    brute_vno  = -2;
    fcvn_vno   = -2;
    gnrc_vno   = -2;
    brute_dist = 999.9;
    fcvn_dist  = 999.9;
    fcv_dist   = 999.9;
    fcvit_dist = 999.9;
    gnrc_dist  = 999.9;

    //------------------------------------
    // brute force
    //------------------------------------
    brute_vno = MRISfindClosestVertex(mris,
                                      probex, probey, probez,
                                      &brute_dist);
    if (brute_vno >=0)
      brute_vtx = &(mris->vertices[brute_vno]);

    // Should distance be returned?
    LegacyAnswerExpected  = (brute_dist <= mhtres     ) ? 1 : 0;
    GenericAnswerExpected = (brute_dist <= searchrange) ? 1 : 0;

    //------------------------------------
    // Using mht
    // Test all functions
    //------------------------------------
    fcv_vtx    = MHTfindClosestVertex(mht, mris, &probe_vtx) ;

    fcvit_vtx  = MHTfindClosestVertexInTable(mht,
                                             mris,
                                             probex, probey, probez,0) ;
    fcvn_vno   = MHTfindClosestVertexNo(mht,
                                        mris,
                                        &probe_vtx,
                                        &fcvn_dist);

    MHTfindClosestVertexGeneric(mht,
                                mris,
                                probex, probey, probez,
                                searchrange,
                                -1,
                                &gnrc_vtx,
                                &gnrc_vno,
                                &gnrc_dist);

    //------------------------------------
    // Calc some distances if possible
    //------------------------------------
    if (fcv_vtx)
      fcv_dist     = sqrt(SQR(probex-fcv_vtx->x)
                          + SQR(probey-fcv_vtx->y)
                          + SQR(probez-fcv_vtx->z));
    if (fcvit_vtx)
      fcvit_dist   = sqrt(SQR(probex-fcvit_vtx->x)
                          + SQR(probey-fcvit_vtx->y)
                          + SQR(probez-fcvit_vtx->z)); ;
    //  fcvn_dist done
    errnum = 0;
    if (LegacyAnswerExpected) {
      if (fcv_vtx   != brute_vtx) { errnum = 1; goto error_found; }
      if (fcvit_vtx != brute_vtx) { errnum = 2; goto error_found; }
      if (fcvn_vno  != brute_vno) { errnum = 3; goto error_found; }
      if (abs(fcvn_dist - brute_dist) > 0.01) {
        errnum = 4; goto error_found; }
    } else {  /* ... did legacy functions return a vertex, because they 
                 should NOT, since code hasn't looked at all vertices 
                 within that distance. */
      if (fcv_vtx)          { errnum = 11; goto error_found; }
      if (fcvit_vtx)        { errnum = 12; goto error_found; }
      if (fcvn_vno > 0)     { errnum = 13; goto error_found; }
    }

    if (GenericAnswerExpected) {
      if (gnrc_vtx != brute_vtx) { errnum = 5; goto error_found; }
      if (gnrc_vno != brute_vno) { errnum = 6; goto error_found; }
      if (abs(gnrc_dist - brute_dist) > 0.01) {
        errnum = 7;  goto error_found; }
    } else {
      if (gnrc_vtx)       { errnum = 15; goto error_found; }
      if (gnrc_vno > -1)  { errnum = 16; goto error_found; }
      if (gnrc_dist < 100) { errnum = 17; goto error_found; }
    }

  error_found:
    if (errnum) 
      rslt = 1;

    switch (purpose) {
    case purpose_Test:

      sprintf(msg, "%5d  %8.4f %8.4f %8.4f %5d %8.4f %8.4f "
              "%8.4f %6d %8.4f %d %d %d %d",
              surfacenum, radius1, radius2, mhtres, probeix,
              probex, probey, probez,
              brute_vno, brute_dist, LegacyAnswerExpected,
              GenericAnswerExpected, errnum, rslt);
      gw_log_message(msg);
      if (errnum || ((probeix & 0x3F) == 0) ) {
        printf("%s\n", msg);
      }
      if (rslt) goto done; // quit on first error if testing
      break;

    case purpose_ListResults:
      strcpy(msg, "------------------------------");
      gw_log_message(msg); printf("%s\n", msg);
      //srf res   ix   prx   pry   prz      LA GA er rslt
      sprintf(msg, "%5d  %8.4f %6d ( %8.4f %8.4f %8.4f ) %5d %d %d",
              surfacenum, mhtres, probeix, probex, probey, probez,
              LegacyAnswerExpected, GenericAnswerExpected, errnum);
      gw_log_message(msg); printf("%s\n", msg);

      sprintf(msg, "brute   %s %6d  %8.4f", "Y", brute_vno, brute_dist);
      gw_log_message(msg); printf("%s\n", msg);
      sprintf(msg, "fcv     %s    ???  %8.4f",
              fcv_vtx  ==brute_vtx ? "Y" : "N", fcv_dist);
      gw_log_message(msg); printf("%s\n", msg);
      sprintf(msg, "fcvit   %s    ???  %8.4f",
              fcvit_vtx==brute_vtx ? "Y" : "N", fcvit_dist);
      gw_log_message(msg); printf("%s\n", msg);
      sprintf(msg, "fcvn    %s %6d  %8.4f"   ,
              fcvn_vno==brute_vno  ? "Y" : "N",  fcvn_vno, fcvn_dist);
      gw_log_message(msg); printf("%s\n", msg);
      sprintf(msg, "gnrc    %s %6d  %8.4f"   ,
              gnrc_vtx==brute_vtx ? "Y" : "N",  gnrc_vno, gnrc_dist);
      gw_log_message(msg); printf("%s\n", msg);

      // keep on going if listing
      break;
    } // switch

  }
  goto done;
 done:

  return rslt;
}
char testhead[] =
"surface rad1 rad2 mhtres rep probex probey probez "
"brute_vno brute_dist LegAns GnrcAns err rslt";
char listhead1[] =
"surface mhtres rep probex probey probez LegAns GnrcAns errnum";
char listhead2[] =
"func correct vno dist";

//-----------------------------------
int main(int argc, char *argv[]) {
//-----------------------------------
  char msg[1000];
  int n;
  char * cp;
  int rslt = 0; // default to OK
  purpose_t purpose = PURPOSE;

  msg[0]=0;

  if (getenv("SKIP_MRISHASH_TEST")) exit(77); // bypass

  Progname = argv[0];
  init_various(Progname);  // and gw_log_init

  gw_log_begin();

  printf("------------------------------\n");
  printf("Program: %s\n", Progname);

  printf("%s\n", msg);
  gw_log_message(msg);

  cp = getenv("FREESURFER_HOME");
  printf("FREESURFER_HOME: %s\n", cp);

  //------------------------------------------
  // Log column heads
  //------------------------------------------
  switch (purpose) {
  case purpose_Test:
    gw_log_message(testhead); printf("%s\n", testhead);
    break;
  case purpose_ListResults:
    gw_log_message(listhead1); printf("%s\n", listhead1);
    gw_log_message(listhead2); printf("%s\n", listhead2);
    break;
  } // switch

    //------------------------------------------
    // Randomize -- may want this or not
    //------------------------------------------
  srand((unsigned int) time((time_t *) NULL) );

  //------------------------------------------
  // Call test routine.
  // Note: change Purpose for different behavior!
  //------------------------------------------
  for (n = 0; n < TestRepetitions; n++) {
    rslt = TestNearestVtxInConcentricIcos
      (n,
       purpose);  // purpose_Test or purpose_ListResults
    if (rslt && (purpose == purpose_Test)) goto done;
  }

 done:

  gw_log_end();

  return rslt;
}

