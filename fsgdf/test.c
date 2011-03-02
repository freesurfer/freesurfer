/**
 * @file  test.c
 * @brief test fsgd stuff using sample data
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:08 $
 *    $Revision: 1.14 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <tcl.h>
#include "fsgdf.h"
#include "fsgdf_wrap.h"

#define Assert(x,s) \
    if( !(x) ) { \
      fprintf (stderr, "Error: %s\n", s ); \
      exit(1); \
    } \

#define AssertTclOK(x,s) \
    if( TCL_OK != (x) ) { \
      fprintf (stderr, "Tcl_Eval returned not TCL_OK:\n"  \
         "Trying: %s\nResult: %s\n", s, interp->result ); \
      exit(1); \
    } \

char *Progname="fsgdf_test_app";

int main (int argc, char** argv) {
  FSGD *gd;
  char *env;
  char fnTestPath[1000];
  char fnTest[1000];
  int  rTcl;

  if ( argc > 1 ) {

    /* Grab the file name from the arg */
    strcpy(fnTest, argv[1]);

  } else {

    /* Build a proper path for our test data. */
    env = getenv("FSDEV_TEST_DATA");
    if (NULL != env) {
      strcpy(fnTestPath, env);
    } else {
      strcpy(fnTestPath, "./test_data");
    }
    sprintf(fnTest, "%s/fsgdf/y-lh.fsgd", fnTestPath);
  }

  /* Just test the read function, first without reading the data and
     then with. */
  gd = gdfRead(fnTest,0);
  if (NULL == gd) {
    printf("ERROR: gdfRead(%s,0) returned NULL",fnTest);
    exit(1);
  }
  gdfFree(&gd);
  gd = gdfRead(fnTest,1);
  if (NULL == gd) {
    printf("ERROR: gdfRead(%s,1) returned NULL",fnTest);
    exit(1);
  }
  gdfFree(&gd);

  printf("C Test successful.\n\n");

  printf("Tcl_CreateInterp..."); fflush(stdout);
  Tcl_Interp* interp = Tcl_CreateInterp();
  Assert( NULL != interp, "Tcl_CreateInterp returned null" );
  printf("done.\n");

  printf("Tcl_Init..."); fflush(stdout);
  rTcl = Tcl_Init( interp );
  AssertTclOK( rTcl, "Tcl_Init returned not TCL_OK" );
  printf("done.\n");

  /* Initialize our Fsgdf functions. This is in fsgdf_wrap.c */
  printf("Fsgdf_Init..."); fflush(stdout);
  rTcl = Fsgdf_Init( interp );
  AssertTclOK( rTcl, "Fsgdf_Init failed" );
  printf("done.\n");

  printf("Tcl_EvalFile( test.tcl )..."); fflush(stdout);
  rTcl = Tcl_EvalFile( interp, "./test.tcl" );
  if( TCL_OK != rTcl ) {
    fprintf (stderr, "Tcl_EvalFile( test.tcl ) returned not TCL_OK:\n" 
             "rTcl=%d\nResult: %s\n", rTcl, interp->result );
    // don't exit 1, as an error will return if a display is not available
  }
  printf("done.\n");

  printf("Tcl Test successful.\n\n");

  return(0);
}
