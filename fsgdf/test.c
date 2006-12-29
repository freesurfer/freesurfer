/**
 * @file  test.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:57 $
 *    $Revision: 1.11 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <tcl.h>
#include "fsgdf.h"
#include "fsgdf_wrap.h"


#define Assert(x,s) \
    if( !(x) ) { \
      fprintf (stderr, "Error: %s\n", s ); \
    } \

#define AssertTclOK(x,s) \
    if( TCL_OK != (x) ) { \
      fprintf (stderr, "Tcl_Eval returned not TCL_OK:\n"  \
         "Trying: %s\nResult: %s\n", s, interp->result ); \
    } \


char *Progname="fsgdf test app";

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


  Tcl_Interp* interp = Tcl_CreateInterp();
  Assert( NULL != interp, "Tcl_CreateInterp returned null" );

  rTcl = Tcl_Init( interp );
  AssertTclOK( rTcl, "Tcl_Init returned not TCL_OK" );

  /* Initialize our Fsgdf functions. This is in fsgdf_wrap.c */
  rTcl = Fsgdf_Init( interp );
  AssertTclOK( rTcl, "Fsgdf_Init failed" );

  rTcl = Tcl_EvalFile( interp, "test.tcl" );
  AssertTclOK( rTcl, "Tcl_EvalFile returned not TCL_OK" );

  printf("Tcl Test successful.\n\n");

  return(0);
}
