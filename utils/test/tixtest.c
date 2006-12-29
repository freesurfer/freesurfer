/**
 * @file  tixtest.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:46 $
 *    $Revision: 1.4 $
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


//
// tixtest.c
//

#include <stdio.h>
#include <stdlib.h>

#include <tcl.h>
#include <tk.h>
#include <tix.h>

extern int Tix_Init ( Tcl_Interp* interp );
extern int Itcl_Init(Tcl_Interp* interp);
extern int Itk_Init(Tcl_Interp* interp);

int main(int argc, char *argv[])
{
  int        eTcl                              = TCL_OK;
  Tcl_Interp *interp=0;

  interp = Tcl_CreateInterp();

  /* read tcl/tk internal startup scripts */
  eTcl = Tcl_Init( interp );
  if ( TCL_OK != eTcl )
  {
    fprintf(stderr, "Tcl_Init returned %d: %s\n", (int)eTcl, interp->result);
    return -1;
  }
  eTcl = Tk_Init(interp);
  if ( TCL_OK != eTcl )
  {
    fprintf(stderr, "Tcl_Init returned %d: %s\n", (int)eTcl, interp->result);
    return -1;
  }
  eTcl = Tix_Init( interp );
  if ( TCL_OK != eTcl )
  {
    fprintf(stderr,
            "Tix_Init returned %d: %s\n", (int)eTcl, interp->result);

    fprintf(stderr,
            "Try adding Itcl_Init() and Itk_Init() before Tix_Init()\n");
    eTcl = Itcl_Init(interp);
    if ( TCL_OK != eTcl )
    {
      fprintf(stderr,
              "Itlc_Init returned %d: %s\n", (int)eTcl, interp->result);
    }
    eTcl = Itk_Init(interp);
    if ( TCL_OK != eTcl )
    {
      fprintf(stderr,
              "Itk_Init returned %d: %s\n", (int)eTcl, interp->result);
    }
    eTcl = Tix_Init( interp );
    if ( TCL_OK != eTcl )
    {
      printf("even after all these tix initialization still failed\n");
    }
    else
      printf("worked fine with Itcl and Itk init before tix init\n");
  }
  else
    printf("no need to use Itcl and Itk\n");

  // cleanup
  Tcl_DeleteInterp(interp);

  return 0;
}
