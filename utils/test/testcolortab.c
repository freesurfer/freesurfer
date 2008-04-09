/**
 * @file  testcolortab.c
 * @brief check the FreeSurferColorLUT.txt for problems
 *
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/04/09 18:45:16 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2008,
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "colortab.h"

int main(int argc, char *argv[])
{
  int errs=ctabDuplicates=0; // CTABreadASCII will increment ctabDuplicates
  fprintf(stdout,"Colortable test...\n");
  char defaultfname[]="../../distribution/FreeSurferColorLUT.txt";
  char* fname=argv[1];
  if (NULL==fname) fname = defaultfname;
  COLOR_TABLE *ct = CTABreadASCII(fname);
  if (NULL==ct)
  {
    fprintf(stderr,"usage:   %s colortable_fname [1]\n",argv[0]);
    fprintf(stderr,"example: %s FreeSurferColorLUT.txt\n",argv[0]);
    fprintf
      (stderr,
       "default: %s ../../distribution/FreeSurferColorLUT.txt\n",
       argv[0]);
    fprintf
      (stderr,
       "option:  if 1 is the second arg, look for duplicate annotations\n");
    exit(1);
  }
  errs=ctabDuplicates; // CTABreadASCII puts duplicate count in ctabDuplicates
  errs+=CTABfindDuplicateNames(ct);
  if ((NULL!=argv[2]) && (0==strcmp(argv[2],"1"))) errs+=CTABfindDuplicateAnnotations(ct);
  if (errs)
  {
    fflush(stderr);
    fflush(stdout);
    fprintf(stderr,"\n!!! FAILURE !!! found %d errors in %s\n",errs,fname);
  }
  return errs;
}
