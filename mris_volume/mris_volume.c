/**
 * @file  mris_volume.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/01/15 23:37:58 $
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
// mris_measure_volume.c
// compute the volume enclosed by a closed triangulated surface
// original author: Xiao Han
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: greve $
// Revision Date  : $Date: 2007/01/15 23:37:58 $
// Revision       : $Revision: 1.4 $
//
////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "matrix.h"
#include "transform.h"
#include "version.h"
#include "label.h"
#include "mrisutils.h"

#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y, v1->z-v0->z)

static int verbose = 0;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;

static void usage_exit(int code) ;

int main(int argc, char *argv[]) {
  char   **av, *in_fname;
  int    ac, nargs;
  MRIS    *mris;
  int    msec, minutes, seconds;
  struct timeb start ;
  double total_volume;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_volume.c,v 1.4 2007/01/15 23:37:58 greve Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    usage_exit(1) ;

  //  printf("command line parsing finished\n");

  /*** Read in the input surfaces ***/
  in_fname = argv[1] ;
  if (verbose) printf("reading %s...\n", in_fname) ;

  mris = MRISread(in_fname) ;
  if(mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s",
              Progname, in_fname) ;

  if(verbose) printf("surface file read in.\n");

  total_volume = MRISvolumeInSurf(mris);

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  if (verbose)
    printf("Volume computation took %d minutes and %d seconds.\n", 
	   minutes, seconds) ;

  if (verbose)
    printf("total volume surrounded by the surface is %g\n", total_volume);
  else
    printf("%lf\n", total_volume);

  MRISfree(&mris);

  exit(0);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  switch (*option) {
  case 'v':
  case 'V':
    verbose = 1;
    break;
  default:
    printf("unknown option %s\n", argv[1]);
    exit(1);
    break;
  }

  return(nargs) ;
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s surface_file_name\n", Progname) ;
  printf("\t This program computes the volume of the given closed surface using a divergence formula \n");
  printf("\t use -v option to output more messages\n");
  exit(code) ;
}



