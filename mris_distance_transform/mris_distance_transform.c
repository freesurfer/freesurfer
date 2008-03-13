/**
 * @file  mris_distance_transform.c
 * @brief program for computing a distance transform on the surface
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2008/03/13 15:23:37 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2004-2007,
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
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "label.h"
#include "MARS_DT_Boundary.h"

static char vcid[] = 
"$Id: mris_distance_transform.c,v 1.1 2008/03/13 15:23:37 fischl Exp $";

int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;


int
main(int argc, char *argv[])
{
  char          **av, *output_fname ;
  int           ac, nargs, msec, mode ;
  LABEL         *area ;
  MRI_SURFACE   *mris ;
  struct timeb  then ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mris_distance_transform.c,v 1.1 2008/03/13 15:23:37 fischl Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;


  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  TimerStart(&then) ;
  mris = MRISread(argv[1]) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s",
              Progname, argv[1]) ;
  area = LabelRead(NULL, argv[2]) ;
  if (area == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read label %s",
              Progname, argv[2]) ;

  if (stricmp(argv[3], "signed") == 0)
    mode = DTRANS_MODE_SIGNED ;
  else if (stricmp(argv[3], "unsigned") == 0)
    mode = DTRANS_MODE_UNSIGNED ;
  else if (stricmp(argv[3], "outside") == 0)
    mode = DTRANS_MODE_OUTSIDE ;
  else
  {
    print_usage() ;
    ErrorExit(ERROR_BADPARM, "unrecognized mode choice %s\n", argv[3]) ;
  }
  output_fname = argv[4] ;

  MRIScomputeMetricProperties(mris) ;
  MRISdistanceTransform(mris, area, mode) ;
  MRISwriteValues(mris, output_fname) ;

  msec = TimerStop(&then) ;
  fprintf(stderr,"distance transform took %2.1f minutes\n", (float)msec/(60*1000.0f));

  exit(0) ;
  return(0) ;  /* for ansi */
}


/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      print_help() ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  printf("usage: %s [options] <surface> <label> <mode> <output file>\n",
         Progname) ;
  printf("where mode is one of 'signed', 'unsigned', 'outside'\n") ;
  
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "This program computes the distance transform of a label on "
          "the surface\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "(see the source code!)\n") ;
  exit(1) ;
}


static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

