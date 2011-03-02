/**
 * @file  mris_remove_intersection.c
 * @brief removes surface intersections
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:32 $
 *    $Revision: 1.6 $
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "utils.h"
#include "timer.h"
#include "version.h"

static char vcid[]=
  "$Id: mris_remove_intersection.c,v 1.6 2011/03/02 00:04:32 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

int main(int argc, char *argv[])
{
  char         **av, *in_surf_fname, *out_fname ;
  int          ac, nargs, msec ;
  MRI_SURFACE  *mris ;
  struct timeb then ;

  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_remove_intersection.c,v 1.6 2011/03/02 00:04:32 nicks Exp $",
   "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mris_remove_intersection.c,v 1.6 2011/03/02 00:04:32 nicks Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Gdiag = DIAG_SHOW ;

  TimerStart(&then) ;
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

  if (argc < 3)
  {
    usage_exit() ;
  }

  in_surf_fname = argv[1] ;
  out_fname = argv[2] ;

  mris = MRISread(in_surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_surf_fname) ;

  MRISaddCommandLine(mris, cmdline) ;

  // MRISsetNeighborhoodSize(mris, 2) ;
  MRISremoveIntersections(mris) ;

  printf("writing corrected surface to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  msec = TimerStop(&then) ;
  fprintf(stderr, "intersection removal took %2.2f hours\n",
          (float)msec/(1000.0f*60.0f*60.0f));

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
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else switch (toupper(*option))
    {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
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

#include "mris_remove_intersection.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_remove_intersection_help_xml,
                mris_remove_intersection_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

