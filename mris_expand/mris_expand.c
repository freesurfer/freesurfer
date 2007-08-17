/**
 * @file  mris_expand.c
 * @brief expand a surface outwards by a specified amount
=======
 * Expands a surface (typically ?h.white) outwards while maintaining smoothness
 * and self-intersection constraints.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/08/17 17:07:43 $
 *    $Revision: 1.8 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;
static INTEGRATION_PARMS parms ;
static int navgs = 0 ;
static int use_thickness = 0 ;

int
main(int argc, char *argv[])
{
  char         **av ;
  int          ac, nargs ;
  char         *in_fname, *out_fname ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  float        mm_out ;
  MRI_SURFACE  *mris ;

  parms.l_spring = .1;
  // parms.l_curv = 1.0 ;
  parms.n_averages = 4 ;
  // parms.l_surf_repulse = .1 ;
  parms.dt = 0.5 ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_expand.c,v 1.8 2007/08/17 17:07:43 fischl Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  in_fname = argv[1] ;
  mm_out = atof(argv[2]) ;
  out_fname = argv[3] ;
  FileNameExtension(out_fname,parms.base_name) ;  // remove hemi (e.g. lh.)

  if (use_thickness)
    printf("expanding surface %s by %2.1f%% of thickness and writing it to %s\n",
           in_fname, 100.0*mm_out, out_fname) ;
  else
    printf("expanding surface %s by %2.1f mm and writing it to %s\n",
           in_fname, mm_out, out_fname) ;
  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: MRISread(%s) failed", Progname, in_fname);
  if (use_thickness)
  {
    printf("reading thickness...\n") ;
    if (MRISreadCurvatureFile(mris, "thickness") != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not load thickness file", Progname) ;
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
    if (MRISreadVertexPositions(mris, "pial") != NO_ERROR)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read pial vertex positions\n",
                Progname) ;
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
  }
  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;
  MRISexpandSurface(mris, mm_out, &parms, use_thickness) ;
  if (navgs > 0)
    MRISaverageVertexPositions(mris, navgs) ;
  printf("writing expanded surface to %s...\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stdout, "surface expansion took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;

  exit(0) ;
  return(0) ;
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
  if (!stricmp(option, "thickness"))
  {
    use_thickness = 1 ;
    printf("using distance as a %% of thickness\n") ;
  }
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      printf("debugging vertex %d\n", Gdiag_no) ;
      break ;
    case 'S':
      parms.l_spring = atof(argv[2]) ;
      nargs = 1 ;
      printf("setting spring term to %2.2f\n", parms.l_spring) ;
      break ;
    case 'T':
      parms.dt = atof(argv[2]) ;
      nargs = 1 ;
      printf("setting dt = %2.2f\n", parms.dt) ;
      break ;
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("writing snapshots of expansion every %d iterations\n",
             parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      printf("smoothing surface with %d iterations after expansion\n", navgs) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code)
{
  printf("Usage: %s [options] <input surface> <mm> <output surface>\n",
         Progname) ;
  exit(code) ;
}
