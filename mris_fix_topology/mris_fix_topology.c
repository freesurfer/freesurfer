
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
#include "icosahedron.h"
#include "mrishash.h"

static char vcid[] = "$Id: mris_fix_topology.c,v 1.2 1999/03/01 00:04:14 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int small = 0 ;
static int huge = 1 ;

static INTEGRATION_PARMS parms ;
static int avgs = 10 ;

#define MAX_VERTICES  0
#define MAX_FACES     0

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *sname, sdir[400], *cp, fname[500] ;
  int           ac, nargs ;
  MRI_SURFACE   *mris, *mris_ico ;
  int           msec ;
  struct timeb  then ;
  float         radius ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.l_tspring = 1 ; parms.l_curv = .0 ; parms.niterations = 0 ;
  parms.tol = 1e-4 ; parms.dt = 0.1f ;
  parms.projection = NO_PROJECTION ;
  parms.integration_type = INTEGRATE_MOMENTUM ;

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

  if (argc < 2)
    usage_exit() ;

  TimerStart(&then) ;
  sname = argv[1] ;
  hemi = argv[2] ;
  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, 
              "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
  strcpy(sdir, cp) ;

  fprintf(stderr, "creating icosahedral representation...\n") ;

  /*  if (huge)*/
    mris_ico = ic163842_make_surface(MAX_VERTICES, MAX_FACES) ;
#if 0
  else
  {
    if (small)
      mris_ico = ic10242_make_surface(MAX_VERTICES, MAX_FACES) ;
    else
      mris_ico = ic40962_make_surface(MAX_VERTICES, MAX_FACES) ;
  }
#endif
  if (!mris_ico)
    ErrorExit(ERROR_NOFILE, "%s: could not create/read MR surface.",Progname) ;
  mris_ico->hemisphere = stricmp(hemi, "lh") ? 
      RIGHT_HEMISPHERE : LEFT_HEMISPHERE ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "qsphere") ;
  fprintf(stderr, "reading input surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read input surface %s",
              Progname, fname) ;
  MRIScomputeMetricProperties(mris) ;
  MRISreadOriginalProperties(mris, "orig") ;
/*
  MRISreadCanonicalCoordinates(mris, "inflated") ;
*/
/*
  MRISreadCanonicalCoordinates(mris, "sphere") ;
*/
  /* save orig. vertex coords */
  radius = MRISaverageRadius(mris_ico) ;
  MRISscaleBrain(mris_ico, mris_ico, 100.0f/radius) ;
  MRISsaveVertexPositions(mris_ico, ORIGINAL_VERTICES) ;

  fprintf(stderr, "mapping icosahedron to cortical surface...\n") ;

  MRISreadOriginalProperties(mris, "smoothwm") ;
  MRISinverseSphericalMap(mris, mris_ico) ;
  MRISfree(&mris) ;
  MRISrestoreVertexPositions(mris_ico, ORIGINAL_VERTICES) ;

#if 0
  MRISaverageVertexPositions(mris_ico, 5) ;
#else
  {
    
    sprintf(mris_ico->fname, 
            "%s/%s/surf/%s.%s", sdir, sname, hemi, "smoothwm_ico") ;
    /*    MRISsetNeighborhoodSize(mris_ico, 2) ;*/
    
    MRIScomputeMetricProperties(mris_ico) ;
    strcpy(parms.base_name, "smoothwm_ico") ;
    mris_ico->status = MRIS_SURFACE ;  /* no longer a sphere */
    MRISstoreMetricProperties(mris_ico) ;
    MRISprintTessellationStats(mris_ico, stderr) ;
    if (parms.niterations > 0)
      MRISintegrate(mris_ico, &parms, 0) ;
  }
#endif
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "smoothwm_ico") ;
  fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
  MRISwrite(mris_ico, fname) ;

  fprintf(stderr, "computing curvature of regenerated surface...\n") ;
  MRISsetNeighborhoodSize(mris_ico, 2) ;
  MRIScomputeSecondFundamentalForm(mris_ico) ;
  MRISuseMeanCurvature(mris_ico) ;
  MRISaverageCurvatures(mris_ico, 10) ;
  MRISwriteCurvature(mris_ico, "curv_ico") ;

/*
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "ico_geo") ;
  fprintf(stderr, "writing output surface to %s...\n", fname) ;
  MRISwrite(mris_ico, fname) ;
*/

  msec = TimerStop(&then) ;
  fprintf(stderr,"topology fixing took %2.1f minutes\n", 
          (float)msec/(60*1000.0f));
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
  case 'W':
    Gdiag |= DIAG_WRITE ;
    sscanf(argv[2], "%d", &parms.write_iterations) ;
    nargs = 1 ;
    fprintf(stderr, "using write iterations = %d\n", parms.write_iterations) ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &parms.niterations) ;
    nargs = 1 ;
    fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
    break ;
  case 'M':
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", (float)parms.momentum) ;
    break ;
  case 'A':
    avgs = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'S':
    small = 1 ;
    break ;
  case 'H':
    huge = 1 ;
    break ;
  case 'B':
    small = 0 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
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
  fprintf(stderr, "usage: %s [options] <subject name> <hemisphere>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program computes a mapping from the unit sphere onto the\n"
          "surface of the cortex from a previously generated approximation\n"
          "of the cortical surface, thus guaranteeing a topologically\n"
          "correct surface.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}




