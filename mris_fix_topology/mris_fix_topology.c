
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

static char vcid[] = "$Id: mris_fix_topology.c,v 1.3 1999/04/20 14:44:13 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;


static char *sphere_name = "qsphere" ;
static char *orig_name = "smoothwm" ;

#define MAX_VERTICES  0
#define MAX_FACES     0

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *sname, sdir[400], *cp, fname[500] ;
  int           ac, nargs ;
  MRI_SURFACE   *mris, *mris_corrected ;
  int           msec, nvert, nfaces, nedges, eno ;
  struct timeb  then ;

  Gdiag |= DIAG_WRITE ;
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

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, sphere_name) ;
  fprintf(stderr, "reading input surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read input surface %s",
              Progname, fname) ;
  strcpy(mris->subject_name, sname) ;
  eno = MRIScomputeEulerNumber(mris, &nvert, &nfaces, &nedges) ;
  fprintf(stderr, "before topology correction, eno=%d (nv=%d, nf=%d, ne=%d,"
          " g=%d)\n", eno, nvert, nfaces, nedges, (2-eno)/2) ;
  MRISprojectOntoSphere(mris, mris, 100.0f) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MRISreadOriginalProperties(mris, orig_name) ;

  fprintf(stderr, "using quasi-homeomorphic spherical map to tessellate "
          "cortical surface...\n") ;

  mris_corrected = MRIScorrectTopology(mris, NULL) ;
  MRISfree(&mris) ;
  eno = MRIScomputeEulerNumber(mris_corrected, &nvert, &nfaces, &nedges) ;
  fprintf(stderr, "after topology correction, eno=%d (nv=%d, nf=%d, ne=%d,"
          " g=%d)\n", eno, nvert, nfaces, nedges, (2-eno)/2) ;

  if (!mris_corrected)  /* for now */
    exit(0) ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "inflated_corrected");
  fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
  MRISwrite(mris_corrected, fname) ;

  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "smoothwm_corrected");
  fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
  MRISwrite(mris_corrected, fname) ;

  MRISrestoreVertexPositions(mris_corrected, CANONICAL_VERTICES) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "sphere_corrected");
  fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
  MRISwrite(mris_corrected, fname) ;

#if 0
  fprintf(stderr, "computing curvature of regenerated surface...\n") ;
  MRISsetNeighborhoodSize(mris_corrected, 2) ;
  MRIScomputeSecondFundamentalForm(mris_corrected) ;
  MRISuseMeanCurvature(mris_corrected) ;
  MRISaverageCurvatures(mris_corrected, 10) ;
  MRISwriteCurvature(mris_corrected, "curv_corrected") ;
#endif

/*
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "ico_geo") ;
  fprintf(stderr, "writing output surface to %s...\n", fname) ;
  MRISwrite(mris_corrected, fname) ;
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
  else if (!stricmp(option, "name"))
  {
    sphere_name = argv[2] ;
    fprintf(stderr,"reading spherical homeomorphism from '%s'\n",sphere_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "orig"))
  {
    orig_name = argv[2] ;
    fprintf(stderr,"reading original coordinates from '%s'\n",orig_name);
    nargs = 1 ;
  }
  else switch (toupper(*option))
  {
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
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




