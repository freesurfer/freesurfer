#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "mrisurf.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "gcsa.h"
#include "transform.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

static char *orig_name = "smoothwm" ;

static char *curv_name = "curv" ;
static char *thickness_name = "thickness" ;
static char *sulc_name = "sulc" ;

static char subjects_dir[STRLEN] ;

int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp,*hemi,
               *canon_surf_name ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MRI_SURFACE  *mris ;
  GCSA         *gcsa ;

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

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment", 
                Progname);
    strcpy(subjects_dir, cp) ;
  }
  if (argc < 5)
    usage_exit(1) ;

  subject_name = argv[1] ; hemi = argv[2] ; canon_surf_name = argv[3] ; 
  out_fname = argv[5] ;

  gcsa = GCSAread(argv[4]) ;
  if (!gcsa)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier from %s",
              Progname, argv[4]) ;

  sprintf(fname, "%s/%s/surf/%s.%s", subjects_dir, subject_name,hemi,orig_name);
  if (DIAG_VERBOSE_ON)
    printf("reading surface from %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s for %s",
              Progname, fname, subject_name) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  if (MRISreadCanonicalCoordinates(mris, canon_surf_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read annot file %s for %s",
              Progname, canon_surf_name, subject_name) ;
  if (gcsa->ninputs > 2)
  {
    if (MRISreadCurvature(mris, thickness_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                Progname, thickness_name, subject_name) ;
    MRIScopyCurvatureToImagValues(mris) ;
  }
  if (gcsa->ninputs > 1)
  {
    if (MRISreadCurvature(mris, sulc_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                Progname, curv_name, subject_name) ;
    MRIScopyCurvatureToValues(mris) ;
    MRIScopyValToVal2(mris) ;
  }
  
  if (MRISreadCurvature(mris, curv_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
              Progname, sulc_name, subject_name) ;
  MRIScopyCurvatureToValues(mris) ;
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

  GCSAlabel(gcsa, mris) ;

  printf("writing output to %s...\n", out_fname) ;
  if (MRISwriteAnnotation(mris, out_fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not write annot file %s for %s",
                  Progname, out_fname, subject_name) ;

  MRISfree(&mris) ;
  GCSAfree(&gcsa) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("classification took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "SDIR"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  }
  else if (!stricmp(option, "ORIG"))
  {
    orig_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as original surface\n", orig_name) ;
  }
  else switch (toupper(*option))
  {
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
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
  printf("usage: %s [options] <subject> <hemi> <canon surf> "
         "<classifier> <output file>\n", Progname) ;
  exit(code) ;
}
