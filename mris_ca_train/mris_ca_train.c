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

static int ninputs = 2 ;  /* curv and sulc */
static int icno = 7 ;

static char *curv_name = "curv" ;
static char *thickness_name = "thickness" ;
static char *sulc_name = "sulc" ;

static char subjects_dir[STRLEN] ;

int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp,*hemi,
               *canon_surf_name, *annot_name ;
  int          ac, nargs, i, train_type ;
  int          msec, minutes, seconds, nsubjects ;
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
  if (argc < 6)
    usage_exit(1) ;

  hemi = argv[1] ; canon_surf_name = argv[2] ; annot_name = argv[3] ;
  out_fname = argv[argc-1] ;
  nsubjects = argc-5 ;

  gcsa = GCSAalloc(ninputs, icno) ;

  for (train_type = 0 ; train_type <= 1 ; train_type++)
  {
    printf("computing %s for %d subject and writing results to %s\n",
           train_type ? "covariances" : "means", nsubjects, out_fname) ;
    for (i = 0 ; i < nsubjects ; i++)
    {
      subject_name = argv[i+4] ;
      printf("processing subject %s, %d of %d...\n", subject_name,i+1,
             nsubjects);
      sprintf(fname, "%s/%s/surf/%s.%s", subjects_dir, subject_name, 
              hemi, orig_name) ;
      if (DIAG_VERBOSE_ON)
        printf("reading surface from %s...\n", fname) ;
      mris = MRISread(fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s for %s",
                  Progname, fname, subject_name) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
      if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read annot file %s for %s",
                  Progname, annot_name, subject_name) ;

      if (MRISreadCanonicalCoordinates(mris, canon_surf_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read annot file %s for %s",
                  Progname, canon_surf_name, subject_name) ;
      if (ninputs > 2)
      {
        if (MRISreadCurvature(mris, thickness_name) != NO_ERROR)
          ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                    Progname, thickness_name, subject_name) ;
        MRIScopyCurvatureToImagValues(mris) ;
      }
      if (ninputs > 1)
      {
        if (MRISreadCurvature(mris, sulc_name) != NO_ERROR)
          ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                    Progname, sulc_name, subject_name) ;
        MRIScopyCurvatureToValues(mris) ;
        MRIScopyValToVal2(mris) ;
      }
      if (MRISreadCurvature(mris, curv_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                  Progname, curv_name, subject_name) ;
      MRIScopyCurvatureToValues(mris) ;

      MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
      MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
      if (train_type == 0)
        GCSAtrainMeans(gcsa, mris) ;
      else
        GCSAtrainCovariances(gcsa, mris) ;
      MRISfree(&mris) ;
    }
    if (train_type == 0)
      GCSAnormalizeMeans(gcsa) ;
    else
      GCSAnormalizeCovariances(gcsa) ;
  }

  printf("writing classifier array to %s...\n", out_fname) ;
  GCSAwrite(gcsa, out_fname) ;
  GCSAfree(&gcsa) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("classifier array training took %d minutes"
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
  else if (!stricmp(option, "IC"))
  {
    icno = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using ico # %d for classifier array\n", icno) ;
  }
  else switch (toupper(*option))
  {
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'N':
    ninputs = atoi(argv[2]) ;
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
  printf("usage: %s [options] <hemi> <canon surf> <annot file> <subject 1> <subject 2> ... <output file>\n",
         Progname) ;
  exit(code) ;
}
