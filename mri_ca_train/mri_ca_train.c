#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "gca.h"
#include "transform.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

static GCA_PARMS parms ;
static char *parc_dir = "parc" ;
static char *orig_dir = "orig" ;
static char *xform_name = "talairach.xfm" ;

static int ninputs = 1 ;  /* T1 intensity */

static char subjects_dir[STRLEN] ;
static char *heq_fname = NULL ;

int
main(int argc, char *argv[])
{
  char   **av, fname[STRLEN], *out_fname, *subject_name, *cp ;
  int    ac, nargs, i ;
  int          msec, minutes, seconds, nsubjects ;
  struct timeb start ;
  GCA          *gca ;
  MRI          *mri_parc, *mri_T1, *mri_eq = NULL ;
  LTA          *lta ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  parms.use_gradient = 0 ;
  parms.spacing = 3.0f ;

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
    if (argc < 3)
      usage_exit(1) ;
  }


  if (heq_fname)
  {
    mri_eq = MRIread(heq_fname) ;
    if (!mri_eq)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not read histogram equalization volume %s", 
                Progname, heq_fname) ;
  }

  out_fname = argv[argc-1] ;
  nsubjects = argc-2 ;
  fprintf(stderr, "training on %d subject and writing results to %s\n",
          nsubjects, out_fname) ;

  gca = GCAalloc(ninputs, parms.spacing, 
                 DEFAULT_VOLUME_SIZE, DEFAULT_VOLUME_SIZE,DEFAULT_VOLUME_SIZE);

  for (i = 0 ; i < nsubjects ; i++)
  {
    subject_name = argv[i+1] ;
    printf("processing subject %s, %d of %d...\n", subject_name,i+1,nsubjects);
    sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, parc_dir) ;
    if (DIAG_VERBOSE_ON)
      fprintf(stderr, "reading parcellation from %s...\n", fname) ;
    mri_parc = MRIread(fname) ;
    if (!mri_parc)
      ErrorExit(ERROR_NOFILE, "%s: could not read parcellation file %s",
                Progname, fname) ;

    sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, orig_dir) ;
    if (DIAG_VERBOSE_ON)
      fprintf(stderr, "reading co-registered T1 from %s...\n", fname) ;
    mri_T1 = MRIread(fname) ;
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE, "%s: could not read T1 data from file %s",
                Progname, fname) ;

    if (mri_eq)
    {
      printf("histogram equalizing input image...\n") ;
      MRIhistoEqualize(mri_T1, mri_eq, mri_T1, 30, 170) ;
    }

    if (xform_name)
    {
      sprintf(fname, "%s/%s/mri/transforms/%s", 
              subjects_dir, subject_name, xform_name) ;
      if (DIAG_VERBOSE_ON)
        fprintf(stderr, "reading transform from %s...\n", fname) ;
      lta = LTAread(fname) ;
      if (!lta)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
                  Progname, fname) ;
    }
    else
      lta = LTAalloc(1, NULL) ;

    GCAtrain(gca, mri_T1, mri_parc, lta) ;
    MRIfree(&mri_parc) ; MRIfree(&mri_T1) ; LTAfree(&lta) ;
  }

  GCAcompleteTraining(gca) ;
  GCAwrite(gca, out_fname) ;
  GCAfree(&gca) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "classifier array training took %d minutes"
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
  if (!stricmp(option, "GRADIENT"))
  {
    parms.use_gradient = 1 ;
    ninputs += 3 ;  /* components of the gradient */
  }
  else if (!stricmp(option, "SPACING"))
  {
    parms.spacing = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "spacing nodes every %2.1f mm\n", parms.spacing) ;
  }
  else if (!stricmp(option, "HEQ"))
  {
    heq_fname = argv[2] ;
    nargs = 1 ;
    printf("reading template for histogram equalization from %s...\n", 
           heq_fname) ;
  }
  else if (!stricmp(option, "T1"))
  {
    orig_dir = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "reading T1 data from subject's mri/%s directory\n",
            orig_dir) ;
  }
  else if (!stricmp(option, "PARC_DIR"))
  {
    parc_dir = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "reading parcellation from subject's mri/%s directory\n",
            parc_dir) ;
  }
  else if (!stricmp(option, "XFORM"))
  {
    xform_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "reading xform from %s\n", xform_name) ;
  }
  else if (!stricmp(option, "NOXFORM"))
  {
    xform_name = NULL ;
    fprintf(stderr, "disabling application of xform...\n") ;
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using %s as subjects directory\n", subjects_dir) ;
  }
  else switch (toupper(*option))
  {
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
  printf("usage: %s [options] <subject 1> <subject 2> ... <output file>\n",
         Progname) ;
  printf(
         "\t-spacing  - spacing of classifiers in canonical space\n");
  printf("\t-gradient - use intensity gradient as input to classifier.\n") ;
  exit(code) ;
}
