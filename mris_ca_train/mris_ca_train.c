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
#include "version.h"

#define MAX_LABELS  1000
#if 0
static int write_ptable(char *fname, int *ptable, int nparcs) ;
#endif
static int find_parc_index(int parc, int *ptable, int nparcs) ;
static int add_to_ptable(MRI_SURFACE *mris, int *ptable, int nparcs) ;
static int *ptable = NULL ;
static int nbrs = 2 ;
static int navgs = 5 ;
static int normalize1_flag = 0 ;
static int normalize2_flag = 0 ;
static int normalize3_flag = 0 ;
static int nparcs = 0 ;
static char *ptable_fname = NULL ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

static char *orig_name = "smoothwm" ;

static int ninputs = 1 ;  /* curv and sulc */
static int icno_priors = 7 ;
static int icno_classifiers = 4 ;

#if 0
static char *curv_name = "curv" ;
#endif
static char *thickness_name = "thickness" ;
static char *sulc_name = "sulc" ;
static int sulconly = 0 ;

static char subjects_dir[STRLEN] ;

int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp, *hemi;
  char         *canon_surf_name, *annot_name ;
  int          ac, nargs, i, train_type ;
  int          msec, minutes, seconds, nsubjects, input1_flags;
  int          input2_flags, input3_flags ;
  struct timeb start ;
  MRI_SURFACE  *mris ;
  GCSA         *gcsa ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, 
  "$Id: mris_ca_train.c,v 1.11 2005/11/16 23:04:41 nicks Exp $", "$Name:  $");
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

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
    {
      cp = getenv("SUBJECTS_DIR") ;
      if (!cp)
        ErrorExit(ERROR_BADPARM, 
                  "%s: SUBJECTS_DIR not defined in environment", 
                  Progname);
      strcpy(subjects_dir, cp) ;
    }
  if (argc < 6)
    usage_exit(1) ;

  hemi = argv[1] ; 
  canon_surf_name = argv[2] ; 
  annot_name = argv[3] ;
  out_fname = argv[argc-1] ;
  nsubjects = argc-5 ;

  gcsa = GCSAalloc(ninputs, icno_priors, icno_classifiers) ;
  input1_flags = input2_flags = input3_flags = 0 ;
  if (normalize1_flag)
    input1_flags |= GCSA_NORMALIZE ;
  if (normalize2_flag)
    input2_flags |= GCSA_NORMALIZE ;
  if (normalize3_flag)
    input3_flags |= GCSA_NORMALIZE ;

  if (sulconly)
    {
      GCSAputInputType(gcsa, 
                       GCSA_INPUT_CURV_FILE, 
                       sulc_name, 
                       0, 
                       0,
                       input1_flags);
    }
  else
    {
      GCSAputInputType(gcsa, GCSA_INPUT_CURVATURE, "mean_curvature", 
                       navgs, input1_flags, 0) ;
      if (ninputs > 1)
        GCSAputInputType(gcsa, 
                         GCSA_INPUT_CURV_FILE, 
                         sulc_name,
                         0,
                         input2_flags,
                         1);
      if (ninputs > 2)
        GCSAputInputType(gcsa, 
                         GCSA_INPUT_CURV_FILE, 
                         thickness_name, 
                         0, 
                         input3_flags, 
                         2);
    }

  for (train_type = 0 ; train_type <= 1 ; train_type++)
    {
      printf("computing %s for %d subject \n",
             train_type ? "covariances" : "means", nsubjects) ;
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
            ErrorExit(ERROR_NOFILE, 
                      "%s: could not read surface file %s for %s",
                      Progname, fname, subject_name) ;
          MRISsetNeighborhoodSize(mris, nbrs) ;
          MRIScomputeSecondFundamentalForm(mris) ;
          MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
          if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
            ErrorExit(ERROR_NOFILE, 
                      "%s: could not read annot file %s for %s",
                      Progname, annot_name, subject_name) ;
          if (ptable)
            nparcs = add_to_ptable(mris, ptable, nparcs) ;

          if (MRISreadCanonicalCoordinates(mris, canon_surf_name) != NO_ERROR)
            ErrorExit(ERROR_NOFILE, 
                      "%s: could not read spherical "
                      "registration file %s for %s",
                      Progname, canon_surf_name, subject_name) ;
          if (ninputs > 2)
            {
              if (MRISreadCurvature(mris, thickness_name) != NO_ERROR)
                ErrorExit(ERROR_NOFILE, 
                          "%s: could not read curv file %s for %s",
                          Progname, thickness_name, subject_name) ;
              if (normalize3_flag)
                MRISnormalizeCurvature(mris) ;
              MRIScopyCurvatureToImagValues(mris) ;
            }
          if (ninputs > 1 || sulconly)
            {
              if (MRISreadCurvature(mris, sulc_name) != NO_ERROR)
                ErrorExit(ERROR_NOFILE, 
                          "%s: could not read curv file %s for %s",
                          Progname, sulc_name, subject_name) ;
              if (normalize2_flag || (sulconly && normalize1_flag))
                MRISnormalizeCurvature(mris) ;
              MRIScopyCurvatureToValues(mris) ;
              MRIScopyValToVal2(mris) ;
            }
          if (!sulconly)
            {
#if 0
              if (MRISreadCurvature(mris, curv_name) != NO_ERROR)
                ErrorExit(ERROR_NOFILE, 
                          "%s: could not read curv file %s for %s",
                          Progname, curv_name, subject_name) ;
#else
              MRISuseMeanCurvature(mris) ;
              MRISaverageCurvatures(mris, navgs) ;
              if (normalize1_flag)
                MRISnormalizeCurvature(mris) ;
#endif
              MRIScopyCurvatureToValues(mris) ;
            }

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

#if 0
  if (ptable)
    write_ptable(fname, ptable, nparcs) ;
#endif
  printf("writing classifier array to %s...\n", out_fname) ;
  gcsa->ptable_fname = ptable_fname ;
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
  else if (!stricmp(option, "nbrs"))
    {
      nbrs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
    }
  else if (!stricmp(option, "ORIG"))
    {
      orig_name = argv[2] ;
      nargs = 1 ;
      printf("using %s as original surface\n", orig_name) ;
    }
  else if (!stricmp(option, "NORM1"))
    {
      printf("normalizing input #1 after reading...\n") ;
      normalize1_flag = 1 ;
    }
  else if (!stricmp(option, "NORM2"))
    {
      printf("normalizing input #2 after reading...\n") ;
      normalize2_flag = 1 ;
    }
  else if (!stricmp(option, "NORM3"))
    {
      printf("normalizing input #3 after reading...\n") ;
      normalize3_flag = 1 ;
    }
  else if (!stricmp(option, "IC"))
    {
      icno_priors = atoi(argv[2]) ;
      icno_classifiers = atoi(argv[3]) ;
      nargs = 2 ;
      printf("using ico # %d for classifier array, and %d for priors\n", 
             icno_classifiers, icno_priors) ;
    }
  else if (!stricmp(option, "SULC") || !stricmp(option, "SULCONLY"))
    {
      printf("using sulc as only input...\n") ;
      sulconly = 1 ;
    }
  else switch (toupper(*option))
    {
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'T':
      ptable_fname = argv[2] ;
      nargs = 1 ;
      ptable = (int *)calloc(MAX_LABELS, sizeof(int)) ;
      break ;
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
  printf("usage: %s [options] <hemi> <canon surf> "
         "<annot file> <subject 1> <subject 2> ... <output file>\n",
         Progname) ;
  exit(code) ;
}
#if 0
static int
write_ptable(char *fname, int *ptable, int nparcs)
{
  FILE   *fp ;
  int    i ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE, 
                 "write_ptable(%s, %d): could not open file",
                 fname, nparcs)) ;

  for (i = 0 ; i < nparcs ; i++)
    fprintf(fp, "%d   %d\n", i, ptable[i]) ;
  fclose(fp) ;
  return(NO_ERROR) ;
}
#endif
static int
add_to_ptable(MRI_SURFACE *mris, int *ptable, int nparcs)
{
  int     vno, i ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      i = find_parc_index(v->annotation, ptable, nparcs) ;
      if (i < 0 || i >= nparcs)
        {
          ptable[i] = v->annotation ;
          i = nparcs++ ;
        }
    }
  return(nparcs) ;
}

static int
find_parc_index(int parc, int *ptable, int nparcs)
{
  int   i ;

  for (i = 0 ; i < nparcs ; i++)
    if (ptable[i] == parc)
      return(i) ;
  return(-1) ;
}

