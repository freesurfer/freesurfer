
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"

static char vcid[] = "$Id: mris_curvature_stats.c,v 1.4 2003/09/05 04:45:40 kteich Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int navgs = 0 ;
static int normalize_flag = 0 ;
static int condition_no = 0 ;
static int stat_flag = 0 ;
static char *label_name = NULL ;
static char *output_fname = NULL ;

#define MAX_FILES 1000
static double means[MAX_FILES] ;

int
main(int argc, char *argv[])
{
  char         **av, surf_name[STRLEN], fname[STRLEN], *sdir ;
  char         *hemi, *subject_name, *curv_fname ;
  int          ac, nargs, i ;
  MRI_SURFACE  *mris ;
  double       mean, sigma, n, total, total_sq ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_curvature_stats.c,v 1.4 2003/09/05 04:45:40 kteich Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  sdir = getenv("SUBJECTS_DIR") ;
  if (!sdir)
    ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in envoronment.\n",Progname);
  ac = argc ;
  av = argv ;
  strcpy(surf_name, "orig") ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  /* parameters are 
     1. subject name
     2. hemisphere
  */
  subject_name = argv[1] ; hemi = argv[2] ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, surf_name) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  if (label_name)
  {
    LABEL *area ;
    area = LabelRead(subject_name, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s",
                Progname, label_name) ;
    LabelRipRestOfSurface(area, mris) ;
    LabelFree(&area) ;
  }
#define START_i  3
  if (label_name)
    fprintf(stdout, "%s: ", label_name) ;
  for (n = total_sq = total = 0.0, i = START_i ; i < argc ; i++)
  {
    curv_fname = argv[i] ;
    if (MRISreadCurvatureFile(mris, curv_fname) != NO_ERROR)
      ErrorExit(ERROR_BADFILE,"%s: could not read curvature file %s.\n",
                Progname, curv_fname);
    if (normalize_flag)
      MRISnormalizeCurvature(mris) ;

    MRISaverageCurvatures(mris, navgs) ;
    mean = MRIScomputeAverageCurvature(mris, &sigma) ;
    fprintf(stdout, "curvature file %s: %2.2f += %2.2f\n", 
            curv_fname, mean, sigma) ;
    means[i-START_i] = mean ;
    total += mean ; total_sq += mean*mean ;
    n += 1.0 ;
  }

  if (n > 1.8)
  {
    mean = total / n ;
    sigma = sqrt(total_sq/n - mean*mean) ;
    fprintf(stdout, "mean across %d subjects: %2.3f +- %2.3f\n",
            (int)n, mean, sigma) ;
  }
  else
    mean = sigma = 0.0 ;
  

  MRISfree(&mris) ;
  if (output_fname)
  {
    FILE *fp ;

    fp = fopen(output_fname, "a") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open output file %s",
                Progname, output_fname) ;
    if (label_name)
      fprintf(fp, "%s: ", label_name) ;
    fprintf(fp, "%2.3f +- %2.3f\n", mean, sigma) ;
    fclose(fp) ;
  }
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
  case 'O':
    output_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "outputting results to %s...\n", output_fname) ;
    break ;
  case 'L':
    label_name = argv[2] ;
    fprintf(stderr, "using label %s...\n", label_name) ;
    nargs = 1 ;
    break ;
  case 'A':
    navgs = atoi(argv[2]) ;
    fprintf(stderr, "averaging curvature %d timesf\n", navgs);
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'S':   /* write out stats */
    stat_flag = 1 ;
    condition_no = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "writing out summary statistics as condition %d\n",
            condition_no) ;
    break ;
  case 'N':
    normalize_flag = 1 ;
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

static void
print_usage(void)
{
  print_help() ;
}

static void
print_help(void)
{
  fprintf(stderr, 
       "\nThis program will compute mean and variances for a set of"
          " curvature files, possibly constrained by a label.\n") ;
  fprintf(stderr, 
          "usage: %s [options] <subject name> <hemi> <curv file> "
          "<curv file> ... \n", Progname) ;
  fprintf(stderr, " curvature files, possibly constrained by a label.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-l <label name>  limit statistics to labelled region\n"
                  "-o <output file> append results to named output file\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

