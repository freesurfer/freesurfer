
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
#include "version.h"

static char vcid[] = "$Id: mris_fix_topology.c,v 1.18 2004/05/19 17:19:48 segonne Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int exit_after_diag = 0 ;

static char *T1_name = "T1" ;
static char *wm_name = "wm" ;
static char *sphere_name = "qsphere" ;
static char *inflated_name = "inflated" ;
static char *orig_name = "orig" ;
static char suffix[STRLEN] = "" ;
static int  add = 1 ;
static int  write_inflated = 0 ;
static int nsmooth = 5 ;

static char sdir[STRLEN] = "" ;
static TOPOLOGY_PARMS parms ;

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *sname, *cp, fname[STRLEN] ;
  int           ac, nargs ;
  MRI_SURFACE   *mris, *mris_corrected ;
  MRI           *mri, *mri_wm ;
  int           msec, nvert, nfaces, nedges, eno ;
  float         max_len ;
  struct timeb  then ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_fix_topology.c,v 1.18 2004/05/19 17:19:48 segonne Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

	parms.max_patches = 10 ;
	parms.max_unchanged = 10 ;
	parms.l_mri = 1 ;
	parms.l_curv = 1 ;
	parms.l_unmri = 1 ;
	parms.niters=-1;
	parms.genetic=0;

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
  if (strlen(sdir) == 0)
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

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

  sprintf(fname, "%s/%s/mri/%s", sdir, sname, T1_name) ;
  printf("reading T1 volume from %s...\n", T1_name) ;
  mri = MRIread(fname) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read T1 volume from %s", Progname, fname) ;

  sprintf(fname, "%s/%s/mri/%s", sdir, sname, wm_name) ;
  printf("reading wm segmentation from %s...\n", wm_name) ;
  mri_wm = MRIread(fname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read T1 volume from %s", Progname, fname) ;

  if (MRISreadOriginalProperties(mris, orig_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read original surface %s",
              Progname, orig_name) ;

  if (MRISreadVertexPositions(mris, inflated_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read inflated surface %s",
              Progname, inflated_name) ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;

  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;

  fprintf(stderr, "using quasi-homeomorphic spherical map to tessellate "
          "cortical surface...\n") ;

  mris_corrected = MRIScorrectTopology(mris, NULL, mri, mri_wm, nsmooth, &parms) ;
  MRISfree(&mris) ;
  eno = MRIScomputeEulerNumber(mris_corrected, &nvert, &nfaces, &nedges) ;
  fprintf(stderr, "after topology correction, eno=%d (nv=%d, nf=%d, ne=%d,"
          " g=%d)\n", eno, nvert, nfaces, nedges, (2-eno)/2) ;

  if (!mris_corrected || exit_after_diag)  /* for now */
    exit(0) ;

  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
  if (add)
    for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
      while (MRISdivideLongEdges(mris_corrected, max_len) > 0)
      {}

  if (write_inflated)
  {
    MRISrestoreVertexPositions(mris_corrected, TMP_VERTICES) ;
    sprintf(fname, "%s/%s/surf/%s.%s%s", sdir,sname,hemi,inflated_name,suffix);
    fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
    MRISwrite(mris_corrected, fname) ;
  }

  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
  sprintf(fname, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, orig_name,suffix);
  fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
  MRISwrite(mris_corrected, fname) ;

  /* compute the orientation changes */
  MRISmarkOrientationChanges(mris_corrected);
  


#if 0
  MRISrestoreVertexPositions(mris_corrected, CANONICAL_VERTICES) ;
  sprintf(fname, "%s/%s/surf/%s.sphere%s", sdir, sname, hemi, suffix);
  fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
  MRISwrite(mris_corrected, fname) ;
#endif

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
  else if (!stricmp(option, "name") || !stricmp(option, "sphere"))
  {
    sphere_name = argv[2] ;
    fprintf(stderr,"reading spherical homeomorphism from '%s'\n",sphere_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "inflated"))
  {
    inflated_name = argv[2] ;
    fprintf(stderr,"reading inflated coordinates from '%s'\n",inflated_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "mri"))
  {
    parms.l_mri = atof(argv[2]) ;
    fprintf(stderr,"setting l_mri = %2.2f\n", parms.l_mri) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "curv"))
  {
    parms.l_curv = atof(argv[2]) ;
    fprintf(stderr,"setting l_curv = %2.2f\n", parms.l_curv) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "unmri"))
  {
    parms.l_unmri = atof(argv[2]) ;
    fprintf(stderr,"setting l_unmri = %2.2f\n", parms.l_unmri) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "patches"))
  {
    parms.max_patches = atoi(argv[2]) ;
    fprintf(stderr,"using %d defect patches/generation...\n", parms.max_patches) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "generations"))
  {
    parms.max_unchanged = atoi(argv[2]) ;
    fprintf(stderr,"terminating evolution after %d generations without change...\n", parms.max_unchanged) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "niters"))
  {
    parms.niters = atoi(argv[2]) ;
    fprintf(stderr,"stopping genetic algorithm after %d iterations\n", parms.niters) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "genetic"))
  {
    parms.genetic=1;
    fprintf(stderr,"using genetic algorithm\n") ;
    nargs = 0 ;
  }
  else if (!stricmp(option, "diag"))
  {
    printf("saving diagnostic information....\n") ;
    Gdiag |= DIAG_SAVE_DIAGS ;
  }
  else if (!stricmp(option, "diagonly"))
  {
    printf("saving diagnostic information and exiting....\n") ;
    Gdiag |= DIAG_SAVE_DIAGS ;
    exit_after_diag = 1 ;
  }
  else if (!stricmp(option, "sdir"))
  {
    strcpy(sdir, argv[2]) ;
    fprintf(stderr,"using %s as subjects_dir\n", sdir);
    nargs = 1 ;
  }
  else if (!stricmp(option, "wi"))
  {
    write_inflated = 1 ;
    fprintf(stderr,"writing fixed inflated coordinates to %s\n",inflated_name);
  }
  else if (!stricmp(option, "suffix"))
  {
    strcpy(suffix, argv[2]) ;
    fprintf(stderr,"adding suffix '%s' to output names\n", suffix) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "add"))
  {
    fprintf(stderr,"adding vertices after retessellation\n") ;
    add = 1 ;
  }
  else if (!stricmp(option, "noadd"))
  {
    fprintf(stderr,"not adding vertices after retessellation\n") ;
    add = 0 ;
  }
  else if (!stricmp(option, "orig"))
  {
    orig_name = argv[2] ;
    fprintf(stderr,"reading original coordinates from '%s'\n",orig_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "seed"))
  {
		setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number genererator to %d\n", atoi(argv[2])) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
  {
  case 'S':
    nsmooth = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing corrected surface for %d iterations\n", nsmooth) ;
    break ;
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




