/**
 * @brief Finds and corrects defects in surface topology.
 *
 * "Topology Correction of Subcortical Segmentation",
 * F. Segonne, E. Grimson and B. Fischl
 * (2003) MICCAI 2003.
 *
 * "Genetic Algorithm for the Topology Correction of Cortical Surfaces",
 * F. Segonne, E. Grimson, B. Fischl
 * (2005) IPMI pp393--405.
 */
/*
 * Original Author: Bruce Fischl
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
#include <sys/time.h>
#include <sys/resource.h>

#include "mri.h"
#include "mrisurf.h"
#include "mrisurf_project.h"
#include "mrishash.h"
#include "icosahedron.h"

#include "error.h"
#include "diag.h"
#include "macros.h"

#include "romp_support.h"
#include "tags.h"
#include "proto.h"
#include "timer.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int noint = 1 ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static void print_parameters(void);

const char *Progname ;

static int exit_after_diag = 0 ;
extern int topology_fixing_exit_after_diag ;

static const char *surf_dir = "surf" ;
static const char *brain_name    = "brain" ;
static const char *wm_name       = "wm" ;
static const char *sphere_name   = "qsphere.nofix" ;
static const char *inflated_name = "inflated.nofix" ;
static const char *orig_name     = "orig.nofix" ;
static const char *out_name      = "orig" ;
const char *defectbasename = "defect";

static char suffix[STRLEN] = "" ;
static int  add = 1 ;
static int  write_inflated = 0 ;
static int nsmooth = 5 ;
static int sphere_smooth = 5;

static char sdir[STRLEN] = "" ;
static TOPOLOGY_PARMS parms ;
static int MGZ = 1; // set to 1 for MGZ

static double nVFMultiplier = 1.1;

int
main(int argc, char *argv[])
{
  ROMP_main
  
  char          **av, *hemi, *sname, *cp, fname[STRLEN] ;
  int           ac, nargs ;
  MRI_SURFACE   *mris, *mris_corrected ;
  MRI           *mri, *mri_wm ;
  int           msec, nvert, nfaces, nedges, eno ;
  float         max_len ;
  Timer then ;


  std::string cmdline = getAllInfo(argc, argv, "mris_fix_topology");

  nargs = handleVersionOption(argc, argv, "mris_fix_topology");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  parms.max_patches = 10 ;
  parms.max_unchanged = 10 ;
  parms.l_mri = 1 ;
  parms.l_curv = 1 ;
  parms.l_qcurv = 1;
  parms.l_unmri = 1 ;
  // default is greedy retessellation
  parms.search_mode=GREEDY_SEARCH;
  // default is to use all vertices
  parms.retessellation_mode=0;
  // default is not to kill vertices
  parms.vertex_eliminate=0;
  //select an initial population of chromosomes
  parms.initial_selection=0;
  //stopping the algo after niters (if niters>0)
  parms.niters=-1;
  //volume resolution for the lodunlikelihood
  parms.volume_resolution=-1;
  //keep all vertices or not
  parms.keep=0;
  //smooth every patch
  parms.smooth=0;
  //match patch onto surface using local intensities
  parms.match=0;
  // don't use the edge table
  parms.edge_table=0;

  //don't save
  parms.save_fname=NULL;
  parms.defect_number=-1;
  //verbose mode
  parms.verbose=0;
  //movie mode
  parms.movie=0;
  // which defect
  parms.correct_defect=-1;
  // self-intersection
  parms.check_surface_intersection=0;
  // use initial mapping only
  parms.optimal_mapping=0;

  //Gdiag |= DIAG_WRITE ;
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
  {
    usage_exit() ;
  }

  print_parameters();

  printf("%s\n",getVersion().c_str());
  printf("  %s\n",getVersion().c_str());
  fflush(stdout);

  then.reset() ;
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

  int req = snprintf(fname, STRLEN, "%s/%s/%s/%s.%s", sdir, sname, surf_dir, hemi, sphere_name) ;
  if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("reading input surface %s...\n", fname) ;
  mris = MRISreadOverAlloc(fname,nVFMultiplier) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read input surface %s",
              Progname, fname) ;
  MRISaddCommandLine(mris, cmdline) ;
  strcpy(mris->subject_name, sname) ;
  eno = MRIScomputeEulerNumber(mris, &nvert, &nfaces, &nedges) ;
  fprintf(stderr, "before topology correction, eno=%d (nv=%d, nf=%d, ne=%d,"
          " g=%d)\n", eno, nvert, nfaces, nedges, (2-eno)/2) ;

  MRISprojectOntoSphere(mris, mris, 100.0f) ;

  MRISsmoothOnSphere(mris,sphere_smooth);

  /* at this point : canonical vertices */
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, brain_name) ;
  if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ)
  {
    strcat(fname,".mgz");
  }
  printf("reading brain volume from %s...\n", brain_name) ;
  mri = MRIread(fname) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read brain volume from %s", Progname, fname) ;

  {
    float fmin, fmax ;
    MRIvalRange(mri, &fmin, &fmax) ;
    if (fmax > 255)
    {
      printf("scaling intensities down by %2.3f\n", 255.0/fmax) ;
      MRIscalarMul(mri, mri, 255/fmax) ;
    }
  }
  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, wm_name) ;
  if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ)
  {
    strcat(fname,".mgz");
  }
  printf("reading wm segmentation from %s...\n", wm_name) ;
  mri_wm = MRIread(fname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read wm volume from %s", Progname, fname) ;

  printf("Reading original properties of %s\n",orig_name);
  if (MRISreadOriginalProperties(mris, orig_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read original surface %s",
              Progname, orig_name) ;

  /* at this point : canonical vertices but orig are in orignal vertices */
  printf("Reading vertex positions of %s\n",inflated_name);
  if (MRISreadVertexPositions(mris, inflated_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read inflated surface %s",
              Progname, inflated_name) ;

  /* at this point : inflated vertices */
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;

  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  /* at this point : canonical vertices */
  MRIScomputeMetricProperties(mris) ;

  fprintf(stderr, "using quasi-homeomorphic spherical map to tessellate "
          "cortical surface...\n") ;

  mris_corrected =
    MRIScorrectTopology(mris, mri, mri_wm, nsmooth, &parms, defectbasename) ;
  /* at this point : original vertices
     real solution in original vertices = corrected smoothed orig vertices */
  MRISfree(&mris) ;

  if (!mris_corrected || exit_after_diag)  /* for now */
    exit(0) ;

  if (noint)
    MRISremoveIntersections(mris_corrected) ;

  eno = MRIScomputeEulerNumber(mris_corrected, &nvert, &nfaces, &nedges) ;
  fprintf(stderr, "after topology correction, eno=%d (nv=%d, nf=%d, ne=%d,"
          " g=%d)\n", eno, nvert, nfaces, nedges, (2-eno)/2) ;

  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
  /* at this point : smoothed corrected orig vertices = solution */
  if (add)
    for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
      while (MRISdivideLongEdges(mris_corrected, max_len) > 0) {}

  MRIScopyVolGeomFromMRI(mris_corrected, mri_wm) ;
  if (write_inflated)
  {
    MRISrestoreVertexPositions(mris_corrected, TMP_VERTICES) ;
    int req = snprintf(fname, STRLEN, "%s/%s/%s/%s.%s%s",
            sdir,sname,surf_dir, hemi,inflated_name,suffix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
    MRISwrite(mris_corrected, fname) ;
  }

  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES) ;
  /* at this point : smoothed corrected orig vertices = solution */
  req = snprintf(fname, STRLEN, "%s/%s/%s/%s.%s%s", sdir, sname, surf_dir, hemi, out_name,suffix);
  if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "writing corrected surface to %s...\n", fname) ;
  MRISwrite(mris_corrected, fname) ;

  /* compute the orientation changes */
  MRISmarkOrientationChanges(mris_corrected);

  /*
    int req = snprintf(fname, STRLEN, "%s/%s/%s/%s.%s", sdir, sname, surf_dir, hemi, "ico_geo") ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "writing output surface to %s...\n", fname) ;
    MRISwrite(mris_corrected, fname) ;
  */

  msec = then.milliseconds() ;
  printf("topology fixing took %2.1f minutes\n", (float)msec/(60*1000.0f));

  // Output formatted so it can be easily grepped
#ifdef HAVE_OPENMP
  int n_omp_threads = omp_get_max_threads();
  printf("FSRUNTIME@ mris_fix_topology %s %7.4f hours %d threads\n",hemi,msec/(1000.0*60.0*60.0),n_omp_threads);
#else
  printf("FSRUNTIME@ mris_fix_topology %s %7.4f hours %d threads\n",hemi,msec/(1000.0*60.0*60.0),1);
#endif

  printf("#VMPC# mris_fix_topology VmPeak  %d\n",GetVmPeak());

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
  else if (!stricmp(option, "name") || !stricmp(option, "sphere"))
  {
    sphere_name = argv[2] ;
    fprintf(stderr,"reading spherical homeomorphism from '%s'\n",
            sphere_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "inflated"))
  {
    inflated_name = argv[2] ;
    fprintf(stderr,"reading inflated coordinates from '%s'\n",inflated_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "wm"))
  {
    wm_name = argv[2] ;
    fprintf(stderr,"reading wm segmentation from '%s'\n",wm_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "verbose"))
  {
    parms.verbose=VERBOSE_MODE_DEFAULT;
    fprintf(stderr,"verbose mode on (default mode)\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "sphere_smooth"))
  {
    sphere_smooth=atoi(argv[2]);
    fprintf
    (stderr,
     "smoothing spherical representation for %d iterations\n",
     sphere_smooth);
    nargs=1;
  }
  else if (!stricmp(option, "rusage"))
  {
    // just ignore
    nargs = 1 ;
  }
  else if (!stricmp(option, "int"))
  {
    noint = 0 ;
    nargs = 0 ;
  }
  else if (!stricmp(option, "verbose_low"))
  {
    parms.verbose=VERBOSE_MODE_LOW;
    fprintf(stderr,"verbose mode on (default+low mode)\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "warnings"))
  {
    parms.verbose=VERBOSE_MODE_MEDIUM;
    fprintf(stderr,"verbose mode on (medium mode): printing warnings\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "errors"))
  {
    parms.verbose=VERBOSE_MODE_HIGH;
    fprintf(stderr,"verbose mode on (high mode): "
            "exiting when warnings appear\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "movie"))
  {
    parms.movie=1;
    fprintf(stderr,"movie mode on\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "intersect"))
  {
    parms.check_surface_intersection=atoi(argv[2]);
    if (parms.check_surface_intersection)
    {
      fprintf(stderr,"check if the final surface self-intersects\n");
    }
    nargs = 1 ;
  }
  else if (!stricmp(option, "mappings"))
  {
    parms.optimal_mapping=atoi(argv[2]);
    if (parms.optimal_mapping)
    {
      fprintf(stderr,"generate several different mappings\n");
    }
    nargs = 1 ;
  }
  else if (!stricmp(option, "correct_defect"))
  {
    parms.correct_defect=atoi(argv[2]);
    fprintf(stderr,"correct defect %d only\n",parms.correct_defect);
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
  else if (!stricmp(option, "qcurv"))
  {
    parms.l_qcurv = atof(argv[2]) ;
    fprintf(stderr,"setting l_qcurv = %2.2f\n", parms.l_qcurv) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "unmri"))
  {
    parms.l_unmri = atof(argv[2]) ;
    fprintf(stderr,"setting l_unmri = %2.2f\n", parms.l_unmri) ;
    if (parms.volume_resolution<1)
    {
      parms.volume_resolution = 2;
    }
    nargs = 1 ;
  }
  else if (!stricmp(option, "vol"))
  {
    if (atoi(argv[2])<1)
    {
      parms.volume_resolution = 2;
    }
    else
    {
      parms.volume_resolution = atoi(argv[2]) ;
    }
    fprintf(stderr,"setting the volume resolution to %d\n",
            parms.volume_resolution) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "patches"))
  {
    parms.max_patches = atoi(argv[2]) ;
    fprintf(stderr,"using %d defect patches/generation...\n",
            parms.max_patches) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "generations"))
  {
    parms.max_unchanged = atoi(argv[2]) ;
    fprintf(stderr,"terminating evolution after %d generations "
            "without change...\n", parms.max_unchanged) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "niters"))
  {
    parms.niters = atoi(argv[2]) ;
    fprintf(stderr,"stopping genetic algorithm after %d iterations\n",
            parms.niters) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "optimize"))
  {
    parms.search_mode=GENETIC_SEARCH;
    parms.keep=0; /* to be discussed */
    parms.vertex_eliminate=1;
    parms.retessellation_mode=0;
    parms.initial_selection=1;
    parms.smooth=2;
    parms.match=1;
    parms.volume_resolution=2;
    parms.l_mri = 1.0f ;
    parms.l_curv = 1.0f ;
    parms.l_qcurv = 1.0f;
    parms.l_unmri = 10.0f ;
    nsmooth=0;
    add=0;
    fprintf(stderr,"using optimized parameters\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "ga"))
  {
    parms.search_mode=GENETIC_SEARCH;
    parms.keep=0;
    parms.vertex_eliminate=1;
    parms.retessellation_mode=0;
    parms.initial_selection=1;
    parms.smooth=2;
    parms.match=1;
    parms.volume_resolution=2;
    parms.l_mri = 1.0f ;
    parms.l_curv = 1.0f ;
    parms.l_qcurv = 1.0f;
    parms.l_unmri = 10.0f ;
    nsmooth=0;
    add=0;
    fprintf(stderr,"using genetic algorithm with optimized parameters\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "match"))
  {
    parms.match=atoi(argv[2]);
    if (parms.match)
    {
      fprintf(stderr,"match patch onto surface using local intensities\n") ;
    }
    else
      fprintf(stderr,"do not match patch onto surface using "
              "local intensities\n") ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "smooth"))
  {
    parms.smooth=atoi(argv[2]);
    if (parms.smooth)
    {
      fprintf(stderr,"smooth patch with mode %d\n",parms.smooth) ;
    }
    else
    {
      fprintf(stderr,"do not smooth patch\n") ;
    }
    nargs = 1 ;
  }
  else if (!stricmp(option, "select"))
  {
    parms.initial_selection=atoi(argv[2]);
    if (parms.initial_selection)
    {
      fprintf(stderr,"using qsphere to infer potential solutions\n") ;
    }
    else
    {
      fprintf(stderr,"random initial selection\n") ;
    }
    nargs = 1 ;
  }
  else if (!stricmp(option, "save"))
  {
    parms.save_fname=argv[2]; //name of the folder
    parms.defect_number=atoi(argv[3]);
    if (parms.defect_number>=0)
      fprintf(stderr,"save results in folder %s for the patch %d\n",
              parms.save_fname,parms.defect_number) ;
    else
      fprintf(stderr,"save results in folder %s for all patches\n",
              parms.save_fname) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "genetic"))
  {
    parms.search_mode=GENETIC_SEARCH;
    if (parms.volume_resolution<=0)
    {
      parms.volume_resolution=2;
    }
    fprintf(stderr,"using genetic algorithm\n") ;
    nargs = 0 ;
  }
  else if (!stricmp(option, "random"))
  {
    parms.search_mode=RANDOM_SEARCH;
    parms.niters=atoi(argv[2]);
    fprintf(stderr,"using random search with %d iterations to "
            "retessellate\n",
            parms.niters);
    nargs = 1 ;
  }
  else if (!stricmp(option, "variable"))
  {
    parms.retessellation_mode=atoi(argv[2]);
    if (parms.retessellation_mode)
    {
      fprintf(stderr,"ordering dependant final number of vertices\n");
    }
    else
    {
      fprintf(stderr,"keeping all vertices in ordered retessellation\n");
    }
    nargs = 1 ;
  }
  else if (!stricmp(option, "eliminate"))
  {
    parms.vertex_eliminate=atoi(argv[2]);
    if (parms.vertex_eliminate)
    {
      fprintf(stderr,"eliminate the less used vertices during search\n");
    }
    else
    {
      fprintf(stderr,"keep all vertices during search\n");
    }
    nargs = 1 ;
  }
  else if (!stricmp(option, "keep"))
  {
    parms.keep=atoi(argv[2]);
    if (parms.keep)
    {
      fprintf(stderr,"keep every vertex in the defect before search\n");
    }
    else
      fprintf(stderr,"select defective vertices in each defect "
              "before search\n");
    nargs = 1 ;
  }
  else if (!stricmp(option, "edge_table"))
  {
    parms.edge_table=atoi(argv[2]);
    if (parms.edge_table)
    {
      fprintf(stderr,"use precomputed edge table\n");
    }
    else
    {
      fprintf(stderr,"do not use precomputed edge table\n");
    }
    nargs = 1 ;
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
    topology_fixing_exit_after_diag = 1 ;
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
    fprintf(stderr,"writing fixed inflated coordinates to %s\n",
            inflated_name);
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
  else if (!stricmp(option, "openmp") || !stricmp(option, "threads") )
  {
    char str[STRLEN] ;
    int req = snprintf(str, STRLEN, "OMP_NUM_THREADS=%d", atoi(argv[2]));
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncatqion on line " << __LINE__ << std::endl;
    }
    putenv(str) ;
#ifdef HAVE_OPENMP
    omp_set_num_threads(atoi(argv[2]));
#else
    fprintf(stderr, "Warning - built without openmp support\n");
#endif
    nargs = 1 ;
    fprintf(stderr, "Setting %s\n", str) ;
  }
  else if (!stricmp(option, "orig"))
  {
    orig_name = argv[2] ;
    fprintf(stderr,"reading original coordinates from '%s'\n",orig_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "defect"))
  {
    defectbasename = argv[2] ;
    printf("defect basename set to '%s'\n",defectbasename);
    nargs = 1 ;
  }
  else if (!stricmp(option, "out"))
  {
    out_name = argv[2] ;
    printf("writing corrected surface to '%s'\n",out_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "brain"))
  {
    brain_name = argv[2] ;
    fprintf(stderr,"using '%s' as brain volume\n",brain_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "surf"))
  {
    surf_dir = argv[2] ;
    fprintf(stderr,"using '%s' instead of surf subdirectory\n",surf_dir);
    nargs = 1 ;
  }
  else if (!stricmp(option, "seed"))
  {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number genererator to %d\n",
            atoi(argv[2])) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "mgz"))
  {
    printf("INFO: assuming .mgz format\n");
    MGZ = 1;
  }
  else if (!stricmp(option, "nomgz"))
  {
    printf("INFO: assuming no .mgz format\n");
    MGZ = 0;
  }
  else switch (toupper(*option))
    {
    case 'S':
      nsmooth = atoi(argv[2]) ;
      nargs = 1 ;
      printf("smoothing corrected surface for %d iterations\n", nsmooth) ;
      break ;
    case '?':
    case 'H':
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

#include "mris_fix_topology.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_fix_topology_help_xml,
                mris_fix_topology_help_xml_len);
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
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static void print_parameters(void)
{
  fprintf(stderr,"\n***********************************"
          "**************************\n");
  fprintf(stderr,"Topology Correction Parameters\n");
  switch (parms.search_mode)
  {
  case GREEDY_SEARCH:
    fprintf(stderr,"retessellation mode :          greedy search\n");
    break;
  case RANDOM_SEARCH:
    fprintf(stderr,"retessellation mode :          random search\n");
    break;
  case GENETIC_SEARCH:
    fprintf(stderr,"retessellation mode:           genetic search\n");
    fprintf(stderr,"number of patches/generation : %d\n",parms.max_patches);
    fprintf(stderr,"number of generations :        %d\n",parms.max_unchanged);
    if (parms.niters>0)
    {
      fprintf(stderr,"stopping after %d iterations\n",parms.niters);
    }
    fprintf(stderr,"surface mri loglikelihood coefficient :         %1.1f\n",
            (float)parms.l_mri);
    fprintf(stderr,"volume mri loglikelihood coefficient :          %1.1f\n",
            (float)parms.l_unmri);
    fprintf(stderr,"normal dot loglikelihood coefficient :          %2.1f\n",
            (float)parms.l_curv);
    fprintf(stderr,"quadratic curvature loglikelihood coefficient : %2.1f\n",
            (float)parms.l_qcurv);
    fprintf(stderr,"volume resolution :                             %d\n",
            parms.volume_resolution);
    fprintf(stderr,"eliminate vertices during search :              %d\n",
            parms.vertex_eliminate);
    fprintf(stderr,"initial patch selection :                       %d\n",
            parms.initial_selection);
    break;
  }
  fprintf(stderr,"select all defect vertices :                    %d\n",
          parms.keep);
  fprintf(stderr,"ordering dependant retessellation:              %d\n",
          parms.retessellation_mode);
  fprintf(stderr,"use precomputed edge table :                    %d\n",
          parms.edge_table);
  fprintf(stderr,"smooth retessellated patch :                    %d\n",
          parms.smooth);
  fprintf(stderr,"match retessellated patch :                     %d\n",
          parms.match);
  fprintf(stderr,"verbose mode :                                  %d\n",
          parms.verbose);
  fprintf(stderr,"\n****************************************"
          "*********************\n");
}




