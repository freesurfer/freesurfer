/**
 * @brief create a parcellation of equal areas
 *
 * make a parcellation where each unit is assigned based on the face 
 * it maps to in a specified spherical surface (usually an icosahedral one)
 * 
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

#include "macros.h"
#include "timer.h"
#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "mrishash.h"
#include "fio.h"
#include "colortab.h"


#define MAX_PARCEL_VERTICES 10000

typedef struct
{
  double *cx, *cy, *cz, *energy, energy_num, energy_den ;
  double *var_within, *var_between ; // for energy_type == ENERGY_VARIANCE
  MRI    *mri_min_timecourse ;
  int    nparcels ;
  int    *num;
  int    *min_vno;
  int    **nbrs ;
} PARCEL_STATS ;

typedef struct
{
  int    write_iterations ;
  char   base_name[STRLEN] ;
  int    max_iterations ;
  double tol ;
  double l_markov ;
  double l_gaussian ;
  double l_eden ;
  double l_border ;
  double l_area ;
  double gaussian_std ;  
  double l_var ;
  PARCEL_STATS stats ;
  MRI    *mri_dmat ;
  MRI    *mri_computed ;
  char   dmat_read_fname[STRLEN] ;
  char   dmat_write_fname[STRLEN] ;
  int    read_dmat ;
  int    write_dmat ;
  int    energy_type ;
  MRI_SURFACE *mris_ico ;
  MRI    *mri_means, *mri_vars, *mri_stats ;
  double avg_area ;
} PARMS ;

#define ENERGY_SIMILARITY    0
#define ENERGY_ETA           1
#define ENERGY_DISTANCE      2
#define ENERGY_VARIANCE      3

int main(int argc, char *argv[]) ;

static int write_annot_correlations(MRI_SURFACE *mris, MRI *mri_cmatrix, PARMS *parms, char *fname) ;
static double compute_distance_energy(MRI_SURFACE *mris, PARMS *parms, MRI *mri_cmatrix, 
                                      double *penergy_num, double *penergy_den)  ;
static double compute_variance_energy(MRI_SURFACE *mris, PARMS *parms, MRI *mri_cmatrix, MRI *mri_means, 
                                      MRI *mri_vars, double *penergy_num, double *penergy_den) ;
static double compute_variance_energy_for_vertex_change(MRI_SURFACE *mris, PARMS *parms, 
                                                        MRI *mri_cmatrix, MRI *mri_means, 
                                                        MRI *mri_vars, double *penergy_num, 
                                                        double *penergy_den, int vno,
                                                        int old_parcel, int new_parcel) ;
static int mark_border_vertices(MRI_SURFACE *mris)  ;
static int segment_similarity_matrix(MRI_SURFACE *mris, MRI_SURFACE *mris_ico,  
                                     MRI *mri_smatrix, PARMS *parms) ;
static MRI *compute_eta_matrix(MRI_SURFACE *mris, MRI *mri_cmatrix, MRI *mri_ematrix, 
                               double thresh, double max_dist) ;
static MRI *compute_similarity_matrix(MRI_SURFACE *mris, MRI *mri_cmatrix, MRI *mri_smatrix, 
                                      double thresh, double max_dist) ;
static int allocate_stats(MRI_SURFACE *mris, PARMS *parms, int nparcels, MRI *mri_cmatrix) ;
static double L1_distance(MRI *mri1, MRI *mri2, int x1, int y1, int z1, int x2, int y2, int z2) ;
static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static int find_parcel_min_timecourses(MRI_SURFACE *mris, MRI *mri_cmatrix, PARMS *parms) ;
static void print_version(void) ;
static int write_snapshot(MRI_SURFACE *mris, PARMS *parms, int n, MRI *mri_cmatrix) ;
static double markov_energy(MRI_SURFACE *mris) ;
static double border_energy(MRI_SURFACE *mris) ;
static double area_energy(MRI_SURFACE *mris, PARMS *parms) ;
static double gaussian_energy(MRI_SURFACE *mris, PARMS *parms) ;
static int adjust_parcellation_boundaries(MRI_SURFACE *mris, MRI_SURFACE *mris_ico, MRI *mri_cmatrix,
                                          PARMS *parms) ;


static double compute_parcellation_energy_change(MRI_SURFACE *mris, 
                                                 PARMS *parms,MRI *mri_cmatrix,
                                                 MRI *mri_means, MRI *mri_vars,
                                                 int vno,int parcel_to_move_to,
                                                 double *den_old,
                                                 double *den_new,
                                                 double *num_old,
                                                 double *num_new) ;
static double compute_parcellation_energy(MRI_SURFACE *mris, MRI *mri_stats, PARMS *parms, double *penergy, MRI *mri_means, MRI *mri_vars, MRI *mri_cmatrix) ;
static double
compute_parcellation_median_energy(MRI_SURFACE *mris, MRI *mri_stats, PARMS *parms, MRI *mri_means, MRI *mri_cmatrix);
static int build_parcellation_border_permutation(MRI_SURFACE *mris, int *vertex_permutation, 
                                                 int *pnborder) ;
static int compute_parcellation_statistics(MRI_SURFACE *mris, MRI *mri_cmatrix, MRI *mri_stats, MRI *mri_means, MRI *mri_vars, PARMS *parms);


static int
update_parcellation_statistics(MRI_SURFACE *mris, int vno, int old_parcel, int new_parcel,
                               MRI *mri_cmatrix, MRI *mri_means, MRI *mri_vars, 
                               MRI *mri_stats, PARMS *parms);
const char *Progname ;

static char *cmatrix_fname = NULL ;
static MRI  *mri_cmatrix ;
static PARMS  parms ;
static char *evaluate_fname = NULL ;
static char *log_fname = NULL ;
static int randomize = 0 ;
static int  do_vertices =1 ;
static char *write_corr_fname = NULL ;
static char *write_annot_fname = NULL ;
char *ctab_fname = NULL ;
int InflatedOK = 0;

int
main(int argc, char *argv[]) {
  char               **av, *in_fname, *ico_fname, *out_fname, path[STRLEN], ico_name[STRLEN] ;
  int                ac, nargs ;
  float              scale ;
  MRI_SURFACE        *mris, *mris_ico ;
  float              radius ;
  int                fno, vno, annot,k ;
  double             fdist ;
  FACE               *face ;
  MHT                *mht ;
  VERTEX             *v ;
  int                msec, minutes, seconds ;
  MRI                *mri_smatrix = NULL ;


  parms.max_iterations = 100 ;
  parms.tol = 1e-6 ;
  parms.l_markov = 0 ;  // 0.001
  parms.l_gaussian = 0 ; // 0.01
  parms.l_var = 1.0 ;
  parms.energy_type = ENERGY_VARIANCE ;
  parms.l_eden = 1.0 ;
  
  Timer start;

  std::string cmdline = getAllInfo(argc, argv, "mris_make_face_parcellation");

  setRandomSeed(1L) ;

  nargs = handleVersionOption(argc, argv, "mris_make_face_parcellation");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if(argc < 4) print_help() ;

  in_fname = argv[1] ;
  ico_fname = argv[2] ;
  out_fname = argv[3] ;
  FileNamePath(out_fname, path) ;
  FileNameOnly(ico_fname, ico_name) ;

  if(strcmp("lh.inflated",fio_basename(in_fname,NULL))==0 ||
     strcmp("rh.inflated",fio_basename(in_fname,NULL))==0){
    if(InflatedOK == 0){
      printf("ERROR: do not use the inflated surface for this function\n");
      printf(" You probably want to use sphere or sphere.reg\n");
      exit(1);
    }
    printf("You have selected inflated surface\n");    
  }

  printf("Reading %s\n",in_fname) ;
  mris = MRISread(in_fname) ;
  if (!mris)  ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",Progname, in_fname) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISresetNeighborhoodSize(mris, 3) ; // reset current size to 1-nbrs

  printf("Reading %s\n",ico_fname) ;
  mris_ico = MRISread(ico_fname) ;
  if (!mris_ico) ErrorExit(ERROR_NOFILE, "%s: could not read ico file %s",Progname, ico_fname) ;
  if (mris_ico->nfaces < 256) scale = 256 / mris_ico->nfaces ;
  else                        scale = 1 ;
  parms.mris_ico = mris_ico ;
  radius = MRISaverageRadius(mris) ;
  MRISscaleBrain(mris_ico, mris_ico, radius / mris_ico->radius) ;
  MRIScomputeMetricProperties(mris_ico) ;
  MRISsaveVertexPositions(mris_ico, CANONICAL_VERTICES) ;

  if (cmatrix_fname)  // read in a correlation matrix for refining the parcellation
  {
    double max_dist ;
    mri_cmatrix = MRIread(cmatrix_fname) ;
    if (mri_cmatrix == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read correlation file %s\n", Progname, cmatrix_fname) ;
    
    if (mri_cmatrix->width != mris->nvertices || mri_cmatrix->nframes != mris->nvertices)
      ErrorExit(ERROR_BADFILE, "%s: cmatrix must be %d x 1 x 1 x %d",Progname,mris->nvertices,mris->nvertices);
    max_dist = 10*radius/(mris_ico->nvertices) ;
    if (parms.energy_type == ENERGY_SIMILARITY)
      mri_smatrix = compute_similarity_matrix(mris, mri_cmatrix, NULL, 0, max_dist) ;
    else if (parms.energy_type == ENERGY_ETA)
      mri_smatrix = compute_eta_matrix(mris, mri_cmatrix, NULL, 0.3, max_dist) ;
  }

  if (write_annot_fname)
  {
    if (MRISreadAnnotation(mris, write_annot_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not load annotation from %s",Progname,write_annot_fname);
    MRIScopyAnnotationsToMarkedIndex(mris) ;
    parms.stats.nparcels = MRISmaxMarked(mris)+1 ;
    write_annot_correlations(mris, mri_cmatrix, &parms, write_corr_fname) ;
    exit(0) ;
  }
  if (evaluate_fname)  // load  a previously computed annotation and compute it's energy
  {
    MRI *mri_stats, *mri_means, *mri_vars ;
    double energy, penergy ;
    int    nparcels ;
    FILE   *fp ;

    if (MRISreadAnnotation(mris, evaluate_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not load annotation %s",Progname,evaluate_fname) ;
    MRIScopyAnnotationsToMarkedIndex(mris) ;
    nparcels = MRISmaxMarked(mris)+1 ;

    allocate_stats(mris, &parms, nparcels, mri_cmatrix) ;
    printf("annotation with %d units read\n", nparcels) ;
    mri_stats = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, 3) ;
    mri_means = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, mri_cmatrix->nframes) ;
    mri_vars = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, mri_cmatrix->nframes) ;
    parms.mri_stats = mri_stats ; parms.mri_means = mri_means ; parms.mri_vars = mri_vars ;
    compute_parcellation_statistics(mris, mri_cmatrix, mri_stats, mri_means,mri_vars,&parms);
    penergy = energy = compute_parcellation_median_energy(mris, mri_stats, &parms, mri_means, mri_cmatrix);
    printf("energy = %2.3f\n", penergy) ;
    if (randomize)
    {
      MRI *mri_rand ;
      double rand_energy, rand_penergy, total_rand_energy, total_rand_penergy ;
      int    n, vnos[300000], tmp, index ;

      setRandomSeed(0L) ;

#define MAX_RAND 10
      MRIScopyMarkedToMarked2(mris) ;
      for (total_rand_energy = total_rand_penergy = 0.0, n = 0 ; n < MAX_RAND ; n++)
      {
        for (vno= 0 ; vno < mris->nvertices ; vno++)
          vnos[vno] = vno ;
        for (vno= 0 ; vno < mris->nvertices ; vno++)
        {
          index = (int)randomNumber(0.0, mris->nvertices-.00001) ;
          tmp = vnos[index] ;
          vnos[index] = vnos[vno] ;
          vnos[vno] = tmp;
        }
        for (vno= 0 ; vno < mris->nvertices ; vno++)
          mris->vertices[vno].marked = mris->vertices[vnos[vno]].marked2 ;
#if 0
        mri_rand = MRIdrand48(mri_cmatrix->width, mri_cmatrix->height, mri_cmatrix->depth, 
                              mri_cmatrix->nframes, 0, 1, NULL) ;
#else
        mri_rand = NULL ;
#endif
        
        if (parms.energy_type == ENERGY_VARIANCE)
          compute_parcellation_statistics(mris, mri_rand, mri_stats,mri_means,mri_vars,&parms);
#if 1
        rand_energy = rand_penergy = compute_parcellation_median_energy(mris, mri_stats, &parms, mri_means, mri_cmatrix);
#else
        rand_energy = compute_parcellation_energy(mris,mri_stats,&parms,&rand_penergy,mri_means,mri_vars, mri_cmatrix);
#endif
        printf("rand_energy = %2.5f (%2.5f)\n", rand_energy, rand_penergy) ;

        total_rand_energy += rand_energy ; total_rand_penergy += rand_penergy ;
        //        MRIfree(&mri_rand) ;
      }
      total_rand_energy /= n ; total_rand_penergy /= n ;
      energy /= total_rand_energy ; penergy /= total_rand_penergy ;
    }
    printf("energy = %2.5f (%2.5f)\n", energy, penergy) ;
    fp = fopen(log_fname, "a") ;
    fprintf(fp, "%d %f %f\n", nparcels, energy, penergy) ;
    fclose(fp) ;
    exit(0) ;
  }

  MRISaddCommandLine(mris, cmdline) ;

  mris->ct = CTABalloc(do_vertices ? mris_ico->nvertices : mris_ico->nfaces) ;
  // Search an additional 100 times for a unique set of RGBs
  CTABunique(mris->ct, 100);
  if(CTABcountRepeats(mris->ct, 1) != 0){
    printf("ERROR: could not find a unique color table\n");
    exit(1);
  }
  printf("parcellating hemisphere into %d units\n", mris->ct->nentries);
  strcpy (mris->ct->fname, ico_fname);
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++) {
    int req = snprintf (mris->ct->entries[vno]->name, STRLEN, "%s_vertex_%d", ico_name, vno);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }
  // Note: there was logic here that attempted to assure that nearby patches did not 
  // have similar color, but it caused some patches to have the same color, and so
  // resolve to the same ROI. Now just uses random ctab
  if(ctab_fname) CTABwriteFileASCII(mris->ct,ctab_fname);

  printf("do_vertices = %d\n",do_vertices);
  if(do_vertices){
    int *nhits;
    nhits = (int*)calloc(mris_ico->nvertices,sizeof(int));
    radius = MRISaverageRadius(mris) ;
    MRISscaleBrain(mris_ico, mris_ico, radius / mris_ico->radius) ;
    MRIScomputeMetricProperties(mris_ico) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      if ((((vno % (mris->nvertices/10))) == 0) && DIAG_VERBOSE_ON)
        printf("%2.1f%% done\n", 100.0*(float)vno / mris->nvertices) ;
      v = &mris->vertices[vno] ;
      fno = MRISfindClosestVertex(mris_ico, v->cx, v->cy, v->cz, NULL, CURRENT_VERTICES) ;
      CTABannotationAtIndex(mris->ct, fno, &annot);
      v->annotation = annot ;
      v->marked = fno ;
      CTABfindAnnotation(mris->ct, annot, &k);
      nhits[k] ++;
    }
    // Check that all have representation
    for(k=0; k<mris_ico->nvertices;k++){
      if(nhits[k] == 0) {
	v = &(mris_ico->vertices[k]);
	printf("Parcellation %d is empty\n",k);
	fflush(stdout);
      }
    }
    free(nhits);
  }
  else{
    mht = MHTcreateFaceTable_Resolution(mris_ico, CURRENT_VERTICES, 1.0);
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      if ((((vno % (mris->nvertices/10))) == 0) && DIAG_VERBOSE_ON)
        printf("%2.1f%% done\n", 100.0*(float)vno / mris->nvertices) ;
      v = &mris->vertices[vno] ;
      MHTfindClosestFaceGeneric(mht, mris, 
                                v->x, v->y, v->z, 
                                mris_ico->avg_vertex_dist/10, 
                                mris_ico->avg_vertex_dist/10, 
                                0,
                                &face, &fno, &fdist) ;
      if (fno < 0)
      {
        fno = MHTBruteForceClosestFace(mris_ico, v->x, v->y, v->z, 
                                       CURRENT_VERTICES, NULL);    
        if (fno  < 0)
        {
          printf("warning: v %d not found in MHT\n", vno) ;
          continue ;
        }
      }
      CTABannotationAtIndex(mris->ct, fno, &annot);
      v->annotation = annot ;
      v->marked = fno ;
    }
  }


  if (mri_smatrix)
  {
    char fname[STRLEN], *cp ;
    FileNameOnly(out_fname, fname) ;
    strcpy(parms.base_name, fname) ;
    cp = strrchr(parms.base_name, '.') ;
    if (cp)
      *cp = 0 ; // take out extension

    allocate_stats(mris, &parms, do_vertices ? mris_ico->nvertices : mris_ico->nfaces, mri_cmatrix) ;
    segment_similarity_matrix(mris, mris_ico, mri_smatrix, &parms) ;
  }
  else if (mri_cmatrix)
  {
    char fname[STRLEN], *cp ;
    if (parms.l_gaussian > 0)
    {
      parms.gaussian_std = 0.1 * sqrt(mris_ico->total_area / (double)mris_ico->nfaces) ;
      printf("setting gaussian std to %f\n", parms.gaussian_std) ;
    }
    FileNameOnly(out_fname, fname) ;
    strcpy(parms.base_name, fname) ;
    cp = strrchr(parms.base_name, '.') ;
    if (cp)
      *cp = 0 ; // take out extension

    allocate_stats(mris, &parms, do_vertices ? mris_ico->nvertices : mris_ico->nfaces, mri_cmatrix) ;
    adjust_parcellation_boundaries(mris, mris_ico, mri_cmatrix, &parms) ;
    if (write_corr_fname)
      write_annot_correlations(mris, mri_cmatrix, &parms, write_corr_fname) ;
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing annotation to %s\n", out_fname) ;
  MRISwriteAnnotation(mris, out_fname) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "parcellation took %d minutes and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))      print_help() ;
  else if (!stricmp(option, "help"))  print_help() ;
  else if (!stricmp(option, "markov")){
    parms.l_markov = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting markov coefficient to %2.6f\n", parms.l_markov) ;
  }
  else if (!stricmp(option, "write_corr")){
    write_corr_fname  = argv[2] ;
    nargs = 1 ;
    printf("writing parcellation correlation matrix to %s\n", write_corr_fname) ;
  }
  else if (!stricmp(option, "wonly")){
    write_corr_fname  = argv[2] ;
    write_annot_fname = argv[3] ;
    nargs = 2 ;
    printf("writing parcellation correlation matrix from %s to %s\n", 
           write_annot_fname, write_corr_fname) ;
  }
  else if (!stricmp(option, "var")){
    parms.l_var = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting variance coefficient to %2.3f\n", parms.l_var) ;
  }
  else if (!stricmp(option, "eden")){
    parms.l_eden = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting eden coefficient to %2.3f\n", parms.l_eden) ;
  }
  else if (!stricmp(option, "read_dmat")){
    parms.read_dmat = 1 ;
    strcpy(parms.dmat_read_fname, argv[2]) ;
    nargs = 1 ;
    printf("reading dmat from %s\n", parms.dmat_read_fname) ;
  }
  else if (!stricmp(option, "write_dmat")){
    parms.write_dmat = 1 ;
    strcpy(parms.dmat_write_fname, argv[2]) ;
    nargs = 1 ;
    printf("writing dmat to %s\n", parms.dmat_write_fname) ;
  }
  else if (!stricmp(option, "tol")){
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting tol to %2.7f\n", parms.tol) ;
  }
  else if (!stricmp(option, "ctab")){
    ctab_fname = argv[2];
    nargs = 1 ;
    printf("Writing ctab to %s\n",ctab_fname);
  }
  else if (!stricmp(option, "-version")){
    print_version() ;
  } else switch (toupper(*option)) {
  case 'C':
    cmatrix_fname = argv[2] ;
    nargs = 1 ;
    printf("reading correlation matrix from %s\n", cmatrix_fname) ;
    break ;
  case 'G':
    parms.l_gaussian = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting gaussian coefficient to %2.6f\n", parms.l_gaussian) ;
    break ;
  case 'F':
    do_vertices = 0 ;
    printf("writing out face parcellation instead of vertex one\n") ;
    break ;
  case 'B':
    parms.l_border = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting border coefficient to %2.6f\n", parms.l_border) ;
    break ;
  case 'A':
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting area coefficient to %2.6f\n", parms.l_area) ;
    break ;
  case 'M':
    parms.max_iterations = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting max iterations to %d\n", parms.max_iterations) ;
    break ;
  case 'I':
    nargs = 0 ;
    InflatedOK = 1;
    printf("OK to use inflated\n");
    break ;
  case 'R':
    randomize = atof(argv[2]) ;
    printf("%susing randomization to normalize energies\n", randomize ? "" : "not ") ;
    nargs = 1 ;
    break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      printf("debugging vertex %d\n", Gdiag_no) ;
      nargs = 1 ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("setting write iterations to %d\n", parms.write_iterations) ;
      break ;
  case 'E':
    evaluate_fname = argv[2] ;
    log_fname = argv[3] ;
    printf("evalauting energy of annotation file %s and logging results to %s\n", evaluate_fname,
           log_fname);
    nargs = 2 ;
    break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
print_usage(void) {
  printf("usage: %s [options] <input surface> <ico file> <output annot>\n\n"
          "example: %s lh.sphere.reg $FREESURFER_HOME/lib/bem/ic3.tri ./lh.ic3.annot\n",
          Progname, Progname) ;
  printf("  surf should be either:\n");
  printf("    sphere: units will be approximately equal size within subject but \n"
	 "            not in correspondence across subjects\n");
  printf("    sphere.reg: units will be different sizes within subject but \n"
	 "            in correspondence across subjects\n");
  printf("  Note: do not use inflated as was suggested by previous versions! \n");
  printf("  \n");
}

static void
print_help(void) {
  print_usage() ;
  printf(
    "\nThis generates a parcellation based on which icosahedral face each vertex maps to.\n"
    "Options: \n"
    "  -ctab colortable.txt\n"
    "  \n"
    );
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}
static int
adjust_parcellation_boundaries(MRI_SURFACE *mris, MRI_SURFACE *mris_ico, MRI *mri_cmatrix,
                               PARMS *parms)
{
  MRI     *mri_stats, *mri_means, *mri_vars ;
  int     *vertex_permutation, parcel, nborder, done, nchanged, iter, index, *nbrs, nparcels,  
    min_parcel,vno, n, nframes;
  double  energy, last_energy, parc_energy, menergy, energy_change, 
    best_energy_change, last_penergy ;

  energy = last_energy = 0.0 ;
  nframes = mri_cmatrix->nframes ;
  nparcels = parms->stats.nparcels ;
  mri_stats = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, 3) ;
  mri_means = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, nframes) ;
  mri_vars = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, nframes) ;
  parms->mri_stats = mri_stats ; parms->mri_means = mri_means ; parms->mri_vars = mri_vars ;

  nbrs = (int *)calloc(nparcels, sizeof(int)) ;
  vertex_permutation = (int *)calloc(mris->nvertices, sizeof(int)) ;
  if (vertex_permutation == NULL)
    ErrorExit(ERROR_NOMEMORY, "adjust_parcellation_boundaries: could not allocated permutation");
  iter = done = 0 ; 
  do
  {
    nchanged = 0 ;
    build_parcellation_border_permutation(mris, vertex_permutation, &nborder) ;

    if (iter == 0)
    {
      compute_parcellation_statistics(mris, mri_cmatrix, mri_stats,mri_means,mri_vars,parms) ;
      last_energy = energy = 
        compute_parcellation_energy(mris, mri_stats, parms, &parc_energy, 
                                    mri_means, mri_vars, mri_cmatrix) ;

      menergy = compute_parcellation_median_energy(mris, mri_stats, parms, 
                                                   mri_means, mri_cmatrix);
      printf("iter %d: nchanged = %d, nborder = %d, energy = %2.5f (%2.5f, %2.5f)\n", 
             iter, nchanged,nborder,energy,parc_energy,menergy);
    }
    else if (iter == Gdiag_no)
      DiagBreak() ;
    if (parms->write_iterations > 0 && ((iter % parms->write_iterations) == 0))
      write_snapshot(mris, parms, iter, mri_cmatrix) ;
    for (index = 0 ; index < nborder ; index++)
    {
      vno = vertex_permutation[index] ;
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
       if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      parcel = v->marked ;
      if (parcel < 0)
        continue ;

      // build list of all neihboring parcels at this point
      memset(nbrs, 0, nparcels*sizeof(int)) ;
      min_parcel = parcel ; best_energy_change = 0 ;
      for (n = 0 ; n < vt->vnum ; n++)
      {
        VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
        if (vt->v[n] == Gdiag_no)
          DiagBreak() ;
        if (vn->marked != parcel && nbrs[vn->marked] != 1) // not already done
        {
          // try changing the parcel and see if it decreases the energy
          nbrs[vn->marked] = 1 ;  // only do it once
          energy_change = compute_parcellation_energy_change(mris, parms, mri_cmatrix, mri_means, mri_vars,
                                                             vno, vn->marked, NULL, NULL, NULL, NULL) ;
          if (energy_change < best_energy_change)
          {
            best_energy_change = energy_change ;
            min_parcel = vn->marked ;
          }
        }
      }
      if (min_parcel != parcel)
      {
        int annot ;
        double den_old, den_new, num_old, num_new ;

        energy_change = 
          compute_parcellation_energy_change(mris, parms, mri_cmatrix, mri_means, mri_vars, vno, 
                                             min_parcel, &den_old, &den_new, 
                                             &num_old, &num_new) ;
        parms->stats.energy_num = parms->stats.energy_num + num_new - num_old ;
        parms->stats.energy_den = parms->stats.energy_den + den_new - den_old ;
        if (parms->energy_type == ENERGY_VARIANCE)
          update_parcellation_statistics(mris, vno, v->marked, min_parcel,
                                         mri_cmatrix, mri_means, mri_vars, mri_stats, parms);
        v->marked = min_parcel ;
        nchanged++ ;
        CTABannotationAtIndex(mris->ct, v->marked, &annot);
        v->annotation = annot ;
        if (vno == Gdiag_no)
        {
          double e, pe ;
          compute_parcellation_statistics(mris, mri_cmatrix, mri_stats,
                                          mri_means,mri_vars,parms) ;
          e = compute_parcellation_energy(mris, mri_stats, parms, &pe, mri_means, mri_vars, mri_cmatrix) ;
          if (e > last_energy)
            DiagBreak() ;
          DiagBreak() ;
        }
      }
      //      energy = compute_parcellation_energy(mris, mri_stats, parms, &parc_energy, mri_means, mri_vars, mri_cmatrix) ;
    }
    last_energy = energy ;
    last_penergy = parc_energy ;
    compute_parcellation_statistics(mris, mri_cmatrix, mri_stats,
                                    mri_means,mri_vars,parms) ;
    energy = compute_parcellation_energy(mris, mri_stats, parms, &parc_energy, mri_means, mri_vars, mri_cmatrix) ;
    done = (nchanged == 0) || (iter++ > parms->max_iterations) ||
      ((last_penergy-parc_energy)/last_penergy < parms->tol) ;
        
    menergy = compute_parcellation_median_energy(mris, mri_stats, parms, mri_means, mri_cmatrix);
    printf("iter %d: nchanged = %d, nborder = %d, energy = %2.5f (%2.5f, %2.5f)\n", 
           iter, nchanged,nborder,energy, parc_energy,menergy);
  }  while (!done) ;
  if (parms->write_iterations > 0 && ((iter % parms->write_iterations) == 0))
    write_snapshot(mris, parms, iter, mri_cmatrix) ;
  MRIfree(&mri_stats) ; MRIfree(&mri_means) ; MRIfree(&mri_vars) ;
  free(nbrs) ; free(vertex_permutation) ;
  return(NO_ERROR) ;
}

static double
compute_parcellation_energy(MRI_SURFACE *mris, MRI *mri_stats, PARMS *parms, 
                            double *penergy, MRI *mri_means, MRI *mri_vars, 
                            MRI *mri_cmatrix)
{
  double  energy = 0.0, parc_energy, markov, gauss ;

  if (parms->l_var > 0)
  {
    double energy_num, energy_den ;

    switch (parms->energy_type)
    {
    case ENERGY_DISTANCE:
      energy = compute_distance_energy(mris, parms, mri_cmatrix, &energy_num, &energy_den) ;
      parms->stats.energy_num = energy_num ;
      parms->stats.energy_den = energy_den ;
      break ;
    case ENERGY_VARIANCE:
      energy = compute_variance_energy(mris, parms, mri_cmatrix, mri_means, mri_vars, &energy_num, &energy_den) ;
      break ;
    default:
      ErrorExit(ERROR_BADPARM, "compute_parcellation_energy: unknown energy type %d", Progname, parms->energy_type);
    }
  }
  else
    energy = 0.0 ;

  if (energy < 0)
    DiagBreak() ;
  energy = parc_energy = parms->l_var * energy ;
  if (parms->l_markov > 0)
  {
    markov = markov_energy(mris) ;
    energy +=  parms->l_markov * markov;
  }
  if (parms->l_border > 0)
  {
    double border ;
    border = border_energy(mris) ;
    energy +=  parms->l_border * border;
  }
  if (parms->l_area > 0)
  {
    double area  ;
    area = area_energy(mris, parms) ;
    energy +=  parms->l_area * area;
  }
  if (parms->l_gaussian > 0)
  {
    gauss = gaussian_energy(mris, parms) ;
    energy += parms->l_gaussian * gauss;
  }
  *penergy = parc_energy ;
  return(energy) ;
}


static double
compute_distance_energy(MRI_SURFACE *mris, PARMS *parms, MRI *mri_cmatrix, double *penergy_num, double *penergy_den) 
{
  double  energy, dist_within, dist_between, dist, dist_between_total ;
  int     parcel, n, nparcels, vno, n_nbrs ;
  VERTEX  *v ;

  nparcels = parms->stats.nparcels ;

  memset(parms->stats.energy, 0, sizeof(double)*nparcels) ;
  for (dist_within = dist_between_total = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    parcel = v->marked ;
    //      dist = L1_distance(mri_cmatrix, parms->stats.mri_min_timecourse, vno, 0, 0, parcel, 0, 0) ;
    dist = MRIgetVoxVal(parms->mri_dmat,vno, parms->stats.min_vno[parcel], 0, 0) ;
    dist_within += dist ;
    for (dist_between = 0.0, n_nbrs = n = 0 ; n < nparcels ; n++)
      if (n != parcel && parms->stats.nbrs[parcel][n] > 0)
      {
        //          dist = L1_distance(mri_cmatrix, parms->stats.mri_min_timecourse, vno, 0, 0, n, 0, 0) ;
        dist = MRIgetVoxVal(parms->mri_dmat, vno, parms->stats.min_vno[n], 0, 0) ;
        dist_between += dist ;
        n_nbrs++ ;
      }
    dist_between_total += dist_between/n_nbrs ;
  }
  dist_between /= n_nbrs ;
  energy = dist_within / (parms->l_eden * dist_between_total) ;
  *penergy_num = dist_within ;
  *penergy_den = dist_between_total ;
  return(energy) ;
}

static int
build_parcellation_border_permutation(MRI_SURFACE *mris, int *vertex_permutation, int *pnborder)
{
  int    vno, nborder, n, tmp, index ;

  for (nborder = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < vt->vnum ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vt->v[n] == Gdiag_no)
        DiagBreak() ;
      if (vn->marked != v->marked)
        break ;
    }
    if (n < vt->vnum)   // is a border vertex
      vertex_permutation[nborder++] = vno ;
  }

  for (n = 0 ; n < nborder ; n++)
  {
    index = randomNumber(0.0, (double)(nborder-0.0001)) ;
    tmp = vertex_permutation[n] ;
    vertex_permutation[n] = vertex_permutation[index] ;
    vertex_permutation[index] = tmp ;
  }

  *pnborder = nborder ;
  return(NO_ERROR) ;
}

/*
  frame 0 of mri_stats will be the mean for the ith parcellation unit
  frame 1 of mri_stats will be the variance for the ith parcellation unit
  frame 2 of mri_stats will be the # of vertices in the ith parcellation unit
*/
static int
compute_parcellation_statistics(MRI_SURFACE *mris, MRI *mri_cmatrix, MRI *mri_stats,
                                MRI *mri_means, MRI *mri_vars, PARMS *parms)
{
  
  int    vno, parcel, nparcels, frame, nframes ;
  VERTEX *v ;
  double  var, mean,val, dof, mean_total, var_total, norm ;

  if (parms->energy_type == ENERGY_DISTANCE)
    find_parcel_min_timecourses(mris, mri_cmatrix, parms) ;
  if (Gdiag_no >= 0) 
    printf("computing parcellation statistics\n") ;
  nframes = mri_cmatrix->nframes ;
  nparcels = mri_stats->width ;
  MRIclear(mri_stats) ; MRIclear(mri_means) ; MRIclear(mri_vars) ;

  memset(parms->stats.num, 0, nparcels*sizeof(int)) ;
  memset(parms->stats.cx, 0, nparcels*sizeof(double)) ;
  memset(parms->stats.cy, 0, nparcels*sizeof(double)) ;
  memset(parms->stats.cz, 0, nparcels*sizeof(double)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    parcel = v->marked ;
    parms->stats.cx[parcel] += v->cx ; 
    parms->stats.cy[parcel] += v->cy ; 
    parms->stats.cz[parcel] += v->cz ; 
    parms->stats.num[parcel]++ ;
    norm = sqrt(SQR(parms->stats.cx[parcel]) + SQR(parms->stats.cy[parcel]) + 
                SQR(parms->stats.cz[parcel]));
#if 0
    parms->stats.cx[parcel] /= norm ; 
    parms->stats.cy[parcel] /= norm ; 
    parms->stats.cz[parcel] /= norm ;
#endif
  }

  parms->avg_area = (double)mris->nvertices / parms->stats.nparcels ;
  for (parcel = 0 ; parcel < nparcels ; parcel++)
    if (parms->stats.num[parcel] > 0)
    {
      parms->stats.cx[parcel] /= parms->stats.num[parcel] ; 
      parms->stats.cy[parcel] /= parms->stats.num[parcel] ; 
      parms->stats.cz[parcel] /= parms->stats.num[parcel] ;
    }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->marked < 0)
      continue ;
    parcel = v->marked ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dof = MRIgetVoxVal(mri_stats, v->marked, 0, 0, 2)+1 ;
    MRIsetVoxVal(mri_stats, v->marked, 0, 0, 2, dof) ;
    if (v->marked == Gdiag_no && Gdiag_no >= 0)
      printf("adding %d to parcel %d, dofs = %d, ", vno, v->marked, (int)dof) ;

    for (frame = 0 ; frame < nframes ; frame++)
    {
      val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
      mean = MRIgetVoxVal(mri_means, v->marked, 0, 0, frame) ;
      var = MRIgetVoxVal(mri_vars, v->marked, 0, 0, frame) ;
      mean += val ;
      var += val*val ;
      if (parcel == 133 && frame == 322 && DIAG_VERBOSE_ON)
      {
        printf("adding vno %d: %f to parcel %d, now %f +- %f (var=%f), dof=%d\n",
               vno, val, parcel, mean, var,
               (var - (mean*mean/dof)) / (dof-1),(int)dof);
      }
      if (dof > 2 && (((var - (mean*mean/dof)) / (dof-1)) < 0))
      {
        double m, v ;
        m = mean / dof ;
        v = (var - dof*m*m) / (dof-1);

        if (v < -1e-6)
          DiagBreak() ;
      }
            
      MRIsetVoxVal(mri_means, v->marked, 0, 0, frame, mean) ;
      MRIsetVoxVal(mri_vars, v->marked, 0, 0, frame, var) ;
      if (v->marked == Gdiag_no && frame == 2 && Gdiag_no >= 0)
        printf("\tval = %2.3f, mean = %2.3f, var = %2.5f\n", val, mean, var) ;
    }
  }

  for (parcel = 0 ; parcel < nparcels ; parcel++)
  {
    dof = MRIgetVoxVal(mri_stats, parcel, 0, 0, 2) ;
    if (dof <= 0)
      continue ;
    for (mean_total = var_total = 0.0, frame = 0 ; frame < nframes ; frame++)
    {
      mean = MRIgetVoxVal(mri_means, parcel, 0, 0, frame) ;
      var = MRIgetVoxVal(mri_vars, parcel, 0, 0, frame) ;
      mean /= dof ;
      if (dof > 1)
        var = (var - dof*mean*mean) / (dof-1);
      else
        var = (var / dof - mean*mean) ;
      if (var < 0 && dof > 2)
      {
        FILE *fp = fopen("val.log", "w") ;

        for (vno = 0 ; vno < mris->nvertices ; vno++)
        {
          v = &mris->vertices[vno] ;
          if (v->marked != parcel)
            continue ;
          val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
          fprintf(fp, "%f\n", val) ;
        }
        fclose(fp) ;
        DiagBreak() ;
      }
      // don't normalize by dofs
      //      MRIsetVoxVal(mri_means, parcel, 0, 0, frame, mean) ;
      //      MRIsetVoxVal(mri_vars, parcel, 0, 0, frame, var) ;
      mean_total += mean ; var_total += var ;
    }
    mean_total /= nframes ;  var_total /= nframes ;
    MRIsetVoxVal(mri_stats, parcel, 0,0, 0, mean_total) ;
    MRIsetVoxVal(mri_stats, parcel, 0,0, 1, var_total) ;
  }
  return(NO_ERROR) ;
}
static int
write_snapshot(MRI_SURFACE *mris, PARMS *parms, int n, MRI *mri_cmatrix)
{
  char fname[STRLEN] ;
  static MRI  *mri_timecourses = NULL ;
  int         vno, parcel, frame, nframes, min_vno ;
  double      val ;

  int req = snprintf(fname, STRLEN, "%s.%3.3d.annot", parms->base_name, n) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  printf("writing snapshot to %s\n", fname) ;
  MRISwriteAnnotation(mris, fname) ;

  if (DIAG_VERBOSE_ON)
  {
    nframes = mri_cmatrix->nframes ;
    if (mri_timecourses == NULL)
      mri_timecourses = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, nframes) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      parcel = mris->vertices[vno].marked ;
      min_vno = parms->stats.min_vno[parcel] ;
      for (frame = 0 ; frame < nframes ; frame++)
      {
        val = MRIgetVoxVal(mri_cmatrix, min_vno, 0, 0, frame) ;
        MRIsetVoxVal(mri_timecourses, frame, 0, 0, vno, val) ;
      }
    }
    int req = snprintf(fname, STRLEN, "%s.%3.3d.mgz", parms->base_name, n) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwrite(mri_timecourses, fname) ;
  }
  return(NO_ERROR) ;
}

static double
markov_energy(MRI_SURFACE *mris)
{
  int    vno, nborders, n ;
  double energy ;

  for (energy = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (nborders = n = 0 ; n < vt->vnum ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->marked != v->marked)
        nborders++ ;
    }
#define EXP_K 10
    energy += exp(EXP_K*(double)nborders/(double)vt->vnum)-1 ;
  }
  energy /= mris->nvertices ;
  return(energy) ;
}

static int
update_parcellation_statistics(MRI_SURFACE *mris, int vno, int old_parcel, int new_parcel,
                               MRI *mri_cmatrix, MRI *mri_means, MRI *mri_vars, 
                               MRI *mri_stats, PARMS *parms)
{
  int    frame, dofs, nframes ;
  double mean, var, val, total_mean, total_var ;

  nframes = mri_means->nframes ;
  dofs = MRIgetVoxVal(mri_stats, old_parcel, 0, 0, 2) - 1 ;
  MRIsetVoxVal(mri_stats, old_parcel, 0, 0, 2, dofs) ;
  if (old_parcel == Gdiag_no)
    printf("removing %d from parcel %d, dofs = %d, ", vno, old_parcel, dofs) ;

  parms->stats.num[old_parcel]-- ;
  parms->stats.num[new_parcel]++ ;
  for (total_mean = total_var = 0.0, frame = 0 ; frame < nframes ; frame++)
  {
    val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;

    mean = MRIgetVoxVal(mri_means, old_parcel, 0, 0, frame) ;
    var = MRIgetVoxVal(mri_vars, old_parcel, 0, 0, frame) ;
    mean -= val ;
    var -= (val*val) ;
    if (old_parcel == 133 && frame == 322 && DIAG_VERBOSE_ON)
      printf("removing vno %d: %f from parcel %d, now %f +- %f (var=%f), dofs=%d\n",
             vno, val, old_parcel, mean, var,
             (var - (mean*mean/dofs)) / (dofs-1),dofs);
            
    MRIsetVoxVal(mri_means, old_parcel, 0, 0, frame, mean) ;
    MRIsetVoxVal(mri_vars, old_parcel, 0, 0, frame, var) ;
    mean /= dofs ;
    if (dofs > 1)
      var = (var - dofs*mean*mean) / (dofs-1);
    else
      var = (var /  dofs - mean*mean) ;
    if (dofs > 2 && var < 0)
      DiagBreak() ;
    total_mean += mean ;
    total_var += var ;
    if (old_parcel == Gdiag_no && frame == 2 && Gdiag_no >= 0)
      printf("\tval = %2.3f, mean = %2.3f, var = %2.5f\n", val, mean, var) ;
  }
  total_mean /= nframes ; total_var /= nframes ;
  MRIsetVoxVal(mri_stats, old_parcel, 0, 0, 0, total_mean) ;
  MRIsetVoxVal(mri_stats, old_parcel, 0, 0, 1, total_var) ;
  
  dofs = MRIgetVoxVal(mri_stats, new_parcel, 0, 0, 2)+1 ;
  MRIsetVoxVal(mri_stats, new_parcel, 0, 0, 2, dofs) ;
  if (new_parcel == Gdiag_no)
    printf("adding %d to parcel %d, dofs = %d, ", vno, new_parcel, dofs) ;
  for (total_mean = total_var = 0.0, frame = 0 ; frame < nframes ; frame++)
  {
    val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;

    mean = MRIgetVoxVal(mri_means, new_parcel, 0, 0, frame) ;
    var = MRIgetVoxVal(mri_vars, new_parcel, 0, 0, frame) ;
    mean += val ;
    var += (val*val) ;
    if (new_parcel == 133 && frame == 322 && DIAG_VERBOSE_ON)
      printf("contributing vno %d: %f to parcel %d, now %f +- %f (var=%f), dofs=%d\n",
             vno, val,new_parcel, mean, var,
             (var - (mean*mean/dofs)) / (dofs-1),dofs);
    MRIsetVoxVal(mri_means, new_parcel, 0, 0, frame, mean) ;
    MRIsetVoxVal(mri_vars, new_parcel, 0, 0, frame, var) ;
    mean /= dofs ;
    if (dofs > 1)
      var = (var - dofs*mean*mean) / (dofs-1);
    else
      var = (var /  dofs - mean*mean) ;
    if (dofs > 2 && var < 0)
      DiagBreak() ;
    total_mean += mean ;
    total_var += var ;
    if (new_parcel == Gdiag_no && frame == 2 && Gdiag_no >= 0)
      printf("\tval = %2.3f, mean = %2.3f, var = %2.5f\n", val, mean, var) ;
  }
  total_mean /= nframes ; total_var /= nframes ;
  MRIsetVoxVal(mri_stats, new_parcel, 0, 0, 0, total_mean) ;
  MRIsetVoxVal(mri_stats, new_parcel, 0, 0, 1, total_var) ;
  MRIsetVoxVal(mri_stats, new_parcel, 0, 0, 2, dofs) ;
  mris->vertices[vno].marked = new_parcel ;
  return(NO_ERROR) ;
}

static double
gaussian_energy(MRI_SURFACE *mris, PARMS *parms)
{
  int    vno, parcel, nparcels ;
  VERTEX *v ;
  double total_energy, radius, circumference, dist, norm, dot ;

  nparcels = parms->stats.nparcels ;
  memset(parms->stats.energy, 0, sizeof(double)*nparcels) ;

  v = &mris->vertices[0] ;
  radius = sqrt(v->cx*v->cx + v->cy*v->cy + v->cz*v->cz) ;
  circumference = 2*M_PI*radius ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    parcel = v->marked ;
    // use euclidean distance for now
    dist = sqrt(SQR(v->cx - parms->stats.cx[parcel]) + 
                SQR(v->cy - parms->stats.cy[parcel]) + 
                SQR(v->cz - parms->stats.cz[parcel])) ;
    norm = sqrt(SQR(v->cx) + SQR(v->cy) + SQR(v->cz)) ;
    dot = fabs(((v->cx)*parms->stats.cx[parcel] + 
               (v->cy)*parms->stats.cy[parcel] + 
                (v->cz)*parms->stats.cz[parcel])/norm) ;
#if 0
    dist = acos(dot) * circumference / (2*M_PI) ;
#endif
    dist /= parms->gaussian_std ;
    parms->stats.energy[parcel] += dist*dist ;
  }
  for (nparcels = 0,total_energy = 0.0, parcel = 0; parcel < parms->stats.nparcels; parcel++)
  {
    if (parcel == Gdiag_no)
      DiagBreak() ;
    if (parms->stats.num[parcel] > 0)
    {
      parms->stats.energy[parcel] /= parms->stats.num[parcel] ;
      total_energy += parms->stats.num[parcel] * parms->stats.energy[parcel] ;
      nparcels++ ;
    }
  }
  total_energy /= mris->nvertices ;

  return(total_energy) ;
}
static double
compute_parcellation_median_energy(MRI_SURFACE *mris, MRI *mri_stats, PARMS *parms, 
                                   MRI *mri_means, MRI *mri_cmatrix)
{
  double  energy, var_within, var_between, mean, mean_nbr, var_between_total, val, 
    var_within_total, var_between_nbr ;
  int     **nbrs, parcel, vno, n, nparcels, n_nbrs, frame, nframes ;

  return(0.0) ;
  nframes = mri_cmatrix->nframes ;
  nparcels = mri_stats->width ;
  nbrs = (int **)calloc(nparcels, sizeof(int *)) ;
  for (parcel = 0 ; parcel < nparcels ; parcel++)
    nbrs[parcel] = (int *)calloc(nparcels, sizeof(int)) ;

  // compute within-cluster L1 distance
  for (var_within_total = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue;
    parcel = v->marked ;
    if (parcel < 0)
    {
      printf("warning, v %d is marked with parcel %d\n", vno, parcel) ;
      continue ;
    }
    for (n = 0 ; n < vt->v[n] ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      nbrs[parcel][vn->marked] = 1 ;
    }
    for (var_within = 0.0, frame = 0 ; frame < nframes ; frame++)
    {
      val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
      mean = MRIgetVoxVal(mri_means, parcel, 0, 0, frame) ;
      var_within += fabs(val-mean) ;
    }
    var_within_total += (var_within / nframes) ;
  }
  
  // compute between-cluster average L1 distance
  for (var_between_total = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX const * const v = &mris->vertices[vno] ; 
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue;
    parcel = v->marked ;
    if (parcel < 0)
    {
      printf("warning, v %d is marked with parcel %d\n", vno, parcel) ;
      continue ;
    }

    for (var_between = 0.0, n_nbrs = n = 0 ; n < nparcels ; n++)
    {
      if (n != parcel && nbrs[parcel][n] > 0)  // a neighboring parcel
      {
        n_nbrs++ ;
        for (var_between_nbr = 0.0, frame = 0 ; frame < nframes ; frame++)
        {
          val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
          mean_nbr = MRIgetVoxVal(mri_means, n, 0, 0, frame) ;
          var_between_nbr += fabs(val-mean) ;
        }
        var_between += var_between_nbr / nframes ;
      }
    }
#if 0
    if (n_nbrs > 0)
      var_between_total += var_between / n_nbrs ;
#else
    var_between_total += var_between ;
#endif
  }
  energy = var_between_total / var_within_total ;

  for (parcel = 0 ; parcel < nparcels ; parcel++)
    free(nbrs[parcel]) ;
  free(nbrs) ; 

  return(energy) ;
}
static int
find_parcel_min_timecourses(MRI_SURFACE *mris, MRI *mri_cmatrix, PARMS *parms)
{
  int     nframes, parcel, vertices[MAX_PARCEL_VERTICES], nvertices, vno1, vno2, n1, n2, 
    min_vno ;
  double  dist, min_dist, avg_dist, dmat[MAX_PARCEL_VERTICES][MAX_PARCEL_VERTICES] ;

  nframes = mri_cmatrix->nframes ;
  for (parcel = 0 ; parcel < parms->stats.nparcels ; parcel++)
  {
    for (nvertices = vno1 = 0 ; vno1 < mris->nvertices ; vno1++)
      if (mris->vertices[vno1].marked == parcel)
        vertices[nvertices++] = vno1 ;    // make a list of all vertices in this parcel

    for (n1 = 0 ; n1 < nvertices ; n1++)
    {
      dmat[n1][n1] = 0 ;
      for (n2 = n1+1 ; n2 < nvertices ; n2++)
      {
        vno1 = vertices[n1] ; vno2 = vertices[n2] ;
        //        dist = L1_distance(mri_cmatrix, mri_cmatrix, vno1, 0,0, vno2, 0, 0) ;
        dist = MRIgetVoxVal(parms->mri_dmat, vno1, vno2, 0, 0) ;
        dmat[n1][n2] = dist ; dmat[n2][n1] = dist ;
      }
    }
    min_dist = mri_cmatrix->nframes * 1e10; min_vno = 0 ;
    for (n1 = 0 ; n1 < nvertices ; n1++)
    {
      for (avg_dist = 0.0, n2 = 0 ; n2 < nvertices ; n2++)
        avg_dist += dmat[n1][n2] ;
      avg_dist /= nvertices ;
      if (avg_dist < min_dist)
      {
        min_dist = avg_dist ;
        min_vno = vertices[n1] ;
      }
    }
    parms->stats.min_vno[parcel] = min_vno ;
  }
      
    
  return(NO_ERROR) ;
}

static int
allocate_stats(MRI_SURFACE *mris, PARMS *parms, int nparcels, MRI *mri_cmatrix)
{
  int      n, vno, parcel, vno2 ;
  double   L1_dist ;

  parms->stats.nparcels = nparcels ;
  parms->stats.cx = (double *)calloc(nparcels, sizeof(double)) ;
  parms->stats.cy = (double *)calloc(nparcels, sizeof(double)) ;
  parms->stats.cz = (double *)calloc(nparcels, sizeof(double)) ;
  parms->stats.energy = (double *)calloc(nparcels, sizeof(double)) ;
  parms->stats.num = (int *)calloc(nparcels, sizeof(int)) ;
  parms->stats.min_vno = (int *)calloc(nparcels, sizeof(int)) ;
  parms->stats.nbrs = (int **)calloc(nparcels, sizeof(int *)) ;
  for (parcel = 0 ; parcel < nparcels ; parcel++)
    parms->stats.nbrs[parcel] = (int *)calloc(nparcels, sizeof(int)) ;
  parms->mri_dmat = MRIalloc(mris->nvertices, mris->nvertices, 1, MRI_FLOAT) ;
  if (parms->mri_dmat == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate complete distance matrix", Progname) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue;
    parcel = v->marked ;
    for (n = 0 ; n < vt->v[n] ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      parms->stats.nbrs[parcel][vn->marked] = 1 ;
    }
  }
  
  if (parms->energy_type == ENERGY_VARIANCE)
  {
    parms->stats.var_within = (double *)calloc(nparcels, sizeof(double)) ;
    if (parms->stats.var_within == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate var within", Progname) ;
    parms->stats.var_between = (double *)calloc(nparcels, sizeof(double)) ;
    if (parms->stats.var_between == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate var within", Progname) ;
  }

  if (parms->read_dmat)
  {
    MRI *mri_norm ;
    double norm1, norm2 ;

    parms->mri_dmat = MRIread(parms->dmat_read_fname) ;
    if (parms->mri_dmat == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read distance matrix volume from %s",
                Progname, parms->dmat_read_fname) ;
    mri_norm = MRIcomputeFrameVectorLength(mri_cmatrix, NULL) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      for (vno2 = vno+1 ; vno2 < mris->nvertices ; vno2++)
      {
        L1_dist = MRIgetVoxVal(parms->mri_dmat, vno, vno2, 0, 0) ;
        norm1 = MRIgetVoxVal(mri_norm, vno, 0, 0, 0) ;
        norm2 = MRIgetVoxVal(mri_norm, vno2, 0, 0, 0) ;
        if (FZERO(norm1))
          norm1 = 1 ;
        if (FZERO(norm2))
          norm2 = 1 ;
        L1_dist /= sqrt(norm1*norm2) ;
        MRIsetVoxVal(parms->mri_dmat, vno, vno2, 0, 0, L1_dist) ;
        MRIsetVoxVal(parms->mri_dmat, vno2, vno, 0, 0, L1_dist) ;
      }
    }
    MRIfree(&mri_norm) ;
  }
  else
  {
    if (parms->energy_type != ENERGY_VARIANCE)
    {
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        for (vno2 = vno+1 ; vno2 < mris->nvertices ; vno2++)
        {
          L1_dist = L1_distance(mri_cmatrix, mri_cmatrix, vno, 0, 0, vno2, 0, 0) ;
          MRIsetVoxVal(parms->mri_dmat, vno, vno2, 0, 0, L1_dist) ;
          MRIsetVoxVal(parms->mri_dmat, vno2, vno, 0, 0, L1_dist) ;
        }
      }
      if (parms->write_dmat)
        MRIwrite(parms->mri_dmat, parms->dmat_write_fname) ;
    }
  }
  return(NO_ERROR) ;
}

static double
L1_distance(MRI *mri1, MRI *mri2, int x1, int y1, int z1, int x2, int y2, int z2)
{
  double   L1, val1, val2 ;
  int      f ;

  for (L1 = 0.0, f = 0 ; f < mri1->nframes ; f++)
  {
    val1 = MRIgetVoxVal(mri1, x1, y1, z1, f) ;
    val2 = MRIgetVoxVal(mri2, x2, y2, z2, f) ;
    L1 += fabs(val1-val2) ;
  }
  
  return(L1/mri1->nframes) ;
}


static double
compute_parcellation_energy_change(MRI_SURFACE *mris, PARMS *parms, 
                                   MRI *mri_cmatrix, MRI *mri_means, MRI *mri_vars, int vno, 
                                   int parcel_to_move_to,
                                   double *dold, double *dnew,
                                   double *nold, double *nnew)
{
  double  dist_within_old, dist_between_old, dist_within_new, dist_between_new, dist, 
    dist_between, markov_old, markov_new, energy_old, energy_new, gauss_old, gauss_new, 
    area_old, area_new, border_old, border_new;
  int     parcel, n, nparcels, n_nbrs, nframes ;
  VERTEX  *v ;

  dist_between_old = dist_between_new = 
    energy_new = energy_old = dist_within_old = dist_within_new = 0.0 ;
  nparcels = parms->stats.nparcels ;

  v = &mris->vertices[vno] ;
  if (vno == Gdiag_no)
    DiagBreak() ;
  parcel = v->marked ;
  nparcels = parms->stats.nparcels ; nframes = mri_cmatrix->nframes ;
  if (parms->l_var > 0)
  {
    if (parms->energy_type == ENERGY_VARIANCE)
    {
#if 0
      int frame, nframes  ;
      double val, mean, old_var_within, new_var_within ;

      // compute new means and vars for these two parcels
      for (frame = 0 ; frame < nframes ; frame++)
      {
        // remove it from old parcel
        val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
        mean = MRIgetVoxVal(mri_means, parcel, 0, 0, frame) ;
        MRIsetVoxVal(mri_means, parcel, 0, 0, mean-val) ;
        var = MRIgetVoxVal(mri_var, parcel, 0, 0, frame) ;
        MRIsetVoxVal(mri_var, parcel, 0, 0, var - (val*val)) ;

        // add it to new parcel
        mean = MRIgetVoxVal(mri_means, parcel_to_move_to, 0, 0, frame) ;
        MRIsetVoxVal(mri_means, parcel_to_move_to, 0, 0, mean+val) ;
        var = MRIgetVoxVal(mri_var, parcel_to_move_to, 0, 0, frame) ;
        MRIsetVoxVal(mri_var, parcel_to_move_to, 0, 0, var + (val*val)) ;
      }

      

      // compute new means and vars for these two parcels
      for (frame = 0 ; frame < nframes ; frame++)
      {
        // remove it from old parcel
        val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
        mean = MRIgetVoxVal(mri_means, parcel, 0, 0, frame) ;
        MRIsetVoxVal(mri_means, parcel, 0, 0, mean-val) ;
        var = MRIgetVoxVal(mri_var, parcel, 0, 0, frame) ;
        MRIsetVoxVal(mri_var, parcel, 0, 0, var - (val*val)) ;

        // add it to new parcel
        mean = MRIgetVoxVal(mri_means, parcel_to_move_to, 0, 0, frame) ;
        MRIsetVoxVal(mri_means, parcel_to_move_to, 0, 0, mean+val) ;
        var = MRIgetVoxVal(mri_var, parcel_to_move_to, 0, 0, frame) ;
        MRIsetVoxVal(mri_var, parcel_to_move_to, 0, 0, var + (val*val)) ;
      }


      for (old_var_within = 0.0, frame = 0 ; frame < nframes ; frame++)
      {
        val = MRIFseq_vox(mri_cmatrix, vno, 0, 0, frame) ;
        mean = MRIFseq_vox(mri_means, parcel, 0, 0, frame)/parms->stats.num[parcel] ;
        old_var_within += SQR(mean-val) ;
      }
      old_var_within /= nframes ;
      for (old_var_between = 0.0, n = 0 ; n < parms->mris_ico->vertices[parcel].vnum ; n++)
      {
        int nbr_parcel ;

        nbr_parcel = parms->mris_ico->vertices[parcel].v[n] ;
        for (var = 0.0, frame = 0 ; frame < nframes ; frame++)
        {
          mean = MRIFseq_vox(mri_means, parcel, 0, 0, frame)/parms->stats.num[parcel] ;
          val = MRIFseq_vox(mri_means, nbr_parcel, 0, 0, frame)/parms->stats.num[nbr_parcel] ;
          var += SQR(mean-val) ;
        }
        old_var_between += (var/frames) ;
      }
      old_var_between /= parms->mris_ico->vertices[parcel].vnum ;


      energy_old = parms->stats.energy_num / parms->stats.energy_den ;
      delta_var = new_within_var - old_within_var ;
      energy_new = (parms->stats.energy_num + delta_var) / parms->stats.nenergy_den ;
#else
      update_parcellation_statistics(mris, vno, parcel, parcel_to_move_to,
                                     mri_cmatrix, mri_means, mri_vars, parms->mri_stats, parms);
      energy_new = 
        compute_variance_energy_for_vertex_change(mris, parms, mri_cmatrix, mri_means, mri_vars, 
                                                  &dist_within_new, &dist_between_new, vno, parcel,
                                                  parcel_to_move_to) ;
      update_parcellation_statistics(mris, vno, parcel_to_move_to, parcel,
                                     mri_cmatrix, mri_means, mri_vars, parms->mri_stats, parms);
      energy_old = 
        compute_variance_energy_for_vertex_change(mris, parms, mri_cmatrix, mri_means, mri_vars, 
                                                  &dist_within_new, &dist_between_new, vno, parcel_to_move_to,
                                                  parcel) ;
#endif
    }
    else if (parms->energy_type == ENERGY_DISTANCE)
    {
      
      // compute energy for this vertex in current parcel
      //    dist = L1_distance(mri_cmatrix, parms->stats.mri_min_timecourse, vno, 0, 0, parcel, 0, 0) ;
      dist = MRIgetVoxVal(parms->mri_dmat, vno,parms->stats.min_vno[parcel],0,0);
      dist_within_old = dist ;
      
      for (dist_between = 0.0, n_nbrs = n = 0 ; n < nparcels ; n++)
        if (n != parcel && parms->stats.nbrs[parcel][n] > 0)
        {
          //        dist = L1_distance(mri_cmatrix, parms->stats.mri_min_timecourse, vno, 0, 0, n, 0, 0) ;
          dist = MRIgetVoxVal(parms->mri_dmat, vno, parms->stats.min_vno[n],0,0);
          dist_between += dist ;
          n_nbrs++ ;
        }
      
      dist_between_old = dist_between/n_nbrs ;
      energy_old = dist_within_old / dist_between_old ;
      
      // compute energy for this vertex in new parcel
      //    dist = L1_distance(mri_cmatrix, parms->stats.mri_min_timecourse, vno, 0, 0, parcel_to_move_to, 0, 0) ;
      dist = MRIgetVoxVal(parms->mri_dmat, vno, 
                          parms->stats.min_vno[parcel_to_move_to], 0, 0) ;
      dist_within_new = dist ;
      
      for (dist_between = 0.0, n_nbrs = n = 0 ; n < nparcels ; n++)
        if (n!= parcel_to_move_to && parms->stats.nbrs[parcel_to_move_to][n] > 0)
        {
          //        dist = L1_distance(mri_cmatrix, parms->stats.mri_min_timecourse, vno, 0, 0, n, 0, 0) ;
          dist = MRIgetVoxVal(parms->mri_dmat, vno, parms->stats.min_vno[n],0,0);
          dist_between += dist ;
          n_nbrs++ ;
        }
      
      dist_between_new = dist_between/n_nbrs ;
      energy_new = dist_within_new / dist_between_new ;
      
#if 1
      energy_new = (parms->stats.energy_num-dist_within_old+dist_within_new) / 
        (parms->l_eden*(parms->stats.energy_den-dist_between_old+dist_between_new)) ;
      energy_old = parms->stats.energy_num / (parms->l_eden * parms->stats.energy_den);
#endif
    }
  }
  else
    energy_old = energy_new = 0.0 ;

  //  energy_old = parms->l_var * energy_old / parms->stats.num[v->marked] ;
  energy_old = parms->l_var * energy_old ;
  if (parms->l_markov > 0)
  {
    markov_old = markov_energy(mris) ;
    energy_old +=  parms->l_markov * markov_old ;
  }
  if (parms->l_border > 0)
  {
    border_old = border_energy(mris) ;
    energy_old +=  parms->l_border * border_old ;
  }
  if (parms->l_area > 0)
  {
    area_old = area_energy(mris, parms) ;
    energy_old +=  parms->l_area * area_old ;
  }
  if (parms->l_gaussian > 0)
  {
    gauss_old = gaussian_energy(mris, parms) ;
    energy_old += parms->l_gaussian * gauss_old ;
  }

  parcel = v->marked; v->marked = parcel_to_move_to ;
  energy_new = parms->l_var * energy_new ;
  if (parms->l_markov > 0)
  {
    markov_new = markov_energy(mris) ;
    energy_new +=  parms->l_markov * markov_new ;
  }
  if (parms->l_border > 0)
  {
    border_new = border_energy(mris) ;
    energy_new +=  parms->l_border * border_new ;
  }
  if (parms->l_area > 0)
  {
    parms->stats.num[parcel]-- ;
    parms->stats.num[parcel_to_move_to]++ ;
    area_new = area_energy(mris, parms) ;
    parms->stats.num[parcel]++ ;
    parms->stats.num[parcel_to_move_to]-- ;
    energy_new +=  parms->l_area * area_new ;
  }
  if (parms->l_gaussian > 0)
  {
    gauss_new = gaussian_energy(mris, parms) ; ;
    energy_new += parms->l_gaussian * gauss_new ;
  }
  v->marked = parcel ;  // restore it

  if (dold)
    *dold = dist_between_old ;
  if (dnew)
    *dnew = dist_between_new ;
  if (nnew)
    *nnew = dist_within_new ;
  if (nold)
    *nold = dist_within_old ;
  return(energy_new-energy_old) ;
}

static MRI *
compute_similarity_matrix(MRI_SURFACE *mris, MRI *mri_cmatrix, MRI *mri_smatrix, 
                          double thresh, double max_dist)
{
  int     vno, n, vno2, vno3 ;
  VERTEX  *v, *v2 ;
  double  count, val1, val2, dist, norm1, norm2 ;

  printf("computing similarity matrix...\n") ;
  mri_smatrix = MRIclone(mri_cmatrix, mri_smatrix) ;
  for (vno = 0 ;  vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    //    for (n = 0; n < v->vtotal ; n++)
    for (n = 0; n < mris->nvertices ; n++)
    {
      //      vno2 = v->v[n] ;
      vno2 = n ;
      if (vno2 > vno)
        continue ;   // only do it once
      v2 = &mris->vertices[vno2] ;
      if (vno2 == Gdiag_no)
        DiagBreak() ;
      dist = sqrt(SQR(v->x-v2->x)+SQR(v->y-v2->y)+SQR(v->z-v2->z));
      if (dist > max_dist)
        continue ;

      // compute the similarity of vno and vno2 in all other maps
      for (norm1 = norm2 = count = 0.0, vno3 = 0 ;  vno3 < mris->nvertices ; vno3++)
      {
#if 0
        val1 = MRIgetVoxVal(mri_cmatrix, vno3, 0, 0, vno) ;
        val2 = MRIgetVoxVal(mri_cmatrix, vno3, 0, 0, vno2) ;
#else
        val1 = MRIFseq_vox(mri_cmatrix, vno3, 0, 0, vno) ;
        val2 = MRIFseq_vox(mri_cmatrix, vno3, 0, 0, vno2) ;
#endif
        if (val1 > thresh && val2 > thresh)
          count += (val1*val2) ;
        norm1 += fabs(val1) ;
        norm2 += fabs(val2) ;
      }
      count = sqrt(count/(norm1*norm2)) ;
      MRIsetVoxVal(mri_smatrix, vno, 0, 0, vno2, count) ;
      MRIsetVoxVal(mri_smatrix, vno2, 0, 0, vno, count) ;
    }
  }

  MRIwrite(mri_smatrix, "s.mgz");
  return(mri_smatrix) ;
}
static MRI *
compute_eta_matrix(MRI_SURFACE *mris, MRI *mri_cmatrix, MRI *mri_ematrix, 
                   double thresh, double max_dist)
{
  int     vno, n, vno2, vno3 ;
  VERTEX  *v, *v2 ;
  double  eta, val1, val2, dist, M, m, num, den ;
  MRI     *mri_mean ;

  printf("computing eta matrix...\n") ;
  mri_mean = MRImeanTimecourse(mri_cmatrix, NULL) ;
  
  mri_ematrix = MRIclone(mri_cmatrix, mri_ematrix) ;
  for (vno = 0 ;  vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    //    for (n = 0; n < v->vtotal ; n++)
    for (n = 0; n < mris->nvertices ; n++)
    {
      //      vno2 = v->v[n] ;
      vno2 = n ;
      if (vno2 > vno)
        continue ;   // only do it once
      v2 = &mris->vertices[vno2] ;
      if (vno2 == Gdiag_no)
        DiagBreak() ;
      dist = sqrt(SQR(v->x-v2->x)+SQR(v->y-v2->y)+SQR(v->z-v2->z));
      if (dist > max_dist)
        continue ;

      // compute the similarity of vno and vno2 in all other maps
      M = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, 0) ;
      M += MRIgetVoxVal(mri_cmatrix, vno2, 0, 0, 0) ;
      M /= 2 ;
      for (num = den = 0.0, vno3 = 0 ;  vno3 < mris->nvertices ; vno3++)
      {
        val1 = MRIFseq_vox(mri_cmatrix, vno3, 0, 0, vno) ;
        val2 = MRIFseq_vox(mri_cmatrix, vno3, 0, 0, vno2) ;
        m = (val1+val2)/2 ;
        num += SQR(val1-m) + SQR(val2-m) ;
        den += SQR(val1-M) + SQR(val2-M) ;
      }
      eta = 1 - num/den ;
      MRIsetVoxVal(mri_ematrix, vno, 0, 0, vno2, eta) ;
      MRIsetVoxVal(mri_ematrix, vno2, 0, 0, vno, eta) ;
    }
  }

  MRIfree(&mri_mean) ;
  MRIwrite(mri_ematrix, "eta.mgz");
  return(mri_ematrix) ;
}

#define MAX_PARCELS 20000
static int
segment_similarity_matrix(MRI_SURFACE *mris, MRI_SURFACE *mris_ico, 
                                     MRI *mri_smatrix, PARMS *parms)
{
  int      nchanged, min_border = 1, vno, vno2, n, max_parcel, iter, nparcels, parcel ;
  double   max_similarity, val, sims[MAX_PARCELS] ;

  nparcels = mris_ico->nvertices ;
  iter = 0 ;
  if (parms->write_iterations > 0)
    write_snapshot(mris, parms, 0, mri_cmatrix) ;
  do
  {
    nchanged = 0 ;

    mark_border_vertices(mris) ;

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->border < min_border)
        continue ;

      // initialize parcel similarities to 0
      sims[v->marked] = 0 ;
      for (n = 0 ; n < vt->vtotal ; n++)
      {
        vno2 = vt->v[n] ;
        parcel = mris->vertices[vno2].marked ;
        sims[parcel] = 0 ;
      }

      // compute neighboring parcel similarities 
      for (n = 0 ; n < vt->vtotal ; n++)
      {
        vno2 = vt->v[n] ;
        parcel = mris->vertices[vno2].marked ;
        val = MRIgetVoxVal(mri_smatrix, vno, 0, 0, vno2) ;
#if 0
        if (val > sims[parcel])
          sims[parcel] = val ;
#else
        sims[parcel] += val ;
#endif
        
      }

      // find max neighboring parcel similarity
      max_similarity = sims[v->marked]; max_parcel = v->marked;
      for (n = 0 ; n < vt->vtotal ; n++)
      {
        vno2 = vt->v[n] ; 
        VERTEX const * const vn = &mris->vertices[vno2] ;
        if (sims[vn->marked] > max_similarity)
        {
          max_similarity = sims[vn->marked] ;
          max_parcel = vn->marked ;
        }
      }
      if (max_parcel != v->marked)
      {
        int annot ;
        if (vno == Gdiag_no)
          DiagBreak() ;
        nchanged++ ;
        CTABannotationAtIndex(mris->ct, v->marked, &annot);
        v->annotation = annot ;
        v->marked = max_parcel ;
      }
    }

    iter++ ;
    printf("iter %3.3d: nchanged = %d\n", iter, nchanged) ;
    if (parms->write_iterations > 0 && ((iter % parms->write_iterations) == 0))
      write_snapshot(mris, parms, iter, mri_cmatrix) ;
  } while (nchanged > 0) ;

  return(NO_ERROR) ;
}

static int
mark_border_vertices(MRI_SURFACE *mris) 
{
  int     vno, n, nborder ;

  for (nborder = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    v->border = 0 ;
    for (n = 0 ; n < vt->vnum ; n++)
      if (mris->vertices[vt->v[n]].marked != v->marked)
      {
        if (v->border == 0)  // only count it once
          nborder++ ;
        v->border++ ;
      }
  }
  return(nborder) ;
}

 static double
 compute_variance_energy(MRI_SURFACE *mris, PARMS *parms, MRI *mri_cmatrix, MRI *mri_means, 
                         MRI *mri_vars, double *penergy_num, double *penergy_den)
{
  double  energy, var_within, var_between, var, var_between_total, var_within_total, val, mean ;
  int     parcel, n, nparcels, n_nbrs, frame, nframes, nbr_parcel ;

  nparcels = parms->stats.nparcels ; nframes = mri_cmatrix->nframes ;

  memset(parms->stats.energy, 0, sizeof(double)*nparcels) ;
  for ( var_within_total = var_between_total = 0.0, parcel = 0 ; parcel < nparcels ; parcel++)
  {
    if (parcel == Gdiag_no)
      DiagBreak() ;
    if (parms->stats.num[parcel] <= 0)
      continue ;
    for (var_within = 0.0, frame = 0 ; frame < nframes ; frame++)
    {
      mean = MRIgetVoxVal(mri_means, parcel, 0, 0, frame) ;
      var =  MRIgetVoxVal(mri_vars, parcel, 0, 0, frame) ;
      mean /= parms->stats.num[parcel] ; 
      if (parms->stats.num[parcel] > 1)
        var = (var - parms->stats.num[parcel] * mean*mean) / (parms->stats.num[parcel]-1)  ;
      else
        var = (var / (parms->stats.num[parcel]) - mean*mean) ;
      var_within += var ;
    }
    parms->stats.var_within[parcel] = parms->stats.num[parcel] * (var_within / nframes) ;
    var_within_total += parms->stats.num[parcel] * (var_within / nframes) ;

    for (var_between = 0.0, n_nbrs = n = 0 ; n < parms->mris_ico->vertices_topology[parcel].vnum ; n++)
    {
      nbr_parcel = parms->mris_ico->vertices_topology[parcel].v[n] ;
      for (var = 0.0, frame = 0 ; frame < nframes ; frame++)
      {
        mean = MRIFseq_vox(mri_means, parcel, 0, 0, frame)/parms->stats.num[parcel] ;
        val = MRIFseq_vox(mri_means, nbr_parcel, 0, 0, frame)/parms->stats.num[nbr_parcel] ;
        var += SQR(mean-val) ;
      }
      var /= nframes ;
      var_between += var ;
    }
    var_between_total += (var_between/ parms->mris_ico->vertices_topology[parcel].vnum ) ;
  }
  var_within_total /= mris->nvertices ;
  if (FZERO(var_between_total))
    var_between_total = 1.0 ;
  energy = var_within_total / (parms->l_eden * var_between_total) ;
  *penergy_num = var_within_total ;
  *penergy_den = var_between_total ;
  return(parms->stats.nparcels * energy) ;
}

static double
border_energy(MRI_SURFACE *mris)
{
  return((double)mark_border_vertices(mris)/(double)mris->nvertices) ;
}

static double
area_energy(MRI_SURFACE *mris, PARMS *parms)
{
  int    parcel ;
  double energy ;

  for (parcel = 0, energy = 0.0 ; parcel < parms->stats.nparcels ; parcel++)
    energy += SQR((parms->stats.num[parcel]-parms->avg_area)/parms->avg_area) ;
  return(energy/parms->stats.nparcels) ;
}

static double
compute_variance_energy_for_vertex_change(MRI_SURFACE *mris, PARMS *parms, 
                                          MRI *mri_cmatrix, MRI *mri_means, 
                                          MRI *mri_vars, double *penergy_num, 
                                          double *penergy_den, int vno,
                                          int old_parcel, int new_parcel)
{
  double  energy, var_within, var_between, var, var_between_total, var_within_total, val, mean ;
  int     parcel, n, nparcels, n_nbrs, frame, nframes, nbr_parcel ;

  nparcels = parms->stats.nparcels ; nframes = mri_cmatrix->nframes ;

  memset(parms->stats.energy, 0, sizeof(double)*nparcels) ;
  for ( var_within_total = var_between_total = 0.0, parcel = 0 ; parcel < nparcels ; parcel++)
  {
    if (parcel == Gdiag_no)
      DiagBreak() ;

    nbr_parcel = 0 ;
    for (n = 0 ; n < parms->mris_ico->vertices_topology[parcel].vnum ; n++)
      if (parms->mris_ico->vertices_topology[parcel].v[n] == old_parcel ||
          parms->mris_ico->vertices_topology[parcel].v[n] == new_parcel)
      {
        nbr_parcel = 1 ;
        break ;
      }
    if (parcel != old_parcel && parcel != new_parcel && nbr_parcel == 0)
    {
      var_within_total += parms->stats.var_within[parcel] ;
      var_between_total += parms->stats.var_between[parcel] ;
      continue ;
    }

    for (var_within = 0.0, frame = 0 ; frame < nframes ; frame++)
    {
      mean = MRIgetVoxVal(mri_means, parcel, 0, 0, frame) ;
      var =  MRIgetVoxVal(mri_vars, parcel, 0, 0, frame) ;
      mean /= parms->stats.num[parcel] ; 
      if (parms->stats.num[parcel] > 1)
        var = (var - parms->stats.num[parcel] * mean*mean) / (parms->stats.num[parcel]-1) ;
      else
        var = (var / parms->stats.num[parcel]) - mean*mean  ;
      var_within += var ;
    }
    var_within_total += parms->stats.num[parcel] * (var_within / nframes) ;
    parms->stats.var_within[parcel] = parms->stats.num[parcel] * (var_within / nframes) ;

    for (var_between = 0.0, n_nbrs = n = 0 ; n < parms->mris_ico->vertices_topology[parcel].vnum ; n++)
    {
      nbr_parcel = parms->mris_ico->vertices_topology[parcel].v[n] ;
      for (var = 0.0, frame = 0 ; frame < nframes ; frame++)
      {
        mean = MRIFseq_vox(mri_means, parcel, 0, 0, frame)/parms->stats.num[parcel] ;
        val = MRIFseq_vox(mri_means, nbr_parcel, 0, 0, frame)/parms->stats.num[nbr_parcel] ;
        var += SQR(mean-val) ;
      }
      var /= nframes ;
      var_between += var ;
    }
    var_between_total += (var_between/ parms->mris_ico->vertices_topology[parcel].vnum ) ;
    parms->stats.var_between[parcel] = (var_between/ parms->mris_ico->vertices_topology[parcel].vnum ) ;
  }
  var_within_total /= mris->nvertices ;
  if (FZERO(var_between_total))
    var_between_total = 1.0 ;
  energy = var_within_total / (parms->l_eden * var_between_total) ;
  *penergy_num = var_within_total ;
  *penergy_den = var_between_total ;
  return(parms->stats.nparcels * energy) ;
}

static int
write_annot_correlations(MRI_SURFACE *mris, MRI *mri_cmatrix, PARMS *parms, char *fname)
{
  MRI  *mri_corr, *mri_counts ;
  int  vno, source_parcel, target_parcel, frame, nframes, nparcels ;
  double corr, total_corr, num ;

  nparcels = parms->stats.nparcels ; nframes = mri_cmatrix->nframes ;
  mri_corr = MRIalloc(nparcels, nparcels, 1, MRI_FLOAT) ;
  mri_counts = MRIalloc(nparcels, nparcels, 1, MRI_FLOAT) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    source_parcel = mris->vertices[vno].marked ;
    for (frame = 0 ; frame < nframes ; frame++)
    {
      target_parcel = mris->vertices[frame].marked ;
      corr = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
      total_corr = MRIgetVoxVal(mri_corr, source_parcel, target_parcel, 0, 0) ;
      MRIsetVoxVal(mri_corr, source_parcel, target_parcel, 0, 0, total_corr + corr) ;
      num = MRIgetVoxVal(mri_counts,source_parcel, target_parcel, 0, 0) ;
      MRIsetVoxVal(mri_counts, source_parcel, target_parcel, 0, 0, num+1) ;
    }
  }

  for (source_parcel = 0 ; source_parcel < nparcels ; source_parcel++)
    for (target_parcel = 0 ; target_parcel < nparcels ; target_parcel++)
    {
      num = MRIgetVoxVal(mri_counts, source_parcel, target_parcel, 0, 0) ;
      if (num > 0)
      {
        corr = MRIgetVoxVal(mri_corr, source_parcel, target_parcel, 0, 0) ;
        MRIsetVoxVal(mri_corr, source_parcel, target_parcel, 0, 0, corr/num) ;
      }
    }

  MRIwrite(mri_corr, fname) ;
  MRIfree(&mri_corr) ;
  return(NO_ERROR) ;
}
