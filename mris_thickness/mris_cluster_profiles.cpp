/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "transform.h"

typedef struct {
  VECTOR *v_mean ;
  MATRIX *m_cov ;
  int    npoints ;
  int    vno ;
}
CLUSTER ;


//int MRISfindClosestVertex(MRI_SURFACE *mris1, MRI_SURFACE *mris2, int vno1) ;
static int clusters_connected(MRI_SURFACE *mris, int vno1, int vno2, int n) ;
static int okay_to_swap(MRI_SURFACE *mris, int vno, int c1, int c2);
static int rip_bad_vertices(MRI_SURFACE *mris, MRI *mri_profiles) ;
static int remove_vertex_from_cluster(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int c, int vno) ;
static int add_vertex_to_cluster(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int c, int vno) ;
static int find_most_likely_unmarked_vertex(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k, int *pcluster_no) ;
static int load_vals(MRI *mri_profile, VECTOR *v, int vno) ;
static int find_most_likely_vertex_to_swap(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k, int *pcluster_no) ;
static CLUSTER *MRISclusterAgglomerative(MRI_SURFACE *mris, MRI *mri_profiles,
    int k, char *start_fname, MRI_SURFACE *mris_ico) ;
static int rip_vertices_out_of_fov(MRI_SURFACE *mris, MRI *mri_profiles) ;
int main(int argc, char *argv[]) ;

static int mark_clusters(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k) ;
static int compute_cluster_statistics(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k) ;
CLUSTER *MRISclusterKMeans(MRI_SURFACE *mris, MRI *mri_profiles, int k,
                           char *start_fname, MRI_SURFACE *mris_ico);
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int initialize_kmeans_centers(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k) ;
static int initialize_cluster_centers(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k) ;
static int initialize_cluster_centers_with_ico(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, MRI_SURFACE *mris_ico) ;

CLUSTER *MRIScluster(MRI_SURFACE *mris, MRI *mri_profiles, int cluster_type,
                     int k, char *start_fname, MRI_SURFACE *mris_ico) ;

char *start_fname = NULL ;

static int max_iterations = 500000 ;
const char *Progname ;
#define MAX_LABELS 10000
static char *label_names[MAX_LABELS] ;
static int nlabels = 0 ;

static int nbhd_size = 2 ;
static char *sdir = NULL ;
static int num_erode = 0 ;
static int navgs = 0 ;
static int write_iters = 500 ;

#define K_MEANS        0
#define AGGLOMERATIVE  1

static int cluster_type = AGGLOMERATIVE ;
static int k = 50 ;

static char *ico_fname = NULL ;

int
main(int argc, char *argv[]) {
  char          **av, *surf_fname, *profile_fname, *seg_fname ;
  int           ac, nargs ;
  MRI_SURFACE   *mris, *mris_ico = NULL ;
  // LABEL         *label = NULL ;
  MRI           *mri_profiles ;
  CLUSTER       *ct ;

  setRandomSeed(10L) ;

  nargs = handleVersionOption(argc, argv, "mris_cluster_profiles");
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

  if (argc < 4)
    usage_exit() ;

  profile_fname = argv[1] ;
  surf_fname = argv[2] ;
  seg_fname = argv[3] ;
  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;
  MRISsetNeighborhoodSizeAndDist(mris, 2) ;

  if (MRISreadAnnotation(mris, "ad_aparc") != NO_ERROR) {
    if (MRISreadAnnotation(mris, "aparc") != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read annotation file ad_aparc", Progname) ;
  }

  printf("reading intensity profile volume from %s...\n", profile_fname) ;
  mri_profiles = MRIread(profile_fname) ;
  if (!mri_profiles)
    ErrorExit(ERROR_NOFILE, "%s: could not read intensity volume %s",
              Progname, profile_fname) ;
  rip_vertices_out_of_fov(mris, mri_profiles) ;
  rip_bad_vertices(mris, mri_profiles) ;
  MRISclearAnnotations(mris) ;

#if 0
  if (nlabels > 0) {
    int l ;
    char label_name[STRLEN] ;
    LABEL *ltotal = NULL ;

    for (l = 0 ; l < nlabels ; l++) {
      sprintf(label_name, "%s/%s/label/%s.%s.label", sdir, sname, hemi,label_names[l]) ;

      label = LabelRead(NULL, label_name) ;
      if (!label)
        ErrorExit(ERROR_NOFILE, "%s: could not read label file %s...\n", Progname,
                  label_name) ;
      if (num_erode > 0) {
        printf("eroding label %d times, npoints went from %d ", num_erode,label->n_points) ;
        LabelErode(label, mris, num_erode) ;
        printf("to %d ", label->n_points) ;
      }
      ltotal = LabelCombine(label, ltotal) ;
    }
    if (nlabels == 0)
      ltotal = LabelInFOV(mris, mri, MIN_BORDER_DIST) ;

    LabelRipRestOfSurfaceWithThreshold(ltotal, mris, thresh) ;
  }
#endif

  if (navgs > 0) {
    printf("smoothing profiles %d times\n", navgs) ;
    MRISsmoothFrames(mris, mri_profiles, navgs) ;
  }
  if (ico_fname) {
    mris_ico = MRISread(ico_fname) ;
    if (mris_ico == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not read icosahedron from %s...\n", Progname, ico_fname) ;
  }
  ct = MRIScluster(mris, mri_profiles, cluster_type, k, start_fname, mris_ico) ;
  printf("writing cortical intensity clusters to %s...\n", seg_fname) ;
  MRISwriteAnnotation(mris, seg_fname) ;
  {
    int    vno ;
    VERTEX *v ;
    int    c, i ;
    char   fname[STRLEN], ext[STRLEN] ;

    // write average profiles into mri_profiles and write it out
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      c = v->curv ;
      if (c < 0)
        continue ;
      for (i = 0 ; i < mri_profiles->nframes ; i++)
        MRIsetVoxVal(mri_profiles, vno, 0, 0, i, VECTOR_ELT(ct[c].v_mean,i+1)) ;
    }
    FileNameExtension(seg_fname, ext) ;
    FileNameRemoveExtension(seg_fname, fname) ;
    strcat(fname, "_cluster_avg.mgz") ;
    printf("writing average cluster profiles to %s...\n", fname) ;
    MRIwrite(mri_profiles, fname) ;
  }

  MRIfree(&mri_profiles) ;

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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "erode")) {
    num_erode = atoi(argv[2]) ;
    fprintf(stderr,  "eroding label %d times\n", num_erode) ;
    nargs = 1 ;
  } else if (!stricmp(option, "sdir")) {
    sdir = argv[2] ;
    nargs = 1 ;
  } else if (!stricmp(option, "ic")) {
    ico_fname = argv[2] ;
    printf("using icosahedral surface %s to initialize clusters\n", ico_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "start")) {
    start_fname = argv[2] ;
    nargs = 1 ;
    printf("using %s as initial clustering...\n", start_fname) ;
  } else switch (toupper(*option)) {
    case 'M':
      max_iterations = atoi(argv[2]) ;
      printf("setting max_iterations = %d\n", max_iterations) ;
      nargs = 1 ;
      break ;
    case 'K':
      k = atof(argv[2]) ;
      nargs = 1 ;
      printf("using k-means clustering with k=%d\n", k) ;
      break ;
    case 'L':
      label_names[nlabels] = argv[2] ;
      nargs = 1 ;
      printf("limiting profile calculation to label %s\n", label_names[nlabels]) ;
      nlabels++ ;
      break ;
    case 'W':
      write_iters = atoi(argv[2]) ;
      printf("setting write iterations to %d\n", write_iters) ;
      nargs = 1 ;
      break ;
    case 'N':
      nbhd_size = atoi(argv[2]) ;
      fprintf(stderr, "using neighborhood size=%d\n", nbhd_size) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1;
      printf("smoothing profiles %d times across space\n", navgs) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input profile file> <spherical surface> <output clustering>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program clusters the intensity profile of the cortical ribbon\n"
          "and writes the resulting measurement into a 'curvature' file\n"
          "<output file>.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-sdir %%s specifies the SUBJECTS_DIR \n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

#define MAX_CLUSTERS 10000
CLUSTER *
MRIScluster(MRI_SURFACE *mris, MRI *mri_profiles, int cluster_type, int k,
            char *start_fname, MRI_SURFACE *mris_ico) {
  CLUSTER  *cluster_table ;

  switch (cluster_type) {
  default:
  case K_MEANS:
    cluster_table = MRISclusterKMeans(mris, mri_profiles, k, start_fname, mris_ico) ;
    break ;
  case AGGLOMERATIVE:
    cluster_table = MRISclusterAgglomerative(mris, mri_profiles, k, start_fname, mris_ico) ;
    break ;
  }

  return(cluster_table) ;
}
static CLUSTER *
MRISclusterAgglomerative(MRI_SURFACE *mris, MRI *mri_profiles, int k,
                         char *start_fname, MRI_SURFACE *mris_ico) {
  int    i, nsamples, iter, done, vno, cluster ;
  char   fname[STRLEN] ;
  CLUSTER *ct ;

  if (start_fname) {
    VERTEX *v ;

    if (MRISreadAnnotation(mris, start_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read initial annotation file %s",
                Progname, start_fname) ;
    k = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->annotation == 0) {
        v->ripflag = 1;
        continue ;
      }
      CTABfindAnnotation(mris->ct, v->annotation, &cluster);
      if (cluster >= k)
        k = cluster+1 ;
      v->curv = cluster ;
    }
    printf("%d clusters detected...\n", k) ;
    ct = calloc(k, sizeof(CLUSTER)) ;
    if (!ct)
      ErrorExit(ERROR_BADPARM, "%s: could not allocate %d clusters", Progname, k) ;
    nsamples = mri_profiles->nframes ;
    for (i = 0 ; i < k ; i++) {
      ct[i].v_mean = VectorAlloc(nsamples, MATRIX_REAL) ;
      ct[i].m_cov = MatrixIdentity(nsamples, NULL) ;
      ct[i].npoints = 0 ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      add_vertex_to_cluster(mris, mri_profiles, ct, v->curv, vno);
    }
  } else {
    if (mris_ico)
      k = mris_ico->nvertices ;
    ct = calloc(k, sizeof(CLUSTER)) ;
    if (!ct)
      ErrorExit(ERROR_BADPARM, "%s: could not allocate %d clusters", Progname, k) ;
    nsamples = mri_profiles->nframes ;
    for (i = 0 ; i < k ; i++) {
      ct[i].v_mean = VectorAlloc(nsamples, MATRIX_REAL) ;
      ct[i].m_cov = MatrixIdentity(nsamples, NULL) ;
      ct[i].npoints = 0 ;
    }

    MRISsetCurvature(mris, -1) ;
    if (mris_ico)
      initialize_cluster_centers_with_ico(mris, mri_profiles, ct, mris_ico) ;
    else
      initialize_cluster_centers(mris, mri_profiles, ct, k) ;

    done = iter = 0 ;
    do {
      vno = find_most_likely_unmarked_vertex(mris, mri_profiles, ct, k, &cluster);
      if (vno < 0)
        break ;
      add_vertex_to_cluster(mris, mri_profiles, ct, cluster, vno);
      if (write_iters > 0 && ((iter % write_iters) == 0)) {
        sprintf(fname, "%s.clusters%6.6d.annot",
                mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", iter) ;
        printf("%6.6d: writing %s\n", iter, fname) ;
        MRISwriteAnnotation(mris, fname) ;
        sprintf(fname, "./%s.clusters%6.6d.indices",
                mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", iter) ;
        //   MRISwriteCurvature(mris, fname) ;
      }
      if (write_iters == 0 && ((iter % 5000) == 0))
        printf("%6.6d of %6.6d\n", iter, mris->nvertices-k) ;
      if (iter++ > mris->nvertices-k || iter > max_iterations)
        done = 1 ;
    } while (!done) ;
  }
  iter = done = 0 ;
  do {
    vno = find_most_likely_vertex_to_swap(mris, mri_profiles, ct, k, &cluster);
    if (vno < 0)
      break ;
    remove_vertex_from_cluster(mris, mri_profiles, ct, mris->vertices[vno].curv, vno) ;
    add_vertex_to_cluster(mris, mri_profiles, ct, cluster, vno);
    if (write_iters > 0 && ((iter % write_iters) == 0)) {
      sprintf(fname, "%s.more_clusters%6.6d.annot",
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", iter) ;
      printf("%6.6d: writing %s\n", iter, fname) ;
      MRISwriteAnnotation(mris, fname) ;
      sprintf(fname, "./%s.more_clusters%6.6d.indices",
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", iter) ;
      //   MRISwriteCurvature(mris, fname) ;
    }
    if (write_iters == 0 && ((iter % 5000) == 0))
      printf("%6.6d of %6.6d\n", iter, mris->nvertices) ;
    if (iter++ > mris->nvertices || iter > max_iterations)
      done = 1 ;
  } while (!done) ;
  return(ct);
}
CLUSTER *
MRISclusterKMeans(MRI_SURFACE *mris, MRI *mri_profiles, int k, char *start_fname, MRI_SURFACE *mris_ico) {
  int    i, nsamples, iter, done, nchanged ;
  char   fname[STRLEN] ;
  CLUSTER *ct ;

  nsamples = mri_profiles->nframes ;
  ct = calloc(k, sizeof(CLUSTER)) ;
  for (i = 0 ; i < k ; i++) {
    ct[i].v_mean = VectorAlloc(nsamples, MATRIX_REAL) ;
    ct[i].m_cov = MatrixIdentity(nsamples, NULL) ;
    ct[i].npoints = 0 ;
  }

  initialize_kmeans_centers(mris, mri_profiles, ct, k) ;

  done = iter = 0 ;
  do {
    nchanged = mark_clusters(mris, mri_profiles, ct, k) ;
    if (nchanged == 0)
      done = 1 ;
    compute_cluster_statistics(mris, mri_profiles, ct, k) ;
    sprintf(fname, "%s.clusters%6.6d.annot",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", iter) ;
    printf("%6.6d: writing %s\n", iter, fname) ;
    MRISwriteAnnotation(mris, fname) ;
    sprintf(fname, "./%s.clusters%6.6d.indices",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", iter) ;
    MRISwriteCurvature(mris, fname) ;
    if (iter++ > max_iterations)
      done = 1 ;
  } while (!done) ;

  return(ct) ;
}
static int
initialize_kmeans_centers(MRI_SURFACE *mris,  MRI *mri_profiles, CLUSTER *ct, int k) {
  int i, j, done, vno, nsamples, vnos[MAX_CLUSTERS] ;
  double dist, min_dist ;
  VERTEX *v, *vn ;

  nsamples = mri_profiles->nframes ;
  min_dist = sqrt(mris->total_area/k) / 2 ;
  for (i = 0 ; i < k ; i++) {
    do {
      do {
        vno = (int)randomNumber(0, mris->nvertices-1) ;
        v = &mris->vertices[vno] ;
      } while (v->ripflag) ;

      // see if any of the other centers are too close to this one
      done = 1 ;
      for (j = 0 ; done && j < i ; j++) {
        if (j == i)
          continue ;
        vn = &mris->vertices[vnos[j]] ;
        dist = sqrt(SQR(vn->x-v->x)+SQR(vn->y-v->y)+SQR(vn->z-v->z)) ;
        if (dist < min_dist) {
          done = 0 ;
          break ;
        }
      }
    } while (!done) ;
    vnos[i] = vno ;
    for (j = 0 ; j < nsamples ; j++)
      VECTOR_ELT(ct[i].v_mean, j+1) = MRIgetVoxVal(mri_profiles, vno, 0, 0, j) ;
  }
  for (i = 0 ; i < k ; i++) {
    mris->vertices[vnos[i]].curv = i ;
    CTABannotationAtIndex(mris->ct, i, &mris->vertices[vnos[i]].annotation) ;
  }
  return(NO_ERROR) ;
}


static int
mark_clusters(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k) {
  int    i, vno, min_i, nsamples, num[MAX_CLUSTERS], nchanged ;
  double min_dist, dist ;
  VECTOR *v1 ;
  VERTEX *v ;

  memset(num, 0, sizeof(num)) ;
  nsamples = mri_profiles->nframes ;
  v1 = VectorAlloc(nsamples, MATRIX_REAL) ;
  for (nchanged = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    for (i = 0 ; i < nsamples ; i++)
      VECTOR_ELT(v1, i+1) = MRIgetVoxVal(mri_profiles, vno, 0, 0, i) ;
    min_dist = MatrixMahalanobisDistance(ct[0].v_mean, ct[0].m_cov, v1) ;
    min_i = 0 ;
    for (i = 1 ; i < k ; i++) {
      if (i == Gdiag_no)
        DiagBreak() ;
      dist = MatrixMahalanobisDistance(ct[i].v_mean, ct[i].m_cov, v1) ;
      if (dist < min_dist) {
        min_dist = dist ;
        min_i = i ;
      }
    }
    CTABannotationAtIndex(mris->ct, min_i, &v->annotation) ;
    if (v->curv != min_i)
      nchanged++ ;
    v->curv = min_i ;
    num[min_i]++ ;
  }
  for (i = 0 ; i < k ; i++)
    if (num[i] == 0)
      DiagBreak() ;
  VectorFree(&v1) ;
  printf("%d vertices changed clusters...\n", nchanged) ;
  return(nchanged) ;
}

// compute means and inverse covariances (in m_covs) for each cluster
#if 1
static int
compute_cluster_statistics(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k) {
  int    i, vno, cluster, nsamples, num[MAX_CLUSTERS];
  VECTOR *v1 ;

  memset(num, 0, sizeof(num)) ;
  nsamples = mri_profiles->nframes ;

  v1 = VectorAlloc(nsamples, MATRIX_REAL) ;

  for (cluster = 0 ; cluster < k ; cluster++)
    VectorClear(ct[cluster].v_mean) ;

  // compute means
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    cluster = mris->vertices[vno].curv ;
    for (i = 0 ; i < nsamples ; i++) {
      VECTOR_ELT(v1, i+1) = MRIgetVoxVal(mri_profiles, vno, 0, 0, i) ;
    }
    num[cluster]++ ;
    VectorAdd(ct[cluster].v_mean, v1, ct[cluster].v_mean) ;
  }

  for (cluster = 0 ; cluster < k ; cluster++)
    if (num[cluster] > 0)
      VectorScalarMul(ct[cluster].v_mean, 1.0/(double)num[cluster], ct[cluster].v_mean) ;

  VectorFree(&v1) ;
  return(NO_ERROR) ;
}
#else
static int
compute_cluster_statistics(MRI_SURFACE *mris, MRI *mri_profiles, MATRIX **m_covs, VECTOR **v_means, int k) {
  int    i, vno, cluster, nsamples, num[MAX_CLUSTERS];
  int    singular, cno_pooled, cno ;
  MATRIX *m1, *mpooled, *m_inv_covs[MAX_CLUSTERS] ;
  VECTOR *v1 ;
  FILE   *fp ;
  double det, det_pooled ;

  memset(num, 0, sizeof(num)) ;
  nsamples = mri_profiles->nframes ;

  v1 = VectorAlloc(nsamples, MATRIX_REAL) ;
  m1 = MatrixAlloc(nsamples, nsamples, MATRIX_REAL) ;
  mpooled = MatrixAlloc(nsamples, nsamples, MATRIX_REAL) ;

  for (cluster = 0 ; cluster < k ; cluster++) {
    VectorClear(v_means[cluster]) ;
    MatrixClear(m_covs[cluster]) ;
  }

  // compute means
  // fp = fopen("co.dat", "w") ;
  fp = NULL ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    cluster = mris->vertices[vno].curv ;
    for (i = 0 ; i < nsamples ; i++) {
      VECTOR_ELT(v1, i+1) = MRIgetVoxVal(mri_profiles, vno, 0, 0, i) ;
      if (cluster == 0 && fp)
        fprintf(fp, "%f ", VECTOR_ELT(v1, i+1));
    }
    if (cluster == 0 && fp)
      fprintf(fp, "\n") ;
    num[cluster]++ ;
    VectorAdd(v_means[cluster], v1, v_means[cluster]) ;
  }

  if (fp)
    fclose(fp) ;
  for (cluster = 0 ; cluster < k ; cluster++)
    if (num[cluster] > 0)
      VectorScalarMul(v_means[cluster], 1.0/(double)num[cluster], v_means[cluster]) ;

  // compute inverse covariances
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    cluster = mris->vertices[vno].curv ;
    for (i = 0 ; i < nsamples ; i++)
      VECTOR_ELT(v1, i+1) = MRIgetVoxVal(mri_profiles, vno, 0, 0, i) ;
    VectorSubtract(v_means[cluster], v1, v1) ;
    VectorOuterProduct(v1, v1, m1) ;
    MatrixAdd(m_covs[cluster], m1, m_covs[cluster]) ;
    MatrixAdd(mpooled, m1, mpooled) ;
  }

  MatrixScalarMul(mpooled, 1.0/(double)mris->nvertices, mpooled) ;
  cno_pooled = MatrixConditionNumber(mpooled) ;
  det_pooled = MatrixDeterminant(mpooled) ;
  for (cluster = 0 ; cluster < k ; cluster++)
    if (num[cluster] > 0)
      MatrixScalarMul(m_covs[cluster], 1.0/(double)num[cluster], m_covs[cluster]) ;


  // invert all the covariance matrices
  MatrixFree(&m1) ;
  singular = 0 ;
  for (cluster = 0 ; cluster < k ; cluster++) {
    m1 = MatrixInverse(m_covs[cluster], NULL) ;
    cno = MatrixConditionNumber(m_covs[cluster]) ;
    det = MatrixDeterminant(m_covs[cluster]) ;
    if (m1 == NULL)
      singular++ ;
    while (cno > 100*cno_pooled || 100*det < det_pooled) {
      if (m1)
        MatrixFree(&m1) ;
      m1 = MatrixScalarMul(mpooled, 0.1, NULL) ;
      MatrixAdd(m_covs[cluster], m1, m_covs[cluster]) ;
      MatrixFree(&m1) ;
      cno = MatrixConditionNumber(m_covs[cluster]) ;
      m1 = MatrixInverse(m_covs[cluster], NULL) ;
      det = MatrixDeterminant(m_covs[cluster]) ;
    }
    m_inv_covs[cluster] = m1 ;
  }

  for (cluster = 0 ; cluster < k ; cluster++) {
    if (m_inv_covs[cluster] == NULL)
      DiagBreak() ;
    else {
      MatrixFree(&m_covs[cluster]) ;
      m_covs[cluster] = m_inv_covs[cluster] ;
      //   MatrixIdentity(m_covs[cluster]->rows, m_covs[cluster]);
    }
  }
  MatrixFree(&mpooled) ;
  VectorFree(&v1) ;
  return(NO_ERROR) ;
}
#endif
static int
rip_vertices_out_of_fov(MRI_SURFACE *mris, MRI *mri_profiles) {
  int  vno, n, i, good, nsamples ;
  VERTEX *v ;

  nsamples = mri_profiles->nframes ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    for (good = i = 0 ; good == 0 && i < nsamples ; i++)
      if (!FZERO(MRIgetVoxVal(mri_profiles, vno, 0, 0, i)))
        good = 1 ;
    if (!good) {
      v->ripflag = 1 ;
      v->annotation = 0 ;
    }
  }

  // rip stuff around the identically 0 ones, as they can't be trusted
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag != 1)
      continue ;
    for (n = 0 ; n < v->vtotal ; n++) {
      mris->vertices[v->v[n]].ripflag = 2 ;
      mris->vertices[v->v[n]].annotation = 0 ;
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag == 2)
      v->ripflag = 1 ;
  }
  return(NO_ERROR) ;
}
static int
find_most_likely_unmarked_vertex(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k, int *pcluster_no) {
  double   dist, min_dist ;
  int      vno, c, min_cluster, min_vno, n, okay ;
  VECTOR   *v_vals ;
  VERTEX   *v, *vn ;

  v_vals = VectorAlloc(mri_profiles->nframes, MATRIX_REAL) ;

  min_dist = -1 ;
  min_vno = -1 ;
  min_cluster = -1 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->curv < 0)   // not in a cluster yet
      continue ;

    c = v->curv ;
    for (n = 0 ; n < v->vnum ; n++) {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->curv >= 0)
        continue ;
      if (vn->ripflag)
        continue ;
      load_vals(mri_profiles, v_vals, v->v[n]) ;
      dist = MatrixMahalanobisDistance(ct[c].v_mean, ct[c].m_cov, v_vals);
      if (min_vno < 0 || dist < min_dist) {
        if (ct[c].npoints == 1)
          okay = 1 ;
        else // check to make sure it has at least 2 nbrs of the right class
        {
          int n2, nc ;
          for (n2 = nc = 0 ; n2 < vn->vnum ; n2++)
            if (mris->vertices[vn->v[n2]].curv == c)
              nc++ ;
          okay = nc > 1 ;
        }
        if (okay) {
          min_cluster = c ;
          min_dist = dist ;
          min_vno = v->v[n] ;
        }
      }
    }
  }

  v = &mris->vertices[min_vno] ;
  *pcluster_no = min_cluster ;
  VectorFree(&v_vals) ;
  return(min_vno) ;
}
static int
find_most_likely_vertex_to_swap(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct,
                                int k, int *pcluster_no) {
  double   dist_current, dist_changed, dist_dec, max_dist_dec;
  int      vno, c1, c2, min_cluster, min_vno, n ;
  static VECTOR   *v_vals = NULL ;
  VERTEX   *v, *vn ;

  if (v_vals == NULL)
    v_vals = VectorAlloc(mri_profiles->nframes, MATRIX_REAL) ;

  // search for the vertex that has the biggest distance diff
  max_dist_dec = -1 ;
  min_vno = -1 ;
  min_cluster = -1 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->curv < 0)   // not in a cluster yet
      continue ;

    c1 = v->curv ;
    load_vals(mri_profiles, v_vals, vno) ;
    dist_current = MatrixMahalanobisDistance(ct[c1].v_mean, ct[c1].m_cov, v_vals);
    for (n = 0 ; n < v->vnum ; n++) {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->curv < 0)
        continue ;
      if (vn->ripflag)
        continue ;
      c2 = vn->curv ;
      if (c2 == c1)
        continue ;
      // distance to putative new cluster
      dist_changed = MatrixMahalanobisDistance(ct[c2].v_mean, ct[c2].m_cov, v_vals);
      dist_dec = (dist_current-dist_changed) ;
      if (min_vno < 0 || dist_dec > max_dist_dec) {
        if (okay_to_swap(mris, vno, c1, c2)) {
          min_cluster = c2 ;  // change it to new cluster
          max_dist_dec = dist_dec ;
          min_vno = vno ;
        }
      }
    }
  }

  v = &mris->vertices[min_vno] ;
  *pcluster_no = min_cluster ;
  return(min_vno) ;
}
static int
load_vals(MRI *mri_profiles, VECTOR *v, int vno) {
  int i ;

  for (i = 0 ; i < mri_profiles->nframes ; i++)
    VECTOR_ELT(v, i+1) = MRIgetVoxVal(mri_profiles, vno, 0, 0, i) ;
  return(NO_ERROR) ;
}

static int
remove_vertex_from_cluster(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct,
                           int c, int vno) {
  int    nsamples;
  static VECTOR *v_vals = NULL ;

  nsamples = mri_profiles->nframes ;

  mris->vertices[vno].curv = c ;
  if (v_vals == NULL)
    v_vals = VectorAlloc(nsamples, MATRIX_REAL) ;

  VectorScalarMul(ct[c].v_mean, ct[c].npoints, ct[c].v_mean) ;
  load_vals(mri_profiles, v_vals, vno) ;
  VectorSubtract(ct[c].v_mean, v_vals, ct[c].v_mean) ;
  ct[c].npoints-- ;
  if (ct[c].npoints == 0)
    ErrorPrintf(ERROR_BADPARM, "%s: empty cluster %d (vertex %d removed)",
                Progname, c, vno) ;
  else
    VectorScalarMul(ct[c].v_mean, 1.0/(double)ct[c].npoints, ct[c].v_mean) ;

  return(NO_ERROR) ;
}
static int
add_vertex_to_cluster(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int c, int vno) {
  int    nsamples;
  static VECTOR *v_vals = NULL ;

  nsamples = mri_profiles->nframes ;

  if (ct[c].npoints == 0)
    ct[c].vno = vno ;
  mris->vertices[vno].curv = c ;
  CTABannotationAtIndex(mris->ct, c, &mris->vertices[vno].annotation);
  if (v_vals == NULL)
    v_vals = VectorAlloc(nsamples, MATRIX_REAL) ;

  VectorScalarMul(ct[c].v_mean, ct[c].npoints, ct[c].v_mean) ;
  load_vals(mri_profiles, v_vals, vno) ;
  VectorAdd(ct[c].v_mean, v_vals, ct[c].v_mean) ;
  ct[c].npoints++ ;
  VectorScalarMul(ct[c].v_mean, 1.0/(double)ct[c].npoints, ct[c].v_mean) ;

  return(NO_ERROR) ;
}

static int
initialize_cluster_centers_with_ico(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, MRI_SURFACE *mris_ico) {
  int             i, j, vno, nsamples, vnos[MAX_CLUSTERS], k ;
  double          r1, r2, res ;
  float           fmin ;
  MRIS_HASH_TABLE *mht ;
  VERTEX          *vico ;

  k = mris_ico->nvertices ;
  MRISstoreRipFlags(mris) ;
  MRISunrip(mris) ;
  r1 = MRISaverageRadius(mris) ;
  r2 = MRISaverageRadius(mris_ico) ;
  MRISscaleBrain(mris_ico,mris_ico, r1/r2);

  res = sqrt(mris->total_area/mris->nvertices) ;
  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 2*res) ;
  nsamples = mri_profiles->nframes ;
  for (i = 0 ; i < mris_ico->nvertices ; i++) {
    vico = &mris_ico->vertices[i] ;
    vno = MRISfindClosestVertex(mris, vico->x, vico->y, vico->z, &fmin) ;
    if (vno < 0)
      continue ;
    vnos[i] = vno ;
    for (j = 0 ; j < nsamples ; j++)
      VECTOR_ELT(ct[i].v_mean, j+1) = MRIgetVoxVal(mri_profiles, vno, 0, 0, j) ;
  }
  mris->ct = CTABalloc(k) ;
  for (i = 0 ; i < k ; i++) {
    mris->vertices[vnos[i]].curv = i ;
    ct[i].npoints++ ;
    ct[i].vno = vnos[i] ;
    CTABannotationAtIndex(mris->ct, i, &mris->vertices[vnos[i]].annotation) ;
  }
  MRISrestoreRipFlags(mris) ;
  return(NO_ERROR) ;
}
static int
initialize_cluster_centers(MRI_SURFACE *mris, MRI *mri_profiles, CLUSTER *ct, int k) {
  int i, j, done, vno, nsamples, vnos[MAX_CLUSTERS], iter ;
  double dist, min_dist ;
  VERTEX *v, *vn ;

  nsamples = mri_profiles->nframes ;
  min_dist = sqrt(mris->total_area/k) / 2 ;
  for (i = 0 ; i < k ; i++) {
    iter = 0 ;
    do {
      do {
        vno = (int)randomNumber(0, mris->nvertices-1) ;
        v = &mris->vertices[vno] ;
      } while (v->ripflag) ;

      // see if any of the other centers are too close to this one
      done = 1 ;
      for (j = 0 ; done && j < i ; j++) {
        if (j == i)
          continue ;
        vn = &mris->vertices[vnos[j]] ;
        dist = sqrt(SQR(vn->x-v->x)+SQR(vn->y-v->y)+SQR(vn->z-v->z)) ;
        if (dist < min_dist) {
          done = 0 ;
          break ;
        }
      }
      if (iter++ > mris->nvertices)
        done = 1 ;
    } while (!done) ;
    vnos[i] = vno ;
    for (j = 0 ; j < nsamples ; j++)
      VECTOR_ELT(ct[i].v_mean, j+1) = MRIgetVoxVal(mri_profiles, vno, 0, 0, j) ;
  }
  mris->ct = CTABalloc(k) ;
  for (i = 0 ; i < k ; i++) {
    mris->vertices[vnos[i]].curv = i ;
    ct[i].npoints++ ;
    ct[i].vno = vnos[i] ;
    CTABannotationAtIndex(mris->ct, i, &mris->vertices[vnos[i]].annotation) ;
  }
  return(NO_ERROR) ;
}
static int
rip_bad_vertices(MRI_SURFACE *mris, MRI *mri_profiles) {
  int    unknown_index, bad, i, vno, index ;
  VERTEX *v ;

  CTABfindName(mris->ct, "Unknown", &unknown_index) ;
  printf("unknown index = %d\n", unknown_index) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    CTABfindAnnotation(mris->ct, v->annotation, &index);
    if (index == unknown_index) {
      v->annotation = 0 ;
      v->ripflag = 1 ;
    } else {
      bad = 1 ;
      for (i = 0 ; i < mri_profiles->nframes ; i++) {
        if (!FZERO(MRIgetVoxVal(mri_profiles, vno, 0, 0, i)))
          bad = 0 ;
      }
      if (bad) {
        v->ripflag = 1 ;
        v->annotation = 0 ;
      }
    }

  }
  MRISsetRipInFacesWithRippedVertices(mris) ;
  MRIScomputeMetricProperties(mris) ;
  return(NO_ERROR) ;
}
static int
okay_to_swap(MRI_SURFACE *mris, int vno, int c1, int c2) {
  int    n, nc1, nc2, n2, okay ;
  VERTEX *v, *vn, *vn2 ;
  float  old_curv ;

  okay = 0 ;
  v = &mris->vertices[vno] ;
  for (nc1 = nc2 = n = 0 ; n < v->vnum ; n++) {
    vn = &mris->vertices[v->v[n]] ;
    if (vn->curv == c1)
      nc1++ ;
    else if (vn->curv == c2)
      nc2++ ;
  }
  if (nc2 > 1) // could be okay, check connectivity.
  {
    // check to make sure the c1 vertices are still connected
    old_curv = v->curv ;
    v->curv = -1 ;  // don't use this for connectivity of c1
    okay = 1 ;
    for (n = 0 ; okay && n < v->vnum ; n++) {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->curv == c1) {
        for (n2 = 0 ; n2 < n ; n2++) {
          vn2 = &mris->vertices[v->v[n2]] ;
          if (vn2->curv != c1)
            continue ;
          if (clusters_connected(mris, v->v[n], v->v[n2], 3) == 0) {
            okay = 0 ;
            break ;
          }
        }
      }
    }
    v->curv = old_curv ;
  }
  return(okay) ;
}

static int
clusters_connected(MRI_SURFACE *mris, int vno1, int vno2, int num) {
  int    n, n2, n3, connected ;
  VERTEX *v, *vn, *vn2, *vn3 ;

  v = &mris->vertices[vno1] ;
  for (connected = n = 0 ; !connected && n < v->vnum ; n++) {
    vn = &mris->vertices[v->v[n]] ;
    if (v->v[n] == vno2)
      connected = 1 ;
    else if (vn->curv == v->curv)  // in same cluster - search its nbhd
    {
      for (n2 = 0 ; !connected && n2 < vn->vnum ; n2++) {
        vn2 = &mris->vertices[vn->v[n2]] ;
        if (vn->v[n2] == vno2)
          connected = 1 ;
        else if (vn2->curv == v->curv)  // in same cluster - search it's nbhd
        {
          for (n3 = 0 ; !connected && n3 < vn2->vnum ; n3++) {
            vn3 = &mris->vertices[vn2->v[n3]] ;
            if (vn2->v[n3] == vno2)
              connected = 1 ;
          }
        }
      }
    }
  }
  return(connected) ;
}

#if 0
int
MRISfindClosestVertex(MRI_SURFACE *mris1, MRI_SURFACE *mris2, int vno1) {
  int    vno, min_vno ;
  double dist, min_dist ;
  VERTEX *v1, *v2 ;

  v1 = &mris1->vertices[vno1] ;
  v2 = &mris2->vertices[0] ;
  min_dist = SQR(v1->x-v2->x)+SQR(v1->y-v2->y)+SQR(v1->z-v2->z);
  min_vno = 0 ;
  for (vno = 1 ; vno < mris2->nvertices ; vno++) {
    v2 = &mris2->vertices[vno] ;
    dist = SQR(v1->x-v2->x)+SQR(v1->y-v2->y)+SQR(v1->z-v2->z);
    if (dist < min_dist) {
      min_dist = dist ;
      min_vno = vno ;
    }
  }
  return(min_vno) ;
}
#endif
