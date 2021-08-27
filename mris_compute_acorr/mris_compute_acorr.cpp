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

#include "mrisurf.h"
#include "mrishash_internals.h"

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "macros.h"
#include "fio.h"
#include "sig.h"
#include "version.h"



/*-------------------------------- CONSTANTS -----------------------------*/

#define BIN_SIZE  1
#define MAX_DIST  20

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
#if 1
static double *MRIScomputeCurvatureAutocorrelation(MRI_SURFACE *mris,
    float *curv, double *acorr,
    double *counts,
    int nbins, float bin_size) ;
#else
static double *MRIScomputeCurvatureAutocorrelation(MRI_SURFACE *mris,
    float bin_size,
    float max_dist, int *pn) ;
#endif



static int   cvector_normalize(float *v, float norm, int num) ;
static int   cvector_accumulate(float *v, float *vtotal, int num) ;
static int   cvector_accumulate_square(float *v, float *vtotal, int num) ;
static int   cvector_subtract(float *v1, float *v2, float *vdst, int num) ;
static int   cvector_compute_variance(float *var,float *mean,int norm,int num);
static int   cvector_compute_t_test(float *c1_mean, float *c1_var,
                                    float *c2_mean, float *c2_var,
                                    int num_class1, int num_class2,
                                    float *pvals, int num) ;
#if 0
static double    cvector_len(float *v, int num) ;
#endif
static double cvector_compute_snr(float *c1_mean, float *c2_mean,
                                  float *verror, float *snr, int num, int *pi);

static int   cvector_mark_low_prob_vertices(float *pvals, float pthresh,
    MRI_SURFACE *mris) ;
static float *cvector_alloc(int num) ;
static int   cvector_clear(float *v, int num) ;
static int   cvector_add_variances(float *c1_var, float *c2_var,
                                   int num_class1, int num_class2,
                                   float *total_var, int nvertices) ;
static int   cvector_track_best_snr(float *vsnr, float *vbest_snr,
                                    float *vbest_avgs, int avgs, int num) ;


static int   write_acorr(const char *fname, double *acorr, double *counts, int n,
                         float bin_size) ;

static int   fill_acorr_holes(double *acorr, double *counts, int nbins) ;

/*-------------------------------- DATA ----------------------------*/

const char *Progname ;

static char *output_subject = NULL ;
static int navgs = 0 ;
static const char *noise_name = "noise_acorr.dat" ;
static const char *signal_name = "signal_acorr.dat" ;
static char *label_name ;
static double pthresh = 0.0 ;

static int max_avgs = 100 ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  MRI_SURFACE  *mris ;
  char         **av, *curv_name, *surf_name, *hemi, fname[STRLEN],
  *cp, *subject_name, subjects_dir[STRLEN],
  **c1_subjects, **c2_subjects ;
  int          ac, nargs, n, num_class1, num_class2, i, nvertices, nbins;
  float        **c1_curvs, **c2_curvs, *curvs, *total_mean, *c1_mean, *c2_mean,
  *class_mean, *c1_var, *c2_var, *class_var, *pvals ;
  double       *noise_acorr, *signal_acorr, *noise_counts, *signal_counts ;
  MRI_SP       *mrisp ;
  LABEL        *area ;
  FILE         *fp ;

  nargs = handleVersionOption(argc, argv, "mris_compute_acorr");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  fp = fopen("scalespace.dat", "w") ;

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

  /* subject_name hemi surface curvature */
  if (argc < 7)
    usage_exit() ;
  if (output_subject == NULL)
    ErrorExit(ERROR_BADPARM,
              "output subject must be specified with -o <subject name>");

  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
              Progname) ;

  strcpy(subjects_dir, cp) ;

  hemi = argv[1] ;
  surf_name = argv[2] ;
  curv_name = argv[3] ;

#define ARGV_OFFSET 4

  /* first determine the number of subjects in each class */
  num_class1 = 0 ;
  n = ARGV_OFFSET ;
  do {
    num_class1++ ;
    n++ ;
    if (argv[n] == NULL || n >= argc)
      ErrorExit(ERROR_BADPARM, "%s: must spectify ':' between class lists",
                Progname) ;
  } while (argv[n][0] != ':') ;

  /* find  # of vertices in output subject surface */
  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s",
		     subjects_dir,output_subject,hemi,surf_name);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  nvertices = mris->nvertices ;
  MRISfree(&mris) ;

  total_mean = (float *)calloc(nvertices, sizeof(float)) ;
  if (!total_mean)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate mean list of %d curvatures",
              Progname, n, nvertices) ;
  c1_mean = (float *)calloc(nvertices, sizeof(float)) ;
  if (!c1_mean)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate c1 mean list of %d curvatures",
              Progname, n, nvertices) ;
  pvals = (float *)calloc(nvertices, sizeof(float)) ;
  if (!pvals)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate pvals",
              Progname, n, nvertices) ;
  c2_mean = (float *)calloc(nvertices, sizeof(float)) ;
  if (!c2_mean)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate c2 mean list of %d curvatures",
              Progname, n, nvertices) ;

  c1_var = (float *)calloc(nvertices, sizeof(float)) ;
  if (!c1_var)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate c1 var list of %d curvatures",
              Progname, n, nvertices) ;
  c2_var = (float *)calloc(nvertices, sizeof(float)) ;
  if (!c2_var)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate c2 var list of %d curvatures",
              Progname, n, nvertices) ;

  num_class2 = 0 ;
  n++ ; /* skip ':' */
  if (n >= argc)
    ErrorExit(ERROR_BADPARM, "%s: class2 list empty", Progname) ;
  do {
    num_class2++ ;
    n++ ;
    if (n >= argc)
      break ;
  } while (argv[n] != NULL) ;

  fprintf(stderr, "%d subjects in class1, %d subjects in class2\n",
          num_class1, num_class2) ;

  c1_subjects = (char **)calloc(num_class1, sizeof(char *)) ;
  c1_curvs = (float **)calloc(num_class1, sizeof(char *)) ;
  c2_subjects = (char **)calloc(num_class2, sizeof(char *)) ;
  c2_curvs = (float **)calloc(num_class2, sizeof(char *)) ;
  for (n = 0 ; n < num_class1 ; n++) {
    c1_subjects[n] = argv[ARGV_OFFSET+n] ;
    c1_curvs[n] = (float *)calloc(nvertices, sizeof(float)) ;
    if (!c1_curvs[n])
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate %dth list of %d curvatures",
                Progname, n, nvertices) ;

    strcpy(c1_subjects[n], argv[ARGV_OFFSET+n]) ;
    /*    fprintf(stderr, "class1[%d] - %s\n", n, c1_subjects[n]) ;*/
  }
  i = n+1+ARGV_OFFSET ;  /* starting index */
  for (n = 0 ; n < num_class2 ; n++) {
    c2_subjects[n] = argv[i+n] ;
    c2_curvs[n] = (float *)calloc(nvertices, sizeof(float)) ;
    if (!c2_curvs[n])
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate %dth list of %d curvatures",
                Progname, n, nvertices) ;
    strcpy(c2_subjects[n], argv[i+n]) ;
    /*    fprintf(stderr, "class2[%d] - %s\n", n, c2_subjects[n]) ;*/
  }

  if (label_name) {
    area = LabelRead(output_subject, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname,
                label_name) ;
  } else
    area = NULL ;

  /* real all the curvatures in for group1 */
  for (n = 0 ; n < num_class1+num_class2 ; n++) {
    /* transform each subject's curvature into the output subject's space */
    subject_name = n < num_class1 ? c1_subjects[n] : c2_subjects[n-num_class1];
    fprintf(stderr, "reading subject %d of %d: %s\n",
            n+1, num_class1+num_class2, subject_name) ;
    int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s",
		       subjects_dir,subject_name,hemi,surf_name);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    if (strchr(curv_name, '/') != NULL) {
      strcpy(fname, curv_name) ;  /* full path specified */
    } else {
      int req = snprintf(fname,STRLEN,"%s/%s/surf/%s.%s",
			 subjects_dir,subject_name,hemi,curv_name);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: could no read curvature file %s",Progname,fname) ;
    MRISaverageCurvatures(mris, navgs) ;
    mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
    MRISfree(&mris) ;

    req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s",
		   subjects_dir,output_subject,hemi,surf_name);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISfromParameterization(mrisp, mris, 0) ;
    if (label_name)
      MRISmaskNotLabel(mris, area) ;
    curvs = (n < num_class1) ? c1_curvs[n] : c2_curvs[n-num_class1] ;
    class_mean = (n < num_class1) ? c1_mean : c2_mean ;
    class_var = (n < num_class1) ? c1_var : c2_var ;
    MRISextractCurvatureVector(mris, curvs) ;
    cvector_accumulate(curvs, total_mean, nvertices) ;
    cvector_accumulate(curvs, class_mean, nvertices) ;
    cvector_accumulate_square(curvs, class_var, nvertices) ;
    MRISPfree(&mrisp) ;
    MRISfree(&mris) ;
  }

  /* compute within-group means, and total mean */
  cvector_normalize(total_mean, num_class1+num_class2, nvertices) ;
  cvector_normalize(c1_mean, num_class1, nvertices) ;
  cvector_normalize(c2_mean, num_class2, nvertices) ;
  cvector_compute_variance(c1_var, c1_mean, num_class1, nvertices) ;
  cvector_compute_variance(c2_var, c2_mean, num_class2, nvertices) ;
  cvector_compute_t_test(c1_mean, c1_var, c2_mean, c2_var,
                         num_class1, num_class2, pvals, nvertices) ;

  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s",
		 subjects_dir,output_subject,hemi,surf_name);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "reading output surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  if (max_avgs) {
    float  *vbest_snr, *vbest_avgs, *vtotal_var, *vsnr ;
    int    avgs, max_snr_avgs, i ;
    double snr, max_snr ;

    vbest_snr = cvector_alloc(nvertices) ;
    vbest_avgs = cvector_alloc(nvertices) ;
    vtotal_var = cvector_alloc(nvertices) ;
    vsnr = cvector_alloc(nvertices) ;

    cvector_add_variances(c1_var, c2_var, num_class1, num_class2,
                          vtotal_var, nvertices) ;
    snr = cvector_compute_snr(c1_mean, c2_mean, vtotal_var,
                              vsnr, nvertices, &i);
    fprintf(stderr, "raw SNR=%2.2f\n", snr/(double)nvertices) ;
    max_snr = snr ;
    max_snr_avgs = 0 ;
    cvector_track_best_snr(vsnr, vbest_snr, vbest_avgs, 0, nvertices) ;

    for (avgs = 1 ; avgs <= max_avgs ; avgs++) {
      fprintf(stderr, "testing %d averages...\n", avgs) ;
      cvector_clear(c1_mean, nvertices) ;
      cvector_clear(c2_mean, nvertices) ;
      cvector_clear(c1_var, nvertices) ;
      cvector_clear(c1_var, nvertices) ;
      cvector_clear(total_mean, nvertices) ;
      for (n = 0 ; n < num_class1 ; n++) {
#if 0
        fprintf(stderr, "processing subject %d of %d: %s\r",
                n+1, num_class1+num_class2, c1_subjects[n]) ;
#endif
        MRISimportCurvatureVector(mris, c1_curvs[n]) ;
        MRISaverageCurvatures(mris, 1) ;
        MRISextractCurvatureVector(mris, c1_curvs[n]) ;
        cvector_accumulate(c1_curvs[n], total_mean, nvertices) ;
        cvector_accumulate(c1_curvs[n], c1_mean, nvertices) ;
        cvector_accumulate_square(c1_curvs[n], c1_var, nvertices) ;
      }
      for (n = 0 ; n < num_class2 ; n++) {
#if 0
        fprintf(stderr, "processing subject %d of %d: %s\r",
                n+1+num_class1, num_class1+num_class2, c2_subjects[n]) ;
#endif
        MRISimportCurvatureVector(mris, c2_curvs[n]) ;
        MRISaverageCurvatures(mris, 1) ;
        MRISextractCurvatureVector(mris, c2_curvs[n]) ;
        cvector_accumulate(c2_curvs[n], total_mean, nvertices) ;
        cvector_accumulate(c2_curvs[n], c2_mean, nvertices) ;
        cvector_accumulate_square(c2_curvs[n], c2_var, nvertices) ;
      }
      cvector_normalize(total_mean, num_class1+num_class2, nvertices) ;
      cvector_normalize(c1_mean, num_class1, nvertices) ;
      cvector_normalize(c2_mean, num_class2, nvertices) ;
      cvector_compute_variance(c1_var, c1_mean, num_class1, nvertices) ;
      cvector_compute_variance(c2_var, c2_mean, num_class2, nvertices) ;
      cvector_add_variances(c1_var, c2_var, num_class1, num_class2,
                            vtotal_var, nvertices) ;
      snr =
        cvector_compute_snr(c1_mean, c2_mean, vtotal_var, vsnr, nvertices,&i);
      fprintf(fp, "%d %2.1f  %2.2f %2.2f %2.2f ",
              avgs, sqrt((float)avgs), sqrt(snr), c1_mean[i]-c2_mean[i],
              sqrt(vtotal_var[i])) ;
      fflush(fp) ;
      for (n = 0 ; n < num_class1 ; n++)
        fprintf(fp, "%2.2f ", c1_curvs[n][i]) ;
      for (n = 0 ; n < num_class2 ; n++)
        fprintf(fp, "%2.2f ", c2_curvs[n][i]) ;
      fprintf(fp, "\n") ;
      if (snr > max_snr) {
        fprintf(stderr,
                "new max SNR found at avgs=%d (%2.1f mm)=%2.2f, n=%2.2f, "
                "d=%2.2f\n",
                avgs, sqrt((float)avgs), sqrt(snr), c1_mean[i]-c2_mean[i],
                sqrt(vtotal_var[i])) ;
        max_snr = snr ;
        max_snr_avgs = avgs ;
      }
      cvector_track_best_snr(vsnr, vbest_snr, vbest_avgs, avgs, nvertices) ;
    }
    printf("max snr=%2.2f at %d averages\n", max_snr, max_snr_avgs) ;
    MRISimportValVector(mris, vbest_snr) ;
    sprintf(fname, "./%s.best_snr", hemi) ;
    MRISwriteValues(mris, fname) ;
    MRISimportValVector(mris, vbest_avgs) ;
    sprintf(fname, "./%s.best_avgs", hemi) ;
    MRISwriteValues(mris, fname) ;
    fclose(fp) ;
    exit(0) ;
  }
  fclose(fp) ;

  nbins = MAX_DIST / BIN_SIZE ;

  MRISclearMarks(mris) ;
  if (pthresh > 0)
    cvector_mark_low_prob_vertices(pvals, pthresh, mris) ;

  noise_acorr = (double *)calloc(nbins, sizeof(double)) ;
  noise_counts = (double *)calloc(nbins, sizeof(double)) ;
  signal_acorr = (double *)calloc(nbins, sizeof(double)) ;
  signal_counts = (double *)calloc(nbins, sizeof(double)) ;
  curvs = cvector_alloc(nvertices) ;
  cvector_subtract(c1_mean, c2_mean, curvs, nvertices) ;
  MRISimportCurvatureVector(mris, curvs) ;
  sprintf(fname, "./%s.c1_c2", hemi) ;
  MRISwriteCurvature(mris, fname) ;

  sprintf(fname, "./%s.c1", hemi) ;
  MRISimportCurvatureVector(mris, c1_mean) ;
  MRISwriteCurvature(mris, fname) ;

  sprintf(fname, "./%s.c2", hemi) ;
  MRISimportCurvatureVector(mris, c2_mean) ;
  MRISwriteCurvature(mris, fname) ;

  fprintf(stderr, "computing autocorrelation functions...\n") ;
  for (n = 0 ; n < num_class1 ; n++) {
    fprintf(stderr, "processing subject %d of %d: %s\n",
            n+1, num_class1+num_class2, c1_subjects[n]) ;
    cvector_subtract(c1_mean, c1_curvs[n], curvs, nvertices) ;
    MRIScomputeCurvatureAutocorrelation(mris,curvs,noise_acorr,noise_counts,
                                        nbins,BIN_SIZE);
#if 0
    cvector_subtract(c2_mean, c1_curvs[n], curvs, nvertices) ;
    MRIScomputeCurvatureAutocorrelation(mris,curvs,signal_acorr,signal_counts,
                                        nbins,BIN_SIZE);
#endif
  }
  for (n = 0 ; n < num_class2 ; n++) {
    fprintf(stderr, "processing subject %d of %d: %s\n",
            n+1+num_class1, num_class1+num_class2, c2_subjects[n]) ;
    cvector_subtract(c2_mean, c2_curvs[n], curvs, nvertices) ;
    MRIScomputeCurvatureAutocorrelation(mris,curvs,noise_acorr,noise_counts,
                                        nbins,BIN_SIZE);
#if 0
    cvector_subtract(c1_mean, c2_curvs[n], curvs, nvertices) ;
    MRIScomputeCurvatureAutocorrelation(mris, curvs,signal_acorr,signal_counts,
                                        nbins,BIN_SIZE);
#endif
  }
  cvector_subtract(c1_mean, c2_mean, curvs, nvertices) ;
  MRIScomputeCurvatureAutocorrelation(mris, curvs,signal_acorr,signal_counts,
                                      nbins,BIN_SIZE);
  for (n = 0 ; n < nbins ; n++) {
    if (signal_counts[n])
      signal_acorr[n] /= signal_counts[n] ;
    if (noise_counts[n])
      noise_acorr[n] /= noise_counts[n] ;
  }
  fill_acorr_holes(signal_acorr, signal_counts, nbins) ;
  fill_acorr_holes(noise_acorr, noise_counts, nbins) ;
  write_acorr(signal_name, signal_acorr, signal_counts, nbins, BIN_SIZE) ;
  write_acorr(noise_name, noise_acorr, noise_counts, nbins, BIN_SIZE) ;
#if 0
  exit(0) ;

  acorr = MRIScomputeCurvatureAutocorrelation(mris, BIN_SIZE, MAX_DIST, &n) ;
  free(acorr) ;
  MRISfree(&mris) ;
#endif

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
  else if (!stricmp(option, "avgs")) {
    max_avgs = atoi(argv[2]) ;
    fprintf(stderr,
            "computing kernel for maximum snr up to %d averages (%2.1f mm)\n",
            max_avgs, sqrt((double)max_avgs)) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'L':
      label_name = argv[2] ;
      fprintf(stderr, "masking label %s\n", label_name) ;
      nargs = 1 ;
      break ;
    case 'N':
      noise_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, "outputting noise autocorrelation to %s...\n", noise_name);
      break ;
    case 'S':
      signal_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr,"outputting signal autocorrelation to %s...\n",signal_name);
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "averaging curvature patterns %d times...\n", navgs) ;
      break ;
    case 'T':
      pthresh = atof(argv[2]) ;
      fprintf(stderr, "using p value threshold %2.3f\n", (float)pthresh) ;
      nargs = 1 ;
      break ;
    case 'O':
      output_subject = argv[2] ;
      fprintf(stderr, "using %s as output subject\n", output_subject) ;
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
          "usage: %s -o <output subject> [options] \n"
          "\t<hemi> <surf> <curv> \n\t<c1_subject1> <c1_subject2>... : \n"
          "\t<c2_subject1> <c2_subject2>...\n",
          Progname) ;
  fprintf(stderr, "where surf must be a spherical surface suitable for "
          "computing geodesics\n") ;
  fprintf(stderr, "The <c1_subject> ... is a list of subjects from one class\n"
          "and the <c2_subject>... is a list of subjects from another "
          "class.\n");
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will compute the autocorrelation function of"
          " a curvature file\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

#if 1
static double *
MRIScomputeCurvatureAutocorrelation(MRI_SURFACE *mris,float*curv,double *acorr,
                                    double *counts, int nbins, float bin_size) {
  static MHT     *mht = NULL ;
  static int nv = 0 ;
  int     vno, n, i, index ;
  VERTEX  *v, *vn ;
  float   x, y, z, radius, dist, max_dist ;
  double  angle, circumference ;
  VECTOR  *v1, *v2 ;

  max_dist = nbins * bin_size ;
  if (nv != mris->nvertices && mht != NULL)
    MHTfree(&mht) ;

  nv = mris->nvertices ;
  if (!mht) {
    fprintf(stderr, "building spatial LUT...\n") ;
    mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 2*max_dist) ;
  }
  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;
  radius = MRISaverageRadius(mris) ;
  circumference = M_PI * 2.0 * radius ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || !v->marked)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
#if 0
    if (!(vno % 10000))
      fprintf(stderr, "%d of %d vertices processed\n", vno, mris->nvertices) ;
#endif
    x = v->x ;
    y = v->y ;
    z = v->z ;
    MHBT* bucket = MHTacqBucket(mht, x, y, z) ;
    VECTOR_LOAD(v1, v->x, v->y, v->z) ;  /* radius vector */
    MHB* bin ;
    for (bin = bucket->bins, i = 0 ; i < bucket->nused ; i++, bin++) {
      n = bin->fno ;
      vn = &mris->vertices[n] ;
      VECTOR_LOAD(v2, vn->x, vn->y, vn->z) ;  /* radius vector */
      angle = fabs(Vector3Angle(v1, v2)) ;
#if 0
      xd = v->x - vn->x ;
      yd = v->y - vn->y ;
      zd = v->z - vn->z ;
      dist = sqrt(xd*xd + yd*yd + zd*zd) ;
#else
      dist = circumference * angle / (2.0 * M_PI) ;
#endif
      if (dist < max_dist) {
        index = (int)((float)dist/bin_size) ;
        counts[index]++ ;
#if 0
        if (!index) {
          double a = (double)(curv[vno] * curv[n]) ;
          DiagBreak() ;
        }
#endif
        acorr[index] += (double)(curv[vno] * curv[n]) ;
      }
    }
    MHTrelBucket(&bucket);
  }

  VectorFree(&v1) ;
  VectorFree(&v2) ;
  return(acorr) ;
}
#else
static double *
MRIScomputeCurvatureAutocorrelation(MRI_SURFACE *mris, float bin_size,
                                    float max_dist, int *pn) {
  MHT     *mht ;
  int     vno, n, i, index, *counts, nbins ;
  VERTEX  *v, *vn ;
  float   x, y, z, radius, dist ;
  double  angle, circumference, *acorr ;
  VECTOR  *v1, *v2 ;

  fprintf(stderr, "building spatial LUT...\n") ;
  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;
  radius = MRISaverageRadius(mris) ;
  circumference = M_PI * 2.0 * radius ;
  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 2*max_dist) ;

  nbins = max_dist/bin_size+1 ;
  counts = (int *)calloc(nbins, sizeof(int)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (!(vno % 10000))
      fprintf(stderr, "%d of %d vertices processed\n", vno, mris->nvertices) ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    MHBT* bucket = MHTacqBucket(mht, x, y, z) ;
    VECTOR_LOAD(v1, v->x, v->y, v->z) ;  /* radius vector */
    MHB* bin ;
    for (bin = bucket->bins, i = 0 ; i < bucket->nused ; i++, bin++) {
      n = bin->fno ;
      vn = &mris->vertices[n] ;
      VECTOR_LOAD(v2, vn->x, vn->y, vn->z) ;  /* radius vector */
      angle = fabs(Vector3Angle(v1, v2)) ;
      dist = circumference * angle / (2.0 * M_PI) ;
      if (dist < max_dist) {
        index = (int)((float)dist/bin_size) ;
        counts[index]++ ;
        acorr[index] += (double)(v->curv * vn->curv) ;
      }
    }
    MHTrelBucket(&bucket);
  }

  MHTfree(&mht) ;

  for (i = 0 ; i < nbins ; i++) {
    if (counts[i]) {
      acorr[i] /= (float)counts[i] ;
      printf("%2.4f  %2.4f  %d\n",
             (float)i*bin_size, acorr[i], counts[i]) ;
    }
#if 0
    else
      printf("0  0  0\n") ;
#endif
  }
  *pn = nbins ;
  free(counts) ;
  VectorFree(&v1) ;
  VectorFree(&v2) ;
  return(acorr) ;
}
#endif

static int
cvector_compute_variance(float *var,float *mean,int norm,int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    var[i] = var[i] / (float)norm - mean[i]*mean[i] ;
  return(NO_ERROR) ;
}
static int
cvector_normalize(float *v, float norm, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    v[i] /= norm ;
  return(NO_ERROR) ;
}

static int
cvector_accumulate_square(float *v, float *vtotal, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    vtotal[i] += v[i]*v[i] ;
  return(NO_ERROR) ;
}
static int
cvector_accumulate(float *v, float *vtotal, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    vtotal[i] += v[i] ;
  return(NO_ERROR) ;
}

static int
cvector_subtract(float *v1, float *v2, float *vdst, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    vdst[i] = v1[i] - v2[i] ;
  return(NO_ERROR) ;
}
static int
cvector_compute_t_test(float *c1_mean, float *c1_var,
                       float *c2_mean, float *c2_var,
                       int num_class1, int num_class2,
                       float *pvals, int num) {
  int    i ;
  double t, numer, denom ;

  for (i = 0 ; i < num ; i++) {
    numer = (c1_mean[i] - c2_mean[i]) ;
    denom = sqrt(c1_var[i]/num_class1) + sqrt(c2_var[i]/num_class2) ;
    t = numer / denom ;
    pvals[i] = sigt(t, num_class1+num_class2-2) ;
  }
  return(NO_ERROR) ;
}

static int
write_acorr(const char *fname,double *acorr,double *counts,int nbins, float bin_size) {
  int  i ;
  FILE *fp ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE, "%s: could not open autocorrelation file %s",
                 Progname, fname)) ;

  for (i = 0 ; i < nbins ; i++) {
    /*    if (counts[i])*/
    {
      fprintf(fp, "%2.4f  %2.4f  %d\n",
              (float)i*bin_size, acorr[i], (int)counts[i]) ;
    }
#if 0
    else
      printf("0  0  0\n") ;
#endif
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

static int
fill_acorr_holes(double *acorr, double *counts, int nbins) {
  int    i, n ;

  for (n = 0 ; n < 5000  ; n++) {
    for (i = 0 ; i < nbins ; i++) {
      if (counts[i] == 0) {
        if (i == 0)
          acorr[i] = (acorr[i] + acorr[i+1]) / 2 ;
        else if (i == nbins-1)
          acorr[i] = (acorr[i] + acorr[i-1]) / 2 ;
        else
          acorr[i] = (acorr[i] + acorr[i-1] + acorr[i+1]) / 3 ;
      }
    }
  }
  return(NO_ERROR) ;
}


static int
cvector_mark_low_prob_vertices(float *pvals, float pthresh, MRI_SURFACE *mris) {
  int    i, num ;

  for (num = i = 0 ; i < mris->nvertices ; i++) {
    if (pvals[i] < pthresh) {
      num++ ;
      mris->vertices[i].marked = 1 ;
    }
  }
  fprintf(stderr, "%d vertices p < %2.2e\n", num, pthresh) ;
  return(NO_ERROR) ;
}

static double
cvector_compute_snr(float *c1_mean, float *c2_mean, float *vvar, float *snr,
                    int num, int *pi) {
  int    i, max_i ;
  float  f, max_snr ;
  double total_snr ;

  max_i = -1 ;
  for (max_snr = total_snr = 0.0, i = 0 ; i < num ; i++) {
    f = (c1_mean[i]-c2_mean[i]) ;
    f *= f ;
    if (!iszero(vvar[i]))
      f /= (vvar[i]) ;

    if (f > max_snr) {
      max_snr = f ;
      max_i = i ;
    }
    snr[i] = f ;
    total_snr += snr[i] ;
  }
  total_snr /= (double)num ;
  *pi = max_i ;
  return(max_snr) ;
}
#if 0
static double
cvector_len(float *v, int num) {
  int    i ;
  double len ;

  for (len = 0.0, i = 0 ; i < num ; i++)
    len += v[i]*v[i] ;
  len /= (double)num ;
  len = sqrt(len) ;
  return(len) ;
}
#endif
static float *
cvector_alloc(int num) {
  float *v ;

  v = (float *)calloc(num, sizeof(float)) ;
  if (!v)
    ErrorExit(ERROR_NOMEMORY, "cvector_alloc(%d): calloc failed", num) ;
  return(v) ;
}

static int
cvector_clear(float *v, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    v[i] = 0 ;
  return(NO_ERROR) ;
}
static int
cvector_add_variances(float *c1_var, float *c2_var, int num_class1,
                      int num_class2, float *vtotal_var, int nvertices) {
  int     i, total_dof ;

  total_dof = num_class1 + num_class2 ;
  for (i = 0 ; i < nvertices ; i++)
    vtotal_var[i] = (c1_var[i]*num_class1 * c2_var[i]*num_class2) / total_dof ;

  return(NO_ERROR) ;
}


static int
cvector_track_best_snr(float *vsnr, float *vbest_snr, float *vbest_avgs,
                       int avgs, int num) {
  int    i ;

  for (i = 0 ; i < num ; i++) {
    if (vsnr[i] > vbest_snr[i]) {
      vbest_snr[i] = vsnr[i] ;
      vbest_avgs[i] = avgs ;
    }
  }

  return(NO_ERROR) ;
}

