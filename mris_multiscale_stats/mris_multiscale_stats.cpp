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

#include "timer.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"
#include "fio.h"
#include "mrishash.h"
#include "sig.h"
#include "version.h"



/*-------------------------------- CONSTANTS -----------------------------*/

#define MIN_LABELS   5

/*-------------------------------- PROTOTYPES ----------------------------*/

static int  cvector_expand_label(LABEL *area, float f, float *vout) ;
static int  extract_thickness_at_best_scale(MRI_SURFACE *mris,
    float **c1_avg_curvs,
    float *vbest_avgs,
    float **c1_curvs,
    int num, int nvectors);
static int  mark_thresholded_vertices(MRI_SURFACE *mris, float *vbest_snr,
                                      float *vbest_avgs, float thresh) ;
static int  segment_and_write_labels(char *subject, char *fname,
                                     MRI_SURFACE *mris, LABEL ***plabels,
                                     int *pnlabels, int offset) ;
int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;



static int cvector_divide(float *v, float div, int num) ;
static int cvector_compute_t_test(float *c1_mean, float *c1_var,
                                  float *c2_mean, float *c2_var,
                                  int num_class1, int num_class2,
                                  float *pvals, int num) ;
static int   cvector_normalize(float *v, float norm, int num) ;
static int   cvector_accumulate(float *v, float *vtotal, int num) ;
static int   cvector_accumulate_square(float *v, float *vtotal, int num) ;
static int   cvector_compute_variance(float *var,float *mean,int norm,int num);
static int   cvector_bonferroni_correct(float *pvalues, float *avgs, int num) ;
static double cvector_average_in_label(float *v, LABEL *area, int num) ;
static int   cvector_extract_best_avg(float *vbest_avgs, float *vsrc,
                                      float *vdst, int navgs, int num) ;
#if 0
static int  extract_thickness_at_best_scale(MRI_SURFACE *mris,
    float **c1_avg_curvs,
    float *vbest_avgs,
    float **c1_curvs,
    int num, int nvectors);
static double    cvector_len(float *v, int num) ;
#endif
static int    cvector_copy(float *v1, float *v2, int num) ;
static double cvector_compute_dist_free_snr(float **c1_curvs, int num_class1,
    float **c2_curvs, int num_class2,
    float *c1_mean, float *c2_mean,
    float *vsnr, int num, int *pi);
static double cvector_compute_snr(float *c1_mean, float *c2_mean,
                                  float *verror, float *snr, int num, int *pi,
                                  float bonferroni);

#if 0
static int   cvector_subtract(float *v1, float *v2, float *vdst, int num) ;
static int   cvector_mark_low_prob_vertices(float *pvals, float pthresh,
    MRI_SURFACE *mris) ;
static int   cvector_track_best_snr(float *vsnr, float *vbest_snr,
                                    float *vbest_avgs, int avgs, int num) ;

#endif

static float *cvector_alloc(int num) ;
static int   cvector_clear(float *v, int num) ;
static int   cvector_add_variances(float *c1_var, float *c2_var,
                                   int num_class1, int num_class2,
                                   float *total_var, int nvertices) ;
static int   cvector_track_best_stats(float *vsnr, float *vbest_snr,
                                      float *vbest_avgs,
                                      float *c1_var, float *c2_var,
                                      float *c1_best_var, float *c2_best_var,
                                      float *c1_mean, float *c2_mean,
                                      float *c1_best_mean, float *c2_best_mean,
                                      int avgs, int num) ;


/*-------------------------------- DATA ----------------------------*/

static float tthresh = 1.0f ;
static int roi_flag = 0 ;
static int cond_no1 = 0 ;
static int cond_no2 = 1 ;

const char *Progname ;

static float min_label_area = 30.0f ;
static int write_flag = 0 ;
static char *output_subject = NULL ;
static char *label_name = NULL ;
static const char *prefix = "" ;

static int max_avgs = 500 ;
static int use_no_distribution = 0 ;

static int bonferroni = 0 ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  MRI_SURFACE  *mris ;
  char         **av, *curv_name, *surf_name, *hemi, fname[STRLEN],
  *cp, *subject_name, subjects_dir[STRLEN],
  **c1_subjects, **c2_subjects ;
  int          ac, nargs, n, num_class1, num_class2, i, nvertices,
  avgs, max_snr_avgs, nlabels ;
  float        **c1_curvs, **c2_curvs, *curvs, *c1_mean, *c2_mean,
  *class_mean, *c1_var, *c2_var, *class_var,*pvals,**c1_avg_curvs,
  *vbest_snr, *vbest_avgs, *vtotal_var, *vsnr, **c2_avg_curvs,
  *vbest_pvalues, *c1_best_var, *c2_best_var, *c1_best_mean,
  *c2_best_mean ;
  MRI_SP       *mrisp ;
  LABEL        *area, **labels = NULL ;
  FILE         *fp = NULL ;
  double       snr, max_snr ;
  double       **c1_thickness, **c2_thickness ;
  Timer start ;
  int          msec, minutes, seconds ;

  nargs = handleVersionOption(argc, argv, "mris_multiscale_stats");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  if (write_flag && DIAG_VERBOSE_ON)
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

  start.reset() ;

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

  pvals = cvector_alloc(nvertices) ;
  c1_mean = cvector_alloc(nvertices) ;
  c2_mean = cvector_alloc(nvertices) ;
  c1_best_mean = cvector_alloc(nvertices) ;
  c2_best_mean = cvector_alloc(nvertices) ;
  c1_var = cvector_alloc(nvertices) ;
  c2_var = cvector_alloc(nvertices) ;
  c1_best_var = cvector_alloc(nvertices) ;
  c2_best_var = cvector_alloc(nvertices) ;

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

  fprintf(stderr, "%d subjects in class 1, %d subjects in class 2\n",
          num_class1, num_class2) ;

  c1_subjects = (char **)calloc(num_class1, sizeof(char *)) ;
  c1_curvs = (float **)calloc(num_class1, sizeof(char *)) ;
  c1_avg_curvs = (float **)calloc(num_class1, sizeof(char *)) ;
  c2_subjects = (char **)calloc(num_class2, sizeof(char *)) ;
  c2_curvs = (float **)calloc(num_class2, sizeof(char *)) ;
  c2_avg_curvs = (float **)calloc(num_class2, sizeof(char *)) ;
  for (n = 0 ; n < num_class1 ; n++) {
    c1_subjects[n] = argv[ARGV_OFFSET+n] ;
    c1_curvs[n] = cvector_alloc(nvertices) ;
    c1_avg_curvs[n] = cvector_alloc(nvertices) ;
    strcpy(c1_subjects[n], argv[ARGV_OFFSET+n]) ;
  }
  i = n+1+ARGV_OFFSET ;  /* starting index */
  for (n = 0 ; n < num_class2 ; n++) {
    c2_subjects[n] = argv[i+n] ;
    c2_curvs[n] = cvector_alloc(nvertices) ;
    c2_avg_curvs[n] = cvector_alloc(nvertices) ;
    strcpy(c2_subjects[n], argv[i+n]) ;
  }

  if (label_name) {
    area = LabelRead(output_subject, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname,
                label_name) ;
  } else
    area = NULL ;

  /* real all the curvatures in for both groups  */
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
    if (area)
      MRISmaskNotLabel(mris, area) ;
    curvs = (n < num_class1) ? c1_curvs[n] : c2_curvs[n-num_class1] ;
    class_mean = (n < num_class1) ? c1_mean : c2_mean ;
    class_var = (n < num_class1) ? c1_var : c2_var ;
    MRISextractCurvatureVector(mris, curvs) ;
    cvector_accumulate(curvs, class_mean, nvertices) ;
    cvector_accumulate_square(curvs, class_var, nvertices) ;
    MRISPfree(&mrisp) ;
    MRISfree(&mris) ;
  }

  /* compute within-group means, and total mean */
  cvector_normalize(c1_mean, num_class1, nvertices) ;
  cvector_normalize(c2_mean, num_class2, nvertices) ;
  cvector_compute_variance(c1_var, c1_mean, num_class1, nvertices) ;
  cvector_compute_variance(c2_var, c2_mean, num_class2, nvertices) ;

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

  if (area)
    MRISripNotLabel(mris, area) ;
  vbest_snr = cvector_alloc(nvertices) ;
  vbest_pvalues = cvector_alloc(nvertices) ;
  vbest_avgs = cvector_alloc(nvertices) ;
  vtotal_var = cvector_alloc(nvertices) ;
  vsnr = cvector_alloc(nvertices) ;

  cvector_add_variances(c1_var, c2_var, num_class1, num_class2,
                        vtotal_var, nvertices) ;
  if (use_no_distribution)
    snr = cvector_compute_dist_free_snr(c1_curvs, num_class1, c2_curvs,
                                        num_class2, c1_mean, c2_mean,
                                        vsnr, nvertices, &i);
  else
    snr = cvector_compute_snr(c1_mean, c2_mean, vtotal_var, vsnr, nvertices,
                              &i, 0.0f);
  fprintf(stderr,
          "raw SNR %2.2f, n=%2.4f, d=%2.4f, vno=%d\n",
          sqrt(snr), c1_mean[i]-c2_mean[i], sqrt(vtotal_var[i]), i) ;
  max_snr = snr ;
  max_snr_avgs = 0 ;
  cvector_track_best_stats(vsnr, vbest_snr, vbest_avgs,
                           c1_var, c2_var, c1_best_var, c2_best_var,
                           c1_mean, c2_mean, c1_best_mean, c2_best_mean,
                           0, nvertices) ;

  /*
    c?_avg_curvs will contain the curvature at the current scale,
    c?_curvs will contain the curvature at the best scale.
    c?_var will contain the the variance at the current scale.
    c?_best_var will contain the variance at the best scale.
  */
  for (n = 0 ; n < num_class1 ; n++)
    cvector_copy(c1_curvs[n], c1_avg_curvs[n], nvertices) ;
  for (n = 0 ; n < num_class2 ; n++)
    cvector_copy(c2_curvs[n], c2_avg_curvs[n], nvertices) ;
  for (avgs = 1 ; avgs <= max_avgs ; avgs++) {
    if (!(avgs % 50))
      fprintf(stderr, "testing %d averages...\n", avgs) ;
    cvector_clear(c1_mean, nvertices) ;
    cvector_clear(c2_mean, nvertices) ;
    cvector_clear(c1_var, nvertices) ;
    cvector_clear(c1_var, nvertices) ;
    for (n = 0 ; n < num_class1 ; n++) {
#if 0
      fprintf(stderr, "processing subject %d of %d: %s\r",
              n+1, num_class1+num_class2, c1_subjects[n]) ;
#endif
      MRISimportCurvatureVector(mris, c1_avg_curvs[n]) ;
      MRISaverageCurvatures(mris, 1) ;
      MRISextractCurvatureVector(mris, c1_avg_curvs[n]) ;
      cvector_accumulate(c1_avg_curvs[n], c1_mean, nvertices) ;
      cvector_accumulate_square(c1_avg_curvs[n], c1_var, nvertices) ;
    }
    for (n = 0 ; n < num_class2 ; n++) {
#if 0
      fprintf(stderr, "processing subject %d of %d: %s\r",
              n+1+num_class1, num_class1+num_class2, c2_subjects[n]) ;
#endif
      MRISimportCurvatureVector(mris, c2_avg_curvs[n]) ;
      MRISaverageCurvatures(mris, 1) ;
      MRISextractCurvatureVector(mris, c2_avg_curvs[n]) ;
      cvector_accumulate(c2_avg_curvs[n], c2_mean, nvertices) ;
      cvector_accumulate_square(c2_avg_curvs[n], c2_var, nvertices) ;
    }
    cvector_normalize(c1_mean, num_class1, nvertices) ;
    cvector_normalize(c2_mean, num_class2, nvertices) ;
    cvector_compute_variance(c1_var, c1_mean, num_class1, nvertices) ;
    cvector_compute_variance(c2_var, c2_mean, num_class2, nvertices) ;
    cvector_add_variances(c1_var, c2_var, num_class1, num_class2,
                          vtotal_var, nvertices) ;
    if (use_no_distribution)
      snr = cvector_compute_dist_free_snr(c1_avg_curvs, num_class1,
                                          c2_avg_curvs, num_class2, c1_mean,
                                          c2_mean, vsnr, nvertices, &i);
    else
      snr =
        cvector_compute_snr(c1_mean, c2_mean, vtotal_var, vsnr, nvertices,&i,
                            bonferroni ? log((double)avgs) : 0.0f);
    if (write_flag && DIAG_VERBOSE_ON) {
      fprintf(fp, "%d %2.1f  %2.2f %2.2f %2.2f ",
              avgs, sqrt((float)avgs), sqrt(snr), c1_mean[i]-c2_mean[i],
              sqrt(vtotal_var[i])) ;
      fflush(fp) ;
      for (n = 0 ; n < num_class1 ; n++)
        fprintf(fp, "%2.2f ", c1_avg_curvs[n][i]) ;
      for (n = 0 ; n < num_class2 ; n++)
        fprintf(fp, "%2.2f ", c2_avg_curvs[n][i]) ;
      fprintf(fp, "\n") ;
    }
    if (snr > max_snr) {
      fprintf(stderr,
              "new max SNR found at avgs=%d (%2.1f mm)=%2.1f, n=%2.4f, "
              "d=%2.4f, vno=%d\n",
              avgs, sqrt((float)avgs), sqrt(snr), c1_mean[i]-c2_mean[i],
              sqrt(vtotal_var[i]), i) ;
      max_snr = snr ;
      max_snr_avgs = avgs ;
    }
    cvector_track_best_stats(vsnr, vbest_snr, vbest_avgs,
                             c1_var, c2_var, c1_best_var, c2_best_var,
                             c1_mean, c2_mean, c1_best_mean, c2_best_mean,
                             avgs, nvertices) ;
  }
  printf("max snr=%2.2f at %d averages\n", max_snr, max_snr_avgs) ;
  cvector_compute_t_test(c1_best_mean, c1_best_var, c2_best_mean, c2_best_var,
                         num_class1, num_class2, vbest_pvalues, nvertices) ;
  if (bonferroni)
    cvector_bonferroni_correct(vbest_pvalues, vbest_avgs, nvertices) ;

  free(c1_mean) ;
  free(c2_mean) ;
  free(c1_var);
  free(c2_var);
  if (roi_flag)  /* collapse across label */
  {
    do {
      int   npos_labels, nneg_labels ;
      LABEL **pos_labels, **neg_labels ;

      MRISclearMarks(mris) ;
      sprintf(fname, "%s-%s_thickness", hemi, prefix ? prefix : "") ;
      mark_thresholded_vertices(mris, vbest_pvalues, vbest_avgs, tthresh) ;
      segment_and_write_labels(output_subject, fname, mris,
                               &pos_labels, &npos_labels, 0) ;
      MRISclearMarks(mris) ;
      mark_thresholded_vertices(mris, vbest_pvalues, vbest_avgs, -tthresh) ;
      segment_and_write_labels(output_subject, fname, mris, &neg_labels,
                               &nneg_labels, npos_labels) ;

      nlabels = nneg_labels + npos_labels ;
      if (nlabels) {
        labels = (LABEL **)calloc(nlabels, sizeof(LABEL *)) ;
        for (i = 0 ; i < npos_labels ; i++)
          labels[i] = pos_labels[i] ;
        for (i = 0 ; i < nneg_labels ; i++)
          labels[i+npos_labels] = neg_labels[i] ;
        free(pos_labels) ;
        free(neg_labels) ;
      }

      if (nlabels < MIN_LABELS) {
        if (FZERO(tthresh))
          break ;
        for (i = 0 ; i < nlabels ; i++)
          LabelFree(&labels[i]) ;
        if (nlabels)
          free(labels) ;
        tthresh *= 0.75 ;
        fprintf(stderr,
                "%d labels found (min %d), reducing threshold to %2.1f\n",
                nlabels, MIN_LABELS, tthresh) ;
      }
    } while (nlabels < MIN_LABELS) ;

    extract_thickness_at_best_scale(mris, c1_avg_curvs, vbest_avgs, c1_curvs,
                                    nvertices, num_class1);
    extract_thickness_at_best_scale(mris, c2_avg_curvs, vbest_avgs, c2_curvs,
                                    nvertices, num_class2);


    /* We have the thickness values at the most powerful scale stored for
       each subject in the c1_avg_curvs and c2_avg_curvs vectors.  Now collapse
       them across each label and compute values for t-test.
    */
    c1_thickness = (double **)calloc(num_class1, sizeof(double *)) ;
    c2_thickness = (double **)calloc(num_class2, sizeof(double *)) ;
    for (n = 0 ; n < num_class1 ; n++)
      c1_thickness[n] = (double *)calloc(nlabels, sizeof(double)) ;
    for (n = 0 ; n < num_class2 ; n++)
      c2_thickness[n] = (double *)calloc(nlabels, sizeof(double)) ;

    fprintf(stderr, "collapsing thicknesses within labels for class 1\n") ;
    for (n = 0 ; n < num_class1 ; n++)
      for (i = 0 ; i < nlabels ; i++)
        c1_thickness[n][i] =
          cvector_average_in_label(c1_avg_curvs[n], labels[i], nvertices) ;
    fprintf(stderr, "collapsing thicknesses within labels for class 2\n") ;
    for (n = 0 ; n < num_class2 ; n++)
      for (i = 0 ; i < nlabels ; i++)
        c2_thickness[n][i] =
          cvector_average_in_label(c2_avg_curvs[n], labels[i], nvertices) ;


    cvector_clear(c1_best_mean, nvertices) ;
    cvector_clear(c2_best_mean, nvertices) ;
    cvector_clear(c1_best_var, nvertices) ;
    cvector_clear(c2_best_var, nvertices) ;
    c1_mean = cvector_alloc(nlabels) ;
    c2_mean = cvector_alloc(nlabels) ;
    c1_var = cvector_alloc(nlabels) ;
    c2_var = cvector_alloc(nlabels) ;
    for (i = 0 ; i < nlabels ; i++)  /* for each row */
    {
      double total, totalsq ;

      total = totalsq = 0.0 ;
      for (n = 0 ; n < num_class1 ; n++)  /* for each column */
      {
        total += c1_thickness[n][i] ;
        totalsq += (c1_thickness[n][i]*c1_thickness[n][i]) ;
      }
      c1_mean[i] = total / (float)num_class1 ;
      c1_var[i] = totalsq / (float)num_class1 - c1_mean[i]*c1_mean[i] ;
      cvector_expand_label(labels[i], c1_mean[i], c1_best_mean) ;
      cvector_expand_label(labels[i], c1_var[i], c1_best_var) ;
    }
    for (i = 0 ; i < nlabels ; i++)  /* for each row */
    {
      double total, totalsq ;

      total = totalsq = 0.0 ;
      for (n = 0 ; n < num_class2 ; n++)  /* for each column */
      {
        total += c2_thickness[n][i] ;
        totalsq += (c2_thickness[n][i]*c2_thickness[n][i]) ;
      }
      c2_mean[i] = total / (float)num_class2 ;
      c2_var[i] = totalsq / (float)num_class2 - c2_mean[i]*c2_mean[i] ;
      cvector_expand_label(labels[i], c2_mean[i], c2_best_mean) ;
      cvector_expand_label(labels[i], c2_var[i], c2_best_var) ;
    }
    cvector_compute_t_test(c1_best_mean, c1_best_var, c2_best_mean,c2_best_var,
                           num_class1, num_class2, vbest_pvalues, nvertices) ;
  }

  MRISimportValVector(mris, vbest_snr) ;
  sprintf(fname, "./%s.%s_best_snr", hemi,prefix) ;
  MRISwriteValues(mris, fname) ;
  MRISimportValVector(mris, vbest_avgs) ;
  sprintf(fname, "./%s.%s_best_avgs", hemi, prefix) ;
  MRISwriteValues(mris, fname) ;
  fclose(fp) ;

  /* write out p values */
  MRISimportValVector(mris, vbest_pvalues) ;
  sprintf(fname, "./%s.%s_sig%d_%dt", hemi,prefix, cond_no1, cond_no2) ;
  MRISwriteValues(mris, fname) ;

  /* write out means */
  MRISimportValVector(mris, c1_best_mean) ;
  sprintf(fname, "./%s.%s_sigavg%d", hemi,prefix, cond_no1) ;
  MRISwriteValues(mris, fname) ;
  MRISimportValVector(mris, c2_best_mean) ;
  sprintf(fname, "./%s.%s_sigavg%d", hemi,prefix, cond_no2) ;
  MRISwriteValues(mris, fname) ;

  /* write out variances */
  cvector_divide(c1_best_var, (float)num_class1, nvertices) ;
  cvector_divide(c2_best_var, (float)num_class2, nvertices) ;
  MRISimportValVector(mris, c1_best_mean) ;
  sprintf(fname, "./%s.%s_sigvar%d", hemi,prefix, cond_no1) ;
  MRISwriteValues(mris, fname) ;
  MRISimportValVector(mris, c2_best_mean) ;
  sprintf(fname, "./%s.%s_sigvar%d", hemi,prefix, cond_no2) ;
  MRISwriteValues(mris, fname) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "multi-scale analysis took %d minutes and %d seconds.\n",
          minutes, seconds) ;
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
  else if (!stricmp(option, "max")) {
    max_avgs = atoi(argv[2]) ;
    fprintf(stderr,
            "computing kernel for maximum snr up to %d averages (%2.1f mm)\n",
            max_avgs, sqrt((double)max_avgs)) ;
    nargs = 1 ;
  } else if (!stricmp(option, "roi")) {
    roi_flag = 1 ;
    fprintf(stderr, "automatically generating regions of interest...\n") ;
  } else switch (toupper(*option)) {
    case 'C':
      cond_no1 = atoi(argv[2]) ;
      cond_no2 = atoi(argv[3]) ;
      nargs = 2 ;
      fprintf(stderr, "writing stats as condition numbers %d and %d\n",
              cond_no1, cond_no2) ;
      break ;
    case 'L':
      label_name = argv[2] ;
      fprintf(stderr, "masking label %s\n", label_name) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'M':
      min_label_area = atof(argv[2]) ;
      fprintf(stderr,
              "discarding labels with surface area smaller than %2.1f mm\n",
              min_label_area) ;
      nargs = 1 ;
      break ;
    case 'P':
      prefix = argv[2] ;
      fprintf(stderr, "using label prefix %s\n", prefix) ;
      nargs = 1 ;
      break ;
    case 'T':
      tthresh = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "using t threshold of %2.2f...\n", tthresh) ;
      break ;
    case 'W':
      write_flag = 1 ;
      break ;
    case 'O':
      output_subject = argv[2] ;
      fprintf(stderr, "using %s as output subject\n", output_subject) ;
      nargs = 1 ;
      break ;
    case 'N':
      use_no_distribution = 1 ;
      fprintf(stderr, "using distribution free estimate of snr...\n") ;
      break ;
    case 'B':
      bonferroni = 1 ;
      fprintf(stderr, "doing bonferroni correction of SNR values...\n") ;
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


static int
cvector_compute_variance(float *var,float *mean,int norm,int num) {
  int   i ;

  for (i = 0 ; i < num ; i++) {
    var[i] = var[i] / (float)norm - mean[i]*mean[i] ;
    if (var[i] < 0)  /* numerical instability */
      var[i] = 0 ;
  }
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
cvector_compute_t_test(float *c1_mean, float *c1_var,
                       float *c2_mean, float *c2_var,
                       int num_class1, int num_class2,
                       float *pvals, int num) {
  int    i ;
  double t, numer, denom, sig ;

  for (i = 0 ; i < num ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    numer = (c1_mean[i] - c2_mean[i]) ;
    denom = sqrt(c1_var[i]/num_class1 + c2_var[i]/num_class2) ;
    t = numer / denom ;
    if (FZERO(denom))
      sig = FZERO(numer) ? 1.0f : 0.0f ;
    else
      sig = sigt(t, num_class1+num_class2-2) ;
    if (numer > 0)
      pvals[i] = -log10(sig) ;
    else
      pvals[i] = log10(sig) ;
    if (!std::isfinite(numer) || !std::isfinite(denom) || !std::isfinite(pvals[i])
        || i == Gdiag_no)
      DiagBreak() ;
  }
  return(NO_ERROR) ;
}


#if 0
static int
cvector_subtract(float *v1, float *v2, float *vdst, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    vdst[i] = v1[i] - v2[i] ;
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
#endif

static double
cvector_compute_dist_free_snr(float **c1_curvs, int num_class1,
                              float **c2_curvs, int num_class2,
                              float *c1_mean, float *c2_mean, float *vsnr,
                              int num, int *pi) {
  int     i, max_i, n ;
  double  max_snr, mean, snr ;

  max_i = -1 ;
  for (max_snr = 0.0, i = 0 ; i < num ; i++) {
    mean = (c1_mean[i] + c2_mean[i])/2 ;
    snr = 0 ;
    if (c1_mean[i] > c2_mean[i]) {
      for (n = 0 ; n < num_class1 ; n++)
        if (c1_curvs[n][i] > mean)
          snr++ ;
      for (n = 0 ; n < num_class2 ; n++)
        if (c2_curvs[n][i] < mean)
          snr++ ;
    } else {
      for (n = 0 ; n < num_class1 ; n++)
        if (c1_curvs[n][i] < mean)
          snr++ ;
      for (n = 0 ; n < num_class2 ; n++)
        if (c2_curvs[n][i] > mean)
          snr++ ;
    }

    vsnr[i] = snr ;
    if (snr > max_snr)
      max_snr = snr ;
  }
  *pi = max_i ;
  return(max_snr) ;
}
static double
cvector_compute_snr(float *c1_mean, float *c2_mean, float *vvar, float *snr,
                    int num, int *pi, float bonferroni) {
  int    i, max_i ;
  float  f, max_snr ;

  max_i = -1 ;
  for (max_snr = i = 0 ; i < num ; i++) {
    f = (c1_mean[i]-c2_mean[i]) ;
    f *= f ;
    if (!iszero(vvar[i]))
      f /= (vvar[i]) ;

    f += bonferroni ;
    if (c2_mean[i] > c1_mean[i])
      f *= -1 ;   /* make it a signed quantity */

    if (fabs(f) > max_snr) {
      max_snr = fabs(f) ;
      max_i = i ;
    }
    snr[i] = f ;
  }
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
    vtotal_var[i] = (c1_var[i]*num_class1 + c2_var[i]*num_class2) / total_dof ;

  return(NO_ERROR) ;
}


#if 0
static int
cvector_track_best_snr(float *vsnr, float *vbest_snr, float *vbest_avgs,
                       int avgs, int num) {
  int    i ;

  for (i = 0 ; i < num ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    if (fabs(vsnr[i]) > fabs(vbest_snr[i])) {
      vbest_snr[i] = vsnr[i] ;
      vbest_avgs[i] = avgs ;
    }
  }

  return(NO_ERROR) ;
}
#endif

static int
cvector_track_best_stats(float *vsnr, float *vbest_snr, float *vbest_avgs,
                         float *c1_var, float *c2_var,
                         float *c1_best_var, float *c2_best_var,
                         float *c1_mean, float *c2_mean,
                         float *c1_best_mean, float *c2_best_mean,
                         int avgs, int num) {
  int    i ;

  for (i = 0 ; i < num ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    if (fabs(vsnr[i]) > fabs(vbest_snr[i])) {
      if (i == Gdiag_no)
        DiagBreak() ;
      if (i == Gdiag_no && avgs == 286)
        DiagBreak() ;
      vbest_snr[i] = vsnr[i] ;
      vbest_avgs[i] = avgs ;
      c1_best_var[i] = c1_var[i] ;
      c2_best_var[i] = c2_var[i] ;
      c1_best_mean[i] = c1_mean[i] ;
      c2_best_mean[i] = c2_mean[i] ;
    }
  }

  return(NO_ERROR) ;
}


static int
cvector_copy(float *v1, float *v2, int num) {
  int i ;

  for (i = 0 ; i < num ; i++)
    v2[i] = v1[i] ;

  return(NO_ERROR) ;
}

#if 0
static int
cvector_extract_best_avg(float *vbest_avgs,
                         float *vsrc, float *vdst, int navgs, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++) {
    if (nint(vbest_avgs[i]) == navgs)
      vdst[i] = vsrc[i] ;
  }
  return(NO_ERROR) ;
}

static int
extract_thickness_at_best_scale(MRI_SURFACE *mris, float **c1_avg_curvs,
                                float *vbest_avgs, float **c1_curvs,
                                int num, int nvectors) {
  int    i, max_avgs, avgs, n ;

  for (max_avgs = i = 0 ; i < num ; i++)
    if (nint(vbest_avgs[i]) >= max_avgs)
      max_avgs = nint(vbest_avgs[i]) ;

  for (avgs = 0 ; avgs <= max_avgs ; avgs++) {
    for (n = 0 ; n < nvectors ; n++) {
      cvector_extract_best_avg(vbest_avgs, c1_curvs[n], c1_avg_curvs[n],
                               avgs-1, num) ;
      MRISimportCurvatureVector(mris, c1_curvs[n]) ;
      MRISaverageCurvatures(mris, 1) ;
      MRISextractCurvatureVector(mris, c1_curvs[n]) ;
    }
  }
  return(NO_ERROR) ;
}
#endif

static int
cvector_bonferroni_correct(float *pvalues, float *avgs, int num) {
  int    i ;
  double correction ;

  for (i = 0 ; i < num ; i++) {
    correction = (double)num / (double)(avgs[i]+1) ;
    pvalues[i] *= correction ;
  }
  return(NO_ERROR) ;
}

static int
cvector_divide(float *v, float div, int num) {
  int  i ;

  for (i = 0 ; i < num ; i++)
    v[i] /= div ;
  return(NO_ERROR) ;
}


static int
mark_thresholded_vertices(MRI_SURFACE *mris,float *vbest_snr,
                          float *vbest_avgs, float thresh) {
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    if (vno == Gdiag_no)
      DiagBreak();
    v = &mris->vertices[vno] ;
    if (thresh > 0 && vbest_snr[vno] > thresh) {
      v->marked = 1 ;
      v->val = vbest_snr[vno] ;
      v->val2 = vbest_avgs[vno] ;
    } else if (thresh < 0 && vbest_snr[vno] < thresh) {
      v->marked = 1 ;
      v->val = vbest_snr[vno] ;
      v->val2 = vbest_avgs[vno] ;
    }
  }
  return(NO_ERROR) ;
}
static int
segment_and_write_labels(char *subject, char *name, MRI_SURFACE *mris,
                         LABEL ***plabels, int *pnlabels, int offset) {
  LABEL **labels, *area ;
  int   nlabels, i ;
  char  fname[STRLEN] ;

  MRISsegmentMarked(mris, &labels, &nlabels, min_label_area) ;

  for (i = 0 ; i < nlabels ; i++) {
    area = labels[i] ;
    strcpy(area->subject_name, subject) ;
    fprintf(stderr, "label %d: %d points, %2.1f mm\n",
            i, area->n_points, LabelArea(area, mris)) ;
    if (write_flag) {
      sprintf(fname, "%s%d", name, i+offset) ;
      fprintf(stderr, "writing label %s\n", fname) ;
      LabelWrite(area, fname) ;
    }
    strcpy(area->name, name) ;
  }
  *pnlabels = nlabels ;
  *plabels = labels ;
  return(NO_ERROR) ;
}

static int
cvector_extract_best_avg(float *vbest_avgs,
                         float *vsrc, float *vdst, int navgs, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++) {
    if (nint(vbest_avgs[i]) == navgs)
      vdst[i] = vsrc[i] ;
  }
  return(NO_ERROR) ;
}
static double
cvector_average_in_label(float *v, LABEL *area, int num) {
  int    i ;
  double avg ;

  for (avg = 0.0, i = 0 ; i < area->n_points ; i++)
    avg += v[area->lv[i].vno] ;
  avg /= (double)area->n_points ;
  return(avg) ;
}
static int
extract_thickness_at_best_scale(MRI_SURFACE *mris, float **c1_avg_curvs,
                                float *vbest_avgs, float **c1_curvs,
                                int num, int nvectors) {
  int    i, max_avgs, avgs, n ;

  for (max_avgs = i = 0 ; i < num ; i++)
    if (nint(vbest_avgs[i]) >= max_avgs)
      max_avgs = nint(vbest_avgs[i]) ;

  for (avgs = 0 ; avgs <= max_avgs ; avgs++) {
    for (n = 0 ; n < nvectors ; n++) {
      cvector_extract_best_avg(vbest_avgs, c1_curvs[n], c1_avg_curvs[n],
                               avgs-1, num) ;
      MRISimportCurvatureVector(mris, c1_curvs[n]) ;
      MRISaverageCurvatures(mris, 1) ;
      MRISextractCurvatureVector(mris, c1_curvs[n]) ;
    }
  }
  return(NO_ERROR) ;
}
static int
cvector_expand_label(LABEL *area, float f, float *vout) {
  int    i, vno ;

  for (i = 0 ; i < area->n_points ; i++) {
    vno = area->lv[i].vno ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vout[vno] = f ;
  }

  return(NO_ERROR) ;
}

