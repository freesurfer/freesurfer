/**
 * @file  mris_twoclass.c
 * @brief computes autocorrelation function of a curvature file
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:34 $
 *    $Revision: 1.12 $
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


#include "timer.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "const.h"
#include "proto.h"
#include "mrisurf.h"
#include "label.h"
#include "macros.h"
#include "fio.h"
#include "mrishash.h"
#include "sig.h"
#include "version.h"

static char vcid[] = "$Id: mris_twoclass.c,v 1.12 2011/03/02 00:04:34 nicks Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

#define STAT_T           0
#define STAT_F           1
#define STAT_MEAN        2
#define STAT_PCT         3
#define STAT_PSEUDO_T    4

/*-------------------------------- PROTOTYPES ----------------------------*/

static int which_norm = NORM_MEAN ;
static float sigma = 0.0f ;
static int stat_type = STAT_T ;
static int true_class = 1 ;
static int condition_0 = -1 ;
static int condition_1 = -1 ;
static int wfile_flag = 0 ;
static double conf = 0.0 ;
static int rectify_flag = 0 ;

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int write_vertex_data(char *fname, int index, float **v, int num) ;

static int   cvector_scalar_mul(float *v, float m, int num) ;
static float find_t_for_confidence_interval(float conf, int dof) ;
static double   cvector_compute_pvalues(float *c1_mean, float *c1_var,
                                        float *c2_mean, float *c2_var,
                                        int num_class1, int num_class2,
                                        float *pvals, int num, int stat_type,
                                        int *pvno) ;
static double cvector_average_in_label(float *v, LABEL *area, int num) ;
static int   cvector_extract_best_avg(float *vbest_avgs, float *vsrc,
                                      float *vdst, int navgs, int num) ;
static int   cvector_normalize(float *v, float norm, int num) ;
static int   cvector_accumulate(float *v, float *vtotal, int num) ;
static int   cvector_accumulate_square(float *v, float *vtotal, int num) ;
static int   cvector_compute_variance(float *var,float *mean,int norm,int num);
#if 0
static int   cvector_combine_variances(float *c1_var, float *c2_var, int num_class1,
                                       int num_class2, int nvertices, float *c_var) ;
#endif
static int   cvector_combine_variances_into_stderrs(float *c1_var, float *c2_var, int num_class1,
    int num_class2, int nvertices, float *c_var) ;
static double   cvector_compute_t_test(float *c1_mean, float *c1_var,
                                       float *c2_mean, float *c2_var,
                                       int num_class1, int num_class2,
                                       float *pvals, int num, int *pvno) ;
static double   cvector_compute_mean_diff(float *c1_mean, float *c2_mean,
    float *pvals, int num, int *pvno) ;
static double   cvector_compute_pct_diff(float *c1_mean, float *c2_mean,
    float *pvals, int num, int *pvno) ;
#if 0
static int  segment_and_write_labels(char *subject, char *fname,
                                     MRI_SURFACE *mris, LABEL ***plabels,
                                     int *pnlabels, int offset,
                                     float min_label_area) ;

static int  mark_thresholded_vertices(MRI_SURFACE *mris, float *vbest_snr,
                                      float *vbest_avgs, float thresh) ;
static int  extract_thickness_at_best_scale(MRI_SURFACE *mris,
    float **c1_avg_thickness,
    float *vbest_avgs,
    float **c1_thickness,
    int num, int nvectors);
static double    cvector_len(float *v, int num) ;
#endif
static int    cvector_copy(float *v1, float *v2, int num) ;
static double cvector_compute_dist_free_snr(float **c1_thickness,
    int num_class1,
    float **c2_thickness,
    int num_class2,
    float *c1_mean, float *c2_mean,
    float *vsnr, int num, int *pi);
static double cvector_compute_snr(float *c1_mean, float *c2_mean,
                                  float *verror, float *snr, int num, int *pi,
                                  float bonferroni, int stat_type);
static double cvector_compute_snr_F(float *c1_mean, float *c2_mean,
                                    float *verror, float *snr, int num, int *pi,
                                    float bonferroni);


#if 0
static int  *cvector_sort(float *vbest_snr, int nvertices) ;

static int   cvector_subtract(float *v1, float *v2, float *vdst, int num) ;
static int   cvector_mark_low_prob_vertices(float *pvals, float pthresh,
    MRI_SURFACE *mris) ;
#endif

static float *cvector_alloc(int num) ;
static int   cvector_clear(float *v, int num) ;
static int   cvector_set(float *v, float val, int num) ;
static int   cvector_add(float *v1, float *v2, float *vdst, int num) ;
static int   cvector_add_variances(float *c1_var, float *c2_var,
                                   int num_class1, int num_class2,
                                   float *total_var, int nvertices) ;
static int   cvector_multiply_variances(float *c1_var, float *c2_var,
                                        int num_class1, int num_class2,
                                        float *total_var, int nvertices) ;
static int   cvector_track_best_snr(float *vsnr, float *vbest_snr,
                                    float *vbest_avgs,
                                    float *c1_mean, float *c2_mean,
                                    float *c1_best_mean, float *c2_best_mean,
                                    float *c1_var, float *c2_var,
                                    float *c1_best_var, float *c2_best_var,
                                    float **c1_avg_thickness,
                                    float **c2_avg_thickness,
                                    float **c1_best_thicknesses, int nc1,
                                    float **c2_best_thicknesses, int nc2,
                                    int avgs, int num,
                                    float fthresh, int *pnum_found) ;

static int   cvector_track_best_stats(float *vsnr, float *vbest_snr,
                                      float *vbest_avgs,
                                      float *c1_mean, float *c2_mean,
                                      float *c1_best_mean, float *c2_best_mean,
                                      float *c1_var, float *c2_var,
                                      float *c1_best_var, float *c2_best_var,
                                      float **c1_avg_thickness,
                                      float **c2_avg_thickness,
                                      float **c1_best_thicknesses, int nc1,
                                      float **c2_best_thicknesses, int nc2,
                                      int avgs, int num,
                                      float fthresh, int *pnum_found) ;



/*-------------------------------- DATA ----------------------------*/

char *Progname ;

static int find_optimal_scale = TRUE ;
static int navgs = 0 ;   /* only used if find_optimal_scale is 0 */
static int nsort = -1 ;
static int use_buggy_snr = 0 ;
static char *write_dir = NULL ;
static char *read_dir = NULL ;

static int write_flag = 0 ;
static char *output_subject = NULL ;
static char *test_subject = NULL ;
static char *label_name = NULL ;
static char *prefix = "" ;

static int max_avgs = 5000 ;
static int use_no_distribution = 0 ;
static int use_stats = 1 ;

/* these do 77.8% on schizophrenic left hemispheres */
static float min_label_area = 25.0f ;
static float fthresh = 2.0 ;

#define MIN_LABELS   5
static int min_labels = MIN_LABELS ;

static int bonferroni = 0 ;
static FILE *labels_fp = NULL ;
static char *out_label_fname = NULL ;
static int normalize_flag = 0 ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  MRI_SURFACE  *mris ;
  char         **av, *curv_name, *surf_name, *hemi, fname[STRLEN],
  *cp, *subject_name, subjects_dir[STRLEN], *out_prefix,
  **c1_subjects, **c2_subjects, *stat_suffix ;
  int          ac, nargs, n, num_class1, num_class2, i, nvertices,
  avgs, max_snr_avgs, nlabels = 0, num_found = 0, total_found,
                                num_above_thresh ;
#if 0
  int          done, vno ;
#endif
  float        **c1_thickness, **c2_thickness, *curvs, *total_mean,
  *c1_mean, *c2_mean, *c1_best_mean, *c2_best_mean, *c1_best_var,
  *c2_best_var,
  *class_mean, *c1_var, *c2_var, *class_var,*pvals,
  **c1_avg_thickness,
  *vbest_snr, *vbest_avgs, *vtotal_var, *vsnr, **c2_avg_thickness,
  *vbest_pvalues ;
#if 0
  float        current_min_label_area, current_fthresh ;
#endif
  MRI_SP       *mrisp ;
  LABEL        *area, **labels = NULL ;
  FILE         *fp = NULL ;
  double       snr, max_snr ;
  struct timeb start ;
  int          msec, minutes, seconds ;
#if 0
  double       **c1_label_thickness, **c2_label_thickness ;
#endif
  int          *sorted_indices = NULL ;
  float        *test_thickness, *test_avg_thickness ;
  double       label_avg ;

  if (write_flag && DIAG_VERBOSE_ON)
    fp = fopen("scalespace.dat", "w") ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_twoclass.c,v 1.12 2011/03/02 00:04:34 nicks Exp $", "$Name:  $");
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

  TimerStart(&start) ;

  /* subject_name hemi surface curvature */
  if (argc < 8)
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
  out_prefix = argv[4] ;

  switch (stat_type) {
  default:
  case STAT_T:
    stat_suffix = "t" ;
    break ;
  case STAT_F:
    stat_suffix = "F" ;
    break ;
  case STAT_MEAN:
    stat_suffix = "m" ;
    break ;
  case STAT_PCT:
    stat_suffix = "p" ;
    break ;
  case STAT_PSEUDO_T:
    stat_suffix = "pt" ;
    break ;
  }
#define ARGV_OFFSET 5

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
  sprintf(fname, "%s/%s/surf/%s.%s",
          subjects_dir,output_subject,hemi,surf_name);
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  nvertices = mris->nvertices ;
  MRISfree(&mris) ;

  c1_best_mean = cvector_alloc(nvertices) ;
  c2_best_mean = cvector_alloc(nvertices) ;
  c1_best_var = cvector_alloc(nvertices) ;
  c2_best_var = cvector_alloc(nvertices) ;
  total_mean = cvector_alloc(nvertices) ;
  c1_mean = cvector_alloc(nvertices) ;
  pvals = cvector_alloc(nvertices) ;
  c2_mean = cvector_alloc(nvertices) ;
  c1_var = cvector_alloc(nvertices) ;
  c2_var = cvector_alloc(nvertices) ;


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
  c1_thickness = (float **)calloc(num_class1, sizeof(char *)) ;
  c1_avg_thickness = (float **)calloc(num_class1, sizeof(char *)) ;
  c2_subjects = (char **)calloc(num_class2, sizeof(char *)) ;
  c2_thickness = (float **)calloc(num_class2, sizeof(char *)) ;
  c2_avg_thickness = (float **)calloc(num_class2, sizeof(char *)) ;

  /* read in subject names for the two classes */
  for (n = 0 ; n < num_class1 ; n++) {
    c1_subjects[n] = argv[ARGV_OFFSET+n] ;
    c1_thickness[n] = (float *)cvector_alloc(nvertices) ;
    c1_avg_thickness[n] = cvector_alloc(nvertices) ;
    strcpy(c1_subjects[n], argv[ARGV_OFFSET+n]) ;
    /*    fprintf(stderr, "class1[%d] - %s\n", n, c1_subjects[n]) ;*/
  }
  i = n+1+ARGV_OFFSET ;  /* starting index */
  for (n = 0 ; n < num_class2 ; n++) {
    c2_subjects[n] = argv[i+n] ;
    c2_thickness[n] = cvector_alloc(nvertices) ;
    c2_avg_thickness[n] = cvector_alloc(nvertices) ;
    strcpy(c2_subjects[n], argv[i+n]) ;
    /*    fprintf(stderr, "class2[%d] - %s\n", n, c2_subjects[n]) ;*/
  }

  /* read in label if limitting calculation to an ROI */
  if (label_name) {
    area = LabelRead(output_subject, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname,
                label_name) ;
  } else
    area = NULL ;

  if (read_dir)  /* read in precomputed data */
  {
    sprintf(fname, "%s/%s/surf/%s.%s",
            subjects_dir,output_subject,hemi,surf_name);
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

    /* real all the curvatures in for group1 */
    for (n = 0 ; n < num_class1+num_class2 ; n++) {
      /* transform each subject's curvature into the output subject's space */
      subject_name = n < num_class1 ? c1_subjects[n]:c2_subjects[n-num_class1];
      fprintf(stderr, "reading subject %d of %d: %s\n",
              n+1, num_class1+num_class2, subject_name) ;
      sprintf(fname, "%s/%s.%s", read_dir,hemi,subject_name);
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(Gerror,
                  "%s: could not read curvature file %s",Progname,fname);
      if (area)
        MRISmaskNotLabel(mris, area) ;
      curvs = (n < num_class1) ? c1_thickness[n] : c2_thickness[n-num_class1] ;
      class_mean = (n < num_class1) ? c1_mean : c2_mean ;
      class_var = (n < num_class1) ? c1_var : c2_var ;
      MRISexportValVector(mris, curvs) ;
      cvector_accumulate(curvs, total_mean, nvertices) ;
      cvector_accumulate(curvs, class_mean, nvertices) ;
      cvector_accumulate_square(curvs, class_var, nvertices) ;
    }
  }
  else {
    /* real all the data for both groups  */
    for (n = 0 ; n < num_class1+num_class2 ; n++) {
      /* transform each subject's curvature into the output subject's space */
      subject_name = n < num_class1 ? c1_subjects[n]:c2_subjects[n-num_class1];
      fprintf(stderr, "reading subject %d of %d: %s\n",
              n+1, num_class1+num_class2, subject_name) ;
      sprintf(fname, "%s/%s/surf/%s.%s",
              subjects_dir,subject_name,hemi,surf_name);
      mris = MRISread(fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, fname) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
      if (strchr(curv_name, '/') != NULL)
        strcpy(fname, curv_name) ;  /* full path specified */
      else
        sprintf(fname,"%s/%s/surf/%s.%s",
                subjects_dir,subject_name,hemi,curv_name);
      if (wfile_flag) {
        if (MRISreadValues(mris, fname) != NO_ERROR)
          ErrorExit(Gerror,"%s: could no read w file file %s",Progname,fname);
      } else {
        if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
          ErrorExit(Gerror,"%s: could no read curvature file %s",
                    Progname,fname);
        if (normalize_flag)
          MRISnormalizeCurvature(mris, which_norm) ;
      }
      if (rectify_flag)
        MRISrectifyCurvature(mris) ;
      mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
      MRISfree(&mris) ;

      sprintf(fname, "%s/%s/surf/%s.%s",
              subjects_dir,output_subject,hemi,surf_name);
      mris = MRISread(fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, fname) ;
      MRISfromParameterization(mrisp, mris, 0) ;
      if (area)
        MRISmaskNotLabel(mris, area) ;
      curvs = (n < num_class1) ? c1_thickness[n] : c2_thickness[n-num_class1] ;
      class_mean = (n < num_class1) ? c1_mean : c2_mean ;
      class_var = (n < num_class1) ? c1_var : c2_var ;
      MRISextractCurvatureVector(mris, curvs) ;
      cvector_accumulate(curvs, total_mean, nvertices) ; /* across class */
      cvector_accumulate(curvs, class_mean, nvertices) ; /* within class */
      cvector_accumulate_square(curvs, class_var, nvertices) ;
      MRISPfree(&mrisp) ;
      MRISfree(&mris) ;
    }
  }

  /* compute within-group means, and total mean */
  cvector_normalize(total_mean, num_class1+num_class2, nvertices) ;
  cvector_normalize(c1_mean, num_class1, nvertices) ;
  cvector_normalize(c2_mean, num_class2, nvertices) ;
  cvector_compute_variance(c1_var, c1_mean, num_class1, nvertices) ;
  cvector_compute_variance(c2_var, c2_mean, num_class2, nvertices) ;

  /* read output surface */
  sprintf(fname, "%s/%s/surf/%s.%s",
          subjects_dir,output_subject,hemi,surf_name);
  fprintf(stderr, "reading output surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  if (area)  /* remove everything else */
    MRISripNotLabel(mris, area) ;
  vbest_snr = cvector_alloc(nvertices) ;
  vbest_pvalues = cvector_alloc(nvertices) ;
  vbest_avgs = cvector_alloc(nvertices) ;
  vtotal_var = cvector_alloc(nvertices) ;
  vsnr = cvector_alloc(nvertices) ;

  if (find_optimal_scale == TRUE) {
    if (read_dir == NULL)  /* recompute everything */
    {
      if (use_buggy_snr)
        cvector_multiply_variances(c1_var, c2_var, num_class1, num_class2,
                                   vtotal_var, nvertices) ;
      else
        cvector_add_variances(c1_var, c2_var, num_class1, num_class2,
                              vtotal_var, nvertices) ;
      if (use_no_distribution)
        snr = cvector_compute_dist_free_snr(c1_thickness, num_class1,
                                            c2_thickness, num_class2,
                                            c1_mean, c2_mean,
                                            vsnr, nvertices, &i);
      else if (use_stats) {
        snr = cvector_compute_pvalues(c1_mean, c1_var, c2_mean, c2_var,
                                      num_class1, num_class2, vsnr, nvertices,
                                      stat_type, &i) ;
        fprintf(stderr,
                "raw pval %2.2f, n=%2.4f, d=%2.4f, vno=%d\n",
                snr, c1_mean[i]-c2_mean[i], sqrt(vtotal_var[i]), i) ;
        cvector_track_best_stats(vsnr, vbest_snr, vbest_avgs,
                                 c1_mean, c2_mean, c1_best_mean, c2_best_mean,
                                 c1_var, c2_var, c1_best_var, c2_best_var,
                                 c1_thickness, c2_thickness,
                                 c1_thickness, num_class1,
                                 c2_thickness, num_class2,
                                 0, nvertices,
                                 fthresh,
                                 &num_found) ;
      } else {
        snr = cvector_compute_snr(c1_mean, c2_mean, vtotal_var, vsnr,nvertices,
                                  &i, 0.0f, stat_type);
        fprintf(stderr,
                "raw SNR %2.2f, n=%2.4f, d=%2.4f, vno=%d\n",
                sqrt(snr), c1_mean[i]-c2_mean[i], sqrt(vtotal_var[i]), i) ;
        cvector_track_best_snr(vsnr, vbest_snr, vbest_avgs,
                               c1_mean, c2_mean, c1_best_mean, c2_best_mean,
                               c1_var, c2_var, c1_best_var, c2_best_var,
                               c1_thickness, c2_thickness,
                               c1_thickness, num_class1,
                               c2_thickness, num_class2,
                               0, nvertices,
                               fthresh/((float)(num_class1+num_class2)/2),
                               &num_found) ;
      }

      max_snr = snr ;
      max_snr_avgs = 0 ;

      /* the c?_avg_thickness will be the thickness at the current scale */
      for (n = 0 ; n < num_class1 ; n++)
        cvector_copy(c1_thickness[n], c1_avg_thickness[n], nvertices) ;
      for (n = 0 ; n < num_class2 ; n++)
        cvector_copy(c2_thickness[n], c2_avg_thickness[n], nvertices) ;

      /* now incrementally average the data, keeping track of the best
         snr at each location, and at what scale it occurred. vbest_avgs
         and vbest_snr will contain the scale and the snr at that scale.
      */
      total_found = num_above_thresh = num_found ;
      for (avgs = 1 ; avgs <= max_avgs ; avgs++) {
        /* c?_avg_thickness is the thickness at the current scale */
        if (!(avgs % 50))
          fprintf(stderr, "testing %d averages...\n", avgs) ;
        cvector_clear(c1_mean, nvertices) ;
        cvector_clear(c2_mean, nvertices) ;
        cvector_clear(c1_var, nvertices) ;
        cvector_clear(c2_var, nvertices) ;
        for (n = 0 ; n < num_class1 ; n++) {
          MRISimportCurvatureVector(mris, c1_avg_thickness[n]) ;
          MRISaverageCurvatures(mris, 1) ;
          MRISextractCurvatureVector(mris, c1_avg_thickness[n]) ;
          cvector_accumulate(c1_avg_thickness[n], c1_mean, nvertices) ;
          cvector_accumulate_square(c1_avg_thickness[n], c1_var, nvertices) ;
        }
        for (n = 0 ; n < num_class2 ; n++) {
          MRISimportCurvatureVector(mris, c2_avg_thickness[n]) ;
          MRISaverageCurvatures(mris, 1) ;
          MRISextractCurvatureVector(mris, c2_avg_thickness[n]) ;
          cvector_accumulate(c2_avg_thickness[n], c2_mean, nvertices) ;
          cvector_accumulate_square(c2_avg_thickness[n], c2_var, nvertices) ;
        }

        /* c?_mean and c?_var are the means and variances at current scale */
        cvector_normalize(c1_mean, num_class1, nvertices) ;
        cvector_normalize(c2_mean, num_class2, nvertices) ;
        cvector_compute_variance(c1_var, c1_mean, num_class1, nvertices) ;
        cvector_compute_variance(c2_var, c2_mean, num_class2, nvertices) ;
        if (use_buggy_snr)
          cvector_multiply_variances(c1_var, c2_var, num_class1, num_class2,
                                     vtotal_var, nvertices) ;
        else
          cvector_add_variances(c1_var, c2_var, num_class1, num_class2,
                                vtotal_var, nvertices) ;
        if (use_no_distribution)
          snr =
            cvector_compute_dist_free_snr(c1_avg_thickness,num_class1,
                                          c2_avg_thickness, num_class2,c1_mean,
                                          c2_mean, vsnr, nvertices, &i);
        else if (use_stats) {
          snr = cvector_compute_pvalues(c1_mean, c1_var, c2_mean, c2_var,
                                        num_class1, num_class2, vsnr,nvertices,
                                        stat_type, &i);
          cvector_track_best_stats(vsnr, vbest_snr, vbest_avgs,
                                   c1_mean, c2_mean, c1_best_mean,c2_best_mean,
                                   c1_var, c2_var, c1_best_var, c2_best_var,
                                   c1_avg_thickness, c2_avg_thickness,
                                   c1_thickness, num_class1,
                                   c2_thickness, num_class2,
                                   avgs, nvertices,
                                   fthresh,
                                   &num_found) ;
        } else {
          snr =
            cvector_compute_snr(c1_mean, c2_mean, vtotal_var,vsnr,nvertices,&i,
                                bonferroni ? log((double)avgs):0.0f,stat_type);
          cvector_track_best_snr(vsnr, vbest_snr, vbest_avgs,
                                 c1_mean, c2_mean, c1_best_mean, c2_best_mean,
                                 c1_var, c2_var, c1_best_var, c2_best_var,
                                 c1_avg_thickness, c2_avg_thickness,
                                 c1_thickness, num_class1,
                                 c2_thickness, num_class2,
                                 avgs, nvertices,
                                 fthresh/((float)(num_class1+num_class2)/2),
                                 &num_found) ;
        }
        total_found += num_found ;
        if (write_flag && DIAG_VERBOSE_ON && !use_stats) {
          fprintf(fp, "%d %2.1f  %2.2f %2.2f %2.2f ",
                  avgs, sqrt((float)avgs), sqrt(snr), c1_mean[i]-c2_mean[i],
                  sqrt(vtotal_var[i])) ;
          fflush(fp) ;
          for (n = 0 ; n < num_class1 ; n++)
            fprintf(fp, "%2.2f ", c1_avg_thickness[n][i]) ;
          for (n = 0 ; n < num_class2 ; n++)
            fprintf(fp, "%2.2f ", c2_avg_thickness[n][i]) ;
          fprintf(fp, "\n") ;
          fclose(fp) ;
        }
        if (fabs(snr) > fabs(max_snr)) {
          if (use_stats)
            fprintf(stderr,
                    "new max pval found at avgs=%d (%2.1f mm)=%2.1f, "
                    "n=%2.1f-%2.1f, d=(%2.3f+%2.3f), vno=%d\n",
                    avgs, sqrt((float)avgs), snr, c1_mean[i],c2_mean[i],
                    sqrt(c1_var[i]), sqrt(c2_var[i]),i) ;
          else
            fprintf(stderr,
                    "new max SNR found at avgs=%d (%2.1f mm)=%2.1f, n=%2.4f, "
                    "d=%2.4f, vno=%d\n",
                    avgs, sqrt((float)avgs), sqrt(snr), c1_mean[i]-c2_mean[i],
                    sqrt(vtotal_var[i]), i) ;
          max_snr = snr ;
          max_snr_avgs = avgs ;
        }
        if (!(avgs%10)) {
          if (use_stats)
            printf("%d vertices found at scale %d, max p = %2.1f,"
                   "vno=%d,n=%2.1f-%2.1f, d=%2.2f+%2.2f\n",
                   total_found,avgs, snr,
                   i,c1_mean[i], c2_mean[i], c1_var[i], c2_var[i]);
          else
            printf("%d vertices found at scale %d, max snr=%2.2f\n",
                   total_found, avgs, snr) ;
          if (!total_found)
            break ;
          total_found = 0 ;
        }
      }

      snr = cvector_compute_pvalues(c1_best_mean, c1_best_var,
                                    c2_best_mean, c2_best_var,
                                    num_class1, num_class2, vbest_pvalues,
                                    nvertices, stat_type, NULL) ;

      for (n = 0 ; n < nvertices ; n++)
        if (vbest_avgs[n] > 0)
          num_above_thresh++ ;
      printf("max snr=%2.2f at %d averages, %d total vertices found "
             "above thresh\n", max_snr, max_snr_avgs, num_above_thresh);
      if (write_flag) {
        MRISimportValVector(mris, vbest_snr) ;
        sprintf(fname, "./%s.%s_best_snr", hemi,prefix) ;
        MRISwriteValues(mris, fname) ;
        MRISimportValVector(mris, vbest_avgs) ;
        sprintf(fname, "./%s.%s_best_avgs", hemi, prefix) ;
        MRISwriteValues(mris, fname) ;
        MRISimportValVector(mris, vbest_pvalues) ;
        sprintf(fname, "./%s.%s_best_pval", hemi,prefix) ;
        MRISwriteValues(mris, fname) ;
      }
    }
    else  /* read from directory containing precomputed optimal values */
    {
      sprintf(fname, "%s/%s.%s_best_snr", read_dir, hemi, prefix) ;
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: MRISreadValues(%s) failed",Progname,fname) ;
      MRISexportValVector(mris, vbest_snr) ;

      sprintf(fname, "%s/%s.%s_best_avgs", read_dir, hemi, prefix) ;
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: MRISreadValues(%s) failed",Progname,fname) ;
      MRISexportValVector(mris, vbest_avgs) ;
    }
  } else   /* use fixed scale smoothing kernel */
  {
    cvector_set(vbest_avgs, (float)navgs, nvertices) ;
    for (n = 0 ; n < num_class1 ; n++) {
      MRISimportCurvatureVector(mris, c1_thickness[n]) ;
      MRISaverageCurvatures(mris, navgs) ;
      MRISextractCurvatureVector(mris, c1_avg_thickness[n]) ;
      cvector_accumulate(c1_avg_thickness[n], c1_best_mean, nvertices) ;
      cvector_accumulate_square(c1_avg_thickness[n], c1_best_var, nvertices) ;
      cvector_copy(c1_avg_thickness[n], c1_thickness[n], nvertices) ;
    }
    for (n = 0 ; n < num_class2 ; n++) {
      MRISimportCurvatureVector(mris, c2_thickness[n]) ;
      MRISaverageCurvatures(mris, navgs) ;
      MRISextractCurvatureVector(mris, c2_avg_thickness[n]) ;
      cvector_accumulate(c2_avg_thickness[n], c2_best_mean, nvertices) ;
      cvector_accumulate_square(c2_avg_thickness[n], c2_best_var, nvertices) ;
      cvector_copy(c2_avg_thickness[n], c2_thickness[n], nvertices) ;
    }
    cvector_normalize(c1_best_mean, num_class1, nvertices) ;
    cvector_normalize(c2_best_mean, num_class2, nvertices) ;
    cvector_compute_variance(c1_best_var, c1_best_mean, num_class1, nvertices);
    cvector_compute_variance(c2_best_var, c2_best_mean, num_class2, nvertices);
    cvector_add_variances(c1_best_var, c2_best_var, num_class1, num_class2,
                          vtotal_var, nvertices) ;
    if (use_no_distribution)
      snr =
        cvector_compute_dist_free_snr(c1_avg_thickness,num_class1,
                                      c2_avg_thickness, num_class2,
                                      c1_best_mean, c2_best_mean,
                                      vsnr, nvertices, &i);
    else
      snr =
        cvector_compute_snr(c1_best_mean, c2_best_mean, vtotal_var,vsnr,
                            nvertices,&i, bonferroni ? log((double)navgs):0.0f,
                            stat_type);
  }

  /* now write out the optimal scale and statistic to out_prefix */
  cvector_compute_pvalues(c1_best_mean, c1_best_var, c2_best_mean,
                          c2_best_var, num_class1, num_class2,
                          pvals, nvertices, stat_type, NULL) ;

#if 0
  MRISimportValVector(mris, vbest_snr) ;
  sprintf(fname, "%s_snr-%s.w", out_prefix, hemi) ;
  MRISwriteValues(mris, fname) ;
#endif
  if (find_optimal_scale) {
    MRISimportValVector(mris, vbest_avgs) ;
    sprintf(fname, "%s_optimal_scale_%s-%s.w", out_prefix, stat_suffix, hemi) ;
    MRISwriteValues(mris, fname) ;
  }
  MRISimportValVector(mris, pvals) ;  /* was vbest_pvalues */
  sprintf(fname, "%s_%s-%s.w", out_prefix, stat_suffix, hemi) ;
  MRISwriteValues(mris, fname) ;
  if (conf > 0) {
    float  *c_stderr, *c_mean_diff, *c_conf_interval, tconf ;

    c_conf_interval = cvector_alloc(nvertices) ;
    c_stderr = cvector_alloc(nvertices) ;
    c_mean_diff = cvector_alloc(nvertices) ;

    cvector_compute_mean_diff(c1_best_mean, c2_best_mean, c_mean_diff, nvertices,NULL) ;
    cvector_combine_variances_into_stderrs(c1_best_var, c2_best_var, num_class1,
                                           num_class2,nvertices,c_stderr);
    tconf = find_t_for_confidence_interval(conf/2, num_class1+num_class2-2) ;
    printf("using t thresh of %2.2f\n", tconf) ;

    /* add sigma stderrs to mean diff */
    cvector_scalar_mul(c_stderr, tconf, nvertices) ;
    cvector_add(c_mean_diff, c_stderr, c_conf_interval, nvertices) ;
    MRISimportValVector(mris, c_conf_interval) ;
    sprintf(fname, "%s_%s_p%2.1fsigma-%s.w", out_prefix, stat_suffix, sigma, hemi) ;
    printf("writing mean diff + %2.1f sigma to %s...\n", sigma, fname) ;
    MRISwriteValues(mris, fname) ;

    /* subtract sigma stderrs to mean diff */
    cvector_scalar_mul(c_stderr, -1, nvertices) ;
    cvector_add(c_mean_diff, c_stderr, c_conf_interval, nvertices) ;
    MRISimportValVector(mris, c_conf_interval) ;  /* was vbest_pvalues */
    sprintf(fname, "%s_%s_m%2.1fsigma-%s.w", out_prefix, stat_suffix, sigma, hemi) ;
    printf("writing mean diff - %2.1f sigma to %s...\n", sigma, fname) ;
    MRISwriteValues(mris, fname) ;
  }

  if (labels_fp) {
    char   line[STRLEN] ;
    FILE   *out_fp ;
    LABEL  *area ;
    int    i ;
    double mean, var, thick ;

    printf("writing label report to %s...\n", out_label_fname) ;
    out_fp = fopen(out_label_fname, "w") ;
    if (!out_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open label report file %s",
                Progname, out_label_fname) ;

    while (fgetl(line, STRLEN, labels_fp)) {
      printf("reading label %s...\n", line) ;
      area = LabelRead(output_subject, line) ;
      if (!area) {
        fprintf(stderr, "%s: could not read label %s", Progname,line) ;
        continue ;
      }
      fprintf(out_fp, "%s  ", line) ;
      for (n = 0 ; n < num_class1 ; n++) {
        mean = var = 0.0 ;
        for (i = 0 ; i < area->n_points ; i++) {
          thick = c1_thickness[n][i] ;
          mean += thick ;
          var += (thick*thick) ;
        }
        mean /= (double)area->n_points ;
        var = var / (double)area->n_points - (mean*mean) ;
        fprintf(out_fp, "%2.3f %2.4f  ", mean, var) ;
      }
      for (n = 0 ; n < num_class2 ; n++) {
        mean = var = 0.0 ;
        for (i = 0 ; i < area->n_points ; i++) {
          thick = c2_thickness[n][i] ;
          mean += thick ;
          var += (thick*thick) ;
        }
        mean /= (double)area->n_points ;
        var = var / (double)area->n_points - (mean*mean) ;
        fprintf(out_fp, "%2.3f %2.4f  ", mean, var) ;
      }
      fprintf(out_fp, "\n") ;
      LabelFree(&area) ;
    }
    fclose(out_fp) ;
    fclose(labels_fp) ;
  }

  if (condition_0 >= 0)  /* write means and variances */
  {
    char  path[STRLEN] ;
    FILE  *fp ;

    FileNamePath(out_prefix, path) ;

    /* sigvar files actually are squared standard errors */
    cvector_scalar_mul(c1_best_var, 1.0f/(float)num_class1, nvertices) ;
    cvector_scalar_mul(c2_best_var, 1.0f/(float)num_class2, nvertices) ;

    /* write out dof files */
    sprintf(fname, "%s/sigavg%d.dof", path, condition_0) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_BADFILE,"%s: could not open dof file %s",Progname,fname);
    fprintf(fp, "%d\n", num_class1) ;
    fclose(fp) ;

    sprintf(fname, "%s/sigavg%d.dof", path, condition_1) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_BADFILE,"%s: could not open dof file %s",Progname,fname);
    fprintf(fp, "%d\n", num_class2) ;
    fclose(fp) ;

    /* write condition 0 stats */
    sprintf(fname, "%s/sigavg%d-%s.w", path, condition_0, hemi) ;
    MRISimportValVector(mris, c1_best_mean) ;
    MRISwriteValues(mris, fname) ;
    sprintf(fname, "%s/sigvar%d-%s.w", path, condition_0, hemi) ;
    MRISimportValVector(mris, c1_best_var) ;
    MRISwriteValues(mris, fname) ;

    /* write condition 1 stats */
    sprintf(fname, "%s/sigavg%d-%s.w", path, condition_1, hemi) ;
    MRISimportValVector(mris, c2_best_mean) ;
    MRISwriteValues(mris, fname) ;
    sprintf(fname, "%s/sigvar%d-%s.w", path, condition_1, hemi) ;
    MRISimportValVector(mris, c2_best_var) ;
    MRISwriteValues(mris, fname) ;
  }

  if (write_dir) {
    sprintf(fname, "%s/%s.%s_best_snr", write_dir, hemi,prefix) ;
    MRISimportValVector(mris, vbest_snr) ;
    if (MRISwriteValues(mris, fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: MRISwriteValues(%s) failed",Progname,fname) ;

    sprintf(fname, "%s/%s.%s_best_avgs", write_dir, hemi, prefix) ;
    MRISimportValVector(mris, vbest_avgs) ;
    if (MRISwriteValues(mris, fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: MRISwriteValues(%s) failed",Progname,fname) ;
  }

#if 0
  if (nsort < -1)
    nsort = mris->nvertices ;

  if (nsort <= 0) {
    nlabels = 0 ;
    current_min_label_area = min_label_area ;
    for (done = 0, current_fthresh = fthresh ;
         !FZERO(current_fthresh) && !done ;
         current_fthresh *= 0.95) {
      int   npos_labels, nneg_labels ;
      LABEL **pos_labels, **neg_labels ;

      for (current_min_label_area = min_label_area ;
           current_min_label_area > 0.5 ;
           current_min_label_area *= 0.75) {
        MRISclearMarks(mris) ;
        sprintf(fname, "%s-%s_thickness", hemi, prefix ? prefix : "") ;
        mark_thresholded_vertices(mris, vbest_snr, vbest_avgs,current_fthresh);
        segment_and_write_labels(output_subject, fname, mris,
                                 &pos_labels, &npos_labels, 0,
                                 current_min_label_area) ;
        MRISclearMarks(mris) ;
        mark_thresholded_vertices(mris, vbest_snr,vbest_avgs,-current_fthresh);
        segment_and_write_labels(output_subject, fname, mris, &neg_labels,
                                 &nneg_labels, npos_labels,
                                 current_min_label_area) ;

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
        done = (nlabels >= min_labels) ;
        if (done)  /* found enough points */
          break ;

        /* couldn't find enough  points - free stuff and try again */
        for (i = 0 ; i < nlabels ; i++)
          LabelFree(&labels[i]) ;
        if (nlabels)
          free(labels) ;
#if 0
        fprintf(stderr,"%d labels found (min %d), reducing constraints...\n",
                nlabels, min_labels) ;
#endif
      }
    }

    printf("%d labels found with F > %2.1f and area > %2.0f\n",
           nlabels, current_fthresh, current_min_label_area) ;
    for (i = 0 ; i < nlabels ; i++)
      fprintf(stderr, "label %d: %d points, %2.1f mm\n",
              i, labels[i]->n_points, LabelArea(labels[i], mris)) ;
  }

  /* read or compute thickness at optimal scale and put it into
     c?_avg_thickness.
  */
  if (!read_dir) {
    fprintf(stderr, "extracting thickness at optimal scale...\n") ;

    /* now build feature vectors for each subject */
    extract_thickness_at_best_scale(mris, c1_avg_thickness, vbest_avgs,
                                    c1_thickness, nvertices, num_class1);
    fprintf(stderr, "extracting thickness for class 2...\n") ;
    extract_thickness_at_best_scale(mris, c2_avg_thickness, vbest_avgs,
                                    c2_thickness, nvertices, num_class2);
  } else  /* read in precomputed optimal thicknesses */
  {
    char fname[STRLEN] ;

    fprintf(stderr, "reading precomputed thickness vectors\n") ;
    for (n = 0 ; n < num_class1 ; n++) {
      sprintf(fname, "%s/%s.%s", read_dir, hemi, argv[ARGV_OFFSET+n]) ;
      fprintf(stderr, "reading thickness vector from %s...\n", fname) ;
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read thickness file %s",
                  Progname,fname) ;
      MRISexportValVector(mris, c1_avg_thickness[n]) ;
    }
    for (n = 0 ; n < num_class2 ; n++) {
      sprintf(fname, "%s/%s.%s", read_dir, hemi,
              argv[n+num_class1+1+ARGV_OFFSET]) ;
      fprintf(stderr, "reading curvature vector from %s...\n", fname) ;
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read thickness file %s",
                  Progname,fname) ;
      MRISexportValVector(mris, c2_avg_thickness[n]) ;
    }
  }

  if (write_dir)   /* write out optimal thicknesses */
  {
    char fname[STRLEN] ;

    for (n = 0 ; n < num_class1 ; n++) {
      sprintf(fname, "%s/%s.%s", write_dir, hemi, argv[ARGV_OFFSET+n]) ;
      fprintf(stderr, "writing curvature vector to %s...\n", fname) ;
      MRISimportValVector(mris, c1_avg_thickness[n]) ;
      MRISwriteValues(mris, fname) ;
    }
    for (n = 0 ; n < num_class2 ; n++) {
      sprintf(fname, "%s/%s.%s", write_dir, hemi,
              argv[n+num_class1+1+ARGV_OFFSET]) ;
      fprintf(stderr, "writing curvature vector to %s...\n", fname) ;
      MRISimportValVector(mris, c2_avg_thickness[n]) ;
      MRISwriteValues(mris, fname) ;
    }
  }


  /* should free c?_thickness here */

  if (nsort <= 0) {
    /* We have the thickness values at the most powerful scale stored for
       each subject in the c1_avg_thickness and c2_avg_thickness vectors.
       Now collapse them across each label and build  feature vector for
       classification.
    */
    c1_label_thickness = (double **)calloc(num_class1, sizeof(double *)) ;
    c2_label_thickness = (double **)calloc(num_class2, sizeof(double *)) ;
    for (n = 0 ; n < num_class1 ; n++)
      c1_label_thickness[n] = (double *)calloc(nlabels, sizeof(double)) ;
    for (n = 0 ; n < num_class2 ; n++)
      c2_label_thickness[n] = (double *)calloc(nlabels, sizeof(double)) ;

    fprintf(stderr, "collapsing thicknesses within labels for class 1\n") ;
    for (n = 0 ; n < num_class1 ; n++)
      for (i = 0 ; i < nlabels ; i++)
        c1_label_thickness[n][i] =
          cvector_average_in_label(c1_avg_thickness[n], labels[i], nvertices) ;
    fprintf(stderr, "collapsing thicknesses within labels for class 2\n") ;
    for (n = 0 ; n < num_class2 ; n++)
      for (i = 0 ; i < nlabels ; i++)
        c2_label_thickness[n][i] =
          cvector_average_in_label(c2_avg_thickness[n], labels[i], nvertices) ;
    sprintf(fname, "%s_%s_class1.dat", hemi,prefix) ;
    fprintf(stderr, "writing class 1 info to %s...\n", fname) ;
    fp = fopen(fname, "w") ;
    for (i = 0 ; i < nlabels ; i++)  /* for each row */
    {
      for (n = 0 ; n < num_class1 ; n++)  /* for each column */
        fprintf(fp, "%2.2f  ", c1_label_thickness[n][i]) ;
      fprintf(fp, "\n") ;
    }
    fclose(fp) ;

    sprintf(fname, "%s_%s_class2.dat", hemi,prefix) ;
    fprintf(stderr, "writing class 2 info to %s...\n", fname) ;
    fp = fopen(fname, "w") ;
    for (i = 0 ; i < nlabels ; i++) {
      for (n = 0 ; n < num_class2 ; n++)
        fprintf(fp, "%2.2f  ", c2_label_thickness[n][i]) ;
      fprintf(fp, "\n") ;
    }
    fclose(fp) ;
  } else {
    sorted_indices = cvector_sort(vbest_snr, nvertices) ;
    vno = sorted_indices[0] ;
    write_vertex_data("c1.dat", vno, c1_avg_thickness,num_class1);
    write_vertex_data("c2.dat", vno, c2_avg_thickness,num_class2);
    printf("sorting complete\n") ;

    /* re-write class means at these locations */
    sprintf(fname, "%s_%s_class1.dat", hemi,prefix) ;
    fprintf(stderr, "writing class 1 info to %s...\n", fname) ;
    fp = fopen(fname, "w") ;
    for (i = 0 ; i < nsort ; i++) {
      for (n = 0 ; n < num_class1 ; n++)
        fprintf(fp, "%2.2f  ", c1_avg_thickness[n][sorted_indices[i]]) ;
      fprintf(fp, "\n") ;
    }
    fclose(fp) ;
    sprintf(fname, "%s_%s_class2.dat", hemi,prefix) ;
    fprintf(stderr, "writing class 2 info to %s...\n", fname) ;
    fp = fopen(fname, "w") ;
    for (i = 0 ; i < nsort ; i++) {
      for (n = 0 ; n < num_class2 ; n++)
        fprintf(fp, "%2.2f  ", c2_avg_thickness[n][sorted_indices[i]]) ;
      fprintf(fp, "\n") ;
    }
    fclose(fp) ;
  }

#endif

  if (test_subject) {
    test_thickness = cvector_alloc(nvertices) ;
    test_avg_thickness = cvector_alloc(nvertices) ;
    MRISfree(&mris) ;
    fprintf(stderr, "reading subject %s\n", test_subject) ;
    sprintf(fname, "%s/%s/surf/%s.%s",
            subjects_dir,test_subject,hemi,surf_name);
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    if (strchr(curv_name, '/') != NULL)
      strcpy(fname, curv_name) ;  /* full path specified */
    else
      sprintf(fname,"%s/%s/surf/%s.%s",
              subjects_dir,test_subject,hemi,curv_name);
    if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
      ErrorExit(Gerror,"%s: could no read curvature file %s",Progname,fname);
    if (normalize_flag)
      MRISnormalizeCurvature(mris, which_norm) ;
    if (rectify_flag)
      MRISrectifyCurvature(mris) ;

    mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
    MRISfree(&mris) ;

    sprintf(fname, "%s/%s/surf/%s.%s",
            subjects_dir,output_subject,hemi,surf_name);
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISfromParameterization(mrisp, mris, 0) ;
    if (area)
      MRISmaskNotLabel(mris, area) ;
    MRISextractCurvatureVector(mris, test_thickness) ;
    for (avgs = 0 ; avgs <= max_avgs ; avgs++) {
      cvector_extract_best_avg(vbest_avgs, test_thickness,test_avg_thickness,
                               avgs-1, nvertices) ;
      MRISimportCurvatureVector(mris, test_thickness) ;
      MRISaverageCurvatures(mris, 1) ;
      MRISextractCurvatureVector(mris, test_thickness) ;
    }

    if (nsort <= 0) {
      sprintf(fname, "%s_%s.dat", hemi,test_subject) ;
      fprintf(stderr, "writing test subject feature vector to %s...\n",
              fname) ;
      fp = fopen(fname, "w") ;
      for (i = 0 ; i < nlabels ; i++)  /* for each row */
      {
        label_avg =
          cvector_average_in_label(test_avg_thickness, labels[i], nvertices) ;
        fprintf(fp, "%2.2f\n", label_avg) ;
      }
      fclose(fp) ;
    } else   /* use sorting instead of connected areas */
    {
      double classification, offset, w ;
      int    total_correct, total_wrong, first_wrong, vno ;


      sprintf(fname, "%s_%s.dat", hemi,test_subject) ;
      fprintf(stderr, "writing test subject feature vector to %s...\n",
              fname) ;
      fp = fopen(fname, "w") ;

      first_wrong = -1 ;
      total_wrong = total_correct = 0 ;
      for (i = 0 ; i < nsort ; i++) {
        vno = sorted_indices[i] ;
        fprintf(fp, "%2.2f\n ", test_avg_thickness[sorted_indices[i]]) ;
        offset = (c1_mean[vno]+c2_mean[vno])/2.0 ;
        w = (c1_mean[vno]-c2_mean[vno]) ;
        classification = (test_avg_thickness[vno] - offset) * w ;

        if (((classification < 0) && (true_class == 1)) ||
            ((classification > 0) && (true_class == 2))) {
          total_wrong++ ;
          if (first_wrong < 0)
            first_wrong = i ;
        } else
          total_correct++ ;
      }
      fclose(fp) ;
      fprintf(stderr, "%d of %d correct = %2.1f%% (first wrong %d (%d)),"
              "min snr=%2.1f\n",
              total_correct, total_correct+total_wrong,
              100.0*total_correct / (total_correct+total_wrong),
              first_wrong, first_wrong >= 0 ? sorted_indices[first_wrong]:-1,
              vbest_snr[sorted_indices[nsort-1]]) ;

      if (first_wrong >= 0) {
        write_vertex_data("c1w.dat", sorted_indices[first_wrong],
                          c1_avg_thickness,num_class1);
        write_vertex_data("c2w.dat", sorted_indices[first_wrong],
                          c2_avg_thickness,num_class2);
      }
    }
  }

  msec = TimerStop(&start) ;
  free(total_mean);
  free(c1_mean) ;
  free(c2_mean) ;
  free(c1_var);
  free(c2_var);
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "cross-subject statistics took %d minutes and %d seconds.\n",
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
  else if (!stricmp(option, "test")) {
    test_subject = argv[2] ;
    fprintf(stderr, "writing test.dat for subject %s\n", test_subject) ;
    nargs = 1 ;
  } else if (!stricmp(option, "sigma")) {
    sigma = atof(argv[2]) ;
    fprintf(stderr, "writing out confidence intervals of +-%2.1f*stderr\n",sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "rectify") || !stricmp(option, "fabs")) {
    rectify_flag = 1 ;
    fprintf(stderr, "rectifying input vectors.\n") ;
  } else if (!stricmp(option, "conf")) {
    conf = atof(argv[2]) ;
    if (conf > 1)
      conf /= 100 ;   /* assume user used %% */
    fprintf(stderr, "writing out confidence intervals of %%%2.1f confidence interval\n",100*conf) ;
    conf = 1-conf ;
    nargs = 1 ;
  } else if (!stricmp(option, "labels")) {
    labels_fp = fopen(argv[2], "r") ;
    if (!labels_fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open label file list '%s'",
                Progname, argv[2]) ;
    out_label_fname = argv[3] ;
    fprintf(stderr, "reading label names out of %s and writing report to %s\n",
            argv[2], out_label_fname) ;
    nargs = 2 ;
  } else if (!stricmp(option, "wfile")) {
    wfile_flag = TRUE ;
    fprintf(stderr, "reading from w files\n") ;
  } else if (!stricmp(option, "optimal")) {
    find_optimal_scale = TRUE ;
    printf("finding optimal smoothing kernel at each cortical location\n") ;
  } else if (!stricmp(option, "fixed")) {
    find_optimal_scale = FALSE ;
    navgs = atoi(argv[2]) ;
    printf("using fixed size smoothing "
           " kernel (%d iterations) at each cortical location\n", navgs) ;
    nargs = 1 ;
  } else if (!stricmp(option, "num")) {
    min_labels = atoi(argv[2]) ;
    fprintf(stderr, "min labels set to %d\n", min_labels) ;
    nargs = 1 ;
  } else if (!stricmp(option, "mean")) {
    stat_type = STAT_MEAN ;
    fprintf(stderr, "computing mean difference between groups\n") ;
  } else if (!stricmp(option, "pct")) {
    stat_type = STAT_PCT ;
    fprintf(stderr, "computing pct difference between groups\n") ;
  } else if (!stricmp(option, "wt") || !stricmp(option, "write")) {
    write_dir = argv[2] ;
    nargs = 1 ;
    fprintf(stderr,"writing out optimal thickness vectors into directory %s\n",
            write_dir) ;
  } else if (!stricmp(option, "rt") || !stricmp(option, "read")) {
    read_dir = argv[2] ;
    nargs = 1 ;
    fprintf(stderr,"reading optimal thickness vectors from directory %s\n",
            read_dir) ;
  } else if (!stricmp(option, "normalize")) {
    normalize_flag = 1 ;
    fprintf(stderr, "normalizing inputs\n") ;
  } else if (!stricmp(option, "max")) {
    int i = max_avgs ;
    max_avgs = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting max_avgs to %d (was %d)\n", max_avgs, i);
  } else if (!stricmp(option, "conditions")) {
    nargs = 2 ;
    condition_0 = atoi(argv[2]) ;
    condition_1 = atoi(argv[3]) ;
    printf("writing summary statistics to condition %d and %d files...\n",
           condition_0, condition_1) ;
  } else switch (toupper(*option)) {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'L':
      label_name = argv[2] ;
      fprintf(stderr, "masking label %s\n", label_name) ;
      nargs = 1 ;
      break ;
    case 'C':
      true_class = atoi(argv[2]) ;
      fprintf(stderr, "generating stats for test subject as class %d\n",
              true_class) ;
      nargs = 1 ;
      break ;
    case 'S':
      nsort = atoi(argv[2]) ;
      fprintf(stderr, "sorting by SNR and using top %d to classify...\n",nsort) ;
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
      fthresh = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "using F snr threshold of %2.2f...\n", fthresh) ;
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
  printf("usage: %s -o <output subject> [options] \n"
         "\t<hemi> <surf> <curv> <out prefix> \n\t<c1_subject1> "
         "<c1_subject2>... : \n"
         "\t<c2_subject1> <c2_subject2>...\n",
         Progname) ;
  printf("where surf must be a spherical surface suitable for "
         "computing geodesics\n") ;
  printf("The <c1_subject> ... is a list of subjects from one class\n"
         "and the <c2_subject>... is a list of subjects from another "
         "class.\n");
  printf("valid options are:\n") ;
  printf("\t-wfile             -  not implemented yet (reading of *.w files)\n"
         "\t-optimal           -  compute optimal smoothing kernel\n"
         "\t-fixed <navgs>     -  apply a fixed smoothing kernel\n"
         "\t-conditions <c0> <c1>- write out sigavg<c?> and sigvar<c?> files\n"
         "\t-t <fthresh>       -  specify F threshold\n"
         "\t-o <subject>       -  specify output subject name\n"
         "\t-b                 -  perform bonferroni correction\n") ;
}

static void
print_help(void) {
  print_usage() ;
  printf(
    "\nThis program will compute the autocorrelation function of"
    " a curvature file\n") ;
  printf( "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  printf( "%s\n", vcid) ;
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

static double
cvector_compute_t_test(float *c1_mean, float *c1_var,
                       float *c2_mean, float *c2_var,
                       int num_class1, int num_class2,
                       float *pvals, int num, int *pvno) {
  int    i ;
  double t, numer, denom, p, max_p ;

  max_p = 0 ;
  for (i = 0 ; i < num ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    numer = (c1_mean[i] - c2_mean[i]) ;
    denom = sqrt((c1_var[i]/num_class1) + (c2_var[i]/num_class2)) ;
    if (FZERO(denom)) {
      t = 0 ;
      p = 0.0f ;
    } else {
      if (num_class1 == 1 || num_class2 == 1) {
        int dof_out = num_class1 + num_class2 ;
        /*
          if one of the means is not a population but a single subject
          then use standard deviation rather than standard error.
        */
        denom = sqrt((dof_out-1)*
                     ((c1_var[i]/num_class1) + (c2_var[i]/num_class2))) ;
      } else
        denom = sqrt((c1_var[i]/num_class1) + (c2_var[i]/num_class2)) ;
      t = numer / denom ;
      p = sigt(t, num_class1+num_class2-2) ;
      p = log10(p) ;
    }
    if (t > 0)
      p *= -1 ;
    pvals[i] = p ;

    if (fabs(p) > fabs(max_p)) {
      if (pvno)
        *pvno = i ;
      max_p = p ;
    }
  }
  return(max_p) ;
}
static double
cvector_compute_mean_diff(float *c1_mean, float *c2_mean,
                          float *vmean_diff, int num, int *pvno) {
  int    i ;
  double max_diff ;

  max_diff = 0 ;
  for (i = 0 ; i < num ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    vmean_diff[i] = (c1_mean[i] - c2_mean[i]) ;
    if (fabs(vmean_diff[i]) > fabs(max_diff)) {
      if (pvno)
        *pvno = i ;
      max_diff = vmean_diff[i] ;
    }
  }
  return(max_diff) ;
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
  printf( "%d vertices p < %2.2e\n", num, pthresh) ;
  return(NO_ERROR) ;
}
#endif

static double
cvector_compute_dist_free_snr(float **c1_thickness, int num_class1,
                              float **c2_thickness, int num_class2,
                              float *c1_mean, float *c2_mean, float *vsnr,
                              int num, int *pi) {
  int     i, max_i, n, correct, total ;
  double  max_snr, mean, snr ;

  max_i = -1 ;
  for (max_snr = 0.0, i = 0 ; i < num ; i++) {
    mean = (c1_mean[i] + c2_mean[i])/2 ;
    snr = 0 ;
    correct = 0 ;
    if (c1_mean[i] > c2_mean[i]) {
      for (n = 0 ; n < num_class1 ; n++)
        if (c1_thickness[n][i] > mean)
          correct++ ;

      for (n = 0 ; n < num_class2 ; n++)
        if (c2_thickness[n][i] < mean)
          correct++ ;
    } else {
      for (n = 0 ; n < num_class1 ; n++)
        if (c1_thickness[n][i] < mean)
          correct++ ;

      for (n = 0 ; n < num_class2 ; n++)
        if (c2_thickness[n][i] > mean)
          snr++ ;
    }

    total = num_class1 + num_class2 ;
    snr = (double)correct / (double)total ;
    vsnr[i] = snr ;
    if (snr > max_snr) {
      max_i = i ;
      max_snr = snr ;
    }
  }
  *pi = max_i ;
  return(max_snr) ;
}
static double
cvector_compute_snr(float *c1_mean, float *c2_mean, float *vvar, float *vsnr,
                    int num, int *pi, float bonferroni, int stat_type) {
  double snr ;

  switch (stat_type) {
  default:
  case STAT_T:
  case STAT_F:
    snr = cvector_compute_snr_F(c1_mean, c2_mean, vvar, vsnr, num,pi,
                                bonferroni);
    break ;
  }
  return(snr) ;
}
static double
cvector_compute_snr_F(float *c1_mean, float *c2_mean, float *vvar, float *snr,
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

static int
cvector_multiply_variances(float *c1_var, float *c2_var, int num_class1,
                           int num_class2, float *vtotal_var, int nvertices) {
  int     i, total_dof ;

  total_dof = num_class1 + num_class2 ;
  for (i = 0 ; i < nvertices ; i++)
    vtotal_var[i] = (c1_var[i]*num_class1 * c2_var[i]*num_class2) / total_dof ;

  return(NO_ERROR) ;
}


static int
cvector_track_best_snr(float *vsnr, float *vbest_snr, float *vbest_avgs,
                       float *c1_mean, float *c2_mean,
                       float *c1_best_mean, float *c2_best_mean,
                       float *c1_var, float *c2_var,
                       float *c1_best_var, float *c2_best_var,
                       float **c1_avg_thickness, float **c2_avg_thickness,
                       float **c1_best_thicknesses, int nc1,
                       float **c2_best_thicknesses, int nc2,
                       int avgs, int num, float fthresh, int *pnum_found) {
  int    i, n ;

  *pnum_found = 0 ;
  for (i = 0 ; i < num ; i++) {
    if (fabs(vsnr[i]) > fabs(vbest_snr[i])) {
      vbest_snr[i] = vsnr[i] ;
      vbest_avgs[i] = avgs ;
      c1_best_mean[i] = c1_mean[i] ;
      c1_best_var[i] = c1_var[i] ;
      c2_best_mean[i] = c2_mean[i] ;
      c2_best_var[i] = c2_var[i] ;
      if (vsnr[i] >= fthresh)
        *pnum_found += 1 ;
      for (n = 0 ; n < nc1 ; n++)
        c1_best_thicknesses[n][i] = c1_avg_thickness[n][i] ;
      for (n = 0 ; n < nc2 ; n++)
        c2_best_thicknesses[n][i] = c2_avg_thickness[n][i] ;
    }
  }

  return(NO_ERROR) ;
}
static int
cvector_track_best_stats(float *vpvals, float *vbest_pvals, float *vbest_avgs,
                         float *c1_mean, float *c2_mean,
                         float *c1_best_mean, float *c2_best_mean,
                         float *c1_var, float *c2_var,
                         float *c1_best_var, float *c2_best_var,
                         float **c1_avg_thickness, float **c2_avg_thickness,
                         float **c1_best_thicknesses, int nc1,
                         float **c2_best_thicknesses, int nc2,
                         int avgs, int num, float fthresh, int *pnum_found) {
  int    i, n ;

  *pnum_found = 0 ;
  for (i = 0 ; i < num ; i++) {
    if (fabs(vpvals[i]) > fabs(vbest_pvals[i])) {
      vbest_pvals[i] = vpvals[i] ;
      vbest_avgs[i] = avgs ;
      c1_best_mean[i] = c1_mean[i] ;
      c1_best_var[i] = c1_var[i] ;
      c2_best_mean[i] = c2_mean[i] ;
      c2_best_var[i] = c2_var[i] ;
      if (fabs(vpvals[i]) >= fthresh)
        *pnum_found += 1 ;
      for (n = 0 ; n < nc1 ; n++)
        c1_best_thicknesses[n][i] = c1_avg_thickness[n][i] ;
      for (n = 0 ; n < nc2 ; n++)
        c2_best_thicknesses[n][i] = c2_avg_thickness[n][i] ;
    }
  }

  return(NO_ERROR) ;
}

#if 0
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
                         LABEL ***plabels, int *pnlabels, int offset,
                         float min_label_area) {
  LABEL **labels, *area ;
  int   nlabels, i ;
  char  fname[STRLEN] ;

  MRISsegmentMarked(mris, &labels, &nlabels, min_label_area) ;

  for (i = 0 ; i < nlabels ; i++) {
    area = labels[i] ;
    strcpy(area->subject_name, subject) ;
#if 0
    printf( "label %d: %d points, %2.1f mm\n",
            i, area->n_points, LabelArea(area, mris)) ;
#endif
    if (write_flag) {
      sprintf(fname, "%s%d", name, i+offset) ;
      printf("writing label %s\n", fname) ;
      LabelWrite(area, fname) ;
    }
    strcpy(area->name, name) ;
  }
  *pnlabels = nlabels ;
  *plabels = labels ;
  return(NO_ERROR) ;
}
#endif

static int
cvector_copy(float *v1, float *v2, int num) {
  int i ;

  for (i = 0 ; i < num ; i++)
    v2[i] = v1[i] ;

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

  for (avg = 0.0, i = 0 ; i < area->n_points ; i++) {
    if (!finite(v[area->lv[i].vno]))
      DiagBreak() ;
    avg += v[area->lv[i].vno] ;
  }
  avg /= (double)area->n_points ;
  return(avg) ;
}

#if 0
static int
extract_thickness_at_best_scale(MRI_SURFACE *mris, float **c1_avg_thickness,
                                float *vbest_avgs, float **c1_thickness,
                                int num, int nvectors) {
  int    i, max_avgs, avgs, n ;

  for (max_avgs = i = 0 ; i < num ; i++)
    if (nint(vbest_avgs[i]) >= max_avgs)
      max_avgs = nint(vbest_avgs[i]) ;

  for (avgs = 0 ; avgs <= max_avgs ; avgs++) {
    for (n = 0 ; n < nvectors ; n++) {
      cvector_extract_best_avg(vbest_avgs, c1_thickness[n],
                               c1_avg_thickness[n], avgs-1, num) ;
      MRISimportCurvatureVector(mris, c1_thickness[n]) ;
      MRISaverageCurvatures(mris, 1) ;
      MRISextractCurvatureVector(mris, c1_thickness[n]) ;
    }
  }
  return(NO_ERROR) ;
}
#endif


#if 0
typedef struct {
  float  snr ;
  int    index ;
}
SORT_ELT ;

static int compare_sort_elts(const void *vse1, const void  *vse2) ;
static int
compare_sort_elts(const void *vse1, const void  *vse2) {
  const SORT_ELT *se1, *se2 ;

  se1 = vse1 ;
  se2 = vse2 ;
  return(se2->snr > se1->snr) ;
}

static int *
cvector_sort(float *vbest_snr, int nvertices) {
  int *sorted_vertices, i ;
  SORT_ELT  *se_table ;

  se_table = (SORT_ELT*)calloc(nvertices, sizeof(SORT_ELT)) ;
  if (!se_table)
    ErrorExit(ERROR_NOMEMORY, "cvector_sort: could not allocated %d vector",
              nvertices) ;
  sorted_vertices = (int *)calloc(nvertices, sizeof(int)) ;
  if (!sorted_vertices)
    ErrorExit(ERROR_NOMEMORY, "cvector_sort: could not allocated %d vector",
              nvertices) ;

  for (i = 0 ; i < nvertices ; i++) {
    se_table[i].index = i ;
    se_table[i].snr = vbest_snr[i] ;
  }

  qsort(se_table, nvertices, sizeof(SORT_ELT), compare_sort_elts) ;
  for (i = 0 ; i < nvertices ; i++)
    sorted_vertices[i] = se_table[i].index ;
  free(se_table) ;
  return(sorted_vertices) ;
}
#endif

static int
write_vertex_data(char *fname, int index, float **v, int nclass) {
  FILE  *fp ;
  int   i ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE, "write_vertex_data: could not open %s",fname));
  for (i = 0 ; i < nclass ; i++) {
    fprintf(fp, "%2.3f\n", v[i][index]) ;
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}
static double
cvector_compute_pvalues(float *c1_mean, float *c1_var, float *c2_mean,
                        float *c2_var, int num_class1, int num_class2,
                        float *pvals, int num, int stat_type, int *pvno) {
  switch (stat_type) {
  default:
  case STAT_T:
    return(cvector_compute_t_test(c1_mean, c1_var, c2_mean, c2_var,
                                  num_class1, num_class2, pvals, num,pvno)) ;
    break ;
  case STAT_MEAN:
    return(cvector_compute_mean_diff(c1_mean, c2_mean, pvals, num,pvno)) ;
  case STAT_PCT:
    return(cvector_compute_pct_diff(c1_mean, c2_mean, pvals, num,pvno)) ;
  }
  return(NO_ERROR) ;
}

static int
cvector_scalar_mul(float *v, float m, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    v[i] *= m ;

  return(NO_ERROR) ;
}

static int
cvector_set(float *v, float val, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++)
    v[i] = val ;

  return(NO_ERROR) ;
}

static double
cvector_compute_pct_diff(float *c1_mean, float *c2_mean, float *vpct_diff, int num, int *pvno) {
  int    i ;
  double max_diff, denom ;

  max_diff = 0 ;
  for (i = 0 ; i < num ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    denom = (fabs(c1_mean[i]) + fabs(c2_mean[i])) / 2 ;
    if (FZERO(denom))
      denom = 1 ;
    vpct_diff[i] = (c1_mean[i] - c2_mean[i]) / denom  ;
    if (fabs(vpct_diff[i]) > fabs(max_diff)) {
      if (pvno)
        *pvno = i ;
      max_diff = vpct_diff[i] ;
    }
  }
  return(max_diff) ;
}

#if 0
static int
cvector_combine_variances(float *c1_var, float *c2_var, int num_class1,
                          int num_class2, int nvertices, float *c_var) {
  int   i ;

  for (i = 0 ; i < nvertices ; i++) {
    c_var[i] = (c1_var[i] * num_class1 + c2_var[i]*num_class2) / (num_class1 + num_class2) ;
  }
  return(NO_ERROR) ;
}
#endif

static int
cvector_combine_variances_into_stderrs(float *c1_var, float *c2_var, int num_class1,
                                       int num_class2, int nvertices, float *c_stderr) {
  int   i ;

  for (i = 0 ; i < nvertices ; i++) {
    c_stderr[i] = sqrt((c1_var[i] * (num_class1-1) + c2_var[i]*(num_class2-1)) /
                       (num_class1 + num_class2 - 2)) ;
    c_stderr[i] = c_stderr[i] * sqrt(1.0f/(float)num_class1 + 1.0f / (float)num_class2) ;
  }
  return(NO_ERROR) ;
}

static int
cvector_add(float *v1, float *v2, float *vdst, int num) {
  int   i ;

  for (i = 0 ; i < num ; i++) {
    vdst[i] = v1[1] + v2[i] ;
  }
  return(NO_ERROR) ;
}

static float
find_t_for_confidence_interval(float conf, int dof) {
  double t, p ;

  for (t = 0.0 ; t += 0.01 ; t++) {
    p = sigt(t, dof) ;
    if (p < conf)
      return(p) ;
  }
  return(0.0) ;
}

