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
#include "macros.h"
#include "fio.h"
#include "mrishash.h"
#include "sig.h"

static char vcid[] = "$Id: mris_classify_thickness.c,v 1.3 2000/05/26 12:26:27 fischl Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

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


static int  cvector_compute_t(float *vbest_snr, float *vbest_pvalue, 
                              int dof, int num) ;
static double cvector_average_in_label(float *v, LABEL *area, int num) ;
static int   cvector_extract_best_avg(float *vbest_avgs, float *vsrc, 
                                      float *vdst, int navgs, int num) ;
static int   cvector_normalize(float *v, float norm, int num) ;
static int   cvector_accumulate(float *v, float *vtotal, int num) ;
static int   cvector_accumulate_square(float *v, float *vtotal, int num) ;
static int   cvector_compute_variance(float *var,float *mean,int norm,int num);
static int   cvector_compute_t_test(float *c1_mean, float *c1_var, 
                                    float *c2_mean, float *c2_var, 
                                    int num_class1, int num_class2, 
                                    float *pvals, int num) ;
#if 0
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
#endif

static float *cvector_alloc(int num) ;
static int   cvector_clear(float *v, int num) ;
static int   cvector_add_variances(float *c1_var, float *c2_var, 
                                   int num_class1, int num_class2,
                                   float *total_var, int nvertices) ;
static int   cvector_multiply_variances(float *c1_var, float *c2_var, 
                                   int num_class1, int num_class2,
                                   float *total_var, int nvertices) ;
static int   cvector_track_best_snr(float *vsnr, float *vbest_snr,
                                    float *vbest_avgs, int avgs, int num) ;



/*-------------------------------- DATA ----------------------------*/

char *Progname ;

static int use_buggy_snr = 0 ;
static char *write_dir = NULL ;
static char *read_dir = NULL ;
static int compute_stats = 0 ;

static int write_flag = 0 ;
static char *output_subject = NULL ;
static char *test_subject = NULL ;
static char *label_name = NULL ;
static char *prefix = "" ;

static int max_avgs = 500 ;
static int use_no_distribution = 0 ;

/* these do 77.8% on schizophrenic left hemispheres */
static float min_label_area = 25.0f ;
static float fthresh = 5.0 ;

#define MIN_LABELS   5

static int bonferroni = 0 ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[])
{
  MRI_SURFACE  *mris ;
  char         **av, *curv_name, *surf_name, *hemi, fname[STRLEN],
               *cp, *subject_name, subjects_dir[STRLEN],
               **c1_subjects, **c2_subjects ;
  int          ac, nargs, n, num_class1, num_class2, i, nvertices,
               avgs, max_snr_avgs, nlabels ;
  float        **c1_curvs, **c2_curvs, *curvs, *total_mean, *c1_mean, *c2_mean,
               *class_mean, *c1_var, *c2_var, *class_var,*pvals,**c1_avg_curvs,
               *vbest_snr, *vbest_avgs, *vtotal_var, *vsnr, **c2_avg_curvs,
               *vbest_pvalues ;
  MRI_SP       *mrisp ;
  LABEL        *area, **labels = NULL ;
  FILE         *fp = NULL ;
  double       snr, max_snr ;
  struct timeb start ;
  int          msec, minutes, seconds ;

  if (write_flag && DIAG_VERBOSE_ON)
    fp = fopen("scalespace.dat", "w") ;

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

  TimerStart(&start) ;

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
  num_class1 = 0 ; n = ARGV_OFFSET ;
  do
  {
    num_class1++ ;
    n++ ;
    if (argv[n] == NULL || n >= argc)
      ErrorExit(ERROR_BADPARM, "%s: must spectify ':' between class lists",
                Progname) ;
  } while(argv[n][0] != ':') ;

  /* find  # of vertices in output subject surface */
  sprintf(fname, "%s/%s/surf/%s.%s", 
          subjects_dir,output_subject,hemi,surf_name);
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

  num_class2 = 0 ; n++ ; /* skip ':' */
  if (n >= argc)
    ErrorExit(ERROR_BADPARM, "%s: class2 list empty", Progname) ;
  do
  {
    num_class2++ ;
    n++ ;
    if (n >= argc)
      break ;
  } while(argv[n] != NULL) ;

  fprintf(stderr, "%d subjects in class 1, %d subjects in class 2\n",
          num_class1, num_class2) ;

  c1_subjects = (char **)calloc(num_class1, sizeof(char *)) ;
  c1_curvs = (float **)calloc(num_class1, sizeof(char *)) ;
  c1_avg_curvs = (float **)calloc(num_class1, sizeof(char *)) ;
  c2_subjects = (char **)calloc(num_class2, sizeof(char *)) ;
  c2_curvs = (float **)calloc(num_class2, sizeof(char *)) ;
  c2_avg_curvs = (float **)calloc(num_class2, sizeof(char *)) ;
  for (n = 0 ; n < num_class1 ; n++)
  {
    c1_subjects[n] = argv[ARGV_OFFSET+n] ;
    c1_curvs[n] = (float *)calloc(nvertices, sizeof(float)) ;
    c1_avg_curvs[n] = (float *)calloc(nvertices, sizeof(float)) ;
    if (!c1_curvs[n] || !c1_avg_curvs[n])
      ErrorExit(ERROR_NOMEMORY, 
                "%s: could not allocate %dth list of %d curvatures",
                Progname, n, nvertices) ;

    strcpy(c1_subjects[n], argv[ARGV_OFFSET+n]) ;
    /*    fprintf(stderr, "class1[%d] - %s\n", n, c1_subjects[n]) ;*/
  }
  i = n+1+ARGV_OFFSET ;  /* starting index */
  for (n = 0 ; n < num_class2 ; n++)
  {
    c2_subjects[n] = argv[i+n] ;
    c2_curvs[n] = (float *)calloc(nvertices, sizeof(float)) ;
    c2_avg_curvs[n] = (float *)calloc(nvertices, sizeof(float)) ;
    if (!c2_curvs[n] || !c2_avg_curvs[n])
      ErrorExit(ERROR_NOMEMORY, 
                "%s: could not allocate %dth list of %d curvatures",
                Progname, n, nvertices) ;
    strcpy(c2_subjects[n], argv[i+n]) ;
    /*    fprintf(stderr, "class2[%d] - %s\n", n, c2_subjects[n]) ;*/
  }
  
  if (label_name)
  {
    area = LabelRead(output_subject, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname,
                label_name) ;
  }
  else
    area = NULL ;

  if (read_dir)
  {
    sprintf(fname, "%s/%s/surf/%s.%s", 
            subjects_dir,output_subject,hemi,surf_name);
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

    /* real all the curvatures in for group1 */
    for (n = 0 ; n < num_class1+num_class2 ; n++)
    {
      /* transform each subject's curvature into the output subject's space */
      subject_name = n < num_class1 ? c1_subjects[n]:c2_subjects[n-num_class1];
      fprintf(stderr, "reading subject %d of %d: %s\n",
              n+1, num_class1+num_class2, subject_name) ;
      sprintf(fname, "%s/%s.%s", read_dir,hemi,subject_name);
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(Gerror,"%s: could no read curvature file %s",Progname,fname);
      if (area)
        MRISmaskNotLabel(mris, area) ;
      curvs = (n < num_class1) ? c1_curvs[n] : c2_curvs[n-num_class1] ;
      class_mean = (n < num_class1) ? c1_mean : c2_mean ;
      class_var = (n < num_class1) ? c1_var : c2_var ;
      MRISexportValVector(mris, curvs) ;
      cvector_accumulate(curvs, total_mean, nvertices) ;
      cvector_accumulate(curvs, class_mean, nvertices) ;
      cvector_accumulate_square(curvs, class_var, nvertices) ;
    }
  }
  else
  {
    
    /* real all the curvatures in for group1 */
    for (n = 0 ; n < num_class1+num_class2 ; n++)
    {
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
      if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
        ErrorExit(Gerror,"%s: could no read curvature file %s",Progname,fname);
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
  }

  /* compute within-group means, and total mean */
  cvector_normalize(total_mean, num_class1+num_class2, nvertices) ;
  cvector_normalize(c1_mean, num_class1, nvertices) ;
  cvector_normalize(c2_mean, num_class2, nvertices) ;
  cvector_compute_variance(c1_var, c1_mean, num_class1, nvertices) ;
  cvector_compute_variance(c2_var, c2_mean, num_class2, nvertices) ;
  cvector_compute_t_test(c1_mean, c1_var, c2_mean, c2_var, 
                         num_class1, num_class2, pvals, nvertices) ;

  sprintf(fname, "%s/%s/surf/%s.%s", 
          subjects_dir,output_subject,hemi,surf_name);
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

  if (read_dir == NULL)
  {
    if (use_buggy_snr)
      cvector_multiply_variances(c1_var, c2_var, num_class1, num_class2,
                            vtotal_var, nvertices) ;
    else
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
    cvector_track_best_snr(vsnr, vbest_snr, vbest_avgs, 0, nvertices) ;
    
    for (n = 0 ; n < num_class1 ; n++)
      cvector_copy(c1_curvs[n], c1_avg_curvs[n], nvertices) ;
    for (n = 0 ; n < num_class2 ; n++)
      cvector_copy(c2_curvs[n], c2_avg_curvs[n], nvertices) ;
    for (avgs = 1 ; avgs <= max_avgs ; avgs++)
    {
      if (!(avgs % 50))
        fprintf(stderr, "testing %d averages...\n", avgs) ;
      cvector_clear(c1_mean, nvertices) ; cvector_clear(c2_mean, nvertices) ;
      cvector_clear(c1_var, nvertices) ; cvector_clear(c2_var, nvertices) ;
      cvector_clear(total_mean, nvertices) ;
      for (n = 0 ; n < num_class1 ; n++)
      {
#if 0
        fprintf(stderr, "processing subject %d of %d: %s\r",
                n+1, num_class1+num_class2, c1_subjects[n]) ;
#endif
        MRISimportCurvatureVector(mris, c1_avg_curvs[n]) ;
        MRISaverageCurvatures(mris, 1) ;
        MRISextractCurvatureVector(mris, c1_avg_curvs[n]) ;
        cvector_accumulate(c1_avg_curvs[n], total_mean, nvertices) ;
        cvector_accumulate(c1_avg_curvs[n], c1_mean, nvertices) ;
        cvector_accumulate_square(c1_avg_curvs[n], c1_var, nvertices) ;
      }
      for (n = 0 ; n < num_class2 ; n++)
      {
#if 0
        fprintf(stderr, "processing subject %d of %d: %s\r",
                n+1+num_class1, num_class1+num_class2, c2_subjects[n]) ;
#endif
        MRISimportCurvatureVector(mris, c2_avg_curvs[n]) ;
        MRISaverageCurvatures(mris, 1) ;
        MRISextractCurvatureVector(mris, c2_avg_curvs[n]) ;
        cvector_accumulate(c2_avg_curvs[n], total_mean, nvertices) ;
        cvector_accumulate(c2_avg_curvs[n], c2_mean, nvertices) ;
        cvector_accumulate_square(c2_avg_curvs[n], c2_var, nvertices) ;
      }
      cvector_normalize(total_mean, num_class1+num_class2, nvertices) ;
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
        snr = cvector_compute_dist_free_snr(c1_avg_curvs,num_class1,
                                            c2_avg_curvs, num_class2, c1_mean,
                                            c2_mean, vsnr, nvertices, &i);
      else
        snr = 
          cvector_compute_snr(c1_mean, c2_mean, vtotal_var, vsnr, nvertices,&i,
                              bonferroni ? log((double)avgs) : 0.0f);
      if (write_flag && DIAG_VERBOSE_ON)
      {
        fprintf(fp, "%d %2.1f  %2.2f %2.2f %2.2f ",
                avgs, sqrt((float)avgs), sqrt(snr), c1_mean[i]-c2_mean[i],
                sqrt(vtotal_var[i])) ;
        fflush(fp) ;
        for (n = 0 ; n < num_class1 ; n++)
          fprintf(fp, "%2.2f ", c1_avg_curvs[n][i]) ;
        for (n = 0 ; n < num_class2 ; n++)
          fprintf(fp, "%2.2f ", c2_avg_curvs[n][i]) ;
        fprintf(fp, "\n") ;
        fclose(fp) ;
      }
      if (snr > max_snr)
      {
        fprintf(stderr, 
                "new max SNR found at avgs=%d (%2.1f mm)=%2.1f, n=%2.4f, "
                "d=%2.4f, vno=%d\n",
                avgs, sqrt((float)avgs), sqrt(snr), c1_mean[i]-c2_mean[i],
                sqrt(vtotal_var[i]), i) ;
        max_snr = snr ;
        max_snr_avgs = avgs ;
      }
      cvector_track_best_snr(vsnr, vbest_snr, vbest_avgs, avgs, nvertices) ;
    }
    if (compute_stats)
      cvector_compute_t(vbest_snr, vbest_pvalues,num_class1+num_class2, 
                      nvertices) ;
    printf("max snr=%2.2f at %d averages\n", max_snr, max_snr_avgs) ;
    if (write_flag)
    {
      MRISimportValVector(mris, vbest_snr) ;
      sprintf(fname, "./%s.%s_best_snr", hemi,prefix) ; 
      MRISwriteValues(mris, fname) ;
      MRISimportValVector(mris, vbest_avgs) ;
      sprintf(fname, "./%s.%s_best_avgs", hemi, prefix) ; 
      MRISwriteValues(mris, fname) ;
      if (compute_stats)
      {
        MRISimportValVector(mris, vbest_pvalues) ;
        sprintf(fname, "./%s.%s_best_pval", hemi,prefix) ; 
        MRISwriteValues(mris, fname) ;
      }
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

  if (write_dir)
  {
    sprintf(fname, "%s/%s.%s_best_snr", write_dir, hemi,prefix) ; 
    MRISimportValVector(mris, vbest_snr) ;
    if (MRISwriteValues(mris, fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: MRISwriteValues(%s) failed",Progname,fname) ;

    sprintf(fname, "%s/%s.%s_best_avgs", write_dir, hemi, prefix) ; 
    MRISimportValVector(mris, vbest_avgs) ;
    if (MRISwriteValues(mris, fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: MRISwriteValues(%s) failed",Progname,fname) ;
  }

  do
  {
    int   npos_labels, nneg_labels ;
    LABEL **pos_labels, **neg_labels ;

    MRISclearMarks(mris) ;
    sprintf(fname, "%s-%s_thickness", hemi, prefix ? prefix : "") ;
    mark_thresholded_vertices(mris, vbest_snr, vbest_avgs, fthresh) ;
    segment_and_write_labels(output_subject, fname, mris, 
                             &pos_labels, &npos_labels, 0) ;
    MRISclearMarks(mris) ;
    mark_thresholded_vertices(mris, vbest_snr, vbest_avgs, -fthresh) ;
    segment_and_write_labels(output_subject, fname, mris, &neg_labels, 
                             &nneg_labels, npos_labels) ;

    nlabels = nneg_labels + npos_labels ;
    if (nlabels)
    {
      labels = (LABEL **)calloc(nlabels, sizeof(LABEL *)) ;
      for (i = 0 ; i < npos_labels ; i++)
        labels[i] = pos_labels[i] ;
      for (i = 0 ; i < nneg_labels ; i++)
        labels[i+npos_labels] = neg_labels[i] ;
      free(pos_labels) ; free(neg_labels) ;
    }

    if (nlabels < MIN_LABELS)
    {
      if (FZERO(fthresh))
        break ;
      for (i = 0 ; i < nlabels ; i++)
        LabelFree(&labels[i]) ;
      if (nlabels)
        free(labels) ;
      fthresh *= 0.75 ;
      fprintf(stderr,"%d labels found (min %d), reducing threshold to %2.1f\n",
              nlabels, MIN_LABELS, fthresh) ;
    }
  } while (nlabels < MIN_LABELS) ;

  if (!read_dir)
  {
    fprintf(stderr, "%d labels found - extracting thickness at optimal "
            "scale...\n", nlabels) ;

    /* now build feature vectors for each subject */
    extract_thickness_at_best_scale(mris, c1_avg_curvs, vbest_avgs, c1_curvs, 
                                    nvertices, num_class1);
    fprintf(stderr, "extracting thickness for class 2...\n") ;
    extract_thickness_at_best_scale(mris, c2_avg_curvs, vbest_avgs, c2_curvs, 
                                    nvertices, num_class2);
  }
  else  /* read in precomputed optimal thicknesses */
  {
    char fname[STRLEN] ;
    
    fprintf(stderr,"%d labels found - reading precomputed thickness vectors\n",
            nlabels) ;
    for (n = 0 ; n < num_class1 ; n++)
    {
      sprintf(fname, "%s/%s.%s", read_dir, hemi, argv[ARGV_OFFSET+n]) ;
      fprintf(stderr, "reading thickness vector from %s...\n", fname) ;
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read thickness file %s",
                  Progname,fname) ;
      MRISexportValVector(mris, c1_avg_curvs[n]) ;
    }
    for (n = 0 ; n < num_class2 ; n++)
    {
      sprintf(fname, "%s/%s.%s", read_dir, hemi, 
              argv[n+num_class1+1+ARGV_OFFSET]) ;
      fprintf(stderr, "reading curvature vector from %s...\n", fname) ;
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read thickness file %s",
                  Progname,fname) ;
      MRISexportValVector(mris, c2_avg_curvs[n]) ;
    }
  }

  if (write_dir)   /* write out optimal thicknesses */
  {
    char fname[STRLEN] ;
    
    for (n = 0 ; n < num_class1 ; n++)
    {
      sprintf(fname, "%s/%s.%s", write_dir, hemi, argv[ARGV_OFFSET+n]) ;
      fprintf(stderr, "writing curvature vector to %s...\n", fname) ;
      MRISimportValVector(mris, c1_avg_curvs[n]) ;
      MRISwriteValues(mris, fname) ;
    }
    for (n = 0 ; n < num_class2 ; n++)
    {
      sprintf(fname, "%s/%s.%s", write_dir, hemi, 
              argv[n+num_class1+1+ARGV_OFFSET]) ;
      fprintf(stderr, "writing curvature vector to %s...\n", fname) ;
      MRISimportValVector(mris, c2_avg_curvs[n]) ;
      MRISwriteValues(mris, fname) ;
    }
  }


  /* We have the thickness values at the most powerful scale stored for
     each subject in the c1_avg_curvs and c2_avg_curvs vectors.  Now collapse
     them across each label and build  feature vector for classification.
  */
  {
    double  **c1_thickness, **c2_thickness ;
    FILE   *fp ;

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
    sprintf(fname, "%s_%s_class1.dat", hemi,prefix) ;
    fprintf(stderr, "writing class 1 info to %s...\n", fname) ;
    fp = fopen(fname, "w") ;
    for (i = 0 ; i < nlabels ; i++)  /* for each row */
    {
      for (n = 0 ; n < num_class1 ; n++)  /* for each column */
        fprintf(fp, "%2.2f  ", c1_thickness[n][i]) ;
      fprintf(fp, "\n") ;
    }
    fclose(fp) ;

    sprintf(fname, "%s_%s_class2.dat", hemi,prefix) ;
    fprintf(stderr, "writing class 2 info to %s...\n", fname) ;
    fp = fopen(fname, "w") ;
    for (i = 0 ; i < nlabels ; i++)
    {
      for (n = 0 ; n < num_class2 ; n++)
        fprintf(fp, "%2.2f  ", c2_thickness[n][i]) ;
      fprintf(fp, "\n") ;
    }
    fclose(fp) ;

    free(total_mean); free(c1_mean) ; free(c2_mean) ;free(c1_var);free(c2_var);
    if (test_subject)
    {
      float *test_thickness, *test_avg_thickness ;
      double label_avg ;

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
        ErrorExit(Gerror, "%s: could no read curvature file %s",Progname,fname);
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
      for (avgs = 0 ; avgs <= max_avgs ; avgs++)
      {
        cvector_extract_best_avg(vbest_avgs, test_thickness,test_avg_thickness,
                                 avgs-1, nvertices) ;
        MRISimportCurvatureVector(mris, test_thickness) ;
        MRISaverageCurvatures(mris, 1) ;
        MRISextractCurvatureVector(mris, test_thickness) ;
      }

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
    }
  }
  
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "classification took %d minutes and %d seconds.\n", 
          minutes, seconds) ;
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
  else if (!stricmp(option, "test"))
  {
    test_subject = argv[2] ;
    fprintf(stderr, "writing test.dat for subject %s\n", test_subject) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "avgs"))
  {
    max_avgs = atoi(argv[2]) ;
    fprintf(stderr, 
            "computing kernel for maximum snr up to %d averages (%2.1f mm)\n",
            max_avgs, sqrt((double)max_avgs)) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "wt") || !stricmp(option, "write"))
  {
    write_dir = argv[2] ;
    nargs = 1 ;
    fprintf(stderr,"writing out optimal thickness vectors into directory %s\n",
           write_dir) ;
  }
  else if (!stricmp(option, "rt") || !stricmp(option, "read"))
  {
    read_dir = argv[2] ;
    nargs = 1 ;
    fprintf(stderr,"reading optimal thickness vectors from directory %s\n",
           read_dir) ;
  }
  else if (!stricmp(option, "stats"))
  {
    compute_stats = 1 ;
    fprintf(stderr, "computing multi-scale p values...\n") ;
  }
  else if (!stricmp(option, "bug"))
  {
    use_buggy_snr = 1 ;
    fprintf(stderr, "using multiplicative variance in snr calculations...\n") ;
  }
  else switch (toupper(*option))
  {
  case 'L':
    label_name = argv[2] ;
    fprintf(stderr, "masking label %s\n", label_name) ;
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
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
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
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program will compute the autocorrelation function of"
          " a curvature file\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


static int
cvector_compute_variance(float *var,float *mean,int norm,int num)
{
  int   i ;

  for (i = 0 ; i < num ; i++)
  {
    var[i] = var[i] / (float)norm - mean[i]*mean[i] ;
    if (var[i] < 0)  /* numerical instability */
      var[i] = 0 ;
  }
  return(NO_ERROR) ;
}
static int
cvector_normalize(float *v, float norm, int num)
{
  int   i ;

  for (i = 0 ; i < num ; i++)
    v[i] /= norm ;
  return(NO_ERROR) ;
}

static int   
cvector_accumulate_square(float *v, float *vtotal, int num)
{
  int   i ;

  for (i = 0 ; i < num ; i++)
    vtotal[i] += v[i]*v[i] ;
  return(NO_ERROR) ;
}
static int
cvector_accumulate(float *v, float *vtotal, int num)
{
  int   i ;

  for (i = 0 ; i < num ; i++)
    vtotal[i] += v[i] ;
  return(NO_ERROR) ;
}

static int
cvector_compute_t_test(float *c1_mean, float *c1_var, 
                       float *c2_mean, float *c2_var, 
                       int num_class1, int num_class2, 
                       float *pvals, int num)
{
  int    i ;
  double t, numer, denom ;

  for (i = 0 ; i < num ; i++)
  {
    numer = (c1_mean[i] - c2_mean[i]) ;
    denom = sqrt(c1_var[i]/num_class1) + sqrt(c2_var[i]/num_class2) ;
    t = numer / denom ;
    if (FZERO(denom))
      pvals[i] = 1.0f ;
    else
      pvals[i] = sigt(t, num_class1+num_class2-2) ;
  }
  return(NO_ERROR) ;
}


#if 0
static int
cvector_subtract(float *v1, float *v2, float *vdst, int num)
{
  int   i ;

  for (i = 0 ; i < num ; i++)
    vdst[i] = v1[i] - v2[i] ;
  return(NO_ERROR) ;
}
static int
cvector_mark_low_prob_vertices(float *pvals, float pthresh, MRI_SURFACE *mris)
{
  int    i, num ;

  for (num = i = 0 ; i < mris->nvertices ; i++)
  {
    if (pvals[i] < pthresh)
    {
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
                              int num, int *pi)
{
  int     i, max_i, n, correct, total ;
  double  max_snr, mean, snr ;

  max_i = -1 ;
  for (max_snr = 0.0, i = 0 ; i < num ; i++)
  {
    mean = (c1_mean[i] + c2_mean[i])/2 ;
    snr = 0 ; correct = 0 ;
    if (c1_mean[i] > c2_mean[i])
    {
      for (n = 0 ; n < num_class1 ; n++)
        if (c1_curvs[n][i] > mean)
          correct++ ;

      for (n = 0 ; n < num_class2 ; n++)
        if (c2_curvs[n][i] < mean)
          correct++ ;
    }
    else
    {
      for (n = 0 ; n < num_class1 ; n++)
        if (c1_curvs[n][i] < mean)
          correct++ ;

      for (n = 0 ; n < num_class2 ; n++)
        if (c2_curvs[n][i] > mean)
          snr++ ;
    }

    total = num_class1 + num_class2 ;
    snr = (double)correct / (double)total ;
    vsnr[i] = snr ;
    if (snr > max_snr)
    {
      max_i = i ;
      max_snr = snr ;
    }
  }
  *pi = max_i ;
  return(max_snr) ;
}
static double
cvector_compute_snr(float *c1_mean, float *c2_mean, float *vvar, float *snr, 
                    int num, int *pi, float bonferroni)
{
  int    i, max_i ;
  float  f, max_snr ;

  max_i = -1 ;
  for (max_snr = i = 0 ; i < num ; i++)
  {
    f = (c1_mean[i]-c2_mean[i]) ; f *= f ;
    if (!iszero(vvar[i]))
      f /= (vvar[i]) ;

    f += bonferroni ;
    if (c2_mean[i] > c1_mean[i])
      f *= -1 ;   /* make it a signed quantity */

    if (fabs(f) > max_snr)
    {
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
cvector_len(float *v, int num)
{
  int    i ;
  double len ;

  for (len = 0.0, i = 0 ; i < num ; i++)
    len += v[i]*v[i] ;
  len /= (double)num ; len = sqrt(len) ;
  return(len) ;
}
#endif
static float *
cvector_alloc(int num)
{
  float *v ;

  v = (float *)calloc(num, sizeof(float)) ;
  if (!v)
    ErrorExit(ERROR_NOMEMORY, "cvector_alloc(%d): calloc failed", num) ;
  return(v) ;
}

static int
cvector_clear(float *v, int num)
{
  int   i ;

  for (i = 0 ; i < num ; i++)
    v[i] = 0 ;
  return(NO_ERROR) ;
}
static int
cvector_add_variances(float *c1_var, float *c2_var, int num_class1, 
                      int num_class2, float *vtotal_var, int nvertices)
{
  int     i, total_dof ;

  total_dof = num_class1 + num_class2 ;
  for (i = 0 ; i < nvertices ; i++)
    vtotal_var[i] = (c1_var[i]*num_class1 + c2_var[i]*num_class2) / total_dof ;

  return(NO_ERROR) ;
}

static int
cvector_multiply_variances(float *c1_var, float *c2_var, int num_class1, 
                      int num_class2, float *vtotal_var, int nvertices)
{
  int     i, total_dof ;

  total_dof = num_class1 + num_class2 ;
  for (i = 0 ; i < nvertices ; i++)
    vtotal_var[i] = (c1_var[i]*num_class1 * c2_var[i]*num_class2) / total_dof ;

  return(NO_ERROR) ;
}


static int
cvector_track_best_snr(float *vsnr, float *vbest_snr, float *vbest_avgs, 
                       int avgs, int num)
{
  int    i ;
  
  for (i = 0 ; i < num ; i++)
  {
    if (fabs(vsnr[i]) > fabs(vbest_snr[i]))
    {
      vbest_snr[i] = vsnr[i] ;
      vbest_avgs[i] = avgs ;
    }
  }

  return(NO_ERROR) ;
}

static int
mark_thresholded_vertices(MRI_SURFACE *mris,float *vbest_snr, 
                          float *vbest_avgs, float thresh)
{
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak(); 
    v = &mris->vertices[vno] ;
    if (thresh > 0 && vbest_snr[vno] > thresh)
    {
      v->marked = 1 ;
      v->val = vbest_snr[vno] ;
      v->val2 = vbest_avgs[vno] ;
    }
    else if (thresh < 0 && vbest_snr[vno] < thresh)
    {
      v->marked = 1 ;
      v->val = vbest_snr[vno] ;
      v->val2 = vbest_avgs[vno] ;
    }
  }
  return(NO_ERROR) ;
}
static int  
segment_and_write_labels(char *subject, char *name, MRI_SURFACE *mris,
                         LABEL ***plabels, int *pnlabels, int offset)
{
  LABEL **labels, *area ;
  int   nlabels, i ;
  char  fname[STRLEN] ;

  MRISsegmentMarked(mris, &labels, &nlabels, min_label_area) ;

  for (i = 0 ; i < nlabels ; i++)
  {
    area = labels[i] ;
    strcpy(area->subject_name, subject) ;
    fprintf(stderr, "label %d: %d points, %2.1f mm\n", 
            i, area->n_points, LabelArea(area, mris)) ;
    if (write_flag)
    {
      sprintf(fname, "%s%d", name, i+offset) ;
      fprintf(stderr, "writing label %s\n", fname) ;
      LabelWrite(area, fname) ;
    }
    strcpy(area->name, name) ;
  }
  *pnlabels = nlabels ; *plabels = labels ;
  return(NO_ERROR) ;
}

static int
cvector_copy(float *v1, float *v2, int num)
{
  int i ;

  for (i = 0 ; i < num ; i++)
    v2[i] = v1[i] ;

  return(NO_ERROR) ;
}

static int
cvector_extract_best_avg(float *vbest_avgs, 
                         float *vsrc, float *vdst, int navgs, int num)
{
   int   i ;

   for (i = 0 ; i < num ; i++)
   {
     if (nint(vbest_avgs[i]) == navgs)
       vdst[i] = vsrc[i] ;
   }
   return(NO_ERROR) ;
}
static double
cvector_average_in_label(float *v, LABEL *area, int num)
{
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
                                int num, int nvectors)
{
  int    i, max_avgs, avgs, n ;

  for (max_avgs = i = 0 ; i < num ; i++)
    if (nint(vbest_avgs[i]) >= max_avgs)
      max_avgs = nint(vbest_avgs[i]) ;

  for (avgs = 0 ; avgs <= max_avgs ; avgs++)
  {
    for (n = 0 ; n < nvectors ; n++)
    {
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
cvector_compute_t(float *vbest_snr, float *vbest_pvalues, int dof, int num)
{
  int    i ;
  double sig ;

  for (i = 0 ; i < num ; i++)
  {
    sig = sigt(sqrt(dof*fabs(vbest_snr[i])), dof-2) ;
    if (vbest_snr[i] < 0)
      sig = log10(sig) ;
    else
      sig = -log10(sig) ;
    vbest_pvalues[i] = sig ;
  }
  return(NO_ERROR) ;
}

