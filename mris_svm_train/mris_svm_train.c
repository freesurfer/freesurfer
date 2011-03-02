/**
 * @file  mris_svm_train.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:34 $
 *    $Revision: 1.6 $
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
#include "macros.h"
#include "fio.h"
#include "mrishash.h"
#include "cvector.h"
#include "svm.h"
#include "version.h"

static char vcid[] = "$Id: mris_svm_train.c,v 1.6 2011/03/02 00:04:34 nicks Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;


/*-------------------------------- DATA ----------------------------*/

char *Progname ;

static char *wfile_name = NULL ;
static char *c1_name = "class1" ;
static char *c2_name = "class2" ;
static char *sdir = NULL ;

static char *output_subject = NULL ;
static char *test_subject = NULL ;
static char *label_name = NULL ;
static char *prefix = "" ;
static double tol = DEFAULT_SVM_TOL ;
static double svm_C = DEFAULT_SVM_C ;
static double momentum = 0.0 ;
static double rbf_sigma = 0.0 ;
static double poly_d = 0.0 ;
static int max_iter = 1000000 ;

static int navgs = 0 ;

#define TEST_SVM 0
#if TEST_SVM
float x_[20][2] = {
                    {
                      -1.0000,    1.0000
                    },
                    {-0.5000,    1.5000},
                    {0,    2.0000},
                    {0.8633,    3.2528},
                    {-2.1805,    2.2737},
                    {-0.6757,    4.0846},
                    {-1.2589,    3.8058},
                    {-0.8013,    3.6031},
                    {-1.3194,    4.3115},
                    {-1.5635,    2.0287},
                    {1.0000,   -1.0000},
#if 0
                    {1.0000,   -2.0000},
#else
                    {1.5000,   -0.5000},
#endif
                    {2.0000,         0},
                    {2.2501,   -1.1316},
                    {3.0240,   -1.0661},
                    {1.0503,   -2.0306},
                    {2.7714,   -1.1275},
                    {2.2957,   -1.9387},
                    {3.2611,   -2.2788},
                    {3.1396,    0.1683}
                  } ;
float y[20] = {
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1
              } ;
#endif

#define MAX_ANNOTATIONS 100
static int annotations[MAX_ANNOTATIONS] ;
static char *anames[MAX_ANNOTATIONS] ;
static int nannotations = 0 ;
static char *annot_name = "aparc" ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  SVM          *svm ;
  MRI_SURFACE  *mris ;
  char         **av, *curv_name, *surf_name, *hemi, fname[STRLEN],
    *cp, *subject_name, subjects_dir[STRLEN],
    **cv_subjects, *out_fname ;
  int          ac, nargs, n, nsubjects, num_class1, num_class2, i, nvertices,
    index ;
  float        **cv_thickness, *cv_outputs, **cv_avg_thickness ;
  MRI_SP       *mrisp ;
  LABEL        *area ;
  struct timeb start ;
  int          msec, minutes, seconds ;



#if TEST_SVM
  {
    float **x ;
    int   ntraining = 20, i, j, k ;
    double /*sum,*/ w[2], out ;
    x = (float **)calloc(ntraining, sizeof(float*)) ;
    for (i = 0 ; i < ntraining ; i++) {
      x[i] = (float *)calloc(ntraining, sizeof(float)) ;
      x[i][0] = x_[i][0] ;
      x[i][1] = x_[i][1] ;
    }

    for (n = 0 ; n < 2 ; n++)
    {
      if (n == 0)
      {
        svm = SVMalloc(2, "c1", "c2") ;
        svm->type = SVM_KERNEL_RBF ;
        svm->sigma = 2 ;
      }
      else
      {
        svm = SVMread("test.svm") ;
      }

      for (k = 0; k < 1 ; k++) {
        srandom(k);
        if (n == 0)
        {
          SVMtrain(svm, x, y, 20, svm_C, tol, 100000) ;
          SVMwrite(svm, "test.svm") ;
        }

#if 0
        sum = SVMconstraint(svm, y, ntraining) ;
        printf("constraint = %f\n", sum) ;
#endif
        for (j = 0 ; j < 2 ; j++) {
          w[j] = 0 ;
          for (i = 0 ; i < svm->nsupport ; i++)
            w[j] += svm->asupport[i] * svm->ysupport[i] * svm->xsupport[i][j] ;
        }
#if 1
        for (i = 0 ; i < ntraining ; i++) {
          out = SVMclassify(svm, x[i]) ;
          printf("%d (%2.1f,%2.1f) (%2.4f): %f (%f) --> %s\n",
                 i, x[i][0], x[i][1], svm->alpha[i],
                 y[i], out, y[i]*out > 0 ? "CORRECT" : "INCORRECT") ;
        }
#endif
        printf("%d support vectors found, weights = %2.3f, %2.3f\n",
               svm->nsupport, w[0], w[1]) ;
      }
    }
    exit(0) ;
  }
#endif


  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_svm_train.c,v 1.6 2011/03/02 00:04:34 nicks Exp $", "$Name:  $");
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
  if (argc < 7)
    usage_exit() ;
  if (output_subject == NULL)
    ErrorExit(ERROR_BADPARM,
              "output subject must be specified with -o <subject name>");

  if (sdir)
    cp = sdir ;
  else
    cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
              Progname) ;

  strcpy(subjects_dir, cp) ;

  out_fname = argv[argc-1] ;
  printf("writing output to %s...\n", out_fname) ;
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
  sprintf(fname, "%s/%s/surf/%s.%s",
          subjects_dir,output_subject,hemi,surf_name);
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  nvertices = mris->nvertices ;
  MRISfree(&mris) ;

  num_class2 = 0 ;
  n++ ; /* skip ':' */
  if (n >= argc)
    ErrorExit(ERROR_BADPARM, "%s: class2 list empty", Progname) ;
  do {
    num_class2++ ;
    n++ ;
    if (n >= (argc-1))  /* output file is last arg */
      break ;
  } while (argv[n] != NULL) ;

  nsubjects = num_class1+num_class2 ;
  printf("%d subjects in class 1, %d subjects in class 2, %d total\n",
         num_class1, num_class2, nsubjects) ;

  cv_subjects = (char **)calloc(nsubjects, sizeof(char *)) ;
  cv_thickness = (float **)calloc(nsubjects, sizeof(char *)) ;
  cv_avg_thickness = (float **)calloc(nsubjects, sizeof(float *)) ;
  cv_outputs = (float *)calloc(nsubjects, sizeof(float)) ;
  for (n = 0 ; n < num_class1 ; n++) {
    cv_outputs[n] = 1 ;
    cv_subjects[n] = argv[ARGV_OFFSET+n] ;
    cv_thickness[n] = (float *)calloc(nvertices, sizeof(float)) ;
    cv_avg_thickness[n] = (float *)calloc(nvertices, sizeof(float)) ;
    if (!cv_thickness[n] || !cv_avg_thickness[n])
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate %dth list of %d curvatures",
                Progname, n, nvertices) ;

    /*    fprintf(stderr, "class1[%d] - %s\n", n, cv_subjects[n]) ;*/
  }
  i = n+1+ARGV_OFFSET ;  /* starting index */
  for (n = 0 ; n < num_class2 ; n++) {
    index = num_class1+n ;
    cv_outputs[index] = -1 ;
    cv_subjects[index] = argv[i+n] ;
    cv_thickness[index] = (float *)calloc(nvertices, sizeof(float)) ;
    cv_avg_thickness[index] = (float *)calloc(nvertices, sizeof(float)) ;
    if (!cv_thickness[index] || !cv_avg_thickness[index] || !cv_subjects[index])
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate %dth list of %d curvatures",
                Progname, n, nvertices) ;
    /*    fprintf(stderr, "class2[%d] - %s\n", n, c2_subjects[n]) ;*/
  }

  if (label_name) {
    area = LabelRead(output_subject, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname,
                label_name) ;
  } else
    area = NULL ;

  /* read all the curvatures in for group1 */
  for (n = 0 ; n < nsubjects ; n++) {
    /* transform each subject's curvature into the output subject's space */
    subject_name = cv_subjects[n] ;
    fprintf(stderr, "reading subject %d of %d: %s\n",
            n+1, nsubjects, subject_name) ;
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
    if (nannotations > 0)
    {
      int vno, a, found ;
      VERTEX *v ;
      
      if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, 
                  "%s: could not read annot file %s for subject %s",
                  Progname, annot_name, subject_name) ;
      if (n == 0)  // first time fill in annot index table
      {
        for (a = 0 ; a < nannotations ; a++)
        {
          int index ;

          CTABfindName(mris->ct, anames[a], &index) ;
          CTABannotationAtIndex(mris->ct, index, &annotations[a]) ;
          printf("mapping annot %s to %d\n",
                 anames[a], annotations[a]) ;
        }
      }
      // rip all vertices that don't have one of the specified annotations
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;
        found = 0 ;
        for (a = 0 ; a < nannotations ; a++)
          if (v->annotation == annotations[a])
            found = 1 ;
        if (found == 0)
          v->ripflag = 1 ;
      }
    }
    if (navgs > 0)
      MRISaverageCurvatures(mris, navgs) ;
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
    MRISextractCurvatureVector(mris, cv_avg_thickness[n]) ;
    MRISPfree(&mrisp) ;
    MRISfree(&mris) ;
  }

  printf("training svm...\n") ;
  svm = SVMalloc(nvertices, c1_name, c2_name) ;
  if (rbf_sigma > 0) {
    svm->type = SVM_KERNEL_RBF;
    svm->sigma = rbf_sigma ;
  } else if (poly_d > 0) {
    svm->type = SVM_KERNEL_POLYNOMIAL;
    svm->sigma = poly_d ;
  }

  SVMtrain(svm, cv_avg_thickness, cv_outputs, nsubjects, svm_C, tol, max_iter);
  
  if (wfile_name) {
    int vno ;
    char  fname[STRLEN] ;
    
    sprintf(fname, "%s/%s/surf/%s.%s",
            subjects_dir,output_subject,hemi,surf_name);
    fprintf(stderr, "reading output surface %s...\n", fname) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    
    if (area)
      MRISripNotLabel(mris, area) ;
    
    for (vno = 0 ; vno < mris->nvertices ; vno++)
      mris->vertices[vno].val = (float)mris->nvertices*svm->w[vno] ;
    sprintf(fname, "./%s.%s.w", hemi, wfile_name) ;
    sprintf(fname, "./%s", wfile_name) ;
    MRISwriteValues(mris, fname) ;
  }
  
  printf("writing trained SVM to %s...\n", out_fname) ;
  SVMwrite(svm, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("svm training took %d minutes and %d seconds.\n",
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
  } else if (!stricmp(option, "sdir")) {
    sdir = argv[2] ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "c1")) {
    c1_name = argv[2] ;
    printf("using %s as name for class 1\n", c1_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "c2")) {
    c2_name = argv[2] ;
    printf("using %s as name for class 2\n", c2_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max")) {
    max_iter = atoi(argv[2]) ;
    printf("using max iter = %d...\n", max_iter) ;
    nargs = 1 ;
  } else if (!stricmp(option, "rbf")) {
    rbf_sigma = atof(argv[2]) ;
    printf("using RBF SVM with sigma = %f\n", rbf_sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "poly")) {
    poly_d = atof(argv[2]) ;
    printf("using polynomial SVM with dimension = %f\n", poly_d) ;
    nargs = 1 ;
  } else if (!stricmp(option, "tol")) {
    tol = atof(argv[2]) ;
    printf("using integration tolerance = %e\n", tol) ;
    nargs = 1 ;
  } else if (!stricmp(option, "aname")) {
    annot_name = argv[2] ;
    printf("using %s as name of annotation file\n", annot_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "annot")) {
    if (nannotations >= MAX_ANNOTATIONS)
      ErrorExit(ERROR_NOMEMORY, "%s: too many annotations specified (%d)",
                Progname, nannotations) ;
    anames[nannotations] = argv[2] ;
    printf("using annotation %s\n", anames[nannotations]) ;
    nannotations++;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'L':
      label_name = argv[2] ;
      fprintf(stderr, "masking label %s\n", label_name) ;
      nargs = 1 ;
      break ;
    case 'C':
      svm_C = atof(argv[2]) ;
      printf("using C=%f for svm slack variable weight\n", svm_C) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      printf("averaging values %d times\n", navgs) ;
      break ;
    case 'P':
      prefix = argv[2] ;
      fprintf(stderr, "using label prefix %s\n", prefix) ;
      nargs = 1 ;
      break ;
    case 'W':
      wfile_name = argv[2] ;
      printf("writing svm to w file %s...\n", wfile_name) ;
      nargs = 1 ;
      break ;
    case 'O':
      output_subject = argv[2] ;
      fprintf(stderr, "using %s as output subject\n", output_subject) ;
      nargs = 1 ;
      break ;
    case 'M':
      momentum = atof(argv[2]) ;
      printf("using momentum = %f for SVM training.\n", momentum) ;
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
          "\nThis program will train an SVM classifier based on thickness (or other curv) "
          " measures\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1)  ;
}
