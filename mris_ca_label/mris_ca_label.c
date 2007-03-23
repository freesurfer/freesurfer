/**
 * @file  mris_ca_label.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/03/23 20:00:07 $
 *    $Revision: 1.23 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "mrisurf.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "gcsa.h"
#include "transform.h"
#include "annotation.h"
#include "icosahedron.h"
#include "version.h"
#include "cma.h"

static char vcid[] =
  "$Id: mris_ca_label.c,v 1.23 2007/03/23 20:00:07 fischl Exp $";

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int postprocess(GCSA *gcsa, MRI_SURFACE *mris) ;

char *Progname ;
static void usage_exit(int code) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int which_norm = NORM_MEAN ;
static double MIN_AREA_PCT = 0.1 ;
static char *read_fname = NULL ;
static int nbrs = 2 ;
static int filter = 10 ;
static char *orig_name = "smoothwm" ;
static MRI  *mri_aseg ;

#if 0
static int normalize_flag = 0 ;
static int navgs = 5 ;
static char *curv_name = "curv" ;
static char *thickness_name = "thickness" ;
static char *sulc_name = "sulc" ;
#endif

static char subjects_dir[STRLEN] ;
extern char *gcsa_write_fname ;
extern int gcsa_write_iterations ;

static int novar = 0 ;
static int refine = 0;

int
main(int argc, char *argv[]) {
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp,*hemi,
  *canon_surf_name ;
  int          ac, nargs, i ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MRI_SURFACE  *mris ;
  GCSA         *gcsa ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
  }
  if (argc < 6)
    usage_exit(1) ;

  subject_name = argv[1] ;
  hemi = argv[2] ;
  canon_surf_name = argv[3] ;
  out_fname = argv[5] ;

  printf("%s\n",vcid);
  printf("  %s\n",MRISurfSrcVersion());
  fflush(stdout);

  printf("reading atlas from %s...\n", argv[4]) ;
  gcsa = GCSAread(argv[4]) ;
  if (!gcsa)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier from %s",
              Progname, argv[4]) ;

  sprintf(fname, "%s/%s/surf/%s.%s", subjects_dir,subject_name,hemi,orig_name);
  if (DIAG_VERBOSE_ON)
    printf("reading surface from %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s for %s",
              Progname, fname, subject_name) ;
  MRISresetNeighborhoodSize(mris, nbrs) ;
  mris->ct = gcsa->ct ; /* hack so that color table
                                   will get written into annot file */

  //set annotation table from the colortable
  set_atable_from_ctable(gcsa->ct);

  // read colortable from the gcsa if not already done
  //  if(gcsa->ct != NULL)
  //  read_named_annotation_table(gcsa->ct->fname);

  MRIScomputeSecondFundamentalForm(mris) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  if (MRISreadCanonicalCoordinates(mris, canon_surf_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read spherical coordinate system from %s for %s",
              Progname, canon_surf_name, subject_name) ;
#if 1
  for (i = gcsa->ninputs-1 ; i >= 0 ; i--) {
    printf("input %d: %s, flags %0x, avgs %d, name %s\n",
           i, gcsa->inputs[i].type == GCSA_INPUT_CURV_FILE ?
           "CURVATURE FILE" : "MEAN CURVATURE",
           gcsa->inputs[i].flags,
           gcsa->inputs[i].navgs, gcsa->inputs[i].fname) ;
    switch (gcsa->inputs[i].type) {
    case GCSA_INPUT_CURV_FILE:
      if (MRISreadCurvature(mris, gcsa->inputs[i].fname) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                  Progname, gcsa->inputs[i].fname, subject_name) ;
      break ;
    case GCSA_INPUT_CURVATURE:
      MRISuseMeanCurvature(mris) ;
      MRISaverageCurvatures(mris, gcsa->inputs[i].navgs) ;
      break ;
    }
    if (gcsa->inputs[i].flags & GCSA_NORMALIZE)
      MRISnormalizeCurvature(mris, which_norm) ;
    MRIScopyCurvatureToValues(mris) ;
    if (i == 2)
      MRIScopyCurvatureToImagValues(mris) ;
    else if (i == 1)
      MRIScopyValToVal2(mris) ;
  }
#else
  if (gcsa->ninputs > 2) {
    if (MRISreadCurvature(mris, thickness_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                Progname, thickness_name, subject_name) ;
    MRIScopyCurvatureToImagValues(mris) ;
  }
  if (gcsa->ninputs > 1 || sulconly) {
    if (MRISreadCurvature(mris, sulc_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                Progname, curv_name, subject_name) ;
    MRIScopyCurvatureToValues(mris) ;
    MRIScopyValToVal2(mris) ;
  }

  if (!sulconly) {
#if 0
    if (MRISreadCurvature(mris, curv_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                Progname, sulc_name, subject_name) ;
#else
  MRISuseMeanCurvature(mris) ;
  MRISaverageCurvatures(mris, navgs) ;
#endif
    MRIScopyCurvatureToValues(mris) ;
  }
#endif


  if (novar)
    GCSAsetCovariancesToIdentity(gcsa) ;

  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

  if (!read_fname) {
    printf("labeling surface...\n") ;
    GCSAlabel(gcsa, mris) ;
    if (Gdiag_no >= 0)
      printf("vertex %d: label %s\n",
             Gdiag_no,
             annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;
    if (mri_aseg)
      GCSArelabelWithAseg(gcsa, mris, mri_aseg) ;
    printf("relabeling using gibbs priors...\n") ;
    GCSAreclassifyUsingGibbsPriors(gcsa, mris) ;
    if (Gdiag_no >= 0)
      printf("vertex %d: label %s\n",
             Gdiag_no,
             annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;
    if (mri_aseg)
      GCSArelabelWithAseg(gcsa, mris, mri_aseg) ;
    postprocess(gcsa, mris) ;
    if (Gdiag_no >= 0)
      printf("vertex %d: label %s\n",
             Gdiag_no,
             annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;
    if (gcsa_write_iterations != 0) {
      char fname[STRLEN] ;
      sprintf(fname, "%s_post.annot", gcsa_write_fname) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRISwriteAnnotation(mris, fname) ;
    }
  } else {

    MRISreadAnnotation(mris, read_fname) ;
    if (refine != 0) {
      GCSAreclassifyUsingGibbsPriors(gcsa, mris) ;
      if (Gdiag_no >= 0)
        printf("vertex %d: label %s\n",
               Gdiag_no,
               annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;
      postprocess(gcsa, mris) ;
      if (Gdiag_no >= 0)
        printf("vertex %d: label %s\n",
               Gdiag_no,
               annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;
      if (gcsa_write_iterations != 0) {
        char fname[STRLEN] ;
        sprintf(fname, "%s_post.annot", gcsa_write_fname) ;
        printf("writing snapshot to %s...\n", fname) ;
        MRISwriteAnnotation(mris, fname) ;
      }

    }
  }

  MRISmodeFilterAnnotations(mris, filter) ;
  if (Gdiag_no >= 0)
    printf("vertex %d: label %s\n",
           Gdiag_no,
           annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;

  printf("writing output to %s...\n", out_fname) ;
  if (MRISwriteAnnotation(mris, out_fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not write annot file %s for %s",
              Progname, out_fname, subject_name) ;

  MRISfree(&mris) ;
  GCSAfree(&gcsa) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("classification took %d minutes and %d seconds.\n",
         minutes, seconds) ;
  exit(0) ;
  return(0) ;
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;
  char *gcsfile, *fsh, *outannot;
  int  icoorder,err;
  char tmpstr[2000];
  GCSA *gcsa;
  MRIS *ico;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "ml-annot")) {
    // Compute most-likely annotation labeling on ico, save, and exit
    // args: gcs icoorder outannot
    if (argc < 5) {
      printf("ERROR: -ml-annot requires 3 args: gcs icoorder outannot\n");
      exit(1);
    }
    gcsfile=argv[2];/* rel to FREESURFER_HOME/average,
                                       ?h.curvature.buckner40.filled.desikan_killiany.gcs */
    sscanf(    argv[3],"%d",&icoorder); // usually 7
    outannot = argv[4];  // absolute path to output
    printf("ML Label: %s %d %s\n",gcsfile,icoorder,outannot);
    ico = ReadIcoByOrder(icoorder, 100);
    if (ico == NULL) exit(1);
    fsh = getenv("FREESURFER_HOME");
    sprintf(tmpstr,"%s/average/%s",fsh,gcsfile);
    printf("Reading gcsa from %s\n",tmpstr);
    gcsa = GCSAread(tmpstr);
    if (gcsa == NULL) exit(1);
    ico->ct = gcsa->ct;
    printf("Building most likely labels\n");
    GCSAbuildMostLikelyLabels(gcsa,ico);
    printf("Filtering labels\n");
    MRISmodeFilterAnnotations(ico, 2);
    err = MRISwriteAnnotation(ico, outannot);
    if (err) exit(1);
    MRISfree(&ico);
    GCSAfree(&gcsa);
    exit(0);
  } else if (!stricmp(option, "SDIR")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  } else if (!stricmp(option, "aseg")) {
    mri_aseg = MRIread(argv[2]) ;
    nargs = 1 ;
    if (mri_aseg == NULL)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s", Progname, argv[2]);
    printf("using %s aseg volume to correct midline\n", argv[2]) ;
  } else if (!stricmp(option, "MINAREA")) {
    MIN_AREA_PCT = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting minimum area threshold for connectivity to %2.2f\n",
           MIN_AREA_PCT) ;
    if (MIN_AREA_PCT < 0 || MIN_AREA_PCT > 1)
      ErrorExit(ERROR_BADPARM, "%s: MIN_AREA_PCT must be in [0,1]\n",Progname);
  } else if (!stricmp(option, "ORIG")) {
    orig_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as original surface\n", orig_name) ;
  } else if (!stricmp(option, "LONG")) {
    refine = 1 ;
    printf("will refine the initial labeling read-in from -R \n") ;
  } else if (!stricmp(option, "NOVAR")) {
    novar = 1 ;
    printf("setting all covariance matrices to the identity...\n") ;
  }
#if 0
  else if (!stricmp(option, "NORM")) {
    normalize_flag = 1 ;
    printf("normalizing sulc after reading...\n") ;
  } else if (!stricmp(option, "SULC")) {
    sulconly = 1 ;
    printf("using sulc as only input....\n") ;
  }
#endif
  else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  } else if (!stricmp(option, "seed")) {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number generator to %d\n",
            atoi(argv[2])) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
#if 0
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
#endif
    case 'F':
      filter = atoi(argv[2]) ;
      nargs = 1 ;
      printf("applying mode filter %d times before writing...\n", filter) ;
      break ;
    case 'T':
      if (read_named_annotation_table(argv[2]) != NO_ERROR)
        ErrorExit(ERROR_BADFILE,
                  "%s: could not read annotation file %s...\n",
                  Progname, argv[2]) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      printf("printing diagnostic information about vertex %d\n", Gdiag_no) ;
      break ;
    case 'W':
      gcsa_write_iterations = atoi(argv[2]) ;
      gcsa_write_fname = argv[3] ;
      nargs = 2 ;
      printf("writing out snapshots of gibbs process every %d "
             "iterations to %s\n"
             ,gcsa_write_iterations, gcsa_write_fname) ;
      break ;
    case 'R':
      read_fname = argv[2] ;
      nargs = 1 ;
      printf("reading precomputed parcellation from %s...\n", read_fname) ;
      break ;
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
print_usage(void) {
  printf("mris_ca_label [options] <subject> <hemi> "
         "<canon surf> <classifier> <output file>\n");
  printf("\n");
  printf("   subject    - freesurfer subject id\n");
  printf("   hemi       - lh or rh\n");
  printf("   canonsurf  - cannoincal surface, usually ?h.sphere.reg\n");
  printf("   classifier - $FREESURFER_HOME/average/?h.curvature."
         "buckner40.filled.desikan_killiany.gcs\n");
  printf("   outputfile - ?h.aparc.annot\n");
  printf("\n");
  printf(" Options:\n");
  printf("\n");
  printf("  -ml-annot gcs icoorder annot : "
         "Compute most-likely annotation labeling on ico, save, and exit\n");
  printf("  -orig orig_name\n");
  printf("  -long\n");
  printf("  -nbrs\n");
  printf("  -a navgs\n");
  printf("  -f filter\n");
  printf("  -t annottable\n");
  printf("  -v diagno\n");
  printf("  -w fname\n");
  printf("  -r fname : precomputed parcellations\n");
  printf("  -seed N  : set random number generator to seed N\n");
  printf("\n");
}

static void
usage_exit(int code) {
  print_usage();
  exit(code) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\n"
          "Automatically assigns a neuroanatomical label to each location \n"
          "on a cortical surface model. This procedure incorporates both \n"
          "geometric information derived from the cortical model (sulcus \n"
          "and curvature), and neuroanatomical convention, as found in a \n"
          "training set (see mris_ca_train).\n\n");
  fprintf(stderr,
          "Required args:\n"
          "--------------\n\n") ;
  fprintf(stderr,
          "  <subject>          the subject id\n\n");
  fprintf(stderr,
          "  <hemi>             hemisphere: rh or lh\n\n");
  fprintf(stderr,
          "  <canonsurf>        canonical surface filename\n\n");
  fprintf(stderr,
          "  <classifier>       classifier array input filename\n\n");
  fprintf(stderr,
          "  <outputfile>       annotated surface output file\n\n");
  fprintf(stderr,
          "Optional args:\n"
          "--------------\n\n");
  fprintf(stderr,
          "  -sdir <directory>  use <directory> as subjects directory \n"
          "                     (default: $SUBJECTS_DIR)\n\n");
  fprintf(stderr,
          "  -orig <filename>   specify filename of original surface \n"
          "                     (default=smoothwm)\n\n");
  fprintf(stderr,
          "  -long              refines the initial labeling read-in \n"
          "                     from -r (default: disabled)\n\n");
  fprintf(stderr,
          "  -r <filename>      file containing precomputed parcellation\n\n");
  fprintf(stderr,
          "  -novar             sets all covariance matrices to the identify\n"
          "                     (default: disabled)\n\n");
  fprintf(stderr,
          "  -nbrs <number>     neighborhood size (default=2)\n\n");
  fprintf(stderr,
          "  -f <number>        applies mode filter <number> times before \n"
          "                     writing output (default=10)\n\n");
  fprintf(stderr,
          "  -t <filename>      specify parcellation table input file \n"
          "                     (default: none)\n\n");
  fprintf(stderr,
          "  -v <number>        diagnostic level (default=0)\n\n");
  fprintf(stderr,
          "  -w <number> <filename>  writes-out snapshots of gibbs process \n"
          "                          every <number> iterations to <filename>\n"
          "                          (default=disabled)\n\n");
  fprintf(stderr,
          "  --help             print help info\n\n");
  fprintf(stderr,
          "  --version          print version info\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
postprocess(GCSA *gcsa, MRI_SURFACE *mris) {
  LABEL       **larray, *area ;
  int         nlabels, i, j, annotation, n, nchanged, niter = 0, deleted ;
  double      max_area, label_area ;

#define MAX_ITER 5

  do {
    deleted = nchanged  = 0 ;
    MRISsegmentAnnotated(mris, &larray, &nlabels, 0) ;
    /*    printf("%d total segments in Gibbs annotation\n", nlabels) ;*/
    MRISclearMarks(mris) ;
    for (i = 0 ; i < nlabels ; i++) {
      area = larray[i] ;
      if (!area)   /* already processed */
        continue ;
      annotation = mris->vertices[area->lv[0].vno].annotation ;

      /* find label with this annotation with max area */
      max_area = LabelArea(area, mris) ;
      for (n = 1, j = i+1 ; j < nlabels ; j++) {
        if (!larray[j])
          continue ;
        if (annotation !=
            mris->vertices[larray[j]->lv[0].vno].annotation)
          continue ;
        n++ ;
        label_area = LabelArea(larray[j], mris) ;
        if (label_area > max_area)
          max_area = label_area ;
      }
#if 0
      printf("%03d: annotation %s (%d): %d segments, max area %2.1f\n",
             niter, annotation_to_name(annotation, NULL),
             annotation, n, max_area) ;
#endif
      for (j = i ; j < nlabels ; j++) {
        if (!larray[j])
          continue ;
        if (annotation !=
            mris->vertices[larray[j]->lv[0].vno].annotation)
          continue ;

        label_area = LabelArea(larray[j], mris) ;
        if (label_area < MIN_AREA_PCT*max_area) {
          /*          printf("relabeling annot %2.1f mm
                area...\n", label_area) ;*/
          nchanged += GCSAreclassifyLabel(gcsa, mris, larray[j]) ;
          deleted++ ;
        }
        LabelFree(&larray[j]) ;
      }
    }
    free(larray) ;
    printf("%03d: %d total segments, %d labels (%d vertices) changed\n",
           niter, nlabels, deleted, nchanged) ;
  } while (nchanged > 0 && niter++ < MAX_ITER) ;
  return(NO_ERROR) ;
}

