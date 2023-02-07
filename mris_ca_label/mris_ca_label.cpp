/**
 * @brief parcellate the cortex based on an atlas
 *
 * "Automatically Parcellating the Human Cerebral Cortex", Fischl et al.
 * (2004). Cerebral Cortex, 14:11-22.

 * "An automated labeling system for subdividing the human cerebral cortex
 * on MRI scans into gyral based regions of interest", Desikan et al.
 * (2006) NeuroImage 31(3):968-80.
 */
/*
 * Original Author: Bruce Fischl
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
#include <math.h>
#include <ctype.h>

#include "macros.h"

#include "mri.h"
#include "mrisurf.h"
#include "mrisurf_project.h"

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


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int postprocess(GCSA *gcsa, MRI_SURFACE *mris) ;

const char *Progname ;
static void usage_exit(int code) ;
static void print_help(void) ;
static void print_version(void) ;

static int which_norm = NORM_MEAN ;
static double MIN_AREA_PCT = 0.1 ;
static char *read_fname = NULL ;
static char *prob_fname = NULL ;
static int nbrs = 2 ;
static int filter = 10 ;
static const char *orig_name = "smoothwm" ;
static MRI  *mri_aseg ;
static const char *surf_dir = "surf" ;

#if 0
static int normalize_flag = 0 ;
static int navgs = 5 ;
static char *curv_name = "curv" ;
static const char *thickness_name = "thickness" ;
static const char *sulc_name = "sulc" ;
#endif

static char subjects_dir[STRLEN] ;
extern char *gcsa_write_fname ;
extern int gcsa_write_iterations ;

static int novar = 0 ;
static int refine = 0;

static LABEL *cortex_label = NULL ;
static int relabel_unknowns_with_cortex_label(GCSA *gcsa,
    MRI_SURFACE *mris,
    LABEL *cortex_label);

int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp,*hemi,
               *canon_surf_name ;
  int          ac, nargs, i ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI_SURFACE  *mris ;
  GCSA         *gcsa ;

  nargs = handleVersionOption(argc, argv, "mris_ca_label");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
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
  {
    usage_exit(1) ;
  }

  subject_name = argv[1] ;
  hemi = argv[2] ;
  canon_surf_name = argv[3] ;
  out_fname = argv[5] ;

  printf("%s\n",getVersion().c_str());
  printf("  %s\n",getVersion().c_str());
  fflush(stdout);

  printf("reading atlas from %s...\n", argv[4]) ;
  gcsa = GCSAread(argv[4]) ;
  if (!gcsa)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier from %s",
              Progname, argv[4]) ;

  int req = snprintf(fname, STRLEN, "%s/%s/%s/%s.%s", subjects_dir,subject_name,surf_dir,hemi,orig_name); 
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (DIAG_VERBOSE_ON)
  {
    printf("reading surface from %s...\n", fname) ;
  }
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
  for (i = gcsa->ninputs-1 ; i >= 0 ; i--)
  {
    printf("input %d: %s, flags %0x, avgs %d, name %s\n",
           i, gcsa->inputs[i].type == GCSA_INPUT_CURV_FILE ?
           "CURVATURE FILE" : "MEAN CURVATURE",
           gcsa->inputs[i].flags,
           gcsa->inputs[i].navgs, gcsa->inputs[i].fname) ;
    switch (gcsa->inputs[i].type)
    {
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
    {
      MRISnormalizeCurvature(mris, which_norm) ;
    }
    MRIScopyCurvatureToValues(mris) ;
    if (i == 2)
    {
      MRIScopyCurvatureToImagValues(mris) ;
    }
    else if (i == 1)
    {
      MRIScopyValToVal2(mris) ;
    }
  }
#else
  if (gcsa->ninputs > 2)
  {
    if (MRISreadCurvature(mris, thickness_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                Progname, thickness_name, subject_name) ;
    MRIScopyCurvatureToImagValues(mris) ;
  }
  if (gcsa->ninputs > 1 || sulconly)
  {
    if (MRISreadCurvature(mris, sulc_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s for %s",
                Progname, curv_name, subject_name) ;
    MRIScopyCurvatureToValues(mris) ;
    MRIScopyValToVal2(mris) ;
  }

  if (!sulconly)
  {
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
  {
    GCSAsetCovariancesToIdentity(gcsa) ;
  }

  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

  if (!read_fname)
  {
    printf("labeling surface...\n") ;
    GCSAlabel(gcsa, mris) ;
    if(Gdiag_no >= 0) printf("vertex %d: label %s\n", Gdiag_no,annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;
    if(mri_aseg) GCSArelabelWithAseg(gcsa, mris, mri_aseg) ;
    printf("relabeling using gibbs priors...\n") ;
    GCSAreclassifyUsingGibbsPriors(gcsa, mris) ;
    if (Gdiag_no >= 0)
      printf("vertex %d: label %s\n",
             Gdiag_no,
             annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;
    if (mri_aseg)
    {
      GCSArelabelWithAseg(gcsa, mris, mri_aseg) ;
    }
    postprocess(gcsa, mris) ;
    if (Gdiag_no >= 0)
      printf("vertex %d: label %s\n",
             Gdiag_no,
             annotation_to_name(mris->vertices[Gdiag_no].annotation, NULL)) ;
    if (gcsa_write_iterations != 0)
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s_post.annot", gcsa_write_fname) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRISwriteAnnotation(mris, fname) ;
    }
  }
  else
  {
    MRISreadAnnotation(mris, read_fname) ;
    if (refine != 0)
    {
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
      if (gcsa_write_iterations != 0)
      {
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

  if (cortex_label)
  {
    relabel_unknowns_with_cortex_label(gcsa, mris, cortex_label) ;
  }

  printf("writing output to %s...\n", out_fname) ;
  if (MRISwriteAnnotation(mris, out_fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not write annot file %s for %s",
              Progname, out_fname, subject_name) ;

  if (NULL != prob_fname){
    printf("WARNING: Saving probabilities does not work\n");
    // The posteriors need to be cached somewhere in the GCSAlabel() function. 
    // Using val2 breaks the labeling so this is no longer valid
    MRI* prob_image = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;
    for (i = 0 ; i < mris->nvertices ; i++)
      MRIsetVoxVal(prob_image, i, 0, 0, 0, mris->vertices[i].val2) ;
    printf("writing vertex label probabilities to %s...\n", prob_fname) ;
    MRIwrite(prob_image, prob_fname) ;
    MRIfree(&prob_image) ;
  }

  MRISfree(&mris) ;
  GCSAfree(&gcsa) ;
  msec = start.milliseconds() ;
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
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  char *gcsfile, *fsh, *outannot;
  int  icoorder,err;
  char tmpstr[2000];
  GCSA *gcsa;
  MRIS *ico;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "ml-annot"))
  {
    // Compute most-likely annotation labeling on ico, save, and exit
    // args: gcs icoorder outannot
    if (argc < 5)
    {
      printf("ERROR: -ml-annot requires 3 args: gcs icoorder outannot\n");
      exit(1);
    }
    gcsfile=argv[2];/* rel to FREESURFER_HOME/average,
                       ?h.curvature.buckner40.filled.desikan_killiany.gcs */
    sscanf(    argv[3],"%d",&icoorder); // usually 7
    outannot = argv[4];  // absolute path to output
    printf("ML Label: %s %d %s\n",gcsfile,icoorder,outannot);
    ico = ReadIcoByOrder(icoorder, 100);
    if (ico == NULL)
    {
      exit(1);
    }
    fsh = getenv("FREESURFER_HOME");
    sprintf(tmpstr,"%s/average/%s",fsh,gcsfile);
    printf("Reading gcsa from %s\n",tmpstr);
    gcsa = GCSAread(tmpstr);
    if (gcsa == NULL)
    {
      exit(1);
    }
    ico->ct = gcsa->ct;
    printf("Building most likely labels\n");
    GCSAbuildMostLikelyLabels(gcsa,ico);
    printf("Filtering labels\n");
    MRISmodeFilterAnnotations(ico, 2);
    err = MRISwriteAnnotation(ico, outannot);
    if (err)
    {
      exit(1);
    }
    MRISfree(&ico);
    GCSAfree(&gcsa);
    exit(0);
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  }
  else if (!stricmp(option, "aseg"))
  {
    mri_aseg = MRIread(argv[2]) ;
    nargs = 1 ;
    if (mri_aseg == NULL)
    {
      ErrorExit(ERROR_BADFILE, "%s: could not open %s", Progname, argv[2]);
    }
    printf("using %s aseg volume to correct midline\n", argv[2]) ;
  }
  else if (!stricmp(option, "surf_dir"))
  {
    surf_dir = argv[2] ;
    nargs = 1 ;
    printf("using %s instead of surf for subdirectory search\n", argv[2]) ;
  }
  else if (!stricmp(option, "MINAREA"))
  {
    MIN_AREA_PCT = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting minimum area threshold for connectivity to %2.2f\n",
           MIN_AREA_PCT) ;
    if (MIN_AREA_PCT < 0 || MIN_AREA_PCT > 1)
    {
      ErrorExit(ERROR_BADPARM, "%s: MIN_AREA_PCT must be in [0,1]\n",Progname);
    }
  }
  else if (!stricmp(option, "ORIG"))
  {
    orig_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as original surface\n", orig_name) ;
  }
  else if (!stricmp(option, "LONG"))
  {
    refine = 1 ;
    printf("will refine the initial labeling read-in from -R \n") ;
  }
  else if (!stricmp(option, "NOVAR"))
  {
    novar = 1 ;
    printf("setting all covariance matrices to the identity...\n") ;
  }
#if 0
  else if (!stricmp(option, "NORM"))
  {
    normalize_flag = 1 ;
    printf("normalizing sulc after reading...\n") ;
  }
  else if (!stricmp(option, "SULC"))
  {
    sulconly = 1 ;
    printf("using sulc as only input....\n") ;
  }
#endif
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  }
  else if (!stricmp(option, "seed"))
  {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number generator to %d\n",
            atoi(argv[2])) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
    {
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
    case 'L':
      cortex_label = LabelRead(NULL, argv[2]) ;
      if (cortex_label == NULL)
      {
        ErrorExit(ERROR_NOFILE, "") ;
      }
      nargs = 1 ;
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
    case 'P':
      prob_fname = argv[2] ;
      nargs = 1 ;
      printf("saving vertex label probability to %s...\n", prob_fname) ;
      printf("ERROR: saving posteriors does not work\n");
      exit(1);
      break ;
    case 'H':
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
usage_exit(int code)
{
  print_help();
  exit(code) ;
}

#include "mris_ca_label.help.xml.h"
static void
print_help(void)
{
  outputHelpXml(mris_ca_label_help_xml,
                mris_ca_label_help_xml_len);
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static int
postprocess(GCSA *gcsa, MRI_SURFACE *mris)
{
  LABEL       **larray, *area ;
  int         nlabels, i, j, annotation, n, nchanged, niter = 0, deleted ;
  double      max_area, label_area ;

#define MAX_ITER 5

  do
  {
    deleted = nchanged  = 0 ;
    MRISsegmentAnnotated(mris, &larray, &nlabels, 0) ;
    /*    printf("%d total segments in Gibbs annotation\n", nlabels) ;*/
    MRISclearMarks(mris) ;
    for (i = 0 ; i < nlabels ; i++)
    {
      area = larray[i] ;
      if (!area)   /* already processed */
      {
        continue ;
      }
      annotation = mris->vertices[area->lv[0].vno].annotation ;

      /* find label with this annotation with max area */
      max_area = LabelArea(area, mris) ;
      for (n = 1, j = i+1 ; j < nlabels ; j++)
      {
        if (!larray[j])
        {
          continue ;
        }
        if (annotation !=
            mris->vertices[larray[j]->lv[0].vno].annotation)
        {
          continue ;
        }
        n++ ;
        label_area = LabelArea(larray[j], mris) ;
        if (label_area > max_area)
        {
          max_area = label_area ;
        }
      }
#if 0
      printf("%03d: annotation %s (%d): %d segments, max area %2.1f\n",
             niter, annotation_to_name(annotation, NULL),
             annotation, n, max_area) ;
#endif
      for (j = i ; j < nlabels ; j++)
      {
        if (!larray[j])
        {
          continue ;
        }
        if (annotation !=
            mris->vertices[larray[j]->lv[0].vno].annotation)
        {
          continue ;
        }

        label_area = LabelArea(larray[j], mris) ;
        if (label_area < MIN_AREA_PCT*max_area)
        {
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
  }
  while (nchanged > 0 && niter++ < MAX_ITER) ;
  return(NO_ERROR) ;
}

#define MARK_RELABEL 2
#define MAX_EXCLUDED 100
static int
relabel_unknowns_with_cortex_label(GCSA *gcsa,
                                   MRI_SURFACE *mris,
                                   LABEL *cortex_label)
{
  int vno, n, annot;
  int nexcluded, exclude_list[MAX_EXCLUDED];
  int num_marked_for_relabel;
  VERTEX  *v ;

  printf("rationalizing unknown annotations with cortex label\n") ;

  nexcluded = 0 ;
  num_marked_for_relabel = 0;

  // Medial_wall label is in Christophe atlas
  annot = CTABentryNameToAnnotation("Medial_wall", mris->ct);
  if (annot >= 0)
  {
    exclude_list[nexcluded++] = annot ;
    printf("relabeling Medial_wall label...\n");
  }

  // 'unknown' label is essentially the medial wall in Desikan atlas
  annot = CTABentryNameToAnnotation("unknown", mris->ct);
  if (annot >= 0)
  {
    exclude_list[nexcluded++] = annot ;
    printf("relabeling unknown label...\n");
  }

  // corpus callosum labeled on medial wall in Desikan atlas
  annot = CTABentryNameToAnnotation("corpuscallosum", mris->ct);
  if (annot >= 0)
  {
    exclude_list[nexcluded++] = annot ;
    printf("relabeling corpuscallosum label...\n");
  }

  MRISclearMarks(mris) ;
  LabelMark(cortex_label, mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;

    v = &mris->vertices[vno] ;
    if (v->marked == 0)  // cortex label says it's not in cortex
    {
      if (vno == Gdiag_no)
	printf("vertex not in cortex.label - changing annotation to unknown\n") ;
      v->annotation = 0;  // replace with empty (transparent) annotation
    }
    else // cortex label says it is in cortex
    {
      if (v->annotation <= 0) // annotation is Unknown or invalid
      {
	if (vno == Gdiag_no)
	  printf("vertex label unknown but in cortex.label - marking for reclassification\n") ;
        v->marked = MARK_RELABEL;
        num_marked_for_relabel++;
      }
      for (n = 0 ; n < nexcluded ; n++)
        if (v->annotation == exclude_list[n])
        {
	  if (vno == Gdiag_no)
	    printf("vertex label in exclude list - marking for reclassification\n") ;
          v->marked = MARK_RELABEL ;
          num_marked_for_relabel++;
          break ;
        }
    }
  }
  printf("%d vertices marked for relabeling...\n", num_marked_for_relabel);
  GCSAreclassifyMarked(gcsa, mris, MARK_RELABEL, exclude_list, nexcluded) ;

  return(NO_ERROR) ;
}

