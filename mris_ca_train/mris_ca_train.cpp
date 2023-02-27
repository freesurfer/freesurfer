/**
 * @brief builds a cortical parcellation atlas (.gcs) file from a training set
 *
 * Creates a cortical parcellation atlas file based on one or more annotated
 * subjects. mris_ca_train builds probabilistic information estimated from a
 * manually labeled training set (of annotated subjects). Note that an
 * "annotation" is synonymous with a "parcellation", and is used for backwards
 * compatibility. The manual labeling can be carried out directly on surface
 * models using drawing tools in tksurfer, or volumetrically, then sampled
 * onto the surfaces using mris_sample_parc. This information is then used by
 * mris_ca_label to automatically assign a neuroanatomical label to each
 * location on a cortical surface model. This procedure incorporates both
 * geometric information derived from the cortical model (sulcus and
 * curvature), and neuroanatomical convention, as found in the training set.
 * The result of mris_ca_train and mris_ca_label is a complete labeling of
 * cortical sulci and gyri.
 *
 * See http://surfer.nmr.mgh.harvard.edu/fswiki/mris_ca_train
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
#include "version.h"
#include "label.h"


#define MAX_LABELS  1000
#if 0
static int write_ptable(char *fname, int *ptable, int nparcs) ;
#endif
static int find_parc_index(int parc, int *ptable, int nparcs) ;
static int add_to_ptable(MRI_SURFACE *mris, int *ptable, int nparcs) ;
static int *ptable = NULL ;
static int nbrs = 2 ;
static int navgs = 5 ;
static int normalize1_flag = 0 ;
static int normalize2_flag = 0 ;
static int normalize3_flag = 0 ;
static int nparcs = 0 ;
static char *ptable_fname = NULL ;
static COLOR_TABLE *ctab = NULL ;
static int which_norm = NORM_MEAN;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static const char *orig_name = "smoothwm" ;
static char *label_name = NULL ;
static int label_index ;

static int ninputs = 1 ;  /* curv and sulc */
static int icno_priors = 7 ;
static int icno_classifiers = 4 ;

#if 0
static char *curv_name = "curv" ;
#endif
static const char *thickness_name = "thickness" ;
static const char *sulc_name = "sulc" ;
static int sulconly = 0 ;

static char subjects_dir[STRLEN] ;
int nfillmax = -1;
int DoFill = 1;

int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp, *hemi;
  char         *canon_surf_name, *annot_name ;
  int          ac, nargs, i, train_type ;
  int          msec, minutes, seconds, nsubjects, input1_flags;
  int          input2_flags, input3_flags ;
  Timer start ;
  MRI_SURFACE  *mris ;
  GCSA         *gcsa ;
  int          unknown_index = -1 ;

  nargs = handleVersionOption(argc, argv, "mris_ca_train");
  if (nargs && argc - nargs == 1)
    exit (0);
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

  if (ptable_fname)
  {
    ctab = CTABreadASCII(ptable_fname) ;
    if (label_name)
    {
      CTABfindName(ctab, label_name, &label_index) ;
      if (label_index < 0)
        ErrorExit(ERROR_BADPARM, 
                  "%s: could not find index for label %s in table %s",
                  Progname, label_name, ptable_fname) ;

      CTABfindName(ctab, "unknown", &unknown_index) ;
      if (unknown_index < 0)
        CTABfindName(ctab, "medial wall", &unknown_index) ;
      printf("label index = %d, unknown index = %d\n", 
             label_index, unknown_index) ;
    }
  }
  else if (label_name)
    ErrorExit
      (ERROR_UNSUPPORTED, 
       "%s: must specify colortable with -t <ctab> when specifying label",
       Progname) ;

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

  hemi = argv[1] ;
  canon_surf_name = argv[2] ;
  annot_name = argv[3] ;
  out_fname = argv[argc-1] ;
  nsubjects = argc-5 ;

  gcsa = GCSAalloc(ninputs, icno_priors, icno_classifiers) ;
  input1_flags = input2_flags = input3_flags = 0 ;
  if (normalize1_flag)    input1_flags |= GCSA_NORMALIZE ;
  if (normalize2_flag)    input2_flags |= GCSA_NORMALIZE ;
  if (normalize3_flag)    input3_flags |= GCSA_NORMALIZE ;
  if(ctab) gcsa->ct = ctab;

  if (sulconly) { // not the default
    GCSAputInputType(gcsa,GCSA_INPUT_CURV_FILE,sulc_name,0,0,input1_flags);
  }
  else{ // default
    GCSAputInputType(gcsa, GCSA_INPUT_CURVATURE, "mean_curvature",navgs, input1_flags, 0) ;
    if(ninputs > 1)
      GCSAputInputType(gcsa,GCSA_INPUT_CURV_FILE,sulc_name,0,input2_flags,1);
    if (ninputs > 2)
      GCSAputInputType(gcsa,GCSA_INPUT_CURV_FILE,thickness_name,0,input3_flags,2);
  }

  for (train_type = 0 ; train_type <= 1 ; train_type++) {
    printf("computing %s for %d subject \n",train_type ? "covariances" : "means", nsubjects) ;
    for (i = 0 ; i < nsubjects ; i++)  {
      subject_name = argv[i+4] ;
      printf("processing subject %s, %d of %d...\n", subject_name,i+1,
             nsubjects);
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", subjects_dir, subject_name,
			 hemi, orig_name) ;    
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (DIAG_VERBOSE_ON)
        printf("reading surface from %s...\n", fname) ;
      mris = MRISread(fname) ;
      if (!mris) ErrorExit(ERROR_NOFILE,"%s: could not read surface file %s for %s",
			   Progname, fname, subject_name) ;
      MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
      if(label_name){ // not the default
        LABEL *area ;
        int   i ;
        VERTEX *v ;
        int req = snprintf(fname, STRLEN, "%s/%s/label/%s.%s",  subjects_dir, subject_name, hemi, annot_name) ; 
	if( req >= STRLEN )
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
        area = LabelRead(subject_name, fname) ;
        if(area == NULL) ErrorExit(ERROR_NOFILE,"%s: could not read label file %s for %s",Progname, fname, subject_name) ;
        for (i = 0 ; i < area->n_points ; i++) {
          if (area->lv[i].vno < 0)            continue ;
          v = &mris->vertices[area->lv[i].vno] ;
          if (area->lv[i].vno == Gdiag_no)
            DiagBreak() ;
          if (v->ripflag)continue ;
          CTABannotationAtIndex(ctab, label_index, &v->annotation) ;
        }
        // mark rest of surface unknown
        for (i = 0 ; i < mris->nvertices ; i++) {
          v = &mris->vertices[i] ;
          if (i == Gdiag_no)
            DiagBreak() ;
          if (v->ripflag || v->annotation) continue ;
          CTABannotationAtIndex(ctab, unknown_index, &v->annotation) ;
        }
        LabelFree(&area) ;
      }// end of label
      else{// use annotation (default)
	printf("Reading annotation %s\n",annot_name);
        if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
          ErrorExit(ERROR_NOFILE,"%s: could not read annot file %s for %s",Progname, annot_name, subject_name) ;
	//if(ctab) mris->ct = CTABdeepCopy(ctab);
	//CTABprintASCII(mris->ct, stdout);
      }
      if(ptable)
        nparcs = add_to_ptable(mris, ptable, nparcs) ;

      if (MRISreadCanonicalCoordinates(mris, canon_surf_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read spherical ""registration file %s for %s",
                  Progname, canon_surf_name, subject_name) ;

      if (ninputs > 2){ // not the default
        if (MRISreadCurvature(mris, thickness_name) != NO_ERROR)
          ErrorExit(ERROR_NOFILE,
                    "%s: could not read curv file %s for %s",
                    Progname, thickness_name, subject_name) ;
        if (normalize3_flag)
          MRISnormalizeCurvature(mris, which_norm) ;
        MRIScopyCurvatureToImagValues(mris) ;
      }

      if (ninputs > 1 || sulconly) { // not the default
        if (MRISreadCurvature(mris, sulc_name) != NO_ERROR)
          ErrorExit(ERROR_NOFILE,"%s: could not read curv file %s for %s",
                    Progname, sulc_name, subject_name) ;
        if(normalize2_flag || (sulconly && normalize1_flag))
          MRISnormalizeCurvature(mris, which_norm) ;
        MRIScopyCurvatureToValues(mris) ;
        MRIScopyValToVal2(mris) ;
      }
      if(!sulconly){ // default
        MRISuseMeanCurvature(mris) ;
        MRISaverageCurvatures(mris, navgs) ;
        if(normalize1_flag)  MRISnormalizeCurvature(mris, which_norm) ;
        MRIScopyCurvatureToValues(mris) ;
      }

      MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
      MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
      if (train_type == 0)
        GCSAtrainMeans(gcsa, mris) ;
      else
        GCSAtrainCovariances(gcsa, mris) ;
      MRISfree(&mris) ;
    }
    if (train_type == 0)
      GCSAnormalizeMeans(gcsa) ;
    else
      GCSAnormalizeCovariances(gcsa) ;
  }

  if(DoFill){
    if(nfillmax > 0) printf("Filling with a maximum of %d iterations\n",nfillmax);
    else             printf("Filling without a limit of  iterations\n");
    GCSAfill_cpn_holes(gcsa,nfillmax);
    GCSAfill_gcsan_holes(gcsa,nfillmax);
  }

  printf("writing classifier array to %s...\n", out_fname) ;
  gcsa->ptable_fname = ptable_fname ;
  GCSAwrite(gcsa, out_fname) ;
  GCSAfree(&gcsa) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("classifier array training took %d minutes"
         " and %d seconds.\n", minutes, seconds) ;
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

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help") || !stricmp(option, "help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  }
  else if (!stricmp(option, "ORIG"))
  {
    orig_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as original surface\n", orig_name) ;
  }
  else if (!stricmp(option, "NORM1"))
  {
    printf("normalizing input #1 after reading...\n") ;
    normalize1_flag = 1 ;
  }
  else if (!stricmp(option, "NORM2"))
  {
    printf("normalizing input #2 after reading...\n") ;
    normalize2_flag = 1 ;
  }
  else if (!stricmp(option, "NORM3"))
  {
    printf("normalizing input #3 after reading...\n") ;
    normalize3_flag = 1 ;
  }
  else if (!stricmp(option, "IC"))
  {
    icno_priors = atoi(argv[2]) ;
    icno_classifiers = atoi(argv[3]) ;
    nargs = 2 ;
    printf("using ico # %d for classifier array, and %d for priors\n",
           icno_classifiers, icno_priors) ;
  }
  else if (!stricmp(option, "SULC") || !stricmp(option, "SULCONLY"))
  {
    printf("using sulc as only input...\n") ;
    sulconly = 1 ;
  }
  else if (!stricmp(option, "gcs-means")){
    // gcsa inputno out.mgz
    GCSA *gcsa = GCSAread(argv[2]);
    if(!gcsa) exit(1);
    int inputno=0;
    sscanf(argv[3],"%d",&inputno);
    printf("inputno %d\n",inputno);
    MRI *mri = GCSAlikelihoodMeans2MRI(gcsa, inputno);
    if(mri==NULL) exit(1);
    printf("Writing to %s\n",argv[4]);
    int err = MRIwrite(mri,argv[4]);
    exit(err);
  }
  else if(!stricmp(option, "gcs-priors")){
    // gcsa out.mgz
    GCSA *gcsa = GCSAread(argv[2]);
    if(!gcsa) exit(1);
    MRI *mri = GCSApriors2MRI(gcsa);
    if(mri==NULL) exit(1);
    int err = MRIwrite(mri,argv[3]);
    exit(err);
  }
  else if(!strcmp(option, "debug-vertex")){
    Gdiag_no = atoi(argv[2]) ;
    printf("Gdiag_no set to %d\n",Gdiag_no);
    nargs = 1;
  }
  else if(!strcmp(option, "nfill")){
    nfillmax = atoi(argv[2]) ;
    printf("Setting fill iterations to %d\n",nfillmax);
    nargs = 1;
  }
  else if(!strcmp(option, "no-fill")){
    DoFill = 0;
    printf("Turning off fill\n");
  }
  else switch (toupper(*option))
    {
    case 'L':
      label_name = argv[2] ;
      nargs = 1 ;
      printf("interpreting inputs as label files for %s "
             "intead of annotations\n", label_name) ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'T':
      ptable_fname = argv[2] ;
      nargs = 1 ;
      ptable = (int *)calloc(MAX_LABELS, sizeof(int)) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'N':
      ninputs = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
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
print_usage(void)
{
  printf("mris_ca_train [options] <hemi> <canonsurf> <annot file> <subject 1> <subject 2> ... <output file>\n");
}

static void
usage_exit(int code)
{
  print_usage();
  exit(code) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "\n"
          "Creates a cortical parcellation atlas file based on one or \n"
          "more annotated subjects. mris_ca_train builds probabilistic \n"
          "information estimated from a manually labeled training set \n"
          "(of annotated subjects). The manual labeling can be carried out \n"
          "directly on surface models using drawing tools in tksurfer, or \n"
          "volumetrically, then sampled onto the surfaces using \n"
          "mris_sample_parc. This information is then used by mris_ca_label\n"
          "to automatically assign a neuroanatomical label to each location \n"
          "on a cortical surface model. This procedure incorporates both \n"
          "geometric information derived from the cortical model (sulcus \n"
          "and curvature), and neuroanatomical convention, as found in the \n"
          "training set. The result of mris_ca_train and mris_ca_label is \n"
          "a complete labeling of cortical sulci and gyri.\n\n");
  printf("Required args:\n") ;
  printf("  <hemi>               hemisphere: rh or lh\n");
  printf("  <canon surf>         canonical surface filename\n");
  printf("  <annot file>         annotation filename\n");
  printf("  <subject 1> <subject 2> ...  the subject id(s)\n");
  printf("  <outputfile>         classifier array output file\n");
  printf("Optional args:-------------\n");
  printf("  -sdir <directory>    use <directory> as subjects directory (default: $SUBJECTS_DIR)\n");
  printf("  -nbrs <number>       neighborhood size (default=2)\n");
  printf("  -orig <filename>     specify filename of original surface   (default=smoothwm)\n");
  printf("  -norm1               GCSA normalize input #1 after reading   (default: disabled)\n");
  printf("  -norm2               GCSA normalize input #2 after reading   (default: disabled)\n");
  printf("  -norm3               GCSA normalize input #3 after reading   (default: disabled)\n");
  printf("  -ic <number_priors> <number_classifiers>   parameters passed to the classifier routine (default: -ic 7 4)""\n");
  printf("  -sulc                 specify sulc as only input    (default: sulcus and curvature)\n");
  printf("  -sulconly             same as -sulc\n");
  printf("  -a <number>           number of averages (default=5)\n");
  printf("  -t <filename>         specify parcellation table input file  (default: none)\n");
  printf("   -n <number>           number of inputs (default=1)\n");
  printf("   -v <number>            diagnostic level (default=0)\n");
  printf("   -debug-vertex <number> diagnostic level (default=0)\n");
  printf("   -gcs-means gcsa inputno means.mgz : stand-alone to extract\n");
  printf("      likelihood means for all classes for given input\n");
  printf("   -gcs-priors gcsa priors.mgz : stand-alone to extract\n");
  printf("      priors for all classes for given input\n");
  printf("   -nfill nfill : set the max number of iterations for filling empty vertices\n");
  printf("   -no-fill : do not fill at all\n");
  printf("   --help                print help info\n");
  printf("   --version             print version info\n");
  exit(1) ;
}

static void
print_version(void)
{
  printf( "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

#if 0
static int
write_ptable(char *fname, int *ptable, int nparcs)
{
  FILE   *fp ;
  int    i ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "write_ptable(%s, %d): could not open file",
                 fname, nparcs)) ;

  for (i = 0 ; i < nparcs ; i++)
    fprintf(fp, "%d   %d\n", i, ptable[i]) ;
  fclose(fp) ;
  return(NO_ERROR) ;
}
#endif
static int
add_to_ptable(MRI_SURFACE *mris, int *ptable, int nparcs)
{
  int     vno, i ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    i = find_parc_index(v->annotation, ptable, nparcs) ;
    if (i < 0 || i >= nparcs)
    {
      i = nparcs++ ;
      ptable[i] = v->annotation ;
    }
  }
  return(nparcs) ;
}

static int
find_parc_index(int parc, int *ptable, int nparcs)
{
  int   i ;

  for (i = 0 ; i < nparcs ; i++)
    if (ptable[i] == parc)
      return(i) ;
  return(-1) ;
}

