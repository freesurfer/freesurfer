
/**
 * @brief Creates a Random Forest classifier for longitudinal data
 *
 * See:
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

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "gca.h"
#include "gcamorph.h"
#include "transform.h"
#include "cma.h"
#include "version.h"
#include "rforest.h"
#include "rfa.h"
#include "gca.h"
#include "talairachex.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include "romp_support.h"

#define MAX_RFA_INPUTS 1000
#define MAX_TIMEPOINTS 20

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

static char *log_file_name = NULL ;
static int conform = 1 ;
static int binarize = 0 ;
static int binarize_in = 0 ;
static int binarize_out = 0 ;
static char *wmsa_fname = NULL ;

static COLOR_TABLE *ctab = NULL ;
const char *Progname ;
static char *mask_fname = NULL ;
static char *insert_fname = NULL ;
static int  insert_label = 0 ;

static float scale = 0 ;
static int force_inputs = 1 ;

static RFA_PARMS parms ;
static const char *seg_dir = "seg_edited.mgz" ; // default name of manual edit file
static char T1_name[STRLEN] = "orig" ;
static char *xform_name = NULL;
static float smooth = -1 ;
static double TRs[MAX_RFA_INPUTS] ;
static double TEs[MAX_RFA_INPUTS] ;
static double FAs[MAX_RFA_INPUTS] ;

static int ninputs = 1 ;  /* T1 intensity */
static int navgs = 0 ;

static char subjects_dir[STRLEN] ;

static char *input_names[MAX_RFA_INPUTS] =
  {
    T1_name
  } ;

static int do_sanity_check = 0;
static int do_fix_badsubjs = 0;
static int sanity_check_badsubj_count = 0;

static int single_classifier_flag = 0 ;
static int only_nbrs ;   // only pick voxels that are on borders of a wmsa to train
static float max_wm_wmsa_ratio = 5.0 ;
static int make_uchar = 1 ;
static char *gca_name = NULL ;
static float wm_thresh = .8 ;  // only consider voxels with a prior at least this high
static int max_steps = 10 ;
static const char *single_classifier_names[] = 
{ "NOT WMSA", "WMSA", "FUTURE WMSA" } ;

#define NCLASSES     3
#define NOT_WMSA     0
#define WMSA         1
#define FUTURE_WMSA  2
#define MAX_SUBJECTS 1000

static MRI *mri_inputs[MAX_SUBJECTS][MAX_TIMEPOINTS] ;
static MRI *mri_segs[MAX_SUBJECTS][MAX_TIMEPOINTS] ;
static TRANSFORM *transforms[MAX_SUBJECTS][MAX_TIMEPOINTS] ;

static int lateralize_hypointensities(MRI *mri_seg) ;
static void usage_exit(int code) ;
static int replaceLabels(MRI *mri_seg) ;
static int check(MRI *mri_seg, char *subjects_dir, char *subject_name) ;
static RANDOM_FOREST *train_rforest(MRI *mri_inputs[MAX_SUBJECTS][MAX_TIMEPOINTS], MRI *mri_segs[MAX_SUBJECTS][MAX_TIMEPOINTS], TRANSFORM *transforms[MAX_SUBJECTS][MAX_TIMEPOINTS], 
				    int nsubjects, GCA *gca, RFA_PARMS *parms, float wm_thresh,
				    int wmsa_whalf, int ntp) ;
static int wmsa_whalf = 0 ;

int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp, *tp1_name, *tp2_name ;
  char         s1_name[STRLEN], s2_name[STRLEN], *sname ;
  int          ac, nargs, i, n, options, max_index ;
  int          msec, minutes, seconds, nsubjects, input ;
  Timer start ;
  MRI          *mri_seg, *mri_tmp, *mri_in ;
  TRANSFORM    *transform ;
//  int          counts ;
  int          t;
  RANDOM_FOREST *rf = NULL ;
  GCA           *gca = NULL ;

  Progname = argv[0] ;

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  parms.width = parms.height = parms.depth = DEFAULT_VOLUME_SIZE ;
  parms.ntrees = 10 ;
  parms.max_depth = 10 ;
  parms.wsize = 1 ;
  parms.training_size = 100 ;
  parms.training_fraction = .5 ;
  parms.feature_fraction = 1 ;

  nargs = handleVersionOption(argc, argv, "mri_rf_long_train");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  // parse command line args
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
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
  }
  if (argc < 3)
    usage_exit(1) ;


  // options parsed.   subjects, tp1 and tp2 and rf name remaining
  out_fname = argv[argc-1] ;
  nsubjects = (argc-2)/3 ;
  for (options = i = 0 ; i < nsubjects ; i++)
  {
    if (argv[i+1][0] == '-')
    {
      nsubjects-- ;
      options++ ;
    }
  }

  printf("training on %d subject and writing results to %s\n",
         nsubjects, out_fname) ;

  // rf_inputs can be T1, PD, ...per subject
  if (parms.nvols == 0)
    parms.nvols = ninputs ;
  /* gca reads same # of inputs as we read
     from command line - not the case if we are mapping to flash */
  n = 0 ;

  //////////////////////////////////////////////////////////////////
  // set up gca direction cosines, width, height, depth defaults

  gca = GCAread(gca_name) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read GCA from %s", Progname, gca_name) ;
  
  
  /////////////////////////////////////////////////////////////////////////
  // weird way options and subject name are mixed here
  
  /////////////////////////////////////////////////////////
  // first calculate mean
  ////////////////////////////////////////////////////////
  // going through the subject one at a time
  max_index = nsubjects+options ;
  nargs = 0 ;
  mri_in = NULL ; 
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  subject_name = NULL ; sname = NULL ; t = 0 ;
//  counts = 0 ;   would be private
  input = 0 ;
  transform = NULL ;
  tp1_name = tp2_name = NULL ;
  mri_tmp = mri_seg = NULL ;
#pragma omp parallel for if_ROMP(experimental) firstprivate(tp1_name, tp2_name, mri_in,mri_tmp, input, xform_name, transform, subjects_dir, force_inputs, conform, Progname, mri_seg, subject_name, s1_name, s2_name, sname, t, fname) shared(mri_inputs, transforms, mri_segs,argv) schedule(static,1)
#endif
  for (i = 0 ; i < max_index ; i++)
  {
    ROMP_PFLB_begin
    subject_name = argv[3*i+1] ;
    tp1_name = argv[3*i+2] ;
    tp2_name = argv[3*i+3] ;
    sprintf(s1_name, "%s_%s.long.%s_base", subject_name, tp1_name, subject_name) ;
    sprintf(s2_name, "%s_%s.long.%s_base", subject_name, tp2_name, subject_name) ;

    //////////////////////////////////////////////////////////////
    printf("***************************************"
	   "************************************\n");
    printf("processing subject %s, %d of %d (%s and %s)...\n", subject_name,i+1-nargs,
	   nsubjects, s1_name,s2_name);

    for (t = 0 ; t < 2 ; t++)
    {
      sname = t == 0 ? s1_name : s2_name;

      // reading this subject segmentation
      sprintf(fname, "%s/%s/mri/%s", subjects_dir, sname, seg_dir) ;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
	fprintf(stderr, "Reading segmentation from %s...\n", fname) ;
      mri_seg = MRIread(fname) ;
      if (!mri_seg)
	ErrorExit(ERROR_NOFILE, "%s: could not read segmentation file %s",
		  Progname, fname) ;

      if ((mri_seg->type != MRI_UCHAR) && (make_uchar != 0))
      {
	MRI *mri_tmp ;
	mri_tmp = MRIchangeType(mri_seg, MRI_UCHAR, 0, 1,1);
	MRIfree(&mri_seg) ;
	mri_seg = mri_tmp ;
      }

      if (wmsa_fname)
      {
	MRI *mri_wmsa ;
	sprintf(fname, "%s/%s/mri/%s", subjects_dir, sname, wmsa_fname) ;
	printf("reading WMSA labels from %s...\n", fname) ;
	mri_wmsa = MRIread(fname) ;
	if (mri_wmsa == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not read WMSA file %s", fname) ;
	MRIbinarize(mri_wmsa, mri_wmsa,  1, 0, WM_hypointensities) ;
	MRIcopyLabel(mri_wmsa, mri_seg, WM_hypointensities) ;
	lateralize_hypointensities(mri_seg) ;
	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON )
	{
	  char s[STRLEN] ;
	  sprintf(s, "%s/%s/mri/seg_%s",
		  subjects_dir, subject_name, wmsa_fname) ;
	  MRIwrite(mri_seg, s) ;
	}
      }
      if (binarize)
      {
	int j ;
	for (j = 0 ; j < 256 ; j++)
	{
	  if (j == binarize_in)
	    MRIreplaceValues(mri_seg, mri_seg, j, binarize_out) ;
	  else
	    MRIreplaceValues(mri_seg, mri_seg, j, 0) ;
	}
      }
      if (insert_fname)
      {
	MRI *mri_insert ;
	
	sprintf(fname, "%s/%s/mri/%s",
		subjects_dir, subject_name, insert_fname) ;
	mri_insert = MRIread(fname) ;
	if (mri_insert == NULL)
	  ErrorExit(ERROR_NOFILE,
		    "%s: could not read volume from %s for insertion",
		    Progname, insert_fname) ;
	
	MRIbinarize(mri_insert, mri_insert, 1, 0, insert_label) ;
	MRIcopyLabel(mri_insert, mri_seg, insert_label) ;
	MRIfree(&mri_insert) ;
      }
      
      replaceLabels(mri_seg) ;
      MRIeraseBorderPlanes(mri_seg, 1) ;

      ////////////////////////////////////////////////////////////
      if (DIAG_VERBOSE_ON)
	fprintf(stderr,
		"Gather all input volumes for the subject %s.\n",
		subject_name);
      // inputs must be coregistered
      // note that inputs are T1, PD, ... per subject (same TE, TR, FA)
      for (input = 0 ; input < ninputs ; input++)
      {
	//////////// set the gca type //////////////////////////////
	// is this T1/PD training?
	// how can we allow flash data training ???????
	// currently checks the TE, TR, FA to be the same for all inputs
	// thus we cannot allow flash data training.
	////////////////////////////////////////////////////////////
	
	sprintf(fname, "%s/%s/mri/%s", subjects_dir, sname,input_names[input]);
	if (DIAG_VERBOSE_ON)
	  printf("reading co-registered input from %s...\n", fname) ;
	fprintf(stderr, "   reading input %d: %s\n", input, fname);
	mri_tmp = MRIread(fname) ;
	if (!mri_tmp)
	  ErrorExit
	    (ERROR_NOFILE,
	     "%s: could not read image from file %s", Progname, fname) ;
	// input check 1
	if (getSliceDirection(mri_tmp) != MRI_CORONAL)
	{
	  ErrorExit
	    (ERROR_BADPARM,
	     "%s: must be in coronal direction, but it is not\n",
	     fname);
	}
	// input check 2
	if (conform &&
	    (mri_tmp->xsize != 1 || mri_tmp->ysize != 1 || mri_tmp->zsize != 1))
	{
	  ErrorExit
	    (ERROR_BADPARM,
	     "%s: must have 1mm voxel size, but have (%f, %f, %f)\n",
	     fname, mri_tmp->xsize, mri_tmp->ysize, mri_tmp->ysize);
	}
	// input check 3 is removed.  now we can handle c_(ras) != 0 case
	// input check 4
	if (i == 0)
	{
	  TRs[input] = mri_tmp->tr ;
	  FAs[input] = mri_tmp->flip_angle ;
	  TEs[input] = mri_tmp->te ;
	}
	else if ((force_inputs == 0) &&
		 (!FEQUAL(TRs[input],mri_tmp->tr) ||
		  !FEQUAL(FAs[input],mri_tmp->flip_angle) ||
		  !FEQUAL(TEs[input], mri_tmp->te)))
	  ErrorExit
	    (ERROR_BADPARM,
	     "%s: subject %s input volume %s: sequence parameters "
	     "(%2.1f, %2.1f, %2.1f)"
	     "don't match other inputs (%2.1f, %2.1f, %2.1f)",
	     Progname, subject_name, fname,
	     mri_tmp->tr, DEGREES(mri_tmp->flip_angle), mri_tmp->te,
	     TRs[input], DEGREES(FAs[input]), TEs[input]) ;
	// first time do the following
	if (input == 0)
	{
	  int nframes = ninputs ;
	  
	  ///////////////////////////////////////////////////////////
	  mri_in = MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth,
				    mri_tmp->type, nframes) ;
	  if (!mri_in)
	    ErrorExit
	      (ERROR_NOMEMORY,
	       "%s: could not allocate input volume %dx%dx%dx%d",
	       mri_tmp->width, mri_tmp->height, mri_tmp->depth,nframes) ;
	  MRIcopyHeader(mri_tmp, mri_in) ;
	}
	// -mask option ////////////////////////////////////////////
	if (mask_fname)
	{
	  MRI *mri_mask ;
	  
	  sprintf(fname, "%s/%s/mri/%s",
		  subjects_dir, subject_name, mask_fname);
	  printf("reading volume %s for masking...\n", fname) ;
	  mri_mask = MRIread(fname) ;
	  if (!mri_mask)
	  ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
		    Progname, fname) ;
	
	  MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
	  MRIfree(&mri_mask) ;
	}
	MRIcopyFrame(mri_tmp, mri_in, 0, input) ;
	MRIfree(&mri_tmp) ;

      }// end of inputs per subject
    
    
      /////////////////////////////////////////////////////////
      // xform_name is given, then we can use the consistent c_(r,a,s) for gca
      /////////////////////////////////////////////////////////
      if (xform_name)
      {
	// we read talairach.xfm which is a RAS-to-RAS
	sprintf(fname, "%s/%s/mri/transforms/%s", subjects_dir, sname, xform_name) ;
	if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
	  printf("INFO: reading transform file %s...\n", fname);
	if (!FileExists(fname))
	{
	  fprintf(stderr,"ERROR: cannot find transform file %s\n",fname);
	  exit(1);
	}
	transform = TransformRead(fname);
	if (!transform)
	  ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
		    Progname, fname);
	
//        modify_transform(transform, mri_in, gca);
	// Here we do 2 things
	// 1. modify gca direction cosines to
	// that of the transform destination (both linear and non-linear)
	// 2. if ras-to-ras transform,
      // then change it to vox-to-vox transform (linear case)
	
      // modify transform to store inverse also
	TransformInvert(transform, mri_in) ;
      }
      else
      {
//        GCAreinit(mri_in, gca);
	// just use the input value, since dst = src volume
	transform = TransformAlloc(LINEAR_VOXEL_TO_VOXEL, NULL) ;
      }
      
      /////////////////////////////////////////////////////////
      if (do_sanity_check)
      {
	// conduct a sanity check of particular labels, most importantly
	// hippocampus, that such labels do not exist in talairach coords
	// where they are known not to belong (indicating a bad manual edit)
	int errs = check(mri_seg, subjects_dir, subject_name);
	if (errs) 
	{
	  printf(
	    "ERROR: mri_ca_train: possible bad training data! subject:\n"
	    "\t%s/%s\n\n", subjects_dir, subject_name);
	  fflush(stdout) ;
	  sanity_check_badsubj_count++;
	}
      }
      
      mri_segs[i][t] = mri_seg ;
      mri_inputs[i][t] = mri_in ;
      transforms[i][t] = transform ;
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  rf = train_rforest(mri_inputs, mri_segs, transforms, nsubjects, gca, &parms, wm_thresh,wmsa_whalf, 2) ;
  printf("writing random forest to %s\n", out_fname) ;
  if (RFwrite(rf, out_fname) != NO_ERROR)
    ErrorExit
      (ERROR_BADFILE, "%s: could not write rf to %s", Progname, out_fname) ;
  
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("classifier array training took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  static int first_input = 1 ;
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "GRADIENT"))
  {
    parms.use_gradient = 1 ;
    ninputs += 3 ;  /* components of the gradient */
  }

  else if (!stricmp(option, "MAX_RATIO"))
  {
    max_wm_wmsa_ratio = atof(argv[2]) ;
    nargs = 1 ;
    printf("using %2.1f as max wm/wmsa ratio for number of training samples\n", max_wm_wmsa_ratio) ;
  }
  else if (!stricmp(option, "CONFORM"))
  {
    conform = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%sassuming input volumes are conformed\n", 
           conform ? "" : "NOT ") ;
  }
  else if (!stricmp(option, "MAKE_UCHAR"))
  {
    make_uchar = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%smaking input volumes UCHAR\n", make_uchar ? "" : "NOT ") ;
  }
  else if (!stricmp(option, "TRAINING_FRACTION"))
  {
    parms.training_fraction = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting training_fraction = %2.6f\n", parms.training_fraction) ;
  }
  else if (!stricmp(option, "WMSA"))
  {
    wmsa_fname = argv[2] ;
    nargs = 1 ;
    printf("reading white matter signal abnormalities from %s\n", wmsa_fname) ;
  }
  else if (!stricmp(option, "INPUT"))
  {
    if (first_input)
    {
      ninputs-- ;
      first_input = 0 ;
    }

    input_names[ninputs++] = argv[2] ;
    nargs = 1 ;
    printf("input[%d] = %s\n", ninputs-1, input_names[ninputs-1]) ;
  }
  else if (!stricmp(option, "BINARIZE"))
  {
    binarize = 1 ;
    binarize_in = atoi(argv[2]) ;
    binarize_out = atoi(argv[3]) ;
    nargs = 2 ;
    printf("binarizing segmentation values, setting input %d to output %d\n",
           binarize_in, binarize_out) ;
  }
  else if (!stricmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!stricmp(option, "DEBUG_NODE"))
  {
    Ggca_x = atoi(argv[2]) ;
    Ggca_y = atoi(argv[3]) ;
    Ggca_z = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Ggca_x,Ggca_y,Ggca_z) ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
  else if (!stricmp(option, "DEBUG_LABEL"))
  {
    Ggca_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging label %s (%d)\n", cma_label_to_name(Ggca_label),
           Ggca_label) ;
  }
  else if (!stricmp(option, "INSERT"))
  {
    insert_fname = argv[2] ;
    insert_label = atoi(argv[3]) ;
    nargs = 2 ;
    printf("inserting non-zero vals from %s as label %d...\n",
           insert_fname,insert_label);
  }
  else if (!stricmp(option, "ctab"))
  {
    printf("reading color table from %s and embedding in .gca file", argv[2]) ;
    ctab = CTABreadASCII(argv[2]) ;
    if (ctab == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read color table from %s", Progname,argv[2]);
    nargs = 1 ;
  }
  else if (!stricmp(option, "T1"))
  {
    strcpy(T1_name, argv[2]) ;
    nargs = 1 ;
    printf("reading T1 data from subject's mri/%s directory\n",
           T1_name) ;
  }
  else if (!stricmp(option, "PARC_DIR") || !stricmp(option, "SEG_DIR") ||
           !stricmp(option, "SEG") || !stricmp(option, "SEGMENTATION"))
  {
    seg_dir = argv[2] ;
    nargs = 1 ;
    printf("reading segmentation from subject's mri/%s directory\n",
           seg_dir) ;
  }
  else if (!stricmp(option, "XFORM"))
  {
    xform_name = argv[2] ;
    nargs = 1 ;
    printf("reading xform from %s\n", xform_name) ;
  }
  else if (!stricmp(option, "NOXFORM"))
  {
    xform_name = NULL ;
    printf("disabling application of xform...\n") ;
  }
  else if (!stricmp(option, "check"))
  {
    do_sanity_check = 1;
    printf("will conduct sanity-check of labels...\n") ;
  }
  else if (!stricmp(option, "check_and_fix"))
  {
    do_sanity_check = 1;
    do_fix_badsubjs = 1;
    printf("will conduct sanity-check of labels and write corrected "
           "volume to seg_fixed.mgz...\n") ;
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  }
  else if (!stricmp(option, "NTREES"))
  {
    parms.ntrees = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d trees in random forest classifier\n", parms.ntrees) ;
  }
  else if (!stricmp(option, "MAX_DEPTH"))
  {
    parms.max_depth = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d max tree depth in random forest classifier\n", parms.max_depth) ;
  }
  else if (!stricmp(option, "WMSA_WHALF"))
  {
    wmsa_whalf = atoi(argv[2]) ;
    nargs = 1 ;
    printf("only examing voxels that wmsa occurred within %d voxels of\n", wmsa_whalf) ;
  }
  else if (!stricmp(option, "SMOOTH"))
  {
    smooth = atof(argv[2]) ;
    if (smooth <= 0 || smooth > 1)
      ErrorExit(ERROR_BADPARM,
                "%s: smoothing parameter %2.1f must be in [0,1]\n",
                Progname, smooth) ;
    nargs = 1 ;
    printf("imposing %2.1f smoothing on conditional statistics\n", smooth) ;
  }
  else if (!stricmp(option, "NBRS"))
  {
    only_nbrs = 1 ;
    printf("training RF using voxels that neighbor both a wmsa and a non-wmsa\n") ;
  }
  else switch (toupper(*option))
    {
    case 'G':
      single_classifier_flag = 1 ;
      gca_name = argv[2] ;
      nargs = 1 ;
      printf("training a single classifier instead of an array using gca %s\n", gca_name) ;
      break ;
    case 'T':
      wm_thresh = atof(argv[2]) ;
      nargs = 1 ;
      printf("thresholding wm priors at %f to build training set\n", wm_thresh) ;
      break ;
    case 'F':
      force_inputs = 1 ;
      printf("forcing use of inputs even if acquisition parameters don't match\n");
      break ;
    case 'S':
      scale = atof(argv[2]) ;
      printf("scaling all volumes by %2.3f after reading...\n", scale) ;
      nargs = 1 ;
      break ;
    case 'L':
      log_file_name = argv[2] ;
      printf("logging out of bag accuracy to %s\n", log_file_name) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      printf("applying %d mean filters to classifiers after training\n",navgs);
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    case 'W':
      parms.wsize = atoi(argv[2]) ;
      nargs = 1 ;
      printf("using window size = %d for RFA\n", parms.wsize) ;
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
  printf("Purpose: %s trains GCA data with (multiple) subject(s)\n", Progname);
  printf("Usage  : %s [options] <subject1> <subject2> ... "
         "<output rfa fname>\n",
         Progname) ;
  printf
  ("where SUBJECTS_DIR env variable must be set.\n"
   "Options are:\n"
   "\t-seg dir        - (required) segmentation volume "
   "(path relative to $subject/mri).\n"
   "\t-xform xform    -  atlas transform (path relative "
   "to $subject/mri/transforms).\n"
   "\t-mask volname   - use volname as a mask "
   "(path relative to $subject/mri.\n"
   "\t-node_spacing   - spacing of classifiers in canonical space\n"
   "\t-prior_spacing  - spacing of class priors in canonical space\n"
   "\t-input name     - specifying training data "
   "(path relative to $subject/mri).\n"
   "                          can specify multiple inputs.  "
   "If not specified, \"orig\" is used\n"
   "\t-check          - conduct sanity-check of labels for obvious edit errors"
   "\n"
  );
  exit(code) ;
}

static int input_labels[] =
  {
    Left_Cerebral_Exterior,
    Right_Cerebral_Exterior,
    Left_Cerebellum_Exterior,
    Right_Cerebellum_Exterior
  } ;
// replace the values above with the following
static int output_labels[] =
  {
    Left_Cerebral_Cortex,
    Right_Cerebral_Cortex,
    Left_Cerebellum_Cortex,
    Right_Cerebellum_Cortex
  } ;

static int
replaceLabels(MRI *mri_seg)
{
  int    i ;

  for (i = 0 ; i < sizeof(output_labels)/sizeof(output_labels[0]) ; i++)
    MRIreplaceValues(mri_seg, mri_seg, input_labels[i], output_labels[i]) ;
  return(NO_ERROR) ;
}

#if 0
// modify transform to vox-to-vox
static void modify_transform(TRANSFORM *transform, MRI *mri_in, GCA *gca)
{
  LTA *lta=0;
  MATRIX *i_to_r=0, *r_to_i=0, *tmpmat=0, *vox2vox;
  MRI *mri_buf = 0;
  GCA_MORPH *gcam = 0;
  static int warned = 0;

  // temp buf to get the transform
  mri_buf = MRIallocHeader(mri_in->width,
                           mri_in->height,
                           mri_in->depth,
                           mri_in->type,1);
  MRIcopyHeader(mri_in, mri_buf);

  //////////////////////////////////////////////////////////////////////////
  // non-linear transform case
  //////////////////////////////////////////////////////////////////////////
  if (transform->type == MORPH_3D_TYPE)
  {
    gcam = (GCA_MORPH *) transform->xform;
    if (gcam->atlas.valid) // means it contains the dst volume information
    {
      mri_buf->c_r = gcam->atlas.c_r;
      mri_buf->c_a = gcam->atlas.c_a;
      mri_buf->c_s = gcam->atlas.c_s;
      if (warned == 0)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf
          (stderr,
           "INFO: modified c_(r,a,s) using the non-linear "
           "transform dst value.\n");
        warned = 1;
      }
    }
    else // this is an old 3d, I should use c_(ras) = 0
    {
      mri_buf->c_r = 0;
      mri_buf->c_a = 0;
      mri_buf->c_s = 0;
      if (warned == 0)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stderr, "INFO: modified c_(r,a,s) = 0.\n");
        warned = 1;
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////
  /// linear transform case
  ////////////////////////////////////////////////////////////////////////////
  else if (transform->type == LINEAR_VOX_TO_VOX)
  {
    lta = (LTA *) (transform->xform);
    // modify using the xform dst
    if (lta->xforms[0].dst.valid)
    {
      mri_buf->c_r = lta->xforms[0].dst.c_r;
      mri_buf->c_a = lta->xforms[0].dst.c_a;
      mri_buf->c_s = lta->xforms[0].dst.c_s;
      if (warned == 0)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stderr, "INFO: modified c_(r,a,s) using the xform dst.\n");
        warned = 1;
      }
    }
    else // keep the old behavior
    {
      mri_buf->c_r = 0;
      mri_buf->c_a = 0;
      mri_buf->c_s = 0;
      if (warned == 0)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stderr, "INFO: modified c_(r,a,s) = 0.\n");
        warned = 1;
      }
    }
  }
  else if (transform->type == LINEAR_RAS_TO_RAS)
  {
    lta = (LTA *) (transform->xform);
    // modify using the xform dst
    if (lta->xforms[0].dst.valid)
    {
      mri_buf->c_r = lta->xforms[0].dst.c_r;
      mri_buf->c_a = lta->xforms[0].dst.c_a;
      mri_buf->c_s = lta->xforms[0].dst.c_s;
      if (warned == 0)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stderr, "INFO: modified c_(r,a,s) using the xform dst\n");
        warned = 1;
      }
    }
    // dst invalid
    else if (getenv("USE_AVERAGE305"))// use average_305 value
      // (usually ras-to-ras comes from MNI transform)
    {
      mri_buf->c_r = -0.095;
      mri_buf->c_a = -16.51;
      mri_buf->c_s =   9.75;
      if (warned == 0)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        {
          fprintf
          (stderr,
           "INFO: modified c_(r,a,s) using average_305 value\n");
          fprintf
          (stderr,
           "INFO: if this is not preferred, set environment "
           "variable NO_AVERAGE305\n");
        }
        warned = 1;
      }
    }
    else // keep old behavior
    {
      mri_buf->c_r = 0;
      mri_buf->c_a = 0;
      mri_buf->c_s = 0;
      if (warned == 0)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf
          (stderr, "INFO: xform.dst invalid thus modified c_(r,a,s) = 0.\n");
        warned = 1;
      }
    }
    /////////////////////////////////////////////////////////////////
    printf("INFO: original RAS-to-RAS transform\n");
    MatrixPrint(stdout, lta->xforms[0].m_L);
    // going from vox->RAS->TalRAS
    i_to_r = extract_i_to_r(mri_in);
    tmpmat = MatrixMultiply(lta->xforms[0].m_L, i_to_r, NULL);
    r_to_i = extract_r_to_i(mri_buf);
    // going from TalRAS -> voxel
    vox2vox = MatrixMultiply(r_to_i, tmpmat,NULL );
    printf("INFO: modified VOX-to-VOX transform\n");
    MatrixPrint(stdout, vox2vox);
    // store it
    MatrixCopy(vox2vox, lta->xforms[0].m_L);
    // now mark it as vox-to-vox
    transform->type = LINEAR_VOX_TO_VOX;
    // free up memory
    MatrixFree(&r_to_i);
    MatrixFree(&i_to_r);
    MatrixFree(&tmpmat);
    MatrixFree(&vox2vox);
  }
  /////////////////////////////////////////////////////////////////
  // OK now we know what the target c_(ras) should be
  // we reset c_(ras) value for GCA node and priors
  GCAreinit(mri_buf, gca);

  MRIfree(&mri_buf);
}
#endif

#define WSIZE 9
#define WHALF ((WSIZE-1)/2)
static int
lateralize_hypointensities(MRI *mri_seg)
{
  int left, right, n, x, y, z, label_counts[MAX_CMA_LABELS], label ;

  for (x = 0 ; x < mri_seg->width ; x++)
  {
    for (y = 0 ; y < mri_seg->height ; y++)
    {
      for (z = 0 ; z < mri_seg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        if (label != WM_hypointensities)
          continue ;
        MRIcomputeLabelNbhd
        (mri_seg, NULL, x, y, z,
         label_counts,  NULL, WHALF, MAX_CMA_LABELS) ;
        for (left = right = n = 0 ; n < MAX_CMA_LABELS ; n++)
        {
          switch (n)
          {
          case Left_Lateral_Ventricle:
          case Left_Cerebral_White_Matter:
          case Left_Caudate:
          case Left_Putamen:
          case Left_Pallidum:
          case Left_Thalamus:
          case Left_Cerebral_Cortex:
          case Left_Hippocampus:
          case Left_Amygdala:
          case Left_Cerebellum_Cortex:
          case Left_Cerebellum_White_Matter:
            left += label_counts[n] ;
            break ;
          case Right_Lateral_Ventricle:
          case Right_Cerebral_White_Matter:
          case Right_Caudate:
          case Right_Putamen:
          case Right_Pallidum:
          case Right_Thalamus:
          case Right_Cerebral_Cortex:
          case Right_Hippocampus:
          case Right_Amygdala:
          case Right_Cerebellum_Cortex:
          case Right_Cerebellum_White_Matter:
            right += label_counts[n] ;
            break ;
          }
        }
        if (left > right)
          MRIsetVoxVal(mri_seg, x, y, z, 0, Left_WM_hypointensities) ;
        else
          MRIsetVoxVal(mri_seg, x, y, z, 0, Right_WM_hypointensities) ;
      }
    }
  }

  return(NO_ERROR) ;
}


/*
 * check
 *
 * conduct a sanity check of particular labels, most importantly
 * hippocampus, that such labels do not exist in talairach coords
 * where they are known not to belong (indicating a bad manual edit)
 */
static int check(MRI *mri_seg, char *subjects_dir, char *subject_name)
{
  MRI *mri_fixed = NULL;
  int errors=0;
  int x, y, z, label=0;
  double xw=0.0, yw=0.0, zw=0.0; // RAS coords
  double xmt=0.0, ymt=0.0, zmt=0.0; // MNI tal coords
  float xt=0.0, yt=0.0, zt=0.0; // 'real' tal coords

  float max_xtal_l_hippo    = -1000;
  float max_xtal_l_caudate  = -1000;
  float max_xtal_l_amygdala = -1000;
  float max_xtal_l_putamen  = -1000;
  float max_xtal_l_pallidum = -1000;
  float min_xtal_r_hippo    =  1000;
  float min_xtal_r_caudate  =  1000;
  float min_xtal_r_amygdala =  1000;
  float min_xtal_r_putamen  =  1000;
  float min_xtal_r_pallidum =  1000;

  printf("checking labels in subject %s...\n",subject_name);

  if (do_fix_badsubjs)
  {
    mri_fixed = MRIcopy(mri_seg,NULL);
  }

  if (NULL == mri_seg->linear_transform)
  {
    ErrorExit(ERROR_BADFILE,
              "ERROR: mri_ca_train: talairach.xfm not found in %s!\n"
              "Run mri_add_xform_to_header to add to volume.\n",
              seg_dir);
  }

  // now conduct checks for voxels in locations where they shouldnt be
  for (x = 0 ; x < mri_seg->width ; x++)
  {
    for (y = 0 ; y < mri_seg->height ; y++)
    {
      for (z = 0 ; z < mri_seg->depth ; z++)
      {
        /*
         * rules:
         * - no label should have a voxel coord < 6 or > 249
         * - no left or right hippo labels with z tal coord > 15
         * - no left hippo, caudate, amydala, putamen or pallidum 
         *   labels with x tal coord > 5 (or right, with x tal coord < -5)
         */
        label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
        int proper_label = label; // used for making corrections

        if (label != Unknown)
        {
          // note: these won't be caught if -mask brainmask.mgz is included
          // on the command-line, since voxels outside of brainmask get 
          // labelled as Unknown.
          if ((x < 6) || (y < 6) || (z < 6) ||
              (x > 249) || (y > 249) || (z > 249))
          {
            printf
              ("ERROR: %s: found label outside of brain: "
               "%d %d %d\n", 
               cma_label_to_name(label),x,y,z);
            fflush(stdout) ;
            errors++;

            proper_label = Unknown;
          }
        }

        if ((label == Left_Hippocampus)  || (label == Right_Hippocampus) ||
            (label == Left_Caudate)      || (label == Right_Caudate) ||
            (label == Left_Amygdala)     || (label == Right_Amygdala) ||
            (label == Left_Putamen)      || (label == Right_Putamen) ||
            (label == Left_Pallidum)     || (label == Right_Pallidum) ||
            (label == Left_Inf_Lat_Vent) || (label == Right_Inf_Lat_Vent))
        {
          // the 'if' statement above spares some cpu cycles in having to
          // calculate the coord transform for every voxel, ie these 3 lines:
          MRIvoxelToWorld(mri_seg, x, y, z, &xw, &yw, &zw) ;
          transform_point(mri_seg->linear_transform, 
                          xw, yw, zw, &xmt, &ymt, &zmt);// get mni tal coords
          FixMNITal(xmt,ymt,zmt, &xt,&yt,&zt); // get 'real' tal coords

          switch (label)
          {
          case Left_Hippocampus:
            if (zt > 15)
            {
              printf
                ("ERROR: %s: "
                 "%d %d %d, tal x=%f, y=%f, *** z=%f > 15 ***\n", 
                 cma_label_to_name(label),x,y,z,xt,yt,zt);
              fflush(stdout) ;
              errors++;
            }
            // no break (check xt)

          case Left_Caudate:
          case Left_Amygdala:
          case Left_Putamen:
          case Left_Pallidum:
          case Left_Inf_Lat_Vent:
            if (xt > 5)
            {
              printf
                ("ERROR: %s: "
                 "%d %d %d, tal *** x=%f > 5 ***, y=%f, z=%f\n", 
                 cma_label_to_name(label), x,y,z,xt,yt,zt);
              fflush(stdout) ;
              errors++;

              if  (label == Left_Hippocampus) proper_label = Right_Hippocampus;
              else if (label == Left_Caudate)  proper_label = Right_Caudate;
              else if (label == Left_Amygdala) proper_label = Right_Amygdala;
              else if (label == Left_Putamen)  proper_label = Right_Putamen;
              else if (label == Left_Pallidum) proper_label = Right_Pallidum;
              else if (label == Left_Inf_Lat_Vent) 
                proper_label = Right_Inf_Lat_Vent;
            }
            break;

          case Right_Hippocampus:
            if (zt > 15)
            {
              printf
                ("ERROR: %s: "
                 "%d %d %d, tal x=%f, y=%f, *** z=%f > 15 ***\n", 
                 cma_label_to_name(label),x,y,z,xt,yt,zt);
              fflush(stdout) ;
              errors++;
            }
            // no break (check xt)

          case Right_Caudate:
          case Right_Amygdala:
          case Right_Putamen:
          case Right_Pallidum:
          case Right_Inf_Lat_Vent:
            if (xt < -5)
            {
              printf
                ("ERROR: %s: "
                 "%d %d %d, tal *** x=%f < -5 ***, y=%f, z=%f\n", 
                 cma_label_to_name(label), x,y,z,xt,yt,zt);
              fflush(stdout) ;
              errors++;

              if  (label == Right_Hippocampus) proper_label = Left_Hippocampus;
              else if (label == Right_Caudate)  proper_label = Left_Caudate;
              else if (label == Right_Amygdala) proper_label = Left_Amygdala;
              else if (label == Right_Putamen)  proper_label = Left_Putamen;
              else if (label == Right_Pallidum) proper_label = Left_Pallidum;
              else if (label == Right_Inf_Lat_Vent) 
                proper_label = Left_Inf_Lat_Vent;
            }
            break;

          default:
            break ;
          }

          /*
           * collect stats on positioning of structures.
           * used to determine reasonable boundaries.
           */
          switch (label)
          {
          case Left_Hippocampus:
            if (xt > max_xtal_l_hippo)    max_xtal_l_hippo = xt;
            break;
          case Left_Caudate:
            if (xt > max_xtal_l_caudate)  max_xtal_l_caudate = xt;
            break;
          case Left_Amygdala:
            if (xt > max_xtal_l_amygdala) max_xtal_l_amygdala = xt;
            break;
          case Left_Putamen:
            if (xt > max_xtal_l_putamen)  max_xtal_l_putamen = xt;
            break;
          case Left_Pallidum:
            if (xt > max_xtal_l_pallidum) max_xtal_l_pallidum = xt;
            break;

          case Right_Hippocampus:
            if (xt < min_xtal_r_hippo)    min_xtal_r_hippo = xt;
            break;
          case Right_Caudate:
            if (xt < min_xtal_r_caudate)  min_xtal_r_caudate = xt;
            break;
          case Right_Amygdala:
            if (xt < min_xtal_r_amygdala) min_xtal_r_amygdala = xt;
            break;
          case Right_Putamen:
            if (xt < min_xtal_r_putamen)  min_xtal_r_putamen = xt;
            break;
          case Right_Pallidum:
            if (xt < min_xtal_r_pallidum) min_xtal_r_pallidum = xt;
            break;

          default:
            break ;
          }
        }

        /* 
         * if -check_and_fix is being used, then mod our fixed volume
         */
        if (do_fix_badsubjs && (label != proper_label))
        {
          MRIsetVoxVal(mri_fixed, x, y, z, 0, proper_label) ;
        }
      }
    }
  }

  // stats used to determine optimal boundaries
  printf("max_xtal_l_hippo    = %4.1f\n",max_xtal_l_hippo);
  printf("max_xtal_l_caudate  = %4.1f\n",max_xtal_l_caudate);
  printf("max_xtal_l_amygdala = %4.1f\n",max_xtal_l_amygdala);
  printf("max_xtal_l_putamen  = %4.1f\n",max_xtal_l_putamen);
  printf("max_xtal_l_pallidum = %4.1f\n",max_xtal_l_pallidum);
  
  printf("min_xtal_r_hippo    = %4.1f\n",min_xtal_r_hippo);
  printf("min_xtal_r_caudate  = %4.1f\n",min_xtal_r_caudate);
  printf("min_xtal_r_amygdala = %4.1f\n",min_xtal_r_amygdala);
  printf("min_xtal_r_putamen  = %4.1f\n",min_xtal_r_putamen);
  printf("min_xtal_r_pallidum = %4.1f\n",min_xtal_r_pallidum);
  
  if ( do_fix_badsubjs && errors)
  {
    char fname[STRLEN];
    sprintf(fname, "%s/%s/mri/seg_fixed.mgz", subjects_dir, subject_name);
    printf("Writing corrected volume to %s\n",fname);
    MRIwrite(mri_fixed,fname);
    MRIfree(&mri_fixed);
  }

  fflush(stdout);

  return(errors) ;
}

static RANDOM_FOREST *
train_rforest(MRI *mri_inputs[MAX_SUBJECTS][MAX_TIMEPOINTS], MRI *mri_segs[MAX_SUBJECTS][MAX_TIMEPOINTS], TRANSFORM *transforms[MAX_SUBJECTS][MAX_TIMEPOINTS], 
	      int nsubjects, GCA *gca, RFA_PARMS *parms, float wm_thresh,
	      int wmsa_whalf, int ntp) 
{
  RANDOM_FOREST  *rf ;
  int            nfeatures, x, y, z, ntraining, n, tvoxel_size, width, height, depth, xt, yt, zt ;
  double         xatlas, yatlas, zatlas ;
  MRI            *mri_in, *mri_seg, *mri_training_voxels, *mri_wmsa_possible ;
  TRANSFORM      *transform ;
  double         **training_data ;
  int            *training_classes, i, label, tlabel, nwmsa, nfuture, nnot, correct, label_time1, label_time2;

  nwmsa = nnot = nfuture = 0 ;

/*
  features are:
  t1 intensity (3 vols)
  3 priors
  # of unknown voxels in the nbhd
  # of neighboring wmsa voxels at t1
*/
  nfeatures = parms->wsize*parms->wsize*parms->wsize*parms->nvols + 5 ; 

  rf = RFalloc(parms->ntrees, nfeatures, NCLASSES, parms->max_depth, const_cast<char**>(single_classifier_names), max_steps) ;
  if (rf == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not allocate random forest", Progname) ;
  rf->min_step_size = 1 ; 

  tvoxel_size=1 ;
  width = (int)ceil((float)mri_segs[0][0]->width/tvoxel_size) ;
  height = (int)ceil((float)mri_segs[0][0]->height/tvoxel_size) ;
  depth = (int)ceil((float)mri_segs[0][0]->depth/tvoxel_size) ;
  mri_wmsa_possible = MRIalloc(width, height, depth, MRI_UCHAR) ;
  GCAcopyDCToMRI(gca, mri_wmsa_possible) ;
  mri_in = mri_inputs[0][0] ;
  mri_training_voxels = MRIallocSequence(mri_in->width,mri_in->height, mri_in->depth,MRI_UCHAR,nsubjects) ;

#if 1
  // update time 1 segmentation based on labels at time1 and time2
  for (n = 0 ; n < nsubjects ; n++)
  {
    mri_seg = mri_segs[n][0] ;
    for (x = 0 ; x < mri_in->width ; x++)
      for (y = 0 ; y < mri_in->height ; y++)
	for (z = 0 ; z < mri_in->depth ; z++)
	{
	  label_time1 = MRIgetVoxVal(mri_segs[n][0], x, y, z, 0) ;
	  label_time2 = MRIgetVoxVal(mri_segs[n][1], x, y, z, 0) ;
	  if (IS_WMSA(label_time1))
	    MRIsetVoxVal(mri_segs[n][0], x, y, z, 0, label_time1) ;
	  else if (IS_WMSA(label_time2))
	    MRIsetVoxVal(mri_segs[n][0], x, y, z, 0, future_WMSA) ;
	}
  }
#endif

  // build map of spatial locations that WMSAs can possibly occur in
  for (n = 0 ; n < nsubjects ; n++)
  {
    mri_in = mri_inputs[n][1] ; transform = transforms[n][1] ; 
    for (x = 0 ; x < mri_in->width ; x++)
      for (y = 0 ; y < mri_in->height ; y++)
	for (z = 0 ; z < mri_in->depth ; z++)
	  if (is_possible_wmsa(gca, mri_in, transform, x, y, z, 0))
	  {
	    TransformSourceVoxelToAtlas(transform, mri_in, x, y, z, &xatlas, &yatlas, &zatlas) ;
	    xt = nint(xatlas/tvoxel_size) ;
	    yt = nint(yatlas/tvoxel_size) ;
	    zt = nint(zatlas/tvoxel_size) ;
	    if (xt == Gx && yt == Gy && zt == Gz)
	      DiagBreak() ;
	    MRIsetVoxVal(mri_wmsa_possible, xt, yt, zt, 0, 1) ;
	  }
  }
  for ( ; wmsa_whalf > 0 ; wmsa_whalf--)
    MRIdilate(mri_wmsa_possible, mri_wmsa_possible) ;

  // now build map of all voxels in training set
  for (nnot = nwmsa = nfuture = ntraining = n = 0 ; n < nsubjects ; n++)
  {
    mri_in = mri_inputs[n][0] ; mri_seg = mri_segs[n][0] ; transform = transforms[n][0] ;
    for (x = 0 ; x < mri_in->width ; x++)
      for (y = 0 ; y <  mri_in->height; y++)
	for (z = 0 ; z <  mri_in->depth; z++)
	{
	  label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
	  
	  TransformSourceVoxelToAtlas(transform, mri_in, x, y, z, &xatlas, &yatlas, &zatlas) ;
	  xt = nint(xatlas/tvoxel_size) ;
	  yt = nint(yatlas/tvoxel_size) ;
	  zt = nint(zatlas/tvoxel_size) ;
	  if (xt == Gx && yt == Gy && zt == Gz)
	    DiagBreak() ;
	  if ((IS_WMSA(label) == 0) &&
	      MRIgetVoxVal(mri_wmsa_possible,xt,yt, zt,0) == 0)
	    continue ;
	  if (NOT_TRAINING_LABEL(label))
	    continue ;
	  ntraining++ ;
	  
	  if (IS_FUTURE_WMSA(label))
	  {
	    label = FUTURE_WMSA;
	    nfuture++ ;
	  }
	  else if (IS_WMSA(label))
	  {
	    label = WMSA ;
	    nwmsa++ ;
	  }
	  else
	  {
	    label = NOT_WMSA ;
	    nnot++ ;
	  }
          // set label to one more than it will be for training so that 0 means this is not a training voxel
	  MRIsetVoxVal(mri_training_voxels, x, y, z, n, label+1) ;
	}
  }

  correct = MRIcountNonzero(mri_training_voxels) ;
  if (correct != ntraining)
    DiagBreak() ;
  printf("total training set size = %2.1fM\n", (float)ntraining/(1024.0f*1024.0f)) ;
  printf("initial training set found with %dK FUTURE WMSA labels, %dK WMSA labels, and %dK non (ratio=%2.1f:%2.1f)\n",
	 nfuture/1000, nwmsa/1000, nnot/1000, (float)nnot/(float)nfuture, (float)nnot/(float)nwmsa) ;

  MRIfree(&mri_wmsa_possible) ;

  if (max_wm_wmsa_ratio*(nfuture+nwmsa) < nnot)   // too many wm labels w.r.t. # of wmsas - remove some wm
  {
    int removed, total_to_remove =  nnot - (max_wm_wmsa_ratio*(nfuture+nwmsa)) ;
    double premove ;

    premove = (double)total_to_remove  / (double)nnot ;
    printf("removing %dK WM indices to reduce training set imbalance (p < %f)\n", total_to_remove/1000, premove) ;

    for (removed = n = 0 ; n < nsubjects ; n++)
      for (x = 0 ; x < mri_in->width ; x++)
	for (y = 0 ; y <  mri_in->height; y++)
	  for (z = 0 ; z <  mri_in->depth; z++)
	  {
	    label = MRIgetVoxVal(mri_training_voxels, x, y, z, n) ;
	    if (label == 1)   // a WM voxel
	    {
	      if (randomNumber(0,1) < premove)
	      {
		removed++ ;
		MRIsetVoxVal(mri_training_voxels, x, y, z, n, 0) ;  // remove it from training set
	      }
	    }
	  }
    ntraining -= removed ;
    printf("%d WM voxels removed, new training set size = %dM (ratio = %2.1f)\n",
	   removed, ntraining/(1024*1024), (double)(nnot-removed)/(double)(nwmsa+nfuture)) ;
  }

  correct = MRIcountNonzero(mri_training_voxels) ;
  if (correct != ntraining)
    DiagBreak() ;
//  if (Gx >= 0)
  {
    int whalf = (parms->wsize-1)/2 ;
    char buf[STRLEN] ;
    rf->feature_names = (char **)calloc(rf->nfeatures, sizeof(char *)) ;
    for (i = 0, x = -whalf ; x <= whalf ; x++)
      for (y = -whalf ; y <= whalf ; y++)
	for (z = -whalf ; z <= whalf ; z++)
	  for (n = 0 ; n < mri_in->nframes ; n++, i++)
	  {
	    switch (n)
	    {
	    default:
	    case 0: sprintf(buf, "T1(%d, %d, %d)", x, y, z) ; break ;
	    case 1: sprintf(buf, "T2(%d, %d, %d)", x, y, z) ; 
	      break ;
	    case 2: sprintf(buf, "FLAIR(%d, %d, %d)", x, y, z) ; 
	      if (x == 0 && y == 0 && z == 0)
		Gdiag_no = i ;
	      break ;
	    }
	    rf->feature_names[i] = (char *)calloc(strlen(buf)+1, sizeof(char)) ;
	    strcpy(rf->feature_names[i], buf) ;
	  }
    printf("FLAIR(0,0,0) = %dth feature\n", Gdiag_no) ;

    sprintf(buf, "CSF voxels in nbhd") ;
    rf->feature_names[i] = (char *)calloc(strlen(buf)+1, sizeof(char)) ;
    strcpy(rf->feature_names[i], buf) ;
    i++ ; sprintf(buf, "gm prior") ;
    rf->feature_names[i] = (char *)calloc(strlen(buf)+1, sizeof(char)) ;
    strcpy(rf->feature_names[i], buf) ;
    i++ ; sprintf(buf, "wm prior") ;
    rf->feature_names[i] = (char *)calloc(strlen(buf)+1, sizeof(char)) ;
    strcpy(rf->feature_names[i], buf) ;
    i++ ; sprintf(buf, "csf prior") ;
    rf->feature_names[i] = (char *)calloc(strlen(buf)+1, sizeof(char)) ;
    strcpy(rf->feature_names[i], buf) ;
    i++ ; sprintf(buf, "WMSA in nbhd") ;
    rf->feature_names[i] = (char *)calloc(strlen(buf)+1, sizeof(char)) ;
    strcpy(rf->feature_names[i], buf) ;

    if (Gdiag & DIAG_WRITE)
    {
      printf("writing training voxels to tv.mgz\n") ;
      MRIwrite(mri_training_voxels, "tv.mgz") ;
    }
  }

  // now build training features and classes
  training_classes = (int *)calloc(ntraining, sizeof(training_classes[0])) ;
  if (training_classes == NULL)
    ErrorExit(ERROR_NOFILE, "train_rforest: could not allocate %d-length training buffers",ntraining);
  training_data = (double **)calloc(ntraining, sizeof(training_data[0])) ;
  if (training_classes == NULL)
    ErrorExit(ERROR_NOFILE, "train_rforest: could not allocate %d-length training buffers",ntraining);
  for (i = n = 0 ; n < nsubjects ; n++)
  {
    mri_in = mri_inputs[n][0] ; mri_seg = mri_segs[n][0] ; transform = transforms[n][0] ;
    for (x = 0 ; x < mri_in->width ; x++)
      for (y = 0 ; y <  mri_in->height; y++)
	for (z = 0 ; z <  mri_in->depth; z++)
	{
	  if ((int)MRIgetVoxVal(mri_training_voxels, x, y, z, n) == 0)
	    continue ;
	  label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
	  TransformSourceVoxelToAtlas(transform, mri_in, x, y, z, &xatlas, &yatlas, &zatlas) ;
	  xt = nint(xatlas/tvoxel_size) ; yt = nint(yatlas/tvoxel_size) ; zt = nint(zatlas/tvoxel_size) ;
	  if (IS_FUTURE_WMSA(label))
	    tlabel = FUTURE_WMSA ;
	  else if (IS_WMSA(label))
	    tlabel = WMSA ;
	  else
	    tlabel= NOT_WMSA ;
	    
	  training_classes[i] =  tlabel ;
	  training_data[i] = (double *)calloc(nfeatures, sizeof(double)) ;
	  if (training_data[i] == NULL)
	    ErrorExit(ERROR_NOMEMORY, "train_rforest: could not allocate %d-len feature vector #%d",
		      nfeatures, i) ;
	  training_classes[i] = training_classes[i] ;
//	  extract_feature(mri_in, parms->wsize, x, y, z, training_data[i], xatlas, yatlas, zatlas) ;
	  extract_long_features(mri_in, mri_seg, transform, gca, parms->wsize, x, y, z, training_data[i]) ;
	  if (training_data[i][Gdiag_no] < 80 && training_classes[i] == 1)
	    DiagBreak() ;
	  i++ ;
	}
    MRIfree(&mri_in) ; MRIfree(&mri_seg) ; TransformFree(&transform) ;
  }

  if (i < ntraining)
  {
    printf("warning!!!! i (%d) < ntraining (%d)! Setting ntraining=i\n", i, ntraining) ;
    ntraining = i ;
  }
  printf("training random forest with %dK FUTURE WMSA labels, %dK WMSA labels, and %dK non (ratio=%2.1f:%2.1f)\n",
	 nfuture/1000, nwmsa/1000, nnot/1000, (float)nnot/(float)nfuture, (float)nnot/(float)nwmsa) ;
  RFtrain(rf, parms->feature_fraction, parms->training_fraction, training_classes, training_data, ntraining);
  correct = RFcomputeOutOfBagCorrect(rf, training_classes, training_data,ntraining);
  printf("out of bag accuracy: %d of %d = %2.2f%%\n", correct, ntraining,
	 100.0*correct/ntraining) ;
  if (log_file_name)
  {
      struct flock fl;
      int    fd;
      char   line[MAX_LINE_LEN] ;

      printf("writing results to train.log file %s\n", log_file_name) ;
      fd = open(log_file_name, O_WRONLY|O_APPEND|O_CREAT, S_IRWXU|S_IRWXG);
      if (fd < 0)
	ErrorExit(ERROR_NOFILE, "%s: could not open test log file %s", 
		  Progname, log_file_name);
      
      fcntl(fd, F_SETLKW, &fl);  /* F_GETLK, F_SETLK, F_SETLKW */
      sprintf(line, "%f %d %d %f\n", 
	      rf->training_fraction,
	      rf->max_depth,
	      rf->ntrees,
	       100.0*correct/ntraining) ;
      write(fd, line, (strlen(line))*sizeof(char)) ;
      fl.l_type   = F_UNLCK;  /* tell it to unlock the region */
      fcntl(fd, F_SETLK, &fl); /* set the region to unlocked */
      close(fd) ;
  }

  for (i = 0 ; i < ntraining ; i++)  // allow for augmenting with other wmsa examples
    free(training_data[i]) ;
  free(training_data) ;
  free(training_classes) ;
  MRIfree(&mri_training_voxels) ;
  return(rf) ;
}


#if 0
static int
is_wmsa_border(MRI *mri_seg, int x, int y, int z)
{
  int found_wmsa, found_non_wmsa, label, xi, yi, zi, xk, yk, zk ;

  for (found_wmsa = found_non_wmsa = 0, xk = -1 ; xk <= 1 ; xk++)
  {
    xi = mri_seg->xi[x+xk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri_seg->yi[y+yk] ;
      for (zk = -1 ; zk <= 1 ; zk++)
      {
	if (abs(xk)+abs(yk)+abs(zk) > 1)  // only 6-connected
	  continue ;
	zi = mri_seg->zi[z+zk] ;
	label = (int)MRIgetVoxVal(mri_seg, xi, yi, zi, 0) ;
	if (IS_WMSA(label))
	  found_wmsa++ ;
	else
	  found_non_wmsa++ ;
	if (found_wmsa && found_non_wmsa)
	  return(1) ;
      }
    }
  }
  return(0) ;
}

#endif
