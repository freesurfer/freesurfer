/**
 * @brief Creates the Gaussian Classifier Array (GCA) atlas from training set
 *
 * See:
 * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
 * Structures in the Human Brain", Fischl et al.
 * (2002) Neuron, 33:341-355.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2011-2017 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/stat.h>
#include <sys/types.h>

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
#include "flash.h"
#include "version.h"
#ifdef _OPENMP
#include "romp_support.h"
#endif


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int replaceLabels(MRI *mri_seg) ;
static void configure_transform(TRANSFORM *transform, MRI *mri, GCA *gca) ;
static int lateralize_hypointensities(MRI *mri_seg) ;
static int check(MRI *mri_seg, char *subjects_dir, char *subject_name) ;

static int conform = 1 ;
static int flash = 0 ;
static int binarize = 0 ;
static int binarize_in = 0 ;
static int binarize_out = 0 ;
static char *wmsa_fname = NULL ;

static int gca_flags = GCA_NO_FLAGS ;

static COLOR_TABLE *ctab = NULL ;
const char *Progname ;
static void usage_exit(int code) ;
static char *mask_fname = NULL ;
static char *insert_fname = NULL ;
static int  insert_label = 0 ;
static char *histo_fname = NULL ;

static float scale = 0 ;
static int force_inputs = 1 ;

static GCA_PARMS parms ;
static const char *seg_dir = "seg_edited.mgz" ; // default name of manual edit file
static char T1_name[STRLEN] = "orig" ;
static char *xform_name = NULL;
static int prune = 0 ;
static float smooth = -1 ;
static int gca_inputs = 0 ;
static double TRs[MAX_GCA_INPUTS] ;
static double TEs[MAX_GCA_INPUTS] ;
static double FAs[MAX_GCA_INPUTS] ;
static int map_to_flash = 0 ;

static int ninputs = 1 ;  /* T1 intensity */
static int navgs = 0 ;

static char subjects_dir[STRLEN] ;
static char *heq_fname = NULL ;

static char *input_names[MAX_GCA_INPUTS] =
  {
    T1_name
  } ;

static int do_sanity_check = 0;
static int do_fix_badsubjs = 0;
static int sanity_check_badsubj_count = 0;
static int AllowMisMatch=0;
int DoSym = 0;
#ifdef _OPENMP
  int n_omp_threads;
#endif
char *DoneFile = NULL;
char *cmdline, cwd[2000];


static void writeDoneFile(char *DoneFile, int errorcode)
{
  if (DoneFile == NULL) return;
  FILE *fp = fopen(DoneFile, "w");
  if (fp == NULL) return;
  fprintf(fp, "%d\n", errorcode);
  fclose(fp);
}


int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp ;
  int          ac, nargs, i, n, noint = 0, options ;
  int          msec, minutes, seconds, nsubjects, input,
  ordering[MAX_GCA_INPUTS], o ;
  Timer start ;
  GCA          *gca, *gca_prune = NULL ;
  MRI          *mri_seg, *mri_tmp, *mri_eq = NULL, *mri_inputs ;
  TRANSFORM    *transform ;
  LTA          *lta;
  int          used[MAX_GCA_INPUTS];
  int          counts;

  Progname = argv[0] ;

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  cmdline = argv2cmdline(argc,argv);
  getcwd(cwd,2000);

  start.reset() ;

  parms.use_gradient = 0 ;
  parms.node_spacing = 4.0f ;
  parms.prior_spacing = 2.0f ;

  nargs = handleVersionOption(argc, argv, "mri_ca_train");

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

  // catch exceptions in order to write a done file
  throwExceptions(true);
  try {

  printf("\n");
  printf("setenv SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  printf("cd %s\n",cwd);
  printf("%s\n\n",cmdline);

#ifdef HAVE_OPENMP
  n_omp_threads = omp_get_max_threads(); 
  printf("\n== Number of threads available to %s for OpenMP = %d == \n",
         Progname, n_omp_threads);
#endif

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
    if (argc < 3)
      usage_exit(1) ;
  }

  if (heq_fname)
  {
    mri_eq = MRIread(heq_fname) ;
    if(!mri_eq)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read histogram equalization volume %s",
                Progname, heq_fname) ;
  }
  // options parsed.   subjects and gca name remaining
  out_fname = argv[argc-1] ;
  nsubjects = argc-2 ;
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

  // gca_inputs can be T1, PD, ...per subject
  if (gca_inputs == 0)
    gca_inputs = ninputs ;
  /* gca reads same # of inputs as we read
     from command line - not the case if we are mapping to flash */
  n = 0 ;
  if (gca_flags & GCA_GRAD)
  {
    int extra = 0 ;
    if (gca_flags & GCA_XGRAD)
      extra += gca_inputs ;
    if (gca_flags & GCA_YGRAD)
      extra += gca_inputs ;
    if (gca_flags & GCA_ZGRAD)
      extra += gca_inputs ;
    gca_inputs += extra ;
  }

  //////////////////////////////////////////////////////////////////
  do
  {
    // set up gca direction cosines, width, height, depth defaults
    gca = GCAalloc(gca_inputs, parms.prior_spacing,
                   parms.node_spacing, DEFAULT_VOLUME_SIZE,
                   DEFAULT_VOLUME_SIZE,DEFAULT_VOLUME_SIZE, gca_flags);

    /////////////////////////////////////////////////////////////////////////
    // weird way options and subject name are mixed here

    /////////////////////////////////////////////////////////
    // first calculate mean
    ////////////////////////////////////////////////////////
    // going through the subject one at a time
    for (nargs = i = 0 ; i < nsubjects+options ; i++)
    {
      subject_name = argv[i+1] ;
      //////////////////////////////////////////////////////////////
      printf("***************************************"
             "************************************\n");
      printf("processing subject %s, %d of %d...\n", subject_name,i+1-nargs,
             nsubjects);

      if (stricmp(subject_name, "-NOINT") == 0)
      {
        printf("not using intensity information for subsequent subjects...\n");
        noint = 1 ;
        nargs++ ;
        continue ;
      }
      else if (stricmp(subject_name, "-INT") == 0)
      {
        printf("using intensity information for subsequent subjects...\n");
        noint = 0 ;
        nargs++ ;
        continue ;
      }
      // reading this subject segmentation
      int req = snprintf(fname, STRLEN, "%s/%s/mri/%s", subjects_dir, subject_name, seg_dir) ; 
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stderr, "Reading segmentation from %s...\n", fname) ;
      mri_seg = MRIread(fname) ;
      if (!mri_seg)
        ErrorExit(ERROR_NOFILE, "%s: could not read segmentation file %s",
                  Progname, fname) ;
      if ((mri_seg->type != MRI_UCHAR) && (mri_seg->type != MRI_FLOAT))
      {
        //ErrorExit(ERROR_NOFILE,"%s: segmentation file %s is not type UCHAR or FLOAT",Progname, fname) ;
	printf("Info: changing type of seg to float\n");
	MRI *mritmp;
	mritmp = MRISeqchangeType(mri_seg, MRI_FLOAT, 0,0,0);
	MRIfree(&mri_seg);
	mri_seg = mritmp;
      }

      if (wmsa_fname)
      {
        MRI *mri_wmsa ;
        int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			   subjects_dir, subject_name, wmsa_fname) ;  
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
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
          int req = snprintf(s, STRLEN, "%s/%s/mri/seg_%s",
			     subjects_dir, subject_name, wmsa_fname) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
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

        int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			   subjects_dir, subject_name, insert_fname) ;  
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
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

      // subjects loop index i
      if (i != 0)  /* not the first image read -
                      reorder it to be in the same order as 1st */
      {
        // initialize the flag array
        for (input =0; input < ninputs; input++) {
          used[input] = 0;
	}

        for (input = 0 ; input < ninputs ; input++)
        {
          int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			     subjects_dir, subject_name, input_names[input]);    
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          mri_tmp = MRIreadInfo(fname) ;
          if (!mri_tmp)
            ErrorExit(ERROR_NOFILE,
             "%s: could not read image from file %s", Progname, fname) ;
	  if (force_inputs)
	  {
	    ordering[input] = input ;
	    used[input] = 1 ;
	  }
	  else for (o = 0 ; o < ninputs ; o++)
            if (FEQUAL(TRs[o],mri_tmp->tr) &&
                FEQUAL(FAs[o],mri_tmp->flip_angle) &&
                FEQUAL(TEs[o],mri_tmp->te))
            {
              // if this o is not used, then use it
              if (used[o] == 0)
              {
                ordering[input] = o ;
                used[o] = 1;
                break;
              }
            }
          MRIfree(&mri_tmp) ;
        }
        // verify whether it has input values are used
        counts = 0;
        for (input = 0; input < ninputs; input++)
          if (used[input] == 1) counts++;
        if (counts != ninputs){
	  if(!AllowMisMatch)
	    ErrorExit(ERROR_BADPARM,
		      "Input TR, TE, FlipAngle for each subjects must match.\n");
	  printf("Input TR, TE, FlipAngle for each subjects do not match, but mismatch allowed.\n");	  
	}
      } else
        for (o = 0 ; o < ninputs ; o++)
          ordering[o] = o ;


      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      {
        printf("ordering images: ") ;
        for (o = 0 ; o < ninputs ; o++)
          printf("%d ", ordering[o]) ;
        printf("\n") ;
      }

      if (flash)
        gca->type = GCA_FLASH ;
      else
      {
        if (ninputs==2) // T1 and PD?
        {
          if ((strstr(input_names[ordering[0]], "T1") &&
               strstr(input_names[ordering[1]], "PD")) ||
              (strstr(input_names[ordering[1]], "T1") &&
               strstr(input_names[ordering[0]], "PD")))
          {
            gca->type=GCA_PARAM;
          }
        }
        else if (ninputs==1) // single input
          gca->type=GCA_NORMAL;
        else
          gca->type=GCA_UNKNOWN;
      }
      if (DIAG_VERBOSE_ON)
        printf("gca->type = %d\n", gca->type) ;

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

        int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			  subjects_dir, subject_name,input_names[ordering[input]]);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
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

          if (gca_flags & GCA_XGRAD)
            nframes += ninputs ;
          if (gca_flags & GCA_YGRAD)
            nframes += ninputs ;
          if (gca_flags & GCA_ZGRAD)
            nframes += ninputs ;
          ///////////////////////////////////////////////////////////
          mri_inputs =
            MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth,
                             mri_tmp->type, nframes) ;
          if (!mri_inputs)
            ErrorExit
            (ERROR_NOMEMORY,
             "%s: could not allocate input volume %dx%dx%dx%d",
             mri_tmp->width, mri_tmp->height, mri_tmp->depth,nframes) ;
          MRIcopyHeader(mri_tmp, mri_inputs) ;
        }
        // -mask option ////////////////////////////////////////////
        if (mask_fname)
        {
          MRI *mri_mask ;

          int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			     subjects_dir, subject_name, mask_fname); 
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }

          printf("reading volume %s for masking...\n", fname) ;
          mri_mask = MRIread(fname) ;
          if (!mri_mask)
            ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                      Progname, fname) ;

          MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
          MRIfree(&mri_mask) ;
        }
        if (mri_eq && !noint)
        {
          printf("histogram equalizing input image...\n") ;
          MRIhistoEqualize(mri_tmp, mri_eq, mri_tmp, 30, 170) ;
        }
	fflush(stdout);fflush(stderr);
        mri_inputs=MRIcopyFrame(mri_tmp, mri_inputs, 0, input) ;
	fflush(stdout);fflush(stderr);
	if(mri_inputs == NULL) exit(1);
        MRIfree(&mri_tmp) ;
      }// end of inputs per subject

      MRIeraseBorderPlanes(mri_inputs, 1) ;
      // when loaded gca->type = GCA_FLASH, these should have been set?????
      if (i == 0 && flash)   /* first subject */
        GCAsetFlashParameters(gca, TRs, FAs, TEs) ;

      /////////////////////////////////////////////////////////
      // xform_name is given, then we can use the consistent c_(r,a,s) for gca
      /////////////////////////////////////////////////////////
      if (xform_name)
      {
        // we read talairach.xfm which is a RAS-to-RAS
        int req = snprintf(fname, STRLEN, "%s/%s/mri/transforms/%s",
			   subjects_dir, subject_name, xform_name) ; 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          printf("INFO: reading transform file %s...\n", fname);
        if (!FileExists(fname))
        {
          fprintf(stderr,"ERROR: cannot find transform file %s\n",fname);
          exit(1);
        }
        transform = TransformRead(fname);
        if (!transform)
          ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s", Progname, fname);

        configure_transform(transform, mri_inputs, gca);
        // Here we do 2 things
        // 1. modify gca direction cosines to
        // that of the transform destination (both linear and non-linear)
        // 2. if ras-to-ras transform,
        // then change it to vox-to-vox transform (linear case)

        // modify transform to store inverse also
	printf("Inverting transform\n");fflush(stdout);
        TransformInvert(transform, mri_inputs) ;
        // verify inverse
        lta = (LTA *) transform->xform;
      }
      else
      {
        GCAreinit(mri_inputs, gca);
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

      /////////////////////////////////////////////////////////
      if (map_to_flash)
      {
        MRI *mri_tmp ;

        mri_tmp = MRIparameterMapsToFlash
                  (mri_inputs, NULL, TRs, TEs, FAs, gca_inputs) ;
        MRIfree(&mri_inputs) ;
        mri_inputs = mri_tmp ;
      }

      /////////////////////////////////////////////////////////
      if (gca_flags & GCA_GRAD)
      {
        MRI *mri_kernel, *mri_smooth, *mri_grad, *mri_tmp ;
        int i, start ;

        mri_kernel = MRIgaussian1d(1.0, 30) ;
        mri_smooth = MRIconvolveGaussian(mri_inputs, NULL, mri_kernel) ;

        if (mri_inputs->type != MRI_FLOAT)
        {
          mri_tmp = MRISeqchangeType(mri_inputs, MRI_FLOAT, 0, 0, 1) ;
          MRIfree(&mri_inputs) ;
          mri_inputs = mri_tmp ;
        }

        start = ninputs ;
        if (gca_flags & GCA_XGRAD)
        {
          for (i = 0 ; i < ninputs ; i++)
          {
            mri_grad = MRIxSobel(mri_smooth, NULL, i) ;
            mri_inputs=MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
	    if(mri_inputs == NULL) exit(1);
            MRIfree(&mri_grad) ;
          }
          start += ninputs ;
        }
        if (gca_flags & GCA_YGRAD)
        {
          for (i = 0 ; i < ninputs ; i++)
          {
            mri_grad = MRIySobel(mri_smooth, NULL, i) ;
            mri_inputs=MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
	    if(mri_inputs == NULL) exit(1);
            MRIfree(&mri_grad) ;
          }
          start += ninputs ;
        }
        if (gca_flags & GCA_ZGRAD)
        {
          for (i = 0 ; i < ninputs ; i++)
          {
            mri_grad = MRIzSobel(mri_smooth, NULL, i) ;
            mri_inputs=MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
	    if(mri_inputs == NULL) exit(1);
            MRIfree(&mri_grad) ;
          }
          start += ninputs ;
        }

        MRIfree(&mri_kernel) ;
        MRIfree(&mri_smooth) ;
      }

      ////////////////////////////////////////////////////////////////////
      // train gca
      ////////////////////////////////////////////////////////////////////
      // segmentation is seg volume
      // inputs       is the volumes of all inputs
      // transform    is for this subject
      // gca_prune    is so far null
      // noint        is whether to use intensity information or not
      GCAtrain(gca, mri_inputs, mri_seg, transform, gca_prune, noint) ;
      GCAcheck(gca) ;
      MRIfree(&mri_seg) ;
      MRIfree(&mri_inputs) ;
      TransformFree(&transform) ;

      // compactify gca
      gca = GCAcompactify(gca);
      if (gca_prune)
        gca_prune = GCAcompactify(gca_prune);
    }
    GCAcompleteMeanTraining(gca) ;

    ///////////////////////////////////////////////////////////////
    if (do_sanity_check && sanity_check_badsubj_count)
    {
      ErrorExit(-9,
                "\nERROR: mri_ca_train check found %d subjects "
                "with bad labels!\n",
                sanity_check_badsubj_count);     
    }

    ///////////////////////////////////////////////////////////////
    /* now compute covariances */
    ///////////////////////////////////////////////////////////////
    for (nargs = i = 0 ; i < nsubjects+options ; i++)
    {
      subject_name = argv[i+1] ;
      if (stricmp(subject_name, "-NOINT") == 0)
      {
        printf("not using intensity information for subsequent subjects...\n");
        noint = 1 ;
        nargs++ ;
        continue ;
      }
      else if (stricmp(subject_name, "-INT") == 0)
      {
        printf("using intensity information for subsequent subjects...\n");
        noint = 0 ;
        nargs++ ;
        continue ;
      }
      if (noint)
      {
        printf("skipping covariance calculation for subject %s...\n",
               subject_name) ;
        continue ;
      }
      printf("computing covariances for subject %s, %d of %d...\n",
             subject_name,i+1-nargs,
             nsubjects);
      int req = snprintf(fname, STRLEN, "%s/%s/mri/%s", subjects_dir, subject_name, seg_dir) ;  
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (DIAG_VERBOSE_ON)
        printf("reading segmentation from %s...\n", fname) ;
      // seg volume
      ///////////////////////////////
      mri_seg = MRIread(fname) ;
      if (!mri_seg)
        ErrorExit(ERROR_NOFILE, "%s: could not read segmentation file %s",
                  Progname, fname) ;
      if (wmsa_fname)
      {
        MRI *mri_wmsa ;
        int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			   subjects_dir, subject_name, wmsa_fname) ; 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
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
          int req = snprintf(s, STRLEN, "%s/%s/mri/seg_%s",
			     subjects_dir, subject_name, wmsa_fname) ; 
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          MRIwrite(mri_seg, s) ;
        }
      }
      if (insert_fname)
      {
        MRI *mri_insert ;

        int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			   subjects_dir, subject_name, insert_fname) ;
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        mri_insert = MRIread(fname) ;
        if (mri_insert == NULL)
          ErrorExit
          (ERROR_NOFILE,
           "%s: could not read volume from %s for insertion",
           Progname, insert_fname) ;

        MRIbinarize(mri_insert, mri_insert, 1, 0, insert_label) ;
        MRIcopyLabel(mri_insert, mri_seg, insert_label) ;
        MRIfree(&mri_insert) ;
      }

      replaceLabels(mri_seg) ;
      MRIeraseBorderPlanes(mri_seg, 1) ;

      // input volume
      /////////////////////////////////
      // inputs are T1, PD, .... per subject
      for (input = 0 ; input < ninputs ; input++)
      {
        int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			   subjects_dir, subject_name,input_names[ordering[input]]);  
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        if (DIAG_VERBOSE_ON)
          printf("reading co-registered input from %s...\n", fname) ;
        mri_tmp = MRIread(fname) ;
        if (!mri_tmp)
          ErrorExit(ERROR_NOFILE, "%s: could not read T1 data from file %s",
                    Progname, fname) ;

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
            (mri_tmp->xsize != 1 || mri_tmp->ysize != 1 ||mri_tmp->zsize != 1))
        {
          ErrorExit
          (ERROR_BADPARM,
           "%s: must have 1mm voxel size, but have (%f, %f, %f)\n",
           fname, mri_tmp->xsize, mri_tmp->ysize, mri_tmp->ysize);
        }
        // input check 3 is removed c_(ras) != 0 can be handled
        if (input == 0)
        {
          int nframes = ninputs ;

          if (gca_flags & GCA_XGRAD)
            nframes += ninputs ;
          if (gca_flags & GCA_YGRAD)
            nframes += ninputs ;
          if (gca_flags & GCA_ZGRAD)
            nframes += ninputs ;
          mri_inputs =
            MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth,
                             mri_tmp->type, nframes) ;
          if (!mri_inputs)
            ErrorExit
            (ERROR_NOMEMORY,
             "%s: could not allocate input volume %dx%dx%dx%d",
             mri_tmp->width, mri_tmp->height, mri_tmp->depth,ninputs) ;
          MRIcopyHeader(mri_tmp, mri_inputs) ;
        }

        if (mask_fname)
        {
          MRI *mri_mask ;

          int req = snprintf(fname, STRLEN, "%s/%s/mri/%s",
			     subjects_dir, subject_name, mask_fname);  
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          printf("reading volume %s for masking...\n", fname) ;
          mri_mask = MRIread(fname) ;
          if (!mri_mask)
            ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                      Progname, fname) ;

          MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
          MRIfree(&mri_mask) ;
        }
        if (mri_eq && !noint)
        {
          printf("histogram equalizing input image...\n") ;
          MRIhistoEqualize(mri_tmp, mri_eq, mri_tmp, 30, 170) ;
        }
        mri_inputs=MRIcopyFrame(mri_tmp, mri_inputs, 0, input) ;
	if(mri_inputs == NULL) exit(1);
        MRIfree(&mri_tmp) ;
      } // end of building inputs
      ///////////////////////////////////////////////////////////
      if (xform_name)
      {
        int req = snprintf(fname, STRLEN, "%s/%s/mri/transforms/%s",
			   subjects_dir, subject_name, xform_name) ; 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          printf("reading transform from %s...\n", fname) ;
        transform = TransformRead(fname) ;
        if (!transform)
          ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
                    Progname, fname) ;
        // change the transform to vox-to-vox
        configure_transform(transform, mri_inputs, gca);
	printf("Inverting transform\n");fflush(stdout);
        TransformInvert(transform, mri_inputs) ;
        if ((transform->type != MORPH_3D_TYPE) &&
            ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON))
        {
          // verify inverse
          lta = (LTA *) transform->xform;
          {
            MATRIX *go = lta->xforms[0].m_L;
            MATRIX *back = lta->inv_xforms[0].m_L;
            MATRIX *seki = MatrixMultiply(back, go, NULL);
            fprintf(stderr,
                    "You should see the unit matrix to verify the inverse.\n");
            MatrixPrint(stderr, seki);
            MatrixFree(&seki);
          }
        }
      }
      else
        transform = TransformAlloc(LINEAR_VOXEL_TO_VOXEL, NULL) ;

      ////////////////////////////////////////////////////////
      if (map_to_flash)
      {
        MRI *mri_tmp ;

        mri_tmp = MRIparameterMapsToFlash
                  (mri_inputs, NULL, TRs, TEs, FAs, gca_inputs) ;
        MRIfree(&mri_inputs) ;
        mri_inputs = mri_tmp ;
      }
      /////////////////////////////////////////////////////////
      if (gca_flags & GCA_GRAD)
      {
        MRI *mri_kernel, *mri_smooth, *mri_grad, *mri_tmp ;
        int i, start ;

        mri_kernel = MRIgaussian1d(1.0, 30) ;
        mri_smooth = MRIconvolveGaussian(mri_inputs, NULL, mri_kernel) ;

        if (mri_inputs->type != MRI_FLOAT)
        {
          mri_tmp = MRISeqchangeType(mri_inputs, MRI_FLOAT, 0, 0, 1) ;
          MRIfree(&mri_inputs) ;
          mri_inputs = mri_tmp ;
        }

        start = ninputs ;
        if (gca_flags & GCA_XGRAD)
        {
          for (i = 0 ; i < ninputs ; i++)
          {
            mri_grad = MRIxSobel(mri_smooth, NULL, i) ;
            mri_inputs=MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
	    if(mri_inputs == NULL) exit(1);
            MRIfree(&mri_grad) ;
          }
          start += ninputs ;
        }
        if (gca_flags & GCA_YGRAD)
        {
          for (i = 0 ; i < ninputs ; i++)
          {
            mri_grad = MRIySobel(mri_smooth, NULL, i) ;
            mri_inputs=MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
	    if(mri_inputs == NULL) exit(1);
            MRIfree(&mri_grad) ;
          }
          start += ninputs ;
        }
        if (gca_flags & GCA_ZGRAD)
        {
          for (i = 0 ; i < ninputs ; i++)
          {
            mri_grad = MRIzSobel(mri_smooth, NULL, i) ;
            mri_inputs=MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
	    if(mri_inputs == NULL) exit(1);
            MRIfree(&mri_grad) ;
          }
          start += ninputs ;
        }

        MRIfree(&mri_kernel) ;
        MRIfree(&mri_smooth) ;
      }
      //////////////////////////////////////////////////////////////////
      GCAtrainCovariances(gca, mri_inputs, mri_seg, transform) ;
      /////////////////////////////////////////////////////////////////
      MRIfree(&mri_seg) ;
      MRIfree(&mri_inputs) ;
      TransformFree(&transform) ;

      // compactify gca
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        printf("compactify gca\n");
      gca = GCAcompactify(gca);
      if (gca_prune)
        gca_prune = GCAcompactify(gca_prune);
    }
    GCAcompleteCovarianceTraining(gca) ;
    if (gca_prune)
      GCAfree(&gca_prune) ;
    gca_prune = gca ;  // gca_prune is non-zero now

  }
  while (n++ < prune) ;
  ////////////////  end of do ////////////////////////////////////////////

  if (smooth > 0)
  {
    printf("regularizing conditional densities with smooth=%2.2f\n", smooth) ;
    GCAregularizeConditionalDensities(gca, smooth) ;
  }
  if (navgs)
  {
    printf("applying mean filter %d times to conditional densities\n", navgs) ;
    GCAmeanFilterConditionalDensities(gca, navgs) ;
  }

  if(DoSym){
    printf("Symmetrizing GCA\n");
    GCA *gcatmp = GCAsymmetrize(gca);
    GCAfree(&gca);
    gca = gcatmp;
    int err = GCAisNotSymmetric(gca);
    if(err){
      printf("Symmetrization failed %d\n",err);
      exit(1);
    }
  }

  printf("writing trained GCA to %s...\n", out_fname) ;
  GCAcheck(gca) ;
  gca->ct = ctab ;  // read in with -ctab option
  if (GCAwrite(gca, out_fname) != NO_ERROR)
    ErrorExit
    (ERROR_BADFILE, "%s: could not write gca to %s", Progname, out_fname) ;

  {
    // Why is this here?
    MRI *mri ;
    mri = GCAbuildMostLikelyVolume(gca, NULL) ;
    MRIfree(&mri) ;
  }

  if (histo_fname)
  {
    FILE *fp ;
    int   histo_counts[10000], xn, yn, zn, max_count ;
    GCA_NODE  *gcan ;

    memset(histo_counts, 0, sizeof(histo_counts)) ;
    fp = fopen(histo_fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open histo file %s",
                Progname, histo_fname) ;

    max_count = 0 ;
    for (xn = 0 ; xn < gca->node_width;  xn++)
    {
      for (yn = 0 ; yn < gca->node_height ; yn++)
      {
        for (zn = 0 ; zn < gca->node_depth ; zn++)
        {
          gcan = &gca->nodes[xn][yn][zn] ;
          if (gcan->nlabels < 1)
            continue ;
          if (gcan->nlabels == 1 && IS_UNKNOWN(gcan->labels[0]))
            continue ;
          histo_counts[gcan->nlabels]++ ;
          if (gcan->nlabels > max_count)
            max_count = gcan->nlabels ;
        }
      }
    }
    max_count = 20 ;
    for (xn = 1 ; xn < max_count ;  xn++)
      fprintf(fp, "%d %d\n", xn, histo_counts[xn]) ;
    fclose(fp) ;
  }

  GCAfree(&gca) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec / 1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("classifier array training took %d minutes and %d seconds.\n", minutes, seconds) ;
  
  } catch(...) {
    writeDoneFile(DoneFile, 1);
    printf("mri_ca_train failed\n");
    exit(1);
  }

  writeDoneFile(DoneFile, 0);
  printf("mri_ca_train done\n");
  return 0;
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
  if (!stricmp(option, "threads") || !stricmp(option, "nthreads")){
#ifdef _OPENMP
    sscanf(argv[2],"%d",&n_omp_threads);
    omp_set_num_threads(n_omp_threads);
#endif
    nargs = 1;
  }
  else if(!stricmp(option, "done")) {
    // This file gets created when process is finished.
    // Text content is either 0 (no error) or 1 (error)
    // Calling process must make sure it does not exist
    DoneFile = argv[2];
    nargs = 1;
  }
  else if (!stricmp(option, "GRADIENT"))
  {
    parms.use_gradient = 1 ;
    ninputs += 3 ;  /* components of the gradient */
  }
  else if (!stricmp(option, "PRIOR_SPACING"))
  {
    parms.prior_spacing = atof(argv[2]) ;
    nargs = 1 ;
    printf("spacing priors every %2.1f mm\n", parms.prior_spacing) ;
  }
  else if (!stricmp(option, "CONFORM"))
  {
    conform = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%sassuming input volumes are conformed\n", 
           conform ? "" : "NOT ") ;
  }
  else if (!stricmp(option, "XGRAD"))
  {
    gca_flags |= GCA_XGRAD ;
    printf("using x gradient information in training...\n") ;
  }
  else if (!stricmp(option, "WMSA"))
  {
    wmsa_fname = argv[2] ;
    nargs = 1 ;
    printf("reading white matter signal abnormalities from %s\n", wmsa_fname) ;
  }
  else if (!stricmp(option, "SYM"))
  {
    // Make atlas symmetric prior to writing out
    DoSym = 1;
    nargs = 0 ;
    printf("Creating symmetric atlas\n");
  }
  else if (!stricmp(option, "MAKESYM"))
  {
    // Read in a GCA and make it symmetric
    GCA *gca;
    gca = GCAread(argv[2]);
    GCA *gcasym = GCAsymmetrize(gca);
    GCAwrite(gcasym,argv[3]);
    int err = GCAisNotSymmetric(gcasym);
    exit(err);
  }
  else if (!stricmp(option, "CHECKSYM"))
  {
    // Read in a GCA and check whether it is symmetric
    GCA *gca = GCAread(argv[2]);
    int err = GCAisNotSymmetric(gca);
    exit(err);
  }
  else if (!stricmp(option, "YGRAD"))
  {
    gca_flags |= GCA_YGRAD ;
    printf("using y gradient information in training...\n") ;
  }
  else if (!stricmp(option, "ZGRAD"))
  {
    gca_flags |= GCA_ZGRAD ;
    printf("using z gradient information in training...\n") ;
  }
  else if (!stricmp(option, "FLASH"))
  {
#if 1
    flash = 1 ;
    printf("setting gca->type to FLASH\n") ;
#else
    int i ;

    map_to_flash = 1 ;
    gca_inputs = atoi(argv[2]) ;
    nargs = 1+3*gca_inputs ;
    printf("mapping T1/PD inputs to flash volumes:\n") ;
    for (i = 0 ; i < gca_inputs ; i++)
    {
      TRs[i] = atof(argv[3+3*i]) ;
      FAs[i] = RADIANS(atof(argv[4+3*i])) ;
      TEs[i] = atof(argv[5+3*i]) ;
      printf("\tvolume %d: TR=%2.1f msec, flip angle %2.1f, TE=%2.1f msec\n",
             i, TRs[i], DEGREES(FAs[i]), TEs[i]) ;
    }
#endif
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
  else if (!stricmp(option, "NODE_SPACING"))
  {
    parms.node_spacing = atof(argv[2]) ;
    nargs = 1 ;
    printf("spacing nodes every %2.1f mm\n", parms.node_spacing) ;
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
  else if (!stricmp(option, "NOMRF"))
  {
    gca_flags |= GCA_NO_MRF ;
    printf("not computing MRF statistics...\n") ;
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
  else if (!stricmp(option, "DEBUG_PRIOR"))
  {
    Gxp = atoi(argv[2]) ;
    Gyp = atoi(argv[3]) ;
    Gzp = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging prior (%d, %d, %d)\n", Gxp,Gyp,Gzp) ;
  }
  else if (!stricmp(option, "DEBUG_LABEL"))
  {
    Ggca_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging label %s (%d)\n", cma_label_to_name(Ggca_label),
           Ggca_label) ;
  }
  else if (!stricmp(option, "DEBUG_NBR"))
  {
    Ggca_nbr_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging nbr label %s (%d)\n",
           cma_label_to_name(Ggca_nbr_label), Ggca_nbr_label) ;
  }
  else if (!stricmp(option, "INSERT"))
  {
    insert_fname = argv[2] ;
    insert_label = atoi(argv[3]) ;
    nargs = 2 ;
    printf("inserting non-zero vals from %s as label %d...\n",
           insert_fname,insert_label);
  }
  else if (!stricmp(option, "PRUNE"))
  {
    prune = atoi(argv[2]) ;
    nargs = 1 ;
    printf("pruning classifier %d times after initial training\n", prune) ;
  }
  else if (!stricmp(option, "HEQ"))
  {
    heq_fname = argv[2] ;
    nargs = 1 ;
    printf("reading template for histogram equalization from %s...\n",
           heq_fname) ;
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
  else if (!stricmp(option, "mismatch"))
  {
    AllowMisMatch=1;
    printf("will allow MR param mismatch\n") ;
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
  else switch (toupper(*option))
    {
    case 'F':
      force_inputs = 1 ;
      printf("forcing use of inputs even if acquisition parameters don't match\n");
      break ;
    case 'S':
      scale = atof(argv[2]) ;
      printf("scaling all volumes by %2.3f after reading...\n", scale) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      printf("applying %d mean filters to classifiers after training\n",navgs);
      nargs = 1 ;
      break ;
    case 'H':
      histo_fname = argv[2] ;
      nargs = 1 ;
      printf("writing histogram of classes/voxel to %s\n", histo_fname) ;
      break;
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
usage_exit(int code)
{
  printf("Purpose: %s trains GCA data with (multiple) subject(s)\n", Progname);
  printf("Usage  : %s [options] <subject1> <subject2> ... "
         "<output gca fname>\n",
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
   "\t -sym - symmetrize the atlas after creation"
   "\t -makesym input.gca symmetrized.gca : symmetrize an already existing atlas "
   "\t -checksym input.gca symmetrized.gca : check the symmetry of an already existing atlas"
   "If not specified, \"orig\" is used\n"
   "\t-check          - conduct sanity-check of labels for obvious edit errors"
   "\t-threads N : specify number of threads to use (also -nthreads)"
   "\t-done DoneFile : create DoneFile when done, (contents: 0=ok, 1=error)"
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


/*
  Reinitializes the GCA geometry with the transform's target geometry. If the transform
  has no destination information, the geometry of the source image is used and CRAS is set to
  the origin. If the transform is a linear RAS->RAS matrix, it is converted to VOX->VOX.
*/
static void configure_transform(TRANSFORM *transform, MRI *src, GCA *gca)
{
  static bool warned = false;

  // get transform target geometry
  VOL_GEOM *geom = nullptr;
  if (transform->type == MORPH_3D_TYPE) {
    geom = &((GCA_MORPH *)transform->xform)->atlas;
  } else {
    geom = &((LTA *)transform->xform)->xforms[0].dst;
  }

  MRI *target = nullptr;
  if (geom->valid) {
    // create target volume with transform target geometry
    target = MRIallocFromVolGeom(geom, src->type, 1, 1);
  } else {
    // if target does not exist, just use source geometry with cras set to origin
    target = MRIallocHeader(src->width, src->height, src->depth, src->type, 1);
    MRIcopyHeader(src, target);
    bool in305space = ((transform->type == LINEAR_RAS_TO_RAS) && (getenv("USE_AVERAGE305")));
    if (in305space && (!warned)) {
      fprintf(stdout, "INFO: Using average_305 CRAS. Disable by unsetting USE_AVERAGE305 env variable.\n");
      warned = true;
    }
    target->c_r = in305space ? -0.095 : 0.0;
    target->c_a = in305space ? -16.51 : 0.0;
    target->c_s = in305space ?   9.75 : 0.0;
  }

  // reinit the GCA with new geometry
  GCAreinit(target, gca);

  // if transform is ras->ras, convert to vox->vox
  if (transform->type == LINEAR_RAS_TO_RAS) {
    LTA *lta = (LTA *)transform->xform;
    // vox -> RAS -> TalRAS
    MATRIX *i_to_r = extract_i_to_r(src);
    MATRIX *tmpmat = MatrixMultiply(lta->xforms[0].m_L, i_to_r, NULL);
    MATRIX *r_to_i = extract_r_to_i(target);
    // TalRAS -> vox
    MATRIX *vox2vox = MatrixMultiply(r_to_i, tmpmat,NULL );
    MatrixCopy(vox2vox, lta->xforms[0].m_L);
    // mark it as vox-to-vox
    transform->type = LINEAR_VOX_TO_VOX;
    // free up memory
    MatrixFree(&r_to_i);
    MatrixFree(&i_to_r);
    MatrixFree(&tmpmat);
    MatrixFree(&vox2vox);
  }

  MRIfree(&target);
}


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
#if __GNUC__  >= 8
	    [[gnu::fallthrough]];
#endif


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
#if __GNUC__  >= 8
	    [[gnu::fallthrough]];
#endif


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
    int req = snprintf(fname, STRLEN, "%s/%s/mri/seg_fixed.mgz", subjects_dir, subject_name);  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("Writing corrected volume to %s\n",fname);
    MRIwrite(mri_fixed,fname);
    MRIfree(&mri_fixed);
  }

  fflush(stdout);

  return(errors) ;
}

