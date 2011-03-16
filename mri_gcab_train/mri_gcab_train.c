/**
 * @file  mri_gcab_train.c
 * @brief train a gca boundary atlas
 *
 * Routines for supporting boundary deformations to refine the
 * exact border of aseg structures.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/03/16 20:23:33 $
 *    $Revision: 1.4 $
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


/***********************************************************************/
/* mri_gcab_train.c                                                      */
/* by Bruce Fischl                                                     */
/*                                                                     */
/* Warning: Do not edit the following four lines.  CVS maintains them. */
/* Revision Author: $Author: fischl $                                  */
/* Revision Date  : $Date: 2011/03/16 20:23:33 $                       */
/* Revision       : $Revision: 1.4 $                                  */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "gcaboundary.h"
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

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int replaceLabels(MRI *mri_seg) ;
static void modify_transform(TRANSFORM *transform, MRI *mri, GCA *gca);

static int binarize = 0 ;
static int binarize_in = 0 ;
static int binarize_out = 0 ;

static int gca_flags = GCA_NO_FLAGS ;

char *Progname ;
static void usage_exit(int code) ;
static char *mask_fname = NULL ;
static char *insert_fname = NULL ;
static int  insert_label = 0 ;
static char *histo_fname = NULL ;

static float scale = 0 ;

static GCA_PARMS parms ;
static char *seg_dir = "seg" ;
static char T1_name[STRLEN] = "orig" ;
static char *xform_name = NULL;
static float smooth = -1 ;
static double TRs[MAX_GCA_INPUTS] ;
static double TEs[MAX_GCA_INPUTS] ;
static double FAs[MAX_GCA_INPUTS] ;

static int navgs = 0 ;

static char subjects_dir[STRLEN] ;

static char *input_names[MAX_GCA_INPUTS] = {
      T1_name
    } ;

static int target_label = Left_Hippocampus ;
static float spacing = 8.0 ;
int
main(int argc, char *argv[]) {
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp ;
  int          ac, nargs, i, n, noint = 0, options ;
  int          msec, minutes, seconds, nsubjects, input ;
  struct timeb start ;
  GCA          *gca ;
  MRI          *mri_seg, *mri_tmp, *mri_inputs ;
  TRANSFORM    *transform ;
  LTA          *lta;
  GCA_BOUNDARY *gcab ;

  Progname = argv[0] ;

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  parms.use_gradient = 0 ;
  spacing = 8 ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mri_gcab_train.c,v 1.4 2011/03/16 20:23:33 fischl Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  // parse command line args
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  printf("reading gca from %s\n", argv[1]) ;
  gca = GCAread(argv[1]) ;
  if (!gca)
    exit(Gerror) ;

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
    if (argc < 4)
      usage_exit(1) ;
  }

  // options parsed.   subjects and gca name remaining
  out_fname = argv[argc-1] ;
  nsubjects = argc-3 ;
  for (options = i = 0 ; i < nsubjects ; i++) {
    if (argv[i+1][0] == '-') {
      nsubjects-- ;
      options++ ;
    }
  }

  printf("training on %d subject and writing results to %s\n",
         nsubjects, out_fname) ;

  n = 0 ;

  gcab = GCABalloc(gca, 8, 0, 30, 10, target_label);
  strcpy(gcab->gca_fname, argv[1]) ;
  // going through the subject one at a time
  for (nargs = i = 0 ; i < nsubjects+options ; i++) {
    subject_name = argv[i+2] ;
    //////////////////////////////////////////////////////////////
    printf("***************************************"
           "************************************\n");
    printf("processing subject %s, %d of %d...\n", subject_name,i+1-nargs,
           nsubjects);

    if (stricmp(subject_name, "-NOINT") == 0) {
      printf("not using intensity information for subsequent subjects...\n");
      noint = 1 ;
      nargs++ ;
      continue ;
    } else if (stricmp(subject_name, "-INT") == 0) {
      printf("using intensity information for subsequent subjects...\n");
      noint = 0 ;
      nargs++ ;
      continue ;
    }
    // reading this subject segmentation
    sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, seg_dir) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "Reading segmentation from %s...\n", fname) ;
    mri_seg = MRIread(fname) ;
    if (!mri_seg)
      ErrorExit(ERROR_NOFILE, "%s: could not read segmentation file %s",
                Progname, fname) ;
    if ((mri_seg->type != MRI_UCHAR) && (mri_seg->type != MRI_FLOAT)) {
      ErrorExit
        (ERROR_NOFILE,
         "%s: segmentation file %s is not type UCHAR or FLOAT",
         Progname, fname) ;
    }

    if (binarize) {
      int j ;
      for (j = 0 ; j < 256 ; j++) {
        if (j == binarize_in)
          MRIreplaceValues(mri_seg, mri_seg, j, binarize_out) ;
        else
          MRIreplaceValues(mri_seg, mri_seg, j, 0) ;
      }
    }
    if (insert_fname) {
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

    for (input = 0 ; input < gca->ninputs ; input++) {
      //////////// set the gca type //////////////////////////////
      // is this T1/PD training?
      // how can we allow flash data training ???????
      // currently checks the TE, TR, FA to be the same for all inputs
      // thus we cannot allow flash data training.
      ////////////////////////////////////////////////////////////

      sprintf(fname, "%s/%s/mri/%s",
              subjects_dir, subject_name,input_names[input]);
      if (DIAG_VERBOSE_ON)
        printf("reading co-registered input from %s...\n", fname) ;
      fprintf(stderr, "   reading input %d: %s\n", input, fname);
      mri_tmp = MRIread(fname) ;
      if (!mri_tmp)
        ErrorExit
          (ERROR_NOFILE,
           "%s: could not read image from file %s", Progname, fname) ;
      // input check 1
      if (getSliceDirection(mri_tmp) != MRI_CORONAL) {
        ErrorExit
          (ERROR_BADPARM,
           "%s: must be in coronal direction, but it is not\n",
           fname);
      }
      // input check 2
      if (mri_tmp->xsize != 1 || mri_tmp->ysize != 1 || mri_tmp->zsize != 1) {
        ErrorExit
          (ERROR_BADPARM,
           "%s: must have 1mm voxel size, but have (%f, %f, %f)\n",
           fname, mri_tmp->xsize, mri_tmp->ysize, mri_tmp->ysize);
      }
      // input check 3 is removed.  now we can handle c_(ras) != 0 case
      // input check 4
      if (i == 0) {
        TRs[input] = mri_tmp->tr ;
        FAs[input] = mri_tmp->flip_angle ;
        TEs[input] = mri_tmp->te ;
      } else if (!FEQUAL(TRs[input],mri_tmp->tr) ||
                 !FEQUAL(FAs[input],mri_tmp->flip_angle) ||
                 !FEQUAL(TEs[input], mri_tmp->te))
        ErrorExit
          (ERROR_BADPARM,
           "%s: subject %s input volume %s: sequence parameters "
           "(%2.1f, %2.1f, %2.1f)"
           "don't match other inputs (%2.1f, %2.1f, %2.1f)",
           Progname, subject_name, fname,
           mri_tmp->tr, DEGREES(mri_tmp->flip_angle), mri_tmp->te,
           TRs[input], DEGREES(FAs[input]), TEs[input]) ;
      // first time do the following
      if (input == 0) {
        int nframes = gca->ninputs ;

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
      MRIcopyFrame(mri_tmp, mri_inputs, 0, input) ;
      MRIfree(&mri_tmp) ;
    }// end of inputs per subject


    /////////////////////////////////////////////////////////
    // xform_name is given, then we can use the consistent c_(r,a,s) for gca
    /////////////////////////////////////////////////////////
    if (xform_name) 
    {
      // we read talairach.xfm which is a RAS-to-RAS
      sprintf(fname, "%s/%s/mri/transforms/%s",
              subjects_dir, subject_name, xform_name) ;
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

      modify_transform(transform, mri_inputs, gca);
      // Here we do 2 things
      // 1. modify gca direction cosines to
      // that of the transform destination (both linear and non-linear)
      // 2. if ras-to-ras transform,
      // then change it to vox-to-vox transform (linear case)

      // modify transform to store inverse also
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


    ////////////////////////////////////////////////////////////////////
    // train gca
    ////////////////////////////////////////////////////////////////////
    // segmentation is seg volume
    // inputs       is the volumes of all inputs
    // transform    is for this subject
    // noint        is whether to use intensity information or not
    GCABtrain(gcab, mri_inputs, mri_seg, transform, target_label) ;
    MRIfree(&mri_seg) ;
    MRIfree(&mri_inputs) ;
    TransformFree(&transform) ;
  }
  GCABcompleteTraining(gcab) ;

  if (smooth > 0) {
    printf("regularizing conditional densities with smooth=%2.2f\n", smooth) ;
    GCAregularizeConditionalDensities(gca, smooth) ;
  }
  if (navgs) {
    printf("applying mean filter %d times to conditional densities\n", navgs) ;
    GCAmeanFilterConditionalDensities(gca, navgs) ;
  }

  printf("writing trained GCAB to %s...\n", out_fname) ;
  if (GCABwrite(gcab, out_fname) != NO_ERROR)
    ErrorExit
      (ERROR_BADFILE, "%s: could not write gca to %s", Progname, out_fname) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRI *mri ;

    mri = GCAbuildMostLikelyVolume(gca, NULL) ;
    MRIwrite(mri, "m.mgz") ;
    MRIfree(&mri) ;
  }

  if (histo_fname) {
    FILE *fp ;
    int   histo_counts[10000], xn, yn, zn, max_count ;
    GCA_NODE  *gcan ;

    memset(histo_counts, 0, sizeof(histo_counts)) ;
    fp = fopen(histo_fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open histo file %s",
                Progname, histo_fname) ;

    max_count = 0 ;
    for (xn = 0 ; xn < gca->node_width;  xn++) {
      for (yn = 0 ; yn < gca->node_height ; yn++) {
        for (zn = 0 ; zn < gca->node_depth ; zn++) {
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
  msec = TimerStop(&start) ;
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
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "SPACING")) {
    spacing = atof(argv[2]) ;
    nargs = 1 ;
    printf("spacing pdfs every %2.1f mm\n", spacing) ;
  } else if (!stricmp(option, "NODE_SPACING")) {
    parms.node_spacing = atof(argv[2]) ;
    nargs = 1 ;
    printf("spacing nodes every %2.1f mm\n", parms.node_spacing) ;
  } else if (!stricmp(option, "BINARIZE")) {
    binarize = 1 ;
    binarize_in = atoi(argv[2]) ;
    binarize_out = atoi(argv[3]) ;
    nargs = 2 ;
    printf("binarizing segmentation values, setting input %d to output %d\n",
           binarize_in, binarize_out) ;
  } else if (!stricmp(option, "NOMRF")) {
    gca_flags |= GCA_NO_MRF ;
    printf("not computing MRF statistics...\n") ;
  } else if (!stricmp(option, "MASK")) {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  } else if (!stricmp(option, "DEBUG_NODE")) {
    Ggca_x = atoi(argv[2]) ;
    Ggca_y = atoi(argv[3]) ;
    Ggca_z = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Ggca_x,Ggca_y,Ggca_z) ;
  } else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  } else if (!stricmp(option, "DEBUG_LABEL")) {
    Ggca_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging label %s (%d)\n", cma_label_to_name(Ggca_label),
           Ggca_label) ;
  } else if (!stricmp(option, "DEBUG_NBR")) {
    Ggca_nbr_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging nbr label %s (%d)\n",
           cma_label_to_name(Ggca_nbr_label), Ggca_nbr_label) ;
  } else if (!stricmp(option, "INSERT")) {
    insert_fname = argv[2] ;
    insert_label = atoi(argv[3]) ;
    nargs = 2 ;
    printf("inserting non-zero vals from %s as label %d...\n",
           insert_fname,insert_label);
  } else if (!stricmp(option, "T1")) {
    strcpy(T1_name, argv[2]) ;
    nargs = 1 ;
    printf("reading T1 data from subject's mri/%s directory\n",
           T1_name) ;
  } else if (!stricmp(option, "PARC_DIR") || !stricmp(option, "SEG_DIR") ||
             !stricmp(option, "SEG") || !stricmp(option, "SEGMENTATION")) {
    seg_dir = argv[2] ;
    nargs = 1 ;
    printf("reading segmentation from subject's mri/%s directory\n",
           seg_dir) ;
  } else if (!stricmp(option, "XFORM")) {
    xform_name = argv[2] ;
    nargs = 1 ;
    printf("reading xform from %s\n", xform_name) ;
  } else if (!stricmp(option, "NOXFORM")) {
    xform_name = NULL ;
    printf("disabling application of xform...\n") ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  } else if (!stricmp(option, "SMOOTH")) {
    smooth = atof(argv[2]) ;
    if (smooth <= 0 || smooth > 1)
      ErrorExit(ERROR_BADPARM,
                "%s: smoothing parameter %2.1f must be in [0,1]\n",
                Progname, smooth) ;
    nargs = 1 ;
    printf("imposing %2.1f smoothing on conditional statistics\n", smooth) ;
  } else switch (toupper(*option)) {
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
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
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
usage_exit(int code) {
  printf("Purpose: %s trains GCA data with (multiple) subject(s)\n", Progname);
  printf("Usage  : %s [options] <subject1> <subject2> ... "
         "<output gca fname>\n",
         Progname) ;
  printf
  ("where SUBJECTS_DIR env variable must be set.\n"
   "Options are:\n"
   "\t-seg dir   - (Required) segmentation volume "
   "(path relative to $subject/mri).\n"
   "\t-xform xform  -  atlas transform (path relative "
   "to $subject/mri/transforms).\n"
   "\t-mask volname   - use volname as a mask "
   "(path relative to $subject/mri.\n"
   "\t-node_spacing   - spacing of classifiers in canonical space\n"
   "\t-prior_spacing  - spacing of class priors in canonical space\n"
   "\t-input name     - specifying training data "
   "(path relative to $subject/mri).\n"
   "                    can specify multiple inputs.  "
   "If not specified, \"orig\" is ued\n"
  );
  exit(code) ;
}

static int input_labels[] = {
                              Left_Cerebral_Exterior,
                              Right_Cerebral_Exterior,
                              Left_Cerebellum_Exterior,
                              Right_Cerebellum_Exterior
                            } ;
// replace the values above with the following
static int output_labels[] = {
                               Left_Cerebral_Cortex,
                               Right_Cerebral_Cortex,
                               Left_Cerebellum_Cortex,
                               Right_Cerebellum_Cortex
                             } ;

static int
replaceLabels(MRI *mri_seg) {
  int    i ;

  for (i = 0 ; i < sizeof(output_labels)/sizeof(output_labels[0]) ; i++)
    MRIreplaceValues(mri_seg, mri_seg, input_labels[i], output_labels[i]) ;
  return(NO_ERROR) ;
}

// modify transform to vox-to-vox
static void modify_transform(TRANSFORM *transform, MRI *mri_inputs, GCA *gca) {
  LTA *lta=0;
  MATRIX *i_to_r=0, *r_to_i=0, *tmpmat=0, *vox2vox;
  MRI *mri_buf = 0;
  GCA_MORPH *gcam = 0;
  static int warned = 0;

  // temp buf to get the transform
  mri_buf = MRIallocHeader(mri_inputs->width,
                           mri_inputs->height,
                           mri_inputs->depth,
                           mri_inputs->type,1);
  MRIcopyHeader(mri_inputs, mri_buf);

  //////////////////////////////////////////////////////////////////////////
  // non-linear transform case
  //////////////////////////////////////////////////////////////////////////
  if (transform->type == MORPH_3D_TYPE) {
    gcam = (GCA_MORPH *) transform->xform;
    if (gcam->atlas.valid) // means it contains the dst volume information
    {
      mri_buf->c_r = gcam->atlas.c_r;
      mri_buf->c_a = gcam->atlas.c_a;
      mri_buf->c_s = gcam->atlas.c_s;
      if (warned == 0) {
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
      if (warned == 0) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stderr, "INFO: modified c_(r,a,s) = 0.\n");
        warned = 1;
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////
  /// linear transform case
  ////////////////////////////////////////////////////////////////////////////
  else if (transform->type == LINEAR_VOX_TO_VOX) {
    lta = (LTA *) (transform->xform);
    // modify using the xform dst
    if (lta->xforms[0].dst.valid) {
      mri_buf->c_r = lta->xforms[0].dst.c_r;
      mri_buf->c_a = lta->xforms[0].dst.c_a;
      mri_buf->c_s = lta->xforms[0].dst.c_s;
      if (warned == 0) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stderr, "INFO: modified c_(r,a,s) using the xform dst.\n");
        warned = 1;
      }
    } else // keep the old behavior
    {
      mri_buf->c_r = 0;
      mri_buf->c_a = 0;
      mri_buf->c_s = 0;
      if (warned == 0) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stderr, "INFO: modified c_(r,a,s) = 0.\n");
        warned = 1;
      }
    }
  } else if (transform->type == LINEAR_RAS_TO_RAS) {
    lta = (LTA *) (transform->xform);
    // modify using the xform dst
    if (lta->xforms[0].dst.valid) {
      mri_buf->c_r = lta->xforms[0].dst.c_r;
      mri_buf->c_a = lta->xforms[0].dst.c_a;
      mri_buf->c_s = lta->xforms[0].dst.c_s;
      if (warned == 0) {
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
      if (warned == 0) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
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
      if (warned == 0) {
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
    i_to_r = extract_i_to_r(mri_inputs);
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
