/**
 * @brief estimates T1 and PD values and transform from a set of FLASH images
 *
 * This program takes an arbitrary # of FLASH images as input, and estimates
 * the T1 and PD values of the data for voxel, as well as a linear transform
 * aligning each of the images. The T1 and PD maps are written into
 * <output directory> together with synthetic volumes names vol?.mgz, one
 * for each of the input volumes. All the output volumes are generated in
 * the common (motion-corrected) space. Note that TR, TE and the flip angle
 * are read directly from the image header. If this information is not
 * available, it can be specified on the command line using -tr <TR in msec>
 * -te <TE in msec> -fa <flip angle in degrees> before each volume.
 * Use -at <xform file name> or -ait <xform file name> to specify
 * transformation for each individual volume. Note only one for each
 * flip-angle is enough. -at will apply the transform to the following
 * volume to align with others.
 */
/*
 * Original Author: Bruce Fischl
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "timer.h"
#include "tags.h"
#include "matrix.h"
#include "transform.h"
#include "version.h"
#include "label.h"
#include "mrinorm.h"
#include "tukey.h"
#include "mrisegment.h"
#include "mriBSpline.h"

static int check_finite(double val)
{
  if (!devFinite(val))
  {
    DiagBreak() ;
    return(0) ;
  }
  return(1) ;
}

//E/ Maybe should update these three to not demand MRI_SHORT -
// but they're never called.
MRI *MRIsadd(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI *MRIsscalarMul(MRI *mri_src, MRI *mri_dst, float scalar) ;
MRI *MRIssqrt(MRI *mri_src, MRI *mri_dst) ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

static int use_tukey = 0 ;
static LTA *Glta = NULL ;
static int invert_flag = 0 ;
static double scale = 1;
static double sigma = 4 ;
static  double base_dt  = 1e-6;
static double momentum = 0.9 ;
static int debug_slice = -1 ;
static int correct_PD = 0 ;
static int synth_flag = 1 ;
static double flip_angle_scale = -1 ;

static double max_T2star = 1000 ;

static int use_brain_mask = 0;
/* compute brain mask and only use brain voxels when computing SSE */

const char *Progname ;

static void usage_exit(int code) ;

static MRI *mri_faf = NULL ;
static int faf_smooth = -1 ;
static float faf_thresh = 0 ;

static int reciprocity = 0 ;
#if 0
static int nfaf = 0 ;        /* # of coefficients in fourier
                                series approximation to flip angle field */
static float *faf_coefs[3][2] ;  /* coefficients (3 spatial dimensions,
                                    and one sin and one cos */
static LABEL *faf_label ;
#endif
static int niter=10 ;
static int conform = 0 ;

static char compress_char = 'z' ;
static int InterpMethod = SAMPLE_TRILINEAR;  /*E* prev default behavior */
static int sinchalfwindow = 3;

static char *residual_name = NULL ;

#define MAX_IMAGES 500

static double FLASHforwardModel(double flip_angle, double TR, double PD,
                                double T1) ;

#if 0
static double faf_coefs_to_scale(double x, double y,  double z,
                                 double w0x, double w0y, double w0z,
                                 float *faf_coefs[3][2],
                                 int nfaf)  ;
static double compute_fa_rms(double T1_wm, double PD_wm,
                             MRI **mri_flash, int nvolumes,
                             float *faf_coefs[3][2],
                             int nfaf, LABEL *faf_label,
                             MRI *mri_T1, MRI *mri_PD)  ;
static int estimate_flip_angle_field(MRI *mri_T1,  MRI *mri_PD,
                                     MRI **mri_flash, int nvolumes, int nfaf,
                                     float *faf_coefs[3][2], LABEL *faf_label,
                                     MRI *mri_faf)  ;
#endif
static double estimate_ms_params(MRI **mri_flash, MRI **mri_flash_synth,
                                 int nvolumes, MRI *mri_T1, MRI *mri_PD,
                                 MRI *mri_sse, MATRIX **M_reg,LTA *lta,
                                 MRI_BSPLINE **mri_flash_bsplines);
static double estimate_ms_params_with_faf(MRI **mri_flash,
    MRI **mri_flash_synth,
    int nvolumes,
    MRI *mri_T1, MRI *mri_PD,
    MRI *mri_sse,
    MATRIX **M_reg,MRI*mri_faf,
    MRI_BSPLINE **mri_flash_bsplines);
#if 0
static double estimate_ms_params_in_label(MRI **mri_flash,
    MRI **mri_flash_synth,
    int nvolumes,
    MRI *mri_T1, MRI *mri_PD,
    MRI *mri_sse, MATRIX **M_reg,
    MRI*mri_faf, LABEL *area);
#endif

#if 0
static double estimate_ms_params_with_kalpha(MRI **mri_flash,
    MRI **mri_flash_synth,
    int nvolumes,
    MRI *mri_T1, MRI *mri_PD,
    MRI *mri_fa, MRI *mri_sse,
    MATRIX **M_reg);
#endif
static void estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_traget,
                                     MATRIX *M_reg, MRI *mri_mask,
                                     MRI_BSPLINE *mri_target_bspline);

static float        tr = 0, te = 0, fa = 0 ;

static int use_outside_reg = 0;
static char *outside_xfm_fname = NULL; //asume all echos need same xform
static int invert_outside_xfm = 0;

static int write_iterations=0;

static int average_volumes_with_different_echo_times(MRI **mri_flash,
    MRI **mri_all_flash,
    int nvolumes_total) ;
static int
average_volumes_with_different_echo_times_and_set_Mreg(MRI **mri_flash,
    MRI **mri_all_flash,
    int nvolumes_total,
    MATRIX **M_reg_orig,
    MATRIX **M_reg);

static MRI *estimate_T2star(MRI **mri_flash, int nvolumes, MRI *mri_PD,
                            MATRIX **Mreg, LTA *lta, MRI *mri_T1,
                            MRI_BSPLINE **mri_all_bsplines) ;
static MRI *compute_T2star_map(MRI **mri_flash, int nvolumes, int *scan_types,
                               MATRIX **Mreg, LTA *lta,
                               MRI_BSPLINE **mri_flash_bsplines) ;

static int findUniqueTETRFA(MRI *mri_flash[], int numvolumes, float *ptr,
                            float *pte, double *pfa);
static int resetTRTEFA(MRI *mri, float tr, float te, double fa);


static int dx = -1, dy, dz, xo, yo, zo ;

static  MRI *mri_mask = NULL;

int
main(int argc, char *argv[])
{
  char   **av, fname[STRLEN] ;
  int    ac, nargs, i ;
  MRI    *mri_T1 = NULL, *mri_PD = NULL, *mri_sse, *mri_T2star ;
  char   *in_fname, *out_dir ;
  int    msec, minutes, seconds, nvolumes, nvolumes_total ;
  Timer start ;
  double rms ;
  float TR=0;
  float TE=0;
  double FA=0;
  int    modified;
  MRI    **mri_flash;
  MRI    **mri_flash_synth;
  MRI    **mri_all_flash;
  MRI_BSPLINE **mri_flash_bsplines = NULL;
  MRI_BSPLINE **mri_all_bsplines = NULL;
  MATRIX **M_reg;
  MATRIX **M_reg_orig;

  /* The following variables are just for finding the brain mask */
  HISTOGRAM *histo;
  MRI_SEGMENTATION *mriseg;
  MRI *mri_tmp = NULL;
  float thresh;
  int b, segno;

  FA = TE = TR = 0;

  std::string cmdline = getAllInfo(argc, argv, "mri_ms_fitparms");

  nargs = handleVersionOption(argc, argv, "mri_ms_fitparms");
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

  if (argc < 2)
  {
    usage_exit(1) ;
  }

  out_dir = argv[argc-1] ;

  mri_flash = (MRI **)malloc(sizeof(mri_flash)*MAX_IMAGES);
  mri_flash_synth = (MRI **)malloc(sizeof(mri_flash)*MAX_IMAGES);
  mri_all_flash = (MRI **)malloc(sizeof(mri_flash)*MAX_IMAGES);
  M_reg = (MATRIX **)malloc(sizeof(mri_flash)*MAX_IMAGES);
  M_reg_orig = (MATRIX **)malloc(sizeof(mri_flash)*MAX_IMAGES);
  mri_flash_bsplines = (MRI_BSPLINE **)malloc(sizeof(mri_flash)*MAX_IMAGES);
  mri_all_bsplines = (MRI_BSPLINE **)malloc(sizeof(mri_flash)*MAX_IMAGES);


  for (i = 0 ; i < MAX_IMAGES ; i++)
  {
    mri_flash[i] = NULL;
    mri_flash_synth[i] = NULL;
    mri_all_flash[i] = NULL;
    M_reg[i] = NULL;
    M_reg_orig[i] = NULL;
    mri_flash_bsplines = NULL;
    mri_all_bsplines = NULL;
  }

/////////////////////////////////////////////////////////////////////
  nvolumes = 0 ;
  for (i = 1 ; i < argc-1 ; i++)
  {
    if (argv[i][0] == '-')
    {
      if (!stricmp(argv[i]+1, "te"))
      {
        te = atof(argv[i+1]) ;
      }
      else if (!stricmp(argv[i]+1, "tr"))
      {
        tr = atof(argv[i+1]) ;
      }
      else if (!stricmp(argv[i]+1, "fa"))
      {
        fa = RADIANS(atof(argv[i+1])) ;
      }
      else if (!stricmp(argv[i]+1, "at"))
      {
        use_outside_reg = 1;
        outside_xfm_fname = argv[i+1];
      }
      else if (!stricmp(argv[i]+1, "ait"))
      {
        use_outside_reg = 1;
        outside_xfm_fname = argv[i+1];
        invert_outside_xfm = 1;
      }
      else
        ErrorExit(ERROR_BADPARM, "%s: unsupported MR parameter %s",
                  Progname, argv[i]+1) ;
      i++ ;  /* skip parameter */
      continue ;
    }

    in_fname = argv[i] ;
    printf("reading %s...", in_fname) ;

    mri_flash[nvolumes] = MRIread(in_fname) ;
    if (mri_flash[nvolumes] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read volume %s",
                Progname, in_fname) ;
    if (mri_mask)
    {
      MRI *mri_tmp ;

      mri_tmp = MRImask(mri_flash[nvolumes], mri_mask, NULL, 0, 0) ;
      MRIfree(&mri_flash[nvolumes]) ;
      mri_flash[nvolumes] = mri_tmp ;
    }
    mri_flash[nvolumes]->dof = 1 ;
    if (dx > 0)   // extract subimage
    {
      MRI *mri_tmp = MRIextract(mri_flash[nvolumes], NULL, xo, yo, zo, dx, dy, dz) ;
      MRIfree(&mri_flash[nvolumes]) ;
      mri_flash[nvolumes] = mri_tmp ;
    }
    if (mri_flash[nvolumes]->type == MRI_UCHAR)
    {
      fprintf(stderr, "WARNING:  input %s is type UCHAR. Make sure it has not been scaled",
	      in_fname) ;
    }

    if (!FEQUAL(scale, 1.0))
    {
      MRIscalarMul(mri_flash[nvolumes], mri_flash[nvolumes], scale) ;
    }
    mri_flash[nvolumes]->register_mat =
      MRIgetVoxelToRasXform(mri_flash[nvolumes]);

    if (flip_angle_scale > 0)
    {
      mri_flash[nvolumes]->flip_angle *= flip_angle_scale ;
    }

    if (tr > 0)
    {
      mri_flash[nvolumes]->tr = tr ;
      //      tr = 0; //reset! without it, the tr will be applied
      // to all subsequent volumes! Well, maybe this is desired, so
      // doesn't need to specify te,tr,fa for each individual echo
    }
    if (te > 0)
    {
      mri_flash[nvolumes]->te = te ;
      //  te = 0;
    }
    if (fa > 0)
    {
      mri_flash[nvolumes]->flip_angle = fa ;
      // fa = 0;
    }
    printf("TE = %2.2f, TR = %2.2f, flip angle = %2.2f\n",
           mri_flash[nvolumes]->te,
           mri_flash[nvolumes]->tr,
           DEGREES(mri_flash[nvolumes]->flip_angle)) ;

    if (outside_xfm_fname != NULL)
    {
      //read in xform and store it in M_reg_orig[i]

      int transform_type =  TransformFileNameType(outside_xfm_fname);
      LTA *lta = 0;

      if (transform_type == MNI_TRANSFORM_TYPE ||
          transform_type == TRANSFORM_ARRAY_TYPE ||
          transform_type == REGISTER_DAT ||
          transform_type == FSLREG_TYPE
         )
      {
        printf("Reading transform ...\n");
        lta = LTAreadEx(outside_xfm_fname) ;
        if (!lta)
          ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
                    Progname, outside_xfm_fname) ;

        if (transform_type == FSLREG_TYPE)
        {
          //currently unsupported, otherwise the usage would be a mess
          // if(mri == 0 || mri_dst == 0){
          //  fprintf(stderr, "ERROR: fslmat does not have information
          // on the src and dst volumes\n");
          //  fprintf(stderr, "ERROR: you must give options '-src' and
          // '-dst' to specify the src and dst volume infos\n");
          // }
          // LTAmodifySrcDstGeom(lta, mri, mri_dst); // add src and
          // dst information
          // LTAchangeType(lta, LINEAR_VOX_TO_VOX);
        }

        if (lta->xforms[0].src.valid == 0 || lta->xforms[0].dst.valid == 0)
        {
          ErrorExit(ERROR_BADPARM,
                    "transform doesn't have valid source or "
                    "destination volume info!");
        }
      }
      else
      {
        ErrorExit(ERROR_BADPARM, "transform is not of MNI, nor "
                  "Register.dat, nor FSLMAT type");
      }

      //change LTA to LINEAR_RAS_TO_RAS type
      if (lta->type != LINEAR_RAS_TO_RAS)
      {
        LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      }

      if (invert_outside_xfm == 0)
      {
        //note need a registration matrix from parameter map to input image
        //thus need to apply inverse when using -at
        M_reg_orig[nvolumes] = MatrixInverse(lta->xforms[0].m_L, NULL) ;
      }
      else
      {
        M_reg_orig[nvolumes] = MatrixCopy(lta->xforms[0].m_L, NULL);
        invert_outside_xfm = 0; //reset
      }
      if (lta)
      {
        LTAfree(&lta);
      }
      outside_xfm_fname = NULL; //clear it
    }

    if (conform)
    {
      MRI *mri_tmp ;

      printf("embedding and interpolating volume\n") ;
      if (InterpMethod!=SAMPLE_TRILINEAR)
        printf("WARNING: interpolation is trilinear in conform!\n") ;
      mri_tmp = MRIconform(mri_flash[nvolumes]) ;
      /*      MRIfree(&mri_src) ;*/
      mri_flash[nvolumes] = mri_tmp ;
    }
    mri_all_flash[nvolumes] =
      mri_flash[nvolumes] ;  /* for multi-echo, mri_flash will
                                be avg across echo time */
    if (FZERO(mri_flash[nvolumes]->tr) ||
        FZERO(mri_flash[nvolumes]->flip_angle))
      ErrorExit(ERROR_BADPARM, "%s: invalid TR or FA for image %d:%s",
                Progname, nvolumes, in_fname) ;

    nvolumes++ ;

    if (nvolumes >= MAX_IMAGES)
    {
      ErrorExit(
        ERROR_BADPARM,
        "ERROR: too many input images!  max of %d!\n",
        MAX_IMAGES);
    }
  }
  ///////////////////////////////////////////////////////////////////////////
  nvolumes_total = nvolumes ;   /* all volumes read in */

  if (faf_smooth > 0)
  {
    MRI *mri_mask, *mri_tmp ;
    if (mri_faf == NULL)
      ErrorExit(ERROR_BADPARM,
                "%s: -fsmooth only applies if FAF image is specified",
                Progname) ;
    mri_mask = MRIbinarize(mri_flash[0], NULL, faf_thresh, 0, CONTROL_MARKED) ;
    mri_tmp = MRIchangeType(mri_mask, MRI_UCHAR, 0, 255, 1) ;
    MRIfree(&mri_mask) ;
    mri_mask = mri_tmp ;
    if (Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_mask, "m.mgz") ;
    }
    mri_tmp = MRIsoapBubbleLabel(mri_faf,
                                 mri_mask,
                                 NULL,
                                 CONTROL_MARKED,
                                 faf_smooth) ;
    MRIwrite(mri_tmp, "faf_smooth.mgz") ;
    MRIfree(&mri_faf) ;
    mri_faf = mri_tmp ;
    MRIfree(&mri_mask) ;
  }

  if (use_outside_reg == 0)
  {
    nvolumes = average_volumes_with_different_echo_times(mri_flash,
               mri_all_flash,
               nvolumes_total) ;
  }
  else
  {
    nvolumes =
      average_volumes_with_different_echo_times_and_set_Mreg(mri_flash,
          mri_all_flash,
          nvolumes_total,
          M_reg_orig,
          M_reg) ;
    niter = 0;
  }

  if (nvolumes <= 2)
  {
    niter = 0 ;
  }  /* don't bother motion-correcting when
                    we only have 2 volumes */
  if (synth_flag > 0)
  {
    for (i = 0 ; i < nvolumes ; i++)
    {
      mri_flash_synth[i] = MRIclone(mri_flash[i], NULL) ;
      mri_flash_synth[i]->register_mat =
        MRIgetVoxelToRasXform(mri_flash[i]);
    }
  }
  else
  {
    for (i = 0 ; i < nvolumes ; i++)
    {
      mri_flash_synth[i] = NULL ;
    }
  }

  {
    int i=0, j;
    for (j = i+1 ; j < nvolumes ; j++)
    {
      if ((mri_flash[i]->width != mri_flash[j]->width) ||
          (mri_flash[i]->height != mri_flash[j]->height) ||
          (mri_flash[i]->depth != mri_flash[j]->depth))
        ErrorExit(ERROR_BADPARM, "%s:\nvolumes %d (type %d) and %d "
                  "(type %d) don't match (%d x %d x %d) vs "
                  "(%d x %d x %d)\n",
                  Progname, i, mri_flash[i]->type, j,
                  mri_flash[j]->type, mri_flash[i]->width,
                  mri_flash[i]->height, mri_flash[i]->depth,
                  mri_flash[j]->width, mri_flash[j]->height,
                  mri_flash[j]->depth) ;
    }

  }

  if (Glta)
  {
    getVolGeom(mri_flash[0], &Glta->xforms[0].src) ;
    LTAchangeType(Glta, LINEAR_VOX_TO_VOX) ;
    if (invert_flag)
    {
      MATRIX *m ;

      printf("inverting transform...\n") ;
      m = MatrixInverse(Glta->xforms[0].m_L, NULL) ;
      MatrixCopy(m, Glta->xforms[0].m_L) ;
      MatrixFree(&m) ;
    }
  }

  
  if (use_outside_reg == 0)
  {
    int j;
    for (j=0; j<nvolumes; j++)
    {
        M_reg[j] = MatrixIdentity(4,(MATRIX *)NULL);
    }
  }

  /* (mr) compute bspline images for interpolation */
  if (InterpMethod == SAMPLE_CUBIC_BSPLINE)
  {
    int j;
    printf("computing flash B-Spline images for cubic interpolation ...\n");
    for (j=0; j<nvolumes; j++)
    {
        mri_flash_bsplines[j] = MRItoBSpline(mri_flash[j], NULL, 3);
    }
  }

  if (nvolumes > 1)
  {
    int j, iter ; /* should exit when M_reg's converge. */
    LTA *lta;

    printf("using %d FLASH volumes (%d total) to estimate tissue parameters.\n",
           nvolumes, nvolumes_total) ;
    for (j = 0 ; j < nvolumes ; j++)
      printf("dof[%d] = %d ... ", j, mri_flash[j]->dof) ;
    printf("\n") ;
    mri_T1 = MRIcloneDifferentType(mri_flash[0], MRI_FLOAT) ;
    mri_PD = MRIcloneDifferentType(mri_flash[0], MRI_FLOAT) ;
    mri_sse = MRIcloneDifferentType(mri_flash[0], MRI_FLOAT) ;


    if (niter == 0)
    {
      if (mri_faf)
        rms =
          estimate_ms_params_with_faf(mri_flash,
                                      synth_flag ? mri_flash_synth : \
                                      NULL, nvolumes, \
                                      mri_T1, mri_PD, mri_sse, \
                                      M_reg, mri_faf, mri_flash_bsplines) ;
      else
        rms = estimate_ms_params(mri_flash, synth_flag ? \
                                 mri_flash_synth : NULL, nvolumes, \
                                 mri_T1, mri_PD, mri_sse, M_reg, NULL, mri_flash_bsplines) ;
      printf("parameter rms = %2.3f\n", rms) ;
    }

    // here is the iterations
    for (iter=0; iter<niter; iter++)
    {
      printf("parameter estimation/motion correction iteration %d of %d\n",
             iter+1, niter) ;
      if (mri_faf)
        rms = estimate_ms_params_with_faf(mri_flash, mri_flash_synth,
                                          nvolumes, mri_T1, mri_PD,
                                          mri_sse, M_reg, mri_faf, mri_flash_bsplines) ;
      else
        rms = estimate_ms_params(mri_flash, mri_flash_synth, nvolumes,
                                 mri_T1, mri_PD, mri_sse, M_reg, NULL, mri_flash_bsplines) ;
      printf("parameter rms = %2.3f\n", rms) ;

      if (use_brain_mask)
      {
        printf("estimate brain mask using PD map \n");
        histo = MRIhistogram(mri_PD, 256);
        for (b=255; b >= 0; b--)
          if (histo->counts[b] > 0)
          {
            break;
          }
        thresh = histo->bins[b]*.25;
        mriseg = MRIsegment(mri_PD, thresh, histo->bins[b]);
        segno = MRIfindMaxSegmentNumber(mriseg);
        mri_tmp = MRIsegmentToImage(mri_PD, NULL, mriseg, segno);
        MRIsegmentFree(&mriseg);
        MRIerode(mri_tmp, mri_tmp);
        mri_mask = MRIbinarize(mri_tmp, NULL, 1, 0, 1);
        mri_tmp = MRIchangeType(mri_mask, MRI_UCHAR, 0,1,1);
        MRIfree(&mri_mask);
        mri_mask =  mri_tmp;
        //  if(iter == (niter -1))
        //   MRIwrite(mri_mask, "brain_mask.mgz") ;
      }

      for (j=0; j<nvolumes; j++)
      {
        printf("estimating rigid alignment for FLASH volume "
               "#%d of %d\n", j+1, nvolumes) ;
        estimate_rigid_regmatrix(mri_flash_synth[j],mri_flash[j],
                                 M_reg[j], mri_mask, mri_flash_bsplines[j]);
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        {
          printf("M_reg[%d]\n",j);
          MatrixPrint(stdout,M_reg[j]);
        }
      }

      if ((write_iterations>0) && (iter%write_iterations==0))
      {
        sprintf(fname,"%s/T1-%d.mg%c",out_dir,iter, compress_char);
        printf("writing T1 estimates to %s...\n", fname) ;
        if (Glta)
        {
          useVolGeomToMRI(&Glta->xforms[0].dst, mri_T1) ;
        }
        MRIaddCommandLine(mri_T1, cmdline) ;
        MRIwrite(mri_T1, fname) ;
        sprintf(fname,"%s/PD-%d.mg%c",out_dir,iter, compress_char);
        printf("writing PD estimates to %s...\n", fname) ;
        if (Glta)
        {
          useVolGeomToMRI(&Glta->xforms[0].dst, mri_PD) ;
        }
        MRIaddCommandLine(mri_PD, cmdline) ;
        MRIwrite(mri_PD, fname) ;
        sprintf(fname,"%s/sse-%d.mg%c",out_dir,iter, compress_char);
        printf("writing residual sse to %s...\n", fname) ;
        if (Glta)
        {
          useVolGeomToMRI(&Glta->xforms[0].dst, mri_sse) ;
        }
        MRIaddCommandLine(mri_sse, cmdline) ;
        MRIwrite(mri_sse, fname) ;

        for (j=0; j<nvolumes; j++)
        {
          sprintf(fname,"%s/vol%d-%d.mg%c",out_dir,j,iter, compress_char);
          printf("writing synthetic images to %s...\n", fname);
          if (Glta)
          {
            useVolGeomToMRI(&Glta->xforms[0].dst, mri_flash_synth[j]) ;
          }
          MRIaddCommandLine(mri_flash_synth[j], cmdline) ;
          MRIwrite(mri_flash_synth[j], fname) ;
          sprintf(fname,"%s/vol%d-%d.lta",out_dir,j,iter);
          printf("writing regisration matrix to %s...\n", fname);
          lta = LTAalloc(1,NULL) ;
          MatrixCopy(M_reg[j],lta->xforms[0].m_L) ;
          // add src and dst information
          getVolGeom(mri_flash_synth[j], &lta->xforms[0].src);
          getVolGeom(mri_flash[j], &lta->xforms[0].dst);
          lta->type = LINEAR_RAS_TO_RAS ;
          LTAwriteEx(lta,fname) ;
        }
      }
    } // iterations end here
    if (Glta != NULL)
      rms = estimate_ms_params(mri_flash, synth_flag ? mri_flash_synth : \
                               NULL, nvolumes, mri_T1, mri_PD, mri_sse, \
                               M_reg, Glta, mri_flash_bsplines) ;
#if 0
    if (nfaf >  0)
    {
      int i ;
      double rms,  pct_change, last_rms ;

      if (!mri_faf)
      {
        mri_faf =
          MRIalloc(mri_T1->width, mri_T1->height, mri_T1->depth,
                   MRI_FLOAT) ;
        if (!mri_faf)
          ErrorExit(ERROR_NOMEMORY, "%s: could  not allocate flip "
                    "angle field image", Progname) ;
        MRIcopyHeader(mri_T1, mri_faf) ;
        MRIsetValues(mri_faf, 1) ;
      }
      printf("estimating flip angle field...\n") ;
      last_rms = estimate_ms_params_in_label(mri_flash, mri_flash_synth,
                                             nvolumes, mri_T1, mri_PD,
                                             mri_sse, M_reg,
                                             NULL,faf_label) ;
      printf("outer iter %03d: parameter rms = %2.3f\n", 0,  last_rms) ;
      for (i = 0  ;  i  < 1000  ; i++)
      {
        estimate_flip_angle_field(mri_T1, mri_PD, mri_flash, nvolumes,
                                  nfaf, faf_coefs, faf_label, mri_faf)  ;
        rms = estimate_ms_params_in_label(mri_flash, mri_flash_synth,
                                          nvolumes, mri_T1, mri_PD,
                                          mri_sse, M_reg, mri_faf,
                                          faf_label) ;
        pct_change = 100.0*(last_rms-rms)/last_rms  ;
        printf("outer iter %03d: parameter rms = %2.3f (%2.3f%%)\n",
               i+1,  rms,  pct_change) ;
        last_rms = rms ;
#if 1
        if (pct_change  < 0.00000001 && i > 5)
        {
          break ;
        }
#endif
      }
      rms = estimate_ms_params_with_faf(mri_flash, mri_flash_synth,
                                        nvolumes, mri_T1, mri_PD,
                                        mri_sse, M_reg, mri_faf) ;
      sprintf(fname,  "%s/faf%d.mg%c", out_dir, iter, compress_char);
      printf("writing flip angle field estimates to %s...\n", fname) ;
      MRIaddCommandLine(mri_faf, cmdline) ;
      MRIwrite(mri_faf, fname);
      MRIfree(&mri_faf);
    }
#endif
    /////////////////////////////////////////
    // write results
    /////////////////////////////////////////
    // modify TR, TE, FA if they are not unique among all input volumes
    modified = findUniqueTETRFA(mri_flash, nvolumes, &TR, &TE, &FA);
    fprintf(stderr, "TR = %.2f, TE = %.2f, Flip_angle = %.2f are used\n",
            TR, TE, DEGREES(FA));

    if (nvolumes > 1)
    {
      resetTRTEFA(mri_PD, TR, TE, FA);
      sprintf(fname,"%s/PD.mg%c",out_dir, compress_char);
      printf("writing PD estimates to %s...\n", fname) ;
      if (Glta)
      {
        useVolGeomToMRI(&Glta->xforms[0].dst, mri_PD) ;
      }
      MRIaddCommandLine(mri_PD, cmdline) ;
      MRIwrite(mri_PD, fname) ;
      resetTRTEFA(mri_T1, TR, TE, FA);
      sprintf(fname,"%s/T1.mg%c",out_dir, compress_char);
      printf("writing T1 estimates to %s...\n", fname) ;
      if (Glta)
      {
        useVolGeomToMRI(&Glta->xforms[0].dst, mri_T1) ;
      }
      MRIaddCommandLine(mri_T1, cmdline) ;
      MRIwrite(mri_T1, fname) ;

      sprintf(fname,"%s/sse.mg%c",out_dir, compress_char);
      printf("writing residual sse to %s...\n", fname) ;
      if (Glta)
      {
        useVolGeomToMRI(&Glta->xforms[0].dst, mri_sse) ;
      }
      MRIaddCommandLine(mri_sse, cmdline) ;
      MRIwrite(mri_sse, fname) ;
      MRIfree(&mri_sse) ;
    }

    if (mri_faf)
    {
      resetTRTEFA(mri_faf, TR, TE, FA);
      sprintf(fname,"%s/faf.mg%c",out_dir, compress_char);
      printf("writing faf map to %s...\n", fname) ;
      MRIaddCommandLine(mri_faf, cmdline) ;
      MRIwrite(mri_faf, fname) ;
    }
    if (synth_flag > 0 && nvolumes > 1)
    {
      for (j=0; j<nvolumes; j++)
      {
#if 0
        sprintf(fname,"%s/vol%d.mg%c",out_dir,j, compress_char);
#else
        sprintf(fname,"%s/flash%d.mg%c",out_dir,
                nint(DEGREES(mri_flash_synth[j]->flip_angle)), compress_char);
#endif
        printf("writing synthetic images to %s...\n", fname);
        if (Glta)
        {
          useVolGeomToMRI(&Glta->xforms[0].dst, mri_flash_synth[j]) ;
        }
        MRIaddCommandLine(mri_flash_synth[j], cmdline) ;
        MRIwrite(mri_flash_synth[j], fname) ;
#if 0
        sprintf(fname,"%s/vol%d.lta",out_dir,j);
#else
        sprintf(fname,"%s/flash%d.lta",out_dir,
                nint(DEGREES(mri_flash_synth[j]->flip_angle)));
#endif
        printf("writing registration matrix to %s...\n", fname);
        lta = LTAalloc(1,NULL) ;
        MatrixCopy(M_reg[j],lta->xforms[0].m_L) ;
        // add src and dst information
        getVolGeom(mri_flash_synth[j], &lta->xforms[0].src);
        getVolGeom(mri_flash[j], &lta->xforms[0].dst);
        lta->type = LINEAR_RAS_TO_RAS ;
        LTAwriteEx(lta,fname) ;
      }
    }
  }
  
  /* (mr) compute bspline images for interpolation */
  if (InterpMethod == SAMPLE_CUBIC_BSPLINE)
  {
    int j;
    printf("computing all B-Spline images for cubic interpolation (for T2star) ...\n");
    for (j=0; j<nvolumes_total; j++)
    {
        mri_all_bsplines[j] = MRItoBSpline(mri_all_flash[j], NULL, 3);
    }
  }
  mri_T2star = estimate_T2star(mri_all_flash, nvolumes_total, mri_PD,
                               M_reg, Glta, mri_T1, mri_all_bsplines) ; 
                               
  if (mri_T1)
    MRIfree(&mri_T1) ;
  if (mri_T2star)
  {
    resetTRTEFA(mri_T2star, TR, TE, FA);
    if  (correct_PD && nvolumes > 1)
    {
      sprintf(fname,"%s/PDcorrected.mg%c",out_dir, compress_char);
      printf("writing corrected PD estimates to %s...\n", fname) ;
      if (Glta)
      {
        useVolGeomToMRI(&Glta->xforms[0].dst, mri_PD) ;
      }
      MRIaddCommandLine(mri_PD, cmdline) ;
      MRIwrite(mri_PD, fname) ;
    }
    sprintf(fname,"%s/T2star.mg%c",out_dir, compress_char);
    printf("writing T2star estimates to %s...\n", fname) ;
    if (Glta)
    {
      useVolGeomToMRI(&Glta->xforms[0].dst, mri_T2star) ;
    }
    MRIaddCommandLine(mri_T2star, cmdline) ;
    MRIwrite(mri_T2star, fname) ;
  }

  for (i = 0 ; i < MAX_IMAGES ; i++)
  {
    if (M_reg[i] != NULL)
    {
      MatrixFree(&M_reg[i]);
    }
    if (M_reg_orig[i] != NULL)
    {
      MatrixFree(&M_reg_orig[i]);
    }
    if (mri_flash_bsplines[i] != NULL)
    {
      MRIfreeBSpline(&mri_flash_bsplines[i]);
    }
    if (mri_all_bsplines[i] != NULL)
    {
      MRIfreeBSpline(&mri_all_bsplines[i]);
    }
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("parameter estimation took %d minutes and %d seconds.\n",
         minutes, seconds) ;
  exit(0);
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int        nargs = 0 ;
  char       *option ;
  TRANSFORM  *transform ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "debug_slice"))
  {
    debug_slice = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging slice %d...\n", debug_slice) ;
  }
  else if (!stricmp(option, "nosynth"))
  {
    synth_flag = 1 ;
    printf("disabling volume synthesis\n") ;
    niter = 0 ;
  }
  else if (!stricmp(option, "reciprocity"))
  {
    reciprocity = 1 ;
    printf("assuming reciprocity in forward modeling (i.e. same coil for transmit and receive\n") ;
  }
  else if (!stricmp(option, "extract"))
  {
    xo = atoi(argv[2]) ;
    yo = atoi(argv[3]) ;
    zo = atoi(argv[4]) ;
    dx = atoi(argv[5]) ;
    dy = atoi(argv[6]) ;
    dz = atoi(argv[7]) ;

    printf("extracting subimages at (%d, %d, %d), size (%d, %d, %d)\n", xo, yo, zo, dx, dy, dz) ;
    nargs = 6 ;
  }
  else if (!stricmp(option, "fa_scale"))
  {
    flip_angle_scale = atof(argv[2]) ;
    nargs = 1 ;
    printf("scaling all flip angles by %2.2f\n", flip_angle_scale) ;
  }
  else if (!stricmp(option, "dt"))
  {
    base_dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting dt = %e\n", base_dt) ;
  }
  else if (!stricmp(option, "tukey"))
  {
    use_tukey = 1 ;
    printf("using tukey biweight of residuals...\n") ;
  }
  else if (!stricmp(option, "max"))
  {
    max_T2star = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting max T2* to %f\n", max_T2star) ;
  }
  else if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)...\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "correct"))
  {
    correct_PD = 1 ;
    printf("correcting PD by T2* estimates\n") ;
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "scale"))
  {
    scale = atof(argv[2]) ;
    nargs = 1 ;
    printf("scaling volumes by %2.3f after reading\n", scale) ;
  }
  else if (!stricmp(option, "window"))
  {
    printf("window option not implemented\n");
    /*E* window_flag = 1 ; */
  }
  else if (!stricmp(option, "nocompress"))
  {
    compress_char = 'h' ;
    printf("not compressing output volumes (saving as .mgh)") ;
    /*E* window_flag = 1 ; */
  }

  /*E* Interpolation method.  Default is trilinear, other options are
    nearest, cubic, sinc.  You can say -foo or -interp foo.  For sinc,
    you can say -interp sinc 3 or -interp sinc -hw 3 or -sinc 3 or
    -sinc -hw 3.  Maybe -hw 3 should imply sinc, but right now it
    doesn't.  */

  else if (!stricmp(option, "st") ||
           !stricmp(option, "sample") ||
           !stricmp(option, "sample_type") ||
           !stricmp(option, "interp"))
  {
    InterpMethod = MRIinterpCode(argv[2]) ;
    nargs = 1;
    if (InterpMethod==SAMPLE_SINC)
    {
      if ((argc<4) ||
          !strncmp(argv[3],"-",1)) /*E* i.e. no sinchalfwindow
                                     value supplied */
      {
        printf("using sinc interpolation (default windowwidth is 6)\n");
      }
      else
      {
        sinchalfwindow = atoi(argv[3]);
        nargs = 2;
        printf("using sinc interpolation with windowwidth of %d\n",
               2*sinchalfwindow);
      }
    }
  }
  else if (!stricmp(option, "sinc"))
  {
    InterpMethod = SAMPLE_SINC;
    if ((argc<3) || !strncmp(argv[2],"-",1)) /*E* i.e. no
                                               sinchalfwindow value
                                               supplied */
    {
      printf("using sinc interpolation (default windowwidth is 6)\n");
    }
    else
    {
      sinchalfwindow = atoi(argv[2]);
      nargs = 1;
      printf("using sinc interpolation with windowwidth of %d\n",
             2*sinchalfwindow);
    }
  }
  else if (!stricmp(option, "sinchalfwindow") ||
           !stricmp(option, "hw"))
  {
    /*E* InterpMethod = SAMPLE_SINC; //? */
    sinchalfwindow = atoi(argv[2]);
    nargs = 1;
    printf("using sinc interpolation with windowwidth of %d\n",
           2*sinchalfwindow);
  }
  else if (!stricmp(option, "trilinear"))
  {
    InterpMethod = SAMPLE_TRILINEAR;
    printf("using trilinear interpolation\n");
  }
  else if (!stricmp(option, "cubic"))
  {
    InterpMethod = SAMPLE_CUBIC_BSPLINE;
    printf("using cubic B-Spline interpolation\n");
  }
  else if (!stricmp(option, "nearest"))
  {
    InterpMethod = SAMPLE_NEAREST;
    printf("using nearest-neighbor interpolation\n");
  }
  else if (!stricmp(option, "tr"))
  {
    tr = atof(argv[2]) ;
    nargs = 1;
  }
  else if (!stricmp(option, "te"))
  {
    te = atof(argv[2]) ;
    nargs = 1;
  }
  else if (!stricmp(option, "fa"))
  {
    fa = RADIANS(atof(argv[2])) ;
    nargs = 1;
  }
  else if (!stricmp(option, "at"))
  {
    use_outside_reg = 1;
    outside_xfm_fname = argv[2];
    nargs = 1;
  }
  else if (!stricmp(option, "ait"))
  {
    use_outside_reg = 1;
    outside_xfm_fname = argv[2];
    invert_outside_xfm = 1;
    nargs = 1;
  }
  else if (!stricmp(option, "-help"))
  {
    usage_exit(0);
  }
  else if (!stricmp(option, "use_brain_mask"))
  {
    use_brain_mask = 1 ;
    printf("Compute a brain mask from PD map and use it when "
           "compute registration\n") ;
  }
  else if (!stricmp(option, "noconform"))
  {
    conform = 0 ;
    printf("inhibiting isotropic volume interpolation\n") ;
  }
  else if (!stricmp(option, "fsmooth"))
  {
    faf_smooth = atoi(argv[2]) ;
    faf_thresh = atof(argv[3]) ;
    printf("smoothing flip angle map for %d iterations of soap "
           "bubble smoothing (thresh=%2.1f)\n",
           faf_smooth, faf_thresh) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "mask"))
  {
    nargs = 1 ;
    printf("reading in mask volume from %s\n", argv[2]) ;
    mri_mask = MRIread(argv[2]) ;
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not read mask volume from %s\n", argv[2]) ;
  }
  else if (!stricmp(option, "afi"))
  {
    double nominal_fa ;

    nargs = 1 ;
    nominal_fa = 60 ;
    printf("using flip angle map %s with nominal value %2.1f degrees\n",
           argv[2],nominal_fa) ;
    mri_faf = MRIread(argv[2]) ;
    if (!mri_faf)
      ErrorExit(ERROR_NOFILE, "%s: could not read flip angle gain "
                "field from %s\n", argv[2]) ;
    if (mri_faf->type != MRI_FLOAT)
    {
      MRI *mri_tmp  ;
      mri_tmp = MRIchangeType(mri_faf, MRI_FLOAT, 0, 1, 1) ;
      MRIcopyHeader(mri_faf, mri_tmp) ;
      MRIfree(&mri_faf) ;
      mri_faf = mri_tmp ;
    }
    MRIscalarMul(mri_faf, mri_faf, 1.0/nominal_fa) ;
  }
  else if (!stricmp(option, "fam"))
  {
    double nominal_fa ;

    nargs = 2 ;
    nominal_fa = atof(argv[3]) ;
    printf("using flip angle map %s with nominal value %2.1f degrees\n",
           argv[2],nominal_fa) ;
    mri_faf = MRIread(argv[2]) ;
    if (!mri_faf)
      ErrorExit(ERROR_NOFILE, "%s: could not read flip angle gain "
                "field from %s\n", argv[2]) ;
    if (mri_faf->type != MRI_FLOAT)
    {
      MRI *mri_tmp  ;
      mri_tmp = MRIchangeType(mri_faf, MRI_FLOAT, 0, 1, 1) ;
      MRIcopyHeader(mri_faf, mri_tmp) ;
      MRIfree(&mri_faf) ;
      mri_faf = mri_tmp ;
    }
    MRIscalarMul(mri_faf, mri_faf, 1.0/nominal_fa) ;
  }
  else if (!stricmp(option, "faf"))
  {
#if 1
    MRI *mri_v, *mri_kernel, *mri_ctrl ;
    double avg_v_size ;

    nargs = 2 ;
    printf("using flip angle map %s with control points specified in %s\n",
           argv[2],argv[3]) ;
    mri_faf = MRIread(argv[2]) ;
    if (!mri_faf)
      ErrorExit(ERROR_NOFILE, "%s: could not read flip angle gain "
                "field from %s\n", argv[2]) ;
    mri_ctrl = MRIread(argv[3]) ;
    if (!mri_ctrl)
      ErrorExit(ERROR_NOFILE, "%s: could not read flip control "
                "point field from %s\n", argv[3]) ;
    MRIbinarize(mri_ctrl, mri_ctrl, 1, 0, 1) ;
    if (mri_ctrl->type != MRI_UCHAR)
    {
      MRI *mri_tmp = MRIchangeType(mri_ctrl, MRI_UCHAR, 0, 255, 1) ;
      MRIfree(&mri_ctrl) ;
      mri_ctrl = mri_tmp ;
    }
    {
      int i ;
      for (i = 0 ; i < 1 ; i++)
      {
        MRIerode(mri_ctrl, mri_ctrl) ;
      }
    }
    mri_v = MRIbuildVoronoiDiagram(mri_faf, mri_ctrl, NULL) ;
#if 0
    MRIsoapBubble(mri_v, mri_ctrl, mri_faf, 25, 1) ;
#else
    avg_v_size = (mri_faf->xsize + mri_faf->ysize + mri_faf->zsize)/3 ;
    mri_kernel = MRIgaussian1d(sigma/avg_v_size, 0) ;
    MRIconvolveGaussian(mri_v, mri_faf, mri_kernel) ;
    MRIfree(&mri_kernel) ;
#endif
    MRIfree(&mri_v) ;
    MRIfree(&mri_ctrl) ;
#else
    int i, j ;

    nfaf = atoi(argv[2]) ;
    nargs = 2 ;
    printf("fitting flip angle field to label %s and removing "
           "distortion\nusing %d coefficient fourier series\n",
           argv[3], nfaf) ;
    faf_label = LabelRead(NULL, argv[3]) ;
    for (i = 0 ; i < 3 ; i++)
      for (j = 0 ; j < 2 ; j++)
      {
        faf_coefs[i][j] = (float *)calloc(nfaf, sizeof(float*)) ;
        if (!faf_coefs[i][j])
          ErrorExit(ERROR_NOMEMORY, "%s: could not allocate FAF "
                    "array %d,%d\n", Progname, i, j) ;
      }
#endif
  }
  else switch (toupper(*option))
    {
    case 'S':
      sigma = atof(argv[2]) ;
      nargs = 1 ;
      printf("smoothing faf field with sigma=%2.3f\n", sigma) ;
      break ;
    case 'N':
      niter = atoi(argv[2]) ;
      printf("performing estimation/motion correction %d times...\n", niter) ;
      nargs = 1 ;
      break ;
    case 'W':
      write_iterations = atoi(argv[2]) ;
      printf("writing out intermediate results every %d iterations...\n",
             write_iterations) ;
      nargs = 1 ;
      break ;
    case 'R':
      residual_name = argv[2] ;
      printf("writing out residuals to %s...\n", residual_name) ;
      nargs = 1 ;
      break ;
    case 'M':
      momentum = atof(argv[2])  ;
      nargs = 1 ;
      printf("setting momentum=%2.2f\n",momentum) ;
      break ;
    case 'T':
    {
      MRI *mri ;
      transform = TransformRead(argv[2]) ;
      if (!transform)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s",
                  Progname, argv[2]) ;

      Glta = (LTA *)(transform->xform) ;
      nargs = 2 ;
      printf("applying transform %s to output volumes...\n", argv[2]) ;
      printf("reading output volume geometry from %s...\n", argv[3]) ;
      mri = MRIread(argv[3]) ;
      if (mri == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read volume %s for "
                  "output volume settings", Progname, argv[3]) ;
      getVolGeom(mri, &Glta->xforms[0].dst) ;
      MRIfree(&mri) ;
    }
    break ;
    case 'I':
      invert_flag = 1 ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
#include "mri_ms_fitparms.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_ms_fitparms_help_xml,mri_ms_fitparms_help_xml_len);
  exit(code) ;
}


static double
FLASHforwardModel(double flip_angle, double TR, double PD, double T1)
{
  double FLASH, E1 ;
  double  CFA, SFA ;


  CFA = cos(flip_angle) ;
  SFA = sin(flip_angle) ;
  E1 = exp(-TR/T1) ;

  FLASH = PD * SFA ;
  if (!DZERO(T1))
  {
    FLASH *= (1-E1)/(1-CFA*E1);
  }
  return(FLASH) ;
}

#define T1_MIN     10.0
#define T1_MAX     50000.0
#define MAX_NVALS  10000
#define MAX_NVOLS  MAX_IMAGES
#define FAk_MAX    2.0
#define FAk_MIN    0.5


static double SignalTableValues[MAX_NVALS][MAX_NVOLS]; // won't fit on stack
static double SignalTableT1[MAX_NVALS], SignalTableNorm[MAX_NVALS];
static double ImageValues[MAX_NVOLS];
static MATRIX *vox2ras[MAX_NVOLS], *ras2vox[MAX_NVOLS];

static double
estimate_ms_params(MRI **mri_flash, MRI **mri_flash_synth, int nvolumes,
                   MRI *mri_T1, MRI *mri_PD, MRI *mri_sse,
                   MATRIX **M_reg, LTA *lta, MRI_BSPLINE** mri_flash_bsplines)
{
  double   total_sse ;
  double   se, best_se, ss, sse, err, val, norm, T1, PD, xf, yf, zf ;
  int      i, j, x, y, z, indx, min_indx, max_indx, best_indx, center_indx,
           stepindx;
  int      width=mri_T1->width, height=mri_T1->height, depth=mri_T1->depth,
    nvalues=MAX_NVALS, nevals, total_dof; 
  int      nstep=11, step[11]= {1024,512,256,128,64,32,16,8,4,2,1};
  MRI      *mri ;
  MATRIX   *voxvec1, *voxvec2, *rasvec1, *rasvec2, *m_xform = NULL ;

  if (lta)
  {
    m_xform = MatrixInverse(lta->xforms[0].m_L, NULL) ;
    printf("using matrix for resampling data:\n") ;
    MatrixPrint(stdout, m_xform) ;
  }

  voxvec1 = MatrixAlloc(4,1,MATRIX_REAL);
  voxvec1->rptr[4][1] = 1.0;
  voxvec2 = MatrixCopy(voxvec1, NULL);
  rasvec1 = MatrixCopy(voxvec1, NULL);
  rasvec2 = MatrixCopy(voxvec1, NULL);

  for (total_dof = j = 0 ; j < nvolumes ; j++)
  {
    vox2ras[j] = MatrixCopy(mri_flash[j]->register_mat, NULL);
    ras2vox[j] = MatrixInverse(vox2ras[j], NULL);
    total_dof += mri_flash[j]->dof ;
  }

  PD = 1;
  for (i=0; i<nvalues; i++)
  {
    T1 = SignalTableT1[i] = T1_MIN+i*(T1_MAX-T1_MIN)/(nvalues-1);
    ss = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash[j] ;
      val =
        SignalTableValues[i][j] =
          FLASHforwardModel(mri->flip_angle, mri->tr, PD, T1) ;
      ss += mri->dof*val*val;
    }
    norm = SignalTableNorm[i] = sqrt(ss);
    if (norm>0)
      for (j = 0 ; j < nvolumes ; j++)
      {
        SignalTableValues[i][j] /= norm;
      }
  }

  total_sse = 0 ;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++)
      {
        ss = 0;
        for (j = 0 ; j < nvolumes ; j++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          mri = mri_flash[j] ;
          voxvec1->data[0]=x;
          voxvec1->data[1]=y;
          voxvec1->data[2]=z;
          if (m_xform)
          {
            MatrixMultiply(m_xform, voxvec1, voxvec2) ;
            MatrixCopy(voxvec2, voxvec1) ;
          }

          MatrixMultiply(vox2ras[j],voxvec1,rasvec1);
          MatrixMultiply(M_reg[j],rasvec1,rasvec2);
          MatrixMultiply(ras2vox[j],rasvec2,voxvec2);
          xf=voxvec2->data[0];
          yf=voxvec2->data[1];
          zf=voxvec2->data[2];
          if (InterpMethod==SAMPLE_SINC)
          {
            MRIsincSampleVolume(mri, xf, yf, zf, sinchalfwindow, &val) ;
          }
          else if (InterpMethod==SAMPLE_CUBIC_BSPLINE)
          {
            MRIsampleBSpline(mri_flash_bsplines[j], xf, yf, zf, 0, &val);
          }
          else
          {
            MRIsampleVolumeType(mri, xf, yf, zf, &val, InterpMethod) ;
          }
	  check_finite(val) ;
	  if (!devFinite(val))
	    DiagBreak() ;
          ImageValues[j] = val;
          ss += mri->dof*val*val;
	  check_finite(ss) ;
        }
        norm = sqrt(ss);
	check_finite(norm) ;
        if (norm>0)
          for (j = 0 ; j < nvolumes ; j++)
          {
            ImageValues[j] /= norm;
          }

        min_indx = best_indx = 0;
        max_indx = nvalues-1;
        best_indx = -1;
        center_indx = -1;
        best_se = 10000000;
        nevals = 0;
        for (stepindx=0; stepindx<nstep; stepindx++)
        {
          for (indx=min_indx; indx<=max_indx; indx+=step[stepindx])
            if (indx!=center_indx)
            {
              se = 0;
              for (j = 0 ; j < nvolumes ; j++)
              {
                err = ImageValues[j]-SignalTableValues[indx][j];
                se += mri_flash[j]->dof*err*err;
              }
              if (se<best_se)
              {
                best_se = se;
                best_indx = indx;
              }
              nevals++;
            }
          min_indx = MAX(best_indx-step[stepindx]/2,1);
          max_indx = MIN(best_indx+step[stepindx]/2,nvalues-1);
          center_indx = best_indx;
        }

	if (center_indx < 0)
	  DiagBreak() ;
        T1 = SignalTableT1[best_indx];
        MRIsetVoxVal(mri_T1, x, y, z, 0, T1);

        PD = norm/SignalTableNorm[best_indx];
#if 0
        if ((short)PD < 0)
        {
          PD = (double)(0x7fff-1) ;
        }
#endif
	if (devFinite(PD) == 0 || devFinite(norm) == 0 || FZERO(SignalTableNorm[best_indx]) || PD > 1e7)
	{
	  printf("PD at (%d, %d, %d) = %f (%f / %f at index %d\n",
		 x,y,z,PD, norm, SignalTableNorm[best_indx], best_indx) ;
	  DiagBreak();
	}
        MRIsetVoxVal(mri_PD, x, y, z, 0, PD);
        if (mri_flash_synth)
        {
          for (j = 0 ; j < nvolumes ; j++)
          {
            mri = mri_flash_synth[j] ;
            MRIsetVoxVal(mri, x, y, z, 0,
                         PD*SignalTableNorm[best_indx]*
                         SignalTableValues[best_indx][j]);
          }
          sse = 0;
          for (j = 0 ; j < nvolumes ; j++)
          {
            err = MRIgetVoxVal(mri_flash_synth[j], x, y, z, 0) -
                  ImageValues[j]*norm;
	    if (!devFinite(err))
	      DiagBreak() ;
            sse += err*err;
	    if (!devFinite(sse))
	      DiagBreak() ;
          }
        }
        else
        {
          sse = 0;
          for (j = 0 ; j < nvolumes ; j++)
          {
            float pred_val ;

            pred_val = PD*SignalTableNorm[best_indx]*
                       SignalTableValues[best_indx][j] ;
            err = pred_val-ImageValues[j]*norm;
	    if (!devFinite(err))
	      DiagBreak() ;
            sse += mri_flash[j]->dof*err*err;
	    if (!devFinite(sse))
	      DiagBreak() ;
          }
        }

        total_sse += (sse/(double)total_dof) ;
	if (!devFinite(total_sse))
	  DiagBreak() ;
        MRIsetVoxVal(mri_sse, x, y, z, 0, sqrt(sse));
	check_finite(sqrt(sse)) ;
        if (T1 >= 4999 && ImageValues[0] > 70 && ImageValues[1] > 70)
        {
          DiagBreak() ;
        }
      }
  if (m_xform)
  {
    MatrixFree(&m_xform) ;
  }
  MatrixFree(&voxvec1) ;
  MatrixFree(&voxvec2) ;
  MatrixFree(&rasvec1) ;
  MatrixFree(&rasvec2) ;
  for (j = 0 ; j < nvolumes ; j++)
  {
    MatrixFree(&vox2ras[j]) ;
    MatrixFree(&ras2vox[j]) ;
  }
  total_sse = sqrt(total_sse / ((double)width*(double)height*(double)depth)) ;
  check_finite(total_sse) ;
  return(total_sse) ;
}

static double
estimate_ms_params_with_faf(MRI **mri_flash, MRI **mri_flash_synth,
                            int nvolumes, MRI *mri_T1, MRI *mri_PD,
                            MRI *mri_sse, MATRIX **M_reg,
                            MRI *mri_faf, MRI_BSPLINE** mri_flash_bsplines)
{
  double   total_sse ;
  double   se, best_se, ss, sse, err, val, norm, T1, PD, xf, yf, zf, inorm ;
  int      j, x, y, z, indx, min_indx, max_indx, best_indx, center_indx, \
  stepindx;
  int      width=mri_T1->width, height=mri_T1->height, depth=mri_T1->depth, \
    nvalues=MAX_NVALS, nevals, total_dof;
  int      nstep=11, step[11]= {1024,512,256,128,64,32,16,8,4,2,1};
  MRI      *mri ;
  MATRIX   *voxvec1, *voxvec2, *rasvec1, *rasvec2, *m_vox2vox;
  double   x0, y0,  z0, w0x, w0y, w0z, faf_scale ;
  VECTOR   *v1, *v2 ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_flash[0], mri_faf) ;
  x0 = width/2 ;
  y0 = height/2 ;
  z0 = depth/2 ;
  w0x = M_PI/x0 ;
  w0y = M_PI/y0 ;
  w0z = M_PI/z0 ;
  voxvec1 = MatrixAlloc(4,1,MATRIX_REAL);
  voxvec1->rptr[4][1] = 1.0;
  voxvec2 = MatrixCopy(voxvec1, NULL);
  rasvec1 = MatrixCopy(voxvec1, NULL);
  rasvec2 = MatrixCopy(voxvec1, NULL);
  for (total_dof = j = 0 ; j < nvolumes ; j++)
  {
    total_dof += mri_flash[j]->dof ;
    vox2ras[j] = MatrixCopy(mri_flash[j]->register_mat, NULL);
    ras2vox[j] = MatrixInverse(vox2ras[j], NULL);
  }

  total_sse = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        ss = 0;
        for (j = 0 ; j < nvolumes ; j++)
        {
          mri = mri_flash[j] ;
          voxvec1->data[0]=x;
          voxvec1->data[1]=y;
          voxvec1->data[2]=z;
          MatrixMultiply(vox2ras[j],voxvec1,rasvec1);
          MatrixMultiply(M_reg[j],rasvec1,rasvec2);
          MatrixMultiply(ras2vox[j],rasvec2,voxvec2);
          xf=voxvec2->data[0];
          yf=voxvec2->data[1];
          zf=voxvec2->data[2];
          if (InterpMethod==SAMPLE_SINC)
          {
            MRIsincSampleVolume(mri, xf, yf, zf, sinchalfwindow, &val) ;
          }
          else if (InterpMethod==SAMPLE_CUBIC_BSPLINE)
          {
            MRIsampleBSpline(mri_flash_bsplines[j], xf, yf, zf, 0, &val);
          }
          else
          {
            MRIsampleVolumeType(mri, xf, yf, zf, &val, InterpMethod) ;
          }
          ImageValues[j] = val;
          ss += mri->dof*val*val;
        }
        inorm = sqrt(ss);
	check_finite(ss) ;
        if (inorm>0)
          for (j = 0 ; j < nvolumes ; j++)
          {
            ImageValues[j] /= inorm;
          }

        min_indx = best_indx = 0;
        max_indx = nvalues-1;
        best_indx = -1;
        center_indx = -1;
        best_se = 10000000;
        nevals = 0;
        V3_X(v1) = x ;
        V3_Y(v1) = y ;
        V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xf = V3_X(v2) ;
        yf = V3_Y(v2) ;
        zf = V3_Z(v2) ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        MRIsampleVolume(mri_faf, xf, yf, zf, &faf_scale) ;
        if (FZERO(faf_scale))
        {
          faf_scale = .1 ;
        }
#if 0
        faf_scale = MRIgetVoxVal(mri_faf,  x, y, z, 0);/*faf_scale = 1 ;*/
#endif
        for (stepindx=0; stepindx<nstep; stepindx++)
        {
          for (indx=min_indx; indx<=max_indx; indx+=step[stepindx])
            if (indx!=center_indx)
            {
              se = 0;

              /* build signal table at  this T1 */
              PD = 1;
              T1 = SignalTableT1[indx] =
                     T1_MIN+indx*(T1_MAX-T1_MIN)/(nvalues-1);
              ss = 0;
              for (j = 0 ; j < nvolumes ; j++)
              {
                mri = mri_flash[j] ;
                val =
                  SignalTableValues[indx][j] =
                    FLASHforwardModel(faf_scale*mri->flip_angle,
                                      mri->tr, PD, T1) ;
                ss += mri->dof*val*val;
              }
              norm = SignalTableNorm[indx] = sqrt(ss);
	      check_finite(norm) ;
              if (norm>0)
                for (j = 0 ; j < nvolumes ; j++)
                {
                  SignalTableValues[indx][j] /= norm;
                }

              for (j = 0 ; j < nvolumes ; j++)
              {
                err = ImageValues[j]-SignalTableValues[indx][j];
                se += mri_flash[j]->dof*err*err;
              }
              if (se<best_se)
              {
                best_se = se;
                best_indx = indx;
              }
              nevals++;
            }
          min_indx = MAX(best_indx-step[stepindx]/2,1);
          max_indx = MIN(best_indx+step[stepindx]/2,nvalues-1);
          center_indx = best_indx;
        }

        T1 = SignalTableT1[best_indx];
        MRIsetVoxVal(mri_T1, x, y, z, 0, T1);
        PD = (inorm/SignalTableNorm[best_indx]);
	if (reciprocity)
	  PD *= faf_scale ;
        if ((mri_PD->type == MRI_SHORT) && ((short)PD < 0))
        {
          PD = (double)(0x7fff-1) ;
        }
        MRIsetVoxVal(mri_PD, x, y, z, 0, PD);
        for (j = 0 ; j < nvolumes ; j++)
        {
          mri = mri_flash_synth[j] ;
          MRIsetVoxVal(mri, x, y, z, 0,
                       FLASHforwardModel(mri->flip_angle,
                                         mri->tr, PD, T1)) ;
        }
        sse = 0;
        for (j = 0 ; j < nvolumes ; j++)
        {
          mri = mri_flash_synth[j] ;
          err = FLASHforwardModel(faf_scale*mri->flip_angle,
                                  mri->tr, PD, T1)-ImageValues[j]*inorm;
          sse += mri_flash[j]->dof*err*err;
	  if (!devFinite(sse))
	    DiagBreak() ;
        }
        total_sse += (sse/(double)total_dof) ;
	if (!devFinite(total_sse))
	  DiagBreak() ;
        MRIsetVoxVal(mri_sse, x, y, z, 0, sqrt(sse/(double)nvolumes));
	check_finite(sqrt(sse/(double)nvolumes)) ;
        if (T1 >= 4999 && ImageValues[0] > 70 && ImageValues[1] > 70)
        {
          DiagBreak() ;
        }
      }
  }
  total_sse = sqrt(total_sse / ((double)width*(double)height*(double)depth*(double)nvolumes)) ;
  check_finite(total_sse) ;
  MatrixFree(&m_vox2vox) ;
  VectorFree(&v1) ;
  VectorFree(&v2) ;
  return(total_sse) ;
}

#if 0
static double
estimate_ms_params_in_label(MRI **mri_flash, MRI **mri_flash_synth,
                            int nvolumes, MRI *mri_T1, MRI *mri_PD,
                            MRI *mri_sse,
                            MATRIX **M_reg, MRI *mri_faf, LABEL *area)
{
  double   total_sse ;
  double   se, best_se, ss, sse, err, val, norm, T1, PD, xf, yf, zf ;
  int      n, i, j, x, y, z, indx, min_indx, max_indx, best_indx, \
  center_indx, stepindx;
  int      width=mri_T1->width, height=mri_T1->height, \
                                       depth=mri_T1->depth, nvalues=MAX_NVALS, nevals;
  int      nstep=11, step[11]= {1024,512,256,128,64,32,16,8,4,2,1};
  MRI      *mri ;
  MATRIX   *voxvec1, *voxvec2, *rasvec1, *rasvec2;
  double   x0, y0,  z0, w0x, w0y, w0z, faf_scale ;

  x0 = width/2 ;
  y0 = height/2 ;
  z0 = depth/2 ;
  w0x = M_PI/x0 ;
  w0y = M_PI/y0 ;
  w0z = M_PI/z0 ;
  voxvec1 = MatrixAlloc(4,1,MATRIX_REAL);
  voxvec1->rptr[4][1] = 1.0;
  voxvec2 = MatrixCopy(voxvec1, NULL);
  rasvec1 = MatrixCopy(voxvec1, NULL);
  rasvec2 = MatrixCopy(voxvec1, NULL);
  for (j = 0 ; j < nvolumes ; j++)
  {
    vox2ras[j] = MatrixCopy(mri_flash[j]->register_mat, NULL);
    ras2vox[j] = MatrixInverse(vox2ras[j], NULL);
  }

  PD = 1;
  for (i=0; i<nvalues; i++)
  {
    T1 = SignalTableT1[i] = T1_MIN+i*(T1_MAX-T1_MIN)/(nvalues-1);
    ss = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash[j] ;
      val = SignalTableValues[i][j] =
              FLASHforwardModel(mri->flip_angle, mri->tr, PD, T1) ;
      ss += val*val;
    }
    norm = SignalTableNorm[i] = sqrt(ss);
    check_finite(norm) ;
    if (norm>0)
      for (j = 0 ; j < nvolumes ; j++)
      {
        SignalTableValues[i][j] /= norm;
      }
  }

  total_sse = 0 ;
  for (n = 0 ; n  < area->n_points ; n++)
  {
    x = area->lv[n].x ;
    y = area->lv[n].y ;
    z = area->lv[n].z ;
    MRIsurfaceRASToVoxel(mri_flash[0], x, y, z, &xf, &yf, &zf)  ;
    x = nint(xf) ;
    y = nint(yf) ;
    z = nint(zf) ;
    if (mri_faf)  /* rebuild signal table  for each point in  space  */
    {
      faf_scale = MRIgetVoxVal(mri_faf,  x, y, z, 0) ;
      PD = 1;
      for (i=0; i<nvalues; i++)
      {
        T1 = SignalTableT1[i] = T1_MIN+i*(T1_MAX-T1_MIN)/(nvalues-1);
        ss = 0;
        for (j = 0 ; j < nvolumes ; j++)
        {
          mri = mri_flash[j] ;
          val =
            SignalTableValues[i][j] =
              FLASHforwardModel(faf_scale*mri->flip_angle,
                                mri->tr, PD, T1) ;
          ss += val*val;
        }
        norm = SignalTableNorm[i] = sqrt(ss);
	check_finite(norm) ;
        if (norm>0)
          for (j = 0 ; j < nvolumes ; j++)
          {
            SignalTableValues[i][j] /= norm;
          }
      }
    }
    else
    {
      faf_scale  = 1 ;
    }

    ss = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash[j] ;
      voxvec1->data[0]=x;
      voxvec1->data[1]=y;
      voxvec1->data[2]=z;
      MatrixMultiply(vox2ras[j],voxvec1,rasvec1);
      MatrixMultiply(M_reg[j],rasvec1,rasvec2);
      MatrixMultiply(ras2vox[j],rasvec2,voxvec2);
      xf=voxvec2->data[0];
      yf=voxvec2->data[1];
      zf=voxvec2->data[2];
      if (InterpMethod==SAMPLE_SINC)
      {
        MRIsincSampleVolume(mri, xf, yf, zf, sinchalfwindow, &val) ;
      }
      else
      {
        MRIsampleVolumeType(mri, xf, yf, zf, &val, InterpMethod) ;
      }
      ImageValues[j] = val;
      ss += val*val;
    }
    norm = sqrt(ss);
    check_finite(norm) ;
    if (norm>0)
      for (j = 0 ; j < nvolumes ; j++)
      {
        ImageValues[j] /= norm;
      }

    min_indx = best_indx = 0;
    max_indx = nvalues-1;
    best_indx = -1;
    center_indx = -1;
    best_se = 10000000;
    nevals = 0;
    for (stepindx=0; stepindx<nstep; stepindx++)
    {
      for (indx=min_indx; indx<=max_indx; indx+=step[stepindx])
        if (indx!=center_indx)
        {
          se = 0;
          for (j = 0 ; j < nvolumes ; j++)
          {
            err = ImageValues[j]-SignalTableValues[indx][j];
            se += err*err;
          }
          if (se<best_se)
          {
            best_se = se;
            best_indx = indx;
          }
          nevals++;
        }
      min_indx = MAX(best_indx-step[stepindx]/2,1);
      max_indx = MIN(best_indx+step[stepindx]/2,nvalues-1);
      center_indx = best_indx;
    }

    T1 = SignalTableT1[best_indx];
    MRIsetVoxVal(mri_T1, x, y, z, 0, T1);

    PD = norm/SignalTableNorm[best_indx];
#if 0
    if ((short)PD < 0)
    {
      PD = (double)(0x7fff-1) ;
    }
#endif
    MRIsetVoxVal(mri_PD, x, y, z, 0, PD);
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash_synth[j] ;
      if  (mri_faf)
        MRIsetVoxVal(mri, x, y, z, 0,
                     FLASHforwardModel(mri->flip_angle, mri->tr, PD, T1));
      else
        MRIsetVoxVal(mri, x, y, z, 0,
                     PD*SignalTableNorm[best_indx]*
                     SignalTableValues[best_indx][j]);
    }
    sse = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash_synth[j] ;
      err = FLASHforwardModel(faf_scale*mri->flip_angle,
                              mri->tr, PD, T1)-ImageValues[j]*norm;
      sse += err*err;
    }
    total_sse += sse ;
    MRIsetVoxVal(mri_sse, x, y, z, 0, sqrt(sse));
    check_finite(sqrt(sse)) ;
    if (T1 >= 4999 && ImageValues[0] > 70 && ImageValues[1] > 70)
    {
      DiagBreak() ;
    }
  }
  total_sse = sqrt(total_sse / (double)(area->n_points*nvolumes)) ;
  check_finite(total_sse) ;
  return(total_sse) ;
}
#endif
#if 0
#define UNITY_FAk_INDEX  nint((1 - FAk_MIN)*(nvalues-1)/(FAk_MAX - FAk_MIN))
static double
estimate_ms_params_with_kalpha(MRI **mri_flash, MRI **mri_flash_synth,
                               int nvolumes, MRI *mri_T1, MRI *mri_PD,
                               MRI *mri_fa, MRI *mri_sse, MATRIX **M_reg)
{
  double   total_sse, \
  SignalTableFAk[MAX_NVALS], FAk, last_sse, last_FAk ;
  double   se, best_se, ss, sse, err, val, norm, \
  T1, PD, xf, yf, zf, last_T1, last_PD ;
  int      i, j, k, x, y, z, T1_indx, min_T1_indx, \
  max_T1_indx, best_T1_indx, center_T1_indx,
               FA_indx, min_FA_indx, max_FA_indx, best_FA_indx, center_FA_indx,
               stepindx;
  int      width=mri_T1->width, height=mri_T1->height, \
                                       depth=mri_T1->depth, nvalues=MAX_NVALS/3,
                                       nevals, niter;
  int      nstep=11, step[11]= {1024,512,256,128,64,32,16,8,4,2,1};
  MRI      *mri ;
  MATRIX   *voxvec1, *voxvec2, *rasvec1, *rasvec2;

  SignalTableValues = (double ***)calloc(nvalues, sizeof(double **)) ;
  SignalTableNorm = (double **)calloc(nvalues, sizeof(double *)) ;
  for (i = 0 ; i < nvalues ; i++)
  {
    SignalTableValues[i] = (double **)calloc(nvalues, sizeof(double *)) ;
    SignalTableNorm[i] = (double *)calloc(nvalues, sizeof(double)) ;
    for (j = 0 ; j < nvalues ; j++)
    {
      SignalTableValues[i][j] =
        (double *)calloc(nvolumes, sizeof(double)) ;
    }
  }

  voxvec1 = MatrixAlloc(4,1,MATRIX_REAL);
  voxvec1->rptr[4][1] = 1.0;
  voxvec2 = MatrixCopy(voxvec1, NULL);
  rasvec1 = MatrixCopy(voxvec1, NULL);
  rasvec2 = MatrixCopy(voxvec1, NULL);
  for (k = 0 ; k < nvolumes ; k++)
  {
    vox2ras[k] = MatrixCopy(mri_flash[k]->register_mat, NULL);
    ras2vox[k] = MatrixInverse(vox2ras[k], NULL);
  }

  PD = 1;
  for (i=0; i<nvalues; i++)
  {
    T1 = SignalTableT1[i] = T1_MIN+i*(T1_MAX-T1_MIN)/(nvalues-1);
    for (j = 0 ; j < nvalues ; j++)
    {
      FAk = SignalTableFAk[j] =
              FAk_MIN+(double)j*(FAk_MAX-FAk_MIN)/(double)(nvalues-1);
      ss = 0;
      for (k = 0 ; k < nvolumes ; k++)
      {
        mri = mri_flash[k] ;
        val = SignalTableValues[i][j][k] =
                FLASHforwardModel(FAk*mri->flip_angle, mri->tr, PD, T1) ;
        ss += val*val;
      }
      norm = SignalTableNorm[i][j] = sqrt(ss);
      check_finite(norm) ;
      if (norm>0)
        for (k = 0 ; k < nvolumes ; k++)
        {
          SignalTableValues[i][j][k] /= norm;
        }
    }
  }

  total_sse = 0 ;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (debug_slice >= 0 && z != debug_slice)
        {
          continue ;
        }
        if (x == 0 && y == 0 && z == 4)
        {
          DiagBreak() ;
        }
        ss = 0;
        for (k = 0 ; k < nvolumes ; k++)
        {
          mri = mri_flash[k] ;
          voxvec1->data[0]=x;
          voxvec1->data[1]=y;
          voxvec1->data[2]=z;
          MatrixMultiply(vox2ras[k],voxvec1,rasvec1);
          MatrixMultiply(M_reg[k],rasvec1,rasvec2);
          MatrixMultiply(ras2vox[k],rasvec2,voxvec2);
          xf=voxvec2->data[0];
          yf=voxvec2->data[1];
          zf=voxvec2->data[2];
          if (InterpMethod==SAMPLE_SINC)
          {
            MRIsincSampleVolume(mri, xf, yf, zf, sinchalfwindow, &val) ;
          }
          else
          {
            MRIsampleVolumeType(mri, xf, yf, zf, &val, InterpMethod) ;
          }
          ImageValues[k] = val;
          ss += val*val;
        }
        if (ImageValues[0] == 76 && ImageValues[1] == 113)
        {
          DiagBreak() ;
        }
        norm = sqrt(ss);
	check_finite(norm) ;
        if (FZERO(norm))
        {
          continue ;  /* no real data */
        }

        if (norm>0)
          for (k = 0 ; k < nvolumes ; k++)
          {
            ImageValues[k] /= norm;
          }

        min_T1_indx = best_T1_indx = 0;
        max_T1_indx = nvalues-1;
        best_T1_indx = -1;
        center_T1_indx = -1;
        min_FA_indx = 0;
        max_FA_indx = nvalues-1;
        best_FA_indx = -1;
        center_FA_indx = -1;
        best_se = 10000000;
        best_FA_indx = UNITY_FAk_INDEX  ;
        nevals = 0;
#if 0
        best_T1_indx = nint((1000.0 - T1_MIN) *
                            (nvalues-1) / (T1_MAX - T1_MIN)) ;
#endif

        sse = -1 ;
        niter = 0 ;
        last_T1 = last_FAk = last_PD = 0 ;  /* for compiler warning */
        do
        {
          for (stepindx=0; stepindx<nstep; stepindx++)
          {
            for (T1_indx=min_T1_indx;
                 T1_indx<=max_T1_indx;
                 T1_indx+=step[stepindx])
            {
              if (T1_indx!=center_T1_indx)
              {
                se = 0;
                for (k = 0 ; k < nvolumes ; k++)
                {
                  err = ImageValues[k]-
                        SignalTableValues[T1_indx][best_FA_indx][k];
                  se += err*err;
                }
                if (se<best_se)
                {
                  best_se = se;
                  best_T1_indx = T1_indx;
                }
                nevals++;
              }
            }
            min_T1_indx = MAX(best_T1_indx-step[stepindx]/2,1);
            max_T1_indx = MIN(best_T1_indx+step[stepindx]/2,nvalues-1);
            center_T1_indx = best_T1_indx;
          }
          for (stepindx=0; stepindx<nstep; stepindx++)
          {
            for (FA_indx=min_FA_indx;
                 FA_indx<=max_FA_indx;
                 FA_indx+=step[stepindx])
            {
              if (FA_indx!=center_FA_indx)
              {
                se = 0;
                for (k = 0 ; k < nvolumes ; k++)
                {
                  err = ImageValues[k]-
                        SignalTableValues[best_T1_indx][FA_indx][k];
                  se += err*err;
                }
                if (se<best_se)
                {
                  best_se = se;
                  best_FA_indx = FA_indx;
                }
                nevals++;
              }
            }
            min_FA_indx = MAX(best_FA_indx-step[stepindx]/2,1);
            max_FA_indx = MIN(best_FA_indx+step[stepindx]/2,nvalues-1);
            center_FA_indx = best_FA_indx;
          }

          T1 = SignalTableT1[best_T1_indx];
          MRIsetVoxVal(mri_T1, x, y, z, 0, T1);

          FAk = MRIFvox(mri_fa, x, y, z) = SignalTableFAk[best_FA_indx];
          PD = norm/SignalTableNorm[best_T1_indx][best_FA_indx];
#if 0
          if ((short)PD < 0)
          {
            PD = (double)(0x7fff-1) ;
          }
#endif
          MRIsetVoxVal(mri_PD, x, y, z, 0, PD);
          for (k = 0 ; k < nvolumes ; k++)
          {
            mri = mri_flash_synth[k] ;
            MRIsetVoxVal
            (mri, x, y, z, 0,
             PD *
             SignalTableNorm[best_T1_indx][best_FA_indx] *
             SignalTableValues[best_T1_indx][best_FA_indx][k]);
          }
          last_sse = sse ;
          sse = 0 ;
          for (k = 0 ; k < nvolumes ; k++)
          {
            err =
              MRIgetVoxVal(mri_flash_synth[k], x, y, z, 0) -
              ImageValues[k]*norm;
            sse += err*err;
          }
          if (last_sse < 0)
          {
            last_sse = sse+1 ;  /* first time */
          }
          MRIsetVoxVal(mri_sse, x, y, z, 0, sqrt(sse));
	  check_finite(sqrt(sse)) ;
          if (sse > last_sse)  /* revert to old values */
          {
            T1 = last_T1 ;
            FAk = last_FAk ;
            PD = last_PD ;
            best_T1_indx = nint((T1 - T1_MIN) *
                                (nvalues-1) / (T1_MAX - T1_MIN)) ;
            best_FA_indx = nint((FAk - FAk_MIN) *
                                (nvalues-1) / (FAk_MAX - FAk_MIN)) ;
            MRIsetVoxVal(mri_PD, x, y, z, 0, PD) ;
            for (k = 0 ; k < nvolumes ; k++)
            {
              mri = mri_flash_synth[k] ;
              MRIsetVoxVal
              (mri, x, y, z, 0,
               PD * SignalTableNorm[best_T1_indx][best_FA_indx]
               * SignalTableValues[best_T1_indx][best_FA_indx][k]);
            }
            sse = last_sse ;
          }
          else
          {
            last_T1 = T1 ;
            last_FAk = FAk ;
            last_PD = PD ;
          }

        }
        while ((sse < last_sse) && (niter++ < 4));
        total_sse += sse ;
      }
  total_sse = sqrt(total_sse / ((double)width*(double)height*(double)depth)) ;
  check_finite(total_sse) ;
  for (i = 0 ; i < nvalues ; i++)
  {
    for (j = 0 ; j < nvalues ; j++)
    {
      free(SignalTableValues[i][j]) ;
    }
    free(SignalTableValues[i]) ;
    free(SignalTableNorm[i]) ;
  }
  free(SignalTableValues) ;
  free(SignalTableNorm) ;
  return(total_sse) ;
}
#endif

MRI *
MRIsadd(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  short   *p1, *p2, *pdst ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = &MRISvox(mri1, 0, y, z) ;
      p2 = &MRISvox(mri2, 0, y, z) ;
      pdst = &MRISvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        *pdst++ = *p1++ + *p2++ ;
      }
    }
  }
  return(mri_dst) ;
}

static void
estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_target,
                         MATRIX *M_reg, MRI *mri_mask, MRI_BSPLINE *mri_target_bspline)
{
  double   xf, yf, zf, tx, ty, tz, ax, ay, az, ca, sa, \
  val1, val2, err, sse, best_sse, dt=0.1, da=RADIANS(0.025), tol=0.00001;
  int      x, y, z, txi, tyi, tzi, axi, ayi, azi, \
  indx, stepindx, changed, pass;
  int      width=mri_source->width, height=mri_source->height, \
                 depth=mri_source->depth, dx=10, dy=10, dz=10, nvalues;
#if 1
  /*
    int      nstep=10, step[10]={512, 256, 128, 64, 32,16,8,4,2,1}, scale;
  */
  int      nstep=8, step[8]= {128, 64, 32, 16,8,4,2,1}, scale;
#else
  int      nstep=1, step[1]= {1}, scale;
#endif
  MATRIX   *vox2ras_source, *ras2vox_source, *vox2ras_target, \
  *ras2vox_target, *vox_s2vox_t;
  MATRIX   *M_reg_bak, *M_reg_opt, *M_tmp, *M_delta, \
  *M_delta1, *M_delta2, *M_delta3, *M_delta4, *M_delta5, *M_delta6;
  MATRIX   *voxmat1, *voxmat2;
  double   *voxval1, *voxval2, tukey_thresh = 100 ;

  vox2ras_source = MatrixCopy(mri_source->register_mat, NULL);
  vox2ras_target = MatrixCopy(mri_target->register_mat, NULL);
  ras2vox_source = MatrixInverse(vox2ras_source, NULL);
  ras2vox_target = MatrixInverse(vox2ras_target, NULL);
  vox_s2vox_t = MatrixIdentity(4,NULL);

  if (use_brain_mask && mri_mask != NULL)
  {
    /* Only count voxels fall within brain mask
       so reduce sampling steps
    */
    dx = 5;
    dy = 5;
    dz = 5;
  }

  nvalues = 0;

  if (use_brain_mask && mri_mask != NULL)
  {
    /* Only count voxels fall within brain mask
       so reduce sampling steps
    */

    for (z = 0 ; z < depth ; z++)
      for (y = 0 ; y < height ; y++)
        for (x = 0 ; x < width ; x++)
        {
          if ((x%dx==0)&&(y%dy==0)&&(z%dz==0))
          {
            if (MRIvox(mri_mask, x, y, z) > 0)
            {
              nvalues++;
            }
          }
        }

  }
  else
  {
    for (z = 0 ; z < depth ; z++)
      for (y = 0 ; y < height ; y++)
        for (x = 0 ; x < width ; x++)
        {
          if ((x%dx==0)&&(y%dy==0)&&(z%dz==0))
          {
            nvalues++;
          }
        }
  }

  voxmat1 = MatrixAlloc(4,nvalues,MATRIX_REAL);
  voxmat2 = MatrixCopy(voxmat1, NULL);
  voxval1 = (double*)calloc(nvalues+1, sizeof(double*)) ;
  voxval2 = (double*)calloc(nvalues+1, sizeof(double*)) ;
  if (!voxmat1 || !voxmat2 || !voxval1 || !voxval2)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate %d value arrays", Progname, nvalues) ;

  indx = 0;

  if (use_brain_mask && mri_mask != NULL)
  {
    for (z = 0 ; z < depth ; z++)
      for (y = 0 ; y < height ; y++)
        for (x = 0 ; x < width ; x++)
        {
          if ((x%dx==0)&&(y%dy==0)&&(z%dz==0) &&(MRIvox(mri_mask,x,y,z) > 0))
          {
            indx++;
            voxmat1->rptr[1][indx] = x;
            voxmat1->rptr[2][indx] = y;
            voxmat1->rptr[3][indx] = z;
            voxmat1->rptr[4][indx] = 1;
            voxval1[indx] = MRIgetVoxVal(mri_source, x, y, z, 0);
          }
        }

  }
  else
  {
    for (z = 0 ; z < depth ; z++)
      for (y = 0 ; y < height ; y++)
        for (x = 0 ; x < width ; x++)
        {
          if ((x%dx==0)&&(y%dy==0)&&(z%dz==0))
          {
            indx++;
            voxmat1->rptr[1][indx] = x;
            voxmat1->rptr[2][indx] = y;
            voxmat1->rptr[3][indx] = z;
            voxmat1->rptr[4][indx] = 1;
            voxval1[indx] = MRIgetVoxVal(mri_source, x, y, z, 0);
          }
        }
  }

  M_delta = MatrixIdentity(4,NULL);
  M_delta1 = MatrixIdentity(4,NULL);
  M_delta2 = MatrixIdentity(4,NULL);
  M_delta3 = MatrixIdentity(4,NULL);
  M_delta4 = MatrixIdentity(4,NULL);
  M_delta5 = MatrixIdentity(4,NULL);
  M_delta6 = MatrixIdentity(4,NULL);

  M_reg_opt = MatrixCopy(M_reg, NULL);
  M_reg_bak = MatrixCopy(M_reg_opt, NULL);
  M_tmp = MatrixCopy(M_reg, NULL);

  if (use_tukey)
  {
    double total_error ;

    MatrixMultiply(M_reg,vox2ras_source,M_tmp);
    MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);

    MatrixMultiply(vox_s2vox_t,voxmat1,voxmat2);
    for (total_error = 0.0, indx=1; indx<=nvalues; indx++)
    {
      xf=voxmat2->rptr[1][indx];
      yf=voxmat2->rptr[2][indx];
      zf=voxmat2->rptr[3][indx];
      if (InterpMethod==SAMPLE_SINC)
        MRIsincSampleVolume(mri_target, xf, yf, zf,
                            sinchalfwindow, &val2) ;
      else if (InterpMethod==SAMPLE_CUBIC_BSPLINE)
      {
          MRIsampleBSpline(mri_target_bspline, xf, yf, zf, 0, &val2);
      }
      else
      {
        MRIsampleVolumeType(mri_target, xf, yf, zf, &val2, InterpMethod) ;
      }
      val1 = voxval1[indx];
      err = val1-val2;
      total_error += abs(err) ;
    }
    tukey_thresh = total_error / (float)nvalues ;
    printf("setting tukey threshold to %2.3f\n", tukey_thresh) ;
  }

  best_sse = 10000000;
  for (stepindx=0; stepindx<nstep; stepindx++)
  {
    scale = step[stepindx];
    changed = 1;
    pass = 0;
    while (changed)
    {
      pass++;
      changed = 0;
      MatrixCopy(M_reg_opt, M_reg_bak);
      for (txi = -1; txi <= 1; txi++)
        for (tyi = -1; tyi <= 1; tyi++)
          for (tzi = -1; tzi <= 1; tzi++)
            for (axi = -1; axi <= 1; axi++)
              for (ayi = -1; ayi <= 1; ayi++)
                for (azi = -1; azi <= 1; azi++)
                {
                  tx = txi*dt*scale;
                  ty = tyi*dt*scale;
                  tz = tzi*dt*scale;
                  ax = axi*da*scale;
                  ay = ayi*da*scale;
                  az = azi*da*scale;
                  M_delta1->rptr[1][4]=tx;
                  M_delta1->rptr[2][4]=ty;
                  M_delta1->rptr[3][4]=tz;
                  ca = cos(ax);
                  sa = sin(ax);
                  M_delta2->rptr[2][2]=ca;
                  M_delta2->rptr[2][3]=-sa;
                  M_delta2->rptr[3][2]=sa;
                  M_delta2->rptr[3][3]=ca;
                  MatrixMultiply(M_delta2,M_delta1,M_delta5);
                  ca = cos(ay);
                  sa = sin(ay);
                  M_delta3->rptr[1][1]=ca;
                  M_delta3->rptr[1][3]=-sa;
                  M_delta3->rptr[3][1]=sa;
                  M_delta3->rptr[3][3]=ca;
                  MatrixMultiply(M_delta3,M_delta5,M_delta6);
                  ca = cos(az);
                  sa = sin(az);
                  M_delta4->rptr[1][1]=ca;
                  M_delta4->rptr[1][2]=-sa;
                  M_delta4->rptr[2][1]=sa;
                  M_delta4->rptr[2][2]=ca;
                  MatrixMultiply(M_delta4,M_delta6,M_delta);
                  MatrixMultiply(M_delta,M_reg_bak,M_reg);

                  MatrixMultiply(M_reg,vox2ras_source,M_tmp);
                  MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);

                  MatrixMultiply(vox_s2vox_t,voxmat1,voxmat2);
                  sse = 0;
                  for (indx=1; indx<=nvalues; indx++)
                  {
                    xf=voxmat2->rptr[1][indx];
                    yf=voxmat2->rptr[2][indx];
                    zf=voxmat2->rptr[3][indx];
                    if (InterpMethod==SAMPLE_SINC)
                      MRIsincSampleVolume(mri_target, xf, yf, zf,
                                          sinchalfwindow, &val2) ;
      							else if (InterpMethod==SAMPLE_CUBIC_BSPLINE)
      							{
          							MRIsampleBSpline(mri_target_bspline, xf, yf, zf, 0, &val2);
      							}
                    else
                      MRIsampleVolumeType(mri_target, xf, yf, zf,
                                          &val2, InterpMethod) ;
                    voxval2[indx] = val2;
                    val1 = voxval1[indx];
                    if (use_tukey)
                    {
                      err = tukey_biweight(val1-val2, tukey_thresh);
                    }
                    else
                    {
                      err = val1-val2;
                    }
                    sse += err*err;
                  }
                  sse /= nvalues;
                  if (sse<best_sse-tol)
                  {
                    best_sse = sse;
                    MatrixCopy(M_reg, M_reg_opt);
                    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
                      printf("%d (%d) %f %f %f %f %f %f sse = "
                             "%f (%f)\n",scale,pass,tx,ty,tz,ax,
                             ay,az,sse,sqrt(sse));
                    /*
                      printf("M_delta\n"); MatrixPrint(stdout,M_delta);
                    */
                    /*
                      printf("M_delta1\n");
                      MatrixPrint(stdout,M_delta1);
                      printf("M_delta2\n");
                      MatrixPrint(stdout,M_delta2);
                      printf("M_delta3\n");
                      MatrixPrint(stdout,M_delta3);
                      printf("M_delta4\n");
                      MatrixPrint(stdout,M_delta4);
                      printf("M_delta5\n");
                      MatrixPrint(stdout,M_delta5);
                      printf("M_delta6\n");
                      MatrixPrint(stdout,M_delta6);
                    */
                    /*
                      printf("M_reg\n");
                      MatrixPrint(stdout,M_reg);
                      printf("vox_s2vox_t\n");
                      MatrixPrint(stdout,vox_s2vox_t);
                    */
                    changed = 1;
                  }
                }

    }
    printf("step %d (dx=%2.1f, da=%2.1f): sse = %f (%f)\n",
           stepindx,tx, DEGREES(ax), sse,sqrt(sse));
  }
  MatrixCopy(M_reg_opt, M_reg);
}

MRI *
MRIssqrt(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, frame ;
  short   *psrc, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri_src->type)
        {
        case MRI_SHORT:
          psrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
	    check_finite(sqrt(*psrc)) ;
            *pdst++ = sqrt(*psrc++) ;
          }
          break ;
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED,
                       "MRIssqrt: unsupported type %d",
                       mri_src->type)) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIsscalarMul(MRI *mri_src, MRI *mri_dst, float scalar)
{
  int     width, height, depth, x, y, z, frame ;
  short   *psrc, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri_src->type)
        {
        case MRI_SHORT:
          psrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            *pdst++ = *psrc++ * scalar ;
          }
          break ;
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED,
                       "MRIsscalarMul: unsupported type %d",
                       mri_src->type)) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
double
dFLASH_dk(MRI *mri_T1, MRI *mri_PD, MRI *mri_fa,
          double TR, double flip_angle, int x, int y, int z)
{
  double T1, PD, dk, k, e_TR_T1_minus_1, cos_ka, numer, denom, e_TR_T1 ;

  T1 = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
  PD = MRIgetVoxVal(mri_PD, x, y, z, 0) ;
  k = MRIFvox(mri_fa, x, y, z) ;

  e_TR_T1 = exp(TR/T1) ;
  e_TR_T1_minus_1 = e_TR_T1 - 1 ;
  cos_ka = cos(k*flip_angle) ;

  numer = PD * flip_angle * e_TR_T1_minus_1 * (e_TR_T1*cos_ka-1) ;
  denom = (e_TR_T1-cos_ka) ;
  denom*= denom ;
  if (FZERO(denom))
  {
    denom = 1 ;
  }
  dk = numer /denom ;
  return(dk) ;
}
#define PARAMETERS_MATCH(mri1, mri2) ((mri1->tr == mri2->tr) && (mri1->flip_angle == mri2->flip_angle))

static int
average_volumes_with_different_echo_times(MRI **mri_flash,
    MRI **mri_all_flash,
    int nvolumes_total)
{
  int i, j, nvolumes, averaged[MAX_IMAGES], navgs, all_tes_equal ;
  MRI    *mri_avg ;

  memset(averaged, 0, sizeof(averaged)) ;

  for (all_tes_equal = 1, i = 1 ; i < nvolumes_total ; i++)
  {
    if (mri_all_flash[i]->te != mri_all_flash[0]->te)
    {
      all_tes_equal = 0 ;
      break ;
    }
  }
  if (all_tes_equal)
  {
    for (i = 0 ; i < nvolumes_total ; i++)
    {
      mri_flash[i] = mri_all_flash[i] ;
    }
    return(nvolumes_total) ;
  }

  for (nvolumes = i = 0 ; i < nvolumes_total ; i++)
  {
    if (averaged[i])
    {
      continue ;
    }
    averaged[i] = 1 ;
    mri_avg = MRIcopy(mri_all_flash[i], NULL) ;
    mri_avg->register_mat = MRIgetVoxelToRasXform(mri_all_flash[nvolumes]);
    mri_avg->dof = 1 ;
    navgs = 1 ;
    for (j = i+1 ; j < nvolumes_total ; j++)
    {
      if (averaged[j])
      {
        continue ;
      }
      if (PARAMETERS_MATCH(mri_all_flash[j], mri_avg) == 0)
      {
        continue ;
      }
      MRIaverage(mri_all_flash[j], navgs, mri_avg) ;
      averaged[j] = 1 ;
      navgs++ ;
    }
    mri_flash[nvolumes] = mri_avg ;
    nvolumes++ ;
  }

  return(nvolumes) ;
}

static int
average_volumes_with_different_echo_times_and_set_Mreg(MRI **mri_flash,
    MRI **mri_all_flash,
    int nvolumes_total,
    MATRIX **M_reg_orig,
    MATRIX **M_reg)
{
  int i, j, nvolumes, averaged[MAX_IMAGES], navgs, all_tes_equal ;
  MRI    *mri_avg ;

  memset(averaged, 0, sizeof(averaged)) ;

  for (all_tes_equal = 1, i = 1 ; i < nvolumes_total ; i++)
  {
    if (mri_all_flash[i]->te != mri_all_flash[0]->te)
    {
      all_tes_equal = 0 ;
      break ;
    }
  }
  if (all_tes_equal)
  {
    for (i = 0 ; i < nvolumes_total ; i++)
    {
      mri_flash[i] = mri_all_flash[i] ;
      if (M_reg_orig[i] == NULL)
      {
        M_reg[i] = MatrixIdentity(4, (MATRIX *)NULL);
      }
      else
      {
        M_reg[i] = MatrixCopy(M_reg_orig[i], (MATRIX *)NULL);
      }
    }
    return(nvolumes_total) ;
  }

  for (nvolumes = i = 0 ; i < nvolumes_total ; i++)
  {
    if (averaged[i])
    {
      continue ;
    }
    averaged[i] = 1 ;
    mri_avg = MRIcopy(mri_all_flash[i], NULL) ;
    if (M_reg_orig[i] == NULL)
    {
      M_reg[nvolumes] = MatrixIdentity(4, (MATRIX *)NULL);
    }
    else
    {
      M_reg[nvolumes] = MatrixCopy(M_reg_orig[i], (MATRIX *)NULL);
    }

    mri_avg->register_mat = MRIgetVoxelToRasXform(mri_all_flash[nvolumes]);
    mri_avg->dof = 1 ;
    navgs = 1 ;
    for (j = i+1 ; j < nvolumes_total ; j++)
    {
      if (averaged[j])
      {
        continue ;
      }
      if (PARAMETERS_MATCH(mri_all_flash[j], mri_avg) == 0)
      {
        continue ;
      }
      MRIaverage(mri_all_flash[j], navgs, mri_avg) ;

      //the following shouldn't happen, user should assign xform
      //to the first echo
      if (M_reg_orig[j] != NULL)
      {
        M_reg[nvolumes] = MatrixCopy(M_reg_orig[j], M_reg[nvolumes]);
      }

      averaged[j] = 1 ;
      navgs++ ;
    }
    mri_flash[nvolumes] = mri_avg ;
    nvolumes++ ;
  }


  return(nvolumes) ;
}

static MRI *
estimate_T2star(MRI **mri_flash, int nvolumes, MRI *mri_PD,
                MATRIX **Mreg, LTA *lta, MRI *mri_T1, MRI_BSPLINE **mri_flash_bsplines)
{
  int    i, j, nechoes, processed[MAX_IMAGES], nprocessed, x, y,\
  z, different_te, width, depth, height, unique_te ;
  MRI    *mri_T2star = NULL ;
  double T2star ;
  double PD = 10 ;

  /* first decide whether T2* can be estimated at all */
  different_te = 0 ;
  for (i = 0 ; i < nvolumes ; i++)
  {
    unique_te = 1 ;
    for (j = 0 ; j < i ; j++)
    {
      if (mri_flash[i]->te == mri_flash[j]->te)
      {
        unique_te = 0 ;
        break ;
      }
    }
    if (unique_te)
    {
      different_te++ ;
    }
  }
  if (different_te <= 1)
  {
    return(NULL) ;  /* can't estimate T2* */
  }

  memset(processed, 0, sizeof(processed)) ;

  for (nprocessed = nechoes = i = 0 ; i < nvolumes ; i++)
  {
    if (processed[i])
    {
      continue ;
    }
    processed[i] = nprocessed+1 ;
    nechoes = 1 ;

    for (j = i+1 ; j < nvolumes ; j++)
    {
      if (processed[j])
      {
        continue ;
      }
      if (PARAMETERS_MATCH(mri_flash[i], mri_flash[j]) == 0)
      {
        continue ;
      }
      processed[j] = nprocessed+1 ;
    }
    nprocessed++ ;
  }
  printf("estimating T2* with %d different acquisitions, "
         "each with %d echoes...\n",
         nvolumes/different_te, different_te) ;
  mri_T2star = compute_T2star_map(mri_flash, nvolumes, processed, Mreg, lta, mri_flash_bsplines) ;

  /* now update PD map to take out T2* component */
  if (correct_PD)
  {
    VECTOR *v1, *v2 ;
    MATRIX *m_vox2vox = NULL ;
    double T1, faf_scale, xf, yf, zf, S, K ;
    MRI    *mri ;
    int    i ;

    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
    if (mri_faf)
    {
      m_vox2vox = MRIgetVoxelToVoxelXform(mri_flash[0], mri_faf) ;
    }
    width = mri_T2star->width ;
    height = mri_T2star->height ;
    depth = mri_T2star->depth ;
    for (x = 0 ; x < width ; x++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (z = 0 ; z < depth ; z++)
        {
          if (mri_PD)
          {
            PD = MRIgetVoxVal(mri_PD, x, y, z, 0);
            if (PD < 100)
            {
              MRIsetVoxVal(mri_T2star, x, y, z, 0, 0) ;
            }
          }
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          T2star = MRIgetVoxVal(mri_T2star, x, y, z, 0) ;
          if (FZERO(T2star))
          {
            continue ;
          }

          T1 = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
          if (mri_faf == NULL)
          {
            faf_scale = 1.0 ;
          }
          else
          {
            V3_X(v1) = x ;
            V3_Y(v1) = y ;
            V3_Z(v1) = z ;
            MatrixMultiply(m_vox2vox, v1, v2) ;
            xf = V3_X(v2) ;
            yf = V3_Y(v2) ;
            zf = V3_Z(v2) ;
            MRIsampleVolume(mri_faf, xf, yf, zf, &faf_scale) ;
          }
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          PD = 0.0 ;
          for (i = 0 ; i < nvolumes ; i++)
          {
            mri = mri_flash[i] ;
            K = FLASHforwardModel(faf_scale * mri->flip_angle, mri->tr, 1, T1)
                * exp(-mri->te/T2star) ;
            S = MRIgetVoxVal(mri, x, y, z, 0) ;
            if (DZERO(K))
            {
              continue ;
            }
            PD += S / K ;
          }
          PD /= nvolumes ;
          if (devFinite(PD) == 0)
          {
            DiagBreak() ;
          }
          if (devIsnan(PD))
          {
            DiagBreak() ;
          }
          if (std::isnan(PD))
          {
            DiagBreak() ;
          }
          if (mri_PD)
          {
            MRIsetVoxVal(mri_PD, x, y, z, 0, PD) ;
          }
        }
      }
    }
    MRIremoveNaNs(mri_PD, mri_PD) ;
    VectorFree(&v1) ;
    VectorFree(&v2) ;
    if (m_vox2vox)
      MatrixFree(&m_vox2vox) ;
  }

  return(mri_T2star) ;
}

static MRI *
compute_T2star_map(MRI **mri_flash, int nvolumes, int *scan_types,
                   MATRIX **Mreg, LTA *lta, MRI_BSPLINE **mri_flash_bsplines)
{
  MATRIX *mX, *mXpinv = NULL, *m_xform ;
  VECTOR *vY, *vParms = NULL, *v_src, *v_dst, *rasvec1, *rasvec2;

  int    x, y, z, e, width, height, depth, nscans, i ;
  MRI    *mri_T2star ;
  float  T2star, cond ;
  double val, xf, yf, zf ;

  if (lta)
  {
    m_xform = MatrixInverse(lta->xforms[0].m_L, NULL) ;
  }
  else
  {
    m_xform = NULL ;
  }

  v_src = VectorAlloc(4, MATRIX_REAL) ;
  v_dst = VectorAlloc(4, MATRIX_REAL) ;
  v_src->rptr[4][1] = 1.0 ;
  v_dst->rptr[4][1] = 1.0 ;

  rasvec1 =  MatrixCopy(v_src, NULL);
  rasvec2 =  MatrixCopy(v_src, NULL);


  for (i = nscans = 0 ; i < nvolumes ; i++)
  {
    if (scan_types[i] > nscans)
    {
      nscans = scan_types[i] ;
    }

    vox2ras[i] = MatrixCopy(mri_flash[i]->register_mat, NULL);
    ras2vox[i] = MatrixInverse(vox2ras[i], NULL);

    //    printf("volume %d belongs to volume %d\n", i, scan_types[i]-1);
    //    printf("registration matrix for it is:\n");
    //    MatrixPrint(stdout, Mreg[scan_types[i]-1]);
  }

  width = mri_flash[0]->width ;
  height = mri_flash[0]->height ;
  depth = mri_flash[0]->depth ;
  mri_T2star = MRIalloc(width, height, depth, MRI_FLOAT) ;
  if (!mri_T2star)
  {
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate T2* map", Progname) ;
  }
  MRIcopyHeader(mri_flash[0], mri_T2star) ;

  mX = MatrixAlloc(nvolumes, nscans+1, MATRIX_REAL) ;
  vY = VectorAlloc(nvolumes, MATRIX_REAL) ;
  vParms = VectorAlloc(nscans+1, MATRIX_REAL) ;
  for (e = 0 ; e < nvolumes ; e++)
  {
    *MATRIX_RELT(mX, e+1, 1) = -mri_flash[e]->te ;
    for (i = 1 ; i <= nscans ; i++)  /* which multi-echo set
                                        does this volume belong to */
    {
      if (scan_types[e] == i)
      {
        *MATRIX_RELT(mX, e+1, i+1) = 1 ;
      }
      else
      {
        *MATRIX_RELT(mX, e+1, i+1) = 0 ;
      }
    }
  }
  mXpinv = MatrixPseudoInverse(mX, mXpinv) ;
  if (!mXpinv)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "%s: could not invert matrix for T2* estimation",
                       Progname)) ;

  cond = MatrixConditionNumber(mX) ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        for (e = 0 ; e < nvolumes ; e++)
        {
          V3_X(v_src) = x ;
          V3_Y(v_src) = y ;
          V3_Z(v_src) = z ;
          if (m_xform)
          {
            /* map from output voxel-space to
               voxel-space of synth (T1, PD) */
            MatrixMultiply(m_xform, v_src, v_dst) ;
            MatrixCopy(v_dst, v_src);
          }

          if (Mreg[scan_types[e]-1])
          {
            MatrixMultiply(vox2ras[e],v_src,rasvec1);
            MatrixMultiply(Mreg[scan_types[e]-1],rasvec1,rasvec2);
            MatrixMultiply(ras2vox[e],rasvec2,v_dst);
          }
          else
          {
            MatrixCopy(v_src, v_dst) ;
          }

          xf = V3_X(v_dst) ;
          yf = V3_Y(v_dst) ;
          zf = V3_Z(v_dst) ;

          if (InterpMethod==SAMPLE_SINC)
            MRIsincSampleVolume(mri_flash[e], xf, yf, zf,
                                sinchalfwindow, &val) ;
          else if (InterpMethod==SAMPLE_CUBIC_BSPLINE)
      		{
          	MRIsampleBSpline(mri_flash_bsplines[e], xf, yf, zf, 0, &val);
      		}
					else
            MRIsampleVolumeType(mri_flash[e], xf, yf, zf,
                                &val, InterpMethod) ;
          if (val <= 0 || !std::isfinite(val))
          {
            val = 1E-6 ;
          }

          VECTOR_ELT(vY, e+1) = log(val) ;
        }
        vParms = MatrixMultiply(mXpinv, vY, vParms) ;
        if (*MATRIX_RELT(vParms, 1, 1) > 0)
        {
          T2star = 1 / *MATRIX_RELT(vParms, 1, 1) ;
        }
        else
        {
          T2star = 0 ;
        }

        if (T2star > 10000 || T2star < -1000)
        {
          DiagBreak() ;
        }
        if (!std::isfinite(T2star))
        {
          T2star = 0 ;
        }
        if (T2star > max_T2star)
        {
          T2star = max_T2star ;
        }
        MRIsetVoxVal(mri_T2star, x, y, z, 0, T2star) ;
      }
    }
  }

  MatrixFree(&mX) ;
  VectorFree(&vY) ;
  MatrixFree(&mXpinv) ;
  VectorFree(&vParms) ;
  if (m_xform)
  {
    MatrixFree(&m_xform) ;

  }
  VectorFree(&v_src) ;
  VectorFree(&v_dst) ;
  VectorFree(&rasvec1) ;
  VectorFree(&rasvec2) ;

  return(mri_T2star) ;
}


#if 0
static int
estimate_flip_angle_field(MRI *mri_T1,  MRI *mri_PD, MRI **mri_flash,
                          int nvolumes, int nfaf,
                          float *faf_coefs[3][2], LABEL *faf_label,
                          MRI *mri_faf)
{
  int    i,  k,  n, done, iter, width, height, depth, nsmall =  0, \
      xv, yv, zv ;
  double rms, dj_dalpha, last_rms, x, y, z, x0, y0, z0, w0x, w0y, w0z,\
  err,faf_scale,
      dalpha_dax, dalpha_day, dalpha_daz, dalpha_dbx, dalpha_dby, \
      dalpha_dbz, alpha,TR, exp_TR_T1,  numer, denom,
      dj_dax, dj_day, dj_daz, dj_dbx, dj_dby, dj_dbz, dt,  \
      pct_change/*, dj_dT1, dj_dPD, last_dT1, last_dPD, dT1, dPD */;
  double val, T1_wm, PD_wm, last_dax, last_day, last_daz, dax, \
  day, daz, dbx, dby, dbz, last_dbx, last_dby, last_dbz,
       Inorm, Snorm, T1, PD ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;
  x0 = width/2 ;
  y0 = height/2 ;
  z0 = depth/2 ;
  w0x = M_PI/x0 ;
  w0y = M_PI/y0 ;
  w0z = M_PI/z0 ;
  done = iter = 0 ;
  dt = base_dt / (double)(faf_label->n_points*nvolumes);
  T1_wm = PD_wm  = 0.0 ;
  for (k = 0 ; k  < faf_label->n_points ; k++)
  {
    x = faf_label->lv[k].x ;
    y = faf_label->lv[k].y ;
    z = faf_label->lv[k].z ;
    MRIsurfaceRASToVoxel(mri_T1, x, y, z, &x, &y, &z)  ;
    xv = nint(x) ;
    yv = nint(y) ;
    zv = nint(z) ;
    val = (double)MRIgetVoxVal(mri_T1, xv, yv, zv, 0) ;
    T1_wm += val ;
    val = (double)MRIgetVoxVal(mri_PD, xv, yv, zv, 0) ;
    PD_wm += val ;
  }

  T1_wm /= (double)faf_label->n_points ;
  PD_wm /= (double)faf_label->n_points ;
  printf("using white matter T1/PD = %2.1f/%2.1f\n", T1_wm, PD_wm) ;
  last_rms = compute_fa_rms(T1_wm, PD_wm, mri_flash,
                            nvolumes, faf_coefs,  nfaf,
                            faf_label, mri_T1, mri_PD)  ;
  printf("iter %03d: rms = %2.3f\n", 0, last_rms) ;

  /*last_dT1  = last_dPD = */
  last_dax = last_day = last_daz =
                          last_dbx = last_dby = last_dbz = 0 ;
  do
  {
    /*  compute current rms */
    last_rms = compute_fa_rms(T1_wm, PD_wm,
                              mri_flash, nvolumes, faf_coefs,
                              nfaf, faf_label, mri_T1,  mri_PD)  ;

    /* compute deltas for each coefficient*/
    for (n = 0  ; n <  nfaf ;  n++)
    {
      dalpha_dax = dalpha_day = dalpha_daz = dalpha_dbx =
          dalpha_dby = dalpha_dbz = dj_dalpha = 0 ;
      dj_dax = dj_day = dj_daz = dj_dbx = dj_dby = dj_dbz  = 0 ;

      /* compute dS/dalpha */
      for (k = 0 ; k  < faf_label->n_points ; k++)
      {
        x = faf_label->lv[k].x ;
        y = faf_label->lv[k].y ;
        z = faf_label->lv[k].z ;
        MRIsurfaceRASToVoxel(mri_T1, x, y, z, &x, &y, &z)  ;
        xv = nint(x) ;
        yv = nint(y) ;
        zv = nint(z) ;
#define USE_ALL_PARMS 0
#define PARAMETER_WT  0.0
#if USE_ALL_PARMS
        T1 = MRIgetVoxVal(mri_T1, xv,  yv, zv, 0) ;
        PD  = MRIgetVoxVal(mri_PD,  xv, yv, zv, 0) ;
        T1 = (1-PARAMETER_WT)*
             T1_wm+PARAMETER_WT*MRIgetVoxVal(mri_T1, xv,  yv, zv, 0) ;
        PD = (1-PARAMETER_WT)*
             PD_wm+PARAMETER_WT*MRIgetVoxVal(mri_PD,  xv, yv, zv, 0);
#else
        T1 = T1_wm ;
        PD = PD_wm  ;
#endif

        faf_scale = faf_coefs_to_scale(x-x0, y-y0,
                                       z-z0,  w0x, w0y, w0z,
                                       faf_coefs, nfaf)  ;
        for (Inorm = Snorm = 0.0, i = 0 ; i <  nvolumes  ;  i++)
        {
          alpha = faf_scale * mri_flash[i]->flip_angle ;
          TR = mri_flash[i]->tr ;
          val = MRIgetVoxVal(mri_flash[i], xv, yv, zv, 0) ;
          Inorm += val*val;
          val = FLASHforwardModel(alpha, TR, PD_wm, T1_wm) ;
          Snorm = val*val ;
        }
#define NORM_VALUES 0
#if NORM_VALUES
        Inorm = sqrt(Inorm) ;
        Snorm = sqrt(Snorm) ;
	check_finite(Inorm) ;
	check_finite(Snorm) ;
        if (FZERO(Inorm))
        {
          Inorm = 1  ;
        }
        if (FZERO(Snorm))
        {
          Snorm = 1 ;
        }
#else
        Snorm = Inorm = 1 ;
#endif
        for (i = 0 ; i <  nvolumes  ;  i++)
        {
          alpha = faf_scale * mri_flash[i]->flip_angle ;
          TR = mri_flash[i]->tr ;
          err = MRIgetVoxVal(mri_flash[i], xv, yv, zv, 0)/Inorm ;
          val = FLASHforwardModel(alpha, TR, PD, T1)/Snorm ;
          err  = (val - err) ;
          dalpha_dax = mri_flash[i]->flip_angle *
                       cos((n+1)*w0x*(x-x0));
          dalpha_day = mri_flash[i]->flip_angle *
                       cos((n+1)*w0y*(y-y0));
          dalpha_daz = mri_flash[i]->flip_angle *
                       cos((n+1)*w0z*(z-z0));

          dalpha_dbx = mri_flash[i]->flip_angle *
                       sin((n+1)*w0x*(x-x0));
          dalpha_dby = mri_flash[i]->flip_angle *
                       sin((n+1)*w0y*(y-y0));
          dalpha_dbz = mri_flash[i]->flip_angle *
                       sin((n+1)*w0z*(z-z0));

          exp_TR_T1 = exp(TR/T1_wm) ;
          numer = PD_wm  * (exp_TR_T1-1)*(exp_TR_T1*cos(alpha)-1) ;
          denom = (exp_TR_T1-cos(alpha)) ;
          denom *= denom;
          dj_dalpha = err*numer/denom ;

          dj_dax += dj_dalpha * dalpha_dax ;
          dj_day += dj_dalpha * dalpha_day ;
          dj_daz += dj_dalpha * dalpha_daz ;

          dj_dbx += dj_dalpha * dalpha_dbx ;
          dj_dby += dj_dalpha * dalpha_dby ;
          dj_dbz += dj_dalpha * dalpha_dbz ;
        }
      }
      dax = momentum*last_dax + dt * (-dj_dax) ;
      day = momentum*last_day + dt * (-dj_day) ;
      daz = momentum*last_daz + dt * (-dj_daz) ;
      dbx = momentum*last_dbx + dt * (-dj_dbx) ;
      dby = momentum*last_dby + dt * (-dj_dby) ;
      dbz = momentum*last_dbz + dt * (-dj_dbz) ;
      faf_coefs[0][0][n] += dax  ;
      faf_coefs[1][0][n] += day ;
      faf_coefs[2][0][n] += daz ;

      faf_coefs[0][1][n] += dbx ;
      faf_coefs[1][1][n] += dby  ;
      faf_coefs[2][1][n] += dbz  ;

      last_dax  = dax ;
      last_day  = day ;
      last_daz  = daz ;
      last_dbx  = dbx ;
      last_dby  = dby ;
      last_dbz  = dbz ;
      rms = compute_fa_rms(T1_wm, PD_wm, mri_flash, nvolumes,
                           faf_coefs,  nfaf, faf_label, mri_T1, mri_PD)  ;
    }


#if  0
    /* now update T1  and  PD estimates */
    dj_dT1 = dj_dPD = 0 ;
    for (k = 0 ; k  < faf_label->n_points ; k++)
    {
      x = faf_label->lv[k].x ;
      y = faf_label->lv[k].y ;
      z = faf_label->lv[k].z ;
      MRIsurfaceRASToVoxel(mri_T1, x, y, z, &x, &y, &z)  ;
      xv = nint(x) ;
      yv = nint(y) ;
      zv = nint(z) ;
#if USE_ALL_PARMS
      T1 = (1-PARAMETER_WT)*T1_wm+PARAMETER_WT*
           MRIgetVoxVal(mri_T1, xv,  yv, zv, 0) ;
      PD = (1-PARAMETER_WT)*PD_wm+PARAMETER_WT*
           MRIgetVoxVal(mri_PD,  xv, yv, zv, 0);
#else
      T1 = T1_wm;
      PD  = PD_wm ;
#endif

      faf_scale = faf_coefs_to_scale(x-x0, y-y0,  z-z0,
                                     w0x, w0y, w0z, faf_coefs, nfaf)  ;
      for (i = 0 ; i <  nvolumes  ;  i++)
      {
        alpha = faf_scale * mri_flash[i]->flip_angle ;
        TR = mri_flash[i]->tr ;
        err = MRIgetVoxVal(mri_flash[i], xv, yv, zv, 0) ;
        val = FLASHforwardModel(alpha, TR, PD, T1) ;
        err  = (val - err) ;

        exp_TR_T1 = exp(TR/T1) ;
        numer = PD  * exp_TR_T1*TR*(cos(alpha)-1)*sin(alpha) ;
        denom = T1*(exp_TR_T1-cos(alpha)) ;
        denom *= denom;
        dj_dT1 += err*numer/denom ;
        dj_dPD += err*val / PD ;
      }
    }


#if 0
    dT1 = momentum*last_dT1 + 1e6*dt * (-dj_dT1) ;
    dPD = momentum*last_dPD + 1e7*dt * (-dj_dPD) ;
    T1_wm += dT1 ;
#if NORM_VALUES
    PD_wm += dPD ;
#endif
#endif
    last_dT1  =  dT1 ;
    last_dPD = dPD ;
#endif

    rms = compute_fa_rms(T1_wm, PD_wm, mri_flash, nvolumes,
                         faf_coefs,  nfaf, faf_label, mri_T1,  mri_PD)  ;
    pct_change  = 100*(last_rms-rms)/last_rms;
    printf("iter %03d: rms = %2.3f (%2.5f%%), T1,PD=(%2.1f,%2.1f)\n",
           iter+1, rms, pct_change, T1_wm, PD_wm) ;
    printf("xb(r) = ") ;
    for (n = 0 ;  n <  nfaf ; n++)
    {
      printf("(%2.4f,%2.4f) ", faf_coefs[0][0][n], faf_coefs[0][1][n]) ;
    }
    printf("\nyb(r) = ") ;
    for (n = 0 ;  n <  nfaf ; n++)
    {
      printf("(%2.4f,%2.4f) ", faf_coefs[1][0][n], faf_coefs[1][1][n]) ;
    }
    printf("\nzb(r) = ") ;
    for (n = 0 ;  n <  nfaf ; n++)
    {
      printf("(%2.4f,%2.4f) ", faf_coefs[2][0][n], faf_coefs[2][1][n]) ;
    }
    printf("\n") ;
    if (pct_change < 0.000001)
    {
      nsmall++ ;
    }
    else
    {
      nsmall = 0 ;
    }
    if (pct_change <  0)  /*  reset all momentum */
    {
      /*last_dT1 = last_dPD = */
      last_dax = last_day = last_daz =
                              last_dbx = last_dby = last_dbz = 0 ;
    }

    done  = (iter++ > 500) || (nsmall > 5)  ;
  }
  while (!done) ;

  for  (x = 0  ; x <  width ;  x++)
  {
    for (y = 0 ; y <  height  ;  y++)
    {
      for (z = 0 ; z <  depth  ;  z++)
      {
        faf_scale = faf_coefs_to_scale(x-x0, y-y0, z-z0,
                                       w0x,  w0y, w0z,
                                       faf_coefs, nfaf) ;
        MRIsetVoxVal(mri_faf, x, y, z, 0, faf_scale);
      }
    }
  }
  return(NO_ERROR)  ;
}

static double
compute_fa_rms(double T1_wm, double PD_wm, MRI **mri_flash, int nvolumes,
               float *faf_coefs[3][2],  int nfaf, LABEL *faf_label,
               MRI *mri_T1,  MRI *mri_PD)
{
  int   k, i, xv, yv, zv ;
  double rms, x, y, z, x0, y0, z0, w0x, w0y, w0z, val, sval,\
  faf_scale, Inorm, Snorm, T1, PD ;

  x0 = mri_flash[0]->width/2 ;
  y0 = mri_flash[0]->height/2 ;
  z0 = mri_flash[0]->depth/2 ;
  w0x = M_PI/x0 ;
  w0y = M_PI/y0 ;
  w0z = M_PI/z0 ;
  for (rms = 0.0, k = 0 ; k  < faf_label->n_points ; k++)
  {
    x = faf_label->lv[k].x ;
    y = faf_label->lv[k].y ;
    z = faf_label->lv[k].z ;
    MRIsurfaceRASToVoxel(mri_flash[0], x, y, z, &x, &y, &z)  ;
    xv = nint(x) ;
    yv = nint(y) ;
    zv = nint(z) ;
#if  USE_ALL_PARMS
    T1 = (1-PARAMETER_WT)*T1_wm+PARAMETER_WT*
         MRIgetVoxVal(mri_T1, xv,  yv, zv, 0) ;
    PD = (1-PARAMETER_WT)*PD_wm+PARAMETER_WT*
         MRIgetVoxVal(mri_PD,  xv, yv, zv, 0);
#else
    T1 = T1_wm  ;
    PD = PD_wm ;
#endif
    faf_scale = faf_coefs_to_scale(x-x0, y-y0, z-z0,
                                   w0x,  w0y, w0z,
                                   faf_coefs, nfaf) ;
    for (Inorm = Snorm = 0.0, i = 0 ; i < nvolumes ; i++)
    {
      sval = FLASHforwardModel(faf_scale*mri_flash[i]->flip_angle,
                               mri_flash[i]->tr, PD_wm, T1_wm) ;
      val = MRIgetVoxVal(mri_flash[i], xv,  yv,  zv, 0) ;
      Snorm  += sval*sval;
      Inorm += val*val ;
    }
#if NORM_VALUES
    Snorm = sqrt(Snorm) ;
    Inorm = sqrt(Inorm)  ;
    check_finite(Inorm) ;
    check_finite(Snorm) ;
    if (FZERO(Snorm))
    {
      Snorm =  1 ;
    }
    if (FZERO(Inorm))
    {
      Inorm = 1 ;
    }
#else
    Snorm =  Inorm = 1  ;
#endif
    for (i = 0 ; i < nvolumes ; i++)
    {
      sval = FLASHforwardModel(faf_scale*mri_flash[i]->flip_angle,
                               mri_flash[i]->tr, 1, T1_wm)/Snorm ;
      val = MRIgetVoxVal(mri_flash[i], xv,  yv,  zv, 0)/Inorm ;
      rms += (sval-val)*(sval-val);
    }
  }
  check_finite(sqrt(rms/(double)(faf_label->n_points*nvolumes))) ;

  return(sqrt(rms/(double)(faf_label->n_points*nvolumes))) ;
}

static double
faf_coefs_to_scale(double x, double y, double z,
                   double w0x, double w0y, double w0z,
                   float *faf_coefs[3][2], int nfaf)
{
  int   n ;
  double xb, yb, zb ;

  for (xb = 1.0, n=1 ; n <= nfaf ; n++)
    xb += faf_coefs[0][0][n-1] * cos(w0x*n*(x)) +
          faf_coefs[0][1][n-1] * sin(w0x*n*(x)) ;
  for (yb = 1.0, n=1 ; n <= nfaf ; n++)
    yb += faf_coefs[1][0][n-1] * cos(w0y*n*(y)) +
          faf_coefs[1][1][n-1] * sin(w0y*n*(y)) ;
  for (zb = 1.0, n=1 ; n <= nfaf ; n++)
    zb += faf_coefs[2][0][n-1] * cos(w0z*n*(z)) +
          faf_coefs[2][1][n-1] * sin(w0z*n*(z)) ;
  return(xb*yb*zb) ;

}
#endif

static int findUniqueTETRFA(MRI *mri[], int numvolumes,
                            float *ptr, float *pte, double *pfa)
{
  float TR, TE;
  double FA;
  int i;
  int flag = 0;
  // get the first volume values
  TR = mri[0]->tr;
  TE = mri[0]->te;
  FA = mri[0]->flip_angle;

  // note that flag = 4 will always overwrite flag =1,2
  //if the last volume has different FA than the first one;
  //anyway, it's not important -xh
  for (i = 1; i < numvolumes; ++i)
  {
    if (!FZERO(TR - mri[i]->tr))
    {
      fprintf(stderr, "non-equal TR found for the volume %d.\n", i);
      flag = 1;
    }
    if (!FZERO(TE - mri[i]->te))
    {
      fprintf(stderr, "non-equal TE found for the volume %d.\n", i);
      flag = 2;
    }
    if (!FZERO(FA - mri[i]->flip_angle))
    {
      fprintf(stderr, "non-equal flip_angle found for the "
              "volume %d.\n", i);
      flag = 4;
    }
  }

  if ((flag & 1) == 1)
  {
    fprintf(stderr, "TR is set to zero.\n");
    TR = 0.;
  }
  if ((flag & 2) == 2)
  {
    fprintf(stderr, "TE is set to zero.\n");
    TE = 0.;
  }
  if ((flag & 4) == 4)
  {
    fprintf(stderr, "Flip_angle is set to zero.\n");
    FA = 0.;
  }
  *ptr = TR;
  *pte = TE;
  *pfa = FA;

  return flag;
}

static int resetTRTEFA(MRI *mri, float tr, float te, double fa)
{
  mri->tr = tr;
  mri->te = te;
  mri->flip_angle = fa;
  return NO_ERROR;
}
