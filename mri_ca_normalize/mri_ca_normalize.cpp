/**
 * @brief Normalize a volume making use of subcortical atlas data
 *
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "tags.h"
#include "cma.h"
#include "mrinorm.h"
#include "version.h"
#include "mri2.h"
#include "fsinit.h"

#define MM_FROM_EXTERIOR  5  // distance into brain mask to go when erasing super bright CSF voxels

#define MIN_WM_BIAS_PCT 0.8
#define MAX_WM_BIAS_PCT 1.2

const char *Progname ;

static double extra_norm_range = 0.0 ;

static int fill_in_sample_means(GCA_SAMPLE *gcas, GCA *gca, int nsamples);
MRI *normalizeChannelFromLabel
(MRI *mri_in, MRI *mri_dst, MRI *mri_seg, double *fas, int input_index);
MRI *normalizeFromLabel
(MRI *mri_in, MRI *mri_dst, MRI *mri_seg, double *fas) ;
static MRI *normalize_from_segmentation_volume
(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, int *structs, int nstructs) ;

static int dilate_mask = 0 ;

static double TRs[MAX_GCA_INPUTS] ;
static double fas[MAX_GCA_INPUTS] ;
static double TEs[MAX_GCA_INPUTS] ;

static int noedit = 0 ;
static int remove_cerebellum = 0 ;
static int remove_lh = 0 ;
static int remove_rh = 0 ;

static int file_only = 0 ;
static char *normalized_transformed_sample_fname = NULL ;
static char *T2_mask_fname = NULL ;
static double T2_thresh = 0 ;
static char *aparc_aseg_fname = NULL ;
static char *mask_fname = NULL ;
static char *sample_fname = NULL ;
static char *ctl_point_fname = NULL ;
static int novar = 0 ;

static double bias_sigma = 4.0 ;
static float min_prior = 0.6 ;
static FILE *diag_fp = NULL ;

static void usage_exit(int code) ;
static int get_option(int argc, char *argv[]) ;
static int copy_ctrl_points_to_volume(GCA_SAMPLE *gcas, int nsamples,
                                      MRI *mri_ctrl, int frame) ;
static GCA_SAMPLE *copy_ctrl_points_from_volume(GCA *gca,TRANSFORM *transform,
    int *pnsamples, MRI *mri_ctrl,
    int frame) ;

static MRI  *mri_aseg = NULL ;
static float aseg_thresh = 0 ;

static char *seg_fname = NULL ;
static char *long_seg_fname = NULL ;
static char *renormalization_fname = NULL ;
static double TR = 0.0, TE = 0.0, alpha = 0.0 ;
static char *tissue_parms_fname = NULL ;
static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;
static double min_region_prior
(GCA *gca, int xp, int yp, int zp, int wsize, int label) ;
static GCA_SAMPLE *find_control_points
(GCA *gca, GCA_SAMPLE *gcas, int total_nsamples,
 int *pnorm_samples, int nregions, int label,
 MRI *mri_in, TRANSFORM *transform, double min_prior,
 double ctrl_point_pct, char *fsample_fname) ;

static GCA_SAMPLE *gcas_concatenate
(GCA_SAMPLE *gcas1, GCA_SAMPLE *gcas2, int n1, int n2);
static int  gcas_bounding_box
(GCA_SAMPLE *gcas, int nsamples, int *pxmin, int *pymin, int *pzmin,
 int *pxmax, int *pymax, int *pzmax, int label) ;
static int  uniform_region
(GCA *gca, MRI *mri, TRANSFORM *transform,
 int x, int y, int z, int wsize, GCA_SAMPLE *gcas, float nsigma) ;
static int  discard_unlikely_control_points
(GCA *gca, GCA_SAMPLE *gcas_struct, int struct_samples,
 MRI *mri_in, TRANSFORM *transform, const char *name) ;

/*
  command line consists of three inputs:

  argv[1]  - input volume
  argv[2]  - atlas (gca)
  argv[3]  - transform (lta/xfm/m3d)
  argv[4]  - output volume
*/

#define DEFAULT_CTL_POINT_PCT   .25
static double ctl_point_pct = DEFAULT_CTL_POINT_PCT ;

static int normalization_structures[] =
{
  Left_Cerebral_White_Matter,
  Right_Cerebral_White_Matter,
  Left_Cerebellum_White_Matter,
  Right_Cerebellum_White_Matter,
  Brain_Stem
} ;

#define NSTRUCTURES (sizeof(normalization_structures) / sizeof(normalization_structures[0]))

static int nregions = 3 ;  /* divide each struct into 3x3x3 regions */

static char *ctrl_point_fname = NULL ;
static char *read_ctrl_point_fname = NULL ;

int
main(int argc, char *argv[])
{
  char         *gca_fname, *in_fname, *out_fname, **av, *xform_fname ;
  MRI          *mri_in = NULL, *mri_norm = NULL, *mri_tmp, *mri_ctrl = NULL ;
  GCA          *gca ;
  int          ac, nargs, nsamples, msec, minutes, seconds;
  int          i, struct_samples, norm_samples = 0, n, input, ninputs ;
  Timer start ;
  GCA_SAMPLE   *gcas, *gcas_norm = NULL, *gcas_struct ;
  TRANSFORM    *transform = NULL ;

  FSinit();
  
  std::string cmdline = getAllInfo(argc, argv, "mri_ca_normalize");

  nargs = handleVersionOption(argc, argv, "mri_ca_normalize");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  setRandomSeed(-1L) ;
  Progname = argv[0] ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    ErrorExit
    (ERROR_BADPARM,
     "usage: %s [<options>] <inbrain1> <inbrain2> ... "
     "<atlas> <transform file> <output1> <output2> ...\n",
     Progname) ;

  ninputs = (argc - 2) / 2 ;
  printf("reading %d input volume%s\n", ninputs, ninputs > 1 ? "s" : "") ;
  in_fname = argv[1] ;
  gca_fname = argv[1+ninputs] ;
  xform_fname = argv[2+ninputs] ;
  out_fname = argv[3+ninputs] ;

  if (read_ctrl_point_fname)
  {
    mri_ctrl = MRIread(read_ctrl_point_fname) ;
    if (mri_ctrl == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read precomputed control points from %s",
                Progname, read_ctrl_point_fname) ;
  }
  start.reset() ;
  printf("reading atlas from '%s'...\n", gca_fname) ;
  fflush(stdout) ;

  if (long_seg_fname) // read in a segmentation and turn it into a ctrl point volume
  {
    MRI *mri_seg, *mri_tmp2 = NULL ;
    int i ;

    mri_seg = MRIread(long_seg_fname) ;
    if (mri_seg == NULL)
    {
      ErrorExit(ERROR_NOFILE, "%s: could not read segmentation volume %s", Progname,long_seg_fname) ;
    }
    mri_tmp = MRIclone(mri_seg, NULL) ;
    mri_ctrl = MRIclone(mri_tmp, NULL) ;
    for (i = 0 ; i < NSTRUCTURES ; i++)
    {
      MRIcopyLabel(mri_seg, mri_tmp, normalization_structures[i]) ;
      MRIbinarize(mri_tmp, mri_tmp, 1, 0, 1) ;
      MRIerode(mri_tmp, mri_tmp) ;
      MRIerode(mri_tmp, mri_tmp) ;
      mri_tmp2 = MRImask(mri_seg, mri_tmp, mri_tmp2, 0, 0) ;

      MRIadd(mri_tmp2 , mri_ctrl, mri_ctrl) ;
      MRIclear(mri_tmp) ;
      MRIclear(mri_tmp2) ;
    }
    MRIfree(&mri_tmp) ;
    MRIfree(&mri_seg) ;
  }

  if (seg_fname == NULL)
  {
    gca = GCAread(gca_fname) ;
    if (gca == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",
                Progname, gca_fname) ;
    printf("reading transform from '%s'...\n", xform_fname) ;
    if (gca->ninputs != ninputs)
      ErrorExit(ERROR_BADPARM, "gca has %d frames, but %d inputs specified\n",
             gca->ninputs,ninputs) ;
    fflush(stdout) ;
    transform = TransformRead(xform_fname) ;
    if (!transform)
      ErrorExit
      (ERROR_BADPARM,
       "%s: could not open xform file %s", Progname,xform_fname) ;

    if (novar)
    {
      GCAunifyVariance(gca) ;
    }

    if (remove_lh)
      GCAremoveHemi(gca, 1) ; // for exvivo contrast
    if (remove_rh)
      GCAremoveHemi(gca, 0) ; // for exvivo contrast
    if (remove_cerebellum)
    {
      GCAremoveLabel(gca, Brain_Stem) ;
      GCAremoveLabel(gca, Left_Cerebellum_Cortex) ;
      GCAremoveLabel(gca, Left_Cerebellum_White_Matter) ;
      GCAremoveLabel(gca, Right_Cerebellum_White_Matter) ;
      GCAremoveLabel(gca, Right_Cerebellum_Cortex) ;
      gca->flags |= GCA_NO_CEREBELLUM ;
    }
  
    if (renormalization_fname)
    {
      FILE   *fp ;
      int    *labels, nlines, i ;
      float  *intensities, f1, f2 ;
      char   *cp, line[STRLEN] ;

      fp = fopen(renormalization_fname, "r") ;
      if (!fp)
        ErrorExit(ERROR_NOFILE, "%s: could not read %s",
                  Progname, renormalization_fname) ;

      cp = fgetl(line, 199, fp) ;
      nlines = 0 ;
      while (cp)
      {
        nlines++ ;
        cp = fgetl(line, 199, fp) ;
      }
      rewind(fp) ;
      printf("reading %d labels from %s...\n",
             nlines,renormalization_fname) ;
      labels = (int *)calloc(nlines, sizeof(int)) ;
      intensities = (float *)calloc(nlines, sizeof(float)) ;
      cp = fgetl(line, 199, fp) ;
      for (i = 0 ; i < nlines ; i++)
      {
        sscanf(cp, "%e  %e", &f1, &f2) ;
        labels[i] = (int)f1 ;
        intensities[i] = f2 ;
        if (labels[i] == Left_Cerebral_White_Matter)
        {
          DiagBreak() ;
        }
        cp = fgetl(line, 199, fp) ;
      }
      GCArenormalizeIntensities(gca, labels, intensities, nlines) ;
      free(labels) ;
      free(intensities) ;
    }
  }
  else
  {
    gca = NULL ;  /* don't need atlas if using segmentation */
  }

  if (long_seg_fname)   // longitudinal analysis - regularize means
  {
    GCAregularizeConditionalDensities(gca, .5) ;
  }

  for (input = 0 ; input < ninputs ; input++)
  {
    in_fname = argv[1+input] ;
    printf("reading input volume from %s...\n", in_fname) ;
    mri_tmp = MRIread(in_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s",
                Progname, in_fname) ;
    MRImakePositive(mri_tmp, mri_tmp) ;
    if (mri_tmp->type != MRI_UCHAR && mri_in && mri_in->type == MRI_UCHAR)  // scale it down to fit
    {
      MRI *mri_changed = MRIchangeType(mri_tmp, MRI_UCHAR, 0, 255, 0) ;
      MRIfree(&mri_tmp) ; mri_tmp = mri_changed ;
    }
      
    if (mri_tmp && ctrl_point_fname && !mri_ctrl)
    {
      mri_ctrl = MRIallocSequence(mri_tmp->width, mri_tmp->height,
                                  mri_tmp->depth,MRI_FLOAT, nregions*2) ; // labels and means
      MRIcopyHeader(mri_tmp, mri_ctrl) ;
    }
    if (alpha > 0)
    {
      mri_tmp->flip_angle = alpha ;
    }
    if (TR > 0)
    {
      mri_tmp->tr = TR ;
    }
    if (TE > 0)
    {
      mri_tmp->te = TE ;
    }

    TRs[input] = mri_tmp->tr ;
    fas[input] = mri_tmp->flip_angle ;
    TEs[input] = mri_tmp->te ;

    if (input == 0)
    {
      mri_in =
        MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth,
                         mri_tmp->type, ninputs) ;
      if (!mri_in)
        ErrorExit(ERROR_NOMEMORY,
                  "%s: could not allocate input volume %dx%dx%dx%d",
                  mri_tmp->width,mri_tmp->height,mri_tmp->depth,ninputs) ;
      MRIcopyHeader(mri_tmp, mri_in) ;
    }

    if (mask_fname)
    {
      MRI *mri_mask ;

      mri_mask = MRIread(mask_fname) ;
      if (!mri_mask)
        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                  Progname, mask_fname) ;

      if (noedit == 0)
	MRIreplaceValues(mri_mask, mri_mask, WM_EDITED_OFF_VAL, 0) ;
      MRIclose(mri_mask, mri_mask) ;
      if (dilate_mask)
      {
	int i ;
	MRI *mri_mask_tmp = MRIclone(mri_mask, NULL) ;
	printf("dilating mask %d times\n", dilate_mask) ;
	for (i = 0 ; i < dilate_mask ; i++)
	{
	  MRIdilate(mri_mask, mri_mask_tmp) ;
	  MRIcopy(mri_mask_tmp, mri_mask) ;
	}
	MRIfree(&mri_mask_tmp) ;
      }


      MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
      MRIfree(&mri_mask) ;
    }
    if (T2_mask_fname)
    {
      MRI *mri_T2, *mri_aparc_aseg = nullptr;

      mri_T2 = MRIread(T2_mask_fname) ;
      if (!mri_T2)
	ErrorExit(ERROR_NOFILE, "%s: could not open T2 mask volume %s.\n",
		  Progname, mask_fname) ;
      if (aparc_aseg_fname)   // use T2 and aparc+aseg to remove non-brain stuff
      {
 	mri_aparc_aseg = MRIread(aparc_aseg_fname) ;
	if (mri_aparc_aseg == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not open aparc+aseg volume %s.\n",
		    Progname, aparc_aseg_fname) ;
      }

      MRImask_with_T2_and_aparc_aseg(mri_tmp, mri_tmp, mri_T2, mri_aparc_aseg, T2_thresh, MM_FROM_EXTERIOR) ;
      MRIfree(&mri_T2) ; MRIfree(&mri_aparc_aseg) ;
    }
    MRIcopyFrame(mri_tmp, mri_in, 0, input) ;
  }
  MRIaddCommandLine(mri_in, cmdline) ;

  if (seg_fname == NULL)
  {
    if (gca->type == GCA_PARAM)
    {
      GCA *gca_tmp ;

      printf("mapping T1/PD atlas into %d-dimensional FLASH space atlas\n",
             mri_in->nframes) ;
      // that means gca->ninputs = nframes
      gca_tmp = GCAcreateFlashGCAfromParameterGCA
                (gca, TRs, fas, TEs, mri_in->nframes, GCA_DEFAULT_NOISE_PARAMETER);
      // now the type is set gca->type = GCA_FLASH
      GCAfree(&gca) ;
      gca = gca_tmp ;
      GCAhistoScaleImageIntensities(gca, mri_in, 1) ;
    }
    else if (gca->type == GCA_FLASH)
    {
      GCA *gca_tmp ;

      int need_map_flag = 0;
      int n;

      if(gca->ninputs != ninputs)
      {
        need_map_flag = 1;
      }
      else
      {
        for (n = 0 ; n < mri_in->nframes; n++)
        {
          if(!FZERO(gca->TRs[n] - TRs[n]))
          {
            need_map_flag = 1;
          }
          if(!FZERO(gca->FAs[n] - fas[n]))
          {
            need_map_flag = 1;
          }
          if(!FZERO(gca->TEs[n] - TEs[n]))
          {
            need_map_flag = 1;
          }
        }
      }

      if(need_map_flag)
      {
        printf("mapping %d-dimensional flash atlas "
               "into %d-dimensional input space\n",
               gca->ninputs, ninputs) ;

        gca_tmp = GCAcreateFlashGCAfromFlashGCA
                  (gca, TRs, fas, TEs, mri_in->nframes) ;
        GCAfree(&gca) ;
        gca = gca_tmp ;
      }

      if (novar)
      {
        GCAunifyVariance(gca) ;
      }

      GCAhistoScaleImageIntensities(gca, mri_in, 1) ;
    }
    else
    {
      GCAhistoScaleImageIntensities(gca, mri_in, 1) ;
    }

    if (example_T1)
    {
      MRI *mri_T1, *mri_seg ;

      mri_seg = MRIread(example_segmentation) ;
      if (!mri_seg)
        ErrorExit
        (ERROR_NOFILE,
         "%s: could not read example segmentation from %s",
         Progname, example_segmentation) ;
      mri_T1 = MRIread(example_T1) ;
      if (!mri_T1)
        ErrorExit(ERROR_NOFILE,"%s: could not read example T1 from %s",
                  Progname, example_T1) ;
      printf("scaling atlas intensities using specified examples...\n") ;
      MRIeraseBorderPlanes(mri_seg, 1) ;
      GCArenormalizeToExample(gca, mri_seg, mri_T1) ;
      MRIfree(&mri_seg) ;
      MRIfree(&mri_T1) ;
    }

    if (tissue_parms_fname)   /* use FLASH forward model */
    {
      GCArenormalizeToFlash(gca, tissue_parms_fname, mri_in) ;
    }
  }

  if (seg_fname || mri_aseg)   /* use segmentation volume to drive normalization */
  {
    MRI *mri_seg ;
    int  structs[MAX_CMA_LABELS], nstructs ;
    if (mri_aseg)
      mri_seg = mri_aseg ;
    else
      mri_seg = MRIread(seg_fname) ;
    if (!mri_seg)
      ErrorExit
      (ERROR_NOFILE,
       "%s: could not read segmentation volume %s...\n",
       Progname, seg_fname);


    nstructs = 0 ;
    structs[nstructs++] = Left_Cerebral_White_Matter ;
    structs[nstructs++] = Right_Cerebral_White_Matter ;
    mri_norm = normalize_from_segmentation_volume
               (mri_in, NULL, mri_seg, structs, nstructs) ;
    MRIcopy(mri_norm, mri_in) ;  /* for next pass through */
    MRIfree(&mri_norm) ;
    MRIfree(&mri_seg) ;
  }

  if (seg_fname == NULL)  // only run this if using an aseg in the previous block, not a manual one (i.e. -aseg)
  {
    int j ;

    gcas = GCAfindAllSamples(gca, &nsamples, NULL, 1) ;
    printf("using %d sample points...\n", nsamples) ;
    GCAcomputeSampleCoords(gca, mri_in, gcas, nsamples, transform) ;
    if (sample_fname)
      GCAtransformAndWriteSamples
      (gca, mri_in, gcas, nsamples, sample_fname, transform) ;

    if (Gx >= 0)
    {
      int xp, yp, zp, xn, yn, zn ;

      GCAsourceVoxelToPrior(gca, mri_ctrl, transform, Gx, Gy, Gz, &xp, &yp, &zp);
      GCAsourceVoxelToNode(gca, mri_ctrl, transform, Gx, Gy, Gz, &xn, &yn, &zn);
      printf("source voxel (%d, %d, %d) maps to prior (%d, %d, %d)\n", Gx, Gy, Gz, xn, yn, zn) ;
    }

    for (j = 0 ; j < 1 ; j++)
    {
      for (n = 1 ; n <= nregions ; n++)
      {
        if (long_seg_fname)
        {
          gcas_norm =
            copy_ctrl_points_from_volume(gca, transform, &norm_samples,
                                         mri_ctrl, 0) ;
          fill_in_sample_means(gcas_norm, gca,norm_samples);
        }
        else if (read_ctrl_point_fname)
        {
          gcas_norm =
            copy_ctrl_points_from_volume(gca, transform, &norm_samples,
                                         mri_ctrl, n-1) ;
        }
        else
        {
          for (norm_samples = i = 0 ; i < NSTRUCTURES ; i++)
          {
            if (normalization_structures[i] == Gdiag_no)
            {
              DiagBreak() ;
            }
            printf("finding control points in %s....\n",
                   cma_label_to_name(normalization_structures[i])) ;
            gcas_struct = find_control_points
                          (gca, gcas, nsamples, &struct_samples, n,
                           normalization_structures[i], mri_in, transform, min_prior,
                           ctl_point_pct, sample_fname) ;
            discard_unlikely_control_points
            (gca, gcas_struct, struct_samples, mri_in, transform,
             cma_label_to_name(normalization_structures[i])) ;
            if (mri_ctrl && ctrl_point_fname) // store the samples
            {
              copy_ctrl_points_to_volume(gcas_struct, struct_samples, mri_ctrl, n-1) ;
            }
            if (i)
            {
              GCA_SAMPLE *gcas_tmp ;
              gcas_tmp = gcas_concatenate
                         (gcas_norm, gcas_struct, norm_samples, struct_samples) ;
              free(gcas_norm) ;
              norm_samples += struct_samples ;
              gcas_norm = gcas_tmp ;
            }
            else
            {
              gcas_norm = gcas_struct ;
              norm_samples = struct_samples ;
            }
          }

        }
        if (norm_samples == 0)
        {
          printf("skipping region %d with no control points detected\n", n) ;
          continue ;
        }
        printf("using %d total control points "
               "for intensity normalization...\n", norm_samples) ;
        if (normalized_transformed_sample_fname)
          GCAtransformAndWriteSamples(gca, mri_in, gcas_norm, norm_samples,
                                      normalized_transformed_sample_fname,
                                      transform) ;
        mri_norm = GCAnormalizeSamplesAllChannels
                   (mri_in, gca, gcas_norm, file_only ? 0 :norm_samples,
                    transform, ctl_point_fname, bias_sigma) ;
        if (Gdiag & DIAG_WRITE)
        {
          char fname[STRLEN] ;
          sprintf(fname, "norm%d.mgz", n) ;
          printf("writing normalized volume to %s...\n", fname) ;
          MRIwrite(mri_norm, fname) ;
          sprintf(fname, "norm_samples%d.mgz", n) ;
          GCAtransformAndWriteSamples(gca, mri_in, gcas_norm, norm_samples,
                                      fname, transform) ;

        }
        MRIcopy(mri_norm, mri_in) ;  /* for next pass through */
        MRIfree(&mri_norm) ;
      }
    }
  }

  for (input = 0 ; input < ninputs ; input++)
  {
    out_fname  = argv[3+ninputs+input] ;
    printf("writing normalized volume to %s...\n", out_fname) ;
    mri_in->tr = TRs[input] ;
    mri_in->flip_angle = fas[input] ;
    mri_in->te = TEs[input] ;
    if (MRIwriteFrame(mri_in, out_fname, input)  != NO_ERROR)
      ErrorExit(ERROR_BADFILE, "%s: could not write normalized volume to %s",
                Progname, out_fname);
  }

  if (ctrl_point_fname)
  {
    printf("writing control points to %s\n", ctrl_point_fname) ;
    MRIwrite(mri_ctrl, ctrl_point_fname) ;
    MRIfree(&mri_ctrl) ;
  }
  MRIfree(&mri_in) ;

#if 1
  printf("freeing GCA...") ;
  if (gca)
  {
    GCAfree(&gca) ;
  }
#endif
  printf("done.\n") ;
  if (mri_in)
  {
    MRIfree(&mri_in) ;
  }
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("normalization took %d minutes and %d seconds.\n",
         minutes, seconds) ;
  if (diag_fp)
  {
    fclose(diag_fp) ;
  }
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
  StrUpper(option) ;
  if (!strcmp(option, "FSAMPLES"))
  {
    sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing control points to %s...\n", sample_fname) ;
  }
  else if (!strcmp(option, "-HELP")||!strcmp(option, "-USAGE"))
  {
    usage_exit(0) ;
  }
  else if (!stricmp(option, "T2MASK"))
  {
    T2_mask_fname = argv[2] ;
    T2_thresh = atof(argv[3]) ;
    nargs = 2 ;
    printf("using T2 volume %s thresholded at %f to mask input volume...\n", 
	   T2_mask_fname, T2_thresh) ;
  }
  else if (!stricmp(option, "extra_norm"))
  {
    extra_norm_range = atof(argv[2]) ;
    nargs = 1 ;
    printf("expanding norm range to [%2.1f --> %2.1f]\n",  (1.0-extra_norm_range) * MIN_WM_BIAS_PCT * DEFAULT_DESIRED_WHITE_MATTER_VALUE,
	   (1.0+extra_norm_range) * MAX_WM_BIAS_PCT* DEFAULT_DESIRED_WHITE_MATTER_VALUE) ;

  }
  else if (!stricmp(option, "AMASK"))
  {
    aparc_aseg_fname = argv[2] ;
    T2_mask_fname = argv[3] ;
    T2_thresh = atof(argv[4]) ;
    nargs = 3 ;
    printf("using aparc+aseg vol %s and T2 volume %s thresholded at %f to mask input volume...\n", 
	   aparc_aseg_fname, T2_mask_fname, T2_thresh) ;
  }
  else if (!strcmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!strcmp(option, "ASEG"))
  {
    mri_aseg = MRIread(argv[2]) ;
    if (mri_aseg == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s", Progname, argv[2]) ;
    aseg_thresh = atof(argv[3]) ;
    nargs = 2 ;
    printf("using top %2.2f%% of white matter aseg volume %s to do first pass normalization...\n", 
	   aseg_thresh*100, argv[2]) ;
  }
  else if (!strcmp(option, "SEG"))
  {
    seg_fname = argv[2] ;
    nargs = 1 ;
    printf("using segmentation volume %s to generate control points...\n",
           seg_fname) ;
  }
  else if (!strcmp(option, "LONG"))
  {
    long_seg_fname = argv[2] ;
    nargs = 1 ;
    printf("using longitudinal segmentation volume %s to generate control points...\n",
           long_seg_fname) ;
  }
  else if (!strcmp(option, "NOEDIT"))
  {
    noedit = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%sremoving edited off voxels in the mask\n", noedit ? "not " : "") ;
  }
  else if (!stricmp(option, "LH"))
  {
    remove_rh = 1  ;
    printf("removing right hemisphere labels\n") ;
  }
  else if (!stricmp(option, "RH"))
  {
    remove_lh = 1  ;
    printf("removing left hemisphere labels\n") ;
  }
  else if (!stricmp(option, "NOCEREBELLUM"))
  {
    remove_cerebellum = 1 ;
    printf("removing cerebellum from atlas\n") ;
  }
  else if (!strcmp(option, "FONLY"))
  {
    ctl_point_fname = argv[2] ;
    nargs = 1 ;
    file_only = 1 ;
    printf("only using control points from file %s\n", ctl_point_fname) ;
  }
  else if (!strcmp(option, "SIGMA"))
  {
    bias_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("smoothing bias field with sigma = %2.1f\n", bias_sigma) ;
  }
  else if (!strcmp(option, "DILATE") || !strcmp(option, "DILATE_MASK"))
  {
    dilate_mask = atof(argv[2]) ;
    nargs = 1 ;
    printf("dilating brain mask %d times\n", dilate_mask) ;
  }
  else if (!strcmp(option, "DIAG"))
  {
    diag_fp = fopen(argv[2], "w") ;
    if (!diag_fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open diag file %s for writing",
                Progname, argv[2]) ;
    printf("opening diag file %s for writing\n", argv[2]) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  }
  else if (!strcmp(option, "DEBUG_NODE"))
  {
    Ggca_x = atoi(argv[2]) ;
    Ggca_y = atoi(argv[3]) ;
    Ggca_z = atoi(argv[4]) ;
    printf("debugging node (%d, %d, %d)\n", Ggca_x, Ggca_y, Ggca_z) ;
    nargs = 3 ;
  }
  else if (!strcmp(option, "TR"))
  {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  }
  else if (!strcmp(option, "EXAMPLE"))
  {
    example_T1 = argv[2] ;
    example_segmentation = argv[3] ;
    printf("using %s and %s as example T1 and segmentations respectively.\n",
           example_T1, example_segmentation) ;
    nargs = 2 ;
  }
  else if (!strcmp(option, "TE"))
  {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  }
  else if (!strcmp(option, "ALPHA"))
  {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  }
  else if (!strcmp(option, "NSAMPLES"))
  {
    normalized_transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing  transformed normalization control points to %s...\n",
           normalized_transformed_sample_fname) ;
  }
  else if (!stricmp(option, "RENORM") || !stricmp(option, "RENORMALIZE"))
  {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s...\n",
           renormalization_fname) ;
  }
  else if (!strcmp(option, "FLASH"))
  {
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  }
  else if (!strcmp(option, "PRIOR"))
  {
    min_prior = atof(argv[2]) ;
    nargs = 1 ;
    printf("using prior threshold %2.2f\n", min_prior) ;
  }
  else if (!stricmp(option, "NOVAR"))
  {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  }
  else switch (*option)
    {
    case 'R':
      read_ctrl_point_fname = argv[2] ;
      nargs = 1 ;
      printf("reading control point volume from %s\n", read_ctrl_point_fname) ;
      break ;
    case 'C':
      ctrl_point_fname = argv[2] ;
      nargs = 1 ;
      printf("writing control point volume to %s\n", ctrl_point_fname) ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'N':
      nregions = atoi(argv[2]) ;
      printf("using %d regions/struct for normalization\n", nregions) ;
      nargs = 1 ;
      break ;
    case 'F':
      ctl_point_fname = argv[2] ;
      nargs = 1 ;
      printf("reading manually defined control points from %s\n",
             ctl_point_fname) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      usage_exit(0) ;
      break ;
    case 'P':
      ctl_point_pct = atof(argv[2]) ;
      nargs = 1 ;
      printf("using top %2.1f%% wm points as control points....\n",
             100.0*ctl_point_pct) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static GCA_SAMPLE *
find_control_points
(GCA *gca, GCA_SAMPLE *gcas_total,
 int total_samples, int *pnorm_samples, int nregions, int label,
 MRI *mri_in, TRANSFORM *transform, double min_prior, double ctl_point_pct,
 char *sample_fname)
{
  int        i, j, *ordered_indices, nsamples,
    xmin, ymin, zmin, xmax, ymax, zmax, xv,yv,zv, nremoved,
             x, y, z, xi, yi, zi, region_samples,
             used_in_region, prior_wsize=5, image_wsize=3, histo_peak, n,
                             nbins ;
  GCA_SAMPLE *gcas, *gcas_region, *gcas_norm ;
  double     means[MAX_GCA_INPUTS], vars[MAX_GCA_INPUTS], val, outlying_nsigma = 3, nsigma ;
  HISTOGRAM  *histo, *hsmooth ;
  GC1D       *gc ;
  float      fmin, fmax ;
  MRI        *mri_T1 = NULL ;
  char       fname[STRLEN], fname_only[STRLEN] ;
  MRI        *mri_ctrl = NULL ;

#if 0
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s/../T1", mri_in->fname) ;
    mri_T1 = MRIread(fname) ;
#if 0
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE, "could not read T1 volume %s...", fname) ;
#endif
  }
#endif

  if (label == Gdiag_no)
    DiagBreak() ;

  if (sample_fname)
  {
    mri_ctrl = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 2) ;
    MRIcopyHeader(mri_in, mri_ctrl) ;
    FileNameRemoveExtension(sample_fname, fname_only) ;
  }
  MRIvalRange(mri_in, &fmin, &fmax) ;
  nbins = (int)(fmax-fmin+1);
  histo = HISTOalloc(nbins) ;
  hsmooth = HISTOalloc(nbins) ;
  for (nsamples = i = 0 ; i < total_samples ; i++)
  {
    if (gcas_total[i].label != label)
      continue ;

    nsamples++ ;
    if (mri_ctrl)
      MRIsetVoxVal(mri_ctrl, gcas_total[i].x, gcas_total[i].y, gcas_total[i].z, 0, 1) ;
  }

  *pnorm_samples = 0 ;
  printf("found %d control points for structure...\n", nsamples) ;

  if (nsamples == 0)
    return(NO_ERROR) ;

  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  gcas_region = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  gcas_norm = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  if (!gcas || !gcas_region || !gcas_norm)
    ErrorExit
    (ERROR_NOMEMORY,
     "find_control_points: could not allocate %d samples\n",nsamples);

  for (j = i = 0 ; i < total_samples ; i++)
  {
    if (gcas_total[i].label != label)
      continue ;

    memmove(&gcas[j], &gcas_total[i], sizeof(GCA_SAMPLE)) ;
    j++ ;
  }
  ordered_indices = (int *)calloc(nsamples, sizeof(int)) ;

  gcas_bounding_box
  (gcas, nsamples, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax, label) ;
  printf("bounding box (%d, %d, %d) --> (%d, %d, %d)\n",
         xmin, ymin, zmin, xmax, ymax, zmax) ;
  if (mri_ctrl)
  {
    for (x = 0 ; x < mri_ctrl->width ; x++)
      for (y = 0 ; y < mri_ctrl->height ; y++)
	for (z = 0 ; z < mri_ctrl->depth ; z++)
	{
	  GCA_PRIOR *gcap = getGCAP(gca, mri_ctrl, transform, x, y, z) ;
	  double prior ;
	  prior = getPrior(gcap, label) ;
	  MRIsetVoxVal(mri_ctrl, x, y, z, 1, prior) ;
	}
  }
  for (x = 0 ; x < nregions ; x++)
  {
    for (y = 0 ; y < nregions ; y++)
    {
      for (z = 0 ; z < nregions ; z++)
      {
        /* only process samples in this region */
        nsigma = 1.0 ;
        do
        {
          for (region_samples = i = 0 ; i < nsamples ; i++)
          {
            xi = (int)(nregions*(gcas[i].x - xmin) / (xmax-xmin+1)) ;
            yi = (int)(nregions*(gcas[i].y - ymin) / (ymax-ymin+1)) ;
            zi = (int)(nregions*(gcas[i].z - zmin) / (zmax-zmin+1)) ;
            if ((xi < 0 || xi >= nregions) ||
                (yi < 0 || yi >= nregions) ||
                (zi < 0 || zi >= nregions))
              DiagBreak() ;

            xv = gcas[i].x ; yv = gcas[i].y ; zv = gcas[i].z ;
            if (xv == Gx && yv == Gy && zv == Gz)
              DiagBreak() ;

            if (sqrt(SQR(xv-Gx)+SQR(yv-Gy)+SQR(zv-Gz)) < 2)
              DiagBreak() ;

            if (xi != x || yi != y || zi != z
                || gcas_getPrior(gcas[i]) < min_prior)
              continue ;

            if (min_region_prior
                (gca, gcas[i].xp, gcas[i].yp, gcas[i].zp,
                 prior_wsize, label) < .5)  // changed to .5
	    {
	      if (mri_ctrl)
		MRIsetVoxVal(mri_ctrl, xv, yv, zv, 0, 2) ;
              continue ;
	    }

	    if (gca->ninputs > 1)  // 1st input is norm.mgz and we can depend on its intensities
	    {
	      double val_T1 ;
	      MRIsampleVolumeFrame(mri_in, xv, yv, zv, 0, &val_T1) ;
	      if (val_T1 < MIN_WM_BIAS_PCT * DEFAULT_DESIRED_WHITE_MATTER_VALUE  ||
		  val_T1 > MAX_WM_BIAS_PCT  * DEFAULT_DESIRED_WHITE_MATTER_VALUE )
	      {
		if (mri_ctrl)
		  MRIsetVoxVal(mri_ctrl, xv, yv, zv, 0, 3) ;
		continue ;
	      }
	    }
	    else  // be more conservative
	    {
	      double val_T1 ;
	      MRIsampleVolumeFrame(mri_in, xv, yv, zv, 0, &val_T1) ;
	      if (val_T1 < .75*MIN_WM_BIAS_PCT * DEFAULT_DESIRED_WHITE_MATTER_VALUE  ||
		  val_T1 > 1.25*MAX_WM_BIAS_PCT  * DEFAULT_DESIRED_WHITE_MATTER_VALUE )
	      {
#if 1
		if (mri_ctrl)
		  MRIsetVoxVal(mri_ctrl, xv, yv, zv, 0, 3) ;
		continue ;
#endif
	      }
	    }
            if (uniform_region(gca, mri_in, transform,
                               xv, yv, zv,
                               image_wsize, &gcas[i], nsigma) == 0)
	    {
	      int xk, yk, zk, found = 0 ;

	      for (xk = -1 ; !found && xk <= 1 ; xk++)
	      {
		for (yk = -1 ; !found && yk <= 1 ; yk++)
		{
		  for (zk = -1 ; !found && zk <= 1 ; zk++)
		  {
		    if (uniform_region(gca, mri_in, transform,
				       xv+xk, yv+yk, zv+zk,
				       image_wsize, &gcas[i], nsigma) != 0)
		    {
		      // a neighboring voxel is uniform - move this sample over
		      found = 1 ;
		      xv += xk ; yv += yk ; zv += yk ;
		      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
			printf("moving control point (%d, %d, %d) to (%d, %d, %d) for uniformity\n",
			       xv, yv, zv, gcas[i].x, gcas[i].y, gcas[i].z) ;
		      gcas[i].x = xv ;  gcas[i].y = yv ;  gcas[i].z = zv ; 
		      break ;
		    }
		  }
		}
	      }

	      if (found == 0)
	      {
		if (mri_ctrl)
		  MRIsetVoxVal(mri_ctrl, xv, yv, zv, 0, 4) ;
		continue ;
	      }
	    }

	    if (mri_ctrl)
	      MRIsetVoxVal(mri_ctrl, xv, yv, zv, 0, 5) ;

            memmove(&gcas_region[region_samples],
                    &gcas[i],
                    sizeof(GCA_SAMPLE)) ;
            region_samples++ ;
            if (gcas[i].x == Gx &&
                gcas[i].y == Gy &&
                gcas[i].z == Gz)
            {
              DiagBreak() ;
            }
          }
          nsigma *= 1.1 ;
        }
        while (region_samples < 8 && nsigma < 3) ;

        if (region_samples < 8)/* can't reliably estimate statistics */
          continue ;

        if (DIAG_VERBOSE_ON)
          printf("\t%d total samples found in region (%d, %d, %d)\n",
                 region_samples,x, y,z) ;
        /* compute mean and variance of label within this region */
        for (n = 0 ; n < gca->ninputs ; n++)
        {
          HISTOclear(histo, histo) ; HISTOinit(histo, histo->nbins, 0, 255) ; HISTOinit(hsmooth, hsmooth->nbins, 0, 255);
          for (means[n] = vars[n] = 0.0, i = 0 ;
               i < region_samples ;
               i++)
          {
            MRIsampleVolumeFrame
            (mri_in,
             gcas_region[i].x,gcas_region[i].y,gcas_region[i].z,
             n, &val) ;
            if (FZERO(val))
            {
              if (i < (region_samples-1))
                memmove(&gcas_region[i],
                        &gcas_region[i+1],
                        (region_samples-(i+1))*sizeof(GCA_SAMPLE));
              i-- ;
              region_samples-- ;
              continue ;
            }
            histo->counts[(int)val]++ ;
            means[n] += val ;
            vars[n] += (val*val) ;
#if 0
            if (mri_T1)
            {
              val = MRIvox(mri_T1,
                           gcas_region[i].x,
                           gcas_region[i].y,
                           gcas_region[i].z) ;
              if (val < 85 || val > 130)
              {
                FILE *fp ;
                fp = fopen("badpoints.log", "a") ;
                fprintf(fp, "%s: (%d, %d, %d): %f\n",
                        mri_in->fname,
                        (int)gcas_region[i].x,
                        (int)gcas_region[i].y,
                        (int)gcas_region[i].z,val) ;
                fclose(fp) ;
                printf("!!!!!!!!!!!!!!!!!!!!!!! "
                       "%s: (%d, %d, %d): %f "
                       "!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
                       mri_in->fname,
                       (int)gcas_region[i].x,
                       (int)gcas_region[i].y,
                       (int)gcas_region[i].z,val) ;
              }
            }
#endif
          }

          HISTOsmooth(histo, hsmooth, 2) ;
          histo_peak =
            HISTOfindHighestPeakInRegion(hsmooth, 1, hsmooth->nbins) ;
          if (histo_peak < 0)   /* couldn't find a valid peak? */
          {
            break ;
          }

          for (means[n] = vars[n] = 0.0, i = 0 ;
               i < region_samples ;
               i++)
          {
            if (gcas_region[i].label < 0)
            {
              continue ;
            }
            MRIsampleVolumeFrame
            (mri_in,
             gcas_region[i].x,
             gcas_region[i].y,
             gcas_region[i].z,
             n, &val) ;
            means[n] += val ;
            vars[n] += (val*val) ;
          }
          means[n] /= (double)region_samples ;
          vars[n] =
            vars[n] / (double)region_samples - means[n]*means[n] ;

          means[n] = histo_peak ;
          if (DIAG_VERBOSE_ON)
            printf("\tlabel %s[%d]: %2.1f +- %2.1f\n",
                   cma_label_to_name(label),
                   n, means[n], sqrt(vars[n])) ;
        }

        /* ignore GCA mean and variance -
           use image instead (otherwise bias field will mess us up) */
        for (i = 0 ; i < region_samples ; i++)
        {
          int r ;

          for (r = 0 ; r < gca->ninputs ; r++)
          {
            gcas_region[i].means[r] = means[r] ;
          }
          /*          gcas_region[i].var = var ;*/
        }

        GCAcomputeLogSampleProbability
        (gca, gcas_region, mri_in, transform, region_samples, DEFAULT_CLAMP) ;
        GCArankSamples
        (gca, gcas_region, region_samples, ordered_indices) ;
#if 0
        /* use detected peak as normalization value for whole region */
        used_in_region = 1 ;
        j = ordered_indices[0] ;
        MRIvox(mri_in,
               gcas_region[j].x,
               gcas_region[j].y,
               gcas_region[j].z) = histo_peak ;
        memmove(&gcas_norm[*pnorm_samples],
                &gcas_region[j],
                sizeof(GCA_SAMPLE)) ;
        (*pnorm_samples)++ ;
#else
#if 1
        nremoved = GCAremoveOutlyingSamples
          (gca, gcas_region, mri_in, transform, region_samples, outlying_nsigma) ;
#endif
        for (used_in_region = i = 0 ; i < region_samples ; i++)
        {
          j = ordered_indices[i] ;
          if (gcas_region[j].label != label)  /* it was an outlier */
          {
	    if (mri_ctrl)
	      MRIsetVoxVal(mri_ctrl, gcas_region[i].x, gcas_region[i].y, gcas_region[i].z, 0, 1) ;
            continue ;
          }
          memmove
          (&gcas_norm[*pnorm_samples],
           &gcas_region[j],
           sizeof(GCA_SAMPLE)) ;
          (*pnorm_samples)++ ;
          used_in_region++ ;
        }
        if ((used_in_region <= 0) && region_samples>0)
        {
          j = ordered_indices[0] ;
          /*          gcas_region[j].label = label ;*/
          printf("forcing use of sample %d @ (%d, %d, %d)\n", j,
                 gcas_region[j].x,
                 gcas_region[j].y,
                 gcas_region[j].z) ;
          memmove(&gcas_norm[*pnorm_samples],
                  &gcas_region[j],
                  sizeof(GCA_SAMPLE)) ;
          (*pnorm_samples)++ ;
          used_in_region++ ;
        }
#endif
        if (DIAG_VERBOSE_ON)
        {
          printf("\t%d samples used in region\n", used_in_region) ;
        }
      }
    }
  }

  if (mri_ctrl)
  {
    int req = snprintf(fname, STRLEN, "%s_init_label%d.labels.mgz", fname_only, label) ;
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing initial sample points for %s to %s\n", cma_label_to_name(label), fname); 
    MRIwriteFrame(mri_ctrl, fname, 0) ;
    req = snprintf(fname, STRLEN, "%s_init_label%d.priors.mgz", fname_only, label) ;
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing initial label priors for %s to %s\n", cma_label_to_name(label), fname); 
    MRIwriteFrame(mri_ctrl, fname, 1) ;
    MRIfree(&mri_ctrl) ;
  }
  /* put gca means back into samples */
  for (i = 0 ; i < *pnorm_samples ; i++)
  {
    gc = GCAfindPriorGC(gca,
                        gcas_norm[i].xp,
                        gcas_norm[i].yp,
                        gcas_norm[i].zp,
                        gcas_norm[i].label) ;
    if (gc)
    {
      int r, c, v ;

      for (v = r = 0 ; r < gca->ninputs ; r++)
      {
        for (c = r ; c < gca->ninputs ; c++, v++)
        {
          gcas_norm[i].means[v] = gc->means[v] ;
          gcas_norm[i].covars[v] = gc->covars[v] ;
        }
      }
    }
  }
  HISTOfree(&histo) ;
  HISTOfree(&hsmooth) ;
  free(gcas_region) ;
  free(gcas) ;
  if (mri_T1)
  {
    MRIfree(&mri_T1) ;
  }
  return(gcas_norm) ;
}

static GCA_SAMPLE *
gcas_concatenate(GCA_SAMPLE *gcas1, GCA_SAMPLE *gcas2, int n1, int n2)
{
  GCA_SAMPLE *gcas ;
  int        i ;

  gcas = (GCA_SAMPLE *)calloc(n1+n2, sizeof(GCA_SAMPLE)) ;
  if (!gcas)
    ErrorExit
    (ERROR_NOMEMORY,
     "gcas_concatenate: could not allocate %d samples",n1+n2) ;

  for (i = 0 ; i < n1 ; i++)
  {
    memmove(&gcas[i], &gcas1[i], sizeof(GCA_SAMPLE)) ;
  }
  for (i = 0 ; i < n2 ; i++)
  {
    memmove(&gcas[i+n1], &gcas2[i], sizeof(GCA_SAMPLE)) ;
  }

  return(gcas) ;
}

static int
gcas_bounding_box(GCA_SAMPLE *gcas, int nsamples,
                  int *pxmin, int *pymin, int *pzmin,
                  int *pxmax, int *pymax, int *pzmax, int label)
{
  int   i, xmin, ymin, zmin, xmax, ymax, zmax ;

  xmax = ymax = zmax = -1 ;
  xmin = ymin = zmin = 1000000 ;
  for (i = 0 ; i < nsamples ; i++)
  {
    if (gcas[i].x < xmin)
    {
      xmin = gcas[i].x ;
    }
    if (gcas[i].y < ymin)
    {
      ymin = gcas[i].y ;
    }
    if (gcas[i].z < zmin)
    {
      zmin = gcas[i].z ;
    }

    if (gcas[i].x > xmax)
    {
      xmax = gcas[i].x ;
    }
    if (gcas[i].y > ymax)
    {
      ymax = gcas[i].y ;
    }
    if (gcas[i].z > zmax)
    {
      zmax = gcas[i].z ;
    }
  }

  *pxmin = xmin ;
  *pymin = ymin ;
  *pzmin = zmin ;
  *pxmax = xmax ;
  *pymax = ymax ;
  *pzmax = zmax ;
  return(NO_ERROR) ;
}

static double
min_region_prior(GCA *gca, int xp, int yp, int zp, int wsize, int label)
{
  int       whalf, xi, yi, zi, xk, yk, zk ;
  double    min_prior, prior ;
  GCA_PRIOR *gcap ;

  gcap = &gca->priors[xp][yp][zp] ;
  min_prior = getPrior(gcap, label) ;
  whalf = (wsize-1)/(gca->prior_spacing*2) ;
  for (xi = -whalf ; xi <= whalf ; xi++)
  {
    xk = xp+xi ;
    if (xk < 0 || xk >= gca->prior_width)
    {
      continue ;
    }
    for (yi = -whalf ; yi <= whalf ; yi++)
    {
      yk = yp+yi ;
      if (yk < 0 || yk >= gca->prior_height)
      {
        continue ;
      }
      for (zi = -whalf ; zi <= whalf ; zi++)
      {
        zk = zp+zi ;
        if (zk < 0 || zk >= gca->prior_depth)
        {
          continue ;
        }
        gcap = &gca->priors[xk][yk][zk] ;
        prior = getPrior(gcap, label) ;
        if (prior < min_prior)
        {
          min_prior = prior ;
        }
      }
    }
  }

  return(min_prior) ;
}

static int
uniform_region(GCA *gca, MRI *mri, TRANSFORM *transform,
               int x, int y, int z,
               int wsize, GCA_SAMPLE *gcas,
               float nsigma)
{
  int   xk, yk, zk, whalf, xi, yi, zi, n ;
  double val0, val, sigma, min_val,max_val, thresh ;
  MATRIX *m ;
  GC1D   *gc ;

  gc = GCAfindSourceGC(gca, mri, transform, x, y, z, gcas->label) ;
  if (!gc)
  {
    return(0) ;
  }
  m = load_covariance_matrix(gc, NULL, gca->ninputs) ;

  whalf = (wsize-1)/2 ;
  for (n = 0 ; n < gca->ninputs ; n++)
  {
    sigma = sqrt(*MATRIX_RELT(m, n+1, n+1)) ;
    MRIsampleVolumeFrame(mri, (double)x, (double)y, (double)z, n, &val0) ;
    if (sigma < 0.05*val0)   /* don't let it be too small */
    {
      sigma = 0.05*val0 ;
    }
    if (sigma > 0.1*val0)    /* don't let it be too big */
    {
      sigma = 0.1*val0 ;
    }
    min_val = max_val = val0 ;
    thresh = nsigma*sigma ;

    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      xi = mri->xi[x+xk] ;
      for (yk = -whalf ; yk <= whalf ; yk++)
      {
        yi = mri->yi[y+yk] ;
        for (zk = -whalf ; zk <= whalf ; zk++)
        {
          zi = mri->zi[z+zk] ;
          MRIsampleVolumeFrame
          (mri, (double)xi, (double)yi, (double)zi, n, &val) ;
          if (val < min_val)
          {
            min_val = val ;
          }
          if (val > max_val)
          {
            max_val = val ;
          }
          if (fabs(val-val0) > thresh ||
              fabs(max_val-min_val) > thresh)
          {
	    MatrixFree(&m) ;
            return(0) ;
          }
        }
      }
    }
  }

  MatrixFree(&m) ;
  return(1) ;
}

static int
discard_unlikely_control_points(GCA *gca, GCA_SAMPLE *gcas, int nsamples,
                                MRI *mri_in, TRANSFORM *transform, const char *name)
{
  int    i, xv, yv, zv, n, peak, start, end, num ;
  HISTO *h, *hsmooth ;
  float  fmin, fmax ;
  double val,  mean_ratio, min_T1, max_T1 ;

  if (nsamples == 0)
    return(NO_ERROR) ;

  min_T1 = (1.0-extra_norm_range) * MIN_WM_BIAS_PCT * DEFAULT_DESIRED_WHITE_MATTER_VALUE  ;
  max_T1 = (1.0+extra_norm_range) * MAX_WM_BIAS_PCT* DEFAULT_DESIRED_WHITE_MATTER_VALUE  ;
  for (num = n = 0 ; n < gca->ninputs ; n++)
  {
    int niter = 0 ;
    MRIvalRangeFrame(mri_in, &fmin, &fmax, n) ;
    h = HISTOinit(NULL, nint(fmax-fmin)+1, fmin, fmax) ;

    for (i = 0 ; i < nsamples ; i++)
    {
      xv = gcas[i].x ; yv = gcas[i].y ; zv = gcas[i].z ;
      if (xv == Gx && yv == Gy && zv == Gz)
        DiagBreak() ;

      MRIsampleVolumeFrame(mri_in, gcas[i].x,gcas[i].y,gcas[i].z, n, &val) ;
      if (FZERO(val))
        DiagBreak() ;

      if (n >= 1)
      {
	double val_T1 ;
	MRIsampleVolumeFrame(mri_in, gcas[i].x,gcas[i].y,gcas[i].z, 0, &val_T1) ;
	if (val_T1 < min_T1 || val_T1 > max_T1)
	  continue ;
      }
      h->counts[nint(val-fmin)]++ ;
    }

    /* check to see  if peak is unlikely */
    hsmooth = HISTOsmooth(h, NULL, 2) ;
    do
    {
      if (gca->ninputs == 1) /*  find  brightest peak as
                                 for  n=1 it is T1  weighted  */
      {
        peak = HISTOfindLastPeak(hsmooth, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
      }
      else
      {
        peak = HISTOfindHighestPeakInRegion(hsmooth, 0, h->nbins-1) ;
      }
      end = HISTOfindEndOfPeak(hsmooth, peak, 0.01) ;
      start = HISTOfindStartOfPeak(hsmooth, peak, 0.01) ;
      for (mean_ratio = 0.0, i = 0 ; i < nsamples ; i++)
      {
        mean_ratio += hsmooth->bins[peak] / gcas[i].means[n];
      }
      mean_ratio /= (double)nsamples ;
      HISTOclearBins(hsmooth, hsmooth, hsmooth->bins[start], hsmooth->bins[end])  ;
      if (niter++ > 5)
      {
        break ;
      }
      if (niter > 1)
      {
        DiagBreak() ;
      }
    }
    while  (mean_ratio  < 0.5 || mean_ratio > 2.0) ;

    if (fmin+start <  min_T1)
      start = min_T1-fmin ;
    if (fmin+end > max_T1)
      end = max_T1-fmin ;

    printf("%s: limiting intensities to %2.1f --> %2.1f\n",
           name, fmin+start, fmin+end) ;
    for (i = 0 ; i < nsamples ; i++)
    {
      xv = gcas[i].x ; yv = gcas[i].y ; zv = gcas[i].z ;
      if (xv == Gx && yv == Gy && zv == Gz)
        DiagBreak() ;

      MRIsampleVolumeFrame(mri_in,gcas[i].x,gcas[i].y,gcas[i].z,n,&val) ;
      if (val-fmin < start || val-fmin > end)
      {
        num++ ;
        gcas[i].label = 0 ;
      }
      else if (gca->ninputs > 1)
      {
	MRIsampleVolumeFrame(mri_in,gcas[i].x,gcas[i].y,gcas[i].z,0,&val) ;
	if (val < min_T1 || val > max_T1)  // check norm.mgz to make sure it is in wm
	{
	  num++ ;
	  gcas[i].label = 0 ;
	}
      }
    }
    HISTOfree(&h) ;
    HISTOfree(&hsmooth) ;
  }

  printf("%d of %d (%2.1f%%) samples deleted\n",
         num, nsamples, 100.0f*(float)num/(float)nsamples) ;
  return(NO_ERROR) ;
}

#include "mri_ca_normalize.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_ca_normalize_help_xml,
                mri_ca_normalize_help_xml_len);
  exit(code) ;
}

static MRI *
normalize_from_segmentation_volume
(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, int *structs, int nstructs)
{
  MRI *mri_bin, *mri_tmp ;
  int i, x, y, z, label, ctrl ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }
  mri_tmp = MRIclone(mri_seg, NULL) ;
  mri_bin = MRIclone(mri_seg, NULL) ;
  for (i = 0 ; i < nstructs ; i++)
  {
    MRIcopyLabel(mri_seg, mri_tmp, structs[i]) ;
    MRIbinarize(mri_tmp, mri_tmp, 1, 0, 1) ;
    MRIerode(mri_tmp, mri_tmp) ;
    MRIerode(mri_tmp, mri_tmp) ;
    MRIadd(mri_tmp, mri_bin, mri_bin) ;
    MRIclear(mri_tmp) ;
  }

  MRIopen(mri_bin, mri_bin) ;
  //  mri_dst = normalizeFromLabel(mri_src, mri_dst, mri_bin, fas) ;
  //normalize each channel separately

  for(i=0; i < mri_src->nframes; i++)
  {
    mri_dst = normalizeChannelFromLabel(mri_src, mri_dst, mri_bin, fas, i);
  }

  MRInormGentlyFindControlPoints(mri_dst, 110, 20, 10, mri_bin, NULL) ;
  // remove control points that don't agree with the seg
  for (x = 0 ; x < mri_dst->width ; x++)
    for (y = 0 ; y < mri_dst->height ; y++)
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        ctrl = (int)nint(MRIgetVoxVal(mri_bin, x, y, z, 0)) ;
        if (ctrl == 0)
        {
          continue ;
        }
        label = (int)nint(MRIgetVoxVal(mri_seg, x, y, z, 0)) ;
        switch (label)
        {
        default:
          MRIsetVoxVal(mri_bin, x, y, z, 0, 0) ;
          break ;
        case Left_Cerebral_White_Matter:
        case Right_Cerebral_White_Matter:
        case Left_Cerebellum_White_Matter:
        case Right_Cerebellum_White_Matter:
        case Brain_Stem:
        case Left_VentralDC:
        case Right_VentralDC:
          break ;   // these are ok
        }
      }



  for(i=0; i < mri_src->nframes; i++)
  {
    normalizeChannelFromLabel(mri_dst, mri_dst, mri_bin, fas, i);
  }

  MRIfree(&mri_bin) ;
  MRIfree(&mri_tmp) ;
  return(mri_dst) ;
}
#include "mrinorm.h"
MRI *
normalizeFromLabel(MRI *mri_in, MRI *mri_dst, MRI *mri_seg, double *fas)
{
  MRI    *mri_ctrl, *mri_bias ;
  int    x, y, z, width, height, depth, num, total, input, T1_index, i ;
  float   bias ;
  double  mean, sigma, max_fa ;
  double  val ;

  max_fa = fas[T1_index = 0] ;
  for (i = 1 ; i < mri_in->nframes ; i++)
    if (fas[i] > max_fa)
    {
      T1_index = i ;
      max_fa = fas[i] ;
    }
  printf("using volume %d as most T1-weighted for normalization\n",T1_index) ;
  width = mri_in->width ;
  height = mri_in->height ;
  depth = mri_in->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_in, NULL) ;
  }
  mri_ctrl = MRIalloc(width, height, depth, MRI_UCHAR) ;
  MRIcopyHeader(mri_in, mri_ctrl);
  mri_bias = MRIalloc(mri_in->width,mri_in->height,mri_in->depth,MRI_SHORT);
  if (!mri_bias)
    ErrorExit
    (ERROR_NOMEMORY,
     "GCAnormalizeSamples: could not allocate (%d,%d,%d,2) bias image",
     mri_in->width,mri_in->height,mri_in->depth) ;
  MRIcopyHeader(mri_in, mri_bias);

#define NO_BIAS  1000

  /* use all non-zero locations in mri_seg as control points */
  MRIbinarize(mri_seg, mri_ctrl, 1, 0, CONTROL_MARKED) ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        MRISvox(mri_bias, x,y,z) = NO_BIAS ;  /* by default */
        if (MRIvox(mri_ctrl, x, y, z) !=
            CONTROL_MARKED)  /* not read from file */
        {
          continue ;
        }

        MRIsampleVolumeFrame(mri_in, x, y, z, T1_index, &val) ;
        bias = NO_BIAS*DEFAULT_DESIRED_WHITE_MATTER_VALUE / val ;
        MRISvox(mri_bias, x, y, z) = (short)nint(bias) ;
      }
    }
  }

  /* now check for and remove outliers */
  mean = sigma = 0.0 ;
  for (num = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED)
        {
          num++ ;
          bias = (double)MRISvox(mri_bias, x, y, z) ;
          mean += bias ;
          sigma += (bias*bias) ;
        }
      }
    }
  }

  if (num > 0)
  {
    mean /= (double)num ;
    sigma  = sqrt(sigma / (double)num - mean*mean) ;
    printf("bias field = %2.3f +- %2.3f\n", mean/NO_BIAS, sigma/NO_BIAS) ;
  }

  /* now check for and remove outliers */
  for (total = num = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED)
        {
          bias = (double)MRISvox(mri_bias, x, y, z) ;
          total++ ;
          if (fabs(bias-mean) > 4*sigma)
          {
            MRIvox(mri_ctrl, x, y, z) = CONTROL_NONE ;
            num++ ;
            MRISvox(mri_bias, x, y, z) = NO_BIAS ;
          }
        }
      }
    }
  }

  printf("%d of %d control points discarded\n", num, total) ;

  MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias) ;
  /*  MRIwrite(mri_bias, "bias.mgh") ;*/
#if 1
  {
    MRI *mri_kernel, *mri_smooth, *mri_down ;
    float sigma = 16.0f ;

    mri_down = MRIdownsample2(mri_bias, NULL) ;
    mri_kernel = MRIgaussian1d(sigma, 100) ;
    mri_smooth = MRIconvolveGaussian(mri_down, NULL, mri_kernel) ;
    MRIfree(&mri_bias) ;
    MRIfree(&mri_kernel) ;
    mri_bias = MRIupsample2(mri_smooth, NULL) ;
    sigma = 2.0f ;
    MRIfree(&mri_down) ;
    MRIfree(&mri_smooth) ;
    mri_kernel = MRIgaussian1d(sigma, 100) ;
    mri_smooth = MRIconvolveGaussian(mri_bias, NULL, mri_kernel) ;
    MRIfree(&mri_bias) ;
    mri_bias = mri_smooth ;
    MRIfree(&mri_kernel) ;
  }
#else
  MRIsoapBubble(mri_bias, mri_ctrl, mri_bias, 10, -1) ;
#endif
  /*  MRIwrite(mri_bias, "smooth_bias.mgh") ;*/

  width = mri_in->width ;
  height = mri_in->height ;
  depth = mri_in->depth ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        bias = (float)MRISvox(mri_bias, x, y, z)/NO_BIAS ;
        if (bias < 0)
        {
          DiagBreak() ;
        }
        for (input = 0 ; input < mri_in->nframes ; input++)
        {
          MRIsampleVolumeFrame(mri_in, x, y, z, input, &val) ;
          val *= bias ;   /* corrected value */
          switch (mri_in->type)
          {
          case MRI_UCHAR:
            if (val < 0)
            {
              val = 0 ;
            }
            else if (val > 255)
            {
              val = 255 ;
            }
            MRIseq_vox(mri_dst, x, y, z, input) =
              (BUFTYPE)nint(val) ;
            break ;
          case MRI_SHORT:
            MRISseq_vox(mri_dst, x, y, z, input) =
              (short)nint(val) ;
            break ;
          case MRI_FLOAT:
            MRIFseq_vox(mri_dst, x, y, z, input) =
              val ;
            break ;
          default:
            ErrorReturn
            (NULL,
             (ERROR_UNSUPPORTED,
              "GCAnormalizeSamples: unsupported input type %d",
              mri_in->type));
            break ;
          }
        }
      }
    }
  }

  MRIfree(&mri_bias) ;
  MRIfree(&mri_ctrl) ;
  return(mri_dst) ;
}


MRI *
normalizeChannelFromLabel(MRI *mri_in, MRI *mri_dst, MRI *mri_seg,
                          double *fas, int input_index)
{
  MRI    *mri_ctrl, *mri_bias ;
  int    x, y, z, width, height, depth, num, total;
  float   bias ;
  double  mean, sigma;
  double  val ;

  width = mri_in->width ;
  height = mri_in->height ;
  depth = mri_in->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_in, NULL) ;
  }
  mri_ctrl = MRIalloc(width, height, depth, MRI_UCHAR) ;
  MRIcopyHeader(mri_in, mri_ctrl);
  mri_bias = MRIalloc(mri_in->width,mri_in->height,mri_in->depth,MRI_SHORT);
  if (!mri_bias)
    ErrorExit
    (ERROR_NOMEMORY,
     "GCAnormalizeSamples: could not allocate (%d,%d,%d,2) bias image",
     mri_in->width,mri_in->height,mri_in->depth) ;
  MRIcopyHeader(mri_in, mri_bias);

#define NO_BIAS  1000

  /* use all non-zero locations in mri_seg as control points */
  MRIbinarize(mri_seg, mri_ctrl, 1, 0, CONTROL_MARKED) ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        MRISvox(mri_bias, x,y,z) = NO_BIAS ;  /* by default */
        if (MRIvox(mri_ctrl, x, y, z) !=
            CONTROL_MARKED)  /* not read from file */
        {
          continue ;
        }

        MRIsampleVolumeFrame(mri_in, x, y, z, input_index, &val) ;
        bias = NO_BIAS*DEFAULT_DESIRED_WHITE_MATTER_VALUE / val ;
        MRISvox(mri_bias, x, y, z) = (short)nint(bias) ;
      }
    }
  }

  /* now check for and remove outliers */
  mean = sigma = 0.0 ;
  for (num = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED)
        {
          num++ ;
          bias = (double)MRISvox(mri_bias, x, y, z) ;
          mean += bias ;
          sigma += (bias*bias) ;
        }
      }
    }
  }

  if (num > 0)
  {
    mean /= (double)num ;
    sigma  = sqrt(sigma / (double)num - mean*mean) ;
    printf("bias field = %2.3f +- %2.3f\n", mean/NO_BIAS, sigma/NO_BIAS) ;
  }

  /* now check for and remove outliers */
  for (total = num = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED)
        {
          bias = (double)MRISvox(mri_bias, x, y, z) ;
          total++ ;
          if (fabs(bias-mean) > 2*sigma)
          {
            MRIvox(mri_ctrl, x, y, z) = CONTROL_NONE ;
            num++ ;
            MRISvox(mri_bias, x, y, z) = NO_BIAS ;
          }
        }
      }
    }
  }

  printf("%d of %d control points discarded\n", num, total) ;

  MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias) ;
  /*  MRIwrite(mri_bias, "bias.mgh") ;*/
#if 1
  {
    MRI *mri_kernel, *mri_smooth, *mri_down ;
    float sigma = 16.0f ;

    mri_down = MRIdownsample2(mri_bias, NULL) ;
    mri_kernel = MRIgaussian1d(sigma, 100) ;
    mri_smooth = MRIconvolveGaussian(mri_down, NULL, mri_kernel) ;
    MRIfree(&mri_bias) ;
    MRIfree(&mri_kernel) ;
    mri_bias = MRIupsample2(mri_smooth, NULL) ;
    sigma = 2.0f ;
    MRIfree(&mri_down) ;
    MRIfree(&mri_smooth) ;
    mri_kernel = MRIgaussian1d(sigma, 100) ;
    mri_smooth = MRIconvolveGaussian(mri_bias, NULL, mri_kernel) ;
    MRIfree(&mri_bias) ;
    mri_bias = mri_smooth ;
    MRIfree(&mri_kernel) ;
  }
#else
  MRIsoapBubble(mri_bias, mri_ctrl, mri_bias, 10, -1) ;
#endif
  /*  MRIwrite(mri_bias, "smooth_bias.mgh") ;*/


  width = mri_in->width ;
  height = mri_in->height ;
  depth = mri_in->depth ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        bias = (float)MRISvox(mri_bias, x, y, z)/NO_BIAS ;
        if (bias < 0)
        {
          DiagBreak() ;
        }
        {
          MRIsampleVolumeFrame(mri_in, x, y, z, input_index, &val) ;
          val *= bias ;   /* corrected value */
          switch (mri_in->type)
          {
          case MRI_UCHAR:
            if (val < 0)
            {
              val = 0 ;
            }
            else if (val > 255)
            {
              val = 255 ;
            }
            MRIseq_vox(mri_dst, x, y, z, input_index) =
              (BUFTYPE)nint(val) ;
            break ;
          case MRI_SHORT:
            MRISseq_vox(mri_dst, x, y, z, input_index) =
              (short)nint(val) ;
            break ;
          case MRI_FLOAT:
            MRIFseq_vox(mri_dst, x, y, z, input_index) =
              val ;
            break ;
          default:
            ErrorReturn
            (NULL,
             (ERROR_UNSUPPORTED,
              "GCAnormalizeSamples: unsupported input type %d",
              mri_in->type));
            break ;
          }
        }
      }
    }
  }

  MRIfree(&mri_bias) ;
  MRIfree(&mri_ctrl) ;
  return(mri_dst) ;
}
static int
copy_ctrl_points_to_volume(GCA_SAMPLE *gcas, int nsamples, MRI *mri_ctrl, int frame)
{
  int   i, xv, yv, zv, nregions ;

  nregions = mri_ctrl->nframes/2 ;
  for (i = 0 ; i < nsamples ; i++)
  {
    xv = gcas[i].x ;
    yv = gcas[i].y ;
    zv = gcas[i].z ;
    if (xv == Ggca_x && yv == Ggca_y && zv == Ggca_z)
    {
      DiagBreak() ;
    }
    MRIsetVoxVal(mri_ctrl, xv, yv, zv, frame, gcas[i].label) ;
    MRIsetVoxVal(mri_ctrl, xv, yv, zv, frame+nregions, gcas[i].means[0]) ;
  }

  return(NO_ERROR) ;
}


static GCA_SAMPLE *
copy_ctrl_points_from_volume(GCA *gca, TRANSFORM *transform, int *pnsamples,
                             MRI *mri_ctrl, int frame)
{
  GCA_SAMPLE *gcas ;
  int        nsamples, x, y, z, label, xp, yp, zp, nregions ;

  for (nsamples = x = 0 ; x < mri_ctrl->width ; x++)
  {
    for (y = 0 ; y < mri_ctrl->height ; y++)
    {
      for (z = 0 ; z < mri_ctrl->depth ; z++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_ctrl, x, y, z, frame) ;
        if (label > 0)
        {
          nsamples++ ;
        }
      }
    }
  }

  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  if (!gcas)
    ErrorExit
    (ERROR_NOMEMORY,
     "copy_ctrl_points_from_volume: could not allocate %d samples",nsamples) ;
  *pnsamples = nsamples ;
  nregions = mri_ctrl->nframes/2 ;
  for (nsamples = x = 0 ; x < mri_ctrl->width ; x++)
  {
    for (y = 0 ; y < mri_ctrl->height ; y++)
    {
      for (z = 0 ; z < mri_ctrl->depth ; z++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_ctrl, x, y, z, frame) ;
        if (label > 0)
        {
          gcas[nsamples].x = x ;
          gcas[nsamples].y = y ;
          gcas[nsamples].z = z ;
          gcas[nsamples].means = (float*)calloc(1, sizeof(float)) ;
          if (gcas[nsamples].means == NULL)
            ErrorExit(ERROR_NOMEMORY,
                      "%s: could not allocated %dth sample mean",
                      Progname, nsamples) ;
          gcas[nsamples].means[0] =
            MRIgetVoxVal(mri_ctrl, x, y, z, frame+nregions) ;
          gcas[nsamples].label = label ;


          if (!GCAsourceVoxelToPrior(gca, mri_ctrl,
                                     transform,
                                     x, y, z, &xp, &yp, &zp))
          {
            gcas[nsamples].xp = xp ;
            gcas[nsamples].yp = yp ;
            gcas[nsamples].zp = zp ;
          }
          nsamples++ ;
        }
      }
    }
  }
  return(gcas) ;
}

static int
fill_in_sample_means(GCA_SAMPLE *gcas, GCA *gca, int nsamples)
{
  int n ;
  GC1D  *gc ;


  for (n = 0 ; n < nsamples ; n++)

  {
    gc = GCAfindPriorGC(gca,
                        gcas[n].xp,
                        gcas[n].yp,
                        gcas[n].zp,
                        gcas[n].label) ;
    if (gc)
    {
      int r, c, v ;

      for (v = r = 0 ; r < gca->ninputs ; r++)
      {
        for (c = r ; c < gca->ninputs ; c++, v++)
        {
          gcas[n].means[v] = gc->means[v] ;
          //          gcas[n].covars[v] = gc->covars[v] ;
        }
      }
    }
  }

  return(NO_ERROR) ;
}

