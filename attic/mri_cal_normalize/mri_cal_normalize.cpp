/**
 * @brief Normalize a set of longituindal volumes making use of subcortical atlas data
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


const char *Progname ;


static int remove_cerebellum = 0 ;
static int remove_lh = 0 ;
static int remove_rh = 0 ;

static char *ctl_point_fname = NULL ;
static char *sample_fname = NULL ;
static char *normalized_transformed_sample_fname = NULL ;

static int file_only = 0 ;
static char *mask_fname = NULL ;
static int novar = 0 ;

static double bias_sigma = 4.0 ;
static float min_prior = 0.6 ;
static FILE *diag_fp = NULL ;

static MRI *scale_all_images(MRI *mri_in, MRI *mri_out) ;
static MRI *normalize_timepoints_with_parzen_window(MRI *mri_in, MRI *mri_out, double cross_time_sigma) ;
//static int normalize_timepoints_with_samples(MRI *mri, GCA_SAMPLE *gcas, int nsamples, int nsoap) ;
static int normalize_timepoints(MRI *mri, double thresh, double cross_time_sigma) ;
static void usage_exit(int code) ;
static int get_option(int argc, char *argv[]) ;
static int copy_ctrl_points_to_volume(GCA_SAMPLE *gcas, int nsamples, 
                                      MRI *mri_ctrl, int frame) ;
static int discard_control_points_with_different_labels(GCA_SAMPLE *gcas, int nsamples, MRI *mri_aseg) ;

static const char *aseg_fname = "aseg.mgz" ;
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
 double ctrl_point_pct) ;

static GCA_SAMPLE *gcas_concatenate(GCA_SAMPLE *gcas1, GCA_SAMPLE *gcas2, int n1, int n2);
static int  gcas_bounding_box(GCA_SAMPLE *gcas, int nsamples, int *pxmin, int *pymin, int *pzmin,
                              int *pxmax, int *pymax, int *pzmax, int label) ;
static int  uniform_region(GCA *gca, MRI *mri, TRANSFORM *transform,
                           int x, int y, int z, int wsize, GCA_SAMPLE *gcas, float nsigma) ;
static int  discard_unlikely_control_points(GCA *gca, GCA_SAMPLE *gcas_struct, int struct_samples,
                                            MRI *mri_in, TRANSFORM *transform, const char *name) ;

/*
  command line consists of these inputs:

  argv[1]  - base tps file
  argv[2]  - input volume name 
  argv[3]  - atlas (gca)
  argv[4]  - transform (lta/xfm/m3d)
  argv[5]  - output volume name
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

static double cross_time_sigma = 1.0 ;

static int nregions = 3 ;  /* divide each struct into 3x3x3 regions */

static char *ctrl_point_fname = NULL ;
static char *read_ctrl_point_fname = NULL ;

static int longinput = 0;

#define MAX_TIMEPOINTS 2000
static char *subjects[MAX_TIMEPOINTS] ;
int
main(int argc, char *argv[])
{
  char         *gca_fname, *in_fname, *out_fname, **av, *xform_fname, fname[STRLEN] ;
  MRI          *mri_in, *mri_norm = NULL, *mri_tmp, *mri_ctrl = NULL, *mri_aseg = NULL ;
  GCA          *gca ;
  int          ac, nargs, nsamples, msec, minutes, seconds;
  int          i, struct_samples, norm_samples = 0, n, input, ninputs ;
  Timer start ;
  GCA_SAMPLE   *gcas, *gcas_norm = NULL, *gcas_struct ;
  TRANSFORM    *transform = NULL ;
  char         line[STRLEN], *cp, sdir[STRLEN], base_name[STRLEN] ;
  FILE         *fp ;

  std::string cmdline = getAllInfo(argc, argv, "mri_cal_normalize");

  nargs = handleVersionOption(argc, argv, "mri_cal_normalize");
  if (nargs && argc - nargs == 1)
    exit (0);
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

  if (argc < 6)
    ErrorExit
      (ERROR_BADPARM,
       "usage: %s [<options>] <longitudinal time point file> <in vol name> <atlas> <transform file> <out vol name> \n",
       Progname) ;
  in_fname = argv[2] ;
  gca_fname = argv[3] ;
  xform_fname = argv[4] ;
  out_fname = argv[5] ;

  transform = TransformRead(xform_fname) ;
  if (transform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, xform_fname) ;
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

  gca = GCAread(gca_fname) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",Progname, gca_fname) ;
  GCAregularizeConditionalDensities(gca, .5) ;

  FileNamePath(argv[1], sdir) ;
  cp = strrchr(sdir, '/') ; 
  if (cp)
  {
    strcpy(base_name, cp+1) ;
    *cp = 0 ;  // remove last component of path, which is base subject name
  }
  ninputs = 0 ;
  fp = fopen(argv[1], "r") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read time point file %s", Progname, argv[1]) ;

  do
  {
    cp = fgetl(line, STRLEN-1, fp) ;
    if (cp != NULL && strlen(cp) > 0)
    {
      subjects[ninputs] = (char *)calloc(strlen(cp)+1, sizeof(char)) ;
      strcpy(subjects[ninputs], cp) ;
      ninputs++ ;
    }
  } while (cp != NULL && strlen(cp) > 0) ;
  fclose(fp) ;
  printf("processing %d timepoints in SUBJECTS_DIR %s...\n", ninputs, sdir) ;
  for (input = 0 ; input < ninputs ; input++)
  {
    printf("reading subject %d of %d: %s \n", input+1, ninputs, subjects[input]) ;
    if (longinput)
      sprintf(fname, "%s/%s.long.%s/mri/%s", sdir, subjects[input], base_name, in_fname) ;
    else
      sprintf(fname, "%s/%s/longtp/%s/%s", sdir, base_name, subjects[input], in_fname) ;
    
    mri_tmp = MRIread(fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s . If you are trying to run this on top of existing longitudinal data, switch on -longinput.",
                Progname, fname) ;
    MRImakePositive(mri_tmp, mri_tmp) ;
    if (mri_tmp && ctrl_point_fname && !mri_ctrl)
    {
      mri_ctrl = MRIallocSequence(mri_tmp->width, mri_tmp->height, 
                                  mri_tmp->depth,MRI_FLOAT, nregions*2) ; // labels and means
      MRIcopyHeader(mri_tmp, mri_ctrl) ;
    }
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
      mri_aseg = MRIclone(mri_in, NULL) ;
    }

    if (mask_fname)
    {
      int i ;
      MRI *mri_mask ;

      mri_mask = MRIread(mask_fname) ;
      if (!mri_mask)
        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                  Progname, mask_fname) ;

      for (i = 1 ; i < WM_MIN_VAL ; i++)
        MRIreplaceValues(mri_mask, mri_mask, i, 0) ;
      MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
      MRIfree(&mri_mask) ;
    }
    MRIcopyFrame(mri_tmp, mri_in, 0, input) ;
    MRIfree(&mri_tmp) ;
    if (longinput)
      sprintf(fname, "%s/%s.long.%s/mri/%s", sdir, subjects[input], base_name, aseg_fname) ;
    else
      sprintf(fname, "%s/%s/longtp/%s/%s", sdir, base_name, subjects[input], aseg_fname) ;     
    mri_tmp = MRIread(fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s", Progname, fname) ;
    MRIcopyFrame(mri_tmp, mri_aseg, 0, input) ;
    MRIfree(&mri_tmp) ;
  }
  MRIaddCommandLine(mri_in, cmdline) ;

//  GCAhistoScaleImageIntensitiesLongitudinal(gca, mri_in, 1) ;

// don't do this (there are still bugs, I think. Instead use robust template with --iscaleonly if necessary at all).
if (0) 
  scale_all_images(mri_in, mri_in) ;


  {
    int j ;

    gcas = GCAfindAllSamples(gca, &nsamples, NULL, 1) ;
    printf("using %d sample points...\n", nsamples) ;
    GCAcomputeSampleCoords(gca, mri_in, gcas, nsamples, transform) ;
    if (sample_fname)
      GCAtransformAndWriteSamples
        (gca, mri_in, gcas, nsamples, sample_fname, transform) ;

    for (j = 0 ; j < 1 ; j++)
    {
      for (n = 1 ; n <= nregions ; n++)
      {
        for (norm_samples = i = 0 ; i < NSTRUCTURES ; i++)
        {
          if (normalization_structures[i] == Gdiag_no)
            DiagBreak() ;
          printf("finding control points in %s....\n",
                 cma_label_to_name(normalization_structures[i])) ;
          gcas_struct = find_control_points(gca, gcas, nsamples, &struct_samples, n,
                                            normalization_structures[i], mri_in, transform, min_prior,
                                            ctl_point_pct) ;
          discard_unlikely_control_points(gca, gcas_struct, struct_samples, mri_in, transform,
                                          cma_label_to_name(normalization_structures[i])) ;
          discard_control_points_with_different_labels(gcas_struct, struct_samples, mri_aseg) ;
          if (mri_ctrl && ctrl_point_fname) // store the samples
            copy_ctrl_points_to_volume(gcas_struct, struct_samples, mri_ctrl, n-1) ;
          if (i)
          {
            GCA_SAMPLE *gcas_tmp ;
            gcas_tmp = gcas_concatenate(gcas_norm, gcas_struct, norm_samples, struct_samples) ;
            free(gcas_norm) ;
            norm_samples += struct_samples ;
            gcas_norm = gcas_tmp ;
          }
          else
          {
            gcas_norm = gcas_struct ; norm_samples = struct_samples ;
          }
        }

	if (norm_samples == 0)
	{
	  printf("could not find control points for region %d\n", n) ;
	  continue ;
	}
        printf("using %d total control points "
                 "for intensity normalization...\n", norm_samples) ;
        if (normalized_transformed_sample_fname)
          GCAtransformAndWriteSamples(gca, mri_in, gcas_norm, norm_samples,
                                      normalized_transformed_sample_fname,
                                      transform) ;
        mri_norm = GCAnormalizeSamplesAllChannels(mri_in, gca, gcas_norm, file_only ? 0 :norm_samples,
                                                  transform, ctl_point_fname, bias_sigma) ;
        if (Gdiag & DIAG_WRITE)
        {
          char fname[STRLEN] ;
          sprintf(fname, "norm%d.mgz", n) ;
          printf("writing normalized volume to %s\n", fname) ;
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

  // now do cross-time normalization to bring each timepoint closer to the mean at each location
#if 0
  // seems to hurt more than it helps in some cases
  if (mri_in->nframes > 1 && norm_samples > 0)
  {
    MRI   *mri_frame1, *mri_frame2, *mri_tmp ;
    double rms_before, rms_after ;
    int    i ;

    mri_tmp = MRIcopy(mri_in, NULL) ;
    mri_frame1 = MRIcopyFrame(mri_in, NULL, 0, 0) ;
    mri_frame2 = MRIcopyFrame(mri_in, NULL, 1, 0) ;
    rms_before = MRIrmsDiff(mri_frame1, mri_frame2) ;
    printf("RMS before = %2.2f\n", rms_before) ;
    MRIfree(&mri_frame1) ; MRIfree(&mri_frame2) ;
    for (i = 50 ; i <= 50 ; i += 25)
    {
      MRIcopy(mri_tmp, mri_in) ;
      normalize_timepoints_with_samples(mri_in, gcas_norm, norm_samples, i) ;
      mri_frame1 = MRIcopyFrame(mri_in, NULL, 0, 0) ;
      mri_frame2 = MRIcopyFrame(mri_in, NULL, 1, 0) ;
      rms_after = MRIrmsDiff(mri_frame1, mri_frame2) ;
      printf("RMS after (%d) = %2.2f\n", i, rms_after) ;
      MRIfree(&mri_frame1) ; MRIfree(&mri_frame2) ;
    }
  }
#endif

  // try to bring the images closer to each other at each voxel where they seem to come from the same distribution
  if (mri_in->nframes > 1)
  {
    MRI   *mri_frame1, *mri_frame2 ;
    double rms_after ;

    if (0)
      normalize_timepoints(mri_in, 2.0, cross_time_sigma) ;
    else
    {
      MRI * mri_temp = normalize_timepoints_with_parzen_window(mri_in, NULL, cross_time_sigma) ;
      MRIfree(&mri_in);
      mri_in = mri_temp;
    }
      
    mri_frame1 = MRIcopyFrame(mri_in, NULL, 0, 0) ;
    mri_frame2 = MRIcopyFrame(mri_in, NULL, 1, 0) ;
    rms_after = MRIrmsDiff(mri_frame1, mri_frame2) ;
    MRIfree(&mri_frame1) ; MRIfree(&mri_frame2) ;
    printf("RMS (first 2 inputs) after intensity cohering  = %2.2f (sigma=%2.2f)\n", rms_after, cross_time_sigma) ;
  }

  for (input = 0 ; input < ninputs ; input++)
  {
    if (longinput)
      sprintf(fname, "%s/%s.long.%s/mri/%s", sdir, subjects[input], base_name, out_fname) ;
    else
      sprintf(fname, "%s/%s/longtp/%s/%s", sdir, base_name, subjects[input], out_fname) ;
    
      printf("writing normalized volume to %s\n", fname) ;
    if (MRIwriteFrame(mri_in, fname, input)  != NO_ERROR)
      ErrorExit(ERROR_BADFILE, "%s: could not write normalized volume to %s",Progname, fname);
  }

  if (ctrl_point_fname)
  {
    printf("writing control points to %s\n", ctrl_point_fname) ;
    MRIwrite(mri_ctrl, ctrl_point_fname) ;
    MRIfree(&mri_ctrl) ;
  }
  MRIfree(&mri_in) ;

  if (gca)
    GCAfree(&gca) ;
  printf("done.\n") ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("normalization took %d minutes and %d seconds.\n",
         minutes, seconds) ;
  if (diag_fp)
    fclose(diag_fp) ;
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
    usage_exit(0) ;
  else if (!strcmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!strcmp(option, "ASEG"))
  {
    aseg_fname = argv[2] ;
    nargs = 1 ;
    printf("using segmentation volume %s to generate control points...\n",
           aseg_fname) ;
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
  else if (!strcmp(option, "NOCEREBELLUM"))
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
  else if (!strcmp(option, "LONGINPUT"))
  {
    longinput = 1;
    printf("reading inputs from longitudinal directories\n") ; 
  }
  else if (!strcmp(option, "SIGMA"))
  {
    bias_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("smoothing bias field with sigma = %2.1f\n", bias_sigma) ;
  }
  else if (!stricmp(option, "cross_time_sigma"))
  {
    cross_time_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("smoothing temporal bias field with sigma = %2.1f\n", cross_time_sigma) ;
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
  else if (!strcmp(option, "RENORM"))
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
find_control_points(GCA *gca, GCA_SAMPLE *gcas_total,
                    int total_samples, int *pnorm_samples, int nregions, int label,
                    MRI *mri_in, TRANSFORM *transform, double min_prior, double ctl_point_pct)
{
  int        i, j, *ordered_indices, nsamples,
    xmin, ymin, zmin, xmax, ymax, zmax, xv,yv,zv,
    x, y, z, xi, yi, zi, region_samples,
    used_in_region, prior_wsize=5, image_wsize=3, histo_peak, n,
    nbins ;
  GCA_SAMPLE *gcas, *gcas_region, *gcas_norm ;
  double     means[MAX_GCA_INPUTS], vars[MAX_GCA_INPUTS], val, nsigma ;
  HISTOGRAM  *histo, *hsmooth ;
  GC1D       *gc ;
  float      fmin, fmax ;
  MRI        *mri_T1 = NULL ;


  if (label == Gdiag_no)
    DiagBreak() ;

  MRIvalRange(mri_in, &fmin, &fmax) ;
  nbins = (int)(fmax-fmin+1);
  histo = HISTOalloc(nbins) ; hsmooth = HISTOalloc(nbins) ;
  for (nsamples = i = 0 ; i < total_samples ; i++)
  {
    if (gcas_total[i].label != label)
      continue ;
    nsamples++ ;
  }

  *pnorm_samples = 0 ;
  printf("found %d control points for structure...\n", nsamples) ;
  if (nsamples == 0)
  {
    DiagBreak() ;
    return(NO_ERROR) ;
  }
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

  gcas_bounding_box(gcas, nsamples, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax, label) ;
  printf("bounding box (%d, %d, %d) --> (%d, %d, %d)\n",
         xmin, ymin, zmin, xmax, ymax, zmax) ;
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
            if (xi != x || yi != y || zi != z
                || gcas_getPrior(gcas[i]) < min_prior)
              continue ;

            if (xv == Gx && yv == Gy && zv == Gz)
              DiagBreak() ;
            if (sqrt(SQR(xv-Gx)+SQR(yv-Gy)+SQR(zv-Gz)) < 2)
              DiagBreak() ;
            if (min_region_prior(gca, gcas[i].xp, gcas[i].yp, gcas[i].zp,prior_wsize, label) < min_prior)
              continue ;
            if (uniform_region(gca, mri_in, transform,
                               xv, yv, zv,
                               image_wsize, &gcas[i], nsigma) == 0)
              continue ;
            memmove(&gcas_region[region_samples],
                    &gcas[i],
                    sizeof(GCA_SAMPLE)) ;
            region_samples++ ;
            if (gcas[i].x == Gx &&
                gcas[i].y == Gy &&
                gcas[i].z == Gz)
              DiagBreak() ;
          }
          nsigma *= 1.1 ;
        } while (region_samples < 8 && nsigma < 3) ;

        if (region_samples < 8)/* can't reliably estimate statistics */
          continue ;
        if (DIAG_VERBOSE_ON)
          printf("\t%d total samples found in region (%d, %d, %d)\n",
                 region_samples,x, y,z) ;
        /* compute mean and variance of label within this region */
        for (n = 0 ; n < mri_in->nframes ; n++)
        {
          HISTOclear(histo, histo) ;
          histo->bin_size = 1 ;
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
          }

          HISTOsmooth(histo, hsmooth, 2) ;
          histo_peak =
            HISTOfindHighestPeakInRegion(hsmooth, 1, hsmooth->nbins) ;
          if (histo_peak < 0)   /* couldn't find a valid peak? */
            break ;

          for (means[n] = vars[n] = 0.0, i = 0 ;
               i < region_samples ;
               i++)
          {
            if (gcas_region[i].label < 0)
              continue ;
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
          vars[n] = vars[n] / (double)region_samples - means[n]*means[n] ;

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
            gcas_region[i].means[r] = means[r] ;
          /*          gcas_region[i].var = var ;*/
        }

        GCAcomputeLogSampleProbability
          (gca, gcas_region, mri_in, transform, region_samples, DEFAULT_CLAMP);
        GCArankSamples
          (gca, gcas_region, region_samples, ordered_indices) ;
        GCAremoveOutlyingSamples
          (gca, gcas_region, mri_in, transform, region_samples, 2.0) ;
        for (used_in_region = i = 0 ; i < region_samples ; i++)
        {
          j = ordered_indices[i] ;
          if (gcas_region[j].label != label)  /* it was an outlier */
            continue ;
          memmove
            (&gcas_norm[*pnorm_samples],
             &gcas_region[j],
             sizeof(GCA_SAMPLE)) ;
          (*pnorm_samples)++ ; used_in_region++ ;
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
          (*pnorm_samples)++ ; used_in_region++ ;
        }
        if (DIAG_VERBOSE_ON)
          printf("\t%d samples used in region\n", used_in_region) ;
      }
    }
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
  HISTOfree(&histo) ; HISTOfree(&hsmooth) ;
  free(gcas_region) ;
  free(gcas) ;
  if (mri_T1)
    MRIfree(&mri_T1) ;
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
    memmove(&gcas[i], &gcas1[i], sizeof(GCA_SAMPLE)) ;
  for (i = 0 ; i < n2 ; i++)
    memmove(&gcas[i+n1], &gcas2[i], sizeof(GCA_SAMPLE)) ;

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
      xmin = gcas[i].x ;
    if (gcas[i].y < ymin)
      ymin = gcas[i].y ;
    if (gcas[i].z < zmin)
      zmin = gcas[i].z ;

    if (gcas[i].x > xmax)
      xmax = gcas[i].x ;
    if (gcas[i].y > ymax)
      ymax = gcas[i].y ;
    if (gcas[i].z > zmax)
      zmax = gcas[i].z ;
  }

  *pxmin = xmin ; *pymin = ymin ; *pzmin = zmin ;
  *pxmax = xmax ; *pymax = ymax ; *pzmax = zmax ;
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
      continue ;
    for (yi = -whalf ; yi <= whalf ; yi++)
    {
      yk = yp+yi ;
      if (yk < 0 || yk >= gca->prior_height)
        continue ;
      for (zi = -whalf ; zi <= whalf ; zi++)
      {
        zk = zp+zi ;
        if (zk < 0 || zk >= gca->prior_depth)
          continue ;
        gcap = &gca->priors[xk][yk][zk] ;
        prior = getPrior(gcap, label) ;
        if (prior < min_prior)
          min_prior = prior ;
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
  int   xk, yk, zk, whalf, xi, yi, zi, n, frame ;
  double val0, val, sigma, min_val,max_val, thresh ;
  MATRIX *m ;
  GC1D   *gc ;

  gc = GCAfindSourceGC(gca, mri, transform, x, y, z, gcas->label) ;
  if (!gc)
    return(0) ;
  m = load_covariance_matrix(gc, NULL, gca->ninputs) ;

  whalf = (wsize-1)/2 ;
  for (n = 0 ; n < gca->ninputs ; n++)
  {
    sigma = sqrt(*MATRIX_RELT(m, n+1, n+1)) ;
    MRIsampleVolumeFrame(mri, (double)x, (double)y, (double)z, n, &val0) ;
    if (sigma < 0.05*val0)   /* don't let it be too small */
      sigma = 0.05*val0 ;
    if (sigma > 0.1*val0)    /* don't let it be too big */
      sigma = 0.1*val0 ;
    min_val = max_val = val0 ;
    thresh = nsigma*sigma ;

    for (frame = 0 ; frame < mri->nframes ; frame++)
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
            (mri, (double)xi, (double)yi, (double)zi, frame, &val) ;
          if (val < min_val)
            min_val = val ;
          if (val > max_val)
            max_val = val ;
          if (fabs(val-val0) > thresh ||
              fabs(max_val-min_val) > thresh)
            return(0) ;
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
  double val,  mean_ratio ;

  if (nsamples == 0)
    return(NO_ERROR) ;

  for (num = n = 0 ; n < mri_in->nframes ; n++)
  {
    int niter = 0 ;
    MRIvalRangeFrame(mri_in, &fmin, &fmax, n) ;
    h = HISTOalloc(nint(fmax-fmin)+1) ;
    h->bin_size = (fmax-fmin)/(float)h->nbins ;
    for (i = 0 ; i < h->nbins ; i++)
      h->bins[i] = (i+1)*h->bin_size+fmin ;

    for (i = 0 ; i < nsamples ; i++)
    {
      xv = gcas[i].x ; yv = gcas[i].y ; zv = gcas[i].z ;
      if (xv == Gx && yv == Gy && zv == Gz)
        DiagBreak() ;
      MRIsampleVolumeFrame(mri_in, gcas[i].x,gcas[i].y,gcas[i].z, n, &val) ;
      if (FZERO(val))
        DiagBreak() ;
      h->counts[nint(val-fmin)]++ ;
    }

    /* check to see  if peak is unlikely */
    hsmooth = HISTOsmooth(h, NULL, 2) ;
    do
    {
      if (gca->ninputs == 1) /*  find  brightest peak as
                                 for  n=1 it is T1  weighted  */
        peak = HISTOfindLastPeak(hsmooth, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
      else
        peak = HISTOfindHighestPeakInRegion(hsmooth, 0, h->nbins-1) ;

      end = HISTOfindEndOfPeak(hsmooth, peak, 0.01) ;
      start = HISTOfindStartOfPeak(hsmooth, peak, 0.01) ;
      if (gca->ninputs == 1) 
      {
	int opeak ;
        opeak = HISTOfindFirstPeak(hsmooth, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
	if (hsmooth->bins[opeak] < start)
	{
	  start = HISTOfindStartOfPeak(hsmooth, opeak, 0.01) ;
	}
      }
      for (mean_ratio = 0.0, i = 0 ; i < nsamples ; i++)
      {
        mean_ratio += hsmooth->bins[peak] / gcas[i].means[0];
      }
      mean_ratio /= (double)nsamples ;
      HISTOclearBins
        (hsmooth, hsmooth, hsmooth->bins[start], hsmooth->bins[end])  ;
      if (niter++ > 5)
        break ;
      if (niter > 1)
        DiagBreak() ;
    } while  (mean_ratio  < 0.5 || mean_ratio > 2.0) ;

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
        num++ ; gcas[i].label = 0 ;
      }
    }
    HISTOfree(&h) ; HISTOfree(&hsmooth) ;
  }

  printf("%d of %d (%2.1f%%) samples deleted\n",
         num, nsamples, 100.0f*(float)num/(float)nsamples) ;
  return(NO_ERROR) ;
}

static void
usage_exit(int code)
{
  //  outputHelp(Progname);
#ifdef GREGT
  printf("usage: %s [<options>] <in volume> <atlas> <transform> "
         "<normalized volume>\n\n",
         Progname) ;
  printf("\t<atlas>                      path/to/some.gca file "
         "(or 'noatlas') \n");
  printf("\t<transform>                  ex. transforms/talairach.lta "
         "(or 'noxform') \n");
  printf("\noptions:\n");
  printf("\t-aseg <filename>             aseg file, to help normalization\n");
  printf("\t-sigma <bias sigma>          smoothing sigma for bias field if control points specified (def=4)\n");
  printf("\t-fsamples <filename>         write control points to filename\n");
  printf("\t-nsamples <filename>         write transformed "
         "normalization control points to filename\n");
  printf("\t-mask <mri_vol>              use mri_vol to mask input\n");
  printf("\t-f <filename>                define control points "
         "from filename\n");
  printf("\t-fonly <filename>            only use control points "
         "from filename\n");
  printf("\t-longinput                   load aseg and nu from <tp>.long.<base> dirs, "
         "instead of <base>/longtp subdirs.\n");
  printf("\t-diag <filename>             write to log file\n");
  printf("\t-debug_voxel <x> <y> <z>     debug voxel\n");
  printf("\t-debug_node <x> <y> <z>      debug node\n");
  printf("\t-tr <float n>                set TR in msec\n");
  printf("\t-te <float n>                set TE in msec\n");
  printf("\t-alpha <float n>             set alpha in radians\n");
  printf("\t-example <mri_vol> <segmentation> use T1 (mri_vol) "
         "and segmentation as example\n");
  printf("\t-novar                       do not use variance estimates\n");
  printf("\t-renorm <mri_vol>            renormalize using "
         "predicted intensity values in mri_vol\n");
  printf("\t-flash                       use FLASH forward model "
         "to predict intensity values\n");
  printf("\t-prior <float t>             use prior threshold t "
         "(default=.6)\n");
  printf("\t-w                           write normalized volume "
         "each nregion iteration to norm(n).mgh(see -n)\n");
  printf("\t-n <int n>                   use n regions/struct "
         "for normalization\n");
  printf("\t-v <int n>                   does nothing as far "
         "as i can tell, but an option\n");
  printf("\t-p <float p>                 use top p percent(default=.25) "
         "white matter points as control points\n");
#endif
  exit(code) ;
}
static int
copy_ctrl_points_to_volume(GCA_SAMPLE *gcas, int nsamples, MRI *mri_ctrl, int frame)
{
  int   i, xv, yv, zv, nregions ;

  nregions = mri_ctrl->nframes/2 ;
  for (i = 0 ; i < nsamples ; i++)
  {
    xv = gcas[i].x ; yv = gcas[i].y ; zv = gcas[i].z ;
    if (xv == Ggca_x && yv == Ggca_y && zv == Ggca_z)
      DiagBreak() ;
    MRIsetVoxVal(mri_ctrl, xv, yv, zv, frame, gcas[i].label) ;
    MRIsetVoxVal(mri_ctrl, xv, yv, zv, frame+nregions, gcas[i].means[0]) ;
  }

  return(NO_ERROR) ;
}

#if 0
static int
normalize_timepoints_with_samples(MRI *mri, GCA_SAMPLE *gcas, int nsamples, int nsoap)
{
  int   frame, i, x, y, z ;
  double target, val ;
  MRI    *mri_ctrl, *mri_bias, *mri_target, *mri_frame ;

  mri_ctrl = MRIcloneDifferentType(mri, MRI_UCHAR) ;
  mri_bias = MRIcloneDifferentType(mri, MRI_FLOAT) ;
  mri_target = MRIcloneDifferentType(mri, MRI_FLOAT) ;

  for (i = 0 ; i < nsamples ; i++)
  {
    if (i == Gdiag_no)
      DiagBreak() ;
    x = nint(gcas[i].x)  ; y = nint(gcas[i].y) ; z = nint(gcas[i].z) ;
    MRIsetVoxVal(mri_ctrl, x, y, z, 0, CONTROL_MARKED) ;
    for (target = 0.0, frame = 0 ; frame < mri->nframes ; frame++)
      target += MRIgetVoxVal(mri, x, y, z, frame) ;
    target /= mri->nframes ;
    MRIsetVoxVal(mri_target, x, y, z, 0, target) ;
  }

  // build a bias correction for each time point (which each has its own frame)
  for (frame = 0 ; frame < mri->nframes ; frame++)
  {
    MRIclear(mri_bias) ; 
    for (i = 0 ; i < nsamples ; i++)
    {
      if (i == Gdiag_no)
        DiagBreak() ;
      x = nint(gcas[i].x)  ; y = nint(gcas[i].y) ; z = nint(gcas[i].z) ;
      target = MRIgetVoxVal(mri_target, x, y, z, 0) ;
      val = MRIgetVoxVal(mri, x, y, z, frame) ;
      if (FZERO(val))
        val = 1.0 ;
      MRIsetVoxVal(mri_bias, x, y, z, 0, target/val) ;
    }
    MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias) ;
    MRIsoapBubble(mri_bias, mri_ctrl, mri_bias, nsoap, 1) ;
    mri_frame = MRIcopyFrame(mri, NULL, frame, 0) ;
    MRImultiply(mri_frame, mri_bias, mri_frame) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      char fname[STRLEN] ;
      sprintf(fname, "frame%d.mgz", frame) ;
      MRIwrite(mri_frame, fname) ;
      sprintf(fname, "bias%d.mgz", frame) ;
      MRIwrite(mri_bias, fname) ;
      sprintf(fname, "target%d.mgz", frame) ;
      MRIwrite(mri_target, fname) ;
    }
    MRIcopyFrame(mri_frame, mri, 0, frame) ;
  }
  MRIfree(&mri_bias) ; MRIfree(&mri_target) ; MRIfree(&mri_ctrl) ;
  return(NO_ERROR) ;
}
#endif
static int
normalize_timepoints(MRI *mri, double thresh, double cross_time_sigma)
{
  int   frame, x, y, z, skip, nvox ;
  double target, val ;
  MRI    *mri_ctrl, *mri_bias, *mri_target, *mri_frame, *mri_kernel ;

  mri_ctrl = MRIcloneDifferentType(mri, MRI_UCHAR) ;
  mri_bias = MRIcloneDifferentType(mri, MRI_FLOAT) ;
  mri_target = MRIcloneDifferentType(mri, MRI_FLOAT) ;
  mri_kernel = MRIgaussian1d(cross_time_sigma, -1) ;

  for (nvox = x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        for (target = 0.0, frame = 0 ; frame < mri->nframes ; frame++)
          target += MRIgetVoxVal(mri, x, y, z, frame) ;
        target /= mri->nframes ;
        if (FZERO(target))
          continue ;  // both vals  0
        skip = 0 ;
        for (frame = 0 ; frame < mri->nframes ; frame++)
        {
          val = MRIgetVoxVal(mri, x, y, z, frame) ;
          if (fabs(val-target) > thresh)
          {
            skip = 1 ;
            break ;
          }
        }
        if (skip)
          continue ;
        nvox++ ;
        MRIsetVoxVal(mri_ctrl, x, y, z, 0, CONTROL_MARKED) ;
        MRIsetVoxVal(mri_target, x, y, z, 0, target) ;
      }

  printf("%d voxels found to base intensity correction on\n", nvox) ;

  // build a bias correction for each time point (which each has its own frame)
  for (frame = 0 ; frame < mri->nframes ; frame++)
  {
    MRIclear(mri_bias) ; 
    for (x = 0 ; x < mri->width ; x++)
      for (y = 0 ; y < mri->height ; y++)
        for (z = 0 ; z < mri->depth ; z++)
        {
          target = MRIgetVoxVal(mri_target, x, y, z, 0) ;
          val = MRIgetVoxVal(mri, x, y, z, frame) ;
          if (FZERO(val))
            val = 1.0 ;
          MRIsetVoxVal(mri_bias, x, y, z, 0, target/val) ;
        }
    MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias) ;
    MRIconvolveGaussian(mri_bias, mri_bias, mri_kernel) ;
    //    MRIsoapBubble(mri_bias, mri_ctrl, mri_bias, nsoap,1) ;
    mri_frame = MRIcopyFrame(mri, NULL, frame, 0) ;
    MRImultiply(mri_frame, mri_bias, mri_frame) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      char fname[STRLEN] ;
      sprintf(fname, "frame%d.mgz", frame) ;
      MRIwrite(mri_frame, fname) ;
      sprintf(fname, "bias%d.mgz", frame) ;
      MRIwrite(mri_bias, fname) ;
      sprintf(fname, "target%d.mgz", frame) ;
      MRIwrite(mri_target, fname) ;
    }
    MRIcopyFrame(mri_frame, mri, 0, frame) ;
  }
  MRIfree(&mri_bias) ; MRIfree(&mri_kernel) ; MRIfree(&mri_target) ; MRIfree(&mri_ctrl) ;
  return(NO_ERROR) ;
}



static MRI *
normalize_timepoints_with_parzen_window(MRI *mri_in, MRI *mri_out, double cross_time_sigma)
{
  int   frame1, frame2, x, y, z ;
  double val0, val, total, g, total_norm ;
  double s;
  // if sigma is zero, skip smoothing
  // copying and returning here avoids numerical problems below (division by zero)
  if (fabs(cross_time_sigma) <= 0.000001)
  {
     mri_out = MRIcopy(mri_in, mri_out);
     return mri_out;
  }
  s = -0.5 / SQR(cross_time_sigma);
  if (mri_out == NULL)
    mri_out = MRIclone(mri_in, NULL) ;
  for (x = 0 ; x < mri_in->width ; x++)
    for (y = 0 ; y < mri_in->height ; y++)
      for (z = 0 ; z < mri_in->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        for (frame1 = 0 ; frame1 < mri_in->nframes ; frame1++)
        {
          val0 = MRIgetVoxVal(mri_in, x, y, z, frame1) ;
          for (total = total_norm = 0.0, frame2 = 0 ; frame2 < mri_in->nframes ; frame2++)
          {
            val = MRIgetVoxVal(mri_in, x, y, z, frame2) ;
            g = exp( SQR(val-val0) * s) ;
            total += g*val ; 
            total_norm += g ;
          }
          total /= total_norm ;
          MRIsetVoxVal(mri_out, x, y, z, frame1, total) ;
        }
      }

  return(mri_out) ;
}


#define HBINS 1000

static MRI *
scale_all_images(MRI *mri_in, MRI *mri_out)
{
  HISTOGRAM *h ;
  int       t, b ;
  MRI       *mri_ratio, *mri_f0, *mri_f ;
  float     scale ;

  if (mri_out == NULL)
    mri_out = MRIclone(mri_in, NULL) ;

  
  mri_f0 = MRIcopyFrame(mri_in, NULL, 0, 0) ;
  MRIthresholdRangeInto(mri_f0, mri_f0, 50, 120) ;

  for (t = 1 ; t < mri_in->nframes ; t++)
  {
    mri_f = MRIcopyFrame(mri_in, NULL, t, 0) ;
    mri_ratio = MRIdivide(mri_f0, mri_f, NULL) ;
    h = MRIhistogram(mri_ratio, HBINS) ;
    b = HISTOfindHighestPeakInRegion(h, 1, h->nbins) ;   // ignore zero bin
    scale = h->bins[b] ;
    printf("scaling image %d by %2.3f\n", t+1, scale) ;
    MRIscalarMul(mri_f, mri_f, scale) ;
    MRIcopyFrame(mri_f, mri_out, 0, t) ;

    MRIfree(&mri_f) ; MRIfree(&mri_ratio) ; HISTOfree(&h) ;
  }

  MRIfree(&mri_f0) ;
  return(mri_out) ;
}

static int
discard_control_points_with_different_labels(GCA_SAMPLE *gcas, int nsamples, MRI *mri_aseg)
{
  int i, x, y, z, f, delete_this_sample, deleted = 0, label, label2 ;

  for (i = 0 ; i < nsamples ; i++)
  {
    x = nint(gcas[i].x) ; y = nint(gcas[i].y) ; z = nint(gcas[i].z);
    label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
    delete_this_sample = IS_HYPO(label) ;
    for (f = 1 ; f < mri_aseg->nframes ; f++)
    {
      label2 = MRIgetVoxVal(mri_aseg, x, y, z, f) ;
      if (label2 != label || IS_HYPO(label2))
      {
	delete_this_sample = 1 ;
	break ;
      }
    }
    if (delete_this_sample)
    {
      gcas[i].label = 0 ; deleted++ ;
    }
  }
  printf("%d control points deleted due to different labels across time\n", deleted) ;
  return(deleted) ;
}

