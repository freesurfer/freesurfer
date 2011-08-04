/**
 * @file  mri_segment.c
 * @brief segments white matter from a brain volume
 *
 * "Cortical Surface-Based Analysis I: Segmentation and Surface
 * Reconstruction", Dale, A.M., Fischl, B., Sereno, M.I.
 * (1999) NeuroImage 9(2):179-194
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/08/04 19:39:41 $
 *    $Revision: 1.41 $
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

const char *MRI_SEGMENT_VERSION = "$Revision: 1.41 $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "utils.h"
#include "proto.h"
#include "classify.h"
#include "mrisegment.h"
#include "mri.h"
#include "tags.h"
#include "mrinorm.h"
#include "timer.h"
#include "version.h"

#define BRIGHT_LABEL         130
#define BRIGHT_BORDER_LABEL  100

static int extract = 0 ;
static int  verbose = 0 ;
static int wsize = 11 ;
static float pct = 0.8 ;
static float pslope = 1.0f ;
static float nslope = 1.0f ;
static float wm_low = 90 ;
static float wm_hi = 125 ;
static float gray_hi = 100 ;
static float gray_low = 30 ;
static int gray_low_set = 0 ;
static int niter = 1 ;
static int gray_hi_set = 0 ;
static int wm_hi_set = 0 ;
static int wm_low_set = 0 ;

static int scan_type = MRI_UNKNOWN ;
static int thickness = 4 ;
static int thicken = 1 ;
static int nsegments = 20 ;
static int fill_bg = 0 ;
static int fill_ventricles = 0 ;

static int keep_edits = 0 ;

static int auto_detect_stats =  1 ;
static int log_stats = 1 ;
static void  usage_exit(int code) ;

#define BLUR_SIGMA 0.25f
static float blur_sigma = BLUR_SIGMA ;

const char *Progname ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *MRIremoveWrongDirection(MRI *mri_src, MRI *mri_dst, int wsize,
                             float low_thresh, float hi_thresh,
                             MRI *mri_labels) ;
MRI *MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm) ;
MRI *MRIfilterMorphology(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIfillBasalGanglia(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIfillVentricles(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIremove1dStructures(MRI *mri_src, MRI *mri_dst, int max_iter,
                           int thresh, MRI *mri_labels) ;

static MRI *MRIrecoverBrightWhite(MRI *mri_T1, MRI *mri_src, MRI *mri_dst,
                                  float wm_low, float wm_hi, float slack,
                                  float pct_thresh) ;
static int is_diagonal(MRI *mri, int x, int y, int z) ;
#if 0
MRI *MRIremoveFilledBrightStuff(MRI *mri_src, MRI *mri_dst, int filled_label,
                                float thresh) ;
static int MRIcheckRemovals(MRI *mri_T1, MRI *mri_dst,
                            MRI *mri_labels, int wsize) ;
#endif

int
main(int argc, char *argv[])
{
  MRI     *mri_src, *mri_dst, *mri_tmp, *mri_labeled, *mri_labels;
  char    *input_file_name, *output_file_name ;
  int     nargs, i, msec ;
  struct timeb  then ;
  float   white_mean, white_sigma, gray_mean, gray_sigma ;

  char cmdline[CMD_LINE_LEN] ;

  TAGmakeCommandLineString(argc, argv, cmdline) ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mri_segment.c,v 1.41 2011/08/04 19:39:41 fischl Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
  {
    usage_exit(1);
  }

  TimerStart(&then) ;
  input_file_name = argv[1] ;
  output_file_name = argv[2] ;

  mri_src = MRIread(input_file_name) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read source volume from %s",
              Progname, input_file_name) ;
  MRIaddCommandLine(mri_src, cmdline) ;
  if (mri_src->type != MRI_UCHAR)
  {
    MRI *mri_tmp ;
    printf("changing input type from %d to UCHAR\n", mri_src->type) ;
    mri_tmp = MRIchangeType(mri_src, MRI_UCHAR, 0, 1000, 1) ;
    MRIfree(&mri_src) ;
    mri_src = mri_tmp ;
  }

  if (thicken > 1)
  {
    mri_dst = MRIcopy(mri_src, NULL) ;
    /*    MRIfilterMorphology(mri_dst, mri_dst) ;*/
    fprintf(stderr, "removing 1-dimensional structures...\n") ;
    MRIremove1dStructures(mri_dst, mri_dst, 10000, 2, NULL) ;
#if 0
    MRIcheckRemovals(mri_src, mri_dst, mri_labels, 5) ;
    fprintf(stderr, "thickening thin strands....\n") ;
    MRIthickenThinWMStrands(mri_src, mri_dst, mri_dst, thickness, nsegments,
                            wm_hi) ;
#endif
    MRIwrite(mri_dst, output_file_name) ;
    exit(0) ;
  }

  mri_labels = MRIclone(mri_src, NULL) ;
  if (auto_detect_stats && !wm_low_set) /* widen range to allow
                                           for more variability */
  {
    wm_low -= 10 ;
  }
  fprintf(stderr, "doing initial intensity segmentation...\n") ;
  mri_tmp = MRIintensitySegmentation(mri_src, NULL, wm_low, wm_hi, gray_hi);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_tmp, "tmp1.mgz") ;
  }
  fprintf(stderr, "using local statistics to label ambiguous voxels...\n") ;
  MRIhistoSegment(mri_src, mri_tmp, wm_low, wm_hi, gray_hi, wsize, 3.0f) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_tmp, "tmp2.mgz") ;
  }

  if (auto_detect_stats)
  {

    fprintf(stderr, "computing class statistics for intensity windows...\n") ;
    MRIcomputeClassStatistics(mri_src, mri_tmp, gray_low, WHITE_MATTER_MEAN,
                              &white_mean, &white_sigma, &gray_mean,
                              &gray_sigma) ;
    if (!finite(white_mean) || !finite(white_sigma) ||
        !finite(gray_mean) || !finite(gray_sigma))
      ErrorExit
      (ERROR_BADPARM,
       "%s: class statistics not finite - check input volume!",
       Progname);

    if (!wm_low_set)
    {
      if (FZERO(gray_sigma))
      {
        wm_low = (white_mean+gray_mean) / 2 ;
      }
      else
      {
        wm_low = gray_mean + gray_sigma ;
      }
    }

    if (!gray_hi_set)
    {
      gray_hi = gray_mean + 2*gray_sigma ;
#if 1
      if (gray_hi >= white_mean)
      {
        gray_hi = white_mean-1 ;
      }
#endif
    }
    fprintf(stderr, "setting bottom of white matter range to %2.1f\n",wm_low);
    fprintf(stderr, "setting top of gray matter range to %2.1f\n", gray_hi) ;

    if (log_stats)
    {
      FILE *fp ;

      fp = fopen("segment.dat", "w") ;
      if (fp)
      {
        fprintf(fp, "WM: %2.1f +- %2.1f\n",white_mean, white_sigma) ;
        fprintf(fp, "GM: %2.1f +- %2.1f\n",gray_mean, gray_sigma) ;
        fprintf(fp, "setting bottom of white matter range to %2.1f\n",wm_low);
        fprintf(fp, "setting top of gray matter range to %2.1f\n", gray_hi) ;
        fclose(fp) ;
      }
    }

    fprintf(stderr, "doing initial intensity segmentation...\n") ;
    mri_tmp = MRIintensitySegmentation(mri_src, NULL, wm_low, wm_hi, gray_hi);

    fprintf(stderr, "using local statistics to label ambiguous voxels...\n") ;
    MRIhistoSegment(mri_src, mri_tmp, wm_low, wm_hi, gray_hi, wsize, 3.0f) ;
  }
  else
  {
    /* just some not-too-dopey defaults - won't really be used */
    white_mean =  110 ;
    white_sigma = 5.0 ;
    gray_mean = 65 ;
    gray_sigma = 12 ;
  }

  fprintf(stderr,
          "using local geometry to label remaining ambiguous voxels...\n") ;
  mri_labeled = MRIcpolvMedianCurveSegment(mri_src, mri_tmp, NULL, 5, 3,
                gray_hi, wm_low);
  fprintf(stderr,
          "\nreclassifying voxels using Gaussian border classifier...\n") ;

  /*
    now use the gray and white matter border voxels to build a Gaussian
    classifier at each point in space and reclassify all voxels in the
    range [wm_low-5,gray_hi].
    */
  for (i = 0 ; i < niter ; i++)
  {
    MRIreclassify(mri_src, mri_labeled, mri_labeled, wm_low-5,gray_hi,wsize);
  }
  MRIfree(&mri_tmp) ;

  mri_dst = MRImaskLabels(mri_src, mri_labeled, NULL) ;
  MRIfree(&mri_labeled) ;
  MRIrecoverBrightWhite(mri_src, mri_dst,mri_dst,wm_low,wm_hi,white_sigma,.33);
  fprintf(stderr,
          "\nremoving voxels with positive offset direction...\n") ;

#if 0
  MRIremoveWrongDirection(mri_dst, mri_dst, 3, wm_low-5, gray_hi, mri_labels) ;
#else
  MRIremoveWrongDirection(mri_dst, mri_dst, 3, wm_low-5, gray_hi, NULL) ;
#endif

  if (thicken)
  {
    /*    MRIfilterMorphology(mri_dst, mri_dst) ;*/
    fprintf(stderr, "removing 1-dimensional structures...\n") ;
    MRIremove1dStructures(mri_dst, mri_dst, 10000, 2, mri_labels) ;
#if 0
    MRIcheckRemovals(mri_src, mri_dst, mri_labels, 5) ;
#endif
    fprintf(stderr, "thickening thin strands....\n") ;
    MRIthickenThinWMStrands(mri_src, mri_dst, mri_dst, thickness, nsegments,
                            wm_hi) ;
  }

  mri_tmp = MRIfindBrightNonWM(mri_src, mri_dst) ;
  MRIbinarize(mri_tmp, mri_tmp, WM_MIN_VAL, 255, 0) ;
  MRImaskLabels(mri_dst, mri_tmp, mri_dst) ;
  MRIfilterMorphology(mri_dst, mri_dst) ;

  if (fill_bg)
  {
    fprintf(stderr, "filling basal ganglia....\n") ;
    MRIfillBasalGanglia(mri_src, mri_dst) ;
  }
  if (fill_ventricles)
  {
    fprintf(stderr, "filling ventricles....\n") ;
    MRIfillVentricles(mri_dst, mri_dst) ;
  }


  MRIfree(&mri_src) ;
  msec = TimerStop(&then) ;
  fprintf(stderr, "white matter segmentation took %2.1f minutes\n",
          (float)msec/(1000.0f*60.0f));
  fprintf(stderr, "writing output to %s...\n", output_file_name) ;
  if (keep_edits)
  {
    MRI *mri_old ;

    mri_old = MRIread(output_file_name) ;
    if (!mri_old)
    {
      ErrorPrintf
      (ERROR_NOFILE, "%s: could not read file %s to preserve edits",
       Progname, output_file_name) ;
      exit(1);
    }
    else
    {
      MRIcopyLabel(mri_old, mri_dst, WM_EDITED_ON_VAL) ;
      MRIcopyLabel(mri_old, mri_dst, WM_EDITED_OFF_VAL) ;
      MRIfree(&mri_old) ;
    }
  }
  MRIwrite(mri_dst, output_file_name) ;

  MRIfree(&mri_dst) ;

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
  if (!strcasecmp(option, "-version"))
  {
    fprintf(stderr, "Version: %s\n", MRI_SEGMENT_VERSION);
    exit(0);
  }
  else if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    usage_exit(0);
  }
  else if (!stricmp(option, "MGH_MPRAGE") || !stricmp(option, "MPRAGE"))
  {
    scan_type = MRI_MGH_MPRAGE;
    printf("assuming input volume is MGH (Van der Kouwe) MP-RAGE\n") ;
    gray_hi = 99 ;
    wm_low = 89 ;
  }
  else if (!stricmp(option, "WASHU_MPRAGE"))
  {
    scan_type = MRI_WASHU_MPRAGE;
    printf("assuming input volume is WashU MP-RAGE (dark GM)\n") ;
    gray_hi = 85 ;
    wm_low = 80 ;
  }
  else if (!stricmp(option, "slope"))
  {
    nslope = pslope = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using curvature slope = %2.2f\n", pslope) ;
  }
  else if (!stricmp(option, "pslope"))
  {
    pslope = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using curvature pslope = %2.2f\n", pslope) ;
  }
  else if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "auto"))
  {
    auto_detect_stats = !auto_detect_stats ;
    fprintf(stderr, "%sautomatically detecting class statistics...\n",
            auto_detect_stats ? "" : "not ") ;
  }
  else if (!stricmp(option, "noauto"))
  {
    auto_detect_stats = 0 ;
    fprintf(stderr, "%sautomatically detecting class statistics...\n",
            auto_detect_stats ? "" : "not ") ;
  }
  else if (!stricmp(option, "log"))
  {
    log_stats = 1 ;
    fprintf(stderr, "logging class statistics and thresholds...\n") ;
  }
  else if (!stricmp(option, "keep"))
  {
    keep_edits = 1 ;
    fprintf(stderr, "preserving editing changes in output volume...\n");
  }
  else if (!stricmp(option, "nslope"))
  {
    nslope = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using curvature nslope = %2.2f\n", nslope) ;
  }
  else if (!stricmp(option, "ghi") || !stricmp(option, "gray_hi"))
  {
    gray_hi = atof(argv[2]) ;
    gray_hi_set = 1 ;
    nargs = 1 ;
    fprintf(stderr, "using gray hilim = %2.1f\n", gray_hi) ;
  }
  else if (!stricmp(option, "glo") || !stricmp(option, "gray_low"))
  {
    gray_low = atof(argv[2]) ;
    gray_low_set = 1 ;
    nargs = 1 ;
    fprintf(stderr, "using gray low limit = %2.1f\n", gray_low) ;
  }
  else if (!stricmp(option, "wlo") || !stricmp(option, "wm_low"))
  {
    wm_low = atof(argv[2]) ;
    wm_low_set = 1 ;
    nargs = 1 ;
    fprintf(stderr, "using white lolim = %2.1f\n", wm_low) ;
  }
  else if (!stricmp(option, "whi") || !stricmp(option, "wm_hi"))
  {
    wm_hi = atof(argv[2]) ;
    wm_hi_set = 1 ;
    nargs = 1 ;
    fprintf(stderr, "using white hilim = %2.1f\n", wm_hi) ;
  }
  else if (!stricmp(option, "nseg"))
  {
    nsegments = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"thickening the %d largest thin strands\n",
            nsegments) ;
  }
  else if (!stricmp(option, "thicken"))
  {
    thicken = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"%sthickening thin strands\n", thicken ? "" : "not ") ;
  }
  else if (!stricmp(option, "thickenonly"))
  {
    thicken = 2 ;
    fprintf(stderr,"%sthickening thin strands\n", thicken ? "" : "not ") ;
  }
  else if (!stricmp(option, "fillbg") || !stricmp(option, "fill_bg"))
  {
    fill_bg = !fill_bg ;
    fprintf(stderr,"%sfilling basal ganglia\n", fill_bg ? "" : "not ") ;
  }
  else if (!stricmp(option, "fillv") ||
           !stricmp(option, "fill_ventricles"))
  {
    fill_ventricles = !fill_ventricles ;
    fprintf(stderr,"%sfilling ventricles\n", fill_ventricles ? "" : "not ") ;
  }
  else switch (toupper(*option))
    {
    case 'B':
      blur_sigma = atof(argv[1]) ;
      nargs = 1 ;
      break ;
    case 'N':
      niter = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "running border classification %d times\n", niter) ;
      break ;
    case 'T':
      thickness = atoi(argv[2]) ;
      fprintf(stderr, "finding signicant strands thinner than %d mm\n",
              thickness) ;
      nargs = 1 ;
      break ;
    case 'V':
      verbose = !verbose ;
      break ;
    case 'P':
      pct = atof(argv[2]) ;
      fprintf(stderr, "using %2.0f%% threshold\n", pct*100.0f) ;
      nargs = 1 ;
      break ;
    case 'X':
      if (sscanf(argv[2], "%d", &extract) != 1)
        ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                  Progname, argv[2]) ;
      nargs = 1 ;
      break ;
    case 'W':
      wsize = atoi(argv[2]) ;
      fprintf(stderr, "using wsize = %d\n", wsize) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      usage_exit(1) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


MRI *
MRIremoveWrongDirection(MRI *mri_src, MRI *mri_dst, int wsize,
                        float low_thresh, float hi_thresh, MRI *mri_labels)
{
  MRI  *mri_kernel, *mri_smooth ;
  int  x, y, z, width, height, depth, val, nchanged, ntested ;
  float dir /*, d2I_dg2*/ ;

  mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
  fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n", blur_sigma) ;
  mri_smooth = MRIclone(mri_src, NULL) ;
  MRIconvolveGaussian(mri_src, mri_smooth, mri_kernel) ;
  MRIfree(&mri_kernel) ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

  nchanged = ntested = 0 ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (z == 101 && y == 133 && x == 152)
        {
          DiagBreak() ;
        }
        if (z == 87 && y == 88 && x == 163)
        {
          DiagBreak() ;
        }
        if (z == 88 && y == 89 && x == 163)
        {
          DiagBreak() ;
        }
        val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
        if (val >= low_thresh && val <= hi_thresh)
        {
          ntested++ ;
          dir = MRIvoxelDirection(mri_smooth, x, y, z, wsize) ;
          if (dir > 0.0)
          {
#if 0
            d2I_dg2 = MRIvoxelGradientDir2ndDerivative(mri_smooth,x,y,z,wsize);
            if (d2I_dg2 < 0)
#endif
            {
              nchanged++ ;
              val = 0 ;
              if (mri_labels)
              {
                MRIsetVoxVal(mri_labels, x, y, z, 0, 255) ;
              }
            }
          }
        }
        MRIsetVoxVal(mri_dst, x, y, z, 0,  val) ;
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "               %8d voxels tested (%2.2f%%)\n",
            ntested, 100.0f*(float)ntested/ (float)(width*height*depth));
    fprintf(stderr, "               %8d voxels changed (%2.2f%%)\n",
            nchanged, 100.0f*(float)nchanged/ (float)(width*height*depth));
  }
  MRIfree(&mri_smooth) ;
  return(mri_dst) ;
}



MRI *
MRIfillBasalGanglia(MRI *mri_src, MRI *mri_dst)
{
  float  low_thresh, hi_thresh ;
  int    total_filled, depth, height, width, x, y, z,
         xi, yi, zi, xk, yk, zk, fill, val0, val, i ;
  MRI    *mri_bg ;
  Real   tx, ty, tz ;
  MRI_SEGMENTATION  *mriseg ;
  MRI_SEGMENT       *mseg ;
  float  dx_left, dx_right, dy, dz, dist_left, dist_right ;

  if (!mri_dst)
  {
    mri_dst = MRIcopy(mri_src, NULL) ;
  }
  mri_bg = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  low_thresh = 85 ;
  hi_thresh =  105 ;
  total_filled = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == 152 && y == 117 && z == 132)  /* 93 */
        {
          DiagBreak() ;
        }
        val0 = MRIgetVoxVal(mri_src, x, y, z, 0) ;
#if 0
        if (val0 >= low_thresh && val0 <= hi_thresh &&
            MRIvox(mri_dst,x,y,z) < WM_MIN_VAL)
#else
        if (val0 >= 85 && val0 <= 105)
#endif
        {
#undef WHALF
#undef WSIZE
#define WSIZE   7
#define WHALF  ((WSIZE-1)/2)
          fill = 1 ;
          for (zk = -WHALF ; fill && zk <= WHALF ; zk++)
          {
            zi = mri_src->zi[z+zk] ;
            for (yk = -WHALF ; fill && yk <= WHALF ; yk++)
            {
              yi = mri_src->yi[y+yk] ;
              for (xk = -WHALF ; fill && xk <= WHALF ; xk++)
              {
                xi = mri_src->xi[x+xk] ;
                val = MRIgetVoxVal(mri_src, xi, yi, zi, 0) ;
                if (val < 85 || val > 110)
                {
                  fill = 0 ;  /* not homogeneous enough */
                }
              }
            }
          }
        }
        else
        {
          fill = 0 ;
        }
        if (fill)
        {
          total_filled++ ;
        }
        if (fill)
        {
          MRIsetVoxVal(mri_bg, x, y, z, 0, BASAL_GANGLIA_FILL) ;
        }
      }
    }
  }

  MRIclose(mri_bg, mri_bg) ;  /* remove small holes */

  /* segment into connected components */
  mriseg = MRIsegment(mri_bg, 1, 255) ;
  fprintf(stderr, "segmenting thick gray regions: %d %d mm segments found\n",
          mriseg->nsegments, WSIZE) ;

  /* dilate into regions that are not on */
  for (i = 0 ; i < 2*WSIZE ; i++)
  {
    MRIsegmentDilateThreshold(mriseg, mri_src, mri_src, 80, 100) ;
  }

  /* fill basal ganglia components */
  MRIclear(mri_bg) ;
  for (total_filled = i = 0 ; i < mriseg->nsegments ; i++)
  {
#define TAL_BG_LEFT_X    -30
#define TAL_BG_RIGHT_X   30
#define TAL_BG_Y         5
#define TAL_BG_Z         5
#define MAX_DIST         25

    mseg = &mriseg->segments[i] ;
    MRIvoxelToTalairach(mri_src, mseg->cx, mseg->cy, mseg->cz,&tx, &ty,&tz);
    dx_left = tx - TAL_BG_LEFT_X ;
    dx_right = tx - TAL_BG_RIGHT_X ;
    dy = ty - TAL_BG_Y ;
    dz = tz - TAL_BG_Z ;
    dist_left = sqrt(dx_left*dx_left+dy*dy+dz*dz) ;
    dist_right = sqrt(dx_right*dx_right+dy*dy+dz*dz) ;
    if (dist_left > MAX_DIST && dist_right > MAX_DIST)
    {
      continue ;
    }
    fprintf(stderr, "filling segment %d with %d voxels\n\tc = "
            "(%2.1f,%2.1f,%2.1f) tal (%2.1f,%2.1f,%2.1f), dist %2.1f,%2.1f\n",
            i, mseg->nvoxels, mseg->cx, mseg->cy, mseg->cz,
            tx, ty, tz, dist_left, dist_right) ;
    MRIsegmentToImage(mri_src, mri_bg, mriseg, i) ;
    total_filled += mseg->nvoxels ;
  }

#if 1
  /*  MRIremoveIslands(mri_bg, mri_bg, 3, .56) ;*/
  MRIbinarize(mri_bg, mri_bg, WM_MIN_VAL, 0, BASAL_GANGLIA_FILL) ;
#endif
  MRIsegmentFree(&mriseg) ;

  MRIunion(mri_dst, mri_bg, mri_dst) ;
  MRIfree(&mri_bg) ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "%d basal ganglia points filled\n", total_filled);
  }

  return(mri_dst) ;
}

MRI *
MRIfilterMorphology(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, xk, yk, zk, xi, yi, zi, nsix, nvox,
          total ;
  BUFTYPE val ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  /*  mri_dst = MRIremove1dStructures(mri_src, mri_dst, 10000, 2) ;*/

  total = 0 ;
  do
  {
    nvox = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
          if (x == 159 && y == 138 && z == 101)
          {
            DiagBreak() ;
          }
          if (val < WM_MIN_VAL)
          {
            continue ;
          }
          nsix = 0 ;
          for (zk = -1 ; zk <= 1 ; zk++)
          {
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              for (xk = -1 ; xk <= 1 ; xk++)
              {
                if ((fabs(xk) + fabs(yk) + fabs(zk)) <= 1)
                {
                  continue ;
                }
                xi = mri_dst->xi[x+xk] ;
                yi = mri_dst->yi[y+yk] ;
                zi = mri_dst->zi[z+zk] ;
                if (xi == 159 && yi == 138 && zi == 101)
                {
                  DiagBreak() ;
                }
                if (MRIgetVoxVal(mri_dst, xi, yi, zi, 0) >= WM_MIN_VAL)
                {
                  if (xk && MRIgetVoxVal(mri_dst, xi, y, z, 0) >= WM_MIN_VAL)
                  {
                    continue ;
                  }
                  if (yk && MRIgetVoxVal(mri_dst, x, yi, z, 0) >= WM_MIN_VAL)
                  {
                    continue ;
                  }
                  if (zk && MRIgetVoxVal(mri_dst, x, y, zi, 0) >= WM_MIN_VAL)
                  {
                    continue ;
                  }
                  if (xk)
                  {
                    if (xi == 141 && y == 132 && z == 30)
                    {
                      DiagBreak() ;
                    }
                    MRIsetVoxVal(mri_dst, xi, y, z, 0, DIAGONAL_FILL) ;
                    if (is_diagonal(mri_dst, xi, y, z))
                    {
                      MRIsetVoxVal(mri_dst, xi, y, z, 0, 0) ;
                    }
                    else
                    {
                      nvox++ ;
                    }
                  }
                  if (yk)
                  {
                    MRIsetVoxVal(mri_dst, x, yi, z, 0, DIAGONAL_FILL) ;
                    if (is_diagonal(mri_dst, x, yi, z))
                    {
                      MRIsetVoxVal(mri_dst, x, yi, z, 0, 0) ;
                    }
                    else
                    {
                      nvox++ ;
                    }
                  }
                  if (zk)
                  {
                    MRIsetVoxVal(mri_dst, x, y, zi, 0, DIAGONAL_FILL) ;
                    if (is_diagonal(mri_dst, x, y, zi))
                    {
                      MRIsetVoxVal(mri_dst, x, y, zi, 0, 0) ;
                    }
                    else
                    {
                      nvox++ ;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    total += nvox ;
  }
  while (nvox > 0) ;
  fprintf(stderr, "%d diagonally connected voxels added...\n", total) ;

  return(mri_dst) ;
}

MRI *
MRIremove1dStructures(MRI *mri_src, MRI *mri_dst, int max_iter, int thresh,
                      MRI *mri_labels)
{
  int     width, height, depth, x, y, z, xk, yk, zk, xi, yi, zi, nsix, nvox,
          total, niter ;
  BUFTYPE val ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  MRIcopy(mri_src, mri_dst) ;

  if (max_iter == 0)
  {
    max_iter = 1000 ;
  }

  niter = total = 0 ;
  do
  {
    nvox = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
          nsix = 0 ;
          if (val < WM_MIN_VAL)
          {
            continue ;
          }

          for (zk = -1 ; zk <= 1 ; zk++)
          {
            zi = z+zk ;
            if (zi < 0 || zi >= depth)
            {
              continue ;
            }
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = y+yk ;
              if (yi < 0 || yi >= height)
              {
                continue ;
              }
              for (xk = -1 ; xk <= 1 ; xk++)
              {
#if 0
                if ((z == 124 || z == 125) && y == 162 && x == 158)
                {
                  DiagBreak() ;
                }
                else
#endif
                  if ((fabs(xk) + fabs(yk) + fabs(zk)) != 1)
                  {
                    continue ;
                  }
                xi = x+xk ;
                if (xi < 0 || xi >= width)
                {
                  continue ;
                }
                if (MRIgetVoxVal(mri_dst, xi, yi, zi,0) >= WM_MIN_VAL)
                {
                  nsix++ ;
                }
              }
            }
          }
          if (nsix < thresh && val >= WM_MIN_VAL)
          {
            if ((z == 124 || z == 125) && y == 162 && x == 158)
            {
              DiagBreak() ;
            }
            nvox++ ;
            MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
            if (mri_labels)
            {
              MRIsetVoxVal(mri_labels,x,y,z,0,255) ;
            }
          }
        }
      }
    }
    total += nvox ;
  }
  while (nvox > 0 || ++niter > max_iter) ;

  if (width == 256)
  {
    fprintf(stderr, "%d sparsely connected voxels removed...\n", total) ;
  }

  return(mri_dst) ;
}

static int
is_diagonal(MRI *mri, int x, int y, int z)
{
  int xk, yk, zk, xi, yi, zi ;

  for (zk = -1 ; zk <= 1 ; zk++)
  {
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        if ((fabs(xk) + fabs(yk) + fabs(zk)) <= 1)
        {
          continue ;
        }
        xi = mri->xi[x+xk] ;
        yi = mri->yi[y+yk] ;
        zi = mri->zi[z+zk] ;
        if (!MRIgetVoxVal(mri, xi, yi, zi, 0))
        {
          continue ;  /* not a diagonal */
        }

        if (xk && MRIgetVoxVal(mri, xi, y, z, 0) >= WM_MIN_VAL)
        {
          continue ;  /* not a diagonal */
        }
        if (yk && MRIgetVoxVal(mri, x, yi, z, 0) >= WM_MIN_VAL)
        {
          continue ;  /* not a diagonal */
        }
        if (zk && MRIgetVoxVal(mri, x, y, zi, 0) >= WM_MIN_VAL)
        {
          continue ;  /* not a diagonal */
        }
        return(1) ;   /* found a diagonal with no 4-connection */
      }
    }
  }
  return(0) ;
}


MRI *
MRIfillVentricles(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, xk, yk, xi, yi, nfilled, total, s ;
  MRI     *mri_filled, *mri_ventricles = NULL ;
  MRI_SEGMENTATION *mriseg ;
  MRI_SEGMENT      *mseg ;
  Real    xt, yt, zt ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  MRIcopy(mri_src, mri_dst) ;
  mri_filled = MRIcopy(mri_src, NULL) ;

  MRIreplaceValues(mri_filled, mri_filled, VENTRICLE_FILL,
                   VENTRICLE_FILL-1) ;

  /* first fill each coronal slice starting from a background seed */
  for (z = 0 ; z < depth ; z++)
  {
    total = 0 ;
    do
    {
      nfilled = 0 ;
      MRIsetVoxVal(mri_filled, 0, 0, z, 0, VENTRICLE_FILL) ;
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (MRIgetVoxVal(mri_filled, x, y, z, 0) == VENTRICLE_FILL)
          {
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_src->yi[y+yk] ;
              for (xk = -1 ; xk <= 1 ; xk++)
              {
                xi = mri_src->xi[x+xk] ;
                if (!MRIgetVoxVal(mri_filled, xi, yi, z, 0))
                {
                  nfilled++ ;
                  MRIsetVoxVal(mri_filled, xi, yi, z, 0, VENTRICLE_FILL) ;
                }
              }
            }
          }
        }
      }
      total += nfilled ;
    }
    while (nfilled > 0) ;
  }

  MRIcomplement(mri_filled, mri_filled) ;
  MRIreplaceValues(mri_filled, mri_filled, 1, VENTRICLE_FILL) ;
  mriseg = MRIsegment(mri_filled, 1, 255) ;
  fprintf(stderr, "%d segments found...\n", mriseg->nsegments) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_filled, "exterior_filled.mgh") ;
  }
  for (s = 0 ; s < mriseg->nsegments ; s++)
  {
    mseg = &mriseg->segments[s] ;
    if (mseg->nvoxels < 100 ||
        mseg->z1-mseg->z0 < 7)
    {
      continue ;
    }
    MRIvoxelToTalairach(mri_src, mseg->cx, mseg->cy, mseg->cz, &xt, &yt, &zt);
    fprintf(stderr, "added segment %d, nvox=%d, bbox [%d:%d, %d:%d, %d:%d]\n"
            "\tc = %2.1f, %2.1f, %2.1f (tal = %2.1f, %2.1f, %2.1f)\n",
            s, mseg->nvoxels, mseg->x0, mseg->x1, mseg->y0,mseg->y1,mseg->z0,
            mseg->z1, mseg->cx, mseg->cy, mseg->cz, xt, yt, zt) ;
    mri_ventricles = MRIsegmentToImage(mri_filled, mri_ventricles, mriseg, s) ;
  }

  MRIfree(&mri_filled) ;
  MRIsegmentFree(&mriseg) ;

  /* remove voxels close to the midline so that the cc can still be found */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        MRIvoxelToTalairach(mri_src, x, y, z, &xt, &yt, &zt);
        if (fabs(xt) < 5)
        {
          MRIsetVoxVal(mri_ventricles, x, y, z, 0, 0) ;
        }
      }
    }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_ventricles, "ventricles.mgh") ;
  }
  MRIunion(mri_ventricles, mri_dst, mri_dst) ;

  MRIfree(&mri_ventricles) ;
  return(mri_dst) ;
}
static MRI *
MRIrecoverBrightWhite(MRI *mri_T1, MRI *mri_src, MRI *mri_dst,
                      float wm_low, float wm_hi,
                      float slack, float pct_thresh)
{
  int     x, y, z, xk, yk, zk, xi, yi, zi, width, height, depth, nwhite,
          ntested, nchanged ;
  BUFTYPE val ;
  float   intensity_thresh, nvox_thresh ;

  if (!mri_dst)
  {
    mri_dst = MRIcopy(mri_src, NULL) ;
  }
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  ntested = nchanged = 0 ;
  intensity_thresh = wm_hi+slack ;
  nvox_thresh = (3*3*3-1)*pct_thresh ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
        if (val > wm_hi && val <= intensity_thresh &&
            MRIgetVoxVal(mri_src, x, y, z, 0) < WM_MIN_VAL)
        {
          ntested++ ;
          nwhite = 0 ;
          for (zk = -1 ; zk <= 1 ; zk++)
          {
            zi = mri_src->zi[z+zk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_src->yi[y+yk] ;
              for (xk = -1 ; xk <= 1 ; xk++)
              {
                xi = mri_src->xi[x+xk] ;
                val = MRIgetVoxVal(mri_T1, xi, yi, zi, 0);
                if (val >= wm_low && val <= wm_hi)
                {
                  nwhite++ ;
                }
              }
            }
          }
          if (nwhite >= nvox_thresh)
          {
            nchanged++ ;
            MRIsetVoxVal(mri_dst, x, y, z, 0, MRIgetVoxVal(mri_T1, x, y, z,0)) ;
          }
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "               %8d voxels tested (%2.2f%%)\n",
            ntested, 100.0f*(float)ntested/ (float)(width*height*depth));
    fprintf(stderr, "               %8d voxels changed (%2.2f%%)\n",
            nchanged, 100.0f*(float)nchanged/ (float)(width*height*depth));
  }
  return(mri_dst) ;
}

#if 0
#define NPTS (sizeof(xpts) / sizeof(xpts[0]))
static int xpts[] =
{
  0,   1, 1, 1, 0, -1, -1, -1
} ;
static int ypts[] =
{
  -1, -1, 0, 1, 1,  1,  0, -1
} ;

static int
MRIcheckRemovals(MRI *mri_T1, MRI *mri_dst, MRI *mri_labels, int wsize)
{
  int    x, y, z, width, depth, height, whalf, ntested, nchanged, on, vertex;
  MRI    *mri_tmp, *mri_region, *mri_plane, *mri_binary_plane ;
  float  min_on ;

  whalf = (wsize-1)/2 ;

  mri_tmp = MRIcopy(mri_dst, NULL) ;
  mri_region = MRIalloc(wsize, wsize, wsize, MRI_UCHAR) ;
  min_on = .1*wsize*wsize ;
  MRIcopyLabel(mri_labels, mri_tmp, 255) ;

  MRIbinarize(mri_tmp, mri_tmp, WM_MIN_VAL, 0, 100) ;
  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  ntested = nchanged = 0 ;
  if (Gdiag == 99)
  {
    MRIwrite(mri_tmp, "tmp.mgh") ;
  }
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (z == 87 && y == 88 && x == 163)  /* test1 cs filled */
        {
          DiagBreak() ;
        }
        if (z == 88 && y == 89 && x == 163)  /* test1 cs filled */
        {
          DiagBreak() ;
        }

        if (z == 101 && y == 133 && x == 152)
        {
          DiagBreak() ;
        }
        if (x == 157 && y == 143 && z == 98)
        {
          DiagBreak() ;
        }
        if (x == 156 && y == 143 && z == 98)
        {
          DiagBreak() ;
        }

        if (x == 154 && y == 167 && z == 128)
        {
          DiagBreak() ;
        }
        if (x == 136 && y == 147 && z == 28)
        {
          DiagBreak() ;
        }

        if (x == 163 && y == 88 && z == 86)
        {
          DiagBreak() ;
        }

        if ((x == 140 && y == 141 && z == 54) ||
            (x == 140 && y == 141 && z == 53) ||
            (x == 140 && y == 142 && z == 53) ||
            (x == 140 && y == 142 && z == 54) ||
            (x == 140 && y == 140 && z == 53))
        {
          DiagBreak() ;  /* test4 cerebellum */
        }
        if (x == 142 && y == 139 && z == 54)   /* test4 */
        {
          DiagBreak() ;
        }

        if (!MRIgetVoxVal(mri_labels, x, y, z, 0))
        {
          continue ;
        }
        ntested++ ;

        MRIextract(mri_tmp, mri_region, x-whalf,y-whalf,z-whalf,
                   wsize, wsize, wsize) ;

        vertex =
          MRIcountCpolvOnAtVoxel(mri_region, whalf, whalf, whalf, wsize, &on) ;

        mri_plane = MRIextractVertexPlane(mri_tmp, NULL, vertex,x,y,z,wsize);
        MRIthreshold(mri_plane, mri_plane, 50) ;
        MRIremove1dStructures(mri_plane,mri_plane, 10000,2,NULL);
        mri_binary_plane = MRIfillFG(mri_plane, NULL, whalf, whalf, 0,
                                     50, 128, &on) ;
        if (on > min_on)
        {
          int  xk, yk, i, ntransitions, i_prev  ;

          /*
             now look at the winding # (number of white-black transitions
             in a circle around the central point
          */

          ntransitions = 0 ;
          for (i = 0 ; i < NPTS ; i++)
          {
            xk = xpts[i] ;
            yk = ypts[i] ;
            i_prev = i-1 ;
            if (i_prev < 0)
            {
              i_prev = NPTS-1 ;
            }
            if (MRIgetVoxVal(mri_binary_plane, whalf+xpts[i], whalf+ypts[i], 0, 0) !=
                MRIgetVoxVal(mri_binary_plane,whalf+xpts[i_prev],
                             whalf+ypts[i_prev],0, 0))
            {
              ntransitions++ ;
            }
          }
          if (ntransitions > 2)   /* not planar */
          {
            nchanged++ ;
            MRIsetVoxVal(mri_dst, x, y, z, 0, MRIgetVoxVal(mri_T1, x, y, z, 0)) ;
          }
        }
        if (Gdiag & DIAG_WRITE)
        {
          MRIwrite(mri_region, "region.mgh") ;
          MRIwrite(mri_plane, "plane.mgh") ;
          MRIwrite(mri_binary_plane, "binary_plane.mgh") ;
        }
        MRIfree(&mri_plane) ;
        MRIfree(&mri_binary_plane) ;
      }
    }
  }
  MRIfree(&mri_tmp) ;
  MRIfree(&mri_region) ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "               %8d voxels tested (%2.2f%%)\n",
            ntested, 100.0f*(float)ntested/ (float)(width*height*depth));
    fprintf(stderr, "               %8d voxels restored (%2.2f%%)\n",
            nchanged, 100.0f*(float)nchanged/ (float)(width*height*depth));
  }
  return(NO_ERROR) ;
}


MRI *
MRIremoveFilledBrightStuff(MRI *mri_T1, MRI *mri_labeled, MRI *mri_dst,
                           int filled_label, float thresh)
{
  int     x, y, z, width, height, depth, nwhite, ntested, nchanged ;
  BUFTYPE val ;
  float   intensity_thresh ;

  if (!mri_dst)
  {
    mri_dst = MRIcopy(mri_labeled, NULL) ;
  }
  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  ntested = nchanged = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
        ntested++ ;
        if (MRIgetVoxVal(mri_labeled, x, y, z, 0) == filled_label &&
            val > thresh)
        {
          MRIsetVoxVal(mri_labeled, x, y, z, 0) ;
          nchanged++ ;
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "               %8d voxels tested (%2.2f%%)\n",
            ntested, 100.0f*(float)ntested/ (float)(width*height*depth));
    fprintf(stderr, "               %8d voxels changed (%2.2f%%)\n",
            nchanged, 100.0f*(float)nchanged/ (float)(width*height*depth));
  }
  return(mri_dst) ;
}

#endif
MRI *
MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm)
{
  int     width, height, depth, x, y, z, nlabeled, nwhite,
          xk, yk, zk, xi, yi, zi;
  BUFTYPE val, wm ;
  MRI     *mri_labeled, *mri_tmp ;

  mri_labeled = MRIclone(mri_T1, NULL) ;
  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
        wm = MRIgetVoxVal(mri_wm, x, y, z, 0) ;

        if (x == 110 && y == 125 && z == 172)  /* T1=148 */
        {
          DiagBreak() ;
        }
        /* not white matter and bright (e.g. eye sockets) */
        if ((wm < WM_MIN_VAL) && (val > 125))
        {
          nwhite = 0 ;
          for (xk = -1 ; xk <= 1 ; xk++)
          {
            xi = mri_T1->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_T1->yi[y+yk] ;
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = mri_T1->zi[z+zk] ;
                if (MRIgetVoxVal(mri_wm, xi, yi, zi, 0) >= WM_MIN_VAL)
                {
                  nwhite++ ;
                }
              }
            }
          }
#define MIN_WHITE  ((3*3*3-1)/2)
          if (nwhite < MIN_WHITE)
          {
            MRIsetVoxVal(mri_labeled, x, y, z, 0, BRIGHT_LABEL) ;
          }
        }
      }
    }
  }

  /* find all connected voxels that are above 115 */
  MRIdilateThreshLabel(mri_labeled, mri_T1, NULL, BRIGHT_LABEL, 10,115);
  MRIclose(mri_labeled, mri_labeled) ;

  /* expand once more to all neighboring voxels that are bright. At
     worst we will erase one voxel of white matter.
  */
  mri_tmp =
    MRIdilateThreshLabel(mri_labeled, mri_T1, NULL, BRIGHT_LABEL,1,100);
  MRIxor(mri_labeled, mri_tmp, mri_tmp, 1, 255) ;
  MRIreplaceValues(mri_tmp, mri_tmp, 1, BRIGHT_BORDER_LABEL) ;
  MRIunion(mri_tmp, mri_labeled, mri_labeled) ;
#if 0
  fprintf(stderr, "selectively smoothing volume....\n") ;
  MRIsoapBubbleLabel(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 200) ;
#endif
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_labeled, "label.mgh") ;
  }
  /*    MRIwrite(mri_tmp, "tmp.mgh") ;*/
  nlabeled = MRIvoxelsInLabel(mri_labeled, BRIGHT_LABEL) ;
  fprintf(stderr, "%d bright non-wm voxels segmented.\n", nlabeled) ;

  MRIfree(&mri_tmp) ;
  return(mri_labeled) ;
}


#include "mri_segment.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_segment_help_xml,
                mri_segment_help_xml_len);
  exit(code);
}
