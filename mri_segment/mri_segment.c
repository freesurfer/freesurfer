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
#include "mrinorm.h"
#include "timer.h"

static int extract = 0 ;
static int  verbose = 0 ;
static int wsize = 11 ;
static float pct = 0.8 ;
static float pslope = 1.0f ;
static float nslope = 1.0f ;
static float wm_low = 90 ;
static float wm_hi = 125 ;
static float gray_hi = 100 ;
static int niter = 1 ;
static int gray_hi_set = 0 ;
static int wm_hi_set = 0 ;
static int wm_low_set = 0 ;

static int thickness = 4 ;
static int thicken = 1 ;
static int nsegments = 20 ;
static int fill_bg = 0 ;
static int fill_ventricles = 0 ;

static int keep_edits = 0 ;

static int auto_detect_stats =  1 ;
static int log_stats = 1 ;

char *Progname ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *MRIremoveWrongDirection(MRI *mri_src, MRI *mri_dst, int wsize,
                             float low_thresh, float hi_thresh) ;
MRI *MRIfilterMorphology(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIfillBasalGanglia(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIfillVentricles(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIremove1dStructures(MRI *mri_src, MRI *mri_dst, int max_iter,
                           int thresh) ;
static int is_diagonal(MRI *mri, int x, int y, int z) ;

int
main(int argc, char *argv[])
{
  MRI     *mri_src, *mri_dst, *mri_tmp, *mri_labeled ;
  char    *input_file_name, *output_file_name ;
  int     nargs, i, msec ;
  struct timeb  then ;

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
    ErrorExit(ERROR_BADPARM,
              "usage: %s <input volume> <output volume>", Progname);

  TimerStart(&then) ;
  input_file_name = argv[1] ;
  output_file_name = argv[2] ;

  mri_src = MRIread(input_file_name) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read source volume from %s",
              Progname, input_file_name) ;

  if (auto_detect_stats) /* widen range to allow for more variability */
    wm_low -= 10 ;  
  fprintf(stderr, "doing initial intensity segmentation...\n") ;
  mri_tmp = MRIintensitySegmentation(mri_src, NULL, wm_low, wm_hi, gray_hi);

  fprintf(stderr, "using local statistics to label ambiguous voxels...\n") ;
  MRIhistoSegment(mri_src, mri_tmp, wm_low, wm_hi, gray_hi, wsize, 3.0f) ;

  if (auto_detect_stats)
  {
    float   white_mean, white_sigma, gray_mean, gray_sigma ;

    fprintf(stderr, "computing class statistics for intensity windows...\n") ;
    MRIcomputeClassStatistics(mri_src, mri_tmp, 30, WHITE_MATTER_MEAN,
                              &white_mean, &white_sigma, &gray_mean,
                              &gray_sigma) ;
    
    if (!wm_low_set)
      wm_low = gray_mean + gray_sigma ;
    else
      wm_low += 10 ;   /* set it back to it's original value */
    if (!gray_hi_set)
      gray_hi = gray_mean + 2*gray_sigma ;
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

  fprintf(stderr, 
          "using local geometry to label remaining ambiguous voxels...\n") ;
  mri_labeled = MRIcpolvMedianCurveSegment(mri_src, mri_tmp, NULL, 5, 3,
                                           gray_hi, wm_low);
  fprintf(stderr, 
          "\nreclassifying voxels using Gaussian border classifier...\n") ;

  /*
    now use the gray and white mattter border voxels to build a Gaussian
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
  fprintf(stderr, 
          "\nremoving voxels with positive offset direction...\n") ;
  MRIremoveWrongDirection(mri_dst, mri_dst, 3, wm_low-5, gray_hi) ;

  if (thicken)
  {
    fprintf(stderr, "thickening thin strands....\n") ;
    /*    MRIfilterMorphology(mri_dst, mri_dst) ;*/
    MRIremove1dStructures(mri_dst, mri_dst, 10000, 2) ;
    MRIthickenThinWMStrands(mri_dst, mri_dst, thickness, nsegments) ;
  }

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
      ErrorPrintf(ERROR_NOFILE, "%s: could not read file %s to preserve edits",
                Progname, output_file_name) ;
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
  if (!stricmp(option, "slope"))
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
    thicken = !thicken ;
    fprintf(stderr,"%sthickening thin strands\n", thicken ? "" : "not ") ;
  }
  else if (!stricmp(option, "fillbg") || !stricmp(option, "fill_bg"))
  {
    fill_bg = !fill_bg ;
    fprintf(stderr,"%sfilling basal ganglia\n", fill_bg ? "" : "not ") ;
  }
  else if (!stricmp(option, "fillv") || !stricmp(option, "fill_ventricles"))
  {
    fill_ventricles = !fill_ventricles ;
    fprintf(stderr,"%sfilling ventricles\n", fill_ventricles ? "" : "not ") ;
  }
  else switch (toupper(*option))
  {
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
    break ;
  case '?':
  case 'U':
    fprintf(stderr,
            "usage: %s <classifier file> <input volume> <output volume>\n",
            Progname) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
#define BLUR_SIGMA 0.25f
MRI *
MRIremoveWrongDirection(MRI *mri_src, MRI *mri_dst, int wsize, 
                        float low_thresh, float hi_thresh)
{
  MRI  *mri_kernel, *mri_smooth ;
  int  x, y, z, width, height, depth, val, nchanged, ntested ;

  mri_kernel = MRIgaussian1d(BLUR_SIGMA, 100) ;
  fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n", BLUR_SIGMA) ;
  mri_smooth = MRIclone(mri_src, NULL) ;
  MRIconvolveGaussian(mri_src, mri_smooth, mri_kernel) ;
  MRIfree(&mri_kernel) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  nchanged = ntested = 0 ;
  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ; 
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIvox(mri_src, x, y, z) ;
        if (val >= low_thresh && val <= hi_thresh)
        {
          ntested++ ;
          if (MRIvoxelDirection(mri_smooth, x, y, z, wsize) > 0)
          {
            nchanged++ ;
            val = 0 ;
          }
        }
        MRIvox(mri_dst, x, y, z) = val ;
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
    mri_dst = MRIcopy(mri_src, NULL) ;
  mri_bg = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;
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
          DiagBreak() ;
        val0 = MRIvox(mri_src, x, y, z) ;
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
                val = MRIvox(mri_src, xi, yi, zi) ;
                if (val < 85 || val > 110)
                  fill = 0 ;   /* not homogeneous enough */
              }
            }
          }
        }
        else
          fill = 0 ;
        if (fill)
          total_filled++ ;
        if (fill)
          MRIvox(mri_bg, x, y, z) = BASAL_GANGLIA_FILL ;
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
    MRIsegmentDilateThreshold(mriseg, mri_src, mri_src, 80, 100) ;

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
      continue ;
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
    fprintf(stderr, "%d basal ganglia points filled\n", total_filled);

  return(mri_dst) ;
}

MRI *
MRIfilterMorphology(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, xk, yk, zk, xi, yi, zi, nsix, nvox,
          total ;
  BUFTYPE val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  mri_dst = MRIremove1dStructures(mri_src, mri_dst, 10000, 2) ;

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
          val = MRIvox(mri_dst, x, y, z) ;
          if (x == 159 && y == 138 && z == 101)
            DiagBreak() ;
          if (val < WM_MIN_VAL)
            continue ;
          nsix = 0 ;
          for (zk = -1 ; zk <= 1 ; zk++)
          {
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              for (xk = -1 ; xk <= 1 ; xk++)
              {
                if ((fabs(xk) + fabs(yk) + fabs(zk)) <= 1)
                  continue ;
                xi = mri_dst->xi[x+xk] ;
                yi = mri_dst->yi[y+yk] ;
                zi = mri_dst->zi[z+zk] ;
                if (xi == 159 && yi == 138 && zi == 101)
                  DiagBreak() ;
                if (MRIvox(mri_dst, xi, yi, zi) >= WM_MIN_VAL)
                {
                  if (xk && MRIvox(mri_dst, xi, y, z) >= WM_MIN_VAL)
                    continue ;
                  if (yk && MRIvox(mri_dst, x, yi, z) >= WM_MIN_VAL)
                    continue ;
                  if (zk && MRIvox(mri_dst, x, y, zi) >= WM_MIN_VAL)
                    continue ;
                  if (xk)
                  {
                    if (xi == 141 && y == 132 && z == 30)
                      DiagBreak() ;
                    MRIvox(mri_dst, xi, y, z) = DIAGONAL_FILL ;
                    if (is_diagonal(mri_dst, xi, y, z))
                      MRIvox(mri_dst, xi, y, z) = 0 ;
                    else
                      nvox++ ;
                  }
                  if (yk)
                  {
                    MRIvox(mri_dst, x, yi, z) = DIAGONAL_FILL ;
                    if (is_diagonal(mri_dst, x, yi, z))
                      MRIvox(mri_dst, x, yi, z) = 0 ;
                    else
                      nvox++ ;
                  }
                  if (zk)
                  {
                    MRIvox(mri_dst, x, y, zi) = DIAGONAL_FILL ;
                    if (is_diagonal(mri_dst, x, y, zi))
                      MRIvox(mri_dst, x, y, zi) = 0 ;
                    else
                      nvox++ ;
                  }
                }
              }
            }
          }
        }
      }
    }
    total += nvox ;
  } while (nvox > 0) ;
  fprintf(stderr, "%d diagonally connected voxels added...\n", total) ;

  return(mri_dst) ;
}

MRI *
MRIremove1dStructures(MRI *mri_src, MRI *mri_dst, int max_iter, int thresh)
{
  int     width, height, depth, x, y, z, xk, yk, zk, xi, yi, zi, nsix, nvox,
          total, niter ;
  BUFTYPE val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  MRIcopy(mri_src, mri_dst) ;

  if (max_iter == 0)
    max_iter = 1000 ;

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
          val = MRIvox(mri_dst, x, y, z) ;
          nsix = 0 ;
          if (val < WM_MIN_VAL)
            continue ;

          for (zk = -1 ; zk <= 1 ; zk++)
          {
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              for (xk = -1 ; xk <= 1 ; xk++)
              {
                if ((fabs(xk) + fabs(yk) + fabs(zk)) != 1)
                  continue ;
                xi = mri_dst->xi[x+xk] ;
                yi = mri_dst->yi[y+yk] ;
                zi = mri_dst->zi[z+zk] ;
                if (MRIvox(mri_dst, xi, yi, zi) >= WM_MIN_VAL)
                  nsix++ ;
              }
            }
          }
          if (nsix < thresh && val >= WM_MIN_VAL)
          {
            nvox++ ;
            MRIvox(mri_dst, x, y, z) = 0 ;
          }
        }
      }
    }
    total += nvox ;
  } while (nvox > 0 || ++niter > max_iter) ;

  fprintf(stderr, "%d sparsely connected voxels removed...\n", total) ;

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
          continue ;
        xi = mri->xi[x+xk] ;
        yi = mri->yi[y+yk] ;
        zi = mri->zi[z+zk] ;
        if (!MRIvox(mri, xi, yi, zi))
          continue ;   /* not a diagonal */

        if (xk && MRIvox(mri, xi, y, z) >= WM_MIN_VAL)
          continue ;   /* not a diagonal */
        if (yk && MRIvox(mri, x, yi, z) >= WM_MIN_VAL)
          continue ;  /* not a diagonal */
        if (zk && MRIvox(mri, x, y, zi) >= WM_MIN_VAL)
          continue ;  /* not a diagonal */
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
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

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
      MRIvox(mri_filled, 0, 0, z) = VENTRICLE_FILL ;
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (MRIvox(mri_filled, x, y, z) == VENTRICLE_FILL)
          {
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_src->yi[y+yk] ;
              for (xk = -1 ; xk <= 1 ; xk++)
              {
                xi = mri_src->xi[x+xk] ;
                if (!MRIvox(mri_filled, xi, yi, z))
                {
                  nfilled++ ;
                  MRIvox(mri_filled, xi, yi, z) = VENTRICLE_FILL ;
                }
              }
            }
          }
        }
      }
      total += nfilled ;
    } while (nfilled > 0) ;
  } 

  MRIcomplement(mri_filled, mri_filled) ; 
  MRIreplaceValues(mri_filled, mri_filled, 1, VENTRICLE_FILL) ;
  mriseg = MRIsegment(mri_filled, 1, 255) ;
  fprintf(stderr, "%d segments found...\n", mriseg->nsegments) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_filled, "exterior_filled.mgh") ;
  for (s = 0 ; s < mriseg->nsegments ; s++)
  {
    mseg = &mriseg->segments[s] ;
    if (mseg->nvoxels < 100 ||
        mseg->z1-mseg->z0 < 7)
      continue ;
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
          MRIvox(mri_ventricles, x, y, z) = 0 ;
      }
    }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_ventricles, "ventricles.mgh") ;
  MRIunion(mri_ventricles, mri_dst, mri_dst) ;

  MRIfree(&mri_ventricles) ;
  return(mri_dst) ;
}
