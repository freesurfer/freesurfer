/**
 * @brief Normalize a set of longituindal volumes 
 *
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



static char *mask_fname = NULL ;

static double bias_sigma = 4.0 ;
static FILE *diag_fp = NULL ;

static int normalize_timepoints_with_parzen_window(MRI *mri, double cross_time_sigma) ;
static int normalize_timepoints(MRI *mri, double thresh, double cross_time_sigma) ;
static void usage_exit(int code) ;
static int get_option(int argc, char *argv[]) ;


/*
  command line consists of these inputs:

  argv[1]  - input volume
  argv[2]  - transform (lta/xfm/m3d)
  argv[3]  - output volume
*/


static double cross_time_sigma = 1.0 ;


#define MAX_TIMEPOINTS 2000
static char *subjects[MAX_TIMEPOINTS] ;
int
main(int argc, char *argv[])
{
  char         *in_fname, *out_fname, **av, *xform_fname, fname[STRLEN] ;
  MRI          *mri_in, *mri_tmp ;
  int          ac, nargs, msec, minutes, seconds;
  int          input, ninputs ;
  Timer start ;
  TRANSFORM    *transform = NULL ;
  char         line[STRLEN], *cp, subject[STRLEN], sdir[STRLEN], base_name[STRLEN] ;
  FILE         *fp ;

  std::string cmdline = getAllInfo(argc, argv, "mri_fuse_intensity_images");

  nargs = handleVersionOption(argc, argv, "mri_fuse_intensity_images");
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

  if (argc < 5)
    ErrorExit
      (ERROR_BADPARM,
       "usage: %s [<options>] <longitudinal time point file> <in vol> <transform file> <out vol> \n",
       Progname) ;
  in_fname = argv[2] ;
  xform_fname = argv[3] ;
  out_fname = argv[4] ;

  transform = TransformRead(xform_fname) ;
  if (transform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, xform_fname) ;
  start.reset() ;

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
    sprintf(subject, "%s.long.%s", subjects[input], base_name) ;
    printf("reading subject %s - %d of %d\n", subject, input+1, ninputs) ;
    sprintf(fname, "%s/%s/mri/%s", sdir, subject, in_fname) ;
    mri_tmp = MRIread(fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s",
                Progname, fname) ;
    MRImakePositive(mri_tmp, mri_tmp) ;
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
  }
  MRIaddCommandLine(mri_in, cmdline) ;

  // try to bring the images closer to each other at each voxel where they seem to come from the same distribution
  {
    MRI   *mri_frame1, *mri_frame2 ;
    double rms_after ;

    mri_frame1 = MRIcopyFrame(mri_in, NULL, 0, 0) ;
    mri_frame2 = MRIcopyFrame(mri_in, NULL, 1, 0) ;
    rms_after = MRIrmsDiff(mri_frame1, mri_frame2) ;
    printf("RMS before intensity cohering  = %2.2f\n", rms_after) ;
    MRIfree(&mri_frame1) ; MRIfree(&mri_frame2) ; 
    if (0)
      normalize_timepoints(mri_in, 2.0, cross_time_sigma) ;
    else
      normalize_timepoints_with_parzen_window(mri_in, cross_time_sigma) ;
      
    mri_frame1 = MRIcopyFrame(mri_in, NULL, 0, 0) ;
    mri_frame2 = MRIcopyFrame(mri_in, NULL, 1, 0) ;
    rms_after = MRIrmsDiff(mri_frame1, mri_frame2) ;
    MRIfree(&mri_frame1) ; MRIfree(&mri_frame2) ;
    printf("RMS after intensity cohering  = %2.2f (sigma=%2.2f)\n", rms_after, cross_time_sigma) ;
  }

  for (input = 0 ; input < ninputs ; input++)
  {
    sprintf(fname, "%s/%s.long.%s/mri/%s", sdir, subjects[input], base_name, out_fname) ;
    printf("writing normalized volume to %s...\n", fname) ;
    if (MRIwriteFrame(mri_in, fname, input)  != NO_ERROR)
      ErrorExit(ERROR_BADFILE, "%s: could not write normalized volume to %s",Progname, fname);
  }

  MRIfree(&mri_in) ;

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
  if (!strcmp(option, "-HELP")||!strcmp(option, "-USAGE"))
    usage_exit(0) ;
  else if (!strcmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!strcmp(option, "SIGMA"))
  {
    bias_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("smoothing bias field with sigma = %2.1f\n", bias_sigma) ;
  }
  else if (!stricmp(option, "cross_time_sigma") || !stricmp(option, "cross-time-sigma"))
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
  else switch (*option)
  {
  case 'W':
    Gdiag |= DIAG_WRITE ;
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
  default:
    printf("unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}


static void
usage_exit(int code)
{
  //  outputHelp(Progname);
#ifdef GREGT
  printf("usage: %s [<options>] <in volume> <transform> <normalized volume>\n", Progname) ;
  printf("\t<transform>                  ex. transforms/talairach.lta "
         "(or 'noxform') \n");
  printf("\noptions:\n");
  printf("\t-seg <filename>              aseg file, to help normalization\n");
  printf("\t-sigma <bias sigma>          smoothing sigma for bias field if control points specified (def=4)\n");
  printf("\t-mask <mri_vol>              use mri_vol to mask input\n");
  printf("\t-diag <filename>             write to log file\n");
  printf("\t-debug_voxel <x> <y> <z>     debug voxel\n");
  printf("\t-debug_node <x> <y> <z>      debug node\n");
  printf("\t-tr <float n>                set TR in msec\n");
  printf("\t-te <float n>                set TE in msec\n");
  printf("\t-alpha <float n>             set alpha in radians\n");
  printf("\t-v <int n>                   does nothing as far "
         "as i can tell, but an option\n");
#endif
  exit(code) ;
}

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
    //    MRIsoapBubble(mri_bias, mri_ctrl, mri_bias, nsoap) ;
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



static int
normalize_timepoints_with_parzen_window(MRI *mri, double cross_time_sigma)
{
  int   frame1, frame2, x, y, z ;
  double val0, val, total, g, norm, total_norm ;

  norm = 1 / sqrt(2 * M_PI * SQR(cross_time_sigma)) ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        for (frame1 = 0 ; frame1 < mri->nframes ; frame1++)
        {
          val0 = MRIgetVoxVal(mri, x, y, z, frame1) ;
          for (total = total_norm = 0.0, frame2 = 0 ; frame2 < mri->nframes ; frame2++)
          {
            val = MRIgetVoxVal(mri, x, y, z, frame2) ;
            g = norm * exp( - SQR(val-val0) / (2 * SQR(cross_time_sigma))) ;
            total += g*val ; 
            total_norm += g ;
          }
          total /= total_norm ;
          MRIsetVoxVal(mri, x, y, z, frame1, total) ;
        }
      }


  return(NO_ERROR) ;
}


