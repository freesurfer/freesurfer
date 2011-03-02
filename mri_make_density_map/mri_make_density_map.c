/**
 * @file  mri_make_density_map.c
 * @brief make a tissue density map from a segmentation
 *
 * apply a transform (optionally jacobian correcting it) to
 * a segmentation with partial volume estimates to build a
 * density map.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:22 $
 *    $Revision: 1.8 $
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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "gca.h"
#include "gcamorph.h"
#include "cma.h"
#include "transform.h"
#include "version.h"

#define UNIT_VOLUME 128

static char vcid[] = "$Id: mri_make_density_map.c,v 1.8 2011/03/02 00:04:22 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static char *out_like_fname = NULL ;

static int nreductions = 0 ;
static char *xform_fname = NULL ;
static LTA  *lta = NULL;
static TRANSFORM *transform = NULL ;
static float sigma = 0 ;
static int surf_flag = 0 ;
static GCA *gca ;
static double resolution = 0.25 ;

int
main(int argc, char *argv[]) {
  TRANSFORM   *transform ;
  char        **av, *seg_fname, *intensity_fname, *out_fname, *xform_fname ;
  int         ac, nargs, i, label ;
  MRI         *mri_seg, *mri_intensity, *mri_out = NULL, *mri_kernel, *mri_smoothed, 
              *mri_orig_intensity, *mri_target ;
  GCA_MORPH   *gcam ;


  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_make_density_map.c,v 1.8 2011/03/02 00:04:22 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    usage_exit() ;

  seg_fname = argv[1] ;
  intensity_fname = argv[2] ;
  xform_fname = argv[3] ;
  out_fname = argv[argc-1] ;

  if (surf_flag)
  {
#define PAD 4
    MRI_SURFACE *mris ;
    MRI         *mri_tmp ;
    mris = MRISread(seg_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface %s",
                Progname, seg_fname) ;
    mri_tmp = MRISfillInterior(mris, resolution, NULL) ;
    //    MRIreplaceValues(mri_tmp, mri_tmp, 1, target_label) ;
    mri_seg = MRIextractRegionAndPad(mri_tmp, NULL, NULL, PAD) ;
    MRIfree(&mri_tmp) ;
    
    MRISfree(&mris) ;
  }
  else
  {
    printf("reading segmentation volume from %s\n", seg_fname) ;
    mri_seg = MRIread(seg_fname) ;
    if (!mri_seg)
      ErrorExit(ERROR_NOFILE, "%s: could not read segmentation volume %s", Progname,
              seg_fname) ;
  }

  printf("reading intensity volume from %s\n", intensity_fname) ;
  mri_intensity = MRIread(intensity_fname) ;
  if (!mri_intensity)
    ErrorExit(ERROR_NOFILE, "%s: could not read intensity volume %s", Progname,
              intensity_fname) ;
  mri_orig_intensity = MRIcopy(mri_intensity, NULL) ;

#if 1
  if (!MRIgeometryMatched(mri_seg, mri_intensity)) // resample intensity to be like seg
  {
    MRI *mri_tmp ;

    mri_tmp = MRIresample(mri_intensity, mri_seg, SAMPLE_TRILINEAR) ;
    MRIfree(&mri_intensity) ;
    mri_intensity = mri_tmp ;
  }
#endif
  transform = TransformRead(xform_fname) ;
  if (transform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, xform_fname) ;

  if (transform->type == MORPH_3D_TYPE)
    gcam = (GCA_MORPH *)transform->xform ;
  else
    gcam = NULL ;
  for (i = 4 ; i < argc-1 ; i++) {
    label = atoi(argv[i]) ;
    if (gca) // make target volume from it
    {
      MRI        *mri_tmp, *mri_template ; 
      MRI_REGION box ;
      double     scale ;

      mri_tmp = MRIclone(mri_orig_intensity, NULL) ;
      mri_tmp->c_r = mri_tmp->c_a = mri_tmp->c_s = 0.0 ;
      MRIreInitCache(mri_tmp) ;
      GCAbuildMostLikelyVolumeForStructure(gca, mri_tmp, label, 0, NULL, NULL) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_tmp, "t.mgz") ;
      MRIboundingBox(mri_tmp, 1, &box) ;
      mri_target = MRIextractRegionAndPad(mri_tmp, NULL, &box, 10) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_target, "t.mgz") ;
      MRIfree(&mri_tmp) ;

      scale = mri_target->xsize / resolution ;
      mri_template = MRIalloc(scale*mri_target->width, scale*mri_target->height, 
                              scale*mri_target->depth, MRI_UCHAR) ;
      MRIcopyHeader(mri_target, mri_template) ;
      MRIsetResolution(mri_template, resolution, resolution, resolution) ;
      mri_tmp = MRIresample(mri_target, mri_template, SAMPLE_TRILINEAR) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_tmp, "thr.mgz") ;
      MRIfree(&mri_target) ;
      mri_target = MRIchangeType(mri_tmp, MRI_FLOAT, 0, 1, 1) ;
      MRIfree(&mri_tmp) ;
    }
    else
      mri_target = NULL ;

    if (surf_flag)
      MRIreplaceValues(mri_seg, mri_seg, 1, label) ;
    printf("extracting label %d (%s)\n", label, cma_label_to_name(label)) ;
    if (gcam)
      mri_out = GCAMextract_density_map(mri_seg, mri_intensity, gcam, label, mri_out) ;
    else
    {
      MRI *mri_tmp, *mri_cor ;
      MATRIX *m_L ;
      LTA    *lta ;
      
      lta = ((LTA *)(transform->xform));
      if (lta->type == LINEAR_VOX_TO_VOX)
        LTAvoxelToRasXform(lta, NULL, NULL) ;
      mri_cor = MRIalloc(mri_seg->width, mri_seg->height, mri_seg->depth,
                         MRI_FLOAT) ;
      LTArasToVoxelXform(lta, mri_seg, mri_cor) ;
      m_L = lta->xforms[0].m_L ;

      mri_tmp = MRImakeDensityMap(mri_seg, mri_intensity, label, NULL, 1.0) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_tmp, "t.mgz") ;
      if (mri_target)  // output has different resolution than gca
      {
        MATRIX *m_vox2vox, *m_tmp ;
        m_vox2vox = MRIgetVoxelToVoxelXform(mri_cor, mri_target) ;

        m_tmp = MatrixMultiply(m_vox2vox, m_L, NULL) ;
        MatrixCopy(m_tmp, m_L) ; MatrixFree(&m_tmp) ;

        mri_out = MRIclone(mri_target, NULL) ;
        MRIfree(&mri_target) ;
        MatrixFree(&m_vox2vox) ;
      }
      else
      {
        mri_out = MRIclone(mri_cor, NULL) ;
      }
      MRIlinearTransformInterp(mri_tmp, mri_out, m_L, SAMPLE_TRILINEAR);
      MRIfree(&mri_tmp) ; MRIfree(&mri_cor) ;
    }
    /*    extract_labeled_image(mri_in, transform, label, mri_out) ;*/
  }
  if (!FZERO(sigma)) {
    printf("smoothing extracted volume...\n") ;
    mri_kernel = MRIgaussian1d(sigma, 10*sigma) ;
    mri_smoothed = MRIconvolveGaussian(mri_out, NULL, mri_kernel) ;
    MRIfree(&mri_out) ;
    mri_out = mri_smoothed ;
  }
  /* removed for gcc3.3
   * vsprintf(out_fname, out_fname, (va_list) &label) ;
   */
  while (nreductions-- > 0) {
    MRI *mri_tmp ;

    mri_tmp = MRIreduce(mri_out, NULL) ;
    MRIfree(&mri_out) ;
    mri_out = mri_tmp ;
  }
  printf("writing output to %s\n", out_fname) ;
  MRIwrite(mri_out, out_fname) ;


  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "surf"))
  {
    surf_flag = 1 ;
    printf("assuming input is a surface\n") ;
  }
  else if (!stricmp(option, "resolution"))
  {
    resolution = atof(argv[2]) ;
    printf("setting output resolution to %2.2f for surface\n", resolution) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  } else if (!stricmp(option, "out_like") || !stricmp(option, "ol")) {
    out_like_fname = argv[2] ;
    nargs = 1 ;
    printf("shaping output to be like %s...\n", out_like_fname) ;
  } else switch (toupper(*option)) {
  case 'A':
    printf("reading atlas from %s\n", argv[2]) ;
    gca = GCAread(argv[2]) ;
    if (gca == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read atlas from %s", Progname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'R':
    nreductions = atoi(argv[2]) ;
    printf("reducing density maps %d times...\n", nreductions);
    nargs = 1 ;
    break ;
  case 'S':
    sigma = atof(argv[2]) ;
    printf("applying sigma=%2.1f smoothing kernel after extraction...\n",sigma) ;
    nargs = 1 ;
    break ;
  case 'T':
    xform_fname = argv[2] ;
    printf("reading and applying transform %s...\n", xform_fname) ;
    nargs = 1 ;
    transform = TransformRead(xform_fname) ;
    if (!transform)
      ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s",
                Progname, xform_fname) ;

    if (transform->type != MORPH_3D_TYPE)
      lta = (LTA *)(transform->xform) ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <aseg volume> <norm volume> <xform> <label 1> <label 2> ... <output name>\n",
          Progname) ;
  fprintf(stderr, "where options are:\n") ;
  fprintf(stderr,
          "\t-s <sigma>\tapply a Gaussian smoothing kernel\n"
          "\t-r <n>\tapply a Gaussian reduction n times\n") ;
  //          "\t-t <xform file>\tapply the transform in <xform file> to extracted volume\n");
}

static void
print_help(void) {
  fprintf(stderr,
          "\nThis program will extract a set of labeled voxels from an image\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


