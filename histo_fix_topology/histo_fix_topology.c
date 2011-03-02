/**
 * @file  histo_fix_topology.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.3 $
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

#include "numerics.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "version.h"
#include "matrix.h"
#include "mrisegment.h"

#define RGB_SIZE 500

static char vcid[] =
  "$Id: histo_fix_topology.c,v 1.3 2011/03/02 00:04:09 nicks Exp $";

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static double thresh = 50000 ;
static double max_thresh = 64000;

static int open = 5 ;

int
main(int argc, char *argv[]) {
  char        **av ;
  int         ac, nargs, max_i, i, j ;
  MRI         *mri_in, *mri_out ;
  MRI_SEGMENTATION  *mriseg ;
  MRI_SEGMENT       *ms ;
  MRI_SEGMENT_VOXEL *msv ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: histo_fix_topology.c,v 1.3 2011/03/02 00:04:09 nicks Exp $", "$Name:  $");
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

  if (argc < 2)
    usage_exit() ;

  mri_in = MRIread(argv[1]) ;
  if (mri_in == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not open input image %s...\n", argv[1]) ;

  mriseg = MRIsegment(mri_in, thresh, max_thresh) ;
  max_i = MRIfindMaxSegmentNumber(mriseg) ;
  
  mri_out = MRIclone(mri_in, NULL) ;
  for (i = 0 ; i < mriseg->nsegments ; i++)
  {
      if (i == max_i)
        continue ;
      ms = &mriseg->segments[i] ;
      for (j = 0 ; j < ms->nvoxels ; j++)
      {
        msv = &ms->voxels[j] ;
        MRIsetVoxVal(mri_out, msv->x, msv->y, msv->z, 0, 1);
      }
  }
  MRIsegmentFree(&mriseg) ;

  if (open > 0 && 0)
  {
    MRI *mri_opened ;
    int i, x, y ;
    
    mri_opened = MRIcopy(mri_in, NULL) ;
    for (i = 0 ; i < open ; i++)
      MRIerode(mri_opened, mri_opened) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_opened, "e.mgz") ;
    for (i = 0 ; i < 2*open ; i++)
      MRIdilate(mri_opened, mri_opened) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_opened, "d.mgz") ;

    mri_out =  MRIclone(mri_in, NULL) ;
    for (x = 0 ; x < mri_in->width ; x++)
    {
      for (y = 0 ; y < mri_in->height ; y++)
      {
        if (x == Gx && y == Gy)
          DiagBreak() ;
        if ((MRIgetVoxVal(mri_in, x, y, 0, 0) > thresh) &&
            (MRIgetVoxVal(mri_opened, x, y, 0, 0) < thresh))
        {
          MRIsetVoxVal(mri_in, x, y, 0, 0, 0) ;
          MRIsetVoxVal(mri_out, x, y, 0, 0, 1) ;
        }
      }
    }
  }

  printf("writing corrected volume to %s\n", argv[2]) ;
  MRIwrite(mri_out, argv[2]) ;
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
  else if (!stricmp(option, "out_like") || !stricmp(option, "ol")) {
#if 0
    out_like_fname = argv[2] ;
    nargs = 1 ;
    printf("shaping output to be like %s...\n", out_like_fname) ;
#endif
  }
  else if (!stricmp(option, "debug_voxel") || !stricmp(option, "debug_pixel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    nargs = 2 ;
    printf("debugging pixel (%d, %d)\n", Gx, Gy) ;
  }
  else switch (toupper(*option))
  {
  case 'O':
    open = atof(argv[2]) ;
    nargs = 1 ;
    printf("opening image %d times to detect tears\n", open) ;
    break ;
  case 'T':
    thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("using threshold %2.1f for determining background\n", thresh) ;
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
  //  print_usage() ; // print_help _calls print_usage
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,"usage: %s [options] <input image> <output image>\n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program corrects for tears in an image\n") ;
  fprintf(stderr, "-out_like <reference volume> - set out_volume parameters\n") ;
  exit(1) ;
}


static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
static int
mriWriteImageView(MRI *mri, char *base_name, int target_size, int view,
                  int slice, MRI *mri_template) {
  char  fname[STRLEN], *prefix ;
  IMAGE *I ;
  float scale ;

  mri = MRIresampleFill(mri, mri_template, SAMPLE_NEAREST, 255) ;

  switch (view) {
  default:
  case MRI_CORONAL:
    prefix = "cor" ;
    break ;
  case MRI_SAGITTAL:
    prefix = "sag" ;
    break ;
  case MRI_HORIZONTAL:
    prefix = "hor" ;
    break ;
  }

  if (slice < 0)    /* pick middle slice */
  {
    switch (view) {
    default:
    case MRI_CORONAL:
      slice = mri->depth/2;
      break ;
    case MRI_SAGITTAL:
      slice = 6*mri->width/10 ;
      break ;
    case MRI_HORIZONTAL:
      slice = mri->height/2 ;
      break ;
    }
  }
  I = MRItoImageView(mri, NULL, slice, view, 0) ;
  if (!I)
    ErrorReturn(Gerror, (Gerror, "MRItoImageView failed")) ;

  scale = (float)target_size / (float)I->rows ;
  if (!FEQUAL(scale, 1.0f)) {
    IMAGE *Itmp ;

    Itmp = ImageRescale(I, NULL, scale) ;
    ImageFree(&I) ;
    I = Itmp ;
  }
  sprintf(fname, "%s_%s.rgb", prefix, base_name) ;
  ImageWrite(I, fname) ;
  ImageFree(&I) ;
  MRIfree(&mri) ;
  return(NO_ERROR) ;
}

static int
write_snapshot(MRI *mri, MRI *mri1, MRI *mri2, MATRIX *m, char *base, int n, int level, DENSITY *pdf, MRI *mri_seg) {
  MRI  *mri_xformed, *mri_ll ;
  char fname[STRLEN] ;


  mri = MRIresampleFill(mri, mri2, SAMPLE_NEAREST, 255) ;
  mri_xformed = mri_apply_slice_xform(mri, NULL, m, 0) ;

  sprintf(fname, "%s%3.3d", base, n) ;
  printf("writing snapshot to %s...\n", fname) ;
  mriWriteImageView(mri_xformed, fname, RGB_SIZE, MRI_CORONAL, -1, mri2) ;
  sprintf(fname, "%s%3.3d.mgz", base, n) ;
  MRIwrite(mri_xformed, fname) ;
  MRIfree(&mri_xformed) ;

  if (pdf) {
    MATRIX *m_inv ;

    mri_ll = DensityLikelihoodImage(mri1, mri2, NULL, m, pdf, mri_seg, 0) ;
    sprintf(fname, "%s_ll_%3.3d.mgz", base, n) ;
    printf("writing density map to %s...\n", fname) ;
    MRIwrite(mri_ll, fname) ;
    MRIfree(&mri_ll) ;

    m_inv = MatrixInverse(m, NULL) ;
    mri_ll = DensityLikelihoodImage(mri2, mri1, NULL, m_inv, pdf, mri_seg, 1) ;
    sprintf(fname, "%s_ll_src_%3.3d.mgz", base, n) ;
    printf("writing density map to %s...\n", fname) ;
    MRIwrite(mri_ll, fname) ;
    MRIfree(&mri_ll) ;
    MatrixFree(&m_inv) ;
  }
  return(NO_ERROR) ;
}

static MRI *
expand_tears(MRI *mri_histo, MRI *mri_src, MRI *mri_dst, int niter, double thresh)
{
  int  x, y, i, xi, yi, xk, yk ;
  MRI  *mri_tmp ;

  mri_dst = MRIcopy(mri_src, mri_dst);

  for (i = 0 ; i < niter ; i++)
  {
    for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
      {
        if (MRIgetVoxVal(mri_src, x, y, 0, 0) > 0) // in a tear
        {
          for (xk = -1 ; xk <= 1 ; xk++)
          {
            xi = mri_src->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_src->yi[y+yk] ;
            }
          }
        }
      }
    }
  }
  return(mri_dst) ;
}

#endif

