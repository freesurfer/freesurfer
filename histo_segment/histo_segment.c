/**
 * @file  histo_stereology.c
 * @brief program to segment neurons from a nissl stain
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.2 $
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
#include "transform.h"
#include "version.h"
#include "matrix.h"
#include "density.h"
#include "mrisegment.h"
#include "mri_circulars.h"

#define RGB_SIZE 500

static char vcid[] =
  "$Id: histo_segment.c,v 1.2 2011/03/02 00:04:09 nicks Exp $";

#if 0
static int write_snapshot(MRI *mri, MRI *mri_src, MRI *mri_dst, MATRIX *m, char *base, int n, int level, DENSITY *density, MRI *mri_seg) ;
int main(int argc, char *argv[]) ;
#endif
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

int
main(int argc, char *argv[]) {
  char        **av ;
  int         ac, nargs ;
  MRI         *mri_nissl, *mri_seg ;
  MRI_SEGMENTATION *mriseg ;
  double      thresh, min_neuron_area ;


  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: histo_segment.c,v 1.2 2011/03/02 00:04:09 nicks Exp $", "$Name:  $");
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

  if (argc < 4)
    usage_exit() ;

  mri_nissl = MRIread(argv[1]) ;
  if (mri_nissl == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not open Nissl image %s...\n", argv[1]) ;

  thresh = atof(argv[2]) ;
  min_neuron_area = atof(argv[3]) ;
  mriseg = MRIsegment(mri_nissl, 0, thresh) ;
  MRIremoveSmallSegments(mriseg, min_neuron_area) ;
  mri_seg = MRIsegmentToImage(mri_nissl, NULL, mriseg, -1) ;
  //  MRIbinarize(mri_seg, mri_seg, 0, 0, 1) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_seg, "s.mgz") ;

  MRIwrite(mri_seg, argv[4]) ;
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
  }else switch (toupper(*option)) {
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
  fprintf(stderr,
          "usage: %s [options] <Nissl image> <max intensity> <max area> <output segmented image>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program align a histological slice with a block face image\n") ;
  fprintf(stderr, "-out_like <reference volume> - set out_volume parameters\n") ;
  fprintf(stderr, "-I                           - invert transform "
          "coordinates\n") ;
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
#endif
