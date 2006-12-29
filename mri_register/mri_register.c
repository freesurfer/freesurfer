/**
 * @file  mri_register.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:08 $
 *    $Revision: 1.20 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
#include "transform.h"
#include "version.h"

char         *Progname ;
static MORPH_PARMS  parms ;

static int get_option(int argc, char *argv[]) ;
static LTA  *register_mri(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms) ;
static unsigned char thresh_low = 90 ;
static unsigned char thresh_hi = 120 ;
static int  MRIsetFrame(MRI *mri, int frame, float val) ;

static int linear = 0 ;
static char *transform_fname = NULL ;
static int unit_variance = 1 ;
static char *var_fname = NULL ;

/*
   command line consists of three inputs:

   argv[1]  - input brain (to be registered)
   argv[2]  - 'canonical' brain (target of registration)
   argv[3]  - transform fle (output)
*/

int
main(int argc, char *argv[]) {
  char         ref_fname[100], *in_fname, *out_fname, fname[100], **av ;
  MRI          *mri_ref, *mri_in ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MORPH_3D     *m3d ;
  MRI          *mri_in_reduced, *mri_ref_reduced ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_register.c,v 1.20 2006/12/29 02:09:08 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  parms.morph_skull = 0 ;
  parms.niterations = 100 ;
  parms.l_intensity = 1.0 ;
  parms.levels = -1 ;   /* use default */
  parms.l_dist = 1.0 ;
  parms.l_area = 0.05 ;
  parms.l_nlarea = 100 ;
  parms.sigma = 4 ;
  parms.dt = .25 ;
  parms.tol = 0.1 /*1.5e-2*/ ;
  parms.exp_k = 25.0 ;
  Progname = argv[0] ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <in brain> <ref brain> <output file name>\n",
              Progname) ;

  in_fname = argv[1] ;
  strcpy(ref_fname, argv[2]) ;
#if 0
  if (strchr(ref_fname, '#') == NULL)
    strcat(ref_fname, "#-2") ;   /* only read in 2 frames - I know, a hack */
#endif
  out_fname = argv[3] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  Gdiag |= DIAG_WRITE ;
  fprintf(stderr, "logging results to %s.log\n", parms.base_name) ;

  TimerStart(&start) ;
  /*
     check to make sure output directory exists, and, if not, create it.
     */
  fprintf(stderr, "reading '%s'...\n", ref_fname) ;
  fflush(stderr) ;
  mri_ref = MRIread(ref_fname) ;
  if (!mri_ref)
    ErrorExit(ERROR_NOFILE, "%s: could not open reference volume %s.\n",
              Progname, ref_fname) ;
  if (var_fname)  /* read in a volume of standard deviations */
  {
    MRI *mri_var, *mri_tmp ;

    fprintf(stderr, "reading '%s'...\n", var_fname) ;
    mri_var = MRIread(var_fname) ;
    if (!mri_var)
      ErrorExit(ERROR_NOFILE, "%s: could not open variance volume %s.\n",
                Progname, var_fname) ;
    mri_tmp = MRIconcatenateFrames(mri_ref, mri_var, NULL) ;
    MRIfree(&mri_var) ;
    MRIfree(&mri_ref) ;
    mri_ref = mri_tmp ;
  }
  if (unit_variance && mri_ref->nframes > 1)
    MRIsetFrame(mri_ref, 1, 1.0) ;
#if 0
  if (mri_ref->nframes > 2)
    MRIfreeFrames(mri_ref, 2) ;
#endif

  fprintf(stderr, "reading '%s'...\n", in_fname) ;
  fflush(stderr) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;

#if 1
  if (mri_in->width > 128) {
    mri_in_reduced = MRIreduceByte(mri_in, NULL) ;
    /*    MRIfree(&mri_in) ;*/
  } else
    mri_in_reduced = mri_in ;
  if (mri_ref->width > 128) {
    mri_ref_reduced = MRIreduceMeanAndStdByte(mri_ref, NULL) ;
    /*    MRIfree(&mri_ref) ;*/
  } else
    mri_ref_reduced = mri_ref ;

  if (linear || !parms.lta)    /* find optimal linear transformation */
  {
    MRI *mri_in_red, *mri_ref_red ;
    char fname[100], *cp ;

    mri_in_red = MRIcopy(mri_in_reduced, NULL) ;
    mri_ref_red = MRIcopy(mri_ref_reduced, NULL) ;
    fprintf(stderr, "computing optimal linear transformation...\n") ;
    parms.lta = register_mri(mri_in_red, mri_ref_red, &parms) ;
    strcpy(fname, out_fname) ;
    cp = strrchr(fname, '.') ;
    if (!cp)
      cp = fname+strlen(fname) ;
    strcpy(cp, ".lta") ;
    fprintf(stderr, "writing output transformation to %s...\n", fname) ;
    // add src and ref information to lta
    getVolGeom(mri_in_red, &parms.lta->xforms[0].src);
    getVolGeom(mri_ref_red, &parms.lta->xforms[0].dst);
    LTAwriteEx(parms.lta, fname) ;
    MRIfree(&mri_in_red) ;
    MRIfree(&mri_ref_red) ;
  }

  parms.mri_ref = mri_ref ;
  parms.mri_in = mri_in ;
  m3d = MRI3Dmorph(mri_in_reduced, mri_ref_reduced, &parms) ;
  MRIfree(&mri_ref_reduced) ;
  MRIfree(&mri_in_reduced) ;
#else
  m3d = MRI3Dmorph(mri_in, mri_ref, &parms) ;
#endif

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing 3d morph transform to %s...\n", out_fname) ;
  MRI3Dwrite(m3d, out_fname) ;
  LTAfree(&parms.lta) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "morphing took %d minutes and %d seconds.\n",
          minutes, seconds) ;
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
  StrUpper(option) ;
  if (!stricmp(option, "TOL")) {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "tol = %2.2e\n", parms.tol) ;
  } else if (!strcmp(option, "DIST") || !strcmp(option, "DISTANCE")) {
    parms.l_dist = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_dist = %2.2f\n", parms.l_dist) ;
  } else if (!strcmp(option, "COMP") || !strcmp(option, "COMPRESS") ||
             !strcmp(option, "COMPRESSION")) {
    parms.l_compression = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_compression = %2.2f\n", parms.l_compression) ;
  } else if (!strcmp(option, "DT")) {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt = %2.2e\n", parms.dt) ;
  } else if (!stricmp(option, "skull")) {
    parms.morph_skull = 1 ;
    fprintf(stderr, "morphing skulls into register before 3d morph\n") ;
  } else if (!strcmp(option, "NOVAR")) {
    unit_variance = 1 ;
    fprintf(stderr, "disabling use of variance in morph\n") ;
  } else if (!strcmp(option, "USEVAR")) {
    unit_variance = 0 ;
    fprintf(stderr, "enabling use of variance in morph\n") ;
  } else if (!strcmp(option, "AREA")) {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_area = %2.2f\n", parms.l_area) ;
  } else if (!strcmp(option, "NLAREA")) {
    parms.l_nlarea = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nlarea = %2.2f\n", parms.l_nlarea) ;
  } else if (!strcmp(option, "NLDIST")) {
    parms.l_nldist = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nldist = %2.2f\n", parms.l_nldist) ;
  } else if (!strcmp(option, "LEVELS")) {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "levels = %d\n", parms.levels) ;
  } else if (!strcmp(option, "INTENSITY") || !strcmp(option, "CORR")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.2f\n", parms.l_intensity) ;
  } else switch (*option) {
#if 0
    case 'S':
      parms.sigma = atof(argv[2]) ;
      fprintf(stderr, "using sigma=%2.3f as upper bound on blurring.\n",
              parms.sigma) ;
      nargs = 1 ;
      break ;
#else
    case 'S':
      parms.sigma = atof(argv[2]) ;
      fprintf(stderr, "blurring gradient using Gaussian with sigma = %2.2f\n",
              parms.sigma) ;
      nargs = 1 ;
      break ;
#endif
    case '?':
    case 'U':
      printf("usage: %s [image file name]\n", argv[0]) ;
      exit(1) ;
      break ;
#if 0
    case 'T':
      thresh_low = atoi(argv[2]) ;
      thresh_hi = atoi(argv[3]) ;
      fprintf(stderr, "thresholds set to %d --> %d\n", thresh_low, thresh_hi) ;
      nargs = 2 ;
      break ;
#endif
    case 'T': {
      transform_fname = argv[2] ;
      fprintf(stderr, "reading optimal linear alignment from %s...\n",
              transform_fname) ;
      nargs = 1 ;
      parms.lta = LTAread(transform_fname) ;
      if (!parms.lta)
        ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                  Progname, transform_fname) ;
    }
    break ;
    case 'A':
      parms.navgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "navgs = %d\n", parms.navgs) ;
      break ;
    case 'L':
      linear = 1 ;
      fprintf(stderr, "finding optimal linear alignment...\n") ;
      break ;
    case 'V':
      var_fname = argv[2] ;
      fprintf(stderr, "reading variance image from %s...\n", var_fname) ;
      nargs = 1 ;
      break ;
    case 'N':
      parms.niterations = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "niterations = %d\n", parms.niterations) ;
      break ;
    case 'K':
      parms.exp_k = atof(argv[2]) ;
      fprintf(stderr, "using exponential k = %2.1f\n", parms.exp_k) ;
      nargs = 1 ;
      break ;
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'M':
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static LTA *
register_mri(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms) {
  MORPH_PARMS lparms ;

  memset(&lparms, 0, sizeof(lparms)) ;
  lparms.write_iterations = parms->write_iterations ;
  lparms.mri_ref = mri_ref ;
  lparms.mri_in = mri_in ;
  lparms.dt = 5e-6 ;
  lparms.tol = 1e-4 ;
  lparms.niterations = 200 ;
  lparms.momentum = 0.8 ;
  lparms.l_intensity = 1.0f ;
  lparms.lta = LTAalloc(1, mri_in) ;
  if (MRIfindNeck(mri_ref, mri_ref, thresh_low, thresh_hi, &lparms, 1,
                  &lparms.ref_np) == NULL) {
    lparms.disable_neck = 1 ;
    ErrorPrintf(ERROR_BADPARM, "%s: registration failed", Progname) ;
  }
  if (MRIfindNeck(mri_in, mri_in, thresh_low, thresh_hi, &lparms,
                  -1,&lparms.in_np)  == NULL) {
    lparms.disable_neck = 1 ;
    ErrorPrintf(ERROR_BADPARM, "%s: registration failed", Progname) ;
  }

  MRIlinearAlign(mri_in, mri_ref, &lparms) ;

  return(lparms.lta) ;
}

static int
MRIsetFrame(MRI *mri, int frame, float val) {
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      psrc = &MRIseq_vox(mri, 0, y, z, frame) ;
      for (x = 0 ; x < width ; x++) {
        *psrc++ = val ;
      }
    }
  }
  return(NO_ERROR) ;
}
