/**
 * @file  mris_transform.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:34 $
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
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"
#include "transform.h"
#include "version.h"

static char vcid[] = "$Id: mris_transform.c,v 1.8 2011/03/02 00:04:34 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static MRI_SURFACE  *mris ;
static int invert = 0 ;

MRI          *mri = 0;
MRI          *mri_dst = 0;

int
main(int argc, char *argv[]) {
  char         **av, *in_fname, *out_fname, *xform_fname ;
  int          ac, nargs ;
  LTA          *lta=0 ;
  int          transform_type;
  TRANSFORM    *transform ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_transform.c,v 1.8 2011/03/02 00:04:34 nicks Exp $", "$Name:  $");
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

  if (argc < 3)
    usage_exit() ;

  in_fname = argv[1] ;
  xform_fname = argv[2] ;
  out_fname = argv[3] ;

  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  // read transform
  transform_type =  TransformFileNameType(xform_fname);
  transform = TransformRead(xform_fname) ;
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
              Progname, xform_fname) ;
  if (transform->type == MNI_TRANSFORM_TYPE ||
      transform->type == TRANSFORM_ARRAY_TYPE ||
      transform->type  == REGISTER_DAT) {
    lta = (LTA *)(transform->xform) ;

    if (mri == 0 && lta->xforms[0].src.valid == 0) {
      fprintf(stderr, "The transform does not have the valid src volume info.\n");
      fprintf(stderr, "Either you give src volume info by option --src or\n");
      fprintf(stderr, "make the transform to have the valid src info.\n");
      ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
    }
    if (mri_dst == 0 && lta->xforms[0].dst.valid == 0) {
      fprintf(stderr, "The transform does not have the valid dst volume info.\n");
      fprintf(stderr, "Either you give src volume info by option --dst or\n");
      fprintf(stderr, "make the transform to have the valid dst info.\n");
      fprintf(stderr, "If the dst was average_305, then you can set\n");
      fprintf(stderr, "environmental variable USE_AVERAGE305 true\n");
      fprintf(stderr, "without giving the dst volume for RAS-to-RAS transform.\n");
      ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
    }
    //
    // depend on the type of transform you have to handle differently
    //
    //       orig              ------>      RAS   c_(ras) != 0
    //        |                              |
    //        |                              | identity
    //        V                              V
    //    conformed vol        ------>      RAS
    //        |                              |
    //        | identity                     |
    //        V                              V
    //    conformed vol        ------>   surfaceRAS
    //
    // given a volume transform you have to create a surfaceRAS transform
    //
    // Note that vertices are given by surfaceRAS coordinates
    //
    // RAS-to-RAS transform
    //
    //    surfaceRAS--->RAS --(ras-to-ras)-->RAS -->surfaceRAS
    //
    // VOX-to-Vox transform
    //
    //    surfaceRAS--->Vox---(vox-to-vox)-->Vox -->surfaceRAS
    //
    //
    if (invert) {
      VOL_GEOM vgtmp;
      LT *lt;
      MATRIX *m_tmp = lta->xforms[0].m_L ;
      lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
      MatrixFree(&m_tmp) ;
      lt = &lta->xforms[0];
      if (lt->dst.valid == 0 || lt->src.valid == 0) {
        fprintf(stderr, "WARNING:***************************************************************\n");
        fprintf(stderr, "WARNING:dst volume infor is invalid.  Most likely produce wrong inverse.\n");
        fprintf(stderr, "WARNING:***************************************************************\n");
      }
      copyVolGeom(&lt->dst, &vgtmp);
      copyVolGeom(&lt->src, &lt->dst);
      copyVolGeom(&vgtmp, &lt->src);
    }
  } else {
    TransformInvert(transform, mri_dst) ;
    if (invert) {}
    //    ErrorExit(ERROR_BADPARM, "transform is not of MNI, nor Register.dat type");
  }
  //
  MRIStransform(mris, mri, transform, mri_dst) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing surface to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;

  if (mri)
    MRIfree(&mri);
  if (mri_dst)
    MRIfree(&mri_dst);

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
  else if (!stricmp(option, "-src")) {
    fprintf(stderr, "Reading src volume...\n");
    mri = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri) {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  } else if (!stricmp(option, "-dst")) {
    fprintf(stderr, "Reading dst volume...\n");
    mri_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_dst) {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  } else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "-invert"))
    invert = 1 ;
  else switch (toupper(*option)) {
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
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input surf> <transform file> <output surf>\n",
          Progname) ;
  fprintf(stderr, "  options: --src <volumename>  src volume\n");
  fprintf(stderr, "                 use this option if the transform is created by MNI mritotal\n");
  fprintf(stderr, "         : --dst <volumename>  dst volume\n");
  fprintf(stderr, "                 use this option if the transform target is <not> average_305\n");
  fprintf(stderr, "         : --invert    apply inverted transform\n");
  fprintf(stderr, "         : --help     print help\n");
  fprintf(stderr, "         : --version  print version\n");
  fprintf(stderr, "         : -?, -U     print usage.\n");
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will transform an MRI surface into Talairach space.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

