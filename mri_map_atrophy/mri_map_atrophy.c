/**
 * @file  mri_map_atrophy.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:22 $
 *    $Revision: 1.4 $
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
#include "transform.h"
#include "version.h"
#include "cma.h"
#include "gcamorph.h"


static char vcid[] = "$Id: mri_map_atrophy.c,v 1.4 2011/03/02 00:04:22 nicks Exp $";


static MRI *make_atrophy_map(MRI *mri_time1, MRI *mri_time2, MRI *mri_dst, TRANSFORM *transform1, TRANSFORM *transform2,
                             int *gray_labels, int ngray, int *csf_labels, int ncsf);

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static char *out_like_fname = NULL ;
static int invert_flag = 0 ;

static int gray_labels[] = {
                             Left_Hippocampus,
                             Left_Amygdala,
                             Left_Caudate,
                             Right_Hippocampus,
                             Right_Amygdala,
                             Right_Caudate
                           } ;
static int csf_labels[] = {
                            Left_Lateral_Ventricle,
                            Left_Inf_Lat_Vent,
                            Right_Lateral_Ventricle,
                            Right_Inf_Lat_Vent,
                            Unknown
                          };
static int ncsf = (sizeof(csf_labels) / sizeof(csf_labels[0]));
static int ngray = (sizeof(gray_labels) / sizeof(gray_labels[0]));

int
main(int argc, char *argv[]) {
  char        **av, *out_vol ;
  int         ac, nargs ;
  MRI         *mri_time1, *mri_time2, *mri_tmp, *mri_atrophy ;
  TRANSFORM   *transform1, *transform2 ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_map_atrophy.c,v 1.4 2011/03/02 00:04:22 nicks Exp $", "$Name:  $");
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

  out_vol = argv[argc-1] ;

  printf("reading volume from %s...\n", argv[1]) ;
  mri_time1 = MRIread(argv[1]) ;
  if (!mri_time1)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, argv[2]) ;
  mri_time2 = MRIread(argv[2]) ;
  if (!mri_time2)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, argv[2]) ;

  transform1 = TransformRead(argv[3]) ;
  if (!transform1)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, argv[3]) ;

  transform2 = TransformRead(argv[4]) ;
  if (!transform2)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, argv[4]) ;

  mri_tmp = TransformApplyType(transform1, mri_time1, NULL, SAMPLE_NEAREST);
  MRIfree(&mri_time1) ;
  mri_time1 = mri_tmp ;
  mri_tmp = TransformApplyType(transform2, mri_time2, NULL, SAMPLE_NEAREST);
  MRIfree(&mri_time2) ;
  mri_time2 = mri_tmp ;
  mri_atrophy = make_atrophy_map(mri_time1, mri_time2, NULL, transform1, transform2, gray_labels, ngray, csf_labels, ncsf) ;

  MRIwrite(mri_atrophy, out_vol) ;

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
    out_like_fname = argv[2] ;
    nargs = 1 ;
    printf("shaping output to be like %s...\n", out_like_fname) ;
  } else switch (toupper(*option)) {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'I':
      invert_flag = 1 ;
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
          "usage: %s [options] <seg time1> <seg time 2> <transform 1> <transform 2> <output file>\n",Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will apply a transform to mri volume and write out the result.  The output volume is by default 256^3 1mm^3 isotropic, or you can specify an -out_like volume.  I think there's a bug in -i behavior if you're specifying multiple transforms.\n");
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


static MRI *
make_atrophy_map(MRI *mri_time1, MRI *mri_time2, MRI *mri_dst, TRANSFORM *transform1, TRANSFORM *transform2,
                 int *gray_labels, int ngray, int *csf_labels, int ncsf) {
  int            x, y, z, label1, label2, n, found, xp, yp, zp, spacing ;
  GCA_MORPH_NODE *gcamn1, *gcamn2 ;
  GCA_MORPH      *gcam1, *gcam2 ;
  float           volume ;

  if (mri_dst == NULL) {
    mri_dst = MRIalloc(mri_time1->width, mri_time1->height, mri_time1->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_time1, mri_dst) ;
  }

  gcam1 = (GCA_MORPH*)transform1->xform ;
  gcam2 = (GCA_MORPH*)transform2->xform ;
  spacing = gcam1->spacing ;

  for (x = 0 ; x < mri_time1->width ; x++) {
    xp = x / spacing;
    for (y = 0 ; y < mri_time1->height ; y++) {
      yp = y / spacing;
      for (z = 0 ; z < mri_time1->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label1 = MRIgetVoxVal(mri_time1, x, y, z, 0) ;
        label2 = MRIgetVoxVal(mri_time2, x, y, z, 0) ;
        if (label1 == label2)
          continue ;

        /* if label1 was one of the gray types and label2 one of the csf, call it atrophy */
        for (found = n = 0 ; n < ngray ; n++)
          if (label1 == gray_labels[n]) {
            found = 1 ;
            break ;
          }
        if (found == 0)
          continue ;
        for (found = n = 0 ; n < ncsf ; n++)
          if (label2 == csf_labels[n]) {
            found = 1 ;
            break ;
          }
        if (found == 0)
          continue ;
        zp = z / spacing;
        gcamn1 = &gcam1->nodes[xp][yp][zp] ;
        gcamn2 = &gcam2->nodes[xp][yp][zp] ;
        volume = 0 ;
        if (FZERO(gcamn1->area) == 0)
          volume += gcamn1->orig_area / gcamn1->area ;
        if (FZERO(gcamn2->area) == 0)
          volume += gcamn2->orig_area / gcamn2->area ;
        MRIsetVoxVal(mri_dst, x, y, z, 0, volume) ;
      }
    }
  }


  return(mri_dst) ;
}

