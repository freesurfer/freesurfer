/*
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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "transform.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;
static char sdir[STRLEN] = "" ;

#define MAX_SUBJECTS 100

static float compute_overlap(MRI *mri_seg1, MRI *mri_seg2, TRANSFORM *transform1,
                             TRANSFORM *transform2) ;
int
main(int argc, char *argv[]) {
  char         **av, *xform_name, *out_fname, fname[STRLEN], *seg_name, *s1, *s2 ;
  int          ac, nargs, i, nsubjects, j, nvoxels ;
  MRI          *mri_seg[MAX_SUBJECTS] ;
  float        overlap, total_overlap ;
  TRANSFORM    *transform1, *transform2 ;
  FILE         *fp ;

  nargs = handleVersionOption(argc, argv, "mri_evaluate_morph");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (strlen(sdir) == 0) {
    char *cp ;
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in envoronment.\n",Progname);
    strcpy(sdir, cp) ;
  }
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit() ;

  xform_name = argv[1] ;
  seg_name = argv[2] ;
  out_fname = argv[argc-1] ;

#define FIRST_SUBJECT 3
  nsubjects = argc-(FIRST_SUBJECT+1) ;
  printf("processing %d subjects...\n", nsubjects) ;

  for (i = FIRST_SUBJECT ; i < argc-1 ; i++) {
    fprintf(stderr, "processing subject %s...\n", argv[i]) ;
    sprintf(fname, "%s/%s/mri/%s", sdir, argv[i], seg_name) ;
    mri_seg[i-FIRST_SUBJECT] = MRIread(fname) ;
    if (!mri_seg[i-FIRST_SUBJECT])
      ErrorExit(ERROR_NOFILE, "%s: could not read segmentation %s",
                Progname, fname) ;
  }

  fp = fopen(out_fname, "w") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not open output file %s...\n", out_fname) ;

  nvoxels = mri_seg[0]->width * mri_seg[0]->height * mri_seg[0]->depth ;
  for (total_overlap = 0.0f, i = 0 ; i < nsubjects ; i++) {
    for (j = i+1 ; j < nsubjects ; j++) {
      s1 = argv[i+FIRST_SUBJECT] ;
      s2 = argv[j+FIRST_SUBJECT] ;
      printf("reading transforms for subjects %s and %s...\n", s1, s2) ;
      sprintf(fname, "%s/%s/mri/transforms/%s", sdir, s1, xform_name) ;
      transform1 = TransformRead(fname) ;
      if (transform1 == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform %s",
                  Progname, fname) ;
      sprintf(fname, "%s/%s/mri/transforms/%s", sdir, s1, xform_name) ;
      transform2 = TransformRead(fname) ;
      if (transform2 == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform %s",
                  Progname, fname) ;
      printf("computing overlap for subjects %s and %s...\n", s1, s2) ;
      overlap = compute_overlap(mri_seg[i], mri_seg[j], transform1, transform2) ;
      total_overlap += overlap ;
      printf("overlap = %2.0f, total = %2.0f\n", overlap, total_overlap) ;
      fprintf(fp, "%s %s %2.0f %2.1f\n", s1, s2, overlap, 100.0f*overlap/(float)nvoxels) ;
      fflush(fp) ;
      TransformFree(&transform1) ;
      TransformFree(&transform2) ;
    }
  }

  total_overlap /= (float)((nsubjects*(nsubjects-1))/2.0f) ;
  printf("overlap/subject pair = %2.0f (%2.1f %%)\n", total_overlap,
         100.0f*total_overlap/(float)nvoxels) ;

  fclose(fp) ;
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
  else if (!stricmp(option, "sdir")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
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
          "usage: %s [options] <xform name> <s1> <s2> ... <output file>\n", Progname);
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will compute the overlap of a set of segmentations for a given morph\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static float
compute_overlap(MRI *mri_seg1, MRI *mri_seg2, TRANSFORM *transform1,
                TRANSFORM *transform2) {
  int    x, y, z, width, height, depth, l1, l2, x2, y2, z2 ;
  float  overlap, x1, y1, z1 ;

  TransformInvert(transform1, mri_seg1) ;

  width = mri_seg1->width ;
  height = mri_seg1->height ;
  depth = mri_seg1->depth ;
  for (overlap = 0.0f, x = 0 ; x < width ; x++) {
    for (y = 0 ; y < height ; y++) {
      for (z = 0 ; z < depth ; z++) {
        l1 = MRIvox(mri_seg1, x, y, z) ;
        TransformSample(transform1, x, y, z, &x1, &y1, &z1) ; /* atlas coords of s1 */
        TransformSampleInverseVoxel(transform2, width, height, depth,
                                    nint(x1), nint(y1), nint(z1), &x2, &y2, &z2) ;
        l2 = MRIvox(mri_seg2, x2, y2, z2) ;
        if (l1 == l2)
          overlap++ ;
      }
    }
  }
  return(overlap) ;
}

