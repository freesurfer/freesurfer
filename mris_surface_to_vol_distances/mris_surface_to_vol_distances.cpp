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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "mrisurf.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "matrix.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int update_histograms(MRI_SURFACE *mris, MRI_SURFACE *mris_avg, float ***histograms, int nbins) ;

const char *Progname ;
static void usage_exit(int code) ;
static char subjects_dir[STRLEN] = "";

static double min_distance = 1 ;
static double max_distance = 20 ;

#define MAX_SURFACE_SCALE 10

int
main(int argc, char *argv[]) {
  char         **av, fname[STRLEN] ;
  int          ac, nargs, i, j, nbins ;
  char         *avg_subject, *cp, *hemi, *subject, *output_prefix ;
  int          msec, minutes, seconds, nsubjects, vno ;
  Timer start ;
  MRI_SURFACE  *mris, *mris_avg ;
  float        ***histograms ;
  FILE         *fp ;

  nargs = handleVersionOption(argc, argv, "mris_surface_to_vol_distances");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  if (strlen(subjects_dir) == 0) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR must be specified on the command line (-sdir) or the env", Progname) ;
    strcpy(subjects_dir, cp) ;
  }

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  avg_subject = argv[1] ;
  hemi = argv[2] ;
  nsubjects = argc-4 ;
  output_prefix = argv[argc-1] ;
  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.sphere", subjects_dir, avg_subject, hemi) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris_avg = MRISread(fname) ;
  if (mris_avg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read spherical surface from %s", Progname, fname) ;

  nbins = nint(max_distance - min_distance) ;
  histograms = (float ***)calloc(mris_avg->nvertices, sizeof(float **)) ;
  if (histograms == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d histogram pointers", Progname, mris_avg->nvertices) ;
  for (vno = 0 ; vno < mris_avg->nvertices ; vno++) {
    histograms[vno] = (float **)calloc(nbins+1, sizeof(float *)) ;
    if (histograms[vno] == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate histogram bins %d", Progname, vno) ;
    for (i = 0 ; i < nbins ; i++) {
      histograms[vno][i] = (float *)calloc(MAX_SURFACE_SCALE*nbins+1, sizeof(float)) ;
      if (histograms[vno][i] == NULL)
        ErrorExit(ERROR_NOMEMORY, "%s: could not allocate histogram bins %d", Progname, vno) ;
    }
  }

  printf("processing %d subjects and writing results to %s*\n", nsubjects, output_prefix) ;

  for (i = 0 ; i < nsubjects ; i++) {
    subject = argv[i+3] ;
    printf("processing subject %s: %d of %d...\n", subject, i+1, nsubjects) ;
    int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.sphere.reg", subjects_dir, subject, hemi) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    mris = MRISread(fname) ;
    if (mris == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read spherical surface from %s", Progname, fname) ;
    req = snprintf(fname, STRLEN, "%s/%s/surf/%s.sphere", subjects_dir, subject, hemi) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (MRISreadCanonicalCoordinates(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read spherical surface from %s", Progname, fname) ;
    req = snprintf(fname, STRLEN, "%s/%s/surf/%s.white", subjects_dir, subject, hemi) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (MRISreadOriginalProperties(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read white surface from %s", Progname, fname) ;
    if (MRISreadCurvatureFile(mris, "thickness") != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read thickness file", Progname) ;
    mrisFindMiddleOfGray(mris) ;
    update_histograms(mris, mris_avg, histograms, nbins) ;
    MRISfree(&mris) ;
  }

  printf("writing log files with prefix %s...\n", output_prefix) ;
  for (vno = 0 ; vno < mris_avg->nvertices ; vno++) {

    int req = snprintf(fname, STRLEN, "%s%7.7d.histo", output_prefix, vno) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fp = fopen(fname, "w") ;
    for (i = 0 ; i < nbins ; i++) {
      for (j = 0 ; j < nbins*MAX_SURFACE_SCALE ; j++) {
        fprintf(fp, "%f ", histograms[vno][i][j]) ;
      }
      fprintf(fp, "\n") ;
    }
    fclose(fp) ;
  }

  /* compute average, store it in vertex 0 histgram, and write it out */
  for (vno = 1 ; vno < mris_avg->nvertices ; vno++) {
    for (i = 0 ; i < nbins ; i++)
      for (j = 0 ; j < nbins*MAX_SURFACE_SCALE ; j++)
        histograms[0][i][j] += histograms[vno][i][j] ;

  }

  sprintf(fname, "%s.average.histo", output_prefix) ;
  fp = fopen(fname, "w") ;

  for (i = 0 ; i < nbins ; i++) {
    for (j = 0 ; j < nbins*MAX_SURFACE_SCALE ; j++) {
      fprintf(fp, "%f ", histograms[0][i][j]) ;
    }
    fprintf(fp, "\n") ;
  }
  fclose(fp) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("distance histogram compilation took %d minutes"
         " and %d seconds.\n", minutes, seconds) ;
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
  if (stricmp(option, "sdir") == 0) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      printf("debugging vertex %d\n", Gdiag_no) ;
      break ;
    case 'D':
      max_distance = atof(argv[2]) ;
      nargs = 1 ;
      printf("setting maximum distance to %2.1f\n", max_distance) ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <average subject> <hemi> <subject 1> <subject 2> ... <output prefix>\n", Progname) ;
  exit(code) ;
}




/*
 white (or orig or pial) should be in orig vertices
 sphere.reg should be in current vertices
 sphere should be in canon vertices
*/

static int
update_histograms(MRI_SURFACE *mris, MRI_SURFACE *mris_avg, float ***histograms, int nbins) {
  int    vno, vno2, vno_avg ;
  double volume_dist, surface_dist, circumference, angle ;
  VERTEX *v1, *v2 ;
  VECTOR *vec1, *vec2 ;
  MHT    *mht ;
  float  **histogram, min_dist ;

  mht = MHTcreateVertexTable_Resolution(mris_avg, CURRENT_VERTICES, 2.0) ;

  vec1 = VectorAlloc(3, MATRIX_REAL) ;
  vec2 = VectorAlloc(3, MATRIX_REAL) ;

  v1 = &mris->vertices[0] ;
  VECTOR_LOAD(vec1, v1->cx, v1->cy, v1->cz) ;  /* radius vector */
  circumference = M_PI * 2.0 * V3_LEN(vec1) ;
  MRISclearMarks(mris_avg) ;
#if 0
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    if ((vno % 1000) ==  0) {
      printf("\r%d of %d    ", vno, mris->nvertices) ;
      fflush(stdout) ;
    }
    v1 = &mris->vertices[vno] ;
    VECTOR_LOAD(vec1, v1->cx, v1->cy, v1->cz) ;  /* radius vector */
    vno_avg = MHTfindClosestVertexNo(mht, mris_avg, v1, &min_dist) ;  /* which histogram to increment */
    if (vno_avg < 0)
      continue ;
    if (vno_avg == Gdiag_no)
      DiagBreak() ;
    histogram = histograms[vno_avg] ;
    mris_avg->vertices[vno_avg].marked = 1 ;

    for (vno2 = 0 ; vno2 < mris->nvertices ; vno2++) {
      if (vno2 == vno)
        continue ;
      v2 = &mris->vertices[vno2] ;
      VECTOR_LOAD(vec2, v2->cx, v2->cy, v2->cz) ;  /* radius vector */
      volume_dist = sqrt(SQR(v1->origx-v2->origx)+SQR(v1->origy-v2->origy)+SQR(v1->origz-v2->origz)) ;
      if (nint(volume_dist) >= nbins || nint(volume_dist) < 0)
        continue ;
      angle = fabs(Vector3Angle(vec1, vec2)) ;
      surface_dist = circumference * angle / (2.0 * M_PI) ;
      if (surface_dist > nbins*MAX_SURFACE_SCALE)
        surface_dist = nbins*MAX_SURFACE_SCALE ;
      if (surface_dist < 1)
        surface_dist = 1 ;

      histogram[nint(volume_dist)][nint(surface_dist)]++ ;

      if (mht->buckets[0][0] != NULL)
        DiagBreak() ;
    }
  }
  MHTfree(&mht) ;
#endif

  /* map back ones that were missed */
  /* printf("\nfilling holes in mapping\n") ;*/
  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 2.0) ;
  for (vno_avg = 0 ; vno_avg < mris_avg->nvertices ; vno_avg++) {
    if (mris_avg->vertices[vno_avg].marked > 0)
      continue ;

    if ((vno_avg % 1000) ==  0) {
      printf("\r%d of %d    ", vno_avg, mris_avg->nvertices) ;
      fflush(stdout) ;
    }

    vno = MHTfindClosestVertexNo2(mht, mris, mris_avg, &mris_avg->vertices[vno_avg], &min_dist) ;
    if (vno < 0)
      continue ;
    v1 = &mris->vertices[vno] ;
    VECTOR_LOAD(vec1, v1->cx, v1->cy, v1->cz) ;  /* radius vector */
    if (vno_avg < 0)
      continue ;
    if (vno_avg == Gdiag_no)
      DiagBreak() ;
    histogram = histograms[vno_avg] ;
    mris_avg->vertices[vno_avg].marked = 1 ;

    for (vno2 = 0 ; vno2 < mris->nvertices ; vno2++) {
      if (vno2 == vno)
        continue ;
      v2 = &mris->vertices[vno2] ;
      VECTOR_LOAD(vec2, v2->cx, v2->cy, v2->cz) ;  /* radius vector */
      volume_dist = sqrt(SQR(v1->origx-v2->origx)+SQR(v1->origy-v2->origy)+SQR(v1->origz-v2->origz)) ;
      if (nint(volume_dist) >= nbins || nint(volume_dist) < 0)
        continue ;
      angle = fabs(Vector3Angle(vec1, vec2)) ;
      surface_dist = circumference * angle / (2.0 * M_PI) ;
      if (surface_dist > nbins*MAX_SURFACE_SCALE)
        surface_dist = nbins*MAX_SURFACE_SCALE ;
      if (surface_dist < 1)
        surface_dist = 1 ;

      histogram[nint(volume_dist)][nint(surface_dist)]++ ;
    }
  }

  MHTfree(&mht) ;
  printf("\n") ;

  VectorFree(&vec1) ;
  VectorFree(&vec2) ;
  return(NO_ERROR) ;
}
