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
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "label.h"
#include "mri_identify.h"
#include "colortab.h"

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static int segment_rois = 1 ;
static  float  thresh =  0.99 ;

int
main(int argc, char *argv[]) {
  char               **av, fname[STRLEN], *surf_name, *wfile_name ;
  int                ac, nargs, vno ;
  MRI_SURFACE        *mris ;
  VERTEX             *v ;

  nargs = handleVersionOption(argc, argv, "mris_find_flat_regions");
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
    print_help() ;

  surf_name = argv[1] ;
  wfile_name = argv[2] ;

  mris = MRISread(surf_name) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  MRIScomputeMetricProperties(mris) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if ((fabs(v->nx)  > thresh) || (fabs(v->ny)  > thresh) || (fabs(v->nz)  > thresh))
      v->val = 1 ;
  }

  if (segment_rois)
  {
    int   nlabels, lno, num, vno ;
    LABEL    **labels ;

    MRISthresholdValIntoMarked(mris, thresh) ;
    MRISsegmentMarked(mris, &labels, &nlabels, 1) ;
    printf("%d labels total\n", nlabels) ;
    for (num = lno = 0 ; lno < nlabels ; lno++)
    {
      if (labels[lno]->n_points >= segment_rois)
      {
	if (num == 0)
	  printf("vertex %d gets annotation %d\n",  labels[lno]->lv[0].vno, num+1);
	for (vno = 0 ; vno < labels[lno]->n_points ; vno++)
	{
	  mris->vertices[labels[lno]->lv[vno].vno].annotation = num+1 ;
	}

	char fname[STRLEN] ;
	sprintf(fname, "flat%3.3d.label", num) ;
	num++ ;
	//LabelWrite(labels[lno], fname) ;
      }
    }
    printf("%d labels found. Writing annotation to %s\n", num, wfile_name) ;
    mris->ct = CTABalloc(num+1) ;
    CTABrandom(mris->ct) ;
    MRISwriteAnnotation(mris, wfile_name) ;
  }
  else if (mri_identify(wfile_name) == MGH_LABEL_FILE)
  {
    LABEL *area ;

    area = LabelFromSurface(mris, VERTEX_VALS, .9) ;
    LabelWrite(area, wfile_name) ;
  }
  else
    MRISwriteValues(mris, wfile_name) ;

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
  else switch (toupper(*option)) {
    case 'S':
      segment_rois  = atoi(argv[2]) ;
      nargs = 1 ;
      printf("segmenting surface into ROIs with at least %d vertices\n", segment_rois) ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'T':
      thresh = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "using threshold = %2.3f\n", thresh) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <surface> <wfile>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program computed regions in  which the surface is almost perpindicular to one\n"
          "of  the cardinal axes, and writes  the results  to  a label file\n") ;
  fprintf(stderr, "-t <thresh>  "
          "specify the threshold to use  (default=%2.3f)\n", thresh)  ;
  exit(1) ;
}


static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;

}
