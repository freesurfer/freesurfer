/**
 * @brief program to compute ROC curves for surface segmentations
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "romp_support.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static char sdir[STRLEN] ;

static float min_area = 10 ;
static int dilate = 0 ;
static int erode = 0 ;

#define MAX_SUBJECTS 100
static int
compute_segmentation_stats(MRI_SURFACE *mris, LABEL *true_label, LABEL **segments, 
			   int nsegments,int *ptp, int *ptn, int *pfp, int *pfn);
static int  write_roc_curve(MRI_SURFACE *mris[MAX_SUBJECTS], LABEL *labels[MAX_SUBJECTS], 
			    MRI *mri_overlays[MAX_SUBJECTS], float min_area, int dilate, 
			    int erode, char *out_fname, int nsubjects) ;

int
main(int argc, char *argv[]) {
  char         **av, fname[STRLEN], *subject ;
  int          ac, nargs, i ;
  const char *hemi;
  char         *out_fname, *cp, *true_label_name, *segmentation_name ;
  int          msec, minutes, seconds, nsubjects ;
  Timer start ;
  LABEL        *labels[MAX_SUBJECTS] ;
  MRI          *mri_overlays[MAX_SUBJECTS] ;
  MRI_SURFACE  *mris[MAX_SUBJECTS] ;

  nargs = handleVersionOption(argc, argv, "mris_segmentation_stats");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

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

  if (argc < 3)
    usage_exit(1) ;

  if (strlen(sdir) == 0) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in env or cmd line",Progname) ;
    strcpy(sdir, cp) ;
  }
  segmentation_name = argv[1] ;
  true_label_name = argv[2] ;
  nsubjects = argc-4 ;
  for (i = 0 ; i < nsubjects ; i++)
  {
    subject = argv[i+3] ;
    sprintf(fname, "%s/%s/label/lh.%s.label", sdir, subject, true_label_name) ;
    if (FileExists(fname) == 0)
    {
      sprintf(fname, "%s/%s/label/rh.%s.label", sdir, subject, true_label_name) ;
    if (FileExists(fname) == 0)
      ErrorExit(ERROR_NOFILE, "%s: subject %s has no training label for either hemisphere", Progname, subject) ;
      hemi = "rh" ;
    }
    else
      hemi = "lh" ;
    printf("processing subject %s, hemi %s: %d of %d\n", subject, hemi,i+1, nsubjects) ;
    labels[i] = LabelRead(NULL, fname) ;
    if (labels[i] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load label from %s", Progname, fname) ;
    sprintf(fname, "%s/%s/surf/%s.white", sdir, subject, hemi) ;
    mris[i] = MRISread(fname) ;
    if (mris[i] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load surface from %s", Progname, fname) ;
    sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, segmentation_name) ;
    mri_overlays[i] = MRIread(fname) ;
    if (mri_overlays[i] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load overlay from %s", Progname, fname) ;
    MRIScopyMRI(mris[i], mri_overlays[i], 0, "val") ;
  }
  out_fname = argv[argc-1] ;

  write_roc_curve(mris, labels, mri_overlays, min_area, dilate, erode,out_fname,nsubjects);
  printf("writing outputs to %s\n", out_fname) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "segmentation statistics took %d minutes and %d seconds.\n",minutes,seconds);
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
  if (!stricmp(option, "sdir")) {
    strcpy(sdir, argv[2]) ;
    printf("using SUBJECTS_DIR=%s\n", sdir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "dilate")) {
    dilate = atof(argv[2]) ;
    printf("dilating segmentations %d times\n", dilate) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "erode")) {
    erode = atoi(argv[2]) ;
    printf("eroding segmentations %d times\n", erode) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "min_area")) {
    dilate = atof(argv[2]) ;
    printf("setting min segment area = %2.1f\n", min_area) ;
    nargs = 1 ;
  }
  else switch (toupper(*option)) {
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
  printf("usage: %s [options] <overlay name> <segmentation label name> <subject 1>... <roc file>",
         Progname) ;
  exit(code) ;
}

#define STEPS 1000


static int
compute_segmentation_stats(MRI_SURFACE *mris, LABEL *true_label, LABEL **segments, 
			   int nsegments, int *ptp, int *ptn, int *pfp, int *pfn)
{
  int    n, vno, tp, fp, tn, fn ;
  VERTEX *v ;

  MRISclearMarks(mris) ;
  LabelMarkSurface(true_label, mris) ;  // true label v->marked = 1 
  for (n = 0 ; n < nsegments ; n++)
    LabelAddToSurfaceMark(segments[n], mris, 2) ; 

  tp = fp = tn = fn = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    switch (v->marked)
    {
    default:
    case 0: tn++ ; break ;
    case 1: fn++ ; break ;
    case 2: fp++ ; break ;
    case 3: tp++ ; break ;
    }      
  }

  *ptp = tp ; *ptn = tn ; *pfp = fp ; *pfn = fn ;
  return(NO_ERROR) ;
}
static int
write_roc_curve(MRI_SURFACE *mris[MAX_SUBJECTS], LABEL *labels[MAX_SUBJECTS], 
		MRI *mri_overlays[MAX_SUBJECTS], float min_area, int dilate, 
		int erode, char *out_fname, int nsubjects)
{
  FILE   *fp ;
  float   min_val, max_val, total_min, total_max, step, thresh ;
  int     n, nlabels=0 ;
  LABEL   **segments=NULL ;
  int    tpos, fpos, tneg, fneg, tps[MAX_SUBJECTS], fps[MAX_SUBJECTS], tns[MAX_SUBJECTS], fns[MAX_SUBJECTS], i ;

  fp = fopen(out_fname, "w") ;

  total_min = 1e10 ; total_max = -total_min ;
  for (n = 0 ; n < nsubjects ; n++)
  {
    MRIvalRange(mri_overlays[n], &min_val, &max_val) ;
    if (min_val < total_min)
      total_min = min_val ;
    if (max_val > total_max)
      total_max = max_val ;
  }
  step = (total_max-total_min) / (STEPS-1) ;
  for (thresh = total_max+step ; thresh >= total_min-step ; thresh -= step)
  {
    printf("setting thresh = %2.4f\n", thresh) ;
    fpos = fneg = tpos = tneg = 0 ;
    i = 0 ;
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) firstprivate(i, thresh, dilate, erode, nlabels, segments, min_area) shared(mris, labels) schedule(static,1)
#endif
    for (n = 0 ; n < nsubjects ; n++)
    {
      ROMP_PFLB_begin
      MRISthresholdValIntoMarked(mris[n], thresh) ;
      MRISdilateMarked(mris[n], dilate) ;
      MRISerodeMarked(mris[n], erode) ;
      MRISsegmentMarked(mris[n], &segments, &nlabels, min_area) ;
      compute_segmentation_stats(mris[n], labels[n], segments, nlabels, 
				 &tps[n], &tns[n], &fps[n], &fns[n]);
      for (i = 0 ; i < nlabels ; i++)
	LabelFree(&segments[i]) ;
      ROMP_PFLB_end
    }
    ROMP_PF_end
    for (n = 0 ; n < nsubjects ; n++)
    {
      fpos += fps[n] ; fneg += fns[n] ; tpos += tps[n] ; tneg += tns[n] ;
    }
    fprintf(fp, "%f %d %d %d %d\n", thresh, tpos, tneg, fpos, fneg) ;
    fflush(fp) ;
  }


  fclose(fp) ;
  return(NO_ERROR) ;
}

