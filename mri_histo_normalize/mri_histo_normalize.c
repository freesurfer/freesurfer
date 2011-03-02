/**
 * @file  mri_histo_normalize.c
 * @brief do linear normalization of a set of histograms
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:19 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "mrinorm.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "transform.h"
#include "timer.h"
#include "version.h"
#include "histo.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;
static int   adaptive_normalize = 0 ;

static char *mask_fname = NULL ;
static char *xform_fname = NULL ;
static char sdir[STRLEN] = "" ;
static double tol = 0.5 ;  // avg rms change in vox intensities to terminate

#define MAX_SUBJECTS 100

int
main(int argc, char *argv[]) {
  char         **av, *out_vol_name, *in_vol_name, 
    *cp, fname[STRLEN], *subject;
  int          ac, nargs, n, done, niter ;
  int          msec, minutes, seconds, nsubjects ;
  struct timeb start ;
  MRI          *mri_inputs[MAX_SUBJECTS], *mri ;
  HISTOGRAM    *histos[MAX_SUBJECTS], *htemplate ;
  double       rms, rms_avg ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mri_histo_normalize.c,v 1.3 2011/03/02 00:04:19 nicks Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(
        ERROR_BADPARM,
        "%s: SUBJECTS_DIR not defined in environment.\n", 
        Progname) ;
    strcpy(sdir, cp) ;
  }

  in_vol_name = argv[1] ;
  out_vol_name = argv[argc-1] ;
  nsubjects = argc-3 ;
  printf(
    "argc=%d, out name %s, nsubjects = %d\n", 
    argc, out_vol_name, nsubjects);

  htemplate = HISTOalloc(256) ;
  htemplate->bin_size = 1 ; htemplate->min = 0 ; htemplate->max = 255 ;
  for (n = 0 ; n < 256 ; n++)
    htemplate->bins[n] = n ;

  for (n = 0 ; n < nsubjects ; n++)
  {
    subject = argv[n+2] ;
    printf("reading subject %s\n", subject) ;
    sprintf(fname, "%s/%s/mri/%s", sdir, subject, in_vol_name) ;
    mri_inputs[n] = MRIread(fname) ;
    if (mri_inputs[n] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
                Progname, fname) ;
    if (mask_fname) {
      MRI *mri_mask ;

      sprintf(fname, "%s/%s/mri/%s", sdir, subject, mask_fname) ;
      mri_mask = MRIread(fname) ;
      if (!mri_mask)
        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                  Progname, fname) ;
      // if mask == 0, then set dst as 0
      MRImask(mri_inputs[n], mri_mask, mri_inputs[n], 0, 0) ;
      MRIfree(&mri_mask) ;
    }
    histos[n] = MRIhistogram(mri_inputs[n], 0) ;
    HISTOclearZeroBin(histos[n]) ;
    HISTOadd(histos[n], htemplate, htemplate) ;
  }

  HISTOscalarMul(htemplate, 1.0/nsubjects, htemplate) ;
  niter = done = 0 ;
  if (Gdiag & DIAG_WRITE)
  {
    HISTOplot(histos[0], "h0.plt") ;
    HISTOplot(histos[1], "h1.plt") ;
    HISTOplot(htemplate, "h.plt") ;
  }
  do
  {
    float scale, offset ;
    for (rms_avg = 0.0, n = 0 ; n < nsubjects ; n++)
    {
      HISTOfindLinearFit(histos[n], htemplate, 
                         .5, 2, -20, 20, &scale, &offset) ;
      HISTOlinearScale(histos[n], histos[n], scale, offset) ;
      mri = MRIlinearScale(mri_inputs[n], NULL, scale, offset, 1) ;
      rms = MRIrmsDifferenceNonzero(mri, mri_inputs[n]) ;
      printf("histo %d: scaling by [%2.3f + %2.1f], rms = %2.2f\n", 
             n, scale, offset, rms) ;
      MRIfree(&mri_inputs[n]) ; mri_inputs[n] = mri ;
      rms_avg += rms ;
      HISTOfree(&histos[n]) ;
      histos[n] = MRIhistogram(mri_inputs[n], 256) ;
      HISTOclearZeroBin(histos[n]) ;
    }
    rms_avg /= nsubjects ;
    HISTOclear(htemplate, htemplate) ;
    for (n = 0 ; n < nsubjects ; n++)
      HISTOadd(histos[n], htemplate, htemplate) ;
    HISTOscalarMul(htemplate, 1.0/nsubjects, htemplate) ;

    if (Gdiag & DIAG_WRITE)
    {
      char *fname=NULL;
      HISTOplot(htemplate, "h.plt") ;
      for (n = 0 ; n < nsubjects ; n++)
      {
        sprintf(fname, "h%d.plt", n) ;
        HISTOplot(histos[n], fname) ;
        sprintf(fname, "m%d.%d.mgz", n, niter+1) ;
        MRIwrite(mri_inputs[n], fname) ;
      }
    }
    if (++niter > 10 || rms_avg < tol)
      done = 1 ;
  } while (!done) ;


  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  for (n = 0 ; n < nsubjects ; n++)
  {
    subject = argv[n+2] ;
    printf("reading subject %s\n", subject) ;
    sprintf(fname, "%s/%s/mri/%s", sdir, subject, out_vol_name) ;
    printf("writing output %s\n", fname) ;
    MRIwrite(mri_inputs[n], fname) ;
  }
  fprintf(stderr, 
          "histogram normalization took %d minutes and %d seconds.\n",
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
  if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "tol")) {
    tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting stopping tolerance to %2.3f rms intensity change\n",
           tol) ;
  } else if (!stricmp(option, "mask")) {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using %s as mask volume\n", mask_fname) ;
  } else switch (toupper(*option)) {
  case 'A':
    adaptive_normalize = 1 ;
    break ;
  case 'T':
    xform_fname = argv[2] ;
    nargs = 1 ;
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
  printf("usage: %s [options] <in vol name> <subject1 <subject2>... "
         "<out vol name>\n",
         Progname) ;
  printf(
    "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
  );
  exit(code) ;
}
