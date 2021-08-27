/**
 * @brief program for averaging "curvature" (scalar fields on the surface).
 *
 * average "curvature" format files (scalar fields on the surface) across subjects
 * using a surface-based registration (typically sphere.reg)
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


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static char sdir[STRLEN] = "" ;
static int which_norm = NORM_MEAN;
static int normalize_flag = 0 ;
static int condition_no = 0 ;
static int stat_flag = 0 ;
static char *output_surf_name = NULL ;
static int navgs = 0 ;
static char *ohemi = NULL ;
static char *osurf = NULL ;

int
main(int argc, char *argv[]) {
  char         **av, *in_fname, *out_fname, *surf_name, fname[STRLEN], *hemi ;
  int          ac, nargs, i, skipped ;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp, *mrisp_total ;

  nargs = handleVersionOption(argc, argv, "mris_average_curvature");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (!strlen(sdir)) 
  {
    char *cp ;
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    usage_exit() ;

  in_fname = argv[1] ;
  hemi = argv[2] ;
  if (!ohemi)
    ohemi = hemi ;
  surf_name = argv[3] ;
  if (!osurf)
    osurf = surf_name ;
  out_fname = argv[argc-1] ;

  mrisp_total = MRISPalloc(1, 3) ;
  skipped = 0 ;
  for (i = 4 ; i < argc-1 ; i++) {
    fprintf(stderr, "processing subject %s...\n", argv[i]) ;
    int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir, argv[i], hemi, surf_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    mris = MRISread(fname) ;
    if (!mris) {
      ErrorPrintf(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, fname) ;
      skipped++ ;
      continue ;
    }
    if (MRISreadCurvatureFile(mris, in_fname) != NO_ERROR) {
      skipped++ ;
      ErrorPrintf(ERROR_BADFILE,"%s: could not read curvature file %s.\n",
                  Progname, in_fname);
      continue ;
    }
    MRISaverageCurvatures(mris, navgs) ;
    if (normalize_flag)
      MRISnormalizeCurvature(mris, which_norm) ;
    mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
#if 0
    if (!FZERO(sigma)) {
      MRI_SP  *mrisp_blur ;

#if 0
      mrisp_blur = MRISPblur(mrisp, NULL, sigma, 0) ;
#else
      mrisp_blur = MRISPconvolveGaussian(mrisp, NULL, sigma,mris->radius,0);
#endif
      MRISPfree(&mrisp) ;
      mrisp = mrisp_blur ;
    }
#endif
    MRISPcombine(mrisp, mrisp_total, 0) ;
    MRISPfree(&mrisp) ;
    if (i < argc-2)
      MRISfree(&mris) ;
  }

  if (output_surf_name) {
    int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir,output_surf_name,ohemi,osurf);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "reading output surface %s...\n", fname) ;
    if (mris)
      MRISfree(&mris) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
  }
  if (mris==NULL) {
    printf("No output surface files found!\n");
    exit(1);
  }

  MRISfromParameterization(mrisp_total, mris, 0) ;
  if (stat_flag)    /* write out summary statistics files */
  {
    int    vno ;
    VERTEX *v ;
    float  dof ;
#if 1
    MRI    *mri ;

    mri = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 2) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      MRIsetVoxVal(mri, vno, 0, 0, 0, v->curv) ;
    }
    MRISfromParameterization(mrisp_total, mris, 1) ;
    /* change variances to squared standard errors */
    dof = *IMAGEFseq_pix(mrisp_total->Ip, 0, 0, 2) ;
    if (!FZERO(dof)) for (vno = 0 ; vno < mris->nvertices ; vno++) {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;
        v->curv /= dof ;   /* turn it into a standard error */
        MRIsetVoxVal(mri, vno, 0, 0, 1, v->curv) ;
      }
    mri->dof = dof ;
    MRIwrite(mri, out_fname) ;
    MRIfree(&mri) ;

#else
    FILE   *fp ;

    sprintf(fname, "%s/sigavg%d-%s.w", out_fname, condition_no, ohemi);
    fprintf(stderr, "writing output means to %s\n", fname) ;
    MRISwriteCurvatureToWFile(mris, fname) ;

    MRISfromParameterization(mrisp_total, mris, 1) ;

    /* change variances to squared standard errors */
    dof = *IMAGEFseq_pix(mrisp_total->Ip, 0, 0, 2) ;
    if (!FZERO(dof)) for (vno = 0 ; vno < mris->nvertices ; vno++) {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;
        v->curv /= dof ;   /* turn it into a standard error */
      }

    sprintf(fname, "%s/sigvar%d-%s.w", out_fname, condition_no, ohemi);
    fprintf(stderr, "writing output variances to %s\n", fname) ;
    MRISwriteCurvatureToWFile(mris, fname) ;
    /* write out dof file */
    sprintf(fname, "%s/sigavg%d.dof", out_fname, condition_no) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open dof file %s\n",
                Progname,fname);
    fprintf(stderr, "writing dof file %s\n", fname) ;
    fprintf(fp, "%d\n", (int)dof) ;
    fclose(fp) ;
#endif

  } else {
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr,"writing blurred pattern to surface to %s\n",out_fname);
    MRISwriteCurvature(mris, out_fname) ;
  }

  if (skipped > 0)
    printf("%d subjects skipped...\n", skipped) ;

  MRISfree(&mris) ;
  MRISPfree(&mrisp_total) ;
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
  else if (!stricmp(option, "ohemi")) {
    ohemi = argv[2] ;
    fprintf(stderr, "output hemisphere = %s\n", ohemi) ;
    nargs = 1 ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "osurf")) {
    osurf = argv[2] ;
    fprintf(stderr, "output surface = %s\n", osurf) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'A':
      navgs = atoi(argv[2]) ;
      fprintf(stderr, "blurring thickness for %d iterations\n",navgs);
      nargs = 1 ;
      break ;
    case 'O':
      output_surf_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, "painting output onto subject %s.\n", output_surf_name);
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'S':   /* write out stats */
      stat_flag = 1 ;
      condition_no = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "writing out summary statistics as condition %d\n",
              condition_no) ;
      break ;
    case 'N':
      normalize_flag = 1 ;
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
          "usage: %s [options] <input curv. file> <hemi> <surface> <subject> ... "
          " <output curv file >\n", Progname) ;
  fprintf(stderr, "the output curvature file will be painted onto the last "
          "subject name specified\non the command line.\n"
          "if the -s flag is specified then the last parameter specifies\n"
          "the directory in which to write the statistical maps.\n"
          "if the -o flag is specified then it overrides the last subject\n"
          "as the output surface.\n") ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-s <cond #>     generate summary statistics and write\n"
          "                them into sigavg<cond #>-<hemi>.w and\n"
          "                sigvar<cond #>-<hemi>.w.\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

