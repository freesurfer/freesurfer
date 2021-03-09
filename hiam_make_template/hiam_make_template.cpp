/**
 * @brief hippocampus and amygdala version of mris_make_templace
 *
 * This code is a modified version of mris_make_template.c for being 
 * applied on hippocampus and amygdala.
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

#include "version.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static int which_norm = NORM_MEAN;

static const char *surface_names[] = {"hippocampus"};
static const char *curvature_names[] = {"hippocampus.curv"};

#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])
#define PARAM_IMAGES         (IMAGES_PER_SURFACE*SURFACES)

static int nbrs = 1 ;
static int navgs = 0 ;
static float scale = 1 ;
static int no_rot = 0 ;
static char subjects_dir[STRLEN] ;

int
main(int argc, char *argv[]) {
  char         **av, surf_fname[STRLEN], *template_fname, *hemi, *sphere_name,
  *cp, *subject, fname[STRLEN] ;
  int          ac, nargs, ino, sno ;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp, *mrisp_aligned, *mrisp_template ;
  INTEGRATION_PARMS parms ;

  memset(&parms, 0, sizeof(parms)) ;
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

  if (argc < 5)
    usage_exit() ;

  if (!strlen(subjects_dir))  /* not specified on command line*/
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment.\n",
                Progname) ;
    strcpy(subjects_dir, cp) ;
  }
  hemi = argv[1] ;
  sphere_name = argv[2] ;
  template_fname = argv[argc-1] ;
  if (1 || !FileExists(template_fname))  /* first time - create it */
  {
    fprintf(stderr, "creating new parameterization...\n") ;
    mrisp_template = MRISPalloc(scale, PARAM_IMAGES);
    if (no_rot)  /* don't do rigid alignment */
      mrisp_aligned = NULL ;
    else
      mrisp_aligned = MRISPalloc(scale, PARAM_IMAGES);
  } else {
    fprintf(stderr, "reading template parameterization from %s...\n",
            template_fname) ;
    mrisp_aligned = NULL ;
    mrisp_template = MRISPread(template_fname) ;
    if (!mrisp_template)
      ErrorExit(ERROR_NOFILE, "%s: could not open template file %s",
                Progname, template_fname) ;
  }

  argv += 3 ;
  argc -= 3 ;
  for (ino = 0 ; ino < argc-1 ; ino++) {
    subject = argv[ino] ;
    fprintf(stderr, "processing subject %s\n", subject) ;
    sprintf(surf_fname, "%s/%s/surf/%s.%s",
            subjects_dir, subject, hemi, sphere_name) ;
    fprintf(stderr, "reading spherical surface %s...\n", surf_fname) ;
    mris = MRISread(surf_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, surf_fname) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    MRISstoreMetricProperties(mris) ;

    for (sno = 0; sno < SURFACES ; sno++) {
      if (curvature_names[sno])  /* read in precomputed curvature file */
      {
        sprintf(surf_fname, "%s/%s/surf/%s.%s",
                subjects_dir, subject, hemi, curvature_names[sno]) ;
        if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
          ErrorExit(Gerror, "%s: could not read curvature file '%s'\n",
                    Progname, surf_fname) ;
      } else                       /* compute curvature of surface */
      {
        sprintf(surf_fname, "%s/%s/surf/%s.%s",
                subjects_dir, subject, hemi, surface_names[sno]) ;
        if (MRISreadVertexPositions(mris, surf_fname) != NO_ERROR)
          ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                    Progname, surf_fname) ;

        if (nbrs > 1)
          MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
        MRIScomputeMetricProperties(mris) ;
        MRIScomputeSecondFundamentalForm(mris) ;
        MRISuseMeanCurvature(mris) ;
        MRISaverageCurvatures(mris, navgs) ;
        MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
        MRISnormalizeCurvature(mris, which_norm) ;
      }
      fprintf(stderr, "computing parameterization for surface %s...\n",
              surf_fname);
      mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      MRISPcombine(mrisp, mrisp_template, sno*3) ;
      MRISPfree(&mrisp) ;
    }
    MRISfree(&mris) ;
  }

  if (mrisp_aligned)  /* new parameterization - use rigid alignment */
  {
    MRI_SP *mrisp_tmp ;

    if (Gdiag & DIAG_WRITE) {
      char *cp1 ;

      FileNameOnly(template_fname, fname) ;
      cp = strchr(fname, '.') ;
      if (cp) {
        cp1 = strrchr(fname, '.') ;
        if (cp1 && cp1 != cp)
          strncpy(parms.base_name, cp+1, cp1-cp-1) ;
        else
          strcpy(parms.base_name, cp+1) ;
      } else
        strcpy(parms.base_name, "template") ;
      sprintf(fname, "%s.%s.out", hemi, parms.base_name);
      INTEGRATION_PARMS_openFp(&parms, fname, "w") ;
      printf("writing output to '%s'\n", fname) ;
    }
    for (ino = 0 ; ino < argc-1 ; ino++) {
      subject = argv[ino] ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms.fp, "processing subject %s\n", subject) ;
      fprintf(stderr, "processing subject %s\n", subject) ;
      sprintf(surf_fname, "%s/%s/surf/%s.%s",
              subjects_dir, subject, hemi, sphere_name) ;
      fprintf(stderr, "reading spherical surface %s...\n", surf_fname) ;
      mris = MRISread(surf_fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, surf_fname) ;
      MRIScomputeMetricProperties(mris) ;
      MRISstoreMetricProperties(mris) ;
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
      sprintf(surf_fname, "%s/%s/surf/%s.%s",
              subjects_dir, subject, hemi, "hippocampus.curv") ;
      if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read curvature file '%s'\n",
                  Progname, surf_fname) ;
      parms.frame_no = 3 ;
      parms.mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      parms.mrisp_template = mrisp_template ;
      parms.l_corr = 1.0f ;

      MRISrigidBodyAlignGlobal(mris, &parms, 4.0, 32.0, 8) ;
      MRISPfree(&parms.mrisp) ;

#if 0
      /* write out rotated surface */
      sprintf(surf_fname, "%s.rot", mris->fname) ;
      fprintf(stderr, "writing out rigidly aligned surface to '%s'\n",
              surf_fname) ;
      MRISwrite(mris, surf_fname) ;
#endif

      /* now generate new parameterization using the optimal alignment */
      for (sno = 0; sno < SURFACES ; sno++) {
        if (curvature_names[sno])  /* read in precomputed curvature file */
        {
          sprintf(surf_fname, "%s/%s/surf/%s.%s",
                  subjects_dir, subject, hemi, curvature_names[sno]) ;
          if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
            ErrorExit(Gerror, "%s: could not read curvature file '%s'\n",
                      Progname, surf_fname) ;
        } else                       /* compute curvature of surface */
        {
          sprintf(surf_fname, "%s/%s/surf/%s.%s",
                  subjects_dir, subject, hemi, surface_names[sno]) ;
          if (MRISreadVertexPositions(mris, surf_fname) != NO_ERROR)
            ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                      Progname, surf_fname) ;

          if (nbrs > 1)
            MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
          MRIScomputeMetricProperties(mris) ;
          MRIScomputeSecondFundamentalForm(mris) ;
          MRISuseMeanCurvature(mris) ;
          MRISaverageCurvatures(mris, navgs) ;
          MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
          MRISnormalizeCurvature(mris, which_norm) ;
        }
        fprintf(stderr, "computing parameterization for surface %s...\n",
                surf_fname);
        mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
        MRISPcombine(mrisp, mrisp_aligned, sno*3) ;
        MRISPfree(&mrisp) ;
      }
      MRISfree(&mris) ;
    }

    if (Gdiag & DIAG_WRITE)
      fclose(parms.fp) ;

    mrisp_tmp = mrisp_aligned ;
    mrisp_aligned = mrisp_template ;
    mrisp_template = mrisp_tmp ;
    MRISPfree(&mrisp_aligned) ;
  }
  fprintf(stderr, "writing updated template to %s...\n", template_fname) ;
  MRISPwrite(mrisp_template, template_fname) ;
  MRISPfree(&mrisp_template) ;
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
  else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size = %d\n", nbrs) ;
  } else if (!stricmp(option, "sdir")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using SUBJECTS_DIR=%s\n", subjects_dir) ;
  } else if (!stricmp(option, "norot")) {
    no_rot = 1 ;
    fprintf(stderr, "not aligning hemispheres before averaging.\n") ;
  } else switch (toupper(*option)) {
    case 'W':
      Gdiag |= DIAG_WRITE ;
      if (isdigit((int)*argv[2]))
        nargs = 1 ;
      break ;
    case 'S':
      scale = atof(argv[2]) ;
      fprintf(stderr, "scaling parameterization by %2.1f\n", scale) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      fprintf(stderr, "averaging curvature patterns %d times.\n", navgs) ;
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
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <hemi> <surface name> <subject> <subject> ... "
          "<output name>\n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}
