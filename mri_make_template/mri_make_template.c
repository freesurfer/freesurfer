/**
 * @file  mri_make_template.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:22 $
 *    $Revision: 1.26 $
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
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "transform.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;

// debug point
int DEBUG_X =62;
int DEBUG_Y =254;
int DEBUG_Z =92;

static float smooth = 0 ;
static int erode = 0 ;
static int open = 0 ;
static int binarize = 0 ;
static void usage_exit(int code) ;
static MRI *MRIcomputePriors(MRI *mri_priors, int ndof, MRI *mri_char_priors);
int MRIaccumulateMeansAndVariances(MRI *mri, MRI *mri_mean, MRI *mri_std) ;
int MRIcomputeMeansAndStds(MRI *mri_mean, MRI *mri_std, int ndof) ;
int MRIcomputeMaskedMeansAndStds(MRI *mri_mean, MRI *mri_std, MRI *mri_dof) ;
MRI *MRIfloatToChar(MRI *mri_src, MRI *mri_dst) ;
int MRIaccumulateMaskedMeansAndVariances(MRI *mri, MRI *mri_mask,MRI *mri_dof,
    float low_val, float hi_val,
    MRI *mri_mean, MRI *mri_std) ;
static MRI *MRIupdatePriors(MRI *mri_binary, MRI *mri_priors) ;
static MRI *MRIcomputePriors(MRI *mri_priors, int ndof, MRI *mri_char_priors);

static char *transform_fname = NULL ;
static char *T1_name = "T1" ;
static char *var_fname = NULL ;
static char *binary_name = NULL ;

/* just for T1 volume */
#define T1_MEAN_VOLUME        0
#define T1_STD_VOLUME         1

static int first_transform = 0 ;

#define BUILD_PRIORS             0
#define ON_STATS                 1
#define OFF_STATS                2


#define NPARMS  12
static MATRIX *m_xforms[NPARMS] ;
static MATRIX *m_xform_mean = NULL ;
static MATRIX *m_xform_covariance = NULL ;
static char *xform_mean_fname = NULL ;
static char *xform_covariance_fname = NULL ;
static int stats_only = 0 ;
static int novar = 0 ;

static char subjects_dir[STRLEN];

static int
check_mri(MRI *mri) {
  int error = 0, x, y, z ;
  float val ;

  for (x = 0 ; x < mri->width ; x++) {
    for (y = 0 ; y < mri->height ; y++) {
      for (z = 0 ; z < mri->depth ; z++) {
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        if (fabs(val) > 1e5 ) {
          error = 1 ;
          DiagBreak() ;
        }
      }
    }
  }

  return(error) ;
}

int
main(int argc, char *argv[]) {
  char   **av, *cp ;
  int    ac, nargs, i, dof, no_transform, which, sno = 0, nsubjects = 0 ;
  MRI    *mri=0, *mri_mean = NULL, *mri_std=0, *mri_T1=0,*mri_binary=0,*mri_dof=NULL,
                             *mri_priors = NULL ;
  char   *subject_name, *out_fname, fname[STRLEN] ;
  /*  LTA    *lta;*/
  MRI *mri_tmp=0 ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_make_template.c,v 1.26 2011/03/02 00:04:22 nicks Exp $", "$Name:  $");
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

  if (!strlen(subjects_dir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,"%s: SUBJECTS_DIR not defined in environment.\n",
                Progname) ;
    strcpy(subjects_dir, cp) ;
  }

  if (argc < 3)  usage_exit(1) ;

  out_fname = argv[argc-1] ;

  no_transform = first_transform ;
  if (binary_name)   /* generate binarized volume with priors and */
  {                  /* separate means and variances */
    for (which = BUILD_PRIORS ; which <= OFF_STATS ; which++) {
      /* for each subject specified on cmd line */
      for (dof = 0, i = 1 ; i < argc-1 ; i++) {
        if (*argv[i] == '-')   /* don't do transform for next subject */
        { no_transform = 1 ;
          continue ;
        }
        dof++ ;
        subject_name = argv[i] ;
        if (which != BUILD_PRIORS) {
          sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, T1_name);
          fprintf(stderr, "%d of %d: reading %s...\n", i, argc-2, fname) ;
          mri_T1 = MRIread(fname) ;
          if (!mri_T1)
            ErrorExit(ERROR_NOFILE,"%s: could not open volume %s",
                      Progname,fname);
        }

        sprintf(fname, "%s/%s/mri/%s",subjects_dir,subject_name,binary_name);
        fprintf(stderr, "%d of %d: reading %s...\n", i, argc-2, fname) ;
        mri_binary = MRIread(fname) ;
        if (!mri_binary)
          ErrorExit(ERROR_NOFILE,"%s: could not open volume %s",
                    Progname,fname);

        /* only count voxels which are mostly labeled */
        MRIbinarize(mri_binary, mri_binary, WM_MIN_VAL, 0, 100) ;
        if (transform_fname && no_transform-- <= 0) {
          sprintf(fname, "%s/%s/mri/transforms/%s",
                  subjects_dir, subject_name, transform_fname) ;

          fprintf(stderr, "reading transform %s...\n", fname) ;
          ////////////////////////////////////////////////////////
#if 1
          {
            TRANSFORM *transform ;
            transform = TransformRead(fname) ;
            if (transform == NULL)
              ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",Progname, fname) ;
            mri_tmp = TransformApply(transform, mri_T1, NULL) ;
            TransformFree(&transform) ;
          }
#else
          lta = LTAreadEx(fname);
          if (lta == NULL)
            ErrorExit(ERROR_NOFILE,
                      "%s: could not open transform file %s\n",
                      Progname, fname) ;
          /* LTAtransform() runs either MRIapplyRASlinearTransform()
          for RAS2RAS or MRIlinearTransform() for Vox2Vox. */
          /* MRIlinearTransform() calls MRIlinearTransformInterp() */
          mri_tmp = LTAtransform(mri_T1, NULL, lta);
          MRIfree(&mri_T1) ;
          mri_T1 = mri_tmp ;
          LTAfree(&lta);
          lta = NULL;
#endif
          if (DIAG_VERBOSE_ON)
            fprintf(stderr, "transform application complete.\n") ;
        }
        if (which == BUILD_PRIORS) {
          mri_priors =
            MRIupdatePriors(mri_binary, mri_priors) ;
        } else {
          if (!mri_mean) {
            mri_dof = MRIalloc(mri_T1->width, mri_T1->height, mri_T1->depth,
                               MRI_UCHAR) ;
            mri_mean =
              MRIalloc(mri_T1->width, mri_T1->height,mri_T1->depth,MRI_FLOAT);
            mri_std =
              MRIalloc(mri_T1->width,mri_T1->height,mri_T1->depth,MRI_FLOAT);
            if (!mri_mean || !mri_std)
              ErrorExit(ERROR_NOMEMORY, "%s: could not allocate templates.\n",
                        Progname) ;
          }

          if (DIAG_VERBOSE_ON)
            fprintf(stderr, "updating mean and variance estimates...\n") ;
          if (which == ON_STATS) {
            MRIaccumulateMaskedMeansAndVariances(mri_T1, mri_binary, mri_dof,
                                                 90, 100, mri_mean, mri_std) ;
            fprintf(stderr, "T1 = %d, binary = %d, mean = %2.1f\n",
                    (int)MRIgetVoxVal(mri_T1, 141,100,127,0),
                    MRIvox(mri_binary, 141,100,127),
                    MRIFvox(mri_mean, 141,100,127)) ;
          } else  /* computing means and vars for off */
            MRIaccumulateMaskedMeansAndVariances(mri_T1, mri_binary, mri_dof,
                                                 0, WM_MIN_VAL-1,
                                                 mri_mean, mri_std) ;
          MRIfree(&mri_T1) ;
        }
        MRIfree(&mri_binary) ;
      }

      if (which == BUILD_PRIORS) {
        mri = MRIcomputePriors(mri_priors, dof, NULL) ;
        MRIfree(&mri_priors) ;
        fprintf(stderr, "writing priors to %s...\n", out_fname) ;
      } else {
        MRIcomputeMaskedMeansAndStds(mri_mean, mri_std, mri_dof) ;
        mri_mean->dof = dof ;

        fprintf(stderr, "writing T1 means with %d dof to %s...\n", mri_mean->dof,
                out_fname) ;
        if (!which)
          MRIwrite(mri_mean, out_fname) ;
        else
          MRIappend(mri_mean, out_fname) ;
        MRIfree(&mri_mean) ;
        fprintf(stderr, "writing T1 variances to %s...\n", out_fname);
        if (dof <= 1)
          MRIreplaceValues(mri_std, mri_std, 0, 1) ;
        mri = mri_std ;
      }

      if (!which)
        MRIwrite(mri, out_fname) ;
      else
        MRIappend(mri, out_fname) ;
      MRIfree(&mri) ;
    }
  }
  else {
    /* for each subject specified on cmd line */

    if (xform_mean_fname) {
      m_xform_mean = MatrixAlloc(4,4,MATRIX_REAL) ;
      /* m_xform_covariance = MatrixAlloc(12,12,MATRIX_REAL) ;*/
    }

    dof = 0;
    for (i = 1 ; i < argc-1 ; i++) {

      if (*argv[i] == '-') {
        /* don't do transform for next subject */
        no_transform = 1 ;
        continue ;
      }
      dof++ ;

      subject_name = argv[i] ;
      sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, T1_name);
      fprintf(stderr, "%d of %d: reading %s...\n", i, argc-2, fname) ;
      mri_T1 = MRIread(fname) ;
      if (!mri_T1)
        ErrorExit(ERROR_NOFILE,"%s: could not open volume %s",Progname,fname);
      check_mri(mri_T1) ;

      if (binarize)
        MRIbinarize(mri_T1, mri_T1, binarize, 0, 1) ;
      if (erode) {
        int i ;
        printf("eroding input %d times\n", erode) ;
        for (i = 0 ; i < erode ; i++)
          MRIerode(mri_T1, mri_T1) ;
      }
      if (open) {
        int i ;
        printf("opening input %d times\n", open) ;
        for (i = 0 ; i < open ; i++)
          MRIerode(mri_T1, mri_T1) ;
        for (i = 0 ; i < open ; i++)
          MRIdilate(mri_T1, mri_T1) ;
      }

      check_mri(mri_T1) ;
      if (transform_fname) {

        sprintf(fname, "%s/%s/mri/transforms/%s",
                subjects_dir, subject_name, transform_fname) ;

        fprintf(stderr, "reading transform %s...\n", fname) ;
        ////////////////////////////////////////////////////////
#if 1
        {
          TRANSFORM *transform ;
          transform = TransformRead(fname) ;
          if (transform == NULL)
            ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",Progname, fname) ;
          mri_tmp = TransformApply(transform, mri_T1, NULL) ;
          if (DIAG_VERBOSE_ON)
            MRIwrite(mri_tmp, "t1.mgz") ;
          TransformFree(&transform) ;
        }
#else
        lta = LTAreadEx(fname);
        if (lta == NULL)
          ErrorExit(ERROR_NOFILE,
                    "%s: could not open transform file %s\n",
                    Progname, fname) ;
        printf("transform matrix -----------------------\n");
        MatrixPrint(stdout,lta->xforms[0].m_L);
        /* LTAtransform() runs either MRIapplyRASlinearTransform()
        for RAS2RAS or MRIlinearTransform() for Vox2Vox. */
        /* MRIlinearTransform() calls MRIlinearTransformInterp() */
        mri_tmp = LTAtransform(mri_T1, NULL, lta);
        printf("----- -----------------------\n");
        LTAfree(&lta);
#endif
        MRIfree(&mri_T1);
        mri_T1 = mri_tmp ; // reassign pointers
        if (DIAG_VERBOSE_ON)
          fprintf(stderr, "transform application complete.\n") ;
      }

      if (!mri_mean) {
        mri_mean =
          MRIalloc(mri_T1->width, mri_T1->height, mri_T1->depth, MRI_FLOAT) ;
        mri_std =
          MRIalloc(mri_T1->width, mri_T1->height, mri_T1->depth, MRI_FLOAT) ;
        if (!mri_mean || !mri_std)
          ErrorExit(ERROR_NOMEMORY, "%s: could not allocate templates.\n",
                    Progname) ;
        // if(transform_fname == NULL){
        if (DIAG_VERBOSE_ON)
          printf("Copying geometry\n");
        MRIcopyHeader(mri_T1,mri_mean);
        MRIcopyHeader(mri_T1,mri_std);
        // }
      }

      check_mri(mri_mean) ;
      if (!stats_only) {
        if (DIAG_VERBOSE_ON)
          fprintf(stderr, "updating mean and variance estimates...\n") ;
        MRIaccumulateMeansAndVariances(mri_T1, mri_mean, mri_std) ;
      }

      check_mri(mri_mean) ;
      if (DIAG_VERBOSE_ON)
        MRIwrite(mri_mean, "t2.mgz") ;
      MRIfree(&mri_T1) ;
      no_transform = 0;
    } /* end loop over subjects */

    if (xform_mean_fname) {
      FILE   *fp ;
      VECTOR *v = NULL, *vT = NULL ;
      MATRIX *m_vvT = NULL ;
      int    rows, cols ;

      nsubjects = sno ;

      fp = fopen(xform_covariance_fname, "w") ;
      if (!fp)
        ErrorExit(ERROR_NOFILE, "%s: could not open covariance file %s",
                  Progname, xform_covariance_fname) ;
      fprintf(fp, "nsubjects=%d\n", nsubjects) ;

      MatrixScalarMul(m_xform_mean, 1.0/(double)nsubjects, m_xform_mean) ;
      printf("means:\n") ;
      MatrixPrint(stdout, m_xform_mean) ;
      MatrixAsciiWrite(xform_mean_fname, m_xform_mean) ;

      /* subtract the mean from each transform */
      rows = m_xform_mean->rows ;
      cols = m_xform_mean->cols ;
      for (sno = 0 ; sno < nsubjects ; sno++) {
        MatrixSubtract(m_xforms[sno], m_xform_mean, m_xforms[sno]) ;
        v = MatrixReshape(m_xforms[sno], v, rows*cols, 1) ;
        vT = MatrixTranspose(v, vT) ;
        m_vvT = MatrixMultiply(v, vT, m_vvT) ;
        if (!m_xform_covariance)
          m_xform_covariance =
            MatrixAlloc(m_vvT->rows, m_vvT->cols,MATRIX_REAL) ;
        MatrixAdd(m_vvT, m_xform_covariance, m_xform_covariance) ;
        MatrixAsciiWriteInto(fp, m_xforms[sno]) ;
      }

      MatrixScalarMul(m_xform_covariance, 1.0/(double)nsubjects,
                      m_xform_covariance) ;
      printf("covariance:\n") ;
      MatrixPrint(stdout, m_xform_covariance) ;
      MatrixAsciiWriteInto(fp, m_xform_covariance) ;
      fclose(fp) ;
      if (stats_only)
        exit(0) ;
    }

    MRIcomputeMeansAndStds(mri_mean, mri_std, dof) ;
    check_mri(mri_mean) ;
    check_mri(mri_std) ;

    mri_mean->dof = dof ;

    if (smooth) {
      MRI *mri_kernel, *mri_smooth ;

      printf("applying smoothing kernel\n") ;
      mri_kernel = MRIgaussian1d(smooth, 100) ;
      mri_smooth = MRIconvolveGaussian(mri_mean, NULL, mri_kernel) ;
      MRIfree(&mri_kernel) ;
      MRIfree(&mri_mean) ;
      mri_mean = mri_smooth ;
    }
    fprintf(stderr, "\nwriting T1 means with %d dof to %s...\n", mri_mean->dof,
            out_fname) ;
    MRIwrite(mri_mean, out_fname) ;
    MRIfree(&mri_mean) ;
    if (dof <= 1) /* can't calculate variances - set them to reasonable val */
    {
      //               src      dst
      MRIreplaceValues(mri_std, mri_std, 0, 1) ;
    }
    if (!novar) {
      // mri_std contains the variance here  (does it?? I don't think so -- BRF)
      if (!var_fname) {
        fprintf(stderr, "\nwriting T1 standard deviations to %s...\n", out_fname);
        MRIappend(mri_std, out_fname) ;
      } else {
        fprintf(stderr, "\nwriting T1 standard deviations to %s...\n", var_fname);
        MRIwrite(mri_std, var_fname) ;
      }
    }
    MRIfree(&mri_std) ;
    if (mri)
      MRIfree(&mri);
  } /* end if binarize */
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
  if (!stricmp(option, "T1")) {
    T1_name = argv[2] ;
    fprintf(stderr,"reading T1 volume from directory '%s'\n",T1_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "statsonly")) {
    stats_only = 1 ;
  } else if (!stricmp(option, "smooth")) {
    smooth = atof(argv[2]) ;
    nargs = 1 ;
    printf("smoothing output with Gaussian sigma=%2.2f\n", smooth) ;
  } else if (!stricmp(option, "sdir")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "binarize")) {
    binarize = atoi(argv[2]) ;
    printf("binarizing input volume with thresh=%d...\n", binarize) ;
    nargs = 1 ;
  } else if (!stricmp(option, "erode")) {
    erode = atoi(argv[2]) ;
    printf("eroding input volume %d times...\n", erode) ;
    nargs = 1 ;
  } else if (!stricmp(option, "open")) {
    open = atoi(argv[2]) ;
    printf("opening input volume %d times...\n", open) ;
    nargs = 1 ;
  } else if (!stricmp(option, "novar")) {
    novar = 1 ;
    printf("disabling writing of variance estimates\n") ;
  } else switch (toupper(*option)) {
    case 'X':
      xform_mean_fname = argv[2] ;
      xform_covariance_fname = argv[3] ;
      printf("writing means (%s) and covariances (%s) of xforms\n",
             xform_mean_fname, xform_covariance_fname) ;
      nargs = 2 ;
      break ;
    case 'S':
    case 'V':
      var_fname = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, "writing variances to %s...\n", var_fname) ;
      break ;
    case 'B':
      binary_name = argv[2] ;
      fprintf(stderr, "generating binary template from %s volume\n",
              binary_name) ;
      nargs = 1 ;
      break ;
    case 'N':
      first_transform = 1 ;  /* don't use transform on first volume */
      break ;
    case 'T':
      transform_fname = argv[2] ;
      fprintf(stderr, "applying transformation %s to each volume\n",
              transform_fname) ;
      nargs = 1 ;
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
  printf("usage: %s <subject> <subject> ... <output volume>\n", Progname) ;
  printf("options are:\n") ;
  printf("-s <std>     - write out std deviations to file <fname>\n") ;
  printf("-t <xform>   - apply <xform> before computing stats\n") ;
  exit(code) ;
}

int
MRIaccumulateMeansAndVariances(MRI *mri, MRI *mri_mean, MRI *mri_std) {
  int    x, y, z, width, height, depth ;
  float  val, *pmean, *pstd ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pmean = &MRIFvox(mri_mean, 0, y, z) ;
      pstd = &MRIFvox(mri_std, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        val = MRIgetVoxVal(mri,x,y,z,0) ;
        if (x == DEBUG_X && y == DEBUG_Y && z == DEBUG_Z)
          DiagBreak() ;
#if 1
        *pmean++ += (float) val ;
        *pstd++ += ((float) val)*((float) val) ;
#else
        MRIFvox(mri_mean,x,y,z) += val ;
        MRIFvox(mri_std,x,y,z) += val*val ;
#endif
      }
    }
  }
  return(NO_ERROR) ;
}

int
MRIcomputeMeansAndStds(MRI *mri_mean, MRI *mri_std, int ndof) {
  int    x, y, z, width, height, depth ;
  float  sum, sum_sq, mean, var ;

  width = mri_std->width ;
  height = mri_std->height ;
  depth = mri_std->depth ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == DEBUG_X && y == DEBUG_Y && z == DEBUG_Z)
          DiagBreak() ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        sum = MRIFvox(mri_mean,x,y,z) ;
        sum_sq = MRIFvox(mri_std,x,y,z) / ndof ;
        mean = MRIFvox(mri_mean,x,y,z) = sum / ndof ;
        var = sum_sq - mean*mean;
        if (fabs(var) > 1e10 || fabs(mean)  > 1e10)
          DiagBreak() ;
        MRIFvox(mri_std,x,y,z) = sqrt(var) ;
      }
    }
  }
  return(NO_ERROR) ;
}
int
MRIcomputeMaskedMeansAndStds(MRI *mri_mean, MRI *mri_std, MRI *mri_dof) {
  int    x, y, z, width, height, depth, ndof ;
  float  sum, sum_sq, mean, var ;

  width = mri_std->width ;
  height = mri_std->height ;
  depth = mri_std->depth ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == DEBUG_X && y == DEBUG_Y && z == DEBUG_Z)
          DiagBreak() ;
        sum = MRIFvox(mri_mean,x,y,z) ;
        ndof = MRIvox(mri_dof, x, y, z) ;
        if (!ndof)   /* variance will be 0 in any case */
          ndof = 1 ;
        sum_sq = MRIFvox(mri_std,x,y,z) / ndof ;
        mean = MRIFvox(mri_mean,x,y,z) = sum / ndof ;
        var = sum_sq - mean*mean;
        MRIFvox(mri_std,x,y,z) = sqrt(var) ;
      }
    }
  }
  return(NO_ERROR) ;
}

MRI *
MRIfloatToChar(MRI *mri_src, MRI *mri_dst) {
  int   width, height, depth/*, x, y, z, out_val*/ ;
  /*  float fmax, fmin ;*/

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR) ;
#if 1
  MRIcopy(mri_src, mri_dst) ;
#else
  MRIvalRange(mri_src, &fmin, &fmax) ;
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {}
    }
  }
#endif
  return(mri_dst) ;
}

#if 0
MRI *
MRIbinarizeEditting(MRI *mri_src, MRI *mri_dst) {
  int     width, height, depth, x, y, z, val ;
  BUFTYPE *

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  mri_dst = MRIalloc(width, height, depth, MRI_UCHAR) ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRIvalRange(mri_src, &fmin, &fmax) ;
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {}
    }
  }

  return(mri_dst) ;
}
#endif
int
MRIaccumulateMaskedMeansAndVariances(MRI *mri, MRI *mri_mask, MRI *mri_dof,
                                     float low_val,
                                     float hi_val,MRI *mri_mean,MRI *mri_std) {
  int     x, y, z, width, height, depth ;
  float   val ;
  BUFTYPE *pmask, mask ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pmask = &MRIvox(mri_mask, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        if (x == DEBUG_X && y == DEBUG_Y && z == DEBUG_Z)
          DiagBreak() ;
        mask = *pmask++ ;
        if (mask >= low_val && mask <= hi_val) {
          val = MRIgetVoxVal(mri,x,y,z,0) ;
          MRIFvox(mri_mean,x,y,z) += val ;
          MRIFvox(mri_std,x,y,z) += val*val ;
          MRIvox(mri_dof,x,y,z)++ ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}
static MRI *
MRIupdatePriors(MRI *mri_binary, MRI *mri_priors) {
  int     width, height, depth, x, y, z ;
  BUFTYPE *pbin ;
  float   prob ;

  width = mri_binary->width ;
  height = mri_binary->height ;
  depth = mri_binary->depth ;
  if (!mri_priors) {
    mri_priors = MRIalloc(width, height, depth, MRI_FLOAT) ;
  }

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pbin = &MRIvox(mri_binary, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        if (x == DEBUG_X && y == DEBUG_Y && z == DEBUG_Z)
          DiagBreak() ;
        prob = *pbin++ / 100.0f ;
        MRIFvox(mri_priors, x, y, z) += prob ;
      }
    }
  }
  return(mri_priors) ;
}

static MRI *
MRIcomputePriors(MRI *mri_priors, int ndof, MRI *mri_char_priors) {
  int     width, height, depth, x, y, z ;
  BUFTYPE *pchar_prior, char_prior ;
  float   *pprior, prior ;

  width = mri_priors->width ;
  height = mri_priors->height ;
  depth = mri_priors->depth ;
  if (!mri_char_priors) {
    mri_char_priors = MRIalloc(width, height, depth, MRI_UCHAR) ;
  }

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pprior = &MRIFvox(mri_priors, 0, y, z) ;
      pchar_prior = &MRIvox(mri_char_priors, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        if (x == DEBUG_X && y == DEBUG_Y && z == DEBUG_Z)
          DiagBreak() ;
        prior = *pprior++ ;
        if (prior > 0)
          DiagBreak() ;
        if (prior > 10)
          DiagBreak() ;
        char_prior = (BUFTYPE)nint(100.0f*prior/(float)ndof) ;
        if (char_prior > 101)
          DiagBreak() ;
        *pchar_prior++ = char_prior ;
      }
    }
  }
  return(mri_char_priors) ;
}
