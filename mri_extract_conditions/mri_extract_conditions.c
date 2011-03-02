/**
 * @file  mri_extract_conditions.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:15 $
 *    $Revision: 1.5 $
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
#include "matrix.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "version.h"

static char vcid[] = "$Id: mri_extract_conditions.c,v 1.5 2011/03/02 00:04:15 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int detrend = 0 ;
char *Progname ;
static int *conditions = NULL ;
static float *timepoints = NULL ;

MRI *MRIextractConditions(MRI *mri_in, int *conditions, MRI *mri_out) ;
MRI *MRIdetrendVolume(MRI *mri_in, int *conditions, MRI *mri_out) ;

int
main(int argc, char *argv[]) {
  char        **av, *in_vol, *out_vol, *paradigm_fname, line[STRLEN] ;
  int         ac, nargs, t ;
  MRI         *mri_in, *mri_out ;
  FILE        *fp ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_extract_conditions.c,v 1.5 2011/03/02 00:04:15 nicks Exp $", "$Name:  $");
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

  if (argc < 4)
    usage_exit() ;

  in_vol = argv[1] ;
  paradigm_fname = argv[2] ;
  out_vol = argv[3] ;

  printf("reading volume from %s...\n", in_vol) ;
  mri_in = MRIread(in_vol) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname,
              in_vol) ;

  printf("reading paradigm file %s with %d timepoints\n", paradigm_fname, mri_in->nframes) ;
  conditions = (int *)calloc(mri_in->nframes, sizeof(int)) ;
  timepoints = (float *)calloc(mri_in->nframes, sizeof(float)) ;
  fp = fopen(paradigm_fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOMEMORY, "%s: could not read paradigm file %s", Progname,paradigm_fname) ;

  for (t = 0 ; t < mri_in->nframes ; t++) {
    if (fgetl(line, STRLEN-1, fp) == NULL)
      ErrorExit(ERROR_BADFILE, "%s: could only read %d points out of %s",
                Progname, t, paradigm_fname) ;
    if (sscanf(line, "%f %d", &timepoints[t], &conditions[t]) != 2)
      ErrorExit(ERROR_BADFILE, "%s: line %d malformed (%s)",
                Progname, t, line) ;
  }

  fclose(fp) ;

  if (detrend)
    MRIdetrendVolume(mri_in, conditions, mri_in) ;

  mri_out = MRIextractConditions(mri_in, conditions, NULL) ;
  printf("writing output to %s.\n", out_vol) ;
  MRIwrite(mri_out, out_vol) ;

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
    case 'D':
      detrend = 1 ;
      printf("detrending...\n") ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
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
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  printf(
    "usage: %s [options] <input volume> <paradigm file> <output file>\n",
    Progname) ;
  printf("where valid options are\n\t-detrend  - detrend input\n") ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program extract conditions from a run and output them into a single file\n") ;
  fprintf(stderr, "-imageoffset <image offset> - set offset to use\n") ;
  fprintf(stderr, "-detrend                    - detrend each voxel\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

MRI *
MRIdetrendVolume(MRI *mri_in, int *conditions, MRI *mri_out) {
  int  x, y, z, t, ntime, nnull, i ;
  MATRIX  *m_X, *m_Xt = NULL, *m_XtX = NULL, *m_XtXinv = NULL, *m = NULL ;
  VECTOR  *v_out, *v_p ;
  float   val, mean, a, b ;

  if (!mri_out)
    mri_out = MRIcopy(mri_in, NULL) ;

  ntime = mri_in->nframes ;
  for (nnull = t = 0 ; t < ntime ; t++)
    if (conditions[t] == 0)
      nnull++ ;

  m_X = MatrixAlloc(nnull, 2, MATRIX_REAL) ;
  v_out = VectorAlloc(nnull, MATRIX_REAL) ;
  v_p = VectorAlloc(2, MATRIX_REAL) ;

  for (x = 0 ; x < mri_out->width ; x++) {
    for (y = 0 ; y < mri_out->height ; y++) {
      for (z = 0 ; z < mri_out->depth ; z++) {
        if (x == 32 && y == 33 && z == 9)
          DiagBreak() ;
        for (mean = 0.0, i = t = 0 ; t < ntime ; t++) {
          if (conditions[t] > 0)   /* detrend based on null condition only */
            continue ;
          switch (mri_in->type) {
          default:
            ErrorExit(ERROR_UNSUPPORTED, "MRIdetrendVolume: input type %d unsupported",
                      mri_in->type) ;
          case MRI_FLOAT:
            val = MRIFseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_SHORT:
            val = MRISseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_INT:
            val = MRIIseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_UCHAR:
            val = MRIseq_vox(mri_in, x, y, z, t) ;
            break ;
          }
          *MATRIX_RELT(v_out, i+1, 1) = val ;
          *MATRIX_RELT(m_X, i+1, 1) = t ;
          *MATRIX_RELT(m_X, i+1, 2) = 1 ;
          mean += val ;
          i++ ;
        }
        mean /= (float)i ;
        m_Xt = MatrixTranspose(m_X, m_Xt) ;
        m_XtX = MatrixMultiply(m_Xt, m_X, m_XtX) ;
        m_XtXinv = MatrixInverse(m_XtX, m_XtXinv) ;
        if (m_XtXinv == NULL)
          ErrorExit(ERROR_BADPARM, "non-invertible matrix at %d,%d,%d",
                    x, y, z) ;
        m = MatrixMultiply(m_XtXinv, m_Xt, m) ;
        v_p = MatrixMultiply(m, v_out, v_p) ;
        a = VECTOR_ELT(v_p, 1) ;
        b = VECTOR_ELT(v_p, 2) ;

        /* detrend by subtracting (a*t-b - mean) from each time point */
        for (t = 0 ; t < ntime ; t++) {
          switch (mri_in->type) {
          default:
            ErrorExit(ERROR_UNSUPPORTED, "MRIdetrendVolume: input type %d unsupported",
                      mri_in->type) ;
          case MRI_FLOAT:
            val = MRIFseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_SHORT:
            val = MRISseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_INT:
            val = MRIIseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_UCHAR:
            val = MRIseq_vox(mri_in, x, y, z, t) ;
            break ;
          }
          val -= a*t+b-mean ;
          switch (mri_in->type) {
          default:
            ErrorExit(ERROR_UNSUPPORTED, "MRIdetrendVolume: input type %d unsupported",
                      mri_in->type) ;
          case MRI_FLOAT:
            MRIFseq_vox(mri_in, x, y, z, t) = val ;
            break ;
          case MRI_SHORT:
            MRISseq_vox(mri_in, x, y, z, t) = (short)nint(val) ;
            break ;
          case MRI_INT:
            MRIIseq_vox(mri_in, x, y, z, t) = nint(val) ;
            break ;
          case MRI_UCHAR:
            MRIseq_vox(mri_in, x, y, z, t) = (BUFTYPE)nint(val) ;
            break ;
          }
        }
      }
    }
  }

  MatrixFree(&m_X) ;
  MatrixFree(&m_Xt) ;
  MatrixFree(&m_XtXinv) ;
  MatrixFree(&m) ;
  MatrixFree(&m_XtX) ;
  VectorFree(&v_out) ;
  VectorFree(&v_p) ;
  return(mri_out) ;
}

#define MAX_CONDITIONS 1000

MRI *
MRIextractConditions(MRI *mri_in, int *conditions, MRI *mri_out) {
  int   x, y, z, t, ncond, ntime, width, height, depth, counts[MAX_CONDITIONS], c ;
  float means[MAX_CONDITIONS], val, vars[MAX_CONDITIONS] ;

  ntime = mri_in->nframes ;
  for (ncond = t = 0 ; t < ntime ; t++)
    if (conditions[t] > ncond)
      ncond = conditions[t] ;
  ncond++ ;

  printf("extracting %d conditions + null\n", ncond-1) ;
  width = mri_in->width ;
  height = mri_in->height ;
  depth = mri_in->depth ;
  mri_out = MRIallocSequence(width, height, depth, MRI_FLOAT, ncond*2) ;
  for (x = 0 ; x < mri_out->width ; x++) {
    for (y = 0 ; y < mri_out->height ; y++) {
      for (z = 0 ; z < mri_out->depth ; z++) {
        if (x == 32 && y == 33 && z == 9)
          DiagBreak() ;
        memset(means, 0, ncond*sizeof(float)) ;
        memset(vars, 0, ncond*sizeof(float)) ;
        memset(counts, 0, ncond*sizeof(int)) ;
        for (t = 0 ; t < ntime ; t++) {
          switch (mri_in->type) {
          default:
            ErrorExit(ERROR_UNSUPPORTED, "MRIextractConditions: input type %d unsupported",
                      mri_in->type) ;
          case MRI_FLOAT:
            val = MRIFseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_SHORT:
            val = MRISseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_INT:
            val = MRIIseq_vox(mri_in, x, y, z, t) ;
            break ;
          case MRI_UCHAR:
            val = MRIseq_vox(mri_in, x, y, z, t) ;
            break ;
          }
          means[conditions[t]] += val ;
          vars[conditions[t]] += val*val ;
          counts[conditions[t]]++ ;
        }
        for (c = 0 ; c < ncond ; c++) {
          if (counts[c] <= 0)
            continue ;
          means[c] /= (float)counts[c] ;
          vars[c] = vars[c] / (float)counts[c] - means[c]*means[c] ;
          MRIFseq_vox(mri_out, x, y, z, c*2) = means[c] ;
          MRIFseq_vox(mri_out, x, y, z, c*2+1) = sqrt(vars[c]) ;
        }
      }
    }
  }

  return(mri_out) ;
}

