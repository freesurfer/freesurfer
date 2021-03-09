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
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "gca.h"
#include "transform.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
#if 0
static int test(MRI *mri1, MRI *mri2, MRI *mri3, MATRIX *m_vol1_to_vol2_ras) ;
#endif

const char *Progname ;
static void usage_exit(int code) ;

static char *histo_parms = NULL ;
static int write_flag = 0 ;
static char *log_fname = NULL ;
static const char *parc_dir = "parc" ;
static const char *T1_name = "flash/T1.mgh" ;
static const char *PD_name = "flash/PD.mgh" ;

static char *xform_name = NULL ;

static char subjects_dir[STRLEN] ;

int
main(int argc, char *argv[]) {
  TRANSFORM    *transform = NULL ;
  char         **av, fname[STRLEN], *gca_fname, *subject_name, *cp ;
  int          ac, nargs, i, n ;
  int          msec, minutes, seconds, nsubjects ;
  Timer start ;
  GCA          *gca ;
  MRI          *mri_parc, *mri_T1, *mri_PD ;
  FILE         *fp ;

  nargs = handleVersionOption(argc, argv, "mri_ca_tissue_parms");
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

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
    if (argc < 3)
      usage_exit(1) ;
  }


  gca_fname = argv[1] ;
  nsubjects = argc-2 ;
  printf("computing average tissue parameters on %d subject\n",
         nsubjects) ;

  n = 0 ;

  printf("reading GCA from %s...\n", gca_fname) ;
  gca = GCAread(gca_fname) ;

  for (i = 0 ; i < nsubjects ; i++) {
    subject_name = argv[i+2] ;
    printf("processing subject %s, %d of %d...\n", subject_name,i+1,
           nsubjects);
    sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, parc_dir) ;
    if (DIAG_VERBOSE_ON)
      printf("reading parcellation from %s...\n", fname) ;
    mri_parc = MRIread(fname) ;
    if (!mri_parc)
      ErrorExit(ERROR_NOFILE, "%s: could not read parcellation file %s",
                Progname, fname) ;

    sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, T1_name) ;
    if (DIAG_VERBOSE_ON)
      printf("reading co-registered T1 from %s...\n", fname) ;
    mri_T1 = MRIread(fname) ;
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE, "%s: could not read T1 data from file %s",
                Progname, fname) ;

    sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, PD_name) ;
    if (DIAG_VERBOSE_ON)
      printf("reading co-registered T1 from %s...\n", fname) ;
    mri_PD = MRIread(fname) ;
    if (!mri_PD)
      ErrorExit(ERROR_NOFILE, "%s: could not read PD data from file %s",
                Progname, fname) ;


    if (xform_name) {
      /*      VECTOR *v_tmp, *v_tmp2 ;*/

      sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, xform_name) ;
      printf("reading xform from %s...\n", fname) ;
      transform = TransformRead(fname) ;
      if (!transform)
        ErrorExit(ERROR_NOFILE, "%s: could not read xform from %s",
                  Progname, fname) ;
#if 0
      v_tmp = VectorAlloc(4,MATRIX_REAL) ;
      *MATRIX_RELT(v_tmp,4,1)=1.0 ;
      v_tmp2 = MatrixMultiply(lta->xforms[0].m_L, v_tmp, NULL) ;
      printf("RAS (0,0,0) -->\n") ;
      MatrixPrint(stdout, v_tmp2) ;
#endif

      if (transform->type == LINEAR_RAS_TO_RAS) {
        MATRIX *m_L ;
        m_L = ((LTA *)transform->xform)->xforms[0].m_L ;
        MRIrasXformToVoxelXform(mri_parc, mri_T1, m_L,m_L) ;
      }
#if 0
      v_tmp2 = MatrixMultiply(lta->xforms[0].m_L, v_tmp, v_tmp2) ;
      printf("voxel (0,0,0) -->\n") ;
      MatrixPrint(stdout, v_tmp2) ;
      VectorFree(&v_tmp) ;
      VectorFree(&v_tmp2) ;
      test(mri_parc, mri_T1, mri_PD, lta->xforms[0].m_L) ;
#endif
    }
    if (histo_parms)
      GCAhistogramTissueStatistics(gca,mri_T1,mri_PD,mri_parc,transform,histo_parms);
#if 0
    else
      GCAaccumulateTissueStatistics(gca, mri_T1, mri_PD, mri_parc, transform) ;
#endif

    MRIfree(&mri_parc) ;
    MRIfree(&mri_T1) ;
    MRIfree(&mri_PD) ;
  }
  GCAnormalizeTissueStatistics(gca) ;

  if (log_fname) {
    printf("writing tissue parameters to %s\n", log_fname) ;
    fp = fopen(log_fname, "w") ;
    for (n = 1 ; n < MAX_GCA_LABELS ; n++) {
      GCA_TISSUE_PARMS *gca_tp ;

      gca_tp = &gca->tissue_parms[n] ;
      if (gca_tp->total_training <= 0)
        continue ;
      fprintf(fp, "%d  %f  %f\n", n, gca_tp->T1_mean, gca_tp->PD_mean) ;
    }
    fclose(fp) ;
  }

  if (write_flag)
    GCAwrite(gca, gca_fname) ;
  GCAfree(&gca) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("tissue parameter statistic calculation took %d minutes"
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
  if (!stricmp(option, "T1")) {
    T1_name = argv[2] ;
    nargs = 1 ;
    printf("reading T1 data from subject's mri/%s directory\n",
           T1_name) ;
  } else if (!stricmp(option, "PD")) {
    PD_name = argv[2] ;
    nargs = 1 ;
    printf("reading PD map from subject's mri/%s directory\n",
           PD_name) ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as SUBJECTS_DIR\n", subjects_dir) ;
  } else if (!stricmp(option, "xform")) {
    xform_name = argv[2] ;
    nargs = 1 ;
    printf("reading and applying xform from %s for each subject\n",
           xform_name) ;
  } else if (!stricmp(option, "PARC_DIR")) {
    parc_dir = argv[2] ;
    nargs = 1 ;
    printf("reading parcellation from subject's mri/%s directory\n",
           parc_dir) ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  } else switch (toupper(*option)) {
    case 'H':
      histo_parms = argv[2] ;
      nargs = 1 ;
      printf("writing T1/PD histograms to %s...\n", histo_parms) ;
      break ;
      break ;
    case 'L':
      log_fname = argv[2] ;
      nargs = 1 ;
      printf("writing T1/PD class info to %s...\n", log_fname) ;
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
  printf("usage: %s [options] <subject 1> <subject 2> ... <output file>\n",
         Progname) ;
  printf(
    "\t-spacing  - spacing of classifiers in canonical space\n");
  printf("\t-gradient - use intensity gradient as input to classifier.\n") ;
  exit(code) ;
}
#if 0
static int
test(MRI *mri1, MRI *mri2, MRI *mri3, MATRIX *m_vol1_to_vol2_ras) {
  VECTOR *v_test, *v_vox ;
  float  x_ras1, y_ras1, z_ras1, x_ras2, y_ras2, z_ras2, x_vox1, y_vox1,
  z_vox1, x_vox2, y_vox2, z_vox2 ;
  MATRIX  *m_vol2_vox2ras, *m_vol2_ras2vox, *m_vol1_ras2vox, *m_vol1_vox2ras,
  *m_vol3_ras2vox, *m_vol3_vox2ras ;
  int     val ;


  v_test = VectorAlloc(4, MATRIX_REAL) ;
  m_vol1_vox2ras = MRIgetVoxelToRasXform(mri1) ;
  m_vol2_vox2ras = MRIgetVoxelToRasXform(mri2) ;
  m_vol1_ras2vox = MRIgetRasToVoxelXform(mri1) ;
  m_vol2_ras2vox = MRIgetRasToVoxelXform(mri2) ;
  m_vol3_vox2ras = MRIgetVoxelToRasXform(mri3) ;
  m_vol3_ras2vox = MRIgetRasToVoxelXform(mri3) ;

  x_ras1 = 126.50 ;
  y_ras1 = -125.500 ;
  z_ras1 = 127.50 ;
  V3_X(v_test) = x_ras1 ;
  V3_Y(v_test) = y_ras1 ;
  V3_Z(v_test) = z_ras1 ;
  *MATRIX_RELT(v_test, 4, 1) = 1.0 ;
  v_vox = MatrixMultiply(m_vol1_ras2vox, v_test, NULL) ;
  x_vox1 = V3_X(v_vox) ;
  y_vox1 = V3_Y(v_vox) ;
  z_vox1 = V3_Z(v_vox) ;
  val = MRISvox(mri1, nint(x_vox1), nint(y_vox1), nint(z_vox1)) ;
  printf("VOL1: ras (%1.1f, %1.1f, %1.1f) --> VOX (%1.1f, %1.1f, %1.1f) = %d\n",
         x_ras1, y_ras1, z_ras1, x_vox1, y_vox1, z_vox1, val) ;


  x_ras2 = 76.5421 ;
  y_ras2 = 138.5352 ;
  z_ras2 = 96.0910 ;
  V3_X(v_test) = x_ras2 ;
  V3_Y(v_test) = y_ras2 ;
  V3_Z(v_test) = z_ras2 ;
  *MATRIX_RELT(v_test, 4, 1) = 1.0 ;
  v_vox = MatrixMultiply(m_vol2_ras2vox, v_test, NULL) ;
  x_vox2 = V3_X(v_vox) ;
  y_vox2 = V3_Y(v_vox) ;
  z_vox2 = V3_Z(v_vox) ;
  val = MRISvox(mri2, nint(x_vox2), nint(y_vox2), nint(z_vox2)) ;
  printf("VOL2: ras (%2.1f, %2.1f, %2.1f) --> VOX (%2.1f, %2.1f, %2.1f) = %d\n",
         x_ras2, y_ras2, z_ras2, x_vox2, y_vox2, z_vox2, val) ;



  MatrixFree(&v_test) ;
  return(NO_ERROR) ;
}

#endif
