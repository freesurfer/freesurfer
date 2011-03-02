/**
 * @file  mri_aseg_edit_train
 * @brief use a previously trained a classifier to correct an aseg
 *
 * use SVMs to correct an aseg.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:13 $
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


#include "timer.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "const.h"
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"
#include "fio.h"
#include "mrishash.h"
#include "cvector.h"
#include "svm.h"
#include "version.h"
#include "voxlist.h"
#include "cma.h"
#include "class_array.h"

static char vcid[] = "$Id: mri_aseg_edit_reclassify.c,v 1.5 2011/03/02 00:04:13 nicks Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;


/*-------------------------------- DATA ----------------------------*/

char *Progname ;

static float sigmas[] = { 0, .5, 1.0, 2.0 } ;
#define NSCALES (sizeof(sigmas) / sizeof(sigmas[0]))

static int target_label = 17 ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  SVM          *svm ;
  MRI          *mri_aseg, *mri_norm, *mri_grad[NSCALES], *mri_kernel, *mri_smooth[NSCALES],
               *mri_laplacian[NSCALES], *mri_dtrans, *mri_dtrans_grad, *mri_2nd_deriv_s[NSCALES] ;
  char         **av, *norm_name, *input_aseg_name, *output_aseg_name, *svm_name ;
  int          ac, nargs, ninputs, i ;
  struct timeb start ;
  int          msec, minutes, seconds, x, y, z, nchanged ;
  VOXEL_LIST   *vl_border ;
  float        *svm_inputs = NULL, svm_out ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_aseg_edit_reclassify.c,v 1.5 2011/03/02 00:04:13 nicks Exp $", "$Name:  $");
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

  TimerStart(&start) ;

  // <aseg in> <norm> <svm> <aseg out>
  if (argc < 5)
    usage_exit() ;

  input_aseg_name = argv[1] ;
  norm_name  = argv[2] ;
  svm_name = argv[3] ;
  output_aseg_name = argv[4] ;
  
  printf("reading aseg: %s, norm: %s, SVM: %s and writing output to %s...\n", 
         input_aseg_name, norm_name, svm_name, output_aseg_name) ;

#define WSIZE 3

  // every scale at each window location for each of the dx,dy,dz grads and the original image
  ninputs = (WSIZE*WSIZE*WSIZE*NSCALES*(3+1)) ;
  svm = SVMread(svm_name) ;
  if (svm == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not load SVM from %s", Progname, svm_name) ;

  mri_aseg = MRIread(input_aseg_name) ;
  if (!mri_aseg)
    ErrorExit(ERROR_NOFILE, "%s: could not read aseg file %s",
              Progname, input_aseg_name) ;

  vl_border = VLSTcreate(mri_aseg, target_label, target_label, NULL, 0, 1) ;

  mri_norm = MRIread(norm_name) ;
  if (!mri_norm)
    ErrorExit(ERROR_NOFILE, "%s: could not read norm volume %s",
              Progname, norm_name) ;
  for (i = 0 ; i < NSCALES ; i++)
  {
    mri_kernel = MRIgaussian1d(sigmas[i], -1) ;
    mri_smooth[i] = MRIconvolveGaussian(mri_norm, NULL, mri_kernel) ;
    mri_grad[i] = MRIsobel(mri_smooth[i], NULL, NULL) ;
    mri_laplacian[i] = MRIlaplacian(mri_smooth[i], NULL) ;
    mri_2nd_deriv_s[i] = MRI2ndDirectionalDerivative(mri_smooth[i], NULL, 0, -1, 0) ;
    MRIfree(&mri_kernel) ;
  }
  mri_dtrans = MRIdistanceTransform(mri_aseg, NULL, target_label, 10, DTRANS_MODE_SIGNED, NULL) ;
  mri_dtrans_grad = MRIsobel(mri_dtrans, NULL, NULL) ;
  
  for (nchanged = i = 0 ; i < vl_border->nvox ; i++)
  {
    x = vl_border->xi[i] ;  y = vl_border->yi[i] ;  z = vl_border->zi[i] ; 
    if  (x == Gx && y == Gy && z == Gz)
      DiagBreak() ;
    svm_inputs = CAbuildInputsAtVoxel(vl_border, i, mri_smooth, mri_grad, 
                                      mri_laplacian, mri_dtrans, mri_dtrans_grad, mri_2nd_deriv_s,
                                      WSIZE, NSCALES, svm_inputs, 0) ;
    svm_out = SVMclassify(svm, svm_inputs) ;
    if  (x == Gx && y == Gy && z == Gz)
      printf("voxel(%d, %d, %d): svm_out = %2.3f\n", Gx, Gy, Gz, svm_out) ;
    if (svm_out < 0)
    {
      nchanged++ ;
      MRIsetVoxVal(mri_aseg, x, y, z, 0, Left_undetermined) ;
    }
  }
  
  printf("%d voxels changed (%2.1f%%)\n", nchanged, 100*nchanged/(float)vl_border->nvox) ;
  printf("writing edited aseg to %s...\n", output_aseg_name) ;
  MRIwrite(mri_aseg, output_aseg_name) ;
  
  MRIfree(&mri_aseg) ; MRIfree(&mri_norm) ;
  for (i = 0 ; i < NSCALES ; i++)
  {
    MRIfree(&mri_smooth[i]) ; MRIfree(&mri_grad[i]) ;
  }
  
  VLSTfree(&vl_border) ; 

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("svm classification took %d minutes and %d seconds.\n",
         minutes, seconds) ;
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
  else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else switch (toupper(*option)) {
    case 'L':
      target_label = atoi(argv[2]) ;
      fprintf(stderr, "resegmenting target label %s (%d)\n", cma_label_to_name(target_label), target_label) ;
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
          "usage: %s -o <output subject> [options] \n"
          "\t\n\t<subject1> <subject2>... <output file> \n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will compute train an SVM classifier to correct an aseg\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

