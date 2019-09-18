/**
 * @file  mri_train_autoencoder.c
 * @brief main program for creating and training a stacked autoencoder for feature extraction.
 *
H.-C. Shin, M. R. Orton, D. J. Collins, S. J. Doran, and M. O. Leach,
"Stacked Autoencoders for
Unsupervised Feature Learning and Multiple Organ Detectionin a Pilot Study
Using 4D Patient Data,"
IEEE Transaction on Pattern Analysis and Machine Intelligence, 2012.

 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2013/11/14 16:17:41 $
 *    $Revision: 1.2 $
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
#include "error.h"
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
#include "autoencoder.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static int encoder_type = NORMAL_AUTOENCODER ;
static int synthesize = 1 ;
static int whalf = 2 ;
static double tol = 1e-4;
static double momentum = .5 ;
static double dt = 1e-2;
static double scale = .5 ;
static char *read_fname = NULL ;

static int x0 = -1 ;
static int x1 ;
// for some reason y0 and y1 are defined in /usr/include/bits/mathcalls.h so have to include annoying _
static int y0_ ;
static int y1_ ;

static int z0 ;
static int z1 ;

SAE_INTEGRATION_PARMS parms ;

#define MAX_PYR_LEVELS 100 
static int nlevels = 4 ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  char   *in_fname, *out_fname ;
  int          msec, minutes, seconds, n ;
  Timer start ;
  MRI          *mri, *mri_orig, *mri_scaled, *mri_pyramid[MAX_PYR_LEVELS] ;
  SAE          *sae ;
  double       mean ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.integration_type  = INTEGRATE_GRADIENT_DESCENT ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_train_autoencoder.c,v 1.2 2013/11/14 16:17:41 fischl Exp $", "$Name:  $");
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

  parms.momentum  = momentum ;
  parms.orig_dt = parms.dt  = dt ;
  parms.tol  = tol ;
  in_fname = argv[1] ;
  parms.out_fname = out_fname = argv[2] ;
  if (argc < 3)
    usage_exit(1) ;
//  setRandomSeed(-1L) ;

  mri = MRIread(in_fname) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume from %s", Progname, in_fname) ;

  mri_orig = MRIcopy(mri, NULL) ;
  if (mri->type == MRI_UCHAR)
  {
    MRI *mri_tmp = MRIcloneDifferentType(mri, MRI_FLOAT) ;
    MRIscalarMul(mri, mri_tmp, 1.0/255.0) ;
//    MRIwrite(mri_tmp, "scaled.mgz") ;
    MRIfree(&mri) ; mri = mri_tmp ;
  }
  mri_scaled = MRIcopy(mri, NULL) ;
  mean = MRImeanFrame(mri_scaled, 0) ;
  MRIaddScalar(mri, mri, -mean) ;

  if (x0 >= 0)
  {
    MRI *mri_tmp ;
    mri_tmp = MRIextract(mri, NULL, x0, y0_, z0, x1-x0+1, y1_-y0_+1, z1-z0+1) ;
//    MRIwrite(mri_tmp, "ex.mgz") ;
    MRIfree(&mri) ;
    mri = mri_tmp ;
  }

  if (read_fname)
  {
    sae = SAEread(read_fname) ;
    if (sae == NULL)
      ErrorExit(Gerror, "") ;
    SAEaddLayer(sae, scale) ;
    nlevels = sae->nlevels ;
    whalf = sae->whalf ;
  }
  else
    sae = SAEalloc(whalf, nlevels, encoder_type, scale) ;
  if (sae == NULL)
    ErrorExit(Gerror, "") ;

  mri_pyramid[0] = mri ;
//  sprintf(fname, "pyr%d.mgz", n=0) ;
//  MRIwrite(mri_pyramid[n], fname) ;
  for (n = 1 ; n < nlevels ; n++)
  {
    mri_pyramid[n] = MRIreduce(mri_pyramid[n-1], NULL) ;
//    sprintf(fname, "pyr%d.mgz", n) ;
//    MRIwrite(mri_pyramid[n], fname) ;
  }

  SAEtrainFromMRI(sae, mri_pyramid, &parms) ;
  printf("writing autoencoder to %s\n", out_fname) ;
  SAEwrite(sae, out_fname) ;

  if (synthesize)
  {
    int x, y, z, ind  ;
    VECTOR *v ;
    char   path[STRLEN], fname[STRLEN] ;
    MRI    *mri_synth ;

    printf("synthesizing output volume using auto-encoder...\n") ;

    ind = (sae->first->v_output->rows+1)/2 ;
    FileNameRemoveExtension(out_fname, path) ;
    mri_synth = MRIclone(mri_scaled, NULL) ;
    for (x = sae->whalf ; x < mri_orig->width-sae->whalf; x++)
      for (y = sae->whalf ; y < mri_orig->height-sae->whalf ; y++)
	for (z = sae->whalf ; z < mri_orig->depth-sae->whalf ; z++)
	{
	  if (x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;
	  if (FZERO(MRIgetVoxVal(mri_orig, x, y, z, 0)))
	    continue ;
	  SAEfillInputVector(mri_pyramid, sae->nlevels, x, y, z, sae->whalf, sae->first->v_input) ;
	  v = SAEactivateNetwork(sae) ;
	  MRIsetVoxVal(mri_synth, x, y, z, 0, sae->first->v_output->rptr[ind][1]) ;
	} 
    MRIaddScalar(mri_synth, mri_synth, mean) ;
    if (mri_orig->type == MRI_UCHAR)
    {
      MRIscalarMul(mri_synth, mri_synth, 255.0f) ;
      mri_synth = MRIchangeType(mri_synth, MRI_UCHAR, 0, 255, 1) ;
    }
    sprintf(fname, "%s.out.mgz", path) ;
    printf("writing synthesized output to %s\n", fname) ;
    MRIwrite(mri_synth, fname) ;
    MRIfree(&mri_synth) ; 
  }

  SAEfree(&sae) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf(" Training took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gvx, Gvy, Gvz) ;
  }
  else switch (toupper(*option)) {
    case 'N':
      nlevels = atof(argv[2]) ;
      printf("using %d levels in Gaussian pyramid\n", nlevels) ;
      nargs = 1 ;
      break ;
    case 'B':
      parms.acceptance_sigma = atof(argv[2]) ;
      parms.proposal_sigma = atof(argv[3]) ;
      parms.integration_type = INTEGRATE_BOLTZMANN_MACHINE ;
      printf("minimizing using Boltzmann machine with acceptance/proposal = %2.2f / %2.2f\n",
	     parms.acceptance_sigma, parms.proposal_sigma) ;
      nargs = 2 ;
      break ;
    case 'C':
      parms.integration_type = INTEGRATE_CONJUGATE_GRADIENT ;
      printf("using Polak-Ribiere conjugate gradient minimization\n") ;
      break ;
    case '1':
      encoder_type = FOCUSED_AUTOENCODER ;
      printf("training focused auto-encoder\n") ;
      break ;
    case 'W':
      whalf = atoi(argv[2]) ;
      fprintf(stderr, "using half window size = %d\n", whalf) ;
      nargs = 1 ;
      break ;
    case 'X':
      x0 = atoi(argv[2]) ;
      x1 = atoi(argv[3]) ;
      y0_ = atoi(argv[4]) ;
      y1_ = atoi(argv[5]) ;
      z0 = atoi(argv[6]) ;
      z1 = atoi(argv[7]) ;
      nargs = 6 ;
      printf("extracting subregion X %d -->%d, y %d --> %d, z %d --> %d\n", x0, x1, y0_, y1_, z0, z1) ;
      break ;
    case 'R':
      read_fname = argv[2] ;
      nargs = 1 ;
      printf("reading previously trained auto-encoded from %s\n", read_fname) ;
      break ;
    case 'D':
      dt = atof(argv[2]) ;
      fprintf(stderr, "using dt = %e\n", dt) ;
      nargs = 1 ;
      break ;
    case 'S':
      scale = atof(argv[2]) ;
      fprintf(stderr, "scaling hidden layer to be %2.2f as large as input layer\n", scale) ;
      nargs = 1 ;
      break ;
    case 'M':
      momentum = atof(argv[2]) ;
      fprintf(stderr, "using momentum = %e\n", momentum) ;
      nargs = 1 ;
      break ;
    case 'T':
      tol = atof(argv[2]) ;
      fprintf(stderr, "using tol = %e\n", tol) ;
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
  printf("usage: %s [options] <inverse operator> <EEG/MEG data file>",
         Progname) ;
  printf(
    "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
  );
  printf("\tn - noise-sensitivity normalize inverse (default=1)") ;
  exit(code) ;
}





