/**
 * @file  mris_ms_surface_CNR.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:33 $
 *    $Revision: 1.4 $
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


//
// mris_ms_surface_CNR.c
// compute CNR along a surface (white) from multiple input volumes
// original author: Xiao Han
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2011/03/02 00:04:33 $
// Revision       : $Revision: 1.4 $
//
////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mri.h"
#include "mrinorm.h"
#include "mri_conform.h"
#include "mri_identify.h"
#include "mrisurf.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "matrix.h"
#include "transform.h"
#include "version.h"
#include "label.h"


#define MAX_IMAGES 200

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

static int use_thickness = 0;
static int conform = 0 ;

char *Progname ;

static char sdir[80] = ""; //SUBJECTS_DIR
static char *sname = NULL;
static char *hemi = NULL;

static char *out_fname = NULL; /* Output CNR filename, curv or paint */

static char *thickness_fname = NULL; /* filename for surface thickness */

char *trgtypestring = "paint";
int trgtype = MRI_VOLUME_TYPE_UNKNOWN;

static int nSmoothSteps = 60;

int debugflag = 0;
int debugvtx = 0;

static void usage_exit(int code) ;

int
main(int argc, char *argv[]) {
  char   **av, *in_fname = NULL, *cp;
  char filename[80];
  int    ac, nargs;
  MRIS    *mris;
  MRI *mri_flash[MAX_IMAGES];
  MRI *mri_gm_profile;
  MRI *mri_wm_profile;
  MRI *mri_gm_profile_orig;
  MRI *mri_wm_profile_orig;
  MRI *mri_gm_cov_components;
  MRI *mri_wm_cov_components;
  int num_of_cov_components;
  MRI *mri_cnr;
  MRI *mri_tmp ;
  MRI *mri_weight[MAX_IMAGES]; //this is the relative weighting for each component
  int    msec, minutes, seconds;
  struct timeb start ;
  int nvolumes, nvolumes_total;
  VERTEX *vertex;
  int index, i, j, vno;
  double cx, cy, cz;
  Real  vx, vy, vz;
  Real value_in, value_out;
  double cnr;
  MATRIX *SW1, *SW2, *SW, *InvSW;
  float *weight, weight_norm;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_ms_surface_CNR.c,v 1.4 2011/03/02 00:04:33 nicks Exp $", "$Name:  $");
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

  if (argc < 2)
    usage_exit(1) ;

  printf("command line parsing finished\n");

  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  if (out_fname == NULL) {
    ErrorExit(ERROR_BADPARM,
              "%s: please specify output CNR file name.\n", Progname) ;
  }
  if (sname == NULL && hemi == NULL) {
    ErrorExit(ERROR_BADPARM,
              "%s: please specify subject and hemisphere names.\n", Progname) ;
  }


  //////////////////////////////////////////////////////////////////////////////////
  /*** Read in the input multi-echo volumes ***/
  nvolumes = 0 ;
  for (i = 1 ; i < argc; i++) {
    in_fname = argv[i] ;
    printf("reading %s...\n", in_fname) ;

    mri_flash[nvolumes] = MRIread(in_fname) ;
    if (mri_flash[nvolumes] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read volume %s",
                Progname, in_fname) ;
    /* conform will convert all data to UCHAR, which will reduce data resolution*/
    printf("%s read in. \n", in_fname) ;
    if (conform) {

      printf("embedding and interpolating volume\n") ;
      mri_tmp = MRIconform(mri_flash[nvolumes]) ;
      MRIfree(&mri_flash[nvolumes]);
      mri_flash[nvolumes] = mri_tmp ;
    }

    /* Change all volumes to float type for convenience */
    if (mri_flash[nvolumes]->type != MRI_FLOAT) {
      printf("Volume %d type is %d\n", nvolumes+1, mri_flash[nvolumes]->type);
      printf("Change data to float type \n");
      mri_tmp = MRIchangeType(mri_flash[nvolumes], MRI_FLOAT, 0, 1.0, 1);
      MRIfree(&mri_flash[nvolumes]);
      mri_flash[nvolumes] = mri_tmp; //swap
    }

    nvolumes++ ;
  }

  printf("All data read in\n");

  ///////////////////////////////////////////////////////////////////////////
  nvolumes_total = nvolumes ;   /* all volumes read in */

  for (i = 0 ; i < nvolumes ; i++) {
    for (j = i+1 ; j < nvolumes ; j++) {
      if ((mri_flash[i]->width != mri_flash[j]->width) ||
          (mri_flash[i]->height != mri_flash[j]->height) ||
          (mri_flash[i]->depth != mri_flash[j]->depth))
        ErrorExit(ERROR_BADPARM, "%s:\nvolumes %d (type %d) and %d (type %d) don't match (%d x %d x %d) vs (%d x %d x %d)\n",
                  Progname, i, mri_flash[i]->type, j, mri_flash[j]->type, mri_flash[i]->width,
                  mri_flash[i]->height, mri_flash[i]->depth,
                  mri_flash[j]->width, mri_flash[j]->height, mri_flash[j]->depth) ;
    }
  }

  /*** Read in the white surface ***/
  sprintf(in_fname, "%s/%s/surf/%s.white",sdir,sname,hemi);
  printf("reading white surface from %s...\n", in_fname) ;
  mris = MRISread(in_fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read white surface %s",
              Progname, in_fname) ;

  printf("white surface file read in.\n");

  if (use_thickness) {
    printf("reading thickness\n");

    if (thickness_fname == NULL) {
      if (MRISreadCurvatureFile(mris, "thickness") != 0) {
        ErrorExit(ERROR_NOFILE, "%s:could not read thickness", Progname);
      }
    } else {
      if (MRISreadCurvatureFile(mris, thickness_fname) != 0) {
        ErrorExit(ERROR_NOFILE, "%s:could not read thickness from file %s", Progname, thickness_fname);
      }
    }
  }

  //sample volume data to surface
  //assume GM is sampled at 0.5*thickness; WM is sampled at 1mm inwards

  mri_gm_profile = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, nvolumes_total) ;
  mri_wm_profile = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, nvolumes_total) ;

  MRIScomputeNormals(mris);

  printf("sample GM and WM values from input volumes ... \n");

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];

    if (!(vno % 25000))
      fprintf(stdout, "%d of %d vertices processed\n", vno, mris->nvertices) ;

    if (debugflag && debugvtx == vno)
      printf("(nx, ny, nz) = (%g, %g, %g) \n", vertex->nx, vertex->ny, vertex->nz);

    /* Take the voxel location at half-thickness as GM */
    if (use_thickness) {
      cx = vertex->x + 0.5*vertex->curv*vertex->nx;
      cy = vertex->y + 0.5*vertex->curv*vertex->ny;
      cz = vertex->z + 0.5*vertex->curv*vertex->nz;
    } else {
      cx = vertex->x + vertex->nx;
      cy = vertex->y + vertex->ny;
      cz = vertex->z + vertex->nz;
    }

    MRIsurfaceRASToVoxel(mri_flash[0], cx, cy, cz, &vx, &vy, &vz);

    if (debugflag && debugvtx == vno)
      printf("(vx, vy, vz) = (%g, %g, %g) \n", vx, vy, vz);

    for (i = 0 ; i < nvolumes ; i++) {
      MRIsampleVolumeType(mri_flash[i], vx,  vy, vz, &value_out, SAMPLE_TRILINEAR) ;
      MRIFseq_vox(mri_gm_profile, vno, 0, 0, i) = value_out ;
    }


    //wm
    cx = vertex->x - vertex->nx;
    cy = vertex->y - vertex->ny;
    cz = vertex->z - vertex->nz;
    MRIsurfaceRASToVoxel(mri_flash[0], cx, cy, cz, &vx, &vy, &vz);

    for (i = 0 ; i < nvolumes ; i++) {
      MRIsampleVolumeType(mri_flash[i], vx,  vy, vz, &value_in, SAMPLE_TRILINEAR) ;
      MRIFseq_vox(mri_wm_profile, vno, 0, 0, i) = value_in ;
    }

  }

  printf("Free original volumes ...\n");
  for (i=0; i < nvolumes_total; i++) {
    MRIfree(&mri_flash[i]);
  }

  //backup a copy since one copy is needed to compute means; seems unnecessary! since CNR will use mean and covariance, not original values
  mri_gm_profile_orig = MRIcopy( mri_gm_profile, NULL);
  mri_wm_profile_orig = MRIcopy( mri_wm_profile, NULL);

  //now compute covariance matrix at every surface location
  //need cross-correlation for every pair of channels
  num_of_cov_components = nvolumes_total*(nvolumes_total+1)/2;

  mri_gm_cov_components = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, num_of_cov_components) ;
  mri_wm_cov_components = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, num_of_cov_components) ;

  printf("compute cross-correlation ...\n");
  //compute the vertex-wise cross-correlation
  for (vno = 0; vno < mris->nvertices; vno++) {
    index = 0;
    for (i=0; i < nvolumes_total ;i++)
      for (j=i; j < nvolumes_total; j++) {
        MRIFseq_vox(mri_gm_cov_components, vno,0,0,index) =  MRIFseq_vox(mri_gm_profile, vno, 0, 0, i) * MRIFseq_vox(mri_gm_profile, vno, 0, 0, j);
        MRIFseq_vox(mri_wm_cov_components, vno,0,0,index) =  MRIFseq_vox(mri_wm_profile, vno, 0, 0, i) * MRIFseq_vox(mri_wm_profile, vno, 0, 0, j);
        index++;
      }
  }


  if (nSmoothSteps > 0) { //this is a weighted mean, not true mean; should be equivalent! Well, weighted one may be better
    printf("compute averages (weighted one) ... \n");
    MRISsmoothFrames(mris,mri_gm_cov_components,nSmoothSteps);
    MRISsmoothFrames(mris,mri_wm_cov_components,nSmoothSteps);
    MRISsmoothFrames(mris,mri_wm_profile,nSmoothSteps);
    MRISsmoothFrames(mris,mri_gm_profile,nSmoothSteps);
  }

  printf("Now compute vertex-wise CNR on the surface ...\n");
  /* Just intend to initialize the array */
  mri_cnr =  MRIcopyMRIS(NULL, mris, 0, "curv");
  for (i = 0 ; i < nvolumes_total; i++) {
    mri_weight[i] = MRIcopy(mri_cnr, NULL);
  }

  SW1 = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
  SW2 = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
  SW = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
  weight = (float *)malloc(nvolumes_total*sizeof(float));
  for (vno = 0; vno < mris->nvertices; vno++) {
    //construct the vertex-wise covariance matrix for WM and GM
    index = 0;
    for (i=0; i < nvolumes_total ;i++)
      for (j=i; j < nvolumes_total; j++) {
        SW1->rptr[i+1][j+1] = MRIFseq_vox(mri_gm_cov_components, vno,0,0,index)  - MRIFseq_vox(mri_gm_profile, vno, 0, 0, i) * MRIFseq_vox(mri_gm_profile, vno, 0, 0, j) ; //sigma^2 = mean(x^2) - mean(x)*mean(x)
        SW2->rptr[i+1][j+1] = MRIFseq_vox(mri_wm_cov_components, vno,0,0,index)  - MRIFseq_vox(mri_wm_profile, vno, 0, 0, i) * MRIFseq_vox(mri_wm_profile, vno, 0, 0, j) ; //sigma^2 = mean(x^2) - mean(x)*mean(x)
        index++;
      }

    //take the average of GM and WM cov and the final
    for (i=1; i <= nvolumes_total ;i++) {
      for (j=i; j <= nvolumes_total; j++) {
        SW->rptr[i][j] = 0.5*(SW1->rptr[i][j] + SW2->rptr[i][j]);
        SW->rptr[j][i] = SW->rptr[i][j]; //symmetric
      }

      SW->rptr[i][i] += 1e-30; //regularize a little bit
    }

    InvSW = MatrixInverse(SW, NULL);

    /* weight would be the optimal weighting, but just intermediate value for CNR computation too */
    for (i = 0; i < nvolumes_total ; i++) {
      weight[i]= 0.0;
      for (j= 0; j < nvolumes_total; j++) {
        weight[i] += InvSW->rptr[i+1][j+1] *(MRIFseq_vox(mri_wm_profile, vno, 0, 0, j) - MRIFseq_vox(mri_gm_profile, vno, 0, 0, j));
      }
    }

    /* compute CNR */
    cnr = 0 ;
    weight_norm = 0;
    for (i=0; i < nvolumes_total; i++) {
      cnr += weight[i]*(MRIFseq_vox(mri_wm_profile, vno, 0, 0, i) - MRIFseq_vox(mri_gm_profile, vno, 0, 0, i));
      weight_norm += weight[i]*weight[i];
    }
    if (cnr < 0) cnr = 0;
    else cnr = sqrt(cnr);
    MRIsetVoxVal(mri_cnr,vno, 0, 0, 0, cnr);

    weight_norm = sqrt(weight_norm);
    for (i=0; i < nvolumes_total; i++) {
      if (weight[i] < 0) weight[i] = -weight[i];
      weight[i] = weight[i]/(weight_norm + 1e-30);
      MRIsetVoxVal(mri_weight[i],vno, 0, 0, 0, weight[i]);
    }
  }


  printf("Output CNR and weighting files ...\n");
  MRIScopyMRI(mris, mri_cnr, 0, "curv");

  sprintf(filename,"%s_cnr.w", out_fname);
  if (!strcmp(trgtypestring,"paint") || !strcmp(trgtypestring,"w")) {
    MRISwriteCurvatureToWFile(mris,filename);
  } else if (!strcmp(trgtypestring,"curv")) {
    MRISwriteCurvature(mris,filename);
  } else
    fprintf(stderr, "ERROR unknown output file format.\n");

  for (i=0; i < nvolumes_total; i++) {
    sprintf(filename,"%s_weight_%d.w", out_fname,i);
    MRIScopyMRI(mris, mri_weight[i], 0, "curv");
    MRISwriteCurvatureToWFile(mris,filename);
  }

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("CNR computation took %d minutes and %d seconds.\n", minutes, seconds) ;

  for (i = 0 ; i < nvolumes_total; i++) {
    MRIfree(&mri_weight[i]);
  }
  MatrixFree(&SW1);
  MatrixFree(&SW2);
  MatrixFree(&SW);
  MatrixFree(&InvSW);
  MRISfree(&mris);
  MRIfree(&mri_wm_profile);
  MRIfree(&mri_wm_profile_orig);
  MRIfree(&mri_gm_profile);
  MRIfree(&mri_gm_profile_orig);
  MRIfree(&mri_gm_cov_components);
  MRIfree(&mri_wm_cov_components);
  MRIfree(&mri_cnr);

  exit(0);
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

  if (!stricmp(option, "out") ||
      !stricmp(option, "out_file") ||
      !stricmp(option, "out_name") ||
      !stricmp(option, "cnr")) {
    out_fname = argv[2];
    printf("Output CNR map to file %s\n", out_fname);
    nargs = 1 ;
  } else if (!stricmp(option, "sdir")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "sname")) {
    sname = argv[2];
    printf("using %s as subject name \n", sname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "hemi")) {
    hemi = argv[2];
    printf("hemisphere = %s \n", hemi) ;
    nargs = 1 ;
  } else if (!stricmp(option, "debug")) {
    debugflag = 1;
    debugvtx = atoi(argv[2]);
    nargs = 1;
  } else if (!stricmp(option, "use_thickness")) {
    use_thickness = 1;
    printf("Use thickness to guide the choice of GM point \n");
  } else if (!stricmp(option, "thickness") ||
             !stricmp(option, "thickness_file") ||
             !stricmp(option, "thickness_fname")) {
    thickness_fname = argv[2];
    printf("Use file %s as thickness map \n", thickness_fname);
    nargs = 1 ;
  } else if (!stricmp(option, "trg_type")) {
    trgtypestring = argv[2];
    trgtype = string_to_type(trgtypestring);
    nargs = 1 ;
  } else if (!stricmp(option, "nsmooth")) {
    nSmoothSteps = atoi(argv[2]);
    nargs = 1;
    printf("Perform %d steps of smoothing of output CNR map\n", nSmoothSteps);
  } else if (!stricmp(option, "debug")) {
    debugflag = 1;
    debugvtx = atoi(argv[2]);
    nargs = 1;
  } else {
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    usage_exit(0) ;
    exit(1) ;
  }



  return(nargs) ;
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s vol1 vol2 ... \n", Progname) ;
  printf("\t This program computes the image CNR sampled along the white surface \n");
  printf("Options includes:\n");
  printf("\t -debug %%d to set the surface vertex to debug \n");
  printf("\t -cnr %%s to set the output CNR map  filename (prefix) \n");
  printf("\t -sdir %%s to set the SUBJECTS_DIR\n");
  printf("\t -sname %%s to set the subject name \n");
  printf("\t -hemi %%s to choose the hemisphere \n");
  printf("\t -thickness %%s to set the filename for thicknessmap\n");
  printf("\t -nsmooth %%d number of smoothing steps\n");
  printf("\t -trg_type  %%s output format (default = paint) \n");
  printf("\t -use_thickness  use thickness to guide the pick of GM sample (by default) \n");

  exit(code) ;
}

