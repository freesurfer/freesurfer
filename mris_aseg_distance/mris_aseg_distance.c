/**
 * @file  mris_aseg_distance.c
 * @brief computes a map of distances between each point in the cortex and
 * a subcortical centroid
 *
 * predictability of subcortical structure locations from cortical
 * points.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:26 $
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
#include "timer.h"
#include "version.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;


char *Progname ;

MRI *MRIdivideAseg(MRI *mri_src, MRI *mri_dst, int label, int nunits);
static void usage_exit(int code) ;
static char sdir[STRLEN] = "" ;

static int notransform = 0 ;
static char surf_name[STRLEN] = "white" ;
static int dist_flag = 0 ;
static int divide = 1 ;
static int dot_flag = 0 ;
static int normalize = 0 ;

int
main(int argc, char *argv[]) {
  char        **av, *subject, *hemi, *out_fname, fname[STRLEN], *cp ;
  int         ac, nargs, label ;
  MRI         *mri_aseg, *mri_out ;
	LTA         *lta ;
	MRI_SURFACE *mris ;
	double      xc, yc, zc, xv, yv, zv ;
	VECTOR      *v1, *v2 ;
	int         vno, l ;
	VERTEX      *v ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_aseg_distance.c,v 1.4 2011/03/02 00:04:26 nicks Exp $", "$Name:  $");
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

  if (argc < 5)
    usage_exit(1) ;

  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
	subject = argv[1] ;
	hemi = argv[2] ;
	label = atoi(argv[3]) ;
	out_fname = argv[4] ;

	sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, surf_name) ;
	mris = MRISread(fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read surface %s\n",Progname,fname);

	sprintf(fname, "%s/%s/mri/aseg.mgz", sdir, subject) ;
  mri_aseg = MRIread(fname) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read aseg volume %s\n", 
							Progname,fname);

  if (divide > 1)
  {
    MRI *mri_tmp ;
    mri_tmp = MRIclone(mri_aseg, NULL) ;
    MRIcopyLabel(mri_aseg, mri_tmp, label) ;
    MRIfree(&mri_aseg) ; mri_aseg = mri_tmp ;
    MRIdivideAseg(mri_aseg, mri_aseg, label, divide) ;
  }
  if (notransform == 0)
  {
    sprintf(fname, "%s/%s/mri/transforms/talairach.lta", sdir, subject) ;
    lta = LTAread(fname) ;
    if (lta == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not read transform %s\n", 
                Progname,fname);
  }
  else
    lta = NULL ;

  if (dist_flag || dot_flag)
    mri_out = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 1) ;  
  else
    mri_out = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3) ;  
  
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
	VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
  for (l = 0 ; l < divide ; l++)
  {
    MRIcomputeLabelCentroid(mri_aseg, label+l, &xc, &yc, &zc) ;
    printf("centroid of label %s (%d+%d) = (%2.3f, %2.3f, %2.3f)\n", 
           cma_label_to_name(label), label, l, xc, yc, zc) ;
    
    if (notransform == 0)
    {
      V3_X(v1) = xc ;  V3_Y(v1) = yc ;  V3_Z(v1) = zc ; 
      MatrixMultiply(lta->xforms[0].m_L, v1, v2) ;
      xc = V3_X(v2) ;  yc = V3_Y(v2) ;  zc = V3_Z(v2) ; 
      printf("centroid in tal coords = (%2.0f, %2.0f, %2.0f)\n",
             xc, yc, zc) ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
      if (notransform == 0)
      {
        V3_X(v1) = xv ;  V3_Y(v1) = yv ;  V3_Z(v1) = zv ; 
        MatrixMultiply(lta->xforms[0].m_L, v1, v2) ;
        xv = V3_X(v2) ;  yv = V3_Y(v2) ;  zv = V3_Z(v2) ; 
      }
      if (dist_flag)
      {
        double dist, dx, dy, dz ;
        dx = xv-xc ; dy = yv-yc ; dz = zv-zc ; 
        dist = sqrt(dx*dx + dy*dy + dz*dz) ;
        MRIsetVoxVal(mri_out, vno, 0, 0, 0, dist) ;
      }
      else if (dot_flag)
      {
        double dot, dx, dy, dz, nx, ny, nz ;
        dx = xv-xc ; dy = yv-yc ; dz = zv-zc ; 
        MRISvertexNormalInVoxelCoords(mris,
                                      mri_aseg, vno, &nx, &ny, &nz) ;
        dot = dx*nx + dy*ny + dz*nz ;
        MRIsetVoxVal(mri_out, vno, 0, 0, 0, dot) ;
      }
      else
      {
        MRIsetVoxVal(mri_out, vno, 0, 0, 0, xv-xc) ;
        MRIsetVoxVal(mri_out, vno, 0, 0, 1, yv-yc) ;
        MRIsetVoxVal(mri_out, vno, 0, 0, 2, zv-zc) ;
      }
    }

    sprintf(fname, out_fname, l) ;
    printf("writing output to %s\n", fname) ;
    MRIwrite(mri_out, fname) ;
  }
	MatrixFree(&v1) ; MatrixFree(&v2) ;
  if (Gdiag_fp)
    fclose(Gdiag_fp) ;
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
  if (!stricmp(option, "log")) {
    Gdiag_fp = fopen(argv[2], "w") ;
    if (Gdiag_fp == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not open log file %s", argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "dot")) {
    dot_flag = 1 ;
    printf("using dot product instead of distances...\n") ;
  } else if (!stricmp(option, "normalize")) {
    normalize = 1 ;
    printf("normalizing distances...\n") ;
  } else if (!stricmp(option, "divide")) {
    divide = atoi(argv[2]) ;
    printf("dividing aseg into %d units\n", divide) ;
    nargs = 1 ;
  }
  else switch (toupper(*option)) {
  case 'D':
    dist_flag = 1 ;
    printf("storing distances instead of full vector\n") ;
    break ;
  case 'N':
    notransform = 1 ;
    printf("not using Talairach transform\n") ;
    break ;
  case 'L':
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
  printf("usage: %s [options] <subject> <hemi> <label> <output map>\n",
				 Progname);
  printf("\tvalid options are:\n") ;
  exit(code) ;
}


#include "matrix.h"

