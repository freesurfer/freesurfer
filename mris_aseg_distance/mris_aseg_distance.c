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
 *    $Author: fischl $
 *    $Date: 2007/05/04 17:32:34 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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

static void usage_exit(int code) ;
static char sdir[STRLEN] = "" ;

static char surf_name[STRLEN] = "white" ;

int MRIcomputeLabelCentroid(MRI *mri_aseg, int label, 
														double *pxc, double *pyc, double *pzc) ;
int
main(int argc, char *argv[]) {
  char        **av, *subject, *hemi, *out_fname, fname[STRLEN], *cp ;
  int         ac, nargs, label ;
  MRI         *mri_aseg, *mri_out ;
	LTA         *lta ;
	MRI_SURFACE *mris ;
	double      xc, yc, zc, xv, yv, zv ;
	VECTOR      *v1, *v2 ;
	int         vno ;
	VERTEX      *v ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_aseg_distance.c,v 1.1 2007/05/04 17:32:34 fischl Exp $", "$Name:  $");
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

	sprintf(fname, "%s/%s/mri/transforms/talairach.lta", sdir, subject) ;
  lta = LTAread(fname) ;
  if (lta == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read transform %s\n", 
							Progname,fname);

	MRIcomputeLabelCentroid(mri_aseg, label, &xc, &yc, &zc) ;
	printf("centroid of label %s (%d) = (%2.3f, %2.3f, %2.3f)\n", 
				 cma_label_to_name(label), label, xc, yc, zc) ;

	mri_out = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3) ;  

	v1 = VectorAlloc(4, MATRIX_REAL) ;
	v2 = VectorAlloc(4, MATRIX_REAL) ;
	VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
	V3_X(v1) = xc ;  V3_Y(v1) = yc ;  V3_Z(v1) = zc ; 
	MatrixMultiply(lta->xforms[0].m_L, v1, v2) ;
	xc = V3_X(v2) ;  yc = V3_Y(v2) ;  zc = V3_Z(v2) ; 
	for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
		V3_X(v1) = xv ;  V3_Y(v1) = yv ;  V3_Z(v1) = zv ; 
		MatrixMultiply(lta->xforms[0].m_L, v1, v2) ;
		xv = V3_X(v2) ;  yv = V3_Y(v2) ;  zv = V3_Z(v2) ; 
		MRIsetVoxVal(mri_out, vno, 0, 0, 0, xv-xc) ;
		MRIsetVoxVal(mri_out, vno, 0, 0, 1, yv-yc) ;
		MRIsetVoxVal(mri_out, vno, 0, 0, 2, zv-zc) ;
	}

	printf("writing output to %s\n", out_fname) ;
	MRIwrite(mri_out, out_fname) ;
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
  }
  else switch (toupper(*option)) {
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

int
MRIcomputeLabelCentroid(MRI *mri_aseg, int label, 
												double *pxc, double *pyc, double *pzc)
{
	int    l, x, y, z, num ;
	double xc, yc, zc ;

	for (xc = yc = zc = 0.0, num = x = 0 ; x < mri_aseg->width ; x++)
	{
		for (y = 0 ; y < mri_aseg->height ; y++)
		{
			for (z = 0 ; z < mri_aseg->depth ; z++)
			{
				l = nint(MRIgetVoxVal(mri_aseg, x, y, z, 0)) ;
				if (l != label)
					continue ;
				num++ ;
				xc += x; yc += y; zc += z;
			}
		}
	}

	if (num>0)
	{
		xc /= num;  yc /= num;  zc /= num; 
	}
	*pxc = xc ; *pyc = yc ; *pzc = zc ;

	return(NO_ERROR) ;
}

