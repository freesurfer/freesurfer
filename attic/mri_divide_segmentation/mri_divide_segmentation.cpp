/**
 * @file  main_template.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2014/05/28 20:27:57 $
 *    $Revision: 1.3 $
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
#include "const.h"
#include "timer.h"
#include "version.h"
#include "cma.h"
#include "colortab.h"
#include "matrix.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;


#define MAX_DIVISIONS 1000

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs, segno, indices[MAX_DIVISIONS] ;
  int          msec, minutes, seconds, nparts, i, num, label, mx, my, mz, x, y, z ;
  Timer start ;
  MRI          *mri ;
  double       cx, cy, cz, min_dist, dist, dx, dy, dz ;
  float        evalues[3], zf, zf_low, zf_high, ez_x, ez_y, ez_z ;
//  double       e1x, e1y, e1z, e2x, e2y, e2z, e3z, e3y, e3z ;
  MATRIX       *m_obs, *m_obs_T, *m_cov, *m_eig ;

  setRandomSeed(-1L) ;
  nargs = handle_version_option (argc, argv, "$Id: mri_divide_segmentation.c,v 1.3 2014/05/28 20:27:57 fischl Exp $", "$Name:  $");
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

  if (argc < 5)
    usage_exit(1) ;

  mri = MRIread(argv[1]) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s\n", Progname, argv[1]) ;
  if (mri->ct == NULL)
  {
    char fname[STRLEN], *fs ;
    fs = getenv("FREESURFER_HOME") ;
    if (fs == NULL)
      ErrorExit(ERROR_NOFILE, "%s: volume %s does not have an embedded color table\n", Progname, argv[1]) ;
    sprintf(fname, "%s/FreeSurferColorLUT.txt", fs) ;
    mri->ct = CTABreadASCII(fname) ;
    if (mri->ct == NULL)
      ErrorExit(ERROR_NOFILE, "%s: volume %s does not have an embedded color table\n", Progname, argv[1]) ;
  }

  segno = atoi(argv[2]) ;
  nparts = atoi(argv[3]) ;
  printf("dividing segmentation %s (%d) into %d parts along its eigen-axis\n", cma_label_to_name(segno), segno, nparts) ;
  for (i = 0 ; i < nparts ; i++)
  {
    char name[STRLEN] ;
    sprintf(name, "%s.div%d", cma_label_to_name(segno), i) ;
    indices[i] = CTABaddUniqueEntry(mri->ct, name, 50) ;
    mri->ct->entries[indices[i]]->TissueType = mri->ct->entries[segno]->TissueType ;
    printf("%s: index %d\n", name, indices[i]) ;
  }


  // compute centroid of annotation
  num = cx = cy = cz = 0 ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth; z++)
      {
	label = MRIgetVoxVal(mri, x, y, z, 0) ;
	if (label != segno)
	  continue ;
	cx += x ; cy += y ; cz += z ;
	num++ ;
      }
  if (num == 0) // unused parcellation
    ErrorExit(ERROR_BADPARM, "%s: could not find any voxels with index %d", Progname, segno) ;

  cx /= num ; cy /= num ; cz /= num ;

  // find vertex in annotation closest to centroid
  min_dist = 100000 ;
  mx = my = mz = -1 ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth; z++)
      {
	label = MRIgetVoxVal(mri, x, y, z, 0) ;
	if (label != segno)
	  continue ;
	dist = sqrt(SQR(x-cx)+SQR(y-cy)+SQR(z-cz));
	if (dist < min_dist) {
	  mx = x ; my = y; mz = z ;
	}
      }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("using voxel (%d, %d, %d) as closest (%2.3f mm) from centroid (%2.1f, %2.1f, %2.1f)\n",
           mx, my, mz, min_dist, cx, cy, cz) ;

  // now compute eigensystem around this vertex
  m_obs = MatrixAlloc(3, num, MATRIX_REAL) ;
  for (num = x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth; z++)
      {
	label = MRIgetVoxVal(mri, x, y, z, 0) ;
	if (label != segno)
	  continue ;
	dx = x - cx ; dy = y - cy ;  dz = z - cz ;
	*MATRIX_RELT(m_obs, 1, num+1) = dx ;
	*MATRIX_RELT(m_obs, 2, num+1) = dy ;
	*MATRIX_RELT(m_obs, 3, num+1) = dz ;
	num++ ;
      }

  m_obs_T = MatrixTranspose(m_obs, NULL) ;
  m_cov = MatrixMultiply(m_obs,m_obs_T, NULL) ;
  m_eig = MatrixEigenSystem(m_cov, evalues, NULL) ;
  ez_x = *MATRIX_RELT(m_eig, 1, 1) ;
  ez_y = *MATRIX_RELT(m_eig, 2, 1) ;
  ez_z = *MATRIX_RELT(m_eig, 3, 1) ;
  for (i = 2 ; i <= 3 ; i++)  // find eigenvector that is closest to z axis
  {
    if (fabs(*MATRIX_RELT(m_eig, 3, i)) > fabs(ez_z))
    {
      ez_x = *MATRIX_RELT(m_eig, 1, i) ;
      ez_y = *MATRIX_RELT(m_eig, 2, i) ;
      ez_z = *MATRIX_RELT(m_eig, 3, i) ;
    }
  }
  if (ez_z  < 0) // orient it anterior/posterior
  {
    ez_x *= -1 ;
    ez_y *= -1 ;
    ez_z *= -1 ;
  }

  //find the bounding box
  zf_low = 10000 ; zf_high = -zf_low ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth; z++)
      {
	label = MRIgetVoxVal(mri, x, y, z, 0) ;
	if (label != segno)
	  continue ;
        zf = (x-mx)*ez_x + (y-my)*ez_y + (z-mz)*ez_z ;
        if (zf < zf_low)
          zf_low = zf ;
        if (zf > zf_high)
          zf_high = zf ;
      }

  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth; z++)
      {
	int subdivision ;

	label = MRIgetVoxVal(mri, x, y, z, 0) ;
	if (label != segno)
	  continue ;
	zf = (x-mx)*ez_x + (y-my)*ez_y + (z-mz)*ez_z ;
	subdivision = floor((zf-zf_low)/((zf_high-zf_low+1)/nparts));
	if (subdivision < 0)
	  subdivision = 0;
	if (subdivision >= nparts)
	  subdivision = nparts-1;
	MRIsetVoxVal(mri, x, y, z, 0, indices[subdivision]); 
      }

  MatrixFree(&m_obs) ;
  MatrixFree(&m_cov) ;
  MatrixFree(&m_obs_T) ;
  printf("writing output to %s\n", argv[4]) ;
  MRIwrite(mri, argv[4]) ;


  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "segmentation subdivisioning took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  switch (toupper(*option)) {
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
  printf("usage: %s [options] <input aseg> <segmentation index> <# of subdivisions> <output aseg>\n",
         Progname) ;
  exit(code) ;
}





