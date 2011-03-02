/**
 * @file  mri_conform.c
 * @brief resample volume to isotropic 1mm^3
 *
 */
/*
 * Original Author: Christian Haselgrove
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.35 $
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
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include "mri.h"
#include "error.h"
#include "histo.h"
#include "mri_conform.h"

extern int errno;

MATRIX *MRIgetConformMatrix(MRI *mri)
{

  MRI *templ;
  MATRIX *m_resample ;

  if (mri->ras_good_flag == 0)
  {
    setDirectionCosine(mri, MRI_CORONAL);
  }

  templ = MRIallocHeader(256, 256, 256, MRI_UCHAR);

  templ->imnr0 = 1;
  templ->imnr1 = 256;
  templ->thick = 1.0;
  templ->ps = 1.0;
  templ->xsize = templ->ysize = templ->zsize = 1.0;
  templ->xstart = templ->ystart = templ->zstart = -128.0;
  templ->xend = templ->yend = templ->zend = 128.0;
  setDirectionCosine(templ, MRI_CORONAL); // sets c_(r,a,s) = 0
  // retain the src c_(r,a,s)
  templ->c_r =  mri->c_r;
  templ->c_a =  mri->c_a;
  templ->c_s =  mri->c_s;
  templ->ras_good_flag = 1; // use c_(r,a,s)
  templ->tr = mri->tr ;
  templ->te = mri->te ;
  templ->flip_angle = mri->flip_angle ;
  templ->ti = mri->ti ;

  m_resample = MRIgetResampleMatrix(mri, templ);

  MRIfree(&templ) ;

  return(m_resample);
}

MRI *MRIconform(MRI *mri)
{

  MRI *templ, *mri2, *res;

  res = MRIcopy(mri, NULL); /* don't mess with the input */

  if (res->ras_good_flag == 0)
  {
    setDirectionCosine(res, MRI_CORONAL);
  }

  templ = MRIallocHeader(256, 256, 256, MRI_UCHAR);

  templ->imnr0 = 1;
  templ->imnr1 = 256;
  templ->thick = 1.0;
  templ->ps = 1.0;
  templ->xsize = templ->ysize = templ->zsize = 1.0;
  templ->xstart = templ->ystart = templ->zstart = -128.0;
  templ->xend = templ->yend = templ->zend = 128.0;
  setDirectionCosine(templ, MRI_CORONAL); // sets c_(r,a,s) = 0
  // retain the c_(r,a,s)
  templ->c_r = res->c_r;
  templ->c_a = res->c_a;
  templ->c_s = res->c_s;
  templ->ras_good_flag = 1; // use c_(r,a,s)
  templ->tr = mri->tr ;
  templ->te = mri->te ;
  templ->flip_angle = mri->flip_angle ;
  templ->ti = mri->ti ;

  /* ----- change type if necessary ----- */
  if (res->type != templ->type)
  {
    mri2 = MRIchangeType(res, templ->type, 0.0, 0.999, FALSE);
    MRIfree(&res);
    if (mri2 == NULL)
      return(NULL);
    res = mri2;
  }

  /* ----- reslice if necessary ----- */
  if (res->xsize != templ->xsize || 
      res->ysize != templ->ysize || 
      res->zsize != templ->zsize ||
      res->width != templ->width || 
      res->height != templ->height || 
      res->depth != templ->depth ||
      res->x_r != templ->x_r || 
      res->x_a != templ->x_a || 
      res->x_s != templ->x_s ||
      res->y_r != templ->y_r || 
      res->y_a != templ->y_a || 
      res->y_s != templ->y_s ||
      res->z_r != templ->z_r || 
      res->z_a != templ->z_a || 
      res->z_s != templ->z_s)
  {
    mri2 = MRIresample(res, templ, SAMPLE_TRILINEAR);
    MRIfree(&res);
    if (mri2 == NULL)
      return(NULL);
    res = mri2;
  }

  return(res);

}  /*  end MRIconform()  */

