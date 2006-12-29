/**
 * @file  xvmri.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:19 $
 *    $Revision: 1.16 $
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


/*----------------------------------------------------------------------
           File Name:  xvmri.h

           Author:

           Description:

           Conventions:

----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           INCLUDE FILES
----------------------------------------------------------------------*/

#ifndef XVMRI_H
#define XVMRI_H

#include "xvutil.h"
#include "mri.h"
#include "region.h"
#include "image.h"
#include "mrisurf.h"

int  mri_event_handler(XV_FRAME *xvf, Event *event,
                       DIMAGE *dimage, int *px, int *py, int *pz) ;
void XVMRIdrawPoint(XV_FRAME *xvf, int which, int view, int depth, MRI *mri,
                    int x, int y, int z, int color) ;
int  XVMRIdraw3DPoint(XV_FRAME *xvf, int which,int x, int y,int z,int color);

void XVMRIdrawRegion(XV_FRAME *xvf, int which, int view, int depth, MRI *mri,
                     MRI_REGION *reg, int color) ;
IMAGE *XVMRIshow(XV_FRAME *xvf, MRI *mri, int which, int slice) ;
IMAGE *XVMRIshowRange(XV_FRAME *xvf, MRI *mri, int which, int slice,
                      float fmin, float fmax) ;
IMAGE *XVMRIshowFrame(XV_FRAME *xvf, MRI *mri,int which,int slice,int frame);
IMAGE *XVMRIredisplayFrame(XV_FRAME *xvf, MRI *mri, int which, int slice,
                           int frame);
int  XVMRIinit(XV_FRAME *xvf, int view_row, int view_col) ;
void XVMRIshowAll(XV_FRAME *xvf) ;
int  XVMRIfree(MRI **pmri, int which) ;
int  XVMRIsetView(XV_FRAME *xvf, int which, int view) ;
int  XVMRIsetPoint(XV_FRAME *xvf, int which, int x, int y, int z) ;
int  XVMRIsetImageName(XV_FRAME *xvf, char *image_name) ;
int  XVMRIfixDepth(XV_FRAME *xvf, int which, int fix) ;
int  XVMRIspecifySurface(XV_FRAME *xvf, MRI_SURFACE *mrisurf) ;
int  XVMRIshowThreeViews(XV_FRAME *xvf) ;



extern int   mri_slices[] ;
extern IMAGE *Idisplay[] ;
extern int   mri_views[] ;
extern int   mri_depths[] ;
extern MRI   *mris[] ;
extern int   mri_frames[] ;
extern int   which_click ;
extern int   x_click ;
extern int   y_click ;
extern int   z_click ;
extern int   brush ;

#endif
