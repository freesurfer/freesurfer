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

void mri_event_handler(XV_FRAME *xvf,int depth,MRI *mri,Event *event,
                       DIMAGE *dimage, int view, int *px, int *py, int *pz) ;
void XVMRIdrawPoint(XV_FRAME *xvf, int which, int view, int depth, MRI *mri,
                    int x, int y, int z, int color) ;

void XVMRIdrawRegion(XV_FRAME *xvf, int which, int view, int depth, MRI *mri,
                MRI_REGION *reg, int color) ;
#endif
