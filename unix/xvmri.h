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

void mri_event_handler(XV_FRAME *xvf,int depth,MRI *mri,Event *event,
                       DIMAGE *dimage, int view) ;

#endif
