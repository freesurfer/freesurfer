/**
 * @file  fgutil.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
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


/*  $Id: fgutil.h,v 1.4 2011/03/02 00:04:09 nicks Exp $   */
#ifndef _FGUTIL_H
#define _FGUTIL_H

#include "h_logz.h"

/* if GRAB_LOGMAP is set, then x is # of rings, and y is # of spokes */
typedef struct
{
  int          sizex, sizey;
  LOGMAP_INFO  *lmi ;
  int          GrabOpts ;     /* bit field of GRAB_ values */
  IMAGE        *Icartesian ;  /* only used if GrabOpts & GRAB_LOGMAP */
  IMAGE        *Ilog ;        /* only used if GrabOpts & GRAB_LOGMAP */
  IMAGE        *Idsp ;        /* image out of DSP */
}
FGINFO;

#define GRAB_DEFAULT 0x00
#define GRAB_FIELD   0x01 /* 640x240 Default is 640x480 */
#define GRAB_LOGMAP  0x02 /* Default is to return real image */
#define GRAB_DBUFF   0x04 /* Grab using double buffer - Default is single */
#define GRAB_EONLY   0x08 /* Grab even field - Default is odd */
#define GRAB_SUBS    0x10 /* Subsample to true field width */

extern int FGInitialize(int GrabOpts, FGINFO *ptr);
extern void FGClose(void);
extern int FGGrab(FGINFO *fgi, char *buffer);
extern int FGsetLogmapType(FGINFO *fgi, int logmap_type) ;

#endif /* _FGUTIL_H */
