/*  $Id: fgutil.h,v 1.1 1997/03/11 18:43:07 fischl Exp $   */
#ifndef _FGUTIL_H
#define _FGUTIL_H

#include "h_logz.h"

/* if GRAB_LOGMAP is set, then x is # of rings, and y is # of spokes */
typedef struct {
  int          sizex, sizey;
  LOGMAP_INFO  *lmi ;
  int          GrabOpts ;     /* bit field of GRAB_ values */
  IMAGE        *Icartesian ;  /* only used if GrabOpts & GRAB_LOGMAP */
  IMAGE        *Ilog ;        /* only used if GrabOpts & GRAB_LOGMAP */
  IMAGE        *Idsp ;        /* image out of DSP */
} FGINFO;

#define GRAB_DEFAULT 0x00
#define GRAB_FIELD   0x01 /* 640x240 Default is 640x480 */
#define GRAB_LOGMAP  0x02 /* Default is to return real image */
#define GRAB_DBUFF   0x04 /* Grab using double buffer - Default is single */
#define GRAB_EONLY   0x08 /* Grab even field - Default is odd */
#define GRAB_SUBS    0x10 /* Subsample to true field width */

extern int FGInitialize(int GrabOpts, FGINFO *ptr);
extern void FGClose(void);
extern int FGGrab(FGINFO *fgi, char *buffer);

#endif /* _FGUTIL_H */
