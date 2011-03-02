/**
 * @file  map.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
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


/*
  @(#)map.h 1.2
  8/10/95
*/
/*------------------------------------------------------------------------
      File Name:  logmap.h

         Author:  Bruce Fischl

        Created:  Jan. 1993

    Description:

------------------------------------------------------------------------*/


/*----------------------------------------------------------------------
             File Name:

                Author:

           Description:

----------------------------------------------------------------------*/


#ifndef MAP_H
#define MAP_H


#include "hips.h"
#include "macros.h"

typedef struct
{
  int   rows ;
  int   cols ;
  UCHAR *image ;
}
CMAP ;

typedef struct
{
  int  rows ;
  int  cols ;
  int  *image ;
}
IMAP ;

/*
  this macro will alow access to the map cell at row,col.  It
  is only used for 2-D maps.  1-D maps just index directly.
*/
#define MAPcell(m, x, y)  (*((m)->image + ((y) * m->cols) + x))

void   MapIHipsWrite(IMAP *map, char *fname) ;
void   MapCHipsWrite(CMAP *map, char *fname) ;
CMAP   *MapClone(CMAP *m) ;
void   MapCopy(CMAP *mSrc, CMAP *mDst) ;
void   MapAdd(CMAP *m1, CMAP *m2, CMAP *mSum) ;
IMAP   *MapIAlloc(int cols, int rows) ;
CMAP   *MapCAlloc(int cols, int rows) ;
void   MapIFree(IMAP *map) ;
void   MapCFree(CMAP *map) ;
void   MapIHipsWrite(IMAP *map, char *fname) ;
void   MapCHipsWrite(CMAP *map, char *fname) ;
void   MapCShowImage(CMAP *map, char *name) ;
void   MapIShowImage(IMAP *map, char *name) ;
void   MapSoapBubbleRelax(IMAP *mVal, CMAP *mVisited) ;
void   MapIPrint(IMAP *map, FILE *fp) ;
void   MapCPrint(CMAP *map, FILE *fp) ;



#endif
