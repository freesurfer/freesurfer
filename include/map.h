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
} CMAP ;

typedef struct
{
  int  rows ;
  int  cols ;
  int  *image ;
} IMAP ;

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
