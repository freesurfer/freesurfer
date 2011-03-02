/**
 * @file  map.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:12 $
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
  @(#)map.c 1.1
  4/4/94
*/
/*----------------------------------------------------------------------
           File Name:
               map.c

           Author:
             Bruce Fischl with algorithms stolen from Rich Wallace.

           Description:
             All log space maps are addressed (ring, spoke).  All cartesian
             maps are addressed (col,row), i.e. (x, y).

           Conventions:
             i always indicates column no.
             j always indicates row no.
             TV images are always addressed (i, j).
             log map images are always addressed (ring, spoke)

----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           INCLUDE FILES
----------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <hipl_format.h>

#include "const.h"
#include "map.h"
#include "diag.h"
#include "image.h"
#include "proto.h"
/*----------------------------------------------------------------------
                           CONSTANTS
----------------------------------------------------------------------*/


/*----------------------------------------------------------------------
                           TYPEDEFS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           STATIC DATA
----------------------------------------------------------------------*/


/*----------------------------------------------------------------------
                           FUNCTION PROTOTYPES
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           FUNCTIONS
----------------------------------------------------------------------*/

IMAP *
MapIAlloc(int cols, int rows) {
  IMAP *map ;

  map = (IMAP *)calloc(1, sizeof(IMAP)) ;
  if (!map)
    return(NULL) ;

  map->rows = rows ;
  map->cols = cols ;
  map->image = (int *)calloc(rows * cols, sizeof(int)) ;
  if (!map->image) {
    free(map) ;
    return(NULL) ;
  }
  return(map) ;
}
/*----------------------------------------------------------------------
            Parameters:
        rows - the # of rows in the map.
        cols - the # of cols in the map.

           Description:
        initialize the runlength lookup tables.
----------------------------------------------------------------------*/
CMAP *
MapCAlloc(int cols, int rows) {
  CMAP *map ;

  map = (CMAP *)calloc(1, sizeof(CMAP)) ;
  if (!map)
    return(NULL) ;

  map->rows = rows ;
  map->cols = cols ;
  map->image = (unsigned char *)calloc(rows * cols, sizeof(unsigned char)) ;
  if (!map->image) {
    free(map) ;
    return(NULL) ;
  }
  return(map) ;
}
/*----------------------------------------------------------------------
            Parameters:
              map - the map to free

           Description:
              free a map.
----------------------------------------------------------------------*/
void
MapIFree(IMAP *map) {
  free(map->image) ;
  free(map) ;
}
/*----------------------------------------------------------------------
            Parameters:
              map - the map to free

           Description:
              free a map.
----------------------------------------------------------------------*/
void
MapCFree(CMAP *map) {
  free(map->image) ;
  free(map) ;
}

/*----------------------------------------------------------------------
            Parameters:
              map -   the map to write
              fname - the name of the file to write to.

           Description:
              write a map out to a file in hips format
----------------------------------------------------------------------*/
void
MapCHipsWrite(CMAP *map, char *fname) {
  IMAGE image ;

  ImageInitHeader(&image, map->cols, map->rows,  PFBYTE, 1) ;
  image.image = map->image ;

  ImageWrite(&image, fname) ;
}
/*----------------------------------------------------------------------
            Parameters:
              map -   the map to write
              fname - the name of the file to write to.

           Description:
              write a map out to a file in hips format
----------------------------------------------------------------------*/
void
MapCShowImage(CMAP *map, char *name) {
  IMAGE    image ;
#if USE_XWIN
  xwindow_type *xwin ;
#endif

  ImageInitHeader(&image, map->cols, map->rows,  PFBYTE, 1) ;
  image.image = map->image ;

#if USE_XWIN
  xwin = DebugShowImage(&image, NULL, 1, name, 1) ;
  readKeyboard(&image, xwin, 1) ;
  DebugEnd(xwin) ;
#endif
}
/*----------------------------------------------------------------------
            Parameters:
              map -   the map to write
              fname - the name of the file to write to.

           Description:
              write a map out to a file in hips format
----------------------------------------------------------------------*/
void
MapIShowImage(IMAP *map, char *name) {
  IMAGE    image ;
#if USE_XWIN
  xwindow_type *xwin ;
#endif

  ImageInitHeader(&image, map->cols, map->rows,  PFINT, 1) ;
  image.image = (unsigned char *)map->image ;

#if USE_XWIN
  xwin = DebugShowImage(&image, NULL, 1, name, 1) ;
  readKeyboard(&image, xwin, 1) ;
  DebugEnd(xwin) ;
#endif
}
/*----------------------------------------------------------------------
            Parameters:
              map -   the map to write
              fname - the name of the file to write to.

           Description:
              write a map out to a file in hips format
----------------------------------------------------------------------*/
void
MapIHipsWrite(IMAP *map, char *fname) {
  IMAGE image ;

  ImageInitHeader(&image, map->cols, map->rows, PFINT, 1) ;

  image.image = (unsigned char *)map->image ;
  ImageWrite(&image, fname) ;
}
/*----------------------------------------------------------------------
            Parameters:
              m - map to clone

           Description:
              allocate and make an exact copy of map m.
----------------------------------------------------------------------*/
CMAP *
MapClone(CMAP *m) {
  CMAP *mClone ;

  mClone = MapCAlloc(m->cols, m->rows) ;
  MapCopy(m, mClone) ;

  return(mClone) ;
}
/*----------------------------------------------------------------------
            Parameters:
              mSrc - map to copy from
              mDst - map to copy to

           Description:
              copy each cell in mSrc into mDst.
----------------------------------------------------------------------*/
void
MapCopy(CMAP *mSrc, CMAP *mDst) {
  int           size ;
  unsigned char *src, *dst ;

  size = mSrc->rows * mSrc->cols ;
  dst = mDst->image ;
  src = mSrc->image ;

  while (size--)
    *dst++ = *src++ ;
}
/*----------------------------------------------------------------------
            Parameters:
              m1   - input map
              m2   - input map
              mSum - output map

           Description:
              add each cell in maps m1 and m2 and put the result in map
              mSum.  All the maps must be the same size.
----------------------------------------------------------------------*/
void
MapAdd(CMAP *m1, CMAP *m2, CMAP *mSum) {
  int           size ;
  unsigned char *src1, *src2, *dst ;

  size = m1->rows * m1->cols ;

  src1 = m1->image ;
  src2 = m2->image ;
  dst = mSum->image ;

  while (size--)
    *dst++ = *src1++ + *src2++ ;
}
/*----------------------------------------------------------------------
            Parameters:
               mVal     - map of values.
               mVisited - map of points that are real values

           Description:
----------------------------------------------------------------------*/
#define DEBUG_ALL   0
#define MIN_CHANGE  1
#define MEAN_X      10
#define MEAN_Y      10

void
MapSoapBubbleRelax(IMAP *mVal, CMAP *mVisited) {
  int    x, y, kx, ky, count, max_change, sum, xmin, xmax, ymin, ymax, val,
  *cell, trialno ;
  UCHAR  *vcell ;

#if DEBUG_ALL
  Gdiag = 1 ;
  printf("0:\n") ;
  MapIPrint(mVal, stdout) ;
#endif

#if 1
  xmin = mVisited->cols ;
  ymin = mVisited->rows ;
  xmax = ymax = 0 ;
  for (x = 0 ; x < mVisited->cols ; x++) {
    for (y = 0 ; y < mVisited->rows ; y++) {
      if (MAPcell(mVisited, x, y)) {
        if (x < xmin)
          xmin = x ;
        if (x > xmax)
          xmax = x ;
        if (y < ymin)
          ymin = y ;
        if (y > ymax)
          ymax = y ;
      }
    }
  }
#else

  ymin = xmin = 0 ;
  ymax = mVisited->rows - 1 ;
  xmax = mVisited->cols - 1 ;
#endif

  /* first run 10 by 10 mean filter on fixed points to speed stuff up */
  for (y = ymin ; y <= ymax ; y++) {
    for (x = xmin ; x <= xmax ; x++) {
      if (MAPcell(mVisited, x, y))
        continue ;

      sum = 0 ;
      count = 0 ;
      for (ky = y - MEAN_Y/2 ; ky <= y + MEAN_Y/2 ; ky++) {
        if ((ky >= ymin) && (ky <= ymax))
          for (kx = x-MEAN_X/2 ; kx <= x+MEAN_X/2 ; kx++) {
            if ((kx >= xmin) && (kx <= xmax) && MAPcell(mVisited, kx, ky)) {
              count++ ;
              cell = &MAPcell(mVal, kx, ky) ;
              sum += *cell ;
            }
          }
      }
      if (count)   /* sum fixed points in range */
      {
        sum /= count ;
        MAPcell(mVal, x, y) = sum ;
      }
    }
  }

  trialno = 0 ;

  do {
    max_change = 0 ;
    for (y = ymin ; y <= ymax ; y++) {
      for (x = xmin ; x <= xmax ; x++) {
        if (MAPcell(mVisited, x, y))
          continue ;                   /* fix points that were given */

        count = sum = 0 ;
        for (ky = y-1 ; ky <= y+1 ; ky++) {
          if ((ky >= ymin) && (ky <= ymax)) for (kx = x-1 ; kx <= x+1 ;kx++) {
              if ((kx >= xmin) && (kx <= xmax)) {
                count++ ;
                cell = &MAPcell(mVal, kx, ky) ;
                if (*cell > 200) {
                  vcell = &MAPcell(mVisited, x, y) ;
                }

                sum += MAPcell(mVal, kx, ky) ;
              }
            }
        }

        if (!count)
          continue ;

        val = MAPcell(mVal, x, y) ;
        sum = nint((double)sum / (double)count) ;
        if (abs(val - sum) > max_change)
          max_change = abs(val - sum) ;

        MAPcell(mVal, x, y) = sum ;
      }
    }
#if DEBUG_ALL
    printf("%d: delta %d\n", ++trialno, max_change) ;
    MapIPrint(mVal, stdout) ;
#endif
  } while (max_change > MIN_CHANGE) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
MapIPrint(IMAP *map, FILE *fp) {
  int  x, y ;

  if (!Gdiag) return ;

  for (y = 0 ; y < map->rows ; y++) {
    for (x = 0 ; x < map->cols ; x++) {
      fprintf(fp, "%4.4d  ", MAPcell(map, x, y)) ;
    }
    fprintf(fp, "\n") ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
MapCPrint(CMAP *map, FILE *fp) {
  int  x, y ;

  if (!Gdiag) return ;

  for (y = 0 ; y < map->rows ; y++) {
    for (x = 0 ; x < map->cols ; x++) {
      fprintf(fp, "%4.4d  ", (int)MAPcell(map, x, y)) ;
    }
    fprintf(fp, "\n") ;
  }
}
