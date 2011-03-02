/**
 * @file  connect.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:11 $
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


/*
  @(#)connect.c 1.1
  3/31/94
*/
/*----------------------------------------------------------------------
           File Name:
             connect.c

           Author:

           Description:
             Routines for building lists of objects with connected
             components.

           Conventions:

----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           INCLUDE FILES
----------------------------------------------------------------------*/

#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <hipl_format.h>

#include "image.h"
#include "h_logz.h"
#include "diag.h"
#include "connect.h"
#include "macros.h"
#include "const.h"
#include "error.h"
#include "proto.h"

/*----------------------------------------------------------------------
                           FUNCTION PROTOTYPES
----------------------------------------------------------------------*/

static int        badRle(RLE *rle) ;
static int        checkRow(RLE *rle) ;
static void       isFree(COBJ *cobj) ;
static void       cTableFree(COBJ_TABLE *ctable) ;
static COBJ_TABLE *cTableAlloc(int max_objects, int max_rows,
                               IMAGE *image, int xo, int yo) ;
static void        CobjAdd(COBJ *cobj, int ring, int spoke) ;
#if USE_XWIN
static void       showTable(COBJ_TABLE *ctable) ;
#endif
static void        addObjectToChain(COBJ *child, COBJ *parent) ;
static void       mergeObjectChain(COBJ *cobj) ;
static void       mergeObjects(COBJ *src, COBJ *dst) ;
static COBJ       *newObject(COBJ_TABLE *ctable, int ring, int spoke) ;
static int        inChain(COBJ *c1, COBJ *c2) ;
void ClipUnion(CLIP *c1, CLIP *c2, CLIP *cdst) ;
void ClipCopy(CLIP *cSrc, CLIP *cDst) ;
int  CobjOverlap(COBJ *cobj, COBJ_TABLE *ct, int dx, int dy) ;
int  ClipIntersect(CLIP *c1, CLIP *c2) ;

void checkObj(COBJ *cobj) ;
static void dumpTable(COBJ_TABLE *ctable) ;
static void cobjClear(COBJ *cobj) ;
static void cobjClear(COBJ *cobj) ;
static int cobjFree(COBJ *cobj) ;
void ConnTableCountRles(COBJ_TABLE *ctable) ;

/*----------------------------------------------------------------------
                                CONSTANTS
----------------------------------------------------------------------*/

#define DIAG_SCALE      10
#define MAX_OBJECTS     1000
#define MAX_ROWS        100
#define COBJ_ROWNO(c,y) (y - c->starty)
#define MAX_DIST        11
#define MAX_AREA        600
#define MIN_AREA        10
#define RADIUS_FUDGE    3.5

#define DIAG_MATCH             0
#define DIAG_SEGMENTS          0
#define DIAG_SEGMENTS_VERBOSE  0
#define DIAG_MERGE             0

/*----------------------------------------------------------------------
                                FUNCTIONS
---------------------------------------------------------------------*/

int RLES = 0, ROWS = 0 ;

static void
isFree(COBJ *cobj) {
  int rowno ;
  RLE *rle ;

  for (rowno = 0 ; rowno < cobj->nrows ; rowno++) {
    rle = cobj->rows[rowno] ;
    if (rle)
      badRle(rle) ;
  }
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
             find all connected black (or white) components in an image.
----------------------------------------------------------------------*/
COBJ_TABLE *
ConnectedComponents(LOGMAP_INFO *lmi, IMAGE *image, int white, int thresh,
                    int xo, int yo) {
  int         ring, nobjects, spoke, nrings, nspokes, color, objno, freeno,
  tsign ;
  COBJ_TABLE *ctable ;
  LOGPIX     *logpix, *nlogpix ;
  UINT       *pix ;
  NEIGHBOR   *neighbor ;
  COBJ       *cobj = NULL, *freeObj ;

  if (image->pixel_format != PFINT)
    ErrorExit(ERROR_BADPARM,
              "ConnectedComponents: bad pixel format %d (must by int)\n",
              image->pixel_format) ;

#if USE_XWIN
  if (debug & DIAG_SEGMENT_VERBOSE)
    DebugNewWindow(image->cols, image->rows, "segmentation", DIAG_SCALE);
#endif

#if 0
  if (Gdiag & DIAG_TIMER)
    DebugTimerStart() ;
#endif
  ctable = cTableAlloc(MAX_OBJECTS, MAX_ROWS, image, xo, yo) ;
  nrings = lmi->nrings ;
  nspokes = lmi->nspokes ;
  logpix = LOG_PIX(lmi, 0, 0) ;
  pix = IMAGEIpix(image, 0, 0) ;
  ctable->white = white ;
  ctable->thresh = thresh ;
  ctable->lmi = lmi ;

  if (white)             /* looking for white pixels, anything above thresh */
    tsign = 1 ;
  else                  /* looking for black pixels, anything above -thresh */
  {
    thresh *= -1 ;
    tsign = -1 ;
  }

  for_all_log_pixels(lmi, ring, spoke) {
    logpix = LOG_PIX(lmi, ring, spoke) ;
    logpix->user = NULL ;               /* reset - not part of an object yet */
    pix = IMAGEIpix(image, ring, spoke) ;
    if (logpix->area > 0) {
      if (*pix * tsign > thresh) {
        cobj = NULL ;      /* object that this pixel belongs to */
        color = WHITE ;
        for_each_neighbor(lmi, neighbor, ring, spoke) {
          nlogpix = neighbor->logpix ;

          /* only consider neighboring pixels that have already been visited */
          if (nlogpix->spoke > spoke)
            break ;
          else if ((nlogpix->spoke == spoke) && (nlogpix->ring > ring))
            break ;

          if (nlogpix->user)    /* previous pixel has already been assigned */
          {
            /*
              if cobj is NULL this pixel has not yet been assigned to an object.  In
              that case, we assign it to the first neighboring object.

              if cobj is not NULL and it is different from the object that one of our
              neighbors is part of, we must mark these two objects for merging, as
              they are part of the same object.  This can happen for objects such
              as a 'V', where the two arms will be classified as separate objects
              until the base of the 'V' is reached, and the two arms must be merged.
            */
            if (!cobj)
              cobj = (COBJ *)nlogpix->user;
            else if (cobj != nlogpix->user)  /* bordering different objects */
            {
              addObjectToChain((COBJ *)nlogpix->user, cobj) ;
            }
          }
        }

        /* if no neighboring pixels are assigned, allocate new object */
        if (!cobj)
          cobj = newObject(ctable, ring, spoke) ;

        logpix->user = (void *)cobj ;
        CobjAdd(cobj, ring, spoke) ;
      } else color = BLACK ;
    } else color = BLACK ;

#if USE_XWIN
    if (Gdiag & DIAG_SEGMENT_VERBOSE) {
      DebugDrawPoint(ring, spoke, color) ;
      DebugFlush() ;
    }
#endif

    logpix++ ;
    pix++ ;
  }

  /* now merge objects together that are adjacent */
  nobjects = ctable->nobjects ;
  ctable->nobjects = 0 ;
  for (objno = 0 ; objno < nobjects ; objno++) {
    cobj = ctable->objects + objno ;
    if (cobj->used) {
      mergeObjectChain(cobj) ;
      ctable->nobjects++ ;
    }
  }

  /*
    this compacting can be done faster, but for now....
  */

  /*
    compact the table.  Note that at this point, nobjects represents the
    total # of original objects found, while ctable->nobjects represents
    the # of objects after merging.
  */
  for (freeno = 0 ; freeno < ctable->nobjects ; freeno++) {
    freeObj = ctable->objects + freeno ;
    if (freeObj->used) continue ;

    /*
      found a free object, now search forward for one in use to copy into
      the free one.
    */
    for (objno = freeno + 1 ; objno < nobjects ; objno++) {
      cobj = ctable->objects + objno ;
      if (cobj->used) break ;
    }

    if (cobj->used)    /* copy used object into free one */
    {
      isFree(freeObj) ;
      CobjCopy(cobj, freeObj, 0, 0, 1) ;
      cobjFree(cobj) ;
      isFree(cobj) ;
      freeObj->used = 1 ;
    }
  }

  ConnTableCalculateCentroids(ctable) ;
  /*
    compact the table again, getting rid of objects that are too big.
    Note that at this point, nobjects represents the
    total # of original objects found, while ctable->nobjects represents
    the # of objects after merging.
  */
  for (freeno = 0 ; freeno < ctable->nobjects ; freeno++) {
    freeObj = ctable->objects + freeno ;

    if ((freeObj->tv_area > MAX_AREA) || (freeObj->tv_area < MIN_AREA))
      cobjFree(freeObj) ;

    /* keep objects in the right area range */
    if (freeObj->used)
      continue ;

    /*
      found a free object, now search forward for one in use to copy into
      the free one.
    */
    for (objno = freeno + 1 ; objno < ctable->nobjects ; objno++) {
      cobj = ctable->objects + objno ;
      if (cobj->used)
        break ;
    }

    /* copy used object into free one */
    if ((objno < ctable->nobjects) && cobj->used) {
      isFree(freeObj) ;
      CobjCopy(cobj, freeObj, 0, 0, 1) ;
      cobjFree(cobj) ;
      isFree(cobj) ;

      freeObj->used = 1 ;
      freeno-- ;               /* force checking of this object */
    } else {
      freeObj->used = 0 ;   /* last one in table */
      isFree(freeObj) ;
    }
    ctable->nobjects-- ;       /* one less object in table */
  }

#if 0
  if (Gdiag & DIAG_TIMER)
    DebugTimerShow("segmentation") ;
#endif

#if USE_XWIN
  if (Gdiag & DIAG_SEGMENT_VERBOSE) {
    DebugKeyboard() ;
    DebugEndWindow() ;
  }
  if (Gdiag & DIAG_SEGMENT)
    showTable(ctable) ;
#endif

  return(ctable) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             merge the list of components belonging to 'co1' into those
             belonging to 'co2'.
----------------------------------------------------------------------*/
void
ConnTableFree(COBJ_TABLE **pctable) {
  COBJ_TABLE *ctable ;

  ctable = *pctable ;
  cTableFree(ctable) ;
  *pctable = NULL ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static COBJ_TABLE *
cTableAlloc(int max_objects, int max_rows, IMAGE *image, int xo, int yo) {
  COBJ_TABLE *ctable ;
  COBJ       *object ;
  int        i ;

  ctable = (COBJ_TABLE *)calloc(1, sizeof(COBJ_TABLE)) ;
  if (!ctable)
    ErrorExit(ERROR_NO_MEMORY, "CtableAlloc: could not allocate ctable\n") ;

  ctable->objects = (COBJ *)calloc(max_objects, sizeof(COBJ)) ;
  if (!ctable->objects)
    ErrorExit(ERROR_NO_MEMORY,"CtableAlloc: could not allocate %d objects\n", max_objects) ;
  ctable->mapped = (int *)calloc(max_objects, sizeof(int)) ;
  if (!ctable->mapped)
    ErrorExit(ERROR_NO_MEMORY,"CtableAlloc: could not allocate %d mappings\n", max_objects) ;

  ctable->max_objects = max_objects ;
  ctable->image = image ;
  ctable->nrows = max_rows ;
  ctable->xo = xo ;
  ctable->yo = yo ;

  /* initialize all object ids */
  for (i = 0 ; i < max_objects ; i++) {
    object = &ctable->objects[i] ;
    object->id = i ;
    object->starty = 0 ;
    object->nrows = max_rows ;
    object->rows = (RLE **)calloc(max_rows, sizeof(RLE *)) ;
    ROWS++ ;
    if (!object->rows)
      ErrorExit(ERROR_NO_MEMORY,"CtableAlloc: could not allocate %d rows for %dth object\n",
                max_rows, i) ;
  }

  return(ctable) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
cTableFree(COBJ_TABLE *ctable) {
  int   i, rowno, count ;
  RLE  *rle, *next ;
  COBJ *cobj ;

  ConnTableCountRles(ctable) ;
  for (i = 0 ; i < ctable->max_objects ; i++) {
    cobj = &ctable->objects[i] ;
#if 0
    if (i < ctable->nobjects)
      fprintf(stderr, "object %d:\n", i) ;
#endif
    for (rowno = 0 ; rowno < cobj->nrows ; rowno++) {
      rle = cobj->rows[rowno] ;
      checkRow(rle) ;
      count = 0 ;
      while (rle) {
        next = rle->next ;
        free(rle) ;
        RLES-- ;
        count++ ;
        rle = next ;
      }
#if 0
      if (count) {
        fprintf(stderr, "row %d: %d freed\n", rowno, count) ;
        if (i >= ctable->nobjects)
          fprintf(stderr, "freeing RLE in row %d, nobjects %d!!\n",
                  i, ctable->nobjects) ;
      }
#endif
    }
    free(cobj->rows) ;
    ROWS-- ;
  }

  free(ctable->objects) ;
  free(ctable->mapped) ;
  free(ctable) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
CobjAdd(COBJ *cobj, int ring, int spoke) {
  int  end, added, rowno, which ;
  RLE  *new, *rle, *prev, *first,sfirst ;
  CLIP newClip ;

  /* update the bounding box */
  added = 1 ;                        /* assume pixel is a new one */
  newClip.x = ring ;
  newClip.y = spoke ;
  newClip.dx = newClip.dy = 1 ;
  ClipUnion(&newClip, &cobj->clip, &cobj->clip) ;

  /* add this pixel to the appropriate place in the row rle list */
  rowno = COBJ_ROWNO(cobj, spoke) ;
  if ((rowno < 0) || (rowno >= cobj->nrows)) {
    fprintf(stderr, "CobjAdd(%d, %d) --> row out of bounds <%d,%d>\n",
            ring, spoke, cobj->starty, cobj->starty+cobj->nrows-1) ;
    return ;
  }
  first = rle = cobj->rows[rowno] ; /* DEBUG */
  if (first)
    sfirst = *first ; /* DEBUG */

  checkRow(rle) ;
  new = NULL ;
  prev = NULL ;
  while (rle && rle->start <= ring) {
    prev = rle ;
    rle = rle->next ;
  }

  /*
    prev is non-null it will point to the rle whose 'start' is before
    the current pixel, and closest to it.
  */
  if (prev) {
    which = 1 ;
    end = prev->start + prev->len - 1 ;
    if (ring <= end)               /* pixel was already in object */
      added = 0 ;
    else if (ring == (end + 1))    /* add to end of run */
      prev->len++ ;
    else                           /* new run, pixel hasn't been added yet */
    {
      which = 2 ;
      new = rle = (RLE *)calloc(1, sizeof(RLE)) ;
      RLES++ ;
      rle->next = prev->next ;
      prev->next = rle ;
      rle->start = ring ;
      rle->len = 1 ;
    }
  } else    /* add to start of list */
  {
    which = 3 ;
    new = rle = (RLE *)calloc(1, sizeof(RLE)) ;
    RLES++ ;
    rle->next = cobj->rows[rowno] ;
    cobj->rows[rowno] = rle ;
    rle->start = ring ;
    rle->len = 1 ;
  }

  checkRow(cobj->rows[rowno]) ;

  if (added)
    cobj->log_area++ ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#if USE_XWIN
static void
showTable(COBJ_TABLE *ctable) {
  IMAGE *image ;
  COBJ      *object ;
  int       i, row, x, y, col, xo, yo ;
  RLE       *rle ;

  image = ctable->image ;
  xo = ctable->xo ;
  yo = ctable->yo ;

  DebugNewWindow(image->cols, image->rows, "segments", DIAG_SCALE);

  for (i = 0 ; i < ctable->nobjects ; i++) {
    object = ctable->objects + i ;

    printf("segment #%d: (%d, %d) --> (%d, %d), log area %d, "
           "tv area %d, centroid (%d, %d)\n",
           i, object->clip.x, object->clip.y, object->clip.dx,object->clip.dy,
           object->log_area, object->tv_area, object->col_cent,
           object->row_cent) ;

    for (row = 0 ; row < ctable->nrows ; row++) {
      /* go through each run in this object */
      for (rle = object->rows[row] ; rle ; rle = rle->next) {
        for (col = rle->start ; col < rle->start + rle->len ; col++)
          DebugDrawPoint(col - xo, row - yo, WHITE) ;
      }
    }

    DebugFlush() ;
    if (i < ctable->nobjects-1)
      getchar() ;
  }

  DebugKeyboard() ;
  DebugEndWindow() ;
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
dumpTable(COBJ_TABLE *ctable) {
  COBJ      *object ;
  int       i ;

  printf("ctable at %d,%d, with %d objects\n", ctable->xo, ctable->yo,
         ctable->nobjects) ;
  for (i = 0 ; i < ctable->nobjects ; i++) {
    object = ctable->objects + i ;

    printf("segment #%d: (%d, %d) --> (%d, %d), log area %d, center %d, %d, "
           "tv area %d, centroid (%d, %d)\n",
           i, object->clip.x, object->clip.y, object->clip.dx,object->clip.dy,
           object->log_area, object->log_col_cent, object->log_row_cent,
           object->tv_area, object->col_cent, object->row_cent) ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
addObjectToChain(COBJ *child, COBJ *parent) {
  COBJ   *firstChild, *lastChild ;

  if (inChain(child, parent))
    return ;                    /* already linked, don't do it again */

  /*
    find the last child in the parent's current chain so that we can
    add 'child' to the end of the list.
  */
for (lastChild = parent ; lastChild->child ; lastChild = lastChild->child) {}

  /*
    Now find the head of the child's parent chain, this will be the
    object that is linked to the end of the parent's chain.
  */
  for (firstChild = child; firstChild->parent; firstChild = firstChild->parent) {}

  lastChild->child = firstChild ;
  firstChild->parent = lastChild ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              merge all objects attached to this one into it.  This
              includes objects above (parents) and below (children)
              it in the chain.  Mark all objects that were merged as
              empty.
----------------------------------------------------------------------*/
static void
mergeObjectChain(COBJ *cobj) {
  COBJ  *merge, *child ;

  /* find root of parent chain */
  for (merge = cobj ; merge->parent ; merge = merge->parent) {}

  /* now travel down the child chain, merging objects */
  for (  child = merge ; merge ; merge = child) {
    child = merge->child ;       /* merge->child won't be valid after merge */
    if (merge != cobj)           /* don't merge object into itself */
    {
      COBJ *next ;

      next = merge + 1 ;
      mergeObjects(merge, cobj) ;
      isFree(merge) ;
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              merge 'src' into 'dst'.  Note that 'src' will be taken out
              of the chain it is in.
----------------------------------------------------------------------*/
static void
mergeObjects(COBJ *src, COBJ *dst) {
  int row, colOffset, y ;
  RLE *rle, *next ;

  src->used = 0 ;              /* no longer in use */

  /* clip is box union of two clips */
  ClipUnion(&src->clip, &dst->clip, &dst->clip) ;

  for (row = 0 ; row < src->nrows ; row++) {
    rle = src->rows[row] ;
    y = row + src->starty ;
    if (rle)                  /* this row has something in it */
    {
      /* add each pixel in the source row to the destination object */
      while (rle) {
        for (colOffset = 0 ; colOffset < rle->len ; colOffset++)
          CobjAdd(dst, rle->start + colOffset, y) ;

        /* now that the rle has been merged, we can free it */
        next = rle->next ;
        free(rle) ;
        RLES-- ;
        rle = next ;
      }

#if 0
      /* I don't think this can ever happen */
      /* now merge adjacent rle's in the destination row that overlap */
      rle = dst->rows[row] ;
      next = rle->next ;
      while (next) {
        rleEnd = rle->start + rle->len - 1 ;
        if (rleEnd >= next->start)  /* there is overlap, merge them */
        {
          nextEnd = next->start + next->len - 1 ;
          rle->len = MAX(nextEnd, rleEnd) - rle->start + 1 ;

          /* now remove 'next' from the chain */
          rle->next = next->next ;
          free(next) ;
        } else rle = next ;
        next = rle->next ;
      }
#endif
      src->rows[row] = NULL ;
    }
  }

  /* remove src from the link-list */
  if (src->parent) src->parent->child = src->child ;
  if (src->child) src->child->parent = src->parent ;
  src->child = src->parent = NULL ;
  isFree(src) ;
}
void
ClipCopy(CLIP *cSrc, CLIP *cDst) {
  *cDst = *cSrc ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             calculate the box union of 'c1' and 'c2', storing the
             result in 'cdst'.
----------------------------------------------------------------------*/
void
ClipUnion(CLIP *c1, CLIP *c2, CLIP *cdst) {
  int x0, y0, xend, yend, xend1, yend1, xend2, yend2 ;

  x0 = MIN(c1->x, c2->x) ;
  y0 = MIN(c1->y, c2->y) ;
  yend1 = c1->y + c1->dy - 1 ;
  xend1 = c1->x + c1->dx - 1 ;
  xend2 = c2->x + c2->dx - 1 ;
  yend2 = c2->y + c2->dy - 1 ;

  xend = MAX(xend1, xend2) ;
  yend = MAX(yend1, yend2) ;

  cdst->x = x0 ;
  cdst->y = y0 ;
  cdst->dx = xend - x0 + 1 ;
  cdst->dy = yend - y0 + 1 ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              Allocate a new (empty) object.
----------------------------------------------------------------------*/
static COBJ *
newObject(COBJ_TABLE *ctable, int ring, int spoke) {
  COBJ  *cobj ;

  if (ctable->nobjects >= MAX_OBJECTS)
    ErrorExit(ERROR_NO_MEMORY,"newObject: too many objects\n") ;
  cobj = &ctable->objects[ctable->nobjects++] ;
  cobjClear(cobj) ;
  cobj->used = 1 ;
  cobj->clip.x = ring ;
  cobj->clip.y = spoke ;
  cobj->clip.dx = cobj->clip.dy = 1 ;
  /*
    assume we are at the middle of the object, so have as many rows
    above as below.
  */
  cobj->starty = spoke - (cobj->nrows/2) ;
  return(cobj) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             calculate the centroid of each object in the table.
----------------------------------------------------------------------*/
void
ConnTableCalculateCentroids(COBJ_TABLE *ctable) {
  int           objno, row, x, endx, weight, total_weight,
  pix_val, white, amount_on, thresh, all_on, y, xo, yo ;
  LOGMAP_INFO   *lmi ;
  COBJ          *cobj ;
  RLE           *rle ;
  LOGPIX        *logpix ;
  IMAGE     *image ;

  xo = ctable->xo ;
  yo = ctable->yo ;
  thresh = ctable->thresh + 1 ;     /* weight should never be 0 */
  image = ctable->image ;
  white = ctable->white ;
  lmi = ctable->lmi ;
  if (white)
    all_on = 255 ;       /* every tv pixel in log pix is 'on' */
  else
    all_on = thresh ;    /* every tv pixel in log pix is 'on' */

  for (cobj = ctable->objects, objno = 0 ; objno < ctable->nobjects ;
       objno++, cobj++) {
    total_weight = 0 ;
    for (row = 0 ; row < cobj->nrows ; row++) {
      y = row + cobj->starty ;
      for (rle = cobj->rows[row] ; rle ; rle = rle->next) {
        x = rle->start ;
        logpix = LOG_PIX(lmi, x, y) ;

        /* x and y are in image/logpix coordinates */
        for (endx = x + rle->len - 1 ; x <= endx ; x++, logpix++) {
          pix_val = *IMAGEIpix(image, x, y) ;
          if (white)
            amount_on = pix_val - thresh ;
          else
            amount_on = thresh - pix_val ;

          weight = amount_on * logpix->area ;

          total_weight += weight ;

          /* centroids should be in world coordinates */
          cobj->col_cent += weight * (logpix->col_cent) ;
          cobj->row_cent += weight * (logpix->row_cent) ;
          cobj->log_col_cent += logpix->ring ;
          cobj->log_row_cent += logpix->spoke ;
          cobj->tv_area += weight ;
        }
      }
    }

    cobj->tv_area /= all_on ;
    cobj->row_cent /= total_weight ;
    cobj->col_cent /= total_weight ;
    cobj->row_cent += yo ;
    cobj->col_cent += xo ;
    cobj->log_col_cent /= cobj->log_area ;
    cobj->log_row_cent /= cobj->log_area ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             Determine if c1 is already in the child or parent chain
             of c2.  If so, return 1, otherwise return 0.
----------------------------------------------------------------------*/
static int
inChain(COBJ *c1, COBJ *c2) {
  COBJ *cobj ;

  /* find root of c1's parent chain */
  for (cobj = c1 ; cobj ; cobj = cobj->parent) {
    if (cobj == c2)      /* c2 is one of c1's ancestors */
      return(1) ;
  }

  /* find end of c1's child chain */
  for (cobj = c1 ; cobj ; cobj = cobj->child) {
    if (cobj == c2)      /* c2 is one of c1's children */
      return(1) ;
  }
  return(0) ;            /* not found */
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             find the best fit match for the offset of table ct2 from
             ct1.  This routine is a real hack...

           Coordinate System Note:
             the values returned from this routine will transform a
             point in the space of ct2 into ct1.  This uses 'normal'
             (for image processing) y values, which increase as they
             get to the bottom of the image.
----------------------------------------------------------------------*/
#define MIN_AREA_RATIO    0.0
#define MAX_AREA_RATIO    4.0
#define MIN_POINTS        4

int
ConnTableMatch(LOGMAP_INFO *lmi, COBJ_TABLE *ct1, COBJ_TABLE *ct2,
               int *px, int *py, int dx_old, int dy_old, int max_change) {
  int     objno, dx, dy, dist, total_dist, min_error, min_anchor1, min_anchor2,
  i, c2index, cx2, cy2, npoints, xdiff, ydiff, avg_dist,
  min_area, max_area, error, dxTotal, dyTotal, anchorNo, var_x, var_y,
  std_dev_x, std_dev_y, nbad, ring ;
  COBJ    *cobj1, *cobj2, *cAnchor1, *cAnchor2 ;
  double  point_match_factor ;

  min_anchor1 = min_anchor2 = -1 ;
  min_error = -1 ;

  /*
    now try and match the anchor object with each object in ct2 space,
    finding the best match of the other objects.  If the match is good
    enough, return.
  */
#if 0
  if (Gdiag & DIAG_MATCH)
    printf("looking for at least %d points\n", MIN_POINTS) ;
#endif
  for (anchorNo = 0 ; anchorNo < ct1->nobjects ; anchorNo++) {
    cAnchor1 = &ct1->objects[anchorNo] ;
    for (objno = 0 ; objno < ct2->nobjects ; objno++) {
      cAnchor2= &ct2->objects[objno] ;
      npoints = 1 ;    /* # of points matched (1 for anchor) */

      memset((char *)ct1->mapped, -1, sizeof(*ct1->mapped) * ct1->nobjects) ;
      memset((char *)ct2->mapped, -1, sizeof(*ct2->mapped) * ct2->nobjects) ;

      /* link the two objects */
      ct1->mapped[anchorNo] = objno ;
      ct2->mapped[objno] = anchorNo ;

      dx = cAnchor1->col_cent - cAnchor2->col_cent ;
      dy = cAnchor1->row_cent - cAnchor2->row_cent ;

      if ((max_change != NO_CONTINUITY) &&
          ((abs(dx_old - dx) > max_change) || (abs(dy_old - dy) > max_change)))
        continue ;

      /*
        now run through every object in ct1 and find the nearest match
        in ct2, keeping track of the summed distance for all objects.
      */
      total_dist = 0 ;
      for (i = anchorNo ; i < ct1->nobjects ; i++) {
        cobj1 = &ct1->objects[i] ;
        min_area = nint((double)cobj1->tv_area * MIN_AREA_RATIO) ;
        max_area = MIN(nint((double)cobj1->tv_area * MAX_AREA_RATIO),MAX_AREA);

        cobj2 =
          CobjFindMatch(lmi, ct2, cobj1->col_cent - dx, cobj1->row_cent - dy,
                        1, &dist, min_area, max_area, MAX_DIST);

        if (cobj2)   /* found an object */
        {
          total_dist += dist ;
          npoints++ ;
          c2index = cobj2 - ct2->objects ;
          ct2->mapped[c2index] = i ;
          ct1->mapped[i] = c2index ;
        }
      }


      avg_dist = total_dist / npoints ;

      /* matches which account for more points should be favored over
         matches with smaller number of points.
         */
      point_match_factor = (double)ct1->nobjects / (double)npoints ;
      point_match_factor *= (point_match_factor * point_match_factor) ;
      error = nint((double)avg_dist * point_match_factor) ;

#if 0
      if ((Gdiag & DIAG_MATCH) && (npoints > MIN_POINTS-2))
        fprintf(stderr,
                "anchor %d: (%d, %d) --> %d: (%d, %d), avg_dist %d, "
                "pmf - %2.3lf (%d), error %d\n",
                anchorNo, cAnchor1->col_cent, cAnchor1->row_cent,
                objno, cAnchor2->col_cent, cAnchor2->row_cent,
                avg_dist, point_match_factor, npoints, error) ;
#endif

      if (npoints < MIN_POINTS)   /* not enough points for a reliable match */
        continue ;

      if ((min_error < 0) || (error < min_error)) {
        min_error = error ;
        min_anchor1 = anchorNo ;
        min_anchor2 = objno ;
        *px = dx ;
        *py = dy ;
      }
    }
  }

  if (min_error < 0)
    return(-1) ;              /* no reliable match */

#if 0
  if (Gdiag & DIAG_MATCH)
    fprintf(stderr, "min_anchor1 %d, min_anchor2 %d, min_error %d\n",
            min_anchor1, min_anchor2, min_error) ;
#endif

  /* redo mappings for smallest error and calculate average deltas */
  dx = *px ;
  dy = *py ;
  memset((char *)ct1->mapped, -1, sizeof(*ct1->mapped) * ct1->nobjects) ;
  memset((char *)ct2->mapped, -1, sizeof(*ct2->mapped) * ct2->nobjects) ;
  cAnchor1 = &ct1->objects[min_anchor1] ;
  cAnchor2 = &ct2->objects[min_anchor2] ;
  ct1->mapped[min_anchor1] = min_anchor2 ;
  ct2->mapped[min_anchor2] = min_anchor1 ;
  npoints = 1 ;

  dxTotal = cAnchor1->col_cent - cAnchor2->col_cent ;
  dyTotal = cAnchor1->row_cent - cAnchor2->row_cent ;

  for (i = 0 ; i < ct1->nobjects ; i++) {
    if (i == min_anchor1)
      continue ;

    cobj1 = &ct1->objects[i] ;

    min_area = nint((double)cobj1->tv_area * MIN_AREA_RATIO) ;
    max_area = MIN(nint((double)cobj1->tv_area * MAX_AREA_RATIO), MAX_AREA) ;
    cobj2 =
      CobjFindMatch(lmi, ct2, cobj1->col_cent - dx, cobj1->row_cent - dy, 1,
                    &dist, min_area, max_area, MAX_DIST) ;
    if (cobj2)   /* found an object */
    {
      c2index = cobj2 - ct2->objects ;
      ct2->mapped[c2index] = i ;
      ct1->mapped[i] = c2index ;

      dxTotal += (cobj1->col_cent - cobj2->col_cent) ;
      dyTotal += (cobj1->row_cent - cobj2->row_cent) ;
      npoints++ ;
    }
  }

  if (npoints < MIN_POINTS)   /* not enough points for a reliable match */
    return(-1) ;

  /* At this point, a valid match has been found.  Now go about the calculation
     of the deltas between the two tables.
     use difference between center of gravitys of the two 'objects' defined
     by the segment tables as the deltas between the two images.
     */
  dx = dxTotal / npoints ;
  dy = dyTotal / npoints ;

  /*
    calculate the variance and the standard deviation to use in discarding
    outlying points.
  */
  var_x = var_y = 0 ;

  for (i = 0 ; i < ct2->nobjects ; i++) {
    if (ct2->mapped[i] >= 0) {
      cobj2 = &ct2->objects[i] ;
      cobj1 = &ct1->objects[ct2->mapped[i]] ;
      xdiff = cobj1->col_cent - cobj2->col_cent ;
      ydiff = cobj1->row_cent - cobj2->row_cent ;
      var_x += ((xdiff - dx) * (xdiff - dx)) ;
      var_y += ((ydiff - dy) * (ydiff - dy)) ;
    }
  }

  if (((double)var_x / (double)npoints) < 0)
    fprintf(stderr, "varx / npoints - %f (%d) / %f < 0!!\n",
            (double)var_x, var_x, (double)npoints) ;
  if (((double)var_y / (double)npoints) < 0)
    fprintf(stderr, "vary / npoints - %f / %f < 0!!\n",
            (double)var_y, (double)npoints) ;

  std_dev_x = nint(sqrt((double)var_x) / (double)npoints) ;
  std_dev_y = nint(sqrt((double)var_y) / (double)npoints) ;
#if 0
  if (Gdiag & DIAG_MATCH)
    fprintf(stderr, "standard deviation (%d, %d), variance (%d, %d)\n",
            std_dev_x, std_dev_y, var_x, var_y) ;
#endif

  /*
    discard points more than 2 standard deviations away for purposes of
    calculating accurate deltas.  This will not affect whether or not
    a match occurs.  Points are only thrown out as long as sufficient
    points to accurately match remain.
  */
#define DISCARD_POINTS 0
#if DISCARD_POINTS
  for (i = 0 ; i < ct2->nobjects ; i++) {
    if (ct2->mapped[i] >= 0) {
      cobj2 = &ct2->objects[i] ;
      cobj1 = &ct1->objects[ct2->mapped[i]] ;
      xdiff = abs(cobj1->col_cent - (cobj2->col_cent + dx)) ;
      ydiff = abs(cobj1->row_cent - (cobj2->row_cent + dy)) ;
      if (Gdiag & DIAG_MATCH)
        fprintf(stderr,
                "%d (%d, %d) --> %d (%d, %d), delta %d, %d\n",
                i, cobj1->col_cent, cobj1->row_cent,
                ct1->mapped[i], cobj2->col_cent, cobj2->row_cent,
                xdiff, ydiff) ;

      if ((npoints > MIN_POINTS) &&
          ((xdiff > 2 * std_dev_x) || (ydiff > 2 * std_dev_y))) {
        if (Gdiag & DIAG_MATCH)
          fprintf(stderr, "discarding...\n") ;
        npoints-- ;
        ct2->mapped[ct1->mapped[i]] = -1 ;
        ct1->mapped[i] = -1 ;
      }
    }
  }
#endif

  /* recalculate new deltas */
  dxTotal = 0 ;
  dyTotal = 0 ;
  for (i = 0 ; i < ct2->nobjects ; i++) {
    if (ct2->mapped[i] >= 0) {
      cobj2 = &ct2->objects[i] ;
      cobj1 = &ct1->objects[ct2->mapped[i]] ;
      xdiff = cobj1->col_cent - cobj2->col_cent ;
      ydiff = cobj1->row_cent - cobj2->row_cent ;
      dxTotal += xdiff ;
      dyTotal += ydiff ;
    }
  }
  *px = dx = dxTotal / npoints ;
  *py = dy = dyTotal / npoints ;

  /* now calculate total error based on new deltas */
  total_dist = 0 ;
  npoints = 0 ;
  for (i = 0 ; i < ct2->nobjects ; i++) {
    if (ct2->mapped[i] >= 0) {
      cobj2 = &ct2->objects[i] ;
      cobj1 = &ct1->objects[ct2->mapped[i]] ;
      xdiff = cobj1->col_cent - (cobj2->col_cent + dx) ;
      ydiff = cobj1->row_cent - (cobj2->row_cent + dy) ;
      dist = nint(hypot((double)xdiff, (double)ydiff)) ;
      total_dist += dist ;
      npoints++ ;
    }
  }


  *px = dx ;
  *py = dy ;
  if (npoints < MIN_POINTS)
    return(-1) ;

  /*
    now do final check.  Barring occlusion, there should be no points in
    image 1 that are mapped to points near the fovea that are not found.
    Only allow some small # of these.  If too many are found, assume the
    match is no good.
  */
  nbad = 0 ;
  for (i = 0 ; i < ct1->nobjects ; i++) {
    if (ct1->mapped[i] < 0) {
      cobj1 = &ct1->objects[i] ;

      /* find where object should have been mapped */
      cx2 = cobj1->col_cent - dx ;
      cy2 = cobj1->row_cent - dy ;

      if ((cx2 >= 0) && (cx2 < lmi->ncols) &&
          (cy2 >= 0) && (cy2 < lmi->nrows)) {
        ring = TV_TO_RING(lmi, cx2, cy2) ;
        if (ring < lmi->nrings/2)   /* not close to the periphery */
          nbad++ ;
      }
    }
  }

  if (nbad >= (npoints/2))  /* too many points that should have been found */
    return(-1) ;

  avg_dist = total_dist / npoints ;

  return(avg_dist) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              Find the object in ctable with centroid nearest to (x, y).
              If use_map is true, only consider objects with an unused
              mapped field.
----------------------------------------------------------------------*/
COBJ *
CobjFindMatch(LOGMAP_INFO *lmi, COBJ_TABLE *ctable, int x, int y,
              int use_map, int *pdist, int min_area, int max_area,
              int max_dist) {
  COBJ *cobj, *cnearest ;
  int  objno, dist, xdiff, ydiff, min_dist ;

#if 0
  /* if the point is out of the image, return failure */
  if ((x < 0) || (y < 0))
    return(NULL) ;
  if ((x >= lmi->ncols) || (y >= lmi->nrows))
    return(NULL) ;
#endif

  cnearest = NULL ;
  min_dist = max_dist ;


  for (objno = 0 ; objno < ctable->nobjects ; objno++) {
    if (use_map && (ctable->mapped[objno] >= 0))
      continue ;

    cobj = &ctable->objects[objno] ;

    xdiff = cobj->col_cent - x ;
    ydiff = cobj->row_cent - y ;
    dist = nint(hypot((double)xdiff, (double)ydiff)) ;

    /* the area doesn't seem to be much help in matching.  Objects that
       are partially out of the image will not be used if area matching
       is on, and they seem to be helpful.  Turning area matching on
       will probably be somewhat better about detecting spurious matches,
       but will throw out some good ones.
       */
#define USE_AREA_IN_MATCH 0


    /* check distance and area match */
    if ((dist < min_dist)
#if USE_AREA_IN_MATCH
        && (cobj->tv_area >= min_area) && (cobj->tv_area <= max_area)
#endif
       ) {
      min_dist = dist ;
      cnearest = cobj ;
    }
  }

  if (pdist) *pdist = min_dist ;
  return(cnearest) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
checkRow(RLE *rle) {
  RLE *first, *next ;
  int firstno, nextno ;

  if (!rle)
    return(1) ;

  firstno = 0 ;

  /* check for duplicates or bad rles */
  for (first = rle ; first ; first = first->next, firstno++) {
    next = first->next ;
    nextno = firstno + 1 ;
    while (next) {
      if (next == first) {

        fprintf(stderr, "first %d and next %d are the same!!\n",
                firstno, nextno) ;
        badRle(first) ;
      }
      next = next->next ;
      nextno++ ;
    }
  }
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
badRle(RLE *rle) {
  fprintf(stderr, "bad RLE %lx\n", (unsigned long)rle) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              Merge ct1 into ct2 and return ctDst.
----------------------------------------------------------------------*/
COBJ_TABLE *
ConnTableMerge(COBJ_TABLE *ctSrc, COBJ_TABLE *ctDst, int dx, int dy) {
  int    cSrcIndex, srcRadius, dstRadius, dist, copy, dstRing,
  srcRing, ringFovea ;
  COBJ   *cSrcObj, *cDstObj ;

  if (Gdiag & DIAG_MERGE) {
    printf("ConnTableMerge(%d, %d)\n", dx, dy) ;
    printf("src table:\n") ;
    dumpTable(ctSrc) ;
    printf("dst table:\n") ;
    dumpTable(ctDst) ;
  }

  ringFovea = ctSrc->lmi->ring_fovea ;

  /* now merge each dst object that is not accounted for into src */
  for (cSrcIndex = 0 ; cSrcIndex < ctSrc->nobjects ; cSrcIndex++) {

    cSrcObj = &ctSrc->objects[cSrcIndex] ;
    srcRadius = nint(INV_SQRTPI * sqrt((double)cSrcObj->tv_area)) ;
    cDstObj =
      CobjFindMatch(ctSrc->lmi, ctDst, cSrcObj->col_cent + dx,
                    cSrcObj->row_cent + dy, 0, &dist, 0, MAX_AREA*10,
                    srcRadius * 20);

    copy = 1 ;
    if (cDstObj)   /* possibly the same object, check overlap */
    {
      dstRadius = nint(INV_SQRTPI * sqrt((double)cDstObj->tv_area)) ;
      if (Gdiag & DIAG_MERGE)
        printf("%d: src <%d, %d> r %d --> dst <%d, %d> r %d: dist %d",
               cSrcIndex, cSrcObj->col_cent, cSrcObj->row_cent, srcRadius,
               cDstObj->col_cent, cDstObj->row_cent, dstRadius,
               dist) ;

      /* check for overlap with a fudge factor */
      if (dist > nint(RADIUS_FUDGE * (double)(srcRadius + dstRadius)))
        cDstObj = newObject(ctDst, 0, 0) ;  /* allocate new object */
      else    /* decide which object is more accurate using foveal distance */
      {
        dstRing = cDstObj->log_col_cent ;
        srcRing = cSrcObj->log_col_cent ;

        /* same object, check to see which is better */
        if (abs(dstRing-ringFovea) <= abs(srcRing-ringFovea))
          copy = 0 ;
        else     /* source obj is closer to fovea, use it */
          cobjClear(cDstObj) ;
      }
    }
    else
      cDstObj = newObject(ctDst, 0, 0) ;

    if (copy)   /* no shadow in src space, add new object to table */
    {

      CobjCopy(cSrcObj, cDstObj, dx, dy, 0) ;

      if (Gdiag & DIAG_MERGE)
        printf(" adding") ;
    }
    if (Gdiag & DIAG_MERGE)
      printf("\n") ;
  }

  return(ctDst) ;
}
void
ConnTableCompact(COBJ_TABLE *ct) {
  int   o1, o2, merged, radius1, radius2, dist, xdist, ydist, ring1, ring2,
  ringFovea, nobjects, freeno, objno ;
  COBJ  *cobj1, *cobj2, *freeObj, *cobj = NULL ;

  ringFovea = ct->lmi->ring_fovea ;

  do {
    merged = 0 ;
    for (o1 = 0 ; o1 < ct->nobjects ; o1++) {
      cobj1 = &ct->objects[o1] ;
      if (!cobj1->used)
        continue ;

      /* assume edges of object have been thresholded out (1.1) */
      radius1 = nint(INV_SQRTPI * sqrt((double)cobj1->tv_area));
      for (o2 = 0 ; o2 < ct->nobjects ; o2++) {
        cobj2 = &ct->objects[o2] ;
        if ((o1 == o2) || (!cobj2->used))
          continue ;

        radius2 = nint(INV_SQRTPI * sqrt((double)cobj2->tv_area));

        xdist = cobj2->col_cent - cobj1->col_cent ;
        ydist = cobj2->row_cent - cobj1->row_cent ;
        dist = (int)hypot((double)xdist, (double)ydist) ;

        if (dist < nint(RADIUS_FUDGE * (double)(radius1 + radius2))) {
          merged++ ;
          ring1 = cobj1->log_col_cent ;
          ring2 = cobj2->log_col_cent ;

          if (abs(ring1 - ringFovea) > abs(ring2 - ringFovea))
            cobjFree(cobj1) ;
          else
            cobjFree(cobj2) ;
        }
      }
    }
  } while (merged) ;


  /*
    compact the table.  Note that at this point, nobjects represents the
    total # of original objects found, while ctable->nobjects represents
    the # of objects after merging.
  */
  for (nobjects = freeno = 0 ; freeno < ct->nobjects ; freeno++) {
    freeObj = ct->objects + freeno ;
    if (freeObj->used) {
      nobjects++ ;
      continue ;
    }

    /*
      found a free object, now search forward for one in use to copy into
      the free one.
    */
    for (objno = freeno + 1 ; objno < ct->nobjects ; objno++) {
      cobj = ct->objects + objno ;
      if (cobj->used) break ;
    }

    /* copy used object into free one */
    if ((objno < ct->nobjects) && cobj->used) {
      isFree(freeObj) ;
      CobjCopy(cobj, freeObj, 0, 0, 1) ;
      cobjFree(cobj) ;
      isFree(cobj) ;
      freeObj->used = 1 ;
      nobjects++ ;
    }
  }

  ct->nobjects = nobjects ;
}

int
CobjOverlap(COBJ *cobj1, COBJ_TABLE *ct, int dx, int dy) {
  int    objno, xdiff, ydiff, dist, x, y, radius1, radius2 ;
  COBJ   *cobj2 ;

  x = cobj1->col_cent + dx ;
  y = cobj1->row_cent + dy ;
  radius1 = nint((double)sqrt(cobj1->tv_area)) ;
  for (objno = 0 ; objno < ct->nobjects ; objno++) {
    cobj2 = &ct->objects[objno] ;
    radius2 = nint((double)sqrt(cobj2->tv_area)) ;
    xdiff = cobj2->col_cent - x ;
    ydiff = cobj2->row_cent - y ;
    dist = nint(hypot((double)xdiff, (double)ydiff)) ;
    if (dist < (radius1 + radius2))
      return(1) ;
  }
  return(0) ;
}

int
ClipIntersect(CLIP *c1, CLIP *c2) {
  int x1start, x1end, y1start, y1end ;
  int x2start, x2end, y2start, y2end ;

  x1start = c1->x ;
  y1start = c1->y ;
  x1end = x1start + c1->dx - 1 ;
  y1end = y1start + c1->dy - 1 ;

  x2start = c2->x ;
  y2start = c2->y ;
  x2end = x2start + c2->dx - 1 ;
  y2end = y2start + c2->dy - 1 ;

  if (((x2start < x1start) || (x2start > x1end)) &&
      ((x1start < x2start) || (x1start > x2end)))
    return(0) ;

  if (((y2start < y1start) || (y2start > y1end)) &&
      ((y1start < y2start) || (y1start > y2end)))
    return(0) ;

  return(1) ;
}

void
CobjCopy(COBJ *cSrcObj, COBJ *cDstObj, int dx, int dy, int copy_rles) {
  int  y, srow, colOffset ;
  RLE  *rle ;

  cDstObj->tv_area = cSrcObj->tv_area ;
  cDstObj->log_area = cSrcObj->log_area ;
  cDstObj->col_cent = cSrcObj->col_cent + dx ;
  cDstObj->row_cent = cSrcObj->row_cent + dy ;
  cDstObj->log_row_cent = cSrcObj->log_row_cent ;
  cDstObj->log_col_cent = cSrcObj->log_col_cent ;
  cDstObj->starty = cSrcObj->starty + dy ;

  ClipCopy(&cSrcObj->clip, &cDstObj->clip) ;
  cDstObj->clip.x += dx ;
  cDstObj->clip.y += dy ;

  if (copy_rles) {
    for (srow = 0 ; srow < cSrcObj->nrows ; srow++) {
      checkRow(cSrcObj->rows[srow]) ;
      rle = cSrcObj->rows[srow] ;
      y = srow + cSrcObj->starty ;
      if (rle)                  /* this row has something in it */
      {
        /* add each pixel in the source row to the destination object */
        while (rle) {
          for (colOffset = 0 ; colOffset < rle->len ; colOffset++)
            CobjAdd(cDstObj, rle->start + colOffset + dx, y + dy) ;

          rle = rle->next ;
        }
      }
      checkRow(cSrcObj->rows[srow]) ;
    }
  }
}
void
checkObj(COBJ *cobj) {
  int row ;

  for (row = 0 ; row < cobj->nrows ; row++)
    checkRow(cobj->rows[row]) ;
}
void
CobjTranslate(COBJ *cobj, int dx, int dy) {}

static int
cobjFree(COBJ *cobj) {
  int  count, rowno, total ;
  RLE  *rle, *next ;

  cobj->used = 0 ;
  total = 0 ;
  for (rowno = 0 ; rowno < cobj->nrows ; rowno++) {
    rle = cobj->rows[rowno] ;
    count = 0 ;
    while (rle) {
      next = rle->next ;
      free(rle) ;
      RLES-- ;
      count++ ;
      rle = next ;
    }
    cobj->rows[rowno] = NULL ;
    total += count ;
  }
  isFree(cobj) ;
  return(total) ;
}

static void
cobjClear(COBJ *cobj) {
  int   rowno ;

  cobj->parent = cobj->child = NULL ;

  cobj->log_area = 0 ;
  cobj->tv_area = 0 ;
  cobj->col_cent = 0 ;
  cobj->row_cent = 0 ;
  cobj->log_row_cent = 0 ;
  cobj->log_col_cent = 0 ;

  for (rowno = 0 ; rowno < cobj->nrows ; rowno++)
    cobj->rows[rowno] = NULL ;
}
void
ConnTableFreeRles(COBJ_TABLE *ctable) {
  int   i, rowno ;
  RLE  *rle, *next ;
  COBJ *cobj ;

  for (i = 0 ; i < ctable->max_objects ; i++) {
    cobj = &ctable->objects[i] ;
    for (rowno = 0 ; rowno < cobj->nrows ; rowno++) {
      rle = cobj->rows[rowno] ;
      checkRow(rle) ;
      while (rle) {
        next = rle->next ;
        free(rle) ;
        RLES-- ;
        rle = next ;
      }
      cobj->rows[rowno] = NULL ;
    }
  }
}

void
ConnTableCountRles(COBJ_TABLE *ctable) {
  int   i, rowno, count ;
  RLE  *rle, *next ;
  COBJ *cobj ;

  return ;

  for (i = 0 ; i < ctable->max_objects ; i++) {
    cobj = &ctable->objects[i] ;
    if (i < ctable->nobjects)
      fprintf(stderr, "object %d:\n", i) ;
    for (rowno = 0 ; rowno < cobj->nrows ; rowno++) {
      rle = cobj->rows[rowno] ;
      count = 0 ;
      while (rle) {
        count++ ;
        next = rle->next ;
        rle = next ;
      }
      if (count) {
        fprintf(stderr, "row %d: %d\n", rowno, count) ;
        if (i >= ctable->nobjects) {
          fprintf(stderr, "%d rles in row %d of object %d > nobj %d!!\n",
                  count, rowno, i, ctable->nobjects) ;
        }
      }
    }
  }
}
