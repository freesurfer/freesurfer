/**
 * @file  connect.h
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
  @(#)connect.h 1.1
  3/31/94
*/
/*------------------------------------------------------------------------
      File Name:  connect.h

         Author:  Bruce Fischl

        Created:  Jan. 1994

    Description:

------------------------------------------------------------------------*/
#ifndef CONNECT_H
#define CONNECT_H

#include <hipl_format.h>

#include "h_logz.h"
#include "image.h"
#include "macros.h"
#include "struct.h"


typedef struct _s_rle
{
  struct _s_rle *next ;   /* next run of continuous pixels */
  int           start ;   /* starting column of run */
  int           len ;     /* length of run */
}
RLE ;

/*
  images are typically in log coordinates, so clips, dx and dy, and row
  #s will all be in log space.  The centroids are in Cartesian space.
*/
/*
  starty defines the y coordinate of rows[0].
*/
typedef struct _s_cobj
{
  int     id ;               /* the id of this object (also its index) */
  int     used ;             /* 1 if an object is actually used */
  CLIP    clip ;             /* bounding box of object (log space) */
  RLE     **rows ;           /* one link-list of runs for each row in object */
  int     nrows ;
  struct  _s_cobj *child ;   /* object that will be merged into this one */
  struct  _s_cobj *parent ;  /* object that this will be merged into */
  int     col_cent ;         /* column centroid of the object (Cartesian) */
  int     row_cent ;         /* row centroid of the object (Cartesian) */
  int     log_area ;         /* # of log space pixels in this object */
  int     tv_area ;          /* # of Cartesian pixels in this object */
  int     starty ;
  int     log_col_cent ;     /* column center in log space */
  int     log_row_cent ;     /* row center in log space */
}
CONNECTED_OBJECT, COBJ ;

typedef struct
{
  int       nobjects ;    /* # of objects actually in the table */
  int       max_objects ; /* maximum size of table */
  int       nrows ;
  COBJ      *objects ;    /* array of objects */
  int       *mapped ;     /* the object this object is mapped to in matching */
  IMAGE     *image ;      /* image that was segmented */
  int       thresh ;      /* threshold for segmentation */
  int       white ;       /* 1 is white connected, 0 if black connected */
  LOGMAP_INFO *lmi ;      /* the logmap info structure used in segmentation */
  int       xo, yo ;      /* origin of table */
}
CONNECTED_OBJECT_TABLE, COBJ_TABLE ;

COBJ_TABLE  *ConnectedComponents(LOGMAP_INFO *lmi, IMAGE *image, int white,
                                 int thresh, int xo, int yo) ;
void        ConnTableFree(COBJ_TABLE **co) ;
void        ConnTableCalculateCentroids(COBJ_TABLE *ctable) ;
int         ConnTableMatch(LOGMAP_INFO *lmi, COBJ_TABLE *ct1, COBJ_TABLE *ct2,
                           int *x, int *y, int dx_old, int dy_old,int use_old);
COBJ *      CobjFindMatch(LOGMAP_INFO *lmi, COBJ_TABLE *ctable, int x, int y,
                          int use_map, int *pdistance, int min_area,
                          int max_area, int max_dist) ;
COBJ_TABLE  *ConnTableMerge(COBJ_TABLE *ct1, COBJ_TABLE *ct2, int dx, int dy) ;
void        CobjCopy(COBJ *cSrcObj,COBJ *cDstObj,int dx,int dy,int copy_rles) ;
void        CobjTranslate(COBJ *cobj, int dx, int dy) ;
void        ConnTableFreeRles(COBJ_TABLE *ctable) ;
void        ConnTableCompact(COBJ_TABLE *ct) ;

#define CONNECT_BLACK     0
#define CONNECT_WHITE     1
#define NO_CONTINUITY     -1

#endif
