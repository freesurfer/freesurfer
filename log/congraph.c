/**
 * @file  congraph.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:11 $
 *    $Revision: 1.7 $
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
  @(#)congraph.c  1.3
  8/2/95
*/
/*----------------------------------------------------------------------
           File Name:
             congraph.c

           Author:
             Bruce Fischl with algorithms stolen from Rich Wallace.

           Description:
             Build the connectivity graph.

           Conventions:

----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           INCLUDE FILES
----------------------------------------------------------------------*/

#define USE_XWIN  0

#include <stdio.h>
#include <math.h>
#include <hipl_format.h>

#include "image.h"
#include "h_logz.h"
#include "congraph.h"
#include "diag.h"
#include "const.h"
#include "proto.h"
#include "error.h"

/*----------------------------------------------------------------------
                            CONSTANTS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           STATIC DATA
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           FUNCTION PROTOTYPES
----------------------------------------------------------------------*/

static NEIGHBOR *newNeighbor(LOGMAP_INFO *lmi) ;
static NEIGHBOR *appendNeighbor(NEIGHBOR *neighbor, NEIGHBOR *next) ;
static int   neighborP(LOGMAP_INFO *lmi, int ring, int spoke, int neigh_ring,
                       int neigh_spoke) ;
static void  unionNeighbors(LOGMAP_INFO *lmi, int ring, int spoke,
                            int neigh_ring, int neigh_spoke);
static void conGraphShow(LOGMAP_INFO *lmi) ;
static void  orderNeighbors(LOGMAP_INFO *lmi, int ring, int spoke) ;
static void initEightConnected(LOGMAP_INFO *lmi) ;
static int findConnectedRing(LOGMAP_INFO *lmi, int ring, int spoke, int dir) ;

#if 0
static void printNeighbors(LOGMAP_INFO *lmi, int ring, int spoke) ;
static int findConnectedSpoke(LOGMAP_INFO *lmi, int ring, int spoke, int dir) ;
#endif

/*----------------------------------------------------------------------
                           FUNCTIONS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:

           Description:
             Initialize the connectivity graph.  Any logmap pixel that
             contains Cartesian pixels is a neighbor of all logmap pixels
             which contain Cartesian pixels that are neighbors in
             Cartesian space.  This will give connectivity across the
             vertical meridian.  In addition, logmap pixels are also
             neighbors if the neighbor in log space.  This is needed
             to account for aliasing.
----------------------------------------------------------------------*/
void
ConGraphInit(LOGMAP_INFO *lmi) {
  int i, j;
  int spoke, ring, neigh_ring, neigh_spoke, ncols, nrows, nrings, nspokes ;

  nrings = lmi->nrings ;
  nspokes = lmi->nspokes ;
  ncols = lmi->ncols ;
  nrows = lmi->nrows ;

  /* allocate more neighbors than we will ever need */
  lmi->max_neighbors = nrings * nspokes * 8 ;
  lmi->neighbors =
    (NEIGHBOR *)calloc(lmi->max_neighbors, sizeof(NEIGHBOR)) ;
  if (!lmi->neighbors) {
    fprintf(stderr, "ConGraphInit: could not allocate %d neighbors\n",
            lmi->max_neighbors) ;
    exit(1) ;
  }

  /*
    all log pixels which contains Cartesian pixels which are neighbors
    are also neighbors in log space.
  */
  for_each_tv_pixel(lmi, i, j) {
    ring = TV_TO_RING(lmi, i, j) ;
    spoke = TV_TO_SPOKE(lmi, i, j) ;
    if (DEFINED(ring, spoke) && LOG_PIX_AREA(lmi, ring, spoke) > 0) {
      if (i > 0) {
        neigh_ring = TV_TO_RING(lmi, i-1, j) ;
        neigh_spoke = TV_TO_SPOKE(lmi, i-1, j) ;
        unionNeighbors(lmi, ring, spoke, neigh_ring, neigh_spoke) ;
      }
      if (i < ncols-1) {
        neigh_ring = TV_TO_RING(lmi, i+1, j) ;
        neigh_spoke = TV_TO_SPOKE(lmi, i+1, j) ;
        unionNeighbors(lmi, ring, spoke, neigh_ring, neigh_spoke) ;
      }
      if (j > 0) {
        neigh_ring = TV_TO_RING(lmi, i, j-1) ;
        neigh_spoke = TV_TO_SPOKE(lmi, i, j-1) ;
        unionNeighbors(lmi, ring, spoke, neigh_ring, neigh_spoke) ;
      }
      if (j < nrows-1) {
        neigh_ring = TV_TO_RING(lmi, i, j+1) ;
        neigh_spoke = TV_TO_SPOKE(lmi, i, j+1) ;
        unionNeighbors(lmi, ring, spoke, neigh_ring, neigh_spoke) ;
      }
      if (i > 0 && j > 0) {
        neigh_ring = TV_TO_RING(lmi, i-1, j-1) ;
        neigh_spoke = TV_TO_SPOKE(lmi, i-1, j-1) ;
        unionNeighbors(lmi, ring, spoke, neigh_ring, neigh_spoke) ;
      }
      if (i < ncols-1 && j > 0) {
        neigh_ring = TV_TO_RING(lmi, i+1, j-1) ;
        neigh_spoke = TV_TO_SPOKE(lmi, i+1, j-1) ;
        unionNeighbors(lmi, ring, spoke, neigh_ring, neigh_spoke) ;
      }
      if (i > 0 && j < nrows-1) {
        neigh_ring = TV_TO_RING(lmi, i-1, j+1) ;
        neigh_spoke = TV_TO_SPOKE(lmi, i-1, j+1) ;
        unionNeighbors(lmi, ring, spoke, neigh_ring, neigh_spoke) ;
      }
      if (i < ncols-1 && j < nrows-1) {
        neigh_ring = TV_TO_RING(lmi, i+1, j+1) ;
        neigh_spoke = TV_TO_SPOKE(lmi, i+1, j+1) ;
        unionNeighbors(lmi, ring, spoke, neigh_ring, neigh_spoke) ;
      }
    }
  }

  /*
    due to aliasing, some log pixels which are neighbors will not contain
    any Cartesian pixels which are neighbors.  Find them to complete the
    connectivity graph.
  */
  for_each_log_pixel(lmi, ring, spoke) {
    if (spoke > 0 &&  LOG_PIX_AREA(lmi, ring, spoke) > 0)
      unionNeighbors(lmi, ring, spoke, ring, spoke-1);
    if (spoke < nspokes-1 &&  LOG_PIX_AREA(lmi, ring, spoke+1) > 0)
      unionNeighbors(lmi, ring, spoke, ring, spoke+1);
    if (ring > 0  &&  LOG_PIX_AREA(lmi, ring-1, spoke) > 0)
      unionNeighbors(lmi, ring, spoke, ring-1, spoke);
    if (ring < nrings-1  &&  LOG_PIX_AREA(lmi, ring+1, spoke) > 0)
      unionNeighbors(lmi, ring, spoke, ring+1, spoke);

    if (spoke > 0 && ring > 0 && LOG_PIX_AREA(lmi, ring-1, spoke-1) > 0)
      unionNeighbors(lmi, ring, spoke, ring-1, spoke-1);
    if (spoke < nspokes-1 && ring > 0 &&
        LOG_PIX_AREA(lmi, ring-1, spoke+1) > 0)
      unionNeighbors(lmi, ring, spoke, ring-1, spoke+1);
    if (spoke> 0 && ring < nrings-1 && LOG_PIX_AREA(lmi,ring+1,spoke-1) >0)
      unionNeighbors(lmi, ring, spoke, ring+1, spoke-1);
    if (spoke < nspokes-1 && ring < nrings-1 &&
        LOG_PIX_AREA(lmi, ring+1, spoke+1) > 0)
      unionNeighbors(lmi, ring, spoke, ring+1, spoke+1);
  }

  /* now order the neighbors of a pixel by spoke then ring */
  for_each_log_pixel(lmi, ring, spoke)
  orderNeighbors(lmi, ring, spoke) ;

  initEightConnected(lmi) ;
  if (Gdiag & DIAG_CONGRAPH) conGraphShow(lmi) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
             create a new NEIGHBOR structure and return it.
----------------------------------------------------------------------*/
static NEIGHBOR *
newNeighbor(LOGMAP_INFO *lmi) {
  if (lmi->n_neighbors >= lmi->max_neighbors)
    ErrorReturn(lmi->neighbors + lmi->n_neighbors,
                (ERROR_NO_MEMORY,
                 "newNeighbor: too many total neighbors (%d)\n",
                 lmi->n_neighbors)) ;


  return(lmi->neighbors + lmi->n_neighbors++) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             Append the neighbor list 'next' to 'neighbor'
----------------------------------------------------------------------*/
static NEIGHBOR *
appendNeighbor(NEIGHBOR *neighbor, NEIGHBOR *next) {
  if (neighbor == NULL) return(next) ;
  else {
    neighbor->next = next ;
    return(neighbor) ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             determine whether the logmap pixel specified by neigh_ring
             and neigh_spoke is a neighbor of the one at ring,spoke.
----------------------------------------------------------------------*/
static int
neighborP(LOGMAP_INFO *lmi,int ring,int spoke, int neigh_ring, int neigh_spoke) {
  NEIGHBOR *neighbor;

  for_each_neighbor(lmi, neighbor, ring, spoke)
  if (LOG_PIX(lmi, neigh_ring, neigh_spoke) == neighbor->logpix)
    return(1);
  return(0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             If neigh_ring,neigh_spoke is not already a neighbor, add
             it to ring,spoke's neighbor list.
----------------------------------------------------------------------*/
static void
unionNeighbors(LOGMAP_INFO *lmi, int ring, int spoke, int neigh_ring,
               int neigh_spoke) {
  NEIGHBOR *neighbor;

  if (!DEFINED(neigh_ring, neigh_spoke)) return ;

  if ((neigh_spoke != spoke || neigh_ring != ring) &&
      /*
      ** speedup trick: a pixel can't be its own neighbor
      */
      LOG_PIX_AREA(lmi, neigh_ring, neigh_spoke) > 0 &&
      neighborP(lmi, ring, spoke, neigh_ring, neigh_spoke) == 0) {
    neighbor = newNeighbor(lmi);
    neighbor->logpix = LOG_PIX(lmi, neigh_ring, neigh_spoke) ;
    LOG_PIX_NEIGHBORS(lmi, ring, spoke) =
      appendNeighbor(neighbor, LOG_PIX_NEIGHBORS(lmi, ring, spoke)) ;
    LOG_PIX_N_NEIGHBORS(lmi, ring, spoke)++ ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             display the connectivity graph in Cartesian and in Log space.
----------------------------------------------------------------------*/
static void
conGraphShow(LOGMAP_INFO *lmi) {
#if USE_XWIN
  int       spoke, ring, area, radius, x0, y0, x1, y1 ;
  NEIGHBOR  *neighbor ;

  DebugNewWindow(lmi->ncols, lmi->nrows,
                 "log pixel centroids, area, and connectivity", 3) ;

  for_each_log_pixel(lmi, ring, spoke) {
    area = LOG_PIX_AREA(lmi, ring, spoke) ;
    radius = nint(sqrt((double)area / PI)) / 2 ;

    x0 = LOG_PIX_COL_CENT(lmi, ring, spoke) ;
    y0 = LOG_PIX_ROW_CENT(lmi, ring, spoke) ;
    if (radius > 1)
      DebugDrawCircle(x0, y0, radius, xWHITE) ;

    for_each_neighbor(lmi, neighbor, ring, spoke) {
      x1 = neighbor->logpix->col_cent ;
      y1 = neighbor->logpix->row_cent ;

      DebugBufferLine(x0, y0, x1, y1, xWHITE) ;
    }
  }

  DebugFlushLines(xWHITE) ;
  DebugKeyboard() ;
  DebugEndWindow() ;
#endif
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             order the neighbors of a pixel, with the lowest spoke and
             lowest ring.
----------------------------------------------------------------------*/
static void
orderNeighbors(LOGMAP_INFO *lmi, int ring, int spoke) {
  NEIGHBOR *nOne, *nTwo, *nPrev ;
  int      swapped ;

  /* do a bubble sort */
  do {
    swapped = 0 ;
    /*     printNeighbors(lmi, ring, spoke) ;*/
    nOne = nTwo = nPrev = NULL ;
    for_each_neighbor(lmi, nOne, ring, spoke) {
      nTwo = nOne->next ;
      if (nTwo) {
        if ((nTwo->logpix->spoke < nOne->logpix->spoke) ||
            ((nTwo->logpix->spoke == nOne->logpix->spoke) &&
             (nTwo->logpix->ring < nOne->logpix->ring))) {
          /* swap nOne and nTwo */
          swapped = 1 ;
          nOne->next = nTwo->next ;
          nTwo->next = nOne ;

          /* update item previous to nOne */
          if (!nPrev)   /* nOne is first in list */
            LOG_PIX_NEIGHBORS(lmi, ring, spoke) = nTwo ;
          else
            nPrev->next = nTwo ;
          nOne = nTwo ;   /* they have swapped positions */
        }
      }
      nPrev = nOne ;
    }
  } while (swapped) ;

  /*  printNeighbors(lmi, ring, spoke) ;*/
}

#if 0
static void
printNeighbors(LOGMAP_INFO *lmi, int ring, int spoke) {
#if 1
  int    k, n_ring, n_spoke ;
  LOGPIX *npix ;

  printf("neighbors(%d, %d)\n", ring, spoke) ;
  for (k = 0 ; k < NBD_SIZE ; k++) {
    n_ring = LOG_PIX_NBD(lmi, ring, spoke, k)->ring ;
    n_spoke = LOG_PIX_NBD(lmi, ring, spoke, k)->spoke ;
    npix = LOG_PIX(lmi, n_ring, n_spoke) ;
    printf("(%d, %d)\n", n_ring, n_spoke) ;
  }
#else
  NEIGHBOR *neighbor ;

  printf("neighbors(%d, %d)\n", ring, spoke) ;
  for_each_neighbor(lmi, neighbor, ring, spoke) {
    printf("0x%lx (%d, %d)\n", (unsigned long)neighbor, neighbor->logpix->ring,
           neighbor->logpix->spoke) ;
  }
#endif
  fflush(stdout) ;
}
#endif

extern int debug(void) ;

static void
initEightConnected(LOGMAP_INFO *lmi) {
  int      k, n_ring ;
  int      spoke, ring, ncols, nrows, nrings, nspokes ;
  LOGPIX  *pix ;

  nrings = lmi->nrings ;
  nspokes = lmi->nspokes ;
  ncols = lmi->ncols ;
  nrows = lmi->nrows ;

  for_each_log_pixel(lmi, ring, spoke) {
    /* use Neumann boundary conditions */
    for (k = 0 ; k < NBD_SIZE ; k++)
      LOG_PIX_NBD(lmi, ring, spoke, k) = LOG_PIX(lmi, ring, spoke) ;

    if (LOG_PIX_AREA(lmi, ring, spoke) <= 0)
      continue ;

    pix = LOG_PIX(lmi, ring, spoke) ;

    if (ring >= nrings/2)  /* in right half-plane */
    {
      /* find north-eastern neighbor */
      n_ring = findConnectedRing(lmi, ring+1, spoke+1, 1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_NE) = LOG_PIX(lmi, n_ring, spoke+1) ;

      /* find north-western neighbor */
      n_ring = findConnectedRing(lmi, ring-1, spoke+1, -1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_NW) = LOG_PIX(lmi, n_ring, spoke+1) ;

      /* find eastern neighbor */
      n_ring = findConnectedRing(lmi, ring+1, spoke, 1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_E) = LOG_PIX(lmi, n_ring, spoke) ;

      /* find northern neighbor */
      n_ring = findConnectedRing(lmi, ring, spoke+1, -1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_N) = LOG_PIX(lmi, n_ring, spoke+1) ;

      /* find southern neighbor */
      n_ring = findConnectedRing(lmi, ring, spoke-1, -1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_S) = LOG_PIX(lmi, n_ring, spoke-1) ;

      /* find south-eastern neighbor */
      n_ring = findConnectedRing(lmi, ring+1, spoke-1, 1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_SE) = LOG_PIX(lmi, n_ring, spoke-1) ;

      /* find western neighbor */
      n_ring = findConnectedRing(lmi, ring-1, spoke, -1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_W) = LOG_PIX(lmi, n_ring, spoke) ;

      /* find southwestern neighbor */
      n_ring = findConnectedRing(lmi, ring-1, spoke-1, -1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_SW) = LOG_PIX(lmi, n_ring, spoke-1) ;
    } else                                          /* in left half-plane */
    {
      /* find north-eastern neighbor */
      n_ring = findConnectedRing(lmi, ring+1, spoke+1, 1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_NE) = LOG_PIX(lmi, n_ring, spoke+1) ;

      /* find north-western neighbor */
      n_ring = findConnectedRing(lmi, ring-1, spoke+1, -1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_NW) = LOG_PIX(lmi, n_ring, spoke+1) ;

      /* find eastern neighbor */
      n_ring = findConnectedRing(lmi, ring+1, spoke, 1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_E) = LOG_PIX(lmi, n_ring, spoke) ;

      /* find northern neighbor */
      n_ring = findConnectedRing(lmi, ring, spoke+1, 1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_N) = LOG_PIX(lmi, n_ring, spoke+1) ;

      /* find southern neighbor */
      n_ring = findConnectedRing(lmi, ring, spoke-1, 1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_S) = LOG_PIX(lmi, n_ring, spoke-1) ;

      /* find south-eastern neighbor */
      n_ring = findConnectedRing(lmi, ring+1, spoke-1, 1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_SE) = LOG_PIX(lmi, n_ring, spoke-1) ;

      /* find western neighbor */
      n_ring = findConnectedRing(lmi, ring-1, spoke, -1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_W) = LOG_PIX(lmi, n_ring, spoke) ;

      /* find southwestern neighbor */
      n_ring = findConnectedRing(lmi, ring-1, spoke-1, -1) ;
      if (n_ring >= 0)
        LOG_PIX_NBD(lmi, ring, spoke, N_SW) = LOG_PIX(lmi, n_ring, spoke-1) ;
    }

#if 0
    if (ring == 43 && spoke == 40)
      printNeighbors(lmi, ring, spoke) ;
#endif
#if 0
    {
      FILE *fp ;

      fp = fopen("nbd.dat", "w") ;
      for (k = 0 ; k < NBD_SIZE ; k++) {
        npix = LOG_PIX_NBD(lmi, ring, spoke, k) ;
        if ((npix->ring < 0) || (npix->row_cent <= 0) || (npix->col_cent <= 0)) {
          fprintf(stderr, "(%d,%d,%d) - (%d, %d, %d, %d), area %d\n",
                  ring, spoke, k, npix->ring, npix->spoke, npix->col_cent,
                  npix->row_cent, npix->area) ;
        } else {
          fprintf(fp, "%d  %d  %d  %d\n", ring, spoke,
                  LOG_PIX_COL_CENT(lmi,ring,spoke),
                  LOG_PIX_ROW_CENT(lmi, ring, spoke)) ;
          fprintf(fp, "%d  %d  %d  %d\n\n", npix->ring, npix->spoke,
                  npix->col_cent, npix->row_cent) ;
        }
      }
      fclose(fp) ;
    }
#endif
  }
}

#if 0
static int
findConnectedSpoke(LOGMAP_INFO *lmi, int ring, int spoke, int dir) {
  int nspokes ;

  nspokes = lmi->nspokes ;

  /* find closest spoke with positive area */

  if ((spoke >= lmi->nspokes) || (spoke < 0) || (ring >= lmi->nrings) ||
      (ring < 0))
    return(-1) ;

  do {
    if (LOG_PIX_AREA(lmi, ring, spoke) > 0)
      return(spoke) ;
    spoke += dir ;          /* search left or right */

  }   while ((spoke < nspokes) && (spoke >= 0)) ;

  return(-1) ;
}
#endif

static int
findConnectedRing(LOGMAP_INFO *lmi, int ring, int spoke, int dir) {
  int nrings ;

  /* find closest ring with positive area */

  if ((spoke >= lmi->nspokes) || (spoke < 0) || (ring >= lmi->nrings) ||
      (ring < 0))
    return(-1) ;

  nrings = lmi->nrings ;
  do {
    if (LOG_PIX_AREA(lmi, ring, spoke) > 10000) {
      fprintf(stderr, "area(%d, %d) = %d!\n", ring, spoke,
              LOG_PIX_AREA(lmi, ring, spoke)) ;
      return(-1) ;
    }

    if (LOG_PIX_AREA(lmi, ring, spoke) > 0)
      return(ring) ;
    ring += dir ;          /* search left or right */

  }   while ((ring < nrings) && (ring >= 0)) ;
  return(-1) ;
}
