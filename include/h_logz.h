/*
  @(#)logmap.h  1.18
  8/10/95
*/
/*------------------------------------------------------------------------
      File Name:  logmap.h

         Author:  Bruce Fischl

        Created:  Jan. 1993

    Description:  

------------------------------------------------------------------------*/


#ifndef MAPINFO_H
#define MAPINFO_H

#include "macros.h"
#include "image.h"

typedef struct _neighbor_list
{
  struct _LOGPIX        *logpix ;
  struct _neighbor_list *next ;
} NEIGHBOR ;

/* 
  constants for eight neighbors and self.  Note that these must match
  the positions in a 3x3 matrix (row major order), where south and west 
  correspond to lower indices (down and to the left).
*/
#define N_SW   0
#define N_S    1
#define N_SE   2
#define N_W    3
#define N_SELF 4
#define N_E    5
#define N_NW   6
#define N_N    7
#define N_NE   8
#define NBD_SIZE 9

typedef struct _LOGPIX
{
  int   area ;           /* # of Cartesian pixels that map to this log pix */
  int   row_cent ;       /* centroid of logmap pixel in Cartesian space */
  int   col_cent ;       /* centroid of logmap pixel in Cartesian space */
  int   spoke ;          /* the spoke # of this pixel */
  int   ring ;           /* the ring # of this pixel */
  int   n_neighbors ;    /* # of neighbors of this logpix */
  NEIGHBOR *neighbors ;  /* a link-list */
  void  *user ;          /* misc. pointer for whoever */
  double weight ;         /* for gradient calculation */
  struct _LOGPIX    *nbd[9] ;    /* 8-connected neighbors and self */
} LOGPIX ;

typedef struct
{
  int     nrings ;       /* # of rings in logmap */
  int     nspokes ;      /* # of spokes in logmap */
  int     ncols ;        /* # of cols in Cartesian image */
  int     nrows ;        /* # of rows in Cartesian image */
  int     nruns ;        /* # of runs in runlength lookup tables */
  double  alpha ;        /* alpha used in complex log mapping */
  LOGPIX  *logPix ;      /* info about each logmap pixel (nrings x nspokes) */
  int     *runNoToRing ; /* maps a run # to the ring it resides in. */
  int     *runNoToSpoke; /* maps a run # to the ring it resides in. */
  char    *runNoToLen ;  /* contains a list of all run lengths. */
  UCHAR   *tvToRing ;    /* maps Cartesian (x, y) to a ring # */
  UCHAR   *tvToSpoke ;   /* maps Cartesian (x, y) to a spoke # */
  NEIGHBOR *neighbors ;  /* total # of nbrs needed for connectivity graph */
  int     max_neighbors ; /* size of neighbors allocation */
  int     n_neighbors ;  /* # of neighbors allocated so far */
  int     ring_fovea ;
  double  maxr ;          /* max radius in domain */
} LOGMAP_INFO ;

/*
  These macros give access to the mapping from Cartesian coordinates to log
  space coordinates.  In the macros r = row (y), c = col (x) 
*/
#define TV_TO_SPOKE(m, c, r)   (*((m)->tvToSpoke + ((r) * m->ncols) + c))
#define TV_TO_RING(m, c, r)    (*((m)->tvToRing + ((r) * m->ncols) + c))

/* 
  these macros deal with obtaining information about a log space pixel.
  In the macros r = ring (x), s = spoke (y)
*/
#define LOG_PIX(m, r, s)         ((m)->logPix + ((s) * m->nrings) + r)
#define LOG_PIX_AREA(m, r, s)    LOG_PIX(m, r, s)->area
#define LOG_PIX_ROW_CENT(m,r,s)  LOG_PIX(m, r, s)->row_cent
#define LOG_PIX_COL_CENT(m,r,s)  LOG_PIX(m, r, s)->col_cent
#define LOG_PIX_SPOKE(m,r,s)     LOG_PIX(m, r, s)->spoke
#define LOG_PIX_RING(m,r,s)      LOG_PIX(m, r, s)->ring
#define LOG_PIX_WEIGHT(m,r,s)    LOG_PIX(m, r, s)->weight
#define LOG_PIX_NEIGHBORS(m,r,s) LOG_PIX(m, r, s)->neighbors
#define LOG_PIX_N_NEIGHBORS(m,r,s) LOG_PIX(m, r, s)->n_neighbors
#define LOG_PIX_NBD(m,r,s,d)     LOG_PIX(m,r,s)->nbd[d]


/*
  These macros allow access to information about run lengths used in
  doing a fast construction of the logmap image.
*/
#define RUN_NO_TO_LEN(m, no)     ((m)->runNoToLen[no])
#define RUN_NO_TO_RING(m, no)    ((m)->runNoToRing[no])
#define RUN_NO_TO_SPOKE(m, no)   ((m)->runNoToSpoke[no])

LOGMAP_INFO *LogMapInit(double alpha,int cols,int rows,int nrings,int nspokes);
int    LogMapForward(LOGMAP_INFO *mi, IMAGE *inImage, IMAGE *outImage);
IMAGE  *LogMapSample(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst) ;
int    LogMapInverse(LOGMAP_INFO *mi, IMAGE *inImage, IMAGE *outImage);
int    LogMapGradient(LOGMAP_INFO *mi, IMAGE *inImage, 
                      IMAGE *gradImage, int doweight) ;
double LogMapDiffuse(LOGMAP_INFO *mi, IMAGE *inImage, IMAGE *outImage, 
                     double k, int niter, int doweight, int which) ;
double LogMapDiffusePerona(LOGMAP_INFO *lmi, IMAGE *inImage, 
                           IMAGE *outImage, double k, int niterations, 
                           int doweight) ;
double LogMapDiffuseCurvature(LOGMAP_INFO *lmi, IMAGE *inImage, 
                              IMAGE *outImage, double A, int niterations, 
                              int doweight) ;

int   LogMapCurvature(LOGMAP_INFO *lmi, IMAGE *inImage, 
                      IMAGE *gradImage, float A, int doweight) ;

void  LogMapPatchHoles(LOGMAP_INFO *lmi, IMAGE *logImage) ;

#if 0
#define for_each_neighbor(lmi, ptr, i, j, r, s) \
    for (ptr = LOG_PIX_NEIGHBORS(lmi, i, j) ; ptr ; ptr = ptr->next)  \
    if (((s = ptr->logpix->spoke)|1) && ((r = ptr->logpix->ring)|1))
#else
#define for_each_neighbor(lmi, ptr, i, j) \
    for (ptr = LOG_PIX_NEIGHBORS(lmi, i, j) ; ptr ; ptr = ptr->next)
#endif

#define for_each_log_pixel(lmi, r, s) \
    for (s = 0; s < lmi->nspokes; s++) \
      for (r = 0; r < lmi->nrings; r++) \
           if (LOG_PIX_AREA(lmi, r, s) > 0)

/* this macro will go through all log pixels, even those with no area */
#define for_all_log_pixels(lmi, r, s) \
    for (s = 0; s < lmi->nspokes; s++) \
      for (r = 0; r < lmi->nrings; r++)

#define for_each_tv_pixel(lmi, i, j) \
    for (j = 0; j < lmi->nrows ; j++)  \
      for (i = 0; i < lmi->ncols ; i++) 
      

#define UNDEFINED       255
#define DEFINED(r, s)   ((r != UNDEFINED) && (s != UNDEFINED))

#endif
