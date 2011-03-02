/**
 * @file  h_logz.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.24 $
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
}
NEIGHBOR ;

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
  float xcent ;          /* centroid of logmap pixel in Cartesian space */
  float ycent ;          /* centroid of logmap pixel in Cartesian space */
  int   spoke ;          /* the spoke # of this pixel */
  int   ring ;           /* the ring # of this pixel */
  int   n_neighbors ;    /* # of neighbors of this logpix */
  NEIGHBOR *neighbors ;  /* a link-list */
  void  *user ;          /* misc. pointer for whoever */
  double weight ;         /* for gradient calculation */
  struct _LOGPIX    *nbd[9] ;    /* 8-connected neighbors and self */
  float  rho ;           /* radial log coordinate */
  float  phi ;           /* angular log coordinate */
  int    ncpix ;         /* # of cart. pix that map to this logpix */
}
LOGPIX ;

#define MAX_PIX  10
typedef struct
{
  int npix ;      /* # of log pix that this cart. pixel contributes to */
  /* array of log pixels that this pixel contributes to */
  LOGPIX *logpix[MAX_PIX] ;
}
CARTESIAN_PIXEL, CP ;

typedef struct
{
  int    rows ;
  int    cols ;
  CP     *pix ;    /* array of LUTs for doing mapping */
}
CARTESIAN_MAP_INFO, CMI ;

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
  float   min_rho ;
  float   max_rho ;
  float   *rhos ;
  int     *start_spoke ;  /* first spoke with real data in each ring */
  int     *end_spoke ;    /* last  spoke with real data in each ring */
  CMI     cmi ;           /* for overlapping logmap construction */
}
LOGMAP_INFO ;

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
#define LOG_PIX_ROW_FCENT(m,r,s)  LOG_PIX(m, r, s)->ycent
#define LOG_PIX_COL_FCENT(m,r,s)  LOG_PIX(m, r, s)->xcent
#define LOG_PIX_SPOKE(m,r,s)     LOG_PIX(m, r, s)->spoke
#define LOG_PIX_RING(m,r,s)      LOG_PIX(m, r, s)->ring
#define LOG_PIX_WEIGHT(m,r,s)    LOG_PIX(m, r, s)->weight
#define LOG_PIX_NEIGHBORS(m,r,s) LOG_PIX(m, r, s)->neighbors
#define LOG_PIX_N_NEIGHBORS(m,r,s) LOG_PIX(m, r, s)->n_neighbors
#define LOG_PIX_NBD(m,r,s,d)     LOG_PIX(m,r,s)->nbd[d]
#define LOG_PIX_RHO(m,r,s)       LOG_PIX(m, r, s)->rho
#define LOG_PIX_PHI(m,r,s)       LOG_PIX(m, r, s)->phi


/*
  These macros allow access to information about run lengths used in
  doing a fast construction of the logmap image.
*/
#define RUN_NO_TO_LEN(m, no)     ((m)->runNoToLen[no])
#define RUN_NO_TO_RING(m, no)    ((m)->runNoToRing[no])
#define RUN_NO_TO_SPOKE(m, no)   ((m)->runNoToSpoke[no])

/* constants for LogMapInitForwardFilter */
#define LMI_FORWARD_FILTER_CIRCLE    0
#define LMI_FORWARD_FILTER_GAUSSIAN  1

void    LogMapFree(LOGMAP_INFO **plmi) ;
LOGMAP_INFO *LogMapInit(double alpha,int cols,int rows,int nrings,
                        int nspokes);
int    LogMapForward(LOGMAP_INFO *mi, IMAGE *inImage, IMAGE *outImage);
IMAGE  *LogMapSample(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst) ;
IMAGE  *LogMapInverseSample(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst) ;
IMAGE  *LogMapNonlocal(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Ismooth,
                       IMAGE *Idst) ;
int    LogMapInverse(LOGMAP_INFO *mi, IMAGE *inImage, IMAGE *outImage);
int    LogMapGradient(LOGMAP_INFO *mi, IMAGE *inImage,
                      IMAGE *gradImage, int doweight,
                      int start_ring, int end_ring) ;
int    LogMapSobel(LOGMAP_INFO *lmi, IMAGE *inImage, IMAGE *gradImage,
                   IMAGE *Ix, IMAGE *Iy, int doweight, int start_ring,
                   int end_ring) ;
IMAGE  *LogMapFilter(LOGMAP_INFO *lmi, int which, int window_size,
                     IMAGE *Isrc, IMAGE *Idst) ;
IMAGE  *LogMapMeanFilter(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst) ;
double LogMapDiffuse(LOGMAP_INFO *mi, IMAGE *inImage, IMAGE *outImage,
                     double k, int niter, int doweight, int which,
                     int time_type) ;
double LogMapDiffusePerona(LOGMAP_INFO *lmi, IMAGE *inImage,
                           IMAGE *outImage, double k, int niterations,
                           int doweight, int time_type) ;
double LogMapDiffuseCurvature(LOGMAP_INFO *lmi, IMAGE *inImage,
                              IMAGE *outImage, double A, int niterations,
                              int doweight, int time_type) ;

int   LogMapCurvature(LOGMAP_INFO *lmi, IMAGE *inImage, IMAGE *gradImage,
                      float A, int doweight, int start_ring, int end_ring) ;

void  LogMapPatchHoles(LOGMAP_INFO *lmi, IMAGE *Itv, IMAGE *Ilog) ;
IMAGE *LogMapNormalize(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst,
                       float low, float hi) ;
int   LogMapInitForwardFilter(LOGMAP_INFO *lmi, int which) ;
IMAGE *LogMapForwardFilter(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst) ;
IMAGE *LogMapInverseFilter(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst);
IMAGE *LogMapInverseBilinear(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst);
IMAGE *LogMapInverseFilterGaussian(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst);


IMAGE *LogMapMedianFilter(LOGMAP_INFO *lmi, IMAGE *Isrc, int wsize,
                          IMAGE *Ioffset, IMAGE *Idst) ;

/* offset filtering in the log plane */
IMAGE *LogMapOffsetOrientation(LOGMAP_INFO *lmi, int wsize, IMAGE *Isrc,
                               IMAGE *Iorient) ;
IMAGE *LogMapOffsetDirection(LOGMAP_INFO *lmi, IMAGE *Iorient, IMAGE *Ioffset);
IMAGE *LogMapOffsetMagnitude(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst,
                             int maxsteps) ;
IMAGE *LogMapApplyOffset(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Ioffset,
                         IMAGE *Idst) ;
IMAGE *LogMapGaussianFilter(LOGMAP_INFO *lmi, IMAGE *Isrc,
                            IMAGE *Igaussian, IMAGE *Idst) ;
IMAGE *LogMapOffset(LOGMAP_INFO *lmi, IMAGE *Isrc, int wsize, IMAGE *Ioffset) ;
IMAGE *LogMapAddSaltNoise(LOGMAP_INFO *lmi, IMAGE *Isrc,IMAGE *Idst,
                          float density);

#if 0
#define for_each_neighbor(lmi, ptr, i, j, r, s) \
    for (ptr = LOG_PIX_NEIGHBORS(lmi, i, j) ; ptr ; ptr = ptr->next)  \
    if (((s = ptr->logpix->spoke)|1) && ((r = ptr->logpix->ring)|1))
#else
#define for_each_neighbor(lmi, ptr, i, j) \
    for (ptr = LOG_PIX_NEIGHBORS(lmi, i, j) ; ptr ; ptr = ptr->next)
#endif

#define for_each_log_pixel(lmi, r, s) \
    for (r = 0; r < lmi->nrings; r++) \
      for (s = 0; s < lmi->nspokes; s++) \
        if (LOG_PIX_AREA(lmi,r,s) > 0)

/* this macro will go through a specified set of rings */
#define for_each_ring(lmi, r, s, start_ring, end_ring) \
    for (r = start_ring; r <= end_ring ; r++) \
      for (s = lmi->start_spoke[r]; s <= lmi->end_spoke[r]; s++)

/* this macro will go through all log pixels, even those with no area */
#define for_all_log_pixels(lmi, r, s) \
    for (s = 0; s < lmi->nspokes; s++) \
      for (r = 0; r < lmi->nrings; r++)

#define for_each_tv_pixel(lmi, i, j) \
    for (j = 0; j < lmi->nrows ; j++)  \
      for (i = 0; i < lmi->ncols ; i++)


#define UNDEFINED       255
#define DEFINED(r, s)   ((r != UNDEFINED) && (s != UNDEFINED))

#define DIFFUSION_TIME_LOG            0  /* just straight diffusion */
#define DIFFUSION_TIME_PERIPHERAL     1  /* base time on peripheral it. */
#define DIFFUSION_TIME_FOVEAL         2  /*     ||        foveal   ||   */
#define DIFFUSION_TIME_CARTESIAN      3  /* run for equivalent cart. it. */


int log_median(LOGMAP_INFO *lmi, IMAGE *hdi, IMAGE *hdo, int size) ;

#define LOGMAP_RLE        0
#define LOGMAP_SAMPLE     1
#define LOGMAP_NONLOCAL   2
#define LOGMAP_FILTER     3


#endif
