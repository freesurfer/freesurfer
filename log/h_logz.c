/*
  @(#)h_logz.c  1.1
  4/4/94
*/
/*----------------------------------------------------------------------
           File Name:
               h_logz.c

           Author:
             Bruce Fischl with algorithms stolen from Rich Wallace.

           Description:

             mRunNoToRing    maps a run # to the ring it resides in.
             mRunNoToSpoke   maps a run # to the ring it resides in.
             mTvToRing       maps cartesian (x, y) to a ring #
             mTvToSpoke      maps cartesian (x, y) to a spoke #
             mLogPixArea     maps a log pixel (ring,spoke) into its area
             mRunLen         contains a list of all run lengths.


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
#include <memory.h>
#include <unistd.h>
#include <hipl_format.h>

#include "const.h"
#include "diag.h"
#include "runfuncs.h"
#include "congraph.h"
#include "h_logz.h"
#include "filter.h"
#include "map.h"
#include "image.h"
#include "error.h"
#include "proto.h"
#include "timer.h"

/*----------------------------------------------------------------------
                           CONSTANTS
----------------------------------------------------------------------*/

#define MAX_RUN_LENGTH 32

/*----------------------------------------------------------------------
                           TYPEDEFS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           STATIC DATA
----------------------------------------------------------------------*/

static fwd_func  forward_func_array[MAX_RUN_LENGTH+1] ;
static inv_func  inverse_func_array[MAX_RUN_LENGTH+1] ;

/*----------------------------------------------------------------------
                           FUNCTION PROTOTYPES
----------------------------------------------------------------------*/

static int lmConvolve1d(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Igaussian, 
                        IMAGE *Idst, int which);
static void  mapCalculateLogmap(LOGMAP_INFO *mi) ;
static void  mapCalculateLogmapQuadrant(LOGMAP_INFO *mi,
                                        CMAP *mTvToSpoke, CMAP *mTvToRing) ;
static void  mapQuadToFull(LOGMAP_INFO *mi, CMAP *mQuadSpoke,
                           CMAP *mQuadRing);
static void  mapInitLogPix(LOGMAP_INFO *mi) ;

static void  mapInitRuns(LOGMAP_INFO *mi, int max_run_len) ;
static void  mapCalculateParms(LOGMAP_INFO *lmi) ;

static int    logFilterNbd(LOGMAP_INFO *lmi, int filter[NBD_SIZE], 
                           IMAGE *inImage, IMAGE *outImage, 
                           int doweight, int start_ring, int end_ring) ;
#define NEW_DIFFUSION 1
#if !NEW_DIFFUSION
static double diffusionCalculateK(LOGMAP_INFO *lmi,IMAGE *image,double k);
#endif
static void logPatchHoles(LOGMAP_INFO *lmi) ;
static void  fnormalize(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst, 
                        float low, float hi) ;
static void  bnormalize(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst, 
                        byte low, byte hi) ;

static void logBuildRhoList(LOGMAP_INFO *lmi) ;
static void logBuildValidSpokeList(LOGMAP_INFO *lmi) ;
static void writeIterations(char *fname, LOGMAP_INFO *lmi, float end_time) ;
static void writeTimes(char *fname, LOGMAP_INFO *lmi, int niter) ;

/*----------------------------------------------------------------------
                           FUNCTIONS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
LogMapFree(LOGMAP_INFO **plmi)
{
  LOGMAP_INFO *lmi ;

  lmi = *plmi ;
  *plmi = NULL ;

  free(lmi->logPix) ;
  free(lmi->runNoToRing) ;
  free(lmi->runNoToSpoke) ;
  free(lmi->runNoToLen) ;
  free(lmi->tvToRing) ;
  free(lmi->tvToSpoke) ;
  free(lmi->rhos) ;
  free(lmi->cmi.pix) ;
  free(lmi) ; 
}
/*----------------------------------------------------------------------
            Parameters:
              alpha  - the alpha parameter for the log(z+alpha) tranform.
              rows   - the # of rows in the TV image.
              cols   - the # of cols in the TV image.

           Description:
----------------------------------------------------------------------*/
LOGMAP_INFO *
LogMapInit(double alpha, int cols, int rows, int nrings, int nspokes)
{
  int       tvsize, logsize ;
  LOGMAP_INFO  *logMapInfo ;

  logMapInfo = (LOGMAP_INFO *)calloc(1, sizeof(LOGMAP_INFO)) ;
  logMapInfo->nrings = nrings ;
  logMapInfo->nspokes = nspokes ;
  logMapInfo->nrows = rows ;
  logMapInfo->ncols = cols ;
  logMapInfo->alpha = alpha ;

  mapCalculateParms(logMapInfo) ;  /* may change alpha and/or nrings */
  nrings = logMapInfo->nrings ;
  alpha = logMapInfo->alpha ;

  tvsize = rows * cols ;
  logsize = nrings * nspokes ;
  logMapInfo->logPix = 
    (LOGPIX *)calloc(logsize, sizeof(*(logMapInfo->logPix))) ;
  logMapInfo->runNoToRing = 
    (int *)calloc(tvsize,sizeof(*(logMapInfo->runNoToRing)));
  logMapInfo->runNoToSpoke =
    (int *)calloc(tvsize,sizeof(*(logMapInfo->runNoToSpoke)));
  logMapInfo->runNoToLen = 
    (char *)calloc(tvsize,sizeof(*(logMapInfo->runNoToLen)));

  logMapInfo->tvToRing = 
    (UCHAR *)calloc(tvsize, sizeof(*(logMapInfo->tvToRing))) ;
  logMapInfo->tvToSpoke = 
    (UCHAR *)calloc(tvsize, sizeof(*(logMapInfo->tvToSpoke))) ;

  /* calculate the logmap lookup tables */
  mapCalculateLogmap(logMapInfo) ;
  mapInitRuns(logMapInfo, MAX_RUN_LENGTH) ;
  mapInitLogPix(logMapInfo) ;
  logMapInfo->ring_fovea = TV_TO_RING(logMapInfo, cols/2, rows/2) ;

  logPatchHoles(logMapInfo) ;
  ConGraphInit(logMapInfo) ;   /* initialize the connectivity graph */
  logBuildRhoList(logMapInfo) ;
  logBuildValidSpokeList(logMapInfo) ;
  LogMapInitForwardFilter(logMapInfo, LMI_FORWARD_FILTER_CIRCLE) ;

  return(logMapInfo) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              use lookup tables and runlengths to convert the cartesian
              image into a logmap image.
              
              for each row in the cartesian image
                find the ring and spoke that that pixel corresponds to
                using the mRunNoToSpoke and mRunNoToRing lookup tables.
                add the entire run of pixels whose length is stored in
                mRunLen.

               mRunLen, mRunNoToSpoke, and mRunNoToRing are all parallel 
               maps.

              mRunLen maps the run # into the length of the run.
              mRunNoToSpoke maps a run # into a spoke coordinate for that 
              run
              mRunNoToRing  maps a run # into a ring  coordinate for that 
              run

              Note that runs always terminate at the end of a row.
----------------------------------------------------------------------*/
/* just use sampling */
IMAGE *
LogMapSample(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst)
{
  int   ring, spoke, crow, ccol ;
  IMAGE *Iout ;
  float  fval ;
  byte   bval ;

  if (!Idst)
    Idst = ImageAlloc(lmi->nspokes,lmi->nrings, Isrc->pixel_format,1);
  else
    ImageClearArea(Idst, 0, 0, Idst->rows, Idst->cols, 0.0f) ;

  if (Idst->pixel_format != Isrc->pixel_format)
    Iout = ImageAlloc(lmi->nspokes, lmi->nrings, Isrc->pixel_format,1);
  else
    Iout = Idst ;

  for_each_log_pixel(lmi, ring, spoke)
  {
    crow = LOG_PIX_ROW_CENT(lmi, ring, spoke) ;
    ccol = LOG_PIX_COL_CENT(lmi, ring, spoke) ;
    if (crow != UNDEFINED && ccol != UNDEFINED)
    {
      switch (Isrc->pixel_format)
      {
      case PFFLOAT:
        fval = *IMAGEFpix(Isrc, ccol, crow) ;
        *IMAGEFpix(Iout, ring, spoke) = fval ;
        break ;
      case PFBYTE:
        bval = *IMAGEpix(Isrc, ccol, crow) ;
        *IMAGEpix(Iout, ring, spoke) = bval ;
        break ;
      default:
        ErrorReturn(Idst, (ERROR_UNSUPPORTED, 
                         "LogMapForward: unsupported pixel format %d",
                         Iout->pixel_format)) ;
        break ;
      }
    }
  }

  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
LogMapForward(LOGMAP_INFO *lmi, IMAGE *tvImage, IMAGE *logImage)
{
  int   j;
  int   area, ringno, spokeno;
  int   nspokes, nrings ;
  int   *spoke, *ring, sring, sspoke ;
  char  *runl ;
  UCHAR *tvEndRowPtr, *tvPtr ;
  IMAGE *Itmp, *Itmp2, *Iltmp, *Ilsave = NULL ;

  if (tvImage->pixel_format != PFBYTE)
  {
    Itmp2 = ImageAlloc(tvImage->rows, tvImage->cols, PFFLOAT, 1) ;
    Itmp = ImageAlloc(tvImage->rows, tvImage->cols, PFBYTE, 1) ;
    ImageScale(tvImage, Itmp2, 0.0f, 255.0f) ;
    ImageCopy(Itmp2, Itmp) ;
    ImageFree(&Itmp2) ;
    tvImage = Itmp ;
  }
  else
    Itmp = NULL ;

  if (logImage->pixel_format != PFFLOAT)
  {
    Iltmp = ImageAlloc(lmi->nspokes, lmi->nrings, PFFLOAT, 1) ;
    Ilsave = logImage ;
    logImage = Iltmp ;
    fprintf(stderr, "LogMapForward: unsupported log pixel format %d\n",
            logImage->pixel_format) ;
    exit(2) ;
  }
  else
    Iltmp = NULL ;

  spoke = &RUN_NO_TO_SPOKE(lmi, 0) ;
  ring  = &RUN_NO_TO_RING(lmi, 0) ;
  nspokes = lmi->nspokes ;
  nrings = lmi->nrings ;
  runl  = &RUN_NO_TO_LEN(lmi, 0) ;

  /* initialize image values */
  for (spokeno = 0 ; spokeno < nspokes ; spokeno++) 
    for (ringno = 0; ringno < nrings; ringno++) 
    {
      if (LOG_PIX_AREA(lmi, ringno, spokeno) > 0)
        *IMAGEFpix(logImage, ringno, spokeno) = 0.0f ;
      else
        *IMAGEFpix(logImage, ringno, spokeno) = 0.0f * UNDEFINED ;
    }

  /* for each row, use run length tables to build logmap image */
  for (j = 0; j < tvImage->rows; j++) 
  {
    tvPtr = IMAGEpix(tvImage, 0, j) ;
    tvEndRowPtr = tvPtr + tvImage->cols ;

    do 
    {
      /* starting ring and spoke of run */
      sring = *ring ;
      sspoke = *spoke ;

      if (sring < 0 || sspoke < 0 || sring >= nrings || sspoke >= nspokes)
        tvPtr += *runl++ ;
      else               /* sum next *runl tv pixels and add to logmap pixel */
        *IMAGEFpix(logImage, sring, sspoke) += 
          (float)(*forward_func_array[(int)(*runl++)])(&tvPtr) ;

      ring++ ;
      spoke++ ;
    }
    while (tvPtr < tvEndRowPtr); 
  }

  for (spokeno = 0; spokeno < nspokes; spokeno++) 
    for (ringno = 0; ringno < nrings; ringno++) 
    {
      if ((area = LOG_PIX_AREA(lmi, ringno, spokeno)) > 0) 
      {
        /* normailze by area */
        *IMAGEFpix(logImage, ringno, spokeno) /= (float)area;  

#if 0
        /* an attempt to solve the colormap bug */
        if (*IMAGEIpix(logImage, ringno, spokeno) < 16) 
          *IMAGEIpix(logImage, ringno, spokeno) = 16;
        if (*IMAGEIpix(logImage, ringno, spokeno) >= 240) 
          *IMAGEIpix(logImage, ringno, spokeno) = 240;
#endif
      }
    }

  LogMapPatchHoles(lmi, tvImage, logImage) ;

  if (Itmp)            /* input image was not in proper format */
    ImageFree(&Itmp) ;

  if (Iltmp)   /* output image was not in proper format */
  {
    ImageCopy(Iltmp, Ilsave) ;
    ImageFree(&Iltmp) ;
  }
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
LogMapInverse(LOGMAP_INFO *lmi, IMAGE *Ilog, IMAGE *Idst)
{
  int   j, ringno, spokeno ;
  int   *spoke, *ring ;
  char  *runl;
  byte  *tvEndRowPtr, *tvPtr, logval ;
  IMAGE *Iin, *Iout, *Itmp ;

  if (Ilog->pixel_format != PFINT)
  {
    Iin = ImageAlloc(Ilog->rows, Ilog->cols, PFINT, 1) ;
    Itmp = ImageScale(Ilog, NULL, 0.0f, (float)(UNDEFINED-1)) ;
    ImageCopy(Itmp, Iin) ;
    ImageFree(&Itmp) ;
  }
  else
    Iin = Ilog ;

  if (Idst->pixel_format != PFBYTE)
    Iout = ImageAlloc(Idst->rows, Idst->cols, PFBYTE, 1) ;
  else
    Iout = Idst ;

  spoke = &RUN_NO_TO_SPOKE(lmi, 0) ;
  ring  = &RUN_NO_TO_RING(lmi, 0) ;
  runl  = &RUN_NO_TO_LEN(lmi, 0) ;
  for (j = 0; j < Iout->rows; j++) 
  {
    tvPtr = IMAGEpix(Iout, 0, j) ;
    tvEndRowPtr = tvPtr + Iout->cols ;
    do 
    {
      ringno = *ring++ ;
      spokeno = *spoke++ ;
      if (ringno >= 0)
        logval = (UCHAR)*IMAGEIpix(Iin, ringno, spokeno) ;
      else
        logval = (UCHAR)UNDEFINED ;

      /* put logval into next *runl tv pixels */
      (*inverse_func_array[(int)(*runl++)])(&tvPtr, logval) ;
    }
    while (tvPtr < tvEndRowPtr); 
  }

  if (Iin != Ilog)
    ImageFree(&Iin) ;

  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }

  return(0) ;
}

/*----------------------------------------------------------------------
                           STATIC FUNCTIONS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:
              rows     - the # of rows in the tv image
              cols     - the # of cols in the cartesian image
              nrings   - the # of rings in the logmap image.
              nspokes  - the # of spokes in the logmap image.
              mTvToSpoke   - one quadrant spoke map.
              mTvToRing    - one quadrant ring map.

           Description:
             calculate the logmap for one quadrant then use symmetry
             to construct the full logmap.

             The quadrant we will build is not any quadrant that will
             appear in the full image.  It is inverted in terms of ring
             numbering (i.e. ring 0 is innermost, nrings/2 is outermost),
             and in terms of spoke numbering (the 90 degree spoke is 
             0, and the 0 degree spoke is nspokes/2).  These coordinates
             will be accounted for in the building of the full log map.
----------------------------------------------------------------------*/
static void  
mapCalculateLogmap(LOGMAP_INFO *lmi)
{
  CMAP *mQuadRing, *mQuadSpoke ;

  /* temporary maps for one quadrant of rings and spokes */
  mQuadRing = MapCAlloc(lmi->ncols/2, lmi->nrows/2) ;
  mQuadSpoke = MapCAlloc(lmi->ncols/2, lmi->nrows/2) ;
  
  /* calculate one quadrant of lookup tables */
  mapCalculateLogmapQuadrant(lmi, mQuadSpoke, mQuadRing) ;

  /* use symmetry to build full mapping */
  mapQuadToFull(lmi, mQuadSpoke, mQuadRing);

  if (Gdiag & DIAG_LOGMAP)   /* write maps to hips file, and show on screen */
  {
    CMAP map ;

    if (Gdiag & DIAG_LOG_QUAD)
    {
      MapCHipsWrite(mQuadRing, "qring.hipl") ;
      MapCHipsWrite(mQuadSpoke, "qspoke.hipl") ;
      MapCShowImage(mQuadRing, "quad ring") ;
      MapCShowImage(mQuadSpoke, "quad spoke") ;
    }

    map.rows = lmi->nrows ;
    map.cols = lmi->ncols ;
    map.image = lmi->tvToRing ;
    MapCHipsWrite(&map, "ring.hipl") ;
    MapCShowImage(&map, "ring") ;

    map.image = lmi->tvToSpoke ;
    MapCHipsWrite(&map, "spoke.hipl") ;
    MapCShowImage(&map, "spoke") ;
  }
  
  MapCFree(mQuadRing) ;
  MapCFree(mQuadSpoke) ;
}
/*----------------------------------------------------------------------
            Parameters:
              rows     - the # of rows in the tv image
              cols     - the # of cols in the cartesian image
              nrings   - the # of rings in the logmap image.
              nspokes  - the # of spokes in the logmap image.
              mTvToSpoke   - one quadrant spoke map.
              mTvToRing    - one quadrant ring map.

           Description:
             calculate the Cartesian to logmap transformation map for 
             one quadrant.

             The quadrant we will build is not any quadrant that will
             appear in the full image.  It is inverted in terms of ring
             numbering (i.e. ring 0 is innermost, nrings/2 is outermost),
             and in terms of spoke numbering (the 90 degree spoke is 
             0, and the 0 degree spoke is nspokes/2).  These coordinates
             will be accounted for in the building of the full log map.
----------------------------------------------------------------------*/
void
mapCalculateLogmapQuadrant(LOGMAP_INFO *lmi, CMAP *mQuadTvToSpoke, 
                           CMAP *mQuadTvToRing) 
{
  double r, x, y, theta, maxlog, minlog, minangle, alpha ; 
  int    ringnum, spokenum, halfrows, halfcols, nrings, nspokes ;
  int    i, j;
  
  halfcols = lmi->ncols / 2 ;
  halfrows = lmi->nrows / 2 ;
  alpha = lmi->alpha ;
  nspokes = lmi->nspokes ;
  nrings = lmi->nrings ;
  
  /* initialize arrays to 0 */
  for (i = 0; i < halfcols ; i++)
    for (j = 0; j < halfrows ; j++) 
      MAPcell(mQuadTvToSpoke, i, j) = MAPcell(mQuadTvToRing, i, j) = 0;
  
  /*
    the min log occurs along the x axis 
    but the max log is along the y axis
    because there are fewer rows than cols 
    */
/* 
  The maximum r comes along the shorter axis so as to contain the outermost
  ring within the Cartesian image, so test rows and cols, and calculate max
  r accordingly - BRF 7/21/93
*/
#if 0
  if (halfrows < halfcols)
    r = hypot((double)halfrows+0.5, alpha+0.5); 
  else
    r = hypot((double)0.5, alpha+0.5+halfcols) ;
#else
  r = MIN(hypot((double)halfrows+0.5, alpha+0.5), 
          hypot((double)0.5, alpha+0.5+halfcols)) ;
#endif

  maxlog = log(r);                           /* max is along the y axis */
  r = hypot(alpha+0.5, 0.5);                 /* min is along the x axis */
  minlog = log(r); 


  /* 
    x and y are inverted here because of the bizarre nature of the
    quadrant we are building (see function header).
    */
  minangle = atan2(alpha+0.5, halfrows-0.5); /* min angle is along y axis */
  
  for (i = 0; i < halfcols; i++) 
  {
    for (j = 0; j < halfrows; j++) 
    {
      /* 
        offset by 1/2 so map computed at pixel centers.  This eliminates the
        singularity at the origin even if alpha = 0.
        origin in upper left corner 
        x ranges from alpha to alpha+halfcols-1+0.5 
        y ranges from 0 to halfrows-1+0.5 
        */
      x = alpha + 0.5 + (double)i;  
      
      y = 0.5+(double)j;  
      r = hypot(y, x);
      if (r > 1.0) 
      {
        /*
          normalize ring number so that the maximum ring will be
          nrings/2, and the minimum ring will be 0.
          */
        ringnum = 
          (int)((double)(nrings/2)*(log(r)-minlog) /(maxlog-minlog));

        /* clip square image to circular region */
        if (ringnum >= nrings/2) 
          ringnum = UNDEFINED;

        /* 
          x and y are inverted here because of the bizarre nature of the
          quadrant we are building (see function header).
          */
        theta = atan2(x, y);

        /*
          normalize spoke numbers so that theta = PI/2 gives maximum
          spoke # nspokes/2.  Recall that x and y were flipped in the
          calculation of theta.  Therefore, theta = PI/2 corresponds
          to the spoke along the x-axis, and theta = minangle corresponds
          to the spoke along the y-axis.
          */
        spokenum = 
          (int)(((double)(nspokes/2)*(theta-minangle))/(PI/2.0-minangle));
        if (ringnum == UNDEFINED) 
          spokenum = UNDEFINED; /* make both undefined */
      }
      else 
      {
        ringnum = 0;
        spokenum = 0;
      }
      MAPcell(mQuadTvToSpoke, i, j) = spokenum;
      MAPcell(mQuadTvToRing, i, j) = ringnum;
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:
              mQuadSpoke  - the quadrant spoke map
              mQuadRing   - the quadrant ring map
              mTvToSpoke      - the full spoke map
              mTvToRing       - the full ring map

           Description:
              use symmetry to calculate the full logmap from one
              quadrant.

              Note the coordinate systems are somewhat bizarre.  (0, 0)
              in the full logmap is in the lower left hand corner.



(0,rows-1)                                            (cols-1,rows-1)
------------------------------------------------------------------
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |(halfcols, halfrows)           |
------------------------------------------------------------------
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
|                                 |                               |
------------------------------------------------------------------
(0,0)                                                            (cols-1,0)




----------------------------------------------------------------------*/
static void
mapQuadToFull(LOGMAP_INFO *lmi, CMAP *mQuadSpoke, CMAP *mQuadRing)
{
  int           i, j, halfrows, halfcols, nrings, nspokes ;
  unsigned char s, r;
  
  nrings = lmi->nrings ;
  nspokes = lmi->nspokes ;
  halfrows = mQuadSpoke->rows ;
  halfcols = mQuadSpoke->cols;
  
  for (i = 0; i < halfcols; i++) 
  {
    for (j = 0; j < halfrows; j++) 
    {
      /* calculate full spoke map */
      s = MAPcell(mQuadSpoke, i, j); 
      if (s != UNDEFINED) 
        s = MAPcell(mQuadSpoke, i, j); /* redundant */

      /* fill in lower left quadrant ( halfcols - 0, halfrows - 0) */
      TV_TO_SPOKE(lmi, halfcols-1-i, halfrows-1-j) = s;

      /* fill in lower right quadrant (halfcols - cols, halfrows - 0) */
      TV_TO_SPOKE(lmi, halfcols+i, halfrows-1-j) = s;

      /* 
        now fill in top two quadrants, mapping half spoke #s to full
         spoke #s.
         */
      s = MAPcell(mQuadSpoke, i, j); 
      if (s != UNDEFINED) 
        s = nspokes-1-MAPcell(mQuadSpoke, i, j);   /* full spoke # */

      /* fill in upper left quadrant (halfcols - 0, halfrows to rows) */
      TV_TO_SPOKE(lmi, halfcols-1-i, halfrows+j) = s;

      /* fill in upper right quadrant (halfcols - cols, halfrows - rows) */
      TV_TO_SPOKE(lmi, halfcols + i, halfrows+j)= s;

      /* calculate full ring map */
      r = MAPcell(mQuadRing, i, j);

      /*
        in the quad map, 0 is the innermost ring, in the full logmap
        0 is the outermost ring in the left half plane.
        */

      /* calculate left half plane (0 - (nrings/2-1)) */

      /* if it is a real ring, reflect it around nrings/4 so quad map
         0 maps to full map nrings/2 -1, and quad map nrings/2 - 1
         maps to 0.
         */
      if (r != UNDEFINED) 
        r = nrings/2 - 1 - MAPcell(mQuadRing, i, j);

      /* lower left quadrant (halfcols - 0, halfrows - 0) */
      TV_TO_RING(lmi, halfcols - 1 - i, halfrows-1-j)= r;

      /* upper left quadrant (halfcols - 0, halfrows - rows) */
      TV_TO_RING(lmi, halfcols - 1 - i, halfrows+j)=   r;

      /* calculate right half plane (nrings/2 - nrings) */
      r = MAPcell(mQuadRing, i, j);

      /*
        reflect quad map ring # around nrings/2 so quad ring 0 maps to
        fullmap ring # nrings/2, and quad map nrings/2 maps to full map
        nrings.
        */
      if (r != UNDEFINED) 
        r = nrings/2 + MAPcell(mQuadRing, i, j);

      /* calculate lower right quadrant (halfcols - cols, halfrows - 0) */
      TV_TO_RING(lmi, halfcols + i, halfrows-1-j)=  r;

      /* calculate upper right quadrant (halfcols - cols, halfrows - rows) */
      TV_TO_RING(lmi, halfcols + i, halfrows+j)= r;
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:
              rows - the # of rows in the tv image.
              cols - the # of cols in the tv image.

           Description:
              initialize the run length lookup tables.
----------------------------------------------------------------------*/
static void
mapInitRuns(LOGMAP_INFO *lmi, int max_run_len)
{
  int i, j, nrings, nspokes ;
  int numruns ;
  int r_index, s_index ;
  int lr_index, ls_index ; 
  int runLen ;
  int runLenSum, rows, cols ;

  nrings = lmi->nrings ;
  nspokes = lmi->nspokes ;
  runLenSum = 0;
  numruns = 0;
  rows = lmi->nrows ;
  cols = lmi->ncols ;

  for (j = 0; j < rows; j++) 
  {
    i = 0;
    s_index = TV_TO_SPOKE(lmi, i, j) ;
    r_index = TV_TO_RING(lmi, i, j) ;
    ls_index = s_index;   /* left spoke index (at start of run) */
    lr_index = r_index;   /* left ring  index (at start of run) */
    runLen = 1;

    for (i = 1; i < cols; i++) 
    {
      s_index = TV_TO_SPOKE(lmi, i, j) ;   /* index of current spoke */
      r_index = TV_TO_RING(lmi, i, j) ;    /* index of current ring */

      /* 
        check to see if we are done with this run, i.e. if we are no
        longer in the same ring and spoke, or the run length is maximal.
       */
      if (s_index != ls_index || r_index != lr_index ||  runLen == max_run_len)
      {
        runLenSum += runLen;
        /*          printf("%d %d %d\n",ls_index,lr_index,runLen);  */
        RUN_NO_TO_LEN(lmi, numruns)   = runLen;  /* store length of run */
        if (ls_index >= nspokes || lr_index >= nrings) /* out of log image */
        {
          RUN_NO_TO_SPOKE(lmi, numruns) = -1; 
          RUN_NO_TO_RING(lmi, numruns++) = -1;  
        }
        else 
        {
          /* store the spoke and ring indices that this run begins at */
          RUN_NO_TO_SPOKE(lmi, numruns) = ls_index;
          RUN_NO_TO_RING(lmi, numruns++) = lr_index;
        }

        /* start a new run beginning at current position */
        runLen = 0;
        ls_index = s_index;
        lr_index = r_index;
      }
      runLen++;
    }

    /* end of row, terminate run */
    runLenSum += runLen;
    /*     printf("%d %d %d\n",ls_index,lr_index,runLen);    */
    RUN_NO_TO_LEN(lmi, numruns) = runLen;
    if (ls_index >= nspokes || lr_index >= nrings)  /* out of log image */
    {
      RUN_NO_TO_SPOKE(lmi, numruns)   = -1;
      RUN_NO_TO_RING(lmi, numruns++) = -1;
    }
    else 
    {
      RUN_NO_TO_SPOKE(lmi, numruns) = ls_index;
      RUN_NO_TO_RING(lmi, numruns++) = lr_index;
    }
  }
  lmi->nruns = numruns ;

#if 0
  printf("num runs = %d\n",numruns);
  MapCHipsWrite(mRunLen, "runlen.hipl") ;
  MapIHipsWrite(mRunNoToSpoke, "spokerun.hipl") ;
  MapIHipsWrite(mRunNoToRing, "ringrun.hipl") ;
#endif
    
  runFuncInit(forward_func_array, inverse_func_array) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
              Initialize the the LOGPIX entry for each logmap pixel.
              This involves computing area, centroids, neighbors, etc....

              AREA:
              Count the # of cartesian pixels that fall within the
              boundaries of each logmap pixel.
----------------------------------------------------------------------*/
#define DEBUG_WEIGHT 0
static void
mapInitLogPix(LOGMAP_INFO *lmi)
{
  int     i, j, ring, spoke, area, rval, nrings, rows, cols ;
  double  weight, max_ring, min_ring, ring_step, log_range, nrings_m1, 
         wscale ;
  float  rho, phi, x, y, min_rho, max_rho ;
  int    halfrows, halfcols ;
#if DEBUG_WEIGHT
  FILE    *fp ;
#endif

  min_rho = 1000.0f ;
  max_rho = 0.0f ;
  nrings = lmi->nrings ;
  for_each_tv_pixel(lmi, i, j)
  {
    spoke = TV_TO_SPOKE(lmi, i, j) ;
    ring = TV_TO_RING(lmi, i, j) ;

    if ((spoke != UNDEFINED) && (ring != UNDEFINED))
    {
      LOG_PIX_AREA(lmi, ring, spoke)++ ;
      LOG_PIX_ROW_CENT(lmi, ring, spoke) += j ;
      LOG_PIX_COL_CENT(lmi, ring, spoke) += i ;
      LOG_PIX_SPOKE(lmi, ring, spoke) = spoke ;
      LOG_PIX_RING(lmi, ring, spoke) = ring ;
    }
  }

  nrings_m1 = (double)lmi->nrings-1.0 ;
  min_ring = log(lmi->alpha) ;
  max_ring = log(lmi->maxr) ;
  log_range = max_ring - min_ring ;
  ring_step = log_range / nrings_m1 ;
  wscale = 1.0 / exp(-log(lmi->alpha)) ;
wscale = 1.0 ;
#if 0
  fprintf(stderr, 
          "maxr %2.3f, min_ring %2.2f, max_ring %2.2f, ring_step %2.3f\n",
          lmi->maxr, min_ring, max_ring, ring_step) ;
#endif

#if DEBUG_WEIGHT
  fp = fopen("weights.dat", "w") ;
#endif

  rows = lmi->nrows ;
  cols = lmi->ncols ;
  for_each_log_pixel(lmi, ring, spoke)
  {
    area = LOG_PIX_AREA(lmi, ring, spoke) ;
    LOG_PIX_ROW_CENT(lmi, ring, spoke) /= area ;
    LOG_PIX_COL_CENT(lmi, ring, spoke) /= area ;
    if (ring <= nrings/2)
      rval = nrings/2-ring ;
    else
      rval = ring - nrings/2 ;

    /* don't left centroid be on border of image -- BRF for nonlocal stuff */
    if (LOG_PIX_ROW_CENT(lmi,ring,spoke) == 0)
      LOG_PIX_ROW_CENT(lmi,ring,spoke) = 1 ;
    else if (LOG_PIX_ROW_CENT(lmi,ring,spoke) == rows-1)
      LOG_PIX_ROW_CENT(lmi,ring,spoke) = rows-2 ;
    else if (LOG_PIX_COL_CENT(lmi,ring,spoke) == 0)
      LOG_PIX_COL_CENT(lmi,ring,spoke) = 1 ;
    else if (LOG_PIX_COL_CENT(lmi,ring,spoke) == cols-1)
      LOG_PIX_COL_CENT(lmi,ring,spoke) = cols-2 ;

#if 0
    weight = (double)rval * ring_step + min_ring ;
#else

    halfrows = lmi->nrows / 2 ;
    halfcols = lmi->ncols / 2 ;
    
    x = (float)(LOG_PIX_COL_CENT(lmi, ring, spoke) - halfcols) ;
    y = (float)(LOG_PIX_ROW_CENT(lmi, ring, spoke) - halfrows);
    if (x < 0)
      x = -x ;

    /* not quite right, but good enough for weight */
    LOG_PIX_RHO(lmi,ring,spoke) = rho = log(sqrt(SQR(x+lmi->alpha)+SQR(y)));
    LOG_PIX_PHI(lmi,ring,spoke) = phi = atan2(y, x+lmi->alpha) ;
    weight = exp(-rho) ;
    
    if (rho > max_rho)
      max_rho = rho ;
    if (rho < min_rho)
      min_rho = rho ;
    
#endif
#if DEBUG_WEIGHT
    if (spoke == lmi->nspokes/2)
      fprintf(fp, "%d  %d  rho %2.3f  weight %2.3f, cent (%2.0f, %2.0f)\n", 
              ring, spoke, rho, weight, x, y) ;
#endif
    LOG_PIX_WEIGHT(lmi, ring, spoke) = weight ;
  }

#if DEBUG_WEIGHT
  fclose(fp) ;
#endif
  fprintf(stderr, "rho: %2.3f --> %2.3f\n", min_rho, max_rho) ;
  lmi->min_rho = min_rho ;
  lmi->max_rho = max_rho ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int isy[NBD_SIZE] = { -1, -2, -1,   0, 0, 0,   1, 2, 1} ;
static int isx[NBD_SIZE] = {-1, 0, 1,  -2, 0, 2,   -1,  0,  1} ;

int
LogMapGradient(LOGMAP_INFO *lmi, IMAGE *inImage, 
               IMAGE *gradImage, int doweight, int start_ring, int end_ring)
{
  static IMAGE *xImage = NULL, *yImage = NULL ;

  if (!ImageCheckSize(inImage, xImage, 0, 0, 0))
  {
    if (xImage)
    {
      ImageFree(&xImage) ;
      ImageFree(&yImage) ;
    }
    xImage = ImageAlloc(inImage->rows, inImage->cols,inImage->pixel_format,1);
    yImage = ImageAlloc(inImage->rows, inImage->cols, inImage->pixel_format,1);
  }
  else
  {
    ImageSetSize(xImage, inImage->rows, inImage->cols) ;
    ImageSetSize(yImage, inImage->rows, inImage->cols) ;
    yImage->pixel_format = xImage->pixel_format = inImage->pixel_format ;
  }

  LogMapSobel(lmi, inImage, gradImage, xImage, yImage, doweight, 
              start_ring, end_ring);
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
LogMapSobel(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *gradImage, 
            IMAGE *Ix, IMAGE *Iy, int doweight, int start_ring, int end_ring)
{
  int              ring, spoke, val = 0, xval = 0, yval = 0 ;
  float            fval = 0.0f, fxval = 0.0f, fyval = 0.0f ;

  if (start_ring < 0)
    start_ring = 0 ;
  if (end_ring <= 0)
    end_ring = lmi->nrings-1 ;

  if (!ImageCheckSize(Isrc, Ix, 0, 0, 0))
    ErrorReturn(-1, (ERROR_NO_MEMORY, "LogMapSobel: x image not big enough"));
  if (!ImageCheckSize(Isrc, Iy, 0, 0, 0))
    ErrorReturn(-1, (ERROR_NO_MEMORY, "LogMapSobel: y image not big enough"));

  ImageSetSize(Ix, Isrc->rows, Isrc->cols) ;
  ImageSetSize(Iy, Isrc->rows, Isrc->cols) ;
  Iy->pixel_format = Ix->pixel_format = Isrc->pixel_format ;

  logFilterNbd(lmi, isx, Isrc, Ix, doweight, start_ring, end_ring) ;
  logFilterNbd(lmi, isy, Isrc, Iy, doweight, start_ring, end_ring) ;

  if (gradImage)
  {
    switch (Isrc->pixel_format)
    {
    case PFINT:
      for_each_ring(lmi, ring, spoke, start_ring, end_ring)
        {
          xval = *IMAGEIpix(Ix, ring, spoke) ;
          yval = *IMAGEIpix(Iy, ring, spoke) ;
          val = (int)sqrt(((double)(xval*xval) + (double)(yval*yval))/4.0) ;
          *IMAGEIpix(gradImage, ring, spoke) = val ;
        }
      break ;
    case PFFLOAT:
      for_each_ring(lmi, ring, spoke, start_ring, end_ring)
        {
          fxval = *IMAGEFpix(Ix, ring, spoke) ;
          fyval = *IMAGEFpix(Iy, ring, spoke) ;
          fval = (float)sqrt((double)(fxval*fxval)+(double)(fyval*fyval))/4.0f;
          *IMAGEFpix(gradImage, ring, spoke) = fval ;
        }
      break ;

    default:
      fprintf(stderr, "LogMapGradient: unsupported log pixel format %d\n",
              Isrc->pixel_format) ;
      exit(2) ;
      break ;
    }
  }

#if DEBUG_FILTER
  if (Gdiag & DIAG_WRITE)
  {
    ImageWrite(Ix, "x.hipl") ;
    ImageWrite(Iy, "y.hipl") ;
    ImageWrite(gradImage, "grad.hipl") ;
  }
#endif
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int mean_kernel[NBD_SIZE] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 } ;

IMAGE *
LogMapMeanFilter(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst)
{
  if (!ImageCheckSize(Isrc, Idst, 0, 0, 0))
  {
    if (Idst)
      ImageFree(&Idst) ;
    Idst = ImageAlloc(Isrc->rows,Isrc->cols,Isrc->pixel_format,
                      Isrc->num_frame);
  }

  logFilterNbd(lmi, mean_kernel, Isrc, Idst, 0, 0, lmi->nrings-1) ;
  
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
LogMapCurvature(LOGMAP_INFO *lmi, IMAGE *inImage, IMAGE *gradImage, 
                float A, int doweight, int start_ring, int end_ring)
{
  int    ring, spoke, val, xval, yval ;
  double Asq, Imag ;
  float  fxval, fyval, fval ;
  static IMAGE *xImage = NULL, *yImage = NULL ;

  if (start_ring < 0)
    start_ring = 0 ;
  if (end_ring <= 0)
    end_ring = lmi->nrings-1 ;

  if (!ImageCheckSize(inImage, xImage, 0, 0, 0))
  {
    if (xImage)
    {
      ImageFree(&xImage) ;
      ImageFree(&yImage) ;
    }
    xImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1) ;
    yImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1) ;
  }
  else
  {
    ImageSetSize(xImage, inImage->rows, inImage->cols) ;
    ImageSetSize(yImage, inImage->rows, inImage->cols) ;
  }

  logFilterNbd(lmi, isx, inImage, xImage, doweight, start_ring, end_ring) ;
  logFilterNbd(lmi, isy, inImage, yImage, doweight, start_ring, end_ring) ;

  Asq = (double)(A * A) ;
  switch (inImage->pixel_format)
  {
  case PFINT:
    for_each_ring(lmi, ring, spoke, start_ring, end_ring)
    {
      xval = *IMAGEIpix(xImage, ring, spoke) ;
      yval = *IMAGEIpix(yImage, ring, spoke) ;
      Imag = (double)((xval * xval) + (yval * yval)) ;
      val = (int)sqrt(1.0 + Asq * Imag) ;
      *IMAGEIpix(gradImage, ring, spoke) = val ;
    }
    break ;
  case PFFLOAT:
    for_each_ring(lmi, ring, spoke, start_ring, end_ring)
    {
      fxval = *IMAGEFpix(xImage, ring, spoke) / 4.0f ;
      fyval = *IMAGEFpix(yImage, ring, spoke) / 4.0f ;
      Imag = (double)((fxval*fxval) + (fyval*fyval)) ;
      fval = (float)sqrt(1.0 + Asq * Imag);
      *IMAGEFpix(gradImage, ring, spoke) = fval ;
    }
    break ;
  default:
    fprintf(stderr, "LogMapCurvature: unsupported pixel format %d\n",
            gradImage->pixel_format) ;
    exit(3) ;
    break ;
  }

#if DEBUG_FILTER
  if (Gdiag & DIAG_WRITE)
  {
    ImageWrite(xImage, "x.hipl") ;
    ImageWrite(yImage, "y.hipl") ;
    ImageWrite(gradImage, "grad.hipl") ;
  }
#endif

  return(0) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
logFilterNbd(LOGMAP_INFO *lmi, int filter[NBD_SIZE], IMAGE *inImage, 
             IMAGE *outImage, int doweight, int start_ring, int end_ring)
{
  register int     k, ring, spoke, val, n_ring, n_spoke, *pf ;
  register float   fval ;
  LOGPIX  *npix ;

  if (outImage->pixel_format != inImage->pixel_format)
    ErrorReturn(-1, 
                (ERROR_UNSUPPORTED,
                 "logFilterNbd: input and output format must be the same\n"));

  switch (inImage->pixel_format)
  {
  case PFINT:
    for_each_ring(lmi, ring, spoke, start_ring, end_ring)
    {
      val = 0 ;
      pf = filter ;
      for (k = 0 ; k < NBD_SIZE ; k++)
      {
        npix = LOG_PIX_NBD(lmi, ring, spoke, k) ;
        n_ring = npix->ring ;
        n_spoke = npix->spoke ;
        val += *IMAGEIpix(inImage, n_ring, n_spoke) * *pf++ ;
      }
      if (doweight)
        *IMAGEIpix(outImage, ring, spoke) = 
          val * LOG_PIX_WEIGHT(lmi,ring,spoke);
      else
        *IMAGEIpix(outImage, ring, spoke) = val ;
    }
    break ;
  case PFFLOAT:
    for_each_ring(lmi, ring, spoke, start_ring, end_ring)
    {
      fval = 0.0f ;
      pf = filter ;
      for (k = 0 ; k < NBD_SIZE ; k++)
      {
        npix = LOG_PIX_NBD(lmi, ring, spoke, k) ;
        n_ring = npix->ring ;
        n_spoke = npix->spoke ;
        fval += *IMAGEFpix(inImage, n_ring, n_spoke) * *pf++ ;
      }
      if (doweight)
        *IMAGEFpix(outImage, ring,spoke) = 
          fval * LOG_PIX_WEIGHT(lmi,ring,spoke);
      else
        *IMAGEFpix(outImage, ring, spoke) = fval ;

    }
    break ;
  default:
    ErrorReturn(-2, 
                (ERROR_UNSUPPORTED,
                 "logFilterNbd: unsupported input image format %d\n",
                 inImage->pixel_format)) ;
    break ;   /* never used */
  }
  return(0) ;
}


/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
double
LogMapDiffuse(LOGMAP_INFO *lmi, IMAGE *inImage, IMAGE *outImage, double k, 
              int niterations, int doweight, int which, int time_type)
{
  static  IMAGE *fSrcImage = NULL, *fDstImage = NULL ;
  IMAGE         *fOut ;

  if (!ImageCheckSize(inImage, fSrcImage, 0, 0, 0))
  {
    if (fSrcImage)
    {
      ImageFree(&fSrcImage) ;
      ImageFree(&fDstImage) ;
    }
    fSrcImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1) ;
    fDstImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1) ;
  }
  else
  {
    ImageSetSize(fDstImage, inImage->rows, inImage->cols) ;
    ImageSetSize(fSrcImage, inImage->rows, inImage->cols) ;
  }

#if DEBUG  
  if (Gdiag & DIAG_TIMER)
    DebugTimerStart() ;
#endif

  ImageCopy(inImage, fSrcImage) ;
  if (outImage->pixel_format == PFFLOAT)
    fOut = outImage ;
  else
    fOut = fDstImage ;

  switch (which)
  {
  case DIFFUSE_PERONA:
    k = LogMapDiffusePerona(lmi, fSrcImage, fOut, k, niterations, 
                            doweight, time_type);
    break ;
  case DIFFUSE_CURVATURE:
    k = LogMapDiffuseCurvature(lmi, fSrcImage, fOut,k,niterations,doweight,
                               time_type);
    break ;
  default:
    fprintf(stderr, "LogMapDiffuse: unsupported diffusion type %d\n",which);
    exit(1) ;
  }

  if (fOut != outImage)
    ImageCopy(fOut, outImage) ;

#if 0
  switch (outImage->pixel_format)
  {
  case PFINT:
    for (spoke = 0 ; spoke < lmi->nspokes ; spoke++) 
      for (ring = 0; ring < lmi->nrings; ring++) 
      {
        if (LOG_PIX_AREA(lmi, ring, spoke) <= 0)
          *IMAGEIpix(outImage, ring, spoke) = 0 * UNDEFINED ;
      }
    break ;
  case PFFLOAT:
    for (spoke = 0 ; spoke < lmi->nspokes ; spoke++) 
      for (ring = 0; ring < lmi->nrings; ring++) 
      {
        if (LOG_PIX_AREA(lmi, ring, spoke) <= 0)
          *IMAGEFpix(outImage, ring, spoke) = 0.0f * UNDEFINED ;
      }
    break ;
  default:
    break ;
  }
#endif

#if DEBUG
  if (Gdiag & DIAG_TIMER)
  {
    char str[100] ;

    DebugTimerStop() ;
    sprintf(str, "diffusion for %d iterations: ", niterations) ;
    DebugTimerShow(str) ;
  }
#endif

  return(k) ;
}

/*
   KERNEL_MUL is 1/2 dt
*/


#define FOUR_CONNECTED 0
#if FOUR_CONNECTED
#define KERNEL_MUL    0.25f
#else
#define KERNEL_MUL    0.125f
#endif
#define C(grad,k)    ((float)exp(-0.5f * SQR((float)fabs((grad))/k)))

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
double
LogMapDiffusePerona(LOGMAP_INFO *lmi, IMAGE *inImage, IMAGE *outImage, 
              double k, int niterations, int doweight, int time_type)

{
  int     ring, spoke = 0, i, rows, cols, n_ring, n_spoke, ci, j, 
          end_ring, start_ring, nspokes, new_start_ring, new_end_ring,
          start_spoke, end_spoke ;
  float   weight, c[NBD_SIZE], fvals[NBD_SIZE], dst_val, rho,time,end_time;
  FILE    *fp ;
  LOGPIX  *pix, **pnpix, *npix ;
  static  IMAGE *tmpImage = NULL, *gradImage = NULL ;
  register float   *pf, *pc, *pcself ;

  rows = inImage->rows ;
  cols = inImage->cols ;

  if ((outImage->pixel_format != inImage->pixel_format) ||
      (outImage->pixel_format != PFFLOAT))
  {
    fprintf(stderr,
        "ImageDiffusePerona: input and output image format must both be "
            "in float format.\n");
    exit(2) ;
  }

  if ((outImage->rows != inImage->rows) ||
      (outImage->cols != inImage->cols))
  {
    fprintf(stderr,
          "ImageDiffusePerona: input and output image sizes must match.\n");
    exit(2) ;
  }

  if (!ImageCheckSize(inImage, tmpImage, 0, 0, 0))
  {
    if (tmpImage)
      ImageFree(&tmpImage) ;
    tmpImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1) ;
  }
  else
    ImageSetSize(tmpImage, inImage->rows, inImage->cols) ;

  if (!ImageCheckSize(inImage, gradImage, 0, 0, 0))
  {
    if (gradImage)
      ImageFree(&gradImage) ;
    gradImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1) ;
  }
  else
    ImageSetSize(gradImage, inImage->rows, inImage->cols) ;

  ImageCopy(inImage, tmpImage) ;

  if (Gdiag & DIAG_WRITE)
    fp = fopen("diffuse.dat", "w") ;
  else
    fp = NULL ;

  end_time = 0.0f ;
  switch (time_type)
  {
  case DIFFUSION_TIME_FOVEAL:
    /* compute ending time for Cartesian space */
    end_time = ((float)niterations+0.9f)*exp(2.0f*lmi->min_rho) *KERNEL_MUL;
    break ;
  case DIFFUSION_TIME_PERIPHERAL:
    /* for peripheral iterations, assume # of iterations is 10 times
       actual number to allow for fractional # of iterations */
    end_time = ((float)niterations + 9.9f) / 10.0f * 
      exp(2.0f*lmi->max_rho) * KERNEL_MUL ;
    niterations = (int)(end_time / (KERNEL_MUL*exp(2.0f*lmi->min_rho))) ;
    break ;
  case DIFFUSION_TIME_CARTESIAN:
    end_time = 2.0f * KERNEL_MUL * (float)niterations ;  
    break ;
  case DIFFUSION_TIME_LOG:
  default:
    end_time = ((float)niterations+.9f)* exp(2.0f*lmi->max_rho) *KERNEL_MUL;
    break ;
  }

  if (Gdiag & DIAG_LOGDIFF)
  {
    if (time_type == DIFFUSION_TIME_LOG)
      writeTimes("etimes.dat", lmi, niterations) ;
    else
      writeIterations("niter.dat", lmi, end_time) ;
    if (Gdiag & DIAG_WRITE)
      fprintf(stderr, "end time = %2.3f, niter = %d\n", end_time,niterations);
  }

  start_ring = 0 ;
  end_ring = lmi->nrings - 1 ;
  nspokes = lmi->nspokes ;
  time = 0.0f ;

#define LOG_MOVIE_FNAME      "log_diffuse_movie.hipl"

  unlink(LOG_MOVIE_FNAME) ;
  ImageAppend(tmpImage, LOG_MOVIE_FNAME) ;
  for (i = 0 ; i < niterations ; i++)
  {
    if (time_type != DIFFUSION_TIME_LOG)
    {
      /* find starting and ending ring for this time */
      new_start_ring = lmi->nrings ;
      new_end_ring = 0 ;

      for (ring = start_ring ; ring <= end_ring ; ring++)
      {
        rho = lmi->rhos[ring] ;
        if (FZERO(rho))
          continue ;   /* no valid pixels in this ring */

        time = (float)i * exp(2.0f*rho) * KERNEL_MUL ;
        if (time <= end_time)  /* this ring still requires integration */
        {
          if (ring < new_start_ring)
            new_start_ring = ring ;
          if (ring > new_end_ring)
            new_end_ring = ring ;
        }
      }
      
      start_ring = new_start_ring ;
      end_ring = new_end_ring ;

      if ((Gdiag & DIAG_LOGDIFF) && (start_ring <= end_ring))
      {
        float start_rho, end_rho, stime, etime ;
        
        start_rho = lmi->rhos[start_ring] ;
        stime = (float)i * exp(2.0f*start_rho) * KERNEL_MUL ;
        end_rho = lmi->rhos[end_ring] ;
        etime = (float)i * exp(2.0f*end_rho) * KERNEL_MUL ;

#if 0        
        fprintf(stderr, "iteration %d, ring=%d (%2.3f) --> %d (%2.3f) "
                "(time=%2.1f --> %2.1f)\n",
                i, start_ring, start_rho, end_ring, end_rho, stime, etime) ;
#endif
        if (Gdiag & DIAG_WRITE)
        {
          int npixels ;

          npixels = 0 ;
          for (ring = start_ring ; ring <= end_ring ; ring++)
            npixels += lmi->end_spoke[ring] - lmi->start_spoke[ring] + 1 ;

          fprintf(fp, "%d %d %d %d %2.3f %2.3f\n", i, npixels, 
                  start_ring, end_ring, stime, etime) ;
        }
      }
    }
    else
    {
      start_ring = 0 ;
      end_ring = lmi->nrings -1 ;
    }
    
    if (start_ring > end_ring)
      break ;

    LogMapGradient(lmi, tmpImage, gradImage, doweight, start_ring, end_ring) ;

    for (ring = start_ring ; ring <= end_ring ; ring++)
    {
      start_spoke = lmi->start_spoke[ring] ;
      end_spoke = lmi->end_spoke[ring] ;
      for (spoke = start_spoke ; spoke <= end_spoke ; spoke++)
      {
        if (LOG_PIX_AREA(lmi, ring, spoke) <= 0)
          continue ;
        pix = LOG_PIX(lmi, ring, spoke) ;
        if (doweight)
          weight = pix->weight ;
        else
          weight = 1.0 ;

        /* pnpix is a pointer to an array of pointers to neighbors */
        pnpix = pix->nbd ;
        for (pc = c, pf = fvals, j = 0 ; j < NBD_SIZE ; j++, pnpix++)
        {
          npix = *pnpix ;
          n_ring = npix->ring ;
          n_spoke = npix->spoke ;

          *pf++ = *IMAGEFpix(tmpImage, n_ring, n_spoke) ;
          *pc++ = KERNEL_MUL * C(*IMAGEFpix(gradImage, n_ring, n_spoke),k);
        }

        pcself = &c[N_SELF] ;
        for (pc = c, c[N_SELF] = 1.0f, ci = 0 ; ci < NBD_SIZE ; ci++, pc++)
          if (ci != N_SELF)
            *pcself -= *pc ;

        
        for (pc = c,pf = fvals,dst_val = 0.0f,ci = 0 ; ci < NBD_SIZE ; ci++)
          dst_val += *pf++ * *pc++ ;
        
        *IMAGEFpix(outImage, ring, spoke) = dst_val ;

      }   /* each spoke */
    }     /* each ring */

    ImageCopy(outImage, tmpImage) ;
    if (Gdiag & DIAG_MOVIE)
        ImageAppend(tmpImage, LOG_MOVIE_FNAME) ;

    /*    k = k + k * slope ;*/
  }

  if (fp)
    fclose(fp) ;

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
double
LogMapDiffuseCurvature(LOGMAP_INFO *lmi,IMAGE *inImage,IMAGE *outImage,
                       double A, int niterations,int doweight,int time_type)
{
/*  LOGPIX    *npix ;*/
  int       ring, spoke, n_ring, n_spoke, i, j, ci ; 
  float     c[NBD_SIZE], fvals[NBD_SIZE], dst_val ;
  FILE      *fp ;
  IMAGE *srcImage, *dstImage, *tmpImagePtr ;
  static  IMAGE *tmpImage = NULL, *gradImage = NULL ;

  if (inImage->pixel_format != PFFLOAT)
  {
    fprintf(stderr, 
            "LogMapDiffuseCurvature: unsupported input image format %d\n",
            inImage->pixel_format) ;
    exit(1) ;
  }

  if (outImage->pixel_format != PFFLOAT)
  {
    fprintf(stderr, 
            "LogMapDiffuseCurvature: unsupported output image format %d\n",
            outImage->pixel_format) ;
    exit(1) ;
  }

#if DEBUG
  if (Gdiag & DIAG_TIMER)
    DebugTimerStart() ;
#endif

  if (!ImageCheckSize(inImage, tmpImage, 0, 0, 0))
  {
    if (tmpImage)
    {
      ImageFree(&tmpImage) ;
      ImageFree(&gradImage) ;
    }

    tmpImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1) ;
    gradImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1) ;
  }
  else
  {
    ImageSetSize(tmpImage, inImage->rows, inImage->cols) ;
    ImageSetSize(gradImage, inImage->rows, inImage->cols) ;
  }

  ImageCopy(inImage, tmpImage) ;

  if (Gdiag & DIAG_WRITE)
    fp = fopen("diffuse.dat", "w") ;
  else
    fp = NULL ;

  srcImage = tmpImage ;
  dstImage = outImage ;

  for (i = 0 ; i < niterations ; i++)
  {
    if (fp)
    {
      fprintf(fp, "iteration #%d\n", i) ;
    }
    LogMapCurvature(lmi, srcImage, gradImage, A, doweight, -1, -1) ;

    for_each_log_pixel(lmi, ring, spoke)  /* do diffusion on each pixel */
    {
      for (j = 0 ; j < NBD_SIZE ; j++)
      {
        n_ring = LOG_PIX_NBD(lmi, ring, spoke, j)->ring ;
        n_spoke = LOG_PIX_NBD(lmi, ring, spoke, j)->spoke ;
        fvals[j] = *IMAGEFpix(srcImage, n_ring, n_spoke) ;
        c[j] = KERNEL_MUL / *IMAGEFpix(gradImage, n_ring, n_spoke) ;
      }

      /* center coefficient is 1 - (sum of the other coefficients) */
      for (c[N_SELF] = 1.0f, ci = 0 ; ci < NBD_SIZE ; ci++)
        if (ci != N_SELF)
          c[N_SELF] -= c[ci] ;

      for (dst_val = 0.0f, ci = 0 ; ci < NBD_SIZE ; ci++)
        dst_val += fvals[ci] * c[ci] ;

      *IMAGEFpix(dstImage, ring, spoke) = dst_val ;

    }
    /* swap them */
    tmpImagePtr = dstImage ;
    dstImage = srcImage ;
    srcImage = tmpImagePtr ;

    if (fp)
      fprintf(fp, "\\n") ;
  }

  if (dstImage != outImage)
    ImageCopy(dstImage, outImage) ;

  if (fp)
    fclose(fp) ;

  return(A) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
mapCalculateParms(LOGMAP_INFO *lmi)
{
  double minlog, maxlog, logrange, maxr, halfrows, halfcols, k ;

  halfrows = (double)lmi->nrows / 2.0 ;
  halfcols = (double)lmi->ncols / 2.0 ;

  if (lmi->alpha == 0.0)    /* calculate alpha */
  {
    lmi->alpha = (double)lmi->nspokes / (PI) ;
    fprintf(stderr, "setting a = %2.3f\n", lmi->alpha) ;
  }

  if (lmi->nspokes == 0)
  {
    lmi->nspokes = nint(lmi->alpha * PI) ;
    fprintf(stderr, "setting nspokes = %d\n", lmi->nspokes) ;
  }

  k = 2 * lmi->nspokes ;
  maxr = MIN(hypot(halfrows, lmi->alpha), 
          hypot(0.5, lmi->alpha+halfcols)) ;
  maxr = MIN(halfrows, halfcols) ;
  lmi->maxr = maxr ;

  minlog = log(lmi->alpha);             /* min is along the x axis */
  maxlog = log(maxr+lmi->alpha) ;       /* max is along the y axis */
  logrange = maxlog - minlog ;
  if (lmi->nrings == 0)
  {
#if 0
    lmi->nrings = nint((logrange * k) / PI) ;
#else
    lmi->nrings = (int)((logrange * k) / PI) ;
#endif
#if 1
    fprintf(stderr,
            "maxr=%2.3f, minlog %2.3f, maxlog %2.3f, logrange %2.3f\n",
            maxr, minlog, maxlog, logrange) ;
    fprintf(stderr, "nspokes=%d, setting nrings = %d (%2.3f)\n", lmi->nspokes,
            lmi->nrings, (logrange * k) / PI) ;
#endif
  }
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#if !NEW_DIFFUSION

#define PERCENTILE  0.8
#define HIST_BINS   8

static double
diffusionCalculateK(LOGMAP_INFO *lmi, IMAGE *image, double percentile)
{
  int    i, ring, spoke, npix, ppix, hist_pix[HIST_BINS], oval, val ;
  LOGPIX *opix ;
  double hist_vals[HIST_BINS], bin_size, k, rval, lval, deltaW, deltaN, 
         fval, min_val, max_val ;
  IMAGE  *Itmp ;

  if (image->pixel_format != PFINT)
  {
    Itmp = ImageAlloc(image->rows, image->cols, PFINT, 1) ;
    ImageCopy(image, Itmp) ;
    image = Itmp ;
  }
  else
    Itmp = NULL ;

  if (percentile == 0.0)
    percentile = PERCENTILE ;

  npix = 0 ;
  memset(hist_pix, 0, HIST_BINS * sizeof(double)) ;
  max_val = -1.0 ;
  min_val = 100000.0 ;

  /* find range and # of log pixels */
  for_each_log_pixel(lmi, ring, spoke)
  {
    val = *IMAGEIpix(image, ring, spoke) ;

    opix = LOG_PIX_NBD(lmi, ring, spoke, N_N) ;
    oval = *IMAGEIpix(image, opix->ring, opix->spoke) ;
    deltaN = (double)abs(oval - val) ;

    opix = LOG_PIX_NBD(lmi, ring, spoke, N_W) ;
    oval = *IMAGEIpix(image, opix->ring, opix->spoke) ;
    deltaW = (double)abs(oval - val) ;

    fval = (deltaN + deltaW) / 2.0 ;
    if (fval > max_val) 
      max_val = fval ;
    if (fval < min_val)
      min_val = fval ;
    npix++ ;
  }

#if 0
  fprintf(stderr, "range %2.3f --> %2.3f\n", min_val, max_val) ;
#endif

  bin_size = (max_val - min_val) / (double)HIST_BINS ;
  for (i = 0 ; i < HIST_BINS ; i++)
    hist_vals[i] = (i+1)*bin_size ; /* value at right end of histogram range */

    
  ppix = (int)((double)npix * percentile) ;

  /* build histogram */
  for_each_log_pixel(lmi, ring, spoke)
  {
    val = *IMAGEIpix(image, ring, spoke) ;

    opix = LOG_PIX_NBD(lmi, ring, spoke, N_N) ;
    oval = *IMAGEIpix(image, opix->ring, opix->spoke) ;
    deltaN = (double)abs(oval - val) ;

    opix = LOG_PIX_NBD(lmi, ring, spoke, N_W) ;
    oval = *IMAGEIpix(image, opix->ring, opix->spoke) ;
    deltaW = (double)abs(oval - val) ;

    fval = (deltaN + deltaW) / 2 ;
    for (i = 0 ; i < HIST_BINS ; i++)
    {
      if (fval < hist_vals[i])
      {
        hist_pix[i]++ ;
        break ;
      }
    }
  }

#if 0
{
  FILE *fp ;
  fp = fopen("histo.dat", "w") ;

  for (npix = 0, i = 0 ; i < HIST_BINS ; i++)
  {
    fprintf(fp, "%d  %d  %f\n", i, hist_pix[i], hist_vals[i]) ;
  }
  fclose(fp) ;
}
#endif

  for (npix = 0, i = 0 ; i < HIST_BINS ; i++)
  {
    npix += hist_pix[i] ;
    if (npix >= ppix)  /* set k to middle of the range of this bin */
    {
      rval = hist_vals[i] ;
      if (!i)
        lval = 0.0 ;
      else
        lval = hist_vals[i-1] ;
      k = (rval + lval) / 2.0 ;
#if 0
      fprintf(stderr,"using bin %d: setting k=(%2.3f + %2.3f)/2 = %2.3f\n", 
              i, rval, lval, k) ;
#endif
      break ;
    }
  }

  if (Itmp)
    ImageFree(&Itmp) ;

  return(k) ;
}
#endif

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
LogMapPatchHoles(LOGMAP_INFO *lmi, IMAGE *Itv, IMAGE *Ilog)
{
  float  fval ;
  int    crow, ccol, ring, spoke ;
  IMAGE  *Iout ;

  if (Ilog->pixel_format != PFFLOAT)
  {
    Iout = ImageAlloc(Ilog->rows, Ilog->cols, PFFLOAT, 1) ;
    ImageCopy(Ilog, Iout) ;
  }
  else
    Iout = Ilog ;

  for_each_log_pixel(lmi, ring, spoke)
  {
    crow = LOG_PIX_ROW_CENT(lmi, ring, spoke) ;
    ccol = LOG_PIX_COL_CENT(lmi, ring, spoke) ;
    if (crow != UNDEFINED && ccol != UNDEFINED)
    {
      fval = *IMAGEFpix(Iout, ring, spoke) ;
      if (fval != (float)UNDEFINED && !FZERO(fval))
        continue ;

      switch (Itv->pixel_format)
      {
      case PFFLOAT:
        fval = *IMAGEFpix(Itv, ccol, crow) ;
        break ;
      case PFBYTE:
        fval = (float)*IMAGEpix(Itv, ccol, crow) ;
        break ;
      default:
        ErrorSet(ERROR_UNSUPPORTED, 
                         "LogMapPatchHoles: unsupported pixel format %d",
                         Itv->pixel_format) ;
        return ;
      }
      *IMAGEFpix(Iout, ring, spoke) = fval ;
    }
  }

  if (Ilog != Iout)
  {
    ImageCopy(Iout, Ilog) ;
    ImageFree(&Iout) ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
logPatchHoles(LOGMAP_INFO *lmi)
{
  int       nvalid, ring, spoke, i ;
  float     xtotal, ytotal ;

  for (i = 0 ; i < 2 ; i++)
  {
    for (spoke = 0 ; spoke < lmi->nspokes ; spoke++)
    {
      for (ring = 0 ; ring < lmi->nrings ; ring++)
      {
        if (LOG_PIX_AREA(lmi, ring, spoke))
          continue ;

        nvalid = 0 ;
        xtotal = 0.0f ;
        ytotal = 0.0f ;
        if (ring < lmi->nrings-1)
        {
          if (LOG_PIX_AREA(lmi, ring+1, spoke))
          {
            xtotal += LOG_PIX_COL_CENT(lmi, ring+1, spoke) ;
            ytotal += LOG_PIX_ROW_CENT(lmi, ring+1, spoke) ;
            nvalid++ ;
          }
          if (spoke < lmi->nspokes-1)
          {
            if (LOG_PIX_AREA(lmi, ring+1, spoke+1))
            {
              xtotal += LOG_PIX_COL_CENT(lmi, ring+1, spoke+1) ;
              ytotal += LOG_PIX_ROW_CENT(lmi, ring+1, spoke+1) ;
              nvalid++ ;
            }
          }
          if (spoke > 0)
          {
            if (LOG_PIX_AREA(lmi, ring+1, spoke-1))
            {
              xtotal += LOG_PIX_COL_CENT(lmi, ring+1, spoke-1) ;
              ytotal += LOG_PIX_ROW_CENT(lmi, ring+1, spoke-1) ;
              nvalid++ ;
            }
          }
        }
        if (ring > 0)
        {
          if (LOG_PIX_AREA(lmi, ring-1, spoke))
          {
            xtotal += LOG_PIX_COL_CENT(lmi, ring-1, spoke) ;
            ytotal += LOG_PIX_ROW_CENT(lmi, ring-1, spoke) ;
            nvalid++ ;
          }
          if (spoke < lmi->nspokes-1)
          {
            if (LOG_PIX_AREA(lmi, ring-1, spoke+1))
            {
              xtotal += LOG_PIX_COL_CENT(lmi, ring-1, spoke+1) ;
              ytotal += LOG_PIX_ROW_CENT(lmi, ring-1, spoke+1) ;
              nvalid++ ;
            }
          }
          if (spoke > 0)
          {
            if (LOG_PIX_AREA(lmi, ring-1, spoke-1))
            {
              xtotal += LOG_PIX_COL_CENT(lmi, ring-1, spoke-1) ;
              ytotal += LOG_PIX_ROW_CENT(lmi, ring-1, spoke-1) ;
              nvalid++ ;
            }
          }
        }
        if (spoke < lmi->nspokes-1)
        {
          if (LOG_PIX_AREA(lmi, ring, spoke+1))
          {
            xtotal += LOG_PIX_COL_CENT(lmi, ring, spoke+1) ;
            ytotal += LOG_PIX_ROW_CENT(lmi, ring, spoke+1) ;
            nvalid++ ;
          }
        }
        if (spoke > 0)
        {
          if (LOG_PIX_AREA(lmi, ring, spoke-1))
          {
            xtotal += LOG_PIX_COL_CENT(lmi, ring, spoke-1) ;
            ytotal += LOG_PIX_ROW_CENT(lmi, ring, spoke-1) ;
            nvalid++ ;
          }
        }
        if (nvalid > 4)
        {
          LOG_PIX_ROW_CENT(lmi, ring, spoke) = nint(ytotal / (float)nvalid);
          LOG_PIX_COL_CENT(lmi, ring, spoke) = nint(xtotal / (float)nvalid);
#if 0
          fprintf(stderr, "patched hole at (%d, %d) -> (%d, %d)\n",
                  ring, spoke, LOG_PIX_COL_CENT(lmi,ring,spoke),
                  LOG_PIX_ROW_CENT(lmi,ring,spoke)) ;
#endif
        }
      }
    }
    /* reset area now. Couldn't do it earlier otherwise would have cascade*/
    for_all_log_pixels(lmi, ring, spoke)
    {
      int crow, ccol ;

      crow = LOG_PIX_ROW_CENT(lmi,ring,spoke) ;
      ccol = LOG_PIX_COL_CENT(lmi, ring, spoke);
      if (!LOG_PIX_AREA(lmi, ring, spoke) && (crow != 0) && (ccol != 0))
      {
        float x, y, rho, phi, weight ;
        int   halfrows, halfcols ;

        LOG_PIX_SPOKE(lmi, ring, spoke) = spoke ;
        LOG_PIX_RING(lmi, ring, spoke) = ring ;
        LOG_PIX_AREA(lmi, ring, spoke) = 1 ;
        TV_TO_RING(lmi, ccol, crow) = ring ;
        TV_TO_SPOKE(lmi, ccol, crow) = spoke ;
        halfrows = lmi->nrows / 2 ;
        halfcols = lmi->ncols / 2 ;
        
        x = (float)(LOG_PIX_COL_CENT(lmi, ring, spoke) - halfcols) ;
        y = (float)(LOG_PIX_ROW_CENT(lmi, ring, spoke) - halfrows);
        if (x < 0)
          x = -x ;
        
        /* not quite right, but good enough for weight */
        LOG_PIX_RHO(lmi,ring,spoke) = rho =
          log(sqrt(SQR(x+lmi->alpha)+SQR(y)));
        LOG_PIX_PHI(lmi,ring,spoke) = phi = atan2(y, x+lmi->alpha) ;
        weight = exp(-rho) ;
        LOG_PIX_WEIGHT(lmi, ring, spoke) = weight ;
      }
    }
  }
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
LogMapNormalize(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst, float low,float hi)
{
  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 1) ;

  switch (Isrc->pixel_format)
  {
  case PFFLOAT:
    fnormalize(lmi, Isrc, Idst, low, hi) ;
    break ;
  case PFBYTE:
    bnormalize(lmi, Isrc, Idst, (byte)low, (byte)hi) ;
    break ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
fnormalize(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst, float low, float hi)
{
  int     ring, spoke ;
  float  val, min_val, max_val, scale ;

  max_val = 0.0f ;
  min_val = 100000.0f ;
  for_each_log_pixel(lmi, ring, spoke)
  {
    val = *IMAGEFpix(Isrc, ring, spoke) ;
    if (val < min_val)
      min_val = val ;
    if (val > max_val)
      max_val = val ;
  }
  scale = (hi - low) / (max_val - min_val) ;
  for_each_log_pixel(lmi, ring, spoke)
  {
    val = *IMAGEFpix(Isrc, ring, spoke) ;
    val = (val - min_val) * scale + low ;
    *IMAGEFpix(Idst, ring, spoke) = val ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
bnormalize(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst, byte low, byte hi)
{
  int     ring, spoke ;
  byte    val, min_val, max_val ;
  float   scale ;

  max_val = 0 ;
  min_val = 255 ;
  for_each_log_pixel(lmi, ring, spoke)
  {
    val = *IMAGEpix(Isrc, ring, spoke) ;
    if (val < min_val)
      min_val = val ;
    if (val > max_val)
      max_val = val ;
  }
  scale = (float)(hi - low) / (float)(max_val - min_val) ;
  for_each_log_pixel(lmi, ring, spoke)
  {
    val = *IMAGEpix(Isrc, ring, spoke) ;
    val = (byte)((float)(val - min_val) * scale) + low ;
    *IMAGEpix(Idst, ring, spoke) = val ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
logBuildRhoList(LOGMAP_INFO *lmi)
{
  int  ring, spoke, nspokes, nrings, nrho ;
  float total, rho, min_rho, max_rho ;

  nrings = lmi->nrings ;
  nspokes = lmi->nspokes ;

  lmi->rhos = (float *)calloc(nrings, sizeof(float)) ;
  for (ring = 0; ring < nrings; ring++) 
  {
    nrho = 0 ;
    total = 0.0f ;
    min_rho = 1000.0f ;
    max_rho = 0.0f ;
    for (spoke = 0 ; spoke < nspokes ; spoke++) 
    {
      if (LOG_PIX_AREA(lmi, ring, spoke) > 0)
      {
        nrho++ ;
        rho = LOG_PIX_RHO(lmi,ring,spoke) ;
        total += rho ;
        if (rho >= max_rho)
          max_rho = rho ;
        if (!FZERO(rho) && (rho <= min_rho))
          min_rho = rho ;
      }
    }
    if (!nrho)
      lmi->rhos[ring] = 0.0f ;
    else
      lmi->rhos[ring] = total / (float)nrho ;

#if 0   
    fprintf(stdout, "ring %d: rho %2.3f --> %2.3f, avg %2.3f, N %d\n",
            ring, min_rho, max_rho,  lmi->rhos[ring], nrho) ;
#endif
  }

#if 0
  /* prevent empty last ring */
  if (FZERO(lmi->rhos[nrings-1]))
    nrings = lmi->nrings = nrings - 1 ;
#endif

  min_rho = 1000.0f ; 
  max_rho = 0.0f ;
  for (ring = 0; ring < nrings; ring++) 
  {
    rho = lmi->rhos[ring] ;
    if (rho >= max_rho)
      max_rho = rho ;
    if (!FZERO(rho) && (rho <= min_rho))
      min_rho = rho ;
  }
  lmi->min_rho = min_rho ;
  lmi->max_rho = max_rho ;
/*  fprintf(stderr, "rho: %2.3f --> %2.3f\n", min_rho, max_rho) ;*/
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              build a list of starting and ending spokes with valid
              data.
----------------------------------------------------------------------*/
static void
logBuildValidSpokeList(LOGMAP_INFO *lmi)
{
  int  ring, spoke, nspokes, nrings, start_spoke, end_spoke ;

  nrings = lmi->nrings ;
  nspokes = lmi->nspokes ;

  lmi->start_spoke = (int *)calloc(nrings, sizeof(int)) ;
  lmi->end_spoke = (int *)calloc(nrings, sizeof(int)) ;
  for (ring = 0; ring < nrings; ring++) 
  {
    start_spoke = 0 ;
    end_spoke = -1 ;
    for (spoke = 0 ; spoke < nspokes ; spoke++) 
    {
      if (LOG_PIX_AREA(lmi, ring, spoke) > 0)
      {
        start_spoke = spoke ;
        break ;
      }
    }
    for (spoke = nspokes-1 ; spoke >= 0 ; spoke--) 
    {
      if (LOG_PIX_AREA(lmi, ring, spoke) > 0)
      {
        end_spoke = spoke ;
        break ;
      }
    }
    lmi->start_spoke[ring] = start_spoke ;
    lmi->end_spoke[ring] = end_spoke ;
/*    fprintf(stdout, "%d: %d --> %d\n", ring, start_spoke, end_spoke) ;*/
  }
}
static void
writeIterations(char *fname, LOGMAP_INFO *lmi, float end_time)
{
  int    ring, spoke ;
  float  rho, niter ;
  FILE   *fp ;

  if (!strcmp(fname, "stdout"))
    fp = stdout ;
  else
    fp = fopen(fname, "w") ;

  for (spoke = 0 ; spoke < lmi->nspokes ; spoke++)
  {
    for (ring = 0 ; ring < lmi->nrings ; ring++)
    {
      if (LOG_PIX_AREA(lmi, ring, spoke) <= 0)
      {
        niter = 0.0f ;
      }
      else
      {
        rho = LOG_PIX_RHO(lmi, ring, spoke) ;
        niter = end_time / (exp(2.0f * rho) * KERNEL_MUL) ;
      }
      fprintf(fp, "%3.3f ", niter) ;
    }
    fprintf(fp, "\n") ;
  }

  if (fp != stdout)
    fclose(fp) ;
}
static void
writeTimes(char *fname, LOGMAP_INFO *lmi, int niter)
{
  int    ring, spoke ;
  float  rho, end_time ;
  FILE   *fp ;

  if (!strcmp(fname, "stdout"))
    fp = stdout ;
  else
    fp = fopen(fname, "w") ;

  for (spoke = 0 ; spoke < lmi->nspokes ; spoke++)
  {
    for (ring = 0 ; ring < lmi->nrings ; ring++)
    {
      if (LOG_PIX_AREA(lmi, ring, spoke) <= 0)
      {
        end_time = 0.0f ;
      }
      else
      {
        rho = LOG_PIX_RHO(lmi, ring, spoke) ;
        end_time = niter * (exp(2.0f * rho) * KERNEL_MUL) ;
      }
      fprintf(fp, "%3.3f ", end_time) ;
    }
    fprintf(fp, "\n") ;
  }

  if (fp != stdout)
    fclose(fp) ;
}

#define FSCALE  100000.0f

#define ISZERO(f)   FZERO(f*10000000.0f)

static int imageOffsetDirection(IMAGE *Ix, IMAGE *Iy, int wsize, 
                                IMAGE *Iorient, int x0,int y0) ;
static int
imageOffsetDirection(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Iorient, 
                     int x0,int y0)
{
  int    rows, cols, x, y, whalf, xc, yc, yoff, off, d ;
  float  *xpix, *ypix, dx, dy, *or_xpix,*or_ypix, dir, ox, oy ;

  rows = Ix->rows ;
  cols = Ix->cols ;

  whalf = (wsize-1)/2 ;
  xpix = IMAGEFpix(Ix, x0, y0) ;
  ypix = IMAGEFpix(Iy, x0, y0) ;
  or_xpix = IMAGEFpix(Iorient, x0, y0) ;
  or_ypix = IMAGEFseq_pix(Iorient, x0, y0, 1) ;

/*
  Now calculate the orientation for this point by averaging local gradient
  orientation within the specified window.

  x and y are in window coordinates, while xc and yc are in image
  coordinates.
*/
  /* calculate orientation vector */
  ox = *or_xpix++ ;
  oy = *or_ypix++ ;
  
  dir = 0.0f ;

  for (y = -whalf ; y <= whalf ; y++)
  {
    /* reflect across the boundary */
    yc = y + y0 ;
    if ((yc < 0) || (yc >= rows))
      continue ;
    
    yoff = y*cols ;
    for (x = -whalf ; x <= whalf ; x++)
    {
      xc = x0 + x ;
      if ((xc < 0) || (xc >= cols))
        continue ;
      
      off = yoff + x ;
      dx = *(xpix+off) ;
      dy = *(ypix+off) ;
      dir += (x*ox + y*oy) * fabs(dx*ox + dy*oy) ;
    }
  }
  
  if (ISZERO(dir))
    d = 0 ;
  else if (dir > 0.0f)   /* flip by 180 */
    d = -1 ;
  else
    d = 1 ;

  return(d) ;
}

#define BR     7
#define BS     22
#define BAD(r,s)  ((r==BR) && (s==BS))

#define WINDOW_SIZE   3
#define WHALF         ((WINDOW_SIZE-1)/2)
#define WSQ           (WINDOW_SIZE*WINDOW_SIZE)
#define MEDIAN_INDEX  (WSQ/2)


#define USE_MEDIAN 0
#if USE_MEDIAN
static int compare_sort_farray(const void *pf1, const void *pf2) ;
static int compare_sort_barray(const void *pb1, const void *pb2) ;
#endif


 IMAGE *Icalculated = NULL, *ImapOffset = NULL, *ImapOrient = NULL,
         *Ix = NULL, *Iy = NULL;
IMAGE *
LogMapNonlocal(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Ismooth, IMAGE *Idst)
{
  int          msec ;
  struct timeb then ;
  int   rows, cols, x, y, ax, ay, sx, sy, x1, y1, dx, dy, odx, ody, d, xn, yn,
        steps, dir, dot, wsize = 3, maxsteps = 6 ;
  float *src_xpix, *src_ypix, *oxpix, *oypix, fdir ;
  byte  *calculated ;
  int   ring, spoke ;
  IMAGE *Iout ;
#if USE_MEDIAN
  int    xc, yc ;
  float  sort_farray[WINDOW_SIZE*WINDOW_SIZE], *fptr, *fpix ;
  byte   sort_barray[WINDOW_SIZE*WINDOW_SIZE], *bptr, *bpix ;
#endif

  rows = Isrc->rows ;
  cols = Isrc->cols ;

  if (!Icalculated || (rows != Icalculated->rows || cols != Icalculated->cols))
  {
    if (Icalculated)
      ImageFree(&Icalculated) ;
    if (ImapOffset)
      ImageFree(&ImapOffset) ;
    if (ImapOrient)
      ImageFree(&ImapOrient) ;
    if (Ix)
      ImageFree(&Ix) ;
    if (Iy)
      ImageFree(&Iy) ;
    Icalculated = ImageAlloc(rows, cols, PFBYTE, 1) ;
    Ix = ImageAlloc(rows, cols, PFFLOAT, 1) ;
    Iy = ImageAlloc(rows, cols, PFFLOAT, 1) ;
    ImapOffset = ImageAlloc(rows, cols, PFFLOAT, 2) ;
    ImapOrient = ImageAlloc(rows, cols, PFFLOAT, 2) ;
  }
  else
    ImageClear(Icalculated) ;

  if (Gdiag & DIAG_TIMER)
    TimerStart(&then) ;
  ImageSobel(Ismooth, NULL, Ix, Iy) ;
  if (Gdiag & DIAG_TIMER)
  {
    msec = TimerStop(&then) ;
    fprintf(stderr, "sobel took              %d msec (%2.1f Hz)\n",
            msec, 1000.0f / ((float)msec)) ;
    TimerStart(&then) ;
  }
#if 0
  ImapOrient = ImageOffsetOrientation(Ix, Iy, 3, ImapOrient) ;
#else
  ImageCopyFrames(Ix, ImapOrient, 0, 1, 0) ;
  ImageCopyFrames(Iy, ImapOrient, 0, 1, 1) ;
#endif
  if (Gdiag & DIAG_TIMER)
  {
    msec = TimerStop(&then) ;
    fprintf(stderr, "orientation took        %d msec (%2.1f Hz)\n",
            msec, 1000.0f / ((float)msec)) ;
  }

  if (!Idst)
    Idst = ImageAlloc(lmi->nspokes,lmi->nrings, Isrc->pixel_format,1);
  else
    ImageClearArea(Idst, 0, 0, Idst->rows, Idst->cols, 0.0f) ;

  if (Idst->pixel_format != Isrc->pixel_format)
    Iout = ImageAlloc(lmi->nspokes, lmi->nrings, Isrc->pixel_format,1);
  else
    Iout = Idst ;

  if (Gdiag & DIAG_TIMER)
    TimerStart(&then) ;
  for_each_log_pixel(lmi, ring, spoke)
  {
    y1 = LOG_PIX_ROW_CENT(lmi, ring, spoke) ;
    x1 = LOG_PIX_COL_CENT(lmi, ring, spoke) ;
    if (y1 != UNDEFINED && x1 != UNDEFINED)
    {
      calculated = IMAGEpix(Icalculated, x1, y1) ;
      src_xpix = IMAGEFpix(ImapOrient, x1, y1) ;
      src_ypix = IMAGEFseq_pix(ImapOrient, x1, y1, 1) ;

      /* do a Bresenham algorithm do find the offset line at this point */
      if (*calculated == 0)
      {
        dir = imageOffsetDirection(Ix, Iy, wsize, ImapOrient, x1, y1) ;
        fdir = (float)dir ;
        *calculated = 1 ;
        *IMAGEFpix(ImapOffset, x1, y1) = *src_xpix * fdir ;
        *IMAGEFseq_pix(ImapOffset, x1, y1, 1) = *src_ypix * fdir ;
      }
      dx = nint(*IMAGEFpix(ImapOffset,x1,y1) * FSCALE) ;
      dy = nint(*IMAGEFseq_pix(ImapOffset,x1,y1,1) * FSCALE) ;
      x = x1 ;
      y = y1 ;
      ax = ABS(dx) << 1 ;
      sx = SGN(dx) ;
      ay = ABS(dy) << 1 ;
      sy = SGN(dy) ;
      
      oxpix = src_xpix++ ;
      oypix = src_ypix++ ;
      
      if (ax > ay)  /* x dominant */
      {
        d = ay - (ax >> 1) ;
        for (steps = 0 ; steps < maxsteps ; steps++)
        {
          if (!*IMAGEpix(Icalculated, x, y))
          {
            dir = imageOffsetDirection(Ix, Iy, wsize, ImapOrient, x, y) ;
            fdir = (float)dir ;
            *IMAGEpix(Icalculated, x, y) = 1 ;
            *IMAGEFpix(ImapOffset, x, y) = *oxpix * fdir ;
            *IMAGEFseq_pix(ImapOffset, x, y, 1) = *oypix * fdir ;
          }
          odx = nint(*IMAGEFpix(ImapOffset,x,y) * FSCALE) ;
          ody = nint(*IMAGEFseq_pix(ImapOffset,x,y,1) * FSCALE) ;
          dot = odx * dx + ody * dy ;
          if (dot <= 0)
            break ;
          if (d >= 0)
          {
            yn = y + sy ;
            if (yn < 0 || yn >= rows)
              break ;
            oxpix += (sy * cols) ;
            oypix += (sy * cols) ;
            d -= ax ;
          }
          else
            yn = y ;
          oxpix += sx ;
          oypix += sx ;
          xn = x + sx ;
          if (xn < 0 || xn >= cols)
            break ;

          x = xn ;
          y = yn ;
          d += ay ;
        }
      }
      else    /* y dominant */
      {
        d = ax - (ay >> 1) ;
        for (steps = 0 ; steps < maxsteps ; steps++)
        {
          if (!*IMAGEpix(Icalculated, x, y))
          {
            dir = imageOffsetDirection(Ix, Iy, wsize, ImapOrient, x, y) ;
            fdir = (float)dir ;
            *IMAGEpix(Icalculated, x, y) = 1 ;
            *IMAGEFpix(ImapOffset, x, y) = *oxpix * fdir ;
            *IMAGEFseq_pix(ImapOffset, x, y, 1) = *oypix * fdir ;
          }
          odx = nint(*IMAGEFpix(ImapOffset,x,y) * FSCALE) ;
          ody = nint(*IMAGEFseq_pix(ImapOffset,x,y,1) * FSCALE) ;
          dot = odx * dx + ody * dy ;
          if (dot <= 0)
            break ;
          if (d >= 0)
          {
            xn = x + sx ;
            if (xn < 0 || xn >= cols)
              break ;
            oxpix += sx ;
            oypix += sx ;
            d -= ay ;
          }
          else
            xn = x ;
          yn = y + sy ;
          if (yn < 0 || yn >= rows)
            break ;

          x = xn ;
          y = yn ;
          oypix += (sy * cols) ;
          oxpix += (sy * cols) ;
          d += ax ;
        }
      }

      dx = (x - x1) ;
      dy = (y - y1) ;
      *IMAGEFpix(ImapOffset, x1, y1) = (float)dx ;        /* for diagnostics */
      *IMAGEFseq_pix(ImapOffset, x1, y1, 1) = (float)dy ; /* for diagnostics */

      switch (Isrc->pixel_format)
      {
      case PFFLOAT:
#if USE_MEDIAN
        /* build array of neighboring pixels - x and y are reused */
        for (fptr = sort_farray, y = -WHALF ; y <= WHALF ; y++)
        {
          /* reflect across the boundary */
          yc = y + y1 + dy ;

          /* don't need border checking as centroid wont be on border */
          
          fpix = IMAGEFpix(Isrc, 0, yc) ;
          for (x = -WHALF ; x <= WHALF ; x++)
          {
            xc = x1 + x + dx ;
            /* don't need border checking as centroid wont be on border */
            
            *fptr++ = *(fpix + xc) ;
          }
        }
        qsort(sort_farray, WSQ, sizeof(float), compare_sort_farray) ;
        *IMAGEFpix(Iout, ring, spoke) = sort_farray[MEDIAN_INDEX] ;
#else
        *IMAGEFpix(Iout, ring, spoke) = *IMAGEFpix(Ismooth,x1+dx, y1+dy) ;
#endif
        break ;
      case PFBYTE:
        /* build array of neighboring pixels - x and y are reused */
#if USE_MEDIAN
        for (bptr = sort_barray, y = -WHALF ; y <= WHALF ; y++)
        {
          /* reflect across the boundary */
          yc = y + y1 + dy ;
          if (yc < 0)
            yc = 0 ;
          else if (yc >= rows)
            yc = rows - 1 ;
          
          bpix = IMAGEpix(Isrc, 0, yc) ;
          for (x = -WHALF ; x <= WHALF ; x++)
          {
            xc = x1 + x + dx ;
            if (xc < 0)
              xc = 0 ;
            else if (xc >= cols)
              xc = cols - 1 ;
            
            *bptr++ = *(bpix + xc) ;
          }
        }
        qsort(sort_barray, WSQ, sizeof(float), compare_sort_barray) ;
        *IMAGEpix(Iout, ring, spoke) = sort_barray[MEDIAN_INDEX] ; ;
#else
        *IMAGEpix(Iout, ring, spoke) = *IMAGEpix(Ismooth,x1+dx, y1+dy) ;
#endif
        break ;
      default:
        ErrorReturn(Idst, (ERROR_UNSUPPORTED, 
                         "LogMapNonlocal: unsupported pixel format %d",
                         Iout->pixel_format)) ;
        break ;
      }
    }
  }

  if (Gdiag & DIAG_TIMER)
  {
    msec = TimerStop(&then) ;
    fprintf(stderr, "nonlocal filtering took %d msec (%2.1f Hz)\n",
            msec, 1000.0f / ((float)msec)) ;
  }
#if 0
ImageWrite(Icalculated, "calc.hipl") ;
ImageWrite(ImapOffset, "offset_bad.hipl") ;
#endif

  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }
  return(Idst) ;
}

#if USE_MEDIAN
static int
compare_sort_farray(const void *pf1, const void *pf2)
{
  register float f1, f2 ;

  f1 = *(float *)pf1 ;
  f2 = *(float *)pf2 ;

/*  return(f1 > f2 ? 1 : f1 == f2 ? 0 : -1) ;*/
  if (f1 > f2)
    return(1) ;
  else if (f1 < f2)
    return(-1) ;

  return(0) ;
}
static int
compare_sort_barray(const void *pb1, const void *pb2)
{
  register byte b1,b2 ;

  b1 = *(float *)pb1 ;
  b2 = *(float *)pb2 ;

/*  return(b1 > b2 ? 1 : b1 == b2 ? 0 : -1) ;*/
  if (b1 > b2)
    return(1) ;
  else if (b1 < b2)
    return(-1) ;

  return(0) ;
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
LogMapInitForwardFilter(LOGMAP_INFO *lmi, int which)
{
  int                 rows, cols, row, col, npix, ring, spoke, nrings,
                      nspokes, ccol, crow ;
  CARTESIAN_MAP_INFO  *cmi ;
  CARTESIAN_PIXEL     *cp ;
  LOGPIX              *lpix ;
  float               radius, dist ;
  int                 overflows = 0, col_end, row_end, ring1, ring2 ;

  rows = lmi->nrows ;
  cols = lmi->ncols ;  /* size of cartesian image */
  npix = rows * cols ;
  cmi = &lmi->cmi ;
  cp = cmi->pix = (CARTESIAN_PIXEL *)calloc(npix, sizeof(CP)) ;
  if (!cmi->pix)
    ErrorReturn(ERROR_NO_MEMORY, (ERROR_NO_MEMORY,
                "LogMapInitForwardFilter: could not allocate LUT")) ;

#if 0
  fprintf(stderr, "building lookup tables...\n") ;
#endif
  nrings = lmi->nrings ;
  nspokes = lmi->nspokes ;
  for (ring = 0 ; ring < nrings ; ring++)
  {
    if (ring < nrings-1)
    {
      ring1 = ring ;
      ring2 = ring+1 ;
      if (FZERO(lmi->rhos[ring2])) /* necessary because of empty rings */
      {
        ring1-- ;
        ring2-- ;
      }
    }
    else
    {
      ring1 = ring-1 ;
      ring2 = ring ;
      if (FZERO(lmi->rhos[ring1]))
      {
        ring1++ ;
        ring2++ ;
      }
    }
    radius = exp(lmi->rhos[ring2]) - exp(lmi->rhos[ring1]) ;
    radius = fabs(radius) ;

#if 0
    fprintf(stderr, "\rring, %d, radius %2.1f ", ring, radius) ;
#endif
    for (spoke = 0 ; spoke < nspokes ; spoke++)
    {
      lpix = LOG_PIX(lmi, ring, spoke) ;
      if (lpix->area <= 0)
        continue ;

      ccol = lpix->col_cent ;
      crow = lpix->row_cent ;
      row_end = MIN(crow+radius, rows-1) ;
      col_end = MIN(ccol+radius, cols-1) ;
      for (row = MAX(0,crow-radius) ; row <= row_end ; row++)
      {
        for (col = MAX(ccol-radius,0) ; col <= col_end ; col++)
        {
          /* check to see if this cart. pix is within sampling distance */
          dist = sqrt((float)((ccol-col)*(ccol-col)+(crow-row)*(crow-row))) ;
          if (dist <= radius)
          {
            cp = cmi->pix + row * cols + col ;
#if 0
            if (col == 128 && row == 231)
              fprintf(stdout, "%d: (%d, %d) <-- (%d, %d) (%d, %d) (%2.2f)\n",
                      cp->npix, col, row, ccol, crow, ring, spoke, dist) ;
#endif
            if (cp->npix >= MAX_PIX)   /* no more room */
            {
              overflows++ ;
#if 0
              fprintf(stdout, "%d: (%d, %d) <-- (%d, %d) (%d, %d) (%2.2f)\n",
                      cp->npix, col, row, ccol, crow, ring, spoke, dist) ;
#endif
              continue ;
            }
            cp->logpix[cp->npix++] = lpix ;
            lpix->ncpix++ ;
          }
        }
      }
    }
  }

#if 0
  col = 128 ;
  row = 231 ;
  cp = cmi->pix + row * cols + col ;
  fprintf(stdout, "(%d, %d) <-- %d\n", col, row, cp->npix) ;
  fflush(stdout) ;

{
  FILE *fp ;

  fp = fopen("lut.dat", "wb") ;
  cp = cmi->pix ;
  for (row = 0 ; row < rows ;  row++)
  {
    for (col = 0 ; col < cols ; col++, cp++)
      fprintf(fp, "%d ", cp->npix) ;
    fprintf(fp, "\n") ;
  }
  fclose(fp) ;
}
#endif

  fprintf(stderr, "LUT built with %d overflows\n", overflows) ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
LogMapForwardFilter(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Idst)
{
  int    rows, cols, row, col, ring, spoke, npix, i ;
  IMAGE  *Iout ;
  LOGPIX *lpix ;
  CP     *cp ;
  CMI    *cmi ;
  float  *spix ;

  if (!Idst)
    Idst = ImageAlloc(lmi->nspokes,lmi->nrings, Isrc->pixel_format,1);
  else
    ImageClearArea(Idst, 0, 0, Idst->rows, Idst->cols, 0.0f) ;

  if (Idst->pixel_format != Isrc->pixel_format)
    Iout = ImageAlloc(lmi->nspokes, lmi->nrings, Isrc->pixel_format,1);
  else
    Iout = Idst ;

  rows = lmi->nrows ;
  cols = lmi->ncols ;
  cmi = &lmi->cmi ;

  cp = cmi->pix ;
  spix = IMAGEFpix(Isrc, 0, 0) ;
  for (row = 0 ; row < rows ; row++)
  {
    for (col = 0 ; col < cols ; col++, cp++, spix++)
    {
      npix = cp->npix ;
      for (i = 0 ; i < npix ; i++)
      {
        lpix = cp->logpix[i] ;
        ring = lpix->ring ;
        spoke = lpix->spoke ;
        *IMAGEFpix(Iout, ring, spoke) += *spix ;
#if 0
        if (ring == 1 && spoke == 0)
          fprintf(stderr, "(1,0) += (%d, %d) %2.3f = %2.3f\n",
                  col, row, *spix, *IMAGEFpix(Iout, ring, spoke)) ;
#endif
      }
    }
  }
  
  /* now normalize logmap by area of each pixel */
  for_all_log_pixels(lmi, ring, spoke)
  {
    lpix = LOG_PIX(lmi, ring, spoke) ;
#if 0
    if (ring == 1 && spoke == 0)
      fprintf(stderr, "normalizing (1,0) = %2.3f by %d\n", 
              *IMAGEFpix(Iout, ring, spoke), lpix->ncpix) ;
#endif
    if (lpix->ncpix > 0)
      *IMAGEFpix(Iout, ring, spoke) /= (float)lpix->ncpix ;
  }

  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }
  return(Idst) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
LogMapFilter(LOGMAP_INFO *lmi, int which, int wsize, IMAGE *Isrc, IMAGE *Idst)
{
  IMAGE  *Iout ;

  if (wsize != 3)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "LogMapFilter: unsupported window size %d", wsize)) ;

  if (!Idst)
    Idst = ImageAlloc(lmi->nspokes,lmi->nrings, Isrc->pixel_format,1);
  else
    ImageClearArea(Idst, 0, 0, Idst->rows, Idst->cols, 0.0f) ;

  if (Idst->pixel_format != Isrc->pixel_format)
    Iout = ImageAlloc(lmi->nspokes, lmi->nrings, Isrc->pixel_format,1);
  else
    Iout = Idst ;

  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
LogMapOffsetOrientation(LOGMAP_INFO *lmi, int wsize,IMAGE *Isrc,IMAGE *Iorient)
{
  static IMAGE *Itmp ;
  IMAGE  *Iout, *Ix, *Iy ;
  int    rows, cols ;

  rows = Isrc->rows ;
  cols = Isrc->cols ;

  if (!ImageCheckSize(Isrc, Itmp, 0, 0, 0))
  {
    if (Itmp)
      ImageFree(&Itmp) ;
    Itmp = ImageAlloc(rows, cols, PFFLOAT, 1) ;
  }
  else
    ImageSetSize(Itmp, rows, cols) ;

  LogMapMeanFilter(lmi, Isrc, Itmp) ;
  Isrc = Itmp ;

  if (!Iorient)
    Iorient = ImageAlloc(lmi->nspokes,lmi->nrings, Isrc->pixel_format,2);
  else
    ImageClearArea(Iorient, 0, 0, Iorient->rows, Iorient->cols, 0.0f) ;

  if (Iorient->pixel_format != Isrc->pixel_format)
    Iout = ImageAlloc(lmi->nspokes, lmi->nrings, Isrc->pixel_format,2);
  else
    Iout = Iorient ;

  Ix = ImageAllocHeader(rows, cols, PFFLOAT, 1) ;
  Iy = ImageAllocHeader(rows, cols, PFFLOAT, 1) ;
  Ix->image = Iorient->image ;
  Iy->image = Iorient->image + Ix->sizeimage ;

  LogMapSobel(lmi, Isrc, NULL, Ix, Iy, 1, 0, lmi->nrings-1) ;

  ImageFree(&Ix) ;
  ImageFree(&Iy) ;
  if (Iout != Iorient)
  {
    ImageCopy(Iout, Iorient) ;
    ImageFree(&Iout) ;
  }
  return(Iorient) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
LogMapOffsetDirection(LOGMAP_INFO *lmi, IMAGE *Iorient, IMAGE *Ioffset)
{
  IMAGE  *Iout ;
  float  phi, rho, sin_phi, cos_phi, e_to_the_rho, ox=0.0f, oy=0.0f, dx, dy ;
  int    x0, y0, rows, cols, x, y, whalf, k, ring, spoke ;
  float  *xpix, *ypix, ldx, ldy, *or_xpix,*or_ypix, *oxpix, *oypix, dir, lox, 
         loy, nrho, nphi = 0.0f, sin_nphi, cos_nphi, e_to_the_nrho, dot ;
  LOGPIX *npix ;


  rows = Iorient->rows ;
  cols = Iorient->cols ;

  if (!Ioffset)
    Ioffset = ImageAlloc(rows, cols, Iorient->pixel_format, 2) ;
  else
    ImageClearArea(Ioffset, 0, 0, rows, cols, 0.0f) ;

  if (Ioffset->pixel_format != Iorient->pixel_format)
    Iout = ImageAlloc(rows, cols, Iorient->pixel_format, 2) ;
  else
    Iout = Ioffset ;

  whalf = (3-1)/2 ;                       /* fix window size to 3x3 */
  xpix = IMAGEFpix(Iorient, 0, 0) ;
  ypix = IMAGEFseq_pix(Iorient, 0, 0, 1) ;
  or_xpix = IMAGEFpix(Iorient, 0, 0) ;
  or_ypix = IMAGEFseq_pix(Iorient, 0, 0, 1) ;
  oxpix = IMAGEFpix(Iout, 0, 0) ;
  oypix = IMAGEFseq_pix(Iout, 0, 0, 1) ;
  for (y0 = 0 ; y0 < rows ; y0++)
  {
    for (x0 = 0 ; x0 < cols ; x0++, xpix++, ypix++, oxpix++, oypix++)
    {
/*
  Now calculate the orientation for this point by averaging local gradient
  orientation within the specified window.

  x and y are in window coordinates, while xc and yc are in image
  coordinates.
*/
      /* calculate orientation vector */
      lox = *or_xpix++ ;   /* vector in log space */
      loy = *or_ypix++ ;
      if (LOG_PIX_AREA(lmi, x0, y0) <= 0)
        continue ;

      rho = LOG_PIX_RHO(lmi, x0, y0) ;
      phi = LOG_PIX_PHI(lmi, x0, y0) ;
#define USE_JACOBIAN 0
#if USE_JACOBIAN
      sin_phi = sin(phi) ;
      cos_phi = cos(phi) ;
      e_to_the_rho = exp(rho) ;
#else
      sin_phi = 0.0f ;
      cos_phi = 1.0f ;
      e_to_the_rho = 1.0f ;
#endif

      /* use Jacobian to convert to Cartesian vector */
      ox = e_to_the_rho * (lox * cos_phi - loy * sin_phi) ;
      oy = e_to_the_rho * (lox * sin_phi + loy * cos_phi) ;
      dir = 0.0f ;
      for (k = 0, y = -whalf ; y <= whalf ; y++)
      {
        for (x = -whalf ; x <= whalf ; x++, k++)
        {
          npix = LOG_PIX_NBD(lmi, x0, y0, k) ;
          ring = npix->ring ; spoke = npix->spoke ;
          nrho = npix->rho ; nphi = npix->phi ;
#if USE_JACOBIAN
          cos_nphi = cos(nphi) ; sin_nphi = sin(nphi) ;
          e_to_the_nrho = exp(nrho) ;
#else
          e_to_the_nrho = 1.0f ;
          cos_nphi = 1.0f ; sin_nphi = 0.0f ;
#endif
          ldx = *IMAGEFpix(Iorient, ring, spoke) ;
          ldy = *IMAGEFseq_pix(Iorient, ring, spoke, 1) ;

          /* use Jacobian of mapping to transform to Cartesian vectors */
          dx = e_to_the_nrho * (ldx * cos_nphi - ldy * sin_nphi) ;
          dy = e_to_the_nrho * (ldx * sin_nphi + ldy * cos_nphi) ;
          dot = dx*ox + dy*oy ;
          if (dot < 0.0f)
            dot = 0.0f ;
          dir += (x*ox + y*oy) * dot ;
        }
      }

      if (ISZERO(dir))
        lox = loy = 0.0f ;
      else if (dir > 0.0f)   /* flip by 180 */
      {
        lox = -lox ;
        loy = -loy ;
      }
      *oxpix = lox ;
      *oypix = loy ;
    }

  }

  if (Iout != Ioffset)
  {
    ImageCopy(Iout, Ioffset) ;
    ImageFree(&Iout) ;
  }
  return(Ioffset) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
LogMapOffsetMagnitude(LOGMAP_INFO *lmi,IMAGE *Isrc,IMAGE *Idst,int maxsteps)
{
  int     rows, cols, x, y, ax, ay, sx, sy, x1, y1, dx, dy, odx, ody, d,
          steps, dot = 0, xoff, yoff, k, xold, yold ;
  float   *src_xpix, *src_ypix, *dst_xpix, *dst_ypix, *oxpix, *oypix ;
  LOGPIX  *pix, *npix = NULL ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols,Isrc->pixel_format,
                      Isrc->num_frame);

  rows = Isrc->rows ;
  cols = Isrc->cols ;

  src_xpix = IMAGEFpix(Isrc, 0, 0) ;
  src_ypix = IMAGEFseq_pix(Isrc, 0, 0, 1) ;
  dst_xpix = IMAGEFpix(Idst, 0, 0) ;
  dst_ypix = IMAGEFseq_pix(Idst, 0, 0, 1) ;
  for (y1 = 0 ; y1 < rows ; y1++)
  {
    for (x1 = 0 ; x1 < cols ; x1++)
    {
      /* do a Bresenham algorithm do find the offset line at this point */
      dx = nint(*src_xpix++ * FSCALE) ;
      dy = nint(*src_ypix++ * FSCALE) ;
      xold = x = x1 ;
      yold = y = y1 ;
      ax = ABS(dx) << 1 ;
      sx = SGN(dx) ;
      ay = ABS(dy) << 1 ;
      sy = SGN(dy) ;
      
      if (x1 == 31 && y1 == 22)
        DiagBreak() ;

      pix = LOG_PIX(lmi, x1, y1) ;
      if (ax > ay)  /* x dominant, move sx in x each time step, check for y */
      {
        d = ay - (ax >> 1) ;
        for (steps = 0 ; steps < maxsteps ; steps++)
        {
          oxpix = IMAGEFpix(Isrc, x, y) ;
          oypix = IMAGEFseq_pix(Isrc, x, y, 1) ;
          odx = nint(*oxpix * FSCALE) ;
          ody = nint(*oypix * FSCALE) ;
          dot = odx * dx + ody * dy ;

          if (dot <= 0)
            break ;
          if (d >= 0)               /* move in x and y */
          {
            yoff = sy ;
            d -= ax ;
          }
          else                      /* move only in x */
            yoff = 0 ;

          k = 4+3*yoff + sx ;         /* find index of appropriate neighbor */
          npix = pix->nbd[k] ;

          xold = x ;
          yold = y ;
          x = npix->ring ;
          y = npix->spoke ;
          d += ay ;
          pix = npix ;
        }
      }
      else    /* y dominant, move sy in y each time step, check for x */
      {
        d = ax - (ay >> 1) ;
        for (steps = 0 ; steps < maxsteps ; steps++)
        {
          oxpix = IMAGEFpix(Isrc, x, y) ;
          oypix = IMAGEFseq_pix(Isrc, x, y, 1) ;
          odx = nint(*oxpix * FSCALE) ;
          ody = nint(*oypix * FSCALE) ;
          dot = odx * dx + ody * dy ;
          if (dot <= 0)
            break ;
          if (d >= 0)   /* move one in x direction */
          {
            xoff = sx ;
            d -= ay ;
          }
          else         /* only move in y */
            xoff = 0 ;

          k = 4+3*(sy) + xoff ;      /* find index of appropriate neighbor */
          npix = pix->nbd[k] ;

          xold = x ;
          yold = y ;
          x = npix->ring ;
          y = npix->spoke ;
          d += ax ;
          pix = npix ;
        }
      }

      if (dot == 0)  /* zero of vector field, not reversal */
      {
        xold = x ;
        yold = y ;
      }

      *dst_xpix++ = (float)(xold - x1) ;
      *dst_ypix++ = (float)(yold - y1) ;
    }
  }

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
LogMapApplyOffset(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Ioffset, IMAGE *Idst)
{
  return(ImageApplyOffset(Isrc, Ioffset, Idst)) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int compare_sort_array(const void *pf1, const void *pf2) ;
IMAGE *
LogMapMedianFilter(LOGMAP_INFO *lmi, IMAGE *Isrc, int wsize, IMAGE *Ioffset, 
                   IMAGE *Idst)
{
  static float *sort_array = NULL ;
  static int   sort_size = 0 ;
  int    x0, y0, rows, cols, whalf, dx, dy, frame,wsq,
         median_index, k, nring, nspoke ;
  float  *sptr, *outPix ;
  LOGPIX   *npix ;
#define OVERLAPPING_FAST_MEDIAN 0
#if OVERLAPPING_FAST_MEDIAN
  float  min_val, max_val ;
  int    ecode ;
  byte   *in_image, *out_image ;
  IMAGE  *Iout, *Iin ;
#endif

  rows = Isrc->rows ;
  cols = Isrc->cols ;

  if (ImageCheckSize(Isrc, Idst, 0, 0, 0))
  {
    if (Idst)
      ImageFree(&Idst) ;
    Idst = ImageAlloc(rows, cols, PFFLOAT, 1) ;
  }
  else
    ImageSetSize(Idst, rows, cols) ;

#if OVERLAPPING_FAST_MEDIAN
  if (!Ioffset)
  {
    /* h_median only takes byte formatted input and output images */
    if (Isrc->pixel_format != PFBYTE)
    {
      Iin = ImageAlloc(rows, cols,PFBYTE, Isrc->num_frame);
      ImageValRange(Isrc, &min_val, &max_val) ;
      ImageScale(Isrc, Isrc, 0.0f, 255.0f) ;
      ImageCopy(Isrc, Iin) ;
      ImageScale(Isrc, Isrc, min_val, max_val) ;  /* restore old image */
    }
    else
      Iin = Isrc ;

    if (Isrc->pixel_format != PFBYTE)
      Iout = ImageAlloc(rows, cols, PFBYTE, Idst->num_frame);
    else
      Iout = Idst ;

    out_image = Iout->image ;
    in_image = Iin->image ;
    for (frame = 0 ; frame < Isrc->num_frame ; frame++)
    {
      ecode = h_median(Iin, Iout, wsize) ;
      if (ecode != HIPS_OK)
        ErrorReturn(NULL, (ecode, 
                           "LogMapMedianFilter: h_median failed (%d)", ecode));

      Iout->firstpix += Iout->sizeimage ;
      Iout->image += Iout->sizeimage ;

      Iin->firstpix += Iin->sizeimage ;
      Iin->image += Iin->sizeimage ;
    }

    Iout->firstpix = Iout->image = out_image ;
    Iin->firstpix = Iin->image = in_image ;
    if (Iin != Isrc)
      ImageFree(&Iin) ;
    if (Idst != Iout)
    {
      ImageCopy(Iout, Idst) ;
      ImageScale(Idst, Idst, min_val, max_val) ;
      ImageFree(&Iout) ;
    }
    return(NO_ERROR) ;
  }
#endif

  median_index = wsize*wsize/2 ;
  wsq = wsize*wsize ;
  whalf = (wsize-1)/2 ;

  /* create a static array for sorting pixels in */
  if (wsize > sort_size)
  {
    sort_size = wsize ;
    if (sort_array)
      sort_array = NULL ;
  }

  if (!sort_array)
    sort_array = (float *)calloc(wsq, sizeof(float)) ;

  for (frame = 0 ; frame < Isrc->num_frame ; frame++)
  {
    outPix = IMAGEFseq_pix(Idst, 0, 0, frame) ;
    for (y0 = 0 ; y0 < rows ; y0++)
    {
      for (x0 = 0 ; x0 < cols ; x0++, outPix++)
      {
        if (LOG_PIX_AREA(lmi, x0, y0) <= 0)
          continue ;

        if (Ioffset)
        {
          dx = nint(*IMAGEFpix(Ioffset, x0, y0)) ;
          dy = nint(*IMAGEFseq_pix(Ioffset, x0, y0, 1)) ;
        }
        else
          dx = dy = 0 ;
        
        for (sptr = sort_array, k = 0 ; k < NBD_SIZE ; k++)
        {
          npix = LOG_PIX_NBD(lmi, x0+dx, y0+dy, k) ;
          if (!npix)
            continue ;   /* should never happen */
          nring = npix->ring ;
          nspoke = npix->spoke ;
          *sptr++ = *IMAGEFpix(Isrc, nring, nspoke) ;
        }
        qsort(sort_array, wsq, sizeof(float), compare_sort_array) ;
        *outPix = sort_array[median_index] ;
      }
    }
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
compare_sort_array(const void *pf1, const void *pf2)
{
  register float f1, f2 ;

  f1 = *(float *)pf1 ;
  f2 = *(float *)pf2 ;

/*  return(f1 > f2 ? 1 : f1 == f2 ? 0 : -1) ;*/
  if (f1 > f2)
    return(1) ;
  else if (f1 < f2)
    return(-1) ;

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              convolve a logmap with a Gaussian kernel using
              2 one dimensional convolutions
----------------------------------------------------------------------*/
IMAGE *
LogMapGaussianFilter(LOGMAP_INFO *lmi, IMAGE *Isrc, 
                       IMAGE *Igaussian, IMAGE *Idst)
{
  int    rows, cols ;
  IMAGE  *Iout, *Itmp ;

  rows = Isrc->rows ;
  cols = Isrc->cols ;
  if (!ImageCheckSize(Isrc, Idst, 0, 0, 0))
  {
    if (Idst)
      ImageFree(&Idst) ;
    Idst = ImageAlloc(lmi->nspokes,lmi->nrings, Isrc->pixel_format,1);
  }
  else
  {
    ImageSetSize(Idst, rows, cols) ;
    ImageClearArea(Idst, 0, 0, Idst->rows, Idst->cols, 0.0f) ;
  }

  if (Idst->pixel_format != Isrc->pixel_format)
    Iout = ImageAlloc(lmi->nspokes, lmi->nrings, Isrc->pixel_format,1);
  else
    Iout = Idst ;

  Itmp = ImageAlloc(rows, cols, PFFLOAT, 1) ;
  lmConvolve1d(lmi, Isrc, Igaussian, Itmp, IMAGE_VERTICAL) ;
  lmConvolve1d(lmi, Itmp, Igaussian, Iout, IMAGE_HORIZONTAL) ;


  ImageFree(&Itmp) ;
  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              convolve a logmap with a Gaussian kernel using
              2 one dimensional convolutions
----------------------------------------------------------------------*/
static int 
lmConvolve1d(LOGMAP_INFO *lmi, IMAGE *Isrc, IMAGE *Igaussian, IMAGE *Idst, 
             int which)
{
  int    rows, cols, ring, spoke, klen, khalf, nx, ny, k ;
  LOGPIX *pix, *npix ;
  float  *fdst, *kernel, out ;

  klen = Igaussian->cols ;
  khalf = (klen-1) / 2 ;
  rows = Isrc->rows ;
  cols = Isrc->cols ;
  pix = LOG_PIX(lmi, 0, 0) ;
  fdst = IMAGEFpix(Idst, 0, 0) ;
  if (which == IMAGE_VERTICAL)
  {
    for (spoke = 0 ; spoke < rows ; spoke++)
    {
      for (ring = 0 ; ring < cols ; ring++, pix++, fdst++)
      {
        if (pix->area <= 0)
          continue ;    /* not a real pixel */
        out = 0.0f ;

        /* first track from central pixel downwards (South) */
        kernel = IMAGEFpix(Igaussian, khalf, 0) ;
        npix = pix ;
        for (k = khalf ; k >= 0 ; k--)
        {
          nx = npix->ring ;
          ny = npix->spoke ;
          out += *IMAGEFpix(Isrc, nx, ny) * *kernel-- ;
          npix = npix->nbd[N_S] ;
        }

        /* now track upwards (North) */
        kernel = IMAGEFpix(Igaussian, khalf+1, 0) ;
        npix = pix->nbd[N_N] ;
        for (k = khalf+1 ; k < klen ; k++)
        {
          nx = npix->ring ;
          ny = npix->spoke ;
          out += *IMAGEFpix(Isrc, nx, ny) * *kernel++ ;
          npix = npix->nbd[N_N] ;
        }

        *fdst = out ;
      }
    }
  }
  else     /* horizontal convolution */
  {
    for (spoke = 0 ; spoke < rows ; spoke++)
    {
      for (ring = 0 ; ring < cols ; ring++, pix++, fdst++)
      {
        if (pix->area <= 0)
          continue ;    /* not a real pixel */

        out = 0.0f ;

        /* first track from central pixel leftwards */
        kernel = IMAGEFpix(Igaussian, khalf, 0) ;
        npix = pix ;
        for (k = khalf ; k >= 0 ; k--)
        {
          nx = npix->ring ;
          ny = npix->spoke ;
          out += *IMAGEFpix(Isrc, nx, ny) * *kernel-- ;
          npix = npix->nbd[N_W] ;
        }

        /* now track rightwards */
        kernel = IMAGEFpix(Igaussian, khalf+1, 0) ;
        npix = pix ;
        for (k = khalf+1 ; k < klen ; k++)
        {
          nx = npix->ring ;
          ny = npix->spoke ;
          out += *IMAGEFpix(Isrc, nx, ny) * *kernel++ ;
          npix = npix->nbd[N_E] ;
        }

        *fdst = out ;
      }

    }
  }

  return(NO_ERROR) ;
}
