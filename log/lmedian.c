/**
 * @file  lmedian.c
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
  MODIFIED BY Bruce Fischl to use Logmap data structure
  */

/*
 * Modified version (HIPS-2):
 * Copyright (c) 1991 Michael Landy
 *
 * Disclaimer:  No guarantees of performance accompany this software,
 * nor is any responsibility assumed on the part of the authors.  All the
 * software has been tested extensively and every effort has been made to
 * insure its reliability.
 */

/*  Original Copyright (c) 1987 Pierre Landau

Disclaimer:  No guarantees of performance accompany this software,
nor is any responsibility assumed on the part of the authors.  All the
software has been tested extensively and every effort has been made to
insure its reliability.   */

/*
 * h_median.c - apply a median filter to an image
 *
 * where size is the length of the side of the neighborhood in which the
 * median is computed.
 *
 * pixel formats: BYTE
 *
 * Mike Landy - 5/28/82
 * median algorithm replaced <Pierre Landau 1/6/87>
 * HIPS 2 - msl - 6/16/91
 */

#include <hipl_format.h>

#include "image.h"
#include "hips.h"
#include "h_logz.h"
#include "error.h"

static int log_sselect(int k, int *lo, int *hi) ;

static int log_median_b(LOGMAP_INFO *lmi, IMAGE *hdi,IMAGE *hdo, int size) ;
static int log_median_B(LOGMAP_INFO *lmi, byte *imagei,byte *imageo,int nr,
                        int nc,int nlpi,int nlpo,int size) ;


int
log_median(LOGMAP_INFO *lmi, IMAGE *hdi, IMAGE *hdo, int size) {
  switch (hdi->pixel_format) {
  case PFBYTE:
    return(log_median_b(lmi, hdi,hdo,size));
#if 0
  case PFINT:
    return(log_median_i(hdi,hdo,size));
  case PFFLOAT:
    return(log_median_f(lmi, hdi,hdo,size));
#endif
  default:
    return(perr(HE_FMTSUBR,"log_median", hformatname(hdi->pixel_format)));
  }
}


static int
log_median_b(LOGMAP_INFO *lmi, IMAGE *hdi,IMAGE *hdo, int size) {
  return(log_median_B(lmi, hdi->firstpix,hdo->firstpix,hdi->rows,hdi->cols,
                      hdi->ocols,hdo->ocols,size));
}

static int
log_median_B(LOGMAP_INFO *lmi, byte *imagei, byte *imageo, int nr, int nc,
             int nlpi, int nlpo,int size) {
  static h_boolean nballoc = FALSE;
  static int nbsize,*nb;
  static h_boolean colalloc = FALSE;
  static int *col,saverows,savecols;

  int sizesq,halfsz,ir,ic;
  int minus,plus,*np,top,bot,left,right,nexi,nexo,nextrow;
  byte  *ip,*op,*nnp;
  register int i,j,ii,jj;   /* loop counters */

  sizesq = size*size;
  halfsz = (sizesq + 1) / 2;
  plus = size / 2;
  minus = plus - size + 1;
  top = -minus;
  bot = nr - plus;
  left = -minus;
  right = nc - plus;
  nextrow = nlpi*minus + minus;
  if (!nballoc || nbsize < size) {
    if (nballoc)
      free(nb);
    nb = (int *) memalloc(size*size,sizeof(int));
    nballoc = TRUE;
    nbsize = size;
  }
  if (!colalloc || saverows < nr) {
    if (colalloc)
      free(col);
    col = (int *) memalloc(nr,sizeof(int));
    colalloc = TRUE;
    savecols = nlpi+1; /* force computation */
  }
  if (saverows != nr || savecols != nlpi) {
    saverows = nr;
    savecols = nlpi;
    for (ic = -nlpi,i=0;i<nr;i++)
      col[i] = (ic += nlpi); /* vector array */
  }
  ip = imagei;
  op = imageo;
  nexi = nlpi-nc;
  nexo = nlpo-nc;
  for (i=0;i<nr;i++) {
    for (j=0;j<nc;j++) {
      if (LOG_PIX_AREA(lmi, j, i) <= 0) {
        ip++ ;
        op++ ;
        continue ;
      }
      if (i<top || i>=bot || j<left || j>=right ) {
        np = nb;
        for (ii=minus;ii<=plus;ii++)
          for (jj=minus;jj<=plus;jj++) {
            ir = i + ii;
            ic = j + jj;
            ir = ir<0?0:(ir>=nr)?nr-1:ir;
            ic = ic<0?0:(ic>=nc)?nc-1:ic;
            *np++ = imagei[col[ir]+ic];
          }
      } else {
        nnp = ip + nextrow;
        np = nb;
        for (ii=minus;ii<=plus;ii++) {
          for (jj=minus;jj<=plus;jj++)
            *np++ = *nnp++;
          nnp += nlpi - size;
        }
      }
      ip++;
      *op++ = log_sselect(halfsz,nb,nb+sizesq-1);
    }
    ip += nexi;
    op += nexo;
  }
  return(HIPS_OK);
}
#define exchi(a,b) {tmpi = a; a = b; b = tmpi;} /* exchange ints */
#define exchp(a,b) {tmpp = a; a = b; b = tmpp;} /* exchange pointers */

/* select the k'th element from the list between lo and hi
 * this is a implementation of R.W.Floyd's improvement to Hoare's
 * original algorithm, as published in CACM, vol 18 no 3 (march 1975)
 */

static int
log_sselect(int k, int *lo, int *hi) {
  register int *i,*j;
  int *val[3],t,*tmpp,tmpi,df;

  while (1) {
    if (hi == lo)
      return(*lo);
    if ((t = hi-lo) <= 2) /* if the sequence is short (n<3) sort it directly */
    {
      val[0] = lo;
      val[1] = lo+1;
      val[2] = hi;
      if (t == 1) {
        return(*val[0] < *val[1] ? *val[k-1] : *val[2-k]);
      } else {
        if (*val[0] > *val[1]) exchp(val[0],val[1])
          if (*val[0] > *val[2]) exchp(val[0],val[2])
            if (*val[1] > *val[2]) exchp(val[1],val[2])
              return (*val[k-1]);
      }
    }
    else {
      t = *lo; /* take first element of list as pivot */
      i = lo;
      j = hi;
      if (*hi > t)
        exchi(*hi,*lo) /* set up for first exchange */
        while (i < j) {
          exchi(*i,*j)
          i++;
          j--;
          while (*i < t) i++;
          /* scan list for pair to exchange */
          while (*j > t)
            j--;
        }
      if (*lo == t)
        exchi(*lo,*j)
        /* put pivot back where it belongs */
        else {
          j++;
          exchi(*j,*hi)
        }

      /* now adjust hi,lo so they surround the subset
         containing the k-l+1th element */

      df = j-lo+1;
      if (df < k) {
        k = k-(df);
        lo = j+1;
      } else {
        if (df == k)
          return (*j);
        else
          hi = j - 1;
      }
    }
  }
}

#if 0
static float log_sselect_f(int k, float *lo, float *hi) ;

static float
log_sselect_f(int k, float *lo, float *hi) {
  register float *i,*j;
  float *val[3],t,*tmpp,tmpi,df;

  while (1) {
    if (hi == lo)
      return(*lo);
    if ((t = hi-lo) <= 2) /* if the sequence is short (n<3) sort it directly */
    {
      val[0] = lo;
      val[1] = lo+1;
      val[2] = hi;
      if (t == 1) {
        return(*val[0] < *val[1] ? *val[k-1] : *val[2-k]);
      } else {
        if (*val[0] > *val[1]) exchp(val[0],val[1])
          if (*val[0] > *val[2]) exchp(val[0],val[2])
            if (*val[1] > *val[2]) exchp(val[1],val[2])
              return (*val[k-1]);
      }
    }
    else {
      t = *lo; /* take first element of list as pivot */
      i = lo;
      j = hi;
      if (*hi > t)
        exchi(*hi,*lo) /* set up for first exchange */
        while (i < j) {
          exchi(*i,*j)
          i++;
          j--;
          while (*i < t) i++;
          /* scan list for pair to exchange */
          while (*j > t)
            j--;
        }
      if (*lo == t)
        exchi(*lo,*j)
        /* put pivot back where it belongs */
        else {
          j++;
          exchi(*j,*hi)
        }

      /* now adjust hi,lo so they surround the subset
         containing the k-l+1th element */

      df = j-lo+1;
      if (df < k) {
        k = k-(df);
        lo = j+1;
      } else {
        if (df == k)
          return (*j);
        else
          hi = j - 1;
      }
    }
  }
}
static int log_median_f(LOGMAP_INFO *lmi, IMAGE *hdi,IMAGE *hdo, int size) ;
static int log_median_F(LOGMAP_INFO *lmi, float *imagei, float *imageo, int nr,
                        int nc, int nlpi, int nlpo,int size) ;

static int log_median_i(IMAGE *hdi,IMAGE *hdo, int size) ;
static int log_median_I(int *imagei, int *imageo, int nr, int nc, int nlpi,
                        int nlpo,int size) ;

static int
log_median_f(LOGMAP_INFO *lmi, IMAGE *hdi,IMAGE *hdo, int size) {
  return(log_median_F(lmi, (float *)hdi->firstpix, (float *)hdo->firstpix,
                      hdi->rows,hdi->cols, hdi->ocols,hdo->ocols,size));
}

static int
log_median_F(LOGMAP_INFO *lmi, float *imagei,float *imageo,int nr,int nc,
             int nlpi,int nlpo,int size) {
  static h_boolean nballoc = FALSE;
  static float *nb;
  static h_boolean colalloc = FALSE;
  static int *col,saverows,savecols, nbsize;

  int sizesq,halfsz,ir,ic;
  int minus,plus,top,bot,left,right,nexi,nexo,nextrow;
  float *ip,*op,*nnp, *np;
  register int i,j,ii,jj;   /* loop counters */

  sizesq = size*size;
  halfsz = (sizesq + 1) / 2;
  plus = size / 2;
  minus = plus - size + 1;
  top = -minus;
  bot = nr - plus;
  left = -minus;
  right = nc - plus;
  nextrow = nlpi*minus + minus;
  if (!nballoc || nbsize < size) {
    if (nballoc)
      free(nb);
    nb = (float *) memalloc(size*size,sizeof(float));
    nballoc = TRUE;
    nbsize = size;
  }
  if (!colalloc || saverows < nr) {
    if (colalloc)
      free(col);
    col = (int *) memalloc(nr,sizeof(int));
    colalloc = TRUE;
    savecols = nlpi+1; /* force computation */
  }
  if (saverows != nr || savecols != nlpi) {
    saverows = nr;
    savecols = nlpi;
    for (ic = -nlpi,i=0;i<nr;i++)
      col[i] = (ic += nlpi); /* vector array */
  }
  ip = imagei;
  op = imageo;
  nexi = nlpi-nc;
  nexo = nlpo-nc;
  for (i=0;i<nr;i++) {
    for (j=0;j<nc;j++) {
      if (i<top || i>=bot || j<left || j>=right) {
        np = nb;
        for (ii=minus;ii<=plus;ii++)
          for (jj=minus;jj<=plus;jj++) {
            ir = i + ii;
            ic = j + jj;
            ir = ir<0?0:(ir>=nr)?nr-1:ir;
            ic = ic<0?0:(ic>=nc)?nc-1:ic;
            *np++ = imagei[col[ir]+ic];
          }
      } else {
        nnp = ip + nextrow;
        np = nb;
        for (ii=minus;ii<=plus;ii++) {
          for (jj=minus;jj<=plus;jj++)
            *np++ = *nnp++;
          nnp += nlpi - size;
        }
      }
      ip++;
      *op++ = log_sselect_f(halfsz,nb,nb+sizesq-1);
    }
    ip += nexi;
    op += nexo;
  }
  return(HIPS_OK);
}

static int
log_median_i(IMAGE *hdi,IMAGE *hdo, int size) {
  return(log_median_I((int *)hdi->firstpix, (int *)hdo->firstpix,
                      hdi->rows,hdi->cols, hdi->ocols,hdo->ocols,size));
}

static int
log_median_I(int *imagei,int *imageo,int nr,int nc,int nlpi,int nlpo,int size) {
  static h_boolean nballoc = FALSE;
  static int nbsize,*nb;
  static h_boolean colalloc = FALSE;
  static int *col,saverows,savecols;

  int sizesq,halfsz,ir,ic;
  int minus,plus,*np,top,bot,left,right,nexi,nexo,nextrow;
  int *ip,*op,*nnp;
  register int i,j,ii,jj;   /* loop counters */

  sizesq = size*size;
  halfsz = (sizesq + 1) / 2;
  plus = size / 2;
  minus = plus - size + 1;
  top = -minus;
  bot = nr - plus;
  left = -minus;
  right = nc - plus;
  nextrow = nlpi*minus + minus;
  if (!nballoc || nbsize < size) {
    if (nballoc)
      free(nb);
    nb = (int *) memalloc(size*size,sizeof(int));
    nballoc = TRUE;
    nbsize = size;
  }
  if (!colalloc || saverows < nr) {
    if (colalloc)
      free(col);
    col = (int *) memalloc(nr,sizeof(int));
    colalloc = TRUE;
    savecols = nlpi+1; /* force computation */
  }
  if (saverows != nr || savecols != nlpi) {
    saverows = nr;
    savecols = nlpi;
    for (ic = -nlpi,i=0;i<nr;i++)
      col[i] = (ic += nlpi); /* vector array */
  }
  ip = imagei;
  op = imageo;
  nexi = nlpi-nc;
  nexo = nlpo-nc;
  for (i=0;i<nr;i++) {
    for (j=0;j<nc;j++) {
      if (i<top || i>=bot || j<left || j>=right) {
        np = nb;
        for (ii=minus;ii<=plus;ii++)
          for (jj=minus;jj<=plus;jj++) {
            ir = i + ii;
            ic = j + jj;
            ir = ir<0?0:(ir>=nr)?nr-1:ir;
            ic = ic<0?0:(ic>=nc)?nc-1:ic;
            *np++ = imagei[col[ir]+ic];
          }
      } else {
        nnp = ip + nextrow;
        np = nb;
        for (ii=minus;ii<=plus;ii++) {
          for (jj=minus;jj<=plus;jj++)
            *np++ = *nnp++;
          nnp += nlpi - size;
        }
      }
      ip++;
      *op++ = log_sselect(halfsz,nb,nb+sizesq-1);
    }
    ip += nexi;
    op += nexo;
  }
  return(HIPS_OK);
}


#endif
