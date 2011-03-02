/**
 * @file  artmap.c
 * @brief Adaptive Resonance Theory map
 *
 * some kind of vestige of the CNS department
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:42 $
 *    $Revision: 1.5 $
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

/*-----------------------------------------------------------------
              INCLUDE FILES
-----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "machine.h"
#include "artmap.h"
#include "diag.h"
#include "error.h"
#include "proto.h"

/*-----------------------------------------------------------------
              MACROS AND CONSTANTS
-----------------------------------------------------------------*/

#define InsCalloc        calloc
#define InsHalloc(a,b)   calloc((int)(a), b)
#define InsFree          free
#define InsHfree         free

#define BETA         .001
#define DELTA_RHO    .0001  /* amount above Tj to set row after mismatch */
#define ALPHA(art, j)  (1.0 / (art->beta + (double)art->ninputs)) - \
                       (.00001 * ((double)j + 1.0))
#ifndef MIN
#define MIN(a,b)      (a < b ? a : b)
#endif

#ifndef MAX
#define MAX(a,b)      (a > b ? a : b)
#endif


/*----------------------------------------------------------------------
                STRUCTURES
----------------------------------------------------------------------*/

/*-----------------------------------------------------------------
                PROTOTYPES
-----------------------------------------------------------------*/

static void    artInitWeights(ARTMAP *artmap) ;
static int     artFeedBack(ARTMAP *artmap, int class) ;
static int     artFeedForward(ARTMAP *artmap) ;
static double  artChoice(ARTMAP *art, int j) ;
static double  norm(double huge *a, int len) ;
static void    intersect(double huge *a,double huge *b,double huge *c,int len);
static void    artFastLearn(ARTMAP *artmap, int j) ;
int      artProcess(ARTMAP *artmap, double *I) ;
int      mapFeedForward(ARTMAP *artmap, int class) ;
void    mapLearn(ARTMAP *artmap, int ja, int jb) ;

/*-----------------------------------------------------------------
                   MACROS
-----------------------------------------------------------------*/

#define Mzj(artmap, i, j)   (artmap->zj + \
                              (((long)j) * (long)artmap->ninputs) + (long)i)
#define Mwj(artmap, j, k)   (artmap->w + \
                              (((long)j) * (long)artmap->noutputs) + (long)k)

/*-----------------------------------------------------------------
                STATIC DATA
-----------------------------------------------------------------*/

/*-----------------------------------------------------------------
                GLOBAL DATA
-----------------------------------------------------------------*/

/*-----------------------------------------------------------------
                FUNCTIONS
-----------------------------------------------------------------*/

/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
ARTMAP *
ArtmapAlloc(int ninputs, int noutputs, double rho_bar, int max_f2)
{
  ARTMAP  *artmap ;

  artmap = (ARTMAP *)InsCalloc(1, sizeof(ARTMAP)) ;

  artmap->ninputs = ninputs ;
  artmap->noutputs = noutputs ;
  artmap->max_f2 = max_f2 ;
  artmap->beta = BETA ;
  artmap->rho = artmap->rho_bar = rho_bar ;
  artmap->scratch = (double *)InsCalloc(ninputs, sizeof(double)) ;
  artmap->f0 = (double *)InsCalloc(ninputs, sizeof(double)) ;
  artmap->f1 = (double *)InsCalloc(ninputs, sizeof(double)) ;
  artmap->f2 = (double *)InsCalloc(max_f2, sizeof(double)) ;
  artmap->flags = (int *)InsCalloc(max_f2, sizeof(int)) ;
  artmap->zj =
    (double huge *)InsHalloc((long)ninputs * (long)max_f2, sizeof(double)) ;
  artmap->w =
    (int huge *)InsHalloc((long)noutputs * (long)max_f2, sizeof(int)) ;

  artInitWeights(artmap) ;

  return(artmap) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
ARTMAP *
ArtmapRead(char *fname)
{
  ARTMAP  *artmap ;
  int     ninputs, noutputs, f2nodes, i, j, k, wj ;
  double  zj ;
  FILE    *fp ;

  double  beta, rho_bar ;

  fp = fopen(fname, "r") ;
  if (!fp)
    return(NULL) ;

  fscanf(fp, "%d %d %d %lf %lf\n",
         &ninputs, &noutputs, &f2nodes, &beta, &rho_bar) ;

  artmap = ArtmapAlloc(ninputs, noutputs, rho_bar, f2nodes) ;
  artmap->ncommitted = artmap->f2nodes = f2nodes ;
  for (j = 0 ; j < artmap->f2nodes ; j++)
  {
    artmap->flags[j] |= ARTMAP_COMMITTED ;
    for (k = 0 ; k < artmap->noutputs ; k++)
    {
      fscanf(fp, "%d\n", &wj) ;
      *Mwj(artmap, j, k) = wj ;
    }
  }

  for (i = 0 ; i < artmap->ninputs ; i++)
  {
    for (j = 0 ; j < artmap->f2nodes ; j++)
    {
      fscanf(fp, "%lf\n", &zj) ;
      *Mzj(artmap, i, j) = zj ;
    }
  }

  fclose(fp) ;
  return(artmap) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
ArtmapWrite(ARTMAP *artmap, char *fname)
{
  FILE  *fp ;
  int    i, j, k, wj ;
  double zj ;

  fp = fopen(fname, "w") ;
  if (!fp)
    return(-1) ;

  fprintf(fp, "%d %d %d %f %f\n",
          artmap->ninputs, artmap->noutputs, artmap->f2nodes, artmap->beta,
          artmap->rho_bar) ;

  for (j = 0 ; j < artmap->f2nodes ; j++)
  {
    for (k = 0 ; k < artmap->noutputs ; k++)
    {
      wj = *Mwj(artmap, j, k) ;
      fprintf(fp, "%d\n", wj) ;
    }
  }

  for (i = 0 ; i < artmap->ninputs ; i++)
  {
    for (j = 0 ; j < artmap->f2nodes ; j++)
    {
      zj = *Mzj(artmap, i, j) ;
      fprintf(fp, "%f\n", zj) ;
    }
  }

  fclose(fp) ;
  return(0) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
ArtmapFree(ARTMAP **artmap)
{
  InsFree((*artmap)->scratch) ;
  InsFree((*artmap)->f0) ;
  InsFree((*artmap)->f1) ;
  InsFree((*artmap)->f2) ;
  InsFree((*artmap)->flags) ;
  InsHfree((*artmap)->zj) ;
  InsHfree((*artmap)->w) ;
  InsFree(*artmap) ;
  *artmap = NULL ;
  return(0) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
ArtmapProcess(ARTMAP *artmap, double *I)
{
  int     class ;

  if (Gdiag & DIAG_WRITE) printf("ArtmapProcess()\n") ;
  class = artProcess(artmap, I) ;

  if (class == 4)
  {
    int k ;

    if (Gdiag & DIAG_WRITE)
      for (k = 0 ; k < artmap->noutputs ; k++)
        printf("Mwj(4, %d) = %d\n", k, *Mwj(artmap, 4, k)) ;
  }
  return(mapFeedForward(artmap, class)) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
artProcess(ARTMAP *artmap, double *I)
{
  int     nclass, i, class = 0 ;

  for (i = 0 ; i < artmap->ninputs ; i++)
    artmap->f0[i] = artmap->f1[i] = I[i] ;

  for (i = 0 ; i < artmap->f2nodes ; i++)
    artmap->flags[i] &= ~ARTMAP_RESET ;

  nclass = 0 ;
  do
  {
    /*
      in test mode, if the first neither of the 1st two classes
      survive feedback, give it up.
    */
    if (!artmap->learn && ++nclass > 2)
    {
      artmap->match = 0.0 ;
      break ;
    }
    artmap->class = class = artFeedForward(artmap) ;
    if (class < 0) break ;    /* no classes left */
  }
  while (!artFeedBack(artmap, class)) ;

  if (Gdiag & DIAG_WRITE) printf("artProcess() --> %d\n", class) ;
  return(class) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
ArtmapLearn(ARTMAP *artmap, double *I, int class)
{
  int   artClass ;

  if (Gdiag & DIAG_WRITE) printf("ArtmapLearn(%d)\n", class) ;

  artmap->learn = 1 ;      /* will allow committment of new nodes */

  artmap->match = artmap->rho_bar - DELTA_RHO ;
  do
  {
    artmap->rho = artmap->match + DELTA_RHO ;  /* match tracking */
    /*    artmap->f2nodes = artmap->max_f2 ;*/
    artClass = ArtmapProcess(artmap, I) ;
    /*    artmap->f2nodes = artmap->ncommitted ;*/
    if (artClass < 0)      /* commit new node */
    {
      if (artmap->f2nodes < artmap->max_f2)
      {
        if (Gdiag & DIAG_WRITE)
          printf("committing new f2 node %d to class %d\n",
                 artmap->f2nodes, class) ;
        artmap->class = artmap->f2nodes++ ;
        artmap->flags[artmap->class] |= ARTMAP_COMMITTED ;
        artmap->ncommitted++ ;
        break ;
      }
      else return(-1) ;    /* out of memory */
    }
    if (artClass != class)
    {
      if (Gdiag & DIAG_WRITE)
        printf("artClass %d != class %d --> match tracking\n",
               artClass, class);
    }
  }
  while (artClass != class) ;

  artFastLearn(artmap, artmap->class) ;
  mapLearn(artmap, artmap->class, class) ;

  artmap->rho = artmap->rho_bar ;  /* relax to baseline vigilance */
  return(artClass) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
artInitWeights(ARTMAP *artmap)
{
  int     i, j, k ;
  double  *zj ;

  for (i = 0 ; i < artmap->ninputs ; i++)
  {
    for (j = 0 ; j < artmap->max_f2 ; j++)
    {
      zj = Mzj(artmap, i, j) ;
      *zj = 1.0 ;
    }
  }
  for (j = 0 ; j < artmap->max_f2 ; j++)
  {
    for (k = 0 ; k < artmap->noutputs ; k++)
    {
      *Mwj(artmap, j, k) = 1 ;
    }
  }
}


/*----------------------------------------------------------------------
    Parameters:

    Description:
      feed the active f2 nodes activation back through the top down weights
      to f1.  Return 1 if resonance is established, 0 if there is a mismatch.

       Resonance is established if:

       rho <= |I ^ zj|
              --------
                |I|

    Returns:
----------------------------------------------------------------------*/
static int
artFeedBack(ARTMAP *artmap, int class)
{
  int    j, ninputs ;
  double match, norm_I_int_zj, norm_I ;

  j = class ;
  ninputs = artmap->ninputs ;
  intersect(artmap->f0, Mzj(artmap, 0, j), artmap->scratch, ninputs) ;
  norm_I_int_zj = norm(artmap->scratch, ninputs) ;
  norm_I = norm(artmap->f0, ninputs) ;
  if (!norm_I)
    match = 1 ;
  else
    match = norm_I_int_zj / norm_I ;
  artmap->match = match ;      /* for use in match tracking */

  if (Gdiag & DIAG_WRITE)
    printf("artFeedBack(%d) --> match %2.3f, n(I) %2.3f, n(I^Zj) %2.3f, "
           "rho %2.3f %s\n", class, match, norm_I, norm_I_int_zj, artmap->rho,
           match >= artmap->rho ? "RESONANCE" : "MISMATCH") ;

  if (match < artmap->rho)   /* mismatch --> do reset */
  {
    artmap->flags[j] |= ARTMAP_RESET ;
    return(0) ;
  }
  else
    return(1) ;        /* resonance */
}


/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static int
artFeedForward(ARTMAP *artmap)
{
  int       max_j, j ;
  double  f2, max_out ;

  max_j = -1 ;
  max_out = -1.0 ;

  for (j = 0 ; j < artmap->f2nodes ; j++)
  {
    if (artmap->flags[j] & ARTMAP_RESET)
      continue ;      /* this node was previous winner */

    f2 = artmap->f2[j] = artChoice(artmap, j) ;
    if (artmap->f2[j] > max_out)
    {
      max_j = j ;
      max_out = artmap->f2[j] ;
    }
  }

  if (Gdiag & DIAG_WRITE)
    printf("artFeedForward() --> %d (%2.3f)\n", max_j, max_out) ;
  return(max_j) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:
       calculate the choice function for an f2 node:

       Tj =    |I ^ zj|  = Zj . I
            ------------
      BETA + |zj|

    Returns:
----------------------------------------------------------------------*/
static double
artChoice(ARTMAP *artmap, int j)
{
  double  Tj, *zj, norm_zj, norm_I_int_zj ;

  zj = Mzj(artmap, 0, j) ;    /* address of jth nodes weights */

  /* calculate I ^ zj */
  if (artmap->flags[j] & ARTMAP_COMMITTED)
  {
    intersect(artmap->f0, zj, artmap->scratch, artmap->ninputs) ;
    norm_I_int_zj = norm(artmap->scratch, artmap->ninputs) ;
    norm_zj = norm(zj, artmap->ninputs) ;

    Tj = (double)norm_I_int_zj / (artmap->beta + (double)norm_zj) ;
  }
  else   /* Tj = |I| * ALPHA(j) */
  {
    Tj = (double)norm(artmap->f0, artmap->ninputs) ;
    Tj *= ALPHA(artmap, j) ;
  }
  return(Tj) ;
}


/*
   take the intersection of vectors a and b, and return it
   in vector c.  All vectors have length 'len'.
*/
static void
intersect(double huge *a, double huge *b, double huge *c, int len)
{
  register int i ;

  for (i = 0 ; i < len ; i++, a++, b++, c++)
    *c = MIN(*a, *b) ;
}


/*
   take the norm (city block) of vector a of length 'len'.
*/
static double
norm(double huge *a, int len)
{
  register int i ;
  double   norm_val ;

  for (norm_val = 0.0, i = 0 ; i < len ; i++, a++)
    norm_val += fabs(*a) ;

  return(norm_val) ;
}


/*
   f1->f2 learning:
     Zj =    (I ^ zj)
          ------------
          |I ^ zj| + B

   f2->f1 learning:
     zj = I ^ zj
*/
static void
artFastLearn(ARTMAP *artmap, int j)
{
  double  norm_I_int_zj ;
  int     i ;

  if (j < 0) return ;   /* not a valid class */

  intersect(artmap->f0, Mzj(artmap, 0, j), artmap->scratch, artmap->ninputs) ;
  norm_I_int_zj = norm(artmap->scratch, artmap->ninputs) ;
  for (i = 0 ; i < artmap->ninputs ; i++)
  {
    *Mzj(artmap, i, j) = artmap->scratch[i] ;
  }
}


/*----------------------------------------------------------------------
    Parameters:

    Description:
      calculate map field activation x = Yb ^ wj

      where Yb is b field activation, and wj are adaptive weights
      between F2a and the map field.

    Returns:
----------------------------------------------------------------------*/
int
mapFeedForward(ARTMAP *artmap, int class)
{
  int   k, apredict ;

  if (class < 0)
    return(-1) ;  /* ARTa did not make a choice */

  apredict = -1 ;
  for (k = 0 ; k < artmap->noutputs ; k++)
  {
    if ((apredict < 0) && *Mwj(artmap, class, k))    /* 1st node */
    {
      apredict = k ;
      break ;
    }
  }
  if (Gdiag & DIAG_WRITE)
    printf("mapFeedForward(%d) --> %d\n", class, apredict) ;
  return(apredict) ;
}


/*----------------------------------------------------------------------
    Parameters:

    Description:
      learn weights at the map field.  The weights obey

      Wjk = 1 for j and k active.  All other weights are 0.

    Returns:
----------------------------------------------------------------------*/
void
mapLearn(ARTMAP *artmap, int aclass, int mapclass)
{
  register int k, *wj ;

  if (Gdiag & DIAG_WRITE)
    printf("mapLearn(%d-->%d)\n", aclass, mapclass) ;

  for (k = 0 ; k < artmap->noutputs ; k++)
  {
    wj = Mwj(artmap, aclass, k) ;
    if (k != mapclass)
    {
      if (Gdiag & DIAG_WRITE)
        printf("deleting map connection %d --> %d\n", aclass, k) ;
      *wj = 0 ;
    }
  }
}

