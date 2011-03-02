/**
 * @file  art1.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:42 $
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
  @(#)art1.c  1.1
  3/28/94
*/
/*----------------------------------------------------------------------

      File Name:   art1.c

      Description:

  $Header: /space/repo/1/dev/dev/utils/art1.c,v 1.3 2011/03/02 00:04:42 nicks Exp $
  $Log: art1.c,v $
  Revision 1.3  2011/03/02 00:04:42  nicks
  ENH: new license header mods

  Revision 1.2  2006/12/29 01:49:30  nicks
  added license header; ran astyle to set to kr and ansi code styling

  Revision 1.1  1997/03/18 18:26:29  fischl
  Initial revision

----------------------------------------------------------------------*/

/*-----------------------------------------------------------------
              INCLUDE FILES
-----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "machine.h"
#include "art1.h"
#include "diag.h"
#include "error.h"
#include "proto.h"

#define InsCalloc        calloc
#define InsHalloc(a,b)   calloc((int)(a), b)
#define InsFree          free
#define InsHfree         free

/*-----------------------------------------------------------------
              MACROS AND CONSTANTS
-----------------------------------------------------------------*/

#define BETA           .001
#define DELTA_RHO      .0001  /* amount above Tj to set row after mismatch */
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

void    artInitWeights(ART1 *art1) ;
int     artFeedBack(ART1 *art1, int class) ;
int     artFeedForward(ART1 *art1) ;
double  artChoice(ART1 *art, int j) ;
double  norm(double huge *a, int len) ;
void    intersect(double huge *a, double huge *b, double huge *c, int len) ;
void    artFastLearn(ART1 *art1, int j) ;

/*-----------------------------------------------------------------
                   MACROS
-----------------------------------------------------------------*/

#define Mzj(art1, i, j)   (art1->zj + \
                              (((long)j) * (long)art1->ninputs) + (long)i)

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
ART1 *
Art1Alloc(int ninputs, int max_f2, double rho)
{
  ART1  *art1 ;

  art1 = (ART1 *)InsCalloc(1, sizeof(ART1)) ;

  art1->ninputs = ninputs ;
  art1->noutputs = 0 ;        /* no nodes committed yet */
  art1->max_f2 = max_f2 ;
  art1->beta = BETA ;
  art1->rho = rho ;
  art1->scratch = (double *)InsCalloc(ninputs, sizeof(double)) ;
  art1->f0 = (double *)InsCalloc(ninputs, sizeof(double)) ;
  art1->f1 = (double *)InsCalloc(ninputs, sizeof(double)) ;
  art1->f2 = (double *)InsCalloc(max_f2, sizeof(double)) ;
  art1->flags = (int *)InsCalloc(max_f2, sizeof(int)) ;
  art1->zj =
    (double huge *)InsHalloc((long)ninputs * (long)max_f2, sizeof(double)) ;

  artInitWeights(art1) ;

  return(art1) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
ART1 *
Art1Read(char *fname)
{
  ART1  *art1 ;
  int     ninputs, noutputs, max_f2, i, j ;
  double  zj, rho, beta ;
  FILE    *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    return(NULL) ;

  fscanf(fp, "%d %d %d %lf %lf\n",
         &ninputs, &noutputs, &max_f2, &beta, &rho) ;

  art1 = Art1Alloc(ninputs, max_f2, rho) ;
  art1->noutputs = art1->noutputs ;
  for (j = 0 ; j < art1->noutputs ; j++)
    art1->flags[j] |= ART1_COMMITTED ;

  for (i = 0 ; i < art1->ninputs ; i++)
  {
    for (j = 0 ; j < art1->noutputs ; j++)
    {
      fscanf(fp, "%lf\n", &zj) ;
      *Mzj(art1, i, j) = zj ;
    }
  }

  fclose(fp) ;
  return(art1) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
Art1Write(ART1 *art1, char *fname)
{
  FILE  *fp ;
  int    i, j ;
  double zj ;

  fp = fopen(fname, "w") ;
  if (!fp)
    return(-1) ;

  fprintf(fp, "%d %d %d %f %f\n",
          art1->ninputs, art1->noutputs, art1->max_f2, art1->beta, art1->rho) ;

  for (i = 0 ; i < art1->ninputs ; i++)
  {
    for (j = 0 ; j < art1->noutputs ; j++)
    {
      zj = *Mzj(art1, i, j) ;
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
Art1Free(ART1 **art1)
{
  InsFree((*art1)->scratch) ;
  InsFree((*art1)->f0) ;
  InsFree((*art1)->f1) ;
  InsFree((*art1)->f2) ;
  InsFree((*art1)->flags) ;
  InsHfree((*art1)->zj) ;
  InsFree(*art1) ;
  *art1 = NULL ;
  return(0) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
Art1Process(ART1 *art1, double *I)
{
  int     nclass, i, class ;

  if (Gdiag & DIAG_WRITE) printf("Art1Process()\n") ;

  for (i = 0 ; i < art1->ninputs ; i++)
    art1->f0[i] = art1->f1[i] = I[i] ;

  for (i = 0 ; i < art1->noutputs ; i++)
    art1->flags[i] &= ~ART1_RESET ;

  nclass = 0 ;
  do
  {
    art1->class = class = artFeedForward(art1) ;
    if (class < 0) break ;    /* no classes left */
  }
  while (!artFeedBack(art1, class)) ;


  if (class < 0)     /* no prior match was enough to establish resonance */
  {
    if (art1->noutputs >= art1->max_f2)   /* no more nodes left */
    {
      /* find the closest feed-forward match and force it's use */
      for (i = 0 ; i < art1->noutputs ; i++)
        art1->flags[i] &= ~ART1_RESET ;
      class = artFeedForward(art1) ;
    }
    else    /* allocate new f2 node for this category */
    {
      class = art1->class = art1->noutputs++ ;
      art1->flags[class] = ART1_COMMITTED ;
    }
  }

  artFastLearn(art1, art1->class) ;

  if (Gdiag & DIAG_WRITE) printf("Art1Process() --> %d\n", class) ;
  return(class) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
void
artInitWeights(ART1 *art1)
{
  int     i, j ;
  double  *zj ;

  for (i = 0 ; i < art1->ninputs ; i++)
  {
    for (j = 0 ; j < art1->max_f2 ; j++)
    {
      zj = Mzj(art1, i, j) ;
      *zj = 1.0 ;
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
int
artFeedBack(ART1 *art1, int class)
{
  int    j, ninputs ;
  double match, norm_I_int_zj, norm_I ;

  j = class ;
  ninputs = art1->ninputs ;
  intersect(art1->f0, Mzj(art1, 0, j), art1->scratch, ninputs) ;
  norm_I_int_zj = norm(art1->scratch, ninputs) ;
  norm_I = norm(art1->f0, ninputs) ;
  if (!norm_I)
    match = 1 ;
  else
    match = norm_I_int_zj / norm_I ;

  if (Gdiag & DIAG_WRITE)
    printf("artFeedBack(%d) --> match %2.3f, n(I) %2.3f, "
           "n(I^Zj) %2.3f, rho %2.3f %s\n",
           class, match, norm_I, norm_I_int_zj, art1->rho,
           match >= art1->rho ? "RESONANCE" : "MISMATCH") ;

  if (match < art1->rho)   /* mismatch --> do reset */
  {
    art1->flags[j] |= ART1_RESET ;
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
int
artFeedForward(ART1 *art1)
{
  int       max_j, j ;
  double  f2, max_out ;

  max_j = -1 ;
  max_out = -1.0 ;

  for (j = 0 ; j < art1->noutputs ; j++)
  {
    if (art1->flags[j] & ART1_RESET)
      continue ;      /* this node was previous winner */

    f2 = art1->f2[j] = artChoice(art1, j) ;
    if (art1->f2[j] > max_out)
    {
      max_j = j ;
      max_out = art1->f2[j] ;
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
double
artChoice(ART1 *art1, int j)
{
  double  Tj, *zj, norm_zj, norm_I_int_zj ;

  zj = Mzj(art1, 0, j) ;    /* address of jth nodes weights */

  /* calculate I ^ zj */
  if (art1->flags[j] & ART1_COMMITTED)
  {
    intersect(art1->f0, zj, art1->scratch, art1->ninputs) ;
    norm_I_int_zj = norm(art1->scratch, art1->ninputs) ;
    norm_zj = norm(zj, art1->ninputs) ;

    Tj = (double)norm_I_int_zj / (art1->beta + (double)norm_zj) ;
  }
  else   /* Tj = |I| * ALPHA(j) */
  {
    Tj = (double)norm(art1->f0, art1->ninputs) ;
    Tj *= ALPHA(art1, j) ;
  }
  return(Tj) ;
}
/*
   take the intersection of vectors a and b, and return it
   in vector c.  All vectors have length 'len'.
*/
void
intersect(double huge *a, double huge *b, double huge *c, int len)
{
  register int i ;

  for (i = 0 ; i < len ; i++, a++, b++, c++)
    *c = MIN(*a, *b) ;
}

/*
   take the norm (city block) of vector a of length 'len'.
*/
double
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
void
artFastLearn(ART1 *art1, int j)
{
  double  norm_I_int_zj ;
  int     i ;

  if (j < 0) return ;   /* not a valid class */

  intersect(art1->f0, Mzj(art1, 0, j), art1->scratch, art1->ninputs) ;
  norm_I_int_zj = norm(art1->scratch, art1->ninputs) ;
  for (i = 0 ; i < art1->ninputs ; i++)
  {
    *Mzj(art1, i, j) = art1->scratch[i] ;
  }
}
int
Art1SetParms(ART1 *art1, double rho)
{
  art1->rho = rho ;
  return(NO_ERROR) ;
}
