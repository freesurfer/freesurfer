/*
  @(#)art1.h  1.3
  8/10/95
*/
/*----------------------------------------------------------------------

      File Name:   art1.h

      Description:  

  $Header: /space/repo/1/dev/dev/include/art1.h,v 1.1 1997/03/18 18:17:34 fischl Exp $
  $Log: art1.h,v $
  Revision 1.1  1997/03/18 18:17:34  fischl
  Initial revision

----------------------------------------------------------------------*/

#ifndef ART1_H
#define ART1_H


#include "machine.h"

typedef struct
{
  int     ninputs ;       /* # of inputs */
  int     noutputs ;      /* # of outputs committed so far */
  int     max_f2 ;        /* size of f2 */
  double  beta ;          /* for weight initialization */
  double  rho ;           /* vigilance */
  int     class ;         /* class of previous recognition */
  double  *scratch ;      /* for scratch calculations of intersection */
  double  *f0 ;           /* points to inputs (provided by caller) */
  double  *f1 ;           /* input vector */
  double  *f2 ; 
  double  huge *zj ;      /* top down weights */
  int     *flags ;        /* is node committed, reset */
} ART1 ;

#define ART1_RESET      0x0001
#define ART1_COMMITTED  0x0002

ART1   *Art1Alloc(int ninputs, int noutputs, double rho) ;
ART1   *Art1Read(char *fname) ;
int    Art1Write(ART1 *art1, char *fname) ;
int    Art1Free(ART1 **art1) ;
int    Art1Process(ART1 *art, double *I) ;
int    Art1Reset(ART1 *art1) ;
int    Art1SetParms(ART1 *art, double rho) ;


#endif
