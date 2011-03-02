/**
 * @file  art1.h
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
  @(#)art1.h  1.3
  8/10/95
*/
/*----------------------------------------------------------------------

      File Name:   art1.h

      Description:

  $Header: /space/repo/1/dev/dev/include/art1.h,v 1.3 2011/03/02 00:04:09 nicks Exp $
  $Log: art1.h,v $
  Revision 1.3  2011/03/02 00:04:09  nicks
  ENH: new license header mods

  Revision 1.2  2006/12/29 02:08:59  nicks
  added license header; ran astyle to set to kr and ansi code styling

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
}
ART1 ;

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
