/**
 * @file  artmap.h
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


/*----------------------------------------------------------------------

      File Name:   artmap.h

      Description:

  $Header: /space/repo/1/dev/dev/include/artmap.h,v 1.3 2011/03/02 00:04:09 nicks Exp $
  $Log: artmap.h,v $
  Revision 1.3  2011/03/02 00:04:09  nicks
  ENH: new license header mods

  Revision 1.2  2006/12/29 02:08:59  nicks
  added license header; ran astyle to set to kr and ansi code styling

  Revision 1.1  1997/03/18 18:17:39  fischl
  Initial revision

----------------------------------------------------------------------*/

#ifndef ARTMAP_H
#define ARTMAP_H


#include "machine.h"

typedef struct
{
  int     ninputs ;       /* # of inputs */
  int     noutputs ;      /* # of outputs (size of map field) */
  int     max_f2 ;        /* size of f2 */
  int     f2nodes ;       /* actual # of f2 nodes committed */
  int     ncommitted ;    /* # of committed nodes */
  int     learn ;         /* in learn mode */
  double  beta ;          /* for weight initialization */
  double  rho ;           /* vigilance */
  double  rho_bar;        /* baseline vigilance */
  int     class ;         /* class of previous recognition */
  double  *scratch ;      /* for scratch calculations of intersection */
  double  *f0 ;           /* points to inputs (provided by caller) */
  double  *f1 ;           /* input vector */
  double  *f2 ;
  double  huge *zj ;      /* top down weights */
  int     *flags ;        /* is node committed, reset */
  int     huge *w ;       /* f2->map weights */
  double  match ;        /* for use in match tracking */
}
ARTMAP ;

#define ARTMAP_RESET      0x0001
#define ARTMAP_COMMITTED  0x0002

ARTMAP   *ArtmapAlloc(int ninputs, int noutputs, double rho_bar, int max_f2) ;
ARTMAP   *ArtmapRead(char *fname) ;
int      ArtmapWrite(ARTMAP *artmap, char *fname) ;
int      ArtmapFree(ARTMAP **artmap) ;
int      ArtmapProcess(ARTMAP *artmap, double *I) ;
int      ArtmapLearn(ARTMAP *artmap, double *I, int class) ;


#endif
