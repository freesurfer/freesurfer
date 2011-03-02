/**
 * @file  lpafile.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
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


#ifndef LPAFILE_H
#define LPAFILE_H


#define NPOINTS   4

typedef struct
{
  int  xp[NPOINTS] ;    /* x coordinates of corners */
  int  yp[NPOINTS] ;    /* y coordinates of corners */
  int  xc ;             /* x coordinate of centroid */
  int  yc ;             /* y coordinate of centroid */
  long fpos ;           /* position in answer file */
}
LP_BOX ;

typedef struct
{
  char    fname[100] ;     /* name of the file */
  FILE    *fp ;
  char    **filelist ;     /* array of pointers */
  int     nfiles ;         /* size of filelist */
  LP_BOX  *coords ;
  int     last_written ;   /* index of last written lp_box */
  int     current ;
  int     flush ;
}
LP_ANSWER_FILE, LPAF ;

LPAF *LPAFcreate(char *out_fname, int argc, char *argv[]) ;
int  LPAFwrite(LPAF *lpaf, int current) ;
int  LPAFread(LPAF *lpaf, int current) ;
int  LPAFset(LPAF *lpaf, int current, int *xp, int *yp, int xc, int yc) ;
int  LPAFwriteImageAnswer(LPAF *lpaf, int current) ;
int  LPAFresetImageAnswer(LPAF *lpaf, int current) ;
int  LPAFreadImageAnswer(LPAF *lpaf, int current) ;

#endif
