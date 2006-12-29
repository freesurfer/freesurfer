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
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
