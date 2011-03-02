/**
 * @file  box.h
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


#ifndef BOX_H
#define BOX_H

typedef struct
{
  int  x0 ;
  int  y0 ;
  int  z0 ;
  int  x1 ;
  int  y1 ;
  int  z1 ;
  int  width ;
  int  height ;
  int  depth ;
}
BOX ;


int   BoxPrint(BOX *box, FILE *fp) ;
int   BoxExpand(BOX *box_src, BOX *box_dst, int dx, int dy, int dz) ;

#endif
