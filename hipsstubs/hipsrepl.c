/**
 * @file  hipsrepl.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:58 $
 *    $Revision: 1.5 $
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


/*
 *       FILE NAME:   hipsrepl.c
 *
 *       DESCRIPTION: replacement routines for some hips functionality
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <fcntl.h>
#include <unistd.h>  /* for SEEK_ constants */

#include "hmem.h"
#include <hipl_format.h>

#include "hips.h"

#include "error.h"

#include "image.h"
int
init_header(IMAGE *I,char *onm,char *snm,int nfr,char *odt,int rw,int cl,int pfmt,int nc,char *desc) {
  int bytes ;

  I->num_frame = nfr ;
  I->orows = I->rows = rw ;
  I->ocols = I->cols = cl ;
  I->pixel_format = pfmt ;
  bytes = rw*cl*nfr ;
  switch (pfmt) {
  default:
  case PFBYTE:
    I->sizepix = sizeof(char) ;
    break ;
  case PFFLOAT:
    I->sizepix = sizeof(float) ;
    break ;
  case PFDOUBLE:
    I->sizepix = sizeof(double) ;
    break ;
  case PFINT:
    I->sizepix = sizeof(int) ;
    break ;
  case PFSHORT:
    I->sizepix = sizeof(short) ;
    break ;
  case PFRGB:
  case PFBGR:
    I->sizepix = 3*sizeof(byte);
    break;
  case PFRGBZ:
  case PFZRGB:
  case PFBGRZ:
  case PFZBGR:
    I->sizepix = 4*sizeof(byte);
    break;
  case PFSTEREO:
    I->sizepix = sizeof(byte);
    break;
  case PFINTPYR:
    I->sizepix = sizeof(int);
    break;
  case PFFLOATPYR:
    I->sizepix = sizeof(float);
    break;
  }
  bytes *= I->sizepix ;
  I->numpix = I->rows * I->cols ;
  I->sizeimage = I->numpix * I->sizepix ;
  I->firstpix = I->image ;
  I->image = (byte *)calloc(bytes, sizeof(char)) ;
  if (!I->image)
    ErrorExit(ERROR_NOMEMORY, "init_header: could not allocate %d bytes",
              bytes) ;

  return(NO_ERROR) ;
}
int
h_copy(IMAGE *Isrc, IMAGE *Idst) {
  int   bytes ;

  bytes = Isrc->numpix * Isrc->sizepix ;
  memmove(Idst->image, Isrc->image, bytes) ;
  return(NO_ERROR) ;
}
int
free_header(IMAGE *I) {
  if (I->image)
    free(I->image) ;
  free(I) ;
  return(0) ;
}

