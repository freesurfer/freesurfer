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
init_header(IMAGE *I,char *onm,char *snm,int nfr,char *odt,int rw,int cl,int pfmt,int nc,char *desc)
{
  int bytes ;

  I->num_frame = nfr ;
  I->rows = rw ;
  I->cols = cl ;
  I->pixel_format = pfmt ;
  bytes = rw*cl*nfr ;
  switch (pfmt)
  {
  default:
  case PFBYTE:
    I->sizepix = sizeof(char) ;
    break ;
  case PFFLOAT:
    bytes *= sizeof(float) ;
    I->sizepix = sizeof(float) ;
    break ;
  case PFDOUBLE:
    bytes *= sizeof(double) ;
    I->sizepix = sizeof(double) ;
    break ;
  case PFINT:
    bytes *= sizeof(int) ;
    I->sizepix = sizeof(int) ;
    break ;
  case PFSHORT:
    bytes *= sizeof(short) ;
    I->sizepix = sizeof(short) ;
    break ;
  }
  I->numpix = I->rows * I->cols ;
  I->image = (char *)calloc(bytes, sizeof(char)) ;
  if (!I->image)
    ErrorExit(ERROR_NOMEMORY, "init_header: could not allocate %d bytes",
              bytes) ;

  return(NO_ERROR) ;
}
int
h_copy(IMAGE *Isrc, IMAGE *Idst)
{
  int   bytes ;

  bytes = Isrc->numpix * Isrc->sizepix ;
  memmove(Idst->image, Isrc->image, bytes) ;
  return(NO_ERROR) ;
}
