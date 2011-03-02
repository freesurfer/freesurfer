/**
 * @file  pnm.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
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


/* pnm.h - header file for libpnm portable anymap library
*/

#ifndef _PNM_H_
#define _PNM_H_

#include "ppm.h"
typedef pixel xel;
typedef pixval xelval;
#define PNM_MAXMAXVAL PPM_MAXMAXVAL
#define PNM_GET1(x) PPM_GETB(x)
#define PNM_ASSIGN1(x,v) PPM_ASSIGN(x,0,0,v)
#define PNM_EQUAL(x,y) PPM_EQUAL(x,y)
#define PNM_FORMAT_TYPE(f) PPM_FORMAT_TYPE(f)

/* Declarations of routines. */

void pnm_init ARGS(( int* argcP, char* argv[] ));

#define pnm_allocarray( cols, rows ) ((xel**) pm_allocarray( cols, rows, sizeof(xel) ))
#define pnm_allocrow( cols ) ((xel*) pm_allocrow( cols, sizeof(xel) ))
#define pnm_freearray( xels, rows ) pm_freearray( (char**) xels, rows )
#define pnm_freerow( xelrow ) pm_freerow( (char*) xelrow )

xel** pnm_readpnm ARGS(( FILE* file, int* colsP, int* rowsP, xelval* maxvalP, int* formatP ));
void pnm_readpnminit ARGS(( FILE* file, int* colsP, int* rowsP, xelval* maxvalP, int* formatP ));
void pnm_readpnmrow ARGS(( FILE* file, xel* xelrow, int cols, xelval maxval, int format ));

void pnm_writepnm ARGS(( FILE* file, xel** xels, int cols, int rows, xelval maxval, int format, int forceplain ));
void pnm_writepnminit ARGS(( FILE* file, int cols, int rows, xelval maxval, int format, int forceplain ));
void pnm_writepnmrow ARGS(( FILE* file, xel* xelrow, int cols, xelval maxval, int format, int forceplain ));

xel pnm_backgroundxel ARGS(( xel** xels, int cols, int rows, xelval maxval, int format ));
xel pnm_backgroundxelrow ARGS(( xel* xelrow, int cols, xelval maxval, int format ));
xel pnm_whitexel ARGS(( xelval maxval, int format ));
xel pnm_blackxel ARGS(( xelval maxval, int format ));
void pnm_invertxel ARGS(( xel* x, xelval maxval, int format ));
void pnm_promoteformat ARGS(( xel** xels, int cols, int rows, xelval maxval, int format, xelval newmaxval, int newformat ));
void pnm_promoteformatrow ARGS(( xel* xelrow, int cols, xelval maxval, int format, xelval newmaxval, int newformat ));

extern xelval pnm_pbmmaxval;
/* This is the maxval used when a PNM program reads a PBM file.  Normally
** it is 1; however, for some programs, a larger value gives better results
*/

#endif /*_PNM_H_*/
