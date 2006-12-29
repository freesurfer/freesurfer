/**
 * @file  pbm.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.2 $
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


/* pbm.h - header file for libpbm portable bitmap library
*/

#ifndef _PBM_H_
#define _PBM_H_

#include "pbmplus.h"

typedef unsigned char bit;
#define PBM_WHITE 0
#define PBM_BLACK 1


/* Magic constants. */

#define PBM_MAGIC1 'P'
#define PBM_MAGIC2 '1'
#define RPBM_MAGIC2 '4'
#define PBM_FORMAT (PBM_MAGIC1 * 256 + PBM_MAGIC2)
#define RPBM_FORMAT (PBM_MAGIC1 * 256 + RPBM_MAGIC2)
#define PBM_TYPE PBM_FORMAT


/* Macro for turning a format number into a type number. */

#define PBM_FORMAT_TYPE(f) ((f) == PBM_FORMAT || (f) == RPBM_FORMAT ? PBM_TYPE : -1)


/* Declarations of routines. */

void pbm_init ARGS(( int* argcP, char* argv[] ));

#define pbm_allocarray( cols, rows ) ((bit**) pm_allocarray( cols, rows, sizeof(bit) ))
#define pbm_allocrow( cols ) ((bit*) pm_allocrow( cols, sizeof(bit) ))
#define pbm_freearray( bits, rows ) pm_freearray( (char**) bits, rows )
#define pbm_freerow( bitrow ) pm_freerow( (char*) bitrow )

bit** pbm_readpbm ARGS(( FILE* file, int* colsP, int* rowsP ));
void pbm_readpbminit ARGS(( FILE* file, int* colsP, int* rowsP, int* formatP ));
void pbm_readpbmrow ARGS(( FILE* file, bit* bitrow, int cols, int format ));
char* pm_read_unknown_size ARGS(( FILE* file, long* buf ));

void pbm_writepbm ARGS(( FILE* file, bit** bits, int cols, int rows, int forceplain ));
void pbm_writepbminit ARGS(( FILE* file, int cols, int rows, int forceplain ));
void pbm_writepbmrow ARGS(( FILE* file, bit* bitrow, int cols, int forceplain ));

#endif /*_PBM_H_*/
