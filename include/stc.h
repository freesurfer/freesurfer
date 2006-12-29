/**
 * @file  stc.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
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


#ifndef STC_H
#define STC_H

#include "matrix.h"

typedef struct
{
  float   epoch_begin_lat ;
  float   sample_period ;
  int     ntimepts ;
  int     nperdip;
  int     ndipoles;
  int     nvertices ;
  int     *vertices ;
  MATRIX  *m_vals ;
}
STC ;

typedef struct
{
  float   epoch_begin_lat ;
  float   sample_period ;
  int     ntimepts ;
  int     nperdip;
  int     ndipoles;
  int     nvertices ;
  int     *vertices ;
  FILE*   file_handle;
}
STC_FILE;

typedef struct
{
  float   epoch_begin_lat ;
  float   sample_period ;
  int     nperdip;
  int     ndipoles;
  int     nvertices;
  MATRIX  *m_vals ;
}
STC_FRAME;

STC *StcRead(char *fname) ;
STC_FILE *StcOpen(char* fname);
void StcClose(STC_FILE* stc_file);
STC_FRAME *StcReadFrame(int fno,STC_FILE* stc_file);
int StcWriteFrame(STC_FILE* stc_file,STC_FRAME* cframe);
int StcWrite(char *fname, MATRIX *m_data, float epoch_begin_lat,
             float sample_period, int *vertices, int nvertices) ;
void StcFree(STC* stc);
#endif
