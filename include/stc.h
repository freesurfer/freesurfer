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
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.6 $
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
