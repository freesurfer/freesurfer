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
} STC ;

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
} STC_FILE;

typedef struct 
{
  float   epoch_begin_lat ;
  float   sample_period ;
  int     nperdip;
  int     ndipoles;
  int     nvertices;
  MATRIX  *m_vals ;
} STC_FRAME;

STC *StcRead(char *fname) ;
STC_FILE *StcOpen(char* fname);
void StcClose(STC_FILE* stc_file);
STC_FRAME *StcReadFrame(int fno,STC_FILE* stc_file);
int StcWriteFrame(STC_FILE* stc_file,STC_FRAME* cframe);
int StcWrite(char *fname, MATRIX *m_data, float epoch_begin_lat,
             float sample_period, int *vertices, int nvertices) ;
void StcFree(STC* stc);
#endif
