#ifndef STC_H
#define STC_H

#include "matrix.h"

typedef struct
{
  float   epoch_begin_lat ;
  float   sample_period ;
  int     ntimepts ;
  int     nvertices ;
  int     *vertices ;
  MATRIX  *m_vals ;
} STC ;

STC *StcRead(char *fname) ;
int StcWrite(char *fname, MATRIX *m_data, float epoch_begin_lat,
             float sample_period, int *vertices, int nvertices) ;

#endif
