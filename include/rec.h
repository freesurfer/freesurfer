#ifndef REC_H
#define REC_H

#include "matrix.h"

typedef struct
{
  int    ntimepts ;
  int    ptime ;
  int    neeg_channels ;
  int    nmeg_channels ;
  float  *latencies ;
  MATRIX *m_data ;
} REC ;


REC *RecRead(char *fname) ;  /* each column is a time course (except 1st)*/

#endif
