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


/* each column is a time course (except 1st)*/
REC *RecRead(char *fname, int iop_neeg, int iop_nmeg) ;  

#endif
