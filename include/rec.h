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

/* flag == 1 - load only EEG channels
   flag == 2 - load only MEG channels
*/  
REC *RecReadPartially(char *fname, int iop_neeg, int iop_nmeg, int flag) ;
#endif
