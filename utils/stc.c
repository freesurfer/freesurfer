#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error.h"
#include "stc.h"
#include "matrix.h"
#include "fio.h"
#include "const.h"
#include "proto.h"

int
StcWrite(char *fname, MATRIX *m_data, float epoch_begin_lat,
             float sample_period, int *vertices, int nvertices)
{
  int   vno, ntime, row, col ;
  FILE  *fp;
  float val;

  ntime = m_data->cols ;

#if 0
  sprintf(fname,"%s.stc",fstem);
#endif
  fp = fopen(fname,"w");
  if (fp==NULL)
  ErrorReturn(ERROR_NOFILE,
              (ERROR_NOFILE, "MRISwriteStc: can't open file %s\n",fname)) ;


  fwriteFloat(epoch_begin_lat,fp);
  fwriteFloat(sample_period,fp);
  fwriteInt(nvertices,fp);
  for (vno=0;vno<nvertices;vno++)
    fwriteInt(vertices[vno],fp);

  fwriteInt(ntime,fp);
  for (col = 1 ; col <= m_data->cols ; col++)
  {
    for (row = 1 ; row <= m_data->rows ; row++)
    {
      /* row is dipole, and col is time */
      val = *MATRIX_RELT(m_data, row, col) ;
      fwriteFloat(val,fp);
    }
  }
  fclose(fp);
  printf("sol timecourse file %s written\n",fname);
  return(NO_ERROR) ;
}

STC *
StcRead(char *fstem)
{
  int   j,vno, ndipoles, dipoles_per_location ;
  char  fname[STRLEN];
  FILE  *fp;
  STC   *stc ;
  long  here, there ;

  stc = (STC *)calloc(1, sizeof(STC)) ;

  sprintf(fname,"%s.stc",fstem);
  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorReturn(NULL, 
                (ERROR_NOFILE, "StcRead: could not open %s", fname)) ;

  stc->epoch_begin_lat = freadFloat(fp);
  stc->sample_period = freadFloat(fp);
  stc->nvertices = freadInt(fp);
  printf("epoch_begin_lat=%f, sample_period=%f, nvertices=%d\n",
          stc->epoch_begin_lat,stc->sample_period,stc->nvertices);
  stc->vertices = (int *)calloc(stc->nvertices, sizeof(int)) ;
  if (!stc->vertices)
    ErrorExit(ERROR_NOMEMORY, "StcRead(%s): could not allocated %d vector",
              fname, stc->nvertices) ;

  for (vno=0;vno<stc->nvertices;vno++)
  {
    stc->vertices[vno] = freadInt(fp);
#if 0
    if (stc->vertices[vno]>=vertex_index)
    {
      printf("### stc->vertices[%d] = %d out of bounds\n",vno,stc->vertices[vno]);
      exit(1);
    }
#endif
  }
  stc->ntimepts = freadInt(fp);
  printf("ntime=%d\n",stc->ntimepts);
  here = ftell(fp) ;
  fseek(fp, 0L, SEEK_END) ;
  there = ftell(fp) ;
  fseek(fp, here, SEEK_SET) ;

  ndipoles = nint((there-here) / (sizeof(float)) / stc->ntimepts) ;
  dipoles_per_location = nint(ndipoles / (float)stc->nvertices) ;
  fprintf(stderr, "%d dipoles detected - %d per location\n",
          ndipoles, dipoles_per_location) ;

  stc->m_vals = 
    MatrixAlloc(dipoles_per_location*stc->nvertices,stc->ntimepts,MATRIX_REAL);
  if (!stc->m_vals)
    ErrorExit(ERROR_NOMEMORY, "StcRead(%s) could not allocate %dx%d matrix",
              fname, stc->nvertices, stc->ntimepts) ;

  for (j=0;j<stc->ntimepts;j++)
  {
    for (vno=0;vno<stc->m_vals->rows;vno++)
    {
      *MATRIX_RELT(stc->m_vals, vno+1,j+1) = freadFloat(fp);
    }
  }
  fclose(fp);
  printf("soltimecourse file %s read\n",fname);
  return(stc) ;
}
