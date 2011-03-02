/**
 * @file  stc.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.5 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

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

STC_FILE *StcOpen(char* fstem)
{
  STC_FILE* stc_file;
  char fname[STRLEN];
  int vno;

  stc_file = (STC_FILE *)calloc(1, sizeof(STC_FILE));

  sprintf(fname,"%s.stc",fstem);
  stc_file->file_handle = fopen(fname,"r");
  if (stc_file->file_handle==NULL)
    ErrorReturn(NULL,
                (ERROR_NOFILE, "StcRead: could not open %s", fname)) ;

  stc_file->epoch_begin_lat = freadFloat(stc_file->file_handle);
  stc_file->sample_period = freadFloat(stc_file->file_handle);
  stc_file->nvertices = freadInt(stc_file->file_handle);
  stc_file->vertices = (int *)calloc(stc_file->nvertices, sizeof(int)) ;
  if (!stc_file->vertices)
    ErrorExit(ERROR_NOMEMORY, "StcRead(%s): could not allocated %d vector",
              fname, stc_file->nvertices) ;

  for (vno=0;vno<stc_file->nvertices;vno++)
  {
    stc_file->vertices[vno] = freadInt(stc_file->file_handle);
  }
  stc_file->ntimepts = freadInt(stc_file->file_handle);

  return stc_file;
}

void StcClose(STC_FILE *stc_file)
{
  if (stc_file->file_handle)
    fclose(stc_file->file_handle);

  if (stc_file->vertices)
  {
    free(stc_file->vertices);
    stc_file->vertices = 0;
    stc_file->nvertices = 0;
  }
}

STC_FRAME *StcReadFrame(int fno, STC_FILE* stc_file)
{
  STC_FRAME *stc;
  int baseoffset,offset;
  int framesize,vno;

  assert(stc_file->file_handle);

  stc = (STC_FRAME *)calloc(1, sizeof(STC_FRAME));

  stc->nperdip = stc_file->nperdip;
  stc->ndipoles = stc_file->ndipoles;
  stc->nvertices = stc_file->nvertices;

  baseoffset = 12 + 4*stc_file->nvertices + 4;
  framesize = 4*stc_file->ndipoles*stc_file->nperdip;
  offset = baseoffset + fno*framesize;
  fseek(stc_file->file_handle,offset,SEEK_SET);

  stc->m_vals =
    MatrixAlloc(stc->nperdip*stc->nvertices,1,MATRIX_REAL);

  if (!stc->m_vals)
    ErrorExit(ERROR_NOMEMORY, "StcReadFrame could not allocate %dx%d matrix",
              stc_file->nperdip*stc_file->nvertices, 1) ;

  for (vno=0;vno<stc->m_vals->rows;vno++)
  {
    *MATRIX_RELT(stc->m_vals, vno+1,1) = freadFloat(stc_file->file_handle);
  }

  return stc;
}

STC *
StcRead(char *fstem)
{
  int   j,vno;
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
  printf("StcRead: epoch_begin_lat=%f, sample_period=%f, nvertices=%d\n",
         stc->epoch_begin_lat,stc->sample_period,stc->nvertices);
  stc->vertices = (int *)calloc(stc->nvertices, sizeof(int)) ;
  if (!stc->vertices)
    ErrorExit(ERROR_NOMEMORY, "StcRead(%s): could not allocated %d vector",
              fname, stc->nvertices) ;

  for (vno=0;vno<stc->nvertices;vno++)
  {
    stc->vertices[vno] = freadInt(fp);
  }
  stc->ntimepts = freadInt(fp);
  printf("ntime=%d\n",stc->ntimepts);
  here = ftell(fp) ;
  fseek(fp, 0L, SEEK_END) ;
  there = ftell(fp) ;
  fseek(fp, here, SEEK_SET) ;

  stc->ndipoles = nint((there-here) / (sizeof(float)) / stc->ntimepts) ;
  stc->nperdip = nint(stc->ndipoles / (float)stc->nvertices) ;
  fprintf(stderr, "%d dipoles detected - %d per location\n",
          stc->ndipoles, stc->nperdip);

  stc->m_vals =
    MatrixAlloc(stc->nperdip*stc->nvertices,stc->ntimepts,MATRIX_REAL);
  if (!stc->m_vals)
    ErrorExit(ERROR_NOMEMORY, "StcRead(%s) could not allocate %dx%d matrix",
              fname, stc->nperdip*stc->nvertices, stc->ntimepts) ;

  for (j=0;j<stc->ntimepts;j++)
  {
    for (vno=0;vno<stc->m_vals->rows;vno++)
    {
      *MATRIX_RELT(stc->m_vals, vno+1,j+1) = freadFloat(fp);
    }
  }
  fclose(fp);
  /* printf("soltimecourse file %s read\n",fname); */
  return(stc) ;
}

void StcFree(STC* stc)
{
  if (stc->vertices)
    free(stc->vertices);
  if (stc->m_vals)
    MatrixFree(&stc->m_vals);

  free(stc);
}
