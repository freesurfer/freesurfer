#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "error.h"
#include "matrix.h"
#include "rec.h"

REC *
RecRead(char *fname, int iop_neeg, int iop_nmeg)
{
  int   i,j,tnchan;
  float f;
  FILE  *fp;
  REC   *rec ;

  printf("read_rec(%s)\n",fname);

  
  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(NULL,
                (ERROR_BADFILE, "RecRead: could not open file %s",fname)) ;

  rec = (REC *)calloc(1, sizeof(REC)) ;
  if (rec==NULL) 
    ErrorExit(ERROR_NOMEMORY, "RecRead: couldn't allocate rec struct") ;

  fscanf(fp,"%*s");
  fscanf(fp,"%d %d %d",&rec->ntimepts,&rec->nmeg_channels,&rec->neeg_channels);
  tnchan = rec->nmeg_channels+rec->neeg_channels ;
  rec->ptime = (int)(2*pow(2.0,ceil(log((float)rec->ntimepts)/log(2.0))));
  printf("ntimepts=%d, ptime=%d\n",rec->ntimepts, rec->ptime);
#if 0
  if (sol_ntime>0 && sol_ntime!=rec->ntimepts)
  {
    printf("ntime does not match rec->ntimepts (%d != %d)\n",sol_ntime,rec->ntimepts);
    exit(0);
  }

  sol_ntime = rec->ntimepts; sol_ptime = tptime;

  if (rec->nmeg_channels>0 && rec->nmeg_channels!=tnmeg)
  {
    printf("nmeg does not match tnmeg (%d != %d)\n",rec->nmeg_channels,tnmeg);
    exit(0);
  }
  if (sol_neeg>0 && sol_neeg!=tneeg)
  {
    printf("neeg does not match tnmeg (%d != %d)\n",sol_neeg,tneeg);
    exit(0);
  }
#endif

  rec->latencies = (float *)calloc(rec->ptime, sizeof(float));
  if (!rec->latencies)
    ErrorExit(ERROR_NOMEMORY, "RecRead(%s): could not allocate latency vector",
              fname) ;
  rec->m_data = MatrixAlloc(tnchan,rec->ptime, MATRIX_REAL);
  for (j=0;j<rec->ntimepts;j++)
  {
    fscanf(fp,"%f",&f);
    rec->latencies[j] = f;
    for (i=0;i<rec->neeg_channels;i++)
    {
      fscanf(fp,"%f",&f);
      if (iop_neeg > 0)
        *MATRIX_RELT(rec->m_data,i+1,j+1) = f;
    }
    for (i=0;i<rec->nmeg_channels;i++)
    {
      fscanf(fp,"%f",&f);
      if (iop_nmeg > 0)
        *MATRIX_RELT(rec->m_data,i+iop_neeg+1,j+1) = f;
    }
  }
  fclose(fp);
  printf("rec file read, sample period %2.4f, starting latency %2.4f\n",
         rec->latencies[1]-rec->latencies[0], rec->latencies[0]);

#if 0  
  sol_dipcmp_val[sol_nrec] = matrix(sol_nnz*sol_nperdip,sol_ntime);
#endif
  return(rec) ;
}

