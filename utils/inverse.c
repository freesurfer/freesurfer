/**
 * @file  inverse.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.4 $
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
#include <string.h>

#include "inverse.h"
#include "error.h"
#include "matrix.h"
#include "const.h"
#include "macros.h"

IOP *
IOPRead(char *fname, int hemi)
{
  IOP   *iop ;
  int   i,j,k,jc,d, dipoles_in_decimation ;
  FILE  *fp;
  char  c,str[STRLEN];
  float f;
  int   z ;

  printf("read_iop(%s,%d)\n",fname,hemi);

  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorReturn(NULL,
                (ERROR_NOFILE, "IOPRead: can't open file %s\n",fname));
  iop = calloc(1, sizeof(IOP)) ;
  if (!iop)
    ErrorReturn(NULL,
                (ERROR_NOMEMORY, "IOPRead: can't allocate struct\n"));
  iop->pthresh = 1000;
  c = getc(fp);
  if (c=='#')
  {
    fscanf(fp,"%s",str);
    if (!strcmp(str,"version"))
      fscanf(fp,"%d",&iop->version);
    printf("iop version = %d\n",iop->version);
    fscanf(fp,"%d %d %d %d",
           &iop->neeg_channels,
           &iop->nmeg_channels,
           &iop->ndipoles_per_location,
           &iop->ndipole_files);

    iop->nchan = iop->neeg_channels+iop->nmeg_channels;
    for (i=1;i<=hemi;i++)
    {
      fscanf(fp,"%d %d",&dipoles_in_decimation,&iop->ndipoles);
      if (i==hemi)
      {
#if 0
        if (iop->ndipoles_per_location != sol_ndec)
        {
          fclose(fp);
          IOPFree(&iop) ;
          ErrorReturn(NULL,
                      (ERROR_BADFILE,
                       "IOPRead: .dec and .iop file mismatch (%d!=%d)\n",
                       sol_ndec,iop->ndipoles_per_location));
        }
#endif
        iop->m_iop =
          MatrixAlloc(iop->ndipoles*iop->ndipoles_per_location,iop->nchan,
                      MATRIX_REAL);
        if (iop->version==1)
          iop->m_forward =
            MatrixAlloc(iop->nchan,iop->ndipoles*iop->ndipoles_per_location,
                        MATRIX_REAL);

#if 0
        sol_M = matrix(iop->nchan,iop->nchan); /* temporary space for xtalk */
        sol_Mi = matrix(iop->nchan,iop->nchan);
        sol_sensvec1 = vector(iop->nchan);
        sol_sensvec2 = vector(iop->nchan);
        sol_sensval = vector(iop->nchan);
#endif

        iop->dipole_normalization = calloc(iop->ndipoles, sizeof(float));
        iop->dipole_vertices = calloc(iop->ndipoles, sizeof(int));
        if (!iop->dipole_vertices)
          ErrorReturn(NULL,
                      (ERROR_NOMEMORY,
                       "IOPRead: could not allocated %d v indices",
                       iop->ndipoles)) ;
        iop->pvals = VectorAlloc(iop->ndipoles, MATRIX_REAL);
        iop->spatial_priors = VectorAlloc(iop->ndipoles, MATRIX_REAL);


        iop->bad_sensors = calloc(iop->nchan, sizeof(int));
        if (!iop->bad_sensors)
          ErrorReturn(NULL,
                      (ERROR_NOMEMORY,
                       "IOPRead: could not allocate bad sensor array",
                       iop->nchan)) ;

        /* initialize bad sensor locations*/
        for (z=0;z<iop->nchan;z++)
          iop->bad_sensors[z] = 0;

      }
      for (j=0;j<iop->ndipoles;j++)
      {
        if (i==hemi)
        {
          fscanf(fp,"%d",&d);
          iop->dipole_vertices[j] = d;
        }
        else
          fscanf(fp,"%*d");
      }
      for (j=0;j<iop->ndipoles;j++)
      {
        if (i==hemi)
        {
          fscanf(fp,"%f",&f);
          *MATRIX_RELT(iop->pvals,j+1,1) = f;
          f = fabs(f);
          if (f<iop->pthresh)
            iop->pthresh = f;
        }
        else
          fscanf(fp,"%*f");
      }
      for (j=0;j<iop->ndipoles;j++)
      {
        if (i==hemi)
        {
          fscanf(fp,"%f",&f);
          *MATRIX_RELT(iop->spatial_priors,j+1,1) = f;
#if 0
          vertex[iop->dipole_vertices[j]].val = f;
#endif
        }
        else
          fscanf(fp,"%*f");
      }
      for (j=0;j<iop->ndipoles;j++)
      {
        for (jc=0;jc<iop->ndipoles_per_location;jc++)
        {
          for (k=0;k<iop->nchan;k++)
          {
            if (i==hemi)
            {
              fscanf(fp,"%f",&f);
              *MATRIX_RELT(iop->m_iop, j*iop->ndipoles_per_location+jc+1,k+1)
              = f;
            }
            else
              fscanf(fp,"%*f");
          }
        }
      }
      if (iop->version==1)
      {
        for (j=0;j<iop->ndipoles;j++)
        {
          for (jc=0;jc<iop->ndipoles_per_location;jc++)
          {
            for (k=0;k<iop->nchan;k++)
            {
              if (i==hemi)
              {
                fscanf(fp,"%f",&f);
                *MATRIX_RELT(iop->m_forward,
                             k+1,
                             j*iop->ndipoles_per_location+jc+1) = f;
              }
              else
                fscanf(fp,"%*f");
            }
          }
        }
      }
    }
  }
  else
  {
    printf("Can't read binary .iop files\n");
  }
  fclose(fp);
  printf("neeg_channels=%d, nmeg_channels=%d, iop->ndipoles_per_location=%d, "
         "iop->ndipole_files=%d\n",
         iop->neeg_channels, iop->nmeg_channels,iop->ndipoles_per_location,
         iop->ndipole_files);
  return(iop) ;
}
int
IOPWrite(IOP *iop, char *fname)
{
  return(NO_ERROR) ;
}
int
IOPFree(IOP **piop)
{
  IOP *iop ;

  iop = *piop ;
  *piop = NULL ;
  if (iop->dipole_vertices)
    free(iop->dipole_vertices) ;
  if (iop->spatial_priors)
    free(iop->spatial_priors) ;
  if (iop->pvals)
    free(iop->pvals) ;
  if (iop->m_iop)
    MatrixFree(&iop->m_iop) ;
  if (iop->m_forward)
    MatrixFree(&iop->m_forward) ;
  if (iop->bad_sensors)
    free(iop->bad_sensors) ;
  return(NO_ERROR) ;
}

int
IOPNormalize(IOP *iop)
{
  int j,k,jc;
  double sum, val;

  for (j=0;j<iop->ndipoles;j++)
  {
    sum = 0;
    for (jc=0;jc<iop->ndipoles_per_location;jc++)
      for (k=0;k<iop->nchan;k++)
      {
        val = *MATRIX_RELT(iop->m_iop,j*iop->ndipoles_per_location+jc+1,k+1);
        sum += SQR(val);
      }

    sum = sqrt(sum);
    if (!DZERO(sum))
    {
      for (jc=0;jc<iop->ndipoles_per_location;jc++)
        for (k=0;k<iop->nchan;k++)
          *MATRIX_RELT(iop->m_iop,j*iop->ndipoles_per_location+jc+1,k+1)/= sum;
    }
  }
  return(NO_ERROR) ;
}
MATRIX *
IOPapplyInverseOperator(IOP *iop, REC *rec, MATRIX *m_sol)
{
  m_sol = MatrixMultiply(iop->m_iop, rec->m_data, NULL) ;
  return(m_sol) ;
}

