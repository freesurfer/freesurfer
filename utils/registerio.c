/**
 * @file  registerio.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:39 $
 *    $Revision: 1.15 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "resample.h"
#include "registerio.h"
#include "fio.h"
#include "macros.h" // FEQUAL
#include "proto.h" // nint


/* ----------------------------------------------------------
  Name: regio_read_register()
  Reads a registration file.

  subject -- name of subject as found in the data base
  inplaneres -- in-plane resolution
  betplaneres -- between-plane resolution
  intensity -- for the register program
  R - matrix to convert from xyz in COR space to xyz in Volume space,
      ie, xyzVol = R*xyzCOR
  float2int - if the regfile has a line after the matrix, the string
      is passed to float2int_code(), the result of which is passed
      back as float2int. If there is no extra line, FLT2INT_TKREG
      is returned (indicating that the regfile was created by
      tkregister).
  -------------------------------------------------------------*/
int regio_read_register(char *regfile, char **subject, float *inplaneres,
                        float *betplaneres, float *intensity,  MATRIX **R,
                        int *float2int)
{
  FILE *fp;
  char tmp[1000];
  int r,c,n;
  float val;

  fp = fopen(regfile,"r");
  if (fp==NULL)
  {
    perror("regio_read_register()");
    fprintf(stderr,"Could not open %s\n",regfile);
    return(1);
  }

  /* subject name */
  n = fscanf(fp,"%s",tmp);
  if (n != 1)
  {
    perror("regio_read_register()");
    fprintf(stderr,"Error reading subject from %s\n",regfile);
    fclose(fp);
    return(1);
  }
  *subject = (char *) calloc(strlen(tmp)+2,sizeof(char));
  sprintf(*subject,"%s",tmp);

  /* in-plane resolution */
  n = fscanf(fp,"%f",inplaneres);
  if (n != 1)
  {
    perror("regio_read_register()");
    fprintf(stderr,"Error reading inplaneres from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  /* between-plane resolution */
  n = fscanf(fp,"%f",betplaneres);
  if (n != 1)
  {
    perror("regio_read_register()");
    fprintf(stderr,"Error reading betplaneres from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  /* intensity*/
  n = fscanf(fp,"%f",intensity);
  if (n != 1)
  {
    perror("regio_read_register()");
    fprintf(stderr,"Error reading intensity from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  *R = MatrixAlloc(4,4,MATRIX_REAL);
  if (*R == NULL)
  {
    fprintf(stderr,"regio_read_register(): could not alloc R\n");
    fclose(fp);
    return(1);
  }

  /* registration matrix */
  for (r=0;r<4;r++)
  {
    for (c=0;c<4;c++)
    {
      n = fscanf(fp,"%f",&val);
      if (n != 1)
      {
        perror("regio_read_register()");
        fprintf(stderr,"Error reading R[%d][%d] from %s\n",r,c,regfile);
        fclose(fp);
        return(1);
      }
      (*R)->rptr[r+1][c+1] = val;
    }
  }

  /* Get the float2int method string */
  n = fscanf(fp,"%s",&tmp[0]);
  fclose(fp);

  if (n == EOF)
    *float2int = FLT2INT_TKREG;
  else
  {
    *float2int = float2int_code(tmp);
    if ( *float2int == -1 )
    {
      printf("ERROR: regio_read_register(): float2int method %s from file %s,"
             " match not found\n",tmp,regfile);
      return(1);
    }
  }

  return(0);
}
/* -------------------------------------------------------------- */
int regio_print_register(FILE *fp, char *subject, float inplaneres,
                         float betplaneres, float intensity, MATRIX *R,
                         int float2int)
{
  int r,c;
  char *f2imethod ;

  if (subject != NULL)
    fprintf(fp,"%s\n",subject);
  else
    fprintf(fp,"subjectunknown\n");

  fprintf(fp,"%f\n",inplaneres);
  fprintf(fp,"%f\n",betplaneres);
  fprintf(fp,"%f\n",intensity);

  for (r=0;r<3;r++)
  {
    for (c=0;c<4;c++)
    {
      fprintf(fp,"%e ",R->rptr[r+1][c+1]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"0 0 0 1\n");

  switch (float2int)
  {
  case FLT2INT_TKREG:
    f2imethod = "tkregister";
    break;
  case FLT2INT_FLOOR:
    f2imethod = "floor";
    break;
  case FLT2INT_ROUND:
    f2imethod = "round";
    break;
  default:
    printf("ERROR: regio_print_register: f2i code = %d, unrecoginized\n",
           float2int);
    return(1);
  }

  fprintf(fp,"%s\n",f2imethod);


  return(0);
}

/* -------------------------------------------------------------- */
int regio_write_register(char *regfile, char *subject, float inplaneres,
                         float betplaneres, float intensity, MATRIX *R,
                         int float2int)
{
  FILE *fp;

  fp = fopen(regfile,"w");
  if (fp==NULL)
  {
    perror("regio_write_register");
    return(1);
  }

  regio_print_register(fp, subject, inplaneres, betplaneres, intensity, R,
                       float2int);
  fclose(fp);

  return(0);
}
/* --------------------------------------------------------------
   regio_read_mincxfm() - reads a 3x4 transform as the last three
   lines of the xfmfile. Blank lines at the end will defeat it.
   If fileinfo != NULL, reads in the "fileinfo". This is the 3rd
   line in the minc xfm file. It will contain information about
   the center of the transform (-center). This is used by
   ltaMNIreadEx() to account for the non-zero center of the MNI
   talairach template. If the center infomation is not there,
   ltaMNIreadEx() assumes that the center is 0 and will produce
   the wrong transform. Thus, if one is going to write out the
   xfm, then one needs to keep track of this information when
   reading it in. regio_write_mincxfm() takes fileinfo as an
   argument. If one is not going to write out the xfm, then
   simply set fileinfo to NULL.
   -------------------------------------------------------------- */
int regio_read_mincxfm(char *xfmfile, MATRIX **R, char **fileinfo)
{
  FILE *fp;
  char tmpstr[1000];
  int r,c,n,nlines;
  float val;

  memset(tmpstr,'\0',1000);

  fp = fopen(xfmfile,"r");
  if (fp==NULL)
  {
    perror("regio_read_mincxfm");
    fprintf(stderr,"Could read %s\n",xfmfile);
    return(1);
  }

  fgetl(tmpstr, 900, fp) ;  /* MNI Transform File */
  if (strncmp("MNI Transform File", tmpstr, 18))
  {
    printf("ERROR: %s does not start as 'MNI Transform File'",xfmfile);
    return(1);
  }

  fgetl(tmpstr, 900, fp) ;   /* fileinfo */
  if (fileinfo != NULL)
  {
    *fileinfo = strcpyalloc(tmpstr);
    printf("\n%s\n\n",*fileinfo);
  }
  else printf("Not reading in xfm fileinfo\n");

  // Close it and open it up again to rewind it
  fclose(fp);
  fp = fopen(xfmfile,"r");

  /* Count the number of lines */
  nlines = 0;
  while (fgets(tmpstr,1000,fp) != NULL) nlines ++;
  rewind(fp);

  /* skip all but the last 3 lines */
  for (n=0;n<nlines-3;n++) fgets(tmpstr,1000,fp);

  *R = MatrixAlloc(4,4,MATRIX_REAL);
  if (*R == NULL)
  {
    fprintf(stderr,"regio_read_mincxfm(): could not alloc R\n");
    fclose(fp);
    return(1);
  }
  MatrixClear(*R);

  /* registration matrix */
  for (r=0;r<3;r++)
  { /* note: upper limit = 3 for xfm */
    for (c=0;c<4;c++)
    {
      n = fscanf(fp,"%f",&val);
      if (n != 1)
      {
        perror("regio_read_mincxfm()");
        fprintf(stderr,"Error reading R[%d][%d] from %s\n",r,c,xfmfile);
        fclose(fp);
        return(1);
      }
      (*R)->rptr[r+1][c+1] = val;
      /*printf("%7.4f ",val);*/
    }
    /*printf("\n");*/
  }
  (*R)->rptr[3+1][3+1] = 1.0;

  fclose(fp);
  return(0);
}
/* --------------------------------------------------------------
   regio_write_mincxfm() - writes a 3x4 transform in something
   like a minc xfm file. See regio_read_mincxfm() for docs on
   fileinfo.
   -------------------------------------------------------------- */
int regio_write_mincxfm(char *xfmfile, MATRIX *R, char *fileinfo)
{
  FILE *fp;
  int r,c;
  time_t time_now;

  fp = fopen(xfmfile,"w");
  if (fp==NULL)
  {
    perror("regio_write_mincxfm");
    printf("Could open %s for writing\n",xfmfile);
    return(1);
  }
  fprintf(fp,"MNI Transform File\n");
  if (fileinfo) fprintf(fp,"%s\n",fileinfo);
  fprintf(fp,"%% This file was created by %s\n",Progname);
  time(&time_now);
  fprintf(fp,"%% %s\n", ctime(&time_now));
  fprintf(fp,"\n");
  fprintf(fp,"Transform_Type = Linear;\n");
  fprintf(fp,"Linear_Transform =\n");

  for (r=0;r<3;r++)
  { /* note: upper limit = 3 for xfm */
    for (c=0;c<4;c++)
      fprintf(fp,"%e ",R->rptr[r+1][c+1]);
    if (r != 2) fprintf(fp,"\n");
    else       fprintf(fp,";\n");
  }

  fclose(fp);
  return(0);
}
/* --------------------------------------------------------------
   regio_read_xfm4() - reads a 4x4 transform as the last four
   lines of the xfmfile. Blank lines at the end will defeat it.
   -------------------------------------------------------------- */
int regio_read_xfm4(char *xfmfile, MATRIX **R)
{
  FILE *fp;
  char tmpstr[1000];
  int r,c,n,nlines;
  float val;

  memset(tmpstr,'\0',1000);

  fp = fopen(xfmfile,"r");
  if (fp==NULL)
  {
    perror("regio_read_xfm4");
    fprintf(stderr,"Could read %s\n",xfmfile);
    return(1);
  }

  /* Count the number of lines */
  nlines = 0;
  while (fgets(tmpstr,1000,fp) != NULL) nlines ++;
  rewind(fp);

  /* skip all but the last 3 lines */
  for (n=0;n<nlines-4;n++) fgets(tmpstr,1000,fp);

  *R = MatrixAlloc(4,4,MATRIX_REAL);
  if (*R == NULL)
  {
    fprintf(stderr,"regio_read_xfm4(): could not alloc R\n");
    fclose(fp);
    return(1);
  }
  MatrixClear(*R);

  /* registration matrix */
  for (r=0;r<3;r++)
  { /* note: upper limit = 3 for xfm */
    for (c=0;c<4;c++)
    {
      n = fscanf(fp,"%f",&val);
      if (n != 1)
      {
        perror("regio_read_xfm4()");
        fprintf(stderr,"Error reading R[%d][%d] from %s\n",r,c,xfmfile);
        fclose(fp);
        return(1);
      }
      (*R)->rptr[r+1][c+1] = val;
      /*printf("%7.4f ",val);*/
    }
    /*printf("\n");*/
  }
  (*R)->rptr[3+1][3+1] = 1.0;

  fclose(fp);
  return(0);
}
/* --------------------------------------------------------------
   regio_read_xfm() - reads a transform. If the file ends in .xfm,
   then the last three lines are read, and the last line of the
   transform is set to 0 0 0 1. Otherwise, the last four lines
   are read. Blank lines at the end will defeat it. This should
   be able to read tlas properly.
   -------------------------------------------------------------- */
int regio_read_xfm(char *xfmfile, MATRIX **R)
{
  char *ext, *fileinfo;
  int err = 0;

  ext = fio_extension(xfmfile);

  if (strcmp(ext,"xfm") == 0)
  {
    err = regio_read_mincxfm(xfmfile,R,&fileinfo);
    free(fileinfo);
  }
  else
    err = regio_read_xfm4(xfmfile,R);

  free(ext);

  return(err);
}
