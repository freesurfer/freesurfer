#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "resample.h"
#include "registerio.h"


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
  if(fp==NULL){
    perror("regio_read_register()");
    fprintf(stderr,"Could not open %s\n",regfile);
    return(1);
  }

  /* subject name */
  n = fscanf(fp,"%s",tmp);
  if(n != 1){
    perror("regio_read_register()");
    fprintf(stderr,"Error reading subject from %s\n",regfile);
    fclose(fp);
    return(1);
  }
  *subject = (char *) calloc(strlen(tmp)+2,sizeof(char));
  sprintf(*subject,"%s",tmp);

  /* in-plane resolution */
  n = fscanf(fp,"%f",inplaneres);
  if(n != 1){
    perror("regio_read_register()");
    fprintf(stderr,"Error reading inplaneres from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  /* between-plane resolution */
  n = fscanf(fp,"%f",betplaneres);
  if(n != 1){
    perror("regio_read_register()");
    fprintf(stderr,"Error reading betplaneres from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  /* intensity*/
  n = fscanf(fp,"%f",intensity);
  if(n != 1){
    perror("regio_read_register()");
    fprintf(stderr,"Error reading intensity from %s\n",regfile);
    fclose(fp);
    return(1);
  }

  *R = MatrixAlloc(4,4,MATRIX_REAL);
  if(*R == NULL){
    fprintf(stderr,"regio_read_register(): could not alloc R\n");
    fclose(fp);
    return(1);
  }

  /* registration matrix */
  for(r=0;r<4;r++){
    for(c=0;c<4;c++){
      n = fscanf(fp,"%f",&val);
      if(n != 1){
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

  if(n == EOF) 
    *float2int = FLT2INT_TKREG;
  else{
    *float2int = float2int_code(tmp);
    if( *float2int == -1 ){
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
  
  if(subject != NULL)
    fprintf(fp,"%s\n",subject);
  else
    fprintf(fp,"subjectunknown\n");

  fprintf(fp,"%f\n",inplaneres);
  fprintf(fp,"%f\n",betplaneres);
  fprintf(fp,"%f\n",intensity);

  for(r=0;r<4;r++){
    for(c=0;c<4;c++){
      fprintf(fp,"%8.5f ",R->rptr[r+1][c+1]);
    }
    fprintf(fp,"\n");
  }

  switch(float2int){
  case FLT2INT_TKREG: f2imethod = "tkregister"; break;
  case FLT2INT_FLOOR: f2imethod = "floor"; break;
  case FLT2INT_ROUND: f2imethod = "round"; break;
  default: 
    printf("ERROR: regio_print_register: f2i code = %d, unrecoginized\n",
     float2int);
    fclose(fp);
    return(1);
  }

  fprintf(fp,"%s\n",f2imethod);
  fclose(fp);


  return(0);
}

/* -------------------------------------------------------------- */
int regio_write_register(char *regfile, char *subject, float inplaneres, 
       float betplaneres, float intensity, MATRIX *R,
       int float2int)
{
  FILE *fp;

  fp = fopen(regfile,"w");
  if(fp==NULL){
    perror("regio_write_register");
    return(1);
  }

  regio_print_register(fp, subject, inplaneres, betplaneres, intensity, R,
           float2int);
  fclose(fp);

  return(0);
}
/* -------------------------------------------------------------- */
int regio_read_mincxfm(char *xfmfile, MATRIX **R)
{
  FILE *fp;
  char tmpstr[1000];
  int r,c,n,nlines;
  float val;

  memset(tmpstr,'\0',1000);

  fp = fopen(xfmfile,"r");
  if(fp==NULL){
    perror("regio_read_xfm");
    fprintf(stderr,"Could read %s\n",xfmfile);
    return(1);
  }

  /* Count the number of lines */
  nlines = 0;
  while(fgets(tmpstr,1000,fp) != NULL) nlines ++;
  rewind(fp);

  /* skip all but the last 3 lines */
  for(n=0;n<nlines-3;n++) fgets(tmpstr,1000,fp);

  *R = MatrixAlloc(4,4,MATRIX_REAL);
  if(*R == NULL){
    fprintf(stderr,"regio_read_xfm(): could not alloc R\n");
    fclose(fp);
    return(1);
  }
  MatrixClear(*R);

  /* registration matrix */
  for(r=0;r<3;r++){ /* note: upper limit = 3 for xfm */
    for(c=0;c<4;c++){
      n = fscanf(fp,"%f",&val);
      if(n != 1){
  perror("regio_read_xfm()");
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
