#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
  float2int - 0=tkreg.
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

  n = fscanf(fp,"%d",float2int);
  if(n == 0) *float2int = 0;

  fclose(fp);

  return(0);
}
/* -------------------------------------------------------------- */
int regio_print_register(FILE *fp, char *subject, float inplaneres, 
			 float betplaneres, float intensity, MATRIX *R,
			 int float2int)
{
  int r,c;
  
  fprintf(fp,"%s\n",subject);
  fprintf(fp,"%f\n",inplaneres);
  fprintf(fp,"%f\n",betplaneres);
  fprintf(fp,"%f\n",intensity);

  for(r=0;r<4;r++){
    for(c=0;c<4;c++){
      fprintf(fp,"%8.5f ",R->rptr[r+1][c+1]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"%d\n",float2int);

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
  int r,c,n;
  float val;

  memset(tmpstr,'\0',1000);

  fp = fopen(xfmfile,"r");
  if(fp==NULL){
    perror("regio_read_xfm");
    fprintf(stderr,"Could read %s\n",xfmfile);
    return(1);
  }

  /* skip 5 lines */
  fgets(tmpstr,1000,fp);
  if(strcasecmp(tmpstr,"MNI Transform File\n") != 0){
    fprintf(stderr,"%s does not appear to be a MNI xfm file\n",xfmfile);
    fprintf(stderr,"--%s--\n",tmpstr);
    fclose(fp);
    return(1);
  }

  fgets(tmpstr,1000,fp);
  fgets(tmpstr,1000,fp);
  fgets(tmpstr,1000,fp);
  fgets(tmpstr,1000,fp);

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
    }
  }
  (*R)->rptr[3+1][3+1] = 1.0;

  fclose(fp);

  return(0);
}
