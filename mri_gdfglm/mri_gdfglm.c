/*
  Name:    mri_gdfglm.c
  Author:  Douglas N. Greve 
  email:   analysis-bugs@nmr.mgh.harvard.edu
  Date:    4/4/03
  Purpose: performs glm analysis given group descirptor file
           and dependent variable table
  $Id: mri_gdfglm.c,v 1.1 2003/04/05 23:04:53 greve Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "error.h"
#include "diag.h"
#include "proto.h"

#include "matrix.h"

typedef struct tagDEPVARTABLE{
  char *dvtfile;
  char **RowNames;
  char **ColNames;
  MATRIX *D;
} DEPVARTABLE, DVT;

DVT *ReadDVT(char *dvtfile);
int DumpDVT(DVT *dvt, FILE *fp);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  isflag(char *flag);
static int  nth_is_arg(int nargc, char **argv, int nth);
static int  singledash(char *flag);
//static int  stringmatch(char *str1, char *str2);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_gdfglm.c,v 1.1 2003/04/05 23:04:53 greve Exp $";
char *Progname = NULL;

typedef struct tagCOVARPRUNE{
  char *Covar;
  float Min, Max;
} COVARPRUNE;

int debug = 0;

char *GDFile;
char *DVTFile; 

int  nClassList = 0;
char *ClassList[100];

int   nwClass = 0;
float wClass[200];

int  nCovarList = 0;
char *CovarList[100];

int  nwCovar = 0;
float wCovar[200];

int  nCovarPrune = 0;
COVARPRUNE CovarPrune[200];
int TestOffset = 0;

int   nDepVarList = 0;
char  *DepVarList[100];

int   nwDepVar = 0;
float wDepVar[200];

char *WLMSDepVar = NULL;

char *OutBase = NULL;

DVT *dvt;

/*---------------------------------------------------------------*/
int main(int argc, char **argv)
{

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  dump_options(stdout);

  return(0);
}
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;

    else if (!strcmp(option, "--gdf")){
      if(nargc < 1) argnerr(option,1);
      GDFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--dvtf")){
      if(nargc < 1) argnerr(option,1);
      DVTFile = pargv[0];
      dvt = ReadDVT(DVTFile);
      DumpDVT(dvt,stdout);
      nargsused = 1;
    }
    else if (!strcmp(option, "--classes")){
      if(nargc < 1) argnerr(option,1);
      while(nth_is_arg(nargc, pargv,  nargsused)){
	ClassList[nClassList] = pargv[nargsused];
	nClassList++;
	nargsused++;
      }
    }
    else if (!strcmp(option, "--covar")){
      if(nargc < 1) argnerr(option,1);
      while(nth_is_arg(nargc, pargv, nargsused)){
	CovarList[nCovarList] = pargv[nargsused];
	nCovarList++;
	nargsused++;
      }
    }
    else if (!strcmp(option, "--covarprune")){
      if(nargc < 3) argnerr(option,3);
      CovarPrune[nCovarPrune].Covar = pargv[0];
      sscanf(pargv[1],"%f",&CovarPrune[nCovarPrune].Min);
      sscanf(pargv[2],"%f",&CovarPrune[nCovarPrune].Max);
      if(CovarPrune[nCovarPrune].Min > CovarPrune[nCovarPrune].Max){
	printf("ERROR: covarprune min > max\n");
	exit(1);
      }
      nCovarPrune++;
      nargsused = 3;
    }
    else if (!strcmp(option, "--depvar")){
      if(nargc < 1) argnerr(option,1);
      while(nth_is_arg(nargc, pargv, nargsused)){
	DepVarList[nDepVarList] = pargv[nargsused];
	nDepVarList++;
	nargsused++;
      }
    }
    else if (!strcmp(option, "--wdepvar")){
      if(nargc < 1) argnerr(option,1);
      while(nth_is_arg(nargc, pargv, nargsused)){
	sscanf(pargv[nargsused],"%f",&wDepVar[nwDepVar]);
	nwDepVar++;
	nargsused++;
      }
    }
    else if (!strcmp(option, "--wlms")){
      if(nargc < 1) argnerr(option,1);
      WLMSDepVar = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--wclass")){
      if(nargc < 1) argnerr(option,1);
      while(nth_is_arg(nargc, pargv, nargsused)){
	sscanf(pargv[nargsused],"%f",&wClass[nwClass]);
	nwClass++;
	nargsused++;
      }
    }
    else if (!strcmp(option, "--wcovar")){
      if(nargc < 1) argnerr(option,1);
      while(nth_is_arg(nargc, pargv, nargsused)){
	sscanf(pargv[nargsused],"%f",&wCovar[nwCovar]);
	nwCovar++;
	nargsused++;
      }
    }
    else if (!strcmp(option, "--testoffset")){
      TestOffset = 1;
    }
    else if (!strcmp(option, "--o")){
      if(nargc < 1) argnerr(option,1);
      OutBase = pargv[0];
      nargsused = 1;
    }
    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
	fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --gdf gdfile : group descriptor file\n");
  printf("   --dvt dvtfile : dependent variable file  \n");

  printf("\n");
  printf("   --classes Class1 Class2 : use subset of classes  \n");
  printf("   --covar Covar1 Corvar2 ... : use subset of covars\n");
  printf("   --covarprune Covar Min Max : exclude when out of range\n");
  printf("   --depvar DepVar1 <DepVar2 ...> : spec dependent variables \n");
  printf("   --wdepvar wdv1 wdv2 ... : weight depvars (default is 1) \n");
  printf("   --wlms DepVar \n");


  printf("\n");
  printf("   --wclass wc1 wc2 ... : Class weights (def 1)\n");
  printf("   --wcovar wcv1 wcvw ... : Covar slope weights  \n");
  printf("   --testoffset : test offset, not covariat slope \n");

  printf("\n");
  printf("   --o output base name\n");

  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(
  "Performs glm analysis given group descirptor file
   and dependent variable table\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if(n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(GDFile == NULL){
    printf("ERROR: no gdf file\n");
    exit(1);
  }
  if(DVTFile == NULL){
    printf("ERROR: no dep var table file\n");
    exit(1);
  }
  if(OutBase == NULL){
    printf("ERROR: no output specified\n");
    exit(1);
  }
  if(nDepVarList == 0){
    printf("ERROR: no dep vars specified\n");
    exit(1);
  }
  if(nwClass != 0 && nClassList != 0 && nwClass != nClassList ){
    printf("ERROR: diff number of classes and class weights\n");
    exit(1);
  }
  if(nwCovar != 0 && nCovarList != 0 && nwCovar != nCovarList ){
    printf("ERROR: diff number of covars and covar weights\n");
    exit(1);
  }
  if(nwDepVar != 0 && nDepVarList != 0 && nwDepVar != nDepVarList ){
    printf("ERROR: diff number of dep vars and weights\n");
    exit(1);
  }
  /* Check that DepVars, Covars, and Classes are not repeated */

  if(TestOffset && nwCovar != 0){
    printf("ERROR: cannot test offset and set covar weights\n");
    exit(1);
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  int n;
  fprintf(fp,"GDFile %s \n",GDFile);
  fprintf(fp,"DVTFile %s \n", DVTFile);

  if(nClassList > 0){
    fprintf(fp,"%d subset of classes selected \n",nClassList);    
    for(n = 0; n < nClassList; n++){
      fprintf(fp,"%2d %s \n",n+1,ClassList[n]);
    }
  }
  else
    fprintf(fp,"Using All Classes\n");    

  if(nwClass > 0){
    fprintf(fp,"Class Weights\n");    
    for(n = 0; n < nwClass; n++){
      fprintf(fp,"%2d %g \n",n+1,wClass[n]);
    }
  }
  else
    fprintf(fp,"All Classes equally weighted.\n");    

  if(nCovarList > 0){
    fprintf(fp,"%d subset of covars selected \n",nCovarList);    
    for(n = 0; n < nCovarList; n++){
      fprintf(fp,"%2d %s \n",n+1,CovarList[n]);
    }
  }
  else
    fprintf(fp,"Using All Covars\n");    

  if(nwCovar > 0){
    fprintf(fp,"Covar Weights\n");    
    for(n = 0; n < nwCovar; n++){
      fprintf(fp,"%2d %g \n",n+1,wCovar[n]);
    }
  }
  else
    fprintf(fp,"All Covars equally weighted.\n");    

  if(nCovarPrune){
    fprintf(fp,"Pruning used on %d Covars\n",nCovarPrune);
    for(n = 0; n < nCovarPrune; n++){
      fprintf(fp,"  %2d  %s  %g  %g\n",n+1,CovarPrune[n].Covar,
	      CovarPrune[n].Min,CovarPrune[n].Max);
    }
  }



  fprintf(fp,"Dependent Variables (%d)\n",nDepVarList);    
  for(n = 0; n < nDepVarList; n++){
    fprintf(fp,"%2d %s ",n+1,DepVarList[n]);
    if(nwDepVar > 0)
      fprintf(fp," %g \n",wDepVar[n]);
    else
      fprintf(fp," 1\n");
  }

  if(TestOffset)
    fprintf(fp,"Testing Offsets \n");

  if(WLMSDepVar != NULL)
    fprintf(fp,"WLMSDepVar: %s \n",WLMSDepVar);

  fprintf(fp,"OutBase: %s \n",OutBase);
  fprintf(fp," \n");
  fprintf(fp," \n");

  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int isflag(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth)
{
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */
  /* This is typically called from parse_commandline(), for example,
     the following code specifies that --psdwin requires two arguments
     but there may be a third. The first two are read in, then
     nth_is_arg() is called to determine whether the third exists
     and whether it is an argument or a flag. It is not neccessary
     to check whether there are enough arguments (ie, if nargc > 2)
     because this is done in nth_is_arg().
     
    else if (stringmatch(option, "--psdwin")){
      if(nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%f",&PSDMin);
      sscanf(pargv[1],"%f",&PSDMax);
      nargsused = 2;
      if(nargc > 2 && nth_is_arg(2,pargv,nargc)){
	sscanf(pargv[2],"%f",&dPSDfTR);
	nargsused++;
      }
      else dPSDfTR = 1.0;
    }
  */
     

  /* check that there are enough args for nth to exist */
  if(nargc <= nth) return(0); 

  /* check whether the nth arg is a flag */
  if(isflag(argv[nth])) return(0);

  return(1);
}
/*------------------------------------------------------------*/
//static int stringmatch(char *str1, char *str2)
//{
//  if(! strcmp(str1,str2)) return(1);
//  return(0);
//}


DVT *ReadDVT(char *dvtfile)
{
  DVT *dvt = NULL;
  FILE *fp;
  int nrows, ncols, r, c, lentmpstr;
  char tmpstr[4001];

  if(!fio_FileExistsReadable(dvtfile)){
    printf("ERROR: %s does not exist or is not readable\n",dvtfile);
    return(NULL);
  }

  nrows = fio_NLines(dvtfile);
  if(nrows == 0){
    printf("ERROR: reading %s, no rows found\n",dvtfile);
    return(NULL);
  }

  fp = fopen(dvtfile,"r");
  fgets(tmpstr,4000,fp);
  ncols = gdfCountItemsInString(tmpstr);
  if(ncols == 0){
    printf("ERROR: reading %s, no cols found\n",dvtfile);
    fclose(fp);
    return(NULL);
  }
  fclose(fp);

  nrows = nrows-1;
  ncols = ncols-1;

  printf("%s  nrows = %d, ncols = %d\n",dvtfile,nrows,ncols);

  dvt = (DVT *)calloc(sizeof(DVT),1);

  dvt->dvtfile = (char *) calloc(sizeof(char),strlen(dvtfile)+1);
  memcpy(dvt->dvtfile,dvtfile,strlen(dvtfile));

  dvt->D = MatrixAlloc(nrows,ncols,MATRIX_REAL);
  dvt->RowNames = (char **)calloc(sizeof(char *),nrows);
  dvt->ColNames = (char **)calloc(sizeof(char *),ncols);

  fp = fopen(dvtfile,"r");

  /* Read the first row = names of columns */
  fscanf(fp,"%s",tmpstr); /* swallow the first one */
  for(c=0; c < ncols; c++){
    fscanf(fp,"%s",tmpstr);
    lentmpstr = strlen(tmpstr);
    dvt->ColNames[c] = (char *) calloc(sizeof(char),lentmpstr+1);
    memcpy(dvt->ColNames[c],tmpstr,lentmpstr);
  }

  for(r=0; r < nrows; r++){
    /* First column is row name */
    fscanf(fp,"%s",tmpstr);
    lentmpstr = strlen(tmpstr);
    dvt->RowNames[r] = (char *) calloc(sizeof(char),lentmpstr+1);
    memcpy(dvt->RowNames[r],tmpstr,lentmpstr);
    printf("r=%d, RowName = %s\n",r,dvt->RowNames[r]);
    for(c=0; c < ncols; c++){
      fscanf(fp,"%f",&dvt->D->rptr[r+1][c+1]);
    }
  }

  fclose(fp);
  return(dvt);
}

/*-----------------------------------------------------*/
int DumpDVT(DVT *dvt, FILE *fp)
{
  int c,r;

  if(fp == NULL) fp = stdout;

  fprintf(fp,"%s \n",dvt->dvtfile);
  fprintf(fp,"nrows = %d\n",dvt->D->rows);
  fprintf(fp,"ncols = %d\n",dvt->D->cols);

  fprintf(fp,"\n");
  for(c = 0; c < dvt->D->cols; c++)
    fprintf(fp,"%3d %s\n",c+1,dvt->ColNames[c]);

  fprintf(fp,"\n");
  for(r = 0; r < dvt->D->rows; r++)
    fprintf(fp,"%3d %s\n",r+1,dvt->RowNames[r]);

  fprintf(fp,"\n");

  return(0);

}
