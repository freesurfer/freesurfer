/*
  fsgdf.c
  Utilities for reading freesurfer group descriptor file format 
  $Id: fsgdf.c,v 1.10 2002/11/13 19:29:05 kteich Exp $

  See:   http://surfer.nmr.mgh.harvard.edu/docs/fsgdf.txt

  1. Tags are NOT case sensitive.
  2. Labels are case sensitive.
  3. When multiple items appear on a line, they can be 
     separated by any white space (ie, blank or tab).
  4. Any line where # appears as the first non-white space 
     character is ignored (ie, a comment).
  5. The Variables line should appear before the first Input line.
  6. All Class lines should appear before the first Input line.
  7. Variable label replications are not allowed.
  8. Class label replications are not allowed.
  9. Subject Id replications are not allowed.
 10. If a class label is not used, a warning is printed out.
 11. The DefaultVariable must be a member of the Variable list.
 12. No error is generated if a tag does not match.
 13. Empty lines are OK.
 14. A class label can optionally be followed by a class marker.
 15. A class marker can optionally be followed by a class color.

  Example of a legal file:
  ------------------------- cut here ------------------
  GroupDescriptorFile 1
  Title MyTitle
  MeasurementName thickness
  RegistrationSubject average7
  PlotFile /space/greve/1/users/greve/f_000.bfloat

  Class Class1 plus blue 
  CLASS Class2 circle green
  SomeTag
  Variables             Age  Weight   IQ
  Input subjid1 Class1   10    100   1000		
  Input subjid2 Class2   20    200   2000
  #Input subjid3 Class2   20    200   2000
  DefaultVariable Age
  ------------------------- cut here ------------------

  NOTE: SomeTag is not a valid tag, so it will be ignored.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "fsgdf.h"
#include "mri2.h"

/* This should be in ctype.h, but the compiler complains */
#ifndef Darwin
#ifndef isblank
int isblank (int c);
#endif
#endif


static FSGD *gdfReadV1(char *gdfname);
static int gdfPrintV1(FILE *fp, FSGD *gd);
static int gdfCountItemsOnLine(FILE *fp);
static int gdfCountItemsInString(char *str);
static int gdfClassNo(FSGD *gd, char *class);
static int gdfCheckVarRep(FSGD *gd);
static int gdfCheckClassRep(FSGD *gd);
static int gdfCheckAllClassesUsed(FSGD *gd);
static int gdfCheckSubjRep(FSGD *gd);
static int gdfGetDefVarLabelNo(FSGD *gd);

/*--------------------------------------------------*/
FSGD *gdfAlloc(int version)
{
  FSGD *gd;
  gd = (FSGD *) calloc(sizeof(FSGD),1);
  gd->version = version;
  return(gd);
}
/*--------------------------------------------------*/
int gdfFree(FSGD **ppgd)
{
  FSGD *gd;
  gd = *ppgd;
  if(gd->data)
    MRIfree(&gd->data);
  free(gd);
  *ppgd = NULL;
  return(0);
}
/*--------------------------------------------------*/
int gdfPrintHeader(FILE *fp, FSGD *gd)
{
  int r;

  switch(gd->version){
  case 1: r = gdfPrintV1(fp,gd); break;
  default:
    printf("ERROR: FSGDF version %d unsupported\n",gd->version);
    return(1);
  }
  return(r);
}

/*--------------------------------------------------*/
int gdfPrintStdout(FSGD *gd)
{
  return gdfPrintHeader(stdout,gd);
}

/*--------------------------------------------------*/
static int gdfPrintV1(FILE *fp, FSGD *gd)
{
  int n,m;

  if(gd->version != 1){
    fprintf(fp,"ERROR: FSGDF version = %d, should be 1 \n",gd->version);
    return(1);
  }
  fprintf(fp,"GroupDescriptorFile 1\n");

  if(strlen(gd->title) > 0)
    fprintf(fp,"Title %s\n",gd->title);
  if(strlen(gd->measname) > 0)
    fprintf(fp,"MeasurementName %s\n",gd->measname);
  if(strlen(gd->tessellation) > 0)
    fprintf(fp,"Tessellation %s\n",gd->tessellation);
  if(strlen(gd->regsubj) > 0)
    fprintf(fp,"RegistrationSubject %s\n",gd->regsubj);
  if(strlen(gd->datafile) > 0)
    fprintf(fp,"PlotFile %s\n",gd->datafile);
  if(strlen(gd->defvarlabel) > 0)
    fprintf(fp,"DefaultVariable %s\n",gd->defvarlabel);
  if(gd->nclasses > 0){
    for(n=0; n < gd->nclasses; n++){
      fprintf(fp,"Class %s",gd->classlabel[n]);
      if(strlen(gd->classmarker[n])>0)
	fprintf(fp," %s",gd->classmarker[n]);
      if(strlen(gd->classcolor[n])>0)
	fprintf(fp," %s",gd->classcolor[n]);
      fprintf(fp,"\n");
    }
  }
      
  if(gd->nvariables > 0){
    fprintf(fp,"Variables ");
    for(n=0; n < gd->nvariables; n++)
      fprintf(fp,"%s ",gd->varlabel[n]);
    fprintf(fp,"\n");
  }
      
  if(gd->ninputs > 0){
    for(n=0; n < gd->ninputs; n++){
      //fprintf(fp,"Input %s %s (id=%d) ",gd->subjid[n],
      //        gd->classlabel[gd->subjclassno[n]],gd->subjclassno[n]+1);
      fprintf(fp,"Input %s %s ",gd->subjid[n],gd->classlabel[gd->subjclassno[n]]);
      if(gd->nvariables > 0){
	for(m=0; m < gd->nvariables; m++)
	  fprintf(fp,"%g ",gd->varvals[n][m]);
      }
      fprintf(fp,"\n");
    }
  }
      
  return(0);
}

/*--------------------------------------------------*/
FSGD *gdfRead(char *gdfname, int LoadData)
{
  FSGD *gd;
  FILE *fp;
  char tmpstr[1000];
  int version=0;
  int nv;
  MRI *mritmp;
  char *dirname;
  char datafilename[1000];

  printf("gdfReadHeader: reading %s\n",gdfname);

  fp = fopen(gdfname,"r");
  if(fp==NULL){
    printf("ERROR: cannot open %s for reading\n",gdfname);
    return(NULL);
  }

  fscanf(fp,"%s",tmpstr);
  if(strcasecmp(tmpstr,"GroupDescriptorFile") != 0){
    printf("ERROR: %s is nore formated properly\n",gdfname);
    return(NULL);
  }

  fscanf(fp,"%d",&version);
  fclose(fp);

  switch(version){
  case 1: gd = gdfReadV1(gdfname); break;
  default:
    printf("ERROR: FSGDF version %d unsupported (%s) \n",version,gdfname);
    return(NULL);
  }

  /* load the MRI containing our raw data. */
  if(LoadData && strlen(gd->datafile) > 0){

    /* start with the datafile name  we got. extract the path from the
       file name of the gdf we  got, then prepend that to the datafile
       name. */
    dirname = (char*)fio_dirname(gdfname);
    if(NULL != dirname)
      {
	sprintf(datafilename,"%s/%s",dirname,gd->datafile);
      }
    else
      {
	printf("ERROR: Couldn't extract dirname from GDF header name.\n");
	strcpy(datafilename,gd->datafile);
      }

    gd->data = MRIread(datafilename);
    if(NULL == gd->data){
      printf("ERROR: Couldn't read raw data at %s \n",gd->datafile);
      gdfFree(&gd);
      return(NULL);
    }
    nv = gd->data->width * gd->data->height * gd->data->depth;
    if(strcmp(gd->tessellation,"surface")==0 && nv != gd->data->width){
      printf("INFO: gdfRead: reshaping\n");
      mritmp = mri_reshape(gd->data, nv, 1, 1, gd->data->nframes);
      MRIfree(&gd->data);
      gd->data = mritmp;
    }
  }

  return(gd);
}

/*--------------------------------------------------*/
static FSGD *gdfReadV1(char *gdfname)
{
  FSGD *gd;
  FILE *fp;
  char tag[1000];
  char tmpstr[1000];
  char class[100];
  int version,r,n,m;

  fp = fopen(gdfname,"r");
  if(fp==NULL){
    printf("ERROR: cannot open %s for reading\n",gdfname);
    return(NULL);
  }

  fscanf(fp,"%s",tag);
  if(strcasecmp(tag,"GroupDescriptorFile") != 0){
    printf("ERROR: %s is not formated properly\n",gdfname);
    return(NULL);
  }

  fscanf(fp,"%d",&version);
  if(version != 1){
    printf("ERROR: version=%d, != 1 (%s)\n",version,gdfname);
    return(NULL);
  }    

  gd = gdfAlloc(1);

  /*------- begin input loop --------------*/
  while(1){

    r = fscanf(fp,"%s",tag);
    if(r==EOF) break;

    if(!strcasecmp(tag,"Title")){
      r = fscanf(fp,"%s",gd->title);
      if(r==EOF) goto formaterror;
      continue;
    }

    if(!strcasecmp(tag,"MeasurementName")){
      r = fscanf(fp,"%s",gd->measname);
      if(r==EOF) goto formaterror;
      continue;
    }

    if(!strcasecmp(tag,"Tessellation")){
      r = fscanf(fp,"%s",gd->tessellation);
      if(r==EOF) goto formaterror;
      continue;
    }

    if(!strcasecmp(tag,"RegistrationSubject")){
      r = fscanf(fp,"%s",gd->regsubj);
      if(r==EOF) goto formaterror;
      continue;
    }

    if(!strcasecmp(tag,"PlotFile")){
      r = fscanf(fp,"%s",gd->datafile);
      if(r==EOF) goto formaterror;
      continue;
    }

    if(!strcasecmp(tag,"DefaultVariable")){
      r = fscanf(fp,"%s",gd->defvarlabel);
      if(r==EOF) goto formaterror;
      continue;
    }

    if(!strcasecmp(tag,"Class")){
      r = fscanf(fp,"%s",gd->classlabel[gd->nclasses]);
      if(r==EOF) goto formaterror;
      fgets(tmpstr,1000,fp);
      r = gdfCountItemsInString(tmpstr);
      if(r == 1) sscanf(tmpstr,"%s",gd->classmarker[gd->nclasses]);
      if(r == 2) sscanf(tmpstr,"%s %s",gd->classmarker[gd->nclasses],
			gd->classcolor[gd->nclasses]);
      gd->nclasses ++;
      continue;
    }

    if(!strcasecmp(tag,"Variables")){
      if(gd->nvariables != 0){
	printf("ERROR: multiple 'Variables' lines found\n");
	goto formaterror;
      }
      r = gdfCountItemsOnLine(fp);
      if(r==0){
	fprintf(stderr,"WARNING: no variables on 'Variables' line found\n");
	continue;
      }
      for(m=0; m < r; m++)
	fscanf(fp,"%s",gd->varlabel[m]);
      gd->nvariables = r;
      r = gdfCheckVarRep(gd);
      if(r != -1){
	printf("ERROR: variable label %s appears multiple times\n",gd->varlabel[r]);
	goto formaterror;
      }
      continue;
    }

    if(!strcasecmp(tag,"Input")){
      if(gd->nclasses == 0){
	printf("FSGDF Format Error: no classes defined before the first input line.\n");
	return(NULL);
      }
      n = gd->ninputs; /* current input number */
      r = fscanf(fp,"%s %s",gd->subjid[n],class);
      if(r==EOF){
	printf("Input line %d: ",n+1); 
	goto formaterror;
      }
      r = gdfClassNo(gd,class);
      if(r < 0){
	printf("Input line %d, subjid = %s, class %s not defined \n",
	       n+1,gd->subjid[n],class); 
	goto formaterror;
      }
      gd->subjclassno[n] = r;

      r = gdfCountItemsOnLine(fp);
      if(r != gd->nvariables){
	printf("ERROR: Input line %d, subjid = %s\n",n+1,gd->subjid[n]); 
	printf("       Found %d variables, expected. %d \n",r,gd->nvariables);
	//printf("%s\n",tmpstr);
	goto formaterror;
      }
      for(m=0; m < gd->nvariables; m++)
	fscanf(fp,"%f",&gd->varvals[n][m]);

      gd->ninputs ++;
      continue;
    }

    if(tag[0] != '#') fprintf(stderr,"INFO: ignoring tag %s \n",tag);
    fgets(tmpstr,1000,fp);

  }/*------- End loop over tags ----------------*/

  r = gdfCheckClassRep(gd);
  if(r != -1){
    printf("ERROR: class label %s appears multiple times\n",gd->classlabel[r]);
    sprintf(tag,"Class");
    goto formaterror;
  }

  r = gdfCheckAllClassesUsed(gd);
  if(r != -1)
    printf("WARNING: class %s is defined but not used.\n",gd->classlabel[r]);

  r = gdfCheckSubjRep(gd);
  if(r != -1){
    printf("ERROR: subject id %s appears multiple times\n",gd->subjid[r]);
    sprintf(tag,"Input");
    goto formaterror;
  }

 
  r = gdfGetDefVarLabelNo(gd);
  if(r == -1){
    printf("ERROR: default varible %s does not exist in list\n",gd->defvarlabel);
    sprintf(tag,"DefaultVariable");
    goto formaterror;
  }

  return(gd);

  formaterror:
  printf("FSGDF Format Error: file = %s, tag=%s\n",gdfname,tag);
  gdfFree(&gd);
  return(NULL);

}
/*--------------------------------------------------
  gdfClassNo() - returns the zero-based class number 
  associated with a class label.
  --------------------------------------------------*/
static int gdfClassNo(FSGD *gd, char *class)
{
  int code;
  for(code=0; code < gd->nclasses; code++)
    if(!strcmp(gd->classlabel[code],class)) return(code);
  return(-1);
}
/*--------------------------------------------------
  gdfCountItemsInString() returns the number of items
  in the given string, where an item is defined as
  one or more contiguous non-blank characters.
  --------------------------------------------------*/
static int gdfCountItemsInString(char *str)
{
  int len, n, nhits;

  len = strlen(str);

  nhits = 0;
  n = 0;
  while(n < len){
    while(isblank(str[n])) n++;
    if(n >= len) break;
    if(str[n] == '\0' || str[n] == '\n' || str[n] == '\r') break;
    while(!isblank(str[n])) n++;
    nhits++;
  }

  //printf("nhits %d\n",nhits);

  return(nhits);
}
/*--------------------------------------------------
  gdfCountItemsOnLine() returns the number of items
  before the next newline, where an item is defined as
  one or more contiguous non-blank characters. The
  file pointer is not changed.
  --------------------------------------------------*/
static int gdfCountItemsOnLine(FILE *fp)
{
  fpos_t now;
  char tmpstr[1000];
  int nitems;

  fgetpos(fp,&now);
  fgets(tmpstr,1000,fp);
  fsetpos(fp,&now);

  nitems = gdfCountItemsInString(tmpstr);

  return(nitems);
}
/*--------------------------------------------------
  gdfCheckClassRep() - checks whether there are any
  repetitions in the class labels. Returns 0 if
  there are no reps.
  --------------------------------------------------*/
static int gdfCheckClassRep(FSGD *gd)
{
  int n,m;

  for(n=0; n < gd->nclasses; n++)
    for(m=n+1; m < gd->nclasses; m++)
      if(strcmp(gd->classlabel[n],gd->classlabel[m])==0)
	return(n);
  return(-1);
}
/*--------------------------------------------------
  gdfCheckVarRep() - checks whether there are any
  repetitions in the variable labels. Returns -1 if
  there are no reps, otherwise returns the index
  of one of the reps.
  --------------------------------------------------*/
static int gdfCheckVarRep(FSGD *gd)
{
  int n,m;

  for(n=0; n < gd->nvariables; n++)
    for(m=n+1; m < gd->nvariables; m++)
      if(strcmp(gd->varlabel[n],gd->varlabel[m])==0)
	return(n);
  return(-1);
}
/*--------------------------------------------------
  gdfCheckSubjRep() - checks whether there are any
  repetitions in the subject labels. Returns -1 if
  there are no reps, otherwise returns the index
  of one of the reps.
  --------------------------------------------------*/
static int gdfCheckSubjRep(FSGD *gd)
{
  int n,m;

  for(n=0; n < gd->ninputs; n++)
    for(m=n+1; m < gd->ninputs; m++)
      if(strcmp(gd->subjid[n],gd->subjid[m])==0)
	return(n);
  return(-1);
}
/*--------------------------------------------------
  gdfCheckAllClassesUsed() - checks whether all classes
  defined are present in the subject data. Returns -1
  if everything is ok, otherwise returns the class
  number that is not used.
  --------------------------------------------------*/
static int gdfCheckAllClassesUsed(FSGD *gd)
{
  int n,m,ok;
  
  for(n=0; n < gd->nclasses; n++){
    ok = 0;
    for(m=0; m < gd->ninputs; m++)
      if(gd->subjclassno[m] == n) ok = 1;
    if(!ok) return(n);
  }

  return(-1);
}
/*--------------------------------------------------
  gdfGetDefVarLabelNo() - returns the label number
  of the default variable. If the default variable
  is not defined, 0 is returned. If there is no
  match, -1 is returned. Otherwise, returns the
  index.
  --------------------------------------------------*/
static int gdfGetDefVarLabelNo(FSGD *gd)
{
  int n;

  if(strlen(gd->defvarlabel) == 0) return(0);

  for(n=0; n < gd->nvariables; n++)
    if(strcmp(gd->varlabel[n],gd->defvarlabel)==0)
	return(n);
  return(-1);
}
/*--------------------------------------------------------
  gdfMatrixDOSS() - creates a design matrix that models each
  class as having a Different Offset, but models having the 
  Same Slope
  ---------------------------------------------------------*/
MATRIX *gdfMatrixDOSS(FSGD *gd, MATRIX *X)
{
  int nrows, ncols, r,v,c;

  nrows = gd->ninputs;
  ncols = gd->nclasses + gd->nvariables;

  X = MatrixZero(nrows,ncols,X);
  if(X==NULL) return(NULL);

  for(r=0; r<nrows; r++){

    c = gd->subjclassno[r];
    X->rptr[r+1][c+1] = 1;
    
    for(v = 0; v < gd->nvariables; v++){
      c = v + gd->nclasses;
      X->rptr[r+1][c+1] = gd->varvals[r][v];
    }
  }

  return(X);
}
/*--------------------------------------------------------
  gdfMatrixDODS() - creates a design matrix that models each
  class as having a Different Offset and Different Slope. The
  number of rows in the matrix will equal the number of 
  subjets. The number of columns will equal the number of 
  classes times the (number of varibles + 1). The first
  Nclasses columns will be the offsets, the next Nclasses
  columns will be the first variable, etc.
  ---------------------------------------------------------*/
MATRIX *gdfMatrixDODS(FSGD *gd, MATRIX *X)
{
  int nrows, ncols, n,r,v,c;

  nrows = gd->ninputs;
  ncols = gd->nclasses * ( gd->nvariables + 1);

  X = MatrixZero(nrows,ncols,X);
  if(X==NULL) return(NULL);

  for(r=0; r < nrows; r++){

    n = gd->subjclassno[r];
    c = n;
    X->rptr[r+1][c+1] = 1;
    
    for(v = 0; v < gd->nvariables; v++){
      c += gd->nclasses;
      X->rptr[r+1][c+1] = gd->varvals[r][v];
    }
  }

  return(X);
}

/*---------------------------------------------------*/
int gdfCheckMatrixMethod(char *gd2mtx_method)
{
  if( strcmp(gd2mtx_method,"doss") == 0 ||
      strcmp(gd2mtx_method,"dods") == 0 ||
      strcmp(gd2mtx_method,"none") == 0  
      ) return(0);

  printf("ERROR: gd2mtx method %s unrecoginzied.\n",gd2mtx_method);
  printf("       Legal values are dods, doss, and none.\n");
  return(1);
}
/*---------------------------------------------------*/
MATRIX *gdfMatrix(FSGD *gd, char *gd2mtx_method, MATRIX *X)
{
  if(gdfCheckMatrixMethod(gd2mtx_method)) return(NULL);

  if(strcmp(gd2mtx_method,"none") == 0){
    printf("ERROR: cannot create matrix when method is none\n");
    return(NULL);
  }

  if(strcmp(gd2mtx_method,"doss") == 0)
    X = gdfMatrixDOSS(gd,X);
  if(strcmp(gd2mtx_method,"dods") == 0)
    X = gdfMatrixDODS(gd,X);

  return(X);
}

/*------------------------------------------------------------
  gdfGetTitle() - copies the title into the output argument.
  ------------------------------------------------------------*/
int gdfGetTitle(FSGD *gd, char *title)
{
  if(NULL == gd)
    return(-1);

  strcpy(title,gd->title);

  return(0);
}

/*------------------------------------------------------------
  gdfGetMeasurementName() - copies the measurment name into
  the output argument.
  ------------------------------------------------------------*/
int gdfGetMeasurementName(FSGD *gd, char *name)
{
  if(NULL == gd)
    return(-1);

  strcpy(name,gd->measname);

  return(0);
}

/*------------------------------------------------------------
  gdfGetSubjectName() - copies the subject name into the 
  output argument.
  ------------------------------------------------------------*/
int gdfGetSubjectName(FSGD *gd, char *name)
{
  if(NULL == gd)
    return(-1);

  strcpy(name,gd->regsubj);

  return(0);
}

/*------------------------------------------------------------
  gdfGetDataFileName() - copies the data file name into the
  output argument.
  ------------------------------------------------------------*/
int gdfGetDataFileName(FSGD *gd, char *filename)
{
  if(NULL == gd)
    return(-1);

  strcpy(filename,gd->datafile);

  return(0);
}

/*------------------------------------------------------------
  gdfGetNumClasses() - returns the number of classes in the
  output argument.
  ------------------------------------------------------------*/
int gdfGetNumClasses(FSGD *gd, int *nclasses)
{
  if(NULL == gd)
    return(-1);

  *nclasses = gd->nclasses;

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthClassLabel() - copies the nth class label into the
  output argument where nclass is from 0 -> nclasses.
  ------------------------------------------------------------*/
int gdfGetNthClassLabel(FSGD *gd, int nclass, char *label)
{
  if(NULL == gd)
    return(-1);
  if(nclass < 0 || nclass >= gd->nclasses)
    return(-1);

  strcpy(label,gd->classlabel[nclass]);

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthClassMarker() - copies the nth class marker into the
  output argument where nclass is from 0 -> nclasses.
  ------------------------------------------------------------*/
int gdfGetNthClassMarker(FSGD *gd, int nclass, char *marker)
{
  if(NULL == gd)
    return(-1);
  if(nclass < 0 || nclass >= gd->nclasses)
    return(-1);

  strcpy(marker,gd->classmarker[nclass]);

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthClassColor() - copies the nth class color into the
  output argument where nclass is from 0 -> nclasses.
  ------------------------------------------------------------*/
int gdfGetNthClassColor(FSGD *gd, int nclass, char *color)
{
  if(NULL == gd)
    return(-1);
  if(nclass < 0 || nclass >= gd->nclasses)
    return(-1);

  strcpy(color,gd->classcolor[nclass]);

  return(0);
}

/*------------------------------------------------------------
  gdfGetNumVariables() - returns the number of variables in the
  output argument.
  ------------------------------------------------------------*/
int gdfGetNumVariables(FSGD *gd, int *nvariables)
{
  if(NULL == gd)
    return(-1);

  *nvariables = gd->nvariables;

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthVariableLabel() - copies the nth variable label into
  the output argument.
  ------------------------------------------------------------*/
int gdfGetNthVariableLabel(FSGD *gd, int nvariable, char *label)
{
  if(NULL == gd)
    return(-1);
  if(nvariable < 0 || nvariable >= gd->nvariables)
    return(-1);

  strcpy(label,gd->varlabel[nvariable]);

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthVariableDefault() - copies of name of
  the default variable to the output argument.
  ------------------------------------------------------------*/
int gdfGetDefaultVariable(FSGD *gd, char *label)
{
  if(NULL == gd)
    return(-1);
  
  strcpy(label,gd->defvarlabel);

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthVariableDefaultIndex() - copies of index of
  the default variable to the output parameter
  ------------------------------------------------------------*/
int gdfGetDefaultVariableIndex(FSGD *gd, int *nvariable)
{
  if(NULL == gd)
    return(-1);
  
  *nvariable = gdfGetDefVarLabelNo(gd);

  return(0);
}

/*------------------------------------------------------------
  gdfGetNumSubjects() - returns the number of subjects in the
  output parameter.
  ------------------------------------------------------------*/
int gdfGetNumSubjects(FSGD *gd, int *nsubjects)
{
  if(NULL == gd)
    return(-1);

  *nsubjects = gd->ninputs;

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthSubjectID() - copies the id of the nth subject into
  the ouput parameter where nsubject is from 0 -> ninputs.
  ------------------------------------------------------------*/
int gdfGetNthSubjectID(FSGD *gd, int nsubject, char *id)
{
  if(NULL == gd)
    return(-1);
  if(nsubject < 0 || nsubject >= gd->ninputs)
    return(-1);

  strcpy(id,gd->subjid[nsubject]);

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthSubjectClass() - returns the index of the nth subject's
  class in the output parameter where nsubject is from 0 -> ninputs.
  ------------------------------------------------------------*/
int gdfGetNthSubjectClass(FSGD *gd, int nsubject, int *class)
{
  if(NULL == gd)
    return(-1);
  if(nsubject < 0 || nsubject >= gd->ninputs)
    return(-1);

  *class = gd->subjclassno[nsubject];

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthSubjectNthValue() - returns the index of the nth
  subject's nth variable value where nsubject is from 0 -> ninputs
  and nvariable is from 0 -> nvariables.
  ------------------------------------------------------------*/
int gdfGetNthSubjectNthValue(FSGD *gd, int nsubject, 
			     int nvariable, float *value)
{
  if(NULL == gd)
    return(-1);
  if(nsubject < 0 || nsubject >= gd->ninputs)
    return(-1);
  if(nvariable < 0 || nvariable >= gd->nvariables)
    return(-1);

  *value = gd->varvals[nsubject][nvariable];

  return(0);
}

/*------------------------------------------------------------
  gdfGetNthSubjectMeasurement() - returns a measurement value
  in the output paramter for the nth subject where nsubject 
  is 0 -> nsubjects, getting the data out of the MRI data. the
  meaning of (x,y,z) is contextual; in a surface x=y=1 and z=vno
  and in a volume it's a normal anatomical coordinate.
  ------------------------------------------------------------*/
int gdfGetNthSubjectMeasurement(FSGD *gd, int nsubject, 
				int x, int y, int z, float *value)
{
  if(NULL == gd)
    return(-1);
  if(nsubject < 0 || nsubject >= gd->ninputs)
    return(-1);
#if 0
  if(x < gd->data->xstart || x > gd->data->xend ||
     y < gd->data->ystart || y > gd->data->yend ||
     z < gd->data->zstart || z > gd->data->zend ||
     nsubject < 0 || nsubject >= gd->data->nframes)
    return(-1);
#endif    

  switch( gd->data->type ) {
    case MRI_UCHAR:
      *value = MRIseq_vox(gd->data,x,y,z,nsubject);
      break;
    case MRI_INT:
      *value = MRIIseq_vox(gd->data,x,y,z,nsubject);
      break;
    case MRI_LONG:
      *value = MRILseq_vox(gd->data,x,y,z,nsubject);
      break;
    case MRI_FLOAT:
      *value = MRIFseq_vox(gd->data,x,y,z,nsubject);
      break;
    case MRI_SHORT:
      *value = MRISseq_vox(gd->data,x,y,z,nsubject);
      break;
    default:
      break ;
    }

  return(0);
}
