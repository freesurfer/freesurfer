/**
 * @brief Utilities for reading freesurfer group descriptor file format
 *
 * See:   http://surfer.nmr.mgh.harvard.edu/docs/fsgdf.txt
 *  1. Tags are NOT case sensitive.
 *  2. Labels are case sensitive.
 *  3. When multiple items appear on a line, they can be
 *     separated by any white space (ie, blank or tab).
 *  4. Any line where # appears as the first non-white space
 *     character is ignored (ie, a comment).
 *  5. The Variables line should appear before the first Input line.
 *  6. All Class lines should appear before the first Input line.
 *  7. Variable label replications are not allowed.
 *  8. Class label replications are not allowed.
 *  9. Subject Id replications are not allowed.
 * 10. If a class label is not used, a warning is printed out.
 * 11. The DefaultVariable must be a member of the Variable list.
 * 12. No error is generated if a tag does not match.
 * 13. Empty lines are OK.
 * 14. A class label can optionally be followed by a class marker.
 * 15. A class marker can optionally be followed by a class color.
 *
 *  Example of a legal file:
 *  ------------------------- cut here ------------------
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
 *   ------------------------- cut here ------------------
 *
 *  NOTE: SomeTag is not a valid tag, so it will be ignored.
 */
/*
 * Original Author: Doug Greve
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <string>
#include <locale>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <float.h>
#include "mri2.h"
#include "fio.h"
#include "matfile.h"
#include "stats.h"
#include "fsenv.h"
#include "utils.h"
#include "proto.h"
#include "diag.h"

#define FSGDF_SRC
#include "fsgdf.h"
#undef FSGDF_SRC

/* This should be in ctype.h, but the compiler complains */
#ifndef Darwin
#ifndef isblank
int isblank (int c);
#endif
#endif

static FSGD *gdfReadV1(const char *gdfname);
static int gdfPrintV1(FILE *fp, FSGD *gd);
static int gdfCheckVarRep(FSGD *gd);
static int gdfCheckClassRep(FSGD *gd);
static int gdfCheckAllClassesUsed(FSGD *gd);
static int gdfCheckSubjRep(FSGD *gd);
static int gdfGetDefVarLabelNo(FSGD *gd);

/* RKT - hack to get the .so to have Progname declared. I hate this. */
#ifdef DECLARE_PROGNAME
const char *Progname = "fsgdf";
#endif

/*--------------------------------------------------*/
FSGD *gdfAlloc(int version) {
  FSGD *gd;
  gd = (FSGD *) calloc(sizeof(FSGD),1);
  gd->version = version;
  gd->ResFWHM = -1;
  gd->LogY = 0;
  return(gd);
}


/*--------------------------------------------------*/
int gdfFree(FSGD **ppgd) {
  FSGD *gd;
  int n;
  gd = *ppgd;
  if (gd->data)  MRIfree(&gd->data);
  if (gd->X) MatrixFree(&gd->X);
  if (gd->T) MatrixFree(&gd->T);
  for (n=0; n < gd->nvarsfromfile; n++) {
    free(gd->tablefile[n]);
    free(gd->varfield[n]);
  }
  free(gd);
  *ppgd = NULL;
  return(0);
}


/*--------------------------------------------------*/
int gdfWrite(const char *gdfname, FSGD *gd) {
  FILE *fp;

  fp = fopen(gdfname,"w");
  if (fp == NULL) {
    printf("ERROR: could not open %s for writing\n",gdfname);
    return(1);
  }
  gdfPrintHeader(fp,gd);
  fclose(fp);
  return(0);
}


/*--------------------------------------------------*/
int gdfPrintHeader(FILE *fp, FSGD *gd) {
  int r;

  switch (gd->version) {
  case 1:
    r = gdfPrintV1(fp,gd);
    break;
  default:
    printf("ERROR: FSGDF version %d unsupported\n",gd->version);
    return(1);
  }
  return(r);
}


/*--------------------------------------------------*/
int gdfPrintStdout(FSGD *gd) {
  return gdfPrintHeader(stdout,gd);
}


/*--------------------------------------------------*/
static int gdfPrintV1(FILE *fp, FSGD *gd) {
  int n,m;

  if (gd->version != 1) {
    fprintf(fp,"ERROR: FSGDF version = %d, should be 1 \n",gd->version);
    return(1);
  }
  fprintf(fp,"GroupDescriptorFile 1\n");

  if (strlen(gd->title) > 0)
    fprintf(fp,"Title %s\n",gd->title);
  if (strlen(gd->measname) > 0)
    fprintf(fp,"MeasurementName %s\n",gd->measname);
  if (strlen(gd->tessellation) > 0)
    fprintf(fp,"Tessellation %s\n",gd->tessellation);
  if (strlen(gd->regsubj) > 0)
    fprintf(fp,"RegistrationSubject %s\n",gd->regsubj);
  if (strlen(gd->datafile) > 0)
    fprintf(fp,"PlotFile %s\n",gd->datafile);
  if (strlen(gd->DesignMatFile) > 0)
    fprintf(fp,"DesignMatFile %s %s\n",gd->DesignMatFile,gd->DesignMatMethod);
  if(!strcasecmp(gd->DesignMatMethod,"DODS") || !strcasecmp(gd->DesignMatMethod,"DOSS"))
    fprintf(fp,"%s\n",gd->DesignMatMethod);
  fprintf(fp,"DeMeanFlag %d\n",gd->DeMean);
  fprintf(fp,"ReScaleFlag %d\n",gd->ReScale);
  fprintf(fp,"ResidualFWHM %lf\n",gd->ResFWHM);
  fprintf(fp,"LogY %d\n",gd->LogY);
  if (strlen(gd->defvarlabel) > 0)
    fprintf(fp,"DefaultVariable %s\n",gd->defvarlabel);
  if (gd->nclasses > 0) {
    for (n=0; n < gd->nclasses; n++) {
      fprintf(fp,"Class %s",gd->classlabel[n]);
      if (strlen(gd->classmarker[n])>0)
        fprintf(fp," %s",gd->classmarker[n]);
      if (strlen(gd->classcolor[n])>0)
        fprintf(fp," %s",gd->classcolor[n]);
      fprintf(fp,"\n");
    }
  }

  for (n=0; n < gd->nvarsfromfile; n++)
    fprintf(fp,"# VariableFromFile %s %s %d %d\n",
            gd->tablefile[n],gd->varfield[n],
            gd->fieldcol[n],gd->datacol[n]);

  if (gd->nvariables > 0) {
    fprintf(fp,"Variables ");
    for (n=0; n < gd->nvariables; n++)
      fprintf(fp,"%s ",gd->varlabel[n]);
    fprintf(fp,"\n");
  }

  if(gd->nContrasts > 0) {
    for (n=0; n < gd->nContrasts; n++){
      if(!gd->IsFContrast[n]){
	fprintf(fp,"Contrast %s ",gd->ContrastName[n]);
	for(m=0; m < gd->C[n]->cols; m++)
	  fprintf(fp,"%g ",gd->C[n]->rptr[1][m+1]);
	fprintf(fp,"\n");
      }
      if(gd->IsFContrast[n]){
	fprintf(fp,"FContrast %s ",gd->ContrastName[n]);
	for(m=0; m < gd->FContrastNSub[n]; m++)
	  fprintf(fp,"%s ",gd->FContrastSub[n][m]);
	fprintf(fp,"\n");
	MatrixPrint(stdout,gd->C[n]);
      }
    }
  }

  if (gd->ninputs > 0) {
    for (n=0; n < gd->ninputs; n++) {
      //fprintf(fp,"Input %s %s (id=%d) ",gd->subjid[n],
      //        gd->classlabel[gd->subjclassno[n]],gd->subjclassno[n]+1);
      fprintf(fp,"Input %s %s ",
              gd->subjid[n],gd->classlabel[gd->subjclassno[n]]);
      if (gd->nvariables > 0) {
        for (m=0; m < gd->nvariables; m++)
          fprintf(fp,"%g ",gd->varvals[n][m]);
      }
      fprintf(fp,"\n");
    }
  }

  return(0);
}


/*--------------------------------------------------*/
FSGD *gdfRead(const char *gdfname, int LoadData) {
  FSGD *gd;
  FILE *fp;
  char tmpstr[1000];
  int version=0;
  int nv;
  MRI *mritmp;
  char *dirname, *basename;
  std::string datafilename;
  MATRIX *Xt,*XtX,*iXtX;

  printf("gdfRead(): reading %s\n",gdfname);

  nv = fio_FileHasCarriageReturn(gdfname);
  if(nv == -1) return(NULL);

  if(nv != 0){
    printf("\n");
    printf("WARNING: carriage returns have been detected in file %s\n",gdfname);
    printf("Was it created on a Windows computer?\n");
    printf("This may cause an error in reading the FSGD file.\n");
    printf("If so, try running:\n");
    printf("    cat %s | sed 's/\\r/\\n/g' > new.%s \n",gdfname,gdfname);
    printf("Then use new.%s \n",gdfname);
    printf("\n");
  }

  fp = fopen(gdfname,"r");
  if (fp==NULL) {
    printf("ERROR: gdfRead: cannot open %s for reading\n",gdfname);
    return(NULL);
  }

  fscanf(fp,"%s",tmpstr);
  if (strcasecmp(tmpstr,"GroupDescriptorFile") != 0) {
    printf("ERROR: gdfRead: %s is not formatted properly.\n",gdfname);
    printf("  The first string is '%s', should be 'GroupDescriptorFile'\n",
           tmpstr);
    return(NULL);
  }

  fscanf(fp,"%d",&version);
  fclose(fp);

  switch (version) {
  case 1:
    gd = gdfReadV1(gdfname);
    break;
  default:
    printf("ERROR: FSGDF version %d unsupported (%s) \n",version,gdfname);
    return(NULL);
  }
  if (gd == NULL) return(NULL);


  /* Extract the path from the gdf file. */
  dirname = (char*)fio_dirname(gdfname);
  if (NULL == dirname)
    printf("WARNING: Couldn't extract dirname from GDF header name %s.\n",
           gdfname);

  /* Load the design matrix, if there */
  if (strlen(gd->DesignMatFile) != 0) {
    /* Look for DesignMatFile first. If doesn't exist, prepend the
       directory from the gdf file. */
    datafilename = gd->DesignMatFile;
    if (!fio_FileExistsReadable(datafilename.c_str())) {
      datafilename = std::string(dirname) + "/" + std::string(gd->DesignMatFile);
      if (!fio_FileExistsReadable(datafilename.c_str())) {

        /* If that doesn't work, try the path from the GDF file and the
           base of the file name. */
        basename = fio_basename(gd->DesignMatFile,NULL);
	datafilename = std::string(dirname) + "/" + std::string(basename);
        free(basename);
      }

      if (!fio_FileExistsReadable(datafilename.c_str())) {
        printf("ERROR: gdfRead: could not find file %s\n",gd->DesignMatFile);
        return(NULL);
      }
    }
    gd->X = ReadMatlabFileVariable(datafilename.c_str(),"X");
    if (gd->X == NULL) {
      printf("ERROR: gdfRead: could not read variable X from %s\n",
             gd->DesignMatFile);
      return(NULL);
    }
    Xt = MatrixTranspose(gd->X,NULL);
    XtX = MatrixMultiply(Xt,gd->X,NULL);  // X'*X
    iXtX = MatrixInverse(XtX,NULL);       // inv(X'*X)
    gd->T = MatrixMultiply(iXtX,Xt,NULL); // T = inv(X'*X)*X'
    MatrixFree(&Xt);
    MatrixFree(&XtX);
    MatrixFree(&iXtX);

    if (stricmp(gd->DesignMatMethod,"none") == 0) {
      printf(
        "\n WARNING: the creation method of the design matrix is unknown.\n"
        " This is OK, however, you will not be able to view regression\n"
        " lines in the scatter plot viewer of tksurfer/tkmedit.\n\n");
    }
  }

  /* load the MRI containing our raw data. */
  if (LoadData && strlen(gd->datafile) > 0) {

    if (fio_FileExistsReadable(gd->datafile)) {
      datafilename = gd->datafile;
    }
    else {
      /* Construct the path of the data file by concat the
         path from the GDF file and the data file name */
      if (NULL != dirname) {
	datafilename = std::string(dirname) + "/" + std::string(gd->datafile);
      }

      /* If that doesn't work, try the path from the GDF file and the
         base of the file name. */
      if (!fio_FileExistsReadable(datafilename.c_str())) {
        basename = fio_basename(gd->datafile,NULL);
	datafilename = std::string(dirname) + "/" + std::string(basename);
        free(basename);
      }
    }

    gd->data = MRIread(datafilename.c_str());
    if (NULL == gd->data) {
      printf("ERROR: gdfRead: Couldn't read raw data at %s \n",gd->datafile);
      gdfFree(&gd);
      return(NULL);
    }
    nv = gd->data->width * gd->data->height * gd->data->depth;
    if (strcmp(gd->tessellation,"surface")==0 && nv != gd->data->width) {
      printf("INFO: gdfRead: reshaping\n");
      mritmp = mri_reshape(gd->data, nv, 1, 1, gd->data->nframes);
      MRIfree(&gd->data);
      gd->data = mritmp;
    }
  }

  if (gd->LogY) {
    printf("gdfRead(): Computing log of input\n");
    MRIlog(gd->data,NULL,0,0,gd->data);
  }

  if (NULL != dirname) free(dirname);

  return(gd);
}


/*--------------------------------------------------*/
static FSGD *gdfReadV1(const char *gdfname) {
  FSGD *gd;
  FSENV *env;
  FILE *fp;
  char *cp, tag[1000], tmpstr[1000], class_name[100];
  int version,r,n,m,k,err,ncols,c;
  double d;

  env = FSENVgetenv();

  fp = fopen(gdfname,"r");
  if (fp==NULL) {
    printf("ERROR: gdfReadV1: cannot open %s for reading\n",gdfname);
    return(NULL);
  }

  fscanf(fp,"%s",tag);
  if (strcasecmp(tag,"GroupDescriptorFile") != 0) {
    printf("ERROR: gdfReadV1: %s is not formated properly\n",gdfname);
    printf("  The first string is '%s', should be 'GroupDescriptorFile'\n",
           tmpstr);
    return(NULL);
  }

  fscanf(fp,"%d",&version);
  if (version != 1) {
    printf("ERROR: gdfReadV1: version=%d, != 1 (%s)\n",version,gdfname);
    return(NULL);
  }

  gd = gdfAlloc(1);
  gd->nvarsfromfile = 0;
  gd->DeMean = -1;
  gd->ReScale = 0;

  /*------- begin input loop --------------*/
  while (1) {

    r = fscanf(fp,"%s",tag);
    if (r==EOF) break;

    if (Gdiag_no > 0) printf("fsgd tag: %s\n",tag);

    if (!strcasecmp(tag,"Title")) {
      // account for possible whitespace in title text
      char c[2];c[0]=0;c[1]=0;
      while (c[0] != '\n') {
        r = fscanf(fp,"%c",&c[0]);
        if (r==EOF) goto formaterror;
        strcat(gd->title, c);
      }
      continue;
    }

    if (!strcasecmp(tag,"MeasurementName")) {
      r = fscanf(fp,"%s",gd->measname);
      if (r==EOF) goto formaterror;
      continue;
    }

    if (!strcasecmp(tag,"Tessellation")) {
      r = fscanf(fp,"%s",gd->tessellation);
      if (r==EOF) goto formaterror;
      continue;
    }

    if (!strcasecmp(tag,"RegistrationSubject")) {
      r = fscanf(fp,"%s",gd->regsubj);
      if (r==EOF) goto formaterror;
      continue;
    }

    if (!strcasecmp(tag,"PlotFile")) {
      r = fscanf(fp,"%s",gd->datafile);
      if (r==EOF) goto formaterror;
      continue;
    }
    /*----------------- DesignMat Line ---------------------*/
    if (!strcasecmp(tag,"DesignMatFile")) {
      r = fscanf(fp,"%s %s",gd->DesignMatFile,gd->DesignMatMethod);
      if (r==EOF) goto formaterror;
      continue;
    }
    /*----------------- ResidualFWHM Line ---------------------*/
    if (!strcasecmp(tag,"ResidualFWHM")) {
      r = fscanf(fp,"%lf",&gd->ResFWHM);
      if (r==EOF) goto formaterror;
      continue;
    }
    /*----------------- LogY Line ---------------------*/
    if (!strcasecmp(tag,"LogY")) {
      r = fscanf(fp,"%d",&gd->LogY);
      if (r==EOF) goto formaterror;
      continue;
    }
    /*----------------- DeMeanFlag ---------------------*/
    if (!strcasecmp(tag,"DeMeanFlag")) {
      r = fscanf(fp,"%d",&gd->DeMean);
      if (r==EOF) goto formaterror;
      continue;
    }
    /*----------------- ReScaleFlag ---------------------*/
    if (!strcasecmp(tag,"ReScaleFlag")) {
      r = fscanf(fp,"%d",&gd->ReScale);
      if (r==EOF) goto formaterror;
      continue;
    }
    /*----------------- DefaultVariable Line ---------------------*/
    if (!strcasecmp(tag,"DefaultVariable")) {
      r = fscanf(fp,"%s",gd->defvarlabel);
      if (r==EOF) goto formaterror;
      continue;
    }
    /*----------------- Class Line ---------------------*/
    if (!strcasecmp(tag,"Class")) {
      r = fscanf(fp,"%s",gd->classlabel[gd->nclasses]);
      if (r==EOF) goto formaterror;
      fgets(tmpstr,1000,fp);
      r = gdfCountItemsInString(tmpstr);
      if (r == 1) sscanf(tmpstr,"%s",gd->classmarker[gd->nclasses]);
      if (r == 2) sscanf(tmpstr,"%s %s",gd->classmarker[gd->nclasses],
                         gd->classcolor[gd->nclasses]);
      gd->nclasses ++;
      continue;
    }
    /*----------- Matrix Creation Method */
    if(!strcasecmp(tag,"DODS")){
      strcpy(gd->gd2mtx_method,"DODS");
      strcpy(gd->DesignMatMethod,"DODS");
      continue;
    }
    if(!strcasecmp(tag,"DOSS")){
      strcpy(gd->gd2mtx_method,"DOSS");
      strcpy(gd->DesignMatMethod,"DOSS");
      continue;
    }
    /*---------------- Contrast ------------------------*/
    if(!strcasecmp(tag,"Contrast")){
      r = fscanf(fp,"%s",tmpstr);
      if (r==EOF) goto formaterror;
      gd->ContrastName[gd->nContrasts] = strcpyalloc(tmpstr);
      gd->IsFContrast[gd->nContrasts] = 0;
      fgets(tmpstr,1000,fp);
      r = gdfCountItemsInString(tmpstr);
      if(r<1) {
	printf("ERROR: Contrast, not enough items, %s\n",tmpstr);
        goto formaterror;
      }
      gd->C[gd->nContrasts] = MatrixAlloc(1,r,MATRIX_REAL);
      for(n=0; n < r; n++){
	cp = gdfGetNthItemFromString(tmpstr, n);
	sscanf(cp,"%f",&(gd->C[gd->nContrasts]->rptr[1][n+1]));
	free(cp);
      }
      gd->nContrasts++;
      continue;
    }
    /*---------------- FContrast ------------------------*/
    if(!strcasecmp(tag,"FContrast")){
      r = fscanf(fp,"%s",tmpstr);
      if (r==EOF) goto formaterror;
      gd->IsFContrast[gd->nContrasts] = 1;
      gd->ContrastName[gd->nContrasts] = strcpyalloc(tmpstr);
      fgets(tmpstr,1000,fp);
      r = gdfCountItemsInString(tmpstr);
      if(r<0) {
	printf("ERROR: FContrast, not enough items, %s\n",tmpstr);
        goto formaterror;
      }
      gd->FContrastNSub[gd->nContrasts] = r;
      gd->FContrastSub[gd->nContrasts] = (char **)calloc(r,sizeof(char*));
      for(n=0; n < r; n++){
	cp = gdfGetNthItemFromString(tmpstr, n);
	gd->FContrastSub[gd->nContrasts][n] = strcpyalloc(cp);
	free(cp);
      }
      gd->nContrasts++;
      continue;
    }
    /*----------------- Variables Line ---------------------*/
    if (!strcasecmp(tag,"Variables") || !strcasecmp(tag,"Variable")) {
      if (gd->nvariables != 0) {
        printf("ERROR: gdfReadV1: multiple 'Variables' lines found\n");
        goto formaterror;
      }
      r = gdfCountItemsOnLine(fp);
      if (r==0) {
        fprintf(
          stderr,
          "WARNING: gdfReadV1: no variables on 'Variables' line found\n");
        continue;
      }
      for (m=0; m < r; m++){
        fscanf(fp,"%s",gd->varlabel[m]);
	if(strcasecmp(gd->varlabel[m],"sex")==0 || strcasecmp(gd->varlabel[m],"gender")==0){
	  printf("WARNING: variable %d is \"%s\" which is often a discrete factor\n",m,gd->varlabel[m]);
	  printf("  The proper way to handle discrete factors is to create classes.\n");
	  printf("  See https://surfer.nmr.mgh.harvard.edu/fswiki/FsgdExamples\n");
	}
      }
      gd->nvariables = r;
      r = gdfCheckVarRep(gd);
      if (r != -1) {
        printf("ERROR: gdfReadV1: variable label %s appears multiple times\n",
               gd->varlabel[r]);
        goto formaterror;
      }
      continue;
    }
    /*----------------- VariableFromFile Line ---------------------*/
    if (!strcasecmp(tag,"VariableFromFile")) {
      // Eg, VariableFromFile stats/aseg.stats Left-Hippocampus 5 4
      // 5 = field name column
      // 4 = data column
      m = gd->nvarsfromfile;
      fscanf(fp,"%s",tmpstr);
      gd->tablefile[m] = strcpyalloc(tmpstr);
      fscanf(fp,"%s",tmpstr);
      gd->varfield[m] = strcpyalloc(tmpstr);
      fscanf(fp,"%d %d",&gd->fieldcol[m],&gd->datacol[m]);
      gd->nvarsfromfile ++;
    }
    /*----------------- VariableFromASeg Line ---------------------*/
    if (!strcasecmp(tag,"VariableFromASeg")) {
      m = gd->nvarsfromfile;
      gd->tablefile[m] = strcpyalloc("stats/aseg.stats");
      fscanf(fp,"%s",tmpstr);
      gd->varfield[m] = strcpyalloc(tmpstr);
      gd->fieldcol[m] = 5;
      gd->datacol[m]  = 4;
      gd->nvarsfromfile ++;
    }
    /*----------------- Input Line ---------------------*/
    if (!strcasecmp(tag,"Input")) {
      if (gd->ninputs > FSGDF_NINPUTS_MAX) {
        printf("ERROR: gdfReadV1: the number of inputs in FSGD file "
               "exceeds the maximum allowed %d\n",
               FSGDF_NINPUTS_MAX);
        return(NULL);
      }
      if (gd->nclasses == 0) {
        printf("FSGDF Format Error: no classes defined before "
               "the first input line.\n");
        return(NULL);
      }
      n = gd->ninputs; /* current input number */
      r = fscanf(fp,"%s %s",gd->subjid[n],class_name);
      if (r==EOF) {
        printf("Input line %d: ",n+1);
        goto formaterror;
      }
      r = gdfClassNo(gd,class_name);
      if (r < 0) {
        printf("Input line %d, subjid = %s, class %s not defined \n",
               n+1,gd->subjid[n],class_name);
        goto formaterror;
      }
      gd->subjclassno[n] = r;

      r = gdfCountItemsOnLine(fp);
      if (r != gd->nvariables) {
        printf("ERROR: gdfReadV1: Input line %d, subjid = %s\n",
               n+1,gd->subjid[n]);
        printf("       Found %d variables, expected. %d \n",r,gd->nvariables);
        //printf("%s\n",tmpstr);
        goto formaterror;
      }
      for (m=0; m < gd->nvariables; m++) {
	char tmpstr[1000];
	fscanf(fp,"%s",tmpstr);
	if(isalpha(tmpstr[0])){
	  printf("ERROR: gdfReadV1: Format Error: Input line %d, subjid = %s\n",n+1,gd->subjid[n]);
	  printf(" Variable %d has character string %s\n",m+1,tmpstr);
	  printf(" Variables should be continuous numbers\n");
	  goto formaterror;
	}
	sscanf(tmpstr,"%f",&gd->varvals[n][m]);
      }
      for (m=0; m < gd->nvarsfromfile; m++) {
        sprintf(tmpstr,"%s/%s/%s",
                env->SUBJECTS_DIR,gd->subjid[n],gd->tablefile[m]);
        err = gdfGetDDataFromTable(tmpstr,
                                   gd->varfield[m],
                                   gd->fieldcol[m],
                                   gd->datacol[m],&d);
        if (err) {
          gdfFree(&gd);
          return(NULL);
        }
        gd->varvals[n][m+gd->nvariables] = d;
      }

      gd->ninputs ++;
      continue;
    }

    if (tag[0] != '#') fprintf(stderr,"INFO: ignoring tag %s \n",tag);
    fgets(tmpstr,1000,fp);

  }/*------- End loop over tags ----------------*/

  for (m=0; m < gd->nvarsfromfile; m++)
    sprintf(gd->varlabel[m+gd->nvariables],"%s",gd->varfield[m]);
  gd->nvariables += gd->nvarsfromfile;

  r = gdfCheckClassRep(gd);
  if (r != -1) {
    printf("ERROR: gdfReadV1: class label %s appears multiple times\n",
           gd->classlabel[r]);
    sprintf(tag,"Class");
    goto formaterror;
  }

  r = gdfCheckAllClassesUsed(gd);
  if (r != -1){
    printf("ERROR: gdfReadV1: class %s is defined but not used.\n",gd->classlabel[r]);
    goto formaterror;
  }

  r = gdfCheckSubjRep(gd);
  if (r != -1) {
    /* See gdfCheckSubjRep() for fsgdf_AllowSubjRep usage */
    printf("ERROR: gdfReadV1: subject id %s appears multiple times\n",
           gd->subjid[r]);
    sprintf(tag,"Input");
    goto formaterror;
  }

  r = gdfGetDefVarLabelNo(gd);
  if (r == -1) {
    printf("ERROR: gdfReadV1: default variable %s does not exist in list\n",gd->defvarlabel);
    sprintf(tag,"DefaultVariable");
    goto formaterror;
  }

  // Convert FContrast spec to contrast matrix
  for(n=0; n < gd->nContrasts; n++){
    if(!gd->IsFContrast[n]) continue;
    ncols = -1;
    for(m = 0; m < gd->FContrastNSub[n]; m++){
      err = 1;
      for(k=0; k < gd->nContrasts; k++){
	if(! strcmp(gd->FContrastSub[n][m],gd->ContrastName[k])){
	  if(gd->IsFContrast[k]){
	    printf("ERROR: FContrast %s references another FContrast %s\n",
		   gd->ContrastName[n],gd->FContrastSub[n][m]);
	    goto formaterror;
	  }
	  if(ncols == -1) ncols = gd->C[k]->cols;
	  if(gd->C[k]->cols != ncols){
	    printf("ERROR: Contrasts have conflicting numbers of columns\n");
	    goto formaterror;
	  }
	  err = 0;
	  break;
	}
      }
      if(err){
	printf("ERROR: cannot find contrast %s needed for FContrast %s\n",
	       gd->FContrastSub[n][m],gd->ContrastName[n]);
	goto formaterror;
      }
    }
    gd->C[n] = MatrixAlloc(gd->FContrastNSub[n],ncols,MATRIX_REAL);
    for(m = 0; m < gd->FContrastNSub[n]; m++){
      for(k=0; k < gd->nContrasts; k++){
	if(strcmp(gd->FContrastSub[n][m],gd->ContrastName[k])) continue;
	for(c=0; c < ncols; c++)
	  gd->C[n]->rptr[m+1][c+1] = gd->C[k]->rptr[1][c+1];
      }
    }
  }

  if (gd->DeMean == -1 && gd->nvariables > 0) {
    printf("INFO: DeMeanFlag keyword not found,"
           " DeMeaning will NOT be done.\n");
    gd->DeMean = 0;
  } else {
    if (gd->nvariables > 0) {
      if (gd->DeMean) printf("INFO: demeaning continuous variables\n");
      else           printf("INFO: NOT demeaning continuous variables\n");
    }
  }

  gdfVarMeans(gd);
  gdfClassVarMeans(gd);

  r = gdfCheckNPerClass(gd);
  if(r) return(NULL);

  return(gd);

 formaterror:
  printf("FSGDF Format Error: file = %s, tag=%s\n",gdfname,tag);
  gdfFree(&gd);
  FSENVfree(&env);
  return(NULL);
}


/*--------------------------------------------------
  gdfGetDataMRIHeader() - returns the MRI header
  information for the data in the given FSGD header
  file. This is only the header info, not the data.
  --------------------------------------------------*/
MRI *gdfReadDataInfo(const char *gdfname) {
  FSGD *gd=NULL;
  MRI *info=NULL;

  /* Read this header file but don't load the data. */
  gd = gdfRead(gdfname, 0);
  if (NULL==gd) {
    printf("ERROR: gdfReadDataInfo: Couldn't read GDF %s\n",gdfname);
    return(NULL);
  }

  /* Now try and read an MRI struct from the datafile file name we
     got. This doesn't load the data, just the header info. */
  info = MRIreadInfo(gd->datafile);
  if (NULL==info) {
    printf("ERROR: gdfReadDataInfo: Couldn't read MRI %s\n",gd->datafile);
    gdfFree(&gd);
    return(NULL);
  }

  /* Free the GDF info. */
  gdfFree(&gd);

  /* Return the MRI header we got. */
  return info;
}


/*--------------------------------------------------
  gdfClassNo() - returns the zero-based class number
  associated with a class label.
  --------------------------------------------------*/
int gdfClassNo(FSGD *gd, char *class_number) {
  int code;
  for (code=0; code < gd->nclasses; code++)
    if (!strcmp(gd->classlabel[code],class_number)) return(code);
  return(-1);
}


/*--------------------------------------------------
  gdfCountItemsInString() returns the number of items
  in the given string, where an item is defined as
  one or more contiguous non-blank characters.
  --------------------------------------------------*/
int gdfCountItemsInString(const char *str) {
  int len, n, nhits;

  len = strlen(str);

  nhits = 0;
  n = 0;
  while (n < len) {
    while (isblank(str[n])) n++;
    if (n >= len) break;
    if (str[n] == '\0' || str[n] == '\n' || str[n] == '\r') break;
    while (!isblank(str[n])) n++;
    nhits++;
  }

  //printf("nhits %d\n",nhits);

  return(nhits);
}


/*-------------------------------------------------------------------
  gdfGetNthItemFromString() - extracts the nth item from a string.
  An item is defined as one or more non-white space chars. If nth
  is -1, then it returns the last item. item is a string that
  must be freed by the caller.
  ------------------------------------------------------------------*/
char *gdfGetNthItemFromString(const char *str, const int nth) {
  char *item;
  int nitems;

  nitems = gdfCountItemsInString(str);
  if (nth >= nitems) {
    printf("ERROR: asking for item %d, only %d items in string\n",nth,nitems);
    printf("%s\n",str);
    return(NULL);
  }

  const std::string src(str);
  std::vector<std::string> items;
  std::string tmpstr;
  bool inItem = !std::isspace(src.at(0));
  for(auto it=src.begin(); it!=src.end(); ++it ) {
    if( std::isspace(*it) ) {
      if( inItem ) {
	// We've just completed the next item
	items.push_back(tmpstr);
	tmpstr.clear();
      } else {
	// Nothing to do; we're just consuming blanks
      }
      inItem = false;
    } else {
      inItem = true;
      // We're inside an item, so accumulate
      tmpstr.push_back(*it);
    }
  }

  if( items.size() != static_cast<size_t>(nitems) ) {
    std::cerr << __FUNCTION__
	      << ": Length of items vector did not match nitems"
	      << std::endl;
    std::cerr << "str: '" << str << std::endl;
    std::cerr << "items: ";
    for( auto it=items.begin(); it!=items.end(); ++it ) {
      std::cerr << (*it) << " -|- ";
    }
    std::cerr << std::endl;
    throw std::logic_error("Incorrect item count");
  }
  
  if( nth < 0 ) {
    item = strcpyalloc(items.back().c_str());
  } else {
    item = strcpyalloc(items.at(nth).c_str());
  }
  return(item);
}


/*--------------------------------------------------
  gdfCountItemsOnLine() returns the number of items
  before the next newline, where an item is defined as
  one or more contiguous non-blank characters. The
  file pointer is not changed.
  --------------------------------------------------*/
int gdfCountItemsOnLine(FILE *fp) {
  fpos_t now;
  char tmpstr[10000];
  int nitems;

  fgetpos(fp,&now);
  fgets(tmpstr,10000,fp);
  fsetpos(fp,&now);

  nitems = gdfCountItemsInString(tmpstr);

  return(nitems);
}


/*--------------------------------------------------
  gdfCheckClassRep() - checks whether there are any
  repetitions in the class labels. Returns 0 if
  there are no reps.
  --------------------------------------------------*/
static int gdfCheckClassRep(FSGD *gd) {
  int n,m;

  for (n=0; n < gd->nclasses; n++)
    for (m=n+1; m < gd->nclasses; m++)
      if (strcmp(gd->classlabel[n],gd->classlabel[m])==0)
        return(n);
  return(-1);
}


/*--------------------------------------------------
  gdfCheckVarRep() - checks whether there are any
  repetitions in the variable labels. Returns -1 if
  there are no reps, otherwise returns the index
  of one of the reps.
  --------------------------------------------------*/
static int gdfCheckVarRep(FSGD *gd) {
  int n,m;

  for (n=0; n < gd->nvariables; n++)
    for (m=n+1; m < gd->nvariables; m++)
      if (strcmp(gd->varlabel[n],gd->varlabel[m])==0)
        return(n);
  return(-1);
}


/*--------------------------------------------------
  gdfCheckSubjRep() - checks whether there are any
  repetitions in the subject labels. Returns -1 if
  there are no reps, otherwise returns the index
  of one of the reps. If the global variable
  fsgdf_AllowSubjRep is set to 1, then is always
  returns -1 (ie, no reps).
  --------------------------------------------------*/
static int gdfCheckSubjRep(FSGD *gd) {
  int n,m;
  extern int fsgdf_AllowSubjRep;

  if (fsgdf_AllowSubjRep) return(-1);

  if (getenv("FSGDF_ALLOW_SUBJ_REP") != NULL) return(-1);

  for (n=0; n < gd->ninputs; n++)
    for (m=n+1; m < gd->ninputs; m++)
      if (strcmp(gd->subjid[n],gd->subjid[m])==0)
        return(n);
  return(-1);
}


/*--------------------------------------------------
  gdfCheckAllClassesUsed() - checks whether all classes
  defined are present in the subject data. Returns -1
  if everything is ok, otherwise returns the class
  number that is not used.
  --------------------------------------------------*/
static int gdfCheckAllClassesUsed(FSGD *gd) {
  int n,m,ok;

  for (n=0; n < gd->nclasses; n++) {
    ok = 0;
    for (m=0; m < gd->ninputs; m++)
      if (gd->subjclassno[m] == n) ok = 1;
    if (!ok) return(n);
  }

  return(-1);
}

/*!
  \fn int gdfCheckNPerClass(FSGD *gd);
  \brief Checks to make sure that each class has the minimum number
  of members given the mtxmethod and number of variables
 */
int gdfCheckNPerClass(FSGD *gd)
{
  int cno;
  if(stricmp(gd->gd2mtx_method,"doss") == 0) return(0);

  for (cno = 0; cno < gd->nclasses; cno++) {
    //printf("Class %s has %d members\n",gd->classlabel[cno],(int)gd->NPerClass[cno]);
    if(gd->NPerClass[cno] < gd->nvariables+1){
      printf("ERROR: Class %s has %d members. With %d variables and using DODS, "
	     "you need at least %d members\n",gd->classlabel[cno],(int)gd->NPerClass[cno],
	     gd->nvariables,gd->nvariables+1);
      return(1);
    }
    if(gd->NPerClass[cno] == gd->nvariables+1){
      printf("WARNING: Class %s has %d members. With %d variables and using DODS, "
	     "This is the bare minimum which may cause problems with the design matrix.\n",
	     gd->classlabel[cno],(int)gd->NPerClass[cno],gd->nvariables);
    }
  }
  return(0);
}


/*--------------------------------------------------
  gdfGetDefVarLabelNo() - returns the label number
  of the default variable. If the default variable
  is not defined, 0 is returned. If there is no
  match, -1 is returned. Otherwise, returns the
  index.
  --------------------------------------------------*/
static int gdfGetDefVarLabelNo(FSGD *gd) {
  int n;

  if (strlen(gd->defvarlabel) == 0) return(0);

  for (n=0; n < gd->nvariables; n++)
    if (strcmp(gd->varlabel[n],gd->defvarlabel)==0)
      return(n);
  return(-1);
}


/*--------------------------------------------------
  gdfGetVarLabelNo() - returns the label number of the given
  variable. If there is no match, -1 is returned. Otherwise,
  returns the index.
  --------------------------------------------------*/
int gdfGetVarLabelNo(FSGD *gd, char *LabelName) {
  int n;

  for (n=0; n < gd->nvariables; n++)
    if (strcmp(gd->varlabel[n],LabelName)==0)
      return(n);
  return(-1);
}


/*--------------------------------------------------------*/
int gdfVarMeans(FSGD *gd) {
  int vno, n;

  // Init
  for (vno = 0; vno < gd->nvariables; vno++){
    gd->VarMeans[vno] = 0;
    gd->VarStds[vno] = 0;
  }

  // Sum over all inputs regardless of class
  for (n=0; n < gd->ninputs; n++) {
    for (vno = 0; vno < gd->nvariables; vno++) {
      gd->VarMeans[vno] += gd->varvals[n][vno];
    }
  }

  // Divide Sum by ninputs
  for (vno = 0; vno < gd->nvariables; vno++)
    gd->VarMeans[vno] /= gd->ninputs;

  // Computes SumSq
  for (n=0; n < gd->ninputs; n++) {
    for (vno = 0; vno < gd->nvariables; vno++) {
      gd->VarStds[vno] += SQR(gd->varvals[n][vno]-gd->VarMeans[vno]);
    }
  }
  // Compute StdDev 
  for (vno = 0; vno < gd->nvariables; vno++){
    gd->VarStds[vno] = sqrt(gd->VarStds[vno]/gd->ninputs);
    if(gd->VarStds[vno] < FLT_MIN) gd->VarStds[vno] = FLT_MIN;
  }

  if (gd->nvariables > 0) {
    printf("Continuous Variable Means (all subjects)\n");
  }
  for (vno = 0; vno < gd->nvariables; vno++)
    printf("%d %s %g %g\n",vno,gd->varlabel[vno],gd->VarMeans[vno],gd->VarStds[vno]);

  return(0);
}


/*--------------------------------------------------------*/
int gdfClassVarMeans(FSGD *gd) {
  int cno, vno, n;

  // Init
  for (cno = 0; cno < gd->nclasses; cno++) {
    gd->NPerClass[cno] = 0;
    for (vno = 0; vno < gd->nvariables; vno++)
      gd->ClassVarMeans[cno][vno] = 0;
  }

  // Sum and count nper class
  for (n=0; n < gd->ninputs; n++) {
    cno = gd->subjclassno[n];
    gd->NPerClass[cno] ++;
    for (vno = 0; vno < gd->nvariables; vno++) {
      gd->ClassVarMeans[cno][vno] += gd->varvals[n][vno];
    }
  }

  // Divide by nper class
  for (cno = 0; cno < gd->nclasses; cno++) {
    for (vno = 0; vno < gd->nvariables; vno++)
      gd->ClassVarMeans[cno][vno] /= gd->NPerClass[cno];
  }

  if (gd->nvariables > 0) {
    printf("Class Size and Means of each Continuous Variable\n");
    for (cno = 0; cno < gd->nclasses; cno++) {
      printf("%d %s %2d ",cno+1,gd->classlabel[cno],(int)gd->NPerClass[cno]);
      for (vno = 0; vno < gd->nvariables; vno++)
        printf("%8.4f ",gd->ClassVarMeans[cno][vno]);
      printf("\n");
    }
  }

  return(0);
}


/*--------------------------------------------------------
  gdfMatrixDOSS() - creates a design matrix that models each
  class as having a Different Offset, but models having the
  Same Slope
  ---------------------------------------------------------*/
MATRIX *gdfMatrixDOSS(FSGD *gd, MATRIX *X) {
  int nrows, ncols, r,v,c;
  double mn;

  nrows = gd->ninputs;
  ncols = gd->nclasses + gd->nvariables;

  X = MatrixZero(nrows,ncols,X);
  if (X==NULL) return(NULL);

  for (r=0; r<nrows; r++) {

    c = gd->subjclassno[r];
    X->rptr[r+1][c+1] = 1;

    for (v = 0; v < gd->nvariables; v++) {
      if (gd->DeMean) mn = gd->VarMeans[v];
      else           mn = 0;
      c = v + gd->nclasses;
      X->rptr[r+1][c+1] = gd->varvals[r][v] - mn;
      if(gd->ReScale) X->rptr[r+1][c+1] /= gd->VarStds[v];
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
MATRIX *gdfMatrixDODS(FSGD *gd, MATRIX *X) {
  int nrows, ncols, n,r,v,c;
  double mn;

  nrows = gd->ninputs;
  ncols = gd->nclasses * ( gd->nvariables + 1);

  X = MatrixZero(nrows,ncols,X);
  if (X==NULL) return(NULL);

  for (r=0; r < nrows; r++) {

    n = gd->subjclassno[r]; // 0-based class number
    c = n;
    X->rptr[r+1][c+1] = 1;

    for (v = 0; v < gd->nvariables; v++) {
      // Demean without regard to class
      //if (gd->DeMean) mn = gd->ClassVarMeans[n][v];
      if (gd->DeMean) mn = gd->VarMeans[v];
      else            mn = 0;
      c += gd->nclasses;
      X->rptr[r+1][c+1] = gd->varvals[r][v] - mn;
      if(gd->ReScale) X->rptr[r+1][c+1] /= gd->VarStds[v];
    }
  }

  return(X);
}


/*---------------------------------------------------*/
int gdfCheckMatrixMethod(const char *gd2mtx_method) {
  if ( stricmp(gd2mtx_method,"doss") == 0 ||
       stricmp(gd2mtx_method,"dods") == 0 ||
       stricmp(gd2mtx_method,"none") == 0
    ) return(0);

  printf("ERROR: gd2mtx method %s unrecoginzied.\n",gd2mtx_method);
  printf("       Legal values are dods, doss, and none.\n");
  return(1);
}


/*---------------------------------------------------*/
MATRIX *gdfMatrix(FSGD *gd, const char *gd2mtx_method, MATRIX *X) {
  if (gdfCheckMatrixMethod(gd2mtx_method)) return(NULL);

  if (stricmp(gd2mtx_method,"none") == 0) {
    printf("ERROR: gdfMatrix: cannot create matrix when method is none\n");
    return(NULL);
  }

  if (stricmp(gd2mtx_method,"doss") == 0)
    X = gdfMatrixDOSS(gd,X);
  if (stricmp(gd2mtx_method,"dods") == 0)
    X = gdfMatrixDODS(gd,X);

  gd->X = X;

  return(X);
}


/*------------------------------------------------------------
  gdfOffsetSlope() - computes the offset and slope regression
  parameters for the given class and variable numbers of
  the given voxel/vertex. The class and variable numbers
  are zero-based. Returns 0 if successfull, 1 otherwise.

  Note: this was desgined to be
  used with tksurfer/tkmedit on a point-and-click basis.
  It will be very slow if you try to compute the regression
  parameters for an entire volume/surface.
  ------------------------------------------------------------*/
int gdfOffsetSlope(FSGD *gd, int classno, int varno,
                   int c, int r, int s,
                   float *offset, float *slope) {
  MATRIX *y, *b;
  int n,nf;
  int nslope=0;

  if (strlen(gd->DesignMatMethod) == 0 ||
      stricmp(gd->DesignMatMethod,"none") == 0) {
    printf("ERROR: gdfOffsetSlope: cannot determine the offset "
           "and slope for the \n"
           "given group descriptor because the design matrix \n"
           "creation method is unknown\n");
    return(1);
  }

  if (classno >= gd->nclasses) {
    if (gd->nclasses) {
      printf("ERROR: gdfOffsetSlope: class number %d exceeds max %d\n",
             classno,gd->nclasses);
    }
    return(1);
  }
  if (varno >= gd->nvariables)  {
    if (gd->nvariables) {
      printf("ERROR: gdfOffsetSlope: variable number %d exceeds max %d\n",
             varno,gd->nvariables);
    }
    return(1);
  }
  if (c < 0 || c >= gd->data->width ||
      r < 0 || r >= gd->data->height ||
      s < 0 || s >= gd->data->depth) {
    printf("ERROR: gdfOffsetSlope: index exceeds data matrix dimension\n");
    return(1);
  }
  if (gd->T->cols != gd->data->nframes) {
    printf("ERROR: gdfOffsetSlope: dimension mismatch.\n");
    return(1);
  }
  nf = gd->T->cols;

  y = MatrixAlloc(nf,1,MATRIX_REAL);
  for (n=0; n<nf; n++) y->rptr[n+1][1] = MRIFseq_vox(gd->data,c,r,s,n);

  b = MatrixMultiply(gd->T,y,NULL);
  if (Gdiag) {
    printf("c=%d, r=%d, s=%d\n",c,r,s);
    printf("b=\n=========================\n");
    MatrixPrint(stdout,b);
    printf("=========================\n");
  }

  *offset = b->rptr[classno+1][1];

  if (stricmp(gd->DesignMatMethod,"doss") == 0)
    nslope = gd->nclasses + varno;
  if (stricmp(gd->DesignMatMethod,"dods") == 0)
    nslope = classno + ((varno+1) * gd->nclasses);

  if (Gdiag) {
    printf("nc = %d, nv = %d, method = %s, classno=%d, varno = %d, n=%d\n",
           gd->nclasses, gd->nvariables, gd->DesignMatMethod,
           classno,varno,nslope);
  }

  *slope = b->rptr[nslope+1][1];

  if (Gdiag) {
    printf("offset = %g, slope = %g\n",*offset,*slope);
  }

  MatrixFree(&y);
  MatrixFree(&b);

  return(0);
}


/*------------------------------------------------------------
  gdfGetTitle() - copies the title into the output argument.
  ------------------------------------------------------------*/
int gdfGetTitle(FSGD *gd, char *title) {
  if (NULL == gd)
    return(-1);

  strcpy(title,gd->title);

  return(0);
}


/*------------------------------------------------------------
  gdfGetMeasurementName() - copies the measurment name into
  the output argument.
  ------------------------------------------------------------*/
int gdfGetMeasurementName(FSGD *gd, char *name) {
  if (NULL == gd)
    return(-1);

  strcpy(name,gd->measname);

  return(0);
}


/*------------------------------------------------------------
  gdfGetFWHM() - returns ResFWHM.
  ------------------------------------------------------------*/
double gdfGetFWHM(FSGD *gd) {
  if (NULL == gd)  return(-1);
  return(gd->ResFWHM);
}


/*------------------------------------------------------------
  gdfGetLogY() - returns LogY flag
  ------------------------------------------------------------*/
int gdfGetlogY(FSGD *gd) {
  if (NULL == gd)  return(-1);
  return(gd->LogY);
}


/*------------------------------------------------------------
  gdfGetSubjectName() - copies the subject name into the
  output argument.
  ------------------------------------------------------------*/
int gdfGetSubjectName(FSGD *gd, char *name) {
  if (NULL == gd)
    return(-1);

  strcpy(name,gd->regsubj);

  return(0);
}


/*------------------------------------------------------------
  gdfGetDataFileName() - copies the data file name into the
  output argument.
  ------------------------------------------------------------*/
int gdfGetDataFileName(FSGD *gd, char *filename) {
  if (NULL == gd)
    return(-1);

  strcpy(filename,gd->datafile);

  return(0);
}


/*------------------------------------------------------------
  gdfGetNumClasses() - returns the number of classes in the
  output argument.
  ------------------------------------------------------------*/
int gdfGetNumClasses(FSGD *gd, int *nclasses) {
  if (NULL == gd)
    return(-1);

  *nclasses = gd->nclasses;

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthClassLabel() - copies the nth class label into the
  output argument where nclass is from 0 -> nclasses.
  ------------------------------------------------------------*/
int gdfGetNthClassLabel(FSGD *gd, int nclass, char *label) {
  if (NULL == gd)
    return(-1);
  if (nclass < 0 || nclass >= gd->nclasses)
    return(-1);

  strcpy(label,gd->classlabel[nclass]);

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthClassMarker() - copies the nth class marker into the
  output argument where nclass is from 0 -> nclasses.
  ------------------------------------------------------------*/
int gdfGetNthClassMarker(FSGD *gd, int nclass, char *marker) {
  if (NULL == gd)
    return(-1);
  if (nclass < 0 || nclass >= gd->nclasses)
    return(-1);

  strcpy(marker,gd->classmarker[nclass]);

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthClassColor() - copies the nth class color into the
  output argument where nclass is from 0 -> nclasses.
  ------------------------------------------------------------*/
int gdfGetNthClassColor(FSGD *gd, int nclass, char *color) {
  if (NULL == gd)
    return(-1);
  if (nclass < 0 || nclass >= gd->nclasses)
    return(-1);

  strcpy(color,gd->classcolor[nclass]);

  return(0);
}


/*------------------------------------------------------------
  gdfGetNumVariables() - returns the number of variables in the output
  argument. Special case: if nvariables is actually 0, we really do
  have a 'variable,' it's just the subject index. So the graphing code
  doesn't have to know about this, return 1.
  ------------------------------------------------------------*/
int gdfGetNumVariables(FSGD *gd, int *nvariables) {
  if (NULL == gd)
    return(-1);

  if (gd->nvariables == 0) {
    *nvariables = 1;
  } else {
    *nvariables = gd->nvariables;
  }

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthVariableLabel() - copies the nth variable label into the
  output argument. Special case: if 0 variables, copy Subject as the
  label.
  ------------------------------------------------------------*/
int gdfGetNthVariableLabel(FSGD *gd, int nvariable, char *label) {
  if (NULL == gd)
    return(-1);
  if (nvariable < 0 ||
      (gd->nvariables != 0 && nvariable >= gd->nvariables) )
    return(-1);

  if (gd->nvariables == 0) {
    strcpy(label,"Subject");
  } else {
    strcpy(label,gd->varlabel[nvariable]);
  }

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthVariableDefault() - copies of name of
  the default variable to the output argument.
  ------------------------------------------------------------*/
int gdfGetDefaultVariable(FSGD *gd, char *label) {
  if (NULL == gd)
    return(-1);

  if (gd->nvariables == 0) {
    strcpy(label,"Subject");
  } else {
    strcpy(label,gd->defvarlabel);
  }

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthVariableDefaultIndex() - copies of index of
  the default variable to the output parameter
  ------------------------------------------------------------*/
int gdfGetDefaultVariableIndex(FSGD *gd, int *nvariable) {
  if (NULL == gd)
    return(-1);

  *nvariable = gdfGetDefVarLabelNo(gd);

  return(0);
}


/*------------------------------------------------------------
  gdfGetNumSubjects() - returns the number of subjects in the
  output parameter.
  ------------------------------------------------------------*/
int gdfGetNumSubjects(FSGD *gd, int *nsubjects) {
  if (NULL == gd)
    return(-1);

  *nsubjects = gd->ninputs;

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthSubjectID() - copies the id of the nth subject into
  the ouput parameter where nsubject is from 0 -> ninputs.
  ------------------------------------------------------------*/
int gdfGetNthSubjectID(FSGD *gd, int nsubject, char *id) {
  if (NULL == gd)
    return(-1);
  if (nsubject < 0 || nsubject >= gd->ninputs)
    return(-1);

  strcpy(id,gd->subjid[nsubject]);

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthSubjectClass() - returns the index of the nth subject's
  class in the output parameter where nsubject is from 0 -> ninputs.
  ------------------------------------------------------------*/
int gdfGetNthSubjectClass(FSGD *gd, int nsubject, int *class_number) {
  if (NULL == gd)
    return(-1);
  if (nsubject < 0 || nsubject >= gd->ninputs)
    return(-1);

  *class_number = gd->subjclassno[nsubject];

  return(0);
}


/*------------------------------------------------------------
  gdfGetNthSubjectNthValue() - returns the index of the nth subject's
  nth variable value where nsubject is from 0 -> ninputs and nvariable
  is from 0 -> nvariables. Special case: if 0 variables, return the
  subject index.
  ------------------------------------------------------------*/
int gdfGetNthSubjectNthValue(FSGD *gd, int nsubject,
                             int nvariable, float *value) {
  if (NULL == gd)
    return(-1);
  if (nsubject < 0 || nsubject >= gd->ninputs)
    return(-1);
  if (nvariable < 0 ||
      (gd->nvariables != 0 && nvariable >= gd->nvariables) )
    return(-1);

  if (gd->nvariables == 0) {
    *value = nsubject;
  } else {
    *value = gd->varvals[nsubject][nvariable];
  }

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
                                int x, int y, int z, float *value) {
  float v;
  int errs=0;

  if (NULL == gd)
    return(-1);
  if (nsubject < 0 || nsubject >= gd->ninputs)
    return(-1);

  // bounds-checks
  if (x < 0) {
    printf("ERROR: gdfGetNthSubjectMeasurement: x=%d < 0\n", x);
    errs++;
  }
  if (x > gd->data->width) {
    printf("ERROR: gdfGetNthSubjectMeasurement: x=%d > gd->data->width=%d\n",
           x, gd->data->width);
    errs++;
  }
  if (y < 0) {
    printf("ERROR: gdfGetNthSubjectMeasurement: y=%d < 0\n", y);
    errs++;
  }
  if (y > gd->data->height) {
    printf("ERROR: gdfGetNthSubjectMeasurement: y=%d > gd->data->height=%d\n",
           y, gd->data->height);
    errs++;
  }
  if (z < 0) {
    printf("ERROR: gdfGetNthSubjectMeasurement: z=%d < 0\n", z);
    errs++;
  }
  if (z > gd->data->depth) {
    printf("ERROR: gdfGetNthSubjectMeasurement: z=%d > gd->data->depth=%d\n",
           z, gd->data->depth);
    errs++;
  }
  if (nsubject < 0) {
    printf("ERROR: gdfGetNthSubjectMeasurement: nsubject=%d < 0\n",
           nsubject);
    errs++;
  }
  if (nsubject >= gd->data->nframes) {
    printf("ERROR: gdfGetNthSubjectMeasurement: "
           "nsubject=%d >= gd->data->nframes=%d\n",
           nsubject, gd->data->nframes);
    errs++;
  }
  if (errs) return -1;

  switch ( gd->data->type ) {
  case MRI_UCHAR:
    v = MRIseq_vox(gd->data,x,y,z,nsubject);
    *value = v;
    break;
  case MRI_INT:
    v = MRIIseq_vox(gd->data,x,y,z,nsubject);
    *value = v;
    break;
  case MRI_LONG:
    v = MRILseq_vox(gd->data,x,y,z,nsubject);
    *value = v;
    break;
  case MRI_FLOAT:
    v = MRIFseq_vox(gd->data,x,y,z,nsubject);
    *value = v;
    break;
  case MRI_SHORT:
    v = MRISseq_vox(gd->data,x,y,z,nsubject);
    *value = v;
    break;
  default:
    break ;
  }

  return(0);
}


/*-------------------------------------------------------
  gdfSubSet() - creates a new FSGD with only the Classes
  and Variables listed. If nClasses == -1, all classes
  are included. If nVars == -1, all variables are included.
  -------------------------------------------------------*/
FSGD *gdfSubSet(FSGD *infsgd, int nClasses, char **ClassList,
                int nVars, char **VarList)
{
  FSGD *fsgd;
  int n, nCUse, nVUse, c, ic, v, iv, ninputs, ok;

  if (nClasses > 0) {
    nCUse = nClasses;
    for (n=0; n < nClasses; n++) {
      if (gdfClassNo(infsgd,ClassList[n]) == -1) {
        printf("ERROR: gdfSubSet: class %s not found\n",ClassList[n]);
        return(NULL);
      }
    }
  } else nCUse = infsgd->nclasses;
  if (nVars >= 0) {
    nVUse = nVars;
    for (n=0; n < nVars; n++) {
      if (gdfGetVarLabelNo(infsgd,VarList[n]) == -1) {
        printf("ERROR: gdfSubSet: var %s not found\n",VarList[n]);
        return(NULL);
      }
    }
  } else nVUse = infsgd->nvariables;

  fsgd = gdfAlloc(1);

  for (c = 0; c < nCUse; c++) {
    if (nClasses > 0) ic = gdfClassNo(infsgd,ClassList[c]);
    else             ic = c;
    strcpy(fsgd->classlabel[c], infsgd->classlabel[ic]);
    strcpy(fsgd->classmarker[c],infsgd->classmarker[ic]);
    strcpy(fsgd->classcolor[c], infsgd->classcolor[ic]);
    //printf("c = %d, ic = %d   %s  %s\n",
    //   c,ic,fsgd->classlabel[c], infsgd->classlabel[ic]);
  }

  fsgd->nclasses = nCUse;
  for (v = 0; v < nVUse; v++) {
    if (nVars > 0) iv = gdfGetVarLabelNo(infsgd,VarList[v]);
    else          iv = v;
    strcpy(fsgd->varlabel[v],infsgd->varlabel[iv]);
    //printf("v = %d, iv = %d   %s  %s\n",
    //   v,iv,fsgd->varlabel[v], infsgd->varlabel[iv]);
  }
  fsgd->nvariables = nVUse;

  ninputs = 0;
  for (n = 0; n < infsgd->ninputs; n++) {

    ic = infsgd->subjclassno[n];
    ok = gdfClassNo(fsgd,infsgd->classlabel[ic]);
    if (ok == -1) continue;

    strcpy(fsgd->subjid[ninputs],infsgd->subjid[n]);
    fsgd->subjclassno[ninputs] = gdfClassNo(fsgd,infsgd->classlabel[ic]);

    v = 0;
    for (iv = 0; iv < infsgd->nvariables; iv++) {
      ok = gdfGetVarLabelNo(fsgd,infsgd->varlabel[iv]);
      if (ok == -1) continue;
      fsgd->varvals[ninputs][v] = infsgd->varvals[n][iv];
      v++;
    }
    ninputs ++;
  }
  fsgd->ninputs = ninputs;

  strcpy(fsgd->title,infsgd->title);
  strcpy(fsgd->measname,infsgd->measname);
  strcpy(fsgd->tessellation,infsgd->tessellation);
  strcpy(fsgd->regsubj,infsgd->regsubj);
  strcpy(fsgd->datafile,infsgd->datafile);
  //strcpy(fsgd->defvarlabel,infsgd->defvarlabel); //need check here
  return(fsgd);
}


/*---------------------------------------------------------
  gdfStringIndex() - gets the 0-based index number of the
  string in the list. If the string is not found, returns -1.
  ---------------------------------------------------------*/
int gdfStringIndex(char *str, char **list, int nlist) {
  int index;

  for (index = 0; index < nlist; index++)
    if (strcmp(str,list[index]) == 0) return(index);
  return(-1);
}


/*---------------------------------------------------------
  gdfCopySubjIdppc - Copy subjid to a pointer to a pointer
  to a char.
  ---------------------------------------------------------*/
char **gdfCopySubjIdppc(FSGD *fsgd) {
  char **ppc;
  int n,len;

  ppc = (char **) calloc(sizeof(char *),fsgd->ninputs);
  for (n=0; n < fsgd->ninputs; n++) {
    len = strlen(fsgd->subjid[n]);
    ppc[n] = (char *) calloc(sizeof(char),len+1);
    memmove(ppc[n],fsgd->subjid[n],len);
    //printf("n=%d, %s\n",n,ppc[n]);
  }

  return(ppc);
}


/*!
\fn MATRIX *gdfContrastDOSS(FSGD *fsgd, float *wClass, float *wCovar) 
\brief creates a contrast matrix for the DODS design
  matrix. wClass are the weights for each class; if NULL, then all
  ones are assumed. wCovar are the weights for each covariate
  (including offset); if NULL, then all ones are assumed. Note:
  wCovar must have nvariables+1 items where the first one
  corresponds to the offset.
*/
MATRIX *gdfContrastDODS(FSGD *fsgd, float *wClass, float *wCovar) 
{
  MATRIX *C;
  float w;
  int c, v, n;

  if(strcasecmp(fsgd->DesignMatMethod,"dods") != 0){
    printf("ERROR: gdfContrastDODS() cannot be used with %s\n",fsgd->DesignMatMethod);
    return(NULL);
  }

  /* Contrast matrix (+1 for offsets) */
  C = MatrixAlloc(1, (fsgd->nvariables+1) * fsgd->nclasses, MATRIX_REAL);

  n = 0;
  for (v=0; v < fsgd->nvariables+1; v++) {
    for (c=0; c < fsgd->nclasses; c++) {
      w = 1;
      if (wClass != NULL) w *= wClass[c];
      if (wCovar != NULL) w *= wCovar[v];
      C->rptr[1][n+1] = w;
      n++;
    }
  }

  return(C);
}

/*!
\fn MATRIX *gdfContrastDOSS(FSGD *fsgd, float *wClass, float *wCovar) 
\brief creates a contrast matrix for the DOSS design matrix. wClass
  are the weights for each class intercept; if NULL, then all ones are
  assumed. wCovar are the weights for each covariate (excluding
  offset); if NULL, then all ones are assumed. wCovar has
  nvariables items
*/
MATRIX *gdfContrastDOSS(FSGD *fsgd, float *wClass, float *wCovar) 
{
  MATRIX *C;
  float w;
  int c, v, n;

  if(strcasecmp(fsgd->DesignMatMethod,"doss") != 0){
    printf("ERROR: gdfContrastDOSS() cannot be used with %s\n",fsgd->DesignMatMethod);
    return(NULL);
  }

  /* Contrast matrix*/
  C = MatrixAlloc(1, fsgd->nvariables+fsgd->nclasses, MATRIX_REAL);

  n = 0;
  for (c=0; c < fsgd->nclasses; c++) {
    w = 1;
    if (wClass != NULL) w *= wClass[c];
    C->rptr[1][n+1] = w;
    n++;
  }
  for (v=0; v < fsgd->nvariables+1; v++) {
    w = 1;
    if (wCovar != NULL) w *= wCovar[v];
    C->rptr[1][n+1] = w;
    n++;
  }
  return(C);
}


/*-------------------------------------------------------------------------
  gdfGetSDataFromTable() - gets a cell from a table and returns as a string.
  The data is obtained from column datacol (1-based). The row is determined
  by matching the string found in fieldcol (1-based) against the field. Any
  line in the table that begins with # is skipped. Returns NULL if an error
  or if it could not find a match.
  ------------------------------------------------------------------------*/
char *gdfGetSDataFromTable(char *tablefile, char *field,
                           int fieldcol, int datacol) {
  FILE *fp;
  int ncols;
  char *pc, line[2000], tmpstr[2000];
  char *sfield, *sdata;

  fp = fopen(tablefile,"r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n",tablefile);
    return(NULL);
  }

  while (1) {
    pc = fgets(line,2000,fp);
    if (pc == NULL) {
      printf("ERROR: Could not find a match for %s in table file %s col %d\n",
             field,tablefile,fieldcol);
      fclose(fp);
      return(NULL);
    }
    if(line[0] == '#') {
      if(strcmp(field,"IntraCranialVol") != 0) continue;
      sscanf(line,"%*s %*s %s",tmpstr);
      if(strcmp(tmpstr,"IntraCranialVol,") == 0){
	sscanf(line,"%*s %*s %*s %*s %*s %*s %s",tmpstr);
	sdata = strcpyalloc(tmpstr);
	//printf("%s %s\n",field,sdata);
	return(sdata);
      }
      continue;
    }
    ncols = gdfCountItemsInString(line);
    if (fieldcol > ncols) {
      printf("ERROR: table %s has %d cols, but field col %d reqested\n",
             tablefile,ncols,fieldcol);
      fclose(fp);
      return(NULL);
    }
    if (datacol > ncols) {
      printf("ERROR: table %s has %d cols, but data col %d reqested\n",
             tablefile,ncols,datacol);
      fclose(fp);
      return(NULL);
    }
    sfield = gdfGetNthItemFromString(line, fieldcol-1);
    if (stricmp(sfield,field) == 0) {
      sdata = gdfGetNthItemFromString(line, datacol-1);
      free(sfield);
      fclose(fp);
      return(sdata);
    }
    free(sfield);
  }
  // should never get here
  fclose(fp);
  printf("Could not find a match for %s in table file %s col %d\n",
         field,tablefile,fieldcol);
  return(NULL);
}


/*-------------------------------------------------------------------------
  gdfGetDDataFromTable() - same as gdfGetSDataFromTable() but converts
  string to double. If an error is encountered, returns 1. Otherwise
  returns 0.
  ------------------------------------------------------------------------*/
int gdfGetDDataFromTable(char *tablefile, char *field,
                         int fieldcol, int datacol, double *data) {
  char *s;
  s = gdfGetSDataFromTable(tablefile, field, fieldcol, datacol);
  if (s == NULL) return(1);
  sscanf(s,"%lf",data);
  free(s);
  return(0);
}
