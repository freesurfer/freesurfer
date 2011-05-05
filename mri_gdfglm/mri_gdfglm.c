/**
 * @file  mri_gdfglm.c
 * @brief performs glm analysis given group descriptor file and dep. var. table
 *
 * Things to do:
 * Class-based/Covar-based partial model fit
 * Spec search ranges for given covariates.
 * Allow spec of subjects to exclude
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2011/05/05 15:28:03 $
 *    $Revision: 1.9 $
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
#include <unistd.h>
#include <string.h>

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "matfile.h"

#include "matrix.h"
#include "fsgdf.h"
#include "fio.h"
#include "sig.h"
#include "version.h"

#ifdef X
#undef X
#endif

typedef struct tagDEPVARTABLE {
  char *dvtfile;
  char **RowNames;
  char **ColNames;
  MATRIX *D;
  int transposed;
}
DEPVARTABLE, DVT;

DVT *DVTalloc(int nrows, int ncols, char *dvtfile);
DVT *DVTread(char *dvtfile);
int  DVTdump(DVT *dvt, FILE *fp);
int  DVTsetName(DVT *dvt, int nth, char *Name, int dimid);
int  DVTfree(DVT **ppdvt);
DVT *DVTtranspose(DVT *dvt);
int DVTgetNameIndex(DVT *dvt, char *Name, int dimid);
MATRIX *DVTgetDepVar(DVT *dvt, int nDepVars, char **DepVarList,
                     float *DepVarWeight);
DVT *DVTpruneRows(DVT *indvt, int nKeep, char **KeepList);

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
static int CheckReps(char **List, int nList);
static double ContrastVMF(MATRIX *X, MATRIX *C);
static int WriteAllClassDat(char *base, FSGD *fsgd,
                            MATRIX *y, MATRIX *yhat,
                            MATRIX *X, MATRIX *beta);

static int WriteClassDat(char *base, char *Class, FSGD *fsgd,
                         MATRIX *y, MATRIX *yhat, MATRIX *X,
                         MATRIX *beta);

//static int  stringmatch(char *str1, char *str2);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_gdfglm.c,v 1.9 2011/05/05 15:28:03 greve Exp $";
const char *Progname = "mri_gdfglm";

typedef struct tagCOVARPRUNE {
  char *Covar;
  float Min, Max;
}
COVARPRUNE;

int debug = 0;

char *GDFile;
char *DVTFile;

int  nClassList = 0;
char *ClassList[100];

int   nwClass = 0;
float *wClass;

int  nCovarList = 0;
char *CovarList[100];

int  nwCovar = 0;
float *wCovar = NULL;

int  nCovarPrune = 0;
COVARPRUNE CovarPrune[200];
int TestOffset = 0;

int   nDepVarList = 0;
char  *DepVarList[100];

int   nwDepVar = 0;
float *wDepVar = NULL;

char *WLMSDepVar = NULL;

char *OutBase = NULL;
int KeepSubjId = 0;

DVT *dvt0, *dvt;
FSGD *fsgd0, *fsgd;

MATRIX *X, *Xnorm, *pinvX;
MATRIX *y,  *beta, *yhat, *r, *ces, *C;
MATRIX *all;
float Xcondition;
double rvar, rmean, dof, tval, sigtval, vmf;

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  int n, v, c;
  FILE *fp;
  char *covarname;
  char SumFile[2000];
  char DatFile[2000];
  char MatFile[2000];
  char OutGDFile[2000];
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_gdfglm.c,v 1.9 2011/05/05 15:28:03 greve Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  printf("\n\n");
  printf("%s ",Progname);
  for (n=0; n < argc; n++) printf("%s ",argv[n]);
  printf("\n\n");
  printf("%s\n\n",vcid);

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  X = gdfMatrixDODS(fsgd,NULL);
  if (X==NULL) exit(1);
  if (debug) MatrixPrint(stdout,X);


  Xnorm = MatrixNormalizeCol(X,NULL,NULL);
  Xcondition = sqrt(MatrixNSConditionNumber(Xnorm));
  MatrixFree(&Xnorm);
  printf("INFO: Normalized Design Matrix Condition Number is %g\n",
         Xcondition);
  if (Xcondition > 100000) {
    printf("ERROR: Design matrix is badly conditioned, check for linear\n"
           "dependency  between columns (ie, two or more columns \n"
           "that add up to another column).\n\n");
    exit(1);
  }

  printf("Extracting DepVar\n");
  y = DVTgetDepVar(dvt,nDepVarList,DepVarList,wDepVar);

  printf("Performing Estimation\n");
  pinvX = MatrixPseudoInverse(X,NULL);
  beta = MatrixMultiply(pinvX,y,NULL);
  yhat = MatrixMultiply(X,beta,NULL);
  r = MatrixSubtract(y,yhat,NULL);
  dof = X->rows-X->cols;
  rvar = VectorVar(r, &rmean);
  rvar = rvar * (X->rows-1)/dof;

  printf("Beta: -----------------\n");
  MatrixPrint(stdout,beta);
  printf("---------------------------------\n\n");
  printf("rvar = %g, rstd = %g\n",rvar,sqrt(rvar));

  C = gdfContrastDODS(fsgd, wClass, wCovar);
  printf("C: -----------------\n");
  MatrixPrint(stdout,C);
  printf("---------------------------------\n\n");
  ces = MatrixMultiply(C,beta,NULL);
  vmf = ContrastVMF(X,C);
  tval = ces->rptr[1][1]/sqrt(rvar*vmf);
  sigtval = sigt(tval, rint(dof));
  printf("ces = %g, vmf = %g, t = %g, sigt = %g\n",
         ces->rptr[1][1],vmf,tval,sigtval);

  sprintf(SumFile,"%s.sum",OutBase);
  fp = fopen(SumFile,"w");
  fprintf(fp,"mri_gdfglm summary file\n\n");
  fprintf(fp,"Group Descriptor File %s\n",GDFile);
  fprintf(fp,"Dependent Variable File %s\n",DVTFile);
  fprintf(fp,"Dependent Variable Weights: ");
  if (wDepVar == NULL)
    fprintf(fp," all 1s\n");
  else {
    fprintf(fp,"\n");
    for (n=0; n < nwDepVar; n++)
      fprintf(fp," %s %g\n",DepVarList[n],wDepVar[n]);
  }

  fprintf(fp,"\n");
  fprintf(fp,"Class Contrast Weights: ");
  if (nwClass == 0)
    fprintf(fp," all 1s\n");
  else {
    fprintf(fp,"\n");
    for (n=0; n < nwClass; n++)
      fprintf(fp," %s %g\n",fsgd->classlabel[n],wClass[n]);
  }
  fprintf(fp,"\n");

  fprintf(fp,"Covar Contrast Weights: ");
  if (nwCovar == 0)
    if (!TestOffset) fprintf(fp," all 1s\n");
    else            fprintf(fp," all 0s\n");
  else {
    fprintf(fp,"\n");
    for (n=0; n < nwCovar; n++)
      fprintf(fp," %s %g",CovarList[n],wCovar[n]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"TestOffset = %d\n",TestOffset);
  fprintf(fp,"\n");

  fprintf(fp,"Parameter Estimates and Contrast Weighting:\n\n");
  n = 0;
  for (v=0; v < fsgd->nvariables+1; v++) {
    if (v==0) covarname = "Offset";
    else     covarname = fsgd->varlabel[v-1];
    for (c=0; c < fsgd->nclasses; c++) {
      fprintf(fp,"%-10s %-10s  %12.5f   %5.2f\n",fsgd->classlabel[c],
              covarname,beta->rptr[n+1][1],C->rptr[1][n+1]);
      n++;
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");

  fprintf(fp,"Residual Variance %g\n",rvar);
  fprintf(fp,"Residual StdDev   %g\n",sqrt(rvar));
  fprintf(fp,"DOF %g\n",dof);
  fprintf(fp,"\n");

  fprintf(fp,"Contrast Effect Size       %g\n",ces->rptr[1][1]);
  fprintf(fp,"Variance Reduction Factor  %g\n",1/vmf);
  fprintf(fp,"t-Ratio                    %g\n",tval);
  fprintf(fp,"Significance               %g\n",sigtval);
  fprintf(fp,"\n");
  fclose(fp);

  /*----------------------------------------*/
  sprintf(DatFile,"%s.dat",OutBase);
  fp = fopen(DatFile,"w");
  for (n=0; n < fsgd->ninputs; n++) {
    fprintf(fp,"%2d ",n);
    if (KeepSubjId) fprintf(fp,"%s",fsgd->subjid[n]);
    for (v=0; v < fsgd->nvariables; v++)
      fprintf(fp," %g",fsgd->varvals[n][v]);
    fprintf(fp," %g %g",y->rptr[n+1][1],yhat->rptr[n+1][1]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  /*----------------------------------------*/
  sprintf(MatFile,"%s.mat",OutBase);
  all = MatrixHorCat(X,y,NULL);
  all = MatrixHorCat(all,yhat,NULL);
  all = MatrixHorCat(all,r,NULL);
  MatlabWrite(all,MatFile,"X");

  /*----------------------------------------*/
  sprintf(OutGDFile,"%s.gdf",OutBase);
  fp = fopen(OutGDFile,"w");
  gdfPrintHeader(fp,fsgd);
  fclose(fp);

  /*----------------------------------------*/
  WriteAllClassDat(OutBase,fsgd,y,yhat,X,beta);


  return(0);
}
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  DVT *dvttmp;
  float ftmp[2000];

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--keepid"))  KeepSubjId = 1;

    else if (!strcmp(option, "--gdf")) {
      if (nargc < 1) argnerr(option,1);
      GDFile = pargv[0];
      fsgd0 = gdfRead(GDFile,0);
      if (fsgd0 == NULL) exit(1);
      nargsused = 1;
    } else if (!strcmp(option, "--dvt")) {
      if (nargc < 1) argnerr(option,1);
      DVTFile = pargv[0];
      dvttmp = DVTread(DVTFile);
      dvt0 = DVTtranspose(dvttmp);
      DVTfree(&dvttmp);
      nargsused = 1;
    } else if (!strcmp(option, "--classes")) {
      if (nargc < 1) argnerr(option,1);
      while (nth_is_arg(nargc, pargv,  nargsused)) {
        ClassList[nClassList] = pargv[nargsused];
        nClassList++;
        nargsused++;
      }
    } else if (!strcmp(option, "--covar")) {
      if (nargc < 1) argnerr(option,1);
      while (nth_is_arg(nargc, pargv, nargsused)) {
        CovarList[nCovarList] = pargv[nargsused];
        nCovarList++;
        nargsused++;
      }
    } else if (!strcmp(option, "--covarprune")) {
      if (nargc < 3) argnerr(option,3);
      CovarPrune[nCovarPrune].Covar = pargv[0];
      sscanf(pargv[1],"%f",&CovarPrune[nCovarPrune].Min);
      sscanf(pargv[2],"%f",&CovarPrune[nCovarPrune].Max);
      if (CovarPrune[nCovarPrune].Min > CovarPrune[nCovarPrune].Max) {
        printf("ERROR: covarprune min > max\n");
        exit(1);
      }
      nCovarPrune++;
      nargsused = 3;
    } else if (!strcmp(option, "--depvar")) {
      if (nargc < 1) argnerr(option,1);
      while (nth_is_arg(nargc, pargv, nargsused)) {
        DepVarList[nDepVarList] = pargv[nargsused];
        nDepVarList++;
        nargsused++;
      }
    } else if (!strcmp(option, "--wdepvar")) {
      if (nargc < 1) argnerr(option,1);
      while (nth_is_arg(nargc, pargv, nargsused)) {
        sscanf(pargv[nargsused],"%f",&ftmp[nwDepVar]);
        nwDepVar++;
        nargsused++;
      }
      wDepVar = (float *)calloc(sizeof(float),nwDepVar);
      memmove(wDepVar,ftmp,sizeof(float)*nwDepVar);
    } else if (!strcmp(option, "--wlms")) {
      if (nargc < 1) argnerr(option,1);
      WLMSDepVar = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--wclass")) {
      if (nargc < 1) argnerr(option,1);
      while (nth_is_arg(nargc, pargv, nargsused)) {
        sscanf(pargv[nargsused],"%f",&ftmp[nwClass]);
        nwClass++;
        nargsused++;
      }
      wClass = (float *)calloc(sizeof(float),nwClass);
      memmove(wClass,ftmp,sizeof(float)*nwClass);
    } else if (!strcmp(option, "--wcovar")) {
      if (nargc < 1) argnerr(option,1);
      while (nth_is_arg(nargc, pargv, nargsused)) {
        sscanf(pargv[nargsused],"%f",&ftmp[nwCovar]);
        nwCovar++;
        nargsused++;
      }
      wCovar = (float *)calloc(sizeof(float),nwCovar+1);
      memmove(&wCovar[1],ftmp,sizeof(float)*nwCovar);
      wCovar[0] = 0.0;
      //Note: the first weight is for offset
    } else if (!strcmp(option, "--testoffset")) TestOffset = 1;

    else if (!strcmp(option, "--o")) {
      if (nargc < 1) argnerr(option,1);
      OutBase = pargv[0];
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --gdf gdfile : group descriptor file\n");
  printf("   --dvt dvtfile : dependent variable file  \n");

  printf("\n");
  printf("   --classes Class1 Class2 : use subset of classes  \n");
  printf("   --covar Covar1 Corvar2 ... : use subset of covars\n");
  //printf("   --covarprune Covar Min Max : exclude when out of range\n");
  printf("   --depvar DepVar1 <DepVar2 ...> : spec dependent variables \n");
  printf("   --wdepvar wdv1 wdv2 ... : weight depvars (default is 1) \n");
  //printf("   --wlms DepVar \n");


  printf("\n");
  printf("   --wclass wc1 wc2 ... : Class weights (def 1)\n");
  printf("   --wcovar wcv1 wcvw ... : Covar slope weights  \n");
  printf("   --testoffset : test offset, not covariate slope \n");

  printf("\n");
  printf("   --o output base name\n");
  printf("   --keepid : print subjid in output.dat\n");

  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf(
    "Performs glm analysis given group descriptor file (GDF) and dependent\n"
    "variable table (DVT).\n"
    "\n"
    "\n"
    "--gdf gdffile\n"
    "\n"
    "Path to the GDF. See http://surfer.nmr.mgh.harvard.edu/docs/fsgdf.txt\n"
    "for more info. This file will have a list of Classes and Variables.\n"
    "Hereafter, the Variables are referred to as Covariates.\n"
    "\n"
    "--dvt dvtfile\n"
    "\n"
    "Path to the dependent variable table (DVT). This is a text file that\n"
    "contains the data table. The first column is the list the names of\n"
    "each row. The first row is the list of the names of each column. There\n"
    "needs to be a text place-holder at the first row/column. The rest\n"
    "of the table is filled with numbers. Each column should be a subject, \n"
    "and each row an observation for that subject. The arguments of --depvar\n"
    "must come from the row names. A DVT is produced by make-segvol-table.\n"
    "\n"
    "--classes Class1 <Class2 ...>\n"
    "\n"
    "Use only the subjects that belong to the specfied list of classes.\n"
    "The class names must come from the classes as specified in the GDF. If\n"
    "unspecfied, all classes are used.\n"
    "\n"
    "--covar Covar1 <Covar2 ...>\n"
    "\n"
    "Use only the variables that belong to the specfied list. The names\n"
    "must come from the variables as specified in the GDF. If unspecfied,\n"
    "all variables are used.\n"
    "\n"
    "--depvar DepVar1 <DepVar2 ...>\n"
    "\n"
    "Select variables from the DVT. The a weighted average of variables\n"
    "will be computed as the final dependent variable. If unspecified, all\n"
    "variables will be used. See --wdepvar.\n"
    "\n"
    "--wdepvar wDepVar1 wDepVar2 ...\n"
    "\n"
    "Set the weights of the dependent variables. The final dependent\n"
    "variable will be computed as a weighted average of the dependent\n"
    "variables. The number of weights must be equal to either the number of\n"
    "DepVars listed in --depvar or (if --depvar is unspecfied) the number\n"
    "of dependent variables in the DVT. If unspecfied, the weights are set\n"
    "to compute a simple average.\n"
    "\n"
    "--wclass WC1 WC2 ...\n"
    "\n"
    "Class weights for establishing a contrast. The number of weights must\n"
    "be equal to the number of classes (ie, the number listed in --classes\n"
    "or the number in the GDF). If unspecified, all weights are set to 1.\n"
    "This applies only to the contrast; if the weight of a class is set to\n"
    "0, that class is still included in the parameter estimation. If\n"
    "positive and negative weights are used, they should sum to the same\n"
    "value.\n"
    "\n"
    "--wcovar WCV1 WCV2 ...\n"
    "\n"
    "Covariate weights for establishing a contrast. The number of weights\n"
    "must be equal to the number of covariates (ie, the number listed in\n"
    "--covar or the number in the GDF). If unspecified, all weights are set\n"
    "to 1.  This applies only to the contrast; if the weight of a covariate\n"
    "is set to 0, that covariate is still included in the parameter\n"
    "estimation. If positive and negative weights are used, they should sum\n"
    "to the same value.\n"
    "\n"
    "--testoffset\n"
    "\n"
    "The offset is like a special covariate.\n"
    "\n"
    "--o basename\n"
    "\n"
    "Base name of output files. There will be four output files created:\n"
    "(1) the summary file (basename.sum), (2) the data file (basename.dat),\n"
    "(3) the GDF (basename.gdf), and a matrix file (basename.mat). The \n"
    "summary file has a list of the parameters used in the analysis\n"
    "as well as the results, including parameter estimates, contrast\n"
    "effect size, and signficance. The data file contains a table of \n"
    "the final data with each subject on a different row. The first\n"
    "column is the subject number, the next nCV are the nCV covariates,\n"
    "the next column is the final dependent variable, and the final column \n"
    "is the best fit of the dependent variable. The GDF is the final\n"
    "GDF; this will be the same as the input GDF if no classes, covariates\n"
    "or subjects have been excluded. The matfile is a matrix in matlab4\n"
    "format that contains the design matrix concatenated with the \n"
    "final dependent variable, the fit of final dependent variable,\n"
    "and the residual.\n"
    "\n"
    "In addition, each class has its own output dat file called \n"
    "outbase-classlabel.dat. The first column is the subject number,\n"
    "the next nCV are the nCV covariates, the next column is the final \n"
    "dependent variable, and the final column is the best fit of the \n"
    "dependent variable. This output is best for creating scatter\n"
    "plots.\n"
    "\n"
    "EXAMPLES:\n"
    "\n"
    "Consider the following Group Descriptor File:\n"
    "------------------------------------------\n"
    "GroupDescriptorFile 1\n"
    "  Title AlbertGroup\n"
    "  Class NormMale   \n"
    "  Class AlzMale    \n"
    "  Class NormFemale \n"
    "  CLASS AlzFemale  \n"
    "  DefaultVariable Age\n"
    "  Variables   Age MMSE\n"
    "Input 003007 AlzMale 75 30 \n"
    "...\n"
    "------------------------------------------\n"
    "\n"
    "(1) Test whether the left hippocampal volume signficantly varies \n"
    "with MMSE for the Alzheimer group regressing out the effect of gender \n"
    "and age:\n"
    "\n"
    "mri_gdfglm --gdf albert.gdf --dvt asegvol.dat --o lhipmmse \n"
    "  --wclass 0.5 -0.5 0.5 -0.5  --wcovar 0 1 \n"
    "  --depvar Left-Hippocampus \n"
    "\n"
    "(2) Test whether the left-right difference in hippocampal volume \n"
    "signficantly varies across gender using only the normal subjects,\n"
    "and regressing out the effect of age but not MMSE.\n"
    "\n"
    "mri_gdfglm --gdf albert.gdf --dvt asegvol.dat --o hip-lrdiff-mvf \n"
    "  --depvar Left-Hippocampus Right-Hippocampus --wdepvar  1 -1 \n"
    "  --classes NormMale NormFemale   --wclass  1 -1  --covar Age \n"
    "  --testoffset \n"
  );
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void) {
  int n,nC,nV;
  char **ppctmp;

  if (GDFile == NULL) {
    printf("ERROR: no gdf file\n");
    exit(1);
  }
  if (DVTFile == NULL) {
    printf("ERROR: no dep var table file\n");
    exit(1);
  }
  if (OutBase == NULL) {
    printf("ERROR: no output specified\n");
    exit(1);
  }
  if (nDepVarList == 0) {
    printf("ERROR: no dep vars specified\n");
    exit(1);
  } else {
    for (n=0; n < nDepVarList; n++) {
      if (DVTgetNameIndex(dvt0,DepVarList[n],2) == -1) {
        printf("ERROR: DepVar %s does not exist\n",DepVarList[n]);
        exit(1);
      }
    }
  }
  if (nwDepVar != 0 && nDepVarList != 0 && nwDepVar != nDepVarList ) {
    printf("ERROR: diff number of dep vars and weights\n");
    exit(1);
  }
  if (nClassList != 0) {
    for (n=0; n < nClassList; n++) {
      if (gdfClassNo(fsgd0,ClassList[n]) == -1) {
        printf("ERROR: Class %s does not exist\n",ClassList[n]);
        exit(1);
      }
    }
    if (CheckReps(ClassList,nClassList)) {
      printf("ERROR: repetition found in class list\n");
      exit(1);
    }
  }
  if (nCovarList != 0) {
    for (n=0; n < nCovarList; n++) {
      if (gdfGetVarLabelNo(fsgd0,CovarList[n]) == -1) {
        printf("ERROR: Covar %s does not exist\n",CovarList[n]);
        exit(1);
      }
    }
    if (CheckReps(CovarList,nCovarList)) {
      printf("ERROR: repetition found in covar list\n");
      exit(1);
    }
  }

  if (CheckReps(DepVarList,nDepVarList)) {
    printf("ERROR: repetition found in dep var list\n");
    exit(1);
  }

  if (WLMSDepVar != NULL) {
    if (DVTgetNameIndex(dvt0,WLMSDepVar,2) == -1) {
      printf("ERROR: WLMS DepVar %s does not exist\n",WLMSDepVar);
      exit(1);
    }
    /* Should have a check here to make sure that wlms depvar
       is not redundant with DepVarList and that DepVarList
       exists. */
  }

  /* Get new fsgd if needed */
  if (nClassList != 0 || nCovarList != 0) {
    if (nClassList == 0) nC = -1;
    else                nC = nClassList;
    if (nCovarList == 0) nV = -1;
    else                nV = nCovarList;
    printf("Getting FSGDF Subset nC=%d, nV=%d\n",nC,nV);
    fflush(stdout);
    fsgd = gdfSubSet(fsgd0,nC,ClassList,nV,CovarList);
    if (fsgd == NULL) exit(1);
  } else fsgd = fsgd0;

  /* Prune DVT if needed */
  if (nClassList != 0) {
    printf("INFO: pruning DVT rows\n");
    ppctmp = gdfCopySubjIdppc(fsgd);
    //DVTdump(dvt0,stdout);
    dvt = DVTpruneRows(dvt0,fsgd->ninputs,ppctmp);
    if (dvt==NULL) exit(1);
  } else dvt = dvt0;

  if (nwClass != 0 && nClassList != 0 && nwClass != nClassList ) {
    printf("ERROR: diff number of classes and class weights\n");
    exit(1);
  }
  if (nwClass != 0 && nwClass != fsgd->nclasses) {
    printf("ERROR: diff number of classes and class weights\n");
    exit(1);
  }
  if (TestOffset) {
    if (nwCovar == 0)
      wCovar = (float *) calloc(sizeof(float),fsgd->nvariables+1);
    wCovar[0] = 1.0;
  } else if (wCovar == NULL) {
    wCovar = (float *) calloc(sizeof(float),fsgd->nvariables+1);
    for (n=1; n < fsgd->nvariables+1; n++) wCovar[n] = 1.0;
    wCovar[0] = 0.0;
  }
  if (debug) {
    for (n=0; n < fsgd->nvariables+1; n++)
      printf("wCovar[%d] = %g\n",n,wCovar[n]);
  }

  if (nwCovar != 0 && nCovarList != 0 && nwCovar != nCovarList ) {
    printf("ERROR: diff number of covars and covar weights\n");
    exit(1);
  }
  if (nwCovar != 0 && nwCovar != fsgd->nvariables ) {
    printf("ERROR: diff number of covars and covar weights\n");
    exit(1);
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  int n;
  fprintf(fp,"GDFile %s \n",GDFile);
  fprintf(fp,"DVTFile %s \n", DVTFile);

  if (nClassList > 0) {
    fprintf(fp,"%d subset of classes selected \n",nClassList);
    for (n = 0; n < nClassList; n++) {
      fprintf(fp,"%2d %s \n",n+1,ClassList[n]);
    }
  } else
    fprintf(fp,"Using All Classes\n");

  if (nwClass > 0) {
    fprintf(fp,"Class Weights\n");
    for (n = 0; n < nwClass; n++) {
      fprintf(fp,"%2d %g \n",n+1,wClass[n]);
    }
  } else
    fprintf(fp,"All Classes equally weighted.\n");

  if (nCovarList > 0) {
    fprintf(fp,"%d subset of covars selected \n",nCovarList);
    for (n = 0; n < nCovarList; n++) {
      fprintf(fp,"%2d %s \n",n+1,CovarList[n]);
    }
  } else
    fprintf(fp,"Using All Covars\n");

  if (nwCovar > 0) {
    fprintf(fp,"Covar Weights\n");
    for (n = 0; n < nwCovar; n++) {
      fprintf(fp,"%2d %g \n",n+1,wCovar[n]);
    }
  } else
    fprintf(fp,"All Covars equally weighted.\n");

  if (nCovarPrune) {
    fprintf(fp,"Pruning used on %d Covars\n",nCovarPrune);
    for (n = 0; n < nCovarPrune; n++) {
      fprintf(fp,"  %2d  %s  %g  %g\n",n+1,CovarPrune[n].Covar,
              CovarPrune[n].Min,CovarPrune[n].Max);
    }
  }

  fprintf(fp,"Dependent Variables (%d)\n",nDepVarList);
  for (n = 0; n < nDepVarList; n++) {
    fprintf(fp,"%2d %s ",n+1,DepVarList[n]);
    if (nwDepVar > 0)
      fprintf(fp," %g \n",wDepVar[n]);
    else
      fprintf(fp," 1\n");
  }

  if (TestOffset)
    fprintf(fp,"Testing Offsets \n");

  if (WLMSDepVar != NULL)
    fprintf(fp,"WLMSDepVar: %s \n",WLMSDepVar);

  fprintf(fp,"OutBase: %s \n",OutBase);
  fprintf(fp," \n");
  fprintf(fp," \n");

  if (debug) {
    fprintf(fp,"----------------------------------------\n");
    fprintf(fp,"Dependent Variable Table\n");
    DVTdump(dvt,stdout);
    fprintf(fp,"----------------------------------------\n");
    if (dvt != dvt0) {
      fprintf(fp,"----------------------------------------\n");
      fprintf(fp,"Dependent Variable Table (Original)\n");
      DVTdump(dvt0,stdout);
      fprintf(fp,"----------------------------------------\n");
    }

    fprintf(fp,"Group Descriptor\n");
    gdfPrintHeader(fp,fsgd);
    fprintf(fp,"----------------------------------------\n");
    if (fsgd != fsgd0) {
      fprintf(fp,"Group Descriptor (Original)\n");
      gdfPrintHeader(fp,fsgd0);
      fprintf(fp,"----------------------------------------\n");
    }
  }
  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int isflag(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth) {
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
  if (nargc <= nth) return(0);

  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) return(0);

  return(1);
}
/*------------------------------------------------------------*/
//static int stringmatch(char *str1, char *str2)
//{
//  if(! strcmp(str1,str2)) return(1);
//  return(0);
//}

/*------------------------------------------------------------*/
static int CheckReps(char **List, int nList) {
  int n1, n2;

  for (n1 = 0; n1 < nList; n1++) {
    for (n2 = n1+1; n2 < nList; n2++) {
      if (strcmp(List[n1],List[n2])==0) return(1);
    }
  }
  return(0);
}

/*------------------------------------------------------------
  ContrastVMF() - computes the variance multiplication factor
  ------------------------------------------------------------*/
static double ContrastVMF(MATRIX *X, MATRIX *C) {
  float vmf;
  MATRIX *Xt, *XtX, *iXtX, *CiXtX, *Ct, *CiXtXCt;

  Xt = MatrixTranspose(X,NULL);
  XtX = MatrixMultiply(Xt,X,NULL);
  iXtX = MatrixInverse(XtX,NULL);
  CiXtX = MatrixMultiply(C,iXtX,NULL);
  Ct = MatrixTranspose(C,NULL);
  CiXtXCt = MatrixMultiply(CiXtX,Ct,NULL);

  vmf = CiXtXCt->rptr[1][1];

  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&CiXtX);
  MatrixFree(&Ct);
  MatrixFree(&CiXtXCt);

  return(vmf);
}
/*------------------------------------------------------------*/
static int WriteAllClassDat(char *base, FSGD *fsgd,
                            MATRIX *y, MATRIX *yhat, MATRIX *X, MATRIX *beta) {
  int nc, err;

  for (nc = 0; nc < fsgd->nclasses; nc++) {
    err = WriteClassDat(base,fsgd->classlabel[nc],fsgd,
                        y,yhat,X,beta);
    if (err) return(1);
  }

  return(0);
}
/*------------------------------------------------------------*/
static int WriteClassDat(char *base, char *Class, FSGD *fsgd,
                         MATRIX *y, MATRIX *yhat, MATRIX *X, MATRIX *beta) {
  FILE *fp;
  int ClassNo;
  char fname[2000];
  int n, nth, v;

  ClassNo = gdfClassNo(fsgd,Class);
  if (ClassNo == -1) {
    printf("ERROR: Class %s not found in fsgd\n",Class);
    return(1);
  }

  /* SubjNo Covar1 Covar2 ... DepVar DepVarHat */
  sprintf(fname,"%s-%s.dat",base,Class);
  fp = fopen(fname,"w");
  nth = 0;
  for (n=0; n < fsgd->ninputs; n++) {
    if (fsgd->subjclassno[n] != ClassNo) continue;
    fprintf(fp,"%2d",n); /* use absolute subj no */
    for (v=0; v < fsgd->nvariables; v++)
      fprintf(fp," %g",fsgd->varvals[n][v]);
    fprintf(fp," %g",y->rptr[n+1][1]);
    fprintf(fp," %g",yhat->rptr[n+1][1]);
    fprintf(fp,"\n");
    nth++;
  }
  fclose(fp);

  return(0);
}




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
DVT *DVTread(char *dvtfile) {
  DVT *dvt = NULL;
  FILE *fp;
  int nrows, ncols, r, c;
  char tmpstr[4001];

  if (!fio_FileExistsReadable(dvtfile)) {
    printf("ERROR: %s does not exist or is not readable\n",dvtfile);
    return(NULL);
  }

  nrows = fio_NLines(dvtfile);
  if (nrows == 0) {
    printf("ERROR: reading %s, no rows found\n",dvtfile);
    return(NULL);
  }

  fp = fopen(dvtfile,"r");
  fgets(tmpstr,4000,fp);
  ncols = gdfCountItemsInString(tmpstr);
  if (ncols == 0) {
    printf("ERROR: reading %s, no cols found\n",dvtfile);
    fclose(fp);
    return(NULL);
  }
  fclose(fp);

  nrows = nrows-1;
  ncols = ncols-1;
  printf("%s  nrows = %d, ncols = %d\n",dvtfile,nrows,ncols);

  dvt = DVTalloc(nrows,ncols,dvtfile);
  if (dvt == NULL) return(NULL);
  dvt->transposed = 0;

  fp = fopen(dvtfile,"r");

  /* Read the first row = names of columns */
  fscanf(fp,"%s",tmpstr); /* swallow the first one */
  for (c=0; c < ncols; c++) {
    fscanf(fp,"%s",tmpstr);
    DVTsetName(dvt,c,tmpstr,2);
  }

  for (r=0; r < nrows; r++) {
    /* First column is row name */
    fscanf(fp,"%s",tmpstr);
    DVTsetName(dvt,r,tmpstr,1);
    //printf("r=%d, RowName = %s\n",r,dvt->RowNames[r]);
    for (c=0; c < ncols; c++) {
      fscanf(fp,"%f",&dvt->D->rptr[r+1][c+1]);
    }
  }

  fclose(fp);
  return(dvt);
}
/*-----------------------------------------------------*/
int DVTdump(DVT *dvt, FILE *fp) {
  int c,r;

  if (fp == NULL) fp = stdout;

  fprintf(fp,"%s \n",dvt->dvtfile);
  fprintf(fp,"nrows = %d\n",dvt->D->rows);
  fprintf(fp,"ncols = %d\n",dvt->D->cols);
  fprintf(fp,"transposed = %d\n",dvt->transposed);

  fprintf(fp,"\n");
  fprintf(fp,"Column Names \n");
  for (c = 0; c < dvt->D->cols; c++)
    fprintf(fp,"%3d %s\n",c+1,dvt->ColNames[c]);

  fprintf(fp,"\n");
  fprintf(fp,"Row Names \n");
  for (r = 0; r < dvt->D->rows; r++)
    fprintf(fp,"%3d %s\n",r+1,dvt->RowNames[r]);

  fprintf(fp,"\n");
  for (r = 0; r < dvt->D->rows; r++) {
    fprintf(fp,"%3d %s",r+1,dvt->RowNames[r]);
    for (c = 0; c < dvt->D->cols; c++) {
      fprintf(fp," %g",dvt->D->rptr[r+1][c+1]);
    }
    fprintf(fp,"\n");
  }

  fprintf(fp,"\n");

  return(0);

}
/*-----------------------------------------------------*/
DVT *DVTalloc(int nrows, int ncols, char *dvtfile) {
  DVT *dvt;

  dvt = (DVT *)calloc(sizeof(DVT),1);

  if (dvtfile != NULL) {
    dvt->dvtfile = (char *) calloc(sizeof(char),strlen(dvtfile)+1);
    memmove(dvt->dvtfile,dvtfile,strlen(dvtfile));
  }

  dvt->D = MatrixAlloc(nrows,ncols,MATRIX_REAL);
  dvt->RowNames = (char **)calloc(sizeof(char *),nrows);
  dvt->ColNames = (char **)calloc(sizeof(char *),ncols);

  return(dvt);
}
/*-----------------------------------------------------
  DVTsetName() - sets the name of either the nth row or
  nth column depending upon dimid (1=row, 2=col). nth
  is 0-based.
  -----------------------------------------------------*/
int DVTsetName(DVT *dvt, int nth, char *Name, int dimid) {
  int len;

  if ( (dimid == 1 && nth >= dvt->D->rows) ||
       (dimid == 2 && nth >= dvt->D->cols) ) {
    printf("ERROR: nth=%d exceeds matrix dimensions\n",nth);
    return(1);
  }

  len = strlen(Name);
  if (dimid == 1) {
    dvt->RowNames[nth] = (char *) calloc(sizeof(char),len+1);
    memmove(dvt->RowNames[nth],Name,len);
  } else {
    dvt->ColNames[nth] = (char *) calloc(sizeof(char),len+1);
    memmove(dvt->ColNames[nth],Name,len);
  }

  return(0);
}
/*-----------------------------------------------------*/
int DVTfree(DVT **ppdvt) {
  DVT *dvt;
  int r,c;

  dvt = *ppdvt;

  if (dvt->dvtfile != NULL) free(dvt->dvtfile);
  for (c = 0; c < dvt->D->cols; c++) free(dvt->ColNames[c]);
  for (r = 0; r < dvt->D->rows; r++) free(dvt->RowNames[r]);
  free(dvt->ColNames);
  free(dvt->RowNames);
  MatrixFree(&(dvt->D));

  free(*ppdvt);

  return(0);
}
/*-----------------------------------------------------*/
DVT *DVTtranspose(DVT *dvt) {
  DVT *dvtt;
  int r,c;

  dvtt = DVTalloc(dvt->D->cols,dvt->D->rows,dvt->dvtfile);

  for (c = 0; c < dvtt->D->cols; c++)
    DVTsetName(dvtt,c,dvt->RowNames[c],2);
  for (r = 0; r < dvtt->D->rows; r++)
    DVTsetName(dvtt,r,dvt->ColNames[r],1);

  MatrixTranspose(dvt->D, dvtt->D);

  dvtt->transposed = !dvt->transposed;

  return(dvtt);
}
/*-----------------------------------------------------
  DVTgetNameIndex() - gets the 0-based index of Name.
  If dimid=1, gets from RowNames, othwise ColNames.
  Returns -1 if Name is not in list.
  -----------------------------------------------------*/
int DVTgetNameIndex(DVT *dvt, char *Name, int dimid) {
  if (dimid == 1)
    return(gdfStringIndex(Name, dvt->RowNames, dvt->D->rows));
  else
    return(gdfStringIndex(Name, dvt->ColNames, dvt->D->cols));
  return(-1); // should never get here.
}
/*--------------------------------------------------------

  --------------------------------------------------------*/
DVT *DVTpruneRows(DVT *indvt, int nKeep, char **KeepList) {
  DVT *dvt;
  int n, index, c,j;

  dvt = DVTalloc(nKeep,indvt->D->cols,indvt->dvtfile);
  for (n=0; n < nKeep; n++) {
    //printf("n=%d, %s\n",n,KeepList[n]);
    index = gdfStringIndex(KeepList[n],indvt->RowNames,indvt->D->rows);
    if (index == -1) {
      printf("ERROR: DVT RowName --%s-- not found\n",KeepList[n]);
      for (j=0; j < dvt->D->rows; j++)
        printf("%2d %s   %s\n",j,KeepList[n],indvt->RowNames[j]);
      return(NULL);
    }
    DVTsetName(dvt,n,KeepList[n],1);
    for (c=0; c < indvt->D->cols; c++)
      dvt->D->rptr[n+1][c+1] = indvt->D->rptr[index][c+1];
  }

  for (c=0; c < indvt->D->cols; c++)
    DVTsetName(dvt,c,indvt->ColNames[c],2);

  return(dvt);
}

/*------------------------------------------------------
  DVTgetDepVar() - extracts dependent variables, weights,
  returns a single column.
  ------------------------------------------------------*/
MATRIX *DVTgetDepVar(DVT *dvt,
                     int   nDepVars,
                     char **DepVarList,
                     float *DepVarWeight) {
  MATRIX *y, *W;
  int index,n;
  float wc;

  W = MatrixConstVal(1, dvt->D->cols, 1, NULL);
  if (nDepVars != 0) {
    for (n=0; n < dvt->D->cols; n++) {
      index = gdfStringIndex(dvt->ColNames[n],DepVarList,nDepVarList);
      if (index == -1) wc = 0;
      else if (DepVarWeight != NULL) wc = DepVarWeight[index];
      else wc = 1.0;
      W->rptr[n+1][1] = wc;
    }
  }
  //printf("W-------------\n");
  //MatrixPrint(stdout,W);
  //printf("W-------------\n");

  y = MatrixMultiply(dvt->D,W,NULL);
  MatrixFree(&W);

  return(y);
}
