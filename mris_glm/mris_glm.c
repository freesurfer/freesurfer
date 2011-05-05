/**
 * @file  mris_glm.c
 * @brief Computes glm inferences on the surface.
 *
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2011/05/05 15:28:03 $
 *    $Revision: 1.55 $
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

/*
Things to do:
  0. Documentation.
  1. More sophisticated derived variable.
  2. Input volume directly.
  3. Add ability to load in y directly

MC Sim:
  1. Cluster area threshold is in mm^2
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "resample.h"
#include "icosahedron.h"
#include "matrix.h"
#include "matfile.h"
#include "mri.h"
#include "MRIio_old.h"
#include "mri_identify.h"
#include "sig.h"
#include "fmriutils.h"
#include "mri2.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "version.h"
#include "pdf.h"
#include "fsgdf.h"
#include "fio.h"
#include "mri_circulars.h"

#ifdef X
#undef X
#endif

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
static int  stringmatch(char *str1, char *str2);
static int checkfmt(char *fmt);
static int getfmtid(char *fname);
static int IsSurfFmt(char *fmt);
MRIS *MRISloadSurfSubject(char *subj, char *hemi, char *surfid,
                          char *SUBJECTS_DIR);

int ReadAsciiMatrixNRows(char *desmtxfname);
int ReadAsciiMatrixSize(char *desmtxfname, int *pnrows, int *pncols);
int ReadDesignMatrix(char *desmtxfname);
MATRIX *ReadAsciiMatrix(char *asciimtxfname);
int CheckDesignMatrix(MATRIX *X);
static int MatrixWriteFmt(MATRIX *M, char *fname, char *fmt);
static char *getstem(char *bfilename);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mris_glm.c,v 1.55 2011/05/05 15:28:03 greve Exp $";
const char *Progname = "mris_glm";

char *hemi        = NULL;
char *desmtxfname = NULL;
char *fsgdfile = NULL;
char *xmatfile = NULL;
int  xmatonly = 0;
char *xmatfmt = "matlab4";
int  nsmooth   = 0;
int  frame     = 0;

char *surfmeasure    = NULL;

char *surfregid   = "sphere.reg";
int  ninputs = 0;
char *inputlist[1000];
char *inputfmt = NULL;
int   inputfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *subjectlistfile;
int  nsubjects = 0;
char *subjectlist[1000];
int  nregressors = 0;
MATRIX *X; /* design matrix */
FSGD *fsgd=NULL;
char  *gd2mtx_method = "none";

char   *conmtxfname;
MATRIX *C; /* contrast vector */

char *betaid  = NULL;
char *betafmt = NULL;
int  betafmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *beta_in_id  = NULL;
char *beta_in_fmt = NULL;
int  beta_in_fmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *cesid  = NULL;
char *cesfmt = NULL;
int  cesfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *eresid  = NULL;
char *eresfmt = NULL;
int   eresfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *yid  = NULL, *yidbase, *yidstem, *yidbasestem, xmatpath[1000];
char *yfmt = NULL;
int   yfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *yhatid  = NULL;
char *yhatfmt = NULL;
int   yhatfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *eresvarid  = NULL;
char *eresvarfmt = NULL;
int   eresvarfmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *eresvar_in_id  = NULL;
char *eresvar_in_fmt = NULL;
int   eresvar_in_fmtid = MRI_VOLUME_TYPE_UNKNOWN;

char *tid  = NULL;
char *tfmt = NULL;
int   tfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *tmaxfile = NULL;
float tmax;
int nsim = 1, nthsim, MCSim = 0;

char *sigid  = NULL;
char *sigfmt = NULL;
int   sigfmtid = MRI_VOLUME_TYPE_UNKNOWN;

int  IcoOrder     = 7;
float IcoRadius = 100.0;
char *regsurf     = "sphere.reg";
char *trgsubject = NULL;

int dof;
int debug = 0;

int SynthPDF = 0;
int SynthSeed = -1;
double  SynthGaussianMean;
double  SynthGaussianStd;
char   *SynthCDFFile;
double *SynthCDF;
double *SynthXCDF;
int     SynthNCDF;

MRI_SURFACE *IcoSurf=NULL, *SurfReg=NULL;
MRI *SrcVals=NULL, *beta=NULL, *yhat=NULL;
MRI *eres=NULL, *eresvar=NULL, *ces=NULL, *t=NULL, *sig=NULL;
MATRIX *T=NULL, *Xt=NULL, *XtX=NULL, *iXtX=NULL, *Q=NULL, *R=NULL;
MRI *tmpmri=NULL, *tmpmri2=NULL, *SrcHits=NULL;
MRI *SrcDist=NULL, *TrgHits=NULL, *TrgDist=NULL;

float DOF;
char *SUBJECTS_DIR;

int Force=0;
int ParseOnly=0;
char tmpstr[1000];

CHT *cht; // Cluster Hit Table -- for simulations
int nth_ithr, nth_sthr, NClusters;
double ithr, sthr;
SURFCLUSTERSUM *scs;
char *chtfile=NULL;
int n_ithr, n_sthr;
double ithr_lo, ithr_hi, sthr_lo, sthr_hi;
char *ithr_sign;
int abs_flag = 0;

int nvoxels;
FILE *fp;
int DoPermute = 0;

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  int vtx,nthsubj;
  char *subject;
  char *inputfname;
  int  nargs,n;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv,
                                 "$Id: mris_glm.c,v 1.55 2011/05/05 15:28:03 greve Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0];

  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  if (ParseOnly) exit(1);
  dump_options(stdout);

  printf("%s\n",vcid);
  printf("setenv SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  printf("%s\n",getenv("PWD"));
  printf("%s ",Progname);
  for (n=1;n<argc;n++) printf("%s ",argv[n]);
  printf("\n");

  if (xmatfile != NULL) MatrixWriteFmt(X,xmatfile,xmatfmt);

  /* X is the design matrix */
  Xt = MatrixTranspose(X,NULL);
  XtX = MatrixMultiply(Xt,X,NULL);
  iXtX = MatrixInverse(XtX,NULL);
  if (iXtX==NULL) {
    printf("ERROR: could not compute psuedo inverse of X\n");
    exit(1);
  }

  /* Q is the matrix that when multiplied by y gives beta */
  Q = MatrixMultiply(iXtX,Xt,NULL);

  /* T is the matrix that when multiplied by y gives the signal estimate */
  T = MatrixMultiply(X,Q,NULL);

  /* R is the matrix that when multiplied by y gives the residual error */
  R = MatrixSubtract(MatrixIdentity(nsubjects,NULL),T,NULL);
  DOF = X->rows - X->cols;

  printf("Design Matrix ------------------------------------\n");
  MatrixPrint(stdout,X);
  //printf("Q ------------------------------------\n");
  //MatrixPrint(stdout,Q);
  //printf("T ------------------------------------\n");
  //MatrixPrint(stdout,T);
  //printf("R ------------------------------------\n");
  //MatrixPrint(stdout,R);
  printf("Design Covariance Matrix ------------------------------------\n");
  MatrixPrint(stdout,XtX);
  if (C != NULL) {
    printf("Contrast Matrix: -----------------------------\n");
    MatrixPrint(stdout,C);
  }
  printf("DOF = %g\n",DOF);
  if (DOF < 1) {
    printf("ERROR: zero (or fewer) degrees of freedom\n");
    exit(1);
  }
  fflush(stdout);

  /* ------- Load the reg surf for target subject --------- */
  if (stringmatch(trgsubject,"ico")) {
    /* Use Icosahedron as target surface */
    printf("INFO: loading ico (order = %d)\n",IcoOrder);
    fflush(stdout);
    IcoSurf = ReadIcoByOrder(IcoOrder, IcoRadius);
  } else {
    /* Use target subject (still called IcoSurf) */
    IcoSurf = MRISloadSurfSubject(trgsubject,hemi,surfregid,SUBJECTS_DIR);
    if (IcoSurf == NULL) {
      printf("ERROR: could not load registration surface\n");
      exit(1);
    }
  }
  MRIScomputeMetricProperties(IcoSurf) ;
  printf("surface nvertices %d\n",IcoSurf->nvertices);
  printf("surface area %f\n",IcoSurf->total_area);


  /*--- Load previously computed beta and res err var ----*/
  if (beta_in_id != NULL) {
    printf("INFO: loading beta_in %s\n",beta_in_id);
    fflush(stdout);
    if (beta_in_fmt != NULL) beta = MRIreadType(beta_in_id,beta_in_fmtid);
    else                    beta = MRIread(beta_in_id);
    if (beta == NULL) {
      printf("ERROR: loading %s\n",beta_in_id);
      exit(1);
    }
    if (beta->nframes != nregressors) {
      printf("ERROR: number of components in input beta does not match\n");
      printf("       the number of regressors in the design matrix\n");
      exit(1);
    }

    printf("INFO: loading var_in %s\n",eresvar_in_id);
    if (eresvar_in_fmt != NULL) {
      eresvar = MRISloadSurfVals(eresvar_in_id,eresvar_in_fmt,
                                 IcoSurf,trgsubject,hemi,SUBJECTS_DIR);
      //eresvar = MRIreadType(eresvar_in_id,eresvar_in_fmtid);
    } else
      eresvar = MRIread(eresvar_in_id);
    if (eresvar == NULL) {
      printf("ERROR: loading %s\n",eresvar_in_id);
      exit(1);
    }
    if (eresvar->nframes != 1) {
      printf("ERROR: number of frames in input var does not equal 1.\n");
      exit(1);
    }

    if (beta->width  != eresvar->width ||
        beta->height != eresvar->height ||
        beta->depth  != eresvar->depth) {
      printf("ERROR: dimension mismatch between input beta and var.\n");
      exit(1);
    }
  }

  /* Alloc enough data for the data for all the subjects */
  if (beta_in_id == NULL)
    SrcVals = MRIallocSequence(IcoSurf->nvertices, 1, 1,MRI_FLOAT,nsubjects);

  /* ----------------- Load real data --------------- */
  if (beta_in_id == NULL && SynthPDF == 0) {

    for (nthsubj = 0; nthsubj < nsubjects; nthsubj++) {

      subject = subjectlist[nthsubj];
      printf("nthsubj = %d/%d, %s ---------------------------------------\n",
             nthsubj+1,nsubjects,subject);
      fflush(stdout);

      /* Read in the spherical registration surface */
      SurfReg = MRISloadSurfSubject(subject,hemi,surfregid,SUBJECTS_DIR);
      if (SurfReg == NULL) {
        printf("ERROR: could not load registration surface\n");
        exit(1);
      }
      strcpy(SurfReg->subject_name,subject);
      //SurfReg->hemi = hemi;

      if (surfmeasure != NULL) {
        if (stringmatch(inputfmt,"paint") || stringmatch(inputfmt,"w") ||
            stringmatch(inputfmt,"wfile")) {
          sprintf(tmpstr,"%s/%s/surf/%s.%s",
                  SUBJECTS_DIR,subject,hemi,surfmeasure);
          inputfname = tmpstr;
        } else  inputfname = surfmeasure;
      } else inputfname = inputlist[nthsubj];

      /* Read in the input for this subject. */
      /* If the input is a volume, then read in the registration file,
      check the subject name, load white surface, load thickness
      for projection, resample to the white surface, save in tmpmri */

      printf("  INFO: loading input %s as %s\n",inputfname,inputfmt);
      if (stringmatch(inputfmt,"curv"))
        tmpmri = MRISloadSurfVals(inputfname,"curv",SurfReg,
                                  subject,hemi,SUBJECTS_DIR);
      else if (stringmatch(inputfmt,"paint") || stringmatch(inputfmt,"w") ||
               stringmatch(inputfmt,"wfile")) {
        if (MRISreadValues(SurfReg,inputfname)) {
          printf("ERROR: reading input %s\n",inputfname);
          exit(1);
        }
        tmpmri = MRIcopyMRIS(NULL, SurfReg, 0, "val");
      } else {
        if (inputfmt != NULL) tmpmri = MRIreadType(inputfname,inputfmtid);
        else                 tmpmri = MRIread(inputfname);
        nvoxels = tmpmri->width * tmpmri->height * tmpmri->depth;
        if (nvoxels != SurfReg->nvertices) {
          printf("ERROR: number of vertices in input file (%d) != number of vertices on surface (%d)\n",
                 nvoxels,SurfReg->nvertices);
          exit(1);
        }
        if (tmpmri->height > 1 || tmpmri->depth > 1) {
          printf("Reshaping\n");
          tmpmri2 = mri_reshape(tmpmri, SurfReg->nvertices, 1, 1,tmpmri->nframes);
          MRIfree(&tmpmri);
          tmpmri = tmpmri2;
        }
      }
      if (tmpmri == NULL) {
        printf("ERROR: could not load %s\n",inputfname);
        exit(1);
      }

      /* Extract frame (Future: compute derived variable) */
      if (tmpmri->nframes <= frame) {
        printf("ERROR: nframes (%d) <= frame (%d)\n",
               tmpmri->nframes,frame);
        exit(1);
      }

      if (frame > 0) {
        /* extract frame */
        tmpmri2 = MRIallocSequence(SurfReg->nvertices,1,1,MRI_FLOAT,1);
        for (vtx=0; vtx < SurfReg->nvertices; vtx++)
          MRIFseq_vox(tmpmri2,vtx,0,0,0) = MRIFseq_vox(tmpmri,vtx,0,0,frame);
        MRIfree(&tmpmri);
        tmpmri = tmpmri2;
      }

      /* Smooth on the native surface */
      if (nsmooth > 0)
        MRISsmoothMRI(SurfReg, tmpmri, nsmooth, NULL, tmpmri);

      /*------- Resample to target subject -------------------*/
      if (!stringmatch(trgsubject,subject)) {
        printf("  INFO: resampling to %s\n",trgsubject);
        fflush(stdout);
        tmpmri2 = surf2surf_nnfr(tmpmri, SurfReg, IcoSurf,
                                 &SrcHits,&SrcDist,&TrgHits,&TrgDist,1,1);
        if (tmpmri2 == NULL) {
          printf("ERROR: could not resample to %s.\n",trgsubject);
          exit(1);
        }
        MRIfree(&SrcHits);
        MRIfree(&SrcDist);
        MRIfree(&TrgHits);
        MRIfree(&TrgDist);
        MRIfree(&tmpmri);
        tmpmri = tmpmri2;
      }

      /* Copy into SrcVals structure */
      for (vtx = 0; vtx < IcoSurf->nvertices; vtx++)
        MRIFseq_vox(SrcVals,vtx,0,0,nthsubj)=MRIFseq_vox(tmpmri2,vtx,0,0,0);

      MRIfree(&tmpmri);
      MRISfree(&SurfReg);

      printf("\n\n");
      fflush(stdout);
    }
  } /* End: load real data */


  if (SynthPDF != 0) {
    printf("INFO: synthesizing, pdf = %d, seed = %d\n",SynthPDF,SynthSeed);
    printf("INFO: simulating over %d loops\n",nsim);
  }

  /* ------- Start simulation loop -- only one pass for non-sim runs -------- */
  for (nthsim = 0; nthsim < nsim; nthsim++) {

    if (MCSim) {
      printf("%d/%d\n",nthsim,nsim);
      fflush(stdout);
    }

    /* ------ Synthesize data ------ */
    if (beta_in_id == NULL && SynthPDF != 0) {
      if (SynthPDF == 1)
        MRIrandn(IcoSurf->nvertices, 1, 1,nsubjects,
                 SynthGaussianMean, SynthGaussianStd, SrcVals);
      if (SynthPDF == 2)
        MRIsampleCDF(IcoSurf->nvertices, 1, 1,nsubjects,
                     SynthXCDF, SynthCDF, SynthNCDF,SrcVals);

      if (nsmooth > 0) MRISsmoothMRI(IcoSurf, SrcVals, nsmooth, NULL, SrcVals);
    } /*End syntheisze */

    if (abs_flag) SrcVals = MRIabs(SrcVals,SrcVals);

    /* Save the input data */
    if (beta_in_id == NULL && yid != NULL && MCSim == 0) {
      if (MRIwriteAnyFormat(SrcVals,yid,yfmt,-1,NULL)) exit(1);
      if (fsgd != NULL) {
        sprintf(tmpstr,"%s.fsgd",getstem(yid));
        if (surfmeasure != NULL) strcpy(fsgd->measname,surfmeasure);
        else                 strcpy(fsgd->measname,"external");

        yidbase     = fio_basename(yid,NULL);
        yidstem     = getstem(yid);
        yidbasestem = getstem(yidbase);
        sprintf(fsgd->datafile,"%s",yidbase);

        sprintf(fsgd->DesignMatFile,"%s.X.mat",yidbasestem);
        sprintf(xmatpath,"%s.X.mat",yidstem);
        MatlabWrite(X,xmatpath,"X");

        fp = fopen(tmpstr,"w");
        gdfPrintHeader(fp,fsgd);
        fprintf(fp,"Creator          %s\n",Progname);
        fprintf(fp,"SmoothSteps      %d\n",nsmooth);
        fprintf(fp,"SUBJECTS_DIR     %s\n",SUBJECTS_DIR);
        fprintf(fp,"SynthSeed        %d\n",SynthSeed);
        fclose(fp);
      }
    }/* End save the input data */


    /*-------- Compute beta ----------- */
    if (beta_in_id == NULL) {

      if (DoPermute) {
        MatrixRandPermRows(X);
        printf("Permuting rows of design matrix (seed = %d)\n",SynthSeed);
        MatrixPrint(stdout,X);
      }

      if (nthsim == 1) printf("INFO: computing beta \n");
      fflush(stdout);
      beta = fMRImatrixMultiply(SrcVals, Q, beta);
      if (betaid != NULL && MCSim == 0)
        if (MRIwriteAnyFormat(beta,betaid,betafmt,-1,NULL)) exit(1);

      if (nthsim == 1) printf("INFO: computing eres \n");
      fflush(stdout);
      eres = fMRImatrixMultiply(SrcVals, R, eres);
      if (eresid != NULL && MCSim == 0)
        if (MRIwriteAnyFormat(eres,eresid,eresfmt,-1,NULL)) exit(1);

      if (yhatid != NULL && MCSim == 0) {
        if (nthsim == 1) printf("INFO: computing yhat \n");
        fflush(stdout);
        yhat = fMRImatrixMultiply(SrcVals, T, yhat);
        if (MRIwriteAnyFormat(yhat,yhatid,yhatfmt,-1,NULL)) exit(1);
        MRIfree(&yhat);
      }

      if (nthsim == 1) printf("INFO: computing var \n");
      fflush(stdout);
      //eresvar = fMRIvariance(eres,DOF,0,eresvar);
      eresvar = fMRIcovariance(eres,0,eres->nframes-DOF,NULL,eresvar);
      if (eresvarid != NULL && MCSim == 0)
        if (MRIwriteAnyFormat(eresvar,eresvarid,eresvarfmt,0,IcoSurf)) exit(1);
      MRIfree(&eres);
    }/* end compute beta */


    /* Compute contrast-effect size */
    if (tid != NULL || sigid != NULL || cesid != NULL ||
        tmaxfile != NULL || SynthPDF != 0) {
      if (nthsim == 1) printf("INFO: computing contrast effect size \n");
      ces = fMRImatrixMultiply(beta,C,ces);
      printf("ces nframes = %d\n",ces->nframes);
      if (cesid != NULL && MCSim == 0) {
        if (IsSurfFmt(cesfmt)) {
          if (IcoSurf == NULL)
            IcoSurf = MRISloadSurfSubject(trgsubject,hemi,surfregid,SUBJECTS_DIR);
          MRIwriteAnyFormat(ces,cesid,cesfmt,0,IcoSurf);
        } else {
          MRIwriteAnyFormat(ces,cesid,cesfmt,-1,NULL);
        }

      }
    }

    /* Compute t-ratio  */
    if (tid != NULL || sigid != NULL || tmaxfile != NULL || SynthPDF != 0) {
      if (nthsim == 1) printf("INFO: computing t \n");
      if (C->rows == 1) t = fMRIcomputeT(ces, X, C, eresvar, t);
      else         t = fMRIcomputeF(ces, X, C, eresvar, t);
      if (tid != NULL && MCSim == 0) {
        if (IsSurfFmt(tfmt) && IcoSurf == NULL)
          IcoSurf = MRISloadSurfSubject(trgsubject,hemi,surfregid,SUBJECTS_DIR);
        if (MRIwriteAnyFormat(t,tid,tfmt,0,IcoSurf)) exit(1);
      }
      if (tmaxfile != NULL) {
        tmax = fabs(MRIFseq_vox(t,0,0,0,0));
        for (vtx = 0; vtx < t->width; vtx++) {
          if (IcoSurf->vertices[vtx].ripflag) continue;
          if (tmax < fabs(MRIFseq_vox(t,vtx,0,0,0)) )
            tmax = fabs(MRIFseq_vox(t,vtx,0,0,0));
        }
        fp = fopen(tmaxfile,"a");
        fprintf(fp,"%f %f\n",tmax,-log10(sigt(tmax,DOF)));
        fclose(fp);
      }
    }

    /* Compute significance of t-ratio  */
    if (sigid != NULL || SynthPDF != 0) {
      if (nthsim == 1) printf("INFO: computing t significance \n");
      if (C->rows == 1)
        sig = fMRIsigT(t, DOF, sig);
      else {
        printf("Computing sigF\n");
        sig = fMRIsigF(t, DOF, C->rows, sig);
      }
      MRIlog10(sig,NULL,sig,1);
      //if(sigfmt != NULL && MCSim == 0){
      if (!MCSim) {
        if (sigid != NULL) {
          if (IsSurfFmt(sigfmt) && IcoSurf == NULL)
            IcoSurf = MRISloadSurfSubject(trgsubject,hemi,surfregid,SUBJECTS_DIR);
          if (MRIwriteAnyFormat(sig,sigid,sigfmt,0,IcoSurf)) exit(1);
        }
      }
    }

    /* Compute number of clusters for simulation */
    if (MCSim) {
      /* Load the log10p-values into the surface */
      for (vtx = 0; vtx < IcoSurf->nvertices; vtx++)
        IcoSurf->vertices[vtx].val = MRIFseq_vox(sig,vtx,0,0,0);

      if (nthsim == 0) {
        /* Set up Cluster Hit Table */
        cht->nsim = 0;
        cht->nvox = 0;
        cht->nsmooth = nsmooth;
        cht->fwhm = 0;
        cht->totsize = IcoSurf->total_area;
        cht->seed = SynthSeed;
        CHTwrite(chtfile,cht); // Immediately create cht file with zeros
      } else cht = CHTread(chtfile);

      printf("%d Searching for Clusters\n",nthsim+1);
      for (nth_ithr=0; nth_ithr < cht->n_ithr; nth_ithr++ ) {
        ithr = cht->ithr[nth_ithr];
        for (nth_sthr=0; nth_sthr < cht->n_sthr; nth_sthr++ ) {
          sthr = cht->sthr[nth_sthr];
          scs = sclustMapSurfClusters(IcoSurf,ithr,-1,cht->ithr_signid,
                                      sthr,&NClusters,NULL);
          printf("  %6.2f %6.2f %3d \n",ithr,sthr,NClusters);
          if (NClusters > 0) cht->hits[nth_ithr][nth_sthr] ++;
          free(scs);
        }
      }
      cht->nsim ++;
      CHTwrite(chtfile,cht);
    }

  } /*End simulation loop*/
  if (nsim > 1) {
    printf("%d \n",nthsim);
    CHTprint(stdout, cht);
  }

  if (SrcVals) MRIfree(&SrcVals);

  printf("done \n");
  return(0);
}
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  extern MATRIX *X;
  int  nargc , nargsused;
  char **pargv, *option ;
  int m, err;
  float fvtmp[1000];

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
    else if (!strcasecmp(option, "--force"))   Force = 1;
    else if (!strcasecmp(option, "--parseonly")) ParseOnly = 1;
    else if (!strcasecmp(option, "--allowsubjrep"))
      fsgdf_AllowSubjRep = 1; /* external, see fsgdf.h */
    else if ( !strcmp(option, "--xmatonly") ) xmatonly = 1;
    else if ( !strcmp(option, "--abs") ) abs_flag = 1;
    else if ( !strcmp(option, "--permute") ) DoPermute = 1;

    else if (!strcmp(option, "--seed")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&SynthSeed);
      nargsused = 1;
    } else if (!strcmp(option, "--gaussian")) {
      if (nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%lf",&SynthGaussianMean);
      sscanf(pargv[1],"%lf",&SynthGaussianStd);
      SynthPDF = 1;
      nargsused = 2;
    } else if (!strcmp(option, "--cdf")) {
      if (nargc < 1) argnerr(option,1);
      SynthCDFFile = pargv[0];
      err = PDFloadCDF(SynthCDFFile, &SynthXCDF, &SynthCDF, &SynthNCDF);
      if (err) exit(1);
      SynthPDF = 2;
      nargsused = 1;
    } else if (!strcmp(option, "--icoorder")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&IcoOrder);
      nargsused = 1;
    } else if (!strcmp(option, "--hemi")) {
      if (nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--sd")) {
      if (nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--trgsubj") || !strcmp(option, "--ts") ) {
      if (nargc < 1) argnerr(option,1);
      trgsubject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--surfmeas")) {
      if (nargc < 1) argnerr(option,1);
      surfmeasure = pargv[0];
      inputfmt = "curv";
      inputfmtid = checkfmt(inputfmt);
      nargsused = 1;
    } else if (!strcmp(option, "--i")) {
      if (nargc < 2) argnerr(option,2);
      nargsused = 0;
      while ( nth_is_arg(nargc, pargv, nargsused) ) {
        inputlist[ninputs] = pargv[nargsused];
        nargsused ++;
        ninputs++;
      }
      printf("INFO: found %d input files on cmdline \n",ninputs);
    } else if (!strcmp(option, "--ifile")) {
      if (nargc < 2) argnerr(option,2);
      nargsused = 1;
      fp = fopen(pargv[0],"r");
      if (fp==NULL) {
        printf("ERROR: could not open %s\n",pargv[0]);
        exit(1);
      }
      while (fscanf(fp,"%s",tmpstr) != EOF) {
        inputlist[ninputs] = (char *) calloc(strlen(tmpstr)+1,sizeof(char));
        memmove(inputlist[ninputs],tmpstr,strlen(tmpstr));
        ninputs++;
      }
      fclose(fp);
      printf("INFO: found %d input files in %s\n",ninputs,pargv[0]);
    } else if (!strcmp(option, "--ifmt")) {
      if (nargc < 1) argnerr(option,1);
      inputfmt = pargv[0];
      inputfmtid = checkfmt(inputfmt);
      nargsused = 1;
    } else if ( !strcmp(option, "--fsgd") ) {
      if (nargc < 1) argnerr(option,1);
      fsgdfile = pargv[0];
      nargsused = 1;
      fsgd = gdfRead(fsgdfile,0);
      if (fsgd==NULL) exit(1);
      strcpy(fsgd->tessellation,"surface");
      if (nth_is_arg(nargc, pargv, 1)) {
        gd2mtx_method = pargv[1];
        nargsused ++;
        if (gdfCheckMatrixMethod(gd2mtx_method)) exit(1);
      } else gd2mtx_method = "dods";
      printf("INFO: gd2mtx_method is %s\n",gd2mtx_method);
      if (!stringmatch(gd2mtx_method,"none")) {
        X = gdfMatrix(fsgd,gd2mtx_method,NULL);
        CheckDesignMatrix(X);
        nsubjects = X->rows;
        nregressors = X->cols;
        for (m=0; m<nsubjects; m++) subjectlist[m] = fsgd->subjid[m];
      }
      strcpy(fsgd->DesignMatMethod,gd2mtx_method);
    } else if ( !strcmp(option, "--design") ) {
      if (nargc < 1) argnerr(option,1);
      desmtxfname = pargv[0];
      nargsused = 1;
      ReadDesignMatrix(desmtxfname);
      CheckDesignMatrix(X);
    } else if ( !strcmp(option, "--xmat") ) {
      if (nargc < 1) argnerr(option,1);
      xmatfile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--xmatfmt") ) {
      if (nargc < 1) argnerr(option,1);
      xmatfmt = pargv[0];
      if (! (!strcmp(xmatfmt,"matlab4") || !strcmp(xmatfmt,"ascii")) ) {
        printf("ERROR: xmatfmt = %s, must be matlab4 or ascii\n",xmatfmt);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--nsmooth")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nsmooth);
      nargsused = 1;
    } else if (!strcmp(option, "--frame")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
    } else if (!strcmp(option, "--mcsim")) {
      if (nargc < 9) argnerr(option,9);
      MCSim = 1; // Monte Carlo Simulation
      sscanf(pargv[0],"%d",&nsim);
      chtfile = pargv[1];
      sscanf(pargv[2],"%d",&n_ithr);
      sscanf(pargv[3],"%lf",&ithr_lo);
      sscanf(pargv[4],"%lf",&ithr_hi);
      ithr_sign = pargv[5];
      sscanf(pargv[6],"%d",&n_sthr); // Area threshold is in mm
      sscanf(pargv[7],"%lf",&sthr_lo);
      sscanf(pargv[8],"%lf",&sthr_hi);
      cht = CHTalloc(n_ithr,ithr_lo,ithr_hi, n_sthr,sthr_lo,sthr_hi);
      if (cht==NULL) {
        printf("ERROR: with mcsim params\n");
        exit(1);
      }
      if (CHTsetSignString(cht, ithr_sign) == -100) {
        printf("ERROR: with mcsim params\n");
        exit(1);
      }
      nargsused = 9;
    } else if (!strcmp(option, "--contrast")) {
      if (nargc < 1) argnerr(option,1);
      conmtxfname = pargv[0];
      C = ReadAsciiMatrix(conmtxfname);
      if (C==NULL) exit(1);
      nargsused = 1;
    } else if (!strcmp(option, "--gcv")) {
      if (nargc < 1) argnerr(option,1);
      nargsused = 0;
      while ( nth_is_arg(nargc, pargv, nargsused) ) {
        sscanf(pargv[nargsused],"%f",&fvtmp[nargsused]);
        nargsused ++;
      }
      printf("INFO: found %d elements in group contrast vector\n",
             nargsused);
      C = MatrixAlloc(1,nargsused,MATRIX_REAL);
      for (m=0; m < nargsused; m++) C->rptr[1][m+1] = fvtmp[m];
      MatrixPrint(stdout,C);
    } else if (!strcmp(option, "--beta")) {
      if (nargc < 1) argnerr(option,1);
      betaid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        betafmt = pargv[1];
        nargsused ++;
        betafmtid = checkfmt(betafmt);
      } else betafmtid = getfmtid(betaid);
      if (IsSurfFmt(betafmt)) {
        printf("ERROR: cannot use curv or paint as output format for beta\n");
        exit(1);
      }
    } else if (!strcmp(option, "--beta_in")) {
      if (nargc < 1) argnerr(option,1);
      beta_in_id = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        beta_in_fmt = pargv[1];
        nargsused ++;
        beta_in_fmtid = checkfmt(beta_in_fmt);
      }
    } else if (!strcmp(option, "--ces")) {
      if (nargc < 1) argnerr(option,1);
      cesid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        cesfmt = pargv[1];
        nargsused ++;
        cesfmtid = checkfmt(cesfmt);
      } else cesfmtid = getfmtid(cesid);
    } else if (!strcmp(option, "--eres")) {
      if (nargc < 1) argnerr(option,1);
      eresid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        eresfmt = pargv[1];
        nargsused ++;
        eresfmtid = checkfmt(eresfmt);
      } else eresfmtid = getfmtid(eresid);
    } else if (!strcmp(option, "--y")) {
      if (nargc < 1) argnerr(option,1);
      yid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        yfmt = pargv[1];
        nargsused ++;
        yfmtid = checkfmt(yfmt);
      } else yfmtid = getfmtid(yid);
      printf("y stem %s\n",getstem(yid));
    } else if (!strcmp(option, "--yhat")) {
      if (nargc < 1) argnerr(option,1);
      yhatid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        yhatfmt = pargv[1];
        nargsused ++;
        yhatfmtid = checkfmt(yhatfmt);
      } else yhatfmtid = getfmtid(yhatid);
    } else if (!strcmp(option, "--var")) {
      if (nargc < 1) argnerr(option,1);
      eresvarid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        eresvarfmt = pargv[1];
        nargsused ++;
        eresvarfmtid = checkfmt(eresvarfmt);
      } else eresvarfmtid = getfmtid(eresvarid);
      if (stringmatch(eresvarfmt,"curv")) {
        printf("ERROR: cannot use curv as output format for var\n");
        exit(1);
      }
    } else if (!strcmp(option, "--var_in")) {
      if (nargc < 1) argnerr(option,1);
      eresvar_in_id = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        eresvar_in_fmt = pargv[1];
        nargsused ++;
        eresvar_in_fmtid = checkfmt(eresvar_in_fmt);
      } else eresvar_in_fmtid = getfmtid(eresvar_in_id);
    } else if (!strcmp(option, "--t")) {
      if (nargc < 1) argnerr(option,1);
      tid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        tfmt = pargv[1];
        nargsused ++;
        tfmtid = checkfmt(tfmt);
      } else tfmtid = getfmtid(tid);
      if (stringmatch(tfmt,"curv")) {
        printf("ERROR: cannot use curv as output format for t\n");
        exit(1);
      }
    } else if (!strcmp(option, "--tmax")) {
      tmaxfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--sigt")) {
      if (nargc < 1) argnerr(option,1);
      sigid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        sigfmt = pargv[1];
        nargsused ++;
        sigfmtid = checkfmt(sigfmt);
      } else sigfmtid = getfmtid(sigid);
      if (stringmatch(sigfmt,"curv")) {
        printf("ERROR: cannot use curv as output format for sigt\n");
        exit(1);
      }
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
  printf("USAGE: mris_glm (formerly mri_surfglm)\n") ;
  printf("\n");
  printf("Raw Data Input Options\n");
  printf("   --fsgd    fname : FreeSurfer Group Descriptor File\n");
  printf("   --design  fname : name of design matrix (ascii)\n");
  printf("   --surfmeas name  : input file or name of surface measure\n");
  printf("   --i input1 input2 ...> : input file list\n");
  printf("   --frame      M     : use 0-based Mth frame (default is 0)\n");
  printf("   --hemi       hemi  : hemisphere (lh or rh) \n");
  printf("   --trgsubj    subject : target subject \n");
  printf("   --icoorder   order : order of icosahedral tesselation (default 7)\n");
  printf("   --nsmooth    N     : number of smoothing iterations\n");
  printf("   --abs   : use absolute value of input after smoothing\n");
  printf("\n");
  printf("Processed Data Input Options\n");
  printf("   --beta_in    name <fmt> : parameter estimates from previous \n");
  printf("   --var_in     name <fmt> : reserr var from previous \n");
  printf("\n");
  printf("Estimation Output Options\n");
  printf("   --y       name <fmt> : input data after resampling and smoothing\n");
  printf("   --beta    name <fmt> : parameter estimates \n");
  printf("   --var     name <fmt> : variance of residual error \n");
  printf("   --yhat    name <fmt> : signal estimate\n");
  printf("   --eres    name <fmt> : residual error \n");
  printf("   --xmat    matfile    : save design matrix in matlab4 format\n");
  printf("\n");
  printf("Group Contrast Options\n");
  printf("   --contrast   fname : file containing group contrast matrix (ascii)\n");
  printf("   --gcv c1 ... cN : group contrast vector specified on cmdline\n");
  printf("   --ces     name <fmt> : contrast effect size  \n");
  printf("   --t       name <fmt> : t-ratio of contrast \n");
  printf("   --sigt    name <fmt> : signficance of t-ratio (ie, t-Test) \n");
  printf("\n");
  printf("Synthesis and Simulation Options\n");
  printf("   --gaussian mean std : synthesize data with guassian\n");
  printf("   --cdf cdffile : synthesize data with given cdf\n");
  printf("   --seed seed : random number seed when synth data (-1 for auto)\n");
  printf("   --tmax fname : append max t (and min p) to fname \n");
  printf("   --mcsim nsim fname nithr ithrlo ithrhi ithrsign nsthr sthrlo sthrhi: "
         "MC Simulation.\n");
  printf("\n");
  printf("   --xmatonly : only compute and save the design matrix\n");
  printf("   --allowsubjrep : allow subject name to be replicated\n");
  printf("   --force : force processing with badly cond X\n");
  printf("   --sd    subjectsdir : default is env SUBJECTS_DIR\n");
  printf("\n");
  printf("   --version : print version and exit\n");
  printf("   --help : a short story.\n");
  printf("\n");
  printf("%s\n",vcid);
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf(
    "SUMMARY\n"
    "\n"
    "NOTE: this program replaces mri_surfglm (this is just a name change).\n"
    "\n"
    "This program performs inter-subject/group averaging and inference on\n"
    "the surface by fitting a GLM model at each vertex. The model consists\n"
    "of subject parameters (eg, age, gender, etc). The model is the same\n"
    "across all vertices, though the fit may be (will be) different. The\n"
    "user must supply a matrix that represents the GLM. While estimation\n"
    "and inference can be performed in a single call to mris_glm, the tasks\n"
    "can also be separated, which can be much more convenient. Inferences\n"
    "are not corrected for multiple comparisons.\n"
    "\n"
    "MATHEMATICAL BACKGROUND\n"
    "\n"
    "The forward model is given by:\n"
    "\n"
    "    y = XB + n\n"
    "\n"
    "where X is the Ns-by-Nb design matrix, y is the Ns-by-Nv raw data set,\n"
    "B is the Nb-by-Nv regression parameters, and n is noise. Ns is the\n"
    "number of subjects, Nb is the number of regressors, and Nv is the \n"
    "number of vertices. y will have been preprocessed in possibly two\n"
    "ways: (1) it will be sampled on the surface of the target subject, \n"
    "and (2) it may be spatially smoothed (prior to resampling) (see\n"
    "--nsmooth). \n"
    "\n"
    "During the estimation stage, the forward model is inverted to\n"
    "solve for B:\n"
    "\n"
    "    B = inv(X'*X)*X'y\n"
    "\n"
    "This is performed at each vertex, the result of which can be saved\n"
    "with the --beta flag.\n"
    "\n"
    "The signal estimate (which can be saved with --yhat) is computed as \n"
    "\n"
    "    yhat = B*X\n"
    "\n"
    "The residual error (which can be saved with --eres) is computed as \n"
    "\n"
    "    eres = y - yhat\n"
    "\n"
    "The noise variance estimate (computed for each vertex) is computed\n"
    "as the sum of the squares of the residual error divided by the DOF.\n"
    "The DOF equals the number of rows of X minus the number of columns.\n"
    "The noise variance can be saved with --var.\n"
    "\n"
    "A contrast vector C has as many elements as columns of X. The \n"
    "contrast effect size (--ces) is then computed as:\n"
    "\n"
    "   G = C*B\n"
    "\n"
    "The t-ratio (--t) for the contrast is then given by:\n"
    "\n"
    "   t = G/sqrt(var * C*inv(X'X)*C')\n"
    "\n"
    "The signifiance of the t-ratio (based on a double-sided t-test\n"
    "and uncorrected for multiple comparisons) can be saved with \n"
    "the --sigt flag.\n"
    "\n"
    "\n"
    "COMMAND-LINE ARGUMENTS\n"
    "\n"
    "--fsgd fname <gd2mtx>\n"
    "\n"
    "Specify the design with a FreeSurfer Group Descriptor File (FSGDF).\n"
    "See http://surfer.nmr.mgh.harvard.edu/docs/fsgdf.txt for more info.\n"
    "The gd2mtx is the method by which the group description is converted\n"
    "into a design matrix. Legal values are doss (Different Offset, Same\n"
    "Slope) and dods (Different Offset, Different Slope). doss will create\n"
    "a design matrix in which each class has it's own offset but forces all\n"
    "classes to have the same slope. dods models each class with it's own\n"
    "offset and slope. In either case, you'll need to know the order of the\n"
    "regressors in order to correctly specify the contrast vector. For\n"
    "doss, the first NClass columns refer to the offset for each class.\n"
    "The remaining columns are for the continuous variables. In dods, the\n"
    "first NClass columns again refer to the offset for each class.\n"
    "However, there will be NClass*NVar more columns (ie, one column for\n"
    "each variable for each class). The first NClass columns are for the\n"
    "first variable, etc. If neither of these models works for you, you\n"
    "will have to specify the design matrix manually. In that case, you can\n"
    "still specify an FSGDF but use 'none' as the gd2mtx. It can be\n"
    "advantageous to specify an FSGDF because this information will be\n"
    "propagated with the output data and can be used by tksurfer to display\n"
    "scatter plots.\n"
    "\n"
    "--design fname\n"
    "\n"
    "File name for design matrix. The design matrix must be in an ASCII\n"
    "file. The first column must be the id assigned to the subject during\n"
    "FreeSurfer reconstruction. The following columns are the the design\n"
    "matrix. It is possible for the design matrix to be ill-conditioned.\n"
    "This means that two columns are identical or that one column is equal\n"
    "to a weighted sum of any of the other columns.  In this case, the\n"
    "matrix cannot be inverted and so the analysis must stop. The matrix is\n"
    "judged to be ill-conditioned if its condition number of the normalized\n"
    "matrix is greater than 100000. The X matrix is normalized such that\n"
    "the length of each column vector is 1. The threshold of 100000 is\n"
    "somewhat arbitrary. It is possible to test the condition of a matrix\n"
    "without having to include all the command-line options needed for a\n"
    "full analysis. Just run:\n"
    "          mris_glm --design fname \n"
    "It will print out an INFO line with the condition number. You can\n"
    "force processing with the given design matrix by PRECEDING the\n"
    "--design flag with --force. \n"
    "\n"
    "--surfmeas name \n"
    "\n"
    "This is one of the two ways the user can specify the raw data (ie, the\n"
    "input to the estimation stage). This method requires that the data\n"
    "file for each subject reside in the FreeSurfer anatomical directory\n"
    "under the surf subdirectory. This frees the user from having to list\n"
    "all the inputs on the command-line. The name can be one of two things,\n"
    "depending upon the format designation (--ifmt). If the format is curv\n"
    "(or unspecified), then mris_glm will construct the name of the\n"
    "input file as hemi.name (where hemi is either lh or rh as specified \n"
    "by --hemi). If the format is anything else, then it looks for a file\n"
    "called name in the surf subdirectory. Only specify a raw data input\n"
    "when performing estimation.\n"
    "\n"
    "--i input1 input2 ...\n"
    "\n"
    "This is second method that the user can specify the raw data (ie, the\n"
    "input to the estimation stage). This method allows the user to specify\n"
    "all of the input files on the command line. There must be as many\n"
    "files listed as there are rows in the design matrix. The format (same\n"
    "for all inputs) is as designated with the --ifmt flag. If the format\n"
    "is unspecified, mris_glm will attempt to determine the format.\n"
    "curv format cannot be used with explicit inputs. Only specify a raw \n"
    "data input when performing estimation.\n"
    "\n"
    "--ifmt format\n"
    "\n"
    "This flag is used to specify the format of the input raw data files.\n"
    "When a surface measure is specified, the default format becomes curv,\n"
    "otherwise the mris_glm will attempt to determine the format from\n"
    "the file name if a format is not explicitly give. It is not possible\n"
    "to use curv format with explicit inputs (ie, --i). Valid formats are:\n"
    "bshort, bfloat, COR, analyze, spm, paint (or w or wfile), and curv.\n"
    "\n"
    "--frame M\n"
    "\n"
    "This allows the user to specify that the Mth frame of the raw data (if \n"
    "the input format supports multiple frames) should be used as input. \n"
    "M is zero-based. Default is 0.\n"
    "\n"
    "--hemi hemisphere\n"
    "\n"
    "Specify that the input data either the left (lh) or right (rh) hemisphere.\n"
    "\n"
    "--trgsubject subject\n"
    "\n"
    "Resample each input data set from the individual's surface to that of the\n"
    "target subject. \n"
    "\n"
    "--icoorder order\n"
    "\n"
    "When the target subject is ico, this specifies the order of the \n"
    "icosahedron. Default is 7 (163842 vertices). \n"
    "\n"
    "--nsmooth N\n"
    "\n"
    "Perform N iterations of nearest-neighbor spatial smoothing. Note: \n"
    "smoothing is  performed on the surface of each source subject before \n"
    "resampling to the target subject.\n"
    "\n"
    "--beta_in betaname <fmt>\n"
    "\n"
    "This flag (with --var_in) allows the user to use a previous estimation \n"
    "as input to the contrast/inference computation. This arguments should \n"
    "be identical to that specified with --beta when doing the estimation.\n"
    "Note: the user must also specify a pointer to the residual error\n"
    "variance using --var_in.\n"
    "\n"
    "--var_in betaname <fmt>\n"
    "\n"
    "This flag (with --beta_in) allows the user to use a previous estimation \n"
    "as input to the contrast/inference computation. This arguments should be\n"
    "identical to that specified with --var when doing the estimation.\n"
    "\n"
    "--y name <fmt>\n"
    "\n"
    "Save the raw data (after resampling and smoothing) into a single 'volume'.\n"
    "fmt is the format (see OUTPUT FORMATS). If an FSGDF has been specified\n"
    "(with --fsgdf), then this file is copied to name.fsgdf. This can be\n"
    "useful for displaying scatter plots in tksurfer.\n"
    "\n"
    "--beta name <fmt>\n"
    "\n"
    "Save the result of estimation (ie, the map regression coefficients). \n"
    "This can be used in subsequent calls to mris_glm to perform.\n"
    "inference. fmt is the format (see OUTPUT FORMATS).\n"
    "\n"
    "--var name <fmt>\n"
    "\n"
    "Save the estimate of the noise variances (estimated as the variance\n"
    "of the residual error). This can be used in subsequent calls to \n"
    "mris_glm to perform inference.\n"
    "\n"
    "--yhat name <fmt>\n"
    "\n"
    "Save the estimate of the signal. This is only good for debugging.\n"
    "fmt is the format (see OUTPUT FORMATS).\n"
    "\n"
    "--eres name <fmt>\n"
    "\n"
    "Save the residual error. This is only good for debugging.\n"
    "fmt is the format (see OUTPUT FORMATS).\n"
    "\n"
    "--xmat name \n"
    "\n"
    "Save the design matrix. By default, it is saved in matlab4 format. \n"
    "It can be saved in ASCII with --xmatfmt ascii.\n"
    "\n"
    "--xmatfmt format \n"
    "\n"
    "Save the design matrix in the specfied format. format can be \n"
    "matlab4 or ascii. Default is matlab4\n"
    "\n"
    "--xmatonly\n"
    "\n"
    "Save design matrix to output file indicated by --xmat and exit.\n"
    "This is a means of simply creating a design matrix that can \n"
    "be used for other purposes. Note that only --fsgd and --xmat are\n"
    "needed with --xmatonly (ie, you do not need to supply a full.\n"
    "mris_glm command-line).\n"
    "\n"
    "--contrast fname\n"
    "\n"
    "Load the group contrast vector from the ascii file fname. The \n"
    "contrast vector should have as many entries as there are columns\n"
    "in the design matrix. The contrast vector can also be specified \n"
    "on the command-line with --gcv.\n"
    "\n"
    "--gcv c1 ... cN\n"
    "\n"
    "Specify the group contrast vector (gcv) on the command-line as the\n"
    "values c1 ... cN. The contrast vector should have as many entries as \n"
    "there are columns in the design matrix. The contrast vector can also \n"
    "be specified in a file with --contrast.\n"
    "\n"
    "--ces name <fmt>\n"
    "\n"
    "Save the contrast effect size from the contrast vector specified\n"
    "with --contrast or --gcv. fmt is the format (see OUTPUT FORMATS).\n"
    "\n"
    "--t name <fmt>\n"
    "\n"
    "Save the t-ratio of the contrast.\n"
    "\n"
    "--tmax filename\n"
    "\n"
    "Append the maximum t value and corresponding -log10(p) in text file \n"
    "filename. Good for simulations. Best when used with --guassian. Vertices\n"
    "with the ripflag set are ignored.\n"
    "\n"
    "--sigt name <fmt>\n"
    "\n"
    "Save the signficance (p-value) of the t-ratio of the contrast. The\n"
    "value is actually the -log10 of the significance with the same\n"
    "sign as the t-ratio from which it was computed. The significance\n"
    "is computed from a double-sided t-test and is NOT corrected for\n"
    "multiple comparisons across space. fmt is the format (see OUTPUT \n"
    "FORMATS).\n"
    "\n"
    "--force\n"
    "\n"
    "Force processing eventhough the design matrix condition number \n"
    "exceeds 10000.\n"
    "\n"
    "\n"
    "--gaussian mean std\n"
    "\n"
    "Synthesize input data using white gaussian noise with the given\n"
    "mean and stddev. Good for debugging and Monte Carlo simulations.\n"
    "\n"
    "--cdf cdffile\n"
    "\n"
    "Synthesize input data by drawing samples from the given CDF. The CDF\n"
    "text file has two columns. The first column is the x at which the cdf is \n"
    "sampled.  The second column is the value of the cdf. cdf[n] is the \n"
    "probability that the random number will be <= x[n]. Good for debugging\n"
    "and Monte Carlo simulations.\n"
    "\n"
    "--seed seed\n"
    "\n"
    "Seed for random number generator to use when synthesizing data. Good \n"
    "for debugging and Monte Carlo simulations. If no seed is given, one\n"
    "will automatically be chosen based on the time-of-day. This option has\n"
    "no effect unless --gaussian or --cdf is used.\n"
    "\n"
    "--sd subjectsdir\n"
    "\n"
    "Look for FreeSurfer reconstructions in subjectsdir. If unspecified,\n"
    "the SUBJECTS_DIR envionment variable will be used.\n"
    "\n"
    "--allowsubjrep\n"
    "\n"
    "Allow repetitions in the subject name in the fsgdf file. This can\n"
    "be usefull when the input is specified with --i.\n"
    "\n"
    "\n"
    "OUTPUT FORMATS:\n"
    "\n"
    "Output formats can be designated by specifying a format string\n"
    "following the the output name. Valid strings include bfloat,\n"
    "bshort, spm, analyze, analyze4d, COR, paint, w, and wfile. Paint,\n"
    "w, and wfile cannot be used for beta, y, yhat, or eres.\n"
    "\n"
    "\n"
    "EXAMPLES:\n"
    "\n"
    "1. Analyze thickness maps based on gender and age for 5 hypothetical\n"
    "subjects: subj1 (m, 22) , subj2 (m, 57), subj3 (f, 33), subj4 (f, 65), \n"
    "subj5 (m, 27). The design matrix would look something like:\n"
    "\n"
    "   subj1  1  0  22\n"
    "   subj2  1  0  57\n"
    "   subj3  0  1  33\n"
    "   subj4  0  1  65\n"
    "   subj5  1  0  27\n"
    "\n"
    "The first column is the name of the subject as it appears in the\n"
    "FreeSurfer SubjectsDir. The second and third columns categorically\n"
    "code gender (first column male, second column female). The last\n"
    "column codes age. Assume this matrix is stored in a file called\n"
    "genage.mtx\n"
    "\n"
    "  mris_glm --design genage.mtx --hemi lh --surfmeas thickness \n"
    "    --trgsubj average7 --nsmooth 50 --beta beta bfloat \n"
    "    --var var bfloat \n"
    "\n"
    "This will read the thickness maps for each of the subjects, smooth\n"
    "it with 50 iterations of nearest-neighbor smoothing, resample to\n"
    "the average7 subject and save the regression coeffients and\n"
    "noise variance, both in bfloat format No inference was performed.\n"
    "\n"
    "2. Test the data in Example 1 for an effect of age:\n"
    "\n"
    "  mris_glm --design genage.mtx --hemi lh --trgsubj average7 \n"
    "      --beta_in beta bfloat --var_in var bfloat \n"
    "      --gcv 0 0 1 --sigt ./age-sigt-lh.w paint\n"
    "\n"
    "3. Test the data in Example 1 for a difference between males \n"
    "and females with age regressed out:\n"
    "\n"
    "  mris_glm --design genage.mtx --hemi lh --trgsubj average7 \n"
    "      --beta_in beta bfloat --var_in var bfloat \n"
    "      --gcv 1 -1 0 --sigt ./gender-sigt-lh.w paint\n"
    "\n"
    "4. Perform the same analysis as done in Example 1, but use \n"
    "values that have been painted onto the surface of each subject\n"
    "(this could have come from an fMRI analysis):\n"
    "\n"
    "  mris_glm --design genage.mtx --hemi lh\n"
    "    --ifmt paint --i ./subj1-data-lh.w ./subj2-data-lh.w \n"
    "    ./subj3-data-lh.w ./subj4-data-lh.w ./subj5-data-lh.w\n"
    "    --trgsubj average7 --nsmooth 50 --beta beta bfloat \n"
    "    --var var bfloat \n"
    "\n"
    "BUGS\n"
    "\n"
    "No correction for multiple comparisons.\n"
    "\n"
    "BUG REPORTING\n"
    "\n"
    "If you want your bug report or question to have a prayer of being\n"
    "answered, make sure to include the following information (send to\n"
    "freesurfer@surfer.nmr.mgh.harvard.edu): (1) version of mris_glm\n"
    "(run with --version), (2) full command-line used, (3) terminal\n"
    "output, (4) description of the problem. \n"
    "\n"


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
  if (SynthSeed < 0) SynthSeed = PDFtodSeed();
  srand48(SynthSeed);

  if (DoPermute) {
    printf("Permuting rows of design matrix (seed = %d)\n",SynthSeed);
    MatrixRandPermRows(X);
    MatrixPrint(stdout,X);
  }

  if (xmatonly) {
    if (fsgdfile == NULL) {
      printf("ERROR: need fsgd file with --matonly\n");
      exit(1);
    }
    if (xmatfile == NULL) {
      printf("ERROR: need xmat file with --matonly\n");
      exit(1);
    }
    //X = gdfMatrix(fsgd,gd2mtx_method,NULL); // X should already exist
    MatrixWriteFmt(X,xmatfile,xmatfmt);
    exit(0);
  }

  if (desmtxfname == NULL && fsgdfile == NULL) {
    printf("ERROR: must specify a design \n");
    exit(1);
  }

  if (desmtxfname != NULL && !stringmatch(gd2mtx_method,"none")) {
    printf("ERROR: cannot specify a design matrix and create matrix from fsgd\n");
    exit(1);
  }

  if ( (beta_in_id != NULL && eresvar_in_id == NULL) ||
       (beta_in_id == NULL && eresvar_in_id != NULL) ) {
    printf("ERROR: if using a previous estimate, must specify both\n");
    printf("       beta and var inputs\n");
    exit(1);
  }

  if (ninputs == 0 && surfmeasure == NULL && beta_in_id == NULL) {
    printf("ERROR: must specify an input, either surface measure or \n");
    printf("       a list of files or previous estimate (beta_in, var_in)\n");
    exit(1);
  }
  if (surfmeasure != NULL && beta_in_id != NULL) {
    printf("ERROR: cannot specify both surface measure and beta\n");
    exit(1);
  }
  if (ninputs != 0 && beta_in_id != NULL) {
    printf("ERROR: cannot specify both --i input and beta\n");
    exit(1);
  }
  if (surfmeasure != NULL && ninputs != 0) {
    printf("ERROR: cannot specify both surface measure and --i input\n");
    exit(1);
  }

  if (eresid != NULL && beta_in_id != NULL) {
    printf("ERROR: cannot specify an eres output with a previous estimate\n");
    printf("       as input.\n");
    exit(1);
  }

  if (yid != NULL && beta_in_id != NULL) {
    printf("ERROR: cannot specify a y output with a previous estimate.\n");
    printf("       as input.\n");
    exit(1);
  }

  if (ces != NULL || tid != NULL || sigid != NULL) {
    if (C == NULL) {
      printf("ERROR: must specify a contrast matrix with ces, t, or sig output\n");
      exit(1);
    }
    if (C != NULL) {
      if (C->cols != nregressors) {
        printf("ERROR: dimension mismatch: ncols in contrast matrix\n");
        printf("       is %d but the number of regressors is %d\n",
               C->cols,nregressors);
        MatrixPrint(stdout,C);
        exit(1);
      }
      if (C->rows != 1 && 0) {
        printf("ERROR: the contrast matrix can only have one row.\n");
        printf("Ask Doug to add F-test.\n");
        exit(1);
      }
    }
  }

  if (beta_in_id != NULL && cesid == NULL && tid == NULL && sigid == NULL) {
    printf("ERROR: nothing to do. You have specified a previous\n");
    printf("       estimate as input but have not specified an output.\n");
    printf("       The output should be ces, t, or sigt\n");
    exit(1);
  }

  if (trgsubject == NULL) {
    printf("ERROR: no target subject specified.\n");
    exit(1);
  }
  if (hemi == NULL) {
    printf("ERROR: no hemisphere specified.\n");
    exit(1);
  }


  if (SUBJECTS_DIR == NULL) {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL) {
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  if (beta_in_id == NULL) {
    if (hemi == NULL) {
      printf("ERROR: must specify a hemisphere.\n");
      exit(1);
    }
    if (strcmp(hemi,"lh") && strcmp(hemi,"rh")) {
      printf("ERROR: hemi = %s, must be either lh or rh\n",hemi);
      exit(1);
    }
  }

  if (ninputs != 0 && ninputs != X->rows) {
    printf("ERROR: number of inputs does not equal the number of rows.\n");
    printf("       in the design matrix\n");
    exit(1);
  }

  if (MCSim && C==NULL) {
    printf("ERROR: contrast vector needed with mcsim\n");
    exit(1);
  }

  if (MCSim && SynthPDF == 0) {
    /* Force it to used gaussian (0,1) */
    SynthPDF = 1;
    SynthGaussianMean = 0;
    SynthGaussianStd  = 1;
    printf("INFO: using gaussian (0,1) for simulation\n");
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
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

  /* check that there are enough args for nth to exist */
  if (nargc <= nth) return(0);

  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) return(0);

  return(1);
}
/*------------------------------------------------------------*/
static int stringmatch(char *str1, char *str2) {
  if (str1 == NULL && str2 != NULL) return(0);
  if (str2 == NULL && str1 != NULL) return(0);
  if (! strcmp(str1,str2)) return(1);
  return(0);
}
/*------------------------------------------------------------*/
static int getfmtid(char *fname) {
  int fmtid;
  fmtid = mri_identify(fname);
  if (fmtid == MRI_VOLUME_TYPE_UNKNOWN) {
    printf("ERROR: cannot determine format of %s\n",fname);
    exit(1);
  }
  return(fmtid);
}

/*------------------------------------------------------------*/
static int checkfmt(char *fmt) {
  int fmtid;

  if (fmt == NULL) return(MRI_VOLUME_TYPE_UNKNOWN);
  if (stringmatch(fmt,"curv") ||
      stringmatch(fmt,"paint") ||
      stringmatch(fmt,"wfile") ||
      stringmatch(fmt,"w")) return(MRI_VOLUME_TYPE_UNKNOWN);

  fmtid = string_to_type(fmt);
  if (fmtid == MRI_VOLUME_TYPE_UNKNOWN) {
    printf("ERROR: format string %s unrecognized\n",fmt);
    exit(1);
  }
  return(fmtid);
}

/*------------------------------------------------------------*/
static int IsSurfFmt(char *fmt) {
  if (fmt == NULL) return(0);
  if (stringmatch(fmt,"curv") ||
      stringmatch(fmt,"paint") ||
      stringmatch(fmt,"wfile") ||
      stringmatch(fmt,"w")) return(1);
  return(0);
}

/*------------------------------------------------------------*/
int ReadDesignMatrix(char *desmtxfname) {
  extern MATRIX *X;
  extern char *subjectlist[1000];
  extern int nsubjects, nregressors;
  int nrows = 0,ncols = 0, r,c;
  char tmpstring[1001];
  FILE *fp;

  ReadAsciiMatrixSize(desmtxfname, &nrows, &ncols);
  nsubjects = nrows;
  nregressors = ncols-1;
  X = MatrixAlloc(nrows,ncols-1,MATRIX_REAL);

  fp = fopen(desmtxfname,"r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n",desmtxfname);
    exit(1);
  }

  for (r=0; r < nrows; r++) {

    fscanf(fp,"%s",tmpstring);
    subjectlist[r] = (char *)calloc(strlen(tmpstring)+1,sizeof(char));
    memmove(subjectlist[r],tmpstring,strlen(tmpstring)+1);
    //printf("%2d %s\n",r+1,subjectlist[r]);

    for (c=0; c < ncols-1; c++)
      fscanf(fp,"%f",&(X->rptr[r+1][c+1]));

  }
  //MatrixPrint(stdout,X);

  fclose(fp);
  return(0);
}
/*------------------------------------------------------------*/
int ReadAsciiMatrixNRows(char *desmtxfname) {
  FILE *fp;
  int nrows;
  char tmpstring[2001];

  fp = fopen(desmtxfname,"r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n",desmtxfname);
    return(-1);
  }

  nrows = 0;
  while (fgets(tmpstring,2000,fp) != NULL)  nrows ++;
  fclose(fp);

  if (nrows == 0) {
    printf("ERROR: no data found in %s\n",desmtxfname);
    return(-1);
  }

  if (nrows < 0) {
    printf("ERROR: reading number of rows from %s\n",desmtxfname);
    return(-1);
  }

  return(nrows);
}
/*-----------------------------------------------------------------*/
int ReadAsciiMatrixSize(char *asciimtxfname, int *pnrows, int *pncols) {
  FILE *fp;
  int nrows, nitems, ncols;
  char tmpstring[2001];

  nrows = ReadAsciiMatrixNRows(asciimtxfname);
  if (nrows < 0) return(1);

  fp = fopen(asciimtxfname,"r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n",asciimtxfname);
    return(1);
  }

  nitems = 0;
  while (fscanf(fp,"%s",tmpstring) != EOF)  nitems ++;
  fclose(fp);

  if (nitems == 0) {
    printf("ERROR: no items found in %s\n",asciimtxfname);
    return(1);
  }

  ncols = nitems/nrows;
  if (ncols*nrows != nitems) {
    printf("ERROR: number of items (%d) not divisible by nrows (%d)\n",
           nitems,nrows);
    return(1);
  }

  if (ncols < -1) {
    printf("ERROR: reading number of cols from %s\n",desmtxfname);
    return(-1);
  }

  *pnrows = nrows;
  *pncols = ncols;

  return(0);
}
/*------------------------------------------------------------*/
MATRIX *ReadAsciiMatrix(char *asciimtxfname) {
  int err, nrows,ncols, r,c, nread;
  FILE *fp;

  err = ReadAsciiMatrixSize(asciimtxfname, &nrows, &ncols);
  if (err) return(NULL);

  C = MatrixAlloc(nrows,ncols,MATRIX_REAL);

  fp = fopen(asciimtxfname,"r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n",asciimtxfname);
    exit(1);
  }

  for (r=0; r < nrows; r++) {
    for (c=0; c < ncols; c++) {
      nread = fscanf(fp,"%f",&(C->rptr[r+1][c+1]));
      if (nread != 1) {
        printf("ERROR: ReadAsciiMatrix: could not read item %d,%d\n",r,c);
        MatrixFree(&C);
        return(NULL);
      }
    }
  }

  fclose(fp);
  return(C);
}
/*-----------------------------------------------*/
int CheckDesignMatrix(MATRIX *X) {
  extern char *xmatfile;
  extern int Force;
  float Xcondition;
  MATRIX *Xnorm;

  if (X->rows <= X->cols) {
    printf("ERROR: Design Matrix: nrows (%d) <= ncols (%d)\n",
           X->rows,X->cols);
    exit(1);
  }

  Xnorm = MatrixNormalizeCol(X,NULL,NULL);
  Xcondition = sqrt(MatrixNSConditionNumber(Xnorm));
  MatrixFree(&Xnorm);
  printf("INFO: Normalized Design Matrix Condition Number is %g\n",
         Xcondition);
  if (xmatfile != NULL && debug) {
    printf("INFO: Writing mat file to %s\n",xmatfile);
    if (MatlabWrite(X,xmatfile,"X")) {
      printf("ERROR: Writing mat file to %s\n",xmatfile);
      exit(1);
    }
  }
  if (Xcondition > 100000 && !Force) {
    printf("ERROR: Design matrix is badly conditioned, check for linear\n"
           "dependency  between columns (ie, two or more columns \n"
           "that add up to another column).\n\n");
    exit(1);
  }

  return(0);
}
/*---------------------------------------------------*/
static char *getstem(char *filename) {
  int filetype;
  char *stem;
  int len;

  filetype = mri_identify(filename);
  if (filetype == MRI_VOLUME_TYPE_UNKNOWN) {
    printf("ERROR: cannot determine type of %s\n",filename);
    exit(1);
  }

  len = strlen(filename);
  stem = (char *) calloc(sizeof(char),len+1);

  switch (filetype) {
  case BFLOAT_FILE:
    memmove(stem,filename,len-11);
    break;
  case MRI_MGH_FILE:
  case MRI_ANALYZE_FILE:
  case MRI_ANALYZE4D_FILE:
  case NIFTI1_FILE:
  case NII_FILE:
    memmove(stem,filename,len-4);
    break;
  default:
    printf("ERROR: cannot determine stem for %s\n",filename);
    exit(1);
    break;
  }

  return(stem);
}
/*-------------------------------------------------------------------*/
static int MatrixWriteFmt(MATRIX *M, char *fname, char *fmt) {
  int err = 0, r, c;
  FILE *fp;

  if (!strcmp(fmt,"matlab4")) {
    err = MatlabWrite(M,fname,"X");
    return(err);
  }

  fp = fopen(fname,"w");
  if (fp == NULL) {
    printf("ERROR: could not open %s for writing\n",fname);
    return(1);
  }

  for (r=1; r <= M->rows; r++) {
    for (c=1; c <= M->cols; c++)
      fprintf(fp,"%g ",M->rptr[r][c]);
    fprintf(fp,"\n");
  }
  fclose(fp);
  return(0);
}
