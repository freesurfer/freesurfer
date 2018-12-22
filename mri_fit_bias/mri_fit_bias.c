/**
 * @file  mri_fit_bias.c
 * @brief Computes spatial intensity bias field by fitting to a parametric model
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/04/30 21:32:21 $
 *    $Revision: 1.3 $
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

/*!
\file mri_fit_bias.c
\brief Computes spatial intensity bias field by fitting to a parametric model
\author Douglas Greve
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "mri2.h"
#include "macros.h"
#include "diag.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "matrix.h"
#ifdef _OPENMP
#include "romp_support.h"
#endif

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;
MRI *MRIpolyfitBiasField(MRI *vol2, int order2, MRI *seg2, MRI *bias2);

static char vcid[] = "$Id: mri_fit_bias.c,v 1.3 2012/04/30 21:32:21 greve Exp $";
const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *subject, *SUBJECTS_DIR;
char *srcfile=NULL, *trgfile=NULL;
char *biasfile=NULL, *yhatfile=NULL, *resfile=NULL;
char *maskfile = "brainmask.mgz";
char *segfile = "aseg.mgz";
MRI *src, *seg, *segerode, *mask, *trg, *wmmask;
int niters = 10, polyorder = 3;
int nthreads=0;

MATRIX *MRIbiasPolyReg(int order, MRI *mask);
MATRIX *MRIbiasXsegs(MRI *seg);
MRI *MRImergeSegs(MRI *seg, int *seglist, int nsegs, int NewSegId, MRI *newseg);
MRI *MRImatchSegs(MRI *seg, int *seglist, int nsegs, int NewSegId, MRI *newseg);
MRI *MRIdctField(MRI *template, int nP[3]);
MRI *MRIfitDCTField(MRI *field, MRI *mask, int nP[3]);
MATRIX *MRIbiasDCTReg(MRI *dctfield, MRI *mask, MATRIX *X);
/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,c,r,s,k,nmask,err;
  int wmsegs[] = {2,41,251,252,253,254,255},nwmsegs; //7,46
  int csfsegs[] = {4,5,14,24,43,44,75,31,63,15},ncsfsegs;
  int ctxsegs[] = {3,42},nctxsegs;
  int exsegs[] = {30,62,77,85,16,7,8,46,47,12,51,13,52,11,50,17,53,18,54,10,49,28,60,26,58},nexsegs;
  char tmpstr[2000];
  double sumval,val, rvarSeg, rvarBias, mres, ywmmean, ZtZcond;
  MATRIX *Zwm, *Zm, *ywm, *f, *alpha0, *alpha, *X;
  MATRIX *Xt,*XtX,*iXtX,*iXtXXt, *ZtZ;
  MATRIX *phat, *y0, *yhat,  *beta, *res, *fhat;
  MRI *mritmp, *aseg, *dct;
  int nP[3];

  nwmsegs = sizeof(wmsegs)/sizeof(int);
  ncsfsegs = sizeof(csfsegs)/sizeof(int);
  nctxsegs = sizeof(ctxsegs)/sizeof(int);
  nexsegs = sizeof(exsegs)/sizeof(int);

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);
  printf("%s:%d\n",__FILE__,__LINE__);

  nP[0] = polyorder;
  nP[1] = polyorder;
  nP[2] = polyorder;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,srcfile);
  printf("Reading %s\n",tmpstr);
  src = MRIread(tmpstr);
  if(src==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,segfile);
  printf("Reading %s\n",tmpstr);
  seg = MRIread(tmpstr);
  if(seg==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,maskfile);
  printf("Reading %s\n",tmpstr);
  mask = MRIread(tmpstr);
  if(mask==NULL) exit(1);

  printf("Creating DCT field\n");fflush(stdout);
  dct = MRIdctField(src, nP);
  printf(" ... done creating DCT field\n");fflush(stdout);

  // keep a copy of the aseg to use later
  aseg = MRIcopy(seg,NULL); 

  // Create mask from wm in aseg
  wmmask = MRImatchSegs(seg,wmsegs, nwmsegs, 1, NULL);
  wmmask = MRIerodeSegmentation(wmmask, NULL,2,0);
  MRIwrite(wmmask,"wmmask.mgh");

  // Modify aseg to merge and exclude some segs
  seg = MRImergeSegs(seg, wmsegs, nwmsegs, 2, seg);
  seg = MRImergeSegs(seg, csfsegs, ncsfsegs, 0, seg);
  seg = MRImergeSegs(seg, ctxsegs, nctxsegs, 3, seg);
  seg = MRImergeSegs(seg, exsegs, nexsegs, 0, seg);
  seg = MRIerodeSegmentation(seg, NULL,1,0);
  MRIwrite(seg,"newseg.mgz");

  // Create initial estimate of bias only from WM
  //Zwm = MRIbiasPolyReg(polyorder, wmmask);
  Zwm = MRIbiasDCTReg(dct, wmmask, NULL);

  // Extract the WM signal from the input
  ywm = MRIvol2mat(src, wmmask, 1, NULL);
  // Compute the mean
  ywmmean = 0;
  for(r=0; r < ywm->rows; r++) ywmmean += ywm->rptr[r+1][1];
  ywmmean /= ywm->rows;
  printf("ywmmean = %g\n",ywmmean);
  // Rescale so that the mean is 1
  MatrixScalarMul(ywm, 1/ywmmean, ywm);

  alpha0 = MatrixGlmFit(ywm, Zwm, &rvarBias, NULL);
  //for(r=0; r < alpha0->rows; r++) alpha0->rptr[r+1][1] = 
  //  alpha0->rptr[r+1][1]/alpha0->rptr[1][1];
  MatrixWriteTxt("alpha0.txt",alpha0);
  MatrixWriteTxt("Zwm.txt",Zwm);
  MatrixWriteTxt("ywm.txt",ywm);
  printf("wrote alpha0\n"); fflush(stdout);
  MatrixPrint(stdout,alpha0);

  // Extend init bias to all segs
  //Zm = MRIbiasPolyReg(polyorder,seg);
  Zm = MRIbiasDCTReg(dct, seg, NULL);
  phat = MatrixMultiply(Zm,alpha0,NULL);
  ZtZ  = MatrixMtM(Zm, NULL);
  ZtZcond = MatrixConditionNumber(ZtZ);
  printf("ZtZcond = %g\n",ZtZcond);

  X    = MRIbiasXsegs(seg);
  Xt   = MatrixTranspose(X,NULL);
  //XtX  = MatrixMtM(X, NULL);
  XtX  = MatrixMultiply(Xt, X, NULL);
  iXtX = MatrixInverse(XtX,NULL);
  iXtXXt = MatrixMultiplyD(iXtX,Xt,NULL);

  // Extract the signal from all segs
  f = MRIvol2mat(src,seg,1,NULL);
  // Rescale so that the mean in WM is 1
  for(r=0; r < f->rows; r++) f->rptr[r+1][1] /= ywmmean;

  y0 = NULL;
  yhat = NULL;
  beta = NULL;
  alpha = NULL;
  rvarSeg = 0;
  for(k=0; k < niters; k++){
    y0 = MatrixDivideElts(f,phat,y0); // remove bias
    beta = MatrixMultiplyD(iXtXXt,y0,beta); // compute seg means
    yhat = MatrixMultiplyD(X,beta,yhat); // estimated input from seg means
    res  = MatrixSubtract(y0, yhat, NULL);
    rvarSeg = VectorVar(res,&mres);
    rvarSeg = rvarSeg*(X->rows-1)/(X->rows-X->cols);

    phat = MatrixDivideElts(f,yhat,phat); // estimate of bias field as input/yhat
    alpha = MatrixGlmFit(phat, Zm, &rvarBias, alpha); // fit bias field
    phat = MatrixMultiplyD(Zm,alpha,phat); // new estimate of bias field
    printf("%2d %10.8f %10.8f\n",k,rvarSeg,rvarBias);
  }
  MatrixWriteTxt("alpha.txt",alpha);
  MatrixPrint(stdout,alpha);

  //MatrixWrite(Xsegs,"Xsegs.mat","Xsegs");
  //MatrixWrite(K,"Z.mat","Z");

  MatrixFree(&alpha0);
  MatrixFree(&y0);
  MatrixFree(&res);
  MatrixFree(&f);
  MatrixFree(&phat);
  MatrixFree(&yhat);
  MatrixFree(&Zm);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&iXtXXt);

  // Apply to all voxels in the brain
  printf("Apply to all voxels in the brain\n");
  //Zm = MRIbiasPolyReg(polyorder,mask);
  Zm = MRIbiasDCTReg(dct, mask, NULL);
  phat = MatrixMultiplyD(Zm,alpha,NULL);
  f = MRIvol2mat(src,mask,1,NULL);

  printf("Computing output\n");
  y0 = MatrixDivideElts(f,phat,NULL);
  trg = MRImat2vol(y0, mask, 1, NULL);
  MRIcopyPulseParameters(src, trg);

  // Normalize mri so that mean in mask is 110
  printf("Rescaling\n");
  sumval = 0;
  nmask = 0;
  for(s=0; s < wmmask->depth; s++){
    for(c=0; c < wmmask->width; c++){
      for(r=0; r < wmmask->height; r++){
	if(MRIgetVoxVal(wmmask,c,r,s,0) < .0001) continue;
	nmask ++;
	sumval += MRIgetVoxVal(trg,c,r,s,0);
      }
    }
  }
  sumval = sumval/nmask; // now it's the mean
  printf("In-mask mean %g  %d\n",sumval,nmask);
  // Rescale all values
  for(s=0; s < wmmask->depth; s++){
    for(c=0; c < wmmask->width; c++){
      for(r=0; r < wmmask->height; r++){
	if(MRIgetVoxVal(mask,c,r,s,0))	val = MRIgetVoxVal(trg,c,r,s,0);
	else val = 0;
	val = 110*val/sumval;
	if(val <   0) val = 0;
	MRIsetVoxVal(trg,c,r,s,0,val);
      }
    }
  }
  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,trgfile);
  printf("Wrting to %s\n",tmpstr);
  err = MRIwrite(trg,tmpstr);
  if(err) exit(1);

  if(biasfile){
    mritmp = MRImat2vol(phat, mask, 1, NULL);
    MRIcopyPulseParameters(src, mritmp);
    printf("Wrting to %s\n",biasfile);
    err = MRIwrite(mritmp,biasfile);
    if(err) exit(1);
    MRIfree(&mritmp);
  }
  MatrixFree(&phat);

  if(yhatfile || resfile){
    //Zm = MRIbiasPolyReg(polyorder,aseg);
    Zm = MRIbiasDCTReg(dct, aseg, NULL);
    phat = MatrixMultiply(Zm,alpha,NULL);
    f = MRIvol2mat(src,aseg,1,NULL);
    y0 = MatrixDivideElts(f,phat,NULL);
    X = MRIbiasXsegs(aseg);
    beta = MatrixGlmFit(y0, X, &rvarSeg, NULL);
    yhat = MatrixMultiplyD(X,beta,NULL);
    if(yhatfile){
      mritmp = MRImat2vol(yhat, aseg, 1, NULL);
      MRIcopyPulseParameters(src, mritmp);
      printf("Wrting to %s\n",yhatfile);
      err = MRIwrite(mritmp,yhatfile);
      if(err) exit(1);
      MRIfree(&mritmp);
    }
    if(resfile){
      fhat = MatrixMultiplyElts(yhat,phat,NULL);
      res = MatrixSubtract(f,fhat,NULL);
      mritmp = MRImat2vol(res, aseg, 1, NULL);
      MRIcopyPulseParameters(src, mritmp);
      printf("Wrting to %s\n",resfile);
      err = MRIwrite(mritmp,resfile);
      if(err) exit(1);
      MRIfree(&mritmp);
    }
    MatrixFree(&beta);
    MatrixFree(&yhat);
    MatrixFree(&res);
  }

  printf("mri_fit_bias done\n");
  exit(0);
}
/* ------------------------------------------------------*/
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--s")) {
      if (nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--src")) {
      if (nargc < 1) CMDargNErr(option,1);
      srcfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--trg")) {
      if (nargc < 1) CMDargNErr(option,1);
      trgfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--bias")) {
      if (nargc < 1) CMDargNErr(option,1);
      biasfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--yhat")) {
      if (nargc < 1) CMDargNErr(option,1);
      yhatfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--res")) {
      if (nargc < 1) CMDargNErr(option,1);
      resfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--seg")) {
      if (nargc < 1) CMDargNErr(option,1);
      segfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--niters")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&niters);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--order")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&polyorder);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* -----------------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -----------------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --s subject : to get subject's dir\n");
  printf("   --src srcfile : input to intensity normalize\n");
  printf("   --trg trgfile : intensity normalized output\n");
  printf("\n");
  printf("   --niters N\n");
  printf("   --order polyorder\n");
  printf("   --bias biasfile\n");
  printf("   --yhat yhatfile\n");
  printf("   --res resfile\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* -----------------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* -----------------------------------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* -----------------------------------------------------------------------*/
static void check_options(void) {
  if(subject == NULL){
    printf("ERROR: need subject \n");
    exit(1);
  }
  if(srcfile == NULL){
    printf("ERROR: need source file\n");
    exit(1);
  }
  if(trgfile == NULL){
    printf("ERROR: need source file\n");
    exit(1);
  }

  return;
}
/* -----------------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"subject  %s\n",subject);

  return;
}

/*----------------------------------------------------*/
MRI *MRIpolyfitBiasField(MRI *vol, int order, MRI *seg, MRI *bias)
{
  int c,r,s, nseg, nX, *segids, nsegids,z,k,segid,*tmplist,nthseg, npoly=3;
  MATRIX *X, *y, *Xt, *XtX, *iXtX, *Xty, *beta;
  double dc,dr,ds,val;
  double nchalf,nrhalf,nshalf;

  if(!bias) {
    bias = MRIallocSequence(vol->width, vol->height, vol->depth,MRI_FLOAT, 1);
    MRIcopyHeader(vol,bias);
  }
  if(MRIdimMismatch(vol, bias, 0)){
    printf("ERROR: MRIpolyfitBiasField(): vol/bias dim mismatch\n");
    return(NULL);
  }
  if(MRIdimMismatch(vol, seg, 0)){
    printf("ERROR: MRIpolyfitBiasField(): vol/seg dim mismatch\n");
    return(NULL);
  }

  // Count number of voxels in seg so can alloc
  nseg = 0;
  for (c=0; c < seg->width; c++)
    for (r=0; r < seg->height; r++)
      for (s=0; s < seg->depth; s++)
	if(MRIgetVoxVal(seg,c,r,s,0) > 0.5) nseg++;
  printf("MRIpolyfitBiasField(): found %d voxels in seg\n",nseg);

  // Get number of unique list segmentation IDs
  segids = MRIsegmentationList(seg, &nsegids);
  // Check whether there is a segmentation 0
  z = 0;
  for(k=0; k<nsegids;k++) if(segids[k] == 0) z = 1;
  if(z){
    tmplist = (int *) calloc(nsegids-1,sizeof(int));
    nthseg = 0;
    for(k=0; k<nsegids;k++) {
      if(segids[k] == 0) continue;
      tmplist[nthseg] = segids[k];
      nthseg ++;
    }
    free(segids);
    segids = tmplist;
    nsegids = nthseg;
  }
  printf("Found %d non-zero segids\n",nsegids);
  //for(k=0; k<nsegids;k++) printf("%2d %5d\n",k,segids[k]);
  
  // Alloc
  if(order == 2) npoly = 9;
  if(order == 3) npoly = 9+10;
  nX = nsegids + npoly; 
  X = MatrixAlloc(nseg,nX,MATRIX_REAL);
  y = MatrixAlloc(nseg,1,MATRIX_REAL);

  //Scale CRS to make X better conditioned
  nchalf = seg->width/2.0;
  nrhalf = seg->height/2.0;
  nshalf = seg->depth/2.0;

  // Set up the matrices to do the estimation
  nseg = 0;
  for(s=0; s < seg->depth; s++){
    ds = (s - nshalf)/nshalf; 
    for(c=0; c < seg->width; c++){
      dc = (c - nchalf)/nchalf;
      for(r=0; r < seg->height; r++){
	dr = (r - nrhalf)/nrhalf;
	segid = MRIgetVoxVal(seg,c,r,s,0);
	if(segid == 0) continue;
	for(k=0; k<nsegids;k++) if(segid == segids[k]) break;
	// Linear
	X->rptr[nseg+1][1] = dc;
	X->rptr[nseg+1][2] = dr;
	X->rptr[nseg+1][3] = ds;
	// Pure Quadratic
	X->rptr[nseg+1][4] = dc*dc;
	X->rptr[nseg+1][5] = dr*dr;
	X->rptr[nseg+1][6] = ds*ds;
	// Quadratic cross-products
	X->rptr[nseg+1][7] = dc*dr;
	X->rptr[nseg+1][8] = dc*ds;
	X->rptr[nseg+1][9] = dr*ds;
	if(order > 2){
	  // Cubic
	  X->rptr[nseg+1][10] = dc*dc*dc;
	  X->rptr[nseg+1][11] = dr*dr*dr;
	  X->rptr[nseg+1][12] = ds*ds*ds;
	  X->rptr[nseg+1][13] = dc*dc*dr;
	  X->rptr[nseg+1][14] = dc*dc*ds;
	  X->rptr[nseg+1][15] = dc*dr*ds;
	  X->rptr[nseg+1][16] = dc*dr*dr;
	  X->rptr[nseg+1][17] = dc*ds*ds;
	  X->rptr[nseg+1][18] = dr*dr*ds;
	  X->rptr[nseg+1][19] = dr*ds*ds;
	}
	// Constant
	X->rptr[nseg+1][npoly+k+1] = 1.0;
	// Input data
	y->rptr[nseg+1][1] = MRIgetVoxVal(vol,c,r,s,0);
	nseg++;
      }
    }
  }

  // Do the estimation
  Xt   = MatrixTranspose(X,NULL);
  XtX  = MatrixMtM(X, NULL);
  iXtX = MatrixInverse(XtX,NULL);
  Xty  = MatrixMultiplyD(Xt,y,NULL);
  beta = MatrixMultiplyD(iXtX,Xty,NULL);

  //MatrixWrite(X,"X.mat","X");
  //MatrixWrite(y,"y.mat","y");
  //MatrixWrite(beta,"beta.mat","beta");
  //MatrixWrite(Xt,"Xt.mat","Xt");
  //MatrixWrite(XtX,"XtX.mat","XtX");
  //MatrixWrite(iXtX,"iXtX.mat","iXtX");

  MatrixFree(&y);
  MatrixFree(&X);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&Xty);

  // Use the estimation to compute bias field at all voxels
  X = MatrixAlloc(1,nX,MATRIX_REAL);
  for(s=0; s < seg->depth; s++){
    ds = (s - nshalf)/nshalf; //norm makes X better conditioned
    for(c=0; c < seg->width; c++){
      dc = (c - nchalf)/nchalf;
      for(r=0; r < seg->height; r++){
	dr = (r - nrhalf)/nrhalf;
	X->rptr[1][1] = dc;
	X->rptr[1][2] = dr;
	X->rptr[1][3] = ds;
	X->rptr[1][4] = dc*dc;
	X->rptr[1][5] = dr*dr;
	X->rptr[1][6] = ds*ds;
	X->rptr[1][7] = dc*dr;
	X->rptr[1][8] = dc*ds;
	X->rptr[1][9] = dr*ds;
	if(order > 2){
	  // Cubic
	  X->rptr[1][10] = dc*dc*dc;
	  X->rptr[1][11] = dr*dr*dr;
	  X->rptr[1][12] = ds*ds*ds;
	  X->rptr[1][13] = dc*dc*dr;
	  X->rptr[1][14] = dc*dc*ds;
	  X->rptr[1][15] = dc*dr*ds;
	  X->rptr[1][16] = dc*dr*dr;
	  X->rptr[1][17] = dc*ds*ds;
	  X->rptr[1][18] = dr*dr*ds;
	  X->rptr[1][19] = dr*ds*ds;
	}
	// Note: all constant terms excluded
	y = MatrixMultiply(X,beta,y);
	val = MRIgetVoxVal(vol,c,r,s,0) - y->rptr[1][1];
	MRIsetVoxVal(bias,c,r,s,0,val);
      }
    }
  }
  MatrixFree(&y);
  MatrixFree(&X);
  free(segids);
  return(bias);

}

/*----------------------------------------------------*/
MATRIX *MRIbiasXsegs(MRI *seg)
{
  int c,r,s,z,k,nseg, *segids, nsegids,segid,*tmplist,nthseg;
  MATRIX *X;

  // Count number of voxels in seg so can alloc
  nseg = 0;
  for (c=0; c < seg->width; c++)
    for (r=0; r < seg->height; r++)
      for (s=0; s < seg->depth; s++)
	if(MRIgetVoxVal(seg,c,r,s,0) > 0.5) nseg++;
  printf("MRIbiasXsegs(): found %d voxels in seg\n",nseg);

  // Get number of unique list segmentation IDs
  segids = MRIsegmentationList(seg, &nsegids);
  // Check whether there is a segmentation 0
  z = 0;
  for(k=0; k<nsegids;k++) if(segids[k] == 0) z = 1;
  if(z){
    tmplist = (int *) calloc(nsegids-1,sizeof(int));
    nthseg = 0;
    for(k=0; k<nsegids;k++) {
      if(segids[k] == 0) continue;
      tmplist[nthseg] = segids[k];
      nthseg ++;
    }
    free(segids);
    segids = tmplist;
    nsegids = nthseg;
  }
  printf("Found %d non-zero segids\n",nsegids);
  X = MatrixAlloc(nseg,nsegids,MATRIX_REAL);

  // Use slice, col, row to be consistent with matlab
  nseg = 0;
  for (s=0; s < seg->depth; s++){
    for (c=0; c < seg->width; c++){
      for (r=0; r < seg->height; r++){
	segid = MRIgetVoxVal(seg,c,r,s,0);
	if(segid == 0) continue;
	for(k=0; k<nsegids;k++) if(segid == segids[k]) break;
	X->rptr[nseg+1][k+1] = 1;
	nseg++;
      }
    }
  }
  return(X);
}
/*----------------------------------------------------*/
MATRIX *MRIbiasPolyReg(int order, MRI *mask)
{
  int c,r,s, nmask, nX, npoly;
  double dc,dr,ds,nchalf,nrhalf,nshalf;
  MATRIX *X;

  // Count number of voxels in mask so can alloc
  nmask = 0;
  for (c=0; c < mask->width; c++)
    for (r=0; r < mask->height; r++)
      for (s=0; s < mask->depth; s++)
	if(MRIgetVoxVal(mask,c,r,s,0) > 0.0001) nmask++;
  printf("MRIbiasPolyReg(): found %d voxels in mask\n",nmask);

  // Alloc
  npoly=4;
  if(order == 2) npoly = 9;
  if(order == 3) npoly = 9+10;
  nX = npoly+1; 
  X = MatrixAlloc(nmask,nX,MATRIX_REAL);

  //Scale CRS to make X better conditioned
  nchalf = mask->width/2.0;
  nrhalf = mask->height/2.0;
  nshalf = mask->depth/2.0;

  // Set up the matrix to do the estimation
  // Use slice, col, row to be consistent with matlab
  nmask = 0;
  for(s=0; s < mask->depth; s++){
    ds = (s - nshalf)/nshalf; 
    for(c=0; c < mask->width; c++){
      dc = (c - nchalf)/nchalf;
      for(r=0; r < mask->height; r++){
	dr = (r - nrhalf)/nrhalf;
	if(MRIgetVoxVal(mask,c,r,s,0) < 0.0001) continue;
	// Linear
	X->rptr[nmask+1][1] = dc;
	X->rptr[nmask+1][2] = dr;
	X->rptr[nmask+1][3] = ds;
	if(order > 1){
	  // Pure Quadratic
	  X->rptr[nmask+1][4] = dc*dc;
	  X->rptr[nmask+1][5] = dr*dr;
	  X->rptr[nmask+1][6] = ds*ds;
	  // Quadratic cross-products
	  X->rptr[nmask+1][7] = dc*dr;
	  X->rptr[nmask+1][8] = dc*ds;
	  X->rptr[nmask+1][9] = dr*ds;
	}
	if(order > 2){
	  // Cubic
	  X->rptr[nmask+1][10] = dc*dc*dc;
	  X->rptr[nmask+1][11] = dr*dr*dr;
	  X->rptr[nmask+1][12] = ds*ds*ds;
	  X->rptr[nmask+1][13] = dc*dc*dr;
	  X->rptr[nmask+1][14] = dc*dc*ds;
	  X->rptr[nmask+1][15] = dc*dr*ds;
	  X->rptr[nmask+1][16] = dc*dr*dr;
	  X->rptr[nmask+1][17] = dc*ds*ds;
	  X->rptr[nmask+1][18] = dr*dr*ds;
	  X->rptr[nmask+1][19] = dr*ds*ds;
	}
	// Constant
	X->rptr[nmask+1][nX] = 1.0;
	nmask++;
      }
    }
  }

  return(X);
}

/*!
  \fn MRI *MRImergeSegs(MRI *seg, int *seglist, int nsegs, int NewSegId, MRI *newseg)
  \brief Merges multiple segmentations into one. Can be done in-place.
  \parameter seg - original segmentation
  \parameter seglist - list of segmentation IDs to merge
  \parameter nsegs - length of list
  \parameter NewSegId - replace values in list with NewSegId
  \parameter newseg - new segmentation (also passed as output)
*/
MRI *MRImergeSegs(MRI *seg, int *seglist, int nsegs, int NewSegId, MRI *newseg)
{
  int c,r,s,n,segid;

  if(newseg == NULL) newseg = MRIcopy(seg,NULL);

  for (c=0; c < seg->width; c++){
    for (r=0; r < seg->height; r++){
      for (s=0; s < seg->depth; s++){
	segid = MRIgetVoxVal(seg,c,r,s,0);
	MRIsetVoxVal(newseg,c,r,s,0, segid);
	for(n=0; n < nsegs; n++){
	  if(segid == seglist[n]){
	    MRIsetVoxVal(newseg,c,r,s,0,NewSegId);
	    break;
	  }
	}
      }
    }
  }
  return(newseg);
}

/*!
  \fn MRI *MRImatchSegs(MRI *seg, int *seglist, int nsegs, int MaskId, MRI *mask)
  \brief Creates a binary mask of voxels that match any of the IDs in the
    segmentations list. Can be done in-place.
  \parameter seg - original segmentation
  \parameter seglist - list of segmentation IDs to merge
  \parameter nsegs - length of list
  \parameter MaskId - replace values in list with MaskId
  \parameter mask - new segmentation (also passed as output)
*/
MRI *MRImatchSegs(MRI *seg, int *seglist, int nsegs, int MaskId, MRI *mask)
{
  int c,r,s,n,segid;

  if(mask == NULL) mask = MRIcopy(seg,NULL);

  for (c=0; c < seg->width; c++){
    for (r=0; r < seg->height; r++){
      for (s=0; s < seg->depth; s++){
	MRIsetVoxVal(mask,c,r,s,0, 0);
	segid = MRIgetVoxVal(seg,c,r,s,0);
	for(n=0; n < nsegs; n++){
	  if(segid == seglist[n]){
	    MRIsetVoxVal(mask,c,r,s,0,MaskId);
	    break;
	  }
	}
      }
    }
  }
  return(mask);
}



MRI *MRIfitDCTField(MRI *field, MRI *mask, int nP[3])
{
  int r;
  MRI *P, *fit;
  MATRIX *X, *y, *XtX, *iXtX, *Xty, *beta, *yhat;
  double XtXcond,yhatmn;

  P = MRIdctField(field, nP);
  //MRIwrite(P,"P.mgh");

  y = MRIvol2mat(field, mask, 1, NULL);
  if(y==NULL) return(NULL);
  X = MRIvol2mat(P, mask, 1, NULL);
  if(X==NULL) return(NULL);
  printf("nvox = %d\n",y->rows);

  XtX = MatrixMtM(X, NULL);
  iXtX = MatrixInverse(XtX, NULL);
  XtXcond = MatrixConditionNumber(XtX);
  printf("MRIfitDCTField(): XtXcond=%g\n", XtXcond);
  if(iXtX == NULL) {
    XtXcond = MatrixConditionNumber(XtX);
    printf("ERROR: MRIfitDCTField(): matrix cannot be inverted, cond=%g\n", XtXcond);
    MatrixWriteTxt("X.mtx",X);
    return(NULL);
  }
  Xty = MatrixAtB(X, y, NULL);
  beta = MatrixMultiplyD(iXtX, Xty, NULL);
  printf("MRIfitDCTField(): beta\n");
  MatrixPrint(stdout,beta); fflush(stdout);
  yhat = MatrixMultiplyD(X, beta, NULL);

  yhatmn = 0;
  for(r=0; r < yhat->rows; r++) yhatmn += yhat->rptr[r+1][1];
  yhatmn = yhatmn/yhat->rows;
  MatrixScalarMul( yhat, 1/yhatmn, yhat);

  fit = MRImat2vol(yhat, mask, 1, NULL);

  MatrixFree(&X);
  MatrixFree(&y);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&Xty);
  MatrixFree(&beta);
  MatrixFree(&yhat);

  return(fit);
} 

MATRIX *MRIbiasDCTReg(MRI *dctfield, MRI *mask, MATRIX *X)
{
  X = MRIvol2mat(dctfield, mask, 1, X);
  if(X==NULL) return(NULL);
  return(X);
} 


MRI *MRIdctField(MRI *template, int nP[3])
{
  int c, r, s, nPtot,k=0, n;
  double *cdm, *rdm, *sdm;
  MRI *P;

  // dm are vectors that go from -0.5 to +0.5
  cdm = (double *) calloc(template->width,sizeof(double));
  for(k=0;k<template->width;k++) cdm[k] = (k - template->width/2.0)/template->width;
  rdm = (double *) calloc(template->height,sizeof(double));
  for(k=0;k<template->height;k++) rdm[k] = (k - template->height/2.0)/template->height;
  sdm = (double *) calloc(template->depth,sizeof(double));
  for(k=0;k<template->depth;k++) sdm[k] = (k - template->depth/2.0)/template->depth;

  nPtot = nP[0]+nP[1]+nP[2] + 1;
  P = MRIallocSequence(template->width, template->height, template->depth, MRI_FLOAT, nPtot);

  // Start with a constant = 1
  k = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for(c=0; c < template->width; c++){
    int r, s;
    for(r=0; r < template->height; r++){
      for(s=0; s < template->depth; s++){
	MRIsetVoxVal(P,c,r,s,k,1);
      }
    }
  }
  k = k + 1;
  // Bases in the col direction
  for(n=1; n <= nP[0]; n++){
    #ifdef HAVE_OPENMP
    #pragma omp parallel for 
    #endif
    for(c=0; c < template->width; c++){
      int r, s;
      double f;
      f = cos(2*M_PI*cdm[c]*n);
      for(r=0; r < template->height; r++){
	for(s=0; s < template->depth; s++){
	  MRIsetVoxVal(P,c,r,s,k,f);
	}
      }
    }
    k = k + 1;
  }
  // Bases in the row direction
  for(n=1; n <= nP[0]; n++){
    #ifdef HAVE_OPENMP
    #pragma omp parallel for 
    #endif
    for(r=0; r < template->height; r++){
      int c, s;
      double f;
      f = cos(2*M_PI*rdm[r]*n);
      for(c=0; c < template->width; c++){
	for(s=0; s < template->depth; s++){
	  MRIsetVoxVal(P,c,r,s,k,f);
	}
      }
    }
    k = k + 1;
  }
  // Bases in the slice direction
  for(n=1; n <= nP[0]; n++){
    #ifdef HAVE_OPENMP
    #pragma omp parallel for 
    #endif
    for(s=0; s < template->depth; s++){
      int c, r;
      double f;
      f = cos(2*M_PI*sdm[s]*n);
      for(c=0; c < template->width; c++){
	for(r=0; r < template->height; r++){
	  MRIsetVoxVal(P,c,r,s,k,f);
	}
      }
    }
    k = k + 1;
  }
  
  free(cdm);
  free(rdm);
  free(sdm);
  return(P);
}
