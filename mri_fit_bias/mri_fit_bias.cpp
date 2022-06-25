/**
 * @brief Computes spatial intensity bias field by fitting to a parametric model
 *
 */
/*
 * Original Author: Douglas N. Greve
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
#include "cma.h"
#include "region.h"
#ifdef _OPENMP
#include "romp_support.h"
#endif

typedef struct {
  MATRIX *offsetRAS; // 0 point in RAS scanner space
  MATRIX *offsetCRS; // 0 point in col, row, slice
  int nfields;
  double *period; // one for each field
  int    *dim;    // one for each field
  MRI *mri_template;
  MRI *mask;
  MRI *field;
  double lpcutoff; // cutoff period in mm
  MRI_REGION *bbfit; // bounding box
} DCT_FIELD_SPEC, DCTFS;

typedef struct {
  MRI *ubermask;
  MRI *inputbc;
  MRI *input;
  MRI *biasfield;
  MRI *bbinput;
  MRI *srcseg;
  MRI *seg;
  MRI *segerode;
  int nerode;
  MRI_REGION *bb;
  MRI *bbseg;
  int logflag;
  int nP[3];
  MATRIX *y, *Xseg, *Xbf;
  MATRIX *beta, *betaseg, *betabf;
  MATRIX *X, *Xt,*XtX,*iXtX,*Xty;
  double XtXcond;
  int *exsegs, nexsegs;
  int *wmsegs, nwmsegs;
  int *ctxsegs, nctxsegs;
  DCT_FIELD_SPEC *dct;
  double wmmean;
  double thresh = 0;
} MRI_BIAS_FIELD, MRIBF;


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;
MRI *MRIpolyfitBiasField(MRI *vol2, int order2, MRI *seg2, MRI *bias2);

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *srcfile=NULL, *outfile=NULL,*biasfile=NULL;
char *maskfile = NULL, *segfile = NULL;
MRI *src, *seg, *mask, *out;
int nthreads=0;
double lpfcutoffmm = 23;
double thresh = 0;
char *Xfile=NULL, *yfile=NULL, *dctfile=NULL,*betafile=NULL;
char *yhatfile=NULL, *resfile=NULL;
char *SUBJECTS_DIR;

MATRIX *MRIbiasPolyReg(int order, MRI *mask);
MATRIX *MRIbiasXsegs(MRI *seg);
MRI *MRIdctField(MRI *mri_template, int nP[3], int offset);
MRI *MRIfitDCTField(MRI *field, MRI *mask, int nP[3]);
MATRIX *MRIbiasDCTReg(MRI *dctfield, MRI *mask, MATRIX *X);
int MRIbiasFieldCorLog(MRIBF *bf);
int DCTfield(DCTFS *dct);
int DCTspecLPF(DCTFS *dct);
DCTFS *DCTalloc(void);
MRI *MRIapplyBiasField(MRI *input, MRI *bf, MRI *seg, MRI *mask, double targwmval, double thresh, MRI *bc);

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs;
  int wmsegs[] = {2,41,251,252,253,254,255},nwmsegs; //7,46
  int csfsegs[] = {4,5,14,24,43,44,75,31,63,15},ncsfsegs;
  int ctxsegs[] = {3,42},nctxsegs;
  int exsegs[] = {30,62,77,85,16,7,8,46,47,12,51,13,52,11,50,17,53,18,54,10,49,28,60,26,58},nexsegs;

  nwmsegs = sizeof(wmsegs)/sizeof(int);
  ncsfsegs = sizeof(csfsegs)/sizeof(int);
  nctxsegs = sizeof(ctxsegs)/sizeof(int);
  nexsegs = sizeof(exsegs)/sizeof(int);

  nargs = handleVersionOption(argc, argv, "mri_fit_bias");
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

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  printf("Reading %s\n",srcfile);
  src = MRIread(srcfile);
  if(src==NULL) exit(1);

  printf("Reading %s\n",segfile);
  seg = MRIread(segfile);
  if(seg==NULL) exit(1);

  printf("Reading %s\n",maskfile);
  mask = MRIread(maskfile);
  if(mask==NULL) exit(1);

  MRIBF *bf;
  bf = (MRIBF *)calloc(sizeof(MRIBF),1);
  bf->input = src;
  bf->srcseg = seg;
  bf->nerode = 1;
  bf->logflag = 1;
  bf->ubermask = mask;
  bf->thresh = thresh;
  bf->dct = DCTalloc();
  bf->dct->mri_template = src;
  bf->dct->lpcutoff = lpfcutoffmm;
  MRIbiasFieldCorLog(bf);

  // Apply the bias field. This also rescales so that the mean in the eroded
  // WM is 110. It also changes the values in the biasfield to account for the
  // rescaling so that bcoutput = input/biasfield, though I think the unrescaled
  // bias field is probably more informative (or scaling it so that wm = 1)
  MRI *bc=MRIapplyBiasField(bf->input, bf->biasfield, bf->segerode, bf->ubermask, 110, thresh, NULL);
  printf("Writing output to %s\n",outfile);fflush(stdout);
  MRIwrite(bc,outfile);
  MRIfree(&bc);

  if(biasfile)
    MRIwrite(bf->biasfield,biasfile);

  if(dctfile)
    MRIwrite(bf->dct->field,dctfile);

  if(Xfile) MatrixWriteTxt(Xfile,bf->X);
  if(yfile) MatrixWriteTxt(yfile,bf->y);
  if(betafile) MatrixWriteTxt(yfile,bf->beta);

  printf("#VMPC# mri_fit_bias VmPeak  %d\n",GetVmPeak());
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

    else if(!strcmp(option, "--sd") || !strcmp(option, "-SDIR")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      srcfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--bias")) {
      if (nargc < 1) CMDargNErr(option,1);
      biasfile = pargv[0];
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
    else if (!strcasecmp(option, "--X")) {
      if (nargc < 1) CMDargNErr(option,1);
      Xfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--y")) {
      if (nargc < 1) CMDargNErr(option,1);
      yfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--beta")) {
      if (nargc < 1) CMDargNErr(option,1);
      betafile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--dct")) {
      if (nargc < 1) CMDargNErr(option,1);
      dctfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--cutoff")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&lpfcutoffmm);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&thresh);
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
  printf("   --i inputvol : input to intensity normalize\n");
  printf("   --cutoff lpfcutoffmm (%f)\n",lpfcutoffmm);
  printf("   --seg segvol : segmentation to define WM and Cortex (eg, aseg.presurf.mgz)\n");
  printf("   --mask maskvol : zero everthing outside of mask (eg,brainmask.mgz)\n");
  printf("   --o outvol : bias corrected output\n");
  printf("   --bias biasfield  : output bias field\n");
  printf("   --dct dctvol : DCT fields file (debugging)\n");
  printf("   --thresh thresh : mask out anything <= thresh\n");
  printf("\n");
  printf("   --threads nthreads\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
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
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* -----------------------------------------------------------------------*/
static void check_options(void) {
  if(srcfile == NULL){
    printf("ERROR: need input file\n");
    exit(1);
  }
  if(segfile == NULL){
    printf("ERROR: need segmentation file\n");
    exit(1);
  }
  if(outfile == NULL){
    printf("ERROR: need output file\n");
    exit(1);
  }

  return;
}
/* -----------------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
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



MRI *MRIfitDCTField(MRI *field, MRI *mask, int nP[3])
{
  int r;
  MRI *P, *fit;
  MATRIX *X, *y, *XtX, *iXtX, *Xty, *beta, *yhat;
  double XtXcond,yhatmn;

  P = MRIdctField(field, nP, 1);
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


MRI *MRIdctField(MRI *mri_template, int nP[3], int offset)
{
  int c, r, s, nPtot,k=0, n;
  double *cdm, *rdm, *sdm;
  MRI *P;

  // dm are vectors that go from -0.5 to +0.5
  cdm = (double *) calloc(mri_template->width,sizeof(double));
  for(k=0;k<mri_template->width;k++) cdm[k] = (k - mri_template->width/2.0)/mri_template->width;
  rdm = (double *) calloc(mri_template->height,sizeof(double));
  for(k=0;k<mri_template->height;k++) rdm[k] = (k - mri_template->height/2.0)/mri_template->height;
  sdm = (double *) calloc(mri_template->depth,sizeof(double));
  for(k=0;k<mri_template->depth;k++) sdm[k] = (k - mri_template->depth/2.0)/mri_template->depth;

  nPtot = nP[0]+nP[1]+nP[2] + offset;
  P = MRIallocSequence(mri_template->width, mri_template->height, mri_template->depth, MRI_FLOAT, nPtot);

  // Start with a constant = 1
  k = 0;
  if(offset){
    #ifdef HAVE_OPENMP
    #pragma omp parallel for 
    #endif
    for(c=0; c < mri_template->width; c++){
      int r, s;
      for(r=0; r < mri_template->height; r++){
        for(s=0; s < mri_template->depth; s++){
    	  MRIsetVoxVal(P,c,r,s,k,1);
        }
      }
    }
    k = k + 1;
  }
  // Bases in the col direction
  for(n=1; n <= nP[0]; n++){
    #ifdef HAVE_OPENMP
    #pragma omp parallel for 
    #endif
    for(c=0; c < mri_template->width; c++){
      int r, s;
      double f;
      f = cos(2*M_PI*cdm[c]*n);
      for(r=0; r < mri_template->height; r++){
	for(s=0; s < mri_template->depth; s++){
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
    for(r=0; r < mri_template->height; r++){
      int c, s;
      double f;
      f = cos(2*M_PI*rdm[r]*n);
      for(c=0; c < mri_template->width; c++){
	for(s=0; s < mri_template->depth; s++){
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
    for(s=0; s < mri_template->depth; s++){
      int c, r;
      double f;
      f = cos(2*M_PI*sdm[s]*n);
      for(c=0; c < mri_template->width; c++){
	for(r=0; r < mri_template->height; r++){
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

int MRIbiasFieldCorLog(MRIBF *bf)
{
  int c,r,nwm;
  MRI *ones;
  MATRIX *bffit;
  double wmsum;

  printf("MRIbiasFieldCorLog() thresh = %g\n",bf->thresh);fflush(stdout);
  printf("VmPeak  %d\n",GetVmPeak());

  bf->seg = MRIcopy(bf->srcseg,NULL); 
  MRIcopyPulseParameters(bf->srcseg, bf->seg);

  // Zero out non-WM and non-cortex, and merge all WM and cortex
  wmsum = 0;
  nwm = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : wmsum, nwm)
  #endif
  for(c=0; c < bf->seg->width; c++){
    int r, s, segid;
    double val;
    for(r=0; r < bf->seg->height; r++){
      for(s=0; s < bf->seg->depth; s++){
        val = MRIgetVoxVal(bf->input,c,r,s,0);
	if(val <= bf->thresh){
          // Exclude voxels that have a value of <=0 
          MRIsetVoxVal(bf->seg,c,r,s,0,0);
          continue;
        }
	segid = MRIgetVoxVal(bf->seg,c,r,s,0);
	switch (segid) {
	case Left_Cerebral_White_Matter: 
	case Right_Cerebral_White_Matter:
	case 251:
	case 252:
	case 253:
	case 254:
	case 255:
	  MRIsetVoxVal(bf->seg,c,r,s,0,Left_Cerebral_White_Matter);
	  wmsum += MRIgetVoxVal(bf->input,c,r,s,0);
	  nwm ++;
	  break;
	case Left_Cerebral_Cortex: 
	case Right_Cerebral_Cortex:
	  MRIsetVoxVal(bf->seg,c,r,s,0,Left_Cerebral_Cortex);
	  break;
	case Left_Cerebellum_Cortex:
	case Right_Cerebellum_Cortex:
	  MRIsetVoxVal(bf->seg,c,r,s,0,Left_Cerebellum_Cortex);
	  break;
        case Left_Cerebellum_White_Matter:
	case Right_Cerebellum_White_Matter:
	  MRIsetVoxVal(bf->seg,c,r,s,0,Left_Cerebellum_White_Matter);
	  break;
	case Left_Caudate:
	case Right_Caudate:
	  MRIsetVoxVal(bf->seg,c,r,s,0,Left_Caudate);
	  break;
	case 258:
	  MRIsetVoxVal(bf->seg,c,r,s,0,0); // 258
	  break;
	case 259:
	  MRIsetVoxVal(bf->seg,c,r,s,0,0); // 259
	  break;
	default:
	  MRIsetVoxVal(bf->seg,c,r,s,0,0);
	  break;
	}
      }
    }
  }
  bf->wmmean = wmsum/nwm;
  printf("nwm = %d, wmsum = %g, wmmean = %g\n",nwm,wmsum,bf->wmmean);

  // Erode the segmentation
  bf->segerode = MRIerodeSegmentation(bf->seg, NULL, bf->nerode, 0);
  //MRIwrite(bf->segerode,"segerode.mgz");

  // Compute the bounding box, pad=1, probably not needed
  bf->bb = REGIONgetBoundingBox(bf->segerode, 1);
  printf("Bounding box: ");
  REGIONprint(stdout, bf->bb);

  // Extract the bounding box from seg and input
  bf->bbseg   = MRIextractRegion(bf->segerode, NULL, bf->bb);
  bf->bbinput = MRIextractRegion(bf->input,    NULL, bf->bb);
  //MRIwrite(bf->bbseg,"bbseg.mgh");

  // Compute log of input
  bf->y = MRIvol2mat(bf->bbinput, bf->bbseg, 1, NULL);
  for(r=1; r <= bf->y->rows; r++){
    for(c=1; c <= bf->y->cols; c++){
      if(bf->y->rptr[r][c] <= 0){
	printf("ERROR: y matrix is less than or equal to 0, v=%g (%d %d)\n",bf->y->rptr[r][c],r,c);
	return(1);
      }
      bf->y->rptr[r][c] = log(bf->y->rptr[r][c]);
    }
  }
  printf("VmPeak  %d\n",GetVmPeak());

  // Extract the design matrix to compute the segmentation means
  bf->Xseg = MRIbiasXsegs(bf->bbseg);

  // construct the DCT fields
  bf->dct->mask = bf->ubermask;
  bf->dct->bbfit = bf->bb;
  DCTspecLPF(bf->dct);
  DCTfield(bf->dct);
  // extract the DCT fields from the bounding box
  MRI *dctfield = MRIextractRegion(bf->dct->field, NULL, bf->bb);

  // construct regressors for the DCT
  bf->Xbf  = MRIbiasDCTReg(dctfield, bf->bbseg, NULL);

  // Compute the GLM. Repackage using fsglm?
  bf->X  = MatrixHorCat(bf->Xbf,bf->Xseg,NULL);
  bf->Xt = MatrixTranspose(bf->X,NULL);
  //bf->XtX  = MatrixMtM(bf->X, NULL); // something wrong with this in parallel
  bf->XtX  = MatrixMultiplyD(bf->Xt, bf->X,NULL);
  bf->XtXcond = MatrixConditionNumber(bf->XtX);
  printf("MRIbiasFieldLog(): XtXcond=%g\n", bf->XtXcond);
  bf->iXtX = MatrixInverse(bf->XtX,NULL);
  if(bf->iXtX == NULL) {
    printf("ERROR: MRIbiasFieldLog(): matrix cannot be inverted\n");
    MatrixWriteTxt("X.mtx",bf->X);
    return(1);
  }
  bf->Xty  = MatrixMultiplyD(bf->Xt,bf->y,NULL);
  bf->beta = MatrixMultiplyD(bf->iXtX,bf->Xty,NULL);
  MatrixFree(&bf->Xt);
  MatrixFree(&bf->y);

  // Extract the biasfield-related parameters
  bf->betabf = MatrixAlloc(bf->Xbf->cols,bf->beta->cols,MATRIX_REAL);
  for(r=1; r <= bf->betabf->rows; r++){
    for(c=1; c <= bf->betabf->cols; c++){
      bf->betabf->rptr[r][c] = bf->beta->rptr[r][c];
    }
  }
  printf("alpha\n");
  MatrixPrint(stdout,bf->betabf); fflush(stdout);

  printf("Transfering to the full FoV\n"); fflush(stdout);
  printf("VmPeak  %d\n",GetVmPeak());fflush(stdout);
  ones = MRIcopy(bf->bbseg,NULL);
  MRIconst(ones->width,ones->height,ones->depth,1,1,ones);
  MatrixFree(&bf->Xbf);
  bf->Xbf = MRIbiasDCTReg(bf->dct->field, bf->input, NULL);
  bffit = MatrixMultiply(bf->Xbf,bf->betabf,NULL);
  for(r=1; r <= bffit->rows; r++){
    for(c=1; c <= bffit->cols; c++){
      bffit->rptr[r][c] = exp(bffit->rptr[r][c]);
    }
  }
  bf->biasfield = MRImat2vol(bffit, bf->input, 1, NULL);
  printf("MRIbiasFieldCorLog() done\n");fflush(stdout);
  printf("VmPeak  %d\n",GetVmPeak());fflush(stdout);

  return(0);
}


int DCTfield(DCTFS *dct)
{
  int n,c,r,s,dim;
  MRI *mri_template = dct->mri_template;
  double v,period,t=0,offset;
  MATRIX *vox2ras=NULL,*ras2vox=NULL; //*crs=NULL,*ras=NULL;

  printf("DCTfield():\n");
  vox2ras = MRIxfmCRS2XYZ(mri_template,0);
  ras2vox = MatrixInverse(vox2ras,NULL);
  //dct->offsetCRS = MatrixMultiplyD(ras2vox,dct->offsetRAS,dct->offsetCRS);
  // Put the 0 voxel at the corner of the bounding box
  dct->offsetCRS->rptr[1][1] = dct->bbfit->x;
  dct->offsetCRS->rptr[2][1] = dct->bbfit->y;
  dct->offsetCRS->rptr[3][1] = dct->bbfit->z;
  dct->offsetRAS = MatrixMultiplyD(vox2ras,dct->offsetCRS,dct->offsetRAS);
  printf("offsetCRS %7.4f %7.4f %7.4f\n",dct->offsetCRS->rptr[1][1],
	 dct->offsetCRS->rptr[2][1],dct->offsetCRS->rptr[3][1]);
  printf("offsetRAS %7.4f %7.4f %7.4f\n",dct->offsetRAS->rptr[1][1],
	 dct->offsetRAS->rptr[2][1],dct->offsetRAS->rptr[3][1]);

  if(dct->field) MRIfree(&dct->field);
  dct->field = MRIallocSequence(mri_template->width, mri_template->height, 
				mri_template->depth, MRI_FLOAT, dct->nfields);
  MRIcopyHeader(mri_template,dct->field);
  MRIcopyPulseParameters(mri_template,dct->field);

  //crs = MatrixAlloc(4,1,MATRIX_REAL);
  //crs->rptr[4][1] = 1;
  for(n=0; n < dct->nfields; n++){
    dim = dct->dim[n];
    period = dct->period[n];
    offset = dct->offsetCRS->rptr[dim][1];
    printf("n = %2d, dim=%d period=%6.2f (%6.4f) offset = %7.3f\n",n,dim,period,1.0/period,offset);
    fflush(stdout);
    for(c=0; c < mri_template->width; c++){
      //crs->rptr[1][1] = c;
      for(r=0; r < mri_template->height; r++){
        //crs->rptr[2][1] = r;
	for(s=0; s < mri_template->depth; s++){
	  if(dct->mask && MRIgetVoxVal(dct->mask,c,r,s,0)<0.5) continue;
	  //crs->rptr[3][1] = s;
	  if(period == 0){
	    MRIsetVoxVal(dct->field,c,r,s,n,1);	    
	    continue;
	  }
	  switch(dim){
	  case 1: t = (c-offset)*dct->field->xsize; break;
	  case 2: t = (r-offset)*dct->field->ysize; break;
	  case 3: t = (s-offset)*dct->field->zsize; break;
	  }
	  //ras = MatrixMultiplyD(vox2ras,crs,ras);
	  //v = cos(2*M_PI*(ras->rptr[dim][1]-offset)/period);
	  v = cos(2*M_PI*t/period);
	  MRIsetVoxVal(dct->field,c,r,s,n,v);
	}
      }
    }
  }

  MatrixFree(&vox2ras);
  MatrixFree(&ras2vox);
  return(0);
}

int DCTspecLPF(DCTFS *dct)
{
  int dim,nvox,nf;
  double res,F, deltaF;

  // alloc way more than needed (hopefully:)
  dct->period = (double*)calloc(sizeof(double),1000);
  dct->dim    = (int*)calloc(sizeof(int),1000);

  dct->nfields = 0;

  // Do not include a DC term with log fitting
  //dct->period[dct->nfields] = 0;
  //dct->dim[dct->nfields] = 1;
  //dct->nfields ++;

  for(dim=1; dim <= 3; dim++){
    /* For the given dim, determine the size of the dim and the resolution
     to determine the frequencies to use. If a bounding box is specified,
     then the dim size is taken from that. */
    switch(dim){
    case 1:
      if(dct->bbfit) nvox = dct->bbfit->dx;
      else           nvox = dct->mri_template->width;
      res = dct->mri_template->xsize;
      break;
    case 2:
      if(dct->bbfit) nvox = dct->bbfit->dy;
      else           nvox = dct->mri_template->height;
      res = dct->mri_template->ysize;
      break;
    case 3:
      if(dct->bbfit) nvox = dct->bbfit->dz;
      else           nvox = dct->mri_template->depth;
      res = dct->mri_template->zsize;
      break;
    }
    // deltaF is the frequency resolution.  Using 0.5 allows for a
    // near linear trend across the image (this is what SPM
    // does). This does mean that the zero point has to be at the
    // corner of the volume or else the design matrix becomes highly
    // ill-conditioned.
    deltaF = 0.5/(nvox*res); // 0.5/FoV
    printf("dim=%d nvox=%d res=%g deltaF=%g\n",dim,nvox,res,deltaF);

    // count the number of fields for this dim
    F = deltaF; // Start at the lowest frequency
    nf = 0;
    while(F < 1.0/dct->lpcutoff){
      printf("  nf = %d  %g  %g\n",dct->nfields,F,1/F);
      dct->period[dct->nfields] = 1.0/F;
      dct->dim[dct->nfields] = dim;
      dct->nfields ++;
      nf++;
      F += deltaF;
    }
    printf("  dim = %d, nf = %d\n",dim,nf);
  }
  printf("nfields = %d\n",dct->nfields);
  fflush(stdout);
  return(0);
}

DCTFS *DCTalloc(void)
{
  DCTFS *dct = (DCTFS*) calloc(sizeof(DCTFS),1);
  dct->offsetRAS = MatrixAlloc(4,1,MATRIX_REAL);
  dct->offsetRAS->rptr[4][1] = 1;
  dct->offsetCRS = MatrixAlloc(4,1,MATRIX_REAL);
  dct->offsetCRS->rptr[4][1] = 1;
  return(dct);
}

MRI *MRIapplyBiasField(MRI *input, MRI *bf, MRI *seg, MRI *mask, double targwmval, double thresh, MRI *bc)
{
  int nwm,c;
  double wmsum, wmmean, wmscale;

  if(bc==NULL){
    bc = MRIallocSequence(input->width, input->height, 
			  input->depth, MRI_FLOAT, input->nframes);
    MRIcopyHeader(input,bc);
    MRIcopyPulseParameters(input,bc);
  }

  wmsum = 0;
  nwm = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : wmsum, nwm)
  #endif
  for(c=0; c < input->width; c++){
    int r, s,segid;
    double val, b;
    for(r=0; r < input->height; r++){
      for(s=0; s < input->depth; s++){
	MRIsetVoxVal(bc,c,r,s,0,0);
        if(mask && MRIgetVoxVal(mask,c,r,s,0) <0.5)continue;
        b = MRIgetVoxVal(bf,c,r,s,0);
	if(b==0) continue;
        val = MRIgetVoxVal(input,c,r,s,0);
	if(val <= thresh) continue; // matches MRIbiasFieldCorLog()
        val = val/(b+FLT_MIN);
        MRIsetVoxVal(bc,c,r,s,0,val);
	if(seg){
          segid = MRIgetVoxVal(seg,c,r,s,0);
	  switch (segid) {
            case Left_Cerebral_White_Matter: 
            case Right_Cerebral_White_Matter:
   	    case 251:
	    case 252:
	    case 253:
	    case 254:
	    case 255:
            wmsum += val;
  	    nwm ++;
	    break;
	  }
        }
      }
    }
  }

  if(seg==NULL) return(bc);

  wmmean = wmsum/nwm;
  wmscale = targwmval/wmmean;
  printf("nwm = %d, wmsum = %g, wmmean = %g, scale = %g\n",nwm,wmsum,wmmean,wmscale);


  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : wmsum, nwm)
  #endif
  for(c=0; c < input->width; c++){
    int r, s;
    double val,b;
    for(r=0; r < input->height; r++){
      for(s=0; s < input->depth; s++){
        if(mask && MRIgetVoxVal(mask,c,r,s,0) <0.5) continue;
        val = MRIgetVoxVal(bc,c,r,s,0);
        MRIsetVoxVal(bc,c,r,s,0,val*wmscale);
        b = MRIgetVoxVal(bf,c,r,s,0);
        MRIsetVoxVal(bf,c,r,s,0,b/wmscale);
      }
    }
  }

  return(bc);

}
