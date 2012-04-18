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
 *    $Date: 2012/04/18 19:14:50 $
 *    $Revision: 1.2 $
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

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;
MRI *MRIpolyfitBiasField(MRI *vol2, int order2, MRI *seg2, MRI *bias2);

static char vcid[] = "$Id: mri_fit_bias.c,v 1.2 2012/04/18 19:14:50 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *TempVolFile=NULL;
char *subject, *SUBJECTS_DIR;
char *srcfile=NULL, *trgfile=NULL;
char *maskfile = "brainmask.mgz";
char *segfile = "aseg.mgz";
MRI *src, *seg, *segerode, *mask, *trg;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,c,r,s,k,nmask,err;
  char tmpstr[2000];
  double sumval,val;

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

  segerode = MRIerodeSegmentation(seg, NULL, 1, 0);
  trg = MRIpolyfitBiasField(src, 3, segerode, NULL);

  // Normalize mri so that mean in mask is 110
  sumval = 0;
  nmask = 0;
  for(s=0; s < segerode->depth; s++){
    for(c=0; c < segerode->width; c++){
      for(r=0; r < segerode->height; r++){
	k = MRIgetVoxVal(segerode,c,r,s,0);
	if(k!=2 && k!=41) continue;
	nmask ++;
	sumval += MRIgetVoxVal(trg,c,r,s,0);
      }
    }
  }
  sumval = sumval/nmask; // now it's the mean
  printf("In-mask mean %g  %d\n",sumval,nmask);
  // Rescale all values
  for(s=0; s < segerode->depth; s++){
    for(c=0; c < segerode->width; c++){
      for(r=0; r < segerode->height; r++){
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

  exit(0);



  return 0;
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
  XtX  = MatrixMultiplyD(Xt,X,NULL);
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



