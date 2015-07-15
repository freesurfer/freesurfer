/**
 * @file  dummy.c
 * @brief Computes whole-brain connectivity
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2015/07/15 21:54:08 $
 *    $Revision: 1.1 $
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
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "mrisurf.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "timer.h"
#include "mrimorph.h"
#include "fmriutils.h"

typedef struct {
  MRI *fvol, *flh, *frh;
  MRIS *lh, *rh, *lhsph, *rhsph;
  LABEL *lhlabel, *rhlabel;
  MRI *volmask, *lhmask, *rhmask;
  int nvolmask, nlhmask, nrhmask, ntot;
  MRI *volcon, *lhcon, *rhcon;
  MRI *coordtype, *vertexno, *xyz;
  MRI *f, *fnorm;
  double rhothresh;
  int DoDist;
  double distthresh;
  MRI *con;
  int DoMat;
  MATRIX *M;
} WBC;

MRI *WholeBrainCon(WBC *wbc);
int WBCfinish(WBC *wbc);
int WBCprep(WBC *wbc);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_wbc.c,v 1.1 2015/07/15 21:54:08 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

typedef struct {
  char *fvol;
  char *volmask;
  char *flh, *lhsurface, *lhlabel, *lhmask;
  char *frh, *rhsurface, *rhlabel, *rhmask;
  double rhothresh, distthresh;
  char *volcon, *lhcon, *rhcon;
} CMDARGS;

CMDARGS *cmdargs;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err;
  WBC *wbc;

  cmdargs = (CMDARGS *)calloc(sizeof(CMDARGS),1);

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

  wbc = calloc(sizeof(WBC),1);
  wbc->rhothresh = cmdargs->rhothresh;
  wbc->distthresh = 10;
  wbc->DoDist = 0;
  wbc->DoMat = 0;

  if(cmdargs->fvol){
    wbc->fvol = MRIread(cmdargs->fvol);
    if(wbc->fvol == NULL) exit(1);
  }
  if(cmdargs->volmask){
    wbc->volmask = MRIread(cmdargs->volmask);
    if(wbc->volmask == NULL) exit(1);
  }

  WBCprep(wbc);
  WholeBrainCon(wbc);
  WBCfinish(wbc);
  err = MRIwrite(wbc->volcon,cmdargs->volcon);
  if(err) exit(1);

  exit(0);
}

/* -------------------------------------------------------- */
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

    else if (!strcasecmp(option, "--fvol")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->fvol = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--volmask")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volmask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lh")) {
      if(nargc < 2) CMDargNErr(option,2);
      cmdargs->flh = pargv[0];
      cmdargs->lhsurface = pargv[1];
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--lhmask")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhmask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhlabel")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhlabel = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rh")) {
      if(nargc < 2) CMDargNErr(option,2);
      cmdargs->frh = pargv[0];
      cmdargs->rhsurface = pargv[1];
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--rhmask")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhmask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhlabel")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhlabel = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rho")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->rhothresh);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--dist")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->distthresh);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--volcon")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volcon = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhcon")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhcon = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhcon")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhcon = pargv[0];
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
/* -------------------------------------------------------- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --fvol volfile : functional volume\n");
  printf("   --volmask volfile : mask for functional volume \n");
  printf("\n");
  printf("   --lh lhffile lh.surface : functional surface lh\n");
  printf("   --lhmask lhfile : mask for lh functional surface\n");
  printf("   --lhlabel lhlabel : label mask for lh functional surface\n");
  printf("\n");
  printf("   --rh rhffile rh.surface: functional surface rh\n");
  printf("   --rhmask rhfile : mask for rh functional surface\n");
  printf("   --rhlabel rhlabel : label mask for rh functional surface\n");
  printf("\n");
  printf("   --rho rhothresh\n");
  printf("   --dist distthresh\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* -------------------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void check_options(void) {
  if(cmdargs->fvol == NULL && cmdargs->flh == NULL && cmdargs->frh == NULL){
    printf("ERROR: no input\n");
    exit(1);
  }

  return;
}
/* -------------------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  return;
}

/*******************************************************************************/
MRI *WholeBrainCon(WBC *wbc)
{
  MRI **conth;
  int nthreads;
  int nperthread0,npairs,*nperthread,ntot;
  int **k1thread,**k2thread,threadno;
  int k1,k2,nthpair; // note: these are redefined in loop
  struct timeb timer;
  double **pf;

  // Load time courses into a pointer structure
  pf = (double **) calloc(sizeof(double *),wbc->ntot);
  for(k1 = 0; k1 < wbc->ntot; k1++){
    int t;
    pf[k1] = (double *) calloc(sizeof(double),wbc->f->nframes);
    for(t=0; t < wbc->f->nframes; t++) 
      pf[k1][t] = MRIgetVoxVal(wbc->fnorm,k1,0,0,t);    
  }

  // This is for testing. In most apps, M will be huge
  if(wbc->DoMat) wbc->M = MatrixAlloc(wbc->ntot,wbc->ntot,MATRIX_REAL);

  nthreads = 1;
  #ifdef _OPENMP
  nthreads = omp_get_max_threads();
  #endif
  npairs = wbc->ntot*(wbc->ntot-1)/2;
  nperthread0 = nint((double)npairs/nthreads - 1);
  printf("ntot = %d, nthreads = %d, npairs = %d, nperthread0 = %d\n",wbc->ntot,nthreads,npairs,nperthread0);

  conth = (MRI **) calloc(sizeof(MRI*),nthreads);
  nperthread = (int *) calloc(sizeof(int),nthreads);
  ntot = 0;
  for(threadno=0; threadno < nthreads; threadno++){
    int nf = 1;
    if(wbc->DoDist) nf = 3;
    conth[threadno] = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, nf);
    if(conth[threadno]==NULL) {
      printf("ERROR: conth %d: could not alloc\n",threadno);
      return(NULL);
    }
    nperthread[threadno] = nperthread0; // number of pairs per thread
    ntot += nperthread0;
  }
  if(ntot != npairs) nperthread[nthreads-1] += npairs-ntot;

  k1thread = (int **)calloc(sizeof(int*),nthreads);
  k2thread = (int **)calloc(sizeof(int*),nthreads);
  for(threadno=0; threadno < nthreads; threadno++){
    k1thread[threadno] = (int *)calloc(sizeof(int),nperthread[threadno]);
    k2thread[threadno] = (int *)calloc(sizeof(int),nperthread[threadno]);
  }

  // This is just a test
  ntot = 0;
  for(threadno=0; threadno < nthreads; threadno++){
    printf("thread %d %d\n",threadno,nperthread[threadno]);
    ntot += nperthread[threadno];
  }
  printf("ntot = %d vs npairs %d\n",ntot,npairs); // should be equal

  // Assign k1,k2 pairs to each thread
  nthpair = 0;
  threadno = 0;
  for(k1 = 0; k1 < wbc->ntot; k1++){
    for(k2 = k1+1; k2 < wbc->ntot; k2++){
      k1thread[threadno][nthpair] = k1;
      k2thread[threadno][nthpair] = k2;
      nthpair = nthpair + 1;
      if(nthpair >= nperthread[threadno]){
	nthpair = 0;
	threadno ++;
      }
    }
  }

  TimerStart(&timer);
  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
  for(threadno = 0; threadno < nthreads; threadno ++){
    int  k1, k2, t, n, thno,nthpair, q, ct1, ct2;
    double rho,dx,dy,dz,dist,*pf1,*pf2,x1,y1,z1,x2,y2,z2;
    thno = threadno;
    #ifdef _OPENMP
    thno = omp_get_thread_num(); // actual thread number
    #endif
    for(nthpair = 0; nthpair < nperthread[thno]; nthpair++){
      k1 = k1thread[thno][nthpair];
      k2 = k2thread[thno][nthpair];

      rho = 0;
      pf1 = pf[k1];
      pf2 = pf[k2];
      for(t=0; t < wbc->f->nframes; t++){
	rho += (*pf1) * (*pf2);
	pf1++;
	pf2++;
      }
      if(wbc->M != NULL){
	wbc->M->rptr[k1+1][k2+1] = rho;
	wbc->M->rptr[k2+1][k1+1] = rho;
      }
      if(fabs(rho) < wbc->rhothresh) continue;

      n = MRIgetVoxVal(conth[thno],k1,0,0,0);
      MRIsetVoxVal(conth[thno],k1,0,0,0,n+1);
      n = MRIgetVoxVal(conth[thno],k2,0,0,0);
      MRIsetVoxVal(conth[thno],k2,0,0,0,n+1);

      if(!wbc->DoDist) continue;

      ct1 = MRIgetVoxVal(wbc->coordtype,k1,0,0,0);
      ct2 = MRIgetVoxVal(wbc->coordtype,k2,0,0,0);
      if(ct1 == ct2){
	x1 = MRIgetVoxVal(wbc->xyz,k1,0,0,0);
	y1 = MRIgetVoxVal(wbc->xyz,k1,0,0,1);
	z1 = MRIgetVoxVal(wbc->xyz,k1,0,0,2);
	x2 = MRIgetVoxVal(wbc->xyz,k2,0,0,0);
	y2 = MRIgetVoxVal(wbc->xyz,k2,0,0,1);
	z2 = MRIgetVoxVal(wbc->xyz,k2,0,0,2);
	dx = x1-x2;
	dy = y1-y2;
	dz = z1-z2;
	dist = sqrt(dx*dx + dy*dy + dz*dz);
	if(dist < wbc->distthresh) q = 1; // short dist
	else                       q = 2; // long dist
      }
      else q = 2; // long dist
      n = MRIgetVoxVal(conth[thno],k1,0,0,q);
      MRIsetVoxVal(conth[thno],k1,0,0,q,n+1);
      n = MRIgetVoxVal(conth[thno],k2,0,0,q);
      MRIsetVoxVal(conth[thno],k2,0,0,q,n+1);

    } // k2
  } // k1

  // Sum up the threads
  for(threadno=0; threadno < nthreads; threadno++)
    MRIadd(wbc->con,conth[threadno],wbc->con);

  // Divide number of connections by total possible
  printf("Scaling by ntot %d\n",wbc->ntot);
  MRImultiplyConst(wbc->con,1.0/wbc->ntot,wbc->con);

  // Clean up
  for(threadno=0; threadno < nthreads; threadno++){
    MRIfree(&conth[threadno]);
    free(k1thread[threadno]);
    free(k2thread[threadno]);
  }
  free(conth);
  free(k1thread);
  free(k2thread);
  for(k1 = 0; k1 < wbc->ntot; k1++) free(pf[k1]);
  free(pf);

  printf(" t = %6.4f, ntot = %d\n",TimerStop(&timer)/1000.0,wbc->ntot);fflush(stdout);
  return(wbc->con);
}

int WBCprep(WBC *wbc)
{
  int nthvox, c, r, s, t, k, vtxno;
  MATRIX *V, *crs, *xyz=NULL;
  double val;

  if(wbc->lhlabel) wbc->lhmask = MRISlabel2Mask(wbc->lh,wbc->lhlabel,NULL);
  if(wbc->rhlabel) wbc->rhmask = MRISlabel2Mask(wbc->rh,wbc->lhlabel,NULL);

  wbc->nvolmask = 0;
  wbc->nlhmask = 0;
  wbc->nrhmask = 0;
  if(wbc->fvol){
    if(wbc->volmask) wbc->nvolmask = MRIcountAboveThreshold(wbc->volmask, 0.5);
    else wbc->nvolmask = wbc->fvol->width * wbc->fvol->height * wbc->fvol->depth;
  }
  if(wbc->flh){
    if(wbc->lhmask) wbc->nlhmask  = MRIcountAboveThreshold(wbc->lhmask, 0.5);
    else            wbc->nlhmask  = wbc->flh->width;
  }
  if(wbc->frh){
    if(wbc->rhmask) wbc->nrhmask  = MRIcountAboveThreshold(wbc->rhmask, 0.5);
    else            wbc->nrhmask  = wbc->frh->width;
  }

  wbc->ntot = wbc->nvolmask + wbc->nlhmask + wbc->nrhmask;
  printf("nmask %d %d %d   %d\n",wbc->nvolmask,wbc->nlhmask,wbc->nrhmask,wbc->ntot);

  wbc->coordtype = MRIalloc(wbc->ntot,1,1,MRI_INT);
  wbc->vertexno  = MRIalloc(wbc->ntot,1,1,MRI_INT);
  wbc->xyz       = MRIallocSequence(wbc->ntot,1,1,MRI_FLOAT,3);
  wbc->f = MRIallocSequence(wbc->ntot,1,1,MRI_FLOAT,wbc->fvol->nframes);

  nthvox = 0;
  if(wbc->fvol){
    V = MRIxfmCRS2XYZtkreg(wbc->fvol);
    crs = MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
    for(s=0; s < wbc->fvol->depth; s++){
      for(c=0; c < wbc->fvol->width; c++){
	for(r=0; r < wbc->fvol->height; r++){
	  if(wbc->volmask && MRIgetVoxVal(wbc->volmask,c,r,s,0) < 0.5) continue;
	  MRIsetVoxVal(wbc->coordtype,nthvox,0,0,0, 0);
	  crs->rptr[1][1] = c;
	  crs->rptr[2][1] = r;
	  crs->rptr[3][1] = s;
	  xyz = MatrixMultiplyD(V,crs,xyz);
	  for(k=0; k < 3; k++) MRIsetVoxVal(wbc->xyz,nthvox,0,0,k, xyz->rptr[k+1][1]);
	  for(t=0; t < wbc->fvol->nframes; t++){
	    val = MRIgetVoxVal(wbc->fvol,c,r,s,t);
	    MRIsetVoxVal(wbc->f,nthvox,0,0,t,val);
	  } // time
	  nthvox ++;
	} // row
      } // col
    } // slice
    MatrixFree(&V);
    MatrixFree(&xyz);
    MatrixFree(&crs);    
  }
  if(wbc->flh){
    for(vtxno = 0; vtxno < wbc->flh->width; vtxno++){
      if(wbc->lhmask && MRIgetVoxVal(wbc->lhmask,vtxno,0,0,0) < 0.5) continue;
      MRIsetVoxVal(wbc->vertexno,nthvox,0,0,0, vtxno);
      MRIsetVoxVal(wbc->coordtype,nthvox,0,0,0, 1);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,0, wbc->lh->vertices[vtxno].x);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,1, wbc->lh->vertices[vtxno].y);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,2, wbc->lh->vertices[vtxno].z);
      for(t=0; t < wbc->flh->nframes; t++){
	val = MRIgetVoxVal(wbc->flh,vtxno,0,0,t);
	MRIsetVoxVal(wbc->f,nthvox,0,0,t,val);
      } // time
    } // vertex
  }
  if(wbc->frh){
    for(vtxno = 0; vtxno < wbc->frh->width; vtxno++){
      if(wbc->rhmask && MRIgetVoxVal(wbc->rhmask,vtxno,0,0,0) < 0.5) continue;
      MRIsetVoxVal(wbc->vertexno,nthvox,0,0,0, vtxno);
      MRIsetVoxVal(wbc->coordtype,nthvox,0,0,0, 2);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,0, wbc->rh->vertices[vtxno].x);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,1, wbc->rh->vertices[vtxno].y);
      MRIsetVoxVal(wbc->xyz,nthvox,0,0,2, wbc->rh->vertices[vtxno].z);
      for(t=0; t < wbc->frh->nframes; t++){
	val = MRIgetVoxVal(wbc->frh,vtxno,0,0,t);
	MRIsetVoxVal(wbc->f,nthvox,0,0,t,val);
      } // time
    } // vertex
  }

  // normalize by stddev
  wbc->fnorm = MRIframeNorm(wbc->f, NULL, NULL);  

  int nf = 1;
  if(wbc->DoDist) nf = 3;
  wbc->con = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, nf);
  if(wbc->con==NULL) {
    printf("ERROR: WBCprep(): could not alloc con\n");
    return(1);
  }

  return(0);
}

int WBCfinish(WBC *wbc)
{
  int nf = 1,nthvox, c,r,s,t,vtxno;
  double val;

  if(wbc->DoDist) nf = 3;

  nthvox = 0;
  if(wbc->fvol){
    wbc->volcon = MRIallocSequence(wbc->fvol->width,wbc->fvol->height,wbc->fvol->depth,MRI_FLOAT,nf);
    MRIcopyHeader(wbc->fvol,wbc->volcon);
    MRIcopyPulseParameters(wbc->fvol,wbc->volcon);
    for(s=0; s < wbc->fvol->depth; s++){
      for(c=0; c < wbc->fvol->width; c++){
	for(r=0; r < wbc->fvol->height; r++){
	  if(wbc->volmask && MRIgetVoxVal(wbc->volmask,c,r,s,0) < 0.5) continue;
	  for(t=0; t < nf; t++){
	    val = MRIgetVoxVal(wbc->con,nthvox,0,0,t);
	    MRIsetVoxVal(wbc->volcon,c,r,s,t,val);
	  } // time
	  nthvox ++;
	} // row
      } // col
    } // slice
  }
  if(wbc->flh){
    MRIcopyHeader(wbc->flh,wbc->lhcon);
    MRIcopyPulseParameters(wbc->flh,wbc->lhcon);
    for(vtxno = 0; vtxno < wbc->flh->width; vtxno++){
      if(wbc->lhmask && MRIgetVoxVal(wbc->lhmask,vtxno,0,0,0) < 0.5) continue;
      for(t=0; t < wbc->flh->nframes; t++){
	val = MRIgetVoxVal(wbc->con,nthvox,0,0,t);
	MRIsetVoxVal(wbc->lhcon,vtxno,0,0,t,val);
      } // time
    } // vertex
  }
  if(wbc->frh){
    MRIcopyHeader(wbc->frh,wbc->rhcon);
    MRIcopyPulseParameters(wbc->frh,wbc->rhcon);
    for(vtxno = 0; vtxno < wbc->frh->width; vtxno++){
      if(wbc->rhmask && MRIgetVoxVal(wbc->rhmask,vtxno,0,0,0) < 0.5) continue;
      for(t=0; t < wbc->frh->nframes; t++){
	val = MRIgetVoxVal(wbc->con,nthvox,0,0,t);
	MRIsetVoxVal(wbc->rhcon,vtxno,0,0,t,val);
      } // time
    } // vertex
  }
  return(0);
}

