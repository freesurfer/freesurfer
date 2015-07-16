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
 *    $Date: 2015/07/16 19:41:39 $
 *    $Revision: 1.4 $
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
  int nvolmask, nlhmask, nrhmask, ntot,nframes;
  MRI *volcon, *lhcon, *rhcon;
  MRI *volconS, *lhconS, *rhconS;
  MRI *volconL, *lhconL, *rhconL;
  MRI *coordtype, *vertexno, *xyz;
  MRI *f, *fnorm;
  double rholist[100];
  int nrholist;
  int DoDist;
  double distthresh;
  MRI *con,*conS,*conL;
  int DoMat;
  MATRIX *M;
} WBC;

MRI *WholeBrainCon(WBC *wbc);
int WBCfinish(WBC *wbc);
int WBCprep(WBC *wbc);
int Index2UpperSubscript(int N, long i, int *r, int *c);
int Index2UpperSubscriptTest(int seed);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_wbc.c,v 1.4 2015/07/16 19:41:39 greve Exp $";
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
  double rholist[100];
  int nrholist;
  double distthresh;
  int DoDist;
  char *volcon, *lhcon, *rhcon;
  char *volconS, *lhconS, *rhconS;
  char *volconL, *lhconL, *rhconL;
  int DoMat;
  char *matfile;
} CMDARGS;

CMDARGS *cmdargs;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err,n;
  WBC *wbc;

  cmdargs = (CMDARGS *)calloc(sizeof(CMDARGS),1);
  cmdargs->distthresh = 10;
  cmdargs->DoDist = 0;
  cmdargs->DoMat = 0;
  cmdargs->nrholist = 0;

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
  wbc->distthresh = cmdargs->distthresh;
  wbc->DoDist = cmdargs->DoDist;
  wbc->DoMat = cmdargs->DoMat;
  for(n=0; n < cmdargs->nrholist; n++)
    wbc->rholist[n] = cmdargs->rholist[n];
  wbc->nrholist = cmdargs->nrholist;

  if(cmdargs->fvol){
    wbc->fvol = MRIread(cmdargs->fvol);
    if(wbc->fvol == NULL) exit(1);
  }
  if(cmdargs->volmask){
    wbc->volmask = MRIread(cmdargs->volmask);
    if(wbc->volmask == NULL) exit(1);
  }
  if(cmdargs->flh){
    wbc->flh = MRIread(cmdargs->flh);
    if(wbc->flh == NULL) exit(1);
    wbc->lh  = MRISread(cmdargs->lhsurface);
    if(wbc->lh == NULL) exit(1);
    if(cmdargs->lhlabel){
      wbc->lhlabel = LabelRead(NULL, cmdargs->lhlabel);
      if(wbc->lhlabel == NULL) exit(1);
    }
    if(cmdargs->lhmask){
      wbc->lhmask = MRIread(cmdargs->lhmask);
      if(wbc->lhmask == NULL) exit(1);
    }
  }
  if(cmdargs->frh){
    wbc->frh = MRIread(cmdargs->frh);
    if(wbc->frh == NULL) exit(1);
    wbc->rh  = MRISread(cmdargs->rhsurface);
    if(wbc->rh == NULL) exit(1);
    if(cmdargs->rhlabel){
      wbc->rhlabel = LabelRead(NULL, cmdargs->rhlabel);
      if(wbc->rhlabel == NULL) exit(1);
    }
    if(cmdargs->rhmask){
      wbc->rhmask = MRIread(cmdargs->rhmask);
      if(wbc->rhmask == NULL) exit(1);
    }
  }

  WBCprep(wbc);
  WholeBrainCon(wbc);
  WBCfinish(wbc);

  if(cmdargs->volcon){
    err = MRIwrite(wbc->volcon,cmdargs->volcon);
    if(err) exit(1);
    if(wbc->DoDist){
      err = MRIwrite(wbc->volconS,cmdargs->volconS);
      if(err) exit(1);
      err = MRIwrite(wbc->volconL,cmdargs->volconL);
      if(err) exit(1);
    }
  }
  if(cmdargs->lhcon){
    err = MRIwrite(wbc->lhcon,cmdargs->lhcon);
    if(err) exit(1);
    if(wbc->DoDist){
      err = MRIwrite(wbc->lhconS,cmdargs->lhconS);
      if(err) exit(1);
      err = MRIwrite(wbc->lhconL,cmdargs->lhconL);
      if(err) exit(1);
    }
  }
  if(cmdargs->rhcon){
    err = MRIwrite(wbc->rhcon,cmdargs->rhcon);
    if(err) exit(1);
    if(wbc->DoDist){
      err = MRIwrite(wbc->rhconS,cmdargs->rhconS);
      if(err) exit(1);
      err = MRIwrite(wbc->rhconL,cmdargs->rhconL);
      if(err) exit(1);
    }
  }

  if(wbc->DoMat){
    err = MatrixWriteTxt(cmdargs->matfile,wbc->M);
    if(err) exit(1);    
  }

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
      sscanf(pargv[0],"%lf",&cmdargs->rholist[cmdargs->nrholist]);
      cmdargs->nrholist++;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--dist")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->distthresh);
      cmdargs->DoDist = 1;
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
    else if (!strcasecmp(option, "--volconS")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volconS = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhconS")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhconS = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhconS")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhconS = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--volconL")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->volconL = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lhconL")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->lhconL = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rhconL")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->rhconL = pargv[0];
      nargsused = 1;
    } 

    else if (!strcasecmp(option, "--mat")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->matfile = pargv[0];
      cmdargs->DoMat = 1;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--uppersub-test")) {
      if(nargc < 1) CMDargNErr(option,1);
      int N, err;
      sscanf(pargv[0],"%d",&N);
      err = Index2UpperSubscriptTest(N);
      exit(err);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      int nthreads;
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
  printf("   --threads nthreads\n");
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
  if(cmdargs->fvol == NULL && cmdargs->volcon != NULL){
    printf("ERROR: cannot have --volcon without --fvol\n");
    exit(1);
  }
  if(cmdargs->fvol != NULL && cmdargs->volcon == NULL){
    printf("ERROR: need --volcon output with --fvol\n");
    exit(1);
  }
  if(cmdargs->flh == NULL && cmdargs->lhcon != NULL){
    printf("ERROR: cannot have --lhcon without --lh\n");
    exit(1);
  }
  if(cmdargs->flh != NULL && cmdargs->lhcon == NULL){
    printf("ERROR: need --lhcon output with --lh\n");
    exit(1);
  }
  if(cmdargs->frh == NULL && cmdargs->rhcon != NULL){
    printf("ERROR: cannot have --rhcon without --rh\n");
    exit(1);
  }
  if(cmdargs->frh != NULL && cmdargs->rhcon == NULL){
    printf("ERROR: need --rhcon output with --rh\n");
    exit(1);
  }
  if(cmdargs->nrholist == 0){
    cmdargs->rholist[0] = 0.2;
    cmdargs->nrholist = 1;
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
  MRI **conth,**conSth,**conLth;
  int nthreads,nthrho;
  long nperthread0,*nperthread,ntot;
  int threadno;
  long npairs,ia,ib;
  int k1; // note: these are redefined in loop
  struct timeb timer;
  double **pf;
  int *k1a,*k1b,*k2a,*k2b;

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
  npairs = (long)wbc->ntot*(wbc->ntot-1)/2;
  nperthread0 = nint((double)npairs/nthreads - 1);
  printf("ntot = %d, nthreads = %d, npairs = %ld, nperthread0 = %ld\n",wbc->ntot,nthreads,npairs,nperthread0);

  nperthread = (long *) calloc(sizeof(long),nthreads);
  ntot = 0;
  for(threadno=0; threadno < nthreads; threadno++){
    nperthread[threadno] = nperthread0; // number of pairs per thread
    ntot += nperthread0;
  }
  if(ntot != npairs) nperthread[nthreads-1] += npairs-ntot;
  // This is just a test
  ntot = 0;
  for(threadno=0; threadno < nthreads; threadno++){
    printf("thread %d %ld\n",threadno,nperthread[threadno]);
    ntot += nperthread[threadno];
  }
  printf("ntotpairs = %ld vs npairs %ld, diff %ld\n",ntot,npairs,ntot-npairs); // should be equal

  conth = (MRI **) calloc(sizeof(MRI*),nthreads);
  if(wbc->DoDist){
    conSth = (MRI **) calloc(sizeof(MRI*),nthreads);
    conLth = (MRI **) calloc(sizeof(MRI*),nthreads);
  }
  for(threadno=0; threadno < nthreads; threadno++){
    conth[threadno] = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
    if(wbc->DoDist){
      conSth[threadno] = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
      conLth[threadno] = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
    }
  }

  k1a = (int *)calloc(sizeof(int),nthreads);
  k1b = (int *)calloc(sizeof(int),nthreads);
  k2a = (int *)calloc(sizeof(int),nthreads);
  k2b = (int *)calloc(sizeof(int),nthreads);
  ia = 0;
  for(threadno=0; threadno < nthreads; threadno++){
    ib = ia + nperthread[threadno];
    Index2UpperSubscript(wbc->ntot, ia, &k1a[threadno], &k2a[threadno]);
    Index2UpperSubscript(wbc->ntot, ib-1, &k1b[threadno], &k2b[threadno]);
    printf("thread %d %12ld %12ld   %7d %7d   %7d %7d\n",threadno,ia,ib,
	   k1a[threadno],k2a[threadno], k1b[threadno],k2b[threadno]);
    ia = ib;
  }

  printf("rho thresholds (%d): ",wbc->nrholist);
  for(nthrho = 0; nthrho < wbc->nrholist; nthrho++)
    printf("%lf ",wbc->rholist[nthrho]);
  printf("\n");
  
  printf("Starting WBC loop\n"); fflush(stdout);
  TimerStart(&timer);
  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
  for(threadno = 0; threadno < nthreads; threadno ++){
    int  k1, k2, t, n, thno, q, ct1, ct2;
    int k1start,k1stop,k2start,k2stop,k2min,k2max,nthrho;
    double rho,dx,dy,dz,dist,*pf1,*pf2,x1,y1,z1,x2,y2,z2;
    double rhothresh;
    MRI *conDth;
    thno = threadno;
    #ifdef _OPENMP
    thno = omp_get_thread_num(); // actual thread number
    #endif

    k1start = k1a[thno];
    k1stop  = k1b[thno];
    k2start = k2a[thno];
    k2stop  = k2b[thno];

    for(k1=k1start; k1<=k1stop; k1++){
      if(k1 == k1start) k2min = k2start;
      else              k2min = k1+1;
      if(k1 == k1stop)  k2max = k2stop;
      else              k2max = wbc->ntot-1;
      for(k2=k2min; k2<=k2max; k2++){

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

	for(nthrho = 0; nthrho < wbc->nrholist; nthrho++){
	  rhothresh = wbc->rholist[nthrho];
	  if(fabs(rho) < rhothresh) continue;
	
	  n = MRIgetVoxVal(conth[thno],k1,0,0,nthrho);
	  MRIsetVoxVal(conth[thno],k1,0,0,nthrho,n+1);
	  n = MRIgetVoxVal(conth[thno],k2,0,0,nthrho);
	  MRIsetVoxVal(conth[thno],k2,0,0,nthrho,n+1);
	} // rho thresh
	
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
	if(q == 1) conDth = conSth[thno];
	if(q == 2) conDth = conLth[thno];

	for(nthrho = 0; nthrho < wbc->nrholist; nthrho++){
	  rhothresh = wbc->rholist[nthrho];
	  if(fabs(rho) < rhothresh) continue;
	  n = MRIgetVoxVal(conDth,k1,0,0,nthrho);
	  MRIsetVoxVal(conDth,k1,0,0,nthrho,n+1);
	  n = MRIgetVoxVal(conDth,k2,0,0,nthrho);
	  MRIsetVoxVal(conDth,k2,0,0,nthrho,n+1);
	} // rho thresh

      } // k2
    } // k1
  } // thread

  // Sum up the threads
  for(threadno=0; threadno < nthreads; threadno++){
    MRIadd(wbc->con,conth[threadno],wbc->con);
    if(wbc->DoDist){
      MRIadd(wbc->conS,conSth[threadno],wbc->conS);
      MRIadd(wbc->conL,conLth[threadno],wbc->conL);
    }
  }

  // Divide number of connections by total possible
  printf("Scaling by ntot %d\n",wbc->ntot);
  MRImultiplyConst(wbc->con,1.0/wbc->ntot,wbc->con);
  if(wbc->DoDist){
    MRImultiplyConst(wbc->conS,1.0/wbc->ntot,wbc->conS);
    MRImultiplyConst(wbc->conL,1.0/wbc->ntot,wbc->conL);
  }

  // Clean up
  for(threadno=0; threadno < nthreads; threadno++){
    MRIfree(&conth[threadno]);
    if(wbc->DoDist){
      MRIfree(&conSth[threadno]);
      MRIfree(&conLth[threadno]);
    }
  }
  free(k1a);
  free(k1b);
  free(k2a);
  free(k2b);
  free(conth);
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

  if(wbc->fvol)     wbc->nframes = wbc->fvol->nframes;
  else if(wbc->flh) wbc->nframes = wbc->flh->nframes;
  else if(wbc->frh) wbc->nframes = wbc->frh->nframes;

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
  wbc->f = MRIallocSequence(wbc->ntot,1,1,MRI_FLOAT,wbc->nframes);

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
	  for(t=0; t < wbc->nframes; t++){
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

  wbc->con = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
  if(wbc->DoDist){
    wbc->conS = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
    wbc->conL = MRIallocSequence(wbc->ntot, 1,1, MRI_FLOAT, wbc->nrholist);
  }

  return(0);
}

int WBCfinish(WBC *wbc)
{
  int nthvox, c,r,s,vtxno, nthrho;
  double val;

  nthvox = 0;
  if(wbc->fvol){
    wbc->volcon = MRIallocSequence(wbc->fvol->width,wbc->fvol->height,wbc->fvol->depth,MRI_FLOAT,wbc->nrholist);
    MRIcopyHeader(wbc->fvol,wbc->volcon);
    MRIcopyPulseParameters(wbc->fvol,wbc->volcon);
    if(wbc->DoDist){
      wbc->volconS = MRIallocSequence(wbc->fvol->width,wbc->fvol->height,wbc->fvol->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->fvol,wbc->volconS);
      MRIcopyPulseParameters(wbc->fvol,wbc->volconS);
      wbc->volconL = MRIallocSequence(wbc->fvol->width,wbc->fvol->height,wbc->fvol->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->fvol,wbc->volconL);
      MRIcopyPulseParameters(wbc->fvol,wbc->volconL);
    }
    for(s=0; s < wbc->fvol->depth; s++){
      for(c=0; c < wbc->fvol->width; c++){
	for(r=0; r < wbc->fvol->height; r++){
	  if(wbc->volmask && MRIgetVoxVal(wbc->volmask,c,r,s,0) < 0.5) continue;
	  for(nthrho=0; nthrho < wbc->nrholist; nthrho++){
	    val = MRIgetVoxVal(wbc->con,nthvox,0,0,nthrho);
	    MRIsetVoxVal(wbc->volcon,c,r,s,nthrho,val);
	    if(wbc->DoDist){
	      val = MRIgetVoxVal(wbc->conS,nthvox,0,0,nthrho);
	      MRIsetVoxVal(wbc->volconS,c,r,s,nthrho,val);
	      val = MRIgetVoxVal(wbc->conL,nthvox,0,0,nthrho);
	      MRIsetVoxVal(wbc->volconL,c,r,s,nthrho,val);
	    }
	  } // rho
	  nthvox ++;
	} // row
      } // col
    } // slice
  }
  if(wbc->flh){
    wbc->lhcon = MRIallocSequence(wbc->flh->width,wbc->flh->height,wbc->flh->depth,MRI_FLOAT,wbc->nrholist);
    MRIcopyHeader(wbc->flh,wbc->lhcon);
    MRIcopyPulseParameters(wbc->flh,wbc->lhcon);
    if(wbc->DoDist){
      wbc->lhconS = MRIallocSequence(wbc->flh->width,wbc->flh->height,wbc->flh->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->flh,wbc->lhconS);
      MRIcopyPulseParameters(wbc->flh,wbc->lhconS);
      wbc->lhconL = MRIallocSequence(wbc->flh->width,wbc->flh->height,wbc->flh->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->flh,wbc->lhconL);
      MRIcopyPulseParameters(wbc->flh,wbc->lhconL);
    }
    for(vtxno = 0; vtxno < wbc->flh->width; vtxno++){
      if(wbc->lhmask && MRIgetVoxVal(wbc->lhmask,vtxno,0,0,0) < 0.5) continue;
      for(nthrho=0; nthrho < wbc->nrholist; nthrho++){
	val = MRIgetVoxVal(wbc->con,nthvox,0,0,nthrho);
	MRIsetVoxVal(wbc->lhcon,vtxno,0,0,nthrho,val);
	if(wbc->DoDist){
	  val = MRIgetVoxVal(wbc->conS,nthvox,0,0,nthrho);
	  MRIsetVoxVal(wbc->lhconS,vtxno,0,0,nthrho,val);
	  val = MRIgetVoxVal(wbc->conL,nthvox,0,0,nthrho);
	  MRIsetVoxVal(wbc->lhconL,vtxno,0,0,nthrho,val);
	}
      } // rho
    } // vertex
  }
  if(wbc->frh){
    wbc->rhcon = MRIallocSequence(wbc->frh->width,wbc->frh->height,wbc->frh->depth,MRI_FLOAT,wbc->nrholist);
    MRIcopyHeader(wbc->frh,wbc->rhcon);
    MRIcopyPulseParameters(wbc->frh,wbc->rhcon);
    if(wbc->DoDist){
      wbc->rhconS = MRIallocSequence(wbc->frh->width,wbc->frh->height,wbc->frh->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->frh,wbc->rhconS);
      MRIcopyPulseParameters(wbc->frh,wbc->rhconS);
      wbc->rhconL = MRIallocSequence(wbc->frh->width,wbc->frh->height,wbc->frh->depth,MRI_FLOAT,wbc->nrholist);
      MRIcopyHeader(wbc->frh,wbc->rhconL);
      MRIcopyPulseParameters(wbc->frh,wbc->rhconL);
    }
    for(vtxno = 0; vtxno < wbc->frh->width; vtxno++){
      if(wbc->rhmask && MRIgetVoxVal(wbc->rhmask,vtxno,0,0,0) < 0.5) continue;
      for(nthrho=0; nthrho < wbc->nrholist; nthrho++){
	val = MRIgetVoxVal(wbc->con,nthvox,0,0,nthrho);
	MRIsetVoxVal(wbc->rhcon,vtxno,0,0,nthrho,val);
	if(wbc->DoDist){
	  val = MRIgetVoxVal(wbc->conS,nthvox,0,0,nthrho);
	  MRIsetVoxVal(wbc->rhconS,vtxno,0,0,nthrho,val);
	  val = MRIgetVoxVal(wbc->conL,nthvox,0,0,nthrho);
	  MRIsetVoxVal(wbc->rhconL,vtxno,0,0,nthrho,val);
	}
      } // rho
    } // vertex
  }
  return(0);
}

/*!
  \fn int Index2UpperSubscript(int N, long i, int *r, int *c)
  \brief Computes the row and col of the ith index in an upper
  triangular matrix ignoring the diagonal and lower triangular
  components. The application for this involves evaluating pairs
  of items in a list excluding self-comparisons (diag) and 
  reverse comparisons (lower). See also Index2UpperSubscriptTest().
  Row, col, and index are all 0-based.
 */
int Index2UpperSubscript(int N, long i, int *r, int *c)
{
  long i1,r1,c1,ir1;
  i1 = i + 1;
  r1 = ceil((-(1-2.0*N) - sqrt( pow((1-2.0*N),2) - 8.0*i1))/2.0);
  ir1 = N*r1 - (r1*(r1+1))/2;
  c1 = N - (ir1-i1);
  *r = r1 - 1;
  *c = c1 - 1;
  return(0);
}

/*!
  \fn int Index2UpperSubscriptTest(int N)
  \brief Test for Index2UpperSubscript()
 */
int Index2UpperSubscriptTest(int N)
{
  int r, c, rt, ct;
  long i,err;

  printf("Index2UpperSubscriptTest(): N = %d\n",N);
  err = 0;
  i = 0;
  for(r=0; r < N-1; r++){
    for(c=r+1; c < N; c++){
      Index2UpperSubscript(N, i, &rt, &ct);
      if(r != rt || c != ct) {
	err++;
	printf("ERROR: %4d %4d %4ld  %4d %4d\n",r,c,i,rt,ct);
      }
      i = i + 1;
    }
  }
  printf("Found %ld errors\n",err);
  if(err) return(1);

  return(0);
}



