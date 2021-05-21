/**
 * @brief Slice-to-volume registation with in-plane bbr
 *
 * Slice-to-volume registation with in-plane bbr
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


/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
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
#include "macros.h"
#include "diag.h"
#include "minc.h"
#include "fio.h"
#include "mrisurf.h"
#include "fsenv.h"
#include "registerio.h"
#include "timer.h"
#include "matrix.h"
#include "matfile.h"
#include "mri2.h"
#include "mrisutils.h"
#include "transform.h"
#include "resample.h"
#include "volcluster.h"
#include "numerics.h"
#include "romp_support.h"

typedef struct {
  char *movfile;
  int slice;
  char *surffile;
  char *initreg;
  char *outreg,*outreginv;
  double ftol,linmintol;
  int nitersmax;
  int optschema;
  double bbrslope, bbrcenter, bbrsign;
  double outdist, indist; // mm
  double params[12];
  int paramset[12];
  int inc,search,search1d,searchnper,search1diters;
  double searchmul;
  char *outsurffile;
} CMDARGS;

typedef struct {
  MRIS *surf;
  MRI *mov, *ref;
  int slice,nhits;
  double outdist, indist; // mm
  double bbrslope, bbrcenter, bbrsign;
  double cost;
  int inc;
  int nparams;
  double params[12],pmin[12],pmax[12],pdelta[12];
  double searchmul;
  int searchnper;
  int nCostEvaluations;
  double tLastEval;
  double ftol,linmintol;
  float fret;
  int nitersmax,niters;
  int startmin;
  int optschema;
  MATRIX *R; // surf->xyz to vol tkreg
  LTA *R0;
  MATRIX *D,*invD,*Surf2Vox,*Q,*invQ,*Q0,*invQ0;
  MATRIX *Vmov,*Vref,*Tmov,*Tref;
  MATRIX *invVref,*invVmov,*invTmov,*invTref;
  int debug;
} SBBR;

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

double SBBRcost(SBBR *sbbr);
double BBRcost(double vctx, double vwm, double slope, 
	       double center, double sign, double *pct);

float SBBRcostPowell(float *pPowel) ;
int SBBRMinPowell();
double *SBBRoptSchema2MatrixPar(SBBR *sbbr, double *par);
int SBBRmatrixPar2OptSchema(SBBR *sbbr, double *mpar);
MATRIX *COREGmatrix(double *p, MATRIX *M);
double *COREGparams9(MATRIX *M9, double *p);
int SBBRnParams(SBBR *sbbr);
int SBBRmatrices(SBBR *sbbr, int init);
int SBBRsearchRange(SBBR *sbbr);
int SBBRbruteSearch(SBBR *sbbr);
int SBBRbruteSearchUpdate(SBBR *sbbr, double p[12], double popt[12], double *mincost);
int SSBmin1D(SBBR *sbbr);
int MRIScomputeFaceNormal(MRIS *surf, int faceno, double snorm[3]);
int NormalizeVect3(double v[3]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *subject, *hemi, *SUBJECTS_DIR;
SBBR *sbbr;
CMDARGS *cmdargs;


/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,n,err;
  LTA *lta,*ltainv;
  double rms;

  //double mpar0[12];
  cmdargs = (CMDARGS *)calloc(sizeof(CMDARGS),1);
  cmdargs->optschema = 1;
  cmdargs->nitersmax = 10;
  cmdargs->ftol = 1e-8;
  cmdargs->linmintol = 1e-8;
  cmdargs->outdist = 2.0;
  cmdargs->indist = 1.0;
  cmdargs->bbrslope = 0.5;
  cmdargs->bbrcenter = 0.0;
  cmdargs->bbrsign = -1; //-1 for T2, and +1 for T1
  cmdargs->inc = 1;
  cmdargs->search = 0;
  cmdargs->search1d = 0;
  cmdargs->searchmul = 1;
  cmdargs->searchnper = 11;
  cmdargs->search1diters = 10;

  nargs = handleVersionOption(argc, argv, "mri_sbbr");
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

  sbbr = (SBBR *) calloc(sizeof(SBBR),1);
  sbbr->slice = cmdargs->slice;
  sbbr->outdist = cmdargs->outdist;
  sbbr->indist = cmdargs->indist;
  sbbr->bbrslope = cmdargs->bbrslope;
  sbbr->bbrcenter =cmdargs->bbrcenter;
  sbbr->bbrsign = cmdargs->bbrsign;
  sbbr->nitersmax = cmdargs->nitersmax;
  sbbr->ftol = cmdargs->ftol;
  sbbr->linmintol = cmdargs->linmintol;
  sbbr->inc = cmdargs->inc;
  sbbr->optschema = cmdargs->optschema;
  sbbr->searchmul = cmdargs->searchmul;
  sbbr->searchnper = cmdargs->searchnper;
  sbbr->debug = debug;
  SBBRnParams(sbbr); // sets number of parameters

  // Initialize parameters appropriately (eg, scale must be 1)
  SBBRmatrixPar2OptSchema(sbbr, NULL); 
  // Set any parameters input from command line
  for(n=0; n < sbbr->nparams; n++) if(cmdargs->paramset[n]) sbbr->params[n] = cmdargs->params[n];
  printf("Init Params: ");
  for(n=0; n < sbbr->nparams; n++) printf(" %6.4f",sbbr->params[n]);
  printf("\n");

  // Read inputs
  sbbr->surf = MRISread(cmdargs->surffile);
  if(sbbr->surf == NULL) exit(1);
  // Create an MRI header for the ref volume the surf was created from
  sbbr->ref = MRIallocFromVolGeom(&(sbbr->surf->vg), MRI_UCHAR, 1, 1);
  sbbr->mov = MRIread(cmdargs->movfile);
  if(sbbr->mov == NULL) exit(1);
  if(cmdargs->initreg){
    printf("Reading %s\n",cmdargs->initreg);
    sbbr->R0 = LTAread(cmdargs->initreg);
    if(sbbr->R0 == NULL) exit(1);
    LTAchangeType(sbbr->R0,LINEAR_RAS_TO_RAS);
    // Need to check whether it goes in the right direction
    //LTAinvert(sbbr->R0,sbbr->R0);
  }

  // Store the xyz coords in origxyz
  MRISsaveVertexPositions(sbbr->surf, ORIGINAL_VERTICES) ;

  // Create the matrices
  SBBRmatrices(sbbr,1);

  // Compute an init cost to make sure things are ok
  SBBRcost(sbbr);
  if(sbbr->nhits == 0){
    printf("ERROR: there is no intersection between the surf and volume\n");
    exit(1);
  }
  printf("init cost = %6.4lf, nhits = %d\n",sbbr->cost,sbbr->nhits);
  SBBRcost(sbbr);

  SBBRsearchRange(sbbr);
  if(cmdargs->search) SBBRbruteSearch(sbbr);
  if(cmdargs->search1d){
    printf("1D Search ----------------------\n");
    for(n=0; n<10; n++) SSBmin1D(sbbr);
    printf("1D Search done ----------------------\n");
  }

  // Optimize
  sbbr->startmin = 1;
  SBBRMinPowell();

  // Done
  SBBRcost(sbbr);
  printf("final cost = %6.4lf, nhits = %d\n",sbbr->cost,sbbr->nhits);
  printf("TotalCostEvals %d\n",sbbr->nCostEvaluations);

  // Create the matrices and LTA
  SBBRmatrices(sbbr,0);
  lta = LTAcreate(sbbr->mov, sbbr->ref, sbbr->R, LINEAR_RAS_TO_RAS);
  lta->fscale = sbbr->R0->fscale;
  strcpy(lta->subject,sbbr->R0->subject);
  err = LTAwrite(lta,cmdargs->outreg);
  if(err) exit(1);

  // Compute rms diff (wont nec agree with mri_coreg because not TKRAS
  rms = RMSregDiffMJ(sbbr->R0->xforms[0].m_L, lta->xforms[0].m_L, 50.0);
  printf("RMS Diff %6.4lf\n",rms);

  printf("tkregisterfv --mov %s --reg %s ...\n",cmdargs->movfile,cmdargs->outreg);

  if(cmdargs->outreginv || cmdargs->outsurffile)  ltainv = LTAinvert(lta,NULL);
  if(cmdargs->outreginv){
    err = LTAwrite(ltainv,cmdargs->outreginv);
    if(err) exit(1);
    printf("tkregisterfv --targ %s --reg %s ...\n",cmdargs->movfile,cmdargs->outreginv);
  }
  if(cmdargs->outsurffile){
    LTAchangeType(lta,REGISTER_DAT);
    MRISrestoreVertexPositions(sbbr->surf, ORIGINAL_VERTICES) ;
    MRISmatrixMultiply(sbbr->surf,lta->xforms[0].m_L);
    if(MatrixDeterminant(lta->xforms[0].m_L) < 0) MRISreverseFaceOrder(sbbr->surf);
    getVolGeom(sbbr->mov, &sbbr->surf->vg);
    MRIScomputeMetricProperties(sbbr->surf);
    MRISwrite(sbbr->surf,cmdargs->outsurffile);
  }

  printf("mri_sbbr done\n");

  return 0;
}
/* end main --------------------------------------------------- */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
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
    else if (!strcasecmp(option, "--t1")) cmdargs->bbrsign = +1;
    else if (!strcasecmp(option, "--t2")) cmdargs->bbrsign = -1;

    else if (!strcasecmp(option, "--mov")) {
      if (nargc < 1) CMDargNErr(option,1);
      cmdargs->movfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--slice")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cmdargs->slice);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--surf")) {
      if (nargc < 1) CMDargNErr(option,1);
      cmdargs->surffile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--out-surf")) {
      if (nargc < 1) CMDargNErr(option,1);
      cmdargs->outsurffile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--init-reg")) {
      if (nargc < 1) CMDargNErr(option,1);
      cmdargs->initreg = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--reg")) {
      if (nargc < 1) CMDargNErr(option,1);
      cmdargs->outreg = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--reg-inv")) {
      if (nargc < 1) CMDargNErr(option,1);
      cmdargs->outreginv = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--opt")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cmdargs->optschema);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--din")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->indist);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--dout")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->outdist);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--slope")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->bbrslope);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ftol")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->ftol);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--inc")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cmdargs->inc);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--search")){
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&cmdargs->searchnper);
      sscanf(pargv[1],"%lf",&cmdargs->searchmul);
      cmdargs->search = 1;
      nargsused = 2;
    }
    else if (!strcasecmp(option, "--search1d")){
      if(nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%d",&cmdargs->search1diters);
      sscanf(pargv[0],"%d",&cmdargs->searchnper);
      sscanf(pargv[1],"%lf",&cmdargs->searchmul);
      cmdargs->search1d = 1;
      nargsused = 3;
    }
    else if (!strcasecmp(option, "--niters-max")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cmdargs->nitersmax);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--linmintol")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->linmintol);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--p")) {
      if(nargc < 2) CMDargNErr(option,2);
      int paramno;
      sscanf(pargv[0],"%d",&paramno);
      sscanf(pargv[1],"%lf",&cmdargs->params[paramno]);
      cmdargs->paramset[paramno] = 1;
      nargsused = 2;
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
    else if (!strcasecmp(option, "--diag")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag-show"))    Gdiag = (Gdiag & DIAG_SHOW);
    else if (!strcasecmp(option, "--diag-verbose")) Gdiag = (Gdiag & DIAG_VERBOSE);
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
/* -------------------------------------------------- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -------------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --mov volfile : template volume \n");
  printf("   --surf surfacefile \n");
  printf("   --init-reg regfile \n");
  printf("   --t1, --t2\n");
  printf("   --opt optno : optimization type 1, 2, or 3, default is 1 (6 dof)\n");
  printf("   --din dinmm : dist in mm into surface %lf\n",cmdargs->indist);
  printf("   --dout doutmm : dist in mm out of surface %lf\n",cmdargs->outdist);
  printf("   --slope bbrslope : default %lf\n",cmdargs->bbrslope);
  printf("   --ftol ftol : default %le\n",cmdargs->ftol);
  printf("   --linmintol linmintol : default %lf\n",cmdargs->linmintol);
  printf("   --niters-max nmax : default %d\n",cmdargs->nitersmax);
  printf("   --search nper mul : brute force through param space\n");
  printf("   --search1d niters nper mul : 1d search through param space\n");
  printf("   --p parno par : set init parameter \n");
  printf("   --inc n : face number increment (default 1)\n");
  printf("   --slice sliceno : defaults to 0\n");
  printf("   --threads nthreads \n");
  printf("   --reg regfile : output registration \n");
  printf("   --reg-inv invregfile : inverted output registration \n");
  printf("   --out-surf surfname : output surface in slice coords\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --diag diagno turn on diagnositcs\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* -------------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  printf("This program is a special implementation of boundary-based registration\n");
  printf("that computes the registration for a single slice.\n");
  exit(1) ;
}
/* -------------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* -------------------------------------------------- */
static void check_options(void) {
  if(cmdargs->movfile == NULL){
    printf("ERROR: must spec --mov\n");
    exit(1);
  }
  if(cmdargs->surffile == NULL){
    printf("ERROR: must spec --surf\n");
    exit(1);
  }
  if(cmdargs->outreg == NULL){
    printf("ERROR: must spec --reg\n");
    exit(1);
  }
  return;
}
/* -------------------------------------------------- */
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

/*-----------------------------------------------------------------*/
double SBBRcost(SBBR *sbbr)
{
  int faceno,nhits,n;
  double cost,sumcost,detSign;
  double coutdistvox, routdistvox, cindistvox, rindistvox;
  MRIS *surf = sbbr->surf;
  MRI *mov = sbbr->mov;
  static double par[12];

  coutdistvox = sbbr->outdist/mov->xsize;
  routdistvox = sbbr->outdist/mov->ysize;
  cindistvox = sbbr->indist/mov->xsize;
  rindistvox = sbbr->indist/mov->ysize;

  SBBRoptSchema2MatrixPar(sbbr, par);
  sbbr->Q = COREGmatrix(par, sbbr->Q);
  sbbr->invQ = MatrixInverse(sbbr->Q,sbbr->invQ);
  // Surf2Vox = inv(D)*inv(Q)*inv(Q0)
  sbbr->Surf2Vox = MatrixMultiplyD(sbbr->invD,sbbr->invQ,sbbr->Surf2Vox);
  sbbr->Surf2Vox = MatrixMultiplyD(sbbr->Surf2Vox,sbbr->invQ0,sbbr->Surf2Vox);
  //printf("Surf2Vox -----------------\n");
  //MatrixPrint(stdout,sbbr->Surf2Vox);

  detSign = SIGN(MatrixDeterminant(sbbr->Surf2Vox));

  MRISrestoreVertexPositions(surf, ORIGINAL_VERTICES) ;
  MRISmatrixMultiply(surf,sbbr->Surf2Vox);

  //FILE *fp;
  //char tmpstr[1000];
  //sprintf(tmpstr,"cost.%03d.dat",nthcall);
  //fp = fopen(tmpstr,"w");
  nhits = 0;
  sumcost = 0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
#pragma omp parallel for if_ROMP(experimental) reduction(+:nhits,sumcost)
#endif
  for(faceno=0; faceno < surf->nfaces; faceno += sbbr->inc) {
    ROMP_PFLB_begin
    
    FACE *face;
    VERTEX *vtx;
    int npos,vtxno,skip,i0=0,i1=0,i2=0,m,oob,n;
    double dz[3],cA,rA,sA,cB,rB,sB;
    double dAB,nABc,nABr,nABs,t,c0[2],r0[2];
    double c0mn,r0mn,d,uc,ur,lc,lr,snc,snr;
    double lval,uval, fcost, pct;
    double snorm[3];

    face = &(surf->faces[faceno]);
    if(face->ripflag) continue;

    // check whether there is an intersection
    npos = 0;
    skip = 0;
    //printf("%d/%d %d %d\n",faceno,surf->nfaces,npos,skip);
    for(n=0; n < 3; n++) dz[n] = 0;    
    for(n=0; n < 3; n++){
      vtxno = face->v[n];
      vtx = &(surf->vertices[vtxno]);
      dz[n] = vtx->z - sbbr->slice;
      if(dz[n] > 0) npos ++;
      // If any vertex is outofbounds, skip this face
      if(vtx->x < 0) skip = 1;
      if(vtx->x >= mov->width) skip = 1;
      if(vtx->y < 0) skip = 1;
      if(vtx->y >= mov->height) skip = 1;
      //if(vtx->z < 0) skip = 1;
      //if(vtx->z >= mov->depth) skip = 1;
    }
    if(npos ==0 || npos == 3 || skip) continue;
    // If it gets here, there must be an intersection and not oob

    // One of the vertices must be on one side of the slice and the other 
    // two must be on the other side. Assign pairs
    if(SIGN(dz[0]) != SIGN(dz[1]) && SIGN(dz[0]) != SIGN(dz[2])){
      i0 = 0;
      i1 = 1;
      i2 = 2;
    }
    if(SIGN(dz[1]) != SIGN(dz[0]) && SIGN(dz[1]) != SIGN(dz[2])){
      i0 = 1;
      i1 = 0;
      i2 = 2;
    }
    if(SIGN(dz[2]) != SIGN(dz[0]) && SIGN(dz[2]) != SIGN(dz[1])){
      i0 = 2;
      i1 = 0;
      i2 = 1;
    }
    //printf("%d/%d %d %d  %d %d %d  ",faceno,surf->nfaces,npos,skip,i0,i1,i2);

    vtxno = surf->faces[faceno].v[i0];
    vtx = &(surf->vertices[vtxno]);
    cA = vtx->x;
    rA = vtx->y;
    sA = vtx->z-sbbr->slice;
    for(m = 0; m < 2; m++){
      if(m==0) vtxno = surf->faces[faceno].v[i1];
      if(m==1) vtxno = surf->faces[faceno].v[i2];
      vtx = &(surf->vertices[vtxno]);
      cB = vtx->x;
      rB = vtx->y;
      sB = vtx->z-sbbr->slice;

      dAB = sqrt((cA-cB)*(cA-cB) + (rA-rB)*(rA-rB) + (sA-sB)*(sA-sB));
      nABc = (cA-cB)/dAB;
      nABr = (rA-rB)/dAB;
      nABs = (sA-sB)/dAB;

      t = -sA/nABs;
      c0[m] = cA + t * nABc;
      r0[m] = rA + t * nABr;
    }

    //compute midpoint between to intersections
    c0mn = (c0[0]+c0[1])/2;
    r0mn = (r0[0]+r0[1])/2;

    // Compute normal for this face
    MRIScomputeFaceNormal(surf, faceno, snorm);
    d = sqrt(snorm[0]*snorm[0] + snorm[1]*snorm[1]);
    snc = snorm[0]/d;
    snr = snorm[1]/d;

    // sample outside the surface (along the normal)
    uc = c0mn + detSign*coutdistvox*snc;
    ur = r0mn + detSign*routdistvox*snr;
    oob = MRIindexNotInVolume(mov,uc,ur,sbbr->slice);
    if(oob) continue;
    uval = MRIgetVoxVal(mov,nint(uc),nint(ur),sbbr->slice,0);

    // sample inside the surface (along the normal)
    lc = c0mn - detSign*coutdistvox*snc;
    lr = r0mn - detSign*routdistvox*snr;
    oob = MRIindexNotInVolume(mov,lc,lr,sbbr->slice);
    if(oob) continue;
    lval = MRIgetVoxVal(mov,nint(lc),nint(lr),sbbr->slice,0);

    nhits++;
    fcost = BBRcost(uval, lval, sbbr->bbrslope, sbbr->bbrcenter,sbbr->bbrsign,&pct);
    sumcost += fcost;

    //printf(   "%7d   %7.4lf %7.4lf %7.4f    %7.4lf %7.4lf %7.4f   %7.4f\n",
    //	      faceno,uc,ur,uval,lc,lr,lval,fcost);
     //    fprintf(fp,"%7d   %7.4lf %7.4lf %7.4f    %7.4lf %7.4lf %7.4f   %7.4f   %6.4f %6.4f %6.4f %6.4f  %6.4f %6.4f  %6.4f %6.4f  %6.4f\n",
    //	    faceno,uc,ur,uval,lc,lr,lval,0.0,c0[0],r0[0],c0[1],r0[1],c0mn,r0mn,snc,snr,fcost);
    ROMP_PFLB_end
  } // end face loop
  ROMP_PF_end
  
  //fclose(fp);

  sbbr->nhits = nhits;
  if(sbbr->nhits > 0) cost = sumcost/sbbr->nhits;
  else                cost = 10e10;
  sbbr->cost = cost;

  if(sbbr->debug){
    printf("I  %5d %5d %6.4f    ",sbbr->nCostEvaluations++,sbbr->nhits,cost);
    for(n=0; n < sbbr->nparams; n++) printf(" %6.4f",sbbr->params[n]);
    printf("\n");
    fflush(stdout);
  }

  sbbr->nCostEvaluations++;
  return(cost);
}

double BBRcost(double vctx, double vwm, double slope, 
	       double center, double sign, double *pct)
{
  double d,a=0,c,s;
  s = (vctx+vwm);
  if(fabs(s) > .000001)  d = 100*(vctx-vwm)/(s/2.0); // percent contrast
  else d = 0;
  if(sign ==  0) a = -fabs(slope*(d-center)); // not sure this is useful
  if(sign == -1) a = -(slope*(d-center));
  if(sign == +1) a = +(slope*(d-center));
  if(sign == -2){
    if(d >= 0) a = -(slope*(d-center));
    else       a = 0;
  }
  c = 1+tanh(a);
  *pct = d;
  return(c);
}


/*--------------------------------------------------------------------------*/
float SBBRcostPowell(float *pPowel) 
{
  extern SBBR *sbbr;
  int n,newmin;
  float curcost;
  static float initcost=-1,mincost=-1,ppmin[100];
  FILE *fp;

  for(n=0; n < sbbr->nparams; n++) sbbr->params[n] = pPowel[n+1];
  
  // compute cost
  curcost = SBBRcost(sbbr);

  newmin = 0;
  if(sbbr->startmin) {
    newmin = 1;
    initcost = curcost;
    mincost = curcost;
    for(n=0; n<sbbr->nparams; n++) ppmin[n] = sbbr->params[n];
    printf("InitialCost %20.10lf \n",initcost);
    sbbr->startmin = 0;
  }

  if(mincost > curcost) {
    newmin = 1;
    mincost = curcost;
    for(n=0; n<sbbr->nparams; n++) ppmin[n] = sbbr->params[n];
  }

  if(0){
    fp = stdout;
    fprintf(fp,"%4d  ",sbbr->nCostEvaluations);
    for(n=0; n<sbbr->nparams; n++) fprintf(fp,"%12.8f ",sbbr->params[n]);
    fprintf(fp,"  %12.10f %12.5f\n",curcost/initcost,curcost);
    fflush(fp);
  }

  if(newmin){
    printf("#@# %4d  ",sbbr->nCostEvaluations);
    for(n=0; n<sbbr->nparams; n++) printf("%7.5f ",ppmin[n]);
    printf("  %9.7f\n",mincost);
    fflush(stdout);
  }

  return((float)curcost);
}

/*---------------------------------------------------------*/
int SBBRMinPowell()
{
  extern SBBR *sbbr;
  float *pPowel, **xi;
  int    r, c, n,dof;
  Timer timer;

  timer.reset();
  dof = sbbr->nparams;

  printf("\n\n---------------------------------\n");
  printf("Init Powel Params dof = %d\n",dof);
  pPowel = vector(1, dof) ;
  for(n=0; n < dof; n++) pPowel[n+1] = sbbr->params[n];

  xi = matrix(1, dof, 1, dof) ;
  for (r = 1 ; r <= dof ; r++) {
    for (c = 1 ; c <= dof ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }
  OpenPowell2(pPowel, xi, dof, sbbr->ftol, sbbr->linmintol, sbbr->nitersmax, 
	      &sbbr->niters, &sbbr->fret, SBBRcostPowell);
  printf("Powell done niters total = %d\n",sbbr->niters);
  printf("PowellOptTimeSec %4.1f sec\n", timer.seconds());
  printf("PowellOptTimeMin %5.2f min\n", timer.minutes());
  fflush(stdout);

  printf("Final parameters ");
  for(n=0; n < sbbr->nparams; n++){
    sbbr->params[n] = pPowel[n+1];
    printf("%12.8f ",sbbr->params[n]);
  }
  printf("\n");

  SBBRcost(sbbr);
  printf("Final cost %20.15lf\n ",sbbr->cost);

  free_matrix(xi, 1, dof, 1, dof);
  free_vector(pPowel, 1, dof);
  printf("\n\n---------------------------------\n");
  return(NO_ERROR) ;
}

int SBBRmatrixPar2OptSchema(SBBR *sbbr, double *mpar)
{
  int n;
  switch(sbbr->optschema){
  case 1: 
    // rigid
    if(mpar) for(n=0; n<6; n++) sbbr->params[n] = mpar[n];
    else for(n=0; n<6; n++) sbbr->params[n] = 0.0;
    break;
  case 2: 
    // rigid + in-plane scaling
    if(mpar) for(n=0; n<8; n++) sbbr->params[n] = mpar[n];
    else{
      for(n=0; n<6; n++) sbbr->params[n] = 0.0;
      sbbr->params[6] = 1.0;
      sbbr->params[7] = 1.0;
    }
    break;
  case 3: 
    // rigid + in-plane scaling + in-plane shear
    if(mpar) {
      for(n=0; n<8; n++) sbbr->params[n] = mpar[n];
      sbbr->params[8] = mpar[9];
    }
    else {
      for(n=0; n<6; n++) sbbr->params[n] = 0.0;
      sbbr->params[6] = 1.0;
      sbbr->params[7] = 1.0;
      sbbr->params[8] = 0.0;
    }
    break;
  }
  return(0);
}

/*----------------------------------------------------------------*/
double *SBBRoptSchema2MatrixPar(SBBR *sbbr, double *mpar)
{
  int n;

  if(mpar == NULL) mpar = (double *) calloc(12,sizeof(double));

  // proper init, scaling must be 1
  for(n=0; n<12; n++) mpar[n] = 0;
  mpar[6] = mpar[7] = mpar[8] = 1; // scaling

  switch(sbbr->optschema){
  case 1: 
    // rigid
    for(n=0; n < 6; n++) mpar[n] = sbbr->params[n];
    sbbr->nparams = 6;
    break;
  case 2: 
    // rigid + in-plane scaling
    for(n=0; n < 8; n++) mpar[n] = sbbr->params[n];
    sbbr->nparams = 8; // excludes thru-plane scale and all shears
    break;
  case 3: 
    // rigid + in-plane scaling + in-plane shear
    for(n=0; n < 8; n++) mpar[n] = sbbr->params[n];
    mpar[9] = sbbr->params[8];
    sbbr->nparams = 9; // excludes thru-plane scale and shears
    break;
  }
  return(mpar);
}
/*----------------------------------------------------------------*/
int SBBRnParams(SBBR *sbbr)
{
  switch(sbbr->optschema){
  case 1: 
    // rigid
    sbbr->nparams = 6;
    break;
  case 2: 
    // rigid + in-plane scaling
    sbbr->nparams = 8; // excludes thru-plane scale and all shears
    break;
  case 3: 
    // rigid + in-plane scaling + in-plane shear
    sbbr->nparams = 10; // excludes thru-plane scale and shear
    break;
  }
  return(sbbr->nparams);
}
/*----------------------------------------------------------------*/
int SBBRsearchRange(SBBR *sbbr)
{
  int n;
  switch(sbbr->optschema){
  case 1: 
    // rigid
    for(n=0; n < 6; n++){
      sbbr->pmin[n] = -4.0*sbbr->searchmul;
      sbbr->pmax[n] = +4.0*sbbr->searchmul;
    }
    break;
  case 2: 
    // rigid + in-plane scaling
    for(n=0; n < 6; n++){
      sbbr->pmin[n] = -4*sbbr->searchmul;
      sbbr->pmax[n] = +4*sbbr->searchmul;
    }
    for(n=6; n < 8; n++){
      sbbr->pmin[n] = 1 - 0.3*sbbr->searchmul;
      sbbr->pmax[n] = 1 + 0.3*sbbr->searchmul;
    }
    break;
  case 3: 
    // rigid + in-plane scaling + in-plane shear
    for(n=0; n < 6; n++){
      sbbr->pmin[n] = -4*sbbr->searchmul;
      sbbr->pmax[n] = +4*sbbr->searchmul;
    }
    for(n=6; n < 8; n++){
      sbbr->pmin[n] = 1 - 0.3*sbbr->searchmul;
      sbbr->pmax[n] = 1 + 0.3*sbbr->searchmul;
    }
    for(n=8; n < 10; n++){
      sbbr->pmin[n] = -0.2*sbbr->searchmul;
      sbbr->pmax[n] = +0.2*sbbr->searchmul;
    }
    break;
  }

  for(n=0; n < sbbr->nparams; n++)
    sbbr->pdelta[n] = (sbbr->pmax[n]-sbbr->pmin[n])/(sbbr->searchnper-1);
  printf("Search ranges\n");
  for(n = 0; n < sbbr->nparams; n++)  
    printf("  %2d %7.4f %7.4f %7.4f \n",n,sbbr->pmin[n],sbbr->pmax[n],sbbr->pdelta[n]);

  return(0);
}

/*!
  \fn MATRIX *COREGmatrix(double *p, MATRIX *M)
  \brief Computes a RAS-to-RAS transformation matrix given
  the parameters. p[0-2] translation, p[3-5] rotation
  in degrees, p[6-8] scale, p[9-11] shear. 
  M = T*R1*R2*R3*SCALE*SHEAR
  Consistent with COREGparams9()
 */
MATRIX *COREGmatrix(double *p, MATRIX *M)
{
  MATRIX *T, *R1, *R2, *R3, *R, *ZZ, *S;
  //int n;
  //printf("p = [");
  //for(n=0; n<np; n++) printf("%10.3lf ",p[n]);
  //printf("];\n");

  // translations
  T = MatrixIdentity(4,NULL);
  T->rptr[1][4] = p[0];
  T->rptr[2][4] = p[1];
  T->rptr[3][4] = p[2];

  // rotations
  R1 = MatrixIdentity(4,NULL);
  R1->rptr[2][2] = cos(p[3]*M_PI/180);
  R1->rptr[2][3] = sin(p[3]*M_PI/180);
  R1->rptr[3][2] = -sin(p[3]*M_PI/180);
  R1->rptr[3][3] = cos(p[3]*M_PI/180);

  R2 = MatrixIdentity(4,NULL);
  R2->rptr[1][1] = cos(p[4]*M_PI/180);
  R2->rptr[1][3] = sin(p[4]*M_PI/180);
  R2->rptr[3][1] = -sin(p[4]*M_PI/180);
  R2->rptr[3][3] = cos(p[4]*M_PI/180);

  R3 = MatrixIdentity(4,NULL);
  R3->rptr[1][1] = cos(p[5]*M_PI/180);
  R3->rptr[1][2] = sin(p[5]*M_PI/180);
  R3->rptr[2][1] = -sin(p[5]*M_PI/180);
  R3->rptr[2][2] = cos(p[5]*M_PI/180);

  R = MatrixMultiplyD(R1,R2,NULL);
  MatrixMultiplyD(R,R3,R);

  // scale, use ZZ because some idiot #defined Z
  ZZ = MatrixIdentity(4,NULL);
  ZZ->rptr[1][1] = p[6];
  ZZ->rptr[2][2] = p[7];
  ZZ->rptr[3][3] = p[8];
  ZZ->rptr[4][4] = 1;

  // shear
  S = MatrixIdentity(4,NULL);
  S->rptr[1][2] = p[9];
  S->rptr[1][3] = p[10];
  S->rptr[2][3] = p[11];

  // M = T*R*ZZ*S
  M = MatrixMultiplyD(T,R,M);
  MatrixMultiplyD(M,ZZ,M);
  MatrixMultiplyD(M,S,M);
  //MatrixPrint(stdout,M);

  MatrixFree(&T);
  MatrixFree(&R1);
  MatrixFree(&R2);
  MatrixFree(&R3);
  MatrixFree(&R);
  MatrixFree(&ZZ);
  MatrixFree(&S);

  return(M);
}


/*!
  \fn double *COREGparams9(MATRIX *M9, double *p)
  \brief Extracts parameter from a 9 dof transformation matrix.
  This is consistent with COREGmatrix(). Still need to figure
  out how to do 12 dof. Note: p will be alloced to 12.
  Angles are in degrees.
 */
double *COREGparams9(MATRIX *M9, double *p)
{
  double sum;
  int c, r;
  MATRIX *R;

  if(p==NULL) p = (double*) calloc(12,sizeof(double));

  // translation
  for(r=0; r < 3; r++) p[r] = M9->rptr[r+1][4];

  // R is the rotation matrix
  R = MatrixAlloc(3,3,MATRIX_REAL);
  for(c=0; c < 3; c++){
    sum = 0;
    for(r=0; r < 3; r++) sum += (M9->rptr[r+1][c+1]*M9->rptr[r+1][c+1]);
    p[c+6] = sqrt(sum); //scale
    for(r=0; r < 3; r++) R->rptr[r+1][c+1] = M9->rptr[r+1][c+1]/sqrt(sum);
  }

  // extract rotation params
  p[3] = atan2(R->rptr[2][3],R->rptr[3][3])*180/M_PI;
  p[4] = atan2(R->rptr[1][3],sqrt(pow(R->rptr[2][3],2) + pow(R->rptr[3][3],2)))*180/M_PI;
  p[5] = atan2(R->rptr[1][2],R->rptr[1][1])*180/M_PI;

  MatrixFree(&R);

  //for(n=0; n < 9; n++) printf("%10.8lf ",p[n]);
  //printf("\n");

  return(p);
}

int SBBRmatrices(SBBR *sbbr, int init)
{
  double mpar[12];
  int nthp;

  sbbr->D = MatrixAlloc(4,4,MATRIX_REAL);
  sbbr->D->rptr[1][1] = sbbr->mov->xsize;
  sbbr->D->rptr[2][2] = sbbr->mov->ysize;
  sbbr->D->rptr[3][3] = sbbr->mov->zsize;
  sbbr->D->rptr[1][4] = -sbbr->mov->xsize*sbbr->mov->width/2.0;
  sbbr->D->rptr[2][4] = -sbbr->mov->ysize*sbbr->mov->height/2.0;
  sbbr->D->rptr[4][4] = 1;

  sbbr->invD = MatrixInverse(sbbr->D,NULL);
  sbbr->Vmov = MRIxfmCRS2XYZ(sbbr->mov,0);
  sbbr->invVmov = MatrixInverse(sbbr->Vmov,NULL);
  sbbr->Vref = MRIxfmCRS2XYZ(sbbr->ref,0);
  sbbr->invVref = MatrixInverse(sbbr->Vref,NULL);
  sbbr->Tmov = MRIxfmCRS2XYZtkreg(sbbr->mov);
  sbbr->Tref = MRIxfmCRS2XYZtkreg(sbbr->ref);
  sbbr->invTref = MatrixInverse(sbbr->Tref,NULL);

  if(sbbr->R0 == NULL)
    sbbr->R0 = LTAcreate(sbbr->mov, sbbr->ref, NULL, LINEAR_RAS_TO_RAS);

  /* When making Q with init=1, ignore any parameters set on the cmd
     line. This only makes a difference for Q0. In the end, it won't
     make a difference in the final reg matrix. However, the
     parameters that it finds will be relative to the init matrix
     without considering the param set on the cmd line.*/
  if(! init){
    SBBRoptSchema2MatrixPar(sbbr, mpar);
    sbbr->Q = COREGmatrix(mpar, sbbr->Q);
  }
  else sbbr->Q = MatrixIdentity(4,NULL);
  sbbr->invQ = MatrixInverse(sbbr->Q,sbbr->invQ);

  if(init){
    // Q0 = Tref * inv(Vref) * R *Vmov *inv(D) * inv(Q)
    sbbr->Q0 = MatrixMultiplyD(          sbbr->Tref,sbbr->invVref, sbbr->Q0);
    sbbr->Q0 = MatrixMultiplyD(sbbr->Q0, sbbr->R0->xforms[0].m_L, sbbr->Q0);
    sbbr->Q0 = MatrixMultiplyD(sbbr->Q0, sbbr->Vmov, sbbr->Q0);
    sbbr->Q0 = MatrixMultiplyD(sbbr->Q0, sbbr->invD, sbbr->Q0);
    sbbr->Q0 = MatrixMultiplyD(sbbr->Q0, sbbr->invQ, sbbr->Q0);
    sbbr->invQ0 = MatrixInverse(sbbr->Q0,sbbr->invQ0);
  }

  // Surf2Vox = inv(D)*inv(Q)*inv(Q0)
  sbbr->Surf2Vox = MatrixMultiplyD(sbbr->invD,sbbr->invQ,sbbr->Surf2Vox);
  sbbr->Surf2Vox = MatrixMultiplyD(sbbr->Surf2Vox,sbbr->invQ0,sbbr->Surf2Vox);

  // R = Vref * inv(Tref) * Q0 * Q * D * inv(Vmov)
  sbbr->R = MatrixMultiplyD(          sbbr->Vref,sbbr->invTref, sbbr->R);
  sbbr->R = MatrixMultiplyD(sbbr->R, sbbr->Q0, sbbr->R);
  sbbr->R = MatrixMultiplyD(sbbr->R, sbbr->Q, sbbr->R);
  sbbr->R = MatrixMultiplyD(sbbr->R, sbbr->D, sbbr->R);
  sbbr->R = MatrixMultiplyD(sbbr->R, sbbr->invVmov, sbbr->R);

  if(sbbr->debug){
    printf("-------------------------------------------\n");
    printf("Params: ");
    for(nthp = 0; nthp < sbbr->nparams; nthp++) printf("%6.4f ",sbbr->params[nthp]);
    printf("\n");

    printf("Q0 -----------------------------\n");
    MatrixPrint(stdout,sbbr->Q0);
    printf("invQ0 -----------------------------\n");
    MatrixPrint(stdout,sbbr->invQ0);
    printf("Q -----------------------------\n");
    MatrixPrint(stdout,sbbr->Q);
    printf("D -----------------------------\n");
    MatrixPrint(stdout,sbbr->D);
    printf("invD -----------------------------\n");
    MatrixPrint(stdout,sbbr->invD);
    printf("Tref -----------------------------\n");
    MatrixPrint(stdout,sbbr->Tref);
    printf("Tmov -----------------------------\n");
    MatrixPrint(stdout,sbbr->Tmov);
    printf("Vref -----------------------------\n");
    MatrixPrint(stdout,sbbr->Vref);
    printf("Vmov -----------------------------\n");
    MatrixPrint(stdout,sbbr->Vmov);
    printf("Surf2Vox -----------------\n");
    MatrixPrint(stdout,sbbr->Surf2Vox);
    printf("R0 -----------------\n");
    MatrixPrint(stdout,sbbr->R0->xforms[0].m_L);
    printf("R -----------------\n");
    MatrixPrint(stdout,sbbr->R);
  }

  return(0);
}

int NormalizeVect3(double v[3]) 
{
  float d;
  d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (d>0) {
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
  }
  return(0);
}


int MRIScomputeFaceNormal(MRIS *surf, int faceno, double snorm[3])
{
  VERTEX  *v0, *v1, *v2 ;
  double v10[3];
  double v20[3];
  double v21[3];
  FACE *face;

  face = &(surf->faces[faceno]);

  v0 = &(surf->vertices[face->v[0]]) ; //v
  v1 = &(surf->vertices[face->v[1]]) ; //vn1
  v2 = &(surf->vertices[face->v[2]]) ; //vn0

  // v0 = v - vn0 = v0 - v2
  v10[0] = v0->x - v2->x;
  v10[1] = v0->y - v2->y;
  v10[2] = v0->z - v2->z;
  NormalizeVect3(v10);

  // v1 = vn1 - v = v1 - v0
  v20[0] = v1->x - v0->x;
  v20[1] = v1->y - v0->y;
  v20[2] = v1->z - v0->z;
  NormalizeVect3(v20);

  v21[0] = v1->x - v2->x;
  v21[1] = v1->y - v2->y;
  v21[2] = v1->z - v2->z;
  NormalizeVect3(v21);

  snorm[0] = -v20[1]*v10[2] + v10[1]*v20[2];
  snorm[1] =  v20[0]*v10[2] - v10[0]*v20[2];
  snorm[2] = -v20[0]*v10[1] + v10[0]*v20[1];
  NormalizeVect3(snorm);

  if(0){
  int n;
  printf("v10 = [");
  for(n=0; n<3; n++) printf("%lf ",v10[n]);
  printf("]';\n");
  printf("v20 = [");
  for(n=0; n<3; n++) printf("%lf ",v20[n]);
  printf("]';\n");
  printf("v21 = [");
  for(n=0; n<3; n++) printf("%lf ",v21[n]);
  printf("]';\n");
  printf("snorm = [");
  for(n=0; n<3; n++) printf("%lf ",snorm[n]);
  printf("]';\n");
  }

  return(0);
}

int SBBRbruteSearch(SBBR *sbbr)
{
  double p[12],mincost,popt[12];
  long totiter,iter;
  int newmin,nthp;
  Timer timer;

  timer.reset();

  totiter = pow(sbbr->searchnper,sbbr->nparams);
  printf("Starting search over %ld iterations (%d,%d)\n",totiter,sbbr->nparams,sbbr->searchnper);

  mincost = SBBRcost(sbbr);
  printf("Starting cost %lf\n",mincost);
  iter = 0;
  for(p[0] = sbbr->pmin[0]; p[0] <= sbbr->pmax[0]; p[0] += sbbr->pdelta[0]){
    for(p[1] = sbbr->pmin[1]; p[1] <= sbbr->pmax[1]; p[1] += sbbr->pdelta[1]){
      for(p[2] = sbbr->pmin[2]; p[2] <= sbbr->pmax[2]; p[2] += sbbr->pdelta[2]){
	for(p[3] = sbbr->pmin[3]; p[3] <= sbbr->pmax[3]; p[3] += sbbr->pdelta[3]){
	  for(p[4] = sbbr->pmin[4]; p[4] <= sbbr->pmax[4]; p[4] += sbbr->pdelta[4]){
	    for(p[5] = sbbr->pmin[5]; p[5] <= sbbr->pmax[5]; p[5] += sbbr->pdelta[5]){

	      if(sbbr->nparams == 6){
		newmin = SBBRbruteSearchUpdate(sbbr, p, popt, &mincost);
		if(newmin || iter == 0 || iter%1000 == 0){
		  printf("%6ld %6.7f ",iter,mincost);
		  for(nthp = 0; nthp < sbbr->nparams; nthp++)  printf("%7.4f ",popt[nthp]);
		  // Time (min) remaining
		  printf("%7.2f \n",(((totiter-iter)*timer.milliseconds()/(iter+1))/1000.0)/60.0);
		  fflush(stdout);
		}
		iter++;
		continue;
	      }

	      for(p[6] = sbbr->pmin[6]; p[6] <= sbbr->pmax[6]; p[6] += sbbr->pdelta[6]){
		for(p[7] = sbbr->pmin[7]; p[7] <= sbbr->pmax[7]; p[7] += sbbr->pdelta[7]){

		  if(sbbr->nparams == 8){
		    newmin = SBBRbruteSearchUpdate(sbbr, p, popt, &mincost);
		    if(newmin || iter == 0){
		      printf("%6ld %6.7f ",iter,mincost);
		      for(nthp = 0; nthp < sbbr->nparams; nthp++)  printf("%7.4f ",popt[nthp]);
		      printf("%7.2f \n",(((totiter-iter)*timer.milliseconds()/(iter+1))/1000.0)/60.0);
		      fflush(stdout);
		    }
		    iter++;
		    continue;
		  }

		  for(p[8] = sbbr->pmin[8]; p[8] <= sbbr->pmax[8]; p[8] += sbbr->pdelta[8]){
		    for(p[9] = sbbr->pmin[9]; p[9] <= sbbr->pmax[9]; p[9] += sbbr->pdelta[9]){

		      if(sbbr->nparams == 10){
			newmin = SBBRbruteSearchUpdate(sbbr, p, popt, &mincost);
			if(newmin || iter == 0){
			  printf("%6ld %6.7f ",iter,mincost);
			  for(nthp = 0; nthp < sbbr->nparams; nthp++)  printf("%7.4f ",popt[nthp]);
			  printf("%7.2f \n",(((totiter-iter)*timer.milliseconds()/(iter+1))/1000.0)/60.0);
			  fflush(stdout);
			}
			iter++;
			continue;
		      }
		      
		    } // p9
		  } // p8

		} // p7
	      } // p6

	    } // p5
	  } // p4
	} // p3
      } // p2
    } // p1
  } // p0

  // update params
  for(nthp = 0; nthp < sbbr->nparams; nthp++)  sbbr->params[nthp] = popt[nthp];

  printf("SearchTimeSec %4.1f sec\n", timer.seconds());
  printf("SearchTimeMin %5.2f min\n", timer.minutes());
  printf("SearchParams %6ld %6.7f ",iter,mincost);
  for(nthp = 0; nthp < sbbr->nparams; nthp++)  printf("%7.4f ",popt[nthp]);
  printf("\n");

  return(iter);
}

int SBBRbruteSearchUpdate(SBBR *sbbr, double p[12], double popt[12], double *mincost)
{
  int nthp,newmin=0;

  for(nthp = 0; nthp < sbbr->nparams; nthp++)  sbbr->params[nthp] = p[nthp];
  SBBRcost(sbbr);
  if(*mincost > sbbr->cost){
    *mincost = sbbr->cost;
    for(nthp = 0; nthp < sbbr->nparams; nthp++)  popt[nthp] = p[nthp];
    newmin=1;
  }

  return(newmin);
}

int SSBmin1D(SBBR *sbbr)
{
  int nthp,newmin,n;
  double p,popt[12],mincost;

  mincost = sbbr->cost;
  for(n = 0; n < sbbr->nparams; n++)  popt[n] = sbbr->params[n];

  for(nthp = 0; nthp < sbbr->nparams; nthp++){
    for(p = sbbr->pmin[nthp]; p <= sbbr->pmax[nthp]; p += sbbr->pdelta[nthp]){
      sbbr->params[nthp] = p;
      newmin = 0;
      SBBRcost(sbbr);
      if(mincost > sbbr->cost){
	mincost = sbbr->cost;
	popt[nthp] = p;
	newmin=1;
      }
      if(newmin || sbbr->debug){
	printf("1D %6d %8.4f %6.7f %6.7f ",nthp,p,mincost,sbbr->cost);
	for(n = 0; n < sbbr->nparams; n++)  printf("%8.4f ",popt[n]);
	printf("\n");
	fflush(stdout);
      }
    }
    // update params
    sbbr->params[nthp] = popt[nthp];
  }
  sbbr->cost = mincost;
  return(0);
}

