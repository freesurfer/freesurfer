/*
  Name:    mri_surfcluster.c
  Author:  Douglas N. Greve 
  email:   analysis-bugs@nmr.mgh.harvard.edu
  Date:    2/27/02
  Purpose: Finds clusters on the surface.
  $Id: mri_surfcluster.c,v 1.19 2005/12/06 21:47:13 greve Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <sys/utsname.h>

#include "error.h"
#include "diag.h"
#include "proto.h"

#include "matrix.h"
#include "mri.h"
#include "mri2.h"
#include "mri_identify.h"
#include "mrisurf.h"
#include "MRIio_old.h"
#include "fio.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "transform.h"
#include "version.h"

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

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_surfcluster.c,v 1.19 2005/12/06 21:47:13 greve Exp $";
char *Progname = NULL;

char *subjectdir = NULL;
char *hemi = NULL;

char *srcid = NULL;
char *srcfmt  = "paint";
int   srcfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *srcsurfid = "white";
char *srcsubjid = NULL;
int   srcframe = 0;
float thmin = -1, thmax = -1;
char *thsign = NULL;
int   thsignid = 0; 
float minarea = -1;

char *maskid   = NULL;
char *maskfmt  = NULL;
int   maskfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *masksubjid  = NULL;
char *masksign  = NULL;
float maskthresh = -1;

char *labelfile = NULL;
int  nth = -1;
char *labelbase = NULL;
char *labelsubjid  = NULL;

char *omaskid = NULL;
char *omaskfmt = "paint";
char *omasksubjid = NULL;

// Constraining label
char  *clabelfile=NULL;
LABEL *clabel=NULL;
int   clabelinv = 0;


char *outid = NULL;
char *outfmt = "paint";
int   outfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *ocnid = NULL;
char *ocnfmt = "paint";
int   ocnfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *ocpvalid = NULL;
char *ocpvalfmt = "paint";
int   ocpvalfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *sumfile  = NULL;

char *outlabelbase = NULL;
char outlabelfile[1000];
int  nthlab = 0;
LABEL *outlabel;
VERTEX *v;

char *synthfunc  = NULL;
char *subjectsdir = NULL;
int debug = 0;

MRI *srcval, *mritmp,  *cnsurf;
MRI_SURFACE *srcsurf;
int  reshape = 1;
int  reshapefactor;

char *xfmfile = "talairach.xfm";
char xfmpath[2000];
MATRIX *XFM;
int FixMNI = 1;

SURFCLUSTERSUM *scs;

int overallmaxvtx,overallminvtx;
float overallmax,overallmin;

CSD *csd=NULL;
char *csdfile;
double pvalLow, pvalHi, ciPct=90, pval, ClusterSize;

MRI *merged;

/*---------------------------------------------------------------*/
int main(int argc, char **argv)
{
  char fname[2000];
  int  n,NClusters,vtx, rt;
  FILE *fp;
  float totarea;
  int nargs;
  struct utsname uts;
  char *cmdline, cwd[2000];

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_surfcluster.c,v 1.19 2005/12/06 21:47:13 greve Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  dump_options(stdout);

  if(clabelfile){
    printf("Loading clabel %s.\n",clabelfile);
    clabel = LabelRead(NULL, clabelfile);
    if(clabel == NULL){
      fprintf(stderr,"ERROR reading %s\n",clabelfile);
      exit(1);
    }
    printf("Found %d points in clabel.\n",clabel->n_points);
  }

  sprintf(xfmpath,"%s/%s/mri/transforms/%s",subjectsdir,srcsubjid,xfmfile);
  //XFM = LoadxfmMatrix(xfmpath);
  XFM = DevolveXFM(srcsubjid, NULL, xfmfile);
  if(XFM == NULL) exit(1);

  printf("------------- XFM matrix (RAS2RAS) ---------------\n");
  printf("%s\n",xfmpath);
  MatrixPrint(stdout,XFM);
  printf("----------------------------------------------------\n");

  sprintf(fname,"%s/%s/surf/%s.%s",subjectsdir,srcsubjid,hemi,srcsurfid);
  printf("Reading source surface %s\n",fname);
  srcsurf = MRISread(fname) ;
  if (!srcsurf)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;
  printf("Done reading source surface\n");

  printf("Computing metric properties\n");
  MRIScomputeMetricProperties(srcsurf);

  /* ------------------ load the source data ----------------------------*/
  printf("Loading source values\n");
  if(!strcmp(srcfmt,"curv")){ /* curvature file */
    sprintf(fname,"%s/%s/surf/%s.%s",subjectsdir,srcsubjid,hemi,srcid);
    printf("Reading curvature file %s\n",fname);
    rt = MRISreadCurvatureFile(srcsurf, fname);
    if(rt != NO_ERROR) exit(1);
    srcval = MRIallocSequence(srcsurf->nvertices, 1, 1,MRI_FLOAT,1);
    for(vtx = 0; vtx < srcsurf->nvertices; vtx++)
      srcsurf->vertices[vtx].val = srcsurf->vertices[vtx].curv;
  }
  else if(!strcmp(srcfmt,"paint") || !strcmp(srcfmt,"w")){
    rt = MRISreadValues(srcsurf,srcid);
    if(rt != NO_ERROR) exit(1);
    srcval = MRIallocSequence(srcsurf->nvertices, 1, 1,MRI_FLOAT,1);
  }
  else { 
    /* ----------- Use MRIreadType -------------------*/
    if(srcfmt == NULL){
      srcval =  MRIread(srcid);
      if(srcval == NULL){
	printf("ERROR: could not read %s\n",srcid);
	exit(1);
      }
    }
    else{
      srcfmtid = string_to_type(srcfmt); 
      srcval =  MRIreadType(srcid,srcfmtid);
      if(srcval == NULL){
	printf("ERROR: could not read %s as type %s\n",srcid,srcfmt);
	exit(1);
      }
    }
    if(srcval->height != 1 || srcval->depth != 1){
      reshapefactor = srcval->height * srcval->depth;
      printf("Reshaping %d\n",reshapefactor);
      mritmp = mri_reshape(srcval, reshapefactor*srcval->width, 
         1, 1, srcval->nframes);
      MRIfree(&srcval);
      srcval = mritmp;
      reshapefactor = 0; /* reset for output */
    }
    
    if(srcval->width != srcsurf->nvertices){
      fprintf(stderr,"ERROR: dimesion inconsitency in source data\n");
      fprintf(stderr,"       Number of surface vertices = %d\n",
        srcsurf->nvertices);
      fprintf(stderr,"       Number of value vertices = %d\n",srcval->width);
      exit(1);
    }

    if(srcval->nframes <= srcframe){
      printf("ERROR: desired frame (%d) exceeds number available (%d)\n",
       srcframe,srcval->nframes);
      exit(1);
    }
    
    for(vtx = 0; vtx < srcsurf->nvertices; vtx++)
      srcsurf->vertices[vtx].val = MRIFseq_vox(srcval,vtx,0,0,srcframe);
    MRIfree(&srcval); 
  }

  // ---------- Constrain with clabel ----------------------- //
  if(clabelfile){
    if(clabelinv){
      // Constrain to be OUTSIDE clabel by setting INSIDE values to 0
      for(n=0; n<clabel->n_points; n++){
	vtx = clabel->lv[n].vno;
	srcsurf->vertices[vtx].val = 0.0;
      }
    }
    else{ 
      // Constrain to be INSIDE clabel by OUTSIDE setting values to 0
      // Trickier -- use val2 to create a mask
      for(vtx = 0; vtx < srcsurf->nvertices; vtx++)
	srcsurf->vertices[vtx].val2 = 1;
      for(n=0; n<clabel->n_points; n++){
	vtx = clabel->lv[n].vno;
	srcsurf->vertices[vtx].val2 = 0;
      }
      for(vtx = 0; vtx < srcsurf->nvertices; vtx++){
	if(srcsurf->vertices[vtx].val2) 
	  srcsurf->vertices[vtx].val = 0.0;
      }
    }
  }

  /* Compute the overall max and min */
  for(vtx = 0; vtx < srcsurf->nvertices; vtx++){
    if(vtx == 0){
      overallmax = srcsurf->vertices[vtx].val;
      overallmaxvtx = vtx;
      overallmin = srcsurf->vertices[vtx].val;
      overallminvtx = vtx;
    }
    else if(overallmax < srcsurf->vertices[vtx].val){
      overallmax = srcsurf->vertices[vtx].val;
      overallmaxvtx = vtx;
    }
    else if(overallmin > srcsurf->vertices[vtx].val){
      overallmin = srcsurf->vertices[vtx].val;
      overallminvtx = vtx;
    }
  }
  printf("Done loading source values (nvtxs = %d)\n",srcsurf->nvertices);
  printf("overall max = %g at vertex %d\n",overallmax,overallmaxvtx);
  printf("overall min = %g at vertex %d\n",overallmin,overallminvtx);
  printf("surface nvertices %d\n",srcsurf->nvertices);

  /* Compute total cortex surface area */
  MRIScomputeMetricProperties(srcsurf) ;
  //printf("surface area %f\n",srcsurf->total_area);

  totarea = 0;
  for(vtx = 0; vtx < srcsurf->nvertices; vtx++)
    totarea += srcsurf->vertices[vtx].area;
  printf("surface area %f\n",totarea);

  //printf("Surface status %d\n",srcsurf->status);
  totarea = srcsurf->total_area;

  totarea = 0;
  for (vtx = 0 ; vtx < srcsurf->nfaces ; vtx++)
    totarea += srcsurf->faces[vtx].area;
  printf("surface area %f\n",totarea);


  /*---------------------------------------------------------*/
  /* This is where all the action is */
  printf("Searching for Clusters ...\n");
  scs = sclustMapSurfClusters(srcsurf,thmin,thmax,thsignid,
			      minarea,&NClusters,XFM);
  printf("Found %d clusters\n",NClusters);
  /*---------------------------------------------------------*/

  if(FixMNI){
    printf("INFO: fixing MNI talairach coordinates\n");
    for(n=0; n < NClusters; n++){
      FixMNITal(scs[n].xxfm, scs[n].yxfm, scs[n].zxfm, 
    &scs[n].xxfm, &scs[n].yxfm, &scs[n].zxfm); 
    }
  }

  if(csd != NULL){
    for(n=0; n < NClusters; n++){
      ClusterSize = scs[n].area;
      pval = CSDpvalClustSize(csd, ClusterSize, ciPct, &pvalLow, &pvalHi);
      scs[n].pval_clusterwise     = pval;
      scs[n].pval_clusterwise_low = pvalLow;
      scs[n].pval_clusterwise_hi  = pvalHi;
    }
  }

  if(debug){
    printf("-------------------------------------\n");
    DumpSurfClusterSum(stdout, scs, NClusters);
    printf("-------------------------------------\n");
  }

  /* ------ Print summary to a file ---------- */
  if(sumfile != NULL){
    fp = fopen(sumfile,"w");
    if(fp == NULL){
      printf("ERROR: could not open %s for writing\n",sumfile);
      exit(1);
    }
    fprintf(fp,"# Cluster Growing Summary (mri_surfcluster)\n");
    fprintf(fp,"# %s\n",vcid);
    fprintf(fp,"# %s\n",MRISurfSrcVersion());
    fprintf(fp,"# CreationTime %s\n",VERcurTimeStamp());
    fprintf(fp,"# cmdline %s\n",cmdline);
    fprintf(fp,"# cwd %s\n",cwd);
    fprintf(fp,"# sysname  %s\n",uts.sysname);
    fprintf(fp,"# hostname %s\n",uts.nodename);
    fprintf(fp,"# machine  %s\n",uts.machine);
    fprintf(fp,"# \n");
    fprintf(fp,"# Input      %s\n",srcid);  
    fprintf(fp,"# Frame Number      %d\n",srcframe);  
    fprintf(fp,"# srcsubj %s\n",srcsubjid);
    fprintf(fp,"# hemi %s\n",hemi);
    fprintf(fp,"# surface %s\n",srcsurfid);
    fprintf(fp,"# SUBJECTS_DIR %s\n",subjectsdir);

    fprintf(fp,"# Minimum Threshold %g\n",thmin);  
    if(thmax < 0) 
      fprintf(fp,"# Maximum Threshold infinity\n");
    else
      fprintf(fp,"# Maximum Threshold %g\n",thmax);
    fprintf(fp,"# Threshold Sign    %s\n",thsign);  

    fprintf(fp,"# Area Threshold    %g mm^2\n",minarea);  
    if(synthfunc != NULL)
      fprintf(fp,"# Synthesize        %s\n",synthfunc);  
    if(clabelfile){
      fprintf(fp,"# clabelfile %s\n",clabelfile);
      fprintf(fp,"# clabelinv  %d\n",clabelinv);
    }

    if(csd != NULL){
      fprintf(fp,"# CSD thresh  %lf\n",csd->thresh);      
      fprintf(fp,"# CSD nreps    %d\n",csd->nreps);      
      fprintf(fp,"# CSD simtype  %s\n",csd->simtype);      
      fprintf(fp,"# CSD contrast %s\n",csd->contrast);      
      fprintf(fp,"# CSD confint  %lf\n",ciPct);      
    }

    fprintf(fp,"# Overall max %g at vertex %d\n",overallmax,overallmaxvtx);
    fprintf(fp,"# Overall min %g at vertex %d\n",overallmin,overallminvtx);
    fprintf(fp,"# NClusters          %d\n",NClusters);  
    fprintf(fp,"# Total Cortical Surface Area %g (mm^2)\n",totarea);
    fprintf(fp,"# FixMNI = %d\n",FixMNI);  
    fprintf(fp,"# \n");  
    fprintf(fp,"# ClusterNo  Max   VtxMax   Size(mm^2)  TalX   TalY   TalZ ");
    if(csd != NULL)  fprintf(fp,"   CWP    CWPLow    CWPHi\n");
    else fprintf(fp,"\n");
    for(n=0; n < NClusters; n++){
      fprintf(fp,"%4d     %8.3f  %6d  %8.2f   %6.1f %6.1f %6.1f",
       n+1, scs[n].maxval, scs[n].vtxmaxval, scs[n].area,
       scs[n].xxfm, scs[n].yxfm, scs[n].zxfm);
    if(csd != NULL)  
      fprintf(fp,"  %7.5lf  %7.5lf  %7.5lf\n",
	      scs[n].pval_clusterwise,scs[n].pval_clusterwise_low,scs[n].pval_clusterwise_hi);
    else fprintf(fp,"\n");
    }
    fclose(fp);
  }

  if(ocpvalid != NULL){
    merged = MRIallocSequence(srcsurf->nvertices, 1, 1,MRI_FLOAT,4);
    // frame 0 will be filled in below with cluster-wise pval
    MRIcopyMRIS(merged, srcsurf, 1, "val"); // original data
    // frame 2 filled in below with thresholded
    MRIcopyMRIS(merged, srcsurf, 3, "undefval"); // cluster numbers
    // More below
  }

  /* --- Save the output as the thresholded input --- */
  sclustZeroSurfaceNonClusters(srcsurf);
  if(outid != NULL){
    printf("Saving thresholded output to  %s\n",outid);
    if(!strcmp(outfmt,"paint") || !strcmp(outfmt,"w"))
      MRISwriteValues(srcsurf,outid);
    else{
      mritmp = MRIcopyMRIS(NULL,srcsurf,0,"val");
      MRIwrite(mritmp,outid);
      MRIfree(&mritmp);
    }
  }
  if(ocpvalid != NULL){
    MRIcopyMRIS(merged, srcsurf, 2, "val"); // non-clusters removed
    // More below
  }

  /* --- Save the cluster number output --- */
  if(ocnid != NULL){
    printf("Saving cluster numbers to %s\n",ocnid);
    sclustSetSurfaceValToClusterNo(srcsurf);
    if(!strcmp(ocnfmt,"paint") || !strcmp(ocnfmt,"w"))
      MRISwriteValues(srcsurf,ocnid);
    else{
      mritmp = MRIcopyMRIS(NULL,srcsurf,0,"undefval");
      MRIwrite(mritmp,ocnid);
      MRIfree(&mritmp);
    }
  }

  /* --- Save the cluster pval --- */
  if(ocpvalid != NULL){
    sclustSetSurfaceValToCWP(srcsurf,scs);
    MRIcopyMRIS(merged, srcsurf, 0, "val"); // cluster-wise pval
    printf("Saving cluster pval %s\n",ocpvalid);
    MRIwrite(merged,ocpvalid);
  }

  /* -- Save output clusters as labels -- */
  if(outlabelbase != NULL){
    for(n=1; n <= NClusters; n++){
      sprintf(outlabelfile,"%s-%04d.label",outlabelbase,n);
      outlabel = LabelAlloc(scs[n-1].nmembers,srcsubjid,outlabelfile);
      outlabel->n_points =  scs[n-1].nmembers;
      nthlab = 0;
      for(vtx = 0; vtx < srcsurf->nvertices; vtx++){
	v = &srcsurf->vertices[vtx];
	if(v->undefval == n){
	  outlabel->lv[nthlab].x = v->x;
	  outlabel->lv[nthlab].y = v->y;
	  outlabel->lv[nthlab].z = v->z;
	  outlabel->lv[nthlab].vno = vtx;
	  outlabel->lv[nthlab].stat = v->val;
	  nthlab++;
	}
      }
      LabelWrite(outlabel,outlabelfile);
      LabelFree(&outlabel);
    } // End loop over clusters
  }

  return(0);
}
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--fixmni"))   FixMNI = 1;
    else if (!strcasecmp(option, "--nofixmni")) FixMNI = 0;
    else if (!strcasecmp(option, "--clabelinv")) clabelinv = 1;

    else if (!strcmp(option, "--hemi")){
      if(nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      if(!stringmatch(hemi,"lh") && !stringmatch(hemi,"rh")){
	printf("ERROR: hemi = %s, must be lh or rh\n",hemi);
	exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--src")){
      if(nargc < 1) argnerr(option,1);
      srcid = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	srcfmt = pargv[1]; nargsused ++;
      }
    }
    else if (!strcmp(option, "--srcsubj")){
      if(nargc < 1) argnerr(option,1);
      srcsubjid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcsurf")){
      if(nargc < 1) argnerr(option,1);
      srcsurfid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcframe")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&srcframe);
      nargsused = 1;
    }
    else if (!strcmp(option, "--xfm")){
      if(nargc < 1) argnerr(option,1);
      xfmfile = pargv[0]; nargsused = 1;
    }
    else if (!strcmp(option, "--csd")){
      if(nargc < 1) argnerr(option,1);
      csdfile = pargv[0]; 
      csd = CSDreadMerge(csdfile,csd);
      if(csd == NULL) exit(1);
      if(strcmp(csd->anattype,"surface")){
	printf("ERROR: csd must have anattype of surface\n");
	exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--thmin")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&thmin);
      if(thmin < 0) {
	printf("ERROR: thmin = %g, must be >= 0\n",thmin);
	printf("       use sign to set sign\n");
	exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--thmax")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&thmax);
      if(thmax < 0) {
	printf("ERROR: thmax = %g, must be >= 0\n",thmax);
	printf("       use sign to set sign\n");
	exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--thsign")){
      if(nargc < 1) argnerr(option,1);
      thsign = pargv[0];
      if(!stringmatch(thsign,"abs") && 
	 !stringmatch(thsign,"pos") && 
	 !stringmatch(thsign,"neg") ){
	printf("ERROR: thsign = %s, must be abs, pos, or neg\n",thsign);
	exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--minarea")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&minarea);
      nargsused = 1;
    }


    else if (!strcmp(option, "--mask")){
      if(nargc < 1) argnerr(option,1);
      maskid = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	maskfmt = pargv[1]; nargsused ++;
      }
    }
    else if (!strcmp(option, "--masksubj")){
      if(nargc < 1) argnerr(option,1);
      masksubjid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--masksign")){
      if(nargc < 1) argnerr(option,1);
      masksign = pargv[0];
      if(!stringmatch(masksign,"abs") && 
	 !stringmatch(masksign,"pos") && 
	 !stringmatch(masksign,"neg") ){
	printf("ERROR: masksign = %s, must be abs, pos, or neg\n",masksign);
	exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--maskthresh")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&maskthresh);
      if(maskthresh < 0) {
	printf("ERROR: maskthresh = %g, must be >= 0\n",maskthresh);
	printf("       use masksign to set sign\n");
	exit(1);
      }
      nargsused = 1;
    }
    
    else if (!strcmp(option, "--clabel")){
      if(nargc < 1) argnerr(option,1);
      clabelfile = pargv[0];
      nargsused = 1;
    }

    else if (!strcmp(option, "--label")){
      if(nargc < 1) argnerr(option,1);
      labelfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--nth")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nth);
      nargsused = 1;
    }
    else if (!strcmp(option, "--labelbase")){
      if(nargc < 1) argnerr(option,1);
      labelbase = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--labelsubj")){
      if(nargc < 1) argnerr(option,1);
      labelsubjid = pargv[0];
      nargsused = 1;
    }

    else if (!strcmp(option, "--omask")){
      if(nargc < 1) argnerr(option,1);
      omaskid = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	omaskfmt = pargv[1]; nargsused ++;
      }
    }
    else if (!strcmp(option, "--omasksubj")){
      if(nargc < 1) argnerr(option,1);
      omasksubjid = pargv[0];
      nargsused = 1;
    }

    else if (!strcmp(option, "--sum")){
      if(nargc < 1) argnerr(option,1);
      sumfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--o")){
      if(nargc < 1) argnerr(option,1);
      outid = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	outfmt = pargv[1]; nargsused ++;
      }
      else {
	outfmtid = mri_identify(outid);
	if(outfmtid != MRI_VOLUME_TYPE_UNKNOWN) 
	  outfmt = type_to_string(outfmtid);
      }
    }
    else if (!strcmp(option, "--ocn")){
      if(nargc < 1) argnerr(option,1);
      ocnid = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	ocnfmt = pargv[1]; nargsused ++;
      }
      else {
	ocnfmtid = mri_identify(ocnid);
	if(ocnfmtid != MRI_VOLUME_TYPE_UNKNOWN) 
	  ocnfmt = type_to_string(ocnfmtid);
      }
    }
    else if (!strcmp(option, "--ocp")){
      if(nargc < 1) argnerr(option,1);
      ocpvalid = pargv[0]; nargsused = 1;
      if(nth_is_arg(nargc, pargv, 1)){
	ocpvalfmt = pargv[1]; nargsused ++;
      }
    }
    else if (!strcmp(option, "--olab")){
      if(nargc < 1) argnerr(option,1);
      outlabelbase = pargv[0]; nargsused = 1;
    }

    else if (!strcmp(option, "--synth")){
      if(nargc < 1) argnerr(option,1);
      synthfunc = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--sd") || !strcmp(option, "--subjectsdir")){
      if(nargc < 1) argnerr(option,1);
      subjectsdir = pargv[0];
      nargsused = 1;
    }
    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
	fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --hemi hemi : lh or rh \n");
  printf("\n");
  printf("   --src      srcid <fmt> : source of surface values    \n");
  printf("   --srcsubj  subjid    : source surface subject (can be ico)\n");
  printf("   --srcsurf  surface   : get coorindates from surface (white)\n");
  printf("   --srcframe frameno   : 0-based frame number\n");
  printf("   --thmin    threshold : minimum intensity threshold\n");
  printf("   --thmax    threshold : maximum intensity threshold\n");
  printf("   --thsign   sign      : <abs>, pos, neg\n");
  printf("   --minarea  area      : area threshold for a cluster (mm^2)\n");
  printf("\n");
  printf("   --clabel labelfile : constrain to be within clabel\n");
  printf("   --clabelinv : constrain to be OUTSIDE clabel\n");
  //  printf("   --mask       maskid <fmt> \n");
  //  printf("   --maskthresh thresh \n");
  //  printf("   --masksign   sign : <abs>, pos, or neg \n");
  //  printf("   --masksubj   subjid : default that of src \n");
  //  printf("\n");
  //  printf("   --omask      omaskid <fmt> : o mask file id\n");
  //  printf("   --omasksubj  subjid  : saved in subject's space\n");
  //  printf("\n");
  //  printf("   --label     labelfile : output nth cluster as a label\n");
  //  printf("   --nth       nth       : nth for --label\n");
  //  printf("   --labelbase labelbase : output all labels\n");
  //  printf("   --labelsubj subjid    : default that of src \n");
  //  printf("\n");
  printf("   --sum sumfile     : text summary file\n");
  printf("   --o outid <fmt>   : input with non-clusters set to 0\n");
  printf("   --ocn ocnid <fmt>    : value is cluster number \n");
  printf("   --ocp ocpvalid <fmt> : value is cluster-wise pvalue (not ready yet)\n");
  printf("   --olab labelbase     : output clusters as labels \n");
  printf("\n");
  printf("   --xfm xfmfile     : talairach transform (def is talairach.xfm) \n");
  printf("   --<no>fixmni      : <do not> fix MNI talairach coordinates\n");
  printf("   --sd subjects_dir : (default is env SUBJECTS_DIR)\n");
  //  printf("   --synth synthfunc : uniform,loguniform,gaussian\n");
  printf("   --help : answers to ALL your questions\n");
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;

  printf("\n%s\n\n",vcid);

  printf(
"This program assigns each vertex on a cortical surface to a cluster \n"
"based on the distribution of intensity values in the source file. \n"
"A vertex must meet two criteria to be part of a cluster. First, \n"
"its intensity must fall within a certain range (the intensity threshold \n"
"criteria). Second, it must be part of a contiguous set of vertices \n"
"that meet the threshold criteria and whose surface area is greater  \n"
"than a minimum.  \n"
"\n"
"There are three types of output. (1) Summary file - simple text file\n"
"with a list of clusters found including maximum value, area, and \n"
"talairach coordinate. (2) Filtered - identical to the input except that\n"
"vertices not assigned to a cluster are set to 0. (3) Cluster number\n"
"surface - this has the same format as the input and filtered output\n"
"except that the value at each vertex is the cluster number assigned\n"
"to it. \n"
"\n"
"COMMAND-LINE ARGUMENTS\n"
"\n"
"--hemi hemi\n"
"\n"
"Specify the cortical hemisphere that the input represents. Valid values\n"
"are lh and rh.\n"
"\n"
"--src srcid <fmt>\n"
"\n"
"This is the input data to the clustering program. Fmt is the format \n"
"specification. Currently, only paint format is supported.\n"
"\n"
"--srcsurf surface\n"
"\n"
"This is the surface to use when computing the talairach coordinagtes.\n"
"Default is white.\n"
"\n"
"--srcframe frameno\n"
"\n"
"Zero-based frame number of the input file. Default is 0. For paint\n"
"format, zero is the only possible value.\n"
"\n"
"--thmin minthresh\n"
"\n"
"Minimum threshold in the intensity threshold criteria. See \n"
"SETTING THE CLUSTER INTENSITY THRESHOLD CRITERIA below.\n"
"\n"
"--thmax maxthresh\n"
"\n"
"Maximum threshold in the intensity threshold criteria. If negative,\n"
"then the maximum threshold is infinity. See SETTING THE CLUSTER \n"
"INTENSITY THRESHOLD CRITERIA below.\n"
"\n"
"--sign threshold sign\n"
"\n"
"This is used to control the sign of the threshold criteria. Legal\n"
"values are pos, neg, and abs. See SETTING THE CLUSTER INTENSITY \n"
"THRESHOLD CRITERIA below.\n"
"\n"
"--minarea area\n"
"\n"
"Minimum surface area (in mm^2) that a set of contiguous vertices\n"
"must achieve in order to be considered a cluster.\n"
"\n"
"--clabel labelfile\n"
"\n"
"Constrain cluster search to be inside or outside clabel. By default,\n"
"it will be constrained to be INSIDE. For OUTSIDE, specify --clabelinv.\n"
"clabel must be a surface-based label.\n"
"\n"
"--clabelinv \n"
"\n"
"Constrain cluster search to be OUTSIDE clabel. \n"
"\n"
"--sum summaryfile\n"
"\n"
"Text file in which to store the cluster summary. See SUMMARY FILE\n"
"OUTPUT below.\n"
"\n"
"--o outputid <fmt>\n"
"\n"
"File in which to store the surface values after setting the\n"
"non-cluster vertices to zero. Fmt is the format (currently, only\n"
"paint format is supported). Note: make sure to put a ./ in front\n"
"of outputid or else it will attempt to write the output to the\n"
"subject's surf directory.\n"
"\n"
"--ocn ocnid <fmt>\n"
"\n"
"File in which to store the cluster number of each vertex. This can be\n"
"useful for determining to which cluster a particular vertex\n"
"belongs. It can be viewed as with any other surface value file. Fmt is\n"
"the format (currently, only paint format is supported). Note: make\n"
"sure to put a ./ in front of outputid or else it will attempt to write\n"
"the output to the subject's surf directory.\n"
"\n"
"--olab outlabelbase\n"
"\n"
"Save output clusters as labels. There will be a label file for each\n"
"cluster. The name will be outlabelbase-XXXX.label, where XXXX is the\n"
"4-digit, zero-paded cluster number. The stat field in the label will\n"
"be the value of the input statistic at that vertex.\n"
"\n"
"--xfm xfmfile\n"
"\n"
"This is a transform file that is used to compute the Talairach \n"
"coordinates of a vertex for the summary file. The file must be\n"
"found in subjects_dir/subjectid/transforms. The default is\n"
"talairach.xfm which is based on the MNI atlas (see --fixmni).\n"
"\n"
"--fixmni --nofixmni\n"
"\n"
"Fix (or do not fix) MNI Talairach coordinates based on Matthew Brett's\n"
"transform. See http://www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html.\n"
"Default is to fix. The user may elect not to fix if the xfm file\n"
"is not talairach.xfm (see --xfm).\n"
"\n"
"--sd subjects_dir\n"
"\n"
"This allows the user to change the FreeSurfer subjects's directory\n"
"from the command-line. If unspecified, it will be set to the \n"
"environment variable SUBJECTS_DIR\n"
"\n"
"SETTING THE CLUSTER INTENSITY THRESHOLD CRITERIA\n"
"\n"
"Vertices on the surface are excluded or included, in part, based on\n"
"their intensity. For a vertex to be included, its value (ie,\n"
"intensity) must be within a certain range (namely, between the minimum\n"
"threshold, set with --thmin, and the maximum threshold, set with\n"
"--thmax). While thmin and thmax must be positive values, the sign of\n"
"the threshold can be set with --thsign. However, if thmax is negative,\n"
"then the maximum threshold will be set to infinity. For example, if\n"
"--thmin 2 --thmax 5 --thsign pos, then all vertices with values\n"
"between (positive) 2 and 5 will be candidates for clustering. However,\n"
"if --thsign abs is used instead, then all vertices between -2 and -5\n"
"as well as between +2 and +5 will be candidates.\n"
"\n"
"SUMMARY FILE OUTPUT\n"
"\n"
"The summary file (the argument of the --sum flag) will contain a \n"
"summary of the result of the clustering as well as a summary of the\n"
"conditions under which the clustering was performed. It will list\n"
"the clusters (1 to N) along with the maximum value found in the\n"
"cluster (Max), the vertex at which this maximum value was found\n"
"(VtxMax), the surface area of the cluster (Size), and the Talaiarach\n"
"coordinates of the maximum (based on talairach.xfm). A sample \n"
"summary file is shown below.\n"
"\n"
"Cluster Growing Summary (mri_surfcluster)\n"
"$Id: mri_surfcluster.c,v 1.19 2005/12/06 21:47:13 greve Exp $\n"
"Input :      minsig-0-lh.w\n"
"Frame Number:      0\n"
"Minimum Threshold: 5\n"
"Maximum Threshold: infinity\n"
"Threshold Sign:    pos\n"
"Area Threshold:    40 mm^2\n"
"NClusters          37\n"
"Total Cortical Surface Area 115576 (mm^2)\n"
"FixMNI = 1\n"
"\n"
"ClusterNo   Max  VtxMax  Size(mm^2)   TalX   TalY   TalZ\n"
"   1       44.6    6370    636.79    -33.6  -69.8   49.2\n"
"   2       40.3   48234    518.50     -2.9  -10.5   62.4\n"
"   3       39.5   54239    103.19    -51.1  -13.3   45.3\n"
"   4       39.5   55350     47.31    -50.1  -11.0   46.8\n"
"\n"
"BUGS\n"
"\n"
"Currently only supports paint (or w) format for input and output.\n"
"\n"
);



  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(csd != NULL){
    if(hemi == NULL) hemi = strcpyalloc(csd->hemi);
    else{
      if(strcmp(hemi,csd->hemi)){
	printf("ERROR: you have specified hemi=%s on cmdline, but\n",hemi);
	printf("CSD file was created with %s\n",csd->hemi);
	exit(1);
      }
    }
    if(srcsubjid == NULL) srcsubjid = strcpyalloc(csd->subject);
    else{
      if(strcmp(srcsubjid,csd->subject)){
	printf("ERROR: you have specified srcsubjid=%s on cmdline, but\n",srcsubjid);
	printf("CSD file was created with %s\n",csd->subject);
	exit(1);
      }
    }
    if(thmin < 0) thmin = csd->thresh;
    else{
      if(thmin != csd->thresh){
	printf("ERROR: you have specified thmin=%f on cmdline, but\n",thmin);
	printf("CSD file was created with %lf\n",csd->thresh);
	exit(1);
      }
    }
    if(minarea > 0){
      printf("ERROR: you cannot specify a minarea with --csd\n");
      exit(1);
    }
    minarea = 0;
  } // end csd != NULL


  if(hemi == NULL){
    printf("ERROR: hemi must be supplied\n");
    exit(1);
  }

  if(srcid == NULL){
    printf("ERROR: srcid must be supplied\n");
    exit(1);
  }
  if(srcsubjid == NULL){
    printf("ERROR: srcsubjid must be supplied\n");
    exit(1);
  }
  if(thmin < 0){
    printf("ERROR: thmin must be supplied\n");
    exit(1);
  }
  if(minarea < 0){
    printf("ERROR: minarea must be supplied\n");
    exit(1);
  }
  if(thsign == NULL) thsign = "abs";
  if(stringmatch(thsign,"pos")) thsignid = +1;
  if(stringmatch(thsign,"abs")) thsignid =  0;
  if(stringmatch(thsign,"neg")) thsignid = -1;


  if(srcframe < 0) srcframe = 0; 

  if(maskid != NULL){
    if(maskthresh < 0){
      printf("ERROR: must set mask thresh when specifying mask\n");
      exit(1);
    }
    if(masksubjid == NULL) masksubjid = srcsubjid;
    if(masksign == NULL) masksign = "abs";
  }

  if(labelsubjid == NULL) labelsubjid = srcsubjid;
  if(omasksubjid == NULL) omasksubjid = srcsubjid;


  /* check that the outputs can be written to */
  if(sumfile != NULL){
    if(! fio_DirIsWritable(sumfile,1)) {
      printf("ERROR: cannot write to %s\n",sumfile);
      exit(1);
    }
  }
  if(outid != NULL){
    if(! fio_DirIsWritable(outid,1)) {
      printf("ERROR: cannot write to %s\n",outid);
      exit(1);
    }
  }
  if(omaskid != NULL){
    if(! fio_DirIsWritable(omaskid,1)) {
      printf("ERROR: cannot write to %s\n",omaskid);
      exit(1);
    }
  }
  if(ocnid != NULL){
    if(! fio_DirIsWritable(ocnid,1)) {
      printf("ERROR: cannot write to %s\n",ocnid);
      exit(1);
    }
  }

  if(subjectsdir == NULL){
    subjectsdir = getenv("SUBJECTS_DIR");
    if(subjectsdir == NULL){
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  if(stringmatch(srcfmt,"paint") && srcframe != 0){
    printf("ERROR: for source format = paint, frame must be 0\n");
    exit(1);
  }

  if(clabelinv && clabelfile == NULL){
    printf("ERROR: must specify a clabel with --clabelinv\n");
    exit(1);
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"version %s\n",vcid);
  fprintf(fp,"hemi           = %s\n",hemi);
  fprintf(fp,"srcid          = %s %s\n",srcid,srcfmt);
  fprintf(fp,"srcsubjid      = %s\n",srcsubjid);
  fprintf(fp,"srcsurf        = %s\n",srcsurfid);
  fprintf(fp,"srcframe       = %d\n",srcframe);
  fprintf(fp,"thsign         = %s\n",thsign);
  fprintf(fp,"thmin          = %g\n",thmin);
  fprintf(fp,"thmax          = %g\n",thmax);
  fprintf(fp,"minarea        = %g\n",minarea);
  fprintf(fp,"xfmfile        = %s\n",xfmfile);
  if(maskid != NULL){
    fprintf(fp,"maskid         = %s %s\n",maskid, maskfmt);
    fprintf(fp,"masksubjid     = %s\n",masksubjid);
    fprintf(fp,"maskthresh     = %g\n",maskthresh);
    fprintf(fp,"masksign       = %s\n",masksign);
  }
  if(clabelfile != NULL){
    fprintf(fp,"clabelfile     = %s\n",clabelfile);
    fprintf(fp,"clabelinv      = %d\n",clabelinv);
  }
  if(labelfile != NULL)  fprintf(fp,"labelfile   = %s\n",labelfile);
  if(labelsubjid != NULL)fprintf(fp,"labelsubjid = %s\n",labelsubjid);
  if(labelbase != NULL)  fprintf(fp,"labelbase   = %s\n",labelbase);
  if(labelbase != NULL)  fprintf(fp,"nth         = %d\n",nth);

  if(outid != NULL)  fprintf(fp,"outid    = %s %s\n",outid, outfmt);

  if(ocnid != NULL)  fprintf(fp,"ocnid    = %s %s\n",ocnid, ocnfmt);

  if(sumfile != NULL) fprintf(fp,"sumfile  = %s\n",sumfile);
  if(omaskid != NULL){
    fprintf(fp,"omaskid    = %s %s\n",omaskid, omaskfmt);
    fprintf(fp,"omasksubj  = %s\n",omasksubjid);
  }
  fprintf(fp,"subjectsdir    = %s\n",subjectsdir);

  fprintf(fp,"FixMNI = %d\n",FixMNI);

  if(synthfunc != NULL) fprintf(fp,"synthfunc = %s\n",synthfunc);

  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int isflag(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth)
{
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if(nargc <= nth) return(0); 

  /* check whether the nth arg is a flag */
  if(isflag(argv[nth])) return(0);

  return(1);
}
/*------------------------------------------------------------*/
static int stringmatch(char *str1, char *str2)
{
  if(! strcmp(str1,str2)) return(1);
  return(0);
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if(n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}

