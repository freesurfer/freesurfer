/**
 * @brief Finds clusters on the surface.
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <sys/utsname.h>

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "matrix.h"
#include "mri.h"
#include "mri2.h"
#include "mri_identify.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "MRIio_old.h"
#include "fio.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "transform.h"
#include "version.h"
#include "annotation.h"
#include "fsenv.h"
#include "randomfields.h"


LABEL *MaskToSurfaceLabel(MRI *mask, double thresh, int sign);

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
static int  stringmatch(const char *str1, const char *str2);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

char *subjectdir = NULL;
char *hemi = NULL;

char *srcid = NULL;
const char *srcfmt  = "";
int   srcfmtid = MRI_VOLUME_TYPE_UNKNOWN;
const char *srcsurfid = "white";
char *srcsubjid = NULL;
int   srcframe = 0;
double thmin = -1, thmax = -1;
const char *thsign = NULL;
int   thsignid = 0;
float minarea = 0;
char *annotname = NULL;

char *maskid   = NULL;
char *maskfmt  = NULL;
int   maskfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *masksubjid  = NULL;
const char *masksignstr  = NULL;
int  masksign  = 0;
float maskthresh = 0.5; // assume binary

int  nth = -1;

char *omaskid = NULL;
const char *omaskfmt = "paint";
char *omasksubjid = NULL;

// Constraining label
char  *clabelfile=NULL;
char  *maskfile=NULL;
LABEL *clabel=NULL;
int   clabelinv = 0;
int   UseCortexLabel = 0;

char *outid = NULL;
const char *outfmt = "paint";
int   outfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *ocnid = NULL;
const char *ocnfmt = "paint";
int   ocnfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *ocpvalid = NULL;
const char *ocpvalfmt = "paint";
int   ocpvalfmtid = MRI_VOLUME_TYPE_UNKNOWN;
char *sumfile  = NULL;
char *pointsetfile  = NULL;

char *outlabelbase = NULL;
char outlabelfile[1000];
int  nthlab = 0;
LABEL *outlabel;
VERTEX *v;

char *outannot = NULL;

char *synthfunc  = NULL;
char *subjectsdir = NULL;
int debug = 0;

MRI *srcval, *mritmp,  *cnsurf;
MRI_SURFACE *srcsurf;
int  reshape = 1;
int  reshapefactor;

const char *xfmfile = "talairach.xfm";
char xfmpath[2000];
MATRIX *XFM;
int FixMNI = 1;

SURFCLUSTERSUM *scs,*scs2;

int overallmaxvtx,overallminvtx;
float overallmax,overallmin;

CSD *csd=NULL;
char *csdfile;
double pvalLow, pvalHi, ciPct=90, pval, ClusterSize;
char *csdpdffile = NULL;
int csdpdfonly = 0;
char *csdoutfile = NULL;

MRI *merged, *mask;

double thminadj, thmaxadj;
int AdjustThreshWhenOneTail=1;

char *voxwisesigfile=NULL;
MRI  *voxwisesig;
char *maxvoxwisesigfile=NULL;
char *maxcwpvalfile=NULL; // save p-value of largest cluster
int ReallyUseAverage7 = 0;
double fwhm = -1;
double fdr = -1;

double cwpvalthresh = -1; // pvalue, NOT log10(p)!
int Bonferroni = 0;
int BonferroniMax = 0;
int ReportCentroid = 0;
int sig2pmax = 0; // convert max value from -log10(p) to p

MRI *fwhmmap=NULL; // map of vertex-wise FWHM for non-stationary correction
char *maxareafile=NULL;

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  char fname[2000];
  int  n,nsearch,NClusters,NPrunedClusters,vtx, rt;
  FILE *fp;
  float totarea;
  int nargs,err;
  struct utsname uts;
  char *cmdline, cwd[2000];
  double cmaxsize,fwhmvtx;

  nargs = handleVersionOption(argc, argv, "mri_surfcluster");
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

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  dump_options(stdout);

  // Handle label: either read from mask or read from label file
  if (maskfile) {
    mask = MRIread(maskfile);
    if (mask == NULL) exit(1);
    clabel = MaskToSurfaceLabel(mask,maskthresh,masksign);
    if(clabel == NULL) {
      printf("No voxels found in mask, so exiting now\n");
      exit(0); // use exit 0 here for fspalm
    }
    printf("Found %d points in clabel.\n",clabel->n_points);
  }
  if (clabelfile) {
    printf("Loading clabel %s.\n",clabelfile);
    clabel = LabelRead(NULL, clabelfile);
    if (clabel == NULL) {
      fprintf(stderr,"ERROR reading %s\n",clabelfile);
      exit(1);
    }
    printf("Found %d points in clabel.\n",clabel->n_points);
  }

  sprintf(xfmpath,"%s/%s/mri/transforms/%s",subjectsdir,srcsubjid,xfmfile);
  //XFM = LoadxfmMatrix(xfmpath);
  XFM = DevolveXFM(srcsubjid, NULL, xfmfile);
  if (XFM == NULL) exit(1);

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

  if (annotname != NULL) {
    sprintf(fname,"%s/%s/label/%s.%s.annot",subjectsdir,srcsubjid,hemi,annotname);
    printf("Reading annotation %s\n",fname);
    err = MRISreadAnnotation(srcsurf, fname);
    if (err) {
      printf("ERROR: reading %s\n",fname);
      exit(1);
    }
    set_atable_from_ctable(srcsurf->ct);
  }

  printf("Computing metric properties\n");
  MRIScomputeMetricProperties(srcsurf);

  /* ------------------ load the source data ----------------------------*/
  printf("Loading source values\n");
  if (!strcmp(srcfmt,"curv")) { /* curvature file */
    sprintf(fname,"%s/%s/surf/%s.%s",subjectsdir,srcsubjid,hemi,srcid);
    printf("Reading curvature file %s\n",fname);
    rt = MRISreadCurvatureFile(srcsurf, fname);
    if (rt != NO_ERROR) exit(1);
    srcval = MRIallocSequence(srcsurf->nvertices, 1, 1,MRI_FLOAT,1);
    for (vtx = 0; vtx < srcsurf->nvertices; vtx++)
      srcsurf->vertices[vtx].val = srcsurf->vertices[vtx].curv;
  } else if (!strcmp(srcfmt,"paint") || !strcmp(srcfmt,"w")) {
    rt = MRISreadValues(srcsurf,srcid);
    if (rt != NO_ERROR) exit(1);
    srcval = MRIallocSequence(srcsurf->nvertices, 1, 1,MRI_FLOAT,1);
  } else {
    /* ----------- Use MRIreadType -------------------*/
    if (srcfmt == NULL) {
      srcval =  MRIread(srcid);
      if (srcval == NULL) {
        printf("ERROR: could not read %s\n",srcid);
        exit(1);
      }
    } else {
      srcfmtid = string_to_type(srcfmt);
      srcval =  MRIreadType(srcid,srcfmtid);
      if (srcval == NULL) {
        printf("ERROR: could not read %s as type %s\n",srcid,srcfmt);
        exit(1);
      }
    }
    if (srcval->height != 1 || srcval->depth != 1) {
      reshapefactor = srcval->height * srcval->depth;
      printf("Reshaping %d\n",reshapefactor);
      mritmp = mri_reshape(srcval, reshapefactor*srcval->width,
                           1, 1, srcval->nframes);
      MRIfree(&srcval);
      srcval = mritmp;
      reshapefactor = 0; /* reset for output */
    }

    if (srcval->width != srcsurf->nvertices) {
      fprintf(stderr,"ERROR: dimension inconsistency in source data\n");
      fprintf(stderr,"       Number of surface vertices = %d\n",
              srcsurf->nvertices);
      fprintf(stderr,"       Number of value vertices = %d\n",srcval->width);
      exit(1);
    }

    if (srcval->nframes <= srcframe) {
      printf("ERROR: desired frame (%d) exceeds number available (%d)\n",
             srcframe,srcval->nframes);
      exit(1);
    }

    for (vtx = 0; vtx < srcsurf->nvertices; vtx++)
      srcsurf->vertices[vtx].val = MRIgetVoxVal(srcval,vtx,0,0,srcframe);
    MRIfree(&srcval);
  }

  // In case of mask, Fill undefval = 1
  for (vtx = 0; vtx < srcsurf->nvertices; vtx++)
    srcsurf->vertices[vtx].undefval = 1;

  // ---------- Constrain with clabel ----------------------- //
  if (clabel) {
    if (clabelinv) {
      // Constrain to be OUTSIDE clabel by setting INSIDE values to 0
      for (n=0; n<clabel->n_points; n++) {
        vtx = clabel->lv[n].vno;
        srcsurf->vertices[vtx].val = 0.0;
        srcsurf->vertices[vtx].undefval = 0;
      }
    } else {
      // Constrain to be INSIDE clabel by OUTSIDE setting values to 0
      // Fill all undefvals with 0
      for (vtx = 0; vtx < srcsurf->nvertices; vtx++)
	srcsurf->vertices[vtx].undefval = 0;
      // Change undefvals in mask to 1
      for (n=0; n<clabel->n_points; n++) {
        vtx = clabel->lv[n].vno;
        srcsurf->vertices[vtx].undefval = 1;
      }
      // Now set the overlay values to 0
      for(vtx = 0; vtx < srcsurf->nvertices; vtx++) {
        if(! srcsurf->vertices[vtx].undefval)
          srcsurf->vertices[vtx].val = 0.0;
      }
    }
  }
  nsearch=0;
  for (vtx = 0; vtx < srcsurf->nvertices; vtx++)
    if(srcsurf->vertices[vtx].undefval) nsearch++;
  printf("number of voxels in search space = %d\n",nsearch);


  /* Compute the overall max and min */
  for (vtx = 0; vtx < srcsurf->nvertices; vtx++) {
    if (vtx == 0) {
      overallmax = srcsurf->vertices[vtx].val;
      overallmaxvtx = vtx;
      overallmin = srcsurf->vertices[vtx].val;
      overallminvtx = vtx;
    } else if (overallmax < srcsurf->vertices[vtx].val) {
      overallmax = srcsurf->vertices[vtx].val;
      overallmaxvtx = vtx;
    } else if (overallmin > srcsurf->vertices[vtx].val) {
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
  printf("metric props tot surface area %f\n",srcsurf->total_area);
  printf("group_avg_vtxarea_loaded %d\n",srcsurf->group_avg_vtxarea_loaded);

  // Need to constrain to mask and take into account average subject
  totarea = 0;
  for(vtx = 0 ; vtx < srcsurf->nvertices ; vtx++){
    if(srcsurf->vertices[vtx].undefval == 0) continue; // mask
    if(! srcsurf->group_avg_vtxarea_loaded) totarea += srcsurf->vertices[vtx].area;
    else                                    totarea += srcsurf->vertices[vtx].group_avg_area;
  }
  printf("masked surface area %f\n",totarea);

  if(voxwisesigfile) {
    double maxmaxsig;
    printf("Computing voxel-wise significance\n");
    srcval = MRIcopyMRIS(NULL,srcsurf,0,"val");
    voxwisesig = CSDpvalMaxSigMap(srcval, csd, NULL, NULL, &maxmaxsig, Bonferroni);
    MRIwrite(voxwisesig,voxwisesigfile);
    MRIfree(&srcval);
    MRIfree(&voxwisesig);
    if(maxvoxwisesigfile){
      fp = fopen(maxvoxwisesigfile,"w");
      fprintf(fp,"%10.5f\n",maxmaxsig);
      fclose(fp);
    }
  }

  if(fdr > 0){
    printf("Setting voxel-wise threshold with FDR = %lf\n",fdr);
    printf("Assuming input map is -log10(p)\n");
    err = MRISfdr2vwth(srcsurf, fdr, thsignid, 1, 1, &thmin);
    if(err){
      printf("ERROR: computing FDR threshold\n");
      exit(1);
    }
    printf("FDR Voxel-wise threshold is %g\n",thmin);
  }

  if(thsignid != 0 && AdjustThreshWhenOneTail) {
    // user has requested a tailed threshold, so adjust threshold
    // to account for a one-tailed test. This requires that the
    // input be -log10(p), where p is computed from a two-tailed
    // test. One could recompute the p-values in the volume to
    // convert to a one-tailed test, but easier to change the threshold
    printf("Adjusting threshold for 1-tailed test.\n");
    printf("If the input is not a -log10(p) volume, re-run with --no-adjust.\n");
    thminadj = thmin - log10(2.0);
    if (thmax > 0) thmaxadj = thmax - log10(2.0);
    else           thmaxadj = thmax;
  } else {
    printf("NOT Adjusting threshold for 1-tailed test\n");
    thminadj = thmin;
    thmaxadj = thmax;
  }
  printf("thminadj = %g\n",thminadj);

  /*---------------------------------------------------------*/
  /* This is where all the action is */
  printf("Searching for Clusters ...\n");
  printf("thmin=%f (%f), thmax=%f (%g), thsignid=%d, minarea=%lf\n",
         thmin,thminadj,thmax,thmaxadj,thsignid,minarea);
  scs = sclustMapSurfClusters(srcsurf,thminadj,thmaxadj,thsignid,
                              minarea,&NClusters,XFM, fwhmmap);
  printf("Found %d clusters\n",NClusters);
  cmaxsize = sclustMaxClusterArea(scs, NClusters);
  printf("Max cluster size %lf\n",cmaxsize);

  /*---------------------------------------------------------*/

  if (FixMNI) {
    printf("INFO: fixing MNI talairach coordinates\n");
    for (n=0; n < NClusters; n++) {
      FixMNITal(scs[n].xxfm, scs[n].yxfm, scs[n].zxfm,
                &scs[n].xxfm, &scs[n].yxfm, &scs[n].zxfm);
    }
  }

  if (csd != NULL) {
    for (n=0; n < NClusters; n++) {
      ClusterSize = scs[n].area;
      pval = CSDpvalClustSize(csd, ClusterSize, ciPct, &pvalLow, &pvalHi);
      scs[n].pval_clusterwise     = pval;
      scs[n].pval_clusterwise_low = pvalLow;
      scs[n].pval_clusterwise_hi  = pvalHi;
    }
  }
  if(fwhm > 0) { // use RFT
    //fwhmvtx = fwhm/sqrt(totarea/nsearch);
    double grfnsearch;
    grfnsearch = nsearch;
    if(thsignid == 0) grfnsearch = grfnsearch/2.0;  // hack for abs
    fwhmvtx = fwhm/srcsurf->avg_vertex_dist; 
    printf("GRF: fwhm = %g, dist = %g, fwhmvtx = %g, thmin = %g, grfnsearch = %g\n",
	   fwhm,srcsurf->avg_vertex_dist,fwhmvtx,thmin, grfnsearch);
    for (n=0; n < NClusters; n++) {
      ClusterSize = scs[n].nmembers;
      pval = RFprobZClusterSigThresh(ClusterSize, thmin, fwhmvtx, grfnsearch, 2);
      if(thsignid == 0) pval = 2*pval; // hack for abs
      //printf("  %2d %g %g\n",n,ClusterSize, pval);
      scs[n].pval_clusterwise     = pval;
      scs[n].pval_clusterwise_low = 0;
      scs[n].pval_clusterwise_hi  = 0;
    }
  }

  if(Bonferroni > 0){
    // Bonferroni correction -- generally for across spaces
    for (n=0; n < NClusters; n++) {
      pval = scs[n].pval_clusterwise;
      pval = 1 - pow((1-pval),Bonferroni);
      scs[n].pval_clusterwise = pval;

      pval = scs[n].pval_clusterwise_low;
      pval = 1 - pow((1-pval),Bonferroni);
      scs[n].pval_clusterwise_low = pval;

      pval = scs[n].pval_clusterwise_hi;
      pval = 1 - pow((1-pval),Bonferroni);
      scs[n].pval_clusterwise_hi = pval;
    }
  }

  /* Sort the clusters by p-value */
  scs2 = SortSurfClusterSum(scs, NClusters);
  free(scs);
  scs = scs2;

  if(maxcwpvalfile != NULL){
    // save sig of largest cluster
    fp = fopen(maxcwpvalfile,"w");
    if(NClusters > 1) fprintf(fp,"%10.5lf\n",scs->pval_clusterwise);
    else              fprintf(fp,"1.0\n");
    fclose(fp);
  }

  /* Remove clusters that do not meet the minimum clusterwise pvalue */
  if(cwpvalthresh > 0 && (fwhm >0 || csd != NULL) ){
    printf("Pruning by CW P-Value %g\n",cwpvalthresh);
    scs2 = sclustPruneByCWPval(scs, NClusters, cwpvalthresh, &NPrunedClusters, srcsurf);
    NClusters = NPrunedClusters;
    scs = scs2;
  }


  if (debug) {
    printf("-------------------------------------\n");
    DumpSurfClusterSum(stdout, scs, NClusters);
    printf("-------------------------------------\n");
  }

  if(pointsetfile)
    sclustSaveAsPointSet(pointsetfile, scs, NClusters, srcsurf);

  /* ------ Print summary to a file ---------- */
  if (sumfile != NULL) {
    fp = fopen(sumfile,"w");
    if (fp == NULL) {
      printf("ERROR: could not open %s for writing\n",sumfile);
      exit(1);
    }
    fprintf(fp,"# Cluster Growing Summary (mri_surfcluster)\n");
    fprintf(fp,"# %s\n",getVersion().c_str());
    fprintf(fp,"# %s\n",getVersion().c_str());
    fprintf(fp,"# CreationTime %s\n",VERcurTimeStamp());
    fprintf(fp,"# cmdline %s\n",cmdline);
    fprintf(fp,"# cwd %s\n",cwd);
    fprintf(fp,"# sysname  %s\n",uts.sysname);
    fprintf(fp,"# hostname %s\n",uts.nodename);
    fprintf(fp,"# machine  %s\n",uts.machine);
    fprintf(fp,"# FixVertexAreaFlag %d\n",MRISgetFixVertexAreaValue());
    fprintf(fp,"# FixSurfClusterArea %d\n",FixSurfClusterArea);
    fprintf(fp,"# \n");
    fprintf(fp,"# Input      %s\n",srcid);
    fprintf(fp,"# Frame Number      %d\n",srcframe);
    fprintf(fp,"# srcsubj %s\n",srcsubjid);
    fprintf(fp,"# hemi %s\n",hemi);
    fprintf(fp,"# surface %s\n",srcsurfid);
    fprintf(fp,"# group_avg_surface_area %g\n",srcsurf->group_avg_surface_area);
    fprintf(fp,"# group_avg_vtxarea_loaded %d\n",srcsurf->group_avg_vtxarea_loaded);
    if (annotname) fprintf(fp,"# annot %s\n",annotname);
    fprintf(fp,"# SUBJECTS_DIR %s\n",subjectsdir);

    if(fdr>0) fprintf(fp,"# FDR %lf\n",fdr);
    fprintf(fp,"# SearchSpace_mm2 %g\n",totarea);
    fprintf(fp,"# SearchSpace_vtx %d\n",nsearch);
    fprintf(fp,"# Bonferroni %d\n",Bonferroni);
    fprintf(fp,"# Minimum Threshold %g\n",thmin);
    if (thmax < 0)
      fprintf(fp,"# Maximum Threshold infinity\n");
    else
      fprintf(fp,"# Maximum Threshold %g\n",thmax);
    fprintf(fp,"# Threshold Sign    %s\n",thsign);
    fprintf(fp,"# AdjustThreshWhenOneTail %d\n",AdjustThreshWhenOneTail);
    if(cwpvalthresh > 0)
      fprintf(fp,"# CW PValue Threshold: %g \n",cwpvalthresh);

    fprintf(fp,"# Area Threshold    %g mm^2\n",minarea);
    if (synthfunc != NULL)
      fprintf(fp,"# Synthesize        %s\n",synthfunc);
    if (clabelfile) {
      fprintf(fp,"# clabelfile %s\n",clabelfile);
      fprintf(fp,"# clabelinv  %d\n",clabelinv);
    }

    if (csd != NULL) {
      fprintf(fp,"# CSD thresh  %lf\n",csd->thresh);
      fprintf(fp,"# CSD nreps    %d\n",csd->nreps);
      fprintf(fp,"# CSD simtype  %s\n",csd->simtype);
      fprintf(fp,"# CSD contrast %s\n",csd->contrast);
      fprintf(fp,"# CSD confint  %lf\n",ciPct);
    }
    if(fwhm > 0) fprintf(fp,"# fwhm %lf\n",fwhm);

    fprintf(fp,"# Overall max %g at vertex %d\n",overallmax,overallmaxvtx);
    fprintf(fp,"# Overall min %g at vertex %d\n",overallmin,overallminvtx);
    fprintf(fp,"# NClusters          %d\n",NClusters);
    //fprintf(fp,"# Total Cortical Surface Area %g (mm^2)\n",totarea); // why needed?
    fprintf(fp,"# FixMNI = %d\n",FixMNI);
    fprintf(fp,"# \n");
    if(FixMNI) fprintf(fp,"# ClusterNo  Max   VtxMax   Size(mm^2)  TalX   TalY   TalZ ");
    else       fprintf(fp,"# ClusterNo  Max   VtxMax   Size(mm^2)  MNIX   MNIY   MNIZ ");

    if(csd != NULL)  fprintf(fp,"   CWP    CWPLow    CWPHi");
    else if(fwhm > 0) fprintf(fp,"  GRFCWP");
    fprintf(fp,"   NVtxs    WghtVtx");
    if (annotname != NULL)  fprintf(fp,"   Annot");
    fprintf(fp,"\n");

    for (n=0; n < NClusters; n++) {
      double maxval = scs[n].maxval;
      if(sig2pmax) {
	maxval = pow(10.0,-fabs(scs[n].maxval));
	if(BonferroniMax > 1) maxval *= BonferroniMax;
      }
      if(ReportCentroid){ // Report Centriod XYZ
	fprintf(fp,"%4d     %9.4f  %6d  %8.2f   %6.1f %6.1f %6.1f",
		n+1, maxval, scs[n].vtxmaxval, scs[n].area,
		scs[n].cxxfm, scs[n].cyxfm, scs[n].czxfm);
      } else { // Report Max XYZ
	fprintf(fp,"%4d     %9.4f  %6d  %8.2f   %6.1f %6.1f %6.1f",
		n+1, maxval, scs[n].vtxmaxval, scs[n].area,
		scs[n].xxfm, scs[n].yxfm, scs[n].zxfm);
      }
      if (csd != NULL)
        fprintf(fp,"  %7.5lf  %7.5lf  %7.5lf", scs[n].pval_clusterwise,
                scs[n].pval_clusterwise_low,scs[n].pval_clusterwise_hi);
      else if (fwhm > 0)
	fprintf(fp,"  %7.5lf",scs[n].pval_clusterwise);
      fprintf(fp,"  %5d  %10.2f",scs[n].nmembers,scs[n].weightvtx);
      if (annotname != NULL)
        fprintf(fp,"  %s",
                annotation_to_name(srcsurf->vertices[scs[n].vtxmaxval].annotation,
                                   NULL));
      fprintf(fp,"\n");

    }
    fclose(fp);
  }

  if(maxareafile){
    fp = fopen(maxareafile,"w");
    fprintf(fp,"%10.4f\n",scs[0].area); // assuming 0 is max
    fclose(fp);
  }

  if (ocpvalid != NULL) {
    merged = MRIallocSequence(srcsurf->nvertices, 1, 1,MRI_FLOAT,4);
    // frame 0 will be filled in below with cluster-wise pval
    MRIcopyMRIS(merged, srcsurf, 1, "val"); // original data
    // frame 2 filled in below with thresholded
    MRIcopyMRIS(merged, srcsurf, 3, "undefval"); // cluster numbers
    // More below
  }

  /* --- Save the output as the thresholded input --- */
  sclustZeroSurfaceNonClusters(srcsurf);
  if (outid != NULL) {
    printf("Saving thresholded output to  %s\n",outid);
    if (!strcmp(outfmt,"paint") || !strcmp(outfmt,"w"))
      MRISwriteValues(srcsurf,outid);
    else {
      mritmp = MRIcopyMRIS(NULL,srcsurf,0,"val");
      MRIwrite(mritmp,outid);
      MRIfree(&mritmp);
    }
  }
  if (ocpvalid != NULL) {
    MRIcopyMRIS(merged, srcsurf, 2, "val"); // non-clusters removed
    // More below
  }

  /* --- Save the cluster number output --- */
  if (ocnid != NULL) {
    printf("Saving cluster numbers to %s\n",ocnid);
    sclustSetSurfaceValToClusterNo(srcsurf);
    if (!strcmp(ocnfmt,"paint") || !strcmp(ocnfmt,"w"))
      MRISwriteValues(srcsurf,ocnid);
    else {
      mritmp = MRIcopyMRIS(NULL,srcsurf,0,"undefval");
      MRIwrite(mritmp,ocnid);
      MRIfree(&mritmp);
    }
  }

  /* --- Save the cluster pval --- */
  if (ocpvalid != NULL) {
    sclustSetSurfaceValToCWP(srcsurf,scs);
    MRIcopyMRIS(merged, srcsurf, 0, "val"); // cluster-wise -log10(pval)
    printf("Saving cluster pval %s\n",ocpvalid);
    MRIwrite(merged,ocpvalid);
  }

  /* -- Save output clusters as labels -- */
  if(outlabelbase != NULL) {
    for (n=1; n <= NClusters; n++) {
      sprintf(outlabelfile,"%s-%04d.label",outlabelbase,n);
      outlabel = LabelAlloc(scs[n-1].nmembers,srcsubjid,outlabelfile);
      outlabel->n_points =  scs[n-1].nmembers;
      nthlab = 0;
      for (vtx = 0; vtx < srcsurf->nvertices; vtx++) {
        v = &srcsurf->vertices[vtx];
        if (v->undefval == n) {
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

  if(outannot != NULL){
    printf("Constructing output annotation\n");
    sclustAnnot(srcsurf, NClusters);
    printf("Writing annotation %s\n",outannot);
    MRISwriteAnnotation(srcsurf, outannot);
  }


  return(0);
}
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  FILE *fp;

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
    else if (!strcasecmp(option, "--fixmni"))   FixMNI = 1;
    else if (!strcasecmp(option, "--nofixmni")) FixMNI = 0;
    else if (!strcasecmp(option, "--no-fixmni")) FixMNI = 0;
    else if (!strcasecmp(option, "--clabelinv")) clabelinv = 1;
    else if (!strcasecmp(option, "--mask-inv")) clabelinv = 1;
    else if (!strcasecmp(option, "--no-adjust")) AdjustThreshWhenOneTail=0;
    else if (!strcmp(option, "--csdpdf-only")) csdpdfonly = 1;
    else if (!strcmp(option, "--cortex")) UseCortexLabel = 1;
    else if (!strcmp(option, "--centroid")) ReportCentroid = 1;
    else if (!strcmp(option, "--sig2p-max")) sig2pmax = 1;
    else if (!strcmp(option, "--no-fix-vertex-area")) {
      printf("Turning off fixing of vertex area\n");
      MRISsetFixVertexAreaValue(0);
    } 
    else if (!strcasecmp(option, "--diag")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      printf("Gdiag_no = %d\n",Gdiag_no);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--bonferroni")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&Bonferroni);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--bonferroni-max")) {
      if (nargc < 1) argnerr(option,1);
      // only applies when --sig2pmax is used
      sscanf(pargv[0],"%d",&BonferroniMax);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--hemi")) {
      if (nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      if (!stringmatch(hemi,"lh") && !stringmatch(hemi,"rh")) {
        printf("ERROR: hemi = %s, must be lh or rh\n",hemi);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--src") || !strcmp(option, "--in")) {
      if (nargc < 1) argnerr(option,1);
      srcid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        srcfmt = pargv[1];
        nargsused ++;
      }
    } else if (!strcmp(option, "--really-use-average7")) ReallyUseAverage7 = 1;
    else if (!strcmp(option, "--srcsubj") || !strcmp(option, "--subject")) {
      if (nargc < 1) argnerr(option,1);
      srcsubjid = pargv[0];
      if (!strcmp(srcsubjid,"average7")) {
        if (!ReallyUseAverage7) {
          printf("\n");
          printf("ERROR: you have selected subject average7. It is recommended that\n");
          printf("you use the fsaverage subject in $FREESURFER_HOME/subjects.\n");
          printf("If you really want to use average7, re-run this program with\n");
          printf("--really-use-average7 as the first argument.\n");
          printf("\n");
          exit(1);
        } else {
          printf("\n");
          printf("INFO: you have selected subject average7 (and REALLY want to use it)\n");
          printf("instead of fsaverage. So I'm going to turn off fixing of vertex area\n");
          printf("to maintain compatibility with the pre-stable3 release.\n");
          printf("\n");
          MRISsetFixVertexAreaValue(0);
        }
      }
      nargsused = 1;
    } else if (!strcmp(option, "--srcsurf") || !strcmp(option, "--surf")) {
      if (nargc < 1) argnerr(option,1);
      srcsurfid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srcframe") || !strcmp(option, "--frame")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&srcframe);
      nargsused = 1;
    } else if (!strcmp(option, "--xfm")) {
      if (nargc < 1) argnerr(option,1);
      xfmfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--csd")) {
      if (nargc < 1) argnerr(option,1);
      csdfile = pargv[0];
      csd = CSDreadMerge(csdfile,csd);
      if (csd == NULL) exit(1);
      if (strcmp(csd->anattype,"surface")) {
        printf("ERROR: csd must have anattype of surface\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--csd-out")) {
      if(nargc < 1) argnerr(option,1);
      csdoutfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--csdpdf")) {
      if (nargc < 1) argnerr(option,1);
      csdpdffile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--fdr")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&fdr);
      nargsused = 1;
    } else if (!strcmp(option, "--thmin")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&thmin);
      if (thmin < 0) {
        printf("ERROR: thmin = %g, must be >= 0\n",thmin);
        printf("       use sign to set sign\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--thmax")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&thmax);
      if (thmax < 0) {
        printf("ERROR: thmax = %g, must be >= 0\n",thmax);
        printf("       use sign to set sign\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--thsign") || !strcmp(option, "--sign")) {
      if (nargc < 1) argnerr(option,1);
      thsign = pargv[0];
      if (!stringmatch(thsign,"abs") &&
          !stringmatch(thsign,"pos") &&
          !stringmatch(thsign,"neg") ) {
        printf("ERROR: thsign = %s, must be abs, pos, or neg\n",thsign);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--minarea")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&minarea);
      nargsused = 1;
    } else if (!strcmp(option, "--mask-sign")) {
      if (nargc < 1) argnerr(option,1);
      masksignstr = pargv[0];
      if (!stringmatch(masksignstr,"abs") &&
          !stringmatch(masksignstr,"pos") &&
          !stringmatch(masksignstr,"neg") ) {
        printf("ERROR: masksign = %s, must be abs, pos, or neg\n",masksignstr);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--mask-thresh")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&maskthresh);
      if (maskthresh < 0) {
        printf("ERROR: maskthresh = %g, must be >= 0\n",maskthresh);
        printf("       use masksign to set sign\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--clabel")) {
      if (nargc < 1) argnerr(option,1);
      clabelfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--mask")) {
      if (nargc < 1) argnerr(option,1);
      maskfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--annot")) {
      if (nargc < 1) argnerr(option,1);
      annotname = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--nth")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nth);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--sum")) {
      if (nargc < 1) argnerr(option,1);
      sumfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--maxareafile")) {
      if (nargc < 1) argnerr(option,1);
      maxareafile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--pointset")) {
      if (nargc < 1) argnerr(option,1);
      pointsetfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--o")) {
      if (nargc < 1) argnerr(option,1);
      outid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        outfmt = pargv[1];
        nargsused ++;
      } else {
        outfmtid = mri_identify(outid);
        if (outfmtid != MRI_VOLUME_TYPE_UNKNOWN)
          outfmt = type_to_string(outfmtid);
      }
    } else if (!strcmp(option, "--ocn")) {
      if (nargc < 1) argnerr(option,1);
      ocnid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        ocnfmt = pargv[1];
        nargsused ++;
      } else {
        ocnfmtid = mri_identify(ocnid);
        if (ocnfmtid != MRI_VOLUME_TYPE_UNKNOWN)
          ocnfmt = type_to_string(ocnfmtid);
      }
    } else if (!strcmp(option, "--cwsig") || !strcmp(option, "--ocp") ) {
      if (nargc < 1) argnerr(option,1);
      ocpvalid = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        ocpvalfmt = pargv[1];
        nargsused ++;
      }
    } else if (!strcmp(option, "--cwpvalthresh") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&cwpvalthresh);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--vwsig")) {
      if(nargc < 1) argnerr(option,1);
      voxwisesigfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--vwsigmax")) {
      if(nargc < 1) argnerr(option,1);
      maxvoxwisesigfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--maxcwpval")) {
      if(nargc < 1) argnerr(option,1);
      maxcwpvalfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--olab")) {
      if (nargc < 1) argnerr(option,1);
      outlabelbase = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--oannot")) {
      if (nargc < 1) argnerr(option,1);
      outannot= pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--synth")) {
      if (nargc < 1) argnerr(option,1);
      synthfunc = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--sd") || !strcmp(option, "--subjectsdir")) {
      if (nargc < 1) argnerr(option,1);
      subjectsdir = pargv[0];
      FSENVsetSUBJECTS_DIR(subjectsdir);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--fwhm")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&fwhm);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--fwhm-map")) {
      if (nargc < 1) argnerr(option,1);
      fwhmmap = MRIread(pargv[0]);
      if(fwhmmap==NULL) exit(1);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--use-avg-vertex-area")) {
      setenv("FS_CLUSTER_USE_AVG_VERTEX_AREA","1",1);
    } 
    else if (!strcmp(option, "--no-use-avg-vertex-area")) {
      setenv("FS_CLUSTER_USE_AVG_VERTEX_AREA","0",1);
    } 
    else if (!strcmp(option, "--fwhmdat")) {
      if(nargc < 1) argnerr(option,1);
      fp = fopen(pargv[0],"r");
      if(fp == NULL){
	printf("ERROR: opening %s\n",pargv[0]);
	exit(1);
      }
      fscanf(fp,"%lf",&fwhm);
      fclose(fp);
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
  printf("\n");
  printf("   --in infile : source of surface values \n");
  printf("\n");
  printf("   --thmin  threshold : minimum intensity threshold\n");
  printf("   --sign   sign      : <abs> or pos/neg for one-sided tests\n");
  printf("   --no-adjust  : do not adjust thresh for one-tailed tests\n");
  printf("   --fdr FDR  : set thmin with False Discovery Rate\n");
  printf("\n");
  printf("   --subject  subjid    : source surface subject (can be ico)\n");
  printf("   --hemi hemi          : lh or rh \n");
  printf("   --surf     surface   : get coorindates from surface (white)\n");
  printf("   --annot    annotname : report annotation for max vertex (eg, aparc)\n");
  printf("   --frame frameno      : 0-based source frame number\n");
  printf("\n");
  printf("   --csd csdfile <--csd csdfile ...>\n");
  printf("   --vwsig vwsig : map of corrected voxel-wise significances\n");
  printf("   --cwsig cwsig : map of cluster-wise significances\n");
  printf("   --maxcwpval maxcwpvalfile.txt : save p-value of the largest (max) cluster (1.0 if no cluster)\n");
  printf("   --bonferroni N : addition correction across N (eg, spaces)\n");
  printf("   --sig2p-max : convert max from sig to p\n");
  printf("   --bonferroni-max N : apply bonf cor to maximum (only applies with --sig2p-max)\n");
  printf("   --csdpdf csdpdffile\n");
  printf("   --csdpdf-only : write csd pdf file and exit.\n");
  printf("   --csd-out out.csd : write out merged csd files as one.\n");
  printf("   --cwpvalthresh cwpvalthresh : clusterwise threshold\n");
  printf("\n");
  printf("   --fwhm    fwhm     :  fwhm in mm2 for GRF\n");
  printf("   --fwhmdat fwhm.dat :  text file with fwhm in mm2 for GRF\n");
  printf("\n");
  printf("   --clabel labelfile : constrain to be within clabel\n");
  printf("   --cortex : set clabel to be subject/label/hemi.cortex.label\n");
  printf("   --mask maskfile : constrain to be within mask\n");
  printf("   --mask-inv : constrain to be OUTSIDE mask or clabel\n");
  printf("   --centroid : report centroid instead of location of maximum stat\n");
  printf("   --sum sumfile     : text summary file\n");
  printf("   --pointset pointsetfile : file that can be read into freeview with -c\n");
  printf("   --maxareafile file : write area of largest cluster to file\n");
  printf("   --o outid        : input with non-clusters set to 0\n");
  printf("   --ocn ocnid      : value is cluster number \n");
  printf("   --olab labelbase : output clusters as labels \n");
  printf("   --oannot annotname : output clusters as an annotation \n");
  printf("\n");
  printf("   --minarea  area      : area threshold for a cluster (mm^2)\n");
  printf("   --xfm xfmfile     : talairach transform (def is talairach.xfm) \n");
  printf("   --<no>fixmni      : <do not> fix MNI talairach coordinates\n");
  printf("   --sd subjects_dir : (default is env SUBJECTS_DIR)\n");
  printf("   --thmax  threshold : maximum intensity threshold (only use if you know what you are doing)\n");
  //  printf("   --synth synthfunc : uniform,loguniform,gaussian\n");
  printf("   --no-fix-vertex-area : turn off fixing of vertex area (for back comapt only)\n");
  printf("   --help : answers to ALL your questions\n");
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf("\n%s\n\n",getVersion().c_str());

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
    "--in inputfile \n"
    "\n"
    "This is the input data to the clustering program.\n"
    "Reads all formats supported by mri_convert.\n"
    "\n"
    "--subject subjectid\n"
    "\n"
    "Surface values are defined on this subject.\n"
    "\n"
    "--hemi hemi\n"
    "\n"
    "Specify the cortical hemisphere that the input represents. Valid values\n"
    "are lh and rh.\n"
    "\n"
    "--surf surface\n"
    "\n"
    "This is the surface to use when computing the talairach coordinagtes.\n"
    "Default is white.\n"
    "\n"
    "--annot annotationname\n"
    "\n"
    "Report the cortical annotation that the maximum in each label falls into. \n"
    "The annotation used will be SUBJECTS_DIR/subject/label/hemi.annotationname.annot.\n"
    "Eg, --annot aparc\n"
    "\n"
    "--frame frameno\n"
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
    "--no-adjust  \n"
    "\n"
    "Do not adjust thresh for one-tailed tests. By default, the threshold\n"
    "will be adjusted when the --sign is pos or neg by subtracting log10(2.0).\n"
    "This assumes several things: (1) the source is a -log10(p) map, and (2)\n"
    "the pvalue was computed using a two-sided t-test. Under these conditions,\n"
    "subtracting log10(2.0) is like dividing the p-value by 2 (ie, changing it\n"
    "from a two-tailed test to a one-tailed test). If the input map does not\n"
    "meet these criteria, then run with --no-adjust.\n"
    "\n"
    "--fdr FDR\n"
    "\n"
    "Set thmin with False Discovery Rate. 0 < FDR < 1. Usually something \n"
    "like .01 or .05. Only when input is -log10(p).\n"
    "\n"
    "--csd csdfile <--csd csdfile>\n"
    "\n"
    "Load one or more CSD files. CSD stands for 'Cluster Simulation Data'. This\n"
    "file is produced by running mri_glmfit with --sim. The the threshold and hemi\n"
    "info are obtained from the CSD file and cannot be specified on the command-\n"
    "line. If more than one CSD is specified, they are merged into one CSD internally.\n"
    "When a CSD file is specified, three more columns are added to the summary table:\n"
    "  1. CWP - cluster-wise pvalue. The pvalue of the cluster corrected for \n"
    "     multiple comparisons\n"
    "  2. CWPLow - lower 90%% confidence limit of CWP based on binomial distribution\n"
    "  3. CWPHi  - upper 90%% confidence limit of CWP based on binomial distribution\n"
    "In addition, the user can specify --cwsig, which saves the sigificance map of \n"
    "the clusters in which the value of each vertex is the -log10(pvalue) of cluster\n"
    "to which the vertex belongs (the cluster-wise significance). The user can also\n"
    "specify that the vertex-wise significance be computed and saved  with --vwsig.\n"
    "The significance is based on the distribution of the maximum significances \n"
    "found during the CSD simulation.\n"
    "\n"
    "--csdpdf csdpdfile\n"
    "\n"
    "Compute PDF/CDF of CSD data and save in csdpdffile. This is mostly good for debugging.\n"
    "\n"
    "--csdpdf-only\n"
    "\n"
    "Only write the csd pdf file.\n"
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
    "--mask-inv \n"
    "\n"
    "Constrain cluster search to be OUTSIDE mask or clabel. \n"
    "\n"
    "--sum summaryfile\n"
    "\n"
    "Text file in which to store the cluster summary. See SUMMARY FILE\n"
    "OUTPUT below.\n"
    "\n"
    "--o outputid \n"
    "\n"
    "File in which to store the surface values after setting the\n"
    "non-cluster vertices to zero. Fmt is the format (currently, only\n"
    "paint format is supported). Note: make sure to put a ./ in front\n"
    "of outputid or else it will attempt to write the output to the\n"
    "subject's surf directory.\n"
    "\n"
    "--ocn ocnid \n"
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
static void check_options(void) {
  char tmpstr[2000];
  int err;

  if(csdoutfile){
    if(csd == NULL){
      printf("ERROR: need --csd with --csd-out \n");
      exit(1);
    }
    printf("Merging CSD files\n");
    err = CSDwrite(csdoutfile, csd);
    if(err) exit(1);
  }

  if(csdpdffile) {
    if (csd == NULL) {
      printf("ERROR: need --csd with --csdpdf");
      exit(1);
    }
    printf("Creating CDFs from CSD files\n");
    CSDpdf(csd,-1);
    CSDwritePDF(csdpdffile,csd);
    if(csdpdfonly) exit(0);
  }

  if (csd != NULL) {
    if (hemi == NULL) hemi = strcpyalloc(csd->hemi);
    else {
      if (strcmp(hemi,csd->hemi)) {
        printf("ERROR: you have specified hemi=%s on cmdline, but\n",hemi);
        printf("CSD file was created with %s\n",csd->hemi);
        exit(1);
      }
    }
    if (srcsubjid == NULL) srcsubjid = strcpyalloc(csd->subject);
    else {
      if (strcmp(srcsubjid,csd->subject)) {
        printf("ERROR: you have specified srcsubjid=%s on cmdline, but\n",srcsubjid);
        printf("CSD file was created with %s\n",csd->subject);
        exit(1);
      }
    }
    if (thmin < 0) thmin = csd->thresh;
    else {
      if (thmin != csd->thresh) {
        printf("ERROR: you have specified thmin=%f on cmdline, but\n",thmin);
        printf("CSD file was created with %lf\n",csd->thresh);
        exit(1);
      }
    }
    if (minarea > 0) {
      printf("ERROR: you cannot specify a minarea with --csd\n");
      exit(1);
    }
    printf("csd->threshsign = %g\n",csd->threshsign);
    if(csd->threshsign >  0.5)      thsign = "pos";
    else if(csd->threshsign < -0.5) thsign = "neg";
    else                            thsign = "abs";
    minarea = 0;
  } // end csd != NULL
  if (voxwisesigfile != NULL && csd == NULL) {
    printf("ERROR: need csd with --vwsig\n");
    exit(1);
  }

  if(thsign == NULL) thsign = "abs";
  if(stringmatch(thsign,"pos")) thsignid = +1;
  if(stringmatch(thsign,"abs")) thsignid =  0;
  if(stringmatch(thsign,"neg")) thsignid = -1;
  printf("thsign = %s, id = %d\n",thsign,thsignid);

  if (hemi == NULL) {
    printf("ERROR: hemi must be supplied\n");
    exit(1);
  }

  if (srcid == NULL) {
    printf("ERROR: srcid must be supplied\n");
    exit(1);
  }
  if (srcsubjid == NULL) {
    printf("ERROR: srcsubjid must be supplied\n");
    exit(1);
  }
  if(thmin < 0 && fdr < 0) {
    printf("ERROR: thmin or fdr must be supplied\n");
    exit(1);
  }
  if (srcframe < 0) srcframe = 0;

  if (maskid != NULL) {
    if (maskthresh < 0) {
      printf("ERROR: must set mask thresh when specifying mask\n");
      exit(1);
    }
    if (masksubjid == NULL) masksubjid = srcsubjid;
    if (masksignstr == NULL) masksignstr = "abs";
    if (stringmatch(masksignstr,"abs")) masksign = 0;
    if (stringmatch(masksignstr,"pos")) masksign = 1;
    if (stringmatch(masksignstr,"neg")) masksign = -1;
  }

  /* check that the outputs can be written to */
  if (sumfile != NULL) {
    if (! fio_DirIsWritable(sumfile,1)) {
      printf("ERROR: cannot write to %s\n",sumfile);
      exit(1);
    }
  }
  if (outid != NULL) {
    if (! fio_DirIsWritable(outid,1)) {
      printf("ERROR: cannot write to %s\n",outid);
      exit(1);
    }
  }
  if (omaskid != NULL) {
    if (! fio_DirIsWritable(omaskid,1)) {
      printf("ERROR: cannot write to %s\n",omaskid);
      exit(1);
    }
  }
  if (ocnid != NULL) {
    if (! fio_DirIsWritable(ocnid,1)) {
      printf("ERROR: cannot write to %s\n",ocnid);
      exit(1);
    }
  }

  if (subjectsdir == NULL) {
    subjectsdir = getenv("SUBJECTS_DIR");
    if (subjectsdir == NULL) {
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  if (stringmatch(srcfmt,"paint") && srcframe != 0) {
    printf("ERROR: for source format = paint, frame must be 0\n");
    exit(1);
  }

  if(UseCortexLabel && clabelfile){
    printf("ERROR: cannot specify both --clabel and --cortex\n");
    exit(1);
  }
  if(UseCortexLabel){
    sprintf(tmpstr,"%s/%s/label/%s.cortex.label",subjectsdir,srcsubjid,hemi);
    clabelfile = strcpyalloc(tmpstr);
  }

  if (clabelfile != NULL && maskfile != NULL) {
    printf("ERROR: cannot specify both --clabel and --mask\n");
    exit(1);
  }
  if (clabelinv && (clabelfile == NULL && maskfile == NULL)) {
    printf("ERROR: must specify a --clabel or --mask with --mask-inv\n");
    exit(1);
  }

  //if(fwhm > 0 && !strcmp(thsign,"abs")){
  //printf("ERROR: you must specify a pos or neg sign with --fwhm\n");
  //exit(1);
  //}

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"version %s\n",getVersion().c_str());
  fprintf(fp,"hemi           = %s\n",hemi);
  fprintf(fp,"srcid          = %s %s\n",srcid,srcfmt);
  fprintf(fp,"srcsubjid      = %s\n",srcsubjid);
  fprintf(fp,"srcsurf        = %s\n",srcsurfid);
  fprintf(fp,"srcframe       = %d\n",srcframe);
  fprintf(fp,"thsign         = %s\n",thsign);
  fprintf(fp,"thmin          = %g\n",thmin);
  fprintf(fp,"thmax          = %g\n",thmax);
  fprintf(fp,"fdr            = %g\n",fdr);
  fprintf(fp,"minarea        = %g\n",minarea);
  if(Bonferroni) fprintf(fp,"Bonferroni      = %d\n",Bonferroni);
  fprintf(fp,"xfmfile        = %s\n",xfmfile);
  if (maskid != NULL) {
    fprintf(fp,"maskid         = %s %s\n",maskid, maskfmt);
    fprintf(fp,"masksubjid     = %s\n",masksubjid);
    fprintf(fp,"maskthresh     = %g\n",maskthresh);
    fprintf(fp,"masksign       = %s\n",masksignstr);
  }
  if (clabelfile != NULL) {
    fprintf(fp,"clabelfile     = %s\n",clabelfile);
    fprintf(fp,"clabelinv      = %d\n",clabelinv);
  }
  fprintf(fp,"nth         = %d\n",nth);

  if (outid != NULL)  fprintf(fp,"outid    = %s %s\n",outid, outfmt);

  if (ocnid != NULL)  fprintf(fp,"ocnid    = %s %s\n",ocnid, ocnfmt);

  if (sumfile != NULL) fprintf(fp,"sumfile  = %s\n",sumfile);
  if (omaskid != NULL) {
    fprintf(fp,"omaskid    = %s %s\n",omaskid, omaskfmt);
    fprintf(fp,"omasksubj  = %s\n",omasksubjid);
  }
  fprintf(fp,"subjectsdir    = %s\n",subjectsdir);

  fprintf(fp,"FixMNI = %d\n",FixMNI);

  if (synthfunc != NULL) fprintf(fp,"synthfunc = %s\n",synthfunc);

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
static int stringmatch(const char *str1, const char *str2) {
  if (! strcmp(str1,str2)) return(1);
  return(0);
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
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

/*--------------------------------------------------------*/
LABEL *MaskToSurfaceLabel(MRI *mask, double thresh, int sign) {
  int v,c,r,s,nhits,ishit;
  LABEL *label;
  double val;

  nhits = 0;
  for (s=0; s < mask->depth; s++) {
    for (r=0; r < mask->height; r++) {
      for (c=0; c < mask->width; c++) { // cols must be fasest
        ishit = 0;
        val = MRIgetVoxVal(mask,c,r,s,0);
        if (sign == 0 && fabs(val) > thresh) ishit = 1;
        if (sign >  0 && val > +thresh)      ishit = 1;
        if (sign <  0 && val < -thresh)      ishit = 1;
        if (ishit) nhits ++;
      }
    }
  }
  printf("Found %d vertices in mask\n",nhits);
  if (nhits == 0) return(NULL);

  label = LabelAlloc(nhits, "unknown", "mask");
  label->n_points = nhits;

  nhits = 0;
  v = 0;
  for (s=0; s < mask->depth; s++) {
    for (r=0; r < mask->height; r++) {
      for (c=0; c < mask->width; c++) { // cols must be fasest
        ishit = 0;
        val = MRIgetVoxVal(mask,c,r,s,0);
        if (sign == 0 && fabs(val) > thresh) ishit = 1;
        if (sign >  0 && val > +thresh)      ishit = 1;
        if (sign <  0 && val < -thresh)      ishit = 1;
        if (ishit) {
          // only assign vertex no because it's a surface label
          label->lv[nhits].vno = v;
          nhits++;
        }
        v++;
      }
    }
  }

  printf("Found %d vertices in mask\n",nhits);
  return(label);
}


