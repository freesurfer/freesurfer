/*
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

  0. Test - never enough

  6. Free function seems to be broken

  ------- Done --------
  1. Prune by distance (in subject's space ?)
  2. MNI2Tal
  3. Output Thresholded volume
  4. Mask volume (absolute)
  5. Label output (native)
  7. Size threshold in mm^3
  8. Synthesize

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri_identify.h"
#include "matrix.h"
#include "registerio.h"
#include "resample.h"
#include "cmdargs.h"
#include "fio.h"
#include "utils.h"
#include "volcluster.h"
#include "version.h"
#include "randomfields.h"
#include "pdf.h"

static MATRIX *LoadMNITransform(char *regfile, int ncols, int nrows,
                                int nslices, MATRIX **ppCRS2FSA,
                                MATRIX **ppFSA2Func,
                                float *colres, float *rowres,
                                float *sliceres);


static MRI *MRIsynthUniform(int ncols, int nrows, int nslices,
                            int nframes, MRI *tvol);
static MRI *MRIsynthLogUniform(int ncols, int nrows, int nslices,
                               int nframes, MRI *tvol);
static double Gaussian01PDF(void);
static MRI *MRIsynthGaussian(int ncols, int nrows, int nslices,
                             int nframes, MRI *tvol);
static MRI *MRIbinarize01(MRI *vol, float thmin, float thmax,
                          const char *thsign, int invert,
                          int lowval, int highval, int *nhits, MRI *binvol);


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static int  singledash(char *flag);
static void dump_options(FILE *fp);
double round(double); // why is this never defined?!?

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

static char tmpstr[2000];

int debug = 0;
int verbose = 0;

char *volid   = NULL;
char *regfile = NULL;
int   frame   = 0;
int   intype = MRI_VOLUME_TYPE_UNKNOWN;
char  *intypestring;

char *maskid   = NULL;
int   masktype = MRI_VOLUME_TYPE_UNKNOWN;
char  *masktypestring;
float maskthresh = 0.5;
const char  *masksignstring = "abs";
int   maskinvert = 0;
int   maskframe  = 0;

char *outmaskid   = NULL;
int   outmasktype = MRI_VOLUME_TYPE_UNKNOWN;
char  *outmasktypestring;

char *outcnid   = NULL;
int   outcntype = MRI_VOLUME_TYPE_UNKNOWN;
char  *outcntypestring;

char *outid   = NULL;
char *synthfunction;
int   outtype = MRI_VOLUME_TYPE_UNKNOWN;
char  *outtypestring;

char *sumfile;

int nlabelcluster = -1;
char *labelfile;
char *labelbase;

float threshmin  = -1.0;
float threshmax  = -1.0;
const char  *signstring = "abs";
int   threshsign =    0;
float sizethresh    = 0.0;
int   sizethreshvox = 0;
float distthresh =   0.0;
int   allowdiag  = 0;
int sig2pmax = 0; // convert max value from -log10(p) to p

MRI *vol, *HitMap, *outvol, *maskvol, *binmask;
VOLCLUSTER **ClusterList, **ClusterList2;
MATRIX *CRS2MNI, *CRS2FSA, *FSA2Func;
LABEL *label;

float colres, rowres, sliceres, voxsize;

FILE *fpsum;

int fixtkreg = 1;
double threshminadj, threshmaxadj;
int AdjustThreshWhenOneTail=1;

double cwpvalthresh = -1; // pvalue, NOT log10(p)!

CSD *csd=NULL;
char *csdfile;
double pvalLow, pvalHi, ciPct=90, pval, ClusterSize;
char *csdpdffile = NULL;
int csdpdfonly = 0;

char *voxwisesigfile=NULL;
char *maxvoxwisesigfile=NULL;
MRI  *voxwisesig, *clustwisesig;
char *clustwisesigfile=NULL;

char *SUBJECTS_DIR=NULL;
char *subject=NULL;
char *segvolfile=NULL;
char *segvolpath=NULL;
MRI *segvol0=NULL;
MRI *segvol=NULL;
char *ctabfile;
COLOR_TABLE *ctab = NULL;
int ctabindex;

double fwhm = -1;
int nmask;
double searchspace;
int FixMNI = 1;
MATRIX *vox2vox;
int UseFSAverage = 0;
struct utsname uts;
char *cmdline, cwd[2000];
char *segctabfile = NULL;
COLOR_TABLE *segctab = NULL;
int Bonferroni = 0;
int BonferroniMax = 0;
int regheader = 0;
char *pointset = NULL;

/*--------------------------------------------------------------*/
/*--------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------*/
int main(int argc, char **argv) {
  int nhits, *hitcol, *hitrow, *hitslc,nargs;
  int col, row, slc;
  int nthhit, n, m, nclusters, nprunedclusters;
  float x,y,z,val,pval;
  char *stem;
  COLOR_TABLE *ct;
  FILE *fp;

  setRandomSeed(53);

  nargs = handleVersionOption(argc, argv, "mri_volcluster");
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

  if (FixMNI) VolClustFixMNI = 1;

  if (debug) dump_options(stdout);

  /* Load the input volume */
  vol = MRIread(volid);
  if (vol == NULL) {
    fprintf(stderr,"ERROR: reading %s\n",volid);
    exit(1);
  }

  if(UseFSAverage){
    if(vol->xsize != 1 && vol->xsize != 2){
      printf("ERROR: voxel size is %g, must be 1 or 2 with --fsaverage\n",vol->xsize);
      exit(1);
    }
    sprintf(tmpstr,"%s/average/mni305.cor.subfov%d.reg",
	    getenv("FREESURFER_HOME"),(int)round(vol->xsize));
    regfile = strcpyalloc(tmpstr);
    printf("Using %s as regfile\n",regfile);
  }

  /* Check that the specified frame is within range */
  if (frame >= vol->nframes) {
    fprintf(stderr,"ERROR: specified frame = %d, >= nframes\n",frame);
    exit(1);
  }

  /* Load the mask volume */
  if (maskid != NULL) {
    printf("INFO: loading mask volume: %s\n",maskid);
    maskvol = MRIread(maskid);
    if (maskvol == NULL) {
      fprintf(stderr,"ERROR: reading %s\n",maskid);
      exit(1);
    }
    if (vol->width != maskvol->width &&
        vol->depth != maskvol->depth &&
        vol->height != maskvol->height) {
      printf("ERROR: input volume and mask volume dimensions are "
             "inconsistent\n");
      exit(1);
    }
    if (maskframe >= maskvol->nframes) {
      printf("ERROR: mask frame = %d >= nframes in mask (%d)\n",
             maskframe,maskvol->nframes);
      exit(1);

    }

    binmask = MRIbinarize01(maskvol, maskthresh, -1, masksignstring,
                            maskinvert, 0, 1, &nmask, NULL);
    if (binmask == NULL) exit(1);

    //if (outmaskid != NULL) MRIwriteType(binmask,outmaskid,outmasktype);
    if(outmaskid != NULL) MRIwrite(binmask,outmaskid);
    MRIfree(&maskvol);
    printf("Found %d voxels in mask\n",nmask);
  } else {
    binmask = NULL;
    nmask = vol->width  * vol->height * vol->depth;
  }

  /* Load the resolution and geometry information from the register.dat */
  if (regfile != NULL)
    CRS2MNI =
      LoadMNITransform(regfile, vol->width,vol->height,vol->depth,
                       &CRS2FSA, &FSA2Func, &colres, &rowres, &sliceres);
  else  {
    CRS2MNI = MatrixIdentity(4,NULL);
    CRS2FSA = MatrixIdentity(4,NULL);
    FSA2Func = MatrixIdentity(4,NULL);
  }

  colres   = vol->xsize;
  rowres   = vol->ysize;
  sliceres = vol->zsize;
  voxsize = colres * rowres * sliceres;
  if (debug) {
    printf("VolumeRes: %g %g %g (%g)\n",colres,rowres,sliceres,voxsize);
    printf("Registration: ---------------\n");
    MatrixPrint(stdout,FSA2Func);
    printf("CRS2FSA: ---------------\n");
    MatrixPrint(stdout,CRS2FSA);
    printf("CRS2MNI: ---------------\n");
    MatrixPrint(stdout,CRS2MNI);
  }
  searchspace = nmask * voxsize;
  printf("Search Space = %g mm3\n",searchspace);

  if (segvolpath != NULL) {
    segvol0 = MRIread(segvolpath);
    if(segvol0 == NULL) exit(1);
    segvol  = MRIcloneBySpace(vol,-1,1);
    vox2vox = MRIvoxToVoxFromTkRegMtx(vol, segvol0, FSA2Func);
    vox2vox = MatrixInverse(vox2vox,NULL);
    MRIvol2Vol(segvol0,segvol,vox2vox,SAMPLE_NEAREST,0);
    MatrixFree(&vox2vox);
    MRIfree(&segvol0);
    if(segctabfile == NULL){
      segctabfile = (char *) calloc(sizeof(char),1000);
      sprintf(segctabfile,"%s/FreeSurferColorLUT.txt",getenv("FREESURFER_HOME"));
    }
    printf("Using ctab %s\n",segctabfile);
    segctab = CTABreadASCII(segctabfile);
    if (segctab == NULL) {
      printf("ERROR: reading %s\n",segctabfile);
      exit(1);
    }
  }

  if (sizethresh == 0) sizethresh = sizethreshvox*voxsize;

  /* Replace data with synthetic if desired */
  if (synthfunction != NULL) {
    printf("INFO: synthsizing with %s\n",synthfunction);
    if (!strcmp(synthfunction,"uniform"))
      MRIsynthUniform(0,0,0,0,vol);
    if (!strcmp(synthfunction,"loguniform"))
      MRIsynthLogUniform(0,0,0,0,vol);
    if (!strcmp(synthfunction,"gaussian"))
      MRIsynthGaussian(0,0,0,0,vol);
  }

  if(voxwisesigfile) {
    double maxmaxsig;
    printf("Computing voxel-wise significance\n");
    voxwisesig = CSDpvalMaxSigMap(vol, csd, binmask, NULL, &maxmaxsig, Bonferroni);
    MRIwrite(voxwisesig,voxwisesigfile);
    if(maxvoxwisesigfile){
      fp = fopen(maxvoxwisesigfile,"w");
      fprintf(fp,"%10.5f\n",maxmaxsig);
      fclose(fp);
    }
  }


  /* Initialize the hit map - this is a map of voxels that have been
     accounted for as either outside of a the threshold range or
     belonging to a cluster. The initialization is for thresholding */
  HitMap = clustInitHitMap(vol, frame, threshminadj, threshmaxadj, threshsign,
                           &nhits, &hitcol, &hitrow, &hitslc,
                           binmask, maskframe);
  if (HitMap == NULL) {
    printf("ERROR: initializing hit map\n");
    if(nhits == 0){
      printf("  No voxels were found that met the threshold criteria");      
      if(binmask) printf(" within the mask");      
      printf(".\n");
    }
    exit(1);
  }
  //MRIwriteType(HitMap,"hitmap",BSHORT_FILE);

  printf("INFO: Found %d voxels in threhold range\n",nhits);

  /* Allocate an array of clusters equal to the number of hits -- this
     is the maximum number of clusters possible */
  ClusterList = clustAllocClusterList(nhits);
  if (ClusterList == NULL) {
    fprintf(stderr,"ERROR: could not alloc %d clusters\n",nhits);
    exit(1);
  }

  nclusters = 0;
  for (nthhit = 0; nthhit < nhits; nthhit ++) {

    /* Determine whether this hit is still valid. It may
       not be if it was assigned to the cluster of a
       previous hit */
    col = hitcol[nthhit];
    row = hitrow[nthhit];
    slc = hitslc[nthhit];
    if (MRIgetVoxVal(HitMap,col,row,slc,0)) continue;

    /* Grow cluster using this hit as a seed */
    ClusterList[nclusters] = clustGrow(col,row,slc,HitMap,allowdiag);

    /* Determine the member with the maximum value */
    clustMaxMember(ClusterList[nclusters], vol, frame, threshsign);

    //clustComputeXYZ(ClusterList[nclusters],CRS2FSA); /* for FSA coords */
    clustComputeTal(ClusterList[nclusters],CRS2MNI); /*"true" Tal coords */

    /* increment the number of clusters */
    nclusters ++;
  }

  printf("INFO: Found %d clusters that meet threshold criteria\n",
         nclusters);

  if (debug)
    clustDumpClusterList(stdout,ClusterList, nclusters, vol, frame);

  /* Remove clusters that do not meet the minimum size requirement */
  ClusterList2 = clustPruneBySize(ClusterList,nclusters,
                                  voxsize,sizethresh,
                                  &nprunedclusters);
  //clustFreeClusterList(&ClusterList,nclusters);/* Free - does not work */
  nclusters = nprunedclusters;
  ClusterList = ClusterList2;

  printf("INFO: Found %d clusters that meet size criteria\n",
         nclusters);

  if (debug)
    clustDumpClusterList(stdout,ClusterList, nclusters, vol, frame);

  /* Remove clusters that do not meet the minimum distance requirement */
  if (distthresh > 0.0) {
    printf("INFO: pruning by distance %g\n",distthresh);
    ClusterList2 = clustPruneByDistance(ClusterList,nclusters,
                                        distthresh, &nprunedclusters);
    //clustFreeClusterList(&ClusterList,nclusters);/* Free - does not work */
    nclusters = nprunedclusters;
    ClusterList = ClusterList2;
  }

  /* Sort Clusters */
  ClusterList2 = clustSortClusterList(ClusterList2,nclusters,NULL);
  //clustFreeClusterList(&ClusterList,nclusters);/* Free - does not work */
  ClusterList = ClusterList2;

  printf("INFO: Found %d final clusters\n",nclusters);
  if (debug)
    clustDumpClusterList(stdout,ClusterList, nclusters, vol, frame);

  if(csd != NULL) {
    for (n=0; n < nclusters; n++) {
      ClusterSize = ClusterList[n]->nmembers * voxsize;
      pval = CSDpvalClustSize(csd, ClusterSize, ciPct, &pvalLow, &pvalHi);
      ClusterList[n]->pval_clusterwise     = pval;
      ClusterList[n]->pval_clusterwise_low = pvalLow;
      ClusterList[n]->pval_clusterwise_hi  = pvalHi;
    }
  }
  if(fwhm > 0) {
    double grfsearchspace;
    int D=0;
    if(vol->depth == 1) {
      D = 2;
      grfsearchspace = nmask * colres * rowres;
    }
    else{
      D = 3;
      grfsearchspace = nmask * colres * rowres * sliceres;
    }
    for (n=0; n < nclusters; n++) {
      if(D==2) ClusterSize = ClusterList[n]->nmembers * colres * rowres;
      if(D==3) ClusterSize = ClusterList[n]->nmembers * colres * rowres * sliceres;
      if(threshsign == 0) grfsearchspace = grfsearchspace/2.0; // This is a hack for abs
      if(AdjustThreshWhenOneTail) {
	pval = RFprobZClusterSigThresh(ClusterSize, threshmin, fwhm, grfsearchspace, D);
	if(threshsign == 0) pval = 2*pval;
      }
      else{
	// When no adjustment is made and GRF is used, assumes that the 
	// input volume is a z-volume. This is a bit messy.
	pval = RFprobZCluster(ClusterSize, threshmin, fwhm, grfsearchspace, D);
      }
      if(threshsign == 0) pval = 2*pval; // This is a hack for abs
      ClusterList[n]->pval_clusterwise = pval;
    }
  }
  if(Bonferroni > 0){
    // Bonferroni correction -- generally for across spaces
    for (n=0; n < nclusters; n++) {
      pval = ClusterList[n]->pval_clusterwise;
      pval = 1 - pow((1-pval),Bonferroni);
      ClusterList[n]->pval_clusterwise = pval;

      pval = ClusterList[n]->pval_clusterwise_low;
      pval = 1 - pow((1-pval),Bonferroni);
      ClusterList[n]->pval_clusterwise_low = pval;

      pval = ClusterList[n]->pval_clusterwise_hi;
      pval = 1 - pow((1-pval),Bonferroni);
      ClusterList[n]->pval_clusterwise_hi = pval;
    }
  }
  /* Sort Clusters */
  ClusterList2 = clustSortClusterList(ClusterList2,nclusters,NULL);
  //clustFreeClusterList(&ClusterList,nclusters);/* Free - does not work */
  ClusterList = ClusterList2;

  /* Remove clusters that do not meet the minimum clusterwise pvalue */
  if(cwpvalthresh > 0 && (fwhm >0 || csd != NULL) ){
    printf("Pruning by CW P-Value %g\n",cwpvalthresh);
    ClusterList2 = clustPruneByCWPval(ClusterList,nclusters,
                                        cwpvalthresh, &nprunedclusters);
    nclusters = nprunedclusters;
    ClusterList = ClusterList2;
  }

  /* Open the Summary File (or set its pointer to stdout) */
  if (sumfile != NULL) {
    fpsum = fopen(sumfile,"w");
    if (fpsum == NULL) {
      printf("ERROR: could not open %s for writing\n",sumfile);
      exit(1);
    }
    printf("INFO: writing summary to %s\n",sumfile);
  } else fpsum = stdout;

  /* Dump summary to file or stdout */
  fprintf(fpsum,"# Cluster Growing Summary (mri_volcluster)\n");
  fprintf(fpsum,"# %s\n",getVersion().c_str());
  fprintf(fpsum,"# cwd %s\n",cwd);
  fprintf(fpsum,"# cmdline %s\n",cmdline);
  if(SUBJECTS_DIR) fprintf(fpsum,"# SUBJECTS_DIR  %s\n",SUBJECTS_DIR);
  fprintf(fpsum,"# sysname  %s\n",uts.sysname);
  fprintf(fpsum,"# hostname %s\n",uts.nodename);
  fprintf(fpsum,"# machine  %s\n",uts.machine);
  fprintf(fpsum,"# user     %s\n",VERuser());

  fprintf(fpsum,"# Input Volume:      %s\n",volid);
  fprintf(fpsum,"# Frame Number:      %d\n",frame);
  fprintf(fpsum,"# VoxSize_mm3 %g\n",voxsize);
  fprintf(fpsum,"# SearchSpace_mm3 %g\n",searchspace);
  fprintf(fpsum,"# SearchSpace_vox %d\n",nmask);
  fprintf(fpsum,"# Minimum Threshold: %g\n",threshmin);
  fprintf(fpsum,"# Bonferroni %d\n",Bonferroni);
  if (threshmax < 0)
    fprintf(fpsum,"# Maximum Threshold: inifinity\n");
  else
    fprintf(fpsum,"# Maximum Threshold: %g\n",threshmax);
  fprintf(fpsum,"# Threshold Sign:    %s\n",signstring);
  fprintf(fpsum,"# AdjustThreshWhenOneTail %d\n",AdjustThreshWhenOneTail);

  if (distthresh > 0)
    fprintf(fpsum,"# Distance Threshold: %g (mm)\n",distthresh);

  if(cwpvalthresh > 0)
    fprintf(fpsum,"# CW PValue Threshold: %g \n",cwpvalthresh);

  fprintf(fpsum,"# Size Threshold:    %g mm^3\n",sizethresh);
  fprintf(fpsum,"# Size Threshold:    %g voxels\n",sizethresh/voxsize);
  fprintf(fpsum,"# Voxel Size:        %g mm^3\n",voxsize);
  if (regfile) fprintf(fpsum,"# Registration:      %s\n",regfile);
  else fprintf(fpsum,"# Registration:      None : Tal Coords invalid\n");
  if (synthfunction != NULL)
    fprintf(fpsum,"# Synthesize:        %s\n",synthfunction);
  if (maskid != NULL) {
    fprintf(fpsum,"# Mask Vol:          %s\n",maskid);
    fprintf(fpsum,"# Mask Thresh:       %f\n",maskthresh);
    fprintf(fpsum,"# Mask Sign:         %s\n",masksignstring);
    fprintf(fpsum,"# Mask Invert:       %d\n",maskinvert);
  }
  fprintf(fpsum,"# AllowDiag:         %d\n",allowdiag);
  fprintf(fpsum,"# NClusters          %d\n",nclusters);
  if (csd != NULL) {
    fprintf(fpsum,"# CSD thresh  %lf\n",csd->thresh);
    fprintf(fpsum,"# CSD nreps    %d\n",csd->nreps);
    fprintf(fpsum,"# CSD simtype  %s\n",csd->simtype);
    fprintf(fpsum,"# CSD contrast %s\n",csd->contrast);
    fprintf(fpsum,"# CSD confint  %lf\n",ciPct);
  }
  if (fwhm > 0) {
    fprintf(fpsum,"# FWHM        %lf\n",fwhm);
    // dLh is an FSL parameter. Include here for comparison
    fprintf(fpsum,"# dLh         %lf\n", pow((fwhm/colres)/sqrt(4.0*log(2.0)),-3) );
  }

  fprintf(fpsum,"# \n");
  if (regfile) {
    if (FixMNI) {
      fprintf(fpsum,"# Reporting Coordinates in Talairach Space\n");
      fprintf(fpsum,"# Cluster   Size(n)   Size(mm^3)     "
              "TalX   TalY    TalZ              Max");
    } else {
      fprintf(fpsum,"# Reporting Coordinates in MNI305 Space\n");
      fprintf(fpsum,"# Cluster   Size(n)   Size(mm^3)     "
              "MNIX   MNIY    MNIZ              Max");
    }
  } else {
    fprintf(fpsum,"# Reporting Coordinates in Voxel Indices\n");
    fprintf(fpsum,"# Cluster   Size(n)   Size(mm^3)     "
            "VoxX    VoxY    VoxZ             Max");
  }

  if (csd != NULL)  fprintf(fpsum,"    CWP    CWPLow    CWPHi\n");
  else if (fwhm > 0) fprintf(fpsum,"     GRFCWP\n");
  else fprintf(fpsum,"\n");

  for (n = 0; n < nclusters; n++) {
    double maxval = ClusterList[n]->maxval;
    if(sig2pmax) {
      maxval = pow(10.0,-fabs(ClusterList[n]->maxval));
      if(BonferroniMax > 1) maxval *= BonferroniMax;
    }
    clustComputeTal(ClusterList[n],CRS2MNI); /* for "true" Tal coords */
    //clustComputeXYZ(ClusterList[n],CRS2FSA); /* for FSA coords */
    //clustComputeXYZ(ClusterList[n],CRS2MNI); /* for MNI coords */
    col = ClusterList[n]->col[ClusterList[n]->maxmember];
    row = ClusterList[n]->row[ClusterList[n]->maxmember];
    slc = ClusterList[n]->slc[ClusterList[n]->maxmember];
    if (regfile) {
      x = ClusterList[n]->x[ClusterList[n]->maxmember];
      y = ClusterList[n]->y[ClusterList[n]->maxmember];
      z = ClusterList[n]->z[ClusterList[n]->maxmember];
    } else {
      x=col;
      y=row;
      z=slc;
    }
    fprintf(fpsum,"%3d        %5d      %8.1f    %7.2f %7.2f %7.2f   %15.5f",
            n+1,ClusterList[n]->nmembers,voxsize*ClusterList[n]->nmembers,
            x,y,z, maxval);
    if (debug) fprintf(fpsum,"  %3d %3d %3d \n",col,row,slc);
    if (csd != NULL)
      fprintf(fpsum,"  %7.5lf  %7.5lf  %7.5lf",
              ClusterList[n]->pval_clusterwise,
              ClusterList[n]->pval_clusterwise_low,
              ClusterList[n]->pval_clusterwise_hi);
    else if (fwhm > 0)
      fprintf(fpsum,"  %9.7lf",ClusterList[n]->pval_clusterwise);
    if(segvolfile){
      ctabindex = MRIgetVoxVal(segvol,col,row,slc,0);
      fprintf(fpsum,"  %s",segctab->entries[ctabindex]->name);
    }
    fprintf(fpsum,"\n");
  }
  if (sumfile != NULL) fclose(fpsum);

  if(pointset != NULL){
    fp = fopen(pointset,"w");
    MATRIX *vox2ras = MRIxfmCRS2XYZ(vol,0);
    MATRIX *crs = MatrixAlloc(4,1,MATRIX_REAL);
    MATRIX *xyz=NULL;
    crs->rptr[4][1] = 1;
    for (n = 0; n < nclusters; n++) {
      crs->rptr[1][1] = ClusterList[n]->col[ClusterList[n]->maxmember];
      crs->rptr[2][1] = ClusterList[n]->row[ClusterList[n]->maxmember];
      crs->rptr[3][1] = ClusterList[n]->slc[ClusterList[n]->maxmember];
      xyz = MatrixMultiply(vox2ras,crs,xyz);
      fprintf(fp,"%7.4f %7.4f %7.4f\n",xyz->rptr[1][1],xyz->rptr[2][1],xyz->rptr[3][1]);
    }
    fprintf(fp,"info\n");
    fprintf(fp,"numpoints %d\n",nclusters);
    fprintf(fp,"UseRealRAS 1\n");
    fclose(fp);
  }

  /* Write clusters values to a volume */
  if (outid != 0) {
    outvol = clustClusterList2Vol(ClusterList, nclusters, vol,frame, 1);
    //MRIwriteType(outvol,outid,outtype);
    MRIwrite(outvol,outid);
    MRIfree(&outvol);
  }

  /* --- Save the cluster pval --- */
  if (clustwisesigfile != NULL) {
    printf("Saving cluster pval %s\n",clustwisesigfile);
    clustwisesig = MRIclone(vol,NULL);
    for (n = 0; n < nclusters; n++) {
      pval = ClusterList[n]->pval_clusterwise;
      if( pval < 10e-30) pval = 1e-30;
      val = -log10(pval)*SIGN(ClusterList[n]->maxval);
      /*printf("%3d %5d %8.7lf %7.3lf  \n", n,ClusterList[n]->nmembers,
      ClusterList[n]->pval_clusterwise,val);*/
      for (m = 0; m < ClusterList[n]->nmembers; m++) {
        col = ClusterList[n]->col[m];
        row = ClusterList[n]->row[m];
        slc = ClusterList[n]->slc[m];
        MRIsetVoxVal(clustwisesig,col,row,slc,0,val);
      }
    }
    MRIwrite(clustwisesig,clustwisesigfile);
  }

  /* Write clusters numbers to a volume, include color LUT */
  if (outcnid != 0) {
    ct = CTABalloc(nclusters+1);
    strcpy(ct->entries[0]->name,"Unknown");
    for (n = 0; n < nclusters; n++)
      sprintf(ct->entries[n+1]->name,"Cluster-%03d",n+1);
    stem = IDstemFromName(outcnid);
    sprintf(tmpstr,"%s.lut",stem);
    CTABwriteFileASCII(ct, tmpstr);
    outvol = clustClusterList2Vol(ClusterList, nclusters, vol,frame, 0);
    if(outvol->ct) CTABfree(&outvol->ct);
    outvol->ct = ct;
    printf("INFO: writing OCN to %s\n",outcnid);
    MRIwrite(outvol,outcnid);
    //MRIwriteType(outvol,outcnid,outcntype);
    MRIfree(&outvol);
    CTABfree(&ct);
    free(stem);
  }

  /* Write the given cluster to a label file */
  /* Warning: the cluster xyz will change */
  if (labelfile != NULL) {
    if (nlabelcluster > nclusters) {
      fprintf(stderr,"ERROR: selected cluster number %d, "
              "but there are only %d clusters\n",nlabelcluster,nclusters+1);
      exit(1);
    }
    printf("Computing Label\n");
    label = clustCluster2Label(ClusterList[nlabelcluster-1], vol, frame,
                               colres, rowres, sliceres, FSA2Func);
    LabelWrite(label, labelfile);
  }

  if (labelbase != NULL) {
    for (nlabelcluster = 0; nlabelcluster < nclusters; nlabelcluster ++) {

      printf("Computing label for cluster %d\n",nlabelcluster);
      sprintf(tmpstr,"%s-%04d.label",labelbase,nlabelcluster+1);
      label = clustCluster2Label(ClusterList[nlabelcluster], vol, frame,
                                 colres, rowres, sliceres, FSA2Func);
      printf("Saving label to %s\n",tmpstr);
      LabelWrite(label, tmpstr);
      LabelFree(&label);
    }
  }

  // Delete temp reg file created by --regheader
  if(regheader) unlink(regfile);

  printf("mri_volcluster: done\n");

  return(0);
  exit(0);
}

/* ----------------------------------------------------------- */
/* --------------->>>>>>.<<<<<<<------------------------------ */
/* ----------------------------------------------------------- */
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
    else if (!strcasecmp(option, "--verbose")) verbose = 1;
    else if (!strcasecmp(option, "--no-adjust")) AdjustThreshWhenOneTail=0;
    else if (!strcasecmp(option, "--allowdiag")) allowdiag = 1;
    else if (!strcmp(option, "--maskinvert"))    maskinvert = 1;
    else if (!strcmp(option, "--nofixtkreg"))    fixtkreg = 0;
    else if (!strcmp(option, "--fixtkreg"))      fixtkreg = 1;
    else if (!strcmp(option, "--csdpdf-only")) csdpdfonly = 1;
    else if (!strcasecmp(option, "--no-fixmni"))  FixMNI = 0;
    else if (!strcasecmp(option, "--fixmni"))     FixMNI = 1;
    else if (!strcasecmp(option, "--fsaverage"))  UseFSAverage = 1;
    else if (!strcmp(option, "--sig2p-max")) sig2pmax = 1;

    else if (!strcasecmp(option, "--diag")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ctab")) {
      if (nargc < 1) CMDargNErr(option,1);
      ctabfile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--seg")) {
      if(nargc < 1) CMDargNErr(option,2);
      subject = pargv[0];
      segvolfile = pargv[1];
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      if (SUBJECTS_DIR == NULL) {
        printf("ERROR: SUBJECTS_DIR not defined in environment\n");
        exit(1);
      }
      sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,segvolfile);
      if (!fio_FileExistsReadable(tmpstr)) {
        printf("ERROR: cannot find %s\n",tmpstr);
        exit(1);
      }
      segvolpath = strcpyalloc(tmpstr);
      nargsused = 2;
    } else if (!strcasecmp(option, "--sd")) {
      if (nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
    }

    /* -------- source volume inputs ------ */
    else if (!strcmp(option, "--i") || !strcmp(option, "--in")) {
      if (nargc < 1) argnerr(option,1);
      volid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--in_type")) {
      if (nargc < 1) argnerr(option,1);
      intypestring = pargv[0];
      intype = string_to_type(intypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--mask")) {
      if (nargc < 1) argnerr(option,1);
      maskid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--mask_type")) {
      if (nargc < 1) argnerr(option,1);
      masktypestring = pargv[0];
      masktype = string_to_type(masktypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--maskframe")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&maskframe);
      if (maskframe < 0) {
        printf("ERROR: negative frame number: frame = %d\n",maskframe);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--masksign")) {
      if (nargc < 1) argnerr(option,1);
      masksignstring = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--maskthresh")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&maskthresh);
      if (maskthresh < 0.0) {
        printf("ERROR: negative mask threshold not"
               "allowed (use -masksign)\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--outmask")) {
      if (nargc < 1) argnerr(option,1);
      outmaskid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--outmask_type")) {
      if (nargc < 1) argnerr(option,1);
      outmasktypestring = pargv[0];
      outmasktype = string_to_type(outmasktypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--synth")) {
      if (nargc < 1) argnerr(option,1);
      synthfunction = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--o") || !strcmp(option, "--out")) {
      if (nargc < 1) argnerr(option,1);
      outid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--out_type")) {
      if (nargc < 1) argnerr(option,1);
      outtypestring = pargv[0];
      outtype = string_to_type(outtypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--cwsig") ) {
      if (nargc < 1) argnerr(option,1);
      clustwisesigfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--cwpvalthresh") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&cwpvalthresh);
      nargsused = 1;
    } else if (!strcmp(option, "--ocn")) {
      if (nargc < 1) argnerr(option,1);
      outcnid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--ocn_type")) {
      if (nargc < 1) argnerr(option,1);
      outcntypestring = pargv[0];
      outcntype = string_to_type(outcntypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--vwsig")) {
      if (nargc < 1) argnerr(option,1);
      voxwisesigfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--vwsigmax")) {
      if(nargc < 1) argnerr(option,1);
      maxvoxwisesigfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--sum")) {
      if (nargc < 1) argnerr(option,1);
      sumfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--labelfile") ||
               !strcmp(option, "--label")) {
      if (nargc < 1) argnerr(option,1);
      labelfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--labelbase") ||
               !strcmp(option, "--labelbase")) {
      if (nargc < 1) argnerr(option,1);
      labelbase = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--nlabelcluster")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nlabelcluster);
      if (nlabelcluster <= 0) {
        fprintf(stderr,"ERROR: non-postive nlabelcluster not allowed\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--thmin")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&threshmin);
      if(threshmin < 0.0) {
        printf("ERROR: %g negative threshold not allowed (use -sign)\n",threshmin);
        exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcmp(option, "--thmax")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&threshmax);
      if (threshmax < 0.0) {
        fprintf(stderr,"ERROR: negative threshold not allowed (use -sign)\n");
        exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcmp(option, "--match")) {
      int matchval;
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&matchval);
      threshmin = matchval - 0.5;
      threshmax = matchval + 0.5;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--sign")) {
      if (nargc < 1) argnerr(option,1);
      signstring = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--csd")) {
      if (nargc < 1) argnerr(option,1);
      csdfile = pargv[0];
      csd = CSDreadMerge(csdfile,csd);
      if (csd == NULL) exit(1);
      if (strcmp(csd->anattype,"volume")) {
        printf("ERROR: csd must have anattype of volume\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--csdpdf")) {
      if (nargc < 1) argnerr(option,1);
      csdpdffile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--frame")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      if (frame < 0) {
        fprintf(stderr,"ERROR: negative frame number: frame = %d\n",frame);
        exit(1);
      }
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
    else if (!strcmp(option, "--minsize")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&sizethresh);
      if (sizethresh < 0) {
        printf("ERROR: %g negative cluster size threshold not allowed\n",
               sizethresh);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--minsizevox")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&sizethreshvox);
      if (sizethreshvox <  0) {
        printf("ERROR: %d negative cluster size threshold not allowed\n",
               sizethreshvox);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--mindist")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&distthresh);
      if (distthresh <= 0) {
        fprintf(stderr,"ERROR: negative distance threshold not allowed\n");
        exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcmp(option, "--reg")) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--regheader")) {
      if (nargc < 1) argnerr(option,1);
      FILE *fp;
      srand48(PDFtodSeed());
      std::string tmpname = getTempFile(".dat");
      sprintf(tmpstr, "%s", tmpname.c_str());
      regfile = strcpyalloc(tmpstr);
      fp = fopen(regfile,"w");
      fprintf(fp,"%s\n",pargv[0]);
      fprintf(fp,"1\n");
      fprintf(fp,"1\n");
      fprintf(fp,"1\n");
      fprintf(fp,"1 0 0 0\n");
      fprintf(fp,"0 1 0 0\n");
      fprintf(fp,"0 0 1 0\n");
      fprintf(fp,"0 0 0 1\n");
      fprintf(fp,"round\n");
      fclose(fp);
      regheader = 1;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--mni152reg")) {
      sprintf(tmpstr,"%s/average/mni152.register.dat",getenv("FREESURFER_HOME"));
      regfile = strcpyalloc(tmpstr);
      nargsused = 0;
    }
    else if (!strcmp(option, "--pointset")) {
      if (nargc < 1) argnerr(option,1);
      pointset = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--fwhm")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&fwhm);
      nargsused = 1;
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
  printf("   --in infile : source of volume values\n");
  printf("\n");
  printf("   --sum file   : text summary file \n");
  printf("   --out      output volid \n");
  printf("   --ocn      output cluster number volid \n");
  printf("   --cwsig    clusterwise sig volid \n");
  printf("   --pointset pointset.dat : create a freeview pointset of the clusters\n");
  printf("\n");
  printf("   --thmin   minthresh : minimum intensity threshold\n");
  printf("   --thmax   maxthresh : maximum intensity threshold\n");
  printf("   --sign    sign      : <abs> or pos/neg for one-sided tests\n");
  printf("   --no-adjust  : do not adjust thresh for one-tailed tests (assumes z for GRF)\n");
  printf("   --match matchval : set thmin=matchval-0.5 and thmax=matchval+0.5\n");
  printf("\n");
  printf("   --cwpvalthresh pval : require clusters to have cwp < thresh\n");
  printf("      with --fwhm or --csd\n");
  printf("\n");
  printf("   --reg     register.dat : for reporting talairach coords\n");
  printf("   --mni152reg : input is in mni152 space\n");
  printf("   --regheader subject : use header registration with subject\n");
  printf("   --fsaverage : assume input is in fsaverage space\n");
  printf("   --frame   frameno <0>\n");
  printf("\n");
  printf("   --csd csdfile <--csd csdfile ...>\n");
  printf("   --cwsig cwsig : map of corrected cluster-wise significances\n");
  printf("   --vwsig vwsig : map of corrected voxel-wise significances\n");
  printf("   --csdpdf csdpdffile : PDF/CDF of cluster and max sig\n");
  printf("   --csdpdf-only : write csd pdf file and exit.\n");
  printf("\n");
  printf("   --fwhm fwhm : fwhm in mm3, forces GRF analysis\n");
  printf("   --fwhmdat fwhm.dat : text file with fwhm in mm3 for GRF\n");
  printf("\n");
  printf("   --minsize    minimum volume (mm^3)\n");
  printf("   --minsizevox minimum volume (voxels)\n");
  printf("   --mindist distance threshold <0>\n");
  printf("   --allowdiag  : define contiguity to include diagonal\n");
  printf("   --bonferroni N : addition correction across N (eg, spaces)\n");
  printf("   --bonferroni-max N : apply bonf cor to maximum (only applies with --sig2p-max)\n");
  printf("   --sig2p-max : convert max from sig to p\n");
  printf("\n");
  printf("   --mask      mask volid (same dim as input)\n");
  printf("   --mask_type file format \n");
  printf("   --maskframe frameno <0> \n");
  printf("   --maskthresh  upper threshold \n");
  printf("   --masksign   <abs>, neg, pos \n");
  printf("   --maskinvert  \n");
  printf("   --outmask      final binary mask\n");
  printf("   --outmask_type file format \n");
  printf("\n");
  printf("\n");
  printf("   --label   label file\n");
  printf("   --nlabelcluster n : save nth cluster in label\n");
  printf("   --labelbase  base : save all clusters under base-NNNN.label\n");
  printf("\n");
  printf("   --synth   synthfunc (uniform,loguniform,gaussian)\n");
  printf("   --diag diagno : set diagnostic level\n");
  printf("   --help    : how to use this program \n");
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  printf("\n");
  print_usage() ;
  printf("\n");

  printf
  (

    "DESCRIPTION:\n"
    "\n"
    "This program will find clusters in a volume. A cluster is a set of\n"
    "contiguous voxels which meet a threshold criteria. The set may also\n"
    "have to reach a certain minimum number of voxels to be considered a\n"
    "cluster. The results can be saved in four ways: (1) a text summary\n"
    "file, (2) a new volume which is same as the input volume but with all\n"
    "the voxels that do not belong to a cluster set to zero, (3) a volume\n"
    "with each voxel's value equal to the cluster number in the summary\n"
    "file to which the voxel belongs, and (4) one cluster can be saved as a\n"
    "label file. The search space within the volume can be restricted to be\n"
    "within a mask. Two voxels are considered contiguous if they share a \n"
    "common row, column, or slice (except for --allowdiag).\n"
    "\n"
    "COMMAND-LINE ARGUMENTS:\n"
    "\n"
    "  --in input-volid : path/name of the input volume to be clustered.\n"
    "\n"
    "  --in_type typename : name of the file format type of the input\n"
    "    volume. See FILE TYPES below.\n"
    "\n"
    "  --frame n : perform the cluster analysis on the nth frame. The first\n"
    "    frame corresponds to n=0.\n"
    "\n"
    "  --reg registration-file : registration file created by either tkregister\n"
    "    or tkmedit. Note: the subject listed in the registration file must \n"
    "    have an entry under $SUBJECTS_DIR.\n"
    "\n"
    "  --thmin minthresh : voxel values must exceed this amount to be considered\n"
    "    as a candidate for a cluster.\n"
    "\n"
    "  --thmax maxthresh : the value of a voxel must be less than this\n"
    "    amount to be considered as a candidate for a cluster. Default is\n"
    "    infinity.\n"
    "\n"
    "  --sign sign-name : the value of a voxel must be this sign to be considered\n"
    "    as a candidate for a cluster. pos = positive, neg = negative, abs = \n"
    "    absolute (ie, ignore sign). Default is abs. When neg is used, the \n"
    "    interpretation of minthreshold and maxthreshold are reversed.\n"
    "--csd csdfile <--csd csdfile>\n"
    "\n"
    "    Load one or more CSD files. CSD stands for 'Cluster Simulation Data'. This\n"
    "    file is produced by running mri_glmfit with --sim. The the threshold and hemi\n"
    "    info are obtained from the CSD file and cannot be specified on the command-\n"
    "    line. If more than one CSD is specified, they are merged into one CSD internally.\n"
    "    When a CSD file is specified, three more columns are added to the summary table:\n"
    "      1. CWP - cluster-wise pvalue. The pvalue of the cluster corrected for \n"
    "         multiple comparisons\n"
    "      2. CWPLow - lower 90%% confidence limit of CWP based on binomial distribution\n"
    "      3. CWPHi  - upper 90%% confidence limit of CWP based on binomial distribution\n"
    "    In addition, the user can specify --ocp, which saves the sigificance map of \n"
    "    the clusters in which the value of each voxel is the -log10(pvalue) of cluster\n"
    "    to which the vertex belongs.\n"
    "\n"
    "--csdpdf csdpdfile\n"
    "\n"
    "    Compute PDF/CDF of CSD data and save in csdpdffile. This is mostly good for debugging.\n"
    "\n"
    "--csdpdf-only\n"
    "\n"
    "    Only write the csd pdf file.\n"
    "\n"
    "  --minsize volume : minimum volume (in mm^3) of contiguous voxels that\n"
    "    meet the threshold criteria needed to become a cluster. See also\n"
    "    --minsizevox.\n"
    "\n"
    "  --minsizevox number : minimum number of contiguous voxels that meet\n"
    "    the threshold criteria needed to become a cluster. See also --minsize.\n"
    "\n"
    "  --mindistance distance : minimum distance (in mm) that the peaks of two\n"
    "    clusters must be separated in order for them to be considered two\n"
    "    distinct clusters. If two clusters are closer than this amount, the\n"
    "    cluster with the lower peak is eliminated.\n"
    "\n"
    "  --allowdiag : (no argument) allow two voxels that share a corner to \n"
    "    be considered contiguous.\n"
    "\n"
    "  --mask mask-volid: path/name of a mask used to restrict the region\n"
    "    over which clusters will be searched. For example, this could be used\n"
    "    to restrict the search space to only the brain (ie, exclude eyeballs).\n"
    "    The mask must have the same dimension as the input volume.\n"
    "\n"
    "  --mask_type typename : name of the file format type of the mask\n"
    "    volume. See FILE TYPES below.\n"
    "\n"
    "  --maskframe n : use the nth frame of the mask volume as the mask,\n"
    "    where n = 0 indicates the first frame. Default is 0.\n"
    "\n"
    "  --maskthresh maskthresh: use only those voxels in the mask whose\n"
    "    value exceeds the given threshold (see also --masksign and\n"
    "    --maskinvert).\n"
    "\n"
    "  --masksign sign-name : threshold the mask based on the sign of the\n"
    "    value at each mask voxel: pos = positive, neg = negative, abs = \n"
    "    absolute (ie, ignore sign). Default is abs. When neg is used, the \n"
    "    a mask voxel value must be less than -maskthresh.\n"
    "\n"
    "  --maskinverse : after determining which voxels are in the mask based\n"
    "    on the threshold and sign, take only voxels that are outside of the\n"
    "    mask.\n"
    "\n"
    "  --outmask outmask-volid: path/name of the final binary mask after\n"
    "    considering thresholding, sign, and inversion. This is mainly useful\n"
    "    for debuggin.\n"
    "\n"
    "  --outmask_type typename : name of the file format type of the outmask\n"
    "    volume. See FILE TYPES below.\n"
    "\n"
    "  --sumfile sumfilename : save a summary of the results in ASCII into \n"
    "    sumfilename. See SUMMARY FILE below.\n"
    "\n"
    "  --out out-volid: path/name of the output volume after\n"
    "    clustering. All voxels that were not part of a cluster have their\n"
    "    values set to zero. Otherwise, their values do not change.\n"
    "\n"
    "  --out_type typename : name of the file format type of the output \n"
    "    volume. See FILE TYPES below.\n"
    "\n"
    "  --ocn ocn-volid: path/name of the output volume after clustering\n"
    "    where the value at each voxel is the cluster number (as found in the\n"
    "    summary file). Voxels that did not belong to a cluster have their\n"
    "    values set to zero. Also creates ocn-volid.lut. This is a color\n"
    "    lookup table so that the OCN can be loaded as a segmentation in\n"
    "    tkmedit with -segmentation ocn-volid.fmt ocn-volid.lut.\n"
    "\n"
    "  --ocn_type typename : name of the file format type of the output \n"
    "    cluster number volume. See FILE TYPES below.\n"
    "\n"
    "  --label label-file : save the nth cluster (see -nlabelcluster) as a\n"
    "    label file in the subject's anatomical space. Note: make sure that the\n"
    "    label file includes a forward slash (/) or the label will be saved\n"
    "    into the subjects anatomical direcotry. For example: ./mylabel.label.\n"
    "    Requires --nlabelcluster.\n"
    "\n"
    "  --nlabelcluster n : save the nth cluster (see -label) as a label file.\n"
    "\n"
    "  --labelbase base : save each cluster in its own label file. The name\n"
    "    of the file will be base-NNNN.label, where NNNN is the four digit,\n"
    "    zero-padded cluster number. All clusters found will be saved.\n"
    "\n"
    "SUMMARY FILE\n"
    "\n"
    "The summary file produced by --sumfile is in ASCII text. The summary\n"
    "file will have a short header indicating the conditions underwhich the\n"
    "clustering was performed followed by rows of data with 7 columns. Each\n"
    "row reprsents a different cluster. Each column has the following\n"
    "interpretation: (1) cluster number, (2) number of voxels in the\n"
    "cluster, (3) cluster volume in mm^3, (4-6) Talairach X, Y, and Z (mm),\n"
    "(7) maximum value inside the cluster. The Talairach coordinates are\n"
    "the 'true' coordinates (not MNI). Part of a summary file is \n"
    "shown below as an example:\n"
    "\n"
    "-----------------------------------------------------------------------\n"
    "Cluster Growing Summary (mri_volcluster)\n"
    "Input Volume:      grp-subj-p1p2vlor/bold/sxa-h20/tal-ffx-rfx/novvfix/sig\n"
    "Frame Number:      0\n"
    "Minimum Threshold: 2\n"
    "Maximum Threshold: inifinity\n"
    "Threshold Sign:    abs\n"
    "Distance Threshold: 10 (mm)\n"
    "Size Threshold:    640 mm^3\n"
    "Size Threshold:    10 voxels\n"
    "Voxel Size:        64 mm^3\n"
    "Registration:      grp-subj-p1p2vlor/bold/sxa-h20/tal-ffx-rfx/register.dat\n"
    "Mask Vol:          talbrain/mask\n"
    "Mask Thresh:       0.500000\n"
    "Mask Sign:         pos\n"
    "Mask Invert:       0\n"
    "AllowDiag:         1\n"
    "NClusters          26\n"
    "\n"
    "Cluster Size(n) Size(mm^3)  TalX   TalY    TalZ     Max\n"
    "  1      348      22272.0   59.40 -66.72  -13.48   5.66192\n"
    "  2       45       2880.0  -39.60  26.79   -8.07   5.45487\n"
    "  3       27       1728.0   55.44  16.60   21.28   4.95684\n"
    "-----------------------------------------------------------------------\n"
    "\n"
    "\n"
    "FILE TYPES/FORMATS:\n"
    "\n"
    "mri_volcluster can read/write any file format that can be read/written\n"
    "by mri_convert. The most relevent ones are: bshort, bfloat, analyze,\n"
    "analyze4d. When specifying bshort/bfloat, the volume id is the\n"
    "stem. Ie, if the volume is f_000.bshort, f_001.bshort, ..., then the\n"
    "volume id is 'f' (no quotes).\n"
    "\n"
    "KNOWN BUGS:\n"
    "\n"
    "When specifying a label file, make sure that the label file includes a\n"
    "forward slash (/) or the label will be saved into the subjects\n"
    "anatomical direcotry. For example: ./mylabel.label.\n"
    "\n"
    "BUG REPORTS:\n"
    "\n"
    "Send bug reports to analysis-bugs@nmr.mgh.harvard.edu. Make sure to\n"
    "include the version number (--version) , the full command-line used,\n"
    "the type of computer operating system, anything printed to the screen,\n"
    "a description of what you think is wrong and why, and anything else\n"
    "that may be helpful. Users at the NMR Center should also include the\n"
    "directory that they ran it from. NOTE: bug reports that do not have\n"
    "sufficient information to diagnose the problem will probably be either\n"
    "ignored or placed at the bottom of the list. BE CLEAR!\n"
    "\n"
    "AUTHOR:\n"
    "\n"
    "Douglas N. Greve, Ph.D; NMR Center, MGH, greve@nmr.mgh.harvard.edu\n"
  ) ;

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
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
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/* --------------------------------------------- */
static void check_options(void) {
  int err;

  err = 0;

  if (csdpdffile) {
    if (csd == NULL) {
      printf("ERROR: need --csd with --csdpdf");
      exit(1);
    }
    CSDpdf(csd,-1);
    CSDwritePDF(csdpdffile,csd);
    if (csdpdfonly) exit(0);
  }
  if(voxwisesigfile != NULL && csd == NULL) {
    printf("ERROR: need csd with --vwsig\n");
    exit(1);
  }
  if (clustwisesigfile != NULL && csd == NULL && fwhm < 0) {
    printf("ERROR: need csd with --cwsig\n");
    exit(1);
  }

  if (volid == NULL) {
    fprintf(stderr,"ERROR: no input volume supplied\n");
    err = 1;
  }

  // Check cluster data file
  if (csd != NULL) {
    if (threshmin < 0) threshmin = csd->thresh;
    else {
      if (threshmin != csd->thresh) {
        printf("ERROR: you have specified thmin=%f on cmdline, but\n",
               threshmin);
        printf("CSD file was created with %lf\n",csd->thresh);
        exit(1);
      }
    }
    if (sizethresh > 0 || sizethreshvox > 0) {
      printf("ERROR: you cannot specify --minsize or "
             "--minsizevox with --csd\n");
      exit(1);
    }
    if (threshmax > 0) {
      printf("ERROR: you cannot specify --thmax  with --csd\n");
      exit(1);
    }
    threshsign = csd->threshsign;
    if (threshsign ==  0) signstring = "abs";
    if (threshsign == +1) signstring = "pos";
    if (threshsign == -1) signstring = "neg";
  } // end csd != NULL
  else {
    if (      !strncmp(signstring,"a",1) ) threshsign = 0;
    else if ( !strncmp(signstring,"p",1) ) threshsign = +1;
    else if ( !strncmp(signstring,"n",1) ) threshsign = -1;
    else {
      fprintf(stderr,"ERROR: sign = %s, must be neg, abs, or pos\n",
              signstring);
      err = 1;
    }
  }

  if (threshmin < 0) {
    fprintf(stderr,"ERROR: no minimum threshold supplied\n");
    err = 1;
  }

  if(threshsign != 0 && AdjustThreshWhenOneTail) {
    // user has requested a tailed threshold, so adjust threshold
    // to account for a one-tailed test. This requires that the
    // input be -log10(p), where p is computed from a two-tailed
    // test. One could recompute the p-values to convert to a 
    // one-sided test, but easier to just adjust the threshold
    printf("Adjusting threshold for 1-tailed test.\n");
    printf("If the input is not a -log10(p) volume, "
           "re-run with --no-adjust.\n");
    threshminadj = threshmin - log10(2.0);
    threshmaxadj = threshmax - log10(2.0);
    printf("New threshold is %lf\n",threshminadj);
  } else {
    printf("NOT Adjusting threshold for 1-tailed test\n");
    threshminadj = threshmin;
    threshmaxadj = threshmax;
  }
  printf("threshmin %g, threshminadj %g\n",threshmin,threshminadj);

  if (synthfunction != NULL) {
    if (strcmp(synthfunction,"uniform")    &&
        strcmp(synthfunction,"loguniform") &&
        strcmp(synthfunction,"gaussian") ) {
      fprintf(stderr,
              "ERROR: synth = %s, must be uniform, loguniform, or gaussian\n",
              synthfunction);
      err = 1;
    }
  }

  if (intype == MRI_VOLUME_TYPE_UNKNOWN) intype = mri_identify(volid);
  if (intype == MRI_VOLUME_TYPE_UNKNOWN) {
    fprintf(stderr,"ERROR: could not determine type of %s\n",volid);
    exit(1);
  }

  if (outid != 0) {
    if (outtype == MRI_VOLUME_TYPE_UNKNOWN) outtype = mri_identify(outid);
    if (outtype == MRI_VOLUME_TYPE_UNKNOWN) {
      fprintf(stderr,"ERROR: could not determine type of %s\n",outid);
      exit(1);
    }
  }

  if (maskid != 0) {
    if (masktype == MRI_VOLUME_TYPE_UNKNOWN) masktype = mri_identify(maskid);
    if (masktype == MRI_VOLUME_TYPE_UNKNOWN) {
      fprintf(stderr,"ERROR: could not determine type of %s\n",maskid);
      exit(1);
    }
  }

  if (outmaskid != 0) {
    if (outmasktype == MRI_VOLUME_TYPE_UNKNOWN)
      outmasktype = mri_identify(outmaskid);
    if (outmasktype == MRI_VOLUME_TYPE_UNKNOWN) {
      fprintf(stderr,"ERROR: could not determine type of %s\n",outmaskid);
      exit(1);
    }
  }

  if (outcnid != 0) {
    if (outcntype == MRI_VOLUME_TYPE_UNKNOWN)
      outcntype = mri_identify(outcnid);
    if (outcntype == MRI_VOLUME_TYPE_UNKNOWN) {
      fprintf(stderr,"ERROR: could not determine type of %s\n",outcnid);
      exit(1);
    }
  }

  if (labelfile != NULL && nlabelcluster < 1) {
    printf("ERROR: --nlabelcluster must be specified with --label\n");
    err = 1;
  }

  if(segvolfile != NULL){
    if(ctabfile == NULL){
      ctabfile = (char *) calloc(sizeof(char),1000);
      sprintf(ctabfile,"%s/FreeSurferColorLUT.txt",getenv("FREESURFER_HOME"));
      printf("Using defalt ctab %s\n",ctabfile);
    }
    ctab = CTABreadASCII(ctabfile);
    if(ctab == NULL) exit(1);
  }

  if(UseFSAverage && regfile != NULL){
    printf("ERROR: cannot --reg and --fsaverage\n");
    exit(1);
  }

  // This no longer applies
  //if(fwhm > 0 && !strcmp(signstring,"abs")){
  //printf("ERROR: you must specify a pos or neg sign with --fwhm\n");
  //exit(1);
  //}


  if (err) exit(1);
  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"volid    %s\n",volid);
  if (intypestring) fprintf(fp,"voltype  %s\n",intypestring);
  fprintf(fp,"frame    %d\n",frame);
  fprintf(fp,"thmin    %g\n",threshmin);
  fprintf(fp,"thmax    %g\n",threshmax);
  fprintf(fp,"sign     %s\n",signstring);
  fprintf(fp,"sizethresh  %g\n",sizethresh);
  fprintf(fp,"distthresh  %g\n",distthresh);
  fprintf(fp,"fixmni  %d\n",FixMNI);
}

/* ------------------------------------------------------------------
   LoadMNITransform() - returns the matrix that transforms col, row,
   and slice in the input volume to MNI x, y, z.  It also returns the
   matrix that transforms col, row, and slice in the input volume to
   FreeSurfer Anatomical x, y, z.
   --------------------------------------------------------------- */
static MATRIX *LoadMNITransform(char *regfile, int ncols, int nrows,
                                int nslices,
                                MATRIX **ppCRS2FSA, MATRIX **ppFSA2Func,
                                float *colres, float *rowres, float *sliceres) {
  extern int fixtkreg;
  extern MRI *vol;
  int float2int;
  char *SUBJECTS_DIR;
  int err;
  char *subject;
  float ipr, bpr, intensity;
  MATRIX *Rtmp, *R, *iR, *T, *iQ;
  MATRIX *CRS2MNI;
  //char talxfmfile[1000];

  err = regio_read_register(regfile, &subject, &ipr, &bpr,
                            &intensity, &R, &float2int);
  if (err) exit(1);
  iR = MatrixInverse(R,NULL);

  if ( (fabs(vol->xsize - ipr) > .001) || fabs(vol->zsize - bpr) > .001) {
    printf("ERROR: Input volume voxel dimensions do not match those \n"
           "in the registration file. If the input volume is in \n"
           "bshort/bfloat format, check that there is an accompanying \n"
           "bhdr file.\n");
    exit(1);
  }

  if (fixtkreg && (float2int == FLT2INT_TKREG)) {
    printf("INFO: Fixing tkregister matrix\n");
    printf("Original Reg Matrix: ----------------\n");
    MatrixPrint(stdout,R);
    Rtmp = MRIfixTkReg(vol,R);
    MatrixFree(&R);
    R = Rtmp;
    printf("New Reg Matrix: ----------------\n");
    MatrixPrint(stdout,R);
  }

  /* get the SUBJECTS_DIR environment variable */
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL) {
    printf("ERROR: environment variable SUBJECTS_DIR undefined "
           "(use setenv)\n");
    exit(1);
  }
  printf("INFO: subject = %s\n",subject);

  /* Load the talairach.xfm */
  T = DevolveXFM(subject, NULL, NULL);
  if (T==NULL) exit(1);
  //sprintf(talxfmfile,"%s/%s/mri/transforms/talairach.xfm",
  //SUBJECTS_DIR,subject);
  //err = regio_read_mincxfm(talxfmfile, &T);
  //if(err) exit(1);

  iQ =   MRIxfmCRS2XYZtkreg(vol);
  printf("Input volume FOV xfm Matrix: ----------------\n");
  MatrixPrint(stdout,iQ);

  *ppCRS2FSA = MatrixMultiply(iR,iQ,NULL);
  CRS2MNI = MatrixMultiply(T,*ppCRS2FSA,NULL);

  MatrixFree(&iR);
  MatrixFree(&iQ);
  MatrixFree(&T);

  //MatrixFree(&R);
  *ppFSA2Func = R;
  *colres   = ipr;
  *rowres   = ipr;
  *sliceres = bpr;

  return(CRS2MNI);
}
/*---------------------------------------------------------------
  ---------------------------------------------------------------*/
static MRI *MRIsynthUniform(int ncols, int nrows, int nslices,
                            int nframes, MRI *tvol) {
  MRI *vol;
  int col, row, slc, frm;

  if (tvol == NULL) {
    vol = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (vol == NULL) {
      fprintf(stderr,"ERROR: MRIsynthUniform: could not alloc mri\n");
      return(NULL);
    }
  } else vol = tvol;

  for (col = 0; col < vol->width;   col++) {
    for (row = 0; row < vol->height;  row++) {
      for (slc = 0; slc < vol->depth;   slc++) {
        for (frm = 0; frm < vol->nframes; frm++) {
          MRIFseq_vox(vol,col,row,slc,frm) = (float) drand48();
        }
      }
    }
  }

  return(vol);
}
/*---------------------------------------------------------------
  ---------------------------------------------------------------*/
static MRI *MRIsynthLogUniform(int ncols, int nrows, int nslices,
                               int nframes, MRI *tvol) {
  MRI *vol;
  int col, row, slc, frm;

  if (tvol == NULL) {
    vol = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (vol == NULL) {
      fprintf(stderr,"ERROR: MRIsynthLogUniform: could not alloc mri\n");
      return(NULL);
    }
  } else vol = tvol;

  for (col = 0; col < vol->width;   col++) {
    for (row = 0; row < vol->height;  row++) {
      for (slc = 0; slc < vol->depth;   slc++) {
        for (frm = 0; frm < vol->nframes; frm++) {
          MRIFseq_vox(vol,col,row,slc,frm) = (float)(-log10(drand48()));
        }
      }
    }
  }

  return(vol);
}
/*---------------------------------------------------------------
  ---------------------------------------------------------------*/
static MRI *MRIsynthGaussian(int ncols, int nrows, int nslices,
                             int nframes, MRI *tvol) {
  MRI *vol;
  int col, row, slc, frm;

  if (tvol == NULL) {
    vol = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (vol == NULL) {
      fprintf(stderr,"ERROR: MRIsynthGaussian: could not alloc mri\n");
      return(NULL);
    }
  } else vol = tvol;

  for (col = 0; col < vol->width;   col++) {
    for (row = 0; row < vol->height;  row++) {
      for (slc = 0; slc < vol->depth;   slc++) {
        for (frm = 0; frm < vol->nframes; frm++) {
          MRIFseq_vox(vol,col,row,slc,frm) = (float)Gaussian01PDF();
        }
      }
    }
  }

  return(vol);
}

/*********************************************************
 * Name:    Gaussian01PDF()
 * Purpose: generates random numbers that obey a gaussian
 *          distribution with zero mean and std dev of 1:
 *              pdf(x) = e^(x^2/2)/sqrt(2pi)
 ************************************************************/
static double Gaussian01PDF(void) {
  double v1,v2,r2;

  do {
    v1 = 2.0 * drand48() - 1.0;
    v2 = 2.0 * drand48() - 1.0;
    r2 = v1*v1 + v2*v2;
  } while ( r2 > 1.0);

  return( v1 * sqrt( -2.0 * log(r2)/r2 ));
}

/*---------------------------------------------------------------
  MRIbinarize01() - returns a binarized volume of type int. The
  volume is binarized based on on the input volume according to
  the following rules:
  1. The sign of the input voxel is changed according to thsign.
  2. Values between thmin and thmax are binarized to highval,
  otherwise they are set to lowval.
  3. If thmax < 0, then thmax = infinity
  4. If invert=1, the lowval and highval are exchanged.
  Notes:
  1. thsign should be abs, pos, or neg
  2. If binvol is non-NULL, it's type must be MRI_INT and it
  must have the same dimensions as vol.
  ---------------------------------------------------------------*/
static MRI *MRIbinarize01(MRI *vol, float thmin, float thmax,
                          const char *thsign, int invert,
                          int lowval, int highval, int *nhits, MRI *binvol) {
  int ncols, nrows, nslices, nframes;
  int col, row, slice, frame;
  int ithsign;
  short r;
  float val=0.0;

  if (strncasecmp(thsign,"neg",3)==0)      ithsign = -1;
  else if (strncasecmp(thsign,"abs",3)==0) ithsign =  0;
  else if (strncasecmp(thsign,"pos",3)==0) ithsign = +1;
  else {
    printf("ERROR:  MRIbinarize01(): sign string = %s\n",thsign);
    return(NULL);
  }

  ncols   = vol->width;
  nrows   = vol->height;
  nslices = vol->depth;
  nframes = vol->nframes;

  if (binvol == NULL) {
    binvol = MRIallocSequence(ncols, nrows, nslices, MRI_INT, nframes);
    if (binvol == NULL) {
      printf("ERROR: MRIbinarize01(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(vol,binvol);
  } else {
    if ( (binvol->width != ncols) || (binvol->height != nrows) ||
         (binvol->depth != nslices) || (binvol->nframes != nframes) ) {
      printf("ERROR: MRIbinarize01(): dimension missmatch\n");
      return(NULL);
    }
    if (binvol->type != MRI_INT) {
      printf("ERROR: MRIbinarize01(): passed binvol "
             "type = %d, must be int (%d)\n",binvol->type,MRI_INT);
      return(NULL);
    }
  }

  *nhits = 0;
  for (col = 0; col < vol->width; col++) {
    for (row = 0; row < vol->height; row++) {
      for (slice = 0; slice < vol->depth; slice++) {
        for (frame = 0; frame < vol->nframes; frame++) {

          switch ( vol->type ) {
          case MRI_UCHAR:
            val = (float)(MRISCseq_vox(vol,col,row,slice,frame));
            break;
          case MRI_INT:
            val = (float)(MRIIseq_vox(vol,col,row,slice,frame));
            break;
          case MRI_LONG:
            val = (float)(MRILseq_vox(vol,col,row,slice,frame));
            break;
          case MRI_FLOAT:
            val = (float)(MRIFseq_vox(vol,col,row,slice,frame));
            break;
          case MRI_SHORT:
            val = (float)(MRISseq_vox(vol,col,row,slice,frame));
            break;
          }

          if (ithsign ==  0) val = fabs(val);
          if (ithsign == -1) val = -val;

          r = 0;
          if (thmax > 0) {
            if (val >= thmin && val <= thmax) r = 1;
          } else {
            if (val >= thmin ) r = 1;
          }

          if (invert) r = !r;

          if (r) MRIIseq_vox(binvol,col,row,slice,frame) = (short)highval;
          else  MRIIseq_vox(binvol,col,row,slice,frame) = (short)lowval;

          if(r) (*((int*)nhits))++;
        }
      }
    }
  }

  printf("INFO: MRIbinarize01(): nhits = %d\n",*nhits);

  return(binvol);
}

