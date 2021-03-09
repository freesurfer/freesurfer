/*
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
BEGINHELP --------------------------------------------------------------

This program computes tables for performing multiple comparisons on
the surface based on a monte carlo simulation using white gaussian
noise smoothed on the surface. mri_mcsim was used to generate the
tables found in $FREESURFER_HOME/average/mult-comp-cor and used by
mri_glmfit-sim with the --cache option. The tables in mult-comp-cor
are for fsaverage and fsaverage_sym for a search space over the whole
cortex. mri_mcsim can be used to create new tables for a new subject
or for fsaverage/fsaverage_sum with a smaller search space.

Example 1: Create tables for a new subject for whole hemisphere

mri_mcsim --o /path/to/mult-comp-cor/newsubject/lh/cortex --base mc-z 
  --save-iter  --surf newsubject lh --nreps 10000

This may take hours (or even days) to run; see below for
parallelizing.  When running mri_glmfit-sim, add --cache-dir
/path/to/mult-comp-cor

Example 2: Create tables for superior temporal gyrus for fsaverage

First, create the label by running

mri_annotation2label --subject fsaverage --hemi lh --outdir labeldir

mri_mcsim --o /path/to/mult-comp-cor/fsaverage/lh/superiortemporal --base mc-z 
  --save-iter  --surf fsaverage lh --nreps 10000 
  --label labeldir/lh.superiortemporal.label

When running mri_glmfit, make sure to use   --label labeldir/lh.superiortemporal.label
When running mri_glmfit-sim, add --cache-dir /path/to/mult-comp-cor --cache-label superiortemporal

Example 3: running simulations in parallel (two jobs, 5000 iterations
each for a total of 10000)

mri_mcsim --o /path/to/mult-comp-cor/fsaverage/lh/superiortemporal --base mc-z.j001 
  --save-iter  --surf fsaverage lh --nreps 5000
  --label labeldir/lh.superiortemporal.label

mri_mcsim --o /path/to/mult-comp-cor/fsaverage/lh/superiortemporal --base mc-z.j002 
  --save-iter  --surf fsaverage lh --nreps 500
  --label labeldir/lh.superiortemporal.label

When those jobs are done, merge the results into a single table with

mri_surfcluster 
  --csd /path/to/mult-comp-cor/fsaverage/lh/superiortemporal/fwhm10/abs/th20/mc-z.j001.csd 
  --csd /path/to/mult-comp-cor/fsaverage/lh/superiortemporal/fwhm10/abs/th20/mc-z.j002.csd 
  --csd-out /path/to/mult-comp-cor/fsaverage/lh/superiortemporal/fwhm10/abs/th20/mc-z.csd 
  --csdpdf  /path/to/mult-comp-cor/fsaverage/lh/superiortemporal/fwhm10/abs/th20/mc-z.cdf --csdpdf-only

Repeat the above command for each FWHM, sign (pos, neg, abs) and threshold

ENDHELP --------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "matrix.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "randomfields.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int SaveOutput(void);
int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *OutTop = NULL;
const char *csdbase = NULL;
const char *subject = NULL;
const char *surfname = "white";
const char *hemi = NULL;
char *MaskFile = NULL;
const char *LabelFile = "cortex";
int nRepetitions = -1;
int SynthSeed = -1;

int nThreshList;
double ThreshList[100];
int nFWHMList;
double FWHMList[100];
int SignList[3] = {-1,0,1}, nSignList=3;
char *DoneFile = NULL;
char *LogFile = NULL;
char *StopFile = NULL;
char *SaveFile = NULL;
int SaveMask = 1;
int UseAvgVtxArea = 0;
int SaveEachIter = 0;

CSD *csdList[100][100][3], *csd;
MRI *mask=NULL;
MRIS *surf;
char tmpstr[2000];
const char *signstr=NULL;
int msecTime, nmask, nmaskout, nthRep;
int *nSmoothsList;
double fwhmmax=30;
int SaveWeight=0;
int FixFSALH = 1;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, n, err,k, *maskoutvtxno;
  char tmpstr[2000], *SUBJECTS_DIR, fname[2000];
  const char *signstr = NULL; // Is this intended to mask the global?
  //char *OutDir = NULL;
  RFS *rfs;
  int nSmoothsPrev, nSmoothsDelta;
  MRI *z, *zabs=NULL, *sig=NULL, *p=NULL;
  int FreeMask = 0;
  int nthSign, nthFWHM, nthThresh;
  double sigmax, zmax, threshadj, csize, csizeavg, cweightvtx, searchspace,avgvtxarea;
  int csizen;
  int nClusters, cmax,rmax,smax;
  SURFCLUSTERSUM *SurfClustList;
  Timer mytimer;
  LABEL *clabel;
  FILE *fp, *fpLog=NULL;
  float **ppVal, **ppSig, **ppVal0, **ppSig0, **ppZ, **ppZ0;

  nargs = handleVersionOption(argc, argv, "mri_mcsim");
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

  if(LogFile){
    fpLog = fopen(LogFile,"w");
    if(fpLog == NULL){
      printf("ERROR: opening %s\n",LogFile);
      exit(1);
    }
    dump_options(fpLog);
  } 

  if(SynthSeed < 0) SynthSeed = PDFtodSeed();
  srand48(SynthSeed);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  // Create output directory
  printf("Creating %s\n",OutTop);
  err = fio_mkdirp(OutTop,0777);
  if(err) exit(1);
  for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++){
    for(nthThresh = 0; nthThresh < nThreshList; nthThresh++){
      for(nthSign = 0; nthSign < nSignList; nthSign++){
	if(SignList[nthSign] ==  0) signstr = "abs"; 
	if(SignList[nthSign] == +1) signstr = "pos"; 
	if(SignList[nthSign] == -1) signstr = "neg"; 
	sprintf(tmpstr,"%s/fwhm%02d/%s/th%02d",
		OutTop,(int)round(FWHMList[nthFWHM]),
		signstr,(int)round(10*ThreshList[nthThresh]));
	sprintf(fname,"%s/%s.csd",tmpstr,csdbase);
	if(fio_FileExistsReadable(fname)){
	  printf("ERROR: output file %s exists\n",fname);
	  if(fpLog) fprintf(fpLog,"ERROR: output file %s exists\n",fname);
          exit(1);
	}
	err = fio_mkdirp(tmpstr,0777);
	if(err) exit(1);
      }
    }
  }

  // Load the target surface
  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
  printf("Loading %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) return(1);

  printf("group_avg_surface_area %g\n",surf->group_avg_surface_area);
  printf("group_avg_vtxarea_loaded %d\n",surf->group_avg_vtxarea_loaded);
  if(fpLog){
    fprintf(fpLog,"group_avg_surface_area %g\n",surf->group_avg_surface_area);
    fprintf(fpLog,"group_avg_vtxarea_loaded %d\n",surf->group_avg_vtxarea_loaded);
  }

  // Handle masking
  if(LabelFile){
    printf("Loading label file %s\n",LabelFile);
    sprintf(tmpstr,"%s/%s/label/%s.%s.label",
	    SUBJECTS_DIR,subject,hemi,LabelFile);
    if(!fio_FileExistsReadable(tmpstr)){
      printf(" Cannot find label file %s\n",tmpstr);
      sprintf(tmpstr,"%s",LabelFile);
      printf(" Trying label file %s\n",tmpstr);
      if(!fio_FileExistsReadable(tmpstr)){
	printf("  ERROR: cannot read or find label file %s\n",LabelFile);
	exit(1);
      }
    }
    printf("Loading %s\n",tmpstr);
    clabel = LabelRead(NULL, tmpstr);
    mask = MRISlabel2Mask(surf, clabel, NULL);
    FreeMask = 1;
  }
  if(MaskFile){
    printf("Loading %s\n",MaskFile);
    mask = MRIread(MaskFile);
    if(mask == NULL) exit(1);
  }
  if(mask && strcmp(subject,"fsaverage")==0 && strcmp(hemi,"lh")==0 && FixFSALH){
    // This is a hack to delete two stray vertices in the fsaverage lh.cortex.label
    // file. They could be deleted in the label file. They have no effect on cluster
    // correction, but they do affect the vertex-wise correction using the maximum
    // stat because they end up being way to big since there is nothing around
    // them to smooth with.
    if(MRIgetVoxVal(mask,102161,0,0,0) > 0.5 || MRIgetVoxVal(mask,102162,0,0,0) > 0.5){
      printf("Removing vertices 102161 and 102162 from mask\n");
      MRIsetVoxVal(mask,102161,0,0,0,0);
      MRIsetVoxVal(mask,102162,0,0,0,0);
    }
  }
  if(mask && (MRIgetVoxVal(mask,102161,0,0,0) > 0.5 || MRIgetVoxVal(mask,102162,0,0,0) > 0.5))
    printf("Not removing vertices 102161 and 102162 from mask\n");

  if(mask && SaveMask){
    sprintf(tmpstr,"%s/mask.mgh",OutTop);
    printf("Saving mask to %s\n",tmpstr);
    err = MRIwrite(mask,tmpstr);
    if(err) exit(1);
  }

  // Compute search space
  searchspace = 0;
  nmask = 0;
  for(n=0; n < surf->nvertices; n++){
    if(!mask) continue;
    if(mask && MRIgetVoxVal(mask,n,0,0,0) < 0.5) continue;
    searchspace += surf->vertices[n].area;
    nmask++;
  }
  printf("Found %d voxels in mask\n",nmask);

  // Make a list of vertex numbers
  maskoutvtxno = (int *) calloc(surf->nvertices-nmask,sizeof(int));
  nmaskout = 0;
  for(n=0; n < surf->nvertices; n++){
    if(mask && MRIgetVoxVal(mask,n,0,0,0) > 0.5) continue;
    maskoutvtxno[nmaskout] = n;
    nmaskout++;
  }
  printf("Found %d voxels in mask\n",nmask);

  if(surf->group_avg_surface_area > 0)
    searchspace *= (surf->group_avg_surface_area/surf->total_area);
  printf("search space %g mm2\n",searchspace);
  avgvtxarea = searchspace/nmask;
  printf("average vertex area %g mm2\n",avgvtxarea);

  // Determine how many iterations are needed for each FWHM
  nSmoothsList = (int *) calloc(sizeof(int),nFWHMList);
  for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++){
    nSmoothsList[nthFWHM] = MRISfwhm2niters(FWHMList[nthFWHM], surf);
    printf("%2d %5.1f  %4d\n",nthFWHM,FWHMList[nthFWHM],nSmoothsList[nthFWHM]);
    if(fpLog) fprintf(fpLog,"%2d %5.1f  %4d\n",nthFWHM,FWHMList[nthFWHM],nSmoothsList[nthFWHM]);
  }
  printf("\n");

  // Allocate the CSDs
  for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++){
    for(nthThresh = 0; nthThresh < nThreshList; nthThresh++){
      for(nthSign = 0; nthSign < nSignList; nthSign++){
	csd = CSDalloc();
	sprintf(csd->simtype,"%s","null-z");
	sprintf(csd->anattype,"%s","surface");
	sprintf(csd->subject,"%s",subject);
	sprintf(csd->hemi,"%s",hemi);
	sprintf(csd->contrast,"%s","NA");
	csd->seed = SynthSeed;
	csd->nreps = nRepetitions;
	csd->thresh = ThreshList[nthThresh];
	csd->threshsign = SignList[nthSign];
	csd->nullfwhm = FWHMList[nthFWHM];
	csd->varfwhm = -1;
	csd->searchspace = searchspace;
	CSDallocData(csd);
	csdList[nthFWHM][nthThresh][nthSign] = csd;
      }
    }
  }

  sprintf(tmpstr,"%s/fwhm2niters.dat",OutTop);
  fp = fopen(tmpstr,"w");
  for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++)
    fprintf(fp,"%5.1f %4d\n",FWHMList[nthFWHM],nSmoothsList[nthFWHM]);
  fclose(fp);

  // Alloc the z and sig maps
  z   = MRIallocSequence(surf->nvertices, 1,1, MRI_FLOAT, 1);
  sig = MRIallocSequence(surf->nvertices, 1,1, MRI_FLOAT, 1);

  // Set up the random field specification
  rfs = RFspecInit(SynthSeed,NULL);
  rfs->name = strcpyalloc("gaussian");
  rfs->params[0] = 0;
  rfs->params[1] = 1;

  printf("Thresholds (%d): ",nThreshList);
  for(n=0; n < nThreshList; n++) printf("%5.2f ",ThreshList[n]);
  printf("\n");
  printf("Signs (%d): ",nSignList);
  for(n=0; n < nSignList; n++)  printf("%2d ",SignList[n]);
  printf("\n");
  printf("FWHM (%d): ",nFWHMList);
  for(n=0; n < nFWHMList; n++) printf("%5.2f ",FWHMList[n]);
  printf("\n");

  // Set up pointer to the val field
  ppVal = (float **) calloc(surf->nvertices,sizeof(float *));
  ppVal0 = ppVal;
  ppSig = (float **) calloc(surf->nvertices,sizeof(float *));
  ppSig0 = ppSig;
  ppZ = (float **) calloc(surf->nvertices,sizeof(float *));
  ppZ0 = ppZ;
  for(k=0; k < surf->nvertices; k++){
    ppVal[k] = &(surf->vertices[k].val);
    ppSig[k] = &(MRIFseq_vox(sig,k,0,0,0));
    ppZ[k] = &(MRIFseq_vox(z,k,0,0,0));
  }

  // Start the simulation loop
  printf("\n\nStarting Simulation over %d Repetitions\n",nRepetitions);
  if(fpLog) fprintf(fpLog,"\n\nStarting Simulation over %d Repetitions\n",nRepetitions);
  mytimer.reset() ;
  for(nthRep = 0; nthRep < nRepetitions; nthRep++){
    msecTime = mytimer.milliseconds() ;
    printf("%5d %7.2f ",nthRep,(msecTime/1000.0)/60);
    fflush(stdout);
    if(fpLog) {
      fprintf(fpLog,"%5d %7.1f ",nthRep,(msecTime/1000.0)/60);
      fflush(fpLog);
    }
    // Synthesize an unsmoothed z map
    RFsynth(z,rfs,mask); 
    nSmoothsPrev = 0;
    
    // Loop through FWHMs
    for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++){
      printf("%d ",nthFWHM);
      fflush(stdout);
      if(fpLog) {
	fprintf(fpLog,"%d ",nthFWHM);
	fflush(fpLog);
      }
      nSmoothsDelta = nSmoothsList[nthFWHM] - nSmoothsPrev;
      nSmoothsPrev = nSmoothsList[nthFWHM];
      // Incrementally smooth z
      MRISsmoothMRIFastFrame(surf, z, 0, nSmoothsDelta, mask);
      // Rescale
      RFrescale(z,rfs,mask,z);
      // Slightly tortured way to get the right p-values because
      //   RFstat2P() computes one-sided, but I handle sidedness
      //   during thresholding.
      // First, use zabs to get a two-sided pval bet 0 and 0.5
      zabs = MRIabs(z,zabs);
      p = RFstat2P(zabs,rfs,mask,0,p);
      // Next, mult pvals by 2 to get two-sided bet 0 and 1
      MRIscalarMul(p,p,2.0);
      sig = MRIlog10(p,NULL,sig,1); // sig = -log10(p)

      for(nthSign = 0; nthSign < nSignList; nthSign++){
	csd = csdList[nthFWHM][0][nthSign]; // just need csd->threshsign

	// If test is not ABS then apply the sign
	if(csd->threshsign != 0) MRIsetSign(sig,z,0);

	// Get the max stats
	sigmax = MRIframeMax(sig,0,mask,csd->threshsign,
			     &cmax,&rmax,&smax);
	zmax = MRIgetVoxVal(z,cmax,rmax,smax,0);
	if(csd->threshsign == 0){
	  zmax = fabs(zmax);
	  sigmax = fabs(sigmax);
	}
	// Mask
	if(mask) {
	  //MRImask(sig,mask,sig,0.0,0.0); // a little slow
	  for(k=0; k < nmaskout; k++) MRIsetVoxVal(sig, maskoutvtxno[k],0,0,0, 0.0);
	}

	// Copy sig to vertexval
	//MRIScopyMRI(surf, sig, 0, "val"); // MRIScopyMRI() is a little slow
	ppVal = ppVal0;
	ppSig = ppSig0;
	for(k=0; k < surf->nvertices; k++) *(*ppVal++) = *(*ppSig++);

	for(nthThresh = 0; nthThresh < nThreshList; nthThresh++){
	  csd = csdList[nthFWHM][nthThresh][nthSign];

	  // Surface clustering
	  // Set the threshold
	  if(csd->threshsign == 0) threshadj = csd->thresh;
	  else threshadj = csd->thresh - log10(2.0); // one-sided test
	  // Compute clusters
	  SurfClustList = sclustMapSurfClusters(surf,threshadj,-1,csd->threshsign,
						0,&nClusters,NULL,NULL);
	  // Actual area of cluster with max area
	  csize  = sclustMaxClusterArea(SurfClustList, nClusters);
	  // Number of vertices of cluster with max number of vertices. 
	  // Note: this may be a different cluster from above!
	  csizen = sclustMaxClusterCount(SurfClustList, nClusters);
	  cweightvtx = sclustMaxClusterWeightVtx(SurfClustList, nClusters, csd->threshsign);
	  // Area of this cluster based on average vertex area. This just scales
	  // the number of vertices.
	  csizeavg = csizen * avgvtxarea;
	  if(UseAvgVtxArea) csize = csizeavg;
	  // Store results
	  csd->nClusters[nthRep] = nClusters;
	  csd->MaxClusterSize[nthRep] = csize;
	  csd->MaxClusterSizeVtx[nthRep] = csizen;
	  csd->MaxClusterWeightVtx[nthRep] = cweightvtx;
	  csd->MaxSig[nthRep] = sigmax;
	  csd->MaxStat[nthRep] = zmax;
	} // Sign
      } // Thresh
    } // FWHM
    printf("\n");
    if(fpLog) fprintf(fpLog,"\n");
    if(SaveEachIter || fio_FileExistsReadable(SaveFile)) SaveOutput();
    if(fio_FileExistsReadable(StopFile)) {
      printf("Found stop file %s\n",StopFile);
      goto finish;
    }
  } // Simulation Repetition

 finish:

  SaveOutput();

  msecTime = mytimer.milliseconds() ;
  printf("Total Sim Time %g min (%g per rep)\n",
	 msecTime/(1000*60.0),(msecTime/(1000*60.0))/nthRep);
  if(fpLog) fprintf(fpLog,"Total Sim Time %g min (%g per rep)\n",
		    msecTime/(1000*60.0),(msecTime/(1000*60.0))/nthRep);

  if(DoneFile){
    fp = fopen(DoneFile,"w");
    fprintf(fp,"%g\n",msecTime/(1000*60.0));
    fclose(fp);
  }
  printf("mri_mcsim done\n");
  if(fpLog){
    fprintf(fpLog,"mri_mcsim done\n");
    fclose(fpLog);
  }
  exit(0);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused, nth;
  char **pargv, *option ;
  double fwhm;

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
    else if (!strcasecmp(option, "--save-weight")) SaveWeight = 1;
    else if (!strcasecmp(option, "--no-fix-fsalh")) FixFSALH = 0;

    else if (!strcasecmp(option, "--test")){
      double fwhm;
      nFWHMList = 0;
      for(fwhm = 3; fwhm <= 6; fwhm++){
	FWHMList[nFWHMList] = fwhm;
	nFWHMList++;
      }
      nThreshList = 2;
      ThreshList[0] = 2.0;
      ThreshList[1] = 3.0; 
      nSignList = 2;
      subject = "fsaverage5";
      hemi = "lh";
      nRepetitions = 2;
      csdbase = "junk";
      SynthSeed = 53;
    }
    else if (!strcasecmp(option, "--o")) {
      if(nargc < 1) CMDargNErr(option,1);
      OutTop = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--done")) {
      if(nargc < 1) CMDargNErr(option,1);
      DoneFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--log")) {
      if(nargc < 1) CMDargNErr(option,1);
      LogFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--stop")) {
      if(nargc < 1) CMDargNErr(option,1);
      StopFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--save")) {
      if(nargc < 1) CMDargNErr(option,1);
      SaveFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--save-iter")) {
      SaveEachIter = 1;
    } 
    else if (!strcasecmp(option, "--no-save-iter")) {
      SaveEachIter = 0;
    } 
    else if (!strcasecmp(option, "--base")) {
      if(nargc < 1) CMDargNErr(option,1);
      csdbase = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--surf") || 
	     !strcasecmp(option, "--surface")) {
      if(nargc < 2) CMDargNErr(option,1);
      subject = pargv[0];
      hemi = pargv[1];
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--surfname")) {
      if(nargc < 1) CMDargNErr(option,1);
      surfname = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--avgvtxarea"))    UseAvgVtxArea = 1;
    else if (!strcasecmp(option, "--no-avgvtxarea")) UseAvgVtxArea = 0;
    else if (!strcasecmp(option, "--seed")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&SynthSeed);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mask")) {
      if(nargc < 1) CMDargNErr(option,1);
      MaskFile = pargv[0];
      LabelFile = NULL;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--label")) {
      if(nargc < 1) CMDargNErr(option,1);
      LabelFile = pargv[0];
      MaskFile = NULL;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--no-label")) {
      LabelFile = NULL;
    }
    else if (!strcmp(option, "--sd")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--no-save-mask")) SaveMask = 0;
    else if (!strcasecmp(option, "--nreps")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nRepetitions);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--fwhm")) {
      if (nargc < 1) CMDargNErr(option,1);
      nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) ){
	sscanf(pargv[nth],"%lf",&fwhm);
	FWHMList[nFWHMList] = fwhm;
	nFWHMList++;
	nth++;
      }
      nargsused = nth;
    } 
    else if(!strcasecmp(option, "--thresh")) {
      double thresh;
      if (nargc < 1) CMDargNErr(option,1);
      nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) ){
	sscanf(pargv[nth],"%lf",&thresh);
	ThreshList[nThreshList] = thresh;
	nThreshList++;
	nth++;
      }
      nargsused = nth;
    } 
    else if(!strcasecmp(option, "--nsign")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nSignList);
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
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("%s \n",Progname) ;
  printf("\n");
  printf("   --o top-output-dir\n");
  printf("   --base csdbase\n");
  printf("   --surface subjectname hemi\n");
  printf("   --nreps nrepetitions\n");
  printf("   --fwhm FWHM <FWHM2 FWHM3 ...>\n");
  printf("   --fwhm-max FWHMMax : sim with fwhm=1:FWHMMax (default %g)\n",fwhmmax);
  printf("   \n");
  printf("   --avgvtxarea : report cluster area based on average vtx area\n");
  printf("   --seed randomseed : default is to choose based on ToD\n");
  printf("   --label labelfile : default is ?h.cortex.label \n");
  printf("   --mask maskfile : instead of label\n");
  printf("   --no-label : do not use a label to mask\n");
  printf("   --no-save-mask : do not save mask to output (good for mult jobs)\n");
  printf("   --surfname surfname : default is white\n");
  printf("\n");
  printf("   --log  LogFile \n");
  printf("   --done DoneFile : will create DoneFile when finished\n");
  printf("   --stop stopfile : default is ourdir/mri_mcsim.stop \n");
  printf("   --save savefile : default is ourdir/mri_mcsim.save \n");
  printf("   --save-iter : save output after each iteration \n");
  printf("   --sd SUBJECTS_DIR\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
printf("\n");
printf("This program computes tables for performing multiple comparisons on\n");
printf("the surface based on a monte carlo simulation using white gaussian\n");
printf("noise smoothed on the surface. mri_mcsim was used to generate the\n");
printf("tables found in $FREESURFER_HOME/average/mult-comp-cor and used by\n");
printf("mri_glmfit-sim with the --cache option. The tables in mult-comp-cor\n");
printf("are for fsaverage and fsaverage_sym for a search space over the whole\n");
printf("cortex. mri_mcsim can be used to create new tables for a new subject\n");
printf("or for fsaverage/fsaverage_sum with a smaller search space.\n");
printf("\n");
printf("Example 1: Create tables for a new subject for whole hemisphere\n");
printf("\n");
printf("mri_mcsim --o /path/to/mult-comp-cor/newsubject/lh/cortex --base mc-z \n");
printf("  --save-iter  --surf newsubject lh --nreps 10000\n");
printf("\n");
printf("This may take hours (or even days) to run; see below for\n");
printf("parallelizing.  When running mri_glmfit-sim, add --cache-dir\n");
printf("/path/to/mult-comp-cor\n");
printf("\n");
printf("Example 2: Create tables for superior temporal gyrus for fsaverage\n");
printf("\n");
printf("First, create the label by running\n");
printf("\n");
printf("mri_annotation2label --subject fsaverage --hemi lh --outdir labeldir\n");
printf("\n");
printf("mri_mcsim --o /path/to/mult-comp-cor/fsaverage/lh/superiortemporal --base mc-z \n");
printf("  --save-iter  --surf fsaverage lh --nreps 10000 \n");
printf("  --label labeldir/lh.superiortemporal.label\n");
printf("\n");
printf("When running mri_glmfit, make sure to use   --label labeldir/lh.superiortemporal.label\n");
printf("When running mri_glmfit-sim, add --cache-dir /path/to/mult-comp-cor --cache-label superiortemporal\n");
printf("\n");
printf("Example 3: running simulations in parallel (two jobs, 5000 iterations\n");
printf("each for a total of 10000)\n");
printf("\n");
printf("mri_mcsim --o /path/to/mult-comp-cor/fsaverage/lh/superiortemporal --base mc-z.j001 \n");
printf("  --save-iter  --surf fsaverage lh --nreps 5000\n");
printf("  --label labeldir/lh.superiortemporal.label\n");
printf("\n");
printf("mri_mcsim --o /path/to/mult-comp-cor/fsaverage/lh/superiortemporal --base mc-z.j002 \n");
printf("  --save-iter  --surf fsaverage lh --nreps 500\n");
printf("  --label labeldir/lh.superiortemporal.label\n");
printf("\n");
printf("When those jobs are done, merge the results into a single table with\n");
printf("\n");
printf("mri_surfcluster \n");
printf("  --csd /path/to/mult-comp-cor/fsaverage/lh/superiortemporal/fwhm10/abs/th20/mc-z.j001.csd \n");
printf("  --csd /path/to/mult-comp-cor/fsaverage/lh/superiortemporal/fwhm10/abs/th20/mc-z.j002.csd \n");
printf("  --csd-out /path/to/mult-comp-cor/fsaverage/lh/superiortemporal/fwhm10/abs/th20/mc-z.csd \n");
printf("  --csdpdf  /path/to/mult-comp-cor/fsaverage/lh/superiortemporal/fwhm10/abs/th20/mc-z.cdf --csdpdf-only\n");
printf("\n");
printf("Repeat the above command for each FWHM, sign (pos, neg, abs) and threshold\n");
printf("\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if(subject == NULL) {
    printf("ERROR: must specify a surface\n");
    exit(1);
  }
  if(OutTop == NULL){
    printf("ERROR: need to spec an output directory\n");
    exit(1);
  }
  if(csdbase==NULL) {
    printf("ERROR: need to specify a csd base\n");
    exit(1);
  }
  if(MaskFile && LabelFile) {
    printf("ERROR: cannot specify both a mask and a label\n");
    exit(1);
  }
  if(nRepetitions < 1) {
    printf("ERROR: need to specify number of simulation repitions\n");
    exit(1);
  }
  if(nFWHMList == 0){
    double fwhm;
    nFWHMList = 0;
    for(fwhm = 1; fwhm <= fwhmmax; fwhm++){
      FWHMList[nFWHMList] = fwhm;
      nFWHMList++;
    }
  }
  if(nThreshList == 0){
    nThreshList = 6;
    ThreshList[0] = 1.3; 
    ThreshList[1] = 2.0;
    ThreshList[2] = 2.3; 
    ThreshList[3] = 3.0; 
    ThreshList[4] = 3.3; 
    ThreshList[5] = 4.0;
  }
  if(StopFile == NULL) {
    sprintf(tmpstr,"%s/mri_mcsim.stop",OutTop);
    StopFile = strcpyalloc(tmpstr);
  }
  if(SaveFile == NULL) {
    sprintf(tmpstr,"%s/mri_mcsim.save",OutTop);
    SaveFile = strcpyalloc(tmpstr);
  }

  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"OutTop  %s\n",OutTop);
  fprintf(fp,"CSDBase  %s\n",csdbase);
  fprintf(fp,"nreps    %d\n",nRepetitions);
  fprintf(fp,"fwhmmax  %g\n",fwhmmax);
  fprintf(fp,"subject  %s\n",subject);
  fprintf(fp,"hemi     %s\n",hemi);
  fprintf(fp,"surfname %s\n",surfname);
  fprintf(fp,"FixVertexAreaFlag %d\n",MRISgetFixVertexAreaValue());
  if(MaskFile) fprintf(fp,"mask     %s\n",MaskFile);
  fprintf(fp,"UseAvgVtxArea %d\n",UseAvgVtxArea);
  fprintf(fp,"SaveFile %s\n",SaveFile);
  fprintf(fp,"StopFile %s\n",StopFile);
  fprintf(fp,"UFSS %s\n",getenv("USE_FAST_SURF_SMOOTHER"));
  fflush(fp);
  return;
}

int SaveOutput(void)
{
  int nthSign, nthFWHM, nthThresh;
  FILE *fp;

  // Save output
  // printf("Saving results nreps = %d\n",nthRep);
  for(nthFWHM=0; nthFWHM < nFWHMList; nthFWHM++){
    for(nthThresh = 0; nthThresh < nThreshList; nthThresh++){
      for(nthSign = 0; nthSign < nSignList; nthSign++){
	csd = csdList[nthFWHM][nthThresh][nthSign];
	if(csd->threshsign ==  0) signstr = "abs"; 
	if(csd->threshsign == +1) signstr = "pos"; 
	if(csd->threshsign == -1) signstr = "neg"; 
	csd->nreps = nthRep;
	sprintf(tmpstr,"%s/fwhm%02d/%s/th%02d/%s.csd",
		OutTop,(int)round(FWHMList[nthFWHM]),
		signstr,(int)round(10*ThreshList[nthThresh]),
		csdbase);
	fp = fopen(tmpstr,"w");
	fprintf(fp,"# ClusterSimulationData 2\n");
	fprintf(fp,"# mri_mcsim\n");
	fprintf(fp,"# %s\n",cmdline);
	fprintf(fp,"# %s\n",getVersion().c_str());
	fprintf(fp,"# hostname %s\n",uts.nodename);
	fprintf(fp,"# machine  %s\n",uts.machine);
	fprintf(fp,"# runtime_min %g\n",msecTime/(1000*60.0));
	fprintf(fp,"# nvertices-total %d\n",surf->nvertices);
	fprintf(fp,"# nvertices-search %d\n",nmask);
	if(mask) fprintf(fp,"# masking 1\n");
	else     fprintf(fp,"# masking 0\n");
	//fprintf(fp,"# FWHM %g\n",FWHMList[nthFWHM]);
	fprintf(fp,"# SmoothLevel %d\n",nSmoothsList[nthFWHM]);
	fprintf(fp,"# UseAvgVtxArea %d\n",UseAvgVtxArea);
	if(!SaveWeight) CSDprint(fp, csd);
	else CSDprintWeight(fp, csd);
	fclose(fp);
      }
    }
  }
  return(0);
}

