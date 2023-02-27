/**
 * @brief binarizes an image
 *
 * Program to binarize a volume (or volume-encoded surface file). Can also
 * be used to merge with other binarizations. Binarization can be done
 * based on threshold or on matched values.
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

Program to binarize a volume (or volume-encoded surface file). Can also
be used to merge with other binarizations. Binarization can be done
based on threshold or on matched values.

--i invol

Input volume to be binarized.

--min min
--max max
--rmin rmin
--rmax rmax

Minimum and maximum thresholds. If the value at a voxel is >= min and
<= max, then its value in the output will be 1 (or --binval),
otherwise it is 0 (or --binvalnot) or the value of the merge volume at
that voxel.  By default, min = -infinity and max = +infinity, but you
must set one of the thresholds. Cannot be used with --match. If --rmin
or --rmax are specified, then min (or max) is computed as rmin (or
rmax) times the global mean of the input.

--pct P

Set min threshold so that the top P percent of the voxels are captured
in the output mask. The percent will be computed based on the number of
voxels in the volume (if not input mask is specified) or within the
input mask.

--fdr fdrthreshold

Set min threshold to achieve a given FDR. By default, it uses the
absolute value but this can be changed with --fdr-pos and
--fdr-neg. If a mask is passed, it will compute the voxel-wise
threshold only with in the places where mask > 0.5.  The mask
threshold will be ignored.

--match matchvalue <matchvalue2 ...>

Binarize based on matching values. Any number of match values can be 
specified. Cannot be used with --min/--max.

--o outvol

Path to output volume.

--count countfile

Save number of voxels that meet match criteria in ascii
countefile. Four numbers are saved: the number of voxels that match
(nhits), the volume of the voxels that match, the total number of
voxels in the volume (nvoxtot), and the percent matching
(100*nhits/nvoxtot).

--binval    binval
--binvalnot binvalnot

Value to use for those voxels that are in the threshold/match
(--binval) or out of the range (--binvalnot). These must be integer
values. binvalnot only applies when a merge volume is NOT specified.

--replace V1 V2

Replace every occurrence of (int) value V1 with value V2. Multiple 
--replace args are possible.

--replaceonly V1 V2

Replace every occurrence of (int) value V1 with value V2. Multiple 
--replaceonly args are possible. Other locations in the source volume
will be propagated to the output (unlike --replace which masks those
locations).

--replace-nn V1 W

Replace every occurrence of (int) value V1 with that of its nearest
neighbor voxel within a window of W voxels. Multiple --replace-nn args
are possible.

--replaceonly-nn V1 W

Replace every occurrence of (int) value V1 with that of its nearest
neighbor voxel within a window of W voxels. Multiple --replaceonly-nn
args are possible. Other locations in the source volume will be propagated 
to the output (unlike --replace-nn which masks those locations).

--frame frameno

Use give frame of the input. 0-based. Default is 0.

--frame-sum

Sum the frames together before applying threshold.

--frame-and

Treat the multi-frame volume as binary 'AND' the frames together. This
takes an intersection of the individual frames. You do not need to 
specify a --min (the min will be set to nframes-0.5).

--copy copyvol

copy values from copyvol into the output, except where replaceval is matched

--merge mergevol

Merge binarization with the mergevol. If the voxel is within the threshold
range (or matches), then its value will be binval. If not, then it will 
inherit its value from the value at that voxel in mergevol. mergevol must 
be the same dimension as the input volume. Combining this with --binval 
allows you to construct crude segmentations.

--mask maskvol
--mask-thresh thresh

Mask input with mask. The mask volume is itself binarized at thresh
(default is 0.5). If a voxel is not in the mask, then it will be assigned
binvalnot or the value from the merge volume.

--zero-edges

Set the first and last planes in all dimensions to 0 (or --binvalnot). This
makes sure that all the voxels on the edge of the imaging volume are 0.

--zero-slice-edges

Same as --zero-edges, but only for slices.

--uchar

Save output file in 'uchar' format.

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <float.h>

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
#include "cma.h"
#include "colortab.h"
#include "region.h"
#include "romp_support.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *InVolFile=NULL;
char *OutVolFile=NULL;
char *MergeVolFile=NULL;
char *CopyVolFile=NULL;
char *MaskVolFile=NULL;
double MinThresh, MaxThresh;
int MinThreshSet=0, MaxThreshSet=0;
double RMinThresh, RMaxThresh;
int RMinThreshSet=0, RMaxThreshSet=0;
char *CountFile = NULL;

int BinVal=1;
int BinValNot=0;
int frame=0;
int DoFrameLoop=1;
int DoAbs=0;
int DoNeg=0;
int ZeroColEdges = 0;
int ZeroRowEdges = 0;
int ZeroSliceEdges = 0;

int DoMatch = 0;
int nMatch = 0;
int MatchValues[1000];
int Matched = 0;

MRI *InVol,*OutVol,*MergeVol,*MaskVol=NULL, *CopyVol;
double MaskThresh = 0.5;

int nErode2d = 0;
int nErode3d = 0;
int nDilate3d = 0;
int DoBinCol = 0;

int mriTypeUchar = 0;
int DoFrameSum = 0;
int DoFrameAnd = 0;
int DoPercent = 0;
double TopPercent = -1;
double FDR;
int DoFDR = 0;
int FDRSign = 0;

int nErodeNN=0, NNType=0;

static int replace_only = 0 ;
int nReplace = 0, SrcReplace[1000], TrgReplace[1000], 
    nReplaceNN = 0, SrcReplaceNN[1000], ReplaceWindowNN[1000];
char *SurfFile=NULL;
int nsmoothsurf=0;

int noverbose = 0;
int DoBB = 0, nPadBB=0;
int DoCount = 1;
int ReverseFaceOrder = 0;
int FillHoles = 0;
int RemoveIslands = 0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, c, nhits, n, mriType,nvox;
  int fstart, fend, nframes;
  double gmean,gstd,gmax,voxvol;
  FILE *fp;
  MRI *mritmp;

  nargs = handleVersionOption(argc, argv, "mri_binarize");
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
  if (noverbose == 0)
    dump_options(stdout);

  // Load the input volume
  InVol = MRIread(InVolFile);
  if (InVol==NULL) exit(1);
  if (frame >= InVol->nframes) {
    printf("ERROR: requested frame=%d >= nframes=%d\n",
           frame,InVol->nframes);
    exit(1);
  }

  // Load the mask volume (if needed)
  if (MaskVolFile) {
    printf("Loading mask %s\n",MaskVolFile);
    MaskVol = MRIread(MaskVolFile);
    if (MaskVol==NULL) exit(1);
    if (MaskVol->width != InVol->width) {
      printf("ERROR: dimension mismatch between input and mask volumes\n");
      exit(1);
    }
    if (MaskVol->height != InVol->height) {
      printf("ERROR: dimension mismatch between input and mask volumes\n");
      exit(1);
    }
    if (MaskVol->depth != InVol->depth) {
      printf("ERROR: dimension mismatch between input and mask volumes\n");
      exit(1);
    }
  }

  if(DoFDR){
    double FDRThresh;
    // 1 = assume -log10(p)
    MRIfdr2vwth(&InVol, 1, &frame, FDR, FDRSign, 1, &MaskVol, &FDRThresh, NULL);
    printf("FDR %g, Sign=%d, Thresh = %g\n",FDR,FDRSign,FDRThresh);
    if(FDRSign == 0) {
      DoAbs = 1;
      MinThresh = FDRThresh;
      MinThreshSet=1;
    }
    if(FDRSign == +1) {
      MinThresh = FDRThresh;
      MinThreshSet=1;
    }
    if(FDRSign == -1) {
      MaxThresh = -FDRThresh;
      MaxThreshSet=1;
    }
  }

  if(DoFrameSum || DoFrameAnd) {
    printf("Summing over %d frames\n",InVol->nframes);
    InVol = MRIframeSum(InVol,InVol);
    frame = 0;
  }
  if(DoFrameAnd) MinThresh = InVol->nframes - 0.5;

  if(DoAbs) {
    printf("Removing sign from input\n");
    MRIabs(InVol,InVol);
  }
  if(DoNeg) {
    printf("Negating input\n");
    MRImultiplyConst(InVol,-1.0,InVol);
  }

  if(RMinThreshSet || RMaxThreshSet){
    printf("Computing global statistics \n");
    RFglobalStats(InVol, NULL, &gmean, &gstd, &gmax);
    printf("mean = %g, std = %g, max = %g\n",gmean, gstd, gmax);
    if(RMinThreshSet)  {
      MinThresh = gmean * RMinThresh;
      MinThreshSet = 1;
      printf("Setting min thresh to %g\n",MinThresh);
    }
    if(RMaxThreshSet) {
      MaxThresh = gmean * RMaxThresh;
      MaxThreshSet = 1;
      printf("Setting max thresh to %g\n",MinThresh);
    }
  }

  // Load the merge volume (if needed)
  if (MergeVolFile) {
    MergeVol = MRIread(MergeVolFile);
    if (MergeVol==NULL) exit(1);
    if (MergeVol->width != InVol->width) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
    if (MergeVol->height != InVol->height) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
    if (MergeVol->depth != InVol->depth) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
  }

  if (CopyVolFile) {
    CopyVol = MRIread(CopyVolFile);
    if (CopyVol==NULL) exit(1);
    if (CopyVol->width != InVol->width) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
    if (CopyVol->height != InVol->height) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
    if (CopyVol->depth != InVol->depth) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
  }

  if(DoPercent) {
    printf("Computing threshold based on top %g percent\n",TopPercent);
    MinThresh = MRIpercentThresh(InVol, MaskVol, frame, TopPercent);
    printf("  Threshold set to %g\n",MinThresh);
  }

  if(DoFrameLoop){
    fstart = 0;
    fend = InVol->nframes-1;
  }
  else {
    fstart = frame;
    fend = frame;
  }
  nframes = fend - fstart + 1;
  if (noverbose == 0)
    printf("fstart = %d, fend = %d, nframes = %d\n",fstart,fend,nframes);

  // Prepare the output volume
  mriType = MRI_INT;
  if (mriTypeUchar) mriType = MRI_UCHAR;
  OutVol = MRIallocSequence(InVol->width,InVol->height,InVol->depth,mriType,nframes);
  if (OutVol == NULL) exit(1);
  MRIcopyHeader(InVol, OutVol);

  nhits = 0;
  if(!replace_only){
    // Binarize
    for(frame = fstart; frame <= fend; frame++){
      // Do openmp for columns because often nframes = 1
      #ifdef HAVE_OPENMP
      printf("Starting parallel %d\n",omp_get_num_threads());
      #pragma omp parallel for if_ROMP(assume_reproducible) reduction(+ : nhits)
      #endif
      for (c=0; c < InVol->width; c++) {
        int r,s;
        double mergeval = BinValNot, maskval, val;
	for (r=0; r < InVol->height; r++) {
	  for (s=0; s < InVol->depth; s++) {
	    if(MergeVol) mergeval = MRIgetVoxVal(MergeVol,c,r,s,frame);
	    
	    // Skip if on the edge
	    if( (ZeroColEdges &&   (c == 0 || c == InVol->width-1))  ||
		(ZeroRowEdges &&   (r == 0 || r == InVol->height-1)) ||
		(ZeroSliceEdges && (s == 0 || s == InVol->depth-1)) ){
	      MRIsetVoxVal(OutVol,c,r,s,frame-fstart,mergeval);
	      continue;
	    }
	    
	    // Skip if not in the mask
	    if(MaskVol) {
	      maskval = MRIgetVoxVal(MaskVol,c,r,s,0);
	      if(maskval < MaskThresh){
		MRIsetVoxVal(OutVol,c,r,s,frame-fstart,mergeval);
		continue;
	      }
	    }
	    
	    // Get the value at this voxel
	    val = MRIgetVoxVal(InVol,c,r,s,frame);
	    
	    if(DoMatch){
	      // Check for a match
              int Matched = 0;
	      for(n=0; n < nMatch; n++){
		if(fabs(val - MatchValues[n]) < 2*FLT_MIN){
		  MRIsetVoxVal(OutVol,c,r,s,frame-fstart,BinVal);
		  Matched = 1;
		  nhits ++;
		  break;
		}
	      }
	      if(!Matched) MRIsetVoxVal(OutVol,c,r,s,frame-fstart,mergeval);
	    }
	    else{
	      // Determine whether it is in range
	      if((MinThreshSet && (val < MinThresh)) ||
		 (MaxThreshSet && (val > MaxThresh))){
		// It is NOT in the Range
		MRIsetVoxVal(OutVol,c,r,s,frame-fstart,mergeval);
	      }
	      else {
		// It is in the Range
		MRIsetVoxVal(OutVol,c,r,s,frame-fstart,BinVal);
		nhits ++;
	      }
	    }
	    
	  } // slice
	} // row
      } // col
    } // frame
  } // if(!replace_only)

  if(nReplace != 0) {
    if (replace_only)
    {
      printf("Replacing %d and propagating source list\n",nReplace);
      OutVol = MRIcopy(InVol, NULL) ;
      for(n=0; n < nReplace; n++) 
      {
	printf("%2d:  %4d %4d\n",n+1,SrcReplace[n],TrgReplace[n]);
	MRIreplaceValues(OutVol, OutVol, SrcReplace[n], TrgReplace[n]);
      }
    }
    else
    {
      printf("Replacing %d\n",nReplace);
      for(n=0; n < nReplace; n++) printf("%2d:  %4d %4d\n",n+1,SrcReplace[n],TrgReplace[n]);
      OutVol = MRIreplaceList(InVol, SrcReplace, TrgReplace, nReplace, MaskVol, NULL);
    }
  }

  if(nReplaceNN != 0) {
    if (replace_only)
    {
      printf("Replacing %d with nearest neighbors and propagating source list\n", nReplaceNN);
      OutVol = MRIcopy(InVol, NULL) ;
    }
    else
      printf("Replacing %d with nearest neighbors\n", nReplaceNN);

    for(frame = fstart; frame <= fend; frame++){
      for (c=0; c < InVol->width; c++) {
        int r,s;
        for (r=0; r < InVol->height; r++) {
          for (s=0; s < InVol->depth; s++) {
            if (!MaskVol || MRIgetVoxVal(MaskVol, c, r, s, 0) > MaskThresh) {
              int segid = MRIgetVoxVal(InVol, c, r, s, frame);

              for (n = 0; n < nReplaceNN; n++) {
                if (segid == SrcReplaceNN[n]) {
                  int whalf = (ReplaceWindowNN[n] - 1) / 2,
                      nn_val = segid;
                  float min_dist = 100000;

                  for (int ds = -whalf; ds <= whalf; ds++) {
                    int s1 = s + ds;

                    if (s1 < 0 || s1 >= InVol->depth)	continue;

                    for (int dr = -whalf; dr <= whalf; dr++) {
                      int r1 = r + dr;

                      if (r1 < 0 || r1 >= InVol->height)	continue;

                      for (int dc = -whalf; dc <= whalf; dc++) {
                        int c1 = c + dc,
                            newsegid, isrep = 0;
                        float dist;

                        if (c1 < 0 || c1 >= InVol->width)	continue;

                        newsegid = MRIgetVoxVal(InVol, c1, r1, s1, frame);

                        for (int n1 = 0; n1 < nReplaceNN; n1++)
                          if (newsegid == SrcReplaceNN[n1]) {
                            isrep = 1;
                            break;
                          }

                        if (isrep)	continue;

                        dist = sqrt((float) dc * dc + dr * dr + ds * ds);

                        if (dist < min_dist) {
                          min_dist = dist;
                          nn_val = newsegid;
                        }
                      }		// ck
                    }		// rk
                  }		// sk

                  MRIsetVoxVal(OutVol, c, r, s, frame-fstart, nn_val);
                  break;
                }
              }	// n
            }	// mask
          }	// s
        }	// r
      }		// c
    }		// frame
  }

  if (noverbose == 0 && replace_only == 0) 
    printf("Found %d values in range\n",nhits);

  if(nDilate3d > 0){
    printf("Dilating %d voxels in 3d\n",nDilate3d);
    for(n=0; n<nDilate3d; n++) MRIdilate(OutVol,OutVol);
  }
  if(nErode3d > 0){
    printf("Eroding %d voxels in 3d\n",nErode3d);
    for(n=0; n<nErode3d; n++) MRIerode(OutVol,OutVol);
  }
  if(nErode2d > 0){
    printf("Eroding %d voxels in 2d\n",nErode2d);
    for(n=0; n<nErode2d; n++) MRIerode2D(OutVol,OutVol);
  }
  if(nErodeNN > 0){
    printf("Eroding %d voxels using %d\n",nErodeNN,NNType);
    mritmp = NULL;
    for(n=0; n<nErodeNN; n++) {
      mritmp = MRIerodeNN(OutVol,mritmp,NNType);
      MRIcopy(mritmp,OutVol);
    }
    MRIfree(&mritmp);
  }

  nhits = -1;
  if(DoCount && !replace_only){
    if(noverbose == 0) printf("Counting number of voxels in first frame\n");
    #ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(assume_reproducible) reduction(+ : nhits)
    #endif
    for (c=0; c < OutVol->width; c++) {
      double val;
      int r,s;
      for (r=0; r < OutVol->height; r++) {
	for (s=0; s < OutVol->depth; s++) {
	  // Get the value at this voxel
	  val = MRIgetVoxVal(OutVol,c,r,s,0);
	  if(fabs(val-BinVal) < .00001) nhits ++;
	} // slice
      } // row
    } // col
    if(noverbose == 0)  printf("Found %d voxels in final mask\n",nhits);
  }

  if(DoBinCol){
    printf("Filling mask with column number\n");
    MRIbinMaskToCol(OutVol, OutVol);
  }

  if(DoBB){
    // reduce volume to a smaller volume by reducing to a bounding box 
    // around non-zero voxels
    printf("Computing bounding box, npad = %d\n",nPadBB);
    MRI_REGION *region;
    region = REGIONgetBoundingBox(OutVol,nPadBB);
    REGIONprint(stdout, region);
    mritmp = MRIextractRegion(OutVol, NULL, region);
    if(mritmp == NULL) exit(1);
    MRIfree(&OutVol);
    OutVol = mritmp;
  }

  // if we didn't binarize, copy any embedded color table from the input
  if (replace_only && (InVol->ct)) OutVol->ct = CTABdeepCopy(InVol->ct);

  if(RemoveIslands){
    printf("Removing Volume Islands\n");fflush(stdout);
    MRI *tmpvol = MRIremoveVolumeIslands(OutVol, BinVal-0.5, BinVal, NULL);
    MRIfree(&OutVol);
    OutVol = tmpvol;
  }
  if(FillHoles){
    printf("Removing Volume Holes\n");fflush(stdout);
    MRI *tmpvol = MRIremoveVolumeHoles(OutVol, BinVal-0.5, BinVal, BinVal, NULL);
    MRIfree(&OutVol);
    OutVol = tmpvol;
  }

  // Save output
  if(OutVolFile) {
    printf("Writing output to %s\n",OutVolFile);
    int err = MRIwrite(OutVol,OutVolFile);
    if(err) exit(1);
  }

  if(SurfFile){
    printf("Creating surface %s\n",SurfFile);
    MRIS *surf;
    surf = MRIStessellate(OutVol,BinVal,0);
    surf->hemisphere = NO_HEMISPHERE;
    if(nsmoothsurf > 0) MRISaverageVertexPositions(surf, nsmoothsurf) ;
    if(ReverseFaceOrder){
      printf("Reversing face order\n");
      MRISreverseFaceOrder(surf);
      MRIScomputeMetricProperties(surf);
    }
    MRISwrite(surf,SurfFile);
    MRISfree(&surf);
  }

  nvox = OutVol->width * OutVol->height * OutVol->depth;
  voxvol = OutVol->xsize * OutVol->ysize * OutVol->zsize;
  if (noverbose == 0 && replace_only == 0)
    printf("Count: %d %lf %d %lf\n",nhits,nhits*voxvol,nvox,(double)100*nhits/nvox);

  if(CountFile){
    fp = fopen(CountFile,"w");
    if(fp == NULL){
      printf("ERROR: could not open %s\n",CountFile);
      exit(1);
    }
    fprintf(fp,"%d %lf %d %lf\n",nhits,nhits*voxvol,nvox,(double)100*nhits/nvox);
    fclose(fp);
  }

  if (noverbose == 0)
    printf("mri_binarize done\n");

  exit(0);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused, nth;
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
    else if (!strcasecmp(option, "--noverbose"))   noverbose = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--abs")) DoAbs = 1;
    else if (!strcasecmp(option, "--neg")) DoNeg = 1;
    else if (!strcasecmp(option, "--bincol")) DoBinCol = 1;
    else if (!strcasecmp(option, "--uchar")) mriTypeUchar = 1;
    else if (!strcasecmp(option, "--no-count")) DoCount = 0;
    else if (!strcasecmp(option, "--reverse")) ReverseFaceOrder = 1;
    else if (!strcasecmp(option, "--zero-edges")){
      ZeroColEdges = 1;
      ZeroRowEdges = 1;
      ZeroSliceEdges = 1;
    }
    else if (!strcasecmp(option, "--zero-slice-edges")) ZeroSliceEdges = 1;
    else if (!strcasecmp(option, "--fill-holes")) FillHoles = 1;
    else if (!strcasecmp(option, "--no-fill-holes")) FillHoles = 0;
    else if (!strcasecmp(option, "--remove-islands")) RemoveIslands = 1;
    else if (!strcasecmp(option, "--no-remove-islands")) RemoveIslands = 0;
    else if (!strcasecmp(option, "--ctx-wm") || !strcasecmp(option, "--wm")){
      MatchValues[nMatch++] =  2;
      MatchValues[nMatch++] = 41;
      MatchValues[nMatch++] = 77;
      MatchValues[nMatch++] = 251;
      MatchValues[nMatch++] = 252;
      MatchValues[nMatch++] = 253;
      MatchValues[nMatch++] = 254;
      MatchValues[nMatch++] = 255;
      DoMatch = 1;
    }
    else if (!strcasecmp(option, "--all-wm")){
      MatchValues[nMatch++] =  2;
      MatchValues[nMatch++] = 41;
      MatchValues[nMatch++] = 77;
      MatchValues[nMatch++] = 251;
      MatchValues[nMatch++] = 252;
      MatchValues[nMatch++] = 253;
      MatchValues[nMatch++] = 254;
      MatchValues[nMatch++] = 255;
      MatchValues[nMatch++] =  7; // Cerebellar WM
      MatchValues[nMatch++] = 46; // Cerebellar WM
      DoMatch = 1;
    }
    else if (!strcasecmp(option, "--ventricles")){
      MatchValues[nMatch++] =  4; // Left-Lateral-Ventricle
      MatchValues[nMatch++] =  5; // Left-Inf-Lat-Vent
      MatchValues[nMatch++] = 14; // 3rd-Ventricle
      MatchValues[nMatch++] = 43; // Right-Lateral-Ventricle
      MatchValues[nMatch++] = 44; // Right-Inf-Lat-Vent
      MatchValues[nMatch++] = 72; // 5th-Ventricle
      MatchValues[nMatch++] = 31; // Left-choroid-plexus 
      MatchValues[nMatch++] = 63; // Right-choroid-plexus 
      //MatchValues[nMatch++] = 15; // 4th-Ventricle. Good to include?
      DoMatch = 1;
    }

    else if (!strcasecmp(option, "--wm+vcsf")){
      MatchValues[nMatch++] =  2;
      MatchValues[nMatch++] = 41;
      MatchValues[nMatch++] = 77;
      MatchValues[nMatch++] = 251;
      MatchValues[nMatch++] = 252;
      MatchValues[nMatch++] = 253;
      MatchValues[nMatch++] = 254;
      MatchValues[nMatch++] = 255;
      MatchValues[nMatch++] =  7; // Cerebellar WM
      MatchValues[nMatch++] = 46; // Cerebellar WM
      MatchValues[nMatch++] =  4; // Left-Lateral-Ventricle
      MatchValues[nMatch++] =  5; // Left-Inf-Lat-Vent
      MatchValues[nMatch++] = 14; // 3rd-Ventricle
      MatchValues[nMatch++] = 43; // Right-Lateral-Ventricle
      MatchValues[nMatch++] = 44; // Right-Inf-Lat-Vent
      MatchValues[nMatch++] = 72; // 5th-Ventricle
      MatchValues[nMatch++] = 31; // Left-choroid-plexus 
      MatchValues[nMatch++] = 63; // Right-choroid-plexus 
      //MatchValues[nMatch++] = 15; // 4th-Ventricle. Good to include?
      DoMatch = 1;
    }
    else if (!strcasecmp(option, "--gm")) {
      // Create a mask of all other stuff and invert
      MatchValues[nMatch++] =  2;
      MatchValues[nMatch++] = 41;
      MatchValues[nMatch++] = 77;
      MatchValues[nMatch++] = 251;
      MatchValues[nMatch++] = 252;
      MatchValues[nMatch++] = 253;
      MatchValues[nMatch++] = 254;
      MatchValues[nMatch++] = 255;
      MatchValues[nMatch++] =  7; // Cerebellar WM
      MatchValues[nMatch++] = 46; // Cerebellar WM
      MatchValues[nMatch++] =  4; // Left-Lateral-Ventricle
      MatchValues[nMatch++] =  5; // Left-Inf-Lat-Vent
      MatchValues[nMatch++] = 14; // 3rd-Ventricle
      MatchValues[nMatch++] = 43; // Right-Lateral-Ventricle
      MatchValues[nMatch++] = 44; // Right-Inf-Lat-Vent
      MatchValues[nMatch++] = 15; // 4th-Ventricle 
      MatchValues[nMatch++] = 72; // 5th-Ventricle
      MatchValues[nMatch++] = 31; // Left-choroid-plexus 
      MatchValues[nMatch++] = 63; // Right-choroid-plexus 
      MatchValues[nMatch++] =  0; // Background
      MatchValues[nMatch++] = 24; // CSF
      DoMatch = 1;
      // Invert the matches above
      BinVal = 0;
      BinValNot = 1;
    }
    else if (!strcasecmp(option, "--scm-lh")){
      // Subcortical Mass
      MatchValues[nMatch++] =  26; //Left-Accumbens-area
      MatchValues[nMatch++] =  11; //Left-Caudate
      MatchValues[nMatch++] =  2; //Left-Cerebral-White-Matter
      MatchValues[nMatch++] =  4; //Left-Lateral-Ventricle
      MatchValues[nMatch++] =  13; //Left-Pallidum
      MatchValues[nMatch++] =  12; //Left-Putamen
      MatchValues[nMatch++] =  10; //Left-Thalamus
      MatchValues[nMatch++] =  28; //Left-VentralDC
      MatchValues[nMatch++] =  30; //Left-vessel
      MatchValues[nMatch++] =  34; //Left-WMCrowns
      //MatchValues[nMatch++] =  77; //WM-hypointensities, hopefully fixed by removing holes
      MatchValues[nMatch++] =  78; //Left-WM-hypointensities
      RemoveIslands = 1;
      FillHoles = 1;
      DoMatch = 1;
    }
    else if (!strcasecmp(option, "--scm-rh")){
      // Subcortical Mass
      MatchValues[nMatch++] =  58; //Right-Accumbens-area
      MatchValues[nMatch++] =  50; //Right-Caudate
      MatchValues[nMatch++] =  41; //Right-Cerebral-White-Matter
      MatchValues[nMatch++] =  43; //Right-Lateral-Ventricle
      MatchValues[nMatch++] =  52; //Right-Pallidum
      MatchValues[nMatch++] =  51; //Right-Putamen
      MatchValues[nMatch++] =  49; //Right-Thalamus
      MatchValues[nMatch++] =  60; //Right-VentralDC
      MatchValues[nMatch++] =  62; //Right-vessel
      MatchValues[nMatch++] =  66; //Right-WMCrowns
      //MatchValues[nMatch++] =  77; //WM-hypointensities, hopefully fixed by removing holes
      MatchValues[nMatch++] =  79; //Right-WM-hypointensities
      RemoveIslands = 1;
      FillHoles = 1;
      DoMatch = 1;
    }
    else if (!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      InVolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutVolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--surf")) {
      if (nargc < 1) CMDargNErr(option,1);
      SurfFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--surf-smooth")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nsmoothsurf); // number of iterations
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--merge")) {
      if (nargc < 1) CMDargNErr(option,1);
      MergeVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--copy")) {
      if (nargc < 1) CMDargNErr(option,1);
      CopyVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      MaskVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--mask-thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MaskThresh);
      nargsused = 1;
    } else if (!strcasecmp(option, "--pct")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&TopPercent);
      MinThreshSet = 1;
      DoPercent = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--min")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MinThresh);
      MinThreshSet = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--max")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MaxThresh);
      MaxThreshSet = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--rmin")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&RMinThresh);
      RMinThreshSet = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--rmax")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&RMaxThresh);
      RMaxThreshSet = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--fdr")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&FDR);
      DoFDR = 1;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      int nthreads=1;
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--fdr-pos")) FDRSign = +1;
    else if (!strcasecmp(option, "--fdr-neg")) FDRSign = -1;
    else if (!strcasecmp(option, "--fdr-abs")) FDRSign =  0; //default

    else if (!strcasecmp(option, "--bb")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nPadBB);
      DoBB = 1;
      nargsused = 1;
    }    
    else if (!strcasecmp(option, "--replace")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&SrcReplace[nReplace]);
      sscanf(pargv[1],"%d",&TrgReplace[nReplace]);
      nReplace++;
      nargsused = 2;
    }    
    else if (!strcasecmp(option, "--replaceonly")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&SrcReplace[nReplace]);
      sscanf(pargv[1],"%d",&TrgReplace[nReplace]);
      replace_only = 1 ;
      nReplace++;
      nargsused = 2;
    }    
    else if (!strcasecmp(option, "--replace-nn")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&SrcReplaceNN[nReplaceNN]);
      sscanf(pargv[1],"%d",&ReplaceWindowNN[nReplaceNN]);
      nReplaceNN++;
      nargsused = 2;
    }    
    else if (!strcasecmp(option, "--replaceonly-nn")) {
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&SrcReplaceNN[nReplaceNN]);
      sscanf(pargv[1],"%d",&ReplaceWindowNN[nReplaceNN]);
      replace_only = 1 ;
      nReplaceNN++;
      nargsused = 2;
    }    
    else if (!strcasecmp(option, "--binval")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&BinVal);
      nargsused = 1;
    } else if (!strcasecmp(option, "--binvalnot")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&BinValNot);
      nargsused = 1;
    } else if (!strcasecmp(option, "--inv")) {
      BinVal = 0;
      BinValNot = 1;
    } 
    else if (!strcasecmp(option, "--match")) {
      if (nargc < 1) CMDargNErr(option,1);
      nth = 0;
      while(CMDnthIsArg(nargc, pargv, nth) ){
	sscanf(pargv[nth],"%d",&MatchValues[nMatch]);
	nMatch ++;
	nth++;
      }
      nargsused = nth;
      DoMatch = 1;
    } 
    else if (!strcasecmp(option, "--match-ctab")) {
      if(nargc < 1) CMDargNErr(option,1);
      COLOR_TABLE *ctab = CTABreadASCII(pargv[0]);
      if(ctab==NULL) exit(1);
      int ntotalsegid,n;
      CTABgetNumberOfTotalEntries(ctab,&ntotalsegid);
      for(n=0; n < ntotalsegid; n++){
	int valid;
	CTABisEntryValid(ctab,n,&valid);
	if(!valid) continue;
	MatchValues[nMatch] = n;
	nMatch ++;
      }
      CTABfree(&ctab);
      nargsused = 1;
      DoMatch = 1;
    } 
    else if (!strcasecmp(option, "--subcort-gm")) {
      MatchValues[nMatch] = Left_Thalamus; nMatch++;
      MatchValues[nMatch] = Left_Caudate; nMatch++;
      MatchValues[nMatch] = Left_Putamen; nMatch++;
      MatchValues[nMatch] = Left_Pallidum; nMatch++;
      MatchValues[nMatch] = Left_Hippocampus; nMatch++;
      MatchValues[nMatch] = Left_Amygdala; nMatch++;
      MatchValues[nMatch] = Left_Accumbens_area; nMatch++;
      MatchValues[nMatch] = Left_VentralDC; nMatch++;
      MatchValues[nMatch] = Left_Substancia_Nigra; nMatch++;
      MatchValues[nMatch] = Left_Cerebellum_Cortex; nMatch++;
      MatchValues[nMatch] = Right_Thalamus; nMatch++;
      MatchValues[nMatch] = Right_Caudate; nMatch++;
      MatchValues[nMatch] = Right_Putamen; nMatch++;
      MatchValues[nMatch] = Right_Pallidum; nMatch++;
      MatchValues[nMatch] = Right_Hippocampus; nMatch++;
      MatchValues[nMatch] = Right_Amygdala; nMatch++;
      MatchValues[nMatch] = Right_Accumbens_area; nMatch++;
      MatchValues[nMatch] = Right_VentralDC; nMatch++;
      MatchValues[nMatch] = Right_Substancia_Nigra; nMatch++;
      MatchValues[nMatch] = Right_Cerebellum_Cortex; nMatch++;
      //MatchValues[nMatch] = Brain_Stem; nMatch++;
      DoMatch = 1;
    } 
    else if (!strcasecmp(option, "--frame")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
      DoFrameLoop = 0;
    } 
    else if (!strcasecmp(option, "--frame-sum")) {
      DoFrameSum = 1;
      DoFrameLoop = 0;
    } 
    else if (!strcasecmp(option, "--frame-and")) {
      DoFrameAnd = 1;
      MinThreshSet = 1;
      DoFrameLoop = 0;
    } 
    else if (!strcasecmp(option, "--dilate")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nDilate3d);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--erode-face")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nErodeNN);
      NNType = NEAREST_NEIGHBOR_FACE;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--erode-edge")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nErodeNN);
      NNType = NEAREST_NEIGHBOR_EDGE;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--erode-corner")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nErodeNN);
      NNType = NEAREST_NEIGHBOR_CORNER;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--erode")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nErode3d);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--erode2d")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nErode2d);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--count")) {
      if (nargc < 1) CMDargNErr(option,1);
      CountFile = pargv[0];
      DoCount = 1;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--dilate-vertex-sum") || !strcasecmp(option, "--dilate-vertex")){
      int doradius = 0;
      if(!strcasecmp(option, "--dilate-vertex-sum") && nargc < 5) {
	// vno surf area target outfile
	printf("USAGE: --dilate-vertex-sum vno surf measure target outmask\n");
	exit(1);
      }
      if(!strcasecmp(option, "--dilate-vertex")){
	if(nargc < 4) {
	  // vno surf radius outfile
	  printf("USAGE: --dilate-vertex-sum vno surf radius outmask\n");
	  exit(1);
	}
	doradius = 1;
      }
      int vno; double targetSum;
      sscanf(pargv[0],"%d",&vno);
      MRIS *surf = MRISread(pargv[1]);
      if(!surf) exit(1);
      MRI *measure=NULL;
      if(!doradius && strcmp(pargv[2],"nofile")!=0 && strcmp(pargv[2],"area")!=0) {
	measure = MRIread(pargv[2]);
	if(!measure) exit(1);
      }
      if(doradius) {
	double radius;
	sscanf(pargv[2],"%lf",&radius);
	targetSum = M_PI*radius*radius;
      }
      else         sscanf(pargv[3],"%lf",&targetSum);
      MRI *tmpmri = MRISdilateVertexToSum(vno, surf, measure, targetSum);
      if(!tmpmri) exit(1);
      int err;
      if(doradius) err = MRIwrite(tmpmri,pargv[3]);
      else         err = MRIwrite(tmpmri,pargv[4]);
      exit(err);
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --i invol  : input volume \n");
  printf("   \n");
  printf("   --min min  : min thresh (def is -inf)\n");
  printf("   --max max  : max thresh (def is +inf)\n");
  printf("   --pct P : set threshold to capture top P%% (in mask or total volume)\n");
  printf("   --rmin rmin  : compute min based on rmin*globalmean\n");
  printf("   --rmax rmax  : compute max based on rmax*globalmean\n");
  printf("   --fdr fdrthresh : compute min based on FDR (assuming -log10(p) input)\n");
  printf("     --fdr-pos, --fdr-neg, --fdr-abs (use only pos, neg, or abs; abs is default)\n");
  printf("   --match matchval <matchval2 ...>  : match instead of threshold\n");
  printf("   --match-ctab colortable  : match all entries in the given color table\n");
  printf("   --replace V1 V2 : replace voxels=V1 with V2\n");
  printf("   --replaceonly V1 V2 : replace voxels=V1 with V2 and propagate other src voxels instead of binarizing\n");
  printf("   --replace-nn V1 W : replace voxels=V1 with their nearest neighbor within a window of W voxels\n");
  printf("   --replaceonly-nn V1 W : replace voxels=V1 with their nearest neighbor within a window of W voxels and propagate other src voxels instead of binarizing\n");
  printf("   --ctx-wm : set match vals to 2, 41, 77, 251-255 (aseg for cerebral WM)\n");
  printf("   --all-wm : set match vals to 2, 41, 77, 251-255, 7, and 46, (aseg for all WM)\n");
  printf("   --ventricles : set match vals those for aseg ventricles+choroid (not 4th)\n");
  printf("   --wm+vcsf : WM and ventricular CSF, including choroid (not 4th)\n");
  printf("   --gm : match for all WM and VCSF and background, then invert\n");
  printf("   --subcort-gm : subcortical gray matter\n");
  printf("   --scm-lh, --scm-rh : subcortical mass (includes filling holes and removing islands)\n");
  printf("   \n");
  printf("   --o outvol : output volume \n");
  printf("   --count countfile : save number of hits in ascii file (hits,ntotvox,pct)\n");
  printf("   --no-count : turn off counting number of voxels in the first frame -- good for large volumes\n");
  printf("   \n");
  printf("   --binval    val    : set vox within thresh to val (default is 1) \n");
  printf("   --binvalnot notval : set vox outside range to notval (default is 0) \n");
  printf("   --inv              : set binval=0, binvalnot=1\n");
  printf("   --frame frameno    : use 0-based frame of input (default is 0) \n");
  printf("   --frame-sum : sum frames together before binarizing\n");
  printf("   --frame-and : take intersection (AND) of frames. No --min needed.\n");
  printf("   --merge mergevol   : merge with mergevolume \n");
  printf("   --mask maskvol       : must be within mask \n");
  printf("   --mask-thresh thresh : set thresh for mask (def is 0.5) \n");
  printf("   --abs : take abs of invol first (ie, make unsigned)\n");
  printf("   --bincol : set binarized voxel value to its column number\n");
  printf("   --zero-edges : zero the edge voxels\n");
  printf("   --zero-slice-edges : zero the edge slice voxels\n");
  printf("   --dilate ndilate: dilate binarization in 3D\n");
  printf("   --erode  nerode: erode binarization in 3D (after any dilation)\n");
  printf("   --erode2d nerode2d: erode binarization in 2D (after any 3D erosion)\n");
  printf("   --erode-face   nerode: erode binarization using 'face' nearest neighbors\n");
  printf("   --erode-edge   nerode: erode binarization using 'edge' nearest neighbors\n");
  printf("   --erode-corner nerode: erode binarization using 'corner' nearest neighbors (same as --erode)\n");
  printf("   --remove-islands, --no-remove-islands : remove islands in the mask\n");
  printf("   --fill-holes, --no-fill-holes : remove holes in the mask (after removing islands if specified)\n");
  printf("   --bb npad : reduce dim of output to the minimum volume of non-zero voxels with npad boundary\n");
  printf("   --surf surfname : create a surface mesh from the binarization\n");
  printf("   --surf-smooth niterations : iteratively smooth the surface mesh\n");
  printf("   --threads nthreads (won't apply to replace)\n");
  printf("   --noverbose (default *verbose*) \n");
  printf("   --dilate-vertex vno surf radius outmask : dilate vertex to a target area =pi*R^2 (stand-alone option)\n");
  printf("\n");
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
printf("Program to binarize a volume (or volume-encoded surface file). Can also\n");
printf("be used to merge with other binarizations. Binarization can be done\n");
printf("based on threshold or on matched values.\n");
printf("\n");
printf("--i invol\n");
printf("\n");
printf("Input volume to be binarized.\n");
printf("\n");
printf("--min min\n");
printf("--max max\n");
printf("--rmin rmin\n");
printf("--rmax rmax\n");
printf("\n");
printf("--noverbose (default is *verbose*) \n");
printf("\n");
printf("Minimum and maximum thresholds. If the value at a voxel is >= min and\n");
printf("<= max, then its value in the output will be 1 (or --binval),\n");
printf("otherwise it is 0 (or --binvalnot) or the value of the merge volume at\n");
printf("that voxel.  By default, min = -infinity and max = +infinity, but you\n");
printf("must set one of the thresholds. Cannot be used with --match. If --rmin\n");
printf("or --rmax are specified, then min (or max) is computed as rmin (or\n");
printf("rmax) times the global mean of the input.\n");
printf("\n");
printf("--pct P\n");
printf("\n");
printf("Set min threshold so that the top P percent of the voxels are captured\n");
printf("in the output mask. The percent will be computed based on the number of\n");
printf("voxels in the volume (if not input mask is specified) or within the\n");
printf("input mask.\n");
printf("\n");
printf("--fdr fdrthreshold\n");
printf("\n");
printf("Set min threshold to achieve a given FDR. By default, it uses the\n");
printf("absolute value but this can be changed with --fdr-pos and\n");
printf("--fdr-neg. If a mask is passed, it will compute the voxel-wise\n");
printf("threshold only with in the places where mask > 0.5.  The mask\n");
printf("threshold will be ignored.\n");
printf("\n");
printf("--match matchvalue <matchvalue2 ...>\n");
printf("\n");
printf("Binarize based on matching values. Any number of match values can be \n");
printf("specified. Cannot be used with --min/--max.\n");
printf("\n");
printf("--o outvol\n");
printf("\n");
printf("Path to output volume.\n");
printf("\n");
printf("--count countfile\n");
printf("\n");
printf("Save number of voxels that meet match criteria in ascii\n");
printf("countefile. Four numbers are saved: the number of voxels that match\n");
printf("(nhits), the volume of the voxels that match, the total number of\n");
printf("voxels in the volume (nvoxtot), and the percent matching\n");
printf("(100*nhits/nvoxtot).\n");
printf("\n");
printf("--binval    binval\n");
printf("--binvalnot binvalnot\n");
printf("\n");
printf("Value to use for those voxels that are in the threshold/match\n");
printf("(--binval) or out of the range (--binvalnot). These must be integer\n");
printf("values. binvalnot only applies when a merge volume is NOT specified.\n");
printf("\n");
printf("--replace V1 V2\n");
printf("\n");
printf("Replace every occurrence of (int) value V1 with value V2. Multiple \n");
printf("--replace args are possible.\n");
printf("\n");
printf("--replaceonly V1 V2\n");
printf("\n");
printf("Replace every occurrence of (int) value V1 with value V2. Multiple\n");
printf("--replaceonly args are possible. Other locations in the source volume\n");
printf("will be propagated to the output (unlike --replace which masks those\n");
printf("locations).\n");
printf("\n");
printf("--replace-nn V1 W\n");
printf("\n");
printf("Replace every occurrence of (int) value V1 with that of its nearest\n");
printf("neighbor voxel within a window of W voxels. Multiple --replace-nn args\n");
printf("are possible.\n");
printf("\n");
printf("--replaceonly-nn V1 W\n");
printf("\n");
printf("Replace every occurrence of (int) value V1 with that of its nearest\n");
printf("neighbor voxel within a window of W voxels. Multiple --replaceonly-nn\n");
printf("args are possible. Other locations in the source volume will be propagated\n");
printf("to the output (unlike --replace-nn which masks those locations).\n");
printf("\n");
printf("--frame frameno\n");
printf("\n");
printf("Use give frame of the input. 0-based. Default is 0.\n");
printf("\n");
printf("--frame-sum\n");
printf("\n");
printf("Sum the frames together before applying threshold.\n");
printf("\n");
printf("--frame-and\n");
printf("\n");
printf("Treat the multi-frame volume as binary 'AND' the frames together. This\n");
printf("takes an intersection of the individual frames. You do not need to \n");
printf("specify a --min (the min will be set to nframes-0.5).\n");
printf("\n");
printf("--merge mergevol\n");
printf("\n");
printf("Merge binarization with the mergevol. If the voxel is within the threshold\n");
printf("range (or matches), then its value will be binval. If not, then it will \n");
printf("inherit its value from the value at that voxel in mergevol. mergevol must \n");
printf("be the same dimension as the input volume. Combining this with --binval \n");
printf("allows you to construct crude segmentations.\n");
printf("\n");
printf("--mask maskvol\n");
printf("--mask-thresh thresh\n");
printf("\n");
printf("Mask input with mask. The mask volume is itself binarized at thresh\n");
printf("(default is 0.5). If a voxel is not in the mask, then it will be assigned\n");
printf("binvalnot or the value from the merge volume.\n");
printf("\n");
printf("--zero-edges\n");
printf("\n");
printf("Set the first and last planes in all dimensions to 0 (or --binvalnot). This\n");
printf("makes sure that all the voxels on the edge of the imaging volume are 0.\n");
printf("\n");
printf("--zero-slice-edges\n");
printf("\n");
printf("Same as --zero-edges, but only for slices.\n");
printf("\n");
printf("--uchar\n");
printf("\n");
printf("Save output file in 'uchar' format.\n");
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
  if (InVolFile == NULL) {
    printf("ERROR: must specify input volume\n");
    exit(1);
  }
  if(OutVolFile == NULL && CountFile == NULL && SurfFile == NULL) {
    printf("ERROR: must specify output volume or output count file\n");
    exit(1);
  }
  if(MinThreshSet == 0  && MaxThreshSet == 0 &&
     RMinThreshSet == 0 && RMaxThreshSet == 0 &&
     !DoMatch && !DoFDR) {
    if(nReplace > 0 || nReplaceNN > 0)
      replace_only = 1;
    else {
      printf("ERROR: must specify minimum and/or maximum threshold or match values\n");
      exit(1);
    }
  }
  if((MinThreshSet || MaxThreshSet || RMinThreshSet || RMaxThreshSet) && DoMatch ) {
    printf("ERROR: cannot specify threshold and match values\n");
    exit(1);
  }
  if(MinThreshSet && RMinThreshSet){
    printf("ERROR: cannot --rmin and --min \n");
    exit(1);
  }
  if(MaxThreshSet && RMaxThreshSet){
    printf("ERROR: cannot --rmax and --max \n");
    exit(1);
  }
  if(!DoMatch){
    if (MaxThreshSet && MinThreshSet && MaxThresh < MinThresh) {
      printf("ERROR: max thresh = %g < min thresh = %g\n",
	     MaxThresh,MinThresh);
      exit(1);
    }
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  int n;

  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"\n");
  fprintf(fp,"input      %s\n",InVolFile);
  fprintf(fp,"frame      %d\n",frame);
  fprintf(fp,"nErode3d   %d\n",nErode3d);
  fprintf(fp,"nErode2d   %d\n",nErode2d);
  fprintf(fp,"output     %s\n",OutVolFile);

  if(!DoMatch){
    fprintf(fp,"Binarizing based on threshold\n");
    if(MinThreshSet)        fprintf(fp,"min        %g\n",MinThresh);
    else if(RMinThreshSet)  fprintf(fp,"rmin       %g\n",RMinThresh);
    else                    fprintf(fp,"min        -infinity\n");
    if(MaxThreshSet)        fprintf(fp,"max        %g\n",MaxThresh);
    else if(RMaxThreshSet)  fprintf(fp,"rmax       %g\n",RMaxThresh);
    else                    fprintf(fp,"max        +infinity\n");
  }
  else {
    fprintf(fp,"Binarizing based on matching values\n");
    fprintf(fp,"nMatch %d\n",nMatch);
    for(n=0; n < nMatch; n++)
      fprintf(fp,"%2d  %4d\n",n,MatchValues[n]);
  }
  fprintf(fp,"binval        %d\n",BinVal);
  fprintf(fp,"binvalnot     %d\n",BinValNot);
  if (MergeVolFile)
    fprintf(fp,"merge      %s\n",MergeVolFile);
  if (MaskVolFile) {
    fprintf(fp,"mask       %s\n",MaskVolFile);
    fprintf(fp,"maskthresh %lf\n",MaskThresh);
  }
  return;
}
