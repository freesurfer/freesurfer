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

#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri_identify.h"
#include "matrix.h"
#include "registerio.h"
#include "resample.h"

#include "volcluster.h"

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
        char *thsign, int invert, 
        int lowval, int highval, MRI *binvol);


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static int  singledash(char *flag);
static void dump_options(FILE *fp);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_volcluster.c,v 1.4 2002/04/02 19:53:40 greve Exp $";
char *Progname = NULL;

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
char  *masksignstring = "abs";
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
char  *signstring = "abs";
int   threshsign =    0;
float sizethresh    = -1.0;
int   sizethreshvox = -1;
float distthresh =   0.0;
int   allowdiag  = 0;

MRI *vol, *HitMap, *outvol, *maskvol, *binmask;
VOLCLUSTER **ClusterList, **ClusterList2;
MATRIX *CRS2MNI, *CRS2FSA, *FSA2Func;
LABEL *label;

float colres, rowres, sliceres, voxsize;

FILE *fpsum;

/*--------------------------------------------------------------*/
/*--------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
  int nhits, *hitcol, *hitrow, *hitslc;
  int col, row, slc;
  int nthhit, n, nclusters, nprunedclusters;
  float x,y,z;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  if(debug) dump_options(stdout);

  /* Load the input volume */
  vol = MRIreadType(volid,intype);
  if(vol == NULL){
    fprintf(stderr,"ERROR: reading %s\n",volid);
    exit(1);
  }

  /* Check that the specified frame is within range */
  if(frame >= vol->nframes){
    fprintf(stderr,"ERROR: specified frame = %d, >= nframes\n",frame);
    exit(1);
  }

  /* Load the mask volume */
  if(maskid != NULL){
    printf("INFO: loading mask volume: %s\n",maskid);
    maskvol = MRIreadType(maskid,masktype);
    if(maskvol == NULL){
      fprintf(stderr,"ERROR: reading %s\n",maskid);
      exit(1);
    }
    if(vol->width != maskvol->width && 
       vol->depth != maskvol->depth &&
       vol->height != maskvol->height){
      printf("ERROR: input volume and mask volume dimensions are "
       "inconsistent\n");
      exit(1);
    }
    if(maskframe >= maskvol->nframes){
      printf("ERROR: mask frame = %d >= nframes in mask (%d)\n",
       maskframe,maskvol->nframes);
      exit(1);

    }

    binmask = MRIbinarize01(maskvol, maskthresh, -1, masksignstring, 
          maskinvert, 0, 1, NULL);
    if(binmask == NULL) exit(1);

    if(outmaskid != NULL) MRIwriteType(binmask,outmaskid,outmasktype);
    MRIfree(&maskvol);
  }
  else binmask = NULL;

  /* Load the resolution and geometry information from the register.dat */
  CRS2MNI = LoadMNITransform(regfile, vol->width,vol->height,vol->depth,
           &CRS2FSA, &FSA2Func, 
           &colres, &rowres, &sliceres);
  voxsize = colres * rowres * sliceres;
  if(debug){
    printf("VolumeRes: %g %g %g (%g)\n",colres,rowres,sliceres,voxsize);
    printf("Registration: ---------------\n");
    MatrixPrint(stdout,FSA2Func);
    printf("CRS2FSA: ---------------\n");
    MatrixPrint(stdout,CRS2FSA);
    printf("CRS2MNI: ---------------\n");
    MatrixPrint(stdout,CRS2MNI);
  }

  if(sizethresh < 0) sizethresh = sizethreshvox*voxsize;

  /* Replace data with synthetic if desired */
  if(synthfunction != NULL){
    printf("INFO: synthsizing with %s\n",synthfunction);
    if(!strcmp(synthfunction,"uniform")) 
      MRIsynthUniform(0,0,0,0,vol);
    if(!strcmp(synthfunction,"loguniform")) 
      MRIsynthLogUniform(0,0,0,0,vol);
    if(!strcmp(synthfunction,"gaussian")) 
      MRIsynthGaussian(0,0,0,0,vol);
  }

  /* Initialize the hit map - this is a map of voxels that have been
     accounted for as either outside of a the threshold range or
     belonging to a cluster. The initialization is for thresholding */
  HitMap = clustInitHitMap(vol, frame, threshmin, threshmax, threshsign, 
         &nhits, &hitcol, &hitrow, &hitslc, 
         binmask, maskframe);
  if(HitMap == NULL) exit(1);
  //MRIwriteType(HitMap,"hitmap",BSHORT_FILE);

  printf("INFO: Found %d voxels in threhold range\n",nhits);

  /* Allocate an array of clusters equal to the number of hits -- this
     is the maximum number of clusters possible */
  ClusterList = clustAllocClusterList(nhits);
  if(ClusterList == NULL){
    fprintf(stderr,"ERROR: could not alloc %d clusters\n",nhits);
    exit(1);
  }

  nclusters = 0;
  for(nthhit = 0; nthhit < nhits; nthhit ++){

    /* Determine whether this hit is still valid. It may 
       not be if it was assigned to the cluster of a 
       previous hit */
    col = hitcol[nthhit];
    row = hitrow[nthhit];
    slc = hitslc[nthhit];
    if(MRIIseq_vox(HitMap,col,row,slc,0)) continue;

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

  if(debug)
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

  if(debug)
    clustDumpClusterList(stdout,ClusterList, nclusters, vol, frame);

  /* Remove clusters that do not meet the minimum distance requirement */
  if(distthresh > 0.0){
    printf("INFO: pruning by distance %g\n",distthresh);
    ClusterList2 = clustPruneByDistance(ClusterList,nclusters,
          distthresh, &nprunedclusters);
    //clustFreeClusterList(&ClusterList,nclusters);/* Free - does not work */
    nclusters = nprunedclusters;
    ClusterList = ClusterList2;
  }

  /* Sort Clusters by MaxValue */
  ClusterList2 = clustSortClusterList(ClusterList2,nclusters,NULL);
  //clustFreeClusterList(&ClusterList,nclusters);/* Free - does not work */
  ClusterList = ClusterList2;

  printf("INFO: Found %d final clusters\n",nclusters);
  if(debug)
    clustDumpClusterList(stdout,ClusterList, nclusters, vol, frame);

  /* Open the Summary File (or set its pointer to stdout) */
  if(sumfile != NULL){
    fpsum = fopen(sumfile,"w");
    if(fpsum == NULL){
      printf("ERROR: could not open %s for writing\n",sumfile);
      exit(1);
    }
    printf("INFO: writing summary to %s\n",sumfile);
  }
  else fpsum = stdout;

  /* Dump summary to file or stdout */
  fprintf(fpsum,"Cluster Growing Summary (mri_volcluster)\n");
  fprintf(fpsum,"Input Volume:      %s\n",volid);  
  fprintf(fpsum,"Frame Number:      %d\n",frame);  
  fprintf(fpsum,"Minimum Threshold: %g\n",threshmin);  
  if(threshmax < 0) 
    fprintf(fpsum,"Maximum Threshold: inifinity\n");
  else
    fprintf(fpsum,"Maximum Threshold: %g\n",threshmax);
  fprintf(fpsum,"Threshold Sign:    %s\n",signstring);  

  if(distthresh > 0) 
    fprintf(fpsum,"Distance Threshold: %g (mm)\n",distthresh);

  fprintf(fpsum,"Size Threshold:    %g mm^3\n",sizethresh);  
  fprintf(fpsum,"Size Threshold:    %g voxels\n",sizethresh/voxsize);  
  fprintf(fpsum,"Voxel Size:        %g mm^3\n",voxsize);  
  fprintf(fpsum,"Registration:      %s\n",regfile);  
  if(synthfunction != NULL)
    fprintf(fpsum,"Synthesize:        %s\n",synthfunction);  
  if(maskid != NULL){
    fprintf(fpsum,"Mask Vol:          %s\n",maskid);  
    fprintf(fpsum,"Mask Thresh:       %f\n",maskthresh);  
    fprintf(fpsum,"Mask Sign:         %s\n",masksignstring);  
    fprintf(fpsum,"Mask Invert:       %d\n",maskinvert);  
  }
  fprintf(fpsum,"AllowDiag:         %d\n",allowdiag);    
  fprintf(fpsum,"NClusters          %d\n",nclusters);  

  fprintf(fpsum,"\n");  
  fprintf(fpsum,"Cluster   Size(n)   Size(mm^3)     TalX   TalY    TalZ              Max\n");  

  for(n = 0; n < nclusters; n++){
    clustComputeTal(ClusterList[n],CRS2MNI); /* for "true" Tal coords */
    //clustComputeXYZ(ClusterList[n],CRS2FSA); /* for FSA coords */
    //clustComputeXYZ(ClusterList[n],CRS2MNI); /* for MNI coords */
    col = ClusterList[n]->col[ClusterList[n]->maxmember];
    row = ClusterList[n]->row[ClusterList[n]->maxmember];
    slc = ClusterList[n]->slc[ClusterList[n]->maxmember];
    x = ClusterList[n]->x[ClusterList[n]->maxmember];
    y = ClusterList[n]->y[ClusterList[n]->maxmember];
    z = ClusterList[n]->z[ClusterList[n]->maxmember];
    fprintf(fpsum,"%3d        %4d      %7.1f    %7.2f %7.2f %7.2f   %15.5f",
     n+1,ClusterList[n]->nmembers,voxsize*ClusterList[n]->nmembers,
      x,y,z, ClusterList[n]->maxval);
    if(debug)
      fprintf(fpsum,"  %3d %3d %3d \n",col,row,slc);
    else
      fprintf(fpsum,"\n");
  }
  if(sumfile != NULL) fclose(fpsum);

  /* Write clusters values to a volume */
  if(outid != 0){
    outvol = clustClusterList2Vol(ClusterList, nclusters, vol,frame, 1);
    MRIwriteType(outvol,outid,outtype);
    MRIfree(&outvol);
  }

  /* Write clusters numbers to a volume */
  if(outcnid != 0){
    outvol = clustClusterList2Vol(ClusterList, nclusters, vol,frame, 0);
    printf("INFO: writing OCN to %s as type %d\n",outcnid,outtype);
    MRIwriteType(outvol,outcnid,outcntype);
    MRIfree(&outvol);
  }

  /* Write the given cluster to a label file */
  /* Warning: the cluster xyz will change */
  if(labelfile != NULL){
    if(nlabelcluster > nclusters){
      fprintf(stderr,"ERROR: selected cluster number %d, "
        "but there are only %d clusters\n",nlabelcluster,nclusters+1);
      exit(1);
    }
    printf("Computing Label\n");
    label = clustCluster2Label(ClusterList[nlabelcluster-1], vol, frame,
             colres, rowres, sliceres, FSA2Func);
    LabelWrite(label, labelfile);
  }

  if(labelbase != NULL){
    for(nlabelcluster = 0; nlabelcluster < nclusters; nlabelcluster ++){
      
      printf("Computing label for cluster %d\n",nlabelcluster);
      sprintf(tmpstr,"%s-%04d.label",labelbase,nlabelcluster+1);
      label = clustCluster2Label(ClusterList[nlabelcluster], vol, frame,
         colres, rowres, sliceres, FSA2Func);
      printf("Saving label to %s\n",tmpstr);
      LabelWrite(label, tmpstr);
      LabelFree(&label);
    }
  }

  printf("mri_volcluster: done\n");

  return(0);
  exit(0);
}

/* ----------------------------------------------------------- */
/* --------------->>>>>>.<<<<<<<------------------------------ */
/* ----------------------------------------------------------- */
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
    else if (!strcasecmp(option, "--verbose")) verbose = 1;

    else if (!strcasecmp(option, "--allowdiag")) allowdiag = 1;
    else if (!strcmp(option, "--maskinvert"))    maskinvert = 1;

    /* -------- source volume inputs ------ */
    else if (!strcmp(option, "--i") || !strcmp(option, "--in")){
      if(nargc < 1) argnerr(option,1);
      volid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--in_type")){
      if(nargc < 1) argnerr(option,1);
      intypestring = pargv[0];
      intype = string_to_type(intypestring);
      nargsused = 1;
    }
    else if (!strcmp(option, "--mask")){
      if(nargc < 1) argnerr(option,1);
      maskid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--mask_type")){
      if(nargc < 1) argnerr(option,1);
      masktypestring = pargv[0];
      masktype = string_to_type(masktypestring);
      nargsused = 1;
    }
    else if (!strcmp(option, "--maskframe")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&maskframe);
      if(maskframe < 0){
  fprintf(stderr,"ERROR: negative frame number: frame = %d\n",maskframe);
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--masksign")){
      if(nargc < 1) argnerr(option,1);
      masksignstring = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--maskthresh")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&maskthresh);
      if(maskthresh < 0.0){
  fprintf(stderr,"ERROR: negative mask threshold not"
    "allowed (use -masksign)\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--outmask")){
      if(nargc < 1) argnerr(option,1);
      outmaskid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--outmask_type")){
      if(nargc < 1) argnerr(option,1);
      outmasktypestring = pargv[0];
      outmasktype = string_to_type(outmasktypestring);
      nargsused = 1;
    }
    else if (!strcmp(option, "--synth")){
      if(nargc < 1) argnerr(option,1);
      synthfunction = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--o") || !strcmp(option, "--out")){
      if(nargc < 1) argnerr(option,1);
      outid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--out_type")){
      if(nargc < 1) argnerr(option,1);
      outtypestring = pargv[0];
      outtype = string_to_type(outtypestring);
      nargsused = 1;
    }
    else if (!strcmp(option, "--ocn")){
      if(nargc < 1) argnerr(option,1);
      outcnid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--ocn_type")){
      if(nargc < 1) argnerr(option,1);
      outcntypestring = pargv[0];
      outcntype = string_to_type(outcntypestring);
      nargsused = 1;
    }
    else if (!strcmp(option, "--sum")){
      if(nargc < 1) argnerr(option,1);
      sumfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--labelfile") ||
       !strcmp(option, "--label")){
      if(nargc < 1) argnerr(option,1);
      labelfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--labelbase") ||
       !strcmp(option, "--labelbase")){
      if(nargc < 1) argnerr(option,1);
      labelbase = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--nlabelcluster")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nlabelcluster);
      if(nlabelcluster <= 0){
  fprintf(stderr,"ERROR: non-postive nlabelcluster not allowed\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--thmin")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&threshmin);
      if(threshmin < 0.0){
  fprintf(stderr,"ERROR: negative threshold not allowed (use -sign)\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--thmax")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&threshmax);
      if(threshmax < 0.0){
  fprintf(stderr,"ERROR: negative threshold not allowed (use -sign)\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--sign")){
      if(nargc < 1) argnerr(option,1);
      signstring = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--frame")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      if(frame < 0){
  fprintf(stderr,"ERROR: negative frame number: frame = %d\n",frame);
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--minsize")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&sizethresh);
      if(sizethresh <= 0){
  fprintf(stderr,"ERROR: negative cluster size threshold not allowed\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--minsizevox")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&sizethreshvox);
      if(sizethreshvox <= 0){
  fprintf(stderr,"ERROR: negative cluster size threshold not allowed\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--mindist")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&distthresh);
      if(distthresh <= 0){
  fprintf(stderr,"ERROR: negative distance threshold not allowed\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--reg")){
      if(nargc < 1) argnerr(option,1);
      regfile = pargv[0];
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
  fprintf(stdout, "USAGE: %s \n",Progname) ;
  fprintf(stdout, "\n");
  fprintf(stdout, "   --in      input volid \n");
  fprintf(stdout, "   --in_type file format \n");
  fprintf(stdout, "   --frame   frameno <0>\n");
  fprintf(stdout, "   --reg     register.dat\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --thmin   minthresh\n");
  fprintf(stdout, "   --thmax   maxthresh (default is infinity)\n");
  fprintf(stdout, "   --sign    <abs>, neg, pos\n");
  fprintf(stdout, "   --minsize    minimum volume (mm^3)\n");
  fprintf(stdout, "   --minsizevox minimum volume (voxels)\n");
  fprintf(stdout, "   --mindist distance threshold <0>\n");
  fprintf(stdout, "   --allowdiag  : define contiguity to include diagonal\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --mask      mask volid (same dim as input)\n");
  fprintf(stdout, "   --mask_type file format \n");
  fprintf(stdout, "   --maskframe frameno <0> \n");
  fprintf(stdout, "   --maskthresh  upper threshold \n");
  fprintf(stdout, "   --masksign   <abs>, neg, pos \n");
  fprintf(stdout, "   --maskinvert  \n");
  fprintf(stdout, "   --outmask      final binary mask\n");
  fprintf(stdout, "   --outmask_type file format \n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --sum file   : text summary file \n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --out      outupt volid \n");
  fprintf(stdout, "   --out_type file format \n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --ocn      output cluster number volid \n");
  fprintf(stdout, "   --ocn_type file format \n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --label   label file\n");
  fprintf(stdout, "   --nlabelcluster n : save nth cluster in label\n");
  fprintf(stdout, "   --labelbase  base : save all clusters under base-NNNN.label\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --synth   synthfunc (uniform,loguniform,gaussian)\n");
  fprintf(stdout, "   --help    : how to use this program \n");
  fprintf(stdout, "\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  printf("\n");
  print_usage() ;
  printf("\n");

  printf("

DESCRIPTION:

This program will find clusters in a volume. A cluster is a set of
contiguous voxels which meet a threshold criteria. The set may also
have to reach a certain minimum number of voxels to be considered a
cluster. The results can be saved in four ways: (1) a text summary
file, (2) a new volume which is same as the input volume but with all
the voxels that do not belong to a cluster set to zero, (3) a volume
with each voxel's value equal to the cluster number in the summary
file to which the voxel belongs, and (4) one cluster can be saved as a
label file. The search space within the volume can be restricted to be
within a mask. Two voxels are considered contiguous if they share a 
common row, column, or slice (except for --allowdiag).

COMMAND-LINE ARGUMENTS:

  --in input-volid : path/name of the input volume to be clustered.

  --in_type typename : name of the file format type of the input
    volume. See FILE TYPES below.

  --frame n : perform the cluster analysis on the nth frame. The first
    frame corresponds to n=0.

  --reg registration-file : registration file created by either tkregister
    or tkmedit. Note: the subject listed in the registration file must 
    have an entry under $SUBJECTS_DIR.

  --thmin minthresh : voxel values must exceed this amount to be considered
    as a candidate for a cluster.

  --thmax maxthresh : the value of a voxel must be less than this
    amount to be considered as a candidate for a cluster. Default is
    infinity.

  --sign sign-name : the value of a voxel must be this sign to be considered
    as a candidate for a cluster. pos = positive, neg = negative, abs = 
    absolute (ie, ignore sign). Default is abs. When neg is used, the 
    interpretation of minthreshold and maxthreshold are reversed.

  --minsize volume : minimum volume (in mm^3) of contiguous voxels that
    meet the threshold criteria needed to become a cluster. See also
    --minsizevox.

  --minsizevox number : minimum number of contiguous voxels that meet
    the threshold criteria needed to become a cluster. See also --minsize.

  --mindistance distance : minimum distance (in mm) that the peaks of two
    clusters must be separated in order for them to be considered two
    distinct clusters. If two clusters are closer than this amount, the
    cluster with the lower peak is eliminated.

  --allowdiag : (no argument) allow two voxels that share a corner to 
    be considered contiguous.

  --mask mask-volid: path/name of a mask used to restrict the region
    over which clusters will be searched. For example, this could be used
    to restrict the search space to only the brain (ie, exclude eyeballs).
    The mask must have the same dimension as the input volume.

  --mask_type typename : name of the file format type of the mask
    volume. See FILE TYPES below.

  --maskframe n : use the nth frame of the mask volume as the mask,
    where n = 0 indicates the first frame. Default is 0.

  --maskthresh maskthresh: use only those voxels in the mask whose
    value exceeds the given threshold (see also --masksign and
    --maskinvert).

  --masksign sign-name : threshold the mask based on the sign of the
    value at each mask voxel: pos = positive, neg = negative, abs = 
    absolute (ie, ignore sign). Default is abs. When neg is used, the 
    a mask voxel value must be less than -maskthresh.

  --maskinverse : after determining which voxels are in the mask based
    on the threshold and sign, take only voxels that are outside of the
    mask.

  --outmask outmask-volid: path/name of the final binary mask after
    considering thresholding, sign, and inversion. This is mainly useful
    for debuggin.

  --outmask_type typename : name of the file format type of the outmask
    volume. See FILE TYPES below.

  --sumfile sumfilename : save a summary of the results in ASCII into 
    sumfilename. See SUMMARY FILE below.

  --out out-volid: path/name of the output volume after
    clustering. All voxels that were not part of a cluster have their
    values set to zero. Otherwise, their values do not change.

  --out_type typename : name of the file format type of the output 
    volume. See FILE TYPES below.

  --ocn ocn-volid: path/name of the output volume after clustering
    where the value at each voxel is the cluster number (as found in the
    summary file). Voxels that did not belong to a cluster have their
    values set to zero.

  --ocn_type typename : name of the file format type of the output 
    cluster number volume. See FILE TYPES below.

  --label label-file : save the nth cluster (see -nlabelcluster) as a
    label file in the subject's anatomical space. Note: make sure that the
    label file includes a forward slash (/) or the label will be saved
    into the subjects anatomical direcotry. For example: ./mylabel.label.
    Requires --nlabelcluster.

  --nlabelcluster n : save the nth cluster (see -label) as a label file.

  --labelbase base : save each cluster in its own label file. The name
    of the file will be base-NNNN.label, where NNNN is the four digit,
    zero-padded cluster number. All clusters found will be saved.

SUMMARY FILE

The summary file produced by --sumfile is in ASCII text. The summary
file will have a short header indicating the conditions underwhich the
clustering was performed followed by rows of data with 7 columns. Each
row reprsents a different cluster. Each column has the following
interpretation: (1) cluster number, (2) number of voxels in the
cluster, (3) cluster volume in mm^3, (4-6) Talairach X, Y, and Z (mm),
(7) maximum value inside the cluster. The Talairach coordinates are
the 'true' coordinates (not MNI). Part of a summary file is 
shown below as an example:

-----------------------------------------------------------------------
Cluster Growing Summary (mri_volcluster)
Input Volume:      grp-subj-p1p2vlor/bold/sxa-h20/tal-ffx-rfx/novvfix/sig
Frame Number:      0
Minimum Threshold: 2
Maximum Threshold: inifinity
Threshold Sign:    abs
Distance Threshold: 10 (mm)
Size Threshold:    640 mm^3
Size Threshold:    10 voxels
Voxel Size:        64 mm^3
Registration:      grp-subj-p1p2vlor/bold/sxa-h20/tal-ffx-rfx/register.dat
Mask Vol:          talbrain/mask
Mask Thresh:       0.500000
Mask Sign:         pos
Mask Invert:       0
AllowDiag:         1
NClusters          26

Cluster Size(n) Size(mm^3)  TalX   TalY    TalZ     Max
  1      348      22272.0   59.40 -66.72  -13.48   5.66192
  2       45       2880.0  -39.60  26.79   -8.07   5.45487
  3       27       1728.0   55.44  16.60   21.28   4.95684
-----------------------------------------------------------------------


FILE TYPES/FORMATS:

mri_volcluster can read/write any file format that can be read/written
by mri_convert. The most relevent ones are: bshort, bfloat, analyze,
analyze4d. When specifying bshort/bfloat, the volume id is the
stem. Ie, if the volume is f_000.bshort, f_001.bshort, ..., then the
volume id is 'f' (no quotes).

KNOWN BUGS:

When specifying a label file, make sure that the label file includes a
forward slash (/) or the label will be saved into the subjects
anatomical direcotry. For example: ./mylabel.label.

BUG REPORTS:

Send bug reports to analysis-bugs@nmr.mgh.harvard.edu. Make sure to
include the version number (--version) , the full command-line used,
the type of computer operating system, anything printed to the screen,
a description of what you think is wrong and why, and anything else
that may be helpful. Users at the NMR Center should also include the
directory that they ran it from. NOTE: bug reports that do not have
sufficient information to diagnose the problem will probably be either
ignored or placed at the bottom of the list. BE CLEAR!

AUTHOR:

Douglas N. Greve, Ph.D; NMR Center, MGH, greve@nmr.mgh.harvard.edu

") ;

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
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
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/* --------------------------------------------- */
static void check_options(void)
{
  int err;

  err = 0;

  if(volid == NULL){
    fprintf(stderr,"ERROR: no input volume supplied\n");
    err = 1;
  }

  if(threshmin < 0) {
    fprintf(stderr,"ERROR: no minimum threshold supplied\n");
    err = 1;
  }

  if(sizethresh <= 0.0 && sizethreshvox <= 0) {
    fprintf(stderr,"ERROR: no cluster size threshold supplied\n");
    err = 1;
  }

  if(      !strncmp(signstring,"a",1) ) threshsign = 0;
  else if( !strncmp(signstring,"p",1) ) threshsign = +1;
  else if( !strncmp(signstring,"n",1) ) threshsign = -1;
  else {
    fprintf(stderr,"ERROR: sign = %s, must be neg, abs, or pos\n",signstring);
    err = 1;
  }

  if(synthfunction != NULL){
    if(strcmp(synthfunction,"uniform")    && 
       strcmp(synthfunction,"loguniform") &&
       strcmp(synthfunction,"gaussian") ) {
      fprintf(stderr,
        "ERROR: synth = %s, must be uniform, loguniform, or gaussian\n",
        synthfunction);
      err = 1;
    }
  }

  if(regfile == NULL){
    fprintf(stderr,"ERROR: must specify a registration file\n");
    err = 1;
  }

  if(intype == MRI_VOLUME_TYPE_UNKNOWN) intype = mri_identify(volid);
  if(intype == MRI_VOLUME_TYPE_UNKNOWN){
    fprintf(stderr,"ERROR: could not determine type of %s\n",volid);
    exit(1);
  }

  if(outid != 0){
    if(outtype == MRI_VOLUME_TYPE_UNKNOWN) outtype = mri_identify(volid);
    if(outtype == MRI_VOLUME_TYPE_UNKNOWN){
      fprintf(stderr,"ERROR: could not determine type of %s\n",outid);
      exit(1);
    }
  }

  if(maskid != 0){
    if(masktype == MRI_VOLUME_TYPE_UNKNOWN) masktype = mri_identify(volid);
    if(masktype == MRI_VOLUME_TYPE_UNKNOWN){
      fprintf(stderr,"ERROR: could not determine type of %s\n",maskid);
      exit(1);
    }
  }

  if(outmaskid != 0){
    if(outmasktype == MRI_VOLUME_TYPE_UNKNOWN) 
      outmasktype = mri_identify(volid);
    if(outmasktype == MRI_VOLUME_TYPE_UNKNOWN){
      fprintf(stderr,"ERROR: could not determine type of %s\n",outmaskid);
      exit(1);
    }
  }

  if(outcnid != 0){
    if(outcntype == MRI_VOLUME_TYPE_UNKNOWN) 
      outcntype = mri_identify(volid);
    if(outcntype == MRI_VOLUME_TYPE_UNKNOWN){
      fprintf(stderr,"ERROR: could not determine type of %s\n",outcnid);
      exit(1);
    }
  }

  if(labelfile != NULL && nlabelcluster < 1){
    printf("ERROR: --nlabelcluster must be specified with --label\n");
    err = 1;
  }


  if(err) exit(1);

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"volid    %s\n",volid);
  if(intypestring) fprintf(fp,"voltype  %s\n",intypestring);
  fprintf(fp,"frame    %d\n",frame);
  fprintf(fp,"thmin    %g\n",threshmin);
  fprintf(fp,"thmax    %g\n",threshmax);
  fprintf(fp,"sign     %s\n",signstring);
  fprintf(fp,"sizethresh  %g\n",sizethresh);
  fprintf(fp,"distthresh  %g\n",distthresh);
}

/* ------------------------------------------------------------------
   LoadMNITransform() - returns the matrix that transforms col, row,
   and slice in the input volume to MNI x, y, z.  It also returns the
   matrix that transforms col, row, and slice in the input volume to
   FreeSurfer Anatomical x, y, z.
   --------------------------------------------------------------- */
static MATRIX *LoadMNITransform(char *regfile, int ncols, int nrows, 
        int nslices, MATRIX **ppCRS2FSA,
        MATRIX **ppFSA2Func,
        float *colres, float *rowres, float *sliceres)
{
  char *SUBJECTS_DIR;
  int err;
  char *subject;
  float ipr, bpr, intensity;
  MATRIX *R, *iR, *T, *Q, *iQ;
  MATRIX *CRS2MNI;
  int float2int;
  char talxfmfile[1000];
    
  err = regio_read_register(regfile, &subject, &ipr, &bpr, 
            &intensity, &R, &float2int);
  if(err) exit(1);
  iR = MatrixInverse(R,NULL);

  /* get the SUBJECTS_DIR environment variable */
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR==NULL){
   printf("ERROR: environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(1);
  }
  printf("INFO: subject = %s\n",subject);

  /* get the MINC talairach.xfm file name */
  sprintf(talxfmfile,"%s/%s/mri/transforms/talairach.xfm",
    SUBJECTS_DIR,subject);
  err = regio_read_mincxfm(talxfmfile, &T);
  if(err) exit(1);

  Q = FOVQuantMatrix(ncols, nrows, nslices, ipr, ipr, bpr); 
  iQ = MatrixInverse(Q,NULL);

  *ppCRS2FSA = MatrixMultiply(iR,iQ,NULL);
  CRS2MNI = MatrixMultiply(T,*ppCRS2FSA,NULL);

  MatrixFree(&iR);
  MatrixFree(&Q);
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
          int nframes, MRI *tvol)
{
  MRI *vol;
  int col, row, slc, frm;

  if(tvol == NULL){
    vol = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if(vol == NULL){
      fprintf(stderr,"ERROR: MRIsynthUniform: could not alloc mri\n");
      return(NULL);
    }
  }
  else vol = tvol;

  for(col = 0; col < vol->width;   col++){
    for(row = 0; row < vol->height;  row++){
      for(slc = 0; slc < vol->depth;   slc++){
  for(frm = 0; frm < vol->nframes; frm++){
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
             int nframes, MRI *tvol)
{
  MRI *vol;
  int col, row, slc, frm;

  if(tvol == NULL){
    vol = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if(vol == NULL){
      fprintf(stderr,"ERROR: MRIsynthLogUniform: could not alloc mri\n");
      return(NULL);
    }
  }
  else vol = tvol;

  for(col = 0; col < vol->width;   col++){
    for(row = 0; row < vol->height;  row++){
      for(slc = 0; slc < vol->depth;   slc++){
  for(frm = 0; frm < vol->nframes; frm++){
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
           int nframes, MRI *tvol)
{
  MRI *vol;
  int col, row, slc, frm;

  if(tvol == NULL){
    vol = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if(vol == NULL){
      fprintf(stderr,"ERROR: MRIsynthGaussian: could not alloc mri\n");
      return(NULL);
    }
  }
  else vol = tvol;

  for(col = 0; col < vol->width;   col++){
    for(row = 0; row < vol->height;  row++){
      for(slc = 0; slc < vol->depth;   slc++){
  for(frm = 0; frm < vol->nframes; frm++){
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
static double Gaussian01PDF(void)
{
  double v1,v2,r2;

  do
    {
      v1 = 2.0 * drand48() - 1.0;
      v2 = 2.0 * drand48() - 1.0;
      r2 = v1*v1 + v2*v2;
    } 
  while( r2 > 1.0);

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
        char *thsign, int invert, 
        int lowval, int highval, MRI *binvol)
{
  int ncols, nrows, nslices, nframes;
  int col, row, slice, frame;
  int ithsign, nhits;
  short r;
  float val=0.0;

  if(strncasecmp(thsign,"neg",3)==0)      ithsign = -1;
  else if(strncasecmp(thsign,"abs",3)==0) ithsign =  0;
  else if(strncasecmp(thsign,"pos",3)==0) ithsign = +1;
  else {
    printf("ERROR:  MRIbinarize01(): sign string = %s\n",thsign);
    return(NULL);
  }

  ncols   = vol->width;
  nrows   = vol->height;
  nslices = vol->depth;
  nframes = vol->nframes;

  if(binvol == NULL){
    binvol = MRIallocSequence(ncols, nrows, nslices, MRI_INT, nframes);
    if(binvol == NULL){
      printf("ERROR: MRIbinarize01(): could not alloc\n");
      return(NULL);
    }
  }
  else{
    if( (binvol->width != ncols) || (binvol->height != nrows) ||
  (binvol->depth != nslices) || (binvol->nframes != nframes) ){
      printf("ERROR: MRIbinarize01(): dimension missmatch\n");
      return(NULL);
    }
    if(binvol->type != MRI_INT){
      printf("ERROR: MRIbinarize01(): passed binvol "
       "type = %d, must be int (%d)\n",binvol->type,MRI_INT);
      return(NULL);
    }
  }

  nhits = 0;
  for(col = 0; col < vol->width; col++){
    for(row = 0; row < vol->height; row++){
      for(slice = 0; slice < vol->depth; slice++){
  for(frame = 0; frame < vol->nframes; frame++){

    switch ( vol->type ){
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

    if(ithsign ==  0) val = fabs(val);
    if(ithsign == -1) val = -val;

    r = 0;
    if(thmax > 0){
      if(val >= thmin && val <= thmax) r = 1;
    }
    else{
      if(val >= thmin ) r = 1;
    }

    if(invert) r = !r;

    if(r) MRIIseq_vox(binvol,col,row,slice,frame) = (short)highval;
    else  MRIIseq_vox(binvol,col,row,slice,frame) = (short)lowval;

    if(r) nhits ++;

  }
      }
    }
  }
  
  printf("INFO: MRIbinarize01(): nhits = %d\n",nhits);

  return(binvol);
}
       
