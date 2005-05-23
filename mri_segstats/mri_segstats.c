#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "version.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);

typedef struct {
  int id;
  char name[1000];
  int nhits;
  float vol;
  float min, max, range, mean, std;
} STATSUMENTRY;

int MRIsegFrameAvg(MRI *seg, int segid, MRI *mri, double *favg);
int *MRIsegIdList(MRI *seg, int *nlist, int frame);
int MRIsegCount(MRI *seg, int id, int frame);
int MRIsegStats(MRI *seg, int segid, MRI *mri,	int frame, 
		float *min, float *max, float *range, 
		float *mean, float *std);
int compare_ints(const void *v1,const void *v2);
int nunqiue_int_list(int *idlist, int nlist);
int *unqiue_int_list(int *idlist, int nlist, int *nunique);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_segstats.c,v 1.1 2005/05/23 23:42:41 greve Exp $";
char *Progname = NULL, *SUBJECTS_DIR = NULL;
char *SegVolFile = NULL;
char *InVolFile = NULL;
char *subject = NULL;
char *StatTableFile = NULL;
char *FrameAvgFile = NULL;
char *FrameAvgVolFile = NULL;
int DoFrameAvg = 0;
int frame = 0;
int synth = 0;
int debug = 0;
long seed = 0;
MRI *seg, *invol, *famri;
int nsegid, *segidlist;
int NonEmptyOnly = 0;

char *ctabfile = NULL;
COLOR_TABLE *ctab = NULL;
STATSUMENTRY *StatSumTable = NULL;

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs, n, nhits, f, nthseg, nsegidrep;
  float voxelvolume;
  float min, max, range, mean, std;
  FILE *fp;
  double  **favg;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if(subject != NULL){
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if(SUBJECTS_DIR==NULL){
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  /* Make sure we can open the output file*/
  fp = fopen(StatTableFile,"w");
  if(fp == NULL){
    printf("ERROR: could not open %s for writing\n",StatTableFile);
    exit(1);
  }
  fclose(fp);

  if(FrameAvgFile != NULL){
    fp = fopen(FrameAvgFile,"w");
    if(fp == NULL){
      printf("ERROR: could not open %s for writing\n",FrameAvgFile);
      exit(1);
    }
    fclose(fp);
  }

  printf("Loading %s\n",SegVolFile);
  seg = MRIread(SegVolFile);
  if(seg == NULL){
    printf("ERROR: loading %s\n",SegVolFile);
    exit(1);
  }

  if(InVolFile != NULL){
    printf("Loading %s\n",InVolFile);
    invol = MRIread(InVolFile);
    if(invol == NULL){
      printf("ERROR: loading %s\n",InVolFile);
      exit(1);
    }
    /* Should check that they are the same dim, etc*/
  }

  voxelvolume = seg->xsize * seg->ysize * seg->zsize;
  printf("Voxel Volume is %g mm^3\n",voxelvolume);

  if(ctabfile != NULL){
    ctab = CTABread(ctabfile);
    if(ctab == NULL){
      printf("ERROR: reading %s\n",ctabfile);
      exit(1);
    }
    nsegid = ctab->nbins;
    StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
    for(n=0; n < nsegid; n++){
      StatSumTable[n].id = ctab->bins[n].index;
      strcpy(StatSumTable[n].name, ctab->bins[n].name);
    }
  }
  else{
    /* Get list of segmentation ids */
    printf("Generating list of segmentation ids\n");
    segidlist = MRIsegIdList(seg, &nsegid,0);
    StatSumTable = (STATSUMENTRY *) calloc(sizeof(STATSUMENTRY),nsegid);
    for(n=0; n < nsegid; n++){
      StatSumTable[n].id = segidlist[n];
      strcpy(StatSumTable[n].name, "\0");
    }
  }
  printf("Found %3d segmentations\n",nsegid);

  printf("Computing statistics for each segmentation\n");
  for(n=0; n < nsegid; n++){
    printf("%3d ",n);
    if(n%20 == 19) printf("\n");
    fflush(stdout);
    nhits = MRIsegCount(seg, StatSumTable[n].id, 0);
    StatSumTable[n].nhits = nhits;
    if(InVolFile != NULL){
      MRIsegStats(seg, StatSumTable[n].id, invol, frame,
		  &min, &max, &range, &mean, &std);
      StatSumTable[n].min   = min;
      StatSumTable[n].max   = max;
      StatSumTable[n].range = range;
      StatSumTable[n].mean  = mean;
      StatSumTable[n].std   = std;
    }
  }
  printf("\n");

  if(debug){
    /* Dump the table to the screen */
    for(n=0; n < nsegid; n++){
      printf("%3d  %8d %10.1f  ", StatSumTable[n].id,StatSumTable[n].nhits,
	     voxelvolume*StatSumTable[n].nhits);
      if(ctabfile != NULL) printf("%-30s ",StatSumTable[n].name);
      if(InVolFile != NULL)
	printf("%10.4f %10.4f %10.4f %10.4f %10.4f ", 
	       StatSumTable[n].min, StatSumTable[n].max, 
	       StatSumTable[n].range, StatSumTable[n].mean, 
	       StatSumTable[n].std);
      printf("\n");
    }
  }

  // Count the number of actual rows to report
  nsegidrep = 0;
  for(n=0; n < nsegid; n++){
    if(NonEmptyOnly && voxelvolume*StatSumTable[n].nhits==0) continue;
    nsegidrep ++;
  }

  /* Print the table to the output file */
  if(StatTableFile != NULL){
    fp = fopen(StatTableFile,"w");
    fprintf(fp,"# Title Segmentation Statistics \n");
    fprintf(fp,"# SegVolFile %s \n",SegVolFile);
    if(ctabfile)  fprintf(fp,"# ColorTable %s \n",ctabfile);
    if(InVolFile) {
      fprintf(fp,"# InVolFile  %s \n",InVolFile);
      fprintf(fp,"# InVolFrame %d \n",frame);
    }
    if(NonEmptyOnly) fprintf(fp,"# Only reporting non-empty segmentations\n");
    fprintf(fp,"# Col 1 Index \n");
    fprintf(fp,"# Col 2 SegId \n");
    fprintf(fp,"# Col 3 NVoxels \n");
    fprintf(fp,"# Col 4 Volume \n");
    n = 5;
    if(ctabfile) {fprintf(fp,"# Col %d SegName \n",n); n++;}
    if(InVolFile) {
      fprintf(fp,"# Col %d Min \n",n);    n++;
      fprintf(fp,"# Col %d Max \n",n);    n++;
      fprintf(fp,"# Col %d Range \n",n);  n++;
      fprintf(fp,"# Col %d Mean \n",n);   n++;
      fprintf(fp,"# Col %d StdDev \n",n); n++;
    }
    fprintf(fp,"# NCols %d \n",n); 
    fprintf(fp,"# NRows %d \n",nsegidrep);  // Not right with nonempty
    nthseg = 1;
    for(n=0; n < nsegid; n++){
      if(NonEmptyOnly && voxelvolume*StatSumTable[n].nhits==0) continue;
      fprintf(fp,"%3d %3d  %8d %10.1f  ", nthseg, StatSumTable[n].id,
	      StatSumTable[n].nhits,
	      voxelvolume*StatSumTable[n].nhits);
      if(ctabfile != NULL) fprintf(fp,"%-30s ",StatSumTable[n].name);
      if(InVolFile != NULL)
	fprintf(fp,"%10.4f %10.4f %10.4f %10.4f %10.4f ", 
	       StatSumTable[n].min, StatSumTable[n].max, 
	       StatSumTable[n].range, StatSumTable[n].mean, 
	       StatSumTable[n].std);
      fprintf(fp,"\n");
      nthseg++;
    }
    fclose(fp);
  }

  if(DoFrameAvg){
    // Average input across space to create a waveform for each
    // Segmentation
    printf("Computing frame average\n");
    favg = (double **) calloc(sizeof(double *),nsegid);
    for(n=0; n < nsegid; n++)
      favg[n] = (double *) calloc(sizeof(double),invol->nframes);
    for(n=0; n < nsegid; n++){
      printf("%3d ",n);
      if(n%20 == 19) printf("\n");
      fflush(stdout);
      if(NonEmptyOnly && voxelvolume*StatSumTable[n].nhits == 0) continue;
      MRIsegFrameAvg(seg, StatSumTable[n].id, invol, favg[n]);
    }
    printf("\n");

    if(FrameAvgFile){
      // Save as a simple text file
      printf("Writing to %s\n",FrameAvgFile);
      fp = fopen(FrameAvgFile,"w");
      for(f=0; f < invol->nframes; f++){
	for(n=0; n < nsegid; n++){
	  if(NonEmptyOnly && voxelvolume*StatSumTable[n].nhits == 0) continue;
	  fprintf(fp,"%10.4lf ",favg[n][f]);
	}
	fprintf(fp,"\n");
      }
      fclose(fp);
    }

    if(FrameAvgVolFile){
      // Save as an MRI "volume"
      printf("Writing to %s\n",FrameAvgVolFile);
      famri = MRIallocSequence(nsegidrep,1,1,MRI_FLOAT,invol->nframes);
      for(f=0; f < invol->nframes; f++){
	nthseg = 0;
	for(n=0; n < nsegid; n++){
	  if(NonEmptyOnly && voxelvolume*StatSumTable[n].nhits == 0) continue;
	  MRIsetVoxVal(famri,nthseg,0,0,f,(float)favg[n][f]);
	  if(n==10 && f==67)
	    printf("MRI: seg 10, tp 67: %g  %g\n",favg[10][67],
		   MRIgetVoxVal(famri,nthseg,0,0,f));
	  nthseg++;
	}
	
      }
      MRIwrite(famri,FrameAvgVolFile);
    }

    printf("seg 10, tp 67: %g\n",favg[10][67]);
    famri = MRIread(FrameAvgVolFile);
    printf("MRI: seg 10, tp 67: %g  %g\n",favg[10][67],
	   MRIgetVoxVal(famri,10,0,0,67));

  }
  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
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
    else if (!strcasecmp(option, "--nonempty")) NonEmptyOnly = 1;
    else if ( !strcmp(option, "--seg") ) {
      if(nargc < 1) argnerr(option,1);
      SegVolFile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--in") ) {
      if(nargc < 1) argnerr(option,1);
      InVolFile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--sum") ) {
      if(nargc < 1) argnerr(option,1);
      StatTableFile = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--frameavg") ) {
      if(nargc < 1) argnerr(option,1);
      FrameAvgFile = pargv[0];
      DoFrameAvg = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--frameavgvol") ) {
      if(nargc < 1) argnerr(option,1);
      FrameAvgVolFile = pargv[0];
      DoFrameAvg = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--ctab") ) {
      if(nargc < 1) argnerr(option,1);
      ctabfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--frame")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&frame);
      nargsused = 1;
    }
    else if (!strcmp(option, "--subject")){
      if(nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--synth")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%ld",&seed);
      synth = 1;
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
  printf("   --seg segvol : segmentation volume path \n");
  printf("   --sum file   : stats summary table file \n");
  printf("\n");
  printf(" Other Options\n");
  printf("   --in invol : report more stats on the input volume\n");
  printf("   --ctab ctabfile : color table file with seg id names\n");
  printf("   --nonempty  only report non-empty segmentations\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;

  printf(
"This program will comute statistics. \n"
);
  exit(1) ;
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
/* --------------------------------------------- */
static void check_options(void)
{
  if(SegVolFile == NULL){
    printf("ERROR: must specify a segmentation volume\n");
    exit(1);
  }
  if(StatTableFile == NULL){
    printf("ERROR: must specify an output table file\n");
    exit(1);
  }

  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
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







/* ----------------------------------------------------------
   MRIsegCount() - returns the number of times the given 
   segmentation id appears in the volume.
   --------------------------------------------------------- */
int MRIsegCount(MRI *seg, int id, int frame)
{
  int nhits, v, c,r,s;
  nhits = 0;
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	v = (int) MRIgetVoxVal(seg,c,r,s,frame);
	if(v == id) nhits ++;
      }
    }
  }
  return(nhits);
}
/* ----------------------------------------------------------
   MRIsegIdList() - returns a list of the unique segmentation ids in
   the volume. The number in the list is *nlist. The volume need not
   be an int or char, but it is probably what it will be.
   --------------------------------------------------------- */
int *MRIsegIdList(MRI *seg, int *nlist, int frame)
{
  int nvoxels,r,c,s,nth;
  int *tmplist = NULL;
  int *segidlist = NULL;

  nvoxels = seg->width * seg->height * seg->depth;
  tmplist = (int *) calloc(sizeof(int),nvoxels);

  // First, load all voxels into a list
  nth = 0;
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	tmplist[nth] = (int) MRIgetVoxVal(seg,c,r,s,frame);
	nth++;
      }
    }
  }

  segidlist = unqiue_int_list(tmplist, nvoxels, nlist);
  free(tmplist);
  //for(nth=0; nth < *nlist; nth++)
  //printf("%3d %3d\n",nth,segidlist[nth]);
  return(segidlist);
}

/* --------------------------------------------- */
int compare_ints(const void *v1,const void *v2)
{
  int i1, i2;

  i1 = *((int*)v1);
  i2 = *((int*)v2);
  
  if(i1 < i2) return(-1);
  if(i1 > i2) return(+1);
  return(0);
}
/* --------------------------------------------------- 
   nunqiue_int_list() - counts the number of unique items
   in a list of integers. The list will be sorted.
   --------------------------------------------------- */
int nunqiue_int_list(int *idlist, int nlist)
{
  int idprev, nunique, n;

  qsort(idlist,nlist,sizeof(int),compare_ints);
  nunique = 1;
  idprev = idlist[0];
  for(n=1; n<nlist; n++){
    if(idprev != idlist[n]){
      nunique++;
      idprev = idlist[n];
    }
  }
  return(nunique);
}
/* --------------------------------------------------- 
   unqiue_int_list() - the returns the unique items
   in a list of integers. The list will be sorted.
   --------------------------------------------------- */
int *unqiue_int_list(int *idlist, int nlist, int *nunique)
{
  int n, *ulist, nthu;

  /* count number of unique elements in the list,
     this also sorts the list */
  *nunique = nunqiue_int_list(idlist, nlist);

  /* alloc the unqiue list */
  ulist = (int *) calloc(sizeof(int),*nunique);

  nthu = 0;
  ulist[nthu] = idlist[0];
  for(n=1; n<nlist; n++){
    if(ulist[nthu] != idlist[n]){
      nthu ++;
      ulist[nthu] = idlist[n];
    }
  }
  return(ulist);
}

/*---------------------------------------------------------
  MRIsegStats() - computes statistics within a given
  segmentation. Returns the number of voxels in the
  segmentation.
  ---------------------------------------------------------*/
int MRIsegStats(MRI *seg, int segid, MRI *mri,int frame, 
		float *min, float *max, float *range, 
		float *mean, float *std)
{
  int id,nvoxels,r,c,s;
  double val, sum, sum2;

  *min = 0;
  *max = 0;
  sum  = 0;
  sum2 = 0;
  nvoxels = 0;
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	id = (int) MRIgetVoxVal(seg,c,r,s,0);
	if(id != segid) continue;
	val =  MRIgetVoxVal(mri,c,r,s,frame);
	nvoxels++;
	if( nvoxels == 1 ){
	  *min = val;
	  *max = val;
	}
	if(*min > val) *min = val;
	if(*max < val) *max = val;
	sum  += val;
	sum2 += (val*val);
      }
    }
  }

  *range = *max - *min;

  if(nvoxels != 0){
    *mean = sum/nvoxels;
    *std = sqrt(((nvoxels)*(*mean)*(*mean) - 2*(*mean)*sum + sum2)/
		(nvoxels-1));
  }
  else {
    *mean = 0.0;
    *std = 0.0;
  }
  return(nvoxels);
}
/*---------------------------------------------------------
  MRIsegFrameAvg() - computes the average time course withing the
  given segmentation. Returns the number of voxels in the
  segmentation. favg must be preallocated to number of
  frames. favg = (double *) calloc(sizeof(double),mri->nframes);
  ---------------------------------------------------------*/
int MRIsegFrameAvg(MRI *seg, int segid, MRI *mri, double *favg)
{
  int id,nvoxels,r,c,s,f;
  double val;

  /* zero it out */
  for(f=0;f<mri->nframes;f++) favg[f] = 0;

  nvoxels = 0;
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	id = (int) MRIgetVoxVal(seg,c,r,s,0);
	if(id != segid) continue;
	for(f=0;f<mri->nframes;f++){
	  val =  MRIgetVoxVal(mri,c,r,s,f);
	  favg[f] += val;
	}
	nvoxels++;
      }
    }
  }

  if(nvoxels != 0)
    for(f=0;f<mri->nframes;f++) favg[f] /= nvoxels;

  return(nvoxels);
}

