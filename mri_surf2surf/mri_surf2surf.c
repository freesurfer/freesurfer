/*----------------------------------------------------------
  Name: mri_vol2surf.c
  $Id: mri_surf2surf.c,v 1.4 2001/09/25 23:09:04 greve Exp $
  Author: Douglas Greve
  Purpose: Resamples data from one surface onto another. If
  both the source and target subjects are the same, this is
  just a format conversion. The source or target subject may
  be ico.  Can handle data with multiple frames.
  -----------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "icosahedron.h"
#include "fio.h"

#include "MRIio.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mri.h"
#include "mri2.h"

#include "bfileio.h"
#include "registerio.h"
#include "resample.h"
#include "selxavgio.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
static int  check_format(char *fmt);
int GetNVtxsFromWFile(char *wfile);
int GetICOOrderFromValFile(char *filename, char *fmt);
int GetNVtxsFromValFile(char *filename, char *fmt);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_surf2surf.c,v 1.4 2001/09/25 23:09:04 greve Exp $";
char *Progname = NULL;

char *surfreg = "sphere.reg";
char *hemi    = NULL;

char *srcsubject = NULL;
char *srcvalfile = NULL;
char *srcfmt     = "bvolume";
MRI  *SrcVals, *SrcHits, *SrcDist;
MRI_SURFACE *SrcSurfReg;
char *SrcHitFile = NULL;
char *SrcDistFile = NULL;
int nSrcVtxs = 0;
int SrcIcoOrder;

char *trgsubject = NULL;
char *trgvalfile = NULL;
char *trgfmt     = "bvolume";
MRI  *TrgVals, *TrgValsSmth, *TrgHits, *TrgDist;
MRI_SURFACE *TrgSurfReg;
char *TrgHitFile = NULL;
char *TrgDistFile = NULL;
int TrgIcoOrder;

char *mapmethod = "nnfr";

int UseHash = 1;
int framesave = 0;
float IcoRadius = 100.0;
int nSmoothSteps = 0;
int nthstep, nnbrs, nthnbr, nbrvtx, frame;

int debug = 0;

char *SUBJECTS_DIR = NULL;
char *MRI_DIR = NULL;
SXADAT *sxa;
FILE *fp;

char tmpstr[2000];

int ReverseMapFlag = 0;
int cavtx = 0; /* command-line vertex -- for debugging */

int main(int argc, char **argv)
{
  int f,vtx,tvtx,svtx,n;
  float *framepower = NULL;
  char fname[2000];
  int nTrg121,nSrc121,nSrcLost;
  int nTrgMulti,nSrcMulti;
  float MnTrgMultiHits,MnSrcMultiHits, val;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR==NULL){
    fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  MRI_DIR = getenv("MRI_DIR") ;
  if(MRI_DIR==NULL){
    fprintf(stderr,"ERROR: MRI_DIR not defined in environment\n");
    exit(1);
  }

  /* --------- Load the registration surface for source subject --------- */
  if(!strcmp(srcsubject,"ico")){
    SrcIcoOrder = GetICOOrderFromValFile(srcvalfile,srcfmt);
    sprintf(fname,"%s/lib/bem/ic%d.tri",MRI_DIR,SrcIcoOrder);
    SrcSurfReg = ReadIcoByOrder(SrcIcoOrder, IcoRadius);
    printf("Source Ico Order = %d\n",SrcIcoOrder);
  }
  else{
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,surfreg);
    printf("Reading source surface reg %s\n",fname);
    SrcSurfReg = MRISread(fname) ;
  }
  if (!SrcSurfReg)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;
  printf("Done\n");

  /* ------------------ load the source data ----------------------------*/
  printf("Loading source data\n");
  if(!strcmp(srcfmt,"curv")){
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,srcvalfile);
    printf("Reading curvature file %s\n",fname);
    MRISreadCurvatureFile(SrcSurfReg, fname);
    SrcVals = MRIallocSequence(SrcSurfReg->nvertices, 1, 1,MRI_FLOAT,1);
    for(vtx = 0; vtx < SrcSurfReg->nvertices; vtx++){
      MRIFseq_vox(SrcVals,vtx,0,0,0) = SrcSurfReg->vertices[vtx].curv;
      if(vtx == cavtx){
  printf("vtx = %d, curv = %g, val = %g\n",vtx,
         SrcSurfReg->vertices[vtx].curv,
         MRIFseq_vox(SrcVals,vtx,0,0,0));
  DiagBreak();
      }
    }
  }
  else if(!strcmp(srcfmt,"paint") || !strcmp(srcfmt,"w")){
    MRISreadValues(SrcSurfReg,srcvalfile);
    SrcVals = MRIallocSequence(SrcSurfReg->nvertices, 1, 1,MRI_FLOAT,1);
    for(vtx = 0; vtx < SrcSurfReg->nvertices; vtx++)
      MRIFseq_vox(SrcVals,vtx,0,0,0) = SrcSurfReg->vertices[vtx].val;
  }
  else {
    SrcVals = mri_load_bvolume(srcvalfile); 
    if(SrcVals->width != SrcSurfReg->nvertices){
      fprintf(stderr,"ERROR: dimesion inconsitency in source data\n");
      fprintf(stderr,"       Number of surface vertices = %d\n",
        SrcSurfReg->nvertices);
      fprintf(stderr,"       Number of value vertices = %d\n",SrcVals->width);
      exit(1);
    }
    if(is_sxa_volume(srcvalfile)){
      printf("INFO: Source volume detected as selxavg format\n");
      sxa = ld_sxadat_from_stem(srcvalfile);
      if(sxa == NULL) exit(1);
      framepower = sxa_framepower(sxa,&f);
      if(f != SrcVals->nframes){
  fprintf(stderr," number of frames is incorrect (%d,%d)\n",
    f,SrcVals->nframes);
  exit(1);
      }
      printf("INFO: Adjusting Frame Power\n");  fflush(stdout);
      mri_framepower(SrcVals,framepower);
    }
  }
  if(SrcVals == NULL){
    fprintf(stderr,"ERROR loading source values from %s\n",srcvalfile);
    exit(1);
  }
  printf("Done\n");

  if(strcmp(srcsubject,trgsubject)){
    /* ------- Source and Target Subjects are different -------------- */
    /* ------- Load the registration surface for target subject ------- */
    if(!strcmp(trgsubject,"ico")){
      sprintf(fname,"%s/lib/bem/ic%d.tri",MRI_DIR,TrgIcoOrder);
      TrgSurfReg = ReadIcoByOrder(TrgIcoOrder, IcoRadius);
    }
    else{
      sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,hemi,surfreg);
      printf("Reading target surface reg %s\n",fname);
      TrgSurfReg = MRISread(fname) ;
    }
    if (!TrgSurfReg)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;
    printf("Done\n");
    
    if(!strcmp(mapmethod,"nnfr")) ReverseMapFlag = 1;
    else                          ReverseMapFlag = 0;
    
    /*-------------------------------------------------------------*/
    /* Map the values from the surface to surface */
    printf("Mapping Source Volume onto Source Subject Surface\n");
    TrgVals = surf2surf_nnfr(SrcVals, SrcSurfReg,TrgSurfReg,
           &SrcHits,&SrcDist,&TrgHits,&TrgDist,
           ReverseMapFlag,UseHash);
    
    
    /* Compute some stats on the mapping number of srcvtx mapping to a 
       target vtx*/
    nTrg121 = 0;
    MnTrgMultiHits = 0.0;
    for(tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++){
      n = MRIFseq_vox(TrgHits,tvtx,0,0,0);
      if(n == 1) nTrg121++;
      else MnTrgMultiHits += n;
    }
    nTrgMulti = TrgSurfReg->nvertices - nTrg121;
    if(nTrgMulti > 0) MnTrgMultiHits = (MnTrgMultiHits/nTrgMulti);
    else              MnTrgMultiHits = 0;
    printf("nTrg121 = %5d, nTrgMulti = %5d, MnTrgMultiHits = %g\n",
     nTrg121,nTrgMulti,MnTrgMultiHits);
    
    /* Compute some stats on the mapping number of trgvtxs mapped from a 
       source vtx*/
    nSrc121 = 0;
    nSrcLost = 0;
    MnSrcMultiHits = 0.0;
    for(svtx = 0; svtx < SrcSurfReg->nvertices; svtx++){
      n = MRIFseq_vox(SrcHits,svtx,0,0,0);
      if(n == 1) nSrc121++;
      else if(n == 0) nSrcLost++;
      else MnSrcMultiHits += n;
    }
    nSrcMulti = SrcSurfReg->nvertices - nSrc121;
    if(nSrcMulti > 0) MnSrcMultiHits = (MnSrcMultiHits/nSrcMulti);
    else              MnSrcMultiHits = 0;
    
    printf("nSrc121 = %5d, nSrcLost = %5d, nSrcMulti = %5d, MnSrcMultiHits = %g\n",
     nSrc121,nSrcLost,nSrcMulti,MnSrcMultiHits);
    
    /* save the Source Hits into a .w file */
    if(SrcHitFile != NULL){
      printf("INFO: saving source hits to %s\n",SrcHitFile);
      for(vtx = 0; vtx < SrcSurfReg->nvertices; vtx++)
  SrcSurfReg->vertices[vtx].val = MRIFseq_vox(SrcHits,vtx,0,0,0) ;
      MRISwriteValues(SrcSurfReg, SrcHitFile) ;
    }
    /* save the Source Distance into a .w file */
    if(SrcDistFile != NULL){
      printf("INFO: saving source distance to %s\n",SrcDistFile);
      for(vtx = 0; vtx < SrcSurfReg->nvertices; vtx++)
  SrcSurfReg->vertices[vtx].val = MRIFseq_vox(SrcDist,vtx,0,0,0) ;
      MRISwriteValues(SrcSurfReg, SrcDistFile) ;
    }
    /* save the Target Hits into a .w file */
    if(TrgHitFile != NULL){
      printf("INFO: saving target hits to %s\n",TrgHitFile);
      for(vtx = 0; vtx < TrgSurfReg->nvertices; vtx++)
  TrgSurfReg->vertices[vtx].val = MRIFseq_vox(TrgHits,vtx,0,0,0) ;
      MRISwriteValues(TrgSurfReg, TrgHitFile) ;
    }
    /* save the Target Hits into a .w file */
    if(TrgDistFile != NULL){
      printf("INFO: saving target distance to %s\n",TrgDistFile);
      for(vtx = 0; vtx < TrgSurfReg->nvertices; vtx++)
  TrgSurfReg->vertices[vtx].val = MRIFseq_vox(TrgDist,vtx,0,0,0) ;
      MRISwriteValues(TrgSurfReg, TrgDistFile) ;
    }
  }
  else{
    /* --- Source and Target Subjects are the same --- */
    printf("INFO: trgsubject = srcsubject\n");
    TrgSurfReg = SrcSurfReg;
    TrgVals = SrcVals;
  }
       
  /* Smooth if desired */
  if(nSmoothSteps > 0){
    nnbrs = TrgSurfReg->vertices[0].vnum;
    printf("INFO: Smoothing, NSteps = %d, NNbrs = %d\n",nSmoothSteps,nnbrs);

    TrgValsSmth = MRIcopy(TrgVals,NULL);

    for(nthstep = 0; nthstep < nSmoothSteps; nthstep ++){
      printf("Step = %d\n",nthstep);
      fflush(stdout);
      for(vtx = 0; vtx < TrgSurfReg->nvertices; vtx++){
  nnbrs = TrgSurfReg->vertices[vtx].vnum;
  //if(nnbrs != 5) printf("%4d %d\n",vtx,nnbrs);
  for(frame = 0; frame < TrgVals->nframes; frame ++){
    val = MRIFseq_vox(TrgVals,vtx,0,0,frame);
    for(nthnbr = 0; nthnbr < nnbrs; nthnbr++){
      nbrvtx = TrgSurfReg->vertices[vtx].v[nthnbr];
      val += MRIFseq_vox(TrgVals,nbrvtx,0,0,frame) ;
    }/* end loop over neighbor */
    MRIFseq_vox(TrgValsSmth,vtx,0,0,frame) = (val/(nnbrs+1));
  }/* end loop over frame */
      } /* end loop over vertex */
      MRIcopy(TrgValsSmth,TrgVals);
    }/* end loop over smooth step */
    MRIfree(&TrgValsSmth);
  }


  /* readjust frame power if necessary */
  if(is_sxa_volume(srcvalfile)){
    printf("INFO: Readjusting Frame Power\n");  fflush(stdout);
    for(f=0; f < TrgVals->nframes; f++) framepower[f] = 1.0/framepower[f];
    mri_framepower(TrgVals,framepower);
    sxa->nrows = 1;
    sxa->ncols = TrgVals->width;
  }

  /* ------------ save the target data -----------------------------*/
  printf("Saving target data\n");
  if(!strcmp(trgfmt,"paint") || !strcmp(trgfmt,"w")){
    for(vtx = 0; vtx < TrgSurfReg->nvertices; vtx++)
      TrgSurfReg->vertices[vtx].val = MRIFseq_vox(TrgVals,vtx,0,0,framesave);
    MRISwriteValues(TrgSurfReg,trgvalfile);
  }
  else {
    mri_save_as_bvolume(TrgVals,trgvalfile,0,BF_FLOAT); 
    if(is_sxa_volume(srcvalfile)) sv_sxadat_by_stem(sxa,trgvalfile);
  }

  return(0);
}
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
    else if (!strcasecmp(option, "--usehash")) UseHash = 1;
    else if (!strcasecmp(option, "--hash")) UseHash = 1;
    else if (!strcasecmp(option, "--dontusehash")) UseHash = 0;
    else if (!strcasecmp(option, "--nohash")) UseHash = 0;

    /* -------- source value inputs ------ */
    else if (!strcmp(option, "--srcsubject")){
      if(nargc < 1) argnerr(option,1);
      srcsubject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcsurfval")){
      if(nargc < 1) argnerr(option,1);
      srcvalfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcfmt")){
      if(nargc < 1) argnerr(option,1);
      srcfmt = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--nsmooth")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nSmoothSteps);
      if(nSmoothSteps < 1){
  fprintf(stderr,"ERROR: number of smooth steps (%d) must be >= 1\n",
    nSmoothSteps);
      }
      nargsused = 1;
    }

    /* -------- target value inputs ------ */
    else if (!strcmp(option, "--trgsubject")){
      if(nargc < 1) argnerr(option,1);
      trgsubject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trgicoorder")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&TrgIcoOrder);
      nargsused = 1;
    }
    else if (!strcmp(option, "--trgsurfval")){
      if(nargc < 1) argnerr(option,1);
      trgvalfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trgfmt")){
      if(nargc < 1) argnerr(option,1);
      trgfmt = pargv[0];
      if(!strcmp(trgfmt,"curv")){
  fprintf(stderr,"ERROR: Cannot select curv as target format\n");
  exit(1);
      }
      nargsused = 1;
    }

    else if (!strcmp(option, "--frame")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&framesave);
      nargsused = 1;
    }
    else if (!strcmp(option, "--cavtx")){
      /* command-line vertex -- for debugging */
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&cavtx);
      nargsused = 1;
    }
    else if (!strcmp(option, "--hemi")){
      if(nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--mapmethod")){
      if(nargc < 1) argnerr(option,1);
      mapmethod = pargv[0];
      if(strcmp(mapmethod,"nnfr") && strcmp(mapmethod,"nnf")){
  fprintf(stderr,"ERROR: mapmethod must be nnfr or nnf\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--srchits")){
      if(nargc < 1) argnerr(option,1);
      SrcHitFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcdist")){
      if(nargc < 1) argnerr(option,1);
      SrcDistFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trghits")){
      if(nargc < 1) argnerr(option,1);
      TrgHitFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trgdist")){
      if(nargc < 1) argnerr(option,1);
      TrgDistFile = pargv[0];
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
  fprintf(stderr, "USAGE: %s \n",Progname) ;
  fprintf(stderr, "\n");
  fprintf(stderr, "   --srcsubject source subject\n");
  fprintf(stderr, "   --srcsurfval path of file with input values \n");
  fprintf(stderr, "   --srcfmt     source format\n");
  fprintf(stderr, "   --trgsubject target subject\n");
  fprintf(stderr, "   --trgicoorder when trgsubject=ico\n");
  fprintf(stderr, "   --trgsurfval path of file in which to store output values\n");
  fprintf(stderr, "   --trgfmt     target format\n");
  fprintf(stderr, "   --surfreg    surface registration (sphere.reg)  \n");
  fprintf(stderr, "   --mapmethod  nnfr or nnf\n");
  fprintf(stderr, "   --hemi       hemisphere (lh or rh) \n");
  fprintf(stderr, "   --frame      save only nth frame (with --ofmt paint)\n");
  fprintf(stderr, "   --nsmooth    number of smoothing steps\n");  

  fprintf(stderr, "\n");
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"srcsubject = %s\n",srcsubject);
  fprintf(fp,"srcval = %s\n",srcvalfile);
  fprintf(fp,"srcfmt = %s\n",srcfmt);
  fprintf(fp,"trgsubject = %s\n",trgsubject);
  fprintf(fp,"trgval = %s\n",trgvalfile);
  fprintf(fp,"trgfmt = %s\n",trgfmt);
  fprintf(fp,"surfreg = %s\n",surfreg);
  fprintf(fp,"hemi = %s\n",hemi);
  fprintf(fp,"frame = %d\n",framesave);
  return;
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  fprintf(stderr, "\nThis program will resample one surface onto another. \n") ;
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
/* --------------------------------------------- */
static void check_options(void)
{
  if(srcsubject == NULL){
    fprintf(stderr,"ERROR: no source subject specified\n");
    exit(1);
  }
  if(srcvalfile == NULL){
    fprintf(stderr,"A source value path must be supplied\n");
    exit(1);
  }
  if(srcfmt == NULL) srcfmt = "bvolume";
  check_format(srcfmt);

  if(trgsubject == NULL){
    fprintf(stderr,"ERROR: no target subject specified\n");
    exit(1);
  }
  if(trgvalfile == NULL){
    fprintf(stderr,"A target value path must be supplied\n");
    exit(1);
  }
  if(trgfmt == NULL) trgfmt = "bvolume";
  check_format(trgfmt);

  if(hemi == NULL){
    fprintf(stderr,"ERROR: no hemifield specified\n");
    exit(1);
  }

  return;
}

/* --------------------------------------------- */
int check_format(char *trgfmt)
{
  if( strcasecmp(trgfmt,"bvolume") != 0 &&
      strcasecmp(trgfmt,"bfile") != 0 &&
      strcasecmp(trgfmt,"bshort") != 0 &&
      strcasecmp(trgfmt,"bfloat") != 0 &&
      strcasecmp(trgfmt,"w") != 0 &&
      strcasecmp(trgfmt,"curv") != 0 &&
      strcasecmp(trgfmt,"paint") != 0 ){
    fprintf(stderr,"ERROR: format %s unrecoginized\n",trgfmt);
    fprintf(stderr,"Legal values are: bvolume, bfile, bshort, bfloat, and \n");
    fprintf(stderr,"                  paint, w\n");
    exit(1);
  }
  return(0);
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
int GetNVtxsFromWFile(char *wfile)
{
  FILE *fp;
  int i,ilat, num, nvertices;
  int *vtxnum;
  float *wval;

  fp = fopen(wfile,"r");
  if (fp==NULL) {
    fprintf(stderr,"ERROR: Progname: GetNVtxsFromWFile():\n");
    fprintf(stderr,"Could not open %s\n",wfile);
    fprintf(stderr,"(%s,%d,%s)\n",__FILE__, __LINE__,__DATE__);
    exit(1);
  }
  
  fread2(&ilat,fp);
  fread3(&num,fp);
  vtxnum = (int *)   calloc(sizeof(int),   num);
  wval   = (float *) calloc(sizeof(float), num);

  for (i=0;i<num;i++){
    fread3(&vtxnum[i],fp);
    wval[i] = freadFloat(fp) ;
  }
  fclose(fp);

  nvertices = vtxnum[num-1] + 1;

  free(vtxnum);
  free(wval);

  return(nvertices);
}
/*---------------------------------------------------------------*/
int GetNVtxsFromValFile(char *filename, char *fmt)
{
  int nVtxs=0;
  int err,nrows, ncols, nslcs, nfrms, endian, type;

  if(!strcmp(fmt,"curv")){
    fprintf(stderr,"ERROR: cannot get nvertices from curv format\n");
    exit(1);
  }

  if(!strcmp(fmt,"paint") || !strcmp(fmt,"w"))
    nVtxs = GetNVtxsFromWFile(filename);

  if(!strcmp(fmt,"bvolume") || !strcmp(fmt,"bfloat") ||
     !strcmp(fmt,"bshort")  || !strcmp(fmt,"bfile")){
    err = bf_getvoldim(filename,&nrows,&ncols,
           &nslcs,&nfrms,&endian,&type);
    if(err) exit(1);
    if(nrows != 1 || nslcs != 1){
      fprintf(stderr,"ERROR: bvolume %s not a surface value file\n",filename);
      exit(1);
    }
    nVtxs = ncols;
  }

  return(nVtxs);
}
/*---------------------------------------------------------------*/
int GetICOOrderFromValFile(char *filename, char *fmt)
{
  int nIcoVtxs,IcoOrder;

  nIcoVtxs = GetNVtxsFromValFile(filename, fmt);

  IcoOrder = IcoOrderFromNVtxs(nIcoVtxs);
  if(IcoOrder < 0){
    fprintf(stderr,"ERROR: number of vertices = %d, does not mach ico\n",
      nIcoVtxs);
    exit(1);

  }
  
  return(IcoOrder);
}
