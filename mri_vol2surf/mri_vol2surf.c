/*----------------------------------------------------------
  Name: vol2surf.c
  $Id: mri_vol2surf.c,v 1.2 2001/02/06 23:27:22 greve Exp $
  Author: Douglas Greve
  Purpose: Resamples a volume onto a surface. The surface
  may be that of a subject other than the source subject.
  This replaces paint.c and is also used to convert functional
  data onto the icosohedron.

  Volume-to-Volume - V2V is a necessary step when converting functional
  data to talairach or painting onto the surface. The model as of 2/4/01
  assumes that the transformation is completely linear from one space to
  another, though there are variables for intermediate transformations
  built in. The four variables are: (1) Quantization matrix Q, (2) Field-
  of-view matrix F, (3) Warping matrix W, and (4), Registration matrix D.

  D - Registration matrix. Converts from Anatomical Space to Unwarpded
      Scanner Space.
  W - Warping matrix. Converts from Unwarpded Scanner Space to Warped
      Scanner Space.
  F - FOV matrix. Converts from Warpded Scanner Space to Field-of-View
      Space (the space created by the axes centered in the FOV and 
      parallel with the edges of the FOV).
  Q - Quantization matrix. Converts from FOV space to ColRowSlice Index
      Space.
  ----------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "icosahedron.h"

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

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_vol2surf.c,v 1.2 2001/02/06 23:27:22 greve Exp $";
char *Progname = NULL;

char *srcvolid   = NULL;
char *srcfmt     = NULL;
char *srcregfile = NULL;
char *srcwarp    = NULL;
int   srcoldreg  = 0;
char *srcsubject = NULL;

char *hemi    = NULL;
char *surfname = "white";
char *trgsubject = NULL;
int  IcoOrder = -1;
float IcoRadius = 100;
char *surfreg = "sphere.reg";
char *thicknessname = "thickness";
float ProjFrac = 0;

MRI_SURFACE *Surf    = NULL;
MRI_SURFACE *SurfOut = NULL;
MRI_SURFACE *SrcSurfReg = NULL;
MRI_SURFACE *TrgSurfReg = NULL;
int UseHash = 1;

char *outfile  = NULL;
char *ofmt     = NULL;
int  outtype;

char *srchitfile = NULL;
char *trghitfile = NULL;

char  *interpmethod_string = "nearest";
int    interpmethod = -1; 
char *mapmethod = "nnfr";

int debug = 0;

MATRIX *Dsrc, *Wsrc, *Fsrc, *Qsrc;
SXADAT *sxa;

char *SUBJECTS_DIR = NULL;
MRI *SrcVol, *SurfVals, *SurfVals2;
MRI *SrcHits, *SrcDist, *TrgHits, *TrgDist;

FILE *fp;

char tmpstr[2000];

int float2int_src;
char  *float2int_string = "round";
int    float2int = 0;
int ReverseMapFlag = 0;

int framesave = 0;

int main(int argc, char **argv)
{
  int n,err, f, vtx, svtx, tvtx;
  int nrows_src, ncols_src, nslcs_src, nfrms, endian, srctype;
  float ipr, bpr, intensity;
  float colres_src, rowres_src, slcres_src;
  float *framepower;
  char fname[2000];
  int nTrg121,nSrc121,nSrcLost;
  int nTrgMulti,nSrcMulti;
  float MnTrgMultiHits,MnSrcMultiHits;

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

  /* get info about the source volume */
  if(!strcmp(srcfmt,"bvolume")){
    err = bf_getvoldim(srcvolid,&nrows_src,&ncols_src,
           &nslcs_src,&nfrms,&endian,&srctype);
    if(err) exit(1);
    /* Dsrc: read the source registration file */
    err = regio_read_register(srcregfile, &srcsubject, &ipr, &bpr, 
            &intensity, &Dsrc, &float2int_src);
    if(err) exit(1);
    colres_src = ipr; /* in-plane resolution */
    rowres_src = ipr; /* in-plane resolution */
    slcres_src = bpr; /* between-plane resolution */
  }
  if(!strcmp(srcfmt,"cor")){
    nrows_src = 256;
    ncols_src = 256;
    nslcs_src = 256;
    endian = 0;
    nfrms = 1;
    srctype = 0;
    colres_src = 1;
    rowres_src = 1;
    slcres_src = 1;
  }
  /* Wsrc: Get the source warping Transform */
  Wsrc = NULL;

  /* Fsrc: Get the source FOV registration matrix */
  Fsrc = NULL;

  /* Qsrc: Compute the quantization matrix for src volume */
  Qsrc = FOVQuantMatrix(ncols_src,  nrows_src,  nslcs_src, 
      colres_src, rowres_src, slcres_src); 

  /* load in the (possibly 4-D) source volume */
  printf("Loading volume %s ...",srcvolid); fflush(stdout);
  if(!strcmp(srcfmt,"bvolume")){
    SrcVol = mri_load_bvolume(srcvolid);
    if(SrcVol == NULL) exit(1);
  }
  else{
    SrcVol = mri_load_cor_as_float(srcvolid);
  }
  printf("done\n");

  /* If this is a statistical volume, raise each frame to it's appropriate
     power (eg, stddev needs to be squared)*/
  if(is_sxa_volume(srcvolid)){
    printf("INFO: Source volume detected as selxavg format\n");
    sxa = ld_sxadat_from_stem(srcvolid);
    if(sxa == NULL) exit(1);
    framepower = sxa_framepower(sxa,&f);
    if(f != SrcVol->nframes){
      fprintf(stderr," number of frames is incorrect (%d,%d)\n",
        f,SrcVol->nframes);
      exit(1);
    }
    printf("INFO: Adjusting Frame Power\n");  fflush(stdout);
    mri_framepower(SrcVol,framepower);
  }

  /* Load the surface for subject indicated by registration file*/
  sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,surfname);
  printf("Reading surface %s\n",fname);
  Surf = MRISread(fname) ;
  if (!Surf)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;
  printf("Done\n");

  /* Load the thickness for projection along the normal*/
  if(ProjFrac > 0){
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,thicknessname);
    printf("Reading thickness %s\n",fname);
    MRISreadCurvatureFile(Surf, fname);
    printf("Done\n");
  }

  /* Map the values from the volume to the surface */
  printf("Mapping Source Volume onto Source Subject Surface\n");
  SurfVals = vol2surf_linear(SrcVol, Qsrc, Fsrc, Wsrc, Dsrc, 
           Surf, ProjFrac, interpmethod, 
           float2int);
  MRIfree(&SrcVol);

  if(trgsubject != NULL && strcmp(trgsubject,srcsubject)){
    /* load in the source subject registration */
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,surfreg);    
    SrcSurfReg = MRISread(fname);
    if(SrcSurfReg==NULL) exit(1);

    /* load in the target subject registration */
    if(strcmp(trgsubject,"ico")){
      sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,hemi,surfreg);    
      TrgSurfReg = MRISread(fname);
    }
    else TrgSurfReg = ReadIcoByOrder(IcoOrder,IcoRadius);
    if(TrgSurfReg==NULL) exit(1);

    if(!strcmp(mapmethod,"nnfr")) ReverseMapFlag = 1;
    else                          ReverseMapFlag = 0;

    printf("Mapping Surfaces (%s -> %s)\n",srcsubject,trgsubject);
    SurfVals2 = surf2surf_nnfr(SurfVals, SrcSurfReg,TrgSurfReg,
             &SrcHits,&SrcDist,&TrgHits,&TrgDist,
             ReverseMapFlag,UseHash);

    /*Compute some stats on mapping number of trgvtxs mapped from a source vtx*/
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
    MRISfree(&SrcSurfReg);

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

    if(!strcasecmp(ofmt,"w") || !strcasecmp(ofmt,"paint"))
      SurfOut = TrgSurfReg;
    else
      MRISfree(&TrgSurfReg);
    
    /* save the Source Hits into a .w file */
    if(srchitfile != NULL){
      for(vtx = 0; vtx < Surf->nvertices; vtx++)
  Surf->vertices[vtx].val = MRIFseq_vox(SrcHits,vtx,0,0,0) ;
      MRISwriteValues(Surf, srchitfile) ;
      MRIfree(&SrcHits);
    }
    /* save the Target Hits into a .w file */
    if(trghitfile != NULL){
      for(vtx = 0; vtx < SurfOut->nvertices; vtx++)
  SurfOut->vertices[vtx].val = MRIFseq_vox(TrgHits,vtx,0,0,0) ;
      MRISwriteValues(SurfOut, trghitfile) ;
      MRIfree(&TrgHits);
    }
  }
  else{
    SurfVals2 = SurfVals;
    SurfOut = Surf;
  }

  /* If this is a statistical volume, lower each frame to it's appropriate
     power (eg, variance needs to be sqrt'ed) */
  if(is_sxa_volume(srcvolid)){
    printf("INFO: Readjusting Frame Power\n");  fflush(stdout);
    for(f=0; f < SurfVals2->nframes; f++) framepower[f] = 1.0/framepower[f];
    mri_framepower(SurfVals2,framepower);
    sxa->nrows = 1;
    sxa->ncols = SurfVals2->width;
  }

  /* save the target volume in an appropriate format */
  if(!strcasecmp(ofmt,"bshort") || !strcasecmp(ofmt,"bfloat") || 
     !strcasecmp(ofmt,"bfile")  || !strcasecmp(ofmt,"bvolume") ){
    /*-------------- bvolume --------------*/
    if(!strcasecmp(ofmt,"bfile")  || !strcasecmp(ofmt,"bvolume") )
      outtype = srctype;
    else if(!strcasecmp(ofmt,"bshort")) outtype = BF_SHORT;
    else if(!strcasecmp(ofmt,"bfloat")) outtype = BF_FLOAT;

    printf("Saving volume to %s as bvolume\n",outfile); fflush(stdout);
    mri_save_as_bvolume(SurfVals2,outfile,endian,outtype); 

    /* for a stat volume, save the .dat file */
    if(is_sxa_volume(srcvolid)) sv_sxadat_by_stem(sxa,outfile);
  }
  if(!strcasecmp(ofmt,"w") || !strcasecmp(ofmt,"paint")){
    /*-------------- paint or .w --------------*/
    for(vtx = 0; vtx < SurfVals2->width; vtx++)
      SurfOut->vertices[vtx].val = MRIFseq_vox(SurfVals2,vtx,0,0,framesave) ;
    MRISwriteValues(SurfOut, outfile) ;
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

    /* -------- source volume inputs ------ */
    else if (!strcmp(option, "--srcvol")){
      if(nargc < 1) argnerr(option,1);
      srcvolid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcfmt")){
      if(nargc < 1) argnerr(option,1);
      srcfmt = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcreg")){
      if(nargc < 1) argnerr(option,1);
      srcregfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcoldreg")){
      srcoldreg = 1;
    }
    else if (!strcmp(option, "--srcwarp")){
      if(nargc < 1) argnerr(option,1);
      srcwarp = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--frame")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&framesave);
      nargsused = 1;
    }
    else if (!strcmp(option, "--surf")){
      if(nargc < 1) argnerr(option,1);
      surfname = pargv[0];
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
    else if (!strcmp(option, "--trgsubject")){
      if(nargc < 1) argnerr(option,1);
      trgsubject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--icoorder")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&IcoOrder);
      printf("IcoOrder = %d, nIcoVtxs = %d\n",IcoOrder,
       IcoNVtxsFromOrder(IcoOrder));
      nargsused = 1;
    }
    else if (!strcmp(option, "--surfreg")){
      if(nargc < 1) argnerr(option,1);
      surfreg = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--projfrac")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&ProjFrac);
      nargsused = 1;
    }
    else if (!strcmp(option, "--thickness")){
      if(nargc < 1) argnerr(option,1);
      thicknessname = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--interp")){
      if(nargc < 1) argnerr(option,1);
      interpmethod_string = pargv[0];
      interpmethod = interpolation_code(interpmethod_string);
      if(interpmethod == -1){
  fprintf(stderr,"ERROR: interpmethod = %s\n",interpmethod_string);
  fprintf(stderr,"  must be either nearest, tli, or sinc\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--float2int")){
      if(nargc < 1) argnerr(option,1);
      float2int_string = pargv[0];
      float2int = float2int_code(float2int_string);
      if(float2int == -1){
  fprintf(stderr,"ERROR: float2int = %s\n",float2int_string);
  fprintf(stderr,"  must be either round, floor, or tkreg\n");
  exit(1);
      }
      nargsused = 1;
    }
    else if (!strcmp(option, "--o")){
      if(nargc < 1) argnerr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srchits")){
      if(nargc < 1) argnerr(option,1);
      srchitfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trghits")){
      if(nargc < 1) argnerr(option,1);
      trghitfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--ofmt")){
      if(nargc < 1) argnerr(option,1);
      ofmt = pargv[0];
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
  fprintf(stderr, "   --srcvol    input volume path \n");
  fprintf(stderr, "   --srcfmt    input volume format \n");
  fprintf(stderr, "   --srcreg    source anat2scanner registration  \n");
  fprintf(stderr, "   --srcoldreg interpret srcreg as old-style reg.dat \n");
  fprintf(stderr, "   --srcwarp   source scanner warp table\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --surf       target surface (white) \n");
  fprintf(stderr, "   --hemi       hemisphere (lh or rh) \n");
  fprintf(stderr, "   --trgsubject target subject (if different than reg)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " Options for use with --trgsubject\n");
  fprintf(stderr, "   --surfreg    surface registration (sphere.reg)  \n");
  fprintf(stderr, "   --nohash flag to keep the hash table from being used. \n");
  fprintf(stderr, "   --icoorder   order of icosahedron when trgsubject=ico\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " Options for projecting along the surface normal:\n");
  fprintf(stderr, "   --projfrac  (0->1) projection along normal \n");  
  fprintf(stderr, "   --thickness thickness file (thickness)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " Options for output\n");
  fprintf(stderr, "   --o         output path\n");
  fprintf(stderr, "   --ofmt      output format\n");
  fprintf(stderr, "   --frame     save only nth frame (with --ofmt paint)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --interp    interpolation method (<nearest>, tli, or sinc)\n");
  fprintf(stderr, "   --float2int float-to-int conversion method "
    "(<round>, floor, or tkreg )\n");

  fprintf(stderr, "\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  fprintf(stderr, "\nThis program will resample one volume into another. \n") ;
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
  if(srcvolid == NULL){
    fprintf(stderr,"A source volume path must be supplied\n");
    exit(1);
  }
  if(interpmethod_string == NULL) interpmethod_string = "nearest";
  interpmethod = interpolation_code(interpmethod_string);
  if(interpmethod == -1){
    fprintf(stderr,"ERROR: interpmethod = %s\n",interpmethod_string);
    fprintf(stderr,"  must be either nearest, tli, or sinc\n");
    exit(1);
  }

  if(interpmethod != INTERP_NEAREST){
    fprintf(stderr,"ERROR: currently only nearest interpolation is supported\n");
    exit(1);
  }

  if(srcfmt == NULL) srcfmt = "bvolume";
  check_format(srcfmt);

  if(ofmt == NULL) ofmt = "bvolume";
  check_format(ofmt);

  if(hemi == NULL){
    fprintf(stderr,"ERROR: no hemifield specified\n");
    exit(1);
  }

  if(trgsubject != NULL && strcmp(trgsubject,"ico") && IcoOrder > -1){
    fprintf(stderr,"ERROR: --icoorder can only be used with '--trgsubject ico'\n");
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
      strcasecmp(trgfmt,"paint") != 0 &&
      strcasecmp(trgfmt,"cor") != 0 ){
    fprintf(stderr,"ERROR: format %s unrecoginized\n",trgfmt);
    fprintf(stderr,"Legal values are: bvolume, bfile, bshort, bfloat, and cor\n");
    fprintf(stderr,"                  paint, w\n");
    exit(1);
  }
  return(0);
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"srcvol = %s\n",srcvolid);
  if(srcfmt != NULL) fprintf(fp,"srcfmt = %s\n",srcfmt);
  else                  fprintf(fp,"srcfmt unspecified\n");
  if(srcregfile != NULL) fprintf(fp,"srcreg = %s\n",srcregfile);
  else                  fprintf(fp,"srcreg unspecified\n");
  fprintf(fp,"srcregold = %d\n",srcoldreg);
  if(srcwarp != NULL) fprintf(fp,"srcwarp = %s\n",srcwarp);
  else                   fprintf(fp,"srcwarp unspecified\n");

  fprintf(fp,"surf = %s\n",surfname);
  fprintf(fp,"hemi = %s\n",hemi);
  if(trgsubject != NULL){
    fprintf(fp,"trgsubject = %s\n",trgsubject);
    fprintf(fp,"surfreg = %s\n",surfreg);
  }
  if(ProjFrac > 0.0){
    fprintf(fp,"ProjFrac = %g\n",ProjFrac);
    fprintf(fp,"thickness = %s\n",thicknessname);
  }

  fprintf(fp,"interp = %s\n",interpmethod_string);
  fprintf(fp,"float2int = %s\n",float2int_string);

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


