/*----------------------------------------------------------
  Name: mri_surf2surf.c
  $Id: mri_surf2surf.c,v 1.15 2004/02/17 00:50:26 greve Exp $
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

#include "mri.h"
#include "icosahedron.h"
#include "fio.h"

#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mri2.h"
#include "mri_identify.h"

#include "bfileio.h"
#include "registerio.h"
//  extern char *ResampleVtxMapFile;
#include "resample.h"
#include "selxavgio.h"
#include "prime.h"
#include "version.h"

MRI *MRISgaussianSmooth(MRIS *Surf, MRI *Src, double GStd, MRI *Targ);
MRI *MRISgaussianSmooth2(MRIS *Surf, MRI *Src, double GStd, MRI *Targ);
int MRISextendedNeighbors(MRIS *SphSurf,int TargVtxNo, int CurVtxNo,
			 double CosThresh, int *XNbrVtxNo, 
			 double *XNbrCos, int *nXNbrs,
			 int nXNbrsMax);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
int GetNVtxsFromWFile(char *wfile);
int GetICOOrderFromValFile(char *filename, char *fmt);
int GetNVtxsFromValFile(char *filename, char *fmt);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_surf2surf.c,v 1.15 2004/02/17 00:50:26 greve Exp $";
char *Progname = NULL;

char *surfreg = "sphere.reg";
char *hemi    = NULL;

char *srcsubject = NULL;
char *srcvalfile = NULL;
char *srctypestring = NULL;
int   srctype = MRI_VOLUME_TYPE_UNKNOWN;
MRI  *SrcVals, *SrcHits, *SrcDist;
MRI_SURFACE *SrcSurfReg;
char *SrcHitFile = NULL;
char *SrcDistFile = NULL;
int nSrcVtxs = 0;
int SrcIcoOrder;

char *trgsubject = NULL;
char *trgvalfile = NULL;
char *trgtypestring = NULL;
int   trgtype = MRI_VOLUME_TYPE_UNKNOWN;
MRI  *TrgVals, *TrgValsSmth, *TrgHits, *TrgDist;
MRI_SURFACE *TrgSurfReg;
char *TrgHitFile = NULL;
char *TrgDistFile = NULL;
int TrgIcoOrder;

MRI  *mritmp;
int  reshape = 1;
int  reshapefactor;

char *mapmethod = "nnfr";

int UseHash = 1;
int framesave = 0;
float IcoRadius = 100.0;
int nSmoothSteps = 0;
int nthstep, nnbrs, nthnbr, nbrvtx, frame;
double fwhm=0, gstd;

int debug = 0;

char *SUBJECTS_DIR = NULL;
char *FREESURFER_HOME = NULL;
SXADAT *sxa;
FILE *fp;

char tmpstr[2000];

int ReverseMapFlag = 0;
int cavtx = 0; /* command-line vertex -- for debugging */

int main(int argc, char **argv)
{
  int f,tvtx,svtx,n;
  float *framepower = NULL;
  char fname[2000];
  int nTrg121,nSrc121,nSrcLost;
  int nTrgMulti,nSrcMulti;
  float MnTrgMultiHits,MnSrcMultiHits;
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_surf2surf.c,v 1.15 2004/02/17 00:50:26 greve Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
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

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR==NULL){
    fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  FREESURFER_HOME = getenv("FREESURFER_HOME") ;
  if(FREESURFER_HOME==NULL){
    fprintf(stderr,"ERROR: FREESURFER_HOME not defined in environment\n");
    exit(1);
  }

  /* --------- Load the registration surface for source subject --------- */
  if(!strcmp(srcsubject,"ico")){ /* source is ico */
    SrcIcoOrder = GetICOOrderFromValFile(srcvalfile,srctypestring);
    sprintf(fname,"%s/lib/bem/ic%d.tri",FREESURFER_HOME,SrcIcoOrder);
    SrcSurfReg = ReadIcoByOrder(SrcIcoOrder, IcoRadius);
    printf("Source Ico Order = %d\n",SrcIcoOrder);
  }
  else{
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,surfreg);
    printf("Reading source surface reg %s\n",fname);
    SrcSurfReg = MRISread(fname) ;
    if(cavtx > 0) 
      printf("cavtx = %d, srcsurfreg: %g, %g, %g\n",cavtx,
	   SrcSurfReg->vertices[cavtx].x,
	   SrcSurfReg->vertices[cavtx].y,
	   SrcSurfReg->vertices[cavtx].z);
  }
  if (!SrcSurfReg)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;
  printf("Done\n");

  /* ------------------ load the source data ----------------------------*/
  printf("Loading source data\n");
  if(!strcmp(srctypestring,"curv")){ /* curvature file */
    if(fio_FileExistsReadable(srcvalfile)){
      memset(fname,0,strlen(fname));
      memcpy(fname,srcvalfile,strlen(srcvalfile));
    }
    else
      sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,srcvalfile);
    printf("Reading curvature file %s\n",fname);
    if(MRISreadCurvatureFile(SrcSurfReg, fname) != 0){
      printf("ERROR: reading curvature file\n");
      exit(1);
    }
    SrcVals = MRIcopyMRIS(NULL, SrcSurfReg, 0, "curv");
  }
  else if(!strcmp(srctypestring,"paint") || !strcmp(srctypestring,"w")){
    MRISreadValues(SrcSurfReg,srcvalfile);
    SrcVals = MRIcopyMRIS(NULL, SrcSurfReg, 0, "val");
  }
  else { /* Use MRIreadType */
    SrcVals =  MRIreadType(srcvalfile,srctype);
    if(SrcVals == NULL){
      printf("ERROR: could not read %s as type %d\n",srcvalfile,srctype);
      exit(1);
    }
    if(SrcVals->height != 1 || SrcVals->depth != 1){
      reshapefactor = SrcVals->height * SrcVals->depth;
      printf("Reshaping %d\n",reshapefactor);
      mritmp = mri_reshape(SrcVals, reshapefactor*SrcVals->width, 
         1, 1, SrcVals->nframes);
      MRIfree(&SrcVals);
      SrcVals = mritmp;
      reshapefactor = 0; /* reset for output */
    }

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
      sprintf(fname,"%s/lib/bem/ic%d.tri",FREESURFER_HOME,TrgIcoOrder);
      TrgSurfReg = ReadIcoByOrder(TrgIcoOrder, IcoRadius);
      reshapefactor = 6;
    }
    else{
      sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,hemi,surfreg);
      printf("Reading target surface reg %s\n",fname);
      TrgSurfReg = MRISread(fname) ;
    }
    if (!TrgSurfReg)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", 
    Progname, fname) ;
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
      if(n == 1)      nSrc121++;
      else if(n == 0) nSrcLost++;
      else MnSrcMultiHits += n;
    }
    nSrcMulti = SrcSurfReg->nvertices - nSrc121;
    if(nSrcMulti > 0) MnSrcMultiHits = (MnSrcMultiHits/nSrcMulti);
    else              MnSrcMultiHits = 0;
    
    printf("nSrc121 = %5d, nSrcLost = %5d, nSrcMulti = %5d, "
     "MnSrcMultiHits = %g\n", nSrc121,nSrcLost,nSrcMulti,
     MnSrcMultiHits);
    
    /* save the Source Hits into a .w file */
    if(SrcHitFile != NULL){
      printf("INFO: saving source hits to %s\n",SrcHitFile);
      MRIScopyMRI(SrcSurfReg, SrcHits, 0, "val");
      //for(vtx = 0; vtx < SrcSurfReg->nvertices; vtx++)
      //SrcSurfReg->vertices[vtx].val = MRIFseq_vox(SrcHits,vtx,0,0,0) ;
      MRISwriteValues(SrcSurfReg, SrcHitFile) ;
    }
    /* save the Source Distance into a .w file */
    if(SrcDistFile != NULL){
      printf("INFO: saving source distance to %s\n",SrcDistFile);
      MRIScopyMRI(SrcSurfReg, SrcDist, 0, "val");
      MRISwriteValues(SrcSurfReg, SrcDistFile) ;
      //for(vtx = 0; vtx < SrcSurfReg->nvertices; vtx++)
      //SrcSurfReg->vertices[vtx].val = MRIFseq_vox(SrcDist,vtx,0,0,0) ;
    }
    /* save the Target Hits into a .w file */
    if(TrgHitFile != NULL){
      printf("INFO: saving target hits to %s\n",TrgHitFile);
      MRIScopyMRI(TrgSurfReg, TrgHits, 0, "val");
      MRISwriteValues(TrgSurfReg, TrgHitFile) ;
      //for(vtx = 0; vtx < TrgSurfReg->nvertices; vtx++)
      //TrgSurfReg->vertices[vtx].val = MRIFseq_vox(TrgHits,vtx,0,0,0) ;
    }
    /* save the Target Hits into a .w file */
    if(TrgDistFile != NULL){
      printf("INFO: saving target distance to %s\n",TrgDistFile);
      MRIScopyMRI(TrgSurfReg, TrgDist, 0, "val");
      MRISwriteValues(TrgSurfReg, TrgDistFile) ;
      //for(vtx = 0; vtx < TrgSurfReg->nvertices; vtx++)
      //TrgSurfReg->vertices[vtx].val = MRIFseq_vox(TrgDist,vtx,0,0,0) ;
    }
  }
  else{
    /* --- Source and Target Subjects are the same --- */
    printf("INFO: trgsubject = srcsubject\n");
    TrgSurfReg = SrcSurfReg;
    TrgVals = SrcVals;
  }
       
  /* Smooth if desired */
  if(nSmoothSteps > 0)
    MRISsmoothMRI(TrgSurfReg, TrgVals, nSmoothSteps, TrgVals);
  if(fwhm > 0){
    printf("Gaussian smoothing with fwhm = %g, std = %g\n",fwhm,gstd);
    MRISgaussianSmooth(TrgSurfReg, TrgVals, gstd, TrgVals);
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
  if(!strcmp(trgtypestring,"paint") || !strcmp(trgtypestring,"w")){
    MRIScopyMRI(TrgSurfReg, TrgVals, framesave, "val");
    MRISwriteValues(TrgSurfReg,trgvalfile);
  }
  else {
    if(reshape){
      if(reshapefactor == 0) 
  reshapefactor = GetClosestPrimeFactor(TrgVals->width,6);
      
      printf("Reshaping %d (nvertices = %d)\n",reshapefactor,TrgVals->width);
      mritmp = mri_reshape(TrgVals, TrgVals->width / reshapefactor, 
         1, reshapefactor,TrgVals->nframes);
      if(mritmp == NULL){
  printf("ERROR: mri_reshape could not alloc\n");
  return(1);
      }
      MRIfree(&TrgVals);
      TrgVals = mritmp;
    }
    MRIwriteType(TrgVals,trgvalfile,trgtype);
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
    else if (!strcasecmp(option, "--noreshape")) reshape = 0;
    else if (!strcasecmp(option, "--reshape"))   reshape = 1;

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
    else if (!strcmp(option, "--srcfmt") ||
       !strcmp(option, "--src_type")){
      if(nargc < 1) argnerr(option,1);
      srctypestring = pargv[0];
      srctype = string_to_type(srctypestring);
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
    else if (!strcmp(option, "--fwhm")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&fwhm);
      gstd = fwhm/sqrt(log(256.0));
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
    else if (!strcmp(option, "--trgfmt") ||
       !strcmp(option, "--trg_type")){
      if(nargc < 1) argnerr(option,1);
      trgtypestring = pargv[0];
      if(!strcmp(trgtypestring,"curv")){
  fprintf(stderr,"ERROR: Cannot select curv as target format\n");
  exit(1);
      }
      trgtype = string_to_type(trgtypestring);
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
    else if (!strcmp(option, "--vtxmap")){
      if(nargc < 1) argnerr(option,1);
      ResampleVtxMapFile = pargv[0];
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
  fprintf(stderr, "   --src_type   source format\n");
  fprintf(stderr, "   --trgsubject target subject\n");
  fprintf(stderr, "   --trgicoorder when trgsubject=ico\n");
  fprintf(stderr, "   --trgsurfval path of file in which to store output values\n");
  fprintf(stderr, "   --trg_type   target format\n");
  fprintf(stderr, "   --hemi       hemisphere (lh or rh) \n");
  fprintf(stderr, "   --surfreg    surface registration (sphere.reg)  \n");
  fprintf(stderr, "   --mapmethod  nnfr or nnf\n");
  fprintf(stderr, "   --frame      save only nth frame (with --trg_type paint)\n");
  fprintf(stderr, "   --nsmooth    number of smoothing steps\n");  
  fprintf(stderr, "   --noreshape  do not reshape output to multiple 'slices'\n");  

  fprintf(stderr, "\n");
  printf("%s\n", vcid) ;
  printf("\n");

}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(

"This program will resample one surface onto another. The source and \n"
"target subjects can be any subject in $SUBJECTS_DIR and/or the  \n"
"icosahedron (ico). The source and target file formats can be anything \n"
"supported by mri_convert. The source format can also be a curvature \n"
"file or a paint (.w) file. The user also has the option of smoothing \n"
"on the surface. \n"
"\n"
"OPTIONS\n"
"\n"
"  --srcsubject subjectname\n"
"\n"
"    Name of source subject as found in $SUBJECTS_DIR or ico for icosahedron.\n"
"    The input data must have been sampled onto this subject's surface (eg, \n"
"    using mri_vol2surf)\n"
"\n"
"  --srcsurfval sourcefile\n"
"\n"
"    Name of file where the data on the source surface is located.\n"
"\n"
"  --src_type typestring\n"
"\n"
"    Format type string. Can be either curv (for FreeSurfer curvature file), \n"
"    paint or w (for FreeSurfer paint files), or anything accepted by \n"
"    mri_convert. If no type string  is given, then the type is determined \n"
"    from the sourcefile (if possible). If curv is used, then the curvature\n"
"    file will be looked for in $SUBJECTS_DIR/srcsubject/surf/hemi.sourcefile.\n"
"\n"
"  --trgsubject subjectname\n"
"\n"
"    Name of target subject as found in $SUBJECTS_DIR or ico for icosahedron.\n"
"\n"
"  --trgicoorder order\n"
"\n"
"    Icosahedron order number. This specifies the size of the\n"
"    icosahedron according to the following table: \n"
"              Order  Number of Vertices\n"
"                0              12 \n"
"                1              42 \n"
"                2             162 \n"
"                3             642 \n"
"                4            2562 \n"
"                5           10242 \n"
"                6           40962 \n"
"                7          163842 \n"
"    In general, it is best to use the largest size available.\n"
"\n"
"  --trgsurfval targetfile\n"
"\n"
"    Name of file where the data on the target surface will be stored.\n"
"    BUG ALERT: for trg_type w or paint, use the full path.\n"
"\n"
"  --trg_type typestring\n"
"\n"
"    Format type string. Can be paint or w (for FreeSurfer paint files) or anything\n"
"    accepted by mri_convert. NOTE: output cannot be stored in curv format\n"
"    If no type string  is given, then the type is determined from the sourcefile\n"
"    (if possible). If using paint or w, see also --frame.\n"
"\n"
"  --hemi hemifield (lh or rh)\n"
"\n"
"  --surfreg registration_surface"
"\n"
"    If the source and target subjects are not the same, this surface is used \n"
"    to register the two surfaces. sphere.reg is used as the default. Don't change\n"
"    this unless you know what you are doing.\n"
"\n"
"  --mapmethod methodname\n"
"\n"
"    Method used to map from the vertices in one subject to those of another.\n"
"    Legal values are: nnfr (neighest-neighbor, forward and reverse) and nnf\n"
"    (neighest-neighbor, forward only). Default is nnfr. The mapping is done\n"
"    in the following way. For each vertex on the target surface, the closest\n"
"    vertex in the source surface is found, based on the distance in the \n"
"    registration space (this is the forward map). If nnf is chosen, then the\n"
"    the value at the target vertex is set to that of the closest source vertex.\n"
"    This, however, can leave some source vertices unrepresented in target (ie,\n"
"    'holes'). If nnfr is chosen, then each hole is assigned to the closest\n"
"    target vertex. If a target vertex has multiple source vertices, then the\n"
"    source values are averaged together. It does not seem to make much difference. \n"
"\n"
"  --nsmooth niterations\n"
"\n"
"    Number of smoothing iterations. Each iteration consists of averaging each\n"
"    vertex with its neighbors. When only smoothing is desired, just set the \n"
"    the source and target subjects to the same subject.\n"
"\n"
"  --frame framenumber\n"
"\n"
"    When using paint/w output format, this specifies which frame to output. This\n"
"    format can store only one frame. The frame number is zero-based (default is 0).\n"
"\n"
"  --noreshape"
"\n"
"    By default, mri_surf2surf will save the output as multiple\n"
"    'slices'; has no effect for paint/w output format. For ico, the output\n"
"    will appear to be a 'volume' with Nv/R colums, 1 row, R slices and Nf \n"
"    frames, where Nv is the number of vertices on the surface. For icosahedrons, \n"
"    R=6. For others, R will be the prime factor of Nv closest to 6. Reshaping \n"
"    is for logistical purposes (eg, in the analyze format the size of a dimension \n"
"    cannot exceed 2^15). Use this flag to prevent this behavior. This has no \n"
"    effect when the output type is paint.\n"
"\n"
"EXAMPLES:\n"
"\n"
"1. Resample a subject's thickness of the left cortical hemisphere on to a \n"
"   7th order icosahedron and save in analyze4d format:\n"
"\n"
"   mri_surf2surf --hemi lh --srcsubject bert \n"
"      --srcsurfval thickness --src_type curv \n"
"      --trgsubject ico --trgicoorder 7 \n"
"      --trgsurfval bert-thickness-lh.img --trg_type analyze4d \n"
"\n"
"2. Resample data on the icosahedron to the right hemisphere of subject bert.\n"
"   Save in paint so that it can be viewed as an overlay in tksurfer. The \n"
"   source data is stored in bfloat format (ie, icodata_000.bfloat, ...)\n"
"\n"
"   mri_surf2surf --hemi rh --srcsubject ico \n"
"      --srcsurfval icodata-rh --src_type bfloat \n"
"      --trgsubject bert \n"
"      --trgsurfval ./bert-ico-rh.w --trg_type paint \n"
"\n"
"BUG REPORTS: send bugs to analysis-bugs@nmr.mgh.harvard.edu. Make sure \n"
"    to include the version and full command-line and enough information to\n"
"    be able to recreate the problem. Not that anyone does.\n"
"\n"
"\n"
"BUGS:\n"
"\n"
"  When the output format is paint, the output file must be specified with\n"
"  a partial path (eg, ./data-lh.w) or else the output will be written into\n"
"  the subject's anatomical directory.\n"
"\n"
"\n"
"AUTHOR: Douglas N. Greve, Ph.D., MGH-NMR Center (greve@nmr.mgh.harvard.edu)\n"
"\n");

  exit(1) ;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"srcsubject = %s\n",srcsubject);
  fprintf(fp,"srcval     = %s\n",srcvalfile);
  fprintf(fp,"srctype    = %s\n",srctypestring);
  fprintf(fp,"trgsubject = %s\n",trgsubject);
  fprintf(fp,"trgval     = %s\n",trgvalfile);
  fprintf(fp,"trgtype    = %s\n",trgtypestring);
  fprintf(fp,"surfreg    = %s\n",surfreg);
  fprintf(fp,"hemi       = %s\n",hemi);
  fprintf(fp,"frame      = %d\n",framesave);
  return;
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

  if(srctypestring == NULL){
    srctypestring = "bfloat";
    srctype = BFLOAT_FILE;
  }
  if( strcasecmp(srctypestring,"w") != 0 &&
      strcasecmp(srctypestring,"curv") != 0 &&
      strcasecmp(srctypestring,"paint") != 0 ){
    if(srctype == MRI_VOLUME_TYPE_UNKNOWN) {
  srctype = mri_identify(srcvalfile);
  if(srctype == MRI_VOLUME_TYPE_UNKNOWN){
    fprintf(stderr,"ERROR: could not determine type of %s\n",srcvalfile);
    exit(1);
  }
    }
  }

  if(trgsubject == NULL){
    fprintf(stderr,"ERROR: no target subject specified\n");
    exit(1);
  }
  if(trgvalfile == NULL){
    fprintf(stderr,"A target value path must be supplied\n");
    exit(1);
  }

  if(trgtypestring == NULL){
    trgtypestring = "bfloat";
    trgtype = BFLOAT_FILE;
  }
  if( strcasecmp(trgtypestring,"w") != 0 &&
      strcasecmp(trgtypestring,"curv") != 0 &&
      strcasecmp(trgtypestring,"paint") != 0 ){
    if(trgtype == MRI_VOLUME_TYPE_UNKNOWN) {
  trgtype = mri_identify(trgvalfile);
  if(trgtype == MRI_VOLUME_TYPE_UNKNOWN){
    fprintf(stderr,"ERROR: could not determine type of %s\n",trgvalfile);
    exit(1);
  }
    }
  }

  if(hemi == NULL){
    fprintf(stderr,"ERROR: no hemifield specified\n");
    exit(1);
  }

  if(fwhm != 0 && nSmoothSteps != 0){
    printf("ERROR: cannot specify --fwhm and --nsmooth\n");
    exit(1);
  }
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
//MRI *MRIreadHeader(char *fname, int type);
/*---------------------------------------------------------------*/
int GetNVtxsFromValFile(char *filename, char *typestring)
{
  //int err,nrows, ncols, nslcs, nfrms, endian;
  int nVtxs=0;
  int type;
  MRI *mri;

  printf("GetNVtxs: %s %s\n",filename,typestring);

  if(!strcmp(typestring,"curv")){
    fprintf(stderr,"ERROR: cannot get nvertices from curv format\n");
    exit(1);
  }

  if(!strcmp(typestring,"paint") || !strcmp(typestring,"w")){
    nVtxs = GetNVtxsFromWFile(filename);
    return(nVtxs);
  }

  type = string_to_type(typestring);
  mri = MRIreadHeader(filename, type);
  if(mri == NULL) exit(1);
  
  nVtxs = mri->width*mri->height*mri->depth;

  MRIfree(&mri);

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


/*-------------------------------------------------------------------
  MRISgaussianSmooth() - perform gaussian smoothing on a spherical 
  surface. The gaussian is defined by stddev GStd and is truncated
  at 4 stddevs. Note: this will change the val2bak of all the vertices.
  -------------------------------------------------------------------*/
MRI *MRISgaussianSmooth(MRIS *Surf, MRI *Src, double GStd, MRI *Targ)
{
  int vtxno1, vtxno2;
  float val;
  MRI *SrcTmp, *GSum, *GSum2;
  VERTEX *vtx1;
  double Radius, Radius2, dmax, GVar2, f, d, costheta, theta, g, dotprod;
  int n, err, nXNbrs, *XNbrVtxNo;
  double *XNbrDotProd, DotProdThresh;

  if(Surf->nvertices != Src->width){
    printf("ERROR: MRISgaussianSmooth: Surf/Src dimension mismatch\n");
    return(NULL);
  }

  if(Targ == NULL){
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth, 
            MRI_FLOAT, Src->nframes);
    if(Targ==NULL){
      printf("ERROR: MRISgaussianSmooth: could not alloc\n");
      return(NULL);
    }
  }
  else{
    if(Src->width   != Targ->width  || 
       Src->height  != Targ->height || 
       Src->depth   != Targ->depth  ||
       Src->nframes != Targ->nframes){
      printf("ERROR: MRISgaussianSmooth: output dimension mismatch\n");
      return(NULL);
    }
    if(Targ->type != MRI_FLOAT){
      printf("ERROR: MRISgaussianSmooth: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  /* Make a copy in case it's done in place */
  SrcTmp = MRIcopy(Src,NULL);

  /* This is for normalizing */
  GSum = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);
  if(GSum==NULL){
    printf("ERROR: MRISgaussianSmooth: could not alloc GSum\n");
    return(NULL);
  }

  GSum2 = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);
  if(GSum2==NULL){
    printf("ERROR: MRISgaussianSmooth: could not alloc GSum2\n");
    return(NULL);
  }

  vtx1 = &Surf->vertices[0] ;
  Radius2 = (vtx1->x * vtx1->x) + (vtx1->y * vtx1->y) + (vtx1->z * vtx1->z);
  Radius  = sqrt(Radius2);
  dmax = 4*GStd; // truncate after 4 stddevs
  GVar2 = 2*(GStd*GStd);
  f = 1/(sqrt(2*M_PI)*GStd);
  DotProdThresh = Radius2*cos(dmax/Radius)*(1.0001);

  printf("Radius = %g, gstd = %g, dmax = %g, GVar2 = %g, f = %g, dpt = %g\n",
	 Radius,gstd,dmax,GVar2,f,DotProdThresh);

  /* Initialize */
  for(vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++){
    MRIFseq_vox(GSum,vtxno1,0,0,0)  = 0;
    MRIFseq_vox(GSum2,vtxno1,0,0,0) = 0;
    for(frame = 0; frame < Targ->nframes; frame ++)
      MRIFseq_vox(Targ,vtxno1,0,0,frame) = 0;
    Surf->vertices[vtxno1].val2bak = -1;
  }
  
  XNbrVtxNo   = (int *) calloc(Surf->nvertices,sizeof(int));
  XNbrDotProd = (double *) calloc(Surf->nvertices,sizeof(double));

  if(0){
    // This will mess up future searches because it sets
    // val2bak to 0
    printf("Starting Search\n");
    err = MRISextendedNeighbors(Surf,0,0,DotProdThresh, XNbrVtxNo, 
			       XNbrDotProd, &nXNbrs, 1000);
    printf("Found %d (err=%d)\n",nXNbrs,err);
    for(n = 0; n < nXNbrs; n++){
      printf("%d %d %g\n",n,XNbrVtxNo[n],XNbrDotProd[n]);
    }
  }

  printf("nvertices = %d\n",Surf->nvertices);
  for(vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++){
    if(vtxno1%10000==0) 
    nXNbrs = 0;
    err = MRISextendedNeighbors(Surf,vtxno1,vtxno1,DotProdThresh, XNbrVtxNo, 
				XNbrDotProd, &nXNbrs, 1000);
    for(n = 0; n < nXNbrs; n++){
      vtxno2  = XNbrVtxNo[n];
      dotprod =  XNbrDotProd[n];
      costheta = dotprod/Radius2;

      // cos theta might be slightly > 1 due to precision
      if(costheta > +1.0) costheta = +1.0;
      if(costheta < -1.0) costheta = -1.0;

      // Compute the angle between the vertices
      theta = acos(costheta);

      /* Compute the distance bet vertices along the surface of the sphere */
      d = Radius * theta;

      /* Compute weighting factor for this distance */
      g = f*exp( -(d*d)/(GVar2) ); /* f not really nec */

      if(vtxno1 == 0 && 0){
	printf("%d %d %g %g %g %g %g\n",
	       vtxno1,vtxno2,dotprod,costheta,theta,d,g);
      }

      MRIFseq_vox(GSum,vtxno1,0,0,0)  += g;
      MRIFseq_vox(GSum2,vtxno1,0,0,0) += (g*g);
      
      for(frame = 0; frame < Targ->nframes; frame ++){
	val = g*MRIFseq_vox(SrcTmp,vtxno2,0,0,frame);
	MRIFseq_vox(Targ,vtxno1,0,0,frame) += val;
      }
      
    } /* end loop over vertex2 */
    
  } /* end loop over vertex1 */
  

  /* Normalize */
  for(vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++){
    vtx1 = &Surf->vertices[vtxno1] ;
    g = MRIFseq_vox(GSum,vtxno1,0,0,0);
    MRIFseq_vox(GSum2,vtxno1,0,0,0) /= (g*g);
    
    for(frame = 0; frame < Targ->nframes; frame ++){
      val = MRIFseq_vox(Targ,vtxno1,0,0,frame);
      MRIFseq_vox(Targ,vtxno1,0,0,frame) = val/g;
    }
  }

  MRIwrite(GSum,"gsum_000.bfloat");
  MRIwrite(GSum2,"gsum2_000.bfloat");

  MRIfree(&SrcTmp);
  MRIfree(&GSum);
  MRIfree(&GSum2);
  free(XNbrVtxNo);
  free(XNbrDotProd);

  return(Targ);
}

/*-------------------------------------------------------------------
  MRISextendedNeighbors() - read everything! Finds the set of
  "extended" neighbors of a given target vertex and a distance
  threshold. Call with CurVtxNo=TargVtxNo. There are other
  recursive calls where CurVtxNo!=TargVtxNothat. 

  The distance metric is the dot product between the vectors
  pointing from the origin to the target and current vertices. For 
  any given target vertex, the dot product must be greater than
  (NOT less than) the DotProdThresh to be included. 

  XNbrVtxNo is the (already allocated) list of vertex numbers of
  vertices that are within threshold.

  XNbrDotProd is the (already allocated) list of dot products of
  vertices that are within threshold.

  nXNbrs is a pointer to the total number of vertices that are 
  within threshold.

  nXNbrsMax is the maximum number allowed (ie, the allocation
  lengths of XNbrVtxNo and XNbrDotProd.

  NOTE: IMPORTANT!
  1. This really only works on the spherical surface.
  2. Assuming a spherical surface, the distance along the
     sphere between the two vertices can be computed as
        costheta = dotprod/(Radius^2);
        theta = acos(costheta);
        d = Radius * theta;
     This also allows you to work backwards to get a 
     dot product threshold from a given distance threshold.
  3. Modifies vertex->val2bak. It is important that this
     be set to -1 for all vertices prior to the first
     call or before calling with the same target vertex
     again.
  4. It is expected that this will be run for each vertex
     on a surface, so care has been taken to keep the 
     coputations light (eg, not allocating XNbrnVtxNo
     and XNbrnDotProd, using a dot product threshold
     instead of a distance threshold, not resetting val2bak).
  -------------------------------------------------------------------*/
int MRISextendedNeighbors(MRIS *SphSurf,int TargVtxNo, int CurVtxNo,
			  double DotProdThresh, int *XNbrVtxNo, 
			  double *XNbrDotProd, int *nXNbrs,
			  int nXNbrsMax)
{
  static int ncalls = 0;
  VERTEX *vtarg,*vcur;
  int nNNbrs, n, NbrVtxNo, err;
  double DotProd;

  // Get the current vertex
  vcur  = &SphSurf->vertices[CurVtxNo] ;

  // Return if this vertex has been hit
  if((int)vcur->val2bak == TargVtxNo) return(0);

  // Keep track of the number of recursive calls
  if(CurVtxNo == TargVtxNo){
    *nXNbrs = 0;
    ncalls = 0;
  }
  ncalls++;

  // Get the target vertex
  vtarg = &SphSurf->vertices[TargVtxNo] ;

  // Compute the dot product between the two
  DotProd = (vtarg->x*vcur->x) + (vtarg->y*vcur->y) + (vtarg->z*vcur->z);
  DotProd = fabs(DotProd);

  //printf("c %d %d %d %g %d\n",ncalls,TargVtxNo,CurVtxNo,DotProd,*nXNbrs);

  // Compare to threshold
  if(DotProd <= DotProdThresh) return(0);

  // Check whether another neigbor can be added
  if(*nXNbrs >= nXNbrsMax-1) return(1);

  // OK, add this vertex as an extended neighbor
  XNbrVtxNo[*nXNbrs] = CurVtxNo;
  XNbrDotProd[*nXNbrs]  = DotProd;
  (*nXNbrs)++;
  vcur->val2bak = TargVtxNo; // record a hit

  // Now, loop over the current nearest neighbors 
  nNNbrs = SphSurf->vertices[CurVtxNo].vnum;
  for(n = 0; n < nNNbrs; n++){
    NbrVtxNo = SphSurf->vertices[CurVtxNo].v[n];
    err = MRISextendedNeighbors(SphSurf, TargVtxNo, NbrVtxNo, DotProdThresh, 
				XNbrVtxNo, XNbrDotProd, nXNbrs, nXNbrsMax);
    if(err) return(err);
  }

  return(0);
}






/*-------------------------------------------------------------------*/
/* brute force approach */
/*-------------------------------------------------------------------*/
MRI *MRISgaussianSmooth2(MRIS *Surf, MRI *Src, double GStd, MRI *Targ)
{
  int vtxno1, vtxno2;
  float val;
  MRI *SrcTmp, *GSum, *GSum2;
  VERTEX *vtx1, *vtx2 ;
  double Radius, Radius2, dmax, GVar2, f, d;
  double costheta=0, theta=0, g, dotprod;
  int err, nXNbrs, XNbrVtxNo[1000];
  double XNbrDotProd[1000], DotProdThresh;

  if(Surf->nvertices != Src->width){
    printf("ERROR: MRISgaussianSmooth: Surf/Src dimension mismatch\n");
    return(NULL);
  }

  if(Targ == NULL){
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth, 
            MRI_FLOAT, Src->nframes);
    if(Targ==NULL){
      printf("ERROR: MRISgaussianSmooth: could not alloc\n");
      return(NULL);
    }
  }
  else{
    if(Src->width   != Targ->width  || 
       Src->height  != Targ->height || 
       Src->depth   != Targ->depth  ||
       Src->nframes != Targ->nframes){
      printf("ERROR: MRISgaussianSmooth: output dimension mismatch\n");
      return(NULL);
    }
    if(Targ->type != MRI_FLOAT){
      printf("ERROR: MRISgaussianSmooth: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }
  //Targ   = MRIcopy(Src,NULL);
  SrcTmp = MRIcopy(Src,NULL);

  /* This is for normalizing */
  GSum = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);
  if(GSum==NULL){
    printf("ERROR: MRISgaussianSmooth: could not alloc GSum\n");
    return(NULL);
  }

  GSum2 = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);
  if(GSum2==NULL){
    printf("ERROR: MRISgaussianSmooth: could not alloc GSum2\n");
    return(NULL);
  }

  vtx1 = &Surf->vertices[0] ;
  Radius2 = (vtx1->x * vtx1->x) + (vtx1->y * vtx1->y) + (vtx1->z * vtx1->z);
  Radius  = sqrt(Radius2);
  dmax = 4*GStd; // truncate after 4 stddevs
  GVar2 = 2*(GStd*GStd);
  f = 1/(sqrt(2*M_PI)*GStd);
  DotProdThresh = Radius2*cos(dmax/Radius);

  printf("Radius = %g, gstd = %g, dmax = %g, GVar2 = %g, f = %g, dpt = %g\n",
	 Radius,gstd,dmax,GVar2,f,DotProdThresh);

  for(vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++){
    MRIFseq_vox(GSum,vtxno1,0,0,0)  = 0;
    MRIFseq_vox(GSum2,vtxno1,0,0,0) = 0;
    for(frame = 0; frame < Targ->nframes; frame ++)
      MRIFseq_vox(Targ,vtxno1,0,0,frame) = 0;
    Surf->vertices[vtxno1].val2bak = -1;
  }
  
  printf("Starting Search\n");
  err = MRISextendedNeighbors(Surf,0,0,DotProdThresh, XNbrVtxNo, 
			     XNbrDotProd, &nXNbrs, 1000);
  printf("Found %d (err=%d)\n",nXNbrs,err);
  for(vtxno1 = 0; vtxno1 < nXNbrs; vtxno1++){
    printf("%d %d %g\n",vtxno1,XNbrVtxNo[vtxno1],XNbrDotProd[vtxno1]);
  }

  for(vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++){
    if(vtxno1%100==0) printf("vtxno1 = %d\n",vtxno1);
    vtx1 = &Surf->vertices[vtxno1] ;

    for(vtxno2 = vtxno1; vtxno2 < Surf->nvertices; vtxno2++){
      vtx2 = &Surf->vertices[vtxno2] ;

      if(vtxno1 != vtxno2){
	/* Compute the angle between the two vectors (dot product) */
	dotprod = (vtx1->x*vtx2->x) + (vtx1->y*vtx2->y) + (vtx1->z*vtx2->z);
	costheta = dotprod/Radius2;
	// cos theta might be slightly > 1 due to precision
	if(costheta > +1.0) costheta = +1.0;
	if(costheta < -1.0) costheta = -1.0;
	theta = acos(costheta);
	/* Compute the distance bet vertices along the surface of the sphere */
	d = Radius * theta;
	/* Compute weighting factor for this distance */
	g = f*exp( -(d*d)/(GVar2) );
      }
      else {d = 0; g = f; dotprod = Radius2;}

      /* Truncate */
      if(d > dmax) continue;

      if(vtxno1 == 0){
	printf("%d %d %g %g %g %g %g\n",
	       vtxno1,vtxno2,dotprod,costheta,theta,d,g);
      }

      MRIFseq_vox(GSum,vtxno1,0,0,0)  += g;
      MRIFseq_vox(GSum2,vtxno1,0,0,0) += (g*g);

      if(vtxno1 != vtxno2){
	MRIFseq_vox(GSum,vtxno2,0,0,0)  += g;
	MRIFseq_vox(GSum2,vtxno2,0,0,0) += (g*g);
      }
      
      for(frame = 0; frame < Targ->nframes; frame ++){
	val = g*MRIFseq_vox(SrcTmp,vtxno2,0,0,frame);
	MRIFseq_vox(Targ,vtxno1,0,0,frame) += val;
	if(vtxno1 != vtxno2){
	  val = g*MRIFseq_vox(SrcTmp,vtxno1,0,0,frame);
	  MRIFseq_vox(Targ,vtxno2,0,0,frame) += val;
	}
      }
      
    } /* end loop over vertex2 */
    
  } /* end loop over vertex1 */
  

  /* Normalize */
  for(vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++){
    vtx1 = &Surf->vertices[vtxno1] ;
    g = MRIFseq_vox(GSum,vtxno1,0,0,0);
    MRIFseq_vox(GSum2,vtxno1,0,0,0) /= (g*g);
    
    for(frame = 0; frame < Targ->nframes; frame ++){
      val = MRIFseq_vox(Targ,vtxno1,0,0,frame);
      MRIFseq_vox(Targ,vtxno1,0,0,frame) = val/g;
    }
  }

  MRIwrite(GSum,"gsum_000.bfloat");
  MRIwrite(GSum2,"gsum2_000.bfloat");

  MRIfree(&SrcTmp);
  MRIfree(&GSum);
  MRIfree(&GSum2);

  return(Targ);
}

