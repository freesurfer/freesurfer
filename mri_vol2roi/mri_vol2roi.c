/*
  Name:    vol2roi.c
  Author:  Douglas N. Greve 
  email:   analysis-bugs@nmr.mgh.harvard.edu
  Date:    1/2/00
  Purpose: averages the voxels within an ROI. The ROI
           can be constrained structurally (with a label file)
           and/or functionally (with a volumetric mask)
  $Id: mri_vol2roi.c,v 1.23 2006/11/01 20:17:48 nicks Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "MRIio_old.h"
#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "label.h"

#include "bfileio.h"
#include "registerio.h"
#include "resample.h"
#include "corio.h"
#include "selxavgio.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "mri_circulars.h"

LABEL   *LabelReadFile(char *labelfile);

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
/*static int  isoptionflag(char *flag);*/

int CountLabelHits(MRI *SrcVol, MATRIX *Qsrc, MATRIX *Fsrc, 
       MATRIX *Wsrc, MATRIX *Dsrc, 
       MATRIX *Msrc2lbl, LABEL *Label, float labelfillthresh,
       int float2int);
int BTypeFromStem(char *stem);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_vol2roi.c,v 1.23 2006/11/01 20:17:48 nicks Exp $";
char *Progname = NULL;

char *roifile    = NULL;
char *roifmt     = "bvolume";
char *roitxtfile = NULL;
int  oldtxtstyle = 0;
int  plaintxtstyle = 0;

char *srcvolid   = NULL;
char *srcfmt     = NULL;
char *srcregfile = NULL;
char *srcwarp    = NULL;
int   srcoldreg  = 0;

char *labelfile   = NULL;
char *src2lblregfile = NULL;

char *mskvolid   = NULL;
char *mskfmt     = NULL;
char *mskregfile = NULL;
char *msk2srcregfile = NULL;
int   msksamesrc  = 1;

float  mskthresh = 0.5;
char  *msktail = "abs";
int    mskinvert = 0;
int    mskframe = 0;

char  *srcmskvolid = NULL;
char  *finalmskvolid = NULL;
char  *finalmskcrs = NULL;

LABEL *Label;

int debug = 0;

MATRIX *Dsrc, *Wsrc, *Fsrc, *Qsrc;
MATRIX *Dmsk, *Wmsk, *Fmsk, *Qmsk;
MATRIX *Msrc2lbl;
MATRIX *Mmsk2src;

SXADAT *sxa;

char *SUBJECTS_DIR = NULL;
char *FS_TALAIRACH_SUBJECT = NULL;
char *srcsubject, *msksubject;
char *regfile = "register.dat";
MRI *mSrcVol, *mROI, *mMskVol, *mSrcMskVol, *mFinalMskVol, *mSrcMskVol;
FILE *fp;
int nmskhits, nlabelhits, nfinalhits;

char tmpstr[2000];
int float2int_src, float2int_msk;
float labelfillthresh = .0000001;

int fixxfm = 0;
int labeltal  = 0;
char *talxfm = "talairach.xfm";

int bfiletype;

int main(int argc, char **argv)
{
  int n,err, f, nhits, r,c,s;
  int nrows_src, ncols_src, nslcs_src, nfrms, endian, srctype;
  int nrows_msk, ncols_msk, nslcs_msk, msktype;
  int roitype;
  float ipr, bpr, intensity;
  float colres_src, rowres_src, slcres_src;
  float colres_msk, rowres_msk, slcres_msk;
  float *framepower=NULL, val;
  LTA *lta;
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_vol2roi.c,v 1.23 2006/11/01 20:17:48 nicks Exp $", "$Name:  $");
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

  printf("--------------------------------------------------------\n");
  getcwd(tmpstr,2000);
  printf("%s\n",tmpstr);
  printf("%s\n",Progname);
  for(n=0;n<argc;n++) printf(" %s",argv[n]);
  printf("\n");
  printf("version %s\n",vcid);
  printf("--------------------------------------------------------\n");

  dump_options(stdout);

  /* ------------ get info about the source volume ------------------*/
  err = bf_getvoldim(srcvolid,&nrows_src,&ncols_src,
         &nslcs_src,&nfrms,&endian,&srctype);
  if(err) exit(1);
  /* Dsrc: read the source registration file */
  if(srcregfile != NULL){
    err = regio_read_register(srcregfile, &srcsubject, &ipr, &bpr, 
            &intensity, &Dsrc, &float2int_src);
    if(err) exit(1);
    colres_src = ipr; /* in-plane resolution */
    rowres_src = ipr; /* in-plane resolution */
    slcres_src = bpr; /* between-plane resolution */
    printf("srcreg Dsrc -------------\n");
    MatrixPrint(stdout,Dsrc);
    printf("----------------------------------\n");
  }
  else{
    Dsrc = NULL;
    colres_src = 1; /* in-plane resolution */
    rowres_src = 1; /* in-plane resolution */
    slcres_src = 1; /* between-plane resolution */
  }
  /* Wsrc: Get the source warping Transform */
  Wsrc = NULL;
  /* Fsrc: Get the source FOV registration matrix */
  Fsrc = NULL;
  /* Qsrc: Compute the quantization matrix for src volume */
  Qsrc = FOVQuantMatrix(ncols_src,  nrows_src,  nslcs_src, 
			colres_src, rowres_src, slcres_src); 
  printf("ras2vox src (tkreg) Qsrc -------------\n");
  MatrixPrint(stdout,Qsrc);
  printf("----------------------------------\n");

  /* ----------- load in the label ----------------- */
  if(labelfile != NULL){
    Label = LabelReadFile(labelfile);
    if(Label == NULL) exit(1);
    /* load in the source-to-label registration */
    if(src2lblregfile != NULL){
      //err = regio_read_xfm(src2lblregfile, &Msrc2lbl);
      //if(err) exit(1);
      lta = LTAread(src2lblregfile);
      if(lta->type == LINEAR_VOX_TO_VOX){
	printf("INFO: converting LTA to RAS\n");
	LTAvoxelTransformToCoronalRasTransform(lta);
      }
      Msrc2lbl = MatrixCopy(lta->xforms[0].m_L,NULL);
    }
    else if(labeltal){
      /* Load the talairach.xfm and make it approp for reg.dat*/
      Msrc2lbl = DevolveXFM(srcsubject,NULL,talxfm);
      if(Msrc2lbl==NULL) exit(1);
    }
    else Msrc2lbl = NULL;
    if(Msrc2lbl != NULL){
      printf("-- Source2Label %s ---- \n",src2lblregfile);
      MatrixPrint(stdout,Msrc2lbl);
      printf("-------------------------------\n");
    }
  }
  else{
    Label = NULL;
    Msrc2lbl = NULL;
  }

  /* -------------- load mask volume stuff -----------------------------*/
  if(mskvolid != NULL){
    /* get mask volume info */
    err = bf_getvoldim(mskvolid,&nrows_msk,&ncols_msk,
           &nslcs_msk,&nfrms,&endian,&msktype);
    if(err) exit(1);

    /* get the mask registration info */
    /* xyzFOV = Dmsk*xyzAnat (in mask space) */
    if(mskregfile != NULL){
      err = regio_read_register(mskregfile, &msksubject, &ipr, &bpr, 
        &intensity, &Dmsk, &float2int_msk);
      if(err) exit(1);
      colres_msk = ipr; /* in-plane resolution */
      rowres_msk = ipr; /* in-plane resolution */
      slcres_msk = bpr; /* between-plane resolution */
    }
    else if(msksamesrc){
      Dmsk = MatrixCopy(Dsrc,NULL);
      colres_msk = colres_src; /* in-plane resolution */
      rowres_msk = rowres_src; /* in-plane resolution */
      slcres_msk = slcres_src; /* between-plane resolution */
    }
    else{
      Dmsk = NULL;
      colres_msk = 1; /* in-plane resolution */
      rowres_msk = 1; /* in-plane resolution */
      slcres_msk = 1; /* between-plane resolution */
    }

    /* Qmsk: Compute the quantization matrix for msk volume */
    /* crsFOV = Qmsk*xyzFOV */
    Qmsk = FOVQuantMatrix(ncols_msk,  nrows_msk,  nslcs_msk, 
        colres_msk, rowres_msk, slcres_msk); 

    /* get the mask2source registration information */
    /* xyzSrc = Mmsk2src * xyzMsk */
    if(msk2srcregfile != NULL){
      err = regio_read_mincxfm(msk2srcregfile, &Mmsk2src, NULL);
      if(err) exit(1);
    }
    else Mmsk2src = NULL;

    /* load the mask volume (single frame) */
    mMskVol = mri_load_bvolume_frame(mskvolid, mskframe);
    if(mMskVol == NULL) exit(1);

    /* convert from Mask Anatomical to Src FOV */
    if(!msksamesrc){
      mSrcMskVol = vol2vol_linear(mMskVol, Qmsk, NULL, NULL, Dmsk, 
          Qsrc, Fsrc, Wsrc, Dsrc, 
          nrows_src, ncols_src, nslcs_src, 
          Mmsk2src, INTERP_NEAREST, float2int_msk);
      if(mSrcMskVol == NULL) exit(1);
    }
    else mSrcMskVol = mMskVol;

    /* binarize the mask volume */
    mri_binarize(mSrcMskVol, mskthresh, msktail, mskinvert,
     mSrcMskVol, &nmskhits);
  }
  else {mSrcMskVol = NULL; nmskhits = 0;}
  /*-------------- Done loading mask stuff -------------------------*/

  /* --------- load in the (possibly 4-D) source volume --------------*/
  printf("Loading volume %s ...",srcvolid); fflush(stdout);
  //mSrcVol = mri_load_bvolume(srcvolid);
  //mSrcVol->xsize = colres_src;
  //mSrcVol->ysize = rowres_src;
  //mSrcVol->zsize = slcres_src;
  bfiletype = BTypeFromStem(srcvolid);
  mSrcVol = MRIreadType(srcvolid,bfiletype);
  
  if(mSrcVol == NULL) exit(1);
  printf("done\n");
  /* If this is a statistical volume, raise each frame to it's appropriate
     power (eg, stddev needs to be squared)*/
  if(is_sxa_volume(srcvolid)){
    printf("INFO: Source volume detected as selxavg format\n");
    sxa = ld_sxadat_from_stem(srcvolid);
    if(sxa == NULL) exit(1);
    framepower = sxa_framepower(sxa,&f);
    if(f != mSrcVol->nframes){
      fprintf(stderr," number of frames is incorrect (%d,%d)\n",
        f,mSrcVol->nframes);
      exit(1);
    }
    printf("INFO: Adjusting Frame Power\n");  fflush(stdout);
    mri_framepower(mSrcVol,framepower);
  }

  /*--------- Prepare the final mask ------------------------*/
  if(Label != NULL){
    mFinalMskVol = label2mask_linear(mSrcVol, Qsrc, Fsrc, Wsrc, 
             Dsrc, mSrcMskVol,
             Msrc2lbl, Label, labelfillthresh, float2int_src, 
             &nlabelhits, &nfinalhits);

    if(mFinalMskVol == NULL) exit(1);
  }
  else {
    mFinalMskVol = mSrcMskVol;
    nfinalhits = nmskhits;
  }
  
  if(!oldtxtstyle){
    /* count the number of functional voxels = 1 in the mask */
    nfinalhits = 0;
    for(r=0;r<mFinalMskVol->height;r++){
      for(c=0;c<mFinalMskVol->width;c++){
	for(s=0;s<mFinalMskVol->depth;s++){
	  val = MRIFseq_vox(mFinalMskVol,c,r,s,0); 
          if(val > 0.5) nfinalhits ++;
	}
      }
    }    
    if(Label != NULL)
      nlabelhits = CountLabelHits(mSrcVol, Qsrc, Fsrc, 
			  Wsrc, Dsrc, Msrc2lbl, 
			  Label, labelfillthresh,float2int_src);
    else  nlabelhits = 0;
  }

  /*-------------------------------------------------------*/
  /*--------- Map the volume into the ROI -----------------*/
  printf("Averging over ROI\n");  fflush(stdout);
  mROI = vol2maskavg(mSrcVol, mFinalMskVol,&nhits);
  if(mROI == NULL) exit(1);
  printf("Done averging over ROI (nhits = %d)\n",nhits);
  /*-------------------------------------------------------*/

  /* ------- Save the final mask ------------------ */
  if(finalmskvolid != 0){
    //mri_save_as_bvolume(mFinalMskVol,finalmskvolid,endian,BF_FLOAT); 
    MRIwriteAnyFormat(mFinalMskVol,finalmskvolid,"bfloat",-1,NULL);
  }

  /* ------- Save CRS of the the final mask ------------------ */
  if(finalmskcrs != NULL){
    fp = fopen(finalmskcrs,"w");
    if(fp==NULL){
      fprintf(stderr,"ERROR: cannot open %s\n",finalmskcrs);
      exit(1);
    }
    for(r=0;r<mFinalMskVol->height;r++){
      for(c=0;c<mFinalMskVol->width;c++){
	for(s=0;s<mFinalMskVol->depth;s++){
	  val = MRIFseq_vox(mFinalMskVol,c,r,s,0); 
          if(val > 0.5){
	    fprintf(fp,"%d %d %d\n",c,r,s);
	  }
	}
      }
    }    
    fclose(fp);
  }

  /* If this is a statistical volume, lower each frame to it's appropriate
     power (eg, variance needs to be sqrt'ed) */
  if(is_sxa_volume(srcvolid)){
    printf("INFO: Readjusting Frame Power\n");  fflush(stdout);
    for(f=0; f < mROI->nframes; f++) framepower[f] = 1.0/framepower[f];
    mri_framepower(mROI,framepower);
  }

  /* save the target volume in an appropriate format */
  roitype = srctype;
  if(!strcasecmp(roifmt,"bshort") || !strcasecmp(roifmt,"bfloat") || 
     !strcasecmp(roifmt,"bfile")  || !strcasecmp(roifmt,"bvolume") ){
    /*-------------- bvolume --------------*/
    if(!strcasecmp(roifmt,"bfile")  || !strcasecmp(roifmt,"bvolume") )
      roitype = srctype;
    else if(!strcasecmp(roifmt,"bshort")) roitype = BF_SHORT;
    else if(!strcasecmp(roifmt,"bfloat")) roitype = BF_FLOAT;

    printf("Saving ROI to %s as bvolume\n",roifile); fflush(stdout);
    mri_save_as_bvolume(mROI,roifile,endian,roitype); 

    /* for a stat volume, save the .dat file */
    if(is_sxa_volume(srcvolid)) {
      sxa->nrows = 1;
      sxa->ncols = 1;
      sv_sxadat_by_stem(sxa,roifile);
    }
  }

  /* save as text */
  if(roitxtfile != NULL){
    fp = fopen(roitxtfile,"w");
    if(fp==NULL){
      fprintf(stderr,"ERROR: cannot open %s\n",roitxtfile);
      exit(1);
    }
    if(oldtxtstyle){
      printf("INFO: saving as old style txt\n");
      fprintf(fp,"%d \n",nmskhits);
    }
    if(! plaintxtstyle ){
      fprintf(fp,"%d \n",nlabelhits);
      fprintf(fp,"%d \n",nfinalhits);
    }
    for(f=0; f < mROI->nframes; f++) 
      fprintf(fp,"%9.4f\n",MRIFseq_vox(mROI,0,0,0,f));
    fclose(fp);
  }

  /* ------- Mask the source and save it  ------------------ */
  if(srcmskvolid != 0){
    for(r=0;r<mFinalMskVol->height;r++){
      for(c=0;c<mFinalMskVol->width;c++){
	for(s=0;s<mFinalMskVol->depth;s++){
	  val = MRIFseq_vox(mFinalMskVol,c,r,s,0); 
          if(val < 0.5){
	    for(f=0; f < mROI->nframes; f++) 
	      MRIFseq_vox(mSrcVol,c,r,s,f) = 0.0;
	  }
	}
      }
    }    
    MRIwrite(mSrcVol,srcmskvolid);
  }

  return(0);
}
/*--------------------------------------------------------------------*/
/* ------ ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^------------ */
/*--------------------------------------------------------------------*/

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

    else if (!strcasecmp(option, "--oldtxtstyle"))    oldtxtstyle = 1;
    else if (!strcasecmp(option, "--plaintxtstyle"))  plaintxtstyle = 1;
    else if (!strcasecmp(option, "--mskinvert"))  mskinvert = 1;
    else if (!strcmp(option, "--fixxfm"))   fixxfm = 1; 
    else if (!strcmp(option, "--nofixxfm")) fixxfm = 0; 

    /* -------- ROI output file ------ */
    else if (!strcmp(option, "--roiavgtxt")){
      if(nargc < 1) argnerr(option,1);
      roitxtfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--roiavg")){
      if(nargc < 1) argnerr(option,1);
      roifile = pargv[0];
      nargsused = 1;
    }

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

    /* -------- label inputs ------ */
    else if(!strcmp(option, "--labelfile") || 
      !strcmp(option, "--label")){
      if(nargc < 1) argnerr(option,1);
      labelfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--labelreg")){
      if(nargc < 1) argnerr(option,1);
      src2lblregfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--labeltal")){
      labeltal = 1;
      fixxfm = 1;
      nargsused = 0;
    }
    else if (!strcmp(option, "--talxfm")){
      if(nargc < 1) argnerr(option,1);
      talxfm = pargv[0];
      labeltal = 1;
      fixxfm = 1;
      nargsused = 1;
    }
    else if (!strcmp(option, "--labelfillthresh")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&labelfillthresh);
      nargsused = 1;
    }

    /* -------- mask volume inputs ------ */
    else if (!strcmp(option, "--mskvol")){
      if(nargc < 1) argnerr(option,1);
      mskvolid = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--mskfmt")){
      if(nargc < 1) argnerr(option,1);
      mskfmt = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--mskreg")){
      if(nargc < 1) argnerr(option,1);
      mskregfile = pargv[0];
      nargsused = 1;
      msksamesrc = 0;
    }
    else if (!strcmp(option, "--msksamesrc")){
      msksamesrc = 1;
      nargsused = 0;
    }
    else if (!strcmp(option, "--msktail")){
      if(nargc < 1) argnerr(option,1);
      msktail = pargv[0];
      nargsused = 1;
      if(strncasecmp(msktail,"abs",3) &&
	 strncasecmp(msktail,"pos",3) &&
	 strncasecmp(msktail,"neg",3)){
	fprintf(stderr,"ERROR: msk tail = %s, must be abs, pos, or neg\n",
		msktail);
	exit(1);
      }
    }

    else if (!strcmp(option, "--mskthresh")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&mskthresh);
      nargsused = 1;
    }

    else if (!strcmp(option, "--mskframe")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&mskframe);
      nargsused = 1;
    }

    else if (!strcmp(option, "--finalmskvol")){
      if(nargc < 1) argnerr(option,1);
      finalmskvolid = pargv[0];
      nargsused = 1;
    }

    else if (!strcmp(option, "--srcmskvol")){
      if(nargc < 1) argnerr(option,1);
      srcmskvolid = pargv[0];
      nargsused = 1;
    }

    else if (!strcmp(option, "--finalmskcrs")){
      if(nargc < 1) argnerr(option,1);
      finalmskcrs = pargv[0];
      nargsused = 1;
    }

    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
	fprintf(stderr,"       Did you really mean -%s ?\n",option);
      //printf("Match %d\n",strcmp(option, "--roiavgtxt"));
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
  printf("   --help : documentation and bug information\n");
  printf("\n");
  printf("   --srcvol    input volume path \n");
  //printf("   --srcfmt    input volume format \n");
  printf("   --srcreg    source registration (SrcXYZ = R*AnatXYZ) \n");
  //printf("   --srcoldreg interpret srcreg as old-style reg.dat \n");
  //printf("   --srcwarp   source scanner warp table\n");
  printf("\n");
  printf("   --label     path to label file \n");
  printf("   --labelreg  label registration (LabelXYZ = L*AnatXYZ) \n");
  printf("                  See warning before using this option.\n");
  printf("   --labeltal  use tal xfm for labelreg\n");
  printf("   --talxfm    use talxfm intead of talairach.xfm\n");
  printf("   --labelfillthresh thresh : fraction of voxel\n");
  printf("\n");
  printf("   --mskvol     mask volume path \n");
  //printf("   --mskfmt     mask volume format \n");
  //printf("   --mskreg     mask registration  (MaskXYZ = M*AnatXYZ)\n");
  printf("\n");
  printf("   --mskthresh threshold (0.5) mask threshold\n");
  printf("   --msktail   <abs>, pos, or neg (mask tail) \n");
  printf("   --mskframe  0-based mask frame <0> \n");
  printf("   --mskinvert : invert the mask \n");
  printf("\n");
  printf("   --roiavgtxt fname :output text file for ROI average\n");
  printf("   --roiavg    stem : output bfloat stem for ROI average\n");
  printf("   --finalmskvol path in which to save final mask\n");
  printf("   --finalmskcrs fname: save col,row,slice in text fname\n");
  printf("   --srcmskvol path in which to save masked source\n");
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf("\n%s\n\n", vcid) ;
  printf(

"This program will extract a region-of-interest (ROI) from a\n"
"volume. The ROI can be defined in one of three ways: (1) as set of\n"
"voxels in a mask volume the same size as the source volume, (2) as a\n"
"set of label points (defined in a label file), or (3) the intersection\n"
"of (1) and (2). mri_vol2roi can also be used to create a binary mask\n"
"volume from a label; see CREATING A BINARY MASK VOLUME FROM A LABEL\n"
"below.\n"
"\n"
"The result is a text file (argument of --roiavgtxt) with the following \n"
"rows of information. (1) The number of source voxels in the label (0\n"
"if no label is used). (2) The number of voxels in the final ROI. (3)\n"
"the average of the ROI voxels in the source volume. If the source\n"
"volume has multiple frames, then there will be a row for each frame.\n"
"\n"
"In addition, a final ROI mask volume can be saved. This volume is\n"
"the same size as the source volume (one frame only). The value of\n"
"a voxel is 1 if it was in the final ROI or 0 if it was out.\n"
"\n"
"--srcvol srcvolstem\n"
"\n"
"Specify the bfloat/bshort stem of the volume from which the ROI is to\n"
"be extracted.\n"
"\n"
"--srcreg regfile\n"
"\n"
"Registration between src volume and subject's anatomical (ie, a\n"
"register.dat). This is only needed when using a label.\n"
"\n"
"--label labelfile\n"
"\n"
"Path to FreeSurfer label file. Make sure to give full path. The label\n"
"file has, among other things, xyz coorindates of the points in the\n"
"label. Each point has an extent of 1mm^3. These xyz coorindates have\n"
"to be mapped into the source volume. This is done with the\n"
"registration matrix (--srcreg) and (possibly) with the label\n"
"registration (see --labelreg and --labeltal).\n"
"\n"
"--labelreg labelregfile\n"
"\n"
"File with matrix that maps label xyz coordinates to subject's\n"
"anatomical coordinates. Note: this matrix should map to the subject's\n"
"anatomical coordinates in 'tkregister' space. IMPORTANT: passing the\n"
"talairach.xfm file with this option can result in a MISINTERPRETATION\n"
"of the transform.\n"
"\n"
"--labeltal\n"
"\n"
"The label is in talairach space. In this case, the talairach.xfm file\n"
"for the subject in the register.dat (--srcreg) is loaded. This matrix\n"
"is modified to make it appropriate for use with the register.dat (ie,\n"
"it will be interpreted correctly). See also --talxfm and --labelreg.\n"
"\n"
"--talxfm xfmfile\n"
"\n"
"Use xfmfile in subject/mri/transforms for talairach transfrom\n"
"instead of talairach.xfm. This forces --labeltal.\n"
"\n"
"--labelfillthresh thresh\n"
"\n"
"Each label point represents 1mm^3, so it is possible that a source\n"
"voxel is not entirely filled by label points. This option allows the\n"
"user to choose the extent to which a voxel must be filled for it to be\n"
"considered part of the label. thresh is a number between 0 and 1 and\n"
"represents the fraction of the volume of a voxel that must be filled\n"
"in order for that voxel to be part of the label. Setting thresh=0\n"
"will force all voxels to be in the label regardless of what's in\n"
"the label file. Default is a small value just above 0.\n"
"\n"
"--mskvol maskstem\n"
"\n"
"Stem of the bshort/bfloat mask volume. This volume should be the\n"
"same dimension as the source volume (the number of frames can differ). \n"
"The mask volume will be thresholded to determine which voxels will\n"
"be included in the ROI. See --mskthresh, --msktail, --mskframe, and\n"
"--mskinvert for thresholding criteria.\n"
"\n"
"--mskthresh mskthresh \n"
"\n"
"Voxels value magnitude must be above mskthresh to be included. Default\n"
"is 0.5.\n"
"\n"
"--msktail tail\n"
"\n"
"Tail is the sign that a voxel must have in order to be included.\n"
"Tail can be abs, pos, or neg. Abs means that the sign should be\n"
"ignored (the default), pos means that the sign must be positive,\n"
"neg means that the sign must be negative.\n"
"\n"
"--mskframe nthframe\n"
"\n"
"Use the nthframe in the mask volume to derive the mask. nthframe is\n"
"0-based (ie, --mskframe 5 indicates the 6th frame). Default is 0\n"
"(ie, the first frame).\n"
"\n"
"--mskinvert \n"
"\n"
"After selecting voxels based on the threshold criteria (ie, mskthresh,\n"
"msktail, mskframe), invert the mask (ie, the non-selected voxels\n"
"become selected and vise versa).\n"
"\n"
"--roiavgtxt fname\n"
"\n"
"Save output in text format. See introduction.\n"
"\n"
"--roiavg stem\n"
"\n"
"Save output as a bfloat 'volume'. This flag is actually necessary\n"
"even if you are not going to use this output.\n"
"\n"
"--finalmskvol finalmaskstem\n"
"\n"
"Save the final set of voxels selected for the ROI in this volume. See\n"
"introduction for more info.\n"
"\n"
"--finalmskcrs fname\n"
"\n"
"Save the column, row, and slice of the voxels in the mask in a text\n"
"file. The indicies are zero-based.\n"
"\n"
"CREATING A BINARY MASK VOLUME FROM A LABEL\n"
"\n"
"mri_vol2roi can be used to create a binary mask volume, ie, a volume\n"
"with 1s where the label is and 0s everywhere else. Use something like: \n"
" \n"
"   mri_vol2roi --label your.label --srcvol f --srcreg register.dat \n"
"     --finalmskvol labelbinmask --roiavg /tmp/not.wanted.dat\n"
" \n"
"This will create a bshort volume with stem labelbinmask. This will\n"
"be the same size/geometry as the srcvol f. register.dat is the registration \n"
"matrix between the label space and srcvol space. mri_vol2roi requires \n"
"an output for the ROI information (--roiavg), but in this application \n"
"the user typically does not want this information (still has to be \n"
"there on the command line). \n"
" \n"
"BUGS\n"
"\n"
"The matrix used with --labelreg must map to the subject's\n"
"anatomical coordinates in 'tkregister' space. IMPORTANT: passing the\n"
"talairach.xfm file with this option can result in a MISINTERPRETATION\n"
"of the transform.\n"
"\n"
"The source, mask, and finalmask volumes must be in bshort/bfloat\n"
"format.\n"
"\n"
"A roiavgstem must be specified even if you do not want this format.\n"
"The roiavgtxt can be specified independently.\n"
) ;
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
  if(srcvolid == NULL){
    fprintf(stderr,"A source volume path must be supplied\n");
    exit(1);
  }
  if(roifile == NULL){
    fprintf(stderr,"A ROI output path must be supplied\n");
    exit(1);
  }

  if(msksamesrc && mskregfile != NULL){
    fprintf(stderr,"ERROR: cannot specify both --mskreg and --msksamesrc\n");
    exit(1);
  }

  if(msksamesrc) mskregfile = srcregfile;

  if(srcfmt == NULL) srcfmt = "bvolume";
  check_format(srcfmt);

  if(src2lblregfile != NULL && labeltal){
    printf("ERROR: cannot specify --labelreg and --labeltal\n");
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
      strcasecmp(trgfmt,"cor") != 0 ){
    fprintf(stderr,"ERROR: format %s unrecoginized\n",trgfmt);
    fprintf(stderr,"Legal values are: bvolume, bfile, bshort, bfloat, and cor\n");
    exit(1);
  }
  return(0);
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"roifile = %s\n",roifile);
  if(roitxtfile != NULL){
    fprintf(fp,"roitxtfile = %s\n",roitxtfile);
    fprintf(fp,"oldtxtstyle = %d\n",oldtxtstyle);
  }

  fprintf(fp,"srcvol = %s\n",srcvolid);
  if(srcfmt != NULL) fprintf(fp,"srcfmt = %s\n",srcfmt);
  else                  fprintf(fp,"srcfmt unspecified\n");
  if(srcregfile != NULL) fprintf(fp,"srcreg = %s\n",srcregfile);
  else                  fprintf(fp,"srcreg unspecified\n");
  fprintf(fp,"srcregold = %d\n",srcoldreg);
  if(srcwarp != NULL) fprintf(fp,"srcwarp = %s\n",srcwarp);
  else                   fprintf(fp,"srcwarp unspecified\n");

  if(labelfile != NULL){
    fprintf(fp,"label file = %s\n",labelfile);
    if(src2lblregfile != NULL) fprintf(fp,"labelreg = %s\n",src2lblregfile);
    else                     fprintf(fp,"labelreg unspecified\n");
    fprintf(fp,"label fill thresh = %g\n",labelfillthresh);
  }

  if(mskvolid != NULL) fprintf(fp,"mskvol = %s\n",mskvolid);
  else                  fprintf(fp,"mskvol unspecified\n");
  if(mskvolid != NULL) {
    if(mskfmt != NULL) fprintf(fp,"mskfmt = %s\n",mskfmt);
    else               fprintf(fp,"mskfmt unspecified\n");
    if(!msksamesrc){
      if(mskregfile != NULL) fprintf(fp,"mskreg = %s\n",mskregfile);
      else               fprintf(fp,"mskreg unspecified\n");
    }
    else fprintf(fp,"msk volume same as source\n");
    fprintf(fp,"msk tail = %s\n",msktail); 
    fprintf(fp,"msk threshold = %f\n",mskthresh); 
    fprintf(fp,"msk invert = %d\n",mskinvert); 
    fprintf(fp,"msk frame = %d\n",mskframe); 
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
#if 0
/*---------------------------------------------------------------*/
static int isoptionflag(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 3) return(0);

  if(flag[0] != '-' || flag[1] != '-') return(0);
  return(1);
}
#endif

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
LABEL   *LabelReadFile(char *labelfile)
{
  LABEL  *area ;
  char   *fname, line[STRLEN], *cp;
  FILE   *fp ;
  int    vno, nlines ;
  float  x, y, z, stat ;

  fname = labelfile;

  area = (LABEL *)calloc(1, sizeof(LABEL)) ;
  if (!area)
    ErrorExit(ERROR_NOMEMORY,"%s: could not allocate LABEL struct.",Progname);

  /* read in the file */
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NOFILE, "%s: could not open label file %s",
                       Progname, fname)) ;

  cp = fgetl(line, 199, fp) ;
  if (!cp)
    ErrorReturn(NULL,
                (ERROR_BADFILE, "%s: empty label file %s", Progname, fname)) ;
  if (!sscanf(cp, "%d", &area->n_points))
    ErrorReturn(NULL,
                (ERROR_BADFILE, "%s: could not scan # of lines from %s",
                 Progname, fname)) ;
  area->max_points = area->n_points ;
  area->lv = (LABEL_VERTEX *)calloc(area->n_points, sizeof(LABEL_VERTEX)) ;
  if (!area->lv)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, labelfile, sizeof(LV)*area->n_points) ;
  nlines = 0 ;
  while ((cp = fgetl(line, 199, fp)) != NULL)
  {
    if (sscanf(cp, "%d %f %f %f %f", &vno, &x, &y, &z, &stat) != 5)
      ErrorReturn(NULL, (ERROR_BADFILE, "%s: could not parse %dth line in %s",
                Progname, area->n_points, fname)) ;
    area->lv[nlines].x = x ;
    area->lv[nlines].y = y ;
    area->lv[nlines].z = z ;
    area->lv[nlines].stat = stat ;
    area->lv[nlines].vno = vno ;
    nlines++ ;
  }

  fclose(fp) ;
  if (!nlines)
    ErrorReturn(NULL,
                (ERROR_BADFILE, 
                 "%s: no data in label file %s", Progname, fname));
  return(area) ;
}

/*------------------------------------------------------------
  int CountLabelHits(): This constructs a mask only from the
  label as a way to count the number of functional voxels in
  the label itself.
  ------------------------------------------------------------*/
int CountLabelHits(MRI *SrcVol, MATRIX *Qsrc, MATRIX *Fsrc, 
       MATRIX *Wsrc, MATRIX *Dsrc, 
       MATRIX *Msrc2lbl, LABEL *Label, float labelfillthresh,
       int float2int)
{
  MRI * LabelMskVol;
  int nlabelhits, nfinalhits;
  int r,c,s;
  float val;

  LabelMskVol = label2mask_linear(mSrcVol, Qsrc, Fsrc, Wsrc, 
          Dsrc, NULL, Msrc2lbl,
          Label, labelfillthresh, float2int, 
          &nlabelhits, &nfinalhits);

  nlabelhits = 0;
  for(r=0;r<LabelMskVol->height;r++){
    for(c=0;c<LabelMskVol->width;c++){
      for(s=0;s<LabelMskVol->depth;s++){
  val = MRIFseq_vox(LabelMskVol,c,r,s,0); 
  if(val > 0.5) nlabelhits ++;
      }
    }
  }    
  MRIfree(&LabelMskVol);
  return(nlabelhits);
}

/*-------------------------------------------------------------
  BTypeFromStem() - determines whether stem is a bshort or bfloat.
  -------------------------------------------------------------*/
int BTypeFromStem(char *stem)
{
  char tmpstr[2000];

  sprintf(tmpstr,"%s_000.bfloat",stem);
  if(fio_FileExistsReadable(tmpstr)) return(BFLOAT_FILE);

  sprintf(tmpstr,"%s_000.bshort",stem);
  if(fio_FileExistsReadable(tmpstr)) return(BSHORT_FILE);

  printf("WARNING: cannot find file for bstem %s\n",stem);

  return(MRI_VOLUME_TYPE_UNKNOWN);
}
