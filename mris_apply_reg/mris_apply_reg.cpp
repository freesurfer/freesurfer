/**
 * @brief Applies multiple surface registrations
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



/*!
\file mris_apply_reg.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/



/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "resample.h"
#include "pdf.h"
#include "icosahedron.h"
#include "mrisutils.h"
#include "mri2.h"
#include "gcamorph.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
void usage_message(FILE *stream);
void usage(FILE *stream);
int main(int argc, char *argv[]) ;

const char *Progname = NULL;
static int center_surface=0 ;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *SrcValFile=NULL;
char *TrgValFile=NULL;
char *SurfRegFile[100];
char *SurfPatchFile[100];
int ReverseMapFlag = 1;
int DoJac = 0;
int UseHash = 1;
int nsurfs = 0;
int npatches = 0;
int DoSynthRand = 0;
int DoSynthOnes = 0;
int SynthSeed = -1;
char *AnnotFile = NULL;
char *LabelFile = NULL;
char *SurfXYZFile = NULL;
int OutputCurvFormat=0;
LABEL *MRISmask2Label(MRIS *surf, MRI *mask, int frame, double thresh);
int ApplyScaleSurf(MRIS *surf, const double scale);
double SourceSurfRegScale = 0;
double TargetSurfRegScale = 0;
int RegScannerRAS = 0;
int XYZScannerRAS = 0; // convert the source surf xyz to scanner ras

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,n,err;
  MRIS *SurfReg[100],*SurfSrc;
  MRI *SrcVal, *TrgVal;
  char *base;
  COLOR_TABLE *ctab=NULL;

  nargs = handleVersionOption(argc, argv, "mris_apply_reg");
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

  // Load in surface registrations
  if(RegScannerRAS) printf("Converting reg surf vertices to scanner RAS\n");
  for(n=0; n<nsurfs;n++){
    printf("%d Loading %s\n",n+1,SurfRegFile[n]);
    base = fio_basename(SurfRegFile[n],".tri");
    if(strcmp(base,"ic7")==0){
      // Have to do it this way to rescale. Need to find a better more robust way.
      printf("   reading as ico 7, rescaling radius to 100\n");
      SurfReg[n] = ReadIcoByOrder(7, 100);
    }
    else {
      SurfReg[n] = MRISread(SurfRegFile[n]);
    }
    if(! SurfReg[n]) exit(1);
    MRISsaveVertexPositions(SurfReg[n], TMP_VERTICES); // In case --vertex-pair is set
    free(base);
    if(SurfReg[n]==NULL) exit(1);
    if(npatches>0){
      printf("Loading patch %s\n",SurfPatchFile[n]);
      MRISreadPatch(SurfReg[n], SurfPatchFile[n]);
    }
    if(RegScannerRAS) MRIStkr2Scanner(SurfReg[n]);
  }

  if(SourceSurfRegScale > 0){
    printf("Scaling first surface by %g\n",SourceSurfRegScale);
    ApplyScaleSurf(SurfReg[0], SourceSurfRegScale);
  }
  if(TargetSurfRegScale > 0){
    printf("Scaling second surface by %g\n",TargetSurfRegScale);
    ApplyScaleSurf(SurfReg[1], TargetSurfRegScale);
  }

  // Load in source data
  SrcVal = NULL;
  if(DoSynthRand) {
    if (SynthSeed < 0) SynthSeed = PDFtodSeed();
    printf("INFO: synthesizing, seed = %d\n",SynthSeed);
    srand48(SynthSeed);
    MRIrandn(SrcVal->width, SrcVal->height, SrcVal->depth,
             SrcVal->nframes,0, 1, SrcVal);
  }
  else if(DoSynthOnes != 0) {
    printf("INFO: filling input with all 1s\n");
    MRIconst(SrcVal->width, SrcVal->height, SrcVal->depth,
             SrcVal->nframes, 1, SrcVal);
  }
  else if(AnnotFile) {
    printf("Loading annotation %s\n",AnnotFile);
    err = MRISreadAnnotation(SurfReg[0], AnnotFile);
    if(err) exit(1);
    SrcVal = MRISannotIndex2Seg(SurfReg[0]);
    ctab = CTABdeepCopy(SurfReg[0]->ct);
  }
  else if(SurfXYZFile) {
    printf("Loading surface xyz %s\n",SurfXYZFile);
    SurfSrc = MRISread(SurfXYZFile);
    if(SurfSrc==NULL)  exit(1);
    if(XYZScannerRAS) {
      printf("Converting input surface to Scanner RAS\n");
      MRIStkr2Scanner(SurfSrc);
    }
    SrcVal = MRIcopyMRIS(NULL, SurfSrc, 2, "z"); // start at z to autoalloc
    MRIcopyMRIS(SrcVal, SurfSrc, 0, "x");
    MRIcopyMRIS(SrcVal, SurfSrc, 1, "y");
  }
  else if(LabelFile) {
    LABEL *srclabel;
    printf("Loading label %s\n",LabelFile);
    srclabel = LabelRead(NULL, LabelFile);
    if(srclabel == NULL) exit(1);
    SrcVal = MRISlabel2Mask(SurfReg[0],srclabel,NULL);
    printf("   %d points in input label\n",srclabel->n_points);
    LabelFree(&srclabel);
  }
  else {
    printf("Loading %s\n",SrcValFile);
    SrcVal = MRIread(SrcValFile);
    if(SrcVal==NULL) exit(1);
    MRI *tmpmri;
    tmpmri = MRIreshape1d(SrcVal,NULL);
    if(tmpmri == NULL) exit(1);
    MRIfree(&SrcVal);
    SrcVal = tmpmri;
    if(SrcVal->type != MRI_FLOAT) {
      printf("Converting source to float\n");
      tmpmri = MRISeqchangeType(SrcVal,MRI_FLOAT,0,0,0);
      if (tmpmri == NULL) {
        printf("ERROR: could change type\n");
        exit(1);
      }
      MRIfree(&SrcVal);
      SrcVal = tmpmri;
    }
  }

  // Apply registration to source
  TrgVal = MRISapplyReg(SrcVal, SurfReg, nsurfs, ReverseMapFlag, DoJac, UseHash);
  if(TrgVal == NULL) exit(1);

  // Save output
  if(AnnotFile){
    printf("Converting to target annot\n");
    err = MRISseg2annot(SurfReg[nsurfs-1],TrgVal,ctab);
    if(err) exit(1);
    printf("Writing %s\n",TrgValFile);
    MRISwriteAnnotation(SurfReg[nsurfs-1], TrgValFile);
  } 
  else if(LabelFile){
    LABEL *label;
    label = MRISmask2Label(SurfReg[nsurfs-1], TrgVal, 0, 10e-5);
    printf("   %d points in output label\n",label->n_points);
    err = LabelWrite(label,TrgValFile);
    if(err){
      printf("ERROR: writing label file %s\n",TrgValFile);
      exit(1);
    }
    LabelFree(&label);
  }
  else if(SurfXYZFile){
    printf("Writing surface to %s\n",TrgValFile);
    // If this is a patch, reload the target surf without the patch
    MRIS *tmpsurf = SurfReg[nsurfs-1];
    if(npatches > 0) tmpsurf = MRISread(SurfRegFile[nsurfs-1]);
    MRIScopyMRI(tmpsurf,TrgVal,0,"x");
    MRIScopyMRI(tmpsurf,TrgVal,1,"y");
    MRIScopyMRI(tmpsurf,TrgVal,2,"z");
    if(center_surface)  MRIScenter(tmpsurf, tmpsurf);
    if(XYZScannerRAS) {
      printf("Converting output surface back to tkRegRAS\n");
      MRISscanner2Tkr(tmpsurf);
    }
    MRISwrite(tmpsurf, TrgValFile);
    if(npatches > 0) MRISfree(&tmpsurf);
  }
  else{
    printf("Writing %s\n",TrgValFile);
    err = 0;
    if(OutputCurvFormat){
      MRIScopyMRI(SurfReg[nsurfs-1], TrgVal, 0, "curv");
      err = MRISwriteCurvature(SurfReg[nsurfs-1],TrgValFile);
    }
    else
      err = MRIwrite(TrgVal,TrgValFile);
    if(err) {
      printf("ERROR: writing to %s\n",TrgValFile);
      exit(1);
    }
  }
  
  printf("mris_apply_reg done\n");
  exit(0);
}
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--norev")) ReverseMapFlag = 0;
    else if (!strcasecmp(option, "--no-rev")) ReverseMapFlag = 0;
    else if (!strcasecmp(option, "--nnf")) ReverseMapFlag = 0;
    else if (!strcasecmp(option, "--nnfr")) ReverseMapFlag = 1;
    else if (!strcasecmp(option, "--no-hash")) UseHash = 0;
    else if (!strcasecmp(option, "--jac")) DoJac = 1;
    else if (!strcasecmp(option, "--no-jac")) DoJac = 0;
    else if (!strcasecmp(option, "--randn")) DoSynthRand = 1;
    else if (!strcasecmp(option, "--ones")) DoSynthOnes = 1;
    else if (!strcasecmp(option, "--curv")) OutputCurvFormat=1;
    else if (!strcasecmp(option, "--center")) center_surface = 1 ;
    else if (!strcasecmp(option, "--reg-scanner-ras")) RegScannerRAS = 1 ;
    else if (!strcasecmp(option, "--xyz-scanner-ras")) XYZScannerRAS = 1 ;
    else if (!strcasecmp(option, "--scanner-ras")) {RegScannerRAS = 1 ;XYZScannerRAS = 1 ;}

    else if (!strcasecmp(option, "--src") || !strcasecmp(option, "--sval") || !strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      SrcValFile = pargv[0];
      if(!fio_FileExistsReadable(SrcValFile)){
	printf("ERROR: %s does not exist or is not readable by you\n",SrcValFile);
	exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--sval-annot") || !strcasecmp(option, "--src-annot")){
      if (nargc < 1) CMDargNErr(option,1);
      AnnotFile = pargv[0];
      if(!fio_FileExistsReadable(AnnotFile)){
	printf("ERROR: %s does not exist or is not readable by you\n",AnnotFile);
	exit(1);
      }
      DoJac = 0;
      ReverseMapFlag = 0;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--src-xyz")){
      if (nargc < 1) CMDargNErr(option,1);
      SurfXYZFile = pargv[0];
      if(!fio_FileExistsReadable(SurfXYZFile)){
	printf("ERROR: %s does not exist or is not readable by you\n",SurfXYZFile);
	exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--src-reg-scale")){
      if(nargc < 1) CMDargNErr(option,1);
      // Scale the coords of the first surface by SourceSurfRegScale. This
      // was implemented to make it easier to use CAT reg surfaces which
      // have a radius of 1.
      sscanf(pargv[0],"%lf",&SourceSurfRegScale);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--trg-reg-scale")){
      if(nargc < 1) CMDargNErr(option,1);
      // Scale the coords of the first surface by TargetSurfRegScale. This
      // was implemented to make it easier to use CAT reg surfaces which
      // have a radius of 1.
      sscanf(pargv[0],"%lf",&TargetSurfRegScale);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--sval-label") || !strcasecmp(option, "--src-label")){
      if (nargc < 1) CMDargNErr(option,1);
      LabelFile = pargv[0];
      if(!fio_FileExistsReadable(LabelFile)){
	printf("ERROR: %s does not exist or is not readable by you\n",LabelFile);
	exit(1);
      }
      DoJac = 0;
      ReverseMapFlag = 0;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--trg") || !strcasecmp(option, "--tval") 
	     || !strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      TrgValFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--streg") || !strcasecmp(option, "--st")) {
      if (nargc < 2) CMDargNErr(option,2);
      SurfRegFile[nsurfs] = pargv[0];
      nsurfs++;
      SurfRegFile[nsurfs] = pargv[1];
      nsurfs++;
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--patch")) {
      if (nargc < 2) CMDargNErr(option,2);
      SurfPatchFile[npatches] = pargv[0];
      npatches++;
      SurfPatchFile[npatches] = pargv[1];
      npatches++;
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--stvpair")) {
      if (nargc < 1) CMDargNErr(option,1);
      setenv("FS_MRISAPPLYREG_STVPAIR",pargv[0],1);
      /* Set up to create a text file with source-target vertex pairs (STVP) where the source
      is the fist surface and the target is the last surface. The format will be
        srcvtxno srcx srcy srcz trgvtxno trgx trgy trgz 
      The actual coordinates will come from TMP_VERTEX v->{tx,ty,tz}, so make sure those are set */
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--lta") || !strcasecmp(option, "--lta-rot")) {
      /*Automatically determines which direction the LTA goes by looking
	at the volume geometries of the LTA and the surface.*/
      if(nargc < 3) CMDargNErr(option,3);
      printf("Reading in %s\n",pargv[0]);
      MRIS *ltasurf = MRISread(pargv[0]);
      if(ltasurf==NULL) exit(1);
      printf("Reading in %s\n",pargv[1]);
      LTA *lta = LTAread(pargv[1]);
      if(lta==NULL) exit(1);
      if(!strcasecmp(option, "--lta-rot")){
	printf("Extracting rotational components\n");
	LTAmat2RotMat(lta);
      }
      int err = MRISltaMultiply(ltasurf, lta);
      if(err) exit(1);
      if (center_surface)
	MRIScenter(ltasurf, ltasurf);
      printf("Writing surf to %s\n",pargv[2]);
      err = MRISwrite(ltasurf,pargv[2]);
      if(err) exit(1);
      nargsused = 3;
      printf("mris_apply_reg done\n");
      exit(0);
    } 
    else if (!strcasecmp(option, "--lta-patch")) {
      // Surf Patch LTA out
      // The application for this option is exvivo reg between a whole
      // hemi and the WM surf of a block. This option allows a
      // registration to be applied to the block flat map to bring it
      // in reg with a flat map from the whole hemi. After that, then
      // the standard surf2surf reg can be applied (using --patch)
      if(nargc < 4) CMDargNErr(option,4);
      // surface needed to read in patch, actual coords not important
      printf("Reading in %s\n",pargv[0]);
      MRIS *ltasurf = MRISread(pargv[0]); 
      if(ltasurf==NULL) exit(1);
      ltasurf->vg.valid = 0; // needed for ltaMult
      // Now read in the (flat) patch
      printf("Reading in %s\n",pargv[1]);
      MRISreadPatch(ltasurf,pargv[1]);
      printf("Reading in %s\n",pargv[2]);
      /* Read in the LTA. The coordinates are a little confusing here.
      2D images of the curv are created with coords in the flat map
      surface coords (ie, tkreg). The reg comes registering these two
      maps with mri_coreg, so from the LTA standpoint, this is a RAS2RAS 
      eventhough secretly they are tkRAS. Because we are applying the reg
      to the surface xyz (tkRAS), we need to change the LTA coord designation
      to tkRAS so that MRISltaMultiply() does not change
      the matrix to tkRAS. To make matters worse, we have to invert the matrix
      beucase tkR matrices point from targ-to-source. It's a mess.*/
      LTA *lta = LTAread(pargv[2]);
      if(lta==NULL) exit(1);
      if(lta->type != LINEAR_RAS_TO_RAS){
	printf("Changing LTA type to RAS2RAS\n");
	LTAchangeType(lta,LINEAR_RAS_TO_RAS);
      }
      // The hack mentioned above
      lta->type = REGISTER_DAT; 
      MatrixInverse(lta->xforms[0].m_L,lta->xforms[0].m_L);
      int err = MRISltaMultiply(ltasurf, lta);
      // might could have just done this and skipped all the drama above
      //int err = MRISmatrixMultiply(ltasurf, lta->xforms[0].m_L); 
      if(err) exit(1);
      if (center_surface) MRIScenter(ltasurf, ltasurf);
      printf("Writing patch to %s\n",pargv[3]);
      err = MRISwritePatch(ltasurf,pargv[3]);
      if(err) exit(1);
      nargsused = 4;
      printf("mris_apply_reg done\n");
      exit(0);
    } 
    else if (!strcasecmp(option, "--reverse")) {
      if(nargc < 3) CMDargNErr(option,3);
      MRIS *revsurf = MRISread(pargv[0]);
      if(revsurf==NULL) exit(1);
      if(strcmp(pargv[1],"nopatch")!=0)	MRISreadPatch(revsurf,pargv[1]);
      MRISreverse(revsurf, REVERSE_X, 1);
      int err;
      if(strcmp(pargv[1],"nopatch")!=0)	err = MRISwritePatch(revsurf,pargv[2]);
      else                              err = MRISwrite(revsurf,pargv[2]);
      if(err) exit(1);
      nargsused = 3;
      printf("mris_apply_reg done\n");
      exit(0);
    }
    else if (!strcasecmp(option, "--map-vertex")) {
      // vertexno srcsurf trgsurf outfile
      if(nargc < 4) CMDargNErr(option,4);
      int srcvno,trgvno;
      sscanf(pargv[0],"%d",&srcvno);
      MRIS *srcsurf = MRISread(pargv[1]);
      if(!srcsurf) exit(1);
      MRIS *trgsurf = MRISread(pargv[2]);
      if(!trgsurf) exit(1);
      if(srcvno < 0 || srcvno >= srcsurf->nvertices){
	printf("srcvno=%d, out of range 0,%d\n",srcvno,srcsurf->nvertices);
	exit(1);
      }
      VERTEX *v = &(srcsurf->vertices[srcvno]);
      float dist;
      trgvno = MRISfindClosestVertex(trgsurf,v->x,v->y,v->z,&dist, CURRENT_VERTICES);
      VERTEX *v2 = &(srcsurf->vertices[trgvno]);
      printf("%d %8.4f %8.4f %8.4f  %8.4f  %d %8.4f %8.4f %8.4f\n",
	     trgvno,v2->x,v2->y,v2->z,dist,srcvno,v->x,v->y,v->z);
      if(strcmp(pargv[3],"nofile")!=0) {
	FILE *fp = fopen(pargv[3],"w");
	if(!fp){
	  printf("ERROR: could not open %s\n",pargv[3]);
	  exit(1);
	}
	fprintf(fp,"%d %8.4f %8.4f %8.4f  %8.4f  %d %8.4f %8.4f %8.4f\n",
		trgvno,v2->x,v2->y,v2->z,dist,srcvno,v->x,v->y,v->z);
	fclose(fp);
      }
      MRISfree(&srcsurf);
      MRISfree(&trgsurf);
      exit(0);
    }
    else if (!strcasecmp(option, "--gcam") || !strcasecmp(option, "--m3z") || !strcasecmp(option, "--inv-m3z") ){
      // Note: behavior is not changed with --inv-m3z now because it is automatically determined
      // whether the gcam needs to be inverted
      if(nargc < 3) CMDargNErr(option,3);
      printf("Reading in %s\n",pargv[0]);
      MRIS *m3zsurf = MRISread(pargv[0]);
      if(m3zsurf==NULL) exit(1);
      printf("Reading in %s\n",pargv[1]);
      GCA_MORPH *gcam;
      gcam = GCAMread(pargv[1]);
      if(gcam==NULL) exit(1);
      int err = GCAMmorphSurf(m3zsurf, gcam);
      if(err) exit(1);
      if(center_surface) MRIScenter(m3zsurf, m3zsurf);
      printf("Writing surf to %s\n",pargv[2]);
      err = MRISwrite(m3zsurf,pargv[2]);
      if(err) exit(1);
      nargsused = 3;
      printf("mris_apply_reg done\n");
      exit(0);
    } 
    else if (!strcasecmp(option, "--bci") || !strcasecmp(option, "--bci-xyz")){
      MRIS *SurfSrc=NULL;
      MRI *in=NULL;
      if(!strcasecmp(option, "--bci")){
	in = MRIread(pargv[0]);
	if(in == NULL) exit(1);
      }
      else {
	SurfSrc = MRISread(pargv[0]);
	if(SurfSrc==NULL)  exit(1);
	in = MRIcopyMRIS(NULL, SurfSrc, 2, "z"); // start at z to autoalloc
	MRIcopyMRIS(in, SurfSrc, 0, "x");
	MRIcopyMRIS(in, SurfSrc, 1, "y");
      }
      MRIS *srcreg = MRISread(pargv[1]);
      if(srcreg == NULL) exit(1);
      MRIS *trgreg = MRISread(pargv[2]);
      if(trgreg == NULL) exit(1);
      MRISsaveVertexPositions(srcreg,CANONICAL_VERTICES);
      MRISsaveVertexPositions(trgreg,CANONICAL_VERTICES);
      MRI *TrgVal = MRISapplyRegBCI(srcreg, trgreg, in);
      if(TrgVal == NULL) exit(1);
      int err=0;
      if(!strcasecmp(option, "--bci"))	err = MRIwrite(TrgVal,pargv[3]);
      else {
	MRIScopyMRI(trgreg,TrgVal,0,"x");
	MRIScopyMRI(trgreg,TrgVal,1,"y");
	MRIScopyMRI(trgreg,TrgVal,2,"z");
	copyVolGeom(&SurfSrc->vg, &trgreg->vg);// once used MRI, now seems silly
	err = MRISwrite(trgreg,pargv[3]);
      }
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
/*--------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*--------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf(" Input specifcation (pick one):\n");
  printf("   --src srcvalfile : source values (surface overlay) Can also use --i\n");
  printf("   --src-annot srcannotfile : source annotation (implies --no-rev)\n");
  printf("   --src-label labelfile : source label (implies --no-rev)\n");
  printf("   --src-xyz surfacefile : use xyz coords from given surface as input\n");
  printf(" Output specifcation (format depends on type of input):\n");
  printf("   --trg trgvalfile : (Can also use --o)\n");
  printf(" Registration specifcation (srcreg1->trgreg1->srcreg2->trgreg2...):\n");
  printf(" Need at least one --streg pair but can have any number\n");
  printf("   --streg srcreg1 trgreg1 : source and target reg files\n");
  printf("   --streg srcreg2 trgreg2 : more source and target reg files ...\n");
  printf("\n");
  printf("   --jac : use jacobian correction\n");
  printf("   --no-rev : do not do reverse mapping\n");
  printf("   --randn : replace input with WGN\n");
  printf("   --ones  : replace input with ones\n");
  printf("   --center  : place the center of the output surface at (0,0,0)\n");
  printf("   --curv  : save output in curv file format (spec full path)\n");
  printf("\n");
  printf("   --lta source-surf ltafile output-surf : apply LTA transform\n");
  printf("     automatically determins direction; other options do not apply to --lta\n");
  printf("   --lta-patch source-surf surfpatch ltafile output-patch : apply LTA transform to patch\n");
  printf("   --reverse surf patchopt output : LR reverse suface. patchopt=patchfile or nopatch\n");
  printf("   --patch srcregNpatch trgregNpatch : patches, one for each --streg\n");
  printf("   --stvpair stvpairfile : save srcvtxno srcx srcy srcz trgvtxno trgx trgy trgz from first and last surfs \n");
  printf("\n");
  printf("   --m3z source-surf m3zfile output-surf : apply m3z transform\n");
  printf("   --inv-m3z source-surf m3zfile output-surf : apply inverse m3z transform\n");
  printf("     except for --center, other options do not apply to --m3z or --inv-m3z\n");
  printf("\n");
  printf("   --src-reg-scale Scale : Scale the coords of the first surface\n");
  printf("   --trg-reg-scale Scale : Scale the coords of the last surface\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/*--------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  usage(stdout);
  exit(1) ;
}
/*--------------------------------------------------------------*/
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/*--------------------------------------------------------------*/
static void check_options(void) {
  int n;
  if(SrcValFile == NULL && AnnotFile == NULL && LabelFile == NULL && SurfXYZFile == NULL){
    printf("ERROR: need to specify source value file\n");
    exit(1);
  }
  if(SrcValFile && AnnotFile){
    printf("ERROR: cannot spec both --src and --src-annot\n");
    exit(1);
  }
  if(SrcValFile && LabelFile){
    printf("ERROR: cannot spec both --src and --src-label\n");
    exit(1);
  }
  if(AnnotFile && LabelFile){
    printf("ERROR: cannot spec both --src-annot and --src-label\n");
    exit(1);
  }
  if(AnnotFile && SurfXYZFile){
    printf("ERROR: cannot spec both --src-annot and --src-xyz\n");
    exit(1);
  }
  if(SrcValFile && SurfXYZFile){
    printf("ERROR: cannot spec both --src and --src-xyz\n");
    exit(1);
  }
  if(TrgValFile == NULL){
    printf("ERROR: need to specify target value file\n");
    exit(1);
  }
  if(nsurfs == 0){
    printf("ERROR: must specify at least one source:target registration pair\n");
    exit(1);
  }
  for(n=0; n<nsurfs;n++){
    if(!fio_FileExistsReadable(SurfRegFile[n])){
      printf("ERROR: %s does not exist or is not readable by you\n",SurfRegFile[n]);
      exit(1);
    }
  }
  if(npatches>0 && npatches != nsurfs){
    printf("ERROR: number of patches must match number of surfaces\n");
    exit(1);
  }
  return;
}
/*--------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"srcvalfile  %s\n",SrcValFile);
  fprintf(fp,"trgvalfile  %s\n",TrgValFile);
  fprintf(fp,"nsurfs  %d\n",nsurfs);
  fprintf(fp,"jac  %d\n",DoJac);
  fprintf(fp,"revmap  %d\n",ReverseMapFlag);
  return;
}

#include "mris_apply_reg.help.xml.h"
void usage(FILE *stream)
{
  outputHelpXml(mris_apply_reg_help_xml,mris_apply_reg_help_xml_len);
} /* end usage() */

/* EOF */

/*!
  \fn LABEL *MRISmask2Label(MRIS *surf, MRI *mask, int frame, double thresh)
  \brief Converts a surface overlay (mask) to a surface label. mask can be
  non-binary.  Values over thresh are used. The mask value is set to
  be the label stat.
 */
LABEL *MRISmask2Label(MRIS *surf, MRI *mask, int frame, double thresh)
{
  LABEL *label;
  int c, r, s, n, vtxno;
  double val;
  VERTEX *v;

  // Count number of points in the label
  n = 0;
  for(s=0; s < mask->depth; s++){
    for(r=0; r < mask->height; r++){
      for(c=0; c < mask->width; c++){
	val = MRIgetVoxVal(mask,c,r,s,frame);
	if(val < thresh) continue;
	n++;
      }
    }
  }

  // Alloc
  label = LabelAlloc(n,surf->subject_name,NULL);
  label->n_points = n;

  // Asign values
  n = 0;
  vtxno = -1;
  for(s=0; s < mask->depth; s++){
    for(r=0; r < mask->height; r++){
      for(c=0; c < mask->width; c++){
	vtxno++;
	val = MRIgetVoxVal(mask,c,r,s,frame);
	if(val < thresh) continue;
	v = &surf->vertices[vtxno];
	label->lv[n].vno = vtxno;
	label->lv[n].x = v->x;
	label->lv[n].y = v->y;
	label->lv[n].z = v->z;
	label->lv[n].stat = val;
	n++;
      }
    }
  }

  return(label);
}

// A simple function to muliply the surf coords
int ApplyScaleSurf(MRIS *surf, const double scale)
{
  int vtxno;
  VERTEX *v;
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    v = &(surf->vertices[vtxno]);
    v->x *= scale;
    v->y *= scale;
    v->z *= scale;
  }
  return(0);
}
