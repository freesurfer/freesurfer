/**
 * @brief utility for sampling volumes onto a surface.
 *
 * Purpose: Resamples a volume onto a surface. The surface
 * may be that of a subject other than the source subject.
 * This replaces paint.c and is also used to convert functional
 * data onto the icosahedron.
 * Volume-to-Volume - V2V is a necessary step when converting functional
 * data to talairach or painting onto the surface. The model as of 2/4/01
 * assumes that the transformation is completely linear from one space to
 * another, though there are variables for intermediate transformations
 * built in. The four variables are: (1) Quantization matrix Q, (2) Field-
 * of-view matrix F, (3) Warping matrix W, and (4), Registration matrix D.
 *
 * D - Registration matrix. Converts from Anatomical Space to Unwarpded
 *     Scanner Space.
 * W - Warping matrix. Converts from Unwarpded Scanner Space to Warped
 *     Scanner Space.
 * F - FOV matrix. Converts from Warpded Scanner Space to Field-of-View
 *     Space (the space created by the axes centered in the FOV and
 *     parallel with the edges of the FOV).
 * Q - Quantization matrix. Converts from FOV space to ColRowSlice Index
 *     Space.
 */
/*
 * Original Author: Doug Greve
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "icosahedron.h"
#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "mri.h"
#include "mri_identify.h"
#include "mri2.h"
#include "prime.h"
#include "fsenv.h"
#include "registerio.h"
#include "resample.h"
#include "selxavgio.h"
#include "version.h"
#include "fmriutils.h"
#include "proto.h" // nint
#include "cmdargs.h"

#ifndef FZERO
#define FZERO(f)     (fabs(f) < 0.0000001F)
#endif

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
int main(int argc, char *argv[]) ;


const char *Progname = NULL;

static char *defaulttypestring;
static int  defaulttype = MRI_VOLUME_TYPE_UNKNOWN;

static char *srcvolid   = NULL;
static char *srctypestring = NULL;
static int   srctype = MRI_VOLUME_TYPE_UNKNOWN;
static char *srcregfile = NULL;
static int  regheader = 0;
static char *srcwarp    = NULL;
static int   srcoldreg  = 0;
static char *srcsubject = NULL;
static char *srcsubjectuse = NULL;

static const char *ref_vol_name = "orig.mgz" ;

static char *srchitvolid   = NULL;
static char *srchittypestring = NULL;
static int   srchittype = MRI_VOLUME_TYPE_UNKNOWN;

static char *hemi    = NULL;
static char const *surfname = "white";
static char *trgsubject = NULL;
static int  IcoOrder = -1;
static float IcoRadius = 100;
static char const *surfreg = "sphere.reg";
static char const *thicknessname = "thickness";
static float ProjFrac = 0;
static int   ProjOpt = 0 ;
static char  *volume_fraction_fname = NULL ;
static int   ProjDistFlag = 0;
static float ProjFracMin=0.0,ProjFracMax=0.0,ProjFracDelta=1.0;

MRI *build_sample_array(MRI_SURFACE *mris, MRI *mri, MATRIX *m, 
                        float din, float dout, int nsamples, 
                        MRI *mri_wm, MRI *mri_gm, MRI *mri_csf);
MRI *estimate_gm_values(MRI *mri_wm, MRI *mri_gm, MRI *mri_csf,
                        MRI *SrcVol,
                        MATRIX *Qsrc, MATRIX *Fsrc, MATRIX *Wsrc, MATRIX *Dsrc,
                        MRI_SURFACE *TrgSurf,
                        int InterpMethod, int float2int, MRI *SrcHitVol) ;


static MRI_SURFACE *Surf    = NULL;
static MRI_SURFACE *SurfOut = NULL;
static MRI_SURFACE *SrcSurfReg = NULL;
static MRI_SURFACE *TrgSurfReg = NULL;
static int UseHash = 1;

static char *outfile  = NULL;
static char *outtypestring = NULL;
static int  outtype = MRI_VOLUME_TYPE_UNKNOWN;

static char *srchitfile = NULL;
static char *trghitfile = NULL;

static const char  *interpmethod_string = "nearest";
static int  interpmethod = -1;
static const char *mapmethod = "nnfr";

static int debug = 0;
static int reshape = 0;
static int reshapefactor = 0;
static int reshapetarget = 20;
static int reshape3d = 0;

static MATRIX *Dsrc, *Dsrctmp, *Wsrc, *Fsrc, *Qsrc, *vox2ras;

static char *SUBJECTS_DIR = NULL;
static MRI *SrcVol, *SurfVals, *SurfVals2, *SurfValsP;
static MRI *SrcHits, *SrcDist, *TrgHits, *TrgDist;
static MRI *mritmp;
static MRI *SrcHitVol;
static MRI *TargVol;

static FILE *fp;

static char *nvoxfile = NULL;

static char tmpstr[2000];

static int  float2int_src;
static const char *float2int_string = "round";
static int  float2int = -1;
static int  fixtkreg = 0;
static int ReverseMapFlag = 0;

static int framesave = -1;

static float fwhm = 0, gstd = 0;
static float surf_fwhm = 0, surf_gstd = 0;

static int  srcsynth = 0;
int srcsynthindex = 0;
static long seed = -1; /* < 0 for auto */
static char *seedfile = NULL;

static double scale = 0;
static int GetProjMax = 0;
int UseCortexLabel = 0;
static char  *mask_label_name = NULL ;
LABEL *area ;

double angles[3] = {0,0,0};
MATRIX *Mrot = NULL;
double xyztrans[3] = {0,0,0};
MATRIX *Mtrans = NULL;

char *vsmfile = NULL;
MRI *vsm = NULL;
int UseOld = 1;
MRI *MRIvol2surf(MRI *SrcVol, MATRIX *Rtk, MRI_SURFACE *TrgSurf, 
		 MRI *vsm, int InterpMethod, MRI *SrcHitVol, 
		 float ProjFrac, int ProjType, int nskip);

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
int main(int argc, char **argv) {
  int n,err, f, vtx, svtx, tvtx, nproj, nSmoothSteps;
  int nrows_src, ncols_src, nslcs_src, nfrms;
  //float ipr, bpr, intensity;
  float colres_src=0, rowres_src=0, slcres_src=0;
  char fname[2000];
  int nTrg121,nSrc121,nSrcLost;
  int nTrgMulti,nSrcMulti;
  float MnTrgMultiHits,MnSrcMultiHits;
  int nargs;
  int r,c,s,nsrchits;
  LTA *lta;

  nargs = handleVersionOption(argc, argv, "mri_vol2surf");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  vg_isEqual_Threshold = 10e-4;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL) {
    fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  /* voxel indices conversion from float to integer */
  if (float2int < 0) float2int = float2int_src;
  else {
    if (float2int != float2int_src) {
      printf("INFO: float2int on the command line (%d) overrides that \n"
             "      in the registration file (%d).\n",float2int,
             float2int_src);
    }
  }
  printf("INFO: float2int code = %d\n",float2int);

  if(srcsynth == 0 && srcsynthindex == 0) {
    /* Load the Source Volume */
    SrcVol =  MRIreadType(srcvolid,srctype);
    if (SrcVol == NULL) {
      printf("ERROR: could not read %s as type %d\n",srcvolid,srctype);
      exit(1);
    }
    if (SrcVol->type != MRI_FLOAT) {
      printf("INFO: changing type to float\n");
      SrcVol = MRISeqchangeType(SrcVol,MRI_FLOAT,0,0,0);
    }
    printf("Done loading volume\n");
  }
  if(srcsynth){
    /* Synth the Source Volume */
    printf("Synthesizing, seed = %ld\n",seed);
    srand48(seed);
    //srcfmtid = checkfmt(srctype);
    mritmp =  MRIreadType(srcvolid,srctype);
    if (mritmp == NULL) {
      printf("ERROR: could not read %s as type %d\n",srcvolid,srctype);
      exit(1);
    }
    SrcVol = MRIrandn(mritmp->width, mritmp->height,
                      mritmp->depth, mritmp->nframes, 0, 1, NULL);
    MRIcopyHeader(mritmp, SrcVol);
    SrcVol->type = MRI_FLOAT;
    MRIfree(&mritmp);
  }
  if(srcsynthindex){
    printf("Synthesizing with index\n");
    mritmp = MRIreadType(srcvolid,srctype);
    if (mritmp == NULL) {
      printf("ERROR: could not read %s as type %d\n",srcvolid,srctype);
      exit(1);
    }
    SrcVol = MRIindexNo(mritmp,NULL);
  }

  if (!regheader) {
    /* Load the registration matrix */
    lta = LTAread(srcregfile);
    if(lta == NULL) exit(1);
    if(debug) LTAprint(stdout,lta);
    srcsubject = lta->subject;
    if(srcsubjectuse) {
      if(strcmp(srcsubject,srcsubjectuse)){
        printf("INFO: overriding regsubject %s with %s\n",srcsubject,srcsubjectuse);
      }
      srcsubject = srcsubjectuse;
    }
    if(lta->xforms[0].src.valid && lta->xforms[0].dst.valid){
      printf("Input reg is LTA\n");
      if(LTAmriIsTarget(lta,SrcVol)) {
	printf("Inverting LTA\n");
	lta = LTAinvert(lta,lta);        
      }
      else if(!LTAmriIsSource(lta,SrcVol)) {
	printf("\n\n");
	LTAprint(stdout,lta);
	printf("\n\n");
	printf("ERROR: source volume is neither source nor target of the registration\n");
	exit(1);
      }
      LTAchangeType(lta, REGISTER_DAT);
    }
    else {
      printf("Input reg is register.dat\n");
      /* check that the resolution of the SrcVol is the same as in reg file */
      if(fabs(lta->xforms[0].src.xsize-SrcVol->xsize) > .001 ||
	 fabs(lta->xforms[0].src.zsize-SrcVol->zsize) > .001){
	printf("WARNING: the voxel resolution in the source volume"
	       " (%g,%g,%g) differs \n",
	       SrcVol->xsize,SrcVol->ysize,SrcVol->zsize);
	printf("         from that listed in the registration file (%g,%g,%g)\n",
	       lta->xforms[0].src.xsize,lta->xforms[0].src.ysize,lta->xforms[0].src.zsize);
	if(fixtkreg) {
	  printf("ERROR: cannot fix tkreg matrix with resolution inconsistency\n");
	  exit(1);
	}
      }
    }
    Dsrc = lta->xforms[0].m_L;
  } else {
    /* compute the registration from the header */
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,srcsubject, ref_vol_name);
    printf("Computing registration from header.\n");
    printf("  Using %s as target reference.\n",tmpstr);
    TargVol = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
    if (TargVol == NULL) exit(1);
    Dsrc = MRItkRegMtx(TargVol,SrcVol,NULL);
    MRIfree(&TargVol);
  }

  if(Mrot){
    printf("Applying rotation matrix (R=M*R)\n");
    printf("Current Reg Matrix is:\n");
    MatrixPrint(stdout,Dsrc);
    printf("  Angles (deg): %lf %lf %lf\n",angles[0]*180/M_PI,angles[1]*180/M_PI,angles[2]*180/M_PI);
    printf("  Angles (rad): %lf %lf %lf\n",angles[0],angles[1],angles[2]);
    printf("  Rotation matrix:\n");
    MatrixPrint(stdout,Mrot);
    Dsrc = MatrixMultiply(Mrot,Dsrc,Dsrc);
  }

  if(Mtrans){
    printf("Applying translation matrix (R=M*R)\n");
    printf("Current Reg Matrix is:\n");
    MatrixPrint(stdout,Dsrc);
    printf("  Trans (mm): %lf %lf %lf\n",xyztrans[0],xyztrans[1],xyztrans[2]);
    printf("  Translation matrix:\n");
    MatrixPrint(stdout,Mtrans);
    Dsrc = MatrixMultiply(Mtrans,Dsrc,Dsrc);
  }

  ncols_src = SrcVol->width;
  nrows_src = SrcVol->height;
  nslcs_src = SrcVol->depth;
  nfrms = SrcVol->nframes;
  colres_src = SrcVol->xsize; /* in-plane resolution */
  rowres_src = SrcVol->ysize; /* in-plane resolution */
  slcres_src = SrcVol->zsize; /* between-plane resolution */

  if (framesave >= nfrms) {
    printf("ERROR: frame = %d, input volume limits to < %d\n",framesave,nfrms);
    exit(1);
  }

  /* Fix tkregister matrix if necessary */
  if(fixtkreg && float2int == FLT2INT_TKREG) {
    printf("INFO: fixing tkregister matrix\n");
    Dsrctmp = MRIfixTkReg(SrcVol,Dsrc);
    printf("-------- original matrix -----------\n");
    MatrixPrint(stdout,Dsrc);
    printf("-------- fixed matrix -----------\n");
    MatrixPrint(stdout,Dsrctmp);
    MatrixFree(&Dsrc);
    Dsrc = Dsrctmp;
    float2int = FLT2INT_ROUND;
  }

  printf("-------- original matrix -----------\n");
  MatrixPrint(stdout,Dsrc);
  printf("-------- original matrix -----------\n");

  if(vsmfile){
    printf("Reading vsm %s\n",vsmfile);
    vsm =  MRIread(vsmfile);
    if(vsm == NULL) {
      printf("ERROR: could not read %s\n",vsmfile);
      exit(1);
    }
    err = MRIdimMismatch(vsm,SrcVol,0);
    if(err){
      printf("ERROR: vsm dimension mismatch %d\n",err);
      exit(1);
    }
  }

  /* Wsrc: Get the source warping Transform */
  Wsrc = NULL;
  /* Fsrc: Get the source FOV registration matrix */
  Fsrc = NULL;
  // Compute vox2ras for source
  vox2ras = MRIxfmCRS2XYZtkreg(SrcVol);
  // Compute ras2vox (Qsrc: the quantization matrix)
  Qsrc = MatrixInverse(vox2ras,NULL);

  if (fwhm > 0) {
    printf("INFO: smoothing volume at fwhm = %g mm (std = %g)\n",fwhm,gstd);
    MRIgaussianSmooth(SrcVol, gstd, 1, SrcVol); /* 1 = normalize */
  }

  if(trgsubject == NULL) trgsubject = srcsubject;
  if(UseCortexLabel) {
    sprintf(tmpstr,"%s/%s/label/%s.cortex.label",SUBJECTS_DIR,trgsubject,hemi);
    mask_label_name = strcpyalloc(tmpstr);
  }
  if(mask_label_name){
    printf("Loading label %s\n",mask_label_name) ;
    area = LabelRead(NULL, mask_label_name) ;
    if(area == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load label file %s", 
		Progname, mask_label_name);
  }

  /* Load the surface for subject indicated by registration file*/
  sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,surfname);
  printf("Reading surface %s\n",fname);
  Surf = MRISread(fname) ;
  if (!Surf)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;
  printf("Done reading source surface\n");
  fflush(stdout);

  /* Load the thickness for projection along the normal*/
  if ((ProjFrac != 0 && ProjDistFlag == 0) || (ProjOpt != 0)) {
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,
            hemi,thicknessname);
    printf("Reading thickness %s\n",fname);
    err = MRISreadCurvatureFile(Surf, fname);
    if (err) exit(1);
    printf("Done\n");
  }

  SrcHitVol = MRIallocSequence(SrcVol->width,SrcVol->height,
                               SrcVol->depth,MRI_FLOAT,1);
  if (SrcHitVol == NULL) {
    printf("ERROR: could not alloc SrcHitVol\n");
    exit(1);
  }
  MRIcopyHeader(SrcVol,SrcHitVol);

  /* Map the values from the volume to the surface */
  printf("Mapping Source Volume onto Source Subject Surface\n");
  fflush(stdout);
  if (ProjOpt != 0)
  {
    char fname[STRLEN] ;
    MRI  *mri_gm, *mri_wm, *mri_csf ;
    MATRIX *Qsrc, *QFWDsrc ;
    
    sprintf(fname, "%s.gm.mgz", volume_fraction_fname) ;
    printf("reading gm volume fraction from %s\n", fname) ;
    mri_gm = MRIread(fname) ;
    if (mri_gm == NULL)
      ErrorExit(ERROR_NOFILE, "could not read gm volume fraction from %s", 
                fname);
    
    sprintf(fname, "%s.wm.mgz", volume_fraction_fname) ;
    printf("reading wm volume fraction from %s\n", fname) ;
    mri_wm = MRIread(fname) ;
    if (mri_wm == NULL)
      ErrorExit(ERROR_NOFILE, "could not read wm volume fraction from %s", 
                fname);
    sprintf(fname, "%s.csf.mgz", volume_fraction_fname) ;
    printf("reading csf volume fraction from %s\n", fname) ;
    mri_csf = MRIread(fname) ;
    if (mri_csf == NULL)
      ErrorExit(ERROR_NOFILE, "could not read csf volume fraction from %s", 
                fname);
    Qsrc = MRIxfmCRS2XYZtkreg(SrcVol);
    Qsrc = MatrixInverse(Qsrc,Qsrc);
    QFWDsrc = ComputeQFWD(Qsrc,Fsrc,Wsrc,Dsrc,NULL);
    SurfVals = build_sample_array(Surf, SrcVol, QFWDsrc, 1.0, 0.1, 10,
                                  mri_wm, mri_gm, mri_csf) ;
    MatrixFree(&Qsrc) ; MatrixFree(&QFWDsrc) ;
  }
  else
  {
    nproj = 0;
    for (ProjFrac=ProjFracMin; 
         ProjFrac <= ProjFracMax; 
         ProjFrac += ProjFracDelta) {
      printf("%2d %g %g %g\n",nproj+1,ProjFrac,ProjFracMin,ProjFracMax);
      if(UseOld){
        printf("using old\n");
        SurfValsP = 
          vol2surf_linear(SrcVol, Qsrc, Fsrc, Wsrc, Dsrc,
                          Surf, ProjFrac, interpmethod, float2int, SrcHitVol,
                          ProjDistFlag, 1);
      }
      else{
        printf("using new\n");
        SurfValsP = 
          MRIvol2surfVSM(SrcVol, Dsrc, Surf, vsm, interpmethod, SrcHitVol, 
                         ProjFrac, ProjDistFlag,1,NULL);
      }
      fflush(stdout);
      if (SurfValsP == NULL) {
        printf("ERROR: mapping volume to source\n");
        exit(1);
      }
      if (nproj == 0) SurfVals = MRIcopy(SurfValsP,NULL);
      else {
        if (!GetProjMax) MRIadd(SurfVals,SurfValsP,SurfVals);
        else            MRImax(SurfVals,SurfValsP,SurfVals);
      }
      MRIfree(&SurfValsP);
      nproj ++;
    }
    if (!GetProjMax) MRImultiplyConst(SurfVals, 1.0/nproj, SurfVals);
  }

  printf("Done mapping volume to surface\n");
  fflush(stdout);
  MRIfree(&SrcVol);

  /* Count the number of source voxels hit */
  nsrchits = 0;
  for (c=0; c < SrcHitVol->width; c++) {
    for (r=0; r < SrcHitVol->height; r++) {
      for (s=0; s < SrcHitVol->depth; s++) {
        if (MRIFseq_vox(SrcHitVol,c,r,s,0) > 0.5) nsrchits++;
      }
    }
  }
  printf("Number of source voxels hit = %d\n",nsrchits);
  if (nvoxfile != NULL) {
    fp = fopen(nvoxfile,"w");
    if (fp == NULL)
      printf("ERROR: could not open nvox file %s\n",nvoxfile);
    else {
      fprintf(fp,"%d\n",nsrchits);
      fclose(fp);
    }
  }

  if (srchitvolid != NULL) {
    printf("Saving src hit volume.\n");
    MRIwriteType(SrcHitVol,srchitvolid,srchittype);
  }
  MRIfree(&SrcHitVol);

  if (trgsubject != NULL && strcmp(trgsubject,srcsubject)) {
    /* load in the source subject registration */
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,surfreg);
    printf("Reading source surface registration \n  %s\n",fname);
    SrcSurfReg = MRISread(fname);
    if (SrcSurfReg==NULL) {
      printf("ERROR: reading %s\n",fname);
      exit(1);
    }
    printf("Done loading source registration surface\n");
    fflush(stdout);

    /* load in the target subject registration */
    if (strcmp(trgsubject,"ico")) { /* not ico */
      sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,hemi,surfreg);
      printf("Reading target registration \n   %s\n",fname);
      TrgSurfReg = MRISread(fname);
      if (TrgSurfReg == NULL) {
        printf("ERROR: could not read %s\n",fname);
        exit(1);
      }
    } else { /* is ico */
      if (reshapefactor == 0) reshapefactor = 6; /* 6 slices for ico target */
      fflush(stdout);
      printf("Reading icosahedron, order = %d, radius = %g\n",
             IcoOrder,IcoRadius);
      TrgSurfReg = ReadIcoByOrder(IcoOrder,IcoRadius);
      if (TrgSurfReg==NULL) {
        printf("ERROR reading icosahedron\n");
        exit(1);
      }
    }
    printf("Done loading target registration surface\n");
    fflush(stdout);

    if (!strcmp(mapmethod,"nnfr")) ReverseMapFlag = 1;
    else                          ReverseMapFlag = 0;

    printf("Mapping Surfaces (%s -> %s)\n",srcsubject,trgsubject);
    fflush(stdout);
    SurfVals2 = surf2surf_nnfr(SurfVals, SrcSurfReg,TrgSurfReg,
                               &SrcHits,&SrcDist,&TrgHits,&TrgDist,
                               ReverseMapFlag,UseHash);
    if(SurfVals2 == NULL) {
      printf("ERROR: mapping surfaces\n");
      exit(1);
    }
    printf("Done mapping surfaces\n");
    fflush(stdout);

    /*Compute some stats on mapping number 
      of trgvtxs mapped from a source vtx */
    nSrc121 = 0;
    nSrcLost = 0;
    MnSrcMultiHits = 0.0;
    for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
      n = MRIFseq_vox(SrcHits,svtx,0,0,0);
      if (n == 1) nSrc121++;
      else if (n == 0) nSrcLost++;
      else MnSrcMultiHits += n;
    }
    nSrcMulti = SrcSurfReg->nvertices - nSrc121;
    if (nSrcMulti > 0) MnSrcMultiHits = (MnSrcMultiHits/nSrcMulti);
    else              MnSrcMultiHits = 0;
    printf("nSrc121 = %5d, nSrcLost = %5d, "
           "nSrcMulti = %5d, MnSrcMultiHits = %g\n",
           nSrc121,nSrcLost,nSrcMulti,MnSrcMultiHits);
    MRISfree(&SrcSurfReg);

    nTrg121 = 0;
    MnTrgMultiHits = 0.0;
    for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
      n = MRIFseq_vox(TrgHits,tvtx,0,0,0);
      if (n == 1) nTrg121++;
      else MnTrgMultiHits += n;
    }
    nTrgMulti = TrgSurfReg->nvertices - nTrg121;
    if (nTrgMulti > 0) MnTrgMultiHits = (MnTrgMultiHits/nTrgMulti);
    else              MnTrgMultiHits = 0;
    printf("nTrg121 = %5d, nTrgMulti = %5d, MnTrgMultiHits = %g\n",
           nTrg121,nTrgMulti,MnTrgMultiHits);

    /* save the Source Hits into a .w file */
    if (srchitfile != NULL) {
      for (vtx = 0; vtx < Surf->nvertices; vtx++)
        Surf->vertices[vtx].val = MRIFseq_vox(SrcHits,vtx,0,0,0) ;
      MRISwriteValues(Surf, srchitfile) ;
      MRIfree(&SrcHits);
    }
    /* save the Target Hits into a .w file */
    if (trghitfile != NULL) {
      for (vtx = 0; vtx < SurfOut->nvertices; vtx++)
        SurfOut->vertices[vtx].val = MRIFseq_vox(TrgHits,vtx,0,0,0) ;
      MRISwriteValues(SurfOut, trghitfile) ;
      MRIfree(&TrgHits);
    }
    SurfOut = TrgSurfReg;
  } else {
    // Source and target subjects are the same
    SurfVals2 = SurfVals;
    SurfOut = Surf;
  }

  if(surf_fwhm > 0){
    nSmoothSteps = MRISfwhm2nitersSubj(surf_fwhm,trgsubject,hemi,surfname);
    if(nSmoothSteps == -1) exit(1);
    printf("Surface smoothing by fwhm = %g (n=%d)\n",surf_fwhm,nSmoothSteps);
    fflush(stdout);
    MRISsmoothMRI(SurfOut, SurfVals2, nSmoothSteps, NULL, SurfVals2);
  }

  if(scale != 0) {
    printf("Rescaling output by %g\n",scale);
    MRImultiplyConst(SurfVals2,scale,SurfVals2);
  }

  if(framesave > 0){
    mritmp = fMRIframe(SurfVals2, framesave, NULL);
    if(mritmp == NULL) exit(1);
    MRIfree(&SurfVals2);
    SurfVals2 = mritmp;
  }

  if(mask_label_name) {
    printf("Masking with %s\n",mask_label_name);
    MRISclearMarks(SurfOut) ;
    LabelMarkSurface(area, SurfOut) ;
    for (vtx = 0; vtx < SurfVals2->width; vtx++){
      if (SurfOut->vertices[vtx].marked == 0){
	for(f=0; f < SurfVals2->nframes; f++)
	  MRIsetVoxVal(SurfVals2, vtx, 0, 0, f, 0.0) ;
      }
    }
  }

  if (outtypestring != NULL &&
      (!strcasecmp(outtypestring,"w") || !strcasecmp(outtypestring,"paint"))) {
    /*-------------- paint or .w --------------*/
    for(vtx = 0; vtx < SurfVals2->width; vtx++)
      SurfOut->vertices[vtx].val = MRIFseq_vox(SurfVals2,vtx,0,0,0) ;
    if(mask_label_name) LabelMaskSurface(area, Surf) ;
    MRISwriteValues(SurfOut, outfile) ;
  } else {
    if (reshape) {
      if (reshapefactor == 0) {
        if (outtype == MRI_ANALYZE4D_FILE || outtype == MRI_ANALYZE_FILE ||
            outtype == NIFTI1_FILE || outtype == NII_FILE)
          reshapefactor = GetClosestPrimeFactorLess(SurfVals2->width,
                                                    32768);
        else
          reshapefactor = GetClosestPrimeFactor(SurfVals2->width,
                                                reshapetarget);
      }
      if (reshapefactor != SurfVals2->width) {
        printf("Reshaping %d (nvertices = %d)\n",
               reshapefactor,SurfVals2->width);
        mritmp = mri_reshape(SurfVals2, SurfVals2->width / reshapefactor,
                             1, reshapefactor,SurfVals2->nframes);
        if (mritmp == NULL) {
          printf("ERROR: mri_reshape could not alloc\n");
          return(1);
        }
        MRIfree(&SurfVals2);
        SurfVals2 = mritmp;
      } else {
        printf("INFO: nvertices is prime, cannot reshape\n");
      }
    }

    if (reshape3d) {
      if(TrgSurfReg->nvertices != 163842) {
	printf("ERROR: subject must have 163842 vertices to 3d reshape\n");
	exit(1);
      }
      printf("Reshape 3d\n");
      mritmp = mri_reshape(SurfVals2, 42,47,83,SurfVals2->nframes);
      if (mritmp == NULL) {
        printf("ERROR: mri_reshape could not alloc\n");
        return(1);
      }
      MRIfree(&SurfVals2);
      SurfVals2 = mritmp;
    }

    printf("Writing to %s\n",outfile);
    printf("Dim: %d %d %d\n",
           SurfVals2->width,SurfVals2->height,SurfVals2->depth);
    err = MRIwriteType(SurfVals2,outfile,outtype);
    if(err){
      printf("ERROR: saving %s\n",outfile);
      exit(1);
    }
  }

  return(0);
}
/* --------------------------------------------- */
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
    else if (!strcasecmp(option, "--usehash")) UseHash = 1;
    else if (!strcasecmp(option, "--hash")) UseHash = 1;
    else if (!strcasecmp(option, "--dontusehash")) UseHash = 0;
    else if (!strcasecmp(option, "--nohash")) UseHash = 0;
    else if (!strcasecmp(option, "--reshape"))    reshape = 1;
    else if (!strcasecmp(option, "--reshape3d")) {
      reshape3d = 1;reshape = 0;
    }
    else if (!strcasecmp(option, "--noreshape"))  reshape = 0;
    else if (!strcasecmp(option, "--no-reshape")) reshape = 0;
    else if (!strcasecmp(option, "--fixtkreg")) fixtkreg = 1;
    else if (!strcasecmp(option, "--nofixtkreg")) fixtkreg = 0;
    else if (!strcasecmp(option, "--v")) 
    {
      Gdiag_no = atoi(pargv[0]);
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--ref")) 
    {
      ref_vol_name = pargv[0];
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--inflated"))  surfname = "inflated";

    else if ( !strcmp(option, "--sd") ) {
      if (nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      FSENVsetSUBJECTS_DIR(SUBJECTS_DIR);
      nargsused = 1;
    } else if ( !strcmp(option, "--default_type") ) {
      if (nargc < 1) argnerr(option,1);
      defaulttypestring = pargv[0];
      defaulttype = string_to_type(defaulttypestring);
      nargsused = 1;
    }
    /* -------- source volume inputs ------ */
    else if (!strcmp(option, "--srcvol") ||
             !strcmp(option, "--src") || !strcmp(option, "--mov")) {
      if (nargc < 1) argnerr(option,1);
      srcvolid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--src_type") ||
               !strcmp(option, "--srcvol_type") ||
               !strcmp(option, "--srcfmt")) {
      if (nargc < 1) argnerr(option,1);
      srctypestring = pargv[0];
      srctype = string_to_type(srctypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--srchitvol") ||
               !strcmp(option, "--srchit")) {
      if (nargc < 1) argnerr(option,1);
      srchitvolid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srchit_type") ||
               !strcmp(option, "--srchitvol_type") ||
               !strcmp(option, "--srchitfmt")) {
      if (nargc < 1) argnerr(option,1);
      srchittypestring = pargv[0];
      srchittype = string_to_type(srchittypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--reg") || !strcmp(option, "--srcreg")) {
      if (nargc < 1) argnerr(option,1);
      srcregfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--regheader")) {
      if (nargc < 1) argnerr(option,1);
      regheader = 1;
      srcsubject = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--vg-thresh")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&vg_isEqual_Threshold);
      nargsused = 1;
    }
    else if (!strcmp(option, "--mni152reg")) {
      sprintf(tmpstr,"%s/average/mni152.register.dat",
	      getenv("FREESURFER_HOME"));
      srcregfile = strcpyalloc(tmpstr);
    } 
    else if (!strcmp(option, "--rot")) {
      if (nargc < 3) argnerr(option,3);
      // Angles are in degrees
      sscanf(pargv[0],"%lf",&angles[0]);
      sscanf(pargv[1],"%lf",&angles[1]);
      sscanf(pargv[2],"%lf",&angles[2]);
      angles[0] *= (M_PI/180);
      angles[1] *= (M_PI/180);
      angles[2] *= (M_PI/180);
      Mrot = MRIangles2RotMat(angles);
      nargsused = 3;
    } else if (!strcmp(option, "--trans")) {
      if (nargc < 3) argnerr(option,3);
      // Translation in mm
      sscanf(pargv[0],"%lf",&xyztrans[0]);
      sscanf(pargv[1],"%lf",&xyztrans[1]);
      sscanf(pargv[2],"%lf",&xyztrans[2]);
      Mtrans = MatrixIdentity(4,NULL);
      Mtrans->rptr[1][4] = xyztrans[0];
      Mtrans->rptr[2][4] = xyztrans[1];
      Mtrans->rptr[3][4] = xyztrans[2];
      nargsused = 3;
    } else if (!strcmp(option, "--srcoldreg")) {
      srcoldreg = 1;
    } else if (!strcmp(option, "--srcwarp")) {
      if (nargc < 1) argnerr(option,1);
      srcwarp = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--frame")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&framesave);
      nargsused = 1;
    } else if (!strcmp(option, "--surf")) {
      if (nargc < 1) argnerr(option,1);
      surfname = pargv[0];
      if(!strcmp(surfname,"inflated")){
	printf("ERROR: you have chosen the inflated surface.\n");
	printf("This is most likely not what you want to do.\n");
	printf("If you really want to use the inflated surface,\n");
	printf("re-run this program with --inflated instead of\n");
	printf("--surf inflated.\n");
	exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--hemi")) {
      if (nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--mapmethod")) {
      if (nargc < 1) argnerr(option,1);
      mapmethod = pargv[0];
      if (strcmp(mapmethod,"nnfr") && strcmp(mapmethod,"nnf")) {
        fprintf(stderr,"ERROR: mapmethod must be nnfr or nnf\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--trgsubject")) {
      if (nargc < 1) argnerr(option,1);
      trgsubject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srcsubject")) {
      if (nargc < 1) argnerr(option,1);
      srcsubjectuse = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--icoorder")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&IcoOrder);
      printf("IcoOrder = %d, nIcoVtxs = %d\n",IcoOrder,
             IcoNVtxsFromOrder(IcoOrder));
      nargsused = 1;
    } else if (!strcmp(option, "--surfreg")) {
      if (nargc < 1) argnerr(option,1);
      surfreg = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--projfrac")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&ProjFrac);
      ProjFracMin=ProjFrac;
      ProjFracMax=ProjFrac;
      ProjFracDelta=1.0;
      nargsused = 1;
    } else if (!strcmp(option, "--projopt")) {
      if (nargc < 1) argnerr(option,1);
      volume_fraction_fname = pargv[0] ;
      ProjOpt = 1 ;
      nargsused = 1;
    } else if (!strcmp(option, "--projfrac-int") ||
               !strcmp(option, "--projfrac-avg")) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%f",&ProjFracMin);
      sscanf(pargv[1],"%f",&ProjFracMax);
      sscanf(pargv[2],"%f",&ProjFracDelta);
      ProjFrac = 0.5; // just make it non-zero
      nargsused = 3;
    } else if (!strcmp(option, "--projfrac-max")) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%f",&ProjFracMin);
      sscanf(pargv[1],"%f",&ProjFracMax);
      sscanf(pargv[2],"%f",&ProjFracDelta);
      ProjFrac = 0.5; // just make it non-zero
      GetProjMax = 1;
      nargsused = 3;
    } else if (!strcmp(option, "--projdist")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&ProjFrac);
      ProjFracMin=ProjFrac;
      ProjFracMax=ProjFrac;
      ProjFracDelta=1.0;
      ProjDistFlag = 1;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--cortex")) UseCortexLabel = 1;
    else if (!strcmp(option, "--mask")) {
      if (nargc < 1) argnerr(option,1);
      mask_label_name = pargv[0] ;
      printf("masking output with label %s\n", mask_label_name) ;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--projdist-int") ||
               !strcmp(option, "--projdist-avg")) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%f",&ProjFracMin);
      sscanf(pargv[1],"%f",&ProjFracMax);
      sscanf(pargv[2],"%f",&ProjFracDelta);
      ProjFrac = 0.5; // just make it non-zero
      ProjDistFlag = 1;
      nargsused = 3;
    } else if (!strcmp(option, "--projdist-max")) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%f",&ProjFracMin);
      sscanf(pargv[1],"%f",&ProjFracMax);
      sscanf(pargv[2],"%f",&ProjFracDelta);
      ProjFrac = 0.5; // just make it non-zero
      ProjDistFlag = 1;
      GetProjMax = 1;
      nargsused = 3;
    } else if (!strcmp(option, "--scale")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&scale);
      if (scale == 0) {
        printf("ERROR: scale = 0\n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--thickness")) {
      if (nargc < 1) argnerr(option,1);
      thicknessname = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--interp")) {
      if (nargc < 1) argnerr(option,1);
      interpmethod_string = pargv[0];
      interpmethod = interpolation_code(interpmethod_string);
      if (interpmethod == -1) {
        fprintf(stderr,"ERROR: interpmethod = %s\n",interpmethod_string);
        fprintf(stderr,"  must be either nearest or trilinear \n");
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--float2int")) {
      if (nargc < 1) argnerr(option,1);
      float2int_string = pargv[0];
      float2int = float2int_code(float2int_string);
      if (float2int == -1) {
        fprintf(stderr,"ERROR: float2int = %s\n",float2int_string);
        fprintf(stderr,"  must be either round, floor, or tkreg\n");
        exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcmp(option, "--o") || !strcmp(option, "--out")) {
      if (nargc < 1) argnerr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--use-new")) {
      UseOld = 0;
    } 
    else if (!strcmp(option, "--vsm")) {
      if (nargc < 1) argnerr(option,1);
      vsmfile = pargv[0];
      UseOld = 0;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--out_type") || !strcmp(option, "--ofmt")) {
      if (nargc < 1) argnerr(option,1);
      outtypestring = pargv[0];
      outtype = string_to_type(outtypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--nvox")) {
      if (nargc < 1) argnerr(option,1);
      nvoxfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--rf")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&reshapefactor);
      reshape = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--rft")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&reshapetarget);
      reshape = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--srchits")) {
      if (nargc < 1) argnerr(option,1);
      srchitfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trghits")) {
      if (nargc < 1) argnerr(option,1);
      trghitfile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--fwhm") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&fwhm);
      gstd = fwhm/sqrt(log(256.0));
      nargsused = 1;
    } else if ( !strcmp(option, "--surf-fwhm") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&surf_fwhm);
      surf_gstd = surf_fwhm/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcmp(option, "--srcsynth")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%ld",&seed);
      srcsynth = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--srcsynth-index")) {
      srcsynthindex = 1;
      interpmethod_string = "nearest";
      interpmethod = interpolation_code(interpmethod_string);
    } else if (!strcmp(option, "--seedfile")) {
      if (nargc < 1) argnerr(option,1);
      seedfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--profile")) {
      if(nargc < 6){
	printf("ERROR: --profile requires 6 args\n");
	printf("USAGE: --profile surf vol dist delta sigma output\n");
	exit(1);
      }
      MRIS *surf  = MRISread(pargv[0]);
      if(surf==NULL) exit(1);
      MRI *mri = MRIread(pargv[1]); //norm
      if(mri==NULL) exit(1);
      double dist,delta,sigma;
      sscanf(pargv[2],"%lf",&dist);
      sscanf(pargv[3],"%lf",&delta);
      sscanf(pargv[4],"%lf",&sigma);
      if(delta <= 0) delta = mri->xsize/2.0;
      printf("dist %g, delta=%g, sigma=%g\n",dist,delta,sigma);
      MRI *mri2 = MRISsampleMRINorm(surf, mri, -dist, +dist, delta, sigma, NULL);
      if(mri2==NULL) exit(1);
      MRIwrite(mri2,pargv[5]);
      printf("freeview -f %s:overlay=%s\n",pargv[0],pargv[5]);
      printf("mri_vol2surf --profile done\n");
      exit(0);
    }
    else if (!strcmp(option, "--vol2surf")) {
      // This is an alternative way to run vol2surf that does not rely on
      // the recon-all directory structure. It generates the same as the standard
      // invocation when using --use-new. It does not have all the functionality,
      // eg, avg or max or mapping xyz. 
      int err;
      if(nargc < 9){
	printf("ERROR: --vol2surf requires 9 args\n");
	printf("USAGE: --vol2surf vol surf projtype projdist projmap reg vsm interp output\n");
	printf(" projtype 0=absdist, 1=frac (if 0, projmap arg will be ignored)\n");
	printf(" projmap : map to get the values for projection (usually thickness)\n");
	printf(" reg : LTA registration file (or 'regheader', surf must have vol geom) \n");
	printf("   LTA files can go in either direction if surf has a valid vol geom\n");
	printf(" vsm : voxel shift map for B0 correction (or 'novsm')\n");
	printf(" interp 0=nearest, 1=trilin, 5=cubicbspline\n");
	exit(1);
      }
      MRI *mri = MRIread(pargv[0]);
      if(mri==NULL) exit(1);
      MRIS *surf  = MRISread(pargv[1]);
      if(surf==NULL) exit(1);
      int projtype;
      sscanf(pargv[2],"%d",&projtype);
      double projdist;
      sscanf(pargv[3],"%lf",&projdist);
      if(projtype == 1){
	err = MRISreadCurvatureFile(surf, pargv[4]);
	if(err) exit(1);
      }
      LTA *lta=NULL;
      MATRIX *RegMat=NULL;
      if(strcmp(pargv[5],"regheader") != 0){ // not regheader
	lta = LTAread(pargv[5]);
	if(lta==NULL) exit(1);
	LTAchangeType(lta,REGISTER_DAT);
	VOL_GEOM srcvg;
	getVolGeom(mri, &srcvg);
	vg_isEqual_Threshold = 10e-3;
	if(!vg_isEqual(&srcvg, &(lta->xforms[0].src))){
	  if(!vg_isEqual(&srcvg, &(lta->xforms[0].dst))){
	    printf("ERRRO: input volume VG does not match LTA source or target VG\n");
	    exit(1);
	  }
	  printf("INFO: input volume VG matches LTA target VG, inverting \n");
	  LTA *lta2 = LTAinvert(lta, NULL);
	  lta = lta2;
	}
	RegMat = lta->xforms[0].m_L;
      }
      else {
	if(!surf->vg.valid){
	  printf("ERROR: volume geometry of input surface is not valid, cannot use regheader\n");
	  exit(1);
	}
	MRI *TargVol=NULL;
	TargVol = MRIallocFromVolGeom(&(surf->vg), MRI_UCHAR, 1,1);
	RegMat = MRItkRegMtx(TargVol,mri,NULL);
	MRIfree(&TargVol);
      }
      MRI *vsm = NULL;
      if(strcmp(pargv[6],"novsm") != 0){
	vsm = MRIread(pargv[6]);
	if(vsm==NULL) exit(1);
      }
      int interpmethod=0;
      sscanf(pargv[7],"%d",&interpmethod);
      printf("projtype %d, projdist %g, interp %d\n",projtype,projdist,interpmethod);
      MRI *sval = MRIvol2surfVSM(mri, RegMat, surf, vsm, interpmethod, NULL, projdist, projtype, 1,NULL);
      err = MRIwrite(sval,pargv[8]);
      printf("mri_vol2surf --volsurf done\n");
      exit(err);
      // done with --vol2surf
    }
    else if (!strcmp(option, "--norm-pointset")) {
      if(nargc < 5){
	printf("ERROR: --norm-pointset requires 5 args\n");
	printf("USAGE: --norm-pointset surf vtxno dist delta output\n");
	exit(1);
      }
      MRIS *surf  = MRISread(pargv[0]);
      if(surf==NULL) exit(1);
      int vtxno;
      double dist,delta;
      sscanf(pargv[1],"%d",&vtxno);
      sscanf(pargv[2],"%lf",&dist);
      sscanf(pargv[3],"%lf",&delta);
      FILE *fp = fopen(pargv[4],"w");
      if(fp==NULL) {
	printf("ERROR: opening %s for writing\n",pargv[4]);
	exit(1);
      }
      printf("vtxno = %d, dist %g, delta=%g\n",vtxno,dist,delta);
      MRISnorm2Pointset(surf, vtxno, -dist, dist, delta, fp);
      fclose(fp);
      printf("mri_vol2surf --norm-pointset done\n");
      exit(0);
    }
    else {
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
  printf("   --mov input volume path (or --src)\n");
  printf("   --ref reference volume name (default=orig.mgz\n");
  printf("   --reg source registration (can be reg.dat or lta) \n");
  printf("   --regheader subject\n");
  printf("   --mni152reg : $FREESURFER_HOME/average/mni152.register.dat\n");
  printf("   --rot   Ax Ay Az : rotation angles (deg) to apply to reg matrix\n");
  printf("   --trans Tx Ty Tz : translation (mm) to apply to reg matrix\n");
  printf("   --float2int float-to-int conversion method "
         "(<round>, tkregister )\n");
  printf("   --fixtkreg : make make registration matrix round-compatible\n");
  printf("   --fwhm fwhm : smooth input volume (mm)\n");
  printf("   --surf-fwhm fwhm : smooth output surface (mm)\n");
  printf("\n");
  printf("   --trgsubject target subject (if different than reg)\n");
  printf("   --hemi hemisphere (lh or rh) \n");
  printf("   --surf target surface (default = white) DO NOT USE 'inflated' \n");
  printf("      If you want to display on the inflated, sample it on \n");
  printf("      the white surface, then display it on any surface, including inflated\n");
  printf("   --srcsubject source subject (override that in reg)\n");
  printf("\n");
  printf(" Options for use with --trgsubject\n");
  printf("   --surfreg    surface registration (sphere.reg)  \n");
  printf("   --icoorder   order of icosahedron when trgsubject=ico\n");
  //printf("   --nohash flag to keep the hash table from being used. \n");
  printf("\n");
  printf(" Options for projecting along the surface normal:\n");
  printf("   --projfrac frac : (0->1)fractional projection along normal \n");
  printf("   --projfrac-avg min max del : average along normal\n");
  printf("   --projfrac-max min max del : max along normal\n");
  printf("   --projdist mmdist : distance projection along normal \n");
  printf("   --projdist-avg min max del : average along normal\n");
  printf("   --projopt <fraction stem> : use optimal linear estimation and previously\n"
         "computed volume fractions (see mri_compute_volume_fractions)\n");
  printf("   --projdist-max min max del : max along normal\n");
  printf("   --mask label : mask the output with the given label file (usually cortex)\n");
  printf("   --cortex : use hemi.cortex.label from trgsubject\n");
  
  //printf("   --thickness thickness file (thickness)\n");
  printf("\n");
  printf(" Options for output\n");
  printf("   --o         output path\n");
  printf("   --out_type  output format\n");
  printf("   --frame   nth :  save only 0-based nth frame \n");
  printf("   --noreshape do not save output as multiple 'slices'\n");
  printf("   --rf R  integer reshaping factor, save as R 'slices'\n");
  printf("   --srchit   volume to store the number of hits at each vox \n");
  printf("   --srchit_type  source hit volume format \n");
  printf("   --nvox nvoxfile : write number of voxels intersecting surface\n");
  printf("\n");
  printf(" Other Options\n");
  printf("   --reshape : so dims fit in nifti or analyze\n");
  printf("   --noreshape : do not reshape (default)\n");
  printf("   --reshape3d : reshape fsaverage (ico7) into 42 x 47 x 83\n");
  printf("   --scale scale : multiply all intensities by scale.\n");
  printf("   --v vertex no : debug mapping of vertex.\n");
  printf("   --srcsynth seed : synthesize source volume\n");
  printf("   --srcsynth-index : synthesize source volume with volume index no\n");
  printf("   --seedfile fname : save synth seed to fname\n");
  printf("   --sd SUBJECTS_DIR \n");
  printf("   --profile surf vol dist delta sigma output\n");
  printf("     Computes intensity profile from -dist:delta:+dist\n");
  printf("     If delta is <= 0, then xsize/2 is used\n");
  printf("     If sigma >= 0, then the gradient is estimated with smoothing parameter sigma\n");
  printf("   --norm-pointset surf vtxno dist delta output\n");
  printf("     Creates a freeview pointset using points along the normal\n");
  printf("\n");
  printf("\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("   --interp    interpolation method (<nearest> or trilinear)\n");
  printf("   --vg-thresh thrshold : threshold for  'ERROR: LTAconcat(): LTAs 0 and 1 do not match'\n");
  printf("\n");
  printf("   --vol2surf vol surf projtype projdist projmap reg vsm interp output\n");
  printf("    This is an alternative way to run vol2surf that does not rely on the recon-all\n");
  printf("    directory structure. Generates the same output as standard invocation as long \n");
  printf("    as --use-new is added. --vol2surf does not have all the functionality. \n");
  printf("      projtype 0=absdist, 1=frac (if 0, projmap arg will be ignored)\n");
  printf("      projmap : map to get the values for projection (usually thickness)\n");
  printf("      reg : LTA registration file (or 'regheader', surf must have vol geom) \n");
  printf("        LTA files can go in either direction if surf has a valid vol geom\n");
  printf("      vsm : voxel shift map for B0 correction (or 'novsm')\n");
  printf("      interp 0=nearest, 1=trilin, 5=cubicbspline\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
  //printf("   --src_type  input volume format \n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf(
    "This program will resample a volume onto a surface of a subject or the \n"
    "sphere. The output can be viewed on the surface "
    "(using tksurfer) or can \n"
    "be using for surface-based intersubject averaging. "
    "This program supersedes\n"
    "paint.\n"
    "\n"
    "OPTIONS\n"
    "\n"
    "  --mov path to input volume (see below). Can also use --src\n"
    "\n"
    "  --reg file : registration file as computed by tkregister,\n"
    "    spmregister, bbregister, etc. This file\n"
    "    can be an LTA (better) or a register.dat, in which case it has the following format:\n"
    "\n"
    "        subjectname\n"
    "        in-plane resolution(mm)\n"
    "        between-plane resolution(mm)\n"
    "        intensity\n"
    "        m11 m12 m13 x0\n"
    "        m21 m22 m33 y0\n"
    "        m31 m22 m33 z0\n"
    "          0   0   0  1\n"
    "        <float2int_method>\n"
    "\n"
    "    where subjectname is the name of the subject as found in \n"
    "    $SUBJECTS_DIR, in-plane resolution is the distance between\n"
    "    adjacent rows or columes, between-plane resolution is the \n"
    "    distance between adjacent slices, intensity is ignored, and \n"
    "    the remainder is the 4x4 affine tranform that converts XYZ in\n"
    "    the input volume to XYZ in the subject's anatomical space.\n"
    "    The volume is mapped onto the surface of subject subjectname\n"
    "    (unless --trgsubject is specified). There may or may not be \n"
    "    another line with the method used to convert voxel indices\n"
    "    from floating to integer during registration. If tkregiser was\n"
    "    used to register the volume, the method should be "
    "blank or 'tkregister'\n"
    "    (no quotes), otherwise it should be 'round'. "
    "This can be overridden \n"
    "    with --float2int.\n"
    "\n"
    "  --regheader subject \n"
    "\n"
    "    Compute registration from header information, ie, assume the\n"
    "    volume geometry between the mov volume and the subject/mri/orig.mgz\n"
    "    align the two volumes. This is the same as in tkregister2.\n"
    "\n"
    "  --float2int method: override float2int method in registration file.\n"
    "    See BUGS.\n"
    "\n"
    "  --fixtkreg\n"
    "\n"
    "    Attempt to convert the registration matrix so that it is round \n"
    "    (or nearest neighbor) compatible. Setting this flag will only have \n"
    "    an effect if the float2int method is tkregister. It will 'fix' the \n"
    "    matrix and change the float2int method to round. "
    "Don't use this flag \n"
    "    unless you know what you are doing. \n"
    "\n"
    "  --fwhm FWHM\n"
    "\n"
    "    Smooth input volume with a gaussian kernal with FWHM mm.\n"
    "\n"
    "  --trgsubject target subject name : resample volume onto this subject \n"
    "    instead of the one found in the registration file. "
    "The target subject \n"
    "    can be either a subject name (as found in $SUBJECTS_DIR) or ico \n"
    "    (to map onto the sphere). If the target subject is not the source \n"
    "    subject, then the surfaces are mapped using each"
    " subject's spherical\n"
    "    surface registration (?h.sphere.reg or that "
    "specified with --surfreg).\n"
    "    If the target subject is ico, then the volume is resampled onto an\n"
    "    icosahedron, which is used to uniformly sample a sphere. "
    "This requires\n"
    "    specifying the icosahedron order (see --icoorder).\n"
    "\n"
    "\n"
    "  --hemi hemisphere : lh = left hemisphere, rh = right hemisphere\n"
    "\n"
    "  --surf surfacename : the surface on which to resample. The default is\n"
    "  white. It will look for "
    "$SUBJECTS_DIR/subjectname/surf/?h.surfacename\n"
    "  DO NOT specify 'inflated'. If you want to display on the inflated surface,\n"
    "  or any other surface, specify 'white' here then display it on any surface,\n"
    "  including the inflated surface.\n"
    "\n"
    "  --surfreg intersubject registration surface : default (sphere.reg).\n"
    "    This is a representation of a subject's cortical surface after it\n"
    "    has been registered/morphed with a template spherical surface.\n"
    "\n"
    "  --icoorder icosahedron order number: this specifies the size of the\n"
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
    "  --projfrac fraction : fraction (0,1) of the cortical thickness \n"
    "    at each vertex to project along the surface normal. Default 0. \n"
    "    When set at 0.5 with the white surface, this should sample in the\n"
    "    middle of the cortical surface. "
    "This requires that a ?h.thickness file \n"
    "    exist for the source subject. Note, the faction can be less than\n"
    "    zero, in which case it will project into the white matter. See also\n"
    "    --projdist.\n"
    "\n"
    "  --projdist mmdist\n"
    "\n"
    "    Same as --projfrac but projects the given distance in mm at all\n"
    "    points of the surface regardless of thickness.\n"
    "\n"
    "  --projfrac-avg min max delta\n"
    "  --projdist-avg min max delta\n"
    "\n"
    "    Same idea as --projfrac and --projdist, but sample at each of the points\n"
    "    between min and max at a spacing of delta. The samples are then averaged\n"
    "    together. The idea here is to average along the normal.\n"
    "\n"
    "  --o output path : location to store the data (see below)\n"
    "  --out_type format of output (see below)\n"
    "\n"
    "  --frame 0-based frame number : sample and save only the given frame \n"
    "    from the source volume (needed when out_type = paint). Default 0.\n"
    "\n"
    "  --reshape : save the output as multiple 'slices'. This is \n"
    "    for logistical purposes (eg, in analyze and nifti formats\n"
    "    the size of a dimension cannot exceed 2^15). This used to be\n"
    "    the default behavior. This has no effect when the output type\n"
    "    is paint.\n"
    "\n"
    "  --rf R\n"
    "\n"
    "    Explicity set the reshaping factor to R. R must be an integer factor \n"
    "    of the number of vertices. \n"
    "\n"
    "  --srchit volid\n"
    "\n"
    "    Save a volume (the same size as the source) in which the value at\n"
    "    each voxel is the number of surface vertices it maps to. The number of\n"
    "    voxels hit at least once is printed to stdout as :\n"
    "       'Number of source voxels hit' \n"
    "\n"
    "  --nvox nvoxfile\n"
    "\n"
    "    Save the number of voxels intersecting the surface in the file\n"
    "    nvoxfile.\n"
    "\n"
    "  --version : print version and exit.\n"
    "\n"
    "\n"
    "  --srcsynth seed\n"
    "\n"
    "    Synthesize the source volume with white gaussian noise. seed\n"
    "    is the seed to the random generator. Use -1 to have the seed\n"
    "    automatically chosen based on time-of-day. See also --seedfile\n"
    "\n"
    "  --seedfile fname\n"
    "\n"
    "    Save random number generator seed in fname. Has no effect without\n"
    "    --srcsynth. This is only useful for keeping track of the distribution\n"
    "    of seeds for debugging purposes.\n"
    "\n"
    "\n"
    "SPECIFYING THE INPUT/OUTPUT PATH and TYPE\n"
    "\n"
    "  mri_vol2surf accepts all input/output types as mri_convert (see \n"
    "  mri_convert --help for more information). In addition, an output \n"
    "  type of 'paint' can be specified. This outputs data in a form\n"
    "  that can be easily read by tksurfer (also known as a '.w file').\n"
    "  See BUGS for more information on paint output.\n"
    "\n"
    "NOTES\n"
    "\n"
    "  The output will be a data set with Nv/R colums, 1 row, R slices, and Nf frames,\n"
    "  where Nv is the number of verticies in the output surface, and Nf is the \n"
    "  number of frames in the input volume (unless the output format is paint, in\n"
    "  which case only one frame is written out). R is the reshaping factor. R is 6 \n"
    "  for the icosaheron. For non-ico, the prime factor of Nv closest to 6 is chosen. \n"
    "  Reshaping can be important for logistical reasons (eg, Nv can easily exceed \n"
    "  the maximum number of elements allowed in the analyze format). R can be forced \n"
    "  to 1 with --noreshape. The user can also explicity set R from the command-line\n"
    "  using --rf. Any geometry information saved with the output file will \n"
    "  be bogus.\n"
    "\n"
    "  When resampling for fixed-effects intersubject averaging, make sure\n"
    "  to resample variance and not standard deviation. This is automatically\n"
    "  accomplished when the input volume has been produced by the FS-FAST \n"
    "  selxavg program.\n"
    "\n"
    "EXAMPLES:\n"
    "\n"
    "  1. To paint the third frame of bfloat volume sig registered with tkregister\n"
    "      onto a the left hemisphere of a surface\n"
    "\n"
    "     mri_vol2surf --mov sig.bhdr --reg register.dat \n"
    "       --hemi lh --o ./sig-lh.mgh \n"
    "\n"
    "     This will create sig-lh.mgh in the current directory, which can then\n"
    "     be loaded into tksurfer\n"
    "\n"
    "  2. To convert an analyze volume onto fsaverage (right hemisphere)\n"
    "\n"
    "     mri_vol2surf --mov sig.img  --reg register.dat \n"
    "       --hemi rh --o ./sig-rh.img --trgsubject fsaverage --icoorder 7\n"
    "\n"
    "BUG REPORTS: send bugs to analysis-bugs@nmr.mgh.harvard.edu. Make sure \n"
    "  to include the version and full command-line.\n"
    "\n"
    "BUGS:\n"
    "\n"
    "  When the output format is paint, the output file must be specified with\n"
    "  a partial path (eg, ./data-lh.w) or else the output will be written into\n"
    "  the subject's anatomical directory.\n"
    "\n"
    "  Currently no support for searching along the surface normal for a maximum\n"
    "  value (as can be done with the paint program)\n"
    "\n"
    "  The ability to put the float2int conversion method in the registration file\n"
    "  is new as of Fall 2001. Consequently, there may be some registration files\n"
    "  that do not have a method string, which will force the method to be that\n"
    "  of tkregister. This is can be overridden with --float2int.\n"
    "\n"
    "\n"
    "AUTHOR: Douglas N. Greve, Ph.D., MGH-NMR Center (greve@nmr.mgh.harvard.edu)\n"
    "\n"
  ) ;

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
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
/* --------------------------------------------- */
static void check_options(void) {
  struct timeval tv;
  FILE *fp;

  if (srcvolid == NULL) {
    fprintf(stderr,"A source volume path must be supplied\n");
    exit(1);
  }
  if (interpmethod_string == NULL) interpmethod_string = "nearest";
  interpmethod = interpolation_code(interpmethod_string);
  if (interpmethod == -1) {
    fprintf(stderr,"ERROR: interpmethod = %s\n",interpmethod_string);
    fprintf(stderr,"  must be either nearest or trilinear\n");
    exit(1);
  }

  if (srcregfile == NULL && !regheader) {
    printf("ERROR: must specify a source registration file or --regheader\n");
    exit(1);
  }

  if (srctype == MRI_VOLUME_TYPE_UNKNOWN) {
    if (defaulttype == MRI_VOLUME_TYPE_UNKNOWN)
      srctype = mri_identify(srcvolid);
    else
      srctype = defaulttype;
  }
  if (srctype == MRI_VOLUME_TYPE_UNKNOWN) {
    fprintf(stderr,"ERROR: could not determine type of %s\n",srcvolid);
    exit(1);
  }

  if (srchitvolid != NULL) {
    if (srchittype == MRI_VOLUME_TYPE_UNKNOWN) {
      if (defaulttype == MRI_VOLUME_TYPE_UNKNOWN)
        srchittype = mri_identify(srchitvolid);
      else
        srchittype = defaulttype;
    }
    if (srchittype == MRI_VOLUME_TYPE_UNKNOWN) {
      printf("ERROR: could not determine type of %s\n",srcvolid);
      exit(1);
    }
  }

  if (outfile==NULL) {
    printf("ERROR: no output file specified\n");
    exit(1);
  }

  if (outtypestring != NULL &&
      (!strcasecmp(outtypestring,"w") || !strcasecmp(outtypestring,"paint"))) {
    printf("INFO: output format is paint\n");
  } else {
    if (outtype == MRI_VOLUME_TYPE_UNKNOWN) {
      if (defaulttype == MRI_VOLUME_TYPE_UNKNOWN) {
        outtype = mri_identify(outfile);
        if (outtype == MRI_VOLUME_TYPE_UNKNOWN) {
          fprintf(stderr,"ERROR: could not determine type of %s\n",outfile);
          exit(1);
        }
      } else outtype = defaulttype;
    }
  }

  if (hemi == NULL) {
    fprintf(stderr,"ERROR: no hemifield specified\n");
    exit(1);
  }

  if (trgsubject != NULL && strcmp(trgsubject,"ico") && IcoOrder > -1) {
    fprintf(stderr,"ERROR: --icoorder can only be used with "
            "--trgsubject ico\n");
    exit(1);
  }

  if (trgsubject != NULL && !strcmp(trgsubject,"ico") && IcoOrder < 0) {
    fprintf(stderr,"ERROR: need to specify --icoorder with "
            "--trgsubject ico\n");
    exit(1);
  }

  if (seed < 0) {
    gettimeofday(&tv, NULL);
    seed = tv.tv_sec + tv.tv_usec;
  }
  if (seedfile != NULL) {
    fp = fopen(seedfile,"w");
    if (fp == NULL) {
      printf("ERROR: cannot open seed file %s\n",seedfile);
      exit(1);
    }
    fprintf(fp,"%ld\n",seed);
    fclose(fp);
  }

  // paint format but framesave has not been set
  if(framesave < 0 && outtypestring != NULL &&
      (!strcasecmp(outtypestring,"w") || !strcasecmp(outtypestring,"paint"))) 
    framesave = 0;

  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"srcvol = %s\n",srcvolid);
  if (srctypestring != NULL) fprintf(fp,"srctype = %s\n",srctypestring);
  if (srcregfile != NULL) fprintf(fp,"srcreg = %s\n",srcregfile);
  else                  fprintf(fp,"srcreg unspecified\n");
  fprintf(fp,"srcregold = %d\n",srcoldreg);
  if (srcwarp != NULL) fprintf(fp,"srcwarp = %s\n",srcwarp);
  else                   fprintf(fp,"srcwarp unspecified\n");

  if (srchitvolid != NULL) {
    fprintf(fp,"srchitvol = %s\n",srchitvolid);
    if (srchittypestring != NULL)
      fprintf(fp,"srchittype = %s\n",srchittypestring);
  }

  fprintf(fp,"surf = %s\n",surfname);
  fprintf(fp,"hemi = %s\n",hemi);
  if (trgsubject != NULL) {
    fprintf(fp,"trgsubject = %s\n",trgsubject);
    fprintf(fp,"surfreg = %s\n",surfreg);
  }
  if (ProjFrac != 0.0) {
    if (!ProjDistFlag) {
      fprintf(fp,"ProjFrac = %g\n",ProjFrac);
      fprintf(fp,"thickness = %s\n",thicknessname);
    } else {
      fprintf(fp,"ProjDist = %g\n",ProjFrac);
    }
  }

  fprintf(fp,"reshape = %d\n",reshape);
  fprintf(fp,"interp = %s\n",interpmethod_string);
  fprintf(fp,"float2int = %s\n",float2int_string);
  fprintf(fp,"GetProjMax = %d\n",GetProjMax);

  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}



/*---------------------------------------------------------------*/
MRI *MRIvol2surf(MRI *SrcVol, MATRIX *Rtk, MRI_SURFACE *TrgSurf, 
		 MRI *vsm, int InterpMethod, MRI *SrcHitVol, 
		 float ProjFrac, int ProjType, int nskip)
{
  MATRIX *ras2vox, *Scrs, *Txyz;
  MRI *TrgVol;
  int   irow, icol, islc; /* integer row, col, slc in source */
  float frow, fcol, fslc; /* float row, col, slc in source */
  float srcval, *valvect, rshift;
  int frm, vtx,nhits, err;
  double rval;
  float Tx, Ty, Tz;

  if(vsm){
    err = MRIdimMismatch(vsm,SrcVol,0);
    if(err){
      printf("ERROR: MRIvol2surf: vsm dimension mismatch %d\n",err);
      exit(1);
    }
  }

  vox2ras = MRIxfmCRS2XYZtkreg(SrcVol);
  ras2vox = MatrixInverse(vox2ras,NULL);
  if(Rtk != NULL) ras2vox = MatrixMultiply(ras2vox,Rtk,NULL);
  MatrixFree(&vox2ras);
  // ras2vox now converts surfacs RAS to SrcVol vox

  /* preallocate the row-col-slc vectors */
  Scrs = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz->rptr[3+1][0+1] = 1.0;

  /* allocate a "volume" to hold the output */
  TrgVol = MRIallocSequence(TrgSurf->nvertices,1,1,MRI_FLOAT,SrcVol->nframes);
  if (TrgVol == NULL) return(NULL);
  MRIcopyHeader(SrcVol,TrgVol);

  /* Zero the source hit volume */
  if(SrcHitVol != NULL){
    MRIconst(SrcHitVol->width,SrcHitVol->height,SrcHitVol->depth,
             1,0,SrcHitVol);
  }

  srcval = 0;
  valvect = (float *) calloc(sizeof(float),SrcVol->nframes);
  nhits = 0;
  /*--- loop through each vertex ---*/
  for(vtx = 0; vtx < TrgSurf->nvertices; vtx+=nskip){

    if(ProjFrac != 0.0){
      if(ProjType == 0)
        ProjNormDist(&Tx,&Ty,&Tz,TrgSurf,vtx,ProjFrac);
      else
        ProjNormFracThick(&Tx,&Ty,&Tz,TrgSurf,vtx,ProjFrac);
    }
    else{
      Tx = TrgSurf->vertices[vtx].x;
      Ty = TrgSurf->vertices[vtx].y;
      Tz = TrgSurf->vertices[vtx].z;
    }

    /* Load the Target xyz vector */
    Txyz->rptr[0+1][0+1] = Tx;
    Txyz->rptr[1+1][0+1] = Ty;
    Txyz->rptr[2+1][0+1] = Tz;

    /* Compute the corresponding Source col-row-slc vector */
    Scrs = MatrixMultiply(ras2vox,Txyz,Scrs);
    fcol = Scrs->rptr[1][1];
    frow = Scrs->rptr[2][1];
    fslc = Scrs->rptr[3][1];

    icol = nint(fcol);
    irow = nint(frow);
    islc = nint(fslc);

    /* check that the point is in the bounds of the volume */
    if (irow < 0 || irow >= SrcVol->height ||
        icol < 0 || icol >= SrcVol->width  ||
        islc < 0 || islc >= SrcVol->depth ) continue;

    if(vsm){
      MRIsampleSeqVolume(vsm, fcol, frow, fslc, &rshift, 0,0);
      frow += rshift;
      irow = nint(frow);
      if(irow < 0 || irow >= SrcVol->height) continue;
    }

    if (Gdiag_no == vtx){
      printf("diag -----------------------------\n");
      printf("vtx = %d  %g %g %g\n",vtx,Txyz->rptr[1][1],
	     Txyz->rptr[2][1],Txyz->rptr[3][1]);
      printf("fCRS  %g %g %g\n",Scrs->rptr[1][1],
             Scrs->rptr[2][1],Scrs->rptr[3][1]);
      printf("CRS  %d %d %d\n",icol,irow,islc);
    }

    /* only gets here if it is in bounds */
    nhits ++;

    /* Assign output volume values */
    if(InterpMethod == SAMPLE_TRILINEAR) {
      MRIsampleSeqVolume(SrcVol, fcol, frow, fslc,
                         valvect, 0, SrcVol->nframes-1) ;
      if(Gdiag_no == vtx) printf("val = %f\n", valvect[0]) ;
      for (frm = 0; frm < SrcVol->nframes; frm++)
        MRIFseq_vox(TrgVol,vtx,0,0,frm) = valvect[frm];
    }
    else {
      for (frm = 0; frm < SrcVol->nframes; frm++)  {
        switch (InterpMethod) {
        case SAMPLE_NEAREST:
          srcval = MRIgetVoxVal(SrcVol,icol,irow,islc,frm);
          break ;
        case SAMPLE_SINC:      /* no multi-frame */
          MRIsincSampleVolume(SrcVol, fcol, frow, fslc, 5, &rval) ;
          srcval = rval;
          break ;
        } //switch
        MRIFseq_vox(TrgVol,vtx,0,0,frm) = srcval;
        if(Gdiag_no == vtx) printf("val[%d] = %f\n", frm, srcval) ;
      } // for
    }// else
    if(SrcHitVol != NULL) MRIFseq_vox(SrcHitVol,icol,irow,islc,0)++;
  }

  MatrixFree(&ras2vox);
  MatrixFree(&Scrs);
  MatrixFree(&Txyz);
  free(valvect);

  //printf("vol2surf_linear: nhits = %d/%d\n",nhits,TrgSurf->nvertices);

  return(TrgVol);
}

#if 0
MRI *
estimate_gm_values(MRI *mri_wm, MRI *mri_gm, MRI *mri_csf,
        MRI *SrcVol,MATRIX *Qsrc, MATRIX *Fsrc, MATRIX *Wsrc, MATRIX *Dsrc,
        MRI_SURFACE *TrgSurf, int InterpMethod, int float2int, MRI *SrcHitVol)
{
  MATRIX *QFWDsrc;
  MATRIX *Scrs, *Txyz;
  MRI *TrgVol, *mri_samples, *mri_voxels ;
  int   irow_src, icol_src, islc_src; /* integer row, col, slc in source */
  float frow_src, fcol_src, fslc_src; /* float row, col, slc in source */
  float srcval, Tx, Ty, Tz ;
  int frm, FreeQsrc=0;
  int vtx,nhits;
  float *valvect;
  double rval;
  int  ProjDistFlag = 1, nskip = 1 ;
  float ProjFrac = 0.5 ;

  if(Qsrc == NULL){
    Qsrc = MRIxfmCRS2XYZtkreg(SrcVol);
    Qsrc = MatrixInverse(Qsrc,Qsrc);
    FreeQsrc = 1;
  }

  /* compute the transforms */
  QFWDsrc = ComputeQFWD(Qsrc,Fsrc,Wsrc,Dsrc,NULL);
  if (Gdiag_no >= 0)
  {
    printf("QFWDsrc: vol2surf: ------------------------------\n");
    MatrixPrint(stdout,QFWDsrc);
    printf("--------------------------------------------------\n");
  }

  /* preallocate the row-col-slc vectors */
  Scrs = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz->rptr[3+1][0+1] = 1.0;

  /* allocate a "volume" to hold the output */
  TrgVol = MRIallocSequence(TrgSurf->nvertices,1,1,MRI_FLOAT,SrcVol->nframes);
  if (TrgVol == NULL) return(NULL);
  MRIcopyHeader(SrcVol,TrgVol);

  /* Zero the source hit volume */
  if (SrcHitVol != NULL)
  {
    MRIconst(SrcHitVol->width,SrcHitVol->height,SrcHitVol->depth,
             1,0,SrcHitVol);
  }

  srcval = 0;
  valvect = (float *) calloc(sizeof(float),SrcVol->nframes);
  nhits = 0;
  /*--- loop through each vertex ---*/
  mri_samples = build_sample_array(TrgSurf, SrcVol, QFWDsrc, 1.0, 0.1, 10,
                                   mri_wm, mri_gm, mri_csf,
                                   &mri_voxels);

  for (vtx = 0; vtx < TrgSurf->nvertices; vtx+=nskip)
  {
    if (ProjFrac != 0.0)
      if (ProjDistFlag)
        ProjNormDist(&Tx,&Ty,&Tz,TrgSurf,vtx,ProjFrac);
      else
        ProjNormFracThick(&Tx,&Ty,&Tz,TrgSurf,vtx,ProjFrac);
    else
    {
      Tx = TrgSurf->vertices[vtx].x;
      Ty = TrgSurf->vertices[vtx].y;
      Tz = TrgSurf->vertices[vtx].z;
    }

    /* Load the Target xyz vector */
    Txyz->rptr[0+1][0+1] = Tx;
    Txyz->rptr[1+1][0+1] = Ty;
    Txyz->rptr[2+1][0+1] = Tz;

    /* Compute the corresponding Source col-row-slc vector */
    MatrixMultiply(QFWDsrc,Txyz,Scrs);
    fcol_src = Scrs->rptr[1][1];
    frow_src = Scrs->rptr[2][1];
    fslc_src = Scrs->rptr[3][1];

    /* nearest neighbor */
    switch (float2int)
    {
    case FLT2INT_ROUND:
      icol_src = nint(fcol_src);
      irow_src = nint(frow_src);
      islc_src = nint(fslc_src);
      break;
    case FLT2INT_FLOOR:
      icol_src = (int)floor(fcol_src);
      irow_src = (int)floor(frow_src);
      islc_src = (int)floor(fslc_src);
      break;
    case FLT2INT_TKREG:
      icol_src = (int)floor(fcol_src);
      irow_src = (int) ceil(frow_src);
      islc_src = (int)floor(fslc_src);
      break;
    default:
      fprintf(stderr,"vol2surf_linear(): unrecoginized float2int code %d\n",
              float2int);
      MRIfree(&TrgVol);
      return(NULL);
      break;
    }

    /* check that the point is in the bounds of the volume */
    if (irow_src < 0 || irow_src >= SrcVol->height ||
        icol_src < 0 || icol_src >= SrcVol->width  ||
        islc_src < 0 || islc_src >= SrcVol->depth ) continue;

    if (Gdiag_no == vtx)
    {
      printf("diag -----------------------------\n");
      printf("vtx = %d  %g %g %g\n",vtx,Tx,Ty,Tz);
      printf("fCRS  %g %g %g\n",Scrs->rptr[1][1],
             Scrs->rptr[2][1],Scrs->rptr[3][1]);
      printf("CRS  %d %d %d\n",icol_src,irow_src,islc_src);
    }


    /* only gets here if it is in bounds */
    nhits ++;

    /* Assign output volume values */
    if (InterpMethod == SAMPLE_TRILINEAR)
    {
      MRIsampleSeqVolume(SrcVol, fcol_src, frow_src, fslc_src,
                         valvect, 0, SrcVol->nframes-1) ;
      if (Gdiag_no == vtx)
        printf("val = %f\n", valvect[0]) ;
      for (frm = 0; frm < SrcVol->nframes; frm++)
        MRIFseq_vox(TrgVol,vtx,0,0,frm) = valvect[frm];
    }
    else
    {
      for (frm = 0; frm < SrcVol->nframes; frm++)
      {
        switch (InterpMethod)
        {
        case SAMPLE_NEAREST:
          srcval = MRIgetVoxVal(SrcVol,icol_src,irow_src,islc_src,frm);
          break ;
        case SAMPLE_SINC:      /* no multi-frame */
          MRIsincSampleVolume(SrcVol, fcol_src, frow_src, fslc_src, 5, &rval) ;
          srcval = rval;
          break ;
        } //switch
        MRIFseq_vox(TrgVol,vtx,0,0,frm) = srcval;
        if (Gdiag_no == vtx)
          printf("val[%d] = %f\n", frm, srcval) ;
      } // for
    }// else
    if (SrcHitVol != NULL)
      MRIFseq_vox(SrcHitVol,icol_src,irow_src,islc_src,0)++;
  }

  MatrixFree(&QFWDsrc);
  MatrixFree(&Scrs);
  MatrixFree(&Txyz);
  free(valvect);
  if(FreeQsrc) MatrixFree(&Qsrc);

  //printf("vol2surf_linear: nhits = %d/%d\n",nhits,TrgSurf->nvertices);

  return(TrgVol);
}
#endif


#define MAX_SAMPLES 1000
#define FOUND_WM  0x01
#define FOUND_GM  0x02
#define FOUND_CSF 0x04

#define MAX_NBRS 10000
MRI *
build_sample_array(MRI_SURFACE *mris, MRI *mri_src, MATRIX *m, 
                   float din, float dout, int nsamples, 
                   MRI *mri_wm, MRI *mri_gm, MRI *mri_csf)
{
  float  fcol_src, frow_src, fslc_src, Tx, Ty, Tz,
         dist, thick, sample_dist ;
  int    vno, vno_nbr, srcval, icol_src, irow_src, islc_src, n, nfound,  index,
    nbr, vnum,nsize, found, done, failed, vlist[MAX_NBRS],
    min_needed = 3*2 ; // 3 parameters estimated with 2x overdetermination
  MATRIX *Scrs, *Txyz, *m_A, *m_p, *m_S, *m_inv;
  MRI    *mri_samples, *mri_voxels_c, *mri_voxels_r, *mri_voxels_s, *mri_sampled ;
  double rms_mean, rms_var ;

  mri_samples = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;
  Scrs = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz->rptr[3+1][0+1] = 1.0;

  mri_sampled = MRIcloneDifferentType(mri_src, MRI_UCHAR);
  mri_voxels_c = MRIallocSequence(mris->nvertices, 1,1,MRI_INT,nsamples) ;
  if (mri_voxels_c == NULL)
    ErrorExit(ERROR_NOMEMORY, 
              "build_sample_array: couldn't allocate %d x %d array",
              nsamples, mris->nvertices)  ;

  mri_voxels_r = MRIallocSequence(mris->nvertices, 1,1,MRI_INT,nsamples) ;
  if (mri_voxels_r == NULL)
    ErrorExit(ERROR_NOMEMORY, 
              "build_sample_array: couldn't allocate %d x %d array",
              nsamples, mris->nvertices)  ;
  mri_voxels_s = MRIallocSequence(mris->nvertices, 1,1,MRI_INT,nsamples) ;
  if (mri_voxels_s == NULL)
    ErrorExit(ERROR_NOMEMORY, 
              "build_sample_array: couldn't allocate %d x %d array",
              nsamples, mris->nvertices)  ;

  m_A = MatrixAlloc(MAX_SAMPLES, 3, MATRIX_REAL) ; // vfract matrix
  m_S = MatrixAlloc(MAX_SAMPLES, 1, MATRIX_REAL) ; // observed signals
  m_p = MatrixAlloc(3, 1, MATRIX_REAL) ;  // parameters to be computed

  rms_mean = rms_var = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX const * const v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    thick = v->curv ;
    dist = -din ;
    sample_dist = (din + thick + dout) / (nsamples-1) ;
    for (n = 0 ; n < nsamples ; n++, dist += sample_dist)
    {
      /* Load the Target xyz vector */
      ProjNormDist(&Tx,&Ty,&Tz, mris, vno, dist);
      Txyz->rptr[0+1][0+1] = Tx;
      Txyz->rptr[1+1][0+1] = Ty;
      Txyz->rptr[2+1][0+1] = Tz;
      
      /* Compute the corresponding Source col-row-slc vector */
      MatrixMultiply(m,Txyz,Scrs);
      fcol_src = Scrs->rptr[1][1];
      frow_src = Scrs->rptr[2][1];
      fslc_src = Scrs->rptr[3][1];
      icol_src = nint(fcol_src) ;
      irow_src = nint(frow_src) ;
      islc_src = nint(fslc_src) ;
      MRIsetVoxVal(mri_voxels_c, vno, 0, 0, n, icol_src) ;
      MRIsetVoxVal(mri_voxels_r, vno, 0, 0, n, irow_src) ;
      MRIsetVoxVal(mri_voxels_s, vno, 0, 0, n, islc_src) ;
    }
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    nfound = index = 0 ; nsize = 1 ; failed = 0 ;
    found = 0 ;  // bitask to make sure we find at least some of each tissue class
    // find enough samples to do optimal linear estimation
    do
    {
      switch (nsize)
      {
      case 1: vnum = vt->vnum ; break ;
      case 2: vnum = vt->v2num ; break ;
      default: 
      case 3: vnum = vt->v3num ; break ;
      }        
      memmove(vlist+1, vt->v, vnum*sizeof(vt->v[0])) ;
      vlist[0] = vno ; vnum++ ;
      for (nbr = 0 ; nbr < vnum ; nbr++)
      {
        float gm, wm, csf ;

        vno_nbr = vlist[nbr] ;
        for (n = 0 ; n < nsamples ; n++)
        {
          icol_src = (int)MRIgetVoxVal(mri_voxels_c, vno_nbr, 0, 0, n);
          irow_src = (int)MRIgetVoxVal(mri_voxels_r, vno_nbr, 0, 0, n);
          islc_src = (int)MRIgetVoxVal(mri_voxels_s, vno_nbr, 0, 0, n);

          if (MRIgetVoxVal(mri_sampled, icol_src, irow_src, islc_src, 0) == vno+1)
            continue ; // already sampled this voxel 
          MRIsetVoxVal(mri_sampled, icol_src, irow_src, islc_src, 0, vno+1) ;
          srcval = MRIgetVoxVal(mri_src, icol_src, irow_src, islc_src, 0) ;
          wm = MRIgetVoxVal(mri_wm, icol_src, irow_src, islc_src, 0) ;
          gm = MRIgetVoxVal(mri_gm, icol_src, irow_src, islc_src, 0) ;
          csf = MRIgetVoxVal(mri_csf, icol_src, irow_src, islc_src, 0) ;
          if ((csf < 0.3) && (wm > 0.5) && (srcval > 5000))
            DiagBreak() ;
          if (index >= MAX_SAMPLES)
            break ;
          index++ ;
          if (!FZERO(gm))
          {
            nfound++ ;
            found |= FOUND_GM ;
          }
#if 0
          else
            continue ; // don't use
#endif
          if (!FZERO(wm))
            found |= FOUND_WM ;
          if (!FZERO(csf))
            found |= FOUND_CSF ;
          *MATRIX_RELT(m_A, index, 1) = gm ;
          *MATRIX_RELT(m_A, index, 2) = wm ;
          *MATRIX_RELT(m_A, index, 3) = csf ;
          *MATRIX_RELT(m_S, index, 1) = srcval ;
        }
      }
#if 0
      done = (nfound >= min_needed) &&
        ((found & FOUND_GM) && (found & FOUND_WM) && (found & FOUND_CSF));
#endif
      done = ((nfound >= (min_needed/2)) && (found & FOUND_GM) && (index >= min_needed) &&
              (nsize >= 2)) ;
      if (nsize++ >= 3 && !done)
      {
        ErrorPrintf(ERROR_UNSUPPORTED, "GM estimation at vno %d couldn't find enough samples",vno) ;
        failed = 1 ;
        break ;
      }
    }  while (!done);

    //    MRIclear(mri_sampled) ;
    if (failed)
    {
      MRIsetVoxVal(mri_samples, vno, 0, 0, 0, -1) ;
      continue ;
    }

    m_A->rows = index ;
    m_S->rows = index ;
    if ((found & FOUND_WM) == 0)  // can't estimate WM intensity
    {
      m_A->cols = m_p->rows = 1 ;
    }
    else if ((found & FOUND_CSF) == 0)  // can't estimate CSF intensity
    {
      m_A->cols = m_p->rows = 2 ;
    }
    m_inv = MatrixSVDPseudoInverse(m_A, NULL) ;
    if (m_inv == NULL)
      DiagBreak() ;
    MatrixMultiply(m_inv, m_S, m_p) ;
    if (*MATRIX_RELT(m_p, 1, 1) < 0)
      DiagBreak() ;
    if (*MATRIX_RELT(m_p, 1, 1) < -10000)
      DiagBreak() ;

    // compute quality of linear fit
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      double rms ;
      MATRIX *m_fit = MatrixMultiply(m_A, m_p, NULL) ;
      MatrixSubtract(m_fit, m_S, m_fit) ;
      MatrixSquareElts(m_fit, m_fit) ;
      rms = sqrt(MatrixSumElts(m_fit) / index) ;
      v->curv = rms ;
      rms_mean += rms ; rms_var += rms*rms ;
      MatrixFree(&m_fit) ;
    }

    MRIsetVoxVal(mri_samples, vno, 0, 0, 0, *MATRIX_RELT(m_p, 1, 1)) ;
    m_A->cols = m_p->rows = 3 ;  // for next time
    MatrixFree(&m_inv) ;
  }

  MRIfree(&mri_sampled) ;
  MRIfree(&mri_voxels_c) ;
  MRIfree(&mri_voxels_r) ;
  MRIfree(&mri_voxels_s) ;
  MatrixFree(&Scrs) ; MatrixFree(&Txyz) ; MatrixFree(&m_A);
  MatrixFree(&m_p) ; MatrixFree(&m_S) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    rms_mean /= mris->nvertices ;
    rms_var = rms_var / mris->nvertices - rms_mean*rms_mean ;
    printf("RMS = %2.2f +- %2.2f\n", rms_mean, sqrt(rms_var)) ;
    MRISwriteCurvature(mris, "rms_fit") ;
  }
  return(mri_samples) ;
}

