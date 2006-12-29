/**
 * @file  mri_vol2surf.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:09 $
 *    $Revision: 1.34 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


/*----------------------------------------------------------
  Name: vol2surf.c
  $Id: mri_vol2surf.c,v 1.34 2006/12/29 02:09:09 nicks Exp $
  Author: Douglas Greve
  Purpose: Resamples a volume onto a surface. The surface
  may be that of a subject other than the source subject.
  This replaces paint.c and is also used to convert functional
  data onto the icosahedron.

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
#include <sys/time.h>

#include "icosahedron.h"

#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mri.h"
#include "mri_identify.h"
#include "mri2.h"
#include "prime.h"
#include "fsenv.h"

//#include "bfileio.h"
#include "registerio.h"
#include "resample.h"
#include "selxavgio.h"
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
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_vol2surf.c,v 1.34 2006/12/29 02:09:09 nicks Exp $";
char *Progname = NULL;

char *defaulttypestring;
int  defaulttype = MRI_VOLUME_TYPE_UNKNOWN;

char *srcvolid   = NULL;
char *srctypestring = NULL;
int   srctype = MRI_VOLUME_TYPE_UNKNOWN;
char *srcregfile = NULL;
int  regheader = 0;
char *srcwarp    = NULL;
int   srcoldreg  = 0;
char *srcsubject = NULL;
char *srcsubjectuse = NULL;

char *srchitvolid   = NULL;
char *srchittypestring = NULL;
int   srchittype = MRI_VOLUME_TYPE_UNKNOWN;

char *hemi    = NULL;
char *surfname = "white";
char *trgsubject = NULL;
int  IcoOrder = -1;
float IcoRadius = 100;
char *surfreg = "sphere.reg";
char *thicknessname = "thickness";
float ProjFrac = 0;
int   ProjDistFlag = 0;
float ProjFracMin=0.0,ProjFracMax=0.0,ProjFracDelta=1.0;

MRI_SURFACE *Surf    = NULL;
MRI_SURFACE *SurfOut = NULL;
MRI_SURFACE *SrcSurfReg = NULL;
MRI_SURFACE *TrgSurfReg = NULL;
int UseHash = 1;

char *outfile  = NULL;
char *outtypestring = NULL;
int  outtype = MRI_VOLUME_TYPE_UNKNOWN;

char *srchitfile = NULL;
char *trghitfile = NULL;

char  *interpmethod_string = "nearest";
int    interpmethod = -1;
char *mapmethod = "nnfr";

int debug = 0;
int reshape = 1;
int reshapefactor = 0;
int reshapetarget = 20;

MATRIX *Dsrc, *Dsrctmp, *Wsrc, *Fsrc, *Qsrc, *vox2ras;
SXADAT *sxa;

char *SUBJECTS_DIR = NULL;
MRI *SrcVol, *SurfVals, *SurfVals2, *SurfValsP;
MRI *SrcHits, *SrcDist, *TrgHits, *TrgDist;
MRI *mritmp;
MRI *SrcHitVol;
MRI *TargVol;

FILE *fp;

char *nvoxfile = NULL;

char tmpstr[2000];

int   float2int_src;
char  *float2int_string = "round";
int    float2int = -1;
int   fixtkreg = 0;
int ReverseMapFlag = 0;

int framesave = 0;

float fwhm = 0, gstd = 0;

int  srcsynth = 0;
long seed = -1; /* < 0 for auto */
char *seedfile = NULL;

double scale = 0;
int GetProjMax = 0;

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
int main(int argc, char **argv) {
  int n,err, f, vtx, svtx, tvtx, nproj;
  int nrows_src, ncols_src, nslcs_src, nfrms;
  float ipr, bpr, intensity;
  float colres_src=0, rowres_src=0, slcres_src=0;
  float *framepower = NULL;
  char fname[2000];
  int nTrg121,nSrc121,nSrcLost;
  int nTrgMulti,nSrcMulti;
  float MnTrgMultiHits,MnSrcMultiHits;
  int nargs;
  int r,c,s,nsrchits;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_vol2surf.c,v 1.34 2006/12/29 02:09:09 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

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

  if (srcsynth == 0) {
    /* Load the Source Volume */
    SrcVol =  MRIreadType(srcvolid,srctype);
    if (SrcVol == NULL) {
      printf("ERROR: could not read %s as type %d\n",srcvolid,srctype);
      exit(1);
    }
    if (SrcVol->type != MRI_FLOAT) {
      printf("INFO: chaning type to float\n");
      SrcVol = MRISeqchangeType(SrcVol,MRI_FLOAT,0,0,0);
    }
    printf("Done loading volume\n");
  } else {
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

  if (!regheader) {
    /* Load the registration matrix */
    err = regio_read_register(srcregfile, &srcsubject, &ipr, &bpr,
                              &intensity, &Dsrc, &float2int_src);
    if (err) exit(1);
    if (srcsubjectuse) {
      if (strcmp(srcsubject,srcsubjectuse)) {
        printf("INFO: overriding regsubject %s with %s\n",
               srcsubject,srcsubjectuse);
      }
      srcsubject = srcsubjectuse;
    }
  } else {
    /* compute the registration from the header */
    sprintf(tmpstr,"%s/%s/mri/orig.mgz",SUBJECTS_DIR,srcsubject);
    printf("Computing registration from header.\n");
    printf("  Using %s as target reference.\n",tmpstr);
    TargVol = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
    if (TargVol == NULL) exit(1);
    Dsrc = MRItkRegMtx(TargVol,SrcVol,NULL);
    MRIfree(&TargVol);
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

  /* check that the resolution of the SrcVol is the same as in reg file */
  if (!regheader &&
      (fabs(SrcVol->xsize-ipr) > .01 || fabs(SrcVol->zsize-bpr) > .01)) {
    printf("WARNING: the voxel resolution in the source volume (%g,%g,%g) differs \n",
           SrcVol->xsize,SrcVol->ysize,SrcVol->zsize);
    printf("         from that listed in the registration file (%g,%g,%g)\n",
           ipr,ipr,bpr);
    if (fixtkreg) {
      printf("ERROR: cannot fix tkreg matrix with resolution inconsistency\n");
      exit(1);
    }
  }

  /* Fix tkregister matrix if necessary */
  if (fixtkreg && float2int == FLT2INT_TKREG) {
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

  /* Wsrc: Get the source warping Transform */
  Wsrc = NULL;
  /* Fsrc: Get the source FOV registration matrix */
  Fsrc = NULL;
  // Compute vox2ras for source
  vox2ras = MRIxfmCRS2XYZtkreg(SrcVol);
  // Compute ras2vox (Qsrc: the quantization matrix)
  Qsrc = MatrixInverse(vox2ras,NULL);

  /* If this is a statistical volume, raise each frame to it's appropriate
     power (eg, stddev needs to be squared)*/
  if (is_sxa_volume(srcvolid)) {
    printf("INFO: Source volume detected as selxavg format\n");
    sxa = ld_sxadat_from_stem(srcvolid);
    if (sxa == NULL) exit(1);
    framepower = sxa_framepower(sxa,&f);
    if (f != SrcVol->nframes) {
      fprintf(stderr," number of frames is incorrect (%d,%d)\n",
              f,SrcVol->nframes);
      exit(1);
    }
    printf("INFO: Adjusting Frame Power\n");
    fflush(stdout);
    mri_framepower(SrcVol,framepower);
  }
  fflush(stdout);

  if (fwhm > 0) {
    printf("INFO: smoothing volume at fwhm = %g mm (std = %g)\n",fwhm,gstd);
    MRIgaussianSmooth(SrcVol, gstd, 1, SrcVol); /* 1 = normalize */
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
  if (ProjFrac != 0) {
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
  nproj = 0;
  for (ProjFrac=ProjFracMin; ProjFrac <= ProjFracMax; ProjFrac += ProjFracDelta) {
    printf("%2d %g %g %g\n",nproj+1,ProjFrac,ProjFracMin,ProjFracMax);
    SurfValsP = vol2surf_linear(SrcVol, Qsrc, Fsrc, Wsrc, Dsrc,
                                Surf, ProjFrac, interpmethod, float2int, SrcHitVol,
                                ProjDistFlag);
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
    if (SurfVals2 == NULL) {
      printf("ERROR: mapping surfaces\n");
      exit(1);
    }
    printf("Done mapping surfaces\n");
    fflush(stdout);

    /*Compute some stats on mapping number of trgvtxs mapped from a source vtx*/
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
    printf("nSrc121 = %5d, nSrcLost = %5d, nSrcMulti = %5d, MnSrcMultiHits = %g\n",
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

    if (outtypestring != NULL &&
        (!strcasecmp(outtypestring,"w") ||
         !strcasecmp(outtypestring,"paint")) )
      SurfOut = TrgSurfReg;
    else
      MRISfree(&TrgSurfReg);

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
  } else {
    SurfVals2 = SurfVals;
    SurfOut = Surf;
  }

  if (scale != 0) {
    printf("Rescaling output by %g\n",scale);
    MRImultiplyConst(SurfVals2,scale,SurfVals2);
  }

  /* If this is a statistical volume, lower each frame to it's appropriate
     power (eg, variance needs to be sqrt'ed) */
  if (is_sxa_volume(srcvolid)) {
    printf("INFO: Readjusting Frame Power\n");
    fflush(stdout);
    for (f=0; f < SurfVals2->nframes; f++) framepower[f] = 1.0/framepower[f];
    mri_framepower(SurfVals2,framepower);
    sxa->nrows = 1;
    sxa->ncols = SurfVals2->width;
    sv_sxadat_by_stem(sxa,outfile);
  }

  if (outtypestring != NULL &&
      (!strcasecmp(outtypestring,"w") || !strcasecmp(outtypestring,"paint"))) {
    /*-------------- paint or .w --------------*/
    for (vtx = 0; vtx < SurfVals2->width; vtx++)
      SurfOut->vertices[vtx].val = MRIFseq_vox(SurfVals2,vtx,0,0,framesave) ;
    MRISwriteValues(SurfOut, outfile) ;
  } else {
    if (reshape) {
      if (reshapefactor == 0) {
        if (outtype == MRI_ANALYZE4D_FILE || outtype == MRI_ANALYZE_FILE ||
            outtype == NIFTI1_FILE || outtype == NII_FILE)
          reshapefactor = GetClosestPrimeFactorLess(SurfVals2->width,32768);
        else
          reshapefactor = GetClosestPrimeFactor(SurfVals2->width,reshapetarget);
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
    printf("Writing to %s\n",outfile);
    printf("Dim: %d %d %d\n",SurfVals2->width,SurfVals2->height,SurfVals2->depth);
    MRIwriteType(SurfVals2,outfile,outtype);
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
    else if (!strcasecmp(option, "--reshape"))   reshape = 1;
    else if (!strcasecmp(option, "--noreshape")) reshape = 0;
    else if (!strcasecmp(option, "--fixtkreg")) fixtkreg = 1;
    else if (!strcasecmp(option, "--nofixtkreg")) fixtkreg = 0;

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
    } else if (!strcmp(option, "--regheader")) {
      if (nargc < 1) argnerr(option,1);
      regheader = 1;
      srcsubject = pargv[0];
      nargsused = 1;
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
    } else if (!strcmp(option, "--projdist-int") ||
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
        fprintf(stderr,"  must be either nearest, tli, or sinc\n");
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
    } else if (!strcmp(option, "--o") || !strcmp(option, "--out")) {
      if (nargc < 1) argnerr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--out_type") || !strcmp(option, "--ofmt")) {
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
      nargsused = 1;
    } else if (!strcmp(option, "--rft")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&reshapetarget);
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
    } else if (!strcmp(option, "--srcsynth")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%ld",&seed);
      srcsynth = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--seedfile")) {
      if (nargc < 1) argnerr(option,1);
      seedfile = pargv[0];
      nargsused = 1;
    } else {
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
  printf("   --reg source registration  \n");
  printf("   --regheader subject\n");
  printf("   --float2int float-to-int conversion method "
         "(<round>, tkregister )\n");
  printf("   --fixtkreg : make make registration matrix round-compatible\n");
  printf("   --fwhm fwhm : smooth input volume (mm)\n");
  printf("\n");
  printf("   --trgsubject target subject (if different than reg)\n");
  printf("   --hemi       hemisphere (lh or rh) \n");
  printf("   --surf       target surface (white) \n");
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
  printf("   --projdist-max min max del : max along normal\n");
  //printf("   --thickness thickness file (thickness)\n");
  printf("\n");
  printf(" Options for output\n");
  printf("   --out       output path\n");
  printf("   --out_type  output format\n");
  printf("   --frame     save only nth frame (with paint format)\n");
  printf("   --noreshape do not save output as multiple 'slices'\n");
  printf("   --rf R  integer reshaping factor, save as R 'slices'\n");
  printf("   --srchit   volume to store the number of hits at each vox \n");
  printf("   --srchit_type  source hit volume format \n");
  printf("   --nvox nvoxfile : write number of voxels intersecting surface\n");
  printf("\n");
  printf(" Other Options\n");
  printf("   --scale scale : multiply all intensities by scale.\n");
  printf("   --srcsynth seed : synthesize source volume\n");
  printf("   --seedfile fname : save synth seed to fname\n");
  printf("   --sd SUBJECTS_DIR \n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  //  printf("   --interp    interpolation method (<nearest>, tli, or sinc)\n");
  printf("%s\n", vcid) ;
  printf("\n");
  //printf("   --src_type  input volume format \n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf(
    "This program will resample a volume onto a surface of a subject or the \n"
    "sphere. The output can be viewed on the surface (using tksurfer) or can \n"
    "be using for surface-based intersubject averaging. This program supersedes\n"
    "paint.\n"
    "\n"
    "OPTIONS\n"
    "\n"
    "  --mov path to input volume (see below). Can also use --src\n"
    "\n"
    "  --reg file : registration file as computed by tkregister,\n"
    "    tkmedit, mri_make_register, or spmregister. This file\n"
    "    has the following format:\n"
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
    "    used to register the volume, the method should be blank or 'tkregister'\n"
    "    (no quotes), otherwise it should be 'round'. This can be overridden \n"
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
    "    matrix and change the float2int method to round. Don't use this flag \n"
    "    unless you know what you are doing. \n"
    "\n"
    "  --fwhm FWHM\n"
    "\n"
    "    Smooth input volume with a gaussian kernal with FWHM mm.\n"
    "\n"
    "  --trgsubject target subject name : resample volume onto this subject \n"
    "    instead of the one found in the registration file. The target subject \n"
    "    can be either a subject name (as found in $SUBJECTS_DIR) or ico \n"
    "    (to map onto the sphere). If the target subject is not the source \n"
    "    subject, then the surfaces are mapped using each subject's spherical\n"
    "    surface registration (?h.sphere.reg or that specified with --surfreg).\n"
    "    If the target subject is ico, then the volume is resampled onto an\n"
    "    icosahedron, which is used to uniformly sample a sphere. This requires\n"
    "    specifying the icosahedron order (see --icoorder).\n"
    "\n"
    "\n"
    "  --hemi hemisphere : lh = left hemisphere, rh = right hemisphere\n"
    "\n"
    "  --surf surfacename : the surface on which to resample. The default is\n"
    "    white. It will look for $SUBJECTS_DIR/subjectname/surf/?h.surfacename\n"
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
    "    middle of the cortical surface. This requires that a ?h.thickness file \n"
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
    "  --out  output path : location to store the data (see below)\n"
    "  --out_type format of output (see below)\n"
    "\n"
    "  --frame 0-based frame number : sample and save only the given frame \n"
    "    from the source volume (needed when out_type = paint). Default 0.\n"
    "\n"
    "  --noreshape : by default, mri_vol2surf will save the output as multiple\n"
    "    'slices'. This is for logistical purposes (eg, in the analyze format\n"
    "    the size of a dimension cannot exceed 2^15). Use this flag to prevent\n"
    "    this behavior. This has no effect when the output type is paint.\n"
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
  printf("%s\n", vcid) ;
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
    fprintf(stderr,"  must be either nearest, tli, or sinc\n");
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


