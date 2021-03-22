/**
 * @brief Converts a label to a segmentation volume.
 *
 * Converts a label or a set of labels into a volume. For a single label,
 * the volume will be binary: 1 where the label is and 0 where it is not.
 * For multiple labels, the volume will be 0 where no labels were found
 * otherwise the value will the the label number. For a voxel to be
 * declared part of a label, it must have enough hits in the voxel and
 * it must have more hits than any other label.
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "error.h"
#include "diag.h"
#include "proto.h"

#include "mri.h"
#include "mrisutils.h"
#include "MRIio_old.h"
#include "mri_identify.h"
#include "mri2.h"
#include "tags.h"
#include "matrix.h"
#include "version.h"
#include "registerio.h"
#include "resample.h"
#include "annotation.h"
#include "macros.h"
#include "colortab.h"
#include "cmdargs.h"
#include "region.h"
#include "resample.h"
#include "fsenv.h"

#define PROJ_TYPE_NONE 0
#define PROJ_TYPE_ABS  1
#define PROJ_TYPE_FRAC 2

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
static int  checkhemi(char *hemi);
static int get_proj_type_id(char *projtype);
static int get_crs(MATRIX *Tras2vox, float x, float y, float z,
                   int *c, int *r, int *s, MRI *vol);
static int is_surface_label(LABEL *label);
static int load_annotation(char *annotfile, MRIS *Surf);
static int *NthLabelMap(MRI *aseg, int *nlabels);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

char *LabelList[100];
int nlabels = 0;

char *TempVolId = NULL;
char *RegMatFile = NULL;
int RegHeader = 0; 
int RegIdentity = 0; 
int InvertMtx = 0;
double FillThresh = 0.0;
double ProjDelta = .1;
double ProjStart = 0;
double ProjStop = 0;
char *ProjType = NULL;
char *OutVolId = NULL;
char *HitVolId = NULL;
char *PVFVolId = NULL;
char *subject = NULL;
char *hemi    = NULL;
char *AnnotFile = NULL;
char *ASegFSpec = NULL;
const char *SurfId  = "white";
char *LabelVolFile = NULL;

MRI_SURFACE *Surf=NULL;
MRI *OutVol=NULL, *TempVol=NULL, *HitVol=NULL, *PVFVol=NULL;
MRI *ASeg=NULL;
MRI *LabelStatVol=NULL;
char *LabelStatVolFSpec=NULL;
int DoLabelStatVol=0;

int debug;

LABEL *srclabel;
MATRIX *R, *Tvox2ras, *Tras2vox;
double ProjDepth;
int   ProjTypeId;
int   DoProj = 0;
int   SurfNeeded = 0;
char  *SUBJECTS_DIR = NULL;
const char  *thicknessname = "thickness";
char  fname[1000];
double TempVoxVol;
double LabelVoxVol = 1;
double nHitsThresh;
int *ASegLabelList;
int LabelCode;
int UseNativeVox2RAS=0;
int UseNewASeg2Vol=0;

int UseAParcPlusASeg = 0;
int DoStatThresh = 0;
double StatThresh = -1;
int LabelCodeOffset = 0;
COLOR_TABLE *ctTissueType=NULL;
int UpsampleFactor = -1;
double resmm=0;
int FillRibbon = 0;
LTA *lta=NULL;
MRI *ToTkrTemplate=NULL;

MRI *MRIsurfaceLabel2VolOpt(MRI *ribbon, MRIS *surf, LABEL **labels, int nlabels, LTA *Q, 
			    int DoStatThresh, double StatThresh, int DoLabelStatVol, 
			    MRI *vollabel);

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  int  nargs, nthlabel, float2int, err, nthpoint, vtxno;
  float ipr,bpr,intensity;
  char *regsubject;
  float x,y,z,voxvol;
  int c,r,s, oob, nhits, nhitsmax, nhitsmax_label;
  MRI *LabelVol;
  FSENV *fsenv;
  MRI *ribbon;
  char tmpstr[2000];
  LABEL **labels;

  std::string cmdline = getAllInfo(argc, argv, "mri_label2vol");

  nargs = handleVersionOption(argc, argv, "mri_label2vol");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0];

  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);
  printf("%s\n",getVersion().c_str());

  if(TempVolId){
    // Load the template volume
    TempVol = MRIreadHeader(TempVolId,MRI_VOLUME_TYPE_UNKNOWN);
    if (TempVol == NULL) {
      printf("ERROR: reading %s header\n",TempVolId);
    exit(1);
    }
    if (!UseNativeVox2RAS)
      Tvox2ras = MRIxfmCRS2XYZtkreg(TempVol);
    else
      Tvox2ras = MRIgetVoxelToRasXform(TempVol);

    Tras2vox = MatrixInverse(Tvox2ras,NULL);
    printf("Template RAS-to-Vox: --------\n");
    MatrixPrint(stdout,Tras2vox);

    TempVoxVol = TempVol->xsize * TempVol->ysize * TempVol->zsize;
    nHitsThresh = FillThresh*TempVoxVol/LabelVoxVol;
    printf("Template Voxel Volume: %g\n",TempVoxVol);
    printf("nHits Thresh: %g\n",nHitsThresh);

    if(RegMatFile != NULL) {
      // Load registration matrix
      printf("Loading registration from %s\n",RegMatFile);
      err = regio_read_register(RegMatFile, &regsubject, &ipr, &bpr,
				&intensity, &R, &float2int);
      if (err) exit(1);
    } 
    else {
      if(RegHeader){    
	// RegHeader
	printf("Computing registration based on header\n");
	LabelVol = MRIreadHeader(LabelVolFile, MRI_VOLUME_TYPE_UNKNOWN);
	if(LabelVol == NULL){
	  printf("ERROR: reading %s\n",LabelVolFile);
	  exit(1);
	}
	R = MRItkRegMtx(LabelVol, TempVol, NULL);
      }
      if(RegIdentity){
	printf("Using Identity Matrix\n");
	R = MatrixIdentity(4,NULL);
      }
    }
    printf("RegMat: --------\n");
    MatrixPrint(stdout,R);
    if (InvertMtx) {
      printf("Inverting matrix\n");
      MatrixInverse(R,R);
      printf("RegMat: --------\n");
      MatrixPrint(stdout,R);
    }
    
    MatrixMultiply(Tras2vox,R,Tras2vox);
    
    printf("Label RAS-to-Vox: --------\n");
    MatrixPrint(stdout,Tras2vox);
  }

  if (SurfNeeded) {
    // Load the surface used for projection
    Surf = MRISloadSurfSubject(subject,hemi,SurfId,SUBJECTS_DIR);
    if (Surf == NULL) {
      printf("ERROR: could not load surface.\n");
      exit(1);
    }
    // Load the thickness for projection along the normal
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,
            hemi,thicknessname);
    printf("Reading thickness %s\n",fname);
    MRISreadCurvatureFile(Surf, fname);

    if (AnnotFile != NULL) {
      // Load annotation file, get number of labels from it
      nlabels = load_annotation(AnnotFile,Surf);
      if (nlabels==0) {
        printf("ERROR: load annot file %s\n",AnnotFile);
        exit(1);
      }
    }
  }
  /*-------------------------------------------------------------*/
  if(UseAParcPlusASeg){
    if(subject == NULL) subject = regsubject;
    sprintf(fname,"%s/%s/mri/aparc+aseg.mgz",SUBJECTS_DIR,subject);
    ASegFSpec = strcpyalloc(fname);
  }

  if(ASegFSpec != NULL) {
    ASeg = MRIread(ASegFSpec);
    if (ASeg == NULL) {
      printf("ERROR: loading aseg %s\n",ASegFSpec);
      exit(1);
    }
    if(! UseNewASeg2Vol) {
      LTA *aseg2vol;
      aseg2vol = TransformRegDat2LTA(ASeg, TempVol, R);
      OutVol = MRIaseg2volMU(ASeg, aseg2vol, FillThresh, &HitVol, UpsampleFactor, ctTissueType);
      //OutVol = MRIaseg2vol(ASeg, R,TempVol, FillThresh, &HitVol, ctTissueType);
      MRIaddCommandLine(OutVol, cmdline) ;
      err=MRIwrite(OutVol,OutVolId);
      if(err) exit(1);
      if (HitVolId != NULL) {
        MRIaddCommandLine(HitVol, cmdline) ;
        MRIwrite(HitVol,HitVolId);
      }
      if(PVFVolId != NULL) {
	printf("PVF %s\n",PVFVolId);
	printf("Computing PVF %g\n",TempVoxVol);
	PVFVol = MRImultiplyConst(HitVol, 1.0/TempVoxVol, NULL);
	MRIwrite(PVFVol,PVFVolId);
      }
      printf("mri_label2vol done\n");
      exit(0);
    }
    else {
      LTA *aseg2vol;
      int *segidlist,nsegs;
      aseg2vol = TransformRegDat2LTA(ASeg, TempVol, R);
      segidlist = MRIsegIdListNot0(ASeg, &nsegs, 0);
      if(resmm == 0) resmm = ASeg->xsize/UpsampleFactor;
      OutVol = MRIseg2SegPVF(ASeg, aseg2vol, resmm, segidlist, nsegs, NULL, 0, ctTissueType,NULL);
      MRIaddCommandLine(OutVol, cmdline) ;
      err=MRIwrite(OutVol,OutVolId);
      if(err) exit(1);
      printf("mri_label2vol done\n");
      exit(0);
    }
    ASegLabelList = NthLabelMap(ASeg, &nlabels);
  }
  printf("nlabels = %d\n",nlabels);

  if(FillRibbon){
    fsenv = FSENVgetenv();
    sprintf(tmpstr,"%s/%s/mri/ribbon.mgz",fsenv->SUBJECTS_DIR,subject);
    printf("Loading %s\n",tmpstr);
    ribbon = MRIread(tmpstr);
    if(ribbon==NULL) exit(1);
    
    labels = (LABEL **) calloc(sizeof(LABEL *),nlabels);
    for(nthlabel = 0; nthlabel < nlabels; nthlabel++) {
      if(AnnotFile == NULL && ASegFSpec == NULL) {
	printf("Loading %s\n",LabelList[nthlabel]);
	labels[nthlabel] = LabelRead(NULL, LabelList[nthlabel]);
	if(labels[nthlabel] == NULL) {
	  printf("ERROR reading %s\n",LabelList[nthlabel]);
	  exit(1);
	}
	if(labels[nthlabel]->coords != LABEL_COORDS_TKREG_RAS){
	  if(ToTkrTemplate == NULL){
	    printf("ERROR: label is not in tkreg coords so you need to specify a"
		   "template with --tkr-template MRI. The MRI should be the volume"
		   "that the label was created on\n");
	    exit(1);
	  }
	  printf("  Converting to TkRegister RAS\n");
	  LabelToSurfaceRAS(srclabel, ToTkrTemplate, srclabel);
	}
      }
      if(AnnotFile != NULL) {
	labels[nthlabel] = annotation2label(nthlabel,Surf);
	if(labels[nthlabel] == NULL) continue;
      }
    }

    printf("Mapping\n");
    OutVol = MRIsurfaceLabel2VolOpt(ribbon, Surf, labels, nlabels, lta,
				    DoStatThresh, StatThresh, 0, NULL);
    MRIaddCommandLine(OutVol, cmdline) ;
    err=MRIwrite(OutVol,OutVolId);
    if(err) exit(err);
    
    if(DoLabelStatVol) {
      printf("Mapping stats\n");
      LabelStatVol = MRIsurfaceLabel2VolOpt(ribbon, Surf, labels, nlabels, lta,
					    DoStatThresh, StatThresh, DoLabelStatVol, NULL);
      err=MRIwrite(LabelStatVol,LabelStatVolFSpec);
      if(err) exit(err);
    }
    exit(0);
  }
  //---------------------------------------------------------------




  // Create hit volume based on template, one frame for each label
  printf("Allocating Hit Volume (%ld) voxels\n",
	 (long)TempVol->width*TempVol->height*TempVol->depth*nlabels );
  HitVol = MRIallocSequence(TempVol->width, TempVol->height,
                            TempVol->depth, MRI_SHORT, nlabels );
  if (HitVol == NULL) {
    printf("ERROR: could not alloc hit volume\n");
    exit(1);
  }
  MRIcopyHeader(TempVol,HitVol);
  HitVol->nframes = nlabels;

  if(DoLabelStatVol){
    LabelStatVol = MRIallocSequence(TempVol->width, TempVol->height,
				    TempVol->depth, MRI_FLOAT, 1);
    if(LabelStatVol == NULL) {
      printf("ERROR: could not alloc stat volume\n");
      exit(1);
    }
    MRIcopyHeader(TempVol,LabelStatVol);
    LabelStatVol->nframes = 1;
  }

  // Go through each label
  for (nthlabel = 0; nthlabel < nlabels; nthlabel++) {
    if (debug) {
      printf("%2d ",nthlabel);
      fflush(stdout);
      if (nthlabel%20 == 19) {
        printf("\n");
        fflush(stdout);
      }
    }

    if (AnnotFile == NULL && ASegFSpec == NULL) {
      printf("Loading %s\n",LabelList[nthlabel]);
      srclabel = LabelRead(NULL, LabelList[nthlabel]);
      if (srclabel == NULL) {
        printf("ERROR reading %s\n",LabelList[nthlabel]);
        exit(1);
      }
      if(srclabel->coords != LABEL_COORDS_TKREG_RAS){
	if(ToTkrTemplate == NULL){
	  printf("ERROR: label is not in tkreg coords so you need to specify a"
		 "template with --tkr-template MRI. The MRI should be the volume"
		 "that the label was created on\n");
	  exit(1);
	}
	printf("  Converting label coords to TkRegister RAS\n");
	LabelToSurfaceRAS(srclabel, ToTkrTemplate, srclabel);
      }
    }
    if (AnnotFile != NULL) {
      srclabel = annotation2label(nthlabel,Surf);
      if (srclabel == NULL) continue;
    }
    if (ASegFSpec != NULL) {
      printf("   %d/%d ASeg Code: %d\n",nthlabel,nlabels,ASegLabelList[nthlabel]);
      srclabel = LabelfromASeg(ASeg, ASegLabelList[nthlabel]);
      //LabelWrite(lb,"./tmp.label");
    }

    if (DoProj && !is_surface_label(srclabel)) {
      printf("ERROR: label %s is not a surface label.\n",
             LabelList[nthlabel]);
      exit(1);
    }

    // Go through each point in the label
    for (nthpoint = 0; nthpoint < srclabel->n_points; nthpoint++) {
      if(DoStatThresh) 
	if(srclabel->lv[nthpoint].stat < StatThresh) continue; 

      if (debug) printf("  nthpoint = %d\n",nthpoint);

      if (DoProj) { // Project along the surface normal
        vtxno = srclabel->lv[nthpoint].vno;
        ProjDepth = ProjStart;

        while (ProjDepth <= ProjStop) {

          if (ProjTypeId == PROJ_TYPE_ABS)
            ProjNormDist(&x,&y,&z,Surf,vtxno,ProjDepth);
          if (ProjTypeId == PROJ_TYPE_FRAC)
            ProjNormFracThick(&x,&y,&z,Surf,vtxno,ProjDepth);
          oob = get_crs(Tras2vox,x,y,z,&c,&r,&s,TempVol);

          if (debug) printf("   ProjDepth %g   %g %g %g (%g)  %d %d %d   %d\n",
                              ProjDepth,x,y,z,Surf->vertices[vtxno].curv,c,r,s,oob);

          // Accumulate hit volume
          if(!oob){
	    MRISseq_vox(HitVol,c,r,s,nthlabel) ++;
	    if(DoLabelStatVol) 
	      MRIFseq_vox(LabelStatVol,c,r,s,0) = srclabel->lv[nthpoint].stat;
	  }

          ProjDepth += ProjDelta;
          if (ProjDelta == 0) break; // only do once

        } // end loop through projection depths
      }// end Do Projection

      else { // Simply compute crs from the label xyz
        x = srclabel->lv[nthpoint].x;
        y = srclabel->lv[nthpoint].y;
        z = srclabel->lv[nthpoint].z;
        oob = get_crs(Tras2vox,x,y,z,&c,&r,&s,TempVol);
        if (debug) printf("   %g %g %g   %d %d %d   %d\n",x,y,z,c,r,s,oob);
        if (oob) continue; // Out of the volume
        MRISseq_vox(HitVol,c,r,s,nthlabel) ++;
	if(DoLabelStatVol) 
	  MRIFseq_vox(LabelStatVol,c,r,s,0) = srclabel->lv[nthpoint].stat;
      }
    } // end loop over label points

    LabelFree(&srclabel) ;
  } // End loop over labels
  printf("\n");

  if(DoLabelStatVol){
    err=MRIwrite(LabelStatVol,LabelStatVolFSpec);
    if(err) exit(1);
  }

  if(HitVolId != NULL) MRIwrite(HitVol,HitVolId);
  printf("PVF %s\n",PVFVolId);
  if(PVFVolId != NULL) {
    voxvol = TempVol->xsize * TempVol->ysize * TempVol->zsize;
    printf("Computing PVF %g\n",voxvol);
    PVFVol = MRImultiplyConst(HitVol, 1/voxvol, NULL);
    MRIwrite(PVFVol,PVFVolId);
    exit(1);
  }

  // Create output volume based on template, but use 1 frame int
  OutVol = MRIallocSequence(TempVol->width, TempVol->height,
                            TempVol->depth, MRI_INT, 1);
  if (OutVol == NULL) {
    printf("ERROR: could not alloc output volume\n");
    exit(1);
  }
  MRIcopyHeader(TempVol,OutVol);
  OutVol->nframes = 1;

  printf("Thesholding hit volume.\n");
  // Threshold hit volumes and set outvol to nthlabel+1
  for (c=0; c < OutVol->width; c++) {
    for (r=0; r < OutVol->height; r++) {
      for (s=0; s < OutVol->depth; s++) {
        nhitsmax = 0;
        nhitsmax_label = -1;
        for (nthlabel = 0; nthlabel < nlabels; nthlabel++) {
          nhits = (int)MRIgetVoxVal(HitVol,c,r,s,nthlabel);
          //nhits = MRIIseq_vox(HitVol,c,r,s,nthlabel);
          if (nhits <= nHitsThresh) continue;
          if (nhitsmax < nhits) {
            nhitsmax = nhits;
            nhitsmax_label = nthlabel;
          }
        }
        if (nhitsmax_label == -1)
          LabelCode = 0; // No hits -- Unknown, dont add Offset
        else if (ASegFSpec != NULL)
          LabelCode = ASegLabelList[nhitsmax_label] + LabelCodeOffset;
        else if (AnnotFile != NULL)
          LabelCode = nhitsmax_label + LabelCodeOffset;//dont +1, keeps consist with ctab
        else
          LabelCode = nhitsmax_label + 1 + LabelCodeOffset;
        MRIIseq_vox(OutVol,c,r,s,0) = LabelCode;

      }
    }
  }

  // Save out volume
  MRIaddCommandLine(OutVol, cmdline) ;
  err=MRIwrite(OutVol,OutVolId);
  if(err) exit(1);

  printf("mri_label2vol done \n");
  return(0);
}
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
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
    else if (!strcasecmp(option, "--native-vox2ras"))  UseNativeVox2RAS = 1;
    else if (!strcasecmp(option, "--invertmtx"))  InvertMtx = 1;
    else if (!strcasecmp(option, "--new-aseg2vol"))  UseNewASeg2Vol = 1;
    else if (!strcasecmp(option, "--no-new-aseg2vol"))  UseNewASeg2Vol = 0;
    else if (!strcasecmp(option, "--aparc+aseg"))  UseAParcPlusASeg = 1;
    else if (!strcasecmp(option, "--ttype")) 
      ctTissueType = TissueTypeSchema(NULL,"default-jan-2014");
    else if (!strcasecmp(option, "--ttype+head")) 
      ctTissueType = TissueTypeSchema(NULL,"default-jan-2014+head");
    else if (!strcasecmp(option, "--fill-ribbon")) FillRibbon = 1;

    else if (!strcmp(option, "--surf")) {
      if (nargc < 1) argnerr(option,1);
      SurfId = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--label")) {
      if (nargc < 1) argnerr(option,1);
      LabelList[nlabels] = pargv[0];
      nlabels++;
      nargsused = 1;
    }
    else if (!strcmp(option, "--label-stat")) {
      if (nargc < 1) argnerr(option,1);
      DoLabelStatVol = 1;
      LabelStatVolFSpec = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--annot")) {
      if (nargc < 1) argnerr(option,1);
      AnnotFile = pargv[0];
      SurfNeeded = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--seg")) {
      if (nargc < 1) argnerr(option,1);
      ASegFSpec = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--upsample")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&UpsampleFactor);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--resmm")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&resmm);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--temp")) {
      if (nargc < 1) argnerr(option,1);
      TempVolId = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--tkr-template")) {
      if (nargc < 1) argnerr(option,1);
      ToTkrTemplate = MRIreadHeader(pargv[0],MRI_VOLUME_TYPE_UNKNOWN);
      if(ToTkrTemplate==NULL) exit(1);
      nargsused = 1;
    }
    else if (!strcmp(option, "--reg")) {
      if (nargc < 1) argnerr(option,1);
      RegMatFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--lta")) {
      if (nargc < 1) argnerr(option,1);
      lta = LTAread(pargv[0]);
      if(lta == NULL) exit(1);
      subject = lta->subject;
      nargsused = 1;
    } else if (!strcmp(option, "--fillthresh")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&FillThresh);
      nargsused = 1;
    } else if (!strcmp(option, "--labvoxvol")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&LabelVoxVol);
      nargsused = 1;
    } else if (!strcmp(option, "--proj")) {
      if (nargc < 4) argnerr(option,4);
      ProjType = pargv[0];
      sscanf(pargv[1],"%lf",&ProjStart);
      sscanf(pargv[2],"%lf",&ProjStop);
      sscanf(pargv[3],"%lf",&ProjDelta);
      ProjTypeId = get_proj_type_id(ProjType);
      if (ProjStart > ProjStop) {
        printf("ERROR: projection start must be <= stop\n");
        exit(1);
      }
      if (ProjDelta == 0 && (ProjStart != ProjStop)) {
        printf("ERROR: cannot spec projection delta = 0.\n");
        exit(1);
      }
      if (ProjDelta < 0) {
        printf("ERROR: projection delta must be > 0\n");
        exit(1);
      }
      DoProj = 1;
      nargsused = 4;
    } 
    else if (!strcmp(option, "--subject")) {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--stat-thresh")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&StatThresh);
      DoStatThresh = 1;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--identity")){
      RegIdentity = 1;
    }
    else if (!strcmp(option, "--regheader")) {
      if(CMDnthIsArg(nargc, pargv, 0)){
	LabelVolFile = pargv[0];
	nargsused = 1;
      } else {
	if(ASegFSpec) LabelVolFile = ASegFSpec;
	else{
	  printf("ERROR: --regheader requires one argument unless you\n");
	  printf("are going to use a --seg. If so specify that before --regheader\n");
	  exit(1);
	}
      }
      RegHeader = 1;
    } 
    else if (!strcmp(option, "--hemi")) {
      if (nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      checkhemi(hemi);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--hits")) {
      if (nargc < 1) argnerr(option,1);
      HitVolId = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--pvf")) {
      if (nargc < 1) argnerr(option,1);
      PVFVolId = pargv[0];
      printf("PVF %s\n",PVFVolId);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--sd")) {
      if (nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--o")) {
      if (nargc < 1) argnerr(option,1);
      OutVolId = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--offset")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&LabelCodeOffset);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--defects")) {
      // surf defectno template offset mergeflag
      if(nargc < 6) {
	printf("usage: --defects surf defectno.mgz template.mgz offset mergeflag\n");
	argnerr(option,6);
      }
      MRIS *defsurf  = MRISread(pargv[0]); // surf
      if(defsurf==NULL) exit(1);
      MRI *defno    = MRIread(pargv[1]); //defects no
      if(defno==NULL) exit(1);
      MRI *deftemp  = MRIread(pargv[2]); //template
      if(deftemp==NULL) exit(1);
      int defmerge, defoffset,deferr;
      sscanf(pargv[3],"%d",&defoffset);
      sscanf(pargv[4],"%d",&defmerge);
      if(defmerge==0){
	MRIconst(deftemp->width,deftemp->height,deftemp->depth,deftemp->nframes,0,deftemp); // set MRI to 0
	if(deftemp->type != MRI_INT) {
	  printf("Changing input type %d to MRI_INT\n",deftemp->type);
	  MRI *defmritmp = MRIchangeType(deftemp,MRI_INT,0,0,1);
	  MRIfree(&deftemp);
	  deftemp = defmritmp;
	}
      }
      printf("Converting defects to volume: offset=%d, merge=%d\n",defoffset,defmerge);
      deferr = MRISdefects2Seg(defsurf, defno, defoffset, deftemp);
      if(deferr) exit(1);
      printf("Writing to %s\n",pargv[5]);
      deferr = MRIwrite(deftemp,pargv[5]);
      if(deferr) exit(1);
      MRISfree(&defsurf);
      MRIfree(&defno);
      MRIfree(&deftemp);
      exit(0);
      nargsused = 6;
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
  printf("USAGE: mri_label2vol\n") ;
  printf("\n");
  printf("   --o volid : output volume\n");
  printf("\n");
  printf("     Possible inputs:\n");
  printf("   --label labelid <--label labelid>  \n");
  printf("   --annot annotfile : surface annotation file  \n");
  printf("   --seg   segpath : segmentation\n");
  printf("     --upsample USFactor: upsample seg by factor\n");
  printf("     --resmm ResMM : upsample resolution (instead of USF)\n");
  printf("     --ttype : use default tissue type info for seg\n");
  printf("     --ttype+head : use default+head tissue type info for seg\n");
  printf("     --new-aseg2vol : use new function (eventually will become default)\n");
  printf("   --aparc+aseg  : use aparc+aseg.mgz in subjectdir as seg\n");
  printf("\n");
  printf("   --temp tempvolid : output template volume\n");
  printf("\n");
  printf("   --reg regmat : VolXYZ = R*LabelXYZ\n");
  printf("   --regheader volid : label template volume (needed with --label or --annot)\n");
  printf("       for --seg, use the segmentation volume\n");
  printf("   --identity : set R=I\n");
  printf("   --invertmtx : Invert the registration matrix\n");
  printf("\n");
  printf("   --fillthresh thresh : between 0 and 1 (def 0)\n");
  printf("   --labvoxvol voxvol : volume of each label point (def 1mm3)\n");
  printf("   --proj type start stop delta\n");
  printf("   --fill-ribbon\n");
  printf("\n");
  printf("   --subject subjectid : needed with --proj or --annot\n");
  printf("   --hemi hemi : needed with --proj or --annot\n");
  printf("   --surf surface  : use surface instead of white\n") ;
  printf("\n");
  printf("   --hits hitvolid : each frame is nhits for a label\n");
  printf("   --label-stat statvol : map the label stats field into the vol\n");
  printf("   --stat-thresh thresh : only use label point where stat > thresh\n");
  printf("   --offset k : add k to segmentation numbers (good when 0 should not be ignored)\n");
  printf("   --defects surf defectno.mgz orig.mgz offset mergeflag (ignores all other inputs)\n");
  printf("     Creates a segmentation volume of the surface defects, eg, \n");
  printf("       mri_label2vol --defects surf/lh.orig surf/lh.defect_labels mri/orig.mgz    1000 0 mri/defects.mgz\n");
  printf("       mri_label2vol --defects surf/rh.orig surf/rh.defect_labels mri/defects.mgz 2000 1 mri/defects.mgz\n");
  printf("\n");
  printf("   --native-vox2ras : use native vox2ras xform instead of tkregister-style\n");
  printf("   --version : print version and exit\n");
  printf("   --help \n");
  printf("\n");
  printf("%s\n",getVersion().c_str());
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf(

    "\n"
    "Help Outline:\n"
    "  - SUMMARY  \n"
    "  - ARGUMENTS \n"
    "  - RESOLVING MULTI-LABEL AMBIGUITIES\n"
    "  - CHECKING YOUR RESULTS\n"
    "  - EXAMPLES\n"
    "  - KNOWN BUGS\n"
    "  - SEE ALSO\n"
    "\n"
    "SUMMARY\n"
    "\n"
    "Converts a label or a set of labels into a volume. For a single label,\n"
    "the volume will be binary: 1 where the label is and 0 where it is not.\n"
    "For multiple labels, the volume will be 0 where no labels were found\n"
    "otherwise the value will the the label number. For a voxel to be\n"
    "declared part of a label, it must have enough hits in the voxel and\n"
    "it must have more hits than any other label.\n"
    "\n"
    "ARGUMENTS\n"
    "\n"
    "--label labelfile <--label labelfile>\n"
    "\n"
    "Enter the name of the label file. For multiple labels, use multiple\n"
    "--label flags. Labels can be created manually with tkmedit and\n"
    "tksurfer or automatically from a subcortical segmentation or cortical\n"
    "annotation. Labels are simple text files. The first line is a header.\n"
    "Each following line contains data with 5 columns. The first column is\n"
    "the vertex number of the label point. The next 3 columns are the X, Y,\n"
    "and Z of the point. The last can be ignored. If the label is not a \n"
    "surface-based label, then the vertex number will be -1. Not with --annot\n"
    "or --seg\n"
    "\n"
    "--annot annotfile\n"
    "\n"
    "FreeSurfer can perform automatic parcellation of the cortical surface,\n"
    "the results of which are stored in an annotation file. This is a\n"
    "binary file that assigns each surface vertex to a cortical area\n"
    "designated by a unique number. These annotations can be converted to\n"
    "separate labels using mri_annotation2label. These labels can then be\n"
    "input to mri_label2vol using --label. Or, the annotation file can be\n"
    "read in directly using --annot. The map of annotation numbers to \n"
    "annotation names can be found at FreeSurferColorLUT.txt \n"
    "in $FREESURFER_HOME. Not with --label or --seg.\n"
    "\n"
    "--seg segpath\n"
    "\n"
    "Path to input segmentation. A segmentation is a volume in which each voxel\n"
    "is assigned a number indicating it's class. The output volume will keep\n"
    "the same numbering. Voting is used to resolve multiple segs mapping to\n"
    "a single output voxel. The registration in this case goes from the seg to\n"
    "the template volume. Not with --label or --annot.\n"
    "\n"
    "--temp tempvolid\n"
    "\n"
    "Template volume. The output volume will have the same size and geometry\n"
    "as the template. Template must have geometry information (ie, direction\n"
    "cosines and voxel sizes). Required.\n"
    "\n"
    "--reg regmatfile\n"
    "\n"
    "tkregister-style registration matrix (see tkregister2 --help) which maps\n"
    "the XYZ of the label to the XYZ of the template volume. If not specified,\n"
    "then the identity is assumed.\n"
    "\n"
    "--regheader labelvolume\n"
    "\n"
    "Construct registration based on the header info between the template volume\n"
    "and the labelvolume. The label should have been defined on labelvolume (or\n"
    "its equivalent).\n"
    "\n"
    "--identity\n"
    "\n"
    "Use the identity matrix as the registration.\n"
    "\n"
    "--fillthresh thresh\n"
    "\n"
    "Relative threshold which the number hits in a voxel must exceed for\n"
    "the voxel to be considered a candidate for membership in the label. A\n"
    "'hit' is when a label point falls into a voxel. thresh is a value\n"
    "between 0 and 1 indicating the fraction of the voxel volume that must\n"
    "be filled by label points. The voxel volume is determined from the\n"
    "template. It is assumed that the each label point represents a voxel\n"
    "1mm3 in size (which can be changed with --labvoxvol). So, the actual\n"
    "number of hits needed to exceed threshold is\n"
    "thresh*TempVoxVol/LabVoxVol. The default value is 0, meaning that any\n"
    "label point that falls into a voxel makes that voxel a candidate for\n"
    "membership in the label.  Note: a label must also have the most hits\n"
    "in the voxel before that voxel will actually be assigned to the label\n"
    "in the volume. Note: the label voxel volume becomes a little ambiguous\n"
    "for surface labels, particularly when they are 'filled in' with\n"
    "projection.\n"
    "\n"
    "--labvoxvol voxvol\n"
    "\n"
    "Volume covered by each label point. Default is 1mm3. This only affects\n"
    "the fill threshold (--fillthresh). Note: the label voxel volume\n"
    "becomes a little ambiguous for surface labels, particularly when they\n"
    "are projected.\n"
    "\n"
    "--proj type start stop delta\n"
    "\n"
    "Project the label along the surface normal. type can be abs or frac. \n"
    "abs means that the start, stop, and delta are measured in mm. frac\n"
    "means that start, stop, and delta are relative to the thickness at\n"
    "each vertex. The label definition is changed to fill in label\n"
    "points in increments of delta from start to stop. Requires subject\n"
    "and hemi in order to load in a surface and thickness. Uses the\n"
    "white surface. The label MUST have been defined on the surface.\n"
    "\n"
    "--subject subjectid\n"
    "\n"
    "FREESURFER subject identifier. Needed when using --proj.\n"
    "\n"
    "--hemi hemi\n"
    "\n"
    "Hemisphere to use loading the surface for --proj. Legal values are\n"
    "lh and rh.\n"
    "\n"
    "--o volid\n"
    "\n"
    "Single frame output volume in which each voxel will have the number of\n"
    "the label to which it is assigned (or 0 for no label). The label\n"
    "number is the order in which it appears on the command-line.  Takes\n"
    "any format accepted by mri_convert (eg, nii, mgh, spm, analyze, bshort).\n"
    "\n"
    "--hits hitvolid\n"
    "\n"
    "Hit volume. This is a multi-frame volume, with one frame for each\n"
    "label. The value at each voxel for a given frame is the number of hits\n"
    "that voxel received for that label. This is mostly good as a debugging\n"
    "tool, but you could use it to implement your own multi-label\n"
    "arbitration routine. Or you could binarize to have each label\n"
    "represented separately. Takes any format accepted by mri_convert (eg,\n"
    "nii, spm, analyze, bshort, mgh). With --seg, this is a single frame volume\n"
    "with the number of hits from the winning seg id.\n"
    "\n"
    "--native-vox2ras\n"
    "\n"
    "Use the 'native' voxel-to-RAS transform instead of the tkregister-style.\n"
    "The 'native' vox2ras is what is stored with the volume (usually scanner)\n"
    "This may be needed depending upon how the label was created (eg, with scuba).\n"
    "\n"
    "--label-stat labelstatvol\n"
    "\n"
    "Create a volume in which the value at a voxel is the label stat for that voxel.\n"
    "Note: uses the value from last label point hit. Not recommended for multiple\n"
    "labels.\n"
    "\n"
    "RESOLVING MULTI-LABEL AMBIGUITIES\n"
    "\n"
    "When there are multiple labels, it is possible that more than one\n"
    "label will map to a single voxel in the output volume. When this\n"
    "happens, the voxel is assigned to the label with the most label\n"
    "points in that voxel. Note that the voxel must still pass the \n"
    "fill threshold test in order to be considered part of the label.\n"
    "\n"
    "CHECKING YOUR RESULTS\n"
    "\n"
    "It is very important to check that the conversion of the label to the\n"
    "volume was done correctly. It may be that it is way off or it could be\n"
    "off around the edges. This is particularly true for surface-based\n"
    "labels or when converting a label to a low-resolution space.\n"
    "To check the result, load the orig volume into tkmedit. The orig\n"
    "volume should be in the label space. Load the mri_label2vol output\n"
    "volume as an overlay; this makes the labeled voxels appear as\n"
    "'activity'.  Finally, load the label itself. You should see the label\n"
    "(in green) sitting on top of the 'activity' of the labeled volume.\n"
    "See EXAMPLE 1 for an example.\n"
    "\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "1. Convert a label into a binary mask in the functional space; require\n"
    "that a functional voxel be filled at least 50%% by the label:\n"
    "\n"
    "mri_label2vol \n"
    "  --label lh-avg_central_sulcus.label \n"
    "  --temp f.nii \n"
    "  --reg register.dat \n"
    "  --fillthresh .5 \n"
    "  --o cent-lh.nii\n"
    "\n"
    "To see how well the label is mapped into the functional volume, run\n"
    "\n"
    "tkmedit bert orig \n"
    "  -overlay ./cent-lh.nii \n"
    "  -overlay-reg ./register.dat -fthresh .5 -fmid 1\n"
    "\n"
    "Then load the label with File->Label->LoadLabel. The label should\n"
    "overlap with the overlay. The overlap will not be perfect but it\n"
    "should be very close.\n"
    "\n"
    "2. Convert a surface label into a binary mask in the functional space.\n"
    "Fill in all the cortical gray matter. Require that a functional voxel\n"
    "be filled at least 30%% by the label:\n"
    "\n"
    "mri_label2vol \n"
    "  --label lh-avg_central_sulcus.label \n"
    "  --temp f.nii \n"
    "  --reg register.dat \n"
    "  --fillthresh .3 \n"
    "  --proj frac 0 1 .1 \n"
    "  --subject bert --hemi lh\n"
    "  --o cent-lh.nii\n"
    "\n"
    "3. Convert a surface label into a binary mask in the functional space.\n"
    "Sample a 1mm ribbon 2mm below the gray/white surface. Do not require a \n"
    "fill threshold:\n"
    "\n"
    "mri_label2vol \n"
    "  --label lh-avg_central_sulcus.label \n"
    "  --temp f.nii \n"
    "  --reg register.dat \n"
    "  --proj abs -3 -2 .1 \n"
    "  --subject bert --hemi lh\n"
    "  --o cent-lh.nii\n"
    "\n"
    "4. Convert two labels into a volume in the same space as the labels:\n"
    "\n"
    "mri_label2vol \n"
    "  --label lh-avg_central_sulcus.label \n"
    "  --label lh-avg_calcarine_sulcus.label \n"
    "  --temp $SUBJECTS_DIR/bert/orig\n"
    "  --identity\n"
    "  --o cent_calc.img\n"
    "\n"
    "The voxels corresponding to lh-avg_central_sulcus.label will have a of \n"
    "value of 1 whereas those assigned to lh-avg_calcarine_sulcus.label will\n"
    "have a value of 2.\n"
    "\n"
    "KNOWN BUGS\n"
    "\n"
    "1. Cannot convert surface labels with different hemispheres.\n"
    "\n"
    "\n"
    "\n"
    "SEE ALSO\n"
    "\n"
    "mri_label2label, mri_cor2label, mri_annotation2label, mri_mergelabels,\n"
    "tkregister2, mri_convert, tkmedit, tksurfer.\n"
    "\n"
    "http://surfer.nmr.mgh.harvard.edu/docs/tkmedit_guide.html\n"
    "http://surfer.nmr.mgh.harvard.edu/docs/tksurfer_doc.html\n"
    "\n"

  );


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
  if(DoProj || FillRibbon) SurfNeeded = 1;

  if(UseAParcPlusASeg && subject == NULL && RegMatFile == NULL){
    printf("ERROR: you must spec --subject or --reg with --aparc+aseg\n");
    exit(1);
  }

  if(nlabels == 0 && AnnotFile == NULL && 
     ASegFSpec == NULL && !UseAParcPlusASeg) {
    printf("ERROR: you must spec at least one label or annot file\n");
    exit(1);
  }
  if(nlabels != 0 && AnnotFile != NULL) {
    printf("ERROR: you cannot specify a label AND an annot file\n");
    exit(1);
  }
  if(nlabels != 0 && ASegFSpec != NULL) {
    printf("ERROR: you cannot specify a label AND an aseg file\n");
    exit(1);
  }
  if(AnnotFile != NULL && ASegFSpec != NULL) {
    printf("ERROR: you cannot specify an annot file AND an aseg file\n");
    exit(1);
  }
  if(OutVolId == NULL) {
    printf("ERROR: no output specified\n");
    exit(1);
  }
  if(SurfNeeded) {
    if(subject == NULL) {
      printf("ERROR: subject needed in order to load surface.\n");
      exit(1);
    }
    if(hemi == NULL) {
      printf("ERROR: hemi needed in order to load surface.\n");
      exit(1);
    }
  }
  if(subject != NULL && !DoProj && AnnotFile==NULL && !UseAParcPlusASeg)
    printf("INFO: subject not needed, igorning.\n");
  if(hemi != NULL && !DoProj && AnnotFile==NULL && !FillRibbon)
    printf("INFO: hemi not needed, igorning.\n");

  if(DoLabelStatVol && nlabels == 0){
    printf("ERROR: %s: must specify a --label with --label-stat\n",Progname);
    exit(1);
  }
  if(DoLabelStatVol && nlabels == 0){
    printf("ERROR: %s: must specify a --label with --stat-thresh\n",Progname);
    exit(1);
  }

  if(!RegIdentity && !RegHeader && !RegMatFile){
    printf("ERROR: you must specify a registration method:\n");
    printf("  Either: --reg, --regheader, or --identity\n");
    exit(1);
  }

  if(RegHeader && RegMatFile){
    printf("ERROR: cannot spec --regheader and --reg\n");
    exit(1);
  }
  if(RegIdentity && RegMatFile){
    printf("ERROR: cannot spec --identity and --reg\n");
    exit(1);
  }
  if(RegIdentity && RegHeader){
    printf("ERROR: cannot spec --identity and --regheader\n");
    exit(1);
  }
  if(SUBJECTS_DIR == NULL){
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL) {
      printf("ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }
  if(TempVolId == NULL && !FillRibbon){
    printf("ERROR: must specify template volume with --temp\n");
    exit(1);
  }
  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  int n;
  fprintf(fp,"Number of labels: %d\n",nlabels);
  if (nlabels > 0) {
    for (n=0;n<nlabels;n++)
      fprintf(fp,"%s\n",LabelList[n]);
  }
  fprintf(fp,"Annot File:      %s\n",AnnotFile);
  fprintf(fp,"Template Volume: %s\n",TempVolId);
  fprintf(fp,"Outut Volume: %s\n",OutVolId);
  fprintf(fp,"Registration File: %s\n",RegMatFile);
  fprintf(fp,"Fill Threshold: %g\n",FillThresh);
  fprintf(fp,"Label Vox Vol:  %g\n",LabelVoxVol);
  fprintf(fp,"ProjType:       %s\n",ProjType);
  fprintf(fp,"ProjTypeId:     %d\n",ProjTypeId);
  fprintf(fp,"ProjStart:      %g\n",ProjStart);
  fprintf(fp,"ProjStop:       %g\n",ProjStop);
  fprintf(fp,"ProjDelta:      %g\n",ProjDelta);
  fprintf(fp,"Subject:  %s\n",subject);
  fprintf(fp,"Hemi:     %s\n",hemi);
  fprintf(fp,"UseNewASeg2Vol:  %d\n",UseNewASeg2Vol);
  fprintf(fp,"DoLabelStatVol  %d\n",DoLabelStatVol);
  fprintf(fp,"LabelCodeOffset  %d\n",LabelCodeOffset);
  fprintf(fp,"setenv SUBJECTS_DIR %s\n",SUBJECTS_DIR);
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
/*------------------------------------------------------------*/
/* Exits if something is wrong with hemi. */
static int checkhemi(char *hemi) {
  if (! strcmp(hemi,"lh")) return(0);
  if (! strcmp(hemi,"rh")) return(0);
  printf("ERROR: hemi = %s, must be lh or rh\n",hemi);
  exit(1);
}
/*------------------------------------------------------------*/
static int get_proj_type_id(char *projtype) {
  if (projtype == NULL)          return(PROJ_TYPE_NONE);
  if (! strcmp(projtype,"none")) return(PROJ_TYPE_NONE);
  if (! strcmp(projtype,"abs"))  return(PROJ_TYPE_ABS);
  if (! strcmp(projtype,"frac")) return(PROJ_TYPE_FRAC);
  printf("ERROR: projtype = %s, must be abs or frac\n",projtype);
  exit(1);
}
/*------------------------------------------------------------*/
// Computes the nearest col, row, slice of the xyz (RAS) point.
// Returns 1 if the point is out of the volume. Should really
// do this in one of the util libraries.
static int get_crs(MATRIX *Tras2vox, float x, float y, float z,
                   int *c, int *r, int *s, MRI *vol) {
  float cf, rf, sf;

  cf = Tras2vox->rptr[1][1]*x +
       Tras2vox->rptr[1][2]*y +
       Tras2vox->rptr[1][3]*z +
       Tras2vox->rptr[1][4];

  rf = Tras2vox->rptr[2][1]*x +
       Tras2vox->rptr[2][2]*y +
       Tras2vox->rptr[2][3]*z +
       Tras2vox->rptr[2][4];

  sf = Tras2vox->rptr[3][1]*x +
       Tras2vox->rptr[3][2]*y +
       Tras2vox->rptr[3][3]*z +
       Tras2vox->rptr[3][4];

  *c = rint(cf);
  *r = rint(rf);
  *s = rint(sf);

  if (vol == NULL) return(0);

  if (*c < 0 || *c >= vol->width)  return(1);
  if (*r < 0 || *r >= vol->height) return(1);
  if (*s < 0 || *s >= vol->depth)  return(1);

  return(0);
}
/*------------------------------------------------------------*/
static int is_surface_label(LABEL *label) {
  int nthpoint,vtxno;

  for (nthpoint = 0; nthpoint < label->n_points; nthpoint++) {
    vtxno = label->lv[nthpoint].vno;
    if (vtxno == -1) return(0);
  }

  return(1);
}
/*------------------------------------------------------------*/
static int load_annotation(char *annotfile, MRIS *Surf) {
  int err, vtxno, annot, annotid, annotidmax;
  VERTEX *vtx;

  printf("Loading annotations from %s\n",annotfile);
  err = MRISreadAnnotation(Surf, annotfile);
  if (err) {
    printf("ERROR: could not load %s\n",annotfile);
    exit(1);
  }

  // Count the number of annotations
  annotidmax = -1;
  for (vtxno = 0; vtxno < Surf->nvertices; vtxno++) {
    vtx = &(Surf->vertices[vtxno]);
    annot = Surf->vertices[vtxno].annotation;
    if (Surf->ct != NULL)
      CTABfindAnnotation(Surf->ct, annot, &annotid);
    else
      annotid = annotation_to_index(annot);
    if (annotidmax < annotid) annotidmax = annotid;
  }
  // add one to account for the 0 id
  annotidmax++;
  printf("annotidmax = %d\n",annotidmax);
  return(annotidmax);
}
/*--------------------------------------------------------------*/
static int *NthLabelMap(MRI *aseg, int *nlabels) {
  int *labelmap, *tmpmap;
  int nthlabel;
  int c,r,s,v;
  int nmax;

  nmax = 1000000;

  tmpmap = (int *) calloc(nmax,sizeof(int));
  for (c=0; c < aseg->width; c++) {
    for (r=0; r < aseg->height; r++) {
      for (s=0; s < aseg->depth; s++) {
        v = (int)MRIgetVoxVal(aseg,c,r,s,0);
        if (v >= 0 && v < nmax-1) tmpmap[v] = 1;
      }
    }
  }

  *nlabels = 0;
  for (v=0;v<nmax-1;v++) if (tmpmap[v]) (*nlabels)++;
  //printf("nlabels = %d\n",*nlabels);

  labelmap = (int *) calloc(*nlabels,sizeof(int));
  nthlabel = 0;
  for (v=0;v<nmax-1;v++) {
    if (tmpmap[v]) {
      //printf("%2d %2d\n",nthlabel,v);
      labelmap[nthlabel] = v;
      nthlabel++;
    }
  }
  free(tmpmap);
  return(labelmap);
}

MRI *MRIsurfaceLabel2VolOpt(MRI *ribbon, MRIS *surf, LABEL **labels, int nlabels, LTA *Q, 
			    int DoStatThresh, double StatThresh, int DoLabelStatVol, 
			    MRI *vollabel)
{
  MRI *overlay;
  int vtxno, nthlabel, nthpoint;
  LABEL *label;
  double val;

  overlay = MRIalloc(surf->nvertices,1,1,MRI_FLOAT);
  for(nthlabel = 0; nthlabel < nlabels; nthlabel++){
    label = labels[nthlabel];
    for (nthpoint = 0; nthpoint < label->n_points; nthpoint++) {
      if(DoStatThresh) 	if(label->lv[nthpoint].stat < StatThresh) continue; 
      vtxno = label->lv[nthpoint].vno;
      if(DoLabelStatVol) val = label->lv[nthpoint].stat;
      else               val = 1.0;
      MRIsetVoxVal(overlay,vtxno,0,0,0,val);
    }
  }  

  vollabel = MRIsurf2VolOpt(ribbon, &surf, &overlay, 1, Q, vollabel);
  MRIfree(&overlay);
  return(vollabel);
}


