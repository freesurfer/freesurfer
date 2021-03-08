/**
 * @brief map a label from one subject to another
 *
 * Purpose: Converts a label in one subject's space to a label
 * in another subject's space using either talairach or spherical
 * as an intermediate registration space.
 *
 *  Example 1: If you have a label from subject fred called
 *   broca-fred.label defined on fred's left hemispherical
 *   surface and you want to convert it to sally's surface, then
 *
 *    mri_label2label --srclabel broca-fred.label  --srcsubject fred
 *                   --trglabel broca-sally.label --trgsubject sally
 *                   --regmethod surface --hemi lh
 *
 *    This will map from fred to sally using sphere.reg. The registration
 *   surface can be changed with --surfreg.
 * 
 *  Example 2: You could also do the same mapping using talairach
 *   space as an intermediate:
 *
 *    mri_label2label --srclabel broca-fred.label  --srcsubject fred
 *                   --trglabel broca-sally.label --trgsubject sally
 *                   --regmethod volume
 *
 *    Note that no hemisphere is specified with -regmethod.
 *
 *  Example 3: You can specify the --usepathfiles flag to read and write
 *   from a tksurfer path file.
 * 
 *    mri_label2label --usepathfiles ...
 *
 *   When mapping from lh to rh:
 *    src reg: lh.sphere.reg
 *    trg reg: rh.lh.sphere.reg
 *
 */
/*
 * Original Author: Douglas Greve
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

#include "icosahedron.h"
#include "fio.h"

#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "mri.h"
#include "label.h"
#include "registerio.h"
#include "mri.h"
#include "mri2.h"
#include "version.h"
#include "path.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
static int  isflag(char *flag);
static int  nth_is_arg(int nargc, char **argv, int nth);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

static int label_erode = 0 ;
static int label_dilate = 0 ;
static int label_open = 0 ;
static int label_close = 0 ;
static int label_ring = 0;

char  *srclabelfile = NULL;
static char  *sample_surf_file = NULL ;
LABEL *srclabel     = NULL;
LABEL *tmplabel     = NULL;
char  *srcsubject   = NULL;
char  *trglabelfile = NULL;
LABEL *trglabel     = NULL;
char  *trgsubject   = NULL;
const char  *trgsurface   = "white";

const char *regmethod  = NULL;
char *hemi       = NULL;
char *srchemi    = NULL;
char *trghemi    = NULL;
const char *surfreg = "sphere.reg";
const char *srcsurfreg = NULL;
const char *trgsurfreg = NULL;
char *srcsurfregfile = NULL; // just spec the file name with hemi
char *trgsurfregfile = NULL;

int srcicoorder = -1;
int trgicoorder = -1;

MRI_SURFACE *SrcSurfReg;
MRI_SURFACE *TrgSurf;
MRI_SURFACE *TrgSurfReg;
MRI_SURFACE *PaintSurf=NULL;
char *PaintSurfName=NULL;
MATRIX *SrcVolReg;
MATRIX *TrgVolReg;
MATRIX *InvTrgVolReg;
MATRIX *Src2TrgVolReg;

float IcoRadius = 100.0;
float hashres = 16;
int usehash = 1;

int debug = 0;

char *SUBJECTS_DIR = NULL;
char *FREESURFER_HOME = NULL;
FILE *fp;

char tmpstr[2000];

char *srcmaskfile, *srcmaskfmt;
const char* srcmasksign = "abs";
int srcmaskframe = 0;
float srcmaskthresh = 0.0;
MRI *SrcMask;

int useprojabs = 0, useprojfrac = 0;
float projabs = 0.0, projfrac = 0.0;
int reversemap = 1;
int usepathfiles = 0;
char *XFMFile = NULL;
char *RegFile = NULL;
int InvertXFM=0;
LTA *lta_transform;

char *OutMaskFile = NULL;
MRI *outmask;

int SrcInv = 0, TrgInv = 0;
int DoPaint = 0;
double PaintMax = 2.0;
int DoRescale = 1;
int DoOutMaskStat = 0;
int DoDirectSurface = 0 ;
static char *src_annot_fname = NULL ;
static char *trg_annot_fname = NULL ;

int UseScannerCoords = 0;
int ToScannerCoords = 0;
MRI *ToScannerTemplate = NULL;
int ToTkrCoords = 0;
MRI *ToTkrTemplate = NULL;

char *DminminFile=NULL;
MRI *mritmp=NULL;

/*-------------------------------------------------*/
int main(int argc, char **argv) {
  int err,m;
  MATRIX *xyzSrc, *xyzTrg;
  MHT *TrgHash, *SrcHash=NULL, *PaintHash=NULL;
  VERTEX *srcvtx, *trgvtx, *trgregvtx;
  struct { float x,y,z; } v;
  int n,srcvtxno,trgvtxno,allzero,nrevhits,srcvtxnominmin;
  float dmin, dminmin, projdist=0.0, dx, dy, dz;
  float SubjRadius, Scale;
  char fname[2000];
  int nSrcLabel, nTrgLabel;
  int nargs;
  int numpathsread;
  PATH** paths;
  PATH* path;

  nargs = handleVersionOption(argc, argv, "mri_label2label");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  printf("\n");

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  /*--- Get environment variables ------*/
  if (SUBJECTS_DIR==NULL) SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL) {
    fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  FREESURFER_HOME = getenv("FREESURFER_HOME") ;
  if (FREESURFER_HOME==NULL) {
    fprintf(stderr,"ERROR: FREESURFER_HOME not defined in environment\n");
    exit(1);
  }
  printf("SUBJECTS_DIR    %s\n",SUBJECTS_DIR);
  printf("FREESURFER_HOME %s\n",FREESURFER_HOME);

  /*--- Load in Source Label ------*/
  if (usepathfiles) {
    printf("INFO: Attempting to read a path file.\n");
    /* Make sure this is a path file. */
    if (!PathIsPathFile(srclabelfile)) {
      fprintf(stderr,"ERROR: %s is not a path file\n",srclabelfile);
      exit(1);
    }
    /* Try to read the path file. */
    err = PathReadMany(srclabelfile, &numpathsread, &paths);
    if (ERROR_NONE!=err) {
      fprintf(stderr,"ERROR reading %s\n",srclabelfile);
      exit(1);
    }
    /* Print a warning if we got more than one. */
    if (numpathsread>0) {
      printf("WARNING: Multiple paths read, only using first one.\n");
    }
    /* Convert the first path. */
    srclabel = NULL;
    err = PathConvertToLabel(paths[0], &srclabel);
    if (ERROR_NONE!=err) {
      fprintf(stderr,"ERROR: Couldn't convert path to label\n");
      exit(1);
    }
  } else {
    printf("Loading source label.\n");
    srclabel = LabelRead(NULL, srclabelfile);
    if (srclabel == NULL) {
      fprintf(stderr,"ERROR reading %s\n",srclabelfile);
      exit(1);
    }
    if(ToScannerCoords){
      if(srclabel->coords != LABEL_COORDS_SCANNER_RAS){
	printf("Converting label to scanner RAS\n");
	LabelToScannerRAS(srclabel, ToScannerTemplate, srclabel);
      }
    }
    if(ToTkrCoords){
      if(srclabel->coords != LABEL_COORDS_TKREG_RAS){
	printf("Converting label to TkReg RAS\n");
	LabelToSurfaceRAS(srclabel, ToTkrTemplate, srclabel);
      }
    }
  }
  printf("Found %d points in source label.\n",srclabel->n_points);
  fflush(stdout);
  fflush(stderr);

  /* Set up vectors */
  xyzSrc = MatrixAlloc(4,1,MATRIX_REAL);
  xyzSrc->rptr[3+1][0+1] = 1.0;
  xyzTrg = MatrixAlloc(4,1,MATRIX_REAL);

  if (DoDirectSurface)
  {
    int     vno ;
    VERTEX  *vtrg ;

    printf("Starting direct surface-based mapping\n");
    if(srcsurfregfile == NULL)
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,
	      srchemi,srcsurfreg);
    else strcpy(tmpstr,srcsurfregfile);
    SrcSurfReg = MRISread(tmpstr);
    if (SrcSurfReg == NULL) {
      fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
      exit(1);
      if (strcmp(srchemi,trghemi)==0)
	sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,
		trghemi,trgsurfreg);
    }

    if (MRISreadAnnotation(SrcSurfReg, src_annot_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read annotation from %s", Progname, src_annot_fname) ;


    if(trgsurfregfile == NULL){
      if (strcmp(srchemi,trghemi)==0)
	sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,
		trghemi,trgsurfreg);
      else
	sprintf(tmpstr,"%s/%s/surf/%s.%s.%s",SUBJECTS_DIR,srcsubject,
		trghemi,srchemi,srcsurfreg);
    }
    else strcpy(tmpstr,trgsurfregfile);
    
    printf("Reading target registration \n %s\n",tmpstr);
    TrgSurfReg = MRISread(tmpstr);
    if (TrgSurfReg == NULL) {
      fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
      exit(1);
    }
    printf("Building source registration hash (res=%g).\n",hashres);
    SrcHash = MHTcreateVertexTable_Resolution(SrcSurfReg, CURRENT_VERTICES,hashres);

    TrgSurfReg->ct = SrcSurfReg->ct ;
    for (vno = 0 ; vno < TrgSurfReg->nvertices ; vno++)
    {
      vtrg = &TrgSurfReg->vertices[vno] ;
      srcvtxno = MHTfindClosestVertexNoXYZ(SrcHash,SrcSurfReg, vtrg->x,vtrg->y,vtrg->z, &dmin);
      if (srcvtxno < 0)
      {
	printf("trg vertex %d could not be mapped!\n", vno) ;
	continue ;
      }
      vtrg->annotation = SrcSurfReg->vertices[srcvtxno].annotation ;
    }
    printf("writing mapped annotation to %s\n", trg_annot_fname);
    MRISwriteAnnotation(TrgSurfReg, trg_annot_fname) ;
    exit(0) ;
  }
  /*--------------------- VOLUMETRIC MAPPING --------------------------*/
  if (!strcmp(regmethod,"volume")) {

    /* -- Allocate the Target Label ---*/
    trglabel = LabelAlloc(srclabel->n_points,trgsubject,trglabelfile);
    trglabel->n_points = srclabel->n_points;

    printf("Starting volumetric mapping %d points\n",trglabel->n_points);

    if (RegFile == NULL) {
      if (XFMFile) {
        printf("Reading in xmf file %s",XFMFile);
        lta_transform = LTAreadEx(XFMFile);
        Src2TrgVolReg = lta_transform->xforms[0].m_L;
      } else {
        /*** Load the Src2Tal registration ***/
        SrcVolReg = DevolveXFM(srcsubject, NULL, NULL);
        if (SrcVolReg == NULL) exit(1);

        /*** Load the Trg2Tal registration ***/
        TrgVolReg = DevolveXFM(trgsubject, NULL, NULL);
        if (TrgVolReg == NULL) exit(1);

        /* Compte the Src-to-Trg Registration */
        InvTrgVolReg = MatrixInverse(TrgVolReg,NULL);
        Src2TrgVolReg = MatrixMultiply(InvTrgVolReg,SrcVolReg,NULL);
      }
    } else {
      char *pc;
      float ipr,bpr, fscale;
      int float2int;
      printf("Reading reg file %s\n",RegFile);
      err = regio_read_register(RegFile, &pc, &ipr, &bpr,
                                &fscale, &Src2TrgVolReg, &float2int);
      if (err) {
        printf("ERROR: reading registration %s\n",RegFile);
        exit(1);
      }
      printf("Inverting reg to make it Src2Trg\n");
      MatrixInverse(Src2TrgVolReg,Src2TrgVolReg);
    }

    printf("Src2TrgVolReg: -----------------\n");
    MatrixPrint(stdout,Src2TrgVolReg);

    if (InvertXFM) {
      printf("Inverting \n");
      MatrixInverse(Src2TrgVolReg,Src2TrgVolReg);
      printf("Inverted Src2TrgVolReg: -----------------\n");
      MatrixPrint(stdout,Src2TrgVolReg);
    }

    /* Loop through each source label and map its xyz to target */
    for (n = 0; n < srclabel->n_points; n++) {

      /* load source label xyz into a vector */
      xyzSrc->rptr[0+1][0+1] = srclabel->lv[n].x;
      xyzSrc->rptr[0+2][0+1] = srclabel->lv[n].y;
      xyzSrc->rptr[0+3][0+1] = srclabel->lv[n].z;

      /* compute xyz location in target space */
      MatrixMultiply(Src2TrgVolReg,xyzSrc,xyzTrg);

      /* unload vector into target label */
      trglabel->lv[n].vno = srclabel->lv[n].vno;
      trglabel->lv[n].x = xyzTrg->rptr[0+1][0+1];
      trglabel->lv[n].y = xyzTrg->rptr[0+2][0+1];
      trglabel->lv[n].z = xyzTrg->rptr[0+3][0+1];
      trglabel->lv[n].stat = srclabel->lv[n].stat;

      if(n<5){
	printf("%3d  %6.4f %6.4f %6.4f    %6.4f %6.4f %6.4f\n",n,
	       srclabel->lv[n].x,srclabel->lv[n].y,srclabel->lv[n].z,
	       trglabel->lv[n].x,trglabel->lv[n].y,trglabel->lv[n].z);
      }
    }

    if(SrcVolReg) MatrixFree(&SrcVolReg) ;
    if(TrgVolReg) MatrixFree(&TrgVolReg) ;
    if(InvTrgVolReg) MatrixFree(&InvTrgVolReg) ;
    MatrixFree(&Src2TrgVolReg) ;
    MatrixFree(&xyzSrc);
    MatrixFree(&xyzTrg);

  }/* done with volumetric mapping */

  /*--------------------- SURFACE-BASED MAPPING --------------------------*/
  if (!strcmp(regmethod,"surface")) {

    printf("Starting surface-based mapping\n");

    /*** Load the source registration surface ***/
    if (strcmp(srcsubject,"ico")) {
      if(srcsurfregfile == NULL)
	sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,
		srchemi,srcsurfreg);
      else strcpy(tmpstr,srcsurfregfile);

      printf("Reading source registration \n %s\n",tmpstr);
      SrcSurfReg = MRISread(tmpstr);
      if (SrcSurfReg == NULL) {
        fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
        exit(1);
      }
      MRISreadWhiteCoordinates(SrcSurfReg, "white") ;
      LabelFillUnassignedVertices(SrcSurfReg, srclabel, WHITE_VERTICES);
      if (DoRescale) {
        printf("Rescaling ... ");
        SubjRadius = MRISaverageRadius(SrcSurfReg) ;
        Scale = IcoRadius / SubjRadius;
        MRISscaleBrain(SrcSurfReg, SrcSurfReg, Scale);
        printf(" original radius = %g\n",SubjRadius);
      }
    } else {
      printf("Reading icosahedron, order = %d, radius = %g\n",
             srcicoorder,IcoRadius);
      SrcSurfReg = ReadIcoByOrder(srcicoorder,IcoRadius);
      if (SrcSurfReg==NULL) {
        printf("ERROR reading icosahedron\n");
        exit(1);
      }
    }

    /*** Load the target surfaces ***/
    if (strcmp(trgsubject,"ico")) {
      /* load target xyz surface */
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,
              trghemi,trgsurface);
      printf("Reading target surface \n %s\n",tmpstr);
      TrgSurf = MRISread(tmpstr);
      if (TrgSurf == NULL) {
        fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
        exit(1);
      }
      /* load target registration surface */
      // same hemi: hemi.sphere.reg
      // diff hemi: trghemi.srchemi.sphere.reg
      // Eg, when mapping from lh to rh: rh.lh.sphere.reg
      if(trgsurfregfile == NULL){
	if (strcmp(srchemi,trghemi)==0)
	  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,
		  trghemi,trgsurfreg);
	else
	  sprintf(tmpstr,"%s/%s/surf/%s.%s.%s",SUBJECTS_DIR,srcsubject,
		  trghemi,srchemi,srcsurfreg);
      }
      else strcpy(tmpstr,trgsurfregfile);

      printf("Reading target registration \n %s\n",tmpstr);
      TrgSurfReg = MRISread(tmpstr);
      if (TrgSurfReg == NULL) {
        fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
        exit(1);
      }
      if(TrgSurf->nvertices != TrgSurfReg->nvertices){
	printf("ERROR: vertex mismatch between target surface and registration\n");
	exit(1);
      }
      if (DoRescale) {
        printf("Rescaling ... ");
        SubjRadius = MRISaverageRadius(TrgSurfReg) ;
        Scale = IcoRadius / SubjRadius;
        MRISscaleBrain(TrgSurfReg, TrgSurfReg, Scale);
        printf(" original radius = %g\n",SubjRadius);
      }
    } else {
      printf("Reading icosahedron, order = %d, radius = %g\n",
             trgicoorder,IcoRadius);
      TrgSurfReg = ReadIcoByOrder(trgicoorder,IcoRadius);
      if (TrgSurfReg==NULL) {
        printf("ERROR reading icosahedron\n");
        exit(1);
      }
      TrgSurf = TrgSurfReg;
    }
    
    if (usehash) {
      printf("Building target registration hash (res=%g).\n",hashres);
      TrgHash = MHTcreateVertexTable_Resolution(TrgSurfReg, CURRENT_VERTICES,hashres);
      printf("Building source registration hash (res=%g).\n",hashres);
      SrcHash = MHTcreateVertexTable_Resolution(SrcSurfReg, CURRENT_VERTICES,hashres);
    }
    if (useprojfrac) {
      sprintf(fname,"%s/%s/surf/%s.thickness",SUBJECTS_DIR,srcsubject,srchemi);
      printf("Reading thickness %s\n",fname);
      MRISreadCurvatureFile(TrgSurf, fname); // is this right?
      printf("Done\n");
    }

    /* handle source mask */
    if (srcmaskfile != NULL) {
      printf("INFO: masking label\n");
      //SrcMask = MRIloadSurfVals(srcmaskfile, srcmaskfmt, NULL,
      //      srcsubject, hemi, NULL);

      SrcMask = MRISloadSurfVals(srcmaskfile, srcmaskfmt, SrcSurfReg,
                                 NULL,NULL,NULL);
      if (SrcMask == NULL) exit(1);
      tmplabel = MaskSurfLabel(srclabel, SrcMask,
                               srcmaskthresh, srcmasksign, srcmaskframe);
      if (tmplabel == NULL) exit(1);
      LabelFree(&srclabel) ;
      srclabel = tmplabel;
      printf("Found %d points in source label after masking.\n",
             srclabel->n_points);
      if (srclabel->n_points == 0) {
        printf("ERROR: no overlap between mask and label\n");
        exit(1);
      }
    }

    /* Invert Source Label */
    if (SrcInv) {
      printf("Inverting source label\n");
      tmplabel = MRISlabelInvert(SrcSurfReg,srclabel);
      LabelFree(&srclabel);
      srclabel = tmplabel;
    }

    /* -- Allocate the Target Label ---*/
    trglabel = LabelAlloc(srclabel->n_points,trgsubject,trglabelfile);
    trglabel->n_points = srclabel->n_points;

    if(DoPaint){
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,
		  srchemi,PaintSurfName);
      printf("Painting onto %s\n",tmpstr);
      PaintSurf = MRISread(tmpstr);
      if (PaintSurf == NULL) {
        printf("ERROR: could not read %s\n",tmpstr);
        exit(1);
      }
      if(usehash)
	PaintHash = MHTcreateVertexTable_Resolution(PaintSurf, CURRENT_VERTICES,hashres);
    }

    /* Loop through each source label and map its xyz to target */
    allzero = 1;
    m = 0;
    dminmin = 10e10;
    srcvtxnominmin = 0;
    for (n = 0; n < srclabel->n_points; n++) {

      /* vertex number of the source label */
      if (DoPaint) {
        v.x = srclabel->lv[n].x;
        v.y = srclabel->lv[n].y;
        v.z = srclabel->lv[n].z;
        if (usehash)
          srcvtxno = MHTfindClosestVertexNoXYZ(PaintHash,PaintSurf,v.x,v.y,v.z,&dmin);
        else
          srcvtxno = MRISfindClosestVertex(PaintSurf,v.x,v.y,v.z,&dmin, CURRENT_VERTICES);
	if(debug) printf("%3d %6d (%5.2f,%5.2f,%5.2f) %g\n",n,srcvtxno,v.x,v.y,v.z,dmin);
        if (dmin > PaintMax) continue;
	if(dmin < dminmin){
	  dminmin = dmin;
	  srcvtxnominmin = srcvtxno;
	}
      } else {
        srcvtxno = srclabel->lv[n].vno;
        if (srcvtxno < 0 || srcvtxno >= SrcSurfReg->nvertices) {
          printf("ERROR: there is a vertex in the label that cannot be \n");
          printf("matched to the surface. This usually occurs when\n");
          printf("the label and surface are from different subjects or \n");
          printf("hemispheres or the surface has been changed since\n");
          printf("the label was created.\n");
          printf("Label point %d: vno = %d, max = %d\n",
                 n,srcvtxno, SrcSurfReg->nvertices);
          exit(1);
        }

      }

      if (srcvtxno != 0) allzero = 0;

      /* source vertex */
      srcvtx = &(SrcSurfReg->vertices[srcvtxno]);

      /* closest target vertex number */
      if (usehash) {
        trgvtxno = MHTfindClosestVertexNo2(TrgHash,TrgSurfReg,SrcSurfReg,srcvtx,&dmin);
        if (trgvtxno < 0) {
          printf("ERROR: trgvtxno = %d < 0\n",trgvtxno);
          printf("srcvtxno = %d, dmin = %g\n",srcvtxno,dmin);
          printf("srcxyz = %g, %g, %g\n",srcvtx->x,srcvtx->y,srcvtx->z);
          exit(1);
        }
      } else {
        trgvtxno = MRISfindClosestVertex(TrgSurfReg,srcvtx->x,srcvtx->y,
                                         srcvtx->z,&dmin, CURRENT_VERTICES);
      }
      /* target vertex */
      trgvtx = &(TrgSurf->vertices[trgvtxno]);

      if (useprojabs || useprojfrac) {
        if (useprojabs)  projdist = projabs;
        if (useprojfrac) projdist = projfrac * trgvtx->curv;
        dx = projdist*trgvtx->nx;
        dy = projdist*trgvtx->ny;
        dz = projdist*trgvtx->nz;
      } else {
        dx = 0.0;
        dy = 0.0;
        dz = 0.0;
      }

      trglabel->lv[m].vno = trgvtxno;
      trglabel->lv[m].x = trgvtx->x + dx;
      trglabel->lv[m].y = trgvtx->y + dy;
      trglabel->lv[m].z = trgvtx->z + dz;
      trglabel->lv[m].stat = srclabel->lv[m].stat;
      m++;
    }
    printf("INFO: found  %d nlabel points\n",m);
    if(DoPaint){
      printf("dminmin = %lf at source vertex %d\n",dminmin,srcvtxnominmin);
      if(DminminFile){
	mritmp = MRIalloc(SrcSurfReg->nvertices,1,1,MRI_INT);
	MRIsetVoxVal(mritmp,srcvtxnominmin,0,0,0, 1);
	err = MRIwrite(mritmp,DminminFile);
	if(err) exit(1);
      }
    }

    /* Do reverse loop here: (1) go through each target vertex
       not already in the label, (2) find closest source vertex,
       (3) determine if source is in the label, (4) if so add
       the target to the label */

    if (reversemap) {
      printf("Performing mapping from target back to the source label %d\n",TrgSurf->nvertices);
      nrevhits = 0;
      for (trgvtxno = 0; trgvtxno < TrgSurf->nvertices; trgvtxno++) {
	trgvtx = &TrgSurf->vertices[trgvtxno] ;
	if(trgvtx->ripflag) continue;

        /* if vertex is already in target label, skip it */
        nTrgLabel = LabelHasVertex(trgvtxno, trglabel);
        if (nTrgLabel != -1) continue;

        trgregvtx = &(TrgSurfReg->vertices[trgvtxno]);
        trgvtx = &(TrgSurf->vertices[trgvtxno]);

        /* Find number of closest source vertex */
        if (usehash) {
          srcvtxno = MHTfindClosestVertexNo2(SrcHash,SrcSurfReg,TrgSurfReg,
                                            trgregvtx,&dmin);
          if (srcvtxno < 0) {
            printf("ERROR: srcvtxno = %d < 0\n",srcvtxno);
            printf("trgvtxno = %d, dmin = %g\n",trgvtxno,dmin);
            printf("trgregxyz = %g, %g, %g\n",
                   trgregvtx->x,trgregvtx->y,trgregvtx->z);
            printf("  This means that a vertex in the target surface could\n");
            printf("  not be mapped to a vertex in the source surface\n");
            printf("  because the xyz of the target is outside of the \n");
            printf("  range of the hash table.\n");
	    srcvtxno = MRISfindClosestVertex(SrcSurfReg,trgregvtx->x,
					     trgregvtx->y,trgregvtx->z,&dmin, CURRENT_VERTICES);
	    printf("dmin = %g\n",dmin);
            exit(1);
          }
        } else {
          srcvtxno = MRISfindClosestVertex(SrcSurfReg,trgregvtx->x,
                                           trgregvtx->y,trgregvtx->z,&dmin, CURRENT_VERTICES);
        }
        srcvtx = &(SrcSurfReg->vertices[srcvtxno]);

        /* Determine whether src vtx is in the label */
        nSrcLabel = LabelHasVertex(srcvtxno, srclabel);
        if (nSrcLabel == -1) continue;

        /* Compute dist to project along normal */
        if (useprojabs || useprojfrac) {
          if (useprojabs)  projdist = projabs;
          if (useprojfrac) projdist = projfrac * trgvtx->curv;
          dx = projdist*trgvtx->nx;
          dy = projdist*trgvtx->ny;
          dz = projdist*trgvtx->nz;
        } else {
          dx = 0.0;
          dy = 0.0;
          dz = 0.0;
        }

        /* Alloc another vertex to the label */
        LabelRealloc(trglabel, trglabel->n_points + 1);
        nTrgLabel = trglabel->n_points;
        trglabel->lv[nTrgLabel].vno = trgvtxno;
        trglabel->lv[nTrgLabel].x = trgvtx->x + dx;
        trglabel->lv[nTrgLabel].y = trgvtx->y + dy;
        trglabel->lv[nTrgLabel].z = trgvtx->z + dz;
        trglabel->lv[nTrgLabel].stat = srclabel->lv[nSrcLabel].stat;
        trglabel->n_points ++;

        if (trgvtxno == 53018 && 0) {
          printf("trgvtxno = %d\n",trgvtxno);
          printf("vtx xyz = %g, %g, %g\n",trgvtx->x,trgvtx->y,trgvtx->z);
          printf("dx = %g, dy = %g, dz = %g\n",dx,dy,dz);
        }

        nrevhits++;
      }
      printf("Number of reverse mapping hits = %d\n",nrevhits);
      if (usehash) MHTfree(&SrcHash);
    }

    if (allzero) {
      printf("---------------------------------------------\n");
      printf("WARNING: all source vertex numbers were zero.\n");
      printf("Make sure that the source label is surface-based.\n");
      printf("---------------------------------------------\n");
    }

    printf("Checking for and removing duplicates\n");
    // Does not actually remove them, just flags them
    LabelRemoveDuplicates(trglabel);

    /* Invert Targ Label */
    if (TrgInv) {
      printf("Inverting target label\n");
      tmplabel = MRISlabelInvert(TrgSurfReg,trglabel);
      LabelFree(&trglabel);
      trglabel = tmplabel;
    }

    if (label_ring > 0)
    {
      LABEL *label_interior;
	
      label_interior = LabelCopy(trglabel, NULL) ;
      LabelDilate(trglabel, TrgSurf, label_ring, CURRENT_VERTICES) ;
      LabelRemoveOverlap(trglabel, label_interior) ;
    }
    if (label_dilate)
      LabelDilate(trglabel, TrgSurf, label_dilate, CURRENT_VERTICES) ;
    if (label_erode)
      LabelErode(trglabel, TrgSurf, label_erode) ;
    if (label_close)
    {
      LabelDilate(trglabel, TrgSurf, label_close, CURRENT_VERTICES) ;
      LabelErode(trglabel, TrgSurf, label_close) ;
    }
    if (label_open)
    {
      LabelErode(trglabel, TrgSurf, label_open) ;
      LabelDilate(trglabel, TrgSurf, label_open, CURRENT_VERTICES) ;
    }

    if (OutMaskFile) {
      printf("Creating output %s\n",OutMaskFile);
      outmask = MRISlabel2Mask(TrgSurfReg,trglabel,NULL);
      if (DoOutMaskStat) {
        printf("Saving output statistic\n");
        for (n = 0; n < trglabel->n_points; n++) {
          MRIsetVoxVal(outmask, 
                       trglabel->lv[n].vno,
                       0,0,0,
                       trglabel->lv[n].stat);
        }
      }
      MRIwrite(outmask,OutMaskFile);
      MRIfree(&outmask);
    }

    MRISfree(&SrcSurfReg);
    MRISfree(&TrgSurfReg);
    if (usehash) MHTfree(&TrgHash);
    if (strcmp(trgsubject,"ico")) MRISfree(&TrgSurf);

  }/*---------- done with surface-based mapping -------------*/


  if(UseScannerCoords) strcpy(trglabel->space,"scanner");

  if (usepathfiles) {
    /* Convert the label to a path. */
    err = PathCreateFromLabel(trglabel,&path);
    if (ERROR_NONE!=err) {
      fprintf(stderr,"ERROR: Couldn't convert label to path\n");
      exit(1);
    }
    /* Set the first path in the array. */
    PathFree(&paths[0]);
    paths[0] = path;
    /* Write the path file. */
    printf("Writing path file %s \n",trglabelfile);
    err = PathWriteMany(trglabelfile, 1, paths);
    if (ERROR_NONE!=err) {
      fprintf(stderr,"ERROR writing %s\n",trglabelfile);
      exit(1);
    }
  } else {
    if (sample_surf_file) {
      MRI_SURFACE *mris ;
      mris = MRISread(sample_surf_file) ;
      if (mris == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not load sampling surface %s",
                  sample_surf_file) ;
      printf("sampling label onto surface %s...\n", sample_surf_file) ;
      LabelUnassign(trglabel) ;
      LabelFillUnassignedVertices(mris, trglabel, CURRENT_VERTICES);
      MRISfree(&mris) ;
    }
    printf("Writing label file %s %d\n",trglabelfile,trglabel->n_points);
    if (LabelWrite(trglabel,trglabelfile))
      printf("ERROR: writing label file\n");
  }

  printf("mri_label2label: Done\n\n");

  return(0);
}
/* --------------------------------------------- */
/* --------------------------------------------- */
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
    else if (!strcasecmp(option, "--hash"))   usehash = 1;
    else if (!strcasecmp(option, "--nohash")) usehash = 0;
    else if (!strcasecmp(option, "--norevmap")) reversemap = 0;
    else if (!strcasecmp(option, "--revmap")) reversemap = 1;
    else if (!strcasecmp(option, "--usepathfiles")) usepathfiles = 1;
    else if (!strcmp(option, "--xfm-invert")) InvertXFM = 1;
    else if (!strcmp(option, "--src-invert")) SrcInv = 1;
    else if (!strcmp(option, "--direct")) 
    {
      if (nargc < 2) argnerr(option,1);
      DoDirectSurface = 1 ;
      src_annot_fname = pargv[0] ;
      trg_annot_fname = pargv[1] ;
      nargsused = 2 ;
    }
    else if (!strcmp(option, "--trg-invert")) TrgInv = 1;
    else if (!strcmp(option, "--scanner")) UseScannerCoords = 1;
    else if (!strcmp(option, "--to-scanner")){
      if(nargc < 1) argnerr(option,1);
      ToScannerTemplate = MRIreadHeader(pargv[0],MRI_VOLUME_TYPE_UNKNOWN);
      if(ToScannerTemplate==NULL) exit(1);
      nargsused = 1;
    }
    else if (!strcmp(option, "--to-tkr")){
      if (nargc < 1) argnerr(option,1);
      ToTkrCoords = 1;
      ToTkrTemplate = MRIreadHeader(pargv[0],MRI_VOLUME_TYPE_UNKNOWN);
      if(ToTkrTemplate==NULL) exit(1);
      nargsused = 1;
    }

    else if (!strcmp(option, "--s")) {
      if (nargc < 1) argnerr(option,1);
      srcsubject = pargv[0];
      trgsubject = pargv[0];
      nargsused = 1;
    }
    /* -------- source inputs ------ */
    else if (!strcmp(option, "--sd")) {
      if (nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--dilate")) {
      if (nargc < 1) argnerr(option,1);
      label_dilate = atoi(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--ring")) {
      if (nargc < 1) argnerr(option,1);
      label_ring = atoi(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--erode")) {
      if (nargc < 1) argnerr(option,1);
      label_erode = atoi(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--open")) {
      if (nargc < 1) argnerr(option,1);
      label_open = atoi(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--close")) {
      if (nargc < 1) argnerr(option,1);
      label_close = atoi(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--srcsubject")) {
      if (nargc < 1) argnerr(option,1);
      srcsubject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srclabel")) {
      if (nargc < 1) argnerr(option,1);
      srclabelfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--sample")) {
      if (nargc < 1) argnerr(option,1);
      sample_surf_file = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srcicoorder")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&srcicoorder);
      nargsused = 1;
    } else if (!strcmp(option, "--srcmask")) {
      if (nargc < 2) argnerr(option,2);
      srcmaskfile = pargv[0];
      sscanf(pargv[1],"%f",&srcmaskthresh);
      nargsused = 2;
      if (nth_is_arg(nargc, pargv, 2)) {
        srcmaskfmt = pargv[2];
        nargsused ++;
      }
    } else if (!strcmp(option, "--srcmaskframe")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&srcmaskframe);
      nargsused = 1;
    } else if (!strcmp(option, "--srcmasksign")) {
      if (nargc < 1) argnerr(option,1);
      srcmasksign = pargv[0];
      nargsused = 1;
      if (strcmp(srcmasksign,"abs") &&
          strcmp(srcmasksign,"pos") &&
          strcmp(srcmasksign,"neg")) {
        printf("ERROR: srcmasksign = %s, must be either "
               "abs, pos, or neg\n", srcmasksign);
        exit(1);
      }
    }
    /* -------- target inputs ------ */
    else if (!strcmp(option, "--trgsubject")) {
      if (nargc < 1) argnerr(option,1);
      trgsubject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trglabel")) {
      if (nargc < 1) argnerr(option,1);
      trglabelfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trgsurface") ||
               !strcmp(option, "--trgsurf")) {
      if (nargc < 1) argnerr(option,1);
      trgsurface = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--surfreg")) {
      if (nargc < 1) argnerr(option,1);
      surfreg = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--srcsurfreg")) {
      if (nargc < 1) argnerr(option,1);
      srcsurfreg = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--trgsurfreg")) {
      if (nargc < 1) argnerr(option,1);
      trgsurfreg = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--srcsurfreg-file")) {
      if (nargc < 1) argnerr(option,1);
      srcsurfregfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--trgsurfreg-file")) {
      if (nargc < 1) argnerr(option,1);
      trgsurfregfile = pargv[0];
      nargsused = 1;
    } 

    else if (!strcmp(option, "--trgicoorder")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&trgicoorder);
      nargsused = 1;
    } else if (!strcmp(option, "--hashres")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&hashres);
      nargsused = 1;
    } else if (!strcmp(option, "--projabs")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&projabs);
      useprojabs = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--projfrac")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&projfrac);
      useprojfrac = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--hemi")) {
      if (nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srchemi")) {
      if (nargc < 1) argnerr(option,1);
      srchemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trghemi")) {
      if (nargc < 1) argnerr(option,1);
      trghemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--regmethod")) {
      if (nargc < 1) argnerr(option,1);
      regmethod = pargv[0];
      if (strcmp(regmethod,"surface") && strcmp(regmethod,"volume") &&
          strcmp(regmethod,"surf") && strcmp(regmethod,"vol")) {
        fprintf(stderr,"ERROR: regmethod must be surface or volume\n");
        exit(1);
      }
      if (!strcmp(regmethod,"surf")) regmethod = "surface";
      if (!strcmp(regmethod,"vol"))  regmethod = "volume";
      nargsused = 1;
    } 
    else if (!strcmp(option, "--paint")) {
      if (nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%lf",&PaintMax);
      DoPaint = 1;
      DoRescale = 0;
      reversemap = 0;
      regmethod = "surface";
      PaintSurfName = pargv[1];
      nargsused = 2;
    } 
    else if (!strcmp(option, "--dminmin")) {
      if (nargc < 1) argnerr(option,1);
      DminminFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--xfm")) {
      if (nargc < 1) argnerr(option,1);
      XFMFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--reg")) {
      if (nargc < 1) argnerr(option,1);
      RegFile = pargv[0];
      regmethod = "volume";
      nargsused = 1;
    } 
    else if (!strcmp(option, "--outmask")) {
      if (nargc < 1) argnerr(option,1);
      OutMaskFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--surf-label2mask")) {
      if(nargc != 3) {
	printf("--surf-label2mask label surf mask\n");
	exit(1);
      }
      srclabel = LabelRead(NULL, pargv[0]);
      if(!srclabel) exit(1);
      SrcSurfReg = MRISread(pargv[1]);
      if(!SrcSurfReg) exit(1);
      outmask = MRISlabel2Mask(SrcSurfReg,srclabel,NULL);
      int err = MRIwrite(outmask,pargv[2]);
      exit(err);
      nargsused = 3;
    } 
    else if (!strcmp(option, "--outstat")) {
      if (nargc < 1) argnerr(option,1);
      OutMaskFile = pargv[0];
      DoOutMaskStat = 1;
      nargsused = 1;
    } 
    else if(!strcmp(option, "--Gdiag_no")){
      if (nargc < 1) argnerr(option,1);
      Gdiag_no = atoi(pargv[0]) ;
      printf("Gdiag_no set to %d\n",Gdiag_no);
      nargsused = 1;
    }
    else if (!strcmp(option, "--label-cortex")) {
      // surf aseg outlabel
      if(nargc < 4){
	printf("ERROR: when using --label-cortex, the usage is:\n");
	printf("  --label-cortex surf aseg KeepHipAmyg01 outlabel \n");
	exit(1);
      }
      MRIS *lsurf = MRISread(pargv[0]);
      if(lsurf==NULL) exit(1);
      MRI *aseg = MRIread(pargv[1]);
      if(aseg==NULL) exit(1);
      MRIScomputeMetricProperties(lsurf);// might not be needed
      int KeepHipAmyg01;
      sscanf(pargv[2],"%d",&KeepHipAmyg01);
      LABEL *lcortex = MRIScortexLabelDECC(lsurf, aseg, 4, 4, -1, KeepHipAmyg01);
      if(lcortex == NULL) exit(1);
      int err = LabelWrite(lcortex, pargv[3]);
      exit(err);
      nargsused = 4;
    } 
    else if (!strcmp(option, "--baryfill")) {
      if (nargc < 4) argnerr(option,4);
      MRIS *mris;
      LABEL *lab,*outlab;
      double delta;
      mris = MRISread(pargv[0]);
      if(mris==NULL) exit(1);
      lab = LabelRead(NULL,pargv[1]);
      if(lab==NULL) exit(1);
      sscanf(pargv[2],"%lf",&delta);
      printf("delta = %lf\n",delta);
      outlab = LabelBaryFill(mris, lab, delta);
      printf("writing to %s\n",pargv[3]);
      LabelWrite(outlab, pargv[3]);
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
  printf("   --srclabel     input label file \n");
  printf("   --trglabel     output label file \n");
  printf("\n");
  printf("   --erode  N     erode the label N times before writing\n");
  printf("   --open   N     open the label N times before writing\n");
  printf("   --close  N     close the label N times before writing\n");
  printf("   --dilate N    dilate the label N times before writing\n");
  printf("   --ring   N    dilate the label N times then remove the original before writing\n");
  printf("   --srcsubject   source subject\n");
  printf("   --trgsubject   target subject\n");
  printf("   --s subject : use for both target and source\n");
  printf("\n");
  printf("   --outmask      maskfile : save output label as a "
         "binary mask (surf only)\n");
  printf("   --outstat      statfile : save output label stat as a "
         "mask (surf only)\n");
  printf("   --sample       output subject surface : sample label "
         "onto surface \n");
  printf("\n");
  printf("   --regmethod    registration method (surface, volume) \n");
  printf("   --usepathfiles read from and write to a path file\n");
  printf("\n");
  printf("   --hemi        hemisphere (lh or rh) (with surface)\n");
  printf("   --srchemi     hemisphere (lh or rh) (with surface)\n");
  printf("   --trghemi     hemisphere (lh or rh) (with surface)\n");
  printf("   --srcicoorder when srcsubject=ico\n");
  printf("   --trgicoorder when trgsubject=ico\n");
  printf("   --direct <src annot> <trg annot>     use the [xyz] coords for src and trg surfaces to do direct lookup\n");
  printf("   --trgsurf     get xyz from this surface (white)\n");
  printf("   --surfreg     surface registration (sphere.reg)  \n");
  printf("   --srcsurfreg  source surface registration (sphere.reg)\n");
  printf("   --trgsurfreg  target surface registration (sphere.reg)\n");
  printf("   --srcsurfreg-file  specify full path to source reg\n");
  printf("   --trgsurfreg-file  specify full path to source reg\n");

  printf("\n");
  printf("   --paint dmax surfname : map to closest vertex on source surfname if d < dmax\n");
  printf("   --dmindmin overlayfile : bin mask with vertex of closest label point when painting\n");
  printf("   --baryfill surf surflabel delta outlabel\n");
  printf("   --label-cortex surface aseg KeepHipAmyg01 outlabel : create a label like ?h.cortex.label\n");
  printf("   --surf-label2mask label surf mask : stand-alone way to convert a label to a binary mask\n");
  printf("\n");
  printf("   --srcmask     surfvalfile thresh <format>\n");
  printf("   --srcmasksign sign (<abs>,pos,neg)\n");
  printf("   --srcmaskframe 0-based frame number <0>\n");
  printf("\n");
  printf("   --xfm xfmfile : use xfm instead of computing tal xfm\n");
  printf("   --reg regfile : use register.dat file instead of computing "
         "tal xfm\n");
  printf("   --xfm-invert : invert xfm, or reg \n");
  printf("\n");
  printf("   --projabs  dist project dist mm along surf normal\n");
  printf("   --projfrac frac project frac of thickness along surf normal\n");
  printf("\n");
  printf("   --sd subjectsdir : default is to use env SUBJECTS_DIR\n");
  printf("   --nohash : don't use hash table when regmethod is surface\n");
  printf("   --norevmap : don't use reverse mapping regmethod is surface\n");
  printf("   --to-scanner template : convert coords to scanner RAS (if needed) prior \n");
  printf("       to other operations. template is the MRI volume that the label was created on\n");
  printf("   --to-tkr template : convert coords to tkregister RAS (if needed) prior \n");
  printf("       to other operations. template is the MRI volume that the label was created on\n");
  printf("   --scanner : set output coordinate type to scanner\n");
  printf("       NOTE: --scanner does nothing more than change a string in the label file\n");
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf(
    "  Purpose: Converts a label in one subject's space to a label\n"
    "  in another subject's space using either talairach or spherical\n"
    "  as an intermediate registration space. \n"
    "\n"
    "  If a source mask is used, then the input label must have been\n"
    "  created from a surface (ie, the vertex numbers are valid). The \n"
    "  format can be anything supported by mri_convert or curv or paint.\n"
    "  Vertices in the source label that do not meet threshold in the\n"
    "  mask will be removed from the label. See Example 2.\n"
    "\n"
    "  Example 1: If you have a label from subject fred called\n"
    "    broca-fred.label defined on fred's left hemispherical \n"
    "    surface and you want to convert it to sally's surface, then\n"
    "\n"
    "    mri_label2label --srclabel broca-fred.label  --srcsubject fred \n"
    "                    --trglabel broca-sally.label --trgsubject sally\n"
    "                    --regmethod surface --hemi lh\n"
    "\n"
    "    This will map from fred to sally using sphere.reg. The registration\n"
    "    surface can be changed with --surfreg.\n"
    "\n"
    "  Example 2: Same as Example 1 but with a mask\n"
    "\n"
    "    mri_label2label --srclabel broca-fred.label  --srcsubject fred \n"
    "                    --trglabel broca-sally.label --trgsubject sally\n"
    "                    --regmethod surface --hemi lh\n"
    "                    --srcmask  fred-omnibus-sig 2 bfloat\n"
    "\n"
    "    This will load the bfloat data from fred-omnibus-sig and create\n"
    "    a mask by thresholding the first frame absolute values at 2.\n"
    "    To change it to only the positive values of the 3rd frame, add\n"
    "         --srcmasksign pos --srcmaskframe 2   \n"
    "\n"
    "\n"
    "  Example 3: You could also do the same mapping using talairach \n"
    "    space as an intermediate:\n"
    "\n"
    "    mri_label2label --srclabel broca-fred.label  --srcsubject fred \n"
    "                    --trglabel broca-sally.label --trgsubject sally\n"
    "                    --regmethod volume\n"
    "\n"
    "    Note that no hemisphere is specified with --regmethod volume.\n"
    "\n"
    "  Example 4: You have a label in the volume and you want to find \n"
    "  the closest surface vertices:\n"
    "\n"
    "   mri_label2label --srclabel your.volume.label --s subject \n"
    "     --trglabel lh.your.volume.on-pial.label --hemi lh --paint 30 pial\n"
    "     --trgsurf pial\n"
    "  This keeps the label on a single subject (but could also map to \n"
    "  another subject). The label is mapped to vertices on the pial surface\n"
    "  that are within 30mm of the label point. The xyz of the output label\n"
    "  takes the coordinates of the pial surface (--trgsurf pial).\n"
    "\n"
    "  Notes:\n"
    "\n"
    "  1. A label can be converted to/from talairach space by specifying\n"
    "     the target/source subject as 'talairach'.\n"
    "  2. A label can be converted to/from the icosahedron by specifying\n"
    "     the target/source subject as 'ico'. When the source or target\n"
    "     subject is specified as 'ico', then the order of the icosahedron\n"
    "     must be specified with --srcicoorder/--trgicoorder.\n"
    "  3. When the surface registration method is used, the xyz coordinates\n"
    "     in the target label file are derived from the xyz coordinates\n"
    "     from the target subject's white surface. This can be changed\n"
    "     using the --trgsurf option.\n"
    "  4. When the volume registration method is used, the xyz coordinates\n"
    "     in the target label file are computed as xyzTrg = "
    "inv(Ttrg)*Tsrc*xyzSrc\n"
    "     where Tsrc is the talairach transform in \n"
    "     srcsubject/mri/transforms/talairach.xfm, and where Ttrg "
    "is the talairach \n"
    "     transform in trgsubject/mri/transforms/talairach.xfm.\n"
    "  5. The registration surfaces are rescaled to a radius of 100 "
    "(including \n"
    "     the ico)\n"
    "  6. Projections along the surface normal can be either negative or\n"
    "     positive, but can only be used with surface registration method.\n"
    "\n"
    "BUGS:\n"
    "\n"
    "When using volume registration method, you cannot specify the "
    "SUBJECTS_DIR\n"
    "on the command-line.\n"
    "\n"
  );

  exit(1) ;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"srclabel = %s\n",  srclabelfile);
  fprintf(fp,"srcsubject = %s\n",srcsubject);
  fprintf(fp,"trgsubject = %s\n",trgsubject);
  fprintf(fp,"trglabel = %s\n",  trglabelfile);
  fprintf(fp,"regmethod = %s\n",regmethod);
  fprintf(fp,"\n");
  if (!strcmp(regmethod,"surface")) {
    fprintf(fp,"srchemi = %s\n",srchemi);
    fprintf(fp,"trghemi = %s\n",trghemi);
    fprintf(fp,"trgsurface = %s\n",trgsurface);
    if(srcsurfregfile == NULL)
      fprintf(fp,"srcsurfreg = %s\n",srcsurfreg);
    else
      fprintf(fp,"srcsurfregfile = %s\n",srcsurfregfile);
    if(trgsurfregfile == NULL)
      fprintf(fp,"trgsurfreg = %s\n",trgsurfreg);
    else
      fprintf(fp,"trgsurfregfile = %s\n",trgsurfregfile);
  }
  if (!strcmp(srcsubject,"ico")) fprintf(fp,"srcicoorder = %d\n",srcicoorder);
  if (!strcmp(trgsubject,"ico")) fprintf(fp,"trgicoorder = %d\n",trgicoorder);
  fprintf(fp,"usehash = %d\n",usehash);

  if (srcmaskfile != NULL) {
    fprintf(fp,"srcmask %s, %s \n",srcmaskfile, srcmaskfmt);
    fprintf(fp,"srcmaskthresh %g %s\n",srcmaskthresh, srcmasksign);
    fprintf(fp,"srcmaskframe %d\n",srcmaskframe);
  }
  printf("Use ProjAbs  = %d, %g\n",useprojabs,projabs);
  printf("Use ProjFrac = %d, %g\n",useprojfrac,projfrac);
  printf("DoPaint %d\n",DoPaint);
  if (DoPaint)  printf("PaintMax %lf\n",PaintMax);

  fprintf(fp,"\n");

  return;
}
/* --------------------------------------------- */
static void print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
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
/*---------------------------------------------------------------*/
static int isflag(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth) {
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if (nargc <= nth) return(0);

  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) return(0);

  return(1);
}
/* --------------------------------------------- */
static void check_options(void) {
  if (srcsubject == NULL) {
    fprintf(stderr,"ERROR: no source subject specified\n");
    exit(1);
  }
  if (srclabelfile == NULL) {
    fprintf(stderr,"ERROR: A source label path must be supplied\n");
    exit(1);
  }
  if (trgsubject == NULL) {
    fprintf(stderr,"ERROR: no target subject specified\n");
    exit(1);
  }
  if (trglabelfile == NULL) {
    fprintf(stderr,"ERROR: A target label path must be supplied\n");
    exit(1);
  }

  if (regmethod == NULL) {
    fprintf(stderr,"ERROR: Must specify a registration method\n");
    exit(1);
  }

  if (!strcmp(regmethod,"surface")) {
    if (srchemi == NULL && trghemi == NULL && hemi == NULL) {
      fprintf(stderr,"ERROR: no hemisphere specified\n");
      exit(1);
    }
    if ((srchemi == NULL && trghemi != NULL) ||
        (srchemi != NULL && trghemi == NULL) ) {
      fprintf(stderr,"ERROR: must specify either --hemi or "
              "both --srchemi and --trghemi\n");
      exit(1);
    }
    if (srchemi == NULL) srchemi = hemi;
    if (trghemi == NULL) trghemi = hemi;
  } else { /* volume */
    if (!strcmp(srcsubject,"ico") || !strcmp(trgsubject,"ico")) {
      fprintf(stderr,"ERROR: cannot use volume registration "
              "method with subject ico\n");
      exit(1);
    }
    if (hemi != NULL) {
      fprintf(stderr,"ERROR: cannot specify hemisphere with vol reg method\n");
      exit(1);
    }
    if (OutMaskFile) {
      printf("ERROR: cannot specify outmask with vol reg method\n");
      exit(1);
    }
    if (SrcInv) {
      printf("ERROR: cannot specify src-invert with vol reg method\n");
      exit(1);
    }
    if (TrgInv) {
      printf("ERROR: cannot specify trg-invert with vol reg method\n");
      exit(1);
    }
  }

  if (!strcmp(srcsubject,"ico") && srcicoorder < 0) {
    fprintf(stderr,"ERROR: must specify src ico order with srcsubject=ico\n");
    exit(1);
  }

  if (!strcmp(trgsubject,"ico") && trgicoorder < 0) {
    fprintf(stderr,"ERROR: must specify trg ico order with trgsubject=ico\n");
    exit(1);
  }

  if (!strcmp(regmethod,"surface") && (!strcmp(srcsubject,"talairach") ||
                                       !strcmp(trgsubject,"talairach"))) {
    fprintf(stderr,"ERROR: cannot use talairach with surface mapping\n");
    exit(1);
  }

  if (useprojabs && useprojfrac) {
    fprintf(stderr,"ERROR: cannot use absolute and fractional projection\n");
    exit(1);
  }

  if ( (useprojabs || useprojfrac) &&  strcmp(regmethod,"surface") ) {
    fprintf(stderr,"ERROR: must use surface regmethod with absolute "
            "or fractional projection\n");
    exit(1);
  }

  if (srcsurfreg == NULL) srcsurfreg = surfreg;
  if (trgsurfreg == NULL) trgsurfreg = surfreg;

  if (XFMFile && RegFile) {
    printf("ERROR: cannot --xfm and --reg\n");
    exit(1);
  }

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
