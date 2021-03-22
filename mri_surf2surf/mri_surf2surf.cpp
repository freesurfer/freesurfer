/**
 * @brief Resamples data from one surface onto another.
 *
 * Purpose: Resamples data from one surface onto another. If
 * both the source and target subjects are the same, this is
 * just a format conversion. The source or target subject may
 * be ico.  Can handle data with multiple frames.
 */
/*
 * Original Author: Douglas Greve
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


/*
BEGINHELP

This program will resample one surface onto another. The source and
target subjects can be any subject in $SUBJECTS_DIR and/or the
icosahedron (ico). The source and target file formats can be anything
supported by mri_convert. The source format can also be a curvature
file or a paint (.w) file. The user also has the option of smoothing
on the surface.

OPTIONS

  --srcsubject subjectname

    Name of source subject as found in $SUBJECTS_DIR or ico for icosahedron.
    The input data must have been sampled onto this subject's surface (eg,
    using mri_vol2surf)

  --sval sourcefile

    Name of file where the data on the source surface is located.

  --sval-xyz     surfname
  --sval-tal-xyz surfname
  --sval-area    surfname
  --sval-nxyz    surfname

    Use measures from the input surface as the source (instead of specifying
    a source file explicitly with --sval). --sval-xyz extracts the x, y, and
    z of each vertex. --sval-tal-xyz is the same as --sval-xyz, but applies
    the talairach transform from srcsubject/mri/transforms/talairach.xfm.
    --sval-area extracts the vertex area. --sval-nxyz extracts the surface
    normals at each vertex. See also --tval-xyz.

  --projfrac surfname frac
  --projabs  surfname dist

    Use xyz from surfname as the input, project it along the normal, and
    save new xyz surface. Eg, to create a new surface halfway between
    the white and the pial:
      mri_surf2surf --s subject --projfrac white +0.5 --tval lh.mid --hemi lh
    saves $SUBJECTS_DIR/subject/lh.mid.

  --sval-annot annotfile

    Map annotation file to the output. The target data will be saved as an
    annotation.

  --sfmt typestring

    Format type string. Can be either curv (for FreeSurfer curvature file),
    paint or w (for FreeSurfer paint files), or anything accepted by
    mri_convert. If no type string  is given, then the type is determined
    from the sourcefile (if possible). If curv is used, then the curvature
    file will be looked for in $SUBJECTS_DIR/srcsubject/surf/hemi.sourcefile.

  --srcicoorder order

    Icosahedron order of the source. Normally, this can be detected based
    on the number of verticies, but this will fail with a .w file as input.
    This is only needed when the source is a .w file.

  --trgsubject subjectname

    Name of target subject as found in $SUBJECTS_DIR or ico for icosahedron.

  --trgicoorder order

    Icosahedron order number. This specifies the size of the
    icosahedron according to the following table:
              Order  Number of Vertices
                0              12
                1              42
                2             162
                3             642
                4            2562
                5           10242
                6           40962
                7          163842
    In general, it is best to use the largest size available.

  --tval targetfile

    Name of file where the data on the target surface will be stored.
    BUG ALERT: for trg_type w or paint, use the full path.

  --tval-xyz volume

    Use this flag to indicate that the output (specified by --tval)
    will be a binary surface file This requires that --sval-xyz or
    --sval-tal-xyz was also specified. volume is a volume file in
    the target space; the volume geometry from this file is imbedded
    into the surface file. This is a good way to map the surface of
    one subject to an average (talairach) subject. Note: it will save
    targetfile as trgsubject/surf/targetfile unless targetfile has a
    path.

  --tfmt typestring

    Format type string. Can be paint or w (for FreeSurfer paint files) or curv
    or anything accepted by mri_convert. If no type string  is given, then the type
    is determined from the sourcefile (if possible). If using paint, w, or curv,
    see also --frame.

  --hemi hemifield (lh or rh)

  --surfreg registration_surface

    If the source and target subjects are not the same, this surface is used
    to register the two surfaces. sphere.reg is used as the default. Don't change
    this unless you know what you are doing.

  --mapmethod methodname

    Method used to map from the vertices in one subject to those of another.
    Legal values are: nnfr (neighest-neighbor, forward and reverse) and nnf
    (neighest-neighbor, forward only). Default is nnfr. The mapping is done
    in the following way. For each vertex on the target surface, the closest
    vertex in the source surface is found, based on the distance in the
    registration space (this is the forward map). If nnf is chosen, then the
    the value at the target vertex is set to that of the closest source vertex.
    This, however, can leave some source vertices unrepresented in target (ie,
    'holes'). If nnfr is chosen, then each hole is assigned to the closest
    target vertex. If a target vertex has multiple source vertices, then the
    source values are averaged together. It does not seem to make much difference.

  --fwhm-src fwhmsrc
  --fwhm-trg fwhmtrg (can also use --fwhm)

    Smooth the source or target with a gaussian with the given fwhm (mm). This is
    actually an approximation done using iterative nearest neighbor smoothing.
    The number of iterations is computed based on the white surface. This
    method is similar to heat kernel smoothing. This will give the same
    results as --nsmooth-{in,out}, but automatically computes the the
    number of iterations based on the desired fwhm.

  --nsmooth-in  niterations
  --nsmooth-out niterations  [note: same as --smooth]

    Number of smoothing iterations. Each iteration consists of averaging each
    vertex with its neighbors. When only smoothing is desired, just set the
    the source and target subjects to the same subject. --smooth-in smooths
    the input surface values prior to any resampling. --smooth-out smooths
    after any resampling. See also --fwhm-src and --fwhm-trg.

  --label-src sourcelabel
  --label-trg targetlabel
  --cortex
  --no-cortex

    Only smooth within the given label. If --cortex is specified, then
    ?h.cortex.label will be used (this is created by automatically be
    recon-all). Even if you do not have a label of interest, it is
    recommended that you only smooth within cortex in order to prevent
    values from the medial wall from being smoothed into the
    surrounding cortical areas. At some point, this will be the
    default at which point you will have to use --no-cortex to turn it
    off.  This documentation will reflect the change. For --label-src
    and --label-trg, if you do not give it a full path, it will look
    in the subjects label dir. There is no need to specify both source
    and target unless you are smoothing on both (which you probably
    should not be doing).

  --frame framenumber

    When using paint/w output format, this specifies which frame to output. This
    format can store only one frame. The frame number is zero-based (default is 0).

  --mul Mul
  --div Div
    Multiply or divide the input by the given value

  --reshape

    Force mri_surf2surf to save the output as multiple 'slices'; has
    no effect for paint/w output format. For ico, the output will
    appear to be a 'volume' with Nv/R colums, 1 row, R slices and Nf
    frames, where Nv is the number of vertices on the surface. For
    icosahedrons, R=6. For others, R will be the prime factor of Nv
    closest to 6 (can be changed with --reshape-factor). Reshaping is
    for logistical purposes (eg, in the analyze/nifti format the size
    of a dimension cannot exceed 2^15). Use this flag to prevent this
    behavior. This has no effect when the output type is paint. At one
    point, it was the default to reshape.

  --reshape-factor Nfactor

    Attempt to reshape to Nfactor 'slices' (will choose closest prime
    factor) Default is 6.

  --reshape3d

    Reshape fsaverage (ico7) into 42 x 47 x 83

  --sd SUBJECTS_DIR

    Set SUBJECTS_DIR on the command line.

EXAMPLES:

1. Resample a subject's thickness of the left cortical hemisphere on to a
   7th order icosahedron and save in analyze4d format:

   mri_surf2surf --hemi lh --srcsubject bert
      --srcsurfval thickness --src_type curv
      --trgsubject ico --trgicoorder 7
      --trgsurfval bert-thickness-lh.img --trg_type analyze4d

2. Resample data on the icosahedron to the right hemisphere of subject bert.
   Note that both the source and target data are stored in mgh format
   as 'volume-encoded suface' data.

   mri_surf2surf --hemi rh --srcsubject ico --srcsurfval icodata-rh.mgh
      --trgsubject bert --trgsurfval ./bert-ico-rh.mgh

3. Convert the surface coordinates of the lh.white of a subject to a
   (talairach) average (ie, a subject created by make_average_subject):

   mri_surf2surf --s yoursubject --hemi lh --sval-tal-xyz white
     --trgsubject youraveragesubject --tval lh.white.yoursubject 
     --tval-xyz $SUBJECTS_DIR/fsaverage/mri/orig.mgz

   This will create youraveragesubject/surf/lh.white.yoursubject

4. Convert the surface coordinates of the lh.white of a subject to 
   the subject's functional space

   mri_surf2surf --reg register.dat template.nii.gz --hemi lh 
      --sval-xyz white --tval-xyz template.nii.gz --tval ./lh.white.func 
      --s yoursubject

   This will create lh.white.func in the current directory. template.nii.gz
   is a volume in the functional space. register.dat is the registration
   file between anatomical (target) and functional (movable) spaces. 
   View result with:  freeview -v template.nii.gz -f lh.white.func

   When using an LTA instead of a register.dat, do not include a target volume

   mri_surf2surf --reg register.lta --hemi lh 
      --sval-xyz white --tval-xyz template.nii.gz --tval ./lh.white.func 
      --s yoursubject

5. Extract surface normals of the white surface and save in a
   volume-encoded file:

   mri_surf2surf --s yoursubject --hemi lh --sval-nxyz white
      --tval lh.white.norm.mgh

   This will create youraveragesubject/surf/lh.white.yoursubject


6. Convert the annotation for one subject to the surface of another

  mri_surf2surf --srcsubject subj1 --trgsubject subj2 --hemi lh \\
    --sval-annot $SUBJECTS_DIR/subj1/label/lh.aparc.annot \\
    --tval       $SUBJECTS_DIR/subj2/label/lh.subj1.aparc.annot

   This will create $SUBJECTS_DIR/subj2/label/lh.subj1.aparc.annot.
   The --sval-annot flag will also change the map method to nnf so that
   the annot indices are not averaged. Note: this is not a substitute
   for running the cortical parcellation! The parcellations that it
   maps to the new subject may not be appropriate for that subject.

BUG REPORTS: send bugs to analysis-bugs@nmr.mgh.harvard.edu. Make sure
    to include the version and full command-line and enough information to
    be able to recreate the problem. Not that anyone does.

BUGS:

  When the output format is paint, the output file must be specified with
  a partial path (eg, ./data-lh.w) or else the output will be written into
  the subject's anatomical directory.

ENDHELP
*/

#undef X
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#include "mri.h"
#include "icosahedron.h"
#include "fio.h"
#include "pdf.h"

#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "mri2.h"
#include "mri_identify.h"

#include "bfileio.h"
#include "registerio.h"
//  extern char *ResampleVtxMapFile;
#include "resample.h"
#include "selxavgio.h"
#include "prime.h"
#include "version.h"
#include "colortab.h"
#include "fsenv.h"
#include "utils.h"
#include "cmdargs.h"
#include "proto.h"
#include "mri_circulars.h"

int DumpSurface(MRIS *surf, char *outfile);
MRI *MRIShksmooth(MRIS *Surf,
                  MRI *Src,
                  double sigma,
                  int nSmoothSteps,
                  MRI *Targ);
MRI *MRISheatkernel(MRIS *surf, double sigma);

double MRISareaTriangle(double x0, double y0, double z0,
                        double x1, double y1, double z1,
                        double x2, double y2, double z2);
int MRIStriangleAngles(double x0, double y0, double z0,
                       double x1, double y1, double z1,
                       double x2, double y2, double z2,
                       double *a0, double *a1, double *a2);
MRI *MRISdiffusionWeights(MRIS *surf);
MRI *MRISdiffusionSmooth(MRIS *Surf, MRI *Src, double GStd, MRI *Targ);
int MRISareNeighbors(MRIS *surf, int vtxno1, int vtxno2);

double MRISdiffusionEdgeWeight(MRIS *surf, int vtxno0, int vtxnonbr);
double MRISsumVertexFaceArea(MRIS *surf, int vtxno);
int MRIScommonNeighbors(MRIS *surf, int vtxno1, int vtxno2,
                        int *cvtxno1, int *cvtxno2);
int MRISdumpVertexNeighborhood(MRIS *surf, int vtxno);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
int GetNVtxsFromWFile(const char *wfile);
int GetICOOrderFromValFile(const char *filename, const char *fmt);
int GetNVtxsFromValFile(const char *filename, const char *fmt);
int dump_surf(char *fname, MRIS *surf, MRI *mri);
MATRIX *MRIleftRightRevMatrix(MRI *mri);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

const char *srcsurfregfile = NULL;
char *srchemi    = NULL;
const char *trgsurfregfile = NULL;
char *trghemi    = NULL;

char *srcsubject = NULL;
char *srcvalfile = NULL;
const char *srctypestring = "";
int   srctype = MRI_VOLUME_TYPE_UNKNOWN;
MRI  *SrcVals, *SrcHits, *SrcDist;
MRI_SURFACE *SrcSurfReg;
char *SrcHitFile = NULL;
char *SrcDistFile = NULL;
int nSrcVtxs = 0;
int SrcIcoOrder = -1;

int UseSurfSrc=0; // Get source values from surface, eg, xyz
const char *SurfSrcName=NULL;
MRI_SURFACE *SurfSrc, *SurfTrg;
MATRIX *XFM=NULL, *XFMSubtract=NULL;
#define SURF_SRC_XYZ     1
#define SURF_SRC_TAL_XYZ 2
#define SURF_SRC_AREA    3
#define SURF_SRC_NXYZ    4 // surface normals
#define SURF_SRC_RIP     5 // rip flag
#define SURF_SRC_ANNOT   6 // surface annotation
int UseSurfTarg=0; // Put Src XYZ into a target surface
int ApplyReg=0;
char *AnnotFile = NULL;

char *trgsubject = NULL;
char *trgvalfile = NULL;
const char *trgtypestring = "";
int   trgtype = MRI_VOLUME_TYPE_UNKNOWN;
MRI  *TrgVals, *TrgValsSmth, *TrgHits, *TrgDist;
MRI_SURFACE *TrgSurfReg;
char *TrgHitFile = NULL;
char *TrgDistFile = NULL;
int TrgIcoOrder;

MRI  *mritmp;
int  reshape = 0;
int  reshapefactor;
int  reshape3d=0;

const char *mapmethod = "nnfr";

int UseHash = 1;
int framesave = 0;
float IcoRadius = 100.0;
int nthstep, nnbrs, nthnbr, nbrvtx, frame;
int nSmoothSteps = 0;
double fwhm=0, gstd;

double fwhm_Input=0, gstd_Input=0;
int nSmoothSteps_Input = 0;
int usediff = 0;
char *LabelFile=NULL, *LabelFile_Input=NULL; // for masking smoothing

int debug = 0;
char *SrcDumpFile  = NULL;
char *TrgDumpFile = NULL;

char *SUBJECTS_DIR = NULL;
char *FREESURFER_HOME = NULL;
SXADAT *sxa;
FILE *fp;

char tmpstr[2000];

int ReverseMapFlag = 0;
int cavtx = 0; /* command-line vertex -- for debugging */
int jac = 0;

MRI *sphdist;

int SynthPDF = 0;
int SynthOnes = 0;
int SynthSeed = -1;
double ProjDepth;
int    DoProj = 0;
int    ProjType;
int reshapefactortarget = 6;

int DoNormVar=0;
int NormVar(MRI *mri, MRI *mask);
int UseCortexLabel = 0;

int OKToRevFaceOrder = 1;
int RevFaceOrder = 0;

char *RMSDatFile = NULL;
char *RMSMaskFile = NULL;
MRI *RMSMask = NULL;
int SplitFrames=0;
int ConvGaussian = 0;

int UseDualHemi = 0; // Assume ?h.?h.surfreg file name, source only
MRI *RegTarg = NULL;
int UseOldSurf2Surf = 1;
char *PatchFile=NULL, *SurfTargName=NULL;
int nPatchDil=0;
struct utsname uts;
char *cmdline, cwd[2000];
char *TrgSurfVolFile=NULL;
MRI  *TrgSurfVol=NULL;
int   prunemask = 0;
float prune_thr = FLT_MIN;
int DoMultiply = 0;
double MultiplyVal = 0;

/*---------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
  int f,tvtx,svtx,n,err;
  float *framepower = NULL,x,y,z;
  char fname[4000];
  int nTrg121,nSrc121,nSrcLost;
  int nTrgMulti,nSrcMulti;
  float MnTrgMultiHits,MnSrcMultiHits, val;
  int nargs;
  FACE *face;
  VERTEX *vtx0,*vtx1,*vtx2;
  double area, a0, a1, a2, d, dmin, dmax, dsum;
  COLOR_TABLE *ctab=NULL;
  LABEL *MaskLabel;
  MRI *inmask = NULL, *outmask=NULL, *mri2=NULL;
  char *stem, *ext;

  nargs = handleVersionOption(argc, argv, "mri_surf2surf");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  // Make sure that MRImask() does not try to compensate for the geometry
  setenv("FS_MRIMASK_ALLOW_DIFF_GEOM","0",1);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  //SrcSurfReg = ReadIcoByOrder(7,100);
  //MRISwrite(SrcSurfReg,"lh.ico7");
  //exit(1);

  if (argc == 0) usage_exit();

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL) {
    fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  FREESURFER_HOME = getenv("FREESURFER_HOME") ;
  if (FREESURFER_HOME==NULL) {
    fprintf(stderr,"ERROR: FREESURFER_HOME not defined in environment\n");
    exit(1);
  }

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if(TrgSurfVolFile){
    TrgSurfVol = MRIreadHeader(TrgSurfVolFile,MRI_VOLUME_TYPE_UNKNOWN);
    if(TrgSurfVol == NULL) exit(1);
  }

  /* --------- Load the registration surface for source subject --------- */
  if (!strcmp(srcsubject,"ico")) { /* source is ico */
    if (SrcIcoOrder == -1) {
      SrcIcoOrder = GetICOOrderFromValFile(srcvalfile,srctypestring);
    }
    sprintf(fname,"%s/lib/bem/ic%d.tri",FREESURFER_HOME,SrcIcoOrder);
    SrcSurfReg = ReadIcoByOrder(SrcIcoOrder, IcoRadius);
    printf("Source Ico Order = %d\n",SrcIcoOrder);
  } 
  else {
    // Set source reg depending on whether hemis are same or diff
    // Changed to this on 11/30/97
    if(!strcmp(srchemi,trghemi) && UseDualHemi == 0) { 
      // hemis are the same
      sprintf(fname,"%s/%s/surf/%s.%s",
              SUBJECTS_DIR,srcsubject,srchemi,srcsurfregfile);
    } 
    else {
      // hemis are the different
      if(UseDualHemi)
	sprintf(fname,"%s/%s/surf/%s.%s.%s",SUBJECTS_DIR,
		srcsubject,srchemi,trghemi,srcsurfregfile);
      else
	sprintf(fname,"%s/%s/surf/%s.%s",
		SUBJECTS_DIR,srcsubject,srchemi,srcsurfregfile);
    }
    printf("Reading source surface reg %s\n",fname);
    SrcSurfReg = MRISread(fname) ;
    if (cavtx > 0)
      printf("cavtx = %d, srcsurfregfile: %g, %g, %g\n",cavtx,
             SrcSurfReg->vertices[cavtx].x,
             SrcSurfReg->vertices[cavtx].y,
             SrcSurfReg->vertices[cavtx].z);
  }
  if (!SrcSurfReg) {
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;
  }

  MRIScomputeMetricProperties(SrcSurfReg);
  for (n = 0; n < SrcSurfReg->nfaces && 0; n++) {
    face = &SrcSurfReg->faces[n];
    vtx0 = &SrcSurfReg->vertices[face->v[0]];
    vtx1 = &SrcSurfReg->vertices[face->v[1]];
    vtx2 = &SrcSurfReg->vertices[face->v[2]];
    area = MRISareaTriangle(vtx0->x,vtx0->y,vtx0->z,
                            vtx1->x,vtx1->y,vtx1->z,
                            vtx2->x,vtx2->y,vtx2->z);
    MRIStriangleAngles(vtx0->x,vtx0->y,vtx0->z,
                       vtx1->x,vtx1->y,vtx1->z,
                       vtx2->x,vtx2->y,vtx2->z,&a0,&a1,&a2);
    printf("n=%d, area = %f, %f %f %f   %f\n",n,face->area,a0,a1,a2,a0+a1+a2);
  }

  /* ------------------ load the source data ----------------------------*/
  printf("Loading source data\n");
  if (!strcmp(srctypestring,"curv")) { /* curvature file */
    if (fio_FileExistsReadable(srcvalfile)) {
      memset(fname,0,strlen(srcvalfile)+1);
      memmove(fname,srcvalfile,strlen(srcvalfile));
    } 
    else sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,srchemi,srcvalfile);
    printf("Reading curvature file %s\n",fname);
    if (MRISreadCurvatureFile(SrcSurfReg, fname) != 0) {
      printf("ERROR: reading curvature file\n");
      exit(1);
    }
    SrcVals = MRIcopyMRIS(NULL, SrcSurfReg, 0, "curv");
  } else if (!strcmp(srctypestring,"paint") || !strcmp(srctypestring,"w")) {
    MRISreadValues(SrcSurfReg,srcvalfile);
    SrcVals = MRIcopyMRIS(NULL, SrcSurfReg, 0, "val");
  } 
  else if (UseSurfSrc) {
    if(PatchFile) SurfSrcName = srcsurfregfile;
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,srchemi,SurfSrcName);
    printf("Reading surface file %s\n",fname);
    SurfSrc = MRISread(fname);
    if(SurfSrc==NULL)  exit(1);

    if(PatchFile){
      if(fio_FileExistsReadable(PatchFile)) strcpy(fname,PatchFile);
      else
	sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,srchemi,PatchFile);
      printf("Reading source patch file %s\n",fname);
      err = MRISreadPatch(SurfSrc,fname);
      if(err) exit(1);
    }

    MRIScomputeMetricProperties(SurfSrc);
    
    MRISfreeDistsButNotOrig(SurfSrc);
      // MRISsetXYZ will invalidate all of these,
      // so make sure they are recomputed before being used again!

    if(UseSurfSrc == SURF_SRC_XYZ || UseSurfSrc == SURF_SRC_TAL_XYZ) {
      if(DoProj) {
        if(ProjType == 2) {
          // Load thickness
          sprintf(tmpstr,"%s/%s/surf/%s.thickness",
                  SUBJECTS_DIR,srcsubject,srchemi);
          printf("Loading thickness %s\n",tmpstr);
          err = MRISreadCurvatureFile(SurfSrc,tmpstr);
          if(err) {
            exit(1);
          }
        }
        // Project
        printf("Projecting surface %g along the normal (%d)\n",
               ProjDepth,ProjType);
        for(n=0; n < SurfSrc->nvertices; n++) {
          if(ProjType == 1) {
            ProjNormDist(&x, &y, &z,SurfSrc, n, ProjDepth);
          }
          if(ProjType == 2) {
            ProjNormFracThick(&x,&y,&z,SurfSrc, n, ProjDepth);
          }
          //printf("%5d (%g,%g,%g) (%g,%g,%g) \n",n,SurfSrc->vertices[n].x,
          //SurfSrc->vertices[n].y,SurfSrc->vertices[n].z,x,y,z);
          MRISsetXYZ(SurfSrc,n, x,y,z);
        }
      }
      if(UseSurfSrc == SURF_SRC_TAL_XYZ) {
        XFM = DevolveXFM(srcsubject, NULL, NULL);
        if (XFM == NULL) {
          exit(1);
        }
        printf("Applying MNI305 talairach transform\n");
        MatrixPrint(stdout,XFM);
        MRISmatrixMultiply(SurfSrc,XFM);
      }
      if(ApplyReg) {
        if(XFMSubtract) {
          // This is primarily for testing (computing rms diffs below)
          printf("Computing diff in registration matrices\n");
          XFM = MatrixSubtract(XFM,XFMSubtract,XFM);
        }
        printf("Applying linear registration transform\n");
        MatrixPrint(stdout,XFM);
        MRISmatrixMultiply(SurfSrc,XFM);
        if(MatrixDeterminant(XFM) < 0.0 && OKToRevFaceOrder) {
          printf("Determinant of linear transform is negative, "
                 "so reversing face order\n");
          RevFaceOrder = 1;
        }
        if(RMSDatFile) {
          // This is primarily for testing
          printf("Computing RMS\n");
          if(RMSMaskFile) {
            printf("Loading RMS Mask %s\n",RMSMaskFile);
            RMSMask = MRIread(RMSMaskFile);
            if(RMSMask == NULL) {
              exit(1);
            }
            if(RMSMask->width != SurfSrc->nvertices) {
              printf("Dimension Mismatch: %d %d\n",
                     RMSMask->width,SurfSrc->nvertices);
              exit(1);
            }
          }
          dmin = 10e10;
          dmax = 0.0;
          dsum = 0.0;
          f = 0;
          for(n=0; n < SurfSrc->nvertices; n++) {
            if(RMSMask && fabs(MRIgetVoxVal(RMSMask,n,0,0,0)) < 1e-6) {
              continue;
            }
            f++;
            x = SurfSrc->vertices[n].x;
            y = SurfSrc->vertices[n].y;
            z = SurfSrc->vertices[n].z;
            d = sqrt(x*x + y*y + z*z);
            dsum += d;
            if(d > dmax) {
              dmax = d;
            }
            if(d < dmin) {
              dmin = d;
            }
          }
          d = dsum/f;
          printf("RMS %13.10lf %13.10lf %13.10lf %d\n",d,dmin,dmax,f);
          fp = fopen(RMSDatFile,"w");
          fprintf(fp,"%13.10lf %13.10lf %13.10lf %d\n",d,dmin,dmax,f);
          fclose(fp);
          exit(0);
        }
      }
      SrcVals = MRIcopyMRIS(NULL, SurfSrc, 2, "z"); // start at z to autoalloc
      MRIcopyMRIS(SrcVals, SurfSrc, 0, "x");
      MRIcopyMRIS(SrcVals, SurfSrc, 1, "y");
    }
    if (UseSurfSrc == SURF_SRC_NXYZ) {
      printf("Extracting surface normals\n");
      SrcVals = MRIcopyMRIS(NULL, SurfSrc, 2, "nz"); // start at z to autoalloc
      MRIcopyMRIS(SrcVals, SurfSrc, 0, "nx");
      MRIcopyMRIS(SrcVals, SurfSrc, 1, "ny");
    }
    if(UseSurfSrc == SURF_SRC_AREA) {
      // Not affected by loading ?h.white.avg.area.mgh because
      // this uses area and not group_avg_area.
      SrcVals = MRIcopyMRIS(NULL, SurfSrc, 0, "area");
      if (SurfSrc->group_avg_surface_area > 0) {
        val = SurfSrc->group_avg_surface_area / SurfSrc->total_area;
        MRIscalarMul(SrcVals,SrcVals,val);
        // Always fix now (4/9/10)
        //if (getenv("FIX_VERTEX_AREA") != NULL) {
        //printf("INFO: Fixing group surface area\n");
        //val = SurfSrc->group_avg_surface_area / SurfSrc->total_area;
        //MRIscalarMul(SrcVals,SrcVals,val);
        //}
      }
    }
    if (UseSurfSrc == SURF_SRC_ANNOT) {
      err = MRISreadAnnotation(SurfSrc, AnnotFile);
      if (err) {
        exit(1);
      }
      SrcVals = MRISannotIndex2Seg(SurfSrc);
      ctab = CTABdeepCopy(SurfSrc->ct);
    }
    if(UseSurfSrc == SURF_SRC_RIP) SrcVals = MRIcopyMRIS(NULL, SurfSrc, 0, "ripflag");
    MRISfree(&SurfSrc);
  } 
  else { /* Use MRIreadType */
    SrcVals =  MRIreadType(srcvalfile,srctype);
    if (SrcVals == NULL) {
      printf("ERROR: could not read %s as type %d\n",srcvalfile,srctype);
      exit(1);
    }
    if (SrcVals->height != 1 || SrcVals->depth != 1) {
      reshapefactor = SrcVals->height * SrcVals->depth;
      printf("Reshaping %d\n",reshapefactor);
      mritmp = mri_reshape(SrcVals, reshapefactor*SrcVals->width,
                           1, 1, SrcVals->nframes);
      MRIfree(&SrcVals);
      SrcVals = mritmp;
      reshapefactor = 0; /* reset for output */
    }

    if (SrcVals->width != SrcSurfReg->nvertices) {
      fprintf(stderr,"ERROR: dimension inconsistency in source data\n");
      fprintf(stderr,"       Number of surface vertices = %d\n",
              SrcSurfReg->nvertices);
      fprintf(stderr,"       Number of value vertices = %d\n",SrcVals->width);
      exit(1);
    }
    if(SrcVals->type != MRI_FLOAT) {
      printf("Converting source to float\n");
      mritmp = MRISeqchangeType(SrcVals,MRI_FLOAT,0,0,0);
      if (mritmp == NULL) {
        printf("ERROR: could change type\n");
        exit(1);
      }
      MRIfree(&SrcVals);
      SrcVals = mritmp;
    }

    if (is_sxa_volume(srcvalfile)) {
      printf("INFO: Source volume detected as selxavg format\n");
      sxa = ld_sxadat_from_stem(srcvalfile);
      if (sxa == NULL) {
        exit(1);
      }
      framepower = sxa_framepower(sxa,&f);
      if (f != SrcVals->nframes) {
        fprintf(stderr," number of frames is incorrect (%d,%d)\n",
                f,SrcVals->nframes);
        exit(1);
      }
      printf("INFO: Adjusting Frame Power\n");
      fflush(stdout);
      mri_framepower(SrcVals,framepower);
    }
  }
  if (SrcVals == NULL) {
    fprintf(stderr,"ERROR loading source values from %s\n",srcvalfile);
    exit(1);
  }
  if(DoMultiply){
    printf("Multiplying by %lf\n",MultiplyVal);
    MRImultiplyConst(SrcVals, MultiplyVal, SrcVals);
  }

  n = SrcVals->width * SrcVals->height * SrcVals->depth;
  if(SrcSurfReg->nvertices != n) {
    printf("ERROR: dimension mismatch between surface reg (%d)\n",
           SrcSurfReg->nvertices);
    printf("and source data (%d)\n",n);
    exit(1);
  }
  if (SynthPDF != 0) {
    if (SynthSeed < 0) {
      SynthSeed = PDFtodSeed();
    }
    printf("INFO: synthesizing, pdf = %d, seed = %d\n",SynthPDF,SynthSeed);
    srand48(SynthSeed);
    MRIrandn(SrcVals->width, SrcVals->height, SrcVals->depth,
             SrcVals->nframes,0, 1, SrcVals);
  }
  if (SynthOnes != 0) {
    printf("INFO: filling input with all 1s\n");
    MRIconst(SrcVals->width, SrcVals->height, SrcVals->depth,
             SrcVals->nframes, 1, SrcVals);
  }

  if (SrcDumpFile != NULL) {
    printf("Dumping input to %s\n",SrcDumpFile);
    DumpSurface(SrcSurfReg, SrcDumpFile);
  }

  /* Smooth input if desired */
  if (fwhm_Input > 0) {
    nSmoothSteps_Input =
      MRISfwhm2nitersSubj(fwhm_Input,srcsubject,srchemi,"white");
    if (nSmoothSteps_Input == -1) {
      exit(1);
    }
    printf("Approximating gaussian smoothing of source with fwhm = %lf,\n"
           "std = %lf, with %d iterations of nearest-neighbor smoothing\n",
           fwhm_Input,gstd_Input,nSmoothSteps_Input);
  }

  if(LabelFile_Input != NULL) {
    printf("Reading source subject label mask %s\n",LabelFile_Input);
    MaskLabel = LabelRead(srcsubject, LabelFile_Input);
    if(MaskLabel == NULL) exit(1);
    inmask = MRISlabel2Mask(SrcSurfReg, MaskLabel, NULL);
  } 
  else inmask = NULL;

  // This removes voxels from the mask that do not have non-zero in every frame
  // To prevent smoothing of those voxels into others
  if(prunemask && nSmoothSteps_Input > 0) {
    printf("Pruning input smoothing mask by thr: %e\n", prune_thr);
    inmask = MRIframeBinarize(SrcVals,FLT_MIN,inmask);
  }

  if(nSmoothSteps_Input > 0) {
    if(! ConvGaussian) {
      printf("NN smoothing input with n = %d\n",nSmoothSteps_Input);
      MRISsmoothMRI(SrcSurfReg, SrcVals, nSmoothSteps_Input, inmask, SrcVals);
    } 
    else {
      printf("Convolving with gaussian\n");
      MRISgaussianSmooth(SrcSurfReg, SrcVals, gstd_Input, SrcVals, 3.5);
      //MRIShksmooth(SrcSurfReg, SrcVals, gstd_Input, nSmoothSteps_Input, SrcVals);
    }
  }
  else {
    if(inmask){
      printf("masking the input\n");
      mritmp = MRImask(SrcVals,inmask,SrcVals,0,0);
    }
  }

  if(strcmp(srcsubject,trgsubject) || strcmp(srchemi,trghemi) ||
     strcmp(srcsurfregfile,trgsurfregfile) ||
     (!strcmp(srcsubject,"ico") && !strcmp(trgsubject,"ico") && SrcIcoOrder != TrgIcoOrder)){
    /* ------- Source and Target Subjects or Hemis are different ------ */
    /* ------- Load the registration surface for target subject ------- */
    if (!strcmp(trgsubject,"ico")) {
      sprintf(fname,"%s/lib/bem/ic%d.tri",FREESURFER_HOME,TrgIcoOrder);
      TrgSurfReg = ReadIcoByOrder(TrgIcoOrder, IcoRadius);
      reshapefactor = 6;
    } else {
      // Use same target regardless of whether hemis are the same or diff
      // Changed to this on 11/30/97
      sprintf(fname,"%s/%s/surf/%s.%s",
              SUBJECTS_DIR,trgsubject,trghemi,trgsurfregfile);
      printf("Reading target surface reg %s\n",fname);
      TrgSurfReg = MRISread(fname) ;
    }
    if (!TrgSurfReg)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface %s",
                Progname, fname) ;
    printf("Done\n");

    if (!strcmp(mapmethod,"nnfr")) {
      ReverseMapFlag = 1;
    } else {
      ReverseMapFlag = 0;
    }

    /*-------------------------------------------------------------*/
    /* Map the values from the surface to surface */
    if(UseOldSurf2Surf){
      printf("Using surf2surf_nnfr()\n");
      if (!jac) {
	printf("Mapping Source Volume onto Source Subject Surface\n");
	TrgVals = surf2surf_nnfr(SrcVals, SrcSurfReg,TrgSurfReg,
				 &SrcHits,&SrcDist,&TrgHits,&TrgDist,
				 ReverseMapFlag,UseHash);
      } else {
	printf("Mapping Source Volume onto Source Subject Surface "
	       "with Jacobian Correction\n");
	TrgVals = surf2surf_nnfr_jac(SrcVals, SrcSurfReg,TrgSurfReg,
				     &SrcHits,&SrcDist,&TrgHits,&TrgDist,
				     ReverseMapFlag,UseHash);
      }
      /* Compute some stats on the mapping number of srcvtx mapping to a
	 target vtx*/
      nTrg121 = 0;
      MnTrgMultiHits = 0.0;
      for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
	n = MRIFseq_vox(TrgHits,tvtx,0,0,0);
	if (n == 1) {
	  nTrg121++;
	} else {
	  MnTrgMultiHits += n;
	}
      }
      nTrgMulti = TrgSurfReg->nvertices - nTrg121;
      if (nTrgMulti > 0) {
	MnTrgMultiHits = (MnTrgMultiHits/nTrgMulti);
      } else {
	MnTrgMultiHits = 0;
      }
      printf("nTrg121 = %5d, nTrgMulti = %5d, MnTrgMultiHits = %g\n",
	     nTrg121,nTrgMulti,MnTrgMultiHits);
      
      /* Compute some stats on the mapping number of trgvtxs mapped from a
	 source vtx*/
      nSrc121 = 0;
      nSrcLost = 0;
      MnSrcMultiHits = 0.0;
      for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
	n = MRIFseq_vox(SrcHits,svtx,0,0,0);
	if (n == 1) {
	  nSrc121++;
	} else if (n == 0) {
	  nSrcLost++;
	} else {
	  MnSrcMultiHits += n;
	}
      }
      nSrcMulti = SrcSurfReg->nvertices - nSrc121;
      if (nSrcMulti > 0) {
	MnSrcMultiHits = (MnSrcMultiHits/nSrcMulti);
      } else {
	MnSrcMultiHits = 0;
      }
      
      printf("nSrc121 = %5d, nSrcLost = %5d, nSrcMulti = %5d, "
	     "MnSrcMultiHits = %g\n", nSrc121,nSrcLost,nSrcMulti,
	     MnSrcMultiHits);
      
      /* save the Source Hits into a .w file */
      if (SrcHitFile != NULL) {
	printf("INFO: saving source hits to %s\n",SrcHitFile);
	MRIScopyMRI(SrcSurfReg, SrcHits, 0, "val");
	//for(vtx = 0; vtx < SrcSurfReg->nvertices; vtx++)
	//SrcSurfReg->vertices[vtx].val = MRIFseq_vox(SrcHits,vtx,0,0,0) ;
	MRISwriteValues(SrcSurfReg, SrcHitFile) ;
      }
      /* save the Source Distance into a .w file */
      if (SrcDistFile != NULL) {
	printf("INFO: saving source distance to %s\n",SrcDistFile);
	MRIScopyMRI(SrcSurfReg, SrcDist, 0, "val");
	MRISwriteValues(SrcSurfReg, SrcDistFile) ;
	//for(vtx = 0; vtx < SrcSurfReg->nvertices; vtx++)
	//SrcSurfReg->vertices[vtx].val = MRIFseq_vox(SrcDist,vtx,0,0,0) ;
      }
      /* save the Target Hits into a .w file */
      if (TrgHitFile != NULL) {
	printf("INFO: saving target hits to %s\n",TrgHitFile);
	MRIScopyMRI(TrgSurfReg, TrgHits, 0, "val");
	MRISwriteValues(TrgSurfReg, TrgHitFile) ;
	//for(vtx = 0; vtx < TrgSurfReg->nvertices; vtx++)
	//TrgSurfReg->vertices[vtx].val = MRIFseq_vox(TrgHits,vtx,0,0,0) ;
      }
      /* save the Target Hits into a .w file */
      if (TrgDistFile != NULL) {
	printf("INFO: saving target distance to %s\n",TrgDistFile);
	err = MRIwrite(TrgDist,TrgDistFile);
	if(err) {
	  printf("ERROR: writing %s\n",TrgDistFile);
	  exit(1);
	}
      }
    }
    else {
      printf("Using MRISapplyReg()\n");
      MRIS *SurfRegList[2];
      SurfRegList[0] = SrcSurfReg;
      SurfRegList[1] = TrgSurfReg;
      TrgVals = MRISapplyReg(SrcVals, SurfRegList, 2, ReverseMapFlag,jac,UseHash);
    }

  } else {
    /* --- Source and Target Subjects are the same --- */
    printf("INFO: trgsubject = srcsubject\n");
    TrgSurfReg = SrcSurfReg;
    TrgVals = SrcVals;
  }

  /* Smooth output if desired */
  if (fwhm > 0) {
    nSmoothSteps = MRISfwhm2nitersSubj(fwhm,trgsubject,trghemi,"white");
    if (nSmoothSteps == -1) {
      exit(1);
    }
    printf("Approximating gaussian smoothing of target with fwhm = %lf,\n"
           " std = %lf, with %d iterations of nearest-neighbor smoothing\n",
           fwhm,gstd,nSmoothSteps);
  }

  if(LabelFile != NULL) {
    printf("Reading target space mask label %s\n",LabelFile);
    MaskLabel = LabelRead(trgsubject, LabelFile);
    if(MaskLabel == NULL) exit(1);
    outmask = MRISlabel2Mask(TrgSurfReg, MaskLabel, NULL);
  } 
  else  outmask = NULL;

  // This removes voxels from the mask that do not have non-zero in every frame
  // To prevent smoothing of those voxels into others
  if(prunemask && nSmoothSteps > 0) {
    printf("Pruning output smoothing mask by thr: %e\n", prune_thr);
    outmask = MRIframeBinarize(SrcVals,FLT_MIN,outmask);
  }

  if(nSmoothSteps > 0) {
    if(! ConvGaussian) {
      printf("NN smoothing output with n = %d\n",nSmoothSteps);
      MRISsmoothMRI(TrgSurfReg, TrgVals, nSmoothSteps, outmask, TrgVals);
    } else {
      printf("Convolving with gaussian (assuming sphere)\n");
      MRISgaussianSmooth(TrgSurfReg, TrgVals, gstd, TrgVals, 3.5);
      //printf("Diffusion Smoothing\n");
      //MRISdiffusionSmooth(TrgSurfReg, TrgVals, gstd, TrgVals);
      //printf("HK Smoothing\n");
      //MRIShksmooth(SrcSurfReg, SrcVals, gstd, nSmoothSteps, SrcVals);
    }
  }
  else {
    if(outmask){
      printf("masking the input\n");
      mritmp = MRImask(TrgVals,outmask,TrgVals,0,0);
    }
  }

  /* readjust frame power if necessary */
  if (is_sxa_volume(srcvalfile)) {
    printf("INFO: Readjusting Frame Power\n");
    fflush(stdout);
    for (f=0; f < TrgVals->nframes; f++) {
      framepower[f] = 1.0/framepower[f];
    }
    mri_framepower(TrgVals,framepower);
    sxa->nrows = 1;
    sxa->ncols = TrgVals->width;
  }

  if (TrgDumpFile != NULL) {
    /* Dump before reshaping */
    printf("Dumping output to %s\n",TrgDumpFile);
    dump_surf(TrgDumpFile,TrgSurfReg,TrgVals);
  }

  /* ------------ save the target data -----------------------------*/
  printf("Saving target data\n");
  if (!strcmp(trgtypestring,"paint") || !strcmp(trgtypestring,"w")) {
    printf("Saving as paint format\n");
    MRIScopyMRI(TrgSurfReg, TrgVals, framesave, "val");
    MRISwriteValues(TrgSurfReg,trgvalfile);
  } 
  else if (!strcmp(trgtypestring,"curv")) {
    MRIScopyMRI(TrgSurfReg, TrgVals, framesave, "curv");
    MRISwriteCurvature(TrgSurfReg,trgvalfile);
  } 
  else if (UseSurfTarg) {
    MRIScopyMRI(TrgSurfReg,TrgVals,0,"x");
    MRIScopyMRI(TrgSurfReg,TrgVals,1,"y");
    MRIScopyMRI(TrgSurfReg,TrgVals,2,"z");
    if(RevFaceOrder) {
      printf("Reversing Face Order\n");
      MRISreverseFaceOrder(TrgSurfReg);
    }
    if(RegTarg) getVolGeom(RegTarg, &TrgSurfReg->vg);
    else if(TrgSurfVol) getVolGeom(TrgSurfVol, &TrgSurfReg->vg);
    MRISwrite(TrgSurfReg, trgvalfile);
  } 
  else if (UseSurfSrc == SURF_SRC_ANNOT) {
    printf("Converting to target annot\n");
    err = MRISseg2annot(TrgSurfReg,TrgVals,ctab);
    if (err) {
      exit(1);
    }
    printf("Saving to target annot %s\n",trgvalfile);
    MRISwriteAnnotation(TrgSurfReg, trgvalfile);
  } 
  else if(PatchFile){
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,trghemi,SurfTargName);
    printf("Reading surface file for output patch %s\n",fname);
    SurfTrg = MRISread(fname);
    if(SurfTrg==NULL)  exit(1);
    printf("Ripping patch\n");
    for(tvtx = 0; tvtx < SurfTrg->nvertices; tvtx++){
      if(MRIgetVoxVal(TrgVals,tvtx,0,0,0) < 0.000001) continue;
      SurfTrg->vertices[tvtx].ripflag = 1;
    }
    printf("Dilating %d\n",nPatchDil);
    MRISdilateRipped(SurfTrg, nPatchDil);
    MRISsetRipInFacesWithRippedVertices(SurfTrg);
    SurfTrg->patch = 1 ;
    SurfTrg->status = MRIS_CUT ;
    for (tvtx = 0 ; tvtx < SurfTrg->nvertices ; tvtx++)
      if (SurfTrg->vertices_topology[tvtx].num == 0 || SurfTrg->vertices_topology[tvtx].vnum == 0)
	SurfTrg->vertices[tvtx].ripflag = 1 ;
    MRISupdateSurface(SurfTrg);
    printf("Writing patch to %s\n", trgvalfile);
    MRISwritePatch(SurfTrg, trgvalfile);
  }
  else {
    if (reshape) {
      if(reshapefactor == 0) {
        if(TrgSurfReg->nvertices == 163842) {
          reshapefactor = 6;
        } else {
          reshapefactor = GetClosestPrimeFactor(TrgVals->width,
                                                reshapefactortarget);
        }
      }
      printf("Reshaping %d (nvertices = %d)\n",reshapefactor,TrgVals->width);
      mritmp = mri_reshape(TrgVals, TrgVals->width / reshapefactor,
                           1, reshapefactor,TrgVals->nframes);
      if (mritmp == NULL) {
        printf("ERROR: mri_reshape could not alloc\n");
        return(1);
      }
      MRIfree(&TrgVals);
      TrgVals = mritmp;
    }
    if (reshape3d) {
      if(TrgSurfReg->nvertices != 163842) {
	printf("ERROR: subject must have 163842 vertices to 3d reshape\n");
	exit(1);
      }
      printf("Reshape 3d\n");
      mritmp = mri_reshape(TrgVals, 42,47,83,TrgVals->nframes);
      if (mritmp == NULL) {
        printf("ERROR: mri_reshape could not alloc\n");
        return(1);
      }
      MRIfree(&TrgVals);
      TrgVals = mritmp;
    }
    if(DoNormVar) {
      NormVar(TrgVals, NULL);
    }
    if(! SplitFrames) {
      printf("Saving to %s\n",trgvalfile);
      err = MRIwriteType(TrgVals,trgvalfile,trgtype);
      if(err) {
        printf("ERROR: writing %s\n",trgvalfile);
        exit(1);
      }
      if (is_sxa_volume(srcvalfile)) {
        sv_sxadat_by_stem(sxa,trgvalfile);
      }
    } else {
      stem = IDstemFromName(trgvalfile);
      ext = IDextensionFromName(trgvalfile);
      printf("Splitting frames, stem = %s, ext = %s\n",stem,ext);
      mri2 = NULL;
      for(f=0; f < TrgVals->nframes; f++) {
        mri2 = MRIcopyFrame(TrgVals, mri2, f, 0);
        sprintf(tmpstr,"%s%04d.%s",stem,f,ext);
        printf("%2d %s\n",f,tmpstr);
        err = MRIwrite(mri2, tmpstr);
        if (err != NO_ERROR) {
          printf("ERROR: failure writing %s\n",tmpstr);
          exit(1);
        }
      }
    }
  }
  return(0);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option, tmp[STRLEN] ;
  float ipr, bpr, intensity;
  int float2int,err;
  char *regsubject;
  MATRIX *M;
  MRI *lrmri;

  if (argc < 1) {
    usage_exit();
  }

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) {
      printf("%d %s\n",nargc,option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help")) {
      print_help() ;
    }

    else if (!strcasecmp(option, "--version")) {
      print_version() ;
    }

    else if (!strcasecmp(option, "--debug")) debug = 1;
    else if (!strcasecmp(option, "--prune"))    prunemask = 1;
    else if (!strcasecmp(option, "--no-prune")) prunemask = 0;
    else if (!strcasecmp(option, "--prune_thr")){
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&prune_thr); 
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--old"))UseOldSurf2Surf = 1;
    else if (!strcasecmp(option, "--new")) UseOldSurf2Surf = 0;
    else if (!strcasecmp(option, "--usehash")) {
      UseHash = 1;
    } else if (!strcasecmp(option, "--hash")) {
      UseHash = 1;
    } else if (!strcasecmp(option, "--dontusehash")) {
      UseHash = 0;
    } else if (!strcasecmp(option, "--nohash")) {
      UseHash = 0;
    } else if (!strcasecmp(option, "--noreshape")) {
      reshape = 0;
    } else if (!strcasecmp(option, "--reshape")) {
      reshape = 1;
    } else if (!strcasecmp(option, "--reshape3d")) {
      reshape = 0;
      reshape3d = 1;
    } else if (!strcasecmp(option, "--usediff")) {
      usediff = 1;
    } else if (!strcasecmp(option, "--nousediff")) {
      usediff = 0;
    } else if (!strcasecmp(option, "--synth")) {
      SynthPDF = 1;
    } else if (!strcasecmp(option, "--ones")) {
      SynthOnes = 1;
    } else if (!strcasecmp(option, "--jac")) {
      jac = 1;
    } else if (!strcasecmp(option, "--norm-var")) {
      DoNormVar = 1;
    } else if (!strcasecmp(option, "--split")) {
      SplitFrames = 1;
    } else if (!strcasecmp(option, "--conv")) {
      ConvGaussian = 1;
    } else if (!strcasecmp(option, "--no-rev-face-order")) {
      OKToRevFaceOrder = 0;
    }

    else if (!strcmp(option, "--seed")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&SynthSeed);
      nargsused = 1;
    } else if (!strcmp(option, "--sd")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      FSENVsetSUBJECTS_DIR(pargv[0]);
      SUBJECTS_DIR = pargv[0] ;
      nargsused = 1;
    } else if (!strcmp(option, "--s")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      srcsubject = pargv[0];
      trgsubject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srcsubject")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      srcsubject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srcsurfval") || !strcmp(option, "--sval")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      srcvalfile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--sval-xyz")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SurfSrcName = pargv[0];
      UseSurfSrc = SURF_SRC_XYZ;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--patch")) {
      if(nargc < 3) argnerr(option,3);
      PatchFile = pargv[0];
      SurfTargName = pargv[1];
      sscanf(pargv[2],"%d",&nPatchDil);
      UseSurfSrc = SURF_SRC_RIP;
      //UseSurfTarg = 1;
      mapmethod = "nnfr";
      nargsused = 3;
    } 
    else if (!strcmp(option, "--projabs")) {
      if(nargc < 2) {
        argnerr(option,2);
      }
      ProjType = 1;
      SurfSrcName = pargv[0];
      sscanf(pargv[1],"%lf",&ProjDepth);
      UseSurfSrc  = SURF_SRC_XYZ;
      UseSurfTarg = 1;
      DoProj = 1;
      nargsused = 2;
    } 
    else if (!strcmp(option, "--projfrac")) {
      if(nargc < 2) {
        argnerr(option,2);
      }
      ProjType = 2;
      SurfSrcName = pargv[0];
      sscanf(pargv[1],"%lf",&ProjDepth);
      UseSurfSrc  = SURF_SRC_XYZ;
      UseSurfTarg = 1;
      DoProj = 1;
      nargsused = 2;
    } 
    else if (!strcmp(option, "--proj-norm")) {
      // --proj-norm sourcesurf projdistmm outsurf
      if(nargc < 3) argnerr(option,2);
      MRIS *surf = MRISread(pargv[0]);
      if(surf == NULL) exit(1);
      double dist;
      sscanf(pargv[1],"%lf",&dist);
      int n;
      for(n=0; n < surf->nvertices; n++) {
	float x,y,z;
	ProjNormDist(&x, &y, &z,surf, n, dist);
	MRISsetXYZ(surf,n, x,y,z);
      }
      int err = MRISwrite(surf,pargv[2]);
      exit(err);
    } 
    else if (!strcmp(option, "--reshape-factor")) {
      if(nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&reshapefactortarget);
      reshape = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--sval-nxyz")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SurfSrcName = pargv[0];
      UseSurfSrc = SURF_SRC_NXYZ;
      nargsused = 1;
    } else if (!strcasecmp(option, "--sval-tal-xyz")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SurfSrcName = pargv[0];
      UseSurfSrc = SURF_SRC_TAL_XYZ;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--reg") || !strcasecmp(option, "--reg-inv")) {
      if(nargc < 1) argnerr(option,1);
      nargsused = 1;
      if(!stricmp(FileNameExtension(pargv[0], tmp), "LTA")) {
	if(CMDnthIsArg(nargc, pargv, 1) ) {
	  printf("ERROR: do not include a template volume with --reg or --reg-inv when using an LTA file\n");
	  exit(1);
	}
        LTA *lta ;
        lta = LTAread(pargv[0]) ;
        if (lta == NULL) return(1) ;
	if(lta->type != REGISTER_DAT){
	  printf("Converting LTA to REGISTER_DAT\n");
	  LTAchangeType(lta,REGISTER_DAT);
	}
        regsubject = (char *) calloc(strlen(lta->subject)+2,sizeof(char));
        strcpy(regsubject, lta->subject) ;
        intensity = lta->fscale ;
        float2int = FLT2INT_ROUND ;
        ipr = lta->xforms[0].src.xsize ;
        bpr = lta->xforms[0].src.zsize ;
        XFM = MatrixCopy(lta->xforms[0].m_L, NULL) ;
        RegTarg = MRIallocHeader(lta->xforms[0].src.width,
                                 lta->xforms[0].src.height,
                                 lta->xforms[0].src.depth,
                                 MRI_UCHAR,1) ;
        MRIcopyVolGeomToMRI(RegTarg, &lta->xforms[0].src) ;
        LTAfree(&lta) ;
        err = 0 ;
      } 
      else{
        err = regio_read_register(pargv[0], &regsubject, &ipr, &bpr,
                                  &intensity, &XFM, &float2int);
	if(err) exit(1);
	if(CMDnthIsArg(nargc, pargv, 1) ) {
	  printf("Reading header for %s\n",pargv[1]);
	  RegTarg = MRIreadHeader(pargv[1],MRI_VOLUME_TYPE_UNKNOWN);
	  if(RegTarg == NULL) exit(1);
	  nargsused ++;
	}
      }
      if(!strcasecmp(option, "--reg-inv")){
	printf("Inverting matrix\n");
	XFM = MatrixInverse(XFM,XFM);
      }
      ApplyReg = 1;
    } 
    else if (!strcasecmp(option, "--reg-diff")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      err = regio_read_register(pargv[0], &regsubject, &ipr, &bpr,
                                &intensity, &XFMSubtract, &float2int);
      if (err) {
        exit(1);
      }
      nargsused = 1;
    } else if (!strcasecmp(option, "--rms")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      RMSDatFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--rms-mask")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      RMSMaskFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--reg-inv-lrrev")) {
      // See docs below on MRIleftRightRevMatrix() for what this
      // is and how it works.
      if (nargc < 1) {
        argnerr(option,1);
      }
      err = regio_read_register(pargv[0], &regsubject, &ipr, &bpr,
                                &intensity, &XFM, &float2int);
      if (err) {
        exit(1);
      }
      XFM = MatrixInverse(XFM,NULL);
      sprintf(tmpstr,"%s/%s/mri/brain.mgz",SUBJECTS_DIR,regsubject);
      printf("%s\n",tmpstr);
      lrmri = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
      if(lrmri == NULL) {
        exit(1);
      }
      M = MRIleftRightRevMatrix(lrmri);
      M->rptr[1][1] = -1.0;
      M->rptr[1][4] = +1.0;
      XFM = MatrixMultiply(XFM,M,XFM);
      MatrixFree(&M);
      MRIfree(&lrmri);
      ApplyReg = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--sval-area")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SurfSrcName = pargv[0];
      UseSurfSrc = SURF_SRC_AREA;
      nargsused = 1;
    } else if (!strcasecmp(option, "--sval-annot")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SurfSrcName = "orig";
      AnnotFile = pargv[0];
      UseSurfSrc = SURF_SRC_ANNOT;
      printf("Setting mapmethod to nnf\n");
      mapmethod = "nnf";
      nargsused = 1;
    } else if (!strcasecmp(option, "--sval-rip")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SurfSrcName = pargv[0];
      UseSurfSrc = SURF_SRC_RIP;
      nargsused = 1;
    } else if (!strcmp(option, "--srcdump")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SrcDumpFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srcfmt") || !strcmp(option, "--sfmt") ||
               !strcmp(option, "--src_type")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      srctypestring = pargv[0];
      srctype = string_to_type(srctypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--srcicoorder")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&SrcIcoOrder);
      nargsused = 1;
    } else if (!strcmp(option, "--nsmooth-in")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&nSmoothSteps_Input);
      if (nSmoothSteps_Input < 1) {
        fprintf(stderr,"ERROR: number of smooth steps (%d) must be >= 1\n",
                nSmoothSteps_Input);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--nsmooth-out") ||
               !strcmp(option, "--nsmooth")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&nSmoothSteps);
      if (nSmoothSteps < 1) {
        fprintf(stderr,"ERROR: number of smooth steps (%d) must be >= 1\n",
                nSmoothSteps);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--fwhm-src")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%lf",&fwhm_Input);
      gstd_Input = fwhm_Input/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcmp(option, "--fwhm") ||
               !strcmp(option, "--fwhm-trg")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%lf",&fwhm);
      gstd = fwhm/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcmp(option, "--label-trg")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      LabelFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--label-src")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      LabelFile_Input = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--cortex")) {
      UseCortexLabel = 1;
    } else if (!strcmp(option, "--no-cortex")) {
      UseCortexLabel = 0;
    }

    /* -------- target value inputs ------ */
    else if (!strcmp(option, "--trgsubject")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      trgsubject = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trgicoorder")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&TrgIcoOrder);
      nargsused = 1;
    } else if (!strcmp(option, "--trgsurfval")  || !strcmp(option, "--tval") ||
	       !strcmp(option, "--trgval") || !strcmp(option, "--o")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      trgvalfile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--tval-xyz")) {
      if(nargc < 1 || CMDnthIsArg(nargc, pargv, 0) == 0){
	printf("ERROR: --tval-xyz flag needs 1 argument\n");
	printf("   FYI: --tval-xyz now requires a volume. See --help\n");
	exit(1);
      }
      TrgSurfVolFile = pargv[0];
      UseSurfTarg = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--trgdump")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      TrgDumpFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trgfmt") ||!strcmp(option, "--tfmt") ||
               !strcmp(option, "--trg_type")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      trgtypestring = pargv[0];
      trgtype = string_to_type(trgtypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--frame")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&framesave);
      nargsused = 1;
    } else if (!strcmp(option, "--cavtx")) {
      /* command-line vertex -- for debugging */
      if (nargc < 1) {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%d",&cavtx);
      nargsused = 1;
    } else if (!strcmp(option, "--hemi") || !strcmp(option, "--h")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      srchemi = pargv[0];
      trghemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srchemi")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      srchemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trghemi")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      trghemi = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--dual-hemi")) {
      UseDualHemi = 1; // Assume ?h.?h.surfreg file name, source only
    } else if (!strcmp(option, "--surfreg")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      srcsurfregfile = pargv[0];
      trgsurfregfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srcsurfreg")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      srcsurfregfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trgsurfreg")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      trgsurfregfile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--mapmethod")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      mapmethod = pargv[0];
      if (strcmp(mapmethod,"nnfr") && strcmp(mapmethod,"nnf")) {
        fprintf(stderr,"ERROR: mapmethod must be nnfr or nnf\n");
        exit(1);
      }
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--mul") ){
      if (nargc < 1)argnerr(option,1);
      if(! isdigit(pargv[0][0]) && pargv[0][0] != '-' && 
	 pargv[0][0] != '+' && pargv[0][0] != '.'){
        printf("ERROR: value passed to the --mul flag must be a number\n");
        printf("       If you want to multiply two images, use fscalc\n");
        exit(1);
      }
      sscanf(pargv[0],"%lf",&MultiplyVal);
      DoMultiply = 1;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--div") ){
      if (nargc < 1)argnerr(option,1);
      if(! isdigit(pargv[0][0]) && pargv[0][0] != '-' && 
	 pargv[0][0] != '+' && pargv[0][0] != '.'){
        printf("ERROR: value passed to the --div flag must be a number\n");
        printf("       If you want to divide two images, use fscalc\n");
        exit(1);
      }
      sscanf(pargv[0],"%lf",&MultiplyVal);
      if(MultiplyVal == 0){
	printf("ERROR: you can't divide by zero\n");
	exit(1);
      }
      MultiplyVal = 1.0/MultiplyVal;
      DoMultiply = 1;
      nargsused = 1;
    }
    else if (!strcmp(option, "--srchits")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SrcHitFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--srcdist")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      SrcDistFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trghits")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      TrgHitFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--trgdist")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      TrgDistFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--vtxmap")) {
      if (nargc < 1) {
        argnerr(option,1);
      }
      ResampleVtxMapFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--proj-surf")) {
      // --proj-surf surf projmagfile scale outsurf
      if(nargc < 3) argnerr(option,3);
      MRIS *surf0 = MRISread(pargv[0]);
      if(surf0==NULL) exit(1);
      MRI *projmag = MRIread(pargv[1]);
      if(projmag==NULL) exit(1);
      if(surf0->nvertices != projmag->width){
	printf("ERROR: dimension mismatch %d %d\n",surf0->nvertices,projmag->width);
	exit(1);
      }
      double projscale;
      sscanf(pargv[2],"%lf",&projscale);
      printf("proj scale %lf\n",projscale);
      int vtxno;
      MRIScomputeMetricProperties(surf0);
      for(vtxno=0; vtxno < surf0->nvertices; vtxno++){
	VERTEX *v = &(surf0->vertices[vtxno]);
	v->x += (projscale * MRIgetVoxVal(projmag,vtxno,0,0,0) * v->nx);
	v->y += (projscale * MRIgetVoxVal(projmag,vtxno,0,0,0) * v->ny);
	v->z += (projscale * MRIgetVoxVal(projmag,vtxno,0,0,0) * v->nz);
      }
      MRISwrite(surf0,pargv[3]);
      exit(0);
      nargsused = 1;
    } 
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option)) {
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      }
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
  printf("   --srcsubject source subject\n");
  printf("   --sval path of file with input values \n");
  printf("   --sval-xyz  surfname : use xyz of surfname as input \n");
  printf("   --projfrac surfname 0.5 : use projected xyz of surfname as input \n");
  printf("   --projabs  surfname 0.5 : use projected xyz of surfname as input \n");
  printf("   --sval-tal-xyz  surfname : use tal xyz of surfname as input \n");
  printf("   --sval-area surfname : use vertex area of surfname as input \n");
  printf("   --sval-annot annotfile : map annotation \n");
  printf("   --sval-nxyz surfname : use surface normals of surfname as input \n");
  printf("   --patch srcpatchfile targsurf ndilations\n");
  printf("   --sfmt   source format\n");
  printf("   --reg register.dat <volgeom> : apply register.dat to sval-xyz\n");
  printf("   --reg-inv register.dat <volgeom> : apply inv(register.dat) to sval-xyz\n");
  printf("   --srcicoorder when srcsubject=ico and src is .w\n");
  printf("   --trgsubject target subject\n");
  printf("   --trgicoorder when trgsubject=ico\n");
  printf("   --tval path of file in which to store output values (or use --o)\n");
  printf("   --tval-xyz volume: save tval as a surface file with source xyz (volume for geometry)\n");
  printf("   --tfmt target format\n");
  printf("   --trgdist distfile : save distance from source to target vtx\n");
  printf("   --s subject : use subject as src and target\n");
  printf("   --hemi       hemisphere : (lh or rh) for both source and targ\n");
  printf("   --srchemi    hemisphere : (lh or rh) for source\n");
  printf("   --trghemi    hemisphere : (lh or rh) for target\n");
  printf("   --dual-hemi  : assume source ?h.?h.surfreg file name\n");
  printf("   --jac  : turn on jacobian correction, needed when applying to area or volume \n");
  printf("   --surfreg    source and targ surface registration (sphere.reg)  \n");
  printf("   --srcsurfreg source surface registration (sphere.reg)  \n");
  printf("   --trgsurfreg target surface registration (sphere.reg)  \n");
  printf("   --mapmethod  nnfr or nnf\n");
  printf("   --frame      save only nth frame (with --trg_type paint)\n");
  printf("   --fwhm-src fwhmsrc: smooth the source to fwhmsrc\n");
  printf("   --fwhm-trg fwhmtrg: smooth the target to fwhmtrg\n");
  printf("   --nsmooth-in N  : smooth the input\n");
  printf("   --nsmooth-out N : smooth the output\n");
  printf("   --cortex : use ?h.cortex.label as a smoothing mask\n");
  printf("   --no-cortex : do NOT use ?h.cortex.label as a smoothing mask (default)\n");
  printf("   --label-src label : source smoothing mask\n");
  printf("   --label-trg label : target smoothing mask\n");

  printf("     --mul Mul : Multiply the input by the given value\n");
  printf("     --div Div : Divide the input by the given value\n");
  printf("   --reshape  reshape output to multiple 'slices'\n");
  printf("   --reshape-factor Nfactor : reshape to Nfactor 'slices'\n");
  printf("   --reshape3d : reshape fsaverage (ico7) into 42 x 47 x 83\n");
  printf("   --split : output each frame separately\n");
  printf("   --synth : replace input with WGN\n");
  printf("   --ones  : replace input with 1s\n");
  printf("   --normvar : rescale so that stddev=1 (good with --synth)\n");
  printf("   --seed seed : seed for synth (default is auto)\n");
  printf("   --prune - remove any voxel that is zero in any time point (for smoothing)\n");
  printf("   --no-prune - do not prune (default)\n");
  printf("   --proj-surf surf projmagfile scale outsurf : project vertices by mag*scale at each vertex\n");
  printf("   --proj-norm sourcesurf distmm outsurf : project vertices by distmm at each vertex\n");
  printf("\n");
  printf("   --reg-diff reg2 : subtract reg2 from --reg (primarily for testing)\n");
  printf("   --rms rms.dat   : save rms of reg1-reg2 (primarily for testing)\n");
  printf("   --rms-mask mask : only compute rms in mask (primarily for testing)\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");

}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
printf("\n");
printf("This program will resample one surface onto another. The source and\n");
printf("target subjects can be any subject in $SUBJECTS_DIR and/or the\n");
printf("icosahedron (ico). The source and target file formats can be anything\n");
printf("supported by mri_convert. The source format can also be a curvature\n");
printf("file or a paint (.w) file. The user also has the option of smoothing\n");
printf("on the surface.\n");
printf("\n");
printf("OPTIONS\n");
printf("\n");
printf("  --srcsubject subjectname\n");
printf("\n");
printf("    Name of source subject as found in $SUBJECTS_DIR or ico for icosahedron.\n");
printf("    The input data must have been sampled onto this subject's surface (eg,\n");
printf("    using mri_vol2surf)\n");
printf("\n");
printf("  --sval sourcefile\n");
printf("\n");
printf("    Name of file where the data on the source surface is located.\n");
printf("\n");
printf("  --sval-xyz     surfname\n");
printf("  --sval-tal-xyz surfname\n");
printf("  --sval-area    surfname\n");
printf("  --sval-nxyz    surfname\n");
printf("\n");
printf("    Use measures from the input surface as the source (instead of specifying\n");
printf("    a source file explicitly with --sval). --sval-xyz extracts the x, y, and\n");
printf("    z of each vertex. --sval-tal-xyz is the same as --sval-xyz, but applies\n");
printf("    the talairach transform from srcsubject/mri/transforms/talairach.xfm.\n");
printf("    --sval-area extracts the vertex area. --sval-nxyz extracts the surface\n");
printf("    normals at each vertex. See also --tval-xyz.\n");
printf("\n");
printf("  --projfrac surfname frac\n");
printf("  --projabs  surfname dist\n");
printf("\n");
printf("    Use xyz from surfname as the input, project it along the normal, and\n");
printf("    save new xyz surface. Eg, to create a new surface halfway between\n");
printf("    the white and the pial:\n");
printf("      mri_surf2surf --s subject --projfrac white +0.5 --tval lh.mid --hemi lh\n");
printf("    saves $SUBJECTS_DIR/subject/lh.mid.\n");
printf("\n");
printf("  --sval-annot annotfile\n");
printf("\n");
printf("    Map annotation file to the output. The target data will be saved as an\n");
printf("    annotation.\n");
printf("\n");
printf("  --sfmt typestring\n");
printf("\n");
printf("    Format type string. Can be either curv (for FreeSurfer curvature file),\n");
printf("    paint or w (for FreeSurfer paint files), or anything accepted by\n");
printf("    mri_convert. If no type string  is given, then the type is determined\n");
printf("    from the sourcefile (if possible). If curv is used, then the curvature\n");
printf("    file will be looked for in $SUBJECTS_DIR/srcsubject/surf/hemi.sourcefile.\n");
printf("\n");
printf("  --srcicoorder order\n");
printf("\n");
printf("    Icosahedron order of the source. Normally, this can be detected based\n");
printf("    on the number of verticies, but this will fail with a .w file as input.\n");
printf("    This is only needed when the source is a .w file.\n");
printf("\n");
printf("  --trgsubject subjectname\n");
printf("\n");
printf("    Name of target subject as found in $SUBJECTS_DIR or ico for icosahedron.\n");
printf("\n");
printf("  --trgicoorder order\n");
printf("\n");
printf("    Icosahedron order number. This specifies the size of the\n");
printf("    icosahedron according to the following table:\n");
printf("              Order  Number of Vertices\n");
printf("                0              12\n");
printf("                1              42\n");
printf("                2             162\n");
printf("                3             642\n");
printf("                4            2562\n");
printf("                5           10242\n");
printf("                6           40962\n");
printf("                7          163842\n");
printf("    In general, it is best to use the largest size available.\n");
printf("\n");
printf("  --tval targetfile\n");
printf("\n");
printf("    Name of file where the data on the target surface will be stored.\n");
printf("    BUG ALERT: for trg_type w or paint, use the full path.\n");
printf("\n");
printf("  --tval-xyz volume\n");
printf("\n");
printf("    Use this flag to indicate that the output (specified by --tval)\n");
printf("    will be a binary surface file This requires that --sval-xyz or\n");
printf("    --sval-tal-xyz was also specified. volume is a volume file in\n");
printf("    the target space; the volume geometry from this file is imbedded\n");
printf("    into the surface file. This is a good way to map the surface of\n");
printf("    one subject to an average (talairach) subject. Note: it will save\n");
printf("    targetfile as trgsubject/surf/targetfile unless targetfile has a\n");
printf("    path.\n");
printf("\n");
printf("  --tfmt typestring\n");
printf("\n");
printf("    Format type string. Can be paint or w (for FreeSurfer paint files) or curv\n");
printf("    or anything accepted by mri_convert. If no type string  is given, then the type\n");
printf("    is determined from the sourcefile (if possible). If using paint, w, or curv,\n");
printf("    see also --frame.\n");
printf("\n");
printf("  --hemi hemifield (lh or rh)\n");
printf("\n");
printf("  --surfreg registration_surface\n");
printf("\n");
printf("    If the source and target subjects are not the same, this surface is used\n");
printf("    to register the two surfaces. sphere.reg is used as the default. Don't change\n");
printf("    this unless you know what you are doing.\n");
printf("\n");
printf("  --mapmethod methodname\n");
printf("\n");
printf("    Method used to map from the vertices in one subject to those of another.\n");
printf("    Legal values are: nnfr (neighest-neighbor, forward and reverse) and nnf\n");
printf("    (neighest-neighbor, forward only). Default is nnfr. The mapping is done\n");
printf("    in the following way. For each vertex on the target surface, the closest\n");
printf("    vertex in the source surface is found, based on the distance in the\n");
printf("    registration space (this is the forward map). If nnf is chosen, then the\n");
printf("    the value at the target vertex is set to that of the closest source vertex.\n");
printf("    This, however, can leave some source vertices unrepresented in target (ie,\n");
printf("    'holes'). If nnfr is chosen, then each hole is assigned to the closest\n");
printf("    target vertex. If a target vertex has multiple source vertices, then the\n");
printf("    source values are averaged together. It does not seem to make much difference.\n");
printf("\n");
printf("  --fwhm-src fwhmsrc\n");
printf("  --fwhm-trg fwhmtrg (can also use --fwhm)\n");
printf("\n");
printf("    Smooth the source or target with a gaussian with the given fwhm (mm). This is\n");
printf("    actually an approximation done using iterative nearest neighbor smoothing.\n");
printf("    The number of iterations is computed based on the white surface. This\n");
printf("    method is similar to heat kernel smoothing. This will give the same\n");
printf("    results as --nsmooth-{in,out}, but automatically computes the the\n");
printf("    number of iterations based on the desired fwhm.\n");
printf("\n");
printf("  --nsmooth-in  niterations\n");
printf("  --nsmooth-out niterations  [note: same as --smooth]\n");
printf("\n");
printf("    Number of smoothing iterations. Each iteration consists of averaging each\n");
printf("    vertex with its neighbors. When only smoothing is desired, just set the\n");
printf("    the source and target subjects to the same subject. --smooth-in smooths\n");
printf("    the input surface values prior to any resampling. --smooth-out smooths\n");
printf("    after any resampling. See also --fwhm-src and --fwhm-trg.\n");
printf("\n");
printf("  --label-src sourcelabel\n");
printf("  --label-trg targetlabel\n");
printf("  --cortex\n");
printf("  --no-cortex\n");
printf("\n");
printf("    Only smooth within the given label. If --cortex is specified, then\n");
printf("    ?h.cortex.label will be used (this is created by automatically be\n");
printf("    recon-all). Even if you do not have a label of interest, it is\n");
printf("    recommended that you only smooth within cortex in order to prevent\n");
printf("    values from the medial wall from being smoothed into the\n");
printf("    surrounding cortical areas. At some point, this will be the\n");
printf("    default at which point you will have to use --no-cortex to turn it\n");
printf("    off.  This documentation will reflect the change. For --label-src\n");
printf("    and --label-trg, if you do not give it a full path, it will look\n");
printf("    in the subjects label dir. There is no need to specify both source\n");
printf("    and target unless you are smoothing on both (which you probably\n");
printf("    should not be doing).\n");
printf("\n");
printf("  --frame framenumber\n");
printf("\n");
printf("    When using paint/w output format, this specifies which frame to output. This\n");
printf("    format can store only one frame. The frame number is zero-based (default is 0).\n");
printf("\n");
printf("  --mul Mul\n");
printf("  --div Div\n");
printf("    Multiply or divide the input by the given value\n");
printf("\n");
printf("  --reshape\n");
printf("\n");
printf("    Force mri_surf2surf to save the output as multiple 'slices'; has\n");
printf("    no effect for paint/w output format. For ico, the output will\n");
printf("    appear to be a 'volume' with Nv/R colums, 1 row, R slices and Nf\n");
printf("    frames, where Nv is the number of vertices on the surface. For\n");
printf("    icosahedrons, R=6. For others, R will be the prime factor of Nv\n");
printf("    closest to 6 (can be changed with --reshape-factor). Reshaping is\n");
printf("    for logistical purposes (eg, in the analyze/nifti format the size\n");
printf("    of a dimension cannot exceed 2^15). Use this flag to prevent this\n");
printf("    behavior. This has no effect when the output type is paint. At one\n");
printf("    point, it was the default to reshape.\n");
printf("\n");
printf("  --reshape-factor Nfactor\n");
printf("\n");
printf("    Attempt to reshape to Nfactor 'slices' (will choose closest prime\n");
printf("    factor) Default is 6.\n");
printf("\n");
printf("  --reshape3d\n");
printf("\n");
printf("    Reshape fsaverage (ico7) into 42 x 47 x 83\n");
printf("\n");
printf("  --sd SUBJECTS_DIR\n");
printf("\n");
printf("    Set SUBJECTS_DIR on the command line.\n");
printf("\n");
printf("EXAMPLES:\n");
printf("\n");
printf("1. Resample a subject's thickness of the left cortical hemisphere on to a\n");
printf("   7th order icosahedron and save in analyze4d format:\n");
printf("\n");
printf("   mri_surf2surf --hemi lh --srcsubject bert\n");
printf("      --srcsurfval thickness --src_type curv\n");
printf("      --trgsubject ico --trgicoorder 7\n");
printf("      --trgsurfval bert-thickness-lh.img --trg_type analyze4d\n");
printf("\n");
printf("2. Resample data on the icosahedron to the right hemisphere of subject bert.\n");
printf("   Note that both the source and target data are stored in mgh format\n");
printf("   as 'volume-encoded suface' data.\n");
printf("\n");
printf("   mri_surf2surf --hemi rh --srcsubject ico --srcsurfval icodata-rh.mgh\n");
printf("      --trgsubject bert --trgsurfval ./bert-ico-rh.mgh\n");
printf("\n");
printf("3. Convert the surface coordinates of the lh.white of a subject to a\n");
printf("   (talairach) average (ie, a subject created by make_average_subject):\n");
printf("\n");
printf("   mri_surf2surf --s yoursubject --hemi lh --sval-tal-xyz white\n");
printf("     --trgsubject youraveragesubject --tval lh.white.yoursubject \n");
printf("     --tval-xyz $SUBJECTS_DIR/fsaverage/mri/orig.mgz\n");
printf("\n");
printf("   This will create youraveragesubject/surf/lh.white.yoursubject\n");
printf("\n");
printf("4. Convert the surface coordinates of the lh.white of a subject to \n");
printf("   the subject's functional space\n");
printf("\n");
printf("   mri_surf2surf --reg register.dat template.nii.gz --hemi lh \n");
printf("      --sval-xyz white --tval-xyz template.nii.gz --tval ./lh.white.func \n");
printf("      --s yoursubject\n");
printf("\n");
printf("   This will create lh.white.func in the current directory. template.nii.gz\n");
printf("   is a volume in the functional space. register.dat is the registration\n");
printf("   file between anatomical (target) and functional (movable) spaces. \n");
printf("   View result with:  freeview -v template.nii.gz -f lh.white.func\n");
printf("\n");
printf("   When using an LTA instead of a register.dat, do not include a target volume\n");
printf("\n");
printf("   mri_surf2surf --reg register.lta --hemi lh \n");
printf("      --sval-xyz white --tval-xyz template.nii.gz --tval ./lh.white.func \n");
printf("      --s yoursubject\n");
printf("\n");
printf("5. Extract surface normals of the white surface and save in a\n");
printf("   volume-encoded file:\n");
printf("\n");
printf("   mri_surf2surf --s yoursubject --hemi lh --sval-nxyz white\n");
printf("      --tval lh.white.norm.mgh\n");
printf("\n");
printf("   This will create youraveragesubject/surf/lh.white.yoursubject\n");
printf("\n");
printf("\n");
printf("6. Convert the annotation for one subject to the surface of another\n");
printf("\n");
printf("  mri_surf2surf --srcsubject subj1 --trgsubject subj2 --hemi lh \\\n");
printf("    --sval-annot $SUBJECTS_DIR/subj1/label/lh.aparc.annot \\\n");
printf("    --tval       $SUBJECTS_DIR/subj2/label/lh.subj1.aparc.annot\n");
printf("\n");
printf("   This will create $SUBJECTS_DIR/subj2/label/lh.subj1.aparc.annot.\n");
printf("   The --sval-annot flag will also change the map method to nnf so that\n");
printf("   the annot indices are not averaged. Note: this is not a substitute\n");
printf("   for running the cortical parcellation! The parcellations that it\n");
printf("   maps to the new subject may not be appropriate for that subject.\n");
printf("\n");
printf("BUG REPORTS: send bugs to analysis-bugs@nmr.mgh.harvard.edu. Make sure\n");
printf("    to include the version and full command-line and enough information to\n");
printf("    be able to recreate the problem. Not that anyone does.\n");
printf("\n");
printf("BUGS:\n");
printf("\n");
printf("  When the output format is paint, the output file must be specified with\n");
printf("  a partial path (eg, ./data-lh.w) or else the output will be written into\n");
printf("  the subject's anatomical directory.\n");
printf("\n");
  exit(1) ;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"\n");
  fprintf(fp,"setenv SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  fprintf(fp,"cd %s\n",cwd);
  fprintf(fp,"%s\n",cmdline);
  fprintf(fp,"\n");
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"srcsubject = %s\n",srcsubject);
  fprintf(fp,"srcval     = %s\n",srcvalfile);
  fprintf(fp,"srctype    = %s\n",srctypestring);
  fprintf(fp,"trgsubject = %s\n",trgsubject);
  fprintf(fp,"trgval     = %s\n",trgvalfile);
  fprintf(fp,"trgtype    = %s\n",trgtypestring);
  fprintf(fp,"srcsurfreg = %s\n",srcsurfregfile);
  fprintf(fp,"trgsurfreg = %s\n",trgsurfregfile);
  fprintf(fp,"srchemi    = %s\n",srchemi);
  fprintf(fp,"trghemi    = %s\n",trghemi);
  fprintf(fp,"frame      = %d\n",framesave);
  fprintf(fp,"fwhm-in    = %g\n",fwhm_Input);
  fprintf(fp,"fwhm-out   = %g\n",fwhm);
  fprintf(fp,"label-src  = %s\n",LabelFile_Input);
  fprintf(fp,"label-trg  = %s\n",LabelFile);
  fprintf(fp,"OKToRevFaceOrder  = %d\n",OKToRevFaceOrder);
  fprintf(fp,"UseDualHemi = %d\n",UseDualHemi);

  return;
}
/* --------------------------------------------- */
static void print_version(void)
{
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1) {
    fprintf(stdout,"ERROR: %s flag needs %d argument\n",option,n);
  } else {
    fprintf(stdout,"ERROR: %s flag needs %d arguments\n",option,n);
  }
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if (srcsubject == NULL) {
    fprintf(stdout,"ERROR: no source subject specified\n");
    exit(1);
  }
  if (srcvalfile == NULL && UseSurfSrc == 0 && PatchFile == NULL) {
    fprintf(stdout,"A source value path must be supplied\n");
    exit(1);
  }

  if(UseSurfSrc == 0 && PatchFile == NULL) {
    if ( strcasecmp(srctypestring,"w") != 0 &&
         strcasecmp(srctypestring,"curv") != 0 &&
         strcasecmp(srctypestring,"paint") != 0 ) {
      if (srctype == MRI_VOLUME_TYPE_UNKNOWN) {
        srctype = mri_identify(srcvalfile);
        if (srctype == MRI_VOLUME_TYPE_UNKNOWN) {
          fprintf(stdout,"ERROR: could not determine type of %s\n",srcvalfile);
          exit(1);
        }
      }
    }
  } else {
    if(srcvalfile != NULL) {
      printf("ERROR: cannot spec both --sval-xyz and --sval\n");
      exit(1);
    }
  }

  if (trgsubject == NULL) {
    fprintf(stdout,"ERROR: no target subject specified\n");
    exit(1);
  }
  if(trgvalfile == NULL && RMSDatFile == NULL) {
    fprintf(stdout,"A target value path must be supplied\n");
    exit(1);
  }

  if (UseSurfTarg == 0 && UseSurfSrc != SURF_SRC_ANNOT && RMSDatFile == NULL && PatchFile == NULL) {
    if ( strcasecmp(trgtypestring,"w") != 0 &&
         strcasecmp(trgtypestring,"curv") != 0 &&
         strcasecmp(trgtypestring,"paint") != 0 ) {
      if (trgtype == MRI_VOLUME_TYPE_UNKNOWN) {
        trgtype = mri_identify(trgvalfile);
        if (trgtype == MRI_VOLUME_TYPE_UNKNOWN) {
          fprintf(stdout,"ERROR: could not determine type of %s\n",trgvalfile);
          exit(1);
        }
      }
    }
  } 
  else {
    if (UseSurfSrc != SURF_SRC_XYZ && UseSurfSrc != SURF_SRC_TAL_XYZ &&
        UseSurfSrc != SURF_SRC_ANNOT && PatchFile == NULL) {
      printf("ERROR: must use --sval-xyz or --sval-tal-xyz with --tval-xyz\n");
      exit(1);
    }
  }

  if(srcsurfregfile == NULL) {
    srcsurfregfile = "sphere.reg";
  } else {
    printf("Source registration surface changed to %s\n",srcsurfregfile);
  }

  if(trgsurfregfile == NULL) {
    trgsurfregfile = "sphere.reg";
  } else {
    printf("Target registration surface changed to %s\n",trgsurfregfile);
  }

  if (srchemi == NULL) {
    fprintf(stdout,"ERROR: no hemifield specified\n");
    exit(1);
  }

  if (fwhm != 0 && nSmoothSteps != 0) {
    printf("ERROR: cannot specify --fwhm-out and --nsmooth-out\n");
    exit(1);
  }
  if (fwhm_Input != 0 && nSmoothSteps_Input != 0) {
    printf("ERROR: cannot specify --fwhm-in and --nsmooth-in\n");
    exit(1);
  }

  if(DoNormVar && !SynthPDF) {
    printf("WARNING: variance normalization turned on but not synthesizing\n");
  }

  if(UseCortexLabel) {
    sprintf(tmpstr,"%s.cortex.label",srchemi);
    LabelFile_Input = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s.cortex.label",trghemi);
    LabelFile = strcpyalloc(tmpstr);
  }

  return;
}

/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if (len < 2) {
    return(0);
  }

  if (flag[0] == '-' && flag[1] != '-') {
    return(1);
  }
  return(0);
}

/*---------------------------------------------------------------*/
int GetNVtxsFromWFile(const char *wfile)
{
  FILE *fp;
  int i,ilat, num, nvertices;
  int *vtxnum;
  float *wval;

  fp = fopen(wfile,"r");
  if (fp==NULL) {
    fprintf(stdout,"ERROR: Progname: GetNVtxsFromWFile():\n");
    fprintf(stdout,"Could not open %s\n",wfile);
    fprintf(stdout,"(%s,%d,%s)\n",__FILE__, __LINE__,__DATE__);
    exit(1);
  }

  fread2(&ilat,fp);
  fread3(&num,fp);
  vtxnum = (int *)   calloc(sizeof(int),   num);
  wval   = (float *) calloc(sizeof(float), num);

  for (i=0; i<num; i++) {
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
int GetNVtxsFromValFile(const char *filename, const char *typestring)
{
  //int err,nrows, ncols, nslcs, nfrms, endian;
  int nVtxs=0;
  int type;
  MRI *mri;

  printf("GetNVtxs: %s %s\n",filename,typestring);

  if (!strcmp(typestring,"curv")) {
    fprintf(stdout,"ERROR: cannot get nvertices from curv format\n");
    exit(1);
  }

  if (!strcmp(typestring,"paint") || !strcmp(typestring,"w")) {
    nVtxs = GetNVtxsFromWFile(filename);
    return(nVtxs);
  }

  type = string_to_type(typestring);
  mri = MRIreadHeader(filename, type);
  if (mri == NULL) {
    exit(1);
  }

  nVtxs = mri->width*mri->height*mri->depth;

  MRIfree(&mri);

  return(nVtxs);
}
/*---------------------------------------------------------------*/
int GetICOOrderFromValFile(const char *filename, const char *fmt)
{
  int nIcoVtxs,IcoOrder;

  nIcoVtxs = GetNVtxsFromValFile(filename, fmt);

  IcoOrder = IcoOrderFromNVtxs(nIcoVtxs);
  if (IcoOrder < 0) {
    fprintf(stdout,"ERROR: number of vertices = %d, does not mach ico\n",
            nIcoVtxs);
    exit(1);

  }

  return(IcoOrder);
}
/*---------------------------------------------------------------*/
int dump_surf(char *fname, MRIS *surf, MRI *mri)
{
  FILE *fp;
  float val;
  int vtxno,nnbrs;
  VERTEX *vtx;

  fp = fopen(fname,"w");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n",fname);
    exit(1);
  }
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    val = MRIFseq_vox(mri,vtxno,0,0,0); //first frame
    if (val == 0.0) {
      continue;
    }
    nnbrs = surf->vertices_topology[vtxno].vnum;
    vtx = &surf->vertices[vtxno];
    fprintf(fp,"%5d  %2d  %8.4f %8.4f %8.4f   %g\n",
            vtxno,nnbrs,vtx->x,vtx->y,vtx->z,val);
  }
  fclose(fp);
  return(0);
}



/*---------------------------------------------------------------
  See also mrisTriangleArea() in mrisurf.c
  ---------------------------------------------------------------*/
double MRISareaTriangle(double x0, double y0, double z0,
                        double x1, double y1, double z1,
                        double x2, double y2, double z2)
{
  double xx, yy, zz, a;

  xx = (y1-y0)*(z2-z0) - (z1-z0)*(y2-y0);
  yy = (z1-z0)*(x2-x0) - (x1-x0)*(z2-z0);
  zz = (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);

  a = 0.5 * sqrt( xx*xx + yy*yy + zz*zz );
  return(a);
}
/*------------------------------------------------------------*/
int MRIStriangleAngles(double x0, double y0, double z0,
                       double x1, double y1, double z1,
                       double x2, double y2, double z2,
                       double *a0, double *a1, double *a2)
{
  double d0, d1, d2, d0s, d1s, d2s;

  /* dN is the distance of the segment opposite vertex N*/
  d0s = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
  d1s = (x0-x2)*(x0-x2) + (y0-y2)*(y0-y2) + (z0-z2)*(z0-z2);
  d2s = (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1);
  d0 = sqrt(d0s);
  d1 = sqrt(d1s);
  d2 = sqrt(d2s);

  /* Law of cosines */
  *a0 = acos( -(d0s-d1s-d2s)/(2*d1*d2) );
  *a1 = acos( -(d1s-d0s-d2s)/(2*d0*d2) );
  *a2 = M_PI - (*a0 + *a1);

  return(0);
}
/*------------------------------------------------------------*/
MRI *MRISdiffusionWeights(MRIS *surf)
{
  MRI *w;
  int nnbrsmax, nnbrs, vtxno, vtxnonbr;
  double area, wtmp;

  /* count the maximum number of neighbors */
  nnbrsmax = surf->vertices_topology[0].vnum;
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++)
    if (nnbrsmax < surf->vertices_topology[vtxno].vnum) {
      nnbrsmax = surf->vertices_topology[vtxno].vnum;
    }
  printf("nnbrsmax = %d\n",nnbrsmax);

  MRIScomputeMetricProperties(surf); /* for area */

  w = MRIallocSequence(surf->nvertices, 1, 1, MRI_FLOAT, nnbrsmax);
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    area = MRISsumVertexFaceArea(surf, vtxno);
    nnbrs = surf->vertices_topology[vtxno].vnum;
    //printf("%d %6.4f   ",vtxno,area);
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      vtxnonbr = surf->vertices_topology[vtxno].v[nthnbr];
      wtmp = MRISdiffusionEdgeWeight(surf, vtxno, vtxnonbr);
      MRIFseq_vox(w,vtxno,0,0,nthnbr) = (float)wtmp/area;
      //printf("%6.4f ",wtmp);
    }
    //printf("\n");
    //MRISdumpVertexNeighborhood(surf,vtxno);
  }

  return(w);
}
/*----------------------------------------------------------------------*/
MRI *MRISdiffusionSmooth(MRIS *Surf, MRI *Src, double GStd, MRI *Targ)
{
  MRI *w, *SrcTmp;
  double FWHM;
  float wtmp,val,val0,valnbr;
  int vtxno, nthnbr, nbrvtxno, Niters, nthiter;
  double dt=1;

  if (Surf->nvertices != Src->width) {
    printf("ERROR: MRISdiffusionSmooth: Surf/Src dimension mismatch\n");
    return(NULL);
  }

  if (Targ == NULL) {
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth,
                            MRI_FLOAT, Src->nframes);
    if (Targ==NULL) {
      printf("ERROR: MRISdiffusionSmooth: could not alloc\n");
      return(NULL);
    }
  } else {
    if (Src->width   != Targ->width  ||
        Src->height  != Targ->height ||
        Src->depth   != Targ->depth  ||
        Src->nframes != Targ->nframes) {
      printf("ERROR: MRISdiffusionSmooth: output dimension mismatch\n");
      return(NULL);
    }
    if (Targ->type != MRI_FLOAT) {
      printf("ERROR: MRISdiffusionSmooth: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  /* Make a copy in case it's done in place */
  SrcTmp = MRIcopy(Src,NULL);

  /* Compute the weights */
  printf("Computing diffusion weights\n");
  w = MRISdiffusionWeights(Surf);

  printf("Starting iterations\n");
  FWHM = GStd*sqrt(log(256.0));
  Niters = (int)(((FWHM*FWHM)/(16*log(2)))/dt);
  printf("Niters = %d, dt=%g, GStd = %g, FWHM = %g\n",Niters,dt,GStd,FWHM);
  for (nthiter = 0; nthiter < Niters; nthiter ++) {
    //printf("Step = %d\n",nthiter); fflush(stdout);

    for (vtxno = 0; vtxno < Surf->nvertices; vtxno++) {
      nnbrs = Surf->vertices_topology[vtxno].vnum;

      for (frame = 0; frame < Targ->nframes; frame ++) {
        val0 = MRIFseq_vox(SrcTmp,vtxno,0,0,frame);
        val = val0;
        //printf("%2d %5d %7.4f   ",nthiter,vtxno,val0);
        for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
          nbrvtxno = Surf->vertices_topology[vtxno].v[nthnbr];
          valnbr = MRIFseq_vox(SrcTmp,nbrvtxno,0,0,frame) ;
          wtmp = dt*MRIFseq_vox(w,vtxno,0,0,nthnbr);
          val += wtmp*(valnbr-val0);
          //printf("%6.4f ",wtmp);
        }/* end loop over neighbor */
        //printf("   %7.4f\n",val);

        MRIFseq_vox(Targ,vtxno,0,0,frame) = val;
      }/* end loop over frame */

    } /* end loop over vertex */

    MRIcopy(Targ,SrcTmp);
  }/* end loop over smooth step */

  MRIfree(&SrcTmp);
  MRIfree(&w);

  return(Targ);
}
/*-------------------------------------------------------------
  MRISareNeighbors() - tests whether two vertices are neighbors.
  -------------------------------------------------------------*/
int MRISareNeighbors(MRIS *surf, int vtxno1, int vtxno2)
{
  int nnbrs, nthnbr, nbrvtxno;
  nnbrs = surf->vertices_topology[vtxno1].vnum;
  for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
    nbrvtxno = surf->vertices_topology[vtxno1].v[nthnbr];
    if (nbrvtxno == vtxno2) {
      return(1);
    }
  }
  return(0);
}
/*-------------------------------------------------------------
  MRIScommonNeighbors() - returns the vertex numbers of the two
  vertices that are common neighbors of the the two vertices
  listed.
  -------------------------------------------------------------*/
int MRIScommonNeighbors(MRIS *surf, int vtxno1, int vtxno2,
                        int *cvtxno1, int *cvtxno2)
{
  int nnbrs, nthnbr, nbrvtxno;

  *cvtxno1 = -1;
  nnbrs = surf->vertices_topology[vtxno1].vnum;
  for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
    nbrvtxno = surf->vertices_topology[vtxno1].v[nthnbr];
    if (nbrvtxno == vtxno2) {
      continue;
    }
    if (MRISareNeighbors(surf,nbrvtxno,vtxno2)) {
      if (*cvtxno1 == -1) {
        *cvtxno1 = nbrvtxno;
      } else {
        *cvtxno2 = nbrvtxno;
        return(0);
      }
    }
  }
  return(0);
}
/*-------------------------------------------------------------------------
  MRISdiffusionEdgeWeight() - computes the unnormalized weight of an
  edge needed for diffusion-based smoothing on the surface. The actual
  weight must be divided by the area of all the traingles surrounding
  the center vertex. See Chung 2003.
  ----------------------------------------------------------------------*/
double MRISdiffusionEdgeWeight(MRIS *surf, int vtxno0, int vtxnonbr)
{
  int cvtxno1=0, cvtxno2=0;
  VERTEX *v0, *vnbr, *cv1, *cv2;
  double a0, a1, btmp, ctmp, w;

  MRIScommonNeighbors(surf, vtxno0, vtxnonbr, &cvtxno1, &cvtxno2);

  v0   = &surf->vertices[vtxno0];
  vnbr = &surf->vertices[vtxnonbr];
  cv1  = &surf->vertices[cvtxno1];
  cv2  = &surf->vertices[cvtxno2];

  MRIStriangleAngles(cv1->x,cv1->y,cv1->z, v0->x,v0->y,v0->z,
                     vnbr->x,vnbr->y,vnbr->z, &a0,&btmp,&ctmp);
  MRIStriangleAngles(cv2->x,cv2->y,cv2->z, v0->x,v0->y,v0->z,
                     vnbr->x,vnbr->y,vnbr->z, &a1,&btmp,&ctmp);
  w = 1/tan(a0) + 1/tan(a1);
  return(w);
}
/*-------------------------------------------------------------------------
  MRISvertexSumFaceArea() - sum the area of the faces that the given
  vertex is part of. Make sure to run MRIScomputeMetricProperties()
  prior to calling this function.
  ----------------------------------------------------------------------*/
double MRISsumVertexFaceArea(MRIS *surf, int vtxno)
{
  int n, nfvtx;
  FACE *face;
  double area;

  area = 0;
  nfvtx = 0;
  for (n = 0; n < surf->nfaces; n++) {
    face = &surf->faces[n];
    if (face->v[0] == vtxno || face->v[1] == vtxno || face->v[2] == vtxno) {
      area += face->area;
      nfvtx ++;
    }
  }

  if (surf->vertices_topology[vtxno].vnum != nfvtx) {
    printf("ERROR: MRISsumVertexFaceArea: number of adjacent faces (%d) "
           "does not equal number of neighbors (%d)\n",
           nfvtx,surf->vertices_topology[vtxno].vnum);
    exit(1);
  }

  return(area);
}
/*-------------------------------------------------------------------------
  MRISdumpVertexNeighborhood()
  ----------------------------------------------------------------------*/
int MRISdumpVertexNeighborhood(MRIS *surf, int vtxno)
{
  int  n, nnbrs, nthnbr, nbrvtxno, nnbrnbrs, nthnbrnbr, nbrnbrvtxno;
  FACE *face;

  VERTEX const * const v0 = &surf->vertices[vtxno];
  nnbrs = surf->vertices_topology[vtxno].vnum;

  printf("  seed vtx %d vc = [%6.3f %6.3f %6.3f], nnbrs = %d\n",vtxno,
         v0->x,v0->y,v0->z,nnbrs);

  for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
    nbrvtxno = surf->vertices_topology[vtxno].v[nthnbr];
    VERTEX const * const v = &surf->vertices[nbrvtxno];
    printf("   nbr vtx %5d v%d = [%6.3f %6.3f %6.3f]    ",
           nbrvtxno,nthnbr,v->x,v->y,v->z);

    nnbrnbrs = surf->vertices_topology[nbrvtxno].vnum;
    for (nthnbrnbr = 0; nthnbrnbr < nnbrnbrs; nthnbrnbr++) {
      nbrnbrvtxno = surf->vertices_topology[nbrvtxno].v[nthnbrnbr];
      if (MRISareNeighbors(surf,vtxno,nbrnbrvtxno)) {
        printf("%5d ",nbrnbrvtxno);
      }
    }
    printf("\n");
  }

  printf("   faces");
  for (n = 0; n < surf->nfaces; n++) {
    face = &surf->faces[n];
    if (face->v[0] == vtxno || face->v[1] == vtxno || face->v[2] == vtxno) {
      printf("  %7.4f",face->area);

    }
  }
  printf("\n");

  return(0);
}

/*--------------------------------------------------------------------------*/
MRI *MRISheatkernel(MRIS *surf, double sigma)
{
  int vtxno, nbrvtxno, nnbrs, nnbrsmax;
  double K, Ksum, two_sigma_sqr, dx, dy, dz, d2;
  MRI *hk;

  two_sigma_sqr = 2*pow(sigma,2);

  // Count the maximum number of neighbors
  nnbrsmax = 0;
  for (vtxno=0; vtxno < surf->nvertices; vtxno++) {
    nnbrs = surf->vertices_topology[vtxno].vnum;
    if (nnbrsmax < nnbrs) {
      nnbrsmax = nnbrs;
    }
  }

  printf("2s2 = %g,s = %g\n",two_sigma_sqr,sigma);
  printf("max no of neighbors = %d\n",nnbrsmax);
  hk = MRIallocSequence(surf->nvertices, 1, 1, MRI_FLOAT, nnbrsmax+1);

  printf("Filling in heat kernel weights\n");
  for (vtxno=0; vtxno < surf->nvertices; vtxno++) {
    VERTEX_TOPOLOGY const * const cvtxt = &surf->vertices_topology[vtxno];
    VERTEX          const * const cvtx  = &surf->vertices         [vtxno];
    nnbrs = cvtxt->vnum;
    Ksum = 0;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      nbrvtxno = cvtxt->v[nthnbr];
      VERTEX const * const nbrvtx = &surf->vertices[nbrvtxno];
      dx = (cvtx->x - nbrvtx->x);
      dy = (cvtx->y - nbrvtx->y);
      dz = (cvtx->z - nbrvtx->z);
      d2 = dx*dx + dy*dy + dz*dz;
      K = exp(-d2/two_sigma_sqr);
      Ksum += K;
      MRIsetVoxVal(hk,vtxno,0,0,nthnbr,K);
    }
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      K = MRIgetVoxVal(hk,vtxno,0,0,nthnbr);
      MRIsetVoxVal(hk,vtxno,0,0,nthnbr,K/(Ksum+1.0)); // +1 for self
    }
    MRIsetVoxVal(hk,vtxno,0,0,nnbrs,1/(Ksum+1.0));  // self
  }

  printf("Done computing heat kernel weights\n");

  //MRIwrite(hk,"hk.mgh");
  return(hk);
}


/*--------------------------------------------------------------------------*/

MRI *MRIShksmooth(MRIS *Surf, MRI *Src, double sigma,
                  int nSmoothSteps, MRI *Targ)
{
  int nnbrs, nthstep, frame, vtx, nbrvtx, nthnbr;
  double val, w;
  MRI *SrcTmp, *hk;

  if (Surf->nvertices != Src->width) {
    printf("ERROR: MRIShksmooth: Surf/Src dimension mismatch\n");
    return(NULL);
  }

  if (Targ == NULL) {
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth,
                            MRI_FLOAT, Src->nframes);
    if (Targ==NULL) {
      printf("ERROR: MRIShksmooth: could not alloc\n");
      return(NULL);
    }
  } else {
    if (Src->width   != Targ->width  ||
        Src->height  != Targ->height ||
        Src->depth   != Targ->depth  ||
        Src->nframes != Targ->nframes) {
      printf("ERROR: MRIShksmooth: output dimension mismatch\n");
      return(NULL);
    }
    if (Targ->type != MRI_FLOAT) {
      printf("ERROR: MRIShksmooth: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  printf("Computing heat kernel\n");
  hk = MRISheatkernel(SrcSurfReg, sigma);
  //MRIwrite(hk,"hk.mgh");

  printf("Heat kernel smoothing with %d steps\n",nSmoothSteps);
  SrcTmp = MRIcopy(Src,NULL);
  for (nthstep = 0; nthstep < nSmoothSteps; nthstep ++) {
    //printf("Step = %d\n",nthstep); fflush(stdout);

    for (vtx = 0; vtx < Surf->nvertices; vtx++) {
      nnbrs = Surf->vertices_topology[vtx].vnum;

      for (frame = 0; frame < Targ->nframes; frame ++) {
        w   = MRIgetVoxVal(hk,vtx,0,0,nnbrs); // weight for center
        val = MRIFseq_vox(SrcTmp,vtx,0,0,frame); // val for center

        if(0 && vtx == 10000) {
          printf("0 v = %g, w = %g\n",val,w);
        }
        val *= w;

        for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
          nbrvtx = Surf->vertices_topology[vtx].v[nthnbr];
          w = MRIgetVoxVal(hk,vtx,0,0,nthnbr);
          val += w*MRIFseq_vox(SrcTmp,nbrvtx,0,0,frame) ;
          if(0 && vtx == 10000) {
            printf("%d v = %g, w = %g\n",nthnbr+1,val,w);
          }

        }/* end loop over neighbor */

        MRIFseq_vox(Targ,vtx,0,0,frame) = val;
      }/* end loop over frame */

    } /* end loop over vertex */

    MRIcopy(Targ,SrcTmp);
  }/* end loop over smooth step */

  MRIfree(&SrcTmp);

  return(Targ);
}


int DumpSurface(MRIS *surf, char *outfile)
{
  FILE *fp;
  int nnbrsmax, vtxno, nnbrs, nbrvtxno;

  printf("Dumping surface to %s\n",outfile);

  // Count the maximum number of neighbors
  nnbrsmax = 0;
  for (vtxno=0; vtxno < surf->nvertices; vtxno++) {
    nnbrs = surf->vertices_topology[vtxno].vnum;
    if (nnbrsmax < nnbrs) {
      nnbrsmax = nnbrs;
    }
  }

  fp = fopen(outfile,"w");
  if (fp == NULL) {
    printf("ERROR: cannot open %s for writing\n",outfile);
    return(1);
  }

  for (vtxno=0; vtxno < surf->nvertices; vtxno++) {
    nnbrs = surf->vertices_topology[vtxno].vnum;
    VERTEX_TOPOLOGY const * const cvtxt = &surf->vertices_topology[vtxno];
    VERTEX          const * const cvtx  = &surf->vertices         [vtxno];
    fprintf(fp,"%6d   %8.3f %8.3f %8.3f   %2d   ",
            vtxno+1,cvtx->x,cvtx->y,cvtx->z,nnbrs);
    for (nthnbr = 0; nthnbr < nnbrsmax; nthnbr++) {
      if (nthnbr < nnbrs) {
        nbrvtxno = cvtxt->v[nthnbr];
      } else {
        nbrvtxno = -1;
      }
      fprintf(fp,"%6d ",nbrvtxno+1);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  printf("Done dumping surface\n");
  return(0);
}


// This function computes a ras2ras matrix which accomplishes
// a left-right voxel flip. This was developed as part of a
// project to evaluate how well registering an anatomical to
// its left-right flipped self would do at inter-hemispheric
// registration. It is called above in conjunction with a
// register.dat matrix computed with something like:
//   fsl_rigid_register -r brain.mgz -i brain.mgz -o brain.lrrev12.mgz
//      -left-right-reverse -regmat lrrev12.reg -subject subjectid
// You can then
//   mri_surf2surf --s subject --hemi lh --sval-xyz white
//     --tval lh.lrrev12.white --tval-xyz --reg-inv-lrrev ../mri/lrrev12.reg
// To test:
//   tkmedit subject brain.lrrev12.mgz lh.lrrev12.white -aux brain.mgz
// The lh.lrrev12.white should be on the right side of the brain and
// aligned with the folds of brain.lrrev12.mgz.
//
// This surface can then be applied to map functional values from the
// right hemisphere to the left by specifying
//   mri_vol2vol --reg register.dat --hemi lh --surf lrrev12.white ...
//
// M = vox2ras * V * inv(vox2ras), where V is the vox2vox matrix
// that accomplishes a left-right reversal
MATRIX *MRIleftRightRevMatrix(MRI *mri)
{
  MATRIX *V, *K, *invK, *M;

  // V is a Vox2Vox matrix that performs a LR flip
  V = MatrixIdentity(4,NULL);
  V->rptr[1][1] = -1;
  V->rptr[1][4] = mri->width - 1;

  // Vox2RAS matrix
  K = MRIxfmCRS2XYZ(mri,0);
  invK = MatrixInverse(K,NULL);

  // M = K*V*inv(K)
  M = MatrixMultiply(K,V,NULL);
  M = MatrixMultiply(M,invK,M);

  MatrixFree(&V);
  MatrixFree(&K);
  MatrixFree(&invK);

  return(M);
}

/*-------------------------------------------------------
  Rescale so that stddev=var=1. Good for simulations.
  -------------------------------------------------------*/
int NormVar(MRI *mri, MRI *mask)
{
  int c,r,s,f;
  long N;
  double v,sum, mean, sum2, var, stddev, m;

  printf("Normalizing variance\n");

  N = 0;
  sum  = 0;
  sum2 = 0;
  for(c=0; c < mri->width; c++) {
    for(r=0; r < mri->height; r++) {
      for(s=0; s < mri->depth; s++) {
        if(mask) {
          m = MRIgetVoxVal(mask,c,r,s,0);
          if(m < 0.5) {
            continue;
          }
        }
        for(f=0; f < mri->nframes; f++) {
          v = MRIgetVoxVal(mri,c,r,s,f);
          sum  += v;
          sum2 += (v*v);
          N ++;
        }
      }
    }
  }
  mean = sum/N;
  var = (N*mean*mean - 2*mean*sum + sum2)/(N-1);
  stddev = sqrt(var);
  printf("sum = %lf, sum2 = %lf\n",sum,sum2);
  printf("n = %ld, mean = %lf, var = %lf, stddev = %lf\n",
         N,mean,var,stddev);
  for(c=0; c < mri->width; c++) {
    for(r=0; r < mri->height; r++) {
      for(s=0; s < mri->depth; s++) {
        for(f=0; f < mri->nframes; f++) {
          v = MRIgetVoxVal(mri,c,r,s,f);
          v /= stddev;
          MRIsetVoxVal(mri,c,r,s,f,v);
        }
      }
    }
  }
  return(0);
}
