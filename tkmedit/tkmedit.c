/*============================================================================
 Copyright (c) 1996 Martin Sereno and Anders Dale
                    and kevin helped
=============================================================================*/
#define TCL
#define TKMEDIT 
/*#if defined(Linux) || defined(sun) || defined(SunOS) */
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
/*#endif */
#include <fcntl.h>
#include "error.h"
#include "diag.h"
#include "utils.h"
#include "const.h"

#include "tkmedit.h"


#define SET_TCL_ENV_VAR 1
#include <tcl.h>
#include <tclDecls.h>
#include <tk.h>
#include <tix.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "MRIio.h"
#include "volume_io.h"
#include "rgb_image.h"
#include "fio.h"
#include "mri_conform.h"

#ifndef OPENGL
#define OPENGL
#endif

#include "proto.h"
#include "macros.h"

#define CURSOR_VAL   255

/* #define SQR(x)       ((x)*(x)) */
#define MATCH(A,B)   (!strcmp(A,B))
#define MATCH_STR(S) (!strcmp(str,S))

#define NUMVALS 256
#define MAXIM 256
#define MAXPTS 10000
#define MAXPARS 10
#define MAPOFFSET 0
#define CORONAL    0
#define HORIZONTAL 1
#define SAGITTAL   2
#define POSTANT   0
#define INFSUP    1
#define RIGHTLEFT 2  /* radiol */
#define NAME_LENGTH  STRLEN
#define MAX_DIR_DEPTH  30
#define TMP_DIR          "tmp"             /* relative to subjectsdir/pname */
#define TRANSFORM_DIR    "mri/transforms"  /* ditto */
#define TALAIRACH_FNAME  "talairach.xfm"   /* relative to TRANSFORM_DIR */
#define DIR_FILE "ic1.tri"
#define MAXCOR 500
#define MAXLEN 100
#define kDefaultWindowLocationX 0
#define kDefaultWindowLocationY 0
#define kWindowBottomBorderHeight 32
#define CVIDBUF 25

#ifdef IRIX
int __new_cfgetospeed () {
  fprintf( stderr, "__new_cfgetospeed\n" );
}
int __new_cfsetospeed () {
  fprintf( stderr, "__new_cfsetospeed\n" );
}
int __new_cfsetispeed () {
  fprintf( stderr, "__new_cfsetispeed\n" );
}
int __new_tcgetattr () {
  fprintf( stderr, "__new_tcgetattr\n" );
}
int __new_tcsetattr () {
  fprintf( stderr, "__new_tcsetattr\n" );
}
#endif


int xnum=256,ynum=256;
int ptype;
float ps,st,xx0,xx1,yy0,yy1,zz0,zz1;
int zf, ozf;
float fsf;
static int xdim,ydim;
unsigned long bufsize;
unsigned char **im_b[MAXIM];
unsigned char **fill[MAXIM];
unsigned char **dummy_im[MAXIM];
unsigned char **sim[6]; 
unsigned char **sim2[6]; 
int second_im_allocated = FALSE;
int dummy_im_allocated = FALSE;
int wmfilter_ims_allocated = FALSE;
int changed[MAXIM];
int imnr0,imnr1,numimg;
int wx0=114,wy0=302;  /* (100,100), (117,90), (556,90) */
int ptsflag = FALSE;
int surfflag = FALSE;
int surfloaded = FALSE;
int editflag = TRUE;
int fieldsignflag = FALSE; /* overrides curvflag */
int surflinewidth = 1;
int curvloaded = FALSE;
int curvflag = FALSE;
int editedimage = FALSE;
int inplaneflag = TRUE;
int linearflag = FALSE;
int bwflag = FALSE;
int truncflag = FALSE;
int second_im_full = FALSE;
int npts = 0;
int ndip = 0;
int dip_spacing = 10; /* voxels */
float tm[4][4];
int jold[MAXPTS],iold[MAXPTS],imold[MAXPTS];
float ptx[MAXPTS],pty[MAXPTS],ptz[MAXPTS];

float par[MAXPARS],dpar[MAXPARS];

int white_lolim = 80;
int white_hilim = 140;
int gray_hilim = 100;
int flossflag = TRUE;
int spackleflag = TRUE;
int lim3=170,lim2=145,lim1=95,lim0=75;
double ffrac3=1.0,ffrac2=1.0,ffrac1=1.0,ffrac0=1.0;


#include "xDebug.h"
#include "xTypes.h"
#include "xUtilities.h"
#include "xVoxel.h"
#include "xList.h"
#include "x3DList.h"
#include "tkmMeditWindow.h"


                                       /* irix compiler doesn't like inlines */
#ifdef IRIX
#define inline
#endif


// ==================================================================== OUTPUT

/* for regular output. this could be remapped to display to a file or some 
   window widget like a status bar. that'd be cool. */
#define InitOutput
#define DeleteOutput
#define OutputPrint            fprintf ( stdout,
#define EndOutputPrint         );

// ===========================================================================

// ============================================================= FILENAME MGMT

char gsSubjectHomeDir[256] = ""; /* is $SUBJECTS_DIR/subject_name if they used
            the "subject image" argument, or the -f
            arg if they used that, or the cwd */
char gsUserHomeDir[256]    = ""; /* cwd from shell */
char gsTclScriptDir[256]   = ""; /* $MRI_DIR */

tBoolean SetSubjectHomeDirFromEnv ( char* isSubject );
void     SetSubjectHomeDir        ( char* isHomeDir );
tBoolean FindUserHomeDir          ();

typedef enum {
  tkm_tFileName_Functional = 0,
  tkm_tFileName_Parcellation,
  tkm_tFileName_HeadPoints,
  tkm_tFileName_Surface,
  tkm_tFileName_Volume,
  tkm_tFileName_Label,
  tkm_tFileName_RGB,
  tkm_tFileName_ControlPoints,
  tkm_tFileName_Edit,
  tkm_tFileName_TclScript,
  tkm_knNumFileNameTypes
} tkm_tFileName;

/* subdirectories local to subject's home dir */
char ksaFileNameSubDirs[tkm_knNumFileNameTypes][256] = {
  "fmri", "", "bem", "surf", "mri", "label", "image/rgb", "tmp", "tmp", "lib/tcl"
};

/* input starts with ., gsUserHomeDir will be prepended. if ~ or nothing,
   gsSubjectHomeDir will be prepended. anything else, first use the subject
   home dir and then the subdir of the proper type, then the proper
   remainder of the input file name. if it starts with a /, nothing will be
   added. */
void MakeFileName ( char*          isInput,
        tkm_tFileName  iType, 
        char*          osCompleteFileName);

/* attempts to extract a name from the data */
void ExtractSubjectName ( char* isDataSource,
        char* osSubjectName );
void ExtractVolumeName  ( char* isDataSource,
        char* osVolumeName );

// ==========================================================================

// ==================================================== COORDINATE CONVERSIONS

#include "mriTransform.h"

mriTransformRef gRASTransform = NULL;

void InitTransformation ();

                                      /* converts ras coords to
                                          voxel coords and back */
void RASToVoxel ( Real x, Real y, Real z,        // incoming ras coords
                  int *xi, int *yi, int *zi );   // outgoing voxel coords

void VoxelToRAS ( int xi, int yi, int zi,        // incoming voxel coords
                  Real *x, Real *y, Real *z );   // outgoing RAS coords

// ===========================================================================

// =========================================================== BOUNDS CHECKING

inline char IsVoxelInBounds ( int x, int y, int z );

// ===========================================================================

// ================================================== SELECTING CONTROL POINTS

/* returns nearest control point on the same plane. */
tBoolean FindNearestCtrlPt ( xVoxelRef        iAnaIdx, 
           mri_tOrientation iPlane,
           xVoxelRef*       opCtrlPt );

/* find the nearest ctrl pt and add or remove its pointer to the selection
   list. note that this does not allocate or delete points. only a ptr to
   ctrl pt lives in the selection list */
void AddNearestCtrlPtToSelection      ( xVoxelRef        iAnaIdx, 
          mri_tOrientation iPlane );
void RemoveNearestCtrlPtFromSelection ( xVoxelRef        iAnaIdx, 
          mri_tOrientation iPlane );
void SelectCtrlPt                     ( xVoxelRef        iCtrlPt );
void DeselectCtrlPt                   ( xVoxelRef        iCtrlPt );

/* deselect all control points */
void DeselectAllCtrlPts               ();

/* gets every ctrl pt in the selection list and passes it to DeleteCtrlPt.
   clears the selection list */
void DeleteSelectedCtrlPts            ();



/* makes a copy of the ctrl pt and puts it in the ctrl pt list */
void NewCtrlPt           ( xVoxelRef iCtrlPt,
         tBoolean  ibWriteToFile );

/* removes ctrl pt from the list and deletes it */
void DeleteCtrlPt        ( xVoxelRef ipCtrlPt );

/* passes the cursor to NewCtrlPt */
void NewCtrlPtFromCursor ();


/* reads the control.dat file, transforms all pts from RAS space 
   to voxel space, and adds them as control pts */
void ProcessCtrlPtFile ();
tBoolean gParsedCtrlPtFile = FALSE;

/* writes all control points to the control.dat file in RAS space */
void WriteCtrlPtFile   ();


x3DListRef gCtrlPtList = NULL;
xListRef gSelectionList = NULL;

/* function for comparing voxels */
xList_tCompare CompareVoxels ( void* inVoxelA, void* inVoxelB );

// ===========================================================================

// ================================================================== SURFACES

#include "mriSurface.h"

mriSurfaceRef gSurface = NULL;

void LoadSurface          ( char*           isName );
void LoadSurfaceVertexSet ( Surf_tVertexSet iSet,
          char*           fname );
void UnloadSurface        ();

// ===========================================================================

// ========================================================= SELECTING REGIONS

/* selecting regions works much like editing. the user chooses a brush shape and size and paints in a region. the tool can be toggled between selecting and unselecting. the selected pixels are kept in a voxel space for optimized retreival in the draw loop. there are the usual functions for adding and removing voxels as well as saving them out to a file. */

x3DListRef gSelectedVoxels;
char isDisplaySelectedVoxels;

void InitSelectionModule ();
void DeleteSelectionModule ();

/* grabs the list of voxels selected and draws them into the buffer. */
void DrawSelectedVoxels ( char * inBuffer, int inPlane, int inPlaneNum );

/* handles clicks. uses the current brush settings to paint or paint selected voxels. */
void AllowSelectionModuleToRespondToClick (xVoxelRef inScreenVoxel );

  /* adds or removes voxels to selections. if a voxel that isn't in the 
     selection is told to be removed, no errors occur. this is called from the 
     brush function. */
void AddVoxelToSelection      ( xVoxelRef  iAnaIdx );
void RemoveVoxelFromSelection ( xVoxelRef  ipAnaIdx );

/* clears the current selection */
void ClearSelection ();

/* are we currently displaying the selection? */
char IsDisplaySelectedVoxels ();

/* write to and read from label files */
void SaveSelectionToLabelFile ( char * inFileName );
void LoadSelectionFromLabelFile ( char * inFileName );

/* send the selected voxels to the functional display module */
void GraphSelectedRegion ();

/* values used in calculating selection color */
#define kSel_IntensifyValue 100
#define kSel_TooClose       10

// ===========================================================================

// ============================================================ EDITING VOXELS

void EditVoxelInRange(xVoxelRef ipVoxel, 
           tVolumeValue inLow, tVolumeValue inHigh, 
           tVolumeValue inNewValue );

// ===========================================================================

// ============================================================= VOLUME ACCESS

#define knNumVolumeValues 256
#define knVolumeSize 256
static int gVolumeDimension;

static tVolumeRef gAnatomicalVolume    = NULL;
static tVolumeRef gAuxAnatomicalVolume = NULL;

/* talairach transforms for each volume */
General_transform gAnaToTalTransform;
General_transform gAuxAnaToTalTransform;
tBoolean          gbAnaToTalTransformLoaded    = FALSE;
tBoolean          gbAuxAnaToTalTransformLoaded = FALSE;

char gSubjectName[256]            = ""; /* i.e. anders */
char gAuxSubjectName[256]         = ""; 
char gSubjectVolumeName[256]      = ""; /* i.e. T1 */
char gAuxSubjectVolumeName[256]   = "";
char gVolumeSrc[256]              = ""; /* i.e. /home/users/../anders/mri/T1 */
char gAuxVolumeSrc[256]           = "";

void          InitVolume            ( tVolumeRef*        ioVolume, 
              int                inDimension );
void          DeleteVolume          ( tVolumeRef*        ioVolume );
tVolumeValue* GetVolumeSlicePtr     ( tVolumeRef         inVolume,
              int                inSlice );
int           ReadVolumeWithMRIRead ( tVolumeRef*        ioVolume, 
              General_transform* iTransform,
              char*              inFileOrPath );

/* access */
inline tVolumeValue GetVoxelValue ( tVolumeRef   inVolume,
            int          x, 
            int          y, 
            int          z );
inline void         SetVoxelValue ( tVolumeRef   inVolume,
            int          x, 
            int          y, 
            int          z, 
            tVolumeValue inValue );

/* color scale management */
#define kfDefaultVolumeThreshold 0.35
#define kfDefaultVolumeSquash    12.0

static float gfVolumeColorThreshold[tkm_knNumVolumeTypes];
static float gfVolumeColorSquash[tkm_knNumVolumeTypes];
static float gfaVolumeColors[tkm_knNumVolumeTypes][knNumVolumeValues];

void SetVolumeColorScale ( tkm_tVolumeType iVolume,
         float           ifThreshold, 
         float           ifSquash );
void GetVolumeColor      ( tkm_tVolumeType iVolume,
         tVolumeValue    iucValue,
         xColor3fRef     oColor );


/* maximum intensity projection for each volume */
static tVolumeRef gAnatomicalMaxIntProj;
static tVolumeRef gAuxAnatomicalMaxIntProj;

void         BuildVolumeMaxIntProj    ( tVolumeRef       iVolume, 
          tVolumeRef*      iopMaxIntProjVolume );
tVolumeValue GetVolumeMaxIntProjValue ( tVolumeRef       iMaxIntProjVolume, 
          mri_tOrientation iOrientation,
          xVoxelRef        pVoxel );

/* snapshot */
static tVolumeRef gSnapshotVolume = NULL;

void SnapshotVolume            ();
void RestoreVolumeFromSnapshot ();

/* flipping */
void ChangeVolumeAxes( tVolumeRef iVolume,
           char*      isNewXAxis,
           char*      isNewYAxis,
           char*      isNewZAxis );

// ===========================================================================

// ========================================================= FUNCTIONAL VOLUME

#include "tkmFunctionalVolume.h"

tkmFunctionalVolumeRef gFunctionalVolume = NULL;

void LoadFunctionalOverlay    ( char* isPathAndStem );
void LoadFunctionalTimeCourse ( char* isPathAndStem );

// ===========================================================================

// ================================================ PARCELLATION and ROI GROUPS

#include "mriColorLookupTable.h"
#include "mriROIGroup.h"

static mriROIGroupRef         gROIGroup    = NULL;
static mriColorLookupTableRef gColorTable  = NULL;

void LoadParcellationVolume ( char* inVolumeDirWithPrefix,
            char* inColorFileName );

void GetROIColorAtVoxel      ( xVoxelRef   inVoxel, 
             xColor3fRef oColor );
void GetROILabel             ( xVoxelRef   ipVoxel, 
             int*        onIndex,
             char*       osLabel );
void SelectCurrentROI        ( int         inIndex );

// ===========================================================================

// ============================================================== EDITING UNDO

#include "xUndoList.h"

                                   /* this is a pretty simple implementation
              of undo that only supports pixel 
              editing. when editing, the pixels that
              were changed are saved along with their
              previous values in a list, one list per
              editing click. */

void InitUndoList ();
void DeleteUndoList ();

                                   /* note that the list is cleared when the
              mouse button 2 or 3 is pressed down
              right from the event handling code. 
              this is a hack, but it's the best we
              can do until we pass events and not
              just coords to ProcessClick. */
void ClearUndoList ();

                                   /* when pixels are editied, they are added
              to the list. if the user hits undo, the
              entire list is drawn to the screen, 
              using the SetVoxelValue() function. at
              the same time, a new list is made
              to save the positions of all the restored
              voxel values. that list becomes the
              new undo list, so you can effectivly
              undo an undo. */
void AddVoxelAndValueToUndoList (xVoxelRef inVoxel, int inValue );
void RestoreUndoList ();

                                   /* we need a struct for the undo list. this
              is what we add to it and what we get
              back when the list is restored. */
typedef struct {
 xVoxelRef mVoxel;
  tVolumeValue mValue;
} UndoEntry, *UndoEntryRef;

void NewUndoEntry           ( UndoEntryRef* outEntry, 
           xVoxelRef inVoxel, tVolumeValue inValue );
void DeleteUndoEntry        ( UndoEntryRef* ioEntry );

                                   /* these are our callback functions for the
              undo list. the first deletes an entry
              and the second actually performs the
              undo action and hands back the undone
              voxel. */
void DeleteUndoEntryWrapper ( xUndL_tEntryPtr* inEntryToDelete );
void UndoActionWrapper      ( xUndL_tEntryPtr  inUndoneEntry, 
            xUndL_tEntryPtr* outNewEntry );
void PrintEntryWrapper      ( xUndL_tEntryPtr  inEntry );

xUndoListRef gUndoList = NULL;

// ==========================================================================

/* =========================================================== HEAD POINTS */

#include "mriHeadPointList.h"

mriHeadPointListRef gHeadPoints = NULL;

void LoadHeadPts ( char* isHeadPtsFile, 
       char* isTransformFile );

void RestoreHeadPoints        ();
void WriteHeadPointsTransform ();
void WriteHeadPointsFile      ();
void RotateHeadPts            ( float ifDegrees, 
        char  isDirection );
void TranslateHeadPts         ( float ifDistance, 
        char  isDirection );

void SetSelectedHeadPointLabel ( char* isNewLabel );

void AlignSelectedHeadPointToAnaIdx ( xVoxelRef iAnaIdx );

/* ======================================================================= */

/* ========================================================== MEDIT WINDOW */

tkmMeditWindowRef gMeditWindow = NULL;
int gDoRedraw = 1;

/* ======================================================================= */


// ====================================================================== MISC


                                   /* where to find the interface script */
char gInterfaceScriptName [256] = "";

/* the script to run in batch mode, i.e from te -tlc switch */
tBoolean gbRunBatchScript = FALSE;
char gsBatchScriptName[256] = "";

                                   /* set and get the tcl interp to send
            the msg to */
void SetTclInterp ( Tcl_Interp * inInterp );
Tcl_Interp * GetTclInterp ();

                                   /* send a tcl command */
void SendTCLCommand ( char * inCommand );


#define kMaxNumCachedCommands 100
char gCachedTclCommands[kMaxNumCachedCommands][256];
int gNumCachedCommands = 0;
void SendCachedTclCommands ();
                                   /* the tcl interpreter */
Tcl_Interp * gTclInterp = NULL;

                                   /* don't start accepting tcl commands
              until setup is complete */
static tBoolean gbAcceptingTclCommands = FALSE;

                                   /* use a limited tcl/tk interface */
static tBoolean gbUseCsurfInterface = FALSE;

// ===========================================================================


/*--------------------- prototypes ------------------------------*/
#ifdef Linux
extern void scale2x(int, int, unsigned char *);
#endif


int Medit(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[]);
int read_second_images(char *imdir2) ;
void WriteVoxelToControlFile ( xVoxelRef inVolumeVox );
void WriteVoxelToEditFile    ( xVoxelRef inVolumeVox );
void mri2pix(float xpt, float ypt, float zpt, int *jpt, int *ipt,int *impt);

void rotate_brain(float a,char c) ;
void translate_brain(float a,char c) ;
void optimize2(void) ;
void optimize(int maxiter) ;
float Error(int p,float dp) ;
int  imval(float px,float py,float pz) ; // get value of im
void UpdateAndRedraw ();
void redraw(void) ; // post a glut redisplay
void pix_to_rgb(char *fname) ; // another method of saving a screenshot
void scrsave_to_rgb(char *fname) ; // use scrsave to save a screen shot
void save_rgb(char *fname) ; // saves a screen shot
//void edit_pixel(int action);
int write_images(char *fpref) ; // saves volume to COR files
int read_images(char *fpref) ; // reads volume from COR files
void write_dipoles(char *fname) ;
void write_decimation(char *fname) ;
void read_hpts(char *fname) ;
void read_htrans(char *fname) ;
void write_htrans(char *fname) ;
void smooth_3d(int niter) ;
void wmfilter_corslice(int imc) ;
void sagnorm_allslices(void) ;
void sagnorm_corslice(int imc) ;
void alloc_second_im(void) ;
void smooth_surface(int niter) ;

char *Progname ;
/*--------------------- end prototypes ------------------------------*/

/*--------------------- twitzels hacks ------------------------------*/

#ifdef USE_LICENSE

extern char *crypt(const char *, const char *) ;
/* Licensing */
void checkLicense(char* dirname)
{
  FILE* lfile;
  char* email;
  char* magic;
  char* key;
  char* gkey;
  char* lfilename;

  lfilename = (char*)malloc(512);
  email = (char*)malloc(512);
  magic = (char*)malloc(512);
  key = (char*)malloc(512);
  gkey = (char*)malloc(1024);

  sprintf(lfilename,"%s/.license",dirname);

  lfile = fopen(lfilename,"r");
  if(lfile) {
    fscanf(lfile,"%s\n",email);
    fscanf(lfile,"%s\n",magic);
    fscanf(lfile,"%s\n",key);
    
    sprintf(gkey,"%s.%s",email,magic);
    if (strcmp(key,crypt(gkey,"*C*O*R*T*E*C*H*S*0*1*2*3*"))!=0) {
      printf("No valid license key !\n");
      exit(-1);
    }
  }
  else {
    printf("License file not found !\n");
    exit(-1);
  }
  free(email);
  free(magic);
  free(key);
  free(gkey);
  free(lfilename);
  return;  
}
#endif

char hacked_map[256];


/*--------------------- twitzels hack end ---------------------------*/

int Medit(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[]) {

  FunV_tErr               eFunctional                 = FunV_tErr_NoError;
  int                     nCurrentArg                 = 0;
  int                     bFatalError                 = FALSE;
  char                    sArg[128]                   = "";
  tBoolean                bSubjectDeclared            = FALSE;
  tBoolean                bUsingMRIRead               = FALSE;
  char                    sSubject[128]               = "";
  char                    sImageDir[128]              = "";
  tBoolean                bSurfaceDeclared            = FALSE;
  char                    sSurface[128]               = "";
  tBoolean                bLocalImageDir              = FALSE;
  tBoolean                bNoEdit                     = FALSE;
  tBoolean                bLoadingAuxVolume           = FALSE;
  char                    sAuxVolume[128]             = "";
  tBoolean                bLoadingOverlay             = FALSE;
  char                    sOverlayPathAndStem[256]    = "";
  tBoolean                bLoadingTimeCourse          = FALSE;
  char                    sTimeCoursePathAndStem[256] = "";
  tBoolean                bLoadingParcellation        = FALSE;
  char                    sParcellationPath[256]      = "";
  char                    sParcellationColorFile[256] = "";
  tBoolean                bThresh                     = FALSE;
  FunV_tFunctionalValue   min                         = 0;
  tBoolean                bMid                        = FALSE;
  FunV_tFunctionalValue   mid                         = 0;
  tBoolean                bSlope                      = FALSE;
  FunV_tFunctionalValue   slope                       = 0;
  tBoolean                bRevPhaseFlag               = FALSE;
  int                     nRevPhaseFlag               = 0;
  tBoolean                bTruncPhaseFlag             = FALSE;
  int                     nTruncPhaseFlag             = 0;
  tBoolean                bUseOverlayCacheFlag        = FALSE;
  int                     nUseOverlayCacheFlag        = 0;
  tBoolean                bLoadingHeadPts             = FALSE;
  tBoolean                bHaveHeadPtsTransform       = FALSE;
  char                    sHeadPts[256]               = "";
  char                    sHeadPtsTransform[256]      = "";

  /* first get the functional threshold so we don't overwrite the defaults */
  eFunctional = FunV_GetThreshold( gFunctionalVolume, &min, &min, &slope );
  if( FunV_tErr_NoError != eFunctional ) {
    DebugPrint "Medit(): Couldn't get functional volume threshold.\n" 
      EndDebugPrint;
  }

    if (argc<2) {

      printf("tkmedit: integrated volumetric data viewer and editor with surface overlay\n\n");
printf("usage: tkmedit {[subject image_type] | [-f path_to_data]} [surface]\n");
printf("       [-interface script]\n");
printf("       [-parcellation path_to_data color_file] \n");
printf("       [-headpts points_file [transformation_file]] \n");
printf("       [-overlay path/stem] [-timecourse path/stem] \n");
printf("       [-fthresh min_overlay_threshold] [-fmid overlay_threshold_midpoint] \n");
printf("       [-fslope overlay_threshold_slope] [-revphaseflag 1|0] \n");
printf("       [-truncphaseflag 1|0] [-overlaycache 1|0]\n\n");
printf("   subject image_type : reads main subject volume as COR- file in \n");
printf("                        $SUBJECTS_DIR/subject/mri/image_type/\n\n");
printf("   f path : reads main subject volume if other than COR- file, or if not\n");
printf("            in normal $SUBJECTS_DIR path.\n\n");
printf("   surface : reads in surface in $SUBJECTS_DIR/subject/surf/surface\n\n");
printf("   interface script : specify tcl script to use for interface\n");
printf("   parcellation path_to_data color_file : load volume in path as\n");
printf("                                          parcellation volume with color file\n\n");
 printf("   headpts points_file [transformation_file] : read a file of head points\n");
 printf("                 and a transform file. if the points file ends in .hpts\n");
 printf("                 you can omit the transform file and tkmedit will look for\n");
 printf("                 a transform file with the same name as the pts file with\n");
 printf("                 the suffix .trans instead.\n\n");
printf("   overlay path/stem : load volume in path with stem as overlay volume. must\n");
printf("                       be bfile format. i.e. if data is /path/h.bfloat, arg\n");
printf("           should be /path/h\n\n");
printf("   timecourse path/stem : load volume in path with stem as time course volume.\n");
printf("                          note that this can be the same or different than\n");
printf("        the overlay data.\n\n");
printf("   fthresh min_overlay_threshold : use specified value as threshold minimum.\n\n");
printf("   fmin mid_overlay_threshold : use specified value as threshold midpoint.\n\n");
printf("   fslope overlay_threshold_slope : use specified value as threshold slope.\n\n");
printf("   revphaseflag 1|0 : if 1, show overlay data with reversed sign.\n\n");
printf("   truncphaseflag 1|0 : if 1, don't show negative overlay data.\n\n");
printf("   overlaycache 1|0 : if 1, build cache for overlay data. changing slices will\n");
printf("          be faster, but there will be a 1-4 minute pause when\n");
printf("                      caching data.\n\n");
printf("examples:\n\n");
printf("to view a subject in $SUBJECTS_DIR:\n\n");
printf("   tkmedit anders T1\n\n");
printf("to view a volume in another directory:\n\n");
printf("   tkmedit -f /the/path/to/volume\n\n");
printf("to view a subject with surface overlayed:\n\n");
printf("   tkmedit anders T1 lh.white\n\n");
printf("to view a subject with parcellation overlay:\n\n");
printf("   tkmedit anders T1 -parcellation /path/to/volume /the/color/file.dat\n\n");
printf("to view a subject with functional data overlayed:\n\n");
printf("   tkmedit anders T1 -overlay /path/to/bfile/data/stem\n\n");
printf("to view a subject with functional data overlayed and a time course graph:\n\n");
printf("   tkmedit anders T1 -overlay /path/to/bfile/overlay/data/stem -timecourse /path/to/bfile/timecourse/data/stem\n\n");


      exit(1);
    }

    // set the zoom factor for historical reasons
    zf = ozf = 2;
    fsf = (float)zf;
    surflinewidth = zf;

    
    /* start parsing args */
    nCurrentArg = 1;
    bFatalError = FALSE;
    while( nCurrentArg < argc
     && FALSE == bFatalError) {

      /* get the arg */
      strcpy( sArg, argv[nCurrentArg] );
      
      /* check for a switch */
      if( '-' == sArg[0] ) {

  if( MATCH( sArg, "-tcl" ) ) {

    /* check for one more. */
    if( argc > nCurrentArg + 1 
        && '-' != argv[nCurrentArg][0] ) {
      
      /* copy the script name in */
      strcpy( gsBatchScriptName, argv[nCurrentArg+1] );
      gbRunBatchScript = TRUE;
      nCurrentArg += 2;
  
    } else {

      /* misuse of that switch */
      OutputPrint "-tcl switch needs an argument, the name of the script to run\n"
        EndOutputPrint;
      nCurrentArg ++;
    }

  } else if( MATCH( sArg, "-o" ) ) {

    /* make sure there are enough args */
    if( argc > nCurrentArg + 2 ) {
      
      /* read the overlay path and stem */
      sprintf( sOverlayPathAndStem, "%s/%s",
         argv[nCurrentArg+1], argv[nCurrentArg+2] );
      bLoadingOverlay = TRUE;
      nCurrentArg += 3;

    } else {

      /* misuse of that switch */
      OutputPrint "-o switch needs two arguments, path and stem\n"
        EndOutputPrint;
      nCurrentArg ++;
    }

    /* check for two more. */
    if( argc > nCurrentArg + 1 
        && '-' != argv[nCurrentArg][0] ) {

      /* read in time course path and stem. */
      sprintf( sTimeCoursePathAndStem,  "%s/%s",
         argv[nCurrentArg], argv[nCurrentArg+1] );
      bLoadingTimeCourse = TRUE;
      nCurrentArg += 2;
        }

  } else if( MATCH( sArg, "-overlay" ) ) {

    /* make sure there are enough args */
    if( argc > nCurrentArg + 1 ) {
      
      /* copy arg into a destructible string */
      strcpy( sOverlayPathAndStem, argv[nCurrentArg+1] );
      bLoadingOverlay = TRUE;
      nCurrentArg += 2;

    } else {

      /* misuse of that switch */
      OutputPrint "-overlay switch needs argument, path with stem\n"
        EndOutputPrint;
      nCurrentArg ++;
    }

  } else if( MATCH( sArg, "-timecourse" ) ) {

    /* make sure there are enough args */
    if( argc > nCurrentArg + 1 ) {
      
      /* copy arg into a destructible string */
      strcpy( sTimeCoursePathAndStem, argv[nCurrentArg+1] );
      bLoadingTimeCourse = TRUE;
      nCurrentArg += 2;

    } else {

      /* misuse of that switch */
      OutputPrint "-timecourse switch needs argument, path with stem\n"
        EndOutputPrint;
      nCurrentArg ++;
    }

  } else if( MATCH( sArg, "-parcellation" ) ) {

    /* make sure there are enough args */
    if( argc > nCurrentArg + 2 ) {
      
      /* copy path and color file */
      strcpy( sParcellationPath, argv[nCurrentArg+1] );
      strcpy( sParcellationColorFile, argv[nCurrentArg+2] );
      bLoadingParcellation = TRUE;
      nCurrentArg += 3;

    } else {

      /* misuse of that switch */
      OutputPrint "-parcellation switch needs two arguments, path of data and color file to use\n"
        EndOutputPrint;
      nCurrentArg ++;
    }

  } else if( MATCH( sArg, "-f" ) ) {

    /* make sure subject is not already declared */
    if( bSubjectDeclared ) {
      OutputPrint "Used -f switch when subject was already declared\n"
        EndOutputPrint;
      bFatalError = TRUE;

    /* check for path */
    } else if( argc > nCurrentArg + 1 ) {

      /* read the path */
      strcpy( sSubject, argv[nCurrentArg+1] );
      bUsingMRIRead = TRUE;
      bSubjectDeclared = TRUE;
      nCurrentArg += 2;

      /* save subject home */
      SetSubjectHomeDir( sSubject );

    } else {

      /* misuse of that switch */
      OutputPrint "-f switch needs one argument, path to volume\n"
        EndOutputPrint;
      nCurrentArg ++;
    }

  } else if( MATCH( sArg, "-fthresh" ) ) {

    /* check for the value following the switch */
    if( argc > nCurrentArg + 1 ) {

      /* get the value */
      min = (FunV_tFunctionalValue) atof( argv[nCurrentArg+1] );
      bThresh = TRUE;
      nCurrentArg +=2 ;

    } else { 
      
      /* misuse of that switch */
      OutputPrint "-fthresh requires one argument, the threshold as a float\n" EndOutputPrint;
      nCurrentArg += 1;
    }

  } else if( MATCH( sArg, "-fmid" ) ) {

    /* check for the value following the switch */
    if( argc > nCurrentArg + 1 ) {

      /* get the value */
      mid = (FunV_tFunctionalValue) atof( argv[nCurrentArg+1] );
      bMid = TRUE;
      nCurrentArg +=2 ;

    } else { 
      
      /* misuse of that switch */
      OutputPrint "-fmid requires one argument, the midpoint as a float\n" EndOutputPrint;
      nCurrentArg += 1;
    }

  } else if( MATCH( sArg, "-fslope" ) ) {

    /* check for the value following the switch */
    if( argc > nCurrentArg + 1 ) {

      /* get the value */
      slope = (FunV_tFunctionalValue) atof( argv[nCurrentArg+1] );
      bSlope = TRUE;
      nCurrentArg +=2 ;

    } else { 
      
      /* misuse of that switch */
      OutputPrint "-fslope requires one argument, the slope as a float\n" EndOutputPrint;
      nCurrentArg += 1;
    }

  } else if( MATCH( sArg, "-revphaseflag" ) ) {

    /* check for the value following the switch */
    if( argc > nCurrentArg + 1
        && '-' != argv[nCurrentArg+1][0] ) {

      /* get the value */
      nRevPhaseFlag = atoi( argv[nCurrentArg+1] );
      bRevPhaseFlag = TRUE;
      nCurrentArg +=2 ;

    } else { 
      
      /* misuse of that switch */
      OutputPrint "-revphaseflag requires one argument, the revphaseflag value as 0 or 1\n" EndOutputPrint;
      nCurrentArg += 1;
    }

  } else if( MATCH( sArg, "-truncphaseflag" ) ) {

    /* check for the value following the switch */
    if( argc > nCurrentArg + 1
        && '-' != argv[nCurrentArg+1][0] ) {

      /* get the value */
      nTruncPhaseFlag = atoi( argv[nCurrentArg+1] );
      bTruncPhaseFlag = TRUE;
      nCurrentArg +=2 ;

    } else { 
      
      /* misuse of that switch */
      OutputPrint "-truncphaseflag requires one argument, the truncphaseflag value as 0 or 1\n" EndOutputPrint;
      nCurrentArg += 1;
    }

  } else if( MATCH( sArg, "-aux" ) ) {

    /* check for the value following the switch */
    if( argc > nCurrentArg + 1
        && '-' != argv[nCurrentArg+1][0] ) {
      
      /* read in the aux file name */
      strcpy( sAuxVolume, argv[nCurrentArg+1] );
      bLoadingAuxVolume = TRUE;
      nCurrentArg += 2;

    } else { 
      
      /* misuse of that switch */
      OutputPrint "-aux requires one argument, the volume to load\n" EndOutputPrint;
      nCurrentArg += 1;
    }

  } else if( MATCH( sArg, "-headpts" ) ) {

    /* check for the value following the switch */
    if( argc > nCurrentArg + 1
        && '-' != argv[nCurrentArg+1][0] ) {
      
      /* read in the head pts and transform file name */
      strcpy( sHeadPts, argv[nCurrentArg+1] );
      
      /* if they gave us a transform file as well... */
      if( argc > nCurrentArg + 2
    && '-' != argv[nCurrentArg+2][0] ) {
        
        /* save that */
        strcpy( sHeadPtsTransform, argv[nCurrentArg+2] );
        bHaveHeadPtsTransform = TRUE;
      } else {
        bHaveHeadPtsTransform = FALSE;
      }

      bLoadingHeadPts = TRUE;
      nCurrentArg += 3;

    } else { 
      
      /* misuse of that switch */
      OutputPrint "-headpts requires two arguments, the head points file and the transform file\n" EndOutputPrint;
      nCurrentArg += 1;
    }

  } else if( MATCH( sArg, "-overlaycache" ) ) {

    /* check for the value following the switch */
    if( argc > nCurrentArg + 1
        && '-' != argv[nCurrentArg+1][0] ) {

      /* get the value */
      nUseOverlayCacheFlag = atoi( argv[nCurrentArg+1] );
      bUseOverlayCacheFlag = TRUE;
      nCurrentArg +=2 ;

    } else { 
      
      /* misuse of that switch */
      OutputPrint "-overlaycache requires one argument, the value as 0 or 1\n" EndOutputPrint;
      nCurrentArg += 1;
    }

  } else if( MATCH( sArg, "-interface" ) ) {
    
    /* check for another value */
    if( argc > nCurrentArg + 1
        && '-' != argv[nCurrentArg+1][0] ) {
      
      /* copy in the interface name */
      strcpy( gInterfaceScriptName, argv[nCurrentArg+1] );
      nCurrentArg += 2;
      
    } else { 
      
      /* misuse of that switch */
      OutputPrint "-interface requires one argument, the filename of the script to use\n" EndOutputPrint;
      nCurrentArg += 1;
    }
    
  } else if( MATCH( sArg, "-" ) ) {

    /* set no edit mode. */
    bNoEdit = TRUE;
    nCurrentArg += 1;

  } else if( MATCH( sArg, "-csurf" ) ) {

    /* set name of interface to tkmedit_csurf.tcl */
    gbUseCsurfInterface = TRUE;
    nCurrentArg ++;

    
  } else {
    
    OutputPrint "Unrecognized switch %s\n",
      sArg EndOutputPrint;
    nCurrentArg ++;
    
  }

  /* check for local keyword */
      } else if ( MATCH( sArg, "local" ) || MATCH( sArg, "." ) ) {
  
    /* make sure subject is not already declared */
    if( bSubjectDeclared ) {
      OutputPrint "Used local mode switch when subject was already declared\n" EndOutputPrint;
      bFatalError = TRUE;

      /* save subject home */
      strcpy( gsSubjectHomeDir, gsUserHomeDir );

    } else {

      /* set local flag */
      bLocalImageDir = TRUE;
      bSubjectDeclared = TRUE;
      nCurrentArg ++;
    }
    
      } else {
  
  /* this is the name/imagedir part. only do it if we haven't
     declared the sunject yet. */

  if( !bSubjectDeclared ) {

    /*make sure we have enough args and they arn't switches */
    if( argc > nCurrentArg+1 ) {
      
      /* make sure the next two args arn't switches or local/. args */
      if( '-' == argv[nCurrentArg+1][0] ) {
        
        /* image dir is missing. */
        OutputPrint "Image directory is missing.\n" EndOutputPrint;
        bFatalError = TRUE;
        
      } else if (  MATCH( argv[nCurrentArg+1], "local" )
      || MATCH( argv[nCurrentArg+1], "." ) ) {

        /* set image dir to a plain dot and use mriread. */
        strcpy( sSubject, argv[nCurrentArg+1] );
        bUsingMRIRead = TRUE;
        bSubjectDeclared = TRUE;
        nCurrentArg += 1;

        /* save subject home */
        SetSubjectHomeDir( sSubject );

      } else {
        
        /* read in subject and image name. */
        strcpy( sSubject, argv[nCurrentArg] );
        strcpy( sImageDir, argv[nCurrentArg+1] );
        bSubjectDeclared = TRUE;
        nCurrentArg += 2;

        /* save subject home */
        if( !SetSubjectHomeDirFromEnv( sSubject ) ) {

    OutputPrint "If using subject / image dir arguments, SUBJECTS_DIR must be set.\n" EndOutputPrint;
    bFatalError = TRUE;
        }
      }
      
      /* check for a surface. if we have enough args... */
      if( argc > nCurrentArg ) {
        
        /* and this one isn't a switch or the local flag... */
        if( '-' != argv[nCurrentArg][0] 
      && !MATCH( argv[nCurrentArg], "local" )
      && !MATCH( argv[nCurrentArg], "." ) ) {
    
    /* read it as a surface. */
    strcpy( sSurface, argv[nCurrentArg] );
    bSurfaceDeclared = TRUE;
    nCurrentArg ++;
        }
      }
    } else {
      OutputPrint "Requires subject name and image type.\nType \'tkmedit\' with no arguments for more details.\n" EndOutputPrint;
      bFatalError = TRUE;
    }

  } else {
    
    /* totally unrecognized */
    OutputPrint "Unrecognized argument %s\n",
      sArg EndOutputPrint;
    nCurrentArg++;
  }
      }
    }

    /* check for fatal error. */
    if( bFatalError ) {
      exit(1);
    }

    /* if using local directory, copy 'local' in subject name for
       historical reasons. */
    if( bLocalImageDir ) {
      strcpy( sSubject, "local" );
      strcpy( sImageDir, "" );
    }


    /* disable editing */
    if( bNoEdit ) {
      editflag = FALSE;
    }

    /* if we loaded a surface */
    if( bSurfaceDeclared ) {
      
      /* set the stupid flag */
      surfflag = TRUE;

    } else {

      /* else make default surface name */
      strcpy( sSurface, "rh.orig" );
    }

    if ( bUsingMRIRead ) {
      /* attempt to read main volume */
      if( read_images ( sSubject ) ) {
  OutputPrint "ERROR: Couldn't load volume %s.\n",
    sSubject EndOutputPrint;
  exit (1);
      }
    } else {
      /* attempt to read main volume */
      if( read_images ( sImageDir ) ) {
  OutputPrint "ERROR: Couldn't load volume %s.\n",
    sImageDir EndOutputPrint;
  exit (1);
      }
    }

    /* now build max intensity projection */
    BuildVolumeMaxIntProj( gAnatomicalVolume, &gAnatomicalMaxIntProj );

    /* if reading in an aux image... */
    if( bLoadingAuxVolume ) {
      read_second_images( sAuxVolume );
      BuildVolumeMaxIntProj( gAuxAnatomicalVolume, &gAuxAnatomicalMaxIntProj );
    }

    /* create our transformation object */
    InitTransformation ();

    /* load surface. trasnsform must be inited first. */
    if ( surfflag ) {
      LoadSurface( sSurface );
    }

    /* load parcellation */
    if( bLoadingParcellation ) {
      LoadParcellationVolume( sParcellationPath, sParcellationColorFile );
    }

    /* load head pts */
    if( bLoadingHeadPts ) {
      if( bHaveHeadPtsTransform ) 
  LoadHeadPts( sHeadPts, sHeadPtsTransform ); 
      else
  LoadHeadPts( sHeadPts, NULL );
    }

    /* load functional data */
    if( bLoadingOverlay ) {
      LoadFunctionalOverlay( sOverlayPathAndStem );
    }

    if( bLoadingTimeCourse ) {
      LoadFunctionalTimeCourse( sTimeCoursePathAndStem );
    }

    /* set functional color scale stuff */
    if( bThresh || bMid || bThresh ) {
      eFunctional = FunV_SetThreshold( gFunctionalVolume, min, mid, slope );
    }
    if( bTruncPhaseFlag ) {
      eFunctional = FunV_SetDisplayFlag( gFunctionalVolume, 
           FunV_tDisplayFlag_Ol_TruncateNegative,
           (tBoolean) nTruncPhaseFlag );
    }
    if( bRevPhaseFlag ) {
      eFunctional = FunV_SetDisplayFlag( gFunctionalVolume, 
           FunV_tDisplayFlag_Ol_ReversePhase,
           (tBoolean) nRevPhaseFlag );
    }
    if( bUseOverlayCacheFlag ) {
      eFunctional = FunV_UseOverlayCache( gFunctionalVolume,
            (tBoolean) nUseOverlayCacheFlag );
    }

    return(0);
}

void save_rgb( char* fname ) {

  GLint swapbytes, lsbfirst, rowlength;
  GLint skiprows, skippixels, alignment;
  RGB_IMAGE *image;
  int y,xloc,yloc,width,height,size;
  unsigned short *r,*g,*b;
  unsigned short  *red, *green, *blue;
  FILE *fp;

  if( NULL == gMeditWindow )
    return;

  MWin_GetWindowSize( gMeditWindow,
          &xloc, &yloc, &width, &height );

  size = width*height;

  red = (unsigned short *)calloc(size, sizeof(unsigned short));
  green = (unsigned short *)calloc(size, sizeof(unsigned short));
  blue = (unsigned short *)calloc(size, sizeof(unsigned short));

  glGetIntegerv(GL_UNPACK_SWAP_BYTES, &swapbytes);
  glGetIntegerv(GL_UNPACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, &rowlength);
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, &skiprows);
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv(GL_UNPACK_ALIGNMENT, &alignment);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  glReadPixels(0, 0, width, height, GL_RED,  
         GL_UNSIGNED_SHORT, (GLvoid *)red);
  glReadPixels(0, 0, width, height, GL_GREEN,
         GL_UNSIGNED_SHORT, (GLvoid *)green);
  glReadPixels(0, 0, width, height, GL_BLUE, 
         GL_UNSIGNED_SHORT, (GLvoid *)blue);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, swapbytes);
  glPixelStorei(GL_UNPACK_LSB_FIRST, lsbfirst);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, rowlength);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, skiprows);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, skippixels);
  glPixelStorei(GL_UNPACK_ALIGNMENT, alignment);

  fp = fopen(fname,"w");
  if (fp==NULL){printf("medit: ### can't create file %s\n",fname); return;}
  fclose(fp);
  image = iopen(fname,"w",UNCOMPRESSED(1), 3, width, height, 3);
  for(y = 0 ; y < height; y++) {
    r = red + y * width;
    g = green + y * width;
    b = blue + y * width;
    putrow(image, r, y, 0);
    putrow(image, g, y, 1);
    putrow(image, b, y, 2);
  }
  iclose(image);
  free(red); free(green); free(blue);

  printf("medit: file %s written\n",fname);
}


void GotoSurfaceVertex ( Surf_tVertexSet iSurface, int inVertex ) {

  Surf_tErr eSurface = Surf_tErr_NoErr;
  MWin_tErr eWindow  = MWin_tErr_NoErr;
  xVoxel    anaIdx;
            
  /* get the vertex */
  eSurface = Surf_GetNthVertex( gSurface, iSurface, inVertex, &anaIdx );
  if( Surf_tErr_NoErr != eSurface ) 
    goto error;

  /* tell the window to go there. */
  eWindow = MWin_SetCursor ( gMeditWindow, -1, &anaIdx );
  if( MWin_tErr_NoErr != eWindow )
    goto error;

  goto cleanup;
  
 error:

  DebugPrint "Error in GotoSurfaceVertex( %d, %d )\n",
    (int)iSurface, inVertex EndDebugPrint;

 cleanup:
  return;
}

void FindNearestSurfaceVertex ( Surf_tVertexSet iSurface ) {

  Surf_tErr eSurface = Surf_tErr_NoErr;
  MWin_tErr eWindow  = MWin_tErr_NoErr;
  xVoxel    cursor;
  xVoxel    anaIdx;
  char      sResult[256];
  char      sSetName[256];

  /* get the cursor */
  eWindow = MWin_GetCursor ( gMeditWindow, &cursor );
  if( MWin_tErr_NoErr != eWindow )
    goto error;

  /* get the vertex */
  eSurface = Surf_GetClosestVertex( gSurface, iSurface, &cursor, &anaIdx,
            sResult);
  if( Surf_tErr_NoErr != eSurface ) 
    goto error;

  /* print the result string */
  Surf_GetSurfaceSetName( iSurface, sSetName );
  OutputPrint "Nearest %s vertex to %d, %d, %d:\n\t%s\n", 
    sSetName, xVoxl_ExpandInt( &cursor ), sResult EndOutputPrint;

  /* tell the window to go there. */
  eWindow = MWin_SetCursor ( gMeditWindow, -1, &anaIdx );
  if( MWin_tErr_NoErr != eWindow )
    goto error;

  goto cleanup;
  
 error:

  DebugPrint "Error in FindNearestSurfaceVertex( %d )\n",
    (int)iSurface EndDebugPrint;

 cleanup:
  return;
}

void WriteVoxelToControlFile ( xVoxelRef inVolumeVox ) {

  char sFileName[256] = "";
  FILE *fp;
  Real theRASX, theRASY, theRASZ;
 xVoxelRef theVoxel;

  // make a new voxel
  xVoxl_New ( &theVoxel );
  
  /* make the file name */
  MakeFileName( "control.dat", tkm_tFileName_ControlPoints, sFileName );

  /* open it */
  fp = fopen( sFileName, "a+" );
  if ( fp == NULL ) {
    printf("medit: ### can't create file %s\n",sFileName);
    return;
  }

  // convert voxel to ras
  VoxelToRAS ( xVoxl_ExpandInt(inVolumeVox), 
         &theRASX, &theRASY, &theRASZ );
  
  // write RAS space pt to file
  fprintf ( fp,"%f %f %f\n", theRASX, theRASY, theRASZ );
  DebugPrint "writing RAS point to %s...\n", sFileName EndDebugPrint;
  
  // close the file
  fclose( fp );

  // free the voxel
  xVoxl_Delete ( &theVoxel );
}

void WriteVoxelToEditFile( xVoxelRef inVolumeVox ) {

  char sFileName[256];
  FILE *fp;
  Real theTalX, theTalY, theTalZ;
  Real theRASX, theRASY, theRASZ;
 xVoxelRef theVoxel;

  // make a new voxel
  xVoxl_New ( &theVoxel );

   /* make the file name */
  MakeFileName( "edit.dat", tkm_tFileName_Edit, sFileName );

  /* open it */
  fp = fopen( sFileName, "w" );
  if (fp==NULL) {
    DebugPrint "Couldn't create file %s\n", sFileName EndDebugPrint;
    OutputPrint "Couldn't create edit file.\n" EndOutputPrint;
    return;
  }

  // convert to ras
  VoxelToRAS ( xVoxl_ExpandInt(inVolumeVox), 
         &theRASX, &theRASY, &theRASZ );
  
  // write RAS space pt to file
  fprintf ( fp,"%f %f %f\n", theRASX, theRASY, theRASZ );
  DebugPrint "writing RAS point to %s...\n", sFileName EndDebugPrint;
  
  // if we have a tal transform for this volume...
  if ( gbAnaToTalTransformLoaded ) {
    
    // ras to tal
    transform_point ( get_linear_transform_ptr( &gAnaToTalTransform ),
          theRASX, theRASY, theRASZ,
                      &theTalX, &theTalY, &theTalZ );
      
    // write tal space point to file
    fprintf(fp, "%f %f %f\n", theTalX, theTalY, theTalZ );
    DebugPrint "writing Tal point to %s...\n", sFileName EndDebugPrint;
  }
  
  // close the file
  fclose(fp);

  // free the voxel
  xVoxl_Delete ( &theVoxel );
}

void ReadCursorFromEditFile () {

  char theFileName[NAME_LENGTH];
  FILE* theFile;
  float theRASX, theRASY, theRASZ;
  int theVoxX, theVoxY, theVoxZ;
  xVoxelRef theVolumeVox;

  /* make the file name. */
  MakeFileName( "edit.dat", tkm_tFileName_Edit, theFileName );

  /* open it. */
  theFile = fopen ( theFileName, "r" );
  if ( NULL == theFile ) {
    DebugPrint "Couldn't find file %s\n", theFileName EndDebugPrint;
    OutputPrint "Couldn't find edit file.\n" EndOutputPrint;
    return;
  }

  /* read the point */
  fscanf ( theFile, "%f %f %f", &theRASX, &theRASY, &theRASZ );
  fclose ( theFile );

  /* convert to volume voxel. */
  RASToVoxel ( theRASX,  theRASY,  theRASZ,
         &theVoxX, &theVoxY, &theVoxZ ); 
 
  /* build and set cursor */
  xVoxl_New ( &theVolumeVox );
  xVoxl_Set ( theVolumeVox, theVoxX, theVoxY, theVoxZ );
  MWin_SetCursor ( gMeditWindow, -1, theVolumeVox );
  MWin_SetZoomCenterToCursor ( gMeditWindow, -1 );
  xVoxl_Delete ( &theVolumeVox );
}

void UpdateAndRedraw () {

  if ( NULL != gMeditWindow ) {
    MWin_RedrawAll( gMeditWindow );
  }
}

void redraw () {

  /* tell medit window to redraw */
  MWin_Redraw( gMeditWindow );
}


void
mri2pix(float xpt, float ypt, float zpt, int *jpt, int *ipt,int *impt)
{
  if (ptype==0) /* Horizontal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((ypt-yy0)/ps+0.5);
    *impt = (int)((zpt-zz0)/st+0.5);
  } else if (ptype==2) /* Coronal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((255.0-(zz1-zpt)/ps)+0.5);
    *impt = (int)((ypt-yy0)/st+0.5);
  } else if (ptype==1) /* Sagittal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((ypt-yy0)/ps+0.5);
    *impt = (int)((zpt-zz0)/st+0.5);
  }
}

int 
imval(float px,float py,float pz)
{
  float x,y,z;
  int j,i,imn;

  x = px*tm[0][0]+py*tm[0][1]+pz*tm[0][2]+tm[0][3];
  y = px*tm[1][0]+py*tm[1][1]+pz*tm[1][2]+tm[1][3];
  z = px*tm[2][0]+py*tm[2][1]+pz*tm[2][2]+tm[2][3];
  mri2pix(x,y,z,&j,&i,&imn);
  if (imn>=0&&imn<numimg&&i>=0&&i<ynum&&j>=0&&j<xnum) {

    return ( GetVoxelValue ( gAnatomicalVolume, 
           j, ynum-1-i, imn ) );
  }
  else
    return 0;
}

float
Error(int p,float dp)
{
  int i,num;
  float mu,error,sum;
  float mu1,mu2,sum1,sum2;

  if (p>=0)
    par[p] += dp;
  mu = mu1 = mu2 = 0;
  num = 0;
  for (i=0;i<npts;i++)
  {
    mu += imval(ptx[i],pty[i],ptz[i]);
    mu1 += imval(ptx[i]*0.9,pty[i]*0.9,ptz[i]*0.9);
    mu2 += imval(ptx[i]*1.05,pty[i]*1.05,ptz[i]*1.05);
    num ++;
  }
  mu /= num;
  mu1 /= num;
  mu2 /= num;
  sum = sum1 = sum2 = 0;
  num = 0;
  for (i=0;i<npts;i++)
  {
    error = imval(ptx[i],pty[i],ptz[i])-mu;
    sum += error*error;
    error = imval(ptx[i]*0.9,pty[i]*0.9,ptz[i]*0.9)-mu1;
    sum1 += error*error;
    error = imval(ptx[i]*1.05,pty[i]*1.05,ptz[i]*1.05)-mu2;
    sum2 += error*error;
    num ++;
  }
  sum = sqrt((sum+sum2)/num);
  if (p>=0)
    par[p] -= dp;
  return sum;
}

void
optimize(int maxiter)
{
  float lambda = 0.03;
  float epsilon = 0.1;
  float momentum = 0.8;
  int iter,p;
  float dE[3];
  float error;

  for (iter=0;iter<maxiter;iter++)
  {
    error = Error(-1,0);
    printf("%d: %5.2f %5.2f %5.2f %7.3f\n",
           iter,par[0],par[1],par[2],error);
    for (p=0;p<3;p++)
    {
      dE[p] = tanh((Error(p,epsilon/2)-Error(p,-epsilon/2))/epsilon);
    }
    for (p=0;p<3;p++)
    {
      par[p] += (dpar[p] = momentum*dpar[p] - lambda*dE[p]);
    }
  }
  error = Error(-1,0);
  printf("%d: %5.2f %5.2f %5.2f %7.3f\n",
         iter,par[0],par[1],par[2],error);
}
void
optimize2(void)
{
  float epsilon = 0.5;
  float p0,p1,p2,p0min,p1min,p2min;
  float error,minerror;

  p0min = p1min = p2min = 0.0f;
  error = Error(-1,0);
  minerror = error;
  printf("%5.2f %5.2f %5.2f %7.3f\n",
         par[0],par[1],par[2],error);
  for (p0 = -10;p0 <= 10;p0 += epsilon)
  for (p1 = -10;p1 <= 10;p1 += epsilon)
  for (p2 = -10;p2 <= 10;p2 += epsilon)
  {
    par[0] = p0;
    par[1] = p1;
    par[2] = p2;
    error = Error(-1,0);
    if (error<minerror)
    {
      printf("%5.2f %5.2f %5.2f %7.3f\n",
             par[0],par[1],par[2],error);
      minerror = error;
      p0min = p0;
      p1min = p1;
      p2min = p2;
    }
  }
  par[0] = p0min;
  par[1] = p1min;
  par[2] = p2min;
  error = Error(-1,0);
  printf("%5.2f %5.2f %5.2f %7.3f\n",
         par[0],par[1],par[2],error);
}

int
read_images( char *isName ) {

  char  sPath[256]         = "";
  char* pEnd               = NULL;
  int   eRead              = 0;
  char  sTclArguments[256] = "";

  /* make a filename */
  MakeFileName( isName, tkm_tFileName_Volume, sPath );

  /* if there is a /COR- at the end, take it off and pass
     the rest to ReadVolumeWithMRIRead */
  pEnd = strstr( sPath, "/COR-" );
  if( NULL != pEnd )
    *pEnd = '\0';

  /* save the dir name as the volume src */
  strcpy( gVolumeSrc, sPath );

  /* attempt to read the volume */
  eRead = ReadVolumeWithMRIRead( &gAnatomicalVolume, 
         &gAnaToTalTransform,
         sPath );
 
  /* if no good, set volume to null and return */
  if( 0 != eRead ) {
    gAnatomicalVolume         = NULL;
    gbAnaToTalTransformLoaded = FALSE;
    return 1;
  }

  /* save the volume size */
  gVolumeDimension = knVolumeSize;

  /* check to see if the transform is loaded */
  if( get_linear_transform_ptr( &gAnaToTalTransform ) ) {
    gbAnaToTalTransformLoaded = TRUE;
  } else {
    gbAnaToTalTransformLoaded = FALSE;

    /* put up an error dlog */
    sprintf( sTclArguments, "\"Couldn't find Talairach transform for main volume. Local RAS coordinates will be displayed instead.\"" );
    tkm_SendTclCommand( tkm_tTclCommand_AlertDlog, sTclArguments ); 
  }

  /* get the subject and volume name */
  ExtractSubjectName( gVolumeSrc, gSubjectName );
  ExtractVolumeName( gVolumeSrc, gSubjectVolumeName );

  /* set data in window */
  if( NULL != gMeditWindow ) {
    MWin_SetVolume( gMeditWindow, -1, 
        gAnatomicalVolume, gVolumeDimension );
  }

  return 0;
}


int
read_second_images( char *isName )
{

  char sPath[256];
  char* pEnd;
  int  eRead;

  /* make a filename */
  MakeFileName( isName, tkm_tFileName_Volume, sPath );

  //  DebugPrint "argument: %s\n after MakeFileName: %s\n",
  //    isName, sPath EndDebugPrint;

  /* if there is a /COR- at the end, take it off and pass
     the rest to ReadVolumeWithMRIRead */
  pEnd = strstr( sPath, "/COR-" );
  if( NULL != pEnd )
    *pEnd = '\0';

  //  DebugPrint "after stripping for /COR-: %s\n", sPath EndDebugPrint;

  /* save the dir name as the volume src */
  strcpy( gAuxVolumeSrc, sPath );

  /* attempt to read the volume */
  eRead = ReadVolumeWithMRIRead( &gAuxAnatomicalVolume, 
         &gAuxAnaToTalTransform,
         sPath );

  /* if no good, set volume to null and return */
  if( 0 != eRead ) {
    gAuxAnatomicalVolume         = NULL;
    gbAuxAnaToTalTransformLoaded = FALSE;
    return 1;
  }

  /* check to see if transform is loaded */
  if( get_linear_transform_ptr( &gAuxAnaToTalTransform ) ) {
    gbAuxAnaToTalTransformLoaded = TRUE;
  } else {
    gbAuxAnaToTalTransformLoaded = FALSE;
  }

  /* get the subject name */
  ExtractSubjectName( gAuxVolumeSrc, gAuxSubjectName );
  ExtractVolumeName( gAuxVolumeSrc, gAuxSubjectVolumeName );


  /* set data in window */
  if( NULL != gMeditWindow ) {
    MWin_SetAuxVolume( gMeditWindow, -1, 
        gAuxAnatomicalVolume, gVolumeDimension );
  }

  return 0;
}

int
write_images(char *fpref)
{
  int k;
  FILE *fptr;
  char fname[STRLEN];
  tVolumeValue * theSlicePtr;
  tVolumeValue* buf;
  int bufsize;
  tBoolean bSaveAll = FALSE;

  /* if the name is not the directory we loaded from, save all slices. else
     just save the ones we edited. */
  if( !MATCH( fpref, gVolumeSrc )) {
    bSaveAll = TRUE;
  }

  bufsize = ((unsigned long)xnum)*ynum;
  buf = (tVolumeValue *)lcalloc(bufsize,sizeof(char));

  if ( !editflag && !bSaveAll ) { 
    OutputPrint "Changes cannot be made because you launched tkmedit with editing disabled. To make changes, reboot tkmedit without the no edit flag.\n" EndOutputPrint;
    return(0);
  }

  /* for all slices...*/
  for ( k = 0; k < numimg; k++ ) {

    /* if this slice is dirty or we're saving in a different direcotry */
    if ( changed[ k ]
   || bSaveAll ) {

      if( !bSaveAll ) {
  changed[k] = FALSE;
      }

      /* make file name and open the file. */
      MakeFileName( fpref, tkm_tFileName_Volume, fname );
      sprintf( fname, "%s/COR-%03d", fname, k+imnr0 );
      fptr = fopen(fname,"wb");

      /* if we got it... */
      if( NULL != fptr ) {
  
  /* write the data and close the file. */
  theSlicePtr = GetVolumeSlicePtr ( gAnatomicalVolume, k );
  memcpy ( buf, theSlicePtr, xnum*ynum );
  fwrite ( buf,sizeof(char),bufsize,fptr);
  fclose(fptr);

  DebugPrint "Wrote file %s\n", fname EndDebugPrint;

      } else {

  OutputPrint "ERROR: Couldn't open file %s\n",
    fname EndOutputPrint;
      }
    }
  }
  editedimage = FALSE;
  return(0);
}

//static int nIdleCount = 0;

void tkm_HandleIdle () {


  //  if ( ++nIdleCount == 2000 ) {

    // just call the tk event handling function
    Tk_DoOneEvent ( TK_ALL_EVENTS | TK_DONT_WAIT );

    //    nIdleCount = 0;
    //  }
}

// ================================================================== SURFACES

void LoadSurface ( char* isName ) {

  Surf_tErr eSurface   = Surf_tErr_NoErr;
  tBoolean  bLoaded    = FALSE;
  char      sName[256] = "";

  /* make file name */
  MakeFileName( isName, tkm_tFileName_Surface, sName );

  /* create the surface */
  eSurface = Surf_New( &gSurface, sName, gRASTransform );
  if( Surf_tErr_NoErr != eSurface ) {
    DebugPrint "Surf error %d in LoadSurface: %s\n",
      eSurface, Surf_GetErrorString( eSurface ) EndDebugPrint;
  }

  /* see if it was loaded */
  eSurface = Surf_IsVertexSetLoaded( gSurface, Surf_tVertexSet_Main, 
             &bLoaded );
  if( Surf_tErr_NoErr != eSurface ) {
    DebugPrint "Surf error %d in LoadSurface: %s\n",
      eSurface, Surf_GetErrorString( eSurface ) EndDebugPrint;
  }

  /* set the medit window surface. */
  MWin_SetSurface( gMeditWindow, -1, gSurface );

  /* enable our loading options */
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceLoadingOptions, 
          bLoaded ? "1" : "0" );

  /* and viewing options */
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceViewingOptions, 
          bLoaded ? "1" : "0" );

  /* turn surface display on or off */
  MWin_SetDisplayFlag( gMeditWindow, -1, DspA_tDisplayFlag_MainSurface, 
           bLoaded ? TRUE : FALSE);

  /* this is a default setting */
  MWin_SetDisplayFlag( gMeditWindow, -1,
           DspA_tDisplayFlag_InterpolateSurfaceVertices,
           bLoaded ? TRUE : FALSE);

  /* load other vertices */
  LoadSurfaceVertexSet( Surf_tVertexSet_Original, "orig" );
  LoadSurfaceVertexSet( Surf_tVertexSet_Pial, "pial" );
}

void LoadSurfaceVertexSet ( Surf_tVertexSet iSet,
          char* fname ) {

  Surf_tErr eSurface = Surf_tErr_NoErr;
  tBoolean  bLoaded  = FALSE;
  DspA_tDisplayFlag flag = DspA_tDisplayFlag_None;
  tkm_tTclCommand command = 0;

  if( !gSurface )
    return;

  eSurface = Surf_LoadVertexSet( gSurface, fname, iSet );
  if( Surf_tErr_NoErr != eSurface ) {
    DebugPrint "Surf error %d in LoadSurfaceVertexSet: %s\n",
      eSurface, Surf_GetErrorString( eSurface ) EndDebugPrint;
  }

  eSurface = Surf_IsVertexSetLoaded( gSurface, iSet, &bLoaded );
  if( Surf_tErr_NoErr != eSurface ) {
    DebugPrint "Surf error %d in LoadSurfaceVertexSet: %s\n",
      eSurface, Surf_GetErrorString( eSurface ) EndDebugPrint;
  }

  /* get command and flag to set */
  if( iSet == Surf_tVertexSet_Pial ) {
    flag = DspA_tDisplayFlag_CanonicalSurface;
    command = tkm_tTclCommand_ShowCanonicalSurfaceViewingOptions;
  } else if( iSet == Surf_tVertexSet_Original ) {
    flag = DspA_tDisplayFlag_OriginalSurface;
    command = tkm_tTclCommand_ShowOriginalSurfaceViewingOptions;
  }

  /* turn flag on or off and enable or disable viewing optiosn */
  MWin_SetDisplayFlag( gMeditWindow, -1, flag, bLoaded );
  tkm_SendTclCommand( command, bLoaded?"1":"0" );

  /* set the surface in the window to purge the cache */
  MWin_SetSurface( gMeditWindow, -1, gSurface );
}

void UnloadSurface () {

  Surf_tErr eSurface = Surf_tErr_NoErr;

  if( !gSurface )
    return;

  /* free the surface. */
  eSurface = Surf_Delete( &gSurface );
  if( Surf_tErr_NoErr != eSurface ) {
    DebugPrint "Surf error %d in UnloadSurface: %s\n",
      eSurface, Surf_GetErrorString( eSurface ) EndDebugPrint;
  }

  /* disable our loading options */
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceLoadingOptions, "0" );

/* update the medit window. */
  MWin_SetSurface( gMeditWindow, -1, gSurface );
}

// ===========================================================================



void
write_dipoles(char *fname)
{
  int i,j,k,di,dj,dk;
  float x,y,z;
  FILE *fptr;
  int ***neighindex,index,nneighbors;

  ndip = 0;
/*
  ndip=2;
*/
  neighindex=(int ***)malloc((numimg/dip_spacing+1)*sizeof(int **));
  for (k=0;k<numimg;k+=dip_spacing)
  { 
    neighindex[k/dip_spacing]=(int **)malloc((ynum/dip_spacing+1)*sizeof(int *));
    for (i=0;i<ynum;i+=dip_spacing)
      neighindex[k/dip_spacing][i/dip_spacing]=(int *)malloc((xnum/dip_spacing+1)*
                                                             sizeof(int));
  }

  for (k=0;k<numimg;k+=dip_spacing)
  for (i=0;i<ynum;i+=dip_spacing)
  for (j=0;j<xnum;j+=dip_spacing)

    if ( GetVoxelValue ( gAnatomicalVolume, j, i, k ) != 0 )
      //  if (im[k][i][j]!=0)

/*
  if (sqrt(0.0+SQR(k-numimg/2)+SQR(i-ynum/2)+SQR(j-xnum/2))<=100)
*/
  {
    neighindex[k/dip_spacing][i/dip_spacing][j/dip_spacing] = ndip;
    ndip++;
  } else
   neighindex[k/dip_spacing][i/dip_spacing][j/dip_spacing] = -1;

  fptr = fopen(fname,"w");
  if (fptr==NULL) {printf("medit: ### can't create file %s\n",fname); return;}
  fprintf(fptr,"#!ascii\n");
  fprintf(fptr,"%d\n",ndip);
/*
  fprintf(fptr,"%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
          30.0,-75.0,40.0,0.0,0.0,0.0);
  fprintf(fptr,"%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
          45.0,10.0,70.0,0.0,0.0,0.0);
*/
  for (k=0;k<numimg;k+=dip_spacing)
  for (i=0;i<ynum;i+=dip_spacing)
  for (j=0;j<xnum;j+=dip_spacing)
  if (neighindex[(k)/dip_spacing][(i)/dip_spacing][(j)/dip_spacing]>=0)
  {
    x = xx1-ps*j;
    y = yy0+st*k;
    z = zz1-ps*i;
    fprintf(fptr,"%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",x,y,z,0.0,0.0,0.0);
  }
/*
  fprintf(fptr,"%d \n",0);
  fprintf(fptr,"%d \n",0);
*/
  for (k=0;k<numimg;k+=dip_spacing)
  for (i=0;i<ynum;i+=dip_spacing)
  for (j=0;j<xnum;j+=dip_spacing)
  if (neighindex[(k)/dip_spacing][(i)/dip_spacing][(j)/dip_spacing]>=0)
  {
    nneighbors = 0;
    for (dk= -dip_spacing;dk<=dip_spacing;dk+=2*dip_spacing)
      if (dk>0)
      if (neighindex[(k+dk)/dip_spacing][(i)/dip_spacing][(j)/dip_spacing]>=0)
        nneighbors++;
    for (di= -dip_spacing;di<=dip_spacing;di+=2*dip_spacing)
      if (di>0)
      if (neighindex[(k)/dip_spacing][(i+di)/dip_spacing][(j)/dip_spacing]>=0)
        nneighbors++;
    for (dj= -dip_spacing;dj<=dip_spacing;dj+=2*dip_spacing)
      if (dj>0)
      if (neighindex[(k)/dip_spacing][(i)/dip_spacing][(j+dj)/dip_spacing]>=0)
        nneighbors++;
    fprintf(fptr,"%d ",nneighbors);
    for (dk= -dip_spacing;dk<=dip_spacing;dk+=2*dip_spacing)
      if (dk>0)
      if ((index=neighindex[(k+dk)/dip_spacing][(i)/dip_spacing][(j)/dip_spacing])>=0)
        fprintf(fptr,"%d ",index);
    for (di= -dip_spacing;di<=dip_spacing;di+=2*dip_spacing)
      if (di>0)
      if ((index=neighindex[(k)/dip_spacing][(i+di)/dip_spacing][(j)/dip_spacing])>=0)
        fprintf(fptr,"%d ",index);
    for (dj= -dip_spacing;dj<=dip_spacing;dj+=2*dip_spacing)
      if (dj>0)
      if ((index=neighindex[(k)/dip_spacing][(i)/dip_spacing][(j+dj)/dip_spacing])>=0)
        fprintf(fptr,"%d ",index);
    fprintf(fptr,"\n");
  }
  fclose(fptr);
  printf("medit: dipole file %s written\n",fname);
}

void
write_decimation(char *fname)
{
  int k;
  FILE *fptr;

  fptr = fopen(fname,"w");
  if (fptr==NULL) {printf("medit: ### can't create file %s\n",fname); return;}
  fprintf(fptr,"#!ascii\n");
  fprintf(fptr,"%d\n",ndip);
  for (k=0;k<ndip;k++)
  {
    fprintf(fptr,"1\n");
  }
  fclose(fptr);
  printf("medit: decimation file %s written\n",fname);
}

void
smooth_3d(int niter)
{
  int iter,i,j,k;
  int i2,j2,k2;
  int n;
  float sum;
  tVolumeValue theIntensity;

  if (!dummy_im_allocated) {
    for (k=0;k<numimg;k++) {
      dummy_im[k] = (tVolumeValue **)lcalloc(ynum,sizeof(char *));
      for (i=0;i<ynum;i++) {
        dummy_im[k][i] = (tVolumeValue *)lcalloc(xnum,sizeof(char));
      }
    }
    dummy_im_allocated = TRUE;
  }

  for (iter=0;iter<niter;iter++) {  /* smooth */
    for (j=1;j<xnum-1;j++) {
      printf(".");
      for (k=1;k<numimg-1;k++)
      for (i=1;i<ynum-1;i++) {
        sum = 0;
        n = 0;
        for(k2=k-1;k2<k+2;k2++)  /* 27 */
        for(i2=i-1;i2<i+2;i2++)
        for(j2=j-1;j2<j+2;j2++) {
    sum += GetVoxelValue ( gAnatomicalVolume, j2, i2, k2 );
    //          sum += im[k2][i2][j2];
          n++;
        }
        /*dummy_im[k][i][j] = (sum+im[k][i][j])/(float)(n+1);*/
        dummy_im[k][i][j] = (sum + 27*GetVoxelValue(gAnatomicalVolume,j,i,k))/
    (float)(n+27);
      }
      fflush(stdout);
    }
  }

  /* update view */
  for (k=0;k<numimg;k++)
    for (i=0;i<ynum;i++)
      for (j=0;j<xnum;j++)
  SetVoxelValue ( gAnatomicalVolume, j, i, k, dummy_im[k][i][j] );

  for (k=0;k<numimg;k++)
    changed[k] = TRUE;
  editedimage = TRUE;

  for (k=0;k<6;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
    sim[k][i][j] = 0;

  for (k=0;k<numimg;k++)
    for (i=0;i<ynum;i++)
      for (j=0;j<xnum;j++) {

  theIntensity = GetVoxelValue ( gAnatomicalVolume, j, i, k ) / 2;

  if ( theIntensity > sim[3][i][j]) 
    sim[3][i][j] = theIntensity;
  
  if ( theIntensity > sim[4][k][j]) 
    sim[4][k][j] = theIntensity;
  
  if ( theIntensity > sim[5][i][k]) 
    sim[5][i][k] = theIntensity;
      }

  for (i=0;i<ynum;i++)
    for (j=0;j<xnum;j++)
      for (k=0;k<3;k++)
  sim[k][i][j] = sim[k+3][i][j];

  printf("\nmedit: finished %d smooth_3d steps\n",niter);
}

void
wmfilter_corslice(int imc)
{
  char *getenv(),*mri_dir,fname[NAME_LENGTH];
  int nver,ncor;
  int i,j,k,k2,di,dj,dk,m,n,u,ws2=2,maxi,mini;
  float xcor[MAXCOR],ycor[MAXCOR],zcor[MAXCOR];
  float fwhite_hilim=140,fwhite_lolim=80,fgray_hilim=100;
  float numvox,numnz,numz;
  float f,f2,a,b,c,x,y,z,s;
  float cfracz = 0.6;
  float cfracnz = 0.6;
  double sum2,sum,var,avg,tvar,maxvar,minvar;
  double sum2v[MAXLEN],sumv[MAXLEN],avgv[MAXLEN],varv[MAXLEN],nv[MAXLEN];
  FILE *fptr;
  mri_tOrientation theOrientation;
  int plane;

  // get the plane.
  MWin_GetOrientation ( gMeditWindow, &theOrientation );
  plane = (int)theOrientation;

  k2 = imc/zf;
  if (plane!=CORONAL) {
    printf("medit: ### can only wmfilter CORONAL slice\n"); return; }
  if(k2<ws2 || k2>numimg-ws2) {
    printf("medit: ### slice too close to edge\n"); return; }
  mri_dir = getenv("MRI_DIR");
  if (mri_dir==NULL) {
    printf("medit: ### env var MRI_DIR undefined (setenv, restart)\n");return;}
  if (!flossflag && !spackleflag) { 
    printf("medit: ### no spackle or floss  ...skipping wmfilter\n"); return; }

  sprintf(fname,"%s/lib/bem/%s",mri_dir,DIR_FILE);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {printf("medit: ### File %s not found\n",fname); return;}
  fscanf(fptr,"%d",&nver);
  for (i=0,ncor=0;i<nver;i++) {
    fscanf(fptr,"%*d %f %f %f",&x,&y,&z);
    for (j=0;j<ncor;j++)
      if ((x==-xcor[j]) && (y==-ycor[j]) && (z==-zcor[j])) goto L1;
    xcor[ncor] = x;
    ycor[ncor] = y;
    zcor[ncor] = z;
    ncor++;
  L1:;
  }
  fwhite_hilim = (float)white_hilim;
  fwhite_lolim = (float)white_lolim;
  fgray_hilim = (float)gray_hilim;

  if (!wmfilter_ims_allocated) {
    for (k=0;k<numimg;k++) {
      fill[k] = (tVolumeValue **)lcalloc(ynum,sizeof(char *));
      im_b[k] = (tVolumeValue **)lcalloc(ynum,sizeof(char *));
      for (i=0;i<ynum;i++) {
        fill[k][i] = (tVolumeValue *)lcalloc(xnum,sizeof(char));
        im_b[k][i] = (tVolumeValue *)lcalloc(xnum,sizeof(char));
      }
    }
    wmfilter_ims_allocated = TRUE;
  }

  if (!second_im_allocated)
    alloc_second_im();

  /*for (k=0;k<numimg-1;k++)*/
  for (k=k2-ws2;k<=k2+ws2;k++)
    for (i=0;i<ynum-1;i++)
      for (j=0;j<xnum-1;j++) {
  im_b[k][i][j] = GetVoxelValue ( gAnatomicalVolume, j, i, k );
  SetVoxelValue ( gAuxAnatomicalVolume, j, i, k,
      GetVoxelValue ( gAnatomicalVolume, j, i, k ) );
  if (im_b[k][i][j]>fwhite_hilim || im_b[k][i][j]<fwhite_lolim)
    im_b[k][i][j] = 0;
  fill[k][i][j] = im_b[k][i][j];
      }

  k = k2;
  printf("medit: %d unique orientations\n",ncor);
  printf("medit: white_lolim = %d   white_hilim = %d   gray_hilim = %d\n",
                 white_lolim,white_hilim,gray_hilim);
  printf("medit: wmfiltering coronal slice %d\n",k2);
  if (flossflag)   printf("medit:   floss => '.'\n");
  if (spackleflag) printf("medit:   spackle => '#'\n");
  printf("medit:    ...\n"); 

  for (i=ws2;i<ynum-1-ws2;i++)
  for (j=ws2;j<xnum-1-ws2;j++) {
    numvox = numnz = numz = 0;
    for (dk = -ws2;dk<=ws2;dk++)
    for (di = -ws2;di<=ws2;di++)
    for (dj = -ws2;dj<=ws2;dj++) {
      f = im_b[k+dk][i+di][j+dj];
      s = dk*c+di*b+dj*a;
      numvox++;
      if (f!=0) numnz++;
      else      numz++;
    }
    if ((im_b[k][i][j]==0 && numnz>=cfracnz*SQR(2*ws2+1)) ||
        (im_b[k][i][j]!=0 && numz>=cfracz*SQR(2*ws2+1))) {
      maxvar = -1000000;
      minvar = 1000000;
      maxi = mini = -1;
      for (m=0;m<ncor;m++) {
        a = xcor[m];
        b = ycor[m];
        c = zcor[m];
        sum = sum2 = n = 0;
        for (u=0;u<2*ws2+1;u++)
          sumv[u] = sum2v[u] = nv[u] = 0;
        for (dk = -ws2;dk<=ws2;dk++)
        for (di = -ws2;di<=ws2;di++)
        for (dj = -ws2;dj<=ws2;dj++) {
          u = ws2+floor(dk*c+di*b+dj*a+0.5);
          u = (u<0)?0:(u>2*ws2+1-1)?2*ws2+1-1:u;
          f = im_b[k+dk][i+di][j+dj];
          sum2v[u] += f*f;
          sumv[u] += f;
          nv[u] += 1;
          sum2 += f*f;
          sum += f;
          n += 1;
        }
        avg = sum/n;
        var = sum2/n-avg*avg;
        tvar = 0;
        for (u=0;u<2*ws2+1;u++) {
          avgv[u] = sumv[u]/nv[u];
          varv[u] = sum2v[u]/nv[u]-avgv[u]*avgv[u];
          tvar += varv[u];
        }
        tvar /= (2*ws2+1);
        if (tvar>maxvar) {maxvar=tvar;maxi=m;}
        if (tvar<minvar) {minvar=tvar;mini=m;}
      }
      a = xcor[mini];
      b = ycor[mini];
      c = zcor[mini];

      numvox = numnz = numz = sum = sum2 = 0;
      for (dk = -ws2;dk<=ws2;dk++)
      for (di = -ws2;di<=ws2;di++)
      for (dj = -ws2;dj<=ws2;dj++) {
        f = im_b[k+dk][i+di][j+dj];
        f2 = GetVoxelValue ( gAuxAnatomicalVolume, j+dj, i+di, k+dk );
        s = dk*c+di*b+dj*a;
        if (fabs(s)<=0.5) {
          numvox++;
          sum2 += f2;
          if (f!=0) { numnz++; sum += f; }
          else      { numz++; }
        }
      }
      if (numnz!=0) sum /= numnz;
      if (numvox!=0) sum2 /= numvox;
      f = im_b[k][i][j];
      f2 = GetVoxelValue ( gAuxAnatomicalVolume, j, i, k );
      if (flossflag && f!=0 && numz/numvox>cfracz && f<=fgray_hilim) {
        f=0; printf("."); }
      else if (spackleflag && f==0 && numnz/numvox>cfracnz) {
        f=sum; printf("#"); }
      fill[k][i][j] = f;
    }
  }
  printf("\n");
  /*for (k=0;k<numimg-1;k++)*/
  k = k2;
  for (i=0;i<ynum-1;i++)
  for (j=0;j<xnum-1;j++) {
    im_b[k][i][j] = fill[k][i][j];
    if (im_b[k][i][j]>fwhite_hilim || im_b[k][i][j]<fwhite_lolim)
      im_b[k][i][j]=0;
  }

  /*for (k=0;k<numimg-1;k++)*/
  k = k2;
  for (k=k2-ws2;k<=k2+ws2;k++)
  for (i=0;i<ynum-1;i++)
  for (j=0;j<xnum-1;j++) {
    if ( k == k2 ) 
      SetVoxelValue ( gAuxAnatomicalVolume, j, i, k, im_b[k][i][j] );
    else       
      SetVoxelValue ( gAuxAnatomicalVolume, j, i, k, 0 );
  }
  printf("medit: wmfiltered slice put in 2nd image set (can't be saved)\n");
  redraw();
}

void
norm_allslices(int normdir)
{
  int i,j,k;
  int x,y;
  float imf[256][256];
  float flim0,flim1,flim2,flim3;
  mri_tOrientation theOrientation;
  int plane;

  // get the plane.
  MWin_GetOrientation ( gMeditWindow, &theOrientation );
  plane = (int)theOrientation;

  x = y = 0;
  if ((plane==CORONAL && normdir==POSTANT) ||
      (plane==SAGITTAL && normdir==RIGHTLEFT) ||
      (plane==HORIZONTAL && normdir==INFSUP)) {
    printf("medit: ### N.B.: norm gradient not visible in this view\n");
  }
  printf("medit: normalizing all slices...\n");
  flim0 = (float)lim0;
  flim1 = (float)lim1;
  flim2 = (float)lim2;
  flim3 = (float)lim3;
  for (k=0;k<numimg;k++) {
    printf(".");
    for (i=0;i<256;i++)
    for (j=0;j<256;j++) {
      if (normdir==POSTANT)   {x=i; y=255-k;}  /* SAG:x */
      if (normdir==INFSUP)    {x=j;     y=i;}  /* COR:y */
      if (normdir==RIGHTLEFT) {x=k; y=255-j;}  /* HOR:x */
      imf[y][x] = GetVoxelValue ( gAnatomicalVolume, j, i, k );
      if ((255-y)<=flim0)
        imf[y][x]*=ffrac0;
      else if ((255-y)<=flim1)
        imf[y][x]*=(ffrac0+((255-y)-flim0)*(ffrac1-ffrac0)/(flim1-flim0));
      else if ((255-y)<=flim2)
        imf[y][x]*=(ffrac1+((255-y)-flim1)*(ffrac2-ffrac1)/(flim2-flim1));
      else if ((255-y)<=flim3)
        imf[y][x]*=(ffrac2+((255-y)-flim2)*(ffrac3-ffrac2)/(flim3-flim2));
      else
        imf[y][x]*=ffrac3;
      if (imf[y][x]>255) imf[y][x]=255;
      SetVoxelValue ( gAnatomicalVolume, j, i, k, floor(imf[y][x]+0.5) );
    }
    fflush(stdout);
  }
  printf("\n");
  printf("medit: done (to undo: quit w/o SAVEIMG)\n");
  for (k=0;k<numimg;k++)
    changed[k] = TRUE;
  editedimage = TRUE;
  redraw();
}

void norm_slice(int imc, int ic,int jc, int normdir)
{
  int k,i,j;
  int x,y;
  int k0,k1,i0,i1,j0,j1;
  int imc0,ic0,jc0;
  float imf[256][256];
  float flim0,flim1,flim2,flim3;
  mri_tOrientation theOrientation;
  int plane;

  // get the plane.
  MWin_GetOrientation ( gMeditWindow, &theOrientation );
  plane = (int)theOrientation;

  x = y = 0;
  k0 = k1 = i0 = i1 = j0 = j1 = 0;
  if ((plane==CORONAL && normdir==POSTANT) ||
      (plane==SAGITTAL && normdir==RIGHTLEFT) ||
      (plane==HORIZONTAL && normdir==INFSUP)) {
    printf("medit: ### N.B.: norm gradient not visible in this view\n");
  }

  if (!second_im_allocated)
    alloc_second_im();

  flim0 = (float)lim0;
  flim1 = (float)lim1;
  flim2 = (float)lim2;
  flim3 = (float)lim3;
  imc0 = imc/zf;
  ic0 = (ydim-1-ic)/zf;
  jc0 = jc/zf;
  if (plane==CORONAL)   {k0=imc0; k1=imc0+1; i0=0;   i1=256;   j0=0;  j1=256;}
  if (plane==SAGITTAL)  {k0=0;    k1=numimg; i0=0;   i1=256;   j0=jc0;j1=jc0+1;}
  if (plane==HORIZONTAL){k0=0;    k1=numimg; i0=ic0; i1=ic0+1; j0=0;  j1=256;}
  /*printf("k0=%d k1=%d i0=%d i1=%d j0=%d j1=%d\n",k0,k1,i0,i1,j0,j1);*/
  for (k=k0;k<k1;k++)
  for (i=i0;i<i1;i++)
  for (j=j0;j<j1;j++) {
    if (normdir==POSTANT)   {x=i; y=255-k;}  /* SAG:x */
    if (normdir==INFSUP)    {x=j;     y=i;}  /* COR:y */
    if (normdir==RIGHTLEFT) {x=k; y=255-j;}  /* HOR:x */
    imf[y][x] = GetVoxelValue ( gAnatomicalVolume, j, i, k );
    if ((255-y)<=flim0)
      imf[y][x]*=ffrac0;
    else if ((255-y)<=flim1)
      imf[y][x]*=(ffrac0+((255-y)-flim0)*(ffrac1-ffrac0)/(flim1-flim0));
    else if ((255-y)<=flim2)
      imf[y][x]*=(ffrac1+((255-y)-flim1)*(ffrac2-ffrac1)/(flim2-flim1));
    else if ((255-y)<=flim3)
      imf[y][x]*=(ffrac2+((255-y)-flim2)*(ffrac3-ffrac2)/(flim3-flim2));
    else
      imf[y][x]*=ffrac3;
    if (imf[y][x]>255) imf[y][x]=255;
    SetVoxelValue ( gAuxAnatomicalVolume, j, i, k, floor(imf[y][x]+0.5) );
  }
  printf("medit: normalized slice put in 2nd image set (can't be saved)\n");
  redraw();
}

void SnapshotVolume () {

  int             nSlice     = 0;
  tVolumeValue * pSrcSlice  = NULL;
  tVolumeValue * pDestSlice = NULL;

  /* if we already have a snapshot volume, delete it */
  if( NULL != gSnapshotVolume ) {
    DeleteVolume( &gSnapshotVolume );
  }

  /* allocate a new volume */
  InitVolume( &gSnapshotVolume, gVolumeDimension );

  /* copy the main volume into it */
  for( nSlice = 0; nSlice < gVolumeDimension; nSlice ++ ) {
    pSrcSlice  = GetVolumeSlicePtr( gAnatomicalVolume, nSlice );
    pDestSlice = GetVolumeSlicePtr( gSnapshotVolume, nSlice );
    memcpy( pDestSlice, pSrcSlice, gVolumeDimension*gVolumeDimension );
  }

  OutputPrint "Click!\n" EndOutputPrint;
}

void RestoreVolumeFromSnapshot () {

  int             nSlice     = 0;
  tVolumeValue * pSrcSlice  = NULL;
  tVolumeValue * pDestSlice = NULL;

  /* make sure we have a snapshot volume */
  if( NULL == gSnapshotVolume ) {
    OutputPrint "Must take a snapshot of a volume before restoring it.\n"
      EndOutputPrint;
    goto cleanup;
  }
  
  /* copy it into the volume */
  for( nSlice = 0; nSlice < gVolumeDimension; nSlice ++ ) {
    pSrcSlice  = GetVolumeSlicePtr( gSnapshotVolume, nSlice );
    pDestSlice = GetVolumeSlicePtr( gAnatomicalVolume, nSlice );
    memcpy( pDestSlice, pSrcSlice, gVolumeDimension*gVolumeDimension );
  }

  /* big redraw */
  MWin_RedrawAll( gMeditWindow );

 cleanup:

  return;
}

void ThresholdVolume( tVolumeRef ipVolume,
          int        inLevel,
          tBoolean   ibAbove,
          int        inNewLevel ) {

  int           nX       = 0;
  int           nY       = 0;
  int           nZ       = 0;
  tVolumeValue ucValue  = 0;

  /* check the volume */
  if( NULL == ipVolume ) 
    goto cleanup;

  /* make sure the level is in bounds */
  if( inLevel < knMinVolumeValue
      || inLevel > knMaxVolumeValue ) {
    OutputPrint "ERROR: The threshold value is out of bounds. Please choose a level between %d and %d.\n", knMinVolumeValue, knMaxVolumeValue EndOutputPrint;
    goto cleanup;
  }

  /* make sure the new level is valid */
  if( inNewLevel < knMinVolumeValue
      || inNewLevel > knMaxVolumeValue ) {
    OutputPrint "ERROR: The new value is out of bounds. Please choose a value between %d and %d.\n", knMinVolumeValue, knMaxVolumeValue EndOutputPrint;
    goto cleanup;
  }

  /* step through the volume... */
  for( nZ = 0; nZ < gVolumeDimension; nZ++ ) {
    for( nY = 0; nY < gVolumeDimension; nY++ ) {
      for( nX = 0; nX < gVolumeDimension; nX++ ) {
    
  ucValue = GetVoxelValue( ipVolume, nX, nY, nZ );
  
  /* if we're going above and this value is above the thresh, or if we're
     going below and this value is below the thresh...*/
  if( (ibAbove && ucValue > inLevel)
      || (!ibAbove && ucValue < inLevel ) ) {

    /* set the value to the new level. */
    SetVoxelValue( ipVolume, nX, nY, nZ, inNewLevel );
  }
      }
    }
  }

  /* redraw everything */
  MWin_RedrawAll( gMeditWindow );

 cleanup:

  return;
}

void ChangeVolumeAxes( tVolumeRef iVolume,
           char*      isNewXAxis,
           char*      isNewYAxis,
           char*      isNewZAxis ) {

  tVolumeRef   tmpVolume                = NULL;
  int          nX                       = 0;
  int          nY                       = 0;
  int          nZ                       = 0;
  char         acNewAxis[tkm_knNumAxes] = { 'x', 'y', 'z' };
  tBoolean     abFlip[tkm_knNumAxes]    = { FALSE, FALSE, FALSE };
  int          anCoord[tkm_knNumAxes]   = { 0, 0, 0 };
  int          nCoord                   = 0;
  tVolumeValue value                    = 0;

  /* parse the arguments */
  if( isNewXAxis[0] == '-' )
    abFlip[tkm_tAxis_X] = TRUE;
  if( isNewYAxis[0] == '-' )
    abFlip[tkm_tAxis_Y] = TRUE;
  if( isNewZAxis[0] == '-' )
    abFlip[tkm_tAxis_Z] = TRUE;
  acNewAxis[tkm_tAxis_X] = isNewXAxis[ strlen(isNewXAxis) - 1 ];
  acNewAxis[tkm_tAxis_Y] = isNewYAxis[ strlen(isNewYAxis) - 1 ];
  acNewAxis[tkm_tAxis_Z] = isNewZAxis[ strlen(isNewZAxis) - 1 ];

  /* allocate a new volume */
  InitVolume( &tmpVolume, gVolumeDimension );

  /* copy the current volume into it */
  memcpy( tmpVolume, iVolume, 
    pow( gVolumeDimension, 3 ) * sizeof( tVolumeValue ) );

  /* step through the volume... */
  for( nZ = 0; nZ < gVolumeDimension; nZ++ ) {
    for( nY = 0; nY < gVolumeDimension; nY++ ) {
      for( nX = 0; nX < gVolumeDimension; nX++ ) {
    
  for( nCoord = 0; nCoord < tkm_knNumAxes; nCoord++ ) {
    switch( acNewAxis[nCoord] ) {
    case 'x':
      anCoord[nCoord] = nX;
      break;
    case 'y':
      anCoord[nCoord] = nY;
      break;
    case 'z':
      anCoord[nCoord] = nZ;
      break;
    }
    if( abFlip[nCoord] )
      anCoord[nCoord] = gVolumeDimension-1 - anCoord[nCoord];
  }
  
  value = GetVoxelValue( tmpVolume, 
             anCoord[0], anCoord[1], anCoord[2] );
  SetVoxelValue( iVolume, nX, nY, nZ, value );
  
  
      }
    }
  }
  
  DeleteVolume( &tmpVolume );

  MWin_RedrawAll( gMeditWindow );
}

#define RotateCoords(inX,inY,onX,onY)         \
  if( ( inX <= nHalf && inY <= nHalf )        \
      || ( inX > nHalf && inY > nHalf ) ) {   \
    onX = inX;                                \
    onY = gVolumeDimension-1 - inY;           \
  } else {                                    \
    onX = gVolumeDimension-1 - inY;           \
    onY = inX;                                \
  }

void RotateVolume ( tVolumeRef iVolume,
        tkm_tAxis  iAxis,
        float      ifRadians ) {

  tVolumeRef   tmpVolume = NULL;
  float        fDegrees = 0;
  int          nX        = 0;
  int          nY        = 0;
  int          nZ        = 0;
  MATRIX*      mTransform = NULL;
  mriTransformRef transform = NULL;
  xVoxel       srcCoord;
  xVoxel       dstCoord;
  tVolumeValue value     = 0;

  /* convert to degrees */
  fDegrees = ifRadians * (PI / 180.0);

  /* allocate a new volume */
  InitVolume( &tmpVolume, gVolumeDimension );

  /* copy the current volume into it */
  memcpy( tmpVolume, iVolume, 
    pow( gVolumeDimension, 3 ) * sizeof( tVolumeValue ) );

  /* get the rot matrix */
  Trns_New( &transform );
  mTransform = MatrixIdentity( 4, NULL );
  Trns_CopyAtoRAS( transform, mTransform );
  Trns_CopyBtoRAS( transform, mTransform );
  MatrixFree( &mTransform );
  mTransform = MatrixAllocRotation( 4, fDegrees, (int)iAxis);
  Trns_CopyARAStoBRAS( transform, mTransform );
  MatrixFree( &mTransform );

  /* step through the volume... */
  for( nZ = 0; nZ < gVolumeDimension; nZ++ ) {
    for( nY = 0; nY < gVolumeDimension; nY++ ) {
      for( nX = 0; nX < gVolumeDimension; nX++ ) {
    
  /* and transform every coord to get a new value */
  xVoxl_Set( &srcCoord, nX, nY, nZ );
  Trns_ConvertAtoB( transform, &srcCoord, &dstCoord );
  value = GetVoxelValue( tmpVolume, xVoxl_ExpandInt( &dstCoord ) );
  SetVoxelValue( iVolume, nX, nY, nZ, value );
      }
    }
    fprintf( stdout, "\rRotating: %.2f%% done", 
       ((float)nZ / (float)gVolumeDimension) * 100.0 );
    fflush( stdout );
  }
  
  fprintf( stdout, "\rRotating: 100%% done     \n" );

  DeleteVolume( &tmpVolume );
  Trns_Delete( &transform );

  MWin_RedrawAll( gMeditWindow );
}

void
alloc_second_im(void)
{
  int i,k;

  InitVolume ( &gAuxAnatomicalVolume, xnum );

  for (k=0;k<6;k++) {
    sim2[k] = (tVolumeValue **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++) {
      sim2[k][i] = (tVolumeValue *)lcalloc(xnum,sizeof(char));
    }
  }
  second_im_allocated = TRUE;
}



/* =========================================================== TCL WRAPPERS */

int TclLoadHeadPts ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadHeadPts points_file transform_file",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    DebugPrint "LoadHeadPts %s %s\n",
      argv[1], argv[2] EndDebugPrint;
    LoadHeadPts ( argv[1], argv[2] ); 
  }

  return TCL_OK;
}

int TclRotateHeadPts ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: RotateHeadPts degrees x,y,z",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    RotateHeadPts ( atof(argv[1]), argv[2][0] ); 
  }

  return TCL_OK;
}

int TclTranslateHeadPts ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: TranslateHeadPts distance x,y,z",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    TranslateHeadPts ( atof(argv[1]), argv[2][0] ); 
  }

  return TCL_OK;
}


int TclRestoreHeadPts ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: RestoreHeadPts",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    RestoreHeadPoints();
  }

  return TCL_OK;
}

int TclWriteHeadPointsTransform ( ClientData inClientData, 
          Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WriteHeadPointsTransform",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    WriteHeadPointsTransform();
  }

  return TCL_OK;
}

int TclWriteHeadPointsFile ( ClientData inClientData, 
          Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WriteHeadPointsFile",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    WriteHeadPointsFile();
  }

  return TCL_OK;
}

int TclSetSelectedHeadPointLabel ( ClientData inClientData, 
           Tcl_Interp* inInterp,
           int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetSelectedHeadPointLabel",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    SetSelectedHeadPointLabel( argv[1] );
  }

  return TCL_OK;
}

int TclAlignSelectedHeadPointToAnaIdx ( ClientData inClientData, 
          Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  xVoxel anaIdx;

 if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: AlignSelectedHeadPointToAnaIdx",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    xVoxl_Set( &anaIdx, atoi( argv[1] ), atoi( argv[2] ), atoi( argv[3] ) );
    AlignSelectedHeadPointToAnaIdx( &anaIdx );
  }

  return TCL_OK;
}

int TclSaveRGB ( ClientData inClientData, Tcl_Interp* inInterp,
     int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveRGB filename:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    save_rgb ( argv[1] );
  }

  return TCL_OK;
}

int TclWMFilterCorSlice ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

 xVoxelRef theCursor = NULL;

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WMFilterCorSlice",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    xVoxl_New ( &theCursor );
    MWin_GetCursor ( gMeditWindow, theCursor );
    wmfilter_corslice ( xVoxl_GetZ(theCursor) );
    xVoxl_Delete ( &theCursor );
  }

  return TCL_OK;
}

int TclNormSlice ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

 xVoxelRef theCursor = NULL;

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: NormSlice {0=PostAnt,1=InfSup,2=LeftRight}",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    xVoxl_New ( &theCursor );
    MWin_GetCursor ( gMeditWindow, theCursor );
    norm_slice ( xVoxl_GetZ(theCursor), 
     xVoxl_GetY(theCursor),
     xVoxl_GetX(theCursor),
     atoi(argv[1]) ); 
    xVoxl_Delete ( &theCursor );
  }

  return TCL_OK;
}

int TclNormAllSlices ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, 
         "wrong # args: NormAllSlices {0=PostAnt,1=InfSup,2=LeftRight}",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    norm_allslices ( atoi(argv[1]) );
  }  

  return TCL_OK;
}

int TclThresholdVolume ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: ThresholdVolume threshold_value {0=below,1=above} new_value",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    ThresholdVolume( gAnatomicalVolume, 
         atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) );
  }  

  return TCL_OK;
}

int TclFlipVolume ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: FlipVolume 1|0 1|0 1|0", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    ChangeVolumeAxes( gAnatomicalVolume, 
          atoi(argv[1])?"-x":"x", 
          atoi(argv[2])?"-y":"y", 
          atoi(argv[3])?"-z":"z" );
  }  

  return TCL_OK;
}

int TclRotateVolume ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {

  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: RotateVolume x|y|z degrees", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    RotateVolume( gAnatomicalVolume, 
      argv[1][0] == 'x' ? tkm_tAxis_X :
      argv[1][0] == 'y' ? tkm_tAxis_Y :
      tkm_tAxis_Z,
      atof( argv[2] ) );
  }  

  return TCL_OK;
}

int TclLoadVolume ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadVolume image_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    
    read_images( argv[1] );
  }

  return TCL_OK;
}

int TclLoadAuxVolume ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadAuxVolume image_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    //    DebugPrint "Called LoadAuxVolume %s\n", argv[1] EndDebugPrint;
    read_second_images ( argv[1] );
  }

  return TCL_OK;
}

int TclSaveVolume ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveVolume",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    write_images( gVolumeSrc );
  }

  return TCL_OK;
}

int TclSaveVolumeAs ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveVolumeAs volume_path",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    write_images ( argv[1] );
  }

  return TCL_OK;
}

int TclSaveAuxVolume ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveAuxVolume",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  Tcl_SetResult ( inInterp, "SaveAuxVolume: not implemented yet",
      TCL_VOLATILE );
  return TCL_OK;
}

int TclSnapshotVolume ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SnapshotVolume",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    SnapshotVolume();
  }

  return TCL_OK;
}

int TclRestoreVolumeFromSnapshot ( ClientData inClientData, 
           Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: RestoreVolumeFromSnapshot",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    RestoreVolumeFromSnapshot();
  }

  return TCL_OK;
}

int TclSetVolumeColorScale ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetVolumeColorScale volume threshold:float squash:float",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    SetVolumeColorScale ( atoi( argv[1] ), atof( argv[2] ), atof( argv[3] ) );
  }

  return TCL_OK;
}

int TclSaveLabel ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveLabel name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    SaveSelectionToLabelFile ( argv[1] );
  }

  return TCL_OK;
}

int TclLoadLabel ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadLabel name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    LoadSelectionFromLabelFile ( argv[1] );
  }

  return TCL_OK;
}

int TclGraphSelectedRegion ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: GraphSelectedRegion",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    GraphSelectedRegion ();
  }

  return TCL_OK;
}

int TclClearSelection ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ClearSelection",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    ClearSelection ();
  }

  return TCL_OK;
}


int TclSendCursor ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

 xVoxelRef theCursor = NULL;

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SendCursor",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) { 
    xVoxl_New ( &theCursor );
    MWin_GetCursor ( gMeditWindow, theCursor );
    WriteVoxelToEditFile ( theCursor );
    xVoxl_Delete ( &theCursor );
  }    

  return TCL_OK;
}

int TclReadCursor ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadCursor",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    ReadCursorFromEditFile ();
  }

  return TCL_OK;
}


int TclUndoLastEdit ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: UndoLastEdit",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    RestoreUndoList ();
  }

  return TCL_OK;
}

int TclNewControlPoint ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

 xVoxelRef theCursor = NULL;

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: NewControlPoint",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    xVoxl_New ( &theCursor );
    MWin_GetCursor ( gMeditWindow, theCursor );
    NewCtrlPt ( theCursor, TRUE );
    xVoxl_Delete ( &theCursor );
  }

  return TCL_OK;
}

int TclDeselectAllControlPoints ( ClientData inClientData, 
          Tcl_Interp* inInterp,
          int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: DeselectAllControlPoints",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    DeselectAllCtrlPts ();
  }

  return TCL_OK;
}

int TclDeleteSelectedControlPoints ( ClientData inClientData,
             Tcl_Interp* inInterp,
             int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: DeleteSelectedControlPoints",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    DeleteSelectedCtrlPts ();
  }  

  return TCL_OK;
}

int TclWriteControlPointFile ( ClientData inClientData, Tcl_Interp* inInterp,
             int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WriteControlPointFile",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    WriteCtrlPtFile ();
  }

  return TCL_OK;
}

int TclLoadMainSurface ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadMainSurface surface_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    LoadSurface ( argv[1] );
  }

  return TCL_OK;
}

int TclLoadCanonicalSurface ( ClientData inClientData, Tcl_Interp* inInterp,
            int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadCanonicalSurface surface_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    LoadSurfaceVertexSet( Surf_tVertexSet_Pial, argv[1] );
  }

  return TCL_OK;
}

int TclLoadOriginalSurface ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadOriginalSurface surface_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    LoadSurfaceVertexSet( Surf_tVertexSet_Original, argv[1] );
  }  

  return TCL_OK;
}

int TclUnloadAllSurfaces ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: UnloadSurface",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    UnloadSurface ();
  }  

  return TCL_OK;
}


int TclGotoMainVertex ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: GotoMainVertex vertex_num:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    GotoSurfaceVertex ( Surf_tVertexSet_Main, atoi( argv[1] ) );
  }  

  return TCL_OK;
}

int TclGotoCanonicalVertex ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: GotoCanonicalVertex vertex_num:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    GotoSurfaceVertex ( Surf_tVertexSet_Pial, atoi( argv[1] ) );
  }  

  return TCL_OK;
}

int TclGotoOriginalVertex ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: GotoOriginalVertex vertex_num:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    GotoSurfaceVertex ( Surf_tVertexSet_Original, atoi( argv[1] ) );
  }  

  return TCL_OK;
}

int TclShowNearestMainVertex ( ClientData inClientData, Tcl_Interp* inInterp,
             int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestMainVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    FindNearestSurfaceVertex ( Surf_tVertexSet_Main );
  }  

  return TCL_OK;
}

int TclShowNearestOriginalVertex ( ClientData inClientData, 
           Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestOriginalVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    FindNearestSurfaceVertex ( Surf_tVertexSet_Original );
  }  

  return TCL_OK;
}

int TclShowNearestCanonicalVertex ( ClientData inClientData, 
            Tcl_Interp* inInterp,
            int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestCanonicalVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    FindNearestSurfaceVertex ( Surf_tVertexSet_Pial );
  }  

  return TCL_OK;
}


int TclLoadParcellationVolume ( ClientData inClientData, 
        Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: LoadParcellationVolume directory_and_prefix:string color_file:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    LoadParcellationVolume ( argv[1], argv[2] );
  }  

  return TCL_OK;
}
int TclQuitMedit ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: QuitMedit", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    tkm_Quit ();
  }

  return TCL_OK;
}

/*=======================================================================*/

/* for tcl/tk */
static void StdinProc _ANSI_ARGS_((ClientData clientData, int mask));
static void Prompt _ANSI_ARGS_((Tcl_Interp *interp, int partial));
//static Tk_Window mainWindow;
static Tcl_Interp *interp;
static Tcl_DString command;
static int tty;

static char sTclEnvVar[256] = "";
static char sTkEnvVar[256] = "";


int main(argc, argv)   /* new main */
int argc;
char **argv;
{
  int code;
  int scriptok=FALSE;
  int cur_arg;
  //  static char *display = NULL;
  char tkmedit_tcl[NAME_LENGTH];
  char script_tcl[NAME_LENGTH];
  char *envptr;
  FILE *fp ;
  x3Lst_tErr e3DList     = x3Lst_tErr_NoErr;
  xList_tErr eList       = xList_tErr_NoErr;
  FunV_tErr  eFunctional = FunV_tErr_NoError;
  MWin_tErr  eWindow     = MWin_tErr_NoErr;
  int        nArg        = 0;
  time_t     theTime;

#ifdef SET_TCL_ENV_VAR
  tBoolean  bChangedEnvVar    = FALSE;
  char*     sTclLib           = NULL;
  char      sSavedTclLib[256] = "";
  int       nTclLength        = 0;
  char      sNewTclLib[256]   = "";
  char*     sTkLib            = NULL;
  char      sSavedTkLib[256]  = "";
  int       nTkLength         = 0;
  char      sNewTkLib[256]    = "";
#endif

  // init our debugging macro code, if any.
  xDbg_Init();
  EnableDebuggingOutput;

  //  xDbg_PrintStatus();

  time( &theTime );
  DebugPrint "tkmedit started: %s\n\t", ctime( &theTime ) EndDebugPrint;
  for( nArg = 0; nArg < argc; nArg++ ) {
    DebugPrint "%s ", argv[nArg] EndDebugPrint;
  }
  DebugPrint "\n\n" EndDebugPrint;

  /* init medit window */
  glutInit            ( &argc, argv );
  glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
  eWindow = MWin_New ( &gMeditWindow, "", 512, 512 );
  if( MWin_tErr_NoErr != eWindow ) {
    OutputPrint "Fatal error: Couldn't make main window.\n" EndOutputPrint;
    exit( 1 );
  }

  // init the selection list 
  eList = xList_New( &gSelectionList );
  if( xList_tErr_NoErr != eList ) {
    OutputPrint "Fatal error: Couldn't initialize selection list.\n"
      EndOutputPrint;
    exit( 1 );
  }
  xList_SetComparator( gSelectionList, CompareVoxels );

  // init our control pt list
  e3DList = x3Lst_New( &gCtrlPtList, 256 );
  if( x3Lst_tErr_NoErr != e3DList ) {
    OutputPrint "Fatal error: Couldn't initialize control point list.\n"
      EndOutputPrint;
    exit( 1 );
  }
  x3Lst_SetComparator( gCtrlPtList, CompareVoxels );

  // init the undo list.
  InitUndoList ();

  // and the selection module.
  InitSelectionModule ();

  /* create functional volume */
  eFunctional = FunV_New( &gFunctionalVolume,
        UpdateAndRedraw, tkm_SendTclCommand, 
        SendTCLCommand );

  if( FunV_tErr_NoError != eFunctional ) {
    OutputPrint "Couldn't initialize functional module.\n" EndOutputPrint;
    exit( 1 );
  }


  /* set windows data sources */
  MWin_SetControlPointsSpace ( gMeditWindow, -1, gCtrlPtList );
  MWin_SetControlPointsSelectionList ( gMeditWindow, -1, gSelectionList );
  MWin_SetSelectionSpace ( gMeditWindow, -1, gSelectedVoxels );


  /* start by disabling a bunch of stuff. if it gets loaded, the options
     will be enabled later */
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceLoadingOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceViewingOptions,"0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowCanonicalSurfaceViewingOptions,"0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowOriginalSurfaceViewingOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowHeadPointLabel, "0" );

  /* get cwd */
  if( !FindUserHomeDir() ) {
    OutputPrint "FATAL ERROR: Couldn't get current working directory!\n"
      EndOutputPrint;
    exit( 1 );
  }

  /* start program, now as function; gl window not opened yet */
  Medit((ClientData) NULL, interp, argc, argv); /* event loop commented out */


  /* get tkmedit tcl startup script location from environment */
  envptr = getenv("MRI_DIR");
  if (envptr==NULL) {
    printf("tkmedit: env var MRI_DIR undefined (use setenv)\n");
    printf("    [dir containing mri distribution]\n");
    exit(1);
  }
#ifdef USE_LICENSE
  checkLicense(envptr);
#endif

  fp = NULL;


  /* if there is no interface name... */
  if( MATCH( gInterfaceScriptName, "" ) ) {

    /* use local */
    sprintf(tkmedit_tcl,"%s","tkmedit.tcl"); 
    fp = fopen ( tkmedit_tcl,"r" );

    /* if not open, try in lib/tcl */
    if ( NULL == fp ) { 
      sprintf(tkmedit_tcl,"%s/lib/tcl/%s",envptr,"tkmedit.tcl"); 
      fp = fopen ( tkmedit_tcl,"r" );
    }

  } else {
    
    /* copy in interface script location */
    strcpy( tkmedit_tcl, gInterfaceScriptName );
    fp = fopen ( tkmedit_tcl,"r" );
  }

  // if file still not found bail out.
  if ( NULL == fp ) {
    printf("tkmedit: script %s not found\n",tkmedit_tcl);
    exit(1);
  }
  
  OutputPrint "Using interface file %s\n", tkmedit_tcl EndOutputPrint;

  // not acutally using it now, just checking for it.
  fclose(fp);

  /* look for script: (1) cwd, (2) MRI_DIR/lib/tcl */
  strcpy(script_tcl,"newscript.tcl");   /* default if no script */
  for (cur_arg=1; cur_arg<argc; cur_arg++) {  /* just look for -tcl arg */
    if (!MATCH(argv[cur_arg],"-tcl"))
      continue;
    cur_arg++;
    if (cur_arg==argc) {
      printf("tkmedit: ### no script name given\n");
      exit(1);
    }
    if (argv[cur_arg][0]=='-') {
      printf("tkmedit: ### option is bad script name: %s\n",argv[cur_arg]);
      exit(1);
    }
    strcpy(script_tcl,argv[cur_arg]);
    fp = fopen(script_tcl,"r");  /* first look in cwd */
    if (fp==NULL) {
      sprintf(script_tcl,"%s/lib/tcl/%s",envptr,argv[cur_arg]);/* then libdir */      fp = fopen(script_tcl,"r");
      if (fp==NULL) {
        printf(
          "tkmedit: (1) File ./%s not found\n",argv[cur_arg]);
        printf(
          "         (2) File %s/lib/tcl/%s not found\n",envptr,argv[cur_arg]);
        exit(1);
      }
    }
    scriptok = TRUE;
  }

  /* process ctrl pts file */
  ProcessCtrlPtFile();

  /* set window's data */
  MWin_SetWindowTitle( gMeditWindow, gSubjectName );

  /* don't accept tcl commands yet. */
  gbAcceptingTclCommands = FALSE;

  /* ============================================================ TCL INIT */

#if SET_TCL_ENV_VAR
  /* get tcl and tk lib env vars */
  sTclLib = getenv( "TCL_LIBRARY" );
  sTkLib  = getenv( "TK_LIBRARY" );

  /* if we got them... */
  if( NULL != sTclLib 
      && NULL != sTkLib ) {

    /* save them */
    strcpy( sSavedTclLib, sTclLib );
    strcpy( sSavedTkLib, sTkLib );

    /* check out major version number */
    strcpy( sNewTclLib, sTclLib );
    nTclLength = strlen( sNewTclLib );
    strcpy( sNewTkLib, sTkLib );
    nTkLength = strlen( sNewTkLib );

    if( sNewTclLib[nTclLength-3] != '8' ) {
      
      /* we changed it */
      bChangedEnvVar = TRUE;

      /* set verseion to 8.3 */
      sNewTclLib[nTclLength-3] = '8';
      sNewTclLib[nTclLength-1] = '3';
      sNewTkLib[nTkLength-3] = '8';
      sNewTkLib[nTkLength-1] = '3';
    
      /* set env variable */
      sprintf( sTclEnvVar, "%s=%s", "TCL_LIBRARY", sNewTclLib );
      if( putenv( sTclEnvVar ) ) {
  OutputPrint "ERROR: Couldn't set TCL_LIBRARY env var.\n" 
    EndOutputPrint;
  exit( 1 );
      }
      sprintf( sTkEnvVar, "%s=%s", "TK_LIBRARY", sNewTkLib );
      if( putenv( sTkEnvVar ) ) {
  OutputPrint "ERROR: Couldn't set TK_LIBRARY env var.\n" 
    EndOutputPrint;
  exit( 1 );
      }
    }

  } else {

    OutputPrint "ERROR: TCL_LIBRARY or TK_LIBRARY environement variable is not set.\n" EndOutputPrint;
    exit ( 1 );
  }

#endif

  /* start tcl/tk; first make interpreter */
  interp = Tcl_CreateInterp();

  // kt - set global interp
  SetTclInterp ( interp );

  /* read tcl/tk internal startup scripts */
  if (Tcl_Init(interp) == TCL_ERROR) {
    fprintf(stderr, "Tcl_Init failed: %s\n", interp->result); }
  if (Tk_Init(interp)== TCL_ERROR) {
    fprintf(stderr, "Tk_Init failed: %s\n", interp->result); }
  if (Tix_Init(interp) == TCL_ERROR ) {
    fprintf(stderr, "Tix_Init failed: %s\n", interp->result); }

#if SET_TCL_ENV_VAR
  /* restore env vars */
  if( bChangedEnvVar ) {
    sprintf( sTclEnvVar, "%s=%s", "TCL_LIBRARY", sSavedTclLib );
    if( putenv( sTclEnvVar ) ) {
      OutputPrint "ERROR: Couldn't restore TCL_LIBRARY env var.\n" 
  EndOutputPrint;
      exit( 1 );
    }
    sprintf( sTkEnvVar, "%s=%s", "TK_LIBRARY", sSavedTkLib );
    if( putenv( sTkEnvVar ) ) {
      OutputPrint "ERROR: Couldn't restore TCL_LIBRARY env var.\n" 
  EndOutputPrint;
      exit( 1 );
    }
  }
#endif

  /* regsiter the medit window tcl commands */
  MWin_RegisterTclCommands ( gMeditWindow, interp );

  /* set the "tcl_interactive" variable */
  tty = isatty(0);
  Tcl_SetVar(interp, "tcl_interactive", (tty) ? "1" : "0", TCL_GLOBAL_ONLY);


  /*=======================================================================*/
  /* register wrapped surfer functions with interpreter */

  Tcl_CreateCommand ( interp, "LoadHeadPts",
          TclLoadHeadPts,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RotateHeadPts",
          TclRotateHeadPts,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "TranslateHeadPts",
          TclTranslateHeadPts,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RestoreHeadPts",
          TclRestoreHeadPts,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "WriteHeadPointsTransform",
          TclWriteHeadPointsTransform,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "WriteHeadPointsFile",
          TclWriteHeadPointsFile,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SetSelectedHeadPointLabel",
          TclSetSelectedHeadPointLabel,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "AlignSelectedHeadPointToAnaIdx",
          TclAlignSelectedHeadPointToAnaIdx,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveRGB",
          TclSaveRGB,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "WMFilterCorSlice",
          TclWMFilterCorSlice,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "NormSlice",
          TclNormSlice,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "NormAllSlices",
          TclNormAllSlices,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ThresholdVolume",
          TclThresholdVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "FlipVolume",
          TclFlipVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RotateVolume",
          TclRotateVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadVolume",
          TclLoadVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadAuxVolume",
          TclLoadAuxVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveVolume",
          TclSaveVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveVolumeAs",
          TclSaveVolumeAs,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveAuxVolume",
          TclSaveAuxVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SnapshotVolume",
          TclSnapshotVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RestoreVolumeFromSnapshot",
          TclRestoreVolumeFromSnapshot,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SetVolumeColorScale",
          TclSetVolumeColorScale,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveLabel",
          TclSaveLabel,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadLabel",
          TclLoadLabel,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "GraphSelectedRegion",
          TclGraphSelectedRegion,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ClearSelection",
          TclClearSelection,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "UndoLastEdit",
          TclUndoLastEdit,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SendCursor",
          TclSendCursor,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadCursor",
          TclReadCursor,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "NewControlPoint",
          TclNewControlPoint,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "DeselectAllControlPoints",
          TclDeselectAllControlPoints,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "DeleteSelectedControlPoints",
          TclDeleteSelectedControlPoints,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "WriteControlPointFile",
          TclWriteControlPointFile,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadMainSurface",
          TclLoadMainSurface,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadCanonicalSurface",
          TclLoadCanonicalSurface,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadOriginalSurface",
          TclLoadOriginalSurface,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "UnloadSurface",
          TclUnloadAllSurfaces,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "GotoMainVertex",
          TclGotoMainVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "GotoCanonicalVertex",
          TclGotoCanonicalVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "GotoOriginalVertex",
          TclGotoOriginalVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ShowNearestMainVertex",
          TclShowNearestMainVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ShowNearestOriginalVertex",
          TclShowNearestOriginalVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ShowNearestCanonicalVertex",
          TclShowNearestCanonicalVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadParcellationVolume",
          TclLoadParcellationVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "QuitMedit",
          TclQuitMedit,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  /*=======================================================================*/


  /* run tcl/tk startup script to set vars, make interface; no display yet */
  code = Tcl_EvalFile( interp, tkmedit_tcl );
  if (*interp->result != 0)  {
    fprintf( stdout, "Error reading tkmedit.tcl\n\t%s\n", interp->result );
    exit( 1 );
  }

  /* always start up command line shell too */
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  if (tty)
    Prompt(interp, 0);
  fflush(stdout);
  Tcl_DStringInit(&command);
  Tcl_ResetResult(interp);

  /*Tk_MainLoop();*/  /* standard */

  /* finish initializing functional volume */
  eFunctional = FunV_InitGraphWindow( gFunctionalVolume, interp );
  if( FunV_tErr_NoError != eFunctional ) {
    DebugPrint "Error initializing functional time course window.\n" EndDebugPrint;
    OutputPrint "Couldn't initialize functional time course interface window. Please check if file tkm_functional.tcl is in the library path.\n" EndOutputPrint;
    exit( 1 );
  }

  eFunctional = FunV_RegisterTclCommands( gFunctionalVolume, interp );
  if( FunV_tErr_NoError != eFunctional ) {
    DebugPrint "Error registering functional TCL commands.\n" EndDebugPrint;
    OutputPrint "Error building interface between functional module and interface window.\n" EndOutputPrint;
    exit( 1 );
  }

  /* if command line script exists, now run as batch job (possibly exiting) */
  if( gbRunBatchScript ) {  
    
    printf( "tkmedit: run tcl script: %s\n", gsBatchScriptName );
    code = Tcl_EvalFile( interp, gsBatchScriptName );
    if (*interp->result != 0)  printf(interp->result);

  }

  /* set data in window */
  MWin_SetOverlayVolume( gMeditWindow, -1, gFunctionalVolume );

  /* if we have a tal transform, hide the ras coords. otherwise hide
     the tal coords. */
  if ( gbAnaToTalTransformLoaded ) {
    tkm_SendTclCommand ( tkm_tTclCommand_ShowTalCoords, "1" );
    tkm_SendTclCommand ( tkm_tTclCommand_ShowRASCoords, "0" );
  } else {
    tkm_SendTclCommand ( tkm_tTclCommand_ShowTalCoords, "0" );
  }

  /* now let the window accept tcl commands. */
  gbAcceptingTclCommands = TRUE;
  MWin_AcceptTclCommands( gMeditWindow, TRUE );

  /* we probably get sent a few tcl commands before now, and we cached
     them, so send them now. */
  SendCachedTclCommands ();

  /* set the volume color scale */
  SetVolumeColorScale ( tkm_tVolumeType_Main,
      kfDefaultVolumeThreshold, kfDefaultVolumeSquash );
  SetVolumeColorScale ( tkm_tVolumeType_Aux,
      kfDefaultVolumeThreshold, kfDefaultVolumeSquash );

  /* if using csurf interface, call the func that hides a bunch of stuff. */
  if( gbUseCsurfInterface ) {
    tkm_SendTclCommand( tkm_tTclCommand_CsurfInterface, "" );
  }

  /* never returns */
  glutMainLoop ();

  /* this will never get called here, but oh well. */
  tkm_Quit ();
  
  return 0;
}

void tkm_Quit () {
  
  /* delete window */
  MWin_Delete( &gMeditWindow );

  // shut down tcl stuff
  SendTCLCommand ( "exit" );

  // delete everything we allocated before
  FunV_Delete( &gFunctionalVolume );
  DeleteSelectionModule ();
  xList_Delete( &gSelectionList );
  x3Lst_Delete( &gCtrlPtList );
  DeleteUndoList ();
  DeleteDebugging;

  exit(0);
}


/*=== from TkMain.c ===================================================*/
static void StdinProc(clientData, mask)
  ClientData clientData;
  int mask;
{
#define BUFFER_SIZE 4000
  char input[BUFFER_SIZE+1];
  static int gotPartial = 0;
  char *cmd;
  int code, count;

  count = read(fileno(stdin), input, BUFFER_SIZE);
  if (count <= 0) {
    if (!gotPartial) {
      if (tty) {Tcl_Eval(interp, "exit"); exit(1);}
      else     {Tk_DeleteFileHandler(0);}
      return;
    }
    else count = 0;
  }
  cmd = Tcl_DStringAppend(&command, input, count);
  if (count != 0) {
    if ((input[count-1] != '\n') && (input[count-1] != ';')) {
      gotPartial = 1;
      goto prompt; }
    if (!Tcl_CommandComplete(cmd)) {
      gotPartial = 1;
      goto prompt; }
  }
  gotPartial = 0;
  Tk_CreateFileHandler(0, 0, StdinProc, (ClientData) 0);
  code = Tcl_RecordAndEval(interp, cmd, TCL_EVAL_GLOBAL);
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  Tcl_DStringFree(&command);
  if (*interp->result != 0)
    if ((code != TCL_OK) || (tty))
      puts(interp->result);
  prompt:
  if (tty)  Prompt(interp, gotPartial);
  Tcl_ResetResult(interp);
}

/*=== from TkMain.c ===================================================*/
static void Prompt(interp, partial)
  Tcl_Interp *interp;
  int partial;
{
  char *promptCmd;
  int code;

  promptCmd = Tcl_GetVar(interp,
      partial ? "tcl_prompt2" : "tcl_prompt1", TCL_GLOBAL_ONLY);
  if (promptCmd == NULL) {
    defaultPrompt:
    if (!partial)
      fputs("% ", stdout);
  }
  else {
    code = Tcl_Eval(interp, promptCmd);
    if (code != TCL_OK) {
      Tcl_AddErrorInfo(interp,
                    "\n    (script that generates prompt)");
      fprintf(stderr, "%s\n", interp->result);
      goto defaultPrompt;
    }
  }
  fflush(stdout);
}



// kt

// =========================================================== READING VOLUMES

int ReadVolumeWithMRIRead ( tVolumeRef*        iopVolume, 
          General_transform* iTransform,
          char*              inFileOrPath ) {

  MRI*          theMRIVolume       = NULL;
  tVolumeRef    theVolume          = NULL;
  int           theSlice           = 0;
  tVolumeValue* theSlicePtr        = NULL;

  // pass the path to MRIRead
  theMRIVolume = MRIread ( inFileOrPath );

  // make sure the result is good.
  if ( NULL == theMRIVolume ) {
    OutputPrint "Couldn't read volume data at %s\n", 
      inFileOrPath EndOutputPrint;
    return 1;
  }

  // conform it.
  theMRIVolume = MRIconform ( theMRIVolume );
  
  // grab all the data we need.
  imnr0 = theMRIVolume->imnr0;
  imnr1 = theMRIVolume->imnr1;
  ptype = theMRIVolume->ptype;
  xnum = theMRIVolume->width;
  ynum = theMRIVolume->height;
  ps = theMRIVolume->ps;
  st = theMRIVolume->thick;
  xx0 = theMRIVolume->xstart;
  xx1 = theMRIVolume->xend;
  yy0 = theMRIVolume->ystart;
  yy1 = theMRIVolume->yend;
  zz0 = theMRIVolume->zstart;
  zz1 = theMRIVolume->zend;

  // if they want transforms, grab the tal transforms.
  if ( NULL != iTransform ) {
    if( NULL != theMRIVolume->linear_transform ) {
      copy_general_transform ( &theMRIVolume->transform, iTransform );
    }
  }

    numimg = imnr1-imnr0+1; // really the number of slices

  // calc window dimensions according to already defined scale factor
  xdim= xnum * zf;
  ydim= ynum * zf;

  InitVolume ( &theVolume, xnum );

  // read in all image data into im[], set changed[] for all slices to nil.
  for ( theSlice = 0; theSlice < numimg; theSlice++ ) {
    theSlicePtr = GetVolumeSlicePtr ( theVolume, theSlice );
    memcpy ( theSlicePtr, *theMRIVolume->slices[theSlice],
       xnum * ynum * sizeof(BUFTYPE) );
    changed[theSlice] = FALSE;
  }

  MRIfree ( &theMRIVolume );

  /* delete the incoming volume if we had one. */
  if( NULL != *iopVolume ) 
    DeleteVolume( iopVolume );

  /* return volume */
  *iopVolume = theVolume;

  return 0;
}


/* ================================================ Control point utilities */

                                       /* reads the control.dat file, 
                                          transforms all pts from RAS space 
                                          to voxel space, and adds them as 
                                          control pts */
void ProcessCtrlPtFile ( ) {

  char theFilename [NAME_LENGTH];
  FILE * theFile;
  float theTempX, theTempY, theTempZ;
  Real theRASX, theRASY, theRASZ;
  int theVoxX, theVoxY, theVoxZ;
  int theNumPtsRead;
  xVoxelRef theVoxel = NULL;
  
  xVoxl_New( &theVoxel );

  // don't parse the file if we already have.
  if ( TRUE == gParsedCtrlPtFile )
    return;

  /* make the file name */
  MakeFileName( "control.dat", tkm_tFileName_ControlPoints, theFilename );
  
  // open for reading. position file ptr at beginning of file.
  theFile = fopen ( theFilename, "r" );
  
  // check for success. if not, print an error and return.
  if ( NULL == theFile ) {
    DebugPrint "ProcessCtrlPtFile: Couldn't open %s for processing.\n",
      theFilename EndDebugPrint;
    return;
  }

  // while we have points left
  while ( !feof(theFile) ) {
    
    // read in some numbers
    theNumPtsRead = fscanf ( theFile, "%f %f %f",
           &theTempX, &theTempY, &theTempZ );

    // if not successful, file is wierd. print error, close file, and return.
    if ( theNumPtsRead < 3 &&
         theNumPtsRead != EOF ) {

      DebugPrint "ProcessCtrlPtFile: Error parsing file, expected three float values but didn't get them.\n" EndDebugPrint;
      fclose ( theFile );
      return;
    }

    // if we got 3 pts, we got a point.
    if ( 3 == theNumPtsRead ) {

      // cast to reals
      theRASX = (Real) theTempX;
      theRASY = (Real) theTempY;
      theRASZ = (Real) theTempZ;

      // transform from ras to voxel
      RASToVoxel ( theRASX, theRASY, theRASZ,
                   &theVoxX, &theVoxY, &theVoxZ );
      
      // add it to our cntrl points space
      xVoxl_Set( theVoxel, theVoxX, theVoxY, theVoxZ );
      NewCtrlPt( theVoxel, FALSE );
    }
  }    

  // close the file.
  fclose ( theFile );

  // mark that we have processed the file, and shouldn't do it again.
  gParsedCtrlPtFile = TRUE;

  xVoxl_Delete ( &theVoxel );
}

                                       /* writes all control points to the
                                          control.dat file in RAS space */
void WriteCtrlPtFile ( ) {

  char       sFilename[256] = "";
  FILE*      pFile          = NULL;
  int        nPlane         = 0;
  Real       rRASX          = 0;
  Real       rRASY          = 0;
  Real       rRASZ          = 0;
  xListRef   list           = NULL;
  xVoxelRef  voxel          = NULL;
  x3Lst_tErr e3DList        = x3Lst_tErr_NoErr;
  xList_tErr eList          = xList_tErr_NoErr;

  /* make the file name */
  MakeFileName( "control.dat", tkm_tFileName_ControlPoints, sFilename );

  // open the file for writing.
  pFile = fopen( sFilename, "w" );

  // check for success.
  if ( NULL == pFile ) {
    OutputPrint "Couldn't create %s for writing control points.\n",
      sFilename EndOutputPrint;
    return;
  }

  OutputPrint "Saving control points... " EndOutputPrint;

  // get the ctrl pts in the list...
  for ( nPlane = 0; nPlane < gVolumeDimension; nPlane++ ) {

    // get the list for this x value.
    e3DList = x3Lst_GetItemsInXPlane( gCtrlPtList, nPlane, &list );
    if( e3DList != x3Lst_tErr_NoErr )
      goto error;

    /* traverse the list */
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&voxel )) 
     != xList_tErr_EndOfList ) {

      if( voxel ) {

  // transform to ras space.
  VoxelToRAS ( xVoxl_ExpandInt(voxel),
         &rRASX, &rRASY, &rRASZ );
  
  // write to the file
  fprintf( pFile, "%f %f %f\n", rRASX, rRASY, rRASZ );
      }
    }

    if( eList != xList_tErr_EndOfList )
      goto error;
  }

  OutputPrint " done.\n" EndOutputPrint;

  goto cleanup;

 error:

  if( eList != xList_tErr_NoErr ) {
    DebugPrint "xList error %d in WriteCtrlPtFile.\n", eList EndDebugPrint;
    OutputPrint "Error saving control points.\n" EndOutputPrint;
  }

  if( e3DList != x3Lst_tErr_NoErr ) {
    DebugPrint "x3Lst error %d in WriteCtrlPtFile.\n", e3DList EndDebugPrint;
    OutputPrint "Error saving control points.\n" EndOutputPrint;
  }

 cleanup:

  /* close file */
  if( pFile )
    fclose( pFile );
}


tBoolean FindNearestCtrlPt ( xVoxelRef        inVolumeVox, 
           mri_tOrientation inPlane,
           xVoxelRef*       outCtrlPt ) {

  tBoolean     bFound           = FALSE;
  unsigned int nDistance        = 0;
  unsigned int nClosestDistance = 0;
  xListRef     list             = NULL;
  xVoxelRef    voxel            = NULL;
  xVoxelRef    closestVoxel     = NULL;
  x3Lst_tErr   e3DList          = x3Lst_tErr_NoErr;
  xList_tErr   eList            = xList_tErr_NoErr;

  /* get the list to search in */
  switch ( inPlane ) {
  case mri_tOrientation_Coronal: 
    e3DList = x3Lst_GetItemsInZPlane( gCtrlPtList, 
              xVoxl_GetZ(inVolumeVox), &list );
    break;
  case mri_tOrientation_Horizontal: 
    e3DList = x3Lst_GetItemsInYPlane( gCtrlPtList, 
              xVoxl_GetY(inVolumeVox), &list );
    break;
  case mri_tOrientation_Sagittal: 
    e3DList = x3Lst_GetItemsInXPlane( gCtrlPtList, 
              xVoxl_GetX(inVolumeVox), &list );
    break;
  default:
    bFound = FALSE;
    goto error;
  }
    
  if ( e3DList != x3Lst_tErr_NoErr )
    goto error;

  /* if we got a list... */
  if ( NULL != list ) {

    /* start with a large distance */
    nClosestDistance = gVolumeDimension * gVolumeDimension; 
    bFound = FALSE;

    /* traverse the list */
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&voxel )) 
     != xList_tErr_EndOfList ) {

      if( voxel ) {

  /* get the distance to the clicked voxel... */
  nDistance =
    ((xVoxl_GetX(inVolumeVox) - xVoxl_GetX(voxel)) * 
     (xVoxl_GetX(inVolumeVox) - xVoxl_GetX(voxel))) +
    ((xVoxl_GetY(inVolumeVox) - xVoxl_GetY(voxel)) * 
     (xVoxl_GetY(inVolumeVox) - xVoxl_GetY(voxel))) +
    ((xVoxl_GetZ(inVolumeVox) - xVoxl_GetZ(voxel)) * 
     (xVoxl_GetZ(inVolumeVox) - xVoxl_GetZ(voxel)));
  
  /* if it's less than our max, mark the distance and copy the vox */
  if ( nDistance < nClosestDistance ) {
    nClosestDistance = nDistance;
    closestVoxel = voxel;
    bFound = TRUE;
  }
      }
    }

    if( eList != xList_tErr_EndOfList )
      goto error;
 
   /* if we found a voxel */
    if ( bFound ) {

      /* return it. */
      *outCtrlPt = closestVoxel;
    }
  }

  goto cleanup;

 error:

  if( eList != xList_tErr_NoErr )
    DebugPrint "xList error %d in FindNearestCtrlPt.\n", eList EndDebugPrint;

  if( e3DList != x3Lst_tErr_NoErr )
    DebugPrint "x3Lst error %d in FindNearestCtrlPt.\n", eList EndDebugPrint;

 cleanup:


  return bFound;
}

void AddNearestCtrlPtToSelection ( xVoxelRef        inVolumeVox, 
           mri_tOrientation inPlane ) {

  xVoxelRef theCtrlPt;

  // if we find a nearest control point...
  if ( FindNearestCtrlPt ( inVolumeVox, inPlane, &theCtrlPt ) ) {

    // add this point to the selection
    SelectCtrlPt( theCtrlPt );
  }  
}

void RemoveNearestCtrlPtFromSelection ( xVoxelRef        inVolumeVox, 
          mri_tOrientation inPlane ) {
  
  xVoxelRef theCtrlPt;

  // if we find a nearest control point...
  if ( FindNearestCtrlPt ( inVolumeVox, inPlane, &theCtrlPt ) ) {

    // remove this point from selection
    DeselectCtrlPt( theCtrlPt );
  }
}

void DeselectAllCtrlPts () {

  xList_tErr eList = xList_tErr_NoErr;
  
  eList = xList_Clear( gSelectionList );
  if ( eList != xList_tErr_NoErr )
    DebugPrint "xList error %d in DeselectAllControlPoints: %s\n",
      eList, xList_GetErrorString( eList ) EndDebugPrint;

  redraw ();
}

void SelectCtrlPt ( xVoxelRef iVoxel ) {

  xList_tErr eList = xList_tErr_NoErr;

  eList = xList_PushItem( gSelectionList, iVoxel );
  if ( eList != xList_tErr_NoErr )
    DebugPrint "xList error %d in SelectCtrlPt: %s\n",
      eList, xList_GetErrorString( eList ) EndDebugPrint;

  redraw ();
}

void DeselectCtrlPt ( xVoxelRef iVoxel ) {

  xList_tErr eList = xList_tErr_NoErr;
  xVoxelRef  pVoxel = NULL;

  /* save a ptr to the voxel to remove. */
  pVoxel = iVoxel;

  /* this will return a ptr to the actual voxel removed. don't delte it here,
     because it still exists in the big ctrl pt list */
  eList = xList_RemoveItem( gSelectionList, (void**)&pVoxel );
  if ( eList != xList_tErr_NoErr
       && xList_tErr_ItemNotInList != eList )
    DebugPrint "xList error %d in DeselectCtrlPt: %s\n",
      eList, xList_GetErrorString( eList ) EndDebugPrint;

  redraw ();
}


void NewCtrlPt ( xVoxelRef iCtrlPt,
     tBoolean  ibWriteToFile ) {

  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  ctrlPt  = NULL;

  /* allocate a copy of the voxel */
  xVoxl_New( &ctrlPt );
  xVoxl_Copy( ctrlPt, iCtrlPt );

  // add the voxel to the ctrl pt space
  e3DList = x3Lst_AddItem( gCtrlPtList, ctrlPt, ctrlPt );
  if( e3DList != x3Lst_tErr_NoErr )
    DebugPrint "x3Lst error %d in NewCtrlPt.\n", e3DList EndDebugPrint;

  /* write it to the control point file. */
  if( ibWriteToFile )
    WriteVoxelToControlFile( ctrlPt );
}

void DeleteCtrlPt ( xVoxelRef ipCtrlPt ) {

  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  ctrlPt  = NULL;

  ctrlPt = ipCtrlPt;

  /* remove the item */
  e3DList = x3Lst_RemoveItem( gCtrlPtList, ctrlPt, (void**)&ctrlPt );
  if( e3DList != x3Lst_tErr_NoErr )
    goto error;

  /* delete it */
  xVoxl_Delete( &ctrlPt );

  goto cleanup;
  
 error:
  
  if( x3Lst_tErr_NoErr != e3DList
      && x3Lst_tErr_ItemNotInSpace != e3DList ) {
    DebugPrint "x3Lst error %d in DeleteCtrlPt.\n", e3DList EndDebugPrint;
  }

 cleanup:
  return;
}

void DeleteSelectedCtrlPts () {

  xList_tErr eList   = xList_tErr_NoErr;
  xVoxelRef  ctrlPt  = NULL;

  /* while we have items left to pop from selection list... */
  while( (eList = xList_PopItem( gSelectionList, (void**)&ctrlPt )) 
   != xList_tErr_EndOfList ) {

    /* remove the voxel from the ctrl pt list */
    DeleteCtrlPt( ctrlPt );
  }

  redraw ();
}

// ===========================================================================

// ========================================================= SELECTING REGIONS

void InitSelectionModule () {

  x3Lst_tErr e3DList     = x3Lst_tErr_NoErr;

  e3DList = x3Lst_New( &gSelectedVoxels, 256 );
  if( x3Lst_tErr_NoErr != e3DList ) {
    OutputPrint "Error: Couldn't initialize selection module.\n"
      EndOutputPrint;
    gSelectionList = NULL;
  }
  x3Lst_SetComparator( gSelectedVoxels, CompareVoxels );

  if( NULL != gSelectedVoxels )
    isDisplaySelectedVoxels = TRUE;
}

void DeleteSelectionModule () {

  if ( NULL == gSelectedVoxels )
    return;
  
  x3Lst_Delete( &gSelectedVoxels );
}

void AllowSelectionModuleToRespondToClick (xVoxelRef inScreenVoxel ) {

}

void AddVoxelToSelection ( xVoxelRef iAnaIdx ) {

  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  anaIdx  = NULL;

  if ( NULL == gSelectedVoxels )
    return;

  /* allocate a copy of the voxel */
  xVoxl_New( &anaIdx );
  xVoxl_Copy( anaIdx, iAnaIdx );

  // add the voxel to the ctrl pt space
  e3DList = x3Lst_AddItem( gSelectedVoxels, anaIdx, anaIdx );
  if( e3DList != x3Lst_tErr_NoErr )
    DebugPrint "x3Lst error %d in AddVoxelToSelection.\n", 
      e3DList EndDebugPrint;
}

void RemoveVoxelFromSelection ( xVoxelRef iAnaIdx ) {

  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  where   = NULL;
  xVoxelRef  voxel   = NULL;

  if ( NULL == gSelectedVoxels )
    return;

  xVoxl_New( &where );
  xVoxl_Copy( where, iAnaIdx );

  voxel = iAnaIdx;

  /* remove the item. we'll get back a ptr to the actual voxel we added. */
  e3DList = x3Lst_RemoveItem( gSelectedVoxels, where, (void**)&voxel );
  if( e3DList != x3Lst_tErr_NoErr )
    goto error;

  /* delete it */
  xVoxl_Delete( &voxel );

  goto cleanup;
  
 error:

  if( x3Lst_tErr_NoErr != e3DList
      && x3Lst_tErr_ItemNotInSpace != e3DList ) {
    DebugPrint "x3Lst error %d in RemoveVoxelFromSelection.\n",
      e3DList EndDebugPrint;
  }

 cleanup:

    xVoxl_Delete( &where );
}


void ClearSelection () {

  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  
  if ( NULL == gSelectedVoxels )
    return;

  e3DList = x3Lst_Clear( gSelectedVoxels );
  if ( e3DList != x3Lst_tErr_NoErr )
    DebugPrint "x3Lst error %d in ClearSelection\n",
      e3DList EndDebugPrint;

  MWin_RedrawAll( gMeditWindow );
}

void SaveSelectionToLabelFile ( char * isFileName ) {

  LABEL*        pLabel     = NULL;
  LABEL_VERTEX* pVertex    = NULL;
  Real          rRASX      = 0;
  Real          rRASY      = 0;
  Real          rRASZ      = 0;
  xVoxelRef     anaVoxel   = NULL;
  xListRef      list       = NULL;
  x3Lst_tErr    e3DList    = x3Lst_tErr_NoErr;
  xList_tErr    eList      = xList_tErr_NoErr;
  int           nNumVoxels = 0;
  int           nListCount = 0;
  int           nVoxel     = 0;
  int           nPlane     = 0;
  int           eLabel     = 0;

  /* get the number of selected voxels we have, or the sum of the voxels
     in a particular plane. */
  for ( nPlane = 0; nPlane < gVolumeDimension; nPlane++ ) {

    e3DList = x3Lst_GetItemsInXPlane( gSelectedVoxels, nPlane, &list );
    if( x3Lst_tErr_NoErr != e3DList )
      goto error;

    eList = xList_GetCount( list, &nListCount );
    if( xList_tErr_NoErr != eList )
      goto error;

    nNumVoxels += nListCount;
  }

  /* allocate a label file with that number of voxels, our subject name,
     and the passed in label file name. */
  //  pLabel = LabelAlloc( nNumVoxels, pname, isFileName );
  pLabel = LabelAlloc( nNumVoxels, NULL, isFileName );
  if ( NULL == pLabel ) {
    DebugPrint "\tCouldn't allocate label.\n" EndDebugPrint;
    OutputPrint "ERROR: Couldn't save label.\n" EndOutputPrint;
    goto error;
  }

  /* set the number of points in the label */
  pLabel->n_points = nNumVoxels;

  /* for each plane */
  nVoxel = 0;
  for ( nPlane = 0; nPlane < gVolumeDimension; nPlane++ ) {

    // get the list for this x value.
    e3DList = x3Lst_GetItemsInXPlane( gSelectedVoxels, nPlane, &list );
    if( e3DList != x3Lst_tErr_NoErr )
      goto error;

    /* traverse the list */
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&anaVoxel )) 
     != xList_tErr_EndOfList ) {

      if( anaVoxel ) {

  /* convert voxel to ras */
  VoxelToRAS ( xVoxl_ExpandInt(anaVoxel),
         &rRASX, &rRASY, &rRASZ );

  /* get a ptr the vertex in the label file. */
  pVertex = &(pLabel->lv[nVoxel]);
      
  /* set the vertex */
  pVertex->x = rRASX;
  pVertex->y = rRASY;
  pVertex->z = rRASZ;
      
  /* set the vno to -1, which is significant somewhere outside the realm
     of tkmedit. set stat value to something decent and deleted to not */
  pVertex->vno = -1;
  pVertex->stat = 0;
  pVertex->deleted = FALSE;

  /* inc our global count. */
  nVoxel++;
      }
    }

    if( eList != xList_tErr_EndOfList )
      goto error;
  }

  // write the file
  eLabel = LabelWrite( pLabel, isFileName );
  if ( NO_ERROR != eLabel ) {
    DebugPrint "Error in LabelWrite().\n" EndDebugPrint;
    OutputPrint "ERROR: Couldn't write label to file.\n" EndOutputPrint;
    goto error;
  }

  goto cleanup;

 error:

  if( eList != xList_tErr_NoErr ) {
    DebugPrint "xList error %d in SaveSelectionToLabelFile.\n", 
      eList EndDebugPrint;
    OutputPrint "Error saving lavel file.\n" EndOutputPrint;
  }

  if( e3DList != x3Lst_tErr_NoErr ) {
    DebugPrint "x3Lst error %d in SaveSelectionToLabelFile.\n",
      e3DList EndDebugPrint;
    OutputPrint "Error saving label file.\n" EndOutputPrint;
  }

 cleanup:

  if( pLabel )
    LabelFree( &pLabel );
}

void LoadSelectionFromLabelFile ( char* isFileName ) {

  char          sFileName[256] = "";
  LABEL*        pLabel         = NULL;
  LABEL_VERTEX* pVertex        = NULL;
  int           nVoxX          = 0;
  int           nVoxY          = 0;
  int           nVoxZ          = 0;
  xVoxelRef     voxel          = NULL;
  int           nNumVoxels     = 0;
  int           nVoxel         = 0;

  xVoxl_New( &voxel );

  /* make file name */
  MakeFileName( isFileName, tkm_tFileName_Label, sFileName );

  /* read in the label */
  // pLabel = LabelRead( pname, sFileName );
  pLabel = LabelRead( NULL, sFileName );
  if ( NULL == pLabel ) {
    DebugPrint "\tError reading label file.\n" EndDebugPrint;
    OutputPrint "ERROR: Couldn't read the label.\n" EndOutputPrint;
    goto error;
    
  }

  /* for each vertex in there... */
  nNumVoxels = pLabel->max_points;
  for( nVoxel = 0; nVoxel < nNumVoxels; nVoxel++ ) {
    
    /* get the vertex. */
    pVertex = &(pLabel->lv[nVoxel]);
    
    /* only process verticies that arn't deleted. */
    if ( !(pVertex->deleted) ) {
      
      /* covert from ras to voxel */
      RASToVoxel( pVertex->x, pVertex->y, pVertex->z,
      &nVoxX, &nVoxY, &nVoxZ );
      
      /* add to the selection */
      xVoxl_Set( voxel, nVoxX, nVoxY, nVoxZ );
      AddVoxelToSelection( voxel );
    }
  }

  /* set the window's selection again to force a redraw. */
  MWin_SetSelectionSpace( gMeditWindow, -1, gSelectedVoxels );

  goto cleanup;

 error:

 cleanup:

  if( pLabel )
    LabelFree( &pLabel );

  xVoxl_Delete( &voxel );
}

void GraphSelectedRegion () {

  xVoxelRef  voxel       = NULL;
  xListRef   list        = NULL;
  FunV_tErr  eFunctional = FunV_tErr_NoError;
  x3Lst_tErr e3DList     = x3Lst_tErr_NoErr;
  xList_tErr eList       = xList_tErr_NoErr;
  int        nPlane      = 0;

  // clear the functional display list.
  eFunctional = FunV_BeginSelectionRange( gFunctionalVolume );
  if( FunV_tErr_NoError != eFunctional )
    goto error;

  // get the ctrl pts in the list...
  for ( nPlane = 0; nPlane < gVolumeDimension; nPlane++ ) {

    // get the list for this x value.
    e3DList = x3Lst_GetItemsInXPlane( gSelectedVoxels, nPlane, &list );
    if( e3DList != x3Lst_tErr_NoErr )
      goto error;

    /* traverse the list */
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&voxel )) 
     != xList_tErr_EndOfList ) {

      if( voxel ) {
  
  /* add it to the functional display list. */
  eFunctional = 
    FunV_AddAnatomicalVoxelToSelectionRange( gFunctionalVolume, 
             voxel );
  if( FunV_tErr_NoError != eFunctional )
    goto error;
      }
    }

    if( eList != xList_tErr_EndOfList )
      goto error;
  }

  /* finish the list */
  eFunctional = FunV_EndSelectionRange( gFunctionalVolume );
  if( FunV_tErr_NoError != eFunctional )
    goto error;

  goto cleanup;

 error:

  if( eList != xList_tErr_NoErr ) {
    DebugPrint "xList error %d in GraphSelectedRegion.\n", eList EndDebugPrint;
    OutputPrint "Error graphing selection.\n" EndOutputPrint;
  }

  if( e3DList != x3Lst_tErr_NoErr ) {
    DebugPrint "x3Lst error %d in GraphSelectedRegion.\n", eList EndDebugPrint;
    OutputPrint "Error graphing selection.\n" EndOutputPrint;
  }

  if( eFunctional != FunV_tErr_NoError ) {
    DebugPrint "FunV error %d in GraphSelectedRegion.\n", 
      eFunctional EndDebugPrint;
    OutputPrint "Error graphing selection.\n" EndOutputPrint;
  }


 cleanup:
  return;
}

// ===========================================================================

// ============================================================ EDITING VOXELS

void EditVoxelInRange(xVoxelRef     ipVoxel, 
           tVolumeValue inLow, 
           tVolumeValue inHigh, 
           tVolumeValue inNewValue ) {
  
  tVolumeValue value    = 0;
  tVolumeValue newValue = 0;

  /* get the value at this point. */
  value    = GetVoxelValue( gAnatomicalVolume, xVoxl_ExpandInt(ipVoxel) );
  newValue = value;

  /* if it's in the range... */
  if( value >= inLow && value <= inHigh ) {

    /* set a new value */
    newValue = inNewValue;
  }

  /* if values are different, set and add to undo list. */
  if( value != newValue ) {
    SetVoxelValue( gAnatomicalVolume, xVoxl_ExpandInt(ipVoxel), newValue );
    AddVoxelAndValueToUndoList( ipVoxel, value );
  }

  /* mark this slice as changed. */
  changed[ xVoxl_GetZ(ipVoxel) ] = TRUE;
  editedimage = imnr0 + xVoxl_GetZ(ipVoxel);
}

// ===========================================================================

// ============================================================= FILENAME MGMT

tBoolean SetSubjectHomeDirFromEnv ( char* isSubject ) {

  char* sEnvVar = NULL;

  sEnvVar = getenv( "SUBJECTS_DIR" );
  if( NULL == sEnvVar ) {
    return FALSE;
  }

  /* make the dir */
  sprintf( gsSubjectHomeDir, "%s/%s", sEnvVar, isSubject );
  
  return TRUE;
}

void SetSubjectHomeDir ( char* isHomeDir ) {

  strcpy( gsSubjectHomeDir, isHomeDir );
}

tBoolean FindUserHomeDir () {

  char* sCurrentWorkingDir = NULL;

#ifdef Linux
  sCurrentWorkingDir = getcwd(NULL,0);
#else
  sCurrentWorkingDir = getenv("PWD");
#endif
  
  if( NULL == sCurrentWorkingDir ) {
    return FALSE;
  }
  
  /* save cwd */
  strcpy( gsUserHomeDir, sCurrentWorkingDir );
  
  return TRUE;
}

void MakeFileName ( char*          isInput,
        tkm_tFileName  iType, 
        char*          osCompleteFileName ) {

  char sFilename[256] = "";

  /* look at the first character */
  switch( isInput[0] ) {

  /* tilde, attach the subject home dir. */
  case '~':
    sprintf( sFilename, "%s%s", gsSubjectHomeDir, &(isInput[1]) );
    break;
    
  /* period, use cwd nd then rest of input */
  case '.':
    sprintf( sFilename, "%s%s", gsUserHomeDir, &(isInput[1]) );
    break;

  /* slash, it's a full path. */
  case '/':
    strcpy( sFilename, isInput );
    break;

  /* else, prepend subject home, then sub dir, then rest of file name. */
  default:
    sprintf( sFilename, "%s/%s/%s", gsSubjectHomeDir, 
       ksaFileNameSubDirs[iType], isInput );
    break;
  }

  strcpy( osCompleteFileName, sFilename );
}

void ExtractSubjectName ( char* isDataSource,
        char* osSubjectName ) {

  int   nChar      = 0;
  int   nWordChar  = 0;
  char* sWord      = NULL;
  char  sName[256] = "";
  int   nLastSlash = 0;

  /* look for 'subjects' in the title */
  sWord = strstr( isDataSource, "subjects/" );
  if( NULL != sWord ) {

    /* we're at the s in subjects now. scoot ahead to the slash. */
    nChar = 0;
    while( sWord[nChar] != '/' &&
     sWord[nChar] != '\0' ) {
      nChar++;
    }

    /* if found, use the next part as the name */
    nWordChar = 0;
    nChar++; /* get past the slash */
    while( sWord[nChar] != '/' &&
     sWord[nChar] != '\0' ) {
      sName[nWordChar] = sWord[nChar];
      nWordChar++;
      nChar++;
    }

  } else {
    
    /* else just use the last part */
    nChar = 0;
    while( isDataSource[nChar] != '\0' ) {
      if( isDataSource[nChar] == '/' )
  nLastSlash = nChar;
      nChar++;
    }

    /* if we got it, use it, else use the whole source. */
    if( isDataSource[nLastSlash] == '/' ) 
      strcpy( sName, &(isDataSource[nLastSlash+1]) );
    else 
      strcpy( sName, isDataSource );
  }

  /* return the result */
  strcpy( osSubjectName, sName );
}

void ExtractVolumeName ( char* isDataSource,
       char* osVolumeName ) {

  int nChar      = 0;
  int nLastSlash = 0;

  /* look for the last / */
  while( isDataSource[nChar] != '\0' ) {
    if( isDataSource[nChar] == '/' )
      nLastSlash = nChar;
    nChar++;
  }

  /* if we got one, everything from there+1 into the subjectname, 
     otherwise just use the whole name */
  if( isDataSource[nLastSlash] == '/' )
    strcpy( osVolumeName, &(isDataSource[nLastSlash+1]) );
  else
    strcpy( osVolumeName, isDataSource );
}

// ===========================================================================

/* ============================================= Coordinate transformations */

void InitTransformation () {

  MATRIX* mTemp = NULL;

  /* create the transform object */
  Trns_New( &gRASTransform );

  /* set the a to ras matrix */
  mTemp  = MatrixAlloc( 4, 4, MATRIX_REAL );
  *MATRIX_RELT(mTemp,1,1) = -1.0;
  *MATRIX_RELT(mTemp,2,3) = 1.0;
  *MATRIX_RELT(mTemp,3,2) = -1.0;
  *MATRIX_RELT(mTemp,1,4) = 128.0;
  *MATRIX_RELT(mTemp,2,4) = -128.0;
  *MATRIX_RELT(mTemp,3,4) = 128.0;
  *MATRIX_RELT(mTemp,4,4) = 1.0;
  Trns_CopyAtoRAS( gRASTransform, mTemp );

  MatrixIdentity( 4, mTemp );
  Trns_CopyBtoRAS( gRASTransform, mTemp );
  Trns_CopyARAStoBRAS( gRASTransform, mTemp );
  

  MatrixFree( &mTemp );
}

inline
char IsVoxelInBounds ( int x, int y, int z ) {

  if ( x < 0 || y < 0 || z < 0 ||
       x >= xnum || y >= ynum || z >= xnum ) {
    return FALSE;
  }

  return TRUE;
}

inline
char IsRASPointInBounds ( Real x, Real y, Real z ) {

  if ( x < xx0 || x > xx1 || y < yy0 || y > yy1 || z < zz0 || z > zz1 ) {
    return FALSE;
  }

  return TRUE;
}

xVoxel gCoordA, gCoordB;

void RASToVoxel ( Real x, Real y, Real z,        // incoming ras coords
                  int *xi, int *yi, int *zi ) {  // outgoing voxel coords

  xVoxl_SetFloat( &gCoordA, (float)x, (float)y, (float)z );

  tkm_ConvertRASToVolume( &gCoordA, &gCoordB );

  *xi = xVoxl_GetX( &gCoordB );
  *yi = xVoxl_GetY( &gCoordB );
  *zi = xVoxl_GetZ( &gCoordB );

  if ( ! IsVoxelInBounds ( *xi, *yi, *zi ) ) {
    *xi = *yi = *zi = 0;
  }
}               

void VoxelToRAS ( int xi, int yi, int zi,        // incoming voxel coords
                  Real *x, Real *y, Real *z ) {  // outgoing RAS coords

  if ( ! IsVoxelInBounds ( xi, yi, zi ) ) {
    *x = *y = *z = 0;
    return;
  }

  xVoxl_Set( &gCoordA, xi, yi, zi );

  tkm_ConvertVolumeToRAS( &gCoordA, &gCoordB );

  *x = xVoxl_GetFloatX( &gCoordB );
  *y = xVoxl_GetFloatY( &gCoordB );
  *z = xVoxl_GetFloatZ( &gCoordB );
}

/* ===================================================== General utilities */

void SendUpdateMessageToTKWindow () {

  Tcl_Eval ( GetTclInterp(), "unzoomcoords; sendupdate;" );
}

void SetTclInterp ( Tcl_Interp * inInterp ) {

  gTclInterp = inInterp;
}

Tcl_Interp * GetTclInterp () {

  return gTclInterp;
}

void SendTCLCommand ( char * inCommand ) {
  
  int theErr;
  Tcl_Interp * theInterp;

  // get the interp and send the command.
  theInterp = GetTclInterp ();

  if ( NULL != theInterp ) {

    theErr = Tcl_Eval ( theInterp, inCommand );

    // print any error msgs.
    if ( TCL_OK != theErr ) {
      DebugPrint "Cmd: %s\n", inCommand EndDebugPrint;
      DebugPrint "\tCommand did not return OK.\n" EndDebugPrint;
      DebugPrint "\tResult: %s\n", theInterp->result EndDebugPrint;
    }

  } else {

    /* cache the command so we can send it later */
    if( gNumCachedCommands < kMaxNumCachedCommands ) {
      strcpy( gCachedTclCommands[gNumCachedCommands++], inCommand );
    } else {
      DebugPrint "Couldn't cache %s\n", inCommand EndDebugPrint;
    }
  }
}

void SendCachedTclCommands () {

  int nCommand = 0;

  /* send all our cached commands */
  for( nCommand = 0; nCommand < gNumCachedCommands; nCommand++ ) {
    SendTCLCommand( gCachedTclCommands[nCommand] );
  }
}


// ============================================================= VOLUME ACCESS

void InitVolume ( tVolumeRef* ioVolume, int inDimension ) {

  tVolumeRef theVolume = *ioVolume;

  /* if already exsits, trash it */
  if( NULL != theVolume ) {
    free( theVolume );
  }

  gVolumeDimension = inDimension;
  theVolume = (tVolumeRef) malloc ( gVolumeDimension * 
            gVolumeDimension *
            gVolumeDimension );

  if ( NULL == theVolume ) {
    DebugPrint "InitVolume(): Fatal error, couldn't allocate volume.\n"
      EndDebugPrint;
    exit(1);
  }

  *ioVolume = theVolume;
}

void DeleteVolume( tVolumeRef* ioVolume ) {

  tVolumeRef pVolume = NULL;

  if( NULL == ioVolume ) {
    DebugPrint "DeleteVolume(): ioVolume was NULL!\n" EndDebugPrint;
    goto cleanup;
  }

  /* try to get volume */
  pVolume = *ioVolume;
  if( NULL == pVolume ) {
    DebugPrint "DeleteVolume(): Volume was NULL!\n" EndDebugPrint;
    goto cleanup;
  }

  /* free volume */
  free( pVolume );
  
  /* set outgoing to null */
  *ioVolume = NULL;

 cleanup:

  return;
}

inline
void SetVoxelValue ( tVolumeRef inVolume,
         int x, int y, int z, tVolumeValue inValue ) {

  inVolume [ (gVolumeDimension * gVolumeDimension * z) +
     (gVolumeDimension * y) + x ] = inValue;
}

inline
tVolumeValue GetVoxelValue ( tVolumeRef inVolume,
            int x, int y, int z ) {

  return inVolume [ (gVolumeDimension * gVolumeDimension * z) +
      (gVolumeDimension * y) + x ];
}

tVolumeValue * GetVolumeSlicePtr ( tVolumeRef inVolume, int inSlice ) {

  tVolumeValue * theSlicePtr = NULL;

  theSlicePtr = 
    &( inVolume[ (gVolumeDimension * gVolumeDimension * inSlice) ]);

  return theSlicePtr;
}

void SetVolumeColorScale ( tkm_tVolumeType iVolume, 
         float           ifThreshold, 
         float           ifSquash ) {

  float thresh = 0;
  float squash = 0;
  float scale = 0;
  float i = 0;
  int   nColor = 0;
  char sTclArguments [128];

  gfVolumeColorSquash[iVolume]    = ifSquash;
  gfVolumeColorThreshold[iVolume] = ifThreshold;

  /* generate the look up table */
  scale = knMaxVolumeValue;
  thresh = gfVolumeColorThreshold[iVolume];
  squash = gfVolumeColorSquash[iVolume];
  for( nColor = 0; nColor < knNumVolumeValues; nColor++ ) { 
    i = (float) nColor;
    gfaVolumeColors[iVolume][nColor] =
      (1.0 / (1.0 + exp(-squash * ((i/scale) - thresh)))) * scale + 0.5;
    if( gfaVolumeColors[iVolume][nColor] > knMaxVolumeValue)
      gfaVolumeColors[iVolume][nColor] = knMaxVolumeValue;
      
  }

  /* tell window to redraw */
  if( NULL != gMeditWindow ) {
    MWin_RedrawAll( gMeditWindow );
  }
  
  /* update tcl window */
  sprintf ( sTclArguments, "%d %2f %2f", iVolume,
      gfVolumeColorThreshold[iVolume], gfVolumeColorSquash[iVolume] );
  tkm_SendTclCommand ( tkm_tTclCommand_UpdateVolumeColorScale, 
           sTclArguments );
}

void GetVolumeColor ( tkm_tVolumeType iVolume, 
          tVolumeValue    iucValue, 
          xColor3fRef     oColor ) {

  oColor->mfRed   = (float)gfaVolumeColors[iVolume][iucValue] / 
    (float)knMaxVolumeValue;
  oColor->mfGreen = (float)gfaVolumeColors[iVolume][iucValue] / 
    (float)knMaxVolumeValue;
  oColor->mfBlue  = (float)gfaVolumeColors[iVolume][iucValue] / 
    (float)knMaxVolumeValue;
}

void BuildVolumeMaxIntProj( tVolumeRef iVolume, 
          tVolumeRef* iopMaxIntProjVolume ) {

  tVolumeRef volume = NULL;
  int nVolumeSize = 0;
  int nOrientation = 0;
  int nVolumeIndex = 0;
  int nX = 0;
  int nY = 0;
  int nSlice = 0;
  tVolumeValue value = 0;

  volume = *iopMaxIntProjVolume;

  /* if not allocated, allocate it. */
  nVolumeSize = gVolumeDimension * gVolumeDimension * mri_knNumOrientations *
    sizeof( tVolumeValue );
  if( NULL == volume ) 
    volume = (tVolumeValue*) malloc ( nVolumeSize );

  /* zero the volume */
  bzero( volume, nVolumeSize );

  /* for every orientation */
  for( nOrientation = 0;
       nOrientation < mri_knNumOrientations; 
       nOrientation++ ) {

    /* get dest ptr in the volume */
    nVolumeIndex = nOrientation * (gVolumeDimension * gVolumeDimension);

    /* for every point in volume... */
    for( nY = 0; nY < gVolumeDimension; nY++ ) {
      for( nX = 0; nX < gVolumeDimension; nX++ ) {
    
  /* for every slice... */
  for( nSlice = 0; nSlice < gVolumeDimension; nSlice++ ) {
    
    /* get the anatomical value */
    switch( nOrientation ) {
    case mri_tOrientation_Coronal:
      value = GetVoxelValue( iVolume, nX, nY, nSlice );
      break;

    case mri_tOrientation_Horizontal:
      value = GetVoxelValue( iVolume, nX, nSlice, nY );
      break;
    case mri_tOrientation_Sagittal:
      value = GetVoxelValue( iVolume, nSlice, nY, nX );
      break;
    default:
      break;
    }

    /* if greater than the value in the max int proj, save it */
    if( value > volume[nVolumeIndex] ) 
      volume[nVolumeIndex] = value;
  }

  /* advance the volume ptr */
  nVolumeIndex++;
      }
    }
  }

  /* return the volume */
  *iopMaxIntProjVolume = volume;
}

tVolumeValue GetVolumeMaxIntProjValue( tVolumeRef iMaxIntProjVolume, 
               mri_tOrientation iOrientation,
              xVoxelRef iVoxel ) {

  int nX = 0;
  int nY = 0;

  switch( iOrientation ) {
  case mri_tOrientation_Coronal:
    nX = xVoxl_GetX( iVoxel );
    nY = xVoxl_GetY( iVoxel );
    break;
  case mri_tOrientation_Horizontal:
    nX = xVoxl_GetX( iVoxel );
    nY = xVoxl_GetZ( iVoxel );
    break;
  case mri_tOrientation_Sagittal:
    nX = xVoxl_GetZ( iVoxel );
    nY = xVoxl_GetY( iVoxel );
    break;
  default:
    break;
  }
      
 return iMaxIntProjVolume[ (iOrientation * 
          (gVolumeDimension * gVolumeDimension)) +
       (nY * gVolumeDimension) + nX ];
}

// ========================================================= FUNCTIONAL VOLUME

void LoadFunctionalOverlay( char* isFilename ) {

  char      sFilename[256] = "";
  FunV_tErr eFunctional    = FunV_tErr_NoError;

  MakeFileName( isFilename, tkm_tFileName_Functional, sFilename );
  
  eFunctional = FunV_LoadOverlay( gFunctionalVolume, sFilename );
  if( FunV_tErr_NoError != eFunctional ) {
    OutputPrint "ERROR: Couldn't load functional overlay.\n" EndOutputPrint;
  }
}

void LoadFunctionalTimeCourse( char* isFilename ) {

  char      sFilename[256] = "";
  FunV_tErr eFunctional    = FunV_tErr_NoError;

  MakeFileName( isFilename, tkm_tFileName_Functional, sFilename );
  
  eFunctional = FunV_LoadTimeCourse( gFunctionalVolume, sFilename );
  if( FunV_tErr_NoError != eFunctional ) {
    OutputPrint "ERROR: Couldn't load functional time course.\n" EndOutputPrint;
  }
}


// ================================================= PARCELLATION and ROI GROUP

void LoadParcellationVolume ( char* isVolumeDirWithPrefix,
            char* isColorFileName ) {

  char      sParcellationFileName[256] = "";
  char      sColorFileName[256]        = "";
  CLUT_tErr eColorTable                = CLUT_tErr_NoErr;
  ROIG_tErr eROIGroup                  = ROIG_tErr_NoErr;

  /* make a file name */
  MakeFileName( isVolumeDirWithPrefix, tkm_tFileName_Parcellation, 
    sParcellationFileName );
  MakeFileName( isColorFileName, tkm_tFileName_Parcellation, 
    sColorFileName );

  /* try to load color table */
  eColorTable = CLUT_New( &gColorTable, isColorFileName );
  if( CLUT_tErr_NoErr != eColorTable ) {
    OutputPrint "Error parsing color lookup table.\n" EndOutputPrint;
    goto error;
  }

  /* free existing roi group if present. */
  if( NULL != gROIGroup ) {
    eROIGroup = ROIG_Delete( &gROIGroup );
    if( ROIG_tErr_NoErr != eROIGroup )
      goto error;
  }

  /* read in parcellation volume */
  eROIGroup = ROIG_New( &gROIGroup );
  if( ROIG_tErr_NoErr != eROIGroup )
    goto error;

  eROIGroup = ROIG_ImportVolume( gROIGroup, sParcellationFileName, 
         isColorFileName );
  if( ROIG_tErr_NoErr != eROIGroup )
    goto error;

  /* set parcellation volume in window */
  if( gMeditWindow ) {
    MWin_SetROIGroup( gMeditWindow, -1, gROIGroup );
  }

  goto cleanup;

 error:

  if( ROIG_tErr_NoErr != eROIGroup ) {
    DebugPrint "LoadParcellationVolume: Couldn't open %s\n",
      sParcellationFileName EndDebugPrint;
    OutputPrint "Error finding or opening parcellation data.\n" EndOutputPrint;
  }

  if( CLUT_tErr_NoErr != eColorTable ) {
    DebugPrint "LoadParcellationVolume: Couldn't open %s\n",
      sColorFileName EndDebugPrint;
    OutputPrint "Error finding or opening color table.\n" EndOutputPrint;
  }

 cleanup:
  return;
}

void GetROIColorAtVoxel (xVoxelRef ipVoxel, xColor3fRef oColor ) {

  int         index        = 0;
  CLUT_tErr   eColorTable  = CLUT_tErr_NoErr;

  /* get the index of this voxel */
  GetROILabel( ipVoxel, &index, NULL );

  /* If 0, skip */
  if( 0 == index ) {
    oColor->mfRed   = 0;
    oColor->mfGreen = 0;
    oColor->mfBlue  = 0;
    goto cleanup;
  }

  /* get the color out of the color map */
  eColorTable = CLUT_GetColorFloat( gColorTable, index, oColor );
  if( CLUT_tErr_NoErr != eColorTable )
    goto error;

  goto cleanup;

 error:

  DebugPrint "Error in GetParcellationColor.\n" EndDebugPrint;

  oColor->mfRed   = 0;
  oColor->mfGreen = 0;
  oColor->mfBlue  = 0;

 cleanup:
  return;
}

void GetROILabel ( xVoxelRef ipVoxel, 
       int*      onIndex,
       char*     osLabel ) {

  int            index        = 0;
  CLUT_tErr      eColorTable  = CLUT_tErr_NoErr;
  ROIG_tErr      eROIGroup    = ROIG_tErr_NoErr;
  mriROIVoxelRef voxel        = NULL;

  /* get the voxel at this location */
  eROIGroup = ROIG_GetVoxel( gROIGroup, ipVoxel, &voxel );
  if( NULL != voxel ) {
    index = voxel->mLabelIndex;
  } else {
    index = 0;
  }

  /* get the label out of the table and return it and in the index */
  if( NULL != osLabel ) {
    if( 0 == index ) {
      strcpy( osLabel, "None" );
    } else {
      eColorTable = CLUT_GetLabel( gColorTable, index, osLabel );
      if( CLUT_tErr_NoErr != eColorTable )
  goto error;
    }
  }

  /* return the index */
  if( NULL != onIndex ) {
    *onIndex = index;
  }

  goto cleanup;

 error:

  DebugPrint "Error in GetParcellationColor.\n" EndDebugPrint;

  if( NULL != onIndex )
    *onIndex = -1;
  if( NULL != osLabel )
    strcpy( osLabel, "None" );

 cleanup:
  return;
}

void SelectCurrentROI ( int inIndex ) {

#if 0
  xVoxel voxel;
  tVolumeValue index = 0;

  OutputPrint "Selecting... 0%%" EndOutputPrint;

  /* go thru the volume */
  for( xVoxl_SetZ( &voxel, 0 );
       xVoxl_GetZ( &voxel ) < gVolumeDimension;
       xVoxl_SetZ( &voxel, xVoxl_GetZ( &voxel ) + 1 ) ) {
    for( xVoxl_SetY( &voxel, 0 );
   xVoxl_GetY( &voxel ) < gVolumeDimension;
   xVoxl_SetY( &voxel, xVoxl_GetY( &voxel ) + 1 ) ) {
      for( xVoxl_SetX( &voxel, 0 );
     xVoxl_GetX( &voxel ) < gVolumeDimension;
     xVoxl_SetX( &voxel, xVoxl_GetX( &voxel ) + 1 ) ) {

  /* get the index. if it's the same, add this voxel to the selection */
  index = GetVoxelValue( gParcellationVolume, xVoxl_ExpandInt(&voxel) );
  if( index == inIndex )
    AddVoxelToSelection( &voxel );
      }
    }

    OutputPrint "\rSelecting... %.2f%%       ",
      ((float)xVoxl_GetZ( &voxel ) / (float)gVolumeDimension) * 100.0
      EndOutputPrint;
  }

  OutputPrint "\rSelecting... done!          \n" EndOutputPrint;

  /* redraw */
  UpdateAndRedraw();
#endif
}

// ============================================================== EDITING UNDO

void InitUndoList () {

  xUndL_tErr theErr;

  // new our list.
  theErr = xUndL_New ( &gUndoList, 
           &UndoActionWrapper, &DeleteUndoEntryWrapper );
  if ( xUndL_tErr_NoErr != theErr ) {
    DebugPrint "InitUndoList(): Error in xUndL_New %d: %s\n",
      theErr, xUndL_GetErrorString ( theErr ) EndDebugPrint;
    gUndoList = NULL;
  }

  // set the print function so we can print if necessary.
  theErr = xUndL_SetPrintFunction ( gUndoList, &PrintEntryWrapper );
  if ( xUndL_tErr_NoErr != theErr ) {
    DebugPrint "InitUndoList(): Error in xUndL_SetPrintFunction %d: %s\n",
      theErr, xUndL_GetErrorString ( theErr ) EndDebugPrint;
    gUndoList = NULL;
  }
}

void DeleteUndoList () {

  xUndL_tErr theErr;

  // delete the list.
  theErr = xUndL_Delete ( &gUndoList );
  if ( xUndL_tErr_NoErr != theErr ) {
    DebugPrint "DeleteUndoList(): Error in xUndL_Delete %d: %s\n",
      theErr, xUndL_GetErrorString ( theErr ) EndDebugPrint;
  }
}

void NewUndoEntry ( UndoEntryRef* outEntry,
       xVoxelRef inVoxel, tVolumeValue inValue ) {

  UndoEntryRef this = NULL;
  
  // assume failure.
  *outEntry = NULL;

  // allocate the entry.
  this = (UndoEntryRef) malloc ( sizeof(UndoEntry) );
  if ( NULL == this ) {
    DebugPrint "NewUndoEntry(): Error allocating entry.\n" EndDebugPrint;
    return;
  }

  // copy the voxel in.
  xVoxl_New ( &(this->mVoxel) );
  xVoxl_Copy ( this->mVoxel, inVoxel );

  // copy the value in.
  this->mValue = inValue;

  *outEntry = this;
}

void DeleteUndoEntry ( UndoEntryRef* ioEntry ) {

  UndoEntryRef this = NULL;

  this = *ioEntry;
  if ( NULL == this ) {
    DebugPrint "DeleteUndoEntry(): Got NULL entry.\n" EndDebugPrint;
    return;
  }
  
  // delete the voxel.
  xVoxl_Delete ( &(this->mVoxel) );
  
  // delete the entry.
  free ( this );

  *ioEntry = NULL;
}

void PrintUndoEntry ( UndoEntryRef this ) {

  if ( NULL == this ) {
    DebugPrint "PrintUndoEntry() : INVALID ENTRY\n" EndDebugPrint;
    return;
  }

  DebugPrint "%p voxel (%d,%d,%d)  value = %d\n", this,
    xVoxl_ExpandInt(this->mVoxel), this->mValue EndDebugPrint;

}

void ClearUndoList () {

  xUndL_tErr theErr;

  // clear the list.
  theErr = xUndL_Clear ( gUndoList );
  if ( xUndL_tErr_NoErr != theErr ) {
    DebugPrint "ClearUndoList(): Error in xUndL_Clear %d: %s\n",
      theErr, xUndL_GetErrorString ( theErr ) EndDebugPrint;
  }
}

void AddVoxelAndValueToUndoList (xVoxelRef inVoxel, int inValue ) {

  UndoEntryRef theEntry = NULL;
  xUndL_tErr theErr;

  // make the entry.
  NewUndoEntry ( &theEntry, inVoxel, (tVolumeValue)inValue );
  if ( NULL == theEntry ) {
    DebugPrint "AddVoxelAndValueToUndoList(): Couldn't create entry.\n"
      EndDebugPrint;
    return;
  }

  // add the entry.
  theErr = xUndL_AddEntry ( gUndoList, theEntry );
  if ( xUndL_tErr_NoErr != theErr ) {
   DebugPrint "AddVoxelAndValueToUndoList(): Error in xUndL_AddEntry %d: %s\n",
      theErr, xUndL_GetErrorString ( theErr ) EndDebugPrint;
  }
}

void RestoreUndoList () {
  
  xUndL_tErr theErr;

  // restore the list.
  theErr = xUndL_Restore ( gUndoList );
  if ( xUndL_tErr_NoErr != theErr ) {
    DebugPrint "RestoreUndoList(): Error in xUndL_Restore %d: %s\n",
      theErr, xUndL_GetErrorString ( theErr ) EndDebugPrint;
  }

  /* force a redraw in the window */
  MWin_RedrawAll( gMeditWindow );
}

void DeleteUndoEntryWrapper ( xUndL_tEntryPtr* inEntryToDelete ) {

  if ( NULL == *inEntryToDelete ) {
    DebugPrint "DeleteUndoEntryWrapper(): Got null entry.\n" EndDebugPrint;
    return;
  }
    
  // just typecast and call.
  DeleteUndoEntry ( (UndoEntryRef*)inEntryToDelete );
  *inEntryToDelete = NULL;
}

void UndoActionWrapper      ( xUndL_tEntryPtr  inUndoneEntry, 
            xUndL_tEntryPtr* outNewEntry ) {

  UndoEntryRef theEntryToUndo, theUndoneEntry;
  tVolumeValue theVoxelValue;

                                   /* we're getting an undo entry that should
              be undone or restored. we also want to
              create a new entry and pass it back, so
              we can undo the undo. */

  // get the entry and check it.
  theEntryToUndo = (UndoEntryRef) inUndoneEntry;
  if ( NULL == theEntryToUndo ) {
    DebugPrint "UndoActionWrapper(): Got null entry.\n" EndDebugPrint;
    return;
  }

  // get the value at this voxel.
  theVoxelValue = GetVoxelValue ( gAnatomicalVolume,
          xVoxl_ExpandInt(theEntryToUndo->mVoxel) );

  // create an entry for it.
  NewUndoEntry ( &theUndoneEntry, theEntryToUndo->mVoxel, theVoxelValue );

  // set the voxel value.
  SetVoxelValue ( gAnatomicalVolume,
      xVoxl_ExpandInt(theEntryToUndo->mVoxel),
      theEntryToUndo->mValue );

  // pass back the new entry.
  *outNewEntry = theUndoneEntry;
}

void PrintEntryWrapper ( xUndL_tEntryPtr inEntry ) {

  PrintUndoEntry ( (UndoEntryRef) inEntry );
}

/* ============================================================ HEAD POINTS */

void LoadHeadPts ( char* isHeadPtsFile, 
       char* isTransformFile ) {
  
  HPtL_tErr eHeadPts            = HPtL_tErr_NoErr;
  char      sHeadPtsFile[256]   = "";
  char      sTransformFile[256] = "";
  char*     spTransformFileArg  = NULL;

  /* make filenames */
  MakeFileName( isHeadPtsFile, tkm_tFileName_HeadPoints, sHeadPtsFile );
  if( NULL != isTransformFile ) {
    MakeFileName( isTransformFile, tkm_tFileName_HeadPoints, sTransformFile );
    spTransformFileArg = sTransformFile;
  } else {
    spTransformFileArg = NULL;
  }

  eHeadPts = HPtL_New( &gHeadPoints,
           sHeadPtsFile,
           spTransformFileArg,
           gRASTransform );
  if ( HPtL_tErr_NoErr != eHeadPts ) {
    OutputPrint "ERROR: Couldn't load head points file.\n" EndOutputPrint;
    DebugPrint "HPtL error %d in LoadHeadPts: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) EndDebugPrint;
    goto cleanup;
  }

  MWin_SetHeadPointList( gMeditWindow, -1, gHeadPoints );

  /* show head points stuff */
  tkm_SendTclCommand( tkm_tTclCommand_ShowHeadPointLabel, "1" );

 cleanup:
  return;
}

void RestoreHeadPoints () {

  HPtL_tErr eHeadPts  = HPtL_tErr_NoErr;

  eHeadPts = HPtL_RestoreTransform( gHeadPoints );
  if ( HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint "HPtL error %d in RestoreHeadPoints: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) EndDebugPrint;
    goto cleanup;
  }

  UpdateAndRedraw();

 cleanup:
  return;
}

void WriteHeadPointsTransform () {

  HPtL_tErr eHeadPts  = HPtL_tErr_NoErr;

  eHeadPts = HPtL_WriteTransform( gHeadPoints, NULL );
  if ( HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint "HPtL error %d in WriteHeadPointsTransform: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) EndDebugPrint;
    goto cleanup;
  }

  OutputPrint "Wrote head points transform.\n" EndOutputPrint;

 cleanup:
  return;
}

void WriteHeadPointsFile () {

  HPtL_tErr eHeadPts  = HPtL_tErr_NoErr;

  eHeadPts = HPtL_WriteHeadPointFile( gHeadPoints, NULL );
  if ( HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint "HPtL error %d in WriteHeadPointsFile: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) EndDebugPrint;
    goto cleanup;
  }

  OutputPrint "Wrote head points file.\n" EndOutputPrint;

 cleanup:
  return;
}

void RotateHeadPts ( float ifRadians, 
         char  isAxis ) {

  MATRIX*   mRotation = NULL;
  HPtL_tErr eHeadPts  = HPtL_tErr_NoErr;
  float     fDegrees  = 0;
  mri_tOrientation orientation;

  /* convert to degrees */
  fDegrees = ifRadians * (PI / 180.0);

  /* get the orientation */
  MWin_GetOrientation ( gMeditWindow, &orientation );

  /* create the proper rotate matrix for this orientation */
  switch( isAxis ) {
  case 'x':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      mRotation = MatrixAllocRotation( 4, fDegrees, X_ROTATION );
      break;
    case mri_tOrientation_Horizontal:
      mRotation = MatrixAllocRotation( 4, fDegrees, X_ROTATION );
      break;
    case mri_tOrientation_Sagittal:
      mRotation = MatrixAllocRotation( 4, -fDegrees, Y_ROTATION );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'y':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      mRotation = MatrixAllocRotation( 4, fDegrees, Z_ROTATION );
      break;
    case mri_tOrientation_Horizontal:
      mRotation = MatrixAllocRotation( 4, -fDegrees, Y_ROTATION );
      break;
    case mri_tOrientation_Sagittal:
       mRotation = MatrixAllocRotation( 4, fDegrees, Z_ROTATION );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'z':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      mRotation = MatrixAllocRotation( 4, -fDegrees, Y_ROTATION );
      break;
    case mri_tOrientation_Horizontal:
      mRotation = MatrixAllocRotation( 4, -fDegrees, Z_ROTATION );
      break;
    case mri_tOrientation_Sagittal:
      mRotation = MatrixAllocRotation( 4, -fDegrees, X_ROTATION );
      break;
    default: goto cleanup;
      break;
    }
    break;
  default:
    goto cleanup;
  }

  /* apply the transformation */
  eHeadPts = HPtL_ApplyTransform( gHeadPoints, mRotation );
  if ( HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint "HPtL error %d in RotateHeadPts: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) EndDebugPrint;
    goto cleanup;
  }

  UpdateAndRedraw();

 cleanup:

  if( NULL != mRotation ) 
    MatrixFree( &mRotation );
}

void TranslateHeadPts ( float ifDistance, 
      char  isDirection ) {

  MATRIX*   mTranslation = NULL;
  HPtL_tErr eHeadPts     = HPtL_tErr_NoErr;
  mri_tOrientation orientation;

  mTranslation = MatrixIdentity( 4, NULL );

  /* get the orientation */
  MWin_GetOrientation ( gMeditWindow, &orientation );

  /* create the proper translation matrix for this orientation */
  switch( isDirection ) {
  case 'x':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      *MATRIX_RELT(mTranslation,1,4) = ifDistance; /* ana x */
      break;
    case mri_tOrientation_Horizontal:
      *MATRIX_RELT(mTranslation,1,4) = ifDistance; /* ana x */
      break;
    case mri_tOrientation_Sagittal:
      *MATRIX_RELT(mTranslation,2,4) = -ifDistance; /* ana y */
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'y':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      *MATRIX_RELT(mTranslation,3,4) = ifDistance; /* ana z */
      break;
    case mri_tOrientation_Horizontal:
      *MATRIX_RELT(mTranslation,2,4) = ifDistance; /* ana y */
      break;
    case mri_tOrientation_Sagittal:
      *MATRIX_RELT(mTranslation,3,4) = ifDistance; /* ana z */
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'z':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      *MATRIX_RELT(mTranslation,2,4) = -ifDistance; /* ana y */
      break;
    case mri_tOrientation_Horizontal:
      *MATRIX_RELT(mTranslation,3,4) = -ifDistance; /* ana z */
      break;
    case mri_tOrientation_Sagittal:
      *MATRIX_RELT(mTranslation,1,4) = -ifDistance; /* ana x */
      break;
    default: goto cleanup;
      break;
    }
    break;
  default:
    goto cleanup;
  }
  
  /* apply the transform */
  eHeadPts = HPtL_ApplyTransform( gHeadPoints, mTranslation );
  if ( HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint "HPtL error %d in TranslateHeadPts: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) EndDebugPrint;
    goto cleanup;
  }

  UpdateAndRedraw();

 cleanup:

  if( NULL != mTranslation ) 
    MatrixFree( &mTranslation );
}

void SetSelectedHeadPointLabel ( char* isNewLabel ) {

  HPtL_tHeadPointRef pHeadPoint = NULL;
  MWin_tErr          eWindow    = MWin_tErr_NoErr;

  /* check label */
  if( NULL == isNewLabel )
    goto error;

  /* get the point from the window */
  eWindow = MWin_GetSelectedHeadPt( gMeditWindow, &pHeadPoint );
  if( MWin_tErr_NoErr != eWindow )
    goto error;

  /* if we got one... */
  if( NULL != pHeadPoint ) {
    
    /* set the name */
    strcpy( pHeadPoint->msLabel, isNewLabel );
  }

  goto cleanup;

 error:

  if( MWin_tErr_NoErr != eWindow ) {
    DebugPrint "MWin error %d in SetSelectedHeadPointLabel: %s\n",
      eWindow, MWin_GetErrorString( eWindow ) EndDebugPrint;
  }

 cleanup:

  return;
}

void AlignSelectedHeadPointToAnaIdx ( xVoxelRef iAnaIdx ) {

  HPtL_tHeadPointRef pHeadPoint = NULL;
  MWin_tErr          eWindow    = MWin_tErr_NoErr;
  HPtL_tErr          eHeadPts   = HPtL_tErr_NoErr;

  /* get the point from the window */
  eWindow = MWin_GetSelectedHeadPt( gMeditWindow, &pHeadPoint );
  if( MWin_tErr_NoErr != eWindow )
    goto error;

  /* if we got one... */
  if( NULL != pHeadPoint ) {
    
    /* align to it */
    eHeadPts = HPtL_AlignPointToClientVoxel( gHeadPoints, 
               pHeadPoint, iAnaIdx );
    if( HPtL_tErr_NoErr != eHeadPts )
      goto error;

    /* redraw */
    UpdateAndRedraw();
  }

  goto cleanup;

 error:

  if( MWin_tErr_NoErr != eWindow ) {
    DebugPrint "MWin error %d in SetSelectedHeadPointLabel: %s\n",
      eWindow, MWin_GetErrorString( eWindow ) EndDebugPrint;
  }

  if(  HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint "HPtL error %d in SetSelectedHeadPointLabel: %s\n",
      eWindow, HPtL_GetErrorString( eHeadPts ) EndDebugPrint;
  }

 cleanup:

  return;
 
}



/* ============================================================ ACCESS FUNCS */

void tkm_ConvertVolumeToRAS (xVoxelRef inVolumeVox,xVoxelRef outRASVox ) {

  Trns_tErr eTransform = Trns_tErr_NoErr;

  xVoxl_SetFloat( &gCoordA, xVoxl_ExpandFloat(inVolumeVox) );
  eTransform = Trns_ConvertAtoB( gRASTransform, &gCoordA, &gCoordB );
  if( Trns_tErr_NoErr != eTransform )
    goto error;
   xVoxl_SetFloat( outRASVox, xVoxl_ExpandFloat(&gCoordB) );

  goto cleanup;

 error:
  
  if( Trns_tErr_NoErr != eTransform )
    DebugPrint "Error %d in tkm_ConvertVolumeToTal: %s\n",
      eTransform, Trns_GetErrorString( eTransform ) EndDebugPrint;

 cleanup:

  return;
}

void tkm_ConvertRASToVolume (xVoxelRef inRASVox,xVoxelRef outVolumeVox ) {

  Trns_tErr eTransform = Trns_tErr_NoErr;

  xVoxl_SetFloat( &gCoordA, xVoxl_ExpandFloat( inRASVox ) );
  eTransform = Trns_ConvertBtoA( gRASTransform, &gCoordA, &gCoordB );
  if( Trns_tErr_NoErr != eTransform )
    goto error;
  xVoxl_SetFloat( outVolumeVox, xVoxl_ExpandFloat(&gCoordB) );

  goto cleanup;

 error:
  
  if( Trns_tErr_NoErr != eTransform ) {
    DebugPrint "Error %d in tkm_ConvertRASToVolume: %s\n",
      eTransform, Trns_GetErrorString( eTransform ) EndDebugPrint;
  }

 cleanup:

  return;
}

void tkm_ConvertVolumeToTal (xVoxelRef inVolumeVox,xVoxelRef outTalVox ) {

  Real theTalX, theTalY, theTalZ;
  Real theRASX, theRASY, theRASZ;

  if ( gbAnaToTalTransformLoaded ) {

    xVoxl_SetFloat( &gCoordA, 
        xVoxl_GetX( inVolumeVox ),
        xVoxl_GetY( inVolumeVox ),
        xVoxl_GetZ( inVolumeVox ) );
    
    tkm_ConvertVolumeToRAS( &gCoordA, &gCoordB );
    
    theRASX = xVoxl_GetFloatX( &gCoordB );
    theRASY = xVoxl_GetFloatY( &gCoordB );
    theRASZ = xVoxl_GetFloatZ( &gCoordB );

    transform_point ( get_linear_transform_ptr( &gAnaToTalTransform ),
          theRASX, theRASY, theRASZ,
                      &theTalX, &theTalY, &theTalZ );

    xVoxl_SetFloat ( outTalVox, 
         (float)theTalX, (float)theTalY, (float)theTalZ );
  } else {
    xVoxl_Set ( outTalVox, 0, 0, 0 );
  }

  goto cleanup;

 cleanup:

  return;
}

void tkm_ConvertTalToVolume ( xVoxelRef iTalVox,
            xVoxelRef oAnaIdx ) {

  Real      theRASX    = 0;
  Real      theRASY    = 0;
  Real      theRASZ    = 0;

  if ( gbAnaToTalTransformLoaded ) {

    transform_point ( get_inverse_linear_transform_ptr( &gAnaToTalTransform ), 
          xVoxl_ExpandFloat( iTalVox ),
                      &theRASX, &theRASY, &theRASZ );

    xVoxl_SetFloat( &gCoordA, 
        theRASX, theRASY, theRASZ );
    
    tkm_ConvertRASToVolume( &gCoordA, &gCoordB );
    
    xVoxl_Copy( oAnaIdx, &gCoordB );
         
  } else {
    xVoxl_Set ( oAnaIdx, 0, 0, 0 );
  }

  goto cleanup;

 cleanup:

  return;
}

tBoolean tkm_IsValidVolumeIdx ( xVoxelRef iAnaIdx ) {

  return IsVoxelInBounds( xVoxl_ExpandInt( iAnaIdx ) );
}
  
tBoolean tkm_IsValidRAS ( xVoxelRef iRAS ) {

  return IsRASPointInBounds( xVoxl_ExpandFloat( iRAS ) );
}

tVolumeValue tkm_GetVolumeValue ( tVolumeRef inVolume,
          xVoxelRef inVoxel ) {

  return GetVoxelValue ( inVolume, xVoxl_ExpandInt(inVoxel) );
}

void tkm_GetAnatomicalVolumeColor( tVolumeRef   iVolume,
           tVolumeValue inValue, 
           xColor3fRef  oColor ) {
  
  if( iVolume == gAnatomicalVolume )
    GetVolumeColor( tkm_tVolumeType_Main, inValue, oColor );
  else
    GetVolumeColor( tkm_tVolumeType_Aux, inValue, oColor );
}

tVolumeValue tkm_GetMaxIntProjValue( tVolumeRef iVolume, 
             mri_tOrientation iOrientation, 
            xVoxelRef iVoxel ) {
  if( iVolume == gAnatomicalVolume ) {
    return GetVolumeMaxIntProjValue( gAnatomicalMaxIntProj,
             iOrientation, iVoxel );
  } else if( iVolume == gAuxAnatomicalVolume ) {
    return GetVolumeMaxIntProjValue( gAuxAnatomicalMaxIntProj,
             iOrientation, iVoxel );
  }
  return 0;
}

void tkm_AddNearestCtrlPtToSelection (xVoxelRef inVolumeVox, 
               mri_tOrientation inPlane ) {

  AddNearestCtrlPtToSelection ( inVolumeVox, inPlane );
}

void tkm_RemoveNearestCtrlPtFromSelection (xVoxelRef inVolumeVox, 
              mri_tOrientation inPlane ) {

  RemoveNearestCtrlPtFromSelection ( inVolumeVox, inPlane );
}

void tkm_NewCtrlPt () {

 xVoxelRef theCursor = NULL;
  xVoxl_New ( &theCursor );
  MWin_GetCursor ( gMeditWindow, theCursor );
  NewCtrlPt ( theCursor, TRUE );
  xVoxl_Delete ( &theCursor );
}

void tkm_DeselectAllCtrlPts () {

  DeselectAllCtrlPts ();
}

void tkm_DeleteSelectedCtrlPts () {

  DeleteSelectedCtrlPts ();
}

void tkm_WriteControlFile () {

  WriteCtrlPtFile( );
}

void tkm_EditVoxelInRange(xVoxelRef     inVolumeVox, 
         tVolumeValue inLow, 
         tVolumeValue inHigh, 
         tVolumeValue inNewValue ) {

  EditVoxelInRange( inVolumeVox, inLow, inHigh, inNewValue );

}

void tkm_SelectVoxel (xVoxelRef inVolumeVox ) {

  AddVoxelToSelection ( inVolumeVox );
}

void tkm_DeselectVoxel ( xVoxelRef inVolumeVox ) {

  RemoveVoxelFromSelection( inVolumeVox );
}

void tkm_ClearSelection () {

  ClearSelection ();
}

void tkm_ClearUndoList () {

  ClearUndoList ();
}

void tkm_RestoreUndoList () {

  RestoreUndoList ();
}

void tkm_GetHeadPoint ( xVoxelRef           iAnaIdx,
      mri_tOrientation    iOrientation,
      HPtL_tHeadPointRef* opPoint ) {

  HPtL_tErr eHeadPts = HPtL_tErr_NoErr;
  HPtL_tHeadPointRef pPoint = NULL;
  HPtL_tIterationPlane plane = HPtL_tIterationPlane_X;

  if( gHeadPoints ) {
    
    switch( iOrientation ) {
    case mri_tOrientation_Coronal:
      plane = HPtL_tIterationPlane_Z;
      break;
    case mri_tOrientation_Horizontal:
      plane = HPtL_tIterationPlane_Y;
      break;
    case mri_tOrientation_Sagittal:
      plane = HPtL_tIterationPlane_X;
      break;
    default:
      goto error;
      break;
    }

    eHeadPts = HPtL_FindNearestPoint( gHeadPoints, plane,
              1.0, iAnaIdx, &pPoint );
    if( HPtL_tErr_NoErr != eHeadPts )
      goto error;
    
    if( NULL == pPoint ) 
      goto error;

    *opPoint = pPoint;

  }

  goto cleanup;
  
 error:
  
  *opPoint = NULL;
  
 cleanup:
  return;
}

void tkm_WriteVoxelToControlFile ( xVoxelRef inVolumeVox ) {

  WriteVoxelToControlFile ( inVolumeVox );
}

void tkm_WriteVoxelToEditFile ( xVoxelRef inVolumeVox ) {

  WriteVoxelToEditFile ( inVolumeVox );
}

void tkm_ReadCursorFromEditFile () {

  ReadCursorFromEditFile ();
}

char* tkm_GetSubjectName() {

  return gSubjectName;
}

char* tkm_GetVolumeName() {

  return gSubjectVolumeName;
}

char* tkm_GetAuxVolumeName() {

  return gAuxSubjectVolumeName;
}


void tkm_GetROIColorAtVoxel ( xVoxelRef inVoxel, xColor3fRef oColor ) {
  
  GetROIColorAtVoxel( inVoxel, oColor );
}

void tkm_GetROILabel ( xVoxelRef inVoxel, 
           int* onIndex, char* osLabel ) {

  GetROILabel( inVoxel, onIndex, osLabel );
}

void tkm_SelectCurrentROI ( int inIndex ) {
  
  SelectCurrentROI( inIndex );
}

char kTclCommands [tkm_knNumTclCommands][256] = {

  "UpdateLinkedCursorFlag",
  "UpdateVolumeCursor",
  "UpdateRASCursor",
  "UpdateTalCursor",
  "UpdateVolumeName",
  "UpdateVolumeValue",
  "UpdateAuxVolumeName",
  "UpdateAuxVolumeValue",
  "UpdateFunctionalCoords",
  "UpdateFunctionalValue",
  "UpdateZoomLevel",
  "UpdateOrientation",
  "UpdateDisplayFlag",
  "UpdateTool",
  "UpdateBrush",
  "UpdateBrushThreshold",
  "UpdateVolumeColorScaleInfo",
  "UpdateROILabel",
  "UpdateHeadPointLabel",
  
  /* display status */
  "ShowVolumeCoords",
  "ShowRASCoords",
  "ShowTalCoords",
  "ShowAuxValue",
  "ShowROILabel",
  "ShowHeadPointLabel",
  "ShowFuncCoords",
  "ShowFuncValue",
  "tkm_SetMenuItemGroupStatus tMenuGroup_OverlayOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_TimeCourseOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_SurfaceLoading",
  "tkm_SetMenuItemGroupStatus tMenuGroup_SurfaceViewing",
  "tkm_SetMenuItemGroupStatus tMenuGroup_OriginalSurfaceViewing",
  "tkm_SetMenuItemGroupStatus tMenuGroup_CanonicalSurfaceViewing",
  "tkm_SetMenuItemGroupStatus tMenuGroup_EditHeadPointLabel",

  /* interface configuration */
  "raise .; wm geometry .",
  "CsurfInterface",
  "ErrorDlog",
  "AlertDlog"
};

void tkm_SendTclCommand ( tkm_tTclCommand inCommand,
        char* inArguments ) {

  char theCommand[256];

  if ( inCommand < 0
       || inCommand >= tkm_knNumTclCommands )
    return;

  sprintf ( theCommand, "%s %s", kTclCommands[inCommand], inArguments );
  //  DebugPrint "[] %s\n", theCommand EndDebugPrint;
  SendTCLCommand ( theCommand );
}


xList_tCompare CompareVoxels ( void* inVoxelA, void* inVoxelB ) {

  xList_tCompare eResult = xList_tCompare_GreaterThan;
  xVoxelRef      voxelA  = NULL;
  xVoxelRef      voxelB  = NULL;

  voxelA = (xVoxelRef) inVoxelA;
  voxelB = (xVoxelRef) inVoxelB;

  if( xVoxl_IsEqualInt( voxelA, voxelB ) ) {
    eResult = xList_tCompare_Match;
  } else {
    eResult = xList_tCompare_GreaterThan;
  }

  return eResult;
}
