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
#include <tcl8.3.h>
#include <tclDecls.h>
#include <tk8.3.h>
// #include <tix.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "MRIio.h"
#include "volume_io.h"
#include "rgb_image.h"
#include "fio.h"
#include "mrisurf.h"
#include "mri_conform.h"

#ifndef OPENGL
#define OPENGL
#endif

#define TCL

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

static MRI_SURFACE *mris ;
int white_lolim = 80;
int white_hilim = 140;
int gray_hilim = 100;
int flossflag = TRUE;
int spackleflag = TRUE;
int lim3=170,lim2=145,lim1=95,lim0=75;
double ffrac3=1.0,ffrac2=1.0,ffrac1=1.0,ffrac0=1.0;

char *subjectsdir;   /* SUBJECTS_DIR */
char *srname;        /* sessiondir--from cwd */
char *pname;         /* name */
char *title_str ;    /* window title */
char *imtype;        /* e.g., T1 */
char *imtype2;        /* e.g., T1 */
char *surface;       /* e.g., rh.smooth */
char *mfname;        /* abs single image stem name */
char *sfname;        /* abs surface name */
char *tfname;        /* (dir!) subjectsdir/name/tmp/ */
char *dipfname;      /* curr: $session/bem/brain3d.dip */
char *decfname;      /* curr: $session/bem/brain3d.dec */
char *hpfname;       /* headpts */
char *htfname;       /* headtrans */
char *sgfname;       /* rgb */
char *fsfname;       /* $session/fs/$hemi.fs */
char *fmfname;       /* $session/fs/$hemi.fm */
char *cfname;        /* $home/surf/hemi.curv */
char *xffname;       /* $home/name/mri/transforms/TALAIRACHg4251_FNAME */
char *rfname;        /* script */

/* Talairach stuff */
General_transform talairach_transform ; /* the next two are from this struct */
Transform         *linear_transform = NULL ;
Transform         *inverse_linear_transform = NULL ;
int               transform_loaded = 0 ;



#include "xDebug.h"
#include "xTypes.h"
#include "tkmVoxel.h"
#include "tkmVoxelList.h"
#include "tkmVoxelSpace.h"
#include "VoxelValueList.h"
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

// ================================================================== MESSAGES

                                   /* put our messages here to make changing
              them easier. */
#define kMsg_InvalidClick "\nInvalid click - out of voxel bounds.\n"         \
      "Although the screen displays areas that are out of voxel bounds in "  \
      "black, you cannot click on them to select them. Please click inside " \
      "voxel bounds, closer to the non-black image.\n\n"


// ===========================================================================

// =========================================================== READING VOLUMES

                                   /* this function takes a path and passes
              it to MRIRead. it then takes the MRI
              struct and grabs all the information
              it needs from it. this is an alternative
              function to read_images and is 
              independent of the SUBJECTS_DIR env
              variable. */
void ReadVolumeWithMRIRead ( char * inFileOrPath );


// ===========================================================================

// ==================================================== COORDINATE CONVERSIONS

#include "mritransform.h"

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

void DrawControlPoints ( char * inBuffer, int inPlane, int inPlaneNum );

                                       /* draw control point. switches on 
                                          the current display style of the 
                                          point. */
void DrawCtrlPt ( char * inBuffer,  // video buffer to draw into
                  int inX, int inY, // x,y location in the buffer
                  long inColor );   // color to draw in

                                   /* returns nearest control point on
              the same plane. returns 1 if found
              and 0 otherwise. */
int FindNearestCtrlPt ( VoxelRef inVolumeVox, tkm_tOrientation inPlane,
      VoxelRef outCtrlPt );

                                   /* add and remove control points from
              the selected ones. */
void AddNearestCtrlPtToSelection ( VoxelRef inVolumeVox, 
           tkm_tOrientation inPlane );
void RemoveNearestCtrlPtFromSelection ( VoxelRef inVolumeVox, 
           tkm_tOrientation inPlane );

                                        /* remove the selected control points 
                                           from the control point space */
void DeleteSelectedCtrlPts ();

                                   /* make the input anatomical voxel a
              control point */
void NewCtrlPt ( VoxelRef inVoxel );

                                   /* convert the cursor to anatomical coords
              and make it a control point. */
void NewCtrlPtFromCursor ();

                                   /* deselect all control points */
void DeselectAllCtrlPts ();

                                       /* reads the control.dat file, 
                                          transforms all pts from RAS space 
                                          to voxel space, and adds them as 
                                          control pts */
void ProcessCtrlPtFile ( char * inDir );

                                       /* writes all control points to the
                                          control.dat file in RAS space */
void WriteCtrlPtFile ( char * inDir );

                                       /* tsia */
void ToggleCtrlPtDisplayStatus ();
void SetCtrlPtDisplayStatus ( char inDisplay );
void ToggleCtrlPtDisplayStyle ();
void SetCtrlPtDisplayStyle ( int inStyle );

                                       /* global storage for ctrl space. */
VoxelSpaceRef gCtrlPtList = NULL;

                                       /* global storage for selected ctrl 
                                          pts. */
VoxelListRef gSelectionList = NULL;

                                       /* flag for displaying ctrl pts. if 
                                          true, draws ctrl pts in draw loop. */
char gIsDisplayCtrlPts;

                                       /* style of control point to draw */
#define kCtrlPtStyle_FilledVoxel              1
#define kCtrlPtStyle_Crosshair                2
char gCtrlPtDrawStyle;

                                       /* whether or not we've added the 
                                          contents of the control.dat file to
                                          our control pt space. we only want
                                          to add it once. */
char gParsedCtrlPtFile;

// ===========================================================================

// ================================================================== SURFACES

                                      /* for different surface types */

#define kSurfaceType_Current              0
#define kSurfaceType_Original             1
#define kSurfaceType_Canonical            2

                                       /* flags for determining whether a
                                          surface is loaded. */
char gIsCurrentSurfaceLoaded = FALSE;
char gIsOriginalSurfaceLoaded = FALSE;
char gIsCanonicalSurfaceLoaded = FALSE;

// ===========================================================================

// ========================================================= SELECTING REGIONS

/* selecting regions works much like editing. the user chooses a brush shape and size and paints in a region. the tool can be toggled between selecting and unselecting. the selected pixels are kept in a voxel space for optimized retreival in the draw loop. there are the usual functions for adding and removing voxels as well as saving them out to a file. */

VoxelSpaceRef gSelectedVoxels;
char isDisplaySelectedVoxels;

void InitSelectionModule ();
void DeleteSelectionModule ();

/* grabs the list of voxels selected and draws them into the buffer. */
void DrawSelectedVoxels ( char * inBuffer, int inPlane, int inPlaneNum );

/* handles clicks. uses the current brush settings to paint or paint selected voxels. */
void AllowSelectionModuleToRespondToClick ( VoxelRef inScreenVoxel );

  /* adds or removes voxels to selections. if a voxel that isn't in the 
     selection is told to be removed, no errors occur. this is called from the 
     brush function. */
void AddVoxelToSelection ( VoxelRef inVoxel );
void RemoveVoxelFromSelection ( VoxelRef inVoxel );

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

void EditVoxelInRange( VoxelRef ipVoxel, 
           tVolumeValue inLow, tVolumeValue inHigh, 
           tVolumeValue inNewValue );

// ===========================================================================

// ============================================================= VOLUME ACCESS

#define knNumVolumeValues 256
#define kfDefaultVolumeThreshold 0.35
#define kfDefaultVolumeSquash    12.0

static int gVolumeDimension;

static tVolumeRef gAnatomicalVolume    = NULL;
static tVolumeRef gAuxAnatomicalVolume = NULL;
static tVolumeRef gSnapshotVolume      = NULL;

static float gfVolumeColorThreshold = kfDefaultVolumeThreshold;
static float gfVolumeColorSquash    = kfDefaultVolumeSquash;
static float gfaVolumeColors [knNumVolumeValues];

void InitVolume ( tVolumeRef* ioVolume, int inDimension );
void DeleteVolume ( tVolumeRef* ioVolume );
inline tVolumeValue GetVoxelValue ( tVolumeRef inVolume,
             int x, int y, int z );

inline void SetVoxelValue ( tVolumeRef inVolume,
          int x, int y, int z, tVolumeValue inValue );
tVolumeValue * GetVolumeSlicePtr ( tVolumeRef inVolume, int inSlice );

void SetVolumeColorScale ( float ifThreshold, float ifSquash );

void GetVolumeColor ( tVolumeValue iucValue,
          unsigned char *outRed, 
          unsigned char *outGreen,
          unsigned char *outBlue );

void SnapshotVolume ();
void RestoreVolumeFromSnapshot ();

// ===========================================================================

// ========================================================= FUNCTIONAL VOLUME

#include "tkmFunctionalVolume.h"

tkmFunctionalVolumeRef gFunctionalVolume = NULL;

// ===========================================================================

// ============================================================== PARCELLATION

typedef struct {
  int mRed, mGreen, mBlue;
} tColorEntry;

static tVolumeRef gParcellationVolume   = NULL;
static tColorEntry* gParcellationColors = NULL;
static int gNumParcellationColors;

void LoadParcellationVolume ( char* inVolumeDirWithPrefix,
            char* inColorFileName );

void GetParcellationColor ( VoxelRef inVoxel,
          unsigned char *outRed, 
          unsigned char *outGreen,
          unsigned char *outBlue );

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
void AddVoxelAndValueToUndoList ( VoxelRef inVoxel, int inValue );
void RestoreUndoList ();

                                   /* we need a struct for the undo list. this
              is what we add to it and what we get
              back when the list is restored. */
typedef struct {
  VoxelRef mVoxel;
  tVolumeValue mValue;
} UndoEntry, *UndoEntryRef;

void NewUndoEntry           ( UndoEntryRef* outEntry, 
            VoxelRef inVoxel, tVolumeValue inValue );
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

/* ========================================================== MEDIT WINDOW */

tkmMeditWindowRef gMeditWindow = NULL;
int gDoRedraw = 1;

void flush_redraws ();

/* ======================================================================= */


// ====================================================================== MISC

                                   /* determines if a number is odd. */
#define isOdd(x) (x%2)

                                   /* set and get the tcl interp to send
            the msg to */
void SetTclInterp ( Tcl_Interp * inInterp );
Tcl_Interp * GetTclInterp ();

                                   /* send a tcl command */
void SendTCLCommand ( char * inCommand );


#define kMaxNumCachedCommands 50
char gCachedTclCommands[kMaxNumCachedCommands][256];
int gNumCachedCommands = 0;
void SendCachedTclCommands ();
                                   /* the tcl interpreter */
Tcl_Interp * gTclInterp = NULL;

                                   /* don't start accepting tcl commands
              until setup is complete */
static tBoolean gbAcceptingTclCommands = FALSE;

static tBoolean gbUseCsurfInterface = FALSE;

#define set3fv(v,x,y,z) {v[0]=x; v[1]=y; v[2]=z;}

// ===========================================================================


/*--------------------- prototypes ------------------------------*/
#ifdef Linux
extern void scale2x(int, int, unsigned char *);
#endif


int Medit(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[]);
int read_second_images(char *imdir2) ;
void WriteVoxelToControlFile ( char* inPathName, VoxelRef inVolumeVox );
void WriteVoxelToEditFile ( char* inPathName, VoxelRef inVolumeVox );
void mark_file_vertices(char *fname) ; // read in a file, mark vert indicies
void unmark_vertices(void) ; // set marked field of all vertices to 0
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
int read_binary_surface(char *fname) ;
int read_surface(char *fname) ;
int read_orig_vertex_positions(char *name) ;
int read_canonical_vertex_positions(char *fname) ;
int read_binary_surf(char *fname) ;
int dump_vertex(int vno) ; // prints info about a vertex
int write_images(char *fpref) ; // saves volume to COR files
int read_images(char *fpref) ; // reads volume from COR files
void write_dipoles(char *fname) ;
void write_decimation(char *fname) ;
void read_hpts(char *fname) ;
void read_htrans(char *fname) ;
void write_htrans(char *fname) ;
void make_filenames(char *lsubjectsdir,char *lsrname, char *lpname, 
                    char *limdir, char *lsurface) ;
void mirror(void) ;
void read_fieldsign(char *fname) ;
void read_fsmask(char *fname) ;
void read_binary_curvature(char *fname) ;
void smooth_3d(int niter) ;
void flip_corview_xyz(char *newx, char *newy, char *newz) ;
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

  int i,j;
  char fname[NAME_LENGTH];
  char *lsubjectsdir;
  FILE *fp;
  char lsrname[NAME_LENGTH];
  char *getenv();
  char *cwd; 
  char *word;
  char path[MAX_DIR_DEPTH][NAME_LENGTH];
  FunV_tErr eFunctional = FunV_tErr_NoError;
    
  int  nCurrentArg = 0;
  int  bFatalError = FALSE;
  char sArg[128];
  
  int bSubjectDeclared = FALSE;
  int bUsingMRIRead = FALSE;
  char sSubject[128];
  char sImageDir[128];
  int bSurfaceDeclared = FALSE;
  char sSurface[128];
  int bLocalImageDir = FALSE;
  int bNoEdit = FALSE;
  int tclscriptflag = 0;
  int bLoadingOverlay = FALSE;
  char sOverlayPath[128];
  char sOverlayStem[128];
  int bLoadingTimeCourse = FALSE;
  char sTimeCoursePath[128];
  char sTimeCourseStem[128];
  char bThresh = FALSE;
  FunV_tFunctionalValue min = 0;
  char bMid = FALSE;
  FunV_tFunctionalValue mid = 0;
  char bSlope = FALSE;
  FunV_tFunctionalValue slope = 0;
  char bRevPhaseFlag = FALSE;
  int nRevPhaseFlag = 0;
  char bTruncPhaseFlag = FALSE;
  int nTruncPhaseFlag = 0;
  char sPathAndStem[256] = "";

  /* first get the functional threshold so we don't overwrite the defaults */
  eFunctional = FunV_GetThreshold( gFunctionalVolume, &min, &min, &slope );
  if( FunV_tErr_NoError != eFunctional ) {
    DebugPrint "Medit(): Couldn't get functional volume threshold.\n" 
      EndDebugPrint;
  }

    if (argc<2) {
printf("\n");
printf("Usage:\n");
printf("  Using SUBJECTS_DIR environment variable and subject name as a volume source:\n" );
printf("       %s name imagedir [relative_surf_name] [-] [-o functional_path stem] [-tcl script] [-s scale_factor]\n",argv[0]);
 printf("  Options:\n" );     
printf("           name  subjectdir name              (relative to $SUBJECTS_DIR)\n");
printf("       imagedir  orig,T1,brain,wm,filled      (relative to $subject/mri)\n");
printf("       surfname  ?h.orig,?h.smoothwm,?h.plump (relative to $subject/surf)\n");
 printf("functional_path  full path to functional data\n");
 printf("           stem  stem of functional data image files\n");

printf("\n  Specifying a specific dir or file as a volume source:\n");
printf("       %s -f file_or_path [full_surf_path] [-] [-o functinal_path stem] [-tcl script] [-s scale_factor]\n",argv[0]);
 printf("  Options:\n" );     
printf("-f file_or_path  look at file_or_path and guess format to read in\n");
printf(" full_surf_path  the full path including file name to a surface file\n");

 printf("\n  Using volume data in the current directory:\n" );
printf("       %s local [-] [-o functional_path stem] [-tcl script] [-s scale_factor]\n",argv[0]);

 printf("\n  General options:\n" );
 printf("   -tcl script  run in script mode with no visual output\n");
printf("      lone dash  disable editing\n");
printf("-s scale_factor  set window scale factor. 2 is default for 512x512.\n");
printf("\n                                         [vers: 272134--OpenGL]\n");
printf("\n");
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

    /* don't really know what this does. just set the flag and
       skip the next param */
    tclscriptflag = TRUE;
    nCurrentArg += 2;

  } else if( MATCH( sArg, "-o" ) ) {

    /* make sure there are enough args */
    if( argc > nCurrentArg + 2 ) {
      
      /* read the overlay path and stem */
      strcpy( sOverlayPath, argv[nCurrentArg+1] );
      strcpy( sOverlayStem, argv[nCurrentArg+2] );
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
      strcpy( sTimeCoursePath, argv[nCurrentArg] );
      strcpy( sTimeCourseStem, argv[nCurrentArg+1] );
      bLoadingTimeCourse = TRUE;
      nCurrentArg += 2;
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
      if( '-' == argv[nCurrentArg+1][0] 
    || MATCH( argv[nCurrentArg+1], "local" )
    || MATCH( argv[nCurrentArg+1], "." ) ) {
        
        /* image dir is missing. */
        OutputPrint "Image directory is missing.\n" EndOutputPrint;
        bFatalError = TRUE;
        
      } else {
        
        /* read in subject and image name. */
        strcpy( sSubject, argv[nCurrentArg] );
        strcpy( sImageDir, argv[nCurrentArg+1] );
        bSubjectDeclared = TRUE;
        nCurrentArg += 2;
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


    /* parse cwd to set session root (guess path to bem) */
    /* %%%%Sean %%%% */
    /* getwd(cwd); */
#ifdef Linux
    cwd = getcwd(NULL,0);
#else
    cwd = getenv("PWD");
#endif

    /* if cwd is null, we couldn't get the current working directory, probably
       because of restrictions. we need it, tho, so ask the user to cd and
       bail. */
    if( NULL == cwd ) {
      OutputPrint "ERROR: Couldn't get current working directory information. Please change to a valid directory and try again.\n" EndOutputPrint;
      exit (1);
    }

    /* if using local directory, copy 'local' in subject name for
       historical reasons. */
    if( bLocalImageDir ) {
      strcpy( sSubject, "local" );
      strcpy( sImageDir, "" );
    }

    /* check for subjects dir env variable. if not defined and we're not
       using mriread or local, bail out. */
    lsubjectsdir = getenv("SUBJECTS_DIR");
    if ( lsubjectsdir == NULL 
   && !bUsingMRIRead
   && !bLocalImageDir ) {
      printf("medit: env var SUBJECTS_DIR undefined (use setenv)\n");
      exit(1);
    }

    /* make the path we should find data in */
    sprintf( fname, "%s/%s", lsubjectsdir, sSubject );

    /* if we're not using local dir and not using MRIRead... */
    if( !bUsingMRIRead && !bLocalImageDir ) {

      /* try to open the dir. */
      if ((fp=fopen (fname,"r"))==NULL) {
        printf("medit: ### can't find subject %s\n",fname); exit(1);}
      else fclose(fp);
    }

    if( bNoEdit ) {
      editflag = FALSE;
    }

    word = strtok(cwd,"/");
    strcpy(path[0],word);
    i = 1;
    while ((word = strtok(NULL,"/")) != NULL) {  /* save,count */
      strcpy(path[i],word);
      i++;
    }
    if (MATCH(path[i-1],"scripts") && 
        MATCH(path[i-2],sSubject)) {
      printf("medit: in subjects \"scripts\" dir\n");
      j = i-1;
    } else if (MATCH(path[i-1],"scripts") &&
               MATCH(path[i-2],"eegmeg") &&
               MATCH(path[i-3],sSubject)) {
      printf("medit: in subjects \"eegmeg/scripts\" dir\n");
      j = i-1;
    } else if (MATCH(path[i-1],"scripts")) {
      printf("medit: in \"scripts\" dir (not subjects,eegmeg)\n");
      j = i-1;
    } else if (MATCH(sSubject,"local")) {  /* local even if in mri */
      printf("medit: local data set => session root is cwd\n");
      j = i;
    } else if (MATCH(path[i-2],"mri")) {
      printf("medit: in subdir of \"mri\" dir\n");
      j = i-1;
    } else {
      printf(
  "medit: not in \"scripts\" dir or \"mri\" subdir => session root is cwd\n");
      j = i;
    }
    sprintf(lsrname,"/%s",path[0]);  /* reassemble absolute */
    for(i=1;i<j;i++)
      sprintf(lsrname,"%s/%s",lsrname,path[i]);
    printf("medit: session root data dir ($session) set to:\n");
    printf("medit:     %s\n",lsrname);

    /* if we loaded a surface */
    if( bSurfaceDeclared ) {
      
      /* set the stupid flag */
      surfflag = TRUE;

    } else {

      /* else make default surface name */
      strcpy( sSurface, "rh.orig" );
    }

    // make filenames.
    make_filenames(lsubjectsdir,lsrname,sSubject,sImageDir,sSurface);

    // if we want to load a surface, do so. however... if we are using MRIRead,
    // they specified a full surface name, so use that (lsurface). else
    // use the partial surface name concatenated with the SUBJECTS_DIR (sfname)
    if ( surfflag ) {

      if ( bUsingMRIRead ) 
  read_surface( sSurface );
      else
  read_surface( sfname );
    }

    // read in the volume using the selected method. sets many of our volume
    // related variables as well as allocates and fills out the volume array.
    if ( bUsingMRIRead ) {
      ReadVolumeWithMRIRead ( sSubject );
    } else {
      read_images ( mfname );
    }

    // now we can set up the bounds on our transforms.
    trans_SetBounds ( xx0, xx1, yy0, yy1, zz0, zz1 );
    trans_SetResolution ( ps, ps, st );

    /* load functional data */
    if( bLoadingOverlay ) {
      sprintf( sPathAndStem, "%s/%s", sOverlayPath, sOverlayStem );
      eFunctional = FunV_LoadOverlay( gFunctionalVolume, sPathAndStem );
      if( FunV_tErr_NoError != eFunctional ) {
  OutputPrint "ERROR: Couldn't load functional data.\n" EndOutputPrint;
      }
    }

    if( bLoadingTimeCourse ) {
      sprintf( sPathAndStem, "%s/%s", sTimeCoursePath, sTimeCourseStem );
      eFunctional = FunV_LoadTimeCourse( gFunctionalVolume, sPathAndStem );
      if( FunV_tErr_NoError != eFunctional ) {
  OutputPrint "ERROR: Couldn't load functional data.\n" EndOutputPrint;
      }
    }

    /* set functional color scale stuff */
    if( bThresh || bMid || bThresh ) {
      eFunctional = FunV_SetThreshold( gFunctionalVolume, min, mid, slope );
    }
    if( bTruncPhaseFlag ) {
      eFunctional = FunV_SetDisplayFlag( gFunctionalVolume, 
           FunV_tDisplayFlag_Ol_TruncateOverlay,
           (tBoolean) nTruncPhaseFlag );
    }
    if( bRevPhaseFlag ) {
      eFunctional = FunV_SetDisplayFlag( gFunctionalVolume, 
           FunV_tDisplayFlag_Ol_ReversePhase,
           (tBoolean) nRevPhaseFlag );
    }

    if (tclscriptflag) {
      /* called from tkmedit.c; do nothing (don't even open gl window) */
      /* wait for tcl interp to start; tkanalyse calls tcl script */
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

void
unmark_vertices () {

  int vno ;

  if (!mris) {
    fprintf(stderr, "no surface loaded.\n") ;
    return ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].marked = 0 ;

  redraw();
}

void
mark_file_vertices ( char *fname ) {

  FILE  *fp ;
  char  line[200], *cp ;
  int   vno, nvertices, nargs ;
  float area ;

  if ( !mris ) {

    fprintf(stderr, "no surface loaded.\n") ;
    return ;
  }

  fp = fopen(fname, "r") ;
  if ( !fp ) {
    fprintf(stderr, "could not open file %s.\n", fname) ;
    return ;
  }

  fgetl(line, 150, fp) ;
  nargs = sscanf(line, "%d %f", &nvertices, &area) ;
  if (nargs == 2)
    fprintf ( stderr, "marking %d vertices, %2.3f mm^2 surface area\n",
        nvertices, area ) ;
  else if (nargs == 1)
    fprintf(stderr, "marking %d vertices\n", nvertices) ;

  while  ((cp = fgetl(line, 150, fp)) != NULL) {

    sscanf(cp, "%d", &vno) ;
    if (vno >= 0 && vno < mris->nvertices) {
      mris->vertices[vno].marked = 1 ;
    }
  }

  fclose(fp) ;


}

void GotoSurfaceVertex ( tkm_tSurfaceType inSurface, int inVertex ) {

  VERTEX* theVertex;
  int theVoxX, theVoxY, theVoxZ;
  VoxelRef theVolumeVox;
  
  /* make sure surface is loaded. */
  if ( !mris ) {
    OutputPrint "Surface is not loaded.\n" EndOutputPrint;
    return;
  }

  /* make sure this value is in bounds */
  if ( inVertex < 0 
       || inVertex >= mris->nvertices ) {

    OutputPrint "%d is an invalid vertex number.\n", inVertex EndOutputPrint;
    return;
  }

  /* get the vertex. */
  theVertex = &(mris->vertices[inVertex]);
  if ( NULL == theVertex ) {
    DebugPrint "GotoSurfaceVertex( %d, %d ): Vertex was NULL.\n", 
      (int)inSurface, inVertex EndDebugPrint;
    OutputPrint "Error getting vertex from surface.\n" EndOutputPrint;
    return;
  }

  /* convert the coords to volume voxel. */
  switch ( inSurface ) {
  case tkm_tSurfaceType_Current:
    trans_RASToVoxelIndex ( theVertex->x, theVertex->y, theVertex->z,
          &theVoxX,     &theVoxY,     &theVoxZ ); 
    break;
  case tkm_tSurfaceType_Original:
    trans_RASToVoxelIndex ( theVertex->origx, theVertex->origy, theVertex->origz,
          &theVoxX,         &theVoxY,         &theVoxZ ); 
    break;
  case tkm_tSurfaceType_Canonical:
    trans_RASToVoxelIndex ( theVertex->cx, theVertex->cy, theVertex->cz,
          &theVoxX,      &theVoxY,      &theVoxZ ); 
    break;
  default:
    DebugPrint "GotoSurfaceVertex( %d, %d ): Invalid surface.\n",
      (int)inSurface, inVertex EndDebugPrint;
    OutputPrint "Error finding vertex.\n" EndOutputPrint;
    return;
  }

  /* build a voxel */
  Voxel_New ( &theVolumeVox );
  Voxel_Set ( theVolumeVox, theVoxX, theVoxY, theVoxZ );
         
  /* tell the window to go there. */
  MWin_SetCursor ( gMeditWindow, -1, theVolumeVox );
  
  /* tell the window to hilite the vertex. */
  MWin_HiliteSurfaceVertex ( gMeditWindow, -1, inSurface, inVertex );

  Voxel_Delete ( &theVolumeVox );
}

void FindNearestSurfaceVertex ( tkm_tSurfaceType inSurface ) {

  VERTEX* theVertex;
  VoxelRef theCursor = NULL;
  Real theRASX, theRASY, theRASZ;
  int theCurVertex, theClosestVertex;
  float dx, dy, dz, theDistance, theMinDistance;

  /* make sure we have a surface */
  if ( !mris ) {
    OutputPrint "Surface is not loaded.\n" EndOutputPrint;
    return;
  }
  /* get cursor and convert to RAS */
  Voxel_New ( &theCursor );
  MWin_GetCursor ( gMeditWindow, theCursor );
  trans_VoxelToRAS ( EXPAND_VOXEL_INT(theCursor), 
         &theRASX, &theRASY, &theRASZ );
  Voxel_Delete ( &theCursor );

  /* nothing yet */
  theMinDistance = 1e9;
  theClosestVertex = -1;

  /* for all vertices... */
  for ( theCurVertex = 0; theCurVertex < mris->nvertices; theCurVertex++ ) {

    /* the vertex */
    theVertex = &(mris->vertices[theCurVertex]);

    /* calc distance to cursor */
    switch ( inSurface ) {
    case tkm_tSurfaceType_Current:
      dx = theVertex->x - theRASX;
      dy = theVertex->y - theRASY;
      dz = theVertex->z - theRASZ;
      break;
    case tkm_tSurfaceType_Original:
      dx = theVertex->origx - theRASX;
      dy = theVertex->origy - theRASY;
      dz = theVertex->origz - theRASZ;
      break;
    case tkm_tSurfaceType_Canonical:
      dx = theVertex->cx - theRASX;
      dy = theVertex->cy - theRASY;
      dz = theVertex->cz - theRASZ;
      break;
    default:
      return;
    }

    theDistance = dx*dx + dy*dy + dz*dz;

    /* if less than min, save it. */
    if ( theDistance < theMinDistance ) {
      theClosestVertex = theCurVertex;
      theMinDistance = theDistance;
    }
  }

  if ( theCurVertex != -1 ) {

    theVertex = &(mris->vertices[theClosestVertex]);
    OutputPrint "Found vertex %d @ (%2.1f, %2.1f, %2.1f), %2.1f mm away\n",
      theClosestVertex, theVertex->x, theVertex->y, theVertex->z,
      sqrt(theMinDistance) EndOutputPrint;
    
    /* tell window to hilite that vertex */
    MWin_HiliteSurfaceVertex ( gMeditWindow, -1, inSurface, theClosestVertex );

  } else {
    
    OutputPrint "Nothing found!\n" EndOutputPrint;
  }
}

void WriteVoxelToControlFile ( char* inPathName, VoxelRef inVolumeVox ) {

  char fname[NAME_LENGTH];
  FILE *fp;
  Real theRASX, theRASY, theRASZ;
  VoxelRef theVoxel;

  // make a new voxel
  Voxel_New ( &theVoxel );
  
  sprintf(fname,"%s/control.dat",inPathName);
  fp=fopen(fname,"a+");
  
  if (fp==NULL) {
    printf("medit: ### can't create file %s\n",fname);
    return;
  }

  // convert voxel to ras
  VoxelToRAS ( EXPAND_VOXEL_INT(inVolumeVox), 
         &theRASX, &theRASY, &theRASZ );
  
  // write RAS space pt to file
  fprintf ( fp,"%f %f %f\n", theRASX, theRASY, theRASZ );
  DebugPrint "writing RAS point to %s...\n", fname EndDebugPrint;
  
  // close the file
  fclose(fp);

  // free the voxel
  Voxel_Delete ( &theVoxel );
}

void WriteVoxelToEditFile ( char* inPathName, VoxelRef inVolumeVox ) {

  char fname[NAME_LENGTH];
  FILE *fp;
  Real theTalX, theTalY, theTalZ;
  Real theRASX, theRASY, theRASZ;
  VoxelRef theVoxel;

  // make a new voxel
  Voxel_New ( &theVoxel );

  sprintf(fname,"%s/edit.dat",inPathName);
  fp=fopen(fname,"w");

  if (fp==NULL) {
    DebugPrint "Couldn't create file %s\n", fname EndDebugPrint;
    OutputPrint "Couldn't create edit file.\n" EndOutputPrint;
    return;
  }

  // convert to ras
  VoxelToRAS ( EXPAND_VOXEL_INT(inVolumeVox), 
         &theRASX, &theRASY, &theRASZ );
  
  // write RAS space pt to file
  fprintf ( fp,"%f %f %f\n", theRASX, theRASY, theRASZ );
  DebugPrint "writing RAS point to %s...\n", fname EndDebugPrint;
  
  // if we have a tal transform for this volume...
  if ( transform_loaded ) {
    
    // ras to tal
    transform_point ( linear_transform, theRASX, theRASY, theRASZ,
                      &theTalX, &theTalY, &theTalZ );
      
    // write tal space point to file
    fprintf(fp, "%f %f %f\n", theTalX, theTalY, theTalZ );
    DebugPrint "writing Tal point to %s...\n", fname EndDebugPrint;
  }
  
  // close the file
  fclose(fp);

  // free the voxel
  Voxel_Delete ( &theVoxel );
}

void ReadCursorFromEditFile ( char* inDir ) {

  char theFileName[NAME_LENGTH];
  FILE* theFile;
  float theRASX, theRASY, theRASZ;
  int theVoxX, theVoxY, theVoxZ;
  VoxelRef theVolumeVox;

  /* make the file name. */
  sprintf ( theFileName, "%s/edit.dat", inDir );

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
  trans_RASToVoxelIndex ( theRASX,  theRASY,  theRASZ,
        &theVoxX, &theVoxY, &theVoxZ ); 
 
  /* build and set cursor */
  Voxel_New ( &theVolumeVox );
  Voxel_Set ( theVolumeVox, theVoxX, theVoxY, theVoxZ );
  MWin_SetCursor ( gMeditWindow, -1, theVolumeVox );
  MWin_SetZoomCenterToCursor ( gMeditWindow, -1 );
  Voxel_Delete ( &theVolumeVox );
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

void
translate_brain(float a,char c)
{
  int i,j,k;
  float m1[4][4],m2[4][4];

  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    m1[i][j] = (i==j)?1.0:0.0;
  if (c=='y')
    m1[1][3] = a;
  else if (c=='x')
    m1[0][3] = a;
  else if (c=='z')
    m1[2][3] = a;
  else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
  {
    m2[i][j] = 0;
    for (k=0;k<4;k++)
      m2[i][j] += tm[i][k]*m1[k][j];
  }
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = m2[i][j];
}

void
rotate_brain(float a,char c)
{
  int i,j,k;
  float m1[4][4],m2[4][4];
  float sa,ca;

  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    m1[i][j] = (i==j)?1.0:0.0;
  a = a*M_PI/1800;
  sa = sin(a);
  ca = cos(a);
  if (c=='y')
  {
    m1[0][0] = m1[2][2] = ca;
    m1[2][0] = -(m1[0][2] = sa);
  } else
  if (c=='x')
  {
    m1[1][1] = m1[2][2] = ca;
    m1[1][2] = -(m1[2][1] = sa);
  } else
  if (c=='z')
  {
    m1[0][0] = m1[1][1] = ca;
    m1[0][1] = -(m1[1][0] = sa);
  } else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
  {
    m2[i][j] = 0;
    for (k=0;k<4;k++)
      m2[i][j] += tm[i][k]*m1[k][j];
  }
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = m2[i][j];
}

int
read_images(char *fpref)
{
  int i,j,k;
  FILE *fptr, *xfptr;
  char fname[STRLEN], char_buf[STRLEN];
  tVolumeValue theIntensity;
  tVolumeValue * theSlicePtr;
  tVolumeValue* buf;
  int bufsize;

  sprintf(fname,"%s.info",fpref);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {
    strcpy(fpref,"COR-");
    sprintf(fname,"%s.info",fpref);
    fptr = fopen(fname,"r");
    if (fptr==NULL) {
      printf("medit: ### File %s not found (tried local dir, too)\n",fname);
      exit(1);
    }
  } 

  fscanf(fptr,"%*s %d",&imnr0);
  fscanf(fptr,"%*s %d",&imnr1);
  fscanf(fptr,"%*s %d",&ptype);
  fscanf(fptr,"%*s %d",&xnum);
  fscanf(fptr,"%*s %d",&ynum);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %f",&ps);
  fscanf(fptr,"%*s %f",&st);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %f",&xx0); /* strtx */
  fscanf(fptr,"%*s %f",&xx1); /* endx */
  fscanf(fptr,"%*s %f",&yy0); /* strty */
  fscanf(fptr,"%*s %f",&yy1); /* endy */
  fscanf(fptr,"%*s %f",&zz0); /* strtz */
  fscanf(fptr,"%*s %f",&zz1); /* endz */
  fscanf(fptr, "%*s %*f") ;   /* tr */
  fscanf(fptr, "%*s %*f") ;   /* te */
  fscanf(fptr, "%*s %*f") ;   /* ti */
  fscanf(fptr, "%*s %s",char_buf);

  /* marty: ignore abs paths in COR-.info */
  xfptr = fopen(xffname,"r");
  if (xfptr==NULL)
    printf("medit: Talairach xform file not found (ignored)\n");
  else {
    fclose(xfptr);
    if (input_transform_file(xffname, &talairach_transform) != OK)
      printf("medit: ### File %s load failed\n",xffname);
    else {
      printf("medit: Talairach transform file loaded\n");
      printf("medit: %s\n",xffname);
      linear_transform = get_linear_transform_ptr(&talairach_transform) ;
      inverse_linear_transform =
                      get_inverse_linear_transform_ptr(&talairach_transform) ;
      transform_loaded = TRUE;
    }
  }

  ps *= 1000;
  st *= 1000;
  xx0 *= 1000;
  xx1 *= 1000;
  yy0 *= 1000;
  yy1 *= 1000;
  zz0 *= 1000;
  zz1 *= 1000;
  fclose(fptr);
  numimg = imnr1-imnr0+1;
  xdim=xnum*zf;
  ydim=ynum*zf;


  bufsize = ((unsigned long)xnum)*ynum;
  buf = (tVolumeValue*) lcalloc( bufsize, sizeof(tVolumeValue));

  /*
  for (k=0;k<numimg;k++) {
    im[k] = (tVolumeValue **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++) {
      im[k][i] = (tVolumeValue *)lcalloc(xnum,sizeof(char));
    }
  }
  */

  InitVolume ( &gAnatomicalVolume, xnum );
    
  for (k=0;k<6;k++)
  {
    sim[k] = (tVolumeValue **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++)
    {
      sim[k][i] = (tVolumeValue *)lcalloc(xnum,sizeof(char));
    }
  }

  for (k=0;k<numimg;k++)
  {
    changed[k] = FALSE;
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"rb");
    if(fptr==NULL){printf("medit: ### File %s not found.\n",fname);exit(1);}
    fread(buf,sizeof(char),bufsize,fptr);

    theSlicePtr = GetVolumeSlicePtr ( gAnatomicalVolume, k );
    memcpy ( theSlicePtr, buf, xnum*ynum ); 
   //    buffer_to_image ( buf, &theSlicePtr, xnum, ynum );

    fclose(fptr);
  }

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

  free (buf);

  return(0);
}


int
read_second_images(char *imdir2)
{
  int i,j,k,n;
  FILE *fptr;
  char fname[NAME_LENGTH], cmd[NAME_LENGTH], fpref[NAME_LENGTH];
  char fnamefirst[NAME_LENGTH], notilde[NAME_LENGTH];
  tVolumeValue theIntensity;
  tVolumeValue * theSlicePtr = NULL;
  tVolumeValue* buf;
  int bufsize;

  strcpy(imtype2, imdir2) ;

  bufsize = ((unsigned long)xnum)*ynum;
  buf = (tVolumeValue *)lcalloc(bufsize,sizeof(char));


  /* TODO: replace w/setfile (add mri/<subdir>); tk.c: omit arg like others */
  if (imdir2[0] == '/')
    sprintf(fpref,"%s/COR-",imdir2);
  else if (imdir2[0] == '~') {
    for(n=1;n<=strlen(imdir2);n++)  notilde[n-1] = imdir2[n];
    sprintf(fpref,"%s/%s/%s/COR-",subjectsdir,pname,notilde);
  }
  else if (MATCH(pname,"local"))
    sprintf(fpref,"%s/%s/COR-",srname,imdir2);
  else
    sprintf(fpref,"%s/%s/mri/%s/COR-",subjectsdir,pname,imdir2);

  sprintf(fnamefirst,"%s.info",mfname);
  sprintf(fname,"%s.info",fpref);
  sprintf(cmd,"diff %s %s",fnamefirst,fname);
  if (system(cmd)!=0)
    printf("medit: ### second COR-.info file doesn't match first--ignored\n");

  if (!second_im_allocated)
    alloc_second_im();
 
  for (k=0;k<numimg;k++) {
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"rb");
    if(fptr==NULL){printf("medit: ### File %s not found.\n",fname);return(0);}
    fread(buf,sizeof(char),bufsize,fptr);

    theSlicePtr = GetVolumeSlicePtr ( gAuxAnatomicalVolume, k );
    memcpy ( theSlicePtr, buf, xnum*ynum );
    //    buffer_to_image ( buf, &theSlicePtr, xnum, ynum );

    fclose(fptr);
  }

  for (k=0;k<numimg;k++)
    for (i=0;i<ynum;i++)
      for (j=0;j<xnum;j++) {

  theIntensity = GetVoxelValue ( gAuxAnatomicalVolume, j, i, k ) / 2;

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
    sim2[k][i][j] = sim2[k+3][i][j];

  /* set the aux volume in the medit window */
  MWin_SetAuxVolume( gMeditWindow, -1, 
         gAuxAnatomicalVolume, gVolumeDimension );

  /* show the aux value in the tk window */
  tkm_SendTclCommand ( tkm_tTclCommand_ShowAuxValue, "1" );

  redraw();

  free (buf);

  return(0);
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
  if( !MATCH( fpref, mfname )) {
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
      file_name(fpref,fname,k+imnr0,"%03d");
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
int
dump_vertex(int vno)
{
  VERTEX *v, *vn ;
  int    n ;
  float  dx, dy, dz, dist ;

  if (!mris)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "%s: surface must be loaded\n", Progname)) ;

  v = &mris->vertices[vno] ;
  fprintf(stderr, 
          "v %d: x = (%2.2f, %2.2f, %2.2f), n = (%2.2f, %2.2f, %2.2f)\n",
          vno, v->x, v->y, v->z, v->nx, v->ny, v->nz) ;
  if (curvloaded)
    fprintf(stderr, "curv=%2.2f\n", v->curv) ;
  fprintf(stderr, "nbrs:\n") ;
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    dx = vn->x - v->x ; dy = vn->y - v->y ; dz = vn->z - v->z ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    fprintf(stderr, "\tvn %d, dist = %2.2f\n", v->v[n], dist) ;
  }
  return(NO_ERROR) ;
}

void tkm_HandleIdle () {

  // just call the tk event handling function
  Tk_DoOneEvent ( TK_ALL_EVENTS | TK_DONT_WAIT );
}


void
smooth_surface(int niter)
{
  if (mris)
    MRISaverageVertexPositions(mris, niter) ;
  redraw() ;
}
int 
read_canonical_vertex_positions(char *fname)
{
  if (!mris)
    return(NO_ERROR) ;

  // kt - save current vertices to temp
  MRISsaveVertexPositions ( mris, TMP_VERTICES );

  // read vertices into current
  if (MRISreadVertexPositions(mris, fname) != NO_ERROR) {

    DebugPrint "read_canonical_vertex_positions ( %s ):\n\tcould not read canonical vertex positions\n", fname EndDebugPrint;
    OutputPrint "Couldn't read canonical vertices from %s.\n", 
      fname EndOutputPrint;

  } else {

    // save current to canonical
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

    // surface is loaded.
    gIsCanonicalSurfaceLoaded = TRUE;
    OutputPrint "Canonical surface loaded.\n" EndOutputPrint;
  }

  // restore current from temp
  MRISrestoreVertexPositions ( mris, TMP_VERTICES );
  
  return(NO_ERROR) ;
}
int
read_orig_vertex_positions(char *name)
{
  if (!mris)
    return(NO_ERROR) ;

  // kt - save current vertices to temp
  MRISsaveVertexPositions ( mris, TMP_VERTICES );

  // read vertices into current
  if (MRISreadVertexPositions(mris, name) != NO_ERROR) {

    DebugPrint "read_orig_vertex_positions ( %s ):\n\tcould not read original vertex positions\n", name EndDebugPrint;
    OutputPrint "Couldn't read original vertices from %s.\n", 
      name EndOutputPrint;

  } else {

    // save current to canonical
    MRISsaveVertexPositions(mris, ORIG_VERTICES) ;
    
    // surface is loaded.
    gIsOriginalSurfaceLoaded = TRUE;
    OutputPrint "Original surface loaded.\n" EndOutputPrint;
  }

  // kt - restore current from temp
  MRISrestoreVertexPositions ( mris, TMP_VERTICES );
  
  return(NO_ERROR) ;
}

int
read_surface(char *name)
{
  char fname[STRLEN], *cp ;
  int theStatus;

  cp = strchr(name, '/') ;
  if (!cp)  /* no path specified - put the path into it */
  {
    strcpy(surface, name) ;
    sprintf(sfname,"%s/%s/surf/%s",subjectsdir,pname,name);
    strcpy(fname, sfname) ;
  }
  else
    strcpy(fname, name) ;

  OutputPrint "Reading surface file %s\n", fname EndOutputPrint;

  theStatus = read_binary_surface(fname);

  // kt - if the return is 0, it seems the call was successful.
  if ( theStatus == 0 ) {

    // mark the surface as loaded.
    gIsCurrentSurfaceLoaded = TRUE;
    
    // print a status msg.
    OutputPrint "Surface loaded.\n" EndOutputPrint;

    // now load the other surfaces.
    read_orig_vertex_positions ( "orig" );
    read_canonical_vertex_positions ( "pial" );

  } else {

    DebugPrint "read_surface( %s ):\n\tread_binary_surface( %s ) failed, returned %d\n", name, fname, theStatus EndDebugPrint;
    OutputPrint "Surface failed to load.\n" EndOutputPrint;
  }

  /* set the medit window surface. */
  //MWin_SetSurface( gMeditWindow, -1, mris );
  
  return theStatus;
}

int
read_binary_surface(char *fname) {
  
  if (mris)
    MRISfree(&mris) ;
  mris = MRISread(fname) ;

  if (!mris) {
    surfflag = FALSE;
    surfloaded = FALSE;
    curvflag = FALSE;
    curvloaded = FALSE;
    fieldsignflag = FALSE;
    gIsCurrentSurfaceLoaded = FALSE;
    return(ERROR_NOFILE) ;
  }

  surfflag = TRUE;
  surfloaded = TRUE;
  curvflag = FALSE;
  curvloaded = FALSE;
  fieldsignflag = FALSE;
  gIsCurrentSurfaceLoaded = TRUE;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  return(0);
}

void UnloadSurface () {

  if ( !mris ) {
    OutputPrint "Surface not loaded.\n" EndOutputPrint;
    return;
  }

  /* free the surface. */
  MRISfree ( &mris );
  mris = NULL;

  /* update the medit window. */
  MWin_SetSurface( gMeditWindow, -1, mris );
}


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
read_hpts(char *fname)      /* sprintf(fname,"%s","../bem/head2mri.hpts"); */
{
  FILE *fp;
  char line[NAME_LENGTH];
  int i;


  fp = fopen(fname,"r");
  if (fp==NULL) {printf("medit: ### File %s not found\n",fname); return;}

  ptsflag = TRUE;
  fgets(line,NAME_LENGTH,fp);
  i = 0;
  while (!feof(fp))
  {
    jold[i] = iold[i] = imold[i] = 0;
    if (line[0]=='2') /* .show format? */
      sscanf(line,"%*s%*s%*s%*s%*s%*s%*s%f%*s%f%*s%f",
             &ptx[i],&pty[i],&ptz[i]);
    else /* .apos format */
    {
      sscanf(line,"%*s%*s%f%f%f",&ptx[i],&pty[i],&ptz[i]);
      ptx[i] *= 1000; /* convert from m to mm */
      pty[i] *= 1000;
      ptz[i] *= 1000;
    }
    fgets(line,NAME_LENGTH,fp);
    i++;
  }
  printf("head points file %s read\n",fname);
  npts = i;
  fclose(fp);
}

void
read_htrans(char *fname)   /* sprintf(fname,"%s","../bem/head2mri.trans"); */
{
  int i,j,k;
  FILE *fp;
  fp = fopen(fname,"r");
  if (fp==NULL) {
    printf("medit: no existing htrans file: making default\n");
    for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      tm[i][j] = (i==j);
  }
  else {
    for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      fscanf(fp,"%f",&tm[i][j]);
    printf("transformation file %s read\n",fname);
  }
  for (k=0;k<MAXPARS;k++)
    par[k] = dpar[k] = 0;
}

void
write_htrans(char *fname)   /* sprintf(fname,"%s","../bem/head2mri.trans"); */
{
  int i,j;
  FILE *fptr;

  fptr = fopen(fname,"w");
  if (fptr==NULL) {printf("medit: can't create file %s\n",fname);return;}
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++)
      fprintf(fptr,"%13.6e ",tm[i][j]);
    fprintf(fptr,"\n");
  }
  fclose(fptr);
  printf("transformation file %s written\n",fname);
}

void
make_filenames(char *lsubjectsdir,char *lsrname, char *lpname, char *limdir,
               char *lsurface)
{
  char tmpsurface[NAME_LENGTH], surfhemi[NAME_LENGTH];
    
  subjectsdir = (char *)malloc(NAME_LENGTH*sizeof(char)); /* malloc for tcl */
  srname = (char *)malloc(NAME_LENGTH*sizeof(char));
  pname = (char *)malloc(NAME_LENGTH*sizeof(char));
  title_str = (char *)malloc(NAME_LENGTH*sizeof(char));
  imtype = (char *)malloc(NAME_LENGTH*sizeof(char));
  imtype2 = (char *)malloc(NAME_LENGTH*sizeof(char));
  surface = (char *)malloc(NAME_LENGTH*sizeof(char));
  mfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  sfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  tfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  dipfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  decfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  hpfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  htfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  sgfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  fsfname = (char *)malloc(NAME_LENGTH*sizeof(char)); /* fscontour */
  fmfname = (char *)malloc(NAME_LENGTH*sizeof(char)); /* fscontour */
  cfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  xffname = (char *)malloc(NAME_LENGTH*sizeof(char));
  rfname = (char *)malloc(NAME_LENGTH*sizeof(char));

  strcpy(subjectsdir,lsubjectsdir);
  strcpy(srname,lsrname);
  strcpy(pname,lpname);
  strcpy(imtype,limdir);
  strcpy(surface,lsurface);
  sprintf(mfname,"%s/%s/mri/%s/COR-",subjectsdir,pname,imtype);
  sprintf(sfname,"%s/%s/surf/%s",subjectsdir,pname,surface);
  sprintf(tfname,"%s/%s/%s",subjectsdir,pname,TMP_DIR);
  /* next relative to session */
  sprintf(dipfname,"%s/bem/brain3d%d.dip",srname,dip_spacing);
  sprintf(decfname,"%s/bem/brain3d%d.dec",srname,dip_spacing);
  sprintf(hpfname,"%s/bem/head2mri.hpts",srname);
  sprintf(htfname,"%s/bem/head2mri.trans",srname);
  sprintf(sgfname,"/tmp/tkmedit.rgb");
  strcpy(tmpsurface,surface);  /* strtok writes in arg it parses! */
  strcpy(surfhemi,strtok(tmpsurface,"."));
  sprintf(fsfname,"%s/fs/%s.fs",srname,surfhemi);
  sprintf(fmfname,"%s/fs/%s.fm",srname,surfhemi);
  sprintf(cfname,"%s/%s/surf/%s.curv",subjectsdir,pname,surfhemi);
  sprintf(xffname,"%s/%s/%s/%s",
                     subjectsdir,pname,TRANSFORM_DIR,TALAIRACH_FNAME);
}

void
mirror(void)
{
  flip_corview_xyz("-x","y","z");
}

void
read_fieldsign(char *fname)  /* fscontour */
{
  int k,vnum;
  float f;
  FILE *fp;

  if (!surfloaded) {printf("medit: ### surface %s not loaded\n",fname); return;}

  fp = fopen(fname,"r");
  if (fp==NULL) {printf("medit: ### File %s not found\n",fname); return;}
  fread(&vnum,1,sizeof(int),fp);
  printf("medit: mris->nvertices = %d, vnum = %d\n",mris->nvertices,vnum);
  if (vnum!=mris->nvertices) {
    printf("medit: ### incompatible vertex number in file %s\n",fname);
    return; }
  for (k=0;k<vnum;k++)
  {
    fread(&f,1,sizeof(float),fp);
    mris->vertices[k].fieldsign = f;
  }
  fclose(fp);
  fieldsignflag = TRUE;
  surflinewidth = 3;
  
}

void
read_fsmask(char *fname)  /* fscontour */
{
  int k,vnum;
  float f;
  FILE *fp;

  if (!surfloaded) {printf("medit: ### surface %s not loaded\n",fname); return;}

  printf("medit: read_fsmask(%s)\n",fname);
  fp = fopen(fname,"r");
  if (fp==NULL) {printf("medit: ### File %s not found\n",fname); return;}
  fread(&vnum,1,sizeof(int),fp);
  if (vnum!=mris->nvertices) {
    printf("medit: ### incompatible vertex number in file %s\n",fname);
    return; }
  for (k=0;k<vnum;k++)
  {
    fread(&f,1,sizeof(float),fp);
    mris->vertices[k].fsmask = f;
  }
  fclose(fp);
  
}


void
read_binary_curvature(char *fname)
{
  float curvmin, curvmax, curv;
  int   k;

  MRISreadCurvatureFile(mris, fname) ;

  curvmin= 1000000.0f ; curvmax = -curvmin;
  for (k=0;k<mris->nvertices;k++)
  {
    curv = mris->vertices[k].curv;
    if (curv>curvmax) curvmax=curv;
    if (curv<curvmin) curvmin=curv;
  }
  printf("medit: curvature read: min=%f max=%f\n",curvmin,curvmax);
  curvloaded = TRUE;
  curvflag = TRUE;
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

  /* args: [-]{x,y,z} [-]{x,y,z} [-]{x,y,z} */
void
flip_corview_xyz(char *newx, char *newy, char *newz)
{
  int fx,fy,fz;
  int lx,ly,lz;
  int xx,xy,xz,yx,yy,yz,zx,zy,zz;
  int i,j,k;
  int ni,nj,nk;
  tVolumeValue theIntensity;

  ni = nj = nk = 0;
  fx=(newx[0]=='-')?1:0;
  fy=(newy[0]=='-')?1:0;
  fz=(newz[0]=='-')?1:0;
  lx = strlen(newx)-1;
  ly = strlen(newy)-1;
  lz = strlen(newz)-1;
  if(newx[lx]==newy[ly] || newx[lx]==newz[lz] || newy[ly]==newz[lz]) {
    printf("medit: ### degenerate flip\n"); return; }
  xx=(newx[lx]=='x')?1:0; xy=(newx[lx]=='y')?1:0; xz=(newx[lx]=='z')?1:0;
  yx=(newy[ly]=='x')?1:0; yy=(newy[ly]=='y')?1:0; yz=(newy[ly]=='z')?1:0;
  zx=(newz[lz]=='x')?1:0; zy=(newz[lz]=='y')?1:0; zz=(newz[lz]=='z')?1:0;

  if (!dummy_im_allocated) {
    for (k=0;k<numimg;k++) {
      dummy_im[k] = (tVolumeValue **)lcalloc(ynum,sizeof(char *));
      for (i=0;i<ynum;i++) {
        dummy_im[k][i] = (tVolumeValue *)lcalloc(xnum,sizeof(char));
      }
    }
    dummy_im_allocated = TRUE;
  }

  for (k=0;k<numimg;k++) {
    printf(".");
    for (i=0;i<ynum;i++)
    for (j=0;j<xnum;j++) {
      if(xx) nj=j; if(xy) nj=i; if(xz) nj=k;
      if(yx) ni=j; if(yy) ni=i; if(yz) ni=k;
      if(zx) nk=j; if(zy) nk=i; if(zz) nk=k;
      if(fx) nj=255-nj;
      if(fy) ni=255-ni;
      if(fz) nk=255-nk;
      dummy_im[k][i][j] = GetVoxelValue ( gAnatomicalVolume, nj, ni, nk );
      // dummy_im[k][i][j] = im[nk][ni][nj];
    }
    fflush(stdout);
  }
  printf("\n");

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

  printf("medit: imageflip done--to overwrite current COR: write_images\n");
  redraw();
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
  tkm_tOrientation theOrientation;
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
  tkm_tOrientation theOrientation;
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
  tkm_tOrientation theOrientation;
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




/*----------------end medit.c -------------------------------*/

/* %%%%Sean %%% */


/* boilerplate wrap function defines for easier viewing */
#define TCL_ARGS ( ClientData clientData, \
                   Tcl_Interp *interp,    \
                   int argc,char *argv[] )

/*=======================================================================*/
/* function wrappers and errors */

int TclRotateBrainX ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {
  
  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: RotateBrainX degrees",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    rotate_brain ( atof(argv[1]), 'x' ); 
  }

  return TCL_OK;
}

int TclRotateBrainY ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: RotateBrainY degrees",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    rotate_brain ( atof(argv[1]), 'y' ); 
  }

  return TCL_OK;
}

int TclRotateBrainZ ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: RotateBrainZ degrees",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    rotate_brain ( atof(argv[1]), 'z' ); 
  }

  return TCL_OK;
}

int TclTranslateBrainX ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: TranslateBrainX distance",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    translate_brain ( atof(argv[1]), 'x' ); 
  }

  return TCL_OK;
}

int TclTranslateBrainY ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: TranslateBrainY distance",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    translate_brain ( atof(argv[1]), 'y' ); 
  }

  return TCL_OK;
}

int TclTranslateBrainZ ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: TranslateBrainZ distance",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    translate_brain ( atof(argv[1]), 'z' ); 
  }

  return TCL_OK;
}

int TclWriteDipoles ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WriteDipoles",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    write_dipoles ( dipfname );
  }

  return TCL_OK;
}

int TclWriteDecimation ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WriteDecimation",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    write_decimation ( decfname );
  }

  return TCL_OK;
}

int TclMarkFileVertices ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: MarkFileVertices filename:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    mark_file_vertices ( argv[1] );
  }

  return TCL_OK;
}

int TclUnmarkVertices ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: UnmarkVertices",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    unmark_vertices ();
  }

  return TCL_OK;
}

int TclReadHPts ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadHPts",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_hpts ( hpfname );
  }

  return TCL_OK;
}

int TclReadHTrans ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadHTrans",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_htrans ( htfname );
  }

  return TCL_OK;
}

int TclWriteHTrans ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WriteHTrans",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    write_htrans ( htfname );
  }

  return TCL_OK;
}

int TclSaveRGB ( ClientData inClientData, Tcl_Interp* inInterp,
     int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveRGB filename:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    save_rgb ( argv[1] );
  }

  return TCL_OK;
}

int TclMirror ( ClientData inClientData, Tcl_Interp* inInterp,
    int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: Mirror",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    mirror ();
  }

  return TCL_OK;
}

int TclReadFieldSign ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadFieldSign",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_fieldsign ( fsfname );
  }

  return TCL_OK;
}

int TclReadFSMask ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadFSMask",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_fsmask ( fmfname );
  }

  return TCL_OK;
}

int TclReadDefaultBinaryCurvature ( ClientData inClientData, 
            Tcl_Interp* inInterp,
            int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadDefaultBinaryCurvature",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_binary_curvature ( cfname );
  }

  return TCL_OK;
}

int TclReadBinaryCurvature ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: ReadBinaryCurvature filename:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_binary_curvature ( argv[1] );
  }

  return TCL_OK;
}

int TclSmooth3D ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: Smooth3D steps:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    smooth_3d ( atoi(argv[1]) );
  }

  return TCL_OK;
}

int TclFlipCorView ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc < 4 ) {
    Tcl_SetResult ( inInterp, 
           "wrong # args: FlipCorView [-]{x,y,z}> <[-]{x,y,z}> <[-]{x,y,z}",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    flip_corview_xyz ( argv[1], argv[2], argv[3] );
  }

  return TCL_OK;
}

int TclReadBinarySurface ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadBinarySurface",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_binary_surface ( sfname );
  }

  return TCL_OK;
}

int TclDumpVertex ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: DumpVertex index:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    dump_vertex ( atoi(argv[1]) );
  }

  return TCL_OK;
}

int TclSmoothSurface ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SmoothSurface steps:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    smooth_surface ( atoi(argv[1]) );
  }

  return TCL_OK;
}

int TclWMFilterCorSlice ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  VoxelRef theCursor = NULL;

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WMFilterCorSlice",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    Voxel_New ( &theCursor );
    MWin_GetCursor ( gMeditWindow, theCursor );
    wmfilter_corslice ( Voxel_GetZ(theCursor) );
    Voxel_Delete ( &theCursor );
  }

  return TCL_OK;
}

int TclNormSlice ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  VoxelRef theCursor = NULL;

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: NormSlice {0=PostAnt,1=InfSup,2=LeftRight}",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    Voxel_New ( &theCursor );
    MWin_GetCursor ( gMeditWindow, theCursor );
    norm_slice ( Voxel_GetZ(theCursor), 
     Voxel_GetY(theCursor),
     Voxel_GetX(theCursor),
     atoi(argv[1]) ); 
    Voxel_Delete ( &theCursor );
  }

  return TCL_OK;
}

int TclNormAllSlices ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc < 2 ) {
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

  if ( argc < 4 ) {
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

int TclLoadAuxVolume ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadAuxVolume image_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_second_images ( argv[1] );
  }

  return TCL_OK;
}

int TclSaveVolume ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveVolume",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    write_images ( mfname );
  }

  return TCL_OK;
}

int TclSaveVolumeAs ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc < 2 ) {
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

  if ( argc < 1 ) {
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

  if ( argc < 1 ) {
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

  if ( argc < 1 ) {
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

  if ( argc < 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetVolumeColorScale threshold:float squash:float",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    SetVolumeColorScale ( atof( argv[1] ), atof( argv[2] ) );
  }

  return TCL_OK;
}

int TclSaveLabel ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc < 2 ) {
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

  if ( argc < 2 ) {
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

  if ( argc < 1 ) {
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
  
  if ( argc < 1 ) {
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

  VoxelRef theCursor = NULL;

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SendCursor",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) { 
    Voxel_New ( &theCursor );
    MWin_GetCursor ( gMeditWindow, theCursor );
    WriteVoxelToEditFile ( tfname, theCursor );
    Voxel_Delete ( &theCursor );
  }    

  return TCL_OK;
}

int TclReadCursor ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadCursor",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    ReadCursorFromEditFile ( tfname );
  }

  return TCL_OK;
}


int TclUndoLastEdit ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {

  if ( argc < 1 ) {
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

  VoxelRef theCursor = NULL;

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: NewControlPoint",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    Voxel_New ( &theCursor );
    MWin_GetCursor ( gMeditWindow, theCursor );
    NewCtrlPt ( theCursor );
    Voxel_Delete ( &theCursor );
  }

  return TCL_OK;
}

int TclDeselectAllControlPoints ( ClientData inClientData, 
          Tcl_Interp* inInterp,
          int argc, char* argv[] ) {

  if ( argc < 1 ) {
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

  if ( argc < 1 ) {
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

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: WriteControlPointFile",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    WriteCtrlPtFile ( tfname );
  }

  return TCL_OK;
}

int TclLoadMainSurface ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadMainSurface surface_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_surface ( argv[1] );
  }

  return TCL_OK;
}

int TclLoadCanonicalSurface ( ClientData inClientData, Tcl_Interp* inInterp,
            int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadCanonicalSurface surface_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_canonical_vertex_positions ( argv[1] );
  }

  return TCL_OK;
}

int TclLoadOriginalSurface ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadOriginalSurface surface_name:string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    read_orig_vertex_positions ( argv[1] );
  }  

  return TCL_OK;
}

int TclUnloadAllSurfaces ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {

  if ( argc < 1 ) {
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

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: GotoMainVertex vertex_num:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    GotoSurfaceVertex ( tkm_tSurfaceType_Current, atoi( argv[1] ) );
  }  

  return TCL_OK;
}

int TclGotoCanonicalVertex ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: GotoCanonicalVertex vertex_num:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    GotoSurfaceVertex ( tkm_tSurfaceType_Canonical, atoi( argv[1] ) );
  }  

  return TCL_OK;
}

int TclGotoOriginalVertex ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {

  if ( argc < 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: GotoOriginalVertex vertex_num:int",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    GotoSurfaceVertex ( tkm_tSurfaceType_Original, atoi( argv[1] ) );
  }  

  return TCL_OK;
}

int TclShowNearestMainVertex ( ClientData inClientData, Tcl_Interp* inInterp,
             int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestMainVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    FindNearestSurfaceVertex ( tkm_tSurfaceType_Current );
  }  

  return TCL_OK;
}

int TclShowNearestOriginalVertex ( ClientData inClientData, 
           Tcl_Interp* inInterp,
           int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestOriginalVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    FindNearestSurfaceVertex ( tkm_tSurfaceType_Original );
  }  

  return TCL_OK;
}

int TclShowNearestCanonicalVertex ( ClientData inClientData, 
            Tcl_Interp* inInterp,
            int argc, char* argv[] ) {

  if ( argc < 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestCanonicalVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if( gbAcceptingTclCommands ) {
    FindNearestSurfaceVertex ( tkm_tSurfaceType_Canonical );
  }  

  return TCL_OK;
}


int TclLoadParcellationVolume ( ClientData inClientData, 
        Tcl_Interp* inInterp,
        int argc, char* argv[] ) {

  if ( argc < 3 ) {
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

  if ( argc < 1 ) {
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
  char theErr;
  FunV_tErr eFunctional = FunV_tErr_NoError;

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
  InitDebugging;
  EnableDebuggingOutput;

  DebugPrint "Debugging output is on.\n" EndDebugPrint;

  // init the selection list 
  theErr = VList_New ( &gSelectionList );
  if ( theErr ) {
    DebugPrint "Error in VList_Init: %s\n", 
      VList_GetErrorString(theErr) EndDebugPrint;
    exit (1);
  }

  // init our control pt list
  theErr = VSpace_New ( &gCtrlPtList );
  if ( theErr != kVSpaceErr_NoErr ) {
    DebugPrint "Error in VSpace_Init: %s\n",
      VSpace_GetErrorString ( theErr ) EndDebugPrint;
    exit (1);
  }

  // init the undo list.
  InitUndoList ();

  // and the selection module.
  InitSelectionModule ();

  // no surfaces are loaded.
  gIsCurrentSurfaceLoaded = FALSE;
  gIsOriginalSurfaceLoaded = FALSE;
  gIsCanonicalSurfaceLoaded = FALSE;

  /* create functional volume */
  eFunctional = FunV_New( &gFunctionalVolume,
        UpdateAndRedraw, tkm_SendTclCommand, 
        SendTCLCommand );
  if( FunV_tErr_NoError != eFunctional ) {
    DebugPrint "Error creating functional volume.\n" EndDebugPrint;
    OutputPrint "Couldn't initialize functional module.\n" EndOutputPrint;
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

  // unles flag is set, try local script first.
  if ( !getenv ( "DONT_USE_LOCAL_TKMEDIT_TCL" ) ) {
    sprintf(tkmedit_tcl,"%s","tkmedit.tcl"); 
    fp = fopen ( tkmedit_tcl,"r" );
  }

  // if file is not open, try the normal place.
  if ( NULL == fp ) { 
    sprintf(tkmedit_tcl,"%s/lib/tcl/%s",envptr,"tkmedit.tcl"); 
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
  ProcessCtrlPtFile( tfname );

  /* init medit window */
  glutInit            ( &argc, argv );
  glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
  MWin_New ( &gMeditWindow, pname, 512, 512 );

  /* set window's data */
  MWin_SetVolume ( gMeditWindow, -1, gAnatomicalVolume, 256 );
  MWin_SetControlPointsSpace ( gMeditWindow, -1, gCtrlPtList );
  MWin_SetControlPointsSelectionList ( gMeditWindow, -1, gSelectionList );
  MWin_SetSelectionSpace ( gMeditWindow, -1, gSelectedVoxels );
  if ( gIsCurrentSurfaceLoaded )
    MWin_SetSurface ( gMeditWindow, 0, mris );

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
  //if (Tix_Init(interp) == TCL_ERROR ) {
  //  fprintf(stderr, "Tix_Init failed: %s\n", interp->result); }

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

  Tcl_CreateCommand ( interp, "RotateBrainX",
          TclRotateBrainX,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RotateBrainY",
          TclRotateBrainY,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RotateBrainZ",
          TclRotateBrainZ,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "TranslateBrainX",
          TclTranslateBrainX,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "TranslateBrainY",
          TclTranslateBrainY,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "TranslateBrainZ",
          TclTranslateBrainZ,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "WriteDipoles",
          TclWriteDipoles,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "WriteDecimation",
          TclWriteDecimation,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "MarkFileVertices",
          TclMarkFileVertices,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "UnmarkVertices",
          TclUnmarkVertices,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadHPts",
          TclReadHPts,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadHTrans",
          TclReadHTrans,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "WriteHTrans",
          TclWriteHTrans,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveRGB",
          TclSaveRGB,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "Mirror",
          TclMirror,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadFieldSign",
          TclReadFieldSign,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadFSMask",
          TclReadFSMask,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadDefaultBinaryCurvature",
          TclReadDefaultBinaryCurvature,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadBinaryCurvature",
          TclReadBinaryCurvature,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "Smooth3D",
          TclSmooth3D,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "FlipCorView",
          TclFlipCorView,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadBinarySurface",
          TclReadBinarySurface,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "DumpVertex",
          TclDumpVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SmoothSurface",
          TclSmoothSurface,
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

  strcpy(rfname,script_tcl);  /* save in global (malloc'ed in Program) */

  /* run tcl/tk startup script to set vars, make interface; no display yet */
  code = Tcl_EvalFile(interp,tkmedit_tcl);
  if (*interp->result != 0)  printf(interp->result);

  /* if command line script exists, now run as batch job (possibly exiting) */
  if (scriptok) {    /* script may or may not open gl window */
    printf("tkmedit: run tcl script: %s\n",script_tcl);
    code = Tcl_EvalFile(interp,script_tcl);
    if (*interp->result != 0)  printf(interp->result);
  } else {
    ; /* tkmedit has already opened gl window if no script */
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

  /* set data in window */
  MWin_SetOverlayVolume( gMeditWindow, -1, gFunctionalVolume );

  /* if we have a tal transform, hide the ras coords. otherwise hide
     the tal coords. */
  if ( transform_loaded ) {
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
  SetVolumeColorScale ( kfDefaultVolumeThreshold, kfDefaultVolumeSquash );

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
  
  char theErr;

  /* delete window */
  MWin_Delete( &gMeditWindow );

  // shut down tcl stuff
  SendTCLCommand ( "exit" );

  // delete everything we allocated before
  FunV_Delete( &gFunctionalVolume );

  DeleteSelectionModule ();

  theErr = VList_Delete ( &gSelectionList );
  if ( theErr )
    DebugPrint "Error in VList_Delete: %s\n",  
      VList_GetErrorString(theErr) EndDebugPrint;

  theErr = VSpace_Delete ( &gCtrlPtList );
  if ( theErr != kVSpaceErr_NoErr )
    DebugPrint "Error in VSpace_Delete: %s\n",
      VSpace_GetErrorString(theErr) EndDebugPrint;

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

void ReadVolumeWithMRIRead ( char * inFileOrPath ) {

  MRI *theVolume;
  int theSlice, theRow, theCol;
  tVolumeValue theIntensity;
  tVolumeValue * theSlicePtr;

  // pass the path to MRIRead
  theVolume = MRIread ( inFileOrPath );

  // make sure the result is good.
  if ( NULL == theVolume ) {
    OutputPrint "Couldn't read volume data at %s\n", 
      inFileOrPath EndOutputPrint;
    exit ( 1 );
  }

  // conform it.
  theVolume = MRIconform ( theVolume );
  
  // grab all the data we need.
  imnr0 = theVolume->imnr0;
  imnr1 = theVolume->imnr1;
  ptype = theVolume->ptype;
  xnum = theVolume->width;
  ynum = theVolume->height;
  ps = theVolume->ps;
  st = theVolume->thick;
  xx0 = theVolume->xstart;
  xx1 = theVolume->xend;
  yy0 = theVolume->ystart;
  yy1 = theVolume->yend;
  zz0 = theVolume->zstart;
  zz1 = theVolume->zend;

  // grab the tal transforms.
  if ( NULL != theVolume->linear_transform ) {
    copy_general_transform ( &theVolume->transform, &talairach_transform );
    linear_transform = get_linear_transform_ptr ( &talairach_transform );
    inverse_linear_transform = 
      get_inverse_linear_transform_ptr ( &talairach_transform );
  }

  // if we got them, make note of it.
  if ( NULL != linear_transform  && NULL != inverse_linear_transform )
    transform_loaded = TRUE;
  else
    transform_loaded = FALSE;

  numimg = imnr1-imnr0+1; // really the number of slices

  // calc window dimensions according to already defined scale factor
  xdim= xnum * zf;
  ydim= ynum * zf;

  /*
  DebugPrint "\timnr0 = %d\n", imnr0 EndDebugPrint;
  DebugPrint "\timnr1 = %d\n", imnr1 EndDebugPrint;
  DebugPrint "\tnumimg = %d\n", numimg EndDebugPrint;
  DebugPrint "\tptype = %d\n", ptype EndDebugPrint;
  DebugPrint "\txnum = %d\n", xnum EndDebugPrint;
  DebugPrint "\tynum = %d\n", ynum EndDebugPrint;
  DebugPrint "\tps = %2.5f\n", ps EndDebugPrint;
  DebugPrint "\tst = %2.5f\n", st EndDebugPrint;
  DebugPrint "\txx0 = %2.5f\n", xx0 EndDebugPrint;
  DebugPrint "\txx1 = %2.5f\n", xx1 EndDebugPrint;
  DebugPrint "\tyy0 = %2.5f\n", yy0 EndDebugPrint;
  DebugPrint "\tyy1 = %2.5f\n", yy1 EndDebugPrint;
  DebugPrint "\tzz0 = %2.5f\n", zz0 EndDebugPrint;
  DebugPrint "\tzz1 = %2.5f\n", zz1 EndDebugPrint;
  DebugPrint "\ttransform_loaded = %d\n", transform_loaded EndDebugPrint;
  */

  InitVolume ( &gAnatomicalVolume, xnum );

  // sim[][][] is a small volume buffer.
  for ( theSlice = 0 ; theSlice < 6; theSlice++ ) {
    sim[theSlice] = (BUFTYPE**)lcalloc ( ynum, sizeof(BUFTYPE*) );
    for ( theRow = 0; theRow < ynum; theRow++ ) {
      sim[theSlice][theRow] = (BUFTYPE*)lcalloc ( xnum, sizeof(BUFTYPE) );
    }
  }

  // read in all image data into im[], set changed[] for all slices to nil.
  for ( theSlice = 0; theSlice < numimg; theSlice++ ) {
    theSlicePtr = GetVolumeSlicePtr ( gAnatomicalVolume, theSlice );
    memcpy ( theSlicePtr, *theVolume->slices[theSlice],
       xnum * ynum * sizeof(BUFTYPE) );
    changed[theSlice] = FALSE;
  }

  for ( theSlice = 0; theSlice < numimg; theSlice++ )
    for ( theRow = 0; theRow < ynum; theRow++ )
      for ( theCol = 0; theCol < xnum; theCol++ ) {

  theIntensity = 
    GetVoxelValue ( gAnatomicalVolume, theCol, theRow, theSlice ) / 2;
  
  if ( theIntensity > sim[3][theRow][theCol]) 
    sim[3][theRow][theCol] = theIntensity;
  
  if ( theIntensity > sim[4][theSlice][theCol]) 
    sim[4][theSlice][theCol] = theIntensity;
  
  if ( theIntensity > sim[5][theRow][theSlice]) 
    sim[5][theRow][theSlice] = theIntensity;
      }

  for ( theRow = 0; theRow < ynum; theRow++ )
    for ( theCol = 0; theCol < xnum; theCol++ )
      for ( theSlice = 0; theSlice < 3; theSlice++ ) {
  sim[theSlice][theRow][theCol] = sim[theSlice+3][theRow][theCol];
      }

  // editing is always disabled in this mode (for the moment)
  editflag = FALSE;

  OutputPrint "NOTE: Editing is disabled. Do not try and enable it or I will crash.\n" EndOutputPrint;

  MRIfree ( &theVolume );
}


/* ================================================ Control point utilities */

                                       /* reads the control.dat file, 
                                          transforms all pts from RAS space 
                                          to voxel space, and adds them as 
                                          control pts */
void ProcessCtrlPtFile ( char * inDir ) {

  char theFilename [NAME_LENGTH];
  FILE * theFile;
  float theTempX, theTempY, theTempZ;
  Real theRASX, theRASY, theRASZ;
  int theVoxX, theVoxY, theVoxZ;
  int theNumPtsRead;
  char theResult;
  VoxelRef theVoxel;

  Voxel_New ( &theVoxel );

  // don't parse the file if we already have.
  if ( TRUE == gParsedCtrlPtFile )
    return;

  // create the filename
  sprintf ( theFilename, "%s/control.dat", inDir );
  
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
    theNumPtsRead = fscanf ( theFile, "%f %f %f", &theTempX, &theTempY, &theTempZ );

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
      Voxel_Set ( theVoxel, theVoxX, theVoxY, theVoxZ );
      theResult = VSpace_AddVoxel ( gCtrlPtList, theVoxel );
      if ( theResult != kVSpaceErr_NoErr )
        DebugPrint "ProcessCtrlPtFile(): Error in VSpace_AddVoxel: %s\n", 
          VSpace_GetErrorString ( theResult ) EndDebugPrint;
    }
  }

  // close the file.
  fclose ( theFile );

  // mark that we have processed the file, and shouldn't do it again.
  gParsedCtrlPtFile = TRUE;

  Voxel_Delete ( &theVoxel );
}

                                       /* writes all control points to the
                                          control.dat file in RAS space */
void WriteCtrlPtFile ( char * inDir ) {

  char theFilename [NAME_LENGTH];
  FILE * theFile;
  int theVoxX, theVoxY, theVoxZ;
  Real theRASX, theRASY, theRASZ;
  char theResult;
  int theIndex;
  int theListCount, theListIndex;
  VoxelListRef theList;
  VoxelRef theVoxel;
  
  Voxel_New ( &theVoxel );

  // create the filename.
  sprintf ( theFilename, "%s/control.dat", inDir );

  // open the file for writing.
  theFile = fopen ( theFilename, "w" );

  // check for success.
  if ( NULL == theFile ) {
    DebugPrint "WriteCtrlPtFile: Couldn't create %s for writing.\n",
      theFilename EndDebugPrint;
    return;
  }

   DebugPrint "WriteCtrlPtFile: Writing control points to %s...\n",
    theFilename EndDebugPrint;
   OutputPrint "Saving control points... " EndOutputPrint;

  // get the ctrl pts in the list...
  for ( theIndex = 0; theIndex < kVSpace_NumListsInPlane; theIndex++ ) {

    // get the list for this x value.
    theResult = VSpace_GetVoxelsInXPlane ( gCtrlPtList, theIndex, &theList );
    if ( theResult != kVSpaceErr_NoErr ) {
      DebugPrint "WriteCtrlPtFile(): Error in VSpace_GetVoxelsInXPlane: %s\n",
        VSpace_GetErrorString ( theResult ) EndDebugPrint;
      OutputPrint "\nError saving a point!\n" EndOutputPrint;
      continue;
    }

    // get the number of voxels in the list.
    theResult = VList_GetCount ( theList, &theListCount );
    if ( theResult != kVListErr_NoErr ) {
      DebugPrint "WriteCtrlPtFile(): Error in VList_GetCount: %s\n",
        VList_GetErrorString ( theResult ) EndDebugPrint;
      OutputPrint "\nError saving a point!\n" EndOutputPrint;
      continue;
    }

    // get each voxel...
    for ( theListIndex = 0; theListIndex < theListCount; theListIndex++ ) {

      theResult = VList_GetNthVoxel ( theList, theListIndex, theVoxel );
      if ( theResult != kVListErr_NoErr ) {
        DebugPrint "WriteCtrlPtFile(): Error in VList_GetNthVoxel: %s\n",
          VList_GetErrorString ( theResult ) EndDebugPrint;
  OutputPrint "\nError saving a point!\n" EndOutputPrint;
        continue;
      }
      
      // unpack it.
      theVoxX = Voxel_GetX ( theVoxel );
      theVoxY = Voxel_GetY ( theVoxel );
      theVoxZ = Voxel_GetZ ( theVoxel );

      // transform to ras space.
      VoxelToRAS ( theVoxX, theVoxY, theVoxZ,
                   &theRASX, &theRASY, &theRASZ );
      
      // write to the file
      fprintf ( theFile, "%f %f %f\n", theRASX, theRASY, theRASZ );
      
      DebugPrint "\t%d %d %d -> %f %f %f\n", 
        theVoxX, theVoxY, theVoxZ,
        theRASX, theRASY, theRASZ EndDebugPrint;
      
    }
    
  }

  // close file
  fclose ( theFile );

  DebugPrint "\tDone.\n" EndDebugPrint;
  OutputPrint " done.\n" EndOutputPrint;

  Voxel_Delete ( &theVoxel );
}


int FindNearestCtrlPt ( VoxelRef inVolumeVox, tkm_tOrientation inPlane,
       VoxelRef outCtrlPt ) {

  int theListCount, theListIndex;
  char theResult;
  VoxelRef theVoxel;
  VoxelListRef theList;
  unsigned int theDistance, theClosestDistance;
  short theClosestIndex;
  int theReturn;

  // assume not found.
  theReturn = 0;

  Voxel_New ( &theVoxel );

  // get the dimension to search in, based on our current plane, and
  // get the list of voxels for this slice.
  switch ( inPlane ) {
  case tkm_tOrientation_Coronal: 
    theResult = VSpace_GetVoxelsInZPlane ( gCtrlPtList, 
             Voxel_GetZ(inVolumeVox), &theList );
    break;
  case tkm_tOrientation_Horizontal: 
    theResult = VSpace_GetVoxelsInYPlane ( gCtrlPtList, 
             Voxel_GetY(inVolumeVox), &theList );
    break;
  case tkm_tOrientation_Sagittal: 
    theResult = VSpace_GetVoxelsInXPlane ( gCtrlPtList, 
             Voxel_GetX(inVolumeVox), &theList );
    break;
  default:
    theResult = 0;
  }
    
  if ( theResult != kVSpaceErr_NoErr ) {
    DebugPrint "FindNearestCtrlPt(): Error in VSpace_GetVoxelsInX/Y/ZPlane: %s\n", 
      VSpace_GetErrorString ( theResult ) EndDebugPrint;
    theList = NULL;
  }

  // if we got a list...
  if ( theList != NULL ) {

    // start with a large distance and no closest ctrl pt
    theClosestDistance = 0 - 1; 
    theClosestIndex = -1;

    // get the number of voxel and search through the list...
    theResult = VList_GetCount ( theList, &theListCount );
    if ( theResult != kVListErr_NoErr ) {
      DebugPrint "FindNearestCtrlPt(): Error in VList_GetCount: %s\n",
        VList_GetErrorString ( theResult ) EndDebugPrint;
      theListCount = 0;
    }
    
    for ( theListIndex = 0; theListIndex < theListCount; theListIndex++ ) {
 
      // grab a voxel
      theResult = VList_GetNthVoxel ( theList, theListIndex, theVoxel );
      if ( theResult != kVListErr_NoErr ) {
        DebugPrint "FindNearestCtrlPt(): Error in VList_GetNthVoxel: %s\n",
          VList_GetErrorString ( theResult ) EndDebugPrint;
        continue;
      }

      // get the distance to the clicked voxel...
      theDistance =
        ((Voxel_GetX(inVolumeVox) - Voxel_GetX(theVoxel)) * 
   (Voxel_GetX(inVolumeVox) - Voxel_GetX(theVoxel))) +
        ((Voxel_GetY(inVolumeVox) - Voxel_GetY(theVoxel)) * 
   (Voxel_GetY(inVolumeVox) - Voxel_GetY(theVoxel))) +
        ((Voxel_GetZ(inVolumeVox) - Voxel_GetZ(theVoxel)) * 
   (Voxel_GetZ(inVolumeVox) - Voxel_GetZ(theVoxel)));

      // if it's less than our max, mark the distance and vox index
      if ( theDistance < theClosestDistance ) {
        theClosestDistance = theDistance;
        theClosestIndex = theListIndex;
      }
    }

    // if we found a voxel
    if ( theClosestIndex != -1 ) {

      // get it back again
      theResult = VList_GetNthVoxel ( theList, theClosestIndex, theVoxel );
      if ( theResult != kVListErr_NoErr ) {
        DebugPrint "SelectCtrlPt(): Error in VList_GetNthVoxel: %s\n",
          VList_GetErrorString ( theResult ) EndDebugPrint;
        return theResult;
      }

      // return it.
      Voxel_Copy ( outCtrlPt, theVoxel );

      theResult = 1;
    }
  }

  return theResult;
}

void AddNearestCtrlPtToSelection ( VoxelRef inVolumeVox, 
           tkm_tOrientation inPlane ) {

  char theResult;
  VoxelRef theCtrlPt;
  
  Voxel_New ( &theCtrlPt );

  // if we find a nearest control point...
  if ( FindNearestCtrlPt ( inVolumeVox, inPlane, theCtrlPt ) ) {

    // add this point to the selection
    theResult = VList_AddVoxel ( gSelectionList, theCtrlPt );
    if ( theResult != kVListErr_NoErr )
      DebugPrint "SelectCtrlPt(): Error in VList_AddVoxel: %s\n",
  VList_GetErrorString ( theResult ) EndDebugPrint;
  }  

  Voxel_Delete ( &theCtrlPt );
}

void RemoveNearestCtrlPtFromSelection ( VoxelRef inVolumeVox, 
          tkm_tOrientation inPlane ) {
  
  char theResult;
  VoxelRef theCtrlPt;
  
  Voxel_New ( &theCtrlPt );

  // if we find a nearest control point...
  if ( FindNearestCtrlPt ( inVolumeVox, inPlane, theCtrlPt ) ) {

    // remove this point from selection
    theResult = VList_RemoveVoxel ( gSelectionList, theCtrlPt );
    if ( theResult != kVListErr_NoErr )
      DebugPrint "SelectCtrlPt(): Error in VList_RemoveVoxel: %s\n",
  VList_GetErrorString ( theResult ) EndDebugPrint;    
  }
  
  Voxel_Delete ( &theCtrlPt );
}

void DeselectAllCtrlPts () {

  char theResult;
  
  theResult = VList_ClearList ( gSelectionList );
  if ( theResult != kVListErr_NoErr )
    DebugPrint "DeselectAllControlPoints(): Error in VList_ClearList: %s\n",
      VList_GetErrorString ( theResult ) EndDebugPrint;

  redraw ();
}

void NewCtrlPt ( VoxelRef inVoxel ) {

  char theResult;
  
  // add the voxel to the ctrl pt space
  theResult = VSpace_AddVoxel ( gCtrlPtList, inVoxel );
  if ( theResult != kVSpaceErr_NoErr )
    DebugPrint "NewControlPoint(): Error in VSpace_AddVoxel: %s\n",
      VSpace_GetErrorString ( theResult ) EndDebugPrint;
    
  OutputPrint "Made control point (%d,%d,%d).\n",
    EXPAND_VOXEL_INT(inVoxel)  EndOutputPrint;

  /* write it to the control point file. */
  WriteVoxelToControlFile( tfname, inVoxel );
}

                                        /* remove the selected control points 
                                           from the control point space */
void DeleteSelectedCtrlPts () {

  char theResult;
  int theCount, theIndex;
  VoxelRef theCtrlPt;

  Voxel_New ( &theCtrlPt );

  // get the number of selected points we have
  theResult = VList_GetCount ( gSelectionList, &theCount );
  if ( theResult != kVListErr_NoErr ) {
    DebugPrint "Error in VList_GetCount: %s\n",
      VList_GetErrorString ( theResult ) EndDebugPrint;
    return;
  }

  // for each one...
  for ( theIndex = 0; theIndex < theCount; theIndex++ ) {

    // get it
    theResult = VList_GetNthVoxel ( gSelectionList, theIndex, theCtrlPt );
    if ( theResult != kVListErr_NoErr ) {
      DebugPrint "Error in VList_GetNthVoxel: %s\n",
        VList_GetErrorString ( theResult ) EndDebugPrint;
      continue;
    }
    
    // set its value in the space to 0
    theResult = VSpace_RemoveVoxel ( gCtrlPtList, theCtrlPt );
    if ( theResult != kVSpaceErr_NoErr ) {
      DebugPrint "DeleteSelectedCtrlPts(): Error in VSpace_RemoveVoxel: %s\n",
        VSpace_GetErrorString ( theResult ) EndDebugPrint;
      continue;
    }        
  }

  // remove all pts from the selection list
  theResult = VList_ClearList ( gSelectionList );
  if ( theResult != kVListErr_NoErr ) {
    DebugPrint "Error in VList_ClearList: %s\n",
      VList_GetErrorString ( theResult ) EndDebugPrint;
    return;
  }

  Voxel_Delete ( &theCtrlPt );

  redraw ();
}

// ===========================================================================

// ========================================================= SELECTING REGIONS

void InitSelectionModule () {

  char theErr;
  theErr = VSpace_New ( &gSelectedVoxels );
  if ( theErr != kVSpaceErr_NoErr ) {
    DebugPrint "InitSelectionModule(): Error in VSpace_Init: %s\n",
      VSpace_GetErrorString ( theErr ) EndDebugPrint;
    gSelectedVoxels = NULL;
  }

  isDisplaySelectedVoxels = TRUE;
}

void DeleteSelectionModule () {

  char theErr;

  if ( NULL == gSelectedVoxels )
    return;
  
  theErr = VSpace_Delete ( &gSelectedVoxels );
  if ( theErr != kVSpaceErr_NoErr )
    DebugPrint "DeleteSelectionModule(): Error in VSpace_Delete: %s\n",
      VSpace_GetErrorString(theErr) EndDebugPrint;
}

void AllowSelectionModuleToRespondToClick ( VoxelRef inScreenVoxel ) {

}

void AddVoxelToSelection ( VoxelRef inVoxel ) {

  char theErr;

  if ( NULL == gSelectedVoxels )
    return;
  
  theErr = VSpace_AddVoxel ( gSelectedVoxels, inVoxel );
  if ( theErr != kVSpaceErr_NoErr )
    DebugPrint "AddVoxelToSelection(): Error in VSpace_AddVoxel: %s\n",
      VSpace_GetErrorString(theErr) EndDebugPrint;
  
}

void RemoveVoxelFromSelection ( VoxelRef inVoxel ) {

  char theErr;

  if ( NULL == gSelectedVoxels )
    return;
  
  theErr = VSpace_RemoveVoxel ( gSelectedVoxels, inVoxel );
  if ( theErr != kVSpaceErr_NoErr &&
       theErr != kVSpaceErr_VoxelNotInSpace )
    DebugPrint "RemoveVoxelFromSelection(): Error in VSpace_RemoveVoxel: %s\n",
      VSpace_GetErrorString(theErr) EndDebugPrint;
}


void ClearSelection () {

  char theErr;

  if ( NULL == gSelectedVoxels )
    return;
  
  theErr = VSpace_Clear ( gSelectedVoxels );
  if ( theErr != kVSpaceErr_NoErr &&
       theErr != kVSpaceErr_VoxelNotInSpace )
    DebugPrint "ClearSelection(): Error in VSpace_Clear: %s\n",
      VSpace_GetErrorString(theErr) EndDebugPrint;

  MWin_RedrawAll( gMeditWindow );
}

void SaveSelectionToLabelFile ( char * inFileName ) {

  LABEL * theLabel;
  LABEL_VERTEX *theVertex;
  int theNumVoxels, theGlobalVoxelIndex, theListVoxelIndex, theListIndex,
    theListCount;
  Real theRASX, theRASY, theRASZ;
  VoxelRef theAnatomicalVoxel;
  VoxelListRef theList;
  int theLabelError;
  char theError;
  
  DebugPrint "SaveSelectionToLabelFile ( %s )\n",
    inFileName EndDebugPrint;

  Voxel_New ( &theAnatomicalVoxel );

  // get the number of selected voxels we have.
  theNumVoxels = 0;
  for ( theListIndex = 0; 
  theListIndex < kVSpace_NumListsInPlane; 
  theListIndex++ ) {
    theError = VSpace_GetVoxelsInXPlane ( gSelectedVoxels, 
            theListIndex, &theList );
    if ( kVSpaceErr_NoErr != theError ) {
      DebugPrint "\tError in VSpace_GetVoxelsInXPlane (%d): %s\n",
  theError, VSpace_GetErrorString ( theError ) EndDebugPrint;
      Voxel_Delete ( &theAnatomicalVoxel );
      return;
    }

    theError = VList_GetCount ( theList, &theListCount );
    if ( kVListErr_NoErr != theError ) {
      DebugPrint "\tError in VList_GetCount (%d): %s\n",
  theError, VList_GetErrorString ( theError ) EndDebugPrint;
      Voxel_Delete ( &theAnatomicalVoxel );
      return;
    }

    theNumVoxels += theListCount;
  }

  // allocate a label file with that number of voxels, our subject name,
  // and the passed in label file name.
  theLabel = LabelAlloc ( theNumVoxels, pname, inFileName );
  if ( NULL == theLabel ) {
    DebugPrint "\tCouldn't allocate label.\n" EndDebugPrint;
    OutputPrint "ERROR: Couldn't save label.\n" EndOutputPrint;
    Voxel_Delete ( &theAnatomicalVoxel );
    return;
  }

  // set the number of points in the label
  theLabel->n_points = theNumVoxels;

  // for every list in a plane of the space... 
  theGlobalVoxelIndex = 0;
  for ( theListIndex = 0; 
  theListIndex < kVSpace_NumListsInPlane; 
  theListIndex++ ) {

    // get the list
    theError = VSpace_GetVoxelsInXPlane ( gSelectedVoxels, 
            theListIndex, &theList );
    if ( kVSpaceErr_NoErr != theError ) {
      DebugPrint "\tError in VSpace_GetVoxelsInXPlane (%d): %s\n",
  theError, VSpace_GetErrorString ( theError ) EndDebugPrint;
      Voxel_Delete ( &theAnatomicalVoxel );
      return;
    }

    // get the num of voxels in the list.
    theError = VList_GetCount ( theList, &theListCount );
    if ( kVListErr_NoErr != theError ) {
      DebugPrint "\tError in VList_GetCount (%d): %s\n",
  theError, VList_GetErrorString ( theError ) EndDebugPrint;
      Voxel_Delete ( &theAnatomicalVoxel );
      return;
    }

    // note that this is only the index of the voxel within a list, not
    // global count. for each voxel in the list...
    for ( theListVoxelIndex = 0;
    theListVoxelIndex < theListCount;
    theListVoxelIndex++ ) {
      
      // get a voxel.
      theError = VList_GetNthVoxel ( theList, 
             theListVoxelIndex, theAnatomicalVoxel );
      if ( kVListErr_NoErr != theError ) {
  DebugPrint "\tError in VList_GetNthVoxel (%d): %s\n",
    theError, VList_GetErrorString ( theError ) EndDebugPrint;
  Voxel_Delete ( &theAnatomicalVoxel );
  return;
      }

      // get a ptr the vertex in the label file. use the global count to 
      // index.
      theVertex = &(theLabel->lv[theGlobalVoxelIndex]);
      
      // convert to ras
      VoxelToRAS ( EXPAND_VOXEL_INT(theAnatomicalVoxel),
       &theRASX, &theRASY, &theRASZ );

      // set the vertex
      theVertex->x = theRASX;
      theVertex->y = theRASY;
      theVertex->z = theRASZ;
      
      // set the vno to -1, which is significant somewhere outside the realm
      // of tkmedit. set stat value to something decent and deleted to not
      theVertex->vno = -1;
      theVertex->stat = 0;
      theVertex->deleted = FALSE;

      // inc our global count.
      theGlobalVoxelIndex ++;
    }
  }

  // write the file
  theLabelError = LabelWrite ( theLabel, inFileName );
  if ( NO_ERROR != theLabelError ) {
    DebugPrint "Error in LabelWrite().\n" EndDebugPrint;
    OutputPrint "ERROR: Couldn't write label to file.\n" EndOutputPrint;
    Voxel_Delete ( &theAnatomicalVoxel );
    return;
  }

  // free it
  LabelFree ( &theLabel );

  Voxel_Delete ( &theAnatomicalVoxel );
}

void LoadSelectionFromLabelFile ( char * inFileName ) {

  LABEL *theLabel;
  LABEL_VERTEX *theVertex;
  int theNumVoxels, theVoxelIndex;
  int theVoxX, theVoxY, theVoxZ;
  VoxelRef theVoxel;

  Voxel_New ( &theVoxel );

  DebugPrint "LoadSelectionFromLabelFile ( %s )\n", inFileName EndDebugPrint;

  // read in the label
  theLabel = LabelRead ( pname, inFileName );
  if ( NULL == theLabel ) {
    DebugPrint "\tError reading label file.\n" EndDebugPrint;
    OutputPrint "ERROR: Couldn't read the label.\n" EndOutputPrint;
    
  } else {
    
    // for each vertex in there...
    theNumVoxels = theLabel->max_points;
    DebugPrint "Loading %d points.\n", theNumVoxels EndDebugPrint;
    for ( theVoxelIndex = 0; theVoxelIndex < theNumVoxels; theVoxelIndex++ ) {
      
      // get the vertex.
      theVertex = &(theLabel->lv[theVoxelIndex]);
      
      // only process verticies that arn't deleted.
      if ( !(theVertex->deleted) ) {
  
  // covert from ras to voxel
  RASToVoxel ( theVertex->x, theVertex->y, theVertex->z,
         &theVoxX, &theVoxY, &theVoxZ );
  
  // add to the selection
  Voxel_Set ( theVoxel, theVoxX, theVoxY, theVoxZ );
  AddVoxelToSelection ( theVoxel );
  
      }
    }

    // dump the label
    LabelFree ( &theLabel );
  }  

  /* set the window's selection again to force a redraw. */
  MWin_SetSelectionSpace( gMeditWindow, -1, gSelectedVoxels );

  Voxel_Delete ( &theVoxel );
}

void GraphSelectedRegion () {

  char theError;
  int theListIndex, theListVoxelIndex, theListCount;
  VoxelRef theAnatomicalVoxel = NULL;
  VoxelListRef theList = NULL;
  FunV_tErr eFunctional = FunV_tErr_NoError;

  Voxel_New ( &theAnatomicalVoxel );

  // clear the functional display list.
  eFunctional = FunV_BeginSelectionRange( gFunctionalVolume );
  if( FunV_tErr_NoError != eFunctional ) {
    DebugPrint "GraphSelectedRegion(): Error accessing functional volume.\n"
      EndDebugPrint;
    goto cleanup;
  }

  // for every list in a plane of the space... 
  for ( theListIndex = 0; 
  theListIndex < kVSpace_NumListsInPlane; 
  theListIndex++ ) {

    // get the list
    theError = VSpace_GetVoxelsInXPlane ( gSelectedVoxels, 
            theListIndex, &theList );
    if ( kVSpaceErr_NoErr != theError ) {
      DebugPrint "\tError in VSpace_GetVoxelsInXPlane (%d): %s\n",
  theError, VSpace_GetErrorString ( theError ) EndDebugPrint;
      goto cleanup;
    }

    // get the num of voxels in the list.
    theError = VList_GetCount ( theList, &theListCount );
    if ( kVListErr_NoErr != theError ) {
      DebugPrint "\tError in VList_GetCount (%d): %s\n",
  theError, VList_GetErrorString ( theError ) EndDebugPrint;
      goto cleanup;
    }

    for ( theListVoxelIndex = 0;
    theListVoxelIndex < theListCount;
    theListVoxelIndex++ ) {
      
      // get a voxel.
      theError = VList_GetNthVoxel ( theList, 
             theListVoxelIndex, theAnatomicalVoxel );
      if ( kVListErr_NoErr != theError ) {
  DebugPrint "\tError in VList_GetNthVoxel (%d): %s\n",
    theError, VList_GetErrorString ( theError ) EndDebugPrint;
  goto cleanup;
      }

      // add it to the functional display list.
      eFunctional = 
  FunV_AddAnatomicalVoxelToSelectionRange( gFunctionalVolume, 
             theAnatomicalVoxel );
      if( FunV_tErr_NoError != eFunctional ) {
  DebugPrint"GraphSelectedRegion(): Error accessing functional volume.\n"
    EndDebugPrint;
  goto cleanup;
      }
    }
  }

  // finish the list
  eFunctional = FunV_EndSelectionRange( gFunctionalVolume );
  if( FunV_tErr_NoError != eFunctional ) {
    DebugPrint "GraphSelectedRegion(): Error accessing functional volume.\n"
      EndDebugPrint;
    goto cleanup;
  }

 cleanup:

  if ( NULL != theAnatomicalVoxel )
    Voxel_Delete ( &theAnatomicalVoxel );
}

// ===========================================================================

// ============================================================ EDITING VOXELS

void EditVoxelInRange( VoxelRef     ipVoxel, 
           tVolumeValue inLow, 
           tVolumeValue inHigh, 
           tVolumeValue inNewValue ) {
  
  tVolumeValue value    = 0;
  tVolumeValue newValue = 0;

  /* get the value at this point. */
  value    = GetVoxelValue( gAnatomicalVolume, EXPAND_VOXEL_INT(ipVoxel) );
  newValue = value;

  /* if it's in the range... */
  if( value >= inLow && value <= inHigh ) {

    /* set a new value */
    newValue = inNewValue;
  }

  /* if values are different, set and add to undo list. */
  if( value != newValue ) {
    SetVoxelValue( gAnatomicalVolume, EXPAND_VOXEL_INT(ipVoxel), newValue );
    AddVoxelAndValueToUndoList( ipVoxel, value );
  }

  /* mark this slice as changed. */
  changed[ Voxel_GetZ(ipVoxel) ] = TRUE;
  editedimage = imnr0 + Voxel_GetZ(ipVoxel);
}

// ===========================================================================


/* ============================================= Coordinate transformations */

inline
char IsVoxelInBounds ( int x, int y, int z ) {

  if ( x < 0 || y < 0 || z < 0 ||
       x >= xnum || y >= ynum || z >= xnum ) {

    DebugPrint "!!! Voxel coords out of bounds: (%d, %d, %d)\n",
      x, y, z EndDebugPrint;

    return FALSE;
  }

  return TRUE;
}

inline
char IsRASPointInBounds ( Real x, Real y, Real z ) {

  if ( x < xx0 || x > xx1 || y < yy0 || y > yy1 || z < zz0 || z > zz1 ) {

    DebugPrint "!!! RAS coords out of bounds: (%2.1f, %2.1f, %2.1f)\n",
      x, y, z EndDebugPrint;
    DebugPrint "    x: %2.1f, %2.1f  y: %2.1f, %2.1f  z: %2.1f, %2.1f\n",
      xx0, xx1, yy0, yy1, zz0, zz1 EndDebugPrint;

    return FALSE;
  }

  return TRUE;
}

void RASToVoxel ( Real x, Real y, Real z,        // incoming ras coords
                  int *xi, int *yi, int *zi ) {  // outgoing voxel coords

  // call the function in mritransform.
  trans_RASToVoxelIndex ( x, y, z, xi, yi, zi );

  // check our bounds...
  if ( ! IsVoxelInBounds ( *xi, *yi, *zi ) ) {

    // try not to crash.
    *xi = *yi = *zi = 0;
  }
}               

void VoxelToRAS ( int xi, int yi, int zi,        // incoming voxel coords
                  Real *x, Real *y, Real *z ) {  // outgoing RAS coords

  // check our bounds...
  if ( ! IsVoxelInBounds ( xi, yi, zi ) ) {

    // try not to crash.
    *x = *y = *z = 0;
    return;
  }

  // call the function in mritransform.
  trans_VoxelIndexToRAS ( xi, yi, zi, x, y, z );
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

  tVolumeRef theVolume = NULL;

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

void SetVolumeColorScale ( float ifThreshold, float ifSquash ) {

  float thresh = 0;
  float squash = 0;
  float scale = 0;
  float i = 0;
  int   nColor = 0;
  char sTclArguments [128];

  gfVolumeColorSquash    = ifSquash;
  gfVolumeColorThreshold = ifThreshold;

  /* generate the look up table */
  scale = knMaxVolumeValue;
  thresh = gfVolumeColorThreshold;
  squash = gfVolumeColorSquash;
  for( nColor = 0; nColor < knNumVolumeValues; nColor++ ) { 
    i = (float) nColor;
    gfaVolumeColors[nColor] =
      (1.0 / (1.0 + exp(-squash * ((i/scale) - thresh)))) * scale + 0.5;
    if( gfaVolumeColors[nColor] > knMaxVolumeValue)
      gfaVolumeColors[nColor] = knMaxVolumeValue;
      
  }

  /* tell window to redraw */
  if( NULL != gMeditWindow ) {
    MWin_RedrawAll( gMeditWindow );
  }
  
  /* update tcl window */
  sprintf ( sTclArguments, "%2f %2f", 
      gfVolumeColorThreshold, gfVolumeColorSquash );
  tkm_SendTclCommand ( tkm_tTclCommand_UpdateVolumeColorScale, 
           sTclArguments );
}

void GetVolumeColor ( tVolumeValue iucValue,
          unsigned char *oucRed, 
          unsigned char *oucGreen,
          unsigned char *oucBlue ) {

  *oucRed   = gfaVolumeColors[iucValue];
  *oucGreen = gfaVolumeColors[iucValue];
  *oucBlue  = gfaVolumeColors[iucValue];
}

// ============================================================== PARCELLATION

void LoadParcellationVolume ( char* isVolumeDirWithPrefix,
            char* isColorFileName ) {

  char  sFileName[256] = "";
  FILE* pFile = NULL;
  int   nSliceBegin = 0;
  int   nSliceEnd = 0;
  int   type = 0;
  int   nXDimension = 0;
  int   nYDimension = 0;
  int   nBufferSize = 0;
  tVolumeValue* pBuffer = NULL;
  int   nSlice = 0;
  tVolumeValue* pSlice = NULL;
  int nColor;
  char bGood = FALSE;
  char sLine[1024] = "";
  int nRed = 0;
  int nGreen = 0;
  int nBlue = 0;
  int nBiggestIndex = 0;

  /* free existing volume. */
  if( NULL != gParcellationVolume ) {
    free( gParcellationVolume );
    gParcellationVolume = NULL;
  }

  /* check to see if the info file exist. */
  sprintf( sFileName, "%s.info", isVolumeDirWithPrefix );
  pFile = fopen( sFileName, "rb" );
  if( NULL == pFile ) {
    DebugPrint "LoadParcellationVolume: Couldn't open %s\n",
      sFileName EndDebugPrint;
    OutputPrint "Error finding parcellation data.\n" EndOutputPrint;
    goto cleanup;
  }

  /* read the info file. only get the xnum and ynum values. */
  fscanf ( pFile, "%*s %d", &nSliceBegin);
  fscanf ( pFile, "%*s %d", &nSliceEnd);
  fscanf ( pFile, "%*s %d", &type);
  fscanf ( pFile, "%*s %d", &nXDimension);
  fscanf ( pFile, "%*s %d", &nYDimension);
 
  /* if dimensions are not equal to volume dimensions, bail out. */
  if( nSliceBegin    != imnr0
      || nSliceEnd   != imnr1
      || type        != ptype
      || nXDimension != xnum
      || nYDimension != ynum ) {
    DebugPrint "LoadParcellationVolume: dimensions didn't match.\n"
      EndDebugPrint;
    OutputPrint "Error matching parcellation volume description.\n"
      EndOutputPrint;
    goto cleanup;
  }

  /* init the volume. */
  InitVolume( &gParcellationVolume, nXDimension );
  if( NULL == gParcellationVolume ) {
    DebugPrint "LoadParcellationVolume: Volume allocation failed.\n"
      EndDebugPrint;
    OutputPrint "Out of memory.\n" EndOutputPrint;
    goto cleanup;
  }

  /* allocate temp storage. */
  nBufferSize = nXDimension * nYDimension;
  pBuffer = (tVolumeValue*) malloc (nBufferSize*sizeof(tVolumeValue));
  if( NULL == pBuffer ) {
    DebugPrint "LoadParcellationVolume: Temp buffer allocation failed.\n"
      EndDebugPrint;
    OutputPrint "Out of memory.\n" EndOutputPrint;
    goto cleanup;
  }

  /* read the parcellation volume into the volume. */
  for( nSlice = 0; nSlice <= nSliceEnd-nSliceBegin; nSlice++ ) {

    /* make file name */
    sprintf( sFileName, "%s%03d", 
       isVolumeDirWithPrefix, nSlice + nSliceBegin );

    /* open file. */
    pFile = fopen( sFileName, "rb" );
    if( NULL == pFile ) {
      DebugPrint "LoadParcellationVolume: Couldn't open %s\n", 
  sFileName EndDebugPrint;
      OutputPrint "Error opening parcellation volume slice file %d.\n",
  nSlice EndOutputPrint;
      goto cleanup;
    }

    /* read it into temp memory */
    fread( pBuffer, sizeof(tVolumeValue), nBufferSize, pFile );

    /* copy to the volume*/
    pSlice = GetVolumeSlicePtr( gParcellationVolume, nSlice );
    memcpy( pSlice, pBuffer, nBufferSize );

    fclose( pFile );
    pFile = NULL;
  }

  /* free temp storage */
  free( pBuffer );
  pBuffer = NULL;

  /* open the index file. */
  pFile = fopen( isColorFileName, "rb" );
  if( NULL == pFile ) {
    DebugPrint "LoadParcellationVolume: Couldn't open %s\n", 
  isColorFileName EndDebugPrint;
    OutputPrint "Error opening color file %s\n", isColorFileName
      EndOutputPrint;
    
    /* trash volume */
    free( gParcellationVolume );
    gParcellationVolume = NULL;

    goto cleanup;
  }

  /* scan to see how many colors */
  nBiggestIndex = 0;
  while( !feof( pFile ) ) {
    fgets( sLine, 1024, pFile );
    bGood = sscanf( sLine, "%d %*s %d %d %d %*s",
        &nColor, &nRed, &nGreen, &nBlue );
    if( !bGood ) {
      DebugPrint "Error reading color file after %d lines\n",
  gNumParcellationColors EndDebugPrint;
      goto cleanup;
    }  
    if( nColor > nBiggestIndex )
      nBiggestIndex = nColor;
  }
  fclose( pFile );
    pFile = NULL;

  /* we actually counted number of entries there, so make one more so we can
     keep the biggest index. */
  gNumParcellationColors = nBiggestIndex + 1;

  /* allocation color storage*/
  if( NULL != gParcellationColors )
    free( gParcellationColors );
  gParcellationColors = (tColorEntry*) malloc ( gNumParcellationColors * 
            sizeof(tColorEntry) );

  /* read the indicies in. */
  pFile = fopen( isColorFileName, "rb" );
  while( !feof( pFile ) ) {

    fgets( sLine, 1024, pFile );
    bGood = sscanf( sLine, "%d %*s %d %d %d %*s",
        &nColor, &nRed, &nGreen, &nBlue );

    gParcellationColors[nColor].mRed   = nRed;
    gParcellationColors[nColor].mGreen = nGreen;
    gParcellationColors[nColor].mBlue  = nBlue;
  }
  fclose( pFile );
  pFile = NULL;

  /* set parcellation volume in window */
  MWin_SetParcellationVolume( gMeditWindow, -1, gParcellationVolume,
            nXDimension );

 cleanup:

  if( NULL != pFile )
    fclose( pFile );

  if( NULL != pBuffer ) 
    free( pBuffer );
}

void GetParcellationColor ( VoxelRef ipVoxel,
          unsigned char* oucRed,
          unsigned char* oucGreen,
          unsigned char* oucBlue ) {

  tVolumeValue ucIndex = 0;

  /* get the index from the volume */
  ucIndex = GetVoxelValue( gParcellationVolume, EXPAND_VOXEL_INT(ipVoxel) );

  /* get the color out of the color map */
  if( ucIndex != 0 && ucIndex <= gNumParcellationColors ) {
    *oucRed   = gParcellationColors[ucIndex].mRed;
    *oucGreen = gParcellationColors[ucIndex].mGreen;
    *oucBlue  = gParcellationColors[ucIndex].mBlue;
  } else {
    *oucRed   = 0;
    *oucGreen = 0;
    *oucBlue  = 0;
  }  
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
        VoxelRef inVoxel, tVolumeValue inValue ) {

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
  Voxel_New ( &(this->mVoxel) );
  Voxel_Copy ( this->mVoxel, inVoxel );

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
  Voxel_Delete ( &(this->mVoxel) );
  
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
    EXPAND_VOXEL_INT(this->mVoxel), this->mValue EndDebugPrint;

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

void AddVoxelAndValueToUndoList ( VoxelRef inVoxel, int inValue ) {

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
          EXPAND_VOXEL_INT(theEntryToUndo->mVoxel) );

  // create an entry for it.
  NewUndoEntry ( &theUndoneEntry, theEntryToUndo->mVoxel, theVoxelValue );

  // set the voxel value.
  SetVoxelValue ( gAnatomicalVolume,
      EXPAND_VOXEL_INT(theEntryToUndo->mVoxel),
      theEntryToUndo->mValue );

  // pass back the new entry.
  *outNewEntry = theUndoneEntry;
}

void PrintEntryWrapper ( xUndL_tEntryPtr inEntry ) {

  PrintUndoEntry ( (UndoEntryRef) inEntry );
}

void tkm_ConvertVolumeToRAS ( VoxelRef inVolumeVox, VoxelRef outRASVox ) {

  Real theRASX, theRASY, theRASZ;

  VoxelToRAS ( EXPAND_VOXEL_INT ( inVolumeVox ), 
         &theRASX, &theRASY, &theRASZ );
  Voxel_SetFloat ( outRASVox, (float)theRASX, (float)theRASY, (float)theRASZ );
}

void tkm_ConvertVolumeToTal ( VoxelRef inVolumeVox, VoxelRef outTalVox ) {

  Real theRASX, theRASY, theRASZ;
  Real theTalX, theTalY, theTalZ;

  if ( transform_loaded ) {
    VoxelToRAS ( EXPAND_VOXEL_INT ( inVolumeVox ), 
     &theRASX, &theRASY, &theRASZ );
    transform_point ( linear_transform, theRASX, theRASY, theRASZ,
                      &theTalX, &theTalY, &theTalZ );
    Voxel_SetFloat ( outTalVox, 
         (float)theTalX, (float)theTalY, (float)theTalZ );
  } else {
    Voxel_Set ( outTalVox, 0, 0, 0 );
  }
}

void tkm_CovertVolumeToTal ( VoxelRef inVolumeVox, VoxelRef outTalVox ) {

}

tVolumeValue tkm_GetVolumeValue ( tVolumeRef inVolume,
           VoxelRef inVoxel ) {

  return GetVoxelValue ( inVolume, EXPAND_VOXEL_INT(inVoxel) );
}

void tkm_GetAnatomicalVolumeColor( tVolumeValue inValue,
           float* outRed,
           float* outGreen, 
           float* outBlue ) {
  
  unsigned char ucRed, ucGreen, ucBlue;

  GetVolumeColor( inValue, &ucRed, &ucGreen, &ucBlue );

  *outRed   = (float)ucRed   / (float)knMaxVolumeValue;
  *outGreen = (float)ucGreen / (float)knMaxVolumeValue;
  *outBlue  = (float)ucBlue  / (float)knMaxVolumeValue;


}

void tkm_AddNearestCtrlPtToSelection ( VoxelRef inVolumeVox, 
               tkm_tOrientation inPlane ) {

  AddNearestCtrlPtToSelection ( inVolumeVox, inPlane );
}

void tkm_RemoveNearestCtrlPtFromSelection ( VoxelRef inVolumeVox, 
              tkm_tOrientation inPlane ) {

  RemoveNearestCtrlPtFromSelection ( inVolumeVox, inPlane );
}

void tkm_NewCtrlPt () {

  VoxelRef theCursor = NULL;
  Voxel_New ( &theCursor );
  MWin_GetCursor ( gMeditWindow, theCursor );
  NewCtrlPt ( theCursor );
  Voxel_Delete ( &theCursor );
}

void tkm_DeselectAllCtrlPts () {

  DeselectAllCtrlPts ();
}

void tkm_DeleteSelectedCtrlPts () {

  DeleteSelectedCtrlPts ();
}

void tkm_WriteControlFile () {

  WriteCtrlPtFile ( tfname );
}

void tkm_EditVoxelInRange( VoxelRef     inVolumeVox, 
         tVolumeValue inLow, 
         tVolumeValue inHigh, 
         tVolumeValue inNewValue ) {

  EditVoxelInRange( inVolumeVox, inLow, inHigh, inNewValue );

}

void tkm_SelectVoxel ( VoxelRef inVolumeVox ) {

  AddVoxelToSelection ( inVolumeVox );
}

void tkm_DeselectVoxel ( VoxelRef inVolumeVox ) {

  RemoveVoxelFromSelection (inVolumeVox );
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

void tkm_WriteVoxelToControlFile ( VoxelRef inVolumeVox ) {

  WriteVoxelToControlFile ( tfname, inVolumeVox );
}

void tkm_WriteVoxelToEditFile ( VoxelRef inVolumeVox ) {

  WriteVoxelToEditFile ( tfname, inVolumeVox );
}

void tkm_ReadCursorFromEditFile () {

  ReadCursorFromEditFile ( tfname );
}

char* tkm_GetSubjectName() {

  return pname;
}

char* tkm_GetVolumeName() {

  return imtype;
}

char* tkm_GetAuxVolumeName() {

  return imtype2;
}

void tkm_GetParcellationColor( VoxelRef inVoxel, 
             float* outRed,
             float* outGreen, 
             float* outBlue ) {
  
  unsigned char ucRed, ucGreen, ucBlue;
  
  GetParcellationColor( inVoxel, &ucRed, &ucGreen, &ucBlue );

  *outRed   = ucRed   / knMaxVolumeValue;
  *outGreen = ucGreen / knMaxVolumeValue;
  *outBlue  = ucBlue  / knMaxVolumeValue;  
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

  /* display status */
  "ShowVolumeCoords",
  "ShowRASCoords",
  "ShowTalCoords",
  "ShowAuxValue",
  "ShowFuncCoords",
  "ShowFuncValue",
  "ShowOverlayOptions",
  "ShowTimeCourseOptions",

  /* interface configuration */
  "wm geometry .",
  "CsurfInterface",
  "ErrorDlog"
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


